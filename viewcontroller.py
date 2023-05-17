import os
import sys
import glob
import json

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import (
      QVBoxLayout, QWidget, QPushButton, QTabWidget, QFormLayout,
      QHBoxLayout, QLineEdit, QLabel, QFileDialog
      )
from matplotlib.backends.backend_qtagg import (        FigureCanvasQTAgg,
        NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
import h5py


from util import index_phase
from model import datamodel 
from stripeview import stripeview
from globalview import globalview
from labeling_engine import labeler
from cif_view import CIFView
from phase_diagram import PhaseDiagramView, PhaseDiagramList
from popup import Popup


class TopLevelWindow(QtWidgets.QMainWindow):
    def __init__(self, h5_path = "data/AL_23F4_Bi-Ti-O_run_01_0_all_1d.h5",
                 cif_path = "/Users/ming/Desktop/Code/SARA.jl/BiTiO/cifs/sticks.csv" ):
        super().__init__()
        # Temperature fixed
        self.stripeview = stripeview()
        self.globalview = globalview()
        self.h5_path = h5_path
        self.cif_path = cif_path
        self.model = datamodel(self.h5_path) 
        self.labeler = labeler(self.cif_path)
        self.cifview = CIFView([phase.name for phase in self.labeler.phases])
        self.popup = Popup() 

        self.ind = 0
        self.update(self.ind)
        self.stripeview.plot(self.model[self.ind])
        self.globalview.plot(
              self.model.dwells, self.model.tpeaks,
              self.model.labeled_dwells, self.model.labeled_tpeaks,
              self.model.current_dwell, self.model.current_tpeak,
              self.model.x, self.model.y,
              self.model.labeled_x, self.model.labeled_y,
              self.model.current_x, self.model.current_y
        )

        self.phase_diagram_view = PhaseDiagramView()
        self.phase_diagram_list = PhaseDiagramList()
        pd_layout = QHBoxLayout()
        pd_layout.addWidget(self.phase_diagram_view)
        pd_layout.addWidget(self.phase_diagram_list)
        pd_layout.setStretch(0, 3)
        pd_layout.setStretch(1, 1)

        self.globalview.picked.connect(self.update)
        self.cifview.checked.connect(self.update_sticks)
        self.cifview.add.connect(self.add_to_phase_diagram)
        self.phase_diagram_list.checked.connect(self.update_pd_plot) # FIXME
        self.popup.set_clicked.connect(self.update_labeler_hyperparams)

        label_button = QPushButton()
        label_button.setText("Label")
        label_button.clicked.connect(self.label_button_clicked)

        save_button = QPushButton()
        save_button.setText("Save")
        save_button.clicked.connect(self.save_button_clicked)

        back_button = QPushButton()
        back_button.setText("Back")
        back_button.clicked.connect(lambda: self.change_ind(-1))

        next_button = QPushButton()
        next_button.setText("Next")
        next_button.clicked.connect(lambda: self.change_ind(1))

        browse_button = QPushButton()
        browse_button.setText("Browse data file")
        browse_button.clicked.connect(self.browse_button_clicked)

        browse_cif_button = QPushButton()
        browse_cif_button.setText("Browse CIF directory")
        browse_cif_button.clicked.connect(self.browse_cif_button_clicked)
 
        save_progress_button = QPushButton()
        save_progress_button.setText("Save Progress")
        save_progress_button.clicked.connect(self.save_progress_clicked)

        load_progress_button = QPushButton()
        load_progress_button.setText("Load Progress")
        load_progress_button.clicked.connect(self.load_progress_clicked)

        labeler_setting_button = QPushButton()
        labeler_setting_button.setText("Labeler Settings")
        labeler_setting_button.clicked.connect(self.labeler_setting_clicked)
       
        browse_button_layout = QHBoxLayout()
        browse_button_layout.addWidget(browse_button)
        browse_button_layout.addWidget(browse_cif_button)
        browse_button_layout.addWidget(save_progress_button)
        browse_button_layout.addWidget(load_progress_button)
        browse_button_layout.addWidget(labeler_setting_button)


        button_layout = QHBoxLayout()
        button_layout.addWidget(label_button)
        button_layout.addWidget(save_button)
        button_layout.addWidget(back_button)
        button_layout.addWidget(next_button)

        main_fig_layout = QVBoxLayout()
        main_fig_layout.addLayout(browse_button_layout)
        main_fig_layout.addWidget(browse_button)
        main_fig_layout.addWidget(self.stripeview)
        main_fig_layout.addWidget(self.globalview)
        main_fig_layout.addLayout(button_layout)

        outer_layout = QHBoxLayout()
        outer_layout.addLayout(main_fig_layout)
        outer_layout.addWidget(self.cifview)
        outer_layout.setStretch(0, 3)
        outer_layout.setStretch(1, 1)

        widget = QWidget()
        widget.setLayout(outer_layout)
        pd_widget = QWidget()
        pd_widget.setLayout(pd_layout)
        self.tabs = QTabWidget()
        self.tabs.addTab(widget, "Labeler")
        self.tabs.addTab(pd_widget, "Phase Map")
        self.tabs.currentChanged.connect(self.update_pd_tab)
        self.setCentralWidget(self.tabs)

    def browse_button_clicked(self):
        self.file_name, _ = QFileDialog.getOpenFileName(None, "Open", "", "")
        if self.file_name.endswith("h5"):
            self.model = datamodel(self.file_name)
            self.ind = 0
            self.update(self.ind)

    def browse_cif_button_clicked(self):
        #self.cif_dir_name = QFileDialog.getExistingDirectory(self, "Select Directory") 
        self.cif_csv_fn, _ = QFileDialog.getOpenFileName(None, "Open", "", "")
        if self.cif_csv_fn.endswith("csv"):
            self.labeler = labeler(self.cif_csv_fn)
            self.cifview.update_cif_list([phase.name for phase in self.labeler.phases])

    def save_progress_clicked(self):
        self.save_fn, _ = QFileDialog.getSaveFileName(self, 'Save File', "", "JSON Files (*.json)")
        storing_ds = {}
        storing_ds["phases_diagram"] = self.model.get_dict_for_phase_diagram()
        storing_ds["phases"] = self.model.phases
        storing_ds["cif_path"] = self.cif_path
        storing_ds["h5_path"] = self.h5_path
        with open(self.save_fn, 'w') as f:
            json.dump(storing_ds, f)

    def load_progress_clicked(self):
        self.load_fn, _ = QFileDialog.getOpenFileName(None, "Open", "", "JSON Files (*.json)")
        if self.load_fn.endswith(".json"):
            with open(self.load_fn, 'r') as f:
                load_meta_data = json.load(f)

            self.model = datamodel(load_meta_data["h5_path"])
            self.labeler = labeler(load_meta_data["cif_path"])
            self.model.phases = load_meta_data["phases"]
            self.ind = 0

    def labeler_setting_clicked(self):
        self.popup.set_default_text(*self.labeler.hyperparams)
        self.popup.show()
        
    def update(self, ind):    
        self.ind = ind

    def update_sticks(self, isChecked_list):
        phases = {}
        for idx, checked in enumerate(isChecked_list):
            if checked:
                sticks = np.zeros((len(self.labeler.phases[idx].peaks), 2))
                name = self.labeler.phases[idx].name

                for j, peak in enumerate(self.labeler.phases[idx].peaks):
                    sticks[j, 0] = peak.q
                    sticks[j, 1] = peak.I
                    
                phases[name] = sticks 
        self.stripeview.plot_cifs(phases)
        
    def update_pd_tab(self, tab_num):
        if tab_num == 1:
            phase_dict = self.model.get_dict_for_phase_diagram()
            self.phase_diagram_view.plot(phase_dict)
            self.phase_diagram_list.show(list(phase_dict))

    def update_pd_plot(self, mask):
        phase_dict = self.model.get_dict_for_phase_diagram()
        self.phase_diagram_view.plot(phase_dict, mask)

    def label_button_clicked(self):
        result = self.labeler.fit(self.stripeview.q, self.stripeview.avg_pattern)
        self.stripeview.replot_spectra(result)  
        #pass

    def save_button_clicked(self):
        filename = self.model.current_filename
        d = np.vstack((self.stripeview.q, self.stripeview.avg_pattern))
        fn, _ = QFileDialog.getSaveFileName(self, 'Save File', filename, "")
        np.save(fn, d)

    def change_ind(self, change):
        self.ind += change
    
    def add_to_phase_diagram(self, isChecked_list):
        phase_names = self.labeler.get_phase_names(isChecked_list)
        self.model.add_to_phase_diagram(phase_names)
        #self.update(self.ind)

    def update_labeler_hyperparams(self, std_noise, mean, std):
        self.labeler.set_hyperparams(std_noise, mean, std) 

    @property
    def ind(self):
        return self._ind

    @ind.setter
    def ind(self, new_ind):
        if new_ind >= self.model.size: 
            new_ind = 0 
        self.model.ind = new_ind
        self._ind = self.model.ind # let model do the cycling

        self.stripeview.avg_pattern = None # Not good
        self.stripeview.clear_figures()
        self.stripeview.plot(self.model.current_data)
        self.globalview.clear_figures()
        self.globalview.plot(self.model.dwells, self.model.tpeaks,
                             self.model.labeled_dwells, self.model.labeled_tpeaks,
                             self.model.current_dwell, self.model.current_tpeak,
                             self.model.x, self.model.y,
                             self.model.labeled_x, self.model.labeled_y,
                             self.model.current_x, self.model.current_y)
        existing_phase_ind = index_phase(self.model.phases[new_ind], self.labeler.phase_names) 
        
        if existing_phase_ind:
            self.cifview.clear()
            self.cifview.check_boxes(existing_phase_ind)

        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    w = 1920; h = 1080
    window = TopLevelWindow()
    window.resize(w, h)
    window.show()

    sys.exit(app.exec())
