import os
import sys
import glob

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


class TopLevelWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        # Temperature fixed
        self.stripeview = stripeview()
        self.globalview = globalview()
        self.model = datamodel("data/AL_23F4_Bi-Ti-O_run_01_0_all_1d.h5") 
        self.labeler = labeler("/Users/ming/Desktop/Code/SARA.jl/BiTiO/cifs/sticks.csv")
        self.cifview = CIFView([phase.name for phase in self.labeler.phases])

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

        self.globalview.picked.connect(self.update)
        self.cifview.checked.connect(self.update_sticks)
        self.cifview.add.connect(self.add_to_phase_diagram)

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
        
        browse_button_layout = QHBoxLayout()
        browse_button_layout.addWidget(browse_button)
        browse_button_layout.addWidget(browse_cif_button)


        button_layout = QHBoxLayout()
        button_layout.addWidget(label_button)
        button_layout.addWidget(save_button)
        button_layout.addWidget(back_button)
        button_layout.addWidget(next_button)

        layout = QVBoxLayout()
        layout.addLayout(browse_button_layout)
        layout.addWidget(browse_button)
        layout.addWidget(self.stripeview)
        layout.addWidget(self.globalview)
        layout.addLayout(button_layout)

        FPLayout = QHBoxLayout()

        outer_layout = QHBoxLayout()
        outer_layout.addLayout(layout)
        outer_layout.addWidget(self.cifview)
        outer_layout.setStretch(0, 3)
        outer_layout.setStretch(1, 1)

        widget = QWidget()
        widget.setLayout(outer_layout)
        #tab2 = QWidget()
        #self.tabs = QTabWidget()
        #self.tabs.addTab(widget, "Tab 1")
        #self.tabs.addTab(tab2, "Tab 2")
        self.setCentralWidget(widget)#self.tabs)

    def browse_button_clicked(self):
        self.file_name, _ = QFileDialog.getOpenFileName(None, "Open", "", "")
        print(self.file_name)
        #if self.file_name[0] != '':
        #    self.FPLineEdit.setText(self.file_name[0])
        self.model = datamodel(self.file_name)
        self.ind = 0
        self.update(self.ind)

    def browse_cif_button_clicked(self):
        #self.cif_dir_name = QFileDialog.getExistingDirectory(self, "Select Directory") 
        self.cif_csv_fn, _ = QFileDialog.getOpenFileName(None, "Open", "", "")
        self.labeler = labeler(self.cif_csv_fn)


    def update(self, ind):    
        self.ind = ind
        #self.model.ind = ind
        #self.stripeview.clear_figures()
        #self.stripeview.avg_pattern = None
        #self.stripeview.plot(self.model.current_data)
        #self.globalview.clear_figures()
        #self.globalview.plot(self.model.dwells, self.model.tpeaks,
        #                     self.model.labeled_dwells, self.model.labeled_tpeaks,
        #                     self.model.current_dwell, self.model.current_tpeak,
        #                     self.model.x, self.model.y,
        #                     self.model.labeled_x, self.model.labeled_y, 
        #                     self.model.current_x, self.model.current_y)

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
        

    def label_button_clicked(self):
        #result = self.labeler.fit(self.stripeview.q, self.stripeview.avg_pattern)
        #self.stripeview.replot_spectra(result)  
        pass

    def save_button_clicked(self):
        filename = self.model.current_filename
        d = np.vstack((self.stripeview.q, self.stripeview.avg_pattern))
        np.save(filename, d)

    def change_ind(self, change):
        self.ind += change
    
    def add_to_phase_diagram(self, isChecked_list):
        phase_names = self.labeler.get_phase_names(isChecked_list)
        self.model.add_to_phase_diagram(phase_names)
        self.update(self.ind)

    @property
    def ind(self):
        return self._ind

    @ind.setter
    def ind(self, new_ind):
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
        print(existing_phase_ind)
        if existing_phase_ind:
            self.cifview.clear()
            self.cifview.check_boxes(existing_phase_ind)

        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    w = 2280; h = 1520
    window = TopLevelWindow()
    window.resize(w, h)
    window.show()

    sys.exit(app.exec())
