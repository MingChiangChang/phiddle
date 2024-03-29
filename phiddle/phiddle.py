import os
import sys
import glob
import json
import logging

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import (
    QVBoxLayout, QWidget, QPushButton, QTabWidget, QFormLayout,
    QHBoxLayout, QLineEdit, QLabel, QFileDialog, QMenu
)
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar, )
from matplotlib.figure import Figure
import h5py
from tqdm import tqdm


from util import index_phase, minmax_norm
from model import datamodel
from stripeview import stripeview
from globalview import globalview
from labeling_engine import labeler
from cif_view import CIFView
from phase_diagram import PhaseDiagramView, PhaseDiagramList
from lattice_param_view import LatticeParamView, LatticeParamList
from popup import Popup
from cif_to_input_file import cif_to_input
from center_finder_asym import get_center_asym
from temp_profile import LaserPowerMing_Spring2024, left_right_width

class TopLevelWindow(QtWidgets.QMainWindow):

    def __init__(self,
                 h5_path=None,  # "data/AL_23F4_Bi-Ti-O_run_01_0_all_1d.h5",
                 csv_path=None):  # "/Users/ming/Desktop/Code/SARA.jl/BiTiO/cifs/sticks.csv" ):
        super().__init__()
        self._createMenuBar()
        self.logger = logging.getLogger(__name__)
        self.stripeview = stripeview()
        self.globalview = globalview()
        self.h5_path = h5_path
        self.csv_path = csv_path
        self.model = datamodel()
        self.labeler = labeler()
        self.cifview = CIFView([])
        self.popup = Popup()

        # For testing
        if h5_path is not None and csv_path is not None:
            self.model.read_h5(h5_path)
            self.ind = 0
            self.update(self.ind)
            self.labeler.read_csv(csv_path)
            self.cifview.update_cif_list(
                [phase.name for phase in self.labeler.phases])

        self.phase_diagram_view = PhaseDiagramView()
        self.phase_diagram_list = PhaseDiagramList()
        self.phase_diagram_list.save_signal.connect(
            self.phase_diagram_view.save_phase_diagram)
        pd_layout = QHBoxLayout()
        pd_layout.addWidget(self.phase_diagram_view)
        pd_layout.addWidget(self.phase_diagram_list)
        pd_layout.setStretch(0, 3)
        pd_layout.setStretch(1, 1)

        self.lattice_param_view = LatticeParamView()
        self.lattice_param_list = LatticeParamList()
        # self.lattice_param_list.save_signal.connect(
        #     self.phase_diagram_view.save_phase_diagram)
        lp_layout = QHBoxLayout()
        lp_layout.addWidget(self.lattice_param_view)
        lp_layout.addWidget(self.lattice_param_list)
        lp_layout.setStretch(0, 3)
        lp_layout.setStretch(1, 1)


        self.globalview.picked.connect(self.update)
        self.cifview.checked.connect(self.update_sticks)
        self.cifview.add.connect(self.add_to_phase_diagram)
        self.phase_diagram_list.checked_signal.connect(self.update_pd_plot)  # FIXME
        self.phase_diagram_list.dim_change_signal.connect(self.phase_diagram_view.change_dim)
        self.phase_diagram_list.axes_signal.connect(self.phase_diagram_view.change_axes)
        self.lattice_param_list.phase_selection_box.currentTextChanged.connect(self.lp_phase_changed)
        self.lattice_param_list.dim_change_signal.connect(self.lattice_param_view.change_dim)
        self.lattice_param_list.axes_signal.connect(self.lattice_param_view.change_axes)

        self.popup.set_clicked.connect(self.update_labeler_hyperparams)

        label_button = QPushButton()
        label_button.setText("Label")
        label_button.clicked.connect(self.label_button_clicked)

        fit_w_phase_button = QPushButton()
        fit_w_phase_button.setText("Fit With Phase")
        fit_w_phase_button.clicked.connect(self.fit_w_phase_button_clicked)

        label_w_phase_button = QPushButton()
        label_w_phase_button.setText("Label With Phase")
        label_w_phase_button.clicked.connect(self.label_w_phase_button_clicked)

        save_residual_button = QPushButton()
        save_residual_button.setText("Save Residual")
        save_residual_button.clicked.connect(self.save_residual)

        save_button = QPushButton()
        save_button.setText("Save")
        save_button.clicked.connect(self.save_button_clicked)

        back_button = QPushButton()
        back_button.setText("Back")
        back_button.clicked.connect(lambda: self.change_ind(-1))

        next_button = QPushButton()
        next_button.setText("Next")
        next_button.clicked.connect(lambda: self.change_ind(1))

        labeler_setting_button = QPushButton()
        labeler_setting_button.setText("Labeler Settings")
        labeler_setting_button.clicked.connect(self.labeler_setting_clicked)

        previous_label_result_button = QPushButton()
        previous_label_result_button.setText("Previous Label Result")
        previous_label_result_button.clicked.connect(
            self.previous_label_result)

        next_label_result_button = QPushButton()
        next_label_result_button.setText("Next Label Result")
        next_label_result_button.clicked.connect(self.next_label_result)

        top_button_layout = QHBoxLayout()
        top_button_layout.addWidget(labeler_setting_button)
        top_button_layout.addWidget(previous_label_result_button)
        top_button_layout.addWidget(next_label_result_button)

        bottom_button_layout = QHBoxLayout()
        bottom_button_layout.addWidget(label_button)
        bottom_button_layout.addWidget(label_w_phase_button)
        bottom_button_layout.addWidget(fit_w_phase_button)
        bottom_button_layout.addWidget(save_residual_button)
        bottom_button_layout.addWidget(save_button)
        bottom_button_layout.addWidget(back_button)
        bottom_button_layout.addWidget(next_button)

        main_fig_layout = QVBoxLayout()
        main_fig_layout.addLayout(top_button_layout)
        main_fig_layout.addWidget(self.stripeview)
        main_fig_layout.addWidget(self.globalview)
        main_fig_layout.addLayout(bottom_button_layout)

        outer_layout = QHBoxLayout()
        outer_layout.addLayout(main_fig_layout)
        outer_layout.addWidget(self.cifview)
        outer_layout.setStretch(0, 3)
        outer_layout.setStretch(1, 1)

        widget = QWidget()
        widget.setLayout(outer_layout)
        pd_widget = QWidget()
        pd_widget.setLayout(pd_layout)
        lp_widget = QWidget()
        lp_widget.setLayout(lp_layout)

        self.tabs = QTabWidget()
        self.tabs.addTab(widget, "Labeler")
        self.tabs.addTab(pd_widget, "Phase Map")
        self.tabs.addTab(lp_widget, "Lattice Param")
        self.tabs.currentChanged.connect(self.update_pd_tab)
        self.tabs.currentChanged.connect(self.update_lp_tab)
        self.setCentralWidget(self.tabs)


    def _createMenuBar(self):
        menuBar = self.menuBar()
        # menuBar.setNativeMenuBar(False)
        # Creating menus using a QMenu object
        fileMenu = QMenu(" &File", self)
        browse_data_file_act = QtGui.QAction("Browse Data File", self)
        browse_data_file_act.triggered.connect(self.browse_button_clicked)
        browse_csv_file_act = QtGui.QAction("Browse CSV Input File", self)
        browse_csv_file_act.triggered.connect(self.browse_csv_button_clicked)
        browse_cif_file_act = QtGui.QAction("Browse CIF Files", self)
        browse_cif_file_act.triggered.connect(self.browse_cif_button_clicked)

        fileMenu.addAction(browse_data_file_act)
        fileMenu.addAction(browse_cif_file_act)
        fileMenu.addAction(browse_csv_file_act)
        fileMenu.addSeparator()
        save_progress_act = QtGui.QAction("Save Progress", self)
        save_progress_act.triggered.connect(self.save_progress_clicked)
        load_progress_act = QtGui.QAction("Load Progress", self)
        load_progress_act.triggered.connect(self.load_progress_clicked)
        fileMenu.addAction(save_progress_act)
        fileMenu.addAction(load_progress_act)



        menuBar.addMenu(fileMenu)
        # Creating menus using a title
        editMenu = menuBar.addMenu(" &Edit")
        editMenu.addAction("test")
        helpMenu = menuBar.addMenu(" &Help")
        helpMenu.addAction("test")

    def browse_button_clicked(self):
        self.h5_path, _ = QFileDialog.getOpenFileName(None, "Open h5", "", "")
        if self.h5_path.endswith("h5"):
            self.model.read_h5(self.h5_path)
            if hasattr(self.model, "cations"):
                self.phase_diagram_list.composition_dim = len(self.model.cations) # WARNING: difficult to sync
                for i, cation in enumerate(self.model.cations):
                    self.phase_diagram_list.comp_str[i] = cation
                self.phase_diagram_list.update_combo_boxes()
                self.lattice_param_list.update_phase_combo_box()
            self.ind = 0
            # self.update(self.ind)
        elif self.h5_path.endswith("udi"):
            self.model.read_udi(self.h5_path)
            self.ind = 0
            self.update(self.ind) # FIXME: This should tell us how many dimension is allowed

    def browse_csv_button_clicked(self):
        self.csv_path, _ = QFileDialog.getOpenFileName(
                                 None, "Open csv", "", "CSV Files (*.csv)"
                                 )
        if self.csv_path.endswith("csv"):
            self.labeler.read_csv(self.csv_path)
            self.cifview.update_cif_list(
                [phase.name for phase in self.labeler.phases])

    def browse_cif_button_clicked(self):
        self.cif_paths, _ = QFileDialog.getOpenFileNames(
            None, "Open cifs", "", "")
        if not self.cif_paths:
            self.logger.error("Error: No file are chosen")
        elif not np.all([path.endswith("cif") for path in self.cif_paths]):
            self.logger.error("Error: Non cif files were included.")
        else:
            self.csv_path, _ = QFileDialog.getSaveFileName(
                None, "Store csv", "", "CSV Files (*.csv)")
            cif_to_input(self.cif_paths, self.csv_path, (10, 80))
            self.labeler.read_csv(self.csv_path)
            self.cifview.update_cif_list(
                [phase.name for phase in self.labeler.phases])

    def save_progress_clicked(self):
        self.save_fn, _ = QFileDialog.getSaveFileName(
            self, 'Save File', "", "JSON Files (*.json)")
        if self.save_fn:
            storing_ds = {}
            storing_ds["phases_diagram"] = self.model.get_dict_for_phase_diagram()
            storing_ds["phases"] = self.model.phases
            storing_ds["csv_path"] = os.path.abspath(self.csv_path)
            storing_ds["h5_path"] = os.path.abspath(self.h5_path)
            with open(self.save_fn, 'w') as f:
                json.dump(storing_ds, f)

    def load_progress_clicked(self):
        self.load_fn, _ = QFileDialog.getOpenFileName(
            None, "Open", "", "JSON Files (*.json)")
        if self.load_fn:
            with open(self.load_fn, 'r') as f:
                load_meta_data = json.load(f)

            if (os.path.isfile(load_meta_data["h5_path"])
                    and os.path.isfile(load_meta_data["csv_path"])):

                self.model.read_h5(load_meta_data["h5_path"])
                if hasattr(self.model, "cations"):
                    self.phase_diagram_list.composition_dim = len(self.model.cations) # WARNING: difficult to sync
                    for i, cation in enumerate(self.model.cations):
                        self.phase_diagram_list.comp_str[i] = cation
                    self.phase_diagram_list.update_combo_boxes()
                self.h5_path = load_meta_data["h5_path"]
                self.labeler.read_csv(load_meta_data["csv_path"])
                self.csv_path = load_meta_data["csv_path"]
                self.cifview.update_cif_list(
                    [phase.name for phase in self.labeler.phases])
                self.model.phases = load_meta_data["phases"]
                self.ind = 0
            else:
                self.logger.error(
                    f'ERROR: File in .json not found! Check if you have moved you file around')

    def labeler_setting_clicked(self):
        self.popup.set_default_text(*self.labeler.params)
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

    def update_tab(self, tab_num):
        if tab_num == 1:
            self.update_pd_tab()
        if tab_num == 2:
            self.update_lp_tab()


    def update_pd_tab(self):
        phase_dict = self.model.get_dict_for_phase_diagram()
        self.phase_diagram_view.plot(phase_dict,
                                     self.phase_diagram_list.get_current_axes())
        self.phase_diagram_list.show(list(phase_dict))

    def update_pd_plot(self, mask):
        phase_dict = self.model.get_dict_for_phase_diagram()
        self.phase_diagram_view.plot(phase_dict,
                                     self.phase_diagram_list.get_current_axes(),
                                     mask)

    def update_lp_tab(self):
        phase_dict = self.model.get_dict_for_phase_diagram()
        # self.lattice_param_view.plot_defualt(phase_dict,
        #                                 self.lattice_param_list.get_current_axes())
        phase_names = list(phase_dict)
        if "Amorphous" in phase_names:
            phase_names.remove("Amorphous")
        phase_names.insert(0, "")
        self.lattice_param_list.update_phase_combo_box(phase_names)

    def lp_phase_changed(self, phase):
        print(phase)
        if phase == "":
            return

        indicies = self.model.get_index_with_phase(phase)
        if self.model.is_refined(phase):
            refined_result_for_plot = self.model.get_refined_lp_of_phase(phase) 
        else:
            phase_ls = self.model.get_labeled_phases(indicies)

            refined_result_for_plot = []
            refined_lp = []
            for i, phases in tqdm(zip(indicies, phase_ls)):
                data = self.model[i]
                left_width, right_width = left_right_width(self.model.tpeaks[i], self.model.dwells[i])
                center = get_center_asym(data['data'], left_width, right_width)
                y, _, _ = minmax_norm(data['data'][:, center])
                q = data['q']
                res = self.labeler.fit_phases(q, y, phases)
                for cp in res.CPs:
                    refined_lp.append([cp.cl.a, cp.cl.b, cp.cl.c, cp.cl.α, cp.cl.β, cp.cl.γ])
                    if cp.name == phase:
                        refined_result_for_plot.append([cp.cl.a, cp.cl.b, cp.cl.c, cp.cl.α, cp.cl.β, cp.cl.γ])

                self.model.update_refined_lp(i, refined_lp)
            
        # FIXME: This should belong to model
        lp_dict = {}
        lp_dict["Tpeak"] = [self.model.tpeaks[i] for i in indicies]
        lp_dict["Dwell"] = [self.model.dwells[i] for i in indicies]
        cations = self.model.get_cations()
        # if len(cations) > 0:
        for j, cation in enumerate(cations):
            lp_dict[cation] = [self.model.fractions[i][j] for i in indicies]
        lp_dict["refined_lps"] = np.array(refined_result_for_plot)
        self.lattice_param_view.plot(lp_dict, self.lattice_param_list.get_current_axes())
            

    def update_lp_plot(self, mask):
        phase_dict = self.model.get_dict_for_phase_diagram()
        self.lattice_param_view.plot(phase_dict,
                                     self.phase_diagram_list.get_current_axes(),
                                     mask)



    def label_button_clicked(self):
        self.labeler.fit(self.stripeview.avg_q, self.stripeview.avg_pattern)
        self.stripeview.plot_n_store_label_result_w_spectra(
            1, self.labeler.t[0], self.labeler.results[0], self.labeler.bg)

    def fit_w_phase_button_clicked(self):
        selected_phase_names = self.cifview.get_checked_phase_names()
        if selected_phase_names:
            self.labeler.fit_phases(self.stripeview.avg_q,
                                    self.stripeview.avg_pattern,
                                    selected_phase_names)
            self.stripeview.plot_n_store_label_result_w_spectra(
                1, self.labeler.t[0], self.labeler.results[0], self.labeler.bg)

    def label_w_phase_button_clicked(self):
        selected_phase_names = self.cifview.get_checked_phase_names()
        if selected_phase_names:
            self.labeler.fit(self.stripeview.avg_q,
                                    self.stripeview.avg_pattern,
                                    selected_phase_names)
            self.stripeview.plot_n_store_label_result_w_spectra(
                1, self.labeler.t[0], self.labeler.results[0], self.labeler.bg)


    def save_residual(self):
        if not self.labeler.has_labeled:
            self.save_button_clicked()
            return
        filename = self.model.current_filename
        d = np.vstack((self.stripeview.avg_q, self.labeler.residual))
        fn, _ = QFileDialog.getSaveFileName(
            self, 'Save Residual File', filename, "")
        if fn.endswith('xy'):
            np.savetxt(fn, d)
        else:
            np.save(fn, d)

    def save_button_clicked(self):
        filename = self.model.current_filename
        d = np.vstack((self.stripeview.avg_q, self.stripeview.avg_pattern))
        fn, _ = QFileDialog.getSaveFileName(self, 'Save File', filename, "")

        if self.labeler.has_labeled:
            datadict = self.labeler.get_dict_for_storing()
            with open(fn + '.json', 'w') as f:
                json.dump(datadict, f)
        else:
            if fn.endswith('xy'):
                np.savetxt(fn, d)
            else:
                np.save(fn, d)

    def change_ind(self, change):
        self.ind += change

    def add_to_phase_diagram(self, isChecked_list):
        phase_names = self.labeler.get_phase_names(isChecked_list)
        self.model.add_to_phase_diagram(phase_names)
        # self.update(self.ind)

    def update_labeler_hyperparams(
            self,
            std_noise,
            mean,
            std,
            max_phase,
            expand_degree,
            background_length,
            max_iter,
            optimize_mode,
            background_option):
        self.labeler.set_hyperparams(
            std_noise,
            mean,
            std,
            max_phase,
            expand_degree,
            background_length,
            max_iter,
            optimize_mode,
            background_option)

    def next_label_result(self):
        if self.labeler.has_labeled:
            ind, confidence, result, fractions, bg = self.labeler.next_label_result()
            self.stripeview.plot_n_store_label_result_w_spectra(ind, confidence, result, bg)
            print("############## Output ################")
            print("")
            print(f"{ind}th most probable result")
            print(result.CPs)
            print("")
            print("Probability: {confidence}")
            print("Fractions:")
            for i, xi  in enumerate(fractions):
                print(f"    {result.CPs[i].name}: {xi}")
            print("")
            print("#####################################")

            

    def previous_label_result(self):
        if self.labeler.has_labeled:
            ind, confidence, result, fractions, bg = self.labeler.previous_label_result()
            self.stripeview.plot_n_store_label_result_w_spectra(ind, confidence, result, bg)
            print("############## Output ################")
            print("")
            print(f"{ind}th most probable result")
            print(result.CPs)
            print("")
            print(f"Probability: {confidence}")
            print("Fractions:")
            for i, xi  in enumerate(fractions):
                print(f"    {result.CPs[i].name}: {xi}")
            print("")
            print("#####################################")

    @property
    def ind(self):
        return self._ind

    @ind.setter
    def ind(self, new_ind):
        if new_ind >= self.model.size:
            new_ind = 0
        elif new_ind < 0:
            new_ind = self.model.size - 1
        self.model.ind = new_ind
        self._ind = self.model.ind  # let model do the cycling
        self.labeler.has_labeled = False

        self.stripeview.avg_pattern = None  # Not good
        self.stripeview.clear_figures()
        self.stripeview.plot_new_data(
            self.model.current_data,
            self.model.current_xx)
        self.globalview.clear_figures()
        self.globalview.plot(
            self.model.dwells,
            self.model.tpeaks,
            self.model.x,
            self.model.y,
            self.model.labeled,
            self.model.current,
            )
        try:
            existing_phase_ind = index_phase(
                self.model.phases[new_ind], self.labeler.phase_names)

            if existing_phase_ind:
                self.cifview.clear()
                self.cifview.check_boxes(existing_phase_ind)
        except AttributeError:
            pass


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    w = 1920
    h = 1080
    window = TopLevelWindow()
    window.resize(w, h)
    window.show()

    sys.exit(app.exec())
