import os
import sys
import traceback
import json
import logging
import argparse
import time

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import (
    QVBoxLayout, QWidget, QPushButton, QTabWidget, QSizePolicy,
    QHBoxLayout, QFileDialog, QMenu, QMessageBox, QSlider, QSpacerItem
)
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar, )
from tqdm import tqdm


from util import index_phase, minmax_norm, remove_back_slash, __version__, __date__ 
from model import datamodel
from stripeview import stripeview
from globalview import globalview
# from labeling_engine import labeler
from label_api import labeler
from cif_view import CIFView
from phase_diagram import PhaseDiagramView, PhaseDiagramList
from lattice_param_view import LatticeParamView, LatticeParamList
from popup import Popup, SetupPopup
from cif_to_input_file import cif_to_input

# TODO: Better error handling, mostly for loading file
# TODO: Stripeview should update once the new center is Set
# FIXME: Remove from phase diagram does not remove the "labeled" label in model 
# TODO: The sliding bar should start at where the center is
# TODO: Add relative x position in stripe view figure title
# TODO: Add spawn thread when trying to label or update model
# TODO: Button to update all labels to currently selected temperature profile
# TODO: Need a global replot button so when some settings got updated, 
# TODO: A setting class that is linked to a setting window
# TODO: Setup shell commend in windows


class TopLevelWindow(QtWidgets.QMainWindow):
    """ The main window of phiddle """

    def __init__(self,
                 h5_path=None,
                 csv_path=None):
        super().__init__()
        self.setWindowTitle(f"Phiddle v{__version__} {__date__}") 
        self._createMenuBar()
        self.logger = logging.getLogger(__name__)
        self.stripeview = stripeview()
        self.globalview = globalview()
        self.load_fn = ""
        self._h5_path = None
        self.h5_path = h5_path
        self.csv_path = csv_path
        self.model = datamodel()
        # TODO: Should be instantiated once a csv is provided
        #       and lives as long as the phase list stays the same
        self.labeler = labeler() 

        self.cifview = CIFView([])
        self.center_slider = QSlider(orientation=QtCore.Qt.Orientation.Horizontal)
        self.center_slider.valueChanged.connect(self.user_moved_slider)
        self.popup = Popup()
        self.setuppopup = SetupPopup()

        # spacer = QSpacerItem(20, 40, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        if h5_path:
            self.h5_path = h5_path
            self.model.read_h5(h5_path)
            self.ind = 0
            self._update(self.ind)
        if csv_path:
            time.sleep(3)
            self.labeler.read_csv(csv_path)
            self.cifview.update_cif_list(self.labeler.phase_names)
            self.popup.default_phase_options = ["", *self.labeler.get_phase_names()]

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
        self.lattice_param_list.save_signal.connect(self.save_lp_diagram)
        self.lattice_param_list.save_lp_signal.connect(self.save_lp)
        lp_layout = QHBoxLayout()
        lp_layout.addWidget(self.lattice_param_view)
        lp_layout.addWidget(self.lattice_param_list)
        lp_layout.setStretch(0, 3)
        lp_layout.setStretch(1, 1)

        self.globalview.picked.connect(self._update)
        self.cifview.checked.connect(self.update_sticks)
        self.cifview.add.connect(self.add_to_phase_diagram)
        self.cifview.remove.connect(self.remove_from_phase_diagram)
        self.stripeview.heatmap_release.connect(self.update_checked_cif_list)

        self.phase_diagram_list.checked_signal.connect(self.update_pd_plot)  # FIXME
        self.phase_diagram_list.dim_change_signal.connect(self.phase_diagram_view.change_dim)
        self.phase_diagram_list.axes_signal.connect(self.phase_diagram_view.change_axes)
        self.lattice_param_list.phase_selection_box.currentTextChanged.connect(self.lp_phase_changed)
        self.lattice_param_list.dim_change_signal.connect(self.lattice_param_view.change_dim)
        self.lattice_param_list.axes_signal.connect(self.lattice_param_view.change_axes)

        self.popup.set_clicked.connect(self.update_params)
        self.setuppopup.signal_to_heatmap_view.connect(self.stripeview.set_heatmap_params)
        self.setuppopup.signal_to_condition_view.connect(self.update_params)
        self.setuppopup.signal_to_qrange.connect(self.stripeview.set_qrange)

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

        set_center_button = QPushButton()
        set_center_button.setText("Set Center") 
        set_center_button.clicked.connect(self.set_center)

        labeler_setting_button = QPushButton()
        labeler_setting_button.setText("Labeler Settings")
        labeler_setting_button.clicked.connect(self.labeler_setting_clicked)

        setup_setting_button = QPushButton()
        setup_setting_button.setText("Plot Settings")
        setup_setting_button.clicked.connect(self.plot_setting_clicked)

        previous_label_result_button = QPushButton()
        previous_label_result_button.setText("Previous Label Result")
        previous_label_result_button.clicked.connect(
            self.previous_label_result)

        next_label_result_button = QPushButton()
        next_label_result_button.setText("Next Label Result")
        next_label_result_button.clicked.connect(self.next_label_result)

        top_button_layout = QHBoxLayout()
        # top_button_layout.addSpacerItem(spacer)
        top_button_layout.addWidget(self.center_slider, stretch=1)
        top_button_layout.addWidget(set_center_button, stretch=1)
        top_button_layout.addWidget(labeler_setting_button, stretch=1)
        top_button_layout.addWidget(setup_setting_button, stretch=1)
        top_button_layout.addWidget(previous_label_result_button, stretch=1)
        top_button_layout.addWidget(next_label_result_button, stretch=1)

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
        self.tabs.currentChanged.connect(self.update_tab)
        self.setCentralWidget(self.tabs)


    def _createMenuBar(self):
        """ Create menubar """
        menuBar = self.menuBar()
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
        """ Browse h5 file paths. Load if a h5 exists """
        self.h5_path, _ = QFileDialog.getOpenFileName(None, "Open h5", "", "")
        if self.h5_path.endswith("h5"):
            self.model.read_h5(self.h5_path)
            if hasattr(self.model, "cations"):
                self.phase_diagram_list.composition_dim = len(self.model.cations)
                self.lattice_param_list.composition_dim = len(self.model.cations)
                for i, cation in enumerate(self.model.cations):
                    self.phase_diagram_list.comp_str[i] = cation
                    self.lattice_param_list.comp_str[i] = cation
            else:
                self.phase_diagram_list.composition_dim = 0
                self.lattice_param_list.composition_dim  = 0

            self.phase_diagram_list.update_combo_boxes()
            self.lattice_param_list.update_axis_combo_boxes()
            self.ind = 0
            # self._update(self.ind)
        elif self.h5_path.endswith("udi"):
            self.model.read_udi(self.h5_path)
            self.ind = 0
        else:
            print("File format not supported")


    def browse_csv_button_clicked(self):
        """ Browse csv. Load if proper csv is selected. """
        self.csv_path, _ = QFileDialog.getOpenFileName(
                                 None, "Open csv", "", "CSV Files (*.csv)"
                                 )
        if self.csv_path.endswith("csv"):
            # FIXME: Repeating actions when update csv, probably should seperate into a function
            # But not sure how to sync things, maybe play with signals
            self.labeler.read_csv(self.csv_path)
            self.cifview.update_cif_list(self.labeler.phase_names)
            self.popup.default_phase_options = ["", *self.labeler.get_phase_names()]
                # [phase.name for phase in self.labeler.phases])

    def browse_cif_button_clicked(self):
        """
        Browse cif file(s) and simulate their XRD and store information in csv file
        Prompt popup windows will show up
        """
        self.cif_paths, _ = QFileDialog.getOpenFileNames(
            None, "Open cifs", "", "")
        if not self.cif_paths:
            self.logger.error("Error: No file are chosen")
        elif not np.all([path.endswith("cif") for path in self.cif_paths]):
            self.logger.error("Error: Non cif files were included.")
        else:
            self.csv_path, _ = QFileDialog.getSaveFileName(
                None, "Store csv", "", "CSV Files (*.csv)")
            cif_to_input(self.cif_paths, self.csv_path, (10, 60)) # TODO: add range to settings
            self.labeler.read_csv(self.csv_path)
            self.cifview.update_cif_list(self.labeler.phase_names)
            self.popup.default_phase_options = ["", *self.labeler.get_phase_names()]
            #    [phase.name for phase in self.labeler.phases])

    def save_progress_clicked(self):
        """
        Save the current label process in json and can be later loaded by opening with
        the -j flag or "Load Progress" option in the GUI
        """
        self.save_fn, _ = QFileDialog.getSaveFileName(
            self, 'Save File', "", "JSON Files (*.json)")

        if self.save_fn:
            storing_ds = {}
            # storing_ds["phases_diagram"] = self.model.get_dict_for_phase_diagram()
            storing_ds["full_phase_diagram"] = self.model.labeldata.serialize_data() #seralize_data()
            storing_ds["center_idx"] = self.model.get_center_idx()

            all_phases = self.model.get_all_phases()
            for phase in all_phases:
                storing_ds[phase] = self.model.get_dict_for_lp_plot(phase) # FIXME: Load this as well
            storing_ds["phases"] = self.model.get_phases()
            storing_ds["csv_path"] = os.path.abspath(self.csv_path) if self.csv_path else ""
            storing_ds["h5_path"] = os.path.abspath(self.h5_path) if self.h5_path else ""

            with open(self.save_fn, 'w') as f:
                json.dump(storing_ds, f)


    def load_progress_clicked(self):
        """ Load the label process stored in json """
        # TODO: All load/browse button, if run successfully, should be updating all tabs
        load_fn, _ = QFileDialog.getOpenFileName( None, "Open", "", "JSON Files (*.json)")

        if not load_fn:
            self.logger.error(f'ERROR: File in .json not found! Check if you have moved you file around')
            return

        self.load(load_fn)

    def load(self, load_fn=None):
        """ load progress from json. Everything will be setup """
        if load_fn is None:
            return
        self.load_fn = load_fn
        with open(self.load_fn, 'r') as f:
            meta_data = json.load(f)

        if (os.path.isfile(meta_data["h5_path"])
                and os.path.isfile(meta_data["csv_path"])):

            self.model.read_h5(meta_data["h5_path"])
            if hasattr(self.model, "cations"):
                self.phase_diagram_list.composition_dim = len(self.model.cations)
                self.lattice_param_list.composition_dim = len(self.model.cations)
                for i, cation in enumerate(self.model.cations):
                    self.phase_diagram_list.comp_str[i] = cation
                    self.lattice_param_list.comp_str[i] = cation
            else:
                self.phase_diagram_list.composition_dim = 0
                self.lattice_param_list.composition_dim  = 0

            self.phase_diagram_list.update_combo_boxes()
            self.lattice_param_list.update_axis_combo_boxes()

            self.h5_path = meta_data["h5_path"]
            self.csv_path = meta_data["csv_path"]
            self.labeler.read_csv(meta_data["csv_path"])
            self.cifview.update_cif_list(self.labeler.phase_names)
            self.popup.default_phase_options = ["", *self.labeler.get_phase_names()]
            self.model.update_phases(meta_data["phases"])

            if "full_phase_diagram" in meta_data:
                self.model.labeldata.load_stored_label_data(meta_data["full_phase_diagram"])
            else:
                self.model.clear_label_data()
            if "center_idx" in meta_data:
                self.model.load_center_idx(meta_data["center_idx"])
            self.ind = 0


    def labeler_setting_clicked(self):
        """ show the setting popup window """
        self.popup.set_default_text(*self.labeler.params)
        self.popup.show()


    def plot_setting_clicked(self):
        """ show the plot setting popup window """
        heatmap_info = self.stripeview.get_heatmap_fig_info()
        qmin, qmax = self.stripeview.get_qrange()
        condition_info = self.globalview.get_condition_info()
        self.setuppopup.set_default_text(heatmap_info["title"],
                                         heatmap_info["xmin"],
                                         heatmap_info["xmax"],
                                         heatmap_info["xlabel"],
                                         heatmap_info["ymin"],
                                         heatmap_info["ymax"],
                                         heatmap_info["ylabel"],
                                         condition_info["xmin"],
                                         condition_info["xmax"],
                                         condition_info["xlabel"],
                                         condition_info["ymin"],
                                         condition_info["ymax"],
                                         condition_info["ylabel"],
                                         qmin, qmax)
        self.setuppopup.show()

    def _update(self, ind):
        """ change the index of the stripe showing in the GUI """
        self.ind = ind

    def update_sticks(self, isChecked_list):
        """
        Scan the check box for phases and plot the stick patterns for those that are checked
        """
        phases = {}
        for idx, checked in enumerate(isChecked_list):
            if checked:
                peaks = self.labeler.get_peaks_at(idx)
                sticks = np.zeros((len(peaks), 2))
                name = self.labeler.phase_names[idx]

                for j, peak in enumerate(peaks):
                    sticks[j, 0] = peak["q"]
                    sticks[j, 1] = peak["I"]

                phases[name] = sticks
        self.stripeview.plot_cifs(phases)

    def update_tab(self, tab_num):
        """ Changing tabs shown in phiddle  """
        if tab_num == 1:
            self.update_pd_plot()
        if tab_num == 2:
            self.update_lp_tab()

    def update_pd_plot(self, phase_list=None):
        """ update the phase diagram plot and list before showing """
        # TODO: have a check box for full or center phase ploting
        phase_dict = self.model.get_dict_for_phase_diagram()
        phase_dict_full = self.model.labeldata.get_dict_for_phase_diagram()
        cations = self.model.get_cations()

        for phase in phase_dict_full:
            if phase not in phase_dict:
                phase_dict[phase] = {}
                phase_dict[phase]['Dwell'] = []
                phase_dict[phase]['Tpeak'] = []
                for cation in cations:
                    phase_dict[phase][cation] = []
            phase_dict[phase]['Dwell'] += phase_dict_full[phase]['Dwell']
            phase_dict[phase]['Tpeak'] += phase_dict_full[phase]['Tpeak']
            for cation in cations:
                phase_dict[phase][cation] += phase_dict_full[phase][cation]

        self.phase_diagram_list._show(list(phase_dict))
        # FIXME: Bug after changing h5s
        if phase_list is None:
            phase_list = self.phase_diagram_list.get_checked_phase_names()
        plot_convex_hull = self.phase_diagram_list.convex_hull
        self.phase_diagram_view.plot(phase_dict,
                                     self.phase_diagram_list.get_current_axes(),
                                     phase_list, plot_convex_hull)

    def update_lp_tab(self):
        """ update the lattice paramter view and list before showing """
        phase_dict = self.model.get_dict_for_phase_diagram()
        # self.lattice_param_view.plot_defualt(phase_dict,
        #                                 self.lattice_param_list.get_current_axes())
        phase_names = list(phase_dict)
        if "Amorphous" in phase_names:
            phase_names.remove("Amorphous")
        if "Melt" in phase_names:
            phase_names.remove("Melt")
        phase_names.insert(0, "")
        self.lattice_param_list.update_phase_combo_box(phase_names)

    def lp_phase_changed(self, phase):
        """ 
        When the phase is selected in lattice parameter view,
        The top window will do the work to collect data from the model and send
        it to the lattice parameter view to plot

        Args:
            phase: str = phase name to be plotted
        """
        if phase == "":
            return

        indices = self.model.get_index_with_phase(phase)
        if self.model.is_refined(phase):
            refined_result_for_plot = self.model.get_refined_lp_of_phase(phase) 
        else:
            phase_ls = self.model.get_labeled_phases(indices)

            refined_result_for_plot = []
            for ind, phases in tqdm(zip(indices, phase_ls)):
                data_dict = self.model.get_lps_update_dict()
                data = self.model[ind]
                center = data["center_idx"]
                if center is None:
                    _, _, center = self.model.get_temp_profile_at(ind)
                y, _, _ = minmax_norm(data['data'][:, int(center)])

                res, uncer = self.labeler.fit_phases(data['q'], y, phases)

                for j, cp in enumerate(res["CPs"]):
                    data_dict['refined_lps'].append(
                            [cp["cl"]["a"], cp["cl"]["b"], cp["cl"]["c"],
                             cp["cl"]["α"], cp["cl"]["β"], cp["cl"]["γ"]]
                            )
                    data_dict['refined_lps_uncer'].append(uncer[j*8:j*8+6].tolist())
                    data_dict['act'].append(cp["act"])
                    data_dict['act_uncer'].append(uncer[j*8+6])
                    data_dict['width'].append(cp["σ"])
                    data_dict['width_uncer'].append(uncer[j*8+7])
                    if cp["name"] == phase:
                        refined_result_for_plot.append(
                            [cp["cl"]["a"], cp["cl"]["b"], cp["cl"]["c"],
                             cp["cl"]["α"], cp["cl"]["β"], cp["cl"]["γ"]]
                             )
                
                self.model.update_ind(ind, data_dict)
            
        lp_dict = self.model.get_dict_for_lp_plot(phase)
        lp_dict['refined_lps'] = np.array(refined_result_for_plot)
        self.lattice_param_view.plot(lp_dict, self.lattice_param_list.get_current_axes())
            
    def save_lp_diagram(self):
        """ Save the lattice paramter diagram figure """
        fn, _ = QFileDialog.getSaveFileName(
            self, 'Save Lattice Parameter Diagram', "", "")
        self.lattice_param_view.save_digram(fn)


    def save_lp(self, phase_name):
        """ 
        Save the lattice paramters 
        Args:
            phase_name: str = name of the phase that will be saved
        """
        fn, _ = QFileDialog.getSaveFileName(
            self, 'Save Lattice Parameter', remove_back_slash(phase_name), "")
        self.model.save_lp(fn, phase_name)



    def update_lp_plot(self):
        """
        Update the lattice paramter plot
        """
        phase_dict = self.model.get_dict_for_phase_diagram()
        self.lattice_param_view.plot(
                phase_dict,
                self.phase_diagram_list.get_current_axes()
            )
        # mask)


    def label_button_clicked(self):
        """
        Label the XRD pattern with all phases taken into considerations 
        Will show the fitted pattern in stripeview
        """
        self.labeler.label(self.stripeview.current_q, self.stripeview.avg_pattern)
        self.stripeview.plot_n_store_label_result_w_spectra(
            1, self.labeler.probs[0], self.labeler.fit_result, self.labeler.bg)

    def fit_w_phase_button_clicked(self):
        """
        Fit the XRD pattern with user-selected sets of phases
        Will show the fitted pattern in stripeview
        """
        selected_phase_names = self.cifview.get_checked_phase_names()
        if selected_phase_names:
            self.labeler.fit_phases(self.stripeview.current_q,
                                    self.stripeview.avg_pattern,
                                    selected_phase_names)
            self.stripeview.plot_n_store_label_result_w_spectra(
                1, self.labeler.probs[0], self.labeler.fit_result, self.labeler.bg)

    def label_w_phase_button_clicked(self):
        """
        Label the XRD pattern with user-selected sets of phases
        Will show the fitted pattern in stripeview
        """
        selected_phase_names = self.cifview.get_checked_phase_names()
        if selected_phase_names:
            self.labeler.label(self.stripeview.current_q,
                                    self.stripeview.avg_pattern,
                                    selected_phase_names)
            self.stripeview.plot_n_store_label_result_w_spectra(
                1, self.labeler.probs[0], self.labeler.fit_result, self.labeler.bg)


    def save_residual(self):
        """ Save the residual patterns of the fit  """
        if not self.labeler.has_labeled:
            self.save_button_clicked()
            return
        filename = self.model.current_filename
        d = np.vstack((self.stripeview.current_q, self.labeler.residual))
        fn, _ = QFileDialog.getSaveFileName(
            self, 'Save Residual File', filename, "")
        if fn.endswith('xy'):
            np.savetxt(fn, d)
        else:
            np.save(fn, d)

    def save_button_clicked(self):
        """ 
        Save the patterns, and the fit data into file format
        designated by the user input extension
        """
        filename = self.stripeview.get_file_name() # self.model.current_filename
        d = np.vstack((self.stripeview.current_q,
                       self.stripeview.avg_pattern))
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
        """ Use to jump serveral indices """
        self.ind += change

    def add_to_phase_diagram(self, phase_names):
        """ 
        Attach the currently selected phases to the current pattern
        (correspond to specific condition and xrd index)

        Args:
            phase_names: List[str]
               list of phase name to be added 
        """
        x_indices = self.stripeview.get_selected_frames()
        if int(np.round(self.model.current_center)) in x_indices:
            self.model.add_to_phase_diagram(phase_names)
        temps = self.stripeview.selected_temperature

        for t, x in zip(temps, x_indices):
            self.model.labeldata.update(t,
                                  self.model.current_dwell,
                                  self.model.current_composition,
                                  phase_names,
                                  self.model.current_ind, 
                                  x)

        labeled_indices = self.model.get_current_labeled_indices()
        if labeled_indices:
            self.stripeview.plot_label_progress(labeled_indices)


    def remove_from_phase_diagram(self):
        """ 
        Remove the labels currently attach to the current pattern
        (correspond to specific condition and xrd index)
        """
        x_indices = self.stripeview.get_selected_frames()

        if self.model.current_center in x_indices:
            self.model.remove_from_phase_diagram()

        for x in x_indices:
            self.model.labeldata.remove(self.model.current_ind, x)

        self.stripeview.replot_heatmap()
        labeled_indices = self.model.get_current_labeled_indices()
        if labeled_indices:
            self.stripeview.plot_label_progress(labeled_indices)


    def update_params(self, std_noise, mean, std, max_phase,
                      expand_degree, background_length, max_iter,
                      optimize_mode, background_option, default_phase_option,
                      year):
        """ Update parameters stored in the labeler """
        
        self.labeler.set_hyperparams(std_noise, mean, std, max_phase,
                                     expand_degree, background_length,
                                     max_iter, optimize_mode, background_option,
                                     default_phase_option)
        self.model.set_temp_profile_params_by_year(year)
        xaxis, temp_profile_func, _ = self.model.get_current_temp_profile()
        self.stripeview.replot_w_new_center(xaxis, temp_profile_func)
         
    def update_checked_cif_list(self, x_min_ind, x_max_ind):
        """
        Collect all of the labeled phases in the selected region
        check the check boxes
        """
        phases = self.model.get_current_phases_bw_x_range(x_min_ind, x_max_ind+1)
        if phases: 
            self.cifview.set_checked_phase_names(phases)


    def next_label_result(self):
        """ Show the next (less probably) label result """
        if self.labeler.has_labeled:
            ind, confidence, result, fractions, bg = self.labeler.next_label_result()
            self.stripeview.plot_n_store_label_result_w_spectra(
                       ind, confidence, result, bg
                    )
            self.labeler.print_current_CPs()

    def previous_label_result(self):
        """ Show the previous (more probably) label result """
        if self.labeler.has_labeled:
            ind, confidence, result, fractions, bg =\
                       self.labeler.previous_label_result()
            self.stripeview.plot_n_store_label_result_w_spectra(
                    ind, confidence, result, bg
                )
            self.labeler.print_current_CPs()


    @property
    def h5_path(self):
        return self._h5_path


    @h5_path.setter
    def h5_path(self, val):
        self._h5_path = val
        if val is not None:
            self.setWindowTitle(f"Phiddle v{__version__} {__date__} {os.path.basename(val)}") 


    @property
    def ind(self):
        return self._ind

    @ind.setter
    def ind(self, new_ind):
        """ 
        Include all the update that the GUI should do when the index of the stripe
        is changed

        1. update indices in the model
        2. Reset cahced avg pattern 
        3. replot temperature profile, heatmap and XRD pattern
        4. reset center_slider
        5. Update global view to properly show the currently chosen point
        6. Update the check box to show currently labeled phases (if any) 

        """
        # All the update and setups that has to be done when switching samples is in there
        if new_ind >= self.model.size:
            new_ind = 0
        elif new_ind < 0:
            new_ind = self.model.size - 1
        self.model.ind = new_ind
        self._ind = self.model.ind  # let model do the cycling
        self.labeler.has_labeled = False

        self.stripeview.avg_pattern = None  # Not good
        self.stripeview.clear_figures()

        if np.ndim(np.squeeze(self.model.current_data.data)) > 1:
            xaxis, temp_profile_func, center = self.model.get_current_temp_profile()
        else:
            xaxis, temp_profile_func, center = [0], lambda x: x, 0

        xlabel = self.model.get_stripe_xlabel()
        self.stripeview.plot_new_data(
            self.model.current_data,
            xaxis=xaxis, #self.model.current_xx,
            temp_profile_func = temp_profile_func, xlabel=xlabel, xx = self.model.current_xx)

        labeled_indices = self.model.get_current_labeled_indices()
        if labeled_indices:
            self.stripeview.plot_label_progress(labeled_indices)

        self.center_slider.setRange(0, len(xaxis)-1)
        self.center_slider.setValue(int(center * (len(xaxis) / len(self.model.current_xx))))

        self.globalview.clear_figures()
        self.globalview.plot(
            self.model.df['Dwell'],
            self.model.df['Tpeak'],
            self.model.df['x'],
            self.model.df['y'],
            self.model.labeled,
            self.model.current,
            )
        try:
            existing_phase_ind = index_phase(
                self.model.df['phases'][new_ind], self.labeler.phase_names)

            if existing_phase_ind:
                self.cifview.clear()
                self.cifview.check_boxes(existing_phase_ind)
        except AttributeError:
            pass

    def user_moved_slider(self, value):
        """ Make vertical bar heatmap move with the slider """
        self.stripeview.slider_moveto(value) 
        self.stripeview.plot_xrd()
        # self.stripeview.update_temp_profile(value)

    def set_center(self):
        """ Set the center of current heatmap and store in the model  """
        c = self.stripeview.transform_real_x_to_data_x(self.center_slider.value())
        self.model.set_current_center(c)
        xaxis = self.model.get_current_xaxis()
        self.stripeview.replot_w_new_center(xaxis)
        self.model.update_temp_profile_for_stored_labels() 

        labeled_indices = self.model.get_current_labeled_indices()
        if labeled_indices:
            self.stripeview.plot_label_progress(labeled_indices)

    def closeEvent(self, event):
        """ Ask for confirmation before closing """
        confirmation = QMessageBox.question(self, "Confirmation",
              "Do you want to save current progress before closing?",
              QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel)

        if confirmation == QMessageBox.StandardButton.No:
            event.accept()  # Close the app
        elif confirmation == QMessageBox.StandardButton.Yes:
            self.save_progress_clicked()
            event.accept()
        else:
            event.ignore()  #

    def keyPressEvent(self, event):
        """ 
        Handles key press event.

        'd' and right arrow key will move 1 right of the XRD pattern 
        'a' and left arrow key will move 1 left of the XRD pattern 

        'D' moves to the next stripe
        'A' moves to the previous stripe
        """
        if event.key() == QtCore.Qt.Key.Key_Left or event.text() == 'a' :
            self.stripeview.move(-1); return 
        if event.key() == QtCore.Qt.Key.Key_Right or event.text() == 'd':
            self.stripeview.move(1); return
        if  event.text() == 'A':
            self.change_ind(-1); return
        if  event.text() == 'D':
            self.change_ind(1); return

def error_handler(etype, value, tb):
    """ 
    the catch all error handler. This is the temperary solution for now.
    Should implement proper error handling in the future
    """
    error_msg = ''.join(traceback.format_exception(etype, value, tb))
    print(error_msg)

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-j', '--json', type=str, required=False, help="JSON progress file to load upon start up")
    parser.add_argument('-h5', '--h5', type=str, required=False, default="", help="Preload h5 files")
    parser.add_argument('-c', '--csv', type=str, required=False, default="", help="Preload csv files")

    return parser.parse_args()


if __name__ == "__main__":
    sys.stdin.reconfigure(encoding='utf-8')
    sys.stdout.reconfigure(encoding='utf-8')
    sys.excepthook = error_handler # This should be able to catch all exceptions

    app = QtWidgets.QApplication(sys.argv)

    w = 1920
    h = 1080

    args = parse_args()
    window = TopLevelWindow(args.h5, args.csv)
    
    window.resize(w, h)
    window.show()

    if args.json:
        time.sleep(3)
        window.load(args.json)

    sys.exit(app.exec())
