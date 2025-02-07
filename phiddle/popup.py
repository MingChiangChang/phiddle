import sys

from PyQt6 import QtCore, QtWidgets
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QWidget, QGridLayout, QPushButton, QLineEdit,
    QLabel, QComboBox
)
import numpy as np


class Popup(QWidget):
    """ 
    Popup window for collecting settings related to phase diagram, such as noise, mean,
    standard deviation, phase, expansion degree, and optimization options.
    
    Attributes:
        set_clicked (pyqtSignal): Signal emitted when the "Set" button is clicked.
    """
    set_clicked = pyqtSignal(float, list, list, int, int, float, int, str, str, str)

    def __init__(self, parent=None):
        """ Initializes the Popup window and sets up UI elements. """

        super(Popup, self).__init__(parent)
        # Make popup window always stay on top
        self.setWindowFlag(QtCore.Qt.WindowType.WindowStaysOnTopHint, True)

        self.std_noise_edit = QLineEdit()
        self.std_noise_prompt = QLabel()
        self.std_noise_prompt.setText("std_noise")

        grid_layout = QGridLayout()
        self.mean_prompt = QLabel()
        self.mean_prompt.setText("mean")
        self.std_prompt = QLabel()
        self.std_prompt.setText("std")
        self.mean_0_edit = QLineEdit()
        self.mean_1_edit = QLineEdit()
        self.mean_2_edit = QLineEdit()

        self.std_0_edit = QLineEdit()
        self.std_1_edit = QLineEdit()
        self.std_2_edit = QLineEdit()

        self.max_phase_prompt = QLabel()
        self.max_phase_prompt.setText("Max phase")
        self.max_phase_edit = QLineEdit()

        self.expand_degree_prompt = QLabel()
        self.expand_degree_prompt.setText("Expand degree")
        self.expand_degree_edit = QLineEdit()

        self.background_length_prompt = QLabel()
        self.background_length_prompt.setText("Background length")
        self.background_length_edit = QLineEdit()

        self.max_iter_prompt = QLabel()
        self.max_iter_prompt.setText("Max Iterations")
        self.max_iter_edit = QLineEdit()

        self.optimize_mode_prompt = QLabel()
        self.optimize_mode_prompt.setText("Optimize Mode")
        self.optimize_mode_dropdown = QComboBox()
        self.optimize_mode_options = ["Simple", "EM", "With Uncertainty"]
        self.optimize_mode_dropdown.addItems(self.optimize_mode_options)

        self.background_option_prompt = QLabel()
        self.background_option_prompt.setText("Background option")
        self.background_dropdown = QComboBox()
        self.background_options = ["MCBL", "Default", "None"]
        self.background_dropdown.addItems(self.background_options)

        self.temp_profile_prompt = QLabel()
        self.temp_profile_prompt.setText("Temp profile year")
        self.temp_profile_dropdown = QComboBox()
        self.temp_profile_year_options = ["2025", "2024", "2023", "2021"] # IDEA: These option can be made into a config file
        self.temp_profile_dropdown.addItems(self.temp_profile_year_options)

        grid_layout.addWidget(self.std_noise_prompt, 0, 0)
        grid_layout.addWidget(self.std_noise_edit, 0, 1)
        grid_layout.addWidget(self.mean_prompt, 1, 0)
        grid_layout.addWidget(self.std_prompt, 2, 0)
        grid_layout.addWidget(self.mean_0_edit, 1, 1)
        grid_layout.addWidget(self.mean_1_edit, 1, 2)
        grid_layout.addWidget(self.mean_2_edit, 1, 3)
        grid_layout.addWidget(self.std_0_edit, 2, 1)
        grid_layout.addWidget(self.std_1_edit, 2, 2)
        grid_layout.addWidget(self.std_2_edit, 2, 3)
        grid_layout.addWidget(self.max_phase_prompt, 3, 0)
        grid_layout.addWidget(self.max_phase_edit, 3, 1)
        grid_layout.addWidget(self.expand_degree_prompt, 4, 0)
        grid_layout.addWidget(self.expand_degree_edit, 4, 1)
        grid_layout.addWidget(self.background_length_prompt, 5, 0)
        grid_layout.addWidget(self.background_length_edit, 5, 1)
        grid_layout.addWidget(self.max_iter_prompt, 6, 0)
        grid_layout.addWidget(self.max_iter_edit, 6, 1)
        grid_layout.addWidget(self.optimize_mode_prompt, 3, 2)
        grid_layout.addWidget(self.optimize_mode_dropdown, 3, 3)
        grid_layout.addWidget(self.background_option_prompt, 4, 2)
        grid_layout.addWidget(self.background_dropdown, 4, 3)
        grid_layout.addWidget(self.temp_profile_prompt, 7, 0)
        grid_layout.addWidget(self.temp_profile_dropdown, 7, 1)

        set_button = QPushButton()
        set_button.setText("Set")
        set_button.clicked.connect(self.set_button_clicked)

        cancel_button = QPushButton()
        cancel_button.setText("Cancel")
        cancel_button.clicked.connect(self.close)

        grid_layout.addWidget(set_button, 7, 3)
        grid_layout.addWidget(cancel_button, 7, 2)

        self.setLayout(grid_layout)

    def set_default_text(self, std_noise: float, mean: list, std: list,
                         max_phase: int, expand_degree: int,
                         background_length: float, max_iter: int,
                         optimize_mode: str, background_option: str):
        """
        Sets default values for the input fields in the popup window.

        Args:
            std_noise (float): Standard deviation of noise.
            mean (list): List of mean values.
            std (list): List of standard deviations.
            max_phase (int): Maximum phase.
            expand_degree (int): Expansion degree.
            background_length (float): Background length.
            max_iter (int): Maximum iterations for optimization.
            optimize_mode (str): Optimization mode.
            background_option (str): Background option (e.g., "MCBL", "Default").
        """
 
        self.std_noise_edit.setText(str(std_noise))
        self.mean_0_edit.setText(str(mean[0]))
        self.mean_1_edit.setText(str(mean[1]))
        self.mean_2_edit.setText(str(mean[2]))
        self.std_0_edit.setText(str(std[0]))
        self.std_1_edit.setText(str(std[1]))
        self.std_2_edit.setText(str(std[2]))
        self.max_phase_edit.setText(str(max_phase))
        self.expand_degree_edit.setText(str(expand_degree))
        self.background_length_edit.setText(str(background_length))
        self.max_iter_edit.setText(str(max_iter))

        self.optimize_mode_dropdown.setCurrentIndex(
            self.optimize_mode_options.index(optimize_mode)
        )
        self.background_dropdown.setCurrentIndex(
            self.background_options.index(background_option)
        )

    def set_button_clicked(self):
        """
        Emits a signal with the collected input values when the 'Set' button is clicked.

        The values are emitted via the `set_clicked` signal to be handled by the parent.
        """
        try:
            std_noise = float(self.std_noise_edit.text())
            mean_0 = float(self.mean_0_edit.text())
            mean_1 = float(self.mean_1_edit.text())
            mean_2 = float(self.mean_2_edit.text())
            std_0 = float(self.std_0_edit.text())
            std_1 = float(self.std_1_edit.text())
            std_2 = float(self.std_2_edit.text())
            max_phase = int(self.max_phase_edit.text())
            expand_degree = int(self.expand_degree_edit.text())
            background_length = float(self.background_length_edit.text())
            max_iter = int(self.max_iter_edit.text())
            optimize_mode = self.optimize_mode_dropdown.currentText()
            background_option = self.background_dropdown.currentText()
            temp_year_option = self.temp_profile_dropdown.currentText()

            self.set_clicked.emit(std_noise,
                                  [mean_0, mean_1, mean_2],
                                  [std_0, std_1, std_2],
                                  max_phase,
                                  expand_degree,
                                  background_length,
                                  max_iter,
                                  optimize_mode,
                                  background_option,
                                  temp_year_option)
            self.close()
        except ValueError:
            # TODO: Add warning window (QDialog?)
            print("The value cannot be convert to numbers.")


class SetupPopup(QWidget):
    """
    Popup window for configuring settings related to heatmap, condition view, and Q range.

    Attributes:
        signal_to_heatmap_view (pyqtSignal): Signal to update heatmap view settings.
        signal_to_qrange (pyqtSignal): Signal to update Q range.
        signal_to_condition_view (pyqtSignal): Signal to update condition view.
    """

    signal_to_heatmap_view = pyqtSignal(str, float, float, str, float, float, str)
    signal_to_qrange = pyqtSignal(float, float)
    signal_to_condition_view = pyqtSignal(float, float, str, float, float, str)
     
    def __init__(self):
        """ Initializes the SetupPopup window and sets up the UI. """
        super().__init__()
        grid_layout = QGridLayout()

        self.heatmap_title_edit = QLineEdit()
        self.heatmap_title_label = QLabel()
        self.heatmap_title_label.setText("Heat map title")

        self.heatmap_xmin_edit = QLineEdit()
        self.heatmap_xmax_edit = QLineEdit()
        self.heatmap_xlim_label = QLabel()
        self.heatmap_xlim_label.setText("Heat map X limits")

        self.heatmap_ymin_edit = QLineEdit()
        self.heatmap_ymax_edit = QLineEdit()
        self.heatmap_ylim_label = QLabel()
        self.heatmap_ylim_label.setText("Heat map Y limits")

        self.heatmap_xlabel_edit = QLineEdit()
        self.heatmap_xlabel_prompt = QLabel()
        self.heatmap_xlabel_prompt.setText("Heat map X label")

        self.heatmap_ylabel_edit = QLineEdit()
        self.heatmap_ylabel_prompt = QLabel()
        self.heatmap_ylabel_prompt.setText("Heat map Y label")

        self.condition_xmin_edit = QLineEdit()
        self.condition_xmax_edit = QLineEdit()
        self.condition_xlim_label = QLabel()
        self.condition_xlim_label.setText("Condition X limits")

        self.condition_ymin_edit = QLineEdit()
        self.condition_ymax_edit = QLineEdit()
        self.condition_ylim_label = QLabel()
        self.condition_ylim_label.setText("Condition Y limits")

        self.condition_xlabel_edit = QLineEdit()
        self.condition_xlabel_prompt = QLabel()
        self.condition_xlabel_prompt.setText("Condition X label")

        self.condition_ylabel_edit = QLineEdit()
        self.condition_ylabel_prompt = QLabel()
        self.condition_ylabel_prompt.setText("Condition Y label")

        self.stripe_qmin_edit = QLineEdit()
        self.stripe_qmax_edit = QLineEdit()
        self.stripe_qlim_prompt = QLabel()
        self.stripe_qlim_prompt.setText("Q range limit(nm-1)")

        set_button = QPushButton()
        set_button.setText("Set")
        set_button.setDefault(True)
        set_button.clicked.connect(self.set_button_clicked)

        cancel_button = QPushButton()
        cancel_button.setText("Cancel")
        cancel_button.clicked.connect(self.close)

        grid_layout.addWidget(self.heatmap_title_label,     0, 0)
        grid_layout.addWidget(self.heatmap_title_edit,      0, 1)
        grid_layout.addWidget(self.heatmap_xlim_label,      1, 0)
        grid_layout.addWidget(self.heatmap_xmin_edit,       1, 1)
        grid_layout.addWidget(self.heatmap_xmax_edit,       1, 2)
        grid_layout.addWidget(self.heatmap_xlabel_prompt,   2, 0)
        grid_layout.addWidget(self.heatmap_xlabel_edit,     2, 1)

        grid_layout.addWidget(self.heatmap_ylim_label,      3, 0)
        grid_layout.addWidget(self.heatmap_ymin_edit,       3, 1)
        grid_layout.addWidget(self.heatmap_ymax_edit,       3, 2)
        grid_layout.addWidget(self.heatmap_ylabel_prompt,   4, 0)
        grid_layout.addWidget(self.heatmap_ylabel_edit,     4, 1)

        grid_layout.addWidget(self.condition_xlim_label,    5, 0)
        grid_layout.addWidget(self.condition_xmin_edit,     5, 1)
        grid_layout.addWidget(self.condition_xmax_edit,     5, 2)
        grid_layout.addWidget(self.condition_xlabel_prompt, 6, 0)
        grid_layout.addWidget(self.condition_xlabel_edit,   6, 1)

        grid_layout.addWidget(self.condition_ylim_label,    7, 0)
        grid_layout.addWidget(self.condition_ymin_edit,     7, 1)
        grid_layout.addWidget(self.condition_ymax_edit,     7, 2)
        grid_layout.addWidget(self.condition_ylabel_prompt, 8, 0)
        grid_layout.addWidget(self.condition_ylabel_edit,   8, 1)

        grid_layout.addWidget(self.stripe_qlim_prompt,      9, 0)
        grid_layout.addWidget(self.stripe_qmin_edit,        9, 1)
        grid_layout.addWidget(self.stripe_qmax_edit,        9, 2)

        grid_layout.addWidget(set_button,                  10, 1)
        grid_layout.addWidget(cancel_button,               10, 2)


        self.setLayout(grid_layout)


    def set_default_text(self, h_title, h_xmin, h_xmax, h_xlabel,
                         h_ymin, h_ymax, h_ylabel,
                         c_xmin, c_xmax, c_xlabel, c_ymin, c_ymax, c_ylabel,
                         qmin, qmax):
        """
        Sets default values for the input fields in the setup popup.
        
        Args:
            heatmap_title (str): Title for the heatmap.
            heatmap_xlim (list): List with xmin and xmax values for the heatmap.
            heatmap_ylim (list): List with ymin and ymax values for the heatmap.
            condition_xlim (list): List with xmin and xmax values for the condition.
            condition_ylim (list): List with ymin and ymax values for the condition.
            xlabel (str): Label for the x-axis.
            ylabel (str): Label for the y-axis.
        """
        self.heatmap_title_edit.setText(h_title)
        self.heatmap_xmin_edit.setText(str(np.round(h_xmin, 5)))
        self.heatmap_xmax_edit.setText(str(np.round(h_xmax, 5)))
        self.heatmap_xlabel_edit.setText(h_xlabel)
        self.heatmap_ymin_edit.setText(str(np.round(h_ymin, 5)))
        self.heatmap_ymax_edit.setText(str(np.round(h_ymax, 5)))
        self.heatmap_ylabel_edit.setText(h_ylabel)

        self.condition_xmin_edit.setText(str(np.round(c_xmin, 5)))
        self.condition_xmax_edit.setText(str(np.round(c_xmax, 5)))
        self.condition_xlabel_edit.setText(c_xlabel)
        self.condition_ymin_edit.setText(str(np.round(c_ymin, 5)))
        self.condition_ymax_edit.setText(str(np.round(c_ymax, 5)))
        self.condition_ylabel_edit.setText(c_ylabel)

        self.stripe_qmin_edit.setText(str(np.round(qmin, 5)))
        self.stripe_qmax_edit.setText(str(np.round(qmax, 5)))


    def set_button_clicked(self):
        """
        Emits a signal with the collected input values when the 'Set' button is clicked.

        The values are emitted via the appropriate signals to update heatmap,
        condition view, and Q-range settings.
        """
        # TODO: modify a setting object instead and use the object to plot everytime
        title = self.heatmap_title_edit.text()
        xmin = float(self.heatmap_xmin_edit.text())
        xmax = float(self.heatmap_xmax_edit.text())
        xlabel = self.heatmap_xlabel_edit.text()
        ymin = float(self.heatmap_ymin_edit.text())
        ymax = float(self.heatmap_ymax_edit.text())
        ylabel = self.heatmap_ylabel_edit.text()
        self.signal_to_heatmap_view.emit(title, xmin, xmax, xlabel, ymin, ymax, ylabel)

        xmin = float(self.condition_xmin_edit.text())
        xmax = float(self.condition_xmax_edit.text())
        xlabel = self.condition_xlabel_edit.text()
        ymin = float(self.condition_ymin_edit.text())
        ymax = float(self.condition_ymax_edit.text())
        ylabel = self.condition_ylabel_edit.text()
        self.signal_to_condition_view.emit(xmin, xmax, xlabel, ymin, ymax, ylabel)

        qmin = float(self.stripe_qmin_edit.text())
        qmax = float(self.stripe_qmax_edit.text())
        self.signal_to_qrange.emit(qmin, qmax)



if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    w = 1920
    h = 1080
    window = Popup()
    window.resize(w, h)
    window.show()

    sys.exit(app.exec())
