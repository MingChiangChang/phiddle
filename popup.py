import sys

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
        QWidget, QVBoxLayout, QGridLayout, QPushButton, QHBoxLayout, QLineEdit,
        QLabel
        )

class Popup(QWidget):

    set_clicked = pyqtSignal(float, list, list, int, int, float, int)

    def __init__(self, parent=None):

        super(Popup, self).__init__(parent)

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
                         background_length: float, max_iter: int):
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

    def set_button_clicked(self):
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

            self.set_clicked.emit(std_noise,
                                  [mean_0, mean_1, mean_2],
                                  [std_0, std_1, std_2],
                                  max_phase,
                                  expand_degree,
                                  background_length,
                                  max_iter)
            self.close()
        except ValueError:
            print("The value cannot be convert to numbers.")
        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    w = 1920; h = 1080
    window = Popup()
    window.resize(w, h)
    window.show()

    sys.exit(app.exec())
