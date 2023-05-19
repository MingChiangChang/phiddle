import sys

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
        QWidget, QVBoxLayout, QGridLayout, QPushButton, QHBoxLayout, QLineEdit,
        QLabel
        )

class Popup(QWidget):

    set_clicked = pyqtSignal(float, list, list)

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

        
        set_button = QPushButton()
        set_button.setText("Set")
        set_button.clicked.connect(self.set_button_clicked)
        
        cancel_button = QPushButton()
        cancel_button.setText("Cancel")
        cancel_button.clicked.connect(self.close)

        grid_layout.addWidget(set_button, 3, 3)
        grid_layout.addWidget(cancel_button, 3, 2)

        self.setLayout(grid_layout)

    def set_default_text(self, std_noise: float, mean: list, std: list):
        self.std_noise_edit.setText(str(std_noise))
        self.mean_0_edit.setText(str(mean[0]))
        self.mean_1_edit.setText(str(mean[1]))
        self.mean_2_edit.setText(str(mean[2]))
        self.std_0_edit.setText(str(std[0]))
        self.std_1_edit.setText(str(std[1]))
        self.std_2_edit.setText(str(std[2]))

    def set_button_clicked(self):
        try:
            std_noise = float(self.std_noise_edit.text())
            mean_0 = float(self.mean_0_edit.text())
            mean_1 = float(self.mean_1_edit.text())
            mean_2 = float(self.mean_2_edit.text())
            std_0 = float(self.std_0_edit.text())
            std_1 = float(self.std_1_edit.text())
            std_2 = float(self.std_2_edit.text())

            self.set_clicked.emit(std_noise,
                                  [mean_0, mean_1, mean_2],
                                  [std_0, std_1, std_2])
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
