import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
      QVBoxLayout, QWidget, QPushButton, QTabWidget, QFormLayout,
      QCheckBox
      )
from matplotlib.backends.backend_qtagg import (        FigureCanvasQTAgg,
        NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure


class CIFView(QWidget):


    checked = pyqtSignal(list)
    add = pyqtSignal(list)

    def __init__(self, cif_list, parent=None):

        super(CIFView, self).__init__(parent)

        layout = QVBoxLayout()

        self.widget_ls = []
        button = QPushButton()
        button.setText("Add to phase diagram")
        button.clicked.connect(self.add_to_phase_diagram)
        for cif in cif_list:
            checkbox = QCheckBox(cif)
            checkbox.clicked.connect(self.update_stick_pattern)
            self.widget_ls.append(checkbox)
            layout.addWidget(checkbox) 

        layout.addWidget(button)
        self.setLayout(layout)


    def update_stick_pattern(self):
        self.checked.emit([checkbox.isChecked() for checkbox in self.widget_ls])


    def add_to_phase_diagram(self):
        self.add.emit([checkbox.isChecked() for checkbox in self.widget_ls])


    def clear(self):
        for checkbox in self.widget_ls:
            checkbox.setChecked(False)


    def check_boxes(self, ind_ls):
        for ind in ind_ls:
            self.widget_ls[ind].setChecked(True)
        self.update_stick_pattern()





