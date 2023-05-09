import os
import sys
import glob

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import (
      QVBoxLayout, QWidget, QPushButton, QTabWidget, QFormLayout,
      QHBoxLayout, QLineEdit, QLabel
      )
from matplotlib.backends.backend_qtagg import (        FigureCanvasQTAgg,
        NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
import h5py

from model import datamodel 
from stripeview import stripeview
from globalview import globalview
from labeling_engine import labeler

class TopLevelWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        # Temperature fixed
        self.stripeview = stripeview()
        self.globalview = globalview()
        self.model = datamodel("data/AL_23F4_Bi-Ti-O_run_01_0_all_1d.h5") 
        self.labeler = labeler("/Users/ming/Desktop/Code/SARA.jl/BiTiO/cifs/sticks.csv")

        self.ind = 0
        self.update()
        self.stripeview.plot(self.model[self.ind])
        self.globalview.plot(self.model.dwells, self.model.tpeaks,
                             self.model.current_dwell, self.model.current_tpeak,
                             self.model.x, self.model.y,
                             self.model.current_x, self.model.current_y)

        self.globalview.picked.connect(self.update)

        label_button = QPushButton()
        label_button.setText("Label")
        label_button.clicked.connect(self.label_button_clicked)

        save_button = QPushButton()
        save_button.setText("Save")
        save_button.clicked.connect(self.save_button_clicked)

        next_button = QPushButton()
        next_button.setText("Next")
        next_button.clicked.connect(self.next_button_clicked)

        layout = QVBoxLayout()
        layout.addWidget(self.stripeview)
        layout.addWidget(self.globalview)
        layout.addWidget(label_button)
        layout.addWidget(save_button)
        layout.addWidget(next_button)

        FPLayout = QHBoxLayout()
        self.label = QLabel()
        self.label.setObjectName("File Path")
        browseButton = QPushButton('Browse')
        browseButton.clicked.connect(self.browse_path)
        self.FPLineEdit = QLineEdit()
        self.FPLineEdit.setPlaceholderText("/path/to/file")
        FPLayout.addWidget(self.label)
        FPLayout.addWidget(self.FPLineEdit)
        FPLayout.addWidget(browseButton)


        outer_layout = QHBoxLayout()
        outer_layout.addLayout(layout)
        outer_layout.addLayout(FPLayout)

        widget = QWidget()
        widget.setLayout(outer_layout)
        #tab2 = QWidget()
        #self.tabs = QTabWidget()
        #self.tabs.addTab(widget, "Tab 1")
        #self.tabs.addTab(tab2, "Tab 2")
        self.setCentralWidget(widget)#self.tabs)

    def browse_path(self):
        self.file_name = QtWidgets.QFileDialog.getOpenFileName(None, "Open", "", "")
        if self.file_name[0] != '':
            self.FPLineEdit.setText(self.file_name[0])

    def update(self, ind):    
        self.model.ind = ind
        self.stripeview.clear_figures()
        self.stripeview.avg_pattern = None
        self.stripeview.plot(self.model.current_data)
        self.globalview.clear_figures()
        self.globalview.plot(self.model.dwells, self.model.tpeaks,
                             self.model.current_dwell, self.model.current_tpeak,
                             self.model.x, self.model.y,
                             self.model.current_x, self.model.current_y)

    def label_button_clicked(self):
        result = self.labeler.fit(self.stripeview.q, self.stripeview.avg_pattern)
        self.stripeview.replot_spectra(result)  

    def save_button_clicked(self):
        filename = self.model.current_filename
        d = np.vstack((self.stripeview.q, self.stripeview.avg_pattern))
        np.save(filename, d)

    def next_button_clicked(self):
        self.model.ind += 1
        self.stripeview.avg_pattern = None # Not good
        self.stripeview.clear_figures()
        self.stripeview.plot(self.model.current_data)
        self.globalview.clear_figures()
        self.globalview.plot(self.model.dwells, self.model.tpeaks,
                             self.model.current_dwell, self.model.current_tpeak,
                             self.model.x, self.model.y, 
                             self.model.current_x, self.model.current_y)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    w = 2280; h = 1520
    window = TopLevelWindow()
    window.resize(w, h)
    window.show()

    sys.exit(app.exec())
