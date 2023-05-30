import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
      QVBoxLayout, QWidget, QPushButton, QTabWidget, QFormLayout,
      QCheckBox, QFileDialog
      )
from matplotlib.backends.backend_qtagg import (        FigureCanvasQTAgg,
        NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure

from util import COLORS


class PhaseDiagramView(FigureCanvasQTAgg):

    def __init__(self, parent=None, xlim=(250, 10000), ylim=(400, 1400)):

        fig = Figure()
        super(PhaseDiagramView, self).__init__(fig)
        self.setParent = parent

        self.phase_diagram = fig.add_subplot()

        self.xlim = xlim
        self.ylim = ylim

        self.draw()

    def plot(self, phase_dict, mask=None):

        phase_name_ls = np.array(list(phase_dict))
        if mask:
            phase_name_ls = phase_name_ls[mask]
        self.phase_diagram.clear()
        for idx, phase in enumerate(phase_name_ls):
            self.phase_diagram.scatter(phase_dict[phase]['dwell'],
                                       phase_dict[phase]['tpeak'],
                                       label=phase,
                                       color=COLORS[idx],
                                       alpha=0.5)

        self.phase_diagram.set_xlim(self.xlim)
        self.phase_diagram.set_xlabel("Dwell (us)")
        self.phase_diagram.set_xscale("log")
        self.phase_diagram.set_ylim(self.ylim)
        self.phase_diagram.set_ylabel("Peak Temperature ($^oC$)")

        self.phase_diagram.legend()
        self.draw()     

    def save_phase_diagram(self):
        fn, _ = QFileDialog.getSaveFileName(self, 'Save Phase Diagram', "", "")
        self.phase_diagram.figure.savefig(fn)


class PhaseDiagramList(QWidget):

    checked = pyqtSignal(list)
    save = pyqtSignal()

    def __init__(self, parent=None):

        super(PhaseDiagramList, self).__init__(parent)

        self.save_button = QPushButton()
        self.save_button.setText("Save Phase Diagram")
        self.save_button.clicked.connect(lambda: self.save.emit())

        self.outer_layout = QVBoxLayout()
        self.layout = QVBoxLayout()
        self.outer_layout.addLayout(self.layout)
        self.outer_layout.addWidget(self.save_button)
        self.setLayout(self.outer_layout)
        self.widget_ls = []

    def show(self, phases):

        for idx, phase in enumerate(phases):
            if idx >= len(self.widget_ls):
                checkbox = QCheckBox(phase)
                checkbox.clicked.connect(self.update_phase_diagram)
                self.widget_ls.append(checkbox)
                self.layout.addWidget(checkbox)
            else:
                checkbox = self.widget_ls[idx] 
                checkbox.setText(phase)
            checkbox.setChecked(True)

        if len(phases) > len(self.widget_ls):
            for i in range(len(phase),len(self.widget_ls)):
                self.layout.removeWidget(self.widget_ls[i])
                del self.widget_ls[i]


    def update_phase_diagram(self):
        self.checked.emit([checkbox.isChecked() for checkbox in self.widget_ls])
