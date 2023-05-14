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

    def plot(self, phase_dicts):

        self.phase_diagram.clear()
        for idx, phase in enumerate(phase_dicts):
            self.phase_diagram.scatter(phase_dicts[phase]['dwell'],
                                       phase_dicts[phase]['tpeak'],
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
