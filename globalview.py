import os
import sys

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt6.QtCore import pyqtSignal


class globalview(FigureCanvasQTAgg):

    picked = pyqtSignal(int)

    def __init__(self, parent=None):
        fig = Figure()
        super(globalview, self).__init__(fig)
        self.setParent = parent

        gs = fig.add_gridspec(1, 2)
        self.condition_map = fig.add_subplot(gs[0, 0], picker=1)
        self.wafer_map = fig.add_subplot(gs[0, 1], picker=1)

        self.cid1 = self.mpl_connect('pick_event', self.on_pick)

    def on_pick(self, event):
        indexes = event.ind  # Indexes of the data point (array).
        self.picked.emit(indexes[0])

    def clear_figures(self):
        self.condition_map.clear()
        self.wafer_map.clear()


    def plot(self, dwells, tpeaks, current_dwell, current_tpeak,
                  x, y, current_x, current_y): # This shows that this should be separated ..
        self.clear_figures()

        self.condition_map.scatter(dwells, tpeaks, color='b', picker=True)
        self.condition_map.scatter(current_dwell, current_tpeak, color='r', picker=True)
        self.condition_map.set_xscale("log")
        self.condition_map.set_xlabel("Dwell (us)")
        self.condition_map.set_ylabel("Peak temperature ($^oC$)")

        self.wafer_map.scatter(x, y, color='b', picker=True)
        self.wafer_map.scatter(current_x, current_y, color='r')
        self.wafer_map.set_xlabel("x (mm)")
        self.wafer_map.set_ylabel("y (mm)")

        self.draw()
