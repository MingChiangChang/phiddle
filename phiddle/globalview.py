import os
import sys

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from PyQt6.QtCore import pyqtSignal


class globalview(FigureCanvasQTAgg):

    picked = pyqtSignal(int)

    def __init__(self, parent=None):
        self.fig = Figure()
        super(globalview, self).__init__(self.fig)
        plt.tight_layout()
        self.setParent = parent

        gs = self.fig.add_gridspec(1, 2)
        self.condition_map = self.fig.add_subplot(gs[0, 0])  # , picker=1)
        self.wafer_map = self.fig.add_subplot(gs[0, 1])  # , picker=1)

        self.cid1 = self.mpl_connect('pick_event', self.on_pick)
        self.cmap = 'inferno'

    def on_pick(self, event):
        indexes = event.ind  # Indexes of the data point (array).
        if event.artist == self.location_artist or event.artist == self.condition_artist:
            self.picked.emit(indexes[0])

    def clear_figures(self):
        self.condition_map.clear()
        self.wafer_map.clear()

    def plot(self, dwells, tpeaks, x, y,
             labeled, current): 
             # labeled_dwells, labeled_tpeaks,
             # current_dwell, current_tpeak,
             # x, y,
             # labeled_x, labeled_y,
             # current_x, current_y):  # This shows that this should be separated ..
        unlabeled = np.logical_not(labeled)
        self.clear_figures()

        self.condition_artist = self.condition_map.scatter(
            dwells, tpeaks, color='b', s=8, picker=True)
        self.condition_map.scatter(dwells[labeled], tpeaks[labeled], s=8,
                                   color='g', picker=0.)
        self.condition_map.scatter(dwells[unlabeled], tpeaks[unlabeled], s=8,
                                   color='b', picker=0.)

        self.condition_map.scatter(
            dwells[current],
            tpeaks[current],
            s=8,
            color='r',
            picker=0.)
        self.condition_map.set_xscale("log")
        self.condition_map.set_xlabel("Dwell (us)")
        self.condition_map.set_ylabel("Peak temperature ($^oC$)")
        dwell_size = np.log(dwells)
        dwell_size -= np.min(dwell_size)
        dwell_size /= np.max(dwell_size)

        self.location_artist = self.wafer_map.scatter(
            x, y, c=tpeaks, s=8, picker=True, alpha=0.,
            vmin=np.min(tpeaks), vmax=np.max(tpeaks))
        im = self.wafer_map.scatter(x[labeled], y[labeled], 
                                    s=20*dwell_size[labeled]+8,
                               c=tpeaks[labeled], picker=0., cmap=self.cmap)
        self.wafer_map.scatter(x[unlabeled], y[unlabeled], s=20*dwell_size[unlabeled]+8,
                               c=tpeaks[unlabeled],
                               marker='x', picker=0., cmap=self.cmap)

        self.wafer_map.scatter(x[current], y[current], s=50,
                               color="r",
                               picker=0.,
                               )
        self.wafer_map.set_xlabel("x (mm)")
        self.wafer_map.set_ylabel("y (mm)")
        self.fig.subplots_adjust(right=0.9)
        cbar_ax = self.fig.add_axes([0.92, 0.15, 0.02, 0.7])
        self.fig.colorbar(im, cax=cbar_ax)

        self.draw()
