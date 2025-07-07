import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from PyQt6.QtCore import pyqtSignal


class globalview(FigureCanvasQTAgg):
    """
    A class that represents a global view for displaying condition and wafer maps with matplotlib in a PyQt application.
    The class inherits from FigureCanvasQTAgg to enable plotting inside a PyQt widget.
    
    Signal:
        picked (pyqtSignal): Signal emitted when a data point is picked on the plot
    """

    picked = pyqtSignal(int)

    def __init__(self):
        """
        Initializes the GlobalView widget with two subplots: a condition map and a wafer map.
        """
        self.fig = Figure()
        super(globalview, self).__init__(self.fig)
        plt.tight_layout()

        gs = self.fig.add_gridspec(1, 2)
        self.condition_map = self.fig.add_subplot(gs[0, 0])  # , picker=1)
        self.wafer_map = self.fig.add_subplot(gs[0, 1])  # , picker=1)
        self.fig.subplots_adjust(right=0.9)
        self.cbar_ax = self.fig.add_axes([0.92, 0.15, 0.02, 0.7])

        self.cid1 = self.mpl_connect('pick_event', self.on_pick)
        self.cmap = 'inferno'

    def on_pick(self, event):
        """
        Handles pick events on the plot. Emits the signal with the index of the picked data point.
        
        Args:
            event: 
               matplotlib event object containing information about the pick.
        """
        indexes = event.ind  # Indexes of the data point (array).
        if event.artist == self.location_artist or event.artist == self.condition_artist:
            self.picked.emit(indexes[0])

    def clear_figures(self):
        """ Clears all the current plots """
        self.condition_map.clear()
        self.wafer_map.clear()
        self.cbar_ax.clear()

    # def _plot_condition(self, ax, dwells, tpeaks, labeled, current,
    #                     xlabel="Dwell (us)", ylabel="Peak temperature ($^oC$)",
    #                     xscale="log"):
    #     unlabeled = np.logical_not(labeled)
    #     self.clear_figures()

    #     self.condition_artist = self.condition_map.scatter(
    #         dwells, tpeaks, color='b', s=8, picker=True)
    #     self.condition_map.scatter(dwells[labeled], tpeaks[labeled], s=8,
    #                                color='g', picker=0.)
    #     self.condition_map.scatter(dwells[unlabeled], tpeaks[unlabeled], s=8,
    #                                color='b', picker=0.)

    #     self.condition_map.scatter(
    #         dwells[current],
    #         tpeaks[current],
    #         s=8,
    #         color='r',
    #         picker=0.)
    #     self.condition_map.set_xscale(xscale)
    #     self.condition_map.set_xlabel(xlabel)
    #     self.condition_map.set_ylabel(ylabel)

    #     dwell_size = np.log(dwells)
    #     dwell_size -= np.min(dwell_size)
    #     dwell_size /= np.max(dwell_size)

    #     self.location_artist = self.wafer_map.scatter(
    #         x, y, c=tpeaks, s=8, picker=True, alpha=0.)
    #     im = self.wafer_map.scatter(x[labeled], y[labeled], 
    #                            s=20*dwell_size[labeled]+8,
    #                            c=tpeaks[labeled], picker=0., cmap=self.cmap,
    #                            vmin=np.min(tpeaks), vmax=np.max(tpeaks))

    #     self.wafer_map.scatter(x[unlabeled], y[unlabeled], s=20*dwell_size[unlabeled]+8,
    #                            c=tpeaks[unlabeled],
    #                            marker='x', picker=0., cmap=self.cmap)

    #     self.wafer_map.scatter(x[current], y[current], s=50,
    #                            color="r",
    #                            picker=0.,
    #                            )
    #     self.wafer_map.set_xlabel("x (mm)")
    #     self.wafer_map.set_ylabel("y (mm)")
    #     # self.fig.subplots_adjust(right=0.9)
    #     # self.cbar_ax = self.fig.add_axes([0.92, 0.15, 0.02, 0.7])
    #     self.fig.colorbar(im, cax=self.cbar_ax)
    #     self.cbar_ax.set_label("$^oC$")

    #     self.draw()


    def plot(self, dwells, tpeaks, x, y, labeled, current): 
        """
        Plots the condition map and wafer map using the provided data.
        It distinguishes between labeled and unlabeled points,
        and updates the appearance of the plots accordingly.
        
        Args:
            dwells: np.ndarray
                Array of dwell times.
            tpeaks: np.ndarray
                Array of temperature peaks.
            x: np.ndarray
               Array of x-coordinates for the wafer map.
            y: np.ndarray
                Array of y-coordinates for the wafer map.
            labeled: np.ndarray
                Boolean array indicating labeled data points.
            current: np.ndarray
                Array of indices of current data points.
        """
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
            x, y, c=tpeaks, s=8, picker=True, alpha=0.)
        im = self.wafer_map.scatter(x[labeled], y[labeled], 
                               s=20*dwell_size[labeled]+8,
                               c=tpeaks[labeled], picker=0., cmap=self.cmap,
                               vmin=np.min(tpeaks), vmax=np.max(tpeaks))

        self.wafer_map.scatter(x[unlabeled], y[unlabeled], s=20*dwell_size[unlabeled]+8,
                               c=tpeaks[unlabeled],
                               marker='x', picker=0., cmap=self.cmap)

        self.wafer_map.scatter(x[current], y[current], s=50,
                               color="r",
                               picker=0.,
                               )
        self.wafer_map.set_xlabel("x (mm)")
        self.wafer_map.set_ylabel("y (mm)")
        # self.fig.subplots_adjust(right=0.9)
        # self.cbar_ax = self.fig.add_axes([0.92, 0.15, 0.02, 0.7])
        self.fig.colorbar(im, cax=self.cbar_ax)
        self.cbar_ax.set_label("$^oC$")

        self.draw()

    def get_condition_info(self):
        """
        Retrieves information about the current condition map axes,
        including axis limits and labels.
        
        Returns:
            dict:
               A dictionary containing the xmin, xmax, ymin, ymax values and axis labels.
        """
        info = {}
        xmin, xmax = self.condition_map.get_xlim()
        ymin, ymax = self.condition_map.get_ylim()
        info["xmin"] = xmin
        info["xmax"] = xmax
        info["ymin"] = ymin
        info["ymax"] = ymax
        info["xlabel"] = self.condition_map.get_xlabel()
        info["ylabel"] = self.condition_map.get_ylabel()

        return info

    def set_condition_params(self, xmin:float, xmax: float, xlabel:str,
                             ymin:float, ymax:float, ylabel:str):
        """
        Sets the parameters for the condition map axes, including axis limits and labels.
        
        Args:
            xmin: float
               Minimum value for the x-axis.
            xmax: float
                Maximum value for the x-axis.
            xlabel: str
                Label for the x-axis.
            ymin: float
                Minimum value for the y-axis.
            ymax: float
                Maximum value for the y-axis.
            ylabel: str
                Label for the y-axis.
        """

        self.condition_map.set_xlim(xmin, xmax)
        self.condition_map.set_ylim(ymin, ymax)
        self.condition_map.set_xlabel(xlabel)
        self.condition_map.set_ylabel(ylabel)
        self.draw()
