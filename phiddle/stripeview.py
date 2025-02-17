from copy import deepcopy, copy
import requests
import sys

import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backend_bases import MouseButton
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines
from PyQt6.QtCore import pyqtSignal

from util import (minmax_norm, minmax_denorm, COLORS, find_first_larger,
                  get_continue_patches, find_first_smaller)
# from pyPhaseLabel import evaluate_obj
from setting import PlotSettings, HeatmapPlotSettings, set_plot


# TODO: This need major refactoring and abstraction
#       1. Remove repetitive code in stripeview to make clear API
#       2. Model should use emit so view controller has less thing to do
# TODO: Allow vmin vmax selection in heatmap (slide bar?)
# FIXME: Size mismatch after zoom-in and zoom-out

class stripeview(FigureCanvasQTAgg):
    """
    Class that handles plotting and interaction of XRD heatmap and
    XRD patterns

    Signals:
        heatmap_release(int, int)
            emits the minimum and maximum x indices of the selected region
    """

    heatmap_release = pyqtSignal(int, int)

    def __init__(self, parent=None):
        """ Initialize graphs """
        fig = Figure()
        super(stripeview, self).__init__(fig)
        self.setParent = parent

        gs = fig.add_gridspec(3, 2)
        self.heatmap = fig.add_subplot(gs[:-1, 0])
        self.temp_profile = fig.add_subplot(gs[-1, 0])
        self.temp_profile.set_box_aspect(1/2)
        self.spectra = fig.add_subplot(gs[:2, 1])
        self.stick_patterns = fig.add_subplot(gs[-1, 1])
        self.heatmap_settings = HeatmapPlotSettings()
        self.temp_settings = PlotSettings()
        self.temp_settings.ylabel = "Tpeak ($^o$C)"
        self.spectra_settings = PlotSettings("", "q ($nm^{-1}$)",
                                             "Avg intensity (a.u.)",
                                             0., 1.,
                                             -0.3, 1.1)
        self.stick_patterns_settings = PlotSettings()

        self.bottomY = -100
        self.topY = 100
        self.temp_bottomY = 0
        self.temp_topY = 1414

        self.spectra_box_y = np.array([1.1, -1.1, -1.1, 1.1, 1.1])
        self.spectra_left_x = 0
        self.spectra_right_x = 1024

        self.avg_pattern = None
        self.avg_q = None
        self.moving = False
        self._min = 0
        self._max = 1

        self.cid1 = self.mpl_connect("button_press_event", self.onclick)
        self.cid2 = self.mpl_connect("button_release_event", self.onrelease)
        self.cid3 = self.mpl_connect("motion_notify_event", self.onmotion)

        self.fit_result = None
        self.sticks = None

    ######## Getters ########
    @property
    def box_x(self):
        """ x positions of the selected box """
        return np.array([self.LeftX,
                         self.LeftX,
                         self.RightX,
                         self.RightX,
                         self.LeftX], dtype=int)

    @property
    def temp_y(self):
        """ y positions of the selected box in temperature plot """
        return np.array([self.temp_bottomY,
                         self.temp_topY,
                         self.temp_topY,
                         self.temp_bottomY,
                         self.temp_bottomY], dtype=int)


    @property
    def box_y(self):
        """ x positions of the selected box """
        return np.array([self.bottomY,
                         self.topY,
                         self.topY,
                         self.bottomY,
                         self.bottomY], dtype=int)

    @property
    def spectra_box_x(self):
        """ x positions of the selected box in xrd spectra plot """
        return np.array([self.spectra_left_x,
                         self.spectra_left_x,
                         self.spectra_right_x,
                         self.spectra_right_x,
                         self.spectra_left_x])


    @property
    def selected_temperature(self):
        """ returns the temperatures that are in the selection box """
        r = self.get_selected_frames()
        x_pos = list(map(self.transform_data_idx_to_x, r))
        return self.temp_profile_func(np.array(x_pos))

    ###### END Getters ######

    def get_selected_frames(self):
        """ Includes the starting and ending frames """
        return range(self.transform_x_to_data_idx(find_first_larger(self.xaxis, self.LeftX)),
                     self.transform_x_to_data_idx(find_first_larger(self.xaxis, self.RightX))+1)

    def set_heatmap_title(self):
        if self.LeftX == self.RightX: # Repeating code
            self.t = self.get_temperature(self.LeftX)
            self.heatmap.set_title(self.get_title(self.t))
        else:
            self.t_left = self.get_temperature(self.LeftX)
            self.t_right = self.get_temperature(self.RightX)
            self.heatmap.set_title(self.get_title(self.t_left, self.t_right))

    def move(self, move_idx):
        """
        Move the plotting idx by move_idx. If a span is already selected, both
        left and right of the span moves

        Args:
            move_idx: str
                Amount the plotting index changes
        """
        self.spectra_left_x = np.min(self.q)
        self.spectra_right_x = np.max(self.q)

        # Let this be here for now
        self.q_min_ind = find_first_larger(self.q, self.bottomY)
        self.q_max_ind = find_first_smaller(self.q, self.topY) + 1

        self.x_min_ind = self.transform_x_to_data_idx(find_first_larger(self.xaxis, self.LeftX)) + move_idx
        self.x_max_ind = self.transform_x_to_data_idx(find_first_larger(self.xaxis, self.RightX)) + move_idx
        self.x_min_ind, self.x_max_ind = self.check_bounds(self.x_min_ind, self.x_max_ind)

        self.LeftX = self.transform_data_idx_to_x(self.x_min_ind) + 1
        self.RightX = self.transform_data_idx_to_x(self.x_max_ind) + 1

        self.set_heatmap_title()
        
        if self.LeftX == self.RightX:
            self.avg_pattern = self.data[self.q_min_ind:self.q_max_ind, self.x_min_ind]
        else:
            self.avg_pattern = np.mean(
                self.data[self.q_min_ind:self.q_max_ind, self.x_min_ind:self.x_max_ind], axis=1)
        self.avg_q = self.q[self.q_min_ind:self.q_max_ind]
        pattern_to_plot, _min, _max = minmax_norm(self.avg_pattern) 

        try:
            self.aspan.remove()
        except BaseException:
            pass

        self.aspan = self.heatmap.axvspan(self.LeftX, self.RightX, 0, 1, color='b', alpha=0.25,)

        self.spectra.clear()
        (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x, self.spectra_box_y, color='r')
        # (self.selection_box, ) = self.heatmap.plot(self.box_x, self.box_y, color='r')

        self.selection_box.set_xdata(self.box_x)
        self.temp_selection_box.set_xdata(self.box_x)
        self.selection_box.set_ydata(self.box_y)

        self.spectra.plot(self.avg_q, pattern_to_plot, color='k', linewidth=2,
                          label="XRD") 
        self.spectra.set_xlim((self.avg_q[0], self.avg_q[-1]))

        if self.fit_result is not None:
            self.plot_scaled_results(_min, _max)

        self.spectra.set_ylim((-0.3, 1.1))
        self.spectra.set_ylabel("Avg intensity (a.u.)")
        self.spectra.legend(fontsize=7, loc="upper right")
        self.draw()
        self.heatmap_release.emit(self.x_min_ind, self.x_max_ind)


    def onclick(self, event):
        """ 
        Handles click events. First check which plot the click is on. Then update
        the x, y positions for plotting the patch
        """

        if event.inaxes in [self.heatmap]:
            self.LeftX = int(event.xdata)
            self.bottomY = 0  # int(event.ydata)
            self.topY = np.max(self.q)
            self.spectra_left_x = np.min(self.q)
            self.spectra_right_x = np.max(self.q)

            try:
                self.aspan.remove()
            except BaseException:
                pass

            self.moving = True

        elif event.inaxes in [self.spectra]:

            self._clicked_x = event.x
            self.bottomY = event.xdata
            self.spectra_left_x = event.xdata
            try:
                self.aspan.remove()
            except BaseException:
                pass

            self.moving = True


    def is_moveless_right_click(self, event):
        """ Detect moveless right click, which reset the scale of the plots"""
        return ((event.inaxes in [self.spectra])
                and (event.button is MouseButton.RIGHT)
                and (event.x == self._clicked_x))


    def onrelease(self, event):
        """
        Handles all the replotting and updates once user releases the mouse click
        Moveless right click: resets the scale of the XRD spectra plot
        Right click and drag (release in XRD spectra): zoom in the XRD spectra
        Left click and drag (release in XRD spectra): select subsection for XRD spectra
        Left click and drag (release in heatmap): Take average of the selected region and 
                                                  plot in XRD spectra
        """
        if self.is_moveless_right_click(event):
            self.spectra.clear()
            self.q_min_ind = 0
            self.q_max_ind = len(self.q) + 1
            self.spectra_left_x = np.min(self.q)
            self.spectra_right_x = np.max(self.q)
            (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x,
                                                            self.spectra_box_y,
                                                            color='r')

            self.avg_q = deepcopy(self.q)
            if self.x_min_ind == self.x_max_ind:
                self.avg_pattern = self.data[:, self.x_min_ind]
            else:
                self.avg_pattern = np.mean(self.data[:, self.x_min_ind:self.x_max_ind], axis=1)

            if self.fit_result is not None:
                self.plot_label_result_w_spectra()
            else:
                (self.avgplot, ) = self.spectra.plot(self.avg_q, 
                                             minmax_norm(self.avg_pattern)[0],
                                             linewidth=2, color='k', label="XRD")
            self.spectra.legend(fontsize=7, loc="upper right")
            self.spectra.set_ylim((-0.3, 1.1))
            self.spectra.set_xlim((self.q[0], self.q[-1]))

            self.stick_patterns.clear()
            self.plot_cifs(self.sticks, (self.q[0], self.q[-1]))
            # self.stick_patterns.set_xlim((self.q[0], self.q[-1]))
            self.draw()
            return

        if event.inaxes in [self.heatmap]:
            self.RightX = int(event.xdata)
            self.bottomY = 0  # int(event.ydata)
            self.topY = np.max(self.q)
        elif event.inaxes in [self.spectra]:
            self.topY = event.xdata
            self.spectra_right_x = event.xdata
            # self.topY = int(event.ydata)

        if event.inaxes in [self.heatmap, self.spectra]:
            self.selection_box.set_xdata(self.box_x)
            self.temp_selection_box.set_xdata(self.box_x)
            self.selection_box.set_ydata(self.box_y)

            self.aspan = self.heatmap.axvspan(self.LeftX, self.RightX, 0, 1, color='b', alpha=0.25,)

            self.moving = False
            if self.LeftX > self.RightX:
                self.LeftX, self.RightX = self.RightX, self.LeftX
            if self.bottomY > self.topY:
                self.bottomY, self.topY = self.topY, self.bottomY
            self.set_heatmap_title()
                                       
            # Let this be here for now
            self.q_min_ind = find_first_larger(self.q, self.bottomY)
            self.q_max_ind = find_first_smaller(self.q, self.topY) + 1

            self.x_min_ind = self.transform_x_to_data_idx(find_first_larger(self.xaxis, self.LeftX))
            self.x_max_ind = self.transform_x_to_data_idx(find_first_larger(self.xaxis, self.RightX))
            
            if self.LeftX == self.RightX:
                self.avg_pattern = self.data[self.q_min_ind:self.q_max_ind, self.x_min_ind]
            else:
                self.avg_pattern = np.mean(
                    self.data[self.q_min_ind:self.q_max_ind, self.x_min_ind:self.x_max_ind],
                    axis=1)

            # TODO: after label, if you right click and drag then click
            #       next label result and then do right click
            #       This causes mismatch in q and pattern
            #       Only reset the fiited_q cause the scaling to be wrong
            self.avg_q = self.q[self.q_min_ind:self.q_max_ind]
            pattern_to_plot, _min, _max = minmax_norm(self.avg_pattern) 

            self.spectra.clear()
            (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x, self.spectra_box_y, color='r')

            self.spectra.plot(self.avg_q, pattern_to_plot, color='k', linewidth=2, label="XRD") 

            if self.fit_result is not None:
                self.plot_scaled_results(_min, _max)
            
            if event.button is MouseButton.LEFT:
                self.spectra.set_xlim((self.q[0], self.q[-1]))
                self.stick_patterns.set_xlim((self.q[0], self.q[-1]))
                self.plot_cifs(self.sticks, (self.q[0], self.q[-1]))
            elif event.button is MouseButton.RIGHT:

                self.plot_cifs(self.sticks, (self.avg_q[0], self.avg_q[-1]))
                self.spectra.set_xlim((self.avg_q[0], self.avg_q[-1]))
                # if self.sticks is not None:
                #     renormalized_sticks = self.renormalize_stick_in_range(qmin=self.avg_q[0],
                #                                                           qmax=self.avg_q[-1])
                #     self.plot_cifs(renormalized_sticks)
                # self.stick_patterns.set_xlim((self.avg_q[0], self.avg_q[-1]))

            self.spectra.set_ylim((-0.3, 1.1))
            self.spectra.set_ylabel("Avg intensity (a.u.)")
            self.spectra.legend(fontsize=7, loc="upper right")
            self.draw()

            if event.inaxes in [self.heatmap]:
                self.heatmap_release.emit(self.x_min_ind, self.x_max_ind)

    def onmotion(self, event):
        """ Update the box/patch in the figure user is interactive with """
        if not self.moving:
            return
        if event.inaxes is None:
            return
        if event.button not in [MouseButton.RIGHT, MouseButton.LEFT]:
            return

        if event.inaxes == self.heatmap:
            self.RightX = int(event.xdata)
            self.selection_box.set_xdata(self.box_x)
            self.temp_selection_box.set_xdata(self.box_x)
            self.selection_box.set_ydata(self.box_y)
            self.draw()
            
        elif event.inaxes == self.spectra:
            self.topY = event.xdata
            self.spectra_right_x = event.xdata
            self.selection_box.set_xdata(self.box_x)
            # self.temp_selection_box.set_xdata(self.box_x)
            self.selection_box.set_ydata(self.box_y)
            self.spectra_select_box.set_xdata(self.spectra_box_x)
            self.draw()

    def slider_moveto(self, data_idx_value):
        """ Move the center indicator with the movement of the slider  """
        x = self.xaxis[data_idx_value]
        self.temp_selection_box.set_xdata([x])
        self.RightX = x
        self.LeftX = x
        self.selection_box.set_xdata(self.box_x)
        self.selection_box.set_ydata(self.box_y)
        self.draw()

    def clear_figures(self):
        """ Reset heatmap, spectra figure and temperature profile figure"""
        self.heatmap.clear()
        self.spectra.clear()
        self.temp_profile.clear()

        self.LeftX = 0
        self.bottomY = -100
        self.RightX = 0
        self.topY = 100

    def plot_scaled_results(self, _min, _max):#, q_min_ind, q_max_ind):
        """
        Plot the rescaled result with the provided new min and max

        Args:
            _min: float
                new minimum of the plot
            _max: float
                new maximum of the plot
        """
        phase_names = list(self.fit_result)
        fit = np.sum(np.array([self.fit_result[phase_name] for phase_name in phase_names]), 
                     axis=0)
        fit += self.bg
        phase_names.remove("background")
        self.bg_in_range = self.bg + np.array(self.fit_result["background"])

        fit = self.rescale(fit, _min, _max)
        self.bg_in_range = self.rescale(self.bg_in_range, _min, _max)
        
        self.spectra.plot(self.fitted_q, fit, label="Fitted")
        self.spectra.plot(self.fitted_q, self.bg_in_range, label="background")

        for phase in phase_names:
            self.spectra.plot(self.fitted_q,
                              np.array(self.fit_result[phase]) * self._max / _max,
                              label=phase)
        
        self.spectra.set_title(f"No. {self.ind} " + ("_".join(phase_names)) + f" {self.confidence:.4f}")
        self.spectra.legend(fontsize=7, loc="upper right")
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        res = self.rescale(self.fitted_pattern, _min, _max)
        res = res - fit
        self.spectra.plot(self.fitted_q, (res - .2), label="residual", c='grey')

    def replot_w_new_center(self, xaxis, temp_profile_func=None):
        """ 
        replot the with center provided by the passed in xaxis

        Args:
            xaxis: np.ndarray
                x axis values
            temp_profile_func: Optional[Function]
                optional temperature profile function (for update)
        """
        if temp_profile_func is not None:
            self.temp_profile_func = temp_profile_func
        self.heatmap.clear() # TODO: probably don't need to replot the heat map
        self.temp_profile.clear()
        self.xaxis = xaxis
        (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x, self.spectra_box_y, color='r')

        self.LeftX, self.RightX = 0, 0
        self.t_left = self.get_temperature(self.LeftX)

        (self.selection_box, ) = self.heatmap.plot(self.box_x, self.box_y, color='r')
        (self.temp_selection_box, ) = self.temp_profile.plot(self.box_x, self.temp_y, color='r') 


        self.heatmap.imshow(self.data, extent=(xaxis[0], xaxis[-1], self.q[-1], self.q[0]),
                            aspect = (xaxis[-1]-xaxis[0])/(self.q[-1]-self.q[0]))
        self.heatmap.set_box_aspect(1)
        self.title = self.get_title(t_left=self.t_left)
        self.heatmap.set_title(self.title)
        self.heatmap.set_xticks([])
        self.heatmap.set_ylabel("q ($nm^{-1}$)")

        self.temp_profile.plot(xaxis, self.temp_profile_func(xaxis))
        self.temp_profile.set_box_aspect(1/2)
        self.temp_profile.set_xlabel(self.xlabel)
        self.temp_profile.set_ylabel("Tpeak ($^o$C)")
        self.temp_profile.set_xlim(np.min(xaxis), np.max(xaxis))
        self.temp_profile.set_ylim(0, self.tpeak)

        self.aspan = self.heatmap.axvspan(
            self.LeftX, self.RightX, color='k', alpha=0)

        self.plot_xrd()

        self.draw()

    def replot_heatmap(self, xaxis=None, temp_profile_func=None):
        """
        Replot the heatmap with currently stored data.
        Mostly for updating what is overlayed on top of it
        """

        if temp_profile_func is not None:
            self.temp_profile_func = temp_profile_func
        if xaxis is not None:
            self.xaxis = xaxis

        self.heatmap.clear()

        (self.selection_box, ) = self.heatmap.plot(self.box_x, self.box_y, color='r')
        self.heatmap.imshow(self.data,
                            extent=(self.xaxis[0], self.xaxis[-1], self.q[-1], self.q[0]),
                            aspect = (self.xaxis[-1]-self.xaxis[0])/(self.q[-1]-self.q[0]))
        self.heatmap.set_box_aspect(1)
        self.title = self.get_title(t_left=self.t_left)
        self.heatmap.set_title(self.title)
        self.heatmap.set_xticks([])
        self.heatmap.set_ylabel("q ($nm^{-1}$)")

        self.aspan = self.heatmap.axvspan(
            self.LeftX, self.RightX, color='k', alpha=0)
        self.draw()


    def plot_label_progress(self, labeled_indices):
        """ 
        Plot a tiny green patch at the bottom of the heatmap for locations that 
        have been labeled 
        Args:
            labeled_indices: List[int]
                indices in the heatmap that is labeled
        """
        height = self.q[-10]-self.q[-1] 
        start_width_ls = get_continue_patches(labeled_indices)
        if not start_width_ls: return 
        recs = []
        for start_idx, width in start_width_ls:
            recs.append(Rectangle([float(self.transform_data_idx_to_x(start_idx-0.5)), self.q[-1]],
                          width*self.get_unit_length_x(), height))
        pc = PatchCollection(recs, facecolor='lime', edgecolor='lime')
        self.heatmap.add_collection(pc)
        # self.heatmap.scatter(1000*np.random.rand(self.data.shape[1])-500,
        #                      np.zeros(self.data.shape[1])+np.mean(self.q))
        # self.heatmap.set_xlim(self.xaxis[0], self.xaxis[-1])
        # self.heatmap.set_ylim(self.q[-1], self.q[0])
        self.draw()
        

    def plot_new_data(self, data, xaxis, temp_profile_func, xlabel,
                      xx=None, stick_patterns=None):
        """ This will include initializing some attributes like self.q """
        self.clear_figures()
        self.q = data['q']
        self.data = data['data']
        self.topY = np.max(data['q'])
        self.x = data['x']
        self.y = data['y']
        self.xlabel = xlabel
        self.fit_result = None
        self.xx = xx
        self.temp_profile_func = temp_profile_func
        self.xaxis = xaxis
        self.tpeak, self.dwell = data['Tpeak'], data['Dwell']
        self.t_left = self.tpeak
        self.temp_topY = self.tpeak

        self.comp_title = ""
        if 'fracs' in data:
            self.fracs = data['fracs']
            self.cations = data['cations']
            for cation, frac in zip(self.cations, self.fracs):
                self.comp_title += f" {cation}:{frac:.3f}"

        self.avg_q = deepcopy(self.q)

        self.spectra_left_x = np.min(self.q)
        self.spectra_right_x = np.max(self.q)
        (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x, self.spectra_box_y, color='r')

        self.LeftX, self.RightX = 0, 0

        (self.selection_box, ) = self.heatmap.plot(self.box_x, self.box_y, color='r')
        (self.temp_selection_box, ) = self.temp_profile.plot(self.box_x, self.temp_y, color='r') 

        xspan = max(1, xaxis[-1] - xaxis[0])

        std = np.std(self.data) 
        median = np.median(self.data)
        # FIXME: Should be posible to merge all calls of imshow
        self.heatmap.imshow(self.data, extent=(xaxis[0], xaxis[-1], self.q[-1], self.q[0]),
                            aspect = xspan/(self.q[-1]-self.q[0]), vmin=median-3*std, vmax=median+3*std)
        self.heatmap.set_box_aspect(1)
        self.title = self.get_title(t_left=self.t_left)
        self.heatmap.set_title(self.title)
        self.heatmap.set_xticks([])
        self.heatmap.set_ylabel("q ($nm^{-1}$)")

        self.temp_profile.plot(xaxis, temp_profile_func(xaxis))
        # self.temp_profile.set_box_aspect(1/2)
        self.temp_profile.set_xlabel(xlabel)
        self.temp_profile.set_ylabel("Tpeak ($^o$C)")
        self.temp_profile.set_xlim(np.min(xaxis), np.max(xaxis))
        self.temp_profile.set_ylim(0, self.tpeak)

        self.aspan = self.heatmap.axvspan(
            self.LeftX, self.RightX, color='k', alpha=0)

        if self.avg_pattern is None:
            self.avg_pattern = self.data[:, self.transform_x_to_data_idx(find_first_larger(xaxis, 0.))] # FIXME: change to center of t profile
        (self.avgplot, ) = self.spectra.plot(self.avg_q,
                                             minmax_norm(self.avg_pattern)[0],
                                             linewidth=2, color='k', label="XRD")
        self.spectra.legend(fontsize=7, loc="upper right")
        self.stick_patterns.set_xlim((self.q[0], self.q[-1]))
        self.spectra.set_xlim((self.q[0], self.q[-1]))
        self.spectra.set_ylim((-0.3, 1.1))
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        self.spectra.set_ylabel("Avg intensity (a.u.)")

        if stick_patterns is not None:
            self.plot_cifs(stick_patterns)

        self.draw()

    def plot_fit_result(self):
        """ Plot the fitting result store in self """
        phase_names = list(self.fit_result)
        fit = np.sum(np.array([self.fit_result[phase_name] for phase_name in phase_names]), 
                     axis=0)
        fit += self.bg
        phase_names.remove("background")
        bg = self.bg + np.array(self.fit_result["background"])

        fit = fit[self.q_min_ind:self.q_max_ind]
        bg = bg[self.q_min_ind:self.q_max_ind]

        self.spectra.plot(self.fitted_q, fit, label="Fitted")
        self.spectra.plot(self.fitted_q, bg, label="background")

        for phase in phase_names:
            self.spectra.plot(self.fitted_q,
                              self.fit_result[phase][self.q_min_ind:self.q_max_ind],
                              label=phase)

        self.spectra.set_title(f"No. {self.ind} " + ("_".join(phase_names)) + f" {self.confidence:.4f}")
        self.spectra.legend(fontsize=7, loc="upper right")
        # FIXME: spectra settings
        # self.spectra_settings.xmin = self.q[0]
        # self.spectra_settings.xmax = self.q[-1]
        # set_plot(self.spectra, self.spectra_settings)

        self.spectra.set_xlim((self.fitted_q[0], self.fitted_q[-1]))
        self.spectra.set_ylim((-0.3, 1.1))
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        self.spectra.set_ylabel("Avg intensity (a.u.)")
        self.spectra.plot(self.fitted_q,
                          (self.fitted_pattern - fit - .2),
                          label="residual", c='grey')

        self.draw()

    # FIXME: Allow "next label result" while zoomed in
    def plot_n_store_label_result_w_spectra(self, ind, confidence, fit_result, bg=None):
        """ Store in labeling result in self and plot """
        self.fitted_q = deepcopy(self.avg_q)
        self.fitted_pattern, self._min, self._max = minmax_norm(self.avg_pattern)
        self.ind = ind
        self.confidence = confidence
        self.fit_result = fit_result
        self.bg = bg

        self.plot_label_result_w_spectra()

    def plot_label_result_w_spectra(self):
        """ plot both spectra and fit result """
        self.spectra.clear()

        if self.avg_pattern is None:
            self.avg_pattern = self.data[:, round(self.data.shape[1] / 2)]
        
        (self.avgplot, ) = self.spectra.plot(self.fitted_q,
                                             self.fitted_pattern,
                                             linewidth=2,
                                             color='k',
                                             label="XRD")

        (self.spectra_select_box, ) = self.spectra.plot(
            self.spectra_box_x, self.spectra_box_y, color='r')

        self.plot_fit_result()

    def plot_cifs(self, sticks, q_range=None):
        """
        Plot the stick patterns of phases
        Args:
            stick_patterns: dict
                {phase: [q, peak_height]}
            q_range: Tuple(float, float)
                specifies minimum and maximum q range
        """
        if sticks is None:
            return
        self.stick_patterns.clear()
        self.sticks = sticks # Save for responding to user manipulations

        if q_range is not None:
            q_min, q_max = q_range
        elif self.avg_q is not None:
            q_min, q_max = self.avg_q[0], self.avg_q[-1]
        else:
            q_min, q_max = self.q[0], self.q[-1]

        proxies = []
        for i, phase in enumerate(sticks):
            stick = sticks[phase]
            mask = np.logical_and(stick[:, 0]>= q_min, stick[:, 0] <= q_max)
            qs = stick[mask, 0]
            Is = stick[mask, 1]
            if len(Is) != 0:
                Is /= np.max(Is)

            pc = self.get_patches(qs, Is, COLORS[i])
            proxy = self.get_legend_proxy(COLORS[i], phase)
            proxies.append(proxy)
            self.stick_patterns.add_collection(pc)

        self.stick_patterns.plot()
        self.stick_patterns.legend(handles=proxies, fontsize=7, loc="upper right")
        self.stick_patterns.set_xlim((q_min, q_max))
        self.stick_patterns.set_ylim((0, 1))
        self.draw()

    def get_patches(self, qs, Is, color):
        """
        Produce patchs with  
        Args:
            qs: iterator
                1D array/list q vector
            Is: iterator
                1D array/list intensity
        Return:
            PatchCollection
        """
        boxes = [Rectangle((q - 0.1, 0), .2, I)
                 for q, I in zip(qs, Is)]
        pc = PatchCollection(boxes, facecolor=color, alpha=1.)
        return pc

    def get_legend_proxy(self, color, phase_name):
        """ Use this to fake legend signs """
        return mlines.Line2D([], [], marker="_", linewidth=1, color=color,
                             markersize=15, label=phase_name)


    def get_tpeak_dwell_from_cond(self, cond):
        """ Helper function to get tpeak and dwell from conditions """
        s = cond.split('_')
        return float(s[3]), float(s[1])

    def get_title(self, t_left, t_right=None): # FIXME:Sloppy,  probably should get from Model
        """ Get the title for spectra figure """

        title = f"Tpeak: {int(self.tpeak)}C Dwell: {int(self.dwell)}us" 
        if t_right is None:
            title += f" T: {int(t_left)}C"
        else:
            _left = "$_{left}$"
            _right = "$_{right}$"
            title += f" T{_left}: {int(t_left)}C T{_right}: {int(t_right)}C"

        title += '\n'
        title += f"x:{self.x:.1f} y:{self.y:.1f}"

        if self.comp_title:
            title += ' '
            title += self.comp_title

        return title


    def get_file_name(self): # FIXME:Sloppy,  probably should get from Model
        """ Get default file name for storing spectra """

        filename = f"tau_{int(self.dwell)}us_" 

        if self.LeftX == self.RightX:
            self.t = self.get_temperature(self.LeftX)
            filename += f"T_{int(self.t)}"
        else:
            self.t_left = self.get_temperature(self.LeftX)
            self.t_right = self.get_temperature(self.RightX)
            filename += f"Tleft_{int(self.t_left)}_Tright_{int(self.t_right)}"

        if hasattr(self, 'cations'):
            for cation, frac in zip(self.cations, self.fracs):
                # filename += "_self."
                filename += f"{cation}_{frac:.3f}"
        return filename

    def get_temperature(self, x_idx):
        """ 
        Get ther temperature 
        Args:
            x_idx: int
                x index value
        """
        # x_idx has different meaning depending on xx
        _x_idx = np.array(x_idx)
        if self.xx is not None:
            return self.temp_profile_func(_x_idx)
        return self.temp_profile_func(_x_idx*10) # 10 um per column


    def transform_x_to_data_idx(self, x):
        """ 
        helper function to transform between different x axis
        matplotlib uses x values but when accessing data, we need to transform
        that into the data axis 

        args:
            x: int
                matplotlib x axis value
        return:
            int: corresponding data x axis value
        """
        data_idx = int(x/len(self.xaxis) * self.data.shape[1])
        if data_idx >= self.data.shape[1]:
            return self.data.shape[1]-1
        return max(0, data_idx)

    def transform_data_idx_to_x(self, x):
        """ 
        helper function to transform between different x axis
        args:
            x: int
                data x axis value
        return:
            int: corresponding matplotlib x axis value
        """
        return int(self.xaxis[0] + x / self.data.shape[1] * len(self.xaxis))

    def transform_real_x_to_data_x(self, x):
        return x/len(self.xaxis) * self.data.shape[1]

    def get_unit_length_x(self):
        return  len(self.xaxis) / self.data.shape[1]

    def check_bounds(self, *args):
        """
        make every element in args to be with in (0, data.shape[1]-1)

        args:
            *args: multiple element that should be put in range
        return:
            list that contains elements in range
        """

        r = []
        for a in args:
            if a >= self.data.shape[1]-1:
                r.append(self.data.shape[1]-1)
            else:
                r.append(max(0, a)) # CAUTION: we assume x axis is always >= 0
        return r

    def rescale(self, d, _min, _max):
        """ Inverse of min max normalization"""
        d = minmax_denorm(d, self._min, self._max)
        return  (d - _min) / _max

    # def renormalize_stick_in_range(self, qmin, qmax):
    #     """ 
    #     Select the sticks that lies between (qmin, qmax) and normalized to have
    #     max of 1 in intensity for each phase

    #     args:
    #         qmin, qmax: float
    #     return:
    #         dict containing

    #     """
    #     if self.sticks is None:
    #         return 

    #     rs = {} # deepcopy(self.sticks)

    #     for phase in self.sticks:
    #         sticks = rs[phase]
    #         _q = sticks[:, 0] 

    #         mask = np.logical_and(_q >= qmin, _q <= qmax)
    #         _q = _q[mask]
    #         _I = sticks[mask, 1]
    #         if len(_q) != 0:
    #             _I /= np.max(_I)

    #         rs[phase] = np.vstack((_q, _I)).T

    #     return rs
            
    def plot_xrd(self):
        """
        Basic replot; assuming nothing changes and everything is read from current state
        """
        self.spectra.clear()
        (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x, self.spectra_box_y, color='r')

        self.q_min_ind = find_first_larger(self.q, self.bottomY)
        self.q_max_ind = find_first_smaller(self.q, self.topY) + 1
        self.x_min_ind = self.transform_x_to_data_idx(find_first_larger(self.xaxis, self.LeftX))
        if self.LeftX == self.RightX:
            self.avg_pattern = self.data[self.q_min_ind:self.q_max_ind, self.x_min_ind]
        else:
            self.avg_pattern = np.mean(
                self.data[self.q_min_ind:self.q_max_ind, self.x_min_ind:self.x_max_ind], axis=1)
        pattern_to_plot, _min, _max = minmax_norm(self.avg_pattern) 
        self.avg_q = self.q[self.q_min_ind:self.q_max_ind]

        self.spectra.plot(self.avg_q, pattern_to_plot, color='k', linewidth=2, label="XRD") 

        if self.fit_result is not None:
            self.plot_scaled_results(_min, _max)#, q_min_ind, q_max_ind)
        
        self.spectra.set_xlim((self.q[0], self.q[-1]))
        self.stick_patterns.set_xlim((self.q[0], self.q[-1]))
        self.plot_cifs(self.sticks, (self.q[0], self.q[-1]))

        self.spectra.set_ylim((-0.3, 1.1))
        self.spectra.set_ylabel("Avg intensity (a.u.)")
        self.spectra.legend(fontsize=7, loc="upper right")
        self.draw()

    def get_heatmap_fig_info(self):
        """ 
        Get current heatmap figure information and return as dictionary
        info included:
            ["title", "xmin", "xmax", "ymin", "ymax", "xlabel", "ylabel"]
        """

        info = {}
        info["title"] = self.heatmap.get_title()
        xmin, xmax = self.heatmap.get_xlim()
        ymin, ymax = self.heatmap.get_ylim()
        info["xmin"] = xmin
        info["xmax"] = xmax
        info["ymin"] = ymin
        info["ymax"] = ymax
        info["xlabel"] = self.heatmap.get_xlabel()
        info["ylabel"] = self.heatmap.get_ylabel()

        return info

    def set_heatmap_params(self, title: str, xmin: float, xmax: float, xlabel: str,
                           ymin: float, ymax: float, ylabel: str):
        """
        Set heatmap parameters
        args:
            title: str
            xmin: float
            xmax: float
            xlabel: str
            ymin: float
            ymax: float
            ylabel: str
        """
        self.heatmap.set_title(title)
        self.heatmap.set_xlim(xmin, xmax)
        self.heatmap.set_ylim(ymin, ymax)
        self.heatmap.set_xlabel(xlabel)
        self.heatmap.set_ylabel(ylabel)
        self.draw()

    def set_qrange(self, qmin:float, qmax:float):
        """ 
        set q range of spectra
        args:
            qmin: float
            qmax: float
        """
        self.spectra.set_xlim(qmin, qmax)
        self.draw()

    def get_qrange(self):
        """ get q range of spectra"""
        return self.spectra.get_xlim()
