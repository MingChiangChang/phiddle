import os
import sys
from copy import deepcopy

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backend_bases import MouseButton
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines

from util import (minmax_norm, minmax_denorm, COLORS, find_first_larger,
                  find_first_smaller, two_lorentz)
from pyPhaseLabel import evaluate_obj
from temp_profile import LaserPowerMing_Spring2024, left_right_width
from center_finder_asym import get_center_asym


class stripeview(FigureCanvasQTAgg):

    def __init__(self, parent=None):
        fig = Figure()
        super(stripeview, self).__init__(fig)
        self.setParent = parent

        gs = fig.add_gridspec(3, 2)
        self.heatmap = fig.add_subplot(gs[:-1, 0])
        self.temp_profile = fig.add_subplot(gs[-1, 0])
        self.spectra = fig.add_subplot(gs[:2, 1])
        self.stick_patterns = fig.add_subplot(gs[-1, 1])

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

    @property
    def x(self):
        return np.array([self.LeftX,
                         self.LeftX,
                         self.RightX,
                         self.RightX,
                         self.LeftX], dtype=int)

    @property
    def temp_y(self):
        return np.array([self.temp_bottomY,
                         self.temp_topY,
                         self.temp_topY,
                         self.temp_bottomY,
                         self.temp_bottomY], dtype=int)


    @property
    def y(self):
        return np.array([self.bottomY,
                         self.topY,
                         self.topY,
                         self.bottomY,
                         self.bottomY], dtype=int)

    @property
    def spectra_box_x(self):
        return np.array([self.spectra_left_x,
                         self.spectra_left_x,
                         self.spectra_right_x,
                         self.spectra_right_x,
                         self.spectra_left_x])

    def onclick(self, event):

        if event.inaxes in [self.heatmap]:
            self.LeftX = int(event.xdata)
            self.bottomY = 0  # int(event.ydata)
            self.topY = self.q.shape[0]
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


    def onrelease(self, event):

        if ((event.inaxes in [self.spectra])  # Detect moveless right clicks 
                and (event.button is MouseButton.RIGHT)
                and (event.x == self._clicked_x)):
            self.spectra.clear()
            self.spectra_left_x = np.min(self.q)
            self.spectra_right_x = np.max(self.q)
            (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x, self.spectra_box_y, color='r')

            if self.fit_result is not None:
                self.plot_label_result_w_spectra()
            else:
                (self.avgplot, ) = self.spectra.plot(self.avg_q, # FIXME: Limited ranges, not the desired behavior
                                             minmax_norm(self.avg_pattern)[0],
                                             linewidth=2, color='k', label="XRD")
            self.spectra.legend(fontsize=7, loc="upper right")

            self.spectra.set_xlim((self.q[0], self.q[-1]))
            self.stick_patterns.set_xlim((self.q[0], self.q[-1]))

            self.avg_q = deepcopy(self.q)
            if self.x_min_ind == self.x_max_ind:
                self.avg_pattern = self.data[:, self.x_min_ind]
            else:
                self.avg_pattern = np.mean(self.data[:, self.x_min_ind:self.x_max_ind], axis=1)
            self.draw()
            return

        if event.inaxes in [self.heatmap]:
            self.RightX = int(event.xdata)
            self.bottomY = 0  # int(event.ydata)
            self.topY = self.q.shape[0]
        elif event.inaxes in [self.spectra]:
            self.topY = event.xdata
            self.spectra_right_x = event.xdata
            # self.topY = int(event.ydata)

        if event.inaxes in [self.heatmap, self.spectra]:
            self.selection_box.set_xdata(self.x)
            self.temp_selection_box.set_xdata(self.x)
            self.selection_box.set_ydata(self.y)

            self.aspan = self.heatmap.axvspan(
                self.LeftX, self.RightX,
                0, 1,
                color='b',
                alpha=0.25,
            )

            self.moving = False
            if self.LeftX > self.RightX:
                self.LeftX, self.RightX = self.RightX, self.LeftX
            if self.bottomY > self.topY:
                self.bottomY, self.topY = self.topY, self.bottomY
            if self.LeftX == self.RightX:
                self.t = self.get_temperature(self.LeftX)
                self.heatmap.set_title(self.get_title(self.t))
            else:
                self.t_left = self.get_temperature(self.LeftX)
                self.t_right = self.get_temperature(self.RightX)
                self.heatmap.set_title(self.get_title(self.t_left, self.t_right))
                                       
            # Let this be here for now
            q_min_ind = find_first_larger(self.q, self.bottomY)
            q_max_ind = find_first_smaller(self.q, self.topY)

            self.x_min_ind = self.transfrom_x_to_data_idx(find_first_larger(self.xaxis, self.LeftX))
            self.x_max_ind = self.transfrom_x_to_data_idx(find_first_larger(self.xaxis, self.RightX))

            
            if self.LeftX == self.RightX:
                self.avg_pattern = self.data[q_min_ind:q_max_ind, self.x_min_ind]
            else:
                self.avg_pattern = np.mean(
                    self.data[q_min_ind:q_max_ind, self.x_min_ind:self.x_max_ind], axis=1)
            self.avg_q = self.q[q_min_ind:q_max_ind]
            pattern_to_plot, _min, _max = minmax_norm(self.avg_pattern) 

            self.spectra.clear()
            (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x, self.spectra_box_y, color='r')

            self.spectra.plot(self.avg_q, pattern_to_plot, color='k', linewidth=2, label="XRD") 

            if self.fit_result is not None:
                self.plot_scaled_results(_min, _max, q_min_ind, q_max_ind)
            
            if event.button is MouseButton.LEFT:
                self.spectra.set_xlim((self.q[0], self.q[-1]))
                self.stick_patterns.set_xlim((self.q[0], self.q[-1]))
            elif event.button is MouseButton.RIGHT:
                self.spectra.set_xlim((self.avg_q[0], self.avg_q[-1]))
                self.stick_patterns.set_xlim((self.avg_q[0], self.avg_q[-1]))

            self.spectra.set_ylim((-0.3, 1.1))
            self.spectra.set_ylabel("Avg intensity (a.u.)")
            self.spectra.legend(fontsize=7, loc="upper right")
            self.draw()

    def onmotion(self, event):
        if not self.moving:
            return
        if event.inaxes is None:
            return
        if event.button not in [MouseButton.RIGHT, MouseButton.LEFT]:
            return

        if event.inaxes == self.heatmap:
            self.RightX = int(event.xdata)

            self.selection_box.set_xdata(self.x)
            self.temp_selection_box.set_xdata(self.x)
            self.selection_box.set_ydata(self.y)
            self.draw()
        elif event.inaxes == self.spectra:
            self.topY = event.xdata
            self.spectra_right_x = event.xdata
            self.selection_box.set_xdata(self.x)
            # self.temp_selection_box.set_xdata(self.x)
            self.selection_box.set_ydata(self.y)
            self.spectra_select_box.set_xdata(self.spectra_box_x)
            self.draw()

    def clear_figures(self):
        self.heatmap.clear()
        self.spectra.clear()
        self.temp_profile.clear()

        self.LeftX = 0
        self.bottomY = -100
        self.RightX = 0
        self.topY = 100

    def plot_scaled_results(self, _min, _max, q_min_ind, q_max_ind):
        self.bg_in_range = deepcopy(self.bg)
        fit = evaluate_obj(self.fit_result.CPs, self.fitted_q)
        fit += self.bg_in_range
        fit = self.rescale(fit, _min, _max)
        self.bg_in_range = self.rescale(self.bg_in_range, _min, _max)
        
        phase_name = []
        self.spectra.plot(self.fitted_q, fit, label="Fitted")
        self.spectra.plot(self.fitted_q, self.bg_in_range, label="background")
        
        for cp in self.fit_result.CPs:
            tmp = evaluate_obj(cp, self.fitted_q)
            tmp *= self._max / _max
            self.spectra.plot(self.fitted_q, tmp, label=cp.name)
            phase_name.append(cp.name)
        
        self.spectra.set_title(f"No. {self.ind} " + ("_".join(phase_name)) + f" {self.confidence:.4f}")
        self.spectra.legend(fontsize=7, loc="upper right")
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        res = self.rescale(self.fitted_pattern, _min, _max)
        res = res - fit# [q_min_ind:q_max_ind]
        self.spectra.plot(self.fitted_q, (res - .2), label="residual", c='grey')


    def plot_new_data(self, data, xx=None, stick_patterns=None):
        """ This will include initializing some attributes like self.q """
        self.clear_figures()
        self.q = data['q']
        self.data = data['data']
        self.cond = data['cond']
        self.fit_result = None
        self.xx = xx
        self.tpeak, self.dwell = self.get_tpeak_dwell_from_cond(self.cond)
        self.t_left = self.tpeak
        self.temp_topY = self.tpeak
        self.left_width, self.right_width = left_right_width(self.dwell, self.tpeak)
        self.center = get_center_asym(self.data, self.left_width, self.right_width) # FIXME: This should be handled a layer higher

        self.comp_title = ""
        if 'fracs' in data:
            self.fracs = data['fracs']
            self.cations = data['cations']
            for cation, frac in zip(self.cations, self.fracs):
                self.comp_title += f" {cation}:{frac[0]:.3f}"

        self.avg_q = deepcopy(self.q)

        self.spectra_left_x = np.min(self.q)
        self.spectra_right_x = np.max(self.q)
        (self.spectra_select_box, ) = self.spectra.plot(
            self.spectra_box_x, self.spectra_box_y, color='r')

        self.LeftX = 0
        self.RightX = 0
        if self.xx is not None:
            xmin = self.xx[0]
            xmax = self.xx[-1]
            self.xaxis = np.arange(xmax-xmin)# (xmax-xmin)/2
            self.xaxis -= self.xaxis[int(len(self.xaxis)*self.center/len(self.xx))]
            xlabel = "Location (um)"
        else:
            xmin = 0
            xmax = self.data.shape[1]
            self.xaxis = np.arange(xmax-xmin) * 10 # 10 um per column
            self.xaxis -= self.xaxis[self.center]
            xlabel = "Index"

        (self.selection_box, ) = self.heatmap.plot(self.x, self.y, color='r')
        (self.temp_selection_box, ) = self.temp_profile.plot(self.x, self.temp_y, color='r') 

        self.temp_profile_func = two_lorentz(self.tpeak,
                                             0.,
                      self.left_width, self.right_width)

        self.heatmap.imshow(self.data,
                            extent=(self.xaxis[0], self.xaxis[-1], self.q[-1], self.q[0]),
                            aspect = (self.xaxis[-1]-self.xaxis[0])/(self.q[-1]-self.q[0]))
        self.heatmap.set_box_aspect(1)
        self.title = self.get_title(t_left=self.t_left)
        self.heatmap.set_title(self.title)
        self.heatmap.set_xticks([])
        self.heatmap.set_ylabel("q ($nm^{-1}$)")

        self.temp_profile.plot(self.xaxis, self.temp_profile_func(self.xaxis))
        self.temp_profile.set_box_aspect(1/2)
        self.temp_profile.set_xlabel(xlabel)
        self.temp_profile.set_ylabel("Tpeak ($^o$C)")
        self.temp_profile.set_xlim(np.min(self.xaxis), np.max(self.xaxis))
        self.temp_profile.set_ylim(0, self.tpeak)

        self.aspan = self.heatmap.axvspan(
            self.LeftX, self.RightX, color='k', alpha=0)

        if self.avg_pattern is None:
            self.avg_pattern = self.data[:, round(self.data.shape[1] / 2)] # FIXME: change to center of t profile
        (self.avgplot, ) = self.spectra.plot(self.avg_q,
                                             minmax_norm(self.avg_pattern)[0],
                                             linewidth=2, color='k', label="XRD")
        self.spectra.legend(fontsize=7, loc="upper right")
        self.spectra.set_xlim((self.q[0], self.q[-1]))
        self.stick_patterns.set_xlim((self.q[0], self.q[-1]))
        self.spectra.set_ylim((-0.3, 1.1))
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        self.spectra.set_ylabel("Avg intensity (a.u.)")

        if stick_patterns is not None:
            self.plot_cifs(cif_patterns)

        self.draw()

    def plot_fit_result(self):
        fit = evaluate_obj(self.fit_result.CPs, self.fitted_q)

        phase_name = []
        if np.sum(self.bg) != 0.:
            fit += self.bg
        else:
            self.bg = evaluate_obj(self.fit_result.background, self.fitted_q)
            fit += self.bg

        self.spectra.plot(self.fitted_q, fit, label="Fitted")
        self.spectra.plot(self.fitted_q, self.bg, label="background")

        for cp in self.fit_result.CPs:
            self.spectra.plot(self.fitted_q, evaluate_obj(cp, self.fitted_q), label=cp.name)
            phase_name.append(cp.name)

        self.spectra.set_title(f"No. {self.ind} " + ("_".join(phase_name)) + f" {self.confidence:.4f}")
        self.spectra.legend(fontsize=7, loc="upper right")
        self.spectra.set_xlim((self.q[0], self.q[-1]))
        self.spectra.set_ylim((-0.3, 1.1))
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        self.spectra.set_ylabel("Avg intensity (a.u.)")
        self.spectra.plot(self.fitted_q, (self.fitted_pattern - fit - .2), label="residual", c='grey')

        self.draw()

    def plot_n_store_label_result_w_spectra(self, ind, confidence, fit_result=None, bg=None):
        self.fitted_q = deepcopy(self.avg_q)
        self.fitted_pattern, self._min, self._max = minmax_norm(self.avg_pattern)
        self.ind = ind
        self.confidence = confidence
        self.fit_result = fit_result
        self.bg = bg

        self.plot_label_result_w_spectra()

    def plot_label_result_w_spectra(self):
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

    def plot_cifs(self, sticks):
        """
        stick_patterns: dictionary: {phase: [q, peak_height]}
        """
        self.stick_patterns.clear()
        proxies = []
        for i, phase in enumerate(sticks):
            stick = sticks[phase]
            qs = stick[:, 0]
            Is = stick[:, 1]

            boxes = [Rectangle((q - 0.1, 0), .2, I)
                     for q, I in zip(qs, Is)]
            pc = PatchCollection(boxes, facecolor=COLORS[i], alpha=1.)
            proxy = mlines.Line2D(
                [],
                [],
                marker="_",
                linewidth=1,
                color=COLORS[i],
                markersize=15,
                label=phase)
            proxies.append(proxy)
            self.stick_patterns.add_collection(pc)

        self.stick_patterns.plot()
        self.stick_patterns.legend(handles=proxies, fontsize=7, loc="upper right")
        self.stick_patterns.set_xlim((self.q[0], self.q[-1]))
        self.stick_patterns.set_ylim((0, 1))
        self.draw()

    def get_tpeak_dwell_from_cond(self, cond):
        s = cond.split('_')
        return float(s[3]), float(s[1])

    def get_title(self, t_left, t_right=None):

        title = f"Tpeak: {int(self.tpeak)}C Dwell: {int(self.dwell)}us" 
        if t_right is None:
            title += f" T: {int(t_left)}C"
        else:
            _left = "$_{left}$"
            _right = "$_{right}$"
            title += f" T{_left}: {int(t_left)}C T{_right}: {int(t_right)}C"

        if self.comp_title:
            title += '\n'
            title += self.comp_title

        return title

    def get_temperature(self, x_idx):
        # x_idx has different meaning depending on xx
        _x_idx = np.array(x_idx)
        if self.xx is not None:
            return self.temp_profile_func(_x_idx)
        return self.temp_profile_func(_x_idx*10) # 10 um per column


    def transfrom_x_to_data_idx(self, idx):
        return int(idx/len(self.xaxis) * self.data.shape[1])

    def rescale(self, d, _min, _max):
        d = minmax_denorm(d, self._min, self._max)
        return  (d - _min) / _max

