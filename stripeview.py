import os
import sys
from copy import deepcopy

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines

from util import minmax_norm, COLORS, find_first_larger, find_first_smaller
from pyPhaseLabel import evaluate_obj

class stripeview(FigureCanvasQTAgg):

    def __init__(self, parent=None):
        fig = Figure()
        super(stripeview, self).__init__(fig)
        self.setParent = parent

        gs = fig.add_gridspec(3, 2)
        self.heatmap = fig.add_subplot(gs[:, 0])
        self.spectra = fig.add_subplot(gs[:2, 1])
        self.stick_patterns = fig.add_subplot(gs[-1, 1])

        self.bottomY = -100
        self.topY = 100

        self.spectra_box_y = np.array([1.1, -0.1, -0.1, 1.1, 1.1])
        self.spectra_left_x = 0
        self.spectra_right_x = 1024

        self.avg_pattern = None
        self.avg_q = None
        self.moving = False

        self.cid1 = self.mpl_connect("button_press_event", self.onclick)
        self.cid2 = self.mpl_connect("button_release_event", self.onrelease)
        self.cid3 = self.mpl_connect("motion_notify_event", self.onmotion)
        
    @property
    def x(self):
        return np.array([self.LeftX,
                         self.LeftX,
                         self.RightX,
                         self.RightX,
                         self.LeftX], dtype=int)

    @property
    def y(self):
        return np.array([self.bottomY,
                         self.topY,
                         self.topY,
                         self.bottomY,
                         self.bottomY], dtype=int)

    @property
    def spectra_box_x(self):
        return np.array([self.spectra_left_x, self.spectra_left_x, self.spectra_right_x,
                         self.spectra_right_x, self.spectra_left_x])

    def onclick(self, event):

        if event.inaxes in [self.heatmap]:
            self.LeftX = int(event.xdata)
            self.bottomY = 0#int(event.ydata)
            self.topY = self.q.shape[0] 
            self.spectra_left_x = np.min(self.q)
            self.spectra_right_x = np.max(self.q)

            try:
                self.aspan.remove()
            except:
                pass

            self.moving = True
        elif event.inaxes in [self.spectra]:
            self.bottomY = event.xdata
            self.spectra_left_x = event.xdata
            try:
                self.aspan.remove()
            except:
                pass

            self.moving = True


    def onrelease(self, event):
        if event.inaxes in [self.heatmap]:
            self.RightX = int(event.xdata)
            self.bottomY = 0#int(event.ydata)
            self.topY = self.q.shape[0]
        elif event.inaxes in [self.spectra]:
            self.topY = event.xdata
            self.spectra_right_x = event.xdata
            #self.topY = int(event.ydata)

        if event.inaxes in [self.heatmap, self.spectra]:
            self.selection_box.set_xdata(self.x)
            self.selection_box.set_ydata(self.y)

            self.aspan = self.heatmap.axvspan(
                    self.LeftX, self.RightX,
                    0,1,
                    color = 'b',
                    alpha = 0.25,
            )

            self.moving = False
            if self.LeftX > self.RightX:
                self.LeftX, self.RightX = self.RightX, self.LeftX
            if self.bottomY > self.topY:
                self.bottomY, self.topY = self.topY, self.bottomY
            # Let this be here for now
            q_min_ind = find_first_larger(self.q, self.bottomY) 
            q_max_ind = find_first_smaller(self.q, self.topY)
            if self.LeftX == self.RightX: 
                self.avg_pattern = self.data[q_min_ind:q_max_ind, self.LeftX] 
            else:
                self.avg_pattern = np.mean(self.data[q_min_ind:q_max_ind,self.LeftX:self.RightX], axis=1)
            self.avg_q = self.q[q_min_ind:q_max_ind]
            self.avg_pattern = minmax_norm(self.avg_pattern)
            self.spectra.clear()
            self.spectra.plot(self.avg_q, self.avg_pattern, color='k', linewidth=2,label="XRD")
            (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x,
                                                        self.spectra_box_y, color='r')
            self.spectra.set_xlim((self.q[0], self.q[-1]))
            self.spectra.set_ylim((-0.1, 1.1))
            self.spectra.set_ylabel("Avg intensity (a.u.)")
            self.spectra.legend(fontsize=7)
            self.draw()

            
    def onmotion(self, event):
        if not self.moving:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return

        if event.inaxes == self.heatmap:
            self.RightX = int(event.xdata)

            self.selection_box.set_xdata(self.x)
            self.selection_box.set_ydata(self.y)
            self.draw()
        elif event.inaxes == self.spectra:
            self.topY = event.xdata
            self.spectra_right_x = event.xdata
            self.selection_box.set_xdata(self.x)
            self.selection_box.set_ydata(self.y)
            self.spectra_select_box.set_xdata(self.spectra_box_x)
            self.draw()

    def clear_figures(self):
        self.heatmap.clear()
        self.spectra.clear()

        self.LeftX = 0
        self.bottomY = -100
        self.RightX = 0
        self.topY = 100

    def plot_new_data(self, data, stick_patterns=None, fit_result=None):
        """ This will include initializing some attributes like self.q """
        self.clear_figures()
        self.q = data['q']
        self.data = data['data']
        self.cond = data['cond']
        self.LeftX = round(self.data.shape[1]/2) 
        self.RightX = round(self.data.shape[1]/2) 
        self.avg_q = deepcopy(self.q)

        self.spectra_left_x = np.min(self.q)
        self.spectra_right_x = np.max(self.q)
        (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x,
                                                        self.spectra_box_y, color='r')

        (self.selection_box, ) = self.heatmap.plot(self.x, self.y, color='r')
        self.heatmap.imshow(self.data,
                       extent=(0, self.data.shape[1], self.q[-1], self.q[0]),
                       aspect=self.data.shape[1]/(self.q[-1]-self.q[0]))
        self.heatmap.set_title(self.cond)
        self.heatmap.set_xlabel("index")
        self.heatmap.set_ylabel("q ($nm^{-1}$)")       


        self.aspan = self.heatmap.axvspan(self.LeftX, self.RightX, color='k', alpha=0)

        if self.avg_pattern is None:
            self.avg_pattern = minmax_norm(self.data[:,round(self.data.shape[1]/2)])
        (self.avgplot, ) = self.spectra.plot(self.avg_q, self.avg_pattern, linewidth=2,
                                             color='k', label="XRD")
        self.spectra.legend(fontsize=7)
        self.spectra.set_xlim((self.q[0], self.q[-1]))
        self.spectra.set_ylim((-0.1, 1.1))
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        self.spectra.set_ylabel("Avg intensity (a.u.)")

        self.draw()
       

        if stick_patterns is not None:
            self.plot_cifs(cif_patterns)

        self.draw()

    def plot_fit_result(self, fit_result, bg=None):
        fit = evaluate_obj(fit_result.phase_model, self.avg_q)

        phase_name = []
        if bg is not None:
            fit += bg
        else:
            bg = evaluate_obj(fit_result.phase_model.background, self.avg_q)

        self.spectra.plot(self.avg_q, fit, label="Fitted")
        self.spectra.plot(self.avg_q, bg, label="background")

        for cp in fit_result.phase_model.CPs:
            self.spectra.plot(self.avg_q, evaluate_obj(cp, self.avg_q), label=cp.name)
            phase_name.append(cp.name)

        self.spectra.set_title("_".join(phase_name))
        self.spectra.legend(fontsize=7)
        self.spectra.set_xlim((self.q[0], self.q[-1]))
        self.spectra.set_ylim((-0.1, 1.1))
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        self.spectra.set_ylabel("Avg intensity (a.u.)")

        self.draw()


    def replot_spectra(self, fit_result=None, bg=None):
        self.spectra.clear()
        if self.avg_pattern is None:
            self.avg_pattern = minmax_norm(self.data[:,round(self.data.shape[1]/2)])

        (self.avgplot, ) = self.spectra.plot(self.avg_q, self.avg_pattern,
                                             linewidth=2, color='k', label="XRD")
        (self.spectra_select_box, ) = self.spectra.plot(self.spectra_box_x,
                                                        self.spectra_box_y, color='r')

        self.plot_fit_result(fit_result, bg)


    def plot_cifs(self, sticks):
        """ 
        stick_patterns: dictionary: {phase: [q, peak_height]}
        """
        self.stick_patterns.clear()
        proxies = []
        for i, phase in enumerate(sticks):
            stick = sticks[phase] 
            qs = stick[:,0]
            Is = stick[:,1]

            boxes = [Rectangle((q-0.1, 0), .2, I)
                               for q, I in zip(qs, Is)]
            pc = PatchCollection(boxes, facecolor=COLORS[i], alpha=1.)
            proxy = mlines.Line2D([], [], marker="_", linewidth=1,color=COLORS[i], markersize=15, label=phase)
            proxies.append(proxy)
            self.stick_patterns.add_collection(pc)
            

        self.stick_patterns.plot()
        self.stick_patterns.legend(handles=proxies, fontsize=7)
        self.stick_patterns.set_xlim((self.q[0], self.q[-1]))
        self.stick_patterns.set_ylim((0, 1))
        self.draw()
