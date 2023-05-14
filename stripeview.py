import os
import sys

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines

from util import minmax_norm, COLORS
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

        self.bottomLeftX = 0
        self.bottomLeftY = -100
        self.topRightX = 0
        self.topRightY = 100

        self.avg_pattern = None
        self.moving = False

        self.cid1 = self.mpl_connect("button_press_event", self.onclick)
        self.cid2 = self.mpl_connect("button_release_event", self.onrelease)
        self.cid3 = self.mpl_connect("motion_notify_event", self.onmotion)
        
    @property
    def x(self):
        return np.array([self.bottomLeftX,
                         self.bottomLeftX,
                         self.topRightX,
                         self.topRightX,
                         self.bottomLeftX], dtype=int)

    @property
    def y(self):
        return np.array([self.bottomLeftY,
                         self.topRightY,
                         self.topRightY,
                         self.bottomLeftY,
                         self.bottomLeftY], dtype=int)


    def onclick(self, event):

        if event.inaxes in [self.heatmap]:
            self.bottomLeftX = int(event.xdata)
            self.bottomLeftY = int(event.ydata)

            try:
                self.aspan.remove()
            except:
                pass

            self.moving = True

    def onrelease(self, event):
        if event.inaxes in [self.heatmap]:
            self.topRightX = int(event.xdata)
            self.topRightY = int(event.ydata)

            self.selection_box.set_xdata(self.x)

            self.aspan = self.heatmap.axvspan(
                    self.bottomLeftX, self.topRightX,
                    0,1,
                    color = 'b',
                    alpha = 0.25,
            )

            self.moving = False
            if self.bottomLeftX > self.topRightX:
                self.bottomLeftX, self.topRightX = self.topRightX, self.bottomLeftX
            # Let this be here for now
            if self.bottomLeftX == self.topRightX: 
                self.avg_pattern = self.data[:, self.bottomLeftX] 
            else:
                self.avg_pattern = np.mean(self.data[:,self.bottomLeftX:self.topRightX], axis=1)
            self.avg_pattern = minmax_norm(self.avg_pattern)
            self.spectra.clear()
            self.spectra.plot(self.q, self.avg_pattern, color='r', label="XRD")
            self.spectra.set_xlim((self.q[0], self.q[-1]))
            self.spectra.set_ylabel("Avg intensity (a.u.)")
            self.spectra.legend()

            self.draw()
    
    def onmotion(self, event):
        if not self.moving:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return

        self.topRightX = int(event.xdata)
        self.topRightT = int(event.ydata)


        self.selection_box.set_xdata(self.x)
        self.draw()

    def clear_figures(self):
        self.heatmap.clear()
        self.spectra.clear()

        self.bottomLeftX = 0
        self.bottomLeftY = -100
        self.topRightX = 0
        self.topRightY = 100

    def plot(self, data, stick_patterns=None, fit_result=None):
        self.clear_figures()
        self.q = data['q']
        self.data = data['data']
        self.cond = data['cond']

        (self.selection_box, ) = self.heatmap.plot(self.x, self.y, color='r')
        self.heatmap.imshow(self.data,
                       extent=(0, self.data.shape[1], self.q[0], self.q[-1]),
                       aspect=self.data.shape[1]/(self.q[-1]-self.q[0]))
        self.heatmap.set_title(self.cond)
        self.heatmap.set_xlabel("index")
        self.heatmap.set_ylabel("q ($nm^{-1}$)")       


        self.aspan = self.heatmap.axvspan(self.bottomLeftX, self.topRightX, color='k', alpha=0)

        if self.avg_pattern is None:
            self.avg_pattern = minmax_norm(self.data[:,round(self.data.shape[1]/2)])
        (self.avgplot, ) = self.spectra.plot(self.q, self.avg_pattern, color='r', label="XRD")

        phase_name = []
        if fit_result is not None:
            
            self.spectra.plot(self.q, evaluate_obj(fit_result.phase_model, self.q), label="Fitted")
            self.spectra.plot(self.q, evaluate_obj(fit_result.phase_model.background, self.q),
                              label="background")
            for cp in fit_result.phase_model.CPs:
                self.spectra.plot(self.q, evaluate_obj(cp, self.q), label=cp.name)
                phase_name.append(cp.name)
        self.spectra.set_title("_".join(phase_name))
        self.spectra.legend()
        self.spectra.set_xlim((self.q[0], self.q[-1]))
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        self.spectra.set_ylabel("Avg intensity (a.u.)")

        if stick_patterns is not None:
            self.plot_cifs(cif_patterns)

        self.draw()

    def replot_spectra(self, fit_result=None):
        if self.avg_pattern is None:
                self.avg_pattern = minmax_norm(self.data[:,round(self.data.shape[1]/2)])
    
        phase_name = []
        if fit_result is not None:
    
            self.spectra.plot(self.q, evaluate_obj(fit_result.phase_model, self.q), label="Fitted")
            self.spectra.plot(self.q, evaluate_obj(fit_result.phase_model.background, self.q),
                               label="background")
            for cp in fit_result.phase_model.CPs:
                self.spectra.plot(self.q, evaluate_obj(cp, self.q), label=cp.name)
                phase_name.append(cp.name)
        self.spectra.set_title("_".join(phase_name))
        self.spectra.legend()
        self.spectra.set_xlim((self.q[0], self.q[-1])) 
        self.spectra.set_xlabel("q ($nm^{-1}$)")
        self.spectra.set_ylabel("Avg intensity (a.u.)")
    
        self.draw()

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
