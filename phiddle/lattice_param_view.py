import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QVBoxLayout, QWidget, QPushButton, QTabWidget, QFormLayout,
    QCheckBox, QFileDialog, QComboBox, QHBoxLayout, QLabel
)
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar, )
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mpl_toolkits

from util import COLORS

class LatticeParamView(FigureCanvasQTAgg):

    def __init__(self, parent=None, xlim=(250, 10000), ylim=(400, 1400)):

        self.fig = Figure()
        super(LatticeParamView, self).__init__(self.fig)
        plt.tight_layout()
        self.setParent = parent

        self.xlabel = "Dwell (μs)"
        self.ylabel = "Tpeak ($^o$C)"
        self.zlabel = None
        self.font_size = 6
        self.xlim = xlim
        self.ylim = ylim

        self.gs = self.fig.add_gridspec(2, 3)
        self._init_plot()
        self._set_title()
        self._set_ticks()
        self._set_labels()
        
        self.dim = 2
        self.lp_dict = {} # Purposely use stateful widget to separate the weight from main view manager

    def _init_plot(self, projection=None):
        self.a_plot = self.fig.add_subplot(self.gs[0, 0], projection=projection)
        self.b_plot = self.fig.add_subplot(self.gs[0, 1], projection=projection)
        self.c_plot = self.fig.add_subplot(self.gs[0, 2], projection=projection)
        self.α_plot = self.fig.add_subplot(self.gs[1, 0], projection=projection)
        self.β_plot = self.fig.add_subplot(self.gs[1, 1], projection=projection)
        self.γ_plot = self.fig.add_subplot(self.gs[1, 2], projection=projection)
        self.plots = [self.a_plot, self.b_plot, self.c_plot, self.α_plot, self.β_plot, self.γ_plot]

        self.cbars = []
        if projection is None:
            for i, plot in enumerate(self.plots):
                div = make_axes_locatable(plot)
                cax = div.append_axes("right", size="5%", pad=0.1)
                if projection == "3d":
                    im = plot.scatter([0.], [0.], [0.], c=[0.])
                else:
                    im = plot.scatter([0.], [0.], c=[0.])
                plot.set_xscale("log")
                plot.set_xlim(*self.xlim)
                plot.set_ylim(*self.ylim)
                if isinstance(im, mpl_toolkits.mplot3d.art3d.Path3DCollection):
                    cbar = self.fig.colorbar(im, cax=cax, ticklocation='top') 
                else:
                    cbar = self.fig.colorbar(im, cax=cax)
                cbar.ax.tick_params(labelsize=self.font_size)
                self.cbars.append(cbar)


    def _set_title(self):
        self.a_plot.set_title('a')
        self.b_plot.set_title('b')
        self.c_plot.set_title('c')
        self.α_plot.set_title('α')
        self.β_plot.set_title('β')
        self.γ_plot.set_title('γ')


    def _set_labels(self):
        self.a_plot.set_ylabel(self.ylabel)
        self.α_plot.set_ylabel(self.ylabel)
        self.α_plot.set_xlabel(self.xlabel)
        self.β_plot.set_xlabel(self.xlabel)
        self.γ_plot.set_xlabel(self.xlabel)
        

    def _set_ticks(self):
        self.b_plot.set_yticks([])
        self.c_plot.set_yticks([])
        self.β_plot.set_yticks([])
        self.γ_plot.set_yticks([])


    def clear_figures(self):
        # for i, plot in enumerate(self.plots):
        for plot, cbar in zip(self.plots, self.cbars):
            for image in plot.images:
                image.colorbar.clear()
            cbar.remove()
            plot.clear()
            plot.remove()

    def _init_figures(self, projection=None):
        self._init_plot(projection)
        self._set_title()
        self._set_ticks()
        self._set_labels()


    # def plot(self, tpeak, dwell, refined_lp, axes = ["Dwell", "Tpeak"], mask=None):
    def plot(self, lp_dict, axes = ["Dwell", "Tpeak"]):
        self.lp_dict = lp_dict
        self.clear_figures()
        # if self.dim == 2:
        self._init_figures()
        xlim, self.xlabel, xscale = self.get_2d_axis_info(axes[0])
        ylim, self.ylabel, yscale = self.get_2d_axis_info(axes[1])

        for i, plot in enumerate(self.plots):
            if i > 2: 
                c = lp_dict["refined_lps"][:, i]/np.pi*180
                im = plot.scatter(lp_dict[axes[0]], lp_dict[axes[1]], c=c)
            else:
                c = lp_dict["refined_lps"][:, i]
                im = plot.scatter(lp_dict[axes[0]], lp_dict[axes[1]], c=c)
            im.set_clim(np.mean(c)-2*np.std(c), np.mean(c)+2*np.std(c))
            im.colorbar = self.cbars[i]
            im.colorbar_cid = im.callbacks.connect('changed', self.cbars[i].update_normal) # So that the removal does not give errors
            self.cbars[i].update_normal(im)

            plot.set_xscale(xscale)
            plot.set_xlim(xlim)
            plot.set_yscale(yscale)
            plot.set_ylim(ylim)
        self._set_labels()
        self._set_ticks()

        self.draw()


    def get_log_axis(self, plots, scales):
        idx = scales.index('log')
        if idx == 0:
            return plots.xaxis
        elif idx == 1:
            return plots.yaxis
        elif idx == 2:
            return plots.zaxis


    def log_tick_formatter(self, val, pos=None):
        return f"{int(np.round(10**val, -1))}"

    def identity(self, x):
        return x

    def get_2d_axis_info(self, axis):
        if axis == "Tpeak":
            lim = (400, 1400)
            label = f"Peak Temperature ($^oC$)"
            scale = "linear"
        elif axis == "Dwell":
            lim = (250, 10000)
            label = f"Dwell (us)"
            scale = "log"
        elif axis =="x" or axis =="y":
            lim = (-48, 48)
            label = (axis + 'mm')
            scale = 'linear'
        else:
            lim = (0, 1)
            label = f"{axis} (%)"
            scale = "linear"
        return lim, label, scale


    def get_3d_axis_info(self, axis):
        if axis == "Tpeak":
            lim = (400, 1400)
            label = f"Peak Temperature ($^oC$)"
            scale = "linear"
        elif axis == "Dwell":
            lim = (np.log10(250), np.log10(10000))
            label = f"Dwell (us)"
            scale = "log"
        else:
            lim = (0, 1)
            label = f"{axis} (%)"
            scale = "linear"
        return lim, label, scale



    def save_phase_diagram(self):
        fn, _ = QFileDialog.getSaveFileName(self, 'Save Phase Diagram', "", "")
        self.phase_diagram.figure.savefig(fn)


    def change_dim(self, value, axes):
        self.dim = value
        self.plot(self.lp_dict, axes)

    def change_axes(self, axes):
        self.plot(self.lp_dict, axes)




class LatticeParamList(QWidget):

    checked_signal = pyqtSignal(list)
    save_signal = pyqtSignal()
    save_lp_signal = pyqtSignal(str)
    axes_signal = pyqtSignal(list)
    dim_change_signal = pyqtSignal(int, list)

    def __init__(self, composition_dim=0, parent=None):

        # Better way will be to create it just before showing it, but
        # Need to figure out how to do it in pyqt

        super(LatticeParamList, self).__init__(parent)

        self.save_button = QPushButton()
        self.save_button.setText("Save Diagram")
        self.save_button.clicked.connect(lambda: self.save_signal.emit())
        self.save_lp_button = QPushButton()
        self.save_lp_button.setText("Save Lattice Parameters")
        self.save_lp_button.clicked.connect(
                lambda: self.save_lp_signal.emit(self.phase_selection_box.currentText())
        )

        self.dim_selection_box = QComboBox()
        self.dim_selection_box.addItems(["2D phase digrams", "3D phase diagrams"])
        self.dim_selection_box.setCurrentIndex(0)
        self.dim_selection_box.currentIndexChanged.connect(self.dim_changed)

        self.comp1_str = "Composition 1"
        self.comp2_str = "Composition 2"
        self.comp3_str = "Composition 3"
        self.comp_str = [self.comp1_str, self.comp2_str, self.comp3_str]
        self.composition_dim = composition_dim

        self.widget_ls = []

        option_ls = self.get_option_ls()
        self.axis1_layout = QHBoxLayout()
        self.axis1_selection_box = QComboBox()
        self.axis1_selection_box.addItems(option_ls)
        self.axis1_layout.addWidget(QLabel("Axis 1"))
        self.axis1_layout.addWidget(self.axis1_selection_box)
        self.axis1_selection_box.currentIndexChanged.connect(self.axes_changed)

        self.axis2_layout = QHBoxLayout()
        self.axis2_selection_box = QComboBox()
        self.axis2_selection_box.addItems(option_ls)
        self.axis2_selection_box.setCurrentIndex(1)
        self.axis2_layout.addWidget(QLabel("Axis 2"))
        self.axis2_layout.addWidget(self.axis2_selection_box)
        self.axis2_selection_box.currentIndexChanged.connect(self.axes_changed)

        self.axis3_layout = QHBoxLayout()
        self.axis3_selection_box = QComboBox()
        self.axis3_selection_box.addItems(option_ls)
        self.axis3_selection_box.setCurrentIndex(2)
        self.axis3_label = QLabel("Axis 3")
        self.axis3_layout.addWidget(self.axis3_label)
        self.axis3_layout.addWidget(self.axis3_selection_box)
        self.axis3_selection_box.currentIndexChanged.connect(self.axes_changed)

        self.phase_selection_layout = QHBoxLayout()
        self.phase_selection_box = QComboBox()
        self.phase_selection_label = QLabel("Phase")
        self.phase_selection_layout.addWidget(self.phase_selection_label)
        self.phase_selection_layout.addWidget(self.phase_selection_box)

        self.outer_layout = QVBoxLayout()
        self.layout = QVBoxLayout()
        self.outer_layout.addWidget(self.dim_selection_box)
        self.outer_layout.addLayout(self.axis1_layout)
        self.outer_layout.addLayout(self.axis2_layout)
        self.outer_layout.addLayout(self.axis3_layout)
        self.outer_layout.addLayout(self.phase_selection_layout)

        self.axis3_label.hide()
        self.axis3_selection_box.hide()

        self.outer_layout.addLayout(self.layout)
        self.outer_layout.addWidget(self.save_button)
        self.outer_layout.addWidget(self.save_lp_button)
        self.setLayout(self.outer_layout)


    def update_axis_combo_boxes(self):
        option_ls = self.get_option_ls()
        self.update_combo_box(self.axis1_selection_box, option_ls)
        self.update_combo_box(self.axis2_selection_box, option_ls)
        self.update_combo_box(self.axis3_selection_box, option_ls)

    def update_combo_box(self, combo_box, option_ls):
        for i, option in enumerate(option_ls):
            if i < combo_box.count():
                combo_box.setItemText(i, option)
            else:
                combo_box.addItem(option)
        while combo_box.count() > len(option_ls):
            combo_box.removeItem(combo_box.count()-1) # Remove the last item

    def update_phase_combo_box(self, phase_names):
        self.phase_selection_box.setCurrentIndex(0)
        self.update_combo_box(self.phase_selection_box, phase_names)

    def update_phase_lists(self, phases):

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
            for i in range(len(phase), len(self.widget_ls)):
                self.layout.removeWidget(self.widget_ls[i])
                del self.widget_ls[i]


    def update_phase_diagram(self):
        self.checked_signal.emit([checkbox.isChecked()
                          for checkbox in self.widget_ls])


    def dim_changed(self):
        if self.dim_selection_box.currentIndex() == 0:
            self.axis3_label.hide()
            self.axis3_selection_box.hide()
            dim = 2
            axes = [box.currentText() for box in [self.axis1_selection_box,
                                                 self.axis2_selection_box]]
        elif self.dim_selection_box.currentIndex() == 1:
            self.axis3_label.show()
            self.axis3_selection_box.show()
            self.axis3_selection_box.setCurrentIndex(2)
            dim = 3
            axes = [box.currentText()
                    for box in [self.axis1_selection_box,
                                self.axis2_selection_box,
                                self.axis3_selection_box]]

        self.dim_change_signal.emit(dim, axes)


    def axes_changed(self):
        axes = self.get_current_axes()
        self.axes_signal.emit(axes)


    def get_current_axes(self):
        axes = [box.currentText() for box in [self.axis1_selection_box,
                                             self.axis2_selection_box]]
        if self.dim_selection_box.currentIndex() == 1:
            axes.append(self.axis3_selection_box.currentText())
        return axes


    def get_option_ls(self):
        ls = ["Dwell", "Tpeak", "x", "y"]
        i = 0
        while i < self.composition_dim:
            ls.append(self.comp_str[i])
            i += 1
        return ls
