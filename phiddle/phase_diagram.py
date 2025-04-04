import numpy as np
from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtWidgets import (
    QVBoxLayout, QWidget, QPushButton, QCheckBox, QFileDialog,
    QComboBox, QHBoxLayout, QLabel
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.ticker as mticker
from scipy.spatial import ConvexHull
from scipy.spatial._qhull import QhullError

from util import COLORS

# TODO: Allow user to choose volume plot and simplex phase regions
# TODO: Combine phase diagram list signal to just one, much easier to sync

class PhaseDiagramView(FigureCanvasQTAgg):
    """
    Phase Digram View that can switch between 2d and 3d phase diagrams

    Attributes:
        xlim: 2Tuple
            x limits of the phase digram
        ylim: 2Tuple
            y limits of the phase digram
        dim: int
            dimension of the phase diagram. Only 2 or 3 is valid
        phase_dict: dict
            Dictionary for phase diagram info. Comes from the data model
        phase_list: list
            Store phases that should be plotted in colors
    """

    def __init__(self, xlim=(250., 10000.), ylim=(400., 1400.)):
        """ 
        Initialize parameters.


        Args:
            xlim: Tuple[float, float]
                x limits of the plot

            ylim: Tuple[float, float]
                x limits of the plot
        """

        self.fig = Figure(tight_layout=True)
        super(PhaseDiagramView, self).__init__(self.fig)

        self.phase_diagram = self.fig.add_subplot() # Default empty plot
        self.xlim = xlim
        self.ylim = ylim
        self.dim = 2
        self.phase_dict = {} # Purposely use stateful widget to separate the weight from main view manager
        self.phase_list = []
        # self.draw()

    def plot(self, phase_dict, axes = ["Dwell", "Tpeak"],
             phase_list=None, plot_convex_hull=False):
        """
        Plot the phase diagram based on the dictionary provided

        Args:
            phase_dict: dict
                Dictionary for phase diagram info. Comes from the data model
            axes: List[string]
                Listing the axes to be plotted
            phase_list: List[string]
                Listing the phases that should be plotted in colors.
                Other phases will be plotted in light gray
            plot_convex_hull: bool
                Flag of whether to plot convex hulls for phases listed in phase_list
        """

        self.phase_dict = phase_dict
        self.phase_list = phase_list
        phase_name_ls = np.array(list(phase_dict))

        others = []
        if phase_list is not None:
            others = [p for p in phase_name_ls if p not in phase_list]
            phase_name_ls = phase_list 

        if self.dim == 2: # For 2D phase digrams
            if self.fig.axes:
                self.fig.gca().remove()
            self.phase_diagram = self.fig.add_subplot()

            xlim, xlabel, xscale = self.get_2d_axis_info(axes[0])
            ylim, ylabel, yscale = self.get_2d_axis_info(axes[1])

            for idx, phase in enumerate(phase_name_ls):
                x = np.array(phase_dict[phase][axes[0]])
                y = np.array(phase_dict[phase][axes[1]])
                self.phase_diagram.scatter(x, y, label=phase,
                                           color=COLORS[((idx+1) % len(COLORS)-1)],
                                           alpha=0.5)
                if plot_convex_hull: # Plot 2D convex hull of the phase region
                    try:
                        hull = ConvexHull(np.vstack((x, y)).T)

                        for simplex in hull.simplices:
                            self.phase_diagram.plot(x[simplex], y[simplex],
                                          color=COLORS[((idx+1) % len(COLORS)-1)],
                                          alpha=0.5)
                        
                        self.phase_diagram.fill(x[hull.vertices], y[hull.vertices],
                                       color=COLORS[((idx+1) % len(COLORS)-1)],
                                       alpha=0.05)
                    except QhullError:
                        print("Linear points for {phase}. Unable to create convex hullregion.")

            if phase_list is not None:
                for idx, phase in enumerate(others):
                    self.phase_diagram.scatter(phase_dict[phase][axes[0]],
                                               phase_dict[phase][axes[1]],
                                               label="Other" if idx==0 else "_other",
                                               color="#DFDFDF", 
                                               alpha=1., 
                                               s=4.)

            self.phase_diagram.set_xlim(xlim)
            self.phase_diagram.set_xlabel(xlabel)
            self.phase_diagram.set_xscale(xscale)
            self.phase_diagram.set_ylim(ylim)
            self.phase_diagram.set_ylabel(ylabel)
            self.phase_diagram.set_yscale(yscale)
            self.phase_diagram.legend(bbox_to_anchor=(1., 1.))

        elif self.dim == 3: # For 3D phase diagrams
            if self.fig.axes:
                self.fig.gca().remove()
                cids = sorted(list(self.callbacks._pickled_cids))
                if len(cids) > 8: # There are 8 default callbacks for some reason
                   for i in range(3):
                       # Remove the callbacks that were connected to previous axes
                       self.callbacks.disconnect(cids[-i-1])

            self.phase_diagram = self.fig.add_subplot(projection='3d')

            xlim, xlabel, xscale = self.get_3d_axis_info(axes[0])
            ylim, ylabel, yscale = self.get_3d_axis_info(axes[1])
            zlim, zlabel, zscale = self.get_3d_axis_info(axes[2]) # there's inconsistnecy here
            scales = [xscale, yscale, zscale]

            transform = []
            for scale in scales:
                if scale == "log":
                    transform.append(np.log10)
                else:
                    transform.append(self.identity)


            for idx, phase in enumerate(phase_name_ls):
                self.phase_diagram.scatter(transform[0](phase_dict[phase][axes[0]]),
                                      transform[1](phase_dict[phase][axes[1]]),
                                      zs = transform[2](phase_dict[phase][axes[2]]),
                                      label=phase,
                                      color=COLORS[((idx+1) % len(COLORS)-1)],
                                      s=50,
                                      alpha=0.5)
            if phase_list is not None:
                for idx, phase in enumerate(others):
                    self.phase_diagram.scatter(transform[0](phase_dict[phase][axes[0]]),
                                               transform[1](phase_dict[phase][axes[1]]),
                                               zs = transform[2](phase_dict[phase][axes[2]]),
                                               label="Other" if idx==0 else "_other",
                                               s=10,
                                               color="#DFDFDF",
                                               alpha=1.)
            # Note: no convexhull plotting for 3D currently

            self.phase_diagram.set_xlim(xlim)
            self.phase_diagram.set_xlabel(xlabel)
            self.phase_diagram.set_ylim(ylim)
            self.phase_diagram.set_ylabel(ylabel)
            self.phase_diagram.set_zlim(zlim)
            self.phase_diagram.set_zlabel(zlabel)

            if 'log' in scales:
                self.format_log_tick(self.get_log_axis(scales))
            self.phase_diagram.legend()
        self.draw()

    def format_log_tick(self, axis):
        axis.set_major_formatter(mticker.FuncFormatter(self.log_tick_formatter))
        axis.set_major_locator(mticker.FixedLocator(
            np.linspace(np.log10(250), np.log10(10000), 10)
            )
         )

    def get_log_axis(self, scales):
        idx = scales.index('log')
        if idx == 0:
            return self.phase_diagram.xaxis
        elif idx == 1:
            return self.phase_diagram.yaxis
        elif idx == 2:
            return self.phase_diagram.zaxis


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
        self.plot(self.phase_dict, axes, self.phase_list)

    def change_axes(self, axes):
        self.plot(self.phase_dict, axes, self.phase_list)




class PhaseDiagramList(QWidget):
    """ 
    Will be displayed next to the phase diagram
    Contains pull down menu for selected dimensionality of the phase diagram,
    x, y, z axis selection and check boxes to choose what phase to show

    Signals:
        checked_signal(list): send what phases are checked once the user check
                              or uncheck a box
        save_signal(): trigger saving of the phase diagram plot
        axes_signal(list): send what the axes of the phase diagram should be
                           triggered when the axes menu is changed by the user
        dim_change_signal(int, list): 
            triggered when dimensionality pull down menu is changed by the user.
            Send the dimension and the axes
    """

    checked_signal = pyqtSignal(list)
    save_signal = pyqtSignal()
    axes_signal = pyqtSignal(list)
    dim_change_signal = pyqtSignal(int, list)

    def __init__(self, composition_dim=0, parent=None):
        """ initialize widget """

        # Better way will be to create it just before showing it, but
        # Need to figure out how to do it in pyqt

        super(PhaseDiagramList, self).__init__(parent)

        self.save_button = QPushButton()
        self.save_button.setText("Save Phase Diagram")
        self.save_button.clicked.connect(lambda: self.save_signal.emit())

        self.dim_selection_box = QComboBox()
        self.dim_selection_box.addItems(["2D phase digrams", "3D phase diagrams"])
        self.dim_selection_box.setCurrentIndex(0)
        self.dim_selection_box.currentIndexChanged.connect(self.dim_changed)

        self.comp1_str = "Composition 1"
        self.comp2_str = "Composition 2"
        self.comp3_str = "Composition 3"
        self.comp_str = [self.comp1_str, self.comp2_str, self.comp3_str]
        self.composition_dim = composition_dim
        # self.option_ls = []

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

        self.plot_selection_layout = QHBoxLayout()
        self.plot_selection_box = QComboBox()
        self.plot_selection_box.addItems(["Scatter plot", "Convex Hull"])
        self.plot_selection_layout.addWidget(QLabel("Plot Type"))
        self.plot_selection_layout.addWidget(self.plot_selection_box)
        self.plot_selection_box.currentIndexChanged.connect(self.plot_changed)

        self.outer_layout = QVBoxLayout()
        self.layout = QVBoxLayout()
        self.outer_layout.addWidget(self.dim_selection_box)
        self.outer_layout.addLayout(self.axis1_layout)
        self.outer_layout.addLayout(self.axis2_layout)
        self.outer_layout.addLayout(self.axis3_layout)
        self.outer_layout.addLayout(self.plot_selection_layout)

        self.check_all_box = QCheckBox("Select All")
        self.check_all_box.stateChanged.connect(self.select_all)

        self.axis3_label.hide()
        self.axis3_selection_box.hide()

        self.outer_layout.addLayout(self.layout)
        self.outer_layout.addWidget(self.check_all_box)
        self.outer_layout.addWidget(self.save_button)
        self.setLayout(self.outer_layout)

    def select_all(self, check_state_int):
        """ make all the checkboxes visible be checed or unchecked  """
        if check_state_int == 2:
            for w in self.widget_ls:
                w.setCheckState(Qt.CheckState.Checked)
        else:
            for w in self.widget_ls:
                w.setCheckState(Qt.CheckState.Unchecked)
        self.update_phase_diagram()
 
    def update_combo_boxes(self):
        """ Update the pull down menu to the proper options  """
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


    def get_option_ls(self):
        """
        Get the propert option list.
        This includes: ["Dwell", "Tpeak", *cations]
        """
        ls = ["Dwell", "Tpeak"]
        i = 0
        while i < self.composition_dim:
            ls.append(self.comp_str[i])
            i += 1
        return ls


    def _show(self, phases):

        ordered_phase = []
        if "Amorphous" in phases:
            ordered_phase.append("Amorphous")
            phases.remove("Amorphous")
        if "Melt" in phases:
            ordered_phase.append("Melt")
            phases.remove("Melt")
        ordered_phase = sorted(phases) + ordered_phase

        for idx, phase in enumerate(ordered_phase):
            if idx >= len(self.widget_ls):
                checkbox = QCheckBox(phase)
                # checkbox.setChecked(True)
                checkbox.clicked.connect(self.update_phase_diagram)
                self.widget_ls.append(checkbox)
                self.layout.addWidget(checkbox)
            else:
                checkbox = self.widget_ls[idx]
                checkbox.setText(phase)

        if len(self.widget_ls) > len(ordered_phase):
            for i in range(len(ordered_phase), len(self.widget_ls)):
                self.layout.removeWidget(self.widget_ls[-1])
                del self.widget_ls[-1]


    def update_phase_diagram(self):
        """ emit what boxes are being checked  """
        self.checked_signal.emit(self.get_checked_phase_names())

    def get_checked_phase_names(self):
        """ return names of phases that were checked """
        return [checkbox.text() for checkbox in self.widget_ls if checkbox.isChecked()]

    def dim_changed(self):
        """ Trigger things when the dimension is being changed  """
        if self.dim_selection_box.currentIndex() == 0:
            self.axis3_label.hide() # Hide the pull down menu for the 3rd axis
            self.axis3_selection_box.hide()
            dim = 2
            axes = [box.currentText() for box in [self.axis1_selection_box,
                                                 self.axis2_selection_box]]
        elif self.dim_selection_box.currentIndex() == 1:
            self.axis3_label.show() # show the pull down menu for the 3rd axis
            self.axis3_selection_box.show()
            self.axis3_selection_box.setCurrentIndex(2)
            dim = 3
            axes = [box.currentText()
                    for box in [self.axis1_selection_box,
                                self.axis2_selection_box,
                                self.axis3_selection_box]]

        self.dim_change_signal.emit(dim, axes)

    @property
    def convex_hull(self):
        return self.plot_selection_box.currentText() == "Convex Hull"

    def axes_changed(self):
        """ collect the current axes options and emit """
        axes = self.get_current_axes()
        self.axes_signal.emit(axes)

    def get_current_axes(self):
        """ collect the current axes options """
        axes = [box.currentText() for box in [self.axis1_selection_box,
                                             self.axis2_selection_box]]
        if self.dim_selection_box.currentIndex() == 1:
            axes.append(self.axis3_selection_box.currentText())
        return axes


    def plot_changed(self):
        pass
