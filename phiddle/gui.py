import os
import sys

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QVBoxLayout, QWidget, QPushButton
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
import h5py

f = h5py.File("data/AL_23F4_Bi-Ti-O_run_01_0_all_1d.h5")


class MplWidget(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        fig = Figure()
        super(MplWidget, self).__init__(fig)
        self.setParent(parent)

        self.h5 = h5py.File(
            "data/AL_23F4_Bi-Ti-O_run_01_0_all_1d.h5",
            'r')['exp']
        self.conds = list(self.h5)
        self.ind = 0
        self.q, self.data = self.collect_data_and_q(
            self.h5, self.conds[self.ind])

        gs = fig.add_gridspec(1, 2)
        self.ax1 = fig.add_subplot(gs[0, 0])
        self.ax2 = fig.add_subplot(gs[0, 1])

        self.bottomLeftX = 0
        self.bottomLeftY = -100
        self.topRightX = 0
        self.topRightY = 100

        self.x = np.array([self.bottomLeftX, self.bottomLeftX,
                          self.topRightX, self.topRightX, self.bottomLeftX])
        self.y = np.array([self.bottomLeftY, self.topRightY,
                          self.topRightY, self.bottomLeftY, self.bottomLeftY])

        (self.heatmap, ) = self.ax1.plot(self.x, self.y, color='r')
        self.ax1.imshow(self.data,
                        extent=(0, self.data.shape[1], self.q[0], self.q[-1]),
                        aspect=self.data.shape[1] / (self.q[-1] - self.q[0]))
        self.ax1.set_title(self.conds[0])
        self.aspan = self.ax1.axvspan(
            self.bottomLeftX, self.topRightX, color='k', alpha=0)

        self.avg_pattern = self.data[:, round(self.data.shape[1] / 2)]
        (self.avgplot, ) = self.ax2.plot(self.q, self.avg_pattern, color='r')

        self.moving = False
        self.plotSnap = 1

        self.cid1 = self.mpl_connect("button_press_event", self.onclick)
        self.cid2 = self.mpl_connect("button_release_event", self.onrelease)
        self.cid3 = self.mpl_connect("motion_notify_event", self.onmotion)

    @property
    def current_cond(self):
        return self.conds[self.ind]

    @property
    def filename(self):
        cond = self.current_cond
        ion = self.h5[cond].attrs['cations'][0]
        frac = self.h5[self.current_cond].attrs['fracs'][0][0]
        return f"{cond}_{ion}_{frac:.3f}"

    def next(self):
        self.ind += 1
        self.q, self.data = self.collect_data_and_q(
            self.h5, self.conds[self.ind])
        self.ax1.clear()
        self.ax2.clear()
        self.bottomLeftX = 0
        self.bottomLeftY = -100
        self.topRightX = 0
        self.topRightY = 100
        self.x = np.array([self.bottomLeftX, self.bottomLeftX,
                          self.topRightX, self.topRightX, self.bottomLeftX])
        self.y = np.array([self.bottomLeftY, self.topRightY,
                          self.topRightY, self.bottomLeftY, self.bottomLeftY])
        (self.heatmap, ) = self.ax1.plot(self.x, self.y, color='r')
        self.ax1.imshow(self.data,
                        extent=(0, self.data.shape[1], self.q[0], self.q[-1]),
                        aspect=self.data.shape[1] / (self.q[-1] - self.q[0]))
        self.ax1.set_title(self.conds[self.ind])
        self.aspan = self.ax1.axvspan(
            self.bottomLeftX, self.topRightX, color='k', alpha=0)

        self.avg_pattern = self.data[:, round(self.data.shape[1] / 2)]
        (self.avgplot, ) = self.ax2.plot(self.q, self.avg_pattern, color='r')
        self.draw()

    @classmethod
    def collect_data_and_q(cls, h5, cond):
        dim2 = len(h5[cond])
        dim1 = len(h5[cond]['0']['integrated_1d'][0])

        q = h5[cond]['0']['integrated_1d'][0]
        arr = np.zeros((dim1, dim2))
        for i in range(dim2):
            arr[:, i] = h5[cond][str(i)]['integrated_1d'][1]

        return q, arr

    # def update(self):
    #    self.q, self.data = self.collect_data_and_q(self.h5, self.conds[self.ind])
    #    self.ax1.imshow(self.data,
    #                   extent=(0, self.data.shape[1], self.q[0], self.q[-1]),
    #                   aspect=self.data.shape[1]/(self.q[-1]-self.q[0]))
    #    #self.draw()

    def setSnapBase(self, base):
        return lambda value: int(base * round(float(value) / base))

    def onclick(self, event):
        if self.plotSnap <= 0:
            self.bottomLeftX = event.xdata
            self.bottomLeftY = event.ydata
        else:
            self.calculateSnapCoordinates = self.setSnapBase(self.plotSnap)
            self.bottomLeftX = self.calculateSnapCoordinates(event.xdata)
            self.bottomLeftY = self.calculateSnapCoordinates(event.ydata)

        try:
            self.aspan.remove()
        except BaseException:
            pass

        self.moving = True

    def onrelease(self, event):
        if self.plotSnap <= 0:
            self.topRightX = event.xdata
            self.topRightY = event.ydata

        else:
            try:
                calculateSnapCoordinates = self.setSnapBase(self.setSnapBase)
                self.topRightX = calculateSnapCoordinates(event.xdata)
                self.topRightY = calculateSnapCoordinates(event.ydata)
            except BaseException:
                pass

        self.x = np.array([self.bottomLeftX, self.bottomLeftX,
                          self.topRightX, self.topRightX, self.bottomLeftX])

        self.heatmap.set_xdata(self.x)

        self.aspan = self.ax1.axvspan(
            self.bottomLeftX, self.topRightX,
            0, 1,
            color='b',
            alpha=0.25,
        )

        self.moving = False
        if self.bottomLeftX > self.topRightX:
            self.bottomLeftX, self.topRightX = self.topRightX, self.bottomLeftX
        self.avg_pattern = np.mean(
            self.data[:, self.bottomLeftX:self.topRightX], axis=1)
        self.avgplot.set_ydata(self.avg_pattern)

        self.draw()

    def onmotion(self, event):
        if not self.moving:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return

        if self.plotSnap <= 0:
            self.topRightX = event.xdata
            self.topRightT = event.ydata
        else:
            self.calculateSnapCoordinates = self.setSnapBase(self.plotSnap)
            self.topRightX = self.calculateSnapCoordinates(event.xdata)
            self.topRightY = self.calculateSnapCoordinates(event.ydata)

        self.x = np.array([self.bottomLeftX, self.bottomLeftX,
                          self.topRightX, self.topRightX, self.bottomLeftX])

        self.heatmap.set_xdata(self.x)
        self.draw()


class TopLevelWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.canvas = MplWidget()

        button1 = QPushButton()
        button1.setText("Save")
        button1.clicked.connect(self.button1_clicked)

        button2 = QPushButton()
        button2.setText("Next")
        button2.clicked.connect(self.button2_clicked)

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(button1)
        layout.addWidget(button2)
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

    def button1_clicked(self):
        filename = self.canvas.filename
        np.save(filename, self.canvas.avg_pattern)

    def button2_clicked(self):
        self.canvas.next()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    w = TopLevelWindow()
    w.show()

    sys.exit(app.exec())
