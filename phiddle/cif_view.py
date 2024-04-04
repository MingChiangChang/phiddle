from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QVBoxLayout, QWidget, QPushButton, QCheckBox
)
# from matplotlib.backends.backend_qtagg import (
#     FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar, )
from matplotlib.figure import Figure


class CIFView(QWidget):

    checked = pyqtSignal(list)
    add = pyqtSignal(list)

    def __init__(self, cif_list, parent=None):

        super(CIFView, self).__init__(parent)

        self.layout = QVBoxLayout()

        self.widget_ls = []
        self.button = QPushButton()
        self.button.setText("Add to phase diagram")
        self.button.clicked.connect(self.add_to_phase_diagram)
        for cif in cif_list:
            checkbox = QCheckBox(cif)
            checkbox.clicked.connect(self.update_stick_pattern)
            self.widget_ls.append(checkbox)
            self.layout.addWidget(checkbox)
        self.amorphous_checkbox = QCheckBox("Amorphous")
        # Should not let amorphous
        self.widget_ls.append(self.amorphous_checkbox)
        self.layout.addWidget(self.amorphous_checkbox)  # Be hidden in the list

        self.layout.addWidget(self.button)
        self.setLayout(self.layout)

    def update_cif_list(self, cif_list):
        """ Reuse checkboxes  """
        self.layout.removeWidget(self.button)
        self.layout.removeWidget(self.amorphous_checkbox)
        del self.widget_ls[-1]  # amorphous checkbox is always the last one
        for idx, cif in enumerate(cif_list):
            if idx >= len(self.widget_ls):
                checkbox = QCheckBox(cif)
                checkbox.clicked.connect(self.update_stick_pattern)
                self.widget_ls.append(checkbox)
                self.layout.addWidget(checkbox)
            else:
                checkbox = self.widget_ls[idx]
                checkbox.setText(cif)
        for idx, widget in enumerate(self.widget_ls):
            if idx >= len(cif_list):
                self.layout.removeWidget(self.widget_ls[idx])
        diff = len(self.widget_ls) - len(cif_list)
        for idx in range(diff):
            del self.widget_ls[-1]  # always remove the last one

        self.widget_ls.append(self.amorphous_checkbox)
        self.layout.addWidget(self.amorphous_checkbox)
        self.layout.addWidget(self.button)
        self.clear()

    def update_stick_pattern(self):
        self.checked.emit([checkbox.isChecked()
                          for checkbox in self.widget_ls[:-1]])

    def add_to_phase_diagram(self):
        self.add.emit([checkbox.isChecked() for checkbox in self.widget_ls])

    def clear(self):
        for checkbox in self.widget_ls:
            checkbox.setChecked(False)

    def check_boxes(self, ind_ls):
        for ind in ind_ls:
            self.widget_ls[ind].setChecked(True)
        self.update_stick_pattern()

    def get_checked_phase_names(self):
        return [checkbox.text() for checkbox in self.widget_ls
                if checkbox.isChecked()
                if checkbox.text() != "Amorphous"]
