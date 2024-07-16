from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QVBoxLayout, QWidget, QPushButton, QCheckBox
)
# from matplotlib.backends.backend_qtagg import (
#     FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar, )
# from matplotlib.figure import Figure

# TODO: 1. Instead of combo box, use a list that user can modify the order of 
#          phases for convenience
#       2. When user choosed the range with phase, autocheck the boxes for those
#          phases and show cifs 

class CIFView(QWidget):

    checked = pyqtSignal(list)
    add = pyqtSignal(list)
    remove = pyqtSignal()

    def __init__(self, cif_list, parent=None):

        super(CIFView, self).__init__(parent)

        self.layout = QVBoxLayout()

        self.widget_ls = []
        self.add_button = QPushButton()
        self.add_button.setText("Add to phase diagram")
        self.add_button.clicked.connect(self.add_to_phase_diagram)

        self.remove_button = QPushButton()
        self.remove_button.setText("Remove from phase diagram")
        self.remove_button.clicked.connect(self.remove_from_phase_diagram)

        for cif in cif_list:
            checkbox = QCheckBox(cif)
            checkbox.clicked.connect(self.update_stick_pattern)
            self.widget_ls.append(checkbox)
            self.layout.addWidget(checkbox)

        self.amorphous_checkbox = QCheckBox("Amorphous")
        self.melt_checkbox = QCheckBox("Melt")
        # Should not let amorphous
        self.widget_ls.append(self.amorphous_checkbox)
        self.widget_ls.append(self.melt_checkbox)
        self.layout.addWidget(self.amorphous_checkbox)  # Be hidden in the list
        self.layout.addWidget(self.melt_checkbox)  # Be hidden in the list

        self.layout.addWidget(self.add_button)
        self.layout.addWidget(self.remove_button)
        self.setLayout(self.layout)

    def update_cif_list(self, cif_list):
        """ Reuse checkboxes  """
        self.layout.removeWidget(self.add_button)
        self.layout.removeWidget(self.remove_button)
        self.layout.removeWidget(self.amorphous_checkbox)
        self.layout.removeWidget(self.melt_checkbox)
        del self.widget_ls[-1]  # amorphous checkbox is always the last one
        del self.widget_ls[-1]  # melt checkbox is second to last 
        for idx, cif in enumerate(cif_list):
            if idx >= len(self.widget_ls): # Save 2 spots for Default options
                checkbox = QCheckBox(cif)
                checkbox.clicked.connect(self.update_stick_pattern)
                self.widget_ls.append(checkbox)
                self.layout.addWidget(checkbox)
            else:
                checkbox = self.widget_ls[idx]
                checkbox.setText(cif)
        for idx, _ in enumerate(self.widget_ls):
            if idx >= len(cif_list):
                self.layout.removeWidget(self.widget_ls[idx])
        diff = len(self.widget_ls) - len(cif_list)
        for idx in range(diff):
            del self.widget_ls[-1]  # always remove the last one

        self.widget_ls.append(self.amorphous_checkbox)
        self.widget_ls.append(self.melt_checkbox)
        self.layout.addWidget(self.amorphous_checkbox)
        self.layout.addWidget(self.melt_checkbox)
        self.layout.addWidget(self.add_button)
        self.layout.addWidget(self.remove_button)
        self.clear()

    def update_stick_pattern(self):
        self.checked.emit([checkbox.isChecked()
                          for checkbox in self.widget_ls[:-2]])

    def add_to_phase_diagram(self):
        self.add.emit([checkbox.isChecked() for checkbox in self.widget_ls])


    def remove_from_phase_diagram(self):
        self.remove.emit()

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
                if checkbox.text() != "Amorphous" and checkbox.text() != "Melt"]

    def set_checked_phase_names(self, phase_names):
        for checkbox in self.widget_ls:
            if checkbox.text() in phase_names:
                checkbox.setChecked(True)
            else:
                checkbox.setChecked(False)

