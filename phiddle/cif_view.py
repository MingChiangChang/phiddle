from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QVBoxLayout, QWidget, QPushButton, QCheckBox
)
# from matplotlib.backends.backend_qtagg import (
#     FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar, )
# from matplotlib.figure import Figure

# TODO: 1. Instead of combo box, use a list that user can modify the order of 
#          phases for convenience (QListWidget is a little bit difficult to use)

class CIFView(QWidget):
    """
    A display widget for showing a list of CIF files with checkboxes, allowing users to add
    or remove selected CIF files in the stripeview.

    signals:
        checked: pyqtSignal(list) = Emitted when a checkbox state changes, passing a list of boolean values.
        add: pyqtSignal(list) = Emitted when "Add to phase diagram" is clicked, passing a list of selected CIF names.
        remove: pyqtSignal() = Emitted when "Remove from phase diagram" is clicked.
    """

    checked = pyqtSignal(list)
    add = pyqtSignal(list)
    remove = pyqtSignal()

    def __init__(self, cif_list, parent=None):
        """
        Initialize the CIFView widget with the provided list of CIF files.

        parameters:
            cif_list: list
                A list of CIF file names to be displayed as checkboxes.
            parent: QWidget, optional
                The parent widget for this CIFView. Defaults to None.
        """

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
            checkbox.stateChanged.connect(self.update_stick_pattern)
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
        """
        Update the CIF list by adding or removing checkboxes for the CIF files.

        parameters:
            cif_list: list
                A list of CIF file names to update the checkboxes for.
        """
        self.layout.removeWidget(self.add_button)
        self.layout.removeWidget(self.remove_button)
        self.layout.removeWidget(self.amorphous_checkbox)
        self.layout.removeWidget(self.melt_checkbox)
        del self.widget_ls[-1]  # amorphous checkbox is always the last one
        del self.widget_ls[-1]  # melt checkbox is second to last 
        for idx, cif in enumerate(cif_list):
            if idx >= len(self.widget_ls): # Save 2 spots for Default options
                checkbox = QCheckBox(cif)
                checkbox.stateChanged.connect(self.update_stick_pattern)
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
        """
        Emit a signal with the current checkbox states, indicating which CIFs are selected.

        This signal is emitted every time a checkbox is toggled, providing the current
        state of all checkboxes (excluding the last two).
        """
        self.checked.emit([checkbox.isChecked()
                          for checkbox in self.widget_ls[:-2]])

    def add_to_phase_diagram(self):
        """
        Emit a signal to add selected CIF files to the phase diagram.

        The signal is emitted when the "Add to phase diagram" button is clicked, passing
        a list of selected CIF file names to the connected slots.
        """
        self.add.emit([checkbox.text() for checkbox in self.widget_ls if checkbox.isChecked()])


    def remove_from_phase_diagram(self):
        """
        Emit a signal to remove CIF files from the phase diagram.

        The signal is emitted when the "Remove from phase diagram" button is clicked.
        """
        self.remove.emit()

    def clear(self):
        """
        Uncheck all checkboxes.

        This function clears the selection by unchecking all the checkboxes in the list.
        """
        for checkbox in self.widget_ls:
            checkbox.setChecked(False)

    def check_boxes(self, ind_ls):
        """
        Check the checkboxes at the specified indices.

        parameters:
            ind_ls: list
                A list of indices corresponding to the checkboxes that should be checked.
        """
        for ind in ind_ls:
            self.widget_ls[ind].setChecked(True)
        self.update_stick_pattern()

    def get_checked_phase_names(self):
        """
        Get the names of the selected phases.

        return:
            list: A list of phase names corresponding to the selected checkboxes,
                  excluding the "Amorphous" and "Melt" options.
        """
        return [checkbox.text() for checkbox in self.widget_ls
                if checkbox.isChecked()
                if checkbox.text() != "Amorphous" and checkbox.text() != "Melt"]

    def set_checked_phase_names(self, phase_names):
        """
        Set the checkboxes based on the provided phase names.

        parameters:
            phase_names: list
                A list of phase names to check. Checkboxes corresponding to these phase names will be checked.
        """
        for checkbox in self.widget_ls:
            if checkbox.text() in phase_names:
                checkbox.setChecked(True)
            else:
                checkbox.setChecked(False)

