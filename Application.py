#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 23:02:26 2021

@author: goharshoukat
"""

from GUI_Design import Ui_Form
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
import Linear_Airy_Wave_Solution as laws
import numpy as np

class Define_Properties(qtw.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.ui.run_button.clicked.connect(self.calculate)
        
    def calculate(self):
        period = float(self.ui.period_input.text())
        amplitude = float(self.ui.amp_input.text())
        epsilon = float(self.ui.epsilon_input.text())
        if epsilon > 0.3:
             qtw.QMessageBox.information(self, 'Error', 'Wave Steepness is high, Please confirm non-linearity')
        depth = float(self.ui.depth_input.text())
        tidal_velocity = float(self.ui.tidal_input.text())
        width = float(self.ui.width_input.text())
        dx = float(self.ui.dx_input.text())
        dz = float(self.ui.dz_input.text())
        k = period * amplitude#float(period) * float(amplitude)
        ans_str = '{0:0.2f}'.format(k)
        qtw.QMessageBox.information(self, 'Success', ans_str)
        x = np.arange(0, width + dx, dx)
        field = laws.Wave_Field(depth, period, tidal_velocity, amplitude, x, dz)
        
        
        



if __name__ == '__main__':
    app = qtw.QApplication([])
    widget = Define_Properties()
    widget.show()
    app.exec()
    