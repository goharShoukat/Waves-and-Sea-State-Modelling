#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 21:07:26 2021

@author: goharshoukat
"""

import GUI
import Linear_Airy_Wave_Solution
import sys

from PyQt5 import QtCore, QtGui, QtWidgets
app = QtWidgets.QApplication(sys.argv)
WaveFieldGUI = QtWidgets.QMainWindow()
ui = GUI.Ui_WaveFieldGUI()
ui.setupUi(WaveFieldGUI)
WaveFieldGUI.show()
sys.exit(app.exec_())
