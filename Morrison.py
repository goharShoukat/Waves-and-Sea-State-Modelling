#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon 15 Feb 2021
@author: goharshoukat
'''

import math
import numpy as np
import Linear_Airy_Wave_Solution as laws

class Morisons:
    
    def __init__(self, field):
        self.field = field
    
    
    
    
        
field = laws.Set_Wave_Field(1,2,0,9)
force_solver = Morisons(field)
