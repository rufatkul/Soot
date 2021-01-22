# Cantera Plug Flow Reactor Code
# Written by: Rufat Kulakhmetov

import cantera as ct
import numpy as np
import math
import csv
import copy
import os.path
from pdb import set_trace as keyboard
from momic_class import * 


SootWSR = Soot('./Mechanisms/R2highT.cti')
SootWSR.initialize(6,0,0)
SootWSR.WSR.set_moments(np.exp([10,11,12,13,14,15]))
SootWSR.WSR.set_inlet_Moments([0,0,0,0,0,0])
SootWSR.WSR.set_reactor_properties(vol=np.pi*(1.17**2/4)*9.2*(.0254)**3)
SootWSR.WSR.set_inlet_gas(300,101325,'POSF5433:1 O2:1',.06)
SootWSR.WSR.set_outlet('PC',500*101325/14.7)
SootWSR.WSR.solve(1e-5,1,convergence_criteria=1e-3)

# Solution is Stored in SootWSR.WSR.sol
Sol = SootWSR.WSR.Sol
keyboard()
