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
SootWSR.WSR.reactor_properties(vol=np.pi*(1.17**2/4)*9.2*(.0254)**3)
SootWSR.WSR.set_inlet_gas(300,101325,'POSF5433:1 O2:1',.06)
SootWSR.WSR.set_outlet('PC',500*101325/14.7)
SootWSR.WSR.solve(1e-5,1,convergence_criteria=1e-3)



#SootWSR.WSR.integrate_coupled(1e-3,1,convergence_criteria=1e-3)

keyboard()


AA.vol = reactor_vol
AA.init_gas(T0,Pc,gas.Y)
AA.def_inlet_gas(T0,Pc,gas.Y)

# 1838, 1799, 1801
tres = 7e-3

AA.def_outlet('PC',Pc)
AA.WSR_tres(7e-3)#,Treactor=1750)
AA.mdot = AA.gas.density*AA.vol/tres
#AA.solve_energy = False

keyboard()
AA.integrate_coupled(1e-5,1,convergence_criteria=1e-3)
#AA.integrate_to_SS(1000)
print(AA.MOMIC.soot_properties(AA.M))

d = AA.MOMIC.soot_properties(AA.M)[0]
d1 = AA.MOMIC.soot_properties(AA.M)[1]
fv = AA.MOMIC.soot_properties(AA.M)[2]
N = AA.M[0]
M1 = AA.M[1]

data = [d,d1,fv,N,M1]
keyboard()
np.save('4_coal_unc.npy',data)

#fv_1bar = np.append(fv_1bar,AA.MOMIC.soot_properties(AA.M)[2])
#np.save('data_1bar.npy',[dat_1atm[:j,0],fv_1bar])
#keyboard()


