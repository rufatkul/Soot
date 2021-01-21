# Cantera Plug Flow Reactor Code
# Written by: Rufat Kulakhmetov

#import cantera as ct
#import numpy as np
#import math
#import csv
#import copy
#import os.path
#from pdb import set_trace as keyboard
from momic_class import * 
import time


SootPFR = Soot('./Mechanisms/R2highT.cti')
SootPFR.initialize(6,0,2)
SootPFR.MOMIC.set_moments(np.exp([10,11,12,13,14,15]))#,np.exp([1,1,1]))
SootPFR.MOMIC.update_regime(2,1)

Ac = np.pi*(1.17*.0254)**2/4
z_pts = np.linspace(0.3,9.2,50)*.0254
mdot = 0.03

SootPFR.PFR.set_reactor_pts(z_pts,Ac)
SootPFR.PFR.set_reactor_properties(Twall=300)
SootPFR.PFR.set_inlet_Moments([0,0,0,0,0,0])#,[0,0,0])
#SootPFR.PFR.pressure_guess(250*101325/14.7)
#SootPFR.solve_soot = False
SootPFR.PFR.set_PFR_outlet('mass',At=np.pi*(.308*.0254)**2/4,P= 101325)# 241*101325/14.7)#500*101325/14.7)
#SootPFR.PFR.set_PFR_outlet('PC',At=np.pi*(.308*.0254)**2/4,P= 250*101325/14.7)# 241*101325/14.7)
SootPFR.PFR.set_inlet_gas(300,101325*250/14.7,'O2:1',mdot,z=0)
SootPFR.PFR.set_inlet_gas(300,101325*250/14.7,'POSF5433:1',mdot,z=0)
#SootPFR.PFR.set_inlet_liq(T=300,gas='POSF5433:1',D0=200*1e-6,mu=1e-4,mdot=0.1,rho=810,LHV=251000,Tboil=350,convection=True)


tic = time.time()
SootPFR.PFR.solve()
toc = time.time() - tic
print(toc)
keyboard()
