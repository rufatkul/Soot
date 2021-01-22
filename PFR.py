# Cantera Plug Flow Reactor Code
# Written by: Rufat Kulakhmetov


from momic_class import * 
import time

# Initialize Soot Object
SootPFR = Soot('./Mechanisms/R2highT.cti')
# Initialize #M Momments, #P Moments, Regime:(0: based on Kndusen,1:FM,2:Tr,3:Co)
SootPFR.initialize(6,0,2)
# Set Initial Moments (Arbitrary but may help with initial convergence)
SootPFR.MOMIC.set_moments(np.exp([10,11,12,13,14,15]))#,np.exp([1,1,1]))

# Reactor Properties
Ac = np.pi*(1.17*.0254)**2/4
z_pts = np.linspace(0.3,9.2,50)*.0254
mdot = 0.03

# Reactor Properties
SootPFR.PFR.set_reactor_pts(z_pts,Ac)
SootPFR.PFR.set_reactor_properties(Twall=300)
# Inlet Moments
SootPFR.PFR.set_inlet_Moments([0,0,0,0,0,0])#,[0,0,0])
# Pressure Outlet => PC: Constant Pressure (Specifiy P), mass: compressible flow calculation (Specify At, P)
SootPFR.PFR.set_PFR_outlet('PC',At=np.pi*(.308*.0254)**2/4,P= 250*101325/14.7)# 241*101325/14.7)

# Inlet Gas or inlet liquid
SootPFR.PFR.set_inlet_gas(300,101325*250/14.7,'O2:1',mdot,z=0)
SootPFR.PFR.set_inlet_gas(300,101325*250/14.7,'POSF5433:1',mdot,z=0)
#SootPFR.PFR.set_inlet_liq(T=300,gas='POSF5433:1',D0=200*1e-6,mu=1e-4,mdot=0.1,rho=810,LHV=251000,Tboil=350,convection=True)


tic = time.time()
SootPFR.PFR.solve()
toc = time.time() - tic
print('Elapsed Time')
print(toc)

# Solution Stored Here using Cantera SolutionArray Object
Sol = SootPFR.PFR.Sol
keyboard()
