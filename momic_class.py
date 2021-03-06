# SOOT Object Oriented Toolbox
# By: Rufat Kulakhmetov

# Developed During my PhD
# Still has some junk code in it that was used for developing debugging
# If a class method is not used it probably was used to debug

import sys
# Add Path To Fortran Momic Module
sys.path.insert(1,'./MOMIC/')

# Import All the Necessary Packages
import numpy as np
import momic_module as momic
from pdb import set_trace as keyboard
import cantera as ct
import scipy.integrate
import copy
import time


class Soot:
# Main Soot Object Class, containing all other objects

	def __init__(self,reaction_mechanism,surface_mechanism=None):

		# Store reaction_mechanism, surface_mechanism
		self.reaction_mechanism = reaction_mechanism
		self.surface_mechanism = surface_mechanism

		# Create Gas Object
		self.gas = ct.Solution(reaction_mechanism)		
		
		# Flags to control shit
	#	self.solve_coupled = True
	#	self.solve_soot = True
	#	self.constant_pressure = False
	#	self.aggregation = True
		#keyboard()


		# Control Flags
		# ** Add any main control flags for any future development here and then reference to them here**
		# ** This makes debugging a lot easier **
		
		self.solve_coupled = True	# Solve Soot Coupled to Reactor Equations? - Need to Initialize MOMIC Module
		self.solve_soot = False		# Initialize MOMIC Module to do soot 
		self.solve_energy = True	# Solves Energy Equation in Reactor
		self.solve_heat_transfer = False


		# Associated Object Classes
		# Add reference to this main Class
		self.MOMIC = MOMIC(gas_object=self.gas, main_class=self)
		self.deposit = None #deposition(main_class=self)	
		self.droplet = None #droplet(main_class=self)
		self.Sol = None 	# Reference to solution for WSR or PFR
		self.WSR = WSR(reaction_mechanism, main_class=self)
		self.PFR = PFR(reaction_mechanism,surface_mechanism,WSR_object=self.WSR,main_class=self)
				
	
	def initialize(self,Num_M_Moments,Num_P_Moments,regime):
		# Use this method to initialize MOMIC module instead of using input file
		# Inputs are:
		# Num_M_Moments: Number of M Moments to Use
		# Num_P_Moments: Number of P Moments to Use
		# regime: Coagulation regime calculation - 
		#			0: Calculate Based on Knudsen Number
		#			1: Free Molecular
		#			2: Transitional
		#			3: Continuum

		# Initialize MOMIC Module	
		self.MOMIC.initialize(Num_M_Moments,Num_P_Moments,regime)

		# Switch Flag to Solve Soot
		self.solve_soot = True

		if Num_P_Moments == 0:
			self.aggregation = False
		else:
			self.aggregation = True
				 

class MOMIC:
# F90 Method of Moments Subroutine

	def __init__(self,M=[],P=[],gas_object=None,main_class=None):
		# If using this in integrator, pass M/P as reference!
		self.M = []
		self.P = []
		self.Num_M_Moments = None
		self.Num_P_Moments = None
		self.C2H2 = []
		self.H = []
		self.H2 = []
		self.H2O = []
		self.O2 = []
		self.OH = []
		self.Temperature = []
		self.Pressure = []
		self.M_moments = []
		self.P_moments = []
		self.aggregation = False

		# Species Production Rates (Defined after calling calculate_source method)
		self.rC2H2 = []
		self.rCO = []
		self.rH = []
		self.rH2 = []
		self.rH2O = []
		self.rO2 = []
		self.rOH = []
		
		# Source Terms (Defined after calling calculate_source method)
		self.W = []
		self.G = []
		self.R = []
		self.Hr = []
		self.Ragg = []
		self.SM = []
		self.SP = []

		# Check if main class is defined, add references if it is
		if main_class is not None:
			self.main = main_class
			self.gas = self.main.gas
		else:
			self.main = None
			self.gas = gas_object

		if self.gas is not None:
			self.wdot_soot = np.zeros(self.gas.n_species)
		else:
			self.wdot_soot = None


		# Run set Moments
		self.set_moments(M, P)

		# Flags to control MOMIC Source Calculation
		self.aggregation = False
		self.regime = 2	# 2:Free-Molecular 1:Transitional 0:Continuum
		self.condition_switch_flag = False

		# Control for when to turn on aggregation
		self.aggregate_Dstar = 10e-9

	def moment_dict(self):
		# Return Moments as Dictionary
		# Useful for using with Cantera SolutionArray
		
		# Make M and P Dictionaries
		M_Moments = {'M{}'.format(i):M for i,M in enumerate(self.M)}
		P_Moments = {'P{}'.format(i):P for i,P in enumerate(self.P)}

		# Add P Dict to M Dict to make a single Dict
		M_Moments.update(P_Moments)

		return M_Moments

	def initialize(self,Num_M_Moments,Num_P_Moments=0,regime=0):
		# Use this method to initialize MOMIC module
		# Inputs are:
		# Num_M_Moments: Number of M Moments to Use
		# Num_P_Moments: Number of P Moments to Use
		# regime: Coagulation regime calculation - 
		#			0: Calculate Based on Knudsen Number
		#			1: Free Molecular
		#			2: Transitional
		#			3: Continuum
		
		momic.momic.initialize(Num_M_Moments,Num_P_Moments,regime)

		# Assign these attributes
		self.Num_M_Moments = Num_M_Moments
		self.Num_P_Moments = Num_P_Moments

		# Assign some random initial M,P Moments if they're none
		
		if self.M or np.size(self.M) != Num_M_Moments:
			self.M = np.exp(np.arange(1,Num_M_Moments+1,1))

		if Num_P_Moments>0:
			self.aggregation = True
		else:
			self.aggregation = False
			
		if self.P or np.size(self.P) != Num_P_Moments:
			if Num_P_Moments>0:
				self.P = np.exp(np.arange(1,Num_P_Moments+1,1))

	def update_regime(self,regime,force_aggregate=False):
		# Update Regime or Turn on/off aggregation after
		# MOMIC has already been initialized
		# Inputs are:
		# regime: Coagulation regime calculation - 
		#			0: Calculate Based on Knudsen Number
		#			1: Free Molecular
		#			2: Transitional
		#			3: Continuum
		# force_aggregate: True/False

		momic.momic.update_regime(regime,force_aggregate)
	
	def set_moments(self,M,P=[]):
		# Set Moments

		self.M = np.array(M)
		self.P = np.array(P)
	
	def set_ct_gas(self,gas_object):
		# If using in integrator, pass in gas_object as reference!
		self.gas = gas_object
		
		# Vector for storing soot rates
		self.wdot_soot = np.zeros(self.gas.n_species)

	def ct_rate(self):
		
		self.Temperature = self.gas.T
		self.Pressure = self.gas.P
		
		# Convert Cantera Concentration from kmol/m^3 to mol/cm^3
		# MOMIC Units are in mol/cm^3
		self.C2H2 = self.gas.concentrations[self.gas.species_index('C2H2')]*1e-3
		self.H = self.gas.concentrations[self.gas.species_index('H')]*1e-3
		self.H2 = self.gas.concentrations[self.gas.species_index('H2')]*1e-3
		self.H2O = self.gas.concentrations[self.gas.species_index('H2O')]*1e-3
		self.O2 = self.gas.concentrations[self.gas.species_index('O2')]*1e-3
		self.OH = self.gas.concentrations[self.gas.species_index('OH')]*1e-3
		
		# Run MOMIC 
		self.calculate_source()
		
		# Define vector of soot consumption rates
		self.wdot_soot[:] = 0
		self.wdot_soot[self.gas.species_index('C2H2')] = self.rC2H2
		self.wdot_soot[self.gas.species_index('CO')] = self.rCO
		self.wdot_soot[self.gas.species_index('H')] = self.rH
		self.wdot_soot[self.gas.species_index('H2')] = self.rH2
		self.wdot_soot[self.gas.species_index('H2O')] = self.rH2O
		self.wdot_soot[self.gas.species_index('O2')] = self.rO2
		self.wdot_soot[self.gas.species_index('OH')] = self.rOH
		
		# Output species rates in mol/cm^3, convert to kmol/m^3 (Cantera units)
		self.wdot_soot*=1e3
		
		return self.SM, self.SP, [self.W,self.G,self.R,self.Hr,self.Ragg], self.wdot_soot

	def calculate_source(self):
		# Calculate FORTRAN MOMIC Code to Calculate MOMIC Rates

		# Call Fortran Code
		Out = momic.momic.calculate_source(self.M,self.P,	\
				self.gas.T, self.gas.P,						\
				self.C2H2, 									\
				self.H, 									\
				self.H2, 									\
				self.H2O, 									\
				self.O2, 									\
				self.OH )									\
		
		# Split MOMIC Rates Outputs
		self.W = Out[0]
		self.G = Out[1]
		self.R = Out[2]
		self.Hr = Out[3]
		self.Ragg = Out[4]
	
		# Total Rates
		self.SM = self.W+self.G+self.R
		self.SP = self.Hr + self.Ragg
		
		# Species Rates
		self.rC2H2 = Out[5]
		self.rCO = Out[6]
		self.rH = Out[7]
		self.rH2 = Out[8]
		self.rH2O = Out[9]
		self.rO2 = Out[10]
		self.rOH = Out[11]
	
		return self.SM, self.SP, [self.W,self.G,self.R,self.Hr,self.Ragg], Out[5:]

	def k_soot(self,T):
		# Calculate soot thermal conductivity
		# From paper on soot measurments
	
		k = 13.0839 - 0.03495*T+3.82e-5*T**2-1.48e-8*T**3

		return np.max([k,1])

	def kn(self):
		# Calculate knudsent number

		D_soot = self.soot_properties(self.M)[0]
		mfp = self.mean_free_path(self.gas.P,self.gas.T)

		return 2*mfp/D_soot
	
	def update_regime_calc(self):

		# Calculate kn #
		kn = self.kn()

		# Check what regime to set 
		if kn <=0.1:
			regime= 3	# Continuum
		elif kn <= 10: 
			regime = 2	# Tr
		elif kn > 10: 
			regime = 1	# FM

		if self.regime_update_approach == 0:
			if regime >= self.regime:
				self.regime = regime 
		elif self.regime_update_approach == 1:
			if regime <= self.regime:
				self.regime = regime 
		else:
			self.regime = regime


		# Now Check Aggregation
		if self.aggregation:
			# Calculate D of average soot particle
			Davg = self.Davg_soot_particle(self.M)
		
			if Davg< self.aggregate_Dstar and not self.aggregate:
				self.aggregate = False
				self.P[:] = self.M[0]
			else:
				if not self.aggregate  and not self.condition_switch_flag:
					self.condition_switch_flag = True
				else:
					self.condition_switch_flag = False

				self.aggregate = True

		keyboard()



	def update_regime_flags(self,reset=False):
		# Update source term Regime flags based on soot size

		# Reset Flags
		if reset:
			self.regime = 2
			self.aggregate = False
			return

		# Calculate kn
		kn = self.kn()



		if self.regime == regime:
			self.condition_switch_flag = False
		else:
			self.condition_switch_flag = True
			self.regime = regime

		# Now Check Aggregation
		if self.aggregation:
			# Calculate D of average soot particle
			Davg = self.Davg_soot_particle(self.M)

			if Davg< self.aggregate_Dstar and not self.aggregate:
				self.aggregate = False
				self.P[:] = self.M[0]
			else:
				if not self.aggregate  and not self.condition_switch_flag:
					self.condition_switch_flag = True
				else:
					self.condition_switch_flag = False

				self.aggregate = True

		else:
			self.aggregate = False



	@staticmethod
	def Davg_soot_particle(M):
		# Calculate Diameter of average soot particle (NOT AVERAGE DIAMETER)

		# Carbon Mass in grams
		mc = 12*1.67e-24 
		# Diameter of Average Soot Particle (in m)
		return (6/np.pi*mc*M[1]/M[0]/1.8)**(1/3)*1e-2

	@staticmethod
	def soot_properties(M):

		# Calculate Diameter of average particle
		Davg = MOMIC.Davg_soot_particle(M)

		# Append Average Soot Particle Diameter to Average properties from moment calc
		return np.append(Davg,momic.momic.soot_properties(M))
	
	@staticmethod
	def mean_free_path(P,T):
		# Return Mean Free Path for air in meters
		return momic.momic.mfp(P,T)*1e-2
	
	@staticmethod
	def soot_mfr(M,rho_gas):
		return momic.momic.soot_y_calc(M,rho_gas)

	def calculate_source_switcher(self):
		# Calculate MOMIC rates, but select coagulation regime and
		# aggregation using this function, rather than within the 
		# MOMIC F90 code
		#
		# Use this function with scipy ode if multiple regimes 
		# could be present

		# Control regime and aggregate with these properties:
		# self.regime
		# self.aggregate

		# AGGREGATION AND REGIEM SET TO 1 and True
		print('AGG AND TRANSITIONAL REGIME TURNED ON')
		self.regime = 1
		self.aggregate = True

		# Call Fortran Code
		Out = momic.momic.calculate_source_switcher(self.M,self.P,	\
				self.gas.T, self.gas.P,						\
				self.C2H2, 									\
				self.H, 									\
				self.H2, 									\
				self.H2O, 									\
				self.O2, 									\
				self.OH,self.regime,self.aggregate )									\
		
		# Split MOMIC Rates Outputs
		self.W = Out[0]
		self.G = Out[1]
		self.R = Out[2]
		self.Hr = Out[3]
		self.Ragg = Out[4]
	
		# Total Rates
		self.SM = self.W+self.G+self.R
		self.SP = self.Hr + self.Ragg
		
		# Species Rates
		self.rC2H2 = Out[5]
		self.rCO = Out[6]
		self.rH = Out[7]
		self.rH2 = Out[8]
		self.rH2O = Out[9]
		self.rO2 = Out[10]
		self.rOH = Out[11]
	
		return self.SM, self.SP, [self.W,self.G,self.R,self.Hr,self.Ragg], Out[5:]

	def ct_rate_switcher(self):
		# Same as ct rate, but uses calculate_source_switcher
		self.Temperature = self.gas.T
		self.Pressure = self.gas.P
		
		# Convert Cantera Concentration from kmol/m^3 to mol/cm^3
		# MOMIC Units are in mol/cm^3
		self.C2H2 = self.gas.concentrations[self.gas.species_index('C2H2')]*1e-3
		self.H = self.gas.concentrations[self.gas.species_index('H')]*1e-3
		self.H2 = self.gas.concentrations[self.gas.species_index('H2')]*1e-3
		self.H2O = self.gas.concentrations[self.gas.species_index('H2O')]*1e-3
		self.O2 = self.gas.concentrations[self.gas.species_index('O2')]*1e-3
		self.OH = self.gas.concentrations[self.gas.species_index('OH')]*1e-3
		
		# Run MOMIC 
		self.calculate_source_switcher()
		
		# Define vector of soot consumption rates
		self.wdot_soot[:] = 0
		self.wdot_soot[self.gas.species_index('C2H2')] = self.rC2H2
		self.wdot_soot[self.gas.species_index('CO')] = self.rCO
		self.wdot_soot[self.gas.species_index('H')] = self.rH
		self.wdot_soot[self.gas.species_index('H2')] = self.rH2
		self.wdot_soot[self.gas.species_index('H2O')] = self.rH2O
		self.wdot_soot[self.gas.species_index('O2')] = self.rO2
		self.wdot_soot[self.gas.species_index('OH')] = self.rOH
		
		# Output species rates in mol/cm^3, convert to kmol/m^3 (Cantera units)
		self.wdot_soot*=1e3
		
		return self.SM, self.SP, [self.W,self.G,self.R,self.Hr,self.Ragg], self.wdot_soot

		
class WSR():

	def __init__(self,reaction_mechanism,MOMIC_object= None,main_class=None):
		
		# Add Reference to main class objects
		if main_class is not None:
			self.main = main_class
			self.reaction_mechanism = self.main.reaction_mechanism
			self.MOMIC = self.main.MOMIC
			self.gas = self.main.gas
			self.droplet = self.main.droplet
			
			# Control Flags
			self.solve_soot = self.main.solve_soot
			self.solve_coupled = self.main.solve_coupled
			self.solve_energy = self.main.solve_energy
			self.solve_heat_transfer = self.main.solve_heat_transfer 

		else:
			self.main = None
			self.reaction_mechanism = reaction_mechanism
			self.MOMIC = MOMIC_object
			self.gas = ct.Solution(self.reaction_mechanism)
			self.droplet = None

			# Control Flags
			self.solve_soot = False
			self.solve_coupled = False
			self.solve_energy = True 
			self.solve_heat_transfer = False
			
			if self.MOMIC is not None:
				self.solve_soot = True
				self.solve_coupled = True

		self.inlet_gas = ct.Solution(self.reaction_mechanism)
		self.inlet_gas_mdot = None
		self.inlet_liq = None
		self.inlet_liq_mdot = None
		self.inlet_liq_D0 = None
		self.outlet_gas = ct.Solution(self.reaction_mechanism)
		self.outlet_BC = None
		self.inlet_M = None
		self.inlet_P = None
		self.vol = None
		self.tres = None
		self.M = self.MOMIC.M
		self.P = self.MOMIC.P
		self.M_Moments = self.MOMIC.Num_M_Moments
		self.P_Moments = self.MOMIC.Num_P_Moments
		
		self.mdot_gas = 0
		self.mdot_liq = 0
		self.mdot_vap = 0
		self.mdot_soot = 0
		self.mdot_eat = 0


		self.count = 0	# Number of Times the Function is called 
		
		# Geometry Information for other calculations
		self.At = None
		self.Ac = None
		self.Twall = None
		self.Dc = None

		#self.rho_mdot = None

		# WSR Flow Velocity
		# (Used in droplet calculation with evaporation)
		self.U = 0

		# Other Controller Flags 
		# (Turned on when liquid inlet is specified)
		self.solve_droplet = False


	def update_control_flags(self,caller=None):
		# Update Control Flags in WSR

		if caller is not None:
			self.solve_coupled = caller.solve_coupled
			self.constant_pressure = caller.constant_pressure
			self.solve_soot = caller.solve_soot
		elif self.main is not None:
			self.solve_coupled = self.main.solve_coupled
			self.constant_pressure = self.main.constant_pressure
			self.solve_soot = self.main.solve_soot
		else:
			pass 

	def init_gas(self,T0,P0,Y0):
		self.gas.TPY = T0, P0, Y0 
	
	def set_inlet_gas(self,T0,P0,Y0,mdot):
		# WSR Gas Inlet

		# Update Inlet Gas Properties
		self.inlet_gas.TPY = T0, P0, Y0
		self.inlet_gas_mdot = mdot
	
	def set_inlet_liq(self,T=None,gas=None,D0=None,mdot=None,rho=None,mu=None,LHV=None,Tboil=None,convection=False):
		# Liquid Gas Inlet
		# Solves Vaporization via D2-law
		# Specify the following inputs:
		# T0: Initial Droplet Temperature
		# gas: Resulting Gas Phase species in reaction mechanism
		# D0: Initial Droplet Diameter
		# mdot: Liquid Mass Flow Rate
		# rho: Liquid Density
		# mu: Liquid Viscosity
		# LHV: Liquid Latent Heat of Vaporization
		# Tboil: Boiling Temperature
		# Convection: Solve with Droplet Convection True/False

		# Error Check Input
		if T is None:
			raise ValueError('Liquid Injection Initial Temperature "T0" Not Entered')
		elif gas is None:
			raise ValueError('Gas Species of Corresponding Liquid Injection "gas" Not Entered')
		elif D0 is None:
			raise ValueError('Initial Droplet Size of Liquid Injection "D0" Not Entered')
		elif mdot is None:
			raise ValueError('Liquid Injection Mass Flow Rate "mdot" Not Entered')
		elif rho is None:
			raise ValueError('Liquid Injection Density "rho" Not Entered')
		elif mu is None:
			raise ValueError('Liquid Injection Viscosity "mu" Not Entered')
		elif LHV is None:
			raise ValueError('Liquid Injection Latent Heat of Vaporization "LHV" Not Entered')
		elif Tboil is None:
			raise ValueError('Liquid Injection Boiling Temperature "Tboil" Not Entered')
				
		# Initialize Droplet Class
		self.droplet = droplet(rho,mu,LHV,Tboil, reaction_mechanism = self.reaction_mechanism)
		# Initial droplet and mass
		self.droplet.initialize(D0,mdot,convection=convection)
		
		# Create liquid injection corresponding gas object
		self.inlet_liq = ct.Solution(self.reaction_mechanism)
		self.inlet_liq.TPY =  T, 101325, gas # Pressure is not used for anything in calculation
		# Store liquid mdot, D0, gas properties
		self.inlet_liq_mdot = mdot
		self.inlet_liq_D0 = D0
				
		# Enable droplet solve
		self.solve_droplet = True 


	def set_moments(self,M,P=[]):

		# Initialize MOMIC Class
		if self.MOMIC is None:
			self.MOMIC = MOMIC(self.M,self.P)
		else:
			self.MOMIC.set_moments(M,P=self.P)

		self.MOMIC.set_ct_gas(self.gas)
		
		
	def set_inlet_Moments(self,M,P=None):
		
		# Weird bug changes inlet Moment value
		# so using Copy to fix it  
		self.inlet_M = np.copy(M)
		self.inlet_P = []
		
		if P is not None:
			self.inlet_P = np.copy(P)
	
	def set_reactor_properties(self,vol=1,Ac=None,Twall=None):

		print('Setting Reactor Volume to: {} m^3'.format(vol))
		self.vol = vol

		if Ac is not None:
			print('Setting Reactor Cross Sectional Area to: {} m^2'.format(Ac))
			self.Ac = Ac
		
		if Twall is not None:
			print('Setting Reactor Wall Temperature to: {} K'.format(Twall))
			self.Twall = Twall

	def friction_factor(self,Re):
		
		A = -2*np.log10(12/Re)
		B = -2*np.log10(2.51*A/Re)
		C = -2*np.log10(2.51*B/Re)
		f = (A-(B-A)**2/(C-2*B+A))**(-2)
		
		return f

	def tau_wall(self,u,D):
		# Calculate Wall Shear Stress
		
		# Calculate Reynolds Number
		Re = self.gas.density*u*D/self.gas.viscosity
		# Calculate Friction Factor
		f = self.friction_factor(Re)
		
		# Calculate Shear Stress
		return self.gas.density*u**2/2*f
		
	
	def total_prop(self,gas,u):
		
		gamma = gas.cp/gas.cv
		R = 8314/gas.mean_molecular_weight
		T = gas.T
		P = gas.P
		
		M = u/(np.sqrt(gamma*R*T))
	
		Pt = P*(1+(gamma-1)/2*M**2)**((gamma)/(gamma-1))
		Tt = T*(1+(gamma-1)/2*M**2)
	
	
		return Pt,Tt
	
	
	def mix_pfr_inlet(self,inlets):
		
		n_inlets = len(inlets)
		
		#gas = [0]*n_inlets
		mix = [0]*n_inlets
		Q = 0 
		
		for i,inlet in enumerate(inlets):
			gas = ct.Solution(self.reaction_mechanism)
			gas.TPY = inlet[0],inlet[1],inlet[2]
			mix[i] = ct.Quantity(gas,mass=inlet[3],constant='HP')
			Q+= inlet[4]*inlet[3]
		
		return np.sum(mix),Q
	
	def ct_WSR_gas(self,vol):
		# Solve WSR using Cantera Models
		# For Solving with gaseous inlet only

		# Create Inlet/Outlet Reservoirs
		inlet_res = ct.Reservoir(self.inlet_gas)
		outlet_res = ct.Reservoir(self.outlet_gas)

		# Create Reactor and equilibrate inlet gas
		self.gas.TPY = self.inlet_gas.TPY
		self.gas.equilibrate('HP')
		reactor = ct.IdealGasReactor(self.gas)
		reactor.volume = vol

		# Create Inlet MFC
		inlet_mfc = ct.MassFlowController(inlet_res,reactor,mdot=self.inlet_gas_mdot)
							
		# Outlet Mass Flow Rate Function
		outletf = lambda dP: self.mdot_out_func()
		# Create Outlet MFC with outlet function
		outlet = ct.Valve(reactor,outlet_res,K=outletf)

		# Create ct Reactor Network
		sim = ct.ReactorNet([reactor])
		# Advance to Steady State
		sim.advance_to_steady_state()
		
		# WSR gas solution is stored in self.gas


	def ct_WSR(self,vol,mdot,Treactor=None,Q=0,Pc=None,Equil_IC=True):
		# Solve WSR properties

		# Define inlet WSR
		#inlet_gas = ct.Solution(self.reaction_mechanism)
		#inlet_gas.TPY = self.inlet_gas.T, self.inlet_gas.P, self.inlet_gas.Y
		
		if self.solve_droplet:
			#self.WSR_droplet_SS_vol(vol)
			print('CLEAN UP WSR!')
			keyboard()
			return

			Y_gas = self.inlet_gas.Y
			Y_liq = self.droplet.gas.Y
			Y_mix = (Y_gas*self.inlet_gas_mdot + Y_liq*self.liquid_mdot)/ \
				(self.inlet_gas_mdot+self.liquid_mdot)
		
			inlet_gas.TPY = self.inlet_gas.T, self.inlet_gas.P, Y_mix
			mdot+= self.vaporized_mdot
			

		if Equil_IC or self.solve_droplet:
			self.gas.TPY = self.inlet_gas.TPY #inlet_gas.TPY
			self.gas.equilibrate('HP')
			self.gas.TPY = self.gas.T, self.inlet_gas.P,self.gas.Y

		if Pc is not None:
			self.gas.TPY = self.gas.T, Pc, self.gas.Y
		
		if self.gas.P<self.outlet_gas.P:
			self.gas.TPY = self.gas.T, self.outlet_gas.P, self.gas.Y
		
		# Create Inlet Reservoir 
		inlet = ct.Reservoir(self.inlet_gas)
		reactor = ct.IdealGasReactor(self.gas)
		reactor.volume = vol
		
		if Treactor is not None:
			self.gas.TP = Treactor,self.gas.P
			reactor.insert(self.gas)
			reactor.energy_enabled = False
		#else:
		#	# Create Ideal Gas Reactor with equilibrium detonation products
		#	self.gas.TPY = self.inlet_gas.T, self.inlet_gas.P,self.inlet_gas.Y
		#	self.gas.equilibrate('UV') 

		self.inlet_gas_mdot = mdot
		# Exhaust
		exhaust = ct.Reservoir(self.gas)

		inlet_mfc = ct.MassFlowController(inlet,reactor,mdot=self.inlet_gas_mdot)
							
		outletf = lambda dP: self.mdot_out_func()
		
		outlet = ct.Valve(reactor,exhaust,K=outletf)#valve_K)
		#outlet = ct.PressureController(reactor,exhaust,master=inlet_mfc,K=0.01)

		if self.solve_droplet or Q != 0:
			Q+= self.droplet.LHV*self.liquid_mdot
			Wall=ct.Wall(reactor,exhaust,Q=Q)#,A=1)
		
		# Simulation
		sim = ct.ReactorNet([reactor])
		# Get SS Solution
		sim.advance_to_steady_state()
		# Calculate Gas Residence Time 
		self.tres = self.gas.density*vol/outletf(1)
		self.mdot = mdot
		self.vol = vol

	def vap_eq_ratio(self,per_vap,return_stoich_per=False):
		# Calculate equivalence ratio for a certain percent
		# of vaporized propellants
		# Can also determine what vaporized percent corresponds
		# to stoichiometric 

		Y_gas = self.inlet_gas.Y
		Y_liq = self.inlet_liq.Y
		
		mdot_gas = self.inlet_gas_mdot
		mdot_vap = per_vap*self.inlet_liq_mdot
		Y_mix = (Y_gas*mdot_gas + Y_liq*mdot_vap)/ (mdot_gas+mdot_vap)
		self.gas.Y = Y_mix

		if return_stoich_per:
			vmax = 1
			vmin = 0

			for run_counter in range(0,100):
				mid = (vmax+vmin)/2

				eq = self.vap_eq_ratio(mid)

				if eq<1:
					vmin = mid
				else:
					vmax = mid
				
				if np.abs(vmax-vmin)<.001:
					return mid
		else:
			return self.gas.get_equivalence_ratio()


	def ct_WSR_droplet(self,vol,per_vap,per_vap_previous=0,count=0,debug=False):
		# Solve WSR properties with droplet

		# Define inlet gas, with droplet properties
		Y_gas = self.inlet_gas.Y
		Y_liq = self.inlet_liq.Y
		Tgas = self.inlet_gas.T
		Tliq = self.inlet_liq.T

		mdot_gas = self.inlet_gas_mdot
		mdot_vap = per_vap*self.inlet_liq_mdot
		
		Y_mix = (Y_gas*mdot_gas + Y_liq*mdot_vap)/ \
				(mdot_gas+mdot_vap)

		Tmix = (Tgas*mdot_gas+Tliq*mdot_vap)/ \
			(mdot_gas+mdot_vap)

		# Create Gas Mixture Object
		self.gas.TPY = Tmix,self.inlet_gas.P,Y_mix

		# Create Inlet/Outlet Reservoirs
		inlet_res = ct.Reservoir(self.gas)
		outlet_res = ct.Reservoir(self.outlet_gas)

		# Create Reactor and equilibrate inlet gas
		self.gas.equilibrate('HP')
		reactor = ct.IdealGasReactor(self.gas)
		reactor.volume = vol
	
		# Create Inlet MFC
		inlet_mfc = ct.MassFlowController(inlet_res,reactor,mdot=mdot_gas+mdot_vap)

		# Outlet Mass Flow Rate Function
		outletf = lambda dP: self.mdot_out_func()
		# Create Outlet MFC with outlet function
		outlet = ct.Valve(reactor,outlet_res,K=outletf)


		# Vaporization enthalpy
		Q= self.droplet.LHV*mdot_vap
		Wall=ct.Wall(reactor,outlet_res,Q=Q)

		# Create ct Reactor Network
		sim = ct.ReactorNet([reactor])
		
		# Advance to Steady State				
		try:
			sim.advance_to_steady_state()
		except:
			if self.gas.T<1.5*Tmix:
				return 0
			elif self.gas.T>9000:
				return 0
			else:
				keyboard()

		self.tres = self.gas.density*reactor.volume/(mdot_gas+mdot_vap)
		self.U = (mdot_gas+mdot_vap)/(self.gas.density*self.Ac)


		# Calculate droplet vaporization percent
		_,_,per_vap2 = self.droplet.ct_gas(self.gas,self.tres,U=self.U)
		

		# Recursion to find % vaporized for this reactor volume
		if np.abs(per_vap-per_vap2) <0.001:
			return per_vap
		elif count==100:
			print('Recursion is taking a long time, maybe something wrong?')
			keyboard()
		elif count==-1:
			return per_vap2
		else:
			if self.gas.T<1.5*Tmix:
				if per_vap > per_vap_previous:
					return 0
				else:
					return -1

			per_vap_new = (per_vap2+per_vap)/2
			per_vap = self.ct_WSR_droplet(vol,per_vap_new,per_vap_previous=per_vap,\
			count=count+1)
			return per_vap


	def WSR_tres(self,tres,Treactor=None):
		# Calculate WSR properties based on residence time
		
		# Create Inlet Reservoir 
		inlet = ct.Reservoir(self.inlet_gas)
		
		
		# Create Ideal Gas Reactor with inlet equilibrium products
		self.gas.TPY = self.inlet_gas.TPY
	
		if Treactor is not None:
			self.gas.TP = Treactor,self.gas.P
		else:
			# Create Ideal Gas Reactor with equilibrium products
			self.gas.equilibrate('HP') 
		
		reactor = ct.IdealGasReactor(self.gas)
		
		if Treactor is not None:
			reactor.energy_enabled = False
		
		
		# Exhaust
		exhaust = ct.Reservoir(self.gas)
		
		mdot = lambda t: self.gas.density*reactor.volume/tres
		
		inlet_mfc = ct.MassFlowController(inlet,reactor,mdot=mdot)
		outlet = ct.PressureController(reactor,exhaust,master=inlet_mfc,K=.01)
		
		# Simulation
		sim = ct.ReactorNet([reactor])
		
		sim.advance_to_steady_state()
		
		# Calculate Gas Residence Time 
		self.tres = self.gas.density*reactor.volume/inlet_mfc.mdot(sim.time)
		self.mdot = mdot(sim.time)
		self.vol = reactor.volume
	
		
	def equilibrium_gas(self):
		self.gas.equilibrate('HP')
	
	def concentration(self,specie):
		# Get Cantera Gas Concentration
		
		# Convert from kmol/m^3 to mol/cm^3
		return self.gas.concentrations[self.gas.species_index(specie)]*1e-3
	
	def Newton_Step(self):
		# OLD
		# Newton Iteration Step to find SS Moment Solution
		# USED IN OLD SOLUTION APPROACH BUT MAYBE SHOULD LOOK BACK INTO IT
		# TO SPEED UP CONVERGENCE WITH NEW APPROACH
		
		# Get dlogM and dlogP rates
		dMdt, dPdt,_,_ = self.dlog_rates()
		
		Moments = np.append(self.M,self.P)
		# Calculate Jacobian
		J = momic.momic.jacobian_log(Moments,np.size(self.M),np.size(self.P), \
				self.gas.T, self.gas.P, 		\
				self.concentration('C2H2'), 			\
				self.concentration('H'), 				\
				self.concentration('H2'),				\
				self.concentration('H2O'), 				\
				self.concentration('O2'), 				\
				self.concentration('OH') )
		
		#if self.debug_flag:
		#	keyboard()
			
		if self.aggregation:
		
			if np.all(J[(np.size(self.M)+1):,(np.size(self.M)+1):]==0):
				dS = -np.dot(np.linalg.inv(J[0:np.size(self.M),0:np.size(self.M)]),dMdt)
				dS = np.hstack([dS,np.zeros(np.size(self.P))])
			else:
				dS = -np.dot(np.linalg.inv(J),np.append(dMdt,dPdt))
			
			#dS_P = dS[np.size(self.P):]
			#lam_P = self.calc_lambda(np.log(self.P), dS_P)[0]
		else:
			dS = -np.dot(np.linalg.inv(J),dMdt)
			lam_P = 9999
		
		dS_M = dS[0:np.size(self.M)]
		# Calculate minimum step
		lam_M = self.calc_lambda(np.log(self.M), dS_M)[0]

		#lam = min([lam_M,lam_P])
		lam = lam_M

		# Limit maximu change in M to < 2
			
		lam = np.min([lam,5/np.max(np.abs(dS_M))])

		
		#if self.debug_flag:
		#	keyboard()

		old_rate = np.average(np.abs(self.dlog_rates()[0]))
		M0 = self.M		
		self.M = self.M*np.exp(lam*dS_M/1e5)
		new_rate = np.average(np.abs(self.dlog_rates()[0]))
	
		j_iter = 0
		while old_rate < new_rate:
			j_iter+= 1
			lam /= 1e2
			self.M = M0*np.exp(lam*dS_M/1e5)
			new_rate = np.average(np.abs(self.dlog_rates()[0]))
			
			if j_iter== 10:
				self.fail+=1
				print('Newton Iteration Failed, Using Previous Solution')
				self.M = M0
				break

		#if self.aggregation:
		#	self.P = self.P*np.exp(lam*dS_P/1e5)
		#	self.P[0] = self.M[0]
	

	
	def mdot_out_func(self,mdot_set_pt=0,K=1e-4):
		
		# Same as Cantera Valve Object
		if self.outlet_BC == 'P':
			# Can result in negative transient outlet mass flow rate 
			return (self.gas.P-self.outlet_gas.P)
			# This implementation below makes the solver blow up sometimes:
			#return np.max([(self.gas.P-self.outlet_gas.P),0])

		# Same as Cantera Pressure Controller
		if self.outlet_BC == 'PC':
			return mdot_set_pt + K*(self.gas.P-self.outlet_gas.P)

		# Assume flow is choked at all times
		if self.outlet_BC == 'mass':
			return self.mdot_choked()
		
		# Calculate based on Compressible
		if self.outlet_BC == 'compressible':
			return self.mdot_compressible(self.outlet_gas.P)
	
	
	def set_outlet(self,outlet_BC_type, P=101325, At=1):
		
		self.outlet_BC = outlet_BC_type
		self.outlet_gas.TPY = 300, P, 'N2:1'
		self.At = At
			
	def mdot_choked(self):#,Dt):
		
		Pt = self.gas.P
		Tt = self.gas.T
		gamma = self.gas.cp/self.gas.cv
		R = 8314/self.gas.mean_molecular_weight
		Ma = 1
		
		return self.At*Pt*np.sqrt(gamma/(R*Tt))*Ma*(1+(gamma-1)/2*Ma**2)**(-(gamma+1)/(2*(gamma-1)))
	
	def mdot_compressible(self,Pe):
		
		Pt = self.gas.P
		Tt = self.gas.T
		gamma = self.gas.cp/self.gas.cv
		R = 8314/self.gas.mean_molecular_weight
		
		if Pt<Pe:
			return 0
		
		M = np.sqrt(2/(gamma-1)*((Pt/Pe)**((gamma-1)/gamma)-1))
		
		if M > 1:
			M = 1
		
		return self.At*Pt*np.sqrt(gamma/(R*Tt))*M*(1+(gamma-1)/2*M**2)**(-(gamma+1)/(2*(gamma-1)))

	@staticmethod
	def mdot_compressible_gas(gas,Pe,At):
		# Static Method to calculate gas mass flow rate
		# Inputs are:
		# CT Gas Object, Static Pressure (in Pa), Exit Area (in m^2)

		Pt = gas.P
		Tt = gas.T
		gamma = gas.cp/gas.cv
		R = 8314/gas.mean_molecular_weight

		# If stagnation pressure is below static pressure, then there is no forward flow
		# (but probably backflow, but I'm not calculating that, so outputting 0 for mdot)
		if Pt<Pe:
			return 0

		# Calculate Mach Number of gas
		M = np.sqrt(2/(gamma-1)*((Pt/Pe)**((gamma-1)/gamma)-1))

		# In a non-CD flow area, the max flow velocity is sonic
		if M > 1:
			M = 1

		# Calculate and return compressible mass flow rate
		# Formulation copied from NASA Glenn Compressible Mass Flow page
		return At*Pt*np.sqrt(gamma/(R*Tt))*M*(1+(gamma-1)/2*M**2)**(-(gamma+1)/(2*(gamma-1)))
						
	
	def mdot_out(self,gas):
	
		At = (.308*.0254)**2/4*np.pi
		Pt = gas.P
		Tt = gas.T
		gamma = gas.cp/gas.cv
		R = 8314/gas.mean_molecular_weight
		Ma = 1
		
		mout = At*Pt*np.sqrt(gamma/(R*Tt))*Ma*(1+(gamma-1)/2*Ma**2)**(-(gamma+1)/(2*(gamma-1)))
		
		return mout	
	
	def blowout(self,precision):
		# Find minimum tres to prevent blowout
	
		# Starting Bounds
		tres_low = 1e-8
		tres_hi = 1e-3
		
		# T0
		T0 = self.inlet_gas.T
		
		# Initialize search
		
		# If tres hi bound results in blowout, increase hi-bound and rerun
		while True:
			self.WSR_tres(tres_hi)
			
			if self.gas.T <1.5*T0:
				tres_hi*=1e1
			else:
				break
		
		# If tres low bound results in ignition, decrease lo-bound and rerun
		while True:
			self.WSR_tres(tres_low)
			
			if self.gas.T >1.5*T0:
				tres_low*=1e-1
			else:
				break
		
		# Run bisection solver for 101 steps
		for run_counter in range(0,1000):
		
			tres_mid = (tres_low+tres_hi)/2
			
			self.WSR_tres(tres_mid)
			
			if self.gas.T > 1.5*T0:
				tres_hi = tres_mid
			else:
				tres_low = tres_mid

			if (tres_hi-tres_low)<precision:
				break
		
		if self.gas.T>1.5*T0:
			return tres_mid
		else:
			self.WSR_tres(tres_hi)
			return tres_hi


	def blowout_V(self,precision,Pc=None):
		# Find minimum tres to prevent blowout
	
		# Starting Bounds
		V_hi = 1e-3
		
		# T0
		T0 = self.inlet_gas.T
		
		WSR_func = self.WSR
		# Initialize search
		
		# If tres hi bound results in blowout, increase hi-bound and rerun
		while True:
			WSR_func(V_hi,self.inlet_gas_mdot,Pc=Pc,Equil_IC=True)

			if self.gas.T >2*T0:
				V_hi*=1e-1
			else:
				V_low = np.copy(V_hi)
				V_hi*=1e1
				break
		
		# Run bisection solver for 1001 steps
		for run_counter in range(0,1000):
		
			V_mid = (V_low+V_hi)/2
			
			WSR_func(V_mid,self.inlet_gas_mdot,\
				Pc=Pc,Equil_IC=True)

			if self.gas.T > 1.5*T0:
				V_hi = V_mid
			else:
				V_low = V_mid

			if (V_hi-V_low)<precision:
				break
		

		if self.gas.T>2*T0:
			return V_mid
		else:
			WSR_func(V_hi,self.mdot,Pc=Pc)
			return V_hi

	def blowout_vol(self,precision,starting_vol=1e-3):
		# Find minimum volume for WSR to not blowout
		# Setup WSR model first before running this

		# T0 to determine blowout
		T0 = self.inlet_gas.T

		# WSR function to use to calculate blowout
		if self.solve_droplet:
			WSR_func = self.ct_WSR_droplet
		else:
			WSR_func = self.ct_WSR_gas


		# Starting hi Bounds
		V_min = starting_vol
		V_hi = 0
		V_lo = 0

		# Find blowout bounds
		hi_bound_found = False
		lo_bound_found = False
		while True:
			WSR_func(V_min)

			if self.gas.T >1.5*T0:
				V_hi = V_min
				hi_bound_found = True
				V_min*= 1e-1
			else:
				V_lo = V_min
				lo_bound_found = True
				V_min*=1e1

			# Once bounds are found, break
			if hi_bound_found and lo_bound_found:
				break
		
		# Now that blowout bounds are bound, lets use bisection
		# because it's easy to implement
		# 		
		# Run bisection solver for 1000 steps
		for run_counter in range(0,1000):
			
			V_mid = (V_lo+V_hi)/2


			WSR_func(V_mid)

			if self.gas.T > 1.5*T0:
				V_hi = V_mid
			else:
				V_lo = V_mid

			if (V_hi-V_lo)<precision:
				break
		

		if self.gas.T>1.5*T0:
			return V_mid
		else:
			WSR_func(V_hi)
			return V_hi
	
	def blowout_vol_drop(self,precision,starting_vol=1e-3,vap_target=.25):
		# Find minimum volume for WSR to not blowout with liquid injection
		# Setup WSR model first before running this


		WSR_func = self.ct_WSR_droplet

		# Find droplet vaperization percent for stoichiometric mixture 
		vap_per_stoich = self.vap_eq_ratio(.1,return_stoich_per=True)
		
		# Check to see if WSR does not blowout at stoich proporations
		# (Otherwise the WSR is too small)
		# Find volume that results in complete droplet vaporization
		vmax = starting_vol
		vmin = 0
		vap_per = WSR_func(vmax,vap_per_stoich,count=-1)

		for run_counter in range(0,10):
			
			if vap_per >0.995:
				break
			elif run_counter==9:
				keyboard()
			else:
				vmax*=1e1
			
			vap_per = WSR_func(vmax,vap_per_stoich,count=-1)
		
		# Flag for controlling search finding
		partially_vaporized_flag = False
		# Run Bisection to find minimum WSR volume
		for run_counter in range(0,100):
			
			# Midpoint Volume
			vmid = (vmin+vmax)/2
			# Calculate Vaporization percent without droplet vaporization
			vap_per = WSR_func(vmid,vap_per_stoich,count=-1)
			
			# If the droplet is partially vaporized then:
			if vap_per<1 and vap_per>0:
				# Calculate Droplet Vaporization coupled to WSR
				vap_per = WSR_func(vmid,vap_per_stoich)
				
				# If the coupled solution predicts no droplet vaporization
				# then adjust bounds and continue
				if vap_per in [0,-1]:
					# If near partially vaporized area and 
					# no droplet vaporization/blowout detected
					# then this is a low bound
					if partially_vaporized_flag:
						vmin = vmid
					else:
						vmax = vmid
						continue
				# Otherwise if the droplet is partially vaporized,
				# Then store the vap percent and volume of partially vaporized
				# Also we're in the partially vaporized regime so turn on flag
				else:
					partially_vaporized_flag = True
					vap_per_stoich = vap_per
					vstore = vmid
					vmax = vmid

			elif vap_per !=0:
				vmax = vmid
			else:
				vmin = vmid


			print('Volume:{}'.format(vmid))
			print('%Vap:{}'.format(vap_per))
			if vap_per<= vap_target:
				return vmid
			if np.abs(vmax-vmin) <1e-8:
				if vap_per == -1:
					return vmax
				else:
					return vmid

	

	def blowout_droplet(self,precision,Pc=None):
		# Find minimum WSR volume to prevent blowout
		# With Droplet calculation
		
		# Starting Bounds
		tres_hi = 1e-3

		# Inlet Temperature T0
		T0 = self.inlet_gas.T
		
		WSR_func = self.WSR_droplet_SS

		# Check if high starting bound produces ignition, if not increase
		# then keep decreasing bounds by 1e-1 until there is no ignition
		# This makes bisection converge faster since bounds are smaller

		break_loop = False

		while True:
			WSR_func(tres_hi,Pc=Pc)
			if self.gas.T>2*T0:
				tres_hi*=1e-1
				break_loop = True
			else:
				tres_low = np.copy(tres_hi)
				tres_hi*=1e1
				if break_loop:
					break

		# Run bisection solver for 1001 steps
		for run_counter in range(0,1000):
		
			tres_mid = (tres_low+tres_hi)/2
			
			WSR_func(tres_mid,Pc=Pc)

			if self.gas.T > 1.5*T0:
				tres_hi = tres_mid
			else:
				tres_low = tres_mid

			if np.abs(tres_hi-tres_low)/tres_mid<precision:
				break
		
		if self.gas.T>2*T0:
			tres_min = tres_mid
		else:
			WSR_func(tres_hi,Pc=Pc)
			tres_min = tres_hi

		# Output minimum WSR volume
		Vmin = tres_min/self.gas.density*(self.vaporized_mdot+self.inlet_gas_mdot)
		
		return Vmin

	def initialize(self):

		# Update Control Flags
		self.update_control_flags()

		# Update References
		self.M = self.MOMIC.M
		self.P = self.MOMIC.P
		self.M_Moments = self.MOMIC.Num_M_Moments
		self.P_Moments = self.MOMIC.Num_P_Moments

		# Check for Errors
		self.error_check()

	def solve(self,dt,tf,convergence_criteria=1e-5,ct_WSR_IC=True):
		
		# Initialize Model
		self.initialize()

		# Integrate 
		self.integrate(dt,tf,convergence_criteria=convergence_criteria,ct_WSR_IC=ct_WSR_IC)

	def update_control_flags(self):
		
		# Update Control Flags from Main Class
		if self.main is not None:
			self.solve_soot = self.main.solve_soot
			self.solve_coupled = self.main.solve_coupled
			self.solve_energy = self.main.solve_energy 

		if self.inlet_liq:
			self.solve_droplet = True
		else:
			self.solve_droplet = False

	def error_check(self):
		# Check for Errors

		if self.solve_soot:
			if np.size(self.inlet_M) != self.MOMIC.Num_M_Moments:
				raise ValueError('The Number of Inlet M-Moments does not match number of initialized Moments')
			if self.MOMIC.aggregation:
				if np.size(self.inlet_P) != self.MOMIC.Num_P_Moments-1:
					raise ValueError('The Number of Inlet P-Moments does not match number of initialized Moments')
			
		if self.vol is None:
			raise ValueError('Reactor Volume was not entered')

		if self.inlet_gas is None:
			raise ValueError('No Inlet Gas Was Defined')
		if self.outlet_gas is None:
			raise ValueError('No Outlet Was Defined')

		if self.solve_droplet:
			# This is probably not necessary anymore because
			# not possible to reach this condition?
			if self.inlet_liq is None:
				raise ValueError('No Droplet Gas was Defined')
			if all(self.inlet_liq.Y == self.inlet_gas.Y):
				raise ValueError('Liquid and Gas Injection Compositions are the Same')

	def integrate(self,dt,tf,convergence_criteria=None,ct_WSR_IC=False):
		# Main Function controlling integration of WSR equations
		# Integrate WSR equations as coupled or decoupled
		# Run a Cantera WSR to SS before starting the calculation
		#keyboard()

		# Initialize Arrays
		self.mdot_vap = 0

		if ct_WSR_IC:
			self.gas.TPY = self.inlet_gas.TPY
			self.gas.equilibrate('HP')
			self.gas.TPY = self.gas.T, self.inlet_gas.P,self.gas.Y
			self.ct_WSR(self.vol,self.inlet_gas_mdot)
		print('CLEAN UP integrate_WSR!!!')

		
		# Run if solve_droplet is on, in order to help with convergence, solve
		# for gas properties with no soot calculation 
		if self.solve_droplet:
			# Run Droplet Calculations if there are still droplets in flow
			if self.droplet.D0 >0:
				# Run Droplet Calculation with ct formulation
				if ct_WSR_IC:
					pass
				else:
					print('Cleanup in integrate_WSr')
					#self.WSR_droplet_SS_vol(self.vol)
			else:
				self.solve_droplet = False

		if self.solve_droplet:
			# Run Droplet Calculate with soot source turned off
			#self.integrate_WSR_droplet(1e-3,2)
			self.ct_WSR_droplet(self.vol,.1)
		
		# Check IC 
		if self.gas.T<500:
			print('IC for Gas Temperature too low')
			if self.solve_soot:
				self.MOMIC.ct_rate()
				self.mdot_eat = 0
			return

		# Create IC and define corresponding indices for data
		y0 = np.hstack([self.gas.T,self.gas.density,self.gas.Y])
		self.indx = {'T':0,'D':1,'Y':np.arange(2,self.gas.n_species+2)}

		if self.solve_coupled:
			sol_name = self.rates_soot_coupled
			self.indx.update({'M':np.arange(self.indx['Y'][-1]+1,self.indx['Y'][-1]+1+self.M_Moments)})
			self.indx.update({'P':np.arange(self.indx['M'][-1]+1,self.indx['M'][-1]+1+self.P_Moments-1)})
			self.indx.update({'mp':self.indx['M'][-1]+self.P_Moments})
			self.indx.update({'Tp':self.indx['mp']+1})
		else:
			# If Uncoupled, first solve WSR for gas propeerties
			# then solve for MOMIC rates
			sol_name = self.rates_gas

			# Save Original solve_soot flag and turn off
			# soot calculation for decoupled case
			solve_soot = np.copy(self.solve_soot)
			self.solve_soot = False
			self.scipy_solver(sol_name,y0,dt,tf,convergence_criteria)
			self.solve_soot = solve_soot
			
			if self.solve_soot:
				sol_name = self.rates_soot_decoupled
				y0 = []
				self.indx = {'M':np.arange(0,self.M_Moments)}
				self.indx.update({'P':np.arange(self.indx['M'][-1]+1,self.indx['M'][-1]+1+self.P_Moments-1)})
				self.indx.update({'mp':self.indx['M'][-1]+self.P_Moments+1})
				self.indx.update({'Tp':self.indx['mp']+1})
			else:
				return
						
		# If solving for soot then add extra Eq
		if self.solve_soot:
			y0 = np.hstack([y0,np.log(self.M),np.log(self.P),0,self.gas.T])
		
		# Add index for storing time points
		self.indx.update({'t':self.indx['Tp']+1})

		#keyboard()
		#sol_name = self.source_terms_coupled
		# Run Scipy ODE Solver
		self.scipy_solver(sol_name,y0,dt,tf,convergence_criteria)


	def scipy_solver(self,sol_name,y0,dt,tf,convergence_criteria=False):
		# Solve ODE 
		
		#if self.solve_droplet:
		#	self.droplet.ct_gas(self.gas,self.tres)
		
		## Set Up ODE
		# Defining source term function and numerical Jacobian function
		solver = scipy.integrate.ode(sol_name)
		#solver.set_integrator('lsoda',method='bdf', with_jacobian=True)
		solver.set_integrator('vode',method='bdf', with_jacobian=True)

		# Set initial Conditions
		solver.set_initial_value(y0, 0.0)
		
		# Total number of time points to save data
		n_pts = int(np.ceil(tf/dt))
		
		# Preallocate Solution Array and assign IC
		self.Sol = np.zeros([n_pts+1,np.size(y0)+1])

		self.Sol[0,:] = np.append(solver.y,0)
		
		# Initialize some variables
		rates = np.zeros(np.size(solver.y))
		Y0 = np.zeros(np.size(solver.y))
		self.time = 0
		count = 1
		i = 0
		
		self.sources_constant()
		#keyboard()
		while solver.t < tf:
			
			# Save Previous Solution
			# This is copied over from old code that was unstable and required reinitailizing
			# decreasing time step if diverge occurs, or if time step is too large. I don't
			# think I need it anymore but saving in case I do
			
			Y0[:] = solver.y
			t0 = solver.t

			# Moved Droplet Calculation outside of source term,
			# so integrating each time step with constant droplet
			# source, I was having solver problems otherwise
			self.sources_constant()				

			# If Soot conditions changed then reinitialize solver
			# it's 3:47 am and my brain shut down so I'm doing this dumb way to reinitialze solver
			# But when you wake up, read through scipy documentation to see how to reinitialize
			if self.MOMIC.condition_switch_flag:
				solver = scipy.integrate.ode(sol_name)
				solver.set_integrator('vode',method='bdf', with_jacobian=True)
				solver.set_initial_value(Y0, t0)

			#keyboard()
			# Integrate in time	

			solver.integrate(solver.t+dt)

			solver_iter = 0
			

			# Sometimes solver can't advance full time step
			# so keep integrating or pause and give terminal
			# access if progress is too slow
			while solver.t<t0+dt:
				solver.integrate(solver.t+dt)
				solver_iter+= 1
				if solver_iter == 100:
					keyboard()
			
			# Rates at current time step
			rates[:] = solver.y#sol_name(solver.t,solver.y)
			
			# Get Solution at current time step
			self.time = solver.t
			self.rates = rates			
			self.state = np.hstack([solver.y,\
									self.time])


			# Calculate flow velocity if Ac is specified
			if self.Ac is not None:
				self.U = self.mdot_gas/(self.gas.density*self.Ac)

			# Check Solution and reduce time step (for n time steps) 
			# if any criteria is met
			if any(np.isnan(solver.y)):
				print('y is nan')
				#solver.set_initial_value(Y0, t0)
				#dt*=.1
				#print('Decreasing time step to:{}'.format(str(dt)))
				#count=-1*dt_original/dt
				keyboard()
			
			elif 'Y' in self.indx.keys():
				if any(solver.y[self.indx['Y']]>1):
					#print('yi>1')
					#solver.set_initial_value(Y0, t0)
					#dt*=.1
					#print('Decreasing time step to:{}'.format(str(dt)))
					#count=-1*dt_original/dt
					pass
					#keyboard()
			
			# Return time step to original dt after integrating for 10 decreased
			# time steps
			#if count==0:
			#	dt = dt_original
			#	print('Reseting time step back to original stepsize')
			
			# Save Solution
			if count>0:
				try:
					self.Sol[i,:] = self.state
				except:
					keyboard()
				i+=1
				
			# Check if Convergence Criteria Specified
			if convergence_criteria:
				if all(np.abs(rates)<convergence_criteria):
					print('Convergence Criteria Met')

					# Return nonzero Solution Array
					self.Sol = self.Sol[~(self.Sol==0).all(axis=1)]
					return

		# The integration time should be plenty long for gas to reach SS, sometimes
		# moments can have time reaching steady, check if their rates are small,
		# if they are then we can assume it's converged I think?

		if convergence_criteria:
			if self.solve_soot:
				if all(np.abs(self.soot_rates()[0])<convergence_criteria):
					print('Convergence Probably Met?')
					# Return nonzero Solution Array
					self.Sol = self.Sol[~(self.Sol==0).all(axis=1)]
					return

		# Return nonzero Solution Array
		self.Sol = self.Sol[~(self.Sol==0).all(axis=1)]
							
	def rates_soot_coupled(self,t,y):
		
		# Check input for propblems
		# Sometimes scipy ode goes crazy
		# and enters stupid scalars which screws up
		# Cantera/Moment functions - so ignore these inputs
		if any(np.isnan(y)):
			pass
		elif y[self.indx['D']]<0:
			pass
		elif y[self.indx['T']]<100:
			pass
		elif y[self.indx['T']] > 5000:
			pass
		elif y[self.indx['D']] > 1000:
			pass
		else:
			# Assign solution to gas object if input looks good
			self.gas.TDY = y[self.indx['T']],y[self.indx['D']],y[self.indx['Y']]	
				
		# Check input moment for problems - same as above
		if self.solve_soot:
			if any(np.exp(y[self.indx['M']]) > 1e200):
				pass
			elif any(np.isnan(y[self.indx['M']])):
				pass
			elif any(np.diff(y[self.indx['M']]<0)):
				pass
			else:
				# Assign current solution to momic object attribute				
				self.M[:] = np.exp(y[self.indx['M']])
				self.P[:] = np.exp(y[self.indx['P']])
				self.mp = y[self.indx['mp']]
				self.Tp = y[self.indx['Tp']]
		
		# Calculate Source Terms
		self.sources()
		
		# WSR gas balance
		self.mdot_gas = self.inlet_gas_mdot + 0 + self.S_mass_t()		
		self.mdot_gas+= self.S_mass

		# Calculate gas residence time in WSR
		self.tres = self.gas.density*self.vol/self.mdot_gas
		
		# Gas Properties	
		u_gas = self.gas.partial_molar_int_energies
		Y_gas_in = self.inlet_gas.Y
		Y_gas = self.gas.Y
		MW = self.gas.molecular_weights

		rho = self.gas.density
		cv = self.gas.cv
		#Q_source = 0
		wdot = self.gas.net_production_rates
		
		mdot_out =  self.mdot_out_func(mdot_set_pt=self.mdot_gas)

		# Energy Equation
		if self.solve_energy:
			dTdt = self.inlet_gas_mdot*(self.inlet_gas.h - np.dot(u_gas/MW,Y_gas_in))-self.gas.P/rho*mdot_out-np.dot(u_gas,wdot)*self.vol #+ Q_source					
			dTdt += self.S_energy + self.S_energy_t() 		
			dTdt /= rho*self.vol*cv				
		else:
			dTdt = 0
		
		# Soot Equations
		if self.solve_soot:
			# Calculate Rates
			dTpdt = 0#(Q_source)#-0.9*5.67e-8*Ap*(self.Tp**4-300**4)+np.dot(u_gas,wdot_soot)*4/3*np.pi*(d_s/2)**3
			m_soot = (self.M[1]*1.9944235e-26)/(self.M[0])
			Cp_p = 840
			dTpdt *= 1/(m_soot*Cp_p)
			dmpdt = -self.mdot_eat - self.mp/self.tres
			
			dlogMdt,dlogPdt = self.soot_rates()

		# Continuity Equation			
		drhodt = (self.inlet_gas_mdot-mdot_out+ self.S_mass +self.S_mass_t())
		drhodt/= self.vol

		# Store rates
		self.drhodt = drhodt
		self.dTdt = dTdt
		self.mdot_out = mdot_out
		
		# Species Equations
		dYdt = wdot*MW/(rho) + (Y_gas_in-Y_gas)*self.inlet_gas_mdot/(rho*self.vol)
		dYdt+= (self.S_species)/(self.gas.density*self.vol)
		dYdt+= (self.S_species_t())/(self.gas.density*self.vol)

		# Rates output vector
		if self.solve_soot:
			rates = np.hstack([dTdt,drhodt,dYdt,dlogMdt,dlogPdt,dmpdt,dTpdt])
		else:
			rates = np.hstack([dTdt,drhodt,dYdt])
		
		
		return rates
		
	def rates_soot_decoupled(self,t,y):
		
		# Check Solution
		if any(np.isnan(y)):
			pass
		else:
			# Assign current solution to momic object attribute
			self.M[:] = np.exp(y[self.indx['M']])
			self.P[:] = np.exp(y[self.indx['P']])		
			self.mp = y[self.indx['mp']]
			self.Tp = y[self.indx['Tp']]

		# Run Soot calculation
		self.sources()

		# Calculate gas residence time in WSR
		self.tres = self.gas.density*self.vol/self.mdot_gas
			
		# Calculate Rates
		dTpdt = 0#(Q_source)#-0.9*5.67e-8*Ap*(self.Tp**4-300**4)+np.dot(u_gas,wdot_soot)*4/3*np.pi*(d_s/2)**3
		m_soot = (self.M[1]*1.9944235e-26)/(self.M[0])
		Cp_p = 840
		dTpdt *= 1/(m_soot*Cp_p)
		dmpdt = -self.mdot_eat - self.mp/self.tres
		
		# Soot Rates
		dlogMdt,dlogPdt = self.soot_rates()
		
		# Rates output vector
		rates = np.hstack([dlogMdt,dlogPdt,dmpdt,dTpdt])
		
		return rates
		
	def rates_gas(self,t,y):
		
		# Check input for propblems
		# Sometimes scipy ode goes crazy
		# and enters stupid scalars which screws up
		# Cantera/Moment functions - so ignore these inputs
		if any(np.isnan(y)):
			pass
		elif y[self.indx['D']]<0:
			pass
		elif y[self.indx['T']]<100:
			pass
		elif y[self.indx['T']] > 5000:
			pass
		elif y[self.indx['D']] > 1000:
			pass
		else:
			# Assign solution to gas object if input looks good
			self.gas.TDY = y[self.indx['T']],y[self.indx['D']],y[self.indx['Y']]	

		# Calculate total inlet mass flow rate (Gas + Vaporized Droplet)
		self.mdot_gas = self.inlet_gas_mdot + self.S_mass_t()
		
		# Calculate Source Terms (Wall/Soot/any interactions specified here)
		self.sources()
		
		# Calculate gas residence time in WSR
		self.tres = self.gas.density*self.vol/self.mdot_gas

		# Gas Info	
		u_gas = self.gas.partial_molar_int_energies
		Y_gas_in = self.inlet_gas.Y
		Y_gas = self.gas.Y
		MW = self.gas.molecular_weights

		rho = self.gas.density
		cv = self.gas.cv

		wdot = self.gas.net_production_rates
		
		mdot_out =  self.mdot_out_func(mdot_set_pt=self.mdot_gas) 

		# Energy Eq
		if self.solve_energy:
			dTdt = self.inlet_gas_mdot*(self.inlet_gas.h - np.dot(u_gas/MW,Y_gas_in))-self.gas.P/rho*mdot_out-np.dot(u_gas,wdot)*self.vol #+ Q_source					
			dTdt += self.S_energy +self.S_energy_t()
			dTdt /= rho*self.vol*cv				
		else:
			dTdt = 0
			
		# Continuity Eq
		drhodt = (self.inlet_gas_mdot-mdot_out+ self.S_mass + self.S_mass_t())
		drhodt/= self.vol

		self.drhodt = drhodt
		self.dTdt = dTdt
		self.mdot_out = mdot_out
		
		# Energy Equation
		dYdt = wdot*MW/(rho) + (Y_gas_in-Y_gas)*self.inlet_gas_mdot/(rho*self.vol)	
		dYdt+= (self.S_species + self.S_species_t())/(self.gas.density*self.vol)

		# Rates output vector
		rates = np.hstack([dTdt,drhodt,dYdt])

		return rates
		
	def soot_rates(self):	
		
		dlogM = self.MOMIC.SM/self.M - 1/self.tres + self.inlet_M/(self.M*self.tres)	
		dlogP = self.P

		if self.MOMIC.aggregation:
			dlogP = self.MOMIC.SP/self.P - 1/self.tres + self.inlet_P/(self.P*self.tres)
		else:
			dlogP = np.array(self.P)*0
		
		return [dlogM, dlogP]
	
	def sources(self):
	
		# Initialize Sources
		S_mass = 0 
		S_species = 0 
		S_energy = 0 
		
		# Droplet Source Term (Handled in external sources now)
		if self.solve_droplet and False:
			#[m_droplet,Df] = self.droplet.ct_gas(self.gas,self.tres)

			# Use Constant Source for each time step
			# Droplet is calculated in scipy_solver
			
			m_gas = self.droplet.m_gas*self.droplet.N_drops
			#m_gas= m_droplet*self.droplet.N_drops 
			#print('Using this for now in sources')
			#m_gas = 0.03
			
			####### Mass #########
			S_mass+= m_gas
			
			####### Species #########
			S_species+= m_gas*(self.inlet_liquid_gas.Y-self.gas.Y)
			
			####### Energy #########
		
			# Vaporization Enthalpy
			S_energy+= -m_gas*self.droplet.LHV
			
			# Vaporized gas energy 
			S_energy+= m_gas*(self.inlet_liquid_gas.h-\
				np.dot(self.gas.partial_molar_int_energies/self.gas.molecular_weights,self.inlet_liquid_gas.Y))
		
			# Balance of liquid/gas mass
			#self.mdot_gas+= S_mass
			#self.mdot_liq = self.liquid_mdot-m_gas
			# This should be decommisioned. Need to modify/delete
			keyboard()

		if self.solve_soot:
			
			# Run Soot calculation
			self.MOMIC.ct_rate()

			self.mdot_eat = np.dot(self.MOMIC.wdot_soot,self.gas.molecular_weights)*self.vol			
			
			# If energy is coupled
			if self.solve_coupled:
				####### Mass #########
				# Soot Consumption Rate			
				
			
				S_mass+= self.mdot_eat
				
				###### Species #######
				S_species += self.MOMIC.wdot_soot*self.gas.molecular_weights*self.vol \
				- self.gas.Y*self.mdot_eat
				
				###### Energy #######
				
				# Need to partially re-work
				Q_h = 1000
				d_s = self.MOMIC.soot_properties(self.M)[0]
				Ap = np.pi*d_s**2
				Q_source = Q_h*Ap*(self.gas.T-self.Tp)
			
				#S_energy+= self.M[0]*self.vol*Q_source*0

		if self.solve_heat_transfer:
		# Next step to implement
			keyboard()
			dz = self.vol/self.Ac
			#q_heat = -(np.pi*self.Dc*dz)*self.convection_h()*(self.gas.T-self.Twall)
			#q_heat = -(np.pi*self.Dc*dz)*700*(self.gas.T-300)

			#S_energy += q_heat


		self.S_mass = S_mass
		self.S_species = S_species
		self.S_energy = S_energy
		
		#return S_mass, S_species, S_energy
		
	def sources_constant(self):
		# Calculate Source Terms outside of time integration to help address stiffness of system
		
		# Now using this to handle droplet evaporation
		
		# Initialize Sources
		S_mass = lambda : 0
		S_species = lambda :0
		S_energy = lambda :0

		# Initial Some Stuff
		# Droplet
		S_mass_droplet = lambda: 0
		S_species_droplet = lambda: 0
		S_energy_droplet = lambda: 0 
		
		# Heat xfer
		S_mass_heat = lambda: 0
		S_species_heat = lambda: 0
		S_energy_heat = lambda: 0


		if self.solve_soot:
			pass
			#self.MOMIC.update_regime_flags()
			
		
		# Droplet Source Term
		if self.solve_droplet:
			m_gas,Df,vap_per = self.droplet.ct_gas(self.gas,self.tres,U=self.U)

			# Use Constant Source for each time step
			# Droplet is calculated in scipy_solver
			#m_gas = self.droplet.m_gas*self.droplet.N_drops

			#m_gas= m_droplet*self.droplet.N_drops 
			#print('Using this for now in sources')
			#m_gas = 0.03
			
			####### Mass #########
			#S_mass = lambda : m_gas
			S_mass_droplet = lambda : m_gas
			
			####### Species #########
			#S_species= lambda : m_gas*(self.inlet_liquid_gas.Y-self.gas.Y)
			S_species_droplet= lambda : m_gas*(self.inlet_liq.Y-self.gas.Y)
			
			####### Energy #########
		
			# Vaporization Enthalpy
			S_energy1 =  lambda : -m_gas*self.droplet.LHV
			
			# Vaporized gas energy 
			S_energy2 = lambda : m_gas*(self.inlet_liq.h-\
				np.dot(self.gas.partial_molar_int_energies/self.gas.molecular_weights,self.inlet_liq.Y))
			
			#S_energy = lambda : S_energy1() + S_energy2()

			S_energy_droplet = lambda : S_energy1() + S_energy2()

			# Balance of liquid/gas mass)
			#self.mdot_gas+= S_mass()
			self.mdot_vap = S_mass_droplet()#S_mass()
			#self.liquid_mdot = self.droplet.mdot-self.vaporized_mdot

		if self.solve_heat_transfer:
			# Next step to implement
			dz = self.vol/self.Ac

			#S_energy3 = lambda : self.Q_source#-(np.pi*self.Dc*dz)*self.convection_h()*(self.gas.T-300)
			#S_energy = lambda : S_energy3()
			S_energy_heat = lambda : self.Q_source#-(np.pi*self.Dc*dz)*self.convection_h()*(self.gas.T-300)
			

	
		#self.S_mass_t = lambda : S_mass()
		#self.S_species_t = lambda : S_species()
		#self.S_energy_t = lambda : S_energy() 
		self.S_mass_t = lambda : S_mass_droplet() 
		self.S_species_t = lambda : S_species_droplet()
		self.S_energy_t = lambda : S_energy_droplet() + S_energy_heat() 
				

	def sources_soot_decoupled(self):
	
		# Initialize Sources
		S_mass = 0
		S_species = 0
		S_energy = 0
		
		# Droplet Source Term
		if self.solve_droplet:
			#[m_droplet,Df] = self.droplet.ct_gas(self.gas,self.tres)
			#keyboard()
			m_gas= self.droplet.m_gas*self.droplet.N_drops 
			#print('Using this for now in sources')
			#m_gas = 0.03
			
			####### Mass #########
			S_mass+= m_gas
			
			####### Species #########
			S_species+= m_gas*(self.inlet_liquid_gas.Y-self.gas.Y)
			
			####### Energy #########
		
			# Vaporization Enthalpy
			S_energy+= -m_gas*self.droplet.LHV
			
			# Vaporized gas energy 
			S_energy+= m_gas*(self.inlet_liquid_gas.h-\
				np.dot(self.gas.partial_molar_int_energies/self.gas.molecular_weights,self.inlet_liquid_gas.Y))
		
			# Balance of liquid/gas mass
			self.mdot_gas+= S_mass
			self.mdot_liq = self.liquid_mdot-m_gas

		if self.solve_soot:
			
			# Run Soot calculation
			self.MOMIC.ct_rate()	
			
			# Soot Consumption Rate			
			self.mdot_eat = np.dot(self.MOMIC.wdot_soot,self.gas.molecular_weights)*self.vol
				
		#if self.solve_heat:
		# Next step to implement
	
		self.S_mass = S_mass
		self.S_species = S_species
		self.S_energy = S_energy
		
		#return S_mass, S_species, S_energy



	
	def integrate_coupled_closed(self,dt,tf,convergence_criteria=None):
			
		Y = self.gas.Y
			
		# Initial Conditions
		y0 = np.hstack([self.gas.T,self.gas.density,Y,self.M,self.P,0,self.gas.T])
		
		# Indices
		self.indx = {'T':0,'D':1,'Y':np.arange(2,self.gas.n_species+2)}
		self.indx.update({'M':np.arange(self.indx['Y'][-1]+1,self.indx['Y'][-1]+1+self.M_moments)})
		self.indx.update({'P':np.arange(self.indx['M'][-1]+1,self.indx['M'][-1]+1+self.P_moments)})
		self.indx.update({'mp':self.indx['M'][-1]+self.P_moments+1})
		self.indx.update({'Tp':self.indx['mp']+1})
		self.indx.update({'t':self.indx['Tp']+1})

		# Source Term ODE
		sol_name = self.source_terms_coupled_closed#_pfr
		
		#sol_name = self.source_terms_coupled_closed
		
			
		## Set Up ODE
		# Defining source term function and numerical Jacobian function. scipy.integrate.ode has an
		# automatic finite-difference jacobian option but it seems to be doing crazy things which 
		# I cant figure out why. 
		
		solver = scipy.integrate.ode(sol_name)#,jac=self.jacobian_coupled)		
		solver.set_integrator('lsoda',method='bdf', with_jacobian=True)
		
		## 
		# Set initial Conditions
		solver.set_initial_value(y0, 0.0)
		# Original input dt
		dt_original = copy.deepcopy(dt)
		
		# Total number of time points to save data
		n_pts = int(np.ceil(tf/dt))*2
		# Preallocate Solution Array
		self.Sol = np.zeros([n_pts,np.size(y0)+1])

		# Initialize time and counters
		i = 0
		count = 0
		self.time = 0

		while solver.t < tf:
		
			# Save Previous Solution
			Y0 = solver.y
			t0 = solver.t
			# Integrate in time	
			solver.integrate(solver.t+dt)

			# Get Solution at current time step
			self.time = solver.t
			self.gas.TDY = solver.y[self.indx['T']],solver.y[self.indx['D']],solver.y[self.indx['Y']]
			self.M = solver.y[self.indx['M']]
			self.P = solver.y[self.indx['P']]
			
			self.state = np.hstack([self.gas.T,self.gas.density,self.gas.Y,\
									self.M,self.P,\
									solver.y[self.indx['mp']],\
									solver.y[self.indx['Tp']],\
									self.time])
			# Update counter
			count+= 1
			
			
			# Check Solution and reduce time step (for n time steps) 
			# if any criteria is met
			if any(np.isnan(solver.y)):
				print('y is nan')
				solver.set_initial_value(Y0, t0)
				dt*=.1
				print('Decreasing time step to:{}'.format(str(dt)))
				count=-1*dt_original/dt
				keyboard()
				
			elif any(solver.y[self.indx['Y']]>1):
				print('yi>1')
				solver.set_initial_value(Y0, t0)
				dt*=.1
				print('Decreasing time step to:{}'.format(str(dt)))
				count=-1*dt_original/dt
				keyboard()
			
			# Return time step to original dt after integrating for 10 decreased
			# time steps
			if count==0:
				dt = dt_original
				print('Reseting time step back to original stepsize')
			
			# Save Solution
			if count>0:
				self.Sol[i,:] = self.state
				i+=1
				
			# Check if Convergence Criteria Specified
			if convergence_criteria:
				if all(np.abs(sol_name(solver.t,solver.y))<convergence_criteria):
					print('Convergence Criteria Met')
					self.Sol = self.Sol[~(self.Sol==0).all(axis=1)]
					print('Solution Trimmed')
					return
			
		# Return nonzero Solution Array
		self.Sol = self.Sol[~(self.Sol==0).all(axis=1)]
	
	def integrate_coupled_P(self,dt,tf,convergence_criteria=None):
			
		Y = self.gas.Y
			
		# Initial Conditions
		y0 = np.hstack([self.gas.T,Y,np.log(self.M),np.log(self.P),0,self.gas.T])
		
		# Indices
		self.indx = {'T':0,'Y':np.arange(1,self.gas.n_species+1)}
		self.indx.update({'M':np.arange(self.indx['Y'][-1]+1,self.indx['Y'][-1]+1+self.M_moments)})
		self.indx.update({'P':np.arange(self.indx['M'][-1]+1,self.indx['M'][-1]+1+self.P_moments)})
		self.indx.update({'mp':self.indx['M'][-1]+self.P_moments+1})
		self.indx.update({'Tp':self.indx['mp']+1})
		self.indx.update({'t':self.indx['Tp']+1})

		# Source Term ODE
		sol_name = self.source_terms_coupled_const_P#_pfr
		
			
		## Set Up ODE
		# Defining source term function and numerical Jacobian function. scipy.integrate.ode has an
		# automatic finite-difference jacobian option but it seems to be doing crazy things which 
		# I cant figure out why. 
		
		solver = scipy.integrate.ode(sol_name)#,jac=self.jacobian_coupled)		
		solver.set_integrator('lsoda',method='bdf', with_jacobian=True)
		
		## 
		# Set initial Conditions
		solver.set_initial_value(y0, 0.0)
		# Original input dt
		dt_original = copy.deepcopy(dt)
		
		# Total number of time points to save data
		n_pts = int(np.ceil(tf/dt))*2
		# Preallocate Solution Array
		self.Sol = np.zeros([n_pts,np.size(y0)+1])

		# Initialize time and counters
		i = 0
		count = 0
		self.time = 0

		while solver.t < tf:
		
			# Save Previous Solution
			Y0 = solver.y
			t0 = solver.t
			# Integrate in time	
			solver.integrate(solver.t+dt)

			# Get Solution at current time step
			self.time = solver.t
			self.gas.TPY = solver.y[self.indx['T']],self.gas.P,solver.y[self.indx['Y']]
			self.M = np.exp(solver.y[self.indx['M']])
			self.P = np.exp(solver.y[self.indx['P']])
			
			self.state = np.hstack([self.gas.T,self.gas.Y,\
									np.log(self.M),np.log(self.P),\
									solver.y[self.indx['mp']],\
									solver.y[self.indx['Tp']],\
									self.time])
			# Update counter
			count+= 1
			
			
			# Check Solution and reduce time step (for n time steps) 
			# if any criteria is met
			if any(np.isnan(solver.y)):
				print('y is nan')
				solver.set_initial_value(Y0, t0)
				dt*=.1
				print('Decreasing time step to:{}'.format(str(dt)))
				count=-1*dt_original/dt
				keyboard()
				
			elif any(solver.y[self.indx['Y']]>1):
				print('yi>1')
				solver.set_initial_value(Y0, t0)
				dt*=.1
				print('Decreasing time step to:{}'.format(str(dt)))
				count=-1*dt_original/dt
				keyboard()
			
			# Return time step to original dt after integrating for 10 decreased
			# time steps
			if count==0:
				dt = dt_original
				print('Reseting time step back to original stepsize')
			
			# Save Solution
			if count>0:
				self.Sol[i,:] = self.state
				i+=1
				
			# Check if Convergence Criteria Specified
			if convergence_criteria:
				if all(np.abs(sol_name(solver.t,solver.y))<convergence_criteria):
					print('Convergence Criteria Met')
					self.Sol = self.Sol[~(self.Sol==0).all(axis=1)]
					print('Solution Trimmed')

					return
			
		# Return nonzero Solution Array
		self.Sol = self.Sol[~(self.Sol==0).all(axis=1)]
	
	
	def integrate_coupled2(self,dt,tf,convergence_criteria=None):
			
		Y = self.gas.Y
			
		# Initial Conditions
		y0 = np.hstack([self.gas.T,self.gas.density,Y,np.log(self.M),np.log(self.P),0,self.gas.T])
		# Indices
		self.indx = {'T':0,'D':1,'Y':np.arange(2,self.gas.n_species+2)}
		self.indx.update({'M':np.arange(self.indx['Y'][-1]+1,self.indx['Y'][-1]+1+self.M_moments)})
		self.indx.update({'P':np.arange(self.indx['M'][-1]+1,self.indx['M'][-1]+1+self.P_moments)})
		self.indx.update({'mp':self.indx['M'][-1]+self.P_moments+1})
		self.indx.update({'Tp':self.indx['mp']+1})
		self.indx.update({'t':self.indx['Tp']+1})

		
		# Source Term ODE
		sol_name = self.source_terms_coupled_pfr
		
			
		## Set Up ODE
		# Defining source term function and numerical Jacobian function. scipy.integrate.ode has an
		# automatic finite-difference jacobian option but it seems to be doing crazy things which 
		# I cant figure out why. 
		
		solver = scipy.integrate.ode(sol_name)#,jac=self.jacobian_coupled)		
		solver.set_integrator('lsoda',method='bdf', with_jacobian=True)
		
		## 
		# Set initial Conditions
		solver.set_initial_value(y0, 0.0)
		# Original input dt
		dt_original = copy.deepcopy(dt)
		
		# Total number of time points to save data
		n_pts = int(np.ceil(tf/dt))*2
		# Preallocate Solution Array
		self.Sol = np.zeros([n_pts,np.size(y0)+1])

		# Initialize time and counters
		i = 0
		count = 0
		self.time = 0
		
		self.mdot_out = self.mdot

		Ac = np.pi*(1.17)**2/4*(.0254)**2
		dz = self.vol/Ac
		self.Ac = Ac
		u1 = self.mdot/(self.inlet_gas.density*Ac)
		
		P1 = 250*101325/14.7
		
		mdot_in = self.mdot
		
		self.P_outlet = self.gas.P
		solver.integrate(solver.t+dt)
		self.source_terms_coupled_pfr(solver.t,solver.y)
		u1 = self.mdot_out/(self.gas.density*Ac)
		P1 = self.gas.P
		
		P1 = 230*101325/14.7

		rho1 = self.gas.density
		
		while solver.t < tf:
		
			# Save Previous Solution
			Y0 = solver.y
			t0 = solver.t
			
			# Integrate in time			
			solver.integrate(solver.t+dt)

		
			#if solver.y[1]<Y0[1]:
			#	keyboard()
			

			# Get Solution at current time step
			self.time = solver.t
			self.gas.TDY = solver.y[self.indx['T']],solver.y[self.indx['D']],solver.y[self.indx['Y']]
			self.M = np.exp(solver.y[self.indx['M']])
			self.P = np.exp(solver.y[self.indx['P']])
			
			self.state = np.hstack([self.gas.T,self.gas.density,self.gas.Y,\
									np.log(self.M),np.log(self.P),\
									solver.y[self.indx['mp']],\
									solver.y[self.indx['Tp']],\
									self.time])
			self.source_terms_coupled_pfr(solver.t,solver.y)
			
			u2 = (self.mdot+self.mdot_eat)/(self.gas.density*Ac)
			
			P2 = P1 -self.gas.density*u2*(u2-u1)/dz-u2*(self.mdot_eat)/Ac#-101325/14.7
			#keyboard()
			
			R = 8314.4621
			MW = self.gas.mean_molecular_weight
			
			self.mdot_out = u2*(self.gas.density)*Ac +(self.gas.P-P2)/self.gas.P#*1e-2
			#keyboard()
			
			#self.mdot_out = #P2/self.gas.P*self.mdot_out
			#self.mdot_out = 0.06+ (self.gas.P-P2)#P2/(R/MW*self.gas.T)*Ac*u2
			
			#self.gas.P/P2*self.mdot_out#np.max([0,(self.gas.P-240*101325/14.7)])*3e-5 #P2/self.gas.P*self.mdot_out
			
			self.P_outlet = P2
			#solver.set_initial_value(solver.y, solver.t)
			#print(self.mdot_out)
			#yboard()

			if count>100:
				print(self.gas.P*14.7/101325,P2*14.7/101325)
				self.Ac = Ac
				self.U = u2
				self.check_mdot()
				keyboard()
			
			#if solver.t>1e-3:
			#	keyboard()
			
			# Update counter
			count +=1
				
			# Check Solution and reduce time step (for n time steps) 
			# if any criteria is met
			if any(np.isnan(solver.y)):
				print('y is nan')
				solver.set_initial_value(Y0, t0)
				dt*=.1
				print('Decreasing time step to:{}'.format(str(dt)))
				count=-1*dt_original/dt
				keyboard()
				
			elif any(solver.y[self.indx['Y']]>1):
				print('yi>1')
				solver.set_initial_value(Y0, t0)
				dt*=.1
				print('Decreasing time step to:{}'.format(str(dt)))
				count=-1*dt_original/dt
				keyboard()
			
			# Return time step to original dt after integrating for 10 decreased
			# time steps
			if count==0:
				dt = dt_original
				print('Reseting time step back to original stepsize')
			
			# Save Solution
			if count>0:
				self.Sol[i,:] = self.state
				i+=1
				
			# Check if Convergence Criteria Specified
			if convergence_criteria:
				if all(np.abs(sol_name(solver.t,solver.y))<convergence_criteria):
					print('Convergence Criteria Met')
					self.Sol = self.Sol[~(self.Sol==0).all(axis=1)]
					print('Solution Trimmed')

					return
			
		# Return nonzero Solution Array
		self.Sol = self.Sol[~(self.Sol==0).all(axis=1)]
	
	
	def calc_source(self):
		Out = momic.momic.calculate_source(self.M,self.P,	\
			self.gas.T, self.gas.P,						\
			self.concentration('C2H2'), 				\
			self.concentration('H'), 					\
			self.concentration('H2'), 					\
			self.concentration('H2O'), 					\
			self.concentration('O2'), 					\
			self.concentration('OH') )					\

		# Split MOMIC Rates Outputs
		W = Out[0]
		G = Out[1]
		R = Out[2]
		Hr = Out[3]
		Ragg = Out[4]
		
		SM = W+G+R
		SP = Hr + Ragg

		
		rC2H2 = Out[5]*1e3
		rCO = Out[6]*1e3
		rH = Out[7]*1e3
		rH2 = Out[8]*1e3
		rH2O = Out[9]*1e3
		rO2 = Out[10]*1e3
		rOH  = Out[11]*1e3
		
	
		gas = self.gas
		#tres = self.tres#rho*self.vol/self.mdot
		######
		#THIS TRES IS WRONG HERE
		
		wdot = gas.net_production_rates
		
		wdot2 = wdot[:]
		wdot2[:] = 0.0
		
		wdot2[gas.species_index('C2H2')] += rC2H2
		wdot2[gas.species_index('CO')] += rCO
		wdot2[gas.species_index('H')] += rH
		wdot2[gas.species_index('H2')] += rH2
		wdot2[gas.species_index('H2O')] += rH2O
		wdot2[gas.species_index('O2')] += rO2
		wdot2[gas.species_index('OH')] += rOH
		#wdot2[:] = 0.0
		
		self.source = wdot2
		self.momic_S = wdot2

	
	def integrate_coupled_source(self,dt,tf):
		
		Y = self.gas.Y
		Y[np.where(Y<=0)] = 1e-300
		
		# Initial Conditions
		y0 = np.hstack([np.log(self.M),np.log(self.P),self.gas.density,self.gas.T,Y])

		# ODE
		solver = scipy.integrate.ode(self.source_terms_const_source)
		
		# Initial Conditions
		#y0 = np.append(np.log(self.M),np.log(self.P))

		# ODE
		#solver = scipy.integrate.ode(self.scipy_ode_uncoupled)
		
		solver.set_integrator('vode', method='bdf', with_jacobian=True)
		solver.set_initial_value(y0, 0.0)
		
		j = 0
		sol_store = y0
		t_store = 0
		dt_reset_indx = 100
		dt_original = copy.deepcopy(dt)
		
		count = 0
		i = 1
		#keyboard()
		while solver.t < tf:
		
			if dt_reset_indx == 100:
				dt = dt_original
				print('reseting time step')
				
			dt_reset_indx+=1
			
			Y0 = solver.y
			solver.integrate(solver.t+dt)
			count +=1
			
		
			#if self.count>19000:
			#	self.couple = False
			
			#if self.count>40000:
			#	keyboard()

			if any(np.isnan(solver.y)):
				keyboard()
			
			if any(solver.y<0):
				solver.set_initial_value(Y0, solver.t)
				dt*=.1
				count=-10
				#keyboard()
			
			#if solver.t > .1*i:
			#	keyboard()
				
			if count>0:
				dt = dt_original
			
			if np.abs(self.dTdt2)<1e-4:
				if np.abs(self.drhodt)<1e-2:
					if all(np.abs(self.dYdt)<1e-3):
						print('Convergence Met')
						return
			
			#keyboard()
			#if self.decrease_dt:
			#	solver.set_initial_value(sol_store,t_store)
			#	dt = dt*1e-1
			#	dt_reset_indx = 0
			#	print('Decreasing time step to:{}'.format(str(dt)))
			#	self.decrease_dt = False



		
		
			sol_store =  solver.y
			t_store = solver.t	
			
		
	
		
	def integrate_to_SS(self,max_iter):
	
		i = 0
		dt = 1e-4
		while i < max_iter:
			i+=1
			j = 0
			#if i == 15:
			#	keyboard()
			# Integrate for 100 time steps
			self.integrate(dt,dt*100)
			
			#while j < 100:
			#	self.step_in_time(1e-6)
			#	j+= 1
				
			#self.integrate(1e-6,1e-6*100)
			
			#keyboard()
			
			# Take Newton Steps
			#self.Newton_Step() 
			
			if self.aggregation:
				if all(np.abs(self.dlog_rates()[0])<1e-4) and \
				all(np.abs(self.dlog_rates()[1][1:])<1e-4):
					print('Convergence Criteria Met')
					return 1
			else:
				if all(np.abs(self.dlog_rates()[0])<1e-4):
					print('Convergence Criteria Met')
					return 1

		# If SS reached return 1, else return 0
		return 0
	
	def integrate_to_SS_coupled(self,max_iter):
		# Integrate Coupled Gas + Momic to steady state
		
		# Doing a Lazy Way Approach => Calculate SS Uncoupled => Use Solution to integrate Steady
		self.integrate_to_SS(max_iter)
		
		# Initial Conditions
		y0 = np.hstack([np.log(self.M),np.log(self.P),self.gas.T,self.gas.Y])

		# ODE
		solver = scipy.integrate.ode(self.source_terms_coupled)
		solver.set_integrator('vode', method='bdf', with_jacobian=True)
		solver.set_initial_value(y0, 0.0)
		dt = 1e-4
		
		
		iM = 0
		iP = np.size(self.M)
		iT = np.size(self.M)+np.size(self.P)
		iY = iT +1
		
		#keyboard()
		# Integrate Coupled for some time
		while solver.t < 0.5:
			solver.integrate(solver.t+dt)
			print(solver.t)
		
		#t_cond = 0.5
		# Integrate Coupled for some time and check if converged 
		while solver.t < 4:
			self.M = np.exp(solver.y[iM:iP])	
			self.P = np.exp(solver.y[iP:iT])	
			self.gas.TPY = solver.y[iT],self.gas.P,solver.y[iY:]
			solver.integrate(solver.t+dt)
			
			#if solver.t > t_cond:
			#	keyboard()
			
			# Check if Converged
			#print(solver.t)
			if all(np.abs(self.dlog_rates()[0])<1e-4):
				print('Convergence Criteria Met')
				return
			
		# If not fully converged
		print('Not Fully Converged')
		keyboard()

	
	def jacobian_mkl(self):
		
		J = momic.momic.jacobian_log(np.append(self.M,self.P),\
			np.size(self.M), np.size(self.P), 									\
			self.gas.T, self.gas.P, 		\
			self.concentration('C2H2'), 			\
			self.concentration('H'), 				\
			self.concentration('H2'),				\
			self.concentration('H2O'), 				\
			self.concentration('O2'), 				\
			self.concentration('OH') )
		
		return J
		
	
	def jacobian_coupled(self,t,x):
		eps = 1e-5 
		fun = lambda x_eps: self.source_terms_coupled_Jac(t,x_eps)
		
		J = scipy.optimize.slsqp.approx_jacobian(x,fun,eps)
		
		#J+=-np.mean(J)
		#J*= 1/np.max(J)
		
		#J = scipy.sparse.csc_matrix(J)
		
		keyboard()
		
		return J#scipy.optimize.slsqp.approx_jacobian(x,fun,eps)
	
	def sol_P(self):
		

		P = np.zeros(np.shape(self.Sol)[0])
		
		gas = ct.Solution(self.reaction_mechanism)

		for i,_ in enumerate(P):
			gas.TDY = self.Sol[i,self.indx['T']], self.Sol[i,self.indx['D']], self.Sol[i,self.indx['Y']]
			P[i] = gas.P
		
		return P

	def sol_concentrations(self):
	
		C = np.zeros([np.shape(self.Sol)[0],self.gas.n_species])
		
		gas = ct.Solution(self.reaction_mechanism)
		
		for i,_ in enumerate(self.Sol[:,0]):
			gas.TDY = self.Sol[i,self.indx['T']], self.Sol[i,self.indx['D']], self.Sol[i,self.indx['Y']]
			
			C[i,:] = gas.concentrations
			
		return C
	
	def sol_soot_properties(self):
		
		# Moments
		M = np.exp(self.Sol[:,self.indx['M']])
		
		# Calculate soot properties
		prop = np.matrix([self.soot_properties(Mi) for Mi in M])
		
		return prop 

	def def_reactor_properties(self,Dc,Twall):
		# Define Reactor Diameter
		# Need to call this function for transport calc

		self.Dc = Dc
		self.Ac = np.pi*Dc**2/4
		self.Twall = 300
		self.solve_heat_transfer = True

	def gas_velocity(self):
		# Calculate Gas Velocity

		return self.mdot/(self.gas.density*self.Ac)


	def Re(self):
		# Calculate Reynolds Number for WSR composition

		return self.gas.density*self.gas_velocity()*self.Dc/ \
		self.gas.viscosity

	def Nu(self):
		# Calculate Turbulent Nusselt Number for WSR composition

		# Reynolds Number 
		Re = self.Re()
		# Prandtl Number
		Pr = self.gas.cp*self.gas.viscosity/self.gas.thermal_conductivity

		A = ((8/Re)**(10)+(Re/36500)**(20))**(0.5)
			
		# From churchill
		A = ((8/Re)**(10)+(Re/36500)**(20))**(0.5)
		B = (2.2*np.log(Re/7))**(10)
		C = (1/A+B)**(1/5)

		f = 8/C

		Nu_t =  5.76+(0.079*Re*np.sqrt(f/8)*Pr)/(1+Pr**(4/5))**(5/6)

		return Nu_t
	
	def convection_h(self):
		return self.Nu()*self.gas.thermal_conductivity/self.Dc

	def brownian_diff(self):
		# Calculate Brownian Diffusion coefficient based on 
		# average soot properties

		kb = 1.38065e-23	# Boltzmann constant J/K

		# Soot Particle Diameter
		dp = self.soot_properties(self.M)[0]

		# Knudsen #
		Kn = self.kn(self.gas,dp)

		# Cunningham correction
		Cc = 1+Kn*(1.2+0.41*exp(-0.88/Kn))

		return kb*Cc*T/(3*np.pi*self.gas.viscosity*dp)

	def Sh(self,DB=None):
		# Calculate Sherwood Number for WSR composition
		return 0.0165*self.Re()**(0.86)*self.Sc(DB)**(0.333)
	
	def Sc(self,DB=None):
		# Calculate Schmidt Number for WSR composition
		
		if DB is None:
			return self.gas.viscosity/(self.gas.density*self.brownian_diff())
		else:
			return self.gas.viscosity/(self.gas.density*DB)



class droplet:
	# Class for calculating droplet properties
	# This was originally part of WSR class but moved to seperate class
	# to enable for easier future development
	
	def __init__(self,liq_rho,liq_mu,LHV,T_boil,reaction_mechanism = None):
	
		self.LHV = LHV
		self.T_boil = T_boil
		self.liq_rho = liq_rho
		self.liq_mu = liq_mu
	
		# Create Gas Object for internal calculations if reaction
		# mechanism specified
		# Need this if using ct_gas method, otherwise need to manually
		# specific gas phase properties
		if reaction_mechanism is not None:
			self.gas = ct.Solution(reaction_mechanism)
			
	def initialize(self,D0,m_liquid,convection=False):
		# Initial droplet size and total mass flow
		self.D0 = D0
		self.m_liq = m_liquid
		self.N_drops = self.m_liq/(np.pi*D0**3/6*self.liq_rho)
		self.m_drop = self.m_liq/self.N_drops
		self.m_gas = 0
		self.convection = convection
		
	def gas_phase_input(self):
		# Manually define properties to calculate rates
		raise ValueError('Not Implemented')
				
	def ct_gas(self,gas_object,tres,U=0):
		# Calculate Droplet vaporization using
		# Cantera for gas property evaluation

		# If droplet size is 0, then it's vaporized, so return 0
		if self.D0 == 0:
			self.Df = 0
			self.m_gas = 0
			return 0,0,1
		
		# If gas T < boil T then no vaporization can occur
		if gas_object.T < self.T_boil:
			self.Df = self.D0
			self.m_gas = 0
			return self.m_gas,self.Df,0

		# Mean Temperature
		Tmean = (gas_object.T+self.T_boil)/2

		# Calculate cp of gas-phase fuel at mean temperature
		self.gas.TPY = Tmean,gas_object.P,self.gas.Y
		
		# Store Values to pass to other calculations
		self.cp_inf = gas_object.cp
		self.cp_mean = self.gas.cp
		self.k_inf = gas_object.thermal_conductivity
		self.k_mean = self.gas.thermal_conductivity
		self.T_inf = gas_object.T
		self.mu_inf = gas_object.viscosity
		self.Pr = self.cp_inf*self.mu_inf/self.k_inf
		self.B = np.log(1+self.cp_mean*(self.T_inf-self.T_boil)/(self.LHV))
		self.gas_U = U
		self.rho_g = gas_object.density
		self.mu_g = gas_object.viscosity
		
		# Droplet with convection or without
		if U==0 or (self.convection != True):
			Df = self.analytic_no_convection(self.D0,tres)
		else:
			# Re# term makes the integration messy
			# so going with a lazy approach with numerical integration
			if not hasattr(self,'integrator'):
				self.integrator = scipy.integrate.ode(self.scipy_rates)
				self.integrator.set_integrator('dopri5')
			
			# Some weird crap is going on with len() of unsized object
			# Should look into this after my defense
			try:
				self.integrator.set_initial_value(self.D0,0)
			except:
				try:
					self.integrator = scipy.integrate.ode(self.scipy_rates)
					self.integrator.set_integrator('dopri5')
					self.integrator.set_initial_value(self.D0,0)
				except:
					pass
				pass

			# Integrate ODE
			self.integrator.integrate(tres)	

			# Df can't be below 0
			Df = np.max(np.append(0,self.integrator.y))


		m_gas = self.liq_rho*np.pi/6*(self.D0**3-Df**3)*self.N_drops

		# Store droplet and gas mass
		self.Df = Df
		self.m_gas = m_gas
		self.vap_per = self.m_gas/self.m_liq
		
		return m_gas, Df, self.vap_per
	
	def analytic_no_convection(self,D0,tres):
		# Analytical function for droplet vaporization
		
		Df2 = D0**2 
		Df2 += -8/self.liq_rho*(self.k_mean/self.cp_mean)*\
			self.B*tres

		if Df2 <=0:
			return 0
		else:
			# Rounding error can make Df2>D0
			return np.min([D0,np.sqrt(Df2)])
	
	def convect_effect(self,D):
		# Correlation from:
		# Ranz, W. E., and Marshall, W. R., Evaporation from
		# drops, Chem. Eng. Prog., Vol. 48, 1952, Part I, pp. 141–146;
		# Part II, pp. 173–180.
		
		Pr = self.cp_inf*self.mu_inf/self.k_inf
		
		#Re = self.liq_rho*self.gas_U**2*D/self.liq_mu
		Re = self.rho_g*self.gas_U**2*D/self.mu_g
		return (1+0.3*Re**0.5*Pr**0.33)

	def droplet_rates(self,D):	
		
		# If droplet is fully vaporized
		if D<=0:
			return [0,0]
	
		dmfdt = -2*np.pi*D*(self.k_mean/self.cp_mean)*	\
					self.B
	
		# Calculate convective term
		dmfdt*= self.convect_effect(D)
		
		dDdt = dmfdt/(np.pi*self.liq_rho/2*D**2)

		return [dDdt, dmfdt]
	
	def scipy_rates(self,t,y):
		return self.droplet_rates(y)[0]



	
class PFR:
	# PFR class that uses Moments Class

	def __init__(self,reaction_mechanism,surface_mechanism=None,constant_pressure=False,WSR_object = None,main_class=None):
		
		if main_class is not None:
			self.main = main_class
			self.reaction_mechanism = self.main.reaction_mechanism
			self.surface_mechanism = self.main.surface_mechanism
			self.deposit = self.main.deposit
			self.WSR = self.main.WSR

		else:
			self.reaction_mechanism = reaction_mechanism
			self.surface_mechanism = surface_mechanism
			self.deposit = deposition(reaction_mechanism,surface_mechanism)
			self.WSR = WSR(reaction_mechanism)
		
		self.inlet_gas = None
		self.inlet_gas_mdot = None
		self.inlet_gas_z = None
		self.inlet_liq = None
		self.inlet_liq_mdot = None
		self.inlet_liq_D0 = None
		self.inlet_liq_z = None

		self.inlet_M = None
		self.inlet_P = None
		self.vol = None
		self.Ac = None
		self.z = None
		self.Sol = []
		
		self.count = 0	# Number of Times the Function is called 
		
		self.outlet_BC = []
		self.outlet_gas = []
		self.outlet_At = []

		# Scipy solver inputs 
		self.wsr_dt = 1e-3
		self.wsr_tf = 2
		self.wsr_convergence_critera = 1e-3
		self.Pc_guess = None

		# Control Flags 
		self.solve_droplet_pfr = False
		self.solve_soot = False

	def pressure_guess(self,P_guess):
		self.Pc_guess = P_guess

	def error_check(self):
		# Check for Errors

		if self.z is None:
			raise ValueError('No z_pts specified, specify with set_reactor_pts method')
		if np.size(self.z) == 1:
			raise ValueError('Only 1 WSR specified for PFR model, use WSR model instead')
	
	def save_solution(self,i,new_Sol=False,**data):
		# Store WSR properties in PFR using
		# Cantera SolutionArray
		
		Dict = {'z':self.z[i],'U':self.WSR.U}
		Dict.update(data)
		if self.solve_soot:
			Dict.update(self.WSR.MOMIC.moment_dict())
			fv, rho_mix = self.soot_gas_properties()
			additional_data = {'fv':fv}
			Dict.update(additional_data)
		if self.solve_droplet_pfr:
			Dict.update({'Drop':self.WSR.droplet.Df})

		if self.Sol is None or new_Sol:
			del(self.Sol)
			self.Sol = ct.SolutionArray(self.WSR.gas,extra=Dict)
			self.Sol.append(self.WSR.gas.state,**Dict)
		else:
			self.Sol.append(self.WSR.gas.state,**Dict)

	def update_control_flags(self):
		# Update Control Flags in PFR/WSR
		if self.main is not None:
			self.solve_coupled = self.main.solve_coupled
			self.solve_soot = self.main.solve_soot
		else:
			self.WSR.solve_coupled = self.coupled
			self.WSR.solve_soot = self.solve_soot

	def set_inlet_gas(self,T,P,Y,mdot,z=0):
		
		if self.inlet_gas is None:
			self.inlet_gas =[[ct.Solution(self.reaction_mechanism)]]
			self.inlet_gas[0][0].TPY = T,P,Y
			self.inlet_gas_mdot = [[mdot]]
			self.inlet_gas_z = [[z]]
		elif z in self.inlet_gas_z:
			z_loc = self.inlet_gas_z.index(z)
			self.inlet_gas[z_loc].append(ct.Solution(self.reaction_mechanism))
			self.inlet_gas[z_loc][-1].TPY = T,P,Y
			self.inlet_gas_mdot[z_loc].append(mdot)
		else:
			self.inlet_gas.append([ct.Solution(self.reaction_mechanism)])
			self.inlet_gas[-1][-1].TPY = T,P,Y
			self.inlet_gas_mdot.append([mdot])
			self.inlet_gas_z.append(z)
			
	def set_inlet_liq(self,T=None,gas=None,D0=None,mdot=None,rho=None,mu=None,LHV=None,Tboil=None,convection=False,z=0):
		
		# Error Check Input
		if T is None:
			raise ValueError('Liquid Injection Initial Temperature "T0" Not Entered')
		elif gas is None:
			raise ValueError('Gas Species of Corresponding Liquid Injection "gas" Not Entered')
		elif D0 is None:
			raise ValueError('Initial Droplet Size of Liquid Injection "D0" Not Entered')
		elif mdot is None:
			raise ValueError('Liquid Injection Mass Flow Rate "mdot" Not Entered')
		elif rho is None:
			raise ValueError('Liquid Injection Density "rho" Not Entered')
		elif mu is None:
			raise ValueError('Liquid Injection Viscosity "mu" Not Entered')
		elif LHV is None:
			raise ValueError('Liquid Injection Latent Heat of Vaporization "LHV" Not Entered')
		elif Tboil is None:
			raise ValueError('Liquid Injection Boiling Temperature "Tboil" Not Entered')

		if self.inlet_liq_z is None:
			self.inlet_liq_z = [z]
			self.inlet_liq = [ct.Solution(self.reaction_mechanism)]
			self.inlet_liq[0].TPY = T, 101325, gas # Pressure is Not used for anything
			self.inlet_liq_properties = [T,gas,D0,mdot,rho,mu,LHV,Tboil,convection]
			self.solve_droplet_pfr = True
		elif z in self.inlet_liq_z:
			# 1+ Droplet Injections not implemented
			raise NotImplementedError('Multiple Liquid Injection Points Not Yet Implemented')
		else:
			raise NotImplementedError('Multiple Liquid Injection Points Not Yet Implemented')
		
	def set_inlet_Moments(self,M,P=[]):
		self.inlet_M = np.array(M)
		self.inlet_P = np.array(P)	
	
	def set_PFR_outlet(self,outlet_BC_type,P=101325,At=1):
		
		self.outlet_BC = outlet_BC_type
		self.outlet_gas = ct.Solution(self.reaction_mechanism)
		self.outlet_gas.TPY = 300,P,'N2:1'
		self.outlet_At = At

		#if self.outlet_BC == 'PC':
		#	self.main.constant_pressure = True


	def set_reactor_pts(self,z,Ac):
		
		# Add Reactor Z and Ac to object class
		self.z = np.array(z)
				
		if np.size(Ac) == 1:
			self.Ac = Ac*np.ones(np.size(self.z))
		
		if np.size(Ac) >1:
			if np.size(Ac) != np.size(Z):
				raise ValueError('Ac Length does not match Z length')
			else:
				self.Ac = Ac
	
		# PFR Volume
		self.vol = np.sum(np.diff(self.z)*Ac)+self.z[0]*self.Ac[0]
	
	def set_reactor_properties(self,Twall=None):

		if Twall is not None:
			self.WSR.Twall = Twall

	def sort_inlets(self):
		# Organize inlets

		# Create inlet attribute for each reactor
		self.reactor_inlet_gas = [None]*np.size(self.z)
		self.reactor_inlet_liq = [None]*np.size(self.z)

		# Add inputs to reactor_inlet_gas list
		for i,z in enumerate(self.inlet_gas_z):			
			# location of reactor where to add reactor gas inlet
			loc = np.where(self.z-z>=0)[0][0]

			if self.reactor_inlet_gas[loc]:
				self.reactor_inlet_gas[loc].append([[self.inlet_gas[i][j],self.inlet_gas_mdot[i][j]] for j,_ in enumerate(self.inlet_gas[i])][0])
			else:
				self.reactor_inlet_gas[loc] = [[self.inlet_gas[i][j],self.inlet_gas_mdot[i][j]] for j,_ in enumerate(self.inlet_gas[i])]
		
		# Combine inlets in reactor_inlet_gas into a single averaged inlet
		for i,inlets in enumerate(self.reactor_inlet_gas):
			if np.size(inlets)>1:
				# Gas Mixture Object
				temp_gas = []
				for inlet in inlets:
					temp_gas.append(ct.Solution(self.reaction_mechanism))
					temp_gas[-1].TPY = inlet[0].TPY
				
				gas = np.sum([ct.Quantity(temp_gas[j],mass = inlets[j][1],constant='HP') for j,_ in enumerate(inlets)])				
				# Total mass 
				mass = np.sum([inlets[j][1] for j,_ in enumerate(inlets)])
								
				gas_object = ct.Solution(self.reaction_mechanism)
				gas_object.TPY = gas.TPY
				
				self.reactor_inlet_gas[i] = [gas,mass]

		if self.solve_droplet_pfr:
			for i,z in enumerate(self.inlet_liq_z):
				# location of reactor where to add reactor liquid inlet
				loc = np.where(self.z-z>=0)[0][0]

				if self.reactor_inlet_liq[loc]:
					# Not implemented
					keyboard()
				else:
					self.reactor_inlet_liq[loc] = [self.inlet_liq_properties]
	
	def initialize_pfr_model(self,calculate_WSR0=True):
		# Initialize arrays used in PFR calculations/calculate WSR and setup PFR calculations 

		# Update Control Flags
		self.update_control_flags()
		# Check for errors
		self.error_check()	
		# Create vector of inlets for each reactor
		self.sort_inlets()	
				

	def soot_gas_properties(self):
		
		if self.solve_soot:
			fv = self.WSR.MOMIC.soot_properties(self.WSR.M)[2]
		else:
			fv = 0			

		if self.solve_coupled:
			rho_mix = 1800*fv+self.WSR.gas.density*(1-fv)
		else:
			rho_mix = self.WSR.gas.density
		
		return fv, rho_mix
			
	
	def set_PFR_solver(dt=1e-3,tf=2,convergence_criteria=1e-3):
		
		# Set scipy solver properties
		self.wsr_dt = dt
		self.wsr_tf = tf
		self.wsr_convergence_critera = convergence_criteria
	
	
	def solve(self,calculate_WSRmin=True,iterate_pressure=True):
		
		# Use WSR results for initial Pressure distribution in PFR
		# We're iterating on pressure afterwards so this is just some starting pt
		#P0_dist = self.WSR.gas.P*np.ones(np.size(self.z))
		if self.outlet_BC == 'PC':
			Pc = self.outlet_gas.P
		elif self.Pc_guess is not None:
			Pc = self.Pc_guess
		else:
			# Store Old Z Pts
			z_input = self.z
			Ac_input = self.Ac
			# Make PFR 2 WSR reactor
			self.z = np.array([self.z[-1]/2,self.z[-1]])
			self.Ac = np.array([self.z[-1],self.z[-1]])
			self.initialize_pfr_model()
			self.solve_pfr_model(101325*20,calculate_WSRmin=False,Use_PFR_BC=True)
			Pc = self.Sol.P[-1]
			self.z = z_input
			self.Ac = Ac_input
		
		# Initialize PFR
		# Had some things in here that I didn't need, maybe don't
		# even need this method?
		self.initialize_pfr_model()

		# Solve PFR Model
		Sol = self.solve_pfr_model(Pc,calculate_WSRmin=calculate_WSRmin)

		# If everything is constant pressure, then no need to iterate
		if self.outlet_BC == 'PC':
			#self.Sol = Sol
			return
		
		# Otherwise lets iterate on pressure
		if iterate_pressure:
			# loop counter and flags
			loop_count = 0
			low_bound_found = False
			hi_bound_found = False
			
			# PFR Pressure iteration loop
			while True:
				# Calculate compressible mdot out with given throat At
				mdot_out = self.outlet_At/self.WSR.At*self.WSR.mdot_choked()
				# Get gas mdot from last WSR in PFR chain
				mdot_gas = self.Sol.mdot_gas[-1]

				# Find uppder and lower bounds for bisection search
				if mdot_out < mdot_gas:
					P_low = Pc
					Pc = P_low*1.2
					low_bound_found = True
				else:
					P_hi = Pc
					Pc = P_hi*.8
					hi_bound_found = True
			
				# Start Bisection once upper and lower bounds are found
				if low_bound_found and hi_bound_found:
					# Bisection Loop
					convergence = .5*101325/14.7	# Convergence to within 0.5 psi
					loop_count = 0	# Reset loop count
					while True:
						# Midsection pressure
						Pmid = (P_hi+P_low)/2
						# Update P0
						Pc = Pmid
											
						# Run PFR again
						self.solve_pfr_model(Pc,calculate_WSRmin=False)
						
						loop_count+=1

						# Calculate compressible mdot out with given throat At
						mdot_out = self.outlet_At/self.WSR.At*self.WSR.mdot_choked()
						# Get gas mdot from last WSR in PFR chain
						mdot_gas = self.Sol.mdot_gas[-1]

						if  mdot_out < mdot_gas:
							P_low = Pmid
						else:
							P_hi = Pmid

						# If converged then break out of loop
						if (P_hi-P_low) < convergence:
							#self.Sol = Sol
							return
						# If too much loops then something may be wrong, or initial bounds are too big?
						elif loop_count>20:
							keyboard()

				# Solve PFR Model
				self.solve_pfr_model(Pc,calculate_WSRmin=False)
				# Dumb comment \/:
				# increment loop count by 1					
				loop_count+=1
				
				# If too many loops taken then pause and give terminal access
				if loop_count>10:
					keyboard()
					
	
	def solve_pfr_model(self,Pc,init_M = [],init_P=[],calculate_WSRmin = True,Use_PFR_BC=False):
		# Solve Plug Flow Reactor as a series of 0-D Well Stirred Reactor
		# This is an interface for the Well Stirred Reactor Class 

		# Gas/Liq Mass Flow and soot Balance Counter
		mdot_gas = 0
		mdot_liq = 0
		mdot_soot = 0

		# Gaseous inlet to 1st WSR
		self.WSR.set_inlet_gas(self.reactor_inlet_gas[0][0].T,\
			Pc, \
			self.reactor_inlet_gas[0][0].Y,\
			mdot = self.reactor_inlet_gas[0][1])

		# Account for this gas inlet
		mdot_gas += self.reactor_inlet_gas[0][1]

		# Liquid inlet to 1st WSR
		if self.solve_droplet_pfr:
			self.WSR.set_inlet_liq( 								\
			T = self.reactor_inlet_liq[0][0][0],					\
			gas = self.reactor_inlet_liq[0][0][1],				\
			D0 = self.reactor_inlet_liq[0][0][2],				\
			mdot = self.reactor_inlet_liq[0][0][3],				\
			rho = self.reactor_inlet_liq[0][0][4],				\
			mu = self.reactor_inlet_liq[0][0][5],				\
			LHV = self.reactor_inlet_liq[0][0][6],				\
			Tboil = self.reactor_inlet_liq[0][0][7],				\
			convection = self.reactor_inlet_liq[0][0][8])		
			mdot_liq =  self.reactor_inlet_liq[0][0][3]
		
		# If solving for soot then add inlet Moments
		if self.solve_soot:
			# Define Inlet Moments
			self.WSR.set_inlet_Moments(self.inlet_M,self.inlet_P)
			# Initial Moments are not really important - it just helps
			# speed up convergence
			if init_M:
				M0 = init_M
			else:
				M0 = self.WSR.MOMIC.M
			if init_P:
				P0 = init_P
			else:
				P0 = self.WSR.MOMIC.P

			# Set Initial Moments
			self.WSR.MOMIC.set_moments(M0,P0)
		
		# Setup constant pressrue outlet to 1st WSR
		if Use_PFR_BC:
			self.WSR.set_outlet('PC',P=Pc)
		else:
			self.WSR.set_outlet('PC',P=Pc)

		# Initial WSR Volume, either from estimate of minimum required WSR Volume
		# or initial input
		dV = self.z[0]*self.Ac[0]

		if calculate_WSRmin:
			if self.solve_droplet_pfr:
				V_min = self.WSR.blowout_vol_drop(5e-8,dV)
			else:
				V_min=self.WSR.blowout_vol(5e-8,dV)
			
			# If minimum volume is smaller than minimum volume, quit
			if V_min>self.vol:
				raise ValueError('PFR Volume is too small to sustain combustion')
			
			z0 = V_min/self.Ac[0]
			less_than_Vmin = self.z<z0

			# Delete smaller than Vmin (These are too small and will blowout)
			self.z = self.z[~less_than_Vmin]
			self.Ac = self.Ac[~less_than_Vmin]

			# Now add z0 and Vmin to first z and Ac indx
			self.z = np.append(z0,self.z)
			self.Ac = np.append(self.Ac[0],self.Ac)
			
			# Do inlet sorting magic
			self.sort_inlets()

			# Set dV to V_min
			dV = V_min
			
		# Set Initial WSR Volume
		self.WSR.vol = dV
		# Set WSR Ac
		self.WSR.Ac = self.Ac[0]

		# Integrate 1st WSR
		self.WSR.solve(self.wsr_dt,self.wsr_tf,\
			convergence_criteria = self.wsr_convergence_critera, \
			ct_WSR_IC = True)

		# Tally mass flows
		mdot_gas+=self.WSR.mdot_vap+self.WSR.mdot_eat
		mdot_liq-=self.WSR.mdot_vap
		mdot_soot-=self.WSR.mdot_eat

		# If the initial WSR has blown out, then let's adjust initial WSR volume
		# to something that will work
		if self.WSR.gas.T < 1.5*self.WSR.inlet_gas.T:
			print('Initial WSR volume in PFR chain results in blowout, ' + \
				'increasing volume to minimum volume to prevent blowout')

			if self.solve_droplet_pfr:
				V_min = self.WSR.blowout_vol_drop(5e-8,dV)
			else:
				V_min=self.WSR.blowout_vol(5e-8,dV)

			# If minimum volume is smaller than minimum volume, quit
			if V_min>self.vol:
				raise ValueError('PFR Volume is too small to sustain combustion')

			z0 = V_min/self.Ac[0]
			less_than_Vmin = self.z<z0

			# Delete smaller than Vmin (These are too small and will blowout)
			self.z = self.z[~less_than_Vmin]
			self.Ac = self.Ac[~less_than_Vmin]

			# Now add z0 and Vmin to first z and Ac indx
			self.z = np.append(z0,self.z)
			self.Ac = np.append(self.Ac[0],self.Ac)
			
			# Do inlet sorting magic
			self.sort_inlets()

			# Set dV to V_min
			dV = V_min
			
			# Set Initial WSR Volume
			self.WSR.vol = dV
			# Set WSR Ac
			self.WSR.Ac = self.Ac[0]
			
			# Rerun
			self.WSR.solve(self.wsr_dt,self.wsr_tf,\
				convergence_criteria = self.wsr_convergence_critera, \
				ct_WSR_IC = True)
		
		# Save the 1st WSR solution using Cantera SolutionArray Class
		self.save_solution(0,mdot_gas=mdot_gas,mdot_liq=mdot_liq,mdot_soot=mdot_soot,new_Sol=True)
		
		# Run Reactor Chain
		for i,_ in enumerate(self.z[1:],start=1):

			# Calculate Volume of next WSR reactor (Average of current and next volume)
			dV = (self.z[i]-self.z[i-1])*(self.Ac[i]+self.Ac[i-1])/2
			# Set new volume
			self.WSR.vol = dV
			# Set new cross sectional area
			self.WSR.Ac = self.Ac[i] 

			# Define Outlet BC
			#if i == np.size(self.z)-1:
			#	self.WSR.set_outlet(self.outlet_BC,P=self.outlet_gas.P,At=self.At)
			if Use_PFR_BC:
				self.WSR.set_outlet(self.outlet_BC,P=self.outlet_gas.P,At=self.outlet_At)
			else:
				# Set all except last WSR in chain to constant pressure
				self.WSR.set_outlet('PC', P=Pc)
								
			if self.reactor_inlet_gas[i]:
				# Need to develop!!
				# This is if there are other inlets(x) to PFR (like Fuel Film)
				raise NotImplementedError('Multiple Injection Points not yet implemented')
				#mixed_gas = ct.Quantity(self.WSR.gas,mass=self.WSR.mdot_gas)
				#inlet_gas = ct.Solution(self.reaction_mechanism)
				#inlet_gas.TPY = self.reactor_inlet_gas[i][0].TPY
				#inlet_gas = ct.Quantity(inlet_gas,mass=self.reactor_inlet_gas[i][1])
				#mixed_gas += inlet_gas
				#self.WSR.def_inlet_gas(mixed_gas.T,mixed_gas.P,mixed_gas.Y,mdot=mixed_gas.mass)
				#self.WSR.inlet_Moments(self.WSR.M,self.WSR.P)
			else:
				self.WSR.set_inlet_gas(self.WSR.gas.T,self.WSR.gas.P,self.WSR.gas.Y,mdot=mdot_gas)#self.WSR.mdot_gas)
				self.WSR.set_inlet_Moments(self.WSR.M,self.WSR.P)
			
			# If droplet vaporization is enabled, then update droplet D0 for next calculation
			if self.solve_droplet_pfr:
				self.WSR.droplet.m_liq += -self.WSR.mdot_vap
				self.WSR.droplet.D0 = np.copy(self.WSR.droplet.Df)
			
			# Run next WSR in the chain
			self.WSR.solve(self.wsr_dt*10,self.wsr_tf,\
				convergence_criteria = self.wsr_convergence_critera, \
				ct_WSR_IC = False)#calculate_WSR_IC)
			
			# Tally mass flows
			mdot_gas+=self.WSR.mdot_vap+self.WSR.mdot_eat
			mdot_liq-=self.WSR.mdot_vap
			mdot_soot-=self.WSR.mdot_eat

			# Save Solution
			self.save_solution(i,mdot_gas=mdot_gas,mdot_liq=mdot_liq,mdot_soot=mdot_soot)

class deposition:
	# Class for calculating deposition/soot growth rates
	
	def __init__(self,WSR_class=None,surf_mechanism = False, reaction_mechanism = False, main_class = None):
		
		if surf_mechanism is not None:
			self.surf_mechanism = surf_mechanism
		
		if reaction_mechanism is not None:
			self.reaction_mechanism = reaction_mechanism
		

		if main_class is not None:
			self.main = main_class
			self.gas = self.main.gas
			self.MOMIC = self.main.MOMIC
			self.WSR = self.main.WSR
			self.surface_mechanism = self.main.surface_mechanism



		# Soot Properties
		self.kp = 0.5

	def rates_from_Sol(self,Sol,Twall,PAH=False,skip_HACA=True):
		

		# Calculate Average Soot Particle Diameter
		dp = np.array([self.MOMIC.soot_properties(M)[0] for M in Sol.M])
		dp2 = np.array([self.MOMIC.soot_properties(M)[1] for M in Sol.M])

		V_vec = []
		#Ys = []
		Ysoot2 = []
		mth = []
		tth = []
		hacas = []
		masses = []
		masses2 = []
		mth2 = []
		mdot_conde = []

		m3s = []
		for i, z in enumerate(Sol.z):
			
			self.gas.TPY = Sol.TPY[i,0], Sol.TPY[i,1],Sol.TPY[i,2:]
			self.d_p = dp[i]

			#self.gas.TP = 1731.0349338953688, 3446428.571428564
			#self.d_p = 2.5140856149207635e-09			
			#Ysoot = 5.88e-8


			Vth = self.V_Th(Twall,Sol.mdot_gas[i])

			#Vth = 0.0132
			#keyboard()


			Ysoot = self.MOMIC.soot_mfr(Sol.M[i],self.gas.density)

			msoot = Sol.mdot_soot
			mgas = Sol.mdot_gas

			Ys = msoot/(msoot+mgas)

			#Ysoot2.append(Ysoot)
			blah = self.gas.density*.1*Vth*np.pi*self.WSR.Dc*(3.94*.0254)


			#keyboard()
			dYdx_th = np.pi*self.WSR.Dc*self.gas.density*Ysoot*Vth/Sol.mdot_gas[i]


			if i  == 0:
				dz = np.sum(Sol.z[1:])
			elif i == 1:
				dz = Sol.z[0]
			else:
				dz = Sol.z[i]-Sol.z[i-1]
			
			m_th = dYdx_th*Sol.mdot_gas[i]*dz

			m1 = 2e-23*Sol.M[i,1]*(1e2)**3*np.pi*(self.WSR.Dc)**2/4*dz

			m1 = 2e-23*Sol.M[i,1]*(1e2)**3*1e-3
			
			#m_th = msoot[i]*Vth/Sol.Um[i]#np.pi*self.WSR.Dc*1800*Vth*dz
			#keyboard()
			m_th = Sol.M[i,1]*2.2e-23*1e-3*(1e2)**3*np.pi*self.WSR.Dc*dz*Vth

			mth2.append(m_th)
			
			Y_grad = self.calc_Y_grad(Ysoot,self.d_p)
			Db = self.D_brownian(self.gas,self.d_p)

			mmass = -np.pi*self.WSR.Dc*self.gas.density*Db*Y_grad*dz

			Y_grad2 = self.calc_Y_grad(Sol.M[i,1]*2.2e-23*1e-3*(1e2)**3,dp[i])
			mmass2 = -np.pi*self.WSR.Dc*dz*Db*Y_grad2
			masses.append(mmass)
			masses2.append(mmass2)
		#	keyboard()

			if skip_HACA:
				pass
			else:
				#self.integrate_gas(self.gas,6e-3*10)
				mdot_haca = self.haca_growth(self.gas,Twall)
				
				print(mdot_haca)
				mdot_haca*=np.pi*self.WSR.Dc*dz
				hacas.append(mdot_haca)
			

			dS = np.pi*self.WSR.Dc*dz

			m_th = m1*Vth*dS
			t_th = m_th/(1800*dS)*1e6/.05


			t_blah = t_th/3*1e-6

			m3 = t_blah*700*dS

			m3s.append(m3)

			#keyboard()
		
			#keyboard()
			tth.append(t_th)
			mth.append(m_th)
			V_vec.append(Vth)
			
			
			mdot_PAH = 0
			# PAH
			if PAH:
				
				mdot_cond = self.condensation('A4')*np.pi*self.WSR.Dc*dz
				mdot_conde.append(mdot_cond)
			
			
			
		return Sol.z, np.array(mth),np.array(hacas),np.array(masses2)#,mdot_conde
		#return Sol.z, np.array(mth),np.array(tth),hacas,mdot_conde


	def condensation(self,species):
		
		indx = self.gas.species_index(species)
		
		A3indx = self.gas.species_index('A3')
		A2indx = self.gas.species_index('A2')
		A1indx = self.gas.species_index('A1')
		
		Y = self.gas.Y[indx] 
		
		Y+= self.gas.Y[A3indx]+self.gas.Y[A2indx]
		
		A = 2.68713
		B = 1086.824
		C = -262.849
		
		T = 400
		
		Pv = 10**(A-(B/(T+C)))*100000
		
		Y_i  = Pv/self.gas.P
		
		Nu = self.WSR.Nu()
		
		
		hg = Nu*self.gas.thermal_conductivity/self.WSR.Dc
		Sc = self.gas.viscosity/(self.gas.density*self.gas.mix_diff_coeffs[indx])
		
		Pr = self.gas.cp*self.gas.viscosity/self.gas.thermal_conductivity
		
		K = (hg/(self.gas.density*self.gas.cp))*(Pr/Sc)**(2/3)
		
		
		mdot = K*self.gas.density*np.log((1-Y_i)/(1-Y))
		
		
		return mdot
	

	def integrate_gas(self,gas,tres):
		# Integrate Gas in Well Stirred Reactor for some extra residence time


		# Create CT reactor
		reactor = ct.IdealGasReactor(gas)

		# Create reactornetwork
		sim = ct.ReactorNet([reactor])

		# Advance in time
		sim.advance(tres)

	def V_Th(self,Twall,mdot):
		# Calculate Thermophoretic Drift Velocity
			
		# Calculate Talbot K
		Kth = self.K_talbot(self.gas,self.d_p)

		# Calculate Temperature Gradient
		Tgrad = self.calc_T_grad(Twall,mdot)

		return -Kth*self.gas.viscosity/self.gas.T*Tgrad
		
	def haca_growth(self,gas,Twall=False):
		
		# Create Surface Object
		soot = ct.Solution(self.surface_mechanism,'soot')
		# Create Interface (where all reaction rates are calculated)
		soot_interface = ct.Interface(self.surface_mechanism,'soot_interface',[gas,soot])

		#iFEstar = soot_interface.species_index('C_FE*')
		#iFE = soot_interface.species_index('C_FE-H')
		#iACstar = soot_interface.species_index('C_AC*')
		#iAC = soot_interface.species_index('C_AC-H')
		iC = soot_interface.species_index('Csoot-*')
		iCstar = soot_interface.species_index('Csoot-H')

		# Set Initial Site Distribution
		#x = np.zeros(4)
		#x[iFEstar] = 0.5
		#x[iACstar] = 0.5
		#x[iFE] = 0
		#x[iAC] = 0
		x = np.zeros(2)
		x[iC] = 0
		x[iCstar] = 1

		#if Twall:
		#	soot_interface.TP = Twall,gas.P
		#else:
		# Set Interface Properties (Assuming TP equal to gas) 
		soot_interface.TP = gas.TP
		soot_interface.X = x

		# Advance surface coverages to basically SS
		soot_interface.advance_coverages(1000)

		# This is Bulk soot index
		iCBCB3 = soot_interface.kinetics_species_index('CB-CB3')
		# Producution rate
		soot_dot = soot_interface.net_production_rates[iCBCB3]

		#keyboard()
		# Output total mass production rate
		return soot_dot*soot.molecular_weights[0]


	def ct_gas(self,gas_object,del_T):
	
		self.T_gas = gas_object.T
		
		#self.Kn =  
		#self.Re =
		#self.Nu = 
	
		# self.kg = 
		# self.kp
	
	def calc_T_grad(self,Twall,mdot):
	# Calculate Boundary Layer temperature gradient

		# Set WSR mdot to input mdot, 
		# The next part uses WSR methods to calculate stuff
		self.WSR.mdot = mdot
		return self.WSR.Nu()*(Twall-self.WSR.gas.T)/self.WSR.Dc
		
	def calc_Y_grad(self,Ysoot,dp):
		
		Db = self.D_brownian(self.gas,dp)

		return -self.WSR.Sh(DB=Db)*Ysoot/self.WSR.Dc

	def K_talbot(self,gas,dp):
	# Calculate Thermophoretic Constant via Talbot Model
	
		# Constants in Talbots Formulation 
		Cs = 1.147
		Cm = 1.146
		Ct = 2.20
		A1 = 1.2
		A2 = 0.41
		A3 = 0.88

		# Soot Particle Diameter
		#dp = self.MOMIC.soot_properties(self.WSR.M)[0]

		# Gas Thermal Conductivity
		kg = gas.thermal_conductivity
		# Calculate Kn
		Kn = 2*self.MOMIC.mean_free_path(gas.P,gas.T)/dp 

		# Cunningham correction factor
		Cc = 1+Kn*(A1+A2*np.exp(-A3/Kn))
		
		# Kth
		Kth = 2*Cs*Cc/(1+3*Cm*Kn)*((kg/self.kp)+Ct*Kn)/(1+2*kg/self.kp+2*Ct*Kn)
	
		return Kth

		#Vth = -Kth*self.nu/self.T_gas*self.del_T
	
	# Re/Nu calc

	def D_brownian(self,gas,dp):
		# Calculate Brownian Diffusion coefficient based on average soot properties

		kb = 1.38065e-23	# Boltzmann constant J/K

		# Knudsen #
		Kn = 2*self.MOMIC.mean_free_path(gas.P,gas.T)/dp 

		# Cunningham correction
		Cc = 1+Kn*(1.2+0.41*np.exp(-0.88/Kn))

		return kb*Cc*gas.T/(3*np.pi*gas.viscosity*dp)



	def mass_diffusion(self):
		pass
	def surf_reaction(self):
		pass
	

class Sol:
	# Having some pickling issues with Cantera SolutionArray

	def __init__(self):
		pass

	def ctSolArray(self,SolArray):

		keyboard()
		TDY = SolArray.TDY
		extra = SolArray._extra_arrays



class Sol_OLD:
	# Class for storing WSR/PFR Solution
	# Similar to cantera SolutionArray, I probably could have stolen that
	# But using this to practice with several different structures
	# This class ended up being an excercise in making trash
	# Not using it but keeping it incase I might need it for something
	
	# Went back to using cantera SolutionArray just to make stuff easier
	
	def __init__(self,gas,M,P,allocate_size = 100,**kwargs):
		# Initialize Sol Class
		
		# Store Moments and Gas T,P,Y

		# Initialize # Data Pts
		self.n_data = 0
		self.data_size = allocate_size

		# Initialize NAN array
		self.M = np.empty([allocate_size,np.size(M)])*np.nan
		self.P = np.empty([allocate_size,np.size(P)])*np.nan
		self.TPY = np.empty([allocate_size,2+np.size(gas.Y)])*np.nan

		# Variable Names 
		self.variables = ['M','P','TPY']

		for key in kwargs:
			setattr(self,key,np.empty([allocate_size,1])*np.nan)
			self.variables.append(key)
			
			#np.array(kwargs[key]))

		# Append Gas Object
		self.gas_object = gas


		self.add(gas,M,P,**kwargs)


		#self.M = np.array(M)
		#self.P = np.array(P)

		#self.gasT = np.array(gas.T)
		#self.gasP = np.array(gas.P)
		#self.gasY = np.array(gas.Y)
		#self.gas_object = gas

		#self.n_pts = 1

		
		
		# Add optional key words to Sol Attribute
		#for key in kwargs:
		#	setattr(self,key,np.array(kwargs[key]))
		#	self.variables.append(key)
	
	def add(self,gas,M,P,**kwargs):
		# Add data to Sol Class

		# If the number of data pts exceeds data size, then
		# extend data size
		if self.n_data > self.data_size:
			self.extend()

		# Add data
		self.M[self.n_data] = M
		self.P[self.n_data] = P
		self.TPY[self.n_data,:] = np.hstack([gas.T,gas.P,gas.Y])

		for key in kwargs:
			getattr(self,key)[self.n_data] = kwargs[key]

		# Increment # Data pts
		self.n_data += 1



	def extend(self,allocate_more = 100):
		# Extend data
		keyboard()
	
	def adjust_size(self,new_size):

		# Reduce size 
		for key in self.variables:
			setattr(self,key,getattr(self,key)[:new_size])
		
		self.data_size = new_size

	
	def append(self,gas,M,P,**kwargs):

		self.add(gas,M,P,**kwargs)
		return
		
		## OLD 
		
		self.M = np.vstack([self.M,M])
		self.P = np.vstack([self.P,P])

		self.gasT = np.vstack([self.gasT,gas.T])
		self.gasP = np.vstack([self.gasP,gas.P])
		self.gasY = np.vstack([self.gasY,gas.Y])
		
		for key in kwargs:
			value = getattr(self,key)
			setattr(self,key,np.vstack([value,kwargs[key]]))
	
	def gas(self,mix_prop):
		# Return gas mixture property
		
		values = np.zeros((np.size(self.gasT),1))

		#gas = lambda T,P,Y: self.gas.TPY = T,P,Y 
		def gas_prop(T,P,Y):
			self.gas_object.TPY = T,P,Y
			return getattr(self.gas_object,mix_prop)

		return np.c_[[gas_prop(self.gasT[i],self.gasP[i],self.gasY[i]) for i,_ in enumerate(self.gasT)]]
	
	def cleanup(self,i_start):
		
		self.adjust_size(self.n_data)

		for key in self.variables:
			getattr(self,key)[i_start:] = np.nan

		self.n_data = i_start

	
	def clean(self,i_start,i_end):
		# Purge Data from Sol Class

		self.cleanup(i_end)
		return
		#keyboard()

		for var in self.variables:
			value = getattr(self,var)
			setattr(self,var,value[i_start:i_end,:])
			
	def save_npy(self,file_name):
		# Save numpy dictionary array to file 
		
		#Create Dictionary of Values
		DATA = {}
		
		for var in self.variables:
			DATA[var] = getattr(self,var) 
		
		np.save(file_name,DATA)	

	@staticmethod
	def load_npy_sol(file_name,gas_object):	
		# Load previously saved npy while and create
		# Sol class

		# Load Data
		data = np.load(file_name,allow_pickle=True)
		dat_dict = data.flat[0]
		
		vars = [key for key in dat_dict.keys()]

		M = dat_dict['M']
		dat_dict.pop('M')
		vars.remove('M')
		P = dat_dict['P']
		vars.remove('P')
		dat_dict.pop('P')
		TPY = dat_dict['TPY']
		dat_dict.pop('TPY')
		vars.remove('TPY')

		gas_object.TPY = TPY[0,0], TPY[0,1], TPY[0,2:]
		additional_dat = {var:dat_dict[var][0] for var in vars}
		sol_object = Sol(gas_object,M[0],P[0],**additional_dat)

		for i,_ in enumerate(M[1:],start=1):
			gas_object.TPY = TPY[i,0], TPY[i,1], TPY[i,2:]
			additional_dat = {var:dat_dict[var][i] for var in vars}
			sol_object.add(gas_object,M[i],P[i],**additional_dat)
		
		sol_object.adjust_size(sol_object.n_data)

		return sol_object	