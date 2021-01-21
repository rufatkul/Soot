# Soot
Soot Object Oriented Toolbox

Written during my PhD to predict soot formation in a fuel rich kerosene rocket combustor

# Files:
- The MOMIC folder contains the main Method of Moments with Interpolative Closure Model written in Fortran and based on: Frenklach, M. “Method of Moments with Interpolative Closure.” Chemical Engineering Science, Vol. 57, No. 12, 2002, pp. 2229–2239. https://doi.org/10.1016/S0009-2509(02)00113-6.
  - This was originally written to be used in a CFD code but was never done, instead the code was used to calculate Moment source terms and integrated with SciPy
- The main folder contains momic_class.py which uses Cantera toolbox and calls the MOMIC module to calculate MOMENT source terms. This contains several classes including wrapper for the MOMENT code, a soot-coupled Well Stirred Reactor Model, and a soot-coupled Plug Flow Reactor Model 
- The main folder also has several examples of WSR and PFR calculations
 
# Prerequisites:
- Python 3.7 with the following packages:
  - Numpy
  - SciPy
  - Cantera 2.4
 - Fortran Compiler that works with f2py (Verified to work fine with Intel or GNU Fortran)
 
 # Instructions on setting up the code
 1. Run the make file in the MOMIC folder to create the MOMIC Module
 2. All set, setup any WSR/PFR calculations
 
 
 # Please Cite:
 @dissertation{Kulakhmetov2020,
author = "Rufat Kulakhmetov",
title = "{MEASUREMENT AND MODELING OF SOOT FORMATION AND DEPOSITION IN FUEL RICH HIGH PRESSURE KEROSENE COMBUSTION}",
year = "2020",
month = "12",
url = "https://hammer.figshare.com/articles/thesis/MEASUREMENT_AND_MODELING_OF_SOOT_FORMATION_AND_DEPOSITION_IN_FUEL_RICH_HIGH_PRESSURE_KEROSENE_COMBUSTION/13366178",
doi = "10.25394/PGS.13366178.v1"
}
 
