MODULE MOMIC_PARAMETERS
	IMPLICIT NONE
	
	INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=14)
	
	REAL(KIND=DBL) , PARAMETER :: ONE_SIXTH = 1.0_DBL/6.0_DBL
	REAL(KIND=DBL) , PARAMETER :: ONE_THIRD = 1.0_DBL/3.0_DBL
	REAL(KIND=DBL) , PARAMETER :: TWO_THIRDS = 2.0_DBL/3.0_DBL
	REAL(KIND=DBL) , PARAMETER :: ONE_HALF = 1.0_DBL/2.0_DBL
	REAL(KIND=DBL) , PARAMETER :: ZERO = 0.0_DBL
	REAL(KIND=DBL) , PARAMETER :: ONE = 1.0_DBL
	REAL(KIND=DBL) , PARAMETER :: TWO = 2.0_DBL

	REAL(KIND=DBL), PARAMETER :: kb = 1.38065D-16	 ! Boltzmann Constant cm^2 g s^-2 K^-1/molecule
	REAL(KIND=DBL), PARAMETER :: kb_SI = 1.38065D-23 ! Boltzmann Constant J/K
    REAL(KIND=DBL), PARAMETER :: PI = 3.14159D0		 ! Pi
	REAL(KIND=DBL), PARAMETER :: C_MASS = 12.0D0*1.67D-24 ! Mass of Carbon (g)
	REAL(KIND=DBL), PARAMETER :: rho_soot = 1.8D0 !1.86_DBL!*1E3 ! Soot Density (g/cm^3)
	REAL(KIND=DBL), PARAMETER :: R_u = 1.987D-3!8.314254D0/4184D0	 ! Universal Gas Constant (kcal/mol-K)
	REAL(KIND=DBL), PARAMETER :: Na = 6.022D23  ! Avogadros Number 
    
	!CD1 Parameter in Ranjans Code
	REAL(KIND=DBL), PARAMETER :: d_C = (6D0*C_MASS/PI/rho_soot)**(ONE_THIRD)	! C diameter in cm
	REAL(KIND=DBL), PARAMETER :: d_C_SI = d_C*1e-2                              ! C diameter in m 
    
	REAL(KIND=DBL), PARAMETER :: M_ref = 1           ! Divide All Equations by M_ref
	REAL(KIND=DBL), PARAMETER :: x_soot = 2.32D15	 ! Number Density of Csoot-H sites 1/cm^2 (Frenklach Wang 1994)
	REAL(KIND=DBL), PARAMETER :: x_soot_SI = 2.3D19	 ! Number Density converted to 1/m^2
	
	
	! Oxidation Parameters
	REAL(KIND=DBL), PARAMETER :: OXID_RAT = 32D0 ! Turn off Oxidation Source Terms if M(1)/M(0) < OXID_RAT

	! Use constant alpha or correlation
	LOGICAL, PARAMETER :: const_alpha = .FALSE.
	!Alpha Value to Use
	REAL(KIND=DBL), PARAMETER :: alpha_const = 1.0D0
	
    ! OH Mass
    REAL(KIND=DBL), PARAMETER :: OH_MASS = 17D0*1.67D-24
	! OH + SOOT oxidation constant
    REAL(KIND=DBL), PARAMETER :: CBOH = d_C**2*DSQRT(PI*kb/(2D0*OH_MASS))*Na
    
    
	! Aggregation Parameters (See Kazakov paper to see what these are)
	REAL(KIND=DBL), PARAMETER :: d_star = 10e-9!10e-9!1e-9!25e-9		 !
	REAL(KIND=DBL), PARAMETER :: D_f = 1.8
	
	! Initialize Moment Parameters from input file or call to initialize_momic
	LOGICAL, PARAMETER :: USE_INPUT_FILE = .FALSE.
	
	
END MODULE MOMIC_PARAMETERS
