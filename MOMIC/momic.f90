MODULE MOMIC
! Method of Moments with Interpolative Closure 
! Follows the formulation of:
! Frenklach, M. “Method of Moments with Interpolative Closure.” 
! Chemical Engineering Science, Vol. 57, No. 12, 2002, pp. 2229–2239.
! https://doi.org/10.1016/S0009-2509(02)00113-6.

! INCLUDES Kazakov Aggregation Model

! Written by: Rufat Kulakhmetov - 6/22/19 V0.1
! Painfully debugged - 7/21/19 - V0.9
! ADDED Free Molecular and Continuum Calculations - 11/7/19 - V0.98
! Fixed issue with Aggregate Coagulation - 3/23/20 - V0.99
! CHECK AGGREGATION COAG SUM
  
! Some Code Borrowed (Stolen) From Ranjan Mehta Moments Module

! CHECK MFP => USE GAS Properties instead of air?
!  - It doesn't seem to affect MOMENT calculations by much so using  air properties for simplicity
! CHECK VISCOSITY => USE GAS PROPERTIES!
!  - Same comment as above 

USE MOMIC_PARAMETERS

IMPLICIT NONE

! These Variables are used in the MOMIC module and need to be initialized
! Initialization Flag
LOGICAL :: isInitialized = .FALSE.
! MOMENT INPUT READ FLAG
LOGICAL :: inputRead = .FALSE.
! Binomial Coefficient
INTEGER, ALLOCATABLE :: binom(:,:)
! Reduced Moment Interpolants
REAL(KIND=DBL), ALLOCATABLE ::  L_pos(:,:), L_neg(:,:)
INTEGER :: mu_LO, mu_HI

! Data from Input File
! Input File Name
CHARACTER(LEN=20) :: input_file_name = 'momic.dat'
! Aggregation Flag (CONTROLLED IN INPUT FILE)
LOGICAL :: AGGREGATION = .FALSE.
! Calculate MOMIC Rates Flag (CONTROLLED IN INPUT FILE)
LOGICAL :: doSOOT = .FALSE.
! Number of M Moments (CONTROLLED IN INPUT FILE)
INTEGER :: M_Moments = 6	
! Number of P Moments (CONTROLLED IN INPUT FILE)
INTEGER :: P_Moments = 0	! Default Number of Aggregates
! Highest M and P Moment (CONTROLLED IN INPUT FILE)
INTEGER :: M_HI = 5 
INTEGER :: P_HI = 0
! Force REGIME CALCULATION
LOGICAL :: FORCE_REGIME = .FALSE.
LOGICAL :: FORCE_AGGREGATE = .FALSE.
LOGICAL :: OXID ! Oxidiation Flag
INTEGER :: REGIME 


! Coagulation Parameters
REAL(KIND=DBL) :: Kc, Kcp, Kf
! Soot Properties
REAL(KIND=DBL) :: D_soot

!Debugging Variables
INTEGER, SAVE :: NumCalls = 1

!PRIVATE

!PRIVATE :: NChooseK, COMPUTE_BINOM, L_fct, INITIALIZE_MOMIC
    CONTAINS
    
    ! READ INPUT FILE Properties
	SUBROUTINE READ_INPUT_FILE(input_file_name)
		IMPLICIT NONE
		CHARACTER(len=20), INTENT(IN) :: input_file_name
		INTEGER :: Number_of_moments, Number_of_aggregates
		INTEGER :: istat
		!REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: M, P 
		REAL(KIND=DBL), DIMENSION(20) :: M, P	! Temporary arrays, set to some arbitrary max dimension in order to read in all IC properly
		
		! NAMELIST FOR MOMIC Property FILE
		NAMELIST /Properties/ doSOOT, AGGREGATION, M_Moments, P_Moments
		! NAMELIST FOR MOMIC IC
		NAMELIST /IC/ M, P 
		! NAMELIST FOR COAGULATION REGIME
		NAMELIST /COAG_REGIME/ FORCE_REGIME, REGIME
		
		! Variable to check if this is first call 
		!LOGICAL, SAVE :: first_call = .TRUE.
	
		IF(inputRead .eqv. .FALSE.) THEN 
		
			OPEN(7,FILE=input_file_name,IOSTAT=istat,STATUS='OLD')
		
			! CHECK IF FILE OPENING WAS SUCCESFUL
			IF (istat /= 0) THEN
				WRITE(*,*) 'ERROR! MOMENT FILE NAME NOT FOUND: ',input_file_name
				STOP
			ELSE
				! READ Properties
				READ(UNIT=7,NML=Properties,IOSTAT=istat)
				WRITE(*,*) 'Number of Moments - Number of Aggregates'
				WRITE(*,*) M_Moments, P_moments
				
				! Set High M Moment
				M_HI = M_Moments-1
				! Set High P Moment
				P_HI = P_moments-1
				
				IF (M_Moments < 2) THEN
					WRITE(*,*) 'ERROR, THE NUMBER OF MOMENTS SPECIFIED IS <2'
					STOP
				END IF
				
				If (P_Moments > 1 .AND. AGGREGATION) THEN
					!ALLOCATE(P_0(1:Number_of_aggregates+1))
				ELSE
					AGGREGATION = .FALSE.
					WRITE(*,*) 'Aggregation Calculations are turned OFF'
				END IF							
		
				! READ COAGULATION REGIME IF Specified
				REWIND(UNIT=7)
				READ(UNIT=7,NML=COAG_REGIME,IOSTAT=istat)
								
				IF (istat /= 0) THEN
					FORCE_REGIME = .FALSE.
				ELSEIF (FORCE_REGIME .eqv. .FALSE.) THEN
					WRITE(*,*) 'COAGULATION REGIME WILL BE CALCULATED BASED ON Kn #'
					CONTINUE
				ELSEIF (FORCE_REGIME .eqv. .TRUE.) THEN
					IF (REGIME == 0) THEN
						WRITE(*,*) 'SETTING COAGULATION REGIME TO CONTINUUM'
					ELSEIF (REGIME == 1) THEN
						WRITE(*,*) 'SETTING COAGULATION REGIME TO TRANSITIONAL'
					ELSEIF (REGIME == 2) THEN
						WRITE(*,*) 'SETTING COAGULATION REGIEM TO FREE MOLECULAR'
					ELSE
						WRITE(*,*) 'UNKNOWN INPUT SPECIFIED FOR COAGULATION REGIME, &
								SETTING REGIME CALCULATIONS TO AUTOMATIC BASED ON Kn #'
						FORCE_REGIME = .FALSE.
					END IF
				ELSE 
					WRITE(*,*) "THERE IS SOME ERROR IN COAG_REGIME NAMELIST IN INPUT FILE"
					STOP
				END IF 
				
			END IF
			
			CLOSE(UNIT=7)
		
			!first_call = .FALSE.
			! Switch flag now that input is read
			inputRead = .True.
			
        ELSE
            WRITE(*,*) 'MOMIC IS ALREADY INITIALIZED'
        END IF 
		
    END SUBROUTINE READ_INPUT_FILE
    
	SUBROUTINE UPDATE_REGIME(SET_REGIME,FORCE_AGGREGATE_FLAG)
		! Update Regime Calculations During Calculation
		
		INTEGER, INTENT(IN) :: SET_REGIME
		INTEGER, INTENT(IN) :: FORCE_AGGREGATE_FLAG
	
		FORCE_REGIME = .TRUE.
		REGIME = SET_REGIME

		IF (REGIME == 3) THEN
			WRITE(*,*) 'SETTING COAGULATION REGIME TO CONTINUUM'
		ELSEIF (REGIME == 2) THEN
			WRITE(*,*) 'SETTING COAGULATION REGIME TO TRANSITIONAL'
		ELSEIF (REGIME == 1) THEN
			WRITE(*,*) 'SETTING COAGULATION REGIME TO FREE MOLECULAR'
		ELSE
			WRITE(*,*) 'COAGULATION REGIME WILL BE CALCULATED BASED ON Kn #'
			FORCE_REGIME = .FALSE.
		END IF
		
		IF(FORCE_AGGREGATE_FLAG==1) THEN
			WRITE(*,*) 'FORCING AGGREGATE REGIME CALCULATION'
			FORCE_AGGREGATE = .TRUE.
		END IF
	
	END SUBROUTINE UPDATE_REGIME
	
	
	SUBROUTINE INITIALIZE(Num_M_Moments,Num_P_Moments,SET_REGIME)
		! Initialize MOMIC with this Subroutine instead of using input file
		
		IMPLICIT NONE
		! Number of M and P Moments
		INTEGER, INTENT(IN) :: Num_M_Moments, Num_P_Moments
		INTEGER, INTENT(IN) :: SET_REGIME
		
		! Set M/P Moments
		M_Moments = Num_M_Moments
		P_Moments = Num_P_Moments
		
		! Set High M Moment
		M_HI = Num_M_Moments-1
		! Set High P Moment
		P_HI = Num_P_moments-1
		
		! Set Soot Calc Flag
		doSOOT = .TRUE.
	
		
		IF (Num_M_Moments < 2) THEN
			WRITE(*,*) 'ERROR, THE NUMBER OF MOMENTS SPECIFIED IS <2'
			STOP
		END IF
		
		IF (Num_P_Moments > 0) THEN
			AGGREGATION = .TRUE.
			WRITE(*,*) 'TURNING ON AGGREGATION CALCULATION'

		ELSE
			AGGREGATION = .FALSE.
			WRITE(*,*) 'AGGREGATION CALCULATIONS ARE TURNED OFF'
		END IF
		
		WRITE(*,*) 'Number of Moments - Number of Aggregates'
		WRITE(*,*) M_Moments, P_moments	
		
		
		!Check if Regime is Specified
		FORCE_REGIME = .TRUE.
		REGIME = SET_REGIME

		IF (REGIME == 3) THEN
			WRITE(*,*) 'SETTING COAGULATION REGIME TO CONTINUUM'
		ELSEIF (REGIME == 2) THEN
			WRITE(*,*) 'SETTING COAGULATION REGIME TO TRANSITIONAL'
		ELSEIF (REGIME == 1) THEN
			WRITE(*,*) 'SETTING COAGULATION REGIME TO FREE MOLECULAR'
		ELSE
			WRITE(*,*) 'COAGULATION REGIME WILL BE CALCULATED BASED ON Kn #'
			FORCE_REGIME = .FALSE.
		END IF
				
		! Change Flag now that input is read
		inputRead = .TRUE.

		! Initialize MOMIC
		CALL INITIALIZE_MOMIC()
		!ALLOCATE(mu(mu_LO:mu_HI))
		
	END SUBROUTINE INITIALIZE
		
	
	  	
    SUBROUTINE CALCULATE_SOURCE(M,P_in,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH, &
        W,G,R,H,Ragg,rC2H2, rCO, rH, rH2, rH2O, rO2, rOH)
        ! Calculate Moment Source Terms 
        IMPLICIT NONE
    
        REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M, P_in
        REAL(KIND=DBL), INTENT(IN) :: TEMPERATURE, PRESSURE
		REAL(KIND=DBL), INTENT(IN) :: C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH
        REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(M)) :: W, G, R
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(P_in)) :: H, Ragg
        REAL(KIND=DBL), INTENT(OUT) :: rC2H2, rCO, rH, rH2, rH2O, rO2, rOH
		REAL(KIND=DBL) :: rC2H2_nuc, rC2H2_surf
		
        ! Local Variables
        REAL(KIND=DBL), ALLOCATABLE :: mu(:)
        REAL(KIND=DBL), DIMENSION(SIZE(P_in)+1) :: P
		!REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: P 
        INTEGER :: i                         
		LOGICAL :: error_flag = .FALSE. ! Error Flag
		REAL(KIND=DBL) :: kn ! Knudsen Number
		INTEGER :: COAGULATION_REGIME  ! Free-Molecular/Transitional/Continuum Regime (Determined by Knudsen #)
		
		REAL(KIND=DBL), DIMENSION(SIZE(M)) :: f
		REAL(KIND=DBL) :: a,b,rsq 
		
        REAL(KIND=DBL), DIMENSION(SIZE(M)) :: W_C2H2, W_O2, W_OH
                 
        ! Initialize Moments Module:
        ! Calculate Binomial Coefficient
        ! Calculate Lagrange Interpolation Stencil
        ! Calculate Number of P and R moments
        ! Calculate lowest fraction and highest fraction moments
        IF(.NOT. isInitialized) THEN
            ! Read Moments Input
            CALL INITIALIZE_MOMIC()
		
			! Check if correct initial conditions are being used
			!IF(ALL(M==M_0)) THEN
			!	CONTINUE
			!ELSE
			!	WRITE(*,*) 'PROBLEM READING INITIAL CONDITIONS, THE M Moment at initial call to CALCULATE_SOURCE &
			!	is different from what is specified in input IC'
			!	WRITE(*,*) 'M at initial Call: ', M
			!	WRITE(*,*) 'M in IC file: ', M_0 
			!	error_flag = .TRUE.
            !END IF
				
            ! Check P moments if it's allocated
            !IF (CALC_AGGREGATE) THEN
			!    IF(ALL(P==P_0)) THEN 
			!	    CONTINUE
			 !   ELSE
			!	    WRITE(*,*) 'PROBLEM READING INITIAL CONDITIONS, THE P Moment at initial call to CALCULATE_SOURCE &
			!	    is different from what is specified in input IC'
			!	    WRITE(*,*) 'P at initial Call: ', P
			!	    WRITE(*,*) 'P in IC file: ', P_0  
			!	    error_flag = .TRUE.
             !   END IF
            !END IF 
			
			IF(error_flag) THEN
				WRITE(*,*) 'QUITING!!!'
				STOP
			END IF			
		
        END IF
    
	
        ! ************ MOMENT ERROR CHECKING ************
		! Check if P(0) = M(0) within some error tolerance
		!IF(ABS(P(1) - M(1))>1e-8 .AND. AGGREGATION_CALC) THEN 
		!	WRITE(*,*) 'ERROR: P(0):',P(1), ' /= M(0):', M(1)
		!	STOP
        !END IF
	    ! Don't need this, just set P(0) to M(0)
		
		
		!ALLOCATE(P(1:P_Moments))
		
		If (ANY(M(1) > P_in)) THEN
			P(:) = M(1)!1.D0
		ELSE
			! P(0) = M(0)
			P(1) = M(1)
			! Input Moments
			P(2:) = P_in
        END IF 
		
		!IF(present(AGG)) THEN
		!	IF(AGG==1) THEN
		!		AGGREGATION = .FALSE.
		!	ELSEIF(AGG==2) THEN
		!		AGGREGATION = .TRUE.
		!	END IF 
		!END IF
		
		
        ! Initialize Nucleation, Coagulation, Surface Rates, And Aggregate Rates
	    W(:) = 0
	    G(:) = 0
	    R(:) = 0
        H(:) = 0
        Ragg(:) = 0
        rC2H2_nuc = 0
		rC2H2_surf = 0
		rC2H2 = 0 
        rCO = 0
        rH = 0
        rH2 = 0
        rH2O = 0
        rO2 = 0
        rOH = 0
		OXID = .TRUE.
		mu(:) = 0
		
		
		
        ! Check if MOMIC Calculation is enabled
        IF (.NOT. doSOOT) THEN
            WRITE(*,*) 'MOMENT CALCULATION TURNED OFF IN MOMIC INPUT FILE'
            RETURN
        END IF
		
		! Check Inputs are correct
		! THIS CAN CAUSE POTENTIAL MEMORY ISSUES IF THEY ARE NOT THE SAME
		! MAYBE LOOK AT EDITING HOW MOMENTS ARE INPUT?
		IF (SIZE(M) /= M_MOMENTS) THEN
			WRITE(*,*) 'ERROR: THE NUMBER OF INPUT M MOMENTS DOES NOT MATCH UP &
			WITH THE NUMBER OF MOMENTS SPECIFIED IN THE INITIALIZATION'
			STOP
		END IF
		
		IF (AGGREGATION .AND. SIZE(P) /= P_MOMENTS) THEN
			WRITE(*,*) 'ERROR: THE NUMBER OF INPUT P MOMENTS DOES NOT MATCH UP &
			WITH THE NUMBER OF MOMENTS SPECIFIED IN THE INITIALIZATION'
			STOP
		END IF
		    
        ! Calculate Fractional MOMENTS
		ALLOCATE(mu(mu_LO:mu_HI))
		
        ! Interpolate reduced moments
		CALL interpolate_mu(M,mu,mu_LO,mu_HI)		
		
		
		! Oxidation Calculation Control Flag
		IF (mu(6) <= OXID_RAT) THEN
			OXID = .FALSE.
		END IF
		
		
        ! Calculate Kc, Kc' and Kf Coagulation Parameters
		CALL COAGULATION_PARAMETERS(TEMPERATURE,PRESSURE)
		
		! Calculate Average soot diameter 
		D_soot = Soot_D_Calc(mu,mu_LO,mu_HI)
		
		! Calculate Knudsen Number
		kn = 2*MFP(PRESSURE,TEMPERATURE)*(1e-2)/D_soot
					
		! Kn number regimes
		IF (kn < 0.1) THEN
			COAGULATION_REGIME = 3 ! Continuum
		ELSEIF (kn < 10) THEN
			COAGULATION_REGIME = 2 ! Transitional
		ELSEIF (kn >= 10) THEN
			COAGULATION_REGIME = 1 ! Free Molecular
		ELSEIF (isNAN(kn)) THEN
			WRITE(*,*) "ERROR IN Temperature or Pressure input"
		ELSE
			WRITE(*,*) 'kn#: =', kn
			WRITE(*,*) "HOW DID YOU GET TO THIS CONDITION??? THIS SHOULD BE PHYISCALLY IMPOSSIBLE!!"
			STOP
		END IF 
		
		  		  
		IF (FORCE_REGIME) THEN
			COAGULATION_REGIME = REGIME
		END IF         

        ! Calculate Primary Particle Nucleation Rates
		CALL NUCLEATION(M,TEMPERATURE,PRESSURE,C_C2H2,rC2H2_nuc,R,M_Moments)
		
		! If Number density is less than 1/cm^3 then only use Nucleation calculation,  skip the rest of calculation
		! Like in Ranjan's Code, this helps with some cases!!
		IF(M(1)<150D0) THEN
            WRITE(*,*) 'Soot Number Density M(0) < 150'
            WRITE(*,*) 'Only Calculating Nucleation Source'
			RETURN
		END IF 
					
		! Check if Aggregation is turned on
		
		!AGGREGATION = .FALSE.
		
		
		
        If(AGGREGATION) THEN
            CALL CALCULATE_SOURCE_AGG
        ELSE
            CALL CALCULATE_SOURCE_NO_AGG
        END IF
		
		! Get Total Rate of Nucleation and Surface Growth C2H2
		rC2H2 = rC2H2_nuc+rC2H2_surf
		

        ! ********************************************************
        ! ************* ERROR CHECKING ***************************
        ! ********************************************************
		! Check if rates are nan
        DO i = 1,M_Moments
			IF( isNAN(G(i)) .OR. isNAN(W(i)) .OR. isNAN(R(i))) THEN
				WRITE(*,*) 'NAN for Output Rates'
				!WRITE(*,*) 'R -- Gr -- Wr -- Rr'
				!WRITE(*,*) i, G(i), W(i), R(i)
				error_flag = .TRUE.
				!EXIT
            END IF
        END DO 
        
        ! Check if M(n) < M(0)
        !IF(ANY(M(2:)<M(1))) THEN
         !   WRITE(*,*) 'Error higher moments are smaller than 0th moment'
          !  error_flag = .True.
        !END IF
        
		CALL VERIFY_DISTRIBUTION(M,error_flag)
		
	
        ! DEBUGGING STUFF
        numCalls = numCalls + 1
		

		! Check Error Flag and Stop if error
		IF(error_flag) THEN
            WRITE(*,*) 'NumCalls = ',NumCalls
			WRITE(*,*) '*******INPUTS********'
			WRITE(*,*) 'M Moments'
			WRITE(*,"(6ES14.3)") M
			WRITE(*,*) 'P Moments'
			WRITE(*,"(6ES14.3)") P
			WRITE(*,*) 'Temperature - Pressure '
			WRITE(*,*) TEMPERATURE, PRESSURE
			WRITE(*,*) '[C2H2] - [H] - [H2]'
			WRITE(*,*) C_C2H2, C_H, C_H2
			WRITE(*,*) '[H2O] - [O2] - [OH]'
			WRITE(*,*) C_H2O, C_O2, C_OH
			WRITE(*,*) '******Calculations******'
			WRITE(*,*) 'mu'
			WRITE(*,*) mu
			WRITE(*,*) '---- Gr Rates --- '
			WRITE(*,*) G
			WRITE(*,*) '---- Wr Rates --- '
			WRITE(*,*) W
			WRITE(*,*) '---- Rr Rates --- '
			WRITE(*,*) R			
        END IF
		
		! f2py is not compiling if STOP is in the above if statement for some reason
		!IF(error_flag) THEN
		!	STOP
		!END IF 
		
        ! Normalize all Equations by M_ref
		! Suggested by Ranjan to help with Convergence
        R = R/M_ref
        G = G/M_ref
        W = W/M_ref
        Ragg = Ragg/M_ref
        H = H/M_ref

        ! Don't need this anymore so throw it away    		 
		DEALLOCATE(mu)
		
        RETURN
		
    CONTAINS
		SUBROUTINE VERIFY_DISTRIBUTION(M,error_flag)
			! Check Realizability of Distribution Function
			
			IMPLICIT NONE
			REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
			LOGICAL, INTENT(OUT) :: error_flag
			
			
			! Initialize error_flag
			error_flag = .FALSE.
			
			IF(M(0)<1) THEN
				WRITE(*,*) 'ERROR IN DISTRIBUTION: M(0) < 1'
				error_flag = .TRUE.
			ELSEIF(M(0)*M(2)/(M(1)*M(1))<1) THEN
				WRITE(*,*) 'ERROR in Distribution: M(0)M(2)/M(1)^2 < 1'
				error_flag = .TRUE.
			END IF
		END SUBROUTINE
	
		SUBROUTINE PRIMARY_COAGULATION(M,mu,mu_LO,mu_HI,G)
		
			IMPLICIT NONE
			REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
			REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:SIZE(M)-1) :: G
		
		
			IF(COAGULATION_REGIME==3) THEN
				CALL PRIMARY_COAGULATION_C(M,mu,mu_LO,mu_HI,G)
				WRITE(*,*) 'Coagulation is being calculated for Continuum Regime'
			ELSEIF(COAGULATION_REGIME==2) THEN
				CALL PRIMARY_COAGULATION_TR(M,mu,mu_LO,mu_HI,G)
				WRITE(*,*) 'Coagulation is being calculated for Transitional Regime'
			ELSEIF(COAGULATION_REGIME==1) THEN
				CALL PRIMARY_COAGULATION_FM(M,mu,mu_LO,mu_HI,G)
				WRITE(*,*) 'Coagulation is being calculated for Free Molecular Regime'
			ELSE
				WRITE(*,*) "WRONG COAGULATION REGIME FLAG, I DID SOMETHTHING WRONG?"
				STOP
			END IF 
		
		END SUBROUTINE PRIMARY_COAGULATION
		
		SUBROUTINE PRIMARY_COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,G)
			! Calculate Primary Coagulation with aggregation depending on COAGULATION_REGIME Flag
			IMPLICIT NONE
			REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M, P 
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
			REAL(KIND=DBL), INTENT(IN), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
			REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(M)) :: G
			
			IF(COAGULATION_REGIME==3) THEN
				CALL PRIMARY_COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,G)
				WRITE(*,*) 'Coagulation is being calculated for Continuum Regime'
			ELSEIF(COAGULATION_REGIME==2) THEN
				CALL PRIMARY_COAGULATION_AGGREGATE_TR(M,P,mu,mu_LO,mu_HI,pii,G)
				WRITE(*,*) 'Coagulation is being calculated for Transitional Regime'
			ELSEIF(COAGULATION_REGIME==1) THEN
				CALL PRIMARY_COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,G)
				WRITE(*,*) 'Coagulation is being calculated for Free Molecular Regime'
			ELSE
				WRITE(*,*) "WRONG COAGULATION REGIME FLAG, I DID SOMETHTHING WRONG"
				STOP
			END IF 	
		
		END SUBROUTINE PRIMARY_COAGULATION_AGGREGATE
		
		SUBROUTINE COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,H)
			! Calculate Aggregate Coagulation depending on COAGULATION_REGIME Flag
			
			IMPLICIT NONE
			
			!REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: M, P 
			REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M,P
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(IN) :: mu 
			REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2), INTENT(IN) :: pii
			REAL(KIND=DBL), DIMENSION(1:SIZE(P)-1), INTENT(OUT) :: H
			
			IF(COAGULATION_REGIME==3) THEN
				CALL COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,H)
				WRITE(*,*) 'Coagulation is being calculated for Continuum Regime'
			ELSEIF(COAGULATION_REGIME==2) THEN
				CALL COAGULATION_AGGREGATE_TR(M,P,mu,mu_LO,mu_HI,pii,H)
				WRITE(*,*) 'Coagulation is being calculated for Transitional Regime'
			ELSEIF(COAGULATION_REGIME==1) THEN
				CALL COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,H)
				WRITE(*,*) 'Coagulation is being calculated for Free Molecular Regime'
			ELSE
				WRITE(*,*) "WRONG COAGULATION REGIME FLAG, I DID SOMETHTHING WRONG"
				STOP
			END IF 		
		
		END SUBROUTINE COAGULATION_AGGREGATE
	
        SUBROUTINE CALCULATE_SOURCE_NO_AGG()
                ! Calculate Source Terms without Aggregation
                IMPLICIT NONE
        
                WRITE(*,*) 'Calculating Moment Rates with out Aggregation'
                
				
				CALL PRIMARY_COAGULATION(M,mu,mu_LO,mu_HI,G)
				
				
				CALL SURFACE_GROWTH(M,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
									rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
				
				
				! If Oxidation Regime then base surface rates on M1, like in Ranjan's code

				IF(W(2) < 0) THEN
					
					WRITE(*,*) '*** IN OXIDATION REGIME, Basing Moments on M(1) ***'
					DO i = 1,SIZE(M)-1
						f(i) = log(mu(6*i))
					END DO
					
					
					CALL linear(SIZE(M)-1,f,a,b,rsq)
					
						
					DO i = 4,6*(SIZE(M)-1)-1,6
						mu(i) = EXP(a+b*i/6.D0)
					END DO 
					
					
					!WRITE(*,*) '*******ORIGINAL W********'
					!WRITE(*,*) W 
					
					
					CALL SURFACE_GROWTH(M,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
				rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
						
				
					DO i =2,SIZE(M)
						W(i) = W(i)*M(i)/M(1)/EXP(a+b*i)
					END DO
					

				
					!WRITE(*,*) '*****New W********'
					!WRITE(*,*) W
					!STOP
				
				END IF 
                
				
				! Turn off Source Terms Cutoff Criteria Like in S.P Roy, Haworth Paper
				! For Strongly Oxidizing Environment
				!IF (W(2)<0 .AND. M(2)<32*M(1)) THEN
				!	W(:) = 0
				!	RETURN
				!END IF 
				
				
            END SUBROUTINE CALCULATE_SOURCE_NO_AGG
            
            SUBROUTINE CALCULATE_SOURCE_AGG()
	            ! CALCULATE SOURCE TERMS FOR METHOD OF MOMENTS USING KAZAKOV AGGREGATION MODEL
	            ! ASSUMING Average Particle Diameter Dsoot > D*
                IMPLICIT NONE
            
                ! Local Variable for reduced P moment
                REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
                
                ! Initialize Variables
                pii(:,:) = 0
                
                WRITE(*,*) 'Calculating Moment Rates with Aggregation Enabled'
				
                IF(D_soot > d_star .OR. FORCE_AGGREGATE) THEN

	                ! Interpolate pi 
	                CALL interpolate_pii(P,pii)
                    
					
                    ! Calculate Aggregate Nucleation
                    Ragg(:) = R(1)
                    
                    ! If P terms are 0 then only use nucleation term and non-aggregate coagulation, skip rest of calculation
                    IF(ANY(P == 0)) THEN
                        CALL PRIMARY_COAGULATION(M,mu,mu_LO,mu_HI,G)
                        RETURN
                    END IF
                    
            
	                ! Calculate Coagulation
	                CALL PRIMARY_COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,G)

					! Calculate Aggregate Coagulation
					CALL COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,H)
			
					! Calculate Surface Growth 
					CALL SURFACE_GROWTH_AGG(M,P,pii,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
											rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
		
		
					
					! If Oxidation Regime then base surface rates on M1, like in Ranjan's code
					IF(W(2)<0) THEN
						DO i = 1,SIZE(M)-1
							f(i) = log(mu(3*i))
						END DO	
						
						CALL linear(SIZE(M)-1,f,a,b,rsq)
							
						DO i = 2,3*(SIZE(M)-1)-1,3
							mu(i) = EXP(a+b*i/3.D0)
						END DO 
						
						
						!WRITE(*,*) '*******ORIGINAL W********'
						!WRITE(*,*) W 
						
						CALL SURFACE_GROWTH_AGG(M,P,pii,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
											rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
					
					
						DO i =2,SIZE(M)
							W(i) = W(i)*M(i)/M(1)/EXP(a+b*i)
						END DO
						
						! Oxidation should also result in reduced number density as particles get oxidized
						! Original Formulation does not account for that, so using this to account for it
						! W0 = W_OXIDATION/Cmin
						! Alternatively maybe set M0 = M1/Cmin?
						

						W(1) = (W_O2(1)+W_OH(1))/10D0
						
						
						!WRITE(*,*) '*****New W********'
						!WRITE(*,*) W
						!STOP
					
					END IF 							
                                        
                ELSE
                    WRITE(*,*) 'Average Soot D < d*'
                    CALL CALCULATE_SOURCE_NO_AGG()
                END IF
        
				RETURN
            END SUBROUTINE CALCULATE_SOURCE_AGG
                        
	END SUBROUTINE CALCULATE_SOURCE
    
    SUBROUTINE INITIALIZE_MOMIC()
	    ! Initialize Moments Module
		
	    IMPLICIT NONE
		
	    ! Local Variables
	    INTEGER :: r_max, i
	    REAL(KIND=DBL), ALLOCATABLE :: M_neg(:), M_pos(:) 
		
	    ! Read Moments Input
		IF(USE_INPUT_FILE) THEN
			WRITE(*,*) "Reading 'momic.dat'"
			CALL read_input_file(input_file_name)
        ELSE
			IF(inputRead .eqv. .FALSE.) THEN
				WRITE(*,*) 'PLEASE INITIALIZE THE MOMENT CLASS with INITIALIZE subroutine'
				STOP
			END IF
		END IF
		
        ! Calculate Highest Moment Order (Both M,P Moments) and highest M Moment
	    r_max = MAX(M_Moments,P_Moments)
        
	    ! Allocate Binomial Coefficient
	    ALLOCATE(binom(0:r_max,0:r_max))
	    ! Calculate Binomial Coefficient
	    binom = COMPUTE_BINOM(r_max)
		
	    ! Reduced Moments range: 
        ! mu_LO/6 => mu_HI/6
	    mu_LO = -4 
	    mu_HI = 6*(M_HI)+1
				
	    ! Allocate Arrays to store L(x) stencil for interpolating reduced moments
	    ALLOCATE(L_neg(3,mu_LO:-1))
	    ALLOCATE(L_pos(0:(M_HI),0:mu_HI))
		
	    ! Allocate Temporary Arrays to store interpolant x values
	    ! USE 1st 3 whole order moments to interpolate negative reduced values
	    ALLOCATE(M_neg(3))
	    M_neg = (/0,1,2/)
		
	    ! Use all whole order moments to interpolate positive reduced values
	    ALLOCATE(M_pos(M_Moments))
	    M_pos = (/(i, i=0,M_HI)/)
	
	    ! Calculate L(x) for positive and negative moments
	    DO i = mu_LO,-1
		    L_neg(:,i) = L_fct(REAL(i,DBL)/REAL(6,DBL),M_neg)
	    END DO 
	
	    DO i = 0,mu_HI
		    L_pos(:,i) = L_fct(REAL(i,DBL)/REAL(6,DBL),M_pos)
        END DO 	
        
	    ! Throw away garbage
	    DEALLOCATE(M_neg,M_pos)
		
	    ! Change Flag to Initialized
	    isInitialized = .True.
		
    END SUBROUTINE INITIALIZE_MOMIC
	
	! Debug Checked
	INTEGER FUNCTION NChooseK(n,k)
	! Calculate Binomial Coefficient
	
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(IN) :: k
	INTEGER :: i
	REAL(KIND=DBL) :: temp  !Store Results in temporary real data, otherwise Binomial Calculation are not correct for some values 
							
	temp = 1 ! SET VALUE OF TEMP TO 1 HERE AND NOT IN DECLARATION, OTHERWISE THERES SOME WEIRD BUG
			 ! The "WEIRD BUG" APPARENTLY IS A FORTRAN FEATURE, where any value defined in DECLARATION AUTOMATICALLY 
			 ! HAS A SAVE OPTION APPLIED
	
	IF (k > n) THEN
		WRITE(*,*) "ERROR in binomial coefficient k>n: ",k,">",n
		STOP
	ELSEIF (k==0) THEN
		NChooseK = 1
		RETURN
	END IF

	DO i = k,1,-1
		temp = temp*(n-i+1)/i
	END DO
	
	! Output Integer result
	NChooseK = temp

	END FUNCTION NChooseK
	
	! Debug Checked
	FUNCTION COMPUTE_BINOM(r_max)
		IMPLICIT NONE
		! Speed up binomial calculation by saving coefficient values in memory
		INTEGER,INTENT(IN) :: r_max 		
		INTEGER, DIMENSION(0:r_max,0:r_max) :: COMPUTE_BINOM
	
		! Local Variables
		INTEGER :: r, k
	
		! Calculate Binomial coefficient the first time the function is run 
		DO r = 0, r_max
			DO k = 0,r_max
				IF (k>r) THEN
					CONTINUE
				ELSE
					COMPUTE_BINOM(r,k) = NChooseK(r,k)
				END IF
			END DO
		END DO
	
	END FUNCTION COMPUTE_BINOM
	
	! Debug Checked
	FUNCTION L_fct(x_pt,x_vals)
		IMPLICIT NONE
		! GENERATE L(x) stencil for Lagrange interpolation
		! Input interpolation pt and x values
		REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: x_vals
		REAL(KIND=DBL), INTENT(IN) :: x_pt
		REAL(KIND=DBL), DIMENSION(size(x_vals)) :: L_fct
		REAL(KIND=DBL) :: L
		INTEGER :: i, j
	
		! Initialize L(x)
		L_fct = ZERO
	
		! Calculate L(x)
		DO i = 1,size(x_vals)
			L = ONE
			DO j = 1,size(x_vals)
				IF (i /= j) L = L*((x_pt-x_vals(j))/(x_vals(i)-x_vals(j)))
            END DO
			L_fct(i) = L
		END DO
			
	END FUNCTION L_fct

	SUBROUTINE linear(npts, y, a,b, rsq)
		IMPLICIT NONE 
		INTEGER, INTENT(IN) :: npts
		REAL(KIND=DBL), INTENT(IN), DIMENSION(npts) :: y
		REAL(KIND=DBL), INTENT(OUT) :: a, b, rsq
		
		INTEGER :: n
		REAL(KIND=DBL) :: summation, x1, x2, y1, y2, D, ymean, ypred, ess, yss
		
		summation = 1.D0*npts
		x1 = 0.D0
		x2 = 0.D0
		y1 = 0.D0
		y2 = 0.D0
		
		DO n = 1,npts
			x1 = x1 + REAL(n,8)
			x2 = x2 + REAL(n,8)*REAL(n,8)	
			y1 = y1 + y(n)
			y2 = y2 + y(n) * REAL(n,8)
		END DO
		D = summation*x2-x1*x1
		a = (y1*x2-x1*y2)/D
		b = (summation*y2-y1*x1)/D
		ymean = y1/npts
		ess = 0.D0
		yss = 0.D0
				
		DO n =1,npts
			ypred = a + b*n
			ess = ess + (ymean-ypred)*(ymean-ypred)
			yss = yss + (y(n) -ymean)*(y(n)-ymean)
		ENDDO
		rsq = ess/yss
	END SUBROUTINE linear
			
	
	FUNCTION bm(mu_r,pi_r)
		! Calculate Binary Moments required for 2-D aggregation particle size distribution
		
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN) :: mu_r
		REAL(KIND=DBL), INTENT(IN) :: pi_r
		REAL(KIND=DBL) :: bm 
		
		! Approximation suggested by Kazakov
		bm = mu_r*pi_r
	
	END FUNCTION bm
	
	SUBROUTINE interpolate_mu(M,mu,mu_LO,mu_HI)
	! Interpolate Fractional Reduced Moments
	! Input: 
	! Vector of MOMENTS
	
	IMPLICIT NONE
	
	REAL(KIND = DBL), INTENT(IN), DIMENSION(0:) :: M
	INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
	REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(OUT) :: mu

	! Local Variables
	INTEGER :: i
	
	! Initialize mu
	mu(:) = 0
    
	! Calculate reduced moments for whole order moments
	DO i = 0, M_HI
		! IF M(1) is 0 then do this to prevent divide by 0
		IF(M(0) == 0) THEN
			WRITE(*,*) 'M(0) = 0, Setting all reduced moments to 0'
			mu(:) = 0.0_DBL 
			RETURN
		END IF
		! Otherwise calculate reduced moment for whole order moments
		mu(6*i) = M(i)/M(0)
    END DO

	! Calculate Reduced Moments
	! Negative Reduced Moments (USE first 3 whole order moments
	DO i = mu_LO, -1
		mu(i) = DEXP(DOT_PRODUCT(L_neg(:,i),LOG(M(0:2)/M(0))))
    END DO
	
	DO i = 0, mu_HI
	! Positive Reduced Momements (USE all Positive moments)	
		mu(i) = DEXP(DOT_PRODUCT(L_pos(:,i),LOG(M(:)/M(0))))
    END DO 
	
	END SUBROUTINE interpolate_mu
	
	
	SUBROUTINE NUCLEATION(M,T,P,Y_C2H2,rC2H2,R,M_Moments)
		IMPLICIT NONE
		REAL(KIND = DBL), INTENT(IN), DIMENSION(:) :: M	! MOMENTS
        INTEGER, INTENT(IN) :: M_Moments            ! NEED TO SPECIFY or f2py complains!
		REAL(KIND = DBL), INTENT(IN) :: T, P, Y_C2H2 ! Temperature, Pressure, Mass Fraction
		REAL(KIND = DBL), INTENT(OUT), DIMENSION(M_Moments) :: R ! Nucleation Source Term
		REAL(KIND = DBL), INTENT(OUT) :: rC2H2 ! C2H2 Consumption rate
		REAL(KIND = DBL) :: C_min, A_i, Ei_R
		INTEGER :: I
		
		! Initialize W
		R(:) = 0
		
		! NUCLEATION VIA C2H2 MODEL FROM:
		! Joint-scalar transported PDF modeling
		! of soot formation and oxidation, R.P. Lindstedt, S.A Louloudi 2005
		! CONSTANTS ARE TAKEN FROM PAPER
		
		C_min = 10.D0   		! Napthalene
		A_i = 0.63D4			! Pre-exponential CONSTANTS
		Ei_R = 21000.D0			! Activation Energy
		
		
		IF (Y_C2H2 <= 0) THEN
			R(:) = 0
			rC2H2 = 0
			WRITE(*,*) 'Y(C2H2) <= 0 so setting nucleation rate to 0' 
			RETURN
        ELSE
			R(1) = 2.D0*(Na/C_min)*(A_i*EXP(-Ei_R/T))*Y_C2H2

			DO i = 2,SIZE(M)
				!R(i) = R(i-1)
				R(i) = R(i-1)*C_min/2.D0 
			END DO
			
			rC2H2 = -R(1)*C_min/Na	! Mol/cm^3-s 
			
        END IF
        
        ! DEBUGGING STUFF
        !IF(NumCalls==501) THEN
        !    WRITE(*,*) R/R(1)
        !    STOP
        !END IF
    
	END SUBROUTINE Nucleation
	
	FUNCTION MFP(Pin,T)
		! Calculate Mean Free Path for air, like in Ranjan's Code:
		! Mean free path for air in cm (from AK, myf 1997)
		
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN) :: Pin, T 
		REAL(KIND=DBL) :: MFP
		REAL(KIND=DBL) :: P
	
	
		! Convert Pa from Pa to atm
		P = Pin/101325D0
	
		! Correlation requirtes P in atm; T in K
		
		!                   This term \/ is(Pa=>Atm)
		MFP = 2.3701D-03 * T / P / 1.103D+05 

	END FUNCTION MFP
	
	FUNCTION L_LOG_INTERP(x_pt,x_vals,y_vals)
		! Calculate Lagrange log interp for a specific point
		
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN) :: x_pt
		REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: x_vals, y_vals
		REAL(KIND=DBL) :: L_LOG_INTERP
		
		L_LOG_INTERP = DEXP(DOT_PRODUCT(L_FCT(x_pt,x_vals),DLOG(y_vals)))
	
	END FUNCTION L_LOG_INTERP
		
	
	SUBROUTINE PRIMARY_COAGULATION_FM(M,mu,mu_LO,mu_HI,G)
		! Calculate coagulation term in method of moments
		! FOR FREE MOLECULAR COAGULATION WITH NO AGGREGATE
		! USING THE APPROACH IN:
		! Frenklach, M. “Method of Moments with Interpolative Closure.” 
		! Chemical Engineering Science, Vol. 57, No. 12, 2002, pp. 2229–2239. 
		! doi:https://doi.org/10.1016/S0009-2509(02)00113-6.

		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:SIZE(M)-1) :: G
		
		! Local Variables
		INTEGER :: r, k, n
		REAL(KIND=DBL) :: coag_sum 
		REAL(KIND=DBL), DIMENSION(4) :: f
		REAL(KIND=DBL) :: crk0, crk2, crk3, crk4A, crk4B, crk5A, crk5B, M02
			
		! Initialize Coagulation TERM
		G(:) = 0_DBL
	
		!! Calculate Grid Function
		
		!! MOMENT 0 COAGULATION 
		DO  n=1,MIN(M_Moments,4)
			f(n) = gridFun(n-1,0,0,mu,mu_LO,mu_HI)
		END DO 
		
		crk0 = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0,3D0/),f)
			
		
		!! MOMENT 1 COAGULATION = 0 
		
		IF(M_MOMENTS > 2) THEN
		
			DO n = 1,MIN(M_Moments,4)
				f(n) = gridFun(n-1,1,1,mu,mu_LO,mu_HI)
			END DO
			
			crk2 = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0,3D0/),f)
		
			IF (M_MOMENTS > 3) THEN
		
				DO n = 1,MIN(M_Moments,4)
					f(n) = gridFun(n-1,1,2,mu,mu_LO,mu_HI)
				END DO
				
				crk3 = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0,3D0/),f)
		
				IF (M_MOMENTS > 4) THEN
					
					DO n = 1,3
						f(n) = gridFun(n-1,1,3,mu,mu_LO,mu_HI)
					END DO
					
					crk4A = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0/),f(1:3))
					
					DO n = 1,4
						f(n) = gridFun(n-1,2,2,mu,mu_LO,mu_HI)
					END DO
					
					crk4B = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0,3D0/),f)
				
					IF (M_MOMENTS > 5) THEN
					
						DO n = 1,2
							f(n) = gridFun(n-1,1,4,mu,mu_LO,mu_HI)
						END DO
						
						crk5A = kf* L_LOG_INTERP(0.5D0, (/0D0, 1D0/),f(1:2))
						
						DO n = 1,3
							f(n) = gridFun(n-1,2,3,mu,mu_LO,mu_HI)
						END DO 
	
						crk5B = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0, 2D0/),f(1:3))
	
						IF (M_MOMENTS > 6) THEN
							WRITE(*,*) 'ERROR: COAGULATION FOR MOMENTS > 6 NOT IMPLEMENTED'
							STOP
						END IF
					END IF 
				END IF
			END IF
		END IF 
		
		M02 = M(0)*M(0)
		
		G(0) = -0.5D0*crk0*M02
		G(1) = 0.D0	
				
		IF (M_MOMENTS > 2) THEN
			G(2) = crk2*M02
			IF (M_MOMENTS>3) THEN
			G(3) = 3.D0*crk3*M02
				IF (M_MOMENTS>4) THEN
				G(4) = (4.D0*crk4A + 3.D0*crk4B)*M02
					IF (M_MOMENTS > 5) THEN
						G(5) = (5.D0*crk5A + 10.D0*crk5B)*M02
					END IF
				END IF
			END IF
		END IF
			
			
		CONTAINS
		
		FUNCTION gridFun(k,n,m,mu,mu_LO,mu_HI)
		!  This function used in calculation of free-molecular collision
		!  rates, which, for any (n,m), are determined by Lagrange
		!  interpolation using two or more gridFun(k,n,m)
		
			REAL(KIND=DBL) :: gridFun
			INTEGER :: k,n,m
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		
			! Local Variables
			INTEGER :: i,j,l
			
			gridFun = 0.D0
			
			DO l=0,k
			
				i = 6*(k-l+n)
				j = 6*(l+m)
				
				gridFun = gridFun + binom(k,l) * ( 	&
				mu(i+1)*mu(j-3) + 2D0*mu(i-1)*mu(j-1) + mu(i-3)*mu(j+1))
			
			END DO 
	
	
		END FUNCTION gridFun
	
	END SUBROUTINE PRIMARY_COAGULATION_FM
		
	
	SUBROUTINE PRIMARY_COAGULATION_TR(M,mu,mu_LO,mu_HI,G)
		! Calculate coagulation term in method of moments
		! FOR TRANSITIONAL COAGULATION WITH NO AGGREGATE
		! USING THE APPROACH IN:
		! Frenklach, M. “Method of Moments with Interpolative Closure.” 
		! Chemical Engineering Science, Vol. 57, No. 12, 2002, pp. 2229–2239. 
		! doi:https://doi.org/10.1016/S0009-2509(02)00113-6.
				
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:SIZE(M)-1) :: G

		! Local Variables
		REAL(KIND=DBL), DIMENSION (0:SIZE(M)-1) :: Gc, Gf
		
		! Initialize G Variable 
		G(:) = 0_DBL
		
		! Calculate Free Molecular Coagulation Rates
		CALL PRIMARY_COAGULATION_FM(M,mu,mu_LO,mu_HI,Gf)
		
		! Calculate Continuum Coagulation Rates
		CALL PRIMARY_COAGULATION_C(M,mu,mu_LO,mu_HI,Gc)
			
		! Calculate Transitional Coagulation Rates
		G(0) = (Gf(0)*Gc(0))/(Gf(0) + Gc(0))
		
		G(2:) = (Gf(2:)*Gc(2:))/(Gf(2:) + Gc(2:))
		
	
	END SUBROUTINE PRIMARY_COAGULATION_TR
		
		
	SUBROUTINE PRIMARY_COAGULATION_C(M,mu,mu_LO,mu_HI,G)
		! Calculate coagulation term in method of moments
		! FOR CONTINUUM ONLY WITH NO AGGREGATE
		! Using the approach in:
		! Frenklach, M. “Method of Moments with Interpolative Closure.” 
		! Chemical Engineering Science, Vol. 57, No. 12, 2002, pp. 2229–2239. 
		! doi:https://doi.org/10.1016/S0009-2509(02)00113-6.
		
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:SIZE(M)-1) :: G
		
		! Local Variables
		INTEGER :: r, k
		REAL(KIND=DBL) :: coag_sum 
		
		! Initialize Coagulation TERM
		G(:) = 0_DBL
		
		! Calculate Kc and Kcp terms in Coagulation Relations
		!Kc = (2.0_DBL*kb*T)/(3.0_DBL*viscosity(T))       
		!Kcp = 2.514_DBL*MFP(P,T)*(PI*rho_soot/(6*m_C))**(1/3) 
		
		
		! Calculate G0 => G(1)
		
		G(0) = -Kc*M(0)*M(0)*(1+mu(2)*mu(-2)+Kcp*(mu(-2)+mu(2)*mu(-4)))		
		
		DO r = 2,M_HI
			coag_sum = 0
                        
			DO k = 1, r-1
				 coag_sum = coag_sum + binom(r,k)*(2*mu(6*k)*mu(6*r-6*k)+mu(6*k+2)*mu(6*r-6*k-2) + &
					mu(6*k-2)*mu(6*r-6*k+2) + Kcp*(mu(6*k-2)*mu(6*r-6*k)+ &
					mu(6*k)*mu(6*r-6*k-2) + mu(6*k+2)*mu(6*r-6*k-4) + &
					mu(6*k-4)*mu(6*r-6*k+2)))
			END DO
			
			! Shift G indice
			G(r) = ONE_HALF*Kc*M(0)*M(0)*coag_sum
        END DO
        
	END SUBROUTINE PRIMARY_COAGULATION_C
	
	SUBROUTINE PRIMARY_COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,G)
		! Calculate coagulation term in method of moments
		! FOR CONTINUUM ONLY WITH AGGREGATION
		! Using the approach in:
		
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M, P 
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(M)) :: G
		
		! Local Variables
		INTEGER :: r, k
		REAL(KIND=DBL) :: pii_2S,agg_sum 
		
		! Initialize Coagulation TERM
		G(:) = ZERO
		
		! Need to Calculate pi_+2S term for coagulate aggregation
		pii_2S = interpolate_pii_term(P,2_DBL*(1.0_DBL/D_f - ONE_THIRD))

		! Calculate G(0) 
		G(1) = -Kc*M(1)*M(1)*(1+bm(mu(2),pii(0,1))*bm(mu(-2),pii(0,-1)) + &
			Kcp*(bm(mu(-2),pii(0,-1))+bm(mu(2),pii(0,1))*bm(mu(-4),pii_2S)))
		
		
		! Calculate Gr
		DO r = 1,SIZE(M)-1
			agg_sum = 0 
			DO k = 1, r-1
				agg_sum = agg_sum + binom(r,k)*((2*mu(6*k)*mu(6*r-6*k) + &
				bm(mu(6*k+2),pii(0,1))*bm(mu(6*r-6*k-2),pii(0,-1)) + &
				bm(mu(6*k-2),pii(0,-1))*bm(mu(6*r-6*k+2),pii(0,1)) + &
				Kcp*(bm(mu(6*k-2),pii(0,-1))*mu(6*r-6*k) + mu(6*k)*bm(mu(6*r-6*k-2),pii(0,-1)) + &
				bm(mu(6*k+2),pii(0,1))*bm(mu(6*r-6*k-4),pii(0,-2)) + &
				bm(mu(6*k-4),pii(0,-2))*bm(mu(6*r-6*k+2),pii(0,1)))))	
			END DO
			G(r+1) = ONE_HALF*M(1)*M(1)*Kc*agg_sum 
		END DO
		
	END SUBROUTINE PRIMARY_COAGULATION_AGGREGATE_C
	
	SUBROUTINE PRIMARY_COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,G)
		! Calculate coagulation term in method of moments
		! FOR FREE MOLECULAR ONLY WITH AGGREGATION
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M, P 
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:SIZE(M)-1) :: G
		REAL(KIND=DBL), DIMENSION(MIN(M_MOMENTS,4)) :: f
		REAL(KIND=DBL) :: crk0, crk2, crk3, crk4A, crk4B, crk5A, crk5B, M02
		
		! Local Variables
		INTEGER :: r, k, n
		!REAL(KIND=DBL) :: pii_2S,agg_sum 
	

		! Initialize Coagulation TERM
		G(:) = ZERO
		
		! Need to Calculate pi_+2S term for coagulate aggregation
		!pii_2S = interpolate_pii_term(P,2*(1.0_DBL/D_f - ONE_THIRD))
		
		DO n = 1, MIN(M_moments,4)
			f(n) = gridFun(n-1,0,0,mu,mu_LO,mu_HI,pii,P_Moments)
		END DO
	
		crk0 = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0,3D0/),f)
		
	
		IF(M_MOMENTS > 2) THEN
			
			DO n = 1,MIN(M_MOMENTS,4)
				f(n) = gridFun(n-1,1,1,mu,mu_LO,mu_HI,pii,P_Moments)
			END DO
			
			crk2 = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0,3D0/),f)

			IF (M_MOMENTS > 3) THEN
				
				DO n = 1,MIN(M_Moments,4)
					f(n) = gridFun(n-1,1,2,mu,mu_LO,mu_HI,pii,P_Moments)
				END DO
				
				crk3 = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0,3D0/),f)

			
				IF (M_MOMENTS > 4) THEN
					
					DO n = 1,3
						f(n) = gridFun(n-1,1,3,mu,mu_LO,mu_HI,pii,P_Moments)
					END DO
					
					crk4A = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0/),f(1:3))

					DO n = 1,4
						f(n) = gridFun(n-1,2,2,mu,mu_LO,mu_HI,pii,P_Moments)
					END DO 
	
					crk4B = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0,2D0,3D0/),f)
					
					IF (M_MOMENTS > 5) THEN
					
						DO n = 1,2
							f(n) = gridFun(n-1,1,4,mu,mu_LO,mu_HI,pii,P_Moments)
						END DO

						crk5A = kf* L_LOG_INTERP(0.5D0, (/0D0, 1D0/),f(1:2))

						DO n = 1,3
							f(n) = gridFun(n-1,2,3,mu,mu_LO,mu_HI,pii,P_Moments)
						END DO 

						crk5B = kf*L_LOG_INTERP(0.5D0, (/0D0, 1D0, 2D0/),f(1:3))

						IF (M_MOMENTS > 6) THEN
							WRITE(*,*) 'ERROR: COAGULATION FOR MOMENTS > 6 NOT IMPLEMENTED'
							STOP
						END IF
					END IF 
				END IF 
			END IF 
		END IF 
		
		M02 = M(0)*M(0)
		
		G(0) = -0.5D0*crk0*M02
		G(1) = 0.D0
		
		IF (M_MOMENTS > 2) THEN
			G(2) = crk2*M02
			IF (M_MOMENTS>3) THEN
			G(3) = 3.D0*crk3*M02
				IF (M_MOMENTS>4) THEN
				G(4) = (4.D0*crk4A + 3.D0*crk4B)*M02
					IF (M_MOMENTS > 5) THEN
						G(5) = (5.D0*crk5A + 10.D0*crk5B)*M02
					END IF
				END IF
			END IF
		END IF
	
	CONTAINS
		FUNCTION gridFun(l,x,y,mu,mu_LO,mu_HI,pii,P_moments)
		! This is the grid function for calculation in free-molecular aggregate collision rates
		! Formulation follows Kazakov 1998
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: l,x,y
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		INTEGER, INTENT(IN) :: P_moments ! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), DIMENSION(0:P_moments-2,-2:2), INTENT(IN) :: pii ! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL) :: gridFun
		
		! Local variables
		INTEGER :: k 
		
		! Initialize gridFun
		gridFun = 0_DBL
		
		DO k = 0,l
		
		gridFun = gridFun + binom(l,k)*( &
		bm(mu(6*x+6*k+1),pii(0,2))*mu(6*y+6*l-6*k-3) &
		+ 2*(bm(mu(6*x+6*k-1),pii(0,1))*bm(mu(6*y+6*l-6*k-1),pii(0,1))) &
		+ mu(6*x+6*k-3)*bm(mu(6*y+6*l-6*k+1),pii(0,2)))
		
		END DO
		
		END FUNCTION gridFun
	
	END SUBROUTINE PRIMARY_COAGULATION_AGGREGATE_FM
	
	SUBROUTINE PRIMARY_COAGULATION_AGGREGATE_TR(M,P,mu,mu_LO,mu_HI,pii,G)
		! Calculate coagulation term in method of moments
		! FOR TRANSITIONAL COAGULATION WITH NO AGGREGATE
		! USING THE APPROACH IN:
		! Frenklach, M. “Method of Moments with Interpolative Closure.” 
		! Chemical Engineering Science, Vol. 57, No. 12, 2002, pp. 2229–2239. 
		! doi:https://doi.org/10.1016/S0009-2509(02)00113-6.
				
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M, P 
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:SIZE(M)-1) :: G
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:SIZE(P)-1,-2:2) :: pii

		! Local Variables
		REAL(KIND=DBL), DIMENSION (0:SIZE(M)-1) :: Gc, Gf
		
		! Initialize G Variable 
		G(:) = 0_DBL

		! Calculate Free Molecular Coagulation Rates
		CALL PRIMARY_COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,Gf)
		
		! Calculate Continuum Coagulation Rates
		CALL PRIMARY_COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,Gc)	
		
		! Calculate Transitional Coagulation Rates
		G(0) = (Gf(0)*Gc(0))/(Gf(0) + Gc(0))
		
		G(2:) = (Gf(2:)*Gc(2:))/(Gf(2:) + Gc(2:))			
	
	END SUBROUTINE PRIMARY_COAGULATION_AGGREGATE_TR
	
	SUBROUTINE COAGULATION_PARAMETERS(T,P)
		! Calculate Kc, Kc', Kcp constants for use in Coagulation
		
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN) :: T, P
		
		! Continuum Coagulation Parameters
		Kc = TWO_THIRDS*(kb*T)/(viscosity(T))  
		Kcp = 2.514_DBL*MFP(P,T)*(PI*rho_soot/(6_DBL*C_MASS))**(ONE_THIRD) 
			
		! Free Molecular Coagulation Parameter
		Kf = 2.2D0*DSQRT(6D0*kb/rho_soot*T)*(3D0/4D0*(C_MASS/(PI*rho_soot)))**ONE_SIXTH
		
        	
	END SUBROUTINE COAGULATION_PARAMETERS
	
    SUBROUTINE SURFACE_GROWTH(M,mu,mu_LO,mu_HI,T,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
		rC2H2,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
		
		IMPLICIT NONE
		
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu 
		REAL(KIND=DBL), INTENT(IN) :: T, C_C2H2, C_H, C_H2, C_H2O, C_O2, C_OH
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:) :: W_C2H2, W_O2, W_OH, W
		REAL(KIND=DBL), INTENT(OUT) :: rC2H2, rCO, rH, rH2, rH2O, rO2, rOH
        
		! Local Variables
		REAL(KIND=DBL) :: RT, k1f, k1r, k2f, k2r, k3f, k4f, k5f, k6, &
							alpha, a, b, Cs, r, x_soot_star, C1, C2, denom
					
        
		! INITIALIZE Variables
		W_C2H2(:) = 0
		W_O2(:) = 0
		W_OH(:) = 0
		W(:) = 0
		rC2H2 = 0
		rCO = 0
		rH = 0
		rH2 = 0
		rH2O = 0
		rO2 = 0
		rOH = 0
		
		
		! Reaction Rate Coefficients from Appel, Bockhorn, Frenklach 1999
		RT = R_u*T
        k1f = 4.17D13*DEXP(-13.D0/(RT))
		k1r = 3.9D12*DEXP(-11.D0/(RT))
		k2f = 1.0D10*T**(0.734D0)*DEXP(-1.43D0/(RT))
		k2r = 3.68D8*T**(1.139D0)*DEXP(-17.1D0/(RT))
		k3f = 2.0D13
		k4f = 8.0D7*T**(1.56D0)*DEXP(-3.8D0/(RT))
		k5f = 2.2D12*DEXP(-7.5D0/(RT))
		k6 = 0.13D0
                
        IF(.NOT. const_alpha) THEN
		    ! Calculate Fraction of Active Sites
		    ! ABF correlation 
		    a = 12.65D0 - 5.63D-3*T
		    b = -1.38D0 + 6.8D-4*T
		  
		    alpha = MAX(0.0D0,MIN(1.0D0,DTANH(a/DLOG10(mu(3))+b)))           
			
        ELSE
            alpha = alpha_const
        END IF
		
        ! ************ Temp Fix*******************
		!IF(alpha<0) THEN
		!	alpha = 0
		!END IF 
		
		
		
		! IF denominator 0 then return
		IF (k1r*C_H2+k2r*C_H2O+k3f*C_H+k4f*C_C2H2+k5f*C_O2 == 0) THEN
			RETURN
		END IF 
		
		! Number density of surface radicals (in 1/m^2)
		x_soot_star = x_soot*((k1f*C_H+k2f*C_OH)/&
								(k1r*C_H2+k2r*C_H2O+k3f*C_H+k4f*C_C2H2+k5f*C_O2))
		
        
		! Carbon diameter
		!Cs = (6*m_C/(PI*rho_soot))**(1/3)
		
		! Initialize Rates
		W_C2H2(0) = 0_DBL
		W_O2(0) = 0_DBL
		W_OH(0) = 0_DBL 

        
		DO r = 1,M_HI
			W_C2H2(r) = k4f*C_C2H2*alpha*x_soot_star*PI*d_C**2*M(0)*&
						surf_sum(mu,mu_LO,mu_HI,INT(r),2D0)
						
            W_O2(r) = k5f*C_O2*alpha*x_soot_star*PI*d_C**2*M(0)*&
						surf_sum(mu,mu_LO,mu_HI,INT(r),-2D0)
		
            W_OH(r) = k6*C_OH*CBOH*DSQRT(T)*M(0)*surf_sum(mu,mu_LO,mu_HI,INT(r),-1.0D0)!*alpha*x_soot_star*PI*d_C**2*M(0)*&
						!surf_sum(mu,mu_LO,mu_HI,INT(r),-1.0_DBL) 
        END DO
		
		
		

		! Calculate terms with or without oxidation or without
		IF(.NOT. OXID) THEN
			W_O2(:) = 0
			W_OH(:) = 0
		END IF 
		
		W = W_C2H2 + W_O2 + W_OH
        
        ! Calculate Gas Rates (UNITS: mol/cm^3 => MAKE SURE TO USE RIGHT UNITS IN CANTERA!)
        ! Define Constants 
        C1 = alpha*x_soot*PI*d_C**2*M(0)
        C2 = alpha*x_soot_star*PI*d_C**2*M(0)
        
        rC2H2 = -W_C2H2(1)/(2D0*Na)
        rO2 = W_O2(1)/(2D0*Na)
        rOH = ((k2r*C_H2O*C2-k2f*C_OH*C1)*mu(4) + W_OH(1))/Na
        rH = ((k1r*C_H2 +k4f*C_C2H2)*C2 - &
            k1f*C_H*C1 - k3f*C_H*C2)*mu(4)/Na     
        rH2 = (k1f*C_H*C1 - k1r*C_H2*C2)*mu(4)/Na
        rH2O = (k2f*C_OH*C1   - k2r*C_H2O*C2)*mu(4)/Na
        rCO = -(W_O2(1) + W_OH(1))/Na 
		 
		 
		! Surface Growth Summation Function
		CONTAINS 
		
			FUNCTION surf_sum(mu,mu_LO,mu_HI,r,delta)
			! Calculate Summation in Coagulate sum
			IMPLICIT NONE
			
			INTEGER, INTENT(IN) :: mu_LO,mu_HI
			REAL(KIND=DBL), INTENT(IN) :: delta
			REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(IN) :: mu
			INTEGER, INTENT(IN) :: r
			REAL(KIND=DBL) :: surf_sum
			
			! Local Variable
			INTEGER :: k
			
			! Initialize sum
			surf_sum = 0
			
			DO k = 0,r-1
				surf_sum = surf_sum + binom(r,k)*delta**(r-k)*mu(6*k+4)
            END DO
            
			END FUNCTION surf_sum
		
	
	END SUBROUTINE SURFACE_GROWTH
	
	SUBROUTINE SURFACE_GROWTH_AGG(M,P,pii,mu,mu_LO,mu_HI,T,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
		rC2H2,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
			
		IMPLICIT NONE
		
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: P
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
		
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu 
		REAL(KIND=DBL), INTENT(IN) :: T, C_C2H2, C_H, C_H2, C_H2O, C_O2, C_OH
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:) :: W_C2H2, W_O2, W_OH, W
		REAL(KIND=DBL), INTENT(OUT) :: rC2H2, rCO, rH, rH2, rH2O, rO2, rOH
        
		! Local Variables
		REAL(KIND=DBL) :: RT, k1f, k1r, k2f, k2r, k3f, k4f, k5f, k6, &
							alpha, a, b, Cs, r, x_soot_star, C1, C2, denom, pi_third
					
        
		! INITIALIZE Variables
		W_C2H2(:) = 0
		W_O2(:) = 0
		W_OH(:) = 0
		W(:) = 0
		rC2H2 = 0
		rCO = 0
		rH = 0
		rH2 = 0
		rH2O = 0
		rO2 = 0
		rOH = 0
		
		
		! Reaction Rate Coefficients from Appel, Bockhorn, Frenklach 1999
		RT = R_u*T
        k1f = 4.17D13*DEXP(-13.D0/(RT))
		k1r = 3.9D12*DEXP(-11.D0/(RT))
		k2f = 1.0D10*T**(0.734D0)*DEXP(-1.43D0/(RT))
		k2r = 3.68D8*T**(1.139D0)*DEXP(-17.1D0/(RT))
		k3f = 2.0D13
		k4f = 8.0D7*T**(1.56D0)*DEXP(-3.8D0/(RT))
		k5f = 2.2D12*DEXP(-7.5D0/(RT))
		k6 = 0.13D0
                
        IF(.NOT. const_alpha) THEN
		    ! Calculate Fraction of Active Sites
		    ! ABF correlation 
		    a = 12.65D0 - 5.63D-3*T
		    b = -1.38D0 + 6.8D-4*T
		  
		    alpha = MIN(1.0D0,DTANH(a/DLOG10(mu(3))+b))           
        ELSE
            alpha = alpha_const
        END IF
		
        ! ************ Temp Fix*******************
		!IF(alpha<0) THEN
		!	alpha = 0
		!END IF 
		
		
		
		! IF denominator 0 then return
		IF (k1r*C_H2+k2r*C_H2O+k3f*C_H+k4f*C_C2H2+k5f*C_O2 == 0) THEN
			RETURN
		END IF 
		
		! Number density of surface radicals (in 1/m^2)
		x_soot_star = x_soot*((k1f*C_H+k2f*C_OH)/&
								(k1r*C_H2+k2r*C_H2O+k3f*C_H+k4f*C_C2H2+k5f*C_O2))
		
        
		! Carbon diameter
		!Cs = (6*m_C/(PI*rho_soot))**(1/3)
		
		! Initialize Rates
		W_C2H2(0) = 0_DBL
		W_O2(0) = 0_DBL
		W_OH(0) = 0_DBL 

        
		DO r = 1,M_HI
			W_C2H2(r) = k4f*C_C2H2*alpha*x_soot_star*PI*d_C**2*M(0)*&
						surf_sum(mu,mu_LO,mu_HI,INT(r),2D0)
						
            W_O2(r) = k5f*C_O2*alpha*x_soot_star*PI*d_C**2*M(0)*&
						surf_sum(mu,mu_LO,mu_HI,INT(r),-2D0)
		
            W_OH(r) = k6*C_OH*CBOH*DSQRT(T)*M(0)*surf_sum(mu,mu_LO,mu_HI,INT(r),-1.0D0)!*alpha*x_soot_star*PI*d_C**2*M(0)*&
						!surf_sum(mu,mu_LO,mu_HI,INT(r),-1.0_DBL) 
        END DO
		
		
		pi_third = interpolate_pii_term(P,ONE_THIRD)
		
				
		! Enhancement due to aggregation
		W_C2H2 = W_C2H2*pi_third!*M(0)
		W_O2 = W_O2*pi_third!*M(0)
		W_OH = W_OH*pi_third!*M(0)
		
		
        W = W_C2H2 + W_O2 + W_OH
        
        
        ! Calculate Gas Rates (UNITS: mol/cm^3 => MAKE SURE TO USE RIGHT UNITS IN CANTERA!)
        ! Define Constants 
        C1 = alpha*x_soot*PI*d_C**2*M(0)
        C2 = alpha*x_soot_star*PI*d_C**2*M(0)
        
        rC2H2 = -W_C2H2(1)/(2D0*Na)
        rO2 = W_O2(1)/(2D0*Na)
        rOH = ((k2r*C_H2O*C2-k2f*C_OH*C1)*mu(4) + W_OH(1))/Na
        rH = ((k1r*C_H2 +k4f*C_C2H2)*C2 - &
            k1f*C_H*C1 - k3f*C_H*C2)*mu(4)/Na     
        rH2 = (k1f*C_H*C1 - k1r*C_H2*C2)*mu(4)/Na
        rH2O = (k2f*C_OH*C1   - k2r*C_H2O*C2)*mu(4)/Na
        rCO = -(W_O2(1) + W_OH(1))/Na 
		
		 
		! Surface Growth Summation Function
		CONTAINS 
		
			FUNCTION surf_sum(mu,mu_LO,mu_HI,r,delta)
			! Calculate Summation in Coagulate sum
			IMPLICIT NONE
			
			INTEGER, INTENT(IN) :: mu_LO,mu_HI
			REAL(KIND=DBL), INTENT(IN) :: delta
			REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(IN) :: mu
			INTEGER, INTENT(IN) :: r
			REAL(KIND=DBL) :: surf_sum
			
			! Local Variable
			INTEGER :: k
			
			! Initialize sum
			surf_sum = 0
			
			DO k = 0,r-1
				surf_sum = surf_sum + binom(r,k)*delta**(r-k)*mu(6*k+4)
            END DO
            
			END FUNCTION surf_sum
		
	
	END SUBROUTINE SURFACE_GROWTH_AGG

	
	
	
	PURE FUNCTION viscosity(T)
		IMPLICIT NONE
		! Stolen from Ranjan's Code
		!  Viscosity of air as a function of temperature (K)
		!  (from AK, myf 1997)
	
		REAL(KIND=DBL), INTENT(IN) :: T
		REAL(KIND=DBL) :: viscosity
	
		viscosity = 14.58D-06 * (T**1.5D0) / (T + 100.4D0)
		
	END FUNCTION viscosity 
	
	FUNCTION SOOT_D_CALC(mu,mu_LO,mu_HI)
		! Calculate Average Soot Diameter (in m)
		
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		REAL(KIND=DBL) :: SOOT_D_CALC
		        
		SOOT_D_CALC = d_C*mu(2)*1e-2
		
    END FUNCTION SOOT_D_CALC
    
    FUNCTION SOOT_fv_CALC(M)
        ! Calculate Soot Volume Fraction
        
        IMPLICIT NONE
        !INTEGER, INTENT(IN) :: mu_LO, mu_HI ! Need to specify these here otherwise f2py complains!
        !REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
        REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M 
        REAL(KIND=DBL) :: SOOT_fv_CALC
        
        SOOT_fv_CALC = M(1)*C_MASS/rho_soot
        
        SOOT_fv_CALC = MAX(0.0D0, SOOT_fv_CALC)
    
    END FUNCTION SOOT_fv_CALC
	
	FUNCTION SOOT_Spherical_SA_CALC(M,mu,mu_LO,mu_HI)
		! Calculate Soot Surface Area Density
		
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: mu_LO, mu_HI
		REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
		REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M
		REAL(KIND=DBL) :: SOOT_Spherical_SA_CALC
		
		SOOT_Spherical_SA_CALC = PI*(6*C_MASS/(rho_soot*PI))**(TWO_THIRDS)*mu(4)*M(1)
				
	END FUNCTION SOOT_Spherical_SA_CALC
	
	
	FUNCTION SOOT_Y_CALC(M,rho_gas)
		! Calculate Soot Mass Fraction
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M
		REAL(KIND=DBL), INTENT(IN) :: rho_gas	! Gas Density in kg/m^3
		!REAL(KIND=DBL) :: fv
		REAL(KIND=DBL) :: SOOT_Y_CALC
		!REAL(KIND=DBL) :: SOOT_Y_CALC2
		
		SOOT_Y_CALC = (C_MASS*M(2))/(rho_gas*1e-3+C_MASS*M(2))
		
		! fv = SOOT_fv_CALC(M)
				
		!SOOT_Y_CALC = (fv*rho_soot/rho_gas)/(fv*rho_soot/rho_gas+1)
		
	END FUNCTION SOOT_Y_CALC
	
	FUNCTION SOOT_Y_CALC_AGG(M,P,rho_gas)
		! Calculate Soot Mass Fraction for Aggregate
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: P
		REAL(KIND=DBL), INTENT(IN) :: rho_gas	! Gas Density in kg/m^3
		!REAL(KIND=DBL) :: fv
		REAL(KIND=DBL) :: SOOT_Y_CALC_AGG
		
		SOOT_Y_CALC_AGG = SOOT_Y_CALC(M,rho_gas)*P(1)/P(0)
	END FUNCTION SOOT_Y_CALC_AGG
		
        
    SUBROUTINE SOOT_PROPERTIES(M,D,fv,SA)
        ! Calculate Soot Properties from moments, which can be output to some other subroutine
        REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: M 
        REAL(KIND=DBL), INTENT(OUT) :: D,fv, SA
        
        ! Local Variables
        REAL(KIND=DBL), DIMENSION(-2:6*SIZE(M)+2) :: mu
        
		
		IF(isInitialized .eqv. .FALSE.) THEN
			WRITE(*,*) 'INITIALIZE MOMIC MODULE FIRST BEFORE CALLING SOOT_PROPERTIES!!!'
			STOP
		END IF
		
				
        CALL interpolate_mu(M,mu,mu_LO,mu_HI)
				
        ! Calculate Soot Average Diameter
        D = SOOT_D_CALC(mu,mu_LO,mu_HI)
		! Calculate Diameter of Average Soot Particle
		!Davg = d_C*M(2)/M(1) 
        ! Calculate Soot Volume Fraction
        fv = SOOT_fv_CALC(M)
        ! Calculate Spherical Surface Area
		SA = SOOT_Spherical_SA_CALC(M,mu,mu_LO,mu_HI)

    END SUBROUTINE SOOT_PROPERTIES
    
	SUBROUTINE SOOT_PROPERTIES_AGG(M,P,D,fv,SA)
        ! Calculate Soot Properties from moments, which can be output to some other subroutine
        REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: M
		REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: P
        REAL(KIND=DBL), INTENT(OUT) :: D, fv, SA		
	
		! Local Variables
		REAL(KIND=DBL), DIMENSION(-2:6*SIZE(M)+2) :: mu
		REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
		
		CALL interpolate_mu(M,mu,mu_LO,mu_HI)
		CALL interpolate_pii(P,pii)
		
		! Calculate Coalescant Soot Properties
		CALL SOOT_PROPERTIES(M,D,fv,SA)
		
		! Increase Due to Aggregation
		D = D*pii(0,1)
		SA = SA*interpolate_pii_term(P,ONE_THIRD)
		
	
	END SUBROUTINE SOOT_PROPERTIES_AGG
	
	
	FUNCTION interpolate_pii_term(P,term)
		! Interpolate pi for specific term from whole order P MOMENTS
		
		REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: P 
		REAL(KIND=DBL), INTENT(IN) :: term
		REAL(KIND=DBL) :: interpolate_pii_term
				
		! Local Variables
		REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1) :: L_interpP, P_x 
		INTEGER :: i 
	
		! P_x term 
		P_x = (/(i, i=0,SIZE(P)-1)/)
		
		! Calculate L_x stencil
		L_interpP = L_fct(REAL(term,DBL),P_x)
		
		! Check if P0 = 0
		IF(P(1)<=0) THEN
			interpolate_pii_term = 0
			RETURN
		END IF 
		
		interpolate_pii_term = EXP(DOT_PRODUCT(L_interpP,LOG(P(:)/P(1))))
	
	END FUNCTION interpolate_pii_term
	
	
	SUBROUTINE interpolate_pii(P,pii)
	! Interpolate Fractional Aggregate Reduced Moments 
	! Input: 
	! Vector of Aggregate Momements WITH THE FOLLOWING INDICES:
	! pii(i,-2)  : r_i - 2*(1/Df-1/3) 
	! pii(i,-1)  : r_i - (1/Df-1/3)
	! pii(i,0): r
	
		IMPLICIT NONE

		REAL(KIND = DBL), INTENT(IN), DIMENSION(:) :: P
		REAL(KIND=DBL), DIMENSION(0:(SIZE(P)-1),-2:2), INTENT(OUT) :: pii
		
		! Local Variables
		INTEGER :: i
		REAL(KIND=DBL) :: S 	!Required moment order shift for pii  
		LOGICAL, SAVE :: FIRST_CALL = .TRUE.
		
		REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: L_interp_P ! Lagrange Interpolating Function
		REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: P_x ! P moments for interpolation
		
		! Calculate Moment Shift Term
		S = 1/D_f-ONE_THIRD
		
		! Initialize by Calculating L(x) stencil and then using the stencil for all other interpolations
		IF (first_call) THEN

            IF (ALLOCATED(L_interp_P)) THEN
                WRITE(*,*) 'L_interp_P Allocated'
            ELSE
                ALLOCATE(L_interp_P(0:SIZE(P)-1,0:SIZE(P)-1,-2:2))
            END IF 
            
			ALLOCATE(P_x(0:SIZE(P)-1))
			
			P_x = (/(i, i=0,SIZE(P)-1)/)
			
            
			! Calculate L(x) Stencil for P interpolation
			DO i = 0,SIZE(P)-1
				!Fractional Moments
				! i - (1/Df-1/3)
				!L_interp(:,2*i) = L_fct(REAL(i,DBL)-S,P_x)
				! i - 2(1/Df-1/3)
				L_interp_P(:,i,-2) = L_fct(REAL(i,DBL)-TWO*S,P_x)
				! i - (1/Df-1/3)
				L_interp_P(:,i,-1) = L_fct(REAL(i,DBL)-S,P_x)
				! i + (1/Df-1/3) 
				!L_interp(:,1+2*i) = L_fct(REAL(i,DBL)+S,P_x)
				L_interp_P(:,i,1) = L_fct(REAL(i,DBL)+S,P_x)
				! i + (2/Df-2/3)
				L_interp_P(:,i,2) = L_fct(REAL(i,DBL)+TWO*S,P_x)
            END DO 
			
			! Set first_call to false
			first_call = .FALSE.		
			
			! Throw away garbage
			DEALLOCATE(P_x)

		END IF
		
		! Initialize pi
		pii(:,:) = 0
		
		! Calculate reduced moments for whole order moments
		! P(1) = P0, P(2) = P1, ...
        
        ! IF P(1) is 0 then do this to prevent divide by 0
        ! THIS CONDITION SHOULD NOT BE REACHED, OTHERWISE SOMETHING ELSE IS WRONG IN THE CODE!!!
		! INCLUDING THIS HERE AS A CHECK
		IF(P(1) == 0) THEN
			WRITE(*,*) 'P(0) = 0, Setting all reduced moments to 0'
			pii(:,:) = 0
			RETURN
        END IF
        
        ! Calculate Reduced Moments
		DO i = 0,(size(P))-1
			! Otherwise calculate reduced moment for whole order moments
			pii(i,0) = P(i+1)/P(1)
        END DO

		! Calculate Reduced Moments
		DO i = 0,(SIZE(P)-1)
            
            ! If P terms are 0 then use regular interpolation, otherwise log interpolation
            ! This 1st condition shouldn't be reached except for very small initial time step
            IF(ANY(P(:)/P(1) == 0)) THEN
                pii(i,-2) = DOT_PRODUCT(L_interp_P(:,i,-2),P(:)/P(1))
			    pii(i,-1) = DOT_PRODUCT(L_interp_P(:,i,-1),P(:)/P(1))
			    pii(i,1) = DOT_PRODUCT(L_interp_P(:,i,1),P(:)/P(1))
				pii(i,2) = DOT_PRODUCT(L_interp_P(:,i,2),P(:)/P(1))
            ELSE
			    pii(i,-2) = EXP(DOT_PRODUCT(L_interp_P(:,i,-2),LOG(P(:)/P(1))))
			    pii(i,-1) = EXP(DOT_PRODUCT(L_interp_P(:,i,-1),LOG(P(:)/P(1))))
			    pii(i,1) = EXP(DOT_PRODUCT(L_interp_P(:,i,1),LOG(P(:)/P(1))))
				pii(i,2) = EXP(DOT_PRODUCT(L_interp_P(:,i,2),LOG(P(:)/P(1))))
			END IF
            
			!pii(3*i+2) = EXP(DOT_PRODUCT(L_interp(:,1+2*i),LOG(P(:)/P(1))))
		END DO 
		
    END SUBROUTINE interpolate_pii
    
	SUBROUTINE COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,H)
	! Calculate Aggregate Coagulation for continuum regime
	
		IMPLICIT NONE
		
		!REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: M, P 
		REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M,P
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(IN) :: mu 
		REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2), INTENT(IN) :: pii
		REAL(KIND=DBL), DIMENSION(1:SIZE(P)-1), INTENT(OUT) :: H
		
		! Local Variables 
		REAL(KIND=DBL) :: temp, M0
		INTEGER :: r,k,r_max
		
		
		! Initialize H 
		H(:) = 0 
		
		!
		M0 = M(1)
		r_max = SIZE(P)-1
		        
		DO r = 2,r_max
			temp = 0
			DO k = 1, r-1 
				temp = temp + binom(r,k)*(2*pii(k,0)*pii(r-k,0) + &
				bm(mu(1),pii(k,1))*bm(mu(-1),pii(r-k,-1)) + &
				bm(mu(-1),pii(k,-1))*bm(mu(1),pii(r-k,1)) + &
				Kcp*(bm(mu(-1),pii(k,-1)))*pii(r-k,0) ) + &
				(pii(k,0)*(bm(mu(-1),pii(r-k,-1))) + &
				bm(mu(1),pii(k,1))*bm(mu(-2),pii(r-k,-2)) + &
				bm(mu(-2),pii(k,-2))*bm(mu(1),pii(r-k,1)))
			END DO 
			
			H(r) = ONE_HALF*Kc*M0*M0*temp
		
		END DO 

	
    END SUBROUTINE COAGULATION_AGGREGATE_C
    
	SUBROUTINE COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,H)
		! Calculate Aggregate Coagulation for Free Molecular regime
	
		IMPLICIT NONE
		
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M,P
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(IN) :: mu 
		REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2), INTENT(IN) :: pii
		REAL(KIND=DBL), DIMENSION(1:SIZE(P)-1), INTENT(OUT) :: H
		
		! Local Variables 
		REAL(KIND=DBL), DIMENSION(4) :: f
		INTEGER :: r, l,n, rr
		REAL(KIND=DBL) :: M02
		
		
		
		! Initialize H 
		H(:) = 0
		H(1) = 0_DBL 
		
		M02 = M(0)*M(0)
		
		DO r = 2,SIZE(P)-1 
			! Calculate Grid Function 
			DO l = 1, 4
				f(l) = gridFun(l-1,r,mu,mu_LO,mu_HI,pii)
			END DO
			
			! Interpolate Grid Function and Calculate H(r)
			H(r) = 0.5_DBL*kf*M02*L_LOG_INTERP(0.5D0,(/0D0, 1D0, 2D0, 3D0/),f)
		
		END DO

		CONTAINS
		
		FUNCTION gridFun(l,r,mu,mu_LO,mu_HI,pii)
			! Calculate Aggregate Grid Function for FM calculation
			
			IMPLICIT NONE
			INTEGER, INTENT(IN) :: l, r
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(IN) :: mu 	
			REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2), INTENT(IN) :: pii
			REAL(KIND=DBL) :: gridFun
			
			! Local Variables
			INTEGER :: k, q
			
			REAL(KIND=DBL) :: sum1
		
			! Initialize gridFun Output and sum1
			gridFun = 0_DBL
			sum1 = 0_DBL
			
			DO k = 0,l
				
				DO q = 1,r-1
					sum1 = sum1 + binom(r,q) * ( &
					bm(mu(6*k+1),pii(q,2))*bm(mu(6*l-6*k-3),pii(r-q,0)) &
					+ 2*bm(mu(6*k-1),pii(q,1))*bm(mu(6*l-6*k-1),pii(r-q,1)) &
					+ bm(mu(6*k-3),pii(q,0))*bm(mu(6*l-6*k+1),pii(r-q,2)))	
				END DO
				
				gridFun = gridFun + binom(l,k)*sum1
			
			END DO 			
			
		END FUNCTION gridFun
	
	END SUBROUTINE COAGULATION_AGGREGATE_FM
	
	SUBROUTINE COAGULATION_AGGREGATE_TR(M,P,mu,mu_LO,mu_HI,pii,H)
		! Calculate Transitional Regime Rates for Aggregate Coagulation
		
		IMPLICIT NONE
		REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M,P
		INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
		REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(IN) :: mu 
		REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2), INTENT(IN) :: pii
		REAL(KIND=DBL), DIMENSION(1:SIZE(P)-1), INTENT(OUT) :: H
	
		! Local Variables
		REAL(KIND=DBL), DIMENSION(1:SIZE(P)-1) :: Hc, Hf
		
		! Initialize H
		H(:) = 0_DBL 
		
		! Calculate Free Molecular Rates
		CALL COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,Hc)
		
		! Calculate Continuum Coagulation Rates
		CALL COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,Hf)
		
		! Calculate Transitional Rates as Harmonic Mean
		H(2:) = (Hf(2:)*Hc(2:))/(Hf(2:)+Hc(2:))
	
	END SUBROUTINE COAGULATION_AGGREGATE_TR
	
	
	!! JACOBIAN CALCULATIONS
	
    SUBROUTINE Jacobian(x,num_M,num_P,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,J)
    ! Calculate Numerical Jacobian using INTEL MKL
    IMPLICIT NONE
    
    ! Inputs
    REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: x
    REAL(KIND=DBL), INTENT(IN) :: TEMPERATURE, PRESSURE, C_C2H2, C_H, C_H2, C_H2O, C_O2, C_OH
	INTEGER, INTENT(IN) :: num_M, num_P
        
    REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(x),SIZE(X)) :: J
    ! Cols: M1, M2, M3,.. T, C2H2, H, H2, H2O, O2, OH
    ! Rows: dSM1,dSM2,...,dSP1,dSP2,..., dC2H2, dCO, dH, dH2, dH2O, dO2, dOH   
 
    ! Local Variable
    INTEGER :: res, m, n
    REAL(KIND=DBL) :: eps = 1e-5    ! Pertubation
          
    ! NEED TO DEFINE THESE FOR Jacobian Subroutine
    external djacobi
    integer djacobi
    integer TR_SUCCESS
    parameter (TR_SUCCESS = 1501)
    integer TR_INVALID_OPTION
    parameter (TR_INVALID_OPTION = 1502)
    integer TR_OUT_OF_MEMORY
    parameter (TR_OUT_OF_MEMORY = 1503)
    
    ! Initialize Jacobian
    J(:,:) = 0.D0
    
    ! Calculate Dimensions
    m = SIZE(x)
    n = SIZE(x)
     
    ! Calculate Jacobian
    res = djacobi(J_eval,n,m,J,x,eps)
 
    
    IF (res /= TR_SUCCESS) THEN
        WRITE(*,*) "THERE'S SOME ISSUE WITH JACOBIAN CALCULATION"
        WRITE(*,*) "AT x:"
        WRITE(*,*) x
    END IF 
   
    
    CONTAINS    
        SUBROUTINE J_eval(m,n,x,f)
        ! Function for evaluating Jacobian, (CONSTANT T,P, Concentrations)    
        
        INTEGER, INTENT(IN) :: m,n
        REAL(KIND=DBL), INTENT(IN), DIMENSION(n) :: x
        REAL(KIND=DBL), INTENT(OUT), DIMENSION(m) :: f
        
        ! Local Variables	
        REAL(KIND=DBL), DIMENSION(num_P) :: H, Ragg
        REAL(KIND=DBL), DIMENSION(num_M) :: W, G, R
        REAL(KIND=DBL), DIMENSION(num_M) :: Mx
		REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: Px
		
        ! Don't need these for anything 
        REAL(KIND=DBL) :: rC2H2, rCO, rH, rH2, rH2O, rO2, rOH
        
        ! Initialize f
        f(:) = 0.D0
       
	    ! Split input into M and P Moments
		Mx = x(1:num_M)
		
		IF (num_P>0) THEN
			ALLOCATE(Px(num_P))
			Px = x(num_M+1:)
		ELSE
			IF (AGGREGATION) THEN
				WRITE(*,*) 'ERROR IN JACOBIAN CALCULATION, AGGREGATION IS TURNED IN INPUT &
				FILE BUT NO P MOMENTS ARE SPECIFIED IN INPUT TO JACOBIAN SUBROUTINE'
				STOP
			ELSE
				! Arbitrary P values (Not used in calculation)
				ALLOCATE(Px(1))
				Px = 0.D0
				
			END IF
		END IF 
				
        !CALL CALCULATE_SOURCE(xM,xP,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
        !W,G,R,H,Ragg,rC2H2, rCO, rH, rH2, rH2O, rO2, rOH)
		
		
		
        CALL CALCULATE_SOURCE(Mx,Px,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
        W,G,R,H,Ragg,rC2H2, rCO, rH, rH2, rH2O, rO2, rOH)
		
		!f = W+G+R
		
        f(1:num_M) = W + G + R
		
		IF (num_P > 0) THEN
			f(num_M+1:num_P) = H + Ragg
		END IF 
		
		! Throw Away
		DEALLOCATE(Px)
                 
        END SUBROUTINE J_eval
            
    END SUBROUTINE Jacobian

    SUBROUTINE Jacobian_log(x,num_M,num_P,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,J)
    ! Calculate Numerical Jacobian using INTEL MKL
    IMPLICIT NONE
    
    ! Inputs
    REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: x
    REAL(KIND=DBL), INTENT(IN) :: TEMPERATURE, PRESSURE, C_C2H2, C_H, C_H2, C_H2O, C_O2, C_OH
	INTEGER, INTENT(IN) :: num_M, num_P
    
    REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(x),SIZE(X)) :: J
    ! Cols: M1, M2, M3,.. T, C2H2, H, H2, H2O, O2, OH
    ! Rows: dSM1,dSM2,...,dSP1,dSP2,..., dC2H2, dCO, dH, dH2, dH2O, dO2, dOH   
 
    ! Local Variable
    INTEGER :: res, m, n, i
    REAL(KIND=DBL) :: eps = 1e-5    ! Perturbation 
		  
    ! NEED TO DEFINE THESE FOR Jacobian Subroutine
    external djacobi
    integer djacobi
    integer TR_SUCCESS
    parameter (TR_SUCCESS = 1501)
    integer TR_INVALID_OPTION
    parameter (TR_INVALID_OPTION = 1502)
    integer TR_OUT_OF_MEMORY
    parameter (TR_OUT_OF_MEMORY = 1503)
    
    ! Initialize Jacobian
    J(:,:) = 0.D0
    
    ! Calculate Dimensions
    m = SIZE(x)
    n = SIZE(x)
    	

    res = djacobi(J_eval,n,m,J,x,eps)
 
    
    IF (res /= TR_SUCCESS) THEN
        WRITE(*,*) "THERE'S SOME ISSUE WITH JACOBIAN CALCULATION"
        WRITE(*,*) "AT x:"
        WRITE(*,*) x
    END IF 
   
    
    CONTAINS    
        SUBROUTINE J_eval(m,n,x,f)
        ! Function for evaluating Jacobian, (CONSTANT T,P, Concentrations)    
        
        INTEGER, INTENT(IN) :: m,n
        REAL(KIND=DBL), INTENT(IN), DIMENSION(n) :: x
        REAL(KIND=DBL), INTENT(OUT), DIMENSION(m) :: f
        
        ! Local Variables	
        REAL(KIND=DBL), DIMENSION(num_P) :: H, Ragg
        REAL(KIND=DBL), DIMENSION(num_M) :: W, G, R
        REAL(KIND=DBL), DIMENSION(num_M) :: Mx
		REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: Px
		
        ! Don't need these for anything 
        REAL(KIND=DBL) :: rC2H2, rCO, rH, rH2, rH2O, rO2, rOH
        
        ! Initialize f
        f(:) = 0.D0
       
	    ! Split input into M and P Moments
		Mx = x(1:num_M)
		
		IF(ANY(Mx < Mx(1))) THEN
			WRITE(*,*) 'Higher Order Moments < M0 for input in Jacobian Calculation'
			WRITE(*,*) 'Setting Higher Order Moment to ~M0 to prevent program stopping'
			WRITE(*,*) 'MKL Jacobian Calculation Affected, Calculated Jacobian is probably not accurate'
			
			DO i = 2,num_M
				IF (Mx(i)<Mx(i-1)) THEN
					Mx(i) = Mx(i-1)+1D0
				END IF 	
			END DO 

		END IF 
		
		
		IF (num_P>0) THEN
			ALLOCATE(Px(num_P))
			Px = x(num_M+1:)
		ELSE
			IF (AGGREGATION) THEN
				WRITE(*,*) 'ERROR IN JACOBIAN CALCULATION, AGGREGATION IS TURNED ON IN INPUT &
				FILE BUT NO P MOMENTS ARE SPECIFIED IN INPUT TO JACOBIAN SUBROUTINE'
				STOP
			ELSE
				! Arbitrary P values (Not used in calculation)
				ALLOCATE(Px(1))
				Px = 0.D0
				
			END IF
		END IF 
				
        !CALL CALCULATE_SOURCE(xM,xP,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
        !W,G,R,H,Ragg,rC2H2, rCO, rH, rH2, rH2O, rO2, rOH)
		
		
		
        CALL CALCULATE_SOURCE(Mx,Px,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
        W,G,R,H,Ragg,rC2H2, rCO, rH, rH2, rH2O, rO2, rOH)
		
		!f = W+G+R
		
        f(1:num_M) = (W + G + R)/Mx
        
        !f(1) = f(1)/Mx(1)
        !f(2) = f(2)/Mx(2)
        !f(3) = f(3)/Mx(3)
        !f(4) = f(4)/Mx(4)
        !f(5) = f(5)/Mx(5)
        !f(6) = f(6)/Mx(6)
		
		IF (num_P > 0) THEN
			f(num_M+1:num_P) = (H + Ragg)/Px
		END IF 
		
		! Throw Away
		DEALLOCATE(Px)
                 
        END SUBROUTINE J_eval
            
    END SUBROUTINE Jacobian_log
    
    
    SUBROUTINE Calc_Jacobian_no_aggregate(M,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,J)
    ! Calculate Numerical Jacobian using INTEL MKL
    
    IMPLICIT NONE
    ! Inputs
    REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M 
    REAL(KIND=DBL), INTENT(IN) :: TEMPERATURE, PRESSURE, C_C2H2, C_H, C_H2, C_H2O, C_O2, C_OH
    
    ! Local Variable
   ! REAL(KIND=DBL), DIMENSION(SIZE(M)) :: Jeval
    
    ! Jacobian Output
    REAL(KIND=DBL), DIMENSION(SIZE(M),SIZE(M)), INTENT(OUT) :: J
    
    REAL(KIND=DBL), DIMENSION(SIZE(M)) :: f
    
    INTEGER :: res
    REAL(KIND=DBL) :: eps = 1e-5    ! Pertubation
    
    ! Jacobian Subroutine
    external djacobi
    integer djacobi
    integer TR_SUCCESS
    parameter (TR_SUCCESS = 1501)
    integer TR_INVALID_OPTION
    parameter (TR_INVALID_OPTION = 1502)
    integer TR_OUT_OF_MEMORY
    parameter (TR_OUT_OF_MEMORY = 1503)
    
    !res = djacobi(J_eval_no_aggregate,SIZE(M),SIZE(M),J,M,eps)
    res = djacobi(J_eval_no_aggregate,SIZE(M),SIZE(M),J,M,eps)
    WRITE(*,*) SIZE(M),SIZE(M),M
    WRITE(*,*) J(1,:)
    WRITE(*,*) J(2,:)
    WRITE(*,*) J(3,:)

	    
    if (res /= TR_SUCCESS) THEN
        WRITE(*,*) "THERE'S SOME ISSUE WITH JACOBIAN CALCULATION"
        WRITE(*,*) "AT M:"
        WRITE(*,*) M
	END IF 
    
    
    CONTAINS    
        SUBROUTINE J_eval_no_aggregate(m,n,x,f)
        ! Function for evaluating Jacobian (NO aggregate)    
        
        INTEGER, INTENT(IN) :: m,n
        REAL(KIND=DBL), INTENT(IN), DIMENSION(n) :: x
        REAL(KIND=DBL), INTENT(OUT), DIMENSION(m) :: f
        
        ! Local Variables
        REAL(KIND=DBL) :: rC2H2, rCO, rH, rH2, rH2O, rO2, rOH
        REAL(KIND=DBL), DIMENSION(0:1) :: P, H, Ragg
        REAL(KIND=DBL), DIMENSION(m) :: W, G, R
        
        CALL CALCULATE_SOURCE(x,P,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
        W,G,R,H,Ragg,rC2H2, rCO, rH, rH2, rH2O, rO2, rOH)
        
        f = W + G + R 
        
        END SUBROUTINE J_eval_no_aggregate
        
    
    
    END SUBROUTINE Calc_Jacobian_no_aggregate
    
    ! DO THIS
    !FUNCTION Calc_Jacobian_aggregate(M)
    !END FUNCTION Calc_Jacobian_aggregate
    
	
   SUBROUTINE CALCULATE_SOURCE2(M,P_in,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,AGG, &
        W,G,R,H,Ragg,rC2H2, rCO, rH, rH2, rH2O, rO2, rOH)
        ! Calculate Moment Source Terms 
        IMPLICIT NONE
    
        REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M, P_in
        REAL(KIND=DBL), INTENT(IN) :: TEMPERATURE, PRESSURE
		REAL(KIND=DBL), INTENT(IN) :: C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH
		INTEGER, OPTIONAL, INTENT(IN) :: AGG
        REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(M)) :: W, G, R
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(P_in)) :: H, Ragg
        REAL(KIND=DBL), INTENT(OUT) :: rC2H2, rCO, rH, rH2, rH2O, rO2, rOH
		REAL(KIND=DBL) :: rC2H2_nuc, rC2H2_surf
		
        ! Local Variables
        REAL(KIND=DBL), ALLOCATABLE :: mu(:)
        REAL(KIND=DBL), DIMENSION(SIZE(P_in)+1) :: P
		!REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: P 
        INTEGER :: i                         
		LOGICAL :: error_flag = .FALSE. ! Error Flag
		REAL(KIND=DBL) :: kn ! Knudsen Number
		INTEGER :: COAGULATION_REGIME  ! Free-Molecular/Transitional/Continuum Regime (Determined by Knudsen #)
		
		REAL(KIND=DBL), DIMENSION(SIZE(M)) :: f
		REAL(KIND=DBL) :: a,b,rsq 
		
        REAL(KIND=DBL), DIMENSION(SIZE(M)) :: W_C2H2, W_O2, W_OH
                 
        ! Initialize Moments Module:
        ! Calculate Binomial Coefficient
        ! Calculate Lagrange Interpolation Stencil
        ! Calculate Number of P and R moments
        ! Calculate lowest fraction and highest fraction moments
        IF(.NOT. isInitialized) THEN
            ! Read Moments Input
            CALL INITIALIZE_MOMIC()
		
			! Check if correct initial conditions are being used
			!IF(ALL(M==M_0)) THEN
			!	CONTINUE
			!ELSE
			!	WRITE(*,*) 'PROBLEM READING INITIAL CONDITIONS, THE M Moment at initial call to CALCULATE_SOURCE &
			!	is different from what is specified in input IC'
			!	WRITE(*,*) 'M at initial Call: ', M
			!	WRITE(*,*) 'M in IC file: ', M_0 
			!	error_flag = .TRUE.
            !END IF
				
            ! Check P moments if it's allocated
            !IF (CALC_AGGREGATE) THEN
			!    IF(ALL(P==P_0)) THEN 
			!	    CONTINUE
			 !   ELSE
			!	    WRITE(*,*) 'PROBLEM READING INITIAL CONDITIONS, THE P Moment at initial call to CALCULATE_SOURCE &
			!	    is different from what is specified in input IC'
			!	    WRITE(*,*) 'P at initial Call: ', P
			!	    WRITE(*,*) 'P in IC file: ', P_0  
			!	    error_flag = .TRUE.
             !   END IF
            !END IF 
			
			IF(error_flag) THEN
				WRITE(*,*) 'QUITING!!!'
				STOP
			END IF
		
        END IF
    
	
        ! ************ MOMENT ERROR CHECKING ************
		! Check if P(0) = M(0) within some error tolerance
		!IF(ABS(P(1) - M(1))>1e-8 .AND. AGGREGATION_CALC) THEN 
		!	WRITE(*,*) 'ERROR: P(0):',P(1), ' /= M(0):', M(1)
		!	STOP
        !END IF
	    ! Don't need this, just set P(0) to M(0)
		
		
		!ALLOCATE(P(1:P_Moments))
		
		If (ANY(M(1) > P_in)) THEN
			P(:) = M(1)!1.D0
		ELSE
			! P(0) = M(0)
			P(1) = M(1)
			! Input Moments
			P(2:) = P_in
        END IF 
		
		!IF(present(AGG)) THEN
		!	IF(AGG==1) THEN
		!		AGGREGATION = .FALSE.
		!	ELSEIF(AGG==2) THEN
		!		AGGREGATION = .TRUE.
		!	END IF 
		!END IF
		
		
        ! Initialize Nucleation, Coagulation, Surface Rates, And Aggregate Rates
	    W(:) = 0
	    G(:) = 0
	    R(:) = 0
        H(:) = 0
        Ragg(:) = 0
        rC2H2_nuc = 0
		rC2H2_surf = 0
		rC2H2 = 0 
        rCO = 0
        rH = 0
        rH2 = 0
        rH2O = 0
        rO2 = 0
        rOH = 0
		OXID = .TRUE.
		
        ! Check if MOMIC Calculation is enabled
        IF (.NOT. doSOOT) THEN
            WRITE(*,*) 'MOMENT CALCULATION TURNED OFF IN MOMIC INPUT FILE'
            RETURN
        END IF
		
		! Check Inputs are correct
		! THIS CAN CAUSE POTENTIAL MEMORY ISSUES IF THEY ARE NOT THE SAME
		! MAYBE LOOK AT EDITING HOW MOMENTS ARE INPUT?
		IF (SIZE(M) /= M_MOMENTS) THEN
			WRITE(*,*) 'ERROR: THE NUMBER OF INPUT M MOMENTS DOES NOT MATCH UP &
			WITH THE NUMBER OF MOMENTS SPECIFIED IN THE INPUT FILE'
			STOP
		END IF
		
		IF (AGGREGATION .AND. SIZE(P) /= P_MOMENTS) THEN
			WRITE(*,*) 'ERROR: THE NUMBER OF INPUT P MOMENTS DOES NOT MATCH UP &
			WITH THE NUMBER OF MOMENTS SPECIFIED IN THE INPUT FILE'
			STOP
		END IF
		    
        ! Calculate Fractional MOMENTS
		ALLOCATE(mu(mu_LO:mu_HI))
		
        ! Interpolate reduced moments
		CALL interpolate_mu(M,mu,mu_LO,mu_HI)		
		
		
		! Oxidation Calculation Control Flag
		IF (mu(6) <= OXID_RAT) THEN
			OXID = .FALSE.
		END IF
		
		
        ! Calculate Kc, Kc' and Kf Coagulation Parameters
		CALL COAGULATION_PARAMETERS(TEMPERATURE,PRESSURE)
		        	
		
		! Calculate Average soot diameter 
		D_soot = Soot_D_Calc(mu,mu_LO,mu_HI)
		
		! Calculate Knudsen Number
		kn = 2*MFP(PRESSURE,TEMPERATURE)*(1e-2)/D_soot
								
		! Kn number regimes
		IF (kn < 0.1) THEN
			COAGULATION_REGIME = 0 ! Continuum
		ELSEIF (kn < 10) THEN
			COAGULATION_REGIME = 1 ! Transitional
		ELSEIF (kn >= 10) THEN
			COAGULATION_REGIME = 2 ! Free Molecular
		ELSEIF (isNAN(kn)) THEN
			WRITE(*,*) "ERROR IN Temperature or Pressure input"
		ELSE
			WRITE(*,*) "HOW DID YOU GET TO THIS CONDITION??? THIS SHOULD BE PHYISCALLY IMPOSSIBLE!!"
			STOP
		END IF 
		  		  
		IF (FORCE_REGIME) THEN
			COAGULATION_REGIME = REGIME
		END IF         
		
        ! Calculate Primary Particle Nucleation Rates
		CALL NUCLEATION(M,TEMPERATURE,PRESSURE,C_C2H2,rC2H2_nuc,R,M_Moments)
		
		! If Number density is less than 1/cm^3 then only use Nucleation calculation,  skip the rest of calculation
		! Like in Ranjan's Code, this helps with some cases!!
		IF(M(1)<150D0) THEN
            WRITE(*,*) 'Soot Number Density M(0) < 150'
            WRITE(*,*) 'Only Calculating Nucleation Source'
			RETURN
		END IF 
					
		! Check if Aggregation is turned on
		
		!AGGREGATION = .FALSE.
		
        If(AGGREGATION) THEN
            CALL CALCULATE_SOURCE_AGG
        ELSE
            CALL CALCULATE_SOURCE_NO_AGG
        END IF
		
		! Get Total Rate of Nucleation and Surface Growth C2H2
		rC2H2 = rC2H2_nuc+rC2H2_surf
		
		
        ! ********************************************************
        ! ************* ERROR CHECKING ***************************
        ! ********************************************************
		! Check if rates are nan
        DO i = 1,M_Moments
			IF( isNAN(G(i)) .OR. isNAN(W(i)) .OR. isNAN(R(i))) THEN
				WRITE(*,*) 'NAN for Output Rates'
				!WRITE(*,*) 'R -- Gr -- Wr -- Rr'
				!WRITE(*,*) i, G(i), W(i), R(i)
				error_flag = .TRUE.
				!EXIT
            END IF
        END DO 
        
        ! Check if M(n) < M(0)
        !IF(ANY(M(2:)<M(1))) THEN
         !   WRITE(*,*) 'Error higher moments are smaller than 0th moment'
          !  error_flag = .True.
        !END IF
        
		CALL VERIFY_DISTRIBUTION(M,error_flag)
		
	
        ! DEBUGGING STUFF
        numCalls = numCalls + 1
		

		! Check Error Flag and Stop if error
		IF(error_flag) THEN
            WRITE(*,*) 'NumCalls = ',NumCalls
			WRITE(*,*) '*******INPUTS********'
			WRITE(*,*) 'M Moments'
			WRITE(*,"(6ES14.3)") M
			WRITE(*,*) 'P Moments'
			WRITE(*,"(6ES14.3)") P
			WRITE(*,*) 'Temperature - Pressure '
			WRITE(*,*) TEMPERATURE, PRESSURE
			WRITE(*,*) '[C2H2] - [H] - [H2]'
			WRITE(*,*) C_C2H2, C_H, C_H2
			WRITE(*,*) '[H2O] - [O2] - [OH]'
			WRITE(*,*) C_H2O, C_O2, C_OH
			WRITE(*,*) '******Calculations******'
			WRITE(*,*) 'mu'
			WRITE(*,*) mu
			WRITE(*,*) '---- Gr Rates --- '
			WRITE(*,*) G
			WRITE(*,*) '---- Wr Rates --- '
			WRITE(*,*) W
			WRITE(*,*) '---- Rr Rates --- '
			WRITE(*,*) R			
        END IF
		
		! f2py is not compiling if STOP is in the above if statement for some reason
		!IF(error_flag) THEN
		!	STOP
		!END IF 
		
        ! Normalize all Equations by M_ref
        R = R/M_ref
        G = G/M_ref
        W = W/M_ref
        Ragg = Ragg/M_ref
        H = H/M_ref

        ! Don't need this anymore so throw it away    		 
		DEALLOCATE(mu)
		!DEALLOCATE(P)
		
        RETURN
		
    CONTAINS
		SUBROUTINE VERIFY_DISTRIBUTION(M,error_flag)
			! Check Realizability of Distribution Function
			
			IMPLICIT NONE
			REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
			LOGICAL, INTENT(OUT) :: error_flag
			
			
			! Initialize error_flag
			error_flag = .FALSE.
			
			IF(M(0)<1) THEN
				WRITE(*,*) 'ERROR IN DISTRIBUTION: M(0) < 1'
				error_flag = .TRUE.
			ELSEIF(M(0)*M(2)/(M(1)*M(1))<1) THEN
				WRITE(*,*) 'ERROR in Distribution: M(0)M(2)/M(1)^2 < 1'
				error_flag = .TRUE.
			END IF
		END SUBROUTINE
	
		SUBROUTINE PRIMARY_COAGULATION(M,mu,mu_LO,mu_HI,G)
		
			IMPLICIT NONE
			REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
			REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:SIZE(M)-1) :: G
		
		
			IF(COAGULATION_REGIME==0) THEN
				CALL PRIMARY_COAGULATION_C(M,mu,mu_LO,mu_HI,G)
				WRITE(*,*) 'Coagulation is being calculated for Continuum Regime'
			ELSEIF(COAGULATION_REGIME==1) THEN
				CALL PRIMARY_COAGULATION_TR(M,mu,mu_LO,mu_HI,G)
				WRITE(*,*) 'Coagulation is being calculated for Transitional Regime'
			ELSEIF(COAGULATION_REGIME==2) THEN
				CALL PRIMARY_COAGULATION_FM(M,mu,mu_LO,mu_HI,G)
				WRITE(*,*) 'Coagulation is being calculated for Free Molecular Regime'
			ELSE
				WRITE(*,*) "WRONG COAGULATION REGIME FLAG, I DID SOMETHTHING WRONG?"
				STOP
			END IF 
		
		END SUBROUTINE PRIMARY_COAGULATION
		
		SUBROUTINE PRIMARY_COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,G)
			! Calculate Primary Coagulation with aggregation depending on COAGULATION_REGIME Flag
			IMPLICIT NONE
			REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M, P 
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
			REAL(KIND=DBL), INTENT(IN), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
			REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(M)) :: G
			
			IF(COAGULATION_REGIME==0) THEN
				CALL PRIMARY_COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,G)
				WRITE(*,*) 'Coagulation is being calculated for Continuum Regime'
			ELSEIF(COAGULATION_REGIME==1) THEN
				CALL PRIMARY_COAGULATION_AGGREGATE_TR(M,P,mu,mu_LO,mu_HI,pii,G)
				WRITE(*,*) 'Coagulation is being calculated for Transitional Regime'
			ELSEIF(COAGULATION_REGIME==2) THEN
				CALL PRIMARY_COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,G)
				WRITE(*,*) 'Coagulation is being calculated for Free Molecular Regime'
			ELSE
				WRITE(*,*) "WRONG COAGULATION REGIME FLAG, I DID SOMETHTHING WRONG"
				STOP
			END IF 	
		
		END SUBROUTINE PRIMARY_COAGULATION_AGGREGATE
		
		SUBROUTINE COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,H)
			! Calculate Aggregate Coagulation depending on COAGULATION_REGIME Flag
			
			IMPLICIT NONE
			
			!REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: M, P 
			REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M,P
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(IN) :: mu 
			REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2), INTENT(IN) :: pii
			REAL(KIND=DBL), DIMENSION(1:SIZE(P)-1), INTENT(OUT) :: H
			
			IF(COAGULATION_REGIME==0) THEN
				CALL COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,H)
				WRITE(*,*) 'Coagulation is being calculated for Continuum Regime'
			ELSEIF(COAGULATION_REGIME==1) THEN
				CALL COAGULATION_AGGREGATE_TR(M,P,mu,mu_LO,mu_HI,pii,H)
				WRITE(*,*) 'Coagulation is being calculated for Transitional Regime'
			ELSEIF(COAGULATION_REGIME==2) THEN
				CALL COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,H)
				WRITE(*,*) 'Coagulation is being calculated for Free Molecular Regime'
			ELSE
				WRITE(*,*) "WRONG COAGULATION REGIME FLAG, I DID SOMETHTHING WRONG"
				STOP
			END IF 		
		
		END SUBROUTINE COAGULATION_AGGREGATE
	
        SUBROUTINE CALCULATE_SOURCE_NO_AGG()
                ! Calculate Source Terms without Aggregation
                IMPLICIT NONE
        
                WRITE(*,*) 'Calculating Moment Rates with out Aggregation'
                
				
				CALL PRIMARY_COAGULATION(M,mu,mu_LO,mu_HI,G)
				
				
				CALL SURFACE_GROWTH(M,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
									rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
				
				
				! If Oxidation Regime then base surface rates on M1, like in Ranjan's code

				IF(W(2) < 0) THEN
					
					WRITE(*,*) '*** IN OXIDATION REGIME, Basing Moments on M(1) ***'
					DO i = 1,SIZE(M)-1
						f(i) = log(mu(6*i))
					END DO
					
					
					CALL linear(SIZE(M)-1,f,a,b,rsq)
					
						
					DO i = 4,6*(SIZE(M)-1)-1,6
						mu(i) = EXP(a+b*i/6.D0)
					END DO 
					
					
					!WRITE(*,*) '*******ORIGINAL W********'
					!WRITE(*,*) W 
					
					
					CALL SURFACE_GROWTH(M,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
				rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
						
				
					DO i =2,SIZE(M)
						W(i) = W(i)*M(i)/M(1)/EXP(a+b*i)
					END DO
					

				
					!WRITE(*,*) '*****New W********'
					!WRITE(*,*) W
					!STOP
				
				END IF 
                
				
				! Turn off Source Terms Cutoff Criteria Like in S.P Roy, Haworth Paper
				! For Strongly Oxidizing Environment
				!IF (W(2)<0 .AND. M(2)<32*M(1)) THEN
				!	W(:) = 0
				!	RETURN
				!END IF 
				
				
            END SUBROUTINE CALCULATE_SOURCE_NO_AGG
            
            SUBROUTINE CALCULATE_SOURCE_AGG()
	            ! CALCULATE SOURCE TERMS FOR METHOD OF MOMENTS USING KAZAKOV AGGREGATION MODEL
	            ! ASSUMING Average Particle Diameter Dsoot > D*
                IMPLICIT NONE
            
                ! Local Variable for reduced P moment
                REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
                
                ! Initialize Variables
                pii(:,:) = 0
                
                WRITE(*,*) 'Calculating Moment Rates with Aggregation Enabled'
		
				IF(present(AGG)) THEN
					IF (AGG==2) THEN
						D_soot =2*d_star+1
						WRITE(*,*) 'FORCING AGGREGATION ON'
					ELSEIF(AGG==1) THEN
						D_soot = d_star/2-1
						WRITE(*,*) 'FORCING AGGREGATION OFF'
					ENDIF
				ENDIF
		
                IF(D_soot > d_star) THEN

	                ! Interpolate pi 
	                CALL interpolate_pii(P,pii)
                    
					
                    ! Calculate Aggregate Nucleation
                    Ragg(:) = R(1)
                    
                    ! If P terms are 0 then only use nucleation term and non-aggregate coagulation, skip rest of calculation
                    IF(ANY(P == 0)) THEN
                        CALL PRIMARY_COAGULATION(M,mu,mu_LO,mu_HI,G)
                        RETURN
                    END IF
                    
            
	                ! Calculate Coagulation
	                CALL PRIMARY_COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,G)

					! Calculate Aggregate Coagulation
					CALL COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,H)
			
					! Calculate Surface Growth 
					CALL SURFACE_GROWTH_AGG(M,P,pii,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
											rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
		
		
					
					! If Oxidation Regime then base surface rates on M1, like in Ranjan's code
					IF(W(2)<0) THEN
						DO i = 1,SIZE(M)-1
							f(i) = log(mu(3*i))
						END DO	
						
						CALL linear(SIZE(M)-1,f,a,b,rsq)
							
						DO i = 2,3*(SIZE(M)-1)-1,3
							mu(i) = EXP(a+b*i/3.D0)
						END DO 
						
						
						!WRITE(*,*) '*******ORIGINAL W********'
						!WRITE(*,*) W 
						
						CALL SURFACE_GROWTH_AGG(M,P,pii,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
											rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
					
					
						DO i =2,SIZE(M)
							W(i) = W(i)*M(i)/M(1)/EXP(a+b*i)
						END DO
						
						! Oxidation should also result in reduced number density as particles get oxidized
						! Original Formulation does not account for that, so using this to account for it
						! W0 = W_OXIDATION/Cmin
						! Alternatively maybe set M0 = M1/Cmin?
						
						W(0) = (W_O2(1)+W_OH(1))/10D0
						
						
						!WRITE(*,*) '*****New W********'
						!WRITE(*,*) W
						!STOP
					
					END IF 							
                                        
                ELSE
                    WRITE(*,*) 'Average Soot D < d*'
                    CALL CALCULATE_SOURCE_NO_AGG()
                END IF
        
				RETURN
            END SUBROUTINE CALCULATE_SOURCE_AGG
                        
	END SUBROUTINE CALCULATE_SOURCE2
	
	
	SUBROUTINE CALCULATE_SOURCE_SWITCHER(M,P_in,TEMPERATURE,PRESSURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH, COAGULATION_REGIME, AGGREGATE, &
        W,G,R,H,Ragg,rC2H2, rCO, rH, rH2, rH2O, rO2, rOH)
        ! Calculate Moment Source Terms 
        IMPLICIT NONE
    
        REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M, P_in
        REAL(KIND=DBL), INTENT(IN) :: TEMPERATURE, PRESSURE
		REAL(KIND=DBL), INTENT(IN) :: C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH
		INTEGER, INTENT(IN) :: COAGULATION_REGIME
		LOGICAL, INTENT(IN) :: AGGREGATE
        REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(M)) :: W, G, R
		REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(P_in)) :: H, Ragg
        REAL(KIND=DBL), INTENT(OUT) :: rC2H2, rCO, rH, rH2, rH2O, rO2, rOH
		REAL(KIND=DBL) :: rC2H2_nuc, rC2H2_surf
		

		
		
        ! Local Variables
        REAL(KIND=DBL), ALLOCATABLE :: mu(:)
        REAL(KIND=DBL), DIMENSION(SIZE(P_in)+1) :: P
		!REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: P 
        INTEGER :: i                         
		LOGICAL :: error_flag = .FALSE. ! Error Flag
		REAL(KIND=DBL) :: kn ! Knudsen Number
		!INTEGER :: COAGULATION_REGIME  ! Free-Molecular/Transitional/Continuum Regime (Determined by Knudsen #)
		
		REAL(KIND=DBL), DIMENSION(SIZE(M)) :: f
		REAL(KIND=DBL) :: a,b,rsq 
		
        REAL(KIND=DBL), DIMENSION(SIZE(M)) :: W_C2H2, W_O2, W_OH
                 
        ! Initialize Moments Module:
        ! Calculate Binomial Coefficient
        ! Calculate Lagrange Interpolation Stencil
        ! Calculate Number of P and R moments
        ! Calculate lowest fraction and highest fraction moments
        IF(.NOT. isInitialized) THEN
            ! Read Moments Input
            CALL INITIALIZE_MOMIC()
		
			! Check if correct initial conditions are being used
			!IF(ALL(M==M_0)) THEN
			!	CONTINUE
			!ELSE
			!	WRITE(*,*) 'PROBLEM READING INITIAL CONDITIONS, THE M Moment at initial call to CALCULATE_SOURCE &
			!	is different from what is specified in input IC'
			!	WRITE(*,*) 'M at initial Call: ', M
			!	WRITE(*,*) 'M in IC file: ', M_0 
			!	error_flag = .TRUE.
            !END IF
				
            ! Check P moments if it's allocated
            !IF (CALC_AGGREGATE) THEN
			!    IF(ALL(P==P_0)) THEN 
			!	    CONTINUE
			 !   ELSE
			!	    WRITE(*,*) 'PROBLEM READING INITIAL CONDITIONS, THE P Moment at initial call to CALCULATE_SOURCE &
			!	    is different from what is specified in input IC'
			!	    WRITE(*,*) 'P at initial Call: ', P
			!	    WRITE(*,*) 'P in IC file: ', P_0  
			!	    error_flag = .TRUE.
             !   END IF
            !END IF 
			
			IF(error_flag) THEN
				WRITE(*,*) 'QUITING!!!'
				STOP
			END IF
		
        END IF
    
	
        ! ************ MOMENT ERROR CHECKING ************
		! Check if P(0) = M(0) within some error tolerance
		!IF(ABS(P(1) - M(1))>1e-8 .AND. AGGREGATION_CALC) THEN 
		!	WRITE(*,*) 'ERROR: P(0):',P(1), ' /= M(0):', M(1)
		!	STOP
        !END IF
	    ! Don't need this, just set P(0) to M(0)
		
		
		!ALLOCATE(P(1:P_Moments))
		
		
		
		If (ANY(M(1) > P_in)) THEN
			P(:) = M(1)!1.D0
		ELSE
			! P(0) = M(0)
			P(1) = M(1)
			! Input Moments
			P(2:) = P_in
        END IF 
		
		!IF(present(AGG)) THEN
		!	IF(AGG==1) THEN
		!		AGGREGATION = .FALSE.
		!	ELSEIF(AGG==2) THEN
		!		AGGREGATION = .TRUE.
		!	END IF 
		!END IF
		
		
        ! Initialize Nucleation, Coagulation, Surface Rates, And Aggregate Rates
	    W(:) = 0
	    G(:) = 0
	    R(:) = 0
        H(:) = 0
        Ragg(:) = 0
        rC2H2_nuc = 0
		rC2H2_surf = 0
		rC2H2 = 0 
        rCO = 0
        rH = 0
        rH2 = 0
        rH2O = 0
        rO2 = 0
        rOH = 0
		OXID = .TRUE.
		
        
        ! Check if MOMIC Calculation is enabled
        IF (.NOT. doSOOT) THEN
            WRITE(*,*) 'MOMENT CALCULATION TURNED OFF IN MOMIC INPUT FILE'
            RETURN
        END IF
		
		! Check Inputs are correct
		! THIS CAN CAUSE POTENTIAL MEMORY ISSUES IF THEY ARE NOT THE SAME
		! MAYBE LOOK AT EDITING HOW MOMENTS ARE INPUT?
		IF (SIZE(M) /= M_MOMENTS) THEN
			WRITE(*,*) 'ERROR: THE NUMBER OF INPUT M MOMENTS DOES NOT MATCH UP &
			WITH THE NUMBER OF MOMENTS SPECIFIED IN THE INPUT FILE'
			STOP
		END IF
		
		IF (AGGREGATION .AND. SIZE(P) /= P_MOMENTS) THEN
			WRITE(*,*) 'ERROR: THE NUMBER OF INPUT P MOMENTS DOES NOT MATCH UP &
			WITH THE NUMBER OF MOMENTS SPECIFIED IN THE INPUT FILE'
			STOP
		END IF
		    
        ! Calculate Fractional MOMENTS
		ALLOCATE(mu(mu_LO:mu_HI))
		
        ! Interpolate reduced moments
		CALL interpolate_mu(M,mu,mu_LO,mu_HI)		
		
		
		! Oxidation Calculation Control Flag
		IF (mu(6) <= OXID_RAT) THEN
			OXID = .FALSE.
		END IF
		
        ! Calculate Kc, Kc' and Kf Coagulation Parameters
		CALL COAGULATION_PARAMETERS(TEMPERATURE,PRESSURE)
		
		! Calculate Average soot diameter 
		D_soot = Soot_D_Calc(mu,mu_LO,mu_HI)
		
		! Calculate Knudsen Number
		kn = 2*MFP(PRESSURE,TEMPERATURE)*(1e-2)/D_soot
		
		! Kn number regimes
		!IF (kn < 0.1) THEN
		!	COAGULATION_REGIME = 0 ! Continuum
		!ELSEIF (kn < 10) THEN
		!	COAGULATION_REGIME = 1 ! Transitional
		!ELSEIF (kn >= 10) THEN
		!	COAGULATION_REGIME = 2 ! Free Molecular
		!ELSEIF (isNAN(kn)) THEN
		!	WRITE(*,*) "ERROR IN Temperature or Pressure input"
		!ELSE
		!	WRITE(*,*) "HOW DID YOU GET TO THIS CONDITION??? THIS SHOULD BE PHYISCALLY IMPOSSIBLE!!"
		!	STOP
		!END IF 
		  		  
		!IF (FORCE_REGIME) THEN
		!	COAGULATION_REGIME = REGIME
		!END IF         
		
        ! Calculate Primary Particle Nucleation Rates
		CALL NUCLEATION(M,TEMPERATURE,PRESSURE,C_C2H2,rC2H2_nuc,R,M_Moments)
		
		! If Number density is less than 1/cm^3 then only use Nucleation calculation,  skip the rest of calculation
		! Like in Ranjan's Code, this helps with some cases!!
		IF(M(1)<150D0) THEN
            WRITE(*,*) 'Soot Number Density M(0) < 150'
            WRITE(*,*) 'Only Calculating Nucleation Source'
			RETURN
		END IF 
					
		! Check if Aggregation is turned on
		
		!AGGREGATION = .FALSE.
		
		
        If(AGGREGATE) THEN
            CALL CALCULATE_SOURCE_AGG
        ELSE
            CALL CALCULATE_SOURCE_NO_AGG
        END IF
		
		! Get Total Rate of Nucleation and Surface Growth C2H2
		rC2H2 = rC2H2_nuc+rC2H2_surf
		

        ! ********************************************************
        ! ************* ERROR CHECKING ***************************
        ! ********************************************************
		! Check if rates are nan
        DO i = 1,M_Moments
			IF( isNAN(G(i)) .OR. isNAN(W(i)) .OR. isNAN(R(i))) THEN
				WRITE(*,*) 'NAN for Output Rates'
				!WRITE(*,*) 'R -- Gr -- Wr -- Rr'
				!WRITE(*,*) i, G(i), W(i), R(i)
				error_flag = .TRUE.
				!EXIT
            END IF
        END DO 
        
        ! Check if M(n) < M(0)
        !IF(ANY(M(2:)<M(1))) THEN
         !   WRITE(*,*) 'Error higher moments are smaller than 0th moment'
          !  error_flag = .True.
        !END IF
        
		CALL VERIFY_DISTRIBUTION(M,error_flag)
		
	
        ! DEBUGGING STUFF
        numCalls = numCalls + 1
		

		! Check Error Flag and Stop if error
		IF(error_flag) THEN
            WRITE(*,*) 'NumCalls = ',NumCalls
			WRITE(*,*) '*******INPUTS********'
			WRITE(*,*) 'M Moments'
			WRITE(*,"(6ES14.3)") M
			WRITE(*,*) 'P Moments'
			WRITE(*,"(6ES14.3)") P
			WRITE(*,*) 'Temperature - Pressure '
			WRITE(*,*) TEMPERATURE, PRESSURE
			WRITE(*,*) '[C2H2] - [H] - [H2]'
			WRITE(*,*) C_C2H2, C_H, C_H2
			WRITE(*,*) '[H2O] - [O2] - [OH]'
			WRITE(*,*) C_H2O, C_O2, C_OH
			WRITE(*,*) '******Calculations******'
			WRITE(*,*) 'mu'
			WRITE(*,*) mu
			WRITE(*,*) '---- Gr Rates --- '
			WRITE(*,*) G
			WRITE(*,*) '---- Wr Rates --- '
			WRITE(*,*) W
			WRITE(*,*) '---- Rr Rates --- '
			WRITE(*,*) R			
        END IF
		
		! f2py is not compiling if STOP is in the above if statement for some reason
		!IF(error_flag) THEN
		!	STOP
		!END IF 
		
        ! Normalize all Equations by M_ref
        R = R/M_ref
        G = G/M_ref
        W = W/M_ref
        Ragg = Ragg/M_ref
        H = H/M_ref

        ! Don't need this anymore so throw it away    		 
		DEALLOCATE(mu)
		!DEALLOCATE(P)
		
        RETURN
		
    CONTAINS
		SUBROUTINE VERIFY_DISTRIBUTION(M,error_flag)
			! Check Realizability of Distribution Function
			
			IMPLICIT NONE
			REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
			LOGICAL, INTENT(OUT) :: error_flag
			
			
			! Initialize error_flag
			error_flag = .FALSE.
			
			IF(M(0)<1) THEN
				WRITE(*,*) 'ERROR IN DISTRIBUTION: M(0) < 1'
				error_flag = .TRUE.
			ELSEIF(M(0)*M(2)/(M(1)*M(1))<1) THEN
				WRITE(*,*) 'ERROR in Distribution: M(0)M(2)/M(1)^2 < 1'
				error_flag = .TRUE.
			END IF
		END SUBROUTINE
	
		SUBROUTINE PRIMARY_COAGULATION(M,mu,mu_LO,mu_HI,G)
		
			IMPLICIT NONE
			REAL(KIND=DBL), INTENT(IN), DIMENSION(0:) :: M
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
			REAL(KIND=DBL), INTENT(OUT), DIMENSION(0:SIZE(M)-1) :: G
		
		
			IF(COAGULATION_REGIME==0) THEN
				CALL PRIMARY_COAGULATION_C(M,mu,mu_LO,mu_HI,G)
				WRITE(*,*) 'Coagulation is being calculated for Continuum Regime'
			ELSEIF(COAGULATION_REGIME==1) THEN
				CALL PRIMARY_COAGULATION_TR(M,mu,mu_LO,mu_HI,G)
				WRITE(*,*) 'Coagulation is being calculated for Transitional Regime'
			ELSEIF(COAGULATION_REGIME==2) THEN
				CALL PRIMARY_COAGULATION_FM(M,mu,mu_LO,mu_HI,G)
				WRITE(*,*) 'Coagulation is being calculated for Free Molecular Regime'
			ELSE
				WRITE(*,*) "WRONG COAGULATION REGIME FLAG, I DID SOMETHTHING WRONG?"
				STOP
			END IF 
		
		END SUBROUTINE PRIMARY_COAGULATION
		
		SUBROUTINE PRIMARY_COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,G)
			! Calculate Primary Coagulation with aggregation depending on COAGULATION_REGIME Flag
			IMPLICIT NONE
			REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M, P 
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), INTENT(IN), DIMENSION(mu_LO:mu_HI) :: mu
			REAL(KIND=DBL), INTENT(IN), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
			REAL(KIND=DBL), INTENT(OUT), DIMENSION(SIZE(M)) :: G
			
			IF(COAGULATION_REGIME==0) THEN
				CALL PRIMARY_COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,G)
				WRITE(*,*) 'Coagulation is being calculated for Continuum Regime'
			ELSEIF(COAGULATION_REGIME==1) THEN
				CALL PRIMARY_COAGULATION_AGGREGATE_TR(M,P,mu,mu_LO,mu_HI,pii,G)
				WRITE(*,*) 'Coagulation is being calculated for Transitional Regime'
			ELSEIF(COAGULATION_REGIME==2) THEN
				CALL PRIMARY_COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,G)
				WRITE(*,*) 'Coagulation is being calculated for Free Molecular Regime'
			ELSE
				WRITE(*,*) "WRONG COAGULATION REGIME FLAG, I DID SOMETHTHING WRONG"
				STOP
			END IF 	
		
		END SUBROUTINE PRIMARY_COAGULATION_AGGREGATE
		
		SUBROUTINE COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,H)
			! Calculate Aggregate Coagulation depending on COAGULATION_REGIME Flag
			
			IMPLICIT NONE
			
			!REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: M, P 
			REAL(KIND=DBL), INTENT(IN), DIMENSION(:) :: M,P
			INTEGER, INTENT(IN) :: mu_LO, mu_HI	! Need to specify these here otherwise f2py complains!
			REAL(KIND=DBL), DIMENSION(mu_LO:mu_HI), INTENT(IN) :: mu 
			REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2), INTENT(IN) :: pii
			REAL(KIND=DBL), DIMENSION(1:SIZE(P)-1), INTENT(OUT) :: H
			
			IF(COAGULATION_REGIME==0) THEN
				CALL COAGULATION_AGGREGATE_C(M,P,mu,mu_LO,mu_HI,pii,H)
				WRITE(*,*) 'Coagulation is being calculated for Continuum Regime'
			ELSEIF(COAGULATION_REGIME==1) THEN
				CALL COAGULATION_AGGREGATE_TR(M,P,mu,mu_LO,mu_HI,pii,H)
				WRITE(*,*) 'Coagulation is being calculated for Transitional Regime'
			ELSEIF(COAGULATION_REGIME==2) THEN
				CALL COAGULATION_AGGREGATE_FM(M,P,mu,mu_LO,mu_HI,pii,H)
				WRITE(*,*) 'Coagulation is being calculated for Free Molecular Regime'
			ELSE
				WRITE(*,*) "WRONG COAGULATION REGIME FLAG, I DID SOMETHTHING WRONG"
				STOP
			END IF 		
		
		END SUBROUTINE COAGULATION_AGGREGATE
	
        SUBROUTINE CALCULATE_SOURCE_NO_AGG()
                ! Calculate Source Terms without Aggregation
                IMPLICIT NONE
        
                WRITE(*,*) 'Calculating Moment Rates with out Aggregation'
                
				
				CALL PRIMARY_COAGULATION(M,mu,mu_LO,mu_HI,G)
				
				
				CALL SURFACE_GROWTH(M,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
									rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
				
				
				! If Oxidation Regime then base surface rates on M1, like in Ranjan's code

				IF(W(2) < 0) THEN
					
					WRITE(*,*) '*** IN OXIDATION REGIME, Basing Moments on M(1) ***'
					DO i = 1,SIZE(M)-1
						f(i) = log(mu(6*i))
					END DO
					
					
					CALL linear(SIZE(M)-1,f,a,b,rsq)
					
						
					DO i = 4,6*(SIZE(M)-1)-1,6
						mu(i) = EXP(a+b*i/6.D0)
					END DO 
					
					
					!WRITE(*,*) '*******ORIGINAL W********'
					!WRITE(*,*) W 
					
					
					CALL SURFACE_GROWTH(M,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
				rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
						
				
					DO i =2,SIZE(M)
						W(i) = W(i)*M(i)/M(1)/EXP(a+b*i)
					END DO
					

				
					!WRITE(*,*) '*****New W********'
					!WRITE(*,*) W
					!STOP
				
				END IF 
                
				
				! Turn off Source Terms Cutoff Criteria Like in S.P Roy, Haworth Paper
				! For Strongly Oxidizing Environment
				!IF (W(2)<0 .AND. M(2)<32*M(1)) THEN
				!	W(:) = 0
				!	RETURN
				!END IF 
				
				
            END SUBROUTINE CALCULATE_SOURCE_NO_AGG
            
            SUBROUTINE CALCULATE_SOURCE_AGG()
	            ! CALCULATE SOURCE TERMS FOR METHOD OF MOMENTS USING KAZAKOV AGGREGATION MODEL
	            ! ASSUMING Average Particle Diameter Dsoot > D*
                IMPLICIT NONE
            
                ! Local Variable for reduced P moment
                REAL(KIND=DBL), DIMENSION(0:SIZE(P)-1,-2:2) :: pii
                
                ! Initialize Variables
                pii(:,:) = 0
                
                WRITE(*,*) 'Calculating Moment Rates with Aggregation Enabled'
		
		
                IF(.TRUE.) THEN !IF(D_soot > d_star) THEN

	                ! Interpolate pi 
	                CALL interpolate_pii(P,pii)
                    
					
                    ! Calculate Aggregate Nucleation
                    Ragg(:) = R(1)
                    
                    ! If P terms are 0 then only use nucleation term and non-aggregate coagulation, skip rest of calculation
                    IF(ANY(P == 0)) THEN
                        CALL PRIMARY_COAGULATION(M,mu,mu_LO,mu_HI,G)
                        RETURN
                    END IF
                    
            
	                ! Calculate Coagulation
	                CALL PRIMARY_COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,G)

					! Calculate Aggregate Coagulation
					CALL COAGULATION_AGGREGATE(M,P,mu,mu_LO,mu_HI,pii,H)
			
					! Calculate Surface Growth 
					CALL SURFACE_GROWTH_AGG(M,P,pii,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
											rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
		
		
					
					! If Oxidation Regime then base surface rates on M1, like in Ranjan's code
					IF(W(2)<0) THEN
						DO i = 1,SIZE(M)-1
							f(i) = log(mu(3*i))
						END DO	
						
						CALL linear(SIZE(M)-1,f,a,b,rsq)
							
						DO i = 2,3*(SIZE(M)-1)-1,3
							mu(i) = EXP(a+b*i/3.D0)
						END DO 
						
						
						!WRITE(*,*) '*******ORIGINAL W********'
						!WRITE(*,*) W 
						
						CALL SURFACE_GROWTH_AGG(M,P,pii,mu,mu_LO,mu_HI,TEMPERATURE,C_C2H2,C_H,C_H2,C_H2O,C_O2,C_OH,&
											rC2H2_surf,rCO,rH,rH2,rH2O,rO2,rOH,W_C2H2,W_O2,W_OH,W)
					
					
						DO i =2,SIZE(M)
							W(i) = W(i)*M(i)/M(1)/EXP(a+b*i)
						END DO
						
						! Oxidation should also result in reduced number density as particles get oxidized
						! Original Formulation does not account for that, so using this to account for it
						! W0 = W_OXIDATION/Cmin
						! Alternatively maybe set M0 = M1/Cmin?
						
						W(0) = (W_O2(1)+W_OH(1))/10D0
						
						
						!WRITE(*,*) '*****New W********'
						!WRITE(*,*) W
						!STOP
					
					END IF 							
                                        
                ELSE
                    WRITE(*,*) 'Average Soot D < d*'
					WRITE(*,*) 'THIS CONDITION SHOULD NOT BE POSSIBLE!, CHECK CALCULATE_SOURCE_AGG'
					STOP
                    CALL CALCULATE_SOURCE_NO_AGG()
                END IF
        
				RETURN
            END SUBROUTINE CALCULATE_SOURCE_AGG
                        
	END SUBROUTINE CALCULATE_SOURCE_SWITCHER
	
	
	
	
	
	
    
    
END MODULE MOMIC
