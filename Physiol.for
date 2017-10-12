$DEBUG:'D' !DEBUG VERSION; ALL LINES BEGINNING WITH 'D' ARE COMPILED IN

C**********************************************************************
C This file contains all physiology calculations - photosynthesis,
C stomatal conductance, transpiration, respiration.
C
c The main subroutines (called externally) are:
C PSTRANSP - calls the other functions, does the leaf temp iteration calcn
C RESP - calculates maintenance respiration rate
C GRESP - calculates growth respiration rate
C CALCRMW - calculates stem maintenance respiration rate per unit mass
C CALCFBIOM - calculates foliar biomass from leaf area & SLA
C CALCWBIOM - calculates woody biomass from height & diameter
C For water balance: not finished
C   ETCAN - calculates canopy transpiration rate
C   PartitionPPT - partitions precip between drainage & canopy storage
C
C Subsidiary subroutines are
C 1. Photosynthesis
C   PHOTOSYN - calculates photosynthesis from the FvC model
C   GAMMAFN - calculates T dependence of Gamma*
C   KMFN - calculates T dependence of Km
C   ARRH - Arrhenius T dependence
C   JMAXTFN - calculates T dependence of Jmax
C   VCMAXTFN - calculates T dependence of Vcmax
C 2. Conductances  & transpiration
C   GBHFREE - calculates conductance to heat through free convection
C   GBHFORCED - calculates conductance to heat through forced convection
C   GRADIATION - calculates radiation conductance
C   GSJARVIS - calculates stomatal conductance using the Jarvis model
C   GBCAN - calculates boundary layer conductance of canopy
C   PENMON - implements the Penman-Monteith equation
C**********************************************************************


C**********************************************************************
      SUBROUTINE PSTRANSP(
     +  RDFIPT,TUIPT,TDIPT,RNET,WIND,
     +  PAR,TAIR,CA,RH,VPD,VMFD,PRESS,SOILMD,
     +  JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN,
     +  THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW,
     +  MODELGS,GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0,
     &  GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,
     +  G0,D0L,GAMMA,G1,
     +  WLEAF,NSIDES,CWATER,ITERMAX,
     +  GSC,ALEAF,RD,ET,EI,FHEAT,TLEAF,GBH
     +  )
C This subroutine calculates leaf photosynthesis and transpiration.
C These may be calculated by
C (1) assuming leaf temperature = air temperature, Cs = Ca and Ds = Da
C (2) using iterative scheme of Leuning et al (1995) (PC&E 18:1183-1200)
C to calculate leaf temp, Cs & Ca.
C Setting ITERMAX = 0 gives (1); ITERMAX > 0 (suggest 100) gives (2).
C If canopy is wet (CWATER > 0) assume TLEAF = TAIR and calculate both
C transpiration (ET) and evaporation (EI). NB At present CWATER is set to 0. 
C**********************************************************************

      INCLUDE 'maestcom'

      INTEGER MODELGS
      REAL JMAX25,I0,LHV

C Set initial values of leaf temp and surface CO2 & VPD
      TLEAF = TAIR
      DLEAF = VPD
      VMLEAF = VMFD
      RHLEAF = RH
      CS = CA

C Following calculations do not depend on TLEAF
C Latent heat of water vapour at air temperature (J mol-1)
      LHV = (H2OLV0 - 2.365E3 * TAIR) * H2OMW
C Const s in Penman-Monteith equation  (Pa K-1)
      SLOPE = (SATUR(TAIR + 0.1) - SATUR(TAIR)) / 0.1
C Radiation conductance (mol m-2 s-1)
      GRADN = GRADIATION(TAIR,RDFIPT,TUIPT,TDIPT)
C Boundary layer conductance for heat - single sided, forced convection
      GBHU = GBHFORCED(TAIR,PRESS,WIND,WLEAF)

C**********************************************************************
      ITER = 0  ! Counter for iterations - finding leaf temperature
100   CONTINUE  ! Return point for iterations

       CALL PHOTOSYN(
     +  PAR,TLEAF,CS,RHLEAF,DLEAF,VMLEAF,SOILMD,
     +  JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN,
     +  THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW,
     +  MODELGS,GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0,
     &  GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,
     +  G0,D0L,GAMMA,G1,
     +  GSC,ALEAF,RD
     +)

C Boundary layer conductance for heat - single sided, free convection
      GBHF = GBHFREE(TAIR,TLEAF,PRESS,WLEAF)
C Total boundary layer conductance for heat
      GBH = GBHU + GBHF

C Total conductance for heat - two-sided
      GH = 2.*(GBH + GRADN)
C Total conductance for water vapour
      GBV = GBVGBH*GBH
      GSV = GSVGSC*GSC
C      GV = NSIDES*(GBV*GSV)/(GBV+GSV) ! already one-sided value
      GV = (GBV*GSV)/(GBV+GSV)

C  Call Penman-Monteith equation
      ET = PENMON(PRESS,SLOPE,LHV,RNET,VPD,GH,GV)

C If canopy is wet, call Penman equation and finish
	IF (CWATER.GT.0.0) THEN
	  GVEVAP = 2.*GBH*GBVGBH	! Conductance to water vapour of wet surface
	  EI = PENMON(PRESS,SLOPE,LHV,RNET,VPD,GH,GVEVAP)
	  GOTO 200
	ELSE
	  EI = 0.0
	END IF

C End of subroutine if no iterations wanted.
      IF (ITERMAX.EQ.0) GOTO 200

C Otherwise, calculate new TLEAF, DLEAF, RHLEAF & CS
      GBC = GBH/GBHGBC
      CS = CA - ALEAF/GBC
      TDIFF = (RNET - ET*LHV) / (CPAIR * AIRMA * GH)
      TLEAF1 = TAIR + TDIFF
      DLEAF = ET * PRESS / GV
      RHLEAF = 1. - DLEAF/SATUR(TLEAF1)
      VMLEAF = DLEAF/PRESS*1E-3

C Check to see whether convergence achieved or failed
      IF (ABS(TLEAF - TLEAF1).LT.TOL) GOTO 200
      IF (ITER.GT.ITERMAX) THEN
        CALL SUBERROR('FAILED CONVERGENCE IN PSTRANSP',IFATAL,0)
      END IF

C Update temperature & do another iteration
      TLEAF = TLEAF1
      ITER = ITER + 1
      GOTO 100

200   ET = ET*1E6  ! Return ET,EI in umol m-2 s-1
	EI = EI*1E6 
      FHEAT = (TLEAF - TAIR)*2.*GBH*CPAIR*AIRMA

      RETURN
      END ! PsTransp


C**********************************************************************
      SUBROUTINE PHOTOSYN(
     +  PAR,TLEAF,CS,RH,VPD,VMFD,SOILMD,
     +  JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN,
     +  THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW,
     +  MODELGS,GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0,
     &  GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,
     +  G0,D0L,GAMMA,G1,
     +  GS,ALEAF,RD
     +)
C This subroutine calculates photosynthesis according to the ECOCRAFT
C agreed formulation of the Farquhar & von Caemmerer (1982) equations.
C Stomatal conductance may be calculated according to the Jarvis,
C Ball-Berry or BB-Leuning models.
C NB ALEAF is NET leaf photosynthesis. 
C NB The effect of soil water content is currently ignored - needs
C to be added as part of soil water balance routines.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL PAR,TLEAF,CS,RH,VPD,VMFD
      REAL JMAX25,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN
      REAL THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW
      REAL GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0
	REAL GSJA,GSJB,T0,TREF,TMAX
      REAL G0,D0L,GAMMA,G1
      REAL GS,ALEAF,RD
      INTEGER MODELGS,IQERROR,IECO

      REAL GAMMASTAR,KM,JMAX,VCMAX,J,VJ
      REAL A,B,C,AC,AJ,GSDIVA,CIC,CIJ
      REAL KMFN,JMAXTFN

C Calculate photosynthetic parameters from leaf temperature.
      GAMMASTAR = GAMMAFN(TLEAF,IECO)                   ! CO2 compensation point, umol mol-1
      KM = KMFN(TLEAF,IECO)                             ! Michaelis-Menten for Rubisco, umol mol-1
      JMAX = JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)      ! Potential electron transport rate, umol m-2 s-1
      VCMAX = VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)   ! Maximum Rubisco activity, umol m-2 s-1
      RD = RESP(RD0,TLEAF,Q10F,RTEMP,DAYRESP,TBELOW)    ! Day leaf respiration, umol m-2 s-1
      J = QUADM(THETA,-(AJQ*PAR+JMAX),AJQ*PAR*JMAX,IQERROR) ! Actual e- transport rate, umol m-2 s-1
      VJ = J/4.0                                        ! RuBP-regen rate, umol m-2 s-1

C Deal with extreme cases
	IF ((JMAX.LE.0.0).OR.(VCMAX.LE.0.0)) THEN
	  ALEAF = -RD
        IF ((MODELGS.NE.2).AND.(MODELGS.NE.3)) THEN		  ! Minimum Gs value, mol m-2 s-1	
	    GS = GSMIN
	  ELSE
	    GS = G0
	  END IF
	  RETURN
	END IF

C Calculation varies according to the stomatal model chosen.

C Jarvis model
      IF ((MODELGS.NE.2).AND.(MODELGS.NE.3)) THEN
        GS = GSJARVIS(MODELGS,PAR,I0,VPD,D0,VK1,VK2,VPD1,VPD2,VMFD,
     +      VMFD0,CS,GSJA,GSJB,SOILMD,SMD1,SMD2,
     +      TLEAF,T0,TREF,TMAX,GSREF,GSMIN)

        IF (GS.LE.GSMIN) THEN
          ALEAF = -RD
	    GS = GSMIN
          RETURN
        END IF

C Photosynthesis when Rubisco is limiting
        A = 1./GS
        B = (RD - VCMAX)/GS - CS - KM
        C = VCMAX * (CS - GAMMASTAR) - RD * (CS + KM)
        AC = QUADM(A,B,C,IQERROR)

        IF (IQERROR.EQ.1) THEN
          GS = GSMIN
          AC = - RD
        END IF

C Photosynthesis when electron transport is limiting
        B = (RD - VJ)/GS - CS - 2*GAMMASTAR
        C = VJ * (CS - GAMMASTAR) - RD * (CS + 2*GAMMASTAR)
        AJ = QUADM(A,B,C,IQERROR)

        IF (IQERROR.EQ.1) THEN
          GS = GSMIN
          AJ = - RD
        END IF

        ALEAF = AMIN1(AC,AJ)        ! Jarvis model solution

      ELSE
C Calculate soil moisture modifying factor
	  FSOIL = 1.0
	  IF (SOILMD.GT.0.0.AND.SMD1.GT.0.0) THEN
	    FSOIL = 1.0 - SMD1*EXP(SMD2*SOILMD)
	    IF (FSOIL.LT.0.0) FSOIL = 0.0
	  END IF

        IF (MODELGS.EQ.2) THEN
C Ball-Berry model
          GSDIVA = G1 * RH / (CS - GAMMA) * FSOIL
        ELSE IF (MODELGS.EQ.3) THEN
C Ball-Berry-Leuning model
          GSDIVA = G1 / (CS - GAMMA) / (1 + VPD/D0L) * FSOIL
        END IF

C Following calculations are used for both BB & BBL models.
C Solution when Rubisco activity is limiting

        A = G0 + GSDIVA * (VCMAX - RD)
        B = (1. - CS*GSDIVA) * (VCMAX - RD) + G0 * (KM - CS)
     +      - GSDIVA * (VCMAX*GAMMASTAR + KM*RD)
        C = -(1. - CS*GSDIVA) * (VCMAX*GAMMASTAR + KM*RD) - G0*KM*CS
        CIC = QUADP(A,B,C,IQERROR)

        IF ((IQERROR.EQ.1).OR.(CIC.LE.0.0).OR.(CIC.GT.CS)) THEN
          AC = 0.0
        ELSE
          AC = VCMAX * (CIC - GAMMASTAR) / (CIC + KM)
        END IF

C Solution when electron transport rate is limiting
        A = G0 + GSDIVA * (VJ - RD)
        B = (1. - CS*GSDIVA) * (VJ - RD) + G0 * (2.*GAMMASTAR - CS)
     +      - GSDIVA * (VJ*GAMMASTAR + 2.*GAMMASTAR*RD)
        C = -(1. - CS*GSDIVA) * GAMMASTAR * (VJ + 2.*RD)
     +      - G0*2.*GAMMASTAR*CS
        CIJ = QUADP(A,B,C,IQERROR)

        AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
	  IF (AJ-RD.LT.1E-6) THEN	  ! Below light compensation point
	    CIJ = CS
          AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
        END IF

        ALEAF = AMIN1(AC,AJ) - RD  ! Solution for Ball-Berry model
        GS = G0 + GSDIVA*ALEAF
        IF (GS.LT.G0) GS = G0

      END IF

      RETURN
      END !Photosyn

C**********************************************************************
      REAL FUNCTION GAMMAFN(TLEAF,IECO)
C This subroutine calculates Gamma(star), or the CO2 compensation point
C in the absence of non-photorespiratory respiration.
C This is the ECOCRAFT-agreed formulation of this function.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL TLEAF

      IF (IECO.EQ.1) THEN 
C Ecocraft fomulation; based on Brooks & Farquhar and von Caemmerer et al. 
C If TLEAF < -1.0 then calculate Gamma for T = -1 (quadratic not applicable)
        IF (TLEAF.LT.-1.0) THEN
          GAMMAFN = 36.9 + 1.88*(-26.0) + 0.036*(-26.0)*(-26.0)
        ELSE
          GAMMAFN = 36.9 + 1.88*(TLEAF-25) + 0.036*(TLEAF-25)*(TLEAF-25)
        END IF

	ELSE	! Bernacchi et al 2001 PCE 24: 253-260
	  GAMMAFN = ARRH(42.75,37830,TLEAF,25.0)
	ENDIF

      RETURN
      END !Gammafn

C**********************************************************************
      REAL FUNCTION KMFN(TLEAF,IECO)
C This subroutine calculates Km, or the effective Michaelis-Menten
C coefficient of Rubisco activity.
C This is the ECOCRAFT-agreed formulation of this function.
C**********************************************************************

      INCLUDE 'maestcom'
      REAL OI,KC25,KO25,KCEA,KOEA,KC,KO

      OI = 205000         ! Oxygen partial pressure (umol mol-1)
      IF (IECO.EQ.1) THEN
C Physiological constants - values agreed by Ecocraft - Badger & Collatz values
        KC25 = 404          ! MM coefft of Rubisco for CO2 (umol mol-1)
        KO25 = 248000       ! MM coefft of Rubisco for O2 (umol mol-1)
        KCEA = 59400        ! Temp. response of Kc (J mol-1)
        KOEA = 36000        ! Temp. response of Ko (J mol-1)
	ELSE !  Bernacchi et al 2001 PCE 24: 253-260
        KC25 = 404.9        ! MM coefft of Rubisco for CO2 (umol mol-1)
        KO25 = 278400       ! MM coefft of Rubisco for O2 (umol mol-1)
        KCEA = 79430        ! Temp. response of Kc (J mol-1)
        KOEA = 36380        ! Temp. response of Ko (J mol-1)
	END IF

C This function is well-behaved for TLEAF < 0.0
	  KC = ARRH(KC25,KCEA,TLEAF,25.0)
	  KO = ARRH(KO25,KOEA,TLEAF,25.0)
        KMFN = KC * (1. + OI/KO)

      RETURN
      END !KmFn


C**********************************************************************
      REAL FUNCTION JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)
C This subroutine calculates the potential electron transport rate
C (Jmax) at the leaf temperature.
C This is the ECOCRAFT-agreed formulation of this function.
C**********************************************************************

      INCLUDE 'maestcom'
      REAL JMAX25,TLEAF,EAVJ,EDVJ,DELSJ

C This function is well-behaved for TLEAF < 0.0
      JMAXTFN = JMAX25 * EXP((TLEAF-25)*EAVJ/(RCONST*TK(TLEAF)*TK(25.)))
     +  * (1.+EXP((DELSJ*TK(25.)-EDVJ)/(RCONST*TK(25.))))
     +  / (1.+EXP((DELSJ*TK(TLEAF)-EDVJ)/(RCONST*TK(TLEAF))))

C Function allowing Vcmax to be forced linearly to zero at low T - 
C introduced for Duke data
	IF (TLEAF.LT.TVJDN) THEN
	  JMAXTFN = 0.0 
      ELSE IF (TLEAF.LT.TVJUP) THEN
	  JMAXTFN = (TLEAF - TVJDN)/(TVJUP - TVJDN)*JMAXTFN
	END IF
	 
      RETURN
      END !JmaxTFn

C**********************************************************************
      REAL FUNCTION VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)
C This subroutine calculates the maximum Rubisco activity
C (Vcmax) at the leaf temperature.
C This is the ECOCRAFT-agreed formulation of this function.
C**********************************************************************

      INCLUDE 'maestcom'
      REAL VCMAX25,TLEAF,EAVC,EDVC,DELSC

C There is still disagreement as to whether this function has an
C optimum or not. Both forms are available here. If EDVC <= 0 (default 0)
C then the no-optimum form is used.
C Both functions are well-behaved for TLEAF < 0.0

      IF (EDVC.LE.0.0) THEN
        VCMAXTFN = VCMAX25 * EXP(EAVC*(TLEAF - 25)
     +  / (TK(25.)*RCONST*TK(TLEAF)))
      ELSE
        VCMAXTFN = VCMAX25
     +  * EXP((TLEAF-25.)*EAVC/(RCONST*TK(TLEAF)*TK(25.)))
     +  * (1.+EXP((DELSC*TK(25.)-EDVC)/(RCONST*TK(25.))))
     +  / (1.+EXP((DELSC*TK(TLEAF)-EDVC)/(RCONST*TK(TLEAF))))
      END IF

C Function allowing Vcmax to be forced linearly to zero at low T - 
C introduced for Duke data
	IF (TLEAF.LT.TVJDN) THEN
	  VCMAXTFN = 0.0 
      ELSE IF (TLEAF.LT.TVJUP) THEN
	  VCMAXTFN = (TLEAF - TVJDN)/(TVJUP - TVJDN)*VCMAXTFN
	END IF
	 
      RETURN
      END !VcmaxTFn


C**********************************************************************
      REAL FUNCTION RESP(RD0,TLEAF,Q10F,RTEMP,DAYRESP,TBELOW)
C This function calculates respiration from temperature
C using a Q10 (exponential) formulation.
C**********************************************************************

      REAL RD0,TLEAF,Q10F,RTEMP,DAYRESP

      IF (TLEAF.GE.TBELOW) THEN
        RESP = RD0 * EXP(Q10F * (TLEAF-RTEMP)) * DAYRESP
      ELSE
        RESP = 0.0
      END IF

      RETURN
      END !Resp

C**********************************************************************
      REAL FUNCTION ARRH(KT,EA,T,TREF)
C The Arrhenius function.
C KT is the value at Tref deg C; Ea the activation energy (J mol-1) and T the temp (deg C). 
C**********************************************************************

      INCLUDE 'MAESTCOM'
	REAL KT,EA,T,TREF

	ARRH = KT*EXP(EA*(T-TREF)/(RCONST*(T-ABSZERO)*(TREF-ABSZERO)))
	RETURN
	END !Arrh


C**********************************************************************
      REAL FUNCTION GSJARVIS(
     +  MODELGS,PAR,I0,VPD,D0,VK1,VK2,VPD1,VPD2,VMFD,VMFD0,CS,GSJA,GSJB,
     &  SOILMD,SMD1,SMD2,
     +  TLEAF,T0,TREF,TMAX,GSREF,GSMIN
     +  )
C Calculate stomatal conductance according to the Jarvis model.
C This model calculates gs by multiplying together functions of several
C environmental variables: light, VPD, CO2 & temperature.
C**********************************************************************

      INTEGER MODELGS
      REAL PAR,I0,VPD,D0,VK1,VK2,VMFD,VMFD0,CS,GSJA,GSJB
      REAL TLEAF,T0,TREF,TMAX,GSREF,GSMIN,P,SOILMD,SMD1,SMD2

C Defaults for the different factors
      FLIGHT = 1.0
      FVPD = 1.0
      FCO2 = 1.0
      FTEMP = 1.0
	FSOIL = 1.0

C Response to incident radiation - PAR & I0 in umol m-2 s-1
      IF (I0.GT.0.0) FLIGHT = PAR / (PAR + I0)

C Hyperbolic decline with VPD (VPD in Pa) * BRAY (ALEX BOSC)
C or Lohammer response to VPD (VPD & D0 in Pa) * ECOCRAFT
C or mole fraction deficit (VMFD in mmol mol-1) * MARK RAYMENT
C or linear decline with VPD (VPD1, VPD2 in Pa) * GROMIT (TIM RANDLE)
      IF (MOD(MODELGS,100).GE.30) THEN
	  IF (VPD.GT.0.0) FVPD = 1./(VK1*VPD**VK2)
	ELSE IF (MOD(MODELGS,100).GE.20) THEN
        IF (VPD.GE.VPD1) FVPD = 1. - (VPD - VPD1)/(VPD2 - VPD1)
      ELSE IF (MOD(MODELGS,100).GE.10) THEN
        IF (VMFD0.GT.0.0) FVPD = 1.0 - VMFD/VMFD0
      ELSE
        IF (D0.GT.0.0) FVPD = 1. / (1. + VPD/D0) 
      END IF
	IF (FVPD.GT.1.0) FVPD = 1.0
      IF (FVPD.LT.0.0) FVPD = 0.0

C Two possible responses to CO2 - linear or non-linear. CO2 in umol mol-1.
      IF (MOD(MODELGS,10).EQ.0) THEN
        IF (GSJA.NE.0.0) FCO2 = 1 - GSJA * (CS - 350.0)
      ELSE
        IF (GSJB.NE.0.0) FCO2 = (GSJB + 350.0) / (GSJB + CS)
      END IF
      IF (FCO2.LT.0.0) FCO2 = 0.0

C Temperature function is optional. To not use it, set TMAX <= 0.0. (Default).
      IF (TMAX.LE.0.0) THEN
        FTEMP = 1.0
      ELSE IF (MOD(MODELGS,1000).GE.100) THEN
        P = (TMAX - TREF)/(TREF - T0)
        FTEMP = (TLEAF - T0)*((TMAX - TLEAF)**P) / 
     &        ((TREF - T0)*((TMAX - TREF)**P))
      ELSE
        FTEMP = (TLEAF - T0) * (2*TMAX - T0 - TLEAF) /
     +          ((TREF - T0) * (2*TMAX - T0 - TREF))
      END IF      
      IF (FTEMP.LT.0.0) FTEMP = 0.0

C Soil moisture function according to Granier & Loustau 1994
C Will only be used if soil moisture data present and SMD1 > 0.
      IF (SOILMD.GT.0.0.AND.SMD1.GT.0.0) THEN
	  FSOIL = 1.0 - SMD1*EXP(SMD2*SOILMD)
	  IF (FSOIL.LT.0.0) FSOIL = 0.0
	END IF

C Multiply factors together.
C Gsref is in mol m-2 s-1; Gsjarvis required in mol m-2 s-1.
      GSJARVIS = (GSREF-GSMIN)*FLIGHT*FVPD*FCO2*FTEMP*FSOIL + GSMIN

      RETURN
      END !GsJarvis


C**********************************************************************
      REAL FUNCTION ETCAN(
     +  WIND,ZHT,Z0HT,ZPD,PRESS,TAIR,RNET,VPD,GSCAN,STOCKING
     +)
C Calculate transpiration by applying Penman-Monteith to whole canopy.
C Returns umol m-2 s-1.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL LHV

C Get boundary layer conductance
      GB = GBCAN(WIND,ZHT,Z0HT,ZPD,PRESS,TAIR)
      GSV = GSCAN*GSVGSC*STOCKING !convert mol CO2/tree/s to mol H2O/m2/s
      RNETM2 = RNET*STOCKING 
      
      IF (GB*GSV.GT.0.0) THEN 
C Latent heat of water vapour at air temperature (J mol-1)
        LHV = (H2OLV0 - 2.365E3 * TAIR) * H2OMW
C Const s in Penman-Monteith equation  (Pa K-1)
        SLOPE = (SATUR(TAIR + 0.1) - SATUR(TAIR)) / 0.1
C Call Penman-Monteith
        GH = GB
        GV = 1./(1./GSV + 1./GB)
        ETCAN = PENMON(PRESS,SLOPE,LHV,RNETM2,VPD,GH,GV)*1E6
      ELSE
        ETCAN = 0.0
      END IF

      RETURN
      END !ETCan


C**********************************************************************
      REAL FUNCTION PENMON(
     +  PRESS,SLOPE,LHV,RNET,VPD,GH,GV
     +)
C This subroutine calculates evapotranspiration by leaves using the Penman-Monteith equation. 
C Inputs:	PRESS atmospheric pressure, Pa
C		SLOPE slope of VPD/T curve, Pa K-1
C		LHV latent heat of water at air T, J mol-1
C		RNET net radiation, J m-2 s-1
C		VPD vapour pressure deficit of air, Pa
C		GH boundary layer conductance to heat (free & forced & radiative components), mol m-2 s-1
C		GV conductance to water vapour (stomatal & bdry layer components), mol m-2 s-1
C Result in mol H2O m-2 s-1.
C**********************************************************************

      INCLUDE 'maestcom'
      REAL LHV

      GAMMA = CPAIR*AIRMA*PRESS/LHV

	IF (GV.GT.0.0) THEN
        ET = (SLOPE * RNET + VPD * GH * CPAIR * AIRMA) /
     &		(SLOPE + GAMMA * GH/GV)
      ELSE
        ET = 0.0
      END IF
      PENMON = ET / LHV

      RETURN
      END !PenMon


C**********************************************************************
      REAL FUNCTION GRADIATION(TAIR,RDFIPT,TUIPT,TDIPT)
C Returns the 'radiation conductance' at given temperature.
C Formula from Ying-Ping's version of Maestro.
C See also Jones (1992) p. 108.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      GRADIATION = 4.*SIGMA*(TK(TAIR)**3.)
     +  * RDFIPT/TDIPT * EMLEAF * (TDIPT + TUIPT) / (CPAIR * AIRMA)

      RETURN
      END !GRadiation


C**********************************************************************
      REAL FUNCTION GBHFORCED(TAIR,PRESS,WIND,WLEAF)
C Boundary layer conductance for heat - single sided, forced convection
C in mol m-2 s-1
C See Leuning et al (1995) PC&E 18:1183-1200 Eqn E1
C**********************************************************************

      INCLUDE 'MAESTCOM'

      CMOLAR = PRESS / (RCONST * TK(TAIR))
      GBHFORCED = 0.003 * SQRT(WIND/WLEAF) * CMOLAR

      RETURN
      END !GBHForced


C**********************************************************************
      REAL FUNCTION GBHFREE(TAIR,TLEAF,PRESS,WLEAF)
C Boundary layer conductance for heat - single sided, free convection
C in mol m-2 s-1
C See Leuning et al (1995) PC&E 18:1183-1200 Eqns E3 & E4
C**********************************************************************

      INCLUDE 'MAESTCOM'

      CMOLAR = PRESS / (RCONST * TK(TAIR))
      IF ((TLEAF-TAIR).NE.0.0) THEN
        GRASHOF = 1.6E8 * ABS(TLEAF-TAIR) * (WLEAF**3.) ! Grashof number
        GBHFREE = 0.5 * DHEAT * (GRASHOF**0.25) / WLEAF * CMOLAR
      ELSE
        GBHFREE = 0.0
      END IF

      RETURN
      END !GBHFree


C**********************************************************************
      REAL FUNCTION GBCAN(WIND,ZHT,Z0HT,ZPD,PRESS,TAIR)
C Canopy boundary layer conductance (from Jones 1992 p 68) 
C in mol m-2 s-1
C**********************************************************************

      INCLUDE 'MAESTCOM'

      IF (Z0HT.GT.0.0) THEN

C Formula from Jones 1992 p 68
        GBCAN = WIND*(VONKARMAN**2)/(LOG((ZHT - ZPD)/Z0HT))**2
C Convert from mm s-1 to mol m-2 s-1
        CMOLAR = PRESS / (RCONST * TK(TAIR))
        GBCAN = GBCAN * CMOLAR 

      ELSE 
        GBCAN = 0.0
      END IF

      RETURN
      END !GBCan


C**********************************************************************
      SUBROUTINE CALCWBIOM(IDAY,HT,DIAM,COEFFT,EXPONT,WINTERC,
     &  WBIOM,WBINC)
C Calculate the woody biomass (kg DW) on the given day from the height
C (m) and diameter (m). Also calculate the increment in woody biomass
C since previous day (g DW). Needed to calculate woody respiration.
C**********************************************************************

      PREVWBIOM = WBIOM
      WBIOM = COEFFT*HT*(DIAM**EXPONT) + WINTERC
      IF (IDAY.EQ.0) PREVWBIOM = WBIOM
      WBINC = (WBIOM - PREVWBIOM)*1E3

      RETURN
      END ! CalcWBiom

C**********************************************************************
      SUBROUTINE CALCFBIOM(IDAY,NOSLADATES,FOLLAY,SLA,PROP,NOLAY,NOAGEP,
     &  FBIOM,FBINC)
C Calculate foliage biomass from SLA and leaf area - done in layers.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      REAL FOLLAY(MAXLAY)
      REAL SLA(MAXLAY,MAXC)
      REAL PROP(MAXC)

      IF (NOSLADATES.GT.0) THEN
        PREVFBIOM = FBIOM
        FBIOM = 0.0
        DO 10 I = 1,NOLAY
          DO 10 J = 1,NOAGEP
            FBIOM = FBIOM + FOLLAY(I)*PROP(J)/SLA(I,J)
10      CONTINUE
        IF (IDAY.EQ.0) PREVFBIOM = FBIOM
        FBINC = (FBIOM - PREVFBIOM)*1E3
      ELSE
        FBIOM = 0.0
        FBINC = 0.0
      END IF

      RETURN
      END !CalcFBiom


C**********************************************************************
      FUNCTION CALCRMW(MODELRW,COLLA,COLLK,STEMSDW,
     &  DIAM,HT,STEMFORM,RMWAREA,WBIOM,RMW)
C Calculate stem respiration rate per unit biomass if necessary
C**********************************************************************

      INCLUDE 'MAESTCOM'

      IF (MODELRW.EQ.1) THEN
        RMW = COLLA*EXP(COLLK*DIAM*100.0)*STEMSDW
	ELSE IF (MODELRW.EQ.2) THEN
	  STEMAREA = STEMFORM*PI*(DIAM**2)*HT
	  RMW = RMWAREA*STEMAREA/WBIOM
      END IF

      CALCRMW = RMW
	RETURN
	END !CALCRMW


C**********************************************************************
      FUNCTION GRESP(BINC,EFFY)
C Calculate the growth respiration. Use increment in biomass
C (g DW tree-1 d-1) and the efficiency of growth (g g-1 C).
C Returns a value in mol CO2 tree-1 d-1.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      IF (BINC.GT.0.0) THEN
        GRESP = BINC * EFFY / GCPERMOL * CPERDW
      ELSE
        GRESP = 0.0
      ENDIF

      RETURN
      END ! GResp



*********************************************************************
      SUBROUTINE PartitionPPT(PPT,PFLA,CANCAPLA,FOLT,CanopyWater)
C Calculate fraction of incoming gross precipitation 
C intercepted by crown of target tree (ignoring stem part)
*********************************************************************

C     INPUTS   PPT,PFLA,CANCAPLA,FOLT,CanopyWater
C     OUTPUTS  CanopyWater

C interception is a proportion, Pf, of gross precipitation, PPT
C Pf is related to leaf area
      PF = PFLA*FOLT
      PPTint = PF*PPT

C crown input, PPTc, increases the depth of water on the canopy, Cc
C up to a maximum value CANCAP, which is related to leaf area
      CanopyWater = AMIN1(CanopyWater + PPTint,CANCAPLA*FOLT)

	RETURN
      END ! PartitionPPT


C************************************************************************
      SUBROUTINE EvTran(Etpot,Eipot,CANCAPLA,FOLT,CanopyWater,Et,Ei)
C Evapo-transpiration: calculate transpiration,evaporation of
C intercepted water and drainage, dependent upon the depth of water
C on the leaves
C ETpot, EIpot in umol tree-1 s-1
C CanopyWater, CanCAP in mmol tree-1
C************************************************************************

      INCLUDE 'MAESTCOM'

      CanopyCapacity = CANCAPLA*FOLT			! Canopy water capacity
      IF (CanopyWater.GT.0.0) THEN
	  PossEvap = CanopyWater/SPERHR*1E3		! Maximum evaporation rate (given water on crown)
        Ratio = CanopyWater/CanopyCapacity	! Ratio is between 0 and 1 - fraction of canopy which is wet
        Ei = AMIN1(Eipot*Ratio,PossEvap)		! Actual evaporation rate
        Et = Etpot*(1-Ratio)
        CanopyWater = AMAX1(CanopyWater - Ei*SPERHR*1E-3,0.0)
      ELSE
        Ei = 0.0
        Et = Etpot
      END IF

      RETURN
      END !EvTran

