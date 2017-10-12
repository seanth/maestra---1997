$DEBUG:'D' !DEBUG VERSION; ALL LINES BEGINNING WITH 'D' ARE COMPILED IN

C**********************************************************************
C RADN.FOR
C This file contains all the canopy structure and radiation interception
C routines. The major routines are:

C POINTS - sets up the grid points
C SUN - calculates the daylength
C ZENAZ - calculates zenith and azimuth angles
C SLOPES - calculates corrections for slope of plot (to soil reflectance)
C EXDIFF, EXBEAM - calculate the extinction coefficients for diffuse and
C   beam radiation, respectively
C TRANSD, TRANSB - calculate the transmittances of diffuse and beam radiation
C EHC - sets up the equivalent horizontal canopy (used in scattering calculations)
C SCATTER - calculates scattered radiation
C ABSRAD - calculates absorbed radiation
C CATEGO - used to set up PAR histogram

C These subroutines call the following additional routines:
C POINTS
C    SEGMT, CUMUL, BETA, SURFACE
C SUN
C   ALUNAR, ANOM, ECCENT, EPSIL, OMEGA, DAYJUL
C EXBEAM, EXDIFF
C   COSDEL
C TRANSB, TRANSD
C   TREDST, DISTIN, DIST, CDIST, WPATH, SHADED, TESTD, POSSIBLE
C EHC
C   ASSIGN, CHART, TDUMIN
C SCATTER
C   ABSTHERM
C**********************************************************************


C**********************************************************************
      SUBROUTINE POINTS(
     &  NUMPNT,JLEAF,JSHAPE,SHAPE,RXntr,RYntr,RZntr,
     &  ZBCntr,DXTntr,DYTntr,DZTntr,FOLntr,PROP,
     &  BPT,NOAGEC,NOAGEP,NOLAY,
     &  XL,YL,ZL,VL,DLT,DLI,LGP,FOLLAY
     &  )
C This subroutine is used to set up to 120 grid points through
C the crown. There are 12 grid points per layer and a minimum of
c 3 layers (36 points) is recommended. It is also recommended that
C the number of grid points used is a multiple of 36.
C The inputs required are:
C NUMPNT: the number of gridpoints
C JLEAF: 0 - no leaf area dist; 1 - vertical only; 2 - horiz. & vert.
C JSHAPE,SHAPE: indicate crown shape
C RX,RY,RZ: radius & green crown length of target crown
C ZBC: height to base of crown in target crown
C DXT,DYT,DZT: x,y,z co-ordinates of target crown
C FOL: leaf area of target crown
C BPT: coefficients of beta distributions of leaf area
C NOAGEC: no of age classes for which beta distributions are specified
C PROP: proportion of leaf area in each age class
C NOLAY: no of layers of crown
C Routine outputs are:
C XL,YL,ZL: the co-ordinates of each grid point
C VL: the volume of crown associated with each grid point
C DLT, DLI: the amount of leaf area associated with each grid point
C LGP: the physiological layer corresponding to each grid point
C**********************************************************************

      INCLUDE 'maestcom'

      REAL BPT(8,MAXC),PROP(MAXC)
      REAL DLT(MAXP),DLI(MAXC,MAXP)
      REAL XL(MAXP),YL(MAXP),ZL(MAXP),VL(MAXP),AX(MAXP/4),CZ(MAXP/4)
      INTEGER LGP(MAXP)
      REAL FOLLAY(MAXLAY)

C Constants for use later in the subroutine
C     Number of grid points within quarter of a layer (PPQ)
      PPQ = PPLAY/4
C     Relative height of each height interval (HTINT)
      HTINT=1.0/(NOLAY*1.)
C     Radius of crown at 45 degrees to axis.
      RX2=RXntr*RXntr
      RY2=RYntr*RYntr
      IF (JSHAPE.EQ.JBOX) THEN
        RXY = SQRT(RX2+RY2)
      ELSE
	  RXY = SQRT((RX2+RY2)/2.)
C        RXY = SQRT(2.0*(RX2*RY2)/(RX2+RY2)) // only correct if RX = RY, and then pointless BM 11/99
      END IF
C     Number of grid points in each quadrant.
      MUMPNT=NUMPNT/4

C 1. Calculate co-ordinates (XL,YL,ZL) and volume (VL) for each grid point.

C (a) for one quadrant of the crown, set relative heights (CZ), radial
C distances (AX) and volume (VL) of each grid point.
      DO 10 LAYER = 1, NOLAY       ! for each layer
        IPQ = (LAYER-1)*PPQ + 1     ! grid point number
C calculate relative height - same for each grid point in the layer
        CZ(IPQ)=(LAYER-0.5)*HTINT   ! 3 grid points in each quadrant
        CZ(IPQ+1)=CZ(IPQ)
        CZ(IPQ+2)=CZ(IPQ)

C calculate relative radial distance - differs in odd & even layers
C The function 'surface' finds relative radial distance to crown surface as function of relative height
        IF(MOD(LAYER,2).EQ.1) THEN
          CORR=SURFACE(CZ(IPQ),JSHAPE)
        ELSE
          CORR=0.5*SQRT(2.0)*SURFACE(CZ(IPQ),JSHAPE) ! Multiply by cos 45°
        ENDIF
        AX(IPQ)=0.288675*CORR       ! Factors come from assuming volume of each gridpoint is equal
        AX(IPQ+1)=0.696923*CORR     
        AX(IPQ+2)=0.908248*CORR
C Calculate volume for each grid point.
        HUP= LAYER*HTINT*RZntr
        HDN= (LAYER-1)*HTINT*RZntr
C        DH = HUP-HDN // BM 11/99 not used
        VLM = SEGMT(JSHAPE,HUP,HDN,RXNTR,RYNTR,RZNTR)
        VL(IPQ)  = VLM/PPLAY
        VL(IPQ+1)= VL(IPQ)
        VL(IPQ+2)= VL(IPQ)
10    CONTINUE

C (b) Now set x,y,z co-ordinates of each grid point, using AX & CZ from above.
      DO 20 IPQ1 = 1,MUMPNT    ! loop over gridpoints in 1st quadrant
        IPQ2 = IPQ1+MUMPNT     ! gridpoints in 2nd quadrant
        IPQ3 = IPQ1+MUMPNT*2   ! ditto 3rd quadrant
        IPQ4 = IPQ1+MUMPNT*3   ! ditto 4th quadrant
        IF ((MOD(IPQ1,6).LE.3).AND.(MOD(IPQ1,6).NE.0)) THEN  ! Odd layers
          XL(IPQ1)= AX(IPQ1)*RXntr+DXTntr
          XL(IPQ2)= 0.0         +DXTntr
          XL(IPQ3)=-AX(IPQ1)*RXntr+DXTntr
          XL(IPQ4)= 0.0         +DXTntr
          YL(IPQ1)= 0.0         +DYTntr
          YL(IPQ2)= AX(IPQ1)*RYntr+DYTntr
          YL(IPQ3)= 0.0         +DYTntr
          YL(IPQ4)=-AX(IPQ1)*RYntr+DYTntr
        ELSE                                             ! Even layers
          XL(IPQ1)= AX(IPQ1)*RXY+DXTntr
          XL(IPQ2)=-AX(IPQ1)*RXY+DXTntr
          XL(IPQ3)= XL(IPQ2)
          XL(IPQ4)= XL(IPQ1)
          YL(IPQ1)= AX(IPQ1)*RXY+DYTntr
          YL(IPQ2)= YL(IPQ1)
          YL(IPQ3)=-AX(IPQ1)*RXY+DYTntr
          YL(IPQ4)= YL(IPQ3)
        ENDIF
        ZL(IPQ1)=CZ(IPQ1)*RZntr+ZBCntr+DZTntr
C Height & volume is the same for every grid point in the one layer.
        ZL(IPQ2)=ZL(IPQ1)
        ZL(IPQ3)=ZL(IPQ1)
        ZL(IPQ4)=ZL(IPQ1)
        VL(IPQ2)=VL(IPQ1)
        VL(IPQ3)=VL(IPQ1)
        VL(IPQ4)=VL(IPQ1)
20    CONTINUE

C 2. Calculate the leaf area density associated with each grid point.
C (a) for uniform leaf area density, this is easy.
      IF (JLEAF.EQ.0) THEN
        AVGSD=FOLntr/(PI*RXntr*RYntr*RZntr*SHAPE)
	  CORFL = 1.0
        DO 30 IPT=1,NUMPNT
          DLT(IPT)=AVGSD
          DO 30 IAGE=1,NOAGEP
            DLI(IAGE,IPT)=DLT(IPT)*PROP(IAGE)
30      CONTINUE

C (b) less easy when LAD is specified with beta-distributions
      ELSE
        DO 40 IAGE = 1,NOAGEC       ! for each age class
          IF (JLEAF.EQ.2) THEN      ! calc. horizontal factors
            HORIZ1 = CUMUL(0.000000,0.577350,2,BPT,IAGE)
            HORIZ2 = CUMUL(0.577350,0.816497,2,BPT,IAGE)
            HORIZ3 = CUMUL(0.816497,1.000000,2,BPT,IAGE)
	      CORFL = TWOPI*FOLntr/4.0
          ELSE
            HORIZ1 = 1.0
            HORIZ2 = 1.0
            HORIZ3 = 1.0
	      CORFL = FOLntr/REAL(PPLAY)
          END IF
C set DLI for gridpoints in 1st quadrant
          DO 50 LAYER = 1,NOLAY      ! for each layer
            IPQ = (LAYER-1)*PPQ + 1   ! grid point number - 1st quadrant
            DOWN=(LAYER-1)*HTINT      ! calc. vertical factors
            UP=LAYER*HTINT
            VERT = CUMUL(DOWN,UP,1,BPT,IAGE)
            DLI(IAGE,IPQ)  = VERT*HORIZ1*PROP(IAGE)/VL(IPQ)
            DLI(IAGE,IPQ+1)= VERT*HORIZ2*PROP(IAGE)/VL(IPQ+1)
            DLI(IAGE,IPQ+2)= VERT*HORIZ3*PROP(IAGE)/VL(IPQ+2)
50      CONTINUE
C now set for gridpoints in other quadrants
        DO 40 IPQ1 = 1,MUMPNT    ! loop over gridpoints in 1st quadrant
          IPQ2 = IPQ1+MUMPNT     ! gridpoints in 2nd quadrant
          IPQ3 = IPQ1+MUMPNT*2   ! ditto 3rd quadrant
          IPQ4 = IPQ1+MUMPNT*3   ! ditto 4th quadrant
          DLI(IAGE,IPQ2) = DLI(IAGE,IPQ1)
          DLI(IAGE,IPQ3) = DLI(IAGE,IPQ1)
          DLI(IAGE,IPQ4) = DLI(IAGE,IPQ1)
40      CONTINUE  ! Finish looping over age classes as well as gridpoints
      END IF ! If JLEAF = 0

C Calculate total leaf area density for each grid point
      DO 60 IPT = 1,NUMPNT
        DLT(IPT) = 0.0
        DO 60 IAGE = 1,NOAGEC
          DLT(IPT) = DLT(IPT) + DLI(IAGE,IPT)
60    CONTINUE


C Normalize the subvolume and foliage of each grid point.

C BM 11/99 Following is unnecessary. 
C      TOTVL=0.0
C      DO 70 IPT=1,NUMPNT
C        TOTVL=VL(IPT)+TOTVL
C70    CONTINUE
C      CORVL=PI*RXntr*RYntr*RZntr*SHAPE/TOTVL
C      DO 80 IPT=1, NUMPNT
C         VL(IPT)=VL(IPT)*CORVL
C80    CONTINUE

C BM 12/99 Calculate CORFL directly - provides check on function
C      TOTFL=0.0
C      DO 90 IPT=1, NUMPNT
C         TOTFL=DLT(IPT)*VL(IPT)+TOTFL
C90    CONTINUE
C      IF (TOTFL.EQ.0.0) THEN
C        CORFL = 0.0
C      ELSE
C        CORFL=FOLntr/TOTFL
C      END IF
     
      DO 100 IPT=1,NUMPNT
        DLT(IPT)=DLT(IPT)*CORFL
C        DO 100 IAGE = 1,NOAGEC
C          DLI(IAGE,IPT)=DLI(IAGE,IPT)*CORFL
100   CONTINUE

C Re-calculate DLI: it must correspond to NOAGEP
C BM 11/99 Moved this part to after correction for total leaf area
      IF (NOAGEP.EQ.1) THEN
        DO 65 IPT = 1,NUMPNT
          DLI(1,IPT) = DLT(IPT)
65      CONTINUE
      ELSE
        DO 66 IPT = 1,NUMPNT
          DO 66 IAGE = 1,NOAGEP
            DLI(IAGE,IPT) = PROP(IAGE)*DLT(IPT)
66      CONTINUE
      END IF

C Set the layer of each grid point (LGP)
      DO 110 IPT = 1,MUMPNT
         LGP(IPT) = NOLAY - INT(CZ(IPT)*REAL(NOLAY))
         IF (LGP(IPT).EQ.0) LGP(IPT)=1
110   CONTINUE
      DO 120 J=1,MUMPNT
      DO 120 IJ=1,3
         IPT=J+IJ*MUMPNT
         LGP(IPT)=LGP(J)
120   CONTINUE

C Calculate leaf area in each layer of the tree
      DO 130 ILAY = 1,MAXLAY
        FOLLAY(ILAY) = 0.0
130   CONTINUE
      DO 140 IPT = 1,NUMPNT
        FOLLAY(LGP(IPT)) = FOLLAY(LGP(IPT)) + DLT(IPT)*VL(IPT)
140   CONTINUE

      RETURN
      END !Points


C**********************************************************************
      FUNCTION SEGMT(JSHAPE,HUP,HDN,RX1,RY1,RZ1)
C This subroutine calculates the volume of a segment of the crown
C between relative crown heights HUP and HDN. 
C Inputs: JSHAPE = shape of crown; RX1,RY1,RZ1 = dimensions of crown. 
C**********************************************************************

      INCLUDE 'maestcom'

      IF (JSHAPE.EQ.JHELIP) THEN
        V1 = PI*RX1*RY1* (HUP- (HUP**3)/ (3.0*RZ1*RZ1))
        V2 = PI*RX1*RY1* (HDN- (HDN**3)/ (3.0*RZ1*RZ1))
        SEGMT = V1 - V2

      ELSE IF (JSHAPE.EQ.JFELIP) THEN
        V1 = PI*RX1*RY1* (HUP-((HUP-RZ1/2.)**3)/(3.0*(RZ1/2.)**2))
        V2 = PI*RX1*RY1* (HDN-((HDN-RZ1/2.)**3)/(3.0*(RZ1/2.)**2))
        SEGMT = V1 - V2

      ELSE IF (JSHAPE.EQ.JCONE) THEN
        V1 = PI*RX1*RY1*(HUP-(HUP**2/RZ1)+(HUP**3)/(3*RZ1**2))
        V2 = PI*RX1*RY1*(HDN-(HDN**2/RZ1)+(HDN**3)/(3*RZ1**2))
        SEGMT = V1 - V2

      ELSE IF (JSHAPE.EQ.JPARA) THEN
        V1 = PI*RX1*RY1*(HUP-(HUP**2)/(2*RZ1))
        V2 = PI*RX1*RY1*(HDN-(HDN**2)/(2*RZ1))
        SEGMT = V1 - V2

      ELSE IF (JSHAPE.EQ.JCYL) THEN
        SEGMT = PI*RX1*RY1*(HUP-HDN)

      ELSE IF (JSHAPE.EQ.JBOX) THEN
        SEGMT = 4.*RX1*RY1*(HUP-HDN)

      END IF

      RETURN
      END !Segmt


C**********************************************************************
      FUNCTION CUMUL(DOWN,UP,JFUN,BPT,IAGE)
C Use the Gaussian numerical integration method to integrate beta function over gridpoint volume
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL BPT(8,MAXC)

C    WRITE (*,*) 'CALLING CUMUL'

      IF (JFUN.EQ.1) THEN
         B1 = BPT(1,IAGE)
         B2 = BPT(2,IAGE)
         B3 = BPT(3,IAGE)
	   B4 = BPT(4,IAGE)
      ELSE
         B1 = BPT(5,IAGE)
         B2 = BPT(6,IAGE) + 1   ! Because of r term in integral 
         B3 = BPT(7,IAGE)
	   B4 = BPT(8,IAGE)
      END IF

      GI = (UP-DOWN)/20.0
      GS1 = DOWN + GI*0.42265
      GS2 = DOWN + GI*1.5774
      CUMUL = 0.0
      DO 100 I = 1,10
        S11 = BETA(B1,B2,B3,B4,GS1)
        S12 = BETA(B1,B2,B3,B4,GS2)
        CUMUL = CUMUL + S11 + S12
        GS1 = GS1 + 2.0*GI
        GS2 = GS2 + 2.0*GI
  100 CONTINUE
      CUMUL = CUMUL*GI 

      RETURN
      END  !Cumul


C**********************************************************************
      FUNCTION SURFACE(H,JSHAPE)
C Find the relative radial distance to the surface of the tree crown 
C as a function of relative tree height.
C H is the relative height in the crown.
C JSHAPE is the crown shape. 
C**********************************************************************

      INCLUDE 'MAESTCOM'

      IF (JSHAPE.EQ.JCONE) THEN
         SURFACE = (1.- H)

      ELSE IF (JSHAPE.EQ.JHELIP) THEN
         SURFACE = SQRT(1. - H**2)

      ELSE IF (JSHAPE.EQ.JPARA) THEN
         SURFACE = SQRT(1. - H)
 
      ELSE IF (JSHAPE.EQ.JFELIP) THEN
         SURFACE = SQRT(1. - ((H-1./2.)**2)/((1./2.)**2))

      ELSE IF (JSHAPE.EQ.JCYL) THEN
         SURFACE = 1.

      ELSE IF (JSHAPE.EQ.JBOX) THEN
         SURFACE = 1.

      END IF

      RETURN
      END  !Surface


C**********************************************************************
      SUBROUTINE EXDIFF(NALPHA,ALPHA,FALPHA,NZEN,DIFZEN,RANDOM,DEXT)
C    calculate the extinction coefficients for the diffuse radiation
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DIFZEN(MAXANG),ALPHA(MAXANG),FALPHA(MAXANG),DEXT(MAXANG)

      DO 10 IZEN = 1,NZEN
         DEXT(IZEN) = 0.0
         SUM = 0.0
         DO 5 IALP = 1,NALPHA
            SUM = SUM +
     +      COSDEL(RANDOM,DIFZEN(IZEN),ALPHA(IALP))*FALPHA(IALP)
5        CONTINUE
         DEXT(IZEN) = SUM
10    CONTINUE

      RETURN !Exdiff
      END


C**********************************************************************
      FUNCTION COSDEL(RANDOM,THETA,ALPHA)
C FUNCTION TO calculate COSINE DELTA WHERE THETA IS SOLAR ZENITH
C ANGLE AND ALPHA IS INCLINATION ANGLE OF DLEAF.
C**********************************************************************

      INCLUDE 'maestcom'

      IF ((THETA+ALPHA-PID2).LE.0.00000) THEN
         COSDEL = COS(ALPHA)*COS(THETA)
      ELSE
         BP0 = ACOS(1./ (TAN(ALPHA)*TAN(THETA)))
         COSDEL = ((PID2-BP0)*COS(ALPHA)*COS(THETA)+
     +            SIN(ALPHA)*SIN(THETA)*SIN(BP0))/PID2
      END IF

      COSDEL = COSDEL*RANDOM

      RETURN
      END  !COSDEL


C**********************************************************************
      SUBROUTINE TRANSD(
     &  IDAY,IOTUTD,NEWCANOPY,IPROG,NT,XSLOPE,YSLOPE,
     &  NZEN,DIFZEN,NAZ,NUMPNT,DEXT,DIFSKY,
     &  XL,YL,ZL,RX,RY,RZ,DXT,DYT,DZT,
     &  XMAX,YMAX,SHADEHT,
     &  FOLT,ZBC,JLEAF,BPT,NOAGEC,PROP,JSHAPE,SHAPE,
     &  NEWTUTD,TU,TD,RELDF
     &  )
C This subroutine calculates the diffuse transmittances, if required.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DXT(MAXT),DYT(MAXT),DZT(MAXT),RX(MAXT),RY(MAXT),
     +     RZ(MAXT),ZBC(MAXT),FOLT(MAXT)
      REAL BPT(8,MAXC),PROP(MAXC)
      REAL DIFZEN(MAXANG),DEXT(MAXANG)
      REAL XL(MAXP),YL(MAXP),ZL(MAXP)
      REAL TU(MAXP),TD(MAXP),RELDF(MAXP)
      LOGICAL SHADED

C If IOTUTD = 0, then the diffuse transmittances are read in from file.
C This only needs to be done on the first day of the simulation.
C After reading in the transmittances, the subroutine ends.

      NEWTUTD = 0
      IF (IOTUTD.EQ.0) THEN
        IF (IDAY.EQ.0) THEN
            REWIND (UTUTD)
            READ (UTUTD,*) (J,TU(I),TD(I),RELDF(I),I =1,NUMPNT)
            NEWTUTD = 1
        END IF
        RETURN

C The transmittances are only calculated every IOTUTD'th day. 
C If the canopy has changed significantly they will also need recalculation.
C Check to see whether this needs to be done, otherwise return.

      ELSE IF ((MOD(IDAY,IOTUTD).NE.0).AND.(NEWCANOPY.NE.1)) THEN
        RETURN
      END IF

C Calculate transmittances (this is the old function TRANSD).

      NEWTUTD = 1
      DA = TWOPI/FLOAT(NAZ)

      DO 900 IPT = 1,NUMPNT
        WNS = 0.
        WDS = 0.
        WNS1= 0.
        WDS1= 0.
        WNG = 0.
        WDG = 0.

        DO 500 J = 1,NZEN
          DZUP = DIFZEN(J)
          DZDN = PI - DZUP
          SINUP = SIN(DZUP)
          COSUP = COS(DZUP)
          SINDN = SIN(DZDN)
          COSDN = COS(DZDN)
          ADDUP = 0.00
          ADDDN = 0.00
        
          DO 400 K = 1,NAZ
            DAZ = DA * (K-0.5)
            SLOPE = ATAN(COS(DAZ)*TAN(XSLOPE)+SIN(DAZ)*TAN(YSLOPE))
            IF ((PID2-DZUP).LT.SLOPE) GO TO 300
            IF (SHADED(XL(IPT),YL(IPT),ZL(IPT),
     &          DZUP,DAZ,XMAX,YMAX,SHADEHT)) GO TO 300

            IFLAG = 1
            CALL TREDST(IFLAG,IPROG,DZUP,DAZ,
     &         XL(IPT),YL(IPT),ZL(IPT),RX,RY,RZ,DXT,DYT,DZT,
     &         FOLT,ZBC,JLEAF,BPT,NOAGEC,PROP,JSHAPE,SHAPE,
     &         NT,SU,SU1)
            ADDUP = ADDUP + EXP(-DEXT(J)*(SU+SU1))

  300       IFLAG = 2
            CALL TREDST(IFLAG,IPROG,DZDN,DAZ,
     &         XL(IPT),YL(IPT),ZL(IPT),RX,RY,RZ,DXT,DYT,DZT,
     &         FOLT,ZBC,JLEAF,BPT,NOAGEC,PROP,JSHAPE,SHAPE,
     &         NT,SL,SL1)
            ADDDN = ADDDN + EXP(-DEXT(J)*(SL+SL1))

  400     CONTINUE

C        TMPVAR=(1.0+DIFSKY*COSUP)/(1.0+2.0*DIFSKY/PI)
          TMPVAR=(1.0+DIFSKY*COSUP)/(1.0+DIFSKY)
          PTHUP = ADDUP*SINUP*TMPVAR/FLOAT(NAZ)  
          WNS = WNS + COSUP*PTHUP
          WNS1 = WNS1 + DEXT(J)*PTHUP
          WDS = WDS + SINUP*COSUP*TMPVAR
          WNG = WNG + ADDDN*SINDN*COSDN/FLOAT(NAZ)
          WDG = WDG + SINDN*COSDN
  500   CONTINUE

C DIFFUSE TRANSMITTANCES FOR UPPER AND LOWER HEMISPHERES
        TD(IPT) = WNS/WDS
        TU(IPT) = WNG/WDG
        RELDF(IPT)= WNS1/WDS
  900 CONTINUE

      RETURN
      END !Transd


C**********************************************************************
      SUBROUTINE TREDST(IFLAG,IPROG,DZ,DAZ,
     &  XPT,YPT,ZPT,RX,RY,RZ,DXT,DYT,DZT,
     &  FOLT,ZBC,JLEAF,BPT,NOAGEC,PROP,JSHAPE,SHAPE,
     &  NT,S,S1)
C this subroutine is used to calculate the weighted pathlength,
C which is also kernal function of radiation penetration
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DXT(MAXT),DYT(MAXT),DZT(MAXT)
      REAL RX(MAXT),RY(MAXT),RZ(MAXT),ZBC(MAXT),FOLT(MAXT)
      REAL BPT(8,MAXC),PROP(MAXC)
      LOGICAL POSSIBLE

      TANAZ = TAN(DAZ)
      SINAZ = SIN(DAZ)

C Zero the pathlengths.
C S1 is distance in target tree, S is distance in other trees.
C Note: in normal program, gridpoints are all in target tree. In
C test program, gridpoints may be outside. Flag IPROG indicates this.
      S = 0.00
      S1 = 0.00

C  Loop over each tree to see if it affects the ray passing through the point.
      DO 200 M = 1,NT

C  choose the largest of the x and y radii for use in a quick test to
C  see if the tree could possibly be in the way.
        PATH = 0.00
        XTPOS = DXT(M) - XPT
        YTPOS = DYT(M) - YPT

        IF ((M.EQ.1).AND.(IPROG.NE.ITEST)) THEN
          CALL DISTIN(JSHAPE,IFLAG,DZ,DAZ,XPT,YPT,ZPT,
     &      RX(M),RY(M),RZ(M),ZBC(M),DXT(M),DYT(M),DZT(M),
     &      PATH,X1,Y1,Z1,X2,Y2,Z2)
        ELSE IF (POSSIBLE(XTPOS,YTPOS,RX(M),RY(M),SINAZ,TANAZ,DAZ)) THEN
          CALL DIST(JSHAPE,IFLAG,DZ,DAZ,XPT,YPT,ZPT,
     &      RX(M),RY(M),RZ(M),ZBC(M),DXT(M),DYT(M),DZT(M),
     &      PATH,X1,Y1,Z1,X2,Y2,Z2)
        END IF

        IF (PATH.EQ.0.00) GO TO 200

        IF (JLEAF.EQ.0) THEN
          AVGDL = FOLT(M)/ (RX(M)*RY(M)*RZ(M)*SHAPE*PI)
          SS = PATH*AVGDL
        ELSE
          CALL WPATH(JSHAPE,SHAPE,JLEAF,BPT,NOAGEC,PROP,
     &      X1,Y1,Z1,X2,Y2,Z2,PATH,
     &      RX(M),RY(M),RZ(M),ZBC(M),DXT(M),DYT(M),DZT(M),FOLT(M),SS)
        END IF

        IF ((M.GT.1).OR.(IPROG.EQ.ITEST)) THEN
              S = S + SS
            ELSE
              S1 = S1 + SS
            END IF

200   CONTINUE

      RETURN
      END  !Tredst


C**********************************************************************
      subroutine WPATH(JSHAPE,SHAPE,JLEAF,BPT,NOAGEC,PROP,
     &  X1,Y1,Z1,X2,Y2,Z2,
     &  PATH,RXQ,RYQ,RZQ,ZBCQ,DXTQ,DYTQ,DZTQ,FOLTQ,SS)
C This subroutine weighs the pathlength: SS is the weighted path.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL BPT(8,MAXC),PROP(MAXC)

      SS = 0.000

C Sample path at 20 points
      DO 100 I = 1,20

C Calculate co-ordinates of the 20 sample points
        ZZ = Z1 + (I-0.5)* (Z2-Z1)/20.0000 - ZBCQ - DZTQ	! Height of sample point
        RCH = ZZ/RZQ										! relative height of sample point
        IF ((RCH.LT.1.01) .AND. (RCH.GT.1.0)) RCH = 1.00	! numerical tidying
        TRD = SURFACE(RCH,JSHAPE)*SQRT(RXQ*RYQ)			! TRD is mean distance to crown surface 

C Find beta in vertical direction
        BETA1 = 0.0
        DO 10 J = 1,NOAGEC								! normalized leaf area density
          BETA1 = BETA1 + 
     &		BETA(BPT(1,J),BPT(2,J),BPT(3,J),BPT(4,J),RCH) * PROP(J)
10      CONTINUE
	  
        DFT = BETA1*FOLTQ/(TRD*TRD*RZQ)	
	  				
	  IF (JLEAF.EQ.1) THEN
	    DFT = DFT/PI

C Find beta in horizontal direction
        ELSE IF (JLEAF.EQ.2) THEN
          XX = X1 + (I-0.5)* (X2-X1)/20.000 - DXTQ
          YY = Y1 + (I-0.5)* (Y2-Y1)/20.000 - DYTQ
          TEMPRD = SQRT((XX)**2+ (YY)**2)
          IF (RXQ.EQ.RYQ) THEN
		  RR = TEMPRD/TRD
		ELSE
	      IF (XX.EQ.0.0.AND.YY.GE.0.0) THEN
		    ANG = PID2
		  ELSE IF (XX.EQ.0.0.AND.YY.LT.0.0) THEN
		    ANG = -PID2
		  ELSE
		    ANG=ATAN(YY/XX)
		  END IF
		  TRDE = SURFACE(RCH,JSHAPE)*
     &	      SQRT((RXQ*COS(ANG))**2+(RYQ*SIN(ANG))**2)
	      RR = TEMPRD/TRDE
		END IF
          IF (RR.GT.1.00 .AND. (RR.LE.1.10)) RR = 1.00
          BETA2 = 0.0
          DO 20 J = 1,NOAGEC
            BETA2 = BETA2 +
     +          BETA(BPT(5,J),BPT(6,J),BPT(7,J),BPT(8,J),RR)*PROP(J)
20        CONTINUE
          DFT = DFT*BETA2
        ENDIF

        SS = SS + DFT*PATH/20.00

  100 CONTINUE

      RETURN
      END  !Wpath


C**********************************************************************
      subroutine DIST(Jshape,IFLAG,DZZ,DAZ,
     &                XPP,YPP,ZPP,RXQ,RYQ,RZQ,ZBCQ,DXTQ,DYTQ,DZTQ,
     &                PATH,X1,Y1,Z1,X2,Y2,Z2)
C this subroutine is used to calculate the pathlength inside all
C trees around the tree we do all these calculations for
C**********************************************************************

      INCLUDE 'MAESTCOM'

C The default path length is zero.
      PATH = 0.00

C Calculate multipliers for converting from polar to rectangular co-ords.
      POL = SIN(DZZ)*COS(DAZ)
      POM = SIN(DZZ)*SIN(DAZ)
      PON = COS(DZZ)

C Check the path does not lie completely above or below the tree of interest.
      IF (IFLAG.EQ.2 .AND. ZPP.LE. (ZBCQ+DZTQ)) GO TO 100
      IF (IFLAG.NE.2 .AND. ZPP.GT. (ZBCQ+RZQ+DZTQ)) GO TO 100

C Find co-ords of the point of intersection between the ray and the bottom
C plane of the tree crown.
      XX = XPP + POL* (ZBCQ+DZTQ-ZPP)/PON
      YY = YPP + POM* (ZBCQ+DZTQ-ZPP)/PON
      ZZ = ZBCQ + DZTQ

C Two situations:
C 1. The point on the bottom plane lies within the tree crown
C 2. The point on the bottom plane lies outwith the tree crown
C For full ellipsoid, case 2. always holds.
      ISITU = 2
      IF (JSHAPE.NE.JFELIP) THEN
        DISTRADIAL = ((XX-DXTQ)/RXQ)**2 + ((YY-DYTQ)/RYQ)**2
        IF (DISTRADIAL.LT.1.00001) ISITU = 1
      END IF

C Find co-efficients of quadratic for distance between points on crown
C surface and grid point, then solve quadratic.
      call CDIST(JSHAPE,XPP,YPP,ZPP,RXQ,RYQ,RZQ,ZBCQ,
     &                 DXTQ,DYTQ,DZTQ,POL,POM,PON,DELTA,R1,R2)
      IF (DELTA.LE.0.00) GO TO 100       ! No intersection of path with crown

C Handle situation where point on bottom plane (XX,YY,ZZ) is within crown
C Points of interest are (XX,YY,ZZ) and one of two other points.
      IF (ISITU.EQ.1) THEN
        X1 = XX
        Y1 = YY
        Z1 = ZZ
        X2 = XPP + POL*R1
        Y2 = YPP + POM*R1
        Z2 = ZPP + PON*R1
        X3 = XPP + POL*R2
        Y3 = YPP + POM*R2
        Z3 = ZPP + PON*R2
        IF (Z2.GT. (Z1+RZQ) .OR. (Z2.LT.Z1)) THEN
           X2 = X3
           Y2 = Y3
           Z2 = Z3
        END IF

C Handle situation where point on bottom plane (XX,YY,ZZ) is without crown
C Two points of interest given by solns to quadratic.
      ELSE
        X1 = XPP + POL*R1
        Y1 = YPP + POM*R1
        Z1 = ZPP + PON*R1
        X2 = XPP + POL*R2
        Y2 = YPP + POM*R2
        Z2 = ZPP + PON*R2
      END IF

C Swap point1(x1,y1,z1) with point2(x2,y2,z2) if necessary to ensure
C that point2 is always above the point1
      IF (Z1.GT.Z2) THEN
         TEMP = Z2
         Z2 = Z1
         Z1 = TEMP
         TEMP = Y2
         Y2 = Y1
         Y1 = TEMP
         TEMP = X2
         X2 = X1
         X1 = TEMP
      END IF

C If the grid point is between the two intersection points, two possibilities:
C 1. The crowns overlap, in which case take the grid point as one point
C 2. The path is completely outwith the crown, in which case pathlength = 0
      IF ((Z1.LT.ZPP).AND.(ZPP.LT.Z2)) THEN
	  IF ((Z1.LT.ZZ).AND.(Z2.GT.ZZ+RZQ)) GOTO 100 !Case 2
        IF (IFLAG.NE.2) THEN
           X1 = XPP
           Y1 = YPP
           Z1 = ZPP
        ELSE
           Z2 = ZPP
           Y2 = YPP
           X2 = XPP
        END IF
      END IF

C Check the path is not completely below the grid point
      IF (IFLAG.EQ.2 .AND. (Z1.GE.ZPP)) GOTO 100
      IF (IFLAG.NE.2 .AND. (Z2.LE.ZPP)) GOTO 100

C Check the points of intersection are not completely above
C or below the tree crown
      IF ((Z2-0.0001).LE.ZZ .OR. (Z1+0.0001).GE.(RZQ+ZZ)) GOTO 100

C Calculate path length!
      PATH = SQRT((X1-X2)**2+ (Y1-Y2)**2+ (Z1-Z2)**2)

  100 RETURN
      END !Dist


C**********************************************************************
      SUBROUTINE DISTIN(JSHAPE,IFLAG,DZZ,DAZ,XPP,YPP,ZPP,
     &                  RXQ,RYQ,RZQ,ZBCQ,DXTQ,DYTQ,DZTQ,
     &                  PATH,X1,Y1,Z1,X2,Y2,Z2)
C this subroutine is used to calculate the pathlength inside the
C the crown of the tree we are concerned with.
C**********************************************************************

C The default path length is zero.
      PATH = 0.00

C Calculate multipliers for converting from polar to rectangular co-ords.
      POL = SIN(DZZ)*COS(DAZ)
      POM = SIN(DZZ)*SIN(DAZ)
      PON = COS(DZZ)

C Find co-efficients of quadratic for distance between points on crown
C surface and point on bottom plane, then solve quadratic.
      call CDIST(JSHAPE,XPP,YPP,ZPP,RXQ,RYQ,RZQ,ZBCQ,
     &                 DXTQ,DYTQ,DZTQ,POL,POM,PON,DELTA,R1,R2)

C Swap point1(x1,y1,z1) with point2(x2,y2,z2) if necessary to ensure
C that point2 is always above the point1
      IF ((PON*R2).LT.(PON*R1)) THEN
        TEMP = R1
        R1   = R2
        R2   = TEMP
      END IF

C Calculate co-ordinates of intersection of path with crown
      X1 = XPP + POL*R1
      Y1 = YPP + POM*R1
      Z1 = ZPP + PON*R1
      X2 = XPP + POL*R2
      Y2 = YPP + POM*R2
      Z2 = ZPP + PON*R2

C If both points are above the grid point, then take the lower point
C and the point of intersection with the bottom plane.
      IF (Z1.GT.ZPP .AND. Z2.GT.ZPP) THEN
         X2 = X1
         Y2 = Y1
         Z2 = Z1
         Z1 = ZBCQ + DZTQ
         Y1 = YPP + POM* (ZBCQ+DZTQ-ZPP)/PON
         X1 = XPP + POL* (ZBCQ+DZTQ-ZPP)/PON
      END IF

C If the lower point is below the bottom of the crown, take the
C point of intersection with the bottom plane as the lower point.
      IF (Z1.LT. (ZBCQ+DZTQ)) THEN
         Z1 = ZBCQ + DZTQ
         Y1 = YPP + POM* (ZBCQ+DZTQ-ZPP)/PON
         X1 = XPP + POL* (ZBCQ+DZTQ-ZPP)/PON
      END IF

C Discard lower point (if direction is from top) or upper point (if
C direction is from bottom).
      IF (IFLAG.NE.2) THEN
        X1 = XPP
        Y1 = YPP
        Z1 = ZPP
      ELSE
        X2 = XPP
        Y2 = YPP
        Z2 = ZPP
      END IF

C Calculate path length!
      PATH = SQRT((X1-X2)**2+ (Y1-Y2)**2+ (Z1-Z2)**2)

      RETURN
      END !Distin


C**********************************************************************
      SUBROUTINE CDIST(JSHAPE,XPP,YPP,ZPP,RXQ,RYQ,RZQ,ZBCQ,
     &                 DXTQ,DYTQ,DZTQ,POL,POM,PON,DELTA,R1,R2)
C this is subroutine to calculate the coefficients of the quadratic
C equation in the distance calculation subroutines, the values of
C these coefficients are dependent on the shape of the crown.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      IF (JSHAPE.EQ.JCONE) THEN
         A = (POL/RXQ)**2 + (POM/RYQ)**2 - (PON/RZQ)**2
         B = 2.00* ((XPP-DXTQ)*POL/ (RXQ**2)+
     +       (YPP-DYTQ)*POM/ (RYQ**2)+(RZQ+DZTQ+ZBCQ-ZPP)*PON/(RZQ**2))
         C = ((XPP-DXTQ)/RXQ)**2 + ((YPP-DYTQ)/RYQ)**2 -
     +       ((RZQ+DZTQ+ZBCQ-ZPP)/RZQ)**2

      ELSE IF (JSHAPE.EQ.JHELIP) THEN
         A = (POL/RXQ)**2 + (POM/RYQ)**2 + (PON/RZQ)**2
         B = 2.00* ((XPP-DXTQ)*POL/ (RXQ**2)+
     +       (YPP-DYTQ)*POM/ (RYQ**2)+(ZPP-DZTQ-ZBCQ)*PON/(RZQ**2))
         C = ((XPP-DXTQ)/RXQ)**2 + ((YPP-DYTQ)/RYQ)**2 +
     +       ((ZPP-DZTQ-ZBCQ)/RZQ)**2 - 1.00

      ELSE IF (JSHAPE.EQ.JPARA) THEN
         A = (POL/RXQ)**2 + (POM/RYQ)**2
         B = 2.00* ((XPP-DXTQ)*POL/ (RXQ**2)+(YPP-DYTQ)*POM/(RYQ**2))
     +      + PON/RZQ
         C = ((XPP-DXTQ)/RXQ)**2 + ((YPP-DYTQ)/RYQ)**2
     +      + (ZPP-DZTQ-ZBCQ)/RZQ - 1.00

      ELSE IF (JSHAPE.EQ.JFELIP) THEN
         A = (POL/RXQ)**2 + (POM/RYQ)**2 + (PON/(RZQ/2.))**2
         B = 2.00* ((XPP-DXTQ)*POL/ (RXQ**2)+
     +      (YPP-DYTQ)*POM/ (RYQ**2)+
     +      (ZPP-DZTQ-ZBCQ-RZQ/2.)*PON/ ((RZQ/2.)**2))
         C = ((XPP-DXTQ)/RXQ)**2 + ((YPP-DYTQ)/RYQ)**2 +
     +       ((ZPP-DZTQ-ZBCQ-RZQ/2.)/(RZQ/2.))**2 - 1.00

      ELSE IF (JSHAPE.EQ.JCYL) THEN
         A = (POL/RXQ)**2 + (POM/RYQ)**2
         B = 2.00* ((XPP-DXTQ)*POL/(RXQ**2) +
     +       (YPP-DYTQ)*POM/(RYQ**2))
         C = ((XPP-DXTQ)/RXQ)**2 + ((YPP-DYTQ)/RYQ)**2 - 1.00

      END IF

      IF (JSHAPE.NE.JBOX) THEN
        DELTA = B**2 - 4.0*A*C
        IF (DELTA.LT.0.00 .AND. (DELTA.GE.-0.01)) DELTA = 0.00
        IF (DELTA.GT.0.00) THEN
          R1 = (-B+SQRT(DELTA))/ (2.00*A)
          R2 = (-B-SQRT(DELTA))/ (2.00*A)
        END IF

      ELSE ! Solve for BOX shape

C Find shortest distance in each direction where line meets planes of box
        DELTA = -1.0
        IF (POL.NE.0.0) THEN
          D = (DXTQ + RXQ - XPP)/POL
          CALL TESTD(D,YPP+D*POM,DYTQ,RYQ,
     &      ZPP+D*PON,DZTQ+ZBCQ+RZQ/2,RZQ/2,
     &      R1,R2,DELTA)
          D = (DXTQ - RXQ - XPP)/POL
          CALL TESTD(D,YPP+D*POM,DYTQ,RYQ,
     &      ZPP+D*PON,DZTQ+ZBCQ+RZQ/2,RZQ/2,
     &      R1,R2,DELTA)
        END IF
        IF (POM.NE.0.0) THEN
          D = (DYTQ + RYQ - YPP)/POM
          CALL TESTD(D,XPP+D*POL,DXTQ,RXQ,
     &      ZPP+D*PON,DZTQ+ZBCQ+RZQ/2,RZQ/2,
     &      R1,R2,DELTA)
          D = (DYTQ - RYQ - YPP)/POM
          CALL TESTD(D,XPP+D*POL,DXTQ,RXQ,
     &      ZPP+D*PON,DZTQ+ZBCQ+RZQ/2,RZQ/2,
     &      R1,R2,DELTA)
        END IF
        IF (PON.NE.0.0) THEN
          D = (DZTQ + ZBCQ + RZQ - ZPP)/PON
          CALL TESTD(D,XPP+D*POL,DXTQ,RXQ,YPP+D*POM,DYTQ,RYQ,
     &      R1,R2,DELTA)
          D = (DZTQ + ZBCQ - ZPP)/PON
          CALL TESTD(D,XPP+D*POL,DXTQ,RXQ,YPP+D*POM,DYTQ,RYQ,
     &      R1,R2,DELTA)
        END IF
        IF ((DELTA.GT.-1.0).AND.(DELTA.LT.1.0)) 
     &    CALL SUBERROR('ERROR: FINDING INTERSECTION PTS OF BOX',
     &    IFATAL,0)

      END IF
        
      RETURN
      END !Cdist


C**********************************************************************
      SUBROUTINE TESTD(DLEN,C1,DQ1,RQ1,C2,DQ2,RQ2,
     &  R1,R2,DELTA)
C To find the intersection points of the ray with the tree crown when
C it is box shape: find the intersection point with each of the six
C planes of the box (done in CDIST) then test to see if the co-ords
C of each intersection lie within the box (done here).
C**********************************************************************

      IF ((ABS(C1-DQ1).LE.RQ1).AND.(ABS(C2-DQ2).LE.RQ2)) THEN
        IF (DELTA.GT.-1.0) THEN
          R2 = DLEN
        ELSE
          R1 = DLEN
        END IF
        DELTA = DELTA + 1.
      END IF
    
      RETURN
      END !TestD


C**********************************************************************
      LOGICAL FUNCTION SHADED(XPT,YPT,ZPT,ZEN,AZ,XMAX,YMAX,SHADEHT)
C This function calculates the height of a ray in direction (ZEN,AZ)
C passing through (XPT,YPT,ZPT) at the edge of the plot. This is used
C to check against the height of shadecloth around the plot (for those
C with shadecloth!) to determine if the ray enters the plot or not.
C**********************************************************************

      INCLUDE 'MAESTCOM'

C If there's not shadecloth, return straightaway
      SHADED = .FALSE.
      IF (SHADEHT.LE.0.0) RETURN

C A default large value
      VLARGE = 1E7

C Calculate distance from gridpoint to edge of plot, assuming rectangular.
      COSAZ = COS(AZ)
      SINAZ = SIN(AZ)

      IF (COSAZ.GT.0.0) THEN
        DX = (XMAX-XPT)/COSAZ
      ELSE IF (COSAZ.LT.0.0) THEN
        DX = XPT/(-COSAZ)
      ELSE
        DX = VLARGE
      END IF

      IF (SINAZ.GT.0.0) THEN
        DY = YPT/SINAZ
      ELSE IF (SINAZ.LT.0.0) THEN
        DY = (YMAX-YPT)/(-SINAZ)
      ELSE
        DY = VLARGE
      END IF

      DEDGE = AMIN1(DX,DY)

C Calculate height of ray at the edge & check it's above shadecloth
      EDGEHT = ZPT + DEDGE/(TAN(ZEN))

      IF (EDGEHT.LT.SHADEHT) SHADED = .TRUE.

      RETURN
      END !EdgeHt


C**********************************************************************
      FUNCTION POSSIBLE(XT,YT,RX,RY,SINAZ,TANAZ,AZ)
C  This function tests to see if the tree at xt,yt could possibly be in
C  the path of the ray at angle az through the point at xp,yp.
C  It is designed to be a  fast filter to remove those trees which
C  are obviously not in the path of the ray.
C  It returns true if the tree is possibly in the path of the ray and
C  false if it is impossible for the tree to be in the path.
C**********************************************************************

      INCLUDE 'maestcom'

      LOGICAL POSSIBLE

      POSSIBLE = .TRUE.
      RMAX = AMAX1(RX,RY)
      RMAXSINAZ = RMAX/SINAZ
      YDTANAZ = YT/TANAZ

      IF ((AZ.GT.0) .AND. (AZ.LT.PI)) THEN
         IF (XT.LT. (YDTANAZ-RMAXSINAZ)) POSSIBLE = .FALSE.
         IF (XT.GT. (YDTANAZ+RMAXSINAZ)) POSSIBLE = .FALSE.
      ELSE IF ((AZ.GT.PI) .AND. (AZ.LT.TWOPI)) THEN
         IF (XT.GT. (YDTANAZ-RMAXSINAZ)) POSSIBLE = .FALSE.
         IF (XT.LT. (YDTANAZ+RMAXSINAZ)) POSSIBLE = .FALSE.
      END IF

      RETURN
      END !Possible


C**********************************************************************
      SUBROUTINE EHC(NUMPNT,TU,TD,
     &  TOTLAI,XSLOPE,YSLOPE,NAZ,NZEN,DIFZEN,DEXT,
     &  DLAI,EXPDIF,LAYER,MLAYER
     &  )
C This subroutine sets up the equivalent horizontal canopy.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL TU(MAXP),TD(MAXP)
      REAL DIFZEN(MAXANG),DEXT(MAXANG)
      REAL RTA(50)
      INTEGER LAYER(MAXP),MLAYER(MAXP)

      CALL TDUMIN(TU,TD,NUMPNT,TAUMIN)

      CALL CHART(XSLOPE,YSLOPE,NAZ,NZEN,DIFZEN,DEXT,
     &  TOTLAI,TAUMIN,
     &  RTA,NLAY,DLAI,EXPDIF)

      CALL ASSIGN(TU,TD,NUMPNT,NLAY,RTA,
     &  LAYER,MLAYER)

      RETURN
      END !Ehc


C**********************************************************************
      subroutine TDUMIN(TD,TU,NUMPNT,TAUMIN)
C The transmittance of diffuse radiation through one elementary layer
C in the EHC is taken as the minimum of the diffuse transmittances of
C all grid points. This subroutine finds this minimum transmittance.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL TD(MAXP),TU(MAXP)

      TAUMIN = 1.
      DO 10 IPT = 1,NUMPNT
         SMALL = AMIN1(TU(IPT),TD(IPT))
         IF (SMALL.LT.TAUMIN) TAUMIN = SMALL
   10 CONTINUE

      RETURN
      END !Tdumin


C**********************************************************************
      SUBROUTINE CHART(XSLOPE,YSLOPE,NAZ,
     &  NZEN,DIFZEN,DEXT,TOTLAI,TAUMIN,
     &  RTA,NLAY,DLAI,EXPDIF)
C This subroutine calculates the number of elemenraty layers in the EHC
C and the transmittance of diffuse radiation through the EHC.
C Subroutine returns:
C NLAY - number of layers - up to 50 (plus one for soil)
C RTA - transmittances of each layer
C DLAI - thickness (in LAI) of the layers
C EXPDIF - transmissivity of EHC layer
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DIFZEN(MAXANG),DEXT(MAXANG)
      REAL RTA(50)

C Initial thickness of canopy layer
C Norman (1979) says this should be always < 0.5 and preferably closer to 0.1.
C Here for LAI < 5, we try DLAI = 0.1
      DLAI = (IFIX(TOTLAI/5.0)+1)*0.1     

C Calculate transmittance through one elementary layer with thickness DLAI
10    KK = 0		! Initialise variable used to calculate NLAY
      WN = 0.0
      WD = 0.0
      DA = TWOPI/FLOAT(NAZ)
      DO 20 I = 1,NZEN
        CZ = COS(DIFZEN(I))
        SZ = SIN(DIFZEN(I))
        DO 20 J = 1,NAZ
          DAZ = DA * (J-0.5)
          SLOPE = ATAN(COS(DAZ)*TAN(XSLOPE)+SIN(DAZ)*TAN(YSLOPE))
          IF ((PID2-DIFZEN(I)).GE.SLOPE) THEN
            EXPDIF = (-DEXT(I)*DLAI/CZ)
            IF (EXPDIF.LT.-180.0) THEN
              EXPDIF = 0.0
            ELSE
              WN = WN + EXP(EXPDIF)*SZ*CZ
            END IF
            WD = WD + SZ*CZ
          END IF
20    CONTINUE
      EXPDIF = WN/WD

C Check to see not too many elementary layers - if so, try again with
C higher DLAI.
30    KK = KK + 1
      IF (KK.GT.50) THEN
        DLAI = DLAI + 0.1
        GO TO 10
      END IF

C Set transmittance for each layer
      IF (KK.EQ.1) THEN			! First layer
        IF (EXPDIF.EQ.0.0) THEN	! Extreme case of no transmittance
          RTA(KK) = 0.0
        ELSE
          RTA(KK) = 1.0			! Normal case: transmittance above first layer = 1
        END IF
      ELSE
        RTA(KK) = RTA(KK-1)*EXPDIF
      END IF

C NLAY is the number of layers such that the smallest transmittance of the
C real canopy is just larger than the transmittance through the EHC.
C When this threshold is reached, the subroutine finishes.
      IF (RTA(KK).GE.TAUMIN) GO TO 30
      NLAY = KK
      RETURN

      END !Chart


C**********************************************************************
      SUBROUTINE ASSIGN(TU,TD,NUMPNT,NLAY,RTA,
     &  LAYER,MLAYER)
C This subroutine matches the grid point to a point in the EHC.
C The ith grid point will correspond to LAYER(I) in an EHC having
C MLAYER(I) layers.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL RTA(50)
      REAL TU(MAXP), TD(MAXP)
      INTEGER LAYER(MAXP), MLAYER(MAXP)

C Z1 = no of layers above the equivalent layer in the EHC
C Z2 = no of layers below the equivalent layer in the EHC
C The point h in the horizontal canopy is found by finding Z1 such that
C tau(Z1) = TD(IPT) and tau(Z2) = TU(IPT). The point is then at depth
C Z1 in a canopy of thickness (Z1 + Z2). (Norman & Welles 1983)

      DO 100 IPT = 1,NUMPNT

        DIFMU = 1.
        DO 200 ILAYER = 1,NLAY
          DIFUP = ABS(TU(IPT)-RTA(ILAYER))
          IF (DIFMU.GE.DIFUP) GO TO 200
          Z1 = ILAYER - 1
          GO TO 290
200     DIFMU = DIFUP
        Z2 = NLAY

290     DIFMD = 1.
        DO 300 ILAYER = 1,NLAY
          DIFDN = ABS(TD(IPT)-RTA(ILAYER))
          IF (DIFMD.GE.DIFDN) GO TO 300
          Z2 = ILAYER - 1
          GO TO 390
300     DIFMD = DIFDN
        Z1 = NLAY

390     MLAYER(IPT) = Z2 + Z1
        LAYER(IPT) = Z1
        IF (LAYER(IPT).EQ.0) LAYER(IPT) = 1

100   CONTINUE    ! End loop over grid points

      RETURN
      END !Assign


C**********************************************************************
      SUBROUTINE SUN(IDOY,ALAT,TTIMD,DEC,EQNTIM,DAYL,SUNSET)
C this is a subroutine to calculate the daylength. it is from
C  a paper by Bruce Barkstrom (1981,'what time does the sun rise
C and set?' Byte.102-114)
C**********************************************************************

      INCLUDE 'maestcom'

      REAL LUNLON,MLON,NUTOBL,NUTLON,MUM,MUN,MUA

C      DOY = REAL(MOD(IDOY*100,36525)/100)
C      DOY = 31047
      DOY = DAYJUL(IDOY)
      T = DOY/36525
C COMPUTE SOLAR ORBIT
      ECC = ECCENT(T)
      RM = ANOM(T,DOY)
      E = RM
      DO 1122 IM = 1,3
 1122 E = E + (RM- (E-ECC*SIN(E)))/ (1-ECC*COS(E))
      V = 2.0*ATAN(SQRT((1+ECC)/ (1-ECC))*TAN(0.5*E))
      IF (V.LT.0.0) V = V + TWOPI
      R = 1 - ECC*COS(E)
      EPS = EPSIL(T)
      OMEG = OMEGA(T,DOY)
C COMPUTE NUTATION TERMS
      LUNLON = ALUNAR(T,DOY)
      NUTOBL = (2.5583E-3+2.5E-7*T)*COS(LUNLON)*PID180
      EPS = EPS + NUTOBL
      NUTLON = - (4.7872E-3+4.7222E-6*T)*SIN(LUNLON)*PID180
C COMPUTE SOLAR DECLINATION
      DEC = ASIN(SIN(EPS)*SIN(V+OMEG))
C COMPUTE EQN OF TIME
      MLON = OMEG + RM
      IF (MLON.LT.0.0) MLON = MLON + TWOPI
      IF (MLON.GT.TWOPI) MLON = MLON - TWOPI*IFIX(MLON/TWOPI)
      Y = TAN(0.5*EPS)
      Y = Y*Y
      Y = (1-Y)/ (1+Y)
      SL = OMEG + NUTLON + V
      IF (SL.LT.0.0) SL = SL + TWOPI
      IF (SL.GT.TWOPI) SL = SL - TWOPI*IFIX(SL/TWOPI)
      AO = ATAN(Y*TAN(SL))
      EQNTIM = AO - MLON
      EQNTIM = EQNTIM - PI*IFIX(EQNTIM/PI)
      IF (ABS(EQNTIM).GT.0.9*PI) EQNTIM = EQNTIM - PI*EQNTIM/ABS(EQNTIM)
      AO = EQNTIM + MLON
      IF (AO.GT.TWOPI) AO = AO - TWOPI*IFIX(AO/TWOPI)
C DAY LENGTH
      MUM = COS(ALAT-DEC)
      MUN = -COS(ALAT+DEC)
      MUA = 0.0
      REFAC = 0.0
      UMN = -MUM*MUN
      IF (UMN.GT.0.0) REFAC = 0.05556/SQRT(UMN)
      IF (MUN.GT.MUA) MUA = MUN
      IF (MUM.LE.MUA) GO TO 1133
      FRACSU = SQRT((MUA-MUN)/ (MUM-MUA))
      FRACSU = 1.0 - 2.0*ATAN(FRACSU)/PI
      SUNSET = HHRS*FRACSU
      SUNRIS = SUNSET
      SUNSET = SUNSET + REFAC + EQNTIM*HHRS/PI
      SUNRIS = SUNRIS + REFAC - EQNTIM*HHRS/PI
      SUNSET = SUNSET + HHRS + TTIMD
      SUNRIS = HHRS - SUNRIS + TTIMD
      EQNTIM = EQNTIM*HHRS/PI
      DAYL = SUNSET - SUNRIS
      RETURN
 1133 STOP
      END !Sun


C**********************************************************************
      FUNCTION DAYJUL(IDATE)
C calculateS JULIAN DAY - NEEDED FOR SUN
C**********************************************************************

      JD = JDATE(IDATE)
      IYEAR = (IDATE - JD)*100/36525 + 1951
      D = REAL(JD)

      NY = IYEAR + 4712
      NLEAP = (NY-1)/4
      DAYJUL = 365*NY + NLEAP + D
      IF (IYEAR.LT.1583) GO TO 10
      DAYJUL = DAYJUL - 10
      DAYJUL = DAYJUL - (IYEAR-1501)/100
      DAYJUL = DAYJUL + (IYEAR-1201)/400
   10 IF (IYEAR.EQ.1582 .AND. D.GE.319.0) DAYJUL = DAYJUL - 10
      IF (4.0* (IYEAR/4).EQ.IYEAR .AND. JD.GE.60) DAYJUL = DAYJUL + 1
      DAYJUL = DAYJUL - 2415020 ! Subtract DAYJUL(1/1/1900)

      RETURN
      END  !DayJul

C**********************************************************************
      FUNCTION ECCENT(T)
C calculateS ECCENTRICITY - USED IN SUN
C**********************************************************************
      ECCENT = 0.01675104 - (4.08E-5+1.26E-7*T)*T
      RETURN
      END  !Eccent

C**********************************************************************
      FUNCTION ANOM(T,D)
C calculateS MEAN ANOMALY IN RADIANS - USED IN SUN
C**********************************************************************

      INCLUDE 'maestcom'

      ANOM = -1.52417 + (1.5E-4+3.0E-6*T)*T*T
      ANOM = ANOM + 0.98560*D
      IF (ANOM.LE.360.0) GO TO 10
      ANOM = ANOM - 360.0*IFIX(ANOM/360.0)
   10 ANOM = ANOM*PID180

      RETURN
      END !Anom


C**********************************************************************
      FUNCTION EPSIL(T)
C calculate OBLIQUITY OF ECLIPTIC - USED IN SUN
C**********************************************************************

      INCLUDE 'maestcom'

      EPSIL = (23.452294- (1.30125E-2+ (1.64E-6-5.03E-7*T)*T)*T)*PID180

      RETURN
      END !Epsil

C**********************************************************************
      FUNCTION OMEGA(T,D)
C calculateS MEAN LONGITUDE OF PERIGEE - USED IN SUN
C**********************************************************************

      INCLUDE 'maestcom'

      OMEGA = (281.22083+ (4.53E-4+3.0E-6*T)*T*T+4.70684E-5*D)*PID180

      RETURN
      END  !Omega


C**********************************************************************
      FUNCTION ALUNAR(T,D)
C COMPUTE LONGITUDE OF ASCENDING NODE OF LUNAR ORBIT  - USED IN SUN
C**********************************************************************

      INCLUDE 'maestcom'

      ALUNAR = (259.1833+ (2.078E-3+2.0E-6*T)*T*T-0.05295*D)*PID180

      RETURN
      END !Alunar

C**********************************************************************
      SUBROUTINE ZENAZ(ALAT,TTIMD,BEAR,DEC,EQNTIM,ZEN,AZ)
C Calculates hourly zenith and azimuth angles.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL ZEN(KHRS),AZ(KHRS)

      DO 10 I = 1,KHRS
        R = I - 0.5
        HI = (R-TTIMD-EQNTIM-HHRS)*PI/HHRS

        ZEN(I) = ACOS(SIN(ALAT)*SIN(DEC)
     &              + COS(ALAT)*COS(DEC)*COS(HI))
        SINAZ = COS(DEC)*SIN(HI)/SIN(ZEN(I))
        IF (SINAZ.GT.1.0) THEN
          SINAZ = 1.0
        ELSE IF (SINAZ.LT.-1.0) THEN
          SINAZ = -1.0
        END IF
        Azimth = ASIN(SINAZ)
        AZ(I) = Azimth-BEAR
        IF (ALAT.LT.0.0) AZ(I)=-AZ(I)
10    CONTINUE

      RETURN
      END


C**********************************************************************
      SUBROUTINE SLOPES(I,TTIMD,EQNTIM,ALAT,DEC,
     &  XSLOPE,YSLOPE,BEAR,ZEN,
     &  BMULT,DMULT2,SOMULT)
C  calculate slope corrections
C  The modifying mulitilipers are got from M.D. Steven and M.H. unsworth
C   (1979,Quart. J. R. Met. Soc. 105,pp 593-602) and M.D. Steven
C      and M. H. Unsworth(1980,Quart.J. R. Met. Soc. 106,pp57-61) with
C      the consideration of reflecting radiation by the slope.
C NB The routine also calculates DMULT & DMULT1 but these aren't used
C anywhere so they are not passed back to the main program. BM May 1997.
C**********************************************************************

      INCLUDE 'maestcom'

C     We assume alebdo of solar radiation of the slope in question is 0.2
C    according to Steven et al, the parameter to describe the
C     intensity of diffuse on the sky is between that of the cloudless
C    and that of the cloudness sky. values they suggested are
      PARMB1 = -0.87
      PARMB2 = 1.23

      IF (XSLOPE.LE.0.0 .AND. YSLOPE.LE.0.0) THEN
         DMULT1 = 1.0
         DMULT2 = 1.0
         SOMULT = 0.0
         SLOPE = 0.00

      ELSE
         IF (ABS(XSLOPE).LT.0.01) GO TO 64
         IF (ABS(YSLOPE).LT.0.01) GO TO 66
         AZSLOP = ATAN(TAN(YSLOPE)/TAN(XSLOPE))
         SLOPE = ATAN(COS(AZSLOP)*TAN(XSLOPE)+SIN(AZSLOP)*TAN(YSLOPE))
         GO TO 70

   64    AZSLOP = PID2
         SLOPE = YSLOPE
         GO TO 70

   66    AZSLOP = 0.0
         SLOPE = XSLOPE
   70    IF (ALAT.LT.0.0) go to 75
C N HEMISPHERE
         IF (SLOPE.LT.0.0) GO TO 73
C    POSITIVE SLOPE
         AZSLOP = PI - AZSLOP + BEAR
         GO TO 80
C  NEGATIVE SLOPE
   73    AZSLOP = TWOPI - AZSLOP + BEAR
         GO TO 80
C S HEMISPHERE
   75    IF (SLOPE.LT.0.0) GO TO 78
C POSITIVE SLOPE
         AZSLOP = TWOPI - AZSLOP - BEAR
         GO TO 80
C NEGATIVE SLOPE
   78    AZSLOP = PI - AZSLOP - BEAR
   80    SLOPE = ABS(SLOPE)
         ALATST = ASIN(COS(SLOPE)*SIN(ALAT)-
     +            SIN(SLOPE)*COS(ALAT)*COS(AZSLOP))
         GH = ASIN(SIN(SLOPE)*SIN(AZSLOP)/COS(ALATST))
         SOMULT = 0.1* (1.-COS(SLOPE))
         DTEMP = (1.0+COS(SLOPE))/2.0
         TEMP = SIN(SLOPE) - SLOPE*COS(SLOPE) -
     +          3.14*0.5* (1.-COS(SLOPE))
         DMULT1 = DTEMP + 2.*PARMB1/ (3.14* (3.+2.*PARMB1))*TEMP
         DMULT2 = DTEMP + 2.*PARMB2/ (3.14* (3.+2.*PARMB2))*TEMP
         DMULT1 = DMULT1*DTEMP
         DMULT2 = DMULT2*DTEMP
      END IF
 90   IF (SLOPE.EQ.0.0000) THEN
          BMULT = 1.0

      ELSE
          Hrang = (I-TTIMD-EQNTIM-HHRS)*PI/HHRS
C Make comparable to TPUMAES.FOR
          CZENST = SIN(ALATST)*SIN(DEC) +
     +             COS(ALATST)*COS(DEC)*COS(Hrang-GH)
          BMULT = CZENST/COS(ZEN)
      END IF
      IF (BMULT.NE.0.00) THEN
          TEMPSP = 1.00
      ELSE
          TEMPSP = 0.00
      END IF
C      DMULT(I) = FSUN(I)*DMULT1 + (1.0-FSUN(I))*DMULT2
   95 CONTINUE
      RETURN
      END  !Slopes


C**********************************************************************
      SUBROUTINE EXBEAM(NALPHA,ALPHA,FALPHA,RANDOM,BZEN,
     &  BEXT,BEXTANG)
C   calculate the extinction coefficients for beam radiation
C**********************************************************************

      INCLUDE 'maestcom'

      REAL FALPHA(MAXANG),ALPHA(MAXANG),BEXTANG(MAXANG)

      SUM = 0.0
      DO 10 IALP = 1,NALPHA
        BEXTANG(IALP) = COSDEL(RANDOM,BZEN,ALPHA(IALP))
        SUM = SUM + BEXTANG(IALP)*FALPHA(IALP)
   10 CONTINUE
      BEXT = SUM

      RETURN
      END  !Exbeam



C**********************************************************************
      SUBROUTINE TRANSB(IHOUR,IPROG,
     &  ZENITH,AZMTH,XSLOPE,YSLOPE,FBEAM,BEXT,
     &  XPT,YPT,ZPT,RX,RY,RZ,DXT,DYT,DZT,
     &  XMAX,YMAX,SHADEHT,
     &  FOLT,ZBC,JLEAF,BPT,NOAGEC,PROP,JSHAPE,SHAPE,NT,
     &  SLA)
C a subroutine to calculate the transmittances of direct radiation
C and also the within and between-tree shading for beam
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DXT(MAXT),DYT(MAXT),DZT(MAXT),RX(MAXT),RY(MAXT),
     +     RZ(MAXT),ZBC(MAXT),FOLT(MAXT)
      REAL BPT(8,MAXC),PROP(MAXC)
      REAL FBEAM(KHRS,3)
      LOGICAL SHADED

C  Set flag so that only upper hemisphere will be considered.
      IFLAG = 1
C Default value
      SLA = 0.0

      SLOPE = ATAN(COS(AZMTH)*TAN(XSLOPE)+SIN(AZMTH)*TAN(YSLOPE))
      IF ((PID2-ZENITH).LT.SLOPE) RETURN
      IF (SHADED(XPT,YPT,ZPT,ZENITH,AZMTH,XMAX,YMAX,SHADEHT)) RETURN

      IF ((FBEAM(IHOUR,1).GT.0.0).OR.(FBEAM(IHOUR,2).GT.0.00)) THEN
        CALL TREDST(IFLAG,IPROG,ZENITH,AZMTH,
     &      XPT,YPT,ZPT,RX,RY,RZ,DXT,DYT,DZT,
     &      FOLT,ZBC,JLEAF,BPT,NOAGEC,PROP,JSHAPE,SHAPE,
     &      NT,S,S1)
C       BPATH = the total weighted beam pathlength
        BPATH = S + S1
C       SHU1= the weighted pathlength in the target tree.
        SHU1 = S
C       SLA is the sunlit leaf area associated with the grid point.
        BTEMP = -BEXT*BPATH
        IF (BTEMP.LT.-180.0) THEN
          SLA = 0.000000
        ELSE
          SLA = EXP(BTEMP)
        END IF
      endif

      RETURN
      END !Transb


C**********************************************************************
      SUBROUTINE SCATTER(IPT,IWAVE,
     &  MLAYERI,LAYERI,DLAI,EXPDIF,ZEN,BEXT,
     &  DMULT2,SOMULT,BMULT,
     &  RADABV,FBEAM,TAIR,TSOIL,
     &  ARHO,ATAU,RHOSOL,
     &  DIFUP,DIFDN,SCLOST)
C This subroutine calculates the scattered radiation using the iterative
C method of Norman (1979). Outputs are
C DIFUP: the upwards scattered flux below gridpoint IPT
C DIFDN: the downwards scattered flux above gridpoint IPT
C SCLOST: the scattered radiation lost from the top of the canopy
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DIFDN(MAXP,3),DIFUP(MAXP,3),SCLOST(MAXP,3)

      DIMENSION U(3,101),D(3,101),SUP(101),SDN(101),ADUM(101),TBEAM(101)
C These are temporary arrays within this subroutine.
C U(K,I) is the relative upwards flux in wavelength k incident on layer I
C D(K,I) is the relative downwards flux in wavelength k incident on layer I
C SUP(I) is the beam flux scattered upwards from layer I
C SDN(I) is the beam flux scattered downwards from layer I
C ADUM(I) is the relative reflected flux upwards from layer I
C TBEAM(I) is the relative beam flux density at layer I

      IF (IWAVE.EQ.3) THEN	! Call separate subroutine for THERMAL radiation
	  CALL ABSTHERM(IPT,MLAYERI,LAYERI,EXPDIF,
     &  RADABV,TAIR,TSOIL,RHOSOL,
     &  DIFUP,DIFDN,SCLOST)
	  GOTO 1300			! End of this subroutine
      END IF

C Jtot is the no. of layers in EHC (+1 for soil layer)
C NB the top layer is JTOT, the bottom layer is 1.
      JTOT = MLAYERI

C Calculate the transmittance of beam radiation for an elementary layer.
      COSZEN = COS(ZEN)
      EXPDIR = EXP(-BEXT*DLAI/COSZEN)

C calculate LAYER PROPERTIES
      RLAYER = (1.-EXPDIF)*ARHO ! prob. of radiation being scattered upwards
      TLAYER = (1.-EXPDIF)*ATAU + EXPDIF ! prob. of radiation penetrating downwards
      TBEAM(JTOT) = FBEAM ! relative beam flux density at the top layer
      SDN(1) = 0. ! downwards flux density reaching the soil surface

C The scattered beam flux density in each elementary layer is calculated
C with the assumption of no flux reflected from the soil surface.
      DO 600 JJ = JTOT-1, 1, -1
	  TBEAM(JJ) = 0.0
        IF (TBEAM(JJ+1).GT.1E-9.AND.EXPDIR.NE.0.0) THEN
          TEMPS = ALOG10(TBEAM(JJ+1)) + ALOG10(EXPDIR)
          IF (TEMPS.GE.-45.0) TBEAM(JJ) = TBEAM(JJ+1)*EXPDIR
	  END IF
        SUP(JJ+1) = (TBEAM(JJ+1)-TBEAM(JJ))*ARHO
        SDN(JJ+1) = (TBEAM(JJ+1)-TBEAM(JJ))*ATAU
600   CONTINUE

C ADUM(1) is the relative reflected upwards flux from the soil surface.
      ADUM(1) = RHOSOL*(DMULT2+SOMULT)
C ADUM(J) is the ratio of upwards to downwards diffuse flux (eqn 3.6.13, Norman 1979)
      TLAY2 = TLAYER*TLAYER
      DO 700 J = 2,JTOT
        ADUM(J) = ADUM(J-1)*TLAY2/ (1.-ADUM(J-1)*RLAYER) + RLAYER
700   CONTINUE

C Calculate boundary conditions for iteration.
C SUP(1) is the reflected beam flux from the soil surface
      SUP(1) = TBEAM(1)*RHOSOL*BMULT
C D(IWAVE,JTOT) is the relative incident diffuse flux above the top layer.
      D(IWAVE,JTOT) = 1. - FBEAM

C Calculate initial solution (cf. eq 3.6.14, Norman 1979)
      DO 800 JJ = JTOT-1, 1, -1
        D(IWAVE,JJ) = D(IWAVE,JJ+1)*TLAYER/(1.-ADUM(JJ)*RLAYER)
     &				+ SDN(JJ+1)
        U(IWAVE,JJ+1) = ADUM(JJ+1)*D(IWAVE,JJ+1) + SUP(JJ+1)
800   CONTINUE
      U(IWAVE,1) = RHOSOL*(D(IWAVE,1)*DMULT2 +
     +             (D(IWAVE,1)+TBEAM(1))*SOMULT) + SUP(1)

C Recursive method to calculate scattered flux density at IPT until
C equilibrium conditions met:
C ABS(DOWN - D(IWAVE,JJ)) < 0.01, ABS(UP - U(IWAVE,JJ)) < 0.01 for all JJ

      ITER = 0
900   IREPT = 0
      ITER = ITER + 1
      DO 1000 J = 2,JTOT
         JJ = JTOT - J + 1
         JJP1 = JJ + 1
         DOWN = TLAYER*D(IWAVE,JJP1) + U(IWAVE,JJ)*RLAYER + SDN(JJP1)
         IF (ABS(DOWN-D(IWAVE,JJ)).GT.0.01) IREPT = 1
1000  D(IWAVE,JJ) = DOWN

      U(IWAVE,1) = (D(IWAVE,1)*DMULT2+ (D(IWAVE,1)+TBEAM(1))*SOMULT+
     +         TBEAM(1)*BMULT)*RHOSOL
      DO 1200 JJ = 2,JTOT
         JJM1 = JJ - 1
         UP = RLAYER*D(IWAVE,JJ) + U(IWAVE,JJM1)*TLAYER + SUP(JJ)
         IF (ABS(UP-U(IWAVE,JJ)).GT.0.01) IREPT = 1 
         U(IWAVE,JJ) = UP
1200  CONTINUE
      IF (IREPT.NE.0) GOTO 900

C Summarise calculations.
C This is a mix of code from previous versions of Maestro. Think it's right!
      IF (EXPDIF.LT.0.0000001) THEN
        SKYDIF = 0.00
      ELSE
        SKYDIF = (EXPDIF**(JTOT-LAYERI))*(1.-FBEAM)
      END IF
      DIFDN(IPT,IWAVE) = (D(IWAVE,LAYERI)-SKYDIF)*RADABV
      DIFUP(IPT,IWAVE) = U(IWAVE,LAYERI)*RADABV

1300  RETURN
      END  !Scatter


C**********************************************************************
      SUBROUTINE ABSTHERM(IPT,MLAYERI,LAYERI,EXPDIF,
     &  RADABV,TAIR,TSOIL,RHOSOL,
     &  DIFUP,DIFDN,SCLOST)
C Calculate long-wave radiation absorption or emission. After Norman (1979). 
C Note that foliage temperature is assumed equal to air temperature. 
C This is appropriate here because the calculated value is part of the Rn
C used in the Penman-Monteith equation, where the leaf is assumed at air T. 
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DIFDN(MAXP,3),DIFUP(MAXP,3),SCLOST(MAXP,3)

      DIMENSION U(101),D(101)
C These are temporary arrays within this subroutine.
C U(I) is the upwards long-wave flux incident on layer I
C D(I) is the  downwards long-wave flux incident on layer I

C Jtot is the no. of layers in EHC (+1 for soil layer)
C NB the top layer is JTOT, the bottom layer is 1.
 	JTOT = MLAYERI	  ! Total number of layers

C Calculate the values of relative downwards and upwards flux densities
C for thermal radiation for the elementary layers.
      D(JTOT) = RADABV ! Thermal flux above the top layer in EHC.
	ELEAF = SIGMA * (TK(TAIR)**4.) 	! Emission from leaf (assume Tleaf = Tair)
      FF = ELEAF * (1.-EXPDIF) ! Relative emissivity of elementary layer
      DO 10 J = JTOT-1, 1, -1
         D(J) = D(J+1)*EXPDIF + FF
10    CONTINUE

      ESOIL=(1.0-RHOSOL)*SIGMA*(TK(TSOIL)**4) ! Upwards flux from soil
      U(1) = ESOIL + RHOSOL*D(1)
C Error corrected from previous version, which here had D(JTOT-1) in lieu of D(1). V. small effect. 
      DO 20 J = 2,LAYERI
         U(J) = U(J-1)*EXPDIF + FF
20    CONTINUE

C Summarise calculations.
      DIFDN(IPT,3) = D(LAYERI) - ELEAF*EMLEAF	! Downward thermal flux
      DIFUP(IPT,3) = U(LAYERI) - ELEAF*EMLEAF	! Upward thermal flux
      SCLOST(IPT,3) = 0.0

	RETURN
	END ! Abstherm


C**********************************************************************
      subroutine ABSRAD(
     &  IPT,IWAVE,
     &  NZEN,DEXT,BEXT,BMULT,RELDF,
     &  RADABV,FBEAM,ZEN,ABSRP,DIFDN,DIFUP,
     &  DFLUX,BFLUX,SCATFX)
C This subroutine ends all radiation calculations and summarises
C all calculated results which are to be output.
C Outputs: BFLUX = beam radiation flux absorbed at IPT
C   DFLUX = diffuse radiation flux absorbed at IPT (incl. scattered)
C   SCATFX = amount of diffuse radiation which comes from scattering
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DEXT(MAXANG)
      REAL DFLUX(MAXP,3), BFLUX(MAXP,3), SCATFX(MAXP,3)

C Calculate beam radiation absorbed by sunlit leaf area.

      IF ((IWAVE.EQ.3).OR.(BMULT.EQ.0.0)) THEN
        BFLUX(IPT,IWAVE) = 0.0
      ELSE
        BFLUX(IPT,IWAVE) = RADABV*FBEAM/COS(ZEN)*ABSRP
      END IF

C Calculate diffuse radiation absorbed directly.

      DFX = RELDF*RADABV*(1.0-FBEAM)

C Calculate radiation absorbed from scattering.
      SCAT = DIFDN + DIFUP

C Calculate mean diffuse extinction coefficient
      DMEAN = 0.0
      do 10 IZEN=1,NZEN
         DMEAN = DMEAN + DEXT(IZEN)
10    continue
      DMEAN = DMEAN/REAL(NZEN)

C BM DFX does not need to be *DMEAN because RELDF takes this into account 
      SCATFX(IPT,IWAVE) = SCAT*DMEAN

      IF(IWAVE.EQ.3) THEN
         DFLUX(IPT,IWAVE) = SCATFX(IPT,IWAVE)
      ELSE
         DFLUX(IPT,IWAVE) = DFX + SCATFX(IPT,IWAVE)
      ENDIF
      DFLUX(IPT,IWAVE)= DFLUX(IPT,IWAVE)*ABSRP
      SCATFX(IPT,IWAVE)= SCATFX(IPT,IWAVE)*ABSRP
      RETURN
      END !Absrad


C**********************************************************************
      SUBROUTINE CATEGO(AREA,APAR,HISTO,BINSIZE)
C This function creates a histogram of PAR frequency in the canopy 
C over the simulation period. Frequency in m2.HR.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL HISTO(MAXHISTO)

      IHISTO = NINT((APAR+BINSIZE/2.)/BINSIZE)
      IF (IHISTO.LE.MAXHISTO) THEN
        HISTO(IHISTO) = HISTO(IHISTO) + AREA
      ELSE
        CALL SUBERROR('HISTOGRAM OVERFLOW',IWARN,0)
      END IF

      RETURN
      END !Catego

