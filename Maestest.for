
C =======================================================
C | MAESTRO: ECOCRAFT version:  Belinda Medlyn, 1997    |
C | Based on: modified version: Ying-Ping Wang, 1996    |
C =======================================================

C------------------------------------------------------------------------
C This main program tests the radiation interception routines.
C The positions of radiation sensors are specified and the program
C calculates radiation transmission to those positions. 
C This file contains the main program.
C All other code is contained in the additional files: -
C   radn.for - calculation of radiation interception
C   physiol.for - physiology subroutines
C   getmet.for - read in met data
C   inout.for - handle input & output data
C   utils.for - utility subroutines.
C------------------------------------------------------------------------

      PROGRAM MAESTEST

C------------------------------------------------------------------------
C  For more information write to Prof P.G.Jarvis,Institute of Ecology and
C  Resource Management, Darwin Building, Mayfield Rd, Edinburgh, Scotland.
C------------------------------------------------------------------------

      INCLUDE 'maestcom'

C Array declarations.
C List of trees for which to do calculations
      INTEGER ITARGETS(MAXT)
C Tree positions and dimensions - all trees, all dates
      REAL DXT1(MAXT),DYT1(MAXT),DZT1(MAXT)
      REAL RXTABLE1(MAXTDATE,MAXT),RYTABLE1(MAXTDATE,MAXT)
      REAL RZTABLE1(MAXTDATE,MAXT),ZBCTABLE1(MAXTDATE,MAXT)
      REAL FOLTABLE1(MAXTDATE,MAXT),DIAMTABLE1(MAXTDATE,MAXT)
C Tree positions and dimensions - sorted trees, all dates
      REAL RXTABLE(MAXTDATE,MAXTT),RYTABLE(MAXTDATE,MAXTT)
      REAL RZTABLE(MAXTDATE,MAXTT),ZBCTABLE(MAXTDATE,MAXTT)
      REAL FOLTABLE(MAXTDATE,MAXTT),TOTLAITABLE(MAXTDATE)
      REAL DIAMTABLE(MAXTDATE,MAXTT)
C Dates for tree dimensions
      INTEGER DATESX(MAXTDATE),DATESY(MAXTDATE),DATESZ(MAXTDATE)
      INTEGER DATEST(MAXTDATE),DATESLA(MAXTDATE),DATESD(MAXTDATE)
C Tree dimensions on simulation date (by interpolation)
      REAL DXT(MAXTT),DYT(MAXTT),DZT(MAXTT)
      REAL RX(MAXTT),RY(MAXTT),RZ(MAXTT),ZBC(MAXTT)
      REAL FOLT(MAXTT),DIAM(MAXTT)
C Positions of grid points, associated volume & leaf area, etc
      REAL XL(MAXP),YL(MAXP),ZL(MAXP),VL(MAXP),DLT(MAXP),DLI(MAXC,MAXP)
      INTEGER LGP(MAXP),LAYER(MAXP),MLAYER(MAXP)
C Met data
      INTEGER METCOLS(MAXMET)
      REAL windah(KHRS),tsoil(KHRS),tair(KHRS),radabv(KHRS,3)
      REAL FBEAM(KHRS,3),RH(KHRS),VPD(KHRS),VMFD(KHRS)
      REAL CA(KHRS),PRESS(KHRS),PPT(KHRS),SOILMD(KHRS)
	REAL DELTAT(12)
C Physiology inputs by layer
      REAL ABSRP(MAXLAY,3),ARHO(MAXLAY,3),ATAU(MAXLAY,3)
      REAL RHOSOL(3)
C Structural data inputs
      REAL BPT(8,MAXC),PROP(MAXC)
      REAL ALPHA(MAXANG),FALPHA(MAXANG)
C Intermediate calculations
      REAL DIFZEN(MAXANG),DEXT(MAXANG),BEXTANG(MAXANG)
      REAL ZEN(KHRS),AZ(KHRS)
      REAL TU(MAXP),TD(MAXP),RELDF(MAXP)
      REAL DIFDN(MAXP,3),DIFUP(MAXP,3),SCLOST(MAXP,3)
      REAL BFLUX(MAXP,3),DFLUX(MAXP,3),SCATFX(MAXP,3)
C Titles of input files
      CHARACTER*80 CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE

C Set program flag
      IPROG = ITEST
	VTITLE = 'MAESTEST: Version September 2001'

C Open input files
      CALL OpenInputf(IODAILY,IOHRLY,IOTUTD,IOHIST,IORESP,IOSOIL,
     &  CTITLE,TTITLE,PTITLE,STITLE,UTITLE)

C Get input from control file
      CALL INPUTCON(
     &  ISTART, IEND, NSTEP,
     &  NUMPNT, NOLAY, NZEN, DIFZEN, NAZ,
     &  MODELGS, MODELJM, MODELRD, MODELSS, MODELRW, ITERMAX,
     &  IOHIST, BINSIZE,
     &  ICC, CO2INC, TINC,
     &  IOTC, TOTC, WINDOTC, PAROTC, FBEAMOTC
     &  )

C Get input from canopy structure file
      CALL InputStr(
     &  JLEAF,BPT,RANDOM,NOAGEC,
     &  JSHAPE,SHAPE,EXTWIND,
     &  NALPHA,ALPHA,FALPHA,
     &  COEFFT,EXPONT,WINTERC,
     &  BCOEFFT,BEXPONT,BINTERC,
     &  RCOEFFT,REXPONT,RINTERC,FRFRAC,
     &  PFLA,CANCAPLA
     &  )

C Get input from physiology file
C For the test program, only need reflectance & transmittance coeffts
      CALL READAGEP(UPHY, NOAGEC, NOAGEP)
      CALL READPROP(UPHY, NOAGEC, NOAGEP, PROP)
      CALL READABSRP(UPHY,NOLAY,ABSRP,ARHO,ATAU,RHOSOL)

C Get input from trees file
      CALL INPUTTREE(
     &  XSLOPE,YSLOPE,BEAR,XMAX,YMAX,STOCKING,
     &  ZHT,Z0HT,ZPD,
     &  NOALLTREES,NOTREES,NOTARGETS,ITARGETS,SHADEHT,
     &  NOXDATES,NOYDATES,NOZDATES,NOTDATES,NOLADATES,NODDATES,
     &  DATESX,DATESY,DATESZ,DATEST,DATESLA,DATESD,
     &  DXT1,DYT1,DZT1,RXTABLE1,RYTABLE1,RZTABLE1,ZBCTABLE1,
     &  FOLTABLE1,TOTLAITABLE,DIAMTABLE1,
     &  IFLUSH,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN
     &  )

C Open met data file (must be done after ISTART & IEND read)
      CALL OPENMETF(ISTART,IEND,CAK,PRESSK,SWMIN,SWMAX,
     &  DIFSKY,ALAT,TTIMD,DELTAT,
     &  MFLAG,METCOLS,NOMETCOLS,MTITLE,MSTART)

C Open files and read information about points
      CALL GETPOINTSF(NUMPNT,XL,YL,ZL,
     &	CTITLE,TTITLE,MTITLE,STITLE,VTITLE)

C Do calculations which are not day-dependent
      CALL EXDIFF(NALPHA,ALPHA,FALPHA,NZEN,DIFZEN,RANDOM,DEXT)

C Sort trees in order of distance from first sensor position
      CALL SORTTREESP(
     &  XL(1),YL(1),NOALLTREES,NOTREES,
     &  DXT1,DYT1,DZT1,RXTABLE1,RYTABLE1,RZTABLE1,
     &  ZBCTABLE1,FOLTABLE1,DIAMTABLE1,
     &  DXT,DYT,DZT,RXTABLE,RYTABLE,RZTABLE,
     &  FOLTABLE,ZBCTABLE,DIAMTABLE
     &  )

C Start met file
      CALL RESTARTMETF(ISTART,MSTART,MFLAG)

C**********************************************************************
C Begin daily loop
      IDAY = 0
10    CONTINUE
      WRITE (*,105) IDAY
105   FORMAT ('  DAY:',I5)

C Calculate zenith angle of sun
      CALL SUN(IDAY+ISTART,ALAT,TTIMD,DEC,EQNTIM,DAYL,SUNSET)
      CALL ZENAZ(ALAT,TTIMD,BEAR,DEC,EQNTIM,ZEN,AZ)

C Get meteorological data
      CALL GETMET(IDAY+ISTART,MFLAG,ZEN,METCOLS,NOMETCOLS,
     &  CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,DAYL,
     &  windah,tsoil,tair,radabv,fbeam,RH,VPD,VMFD,CA,PRESS,PPT,SOILMD)
      IF (ICC.EQ.0) CALL ALTERMETCC(CA,TAIR,TSOIL,RH,VPD,VMFD,PRESS,
     &  CO2INC,TINC)
      IF (IOTC.EQ.0) CALL ALTERMETOTC(TOTC,WINDOTC,PAROTC,FBEAMOTC,
     &  TAIR,TSOIL,WINDAH,RADABV,FBEAM,RH,VPD,VMFD,PRESS)

C Interpolate to get daily values of parameters - structural only
      CALL INTERPOLATET(IDAY,ISTART,
     &  NOXDATES,DATESX,RXTABLE,NOYDATES,DATESY,RYTABLE,
     &  NOZDATES,DATESZ,RZTABLE,NOTDATES,DATEST,ZBCTABLE,
     &  NODDATES,DATESD,DIAMTABLE,
     &  NOLADATES,DATESLA,FOLTABLE,TOTLAITABLE,NOTREES,
     &  RX,RY,RZ,ZBC,FOLT,TOTLAI,DIAM,STOCKING,
     &  IFLUSH,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN,
     &  NEWCANOPY)

C Calculate diffuse transmittances
      CALL TRANSD(
     &  IDAY,IOTUTD,NEWCANOPY,IPROG,NOTREES,XSLOPE,YSLOPE,
     &  NZEN,DIFZEN,NAZ,NUMPNT,DEXT,DIFSKY,
     &  XL,YL,ZL,RX,RY,RZ,DXT,DYT,DZT,
     &  XMAX,YMAX,SHADEHT,
     &  FOLT,ZBC,JLEAF,BPT,NOAGEC,PROP,JSHAPE,SHAPE,
     &  NEWTUTD,TU,TD,RELDF
     &  )

      IF (NEWTUTD.EQ.1) THEN
C If the diffuse transmittances have changed, must set up the EHC
        CALL EHC(NUMPNT,TU,TD,
     &    TOTLAI,XSLOPE,YSLOPE,NAZ,NZEN,DIFZEN,DEXT,
     &    DLAI,EXPDIF,LAYER,MLAYER
     &  )
      END IF

C**********************************************************************
C Begin hourly loop
      DO 100 IHOUR = 1,KHRS

C Test to see if daylight hours or if any foliage
        IF ((ABS(ZEN(IHOUR)).LE.1.57).AND.
     &    (RADABV(IHOUR,1).GT.0.1).AND.
     &    (FOLT(1).GT.0.0)) THEN

C Comments for this hour
        PAR = RADABV(IHOUR,1)*UMOLPERJ
        WRITE (UPOINTSO,980) IHOUR,PAR
980     FORMAT ('HOUR: ',I2,'  INCIDENT PAR: ',F8.2,' UMOL M-2 S-1')

C Get slope correction factor
          call SLOPES(IHOUR,TTIMD,EQNTIM,ALAT,DEC,
     &      XSLOPE,YSLOPE,BEAR,ZEN(IHOUR),
     &      BMULT,DMULT2,SOMULT)

C Get beam extinction coefficient
          CALL EXBEAM(NALPHA,ALPHA,FALPHA,RANDOM,ZEN(IHOUR),
     &      BEXT,BEXTANG)

C Loop over grid points
          DO 80 IPT = 1,NUMPNT

C Calculate the weighted pathlengths for beam radiation.
            CALL TRANSB(IHOUR,IPROG,
     &        ZEN(IHOUR),AZ(IHOUR),XSLOPE,YSLOPE,FBEAM,BEXT,
     &        XL(IPT),YL(IPT),ZL(IPT),RX,RY,RZ,DXT,DYT,DZT,
     &        XMAX,YMAX,SHADEHT,
     &        FOLT,ZBC,JLEAF,BPT,NOAGEC,PROP,JSHAPE,SHAPE,NOTREES,
     &        SUNLA)

C Just use PAR wavelength
            IWAVE = 1

C Calculate the scattered radiation
            CALL SCATTER(IPT,IWAVE,
     &          MLAYER(IPT),LAYER(IPT),DLAI,EXPDIF,
     &          ZEN(IHOUR),BEXT,
     &          DMULT2,SOMULT,BMULT,
     &          RADABV(IHOUR,IWAVE),FBEAM(IHOUR,IWAVE),
     &          TAIR(IHOUR),TSOIL(IHOUR),
     &          ARHO(1,IWAVE),ATAU(1,IWAVE),RHOSOL(IWAVE),
     &          DIFUP,DIFDN,SCLOST)

C Output transmittances
      TBEAM = FBEAM(IHOUR,IWAVE)*PAR*SUNLA
      TDIFF = (1-FBEAM(IHOUR,IWAVE))*PAR*TD(IPT)
      TSCAT = DIFDN(IPT,IWAVE)
      TTOT = TBEAM + TDIFF + TSCAT
      WRITE (UPOINTSO,500) IPT,TBEAM,TDIFF,TSCAT,TTOT  
500   FORMAT(I3,1X,4(F12.5,1X))

90          CONTINUE ! End loop over wavelengths

80        CONTINUE ! End loop over grid points

        END IF ! If day or night

C Output hourly totals

C**********************************************************************
100   CONTINUE
C End hourly loop

C Output daily totals

C Go to next day of simulation
      IDAY = IDAY + NSTEP
      IF ((ISTART+IDAY).LE.IEND) THEN
        CALL SKIPMET(MFLAG,NSTEP)
        GOTO 10
      END IF

C**********************************************************************
110   CONTINUE
C End daily loop

C Write diffuse transmittances to file
      REWIND (UTUTD)
      WRITE (UTUTD,1313) (IPT,TU(IPT),TD(IPT),RELDF(IPT),IPT=1,NUMPNT)
1313  FORMAT (1X,I3,5X,F8.4,1X,F8.4,1X,F8.4)

      CALL Closef(IODAILY,IOHRLY,IOHIST,IORESP,IOSOIL)

      STOP
      END !MAIN
C**********************************************************************


C**********************************************************************
      SUBROUTINE GETPOINTSF(NUMPNT,XL,YL,ZL,
     &  CTITLE,TTITLE,MTITLE,STITLE,VTITLE)
c Subroutine for testing radiation interception routines.
C Open input & output files and read information about sensor positions.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      CHARACTER*80 CTITLE, TTITLE, PTITLE, STITLE, MTITLE, VTITLE
      REAL XL(MAXP),YL(MAXP),ZL(MAXP),COORDS(MAXP*3)

      NAMELIST /CONTROL/ NOPOINTS,INPUTTYPE
      NAMELIST /XYZ/ COORDS
      NAMELIST /TRANSECT/ ANGLE,SPACING,ZHEIGHT

C Open input file
      OPEN (UPOINTSI, FILE = 'points.dat', STATUS = 'OLD')
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('ERROR: POINTS INPUT FILE DOES NOT EXIST',
     &  IFATAL,IOERROR)

C Read title from input file
990   FORMAT (A60)     ! For reading titles in input files.
      READ (UPOINTSI, 990) PTITLE

C Default values
      NOPOINTS = 0
      INPUTTYPE = 1

C Read control flags: no of points and type of input
      READ (UPOINTSI, CONTROL, IOSTAT = IOERROR)
      IF ((IOERROR.NE.0).OR.(NOPOINTS.EQ.0))
     &  CALL SUBERROR('ERROR: MISSING CONTROL INFO IN POINTS FILE',
     &  IFATAL,IOERROR)
      IF (NOPOINTS.GT.MAXP) THEN
        CALL SUBERROR('WARNING: TOO MANY TEST POINTS SPECIFIED',
     &  IWARN,IOERROR)
        NUMPNT = MAXP
      ELSE 
        NUMPNT = NOPOINTS
      END IF

C Read in list of points
      IF (INPUTTYPE.EQ.1) THEN
        READ (UPOINTSI, XYZ, IOSTAT = IOERROR)
        IF (IOERROR.NE.0)
     &  CALL SUBERROR('ERROR READING GRID POINTS',
     &  IFATAL,IOERROR)
        DO 100 N = 1,NUMPNT
          XL(N) = COORDS((N-1)*3 + 1)
          YL(N) = COORDS((N-1)*3 + 2)
          ZL(N) = COORDS(N*3)
100     CONTINUE

C Read in details of transect & construct points
      ELSE IF (INPUTTYPE.EQ.2) THEN
        READ (UPOINTSI, TRANSECT, IOSTAT = IOERROR)
        IF (IOERROR.NE.0)
     &    CALL SUBERROR('ERROR READING TRANSECT DETAILS',
     &    IFATAL,IOERROR)
        ANGLE = ANGLE*PID180
        COSANG = COS(ANGLE)
        SINANG = SIN(ANGLE)
        DIST = SPACING/2.0
        DO 200 N = 1,NUMPNT
          XL(N) = DIST*COSANG
          YL(N) = DIST*SINANG
          ZL(N) = ZHEIGHT
          DIST =  DIST + SPACING
200     CONTINUE
      END IF
  
C Open output file
      OPEN (UPOINTSO, FILE = 'testflx.dat', STATUS = 'UNKNOWN')
C Write headings to output file
991   FORMAT (A12,A60) ! For writing comments to output files.
992   FORMAT (A3,4(A12,1X))
      WRITE (UPOINTSO, 991) 'Program:    ', VTITLE
      WRITE (UPOINTSO, 991) 'Control:    ', CTITLE
      WRITE (UPOINTSO, 991) 'Trees:      ', TTITLE
      WRITE (UPOINTSO, 991) 'Structure:  ', STITLE
      WRITE (UPOINTSO, 991) 'Points:     ', PTITLE
      WRITE (UPOINTSO, 991) 'Met data:   ', MTITLE
      WRITE (UPOINTSO, *)
      WRITE (UPOINTSO, 992) 'PT','BEAM','DIFF','SCAT','TOTAL'
      WRITE (UPOINTSO, *)

      RETURN
      END ! GetPointsF

