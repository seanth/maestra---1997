
C =======================================================
C | MAESTRO: ECOCRAFT version:  Belinda Medlyn, 1997    |
C | Based on: modified version: Ying-Ping Wang, 1996    |
C =======================================================

C------------------------------------------------------------------------
C This file contains the main program.
C It also contains subroutines to zero arrays (ZEROHR, ZEROHIST, ZEROD)
C and to sum arrays (SUMHR, SUMDAILY).
C All other code is contained in the additional files: -
C   radn.for - calculation of radiation interception
C   physiol.for - physiology subroutines
C   getmet.for - read in met data
C   inout.for - handle input & output data
C   utils.for - utility subroutines.
C------------------------------------------------------------------------

      PROGRAM MAESTRO

C-----------------------------------------------------------------------
C     An array model to predict radiation regime inside a plant
C  community. Especially suited to forest stands where crown shape
C  is described using a cone or elipsoid.
C     This program consists of three main parts, i.e.
C  1: radiation in the atmosphere , daylength and the position of
C     the sun at different hours of the day;
C  2: plant canopy architecture.3-D leaf area distribution and random
C     and nonrandom leaf area dispersion inside crown has been incorporated
C     into the radiation model with the spherical leaf angular distribution
C     assumption;radiation penetration and scattering processes.
C     A multi-scattering model proposed by John Norman is used.
C  3: both photosynthesis and transpiration are calculated using the
C     combined gs-A-E-H model.
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
      REAL FOLLAY(MAXLAY),WINDLAY(MAXLAY)
C Met data
      INTEGER METCOLS(MAXMET)
      REAL windah(KHRS),tsoil(KHRS),tair(KHRS),radabv(KHRS,3)
      REAL FBEAM(KHRS,3),RH(KHRS),VPD(KHRS),VMFD(KHRS)
      REAL CA(KHRS),PRESS(KHRS),PPT(KHRS),SOILMD(KHRS)
C Physiology inputs by layer
      REAL ABSRP(MAXLAY,3),ARHO(MAXLAY,3),ATAU(MAXLAY,3)
      REAL RHOSOL(3)
      REAL JMAXTABLE(MAXPDATE,MAXLAY,MAXC)
      REAL VCMAXTABLE(MAXPDATE,MAXLAY,MAXC)
      REAL RDTABLE(MAXPDATE,MAXLAY,MAXC)
      REAL SLATABLE(MAXPDATE,MAXLAY,MAXC)
	REAL AJQTABLE(MAXPDATE,MAXLAY,MAXC)
      INTEGER DATESJ(MAXPDATE), DATESV(MAXPDATE)
      INTEGER DATESRD(MAXPDATE), DATESSLA(MAXPDATE)
	INTEGER DATESA(MAXPDATE)
      REAL JMAX25(MAXLAY,MAXC),VCMAX25(MAXLAY,MAXC)
      REAL RD0(MAXLAY,MAXC),SLA(MAXLAY,MAXC),AJQ(MAXLAY,MAXC)
C Structural data inputs
      REAL BPT(8,MAXC),PROP(MAXC)
      REAL ALPHA(MAXANG),FALPHA(MAXANG)
C Intermediate calculations
      REAL DIFZEN(MAXANG),DEXT(MAXANG),BEXTANG(MAXANG)
      REAL ZEN(KHRS),AZ(KHRS)
      REAL TU(MAXP),TD(MAXP),RELDF(MAXP)
      REAL DIFDN(MAXP,3),DIFUP(MAXP,3),SCLOST(MAXP,3)
      REAL BFLUX(MAXP,3),DFLUX(MAXP,3),SCATFX(MAXP,3)
C Outputs for each tree
      REAL THRAB(KHRS,3),TDYAB(3),TCAN(KHRS)
      REAL FCO2(KHRS),FH2OT(KHRS),FH2OE(KHRS)
	REAL GSCAN(KHRS),FHEAT(KHRS),FH2OCAN(KHRS)
      REAL FRESPF(KHRS),FRESPW(KHRS),FRESPB(KHRS)
      REAL FRESPFR(KHRS),FRESPCR(KHRS)
      REAL PPAR(MAXLAY,KHRS),PPS(MAXLAY,KHRS),PTRANSP(MAXLAY,KHRS)
C Soil characteristics ! for water balance - unfinished
C	REAL SOILWATER(MAXS),SOILCAP(MAXS),SOILDEPTH(MAXS)
C Titles of input files
      CHARACTER*80 CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE
C Outputs for PAR histogram
      REAL HISTO(MAXHISTO)

C Set program flag
      IPROG = INORMAL
	VTITLE = 'MAESTRA: Version September 2001'

C Open input files
      CALL OpenInputf(IODAILY,IOHRLY,IOTUTD,IOHIST,IORESP,IOSOIL,
     &  CTITLE,TTITLE,PTITLE,STITLE,UTITLE)
C      OPEN (100, FILE = 'tmp.out', STATUS = 'UNKNOWN')	! for water balance - unfinished

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
      CALL InputPhy(
     &  MODELJM,MODELRD,MODELGS,MODELRW,
     &  NOLAY,NOAGEC,NOAGEP,PROP,
     &  ABSRP,ARHO,ATAU,RHOSOL,
     &  JMAXTABLE,DATESJ,NOJDATES,IECO,EAVJ,EDVJ,DELSJ,THETA,
     &  VCMAXTABLE,DATESV,NOVDATES,EAVC,EDVC,DELSC,TVJUP,TVJDN,
     &  SLATABLE,DATESSLA,NOSLADATES,NOADATES,DATESA,AJQTABLE,
     &  RDTABLE,DATESRD,NORDATES,Q10F,RTEMP,DAYRESP,TBELOW,
     &  EFFYRW,RMW,RTEMPW,Q10W,COLLA,COLLK,STEMSDW,RMWAREA,STEMFORM,
     &  RMFR,RMCR,Q10R,RTEMPR,EFFYRF,   
     &  RMB,Q10B,RTEMPB,  
     &  GSREF,GSMIN,PAR0,D0,VK1,VK2,VPD1,VPD2,VMFD0,
     &  GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,
     &  G0,D0L,GAMMA,G1,WLEAF,NSIDES)

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

C Get input from soil/understorey file ! for water balance - unfinished
C	IF (IOSOIL.EQ.1) CALL INPUTSOIL(
C	&  ISTART,
C	&  PF,PS,CANCAP,STEMCAP,STEMCANERAT,
C	&  CANWATER,STEMWATER,
C	&  NSOILLAY,SOILWATER,SOILCAP,SOILDEPTH
C	&  )

C Open output files
      CALL OPENOUTPUTF(IODAILY,IOHRLY,IOHIST,IORESP,IOSOIL,
     &  CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE)

C Do calculations which are not day-dependent
      CALL EXDIFF(NALPHA,ALPHA,FALPHA,NZEN,DIFZEN,RANDOM,DEXT)

C**********************************************************************
C Loop over all chosen trees within subplot
      DO 50 ITAR = 1,NOTARGETS

      ITREE = ITARGETS(ITAR)
      WRITE (*,106) ITREE
 106  FORMAT ('  TREE:',I5)

      CALL SORTTREES(
     &  NOALLTREES,NOTREES,ITREE,
     &  DXT1,DYT1,DZT1,RXTABLE1,RYTABLE1,RZTABLE1,
     &  ZBCTABLE1,FOLTABLE1,DIAMTABLE1,
     &  DXT,DYT,DZT,RXTABLE,RYTABLE,RZTABLE,
     &  FOLTABLE,ZBCTABLE,DIAMTABLE
     &  )

      CALL RESTARTMETF(ISTART,MSTART,MFLAG)

      CALL ZEROTREE(HISTO,CANWATER,STEMWATER)

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

C Interpolate to get daily values of parameters
      CALL INTERPOLATEP(IDAY,ISTART,NOJDATES,DATESJ,JMAXTABLE,
     &  NOVDATES,DATESV,VCMAXTABLE,NORDATES,DATESRD,RDTABLE,
     &  NOSLADATES,DATESSLA,SLATABLE,NOADATES,DATESA,AJQTABLE,
     &  NOLAY,NOAGEP,JMAX25,VCMAX25,RD0,SLA,AJQ)
      CALL INTERPOLATET(IDAY,ISTART,
     &  NOXDATES,DATESX,RXTABLE,NOYDATES,DATESY,RYTABLE,
     &  NOZDATES,DATESZ,RZTABLE,NOTDATES,DATEST,ZBCTABLE,
     &  NODDATES,DATESD,DIAMTABLE,
     &  NOLADATES,DATESLA,FOLTABLE,TOTLAITABLE,NOTREES,
     &  RX,RY,RZ,ZBC,FOLT,TOTLAI,DIAM,STOCKING,
     &  IFLUSH,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN,
     &  NEWCANOPY)

      IF (NEWCANOPY.EQ.1) THEN
C If the canopy has changed (or on the first day) the grid points must be set up
        CALL Points(
     &    NUMPNT,JLEAF,JSHAPE,SHAPE,RX(1),RY(1),RZ(1),
     &    ZBC(1),DXT(1),DYT(1),DZT(1),FOLT(1),
     &    PROP,BPT,NOAGEC,NOAGEP,NOLAY,
     &    XL,YL,ZL,VL,DLT,DLI,LGP,FOLLAY
     &  )
      END IF

C Following functions need FOLLAY. 
      CALL GETWIND(FOLLAY,FOLT(1),TOTLAI,EXTWIND,WINDLAY)
C Calculate woody biomass and woody biomass increment
      CALL CALCWBIOM(IDAY,RZ(1)+ZBC(1),DIAM(1),COEFFT,EXPONT,WINTERC,
     &  WBIOM,WBINC)
      CALL CALCWBIOM(IDAY,RZ(1)+ZBC(1),DIAM(1),BCOEFFT,BEXPONT,BINTERC,
     &  BBIOM,BBINC)
      CALL CALCWBIOM(IDAY,RZ(1)+ZBC(1),DIAM(1),RCOEFFT,REXPONT,RINTERC,
     &  RBIOM,RBINC)
C Calculate foliar biomass and increment
      CALL CALCFBIOM(IDAY,NOSLADATES,FOLLAY,SLA,PROP,NOLAY,NOAGEP,
     &  FBIOM,FBINC)
C Calculate stem respiration rate per unit biomass
	RMW = CALCRMW(MODELRW,COLLA,COLLK,STEMSDW,
     &  DIAM(1),RZ(1)+ZBC(1),STEMFORM,RMWAREA,WBIOM,RMW)

C Output information to layer flux file if required
      IF (IOHRLY.GT.1)
     &  CALL OUTPUTLAY(ULAY,FOLLAY,JMAX25,VCMAX25,NOLAY)

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

C Zero daily fluxes
      CALL ZEROD(TDYAB,TOTCO2,TOTRESPF,TOTRESPWM,
     &  TOTRESPB,TOTRESPCR,TOTRESPFR,TOTH2O,TOTHFX,TOTDRAIN)
C Zero hourly fluxes
      CALL ZEROHR(THRAB,FCO2,FRESPF,FRESPW,FRESPB,FRESPFR,FRESPCR,
     &  FH2OT,FH2OE,GSCAN,FHEAT,PPAR,PPS,PTRANSP,TCAN)

C**********************************************************************
C Begin hourly loop
      DO 100 IHOUR = 1,KHRS

C Partition rainfall ! for water balance - unfinished
C      CALL PartitionPPT(PPT(IHOUR),PFLA,CANCAPLA,FOLT(1),CANWATER)

C Test to see if daylight hours or if any foliage
        IF ((ABS(ZEN(IHOUR)).LE.1.57).AND.
     &    (RADABV(IHOUR,1).GT.0.1).AND.
     &    (FOLT(1).GT.0.0)) THEN

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

C Loop over the 3 wavelengths
            DO 90 IWAVE = 1,3

C Calculate the scattered radiation
              CALL SCATTER(IPT,IWAVE,
     &          MLAYER(IPT),LAYER(IPT),DLAI,EXPDIF,
     &          ZEN(IHOUR),BEXT,
     &          DMULT2,SOMULT,BMULT,
     &          RADABV(IHOUR,IWAVE),FBEAM(IHOUR,IWAVE),
     &          TAIR(IHOUR),TSOIL(IHOUR),
     &          ARHO(LGP(IPT),IWAVE),ATAU(LGP(IPT),IWAVE),RHOSOL(IWAVE),
     &          DIFUP,DIFDN,SCLOST)

C Calculate absorbed radiation
              CALL ABSRAD(IPT,IWAVE,
     &          NZEN,DEXT,BEXT,BMULT,RELDF(IPT),
     &          RADABV(IHOUR,IWAVE),FBEAM(IHOUR,IWAVE),ZEN(IHOUR),
     &          ABSRP(LGP(IPT),IWAVE),DIFDN(IPT,IWAVE),DIFUP(IPT,IWAVE),
     &          DFLUX,BFLUX,SCATFX)

90          CONTINUE ! End loop over wavelengths

C Calculation of photosynthesis may be done for sunlit & shaded leaves
C separately, or by averaging PAR over the total leaf area

            IF ((MODELSS.EQ.0).AND.(FBEAM(IHOUR,1).NE.0.0)) THEN 
C Do calculations separately for sunlit & shaded leaves

            DO 95 ISUNLIT = 1,2 ! Loop over sunlit & shaded leaves

              IF (ISUNLIT.EQ.1) THEN
                APAR = (BFLUX(IPT,1)*BEXT + DFLUX(IPT,1))*UMOLPERJ
                ANIR = BFLUX(IPT,2)*BEXT + DFLUX(IPT,2)
                FAREA = SUNLA
              ELSE
                APAR = DFLUX(IPT,1)*UMOLPERJ
                ANIR = DFLUX(IPT,2)
                FAREA = 1.0 - SUNLA
              END IF

              ATHR = DFLUX(IPT,3)
              RNET = APAR/UMOLPERJ + ANIR + ATHR

              DO 95 IAGE = 1,NOAGEP ! Loop over age classes
                AREA = FAREA * DLI(IAGE,IPT) * VL(IPT) ! m2
                IF (IOHIST.EQ.1) CALL CATEGO(AREA,APAR,HISTO,BINSIZE)

C Call physiology routine
                CALL PSTRANSP(RELDF(IPT),TU(IPT),TD(IPT),
     +  RNET,
     +  WINDAH(IHOUR)*WINDLAY(LGP(IPT)),APAR,TAIR(IHOUR),CA(IHOUR),
     +  RH(IHOUR),VPD(IHOUR),VMFD(IHOUR),PRESS(IHOUR),SOILMD(IHOUR),
     +  JMAX25(LGP(IPT),IAGE),IECO,EAVJ,EDVJ,DELSJ,
     +  VCMAX25(LGP(IPT),IAGE),EAVC,EDVC,DELSC,TVJUP,TVJDN,
     +  THETA,AJQ(LGP(IPT),IAGE),RD0(LGP(IPT),IAGE),
     &  Q10F,RTEMP,DAYRESP,TBELOW,
     +  MODELGS,GSREF,GSMIN,PAR0,D0,VK1,VK2,VPD1,VPD2,VMFD0,
     &  GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,
     +  G0,D0L,GAMMA,G1,WLEAF,NSIDES,CANWATER,ITERMAX,
     +  GSC,ALEAF,RD,ET,EI,HFX,TLEAF,GBH)

C Sum outputs for the hour
              CALL SUMHR(APAR/UMOLPERJ,ANIR,ATHR,
     +  ALEAF,RD,GSC,ET,EI,HFX,TLEAF,
     +  AREA,IHOUR,LGP(IPT),
     +  PPAR,PPS,PTRANSP,
     +  THRAB,FCO2,FRESPF,GSCAN,FH2OT,FH2OE,FHEAT,TCAN)

95          CONTINUE ! End loop over sunlit / shaded leaves

            ELSE IF ((MODELSS.EQ.1).OR.(FBEAM(IHOUR,1).EQ.0.0)) THEN 
C Do calculations for PAR averaged over all leaf area

              APAR = (BFLUX(IPT,1)*BEXT*SUNLA + DFLUX(IPT,1))*UMOLPERJ
              ANIR = BFLUX(IPT,2)*BEXT*SUNLA + DFLUX(IPT,2)
              ATHR = DFLUX(IPT,3)
              RNET = APAR/UMOLPERJ + ANIR + ATHR

              DO 125 IAGE = 1,NOAGEP ! Loop over age classes
                AREA = DLI(IAGE,IPT) * VL(IPT)
                IF (IOHIST.EQ.1) CALL CATEGO(AREA,APAR,HISTO,BINSIZE)

C Call physiology routine
                CALL PSTRANSP(RELDF(IPT),TU(IPT),TD(IPT),RNET,
     +  WINDAH(IHOUR)*WINDLAY(LGP(IPT)),APAR,TAIR(IHOUR),CA(IHOUR),
     +  RH(IHOUR),VPD(IHOUR),VMFD(IHOUR),PRESS(IHOUR),SOILMD(IHOUR),
     +  JMAX25(LGP(IPT),IAGE),IECO,EAVJ,EDVJ,DELSJ,
     +  VCMAX25(LGP(IPT),IAGE),EAVC,EDVC,DELSC,TVJUP,TVJDN,
     +  THETA,AJQ(LGP(IPT),IAGE),RD0(LGP(IPT),IAGE),
     &  Q10F,RTEMP,DAYRESP,TBELOW,
     +  MODELGS,GSREF,GSMIN,PAR0,D0,VK1,VK2,VPD1,VPD2,VMFD0,
     &  GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,
     +  G0,D0L,GAMMA,G1,WLEAF,NSIDES,CANWATER,ITERMAX,
     +  GSC,ALEAF,RD,ET,EI,HFX,TLEAF,GBH)

C Sum outputs for the hour
              CALL SUMHR(APAR/UMOLPERJ,ANIR,ATHR,
     +  ALEAF,RD,GSC,ET,EI,HFX,TLEAF,
     +  AREA,IHOUR,LGP(IPT),
     +  PPAR,PPS,PTRANSP,
     +  THRAB,FCO2,FRESPF,GSCAN,FH2OT,FH2OE,FHEAT,TCAN)

125          CONTINUE ! End loop over age classes

            ELSE IF ((MODELSS.EQ.2).AND.(FBEAM(IHOUR,1).NE.0.0)) THEN 
! Do calculations separately for sunlit & shaded. Further separate sunlit
! into leaf angle classes.

            DO 135 ISUNLIT = 1,NALPHA+1

              IF (ISUNLIT.GT.NALPHA) THEN
                APAR = DFLUX(IPT,1)*UMOLPERJ
                ANIR = DFLUX(IPT,2)
                FAREA = 1.0 - SUNLA
              ELSE
                APAR = (BFLUX(IPT,1)*BEXTANG(ISUNLIT)
     &               + DFLUX(IPT,1))*UMOLPERJ
                ANIR = BFLUX(IPT,2)*BEXTANG(ISUNLIT) + DFLUX(IPT,2)
                FAREA = SUNLA*FALPHA(ISUNLIT)
              END IF

              ATHR = DFLUX(IPT,3)
              RNET = APAR/UMOLPERJ + ANIR + ATHR

              DO 135 IAGE = 1,NOAGEP ! Loop over age classes
                AREA = FAREA * DLI(IAGE,IPT) * VL(IPT) ! m2
                IF (IOHIST.EQ.1) CALL CATEGO(AREA,APAR,HISTO,BINSIZE)

C Call physiology routine
                CALL PSTRANSP(RELDF(IPT),TU(IPT),TD(IPT),
     +  RNET,
     +  WINDAH(IHOUR)*WINDLAY(LGP(IPT)),APAR,TAIR(IHOUR),CA(IHOUR),
     +  RH(IHOUR),VPD(IHOUR),VMFD(IHOUR),PRESS(IHOUR),SOILMD(IHOUR),
     +  JMAX25(LGP(IPT),IAGE),IECO,EAVJ,EDVJ,DELSJ,
     +  VCMAX25(LGP(IPT),IAGE),EAVC,EDVC,DELSC,TVJUP,TVJDN,
     +  THETA,AJQ(LGP(IPT),IAGE),RD0(LGP(IPT),IAGE),
     &  Q10F,RTEMP,DAYRESP,TBELOW,
     +  MODELGS,GSREF,GSMIN,PAR0,D0,VK1,VK2,VPD1,VPD2,VMFD0,
     &  GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,
     +  G0,D0L,GAMMA,G1,WLEAF,NSIDES,CANWATER,ITERMAX,
     +  GSC,ALEAF,RD,ET,EI,HFX,TLEAF,GBH)

C Sum outputs for the hour
              CALL SUMHR(APAR/UMOLPERJ,ANIR,ATHR,
     +  ALEAF,RD,GSC,ET,EI,HFX,TLEAF,
     +  AREA,IHOUR,LGP(IPT),
     +  PPAR,PPS,PTRANSP,
     +  THRAB,FCO2,FRESPF,GSCAN,FH2OT,FH2OE,FHEAT,TCAN)

135          CONTINUE ! End loop over sunlit / shaded leaves

          END IF ! Separating sunlit & shaded foliage, or not

80        CONTINUE ! End loop over grid points

C Calculate transpiration by applying Penman-Monteith to canopy
          FH2OCAN(IHOUR) = ETCAN(WINDAH(IHOUR),ZHT,Z0HT,ZPD,
     +      PRESS(IHOUR),TAIR(IHOUR),
     +      THRAB(IHOUR,1)+THRAB(IHOUR,2)+THRAB(IHOUR,3),
     +      VPD(IHOUR),GSCAN(IHOUR),STOCKING)
C Normalise to get average foliage temperature
		TCAN(IHOUR) = TCAN(IHOUR)/FOLT(1)

        ELSE ! Night-time

C Calculate night-time foliage respiration
          DO 120 IPT = 1,NUMPNT
            DO 120 IAGE = 1,NOAGEP
              AREA = DLI(IAGE,IPT) * VL(IPT)
              RESPF =   AREA *
     +  RESP(RD0(LGP(IPT),IAGE),TAIR(IHOUR),Q10F,RTEMP,1.0,TBELOW)
              FCO2(IHOUR) = FCO2(IHOUR) - RESPF
              FRESPF(IHOUR) = FRESPF(IHOUR) + RESPF
120       CONTINUE

C Calculate night-time evaporation of rainfall and dewfall. (To be added, BM 04/00)

C Canopy T at nighttime is same as air temperature. 
		TCAN(IHOUR) = TAIR(IHOUR)

        END IF ! If day or night

C Calculate non-foliage maintenance respiration (in umol tree-1 s-1)
        FRESPW(IHOUR) = RESP(RMW,TAIR(IHOUR),Q10W,RTEMPW,1.0,TBELOW) 
     &    * WBIOM
        FRESPB(IHOUR) = RESP(RMB,TAIR(IHOUR),Q10B,RTEMPB,1.0,TBELOW) 
     &    * BBIOM
        FRESPFR(IHOUR) = RESP(RMFR,TSOIL(IHOUR),Q10R,RTEMPR,1.0,TBELOW) 
     &    * RBIOM * FRFRAC
        FRESPCR(IHOUR) = RESP(RMCR,TSOIL(IHOUR),Q10R,RTEMPR,1.0,TBELOW) 
     &    * RBIOM * (1. - FRFRAC)

C Calculate evaporation and actual transpiration: ignore stem component 
C for water balance - unfinished
C	    EPOT = FH2OE(IHOUR)
C		TPOT = FH2OT(IHOUR)
C		CALL EvTran(Tpot,Epot,CANCAPLA,FOLT(1),CanWater,
C     &	  FH2OT(IHOUR),FH2OE(IHOUR))
c      WRITE (100,140) IDAY,IHOUR,PPT(IHOUR),EPOT,FH2OE(IHOUR),
c     &	TPOT,FH2OT(IHOUR),CANWATER
c140   FORMAT (I2,1X,I2,6(1X,F12.3))

C Output hourly totals
        CALL OUTPUTHR(IOHRLY,IDAY+1,IHOUR,ITREE,TCAN,
     +    NOLAY,PPAR,PPS,PTRANSP,FOLLAY,
     +    THRAB,FCO2,FRESPF,FRESPW,FRESPB,
     +    FH2OT,FH2OE,GSCAN,FH2OCAN,FHEAT)

C**********************************************************************
100   CONTINUE
C End hourly loop

C Calculate daily growth respiration
      TOTRESPWG = GRESP(WBINC,EFFYRW)
      TOTRESPBG = GRESP(BBINC,EFFYRW)
      TOTRESPFRG = GRESP(RBINC*FRFRAC,EFFYRF)
      TOTRESPCRG = GRESP(RBINC*(1.-FRFRAC),EFFYRW)
      TOTRESPFG = GRESP(FBINC,EFFYRF)

C Output daily totals
      CALL SUMDAILY(THRAB,FCO2,FRESPF,FRESPW,FRESPB,FRESPCR,FRESPFR,
     &  FH2OT,FH2OE,FH2OCAN,FHEAT,
     &  TDYAB,TOTCO2,TOTRESPF,TOTRESPWM,TOTRESPB,TOTRESPCR,TOTRESPFR,
     &  TOTH2O,TOTH2OCAN,TOTHFX)
      CALL OUTPUTDY(IDAY+1,ITREE,IODAILY,TDYAB,TOTCO2,
     &  TOTRESPF,TOTRESPWM,TOTRESPWG,
     &  TOTH2O,TOTH2OCAN,TOTHFX,
     &  IORESP,TOTRESPCR,TOTRESPFR,TOTRESPFRG,TOTRESPCRG,TOTRESPFG,
     &  TOTRESPB,TOTRESPBG)

C Go to next day of simulation
      IDAY = IDAY + NSTEP
      IF ((ISTART+IDAY).LE.IEND) THEN
        CALL SKIPMET(MFLAG,NSTEP)
        GOTO 10
      END IF

C**********************************************************************
110   CONTINUE
C End daily loop

      IF (IOHIST.EQ.1) CALL OUTPUTHIST(UHIST,ITREE,HISTO,BINSIZE)

C**********************************************************************
50    CONTINUE 
C End loop over trees
C**********************************************************************

C Write diffuse transmittances to file
      REWIND (UTUTD)
      WRITE (UTUTD,1313) (IPT,TU(IPT),TD(IPT),RELDF(IPT),IPT=1,NUMPNT)
1313  FORMAT (1X,I3,5X,F8.4,1X,F8.4,1X,F8.4)

      CALL Closef(IODAILY,IOHRLY,IOHIST,IORESP,IOSOIL)

      STOP
      END !MAIN
C**********************************************************************


C**********************************************************************
      SUBROUTINE ZEROTREE(HISTO,CANWATER,STEMWATER)
C Set values particular to each tree to zero. 
C**********************************************************************

      INCLUDE 'maestcom'
      REAL HISTO(MAXHISTO)

      DO 10 I = 1,MAXHISTO
        HISTO(I) = 0.0
10    CONTINUE
      CANWATER = 0.0
	STEMWATER = 0.0

      RETURN
      END !ZeroTree

C**********************************************************************
      SUBROUTINE ZEROD(
     &  TDYAB,TOTCO2,TOTRESPF,TOTRESPWM,TOTRESPB,
     &  TOTRESPCR,TOTRESPFR,TOTH2O,TOTHFX,TOTDRAIN
     &  )
C This is subroutine to set the initial values of daily total variables
C to zero.
C**********************************************************************

      INCLUDE 'maestcom'
      REAL TDYAB(3)

      TDYAB(1) = 0.0
      TDYAB(2) = 0.0
      TDYAB(3) = 0.0
      TOTCO2 = 0.0
      TOTRESPF = 0.0
      TOTRESPWM = 0.0
      TOTRESPB = 0.0
      TOTRESPFR = 0.0
      TOTRESPCR = 0.0
      TOTH2O = 0.0
      TOTHFX = 0.0
	TOTDRAIN = 0.0

      RETURN
      END !ZeroD


C**********************************************************************
      SUBROUTINE ZEROHR(
     &  THRAB, FCO2, FRESPF, FRESPW, FRESPB, FRESPFR, FRESPCR, 
     &  FH2OT, FH2OE, GSCAN, FHEAT, PPAR, PPS, PTRANSP, TCAN
     &  )
C This is subroutine to set the initial values of hourly total variables
C to zero.
C**********************************************************************

      INCLUDE 'maestcom'
      REAL THRAB(KHRS,3),TCAN(KHRS)
      REAL FCO2(KHRS),FRESPF(KHRS),FRESPW(KHRS),FRESPB(KHRS)
      REAL FRESPCR(KHRS),FRESPFR(KHRS)
      REAL GSCAN(KHRS),FH2OT(KHRS),FH2OE(KHRS),FHEAT(KHRS)
      REAL PPAR(MAXLAY,KHRS),PPS(MAXLAY,KHRS),PTRANSP(MAXLAY,KHRS)

      DO 10 I = 1,KHRS
        FCO2(I) = 0.0
        FRESPF(I) = 0.0
        FRESPW(I) = 0.0
        FRESPB(I) = 0.0
        FRESPFR(I) = 0.0
        FRESPCR(I) = 0.0
        FH2OT(I) = 0.0
        FH2OE(I) = 0.0
        GSCAN(I) = 0.0
        FHEAT(I) = 0.0
	  TCAN(I) = 0.0
        DO 20 J = 1,3
          THRAB(I,J) = 0.0
20      CONTINUE
        DO 30 J = 1,MAXLAY
          PPAR(J,I) = 0.0
          PPS(J,I) = 0.0
          PTRANSP(J,I) = 0.0
30      CONTINUE

10    CONTINUE

      RETURN
      END !ZeroHr

C**********************************************************************
      SUBROUTINE SUMHR(APAR,ANIR,ATHR,
     +  ALEAF,RD,GSC,ET,EI,HFX,TLEAF,AREA,IHOUR,ILAY,
     +  PPAR,PPS,PTRANSP,
     +  THRAB,FCO2,FRESPF,GSCAN,FH2OT,FH2OE,FHEAT,TCAN)
C Sum fluxes from each point to give hourly fluxes.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      REAL THRAB(KHRS,3),FCO2(KHRS),FRESPF(KHRS),TCAN(KHRS)
      REAL GSCAN(KHRS),FH2OT(KHRS),FH2OE(KHRS),FHEAT(KHRS)
      REAL PPAR(MAXLAY,KHRS),PPS(MAXLAY,KHRS),PTRANSP(MAXLAY,KHRS)

C Sum PAR, photosynthesis, & transpiration by layer
      PPAR(ILAY,IHOUR) = PPAR(ILAY,IHOUR) +  APAR*UMOLPERJ*AREA
      PPS(ILAY,IHOUR) = PPS(ILAY,IHOUR) + ALEAF*AREA
      PTRANSP(ILAY,IHOUR) = PTRANSP(ILAY,IHOUR) + ET*AREA

C Sum all fluxes for the hour
      THRAB(IHOUR,1) = THRAB(IHOUR,1) + APAR*AREA
      THRAB(IHOUR,2) = THRAB(IHOUR,2) + ANIR*AREA
      THRAB(IHOUR,3) = THRAB(IHOUR,3) + ATHR*AREA
      FCO2(IHOUR) = FCO2(IHOUR) + ALEAF*AREA
C Foliage respiration in umol tree-1 s-1
      FRESPF(IHOUR) = FRESPF(IHOUR) + RD*AREA
C Transpiration in umol tree-1 s-1
      FH2OT(IHOUR) = FH2OT(IHOUR) + ET*AREA
C Evaporation in umol tree-1 s-1
      FH2OE(IHOUR) = FH2OE(IHOUR) + EI*AREA
C Stom cond in mol tree-1 s-1
      GSCAN(IHOUR) = GSCAN(IHOUR) + GSC*AREA
C Heat flux in mol tree-1 s-1
      FHEAT(IHOUR) = FHEAT(IHOUR) + HFX*AREA
C Average leaf temperature - will be divided by total leaf area later. 
	TCAN(IHOUR) = TCAN(IHOUR) + TLEAF*AREA

      RETURN
      END !SumHR


C**********************************************************************
      SUBROUTINE SUMDAILY(
     +  THRAB,FCO2,FRESPF,FRESPW,FRESPB,FRESPCR,FRESPFR,
     +  FH2OT,FH2OE,FH2OCAN,FHEAT,
     +  TDYAB,TOTCO2,TOTRESPF,TOTRESPWM,TOTRESPB,TOTRESPCR,TOTRESPFR,
     +  TOTH2O,TOTH2OCAN,TOTHFX)
C Sum hourly fluxes to give daily ones
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL THRAB(KHRS,3),FCO2(KHRS),FRESPF(KHRS),FRESPW(KHRS)
      REAL FRESPB(KHRS),FRESPCR(KHRS),FRESPFR(KHRS)
      REAL FH2OT(KHRS),FH2OE(KHRS),FH2OCAN(KHRS),FHEAT(KHRS)
      REAL TDYAB(3)

      DO 10 IHOUR = 1,KHRS
        TOTCO2 = TOTCO2 + FCO2(IHOUR)
        TOTRESPF = TOTRESPF + FRESPF(IHOUR)
        TOTRESPWM = TOTRESPWM + FRESPW(IHOUR)
        TOTRESPB  = TOTRESPB  + FRESPB(IHOUR)
        TOTRESPCR = TOTRESPCR + FRESPCR(IHOUR)
        TOTRESPFR = TOTRESPFR + FRESPFR(IHOUR)
        TOTH2O = TOTH2O + FH2OT(IHOUR) + FH2OE(IHOUR)
        TOTH2OCAN = TOTH2OCAN + FH2OCAN(IHOUR)
        TOTHFX = TOTHFX + FHEAT(IHOUR)
        DO 10 J = 1,3
          TDYAB(J) = TDYAB(J) + THRAB(IHOUR,J)
10    CONTINUE

      CONVERT = SPERHR*1E-6

      TOTCO2 = TOTCO2*CONVERT
      TOTRESPF = TOTRESPF*CONVERT
      TOTRESPB = TOTRESPB*CONVERT
      TOTRESPWM = TOTRESPWM*CONVERT
      TOTRESPCR = TOTRESPCR*CONVERT
      TOTRESPFR = TOTRESPFR*CONVERT
      TOTH2O = TOTH2O*CONVERT
      TOTH2OCAN = TOTH2OCAN*CONVERT
      TOTHFX = TOTHFX*CONVERT
      DO 20 J = 1,3
        TDYAB(J) = TDYAB(J)*CONVERT
20    CONTINUE

      RETURN
      END !SumDaily


