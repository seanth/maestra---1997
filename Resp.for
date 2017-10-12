
C**********************************************************************
      PROGRAM RESPIRN
C This program is designed to be run in conjunction with MAESTRA.
C It uses the same input files but it only calculates respiration rates.
C It writes to the file 'resp.dat'.
C**********************************************************************

      INCLUDE 'maestcom'

C Array declarations.
C List of trees for which to do calculations
      INTEGER ITARGETS(MAXT)
C Tree positions and dimensions - all trees, all dates
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
	REAL DELTAT(12)
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
      REAL DIFZEN(MAXANG),ZEN(KHRS),AZ(KHRS)
      REAL TU(MAXP),TD(MAXP),RELDF(MAXP)
      REAL DIFDN(MAXP,3),DIFUP(MAXP,3),SCLOST(MAXP,3)
      REAL BFLUX(MAXP,3),DFLUX(MAXP,3),SCATFX(MAXP,3)
      REAL TDYAB(3)
C Titles of input files
      CHARACTER*80 CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE
C Biomass for each tree
      REAL FBIOM(MAXT),WBIOM(MAXT),RBIOM(MAXT),BBIOM(MAXT)
	REAL RMWVOL(MAXT)

C**********************************************************************
C Program starts.
      VTITLE = 'Program: Maestra Respiration Version: Sep 2001'

C For this program MAXTT must = MAXT: because I'm not sorting the trees.
C Make sure this is OK when compiling!
      IF (MAXTT.NE.MAXT)
     &  CALL SUBERROR('ERROR: MAXTT MUST = MAXT: RECOMPILE',
     &    IFATAL,0)

C Open input files
      CALL OpenInputf(IODAILY,IOHRLY,IOTUTD,IOHIST,IORESP,IOSOIL,
     &  CTITLE,TTITLE,PTITLE,STITLE,UTITLE)
C Set flags determining which input files are needed
      IORESP = 1
      IODAILY = 0
      IOHRLY = 0
      IOHIST = 0

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
     &  DXT,DYT,DZT,RXTABLE,RYTABLE,RZTABLE,ZBCTABLE,
     &  FOLTABLE,TOTLAITABLE,DIAMTABLE,
     &  IFLUSH,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN
     &  )

C Open met data file (must be done after ISTART & IEND read)
      CALL OPENMETF(ISTART,IEND,CAK,PRESSK,SWMIN,SWMAX,
     &  DIFSKY,ALAT,TTIMD,DELTAT,
     &  MFLAG,METCOLS,NOMETCOLS,MTITLE,MSTART)
      CALL RESTARTMETF(ISTART,MSTART,MFLAG)

C Open output files
      CALL OPENOUTPUTF(IODAILY,IOHRLY,IOHIST,IORESP,IOSOIL,
     &  CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE)

C**********************************************************************
C Begin daily loop
      IDAY = 0
10    CONTINUE
      WRITE (*,105) IDAY
105   FORMAT ('  DAY:',I5)

C Zero daily totals
      CALL ZEROD(DRESPFM,DRESPFG,DRESPWM,DRESPWG,
     &  DRESPBM,DRESPBG,
     &  DRESPCRM,DRESPCRG,DRESPFRM,DRESPFRG,
     &  TFBIOM,TWBIOM,TBBIOM,TRBIOM)

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
     &  NOLADATES,DATESLA,FOLTABLE,TOTLAITABLE,NOALLTREES,
     &  RX,RY,RZ,ZBC,FOLT,TOTLAI,DIAM,STOCKING,
     &  IFLUSH,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN,
     &  NEWCANOPY)

C**********************************************************************
C Loop over all chosen trees within subplot
      DO 50 ITAR = 1,NOTARGETS

      ITREE = ITARGETS(ITAR)
      WRITE (*,106) ITREE
 106  FORMAT ('  TREE:',I5)

C Calculate foliage area in each layer (will be converted to biomass later)
      IF (FOLT(ITREE).GT.0.0) THEN
        CALL Points(
     &    NUMPNT,JLEAF,JSHAPE,SHAPE,RX(ITREE),RY(ITREE),RZ(ITREE),
     &    ZBC(ITREE),DXT(ITREE),DYT(ITREE),DZT(ITREE),FOLT(ITREE),
     &    PROP,BPT,NOAGEC,NOAGEP,NOLAY,
     &    XL,YL,ZL,VL,DLT,DLI,LGP,FOLLAY
     &  )
      ELSE
        DO 60 ILAY = 1,NOLAY
          FOLLAY(ILAY) = 0.0
60      CONTINUE
      END IF

C Calculate foliar biomass and increment
      CALL CALCFBIOM(IDAY,NOSLADATES,FOLLAY,SLA,PROP,NOLAY,NOAGEP,
     &  FBIOM(ITREE),FBINC)
      TFBIOM = TFBIOM + FBIOM(ITREE)

C Calculate woody biomass and woody biomass increment
      CALL CALCWBIOM(IDAY,RZ(ITREE)+ZBC(ITREE),DIAM(ITREE),
     &  COEFFT,EXPONT,WINTERC,
     &  WBIOM(ITREE),WBINC)
      CALL CALCWBIOM(IDAY,RZ(ITREE)+ZBC(ITREE),DIAM(ITREE),
     &  BCOEFFT,BEXPONT,BINTERC,
     &  BBIOM(ITREE),BBINC)
      CALL CALCWBIOM(IDAY,RZ(ITREE)+ZBC(ITREE),DIAM(ITREE),
     &  RCOEFFT,REXPONT,RINTERC,
     &  RBIOM(ITREE),RBINC)
      TWBIOM = TWBIOM + WBIOM(ITREE)
      TBBIOM = TBBIOM + BBIOM(ITREE)
      TRBIOM = TRBIOM + RBIOM(ITREE)

C Calculate stem respiration rate per unit biomass
	RMWVOL(ITREE) = CALCRMW(MODELRW,COLLA,COLLK,STEMSDW,
     &  DIAM(ITREE),RZ(ITREE)+ZBC(ITREE),
     &  STEMFORM,RMWAREA,WBIOM(ITREE),RMW)

C**********************************************************************
C Begin hourly loop
      DO 100 IHOUR = 1,KHRS

        IF (NIGHT(ZEN(IHOUR),RADABV(IHOUR,1)).EQ.0) THEN
          RESPL = DAYRESP
        ELSE
          RESPL = 1.0
        END IF 

C Calculate foliage respiration
C        DO 120 IPT = 1,NUMPNT
C          DO 120 IAGE = 1,NOAGEP
C            AREA = DLI(IAGE,IPT) * VL(IPT)
C            RESPF =   AREA *
C     +  RESP(RD0(LGP(IPT),IAGE),TAIR(IHOUR),Q10F,RTEMP,RESPL,TBELOW)
C            DRESPFM = DRESPFM + RESPF * SPERHR * 1E-6
C120     CONTINUE
C Short-cut since all foliage has same resp. rate
        DRESPFM = DRESPFM 
     &    + RESP(RD0(1,1),TAIR(IHOUR),Q10F,RTEMP,RESPL,TBELOW) 
     &    * FOLT(ITREE) * SPERHR * 1E-6

C Calculate non-foliage maintenance respiration (in umol tree-1 s-1)
        DRESPWM = DRESPWM + RESP(RMWVOL(ITREE),TAIR(IHOUR),Q10W,
     &    RTEMPW,1.0,TBELOW) * WBIOM(ITREE) * SPERHR * 1E-6
        DRESPBM = DRESPBM + RESP(RMB,TAIR(IHOUR),Q10B,RTEMPB,1.0,TBELOW) 
     &    * BBIOM(ITREE) * SPERHR * 1E-6
        DRESPFRM = DRESPFRM + 
     &    RESP(RMFR,TSOIL(IHOUR),Q10R,RTEMPR,1.0,TBELOW) 
     &    * RBIOM(ITREE) * FRFRAC * SPERHR * 1E-6
        DRESPCRM = DRESPCRM 
     &    + RESP(RMCR,TSOIL(IHOUR),Q10R,RTEMPR,1.0,TBELOW) 
     &    * RBIOM(ITREE) * (1. - FRFRAC) * SPERHR * 1E-6
C**********************************************************************
100   CONTINUE
C End hourly loop

C Calculate daily growth respiration
      DRESPWG = DRESPWG + GRESP(WBINC,EFFYRW)
      DRESPBG = DRESPBG + GRESP(BBINC,EFFYRW)
      DRESPFRG = DRESPFRG + GRESP(RBINC*FRFRAC,EFFYRF)
      DRESPCRG = DRESPCRG + GRESP(RBINC*(1.-FRFRAC),EFFYRW)
      DRESPFG = DRESPFG + GRESP(FBINC,EFFYRF)

C**********************************************************************
50    CONTINUE 
C End loop over trees

C Convert totals to mol tree-1 s-1 
      CALL PERTREE(DRESPFM,DRESPWM,DRESPWG,DRESPCRM,DRESPFRM,
     &  DRESPFRG,DRESPCRG,DRESPFG,DRESPBM,DRESPBG,NOTARGETS)

C Output daily totals
      CALL OUTPUTDY(IDAY+1,0,IODAILY,TDYAB,TOTCO2,
     &  DRESPFM,DRESPWM,DRESPWG,
     &  TOTH2O,TOTH2OCAN,TOTHFX,
     &  IORESP,DRESPCRM,DRESPFRM,DRESPFRG,DRESPCRG,DRESPFG,
     &  DRESPBM,DRESPBG)

C Go to next day of simulation
      IDAY = IDAY + NSTEP
      IF ((ISTART+IDAY).LE.IEND) THEN
        CALL SKIPMET(MFLAG,NSTEP)
        GOTO 10
      END IF

C**********************************************************************
110   CONTINUE
C End daily loop

C**********************************************************************
      CALL Closef(IODAILY,IOHRLY,IOHIST,IORESP,IOSOIL)

      STOP
      END !MAIN
C**********************************************************************


C**********************************************************************
      SUBROUTINE ZEROD(
     &  DRESPFM,DRESPFG,DRESPWM,DRESPWG,DRESPBM,DRESPBG,
     &  DRESPCRM,DRESPCRG,DRESPFRM,DRESPFRG,
     &  TFBIOM, TWBIOM, TBBIOM, TRBIOM)
C This is subroutine to set the initial values of daily total variables
C to zero.
C**********************************************************************

      DRESPFM = 0.0
      DRESPFG = 0.0
      DRESPWM = 0.0
      DRESPWG = 0.0
      DRESPBM = 0.0
      DRESPBG = 0.0
      DRESPCRM = 0.0
      DRESPCRG = 0.0
      DRESPFRM = 0.0
      DRESPFRG = 0.0
      TFBIOM = 0.0
      TWBIOM = 0.0
      TBBIOM = 0.0
      TRBIOM = 0.0

      RETURN
      END !ZeroD


C**********************************************************************
      SUBROUTINE PERTREE(DRESPFM,DRESPWM,DRESPWG,DRESPCRM,DRESPFRM,
     &  DRESPFRG,DRESPCRG,DRESPFG,DRESPBM,DRESPBG,NOTARGETS)
C Convert totals (for all trees) to per tree basis - divide by no of trees
C**********************************************************************

      DRESPFM  =   DRESPFM / REAL(NOTARGETS)
      DRESPWM  =   DRESPWM / REAL(NOTARGETS)
      DRESPWG  =   DRESPWG / REAL(NOTARGETS)
      DRESPBM  =   DRESPBM / REAL(NOTARGETS)
      DRESPBG  =   DRESPBG / REAL(NOTARGETS)
      DRESPCRM =   DRESPCRM / REAL(NOTARGETS)
      DRESPFRM =   DRESPFRM / REAL(NOTARGETS)
      DRESPFRG =   DRESPFRG / REAL(NOTARGETS)
      DRESPCRG =   DRESPCRG / REAL(NOTARGETS)
      DRESPFG  =   DRESPFG  / REAL(NOTARGETS)

      RETURN
      END ! PerTree
