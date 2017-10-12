$DEBUG:'D' !DEBUG VERSION; ALL LINES BEGINNING WITH 'D' ARE COMPILED IN

C**********************************************************************
c INOUT.FOR
C This file contains all the subroutines for reading and writing to
C input and output files.
C The main subroutines (called externally) are:
C OPENINPUTF - opens the input files
C OPENOUTPUTF - opens the output files
C CLOSEF - closes the files
C INPUTSTR - reads the canopy structure file str.dat
C INPUTPHY - reads the physiology file phy.dat
C INPUTTREE - reads the trees.dat file
C INPUTCON - reads the control file confile.dat
C INPUTSOIL - reads the soil/understorey file soil.dat
C OUTPUTDY - outputs daily fluxes
C OUTPUTHR - outputs hourly fluxes
C OUTPUTLAY - outputs to layer flux file
C OUTPUTHIST - outputs PAR histogram
C SORTTREES - sorts the trees into order of distance from the target tree
C INTERPOLATEP - calls the daily interpolation routines for physiology
C INTERPOLATET - calls the daily interpolation routines for tree dimensions
C NB The subroutines to do with reading the met file are in getmet.for.

C Subsidiary subroutines are:
C INPUTSTR
C   READCROWN - read crown shape 
C   READAERO - read aerodynamic properties (e.g. extinction of wind) 
C   READBETA - read leaf area density distribution
C   READALLOM - read parameters for biomass vs height & diameter
C   READLIA - read leaf incidence angle 
C READLIA uses the following
C   ANGLE - calculate the fraction of leaf area in each angle class
C   AVGLIA - calculate the average LIA for a given ellipsoidal parameter
C   INTEG - integrate leaf angle dist over specified range
C   FANG - ellipsoidal leaf angle density function
C   CALCELP - calculate ellipsoidal parameter from average LIA
C INPUTCON
C   READDATES - read in start & end dates of simulation
C   READMODEL - read in options for which model to use
C   READZEN - read in number of angle classes & layers to use
C   READHIST - read in details of PAR histogram required
C   READCCSCEN - read in climate change scenario
C   READOTC - read in effects of OTC on met data
C INPUTPHY
c   READAGEP - read number of age classes for physiology parameters
C   READPROP - read proportions of leaf area in each age class
C   READABSRP - get leaf absorptance, reflectance & transmittance arrays
C   READGS - read parameters of stomatal model
C   READLEAFN - read in leaf nitrogen concentrations
C   READJMAX - read in Jmax and Vcmax values
C   READRD - read in parameters of leaf respiration
C   READPHYARRAY - read array of physiological parameters 
C     (called by READLEAFN, READJMAX, READRD)
C   READRW - read in parameters for woody respiration
C   READRR - read in parameters for root respiration
C   READRB - read in parameters for branch respiration
C   RINTEG - convert physiology parameters spec. in layers to layers used
C     (called by READBASRP, READPHYARRAY)
C INPUTTREE
C   READPLOT - read dimensions of plot
C   READXYZ - read co-ordinates of trees
C   READTREEARRAY - read arrays of tree dimensions eg height, diameter
C   GETLEAFAREA - A subroutine to read in leaf area array.
C   CALCLAI - calculate plot LAI
C   READZPD - read canopy dimensions (roughness length etc)
C   READCONTREES - get list of target trees
C   INEDGES - determine whether given tree is in plot edges or not
C INTERPOLATEP
C   PHYINTERP - do interpolation of physiological parameters
C INTERPOLATET
C   PHENOL - calculate leaf area from phenology parameters
C   TREEINTERP - do interpolation of tree dimensions
C Unused: NOTREEGET, GETTREENO, READRSOIL
C**********************************************************************


C**********************************************************************
      SUBROUTINE OPENINPUTF(IODAILYI,IOHRLYI,IOTUTDI,IOHISTI,IORESPI,
     &  IOSOIL,CTITLE,TTITLE,PTITLE,STITLE,UTITLE)
C This routine opens the input files.
C The filenames are defined in this routine.
C**********************************************************************

      INCLUDE 'maestcom'

      CHARACTER*80 CTITLE, TTITLE, PTITLE, STITLE, UTITLE

      NAMELIST /CONTROL/ IOHRLY,IOTUTD,IOHIST,IORESP

990   FORMAT (A60)     ! For reading titles in input files.

C Output file for errors and warnings
      OPEN (UERROR, FILE = 'Maeserr.dat', STATUS = 'UNKNOWN')

C Input file with control switches
      OPEN (UCONTROL, FILE = 'confile.dat', STATUS = 'OLD')
C Default values of control switches
      IOHRLY = 0   ! Controls daily, hourly, and/or layer output
      IOTUTD = 1   ! Controls transmittance file output
      IOHIST = 0   ! Controls histogram output
      IORESP = 0   ! Controls respiration output
      IODAILY = 1  ! Controls daily output: FIXED HERE 
	IOSOIL = 1   ! Calculate water balance? Will depend on whether file exists   

C Read control file
      READ (UCONTROL, 990) CTITLE
      READ (UCONTROL, CONTROL, IOSTAT = IOERROR)
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('WARNING: USING DEFAULT VALUES FOR CONTROL FILE',
     &  IWARN,IOERROR)

      IOHRLYI = IOHRLY
      IOTUTDI = IOTUTD
      IOHISTI = IOHIST
      IORESPI = IORESP
      IODAILYI = IODAILY

C Input file with data on tree position and size
      OPEN (UTREES, FILE = 'trees.dat', STATUS='OLD')
C Input file with data on canopy structure
      OPEN (USTR, FILE = 'str.dat', STATUS='OLD')
C Input file with data on physiology
      OPEN (UPHY, FILE = 'phy.dat', STATUS='OLD')
C Input file with data on soil and understory
	OPEN (USOILI, FILE = 'soil.dat', STATUS='OLD',IOSTAT = IOERROR)
	IF (IOERROR.NE.0) IOSOIL = 0  ! If file doesn't exist, don't do water balance

C Input/output file with diffuse transmittances
      OPEN (UTUTD, FILE = 'tutd.dat', STATUS='UNKNOWN')

C Read titles from input files
      READ (UTREES, 990) TTITLE
      READ (USTR, 990) STITLE
      READ (UPHY, 990) PTITLE
	IF (IOSOIL.EQ.1) READ (USOILI, 990) UTITLE

      RETURN
      END !OpenInputF


C**********************************************************************
      SUBROUTINE OPENOUTPUTF(IODAILY,IOHRLY,IOHIST,IORESP,IOSOIL,
     &  CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE)
C This routine opens the output files.
C The filenames are defined in this routine.
C It also writes initial comments to the output files.
C**********************************************************************

      INCLUDE 'maestcom'

      CHARACTER*80 CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE

990   FORMAT (A80)
991   FORMAT (A12,A60) ! For writing comments to output files.

C Output file with daily fluxes
      IF (IODAILY.GT.0) THEN
        OPEN (UDAILY,FILE = 'Dayflx.dat', STATUS='UNKNOWN')
      END IF

C Output file with hourly fluxes (if required).
      IF (IOHRLY.GT.0) THEN
        OPEN (UHRLY, FILE = 'hrflux.dat', STATUS='UNKNOWN')
      ENDIF
C Output file with layer fluxes (if required).
      IF (IOHRLY.GT.1) THEN
        OPEN (ULAY, FILE = 'layflx.dat', STATUS='UNKNOWN')
      ENDIF

C Output file with histogram (if required).
      IF (IOHIST.EQ.1) THEN
        OPEN (UHIST, FILE = 'histo.dat', STATUS = 'UNKNOWN')
      END IF

C Output file for respiration (if required).
      IF (IORESP.EQ.1) THEN
        OPEN (URESP, FILE = 'resp.dat', STATUS = 'UNKNOWN')
      END IF

C Output file for soil water balance (if required).
      IF (IOSOIL.EQ.1) THEN
        OPEN (USOILO, FILE = 'soilout.dat', STATUS = 'UNKNOWN')
      END IF

C Write headings to daily flux files.
      IF (IODAILY.GT.0) THEN
        WRITE (UDAILY, 991) 'Program:    ', VTITLE
        WRITE (UDAILY, 991) 'Control:    ', CTITLE
        WRITE (UDAILY, 991) 'Trees:      ', TTITLE
        WRITE (UDAILY, 991) 'Structure:  ', STITLE
        WRITE (UDAILY, 991) 'Physiology: ', PTITLE
        WRITE (UDAILY, 991) 'Met data:   ', MTITLE
        IF (IOSOIL.EQ.1) THEN
	    WRITE (UDAILY, 991) 'Soil data:  ', UTITLE
	  END IF
        WRITE (UDAILY, 990) ' '
        WRITE (UDAILY,501)
        WRITE (UDAILY,502)
        WRITE (UDAILY,503)
        WRITE (UDAILY,504)
        WRITE (UDAILY,505)
        WRITE (UDAILY,506)
        WRITE (UDAILY,507)
        WRITE (UDAILY,508)
        WRITE (UDAILY,509)
        WRITE (UDAILY,510)
        WRITE (UDAILY,511)
        WRITE (UDAILY, 990) ' '
        WRITE (UDAILY,512)
      END IF

C Comments to hourly output file (if required).
      if (IOHRLY.gt.0) then
        WRITE (UHRLY, 991) 'Program:    ', VTITLE
        WRITE (UHRLY, 991) 'Control:    ', CTITLE
        WRITE (UHRLY, 991) 'Trees:      ', TTITLE
        WRITE (UHRLY, 991) 'Structure:  ', STITLE
        WRITE (UHRLY, 991) 'Physiology: ', PTITLE
        WRITE (UHRLY, 991) 'Met data:   ', MTITLE
        WRITE (UHRLY, 990) ' '
        WRITE (UHRLY,701)
        WRITE (UHRLY,702)
        WRITE (UHRLY,703)
        WRITE (UHRLY,704)
        WRITE (UHRLY,705)
        WRITE (UHRLY,706)
        WRITE (UHRLY,707)
        WRITE (UHRLY,708)
        WRITE (UHRLY,709)
        WRITE (UHRLY,710)
        WRITE (UHRLY,711)
        WRITE (UHRLY,712)
        WRITE (UHRLY,713)
        WRITE (UHRLY,714)
        WRITE (UHRLY,715)
        WRITE (UHRLY, 990) ' '
        WRITE (UHRLY,716)
      END IF

C Comments to respiration output file (if required).
      if (IORESP.gt.0) then
        WRITE (URESP, 991) 'Program:    ', VTITLE
        WRITE (URESP, 991) 'Control:    ', CTITLE
        WRITE (URESP, 991) 'Trees:      ', TTITLE
        WRITE (URESP, 991) 'Structure:  ', STITLE
        WRITE (URESP, 991) 'Physiology: ', PTITLE
        WRITE (URESP, 991) 'Met data:   ', MTITLE
        WRITE (URESP, 990) ' '
        WRITE (URESP,601)
        WRITE (URESP,602)
        WRITE (URESP,603)
        WRITE (URESP,604)
        WRITE (URESP,605)
        WRITE (URESP,606)
        WRITE (URESP,607)
        WRITE (URESP,608)
        WRITE (URESP,609)
        WRITE (URESP,610)
        WRITE (URESP,611)
        WRITE (URESP, 990) ' '
        WRITE (URESP,612)
      END IF

C Write comments to layer output file (if required)
      IF (IOHRLY.GT.1) THEN
        WRITE (ULAY, 991) 'Program:    ', VTITLE
        WRITE (ULAY, 801)
        WRITE (ULAY, 802)
        WRITE (ULAY, 803)
        WRITE (ULAY, 804)
        WRITE (ULAY, *)
      END IF

C Write comments to soil water output file (if required)
      IF (IOSOIL.EQ.1) THEN
        WRITE (USOILO, 991) 'Program:    ', VTITLE
        WRITE (USOILO, 991) 'Control:    ', CTITLE
        WRITE (USOILO, 991) 'Trees:      ', TTITLE
        WRITE (USOILO, 991) 'Structure:  ', STITLE
        WRITE (USOILO, 991) 'Physiology: ', PTITLE
        WRITE (USOILO, 991) 'Met data:   ', MTITLE
        WRITE (USOILO, 991) 'Soil data:  ', UTITLE
	  WRITE (USOILO, *)
	  WRITE (USOILO, 401)
	  WRITE (USOILO, 402)
	  WRITE (USOILO, 403)
	  WRITE (USOILO, 404)
	  WRITE (USOILO, 405)
	  WRITE (USOILO, 406)
	  WRITE (USOILO, 407)
	  WRITE (USOILO, *)
	  WRITE (USOILO, 408)
	END IF

501   format('DOY: simulation date')
502   FORMAT('Tree: tree number')
503   format('absPAR:   absorbed PAR              MJ tree-1 d-1')
504   format('absNIR:   absorbed NIR              MJ  tree-1 d-1')
505   format('absTherm: absorbed thermal          MJ tree-1 d-1')
506   format('totPs: gross photosynthesis         mol tree-1 d-1')
507   format('totRf: daily foliar respiration     mol tree-1 d-1')
508   format('netPs: photosyn. net of foliar resp   mol tree-1 d-1')
509   format('totLE1: daily transpiration         mol H2O tree-1 d-1')
510   format('totLE2: daily transpirn (CANOPY calc) mol H2O m-2 d-1')
511   format('totH:  daily sensible heat flux     MJ tree-1 d-1')
512   format('Columns: DOY, Tree, absPAR, absNIR, absTherm,',
     &  ' totPs, totRf, netPs,',
     +  ' totLE1, totLE2, totH')

701   format('DOY: simulation date')
702   FORMAT('Tree: tree number')
703   format('Hour:  hour of the day')
704   format('hrPAR: absorbed PAR              umol tree-1 s-1')
705   format('hrNIR: absorbed NIR              W tree-1')
706   format('hrTHM: absorbed thermal          W tree-1')
707   format('hrPS: photosynthesis (net of leaf resp) umol tree-1 s-1')
708   format('hrRf:  hourly leaf respiration   umol tree-1 s-1')
709   format('hrRmW: hourly stem + branch Rm   umol tree-1 s-1')
710   format('hrLE:  hourly transpiration      mmol tree-1 s-1')
711   format('hrEV:  hourly evaporation      mmol tree-1 s-1')
712   format('LECAN: hourly transpirn: CANOPY calc : mmol H2O m-2 s-1')
713   format('Gscan: canopy stomatal conductance : mol CO2 tree-1 s-1')
714   format('hrH:   hourly sensible heat flux:  MJ tree-1 s-1')
715	FORMAT('TCAN: Average foliage temperature (deg C)')
716   format('Columns: DOY, Tree, HOUR, hrPAR, hrNIR, hrTHM,',
     +  ' hrPs, hrRf, hrRmW, hrLE, hrEV,',
     +  ' LECAN, Gscan, hrH, TCAN')

801   FORMAT(' Fluxes for each layer on an hourly basis')
802   FORMAT(' Rows: absorbed PAR (umol m-2 leaf s-1) ')
803   FORMAT('       photosynthesis net of Rleaf (umol m-2 leaf s-1) ')
804   FORMAT('       transpiration (umol m-2 leaf s-1) ')

601   FORMAT('Daily maintenance and growth respiration components')
602   FORMAT('Rmf: Foliage maintenance resp.    mol tree-1 d-1')
603   FORMAT('Rmw: Stem maintenance resp.       mol tree-1 d-1')
604   FORMAT('RmB: Branch maintenance resp.     mol tree-1 d-1')
605   FORMAT('Rmcr: Coarse root maintenance resp. mol tree-1 d-1')
606   FORMAT('Rmfr: Fine root maintenance resp. mol tree-1 d-1')
607   FORMAT('Rgf: Foliage growth resp.         mol tree-1 d-1')
608   FORMAT('Rgw: Stem growth resp.            mol tree-1 d-1')
609   FORMAT('Rgb: Branch growth resp.          mol tree-1 d-1')
610   FORMAT('Rgcr: Coarse root growth resp.    mol tree-1 d-1')
611   FORMAT('Rgfr: Fine root growth resp.      mol tree-1 d-1')
612   FORMAT('Columns: DOY, Tree, Rmf, Rmw, Rmb, Rmcr, Rmfr, Rgf, Rgw, 
     &  Rgb, Rgcr, Rgfr')

401   FORMAT('Daily soil water balance components')
402   FORMAT('PPT: Daily precipitation		mm')
403   FORMAT('ET:  Evaporation/Transpiration  mol tree-1 d-1')
404   FORMAT('DR:  Drainage from soil			mm')
405   FORMAT('CAN: Canopy water storage at end of day	mm')
406   FORMAT('ST:  Stem water storage at end of day	mm')
407   FORMAT('SO:  Total soil water storage at end of day	%vol')
408   FORMAT('Columns: DOY, PPT, ET, DR, CAN, ST, SO')

      RETURN
      END !OpenOutputF


C**********************************************************************
      SUBROUTINE CLOSEF(IODAILY,IOHRLY,IOHIST,IORESP,IOSOIL)
C This routine closes the open files.
C**********************************************************************

      INCLUDE 'maestcom'

      CLOSE(UCONTROL)
      CLOSE(UTREES)
      CLOSE(USTR)
      CLOSE(UPHY)
      CLOSE(UMET)
      CLOSE(UTUTD)
      CLOSE(UERROR)
      IF (IODAILY.GT.0) CLOSE(UDAILY)
      IF (IOHRLY.GT.0) CLOSE(UHRLY)
      IF (IOHRLY.GT.1) CLOSE(ULAY)
      IF (IOHIST.EQ.1) CLOSE(UHIST)
      IF (IORESP.EQ.1) CLOSE(URESP)
	IF (IOSOIL.EQ.1) CLOSE(USOILI)
	IF (IOSOIL.EQ.1) CLOSE(USOILO)

      RETURN
      END !Closef


C**********************************************************************
      SUBROUTINE INPUTCON(
     &  ISTART, IEND, NSTEP,
     &  NUMPNT,NOLAY,NZEN,DIFZEN,NAZ,
     &  MODELGS, MODELJM, MODELRD, MODELSS, MODELRW, ITERMAX,
     &  IOHIST, BINSIZE,
     &  ICC, CO2INC, TINC,
     &  IOTC, TOTC, WINDOTC, PAROTC, FBEAMOTC
     &  )
C Read in the information from the control file. 
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL DIFZEN(MAXANG)

      CALL READDATES(UCONTROL, ISTART, IEND, NSTEP)

      CALL READZEN(UCONTROL, NUMPNT, NOLAY, NZEN, NAZ, DIFZEN)

      CALL READMODEL(UCONTROL, MODELGS, MODELJM, MODELRD, 
     &  MODELSS, MODELRW, ITERMAX)

      CALL READHIST(UCONTROL,IOHIST,BINSIZE)

      CALL READCCSCEN(UCONTROL,ICC,CO2INC,TINC)

      CALL READOTC(UCONTROL,IOTC,TOTC,WINDOTC,PAROTC,FBEAMOTC)

      RETURN
      END  ! InputCon


C**********************************************************************
      SUBROUTINE InputStr(
     &  JLEAF,BPT,RANDOM,NOAGEC,
     &  JSHAPE,SHAPE,EXTWIND,
     &  NALPHA,ALPHA,FALPHA,
     &  COEFFT,EXPONT,WINTERC,
     &  BCOEFFT,BEXPONT,BINTERC,
     &  RCOEFFT,REXPONT,RINTERC,FRFRAC,
     &  PFLA,CANCAPLA
     &  )
C This routine reads in input data on canopy structure (from USTR)
C which is required to begin the simulation.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL ALPHA(MAXANG),FALPHA(MAXANG)
      REAL BPT(8,MAXC)

      CALL READCROWN(USTR, JSHAPE, SHAPE)

      CALL READAERO(USTR, EXTWIND)

      CALL READLIA(USTR, NALPHA, ALPHA, FALPHA)

      CALL READBETA(USTR, NOAGEC, JLEAF, BPT, RANDOM)

      CALL READALLOM(USTR, COEFFT, EXPONT, WINTERC,
     &  BCOEFFT, BEXPONT, BINTERC, 
     &  RCOEFFT, REXPONT, RINTERC, FRFRAC)

      CALL READPPT(USTR,PFLA,CANCAPLA)

      RETURN
      END !InputStr


C**********************************************************************
      SUBROUTINE InputPhy(
     &  MODELJM,MODELRD,MODELGS,MODELRW,
     &  NOLAY,NOAGEC,NOAGEP,PROP,
     &  ABSRP,REFLEC,TRANS,RHOSOL,
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
C This routine reads in input data on physiology (from UPHY).
C Some parameters can change with time: JMAX25, VCMAX25, RD0
C - for these, arrays of dates and values are read in, and must
c be interpolated.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL ABSRP(MAXLAY,3),REFLEC(MAXLAY,3),TRANS(MAXLAY,3)
      REAL RHOSOL(3), PROP(MAXC)
      REAL LEAFN(MAXPDATE,MAXLAY,MAXC)
      REAL JMAXTABLE(MAXPDATE,MAXLAY,MAXC)
      REAL VCMAXTABLE(MAXPDATE,MAXLAY,MAXC)
      REAL RDTABLE(MAXPDATE,MAXLAY,MAXC)
      REAL SLATABLE(MAXPDATE,MAXLAY,MAXC)
	REAL AJQTABLE(MAXPDATE,MAXLAY,MAXC)
      INTEGER DATESN(MAXPDATE), DATESJ(MAXPDATE)
      INTEGER DATESV(MAXPDATE), DATESRD(MAXPDATE)
      INTEGER DATESSLA(MAXPDATE), DATESA(MAXPDATE)

      CALL READAGEP(UPHY, NOAGEC, NOAGEP)

      CALL READPROP(UPHY, NOAGEC, NOAGEP, PROP)

      CALL READABSRP(UPHY,NOLAY,ABSRP,REFLEC,TRANS,RHOSOL)

      CALL READGS(UPHY,MODELGS,
     &  GSREF, GSMIN, PAR0, D0, VK1, VK2, VPD1, VPD2, VMFD0, 
     &  GSJA, GSJB, T0, TREF, TMAX, SMD1, SMD2,
     &  G0, D0L, GAMMA, G1, WLEAF, NSIDES)

      CALL READLEAFN(UPHY, MODELJM, MODELRD, NOLAY, NOAGEP,
     &  NONDATES, DATESN, LEAFN)

      CALL READJMAX(UPHY, MODELJM, NOLAY, NOAGEP, NONDATES,
     &  DATESN, LEAFN, NOJDATES, DATESJ, JMAXTABLE,
     &  NOVDATES, DATESV, VCMAXTABLE,
     &  NOADATES, DATESA, AJQTABLE,
     &  IECO, EAVJ, EDVJ, DELSJ, EAVC, EDVC, DELSC, TVJUP, TVJDN, THETA)

C Read in SLA array directly.
      CALL READPHYARRAY(UPHY,5,NOLAY,NOAGEP,
     &    NOSLADATES,DATESSLA,SLATABLE)

      CALL READRD(UPHY, MODELRD, NOLAY, NOAGEP, NONDATES,
     &  DATESN, LEAFN, NORDATES, DATESRD, RDTABLE,
     &  Q10F, RTEMP, DAYRESP, EFFYRF, TBELOW)

      CALL READRW(UPHY,MODELRW,EFFYRW,RMW,RTEMPW,Q10W,
     &  COLLA,COLLK,STEMSDW,RMWAREA,STEMFORM)

      CALL READRR(UPHY,RMCR,RMFR,Q10R,RTEMPR)
      CALL READRB(UPHY,RMB,Q10B,RTEMPB)

      RETURN
      END !InputPhy

		
C**********************************************************************
      SUBROUTINE INPUTTREE(
     &  XSLOPE,YSLOPE,BEAR,XMAX,YMAX,STOCKING,
     &  ZHT,Z0HT,ZPD,
     &  NOALLTREES,NOTREES,NOTARGETS,ITARGETS,SHADEHT,
     &  NOXDATES,NOYDATES,NOZDATES,NOTDATES,NOLADATES,NODDATES,
     &  DATESX,DATESY,DATESZ,DATEST,DATESLA,DATESD,
     &  DX,DY,DZ,R1,R2,R3,TRUNK,FLT,TOTLAI,DIAMA,
     &  IFLUSH,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN
     &  )
C This subroutine should read in the data from the trees.dat file on
C crown positions and dimensions. Some variables can change with time:
C radii, height, diameter & leaf area - for these, arrays of dates & values at
C those dates may be read in, & interpolated during the program.
C x- and y- co-ordinates of the crowns may be read in or calculated from
C plot size and stocking density.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DX(MAXT),DY(MAXT),DZ(MAXT)
      REAL R1(MAXTDATE,MAXT),R2(MAXTDATE,MAXT),R3(MAXTDATE,MAXT)
      REAL TRUNK(MAXTDATE,MAXT),FLT(MAXTDATE,MAXT),TOTLAI(MAXTDATE)
      REAL DIAMA(MAXTDATE,MAXT)
      INTEGER DATESX(MAXTDATE),DATESY(MAXTDATE),DATESZ(MAXTDATE)
      INTEGER DATEST(MAXTDATE),DATESLA(MAXTDATE),DATESD(MAXTDATE)
      INTEGER ITARGETS(MAXT)

C Read in number of trees & number of target tree
      CALL READPLOT(UTREES, X0, Y0, XMAX, YMAX, NOALLTREES,
     &  XSLOPE, YSLOPE, BEAR, SHADEHT)
      STOCKING = NOALLTREES/(XMAX*YMAX)

C Read in aerodynamic properties of canopy
      CALL READZPD(UTREES,ZHT,Z0HT,ZPD)

C Get x, y, z co-ords of each tree
      CALL READXYZ(UTREES,NOALLTREES,X0,Y0,XMAX,YMAX,XSLOPE,YSLOPE,
     &  DX,DY,DZ)

C Get radii in x & y directions of each tree
      CALL READTREEARRAY(UTREES,1,NOALLTREES,NOXDATES,DATESX,R1)
      CALL READTREEARRAY(UTREES,2,NOALLTREES,NOYDATES,DATESY,R2)
C Get green crown height of each tree
      CALL READTREEARRAY(UTREES,3,NOALLTREES,NOZDATES,DATESZ,R3)
C Get trunk length of each tree
      CALL READTREEARRAY(UTREES,4,NOALLTREES,NOTDATES,DATEST,TRUNK)

C Get leaf area parameters 
      CALL GETLEAFAREA(UTREES,IFLUSH,DT1,DT2,DT3,DT4,
     &  EXPTIME,APP,EXPAN,NOALLTREES,NOLADATES,DATESLA,FLT)

C Get diameter of each tree
      CALL READTREEARRAY(UTREES,6,NOALLTREES,NODDATES,DATESD,DIAMA)

C Calculate total LAI
      CALL CALCLAI(NOLADATES,FLT,NOALLTREES,XMAX,YMAX,XSLOPE,YSLOPE,
     &  TOTLAI)

C Read in how many of the trees form the subplot (from confile)
      CALL READCONTREES(UCONTROL,NOALLTREES,DX,DY,XMAX,YMAX,
     &  NOTREES,NOTARGETS,ITARGETS)

      RETURN
      END !InputTree


C**********************************************************************
      SUBROUTINE SORTTREES(
     &  NOALLTREES,NOTREES,NOTARGET,
     &  DX,DY,DZ,R1,R2,R3,TRUNK,FLT,DIAMA,
     &  DXT,DYT,DZT,RX,RY,RZ,FOLT,ZBC,DIAM
     &  )
C This routine selects the 'NOTREES' trees closest to the target tree.
C It sorts the information about each tree into order of distance from
C the target tree. Does this simply by calling SORTTREESP.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DX(MAXT),DY(MAXT),DZ(MAXT)
      REAL R1(MAXTDATE,MAXT),R2(MAXTDATE,MAXT),R3(MAXTDATE,MAXT)
      REAL TRUNK(MAXTDATE,MAXT),FLT(MAXTDATE,MAXT)
      REAL DXT(MAXTT),DYT(MAXTT),DZT(MAXTT)
      REAL RX(MAXTDATE,MAXTT),RY(MAXTDATE,MAXTT),RZ(MAXTDATE,MAXTT)
      REAL ZBC(MAXTDATE,MAXTT),FOLT(MAXTDATE,MAXTT)
      REAL DIAMA(MAXTDATE,MAXT),DIAM(MAXTDATE,MAXTT)

      X = DX(NOTARGET)
      Y = DY(NOTARGET)
      CALL SORTTREESP(
     &  X,Y,NOALLTREES,NOTREES,
     &  DX,DY,DZ,R1,R2,R3,TRUNK,FLT,DIAMA,
     &  DXT,DYT,DZT,RX,RY,RZ,FOLT,ZBC,DIAM
     &  )

      RETURN
      END !SortTrees


C**********************************************************************
      SUBROUTINE SORTTREESP(
     &  X,Y,NOALLTREES,NOTREES,
     &  DX,DY,DZ,R1,R2,R3,TRUNK,FLT,DIAMA,
     &  DXT,DYT,DZT,RX,RY,RZ,FOLT,ZBC,DIAM
     &  )
C This routine selects the 'NOTREES' trees closest to the point (x,y,z).
C It sorts the information about each tree into order of distance from
C this point.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL DX(MAXT),DY(MAXT),DZ(MAXT)
      REAL R1(MAXTDATE,MAXT),R2(MAXTDATE,MAXT),R3(MAXTDATE,MAXT)
      REAL TRUNK(MAXTDATE,MAXT),FLT(MAXTDATE,MAXT)
      REAL DXT(MAXTT),DYT(MAXTT),DZT(MAXTT)
      REAL RX(MAXTDATE,MAXTT),RY(MAXTDATE,MAXTT),RZ(MAXTDATE,MAXTT)
      REAL ZBC(MAXTDATE,MAXTT),FOLT(MAXTDATE,MAXTT)
      REAL DIAMA(MAXTDATE,MAXT),DIAM(MAXTDATE,MAXTT)
      REAL DNT(MAXT)
      INTEGER IT(MAXT)

C We select NT trees which are closest to the tree (NOTARGET) we
C are concerned with.
      DO 10 N = 1,NOALLTREES
        DNT(N) = SQRT((X-DX(N))**2 + (Y-DY(N))**2)
        IT(N) = N
10    CONTINUE

C Perform a sort on the trees to find the NOTREES closest ones to NOTARGET.
      MT1 = NOALLTREES - 1
      DO 20 I = 1,MT1
         MT2 = I + 1
         DO 20 J = MT2,NOALLTREES
            IF (DNT(I).LE.DNT(J)) GO TO 20
            TEMP = DNT(I)
            DNT(I) = DNT(J)
            DNT(J) = TEMP
            ITEMP = IT(I)
            IT(I) = IT(J)
            IT(J) = ITEMP
20    CONTINUE

C Produce new arrays containing the information about the closest trees.
      DO 30 I = 1,NOTREES
        MM = IT(I)
        DXT(I) = DX(MM)
        DYT(I) = DY(MM)
        DZT(I) = DZ(MM)
        DO 30 IDATE = 1,MAXTDATE
          RX(IDATE,I) =  R1(IDATE,MM)
          RY(IDATE,I) = R2(IDATE,MM)
          RZ(IDATE,I) = R3(IDATE,MM)
          ZBC(IDATE,I) = TRUNK(IDATE,MM)
          FOLT(IDATE,I) = FLT(IDATE,MM)
          DIAM(IDATE,I) = DIAMA(IDATE,MM)
30    CONTINUE

      RETURN
      END !SortTrees


C**********************************************************************
      SUBROUTINE ANGLE(ELP,NALPHA,FALPHA)
C This routine calculates the fraction of leaf area in each leaf
C angle class.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL FALPHA(MAXANG)

      DALP = PID2/FLOAT(NALPHA)
      IF (ELP.NE.1) THEN
         IF (ELP.LT.1) THEN
            ECTR = SQRT(1.0-ELP*ELP)
            TEMP = 1.0 + ASIN(ECTR)/ (ECTR*ELP)
         END IF

         IF (ELP.GT.1) THEN
            ECTR = SQRT(1.0- (1.0/ELP)**2)
            TEMP = 1.0 + LOG((1.0+ECTR)/ (1.0-ECTR))/ (2.0*ELP*ELP*ECTR)
         END IF

         TOTAL = 0.000
         DO 100 N = 1,NALPHA
            ALP1 = (N-1)*DALP
            ALP2 = N*DALP
            SUM = 0.000
            CALL INTEG(ALP1,ALP2,ELP,TEMP,SUM)
            FALPHA(N) = SUM
            TOTAL = TOTAL + SUM
  100    CONTINUE
         DO 110 N = 1,NALPHA
            FALPHA(N) = FALPHA(N)/TOTAL
  110    CONTINUE
      ELSE
         DO 200 N = 1,NALPHA
            FALPHA(N) = COS((N-1)*DALP) - COS(N*DALP)
  200    CONTINUE
      END IF

      RETURN
      END ! Angle

C**********************************************************************
      SUBROUTINE INTEG(ALP1,ALP2,ELP,TEMP,SUM)
C     this is subroutine to integrate the leaf area density
C     function over a specified angle range
C**********************************************************************

      H = (ALP2-ALP1)/50.0
      UP = ALP1 + H*0.42265
      DOWN = ALP1 + H*1.57735
      SUM = FANG(ELP,TEMP,UP) + FANG(ELP,TEMP,DOWN)
      DO 100 I = 1,24
         UP = UP + 2.0*H
         DOWN = DOWN + 2.0*H
  100 SUM = SUM + FANG(ELP,TEMP,UP) + FANG(ELP,TEMP,DOWN)
      SUM = H*SUM

      RETURN
      END !Integ


C**********************************************************************
      FUNCTION AVGLIA(ELP)
C     this function is the relation between the parameter ELP and the
C     average leaf angle. (see notes for details)
C**********************************************************************

      IF (ELP.GT.1.0) THEN
         AVGLIA = 1.0/ (0.5901*ELP+0.3037)
      ELSE
         AVGLIA = 1.0/ (0.3782*ELP+0.6131)
      END IF

      RETURN
      END  !AVGLIA


C**********************************************************************
      FUNCTION CALCELP(AVGANG)
C     this function is the relation between the average leaf angle and
C     the parameter ELP. This is the inverse of function AVGLIA.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      AVGANG = AVGANG*PID180    !Convert to radians
      IF (AVGANG.LT.1.0) THEN
         CALCELP = (1./AVGANG - 0.3037)/0.5901
      ELSE
         CALCELP = (1./AVGANG - 0.6131)/0.3782
      END IF

      RETURN
      END  !CalcElp


C**********************************************************************
      FUNCTION FANG(ELP,TEMP,ALP)
C     this is the ellipsoidal leaf angular density function
C     reference: Campbell,G.S. (personal comm.)
C**********************************************************************

      FANG = 2.0*ELP*ELP*SIN(ALP)/ (TEMP*
     +       (((COS(ALP))**2.0+ (ELP*SIN(ALP))**2)**2))

      RETURN
      END  !Fang



C**********************************************************************
      SUBROUTINE OUTPUTHR(IOHRLY,IDAY,IHOUR,ITREE,TCAN,
     &  NOLAY,PPAR,PPS,PTRANSP,FOLLAY,
     &  THRAB,FCO2,FRESPF,FRESPW,FRESPB,
     &  FH2OT,FH2OE,GSCAN,FH2OCAN,FHEAT)
C Output the hourly totals
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL THRAB(KHRS,3),FCO2(KHRS),FRESPF(KHRS)
      REAL FRESPW(KHRS),FRESPB(KHRS)
      REAL GSCAN(KHRS),FH2OT(KHRS),FH2OE(KHRS),FH2OCAN(KHRS),FHEAT(KHRS)
      REAL PPAR(MAXLAY,KHRS),PPS(MAXLAY,KHRS),PTRANSP(MAXLAY,KHRS)
      REAL FOLLAY(MAXLAY),TCAN(KHRS)

      IF (IOHRLY.GE.1) THEN
        WRITE (UHRLY,500) IDAY,ITREE,IHOUR,THRAB(IHOUR,1)*UMOLPERJ,
     &    THRAB(IHOUR,2),THRAB(IHOUR,3),
     &    FCO2(IHOUR),FRESPF(IHOUR),FRESPW(IHOUR)+FRESPB(IHOUR),
     &    FH2OT(IHOUR)*1e-3,FH2OE(IHOUR)*1e-3,
     &    FH2OCAN(IHOUR)*1E-3,GSCAN(IHOUR),
     &    FHEAT(IHOUR)*1E-3,TCAN(IHOUR)
      END IF
500   FORMAT (I7,1X,2(I4,1X),3(F12.5,1X),9(F12.7,1X))

      IF (IOHRLY.GE.2) THEN
        WRITE (ULAY,610) 'DAY',IDAY,'HOUR',IHOUR
        IF (FOLLAY(1).GT.0.0) THEN 
          WRITE (ULAY,600) (PPAR(I,IHOUR)/FOLLAY(I),I=1,NOLAY)
          WRITE (ULAY,600) (PPS(I,IHOUR)/FOLLAY(I),I=1,NOLAY)
          WRITE (ULAY,600) (PTRANSP(I,IHOUR)/FOLLAY(I),I=1,NOLAY)
        ELSE
          WRITE (ULAY,*) 'NO FOLIAGE AT THIS TIME'
        END IF
      END IF
600   FORMAT (10(F10.2,1X))
610   FORMAT (A5,I5,A5,I5)

      RETURN
      END !OutputHr


C**********************************************************************
      SUBROUTINE OUTPUTDY(IDAY,ITREE,IODAILY,TDYAB,TOTCO2,TOTRESPF,
     &  TOTRESPW,TOTRESPWG,TOTH2O,TOTH2OCAN,TOTHFX,
     &  IORESP,TOTRESPCR,TOTRESPFR,TOTRESPFRG,TOTRESPCRG,TOTRESPFG,
     &  TOTRESPB,TOTRESPBG)
C Output the daily totals
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL TDYAB(3)

      IF (IODAILY.EQ.1) THEN
        WRITE (UDAILY,500) IDAY,ITREE,TDYAB,TOTCO2+TOTRESPF,TOTRESPF,
     &    TOTCO2,TOTH2O,TOTH2OCAN,TOTHFX
500   FORMAT (I7,1X,I4,1X,9(F12.5,1X))
      END IF

      IF (IORESP.EQ.1) THEN
        WRITE (URESP,510) IDAY,ITREE,TOTRESPF,TOTRESPW,TOTRESPB,
     &    TOTRESPCR,TOTRESPFR,TOTRESPFG,TOTRESPWG,TOTRESPBG,
     &    TOTRESPCRG,TOTRESPFRG 
      END IF
510   FORMAT (I7,1X,I4,1X,10(F12.5,1X))

      RETURN
      END !OutputDy


C**********************************************************************
      SUBROUTINE READDATES(UFILE, ISTART, IEND, NSTEPI)
C The routine must return start and end dates for the simulation,
C in days-since-1950 format. The function IDATE50 converts a DD/MM/YY
C date to this format. The routine must also return NSTEP where the
C program is to use met data for every NSTEP'th day (default 1).
C**********************************************************************

      
      INTEGER UFILE
      CHARACTER*10 START,END
	NAMELIST /DATES/ START,END,NSTEP

      START = '01/01/50'
      END = '01/01/50'
      NSTEP = 1
      REWIND (UFILE)
      READ (UFILE, DATES, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING: MISSING DATES: USING ALL MET DATA',
     &    IWARN,IOERROR)
        ISTART = 0
        IEND = 0
      ELSE
        ISTART = IDATE50(START)
        IEND = IDATE50(END)
      END IF
      NSTEPI = NSTEP

      RETURN
      END !ReadDates


C**********************************************************************
      SUBROUTINE READPLOT(UFILE, X0I, Y0I, XMAXI, YMAXI, NOALLTREESI,
     &  XSLOPEI, YSLOPEI, BEARI, SHADEHTI)
C Read in plot details. Subroutine must return:
C XMAX, YMAX - dimensions of plot, in m
c NOALLTREES - no of trees in plot
C XSLOPE, YSLOPE - slope of plot, in radians
C SHADEHT - height of shadecloth surrounding plot, if any
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /PLOT/ X0,Y0,XMAX,YMAX,NOTREES,
     &  XSLOPE,YSLOPE,BEARING,SHADEHT
      INTEGER UFILE

      SHADEHT = 0.0
	X0 = 0.0
	Y0 = 0.0

      REWIND (UFILE)
      READ (UFILE, PLOT, IOSTAT = IOERROR)
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('ERROR READING PLOT DETAILS',IFATAL,IOERROR)
	X0I = X0
	Y0I = Y0
      XMAXI = XMAX
      YMAXI = YMAX
      XSLOPEI = XSLOPE*PID180
      YSLOPEI = YSLOPE*PID180
      BEARI = BEARING*PID180
      SHADEHTI = SHADEHT

      IF (NOTREES.GT.MAXT) THEN
        CALL SUBERROR(
     &  'WARNING: NO OF TREES IN TREES FILE EXCEEDED MAXIMUM',
     &  IWARN,IOERROR)
        NOTREES = MAXT
      END IF
      NOALLTREESI = NOTREES

      RETURN
      END !ReadPlot


C**********************************************************************
      SUBROUTINE READZPD(UFILE,ZHTI,Z0HTI,ZPDI)
C Read in z, z0 and d for the canopy boundary layer conductance calcn.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /AERODYN/ ZHT,Z0HT,ZPD
      INTEGER UFILE

      Z0HT = 0.0
      REWIND (UFILE)
      READ (UFILE, AERODYN, IOSTAT = IOERROR)
      IF ((IOERROR.NE.0).OR.(Z0HT.EQ.0.0)) THEN
        CALL SUBERROR('WARNING: NOT CALCULATING CANOPY BDRY LAYER COND',
     &  IWARN,IOERROR)
      END IF
      ZHTI = ZHT
      Z0HTI = Z0HT
      ZPDI = ZPD

      RETURN
      END !ReadZPD
      

C**********************************************************************
      SUBROUTINE READCROWN(UFILE, JSHAPE, SHAPE)
C Read in crown shape parameters.
C JSHAPE indicates the shape of the crowns:
C   JCONE = conical,
C   JHELIP = half-ellipsoid,  (default)
C   JPARA = paraboloid,
C   JFELIP = full ellipsoid,
C   JCYL = cylinder,
C   JBOX = box.
C SHAPE is the factor needed to calculate volume for that crown shape:
C   VOLUME = PI * R2 * H * SHAPE (1/3, 2/3, 1/2, 2/3, 1, 4/PI).
C**********************************************************************

      INCLUDE 'MAESTCOM'
      
      INTEGER UFILE
      CHARACTER*5 CSHAPE
	NAMELIST /CANOPY/ CSHAPE

      CSHAPE = 'ELIP' !Default value
      REWIND (UFILE)
      READ (UFILE, CANOPY, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING: USING DEFAULT VALUE FOR CROWN SHAPE',
     &  IWARN,IOERROR)
      END IF

      IF (cshape.EQ.'CONE') THEN
         JSHAPE = JCONE
         SHAPE = 1.0/3.0

      ELSE IF (cshape.EQ.'ELIP') THEN
         JSHAPE = JHELIP
         SHAPE = 2.0/3.0

      ELSE IF (cshape.EQ.'PARA') THEN
         JSHAPE = JPARA
         SHAPE = 0.5

      ELSE IF (cshape.EQ.'ROUND') THEN
         JSHAPE = JFELIP
         SHAPE = 2.0/3.0

      ELSE IF (cshape.EQ.'CYL') THEN
         JSHAPE = JCYL
         SHAPE = 1.0

      ELSE IF (cshape.EQ.'BOX') THEN
         JSHAPE = JBOX
         SHAPE = 4.0/PI

      END IF

      RETURN
      END !ReadCrown


C**********************************************************************
      SUBROUTINE READAGEP(UFILE, NOAGEC, NOAGEPI)
C Age classes:
C NOAGEC is the number of age classes of foliage for which distributions
C of leaf area density are provided (read from str.dat).
C NOAGEP is the number of age classes
C of foliage for which (some) physiological parameters are provided.
C If neither NOAGEC nor NOAGEP = 1, then NOAGEC must = NOAGEP.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /NOAGES/ NOAGEP
      INTEGER UFILE

      NOAGEP = 1  !Default value

      REWIND (UFILE)
      READ (UFILE, NOAGES, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &  'WARNING: DEFAULT VALUE: NOAGEP=1',
     &  IWARN,IOERROR)
      END IF

C Check value does not exceed maximum.
      IF (NOAGEP.GT.MAXC) THEN
        CALL SUBERROR(
     &    'WARNING: NOAGEP EXCEEDED MAXIMUM AGE CLASSES',IWARN,0)
        NOAGEP = MAXC
      END IF

C Check same as no of structural age classes
      IF ((NOAGEC.NE.1).AND.(NOAGEP.NE.1).AND.(NOAGEC.NE.NOAGEP)) THEN
        CALL SUBERROR(
     &  'ERROR IN SPECIFICATION OF NO OF AGE CLASSES',
     &  IFATAL,IOERROR)
      END IF

      NOAGEPI = NOAGEP

      RETURN
      END !ReadAGEP


C**********************************************************************
      SUBROUTINE READAERO(UFILE, EXTWINDI)
C Read in aerodynamic info about canopy
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /AERO/ EXTWIND
      INTEGER UFILE

      EXTWIND = 0.0
      REWIND (UFILE)
      READ (UFILE,AERO,IOSTAT = IOERROR)
C      IF (IOERROR.NE.0) THEN
C        CALL SUBERROR(
C     &  'WARNING: USING DEFAULT VALUES FOR LEAF INCIDENCE ANGLE',
C     &  IWARN,IOERROR)
C      END IF

      EXTWINDI = EXTWIND

      RETURN
      END !ReadAero
      

C**********************************************************************
      SUBROUTINE READLIA(UFILE, NALPHAI, ALPHA, FALPHAI)
C Read in leaf incidence angles.
C Must return: NALPHA = number of angle classes,
C              ALPHA = angle of each angle class,
C              FALPHA = fraction of leaf area in each angle class.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /LIA/ ELP,NALPHA,FALPHA,AVGANG
      INTEGER UFILE
      REAL FALPHA(MAXANG),FALPHAI(MAXANG),ALPHA(MAXANG)

C Alternatives: (a) AVGANG > 0
C The mean leaf inclination angle is given by AVGANG. This is used to
C generate the LIA distribution assuming an elliptical distribution.
C               (b) ELP > 0.0
C    NALPHA = 1: there is just one leaf angle class, average angle = AVGLIA(ELP)
C    NALPHA > 1 (max 9): there are NALPHA leaf angle classes, the distribution
C of angles is elliptical with parameter ELP.
C               (c) ELP < 0.0
C The proportion of leaf area in each angle class is read in. Number of angle
C classes is given by NALPHA.
C The distribution of angles is read into array FALPHA(MAXANG).

C Default values
      AVGANG = -1.0
      NALPHA = 1
      ELP = 1.0
      FALPHA(1) = 1.0

C Read file
      REWIND (UFILE)
      READ (UFILE,LIA,IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &  'WARNING: USING DEFAULT VALUES FOR LEAF INCIDENCE ANGLE',
     &  IWARN,IOERROR)
      END IF

      IF (NALPHA.EQ.1) THEN
        IF (AVGANG.GT.0.0) THEN
          ALPHA(1) = AVGANG*PID180
        ELSE IF (ELP.GT.0.0) THEN
          ALPHA(1) = AVGLIA(ELP)
        END IF
        FALPHA(1) = 1.0

      ELSE
        IF (AVGANG.GT.0.0) ELP = CALCELP(AVGANG)
        IF (ELP.GT.0.00) THEN
          DALPHA = PID2/FLOAT(NALPHA)
          DO 60 IALP = 1,NALPHA
             ALPHA(IALP) = (IALP-0.5)*DALPHA
60        CONTINUE
          CALL ANGLE(ELP,NALPHA,FALPHA)
        END IF
      END IF

      NALPHAI = NALPHA
      DO 70 IANG = 1,MAXANG
        FALPHAI(IANG) = FALPHA(IANG)
70    CONTINUE

      RETURN
      END !ReadLIA


C**********************************************************************
      SUBROUTINE READZEN(UFILE,NUMPNT,NOLAYI,NZENI,NAZI,DIFZEN)
C Read in number of layers and angles to consider.
C Vertical layers:
C   NOLAY is the number of vertical layers used to calculate radiation
C interception (NUMPNT = NOLAY * points per layer).
C Angles for diffuse radiation:
C   NZEN = no. of zenith angles (default 5)
C   NAZ = no. of azimuth angles (default 11) - no maximum enforced
C   DIFZEN = array containing the NZEN zenith angles (radians).
C**********************************************************************

      INCLUDE 'maestcom'
      NAMELIST /DIFFANG/ NOLAY,NZEN,NAZ
      INTEGER UFILE
      REAL DIFZEN(MAXANG)

C Default values
      NOLAY = 6
      NAZ = 11
      NZEN = 5

C Read file
      REWIND (UFILE)
      READ (UFILE, DIFFANG, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &  'WARNING: USING DEFAULT VALUES FOR NOLAY,NAZ & NZEN',
     &  IWARN,IOERROR)
      END IF

C Check values do not exceed maxima.
      IF (NOLAY.GT.MAXLAY) THEN
        CALL SUBERROR(
     &  'WARNING: EXCEEDED MAXIMUM NO OF LAYERS',IWARN,0)
        NOLAY = MAXLAY
      END IF
      IF (NZEN.GT.MAXANG) THEN
        CALL SUBERROR(
     &  'WARNING: EXCEEDED MAXIMUM NO OF ZENITH ANGLES',IWARN,0)
        NZEN = MAXANG
      END IF

C Set zenith angle array
      DFZ = PID2/FLOAT(NZEN)
      DO 90 I = 1,NZEN
         DIFZEN(I) = (I-0.5)*DFZ
   90 CONTINUE

      NUMPNT = NOLAY * PPLAY
      NOLAYI = NOLAY
      NZENI = NZEN
      NAZI = NAZ

      RETURN
      END !ReadZen


C**********************************************************************
      SUBROUTINE READALLOM(UFILE, COEFFTI, EXPONTI, WINTERCI,
     &  BCOEFFTI, BEXPONTI, BINTERCI,
     &  RCOEFFTI, REXPONTI, RINTERCI, FRFRACI)
C Read in coefficients for allometric relation between stem
C height, diameter, and stem mass.
C**********************************************************************

      INCLUDE 'maestcom'
      NAMELIST /ALLOM/ COEFFT, EXPONT, WINTERC
      NAMELIST /ALLOMB/ BCOEFFT, BEXPONT, BINTERC
      NAMELIST /ALLOMR/ RCOEFFT, REXPONT, RINTERC, FRFRAC
      INTEGER UFILE

C Default values
      COEFFT = 0.0
      WINTERC = 0.0
      EXPONT = 2.
      RCOEFFT = 0.0
      RINTERC = 0.0
      REXPONT = 2.
      FRFRAC = 1.0
      BCOEFFT = 0.0
      BEXPONT = 2.

C Read file
      REWIND (UFILE)
      READ (UFILE, ALLOM, IOSTAT = IOERROR)
      IF (COEFFT.GT.0.0) THEN
        CALL SUBERROR('CALCULATING WOODY RESPIRATION',
     &  IWARN,IOERROR)
      END IF

      REWIND (UFILE)
      READ (UFILE, ALLOMB, IOSTAT = IOERROR)
      IF (RCOEFFT.GT.0.0) THEN
        CALL SUBERROR('CALCULATING BRANCH RESPIRATION',
     &  IWARN,IOERROR)
      END IF

      REWIND (UFILE)
      READ (UFILE, ALLOMR, IOSTAT = IOERROR)
      IF (RCOEFFT.GT.0.0) THEN
        CALL SUBERROR('CALCULATING ROOT RESPIRATION',
     &  IWARN,IOERROR)
      END IF

      COEFFTI = COEFFT
      EXPONTI = EXPONT
      WINTERCI = WINTERC
      BCOEFFTI = BCOEFFT
      BEXPONTI = BEXPONT
      BINTERCI = BINTERC
      RCOEFFTI = RCOEFFT
      REXPONTI = REXPONT
      RINTERCI = RINTERC
      FRFRACI = FRFRAC

      RETURN
      END !ReadAllom


C**********************************************************************
      SUBROUTINE READPPT(UFILE,PFLAI,CANCAPLAI)
C Read in parameters for canopy rainfall interception model.
C**********************************************************************

      INCLUDE 'maestcom'
      NAMELIST /INTERCEP/ PFLA,CANCAPLA
      INTEGER UFILE

C Default values
      PFLA = 0.0
      CANCAPLA = 0.0

C Read file
      REWIND (UFILE)
      READ (UFILE, INTERCEP, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING:CANOPY INTERCEPTION NOT CALCULATED',
     &  IWARN,IOERROR)
      END IF

      PFLAI = PFLA
C      CANCAPLAI = CANCAPLA*1E6/18.
      CANCAPLAI = CANCAPLA

      RETURN
      END !ReadPPT


C**********************************************************************
      SUBROUTINE READPROP(UFILE, NOAGEC, NOAGEP, PROPI)
C Read in the proportion of leaf area in each age class. The number
C of age classes is given by NOAGEC (for beta distributions) and
C NOAGEP (for physiological parameters).
C**********************************************************************

      INCLUDE 'MAESTCOM'

      NAMELIST /PHENOL/ PROP
      REAL PROP(MAXC),PROPI(MAXC)
      INTEGER UFILE

      DO 10 I = 1,MAXC
        PROP(I) = 0.0
10    CONTINUE

C Find number of age classes
      NAGE = NOAGEC
      IF (NOAGEP.GT.NAGE) NAGE = NOAGEP

C Only need to read in proportions if there is > 1 age class.
      IF (NAGE.GT.1) THEN
        REWIND (UFILE)
        READ (UFILE, PHENOL, IOSTAT = IOERROR)
        IF (IOERROR.NE.0) THEN
          CALL SUBERROR('ERROR: NEED DATA ON AGE CLASS PROPORTIONS',
     &    IFATAL,IOERROR)
        END IF
      ELSE
        PROP(1) = 1.0
      END IF

C Check proportions sum to one
      TOTAL = 0.0
      DO 20 I = 1,NAGE
        TOTAL = TOTAL + PROP(I)
20    CONTINUE
      IF ((TOTAL.LT.0.9999).OR.(TOTAL.GT.1.0001))
     &  CALL SUBERROR('ERROR: AGE CLASS PROPORTIONS DO NOT SUM TO ONE',
     &  IFATAL,0)

      DO 30 IAGE = 1,MAXC
        PROPI(IAGE) = PROP(IAGE)
30    CONTINUE

      RETURN
      END !ReadProp


C**********************************************************************
      SUBROUTINE READBETA(UFILE, NOAGECI, JLEAFI, BPTI, RANDOMI)
C Read in beta distributions for leaf area. The number of age classes
C for which beta distributions are specified is given by NOAGEC.
C Function returns:
C switch JLEAF:
C   JLEAF = 0: Uniform leaf area density assumed
C   JLEAF = 1: Leaf area density is variable in vertical direction
C   JLEAF = 2: Leaf area density is variable in both horizontal & vertical
C array BPT: gives the coefficients of the beta distributions
C the clumping factor RANDOM (= shoot projected area: leaf projected area).
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /LADD/ NOAGEC,JLEAF,BPT,BPTEXT,RANDOM
      REAL BPT(6,MAXC),BPTEXT(2,MAXC),BPTI(8,MAXC)
      INTEGER UFILE
C BM 11/99 added another parameter to BETA function - input via BPTEXT - default value 1

C Default values
      NOAGEC = 1
      JLEAF = 0
      RANDOM = 1.0
      DO 10 J = 1,MAXC
	  BPTEXT(1,J) = 1.0
	  BPTEXT(2,J) = 1.0
        DO 10 I = 1,6
          BPT(I,J) = 0.0
10    CONTINUE

C Read file
      REWIND (UFILE)
      READ (UFILE, LADD, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &  'WARNING: DEFAULT VALUE: UNIFORM LEAF AREA DENSITY',
     &  IWARN,IOERROR)
      END IF
      NOAGECI = NOAGEC
      JLEAFI = JLEAF
      RANDOMI = RANDOM
      DO 20 J = 1,MAXC
	  DO 20 I = 0,1
	  BPTI(I*4+1,J) = BPT(I*3+1,J)
	  BPTI(I*4+2,J) = BPT(I*3+2,J)
	  BPTI(I*4+3,J) = BPT(I*3+3,J)
	  BPTI(I*4+4,J) = BPTEXT(I+1,J)
20    CONTINUE

C Elementary error checking:
C If not using uniform leaf area density
      IF (JLEAF.NE.0) THEN
C If the first coefficient of the first age class = 0, or
        IF ((BPT(1,1).EQ.0.0).OR.
C If there's >1 age class and the first coefft of the second age class = 0
     &      ((NOAGEC.GT.1).AND.(BPT(1,2).EQ.0.0))) THEN
          CALL SUBERROR(
     &    'ERROR: MISSING BETA FUNCTION COEFFICIENTS',
     &    IFATAL,IOERROR)
        END IF
      END IF

      RETURN
      END !ReadBeta


C**********************************************************************
      SUBROUTINE READMODEL(UFILE, GSMODI, JMMODI, RDMODI, 
     &  SSMODI, RWMODI, ITERMAXI)
C Read in flags which control the physiological models used.
C MODELGS - controls model of stomatal conductance
C MODELJM - whether JMAX,VCMAX read in (0) or calculated from leaf N (1)
C MODELRD - whether RD is a function of leaf area (0) or leaf N (1)
C MODELSS - whether photosynthesis is calculated for sunlit & shaded leaves
C   together or separately
C ITERMAX - no. of iterations to be used in solving for leaf temperature
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /MODEL/ MODELGS,MODELJM,MODELRD,MODELRW,MODELSS,ITERMAX
      INTEGER UFILE,GSMODI,JMMODI,RDMODI,SSMODI,RWMODI,ITERMAXI

C Default values
      MODELGS = 0     ! The Jarvis model of stom cond
      MODELJM = 0     ! Jmax & Vcmax read in directly
      MODELRD = 0     ! Rd0 read in directly
      MODELRW = 0     ! RW values read in directly
      MODELSS = 0     ! sunlit & shade calculations separate
      ITERMAX = 0     ! The leaf temperature is not calculated

C Read file
      REWIND (UFILE)
      READ (UFILE, MODEL, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &  'WARNING: DEFAULT VALUES USED FOR PHYSIOLOGICAL CONTROL',
     &  IWARN,IOERROR)
      END IF
      GSMODI = MODELGS
      RDMODI = MODELRD
      RWMODI = MODELRW
      JMMODI = MODELJM
      SSMODI = MODELSS
      ITERMAXI = ITERMAX

      RETURN
      END !ReadModel


C**********************************************************************
      SUBROUTINE READHIST(UFILE, IOHIST, BINSIZEI)
C Read in information for the PAR histogram if required.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /HISTO/ BINSIZE
      INTEGER UFILE

      IF (IOHIST.EQ.1) THEN

        REWIND (UFILE)
        READ (UFILE, HISTO, IOSTAT = IOERROR)
        IF (IOERROR.NE.0) THEN
          CALL SUBERROR(
     &    'ERROR: MISSING INFO FOR HISTOGRAM',
     &    IFATAL,IOERROR)
        END IF
        BINSIZEI = BINSIZE

      END IF
      RETURN
      END ! ReadHist


C**********************************************************************
      SUBROUTINE READCCSCEN(UFILE, ICC, CO2INCI, TINCI)
C Read in details of climate change scenario.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /CCSCEN/ CO2INC, TINC
      INTEGER UFILE

C Default values
      ICC = 0
      CO2INC = 0.0
      TINC = 0.0
C Read file
      REWIND (UFILE)
      READ (UFILE, CCSCEN, IOSTAT = ICC)
      IF (ICC.EQ.0) 
     &  CALL SUBERROR('APPLYING CLIMATE CHANGE SCENARIO',IWARN,0)
      CO2INCI = CO2INC
      TINCI = TINC

      RETURN
      END ! ReadCCScen


C**********************************************************************
      SUBROUTINE READOTC(UFILE, IOTC, TOTCI, WINDOTCI,
     &  PAROTCI, FBEAMOTCI)
C Subroutine to read in parameters describing effect of OTC on met data.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /OTC/ TOTC,WINDOTC,PAROTC,FBEAMOTC
      INTEGER UFILE

C Default values
      IOTC = 0
      TOTC = 0
      WINDOTC = -1
      PAROTC = 1.0
      FBEAMOTC = 1.0
C Read file
      REWIND (UFILE)
      READ (UFILE, OTC, IOSTAT = IOTC)
      IF (IOTC.EQ.0) 
     &  CALL SUBERROR('APPLYING EFFECTS OF OTC ON MET DATA',IWARN,0)
      TOTCI = TOTC
      WINDOTCI = WINDOTC
      PAROTCI = PAROTC
      FBEAMOTCI = FBEAMOTC

      RETURN
      END ! ReadOTC

C**********************************************************************
      SUBROUTINE READABSRP(UFILE,NOLAY,ABSRP,REFLEC,TRANS,RHOSOLI)
C Read leaf absorptances and reflectances for 3 wavelengths.
C Required input: File unit (UFILE), No of layers (NOLAY).
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /ABSORP/ NOLAYERS,RHOSOL,ATAU,ARHO
      REAL ABSRP(MAXLAY,3),REFLEC(MAXLAY,3),TRANS(MAXLAY,3)
      REAL ARHO(MAXLAY*3),ATAU(MAXLAY*3)
      REAL RHOSOL(3),RHOSOLI(3),X(MAXLAY+1)
      INTEGER UFILE

C Read file
      REWIND (UFILE)
      READ (UFILE,ABSORP,IOSTAT = IOERROR)
      IF (IOERROR.NE.0)
     &  CALL SUBERROR(
     &  'ERROR READING REFLECTANCE / TRANSMITTANCE',
     &  IFATAL,IOERROR)

      DO 10 IP = 1,NOLAYERS+1
        X(IP) = (REAL(IP) - 1.0)/REAL(NOLAYERS)
10    CONTINUE

C Check values & set absorptance array.
      DO 30 J = 1,3
        RHOSOLI(J) = RHOSOL(J)
        DO 20 I = 1,NOLAY
          D1 = (REAL(I) - 1)/REAL(NOLAY)
          D2 = (REAL(I))/REAL(NOLAY)
          REFLEC(I,J) = RINTEG(D1,D2,X,ARHO,(J-1)*NOLAYERS,NOLAY)
          TRANS(I,J) = RINTEG(D1,D2,X,ATAU,(J-1)*NOLAYERS,NOLAY)
          ABSRP(I,J) = 1.0 - REFLEC(I,J) - TRANS(I,J)
20    CONTINUE

30    CONTINUE

100   FORMAT(10(1X,F8.3))

      RETURN
      END !ReadAbsrp


C**********************************************************************
      SUBROUTINE READGS(UFILE,MODELGS,
     &  GSREFI,GSMINI,PAR0I,D0I,VK1I,VK2I,VPD1I,VPD2I,VMFD0I,
     &  GSJAI,GSJBI,T0I,TREFI,TMAXI,SMD1I,SMD2I,
     &  G0I,D0LI,GAMMAI,G1I,WLEAFI,NSIDESI)
C Get stomatal conductance parameters.
C Required input: File unit (UFILE), Stom cond model (MODELGS).
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /JARGS/ GSREF,GSMIN,PAR0,D0,VPD1,VPD2,VMFD0,GSJA,GSJB,
     &  T0,TREF,TMAX,SMD1,SMD2,VK1,VK2,WLEAF,NSIDES
      NAMELIST /BBGS/ G0,G1,GAMMA,WLEAF,NSIDES,SMD1,SMD2
      NAMELIST /BBLGS/ G0,G1,D0L,GAMMA,WLEAF,NSIDES,SMD1,SMD2
      INTEGER UFILE

C Defaults
      GSMIN = 0.001
	SMD1 = 0.0

      REWIND(UFILE)

      IF (MODELGS.EQ.2) THEN    ! Ball-Berry model parameters
        READ (UFILE,BBGS,IOSTAT = IOERROR)
        G0I = G0
        G1I = G1
        GAMMAI = GAMMA
	  SMD1I = SMD1
	  SMD2I = SMD2

      ELSE IF (MODELGS.EQ.3) THEN   ! Ball-Berry Leuning model parameters
        READ (UFILE, BBLGS,IOSTAT = IOERROR)
        G0I = G0
        G1I = G1
        GAMMAI = GAMMA
        D0LI = D0L
	  SMD1I = SMD1
	  SMD2I = SMD2
        IF (D0L.LT.0.0) 
     &    CALL SUBERROR('ERROR IN GS PARAMETERS: D0L MUST BE > 0',
     &    IFATAL,0)

      ELSE    ! Jarvis model parameters
	  I0 = 0.0
	  TMAX = 0.0
	  GSJA = 0.0
	  GSJB = 0.0
	  SMD1 = 0.0
        READ (UFILE,JARGS,IOSTAT = IOERROR)
        GSREFI = GSREF
        GSMINI = GSMIN
        PAR0I = PAR0
        D0I = D0
        VPD1I = VPD1
        VPD2I = VPD2
	  VK1I = VK1
	  VK2I = VK2
        VMFD0I = VMFD0
        GSJAI = GSJA
        GSJBI = GSJB
        T0I = T0
        TREFI = TREF
        TMAXI = TMAX
	  SMD1I = SMD1
	  SMD2I = SMD2

      END IF

      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &  'ERROR READING STOMATAL CONDUCTANCE PARAMETERS',
     &  IFATAL,IOERROR)
      END IF
      WLEAFI = WLEAF
      NSIDESI = NSIDES

      RETURN
      END !ReadGS


C**********************************************************************
      SUBROUTINE READLEAFN(UFILE,MODELJM,MODELRD,NOLAY,NOAGEP,
     &  NONDATES,DATESN,VALUESN)
C If MODELJM or MODELRD = 1, then leaf N contents are required to
C calculate Jmax or Rd. This subroutine then reads in the leaf N.
C Leaf N may be specified for a maximum of MAXPDATE dates -
C linear interpolation is used between those dates.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      REAL VALUESN(MAXPDATE,MAXLAY,MAXC)
      INTEGER DATESN(MAXPDATE)
      INTEGER UFILE

      IF ((MODELJM.EQ.1).OR.(MODELRD.EQ.1)) THEN  !Leaf n contents reqd.
        CALL READPHYARRAY(UFILE,1,NOLAY,NOAGEP,
     &  NONDATES,DATESN,VALUESN)
      END IF

      RETURN
      END !ReadLeafN


C**********************************************************************
      SUBROUTINE READPHYARRAY(UFILE,NARRAY,NOLAY,NOAGEP,
     &  NDATE,IDATES,VALUESI)
C Read in an array of physiological parameters from UFILE.
C NARRAY is the number of the array to be read (1 = LEAFN; 2 = JMAX;
C 3 = VCMAX; 4 = RD; 5 = SLA; 6 = ALPHAJ)
C Parameters are given for up to MAXPDATE dates, up to MAXLAY layers,
C and up to MAXC age classes.
C**********************************************************************

      INCLUDE 'MAESTCOM'

	REAL VALUES(MAXPDATE*MAXLAY*MAXC),VALUESI(MAXPDATE,MAXLAY,MAXC)
      REAL X(MAXLAY+1)
      CHARACTER*10 DATES(MAXPDATE)
      INTEGER IDATES(MAXPDATE)
      INTEGER UFILE

      NAMELIST /NFOLCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /NFOL/ DATES,VALUES
      NAMELIST /JMAXCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /JMAX/ DATES,VALUES
      NAMELIST /VCMAXCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /VCMAX/ DATES,VALUES
      NAMELIST /RDCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /RD/ DATES,VALUES
      NAMELIST /SLACON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /SLA/ DATES,VALUES
      NAMELIST /AJQCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /AJQ/ DATES,VALUES

C Default values
      NODATES = 1
      NOLAYERS = 1
      NOAGES = 1
      DATES(1) = '01/01/50'
      DO 10 I = 1,MAXPDATE*MAXLAY*MAXC
        VALUES(I) = -1.0
10    CONTINUE

C Read array size from file
      REWIND(UFILE)
      IF (NARRAY.EQ.1) THEN
        READ(UFILE,NFOLCON,IOSTAT=IOERROR)
        CALL SUBERROR('LEAF N ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.2) THEN
        READ(UFILE,JMAXCON,IOSTAT=IOERROR)
        CALL SUBERROR('JMAX ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.3) THEN
        READ(UFILE,VCMAXCON,IOSTAT=IOERROR)
        CALL SUBERROR('VCMAX ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.4) THEN
        READ(UFILE,RDCON,IOSTAT=IOERROR)
        CALL SUBERROR('RD ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.5) THEN
        READ(UFILE,SLACON,IOSTAT=IOERROR)
        IF (IOERROR.EQ.0) THEN
          CALL SUBERROR('SLA ARRAY:',IWARN,0)
        ELSE
          CALL SUBERROR(
     &      'NO SLA VALUES: FOLIAGE GROWTH RESP NOT CALCULATED',
     &      IWARN,IOERROR)
          NODATES = 0
          RETURN ! No SLA values - end subroutine here
        END IF
      ELSE IF (NARRAY.EQ.6) THEN
        READ(UFILE,AJQCON,IOSTAT=IOERROR)
        IF (IOERROR.EQ.0) THEN
          CALL SUBERROR('AJQ ARRAY:',IWARN,0)
        ELSE
          CALL SUBERROR(
     &      'NO AJQ ARRAY; DEFAULT VALUE USED',
     &      IWARN,IOERROR)
          NODATES = 0
          RETURN ! No AJQ values - end subroutine here
        END IF
      END IF

C Error handling
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &  'WARNING: PROBLEM READING ARRAY SIZE',IWARN,IOERROR)
      END IF
      IF ((NODATES.GT.MAXPDATE).OR.(NOLAYERS.GT.MAXLAY).OR.
     &    (NOAGES.GT.MAXC)) THEN
        CALL SUBERROR(
     &  'WARNING: ARRAY EXCEEDED BOUNDS: SOME DATA LOST',
     &  IWARN,IOERROR)
        IF (NODATES.GT.MAXPDATE) NODATES = MAXPDATE
        IF (NOLAYERS.GT.MAXLAY) NOLAYERS = MAXLAY
        IF (NOAGES.GT.MAXC) NOAGES = MAXC
      END IF
      IF ((NOAGES.GT.1).AND.(NOAGEP.NE.NOAGES)) THEN
        CALL SUBERROR(
     &  'ERROR: PHYSIOLOGY CLASSES & PARAMETER CLASSES DO NOT COINCIDE',
     &  IWARN,IOERROR)
      END IF

C Read arrays from file
      REWIND(UFILE)
      IF (NARRAY.EQ.1) THEN
        READ(UFILE,NFOL,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.2) THEN
        READ(UFILE,JMAX,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.3) THEN
        READ(UFILE,VCMAX,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.4) THEN
        READ(UFILE,RD,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.5) THEN
        READ(UFILE,SLA,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.6) THEN
        READ(UFILE,AJQ,IOSTAT=IOERROR)
      END IF

C Error handling
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &  'ERROR: PROBLEM READING ARRAY',IFATAL,IOERROR)
      END IF

C Check values
      DO 20 I = 1,NOAGES*NOLAYERS*NODATES
        IF (VALUES(I).LT.0.0)
     &    CALL SUBERROR('ERROR: VALUES MISSING FROM ARRAY',
     &    IFATAL,0)
20    CONTINUE

C Set up x array
      DO 30 IP = 1,NOLAYERS+1
        X(IP) = (REAL(IP) - 1.0)/REAL(NOLAYERS)
30    CONTINUE

C Assign arrays to data
      DO 40 IDATE = 1,NODATES
        IDATES(IDATE) = IDATE50(DATES(IDATE))
40    CONTINUE
      DO 50 IDATE = 1,NODATES
        DO 50 IAGE = 1,NOAGES
          DO 50 I = 1,NOLAY
            D1 = (REAL(I) - 1)/REAL(NOLAY)
            D2 = (REAL(I))/REAL(NOLAY)
            VALUESI(IDATE,I,IAGE) =
     &  RINTEG(D1,D2,X,VALUES,
     &  (IDATE-1)*NOAGES*NOLAYERS+(IAGE-1)*NOLAYERS,NOLAY)
50    CONTINUE

C Fill in array in case no of ages is less than NOAGEP
	  IF (NOAGES.LT.NOAGEP) THEN 
        DO 60 IDATE = 1,NODATES
          DO 60 ILAY = 1,NOLAY
	      DO 60 I = NOAGES+1,NOAGEP
	        VALUESI(IDATE,ILAY,I) = VALUESI(IDATE,ILAY,1)
60	    CONTINUE 
	  END IF

      NDATE = NODATES

      RETURN
      END !ReadPhyArray


C**********************************************************************
      SUBROUTINE READJMAX(UFILE, MODELJM, NOLAY, NOAGEP, 
     &  NONDATES, DATESN, LEAFN, 
     &  NOJDATES, DATESJ, JMAXTABLE,
     &  NOVDATES, DATESV, VCMAXTABLE, 
     &  NOADATES, DATESA, AJQTABLE,
     &  IECOI, EAVJI, EDVJI, DELSJI, EAVCI, EDVCI, DELSCI, 
     &  TVJUPI, TVJDNI, THETAI)
C Read in all parameters to do with Jmax and Vcmax.
C If MODELJM = 1, use the leaf N contents to calculate Jmax and Vcmax;
C otherwise read in arrays directly.
C**********************************************************************

      INCLUDE 'MAESTCOM'

	REAL LEAFN(MAXPDATE,MAXLAY,MAXC)
      REAL JMAXTABLE(MAXPDATE,MAXLAY,MAXC)
      REAL VCMAXTABLE(MAXPDATE,MAXLAY,MAXC)
	REAL AJQTABLE(MAXPDATE,MAXLAY,MAXC)
      REAL JMAXA,JMAXB
      INTEGER DATESN(MAXPDATE), DATESJ(MAXPDATE)
	INTEGER DATESV(MAXPDATE), DATESA(MAXPDATE)
      INTEGER UFILE

      NAMELIST /JMAXPARS/ THETA,EAVJ,EDVJ,DELSJ,AJQ,IECO
      NAMELIST /VCMAXPARS/ EAVC,EDVC,DELSC,TVJUP,TVJDN
      NAMELIST /JMAXN/ JMAXA,JMAXB,VCMAXA,VCMAXB

      IF (MODELJM.EQ.1) THEN   ! Calculate Jmax, Vcmax from leaf N
        REWIND(UFILE)
        READ (UFILE,JMAXN,IOSTAT = IOERROR)
        IF (IOERROR.NE.0)
     &    CALL SUBERROR('ERROR: MISSING CONSTANTS FOR JMAX/N',
     &    IFATAL,IOERROR)
        NOJDATES = NONDATES
        NOVDATES = NONDATES
        DO 10 IDATE = 1,NONDATES
          DATESJ(IDATE) = DATESN(IDATE)
          DATESV(IDATE) = DATESN(IDATE)
          DO 10 ILAY = 1,NOLAY
            DO 10 IAGE = 1,NOAGEP
              JMAXTABLE(IDATE,ILAY,IAGE) =
     &          JMAXA*LEAFN(IDATE,ILAY,IAGE) + JMAXB
              VCMAXTABLE(IDATE,ILAY,IAGE) =
     &          VCMAXA*LEAFN(IDATE,ILAY,IAGE) + VCMAXB
10      CONTINUE

      ELSE                     ! Read in Jmax, Vcmax arrays
        CALL READPHYARRAY(UFILE,2,NOLAY,NOAGEP,
     &    NOJDATES,DATESJ,JMAXTABLE)
        CALL READPHYARRAY(UFILE,3,NOLAY,NOAGEP,
     &    NOVDATES,DATESV,VCMAXTABLE)
      END IF

C Read in additional parameters
      REWIND (UFILE)
	  AJQ = ALPHAQ !Default values
	  EDVC = 0.0
	  DELSC = 0.0
	  IECO = 1 ! Ecocraft formulation of T-deps of Km and Gamma. For Montpied formulation, put 0. 
	  TVJUP = -100.0
	  TVJDN = -100.0

      READ (UFILE,JMAXPARS,IOSTAT = IOERROR)
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('INPUT ERROR: MISSING JMAXPARS',IFATAL,IOERROR)
      REWIND (UFILE)
      READ (UFILE,VCMAXPARS,IOSTAT = IOERROR)
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('INPUT ERROR: MISSING VCMAXPARS',IFATAL,IOERROR)

      EAVCI = EAVC
      EDVCI = EDVC
      DELSCI = DELSC
      EAVJI = EAVJ
      EDVJI = EDVJ
      DELSJI = DELSJ
      THETAI = THETA
	IECOI = IECO
	TVJUPI = TVJUP
	TVJDNI = TVJDN

C Read in quantum yield array
      CALL READPHYARRAY(UFILE,6,NOLAY,NOAGEP,
     &    NOADATES,DATESA,AJQTABLE)
      IF (NOADATES.EQ.0) THEN
	  NOADATES = 1
        DO 20 ILAY = 1,NOLAY
          DO 20 IAGE = 1,NOAGEP
	      AJQTABLE(1,ILAY,IAGE) = AJQ
20	  CONTINUE
	END IF

      RETURN
      END !ReadJmax


C**********************************************************************
      SUBROUTINE READRD(UFILE, MODELRD, NOLAY, NOAGEP, NONDATES,
     &  DATESN, LEAFN, NORDATES, DATESRD, RDTABLE,
     &  Q10FI, RTEMPI, DAYRESPI, EFFYRFI, TBELOWI)
C Read in all parameters to do with leaf respiration rate.
C If MODELRD = 1, use the leaf N contents to calculate Rd0;
C otherwise read in array directly.
C**********************************************************************

      INCLUDE 'MAESTCOM'


      REAL LEAFN(MAXPDATE,MAXLAY,MAXC)
      REAL RDTABLE(MAXPDATE,MAXLAY,MAXC)
      INTEGER DATESN(MAXPDATE), DATESRD(MAXPDATE)
      INTEGER UFILE

	  NAMELIST /RDPARS/ Q10F, RTEMP, TBELOW, DAYRESP, EFFYRF
      NAMELIST /RDN/ RDA, RDB

      IF (MODELRD.EQ.1) THEN   ! Calculate Jmax, Vcmax from leaf N
        REWIND(UFILE)
        READ (UFILE,RDN,IOSTAT = IOERROR)
        IF (IOERROR.NE.0)
     &    CALL SUBERROR('ERROR: MISSING CONSTANTS FOR RD/N',
     &    IFATAL,IOERROR)
        NORDATES = NONDATES
        DO 10 IDATE = 1,NONDATES
          DATESRD(IDATE) = DATESN(IDATE)
          DO 10 ILAY = 1,NOLAY
            DO 10 IAGE = 1,NOAGEP
              RDTABLE(IDATE,ILAY,IAGE) =
     &          RDA*LEAFN(IDATE,ILAY,IAGE) + RDB
10      CONTINUE

      ELSE                     ! Read in RD0 arrays
        CALL READPHYARRAY(UFILE,4,NOLAY,NOAGEP,
     &    NORDATES,DATESRD,RDTABLE)
      END IF

C Read in additional parameters
      REWIND (UFILE)
      DAYRESP = 1.0            ! Default value - no effect of light
      RTEMP = 0.0              ! Default value - Rd at zero degrees
      EFFYRF = 0.0             ! Default value - don't calc growth respn
      Q10F = 0.0               ! Indicates missing value
      TBELOW = -100.0          ! Default value - Rd continues at all Ts
      READ (UFILE,RDPARS,IOSTAT = IOERROR)
      IF (Q10F.EQ.0.0)
     &  CALL SUBERROR('INPUT ERROR: MISSING Q10F',IFATAL,IOERROR)
      Q10FI = Q10F
      RTEMPI = RTEMP
      DAYRESPI = DAYRESP
      EFFYRFI = EFFYRF
      TBELOWI = TBELOW

      RETURN
      END !ReadRD


C**********************************************************************
      SUBROUTINE READRW(UFILE,MODELRW,EFFYRWI,RMWI,RTEMPWI,Q10WI,
     &  COLLA,COLLK,STEMSDWI,RMAI,STEMFORMI)
C Read in or calculate parameters to do with woody respiration rate.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /WRESP/ EFFY,RM,RMA,RTEMP,Q10W,STEMFORM
      NAMELIST /COLLWRESP/ CA,CK,STEMSDW
      INTEGER UFILE

      REWIND (UFILE)
      EFFY = 0.0        ! Missing value - growth respn not calculated
      RM = 0.0          ! Missing value - maint respn not calculated
      RMA = 0.0         ! Missing value - maint respn not calculated
      RTEMP = 0.0       ! Default temp at which maint resp specified
      Q10W = 0.0        ! Missing value - will cause error if RM spec
      READ (UFILE,WRESP,IOSTAT = IOERROR)
      IF (EFFY.EQ.0.0)
     &  CALL SUBERROR('WARNING: WOODY GROWTH RESP NOT CALCULATED',
     &  IWARN,IOERROR)
      IF ((MODELRW.NE.1).AND.RM.EQ.0.0.AND.RMA.EQ.0.0) THEN
        CALL SUBERROR('WARNING: WOODY MAINT RESP NOT CALCULATED',
     &  IWARN,IOERROR)
      ELSE IF (Q10W.EQ.0.0) THEN
        CALL SUBERROR('INPUT ERROR: MISSING Q10W',IFATAL,IOERROR)
      END IF

      Q10WI = Q10W
      RTEMPWI = RTEMP
      EFFYRWI = EFFY
      RMWI = RM
	  RMAI = RMA
	  STEMFORMI = STEMFORM

      IF (MODELRW.EQ.1) THEN ! Collelongo model
        REWIND(UFILE)        
        READ(UFILE,COLLWRESP,IOSTAT = IOERROR)
        IF (IOERROR.NE.0)
     &    CALL SUBERROR('ERROR: MISSING INFO FOR COLLELONGO RW MODEL',
     &    IFATAL,IOERROR)
        COLLA = CA
        COLLK = CK
        STEMSDWI = STEMSDW
      END IF

	  IF (RM.EQ.0.0.AND.RMA.GT.0.0) MODELRW = 2

      RETURN
      END !ReadRW


C**********************************************************************
      SUBROUTINE READRR(UFILE,RMCRI,RMFRI,Q10RI,RTEMPRI)
C Read in or calculate parameters to do with root respiration rate.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /RRESP/ RMCR,RMFR,Q10R,RTEMPR
      INTEGER UFILE

      REWIND (UFILE)
      RMCR = 0.0        ! Missing value - maint respn not calculated
      RMFR = 0.0        ! Missing value - maint respn not calculated
      RTEMPR = 0.0      ! Default temp at which maint resp specified
      Q10R = 0.0        ! Missing value - will cause error if RM spec
      READ (UFILE,RRESP,IOSTAT = IOERROR)
      IF ((RMCR.EQ.0.0).AND.(RMFR.EQ.0.0)) THEN
        CALL SUBERROR('WARNING: ROOT MAINT RESP NOT CALCULATED',
     &  IWARN,IOERROR)
      ELSE IF (Q10R.EQ.0.0) THEN
        CALL SUBERROR('INPUT ERROR: MISSING Q10R',IFATAL,IOERROR)
      END IF

      Q10RI = Q10R
      RTEMPRI = RTEMPR
      RMCRI = RMCR
      RMFRI = RMFR

      END !ReadRR


C**********************************************************************
      SUBROUTINE READRB(UFILE,RMBI,Q10BI,RTEMPBI)
C Read in or calculate parameters to do with branch respiration rate.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /BRESP/ RMB,Q10B,RTEMPB
      INTEGER UFILE

      REWIND (UFILE)
      RMB = 0.0        ! Missing value - maint respn not calculated
      RTEMPB = 0.0     ! Default temp at which maint resp specified
      Q10B = 0.0       ! Missing value - will cause error if RM spec
      READ (UFILE,BRESP,IOSTAT = IOERROR)
      IF (RMB.EQ.0.0) THEN
        CALL SUBERROR('WARNING: BRANCH MAINT RESP NOT CALCULATED',
     &  IWARN,IOERROR)
      ELSE IF (Q10B.EQ.0.0) THEN
        CALL SUBERROR('INPUT ERROR: MISSING Q10B',IFATAL,IOERROR)
      END IF

      Q10BI = Q10B
      RTEMPBI = RTEMPB
      RMBI = RMB

      END !ReadRB


C**********************************************************************
C CURRENTLY UNUSED 3/6/98 BEM
      SUBROUTINE READRSOIL(UFILE,RSI,RTEMPSI,Q10SI)
C Read in or calculate parameters to do with soil respiration rate.
C Initially assuming Q10 relationship - to be replaced by Hanson et al 1993.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /SRESP/ RSOIL, RTEMPS, Q10S
      INTEGER UFILE

      REWIND (UFILE)
      RSOIL = 0.0       ! Missing value - soil respn not calculated
      RTEMPS = 0.0      ! Default temp at which soil resp specified
      Q10S = 0.0        ! Missing value - will cause error if RS spec
      READ (UFILE,SRESP,IOSTAT = IOERROR)
      IF (RSOIL.EQ.0.0) THEN
        CALL SUBERROR('WARNING: SOIL RESPIRATION NOT CALCULATED',
     &  IWARN,IOERROR)
      ELSE IF (Q10S.EQ.0.0) THEN
        CALL SUBERROR('INPUT ERROR: MISSING Q10S',IFATAL,IOERROR)
      END IF

      Q10SI = Q10S
      RTEMPSI = RTEMPS
      RSI = RSOIL

      RETURN
      END !ReadRSoil


C**********************************************************************
      SUBROUTINE PHYINTERP(IDATE,NODATES,IDATEARR,PARAMTABLE,
     &  NOLAY,NOAGEP,PARAMS)
C Interpolate physiological parameters for this date from the
C date arrays read in from physiology file.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      INTEGER IDATEARR(MAXPDATE)
      REAL PARAMTABLE(MAXPDATE,MAXLAY,MAXC)
      REAL PARAMS(MAXLAY,MAXC)


C If only one date, or before the starting date, take first value
      IF ((NODATES.EQ.1).OR.IDATE.LE.IDATEARR(1)) THEN
        DO 10 IAGE = 1,NOAGEP
          DO 10 ILAY = 1,NOLAY
            PARAMS(ILAY,IAGE) = PARAMTABLE(1,ILAY,IAGE)
10      CONTINUE

C If after the final date, take last value
      ELSE IF (IDATE.GT.IDATEARR(NODATES)) THEN
        DO 20 IAGE = 1,NOAGEP
          DO 20 ILAY = 1,NOLAY
            PARAMS(ILAY,IAGE) = PARAMTABLE(NODATES,ILAY,IAGE)
20      CONTINUE

C Otherwise have to interpolate
      ELSE
        INDEX = 1
        DO WHILE (IDATE.GT.IDATEARR(INDEX))
          INDEX = INDEX + 1
        END DO
        SLOPE = REAL((IDATE - IDATEARR(INDEX-1)))/
     &      REAL((IDATEARR(INDEX) - IDATEARR(INDEX-1)))
        DO 30 IAGE = 1,NOAGEP
          DO 30 ILAY = 1,NOLAY
            Y1 = PARAMTABLE(INDEX-1, ILAY, IAGE)
            Y2 = PARAMTABLE(INDEX, ILAY, IAGE)
            PARAMS(ILAY,IAGE) = Y1 + SLOPE*(Y2 - Y1)
30      CONTINUE

      END IF

      RETURN
      END !PhyInterp


C**********************************************************************
      FUNCTION RINTEG(D1,D2,X,VALUES,IOFFSET,NOLAY)
C Estimate physiological parameters for each canopy layer by averaging.
C At the moment, this is done by integrating the physiological parameter
C over the height of the layer, assuming foliage is evenly distributed.
C Strictly speaking, should probably take the beta-distribution of
C foliage into account - but since the specification of physiology for
C different layers is woolly anyway, this is a good approx.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL X(MAXLAY+1),VALUES(MAXLAY*MAXC*MAXPDATE)
C NB The array VALUES can be of variable size - since this routine is used
C to handle reflectance and transmittance, and leaf N, Jmax, Vcmax and Rd
C arrays (as of 8/97). Here it should be set to the maximum size of
C any array to be passed.

      RINTEG = 0.0
      IP = 1

      DO WHILE (X(IP+1).LT.D1)
        IP = IP+1
      END DO

      DO WHILE (D2.GT.X(IP+1))
        RINTEG = RINTEG + NOLAY*
     &      VALUES(IP+IOFFSET)*(X(IP+1)-AMAX1(D1,X(IP)))
        IP = IP+1
      END DO

      RINTEG = RINTEG + NOLAY*
     &    VALUES(IP+IOFFSET)*(D2-AMAX1(D1,X(IP)))

100   FORMAT(4(1X,F8.3))

      RETURN
      END !RInteg


C**********************************************************************
      SUBROUTINE READCONTREES(UFILE,NOALLTREES,DX,DY,XMAX,YMAX,
     &  NOTREESI,NOTARGETSI,ITARGETSI)
C Read in controls about tree numbers: need 
C (1) no of trees to be included in the shading calculations
C (2) the numbers of the target trees. 
C (1) NOTREES is the number of trees to be considered in the calculation.
C (Default: NOTREES = all trees in plot).
C (2) ITARGETS is an array with the numbers of the target trees. 
C Just one can be specified with NOTARGET. 
C If NORANDOM is defined, then calculations are done for NORANDOM
C randomly-chosed target trees.
C If none of NORANDOM, NOTARGET and ITARGETS is defined, calculations 
C are done for all trees.
C Trees within EDGEDIST m of the plot edges are exempted. 
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /TREESCON/ NOTREES,NOTARGET,ITARGETS,NORANDOM,
     &  EDGEDIST
      REAL DX(MAXT),DY(MAXT)
      INTEGER ITARGETS(MAXT),ITARGETSI(MAXT)
      INTEGER UFILE

C Default values
      NOTREES = 0
      NOTARGET = 0
      NORANDOM = 0
      EDGEDIST = 0
      DO 10 ITAR = 1,MAXT
        ITARGETS(ITAR) = 0
10    CONTINUE

C Read file
      REWIND (UFILE)
      READ (UFILE,TREESCON,IOSTAT = IOERROR)

C Check the namelist was there
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('INPUT ERROR: MISSING TREES IN CONTROL FILE',
     &  IFATAL,IOERROR)

C Set number of trees used in calculations
      IF (NOTREES.EQ.0) NOTREES = NOALLTREES
      IF (NOTREES.GT.MAXTT) THEN
        CALL SUBERROR(
     &    'WARNING: NOTREES IN CONTROL FILE EXCEEDED MAXIMUM',
     &    IWARN,0)
        NOTREES = MAXTT
      END IF

C Set up the array of target tree numbers ITARGETS
C Case 1: One target tree specified
      IF (NOTARGET.GT.0) THEN
        IF (NOTARGET.GT.NOALLTREES)
     &  CALL SUBERROR('INCORRECT TARGET TREE NUMBER',
     &  IFATAL,0)
        NOTARGETS = 1
        ITARGETS(1) = NOTARGET  

C Case 2: Trees to be chosen randomly
      ELSE IF (NORANDOM.GT.0) THEN
        IF (NORANDOM.GT.NOALLTREES)
     &    CALL SUBERROR('TOO MANY TARGET TREES SPECIFIED',
     &    IFATAL,0)
        NOTARGETS = NORANDOM
        DO 20 IRAN = 1,NORANDOM
30        CALL RANDOM(RANVAL)
          RANVAL = RANVAL*REAL(NOALLTREES+1)
          ITREE = NINT(RANVAL)
          IF (ITREE.EQ.0) ITREE = 1
          IF (INEDGES(DX(ITREE),DY(ITREE),XMAX,YMAX,EDGEDIST).LT.0)
     &      GOTO 30          
          ITARGETS(IRAN) = ITREE
20      CONTINUE
          
C Case 3: All trees are to be used, except those in the edges
      ELSE IF (ITARGETS(1).EQ.0) THEN
        NOTARGETS = NOALLTREES
        DO 40 ITAR = 1,NOALLTREES
          IF (INEDGES(DX(ITAR),DY(ITAR),XMAX,YMAX,
     &    EDGEDIST).LT.0) THEN
            NOTARGETS = NOTARGETS - 1
          ELSE          
            ITARGETS(ITAR + NOTARGETS - NOALLTREES) = ITAR
          END IF
40      CONTINUE

      ELSE
C Case 4: A series of target trees is given.
        NOTARGETS = 0
        DO WHILE (ITARGETS(NOTARGETS+1).GT.0)
          NOTARGETS = NOTARGETS + 1
        END DO
      END IF

      NOTREESI = NOTREES
      NOTARGETSI = NOTARGETS
      DO 50 ITAR = 1,MAXT
        ITARGETSI(ITAR) = ITARGETS(ITAR)
50    CONTINUE

      RETURN
      END !ReadConTrees


C**********************************************************************
      SUBROUTINE READXYZ(UFILE,NOALLTREES,X0,Y0,XMAX,YMAX,XSLOPE,YSLOPE,
     &  DX,DY,DZ)
C Read in the X, Y co-ordinates of each crown.
C Calculate the Z co-ordinates of each crown from slope.
C If they are not in the file, then the co-ordinates should be
C calculated from the stocking and the size of the plot.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /XY/ XYCOORDS
      INTEGER UFILE
      REAL DX(MAXT),DY(MAXT),DZ(MAXT),XYCOORDS(MAXT*2)

      REWIND (UFILE)
      READ (UFILE,XY,IOSTAT = IOERROR)

      IF (IOERROR.NE.0) THEN
C Must calculate x, y co-ordinates
        CALL SUBERROR('WARNING: CALCULATING X, Y CO-ORDINATES',
     &  IWARN,IOERROR)
        SPACING = SQRT(XMAX*YMAX/NOALLTREES)
        NOX = XMAX/SPACING + 1
        DO 10 I = 1,NOALLTREES
          DX(I) = (I-1)/NOX * SPACING
          DY(I) = MOD(I-1,NOX) * SPACING
10      CONTINUE

      ELSE
C Read in x, y co-ordinates
        DO 20 I = 1,NOALLTREES
          DX(I) = XYCOORDS(2*I - 1) - X0
          DY(I) = XYCOORDS(2*I) - Y0
20      CONTINUE
      END IF

C Calculate the z co-ordinates from slopes
      ZADD=0.0
      IF (XSLOPE.LT.0.0) ZADD=XMAX*SIN(XSLOPE)
      IF (YSLOPE.LT.0.0) ZADD=ZADD+YMAX*SIN(YSLOPE)
      ZADD=ABS(ZADD)
C  X and Y distances are measured on the slope so the height is based
C  on the sin(slope).
      DO 30 I = 1,NOALLTREES
        DZ(I) = DX(I)*SIN(XSLOPE) + DY(I)*SIN(YSLOPE) + ZADD
30    CONTINUE

      RETURN
      END ! ReadXY


C**********************************************************************
      SUBROUTINE GETLEAFAREA(UFILE,IFLUSH,DT1I,DT2I,DT3I,DT4I,
     &  EXPTIMEI,APP,EXPAN,NOALLTREES,NOLADATES,DATESLA,FLT)
C A subroutine to read in leaf area array.
C First tries to calculate it from phenology (Wang, Rey & Jarvis 1998).
C Otherwise reads in array directly.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL FLT(MAXTDATE,MAXT)
      INTEGER DATESLA(MAXTDATE)
      CHARACTER*10 FLUSHDATE
      INTEGER UFILE

      NAMELIST /PHENOLOGY/ FLUSHDATE,DT1,DT2,DT3,DT4,
     &  EXPTIME,MAXLEAVES,SIZELEAF

      REWIND(UFILE)
      IFLUSH = 0
      READ(UFILE,PHENOLOGY,IOSTAT=IOERROR)
      IF (IOERROR.EQ.0) THEN
        CALL SUBERROR('LEAF AREA FROM PHENOLOGY:',IWARN,0)
        NOLADATES = 0
        IFLUSH = IDATE50(FLUSHDATE)
        DT1I = DT1
        DT2I = DT2
        DT3I = DT3
        DT4I = DT4
        EXPTIMEI = EXPTIME
C Rate of leaf appearance (leaf per day)
        APP = 2. * MAXLEAVES / (DT1*(2*DT2-DT1))
C Rate of leaf area expansion (m2 leaf per day)
        EXPAN = SIZELEAF / EXPTIME
      ELSE
        CALL READTREEARRAY(UFILE,5,NOALLTREES,NOLADATES,DATESLA,FLT)
      END IF

      RETURN
      END ! GetLeafArea


C**********************************************************************
        SUBROUTINE PHENOL(IDAY,ISTART,IFLUSH,DT1,DT2,DT3,DT4,
     &    EXPTIME,APP,EXPAN,STOCKING,NOTREES,FOLT,TOTLAI,NEWCANOPY)
C Implements phenology routine of Wang et al (1998) GCB to appear.
C Parameters:
C   IFLUSH = date of flushing (in days-since-1950 format)
C   DT1 = time from flushing to end of first flush (d)
C   DT2 = time from flushing to end of second flush (d)
C   DT3 = time from flushing to beginning of senescencs (d)
C   DT4 = time from flushing to leaf fall (d)
C   EXPTIME = time a leaf takes to expand (d)
c   APP = rate of leaf appearance (leaves / day)
C   EXPAN = rate of leaf expansion (m2 leaf / day)
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL FOLT(MAXT)

C Apply 7-part formula in Wang et al 1998
      I = IDAY + ISTART - IFLUSH
      IF (I.LE.0) THEN
        FOLT(1) = 0.0
      ELSE IF (I.LE.EXPTIME) THEN
        FOLT(1) = APP*EXPAN*(I**3)/6.
        NEWCANOPY = 1
      ELSE IF (I.LE.DT1) THEN
        FOLT(1) = APP*EXPAN*(I**3 - (I - EXPTIME)**3)/6.
        NEWCANOPY = 1
      ELSE IF (I.LE.(DT1+EXPTIME)) THEN
        FOLT(1) = APP*EXPAN*(DT1*(3.*I**2 + DT1**2 - 3.*DT1*I) - 
     &    (I - EXPTIME)**3)/6.
        NEWCANOPY = 1
      ELSE IF (I.LE.DT2) THEN
        FOLT(1) = APP*EXPAN*DT1*EXPTIME*(I-0.5*(DT1+EXPTIME))
        NEWCANOPY = 1
      ELSE IF (I.LE.(DT2+EXPTIME)) THEN
        FOLT(1) = 0.5*APP*EXPAN*DT1*(2.*DT2*I-DT2**2-DT1*EXPTIME
     &    - (EXPTIME - I)**2)
        NEWCANOPY = 1
      ELSE IF (I.LE.DT3) THEN
        FOLT(1) = APP*EXPAN*DT1*EXPTIME*(DT2-0.5*DT1)
      ELSE IF (I.LE.DT4) THEN
        FOLT(1) = (DT4-I)*(DT2-0.5*DT1)*APP*EXPAN*DT1*EXPTIME
     &    / (DT4 - DT3)
        NEWCANOPY = 1
      ELSE 
        FOLT(1) = 0.0
      END IF

      DO 10 ITREE = 2,NOTREES
        FOLT(ITREE) = FOLT(ITREE)
10    CONTINUE
      TOTLAI = FOLT(1)*STOCKING

      RETURN
      END !Phenol

C**********************************************************************
      SUBROUTINE READTREEARRAY(UFILE,NARRAY,NOALLTREES,
     &  NDATE,IDATES,VALUESI)
C Read in an array of tree parameters from UFILE.
C NARRAY is the number of the array to be read (1 = RADX; 2 = RADY;
C 3 = HTCROWN; 4 = HTTRUNK; 5 = AREALEAF; 6 = DIAM; 7 = LEAFN (for understorey))
C Either read in values for all trees (NOALLTREES, up to MAXT trees) or
C read in average values. All values can be given for a series of dates.
C**********************************************************************

      INCLUDE 'MAESTCOM'

	  REAL VALUES(MAXTDATE*MAXT),VALUESI(MAXTDATE,MAXT)
      CHARACTER*8 DATES(MAXTDATE)
      INTEGER IDATES(MAXTDATE)
      INTEGER UFILE

      NAMELIST /INDIVRADX/ NODATES, DATES, VALUES
      NAMELIST /ALLRADX/ NODATES, DATES, VALUES
      NAMELIST /INDIVRADY/ NODATES, DATES, VALUES
      NAMELIST /ALLRADY/ NODATES, DATES, VALUES
      NAMELIST /INDIVHTCROWN/ NODATES, DATES, VALUES
      NAMELIST /ALLHTCROWN/ NODATES, DATES, VALUES
      NAMELIST /INDIVHTTRUNK/ NODATES, DATES, VALUES
      NAMELIST /ALLHTTRUNK/ NODATES, DATES, VALUES
      NAMELIST /INDIVLAREA/ NODATES, DATES, VALUES
      NAMELIST /ALLLAREA/ NODATES, DATES, VALUES
      NAMELIST /INDIVDIAM/ NODATES, DATES, VALUES
      NAMELIST /ALLDIAM/ NODATES, DATES, VALUES
      NAMELIST /INDIVFOLN/ NODATES, DATES, VALUES
      NAMELIST /ALLFOLN/ NODATES, DATES, VALUES

C Default values
      NODATES = 0
      DATES(1) = '01/01/50'  ! If only one date, doesn't matter what it is
      DO 10 I = 1,MAXTDATE*MAXT
        VALUES(I) = -1.0
10    CONTINUE

C Try to read arrays for individual trees first
      REWIND(UFILE)
      IF (NARRAY.EQ.1) THEN
        READ(UFILE,INDIVRADX,IOSTAT=IOERROR)
        CALL SUBERROR('X RADII ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.2) THEN
        READ(UFILE,INDIVRADY,IOSTAT=IOERROR)
        CALL SUBERROR('Y RADII ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.3) THEN
        READ(UFILE,INDIVHTCROWN,IOSTAT=IOERROR)
        CALL SUBERROR('CROWN HEIGHT ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.4) THEN
        READ(UFILE,INDIVHTTRUNK,IOSTAT=IOERROR)
        CALL SUBERROR('TRUNK HEIGHT ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.5) THEN
        READ(UFILE,INDIVLAREA,IOSTAT=IOERROR)
        CALL SUBERROR('LEAF AREA ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.6) THEN
        READ(UFILE,INDIVDIAM,IOSTAT=IOERROR)
        CALL SUBERROR('DIAMETER ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.7) THEN
        READ(UFILE,INDIVFOLN,IOSTAT=IOERROR)
        CALL SUBERROR('UNDERSTOREY N ARRAY:',IWARN,0)
      END IF

      IF (IOERROR.NE.-1) THEN
C Process arrays, if read in
        IF (NODATES.GT.MAXTDATE) THEN
          CALL SUBERROR(
     &    'WARNING: TOO MANY DATES: SOME DATA LOST',
     &    IWARN,IOERROR)
          NODATES = MAXTDATE
        END IF
        INDEX = 1
        DO 20 IDATE = 1,NODATES
          IDATES(IDATE) = IDATE50(DATES(IDATE))
20      CONTINUE
        DO 30 ITREE = 1,NOALLTREES
          DO 30 IDATE = 1,NODATES
            IF (VALUES(INDEX).LT.0)
     &        CALL SUBERROR('MISSING DATA',IFATAL,0)
            VALUESI(IDATE,ITREE) = VALUES(INDEX)
            INDEX = INDEX+1
30      CONTINUE

      ELSE
C Read in values for all trees
        REWIND(UFILE)
        IF (NARRAY.EQ.1) THEN
          READ(UFILE,ALLRADX,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.2) THEN
          READ(UFILE,ALLRADY,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.3) THEN
          READ(UFILE,ALLHTCROWN,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.4) THEN
          READ(UFILE,ALLHTTRUNK,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.5) THEN
          READ(UFILE,ALLLAREA,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.6) THEN
          READ(UFILE,ALLDIAM,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.7) THEN
          READ(UFILE,ALLFOLN,IOSTAT=IOERROR)
        END IF

        IF ((IOERROR.NE.0).AND.(NARRAY.LT.6))   ! Missing diam, leaf N is OK
     &    CALL SUBERROR('ERROR: MISSING DATA',IFATAL,IOERROR)
        IF (NODATES.GT.MAXTDATE) THEN
          CALL SUBERROR(
     &    'WARNING: TOO MANY DATES: SOME DATA LOST',
     &    IWARN,IOERROR)
          NODATES = MAXTDATE
        END IF
C Assign arrays to data
        DO 40 IDATE = 1,NODATES
          IDATES(IDATE) = IDATE50(DATES(IDATE))
40      CONTINUE
        DO 50 ITREES = 1,NOALLTREES
          DO 50 IDATE = 1,NODATES
            VALUESI(IDATE,ITREES) = VALUES(IDATE)
50      CONTINUE

      END IF

      NDATE = NODATES
      RETURN
      END !ReadTreeArray


C**********************************************************************
      SUBROUTINE CALCLAI(NOLADATES,FLT,NOALLTREES,
     &  XMAX,YMAX,XSLOPE,YSLOPE,TOTLAI)
C Calculate total LAI of the plot - based on horizontal plot area
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL FLT(MAXTDATE,MAXT),TOTLAI(MAXTDATE)

      GROUNDAREA = XMAX*COS(XSLOPE)*YMAX*COS(YSLOPE)
      DO 10 IDATE = 1,NOLADATES
        TOTLAI(IDATE) = 0.0
        DO 20 ITREE = 1,NOALLTREES
          TOTLAI(IDATE) = TOTLAI(IDATE) + FLT(IDATE,ITREE)
20      CONTINUE
        TOTLAI(IDATE) = TOTLAI(IDATE)/GROUNDAREA
10    CONTINUE

      RETURN
      END ! CalcLAI


C**********************************************************************
      SUBROUTINE TREEINTERP(IDAY,ISTART,NODATES,IDATEARR,PARAMTABLE,
     &  NOTREES,NEWCANOPY,PARAMS)
C Interpolate crown dimensions for this date from the
C date arrays read in from the trees file.
C Parameter NEWCANOPY indicates whether crown has changed - in which
C case grid points need to be reassigned.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      INTEGER IDATEARR(MAXTDATE)
      REAL PARAMTABLE(MAXTDATE,MAXTT)
      REAL PARAMS(MAXTT)

      IDATE = IDAY + ISTART
      IF (IDAY.EQ.0) NEWCANOPY = 1

C If no dates, take 0.0
      IF (NODATES.EQ.0) THEN
        DO 10 ITREE = 1,NOTREES
          PARAMS(ITREE) = 0.0
10      CONTINUE

C If only one date, or before the starting date, take first value
      ELSE IF ((NODATES.EQ.1).OR.IDATE.LE.IDATEARR(1)) THEN
        DO 20 ITREE = 1,NOTREES
          PARAMS(ITREE) = PARAMTABLE(1,ITREE)
20      CONTINUE

C If after the final date, take last value
      ELSE IF (IDATE.GT.IDATEARR(NODATES)) THEN
        DO 30 ITREE = 1,NOTREES
          PARAMS(ITREE) = PARAMTABLE(NODATES,ITREE)
30      CONTINUE

C Otherwise have to interpolate
      ELSE
        INDEX = 1
        DO WHILE (IDATE.GT.IDATEARR(INDEX))
          INDEX = INDEX + 1
        END DO
        SLOPE = REAL((IDATE - IDATEARR(INDEX-1)))/
     &      REAL((IDATEARR(INDEX) - IDATEARR(INDEX-1)))
        IF ((SLOPE.LT.-1E-6).OR.(SLOPE.GT.1E-6))
     &      NEWCANOPY = 1
        DO 40 ITREE = 1,NOTREES
            Y1 = PARAMTABLE(INDEX-1, ITREE)
            Y2 = PARAMTABLE(INDEX, ITREE)
            PARAMS(ITREE) = Y1 + SLOPE*(Y2 - Y1)
40      CONTINUE

      END IF

      RETURN
      END !InterpTree


C**********************************************************************
      SUBROUTINE INTERPOLATEP(IDAY, ISTART,
     &  NOJDATES,DATESJ,JMAXTABLE,
     &  NOVDATES,DATESV,VCMAXTABLE,
     &  NORDATES,DATESRD,RDTABLE,
     &  NOSLADATES,DATESSLA,SLATABLE,
     &  NOADATES,DATESA,AJQTABLE,
     &  NOLAY,NOAGEP,
     &  JMAX25,VCMAX25,RD0,SLA,AJQ)
C Controls the calling of the interpolation routines to get daily values
C of Jmax, Vcmax, SLA, Rd.
C**********************************************************************

      INCLUDE 'MAESTCOM'

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

C Interpolate to get daily values of physiological params
      CALL PHYINTERP(IDAY+ISTART,NOJDATES,DATESJ,JMAXTABLE,
     &  NOLAY,NOAGEP,JMAX25)
      CALL PHYINTERP(IDAY+ISTART,NOVDATES,DATESV,VCMAXTABLE,
     &  NOLAY,NOAGEP,VCMAX25)
      CALL PHYINTERP(IDAY+ISTART,NORDATES,DATESRD,RDTABLE,
     &  NOLAY,NOAGEP,RD0)
      IF (NOSLADATES.GT.0) 
     &  CALL PHYINTERP(IDAY+ISTART,NOSLADATES,DATESSLA,SLATABLE,
     &  NOLAY,NOAGEP,SLA)
      CALL PHYINTERP(IDAY+ISTART,NOADATES,DATESA,AJQTABLE,
     &  NOLAY,NOAGEP,AJQ)

      RETURN
      END !InterpolateP


C**********************************************************************
      SUBROUTINE INTERPOLATET(IDAY, ISTART,
     &  NOXDATES,DATESX,RXTABLE,
     &  NOYDATES,DATESY,RYTABLE,
     &  NOZDATES,DATESZ,RZTABLE,
     &  NOTDATES,DATEST,ZBCTABLE,
     &  NODDATES,DATESD,DIAMTABLE,
     &  NOLADATES,DATESLA,FOLTABLE,TOTLAITABLE,NOTREES,
     &  RX,RY,RZ,ZBC,FOLT,TOTLAI,DIAM,STOCKING,
     &  IFLUSH,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN,
     &  NEWCANOPY)
C Controls the calling of the interpolation routines to get daily values
C of crown heights and radii and leaf area.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      REAL RX(MAXTT),RY(MAXTT),RZ(MAXTT),ZBC(MAXTT),FOLT(MAXTT)
      REAL RXTABLE(MAXTDATE,MAXTT),RYTABLE(MAXTDATE,MAXTT)
      REAL RZTABLE(MAXTDATE,MAXTT),ZBCTABLE(MAXTDATE,MAXTT)
      REAL FOLTABLE(MAXTDATE,MAXTT),TOTLAITABLE(MAXTDATE)
      REAL DIAMTABLE(MAXTDATE,MAXTT),DIAM(MAXTT)
      INTEGER DATESX(MAXTDATE),DATESY(MAXTDATE),DATESZ(MAXTDATE)
      INTEGER DATEST(MAXTDATE),DATESLA(MAXTDATE),DATESD(MAXTDATE)

C Interpolate to get daily values of crown dimensions
      NEWCANOPY = 0
      CALL TREEINTERP(IDAY,ISTART,NOXDATES,DATESX,RXTABLE,NOTREES,
     &  NEWCANOPY,RX)

      CALL TREEINTERP(IDAY,ISTART,NOYDATES,DATESY,RYTABLE,NOTREES,
     &  NEWCANOPY,RY)

      CALL TREEINTERP(IDAY,ISTART,NOZDATES,DATESZ,RZTABLE,NOTREES,
     &  NEWCANOPY,RZ)

      CALL TREEINTERP(IDAY,ISTART,NOTDATES,DATEST,ZBCTABLE,NOTREES,
     &  NEWCANOPY,ZBC)

      IF (NOLADATES.GT.0) THEN
        CALL TREEINTERP(IDAY,ISTART,NOLADATES,DATESLA,FOLTABLE,NOTREES,
     &    NEWCANOPY,FOLT)
        CALL TREEINTERP(IDAY,ISTART,NOLADATES,DATESLA,TOTLAITABLE,1,
     &    NEWCANOPY,TOTLAI)
      ELSE
        CALL PHENOL(IDAY,ISTART,IFLUSH,DT1,DT2,DT3,DT4,
     &    EXPTIME,APP,EXPAN,STOCKING,NOTREES,FOLT,TOTLAI,NEWCANOPY)
      END IF

      IF (NODDATES.GT.0)
     &  CALL TREEINTERP(IDAY,ISTART,NODDATES,DATESD,DIAMTABLE,NOTREES,
     &  NEWCANOPY,DIAM)

      RETURN
      END !InterpolateT


C**********************************************************************
      SUBROUTINE OUTPUTLAY(UFILE,FOLLAY,JMAX25,VCMAX25,NOLAY)
C Daily output to layer flux file.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      INTEGER UFILE
      REAL FOLLAY(MAXLAY)
      REAL JMAX25(MAXLAY,MAXC),VCMAX25(MAXLAY,MAXC)

      WRITE(UFILE,*) 'LEAF AREA OF TARGET TREE IN EACH LAYER IN M2'
      WRITE(UFILE,500) (FOLLAY(I),I=1,NOLAY)
      WRITE(UFILE,*)
      WRITE(UFILE,*) 'JMAX (CURRENT) IN EACH LAYER IN UMOL M-2 S-1'
      WRITE(UFILE,500) (JMAX25(I,1),I=1,NOLAY)
      WRITE(UFILE,*)
      WRITE(UFILE,*) 'VCMAX (CURRENT) IN EACH LAYER IN UMOL M-2 S-1'
      WRITE(UFILE,500) (VCMAX25(I,1),I=1,NOLAY)
      WRITE(UFILE,*)

500   FORMAT(10(F10.5,1X))

      RETURN
      END !OutputLay


C**********************************************************************
      SUBROUTINE OUTPUTHIST(UFILE,ITREE,HISTO,BINSIZE)
C Write PAR histogram to file.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL HISTO(MAXHISTO)
      INTEGER UFILE

      WRITE (UFILE,500) ITREE
      WRITE (UFILE,510) 
      DO 10 IBIN = 1,MAXHISTO
        WRITE (UFILE,520) (IBIN-0.5)*BINSIZE,HISTO(IBIN)
10    CONTINUE
500   FORMAT ('TREE NO: ',I6)
510   FORMAT ('  PAR:         FREQUENCY (M^2.HR): ')
520   FORMAT (F8.2,1X,F12.6)

      RETURN
      END !OutputHist


C**********************************************************************
	SUBROUTINE OUTPUTSOIL(UFILE,IDAY,PPT,TOTH2O,DRAINAGE,
     &	CANWATER,STEMWATER,SOILWATER)
C Write out components of soil water balance
C Need to think about how to do this for more than one tree .. 
C**********************************************************************

      INCLUDE 'maestcom'

      INTEGER UFILE
	  REAL PPT(KHRS), SOILWATER(MAXS)

	  TOTPPT = 0.0
	  DO 10 I = 1,KHRS
	    TOTPPT = PPT(I) + TOTPPT
10	  CONTINUE
      TOTSW = 0.0
	  DO 20 I = 1,MAXS
	    TOTSW = TOTSW + SOILWATER(I)
20	  CONTINUE

      WRITE (UFILE,500) IDAY,TOTPPT,TOTH2O,DRAINAGE,
     &	CANWATER,STEMWATER,TOTSW
500   FORMAT (I5,6(1X,F12.6))

	  RETURN
	  END !OutputSoil


C**********************************************************************
	  SUBROUTINE OUTPUTSOILHR(UFILE,IDAY,IHOUR,PPT,
     &	THRUFALL,STEMFLOW,CANOPYDRAIN,STEMEVAPN,
     &	FH2OT,FH2OE,TOTDRAIN,
     &	CANWATER,STEMWATER,SOILWATER)
C**********************************************************************

	  INCLUDE 'maestcom'

      INTEGER UFILE
	  REAL SOILWATER(MAXS)

      TOTSW = 0.0
	  DO 20 I = 1,MAXS
	    TOTSW = TOTSW + SOILWATER(I)
20	  CONTINUE

      WRITE (UFILE,500) IDAY,IHOUR,PPT,
     &	THRUFALL,STEMFLOW,CANOPYDRAIN,STEMEVAPN,
     &	FH2OT,FH2OE,TOTDRAIN,
     &	CANWATER,STEMWATER,TOTSW
500   FORMAT (2(I5,1X),11(F8.2,1X))

	RETURN
	END !OutputSoil


C**********************************************************************
C THIS FUNCTION NO LONGER USED 14/1/98!!!
      FUNCTION GETTREENO(ICALL,NOTARGET,NORANDOM,NOALLTREES,
     &  DXT,DYT,XMAX,YMAX,EDGEDIST,
     &  ITREE)
C Either return the number of the target tree, or a randomly-chosen
C tree number, or the next tree number (if doing them all). 
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL DXT(MAXT),DYT(MAXT)

      GETTREENO = 1.0

C Target tree specified
      IF (NOTARGET.GT.0) THEN
        IF (ICALL.EQ.1) THEN
          ITREE = NOTARGET
        ELSE
          GETTREENO = -1.0
        END IF

C Choosing trees randomly
      ELSE IF (NORANDOM.GT.0) THEN
        IF (ICALL.GT.NORANDOM) THEN
          GETTREENO = -1.0
        ELSE        
10        CALL RANDOM(RANVAL)
          RANVAL = RANVAL*REAL((NOALLTREES+1))
          ITREE = NINT(RANVAL)
          IF (ITREE.EQ.0) ITREE = 1
          IF (INEDGES(DXT(ITREE),DYT(ITREE),XMAX,YMAX,EDGEDIST).LT.0)
     &      GOTO 10          
        END IF

C Working through all the trees
      ELSE
20      ITREE = ITREE + 1
        IF (ITREE.GT.NOALLTREES) THEN
          GETTREENO = -1.0
        ELSE IF (INEDGES(DXT(ITREE),DYT(ITREE),XMAX,YMAX,EDGEDIST).
     &  LT.0) THEN
          GOTO 20
        END IF          
      END IF

      RETURN
      END

C**********************************************************************
C THIS FUNCTION NO LONGER USED 14/1/98!!!
      FUNCTION NOTREEGET(NOTARGET,NORANDOM,NOALLTREES,ITREE)
C Either return the number of the target tree, or a randomly-chosen
C tree number, or the next tree number (if doing them all). 
C**********************************************************************

      IF (NOTARGET.GT.0) THEN
        NOTREEGET = NOTARGET
      ELSE IF (NORANDOM.GT.0) THEN
        CALL RANDOM(RANVAL)
        RANVAL = RANVAL*REAL((NOALLTREES+1))
        NOTREEGET = NINT(RANVAL)
        IF (NOTREEGET.EQ.0) NOTREEGET = 1
      ELSE
        NOTREEGET = ITREE
      END IF

      RETURN
      END


C**********************************************************************
      FUNCTION INEDGES(DX,DY,XMAX,YMAX,EDGEDIST)
C Checks to see if the chosen tree is within EDGEDIST m of the edge 
C of the plot. Returns -1 if yes and 1 if no.
C**********************************************************************

      IF ((DX.LT.EDGEDIST).OR.(DY.LT.EDGEDIST).OR.
     &  (XMAX-DX.LT.EDGEDIST).OR.(YMAX-DY.LT.EDGEDIST)) THEN
        INEDGES = -1
      ELSE
        INEDGES = 1
      END IF

      RETURN
      END ! InEdges


