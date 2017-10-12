C------------------------------------------------------------------------
C CONTENTS:
C
C OPENINFILE - Open the files from which data are to come
C OPENOUTFILE - Open the output files for the totals. 
C GETCOMM - Get the commands to control the summing process.
C GETSTOCKING - Calculate stocking from plot details in trees file.
C READDAY - Read in daily values from file and sum.
C INLIST - Find out if given tree is in list of targets.
C ZEROT - set the initial values of daily & hourly total variables to zero.
C OUTPUTTDAY - Output daily totals to file.
C READHR - Read in hourly values from file and sum.
C OUTPUTTHR - Output hourly totals to file.
C READHIST - Read in and sum histograms.
C OUTPUTTHIST - Write total PAR histogram to file.
C CLOSEF - Close files.
C------------------------------------------------------------------------

      PROGRAM SUMTREES

C------------------------------------------------------------------------
C A program to sum the data in files dayflx.dat, hrflux.dat, histo.dat
C to get weighted totals for the stand.
C------------------------------------------------------------------------

      INCLUDE 'MAESTCOM'

C Declare arrays holding sums
      REAL TOTHRAB(KHRS*MAXDAY,3),TOTFCO2(KHRS*MAXDAY)
      REAL TOTFH2OT(KHRS*MAXDAY),TOTFH2OE(KHRS*MAXDAY)
	REAL TOTFH2OCAN(KHRS*MAXDAY)
      REAL TOTFRF(KHRS*MAXDAY),TOTFRG(KHRS*MAXDAY)
      REAL TOTFHT(KHRS*MAXDAY),TOTGSCAN(KHRS*MAXDAY)
	REAL AVTCAN(KHRS*MAXDAY)
      REAL TOTDYAB(MAXDAY,3),TOTDCO2(MAXDAY),TOTDCO2NET(MAXDAY)
      REAL TOTDH2O(MAXDAY),TOTDH2OCAN(MAXDAY)
      REAL TOTDRF(MAXDAY),TOTDHT(MAXDAY)
      REAL PAR(MAXHISTO),HISTO(MAXHISTO)
      INTEGER ITARGETS(MAXT),ICHECK(MAXDAY)
      REAL WEIGHTS(MAXT)
      CHARACTER*80 CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE

C Open input & output files
      CALL OPENINFILE(IOHRLY,IOHIST,IOCON,
     &  CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE)
      CALL OPENOUTFILE(IOHRLY,IOHIST,
     &  CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE)

C Find out which trees are required
      CALL GETSTOCKING(UTREES,STOCKING)
      CALL GETCOMM(IOCON,ITARGETS,NOTARGETS,WEIGHTS)

C Zero the arrays holding the sums
      CALL ZEROT(TOTHRAB,TOTFCO2,TOTFH2OT,TOTFH2OE,TOTFHT,
     &  TOTFH2OCAN,TOTFRF,TOTFRG,TOTGSCAN,
     &  TOTDYAB,TOTDCO2,TOTDCO2NET,TOTDH2O,TOTDHT,
     &  TOTDH2OCAN,TOTDRF,AVTCAN,
     &  HISTO,ICHECK)

C Read the input files and sum
      CALL READDAY(UDAILY,ITARGETS,NOTARGETS,WEIGHTS,
     &  TOTDYAB,TOTDCO2,TOTDCO2NET,TOTDH2O,TOTDHT,
     &  TOTDH2OCAN,TOTDRF,ICHECK,NTREES)

      IF (IOCON.GE.0) NTREES = 1
      CALL OUTPUTTDAY(ICHECK,STOCKING,NTREES,
     &  TOTDYAB,TOTDCO2,TOTDRF,TOTDCO2NET,
     &  TOTDH2O,TOTDH2OCAN,TOTDHT)
    
      IF (IOHRLY.EQ.1) THEN
        CALL READHR(UHRLY,ITARGETS,NOTARGETS,WEIGHTS,
     &    TOTHRAB,TOTFCO2,TOTFH2OT,TOTFH2OE,TOTFHT,TOTGSCAN,
     &    TOTFH2OCAN,TOTFRF,TOTFRG,AVTCAN)
        CALL OUTPUTTHR(ICHECK,STOCKING,NTREES,
     &    TOTHRAB,TOTFCO2,TOTFRF,TOTFRG,
     &    TOTFH2OT,TOTFH2OE,TOTFH2OCAN,TOTFHT,TOTGSCAN,AVTCAN)
      END IF

      IF (IOHIST.EQ.1) THEN
        CALL READHIST(UHIST,ITARGETS,NOTARGETS,WEIGHTS,
     &    PAR,HISTO)
        CALL OUTPUTTHIST(UTHIST,PAR,HISTO,STOCKING,NTREES)
      END IF

      CALL CLOSEF(IOHRLY,IOHIST,IOCON)

      STOP
      END

C------------------------------------------------------------------------


C**********************************************************************
      SUBROUTINE OPENINFILE(IOHRLY,IOHIST,IOCON,
     &  CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE)
C Open the files from which data are to come, checking they are available.
C**********************************************************************

      INCLUDE 'maestcom'

      CHARACTER*80 CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE

C Output file for errors and warnings
      OPEN (UERROR, FILE = 'Maeserr.dat', STATUS = 'UNKNOWN')

      IOHRLY = 0
      IOHIST = 0

C Command file
      OPEN(UTCON,FILE = 'Sumcon.dat',STATUS = 'OLD',
     &  IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        IOCON = -1
      ELSE
        IOCON = UTCON
      END IF

C Original trees file (for stocking)
      OPEN (UTREES, FILE = 'trees.dat', STATUS='OLD')

C Output file with daily fluxes
      OPEN (UDAILY,FILE = 'Dayflx.dat', STATUS='OLD',
     &  IOSTAT = IOERROR)
      IF (IOERROR.NE.0) 
     &  CALL SUBERROR('ERROR: COULDNT FIND DAYFLX.DAT',
     &  IFATAL,0)

C Read titles of input files
990   FORMAT (A80)     ! For reading titles in input files.
      READ (UDAILY, 990) VTITLE
      READ (UDAILY, 990) CTITLE
      READ (UDAILY, 990) TTITLE
      READ (UDAILY, 990) STITLE
      READ (UDAILY, 990) PTITLE
      READ (UDAILY, 990) MTITLE
      READ (UDAILY, 990) UTITLE

C Output file with hourly fluxes (if available).
      OPEN (UHRLY, FILE = 'hrflux.dat', STATUS='OLD',
     &  IOSTAT = IOERROR)
      IF (IOERROR.EQ.0) IOHRLY = 1

C Output file with histogram (if available).
      OPEN (UHIST, FILE = 'histo.dat', STATUS = 'OLD',
     &  IOSTAT = IOERROR)
      IF (IOERROR.EQ.0) IOHIST = 1

      RETURN
      END !OpenInFile


C**********************************************************************
      SUBROUTINE OPENOUTFILE(IOHRLY,IOHIST,
     &  CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE)
C Open the output files for the totals. 
C**********************************************************************

      INCLUDE 'maestcom'

      CHARACTER*80 CTITLE,TTITLE,PTITLE,STITLE,MTITLE,VTITLE,UTITLE

C Output file with daily fluxes
      OPEN (UTDAY,FILE = 'Daytot.dat', STATUS='UNKNOWN')

C Output file with hourly fluxes (if required).
      IF (IOHRLY.GT.0) THEN
        OPEN (UTHR, FILE = 'hrtot.dat', STATUS='UNKNOWN')
      ENDIF

C Output file with histogram (if required).
      IF (IOHIST.EQ.1) THEN
        OPEN (UTHIST, FILE = 'tothist.dat', STATUS = 'UNKNOWN')
      END IF

C Write headings to daily flux files.
      WRITE (UTDAY, 990) VTITLE
      WRITE (UTDAY, 990) CTITLE
      WRITE (UTDAY, 990) TTITLE
      WRITE (UTDAY, 990) STITLE
      WRITE (UTDAY, 990) PTITLE
      WRITE (UTDAY, 990) MTITLE
      WRITE (UTDAY, 990) UTITLE
      WRITE (UTDAY, 990) ' '
      WRITE (UTDAY,601)
      WRITE (UTDAY,602)
      WRITE (UTDAY,603)
      WRITE (UTDAY,604)
      WRITE (UTDAY,605)
      WRITE (UTDAY,606)
      WRITE (UTDAY,607)
      WRITE (UTDAY,608)
      WRITE (UTDAY,609)
      WRITE (UTDAY,610)
      WRITE (UTDAY, 990) ' '
      WRITE (UTDAY,611)

C Comments to hourly output file (if required).
      if (IOHRLY.gt.0) then
        WRITE (UTHR, 990) VTITLE
        WRITE (UTHR, 990) CTITLE
        WRITE (UTHR, 990) TTITLE
        WRITE (UTHR, 990) STITLE
        WRITE (UTHR, 990) PTITLE
        WRITE (UTHR, 990) MTITLE
        WRITE (UTHR, 990) UTITLE
        WRITE (UTHR, 990) ' '
        WRITE (UTHR,301)
        WRITE (UTHR,302)
        WRITE (UTHR,303)
        WRITE (UTHR,304)
        WRITE (UTHR,305)
        WRITE (UTHR,306)
        WRITE (UTHR,307)
        WRITE (UTHR,308)
        WRITE (UTHR,309)
        WRITE (UTHR,310)
        WRITE (UTHR,311)
        WRITE (UTHR,312)
        WRITE (UTHR,313)
        WRITE (UTHR,314)
        WRITE (UTHR, 990) ' '
        WRITE (UTHR,315)
      END IF

990   FORMAT (A80)

601   format('DOY:      day of the year')
602   format('absPAR:   absorbed PAR              MJ m-2 d-1')
603   format('absNIR:   absorbed NIR              MJ  m-2 d-1')
604   format('absTherm: absorbed thermal          MJ m-2 d-1')
605   format('totPs: gross photosynthesis         mol m-2 d-1')
606   format('totRf: daily foliar respiration     mol m-2 d-1')
607   format('netPs: photosyn. net of foliar resp   mol m-2 d-1')
608   format('totLE1: daily transpiration          mol H2O m-2 d-1')
609   format('totLE2: daily transpirn (CANOPY calc) mol H2O m-2 d-1')
610   format('totH:  daily sensible heat flux     MJ m-2 d-1')
611   format('Columns: DOY, absPAR, absNIR, absTherm,',
     &  ' totPs, totRf, netPs,',
     +  ' totLE1, totLE2, totH')

301   format('DOY:   day of the year')
302   format('Hour:  hour of the day')
303   format('hrPAR: absorbed PAR              umol m-2 s-1')
304   format('hrNIR: absorbed NIR              W m-2')
305   format('hrTHM: absorbed thermal          W m-2')
306   format('hrPS: photosynthesis (net of leaf resp) umol m-2 s-1')
307   format('hrRf:  hourly leaf respiration   umol m-2 s-1')
308   format('hrRmW: hourly woody Rm           umol m-2 s-1')
309   format('hrLE:  hourly transpiration      mmol m-2 s-1')
310   format('hrEV:  hourly evaporation        mmol m-2 s-1')
311   format('LECAN: hourly transpiration: CANOPY calc      ')
312   FORMAT('Gscan: canopy stomatal conductance: mol CO2 m-1 s-1')
313   format('hrH:   hourly sensible heat flux  J m-2 s-1')
314	FORMAT('TCAN: Average foliage temperature deg C')
315   format('Columns: DOY, HOUR, hrPAR, hrNIR, hrTHM,',
     +  ' hrPs, hrRf, hrRmW, hrLE, hrEV,',
     +  ' LECAN, Gscan, hrH, TCAN')

      RETURN
      END !OpenOutFile


C**********************************************************************
      SUBROUTINE GETCOMM(IOCON,ITARGETSI,NOTARGETS,WEIGHTSI)
C Get the commands to control the summing process.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      NAMELIST /TARGETS/ ITARGETS,WEIGHTS
      REAL WEIGHTS(MAXT),WEIGHTSI(MAXT)
      INTEGER ITARGETS(MAXT),ITARGETSI(MAXT)

C Default values
      WEIGHTS(1) = -1.0
      TOTWEIGHT = 0.0
      DO 10 ITAR = 1,MAXT
        ITARGETS(ITAR) = 0
10    CONTINUE

      IF (IOCON.LT.0) THEN
        NOTARGETS = -1
        WEIGHTSI(1) = 1.0
        RETURN
      END IF

C Read file
      REWIND (IOCON)
      READ (IOCON,TARGETS,IOSTAT = IOERROR)

C Check the namelist was there
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('INPUT ERROR: MISSING COMMANDS IN CONTROL FILE',
     &  IFATAL,IOERROR)

C Work out how many target trees specified
      NOTARGETS = 0
      DO WHILE (ITARGETS(NOTARGETS+1).GT.0)
        NOTARGETS = NOTARGETS + 1
      END DO
C Assign equal weights, if weights not given
      IF (WEIGHTS(1).LT.0.0) THEN 
        DO 80 ITAR = 1,NOTARGETS
          WEIGHTSI(ITAR) = 1./NOTARGETS
80      CONTINUE
      ELSE
C Or, read in given weights & check they sum to 1
        DO 60 ITAR = 1,MAXT
          WEIGHTSI(ITAR) = WEIGHTS(ITAR)
          TOTWEIGHT = TOTWEIGHT + WEIGHTS(ITAR)
60      CONTINUE
        IF ((TOTWEIGHT.LE.0.99).OR.(TOTWEIGHT.GT.1.01)) THEN
        CALL SUBERROR(
     &    'ERROR: WEIGHTS FOR TARGET TREES MUST SUM TO 1',
     &    IFATAL,0)
        END IF
      END IF

      DO 50 ITAR = 1,MAXT
        ITARGETSI(ITAR) = ITARGETS(ITAR)
50    CONTINUE

      RETURN
      END !GetComm


C**********************************************************************
      SUBROUTINE GETSTOCKING(UFILE,STOCKING)
C Calculates stocking from plot details in trees file.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /PLOT/ XMAX,YMAX,NOTREES,XSLOPE,YSLOPE,BEARING,SHADEHT,
     &  X0,Y0
      INTEGER UFILE

      REWIND (UFILE)
      READ (UFILE, PLOT, IOSTAT = IOERROR)
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('ERROR READING PLOT DETAILS',IFATAL,IOERROR)
      STOCKING = NOTREES/(XMAX*YMAX)

      RETURN
      END !GetStocking


C**********************************************************************
      SUBROUTINE READDAY(UFILE,ITARGETS,NOTARGETS,WEIGHTS,
     &  TOTDYAB,TOTDCO2,TOTDCO2NET,TOTDH2O,TOTDHT,
     &  TOTDH2OCAN,TOTDRF,ICHECK,NTREES)
C Read in daily values from file and sum.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      REAL TOTDYAB(MAXDAY,3),TOTDCO2(MAXDAY),TOTDCO2NET(MAXDAY)
      REAL TOTDH2O(MAXDAY),TOTDH2OCAN(MAXDAY)
      REAL TOTDRF(MAXDAY)
      REAL TOTDHT(MAXDAY),TDYAB(3)
      REAL WEIGHTS(MAXT)
      INTEGER ITARGETS(MAXT),ICHECK(MAXDAY)
      INTEGER UFILE
      CHARACTER*80 TMPSTR

C Skip comments at start of file
990   FORMAT (A80)
      REWIND (UFILE)
      DO 10 I = 1,21
        READ (UFILE,990) TMPSTR
10    CONTINUE

      IPTREE = 0
      NTREES = 0

C Read one day at a time. (Program returns to here).
500   FORMAT (I7,1X,I4,1X,9(F12.5,1X))
20    READ (UFILE,500,IOSTAT=IOERROR) 
     &  IDAY,ITREE,TDYAB(1),TDYAB(2),TDYAB(3),
     &  TOTCO2,TOTRESPF,CO2NET,
     &  TOTH2O,TOTH2OCAN,TOTHFX
C If reached the end of the file, return.
      IF (IOERROR.NE.0) RETURN
C Is this a new tree? If so check if in list.
      IF (ITREE.NE.IPTREE) THEN
        ITAR = INLIST(ITREE,ITARGETS,NOTARGETS)
        NTREES = NTREES + 1
        IPTREE = ITREE
      END IF
C If the tree was in the list then add numbers to the totals.
      IF (ITAR.GT.0) THEN
        DO 50 I = 1,3
          TOTDYAB(IDAY,I) = TOTDYAB(IDAY,I) + TDYAB(I)*WEIGHTS(ITAR)
50      CONTINUE
        TOTDCO2(IDAY) = TOTDCO2(IDAY) + TOTCO2*WEIGHTS(ITAR)
        TOTDRF(IDAY) = TOTDRF(IDAY) + TOTRESPF*WEIGHTS(ITAR)
        TOTDCO2NET(IDAY) = TOTDCO2NET(IDAY) + CO2NET*WEIGHTS(ITAR)
        TOTDH2O(IDAY) = TOTDH2O(IDAY) + TOTH2O*WEIGHTS(ITAR)
        TOTDH2OCAN(IDAY) = TOTDH2OCAN(IDAY) + TOTH2OCAN*WEIGHTS(ITAR)
        TOTDHT(IDAY) = TOTDHT(IDAY) + TOTHFX*WEIGHTS(ITAR)
C ICHECK is to know whether this day had values or not.
        ICHECK(IDAY) = ICHECK(IDAY) + 1
      END IF
      GOTO 20

      RETURN
      END !ReadDay


C**********************************************************************
      INTEGER FUNCTION INLIST(ITREE,ITARGETS,NOTARGETS)
C Find out if given tree is in list of targets.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      INTEGER ITARGETS(MAXT)

      INLIST = 0
      IF (NOTARGETS.LT.0) THEN
        INLIST = 1
      ELSE
        DO 10 I = 1,NOTARGETS
          IF (ITARGETS(I).EQ.ITREE) INLIST = I
10      CONTINUE
      END IF

      RETURN
      END !InList


C**********************************************************************
      SUBROUTINE ZEROT(
     &  TOTHRAB,TOTFCO2,TOTFH2OE,TOTFH2OT,TOTFHT,
     &  TOTFH2OCAN,TOTFRF,TOTFRG,
     &  TOTGSCAN,
     &  TOTDYAB,TOTDCO2,TOTDCO2NET,TOTDH2O,TOTDHT,TOTDH2OCAN,TOTDRF,
     &  AVTCAN,HISTO,ICHECK)
C This is subroutine to set the initial values of daily & hourly 
C total variables (for all trees) to zero.
C**********************************************************************

      INCLUDE 'maestcom'
      REAL TOTHRAB(KHRS*MAXDAY,3),TOTFCO2(KHRS*MAXDAY)
      REAL TOTFH2OT(KHRS*MAXDAY),TOTFH2OE(KHRS*MAXDAY)
	REAL TOTFH2OCAN(KHRS*MAXDAY)
      REAL TOTFRF(KHRS*MAXDAY),TOTFRG(KHRS*MAXDAY)
      REAL TOTFHT(KHRS*MAXDAY),TOTGSCAN(KHRS*MAXDAY)
      REAL TOTDYAB(MAXDAY,3),TOTDCO2(MAXDAY),TOTDCO2NET(MAXDAY)
      REAL TOTDH2O(MAXDAY),TOTDH2OCAN(MAXDAY)
      REAL TOTDRF(MAXDAY),TOTDHT(MAXDAY),AVTCAN(KHRS*MAXDAY)
      REAL HISTO(MAXHISTO)
      INTEGER ICHECK(MAXDAY)

      DO 100 IDY = 1,MAXDAY

        DO 10 IHR = 1,KHRS
          I = (IDY-1)*KHRS+IHR
          TOTFCO2(I) = 0.0
          TOTFRF(I) = 0.0
          TOTFRG(I) = 0.0
          TOTFH2OT(I) = 0.0
          TOTFH2OE(I) = 0.0
          TOTFH2OCAN(I) = 0.0
          TOTGSCAN(I) = 0.0
          TOTFHT(I) = 0.0
	    AVTCAN(I) = 0.0
          DO 20 J = 1,3
            TOTHRAB(I,J) = 0.0
20        CONTINUE
10      CONTINUE

        DO 30 J = 1,3
          TOTDYAB(IDY,J) = 0.0
30      CONTINUE

        ICHECK(IDY) = 0
        TOTDCO2(IDY) = 0.0
        TOTDCO2NET(IDY) = 0.0
        TOTDRF(IDY) = 0.0
        TOTDH2O(IDY) = 0.0
        TOTDH2OCAN(IDY) = 0.0
        TOTDHT(IDY) = 0.0
100   CONTINUE

      DO 40 I = 1,MAXHISTO
        HISTO(I) = 0.0
40    CONTINUE

      RETURN
      END !ZeroT


C**********************************************************************
      SUBROUTINE OUTPUTTDAY(ICHECK,STOCKING,NTREES,
     &  TOTDYAB,TOTDCO2,TOTDRF,TOTDCO2NET,
     &  TOTDH2O,TOTDH2OCAN,TOTDHT)
C Output daily totals to file.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      INTEGER ICHECK(MAXDAY)
      REAL TOTDYAB(MAXDAY,3),TOTDCO2(MAXDAY),TOTDCO2NET(MAXDAY)
      REAL TOTDH2O(MAXDAY),TOTDH2OCAN(MAXDAY)
      REAL TOTDRF(MAXDAY),TOTDHT(MAXDAY)

      DO 10 ID = 1,MAXDAY
        IF (ICHECK(ID).GT.0) THEN
          WRITE (UTDAY,600) ID,TOTDYAB(ID,1)*STOCKING/NTREES,
     &      TOTDYAB(ID,2)*STOCKING/NTREES,TOTDYAB(ID,3)*STOCKING/NTREES,
     &      TOTDCO2(ID)*STOCKING/NTREES,TOTDRF(ID)*STOCKING/NTREES,
     &      TOTDCO2NET(ID)*STOCKING/NTREES,
     &      TOTDH2O(ID)*STOCKING/NTREES,TOTDH2OCAN(ID)/NTREES,
     &      TOTDHT(ID)*STOCKING/NTREES
        END IF
10    CONTINUE
600   FORMAT (I4,1X,10(F12.5,1X))

      RETURN
      END !OutputTday


C**********************************************************************
      SUBROUTINE CLOSEF(IOHRLY,IOHIST,IOCON)
C Close files.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      IF (IOCON.GT.0) CLOSE(UTCON)
      CLOSE(UTREES)
      CLOSE(UDAILY)
      CLOSE(UERROR)
      IF (IOHRLY.GT.0) CLOSE(UHRLY)
      IF (IOHIST.EQ.1) CLOSE(UHIST)
      CLOSE(UTDAY)
      IF (IOHRLY.GT.0) CLOSE(UTHR)
      IF (IOHIST.EQ.1) CLOSE(UTHIST)

      RETURN
      END !Closef


C**********************************************************************
      SUBROUTINE READHR(UFILE,ITARGETS,NOTARGETS,WEIGHTS,
     &    TOTHRAB,TOTFCO2,TOTFH2OT,TOTFH2OE,TOTFHT,TOTGSCAN,
     &    TOTFH2OCAN,TOTFRF,TOTFRG,AVTCAN)
C Read in hourly values from file and sum.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL TOTHRAB(KHRS*MAXDAY,3),TOTFCO2(KHRS*MAXDAY)
      REAL TOTFH2OT(KHRS*MAXDAY),TOTFH2OE(KHRS*MAXDAY)
	REAL TOTFH2OCAN(KHRS*MAXDAY)
      REAL TOTFRF(KHRS*MAXDAY),TOTFRG(KHRS*MAXDAY)
      REAL TOTFHT(KHRS*MAXDAY),TOTGSCAN(KHRS*MAXDAY)
	REAL AVTCAN(KHRS*MAXDAY)
      REAL VALUES(12)
      REAL WEIGHTS(MAXT)
      INTEGER ITARGETS(MAXT)
      INTEGER UFILE

C Skip comments at start of file
990   FORMAT (A80)
      REWIND (UFILE)
      DO 10 I = 1,24
        READ (UFILE,990)
10    CONTINUE

C Read hours. If tree is in list, add to sum. 
500   FORMAT (I7,1X,2(I4,1X),3(F12.5,1X),9(F12.7,1X))
      IPTREE = 0
20    READ (UFILE,500,IOSTAT=IOERROR) 
     &  IDAY,ITREE,IHOUR,(VALUES(I),I=1,12)
      IF (IOERROR.NE.0) RETURN
      IF (ITREE.NE.IPTREE) ITAR = INLIST(ITREE,ITARGETS,NOTARGETS)
      IPTREE = ITREE
      IF (ITAR.GT.0) THEN
        IHR = (IDAY-1)*KHRS + IHOUR
        DO 50 I = 1,3
          TOTHRAB(IHR,I) = TOTHRAB(IHR,I) + VALUES(I)*WEIGHTS(ITAR)
50      CONTINUE
        TOTFCO2(IHR) = TOTFCO2(IHR) + VALUES(4)*WEIGHTS(ITAR)
        TOTFRF(IHR) = TOTFRF(IHR) + VALUES(5)*WEIGHTS(ITAR)
        TOTFRG(IHR) = TOTFRG(IHR) + VALUES(6)*WEIGHTS(ITAR)
        TOTFH2OT(IHR) = TOTFH2OT(IHR) + VALUES(7)*WEIGHTS(ITAR)
        TOTFH2OE(IHR) = TOTFH2OE(IHR) + VALUES(8)*WEIGHTS(ITAR)
        TOTFH2OCAN(IHR) = TOTFH2OCAN(IHR) + VALUES(9)*WEIGHTS(ITAR)
        TOTGSCAN(IHR) = TOTGSCAN(IHR) + VALUES(10)*WEIGHTS(ITAR)
        TOTFHT(IHR) = TOTFHT(IHR) + VALUES(11)*WEIGHTS(ITAR)
	  AVTCAN(IHR) = AVTCAN(IHR) + VALUES(12)*WEIGHTS(ITAR)
      END IF
      GOTO 20

      RETURN
      END !ReadHr


C**********************************************************************
      SUBROUTINE OUTPUTTHR(ICHECK,STOCKING,NTREES,
     &    TOTHRAB,TOTFCO2,TOTFRF,TOTFRG,
     &    TOTFH2OT,TOTFH2OE,TOTFH2OCAN,TOTFHT,TOTGSCAN,AVTCAN)
C Output hourly totals to file.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL TOTHRAB(KHRS*MAXDAY,3),TOTFCO2(KHRS*MAXDAY)
      REAL TOTFH2OT(KHRS*MAXDAY),TOTFH2OE(KHRS*MAXDAY)
	REAL TOTFH2OCAN(KHRS*MAXDAY),AVTCAN(KHRS*MAXDAY)
      REAL TOTFRF(KHRS*MAXDAY),TOTFRG(KHRS*MAXDAY)
      REAL TOTFHT(KHRS*MAXDAY),TOTGSCAN(KHRS*MAXDAY)
      INTEGER ICHECK(MAXDAY)

      DO 10 IDAY = 1,MAXDAY
        IF (ICHECK(IDAY).GT.0) THEN
          DO 20 IHOUR = 1,KHRS
            IHR = (IDAY-1)*KHRS+IHOUR
            WRITE (UTHR,500) IDAY,IHOUR,
     &        TOTHRAB(IHR,1)*STOCKING/NTREES,
     &        TOTHRAB(IHR,2)*STOCKING/NTREES,
     &        TOTHRAB(IHR,3)*STOCKING/NTREES,
     &        TOTFCO2(IHR)*STOCKING/NTREES,
     &        TOTFRF(IHR)*STOCKING/NTREES,TOTFRG(IHR)*STOCKING/NTREES,
     &        TOTFH2OT(IHR)*STOCKING/NTREES,
     &        TOTFH2OE(IHR)*STOCKING/NTREES,
     &		TOTFH2OCAN(IHR)/NTREES,
     &        TOTGSCAN(IHR)*STOCKING/NTREES,
     &		TOTFHT(IHR)*STOCKING/NTREES*1000.0,
     &		AVTCAN(IHR)/NTREES
20        CONTINUE
        END IF
10    CONTINUE
500   FORMAT (2(I4,1X),3(F12.5,1X),9(F12.7,1X))

      RETURN
      END !OutputTHr


C**********************************************************************
      SUBROUTINE READHIST(UFILE,ITARGETS,NOTARGETS,WEIGHTS,
     &    PAR,HISTO)
C Read in and sum histograms.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      INTEGER ITARGETS(MAXT)
      REAL WEIGHTS(MAXT)
      REAL PAR(MAXHISTO),HISTO(MAXHISTO)
      CHARACTER*10 TMPSTR
      INTEGER UFILE

500   FORMAT (A9,I6)
520   FORMAT (F8.2,1X,F10.6)
990   FORMAT (A80)
20    READ (UFILE,500,IOSTAT = IOERROR) TMPSTR,ITREE
      IF (IOERROR.NE.0) RETURN
      ITAR = INLIST(ITREE,ITARGETS,NOTARGETS)
      READ (UFILE,990)
      DO 10 I = 1,MAXHISTO
        READ (UFILE,520) PAR(I),HISTMP
        IF (ITAR.GT.0) THEN
          HISTO(I) = HISTO(I) + HISTMP*WEIGHTS(ITAR)
        END IF
10    CONTINUE
      GOTO 20
      RETURN
      END !ReadHist

C**********************************************************************
      SUBROUTINE OUTPUTTHIST(UFILE,PAR,HISTO,STOCKING,NTREES)
C Write total PAR histogram to file.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL PAR(MAXHISTO),HISTO(MAXHISTO)
      INTEGER UFILE

      WRITE (UFILE,510) 
      DO 10 IBIN = 1,MAXHISTO
        WRITE (UFILE,520) PAR(IBIN),HISTO(IBIN)*STOCKING/NTREES
10    CONTINUE
510   FORMAT ('  PAR:         FREQUENCY (M^2.S): ')
520   FORMAT (F8.2,1X,F10.6)

      RETURN
      END !OutputTHist


