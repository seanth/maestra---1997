$DEBUG:'D' !DEBUG VERSION; ALL LINES BEGINNING WITH 'D' ARE COMPILED IN

C**********************************************************************
C GETMET.FOR
C
C This file contains all the subroutines required to read the met data
C file. The main routines (called externally to this file) are:
C OPENMETF - opens the met file & reads the format information
C RESTARTMETF - finds the start of the actual met data in the met file
C GETMET - calls the appropriate subroutine to read the met data, based
C          on the flag MFLAG describing the format
C SKIPMET - skips any unwanted met data
C ALTERMETOTC - alters met data to account for effects of OTC
C ALTERMETCC - alters met data in accordance with climate change scenario
C GETWIND - calculates wind speed through the canopy using exponential formula
C
C Subsidiary routines are:
C 1. Subroutines to read in data
C GETMETDAY - reads the standard daily met data format
C GETMETHR - reads the standard hourly met data format
C READENVIRON - reads in constant met data values
C READLAT - reads latitude and longitude
C READDELTAT - reads in information needed to calculate incident radiation, if needed
C 2. Subroutines to process temperatures & rainfall
C CALCTHRLY - calculates hourly air temps from Tmax & Tmin
C CALCTSOIL - sets hourly soil temp equal to average daily air temp
C ASSIGNRAIN - assigns daily rainfall to hours in the day
C 3. Subroutines to process VPD/RH
C CALCRH - calculates hourly RH from hourly air temp & the minimum temp
C RHTOVPD - converts RH to VPD
C VPDTORH - converts VPD to RH
C TDEWTORH - calculates RH from air temperature & dewpoint temperature
C MFDTORH - converts from water vapour mole fraction to RH
C VPDTOMFD - converts from VPD to water vapour mole fraction
C CALCAH - converts relative humidity to absolute humidity (called by ALTERMETCC)
C AHTORH - converts absolute humidity to relative humidity (called by ALTERMETCC)
C 4. Subroutines to process radiation
C ESTPARIN - estimates daily incident radiation from temperature
C CALPARHHRLY - calculates hourly incident PAR from daily totals of beam & diffuse
C CALCFBMD - calculates diffuse fraction of daily incident radiation
C CALCFBMH - calculates diffuse fraction of hourly incident radiation
C CALCFSUN - calculate sunlit fraction of hour (presently unused)
C CALCNIR - calculates incident NIR from incident PAR
C THERMAL - calculates incident thermal radiation from air temperature
C CALCFBMWN - an alternative method to calculate the diffuse fraction (currently unused)
C 5. Additionally, the following functions are provided:
C SATUR - calculate water vapour saturation pressure at given T
C TK - calculate T in Kelvin given T in Celsius
C ETRAD - calculate radiation incident on the atmosphere
C**********************************************************************

C**********************************************************************
      SUBROUTINE OPENMETF(ISTART,IEND,CAK,PRESSK,SWMIN,SWMAX,
     &  DIFSKY,ALAT,TTIMD,DELTAT,
     &  MFLAG,METCOLS,NOMETCOLS,MTITLE,MSTART)
C This subroutine opens the meteorological data file and reads the
C data describing it. It checks that it covers the dates of the simulation
C and returns the format of the met file.
C INPUTS:
C ISTART, IEND - Desired start and end of simulation (in days-since-1950) - from confile. If not defined,
C	program uses all met data. 
C OUTPUTS:
C CAK - Constant CO2 concentration (umol mol-1), if specified; 0 otherwise
C PRESSK - Constant air pressure (Pa), if specified; PATM otherwise
C SWMIN - Minimum soil water content
C SWMAX - Maximum soil water content - in same units as soil water content in met data
C DIFSKY - Parameter indicating light distribution over cloudy sky, if specified; default 0
C ALAT - latitude, in radians (-ve for S hemisphere)
C TTIMD - time difference between longitude of plot & longitude of time zone, in hours
C MFLAG - indicates whether daily (0) or hourly (1) met data
C METCOLS - array containing the number of the column corresponding to each met data element used by program
C NOMETCOLS - number of columns of met data
C MTITLE - title of met data file
C MSTART - start date of met data (in days-since-1950 format)
C DELTAT - monthly mean daily temperature amplitude
C**********************************************************************

      INCLUDE 'MAESTCOM'
      INCLUDE 'METCOM'

      
      CHARACTER*8 COLUMNS(MAXMET)
      CHARACTER*10 START,END
      CHARACTER*80 MTITLE
	REAL DELTAT(12)
      INTEGER METCOLS(MAXMET)
	NAMELIST /FORMAT/ DAYORHR, NOCOLUMNS, COLUMNS, START, END

990   FORMAT (A60)     ! For reading titles in input files.
991   FORMAT (A12,A60) ! For writing comments to output files.

C Open input file with met data
      OPEN (UMET, FILE = 'met.dat', STATUS='OLD')
      READ (UMET, 990) MTITLE

C Read in unchanging met data
      CALL READENVIRON(UMET, CAK, PRESSK, DIFSKY, SWMIN, SWMAX)

C Read in latitude & longitude
      CALL READLAT(UMET, ALAT, TTIMD)

C Read format
      REWIND(UMET)
      READ (UMET,FORMAT,IOSTAT = IOERROR)
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('ERROR IN MET FILE: MISSING FORMAT INFORMATION',
     &  IFATAL,IOERROR)

C Check dates of simulation are covered
      MSTART = IDATE50(START)
      MEND = IDATE50(END)
      IF ((ISTART.EQ.0).AND.(IEND.EQ.0)) THEN
        ISTART = MSTART
        IEND = MEND
      ELSE
        IF ((ISTART.LT.MSTART).OR.(IEND.GT.MEND))
     &    CALL SUBERROR('ERROR: MET FILE DOES NOT COVER DATES',
     &    IFATAL,IOERROR)
      END IF

C Process columns - assign column numbers to METCOLS array
      IF (NOCOLUMNS.GT.MAXMET) THEN
        CALL SUBERROR('ERROR: TOO MANY MET DATA COLUMNS: DATA LOST',
     &  IWARN,IERROR)
        NOCOLUMNS = MAXMET
      END IF

C Initialise array of met columns
      DO 10 ICOL = 1,MAXMET
        METCOLS(ICOL) = MISSING
10    CONTINUE

      MFLAG = DAYORHR

C If MFLAG = 0, look for daily values; if MFLAG = 1, hourly values
      IF (MFLAG.EQ.0) THEN ! Daily values
        DO 20 I = 1,NOCOLUMNS
          ICOL = MISSING
          IF (COLUMNS(I).EQ.'DATE')     THEN
            ICOL = MDDATE
          ELSEIF (COLUMNS(I).EQ.'WIND') THEN
            ICOL = MDWIND
          ELSEIF (COLUMNS(I).EQ.'TMIN') THEN
            ICOL = MDTMIN
          ELSEIF (COLUMNS(I).EQ.'TMAX') THEN
            ICOL = MDTMAX
          ELSEIF (COLUMNS(I).EQ.'PAR')  THEN
            ICOL = MDPAR
          ELSEIF (COLUMNS(I).EQ.'SI')  THEN
            ICOL = MDSI
          ELSEIF (COLUMNS(I).EQ.'FBEAM')  THEN
            ICOL = MDFBEAM
          ELSEIF (COLUMNS(I).EQ.'CA')  THEN
            ICOL = MDCA
          ELSEIF (COLUMNS(I).EQ.'PRESS')  THEN
            ICOL = MDPRESS
          ELSEIF (COLUMNS(I).EQ.'PPT')  THEN
            ICOL = MDPPT
          ELSEIF (COLUMNS(I).EQ.'SW')  THEN
            ICOL = MDSW
          ENDIF
          IF (ICOL.NE.MISSING) METCOLS(ICOL) = I
20      CONTINUE

C Check to see if any of the essential information is missing.
        IF (METCOLS(MDTMIN).EQ.MISSING)
     &    CALL SUBERROR('ERROR: NEED VALUES OF TMIN IN MET FILE',
     &    IFATAL,0)
        IF (METCOLS(MDTMAX).EQ.MISSING)
     &    CALL SUBERROR('ERROR: NEED VALUES OF TMAX IN MET FILE',
     &    IFATAL,0)
        IF ((METCOLS(MDPAR).EQ.MISSING).AND.(METCOLS(MDSI).EQ.MISSING))
     &    CALL READDELTAT(UMET,DELTAT,IOERROR)
        IF ((CAK.EQ.0).AND.(METCOLS(MDCA).EQ.MISSING))
     &    CALL SUBERROR('ERROR: NO VALUE FOR CA IN MET FILE',
     &    IFATAL,0)

      ELSE                 ! Hourly values
        DO 30 I = 1,NOCOLUMNS
          ICOL = MISSING
          IF (COLUMNS(I).EQ.'DATE')     THEN
            ICOL = MHDATE
          ELSEIF (COLUMNS(I).EQ.'WIND') THEN
            ICOL = MHWIND
          ELSEIF (COLUMNS(I).EQ.'TAIR') THEN
            ICOL = MHTAIR
          ELSEIF (COLUMNS(I).EQ.'TSOIL') THEN
            ICOL = MHTSOIL
          ELSEIF (COLUMNS(I).EQ.'RH')  THEN
            ICOL = MHRH
          ELSEIF (COLUMNS(I).EQ.'RH%')  THEN
            ICOL = MHRHP
          ELSEIF (COLUMNS(I).EQ.'VPD')  THEN
            ICOL = MHVPD
          ELSEIF (COLUMNS(I).EQ.'PAR')  THEN
            ICOL = MHPAR
          ELSEIF (COLUMNS(I).EQ.'RAD')  THEN
            ICOL = MHRAD
          ELSEIF (COLUMNS(I).EQ.'FBEAM')  THEN
            ICOL = MHFBEAM
          ELSEIF (COLUMNS(I).EQ.'CA')  THEN
            ICOL = MHCA
          ELSEIF (COLUMNS(I).EQ.'MFD')  THEN
            ICOL = MHMFD
          ELSEIF (COLUMNS(I).EQ.'PRESS')  THEN
            ICOL = MHPRESS
          ELSEIF (COLUMNS(I).EQ.'TDEW')  THEN
            ICOL = MHTDEW
          ELSEIF (COLUMNS(I).EQ.'PPT')  THEN
            ICOL = MHPPT
          ELSEIF (COLUMNS(I).EQ.'SW')  THEN
            ICOL = MHSW
          ENDIF
          IF (ICOL.NE.MISSING) METCOLS(ICOL) = I
30      CONTINUE

C Check to see if any of the essential information is missing.
        IF (METCOLS(MHTAIR).EQ.MISSING)
     &    CALL SUBERROR('ERROR: NEED VALUES OF TAIR IN MET FILE',
     &    IFATAL,0)
        IF ((METCOLS(MHPAR).EQ.MISSING).AND.(METCOLS(MHRAD).EQ.MISSING)) 
     &    CALL READDELTAT(UMET,DELTAT,IOERROR)
        IF ((CAK.EQ.0).AND.(METCOLS(MHCA).EQ.MISSING))
     &    CALL SUBERROR('ERROR: NO VALUE FOR CA IN MET FILE',
     &    IFATAL,0)

      END IF

      NOMETCOLS = NOCOLUMNS

      RETURN
      END !OpenMetF


C**********************************************************************
      SUBROUTINE RESTARTMETF(ISTART,MSTART,MFLAG)
C Position met file at start of met data
C BM Changed 10/00: Start of data must be directly specified
C INPUTS:
C ISTART - Date for which met data is wanted
C MSTART - Date on which met data starts 
C MFLAG - indicates whether daily or hourly data
C**********************************************************************

      INCLUDE 'MAESTCOM'

      CHARACTER*60 TMPTXT, DATASTART

      DATASTART = 'DATA STARTS'
      REWIND(UMET)

C Read through format information to find start of data
990   FORMAT (A60)
30    READ (UMET,990,IOSTAT = IOERROR) TMPTXT
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('ERROR: COULD NOT FIND START OF MET DATA',
     &  IFATAL,IOERROR)
	IF (TMPTXT.NE.DATASTART) GOTO 30

C Read data until the start date
      IF (MFLAG.EQ.0) THEN
        LINESTOSKIP = ISTART-MSTART
      ELSE
        LINESTOSKIP = KHRS*(ISTART-MSTART)
      END IF

      DO 40 I = 1,LINESTOSKIP
        READ (UMET,990,IOSTAT = IOERROR) TMPTXT
        IF (IOERROR.NE.0)
     &    CALL SUBERROR('READ PAST EOF IN RESTARTMETF',
     &    IFATAL,IOERROR)
40    CONTINUE

      RETURN
      END !RestartMetF


C**********************************************************************
      SUBROUTINE READENVIRON(UFILE,CAK,PRESSK,DIFSKYI,SWMINI,SWMAXI)
C Read in environmental conditions which do not change with time.
C These must be at the start of the met.dat file.
C INPUT:
C UFILE - file number of met data file
C OUTPUTS:
C CAK - Constant CO2 concentration (umol mol-1), if specified; 0 otherwise
C PRESSK - Constant air pressure (Pa), if specified; PATM otherwise
C DIFSKYI - Parameter indicating light distribution over cloudy sky, if specified; default 0
C SWMINI - Minimum soil water content
C SWMAXI - Maximum soil water content - in same units as soil water content in met data
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /ENVIRON/ CA,PRESS,DIFSKY, SWMIN, SWMAX
      INTEGER UFILE

C Default values
      PRESS = PATM
      CA = 0
      DIFSKY = 0.0
	SWMIN = 0.0
	SWMAX = 1.0

C Read namelist
      READ (UFILE, ENVIRON, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &  'WARNING: DEFAULT VALUES: PRESS = 101.25 kPa, DIFSKY = 0',
     &  IWARN,IOERROR)
      END IF
      CAK = CA
      PRESSK = PRESS
      DIFSKYI = DIFSKY
	SWMINI = SWMIN
	SWMAXI = SWMAX

      RETURN
      END !ReadEnviron

C**********************************************************************
      SUBROUTINE READDELTAT(UFILE,DELTATI,IOERROR)
C Read in mean monthly daily temperature amplitudes
C Needed to calculate incident PAR (if not specified) using routine of Bristow & Campbell 1984
C If not present, program will stop (only called if PAR not specified). 
C**********************************************************************

      INCLUDE 'MAESTCOM'
      NAMELIST /BRISTO/ DELTAT
	REAL DELTAT(12),DELTATI(12)
      INTEGER UFILE

C Read namelist
	REWIND(UFILE)
      READ (UFILE, BRISTO, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR(
     &    'ERROR: NEED VALUES OF PAR IN MET FILE OR DELTAT',
     &    IFATAL,0)
      END IF

	DO 10 I = 1,12
	  DELTATI(I) = DELTAT(I)
10	CONTINUE

      RETURN
      END !ReadDeltat

C**********************************************************************
      SUBROUTINE GETMET(IDATE,MFLAG,ZEN,METCOLS,NOMETCOLS,
     &  CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,DAYL,
     &  windah,tsoil,tair,radabv,fbeam,RH,VPD,VMFD,CA,PRESS,PPT,SOILMD)

C According to the value of MFLAG, this subroutine calls the appropriate
C subroutine to read in the meteorological data for use in MAESTRO.
C It is important that all handling of meteorological data occurs within
C this function. In particular, the flag MFLAG should not be referred to
C within the main body of the program. This maintains modularity & makes it
C easier to add new met data formats.

C The values passed to this function are:
C   IDATE - date of simulation, in years-since-1950 format
C   MFLAG - flag indicating the format of the input file
C   ZEN - the zenith angle of the sun
C   NOMETCOLS - the number of columns in the met file
C   METCOLS - the format of the columns in the met file
C   CAK - Constant CO2 concentration (umol mol-1), if specified; 0 otherwise
C   PRESSK - Constant air pressure (Pa), if specified; PATM otherwise
C   SWMINI - Minimum soil water content
C   SWMAXI - Maximum soil water content - in same units as soil water content in met data


C The values which must be returned are:
C arrays of KHRSly values of the following met variables:
C   WINDAH - wind speed in m s-1
C   TSOIL - soil surface temperature (degrees C)
C   TAIR - air temperature (degrees C)
C   RADABV(3) - incident radiation in three wavelengths (J m-2 s-1)
C   FBEAM(3) - beam fractions of radiation in 3 wavelengths (fraction)
C   RH - relative humidity of air (fraction)
C   VPD - vapour pressure deficit of air (Pa)
C   VMFD - mole fraction deficit (mmol mol-1)
C   CA - the atmospheric CO2 concentration (umol mol-1)
C   PRESS - the atmospheric pressure (Pa)
C   PPT - rainfall over the time period (mmol H2O m-2)
C   SOILMD - soil moisture deficit (dimensionless)
C**********************************************************************

      INCLUDE 'maestcom'

      REAL ZEN(KHRS)
      REAL windah(KHRS),tsoil(KHRS),tair(KHRS)
      REAL radabv(KHRS,3),fbeam(KHRS,3)
      REAL RH(KHRS),VPD(KHRS),CA(KHRS),VMFD(KHRS),PRESS(KHRS)
	REAL PPT(KHRS), SOILMD(KHRS)
	REAL DELTAT(12)
      INTEGER METCOLS(MAXMET)

C Call appropriate function to read met file.

      IF (MFLAG.EQ.0) THEN
        CALL GETMETDAY(IDATE,ZEN,NOMETCOLS,METCOLS,
     &  CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,DAYL,
     &  windah,tsoil,tair,radabv,fbeam,RH,VPD,VMFD,CA,PRESS,PPT,SOILMD)

      ELSE IF (MFLAG.EQ.1) THEN
        CALL GETMETHR(IDATE,ZEN,NOMETCOLS,METCOLS,
     &  CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,DAYL,
     &  windah,tsoil,tair,radabv,fbeam,RH,VPD,VMFD,CA,PRESS,PPT,SOILMD)

c      ELSE
C insert other formats if required
      END IF

      RETURN
      END !GetMet


C**********************************************************************
      subroutine GETMETDAY(IDATE,ZEN,NOMETCOLS,METCOLS,
     &  CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,DAYL,
     &  windah,tsoil,tair,radabv,fbeam,RH,VPD,VMFD,CA,PRESS,PPT,SOILMD)
C Read daily met data: see function GETMET for parameter definitions
C**********************************************************************

      INCLUDE 'maestcom'
      INCLUDE 'METCOM'

      REAL ZEN(KHRS)
      REAL DATAIN(MAXMET)
      REAL windah(KHRS),tsoil(KHRS),tair(KHRS),PPT(KHRS),SOILMD(KHRS)
      REAL radabv(KHRS,3),fbeam(KHRS,3),fsun(KHRS)
      REAL RH(KHRS),VPD(KHRS),VMFD(KHRS),CA(KHRS),PRESS(KHRS)
	REAL DELTAT(12)
      INTEGER METCOLS(MAXMET)

C Read in day's data
      READ (UMET,*,IOSTAT = IOERROR) (DATAIN(I), I = 1,NOMETCOLS)
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('ERROR READING MET DATA',IFATAL,IERROR)

C Set up wind speed array
      IF (METCOLS(MDWIND).EQ.MISSING) THEN
        WIND = DEFWIND  ! Very default value - CHECK!!
      ELSE IF (DATAIN(METCOLS(MDWIND)).LE.0) THEN
        WIND = DEFWIND
      ELSE
        WIND = DATAIN(METCOLS(MDWIND))
      END IF
      DO 30 I = 1,KHRS
        WINDAH(I) = WIND
30    CONTINUE

C Set up pressure array
      IF (METCOLS(MDPRESS).NE.MISSING) THEN
        DO 35 IHR = 1,KHRS
          PRESS(IHR) = DATAIN(METCOLS(MDPRESS))
35      CONTINUE
      ELSE
        DO 37 IHR = 1,KHRS
          PRESS(IHR) = PRESSK
37      CONTINUE
      END IF

C Calculate hourly temperatures
      TMAX = DATAIN(METCOLS(MDTMAX))
      TMIN = DATAIN(METCOLS(MDTMIN))
      CALL CALCTHRLY(TMAX,TMIN,DAYL,TAIR)

C Calculate soil temperatures
      CALL CALCTSOIL(TAIR,TSOIL)

C Calculate relative humidities
      CALL CALCRH(TMIN,TAIR,RH)

C Calculate VPDs
      CALL RHTOVPD(RH,TAIR,VPD)

C Calculate VMFDs
      CALL VPDTOMFD(VPD,PRESS,VMFD)

C Calculate incident PAR
      IF (METCOLS(MDSI).EQ.MISSING) THEN
     	  IF (METCOLS(MDPAR).EQ.MISSING) THEN
	    PRECIP = 0.0
          IF (METCOLS(MDPPT).NE.MISSING) PRECIP = DATAIN(METCOLS(MDPPT))
	    CALL BRISTO(IDAY,TMAX,TMIN,PRECIP,DELTAT,ALAT,DEC,DAYL,PAR)
	  ELSE
	    PAR = DATAIN(METCOLS(MDPAR))
	  END IF
	END IF

C Calculate diffuse fraction
      IF (METCOLS(MDFBEAM).EQ.MISSING) THEN
        IF (METCOLS(MDSI).EQ.MISSING) THEN
          CALL CALCFBMD(IDATE,ZEN,PAR,FBM)
        ELSE 
          CALL CALCFBMD(IDATE,ZEN,DATAIN(METCOLS(MDSI))*FPAR,FBM)
        END IF
      ELSE
        FBM = DATAIN(METCOLS(MDFBEAM))
      END IF

C Calculate PAR
      IF (METCOLS(MDSI).EQ.MISSING) THEN
        RADBM = PAR*FBM
        RADDF = PAR*(1.-FBM)
      ELSE
        RADBM = DATAIN(METCOLS(MDSI))*FBM*FPAR
        RADDF = DATAIN(METCOLS(MDSI))*(1.-FBM)*FPAR
      END IF
      CALL CALCPARHRLY(RADBM,RADDF,ZEN,RADABV,FBEAM)

C Calculate NIR
      CALL CALCNIR(RADABV,FBEAM)

C Calculate FSUN
      CALL CALCFSUN(FBEAM,FSUN)

C Calculate thermal radiation
      CALL THERMAL(TAIR,VPD,FSUN,RADABV)

C Read in value of CA if present - NB could add daily variation if wished?
      IF (METCOLS(MDCA).NE.MISSING) THEN
        DO 40 IHR = 1,KHRS
          CA(IHR) = DATAIN(METCOLS(MDCA))
40      CONTINUE
      ELSE
        DO 45 IHR = 1,KHRS
          CA(IHR) = CAK
45      CONTINUE
      END IF

C Read in rainfall data
	DO 50 IHR = 1,KHRS
	    PPT(IHR) = 0.0  ! Initially set all values to zero
50    CONTINUE
      IF (METCOLS(MDPPT).NE.MISSING) THEN
	  CALL ASSIGNRAIN(DATAIN(METCOLS(MDPPT)),PPT)
	END IF

C Calculate soil moisture deficit - assume unchanging over course of day
      IF (METCOLS(MDSW).NE.MISSING) THEN
        SMD = (SWMAX - DATAIN(METCOLS(MDSW)))/(SWMAX - SWMIN)
	ELSE
	  SMD = 0.0
	END IF
	DO 60 IHR = 1,KHRS
	  SOILMD(IHR) = SMD
60    CONTINUE

      RETURN
      END !GetMetDay


C**********************************************************************
      subroutine GETMETHR(IDATE,ZEN,NOMETCOLS,METCOLS,
     &  CAK,PRESSK,SWMIN,SWMAX,DELTAT,ALAT,DEC,DAYL,
     &  windah,tsoil,tair,radabv,fbeam,RH,VPD,VMFD,CA,PRESS,PPT,SOILMD)
C Read hourly met data: see function GETMET for parameter definitions
C**********************************************************************

      INCLUDE 'maestcom'
      INCLUDE 'METCOM'

      REAL DATAIN(KHRS,MAXMET)
      REAL windah(KHRS),tsoil(KHRS),tair(KHRS)
      REAL radabv(KHRS,3),fbeam(KHRS,3),fsun(KHRS)
      REAL RH(KHRS),VPD(KHRS),VMFD(KHRS),TDEW(KHRS)
      REAL CA(KHRS),PRESS(KHRS),ZEN(KHRS),PPT(KHRS),SOILMD(KHRS)
	REAL DELTAT(12)
      INTEGER METCOLS(MAXMET)

C Read in one day's worth of data at a time.
      DO 20 IHR = 1,KHRS
        READ (UMET,*,IOSTAT = IOERROR) (DATAIN(IHR,I), I = 1,NOMETCOLS)
        IF (IOERROR.NE.0)
     &    CALL SUBERROR('ERROR READING MET DATA',IFATAL,IERROR)
20    CONTINUE

C Set up pressure array
      IF (METCOLS(MHPRESS).NE.MISSING) THEN
        DO 35 IHR = 1,KHRS
          PRESS(IHR) = DATAIN(IHR,METCOLS(MHPRESS))
35      CONTINUE
      ELSE
        DO 37 IHR = 1,KHRS
          PRESS(IHR) = PRESSK
37      CONTINUE
      END IF

C Set up wind speed array
      DO 30 IHR = 1,KHRS
        IF (METCOLS(MHWIND).EQ.MISSING) THEN
          WINDAH(IHR) = DEFWIND     ! Very default value!!
        ELSE IF (DATAIN(IHR,METCOLS(MHWIND)).LE.0) THEN
          WINDAH(IHR) = DEFWIND
        ELSE
          WINDAH(IHR) = DATAIN(IHR,METCOLS(MHWIND))
        END IF
30    CONTINUE

C Set hourly air temperatures
      TMIN = DATAIN(1,METCOLS(MHTAIR))
      TMAX = DATAIN(1,METCOLS(MHTAIR))
      DO 40 IHR = 1,KHRS
        TAIR(IHR) = DATAIN(IHR,METCOLS(MHTAIR))
        IF (TAIR(IHR).LT.TMIN) TMIN = TAIR(IHR)
        IF (TAIR(IHR).GT.TMAX) TMAX = TAIR(IHR)
40    CONTINUE

C Set hourly soil temperatures
      IF (METCOLS(MHTSOIL).EQ.MISSING) THEN
        CALL CALCTSOIL(TAIR,TSOIL)
      ELSE
        DO 50 IHR = 1,KHRS
          TSOIL(IHR) = DATAIN(IHR,METCOLS(MHTSOIL))
50      CONTINUE
      END IF

C Read in RH, VPD, VMFD, Tdew
      IF (METCOLS(MHRH).NE.MISSING) THEN
        DO 60 IHR = 1,KHRS
          RH(IHR) = DATAIN(IHR,METCOLS(MHRH))
60      CONTINUE
      ELSE IF (METCOLS(MHRHP).NE.MISSING) THEN
        DO 70 IHR = 1,KHRS
          RH(IHR) = DATAIN(IHR,METCOLS(MHRHP))/100.0
70      CONTINUE
      END IF
      IF (METCOLS(MHVPD).NE.MISSING) THEN
        DO 80 IHR = 1,KHRS
          VPD(IHR) = DATAIN(IHR,METCOLS(MHVPD))
80      CONTINUE
      END IF
      IF (METCOLS(MHTDEW).NE.MISSING) THEN
        DO 82 IHR = 1,KHRS
          TDEW(IHR) = DATAIN(IHR,METCOLS(MHTDEW))
82      CONTINUE
      END IF
      IF (METCOLS(MHMFD).NE.MISSING) THEN
        DO 85 IHR = 1,KHRS
          VMFD(IHR) = DATAIN(IHR,METCOLS(MHMFD))
85      CONTINUE
      END IF

C Calculate RH if not in file
      IF ((METCOLS(MHRH).EQ.MISSING).
     &  AND.(METCOLS(MHRHP).EQ.MISSING)) THEN
        IF (METCOLS(MHVPD).NE.MISSING) THEN
          CALL VPDTORH(VPD,TAIR,RH)
        ELSEIF (METCOLS(MHTDEW).NE.MISSING) THEN
          CALL TDEWTORH(TDEW,TAIR,RH)
        ELSEIF (METCOLS(MHMFD).NE.MISSING) THEN
          CALL MFDTORH(VMFD,PRESS,TAIR,RH)
        ELSE
          CALL CALCRH(TMIN,TAIR,RH)
        END IF
      END IF
C Calculate VPD if not in file
      IF (METCOLS(MHVPD).EQ.MISSING) THEN
        CALL RHTOVPD(RH,TAIR,VPD)
      END IF
C Calculate VMFD if not in file
      IF (METCOLS(MHMFD).EQ.MISSING) THEN
        CALL VPDTOMFD(VPD,PRESS,VMFD)
      END IF

C Read in rainfall if present. Convert mm to mmol m-2.
	DAYPPT = 0.0
      IF (METCOLS(MHPPT).NE.MISSING) THEN
        DO 140 IHR = 1,KHRS
          PPT(IHR) = DATAIN(IHR,METCOLS(MHPPT))*1E6/18.
	    DAYPPT = DAYPPT+PPT(IHR)
140    CONTINUE
      ELSE
        DO 150 IHR = 1,KHRS
          PPT(IHR) = 0.0
150     CONTINUE
      END IF

C Must have either PAR (umol m-2 s-1) or global radiation (W m-2)
	IF (METCOLS(MHPAR).NE.MISSING) THEN
        DO 90 IHR = 1,KHRS
          RADABV(IHR,1) = DATAIN(IHR,METCOLS(MHPAR)) / UMOLPERJ
90      CONTINUE
	ELSE IF(METCOLS(MHSW).NE.MISSING) THEN
        DO 95 IHR = 1,KHRS
          RADABV(IHR,1) = DATAIN(IHR,METCOLS(MHRAD)) * FPAR
95      CONTINUE
	ELSE
	  CALL BRISTO(IDAY,TMAX,TMIN,DAYPPT,DELTAT,ALAT,DEC,DAYL,PAR)
        CALL CALCFBMD(IDATE,ZEN,PAR,FBM)
        RADBM = PAR*FBM
        RADDF = PAR*(1.-FBM)
        CALL CALCPARHRLY(RADBM,RADDF,ZEN,RADABV,FBEAM)
	END IF

C Calculate beam fractions
      IF (METCOLS(MHFBEAM).NE.MISSING) THEN
        DO 100 IHR = 1,KHRS
          FBEAM(IHR,1) = DATAIN(IHR,METCOLS(MHFBEAM))
100     CONTINUE
      ELSE
        DO 110 IHR = 1,KHRS
          FBEAM(IHR,1) = CALCFBMH(IDATE,ZEN(IHR),RADABV(IHR,1))
C 29/3     FBEAM(IHR,1) = CALCFBMH(IDATE,IHR,ALAT,DEC,RADABV(IHR,1))
C          FBEAM(IHR,1) = CALCFBMWN(IDATE,IHR,ZEN(IHR),RADABV(IHR,1))
110     CONTINUE
      END IF

C Calculate NIR
      CALL CALCNIR(RADABV,FBEAM)

C Calculate FSUN
      CALL CALCFSUN(FBEAM,FSUN)

C Calculate thermal radiation
      CALL THERMAL(TAIR,VPD,FSUN,RADABV)

C Read in values of CA if present
      IF (METCOLS(MHCA).NE.MISSING) THEN
        DO 120 IHR = 1,KHRS
          CA(IHR) = DATAIN(IHR,METCOLS(MHCA))
120    CONTINUE
      ELSE
        DO 130 IHR = 1,KHRS
          CA(IHR) = CAK
130     CONTINUE
      END IF

C Calculate soil moisture deficit, if possible.
      IF (METCOLS(MHSW).NE.MISSING) THEN
	  DO 160 IHR = 1,KHRS
	    SOILMD(IHR) = (SWMAX - DATAIN(IHR,METCOLS(MHSW)))
     &    /(SWMAX - SWMIN)
160     CONTINUE
      ELSE
        DO 170 IHR = 1,KHRS
          SOILMD(IHR) = 0.0
170     CONTINUE      
      END IF

      RETURN
      END !GetMetHr


C**********************************************************************
      SUBROUTINE CALCTHRLY(TMAX,TMIN,DAYL,TAIR)
C Calculate a daily variation in temperature from max & min temperatures.
C Temp varies linearly between sunset & sunrise, and sinusoidally during the day.
C INPUTS:
C TMIN, TMAX - minimum and maximum daily temperatures, °C
C DAYL - daylength, hours
C OUTPUTS:
C TAIR - array of hourly air temperatures, °C
C**********************************************************************

      INCLUDE 'maestcom'

      REAL TAIR(KHRS)

      TAV = (TMAX+TMIN)/2.0
      TAMPL = (TMAX-TMIN)/2.0

      DO 105 I=1,KHRS
        HRTIME = I-0.5
        TIME = HRTIME + DAYL*0.5 - (KHRS/2.0)
        IF (TIME.LT.0.0.OR.TIME.GT.DAYL) THEN
          IF (TIME.LT.0.0) HRTIME=HRTIME+KHRS
          TAIR(I) = TAV - (TAV-TMIN)*(HRTIME-DAYL*0.5-(KHRS/2.0))/
     &      (KHRS-DAYL)
        ELSE
          TAIR(I) = TAV - TAMPL*COS(1.5*PI*TIME/DAYL)
        ENDIF

105   CONTINUE

      RETURN
      END !CalcTHrly


C**********************************************************************
      SUBROUTINE CALCTSOIL(TAIR,TSOIL)
C Calculate soil temperatures.
C Set equal to average daily air temperature.
C INPUTS:
C TAIR - array of hourly air temperatures, °C
C OUTPUTS:
C TSOIL - array of hourly soil temperatures, °C
C**********************************************************************

      INCLUDE 'maestcom'

      REAL TSOIL(KHRS),TAIR(KHRS)

      TAV = 0.0
      DO 10 I = 1,KHRS
        TAV = TAV + TAIR(I)
10    CONTINUE
      TAV = TAV/REAL(KHRS)

      DO 20 I = 1,KHRS
        TSOIL(I) = TAV
20    CONTINUE

      RETURN
      END !CalcTSoil


C**********************************************************************
	SUBROUTINE ASSIGNRAIN(TOTAL,PPT)
C Given daily total of PPT, assign to hours through the day.
C Use algorithm from GRAECO (model of D. Loustau). 
C**********************************************************************

	USE MSFLIB
      INCLUDE 'maestcom'
      REAL(4)  ran
      REAL PPT(KHRS) ! Hourly rainfall: assume already set to zero


      CALL SEED(1995) ! For testing - same random numbers - change to 0 later. 

C 1. All rain falls in one hour for light storms (<2 mm)
      IF (TOTAL.LE.2.) THEN
      CALL RANDOM(ran)   ! randomly select hour
	  IRAIN = INT(ran*KHRS)+1
	  PPT(IRAIN) = TOTAL

C 2. All rain falls in 24 hours for storms >46 mm 
      ELSE IF (TOTAL.GT.46.) THEN
	  RAIN = TOTAL/REAL(KHRS)
	  DO 10 IHR = 1,KHRS
		PPT(IHR) = RAIN
10	  CONTINUE

C 3. All rain falls at 2mm/hour at a random time of the day */
      ELSE
	  IHRSWITHRAIN = MIN(INT((TOTAL/2)*KHRS/24),KHRS)
	  RATE = TOTAL/REAL(IHRSWITHRAIN)
	  DO 20 I = 1, IHRSWITHRAIN
        CALL RANDOM(ran)   ! randomly select hours
	    IRAIN = INT(ran*KHRS)+1
	    PPT(IRAIN) = PPT(IRAIN)+RATE
20	  CONTINUE
      END IF

C Convert from mm to mmol m-2
      DO 30 IHR = 1, KHRS
	  PPT(IHR) = PPT(IHR)*1E6/18.
30    CONTINUE

      RETURN
	END !AssignRain


C**********************************************************************
      SUBROUTINE CALCRH(TMIN,TAIR,RH)
C Calculate hourly relative humidity from air temperature and minimum
C daily temperature (assumed to be the dewpoint).
C INPUTS:
C TMIN - daily minimum temperature, °C
C TAIR - array of hourly air temperatures, °C
C OUTPUTS:
C RH - array of hourly relative humidity, fraction
C**********************************************************************

      INCLUDE 'maestcom'

      REAL TAIR(KHRS),RH(KHRS)

      HUMABS = SATUR(TMIN)
      DO 10 I = 1,KHRS
        RH(I) = HUMABS/SATUR(TAIR(I))
10    CONTINUE

      RETURN
      END !CalcRH


C**********************************************************************
      FUNCTION SATUR(TAC)
C Calculate saturated water vapour pressure (Pa) at temperature TAC (Celsius)
C from Jones 1992 p 110 (note error in a - wrong units)
C**********************************************************************

      SATUR = 613.75*EXP(17.502*TAC/(240.97+TAC))

      RETURN
      END !Satur


C**********************************************************************
      SUBROUTINE RHTOVPD(RH,TAIR,VPD)
C Convert from RH to VPD.
C Inputs: RH(KHRS) - hourly relative humidity, fraction
C   TAIR(KHRS) - hourly air temperature, degrees C
C Outputs: VPD(KHRS) - hourly VPD, Pa
C**********************************************************************

      INCLUDE 'maestcom'

      REAL TAIR(KHRS),RH(KHRS),VPD(KHRS)

      DO 10 I = 1,KHRS
         VPD(I)=(1.0-RH(I))*SATUR(TAIR(I))
10    CONTINUE

      RETURN
      END !RHToVPD

C**********************************************************************
      SUBROUTINE VPDTORH(VPD,TAIR,RH)
C Convert from VPD to RH.
C Inputs:
C   VPD(KHRS) - hourly VPD, Pa
C   TAIR(KHRS) - hourly air temperature, degrees C
C Outputs:
C   RH(KHRS) - hourly relative humidity, fraction
C**********************************************************************

      INCLUDE 'maestcom'

      REAL TAIR(KHRS),RH(KHRS),VPD(KHRS)

      DO 10 I = 1,KHRS
         RH(I) = 1.0 - VPD(I)/SATUR(TAIR(I))
10    CONTINUE

      RETURN
      END !VPDToRH

C**********************************************************************
      SUBROUTINE TDEWTORH(TDEW,TAIR,RH)
C Convert from TDEW to RH.
C Inputs:
C   TDEW(KHRS) - hourly dewpoint temperature, degrees C
C   TAIR(KHRS) - hourly air temperature, degrees C
C Outputs:
C   RH(KHRS) - hourly relative humidity, fraction
C**********************************************************************

      INCLUDE 'maestcom'

      REAL TAIR(KHRS),RH(KHRS),TDEW(KHRS)

      DO 10 I = 1,KHRS
         RH(I) = SATUR(TDEW(I))/SATUR(TAIR(I))
10    CONTINUE

      RETURN
      END !TDEWToRH


C**********************************************************************
      SUBROUTINE MFDTORH(VMFD,PRESS,TAIR,RH)
C Convert from mole fraction deficit to relative humidity.
C Inputs:
C   VMFD(KHRS) - hourly mole fraction deficit, mmol mol-1
C   TAIR(KHRS) - hourly air temperature, degrees C
C   PRESS(KHRS) - hourly pressure, Pa
C Outputs:
C   RH(KHRS) - hourly relative humidity, fraction
C**********************************************************************

      INCLUDE 'maestcom'

      REAL VMFD(KHRS),PRESS(KHRS),RH(KHRS),TAIR(KHRS)

      DO 10 I = 1,KHRS
        RH(I) = 1 - VMFD(I)*PRESS(I)*1E-3/SATUR(TAIR(I))
10    CONTINUE

      RETURN
      END !MFDToRH


C**********************************************************************
      SUBROUTINE VPDTOMFD(VPD,PRESS,VMFD)
C Convert from VPD to VMFD.
C Inputs:
C   VPD(KHRS) - hourly VPD, kPa
C   PRESS(KHRS) - hourly pressure, Pa
C Outputs:
C   VMFD(KHRS) - hourly mole fraction deficit, mmol mol-1
C**********************************************************************

      INCLUDE 'maestcom'

      REAL VPD(KHRS),PRESS(KHRS),VMFD(KHRS)

      DO 10 I = 1,KHRS
         VMFD(I)= VPD(I)/PRESS(I)*1E3
10    CONTINUE

      RETURN
      END ! VPDToMFD


C**********************************************************************
	SUBROUTINE BRISTO(IDAY,TMAX,TMIN,PRECIP,DELTAT,ALAT,DEC,DAYL,PAR)
C Subroutine implements Bristow & Campbell (1984) Ag For Met 31:159-166
C Calculates incident daily radiation from temperature data
C**********************************************************************

	INCLUDE 'MAESTCOM'
	REAL DELTAT(12)

	DELT = TMAX - TMIN					!Daily temperature amplitude - should include prev. day but ..
	IF (PRECIP.GT.10.) DELT = DELT*0.75	!Reduce if rainy
	IMON = MONTH(IDAY)					!Month
	BRISTOK = 0.036*EXP(-0.154*DELTAT(IMON))	
	TRANSM = 1 - EXP(-BRISTOK*(DELT)**2.4)	!Transmittance

	ANGDL=DAYL*PI/KHRS					! Angular daylength
      RADCLR=0.0864*SOLARC/PI*TAU* (ANGDL*SIN(ALAT)*SIN(DEC)
     &   + COS(ALAT)*COS(DEC)*SIN(ANGDL))	!Incident radiation on clear day
	RADTOT = RADCLR*TRANSM				!Estimated solar radiation
	PAR = RADTOT*FPAR					!Estimated PAR

	RETURN
	END !BristoPAR


C**********************************************************************
      SUBROUTINE CALCPARHRLY(RADBM,RADDF,ZEN,RADABV,FBEAM)
C Calculate daily course of incident PAR from daily totals of
C beam and diffuse PAR (in MJ m-2 d-1).
C INPUTS:
C RADBM - daily total beam PAR (MJ m-2 d-1)
C RADDF - daily total diffuse PAR (MJ m-2 d-1)
C ZEN - zenith angle (radians)
C OUTPUTS:
C RADABV - array of hourly incident PAR (J m-2 s-1)
C FBEAM - array of hourly beam fractions
C**********************************************************************

      INCLUDE 'maestcom'

      REAL ZEN(KHRS),RADABV(KHRS,3),FBEAM(KHRS,3)
      REAL COSBM(KHRS),COSDF(KHRS)

      SUMBM = 0.0
      SUMDF = 0.0

      DO 105 I=1,KHRS
        COSBM(I) = 0.0
        COSDF(I) = 0.0
        HRTIME = I-0.5
        COSZEN = COS(ZEN(I))
        IF (COSZEN.GT.0.0) THEN
          IF (ZEN(I).LT.80*PID180) THEN  !Set FBM = 0.0 for ZEN > 80 degrees
            COSBM(I)=COSZEN*TAU**(1.0/COSZEN)
          ELSE
            COSBM(I)=0.0
          END IF
          COSDF(I)=COSZEN
          SUMBM=SUMBM+COSBM(I)
          SUMDF=SUMDF+COSDF(I)
        ENDIF
105   CONTINUE

      DO 205 I=1,KHRS
        IF (SUMBM.GT.0.0) THEN
          RDBM = RADBM*COSBM(I)/SUMBM
        ELSE
          RDBM = 0.0
        END IF
        IF (SUMDF.GT.0.0) THEN
          RDDF = RADDF*COSDF(I)/SUMDF
        ELSE
          RDDF = 0.0
        END IF
C Convert from MJ hr-1 to J s-1
        RADABV(I,1)=(RDDF+RDBM)*1e6/SPERHR
        IF ((RDBM+RDDF).GT.0.0) THEN
          FBEAM(I,1) = RDBM/(RDBM+RDDF)
        ELSE
          FBEAM(I,1)= 0.00
        END IF
205   CONTINUE

      RETURN
      END !CalcPARHrly


C**********************************************************************
      SUBROUTINE CALCFBMD(IDATE,ZEN,PAR,FBM)
C Calculate the beam fraction from the total daily incident radiation.
C Use the formula of Spitters et al. (1986) Agric For Met 38:217-229.
C INPUTS:
C IDATE - date in days-since-1950 format
C ZEN - array of hourly sun zenith angle (radians)
C PAR - daily total incident PAR (MJ m-2 d-1)
C OUTPUTS:
C FBM - daily beam fraction
C**********************************************************************

      INCLUDE 'MAESTCOM'
	REAL ZEN(KHRS)

C Calculate extra-terrestrial radiation
      S0 = 0.0
      DO 10 IHR = 1,KHRS
	  SINB = SIN(PID2-ZEN(IHR))
        S0 = S0 + ETRAD(IDATE,SINB)*SPERHR/1E6
10    CONTINUE

C Spitter's formula
      TRANS = (PAR/FPAR) / S0
      IF (TRANS.LT.0.07) THEN
        FDIF = 1.
      ELSE IF (TRANS.LT.0.35) THEN
        FDIF = 1. - 2.3*(TRANS-0.07)**2
      ELSE IF (TRANS.LT.0.75) THEN
        FDIF = 1.33 - 1.46*TRANS
      ELSE
        FDIF = 0.23
      END IF
      FBM = 1. - FDIF

      RETURN
      END !CalcFBMD

C**********************************************************************
      FUNCTION CALCFBMH(IDATE,ZEN,RADABV)
C Calculate the beam fraction from the hourly incident radiation.
C Use the formula of Spitters et al. (1986) Agric For Met 38:217-229.
C INPUTS:
C IDATE - date in days-since-1950 format
C ZEN - hourly sun zenith angle (radians)
C RADABV - total incident PAR (J m-2 s-1)
C RETURNS:
C CALCFBMH - hourly beam fraction
C**********************************************************************

      INCLUDE 'MAESTCOM'

C SINB is the sine of the elevation of the sun above the horizon
C (= SIN(90-ZEN)). For zenith angles > 80 degrees will put FBEAM = 0.
	SINB=SIN(PID2-ZEN)

C Calculate extra-terrestrial radiation
      S0 = ETRAD(IDATE,SINB)

      IF (SINB.GT.0.17) THEN
C Spitter's formula
        TRANS = (RADABV/FPAR) / S0
        SPITR = 0.847 - 1.61*SINB + 1.04*SINB**2
        SPITK = (1.47 - SPITR)/1.66
        IF (TRANS.LE.0.22) THEN
          FDIF = 1.
        ELSE IF (TRANS.LE.0.35) THEN
          FDIF = 1. - 6.4*(TRANS-0.22)**2
        ELSE IF (TRANS.LE.SPITK) THEN
          FDIF = 1.47 - 1.66*TRANS
        ELSE
          FDIF = SPITR
        END IF
        CALCFBMH = 1. - FDIF
      ELSE
        CALCFBMH = 0.0
      END IF
      RETURN
      END !CalcFBMH


C**********************************************************************
      SUBROUTINE CALCNIR(RADABV,FBEAM)
C From incident PAR, calculate incident NIR.
C Parameters:
C RADABV - hourly incident radiation in 3 wavelengths (J m-2 s-1)
C FBEAM - beam fraction of incident radiation for 3 wavelengths
C**********************************************************************

      INCLUDE 'maestcom'

      REAL RADABV(KHRS,3),FBEAM(KHRS,3)

C Fbeam is assumed the same for NIR as PAR. Amount of NIR is calculated from
C amount of PAR & the fraction of total solar that is PAR (FPAR).
      DO 10 I = 1,KHRS
        RADABV(I,2) = RADABV(I,1)*(1./FPAR - 1.0)
        FBEAM(I,2) = FBEAM(I,1)
10    CONTINUE

      RETURN
      END !CalcNIR


C**********************************************************************
      SUBROUTINE THERMAL(TAIR,VPD,FSUN,RADABV)
C Calculate incident thermal radiation, if it has not been measured. 
C Several different formulae exist:
C Ying-Ping had the following formula, from MH Unsworth & JL Monteith
C (1975) Quart. J. R. Met. Soc. 101. pp13-24, pp25-34
C        RADABV(I,3) = 1.06*SIGMA* (TK(TAIR(I))**4) - 119.
C I replaced this with a formula from Brutsaert (1975) Water Resources
C Research 11: 742-744. Leuning et al (1995) PC&E 18: 1183-1200 say that
C Hatfield et al (1983) Water Resources Research 19: 285-288 tested
C a range of alternative formulae & found Brutsaert's to be the best.
C BM 22/5/98.
C This formula is appropriate for clear skies only. The thermal radiation
C from a cloudy sky will be larger since clouds radiate more effectively. 
C Changed to the formula taking this into account, given by Monteith & Unsworth
C 1990 Principles of Environmental Physics 2nd ed p53. 
C**********************************************************************

      INCLUDE 'maestcom'

      REAL RADABV(KHRS,3),TAIR(KHRS),VPD(KHRS),FSUN(KHRS)

      DO 10 I = 1,KHRS
C        EA = SATUR(TAIR(I)) - VPD(I) !Old formula - see comments
C        EMSKY = 0.642*(EA/TK(TAIR(I)))**(1./7.)
C        RADABV(I,3) = EMSKY*SIGMA*(TK(TAIR(I))**4)
	  SIGT4 = SIGMA*(TK(TAIR(I))**4)
	  EMCLEAR = 1.06*SIGT4 - 119.
	  RADABV(I,3) = FSUN(I)*EMCLEAR + (1.-FSUN(I))*(0.84+0.16*EMCLEAR)
10    CONTINUE

      RETURN
      END !Thermal


C**********************************************************************
      SUBROUTINE CALCFSUN(FBEAM,FSUN)
C Calculates the sunlit time from the beam fraction.
C**********************************************************************

      INCLUDE 'maestcom'

      REAL FBEAM(KHRS,3),FSUN(KHRS)

      DO 10 I = 1,KHRS
        FSUN(I) = 1.16*FBEAM(I,1)-0.01
        IF (FSUN(I).LT.0.00) FSUN(I)=0.00
        IF (FSUN(I).GT.1.00) FSUN(I)=1.00
10    CONTINUE

      RETURN
      END !CalcFsun


C**********************************************************************
      FUNCTION ETRAD(IDATE,SINB)
C Calculate the radiation incident on the atmosphere.
C Using formulae from Spitters et al (1986) Agric For Met 38:217
C Returns value in J m-2 s-1.
C**********************************************************************

      INCLUDE 'MAESTCOM'

      JDAY = JDATE(IDATE)

C Spitters' formula
      IF (SINB.GT.0.0) THEN
        ETRAD = SOLARC * (1 + 0.033*COS(REAL(JDAY)/365.0*TWOPI)) * SINB
      ELSE
        ETRAD = 0.0
      END IF

      RETURN
      END !ETRad


C**********************************************************************
      FUNCTION TK(TCELSIUS)
C Converts Celsius temperature to Kelvin.
C**********************************************************************

      INCLUDE 'maestcom'

      TK = TCELSIUS - ABSZERO
      RETURN
      END !TCelsius


C**********************************************************************
      SUBROUTINE SKIPMET(MFLAG,NSTEP)
C Read data from the met file until in the right position
C**********************************************************************

      INCLUDE 'MAESTCOM'

      IF (MFLAG.EQ.0) THEN
        LINESTOSKIP = NSTEP - 1
      ELSE
        LINESTOSKIP = KHRS*(NSTEP - 1)
      END IF

      DO 40 I = 1,LINESTOSKIP
        READ (UMET, 990) MTITLE
40    CONTINUE
990   FORMAT (A80)

      RETURN
      END !SkipMet


C**********************************************************************
      SUBROUTINE READLAT(UFILE, ALAT, TTIMD)
C Read in Lat & Long details - from met file:
C ALAT - latitude, in radians (-ve for S hemisphere)
C TTIMD - time difference between longitude of plot & longitude of time
C   zone, in hours
C**********************************************************************

      INCLUDE 'MAESTCOM'
      
      INTEGER UFILE
      REAL LAT(3),LONG(3)
      CHARACTER*1 LATHEM,LONHEM
	NAMELIST /LATLONG/ LATHEM,LAT,LONHEM,LONG,TZLONG

      REWIND (UFILE)
      READ (UFILE, LATLONG, IOSTAT = IOERROR)
      IF (IOERROR.NE.0)
     &  CALL SUBERROR('ERROR READING LATITUDE: IT MUST BE IN MET.DAT',
     &  IFATAL,IOERROR)
      ALAT = (LAT(1)+LAT(2)/60.0+LAT(3)/3600.0)*PID180
      IF (LATHEM.EQ.'S') ALAT = -ALAT
      ALONG = LONG(1) + LONG(2)/60.0 + LONG(3)/3600.0
      IF (LONHEM.EQ.'E') THEN
        ALONG = 360.0 - ALONG
        TZLONG = 360.0 - TZLONG
      END IF
      ALONG = ALONG*PID180
      TZLONG = TZLONG*PID180
      TIMDIF = KHRS/2.0*ALONG/PI
      TTIMD = KHRS/2.0/PI*(ALONG - TZLONG)

      RETURN
      END !ReadLat


C**********************************************************************
      SUBROUTINE GETWIND(FOLLAY,FOLntr,TOTLAI,EXTWIND,WINDLAY)
C Calculate decrease in wind speed with increasing canopy depth
C Uses exponential decline. Does not take into account different sizes of trees.
C INPUTS:
C FOLLAY - amount of foliage in each layer of target tree
C FOLNTR - total foliage of target tree
C TOTLAI - total LAI of forest
C EXTWIND - extinction coefficient of wind
C OUTPUTS:
C WINDLAY - windspeed in each layer of the canopy (m s-1)
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL FOLLAY(MAXLAY),FOLABV(MAXLAY),WINDLAY(MAXLAY)

C Calculate approximately LAI above each point (for wind speed)
      FOLABV(1) = FOLLAY(1)/2.
      DO 10 ILAY = 2,MAXLAY
        FOLABV(ILAY) = FOLABV(ILAY-1) 
     &      + (FOLLAY(ILAY-1) + FOLLAY(ILAY))/2.
10    CONTINUE

C Calculate decrease in wind speed with canopy depth
      DO 20 ILAY = 1,MAXLAY
        IF (FOLNTR.GT.0.0) THEN
          FOLABV(ILAY) = FOLABV(ILAY)/FOLntr * TOTLAI
        ELSE
          FOLABV(ILAY) = 0.0
        END IF
        WINDLAY(ILAY) = EXP(-EXTWIND*FOLABV(ILAY))
20    CONTINUE

      RETURN
      END !GetWind


C**********************************************************************
      SUBROUTINE ALTERMETCC(CA,TAIR,TSOIL,RH,VPD,VMFD,PRESS,
     &  CO2INC,TINC)
C Subroutine to change met data according to climate change scenario.
C Currently can increase CO2 by fixed amount (CO2INC, ppm)
C and/or T by a fixed amount (TINC, deg C).
C The absolute humidity (in g m-3) is maintained constant.
C VPD, RH and VMFD are re-calculated if temperature is changed.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL CA(KHRS),TAIR(KHRS),TSOIL(KHRS)
      REAL VPD(KHRS),RH(KHRS),VMFD(KHRS),PRESS(KHRS)

      DO 10 IHR = 1,KHRS
        CA(IHR) = CA(IHR) + CO2INC
        ABSHUM = CALCAH(RH(IHR),TAIR(IHR))
        TAIR(IHR) = TAIR(IHR) + TINC
        TSOIL(IHR) = TSOIL(IHR) + TINC
        RH(IHR) = AHTORH(ABSHUM,TAIR(IHR))
10    CONTINUE
      IF (TINC.NE.0.0) THEN
        CALL RHTOVPD(RH,TAIR,VPD)
        CALL VPDTOMFD(VPD,PRESS,VMFD)
      END IF

      RETURN
      END !AlterMetCC


C**********************************************************************
      FUNCTION CALCAH(RH,TAIR)
C Calculate absolute humidity (g m-3) from relative humidity & T.
C Conversion from Jones (1992) pp. 109-110.
C**********************************************************************

      CALCAH = RH*SATUR(TAIR)*2.17/TK(TAIR)

      RETURN
      END  ! CalcAH


C**********************************************************************
      FUNCTION AHTORH(ABSHUM,TAIR)
C Calculate relative humidity from absolute humidity & T.
C Conversion from Jones (1992) pp. 109-110.
C**********************************************************************

      AHTORH = ABSHUM*TK(TAIR)/2.17/SATUR(TAIR)

      RETURN
      END  ! CalcAH


C**********************************************************************
      SUBROUTINE ALTERMETOTC(TOTC,WINDOTC,PAROTC,FBEAMOTC,
     &  TAIR,TSOIL,WINDAH,RADABV,FBEAM,RH,VPD,VMFD,PRESS)
C Subroutine to change met data according to effects of OTC.
C Currently can increase T by fixed amount (deg C),
C set wind speed to constant amount (m s-1),
C reduce radiation by a constant fraction,
C and/or alter the beam fraction of radiation.
C The absolute humidity (in g m-3) is maintained constant.
C VPD, RH and VMFD are re-calculated if temperature is changed.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL windah(KHRS),tsoil(KHRS),tair(KHRS)
      REAL radabv(KHRS,3),fbeam(KHRS,3)
      REAL VPD(KHRS),RH(KHRS),VMFD(KHRS),PRESS(KHRS)

      DO 10 IHR = 1,KHRS
        ABSHUM = CALCAH(RH(IHR),TAIR(IHR))
        TAIR(IHR) = TAIR(IHR) + TOTC
        TSOIL(IHR) = TSOIL(IHR) + TOTC
        RH(IHR) = AHTORH(ABSHUM,TAIR(IHR))
        IF (WINDOTC.GT.0.0) WINDAH(IHR) = WINDOTC
        RADABV(IHR,1) = RADABV(IHR,1)*PAROTC
        FBEAM(IHR,1) = FBEAM(IHR,1)*FBEAMOTC
10    CONTINUE
      IF (TOTC.NE.0.0) THEN
        CALL RHTOVPD(RH,TAIR,VPD)
        CALL VPDTOMFD(VPD,PRESS,VMFD)
      END IF

      RETURN
      END !AlterMetOTC


C**********************************************************************
      FUNCTION CALCFBMWN(IDATE,IHR,ZEN,RADABV)
C Calculate the beam fraction of PAR using the formula of Weiss &
C Norman (1985) Agric For Met 34:205-214.
C An alternative expression to that of Spitters et al (1986)
C (see CALCFBMD and CALCFBMH).
C Assume that pressure = sea level pressure (101.325kPa)
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL FBEAM(KHRS,3)

      PARMA = 0.90
      PARMB = 0.70
      PARMC = 0.88
      PARMD = 0.68

c      DO 10 IWAVE = 1,3
c        FBEAM(IHR,IWAVE) = 0.00
c10    CONTINUE
      CALCFBMWN = 0.0

      IF (NIGHT(ZEN,RADABV).EQ.1) GO TO 100

      COSZEN = COS(ZEN)
      AIRMASS = 1./COSZEN

c Potential direct-beam PAR on a horizontal surface can be approx. by
      RBV = 600.*EXP(-0.185*AIRMASS)*COSZEN
C Potential visible diffuse radiation is approx. by
C BK 2/95 Changed after advice by F. Villalobos, Cordoba
      RDV = 0.4*(600.*COSZEN - RBV)

      AXLOG = ALOG10(AIRMASS)
      ADUM = (-1.195+.4495*AXLOG-.0345*AXLOG*AXLOG)
      WATABS = 1320.*(10.**(ADUM))
C Potential direct-beam NIR
      RBN = (720.*EXP(-0.06*AIRMASS)-WATABS)*COSZEN
C Potential diffuse NIR
C BK 2/95 Changed after advice by F. Villalobos, Cordoba
      RDN = 0.6* (720.*COSZEN-RBN-WATABS*COSZEN)

C Total potential visible radiation (estimated)
      RTOTV = RBV + RDV
C Total potential NIR radiation (estimated)
      RTOTN = RBN + RDN
C Total incoming radiation (measured)
      RMEAS = RADABV/FPAR
C Ratio of measured to potential solar radiation
      RATIO = RMEAS/(RTOTV + RTOTN)

      IF (RATIO.GT.PARMA) THEN
C         FBEAM(IHR,1) = RBV/RTOTV
         CALCFBMWN = RBV/RTOTV
      ELSE
C         FBEAM(IHR,1) = RBV/RTOTV* (1.- ((PARMA-RATIO)/PARMB)** (2./3.))
         CALCFBMWN = RBV/RTOTV* (1.- ((PARMA-RATIO)/PARMB)** (2./3.))
      END IF

C      IF (RATIO.GT.PARMC) THEN
C         FBEAM(IHR,2) = RBN/RTOTN
C      ELSE
C         FBEAM(IHR,2) = RBN/RTOTNN* (1.- ((PARMC-RATIO)/PARMD)** (2./3.))
C      END IF

       IF (CALCFBMWN.LT.0.000) CALCFBMWN = 0.00
C      IF (FBEAM(IHR,1).LT.0.000) FBEAM(IHR,1) = 0.00
C      IF (FBEAM(IHR,2).LT.0.000) FBEAM(IHR,2) = 0.00

C    *****************************************
C The following equation is from the Ph.D thesis of M.H. Steven,
C Univ. of Nottingham, 1977, entitled  'Angular distribution and
C interception of diffuse solar radiation'.
C   diffuse/global=-0.86*sunshine duration+0.99
         FSUN = 1.16*CALCFBMWN - 0.01
C         FSUN = 1.16*FBEAM(IHR,1) - 0.01
C Add the circumsolar diffuse components to the direct radiation
C for the sake of geometric simplicity.
C According to M.H. Steven this component is 51.0 per cent of the total
C diffuse.
      IF (FSUN.GT.0.8) THEN
         CALCFBMWN = CALCFBMWN + 0.51* (1.0-CALCFBMWN)*FSUN
C         FBEAM(IHR,1) = FBEAM(IHR,1) + 0.51* (1.0-FBEAM(IHR,1))*FSUN
C         FBEAM(IHR,2) = FBEAM(IHR,2) + 0.51* (1.0-FBEAM(IHR,2))*FSUN
      END IF

100   RETURN
      END  !CalcFBMWH

***************************************************************************

