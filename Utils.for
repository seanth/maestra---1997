$DEBUG:'D' !DEBUG VERSION; ALL LINES BEGINNING WITH 'D' ARE COMPILED IN

C**********************************************************************
C UTILS.FOR
C
C This file contains all the 'utility' functions which are called throughout
C the other files. I have set the files up so that, when writing a test program
C which calls one of the program subroutines in file xxx.for, the program
C can be compiled by compiling test.for, xxx.for, and utils.for. 
C
C The functions are:
C IDATE50 - translates string-format date in days-since-1950 format
C IDATEYR - given string-format date, finds Julian day
C JDATE - given days-since-1950-format date, finds Julian day
C SUBERROR - error-handling subroutine
C BETA - implements beta function
C QUADM, QUADP - solve quadratic
C MEANVAL - find mean of array of values
C NIGHT - whether an hour qualifies as night-time or not
C**********************************************************************

C**********************************************************************
      FUNCTION IDATE50(STRDATE)
C This function translates a string-format date (DD/MM/YY) into the number 
C of days since 1/1/1950. Provides a numerical way of handling dates.
C**********************************************************************

      CHARACTER*10 STRDATE

C Get the day of year and the year number from function IDATEYR.
      IDAY = IDATEYR(STRDATE,IYEAR)

C Calculate how many days in full years have passed since 1950.
      IYRD = 365*(IYEAR - 50)
      IYRD = IYRD + (IYEAR - 49)/4
      IDATE50 = IYRD + IDAY

      RETURN
      END !Idate50


C**********************************************************************
      FUNCTION IDATEYR(STRDATE,IYEAR)
C Given a date strdate, this function calculates the number of the day 
C from 1 on 1st Jan to 365/366 on 31st Dec.
C**********************************************************************

      CHARACTER*10 STRDATE
      LOGICAL LEAPYR
      INTEGER IFD(12)
      DATA IFD/0,31,59,90,120,151,181,212,243,273,304,334/

      READ (STRDATE,10) IDAY,IMON,IYEAR
10    FORMAT(T1,I2,T4,I2,T7,I2)
C Code to handle passage to 2000: note this will only work until 2050
C I hope this model has been superseded by then!
	IF (IYEAR.LT.50) IYEAR = IYEAR+100

      LEAPYR = .FALSE.
      IF (4.0* (IYEAR/4).EQ.IYEAR) LEAPYR = .TRUE.

      IDATEYR = IFD(IMON) + IDAY
      IF (LEAPYR .AND. (IMON.GE.3)) IDATEYR = IDATEYR + 1

      RETURN
      END !Idateyr


C**********************************************************************
      FUNCTION MONTH(IDATE50)
C This function returns the month given a date in days-since-1950 format.
C**********************************************************************

	JD = JDATE(IDATE50)
	MONTH = JD/12 + 1

	RETURN
	END !Mon


C**********************************************************************
      FUNCTION JDATE(IDATE50)
C This function returns the julian day given a date in days-since-1950
C format.
C**********************************************************************

      IYR = IDATE50 / 365
      IYRD = 365*IYR
      IYRD = IYRD + (IYR - 1)/4 
      JDATE = IDATE50 - IYRD + 1

      RETURN
      END !Jdate


C**********************************************************************
      SUBROUTINE SUBERROR(MESSAGE,IFLAG,IOERROR)
C The error-handling subroutine. When called, writes MESSAGE to the
C error file (UERROR). Uses IFLAG to determine whether to terminate program
C or not. IOERROR is a FORTRAN error number, reported in the error file. 
C**********************************************************************

      INCLUDE 'MAESTCOM'
      CHARACTER MESSAGE *(*)

      WRITE(UERROR,*) MESSAGE
      IF (IOERROR.NE.0) WRITE(UERROR,20) IOERROR
10    FORMAT (A80)
20    FORMAT ('FORTRAN ERROR CODE NO: ',I10)

      IF (IFLAG.EQ.IFATAL) THEN
        STOP 'FATAL ERROR - SEE MAESERR.DAT FOR DETAILS.'
      END IF
      
      RETURN
      END !SubError


C**********************************************************************
      FUNCTION BETA(BC1,BC2,BC3,BC4,RP)
C This a function to define the vertical and horizontal leaf area
C density distribution.
C BM 11/99 added fourth parameter BC4
C**********************************************************************

      IF (RP.EQ.0.0) RP=0.0000001
      IF (RP.EQ.1.0) RP=0.9999999
	IF (RP.GE.BC4) THEN
	  BETA = 0.0
	ELSE
        BETA = BC1* (RP**BC2)* ((BC4-RP)**BC3)
	END IF

      RETURN
      END !Beta


C**********************************************************************
      REAL FUNCTION QUADM(A,B,C,IQERROR)
C Solves the quadratic equation - finds smaller root.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL A,B,C

      IQERROR = 0

      IF ((B*B - 4.*A*C).LT.0.0) THEN
        CALL SUBERROR('WARNING:IMAGINARY ROOTS IN QUADRATIC',
     &  IWARN,0)
        IQERROR = 1
        QUADM = 0.0
      ELSE
        IF (A.EQ.0.0) THEN
          IF (B.EQ.0.0) THEN
            QUADM = 0.0
            IF (C.NE.0.0) 
     &        CALL SUBERROR('ERROR: CANT SOLVE QUADRATIC',IWARN,0)
          ELSE
            QUADM = -C/B
          END IF
        ELSE
          QUADM = (- B - SQRT(B*B - 4*A*C)) / (2.*A)
        END IF
      END IF

      RETURN
      END !QuadM


C**********************************************************************
      REAL FUNCTION QUADP(A,B,C,IQERROR)
C Solves the quadratic equation - finds larger root.
C**********************************************************************

      INCLUDE 'MAESTCOM'
      REAL A,B,C

      IQERROR = 0

      IF ((B*B - 4.*A*C).LT.0.0) THEN
        CALL SUBERROR('WARNING:IMAGINARY ROOTS IN QUADRATIC',
     &  IWARN,0)
        IQERROR = 1
        QUADP = 0.0
      ELSE
        IF (A.EQ.0.0) THEN
          IF (B.EQ.0.0) THEN
            QUADP = 0.0
            IF (C.NE.0.0) 
     &        CALL SUBERROR('ERROR: CANT SOLVE QUADRATIC',IWARN,0)
          ELSE
            QUADP = -C/B
          END IF
        ELSE
          QUADP = (- B + SQRT(B*B - 4*A*C)) / (2.*A)
        END IF
      END IF

      RETURN
      END !QuadP

         
C**********************************************************************
      REAL FUNCTION MEANVAL(ARR,NUMVAL)
C Finds mean value of an array of NUMVAL values.
C**********************************************************************

      INCLUDE 'maestcom'
	REAL ARR(MAXP)

      MEANVAL = 0.0
	DO 10 I = 1,NUMVAL
	  MEANVAL = MEANVAL + ARR(I)
10    CONTINUE
      MEANVAL = MEANVAL/REAL(NUMVAL)

	RETURN
	END ! MeanVal


C**********************************************************************
      FUNCTION NIGHT(ZEN,PAR)
C Function to determine whether a particular hour is at night or not
C**********************************************************************

      IF ((ABS(ZEN).LE.1.57).AND.(PAR.GT.0.1)) THEN
        NIGHT = 0
      ELSE
        NIGHT = 1
      END IF
    
      RETURN
      END !Night


