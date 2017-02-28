      subroutine gday(idd,imm,iyy,icc,kd,ierr)
C!
C!  GIVEN DAY,MONTH,YEAR AND CENTURY(EACH 2 DIGITS), GDAY RETURNS
C!  THE DAY#, KD BASED ON THE GREGORIAN CALENDAR.
C!  THE GREGORIAN CALENDAR, CURRENTLY 'UNIVERSALLY' IN USE WAS
C!  INITIATED IN EUROPE IN THE SIXTEENTH CENTURY. NOTE THAT GDAY
C!  IS VALID ONLY FOR GREGORIAN CALENDAR DATES.
C
C   KD=1 CORRESPONDS TO JANUARY 1, 0000
C       
c       Note that the Gregorian reform of the Julian calendar 
c       omitted 10 days in 1582 in order to restore the date
c       of the vernal equinox to March 21 (the day after
c       Oct 4, 1582 became Oct 15, 1582), and revised the leap 
c       year rule so that centurial years not divisible by 400
c       were not leap years.
c
C   THIS ROUTINE WAS WRITTEN BY EUGENE NEUFELD, AT IOS, IN JUNE 1990.
C
c       STOP statements replaced by error returns.
c       Mar 1999, Bruce Hackett DNMI

      INTEGER NDP(13)
      INTEGER NDM(12)
      integer     idd,imm,iyy,icc,kd,lp,kkd,jfh,jcc,jfy,kk,jyy,jyyy
      integer     i,j,k,l
      integer     ierr
      DATA NDP/0,31,59,90,120,151,181,212,243,273,304,334,365/
      DATA NDM/31,28,31,30,31,30,31,31,30,31,30,31/
C!
      LP = 6
      ierr = 0
C!  TEST FOR INVALID INPUT:
      IF(ICC.LT.0)THEN
         WRITE(LP,5000)ICC
         ierr = -1
         return
      ENDIF
      IF(IYY.LT.0.OR.IYY.GT.99)THEN
         WRITE(LP,5010)IYY
         ierr = -1
         return
      ENDIF
      IF(IMM.LE.0.OR.IMM.GT.12)THEN
         WRITE(LP,5020)IMM
         ierr = -1
         return
      ENDIF
      IF(IDD.LE.0)THEN
         WRITE(LP,5030)IDD
         ierr = -1
         return
      ENDIF
      IF(IMM.NE.2.AND.IDD.GT.NDM(IMM))THEN
         WRITE(LP,5030)IDD
         ierr = -1
         return
      ENDIF
      IF(IMM.EQ.2.AND.IDD.GT.29)THEN
         WRITE(LP,5030)IDD
         ierr = -1
         return
      ENDIF
      IF(IMM.EQ.2.AND.IDD.GT.28.AND.((IYY/4)*4-IYY.NE.0.OR.(IYY.EQ.0.AND
     .    .(ICC/4)*4-ICC.NE.0)))THEN
         WRITE(LP,5030)IDD
         ierr = -1
         return
      ENDIF
5000  FORMAT(' INPUT ERROR. ICC = ',I7)
5010  FORMAT(' INPUT ERROR. IYY = ',I7)
5020  FORMAT(' INPUT ERROR. IMM = ',I7)
5030  FORMAT(' INPUT ERROR. IDD = ',I7)
C!
C!  CALCULATE DAY# OF LAST DAY OF LAST CENTURY:
      KD = ICC*36524 + (ICC+3)/4
C!
C!  CALCULATE DAY# OF LAST DAY OF LAST YEAR:
      KD = KD + IYY*365 + (IYY+3)/4
C!
C!  ADJUST FOR CENTURY RULE:
C!  (VIZ. NO LEAP-YEARS ON CENTURYS EXCEPT WHEN THE 2-DIGIT
C!  CENTURY IS DIVISIBLE BY 4.)
      IF(IYY.GT.0.AND.(ICC-(ICC/4)*4).NE.0) KD=KD-1
C!  KD NOW TRULY REPRESENTS THE DAY# OF THE LAST DAY OF LAST YEAR.
C!
C!  CALCULATE DAY# OF LAST DAY OF LAST MONTH:
      KD = KD + NDP(IMM)
C!
C!  ADJUST FOR LEAP YEARS:
      IF(IMM.GT.2.AND.((IYY/4)*4-IYY).EQ.0.AND.((IYY.NE.0).OR.
     .   (((ICC/4)*4-ICC).EQ.0)))   KD=KD+1
C!  KD NOW TRULY REPRESENTS THE DAY# OF THE LAST DAY OF THE LAST
C!  MONTH.
C!
C!  CALCULATE THE CURRENT DAY#:
      KD = KD + IDD
      RETURN
C!
C!
      ENTRY DMY(IDD,IMM,IYY,ICC,KD,ierr)
C!
C!  GIVEN THE (GREGORIAN) DAY#, KD, AS CALCULATED ABOVE IN THIS ROUTINE,
C!  ENTRY DMY RETURNS THE (GREGORIAN) DAY, MONTH, YEAR AND CENTURY.
C!
      lp = 6

C!  TEST FOR VALID INPUT:
      IF(KD.LE.0) WRITE(LP,5040)KD
5040  FORMAT(' DMY: KD = ',I7,'  INVALID INPUT.')
C!
C!  SAVE KD
      KKD=KD
C!  CALCULATE ICC AND SUBTRACT THE NUMBER OF DAYS REPRESENTED BY ICC
C!  FROM KKD
C!  JFH IS THE NUMBER OF 400 YEAR INTERVALS UP TO KKD
C!  JCC IS THE NUMBER OF ADDITIONAL CENTURIES UP TO KKD
      JFH = KKD/146097
      KKD = KKD - JFH*146097
      IF(KKD.LT.36525)THEN
         JCC = 0
      ELSE
         KKD = KKD - 36525
         JCC = 1 + KKD/36524
         KKD = KKD - (JCC-1)*36524
      END IF
      ICC = 4*JFH + JCC
      IF(KKD.EQ.0)THEN
         ICC = ICC-1
         IYY = 99
         IMM = 12
         IDD = 31
         RETURN
      ENDIF
C!
C!  CALCULATE IYY. JFY IS THE NUMBER OF FOUR YEAR INTERVALS IN THE
C!  CURRENT CENTURY. THE FIRST FOUR YEAR INTERVAL IS SHORT (1460 DAYS
C!  RATHER THAN 1461)IF THE CURRENT CENTURY IS NOT DIVISIBLE BY 4, AND
C!  IN THIS CASE JCC.NE.0 AS CALCULATED ABOVE.
C!
C!  CALCULATE JFY:
      JFY = 0
      IF(JCC.EQ.0)GOTO 10
      IF(KKD.LT.1460)GOTO 10
      JFY = 1
      KKD = KKD - 1460
10    KK = KKD/1461
      JFY = JFY + KK
      KKD = KKD - KK*1461
C!
C!  CALCULATE JYY, THE REMAINING YEARS OF THE CURRENT CENTURY UP TO THE
C!  CURRENT DAY:
      JYY = 0
C!  THE NEXT YEAR IS NOT A LEAP YEAR IF JFY=0 AND JCC.NE.0.
      IF(JFY.EQ.0.AND.JCC.NE.0)GOTO 20
      IF(KKD.LT.366)GOTO 30
      JYY = 1
      KKD = KKD - 366
20    JYYY = KKD/365
      JYY = JYY + JYYY
      KKD = KKD - JYYY*365
30    IYY = 4*JFY + JYY
      IF(KKD.EQ.0) THEN
         IYY=IYY-1
         IMM=12
         IDD=31
         RETURN
      END IF
C!
C!  SET L=1 IF WE HAVE A LEAP YEAR.
      L=0
      IF(IYY-(IYY/4)*4.NE.0)GOTO 40
      IF(IYY.EQ.0.AND.(ICC-(ICC/4)*4).NE.0)GOTO 40
      L=1
C!
C!  CALCULATE IMM AND IDD
40    IF(KKD.GT.31) GOTO 50
      IMM=1
      IDD=KKD
      RETURN
C!
50    IF(KKD.GT.59)GOTO 60
      IMM = 2
      IDD = KKD-31
      RETURN
C!
60    IF(KKD.GT.60)GOTO 70
      IF(L.EQ.0)GOTO 70
      IMM = 2
      IDD = 29
      RETURN
C!
70    IF(L.EQ.1) KKD=KKD-1
      DO 80 I=4,13
         IF(KKD.GT.NDP(I))GOTO 80
         IMM = I-1
         IDD = KKD - NDP(I-1)
         RETURN
C!
80    CONTINUE
90    WRITE(LP,5050)
5050  FORMAT(' ERROR IN DMY.')
      ierr = -1
      return
      END
