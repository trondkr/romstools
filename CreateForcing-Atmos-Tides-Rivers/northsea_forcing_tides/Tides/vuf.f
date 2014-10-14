      subroutine vuf (konx,vux,fx,ierr)
C
C***********************************************************************
C*  THIS SUBROUTINE CALCULATES V (ASTRONOMICAL PHASE ARGUMENT), U AND F
C*  (NODAL MODULATION PHASE AND AMPLITUDE CORRECTIONS) FOR ALL CONSTITU-
C*  ENTS.
c*
C*	This October 1992 version also recalculates the constituent
c	frequencies for the middle of the analysis period

c       Moved variables in common to arguments in order to make 
c       creation of MEX-file easier.
c       Aug 1998, Bruce Hackett DNMI
c       STOP statements replaced with error returns.
c       Mar 1999, Bruce Hackett DNMI
c       Divided into individual subroutines
c       May 1999, Bruce Hackett DNMI
C
      integer     ierr
      integer     maxfreq
      character   konx*5, name*5
      real        vux, fx, xlat, freq

      integer     ntotal
      integer     ii(50),jj(50),kk(50),ll(50),mm(50),nn(50),nj(170)
      real        semi(50)
      character*5 kon(170)
      common      /vufmain/ ntotal,kon,ii,jj,kk,ll,mm,nn,semi,nj

      real        coef(320),f(170),vu(170)
      character*5 konco(320)
      common      /vufshal/ konco,coef,vu,f

      integer     k, lp
C
C***********************************************************************
C*  THE DIMENSION OF KON, VU, F, AND NJ SHOULD BE AT LEAST EQUAL TO THE
C*  TOTAL NUMBER OF POSSIBLE CONSTITUENTS (PRESENTLY 146), THE DIMENSION
C*  OF II, JJ, KK, LL, MM, NN AND SEMI SHOULD BE AT LEAST EQUAL TO THE
C*  NUMBER OF MAIN CONSTITUENTS (PRESENTLY 45), THE DIMENSION OF EE,
C*  LDEL, MDEL, NDEL, IR, AND PH SHOULD BE AT LEAST EQUAL TO THE TOTAL
C*  NUMBER OF SATELLITES TO ALL THE MAIN CONSTITUENTS PLUS THE NUMBER
C*  OF CONSTITUENTS WITH NO SATELLITES (PRESENTLY 162+8),
C*  AND THE DIMENSION OF KONCO, AND COEFF SHOULD BE AT LEAST EQUAL TO
C*  THE SUM FOR ALL SHALLOW WATER CONSTITUENTS OF THE NUMBER OF MAIN
C*  CONSTITUENTS FROM WHICH EACH IS DERIVED (PRESENTLY 251).
C***********************************************************************
C* GIVEN CONSTITUENT KONX , THE NODAL CORRECTIONS V+U AND F ARE RETURNED
C
      ierr = 0
      lp   = 6

      DO 20 K=1,NTOTAL
      IF(KON(K).eq.KONX) go to 40
20    CONTINUE
      WRITE(LP,30)KONX
   30 FORMAT(' VUF: Error on ',A5)
      ierr = -1
      return
40    VUX=VU(K)
      FX=F(K)
      RETURN
      end
C
C***********************************************************************
C***********************************************************************

      subroutine OPNVUF(KONX,VUX,FX,kr,ierr)

      integer     kr,ierr
      real        vux,fx
      character*5 konx

      integer     ntotal
      integer     ii(50),jj(50),kk(50),ll(50),mm(50),nn(50),nj(170)
      real        semi(50)
      character*5 kon(170)
      common      /vufmain/ ntotal,kon,ii,jj,kk,ll,mm,nn,semi,nj

      integer     ntidal,ldel(180),mdel(180),ndel(180),ir(180)
      real        ee(180),ph(180)
      common      /vufsat/  ntidal,ee,ldel,mdel,ndel,ir,ph

      real        coef(320),f(170),vu(170)
      character*5 konco(320)
      common      /vufshal/ konco,coef,vu,f

      integer     lp,j,k
      real*8      s0,h0,p0,enp0,pp0
      real*8      d1,h,pp,s,p,enp,dh,dpp,ds,dp,dnp,hh,tau,dtau
      integer     jbase,j1,jl,j4,k1
      character*5 kblank

      DATA KBLANK/'     '/

      lp   = 6
      ierr = 0

C***********************************************************************
C*  THE ASTRONOMICAL ARGUMENTS AND THEIR RATES OF CHANGE,
C*  S0,H0,P0,ENP0,PP0,DS,DH,DP,DNP,DPP,  ARE READ FROM TWO RECORDS IN
C*  THE FORMAT(5F13.10):
C*     S0  = MEAN LONGITUDE OF THE MOON (CYCLES) AT 000 ET 1/1/1976.
C*     H0  = MEAN LONGITUDE OF THE SUN.
C*     P0  = MEAN LONGITUDE OF THE LUNAR PERIGEE.
C*     ENP0= NEGATIVE OF THE MEAN LONGITUDE OF THE ASCENDING NODE.
C*     PP0 = MEAN LONGITUDE OF THE SOLAR PERIGEE (PERIHELION).
C*     DS,DH,DP,DNP,DPP ARE THEIR RESPECTIVE RATES OF CHANGE OVER A 365
C*     DAY PERIOD AS OF 000 ET 1/1/1976.
C
c	These values are no longer used though they are still
c	read in. More accurate polynomial approximations are 
c	now employed.
      READ(KR,50)S0,H0,P0,ENP0,PP0,DS,DH,DP,DNP,DPP
 50   FORMAT(5F13.10)
C
C***********************************************************************
C*  HERE THE MAIN CONSTITUENTS AND THEIR DOODSON NUMBERS ARE READ IN
C*  FORMAT (6X,A5,1X,6I3,F5.2,I4). THE VALUES ARE RESPECTIVELY
C*     KON    = CONSTITUENT NAME
C*  II,JJ,KK,LL,MM,NN = THE SIX DOODSON NUMBERS
C*     SEMI   = PHASE CORRECTION
C*     NJ     = THE NUMBER OF SATELLITES FOR THIS CONSTITUENT.
C*  THE END OF ALL MAIN CONSTITUENTS IS DENOTED BY A BLANK CARD.
c   There are 45 main constituents in the original package.   - bh
C
      JBASE=0
      DO 90 K=1,1000
      READ(KR,60)KON(K),II(K),JJ(K),KK(K),LL(K),MM(K),NN(K),SEMI(K),
     2 NJ(K)
   60 FORMAT(6X,A5,1X,6I3,F5.2,I4)
      IF(KON(K).eq.KBLANK) go to 100
      J1=JBASE+1
      IF(NJ(K).GE.1) GO TO 75
      NJ(K)=1
      JL=J1
      PH(J1)=0.
      EE(J1)=0.
      LDEL(J1)=0
      MDEL(J1)=0
      NDEL(J1)=0
      IR(J1)=0
      GO TO 90
   75 JL=JBASE+NJ(K)
C
C***********************************************************************
C*  IF NJ>0, INFORMATION ON THE SATELLITE CONSTITUENTS IS READ , THREE
C*  SATELLITES PER CARD, IN THE FORMAT (11X,3(3I3,F4.2,F7.4,1X,I1,1X)).
C*  FOR EACH SATELLITE THE VALUES READ ARE
C*     LDEL,MDEL,NDEL = THE CHANGES IN THE LAST THREE DOODSON NUMBERS
C*                      FROM THOSE OF THE MAIN CONSTITUENT.
C*     PH     = THE PHASE CORRECTION
C*     EE     = THE AMPLITUDE RATIO OF THE SATELLITE TIDAL POTENTIAL TO
C*              THAT OF THE MAIN CONSTITUENT.
C*     IR     = 1 IF THE AMPLITUDE RATIO HAS TO BE MULTIPLIED BY THE
C*                LATITUDE CORRECTION FACTOR FOR DIURNAL CONSTITUENTS
C*              2 IF THE AMPLITUDE RATIO HAS TO BE MULTIPLIED BY THE
C*                LATITUDE CORRECTION FACTOR FOR SEMI-DIURNAL CONSTI-
C*                TUENTS.
C*              OTHERWISE IF NO CORRECTION IS REQUIRED TO THE AMPLITUDE
C*                RATIO.
C
      READ(KR,80)(LDEL(J),MDEL(J),NDEL(J),PH(J),EE(J),IR(J),J=J1,JL)
   80 FORMAT((11X,3(3I3,F4.2,F7.4,1X,I1,1X)))
90    JBASE=JL
100   NTIDAL=K-1
C
C***********************************************************************
C*  THE SHALLOW WATER CONSTITUENTS AND THE MAIN CONSTITUENTS FROM WHICH
C*  THEY ARE DERIVED ARE READ IN HERE WITH THE FORMAT
C*  (6X,A5,I1,2X,4(F5.2,A5,5X)). THE VALUES ARE RESPECTIVELY
C*     KON    = NAME  OF THE SHALLOW WATER CONSTITUENT
C*     NJ     = NUMBER OF MAIN CONSTITUENTS FROM WHICH IT IS DERIVED.
C*     COEF,KONCO = COMBINATION NUMBER AND NAME OF THESE MAIN
C*                  CONSTITUENTS.
C*  THE END OF THESE CONSTITUENTS IS DENOTED BY A BLANK CARD.
c   There are 101 shallow water constituents in the original package.
c   Together with the main constituents we get 146 altogether.
C
      JBASE=0
      K1=NTIDAL+1
      DO 160 K=K1,1000
      J1=JBASE+1
      J4=J1+3
      READ(KR,130)KON(K),NJ(K),(COEF(J),KONCO(J),J=J1,J4)
  130 FORMAT(6X,A5,I1,2X,4(F5.2,A5,5X))
      IF(KON(K).eq.KBLANK) go to 170
160   JBASE=JBASE+NJ(K)
170   NTOTAL=K-1
      RETURN
      end
C
C***********************************************************************
C***********************************************************************

      subroutine SETVUF(KH,KONX,VUX,FX,XLAT,ierr)

      integer     kh,ierr
      real        vux,fx,xlat
      character*5 konx

      integer     mtot
      real        freq(170)
      character*5 name(170)
      common      /const/ name,freq,mtot

      integer     ntotal
      integer     ii(50),jj(50),kk(50),ll(50),mm(50),nn(50),nj(170)
      real        semi(50)
      character*5 kon(170)
      common      /vufmain/ ntotal,kon,ii,jj,kk,ll,mm,nn,semi,nj

      integer     ntidal,ldel(180),mdel(180),ndel(180),ir(180)
      real        ee(180),ph(180)
      common      /vufsat/  ntidal,ee,ldel,mdel,ndel,ir,ph

      real        coef(320),f(170),vu(170)
      character*5 konco(320)
      common      /vufshal/ konco,coef,vu,f

      integer     lp,j,k,l,lk,km1,int24,intdys
      real        pi,twopi,slat,vdbl,v,sumc,sums,rr,uudbl,uu
      real*8      s0,h0,p0,enp0,pp0
      real*8      d1,h,pp,s,p,enp,dh,dpp,ds,dp,dnp,hh,tau,dtau
      integer     jbase,j1,jl,k1,iv,iuu,iflag,kd0,ier
      integer     indx(170)
      character*5 kblank

C***********************************************************************
C*  NTIDAL IS THE NUMBER OF MAIN CONSTITUENTS
C*  NTOTAL IS THE NUMBER OF CONSTITUENTS (MAIN + SHALLOW WATER)
C*  FOR  THE GIVEN TIME KH, THE TABLE OF F AND V+U VALUES IS
C*  CALCULATED FOR ALL THE CONSTITUENTS.
C*     F IS THE NODAL MODULATION ADJUSTMENT FACTOR FOR AMPLITUDE
C*     U IS THE NODAL MODULATION ADJUSTMENT FACTOR FOR PHASE
C*     V IS THE ASTRONOMICAL ARGUMENT ADJUSTMENT FOR PHASE.
C
      ierr = 0
      PI=3.1415926536
      TWOPI=2.*PI
      SLAT=SIN(PI*XLAT/180.)
c      CALL GDAY(1,1,76,19,KD)
c      YEARS=(KH/24.D0-KD)/365.00D0
C
C***********************************************************************
C*  THE ASTRONOMICAL ARGUMENTS ARE CALCULATED BY LINEAR APPROXIMATION
C*  AT THE MID POINT OF THE ANALYSIS PERIOD.
C
c      S=S0+YEARS*DS
c      H=H0+YEARS*DH
c      P=P0+YEARS*DP
c      ENP=ENP0+YEARS*DNP
c      PP=PP0+YEARS*DPP
c	day number measured from January 0.5 1900 (i.e.,
c	1200 UT December 31, 1899
      d1=kh/24.d0
      call gday(31,12,99,18,kd0,ier)
      if (ier.lt.0) then
         ierr = -1
         return
      end if
      d1=d1-dfloat(kd0)-0.5d0
      call astr(d1,h,pp,s,p,enp,dh,dpp,ds,dp,dnp)
      INT24=24
      INTDYS=KH/INT24
      HH=dfloat(KH-INTDYS*INT24)
      TAU=HH/24.D0+H-S
	dtau=365.d0+dh-ds
C
C***********************************************************************
C*  ONLY THE FRACTIONAL PART OF A SOLAR DAY NEED BE RETAINED FOR COMPU-
C*  TING THE LUNAR TIME TAU.
C
      JBASE=0
      DO 210 K=1,NTIDAL
	do 209 l=1,mtot
	if(kon(k).eq.name(l)) then
      FREQ(l)=(II(K)*DTAU+JJ(K)*DS+KK(K)*DH+LL(K)*DP+MM(K)*DNP+
     1NN(K)*DPP)/(365.*24.)
	indx(k)=l
	end if
209	continue
      VDBL=II(K)*TAU+JJ(K)*S+KK(K)*H+LL(K)*P+MM(K)*ENP+NN(K)*PP+SEMI(K)
      IV=VDBL
      IV=(IV/2)*2
      V=VDBL-IV
      J1=JBASE+1
      JL=JBASE+NJ(K)
      SUMC=1.
      SUMS=0.
      DO 200 J=J1,JL
C
C***********************************************************************
C*  HERE THE SATELLITE AMPLITUDE RATIO ADJUSTMENT FOR LATITUDE IS MADE
C
      RR=EE(J)
      L=IR(J)+1
      GO TO (901,902,903),L
  902 RR=EE(J)*0.36309*(1.-5.*SLAT*SLAT)/SLAT
      GO TO 901
  903 RR=EE(J)*2.59808*SLAT
  901 CONTINUE
      UUDBL=LDEL(J)*P+MDEL(J)*ENP+NDEL(J)*PP+PH(J)
      IUU=UUDBL
      UU=UUDBL-IUU
      SUMC=SUMC+RR*COS(UU*TWOPI)
      SUMS=SUMS+RR*SIN(UU*TWOPI)
  200 CONTINUE
      F(K)=SQRT(SUMC*SUMC+SUMS*SUMS)
      VU(K)=V+ATAN2(SUMS,SUMC)/TWOPI
210   JBASE=JL
C
C***********************************************************************
C*  HERE F AND V+U OF THE SHALLOW WATER CONSTITUENTS ARE COMPUTED FROM
C*  THE VALUES OF THE MAIN CONSTITUENT FROM WHICH THEY ARE DERIVED.
C
      JBASE=0
      K1=NTIDAL+1
      IF(K1.GT.NTOTAL) RETURN
      DO 270 K=K1,NTOTAL
      F(K)=1.0
      VU(K)=0.0
      iflag=0
      do 269 lk=1,mtot
         if(kon(k).eq.name(lk)) then
            FREQ(lk)=0.
            iflag=1
            go to 268
         end if
 269  continue
 268  J1=JBASE+1
      JL=JBASE+NJ(K)
      DO 260 J=J1,JL
      KM1=K-1
      DO 240 L=1,KM1
      IF(KON(L).eq.KONCO(J)) go to 250
240   CONTINUE
      WRITE(LP,241)KONCO(J)
  241 FORMAT('   SETVUF: Error on ',A5)
      ierr = -1
      return
250   F(K)=F(K)*F(L)**ABS(COEF(J))
      VU(K)=VU(K)+COEF(J)*VU(L)
      if(iflag.eq.1) FREQ(lk)=FREQ(lk)+COEF(J)*FREQ(indx(L))
 260  continue
270   JBASE=JL
      RETURN
      END
