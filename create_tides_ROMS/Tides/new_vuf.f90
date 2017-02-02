subroutine new_vuf(vx,ux,fx,dd,xlat,constituent)
  implicit none
      ! Parameters
      integer, parameter   ::   maxfreq=170
      real,parameter       ::   pi=3.141592653589793      
      complex,parameter    ::   i=(0.0,1.0)

      ! Outout variables
      real*8      ::   fx,vx,ux

      ! Input variables        
      real*8      ::   xlat, dd 
      character(len=5)     ::   constituent

      ! Variables for main constituents
      integer,dimension(maxfreq)            ::   ii=0,jj=0,kk=0,ll=0,mm=0,nn=0 
      integer,dimension(maxfreq)            ::   ikmpr=0,nsat=0
      integer,dimension(maxfreq)            ::   nshallow=0,ishallow=0
      real,dimension(maxfreq)               ::   freq=0,df=0,semi=0
      character(len=5),dimension(maxfreq)   ::   name,kmpr,kon

      ! Variables for satellite data
      integer,dimension(180)                ::   ldel=0,mdel=0,ndel=0
      integer,dimension(180)                ::   ilatfac=0,iconst=0
      real,dimension(180)                   ::   amprat=0,phcorr=0

      ! Variables for shallow water constituents 
      integer,dimension(320)                ::   shallow_iconst=0,shallow_iname=0
      real,dimension(320)                   ::   shallow_coef=0


      ! Temporary variables for file input
      integer                          ::   iii,jjj,kkk,lll,mmm,nnn,nnj 
      real                             ::   nsemi
      real,dimension(4)                ::   ncoef
      character(len=5)                 ::   nkon
      character(len=5),dimension(4)    ::   nkonco

      ! Variables for calculation of output arguments
      real,dimension(146)              ::    f=0,u=0,v=0
      real,dimension(180)              ::    rr=0,uu=0
      complex,dimension(146)           ::    fsum=0
      ! Temporary variables used in the calculations
      real                             ::    dum,dumm,tmp2
      complex                          ::    tmp
      ! Loop variables
      integer   ::   t,j,k
      integer   ::   ntotal,ntidal,nsatellite,m
      integer   ::   jbase,j1,jend,jstart

      ! The astronomical arguments
      real,dimension(6)   ::   astro,ader

      ! Miscellaneous
      integer             ::   konx,kr,ierr,ier
      real                ::   slat
      character(len=5)    ::   kblank
      

!***************************************************************************
!***************************************************************************
 
      xlat=55
      kr=82
      DATA kblank/'     '/
      if (kr.gt.0) then
         open(unit=kr,file='constits.dat'                                &
              ,status='old'                                               &
              ,form='formatted')
      end if

      do k=1,maxfreq
         READ(kr,1010) name(k),freq(k),kmpr(k)
1010     FORMAT(4X,A5,3X,F13.10,4X,A5)

         if (name(k).eq.constituent) then
            konx=k
         elseif(name(k).eq.kblank) then
            m=k-1
            exit
         end if

      end do

      do k=1,m
         do t=1,m
           
            if (kmpr(k)==name(t))then
               ikmpr(k)=t
               df(k)=abs(freq(k)-freq(t))
               exit
            end if

         end do
      end do
      df(1)=0

      ierr = 0

!***********************************************************************
!  THE ASTRONOMICAL ARGUMENTS AND THEIR RATES OF CHANGE,
!  S0,H0,P0,ENP0,PP0,DS,DH,DP,DNP,DPP,  ARE READ FROM TWO RECORDS IN
!  THE FORMAT(5F13.10):
!
!	These values are no longer used though they are still
!	read in. More accurate polynomial approximations are 
!	now employed.
      read(kr,50)dum,dum,dum,dum,dum,dum,dum,dum,dum,dum
50    FORMAT(5F13.10)

!***********************************************************************
!  HERE THE MAIN CONSTITUENTS AND THEIR DOODSON NUMBERS ARE READ IN
!  FORMAT (6X,A5,1X,6I3,F5.2,I4). THE VALUES ARE RESPECTIVELY
!     kon    = CONSTITUENT NAME
!  ii,jj,kk,ll,mm,nn = THE SIX DOODSON NUMBERS
!     semi   = PHASE CORRECTION
!     nj     = THE NUMBER OF SATELLITES FOR THIS CONSTITUENT.
!  THE END OF ALL MAIN CONSTITUENTS IS DENOTED BY A BLANK CARD.
!  There are 45 main constituents in the original package.   - bh
!
      nsatellite=0

      jbase=1
      jend=0

      do k=1,1000
         read(kr,60)nkon,iii,jjj,kkk,lll,mmm,nnn,nsemi,nnj
60       FORMAT(6X,A5,1X,6I3,F5.2,I4)
        
         if(nkon.eq.kblank)then
            ntidal=k-1
            exit
         end if
         
         do t=1,m
            if (nkon.eq.name(t)) then
               j1=t
               exit
            end if
         end do

         kon(j1)=nkon
         ii(j1)=iii
         jj(j1)=jjj
         kk(j1)=kkk
         ll(j1)=lll
         mm(j1)=mmm
         nn(j1)=nnn
         semi(j1)=nsemi
         nsat(j1)=nnj
  
         if(nsat(j1).ge.1)then 
            nsatellite=nsatellite+nsat(j1)
            jend=jend+nsat(j1)
            !******************************************************************
            !  IF NJ>0, INFORMATION ON THE SATELLITE CONSTITUENTS IS READ , 
            !  THREE SATELLITES PER CARD, IN THE FORMAT 
            ! (11X,3(3I3,F4.2,F7.4,1X,I1,1X)).
            !  FOR EACH SATELLITE THE VALUES READ ARE
            !  LDEL,MDEL,NDEL = THE CHANGES IN THE LAST THREE DOODSON NUMBERS
            !                   FROM THOSE OF THE MAIN CONSTITUENT.
            !     ph     = THE PHASE CORRECTION
            !     ee     = THE AMPLITUDE RATIO OF THE SATELLITE TIDAL POTENTIAL
            !              TO THAT OF THE MAIN CONSTITUENT.
            !     ir     = 1 IF THE AMPLITUDE RATIO HAS TO BE MULTIPLIED BY THE
            !              LATITUDE CORRECTION FACTOR FOR DIURNAL CONSTITUENTS
            !              2 IF THE AMPLITUDE RATIO HAS TO BE MULTIPLIED BY THE
            !              LATITUDE CORRECTION FACTOR FOR SEMI-DIURNAL CONSTI-
            !              TUENTS.
            !              OTHERWISE IF NO CORRECTION IS REQUIRED TO THE 
            !              AMPLITUDE RATIO
            read(kr,80)(ldel(j),mdel(j),ndel(j),phcorr(j),amprat(j),ilatfac(j),j=jbase,jend)
80          FORMAT((11X,3(3I3,F4.2,F7.4,1X,I1,1X)))
            do t=jbase,jend
               iconst(t)=j1     
            end do
            jbase=jend+1
         end if         
      end do

!
!***********************************************************************
!  THE SHALLOW WATER CONSTITUENTS AND THE MAIN CONSTITUENTS FROM WHICH
!  THEY ARE DERIVED ARE READ IN HERE WITH THE FORMAT
!  (6X,A5,I1,2X,4(F5.2,A5,5X)). THE VALUES ARE RESPECTIVELY
!     KON    = NAME  OF THE SHALLOW WATER CONSTITUENT
!     NJ     = NUMBER OF MAIN CONSTITUENTS FROM WHICH IT IS DERIVED.
!     COEF,KONCO = COMBINATION NUMBER AND NAME OF THESE MAIN
!                  CONSTITUENTS.
!  THE END OF THESE CONSTITUENTS IS DENOTED BY A BLANK CARD.
!   There are 101 shallow water constituents in the original package.
!   Together with the main constituents we get 146 altogether.
!
      jbase=0
      do k=1,1000
         read(kr,130)nkon,nnj,(ncoef(j),nkonco(j),j=1,4)
130      FORMAT(6X,A5,I1,2X,4(F5.2,A5,5X))
      
         if(nkon.eq.kblank) then
            ntotal=k-1
           exit
           
        end if
        
        do t=1,m
           if (nkon.eq.name(t)) then
              j1=t
              exit
           end if
        end do

        nshallow(j1)=nnj
        ishallow(j1)=jbase+1

        do t=1,nnj
           jbase=jbase+1
           shallow_coef(jbase)=ncoef(t)
           shallow_iconst(jbase)=j1

           do kkk=1,m

              if (nkonco(t).eq.name(kkk)) then
                 shallow_iname(jbase)=kkk
                 exit
              end if

           end do

        end do

     end do
     close(unit=KR)


! Calculate astronomical arguments at mid-point of data time series

     call new_astro(dd,astro,ader)
  
! Phase relative to Greenwich (in units of cycles, presumeably).
! (This only returns values when we have doodson#s, i.e., not for the 
! shallow water components, but these will be computed later.)
!v=rem( const.doodson*astro+const.semi, 1);

     do k=1,m
        dum=ii(k)*astro(1) + jj(k)*astro(2) + kk(k)*astro(3) +         &
             ll(k)*astro(4) + mm(k)*astro(5) + nn(k)*astro(6) +         &
             semi(k)
        v(k)=dum-floor(dum)
     end do

!     call new_astro(dd,astro,ader)
     ! Apparently the second-order terms in the tidal potential go to zero
     ! at the equator, but the third-order terms do not. Hence when trying
     ! to infer the third-order terms from the second-order terms, the 
     ! nodal correction factors blow up. In order to prevent this, it is 
     ! assumed that the equatorial forcing is due to second-order forcing 
     ! OFF the equator, from about the 5 degree location. Latitudes are 
     ! hence (somewhat arbitrarily) forced to be no closer than 5 deg to 
     ! the equator.
     

     if (abs(xlat).lt.5.0) then
        if (xlat<0) then
           xlat=-5.0
        else
           xlat=5.0
        end if
     end if


     slat=sin(pi*xlat/180.0)

     ! satellite amplitude ratio adjustment for latitude

     rr=amprat  ! no amplitude correction
    
     do k=1,nsatellite
        if (ilatfac(k).eq.1) then
           ! latitude correction for diurnal constituents
           rr(k)=rr(k)*0.36309*(1.0-5.0*slat*slat)/slat
        elseif (ilatfac(k).eq.2) then
           ! latitude correction for semi-diurnal constituents
           rr(k)=rr(k)*2.59808*slat
        end if
     end do

     ! calculate nodal amplitude and phase corrections
     do k=1,nsatellite
        dum=ldel(k)*astro(4)+mdel(k)*astro(5)+ndel(k)*astro(6)+phcorr(k)
        uu(k)=dum-floor(dum)
     end do


     ! Sum up all of the satellite factors for all satellites
     fsum = 0  ! emsures that all elements of fsum equals zero before we start
     do k=1,nsatellite
        fsum(iconst(k))=fsum(iconst(k))+rr(k)*exp(i*2.0*pi*uu(k))
     end do

     fsum=fsum+(1,0)

     f=abs(fsum)

     do k=1,m
        tmp2=aimag(fsum(k))
        tmp=tmp2*i
        tmp=fsum(k)-tmp
        u(k)=atan2(tmp2,abs(tmp))/(2*pi)
     end do

     ! Compute amplitude and phase corrections for shallow water constituents
     do k=1,m
 
        if(nshallow(k).gt.0) then
           jbase=ishallow(k)
           jend=jbase+nshallow(k)-1
           tmp2=1
           dum=0
           dumm=0
           do t=jbase,jend
              tmp2=f(shallow_iname(t))**abs(shallow_coef(t))*tmp2
              dum=u(shallow_iname(t))*shallow_coef(t)+dum
              dumm=v(shallow_iname(t))*shallow_coef(t)+dumm
           end do
           f(k)=tmp2
           u(k)=dum      
           v(k)=dumm
        end if

     end do
     fx=f(konx)
     ux=u(konx)
     vx=v(konx)

   end subroutine new_vuf
