 subroutine shiftphase(tide_Eamp,tide_Ephase,uamp, uphase,vamp,          &
                       vphase,imax,jmax,constit,reflat)

implicit none

real, dimension(imax+2,jmax+2) :: tide_Eamp
real, dimension(imax+2,jmax+2) :: tide_Ephase
real, dimension(imax+2,jmax+2) :: uamp
real, dimension(imax+2,jmax+2) :: uphase
real, dimension(imax+2,jmax+2) :: vamp
real, dimension(imax+2,jmax+2) :: vphase
integer :: imax, jmax,ierr,ierror,i,j
integer,parameter :: mtime=5
integer, dimension(mtime)::startime, stoptime
character(len=5) :: constit
real*8 ::vx,ux,fx,vp,up,fp
real*8 :: reflat,khma,khmb,khmm,nobsum
integer  ::   ihha,idda,imma,iyya,icca,ihhb,iddb,immb,iyyb,iccb,kha,imia,imib
integer  ::   kd,khb
integer  ::   khm
integer :: yy_start, mm_start, dd_start, hh_start, mi_start
integer ::  yy_end, mm_end, dd_end, hh_end, mi_end


call readsup(yy_start, mm_start, dd_start, hh_start, mi_start,    &
                     yy_end, mm_end, dd_end, hh_end, mi_end)

startime(1) = yy_start
startime(2) = mm_start
startime(3) = dd_start
startime(4) = hh_start
startime(5) = mi_start
stoptime(1) = yy_end
stoptime(2) = mm_end
stoptime(3) = dd_end
stoptime(4) = hh_end
stoptime(5) = mi_end


!C***********************************************************************
!C*  CASE INITIALIZATION
!C
!c       Decode the start and stop times

      icca = startime(1) / 100                  
      iyya = startime(1) - icca*100
      imma = startime(2)
      idda = startime(3)
      ihha = startime(4)
      imia = startime(5)
      iccb = stoptime(1) / 100
      iyyb = stoptime(1) - iccb*100
      immb = stoptime(2)
      iddb = stoptime(3)
      ihhb = stoptime(4)
      imib = stoptime(5)
      if ((idda.eq.0) .or. (iddb.eq.0)) then 
         ierror = -1
         write (*,*)                                          &
                '**shiftphase error: invalid time stamp'
         return
      end if

!C***********************************************************************
!C*   HERE THE NUMBER OF USEABLE OBSERVATIONS (NOBSU) & THE TIME OF THE
!C*   MIDDLE OBSERVATION ARE FOUND.
!C
      CALL GDAY(IDDA,IMMA,IYYA,ICCA,KD,ierr)
      if (ierr.lt.0) then
         ierror = -3
         return
      end if
      KHA=KD*24+IHHA
      ! khma er starttidspunktet i desimaldager
      khma=real(KHA)/24.0
      khma=khma+(imia/60.0)/24.0
      khma=khma*24.0
      CALL GDAY(IDDB,IMMB,IYYB,ICCB,KD,ierr)
      if (ierr.lt.0) then
         ierror = -3
         return
      end if

      KHB=KD*24+IHHB
      ! khmb er sluttidspunktet i desimaldager
      khmb=real(KHB)/24.0
      khmb=khmb+(imib/60.0)/24.0
      khmb=khmb*24.0


      nobsum=((khmb-khma-1-0)/2.0)*2.0+1.0
      ! khmm er midtidspunktet i desimaldager
      khmm=khma+NOBSUM/2
      khmm=khmm/24.0
      khma=khma/24.0


      call new_vuf(vx,ux,fx,khma,reflat,constit)
      call new_vuf(vp,up,fp,khmm,reflat,constit) ! only used for nodal corrections
      ! u and v are returned in cycles, so multiplying by 360 to get degrees
      vx=vx*360.0
      ux=ux*360.0
      vp=vp*360.0
      up=up*360.0

      do i=1,imax+2
         do j=1,jmax+2
            tide_Ephase(i,j)=tide_Ephase(i,j)-up-vx
            uphase(i,j)=uphase(i,j)-up-vx
            vphase(i,j)=vphase(i,j)-up-vx
            tide_Eamp(i,j)=tide_Eamp(i,j)*fp
            uamp(i,j)=uamp(i,j)*fp
            vamp(i,j)=vamp(i,j)*fp
            tide_Ephase(i,j)=mod(tide_Ephase(i,j),360.0)
            uphase(i,j)=mod(uphase(i,j),360.0)
            vphase(i,j)=mod(vphase(i,j),360.0)
         end do
      end do
return
 909  continue
      write (*, *) ' HT_ANA Error: constituent file not found'
      ierror = -4
      return
close(333)
end subroutine shiftphase


