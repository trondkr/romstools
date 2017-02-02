
!***********************************************************************
!
! Author:  Jens Debernard
! Date:    20080606
! Purpose: Extract defined(?) values from a field into undef points. 
!          Undef points are replaced by the weighted mean of the defined 
!          neighbour points.
!
!***********************************************************************

! Arguments

! Local variabels

! Initialisation and start

subroutine creep_fill(finn,mask,ilo,ihi,jlo,jhi,maxit)
  use fill_mod
  implicit none


  integer, intent(in) :: ilo,jlo,ihi,jhi,maxit
  real, intent(in),    dimension(ilo:ihi,jlo:jhi) :: mask
  real, intent(inout), dimension(ilo:ihi,jlo:jhi) :: finn

!jd  real, intent(in),    dimension(:,:) :: mask
!jd  real, intent(inout), dimension(:,:) :: f
  
  real, dimension(ilo-1:ihi+1,jlo-1:jhi+1) :: f

  integer :: itr, pt, i,j, k,l
  real :: sumw, wginv

  f = 0.0
  f(ilo:ihi,jlo:jhi) = finn

  call fill_init(ilo,ihi,jlo,jhi,mask)

  do itr = 1,maxit
     npo = npn

     nm1 = n
     n = 3-nm1


     do pt=1,npo
        i = iind(pt,nm1); j=jind(pt,nm1); 
        sumw=0
        do k=-1,1
           do l=-1,1
              mwg(l,k) = rmask(i+l,j+k)*wg(l,k)
              sumw = sumw + mwg(l,k)
           end do
        end do
!        write(*,*) i,j, sumw


        imask(pt) = 0
        fval(pt) = 0.0
        if (sumw > 0.1) then
           sumw = 1.0/sumw
           imask(pt) = 1
           do k=-1,1
              do l=-1,1
                 wginv = mwg(l,k)*sumw
                 fval(pt) = fval(pt) + wginv*f(i+l,j+k)
!                 write(*,*) pt, wginv, sumw, fval(pt)
              end do
           end do
        end if
     end do

     npn=0
     do pt=1,npo
        i = iind(pt,nm1); j=jind(pt,nm1); 
        if ( imask(pt) /= 0 ) then
           f(i,j) = fval(pt)
           rmask(i,j) = 1.0
        else
           npn = npn + 1
           rmask(i,j) = 0.0
           iind(npn,n) = i
           jind(npn,n) = j
        end if
     end do

!     call printf(rmask,ilo,ihi,jlo,jhi,'rmask') 
!     call printf(f,ilo,ihi,jlo,jhi,'f-creep') 

!jd     write(*,*) itr, npn
     if (npn == 0) exit
  end do

  finn = f(ilo:ihi,jlo:jhi) 


  write(*,*) 'Itr, maxit = ', itr, maxit

end subroutine creep_fill
