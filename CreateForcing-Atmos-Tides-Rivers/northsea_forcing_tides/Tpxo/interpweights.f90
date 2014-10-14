subroutine interpweights(lon_i, lat_i,nx,ny, lon_o, lat_o,imax, jmax, r1, r2, r3, r4, i1, j1,topo)
  implicit none
  integer, intent(in) :: nx, ny, imax, jmax
  real, intent(in), dimension(ny,nx) :: lon_i, lat_i
  integer, intent(in), dimension(ny,nx) ::  topo
  real, intent(in), dimension(imax,jmax) :: lon_o, lat_o
  integer, intent(inout), dimension(imax,jmax) :: i1, j1
  real, intent(inout), dimension(imax,jmax) :: r1, r2, r3, r4
  real :: rsum, x,y, rx, ry
  integer :: i,j

  do j=1,jmax
     do i=1, imax

        ! For h-points
        x=(lon_o(i,j) - lon_i(1,1))/(lon_i(1,2)-lon_i(1,1)) + 1.0

        y=(lat_o(i,j) - lat_i(1,1))/(lat_i(2,1)-lat_i(1,1)) + 1.0

        ! Beware of negative x values! (Longitude wrap in global grids)
        if (x.lt.0.0) x=x+nx
!!$
        i1(i,j)=ifix(x)
        j1(i,j)=ifix(y)

        rx=(x-ifix(x))
        ry=(y-ifix(y))


        r1(i,j) = (1.0-rx)*(1.0-ry)
        r2(i,j) =      rx *(1.0-ry)
        r3(i,j) =      rx *     ry
        r4(i,j) = (1.0-rx)*     ry 

!!$        ! Check for landpoints!
!!$        r1(i,j) = r1(i,j)*real(topo(j1(i,j),i1(i,j)))
!!$        r2(i,j) = r2(i,j)*real(topo(j1(i,j),i1(i,j)+1))
!!$        r3(i,j) = r3(i,j)*real(topo(j1(i,j)+1,i1(i,j)+1))
!!$        r4(i,j) = r4(i,j)*real(topo(j1(i,j)+1,i1(i,j)))
!!$
!!$        rsum = r1(i,j) + r2(i,j) + r3(i,j) + r4(i,j)
!!$        if (rsum.gt.0.0) then
!!$           rsum=1.0/rsum
!!$           r1(i,j) = r1(i,j)*rsum
!!$           r2(i,j) = r2(i,j)*rsum
!!$           r3(i,j) = r3(i,j)*rsum
!!$           r4(i,j) = r4(i,j)*rsum
!!$        end if

     end do
  end do

end subroutine interpweights
