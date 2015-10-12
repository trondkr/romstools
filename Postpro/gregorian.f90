!=======================================================================
! Subroutine to find Gregorian date from julian day number
SUBROUTINE gregorian (julian, year, month, day, hour)

IMPLICIT NONE

real,    intent(in)  :: julian
integer, intent(out) :: year, month, day, hour

integer :: i, j, k, l, n

l = julian + 68569
n = 4*l/146097
l = l - (146097*n+3)/4
i = 4000*(l+1)/1461001
l = l - 1461*i/4 + 31
j = 80*l/2447
k = l - 2447*j/80
l = j/11
j = j+2-12*l
i = 100*(n-49) + i + l

year = i
month = j
day = k
hour = NINT((julian - INT(julian))*24)

RETURN

END SUBROUTINE gregorian

