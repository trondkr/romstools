SUBROUTINE fill(nx,ny,i1,i2,j1,j2,za,tx,critx,cor,mxs,rmask,error,nvalue)

! ** This routine originates from the USA.
! ** It was delivered (by 3.5" diskette)
! ** to DNMI in nov. 1990.

! ***** Re-edited at DNMI in January '91 by H.Engedahl.
! ***** Again re-edited at DNMI in March '92 by H.Engedahl to fit the
! ***** ECOM3D output configuration.
! ***** DNMI/FoU 25.08.1993 A.Foss : Recoded for speedup.

!// Solves Laplace's equation with Neumann boundary conditions
!// (dA/dn = 0) in RECTANGULAR coordinates by an iterative method to
!// fill in reasonable values at gridpoints containing values like
!// "undef".

!// NOTE: It is impossibel to make a really parallel (MPP) version of this
!//       routine. One MPP node should work on a complete field.
!//       Parallelization is done on a 'higher' level.  (A.Foss '98)

!-----------------------------------------------------------------------
!//     nx   = Array x dimension
!//     ny   = Array y dimension
!//   i1,i2  = Subarea in x direction to be filled
!//   j1,j2  = Subarea in y direction to be filled
!//     za   = Array to be filled (REAL)
!//     tx   = All values in array A GREATER than Tx are filled (REAL)
!//    critx = Criteria for relaxation, DEL**2 = CRIT
!//            (Usually 4 orders of magnitude DOWN from data in A)
!//     cor  = Coef. of overrelaxation, between +1.2 and +2.0
!//     mxs  = Max. allowed no. of scans in relaxation procedure.
!//   rmask  = Work array
!//   error  = Array containing the errors in the relaxation procedure.
!//   nvalue = No. of gridpoints with value (possibly 0) ... output
!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN)  :: nx, ny, i1, i2, j1, j2, mxs
REAL,    INTENT(IN)  :: tx,critx,cor
INTEGER, INTENT(OUT) :: nvalue

REAL                 :: za(nx,ny)
REAL                 :: rmask(nx,ny)
REAL                 :: error(nx,ny)
INTEGER              :: n, j, i, i1p1, i2m1, j1p1, j2m1, nnn, nbad
REAL                 :: suma, asuma, crit, crtest

n = 0
suma = 0.

DO j=j1,j2
DO i=i1,i2
  IF (za(i,j) < tx) THEN
    suma = suma + za(i,j)
    n = n+1
  END IF
END DO
END DO

nvalue = n

IF (n < 1) THEN
  PRINT *
  PRINT *,"******************  WARNING  *******************"
  PRINT *,"SUBROUTINE FILL : NO USEFUL DATA IN THE FIELD"
  PRINT *,"ALL DATA IN THE INTERIOR (SUBAREA) WERE > ",tx
  PRINT *,"******************  WARNING  *******************"
  PRINT *
  RETURN
END IF

suma = suma/n
asuma = 0.

DO j=j1,j2
DO i=i1,i2
  IF (za(i,j) < tx) THEN
    asuma = asuma + ABS(za(i,j)-suma)
    rmask(i,j) = 0.
  ELSE
    za(i,j) = suma
    rmask(i,j) = 1.
  END IF
END DO
END DO

asuma = asuma/n

crit = critx*asuma  ! af:  Using asuma not suma

!// The value of suma ( i.e. the MEAN value of the array A ) is filled in
!// the array A at all points with an initial value >= t. 
!// za(i,j) = suma may be regarded as the "first guess" in the iterative
!// method.

i1p1 = i1 + 1
i2m1 = i2 - 1
j1p1 = j1 + 1
j2m1 = j2 - 1

DO j=j1p1,j2m1
DO i=i1p1,i2m1
  rmask(i,j) = cor*rmask(i,j)
END DO
END DO

DO nnn=1,mxs

  DO j=j1p1,j2m1
    DO i=i1p1,i2m1
      error(i,j) = (za(i+1,j)+za(i,j+1)+za(i,j-1))*0.25-za(i,j)
    END DO
    DO i=i1p1,i2m1
      error(i,j) = error(i,j) + za(i-1,j)*0.25
      za(i,j) = za(i,j) + error(i,j)*rmask(i,j)
    END DO
  END DO

  !// Test convergence now and then (slow test loop)
  IF (nnn < mxs-5 .AND. MOD(nnn,10) == 0) THEN
    crtest = crit*cor
    nbad = 0
    j = j1
    DO WHILE (nbad == 0 .AND. j < j2m1)
      j=j+1
      DO i=i1p1,i2m1
        IF (ABS(error(i,j))*rmask(i,j) > crtest) nbad = 1
      END DO
    END DO
    IF (nbad == 0) GOTO 300
  END IF

  DO j=j1p1,j2m1
    za(i1,j) = za(i1,j) + (za(i1p1,j)-za(i1,j))*rmask(i1,j)
    za(i2,j) = za(i2,j) + (za(i2m1,j)-za(i2,j))*rmask(i2,j)
  END DO
  DO i=i1,i2
     za(i,j1) = za(i,j1) + (za(i,j1p1)-za(i,j1))*rmask(i,j1)
     za(i,j2) = za(i,j2) + (za(i,j2m1)-za(i,j2))*rmask(i,j2)
  END DO

ENDDO

300 CONTINUE

RETURN

END
