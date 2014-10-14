    
!C                                                                       
!C ********************************************************************
!C
subroutine fill(nx,ny,i1,i2,j1,j2,za,tx,critx,cor,mxs, &
                 rmask,error,nvalue)
  !!
  !! ** This routine originates from the USA.
  !! ** It was delivered (by 3.5" diskette)
  !! ** to DNMI in nov. 1990.
  !!
  !! ***** Re-edited at DNMI in January '91 by H.Engedahl.
  !! ***** Again re-edited at DNMI in March '92 by H.Engedahl to fit the
  !! ***** ECOM3D output configuration.
  !! ***** DNMI/FoU 25.08.1993 A.Foss : Recoded for speedup.
  !!
  !! Solves Laplace's equation with Neumann boundary conditions
  !! (dA/dn = 0) in RECTANGULAR coordinates by an iterative method to
  !! fill in reasonable values at gridpoints containing values like
  !! "undef".
  !!
  !! NOTE: It is impossibel to make a really parallel (MPP) version of this
  !!       routine. One MPP node should work on a complete field.
  !!       Parallelization is done on a 'higher' level.  (A.Foss '98)
  !!
  !!-----------------------------------------------------------------------
  !!     NX   = Array x dimension
  !!     NY   = Array y dimension
  !!   I1,I2  = Subarea in x direction to be filled
  !!   J1,J2  = Subarea in y direction to be filled
  !!     ZA   = Array to be filled (REAL)
  !!     Tx   = All values in array A GREATER than Tx are filled (REAL)
  !!    CRITX = Criteria for relaxation, DEL**2 = CRIT
  !!            (Usually 4 orders of magnitude DOWN from data in A)
  !!     COR  = Coef. of overrelaxation, between +1.2 and +2.0
  !!     MXS  = Max. allowed no. of scans in relaxation procedure.
  !!   RMASK  = Work array
  !!   ERROR  = Array containing the errors in the relaxation procedure.
  !!   NVALUE = no. of gridpoints with value (possibly 0) ... output
  !!-----------------------------------------------------------------------
  !!
  implicit none
  !c
  INTEGER NX,NY,I1,I2,J1,J2,MXS,NVALUE
  REAL    Tx,CRITX,COR
  !c
  REAL    ZA(NX,NY)
  REAL    RMASK(NX,NY)
  real    ERROR(NX,NY)
  integer K,N,J,I,I1P1,I2M1,J1P1,J2M1,NNN,NBAD
  real    SUMA,ASUMA,CRIT,CRTEST,UNDEF
  !C
  N=0
  SUMA=0.
  !C
  DO J=J1,J2
     DO I=I1,I2
        IF(ZA(I,J).LT.Tx) THEN
           SUMA=SUMA+ZA(I,J)
           N=N+1
        END IF
     END DO
  END DO

  NVALUE=N
  IF(N.LT.1) THEN
     !CC       WRITE(6,*)
     !CC       WRITE(6,*)'******************  WARNING  *******************'
     !CC       WRITE(6,*)'SUBROUTINE FILL : NO USEFUL DATA IN THE FIELD'
     !CC       WRITE(6,*)'ALL DATA IN THE INTERIOR (SUBAREA) WERE > ',Tx
     !CC       WRITE(6,*)'******************  WARNING  *******************'
     !CC       WRITE(6,*)
     RETURN
  END IF
  !C
  SUMA=SUMA/N
  ASUMA=0.
  !c
  DO J=J1,J2
     DO I=I1,I2
        IF(ZA(I,J).LT.Tx) THEN
           ASUMA=ASUMA+ABS(ZA(I,J)-SUMA)
           RMASK(I,J)=0.
        ELSE
           ZA(I,J)=SUMA
           RMASK(I,J)=1.
        END IF
     END DO
  END DO

  ASUMA=ASUMA/N
  !!Caf  Using ASUMA not SUMA
  CRIT=CRITX*ASUMA
  !C
  !C The value of SUMA ( i.e. the MEAN value of the array A ) is filled in
  !C the array A at all points with an initial value GREATER than or eq. to
  !C T. ZA(I,J) = SUMA may be regarded as the "first guess" in the iterative
  !C method.
  !C
  I1P1=I1+1
  I2M1=I2-1
  J1P1=J1+1
  J2M1=J2-1
  !C
  DO J=J1P1,J2M1
     DO I=I1P1,I2M1
        RMASK(I,J)=COR*RMASK(I,J)
     END DO
  END DO
  !C
  DO 200 NNN=1,MXS
     !C
     DO J=J1P1,J2M1
        !C  This is a fast vector loop, the next a slow scalar loop
        DO I=I1P1,I2M1
           ERROR(I,J)=(ZA(I+1,J)+ZA(I,J+1)+ZA(I,J-1))*0.25-ZA(I,J)
        END DO
        DO I=I1P1,I2M1
           !CCC  see loop above
           !CCC            ERROR(I,J)=(ZA(I+1,J)+ZA(I-1,J)+ZA(I,J+1)+ZA(I,J-1))
           !CCC                       *0.25 - ZA(I,J)
           ERROR(I,J)=ERROR(I,J)+ZA(I-1,J)*0.25
           ZA(I,J)=ZA(I,J)+ERROR(I,J)*RMASK(I,J)
        END DO
     END DO
     !C
     !C Test convergence now and then (slow test loop)
     IF(NNN.LT.MXS-5 .AND. MOD(NNN,10).EQ.0) THEN
        CRTEST=CRIT*COR
        NBAD=0
        J=J1
        DO WHILE (NBAD.EQ.0 .AND. J.LT.J2M1)
           J=J+1
           DO I=I1P1,I2M1
              IF(ABS(ERROR(I,J))*RMASK(I,J).GT.CRTEST) NBAD=1
           END DO
        END DO
        IF(NBAD.EQ.0) GOTO 300
     END IF
     !C
     DO J=J1P1,J2M1
        ZA(I1,J)=ZA(I1,J)+(ZA(I1P1,J)-ZA(I1,J))*RMASK(I1,J)
        ZA(I2,J)=ZA(I2,J)+(ZA(I2M1,J)-ZA(I2,J))*RMASK(I2,J)
     END DO
     DO I=I1,I2
        ZA(I,J1)=ZA(I,J1)+(ZA(I,J1P1)-ZA(I,J1))*RMASK(I,J1)
        ZA(I,J2)=ZA(I,J2)+(ZA(I,J2M1)-ZA(I,J2))*RMASK(I,J2)
     END DO
     !C
200  CONTINUE
     !C
300  CONTINUE
     !C
     RETURN
 
end subroutine fill
