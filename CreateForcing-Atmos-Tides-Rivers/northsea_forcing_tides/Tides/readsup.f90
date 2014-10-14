  subroutine readsup(yy_start, mm_start, dd_start, hh_start, mi_start,    &
                     yy_end, mm_end, dd_end, hh_end, mi_end)
! -----------------

!  logical, parameter, :: DEBUG = .TRUE.
  integer, parameter  :: unit = 17
  character(len=20) :: supfile 
!!$  logical,:: EOF = .FALSE.
  integer :: yy_start, mm_start, dd_start, hh_start, mi_start
  integer ::  yy_end, mm_end, dd_end, hh_end, mi_end
  character(len=20) :: line
!  character(len=80)::line
!
! --------------------

supfile="tidenc2roms.sup"
!!$  if (DEBUG) write(*,*) 'readsup: starting'
!
  open(unit, file=supfile, status='old', action='read')
!  if (DEBUG) write(*,*) 'readsup: setup-file is opened'
!


  call readln(line,unit)
  read(line,*)yy_start

  
  call readln(line,unit)
  read(line,*)mm_start

  call readln(line,unit)
  read(line,*)dd_start

  call readln(line,unit)
  read(line,*)hh_start 

  call readln(line,unit)
  read(line,*)mi_start 
  
  call readln(line,unit)
  read(line,*)yy_end

  call readln(line,unit)
  read(line,*)mm_end

  call readln(line,unit)
  read(line,*)dd_end

  call readln(line,unit)
  read(line,*)hh_end

  call readln(line,unit)
  read(line,*)mi_end
!!$!
!!$! --- Clean up and finish.
!!$!
    close(unit=unit)
!!$    if (DEBUG) write(*,*) "readsup: finished"

  end subroutine readsup
