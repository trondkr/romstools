  subroutine readln(Line,unit)
! -----------------------------
  !
  ! --------------------------------------
  !  Reads a line from a file, skipping
  !  comments and blank lines.
  !
  !  Comments starts with a character from
  !  COMCHAR and continues to the end of line.
  ! 
  !  Readln reads a line.
  !  If the line contains a comment, the comment is removed.
  !  If thereafter the line is blank, the line is skipped
  !  and the procedure repeated with a new line.
  !
  !  Bjoern aadlandsvik,
  !  IMR, October 1997
  ! --------------------------------------

  !
  ! --- Arguments ---
   character(len=20) :: Line
    integer :: lline
    integer :: unit
  !   First non comment line read
  !
  ! --- Local constants
  !
    character(len=*), parameter :: COMCHAR = "*!#"
  !   Comment starting characters
  !
  ! --- Local variables
  !
    integer :: ios
    integer :: ipos
  !   Start position for comment
  !
  ! --------------------------------
  !
    do
    !
    ! --- Read a line
    !
      read(unit=unit, fmt="(A)", iostat=ios) Line

      if (ios /= 0) then
!!!        EOF = .TRUE.
        return
      end if
    !
    ! --- Scan for comments
    !
      ipos = scan(Line, COMCHAR)
    !
    ! --- Remove any comments 
    !
      if (ipos /= 0) then  
        Line = Line(:ipos-1)
      end if
    !
    ! --- Skip the result if it is blank
    !
      if (len_trim(Line) == 0) then  
        cycle
      end if

      exit
    end do

    return

  end subroutine readln
! ------------------------------
! ******************************



