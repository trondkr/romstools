SUBROUTINE getflt(mode,iunit,itime,ident,ldata,fdata,igtype,gparam,nfskip,ierr)

!// 16 bit input and unpack (not handling 'undefined')

!// input:
!//    mode:   0 = read field
!//            1 = read field, skip fields until time > itime
!//            2 = read field if field time = itime
!//                (otherwise next read starts at the same field)
!//          100 = read field identification
!//                (next read starts at the same field)
!//          101 = read field identification, skip fields
!//                until time > itime
!//                (next read starts at the same (last) field)
!//          102 = read field identification, skip fields
!//                until time >= itime
!//                (next read starts at the same (last) field)
!//          200 = scan rest of the file and read field with
!//                matching identification, specified identification
!//                input in ident(1:20) where -32767 means any value
!//          201 = scan the whole file and read field with
!//                matching identification, specified identification
!//                input in ident(1:20) where -32767 means any value
!//           -1 = clean up after a file is closed, and the same
!//                file unit no. is used for another file.
!//    iunit: file unit no.
!//    itime: used if mode is 1, 2, 101 or 102.
!//           itime(1) - year
!//                (2) - month
!//                (3) - day
!//                (4) - hour (utc)
!//                (5) - forecast time in hours
!//                      (added to previous date and time)
!//    ident:  used if mode is 200 or 201
!//    ldata:  length of fdata (max field size)

!// output:
!//    ident(20): field identification
!//    fdata(..): field (unscaled, according to identification)
!//    igtype:    grid type (1=polar. 2=geo. 3+rot.sph. ...)
!//    gparam(6): grid parameters
!//    nfskip:    no. of fields skipped (for mode=2,101,102)
!//    ierr = 0:  read o.k.
!//           1:  read error
!//           2:  read error, end_of_file
!//           3:  no field returned due to mode=2 and itime spec.

!// note: mode=100,101,102 and possibly mode=2 will only read
!//       a field's identification and not the field (data).
!//       to avoid (cray) backspace problems, the identification
!//       is stored until the next read.

!// warning: using file unit no. (not file name) to identify
!//          files when storing field identification.
!//          If more than one file is opened with the same
!//          unit, use 'call getflt(-1,...)' after closing
!//          a file to avoid errors.

!// no computer dependant i/o methodes
!//                 ------
!// Does also work with CRAY's from PrgEnv 3.0.2.1
!// and above, provided an assign -F f77 -N ieee u:iunit
!// is applied somewhere.
!// NTNU/ITEA 1998-06-23 Jorn Amundsen
!//                 ------

!// DNMI/FoU  12.04.1992  Anstein Foss
!// DNMI/FoU  19.08.1993  Anstein Foss
!// DNMI/FoU  16.04.1998  Anstein Foss ... cray.t3e
!// DNMI/FoU  02.07.1998  Anstein Foss ... integer*2 also for Cray

IMPLICIT NONE

INTEGER, PARAMETER :: mfsize=1200000  ! Maximum field size, input and output

!// input/output

INTEGER            :: mode
INTEGER            :: iunit
INTEGER            :: ldata
INTEGER            :: igtype
INTEGER            :: nfskip         
INTEGER            :: ierr
INTEGER            :: itime(5)
INTEGER            :: ident(20)
REAL               :: fdata(ldata)
REAL               :: gparam(6)

!// local

INTEGER*2          :: ident2(20),idata2(mfsize)

INTEGER            :: itimef(5),idents(20)

INTEGER, PARAMETER :: maxsav=5, maxs20=maxsav*20
INTEGER            :: iusave(maxsav),idsave(20,maxsav)

INTEGER, PARAMETER :: mlgeom=50
INTEGER            :: lgtype
INTEGER*2          :: idgeom(20+1+mlgeom)
REAL               :: gplast(6)

INTEGER            :: i,nsave,isave,ios,iexit,iskip,n,ihours,ierr1,ierr2
INTEGER            :: nxin,nyin,nword,lgeom,nw,iscale,newgrid,ix,iy,ierror
REAL               :: scale

DATA iusave/maxsav*0/
DATA idsave/maxs20*0/

DATA lgtype/0/
DATA idgeom/71*0/
DATA gplast/6*999999./

IF (iunit < 1 .OR. iunit > 99 .OR.&
   (mode /= 0   .AND. mode /= 1   .AND. mode /= 2   .AND.&
    mode /= 100 .AND. mode /= 101 .AND. mode /= 102 .AND.&
    mode /= 200 .AND. mode /= 201 .AND. mode /= -1)) THEN
  PRINT *," **getflt** error ** mode,iunit: ",mode,iunit
  STOP
ENDIF

IF (mode == 200 .OR. mode.eq.201) THEN
  DO i=1,20
    idents(i)=ident(i)
  ENDDO
ENDIF

IF (mode == 201) REWIND(iunit)

!//      a) read field identification (one record)
!//      b) read field data           (one record)

ierr = 0
nfskip = 0

nsave = 0
isave = 0
DO i=1,maxsav
  IF (iusave(i) > 0)      nsave = nsave+1
  IF (iusave(i) == iunit) isave = i
ENDDO

IF (mode == -1) THEN
  IF (isave > 0) iusave(isave) = 0
  GOTO 990
ENDIF

IF (mode == 201 .AND. isave > 0) THEN
  iusave(isave) = 0
  isave = 0
  nsave = nsave-1
ENDIF

IF (isave > 0) THEN
  DO i=1,20
    ident(i) = idsave(i,isave)
  ENDDO
  iusave(isave) = 0
  nsave = nsave-1
  GOTO 30
ENDIF

20 CONTINUE

READ(iunit,IOSTAT=ios,ERR=910,END=920) (ident2(i),i=1,20)
DO i=1,20
  ident(i) = ident2(i)
ENDDO
!IF (ident(9) == 1010) ident(9) = 1

30 CONTINUE

!// ..iexit set to 1 if the next field is not read
!// ..iskip set to 1 if the next field is skipped
iexit = 0
iskip = 0

IF (mode == 100) THEN
  iexit = 1
ELSEIF (mode == 200 .OR. mode == 201) THEN
  n = 0
  DO i=1,20
    IF (idents(i) == -32767 .OR. idents(i) == ident(i)) n = n+1
  ENDDO
  IF (n /= 20) iskip = 1
ELSEIF (mode /= 0) THEN
  itimef(1) = ident(12)
  itimef(2) = ident(13)/100
  itimef(3) = ident(13)-(ident(13)/100)*100
  itimef(4) = ident(14)/100
  itimef(5) = ident(4)
  CALL hrdiff(0,0,itime,itimef,ihours,ierr1,ierr2)
  IF (mode == 1   .AND. ihours <=  0) iskip = 1
  IF (mode == 2   .AND. ihours /= 0)  iexit = 1
  IF (mode == 101 .AND. ihours <=  0) iskip = 1
  IF (mode == 101 .AND. ihours > 0)   iexit = 1
  IF (mode == 102 .AND. ihours < 0)   iskip = 1
  IF (mode == 102 .AND. ihours >=  0) iexit = 1
ENDIF

IF (iexit == 1) THEN
  !// Return without reading the field (data). save identification.
  IF (nsave == maxsav) iusave(1) = 0
  isave = 0
  DO n=1,maxsav
    IF (iusave(n) > 0) THEN
      isave = isave+1
      iusave(isave) = iusave(n)
      DO i=1,20
        idsave(i,isave) = idsave(i,n)
      ENDDO
    ENDIF
  ENDDO
  isave = maxsav
  iusave(isave) = iunit
  DO i=1,20
    idsave(i,isave) = ident(i)
  ENDDO
ENDIF

IF (iexit == 1) THEN
  ierr = 0
  IF (mode == 2) ierr=3
  GOTO 990
ENDIF

nxin = ident(10)
nyin = ident(11)
nword = nxin*nyin

lgeom = 0
IF (ident(9) >=  1000) lgeom=ident(9)-(ident(9)/1000)*1000

IF (nword > ldata .AND. iskip == 0) THEN
  PRINT *," **getflt** field length too big, (input ldata too small)"
  PRINT *," **ldata = ",ldata
  PRINT *," **ident: ",(ident(i),i=1,11)
  PRINT *," **       ",(ident(i),i=12,20)
  PRINT *," **nx,ny,nx*ny: ",nxin,nyin,nword
  ierr = 1
  GOTO 990
ENDIF

IF (nword+lgeom > mfsize) THEN
  PRINT *," **getflt** field length too big, (buffer too small)"
  PRINT *," **mfsize = ",mfsize
  PRINT *," **ident: ",(ident(i),i=1,11)
  PRINT *," **       ",(ident(i),i=12,20)
  PRINT *," **nx,ny,nx*ny,lgeom: ",nxin,nyin,nword,lgeom
  ierr = 1
  GOTO 990
ENDIF

nw = nword+lgeom
READ(iunit,IOSTAT=ios,ERR=910,END=920) (idata2(i),i=1,nw)       

IF (iskip == 0) THEN
  iscale = ident(20)
  scale = 10.**iscale
  !// Following (fast) loop does not handle 'undefined' values
  DO i=1,nword
    fdata(i) = scale*idata2(i)
  ENDDO
!//.............................................................
  !// The next (slow) loop handles 'undefined' values
!  DO i=1,nword
!    IF (idata2(i) /= -32767) THEN
!      fdata(i) = scale*idata2(i)
!    ELSE
!      fdata(i)=undef
!    ENDIF
!  ENDDO
ENDIF

IF (iskip /= 0) THEN
  nfskip = nfskip+1
  GOTO 20
ENDIF

IF (ident(9) > 0) THEN
  newgrid = 0
  IF (ident(9) /= idgeom(9)) THEN
    newgrid = 1
  ELSE
    DO i=15,18
      IF (ident(i) /= idgeom(i)) newgrid = 1
    ENDDO
    DO i=1,lgeom
      IF (idata2(nword+i) /= idgeom(21+i)) newgrid = 1
    ENDDO
  ENDIF
  IF (newgrid == 1) THEN
    DO i=1,20
      idgeom(i) = ident(i)
    ENDDO
    idgeom(10) = 1
    idgeom(11) = 1
    idgeom(21) = 0
    DO i=1,lgeom
      idgeom(21+i) = idata2(nword+i)
    ENDDO
    CALL gridpars(+1,20+1+mlgeom,idgeom,lgtype,ix,iy,gplast,ierror)
    IF (ierror /= 0) THEN
      PRINT *,"GETFLT: GRIDPARS ERROR ",ierror
      idgeom(9) = -9999
      ierr = 2
      GOTO 990
    ENDIF
  ENDIF
  igtype = lgtype
  DO i=1,6
    gparam(i) = gplast(i)
  ENDDO
ELSE
  igtype    = ident(9)
  gparam(1) = ident(15)
  gparam(2) = ident(16)
  gparam(3) = ident(17)
  gparam(4) = ident(18)
  gparam(5) = 0.
  gparam(6) = 0.
ENDIF

ierr = 0
GOTO 990

910 CONTINUE
  ierr = 1
  PRINT *," **getflt** read error. file,iostat: ",iunit,ios
  IF (mode == 200 .OR. mode == 201) GOTO 940
  GOTO 990

920 CONTINUE
  ierr = 2
  PRINT *," **getflt** end_of_file. file,iostat: ",iunit,ios
  IF (mode == 200 .OR. mode == 201) GOTO 940
  GOTO 990

940 CONTINUE
  DO i=1,20
    ident(i) = idents(i)
  ENDDO

990 CONTINUE

RETURN

END


!*******************************************************************

SUBROUTINE putflt(iunit,ident,igtype,gparam,ldata,fdata,scale,ierr)

!// 16 bit pack and output (handling 'undefined')

!// input:
!//    iunit: file unit
!//    ident: field identification, (set ident(20)=-32767 for automatic scaling)
!//    igtype: grid type
!//    gparam: grid specification parameters
!//    ldata: length of the field
!//           (field size given in identification used when writing the field)
!//    fdata: the field, unscaled
!//    scale: additional scaling of the field before output, usually 1
!//           (this is in addition to the scaling given in the identification, 
!//            for changing the basic unit)

!// output:
!//    ierr = 0:  write o.k.
!//           1:  write error

!// no computer dependant i/o methodes.
!//                 ------
!// Does also work with CRAY's from PrgEnv 3.0.2.1
!// and above, provided an assign -F f77 -N ieee u:iunit
!// is applied somewhere.
!// NTNU/ITEA 1998-06-23 Jorn Amundsen
!//                 ------

!// DNMI/FoU  12.04.1992  Anstein Foss
!// DNMI/FoU  19.08.1993  Anstein Foss
!// DNMI/FoU  10.10.1997  Anstein Foss ... added extra geometry ident.
!// DNMI/FoU  16.04.1998  Anstein Foss ... cray.t3e
!// DNMI/FoU  02.07.1998  Anstein Foss ... integer*2 also for Cray

IMPLICIT NONE

INTEGER, PARAMETER :: mfsize=1200000  ! Maximum field size, input and output
REAL,    PARAMETER :: undef=+1.e+35
REAL,    PARAMETER :: udef=undef*0.9

INTEGER, PARAMETER :: mldata=mfsize+50

!// input/output
INTEGER            :: iunit
INTEGER            :: igtype
INTEGER            :: ident(20)
INTEGER            :: ldata
INTEGER            :: ierr
REAL               :: gparam(6)
REAL               :: fdata(ldata)
REAL               :: scale

!// local
INTEGER            :: idento(20)
INTEGER*2          :: ident2(20),idata2(mldata)

INTEGER, PARAMETER :: mlgeom=50
INTEGER            :: lgtype,llgeom
INTEGER*2          :: idgeom(20+1+mlgeom)
REAL               :: gplast(6)

INTEGER            :: i,newgrid,ix,iy,ierror,lgeom,nxout,nyout,nword
INTEGER            :: sc,iscale,ifmax,ios
REAL               :: fmax

DATA lgtype,llgeom/-9999,0/
DATA idgeom/71*0/
DATA gplast/6*999999./

print *, fdata(138672)
DO i=1,20
  idento(i) = ident(i)
ENDDO

IF (igtype > 0) THEN
  newgrid = 0
  IF (igtype /= lgtype) newgrid = 1
  DO i=1,6
    IF (gparam(i) /= gplast(i)) newgrid = 1
  ENDDO
  IF (newgrid == 1) THEN
    lgtype = igtype
    DO i=1,6
      gplast(i) = gparam(i)
    ENDDO
    ix = 1
    iy = 1
    CALL gridpars(-1,20+1+mlgeom,idgeom,igtype,ix,iy,gparam,ierror)
    IF (ierror /= 0) THEN
      PRINT *,"PUTFLT: GRIDPARS ERROR ",ierror
      lgtype = -999
      ierr = 2
      GOTO 990
    ENDIF
    llgeom = 0
    IF (idgeom(9) >=  1000) llgeom = idgeom(9)-(idgeom(9)/1000)*1000
  ENDIF
  idento(9)  = idgeom( 9)
  idento(15) = idgeom(15)
  idento(16) = idgeom(16)
  idento(17) = idgeom(17)
  idento(18) = idgeom(18)
  lgeom = llgeom
ELSE
  !// in case of "something" stored as if it were a field, "geometry"
  !// identification is left unchanged (extended "geometry" not possible)
  lgeom = 0
ENDIF

nxout = ident(10)
nyout = ident(11)
nword = nxout*nyout

ierr = 0

IF (nword+lgeom > mldata) THEN
  PRINT *," **putflt** field length too big, (buffer too small)"
  PRINT *," **mfsize = ",mfsize
  PRINT *," **mldata = ",mldata
  PRINT *," **ident: ",(idento(i),i=1,11)
  PRINT *," **       ",(idento(i),i=12,20)
  PRINT *," **nx,ny,nx*ny,lgeom: ",nxout,nyout,nword,lgeom
  ierr = 1
  GOTO 990
ENDIF

sc = scale
IF (sc == 0.) sc = 1.
iscale = ident(20)

IF (iscale == -32767) THEN
  !// ..automatic scaling
  fmax = 0.
  !// Fast code, not handling undefines values
!  DO i=1,nword
!    fmax = MAX(fmax,ABS(fdata(i)))
!  ENDDO
!//.............................................................
  !//  Slow code, handling undefines values
  DO i=1,nword
    IF (fdata(i) < udef) fmax = MAX(fmax,ABS(fdata(i)))
  ENDDO
  IF (fmax > 0.) THEN
    fmax = fmax*ABS(sc)
    iscale = LOG10(fmax)-4.
    ifmax = NINT(fmax*10.**(-iscale))
    IF (ifmax < 3278) THEN
      iscale = iscale-1
      ifmax = NINT(fmax*10.**(-iscale))
    ENDIF
    IF (ifmax > 32766) iscale = iscale+1
      iscale = MAX(iscale,-30)
  ELSE
    iscale = 0
  ENDIF
  idento(20) = iscale
ENDIF

!// Scale and 'pack' data
sc = sc*10.**(-iscale)

!// Fast code not handling undefines values
!DO i=1,nword
!  idata2(i)= NINT(sc*fdata(i))
!ENDDO
!//.............................................................
!// Slow code handling undefines values
DO i=1,nword
  IF (fdata(i) < udef) THEN
    idata2(i) = NINT(sc*fdata(i))
  ELSE
    idata2(i) = -32767
  ENDIF
ENDDO
!//.............................................................
!//Fixed scaling, check integer*2 value range
!IF (ident(20) /= -32767) THEN
!  DO i=1,nword
!    idata2(i) = MAX(MIN(idata2(i),+32767),-32768)
!  ENDDO
!ENDIF
!//.............................................................
IF (lgeom > 0) THEN
  DO i=1,lgeom
    idata2(nword+i) = idgeom(21+i)
  ENDDO
  nword = nword+lgeom
ENDIF

DO i=1,20
  ident2(i) = idento(i)
ENDDO

!// Write field identification and field data

WRITE(iunit,IOSTAT=ios,ERR=900) (ident2(i),i=1,20)
WRITE(iunit,IOSTAT=ios,ERR=900) (idata2(i),i=1,nword)

GOTO 990

900 CONTINUE
  ierr = 1
  PRINT *," **putflt** write error. file,iostat: ",iunit,ios

990 CONTINUE

RETURN

END


!*******************************************************************

SUBROUTINE hrdiff(iup1,iup2,itime1,itime2,ihours,ierr1,ierr2)

!// Calculate interval in hours between 'itime1' and 'itime2',
!// ihours = 'itime2' - 'itime1'  (positive,zero or negative)

!// itime: year,month,day,time in hours, forecast length in hours
!// (positive,zero or negative)

!// If iup1=1: 'itime1' is updated to give valid date,time
!// If iup2=1: 'itime2' is updated to give valid date,time
!//    (forecast length = 0)

!// ierr1,ierr2: 0 = o.k. itime1,itime2
!//              1 = not o.k. itime1,itime2

!// DNMI/FoU  xx.xx.1990  Anstein Foss

INTEGER :: itime1(5)
INTEGER :: itime2(5)

INTEGER :: mdays(12)
INTEGER :: it(5,2)

DATA mdays/31,28,31,30,31,30,31,31,30,31,30,31/

DO i=1,5
  it(i,1) = itime1(i)
  it(i,2) = itime2(i)
ENDDO

!// Put prog.time into year,month,day,time

CALL vtime(it(1,1),ierr1)
CALL vtime(it(1,2),ierr2)

!//  'it' now gives "veryfing" time of 'itime1' and 'itime2'

IF (iup1 == 1 .AND. ierr1 == 0) THEN
  DO i=1,5
    itime1(i) = it(i,1)
  ENDDO
ENDIF
IF (iup2 == 1 .AND. ierr2 == 0) THEN
  DO i=1,5
    itime2(i) = it(i,2)
  ENDDO
ENDIF

IF (ierr1 /= 0 .OR. ierr2 /= 0) THEN
  ihours = -32767
  RETURN
ENDIF

DO i=1,4
  IF (it(i,1) /= it(i,2)) GOTO 55
ENDDO

!// No time difference
ihours = 0
RETURN

55 CONTINUE

nt1 = 1
nt2 = 2
IF (it(i,1) > it(i,2)) THEN
  nt1 = 2
  nt2 = 1
ENDIF
nhh = 0

IF (it(1,nt1) == it(1,nt2)) GOTO 70
iy = it(1,nt1)
!// Remaining hours first year:
DO im=it(2,nt1),12
  md = mdays(im)
  IF (im == 2) THEN
    IF (iy/4*4 == iy) md = 29
    IF (iy/100*100 == iy .AND. iy/400*400 /= iy) md = 28
  ENDIF
  nhh = nhh+md*24
ENDDO
nhh = nhh-(it(3,nt1)-1)*24-it(4,nt1)
!// One year steps
DO iy=it(1,nt1)+1,it(1,nt2)-1
  nd = 365
  IF (iy/4*4 == iy) nd = 366
  IF (iy/100*100 == iy .AND. iy/400*400 /= iy) nd = 365
  nhh = nhh+nd*24
ENDDO
it(1,nt1) = it(1,nt2)
it(2,nt1) = 1
it(3,nt1) = 1
it(4,nt1) = 0

70 CONTINUE

IF (it(2,nt1) == it(2,nt2)) GOTO 80
!// Remaining hours first month
iy = it(1,nt1)
im = it(2,nt1)
md = mdays(im)
IF (im == 2) THEN
  IF (iy/4*4 == iy) md = 29
  IF (iy/100*100 == iy .AND. iy/400*400 /= iy) md = 28
ENDIF
nhh = nhh+(md-it(3,nt1)+1)*24-it(4,nt1)
!// One month steps
DO im=it(2,nt1)+1,it(2,nt2)-1
  md = mdays(im)
  IF (im == 2) THEN
    IF (iy/4*4 == iy) md = 29
    IF (iy/100*100 == iy .AND. iy/400*400 /= iy) md = 28
  ENDIF
  nhh = nhh+md*24
ENDDO
it(2,nt1) = it(2,nt2)
it(3,nt1) = 1
it(4,nt1) = 0

80 CONTINUE

IF (it(3,nt1) /= it(3,nt2)) THEN
  nhh = nhh+(it(3,nt2)-it(3,nt1))*24-it(4,nt1)
  it(3,nt1) = it(3,nt2)
  it(4,nt1) = 0
ENDIF
!// Hours last day
nhh = nhh+it(4,nt2)-it(4,nt1)
it(4,nt1) = it(4,nt2)

IF (nt1 == 1) THEN
  ihours = nhh
ELSE
  ihours = -nhh
ENDIF

RETURN

END


!*******************************************************************

SUBROUTINE vtime(itime,ierror)

!// 'itime' is updated to give "valid" date,time (with prog.time = 0)

!// input:  itime(5) - itime(1): year
!//  itime(2): month (1-12)
!//  itime(3): day (1-28/29/30/31)
!//  itime(4): time in hours (00-23)
!//  itime(5): time in hours of prognosis
!//            (negative, zero or positive)
!// output: itime(5) -  as above, itime(5) = 0
!//   ierror   -  0 = o.k. input date/time
!//               1 = not o.k. input date/time
!// ('itime' not changed)

!// DNMI/FoU  xx.xx.1992  Anstein Foss

INTEGER :: itime(5)

INTEGER :: mdays(12)

DATA mdays/31,28,31,30,31,30,31,31,30,31,30,31/

iy = itime(1)
im = itime(2)
id = itime(3)
ih = itime(4)

!// Test input time
ierror = 0
IF (im < 1 .OR. im > 12) THEN
  ierror = 1
ELSE
  md = mdays(im)
  IF (im == 2) THEN
    IF (iy/4*4 == iy) md = 29
    IF (iy/100*100 == iy .AND. iy/400*400 /= iy) md = 28
  ENDIF
  IF (id < 1 .OR. id > md) ierror = 1
ENDIF
IF (ih < 00 .OR. ih > 23) ierror = 1

IF (ierror /= 0) RETURN

ih = ih+itime(5)
IF (ih >=  0 .AND. ih <=  23) GOTO 50
IF (ih < 0) GOTO 30

nd = ih/24
ih = ih-24*nd
DO n=1,nd
  id = id+1
  IF (id > md) THEN
    im = im+1
    IF (im > 12) THEN
      iy = iy+1
      im = 1
      md = mdays(im)
    ELSEIF (im == 2) THEN
      md = mdays(im)
      IF (iy/4*4 == iy) md = 29
      IF (iy/100*100 == iy .AND. iy/400*400 /= iy) md = 28
    ELSE
      md = mdays(im)
    ENDIF
    id = 1
  ENDIF
ENDDO

GOTO 50

30 CONTINUE

nd = (-ih+23)/24
ih = ih+24*nd
DO n=1,nd
  id = id-1
  IF (id < 1) THEN
    im = im-1
    IF (im < 1) THEN
      iy = iy-1
      im = 12
      md = mdays(im)
    ELSEIF(im == 2) THEN
      md = mdays(im)
      IF (iy/4*4 == iy) md = 29
      IF (iy/100*100 == iy .AND. iy/400*400 /= iy) md = 28
    ELSE
      md = mdays(im)
    ENDIF
    id = md
  ENDIF
ENDDO

50 CONTINUE

itime(1) = iy
itime(2) = im
itime(3) = id
itime(4) = ih
itime(5) = 0

RETURN

END


!*******************************************************************

SUBROUTINE gridpars(icall,ldata,idata,igtype,nx,ny,grid,ierror)

!// NAME:
!//    gridpar

!// PURPOSE:
!//    Conversion between integer*2 field identification and
!//    variables used in programs. Note that some gridtypes also
!//    has extended (geometry) identification behind the field data.

!// SYNOPSIS:
!//    SUBROUTINE gridpar(icall,ldata,idata,igtype,nx,ny,grid,ierror)
!//    integer   icall,ldata,igtype,nx,ny,ierror
!//    integer*2 idata(ldata)
!//    realgrid(6)

!// INPUT:
!//    icall  - +1 : from idata to igtype,nx,ny,grid
!//             -1 : from igtype,nx,ny,grid to idata
!//    ldata  - length of the idata array

!// INPUT/OUTPUT:
!//    idata  - field identification, field data (not touched)
!//             and possibly extra geometry specifications
!//             (max 20 words in this version)
!//    igtype - grid type, 1 = polarstereographic grid (true at 60N)
!//                        2 = geographic
!//                        3 = spherical rotated grid
!//                        4 = polarstereographic grid
!//                        5 = mercator grid (unrotated)
!//                        * = unknown grid type (-32767 to +32 accepted)
!//             (note that 'codes' are used in the field identifcation
!//             when extra geometry identification is stored after
!//             the field data, the reason for the +32 limit)
!//    nx     - no. of gridpoints in x direction (1 - 32767)
!//    ny     - no. of gridpoints in y direction (1 - 32767)
!//    grid   - grid parameters

!//    igtype=1,4, polarstereographic grid:
!//      grid(1) = x position of north pole (xp)
!//      grid(2) = y position of north pole (yp)
!//      grid(3) = no. of grid units between North Pole and Equator
!//      grid(4) = grid rotation angle (degrees), positive east, negative west
!//      grid(5) = projection latitude (degrees), standard is 60 (60 deg. N)
!//      grid(6) = 0. (not used)

!//    igtype=2,3, geographic or spherical rotated grid:
!//      grid(1) = western boundary (degrees) (longitude for x=1)
!//      grid(2) = southern boundary (degrees) (latitude for y=1)
!//      grid(3) = longitude increment (degrees)
!//      grid(4) = latitude  increment (degrees)
!//      grid(5) = longitude position of rotated equator (degrees) 
!//                (0 if geographic grid)
!//      grid(6) = latitude position of rotated equator (degrees)
!//                (0 if geographic grid)

!//    igtype=5, mercator (unrotated) grid:
!//      grid(1) = western boundary (degrees) (longitude for x=1)
!//      grid(2) = southern boundary (degrees) (latitude for y=1)
!//      grid(3) = x (longitude) increment (km)
!//      grid(4) = y (latitude)  increment (km)
!//      grid(5) = reference (construction) latitude (degrees)
!//      grid(6) = 0. (not used)

!//    igtype=*, unknown grid type, only use grid type less than 1 if the
!//              grid parameters have no meaning:
!//      grid(1:6) : unknown grid parameters

!// OUTPUT:
!//    ierror - error status: 0 = no error
!//                           1 = bad value in input identification
!//                               or in input grid parameters etc.
!//                           2 = ldata too small
!//                           3 = unknown icall
!//                           4 = ldata too small for needed extra
!//                               geometry identification (icall=-1),
!//                               but the best possible identification
!//                               is done (see NOTES below)

!// DOCUMENTATION:
!//    Basic document on felt files:
!//   FILE STRUKTUR FOR "SANNTIDS" LAGRING AV GRID-DATA
!//  Forskningsavdeling DNMI, oktober 1982
!//    See also /usr/local/doc/felt.doc (at DNMI)

!// NOTES:
!//  - This routine maintain compability with old formats
!//    (see comments in the source code),
!//    new formats added when required for some reason.
!//  - Specyfing ldata too short to hold extra geometry identification
!//    (for icall=-1) will restrict this routine from making this even
!//    when it seems appropriate. In order to be compatible with
!//    old unextended (or less extended) formats you then may
!//    disregard the returned ierror=4 (in special cases).
!//  - Avoid calling this routine more often than necessary,
!//    i.e. only after reading the first field or only before output
!//    of the first field.
!//    In models with several calling routines you may do a 'setup' call
!//    and keep (part of) the first 20 words of identification and
!//    possibly extra geometry specification locally (see note below).
!//  - The value of nx and ny (in idata if icall=+1) does not influence
!//    conversion of grid parameters.
!//    You may then use this routine to convert between field
!//    identification and grid parameters with nx=ny=1 and possibly
!//    find extra identification in word 22,23,... in idata.
!//    (see DOCUMENTATION to find format description).

!// DNMI/FoU  05.05.1995  Anstein Foss
!// DNMI/FoU  09.06.1995  Anstein Foss
!// DNMI/FoU  14.05.1996  Anstein Foss ... mercator (unrotated)
!// DNMI/FoU  02.09.1996  Anstein Foss ... gridtype 2012 -> 2 (bad test)
!// DNMI/FoU  15.10.1996  Anstein Foss ... even better scaling when needed
!// DNMI/FoU  17.02.1997  Anstein Foss ... and again (for 'image' fields)

IMPLICIT NONE

!// Input/output
INTEGER   :: icall
INTEGER   :: ldata
INTEGER   :: igtype
INTEGER   :: nx,ny
INTEGER   :: ierror
INTEGER*2 :: idata(ldata)
REAL      :: grid(6)

!// Local
INTEGER   :: ngw,lgeom,i,igr,ig1,ig2,ld,igscale,kgeom,lgeom1,lgeom2
INTEGER*2 :: igeom1(12),igeom2(20)
REAL      :: gscale,glog,gws,dglim,dgmax,dgprv,grx
REAL      :: gw(6),gr(6)

ierror = 0

IF (icall == +1) THEN   ! idata(ldata) -> igtype,nx,ny,grid(6)

  igtype = 0
  nx = 0
  ny = 0
  DO i=1,6
    grid(i) = 0.
  ENDDO
  IF (ldata < 20) THEN
    ierror = 2
    RETURN
  ENDIF
  igtype = idata(9)
  lgeom = 0
  IF (igtype > 999) THEN
    i = igtype
    igtype = igtype/1000
    lgeom = i-igtype*1000
  ENDIF
  nx = idata(10)
  ny = idata(11)
  IF (nx < 1 .OR. ny < 1) THEN
    ierror = 1
    RETURN
  ENDIF

  IF (igtype == 1 .OR. igtype == 4) THEN

    !// Polarstereographic (type 1 always true at 60 degrees North)

    IF (idata(17) > 0) THEN
      !// The standard
      grid(1) = idata(15)*0.01
      grid(2) = idata(16)*0.01
      grid(3) = idata(17)*0.1
      grid(4) = idata(18)
    ELSE
      !// An old extension
      grid(1) = idata(15)
      grid(2) = idata(16)
      grid(3) = -idata(17)*0.1
      grid(4) = idata(18)
    ENDIF
    grid(5) = 60.

  ELSEIF (igtype == 2 .OR. igtype == 3) THEN

    !// Geographic (2) or spherical rotated grid (3)

    grid(1) = idata(16)*0.01
    grid(2) = idata(15)*0.01
    grid(3) = idata(18)*0.01
    grid(4) = idata(17)*0.01

  ELSEIF (igtype == 5) THEN

    !// Mercator grid (5)

    grid(1) = idata(15)*0.01
    grid(2) = idata(16)*0.01
    grid(3) = idata(17)*0.1
    grid(4) = idata(18)*0.1

  ELSE

    !// Unknown/undefined grid type

    grid(1) = idata(15)
    grid(2) = idata(16)
    grid(3) = idata(17)
    grid(4) = idata(18)

  ENDIF

  IF (lgeom > 0) THEN
    IF (igtype == 1 .OR. igtype == 4) THEN
      gscale = 100.
      ngw = 5
    ELSEIF (igtype == 2 .OR. igtype == 3) THEN
      gscale = 10000.
      ngw = 6
    ELSEIF (igtype == 5) THEN
      gscale = 10000.
      ngw = 5
    ELSE
      gscale = 100.
      ngw = 6
    ENDIF

    ld=20+nx*ny

    IF (lgeom == ngw*2 .AND. ld+lgeom <=  ldata) THEN

      !// First extended method
      DO i=1,ngw
        ig1 = idata(ld+1)
        ig2 = idata(ld+2)
        ld = ld+2
        grid(i) = FLOAT(ig1*10000+ig2)/gscale
      ENDDO

    ELSEIF (lgeom == 2+ngw*3 .AND. ld+lgeom <=  ldata) THEN

      !// Second extended method
      IF (idata(ld+1) == ngw .AND. idata(ld+2) == 3) THEN
        ld=ld+2
        DO i=1,ngw
          igscale = idata(ld+1)
          ig1 = idata(ld+2)
          ig2 = idata(ld+3)
          ld = ld+3
          gscale = 10.**igscale
          grid(i) = FLOAT(ig1*10000+ig2)*gscale
        ENDDO
      ELSE
        ierror = 2
      ENDIF

    ELSE

      ierror = 2

    ENDIF

  ENDIF

  IF (ierror == 0) THEN
    IF (igtype == 1 .OR. igtype == 4) THEN
      !// The DNMI "standard" (150km grid => an=grid(3)=79.)
      IF (grid(3) /= 0.) grid(3) = 79.*150./grid(3)
      IF (grid(3) == 0.  .OR. grid(5) == 0. .OR.&
          grid(5) < -90. .OR. grid(5) > +90.) ierror = 1
    ELSEIF (igtype == 2 .OR. igtype == 3) THEN
      IF (grid(3) == 0. .OR. grid(4) == 0.) ierror = 1
    ELSEIF (igtype == 5) THEN
      IF (grid(3) == 0. .OR. grid(4) == 0.) ierror = 1
    ENDIF
  ENDIF

ELSEIF (icall == -1) THEN   ! igtype,nx,ny,grid(6) -> idata(ldata)

  IF (ldata < 20) THEN
    ierror = 2
    RETURN
  ENDIF

  IF (igtype > 32 .OR. nx < 1 .OR. nx > 32767 .OR. ny < 1 .OR. ny > 32767) THEN
    ierror = 1
    RETURN
  ENDIF

  idata(9) = igtype
  idata(10) = nx
  idata(11) = ny
  DO i=1,6
    gw(i) = grid(i)
  ENDDO

  IF (igtype == 1 .OR. igtype == 4) THEN

    !// Polarstereographic (type 1 always true at 60 degrees North)

    IF (gw(3) == 0.) THEN
      ierror = 1
      RETURN
    ENDIF

    !// The DNMI "standard" (150km grid => an=grid(3)=79.)
    gw(3) = 150.*79./gw(3)

    IF (ABS(gw(1)) < 327.66 .AND. ABS(gw(2)) < 327.66) THEN
      !// The standard
      idata(15) = NINT(gw(1)*100.)
      idata(16) = NINT(gw(2)*100.)
      idata(17) = NINT(gw(3)*10.)
      idata(18) = NINT(gw(4))
      gr(1) = FLOAT(idata(15))*0.01
      gr(2) = FLOAT(idata(16))*0.01
      gr(3) = FLOAT(idata(17))*0.1
      gr(4) = FLOAT(idata(18))
      gr(5) = 60.
    ELSEIF(ABS(gw(1)) < 32766. .AND. ABS(gw(2)) < 32766.) THEN
      !// An old extension
      idata(15) = NINT(gw(1))
      idata(16) = NINT(gw(2))
      idata(17) = -NINT(gw(3)*10.)
      idata(18) = NINT(gw(4))
      gr(1) = FLOAT(idata(15))
      gr(2) = FLOAT(idata(16))
      gr(3) = -FLOAT(idata(17))*0.1
      gr(4) = FLOAT(idata(18))
      gr(5) = 60.
    ELSE
      !// Old impossible case
      idata(15) = 32767
      idata(16) = 32767
      idata(17) = 32767
      idata(18) = 32767
      gr(1) = 327.67
      gr(2) = 327.67
      gr(3) = 3276.7
      gr(4) = 32767.
      gr(5) = 60.
    ENDIF

    gscale = 100.
    ngw = 5

  ELSEIF (igtype == 2 .OR. igtype == 3) THEN

    !// Geographic (2) or spherical rotated grid (3)

    IF (gw(3) == 0. .OR. gw(4) == 0.) THEN
      ierror = 1
      RETURN
    ENDIF

    IF (gw(1) > +180.) gw(1) = gw(1)-360.
    IF (gw(1) < -180.) gw(1) = gw(1)+360.

    idata(15) = NINT(gw(2)*100.)
    idata(16) = NINT(gw(1)*100.)
    idata(17) = NINT(gw(4)*100.)
    idata(18) = NINT(gw(3)*100.)
    gr(1) = FLOAT(idata(16))*0.01
    gr(2) = FLOAT(idata(15))*0.01
    gr(3) = FLOAT(idata(18))*0.01
    gr(4) = FLOAT(idata(17))*0.01
    gr(5) = 0.
    gr(6) = 0.

    gscale = 10000.
    ngw = 6

  ELSEIF (igtype == 5) THEN
 
    !// Mercator grid (5)

    IF (gw(3) == 0. .OR. gw(4) == 0.) THEN
      ierror = 1
      RETURN
    ENDIF

    idata(15) = NINT(gw(1)*100.)
    idata(16) = NINT(gw(2)*100.)
    idata(17) = NINT(gw(3)*10.)
    idata(18) = NINT(gw(4)*10.)
    gr(1) = FLOAT(idata(15))*0.01
    gr(2) = FLOAT(idata(16))*0.01
    gr(3) = FLOAT(idata(17))*0.1
    gr(4) = FLOAT(idata(18))*0.1
    gr(5) = 0.

    gscale = 10000.
    ngw = 5

  ELSE

    !// Unknown/undefined grid type

    idata(15) = NINT(gw(1))
    idata(16) = NINT(gw(2))
    idata(17) = NINT(gw(3))
    idata(18) = NINT(gw(4))
    gr(1) = FLOAT(idata(15))
    gr(2) = FLOAT(idata(16))
    gr(3) = FLOAT(idata(17))
    gr(4) = FLOAT(idata(18))
    gr(5) = 0.
    gr(6) = 0.

    gscale = 100.
    ngw = 6

  ENDIF

  !// Check if the standard packing above was good enough
  !// or if the first or second extended method should be used,
  !// for compability with old formats the least extended is preferred

  !// A limit to avoid high precision complications
  dglim = 1.e-8

  dgmax = 0.

  IF (igtype > 0) THEN
    DO i=1,ngw
      IF (gw(i) /= 0.) dgmax = MAX(dgmax,ABS((gr(i)-gw(i))/gw(i)))
    ENDDO
  ENDIF

  IF (dgmax > dglim) THEN

    kgeom = 0
    ld = 20+nx*ny

  !// First extended method, same scaling for all grid parameters
  !// but don't use it unless it gives better results
  !// (kept in the code for compability with old formats, possibly
  !//  in models etc. not using this routine)

    dgprv = dgmax
    dgmax = 0.
    lgeom1 = 0

    DO i=1,ngw
      gws = gw(i)*gscale
      !// Check overflow
      IF (gws /= 0. .AND. ABS(gws) < 3.e+8) THEN
        igr = NINT(gws)
        ig1 = igr/10000
        ig2 = igr-ig1*10000
        igeom1(lgeom1+1) = ig1
        igeom1(lgeom1+2) = ig2
        grx = FLOAT(ig1*10000+ig2)/gscale
        dgmax = MAX(dgmax,ABS((grx-gw(i))/gw(i)))
      ELSE
        !// Zero or value not possible for this method
        igeom1(lgeom1+1) = 0
        igeom1(lgeom1+2) = 0
        IF (gws /= 0.) dgmax = 1.
      ENDIF
      lgeom1 = lgeom1+2
    ENDDO

    IF (dgmax < dgprv) THEN
      IF (ldata >= ld+lgeom1) THEN
        kgeom = 1
      ELSE
        ierror = 4
      ENDIF
    ELSE
      dgmax = dgprv
    ENDIF

    IF (dgmax > dglim .AND. ierror == 0) THEN

      !// Second extended method, good enough for any real*4 precision value
      !// but don't use it unless it gives better results

      dgprv = dgmax
      dgmax = 0.

      igeom2(1) = ngw
      igeom2(2) = 3
      lgeom2 = 2
      DO i=1,ngw
        IF (gw(i) /= 0.) THEN
          !// 8 decimals precision (more than enough for real*4)
          glog = LOG10(ABS(gw(i)))-8.
          igscale = INT(glog)
          IF (glog > 0.) igscale = igscale+1
          !// Keep scaling within real*4 range
          igscale = MAX(MIN(igscale,+25),-25)
          gscale = 10.**(-igscale)
          igr = NINT(gw(i)*gscale)
          ig1 = igr/10000
          ig2 = igr-ig1*10000
          igeom2(lgeom2+1) = igscale
          igeom2(lgeom2+2) = ig1
          igeom2(lgeom2+3) = ig2
          grx = FLOAT(ig1*10000+ig2)/gscale
          dgmax = MAX(dgmax,ABS((grx-gw(i))/gw(i)))
        ELSE
          igeom2(lgeom2+1) = 0
          igeom2(lgeom2+2) = 0
          igeom2(lgeom2+3) = 0
        ENDIF
        lgeom2 = lgeom2+3
      ENDDO

      IF (dgmax < dgprv) THEN
        IF (ldata >= ld+lgeom2) THEN
          kgeom = 2
        ELSE
          ierror = 4
        ENDIF
!      ELSE
!        dgmax=dgprv
      ENDIF

    ENDIF

    IF (kgeom == 1) THEN
      idata(9) = igtype*1000+lgeom1
      DO i=1,lgeom1
        idata(ld+i) = igeom1(i)
      ENDDO
    ELSEIF (kgeom == 2) THEN
      idata(9) = igtype*1000+lgeom2
      DO i=1,lgeom2
        idata(ld+i) = igeom2(i)
      ENDDO
    ENDIF

  ENDIF


ELSE

  !// Wrong icall

  ierror = 3

ENDIF

RETURN

END


!*******************************************************************

SUBROUTINE setbilin(igtype0,gparam0,igtype,gparam,Lp,Mp,lon_rho,lat_rho, &
                    i1,j1,r1,r2,r3,r4)

!// set bi-linear interpolation factors

!// NOTES:
!//  below, igtype0 must be initialized to 0;
!//         igtype0, gparam0 must NOT be changed outside setbilin

!// input:
!//    igtype0: grid type from previous call (set this to 0 initially!)
!//    gparam0: grid specification parameters from previous call
!//    igtype:  present grid type
!//    gparam:  present grid specification parameters
!//    Lp, Mp:  x,y dimensional size of output grid
!//    lon_rho: longitude array for output grid
!//    lat_rho: latitude array for output grid

!// output:
!//    i1,j1:   position in input grid to the "lower left" of the
!//              corresponding output grid position
!//    r1,r2,r3,r4:
!//             interpolation weights (sum 1) for
!//              lower left, lower right, upper right, upper left, respectively

!// met.no/FoU  29.12.2008  Arne Melsom (from software by B. Ådlandsvik, IMR)


IMPLICIT NONE

!// input/output
INTEGER            :: igtype0,igtype
REAL               :: gparam0(6),gparam(6)
INTEGER            :: Lp,Mp
REAL               :: lon_rho(Lp,Mp),lat_rho(Lp,Mp)
REAL               :: r1(Lp,Mp),r2(Lp,Mp),r3(Lp,Mp),r4(Lp,Mp)
INTEGER            :: i1(Lp,Mp),j1(Lp,Mp)

!// local
INTEGER            :: i,j,n
REAL               :: pi,rad,vxr,alfa,beta,rr,xr,yr,rx,ry,x,y
LOGICAL            :: samegrid


! Same grid as in previous call to this subroutine?
samegrid = (igtype0 == igtype)
n = 0
DO WHILE (samegrid .AND. n < 6)
   n = n + 1
   samegrid = (gparam0(n) == gparam(n))
ENDDO
IF (.NOT. samegrid .AND. igtype0 /= 0) &
  STOP "Grid changed, possible conflicts wrt. input depths and/or vector components!"


! (Re-)Set interpolation weights?

IF (.NOT. samegrid) THEN

   ! Re-set grid spec. to the present spec.:
   igtype0 = igtype
   DO n = 1,6
      gparam0(n) = gparam(n)
   ENDDO

   ! Set constants

   pi = 4.*ATAN(1.)
   rad = pi/180.

   ! Set some interpolation factors for polar stereographic grids:
   IF (igtype == 1 .OR. igtype == 4) THEN
      vxr = (90.+gparam(4))*rad
      beta = SIN(vxr)
      alfa = COS(vxr)
   ENDIF

   ! Horizontal bilinear interpolation to ROMS grid
   DO j = 1,Mp
      DO i = 1,Lp

         !// Find gridpoint in input grid
         IF (igtype == 1 .OR. igtype == 4) THEN
            rr = gparam(3)*COS(lat_rho(i,j)*rad)/(1.+SIN(lat_rho(i,j)*rad))
            xr = +rr*SIN(lon_rho(i,j)*rad)
            yr = -rr*COS(lon_rho(i,j)*rad)
            x  = xr*beta - yr*alfa + gparam(1)
            y  = yr*beta + xr*alfa + gparam(2)
         ELSEIF (igtype == 2) THEN
            x = (lon_rho(i,j) - gparam(1))/gparam(3) + 1.
            y = (lat_rho(i,j) - gparam(2))/gparam(4) + 1.
         ELSE
            STOP "Program not adjusted to this grid type!"
         ENDIF
         i1(i,j) = ifix(x)
         j1(i,j) = ifix(y)

         ! Compute, store interpolation weights:
         rx = (x-ifix(x))
         ry = (y-ifix(y))
         r1(i,j) = (1.-rx)*(1.-ry)
         r2(i,j) =     rx *(1.-ry)
         r3(i,j) =     rx *    ry
         r4(i,j) = (1.-rx)*    ry

      ENDDO ! i
   ENDDO  ! j

ENDIF

RETURN

END SUBROUTINE setbilin


SUBROUTINE LLangle(igtypev,gparamv,igtype,gparam, &
                   Lp,Mp,lon_rho,lat_rho,uu,uv,vu,vv,ierror)

!// NOTES:
!//  below, igtypev must be initialized to 0;
!//         igtypev, gparamv must NOT be changed outside LLangle

!// input:
!//    igtypev: grid type from previous call (set this to 0 initially!)
!//    gparamv: grid specification parameters from previous call
!//    igtype:  present grid type
!//    gparam:  present grid specification parameters
!//    Lp, Mp:  x,y dimensional size of output grid



IMPLICIT NONE

!// input/output
INTEGER    :: igtypev,igtype,ierror
REAL       :: gparamv(6),gparam(6)
INTEGER    :: Lp,Mp
REAL       :: lon_rho(Lp,Mp),lat_rho(Lp,Mp)
REAL       :: uu(Lp,Mp),uv(Lp,Mp),vu(Lp,Mp),vv(Lp,Mp)

!// local
INTEGER    :: i,j,n,irot
LOGICAL    :: samegrid
REAL       :: zpir18,dxr,dyr,xcr,ycr,xca,yca,zsyca,zcyca,xwr,ysr
REAL       :: x2,y2,x3,y3
REAL       :: zsxsph,zcxsph,zsysph,zcysph,zxmxc,zsxmxc,zcxmxc
REAL       :: zsxrot,zcxrot,zsyrot,zcyrot,fia,zfia,zns,zsfi,zcfi


ierror = 0

! Same grid as in previous call to this subroutine?
samegrid = (igtypev == igtype)
IF (.NOT. samegrid .AND. igtypev /= 0) &
  STOP "Input grid changed, terminates due to possible conflicts wrt. input depths!"


! (Re-)Set interpolation weights?

IF (.NOT. samegrid) THEN

  ! Re-set grid spec. to the present spec.:
  igtypev = igtype
  DO n = 1,6
    gparamv(n) = gparam(n)
  ENDDO

  zpir18 = 2.0*asin(1.0)/180.

  dxr = zpir18
  dyr = zpir18
  xcr = 0.
  ycr = 0.

  IF (igtype == 2) THEN

    ! Nothing to do if grid is already spherical:
    uu(:,:) = 1.
    uv(:,:) = 0.
    vu(:,:) = 0.
    vv(:,:) = 1.

  ELSE IF (igtype == 3) THEN

    !..rot.sph->sph

    xca = gparam(5)*zpir18
    yca = gparam(6)*zpir18
    if (xca /= 0. .or. yca /= 0.) then
      irot=1
    else
      irot=0
    end if

    zsyca = sin(yca)
    zcyca = cos(yca)

    do j = 1, Mp
      do i = 1, Lp

        if (irot == 1) then

          xwr = lon_rho(i,j)*zpir18
          ysr = lat_rho(i,j)*zpir18

          x2 = xwr + dxr
          y2 = ysr + dyr

          x3 = x2
          y3 = y2
          call ROMSsph2rot(x3,y3,xca,yca)

          zsxsph = sin(x2)
          zcxsph = cos(x2)
          zsysph = sin(y2)
          zcysph = cos(y2)
          zxmxc  = x2 - xca
          zsxmxc = sin(zxmxc)
          zcxmxc = cos(zxmxc)
          zsxrot = sin(x3)
          zcxrot = cos(x3)
          zsyrot = sin(y3)
          zcyrot = cos(y3)
          uu(i,j) = zcxmxc*zcxrot + zcyca*zsxmxc*zsxrot
          uv(i,j) = zcyca*zsxmxc*zcxrot*zsyrot + zsyca*zsxmxc*zcyrot - &
                    zcxmxc*zsxrot*zsyrot
          vu(i,j) =-zsyca*zsxrot/zcysph
          vv(i,j) = (zcyca*zcyrot - zsyca*zcxrot*zsyrot)/zcysph

        else

          uu(i,j) = 1.
          uv(i,j) = 0.
          vu(i,j) = 0.
          vv(i,j) = 1.

        endif

      enddo ! i
    enddo ! j

  ELSEIF (igtype == 1 .or. igtype == 4) THEN

    !..pol->sph

    fia = gparam(4)*zpir18

    zfia = fia
    zns  = 1.
    if (gparam(5) < 0.) then
      zfia = -fia
      zns  = -1.
    end if

    zsfi = sin(zfia)
    zcfi = cos(zfia)

    do j = 1, Mp

      do i = 1, Lp

        xwr = lon_rho(i,j)*zpir18
        ysr = lat_rho(i,j)*zpir18

        x2 = xwr + dxr
        y2 = ysr + dyr

        zsxsph = sin(x2)
        zcxsph = cos(x2)
        zsysph = sin(y2)
        zcysph = cos(y2)
        uu(i,j) = zsxsph*zsfi + zcxsph*zcfi
        uv(i,j) = zsxsph*zcfi - zcxsph*zsfi
        vu(i,j) = zcxsph*zsfi - zsxsph*zcfi
        vv(i,j) = zcxsph*zcfi + zsxsph*zsfi

      enddo ! i
    enddo ! j

  ELSE

    write(*,*)'Input grid type in subroutine LLangle not recognized; type:', &
              igtype
    ierror = 1

  ENDIF

ENDIF


RETURN

END SUBROUTINE LLangle


subroutine ROMSsph2rot(x,y,xcen,ycen)
!
!  conversion between spherical (x=xsph,y=ysph) and spherical rotated
!  (x=xrot,y=yrot) coordinates. (xcen,ycen) is the position of the
!  rotated equator/greenwich in terms of (longitude,latitude).
!  all values are given in radians.
!
!----------------------------------------------------------------------
!   DNMI/FoU  xx.xx.1995  Jan Erik Haugen ... Hirlam code (?)
!   DNMI/FoU  xx.05.1995  Anstein Foss ...... DNMI library version
! met.no/FoU  22.06.2006  Anstein Foss ...... double precision computations
! met.no/FoU  02.01.2009  Arne Melsom  ...... single node "ROMS version"
!----------------------------------------------------------------------

implicit none

real    x, y, xcen, ycen

double precision zsycen,zcycen,xsph,ysph,zxmxc,zsxmxc,zcxmxc, &
                 zsysph,zcysph,zsyrot,yrot,                   &
                 zcyrot,zcxrot,zsxrot,xrot,zxcen,zycen

zxcen= xcen
zycen= ycen
zsycen = sin(zycen)
zcycen = cos(zycen)

!  compute spherical rotated coordinates as function of
!  spherical coordinates

xsph = x
ysph = y
zxmxc  = xsph - zxcen
zsxmxc = sin(zxmxc)
zcxmxc = cos(zxmxc)
zsysph = sin(ysph)
zcysph = cos(ysph)
zsyrot = zcycen*zsysph - zsycen*zcysph*zcxmxc
zsyrot = max(zsyrot,-1.0d0)
zsyrot = min(zsyrot,+1.0d0)
yrot   = asin(zsyrot)
zcyrot = cos(yrot)
zcxrot = (zcycen*zcysph*zcxmxc + zsycen*zsysph)/zcyrot
zcxrot = max(zcxrot,-1.0d0)
zcxrot = min(zcxrot,+1.0d0)
zsxrot = zcysph*zsxmxc/zcyrot
xrot   = acos(zcxrot)
if (zsxrot < 0.0d0) xrot = -xrot
x = xrot
y = yrot

return
end subroutine ROMSsph2rot
