!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! read coordinates
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
module qcxms_read_coordinates
  use readcommon
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_convert
  use xtb_mctc_symbols, only: toNumber
  implicit none
  
  contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initial read routine for defining number of atoms (nuc)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine rd0(fname,nuc)
!  
!    integer  :: nuc,nn
!    integer  :: io_coord, io_error
!  
!    real(wp) :: xx(10)
!  
!    character(len=:),allocatable :: fname
!    character(len=:),allocatable :: xyzfile
!
!    logical :: ex
!
!    inquire(file='coord',exist=ex)
!    if ( ex ) then
!      fname = 'coord'
!!      tmcoord = .true.
!    else
!      write(*,*) '- Which input file ? - '
!!      call get_command_line_argument(1, xyzfile)
!      if (xyzfile /= ' ') then 
!!        tmcoord = .false.
!        fname = xyzfile
!        write(*,*) ' Reading input file : ', fname
!      else
!        stop 'Input file is empty'
!      endif
!    endif
!
!    open(file=fname, newunit=io_coord, status='old', iostat=io_error)
!
!    if (io_error > 0) then
!      write (*,*)
!      stop ' - Missing coordinate file - E X I T - '
!      write (*,*)
!    endif
!
!    nuc = 0
!    do
!      read(io_coord,'(a)',iostat=iocheck)line
!      if (iocheck < 0) exit
!      if (iocheck > 0) stop 'Error in rd0'
!      if(index(line,'$redu').ne.0) exit 
!      if(index(line,'$user').ne.0) exit
!  
!      call readl(line,xx,nn)
!  
!      if ( nn < 3) cycle
!  
!      nuc = nuc + 1
!  
!    enddo
!    close(io_coord)
!
!  end subroutine rd0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read the coordinates in 'coord' style
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rd(init,fname,nuc,xyz,iat)
  
    intrinsic :: get_command_argument

    integer  :: nuc,nn, k
    integer,optional  :: iat(nuc)
    integer  :: io_coord, io_error
    integer :: check
    
    real(wp),optional :: xyz(3,nuc)
    real(wp) :: xx(10)
    
    character(len=:), allocatable :: fname
    character(len=30) :: xyzfile

    logical :: tmcoord
    logical :: init
    logical :: ex

    tmcoord = .false.

    ! this is not yet finished
!start: if ( init ) then
!      inquire(file='coord',exist=ex)
!
!exst: if ( ex ) then
!        fname = 'coord'
!
!      else ! get argument line input file or give file name
!
!        call get_command_argument(1, value=xyzfile, status=check)
!        fname = xyzfile
!
!        if (check > 0 ) then
!          write(*,*) '- Which input file ? - '
!          read(*,'(a)') xyzfile
!
!          if (xyzfile /= ' ') then 
!            tmcoord = .false.
!            fname = xyzfile
!            write(*,*) ' Reading input file : ', fname
!          else
!            stop '- No input provided - E X I T - '
!          endif
!
!        endif 
!
!      endif exst
!
!    endif start

!    if (fname == 'coord') tmcoord = .true.

    inquire(file='coord',exist=ex)
exst: if ( ex ) then
        fname = 'coord'

      else ! get argument line input file or give file name
        stop ' - No coord file available - E X I T - '
      endif exst

    open(file=fname, newunit=io_coord, status='old', iostat=io_error)

    if (io_error > 0) then
      write (*,*)
      stop ' - Missing coordinate file - E X I T - '
      write (*,*)
    endif

    
    nuc=0
  
    ! either turbomole coord file or xyz file

!    if (tmcoord) then
      do
        read(io_coord,'(a)',iostat=iocheck)line
        if (iocheck < 0) exit
        if (iocheck > 0) stop 'Error in rd'
      
        if(index(line,'$redu').ne.0)return
        if(index(line,'$user').ne.0)return
      
        call readl(line,xx,nn)

        if (nn < 3) cycle 

        nuc=nuc+1

        if (.not. init) then
          xyz(1,nuc)=xx(1)
          xyz(2,nuc)=xx(2)
          xyz(3,nuc)=xx(3)
          iat(nuc)=toNumber(line)
        endif

      enddo

!    else ! no coord file (here we take xyz)
!      read(io_coord,*) nuc
!      read(io_coord,*)
!
!      if (.not. init) then
!
!        nuc = 0
!        do
!          read(io_coord, '(a)',iostat=iocheck) line
!          if (iocheck > 0) stop 'error in xyz file'
!          if (iocheck < 0) exit ! EOF
!
!          nuc=nuc+1
!
!          call readl(line,xx,nn)
!          iat(nuc)=toNumber(line)
!          xyz(1,nuc)=xx(1)*aatoau
!          xyz(2,nuc)=xx(2)*aatoau
!          xyz(3,nuc)=xx(3)*aatoau
!
!        end do
!      endif
!    endif
    
    close(io_coord)

  end subroutine rd
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read qcxms.start file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rdstart(itrj,nat,xyz,velo,velof,tadd,eimp)
     use xtb_mctc_accuracy, only: wp
     implicit none
  
     integer  :: itrj,nat
     integer  :: j,ndum
     integer  :: io_start, ierror

     real(wp) :: xyz (3,nat)
     real(wp) :: velo(3,nat)
     real(wp) :: velof   (nat)
     real(wp) :: tadd,eimp
  
     character(len=80) :: fname
  
     fname='qcxms.start'
  
     open(file=fname,newunit=io_start,status='old', &
       action='read',iostat=ierror)

     if(ierror > 0) stop ' - Missing qcxms.start file! -'

     read (io_start,*) itrj,ndum
     read (io_start,'(2D22.14)') eimp,tadd
  
     if (ndum /= nat) stop '- error in rdstart -'
  
     do j = 1, nat
        read (io_start,'(7D22.14)') xyz(1:3,j), velo(1:3,j), velof(j)
     enddo
  
     close(io_start)
  
  end subroutine rdstart


end module qcxms_read_coordinates
   
