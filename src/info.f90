module qcxms_info
  use cidcommon
  use common1
  use newcommon
  use qcxms_mdinit, only: ekinet
  use xtb_mctc_accuracy, only : wp
  use xtb_mctc_symbols, only: toSymbol 
  use xtb_mctc_convert
  use xtb_mctc_constants, only: kB

  implicit none

  contains


  subroutine info_main(ntraj, tstep, tmax, Tinit, trelax, eimp0, &
      & ieeatm, iee_a, iee_b, btf, fimp, hacc, ELAB, ECOM, MaxColl, CollNo, CollSec,  &
      & ESI, tempESI, eTempin, maxsec, betemp, nfragexit, iseed, iprog, edistri, legacy)
      
  integer  :: ntraj,iseed(1)
  integer  :: MaxColl
  integer  :: CollSec(3),CollNo(3)
  integer  :: maxsec
  integer  :: nfragexit
  integer  :: dumprint
  integer  :: iprog
  integer  :: edistri
  integer  :: i

  real(wp) :: tstep,tmax,etempin,betemp
  real(wp) :: Tinit,trelax
  real(wp) :: eimp0
  real(wp) :: iee_a,iee_b,hacc,btf,ieeatm
  real(wp) :: fimp
  real(wp) :: ELAB, ECOM
  real(wp) :: ESI,tempESI

  logical  :: legacy

  character(len=20) :: line
  character(len=20) :: line2

  write(*,'(5(''-''),(a),5(''-''))') ' Internal program parameters '
  write(*,*)
!  if( verbose ) write(*,'('' MS Method             : '',a8)')method

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! QC PROGRAMS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print the QC Protgram used 
  call qcstring(prog,line,line2) 
  write(*,'('' QC Program            : '',(a))') trim(line)

  ! print out extra information on SQM level or Basis/functional
  write(*,'('' QC Level              : '',(a))') trim(line2)

  call qcstring(iprog,line,line2) 
  if ( method == 2 .or. method == 4 ) then
     write(*,'('' QC Prog. for EAs      : '',a)')line 
  else
     write(*,'('' QC Level for IPs      : '',a)')line2 
  endif

  if ( verbose .and. ax > 0 )then
    write(*,'('' "Fock"-exchange ax    : '',F8.3)') ax
    write(*,*)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DISPERSION
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! printout of Disperion setting
  if ( prog == 8 .and. gfnver == 3 ) then
    write(*,'('' Dispersion            : '',a2)') 'D4'
  else
    write(*,'('' Dispersion            : '',a2)') 'D3'
  endif

  ! optional printout on disperion parameters
  if ( verbose ) then
    if ( a2 /= 0 ) then
      write(*,*)
      write(*,'('' D3(BJ)             a1 : '',F8.3)')a1
      write(*,'('' D3(BJ)             a2 : '',F8.3)')a2
      write(*,'('' D3(BJ)             s8 : '',F8.3)')s8
      write(*,*)
    endif   
  endif   

  ! optional output on gCP setting (see: input.f90)
  if ( verbose ) then
    if ( gcp ) then
      write(*,*) 'gCP                   : on '
    else
      write(*,*) 'gCP                   : off '
    endif
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MO spectrum calculations
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ( method /= 3 .and. method /= 4 ) then
    dumprint = 3 ! ORCA
    if ( XTBMO ) dumprint = 7 ! XTB
    call qcstring(dumprint,line,line2) 
    write(*,'('' MO spectrum with      : '',a8)  ')line !pname(dumprint+1)
  endif

  write(*,*)
  write(*,'('' total traj.   (ntraj) : '',i4  )')ntraj
  write(*,'('' time steps    (tstep) : '',f7.2,'' fs'')')tstep 
  write(*,'('' max. sim. time (tmax) : '',f7.2,'' ps'')')tmax/1000.0_wp
  write(*,'('' Initial temp. (tinit) : '',f7.2,'' K'')')Tinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHOSE:if ( method /= 3 .and. method /= 4 ) then ! not CID

    write(*,*)
    write(*,'(12(''-''),(a),11(''-''))') ' EI settings '
    write(*,'('' E(e-impact)   (eimp0) : '',f7.2,'' eV'')'     )eimp0
    if (ieetemp > 0.0_wp ) write(*,'('' Iee-temp.   (ieetemp) : '',f7.2,'' K/Eh'')')&
    &   ieetemp
    write(*,'('' Iee/atom     (ieeatm) : '',f7.2,'' eV/atom'')')ieeatm
    if ( legacy ) write(*,*) ' --- IEE distr. pseudo-random! --- '
    write(*,'('' relax. time  (trelax) : '',f7.2,'' fs'')'     )trelax
    if( verbose ) then
      if ( edistri == 0 ) then
         write(*,'('' G iee_a               : '',f7.2,'' eV⁻¹'')')iee_a
         write(*,'('' G iee_b               : '',f7.2,'' eV'')')iee_b
      endif
      if ( edistri == 1 ) then
         write(*,'('' P iee_a               : '',f7.2,'' eV'')')iee_a
         write(*,'('' P iee_b               : '',f7.2,'' eV'')')iee_b
      endif
      write(*,'('' btf                   : '',f7.2)')btf
      write(*,'('' fimp                  : '',f7.2)')fimp
      write(*,'('' hacc                  : '',f7.2)')hacc
    endif

  ! Print CID settings
  else  
    write(*,*)
    if ( .not. TempRun ) then
      write(*,'(11(''-''),(a),11(''-''))') ' CID settings '
      if (gas%IndAtom == 7 ) then
        write(*,'('' Collision Gas         : '',a,''2'')') trim(toSymbol(gas%IndAtom))
      else
        write(*,'('' Collision Gas         : '',a)') trim(toSymbol(gas%IndAtom))
      endif

      if ( ECOM > 0.0_wp ) then
        write(*,'('' E (COM)               : '',f7.2,'' eV'')') ECOM
      else
        write(*,'('' E (LAB)               : '',f7.2,'' eV'')') ELAB
      endif

    endif

    ! General run-type
    if ( Fullauto) then
       write(*,*)
       write(*,'('' Activation Run - Type : '',(a) )') 'General'
       write(*,'('' Gas Pressure   (PGas) : '',f8.3,'' Pa'')') cell%PGas
       write(*,'('' Gas Temp.      (TGas) : '',f7.2,''  K'')') cell%TGas
       write(*,'('' Cell length  (lchamb) : '',f8.3,'' m'')' ) cell%lchamb
    endif

    ! Forced -> only M+ Gas coll.
    if ( MaxColl > 0 ) then
       write(*,'('' Run - Type            : '',(a) )') 'Manual'
       if( verbose ) write(*,*) ' Collisions between M+ and Gas '
       write(*,'('' Maximum collisions    : '', i2)') MaxColl
    endif

    ! Forced -> collisions for runs (1,2,3) set (see main.f90)
    if ( CollNo(1) > 0 ) then
       write(*,'('' Run - Type            : '',(a) )') 'Manual (FGC)'
       if( verbose ) write(*,*) ' ALL Collisions (incl. fgc) '
       write(*,'('' Total collisions 1    : '', i2)') CollNo(1)
    endif
    if ( CollNo(2) > 0 ) write(*,'('' Total collisions 2    : '', i2)') CollNo(2)
    if ( CollNo(3) > 0 ) write(*,'('' Total collisions 3    : '', i2)') CollNo(3)

    ! Fragmentation -> set amount of fragmentations after it finishes
    if ( CollSec(1) > 0 ) then
       write(*,'('' Run - Type            : '',(a) )') 'Fragmenting'
       if( verbose ) write(*,*) ' Until max. <CollSec> framentations '
       write(*,'('' Maximum no. fragm. 1  : '', i2)') CollSec(1)
    endif
    if ( CollSec(2) > 0 ) write(*,'('' Maximum no. fragm. 2  : '', i2)') CollSec(2)
    if ( CollSec(3) > 0 ) write(*,'('' Maximum no. fragm. 3  : '', i2)') CollSec(3)

    ! Thermal activation -> only scale temp. (comparable to EI run-rype)
    if ( Temprun ) then 
       write(*,'('' Run - Type            : '',(a) )') 'Thermal'
       if( verbose ) write(*,*) ' Scale internal energy/temperature '
       if ( ESI > 0 )write(*,'('' Scaling to E(int)     : '', f7.2, '' eV'')') ESI
       if ( tempESI > 0 )write(*,'('' Scaling to Temp.      : '', f7.2,'' K'')') tempESI
    endif
        
  endif CHOSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)

  ! Give extra information
  if ( verbose ) then
    if (eTempin == -1 ) then
      write(*,'('' eTemp                 : '',(a))') ' set automatically!' 
    else
      write(*,'('' eTemp                 : '',f7.2, '' K'')') eTempin
    endif

    write(*,'('' base eTemp            : '',f7.2,'' K'')')betemp
    write(*,'('' # ion tracks          : '',i4  )')maxsec
    write(*,'('' nfragexit             : '',i4  )')nfragexit
    write(*,'('' iseed                 : '',i4  )')iseed
  endif

  write(*,'(44(''-''))')
  write(*,'('' qc path              '',a   )')trim(path)
  write(*,'(44(''-''))')

  ! Set xTB environment
  if ( prog == 6 .or. prog == 7 .or. prog == 8 ) then
    xtbhome=''
    call getenv('XTBHOME',xtbhome)
    if(xtbhome.eq.'') xtbhome='~/.XTBPARAM/'
    i=len(trim(xtbhome))
    if(xtbhome(i:i).ne.'/') xtbhome(i+1:i+1)='/'
    write(*,'('' xtbhome directory    '',a)')trim(xtbhome)
    write(*,'(44(''-''))')
  endif
  write(*,*)

  end subroutine info_main

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Sum up information after run is concluded (main.f90)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine info_sumup(ntraj, tstep, tmax, Tinit, trelax, eimp0, &
      & ieeatm, iee_a, iee_b, ELAB, ECOM, ESI, tempESI,  & 
      & nfragexit, iprog, nuc, velo, mass)

  integer  :: ntraj
  integer  :: iprog
  integer  :: edistri
  integer  :: nfragexit
  integer,intent(in) :: nuc
  integer :: io_log

  real(wp) :: tstep,tmax
  real(wp) :: Tinit,trelax
  real(wp) :: eimp0
  real(wp) :: iee_a,iee_b,ieeatm
  real(wp) :: ELAB, ECOM
  real(wp) :: ESI,tempESI
  real(wp), intent(in)  :: velo(3,nuc),mass(nuc)
  real(wp) :: Ekin,Temp, E_int

  character(len=20) :: line
  character(len=20) :: line2

  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'if this is not a restart (ie if you dont'
  write(*,*)'want to improve statistics), delete <qcxms.res>!'
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)
  write(*,*)'re-call of important run parameters:'
  write(*,*)

  call qcstring(prog,line,line2) 
  write(*,'('' QC Program            : '',(a))') trim(line)

  ! print out extra information on SQM level or Basis/functional
  write(*,'('' QC Level              : '',(a))') trim(line2)

  if ( prog /= iprog ) then
    call qcstring(iprog,line,line2) 
    if ( method == 2 .or. method == 4 ) then
       write(*,'('' QC Prog. for EAs      : '',a)')line 
    else
       write(*,'('' QC Level for IPs      : '',a)')line2 
    endif
  end if
  write(*,*)

  write(*,'('' total traj.   (ntraj) : '',i4  )')ntraj
  write(*,'('' time steps    (tstep) : '',f7.2,'' fs'')')tstep * autofs
  write(*,'('' max. sim. time (tmax) : '',f7.2,'' ps'')')tmax/1000.0_wp
  write(*,'('' Initial temp. (tinit) : '',f7.2,'' K'')')Tinit

  !write out informations
info: if ( method /= 3 .and. method /= 4 )then

    write(*,*)
    write(*,'(12(''-''),(a),11(''-''))') ' EI settings '
    write(*,'('' E(e-impact)   (eimp0) : '',f7.2,'' eV'')'     )eimp0 * autoev
    if (ieetemp > 0.0_wp ) write(*,'('' Iee-temp.   (ieetemp) : '',f7.2,'' K/Eh'')')&
    &   ieetemp
    write(*,'('' Iee/atom     (ieeatm) : '',f7.2,'' eV/atom'')')ieeatm
    write(*,'('' relax. time  (trelax) : '',f7.2,'' fs'')'     )trelax
    if( verbose ) then
      if ( edistri == 0 ) then
         write(*,'('' G iee_a               : '',f7.2,'' eV⁻¹'')')iee_a
         write(*,'('' G iee_b               : '',f7.2,'' eV'')')iee_b
      endif
      if ( edistri == 1 ) then
         write(*,'('' P iee_a               : '',f7.2,'' eV'')')iee_a
         write(*,'('' P iee_b               : '',f7.2,'' eV'')')iee_b
      endif
    endif

  ! Print CID settings
  else  
    if ( .not. TempRun ) then
      write(*,*)
      write(*,'(11(''-''),(a),11(''-''))') ' CID settings '
      write(*,'('' Collision Gas         : '',a2)')      toSymbol(gas%IndAtom)
      if (ECOM > 0.0_wp) then
        write(*,'('' E (COM)               : '',f7.2,'' eV'')') ECOM
      else
        write(*,'('' E (LAB)               : '',f7.2,'' eV'')') ELAB
      endif
    endif

    ! General run-type
    if ( Fullauto) then
       write(*,*)
       write(*,'('' Activation Run - Type : '',(a) )') 'General'
       write(*,'('' Gas Pressure   (PGas) : '',f8.3,'' Pa'')') cell%PGas
       write(*,'('' Gas Temp.      (TGas) : '',f7.2,''  K'')') cell%TGas
       write(*,'('' Cell length  (lchamb) : '',f8.3,'' m'')' ) cell%lchamb
    endif

    write(*,*)

    call ekinet(nuc,velo,mass,ekin,temp)
    E_int = (temp * (0.5 * 3 * nuc * kB)) * autoev

    write(*,*)
    write(*,'('' internal Energy        : '',F6.2,'' eV'',a3)')E_int
    if (ESI > 0)     write(*,'('' Scaling to             : '',F6.2, '' eV'')')ESI
    if (tempESI > 0) write(*,'('' Scaling to             : '',F6.2, '' K'')')tempESI

  endif info

  if (verbose) write(*,'('' nfragexit             : '',i8  )')nfragexit


  !open(file='qcxms.log',unit=io_log,status='old')
  !! write into qcxms.log - no information I consider important till now
  !if ( method /= 3 .and. method /= 4 ) then
  !  write ( io_log, '(''prog          '',a6)  ' ) line
  !  write ( io_log, '(''ntraj         '',i8  )' ) ntraj
  !  write ( io_log, '(''tstep  (fs)   '',f8.2)' ) tstep*autofs
  !  write ( io_log, '(''eimp0  (eV)   '',f8.2)' ) eimp0*autoev
  !  write ( io_log, '(''ieeatm (eV)   '',f8.2)' ) ieeatm
  !  write ( io_log, '(''iee_a         '',f8.2)' ) iee_a
  !  write ( io_log, '(''iee_b         '',f8.2)' ) iee_b
  !  write ( io_log, '(''nfragexit     '',i8  )' ) nfragexit
  !  write ( io_log, '(''QC method: '',a)')line2
  !  !call version(2)
  !  close ( io_log )
  !endif

  end subroutine info_sumup

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check if CID inputs are good
  ! Fancy output, better for users to understand whats wrong
  ! Also, QCxMS run #1 is not influenced
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check if input makes sense before ESI calc. starts
  subroutine cidcheck(MaxColl, CollNo, CollSec)

  integer  :: MaxColl
  integer  :: CollNo(3)
  integer  :: CollSec(3)
  integer  :: i

  logical  :: errorCIDinfo
  logical  :: logical_CollNo 
  logical  :: logical_CollSec

  errorCIDinfo = .false.
  logical_CollNo = .false. 
  logical_CollSec = .false. 

  ! Check manual settings
  do i = 1, 3
    if ( CollNo(i)  > 0 ) logical_CollNo = .true.
    if ( CollSec(i) > 0 ) logical_CollSec = .true.
  enddo

  ! Check if both CollNo and CollSec are set
  if ( logical_CollNo .and. logical_CollSec ) then
    errorCIDinfo = .true. 
    write(*,*) 'S T O P - CollNo and CollSec are mutually exclusive!'
  endif

  ! Check CollNo or CollSet NOT with MaxColl
  if ( (logical_CollNo .or. logical_CollSec) .and. Maxcoll /= 0 ) then 
    errorCIDinfo = .true. 
    write(*,*) 'S T O P - Maxcoll and CollNo/CollSec are mutually exclusive!'
  endif

  ! Check if run-types are simulatniousely set  
  if ( CollAuto .and. FullAuto ) then
    errorCIDinfo = .true. 
    write(*,*) 'S T O P - Cannot set CollAuto and FullAuto !'
    write(*,*) ' Choose one of the two automatic procedures!'
  endif

  if ( (CollAuto .or. FullAuto) .and. TempRun ) then
    errorCIDinfo = .true. 
    write(*,*) 'S T O P - Cannot set Auto and TempRun!'
    write(*,*) ' Choose one of the automatic procedures!'
  endif

  ! Check Collauto settings
  if ( Collauto .and.( Maxcoll /= 0 .or. logical_CollNo .or. logical_CollSec ))then
    errorCIDinfo = .true. 
    write(*,*) 'S T O P - CollAuto is an automatic procedure!'
    write(*,*) 'Change the number of collisions in CollAuto with &
      &SetColl instead !'
  endif

  ! Check Fullauto settings
  if ( FullAuto .and. (Maxcoll /= 0 .or. logical_CollNo .or. logical_CollSec ))then
    errorCIDinfo = .true. 
    write(*,*) 'S T O P - FullAuto is an automatic procedure!'
    write(*,*) 'The number of collisions is dependend on the collision &
      &cell settings (PGas, TGas, lchamb) !'
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Stop the program and print out exiting the program
  if (errorCIDinfo)  stop  '   - - - E  X  I  T - - -  '
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine cidcheck


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write program strings
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine qcstring(iprog,line,line2)
  
  integer :: iprog
  
  character(len=20) :: line
  character(len=20) :: line2
  character(len=30) :: basis(14),df(0:15)
  character(len=20) :: lab
  
  basis(1)  = 'SV'
  basis(2)  = 'SVx'
  basis(3)  = 'SV(P)'
  basis(4)  = 'SVP'
  basis(5)  = 'TZVP'
  basis(6)  = 'ma-def2-SVP'
  basis(7)  = 'ma-def2-TZVP'
  basis(8)  = 'ma-def2-TZVPP'
  basis(9)  = 'def2-SV(P)'
  basis(10) = 'def2-SVP'
  basis(11) = 'def2-TZVP'
  basis(12) = 'def2-QZVP'
  basis(13) = 'QZVP'
  basis(14) = 'ma-def2-QZVP'
  df   (0)  = 'PBE0'
  df   (1)  = 'PBE38'
  df   (2)  = 'PBE12'
  df   (3)  = 'LDA12'
  df   (4)  = 'M062X'
  df   (5)  = 'PBE'
  df   (6)  = 'B97D'
  df   (7)  = 'B3LYP'
  df   (8)  = 'PW6B95'
  df   (9)  = 'B3PW91'
  df   (10) = 'BLYP'
  df   (11) = 'BP86'
  df   (12) = 'TPSS'
  df   (13) = 'REVPBE'
  df   (14) = 'PBEh-3c'
  df   (15) = 'BHLYP'
  
  line  = ''
  line2 = ''
  lab   = ''

  ! D3
  if ( a2 > 0 ) lab = '-D3'
  
  !DFTB
  if (iprog == 0 ) line = 'DFTB'

  ! MOPAC
  if (iprog == 1 ) then
    line = 'MOPAC'
    if ( ihamilt == 1  ) line2 = 'AM1'
    if ( ihamilt == 2  ) line2 = 'PM3'
    if ( ihamilt == 3  ) line2 = 'RM1'
    if ( ihamilt == 4  ) line2 = 'PM6'
  endif   

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TURBOMOLE
  if (iprog == 2 ) line = 'TURBOMOLE'

  ! ORCA
  if (iprog == 3 .and. orca_version == 5 ) line = 'ORCA 5'
  if (iprog == 3 .and. orca_version == 4 ) line = 'ORCA 4'

  ! DFT functionals and basis output 
  if ( iprog == 2 .or. iprog == 3 ) then
     write(line2,'(a,''/'',a)')trim(df(func)),trim(basis(bas))
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! MSINDO
  if (iprog == 4 ) then
    line = 'MSINDO'
    line2 = 'NDDO'
  endif

  ! MNDO
  if (iprog == 5 ) then
    line = 'MNDO'
    if ( ihamilt == 1  ) line2 = 'AM1'
    if ( ihamilt == 2  ) line2 = 'PM3'
    if ( ihamilt == 6  ) line2 = 'OM2'
    if ( ihamilt == 8  ) line2 = 'OM3'
    if ( ihamilt == 10 ) line2 = 'MNDO/d'
  endif

  ! XTB
  if (iprog == 6 ) line = '(CALL) XTB'
  if (iprog == 7 .or. iprog == 8 ) then
    line = 'xTB'
    if ( iprog == 7 ) line2 = 'GFN1-xTB'
    if ( iprog == 8 ) line2 = 'GFN2-xTB'
  endif

  if ( ihamilt /= 0) write(line2,'(2a)')trim(line2),trim(lab)
  
  
  end subroutine qcstring

end module qcxms_info
