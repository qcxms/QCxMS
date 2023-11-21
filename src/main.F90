!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Copyright (C) 2009 - 2021 Stefan Grimme, University of Bonn, Germany
!
! the code has 4 main modi of operation:
!
! call as  qcxms  in a clean dir              : generates <qcxms.gs>
!                                                neutral M trajectory
! call as  qcxms  when <qcxms.gs> exists     : create inputs for production
!                                                runs, qcxms.in is read and used
! call as  qcxms -prod                        : do paralell fragmentation runs
!                                                as called by pmsprep script
! the code can be run with any of the following QC programs:
! TURBOMOLE, ORCA, MOPAC2012, DFTB+, MNDO99 or MSINDO
!
! Alternatively, the xTB TBlite library has been implemented (subprojects) and
! set as default (running GFN1- and GFN2-xTB for grad. calculations).
! The library has to be updated separately. Go into /subprojects/tblite/ and pull 
! the newest updates via 'git pull' or ' git fetch' + 'git merge'
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program QCxMS
use cidcommon
use common1
use newcommon
use readcommon
use get_version
use qcxms_boxmuller, only: vary_collisions, vary_energies
use qcxms_cid_routine, only: cid, collision_setup, calc_ECOM
use qcxms_mo_spec, only: getspec
use qcxms_fragments
use qcxms_iee, only: getieeab, getmaxiee, gauss0, gauss1, poiss0
use qcxms_impact, only: calctrelax
use qcxms_iniqm, only: iniqm
use qcxms_input, only: input, read_struc_commandline
use qcxms_mdinit, only: mdinitu, ekinet
use qcxms_molecular_dynamics, only: md
use qcxms_info, only: info_main, info_sumup, cidcheck, start_info
use qcxms_use_orca, only: copyorc
use qcxms_use_turbomole, only: copytm, copytm_ip
use qcxms_utility
use qcxms_write_fragments, only: manage_fragments 
use xtb_mctc_accuracy, only : wp
use xtb_mctc_convert
use xtb_mctc_constants
!use xtb_mctc_symbols, only: toSymbol 
use mctc_env, only : error_type, get_argument
use mctc_io, only : structure_type, new, to_symbol

!use cid_module, only: cid_outer

use get_settings

implicit none

integer  :: nmax,nmax0,dumpstep
integer  :: nfrag,ndumpGS
integer  :: ntraj
integer  :: i,j,k,m
integer  :: nrun,nfragexit
integer  :: mspin,iprog,tcont
integer  :: mchrg, mchrg_prod
integer  :: maxsec,idum,edistri
integer  :: ihomo,nb,nuc,mo1,mo2,fragstate
integer  :: num_frags
integer  :: itrj,icoll,isec !,prestep
integer  :: fragat(200,10)
integer  :: imassn(10000)
integer  :: io_res, io_gs, io_log, io_eimp
integer  :: CollSec(3),CollNo(3)
integer  :: collisions
integer  :: minmass,numb
integer  :: frag_counter,new_counter,save_counter
integer  :: rand_int,dep,convetemp
integer  :: manual_dist
integer  :: arg_list
! GBSA
!integer  :: gsolvstate

integer,allocatable  :: list(:)
integer,allocatable  :: iat (:)
integer,allocatable  :: icalc(:)
integer,allocatable  :: imass(:)

real(wp) :: etemp
real(wp) :: tstep,edum,tmax,etemp_in,exc,dums,betemp
real(wp) :: t,Tsoll,Tdum,Tinit,trelax,ttime,ehomo
real(wp) :: Epav,Ekav,Tav,dum,t1,t2,w1,w2,eimpw,pmax,tta
real(wp) :: tadd, pretadd
real(wp) :: fimp,aTlast
real(wp) :: chrgcont
real(wp) :: cema(3),eimp0,eimp,dtime
real(wp) :: iee_a,iee_b,hacc,ieeel,btf,ieeatm,etempGS
real(wp) :: vScale,new_velo,MinPot,ESI,tempESI
real(wp) :: tScale,new_temp
real(wp) :: summass,beta
real(wp) :: direc(3),cm(3),cm1(3),cm2(3)
!real(wp) :: cm_out
real(wp) :: randx
real(wp) :: x
real(wp) :: lowerbound,upperbound
real(wp) :: ELAB, ECOM
real(wp) :: gaus(1000)
real(wp) :: a,b
real(wp) :: rtot,cross,mfpath
real(wp) :: E_KIN,E_COM
real(wp) :: calc_collisions,E_Scale,ENe(10)
! Ion tracking arrays
real(wp) :: fragm(10) !10 fragments max per fragmentation,
real(wp) :: xyzn(3,10000),velon(3,10000),iatn(10000),massn(10000) !10000 atoms max
real :: snorm
real(wp), allocatable :: ergebnis(:)

! Allocatables
real(wp),allocatable :: xyz (:,:)
real(wp),allocatable :: axyz(:,:)
real(wp),allocatable :: grad(:,:)
real(wp),allocatable :: velo(:,:)
real(wp),allocatable :: velof(:)
real(wp),allocatable :: mass(:)
real(wp),allocatable :: atm_charge(:)
real(wp),allocatable :: spin(:)
real(wp),allocatable :: emo  (:)
real(wp),allocatable :: mopop(:,:)
real(wp),allocatable :: modum(:)

! Strings
character(len=:), allocatable :: res_name
character(len=80)    :: fragf(10)
character(len=120)   :: asave
character(len=:), allocatable :: inp_fname
!character*80 solvent

! logicals
logical :: ex
logical :: md_ok
logical :: scani = .false.
logical :: nfrag_ok = .true.
logical :: check,eonly,eonly0,eonly1,noeq
logical :: prod,iniok
logical :: unity ! scales velocities uniformly 
logical :: metal3d ! switches on some features for transition metal compounds 
logical :: ECP !switch on ECP for orca calculations
logical :: stopcid !if error comes up in CID module - code is killed
logical :: noecp,nometal !forbid checking of ECP and METAL
logical :: ConstVelo,eExact
logical :: small,littlemass
logical :: No_ESI,NoScale
logical :: starting_md
logical :: legacy
!logical gbsa ! set solvation model

intrinsic :: get_command_argument
external  :: system

type(structure_type) :: mol
type(error_type), allocatable :: error
type(run_settings) :: set
type(collision_type) :: coll
!type(charge_type) :: chrg


!+*************************
integer :: samples, n, ind
integer(wp) :: bin(-50:50) = 0

real(wp) :: rnd1, rnd2
real(wp) :: s
real(wp) :: z0, z1
real(wp) :: sumn, sumnsq

real(wp) :: mean, stddev
real(wp) :: Tcheck

!+*************************


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start the program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! use/uncomments this piece of code for the lfc compiler which generates other system calls
!     call execute_command_line('echo $SHELL > .tmpqcxms')
!     open(unit=1,file='.tmpqcxms')
!     read(1,'(a)')adum
!     close(1,status='delete')
!     if(index(adum,'tcsh') /= 0) shell=1
!     if(index(adum,'bash') /= 0) shell=2

call start_info
shell=2
!write(*,*)
!write(*,*) 'Shell ',shell
!write(*,*)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! some initialization:

! MS method
method = 0
! logical for too small fragments
small=.false.
littlemass=.false.
! no checking of input etc
check=.false.
! more infos
verbose = .false.
! no production
prod   =.false.
! no equilibration 
noeq   =.true.
! test calls
eonly  =.false.
eonly0 =.false.
eonly1 =.false.
! uniform velocity scaling
unity=.false.
! check if cid was OK
stopcid = .false.
starting_md=.false.
No_ESI = .false.
! HS-UHF ini only for frag runs
iniok  =.true.
! dump every dumpstep MD steps for MOLDEN movie (=4 fs as default)
dumpstep=4
! counts the number of QC calls
calls=0
! set scaling temp to 0
tscale = 0.0_wp

! GBSA Solvation Model
!solvent='none'
!gsolvstate=0

new_velo = 0.0d0
! total MD time including ion tracks
ttime=0
dtime = 0.0d0
! undefined spin state
mspin=0
! neutral M
mchrg=0
! is fragmented?
fragstate=0      

t1 = 0.0_wp
t2 = 0.0_wp
w1 = 0.0_wp
w2 = 0.0_wp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! R   E   A   D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the INPUT file
! input, all defaults are set here
call input(tstep,tmax,ntraj,etemp_in,Tinit, mchrg_prod,           &
&          iee_a,iee_b,eimp0,eimpw,fimp,iprog,                            &
&          trelax,hacc,nfragexit,maxsec,edistri,btf,ieeatm,               &
&          scani,lowerbound,upperbound,                                   &
&          ELAB,ECOM,eExact,ECP,unity,noecp,nometal,                      &
&          vScale,CollNo,CollSec,ConstVelo,                               &
&          minmass,etempGS,coll,               &
&          MinPot,ESI,tempESI,No_ESI,NoScale,manual_dist,legacy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call read_struc_commandline(mol, check, prod, noeq, eonly0, eonly1, eonly, inp_fname)

nuc = mol%nat !      nuc : Number of atoms in molecule

! allocate memory for molecular properties
allocate(xyz (3,nuc),      &
&        axyz(3,nuc),      &
&        grad(3,nuc),      &
&        velo(3,nuc),      &
&        velof(nuc),       &
&        emo(10*nuc),      &
&        mopop(5*nuc,nuc), &
&        modum(5*nuc),     &
&        atm_charge(nuc),  & !      atm_charge(nuc) : Charge of atom
&        spin(nuc),        &
&        iat (nuc),        &
&        list(nuc),        &
&        imass(nuc),       &
&        mass(nuc))           !      mass(nuc) : Mass of the atom

xyz = mol%xyz
iat = mol%num(mol%id) !      iat(nuc) : typ of atom (index)

! Description:
!      ncao / nbf : number of cartesian AOs
!      nel : number of valence electrons (-1)
!      ndim / : numer of AOs
!      ihomo : integer of HOMO orbital
!      edum : Energy of the system (3/2*k_B*T*nuc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Choose the MS method (DEA -> method = 2 for legacy support)
write(*,*)
if (method  ==  0) then
  write(*,'(14x,60(''*''))')
  if ( mchrg_prod > 0 ) write(*,'(22x,'' Mode:          Electron Impact (EI)         '')')
  if ( mchrg_prod < 0 ) write(*,'(22x,'' Mode: Dissociative Electron Attachment (DEA)'')')
  if ( mchrg_prod == 0 ) then
    write(*,*)
    stop " -- Can't use no charge as input! Select charge! -- "
  endif
  write(*,'(14x,60(''*''))')
elseif (method  ==  1) then !This will be closed shell with EI
  write(*,'(14x,60(''*''))')
  write(*,'(22x,'' Mode:         Closed Shell Cation (CSC)     '')')
  write(*,'(14x,60(''*''))')
!  elseif (method  ==  2) then
!    write(*,'(14x,60(''*''))')
!    write(*,'(22x,'' Mode: Dissociative Electron Attachment (DEA)'')')
!    write(*,'(14x,60(''*''))')
elseif (method  ==  3) then
  write(*,'(14x,60(''*''))')
  write(*,'(22x,'' Mode: Collision Induced Dissociation (CID)  '')')
  if ( mchrg_prod > 0 ) write(*,'(22x,''           + Positive Ion mode +'')')
  if ( mchrg_prod < 0 ) write(*,'(22x,''           - Negative Ion mode -'')')
  if ( mchrg_prod == 0 ) then
    write(*,*)
    stop " -- Can't use no charge as input! Select charge! -- "
  endif
  write(*,'(14x,60(''*''))')
!elseif (method  ==  4) then
!   write(*,'(14x,60(''*''))')
!   write(*,'(22x,'' Mode: Collision Induced Dissociation (CID)  '')')
!        write(*,'(14x,60(''*''))')
else
   write(*,'(22x,''**************************'')')
   write(*,'(22x,''* Error Error Error      *'')')
   write(*,'(22x,''* Method Option failed?  *'')')
   write(*,'(22x,''**************************'')')
   stop 'Something went terribly wrong!'
endif
write(*,*)

! initialize D3
if (prog /= 8) call copyc6

! set traj automatically
if(ntraj <= 0) ntraj = nuc*25

! get for the QC method used the eTemp as info
call setetemp(1,-1.0d0,betemp)


  !> inialize random numbers
  !>> if not exlicitly set, use true random
#ifdef __INTEL_COMPILER
  if (iseed(1) == 0) then
    call random_seed()
  else
  !>> if explicitly set, use seed number from input
    call random_seed(put=iseed)
  endif
  call random_seed (put=iseed)
#else
  ! FIXME: this is probably not the right way to do this
  call random_seed (put=spread(iseed(1), 1, 8))
#endif


call random_number(randx)

! for GS every n steps a struc. is used later for production, i.e.
! 100 fragmentation runs require about 5000 GS MD steps
! this ensures uncorrelated fragmentations
nmax0 = ntraj * 50


! Test etemp for CID - hardcoded for XTB for now
if (method == 3 .and. prod)then
   if(etemp_in <= 0)then
      betemp = 5000.0_wp 
      !betemp = 2500.0_wp !5000
      !etemp_in = 2500.0_wp
   else
      betemp = etemp_in
   endif
endif


! change time to fs      
tmax = tmax*1000.0_wp
! # MD steps in a frag run
nmax = tmax / tstep
!! for CID make dependend on atom size
!if(method == 3) then
!  nmax = nuc * 100
!  tmax = nmax * tstep
!endif

! timesteps in au
tstep = tstep * fstoau

! the "70" eV in a.u. for EI
eimp0 = eimp0 * evtoau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! printing runtype information and chosen parameters
call info_main(ntraj, tstep/fstoau, tmax, Tinit, trelax, eimp0*autoev, mchrg, &
  mchrg_prod, ieeatm, iee_a, iee_b, btf, fimp, hacc, ELAB, ECOM, coll%max_coll, &
  CollNo, CollSec, ESI, tempESI, etemp_in, maxsec, betemp, nfragexit, iprog, &
  edistri, legacy)

!if (prog == 2 .or. iprog == 2) call wrcoord(nuc, xyz, iat)

! check for 3d transition metals and 4s-4d/5s-5d elements
if (.not. noecp)then
   do i=1,nuc
      if(iat(i) <= 86.and.iat(i) >= 37) ECP=.true.
   enddo
endif

! convert SV,SV(P),SVP,TZVP to same-size def2 basis set if ECP is on
! for both orca and turbomole.
if ( prog == 2 .or. prog == 3 ) then
   if ( ECP ) then
      write(*,*)
      write(*,*)'Found 5s-4d-5p/6s-5d-6p elements.'
      write(*,*)'ECPs will be used!'
      write(*,*)
      if(bas == 1 .or. bas == 2 .or. bas == 3 .or. bas == 6)then
         write(*,*) 'basis set is changed to def2-SV(P)'
         bas = 9
      elseif(bas == 4 .or. bas == 10)then
         write(*,*) 'basis set is changed to def2-SVP'
         bas = 10
      elseif(bas == 5 .or. bas == 11 .or. bas == 12)then
         write(*,*) 'basis set is changed to def2-TZVP'
         bas = 11
      endif
   endif
endif

! read file qcxms.in for isotopes: after keyword
! isotope <atom number> <mass>  <atom number> <mass> <atom number> <mass> ...
call read_isotopes(nuc,imass)

! assign masses
call setmass(nuc,iat,mass,imass)

if (.not.prod .and. verbose ) then
   write(*,*)
   write(*,*)'molecule:'

   do i=1,nuc
      write(*,'(i3,3f13.6,4x,a2,f9.4)') &
      &     i,xyz(1,i),xyz(2,i),xyz(3,i),to_symbol(iat(i)),mass(i) * autoamu
   enddo

   write(*,*)
endif


! this is to test the QM in SP calcs
if(eonly.or.eonly0.or.eonly1)then
   if(eonly)then
      inquire(file='.CHRG',exist=ex)
      if(ex)then
         open(unit=11,file='.CHRG')
         read(11,*)dum
         close(11)
         mchrg=int(dum)
      endif
      inquire(file='.UHF',exist=ex)
      if(ex)then
         open(unit=11,file='.UHF')
         read(11,*)dum
         close(11)
         mspin=int(dum)+1
      endif
   endif
   if(eonly0) mchrg=0
   if(eonly1) mchrg=1
   idum=0
   call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)
   if(.not.iniok) stop 'fatal QC error. Must stop!'
   write(*,*) 'test energy written to <.ENERGY>'
   open(unit=11,file='.ENERGY')
   write(11,*)edum
   close(11)
   call version(2)
   stop 'normal termination of QCxMS'
endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! skip if production runs are selected (parallel)
prun: if(.not.prod) then 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set the charge of the chosen method 
  ! 0 = EI , 1 = CSC,  3 = CID
  if     (method  ==  0 ) then 
    mchrg = 0
  elseif (method  ==  1 .or. method  ==  3) then 
    mchrg = mchrg_prod 
  endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! neutral GS trajectory to generate many
! different fragmentation trajectories
! initial velo (uniform) distribution on all atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  inquire(file='qcxms.gs',exist=ex)
GS: if(.not.ex)then

    velo = 0.0_wp

    ! Calc. inner Energy at given energy Tinit (500K default)
    ! 3/2*n*k_B*T (* number of atoms) (in Hartree):
    edum =   3.0_wp * 0.5_wp * kB * Tinit * nuc

    ! initilize the velocities, distribute among the atoms
    call mdinitu(nuc,velo,mass,edum)

    ! init the QM code  ! Here edum is total energy, not inner energy anymore!
    write(*,*) 'Groundstate eTemp set to',etempGS
    call iniqm(nuc,xyz,iat,mchrg,mspin,etempGS,edum,iniok,ECP)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! equilibrate GS
    if(.not.iniok) stop 'fatal QC error. Must stop!'
    
    write(*,*)'M trajactory, equilibration ...',nmax0/2

       call md(-1,0,0,nuc,nmax0/2,xyz,iat,mass,imass,mchrg,grad,     &
       &       velo,velof,list,tstep,ndumpGS,99,                     &
       &       fragm,fragf,fragat,dumpstep,etempGS,                  &
       &       md_ok,atm_charge,spin,axyz,                             &
       &       Tinit,0.0d0,0.0d0,.false.,Tav,Epav,Ekav,ttime,aTlast, &
       &       fragstate,dtime,ECP,.false.,0.0d0)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do MD GS
    ndumpGS=0

    write(*,*)
    write(*,*)'M trajactory, sampling ...',nmax0

       call md(0,0,0,nuc,nmax0,xyz,iat,mass,imass,mchrg,             &
       &       grad,velo,velof,list,tstep,ndumpGS,99,                &
       &       fragm,fragf,fragat,dumpstep,etempGS,                  &
       &       md_ok,atm_charge,spin,axyz,                             &
       &       Tinit,0.0d0,0.0d0,.false.,Tav,Epav,Ekav,ttime,aTlast, &
       &       fragstate,dtime,ECP,.false.,0.0d0)


    write(*,'('' average Ekin '',F12.6)')Ekav
    write(*,'('' average Epot '',F12.6)')Epav
    write(*,'('' average Etot '',F12.6)')Epav+Ekav
    write(*,'('' average T    '',F12.1)')Tav

    if (method == 3) then
       call ekinet(nuc,velo,mass,edum,t)
       ENe(3) =  ( t * (0.5_wp * 3 * nuc * kB)) * autoev 
       write(*,'('' Inner E (eV) '',F12.4)')ENe(3)
    endif

    write(*,*)'structures for fragmentation runs written'
    call version(2)
    stop 'normal termination of QCxMS'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! END QCxMS groundstate call (1)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  else ! to re-calculate .gs file if not want to calc TMPQCXMS files
      ! no idea why this is done, but still keeping it (for nostalgia or smth, idk)

    write(*,*) 'reading qcxms.gs ...'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ignored if not used (default)
    ! do MD using previous GS traject
    if ( .not. noeq ) then
       ! Read the groundsate traj. file
       write(*,*) 'Ground-State trajactory, sampling ...'
       open(file='qcxms.gs',unit=io_gs,status='old')
       read(io_gs,*) ndumpGS

       do k=1,nuc
          read(io_gs,*)(xyz(j,k),j=1,3),(velo(j,k),j=1,3)
       enddo

       close(io_gs,status='delete')
       call md(0,0,0,nuc,nmax0,xyz,iat,mass,imass,mchrg,             &
       &       grad,velo,velof,list,tstep,ndumpGS,99,          &
       &       fragm,fragf,fragat,dumpstep,etempGS,                  &
       &       md_ok,atm_charge,spin,axyz,                             &
       &       Tinit,0.0d0,0.0d0,.false.,Tav,Epav,Ekav,ttime,aTlast, &
       &       fragstate,dtime,ECP,.false.,0.0d0)


       write(*,'(/,'' average Ekin '',F12.6)')Ekav
       write(*,'('' average Epot '',F12.6)')Epav
       write(*,'('' average Etot '',F12.6)')Epav+Ekav
       write(*,'('' average T    '',F12.1)')Tav
       write(*,*)'structures for fragmentation runs written'
       call execute_command_line('rm -rf fort.11 fort.15')
       call execute_command_line('rm -rf charges.bin detailed.outi dftb_in.hsd')
       call execute_command_line('rm -rf dftb_pin.hsd band.out')
       call version(2)
       stop 'normal termination of QCxMS'
    endif
    ! ignored if not used (default)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  endif GS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! START QCxMS set-up call (2) 
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call timing(t1,w1)

  !> Check if CID input makes sense before creating the dirs.
  if ( method == 3 ) then 
    call cidcheck(coll%max_coll, CollNo, CollSec )
  endif

  !> set the charge to input value (default: +1)
  mchrg = mchrg_prod 

  !> Test SP energy calc. of the used QM method on the input geom with charges (+1,-1)
  !> (It is important to do for calcs with TM and ORCA)
  write(*,'(''--- Checking QC method for ions ---'')')
  if (etemp_in < 0 ) call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)
  if (etemp_in > 0 ) call iniqm(nuc,xyz,iat,mchrg,mspin,etemp_in,edum,iniok,ECP)
  if (iniok) then 
    write(*,'(''--- QC method okay ---'')')
  else
    stop '>>> fatal QC error. Must stop! <<<'
  endif

  call timing(t2,w2)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !generate the starting points for production run
  write(*,*)
  write(*,*)'*******************************************'
  write(*,*)'      generating randomized ensemble '
  write(*,*)'*******************************************'
  write(*,*)'reading inital coord/velo from qcxms.gs ...'

  open (file='qcxms.gs',newunit=io_gs, status='old')
  read(io_gs,*) ndumpGS

  ! allocating needed variables
  allocate(set%xyzr  (3,nuc,ntraj),       &
  &        set%velor (3,nuc,ntraj),       &
  &        set%velofr(  nuc,ntraj),       &
  &        set%eimpr(ntraj), &
  &        set%taddr(ntraj), &
  &        icalc(ndumpGS))


  ! check if GS was sampled long enough
  if ( ndumpGS <= 2*ntraj ) stop 'Error: compute longer GS trajectory'

  icalc(1:ndumpGS) = 0
  k = 0

  ! determine which structures are taken, randomly on the GS trj
  ! do not take structures equidistant steps to avoid correlations
  ! irand -> function: (0 - ndumpGS) numbers (utilitly.f90)
  do while (k <= ntraj)
    j = irand(ndumpGS)
    if ( icalc(j) == 0) then
      k = k + 1
      icalc(j) = 1
    endif
  enddo

  gaus = 0
  list = 0
  nrun = 0
  Tav  = 0
  tta  = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  ! Calculate or set IEE distribution for EI calculations 
  ! This is only important for the "standard values"
  ! Only pmax is used further as the maximum prob. for a given IEE energy (see below)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  ! NOT needed for CID, no need for IEE
iee0: if (method /= 3 ) then 

    emo = 0

    ! mchrg = 0 to get the values for ground-state
    !call getspec(.true.,nuc,iat,xyz,0,eps,ehomo,mopop,ihomo,nb,ECP) !eps has to be changed in geteps 
    call getspec(.true.,nuc,iat,xyz,0,emo,ehomo,mopop,ihomo,nb,ECP)

    ! adjust IEE distr. if necessary (max in intervall)
    if (scani) then
       write(*,'(/,'' IEE scanning feature  ...'')')
    else
       write(*,'(/,'' preparing the IEE distribution ...'')')
    endif

    exc   = (eimp0 - ehomo) * autoev
    ieeel = dble (ihomo + nb)

    ! user input
iee1: if ( iee_a > 0 .and. iee_b > 0 ) then
      call getmaxiee(iee_a,iee_b,ieeel,edistri,exc,dum,pmax,dums)
      write(*,'('' IEE a (eV) [set]               : '',f8.3)')  iee_a
      write(*,'('' IEE b (eV) [set]               : '',f8.3)')  iee_b

    ! automatic
    else
       call getieeab (iee_a,iee_b,ieeel,edistri,exc,nuc,ieeatm)
       call getmaxiee(iee_a,iee_b,ieeel,edistri,exc,dum,pmax,dums)
       write(*,'('' IEE a (eV) [determined]        : '',f8.3)') iee_a
       write(*,'('' IEE b (eV) [determined]        : '',f8.3)') iee_b

    endif iee1 ! EI

    if (scani) then
       write(*,'('' minimum IEE (eV)               : '',f8.3)') lowerbound
       write(*,'('' maximum IEE (eV)               : '',f8.3)') upperbound
    else
       write(*,'('' maximum IEE (eV)               : '',f8.3)') exc
       write(*,'('' maximum of P(E) at (eV)        : '',f8.3)') dum
       write(*,'('' average E for P(E) (eV)        : '',f8.3)') dums
    endif
    write(*,*)

    !> For the different structures, use the above settings
iee2:do i = 1, ndumpGS
      !> read structure and velo from previous GS MD
      do k = 1, nuc
        read(io_gs,*, iostat = iocheck)(xyz(j,k),j=1,3),(velo(j,k),j=1,3)
        if (iocheck>0)then     !Fail
          write(*,*) 'Something is wrong in the input structure.'
          write(*,*) ' --- Exiting --- '
          stop
         ! End-of-file
        elseif (iocheck<0)then !EOF
          exit
        endif
      enddo

      if (icalc(i) == 0) cycle
      if (nrun == ntraj) exit

      !> 1. generate an e- that can ionize any MO and vary it by boxuller
      do
        !> Legacy support. Do runs with SAME random seed for each run (pseudo-random)
        if ( legacy ) then
          Edum = eimp0 + eimpw * eimp0 * snorm() 

        !> DEFAULT: vary by normal distributed boxmuller random number
        else
          Edum = vary_energies(eimp0,eimpw)
        endif

        if ( Edum >= ehomo ) exit
      enddo

      ! 2. subtract E(HOMO) to get IEE
      Edum = Edum - ehomo

      ! the random IEE in the possible interval
      do 
        call random_number(randx)
        exc = randx * Edum * autoev

        if(edistri == 0)then
          call gauss0(iee_a,iee_b,ieeel,exc,dum)
        else
          call poiss0(iee_a,iee_b,ieeel,exc,dum)
        endif
      
        ! the probability for this energy
        dum = dum / pmax
        call random_number(randx)
        if(dum >= randx)exit
      enddo

      !> make the IEE range dependend on ntraj and nrun
      if (scani) then
        exc = lowerbound + (((upperbound-lowerbound)/dble(ntraj))*dble(nrun))
      endif

      ! this is the IEE energy used -> Edum = exc
      ! fimp is a scaling factor (default = 1.0)
      ! exc is now randomly chosen, but at least larger than ehomo
      edum = fimp * exc * evtoau 

      ! map IEE to MO number, shake-up possibility (ie MO2>0)
      call momap(ihomo, emo, edum+ehomo, mo1, mo2)

      ! Temperature of the molecule
      Tsoll = Edum / (0.5_wp * 3.0_wp * nuc * kB)

      ! Avg Temperature count
      Tav = Tav + Tsoll

      ! plot (for checkin purposes)
      call gauss(gaus,edum*autoev)

      ! heating scale factors from MO populations 
      modum(1:nuc) = mopop(mo1,1:nuc)
      if (mo2 > 0) modum(1:nuc) = modum(1:nuc) + mopop(mo2,1:nuc)

      ! more H move
      do m = 1, nuc
         if (iat(m) == 1) modum(m) = modum(m) * hacc
      enddo

      ! normalize max to unity, actual v-scaling determined in impactscale
      dum = 0
      do m = 1, nuc
         if (modum(m) > dum) dum = modum(m)
      enddo
      !velof(1:nuc)=modum(1:nuc)/dum !some atoms becomes artificially slow based on MO
      velof(1:nuc) = 1 + modum(1:nuc) / dum !some atoms become artificially fast based on MO

      ! for 'big' systems scale velocities uniformly (unity) or if specified as command line argument
      if (nuc > 35) then !if# atoms is greater than 35 - uniform velocity scaling is employed by default
        velof(1:nuc) = 1.0_wp
      elseif (unity) then
        velof(1:nuc) = 1.0_wp
      else
        velof(1:nuc) = modum(1:nuc) / dum
      endif

      ! ad hoc parameter trelax to relate impact energy to relaxation time, randomly broadened by 10 %
      call calctrelax(nuc,emo,ihomo,mo1,trelax,tadd)
      dum = tadd
      if(mo2 > 0)then
        call calctrelax(nuc,emo,ihomo,mo2,trelax,tadd)
        dum = dum+tadd
      endif

      dum = dum / (ihomo + nb)
      ! not less than
      dum = max(dum, trelax / 10.0_wp)

      !> heating time 
      tta = tta + dum

      nrun = nrun + 1

      if (.not.check) then
        call setetemp(1,Edum,etemp)
        write(*,'('' Run: '',i4,'', step on M trj: '',i5,'' MOs '',2i3,'' IEE &
          & (eV)  = '',F6.1,'' heating time ='',F6.0,'' eTemp ='',F6.0)')     &
          & nrun,i,mo1,mo2,Edum*autoev,dum,etemp
      endif


      set%xyzr (1:3,1:nuc,nrun) = xyz (1:3,1:nuc)
      set%velor(1:3,1:nuc,nrun) = velo(1:3,1:nuc)
      set%velofr(   1:nuc,nrun) = velof(   1:nuc)
      set%eimpr (         nrun) = Edum
      set%taddr (         nrun) = dum*fstoau

    enddo iee2

    close(io_gs)
    !----

    dum  = 0
    dums = Tav * 0.5 * 3 * nuc * kB / dble(nrun)

    do i = 1, nrun
       dum = dum + (set%eimpr(i) - dums)**2
    enddo

    write(*,*)
    write(*,*)nrun,' done.'
    write(*,*)

    edum = Tav * 0.5_wp * 3 * nuc * kB / dble(nrun)

    ! output
    !open(file='qcxms.log',newunit=io_log, status='replace')
    write(*,*)
    write(*,'('' average IEE (eV) +/-  '',2F9.2)')edum*autoev,sqrt(dum/float(nrun-1))*autoev
    write(*,'('' average T (K)          '',F9.1)')Tav/dble(nrun)
    write(*,'('' average heating t (fs) '',F9.1)')tta/dble(nrun)
    write(*,'('' tmax  ps               '',f9.1)')tmax/1000.
    write(*,'('' # valence electrons    '',i4  )')ihomo+nb
    !write(io_log,'('' average IEE (eV) +/-  '',2F9.2)')edum*autoev,autoev*sqrt(dum/float(nrun-1))
    !write(io_log,'('' average T (K)         '',F9.1)')Tav/dble(nrun)
    !write(io_log,'('' average heating t (fs)'',F9.1)')tta/dble(nrun)
    !write(io_log,'('' tmax  ps              '',F9.1)')tmax/1000.
    !write(io_log,'('' # valence electrons   '',i4  )')ihomo+nb

    ! Check the IEE energy
    if ( check ) then
      open(file='eimp.dat', newunit= io_eimp)

      do i = 1, 1000
        write(io_eimp,*) i*70.0_wp/1000.0_wp,gaus(i)
      enddo

      close(io_eimp)

      write(*,*)
      write(*,*)'IEE distribution written to file <eimp.dat>'
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---- CID part --- 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> other kind of distribution, because we don't know timings etc
    elseif ( method == 3 ) then 

      if (No_ESI .and. ESI > 0) error stop 'Either ESI or not ...'

      !> do the ESI temperature scaling
      if (.not. No_ESI) then
        sumn = 0
        sumnsq = 0
        !>> empirical temperature (and energy) to which the ion is scaled 
        !>> (we do NOT want fragmentation by this scaling) 
        new_temp = 6.0_wp * nuc + 600_wp
        s = new_temp !for write-out later

        !>> calculate user input temp. scaling value
        if (ESI > 0)     tScale = ESI*evtoau / (0.5_wp * 3.0_wp * nuc * kB) 
        if (tempESI > 0) tScale = tempESI
        if (ESI > 0 .and. tempESI > 0) error stop 'Cannot provide ESI and tempESI!'

        !>> add user temp. to scaling 
        new_temp = new_temp + tScale
        
        !>> calculate the energy for Box-Muller distribution
        Edum = new_temp * (0.5_wp * 3.0_wp * nuc * kB) 

        allocate (ergebnis(ntraj))

        !>> vary energy for each traj by Box-Muller and calc mean
        do i = 1, ntraj
          ergebnis(i) = vary_energies(Edum, eimpw)
          Tcheck = ergebnis(i) / (0.5_wp * 3.0_wp * nuc * kB)
          if(verbose) write(*,*) i, Tcheck, ergebnis(i)*autoev
          sumn = sumn + ergebnis(i)
          sumnsq = sumnsq + (ergebnis(i)*ergebnis(i))
        enddo

        write(*,*) 

        mean = sumn / ntraj
        stddev = sqrt(sumnsq/ntraj - mean*mean)

        Edum = s * (0.5_wp * 3.0_wp * nuc * kB) 
        !>> printout
        write(*,'(41(''-''))')
        write(*,"(a,f10.2)")           "Starting in. Energy (eV)  : ", Edum * autoev 
        if(ESI>0) write(*,"(a,f10.2)") "Additional Energy (eV)   : ", ESI
        if(tempESI>0) write(*,"(a,f10.2)") "Additional Temp. (K)     : ", tScale
        write(*,"(a,f10.2)")           "Energy spread by... (eV) : ", eimpw 
        write(*,'(41(''-''))')
        Tcheck = mean / (0.5_wp * 3.0_wp * nuc * kB)
        write(*, "(a, f14.6)")       "Mean Energy (eV)         : ", mean*autoev
        write(*, "(a, f14.6)")       "Mean Temperature (eV)    : ", Tcheck
        write(*, *)

        if (verbose) then
          Tcheck = stddev / (0.5_wp * 3.0_wp * nuc * kB)
          write(*, "(a, f14.6)")     "Std deviation (eV)       : ", stddev*autoev
          write(*, "(a, f14.6)")     "Std deviation (K)        : ", Tcheck
          write(*, * )
        endif

        Tcheck = maxval(ergebnis) / (0.5_wp * 3.0_wp * nuc * kB)
        write(*, "(a, f14.6)")       "Highest Energy (eV)      : ", maxval(ergebnis) &
          *autoev
        write(*, "(a, f14.6)")       "Highest Temperature (K)  : ", Tcheck
        write(*, *)
        Tcheck = minval(ergebnis) / (0.5_wp * 3.0_wp * nuc * kB)
        write(*, "(a, f14.6)")       "Lowest Energy (eV)       : ", minval(ergebnis) &
          * autoev
        write(*, "(a, f14.6)")       "Lowest Temperature (K)   : ", Tcheck
        write(*,'(41(''-''))')

        if ( minval(ergebnis) <= 0.0_wp) error stop 'energy spread too large! Reduce eimpw values and try again'

      endif



      do i = 1, ndumpGS

        ! read structure and velo from previous GS MD
        do k = 1, nuc
            read ( io_gs, *) ( xyz(j,k), j=1,3 ), ( velo(j,k), j=1,3 )
        enddo

        if ( icalc(i) == 0 ) cycle
        if ( nrun == ntraj ) exit

        nrun=nrun+1
        set%xyzr (1:3,1:nuc,nrun)=xyz (1:3,1:nuc)
        set%velor(1:3,1:nuc,nrun)=velo(1:3,1:nuc)
        set%velofr(   1:nuc,nrun)=velof(   1:nuc)

        if (.not. no_ESI) then
          set%eimpr (         nrun)= ergebnis(nrun)
        else
          set%eimpr (         nrun)= 0.0_wp
        endif
        set%taddr (         nrun)= 0.0_wp 

      enddo

      close(io_gs)
      if (.not. no_ESI) deallocate (ergebnis)

      write(*,*)
      write(*,*) nrun,' done.'
      write(*,*)

      !Take out the taddr and eimp for CID
      !taddr = 0.0d0
      !eimpr = 0.0d0

  endif iee0 !CID
  !----



  ! assume 5000 steps per trj
  write(*,'(/,'' Estimated single core run time (h/d)'', &
  &            F8.1,F8.2,/)') &
  &     (w2-w1)*5000*ntraj/3600.0_wp, &
  &     (w2-w1)*5000*ntraj/3600.0_wp/24.0_wp

  if ( nrun /= ntraj ) then
     write(*,*) nrun,ntraj
     stop 'nrun <> ntraj'
  endif

  if ( check ) stop 'normal termination of QCxMS'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  create the directories (but do checks first)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  create the directories 
  call execute_command_line('rm -rf TMPQCXMS')
  call execute_command_line('mkdir  TMPQCXMS')

  write(*,*) 'generating temporary dirs in TMPQCXMS ...'
  write(*,*)

  do i = 1, ntraj ! copy all needed information to temporary dirs
    call mdtmpdir(i)
    if ( prog == 8 ) call copyorc(i) !for xtb2 - orca setup is fine
    if ( prog == 7 ) call copyorc(i) !for xtb  - orca setup is fine
    if ( prog == 6 ) call copyorc(i) !for cal lxtb - orca setup is fine
    if ( prog == 5 ) call copyorc(i)
    if ( prog == 4 ) call copymsindo(i)
    if ( prog == 3 ) call copyorc(i)
    if ( prog == 2 ) call copytm (i)
    if ( prog /= 2 .and. iprog == 2 ) call copytm_ip(i)
    if ( prog == 1 ) call copymop(i)
    if ( prog == 0 ) call copytb (i)

    ! write coords, velo and counter to starting info file <start>
    call wrstart(i,mol,nuc, set%xyzr(1,1,i), set%velor(1,1,i), set%velofr(1,i), & 
    &        set%eimpr(i), set%taddr(i))
  enddo

  ! sum up the information at the end of creating tmp directories (call #2)
  call info_sumup(ntraj, mchrg_prod, tstep, tmax, Tinit, trelax, eimp0, &
    & ieeatm, iee_a, iee_b, ELAB, ECOM, ESI, tempESI,  & 
    & nfragexit, iprog, nuc, velo, mass)


  write(*,*)
  call version(2)
  write(*,*)
  stop ' --- normal termination of QCxMS --- '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END QCxMS set-up call (2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!           do production runs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else !prun - if production run == .true.

  ! check whether methods are suitable with codes in input
  if ( method == 0 .and. mchrg < 0 ) then
    if ( prog == 0 ) then
      stop 'Exit: DEA is not suitable with DFTB'
    elseif ( prog == 1 ) then
      stop 'Exit: DEA is not suitable with MOPAC'
    elseif ( prog == 4 ) then
       stop 'Exit: DEA is not suitable with MSINDO'
    elseif ( prog == 5 ) then
       stop 'Exit: DEA is not suitable with MNDO99'
    elseif ( prog == 6 ) then
       stop 'Exit: DEA is not suitable with (call) XTB'
    elseif ( prog == 7 ) then
       !stop 'Exit: DEA is not suitable with XTB'
    elseif ( prog == 8 ) then
       !stop 'Exit: DEA is not suitable with XTB2'
    elseif ( prog == 2 ) then
       stop 'Exit: DEA is not suitable with TURBOMOLE'
    endif

  ! CID
  elseif ( method ==  3 ) then 
     if ( prog == 0 ) then
        stop 'Exit: CID is not suitable with DFTB'
     elseif ( prog == 5 .and. gas%Iatom /= 6 .and. gas%Iatom /= 0 ) then !MNDO99 
        stop 'Exit: CID with MNDO99 only with N2 or He'
     !elseif (prog == 4) then
     !   stop 'Exit: CID is not suitable with MSINDO'
     endif
  endif

  ! save the coord file (using TMOL) into coord.original, because it will be overwritten
  if (iprog == 2) call execute_command_line('cp coord coord.original')

  ttime = 0
  nmax0 = nmax
  icoll = 0
  isec  = 0
  asave ='NOT USED'
  tcont = 0
  num_frags = 0
  frag_counter = 0
  save_counter = 0

  call timing(t1,w1)

  !> if the program did not rin before, there is no "ready" file 
  !> (see end of code), so for automatic runs, it checks if running in 
  !> normal mode is wanted.
  ex = .false.
  inquire(file='ready', exist=ex)
  !if (.not. ex ) then
    !> Read the qcxms.start file
    if (method < 3)                      call rdstart(itrj,nuc,velo,velof,tadd,eimp)
    if (method == 3 .and.       TempRun) call rdstart(itrj,nuc,velo,velof,tadd,eimp)
    if (method == 3 .and. .not. TempRun) call rdstart(itrj,nuc,velo,velof,tadd,eimp)
  !else
  !  !>> in this case, the run conditions have to come from other settings
  !  write(*,*) ' - Not using qcxms.start conditions, but settings in qcxms.in -'
  !  !>>> For EI I do not see the need for this 
  !  if (method == 0 .or. method == 1) stop 'Not useful for EI'
  !endif


  if ( method == 3 ) then 
    res_name = 'qcxms_cid.res'
  else
    res_name = 'qcxms.res'
  endif

  open(file = res_name, newunit = io_res, status='replace')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!      !!!!!!!!!!                  !!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!      !!!!!!!!!!                  !!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!      !!!!!!!!!!                  !!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!      !!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!      !!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!      !!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!!!!      !!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!!!!      !!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!!!!      !!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!!!!      !!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!!!!      !!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!!!!      !!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!!!!      !!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!!!!      !!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!!!!      !!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!      !!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!      !!!!!!!!
!!!!!!!                   !!!!!!!!!!      !!!!!!!!!!      !!!!!!!!!      !!!!!!!!
!!!!!!!                   !!!!!!!!!!      !!!!!!!!!!                  !!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!      !!!!!!!!!!                  !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize QM code for production runs
! The cid module is called dependend on the number of collisions it
! undergoes.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Go into CID module
mCID:if ( method == 3 ) then 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !set%tadd  = tadd !should be 0 anyway.
    velof = 1.0_wp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> set charge 
    mchrg =  mchrg_prod

    chrgcont  = real(mchrg,wp) 

    collisions = 0

    edum = 0

    !> Check if input (qcxms.in) makes sense
    call cidcheck(coll%max_coll, CollNo, CollSec)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Do ESI MD if ESI larger than 0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
noESI: if (.not. No_ESI )then

      write(*,*) ''
      write(*,*) '# ## ### ## # ## ### ## # ##'
      write(*,*) '# ESI simulation turned ON #'
      write(*,*) '# ## ### ## # ## ### ## # ##'
      write(*,*) ''

      !> get internal energy
      call ekinet(nuc,velo,mass,edum,t)
      edum = t * (0.5 * 3 * nuc *kB )
      !>> starting energy
      ENe(1) = edum * autoev

      !>> scaling energy
      ENe(2) = eimp * autoev

      write(*,'('' Randomly Scaling internal energy: '')')
      write(*,'(50(''-''))') 
      write(*,'('' Inner Energy before          : '', f14.6, a3)') ENe(1), ' eV'
      write(*,'('' Temperature  before          : '', f14.6, a3)')      T, ' K'
      write(*,'(50(''-''))') 
      write(*,'('' Wanted value Energy          : '', f7.4, a3)'  ) Ene(2), ' eV'
      write(*,'(50(''-''))') 


      !>> scale depending on random number (or not, if value too low)
      !edum = rand_int - edum
      edum = ENe(2) - ENe(1)
      if (edum > 0) then
        E_Scale = Ene(2)
      else
        write(*,*) ' ! No Scaling ! '
        E_Scale = 0.0_wp
      endif

      !> stop if both values were provided (to circumvent errors in input)
      if ( tempESI > 0 .and. ESI > 0 ) stop 'Cannot provide both, Energy and Temp.!'

      !> make sure nothing strange is scaled
      if (E_Scale <= 0) E_Scale = 0

      !> convert values for output
      tScale = (E_Scale * 2.0d0/3.0d0) /(nuc * kB * autoev)
      ENe(4) = tScale - t

      if (tScale > 0) then
        write(*,'('' Scaling to Temperature       : '', f14.6,a3)') tScale,' K'
        write(*,*) ' '
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Do the ESI pre-MD to scale to the energy given
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
doESI:  if ( E_Scale > 0 ) then

ESI_loop: do
          isec = isec + 1

          write(*,*)'Trajectory: ', itrj
          write(*,'(/,80(''=''))')
          write(*,*)'Heating  trajectory            ',isec
          write(*,'(80(''=''),/)')

          if ( verbose ) then
            write(*,*)'initial Cartesian coordinates:'

            do i=1,nuc
              write(*,'(i3,3f10.5,4x,a2,f10.3)')&
              &  i,xyz(1,i),xyz(2,i),xyz(3,i),to_symbol(iat(i)),mass(i)/amutoau
            enddo
          endif

          write(*,'(/,'' statistical charge  = '',F10.4,/)')chrgcont

          ! HS-UHF ini for closed-shells allowed
          mspin=0
          call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)

          ! a second attempt if this fails
          if ( .not. iniok ) then
            iniok=.false.
            call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)
            if(.not.iniok) stop 'fatal QC error. Must stop!'
          endif

          if (.not. TempRun)then

            if (isec == 1)then
              nmax = nint(ENe(4) * 2.5 )
              !prestep = prestep / (4*nuc/10)  !200
              pretadd = float(nmax)*fstoau *(tstep*autofs) * 0.75_wp
              starting_MD = .true.
            else
              pretadd = 0.0_wp
              starting_MD = .false.
              if(fragstate == 2) nmax = nmax0 * 0.75
            endif

          elseif (TempRun .and. isec == 1) then
            pretadd = float(nmax)*fstoau *(tstep*autofs) * 0.75_wp
            starting_MD = .true.
          elseif (TempRun .and. isec > 1) then
            starting_MD = .false.
            pretadd = 0.0_wp
            if(fragstate == 2) nmax = nmax0 * 0.75
          endif

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! Pre-CID MD Loop => ESI MD
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call md(itrj,icoll,isec,nuc,nmax,xyz,iat,mass,imass,mchrg,     &
          & grad,velo,velof,list,tstep,j,nfragexit,fragm,fragf,             &
          & fragat,dumpstep,etemp_in,md_ok,atm_charge,spin,axyz,tscale,pretadd,&
          & E_Scale*evtoau, .false.,Tav,Epav,Ekav,ttime,aTlast,fragstate,dtime,      &
          & ECP,starting_MD,0.0d0)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          write(*,'(/10x,''Results'',/,'' average Ekin '',F12.6)')Ekav
          write(*,'('' average Epot  '',F12.6)')Epav
          write(*,'('' average Etot  '',F12.6)')Epav+Ekav
          write(*,'('' average T (acc.)'',F12.1)')Tav
          !write(*,'('' Temp without acc '',F12.1)')new_Temp
          write(*,'('' average last T'',F12.1)')aTlast
          call ekinet(nuc,velo,mass,edum,t)
          ENe(3) = (t*(0.5*3*nuc*kB))/evtoau
          write(*,'('' final inner E'',F12.6)')ENe(3)


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! printout and count
          if(md_ok)then

            !> do IP calc and write out fragment files
            call manage_fragments(nuc, iat, xyz, axyz, mchrg, atm_charge, spin, mass, &
              &  imass, iprog, aTlast, itrj, icoll, isec, list, chrgcont, nfrag_ok,&
              &  tcont, nfrag, nometal, ECP, btf, maxsec, dtime, asave, io_res )

            write(*,*)

          !If MD fails
          else
            write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*)'something went wrong! Dont worry, the run'
            write(*,*)'is just not further counted.'
            write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            !if(index(asave,'NOT USED') == 0)write(42,'(a)')asave
            !fail = .true.
            !exit
            stop
          endif

          if (.not. nfrag_ok) exit ESI_loop !exit if too many frags


          ! Save the new coordinates
          k=0
          dum=0
          cema=0
          if (tcont > 0) then
            do i = 1, nuc
              j = list(i)
              if(j == tcont) then
                k = k + 1     
                xyzn (1:3,k) = xyz (1:3,i)
                velon(1:3,k) = velo(1:3,i)
                iatn (    k) = iat (    i)  
                massn(    k) = mass(    i)  
                imassn(   k) = imass(   i)  
                dum          = dum + iatn(k)
                cema(1:3)    = cema(1:3) + xyzn(1:3,k) * iatn(k)
              endif   
            enddo

            ! reset number of atoms
            nuc=k
            ! count number of fragmentations
            frag_counter = frag_counter + 1
            write(*,'(40(''!''))')
            write(*,'(''!! Fragmentation in ESI MD !!'')')
            write(*,'(''-- No of overall fragmentations: '',i3, '' --'')')&
              frag_counter
            write(*,'(40(''!''))')
            write(*,*)

          elseif (tcont == 0) then
            do i = 1, nuc
              xyzn (1:3,i) = xyz (1:3,i)
              velon(1:3,i) = velo(1:3,i)
              iatn (    i) = iat (    i)  
              massn(    i) = mass(    i)  
              imassn(   i) = imass(   i)  
              dum          = dum + iatn(i)
              cema(1:3)    = cema(1:3) + xyzn(1:3,i) * iatn(i)
            enddo
          endif

          ! move to Center-of-mass
          cema(1:3) = cema(1:3) / dum


          ! reduce the simulation time
          if ( nuc < 12 ) nmax = nmax0 / 4 

          ! do not continue with small fragments
          if ( nuc <= 7 ) then
            small = .true.
            exit
          endif

          ! do not continue with low masses/resolution of instrument (user)
          !if ( sum(mass(1:nuc)) / amutoau <=  minmass ) then
          !  littlemass = .true.
          !  exit
          !endif

          ! reallocate the variables, as they change for smaller systemsizes
          deallocate(xyz,axyz,grad,velo,velof,atm_charge,spin,iat,list,mass,imass)

          allocate(xyz (3,nuc),        &
                   axyz(3,nuc),        &
                   grad(3,nuc),        &
                   velo(3,nuc),        &
                   velof(nuc),         &
                   atm_charge(nuc),    &
                   spin(nuc),          &
                   iat (nuc),          &
                   list(nuc),          &
                   imass(nuc),         &
                   mass(nuc))

          do i=1,3
            xyz(i,1:nuc) = xyzn (i,1:nuc) - cema(i)
          enddo

          velo(1:3,1:nuc) = velon(1:3,1:nuc)
          iat (    1:nuc) = iatn (    1:nuc)
          mass(    1:nuc) = massn(    1:nuc)
          imass(   1:nuc) = imassn(   1:nuc)


          ! Check fragmentation state and chose what to do 
          if ( tcont > 0 ) then
            cycle ESI_loop ! CYCLE the MD loop if fragmentation happend
          else
            exit ESI_loop  ! EXIT the MD loop if nothing happend
          endif

        enddo ESI_loop ! ENDDO loop the MD module

        if ( small .or. littlemass ) then
          if( index(asave,'NOT USED') == 0 ) then
            write(io_res,'(a)')asave
          endif
        endif

      endif doESI ! ENDIF E_Scale gt 0

    endif noESI ! ENDIF No_ESI

    tadd  = 0.0_wp !should be 0 anyway.
    eimp  = 0.0_wp !should be 0 anyway.
    velof = 1.0_wp

    !Give info for not ESI runs and then no Coll. either
    if (.not. TempRun .and. small ) write(*,*) &
      ' Simulation stopped - too small molecule for collisions '

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Initial caclulation of number of collisions 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cnt:  if (.not. TempRun .and. .not. small) then

      starting_md = .false.
      new_counter = 0

      ! Full Auto = Total Program control depending on mol. size, pressure cham. length etc.
auto:   if ( FullAuto )then

        write(*,*) '--- No. of collisions automatically determined --'

        ! Calculate radius as distance between COM and Atom that is the farest away from the COM
        call collision_setup(nuc,iat,xyz,mass,rtot,cross,mfpath,calc_collisions)

        ! distribute the no. of coll. to randomize a bit
        collisions = vary_collisions(calc_collisions)

        !write(*,'(/,80(''-''))')

      ! Coll Auto = Other run modes
      elseif ( CollAuto ) then
        collisions = coll%set_coll ! starting coll are user set (default = 10)
        !write(*,'(/,80(''-''))')
        !write(*,*) 'Max. No of collisions : ', collisions

      elseif ( Manual ) then

        ! Users choose the fragmentation amount for all mod(itrj,x) runs
        ! The entire spectrum is pieced together from 3 different amount of fragmentations
noauto:   if ( CollSec(1) /= 0 ) then
          write(*,'(/,80(''-''))')
          write(*,*) '!!! Number of fragmentations are user set !!!'
          if (mod(itrj,20)  ==  0)then !all 20 runs are stopped after CollSec(3) fragmentations
             collisions = coll%set_coll
             new_counter = CollSec(3)
             write(*,*) ' - Fragmentations this run: ', new_counter
          elseif (mod(itrj,3)  ==  0)then !all 3 runs are stopped after CollSec(3) fragmentations
             collisions = coll%set_coll
             new_counter = CollSec(2)
             write(*,*) ' - Fragmentations this run: ', new_counter
          else
             collisions = coll%set_coll
             new_counter = CollSec(1) ! The rest is stopped after CollSec(1) fragmentation
             write(*,*) ' - Fragmentations this run: ', new_counter
          endif

        ! Users chose the amount of collisions for different amount (percent) of runs
        ! The entire spectrum is pieced together from 3 different amounts of collisions
        elseif(CollNo(1) /= 0)then
          write(*,'(/,80(''-''))')
          write(*,*) '!!! Amount of entire collisions are set !!!'
          write(*,*) '!!! i.e. M+ AND fragment-gas-coll (fgc) !!!'

          if     (mod(itrj,10)  ==  0) then ! every 10th run has CollNo(3) number of collisions
            collisions  = CollNo(3)
          elseif (mod(itrj,3)  ==  0) then ! every 3rd run has CollNo(2) number of collisions
            collisions  = CollNo(2)
          else
            collisions  = CollNo(1) ! The rest is stopped after CollNo(1) collisions
          endif

          new_counter = 0 ! Set to 0 so exit condition will be skipped

        ! Users chose the maximum amount of collisions between M+ and Gas
        ! There will be NO (fgc) collisions as soon as a fragmentation occurs !!!
        elseif ( coll%max_coll /= 0 ) then
          write(*,'(/,80(''-''))')
          write(*,*) '!!! M+ collisions are user set !!!'
          collisions  = coll%max_coll
          new_counter = 1

        endif noauto

      endif auto


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start Loop for CID collision and subsequent MD simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cidlp:  do
        if(icoll /= 0)then
          write(*,*) ' '
          write(*,*) ' '
        endif

        isec      = 1
        icoll     = icoll + 1
        fragstate = 0


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! First collision
        write(*,*)
        write(*,'(80(''=''))')
        write(*,'(a10,2x,i4,6x,a9,2x,i2,1x,a1,i2)')&
          & 'trajectory ',itrj, 'collision ',icoll,'/',collisions
        write(*,'(80(''-''))')
        write(*,*)
        write(*,'('' total charge          : '',i3)')       mchrg
        write(*,'('' statistical charge    : '',f8.5)')  chrgcont
        write(*,'(80(''=''))')

        ! If verbose, write Coordinates
        if ( verbose ) then
          write(*,*)'initial Cartesian coordinates:'

          do i = 1, nuc
            write(*,'(i3,3f10.5,4x,a2,f10.3)')&
            & i,xyz(1,i),xyz(2,i),xyz(3,i),to_symbol(iat(i)),mass(i)/amutoau
          enddo
        endif


        ! set the direction for the CID module after the first coll
        if ( icoll /=  1 ) then
          direc = cm2(:) - cm1(:)
          direc = direc/norm2(direc)
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Call CID module
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call cid(nuc, iat, mass, xyz, velo, tstep, mchrg, etemp_in,    &
        & stopcid, ELAB, ECOM, axyz, ttime, eExact, ECP, manual_dist,  &
        & vScale, MinPot, ConstVelo, cross,                            &
        & mfpath, rtot, atm_charge, icoll, collisions, direc, new_velo,      &
        & aTlast, calc_collisions, imass)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        new_temp = aTlast

        ! Grad failure (too big steps)
        if ( stopcid ) write(*,*) 'Error occured in the CID module'

        ! END STUFF HERE
        if ( stopcid ) exit

        !> do IP calc and write out fragment files
        call manage_fragments(nuc, iat, xyz, axyz, mchrg, atm_charge, spin, mass, &
          &  imass, iprog, aTlast, itrj, icoll, isec, list, chrgcont, nfrag_ok, &
          &  tcont, nfrag, nometal,  ECP, btf, maxsec, dtime, asave, io_res)

        if (.not. nfrag_ok) exit cidlp !exit if too many frags

        ! still have to figure out how to calculate the avg. xyz coords
         !write(*,*) 'AVG XYZ'
         !write(*,*) nuc
         !write(*,*) 
         !do i = 1,nuc
         !   write(*,*) trim(asym(iat(i))) &
         !        ,' ',axyz(1,i)* autoaa     &
         !        ,' ',axyz(2,i)* autoaa     &
         !        ,' ',axyz(3,i)* autoaa 
         !end do        


        k    = 0
        dum  = 0
        cema = 0

        if ( tcont > 0 ) then
          do i = 1, nuc
            j = list(i)
            if  ( j == tcont ) then
              k = k + 1     
              xyzn (1:3,k) = xyz (1:3,i)
              velon(1:3,k) = velo(1:3,i)
              iatn (    k) = iat (    i)  
              massn(    k) = mass(    i)  
              imassn(   k) = imass(   i)  
              dum          = dum + iatn(k)
              cema(1:3)    = cema(1:3) + xyzn(1:3,k) * iatn(k)
            endif   
          enddo

          ! reset number of atoms
          nuc=k
          ! count number of fragmentations
          frag_counter = frag_counter + 1
          write(*,'(''-- No of overall fragmentations: '',i3, '' --'')')&
            frag_counter

        elseif (tcont == 0) then
          do i = 1, nuc
            xyzn (1:3,i) = xyz (1:3,i)
            velon(1:3,i) = velo(1:3,i)
            iatn (    i) = iat (    i)  
            massn(    i) = mass(    i)  
            imassn(   i) = imass(   i)  
            dum          = dum+iatn(i)
            cema(1:3)    = cema(1:3) + xyzn(1:3,i) * iatn(i)
          enddo
        endif

        ! move to Center-of-mass
        cema(1:3) = cema(1:3) / dum

        ! do not continue with small fragments
        if ( nuc <= 7 ) then
          small = .true.
          exit
        endif

        ! do not continue with low masses/resolution of instrument (user)
        if ( sum(mass(1:nuc)) / amutoau <=  minmass ) then
          littlemass = .true.
          exit
        endif


        ! reallocate the variables, as they change for smaller systemsizes
        ! and for better visibility in the MD visualization 
        deallocate(xyz,axyz,grad,velo,velof,atm_charge,spin,iat,list,mass,imass)

        allocate(xyz (3,nuc), &
        &        axyz(3,nuc), &
        &        grad(3,nuc), &
        &        velo(3,nuc), &
        &        velof(nuc),  &
        &        atm_charge(nuc),   &
        &        spin(nuc),   &
        &        iat (nuc),   &
        &        list(nuc),   &
        &        imass(nuc),  &
        &        mass(nuc))

        do i = 1, 3
          xyz(i,1:nuc) = xyzn (i,1:nuc) - cema(i)
        enddo

        velo(1:3,1:nuc) = velon(1:3,1:nuc)
        iat (    1:nuc) = iatn (    1:nuc)
        mass(    1:nuc) = massn(    1:nuc)
        imass(   1:nuc) = imassn(   1:nuc)



        if ( small ) exit
        if ( littlemass ) exit

        ! We do not cycle the CID module, because we want to calculate
        ! the mean-free-path in between collisions


        write(*,'(80(''=''))')
        write(*,'(a)') ' - Entering Mean-Free-Path simulation - '
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!        Mean-free-path MD                                       !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Loop MD simulations if fragmentation occurs until nothing happens
MFPloop:  do
          isec = isec + 1
          Tdum = 0

          ! calculate the center-of-mass and reset for correct collision sim
          call center_of_mass(nuc,mass,xyz,cm)
          cm1(:) = cm(:) 

          !write(*,'(80(''-''),/)')
          !write(*,'('' MFP trajectory               : '',i2)') isec - 1
          write(*,'(/,80(''-''))')
          write(*,'(a,i4,6x,a9,2x,i2,1x,a1,i2,a10,1x,i4)') &
            & 'Run #',itrj, 'collision ',icoll,'/',collisions, 'MFP traj.', isec -1
          write(*,'('' total charge                 : '',i3)')mchrg
          write(*,'('' statistical charge           : '',F10.4)')chrgcont
          write(*,*)
          write(*,'(80(''=''),/)')
          if ( verbose ) then
            write(*,'(''initial Cartesian coordinates :'')')
            write(*,*)

            do i = 1, nuc
              write(*,'(i3,3f10.5,4x,a2,f10.3)')i,xyz(1,i),xyz(2,i),xyz(3,i), &
              &    to_symbol(iat(i)) , mass(i) * autoamu
            enddo
          endif


          ! HS-UHF ini for closed-shells allowed
          mspin=0

          call iniqm(nuc,xyz,iat,mchrg,mspin,etemp_in,edum,iniok,ECP)

          ! a second attempt if this fails
          if ( .not. iniok ) then
            iniok=.false.
            call iniqm(nuc,xyz,iat,mchrg,mspin,etemp_in,edum,iniok,ECP)
            if ( .not. iniok ) stop 'fatal QC error. Must stop!'
          endif


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !> reduce the simulation timings for performance

          !> change MFP times to reduce timings (empirical values)
          nmax = nuc * 100

          !> reduce the MD time if fragmentation in MFP occurs
          !> even if manually set
          if (isec == 3) nmax =int(nmax * 0.75_wp)
          if (isec == 4) nmax =int(nmax * 0.6_wp )
          if (isec >= 5) nmax =int(nmax * 0.5_wp )

          !>> not too short/long simulations
          if ( nmax < 1000   ) nmax = 1000
          if ( nmax > 10000  ) nmax = 10000

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! Do Mean-free-path (MFP) MD with nmax timesteps
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call md(itrj,icoll,isec,nuc,nmax,xyz,iat,mass,imass,mchrg,grad,&
          &       velo,velof,list,tstep,j,nfragexit,                      &
          &       fragm,fragf,fragat,dumpstep,etemp_in,                   &
          &       md_ok,atm_charge,spin,axyz,                             &
          &       Tdum,tadd,eimp,.false.,Tav,Epav,Ekav,ttime,aTlast,      &
          &       fragstate,dtime,ECP,.false.,new_velo)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !> print some energy results if wanted
          if ( verbose ) then
            write(*,'(/10x,''Results'')')
            write(*,'(20(''-''))')
            write(*,'('' average Ekin     '',F12.6)')Ekav
            write(*,'('' average Epot     '',F12.6)')Epav
            write(*,'('' average Etot     '',F12.6)')Epav+Ekav
            write(*,'('' average T (acc.) '',F12.1)')Tav
            write(*,'('' Temp without acc '',F12.1)')new_Temp
            write(*,'('' average last T   '',F12.1)')aTlast
          endif !verbose

          ! calculate the new center-of-mass as reference
          call center_of_mass(nuc,mass,xyz,cm)
          cm2(:) = cm(:)
          !!!

          ! cm_out = sqrt(cm2(1)**2+cm2(2)**2+cm2(3)**2)
          !.       - sqrt(cm1(1)**2+cm1(2)**2+cm1(3)**2)
          !
          ! new_velo = (cm_out / (ttime*fstoau)) / mstoau
          ! new_velo = (cm_out / (2000*fstoau)) / mstoau
          ! ttime=0 ! Set the md time to 0 for correct simulation time (maybe change)


          ! printout and count
          if(md_ok)then
            !> do IP calc and write out fragment files
            call manage_fragments(nuc, iat, xyz, axyz, mchrg, atm_charge, spin, mass, &
              &  imass, iprog, aTlast, itrj, icoll, isec, list, chrgcont, nfrag_ok, &
              &  tcont, nfrag, nometal, ECP, btf, maxsec, dtime, asave, io_res )

             !If MD fails
            write(*,*)
          else
            write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*)'something went wrong! Dont worry, the run'
            write(*,*)'is just not further counted.'
            write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            !stop
            exit cidlp 
          endif

          if (.not. nfrag_ok) exit cidlp !exit if too many frags

          ! Save the new coordinates
          k    = 0
          dum  = 0
          cema = 0

          if ( tcont > 0 ) then
            do i = 1, nuc
              j = list(i)
              if ( j == tcont ) then
                k = k + 1     
                xyzn (1:3,k) = xyz (1:3,i)
                velon(1:3,k) = velo(1:3,i)
                iatn (    k) = iat (    i)  
                massn(    k) = mass(    i)  
                imassn(   k) = imass(   i)  
                dum          = dum+iatn(k)
                cema(1:3)    = cema(1:3) + xyzn(1:3,k) * iatn(k)
              endif   
            enddo

            ! reset number of atoms
            nuc = k
            ! count number of fragmentations
            frag_counter = frag_counter + 1
            write(*,'(''-- No of overall fragmentations: '',i3, '' --'')')&
              frag_counter

          elseif ( tcont == 0 ) then
            do i = 1, nuc
              xyzn (1:3,i) = xyz (1:3,i)
              velon(1:3,i) = velo(1:3,i)
              iatn (    i) = iat (    i)  
              massn(    i) = mass(    i)  
              imassn(   i) = imass(   i)  
              dum          = dum+iatn(i)
              cema(1:3)    = cema(1:3)+xyzn(1:3,i)*iatn(i)
            enddo
          endif

          ! move to Center-of-mass
          cema(1:3) = cema(1:3) / dum

          ! do not continue with small fragments
          if ( nuc <= 7 ) then
            small = .true.
            exit
          endif

          ! do not continue with low masses/resolution of instrument (user)
          if ( sum(mass(1:nuc)) / amutoau <=  minmass ) then
            littlemass = .true.
            exit
          endif

          ! reallocate the variables, as they change for smaller systemsizes
          deallocate(xyz,axyz,grad,velo,velof,atm_charge,spin,iat,list,mass,imass)

          allocate(xyz (3,nuc), &
          & axyz(3,nuc),        &
          & grad(3,nuc),        &
          & velo(3,nuc),        &
          & velof(nuc),         &
          & atm_charge(nuc),    &
          & spin(nuc),          &
          & iat (nuc),          &
          & list(nuc),          &
          & imass(nuc),         &
          & mass(nuc))

          do i = 1, 3
            xyz(i, 1:nuc) = xyzn (i,1:nuc) - cema(i)
          enddo

          velo(1:3,1:nuc) = velon(1:3,1:nuc)
          iat (    1:nuc) = iatn (    1:nuc)
          mass(    1:nuc) = massn(    1:nuc)
          imass(   1:nuc) = imassn(   1:nuc)


          ! chose what to do
          if ( tcont > 0 ) then
            cycle MFPloop ! CYCLE the MD loop if fragmentation happend
          else
            exit MFPloop  ! EXIT the MD loop if nothing happend
          endif

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          enddo MFPloop ! ENDDO loop the MD module
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Define some stopping criteria for the collisions
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! 1.) If the fragments become too small
          if ( small ) exit
          if ( littlemass ) exit


          ! 2.) If the velocity or the !COM Energy! get too low
          summass = 0

          do i = 1, nuc
            summass = summass +  mass(i)
          enddo

          E_KIN = 0.5_wp * summass * (( new_velo * mstoau )**2)
          beta = gas%mIatom / (gas%mIatom + summass)
          E_COM = calc_ECOM(beta,E_KIN) !(beta * E_KIN) * autoev
          !write(*,*) 'BETA      :    ',beta
          !write(*,*) 'ECOM (eV) :    ', E_COM
          !write(*,*)'NEW VELO MD:', new_velo

          ! set end-conditions for the CID module
          if ( new_velo <= 800 .or. E_COM <= 0.85_wp .and. MinPot == 0 ) then
            if(index(asave,'NOT USED') == 0)then
              write(io_res,'(a)')asave
            endif
            write(*,'(40(''!''))')
            write(*,'(''!   Low velocity and thus low E(COM)  !'')')
            write(*,'('' -> The velocity is now: '',f12.4)') new_velo
            write(*,'('' -> E(COM)       is now: '',f12.4)') E_COM
            write(*,'(40(''!''))')
            exit
          endif

          !--------------------------------------------------------------
          ! 3.1.) If the number of fragmentaion steps are exceeded
          !Collauto  == false
          if ( Manual .and. new_counter > 0 ) then
            if ( frag_counter >= new_counter ) then
              if ( index(asave,'NOT USED') == 0) then
                write(io_res,'(a)')asave
              endif
              write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              write(*,*) ' Maximum number of fragmentations done. Exit.'
              write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              exit !Make this value depnendent on first nmbr of nucs
            endif
          endif

          !--------------------------------------------------------------
          ! 3.2.) If fragmentation occured, check the rest of the collisions
          !       Fullauto run-type!
Full:     if (FullAuto .and. frag_counter > save_counter) then

            save_counter = frag_counter

            ! Set-up collision number and vary the amount (cid.f90)
            call collision_setup(nuc,iat,xyz,mass,rtot,cross,mfpath, &
            &    calc_collisions)

            ! set varied number of collisions by BoxMuller distr.
            collisions = vary_collisions(calc_collisions)

            ! Check if molecule is not too small to make collisions
            if ( collisions == 0 .and. icoll /= 1 ) then
              write(*,*) 'Molecule too small for further collisions'

              if(index(asave,'NOT USED') == 0)then
                write(io_res,'(a)')asave
              endif

              exit ! Finished 

            endif

            !write(*,*) ' Radius Ion (m)               : ', rtot
            !write(*,*) ' Radius Gas Atom (m)          : ', (rIatom/aatoau) * 1E-10 
            write(*,*) ' Cross-section (M+) (m²)      : ', cross
            write(*,*) ' Mean free path (m)           : ', mfpath
            write(*,*) ' chamber length (m)           : ', cell%lchamb
            !write(*,*) ' Number of collisions varied  : ', collisions

          endif Full

          !--------------------------------------------------------------
          ! 3.3) For Collauto (NOT Fullauto)
          ! Vary the amount of collisions depending on the number of atoms
CollLp:   if (CollAuto .and. frag_counter > save_counter) then

            ! Set-up collision number and vary the amount (cid.f90)
            call collision_setup(nuc,iat,xyz,mass,rtot,cross,mfpath, &
            &    calc_collisions)

            if ( collisions > 0 )then
              !> set random seed and number 
              !> else the random seed from start is taken
              if (iseed(1) == 0) then
                !call random_seed()
                call random_seed(numb)
              endif

              call random_number(a)

              b = nuc / 10.0
              dep = nint(b)

              rand_int =  FLOOR((dep+1)*a) ! 0-dep
              write(*,'('' Random integer        : '',I2)') rand_int

              collisions = icoll + rand_int !+ dep
              numb=numb+1
              write(*,'('' Total no. collisions  : '',i2)') collisions

            elseif( collisions == 0 .and. icoll /= 1)then
              write(*,*) 'Molecule too small for further collisions'

              if(index(asave,'NOT USED') == 0)then
                write(io_res,'(a)')asave
              endif

              exit ! Finished 

            endif
          endif CollLp

          ! 4.) If the max number of collisions is reached
          if ( icoll >= collisions ) then
            write(*,*)'-------------------------------------------------'
            write(*,*)' Last MD is finished. All is written out'
            write(*,*)'-------------------------------------------------'
            write(io_res,'(a)')asave

            if(num_frags == 0)then
               write(*,*)''
               write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(*,*)'    No fragmentation in the simulation!   '
               write(*,*)'   Increase energy or time of sampling.   '
               write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            endif

            exit cidlp ! End the collision routine

          elseif ( icoll /= collisions .and. tcont == 0 ) then

            write(*,*)''
            write(*,*)'--------------------------------------------------'
            write(*,*)'-- NOTHING HAPPENED - Next collision simulation --'
            write(*,*)'--------------------------------------------------'

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo cidlp! ENDDO Collision routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Write out too small molecules
        if(small)then
           if(index(asave,'NOT USED') == 0)then
              write(io_res,'(a)')asave
           endif
           write(*,'(60(''!''))')
           write(*,*)'     run not continued because of small fragment(s)'
           write(*,'(60(''!''))')
           write(*,*)''
        endif

        ! Write out too light molecules
        if(littlemass)then
           if(index(asave,'NOT USED') == 0)then
              write(io_res,'(a)')asave
           endif
           write(*,'(60(''!''))')
           write(*,'(''run not continued because of resolution'')')
           write(*,'(60(''!''))')
           write(*,'(''Threshold  : '',4x,i4)') minmass
           write(*,'(''Frag. mass : '',4x,f10.6)') sum(mass(1:nuc)) * autoamu
           write(*,*)
        endif

        ! ERROR
        if(stopcid)then
           if(index(asave,'NOT USED') == 0)write(io_res,'(a)')asave
           write(*,'(60(''!''))')
           write(*,*)' --- run aborted, last structure saved! --- '
           write(*,'(60(''!''))')
           write(*,*)
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif cnt ! ENDIF not temprun 
    endif mCID ! ENDIF CID MODULE (Method == 3)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start Electron ionization Routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   !!!!!!!!!!!!!!!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start Electron ionization Routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ei: if (method /= 3 ) then !.and. method /= 4)then

      !> set the starting overall charge to user input (or default)
      mchrg =  mchrg_prod

      !> set the statistical charge to user input (or default)
      chrgcont  = real(mchrg,wp) 

loop: do
        !> increase the count of runs
        isec = isec + 1

        !> write out info
        write(*,'(/,80(''=''))')
        write(*,*)'                      trajectory ',itrj,isec
        write(*,'(80(''=''),/)')
        if ( verbose ) then
          write(*,*)'initial Cartesian coordinates:'

          !> write out coordinates
          do i = 1, nuc
             write(*,'(i3,3f10.5,4x,a2,f10.3)')&
             & i,xyz(1,i),xyz(2,i),xyz(3,i),to_symbol(iat(i)),mass(i)/amutoau
          enddo
        endif

        !> write out statistical charge
        write(*,'(/,'' statistical charge  = '',F10.4,/)')chrgcont

        ! re- initialize the GBSA module (if ever becomes important)
        !         if (lgbsa)then
        !           call ekinet(nuc,velo,mass,edum,t)
        !           call init_gbsa(nuc,iat,solvent,gsolvstate,t)
        !         endif
       
        !> initialize the QC code
        !> HS- UHF ini for closed-shells allowed
        mspin = 0
        call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)

        !> a second attempt if this fails
        if ( .not. iniok ) then
           iniok =.false.
           call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)
           if ( .not. iniok) stop 'fatal QC error. Must stop!'
        endif

        !> if the loop runs more times (in EI i.e. consecutive 
        !  fragmentation, set nfragexit = 2
        if ( isec > 1 ) nfragexit = 2

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>  do production MD run
        Tdum=0
        call md(itrj,0,isec,nuc,nmax,xyz,iat,mass,imass,mchrg,grad, &
        &       velo,velof,list,tstep,j,nfragexit,                  &
        &       fragm,fragf,fragat,dumpstep,etemp_in,                &
        &       md_ok,atm_charge,spin,axyz,                                &
        &       Tdum,tadd,eimp,.false.,Tav,Epav,Ekav,ttime,aTlast,  &
        &       fragstate,dtime,ECP,.false.,0.0_wp)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !> write out the energies
        write(*,'(/10x,''Results'',/,'' average Ekin '',F12.6)')Ekav
        write(*,'('' average Epot  '',F12.6)')Epav
        write(*,'('' average Etot  '',F12.6)')Epav+Ekav
        write(*,'('' average T     '',F12.1)')Tav
        write(*,'('' average last T'',F12.1)')aTlast

        ! printout and count
okmd:   if(md_ok)then
          !> do IP calc and write out fragment files
          call manage_fragments(nuc, iat, xyz, axyz, mchrg, atm_charge, spin, mass, imass, &
          &  iprog, aTlast, itrj, icoll, isec, list, chrgcont, nfrag_ok,  &
          &  tcont, nfrag, nometal, ECP, btf, maxsec, dtime, asave, io_res )

          if (.not. nfrag_ok) exit !exit if too many frags

          ! this is likely to be wrong (too many fragmnets for EI)
          !if ( nfrag > 5 ) then
          !  write(*,*)'    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          !  write(*,*)'    !! something went wrong! Dont worry, the run !!'
          !  write(*,*)'    !! is just not counted.                      !!'
          !  write(*,*)'    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          !  exit
          !endif

frag:     if(tcont > 0 .and. maxsec > 0)then
            ! because IEE decreases, we can reduce the max length of
            ! the trj by 20% for each generation
            ! but not less than nmax(start)/2
            nmax=nmax-nmax0/5
            nmax=max(nmax,nmax0/2)
            ! also if frag is small (see below)

            ! if nothing happened already in the last MD, make it short
            if(fragstate == 2) nmax=nmax0/3

            k=0
            dum=0
            cema=0
            do i=1,nuc
               j=list(i)
               if(j == tcont)then
                  k=k+1
                  xyzn (1:3,k)=xyz (1:3,i)
                  velon(1:3,k)=velo(1:3,i)
                  iatn (    k)=iat (    i)
                  massn(    k)=mass(    i)
                  imassn(   k)=imass(   i)
                  dum         =dum+iatn(k)
                  cema(1:3)   =cema(1:3)+xyzn(1:3,k)*iatn(k)
               endif
            enddo
            ! move to center-of-mass (CEMA)
            cema(1:3) = cema(1:3) / dum
            nuc=k
            if(nuc <= 10) nmax=nmax0/10

            ! do not continue with small fragments
            if(nuc <= 3)then
               if(index(asave,'NOT USED') == 0)write(io_res,'(a)')asave
               write(*,*)'     run not continued because of small fragment(s)'
               exit
            endif

            deallocate(xyz,axyz,grad,velo,velof,atm_charge,spin,iat,list,mass,imass)

            allocate(xyz (3,nuc), &
            &  axyz(3,nuc),       &
            &  grad(3,nuc),       &
            &  velo(3,nuc),       &
            &  velof(nuc),        &
            &  atm_charge(nuc),         &
            &  spin(nuc),         &
            &  iat (nuc),         &
            &  list(nuc),         &
            &  imass(nuc),        &
            &  mass(nuc))

            do i=1,3
               xyz(i,1:nuc)=xyzn (i,1:nuc)-cema(i)
            enddo
            velo(1:3,1:nuc)=velon(1:3,1:nuc)
            iat (    1:nuc)=iatn (    1:nuc)
            mass(    1:nuc)=massn(    1:nuc)
            imass(   1:nuc)=imassn(   1:nuc)
            velof=0
            tadd=0
            eimp=0

            cycle loop !if fragmented

          else ! if No fragmentation

            exit loop ! exit the EI- MD loop

          endif frag

        else !If MD failed
          write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write(*,*)'something went wrong! Dont worry, the run'
          write(*,*)'is just not further counted.'
          write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          if(index(asave,'NOT USED') == 0)write(io_res,'(a)')asave

          exit loop

        endif okmd
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo loop  ! End the prod. run loop 
     
    endif ei ! END IF the EI/DEA Module

  endif prun ! End the production run tree

!ENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDEND

  ! clear the node even when mdok=.F. (<ready> is read by script)
  call execute_command_line('touch ready')
  close(io_res)
  close(io_log)


  call timing(t2,w2)
  write(*,'(/,'' wall time (min)'',F10.2  )')(w2-w1)/60.0_wp
  !write(*,'(  '' # of QC calls  '',I10  ,/)')calls

  call execute_command_line('date')
  call version(2)
  write(*,*)'normal termination of QCxMS'

end program QCxMS
