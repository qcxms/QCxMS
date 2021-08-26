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
  use get_version
  use qcxms_boxmuller, only: vary_collisions, vary_energies
  use qcxms_mo_spec, only: getspec
  use qcxms_fragments
  use qcxms_iee
  use qcxms_impact, only: calctrelax
  use qcxms_iniqm, only: iniqm
  use qcxms_mdinit, only: mdinitu, ekinet
  use qcxms_info, only: info_main, info_sumup, cidcheck
  use qcxms_read_coordinates
  use qcxms_use_orca, only: copyorc
  use qcxms_use_turbomole, only: copytm, copytm_ip
  use qcxms_utility
  use qcxms_write_fragments, only: manage_fragments 
  use xtb_mctc_accuracy, only : wp
  use xtb_mctc_convert
  use xtb_mctc_constants
  use xtb_mctc_symbols, only: toSymbol 

  implicit none

  integer  :: nmax,nmax0,dumpstep
  integer  :: nfrag,ndumpGS
  integer  :: ntraj,iseed(1)
  integer  :: i,j,k,m
  integer  :: nrun,nfragexit
  integer  :: mspin,mchrg,iprog,tcont
  integer  :: maxsec,idum,edistri
  integer  :: ihomo,nb,nuc,mo1,mo2,fragstate
  integer  :: num_frags
  integer  :: itrj,icoll,isec,prestep
  integer  :: scani
  integer  :: fragat(200,10)
  integer  :: imassn(10000)
  integer  :: io_res, io_gs, io_log, io_eimp
  integer  :: CollSec(3),CollNo(3)
  integer  :: collisions,set_coll,MaxColl
  integer  :: minmass,numb
  integer  :: coll_counter,new_counter,frag_counter
  integer  :: simMD,save_MD
  integer  :: rand_int,dep,convetemp
  ! GBSA
  !integer  :: gsolvstate

  integer,allocatable  :: list(:)
  integer,allocatable  :: iat (:)
  integer,allocatable  :: icalc(:)
  integer,allocatable  :: imass(:)

  real(wp) :: etemp
  real(wp) :: tstep,edum,tmax,etempin,exc,dums,betemp
  real(wp) :: t,Tsoll,Tdum,Tinit,trelax,ttime,ehomo
  real(wp) :: Epav,Ekav,Tav,dum,t1,t2,w1,w2,eimpw,pmax,tta
  real(wp) :: tadd,fimp,aTlast
  real(wp) :: chrgcont,cema(3),eimp0,eimp,dtime
  real(wp) :: iee_a,iee_b,hacc,ieeel,btf,ieeatm,etempGS
  real(wp) :: vScale,new_velo,avgT,MinPot,ESI,tempESI
  real(wp) :: tScale,new_temp
  real(wp) :: summass,beta
  real(wp) :: direc(3),cm(3),cm1(3),cm2(3)
  !real(wp) :: cm_out
  real(wp) :: randx
  real(wp) :: lowerbound,upperbound
  real(wp) :: Eimpact
  real(wp) :: gaus(1000)
  real(wp) :: a,b
  real(wp) :: rtot,cross,mfpath
  real(wp) :: E_KIN,E_COM
  real(wp) :: calc_collisions,E_Scale,ENe(10)
  ! Ion tracking arrays
  real(wp) :: fragm(10) !10 fragments max per fragmentation,
  real(wp) :: xyzn(3,10000),velon(3,10000),iatn(10000),massn(10000) !10000 atoms max


  ! Allocatables
  real(wp),allocatable :: xyz (:,:)
  real(wp),allocatable :: axyz(:,:)
  real(wp),allocatable :: xyzr(:,:,:)
  real(wp),allocatable :: grad(:,:)
  real(wp),allocatable :: velo(:,:)
  real(wp),allocatable :: velor(:,:,:)
  real(wp),allocatable :: velof(:)
  real(wp),allocatable :: eimpr(:)
  real(wp),allocatable :: taddr(:)
  real(wp),allocatable :: velofr(:,:)
  real(wp),allocatable :: mass(:)
  real(wp),allocatable :: chrg(:)
  real(wp),allocatable :: spin(:)
  real(wp),allocatable :: emo  (:)
  real(wp),allocatable :: mopop(:,:)
  real(wp),allocatable :: modum(:)

  ! Strings
  character(len=:), allocatable :: res_name
  character(len=:), allocatable :: file_name
  character(len=80)    :: fragf(10)
  character(len=120)   :: asave
  character(len=80)    :: arg(10)
  !character*80 solvent

  ! logicals
  logical :: ex,mdok,check,takempop
  logical :: prod,noeq,iniok
  logical :: eonly,eonly0,eonly1
  logical :: unity ! scales velocities uniformly 
  logical :: metal3d ! switches on some features for transition metal compounds 
  logical :: ECP !switch on ECP for orca calculations
  logical :: stopcid !if error comes up in CID module - code is killed
  logical :: noecp,nometal !forbid checking of ECP and METAL
  logical :: ConstVelo,eExact
  logical :: small,littlemass
  logical :: No_ESI,NoScale
  logical :: starting_md
  !logical gbsa ! set solvation model

  external system

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Start the program
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call system('date')
  write(*,'(//&
  &          22x,''*********************************************'')')
  write(*,'(22x,''*                                           *'')')
  write(*,'(22x,''*            Q   C   x   M   S              *'')')
  write(*,'(22x,''*                                           *'')')
  call version(0)
  write(*,'(22x,''*                                           *'')')
  write(*,'(22x,''*                S. Grimme                  *'')')
  write(*,'(22x,''* Mulliken Center for Theoretical Chemistry *'')')
  write(*,'(22x,''*             Universitaet Bonn             *'')')
  write(*,'(22x,''*                  2008-21                  *'')')
  call version(1)
  write(*,'(22x,''*                                           *'')')
  write(*,'(22x,''*********************************************'')')
  write(*,*)
  write(*,'('' QCxMS is free software: you can redistribute it and/or &
           & modify it under'')')
  write(*,'('' the terms of the GNU Lesser General Public License as &
           & published by '')') 
  write(*,'(''the Free Software Foundation, either version 3 of the License, or '')') 
  write(*,'(''(at your option) any later version.'')')
  write(*,*)
  write(*,'('' QCxMS is distributed in the hope that it will be useful, '')')
  write(*,'('' but WITHOUT ANY WARRANTY; without even the implied warranty of '')') 
  write(*,'('' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the '')')
  write(*,'('' GNU Lesser General Public License for more details.'')')
  write(*,*)
  write(*,'(''Cite this work as:'')')
  write(*,'(''S.Grimme, Angew.Chem.Int.Ed. 52 (2013) 6306-6312.'')')
  write(*,*)
  write(*,'(''for the CID module:'')')
  write(*,'(''J. Koopman, S. Grimme, J. Am. Soc. Mass Spectrom., (2021), &
           & DOI: 10.1021/jasms.1c00098 '')')
  write(*,*)
  write(*,'(''for the GFN1-xTB implementation:'')')
  write(*,'(''V. Asgeirsson, C.Bauer, S. Grimme, Chem. Sci. 8 (2017) 4879'')')
  write(*,*)
  write(*,'(''for the GFN2-xTB implementation:'')')
  write(*,'('' J. Koopman, S. Grimme, ACS Omega 4 (12) (2019) 15120-15133, &
           & DOI: 10.1021/acsomega.9b02011 '')')
  write(*,*)
  write(*,'(50(''-''))')
  write(*,'(''Current Dev. : J. Koopman '')')
  write(*,'(''Former Dev.  : V.Asgeirsson, C.Bauer '')')
  write(*,'(50(''-''))')
  write(*,*)

  ! use/uncomments this piece of code for the lfc compiler which generates other system calls
  !     call system('echo $SHELL > .tmpqcxms')
  !     open(unit=1,file='.tmpqcxms')
  !     read(1,'(a)')adum
  !     close(1,status='delete')
  !     if(index(adum,'tcsh') /= 0) shell=1
  !     if(index(adum,'bash') /= 0) shell=2


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
  noeq   =.false.
  ! test calls
  eonly  =.false.
  eonly0 =.false.
  eonly1 =.false.
  takempop=.false.
  ! uniform velocity scaling
  unity=.false.
  ! check if cid was OK
  stopcid = .false.
  starting_md=.false.
  ! HS-UHF ini only for frag runs
  iniok  =.true.
  ! dump every dumpstep MD steps for MOLDEN movie (=4 fs as default)
  dumpstep=4
  ! counts the number of QC calls
  calls=0
  ! GS Etemp (to converge radicals etc)
  etempGS=298.15 ! normal ! Maybe make this input relevant
  convetemp=0

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
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read the INPUT file
  ! input, all defaults are set here
  call input(tstep,tmax,ntraj,iseed(1),etempin,Tinit, &
  &          iee_a,iee_b,eimp0,eimpw,fimp,iprog,                            &
  &          trelax,hacc,nfragexit,maxsec,edistri,btf,ieeatm,               &
  &          scani,lowerbound,upperbound,metal3d,                    &
  &          Eimpact,eExact,ECP,unity,noecp,nometal,   &
  &          vScale,CollNo,CollSec,ConstVelo,   &
  &          minmass,simMD,convetemp,set_coll,MaxColl,          &
  &          MinPot,ESI,tempESI,No_ESI,NoScale)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Choose the MS method
  write(*,*)
  if (method  ==  0) then
     write(*,'(14x,60(''*''))')
     write(*,'(22x,'' Mode:          Electron Impact (EI)         '')')
     write(*,'(14x,60(''*''))')
  elseif (method  ==  1) then !This will be closed shell with EI
     write(*,'(14x,60(''*''))')
     write(*,'(22x,'' Mode:         Closed Shell Cation (CSC)     '')')
     write(*,'(14x,60(''*''))')
  elseif (method  ==  2) then
     write(*,'(14x,60(''*''))')
     write(*,'(22x,'' Mode: Dissociative Electron Attachment (DEA)'')')
     write(*,'(14x,60(''*''))')
  elseif (method  ==  3) then
     write(*,'(14x,60(''*''))')
     write(*,'(22x,'' Mode: Collision Induced Dissociation (CID)  '')')
     write(*,'(22x,''             Positive Ion mode               '')')
     write(*,'(14x,60(''*''))')
  elseif (method  ==  4) then
     write(*,'(14x,60(''*''))')
     write(*,'(22x,'' Mode: Collision Induced Dissociation (CID)  '')')
     write(*,'(22x,''             Negative Ion mode  '')')
     write(*,'(14x,60(''*''))')
  else
     write(*,'(22x,''**************************'')')
     write(*,'(22x,''* Error Error Error      *'')')
     write(*,'(22x,''* Method Option failed?  *'')')
     write(*,'(22x,''**************************'')')
     stop 'Something went terribly wrong!'
  endif

  ! initialize D3
  if (prog /= 8) call copyc6



  ! how many atoms? Initialize to provide setup step
  !call rd0('coord',nuc)
  call rd(.true.,file_name,nuc)

  ! set traj automatically
  if(ntraj <= 0) ntraj = nuc*25

  ! get for the QC method used the eTemp as info
  call setetemp(1,-1.0d0,betemp)

  ! ini RNDs
  call random_seed (put=iseed)
  call random_number(randx)

  ! for GS every n steps a struc. is used later for production, i.e.
  ! 100 fragmentation runs require about 5000 GS MD steps
  ! this ensures uncorrelated fragmentations
  nmax0 = ntraj * 50

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! READ command line options
  do i=1,10
     call getarg(i,arg(i))
  enddo
  do i=1,10
     if(index(arg(i),'-check') /= 0)check  =.true. !Check IEE settings
     if(index(arg(i),'-c'    ) /= 0)check  =.true. !Check IEE settings
     if(index(arg(i),'-prod' ) /= 0)prod   =.true. !Do production run
     if(index(arg(i),'-p'    ) /= 0)prod   =.true. !Do production run
     if(index(arg(i),'-noeq' ) /= 0)noeq   =.true. !Skip equilibration MD
     if(index(arg(i),'-e0'   ) /= 0)eonly0 =.true. !Only calc. energy chrg = 0
     if(index(arg(i),'-e1'   ) /= 0)eonly1 =.true. !Only calc. energy chrg = 1
     if(index(arg(i),'-eonly') /= 0)eonly  =.true. !Only calc. energy chrg = .CHRG file
     if(index(arg(i),'-v'    ) /= 0)verbose  =.true.   ! more infos 
     if(index(arg(i),'-unity') /= 0)unity  =.true. !Set velocity scaling to unity
     if(index(arg(i),'-qcp'  ) /= 0)path=arg(i+1)  !Use GCP model (not implemented)
     !if(index(arg(i),'-gbsa' )  /= 0)then
     !  lgbsa =.true.
     !  write(*,*)'Solvation activated'
     !  if(index(atmp,'-') == 0.and.atmp(1:1) /= ' ')then
     !     solvent=arg(i+1)
     !  endif
     !endif
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Test etemp for CID - hardcoded for XTB for now
  if (method == 3 .and. prod)then
     if(etempin <= 0)then
        betemp = 5000.0_wp
     else
        betemp = etempin
     endif
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! printing runtype information and chosen parameters
  call info_main(ntraj, tstep, tmax, Tinit, trelax, eimp0, &
      & ieeatm, iee_a, iee_b, btf, fimp, hacc, eimpact, MaxColl, CollNo, CollSec,  &
      & ESI, tempESI, eTempin, maxsec, betemp, nfragexit, iseed, iprog, edistri)


  ! # MD steps in a frag run
  nmax = tmax / tstep
  ! timesteps in au
  tstep = tstep * fstoau
  ! the "70" eV in a.u.
  eimp0 = eimp0 * evtoau


  ! Check if input provides reasonable molecule
  if(nuc < 1)     stop 'no reasonable molecule (coord) found!'
  if(nuc > 10000) stop 'too many atoms'

  ! allocate memory for molecular properties
  allocate(xyz (3,nuc),      &
  &        axyz(3,nuc),      &
  &        grad(3,nuc),      &
  &        velo(3,nuc),      &
  &        velof(nuc),       &
  &        emo(10*nuc),      &
  &        mopop(5*nuc,nuc), &
  &        modum(5*nuc),     &
  &        chrg(nuc),        &
  &        spin(nuc),        &
  &        iat (nuc),        &
  &        list(nuc),        &
  &        imass(nuc),       &
  &        mass(nuc))

  ! Description:
  !      nuc : Number of atoms in molecule
  !      mass(nuc) : Mass of the atom
  !      chrg(nuc) : Charge of atom
  !      iat(nuc) : typ of atom (index)
  !      ncao / nbf : number of cartesian AOs
  !      nel : number of valence electrons (-1)
  !      ndim / : numer of AOs
  !      ihomo : integer of HOMO orbital
  !      edum : Energy of the system (3/2*k_B*T*nuc)

  ! read coord in TM format
!  call rd('coord',nuc,xyz,iat)
  call rd(.false.,file_name,nuc,xyz,iat)

  ! check for 3d transition metals and 4s-4d/5s-5d elements
  if (.not. noecp)then
     do i=1,nuc
        if(iat(i) <= 86.and.iat(i) >= 37) ECP=.true.
     enddo
  endif
  if (.not. nometal)then
     do i=1,nuc
        if(iat(i) <= 30.and.iat(i) >= 22) metal3d=.true.
     enddo
  endif
  if(metal3d) then
     write(*,*)
     write(*,*)'* 3d metal found, check for multiplicities:&
     & ',metal3d,'*'
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

  if(.not.prod)then
     write(*,*)
     write(*,*)'molecule:'

     do i=1,nuc
        write(*,'(i3,3f13.6,4x,a2,f9.4)') &
        &     i,xyz(1,i),xyz(2,i),xyz(3,i),toSymbol(iat(i)),mass(i) * autoamu
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

  ! Set the charge of the chosen method 
  ! 0 = EI , 1 = CSC, 2 = DEA, 3 = CID+, 4 = CID-
  if     (method  ==  0 .or. method  ==  2) then
    mchrg = 0
  elseif (method  ==  1 .or. method  ==  3) then
    mchrg = 1
  elseif (method  ==  4) then
    mchrg = -1
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! skip if production runs are selected (parallel)
prun: if(.not.prod) then 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

      if(convetemp /= 0)then
         etempGS=convetemp
         write(*,*) 'Groundstate eTemp set to',etempGS
      endif


      ! init the QM code  ! Here edum is total energy, not inner energy anymore!
      call iniqm(nuc,xyz,iat,mchrg,mspin,etempGS,edum,iniok,ECP)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! equilibrate GS
      if(.not.iniok) stop 'fatal QC error. Must stop!'
      
      write(*,*)'M trajactory, equilibration ...',nmax0/2

         call md(-1,0,0,nuc,nmax0/2,xyz,iat,mass,imass,mchrg,grad,     &
         &       velo,velof,list,tstep,ndumpGS,99,                     &
         &       fragm,fragf,fragat,dumpstep,etempGS,                  &
         &       mdok,chrg,spin,axyz,                             &
         &       Tinit,0.0d0,0.0d0,.false.,Tav,Epav,Ekav,ttime,aTlast, &
         &       fragstate,dtime,ECP,.false.,0.0d0)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do MD GS
      ndumpGS=0

      write(*,*)'M trajactory, sampling ...',nmax0

         call md(0,0,0,nuc,nmax0,xyz,iat,mass,imass,mchrg,             &
         &       grad,velo,velof,list,tstep,ndumpGS,99,                &
         &       fragm,fragf,fragat,dumpstep,etempGS,                  &
         &       mdok,chrg,spin,axyz,                             &
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END QCxMS groundstate call (1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*)'structures for fragmentation runs written'
      call version(2)
      stop 'normal termination of QCxMS'


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! START QCxMS set-up call (2) 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else ! if qcxms.gs is NOT yet created

      ! Check if CID input makes sense before creating the dirs.
      if ( method == 3 .or. method == 4 ) then
        call cidcheck(MaxColl, CollNo, CollSec )
      endif

      write(*,*) 'reading qcxms.gs ...'

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ignored if not used (default)
      ! do MD using previous GS traject
      if ( noeq ) then
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
         &       mdok,chrg,spin,axyz,                             &
         &       Tinit,0.0d0,0.0d0,.false.,Tav,Epav,Ekav,ttime,aTlast, &
         &       fragstate,dtime,ECP,.false.,0.0d0)


         write(*,'(/,'' average Ekin '',F12.6)')Ekav
         write(*,'('' average Epot '',F12.6)')Epav
         write(*,'('' average Etot '',F12.6)')Epav+Ekav
         write(*,'('' average T    '',F12.1)')Tav
         write(*,*)'structures for fragmentation runs written'
         call system('rm -rf fort.11 fort.15')
         call system('rm -rf charges.bin detailed.outi dftb_in.hsd')
         call system('rm -rf dftb_pin.hsd band.out')
         call version(2)
         stop 'normal termination of QCxMS'
      endif
      ! ignored if not used (default)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    endif GS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Test SP energy calc. of the used QM method on the input geom with charges (+1,-1)
    ! (It is important to do for calcs with TM and ORCA)
    call timing(t1,w1)
    if     (method  ==  0 .or. method  ==  1 .or. method == 3) then !cations
      mchrg = 1
    elseif (method  ==  2 .or. method == 4) then ! anions
      mchrg = -1
    endif

    write(*,'(''--- Checking QC method for ions ---'')')
    call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)
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
    allocate(xyzr  (3,nuc,ntraj),       &
    &        velor (3,nuc,ntraj),       &
    &        velofr(  nuc,ntraj),       &
    &        eimpr(ntraj),taddr(ntraj), &
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
iee0: if (method /= 3 .and. method /= 4)then ! Not important for CID

      emo = 0

      ! mchrg = 0 to get the values for ground-state
      !call getspec(.true.,nuc,iat,xyz,0,eps,ehomo,mopop,ihomo,nb,ECP) !eps has to be changed in geteps 
      call getspec(.true.,nuc,iat,xyz,0,emo,ehomo,mopop,ihomo,nb,ECP)

      ! adjust IEE distr. if necessary (max in intervall)
      if (scani  ==  0 ) then
         write(*,'(/,'' preparing the IEE distribution ...'')')
      elseif (scani  ==  1) then
         write(*,'(/,'' IEE scanning feature  ...'')')
      endif

      exc   = (eimp0 - ehomo) * autoev
      !exc   = (eimp0 - ehomo) * autoev TEST
      ieeel = dble (ihomo + nb)

      ! user input
iee1: if ( iee_a > 0 .and. iee_b > 0 ) then
        call getmaxiee(iee_a,iee_b,ieeel,edistri,exc,dum,pmax,dums)
        write(*,'('' IEE a (eV) [set]               : '',f8.2)')iee_a
        write(*,'('' IEE b (eV) [set]               : '',f8.2)')iee_b

      ! automatic
      else
         call getieeab (iee_a,iee_b,ieeel,edistri,exc,nuc,ieeatm)
         call getmaxiee(iee_a,iee_b,ieeel,edistri,exc,dum,pmax,dums)
         write(*,'('' IEE a (eV) [determined]        : '',f8.2)')iee_a
         write(*,'('' IEE b (eV) [determined]        : '',f8.2)')iee_b

      endif iee1 ! EI

      if (scani  ==  0) then
         write(*,'('' maximum IEE (eV)               : '',f8.2)')exc
         write(*,'('' maximum of P(E) at (eV)        : '',f8.2)')dum
         write(*,'('' average E for P(E) (eV)        : '',f8.2)')dums
      elseif (scani  ==  1) then
         write(*,'('' minimum IEE (eV)               : '',f8.2)')lowerbound
         write(*,'('' maximum IEE (eV)               : '',f8.2)')upperbound
      endif
      write(*,*)

      ! For the different structures, use the above settings
iee2:  do i = 1, ndumpGS
        ! read structure and velo from previous GS MD
        do k = 1, nuc
          read(io_gs,*)(xyz(j,k),j=1,3),(velo(j,k),j=1,3)
        enddo

        if (icalc(i) == 0) cycle
        if (nrun == ntraj) exit

        ! 1. generate an e- that can ionize any MO and vary it by boxuller
        do 
          ! vary by normal distributed boxmuller random number
          Edum = vary_energies(eimp0,eimpw)
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

        ! this is the IEE energy used -> Edum = exc
        ! fimp is a scaling factor (default = 1.0)
        ! exc is now randomly chosen, but at least larger than ehomo
        if (scani  ==  0) then
          edum = fimp * exc / autoev
        elseif (scani  ==  1) then
          exc = lowerbound + (((upperbound-lowerbound)/dble(ntraj))*dble(nrun))
          edum = fimp * exc / autoev
        endif

        ! map IEE to MO number, shake-up possibility (ie MO2>0)
        call momap(ihomo, emo, edum+ehomo, mo1, mo2)

        ! Temperature of the molecule
        Tsoll = Edum / (0.5_wp * 3.0_wp * nuc * kB)

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

        nrun = nrun + 1

        if (.not.check) then
           call setetemp(1,Edum,etemp)
           write(*,'('' Run: '',i4,'', step on M trj: '',i5,'' MOs '',2i3,'' IEE (eV)  = '',F6.1,    &
             & '' heating time ='',F6.0,'' eTemp ='',F6.0)')nrun,i,mo1,mo2,Edum*autoev,dum,etemp
        endif
        Tav = Tav + Tsoll
        tta = tta + dum

        xyzr (1:3,1:nuc,nrun) = xyz (1:3,1:nuc)
        velor(1:3,1:nuc,nrun) = velo(1:3,1:nuc)
        velofr(   1:nuc,nrun) = velof(   1:nuc)
        eimpr (         nrun) = Edum
        taddr (         nrun) = dum*fstoau

      enddo iee2

      close(io_gs)
      !----

      dum  = 0
      dums = Tav * 0.5 * 3 * nuc * kB / dble(nrun)

      do i = 1, nrun
         dum = dum + (eimpr(i) - dums)**2
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


      !---- CID part --- 
      !no need for EIMP or TADD
      elseif ( method == 3 .or. method == 4 )then ! iee0

        do i = 1, ndumpGS

          ! read structure and velo from previous GS MD
          do k = 1, nuc
              read ( io_gs, *) ( xyz(j,k), j=1,3 ), ( velo(j,k), j=1,3 )
          enddo

          if ( icalc(i) == 0 ) cycle
          if ( nrun == ntraj ) exit

          nrun=nrun+1
          xyzr (1:3,1:nuc,nrun)=xyz (1:3,1:nuc)
          velor(1:3,1:nuc,nrun)=velo(1:3,1:nuc)
          velofr(   1:nuc,nrun)=velof(   1:nuc)
          eimpr (         nrun)=0.0_wp
          taddr (         nrun)=0.0_wp 
        enddo

        close(io_gs)

        write(*,*)
        write(*,*) nrun,' done.'
        write(*,*)

        !Take out the taddr and eimp for CID
        taddr = 0.0d0
        eimpr = 0.0d0

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
    call system('rm -rf TMPQCXMS')
    call system('mkdir  TMPQCXMS')

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
      call wrstart(i,nuc,xyzr(1,1,i),velor(1,1,i),velofr(1,i),eimpr(i),taddr(i))
    enddo

    ! sum up the information at the end of creating tmp directories (call #2)
    call info_sumup(ntraj, tstep, tmax, Tinit, trelax, eimp0, &
      & ieeatm, iee_a, iee_b, eimpact, ESI, tempESI,  & 
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
    if ( method == 2 ) then
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

    !method .eq. 3 = CID
    elseif ( method ==  3 .or. method == 4 ) then
       if ( prog == 0 ) then
          stop 'Exit: CID is not suitable with DFTB'
       elseif ( prog == 5 .and. gas%Iatom /= 6 .and. gas%Iatom /= 0 ) then !MNDO99 
          stop 'Exit: CID with MNDO99 only with N2 or He'
       !elseif (prog == 4) then
       !   stop 'Exit: CID is not suitable with MSINDO'
       endif
    endif

    ! save the coord file (using TMOL) into coord.original, because it will be overwritten
    if (iprog == 2) call system('cp coord coord.original')

    ttime = 0
    nmax0 = nmax
    icoll = 0
    isec  = 0
    asave ='NOT USED'
    tcont = 0
    chrgcont  = 1.0
    num_frags = 0
    coll_counter = 0
    frag_counter = 0

    call timing(t1,w1)

    ! Read the qcxms.start file
    call rdstart(itrj,nuc,xyz,velo,velof,tadd,eimp)

    if ( method == 3 .or. method == 4 ) then
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Initialize QM code for production runs
! Note; if method=3 (CID) then the initalizing is done within cid.f
! instead of inside the main. The reason is to avoid
! declaring the impact-atom related variables inside main.
! The cid module is called dependend on the number of collisions it
! undergoes.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Go into CID module
mCID:if ( method == 3 .or. method == 4 ) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tadd  = 0.0_wp !should be 0 anyway.
      eimp  = 0.0_wp !should be 0 anyway.
      velof = 0.0_wp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( method == 3 ) mchrg =  1
      if ( method == 4 ) mchrg = -1

      collisions = 0

      edum = 0

      ! Check if input (qcxms.in) makes sense
      call cidcheck(MaxColl, CollNo, CollSec)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Do ESI MD if ESI larger than 0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
noESI: if (.not. No_ESI )then

        write(*,*) ''
        write(*,*) '# ## ### ## # ## ### ## # ##'
        write(*,*) '# ESI simulation turned ON #'
        write(*,*) '# ## ### ## # ## ### ## # ##'
        write(*,*) ''

        call ekinet(nuc,velo,mass,edum,t)
        edum = t * (0.5 * 3 * nuc *kB * autoev)
        ENe(1) = edum

        if (ESI == 0 .and. tempESI == 0 ) then
          call random_seed()
          call random_number(a)

          b    = nuc / 10.0
          dep  = nint(b)
          if ( b < 2.0_wp ) dep = 1 ! hard coded for small molecules

          ! Make random energy depending on molecular size
          if (dep == 1) rand_int = 1 + FLOOR(2*a)     ! 1-2
          if (dep == 2) rand_int = 2 + FLOOR(2*a)     ! 2-3
          if (dep == 3) rand_int = 3 + FLOOR(3*a)     ! 3-5 !höher für anderen bereich
          if (dep == 4) rand_int = 3 + FLOOR(4*a)     ! 3-6
          if (dep >= 5) rand_int = 4 + FLOOR(5*a)     ! 4-8

          !write(*,*) 'DEP', dep

          write(*,'('' Randomly Scaling internal energy: '')')
          write(*,'(50(''-''))') 
          write(*,'('' Inner Energy before          : '', f14.6, a3)') ENe(1), ' eV'
          write(*,'('' Temperature  before          : '', f14.6, a3)')      T, ' K'
          write(*,'(50(''-''))') 
          write(*,'('' Wanted value Energy          : '', i4, a3)'  ) rand_int, ' eV'
          write(*,'(50(''-''))') 

          edum = rand_int - edum
          if (edum > 0) then
            E_Scale = vary_energies(edum, 0.1_wp)
            if (rand_int  == 0) E_Scale = 0
            write(*,'('' Scaling to inner Energy   : '', f14.6,a3)') E_Scale+ENe(1),' eV'
          else
            write(*,*) ' ! No Scaling ! '
            E_Scale = 0.0_wp
          endif

        elseif ( ESI > 0.0_wp ) then
          edum = (ESI - edum)
          if ( edum <= 0.4_wp ) then !do not scale if energy is low 
            write(*,'(''No Scaling required ...'')')
            E_Scale = 0.0_wp
          else

            write(*,'(''! Scaling internal energy !'')')
            write(*,'(50(''-''))') 
            write(*,'('' Inner Energy before          : '',f14.6,a3)') ENe(1), ' eV'
            write(*,'('' Temperature  before          : '',f14.6,a3)')    T,   ' K'
            write(*,'(50(''-''))') 
            write(*,'(''Wanted value Energy           : '',f14.6,a3)') ESI,    ' eV'
            write(*,'(50(''-''))') 

            if ( NoScale ) then ! No distribution of ESI 
              E_Scale = edum
            else
              E_Scale = vary_energies(edum, 0.2_wp)
            endif
            write(*,'('' Scaling to inner Energy      : '', f14.6,a3)') E_Scale+ENe(1),' eV'
          endif

        elseif ( tempESI > 0 ) then
          ESI = 0
          write(*,'(''! Scaling Temperature  !'')')
          write(*,'(50(''-''))') 
          write(*,'('' Inner Energy before          : '',f14.6,a3)') ENe(1), ' eV'
          write(*,'('' Temperature  before          : '',f14.6,a3)')    T,   ' K'
          write(*,'(50(''-''))') 
          write(*,'('' Wanted value Temperature     : '', f14.6,a3)') tempESI,'K'
          write(*,'(50(''-''))') 
          edum = tempESI - t
          edum =(3 * 0.5 * edum * nuc* kB)/evtoau
          E_Scale = vary_energies(edum, 0.2_wp)
          write(*,'(''Scaling to inner Energy       : '',f14.6,a3)') E_Scale+ENe(1),' eV'
        endif

        if ( tempESI > 0 .and. ESI > 0 ) then
          write(*,*) 'Cannot provide both, Energy and Temp.!'
          stop
        endif

        if (E_Scale <= 0) E_Scale = 0

        tScale = (E_Scale * 2.0d0/3.0d0) /(nuc * kB * autoev)
        ENe(4) = tScale
        tScale = tScale + T
        write(*,'('' Scaling to Temperature       : '', f14.6,a3)') tScale,' K'
        write(*,*) ' '

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
            write(*,*)'initial Cartesian coordinates:'

            do i=1,nuc
              write(*,'(i3,3f10.5,4x,a2,f10.3)')&
              &  i,xyz(1,i),xyz(2,i),xyz(3,i),toSymbol(iat(i)),mass(i)/amutoau
            enddo

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
                prestep = nint(ENe(4)/100)
                prestep = prestep*300
                starting_MD = .true.
              else
                prestep=simMD
                starting_MD = .false.
                if(fragstate == 2) prestep=simMD/2
              endif

            elseif (TempRun .and. isec == 1) then
              prestep=simMD
              starting_MD = .true.
            elseif (TempRun .and. isec > 1) then
              starting_MD = .false.
              if(fragstate == 2) prestep = prestep/2
            endif


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Pre-CID MD Loop => ESI MD
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call md(itrj,icoll,isec,nuc,prestep,xyz,iat,mass,imass,mchrg,   &
            & grad,velo,velof,list,tstep,j,nfragexit,fragm,fragf,     &
            & fragat,dumpstep,etempin,mdok,chrg,spin,axyz,tscale,tadd, &
            & eimp,.false.,Tav,Epav,Ekav,ttime,aTlast,fragstate,dtime,      &
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
            if(mdok)then

              call manage_fragments(nuc, iat, xyz, axyz, chrg, spin, mass, imass, &
                &  iprog, aTlast, itrj, icoll, isec, list, chrgcont,  &
                &  tcont, nfrag, metal3d, ECP, btf, maxsec, dtime, asave, io_res )

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
              coll_counter = coll_counter + 1
              write(*,*) 'Number of frag Processes',coll_counter
              write(*,*) 'Fragmentation in ESI MD'

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
            if ( nuc <= 6 ) nmax = nmax0 / 4 !CID may fragment more

            ! do not continue with small fragments
            if ( nuc <= 5 ) then
               small = .true.
               exit
            endif

            ! do not continue with low masses/resolution of instrument (user)
            if ( sum(mass(1:nuc)) / amutoau <=  minmass ) then
              littlemass = .true.
              exit
            endif

            ! reallocate the variables, as they change for smaller systemsizes
            deallocate(xyz,axyz,grad,velo,velof,chrg,spin,iat,list,mass,imass)

            allocate(xyz (3,nuc), &
            & axyz(3,nuc),        &
            & grad(3,nuc),        &
            & velo(3,nuc),        &
            & velof(nuc),         &
            & chrg(nuc),          &
            & spin(nuc),          &
            & iat (nuc),          &
            & list(nuc),          &
            & imass(nuc),         &
            & mass(nuc))

            do i=1,3
              xyz(i,1:nuc) = xyzn (i,1:nuc) - cema(i)
            enddo

            velo(1:3,1:nuc) = velon(1:3,1:nuc)
            iat (    1:nuc) = iatn (    1:nuc)
            mass(    1:nuc) = massn(    1:nuc)
            imass(   1:nuc) = imassn(   1:nuc)


            ! Check fragmentation state and chose what to do 
            if (tcont > 0) then
              cycle ESI_loop ! CYCLE the MD loop if fragmentation happend
            else
              exit ESI_loop  ! EXIT the MD loop if nothing happend
            endif

          enddo ESI_loop ! ENDDO loop the MD module

          if ( TempRun .or. small .or. isec > 2 ) then
            if( index(asave,'NOT USED') == 0 ) then
              write(io_res,'(a)')asave
            endif
          endif

        endif doESI ! ENDIF E_Scale gt 0

      endif noESI ! ENDIF No_ESI


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Initial caclulation of number of collisions 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cnt:  if (.not. TempRun .and. .not. small .and. isec < 3) then

        starting_md = .false.
        save_MD = simMD
        new_counter = 0

        ! Full Auto = Total Program control depending on mol. size, pressure cham. length etc.
auto:   if ( FullAuto )then

          write(*,*) 'Collisions set to automatic'

          ! Calculate radius as distance between COM and Atom that is the farest away from the COM
          call collision_setup(nuc,iat,xyz,mass,rtot,cross,mfpath,calc_collisions)

          ! distribute the no. of coll. to randomize a bit
          collisions = vary_collisions(calc_collisions)

          !write(*,'(/,80(''-''))')

        ! Coll Auto = Other run modes
        elseif ( CollAuto ) then
          collisions = set_coll ! starting coll are user set (default = 10)
          !write(*,'(/,80(''-''))')
          !write(*,*) 'Max. No of collisions : ', collisions

        elseif ( Manual ) then

          ! Users choose the fragmentation amount for all mod(itrj,x) runs
          ! The entire spectrum is pieced together from 3 different amount of fragmentations
noauto:   if ( CollSec(1) /= 0 ) then
            write(*,'(/,80(''-''))')
            write(*,*) '!!! Number of fragmentations are user set !!!'
            if (mod(itrj,20)  ==  0)then !all 20 runs are stopped after CollSec(3) fragmentations
               collisions = set_coll
               new_counter = CollSec(3)
               write(*,*) ' - Fragmentations this run: ', new_counter
            elseif (mod(itrj,3)  ==  0)then !all 3 runs are stopped after CollSec(3) fragmentations
               collisions = set_coll
               new_counter = CollSec(2)
               write(*,*) ' - Fragmentations this run: ', new_counter
            else
               collisions = set_coll
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
          elseif ( MaxColl /= 0 ) then
            write(*,'(/,80(''-''))')
            write(*,*) '!!! M+ collisions are user set !!!'
            collisions  = MaxColl
            new_counter = 1

          endif noauto

        endif auto


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start Loop for CID collision and subsequent MD simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cidlp:  do
          if(icoll /= 0)then
            write(*,*) ' '
            write(*,*) ' '
          endif

          isec = 1
          icoll=icoll+1
          fragstate = 0
          simMD = save_MD


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Write Coordinates
          write(*,'(80(''=''))')
          write(*,'(6x,a10,2x,i4,6x,a9,2x,i2,1x,a1,i2)')&
          & 'trajectory ',itrj, 'collision ',icoll,'/',collisions
          write(*,'(80(''=''),/)')
          write(*,*)'initial Cartesian coordinates:'

          do i = 1, nuc
            write(*,'(i3,3f10.5,4x,a2,f10.3)')&
            & i,xyz(1,i),xyz(2,i),xyz(3,i),toSymbol(iat(i)),mass(i)/amutoau
          enddo


          ! First collision
          write(*,'(/,'' statistical charge  = '',F10.4,/)')chrgcont

          ! set the direction for the CID module after the first coll
          if ( icoll /=  1 ) then
            direc = cm2(:) - cm1(:)
            direc = direc/norm2(direc)
          endif

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Call CID module
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call cid(nuc,iat,mass,xyz,velo,tstep,mchrg,eTempin,               &
          &  stopcid,Eimpact,axyz,ttime,eExact,ECP,           &
          &  vScale,MinPot,ConstVelo,cross,                   &
          &  mfpath,rtot,chrg,icoll,collisions,direc,new_velo,aTlast,         &
          &  calc_collisions,imass)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          new_temp = aTlast

          ! Grad failure (too big steps)
          if ( stopcid ) write(*,*) 'Error occured in the CID module'

          ! END STUFF HERE
          if ( stopcid ) exit

          call manage_fragments(nuc, iat, xyz, axyz, chrg, spin, mass, imass, &
            &  iprog, aTlast, itrj, icoll, isec, list, chrgcont,  &
            &  tcont, nfrag, metal3d, ECP, btf, maxsec, dtime, asave, io_res)


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
            coll_counter = coll_counter + 1
            write(*,*) 'Number of frag Processes',coll_counter

          elseif (tcont == 0) then
            do i = 1, nuc
              xyzn (1:3,i)=xyz (1:3,i)
              velon(1:3,i)=velo(1:3,i)
              iatn (    i)=iat (    i)  
              massn(    i)=mass(    i)  
              imassn(   i)=imass(   i)  
              dum         =dum+iatn(i)
              cema(1:3)   =cema(1:3)+xyzn(1:3,i)*iatn(i)
            enddo
          endif

          ! move to Center-of-mass
          cema(1:3) = cema(1:3) / dum

          ! do not continue with small fragments
          if ( nuc <= 5 ) then
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
          deallocate(xyz,axyz,grad,velo,velof,chrg,spin,iat,list,mass,imass)

          allocate(xyz (3,nuc), &
          &        axyz(3,nuc), &
          &        grad(3,nuc), &
          &        velo(3,nuc), &
          &        velof(nuc),  &
          &        chrg(nuc),   &
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


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!        Mean-free-path MD                                       !!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! Loop MD simulations if fragmentation occurs until nothing happens
MFPloop:  do
            isec = isec + 1

            if ( fragstate == 2 ) simMD = simMD / 2

            ! calculate the center-of-mass and reset for correct collision sim
            call cofmass(nuc,mass,xyz,cm)
            cm1(:) = cm(:) 

            write(*,'(/,80(''=''))')
            write(*,*)'  Collision ',icoll,' MD trajectory            ',isec
            write(*,'(80(''=''),/)')
            write(*,*)'initial Cartesian coordinates:'

            do i = 1, nuc
              write(*,'(i3,3f10.5,4x,a2,f10.3)')i,xyz(1,i),xyz(2,i),xyz(3,i), &
              &    toSymbol(iat(i)) , mass(i) * autoamu
            enddo

            write(*,'(/,'' statistical charge  = '',F10.4,/)')chrgcont

            ! HS-UHF ini for closed-shells allowed
            mspin=0

            call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)

            ! a second attempt if this fails
            if ( .not. iniok ) then
              iniok=.false.
              call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)
              if ( .not. iniok ) stop 'fatal QC error. Must stop!'
            endif

            ! do production MD
            Tdum=0

            ! scale the sim MD if fragmentation occurs
            ! if (isec == 2) simMD =int(simMD/2)
            ! if (isec == 3) simMD =int(simMD/3)
            ! if (isec == 4) simMD =int(simMD/4)
            ! if (isec == 3) simMD =int(simMD/2)
            ! if (isec == 4) simMD =int(simMD/3)


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Do Mean-free-path (MFP) MD with simMD timesteps
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call md(itrj,icoll,isec,nuc,simMD,xyz,iat,mass,imass,mchrg,grad,&
            &       velo,velof,list,tstep,j,nfragexit,                      &
            &       fragm,fragf,fragat,dumpstep,etempin,                    &
            &       mdok,chrg,spin,axyz,                                    &
            &       Tdum,tadd,eimp,.false.,Tav,Epav,Ekav,ttime,aTlast,      &
            &       fragstate,dtime,ECP,.false.,new_velo)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            write(*,'(/10x,''Results'',/,'' average Ekin '',F12.6)')Ekav
            write(*,'('' average Epot  '',F12.6)')Epav
            write(*,'('' average Etot  '',F12.6)')Epav+Ekav
            write(*,'('' average T (acc.)'',F12.1)')Tav
            write(*,'('' Temp without acc '',F12.1)') new_Temp
            write(*,'('' average last T'',F12.1)')aTlast

            ! if (isec == 2) simMD = simMD*2
            ! if (isec == 3) simMD = simMD*3
            ! if (isec == 4) simMD = simMD*4
            ! if (isec == 3) simMD =simMD*2
            ! if (isec == 4) simMD =simMD*3

            ! calculate the new center-of-mass as reference
            call cofmass(nuc,mass,xyz,cm)
            cm2(:) = cm(:)
            !!!

            ! cm_out = sqrt(cm2(1)**2+cm2(2)**2+cm2(3)**2)
            !.       - sqrt(cm1(1)**2+cm1(2)**2+cm1(3)**2)
            !
            ! new_velo = (cm_out / (ttime*fstoau)) / mstoau
            ! new_velo = (cm_out / (2000*fstoau)) / mstoau
            ! ttime=0 ! Set the md time to 0 for correct simulation time (maybe change)


            ! printout and count
            if(mdok)then
              call manage_fragments(nuc, iat, xyz, axyz, chrg, spin, mass, imass, &
                &  iprog, aTlast, itrj, icoll, isec, list, chrgcont,  &
                &  tcont, nfrag, metal3d, ECP, btf, maxsec, dtime, asave, io_res )

               !If MD fails
              write(*,*)
            else
              write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              write(*,*)'something went wrong! Dont worry, the run'
              write(*,*)'is just not further counted.'
              write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              stop
            endif

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
              coll_counter = coll_counter + 1
              write(*,*) 'Number of frag Processes',coll_counter

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

            ! reduce the simulation time
            if ( nuc <= 6 ) nmax = nmax0 / 4 !CID may fragment more

            ! do not continue with small fragments
            if ( nuc <= 5 ) then
              small = .true.
              exit
            endif

            ! do not continue with low masses/resolution of instrument (user)
            if ( sum(mass(1:nuc)) / amutoau <=  minmass ) then
              littlemass = .true.
              exit
            endif

            ! reallocate the variables, as they change for smaller systemsizes
            deallocate(xyz,axyz,grad,velo,velof,chrg,spin,iat,list,mass,imass)

            allocate(xyz (3,nuc), &
            & axyz(3,nuc),        &
            & grad(3,nuc),        &
            & velo(3,nuc),        &
            & velof(nuc),         &
            & chrg(nuc),          &
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
          E_COM = (beta * E_KIN) * autoev
          !write(*,*) 'BETA      :    ',beta
          !write(*,*) 'ECOM (eV) :    ', E_COM
          !write(*,*)'NEW VELO MD:', new_velo


          if ( new_velo <= 800 .or. E_COM <= 0.4 .and. MinPot == 0 ) then
            if(index(asave,'NOT USED') == 0)then
              write(io_res,'(a)')asave
            endif
            write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*) ' ! Low velocity and thus low E(COM) !'
            write(*,*) '   -> The velocity',new_velo, '<-    !'
            write(*,*) '   -> E_COM       ',E_COM, '<-    !'
            write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            exit
          endif

          !--------------------------------------------------------------
          ! 3.1.) If the number of fragmentaion steps are exceeded
          !Collauto  == false
          if ( Manual .and. new_counter > 0 ) then
            if ( coll_counter >= new_counter ) then
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
Full:     if (FullAuto .and. coll_counter > frag_counter) then

            frag_counter = coll_counter

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
Coll:     if (CollAuto .and. coll_counter > frag_counter) then

            ! Set-up collision number and vary the amount (cid.f90)
            call collision_setup(nuc,iat,xyz,mass,rtot,cross,mfpath, &
            &    calc_collisions)

            if ( collisions > 0 )then
               call random_seed(numb)
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
          endif Coll

          ! 4.) If the max number of collisions is reached
          if ( icoll >= collisions ) then
            write(*,*)'-------------------------------------------------'
            write(*,*)' Last MD is finished. All is written out'
            write(*,*)'-------------------------------------------------'
            write(io_res,'(a)')asave

            !if(num_frags == 0)then
            !   write(*,*)''
            !   write(*,*)'------------------------------------------'
            !   write(*,*)'    No fragmentation in the simulation!   '
            !   write(*,*)'   Increase energy or time of sampling.   '
            !   write(*,*)'------------------------------------------'
            !endif

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
           write(*,*)'     run not continued because of small fragment(s)'
           write(*,*)''
        endif

        ! Write out too light molecules
        if(littlemass)then
           if(index(asave,'NOT USED') == 0)then
              write(io_res,'(a)')asave
           endif
           write(*,*)'     run not continued because of resolution'
           write(*,*) 'Input      : ', minmass
           write(*,*) 'Actual mass: ', sum(mass(1:nuc))/amutoau
           write(*,*)''
        endif

        ! ERROR
        if(stopcid)then
           if(index(asave,'NOT USED') == 0)write(io_res,'(a)')asave
           write(*,*)'     run aborted, last structure saved '
           write(*,*)''
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

ei: if (method /= 3 .and. method /= 4)then
      if (method  ==  0 .or. method  ==  1 ) mchrg =  1 ! EI
      if (method  ==  2 )                    mchrg = -1 ! DEA

loop: do 
        isec=isec+1

        write(*,'(/,80(''=''))')
        write(*,*)'                      trajectory ',itrj,isec
        write(*,'(80(''=''),/)')
        write(*,*)'initial Cartesian coordinates:'

        do i = 1, nuc
           write(*,'(i3,3f10.5,4x,a2,f10.3)')&
           & i,xyz(1,i),xyz(2,i),xyz(3,i),toSymbol(iat(i)),mass(i)/amutoau
        enddo

        write(*,'(/,'' statistical charge  = '',F10.4,/)')chrgcont

        ! re- initialize the GBSA module
        !         if (lgbsa)then
        !           call ekinet(nuc,velo,mass,edum,t)
        !           call init_gbsa(nuc,iat,solvent,gsolvstate,t)
        !         endif
        
        ! HS- UHF ini for closed-shells allowed
        mspin = 0
        call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)

        ! a second attempt if this fails
        if ( .not. iniok ) then
           iniok =.false.
           call iniqm(nuc,xyz,iat,mchrg,mspin,betemp,edum,iniok,ECP)
           if ( .not. iniok) stop 'fatal QC error. Must stop!'
        endif

        if ( isec > 1 ) nfragexit = 2

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  do production MD
        Tdum=0
        call md(itrj,0,isec,nuc,nmax,xyz,iat,mass,imass,mchrg,grad, &
        &       velo,velof,list,tstep,j,nfragexit,                  &
        &       fragm,fragf,fragat,dumpstep,etempin,                &
        &       mdok,chrg,spin,axyz,                                &
        &       Tdum,tadd,eimp,.false.,Tav,Epav,Ekav,ttime,aTlast,  &
        &       fragstate,dtime,ECP,.false.,0.0_wp)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        write(*,'(/10x,''Results'',/,'' average Ekin '',F12.6)')Ekav
        write(*,'('' average Epot  '',F12.6)')Epav
        write(*,'('' average Etot  '',F12.6)')Epav+Ekav
        write(*,'('' average T     '',F12.1)')Tav
        write(*,'('' average last T'',F12.1)')aTlast

        ! printout and count
okmd:   if(mdok)then
          call manage_fragments(nuc, iat, xyz, axyz, chrg, spin, mass, imass, &
          &  iprog, aTlast, itrj, icoll, isec, list, chrgcont,  &
          &  tcont, nfrag, metal3d, ECP, btf, maxsec, dtime, asave, io_res )


          ! this is likely to be wrong (too many fragmnets for EI)
          if ( nfrag > 5 ) then
            write(*,*)'    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*)'    !! something went wrong! Dont worry, the run !!'
            write(*,*)'    !! is just not counted.                      !!'
            write(*,*)'    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            exit
          endif

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

            deallocate(xyz,axyz,grad,velo,velof,chrg,spin,iat,list,mass,imass)

            allocate(xyz (3,nuc), &
            &  axyz(3,nuc),       &
            &  grad(3,nuc),       &
            &  velo(3,nuc),       &
            &  velof(nuc),        &
            &  chrg(nuc),         &
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
  call system('touch ready')
  close(io_res)
  close(io_log)


  call timing(t2,w2)
  write(*,'(/,'' wall time (min)'',F10.2  )')(w2-w1)/60.0_wp
  write(*,'(  '' # of QC calls  '',I10  ,/)')calls

  call system('date')
  call version(2)
  write(*,*)'normal termination of QCxMS'

end program QCxMS
