module qcxms_input
  use readcommon
  use cidcommon
  use common1 
  use newcommon
  use qcxms_info, only: qcstring
  use mctc_env, only : error_type, get_argument!, fatal_error
  use mctc_io, only : structure_type, read_structure, write_structure, &
      & filetype, get_filetype, to_symbol
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_convert

  use get_settings, only: collision_type
  implicit none

  contains

subroutine input(tstep,tmax,ntraj,etemp_in,Tinit, mchrg_prod,                  &
        iee_a,iee_b,eimp0,eimpw,fimp,iprog,trelax,hacc,nfragexit,maxsec,       &
        edistri,btf,ieeatm,                                                    &
        scanI,lowerbound,upperbound,ELAB,ECOM, eExact,ECP,unity,noecp,         &
        nometal,vScale,CollNo,CollSec,ConstVelo,                               & 
        minmass,etempGS, coll,                       & 
        MinPot,ESI,tempESI,No_ESI,NoScale,manual_dist, legacy)
!  use gbobc, only: lgbsa
    
  integer  :: iprog
  integer  :: ntraj
  integer  :: nfragexit
  integer  :: maxsec      
  integer  :: edistri     
  integer  :: nn,i,j
  integer  :: io_in
  integer  :: error
  integer  :: mchrg_prod
  
  real(wp) :: tstep
  real(wp) :: tmax 
  real(wp) :: etemp_in
  real(wp) :: Tinit
  real(wp) :: eimp0
  real(wp) :: eimpw
  real(wp) :: fimp 
  real(wp) :: iee_a,iee_b       
  real(wp) :: ieeatm            
  real(wp) :: trelax
  real(wp) :: hacc  
  real(wp) :: btf              
  real(wp) :: xx(10),axi
  real(wp) :: lowerbound
  real(wp) :: upperbound
  real(wp) :: eTempGS
  
  
!  interface
!    function to_upper(strIn) result(strOut)
!      implicit none
!      character(len=*), intent(in) :: strIn
!      character(len=len(strIn)) :: strOut
!    end function
!  end interface
  
  logical :: ex
  logical :: ECP
  logical :: unity
  logical :: scanI 
  logical :: noecp,nometal
  logical :: Plasma
  logical :: legacy
  ! logical gbsa
  
  !-----------------------------------
  !----- CID related parameters ------
  
  !integer  :: MaxColl
  integer  :: CollNo(3)
  integer  :: CollSec(3)
  integer  :: minmass
  !integer  :: set_coll
  integer  :: manual_dist
  
  real(wp) :: ELAB,ECOM
  real(wp) :: vScale,ESI,tempESI
  real(wp) :: MinPot
  
  logical  :: ConstVelo
  logical  :: No_ESI
  logical  :: NoScale
  logical  :: eExact

  type(collision_type) :: coll

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! DEFAULTS: !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!
  ! set defaults to 0

  !!!
  ! method = 0 [EI], = 1 [CSC], = 2[DEA] , = 3 [CID]
  method = 0 
  ! SQM hamiltonians
  ihamilt = 0
  ! basis: SV(P)
  bas = 3      
  ! func (default: PBE0)
  func = 0      

  ! prog: 0=dftb, 1=mopac, 2=TM, 3=orca, 4=msindo, 5=mndo99, 6=callxtb, 7=xtb, 8=xtb2
  ! default: GFN2-XTB
  prog = 8 ! GFN2-xTB
  gfnver = 3 

  ! which prog for IP:  GFN2-xTB
  iprog = 8

  ! Calculate MO spectrum with XTB (default = True)
  ! IF =False then use ORCA
  XTBMO=.True.

  ! make charge dependend on input (in Prods min. = 1)
  mchrg_prod = 1
  
  ! maximum simulation time
  tmax = 5
  ! MD time step in fs (0.5-1 fs) 
  tstep = 0.5_wp
  ! # of frag traject 
  ntraj = -1  
  ! rand ini
  !iseed = 42
  iseed = 0
  ! electronic temperature (Fermi smearing) and used for IP/Boltzmann
  ! -1 means automatic      
  etemp_in = -1        
  ! neutral MD temperature 
  Tinit = 500
  !> switch on etemp 
  No_eTemp = .false.
  ! GS Etemp (to converge radicals etc)
  etempGS=298.15 ! normal 

  !!!             !!!! 
  !!! QC Settings !!!!
  !!!             !!!! 
  ! gcp for DFT
  gcp = .False. 
  ! conventional integral vs direct [in orca only]
  noconv = .false.
  ! slowconv for convergence issues [in orca only]
  slowconv = .false.
  ! gridsize
  grid_orca = 2      
  !grid_tmol = 4   
  ! DFTB3
  hhmod = .true. 
  tbcorr = .false.
  ! tunneling correction
  tun = .false.
  ! dftb corr pot
  cab = 0
  alpc = 1.5
  ! D3 parameter
  a1 = 0
  s8 = 0
  a2 = 0
  axi=-1
  ! default parameters for scan
  lowerbound = 0.0d0
  upperbound = 20.0d0
  ! ECP
  ECP = .false.
  
  ! Do not allow checking for ECP or Metals (DEFAULT=FALSE)
  noecp = .False.
  nometal = .False.
  ! Print out neutrals as well (default: no)
  Plasma = .False.
  ! GBSA
!  lgbsa=.false.
  ! memory for QC prog in Mb, input in Gb
  qcmem = 5000      
  ! ORCA version
  orca_version = 0
  ! ORCA nprocs
  nproc_orca = 0
  
  !!!           !!!! 
  !!!!!!! EI !!!!!!!
  !!!           !!!! 
  ! 70 eV standard impact
  eimp0 = 70.0_wp
  ! width (rel)
  eimpw = 0.10_wp 
  ! average IEE per atom (eV)
  ieeatm = 0.6_wp
  ! relate IEE to eTemp in SCF
  ieetemp = 0.0_wp
  ! automatically det.
  iee_a = -99
  iee_b = -99
  ! scaling factor for Boltzmann/IP temp.
  btf = 1.0_wp
  ! scaling of comuted impact energy for fragmentation 
  fimp = 1.0_wp
  ! special H atom acceleration factor
  hacc = 3.0_wp
  ! poisson IEE distribution =1 (=0 Gauss)
  edistri = 1
  ! average heating time per eV Egap
  trelax = 2000.0_wp
  ! how many secondary runs + 1 max?
  maxsec = 7
  ! exit when nfragexit fragments are produced
  nfragexit = 3     
  ! unity scaling
  unity = .False.
  ! Legacy IEE pseudo-randomized (default: NO)
  legacy = .false.
  
  !!!           !!!! 
  !!!!!! CID !!!!!!!
  !!!           !!!!
  ELAB       =  40.0_wp  ! The laboratory energy frame 
  ECOM       =  0.0_wp  ! The center-of-mass energy frame 
  gas%Iatom  = 0        ! Index of collision atom
  manual_dist  = 0

  !!!                 !!!
  ! Different run types !

  !--- The automated run-types (see paper) --- !

  ! 1) Forced activation run-type
  CollAuto   = .False.    ! Coll. until Fragmentation
  coll%set_coll   = 10         ! set_coll max. collisions 

  ! 2) General activation run-type using coll. cell parameters
  FullAuto        = .False.     
  cell%TGas       = 300.0 ! (K) Temperatur of Gas 
  cell%PGas       = 0.132 ! (Pa = 1mTorr) Pressure of Gas
  cell%lchamb     = 0.200 ! (m = 20,0 cm) coll. cell. length
    
  ! 3) Thermal activation run-type
  TempRun    = .False.  
  tempESI    = 0.0_wp     ! Temperature for scaling
  ESI        = 0.0_wp     ! Internal E for scaling (also pre-MD)

  !----------------------------------------------------------------------!
  ! More control / different run-types 
  ! Manual run-types
  Manual = .false.
  coll%max_coll    = 0        ! Max coll.          -> only M+ + Gas (no fgc)

  ! Max overall coll.  -> all coll. (fgc) 
  CollNo(1)  = 0        ! 3 Values can be set, that determine how many 
  CollNo(2)  = 0        ! runs undergo how many coll. (see below)
  CollNo(3)  = 0        

  ! Number of fragmentations until stop
  CollSec(1) = 0        ! 3 Values can be set, that determine how many  
  CollSec(2) = 0        ! runs stop after CollSec fragm. (see below)
  CollSec(3) = 0        

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Extra settings
  eExact       = .False.  ! switch off ELAB velocity scaling (exact velocities)
  minmass      = 45       ! set resolution (lower masses are cut)
  No_ESI       = .false.  ! Don't do pre-scaling of int. Energy (if true)
  NoScale      = .false.  ! no distributing ESI energy (if true)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!                 !!!
  ! Experimental stuff  !

  ConstVelo  = .False.  ! scale to const. velo after collision 
  MinPot     = 0        ! velocity scaling as well
  vScale     = 0.00     ! velocity scaling as well
  
  
  !--------------------------------------------------------------
  !-   Read the qcxms.in file or create it when not existent   -
  !--------------------------------------------------------------
  
  inquire(file='qcxms.in',exist=ex)
  ! create dummy for copy      
  if(.not.ex)then
     open(file='qcxms.in', newunit = io_in, status='new')
     write(io_in,'(''# EMPTY INPUT FILE'')')
     close(io_in)
  else
     write(*,*) 'changed by input:'
     open ( file='qcxms.in', newunit = io_in, status = 'old' , iostat = error)
     if ( error > 0 ) error stop '!!! error while opening qcxms.in !!!'
  
     ! Read line-by-line
     do
       read(io_in, '(a)', iostat = iocheck)line    

       ! Check for errors in the input 
       if (iocheck>0)then     !Fail
         write(*,*) 'Something is wrong in the input. Exiting...'
         stop

       ! End-of-file
       elseif (iocheck<0)then !EOF
         exit

       ! Read the file
       else
        write(*,'(''>'',a)') line
  
        ! USE TO_UPPER TO MAKE THE STRING CAPITALIZED
        line = to_upper(line)
        ! TRIM THE LINE
        line = trim(line)
  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! METHODS 
        if(line == 'EI') then
          method=0
          ! length of production runs in ps
        endif
        if(line == 'CSC') then
          method=1
        endif
        if(line == 'DEA') then
          mchrg_prod = -1
        endif

        if(line == 'CID')then
          method=3
          gas%Iatom = 3 !Argon
        endif
  
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! PROGRAMS 

         ! DFTB
         if(line == 'DFTB')         prog=0

         ! MOPAC
         if(line == 'MOPAC') then
                                    prog=1
                                    ihamilt = 4
         endif

         ! TURBOMOLE
         if(line == 'TMOL' .or. &
            line == 'TURBOMOLE')then
                                    prog  = 2
                                    iprog = 2
         endif

         ! ORCA !!!!!
         if(line == 'ORCA' .or. line == 'ORCA5')then
                                    prog  = 3
                                    iprog = 3
                                    orca_version = 5
         endif
         if(line == 'ORCA4')then
                                    prog  = 3
                                    iprog = 3
                                    orca_version = 4
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! MSINDO
         if(line == 'MSINDO')       prog=4

         ! MNDO
         if(line == 'MNDO' .or. &
            line == 'MNDO99') then
                                    prog=5
                                    iprog=5
                                    ihamilt = 8
         endif

         ! Call upon XTB (OLD, not used anymore)
         if(line == 'CALL_XTB')     prog=6

         ! GFN1-xTB
         if(line == 'XTB') then
                                    prog=7 
                                    gfnver =1
                                    iprog =7  !calculate IP with XTB
         endif

         ! GFN2-xTB with D3-Disp
         if(line == 'XTB2 D3') then
                                    prog=8
                                    gfnver=2  !Use D3-Method
                                    iprog=7   !calculate IP with IPEA
         endif

         ! GFN2-xTB with D4-Disp
         if(line == 'XTB2') then
                                    prog=8
                                    gfnver=3  !Use D4-Method
                                    iprog=8   !calculate IP with XTB2 
         endif
  
  
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! MO SPECTRUM - only two options [XTB/ORCA]
         if ( line == 'MO-XTB')       XTBMO=.True.
         if ( line == 'MO-ORCA')      XTBMO=.False.

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! IP CALCULATIONS
         if ( line == 'IP-MOPAC')     iprog=1
         if ( line == 'IP-TMOL')      iprog=2
         if ( line == 'IP-MNDO')      iprog=5  !OM2 standard
         if ( line == 'IP-XTB')       iprog=7  !Use XTB  
         if ( line == 'IP-XTB2')      iprog=8  !Use XTB2 (default)
         if ( (line == 'IP-ORCA' .or. line == 'IP-ORCA5') &
           .and. orca_version == 0) then
                                      iprog = 3
                                      orca_version = 5
         elseif ( line == 'IP-ORCA4' .and. orca_version == 0 ) then
                                      iprog = 3
                                      orca_version = 4
         endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! SEMI-EMPIRCAL HAMILTONIANS
         if ( line == 'AM1')          ihamilt=1
         if ( line == 'PM3')          ihamilt=2
         if ( line == 'RM1')          ihamilt=3
         if ( line == 'PM6')          ihamilt=4
         if ( line == 'OM2')          ihamilt=6
         if ( line == 'OM3')          ihamilt=8
         if ( line == 'MNDOD')        ihamilt=10

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BASIS SET FOR ORCA/TMOL 
         ! Ahlrichs        !
         if ( line == 'SV')           bas=1
         if ( line == 'SVX')          bas=2
         if ( line == 'SV(P)')        bas=3
         if ( line == 'SVP')          bas=4
         if ( line == 'TZVP')         bas=5
         if ( line == 'QZVP')         bas=13
         !def2-Ahlrichs    
         if ( line == 'DEF2-SV(P)')   bas=9 !not available
         if ( line == 'DEF2-SVP')     bas=10
         if ( line == 'DEF2-TZVP')    bas=11
         if ( line == 'DEF2-QZVP')    bas=12
         !Minimally augmented (MA) Ahlrichs [THIS IS NOT AVAILABLE FOR TMOL]
         if ( line == 'MA-DEF2-SVP')  bas=6
         if ( line == 'MA-DEF2-TZVP') bas=7
         if ( line == 'MA-DEF2-TZVPP')bas=8
         if ( line == 'MA-DEF2-QZVP') bas=14
         !XC-Functionals
         if ( line == 'PBE0')         func=0
         if ( line == 'PBE38')        func=1
         if ( line == 'PBE12')        func=2
         if ( line == 'LDA12')        func=3
         if ( line == 'M062X')        func=4
         if ( line == 'PBE')          func=5
         if ( line == 'B97D')         func=6
         if ( line == 'B3LYP')        func=7 !there was a comment saying this is only for orca, is this true?
         if ( line == 'PW6B95')       func=8
         if ( line == 'B3PW91')       func=9
         if ( line == 'BLYP')         func=10
         if ( line == 'BP86')         func=11
         if ( line == 'TPSS')         func=12
         if ( line == 'REVPBE')       func=13
         if ( line == 'PBEH3C')       func=14
         if ( line == 'BHLYP')        func=15  
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! DFTB3
         if( line == 'NO-HHMOD')     hhmod=.False.
         if( line == 'TBCORR')       tbcorr=.True.
  
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! MISC Logicals
         if ( line == 'TUN')          tun=.true.      !TUNNELING EFFECTS
         if ( line == 'GCP')          gcp=.true.      !USE GCP
         if ( line == 'NO-GCP')       gcp=.False.     !DO NOT USE GCP
         if ( line == 'ECP')          ECP=.true.      !USE ECP (FOR ORCA/TMOL ONLY)
         if ( line == 'NOCONV')       noconv=.true.   !USE NOCONV (FOR ORCA)
         if ( line == 'SLOWCONV')     slowconv=.true. !USE SLOWCONV ( FOR ORCA)
         if ( line == 'UNITY')        unity=.true.    !USE UNITY VELOCITY SCALING
         if ( line == 'NO-ECP')       noecp=.true.    !Do not check for ECP
         if ( line == 'NO-METAL')     nometal=.true.  !Do not check for metal
         if ( line == 'PLASMA')       Plasma = .True.  ! switch off ESI energy distribution
         if ( line == 'LEGACY')       Legacy = .True. ! Legacy support (for IEE dist.)
         if ( line == 'NOETEMP')      No_eTemp = .True. ! Shut off etemp (test for CID)
  
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! CID Logicals 
          if ( line == 'COLLAUTO' ) then              ! 'forced activation' mode
            CollAuto = .True.                          
          endif
          if ( line == 'FULLAUTO' ) then              ! 'general activation' mode
            FullAuto = .True.
          endif
          if ( line == 'TEMPRUN' ) then               !'Thermal activation' mode
             TempRun=.true.
             gas%Iatom = 0
          endif
          if ( line == 'EEXACT')      eExact = .True.   ! switch off CE distribution
          if ( line == 'CONSTVELO')   ConstVelo = .True.! constant velo after collision 
          if ( line == 'NOESI')       No_ESI = .True.   ! switch off auto ESI 
          if ( line == 'NOSCALE')     NoScale = .True.  ! switch off ESI energy distribution
  
  
         ! CID Atoms !change to right iat etc.
  
         if(line == 'IATOM HE') then
           gas%Iatom   = 1  ! Helium
           gas%IndAtom = 2               ! 10 is Neon 
           gas%mIatom  = 4.002 * amutoau   !mass of atom
           gas%rIatom  = 2.64560263       ! bohr vdw-radius of the atom
         endif
  
         !ATOM = Neon
         if(line == 'IATOM NE') then
           gas%Iatom   = 2  ! Helium
           gas%IndAtom = 10               ! 10 is Neon 
           gas%mIatom  = 20.18 * amutoau   !mass of atom
           gas%rIatom  = 2.91016289       ! bohr vdw-radius of the atom
         endif
  
         !ATOM = Argon <== DEFAULT 
         if(line == 'IATOM '.or. line == 'IATOM AR'.or. gas%Iatom == 3)then
           gas%Iatom   = 3  ! default = Argon
           gas%IndAtom = 18               ! 18 is Argon (index of atom)
           gas%mIatom  = 39.948 * amutoau  ! mass of index atom
           gas%rIatom  = 3.55266638       ! bohr vdw-radius of the atom
         endif
  
         !ATOM = Krypton
         if(line == 'IATOM KR')then
          gas%Iatom   = 4  ! Krypton
          gas%IndAtom = 36               ! 36 is Krypton 
          gas%mIatom  = 83.798 * amutoau  ! mass of atom
          gas%rIatom  = 3.81722665       ! bohr vdw-radius of the atom
         endif
  
         !ATOM = Xenon
         if(line == 'IATOM XE')then
           gas%Iatom   = 5  ! Xenon
           gas%IndAtom = 54               ! 54 is Xenon
           gas%mIatom  = 131.29 * amutoau  ! mass of atom
           gas%rIatom  = 4.081786915      ! bohr vdw-radius of the atom
         endif
  
         !ATOMS = Nitrogen(x2)
         if(line == 'IATOM N2')then
           gas%Iatom   = 6  ! Nitrogen
           gas%IndAtom = 7                ! 7 is Nitrogen. 
           gas%mIatom  = 14.007 * amutoau  ! mass of atom
           !  gas%rIatom = 1.341698476 * 2   ! bohr cov-radius of the molecule -> must be fixed for real radius!
           gas%rIatom  = 3.64              !  bohr kinetic-radius of the molecule !fast fix
         endif
         
         !DUMMY
         if(line == 'IATOM DUM')then
            gas%Iatom   = 0  
            gas%IndAtom = 0                 
            gas%mIatom  = 0 
            gas%rIatom  = 0 
         endif
  
  
  !       GBSA (solvation model...most prob unimportant for CID)
  !        if(index(line,'GBSA') /= 0)then 
  !           lgbsa = .true.
  !           solvent='h2o' ! standard H2O, rest not implemented 
  !           if(index(atmp,'-') == 0.and.atmp(1:1) /= ' ')then
  !                solvent=line(i+1)
  !           endif
  !        endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! START READING ACTUAL SIMULATION VALUES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  !       TIMESTEP
          if(index(line,'TSTEP') /= 0)then            
             call readl(line,xx,nn)
             tstep=xx(1)
          endif
  !       MAX SIM. TIME
          if(index(line,'TMAX' ) /= 0)then            
             call readl(line,xx,nn)
             tmax=int(xx(1))
          endif
  !       INITIAL TEMP
          if(index(line,'TINIT') /= 0)then            
             call readl(line,xx,nn)
             Tinit=xx(1)
          endif
  !       NUMBER OF TRAJ.
          if(index(line,'NTRAJ') /= 0)then            
             call readl(line,xx,nn)
             ntraj=int(xx(1))
          endif

  !       CHARGE IN PROD. RUNS
          if(index(line,'CHARGE') /= 0)then            
             call readl(line,xx,nn)
             mchrg_prod=int(xx(1))
          endif

  !       IMPACT ENERGY SREAD (also scaling for CID-ESI)
          if(index(line,'EIMPW') /= 0)then            
             call readl(line,xx,nn)
             eimpw = xx(1)
          endif

  !       NUMBER OF FRAGMENTS FOR EXIT
          if(index(line,'NFRAGEXIT') /= 0)then            
             call readl(line,xx,nn)
             nfragexit=int(xx(1))
          endif
  !       MAX. NUMBER OF CASCADING RUNS
          if(index(line,'MAXSEC') /= 0)then            
             call readl(line,xx,nn)
             maxsec=int(xx(1))
          endif
  !       SEED FOR RANDOM NO. GENERATOR
          if(index(line,'ISEED') /= 0)then            
             call readl(line,xx,nn)
             iseed=int(xx(1))
          endif
  !       HYDROGEN (additional) VELOCITY ACC.
          if(index(line,'HACC' ) /= 0)then            
             call readl(line,xx,nn)
             hacc=xx(1)
          endif
  !       AMOUNT OF EXACT EXCHANGE 
          if(index(line,'FOCK' ) /= 0)then            
             call readl(line,xx,nn)
             axi=xx(1)
          endif
  !       MEMORY
          if(index(line,'MEM' ) /= 0)then            
             call readl(line,xx,nn)
             qcmem=int(1000*xx(1))
          endif
  !       ELECTRONIC TEMP. (Fermi Smearing) in prod runs
          if(index(line,'ETEMP') /= 0)then            
             call readl(line,xx,nn)
             if(nn.gt.0)  etemp_in=xx(1)
          endif
          ! set etemp in equilibration MD
          if(index(line,'ETEMPGS') /= 0)then            
             call readl(line,xx,nn)
             etempGS=xx(1)
          endif

          ! set grid for TMOL (not yet implemnted)
          !if ( prog == 2 .and. index(line,'GRID') /= 0 )then            
          !     call readl(line,xx,nn)
          !     grid_tmol=xx(1)
          !endif

          ! set grid for ORCA
          if ( prog == 3 .and. index(line,'GRID') /= 0)then            
               call readl(line,xx,nn)
               grid_orca=xx(1)
          endif

          !! set numthreads orca
          !if ( prog == 3 .and. index(line,'NPROC') /= 0)then            
          !     call readl(line,xx,nn)
          !     nproc_orca = xx(1)
          !     write(*,*) 'ORCA PROCS', nproc_orca
          !endif

  
  !       Boltzmann factor/weight
          if(index(line,'BTF') /= 0)then            
             call readl(line,xx,nn)
             btf=xx(1)
          endif
  
  ! ---------------------------------------------------------------
  !!! EI PARAMETERS !!! 
  ! ---------------------------------------------------------------
  
  !       ELECTRON IMPACT ENERGY 
          if(index(line,'EIMP0') /= 0)then            
             call readl(line,xx,nn)
             eimp0=xx(2)
          endif
  !       SCALING FACTOR FOR IMPACT ENERGY
          if(index(line,'FIMP' ) /= 0)then 
             call readl(line,xx,nn)
             fimp=xx(1)
          endif
  !       RELAXATION TIME
          if(index(line,'TRLEAX' ) /= 0)then            
             call readl(line,xx,nn)
             trelax=xx(1)
          endif
  
  ! IEE settings !
          if(index(line,'IEEATM') /= 0.or.index(line,'IEEBND') /= 0)then !IEE PER ATOM
             call readl(line,xx,nn)
             ieeatm=xx(1)
          endif
          !     POISSON OR GAUSSIAN DISTRI.  !
          if(index(line,'POISSON') /= 0)   then
             edistri=1                                       
             iee_a=0.2
             iee_b=1.
             call readl(line,xx,nn)
             if(nn == 2)then
                iee_a=xx(1)
                iee_b=xx(2)
             endif
          endif
          if(index(line,'GAUSS') /= 0)   then
             edistri=0   
             iee_a=0.2
             iee_b=0.5
             call readl(line,xx,nn)
             if(nn == 2)then
                iee_a=xx(1)
                iee_b=xx(2)
             endif
          endif
          !  IEETEMP
          if(index(line,'IEETEMP' ) /= 0)then            
             call readl(line,xx,nn)
             ieetemp=xx(1)
          endif
  !       SCAN FUNCTION 
          if(index(line,'UPPER') /= 0)then            
             call readl(line,xx,nn)
             upperbound=xx(1)
             scanI = .true.
          endif
          if(index(line,'LOWER') /= 0)then            
             call readl(line,xx,nn)
             lowerbound=xx(1)
             scanI = .true.
          endif
  ! ---------------------------------------------------------------
  !!! CID PARAMETERS !!! 
  ! ---------------------------------------------------------------
  ! Collision energy (eImpact)     
          !> LAB frame (DEFAULT: OFF)
          if(index(line,'ELAB') /= 0)then            
             call readl(line,xx,nn)
             ELAB=xx(1)
          endif
          !> COM frame (DEFAULT: ON)
          if(index(line,'ECOM') /= 0)then            
             call readl(line,xx,nn)
             ECOM=xx(1)
          endif
          ! Set a maximum number collisions till fragmentation 
          ! This is for the CollAuto run, i.e. collision until
          ! fragmentation with maximum of 'SetColl' collisions
          if(index(line,'SETCOLL') /= 0)then           
             call readl(line,xx,nn)
             coll%set_coll=xx(1)
          endif
          !  Vary the different collision number parameters
          if(index(line,'COLLNO') /= 0)then           
            Manual = .True.
            call readl(line,xx,nn)
            CollNo(1)=xx(1)
            CollNo(2)=xx(2)
            CollNo(3)=xx(3)
          endif
          ! collisions until CollSec fragmentations       
          if(index(line,'COLLSEC') /= 0)then           
            Manual = .True.
            call readl(line,xx,nn)
            CollSec(1)=xx(1)
            CollSec(2)=xx(2)
            CollSec(3)=xx(3)
           
            ! Always enable one fragmentation
            if (CollSec (1) == 0) CollSec(1) = 1
            if (CollSec (2) == 0) CollSec(2) = 1
            if (CollSec (3) == 0) CollSec(3) = 1
          endif
          ! maximum number of Collisions (in FullColl)
          if(index(line,'MAXCOLL') /= 0)then            
             call readl(line,xx,nn)
             coll%max_coll=xx(1)
             Manual = .True.
          endif
          ! Minimum e-field potential for re-acceleration 
          if(index(line,'MINPOT') /= 0)then            
             call readl(line,xx,nn)
             MinPot=xx(1)
          endif
          ! scaling velocity after collision
          if(index(line,'VSCALE') /= 0)then            
             call readl(line,xx,nn)
             vScale=xx(1)
          endif
          ! set ESI internal scaling energy
          if(index(line,'ESI') /= 0)then            
             call readl(line,xx,nn)
             ESI=xx(1)
          endif
          ! set temperature for ESI internal energy
          ! (only choose one of the two)
          if(index(line,'TSCALE') /= 0)then            
             call readl(line,xx,nn)
             tempESI=xx(1)
          endif
  
          ! CollGas Pressure
          if(index(line,'PGAS') /= 0)then            
             call readl(line,xx,nn)
             cell%PGas=xx(1)
          endif
  
          ! CollGas Temperature
          if(index(line,'TGAS') /= 0)then            
             call readl(line,xx,nn)
             cell%TGas=xx(1)
          endif
  
          ! CollGasChamber length
          if(index(line,'LCHAMB') /= 0)then            
             call readl(line,xx,nn)
             cell%lchamb=xx(1)
          endif
  
          ! minimum mass when calcs are stopped
          if(index(line,'MINMASS') /= 0)then            
             call readl(line,xx,nn)
             minmass=xx(1)
          endif
  
          ! set number of steps until collision (circa)
          if(index(line,'DIST') /= 0)then            
             call readl(line,xx,nn)
             manual_dist=int(xx(1))
          endif
  
        endif ! end if CHECK
  
     enddo ! end read loop 
     close(io_in)
   endif  ! End read input (qcxms.in)

   ! Set CID to FullAuto if nothing else was specified
  if (.not. Manual .and. .not. CollAuto .and. .not. TempRun) then
    FullAuto = .true.
  endif
  
  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! Specifics for programs
  
  ! MNDO
  if(prog == 5)then
    ax=0.6
    if(ihamilt == 6)then
    ! fitted for OM2 on S66x8, MAD(S66)=0.77
      a1=0.690
      s8=0.531
      a2=3.446            
      ax=0.8
    endif
    if(ihamilt == 8)then
    ! fitted for OM3 on S66x8, MAD(S66)=0.73
      a1=0.613
      s8=0.501
      a2=3.258            
      ax=0.8
    endif
  
  ! PM3
    if(ihamilt == 2)then
  ! fitted on S66x8 no HB: MAD(S66x8 noHB)=1.15 kcal/mol, MAD(full S77x8)=1.85     
      a1=0
      a2=6.381
      s8=1.668
    endif
  endif
  
  ! DFTB
  if(prog == 0)then
  ! new (10/13) from Gerit:
    a1 = 0.5719
    s8 = 0.5883
    a2 = 3.6017
    ax=0
     write(*,'('' DFTB3 HH mod potential           '',l)')hhmod
  endif
  
  ! MSINDO
  if(prog == 4)then
  ! fitted for MSINDO(NDDO) on S66              
    a1=0.0
    s8=0.7
    a2=3.5               
    ax=0.8
  endif
  
  !      if(prog == 3.or.prog == 2)then
  ! pbe12/sv+gcp on S66x8: MAD(S66x8)=1.05 kcal/mol
  !      if(bas == 1.and.func == 2)then
  !         a1=0.117
  !         s8=2.975
  !         a2=7.150             
  !         ax=0.5
  !      endif
  ! pbe12/svx+gcp on S66x8: MAD(S66x8)=0.71 kcal/mol
  !      if(bas == 2.and.func == 2)then
  !         a1=0.000
  !         s8=2.570
  !         a2=7.340             
  !         ax=0.5
  !      endif
  ! pbe12/sv(p)+gcp on S66x8: MAD(S66x8)=0.71 kcal/mol
  !      if(bas == 3.and.func == 2)then
  !         a1=0.000
  !         s8=1.476
  !         a2=6.911             
  !         ax=0.5
  !      endif
  ! pbe38/svx+gcp on S66x8: MAD(S66x8)=0.76 kcal/mol
  !      if(bas == 2.and.func == 1)then
  !         a1=0.000
  !         s8=0.708
  !         a2=6.270             
  !         ax=0.375
  !      endif
  ! pbe0/svx+gcp on S66x8: MAD(S66x8)=0.76 kcal/mol
  !      if(bas == 2.and.func == 0)then
  !         a1=0.144
  !         s8=0.701
  !         a2=5.431             
  !         ax=0.25
  !      endif
  ! pbe0/sv(p)+gcp on S66x8: MAD(S66x8)=0.84 kcal/mol
  !      if(bas == 3.and.func == 0)then
  !         a1=0.148
  !         s8=0.975
  !         a2=5.725             
  !         ax=0.25
  !      endif
  ! pbe/svp+gcp on S66x8: MAD(S66x8)=1.01 kcal/mol
  !      if(bas == 4.and.func == 5)then
  !         a1=0.174
  !         s8=0.674
  !         a2=5.471             
  !         ax=0.0
  !      endif
  ! pbe0/svp+gcp on S66x8: MAD(S66x8)=0.81 kcal/mol
  !      if(bas == 4.and.func == 0)then
  !         a1=0.153
  !         s8=0.707
  !         a2=5.466             
  !         ax=0.25
  !      endif
  ! pbe38/svp+gcp on S66x8: MAD(S66x8)=0. kcal/mol
  !      if(bas == 4.and.func == 1)then
  !         a1=0.140
  !         s8=0.986
  !         a2=5.716          
  !         ax=0.375
  !      endif
  !      endif
  
  ! MOPAC
  ! mopac NDDO: the STRB52 performance is only very weakly dep. on eTemp
  ! but for large values convergence becomes bad     
  
  if(prog == 1)then
     ax=0.60 
  endif
  
  ! set all fock exchanges regardless of basis set
  if(prog == 2 .or. prog == 3)then
     if (func == 0)  ax = 0.25  !PBE0
     if (func == 1)  ax = 0.375 !PBE38
     if (func == 2)  ax = 0.5   !PBE12
     if (func == 3)  ax = 0.0   !LDA12 (unused)
     if (func == 4)  ax = 0.54  !M062X
     if (func == 5)  ax = 0.0   !PBE
     if (func == 6)  ax = 0.0   !B97-D3
     if (func == 7)  ax = 0.20  !B3LYP
     if (func == 8)  ax = 0.28  !PW6B95
     if (func == 9)  ax = 0.20  !B3PW91
     if (func == 10) ax = 0.0   !BLYP
     if (func == 11) ax = 0.0   !BP86
     if (func == 12) ax = 0.0   !TPSS
     if (func == 13) ax = 0.0   !REVPBE
  endif
  
  ! overwrite ax if user specified in input
  if(axi.gt.-0.001) ax=axi
  
!  write(*,*)
!  call qcstring(prog,line) 
!  write(*,'('' QC method: '',a)')trim(line)               
!  write(*,'('' "Fock"-exchange ax ='',F6.3)')ax
!  
!  if(a2 /= 0)then
!     write(*,*)
!     write(*,'('' D3(BJ)   a1='',F6.3)')a1
!     write(*,'('' D3(BJ)   a2='',F6.3)')a2
!     write(*,'('' D3(BJ)   s8='',F6.3)')s8
!  endif   
  
  if(tun)write(*,'('' tunneling correction switched on'')')
  
  ! DFTB3 OH/CH bond correction parameters, cab > 0 means attractive pot.      
  if(prog == 0.and.tbcorr)then
     cab(1,8)= 0.04   
     cab(1,6)=-0.007
     do i=1,100
        do j=1,i
           cab(i,j)=cab(j,i)
        enddo   
     enddo   
     write(*,'('' DFTB3 correction potential, alpha'',F9.4)')alpc
     do i=1,100
       do j=1,i
          if(abs(cab(i,j)).gt.1.d-6)then
             write(*,'('' DFTB3 correction potential kAB'',2i3,'' ='',F9.4)')i,j,cab(i,j)
          endif
       enddo      
     enddo
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! correct for nonsense options
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ( prog /= 3 .and. prog /= 2 ) gcp = .false.      
  
  if ( gcp ) then
  ! if basis = TZVP or QZVP
     if ( bas == 5 ) then
        write(*,'(''* gcp turned off *'')')
        gcp=.false. 
     endif
     if ( bas == 13 ) then
        write(*,'(''* gcp turned off *'')')
        gcp=.false. 
     endif
  ! if basis = any minimally augmented (SVP,TZVP,TZVPP,QZVP)
     if ( bas == 6 .or. bas == 7 .or. bas == 8 .or. bas == 14 ) then
        write(*,'(''* gcp turned off *'')')
        gcp=.false.
     endif
  ! bas=11 -> def2-TZVP , bas=12 -> def2-QZVP         
     if ( bas == 11 .or. bas == 12 ) then
        write(*,'(''* gcp turned off *'')')
        gcp=.false.
     endif
  endif
  
  
end subroutine input
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

subroutine read_struc_commandline(mol, check, prod, noeq, eonly0, eonly1, &
    eonly, input_file)


  integer :: iarg, narg
  integer, allocatable :: input_format

  character(len=:), allocatable :: arg
  character(len=:), allocatable :: input_file
  logical :: check,eonly,eonly0,eonly1,noeq, prod
  logical :: ex = .false.

  type(error_type), allocatable :: error
   type(structure_type) :: mol

  iarg = 0
  prod = .false.
  narg = command_argument_count()


  do while(iarg < narg)
    iarg = iarg + 1
    call get_argument(iarg, arg)

    select case(arg)
    case("--check","-c")
      check = .true.
    case("--prod","-p")
      prod = .true.
    case("--verbose","-v")
      verbose = .true.
    case("-e0")
      eonly0 = .true.
    case("-e1")
      eonly1 = .true.
    case("-eonly")
      eonly = .true.
    case("-qcp", "--qcpath")
      iarg = iarg + 1
      call get_argument(iarg, arg)
      if (.not.allocated(arg)) then
        write(*,*)
        write(*,'(60(''-''))')
        write(*,'("No QC-Path provided! - Using default: ", a)') &
        & , path
        write(*,'(60(''-''))')
        write(*,*)
        exit
      end if
      path = arg 

    case("-i", "--input")
      iarg = iarg + 1
      call get_argument(iarg, arg)
      if (.not.allocated(arg)) then
        write(*,*)
        write(*,'(60(''-''))')
        write(*,'("No input file provided")')
        write(*,'(60(''-''))')
        write(*,*)
        stop
      else
        input_format = get_filetype("."//arg)
        call read_structure(mol, arg, error, input_format)
        input_file = arg
      end if

    case default
      error stop ' -- Unrecognized Keyword -- '
    end select
  enddo

  !> If production run, search for start.xyz
  if (mol%nat < 1 .and. prod ) then
    call read_structure(mol, 'start.xyz', error, filetype%xyz)

    if (mol%nat < 1) then
      inquire(file='start.xyz', exist=ex)

      if (ex) error stop ' -- something wrong in start.xyz - aborting! -- '

      if (.not. ex) then
        write(*,*) ' -- no reasonable molecule found (searched for start.xyz)! --'
        error stop ' -- provide strucutre with qcxms -i <structure.xyz> -- ' 
      endif

    endif
     
    if(mol%nat > 10000) error stop ' -- too many atoms! (exceeding 10000 atms) --'


  !> If no input is provided, check for coord  
  elseif (mol%nat < 1 .and. .not. prod ) then


    call read_structure(mol, 'coord', error, filetype%tmol )

    !>> check if input provides reasonable molecule
    if (mol%nat < 1) then
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '--- No structure provided! ---'
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) 'Provide structure with -i flag or &
      (turbomole-style) coord file and check in number of atoms is correct'

      write(*,*) 
      error stop 'no reasonable molecule (coord) found! - EXIT'
    endif

    if(mol%nat > 10000) error stop 'too many atoms'

  endif

end subroutine read_struc_commandline
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)                                         
! Original author: Clive Page 

!implicit none

character(len=*), intent(in) :: strIn
character(len=len(strIn)) :: strOut
integer :: i,j

do i = 1, len(strIn)
   j = iachar(strIn(i:i))
   if (j>= iachar("a") .and. j<=iachar("z") ) then
      strOut(i:i) = achar(iachar(strIn(i:i))-32)
   else
      strOut(i:i) = strIn(i:i)
   end if
end do

end function to_upper

end module qcxms_input
