subroutine cid( nuc, iat, mass, xyz, velo, time_step, mchrg, etemp, &
    stopcid, ELAB, ECOM, axyz, ttime, eExact, ECP, manual_dist,     &
    vScale, MinPot, ConstVelo, cross,                               &
    mfpath, r_mol, achrg, icoll, collisions, direc, velo_cm,        &
    aTlast, calc_collisions,imass)

  use common1 
  use cidcommon
  use rmsd_ls, only : get_rmsd
  use qcxms_analyse, only: avg_frag_struc 
  use qcxms_boxmuller, only: vary_energies
  use qcxms_cid_rotation
  use qcxms_fragments
  use qcxms_iniqm, only: iniqm
  use qcxms_mdinit, only: ekinet
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_convert
  use xtb_mctc_constants, only:pi,kB
  use xtb_mctc_symbols, only: toSymbol 
  implicit none
  
  ! note that xyz from main code is in bohr
  real(wp),parameter :: kB_eV  = kB*autoev
  real(wp),parameter :: kB_J   = 1.38064852E-23
  
  integer  :: nuc,spin,iat(nuc)
  integer  :: nstp
  integer  :: dumpavg,dumpxyz,dumpscreen,dumpcoord,dumpdist
  integer  :: average_dump,xyzavg_dump,screen_dump,coord_dump,distance_dump
  integer  :: i,j,k,ind,m
  integer  :: time_step_count
  integer  :: icoll
  integer  :: mchrg
  integer  :: collisions
  integer  :: step_counter
  integer  :: gsolvstate
  integer  :: nuc0,ntot
  integer  :: xtra
  integer  :: list(nuc),nfrag
  integer  :: fragat(200,10)  
  integer  :: imass(nuc)
  integer  :: istep, step_dist, manual_dist
  integer  :: io_cid, io_rotate, io_test
  integer  :: morestep , fconst

  !
  integer :: iatf(nuc,10)
  integer :: idum(nuc,10)
  real(wp) :: xyzf(3,nuc,10)
  real(wp) :: dum (3,nuc,10)
  integer :: natf(10)
  integer :: cnt 
  integer :: check_fragmented


  !

  real(wp) :: calc_collisions
  real(wp) :: etemp,E,ke
  real(wp) :: xyz(3,nuc),velo(3,nuc),grad(3,nuc)
  real(wp), intent(out) :: axyz(3,nuc)
  real(wp) :: avxyz(3,nuc)
  real(wp) :: avxyz2(3,nuc)
  real(wp) :: store_avxyz(3,nuc)
  real(wp) :: ttime
  real(wp) :: mass(nuc)
  real(wp) :: new_velo,velo_diff,velo_cm
  real(wp) :: achrg(nuc) 
  real(wp) :: time_step
  real(wp) :: Edum
  real(wp) :: cm(3)
  real(wp) :: highestCOM,lowestCOM
  real(wp) :: old_cm(3)
  real(wp) :: distcrit !distance criterion
  real(wp) :: start_dist,new_dist
  real(wp) :: Tinit
  real(wp) :: E_velo,E_kin_diff, E_Therm
  real(wp) :: E_COM_calc,beta !,ny,yamma
  real(wp) :: Ekin,T,Tav,avgT
  real(wp) :: new_temp
  real(wp) :: T_GBSA
  real(wp) :: r_mol
  real(wp) :: mfpath
  real(wp) :: Eimpact, ECOM, ELAB
  real(wp) :: E_Scale,MinPot,vScale,W,scale_MinPot
  real(wp) :: summass,cross
  real(wp) :: direc(3),scale_velo(3)
  real(wp) :: fasti
  real(wp) :: diff_cm(3),cm_out
  real(wp) :: normalize_velo
  real(wp) :: xyzAr(3),xyz_start(3) 
  real(wp) :: lowestx,highestx,lowesty,highesty,lowx,lowy,highx,highy,diff1,diff2
  real(wp) :: fragm(10)     
  real(wp) :: E_int(10)
  real(wp) :: fragT(10)     
  real(wp) :: E_Distr
  real(wp) :: aTlast

  real(wp) :: save_grad!(5000)
  real(wp) :: set_grad
  
  !Rotation parameters
!  real(wp) :: a,b,c
  real(wp) :: d,f,g,lmin,lpos
  real(wp) :: velo_rot(3,nuc), E_Rot

  ! Allocatables
  real(wp), allocatable :: xyz0(:,:),velo0(:,:),grad0(:,:)
  real(wp), allocatable :: mass0(:),achrg0(:),aspin0(:)
  integer,  allocatable :: iat0(:)


  !Characters
  character(len=20) :: fname
  character(len=80) :: solvent
  character(len=80) :: fragf(nuc)     
  
  !Logicals
  logical :: ConstVelo
  logical :: ECP
  logical :: iniok
  logical :: stopcid
  logical :: eExact
  logical :: gradfail
  logical :: xstps = .false.
  logical :: avg_struc 

  interface
    function calc_ECOM(beta,e_kin) result(E_COM)
      !use cidcommon
      use xtb_mctc_accuracy, only: wp
      !use xtb_mctc_convert
      implicit none
    
      real(wp) :: beta, e_kin, E_COM
    end function calc_ECOM
  end interface

  ! initiate random numbers
  if (iseed(1) == 0) then
    call random_seed()
  endif

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Allocate molecule+Impactatom
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(gas%Iatom /= 6)nuc0 = nuc+1
  if(gas%Iatom == 6)nuc0 = nuc+2 !N2, but has to be fixed for more atoms
  allocate(xyz0(3,nuc0))
  allocate(velo0(3,nuc0))
  allocate(grad0(3,nuc0))
  allocate(mass0(nuc0))
  allocate(iat0(nuc0))
  allocate(achrg0(nuc0))
  allocate(aspin0(nuc0))
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !OPEN FILES
  if (icoll == 1)then
  
     open(file='rotate.xyz',newunit=io_rotate,status='replace')
     write(io_rotate,*) nuc
     write(io_rotate,*) ' '
     do i=1,nuc
        write(io_rotate,*) toSymbol(iat(i)),' ',xyz(1,i)/aatoau &
            ,' ',xyz(2,i)/aatoau &
            ,' ',xyz(3,i)/aatoau
     end do

   endif

  fname='CID.xyz'

 ! elseif(icoll > 1)then
 !                          write(fname,'(''CID'',i4,''.xyz'')')icoll
  if(icoll < 1000)write(fname,'(''CID'',i3,''.xyz'')')icoll
  if(icoll < 100) write(fname,'(''CID'',i2,''.xyz'')')icoll
  if(icoll < 10)  write(fname,'(''CID'',i1,''.xyz'')')icoll

  open(file=fname,newunit=io_cid,status='replace')
  
  !open(file='avgxyz.xyz',newunit=io_test,status='replace')
  
  
  !---------------------------------------------------------
  !logical
  
  iniok      = .True.
  stopcid    = .False.
  gradfail   = .False.
  avg_struc =  .False.
  
  !------------------------------------------------------     
  
  !------------------------------------------------------     
  !Initial hardcoded parameter block
  !start_dist = 20 * aatoau !distance for Ar-M 
  
  dumpscreen = 100   !interval for screen dumping
  dumpcoord  = 4     !interval for coordinate dumping
  dumpdist   = 10    !interval for distance dumping
  dumpavg    = 50
  ntot       = 15000 !maximum number of steps 
  
  new_velo   = 0.0d0
  Tinit      = 0.0d0
  
  time_step_count = 0
  step_counter = 0
  
  avxyz = 0
  avxyz2 = 0

  cnt = 0


  set_grad = huge(0.0_wp)

  !velo_rot(3,3) = 0.0d0
  
  T_GBSA=300
  solvent='h2o'
  gsolvstate=0
  !------------------------------------------------------     
  nfrag = 1 
  !!!!!!!!!!!
  !> count steps after frag
  morestep = 0
  fconst = 0
  check_fragmented = 1
  !!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! set eletronic temperature 
  if(etemp <= 0)then    ! Fermi-smearing levels
    etemp = 5000.0d0   !+ 20000.0d0 
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! rotate molecule and give angular momentum
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (icoll == 1)then

    !> Calculate center of mass
    call cofmass(nuc,mass,xyz,cm)

    !> translate center of mass to origin
    xyz(1,:) = xyz(1,:)-cm(1)
    xyz(2,:) = xyz(2,:)-cm(2)
    xyz(3,:) = xyz(3,:)-cm(3)


    !> Do Euler Rotation
    call euler_rotation(nuc, iat, xyz, velo, io_rotate)
  
    !> Give angular momentums 
    call rotation_velo(xyz, nuc, mass, velo, velo_rot, E_rot)
  
    ! Calc. ratio of energies (not important)
    !jpoint(1)=  tinit*0.5*3*nuc* kB
    !jpoint(2) = E_Rot
    !jpoint(3) = jpoint(1)/jpoint(2)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i =1,nuc
      velo(1,i) = velo(1,i)  + velo_rot(1,i)
      velo(2,i) = velo(2,i)  + velo_rot(2,i)
      velo(3,i) = velo(3,i)  + velo_rot(3,i)
    end do
  
  endif

  !########################################################################
  ! Set starting conditions
  !########################################################################

  !> Calculate masses
  summass = 0.0d0
  do i = 1,nuc
     summass = summass +  mass(i) 
  enddo

  beta  = gas%mIatom/(gas%mIatom + summass)
  
  !> set the kinetic energy 
  !>> For first collision 
  if (icoll == 1)then

    !> determine if start E is LAB or COM 
    if ( ELAB > 0.0_wp ) then
      Eimpact = ELAB
    else
      Eimpact = ECOM / beta
    endif

    !1. Get Impact Energy! Than vary it via box-muller distribution
    if (.not. eExact) then
      E_Distr = 0.1_wp
      E_velo = vary_energies(Eimpact, E_Distr) * evtoau
    !2. Or set Impact Energy as given value
    else
      E_velo = Eimpact * evtoau !set Energy constant in all prod runs 
    endif

  !>> For consecutive collisions
  else
    ! set this, because it is set in main to au units 
    E_velo = 0.5_wp * summass * ((velo_cm * mstoau)**2) 
  endif


  !> calculate thermal energy molecule 
  call ekinet(nuc,velo,mass,Ekin,Tinit)
  E_Therm = tinit * 0.5_wp * 3.0_wp * nuc * kB
  E_COM_calc = calc_ECOM(beta,E_velo) ! beta * E_velo 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Stuff for scaling velos
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !ny    = 2 
  !yamma = ((2 * summass) + (ny * gas%mIatom)) / (summass + gas%mIatom)
  


  !write(*,*) 'BETA',beta
  !write(*,*) 'GAMMA',yamma
 

  !########################################################################
  !### BEGIN CID MODULE  ################################################### 
  !########################################################################
  
  if (icoll == 1)then
     !set Energy and velocity
     fasti = sqrt(2 * E_velo / summass)

     write(*,*) ' '
     write(*,'(80(''=''))')
     write(*,*) '    CID settings:   '
  
     if ( mchrg > 0 ) then
       write(*,*) '    + Positive Ion mode +   '
     elseif ( mchrg < 0 ) then
       write(*,*) '    - Negative Ion mode -   '
     endif
  
     write(*,'(80(''=''),/)')
     write(*,*) ' '
     write(*,*) '----- Precursor and Gas weights --------------'
     write(*,'('' Collision Gas                :'',  A8)'  )  toSymbol(gas%IndAtom)
     write(*,'('' Number of Precursor atoms    : '',i7)'  ) nuc
     write(*,'('' Mass Collsion Gas  (Da)      : '',f14.6)') gas%mIatom * autoamu
     write(*,'('' Mass Precursor Ion (Da)      : '',f14.6)') summass* autoamu 
!     write(*,'('' Mass Precursor Ion (kg)      : '',e10.2)')summass * amutokg
     write(*,*)
     write(*,*) ' ----- Collision Energetics  -----------------'
     write(*,'('' Kinetic Energy eV (LAB)      : '',f14.6)')  E_velo / evtoau
     write(*,'('' Kinetic Energy au (LAB)      : '',f14.6)')  E_velo 
     write(*,'('' Kinetic Energy eV (COM)      : '',f14.6)')  E_COM_calc
     write(*,'('' Kinetic Energy au (COM)      : '',f14.6)')  E_COM_calc * evtoau
     write(*,*) ' '
     write(*,'('' Velocity (au)                : '',e14.2)')  fasti
     write(*,'('' Velocity (m/s)               : '',f14.6)')  fasti / mstoau 
     write(*,*) ' '
     write(*,*) ' ----- Molecular ion properties --------------'
     write(*,'('' Molecular temp. (K)          : '',f14.6)') Tinit
     write(*,'('' Rotational Energy (eV)       : '',f14.6)') E_Rot   / evtoau
     write(*,'('' Thermal Energy (eV)          : '',f14.6)') E_Therm / evtoau 
     write(*,'(32x,15(''=''))') 
     write(*,'('' Internal Energy (eV)         : '',f14.6)') (E_Rot + E_Therm) / evtoau 
     write(*,*) ' '
     write(*,*) ' ----- Collision calculations ----------------'
     write(*,'('' Radius Ion (m)               : '',e14.2)') r_mol   
     write(*,'('' Radius Gas Atom (m)          : '',e14.2)') (gas%rIatom*autoaa) * 1E-10 ! m 
     write(*,'('' Cross-section (m²)           : '',e14.2)') cross  
     write(*,'('' Temperature coll.Gas (K)     : '',f14.6)') cell%TGas
     write(*,'('' Pressure coll.Gas (Pa)       : '',f14.6)') cell%PGas
     write(*,'('' Collision cell length (m)    : '',f14.6)') cell%lchamb
     write(*,'(32x,15(''-''))')
     write(*,'('' Mean free path (m)           : '',f14.6)') mfpath
     write(*,'('' Calculated collisions        : '',f14.6)') calc_collisions 
     write(*,'(32x,15(''=''))')
     write(*,'('' Number of collisions varied  : '',i7)'   ) collisions 
     write(*,*)
     write(*,'(47(''-''),/)')
  else
     write(*,*) '----- Ion properties -------------------------'
     write(*,'('' Number of atoms              : '',i7)'  ) nuc
     write(*,'('' Mass collision Ion           : '',f14.6,a3)') summass/amutoau, ' Da'
     write(*,'('' Mass Collsion Gas            : '',f14.6,a3)') gas%mIatom/amutoau,  ' Da'
     write(*,*) ' '
     write(*,*) ' ----- Collision Energetics  -----------------'
     write(*,'('' Kinetic Energy (LAB)         : '',f14.6,a3)') E_velo / evtoau,  ' eV'
     write(*,'('' Kinetic Energy (COM)         : '',f14.6,a3)') E_COM_calc     ,  ' eV'
     write(*,'('' Velocity now                 : '',f14.6,a4)') velo_cm        , ' m/s'
     write(*,*) ' '
     write(*,*) ' ----- Collision calculations ----------------'
     write(*,'('' Radius Ion                   : '',e14.2,a2)') r_mol   ,' m'
     write(*,'('' Radius Gas Atom              : '',e14.2,a2)') (gas%rIatom/aatoau) * 1E-10 , ' m'
     write(*,'('' Cross-section                : '',e14.2,a4)') cross ,'m²' 
     write(*,'(32x,15(''-''))')
     write(*,'('' Mean free path               : '',f14.6,a2)') mfpath, ' m'
     write(*,'('' Calculated collisions        : '',f14.6)') calc_collisions 
     write(*,'(32x,15(''=''))')
     write(*,'('' Number of collisions varied  : '',i7)'   ) collisions 
     write(*,*) ' '
     write(*,'(32x,15(''-''))')

     if (MinPot > 0 .and. E_COM_calc < MinPot)then
       scale_MinPot = MinPot / beta
       W = mfpath * scale_MinPot 
       vScale = sqrt((2*W)/summass)
       write(*,*) 'scale MinPot   :',scale_MinPot , 'au'
       write(*,*) 'scale MinPot eV:',scale_MinPot /evtoau , 'eV'
       write(*,*) 'scale   :',vScale/evtoau , 'eV'
       write(*,*) 'scale   :',vScale/mstoau , 'm/s'
       vScale = abs(vScale-(velo_cm*mstoau))
       write(*,*) 'scale diff   :',vScale/mstoau , 'm/s'
       write(*,*) ' '
     endif
  endif
  
  if (ConstVelo .and. icoll  >  1)then
     W = sqrt(2*E_velo/summass) !take the FIRST energy
     W = W - (velo_cm * mstoau)
     write(*,*) 'W scale m/s',velo_cm+ W/mstoau

     if (W  >  0)then
       vScale= W
       write(*,*) 'Scale Velocity!'
     endif
  
  
     write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) ' ATTENTION: Velocity scaled towards constant Velocity'
     write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  endif
  write(*,*) ' '
  
  write(*,*) '============================ '
  write(*,*) 'Starting Collision Dynamics: '
  write(*,*) '============================ '
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Argon placement and initial velocity determination
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
  call cofmass(nuc,mass,xyz,cm)
  
  
  !if (icoll == 1)then
  !   call random_number(a)
  !   call random_number(b)
  !   call random_number(c)
  !endif
  call random_number(d)
  call random_number(f)
  call random_number(g)
  call random_number(lmin) ! Pos/Neg radomizer
  call random_number(lpos) ! Pos/Neg radomizer
  
  !! Randomize the xyz direction in which the molecule is accelerated to ensure different collision angles 
  !   a=0
  !   b=0
  !   c=1  ! set to z-axis, because Euler-rotation should account for collision angles
  !!   write(*,*)  xyz
  !if    (a >= b.and.a >= c)then
  !   a=1
  !   b=0
  !   c=0
  !   g=0
  !   ! Random d,f
  !   if(lmin < 0.5)f= (f * -1.0d0)
  !   if(lpos < 0.5)d= (d * -1.0d0)
  !   write(*,*) 'X-Axis'
  !elseif(b >= a.and.b >= c)then
  !   a=0
  !   b=1
  !   c=0
  !   f=0
  !   ! Random d,g
  !   if(lmin < 0.5)d= (d * -1.0d0)
  !   if(lpos < 0.5)g= (g * -1.0d0)
  !   write(*,*) 'Y-Axis'
  !elseif(c >= a.and.c >= b)then
  !   a=0
  !   b=0
  !   c=1
  !   d=0 
  !   ! Random f,g
  !   if(lmin < 0.5)f= (f * -1.0d0) 
  !   if(lpos < 0.5)g= (g * -1.0d0) 
  !   write(*,*) 'Z-Axis'

  !Vary the collision angle depending on the maximum radius of the molecule
  lowestx  =  huge(0.0_wp)
  lowesty  =  huge(0.0_wp)
  highestx = -huge(0.0_wp)
  highesty = -huge(0.0_wp)
  do i=1,nuc
    lowx  = xyz(1,i)
    lowy  = xyz(2,i)
    highx = xyz(1,i)
    highy = xyz(2,i)
    if (lowx  < lowestx)  lowestx  = lowx
    if (highx > highestx) highestx = highx
    if (lowy  < lowesty)  lowesty  = lowy
    if (highy > highesty) highesty = highy
  enddo
  
  
     if(lmin < 0.5)then
       diff1 = (lowesty ) * f 
     else
       diff1 =  (highesty )* f 
     endif
     if(lpos < 0.5)then
       diff2 = (lowestx )* g 
     else
       diff2 = highestx * g 
     endif
  
     if (g > 0.85) f = f*0.5 !reduce the amount of x-axis if y is very large
     if (f > 0.85) g = g*0.5
  !endif
 

  
  if(icoll == 1)then

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> set (step) distance automatically
    if ( manual_dist == 0 ) then
      step_dist  =   2 * nuc * 10            ! set number of steps until collision depending on the size of the fragment/ion
      if (step_dist > 800 ) step_dist = 800  ! set upper limit to reduce cost
    !> or manually
    else
      step_dist = manual_dist
    endif

    start_dist = ( fasti ) * ( 2 * time_step ) ! convert velo. into distance
    start_dist =  step_dist * start_dist * autoaa ! re-convert into angstrom

    if ( start_dist < 10.0_wp ) start_dist = 10.0_wp ! set lower limit to account for too short distances
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Set the initial direction of the hit depending on the axis-of-flight
    xyz_start(1) = cm(1) !+ start_dist * a 
    xyz_start(2) = cm(2) ! + start_dist * b 
    xyz_start(3) = cm(3) + start_dist 
  
    direc(:) = xyz_start(:) 
  
    !change direc to unity vector
    do i = 1,3
       direc(i) = direc(i) / (sqrt(direc(1)**2 + direc(2)**2 + direc(3)**2))
    end do
  
  
  !  write(*,*) 'Radius', rtot/aatoau
    ! Set the collision gas atom away from the COM by random factors depending on the axis-of-flight
    xyzAr(1) = xyz_start(1) + diff2 *0.8 !(g * rtot)*0.6 
    xyzAr(2) = xyz_start(2) + diff1 *0.8 !(f * rtot)*0.6
    xyzAr(3) = xyz_start(3) !+  !(d * rtot)*0.6
  !  write(*,*) 'Ar XYZ', xyzAr/aatoau
  !  write(*,*) 'randoms g=x',  g 
  !  write(*,*) 'randoms f=y',  f 
  !  write(*,*) 'randoms d=z',  d 
  
    scale_velo = (direc * fasti)
  
  else
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> set (step) distance automatically
    if ( manual_dist == 0 ) then
      step_dist  =   2 * nuc * 10            ! set number of steps until collision depending on the size of the fragment/ion
      if (step_dist > 800 ) step_dist = 800  ! set upper limit to reduce cost
    !> or manually
    else
      step_dist = manual_dist
    endif

    ! Make collision distance dependend on mol. velo.
    start_dist = ( velo_cm * mstoau ) * ( 2 * time_step ) ! convert velo. into distance
    start_dist =  step_dist * start_dist * autoaa ! re-convert into angstrom

    if ( start_dist < 17.0_wp ) start_dist = 17.0_wp ! set lower limit to account for too short distances
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Re-initialize the coordinates
    xyz(1,:) = xyz(1,:) - cm(1)
    xyz(2,:) = xyz(2,:) - cm(2)
    xyz(3,:) = xyz(3,:) - cm(3)
    ! Set the collision gas atom away from the COM by random factors depending on the axis-of-flight
    xyzAr(1) = cm(1) + (direc(1) * start_dist) + diff2 * 0.7 !(g * rtot)*0.90 
    xyzAr(2) = cm(2) + (direc(2) * start_dist) + diff1 * 0.7 !(f * rtot)*0.90 
    xyzAr(3) = cm(3) + (direc(3) * start_dist) !+ (d * rtot)*0.90 
  !  xyzAr(:) = cm(:) + direc * start_dist 
    scale_velo = 0
  endif
  
  
  ! Scale velocities for fragmentation product
  if (icoll > 1 .and. vScale > 0)then
     if (MinPot > 0 .or. ConstVelo)then
       scale_velo=vScale * direc !*mstoau
     else
       scale_velo=vScale * mstoau
     endif
  
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  old_cm(:) = cm(:)
  
  call cofmass(nuc,mass,xyz,cm)
  
  ! Set the molecule, atom, grads and velos        
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i=1,nuc
     velo0(1:3,i) =  velo(1:3,i) + scale_velo 
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  xyz0(1:3,1:nuc)   = xyz(1:3,1:nuc)
  grad0(1:3,1:nuc)  = grad(1:3,1:nuc)
  mass0(1:nuc)      = mass(1:nuc)
  iat0(1:nuc)       = iat(1:nuc)
  
  xyz0(1:3,nuc0)    = xyzAr(1:3) 
  velo0(1:3,nuc0)   = 0.0d0 
  grad0(1:3,nuc0)   = 0.0d0
  mass0(nuc0)       = gas%mIatom
  iat0(nuc0)        = gas%IndAtom

  if ( gas%Iatom ==  6 ) then !N2 collision gas
     xyz0(1:2,nuc0-1)   = xyzAr(1:2)
     xyz0(3,nuc0-1)     = xyzAr(3)+ 1.09 * aatoau
     velo0(1:3,nuc0-1)  = 0.0d0 
     grad0(1:3,nuc0-1)  = 0.0d0
     mass0(nuc0-1)      = gas%mIatom
     iat0(nuc0-1)       = gas%IndAtom
  endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Initialize the QM code
  !It is done within cid.f90 - to avoid declaring all variables inside main
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!  if (lgbsa)then
!     write(*,*) 'Solvation model switched on!'
!  !   call ekinet(nuc0,velo0,mass0,edum,t)
!     call init_gbsa(nuc0,iat0,solvent,gsolvstate,T_GBSA)
!  endif

  if ( prog == 3 ) then
    call system('rm ORCA.INPUT.*')
  endif
  
  spin = 0
  call iniqm(nuc0,xyz0,iat0,mchrg,spin,etemp,edum,iniok,ECP)
  if(.not.iniok) then
     iniok=.false.
     call iniqm(nuc0,xyz0,iat0,mchrg,spin,etemp,edum,iniok,ECP)
     if(.not.iniok)then
  
        write(*,*) 'Attention:'
        write(*,*) '--------------------------------------------'
        write(*,*) 'There was a problem with the QC code        '
        write(*,*) 'initialization.'
        write(*,*) '--------------------------------------------'
        stopcid = .True.
        stop 'E X I T - because of failure'
     endif
  endif
  
  
  !Initialize the Epot and forces
  call egrad(.True.,nuc0,xyz0,iat0,mchrg,spin,etemp,E,grad0, &
             achrg0,aspin0,ECP,gradfail)
  
  write(*,*)
  write(*,'('' Est. no. steps (collision): '',i4,3x,'' Distance (A): '',f12.6)') &
    & step_dist, start_dist
  write(*,*)
  write(*,'(''   step  time [fs]'',4x,''Epot'',7x,''Ekin'',7x, &
    &    ''Etot'',6x,''eTemp'',4x,''Avg. Temp.'')')
  screen_dump = 0
  coord_dump = 0
  distance_dump = 0
  xyzavg_dump = 0 
  average_dump = 0 
  nstp  = 0
  Tav   = 0.0d0
  m     = 0
  xtra  = 0
 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MD loops
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Calculate the distance between COM and Ar atom; set lowest distance as 
  !criteri for rest calculations
  
  call cofmass(nuc,mass,xyz,cm)
  call distArCOM(nuc0,xyz0,cm,new_dist)
  call minmaxCOM(new_dist,highestCOM,lowestCOM)

  lowestCOM = new_dist
  distcrit = 1.05d0 * lowestCOM  !say at which point the sim ends (by distance)
 
  !> start the leapfrog algorithm loop
  do istep = 1, ntot

    !> re-set xyz_avg structures
    if (xyzavg_dump == dumpavg) then
      xyzavg_dump = 0
      avxyz = 0
    endif

    ! increase the step counters
    ttime = ttime + time_step * autofs 
    screen_dump   = screen_dump   + 1 ! counter for screen dump
    coord_dump    = coord_dump    + 1 ! counter for coord dump
    distance_dump = distance_dump + 1 ! counter for distance criterion
    xyzavg_dump   = xyzavg_dump   + 1 ! counter for average xyz structure
    average_dump   = average_dump   + 1 ! counter for average 
 

    ! Do the leap frog step 
    call leapfrog(nuc0,grad0,mass0,time_step,xyz0,velo0,ke,nstp)

    ! Calculate Grad for entire system
    call egrad(.False.,nuc0,xyz0,iat0,mchrg,spin,etemp,E,grad0,achrg0,aspin0,ECP,gradfail)
  
    if (gradfail)then !If 'grad seems to be bogus!' in egrad
       stopcid=.true.
       exit
    endif

    ! Calculate center-of-mass
    call cofmass(nuc,mass0(:nuc),xyz0(:,:nuc),cm)
  
    diff_cm(:) = cm(:) - old_cm(:)
    cm_out = sqrt(diff_cm(1)**2 + diff_cm(2)**2 +  diff_cm(3)**2)
    old_cm(:) = cm(:)
  
    ! compute new velocity and convert to m/s
    if (istep /= 1)then
      new_velo = (cm_out / time_step ) / mstoau 
    else !i.e. in the first step: 
      new_velo = 0
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate different energies/temperatures 
    call ekinet(nuc,velo0(:,:nuc),mass0(:nuc),Ekin,T)
  
    E_velo = 0.5 * summass * ((new_velo*mstoau)**2) 
    E_kin_diff = Ekin - E_velo
  
    new_temp = (2*E_kin_diff) / (3 * kB * nuc)
    if (istep == 1) new_temp = tinit
    Tav = Tav + new_temp
  
    m    = m + 1
    avgT = Tav / m

    !set collision energy (eV)
    E_COM_calc = calc_ECOM(beta,E_velo)  !(beta * E_velo)/evtoau

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! increase number of steps to allow longer collision MD if fragmentation occured
    call fragment_structure(nuc,iat0(:nuc),xyz0(:,:nuc),3.0d0,1,0,list)
    call fragmass(nuc,iat0(:nuc),list,mass0(:nuc),imass(:nuc),nfrag,fragm,fragf,fragat)
 
     !> get the average structure of the ion/fragments
    avxyz  = avxyz  + xyz0(:,:nuc)

    if (nfrag > check_fragmented) then
      cnt = cnt + 1
      avxyz2  = avxyz2  + xyz0(:,:nuc)
      if (cnt == 50) then
        store_avxyz  = avxyz2 / cnt
        !avxyz2  = avxyz / xyzavg_dump
        check_fragmented = nfrag
        write(*,*) 'Count', cnt
        cnt = 0
        avxyz2 = 0
      endif
    endif

    if(nfrag > 3) then ! nfragexit) then
       !write(*,8000)nstep,ttime,Epot,Ekin,Epot+Ekin,Eerror,nfrag,etemp,fragT(1:nfrag)
       !write(*,9001)
       !fragstate=1
       write(*,*) 'Too many frags'
       exit
    endif   

    ! if fragmented, start counting steps
    if(nfrag >= 2)then
       fconst=fconst+1
    else
       fconst=0          
    endif
    ! exit if nfrag=2 is constant for some time  
    if(fconst > 1000) then
    !   write(*,8000)nstep,nstep*tstep/fstoau,Epot,Ekin,Epot+Ekin,Eerror,nfrag,etemp,fragT(1:nfrag)
    !   write(*,9002)
       write(*,*) 'NOTHING HAPPNENS'
       !fragstate=2       
       exit
    endif   

    ! add a few more cycles because fragmentation can directly proceed further and we don't want to miss this         
    if(nfrag >= 3 ) then !nfragexit) then
       morestep=morestep+1
       if(morestep > 250) then !more) then
          !write(*,8000)nstep,ttime,Epot,Ekin,Epot+Ekin,Eerror,nfrag,etemp,fragT(1:nfrag)
          !write(*,9003)
          write(*,*) 'three or more frags'
          !fragstate=1
          exit
       endif 
    endif        

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! D U M P s
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !DUMP coordinates to CID.xyz and DUMP text to screen                  
    
    if (coord_dump == dumpcoord .or. istep == 1)then
       coord_dump = 0
       write(io_cid,*) nuc0
       write(io_cid,*) 'E=',ke+E
       do ind = 1,nuc0
          write(io_cid,*) trim(toSymbol(iat0(ind))) &
               ,' ',xyz0(1,ind) * autoaa     &
               ,' ',xyz0(2,ind) * autoaa     &
               ,' ',xyz0(3,ind) * autoaa 
       end do        
       if (istep == 1)coord_dump=1
    end if
    
    ! Screen dump
    if (screen_dump == dumpscreen)then
       screen_dump = 0
       
       write(*,'(i7,f8.0,3x,F11.4,2x,F8.4,F11.4,2x,F8.0,F12.3)') &
         &      nstp,ttime,E/evtoau,ke/evtoau,E+ke/evtoau,etemp,avgT
    end if

    if (average_dump == dumpavg)then
      average_dump = 0
      aTlast = 0
    end if

    aTlast = avgT    


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
    !Check whether to STOP the CID module 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (distance_dump == dumpdist)then    ! every 10 steps (=dumpdist)
       distance_dump = 0
    
       ! Calculate the new velocity 
    
       !1.Stop MD loop if shortest distance between Ar and molecule is more 
       !  than initial dist or if after impact a certain time has passed
       call cofmass(nuc,mass0(:nuc),xyz0(:,:nuc),cm)
       call distArCOM(nuc0,xyz0,cm,new_dist)
       call minmaxCOM(new_dist,highestCOM,lowestCOM)
  
       if (lowestCOM < new_dist)then
          step_counter = step_counter + 1
       else 
          step_counter = 0
       endif

       if (step_counter == 10)then !100 steps (see dumpdist)
          write(*,*)
          write(*,'(''COLLISION after '',i6,a6)') nstp,' steps'
  
          time_step_count = nstp + nint(800.0_wp * (2 * time_step * autofs)) 

          write(*,'(''STOP      after '',i6,a6)') time_step_count,' steps'
          write(*,*) 
          Tav = 0  !Start counting again to ensure only the
          m   = 0  !temp after the coll. is taken
       endif

       if (nfrag > 1 .and.step_counter == 50)then !500 steps (we need some time after frag)
          xtra = 600  !should be made dependend on the velocity of fragment
                
          write(*,'('' FRAGMENTATION occured!'')')
          write(*,'('' Do '',i3,a12)') xtra, ' extra steps'


       endif
  
      ! Stop if either the dist-criterion...
       if (lowestCOM > distcrit)then
          stopcid = .False.
          write(*,*) 'The collision MD finished in',nstp,' steps'
          write(*,*) 'Collision MD finished due to DISTANCE!'
  !        write(*,*) ''
  !        write(*,*) 'The LAB velocity:',  new_velo, 'm/s'
          velo_diff = velo_cm - new_velo
          time_step_count = 0
          exit
       endif

     ! ...or the step-criterion is fullfilled
       if (nstp == time_step_count+xtra)then  ! 800 Steps + xtra Steps
          stopcid = .False.
          write(*,*) 'Collision MD finished due to STEPS!'
          write(*,*) ''
          velo_diff = velo_cm - new_velo
          E_velo = (0.5_wp * summass * ((new_velo)**2)) !/ evtoau
          !write(*,*) ''
          !write(*,*) 'Velo Difference' , velo_diff
          !write(*,*) 'Velocity now  :',  new_velo, 'm/s'
          !write(*,*) ''
          time_step_count = 0
          exit
       endif
    endif
  
  
  !2.Stop simulation because the max nsteps are reached. 
     if (istep == ntot)then
        write(*,*) 'Attention:'
        write(*,*) '-----------------------------------------------'
        write(*,*) 'Maximum number of steps has been reached, but  '
        write(*,*) 'probably no collision occured.                 '
        write(*,*) '-----------------------------------------------'
        stopcid = .False.
     endif
  
  
  end do
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Save new coordinates and velocities
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ConstVelo .and. velo_diff  /=  0)then
     normalize_velo = velo_diff * mstoau
  write(*,*) 'Velocity increased by ', normalize_velo/mstoau, 'm/s'
  endif
  
  call ekinet(nuc,velo0(:,:nuc),mass0(:nuc),Ekin,T)
           
  E_velo = 0.5 * summass * ((new_velo*mstoau)**2)
  
!  e_int0 = avgT*(0.5*3*nuc*kB)


  write(*,'(50(''#''))')
  write(*,*) 'Energetics after collision: '
  write(*,*)

  if (nfrag > 1)then
     do j = 1, nfrag
        call intenergy(nuc,list,mass0(:nuc),velo0(:,:nuc),nfrag,fragT,E_int)
        write(*,'('' Internal Energy fragment #'',i2,'' : '',f14.6,a3)') j,E_int(j), ' eV' 
        write(*,'('' Internal Temp fragment #'',i2,'' : '',f14.6,a3)') j,fragT(j), ' eV' 
     enddo
  else
     call intenergy(nuc,list,mass0(:nuc),velo0(:,:nuc),nfrag,fragT,E_int)
     write(*,'('' Internal Energy               : '',f14.6,a3)')E_int(1), ' eV' 
     write(*,'('' Internal Temp                 : '',f14.6,a3)')fragT(1), ' eV' 
  endif
  write(*,'('' Kinetic Energy  (COM)         : '',f14.6,a3)') E_COM_calc           ,' eV' 
  write(*,'(50(''#''))')
  write(*,*)
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Finish due to fragmentation in CID module
  write(*,*) 'E X I T - the CID module is finished'
  write(*,*) ''
  xyz (1:3,1:nuc) = xyz0  (1:3,1:nuc) 
  velo(1:3,1:nuc) = velo0(1:3,1:nuc)
  grad(1:3,1:nuc) = grad0(1:3,1:nuc)
  achrg(1:nuc)    = achrg0(1:nuc)
  velo_cm         = new_velo

  if (check_fragmented > 1 ) then 
    !axyz (1:3,1:nuc) = avxyz2  (1:3,1:nuc) 
    axyz (1:3,1:nuc) = store_avxyz  (1:3,1:nuc) 
  else
   axyz (1:3,1:nuc) = avxyz  (1:3,1:nuc) / xyzavg_dump
  endif

  deallocate(xyz0,velo0,grad0,mass0,iat0,achrg0,aspin0)
  
  close(io_cid)
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Routines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !Calculates center of mass and returns it in variable cm      
  subroutine cofmass(nuc,mass,xyz,cm)
  use xtb_mctc_accuracy, only: wp
  implicit none      

  integer  :: i, nuc

  real(wp) :: totmass
  real(wp),intent(in) :: xyz(3,nuc), mass(nuc)
  real(wp),intent(out) :: cm(3)


  totmass = 0.0d0
  cm = 0.0d0
  do i = 1,nuc
     cm(1) = cm(1) + mass(i) * xyz(1,i)
     cm(2) = cm(2) + mass(i) * xyz(2,i)
     cm(3) = cm(3) + mass(i) * xyz(3,i)
     totmass  = totmass + mass(i)
  end do
  cm(1) = cm(1) / totmass
  cm(2) = cm(2) / totmass
  cm(3) = cm(3) / totmass
  end subroutine
 

  !############################################################################
  !Distance coll atom to COM - N2 doesnt matter
  
  subroutine distArCOM(nuc0,xyz0,cm,new_dist)
  use xtb_mctc_accuracy, only: wp
  implicit none
  integer  :: j
  integer  :: nuc0
  real(wp), intent(in) :: xyz0(3,nuc0),cm(3)
  real(wp) :: dist(3) 
  real(wp),intent(out) :: new_dist
  
  do j=1,3
    dist(j) = xyz0(j,nuc0) - cm(j)
  enddo
  
  new_dist = sqrt(dist(1)**2 + dist(2)**2 + dist(3)**2)
  
  end subroutine
  !############################################################################
  
  !Calculate the  
  subroutine minmaxCOM(new_dist,highestCOM,lowestCOM)
  use xtb_mctc_accuracy, only: wp
  implicit none
  
  real(wp) :: new_dist
  real(wp) :: lowestCOM,highestCOM
  
  if (new_dist>highestCOM) then
     highestCOM = new_dist
  endif
  if (new_dist<lowestCOM)then
     lowestCOM = new_dist
  endif
  end subroutine



!#########################################################################
! Calcualte the molecular radius, cross section and mean-free-path
!#########################################################################
  subroutine collision_setup(nuc,iat,xyz,mass,r_mol,cross,mfpath,calc_collisions)
  use cidcommon
  use covalent_radii, only: Rad
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_convert
  use xtb_mctc_constants, only: pi
  implicit none

  integer  :: nuc
  integer  :: iat(nuc)
  integer  :: i

  real(wp),intent(in)  :: xyz(3,nuc),mass(nuc)
!  real(wp),intent(in)  :: TGas,PGas,lchamb
  real(wp),intent(out) :: r_mol,cross,mfpath
  real(wp) :: cm(3),mol_rad,rtot
  real(wp) :: r_atom,calc_collisions

  real(wp),parameter :: kB_J   = 1.38064852E-23

  mol_rad = 0
  rtot    = 0

  call cofmass(nuc,mass,xyz,cm)

  ! Determine the radii by checking the largest distance from the center of mass
  do i=1,nuc
      mol_rad   = sqrt((xyz(1,i)-cm(1))**2 + (xyz(2,i)-cm(2))**2 + (xyz(3,i)-cm(3))**2) + Rad(iat(i))

      if (mol_rad > rtot) rtot = mol_rad 
  enddo

  ! Calculate radii in meters
  r_atom = (gas%rIatom/aatoau) * 1E-10 ! m
  r_mol  =  (rtot  /aatoau) * 1E-10 ! m

  !write(*,*) 'R TOT', rtot
  !write(*,*) 'Radius Molecule', r_mol

  ! calculate collisional cross section and mean free path
  cross  = pi * ((r_mol + r_atom)**2)
  mfpath = (kB_J * cell%TGas) / (cross * cell%PGas) 

  calc_collisions = cell%lchamb/mfpath

  end subroutine
  

!#########################################################################
! Calculate the Center-of-mass energy
!#########################################################################

function calc_ECOM(beta,e_kin) result(E_COM)
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_convert
  implicit none

  real(wp) :: beta, e_kin, E_COM

  E_COM = (beta * E_KIN) * autoev

end function calc_ECOM
