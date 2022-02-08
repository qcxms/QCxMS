!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
! main MD part of the code with various subtleties (a big bug can
! occur in every line, be carefull!)
! the variable it determines the mode:
!                          -1 : M equilibration      
!                           0 : M sampling           
!                          >0 : M+ fragmentation       (most important)
!                        9999 : M+ fragmentation trials(not used anymore)
!
! mdok is set true if the result should be taken
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
subroutine md(it,icoll,isec,nuc,nmax,xyz,iat,mass,imass,mchrg,grad, &
              velo,velof,list,tstep,totdump,nfragexit,              &
              fragm,fragf,fragat,dumpstep,etempin,                  &
              mdok,achrg,aspin,axyz,                                &
              Tsoll,tadd,eimp,restart,Tav,Epav,Ekav,ttime,aTlast,   &
              fragstate,dtime,ECP,starting_md,new_velo)
  use common1
  use cidcommon
  use qcxms_impact, only: impactscale
  use qcxms_fragments
  use qcxms_mdinit, only: ekinet
  use qcxms_utility, only: setetemp
  use xtb_mctc_accuracy, only: wp            
  use xtb_mctc_convert
  use xtb_mctc_constants
  use xtb_mctc_symbols, only: toSymbol 
  implicit none
  
  integer :: nuc,iat(nuc),list(nuc)
  integer :: nmax,it,k,totdump,icoll,isec
  integer :: dumpstep,nfragexit,fragstate
  integer :: imass(nuc)
  integer :: fragat(200,10) ! 100 elements+100 for isotopes, 10 fragments max per fragmentation      
  integer :: nstep,ndump,mdump,nfrag,i,j,avdump,kdump,vdump
  integer :: screendump,nadd,morestep,more,spin,fconst
  integer :: mchrg
  integer :: io_GS, io_OUT
  integer :: step

  integer :: check_fragmented
  integer :: cnt
  
  real(wp) :: xyz (3,nuc)
  real(wp) :: grad(3,nuc)
  real(wp) :: velo(3,nuc)
  real(wp) :: axyz(3,nuc)
  real(wp) :: velof(nuc)
  real(wp) :: aspin(nuc)
  real(wp) :: achrg(nuc)
  real(wp) :: mass (nuc)
  real(wp) :: fragm(10)     
  real(wp) :: fragT(10)     
  real(wp) :: tstep,tadd,eimp,aTlast,dtime
  real(wp) :: T,Tav,Epav,Ekav,etempin,Tsoll,ttime
  real(wp) :: Epot,Eerror,Ekin,Edum,dum,etemp,Eav,fadd
  real(wp) :: Ekinstart
  real(wp) :: avspin(nuc),avchrg(nuc),avxyz(3,nuc)
  real(wp) ::  sca
  !!! CID Stuff
  real(wp) :: diff_cm(3),cm_out
  real(wp) :: cm(3),old_cm(3)
  real(wp) :: new_velo,new_temp
  real(wp) :: E_kin,E_kin_diff,summass
  real(wp) :: Ekin_new,tinit
  real(wp) :: E_int(10)
!  real(wp), allocatable :: E_int(:)
  real(wp) :: avxyz2(3,nuc)
  
  character(len=20) :: fname
  character(len=80) :: fragf(10)     
  
  logical err1,err2
  logical ECP
  logical gradfail
  logical starting_md
  logical mdok,restart
 
  !!!!!!!!!!!
  !> count steps after frag
  step = 1
  !!!!!!!!!!!
  
  
  write(*,'(/13x''E N T E R I N G   M D   M O D U L E'',/)')
  
  mdok =.false.
  nfrag=1
  ! status of fragments when run was ok (=0 is undefined)
  ! 1: normal
  ! 2: nothing happend for nfrag=2 for some time      
        fragstate=0       
  
  ! it=-1 : equilibration GS, no dump
  ! it= 0 : GS for generating starting points
  ! it> 1 : frag. MD
  
  ! MOLDEN file
  if(it == 0) then
                      fname='trjM'
  elseif(it > 0.and.it < 9999.and.icoll < 10)then
                  write(fname,'(''MDtrj.'',i4,''.'',i1,''.'',i1)')it,icoll,isec
    if(it < 1000)write(fname,'(''MDtrj.'',i3,''.'',i1,''.'',i1)')it,icoll,isec
    if(it < 100) write(fname,'(''MDtrj.'',i2,''.'',i1,''.'',i1)')it,icoll,isec
    if(it < 10)  write(fname,'(''MDtrj.'',i1,''.'',i1,''.'',i1)')it,icoll,isec
  elseif(it > 0.and.it < 9999.and.icoll >= 10)then
                  write(fname,'(''MDtrj.'',i4,''.'',i2,''.'',i1)')it,icoll,isec
    if(it < 1000)write(fname,'(''MDtrj.'',i3,''.'',i2,''.'',i1)')it,icoll,isec
    if(it < 100) write(fname,'(''MDtrj.'',i2,''.'',i2,''.'',i1)')it,icoll,isec
    if(it < 10)  write(fname,'(''MDtrj.'',i1,''.'',i2,''.'',i1)')it,icoll,isec
  endif
  if(it >= 0.and.it < 9999)  open (file = fname, newunit = io_OUT)
  
  ! this file is used to get the starting points = snapshots
  if(it == 0)then
     open (file = 'qcxms.gs', newunit = io_GS)
     write(io_GS,*)nmax
  endif   
  
  ! ini Ekin
  call ekinet(nuc,velo,mass,Ekin,T)
  
  if(it == 9999)Tsoll=T
  !if(method == 3.or.method == 4)Tsoll=Tsoll
  if(method == 3)Tsoll=Tsoll
  Ekinstart=Ekin
  
  !if (method == 3 .and.it > 0.and.it < 9999)then
  !  Ekin=esub
  !endif !method /= 3
  
  ! in frag runs the electronic temp is set 
  if( it > 0 .and. it < 9999 .and. etempin < 0)then
     call setetemp(nfrag,eimp,etemp)
  else
     etemp=etempin
  endif
  if (No_etemp) etemp=etempin

  ! CID mopac-pm6 behaves better without huge Electronic Temp.
  !if(prog == 1.and.method == 3.and.temprun == .false.)then
  !   write(*,*) 'NOTE:In MOPAC CID runs the electronic T. is set to constant 300 K' 
  !   etemp=300.0d0
  !endif
  
  ! ini Epot      
  spin=0
  call egrad(.true.,nuc,xyz,iat,mchrg,spin,etemp,Epot,grad,achrg,aspin,ECP,gradfail)
  if(gradfail)write(*,*)'Gradient fails in MD 1' !Prob wont happen. Otherwise go into CID and save frags
  if(Epot == 0) return
  
  ! do more more steps in a fragmentation run when nfragexit frag. 
  ! are already there (i.e. after forming 2 frags, a third often occurs fast)
  more    =250
  avdump  =50
  screendump=100
  Tav     =0
  Epav    =0
  Ekav    =0
  Edum    =0
  Eerror  =0
  nstep   =0
  avchrg  =0
  avspin  =0
  aTlast  =0
  morestep=0
  fconst  =0
  dtime   =0
  new_temp=0
  summass =0

  check_fragmented = 1
  cnt = 0
  avxyz2 = 0
  
  ndump=dumpstep
  kdump=avdump     
  vdump=avdump   

  err1 = .false.
  err2 = .false.
  
  ! Stuff for velocity
  call cofmass(nuc,mass,xyz,old_cm)
  do i = 1,nuc
     summass = summass +  mass(i) 
  enddo
  E_kin = 0.5 * summass * ((new_velo*mstoau)**2) 
  E_kin_diff = Ekin - E_kin
  new_temp = (2*E_kin_diff) / (3 *kB* nuc)
  if(method == 3.and.icoll > 0)Ekin = E_kin_diff
  !if(method == 4.and.icoll > 0)Ekin = E_kin_diff
  
  ! no additional calcs in ftemp (trial) runs      
  if(it == 9999)morestep=more+1
  
  ! add the e-hit energy in the first MD steps linearly     
  if(it > 0.and.it < 9999)then
     write(*,'(''Eimp (eV) = '',F6.1,5x,''tauIC (fs) = '',F6.0,6x,''nstep = '',i7,/)')eimp*autoev,tadd/fstoau,nmax

     nadd=(tadd+tstep)/tstep-1   
     fadd=tstep/(tadd+tstep)
     

     write(*,'(''avcycle = '',i4,''    more = '',i4,/)')avdump,more
  else   
     nadd =0
     velof=1.0d0
     screendump=500
  endif
  mdump=screendump
  
  
  write(*,'(''step   time [fs]'',4x,''Epot'',7x,''Ekin'',7x,''Etot'',4x,''error'',2x,''#F   eTemp   frag. T'')')
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MD loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do k=1,nmax
  
     T=Ekin/(0.5*3*nuc*kB)
     ! in frag runs the electronic temp is set 
     if(it <= 0.or.it == 9999)then
        etemp=etempin
        !in mopac cid the elc temp should be 300 ! True? And Negative ion?=
           if(prog == 1.and.method == 3)then
              etemp=300.0d0
           endif
     endif
  
     Tav =Tav+T
     Epav=Epav+Epot
     Ekav=Ekav+Ekin
     !> if no. of steps exceed ~2x tadd (which is nadd roughly)
     !> the average energy is divided by nstep - nadd 
     if(nstep > nadd)then
        Edum=Edum+Epot+Ekin
        Eav =Edum/float(nstep-nadd)
     else
        Eav=Epot+Ekin
     endif 
     Eerror=Eav-Epot-Ekin
  
     ! Avg. Temp.
     if(nfrag == 1)fragT(1)=Tav/k
  
     ! check av. Etot
  !   Eerror=Eav-Epot-Ekin
     if(it == 9999.or.it < 0) Eerror = 0
     if(starting_md ) Eerror = 0

     ! an error ocurred (normally failed to achieve SCF)      
     if (Epot == 0)  err1 = .true. 
     if(method /= 3) then !.or.method /= 4)then
       if (abs(Eerror) > 0.1.and.it /= 9999) err2 = .true.
     else
       if (abs(Eerror) > 0.2.and.it /= 9999) err2 = .true.
     endif
  
     if(err1.or.err2) then
     ! in case of errors take the traj anyway if fragments have been produced              
        if( (nfrag > 1.and.nfrag <= 4).or. isec > 1 ) then
           write(*,8000)nstep,nstep*tstep/fstoau,Epot,Ekin,Epot+Ekin,Eerror,nfrag,etemp,fragT(1:nfrag)
           write(*,*)'EXIT DUE TO SCF CONVERGENCE PROBLEM'
           write(*,*)'OR LARGE MD INTEGRATION ERROR.'
           write(*,*)'because already nfrag > 1 result is taken'
           write(*,*)'(error often occurs for large inter-fragment distances)'
           write(*,*)' and is therefore not very meaningfull)'
           mdok=.true.
        else   
           mdok=.false.
           write(*,*) 'E,error ',Epot,Eerror
        endif   
        !goto 1000
        exit
     endif      
  
     ! average over last avdump steps              
     if(kdump > avdump-1)then
        kdump =0
        avspin=0
        avchrg=0
        avxyz =0
        aTlast=0
     endif

     avchrg = avchrg + achrg
     avspin = avspin + aspin
     avxyz  = avxyz  + xyz 
  
     if(method == 3.and. icoll > 0)   aTlast = aTlast + new_temp    
     if(method /= 3 .or. icoll == 0)  aTlast = aTlast + T    
  
     ! print out every screendump steps
     if(mdump > screendump-1)then
  !      if(method == 3)then
  !        e_int = (t*(0.5*3*nuc*0.316681534524639E-05))/0.03674932217565499
  !        write(*,8000)nstep,ttime,Epot,Ekin,Epot+Ekin,Eerror,nfrag,etemp,fragT(1:nfrag),e_int
  !      else
          write(*,8000)nstep,ttime,Epot,Ekin,Epot+Ekin,Eerror,nfrag,etemp,fragT(1:nfrag)
  !      endif
        mdump=0
     endif
  
     ! dump for MOLDEN
     if(ndump > dumpstep-1.and.it >= 0.and.it < 9999)then
        ndump=0
        write(io_OUT,*)nuc
        write(io_OUT,*)Epot
        do i=1,nuc  
        write(io_OUT,'(1x,a2,1x,3F14.6)')toSymbol(iat(i)),(xyz(j,i)/aatoau,j=1,3)
        enddo
     endif
  
     ! no fragment run, dump to qcxms.gs
     if(it == 0) then
        totdump=totdump+1
        do i=1,nuc  
           write(io_GS,'(6d16.8)')(xyz(j,i),j=1,3),(velo(j,i),j=1,3) 
        enddo
     endif
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! do MD step
     call leapfrog(nuc,grad,mass,tstep,xyz,velo,Ekin,nstep)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! total time including secondaries      
     ttime=ttime+tstep/fstoau
  ! calc E and forces
     call egrad(.false.,nuc,xyz,iat,mchrg,spin,etemp,Epot,grad,achrg,aspin,ECP,gradfail) 
     if(gradfail)write(*,*)'Gradient fails in MD' !Prob wont happen. Otherwise go into CID and save frags
  
     ndump=ndump+1
     mdump=mdump+1
     kdump=kdump+1
     vdump=vdump+1
  
     ! rescale to get Tsoll (NVT ensemble) for equilibration
     dum=100.*abs(Tav/k-Tsoll)/Tsoll

     ! the GS sampling is done in NVE      
     if(dum > 5.0d0.and.it < 0.and.k > 50)then
        velo = velo/sqrt(Tav/k/Tsoll)
     endif
     
     ! Etemp for mopac is not performing well! Set it to fixed 300
     if(prog == 1.and.method == 3)then
        etemp=5000.0d0
     endif

     ! if it oscillates between bonded and not, reset more
     if(nfrag == 1) morestep=0
  
     if(nfrag > 1.and.dtime < 1.d-6) dtime=ttime/1000.
  
     ! check for fragmentation, EXIT section for frag runs
ifit:if(it > 0)then

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Berendsen Thermostat ! Only if not fragmented
        if(method == 3.and.k > 10.and.starting_md.and. nfrag == 1)then
            sca = dsqrt(1.0_wp + ((tstep/fstoau) / 200)*(Tsoll/T-1.0d0))
            velo= sca * (velo) 
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
        ! add the IEE but only if not already fragmented
        if(method /= 3 ) then 
          if(nstep <= nadd.and.nfrag == 1)then
             call impactscale(nuc,velo,mass,velof,eimp,fadd*nstep,Ekinstart)
          endif        
  
          ! reduce eTemp when system is heated up (but only for nfrag=1)        
          dum=eimp-eimp*float(nstep)/float(nadd)
          call setetemp(nfrag,dum,etemp)
        endif        


        !> check if structure is fragmented
        call fragment_structure(nuc,iat,xyz,3.0d0,1,0,list)
        call fragmass(nuc,iat,list,mass,imass,nfrag,fragm,fragf,fragat)
        if(nfrag > 1) then
          call intenergy(nuc,list,mass,velo,nfrag,fragT,E_int)
        endif

  
        ! probably an error
        if(nfrag > 6)then
          write(*,*) 'More than 6 fragments. Error?'
          !goto 1000
          exit
        endif
  
        !> it is important to get the RIGHT kinetic energy inside the molecule
        !> therefor the kinetic energy has to be subtracted from the velocity
        !if (method == 3.and.icoll >= 1.or.method == 4.and.icoll >= 1 &
        !  & .and. .not. Temprun)then 
        if (method == 3.and.icoll >= 1 .and. .not. Temprun)then 

          !> get the length of COM difference
          call cofmass(nuc,mass,xyz,cm)
          diff_cm(:) = cm(:) - old_cm(:)
          cm_out = sqrt(diff_cm(1)**2 + diff_cm(2)**2 + diff_cm(3)**2)

          !> calculate the velocity for this length in steps
          new_velo = (cm_out / tstep) /mstoau !in ms
          old_cm(:) = cm(:)
 
          !> calculate the kinetic energy diff, for new velo
          E_kin = 0.5_wp * summass * ((new_velo*mstoau)**2) 

          !> check the kin. E difference
          E_kin_diff = Ekin - E_kin
          new_temp = (2.0_wp * E_kin_diff) / (3.0_wp * kB * nuc)
          ! write(*,*) new_temp
          vdump = 0

          !> Now set the new Ekin to a value that does not include the difference
          !> by the velocity of the ion 
          Ekin = E_kin_diff
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Check out the fragments. If > 2, do 1000 steps. If more, do 250 steps. 
        ! If nfragexit (can be user-set), exit immidiately
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (nfrag > check_fragmented) then
          cnt = cnt + 1
          !avxyz2  = avxyz2  + xyz0(:,:nuc)
          if (cnt == 50) then
            avxyz2  = avxyz / kdump !cnt
            check_fragmented = nfrag
            !avg_struc = .true.
            cnt = 0
            !avxyz  = 0
          endif
        endif

        ! exit always immidiately if we have > nfragexit frags
        if(nfrag > nfragexit) then
           write(*,8000)nstep,ttime,Epot,Ekin,Epot+Ekin,Eerror,nfrag,etemp,fragT(1:nfrag)
           write(*,9001)
           fragstate=1
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
           write(*,8000)nstep,nstep*tstep/fstoau,Epot,Ekin,Epot+Ekin,Eerror,nfrag,etemp,fragT(1:nfrag)
           write(*,9002)
           fragstate=2       
           exit
        endif   

        ! add a few more cycles because fragmentation can directly proceed further and we don't want to miss this         
        if(nfrag >= nfragexit) then
           morestep=morestep+1
           if(morestep > more) then
              write(*,8000)nstep,ttime,Epot,Ekin,Epot+Ekin,Eerror,nfrag,etemp,fragT(1:nfrag)
              write(*,9003)
              fragstate=1
              exit
           endif 
        endif        
      endif ifit
  enddo ! END LOOP
  
  ! all done and nice
  mdok=.true.
  if(k >= nmax)then
     write(*,9004)
     fragstate=1
  endif

! error exit (Skip if all okay)
!1000  if( it >= 0.and.it < 9999.and.(.not.restart) )close(io_OUT)
  if( it >= 0.and.it < 9999.and.(.not.restart) )close(io_OUT)
  if( it == 0) close(io_GS)

! printed or used in main      
  Tav    = Tav    / k
  Epav   = Epav   / k
  Ekav   = Ekav   / k
  aspin  = avspin / kdump
  achrg  = avchrg / kdump
  !axyz   = avxyz  / kdump
  if (check_fragmented > 1 ) then 
    axyz (1:3,1:nuc) = avxyz2  (1:3,1:nuc) 
  else
   axyz (1:3,1:nuc) = avxyz  (1:3,1:nuc) / kdump
  endif
  aTlast = aTlast / kdump

8000  format(i7,f8.0,F13.5,F9.4,F12.5,F8.4,1x,I1,F9.0,2x,5F7.0,5F3.3) 
9001  format(/8x,'E X I T   M D  because of multiple fragments') 
9002  format(/8x,'E X I T   M D  because nothing happens here anymore') 
9003  format(/8x,'E X I T   M D  because of nfrag=nfragexit') 
9004  format(/8x,'E X I T   M D  because of tmax has been reached') 

      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the internal energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine intenergy(nuc,list,mass,velo,nfrag,T,E_int)
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   use xtb_mctc_constants
   implicit none

   integer  :: i,j,n(10)
   integer  :: nuc,list(nuc),nfrag

   real(wp) :: mass(nuc),velo(3,nuc),T(nfrag)
   real(wp) :: dum
   real(wp) :: E_int(10)

   n=0
   E_int = 0

   do i = 1, nuc
      j = list(i)
      E_int(j) = E_int(j) + 0.5 * mass(i) * (velo(1,i)**2 + velo(2,i)**2 + velo(3,i)**2)
      n(j) = n(j) + 1
   enddo

   dum=sum(E_int(1:nfrag))

   do i=1,nfrag
      T(i) = E_int(i) / (0.5 * 3 * n(i) * kB)
!        write(*,'(''fragment '',i3,'' Ekin = '',F8.4,'' in %'',
!   -               F6.1,5x,''T='',F8.0)')i,e(i),100.*e(i)/dum,T(i)
   enddo

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate MD step 
! Verlet leapfrog algorithm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine leapfrog(nat,grad,amass,tstp,xyz,vel,ke,nstp)
   use common1
   use xtb_mctc_accuracy, only: wp
   implicit none

   integer  :: nat,nstp
   integer  :: j,k

   real(wp) :: grad(3,nat), amass(nat)
   real(wp) :: xyz (3,nat), vel(3,nat)
   real(wp) :: tstp,ke
   real(wp) :: velold, velavg,mass

   ke = 0.0_wp
   do k = 1,nat
      mass = amass(k)
      do j=1,3
         velold   = vel(j,k)
         vel(j,k) = vel(j,k) - (tstp * grad(j,k) / mass)
         velavg   = 0.5_wp * (velold + vel(j,k))
         xyz(j,k) = xyz(j,k) + (tstp * vel(j,k))
         ke       = ke + 0.5_wp * (mass * velavg * velavg)
      enddo
   enddo

   nstp = nstp + 1

end subroutine leapfrog
