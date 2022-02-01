module cid_module

subroutine cid_outer(set, mchrg_prod,

  integer :: mchrg, mchrg_prod

  real(wp) :: chrgcont  

  type(run_settings)   :: set
  type(collision_type) :: coll

!mCID:if ( method == 3 ) then 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      set%tadd  = 0.0_wp !should be 0 anyway.
      set%eimp  = 0.0_wp !should be 0 anyway.
      set%velof = 0.0_wp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !if ( method == 3 ) mchrg =  mchrg_prod
      !if ( method == 4 ) mchrg = -1 * mchrg_prod

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
        edum = t * (0.5 * 3 * nuc *kB * autoev)
        ENe(1) = edum

        !> if not set manually, determine automatically and distribute
        if (ESI == 0 .and. tempESI == 0 ) then
          !> set random seed and number 
          !> else the random seed from start is taken
          if (iseed(1) == 0) then
            call random_seed()
          endif

          call random_number(a)

          b    = nuc / 10.0
          dep  = nint(b)
          if ( b < 2.0_wp ) dep = 1 ! hard coded for small molecules

          !>> Make random energy depending on molecular size
          if (dep == 1) rand_int = 1 + FLOOR(2*a)     ! 1-2
          if (dep == 2) rand_int = 2 + FLOOR(2*a)     ! 2-3
          if (dep == 3) rand_int = 3 + FLOOR(3*a)     ! 3-5 !höher für anderen bereich
          if (dep == 4) rand_int = 3 + FLOOR(4*a)     ! 3-6
          if (dep >= 5) rand_int = 4 + FLOOR(5*a)     ! 4-8


          write(*,'('' Randomly Scaling internal energy: '')')
          write(*,'(50(''-''))') 
          write(*,'('' Inner Energy before          : '', f14.6, a3)') ENe(1), ' eV'
          write(*,'('' Temperature  before          : '', f14.6, a3)')      T, ' K'
          write(*,'(50(''-''))') 
          write(*,'('' Wanted value Energy          : '', i4, a3)'  ) rand_int, ' eV'
          write(*,'(50(''-''))') 

          !>> scale depending on random number (or not, if value too low)
          edum = rand_int - edum
          if (edum > 0) then
            E_Scale = vary_energies(edum, 0.1_wp)
            if (rand_int  == 0) E_Scale = 0
            write(*,'('' Scaling to inner Energy      : '', f14.6,a3)') E_Scale+ENe(1),' eV'
          else
            write(*,*) ' ! No Scaling ! '
            E_Scale = 0.0_wp
          endif

        !> if ESI value was set manually
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

            !>> decide if manual value is distributed randomly or not
            if ( NoScale ) then ! No distribution of ESI 
              E_Scale = edum
            else
              E_Scale = vary_energies(edum, 0.2_wp) !distribute width 0.2 (see boxmuller)
            endif
            write(*,'('' Scaling to inner Energy      : '', f14.6,a3)') E_Scale+ENe(1),' eV'
          endif

        !> if temperature is scaled instead of energy
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

        !> stop if both values were provided (to circumvent errors in input)
        if ( tempESI > 0 .and. ESI > 0 ) stop 'Cannot provide both, Energy and Temp.!'

        !> make sure nothing strange is scaled
        if (E_Scale <= 0) E_Scale = 0

        !> convert values for output
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
                prestep = nint(ENe(4)/100)
                prestep = prestep*200
                starting_MD = .true.
              else
                prestep=simMD
                starting_MD = .false.
                if(fragstate == 2) prestep=simMD * 0.75
              endif

            elseif (TempRun .and. isec == 1) then
              prestep = simMD
              starting_MD = .true.
            elseif (TempRun .and. isec > 1) then
              starting_MD = .false.
              if(fragstate == 2) prestep = prestep * 0.75
            endif


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Pre-CID MD Loop => ESI MD
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call md(itrj,icoll,isec,nuc,prestep,xyz,iat,mass,imass,mchrg,     &
            & grad,velo,velof,list,tstep,j,nfragexit,fragm,fragf,             &
            & fragat,dumpstep,etemp_in,md_ok,chrg,spin,axyz,tscale,tadd,      &
            & eimp,.false.,Tav,Epav,Ekav,ttime,aTlast,fragstate,dtime,        &
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
              call manage_fragments(nuc, iat, xyz, axyz, mchrg, chrg, spin, mass, &
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
              coll_counter = coll_counter + 1
              write(*,'(40(''!''))')
              write(*,'(''!! Fragmentation in ESI MD !!'')')
              write(*,'(''-- No of overall fragmentations: '',i3, '' --'')')&
                coll_counter
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
            if ( nuc <= 6 ) nmax = nmax0 / 4 !CID may fragment more

            ! do not continue with small fragments
            if ( nuc <= 5 ) then
              small = .true.
              exit
            endif

            ! do not continue with low masses/resolution of instrument (user)
            !if ( sum(mass(1:nuc)) / amutoau <=  minmass ) then
            !  littlemass = .true.
            !  exit
            !endif

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

      !Give info for not ESI runs and then no Coll. either
      if (.not. TempRun .and. small ) write(*,*) &
        ' Simulation stopped - too small molecule for collisions '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Initial caclulation of number of collisions 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cnt:  if (.not. TempRun .and. .not. small .and. isec < 3) then

        starting_md = .false.
        save_simMD = simMD
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
          simMD     = save_simMD


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
          & mfpath, rtot, chrg, icoll, collisions, direc, new_velo,      &
          & aTlast, calc_collisions, imass)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          new_temp = aTlast

          ! Grad failure (too big steps)
          if ( stopcid ) write(*,*) 'Error occured in the CID module'

          ! END STUFF HERE
          if ( stopcid ) exit

          !> do IP calc and write out fragment files
          call manage_fragments(nuc, iat, xyz, axyz, mchrg, chrg, spin, mass, &
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
            coll_counter = coll_counter + 1
            write(*,'(''-- No of overall fragmentations: '',i3, '' --'')')&
              coll_counter

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
            call cofmass(nuc,mass,xyz,cm)
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

            !> store the old value to reset later
            save_simMD = simMD

            !> change MFP times to reduce timings (empirical values)
            !> but only if not set manually
            if ( manual_simMD == 0 ) then
              simMD = icoll * 0.6 * nuc * 100
              if ( simMD > 8000 .and. coll_counter <= 2 ) simMD = 8000
              
              !>> make some timing adjustments
              if ( simMD > 8000 .and. coll_counter > 2  ) then
                simMD = 8000 * 0.75_wp
              elseif ( simMD > 8000 .and. coll_counter > 3  ) then
                simMD = 8000 * 0.6_wp
              elseif ( simMD > 8000 .and. coll_counter >= 4  ) then
                simMD = 8000 * 0.5_wp
              endif

              !>> not too short simulations
              if ( simMD < 2000 ) simMD = 2000
            endif

            !> reduce the MD time if fragmentation in MFP occurs
            !> even if manually set
            if (isec == 3) simMD =int(simMD * 0.75_wp)
            if (isec == 4) simMD =int(simMD * 0.6_wp )
            if (isec >= 5) simMD =int(simMD * 0.5_wp )
            !if ( fragstate == 2 ) simMD = simMD / 2


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Do Mean-free-path (MFP) MD with simMD timesteps
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call md(itrj,icoll,isec,nuc,simMD,xyz,iat,mass,imass,mchrg,grad,&
            &       velo,velof,list,tstep,j,nfragexit,                      &
            &       fragm,fragf,fragat,dumpstep,etemp_in,                    &
            &       md_ok,chrg,spin,axyz,                                    &
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

            ! reset sim MD time
            simMD = save_simMD

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
            if(md_ok)then
              !> do IP calc and write out fragment files
              call manage_fragments(nuc, iat, xyz, axyz, mchrg, chrg, spin, mass, &
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
              coll_counter = coll_counter + 1
              write(*,'(''-- No of overall fragmentations: '',i3, '' --'')')&
                coll_counter

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
          endif Coll

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

end subroutine cid_outer
end module cid_module
