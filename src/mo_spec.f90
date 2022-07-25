module qcxms_mo_spec
   use common1
   use newcommon
   use qcxms_eps, only: geteps
   use qcxms_mo_energy, only: readpopmo
   use qcxms_iniqm, only: eqm
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   implicit none

   contains
     subroutine getspec(echo,nat,iat,xyz,mchrg,emo,ehomo,mopop,ihomo,nb,ECP)

     integer  :: nat,iat(nat),mchrg,ihomo
     integer  :: i,j,k,ncore
     integer  :: nel,nb,nao

     real(wp) :: eps (2,5*nat)
     real(wp) :: emo (5*nat)
     real(wp) :: mopop(5*nat,nat)
     real(wp) :: xyz(3,nat)
     real(wp) :: pmo(5*nat),energy,dum
     real(wp) :: focc(5*nat)
     real(wp) :: ehomo

     logical  :: echo
     logical  :: ECP
     logical  :: ex
     logical  :: spec_calc = .true.
     logical  :: ipok


     mopop = 0.0d0

     ipok = .true.

     ! save original level
      j=func
      k=bas

     ! for PBE12/SVx:
     bas  = 2
     func = 2
     if(ecp)then
        bas = 9
     endif
     inquire(file='qcxms.Mspec',exist=ex)
     if(.not.ex)then
        if(.not. XTBMO)then

          ! USE ORCA
           write(*,*) 'Not yet working'
           ! eps has to be rewritten; changed it in readpopmo to emo and focc
           error stop
           write(*,*) 'doing DFT for MO spectrum for M ...'
           call eqm(3,nat,xyz,iat,mchrg,-1,0.0_wp,.false.,ipok,energy,nel,nb,ECP,spec_calc)
           call execute_command_line('rm -f ORCA.INPUT.prop')
           call execute_command_line('mv neutral.out qcxms.Mspec')
           if(abs(energy).lt.1.d-8.or. .not. ipok ) error stop 'QC initialization error'
           call geteps('qcxms.Mspec',nat,iat,mchrg,eps,mopop,ncore,ihomo,nb)
           emo (1:ihomo) = eps(1,1:ihomo)

        ! USE XTB
        elseif(prog.eq.7)then
           write(*,*) 'doing XTB for MO spectrum for M ...'
           call eqm(7,nat,xyz,iat,mchrg,-1,0.0_wp,.false.,ipok,energy,nel,nb,ECP,spec_calc)
           call execute_command_line('mv neutral.out qcxms.Mspec.xtb')
           if(abs(energy).lt.1.d-8.or. .not. ipok) error stop 'QC initialization error'
           call readpopmo(nat,nao,ihomo,emo,focc,mopop)
           
        ! USE XTB2 (default for everything else)
        else
           write(*,*) 'doing XTB2 for MO spectrum for M ...'
           call eqm(8,nat,xyz,iat,mchrg,-1,0.0_wp,.false.,ipok,energy,nel,nb,ECP,spec_calc)
           call execute_command_line('mv neutral.out qcxms.Mspec.xtb2')
           if(abs(energy).lt.1.d-8.or. .not. ipok) error stop 'QC initialization error'
           call readpopmo(nat,nao,ihomo,emo,focc,mopop)
        endif

     else
        ! For xtb this will be ignored, it is ORCA formatted 
        call geteps('qcxms.Mspec',nat,iat,mchrg,eps,mopop,ncore,ihomo,nb)
        emo (1:ihomo) = eps(1,1:ihomo)
     endif

      ! RESTORE level of theory to original for future calculations
      func=j
      bas =k

!       if(ihomo.ne.nb) stop 'ihomo<>nb'

     if(echo)then
        if (.not. XTBMO) then
           write(*,*) 'PBE12/SVx levels (eV):'
        else
           write(*,*) 'XTB-SCC levels (eV):'
        endif
!        write(*,'('' M alpha/beta '',10F8.2)') eps(1,1:ihomo)
        write(*,'('' M alpha/beta '',10F8.2)') emo(1:ihomo)
        do j=1,ihomo
           dum=0
           do i=1,nat
              dum=dum+mopop(j,i)**2
           enddo
           pmo(j) = 1.0_wp / (dum+1.0e-8_wp)
        enddo
        write(*,*) 'MO localization degree (number of centers)'
        write(*,'('' alpha/beta '',10F8.2)')pmo(1:ihomo)
     endif

     if (emo(ihomo)-emo(1).lt.1.0) error stop 'weird epsilons'

     emo =(-1.0_wp) * emo * evtoau  ! convert to neg. and au

     ehomo = emo(ihomo) ! get the HOMO orbital energy

     !return

     end subroutine getspec

end module qcxms_mo_spec
