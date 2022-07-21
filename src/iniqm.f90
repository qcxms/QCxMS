module qcxms_iniqm
  use common1
  use dftd4, only: gam, gam3
  use newcommon
  use qcxms_tblite, only: get_xtb_egrad, gfn1_xtb, gfn2_xtb, ipea1_xtb
  use qcxms_use_mndo, only: mndoout, mndograd 
  use qcxms_use_msindo, only: getmsindograd 
  use qcxms_use_orca, only: orcaout, rdorcagrad, rdorcaen 
  use qcxms_use_turbomole, only: initm, rdtmgrad, wrcoord, setfermi, rdtmenergy
  use qcxms_utility, only: getspin, qccall
  use xtb_mctc_accuracy, only: wp
  implicit none


  contains

    function get_core_e(iat) result (ncore)
       integer :: iat, ncore
    
       if(iat <= 2)then
          ncore=0
       elseif(iat <= 10)then
          ncore=2
       elseif(iat <= 18)then
          ncore=10
       elseif(iat <= 29)then   !zn
          ncore=18
       elseif(iat <= 36)then
          ncore=28
       elseif(iat <= 47)then
          ncore=36
       elseif(iat <= 54)then
          ncore=46
       elseif(iat <= 71)then
          ncore=54
       elseif(iat <= 79)then
          ncore=68
       elseif(iat <= 86)then
          ncore=78
       endif
    end function get_core_e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine electrons_amount(nat, iat, chrg, nel, nb, z)

       integer, intent(in) :: nat, iat(nat), chrg
       integer, intent(out) :: nel, nb
       integer :: i
       real(wp)  :: z(nat)
      

    ! Set z as valence electrons
       do i=1,nat
          z(i) = iat(i) - get_core_e( iat(i) )
       enddo
    
       nel = idint(sum(z))
       nel=nel-chrg
    
       nb = nel/2          !nb= half the electrons, for uneven e-

    end subroutine electrons_amount

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine iniqm(nat,xyz,iat,chrg,spin,etemp,edum,ok,ECP)

       !!! GCP !!!
    
       integer  :: nat,iat(nat),chrg,spin
       integer  :: idum,jdum
       integer  :: stat
    
       real(wp) :: xyz(3,nat),etemp,edum
       real(wp) :: grad(3,nat),dum(nat),dum2(nat),ec,disp
       real(wp) :: d3a1,d3a2,d3s8,d3atm
       real(wp) :: cn(nat),dcn(3,nat,nat)
    
       logical :: ok
       logical :: ECP
    
       character(len=20) :: atmp
    
       stat=0
       edum=0
    
    ! save info
       jdum=spin
       call getspin(nat,iat,chrg,idum)
    ! sym breaking ini of S=0 but not for dftb+ and only in a
    ! frag run (ok=.false.)
       if(.not.ok.and.idum == 1.and.chrg == 1.and.prog > 0)then
          spin=3
       endif
    
       ec=0
       disp=0
       grad =0
    
       write(*,*)

       ! DFTB+
       if(prog == 0)then
          write(*,*)'initializing dftb+ ...'
          call execute_command_line('rm -f charges.bin')
          call dftbout(nat,xyz,iat,chrg,spin,.false.,etemp,1.d-5)
          call qccall(0,'job.last')
          call dftbenergy(edum)
         ! one more attempt with different diis settings and etemp>0
          if(abs(edum) < 1.d-10)then
             call dftbout(nat,xyz,iat,chrg,spin,.false.,5000.d0,1.d-3)
             call qccall(0,'job.last')
             call dftbenergy(edum)
          endif
          call execute_command_line('rm -f band.out dftb_*in.hsd tmp-broyden*')
       endif

       ! MOPAC
       if(prog == 1)then
          write(*,*)'initializing MOPAC ...'
          call getmopgrad(nat,iat,xyz,grad,chrg,etemp,.true.,&
            edum,dum,dum2)
       endif

       ! TURBOMOLE
       if(prog == 2)then
          call wrcoord(nat,xyz,iat)
          write(*,*)'initializing TM ...'
          call execute_command_line('rm -f control gradient energy')
          call initm(nat,iat,chrg,spin)
          if(etemp > 0) call setfermi(etemp)
          call execute_command_line('kdg scfdump')
          call qccall(2,'job.last')
          call rdtmenergy('energy',edum)
       endif
    
       ! ORCA
       if(prog == 3)then
          write(*,*)'initializing ORCA ...'
          call execute_command_line('rm -f ORCA.INPUT.gbw')
          call orcaout(nat,xyz,iat,chrg,spin,etemp,.false.,ECP)
          call qccall(3,'job.last')
          call rdorcagrad('job.last',nat,grad,dum,dum2,edum)
          call execute_command_line('rm -f ORCA.INPUT.prop')
       endif

       ! MSINDO
       if(prog == 4)then
          write(*,*)'initializing MSINDO ...'
          call getmsindograd(.true.,nat,iat,chrg,spin,xyz,&
          &etemp,14,grad,edum,dum,dum2)
          call execute_command_line('rm -f fort.*')
       endif
    
       ! MNDO
       if(prog == 5)then
          write(*,*)'initializing MNDO99 ...'
          call execute_command_line('rm -f fort.11 fort.15')
          call mndoout(nat,xyz,iat,chrg,spin,etemp,4,10)
          call qccall(5,'job.last')
          call mndograd('job.last',nat,grad,dum,dum2,edum)
         ! two more attempts with different diis settings and etemp>0
          if(abs(edum) < 1.d-10)then
             call mndoout(nat,xyz,iat,chrg,spin,15000.d0,4,20)
             call qccall(5,'job.last')
             call mndograd('job.last',nat,grad,dum,dum2,edum)
          endif
          if(abs(edum) < 1.d-10)then
             call mndoout(nat,xyz,iat,chrg,spin,20000.d0,4,2)
             call qccall(5,'job.last')
             call mndograd('job.last',nat,grad,dum,dum2,edum)
          endif
          call execute_command_line('rm -f fort.15 fort.11')
       endif
    
    
       ! Call-xTB
       if(prog == 6)then
          write(*,*)'initializing (call) XTB ...'
          ! idum=spin
          call callxtb(nat,xyz,iat,chrg,idum,etemp,edum,grad,dum,dum2)
       endif
    
       ! GFN1-xTB
       if(prog == 7)then
          write(*,*)'initializing GFN1-xTB ...'
          call get_xtb_egrad(iat, xyz, chrg, idum, gfn1_xtb, etemp, &
                & "xtb.last", dum, edum, grad, stat, .false.)
       endif
    
       ! GFN2-xTB
       if(prog == 8)then
          write(*,*)'initializing GFN2-xTB ...'
          call get_xtb_egrad(iat, xyz, chrg, idum, gfn2_xtb, etemp, &
                & "xtb.last", dum, edum, grad, stat, .false.)
       endif
    
    
   
       ! Do counter-poise correction
       if(gcp)then
          write(*,*)'initializing gCP ...'
          if(bas == 1)atmp='dft/sv'
          if(bas == 2)atmp='dft/svx'
          if(bas == 3.or.bas == 9)atmp='dft/sv(p)'
          if(bas == 4.or.bas == 10)atmp='dft/svp'
          if(bas <= 4)call calcgcp(nat,iat,xyz,.false.,atmp,ec,grad)
       endif
    
       if(abs(edum) > 1.d-10)then
          edum = edum + ec
    ! D3 is built into xtb, orca, tm and should not be double counted!
          if(prog /= 7 .and. prog /= 6 .and. prog /= 3 .and. prog /=  2 .and. prog /= 8)then
          call gdisp(nat,iat,xyz,d3a1,d3a2,d3s8,d3atm,disp,grad,cn,dcn,.false.)
    !         call gdispxtb(nat,xyz,iat,.true.,grad,ec)
             edum = edum + disp
          endif
    
    
          ok=stat==0
          if (ok) &
          write(*,*)'initialization successful!'
          write(*,'('' energy = '',F14.5)')edum
          write(*,'('' charge = '',I14  )')chrg
          write(*,'('' mult   = '',I14  )')idum
          write(*,'('' etemp  = '',F14.1)')etemp
          if(idum == 1.and.chrg == 1.and.prog > 0)write(*,*)'initializing S=0 state using HS-UHF.'
       else
          write(*,*) 'QM code initialization failure'
          ok=.false.
       endif
    
       spin=jdum
    
    end subroutine iniqm
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    ! compute energies only as used in fragment IP calc
    subroutine eqm(prog,nat,xyz,iat,mchrg,dumspin, &
        & etemp,newcalc,ipok,energy,nel,nb,ECP,spec_calc)
    
       integer  :: prog,nat,iat(nat),mchrg
       integer  :: spin
       integer  :: dumspin
       integer  :: nb
       integer  :: nel
    
       real(wp) :: xyz(3,nat),energy,etemp
       real(wp) :: gradient(3,nat),qat(nat),dum2(nat)
       real(wp) :: z(nat)
    
       logical  :: newcalc, spec_calc
       logical  :: ECP
       logical  :: ipok
    
       integer :: stat
       character(len=:), allocatable :: output_name
    
       call electrons_amount(nat, iat, mchrg, nel, nb, z)

    ! spec.f and regular IP/EA evaluations will use this input option
    ! where the algorithm assigns the spin automatically
    ! However, because metals need to vary the spin, it needs to be an
    ! input variable as well.
    ! NOTE: IF SPIN IS INPUT VARIABLE IT SHOULD COME IN AS MULTIPLICITY?
   
       if(dumspin == -1)then
          call getspin(nat,iat,mchrg,spin)
       else
          spin = dumspin
       endif
       stat=0
       energy=0
      ! no electrons 
      ! we do not want this for XTB since proton has 0.2 Ha energy in XTB
       if(prog < 6)then
          if(spin < 1) then
             energy=1.d-6
             return
          endif
       endif
    
       gradient=0

       ! DFTB+
       if(prog == 0)then
          call execute_command_line('rm -f charges.bin')
          call dftbout(nat,xyz,iat,mchrg,spin,.false.,etemp,1.d-5)
          if(mchrg == 0)then
             call qccall(0,'neutral.out')
          else
             call qccall(0,'ion.out')
          endif
          call dftbenergy(energy)
       endif

       ! MOPAC
       if(prog == 1)then
          call getmopgrad(nat,iat,xyz,gradient,mchrg,etemp, &
            .true.,energy,qat,dum2)
       endif

       ! TURBOMOLE
       if(prog == 2)then
          call execute_command_line('rm -f control gradient energy')
          call wrcoord(nat,xyz,iat)
    
          if(mchrg == 0)then
             call initm(nat,iat,mchrg,spin)
             if(etemp > 0) call setfermi(etemp)
             call execute_command_line('kdg scfdump')
             call qccall(2,'neutral.out')
             call execute_command_line('cp energy energy.neutral')
             call execute_command_line('cp gradient gradient.neutral')
             call rdtmenergy('energy.neutral',energy)
          else
             call initm(nat,iat,mchrg,spin)
             if(etemp > 0) call setfermi(etemp)
             call execute_command_line('kdg scfdump')
             call qccall(2,'ion.out')
             call execute_command_line('cp energy energy.ion')
             call execute_command_line('cp gradient gradient.ion')
             call rdtmenergy('energy.ion',energy)
          endif
       endif

       ! ORCA
       if(prog == 3)then
          if(newcalc) call execute_command_line('rm -f ORCA.INPUT.gbw')
          call orcaout(nat,xyz,iat,mchrg,spin,etemp,.false.,ECP)
          if(mchrg == 0)then
             call qccall(3,'neutral.out')
             call rdorcaen('neutral.out',energy)
          else
             call qccall(3,'ion.out')
             call rdorcaen('ion.out',energy)
          endif
          if(abs(energy) < 1.d-10)then
             if(mchrg == 0)then
                call qccall(3,'neutral.out')
                call rdorcaen('neutral.out',energy)
             else
                call qccall(3,'ion.out')
                call rdorcaen('ion.out',energy)
             endif
          endif
       endif

       ! MSINDO
       if(prog == 4)then
          call getmsindograd(.true.,nat,iat,mchrg,spin,xyz,etemp,14, &
            gradient,energy,qat,dum2)
          call execute_command_line('rm -f fort.*')
       endif
    
       ! MNDO
       if(prog == 5)then
          call execute_command_line('rm fort.11 fort.15')
          call mndoout(nat,xyz,iat,mchrg,spin,etemp,4,10)
          if(mchrg == 0)then
             call qccall(5,'neutral.out')
             call mndograd('neutral.out',nat,gradient,qat,dum2,energy)
          else
             call qccall(5,'ion.out')
             call mndograd('ion.out',nat,gradient,qat,dum2,energy)
          endif
       endif
    
    
       ! Call-xTB
       if (prog == 6) then
          call callxtb(nat,xyz,iat,mchrg,spin,etemp,energy,gradient,qat,dum2)
          if( mchrg == 0)then
             call execute_command_line('mv xtb.last neutral.out')
          else
             call execute_command_line('mv xtb.last ion.out')
          endif
       end if
    
       ! GFN1-xTB
       if (prog == 7) then
          if (mchrg == 0) then
             output_name = "neutral.out"
          else
             output_name = "ion.out"
          end if
          if(nel == 0)then
            call eself(nat,iat,z,energy)
            qat    = 0.0_wp
          else
            call get_xtb_egrad(iat, xyz, mchrg, spin, ipea1_xtb, etemp, &
                   & output_name, qat, energy, gradient, stat, spec_calc)
          end if
       end if
    
       ! GFN2-xTB
       if (prog == 8) then
          if (mchrg == 0) then
             output_name = "neutral.out"
          else
             output_name = "ion.out"
          end if
          if (nel == 0) then
            call eself(nat,iat,z,energy)
            qat    = 0.0_wp
          else
            call get_xtb_egrad(iat, xyz, mchrg, spin, gfn2_xtb, etemp, &
                &    output_name, qat, energy, gradient, stat, spec_calc)
            if (stat /= 0) then
              error stop "[Fatal] Calculation in tblite library failed"
            end if
          end if
       endif
    
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (abs(energy) < 1.d-10 .or. stat/=0) then
          write(*,*)'QM code failure'
          write(*,*)'prog   = ',prog
          write(*,*)'charge = ',mchrg
          write(*,*)'mult   = ',spin
          ipok = .false.
       endif

    end subroutine eqm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine eself(n,at,z,edum)
    
       integer  :: n,at(n)
       integer  :: i,ii
    
       real(wp) :: z(n)
       real(wp) :: e
       real(wp), intent(out) :: edum
    
       edum=0.0d0
       e=0
       do i=1,n
          ii = at(i)
          e = e + 0.5_wp * z(i)**2 * gam(ii) + z(i)**3 * gam3(ii) / 3.0_wp
       enddo
    
    !   write(*,*)'NO ELECTRONS!'
    !   write(*,'(''molecular/atomic self energy :'',F10.6)') e
       edum = e
    !      call wren(e)
    
    end subroutine eself

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! get energy, gradient, charge and spin density on atoms from QC
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  subroutine egrad(first,nuc,xyz,iat,mchrg,spin,etemp,E,grad,qat,aspin,ECP,gradfail)
     use common1
     use qcxms_tblite, only: get_xtb_egrad, gfn1_xtb, gfn2_xtb
     use qcxms_use_mndo, only: mndoout, mndograd 
     use qcxms_use_msindo, only: getmsindograd 
     use qcxms_use_orca, only: orcaout, rdorcagrad 
     use qcxms_use_turbomole, only: initm, rdtmgrad, setfermi
     use qcxms_utility, only: qccall, getspin
     use xtb_mctc_accuracy, only: wp
     use xtb_mctc_convert
     use xtb_mctc_symbols, only: toSymbol 
     implicit none
  
     integer :: nuc
     integer :: mchrg
     integer :: spin
     integer :: iat(nuc)
     integer :: i,idum,stat
  
     real(wp) :: xyz (3,nuc)
     real(wp) :: grad(3,nuc)
     real(wp) :: qat(nuc) ! atomic charge (q-at)
     real(wp) :: aspin(nuc)
     real(wp) :: E,etemp
     real(wp) :: d3a1,d3a2,d3s8,d3atm
     real(wp) :: cn(nuc),dcn(3,nuc,nuc)
     real(wp) :: egcp,ggcp(3,nuc),disp
  
     character(len=20) :: atmp
  
     logical :: first
     logical :: ECP
     logical :: ok,gradfail
  
     disp   = 0
     E    = 0
     grad = 0
     aspin= 0
     qat= 0
     gradfail=.false.
  
    !ccccccccccccccccccccccccccc
    ! dftb+ (DFTB3-D3(BJ))
    !ccccccccccccccccccccccccccc
  
     if(prog == 0)then
        call dftbout(nuc,xyz,iat,mchrg,spin,.true.,etemp,1.d-5)
        call qccall(0,'job.last')
        call dftbgrad(nuc,grad,qat,aspin,E)
        call checkqc(nuc,E,grad,qat,ok,mchrg)
        if(.not.ok) then
           call dftbout(nuc,xyz,iat,mchrg,spin,.false.,etemp,1.d-4)
           call qccall(0,'job.last2')
           call dftbgrad(nuc,grad,qat,aspin,E)
           call checkqc(nuc,E,grad,qat,ok,mchrg)
           if(ok)write(*,*)'healed!'
        endif
     endif
  
  !ccccccccccccccccccccccccccc
  ! mopac
  !ccccccccccccccccccccccccccc
     if(prog == 1)then
        call getmopgrad(nuc,iat,xyz,grad,mchrg,etemp,.false.,E,qat,aspin)
        call checkqc(nuc,E,grad,qat,ok,mchrg)
     endif
  
  !ccccccccccccccccccccccccccc
  ! TM
  !ccccccccccccccccccccccccccc
  
     if(prog == 2)then
        call execute_command_line('rm -rf gradient energy dscf_problem')
        open(unit=87,file='coord')
        write(87,'(a)')'$coord'
        do i=1,nuc
           write(87,'(3F16.10,10x,a2)')&
           &xyz(1,i),xyz(2,i),xyz(3,i),toSymbol(iat(i))
        enddo
        write(87,'(a)')'$end'
        close(87)
        if(etemp.gt.0) call setfermi(etemp)
        call qccall(2,'job.last')
        call rdtmgrad (nuc,grad,qat,aspin,E)
        call checkqc(nuc,E,grad,qat,ok,mchrg)
     endif
  
  !ccccccccccccccccccccccccccc
  ! ORCA
  !ccccccccccccccccccccccccccc
  
     if(prog == 3)then
        call orcaout(nuc,xyz,iat,mchrg,spin,etemp,.true.,ECP)
        call qccall(3,'job.last')
        call rdorcagrad('job.last',nuc,grad,qat,aspin,E)
  ! second attempt
        call checkqc(nuc,E,grad,qat,ok,mchrg)
        if(.not.ok) then
           call qccall(3,'job.last2')
           call rdorcagrad('job.last2',nuc,grad,qat,aspin,E)
           call checkqc(nuc,E,grad,qat,ok,mchrg)
           if(ok)write(*,*)'healed!'
        endif
     endif
  
  !ccccccccccccccccccccccccccc
  ! MSINDO
  !ccccccccccccccccccccccccccc
  
     if(prog == 4)then
        call getmsindograd(first,nuc,iat,mchrg,spin,xyz,etemp,4,grad,E,qat,aspin)
 
        !> for msindo, this was originally checkqc2, but no warning was triggered
        call checkqc(nuc,E,grad,qat,ok,mchrg)
        if(.not.ok) call getmsindograd(first,nuc,iat,mchrg,spin,xyz,etemp,14,grad,E,qat,aspin)
  
        call checkqc(nuc,E,grad,qat,ok,mchrg)
        if(.not.ok) call getmsindograd(first,nuc,iat,mchrg,spin,xyz,etemp,16,grad,E,qat,aspin)
  
        call checkqc(nuc,E,grad,qat,ok,mchrg)
     endif
  
  !ccccccccccccccccccccccccccc
  ! MNDO99
  !ccccccccccccccccccccccccccc
  
     if(prog == 5)then
        call mndoout(nuc,xyz,iat,mchrg,spin,etemp,5,10)
        call qccall(5,'job.last')
        call mndograd('job.last',nuc,grad,qat,aspin,E)
        call checkqc(nuc,E,grad,qat,ok,mchrg)
        if(.not.ok) then
           call mndoout(nuc,xyz,iat,mchrg,spin,etemp,5,20)
           call qccall(5,'job.last2')
           call mndograd('job.last2',nuc,grad,qat,aspin,E)
           call checkqc(nuc,E,grad,qat,ok,mchrg)
           if(ok)write(*,*)'healed!'
        endif
        if(.not.ok) then
           call mndoout(nuc,xyz,iat,mchrg,spin,etemp,4,2)
           call qccall(5,'job.last2')
           call mndograd('job.last2',nuc,grad,qat,aspin,E)
           call checkqc(nuc,E,grad,qat,ok,mchrg)
           if(ok)write(*,*)'healed!'
        endif
     endif
  
  
  !cccccccccccccccccccccccccc
  !  XTB
  !ccccccccccccccccccccccccc
     if(prog == 6)then
  ! idum = spin
        call getspin(nuc,iat,mchrg,idum)
        call callxtb(nuc,xyz,iat,mchrg,idum,etemp,E,grad,qat,aspin)
        call checkqc(nuc,E,grad,qat,ok,mchrg)
        if(.not.ok)then
           gradfail=.true.
           write(*,*) 'GRAD failed!'
        endif
     endif
  
  !cccccccccccccccccccccccccc
  ! GFN1-xTB
  !ccccccccccccccccccccccccc
     if(prog == 7)then
        call getspin(nuc,iat,mchrg,idum)
  !      if (with_tblite) then
           call get_xtb_egrad(iat, xyz, mchrg, idum, gfn1_xtb, etemp, &
              & "xtb.last", qat, E, grad, stat, .false.)
           ok = stat == 0
  !      else
  !         call qcxtb2(nuc,xyz,iat,mchrg,idum,E,etemp,grad,qat,aspin,.false.)
  !      end if
        call checkqc(nuc,E,grad,qat,ok,mchrg)
        if(.not.ok)then
           gradfail=.true.
           write(*,*) 'QC calc failed!'
        endif
     endif
  
  !cccccccccccccccccccccccccc
  ! GFN2-xTB
  !ccccccccccccccccccccccccc
     if(prog == 8)then
        call getspin(nuc,iat,mchrg,idum)
  !      if (with_tblite) then
           call get_xtb_egrad(iat, xyz, mchrg, idum, gfn2_xtb, etemp, &
              &               "xtb.last", qat, E, grad, stat, .false.)
           ok = stat == 0
  !      else
  !         call qcxtb2(nuc,xyz,iat,mchrg,idum,E,etemp,grad,qat,aspin,.false.)
  !      end if
        call checkqc(nuc,E,grad,qat,ok,mchrg)
        if(.not.ok)then
           gradfail=.true.
           write(*,*) 'QC calc failed!'
        endif
     endif
  
     if(.not.ok) return
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! gcp and D3
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(gcp) then
        atmp='dft/sv'
        if(bas == 2)atmp='dft/svx'
        if(bas == 3.or.bas == 9)atmp='dft/sv(p)'
        if(bas == 4.or.bas == 10)atmp='dft/svp'
        call calcgcp(nuc,iat,xyz,.false.,atmp,egcp,ggcp)
        grad(1:3,1:nuc)=grad(1:3,1:nuc)+ggcp(1:3,1:nuc)
        E=E+egcp
     endif
  
    ! return with ec=0 if parameters are zero (see input.f)
    ! make d3 standalone only available for dftb(0), msindo(4) or mndo99(5)
     if (prog == 0 .or. prog ==  4 .or. prog ==  5)then
        call gdisp(nuc,iat,xyz,d3a1,d3a2,d3s8,d3atm,disp,grad,cn,dcn,.false.)
        E = E + disp
     endif
  
  end subroutine egrad
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! is E, grad and charge reasonable ?
  subroutine checkqc(nuc,e,grad,qat,ok,mchrg)
  
     integer  :: nuc,mchrg
     real(wp) :: grad(3,nuc),e,qat(nuc)
     real(wp) :: gn
     logical  :: ok
  
  
     ok=.false.
  
     if(abs(E).lt.1.d-8) then
        write(*,*)'energy seems to be bogus!'
        return
     endif
  
     gn=gnorm(nuc,grad)
     if(gn.lt.1.d-8.or.gn.gt.20.0) then
        write(*,*)'grad seems to be bogus!'
        E=0
        return
     endif
  
     if ( mchrg > 0 ) then
        if(abs(maxval(qat)).lt.1.d-5) then
           write(*,*)'WF seems to be bogus!',sum(qat)
           E=0
           return
        endif
     endif
     ok=.true.
  
  end subroutine checkqc
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute norm of electronic gradient for dftd3
  function gnorm(nuc,grad)
     integer  :: nuc,i
     real(wp) :: gnorm 
     real(wp) :: grad(3,nuc)
     real(wp) :: gn
  
     gn=0
     do i=1,nuc
        gn=gn+grad(1,i)**2+grad(2,i)**2+grad(2,i)**2
     enddo
  
     gnorm=sqrt(gn)
  
  end function gnorm

end module qcxms_iniqm
