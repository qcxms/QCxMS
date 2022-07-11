! This file is part of qcxms.

!> This module defines the bridge between QCxMS and the tblite library.
module qcxms_tblite
   use mctc_env, only : error_type
   use mctc_io, only : structure_type, new
   use tblite_context_type, only : context_type
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_integral_overlap !, only : get_overlap_lat
   use tblite_cutoff !, only : get_lattice_points
   use tblite_basis_type, only : get_cutoff
   use qcxms_mo_energy, only: write_qmo 
   use xtb_mctc_convert
   implicit none
   private

   public :: get_xtb_egrad
   public :: gfn1_xtb, gfn2_xtb, ipea1_xtb

   !> Required precision for the library
   integer, parameter :: wp = selected_real_kind(15)

   !> Type to implement an enumerated method selector
   type :: method_selector
      integer :: id
   end type method_selector

   !> Selector for GFN1-xTB Hamiltonian
   type(method_selector), parameter :: gfn1_xtb = method_selector(1)

   !> Selector for GFN2-xTB Hamiltonian
   type(method_selector), parameter :: gfn2_xtb = method_selector(2)

   !> Selector for IPEA1-xTB Hamiltonian
   type(method_selector), parameter :: ipea1_xtb = method_selector(11)

   !> Conversion factor from Kelvin to Hartree
   real(wp), parameter :: ktoau = 3.166808578545117e-06_wp

   !> Calculation accuracy
   real(wp), parameter :: accuracy = 1.0_wp

   !> Label for sending informative output back
   character(len=*), parameter :: info_label = "[Info] "

   !> Label for sending errors back
   character(len=*), parameter :: error_label = "[Fatal] "

   !> Return status in case of unknown error
   integer, parameter :: stat_fatal = -1

   !> Return status in case an unknown method has been requested
   integer, parameter :: stat_unknown_method = 5


contains


!> Entry point for QCxMS to request calculations from the tblite library
subroutine get_xtb_egrad(num, xyz, charge, multiplicity, method, etemp, &
      & output_file, qat, energy, gradient, stat, spec_calc)
   !> Atomic numbers for each atom
   integer, intent(in) :: num(:)
   !> Cartesian coordinates for each atom in Bohr
   real(wp), intent(in) :: xyz(:, :)
   !> Total atomic charge of the system
   integer, intent(in) :: charge
   !> Total multiplicity of the system
   integer, intent(in) :: multiplicity
   !> Selected method
   type(method_selector), intent(in) :: method
   !> Electronic temperature
   real(wp), intent(in) :: etemp
   !> Name of the output file to write to
   character(len=*), intent(in) :: output_file
   !> Atomic partial charges
   real(wp), intent(out) :: qat(:)
   !> Total energy
   real(wp), intent(out) :: energy
   !> Molecular gradient
   real(wp), intent(out) :: gradient(:, :)
   !> Error status
   integer, intent(out) :: stat
   !> HOMO
   integer  :: ihomo
   real(wp) :: ehomo
   !> Calculate MO ?
   logical  :: spec_calc

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(error_type), allocatable :: error
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: sigma(3, 3)

   real(wp), allocatable :: overlap(:, :), trans(:, :)
   real(wp) ::  cutoff

   stat = 0

   ! Associate output file with calculation context, all printout will be redirected
   open(newunit=ctx%unit, file=output_file)

   ! Create a new molecular structure data container from the provided input
   call new(mol, num, xyz, charge=real(charge, wp), uhf=min(multiplicity-1, 0))

   ! Expand the method selector to an actual calculator
   select case(method%id)
   case default
      call ctx%message(error_label//"Unknown method requested for calculation in tblite")
      stat = stat_unknown_method

      close(ctx%unit)
      return
   case(gfn2_xtb%id)
      call ctx%message(info_label//"Starting GFN2-xTB calculation in tblite")
      call new_gfn2_calculator(calc, mol)
   case(gfn1_xtb%id)
      call ctx%message(info_label//"Starting GFN1-xTB calculation in tblite")
      call new_gfn1_calculator(calc, mol)
   case(ipea1_xtb%id)
      call ctx%message(info_label//"Starting IPEA1-xTB calculation in tblite")
      call new_ipea1_calculator(calc, mol)
   end select

   ! Create a new wavefunction for every calculation
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, etemp * ktoau)

   ! Perform the actual calculation
   call xtb_singlepoint(ctx, mol, calc, wfn, accuracy, energy, gradient, sigma, 1)


   ! Check the calculation context for errors
   if (ctx%failed()) then
      ! Tear down the error stack to send the actual error messages back
      call ctx%message(error_label//"Singlepoint calculation failed")
      do while(ctx%failed())
         call ctx%get_error(error)
         call ctx%message("-> "//error%message)
      end do
      stat = stat_fatal
   end if

!   ! Copy results to output variables
   qat(:) = wfn%qat(:, 1)

   if (spec_calc)then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(overlap(calc%bas%nao, calc%bas%nao))
      cutoff = get_cutoff(calc%bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)
      call get_overlap(mol, trans , cutoff, calc%bas, overlap)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ihomo  = max(wfn%homo(1), 1)
      ehomo  = wfn%emo(ihomo, 1)

      call write_qmo(mol%nat, calc%bas%nao, calc%bas%ao2at, wfn%coeff(:, :, 1), &
         & overlap, wfn%emo(:, 1), wfn%focc(:, 1), ihomo)
   endif
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   ! Cleanup, close all opened files
   close(ctx%unit)

end subroutine get_xtb_egrad

end module qcxms_tblite
