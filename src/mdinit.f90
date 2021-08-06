module qcxms_mdinit
   use common1
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_constants
   implicit none 
  
   contains

! initialize velocities uniformly
   subroutine mdinitu(nat,velo,mass,E_Kin_in)
      use xtb_mctc_accuracy, only: wp
      implicit none

      integer  :: nat,i
      real(wp),intent(inout) :: velo(3,nat)
      real(wp) :: mass(nat)
      real(wp) :: E_Kin_in
      real(wp) :: x
      real(wp) :: eperat,v,f,edum
      real(wp) :: Temp
   
   ! nat: number of atoms
   ! E_kin_in: Edum in main program (inner energy)
   ! eperat : Energy per atom
      eperat = E_Kin_in / (3.0_wp * nat)
   
   ! Distribute for all atoms among all 3 spatial coords
      do i=1,nat
         v = sqrt( 2 * eperat / mass(i))
         
         f = 1.0_wp
         call random_number(x)
         if (x > 0.5) f = -1.0_wp
         velo(1,i) = velo(1,i) + v * f

         f = 1.0_wp
         call random_number(x)
         if (x > 0.5) f = -1.0_wp
         velo(2,i) = velo(2,i) + v * f

         f = 1.0_wp
         call random_number(x)
         if (x > 0.5) f = -1.0_wp
         velo(3,i) = velo(3,i) + v * f
      enddo
   
      call ekinet(nat,velo,mass,edum,Temp)
   
      write(*,'('' Generating uniform velocity distribution, T='',&
      &          f10.1)')Temp
   
   end subroutine mdinitu
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! computed nuclear kinetic energy
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ekinet(nat,velo,mass,E_Kin,Temp)
      implicit none
   
      integer :: i
      integer,intent(in) :: nat
      real(wp) :: ms
      real(wp), intent(in)  :: velo(3,nat),mass(nat)
      real(wp), intent(out) :: E_Kin,Temp
   
      E_Kin = 0.0_wp

      do i=1,nat
         ms=mass(i)
         E_Kin = E_Kin + ms * (velo(1,i)**2+velo(2,i)**2+velo(3,i)**2)
      enddo

      E_Kin = E_Kin * 0.5_wp
      Temp  = E_Kin / (0.5_wp * 3.0_wp * nat * kB)

   end subroutine ekinet

end module qcxms_mdinit
