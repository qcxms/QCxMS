
! set masses in vector mass
subroutine setmass(n,iat,mass,imass)
   use xtb_atomic_masses
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   implicit none

   integer  :: n,iat(n),i
   integer  :: imass(n)

   real(wp) :: mass(n)

   do i=1,n
      if(imass(i).lt.0) then
         mass(i)=atomic_mass(iat(i))
      else
         mass(i)=imass(i) * amutoau
      endif
   enddo

end

! read mass from qcxms.in
subroutine read_isotopes(n,mass)
!   use io_reader
   use readcommon
   use xtb_mctc_accuracy, only: wp
   implicit none

   integer :: n
   integer :: nn,i,k,iat
   integer :: mass(n)
   integer :: io_in

   real(wp) :: xx(10)

   mass = -1

   open(file='qcxms.in',newunit=io_in, status='old')

   do
      read(io_in,'(a)',iostat=iocheck)line
      if (iocheck < 0) exit !EOF
      if (iocheck > 0) stop 'ERROR in read_isotopes'

      if(index(line,'isotop').ne.0)then
         write(*,*) 'reading isotopic masses ...'
         call readl(line,xx,nn)
         k = 0
         do i = 1, nn/2
            k = k + 1
            iat = idint(xx(k))
            k = k + 1
            mass(iat) = idint(xx(k))
         enddo
      endif
   enddo

   close(io_in)

end
