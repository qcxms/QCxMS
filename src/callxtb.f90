subroutine callxtb(nuc,xyz,iat,chrg,spin,etemp,E,g,achrg,aspin)
!   use io_reader
   use newcommon
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_symbols, only: toSymbol 
   implicit none

   integer :: j
   integer :: nn,i
   integer :: nue
   integer :: nuc
   integer :: iat(nuc)
   integer :: chrg,spin !takes in multiplicity
   integer :: iocheck

   real(wp) :: achrg(nuc),aspin(nuc)
   real(wp) :: xyz(3,nuc)
   real(wp) :: etemp
   real(wp) :: xx(100),edum
   real(wp) :: E,g(3,nuc)

   character(len=2 ) :: asym
   character(len=90) :: atmp

   logical :: ex


! calculate number of electrons
   nue=int(dble(spin)/2.0d0)

! aspin (for now - xtb doesnt print it)
   aspin = 0.0d0

! check if the file exist
   inquire(file='inp',exist=ex)
   if(ex)call system('rm inp')

   inquire(file='gradient',exist=ex)
   if(ex)call system('rm gradient')

   inquire(file='charges',exist=ex)
   if(ex)call system('rm charges')

! xtb reads coord file in tm format
! thus write inp file for xtb in AU
   open(unit=3,file='inp')
   write(3,*) '$coord'
   do j=1,nuc
      write(3,'(3x,3(f20.14),3x,a2)')xyz(1,j),xyz(2,j),xyz(3,j),toSymbol(iat(j))
   enddo
   write(3,*)'$end'
   close(3)
!ccccccccccccccccccccccccccc
! CALCULATION
!cccccccccccccccccccccccccc
! system input line
   write(atmp,'(a,''xtb inp -scc -uhf '',i2,'' -chrg '',i2,'' -etemp '',f8.2,''      -grad  > job.last'')')&
     & trim(xtbpath),nue,chrg,etemp
!      write(*,*) atmp

!  call xtb
   call system(atmp)

!cccccccccccccccccccccccccccc
!      READ OUTPUT
!cccccccccccccccccccccccccccc

   open(3,file='gradient')


! read energies and forces
   do
   read(3,'(a)',iostat=iocheck)atmp
     if (iocheck>0)then     !Fail
       write(*,*) 'Something is wrong in the input. Exiting...'
       stop
     elseif (iocheck<0)then !EOF
       exit
     else
       if(index(atmp,'SCF energy =').ne.0)then
          call readl(atmp,xx,nn)
          edum=xx(nn-1)
          do i=1,nuc
             read(3,'(a)')atmp
          end do
          do i=1,nuc
             read(3,*)g(1:3,i)
          enddo
          exit
       endif
     endif
   enddo
   close (3)
   E = edum

! move gradient to gradient.last - to be sure that the file gets overwritten
! xtb works in a way where it adds subsequent calculations to grad
   CALL RENAME('gradient','gradient.last')

! read atomic charges
   open(3,file='charges')
   do i=1,nuc
      read(3,*)achrg(i)
   enddo
   CALL RENAME('charges','charges.last')



end subroutine
