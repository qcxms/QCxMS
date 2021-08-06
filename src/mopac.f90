! get gradient with actual coords
subroutine getmopgrad(nat,ic,xyz,g,chrg,etemp, &
  first,edum,achrg,spin)
   use common1
   use newcommon
   use readcommon
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only: toSymbol 
   implicit none

   integer :: chrg
   integer :: nat,i,j,ic(nat),idum,nn

   real(wp) :: g(3,nat),xx(100),edum
   real(wp) :: xyz(3,nat),achrg(nat),spin(nat),etemp
   real(wp) :: dum

   character(len=148) :: btmp
   character(len=30 ) :: methode

   logical :: first

   dum=etemp
   if(dum.lt.0) dum=0.0
   if(ihamilt.eq.1) write(methode,'(''am1 fermi='',f6.0)') dum
   if(ihamilt.eq.2) write(methode,'(''pm3 fermi='',f6.0)') dum
   if(ihamilt.eq.3) write(methode,'(''rm1 fermi='',f6.0)') dum
   if(ihamilt.eq.4) write(methode,'(''pm6-dh2 fermi='',f6.0)') dum
   if(ihamilt.gt.4) stop 'am1 or pm3 or rm1 or pm6 for mopac'

   if(first)then
      if(chrg.eq.1)then
         write(btmp,'(''1scf gradients aux(42,PRECISION=9,MOS=-99996,COMP)'',&
         &'' uhf denout scfcrt=1.d-5 geo-ok xyz spin charge=1 '',a)')trim(methode)
      else
         write(btmp,'(''1scf gradients aux(42,PRECISION=9,MOS=-99999,COMP)'',&
         &'' uhf denout scfcrt=1.d-5 geo-ok xyz spin '',a)')trim(methode)
      endif
   else
      if(chrg.eq.1)then
         write(btmp,'(''1scf gradients aux(42,PRECISION=9,MOS=-99999,COMP)'',&
         &'' uhf denout oldens scfcrt=1.d-5 xyz spin geo-ok charge=1 '',a)')&
         &trim(methode)
      else
         write(btmp,'(''1scf gradients aux(42,PRECISION=9,MOS=-99999,COMP)'',&
         &'' uhf denout oldens scfcrt=1.d-5 xyz spin geo-ok '',a)')trim(methode)
      endif
   endif

   idum=1
   open(unit=3,file='inp')
   write(3,'(a)') trim(btmp)
   write(3,*)
   write(3,*)
   do j=1,nat
      write(3,'(3x,a2,3(f20.14,i5))')toSymbol(ic(j)),&
      &                              autoaa*xyz(1,j),idum,&
      &                              autoaa*xyz(2,j),idum,&
      &                              autoaa*xyz(3,j),idum
   enddo
   close(3)

   calls = calls + 1
   write(line,'(a,''mopac inp 2>/dev/null'')')trim(path)
   call system(line)

   open(unit=3,file='inp.aux')
!3333 read(3,'(a)',end=3334)atmp
   do 
      read(3,'(a)',iostat=iocheck) line
      if (iocheck < 0) exit 
      if(index(line,'HEAT_OF_FORMATION').ne.0)then
         call readl(line,xx,nn)
         edum=xx(nn) * kcaltoau 
      endif

      if(index(line,'ATOM_CHARGES').ne.0)then
         read(3,'(10F18.15)')(achrg(j),j=1,nat)
      endif

      if(index(line,'GRADIENTS:KCAL/MOL/ANGSTROM').ne.0)then
         read(3,'(10F18.15)')((g(i,j),i=1,3),j=1,nat)
         !goto 3334
         exit
      endif
   enddo
!   goto 3333
!3334 close (3)
   close (3)

   g = g/ (autokcal/autoaa) ! / (627.509541d0/0.529177260d0) autokcal * aatoau

   open(unit=4,file='inp.out')
!4444 read(3,'(a)',end=4445)atmp
   do
     read(4,'(a)',iostat=iocheck) line
     if (iocheck < 0) exit 
     if(index(line,'ATOMIC ORBITAL SPIN POPULATIONS').ne.0)then
        read(4,'(a)') line
        read(4,'(a)') line
        read(4,'(a)') line
        do i=1,nat
           read(4,'(a)') line
           call readl(line,xx,nn)
           spin(i) = xx(2)
        enddo
        exit
!      goto 4445
     endif
   enddo
!   goto 4444
!4445 close (4)
   close (4)

  return
end

