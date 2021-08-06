! read DFTB+ E and forces/gradient
subroutine dftbgrad(nat,g,chrg,spin,edum)
   use readcommon
   use xtb_mctc_accuracy, only: wp
!      use io_reader

   integer nat,j,nn

   real(wp) :: g(3,nat),xx(5),edum,chrg(nat),spin(nat)


   chrg = 0
   spin = 0
   open(unit=33,file='detailed.out')
!10 read(33,'(a)',end=100)line
   do
      read(33,'(a)',iostat=iocheck)line
      if (iocheck < 0) exit

      if(index(line,'Total Mermin free energy:').ne.0)then
         call readl(line,xx,nn)
         edum =xx(1)
      endif

      if(index(line,' Total Forces').ne.0)then
         do j=1,nat
            read(33,*)g(1,j),g(2,j),g(3,j)
         enddo
         exit
         !goto 100
      endif

      if(index(line,'Net atomic charges').ne.0)then
         read(33,'(a)',iostat=iocheck)line
         if (iocheck < 0) exit
         do j=1,nat
         !   read(33,'(a)',end=100)line
            read(33,'(a)',iostat=iocheck)line
            if (iocheck < 0) exit
            call readl(line,xx,nn)
            chrg(j)=xx(nn)
            spin(j)=0
         enddo
         if (iocheck < 0) exit
      endif
      if(index(line,'Atom populations (up').ne.0)then
         !read(33,'(a)',end=100)line
         read(33,'(a)',iostat=iocheck)line
         if (iocheck < 0) exit
         do j=1,nat
         !   read(33,'(a)',end=100)line
            read(33,'(a)',iostat=iocheck)line
            if (iocheck < 0) exit
            call readl(line,xx,nn)
            chrg(j)=xx(nn)
            spin(j)=xx(nn)
         enddo
         if (iocheck < 0) exit
      endif
      if(index(line,'Atom populations (down').ne.0)then
         !read(33,'(a)',end=100)line
         read(33,'(a)',iostat=iocheck)line
         if (iocheck < 0) exit
         do j=1,nat
         !   read(33,'(a)',end=100)line
            read(33,'(a)',iostat=iocheck)line
            if (iocheck < 0) exit
            call readl(line,xx,nn)
            chrg(j)=chrg(j)+xx(nn)
            spin(j)=spin(j)-xx(nn)
         enddo
         if (iocheck < 0) exit
      endif
   enddo
!   goto 10
!100 continue
   close(33)
! grad
   g = - g
   return
end

subroutine dftbenergy(edum)
   use xtb_mctc_accuracy, only: wp
   use readcommon
!      use io_reader
   integer nn
   real(wp) :: edum,xx(10)
   edum = 0
   open(unit=33,file='detailed.out')
!10 read(33,'(a)',end=100)line
   do
      read(33,'(a)',iostat=iocheck)line
      if (iocheck < 0) exit
      if (index(line,'Total Mermin').ne.0)then
         call readl(line,xx,nn)
         edum =xx(1)
         exit
         !goto 100
      endif
   enddo
!   goto 10
!100 close(33)
   close(33)
end

subroutine dftbout(n,xyz,ic,chrg,isp,rest,etemp,conv)
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only: toSymbol 
   use common1
   use qcxms_utility, only: getspin
   implicit none

   integer  :: idum(100),idum2(100),i,j,k,m,ich
   integer  :: chrg,n,ic(*),isp

   real(wp) :: xyz(3,*),etemp,conv

   character(len=80) :: s,p,a

   logical uhf
   logical rest


   if(isp.eq.0) call getspin(n,ic,chrg,isp)

   uhf=.false.
   if(isp.gt.1.or.chrg.gt.0)uhf=.true.

   ich=143

   open(unit=ich,file='dftb_in.hsd')
   write(ich,'(''Geometry = GenFormat {'')')
   write(ich,'(i4,'' C'')')n
   idum = 0
   do i=1,n
      idum(ic(i))=idum(ic(i))+1
   enddo
   k=2
   m=0
   s='                                       '
   do i=1,100
      if(idum(i).ne.0)then
         m=m+1
         idum2(i)=m
         s(k:k+2)=toSymbol(i)
         k=k+3
      endif
   enddo
   write(ich,'(A)')s

   do 200 i=1,n
      write(ich,'(2i4,3F22.12)')&
      &i,idum2(ic(i)),autoaa*xyz(1,i),autoaa*xyz(2,i),autoaa*xyz(3,i)
200 continue

   j=ich
   write(j,*)'}'
   write(j,*)'Options = {'
   write(j,*)'CalculateForces = Yes'
   write(j,*)'}'


   write(j,*)'Hamiltonian = DFTB {'
   if(etemp.gt.0)then
      write(j,*)'Filling = Fermi {'
      write(j,*)'Temperature [K] = ',idint(etemp),' }'
   endif
   write(j,*)'Charge = ',chrg
   if(uhf)&
   &write(j,'('' SpinPolarisation=Colinear{UnpairedElectrons ='',F4.1,&
   &       ''}'')')float(isp-1)
   write(j,*)'SCC = Yes'
   write(j,'(''SCCTolerance ='',E12.5)')conv
   write(j,*)'ThirdOrderFull = Yes'
   write(j,*)'DampXH = Yes'
   write(j,*)'DampXHExponent = 4.20'
   write(j,*)'HubbardDerivs {'
! value from M. Kubillus (SG, 10/13)
! use -0.07 for P and S in case of SCC problems
   if(idum(1).ne.0)write(j,*)'H =   -0.1857'
   if(idum(5).ne.0)write(j,*)'B =   -0.1450'
   if(idum(6).ne.0)write(j,*)'C =   -0.1492'
   if(idum(7).ne.0)write(j,*)'N =   -0.1535'
   if(idum(8).ne.0)write(j,*)'O =   -0.1575'
   if(idum(9).ne.0)write(j,*) 'F =  -0.1623'
   if(idum(11).ne.0)write(j,*)'Na=  -0.0350'
   if(idum(12).ne.0)write(j,*)'Mg=  -0.0200'
   if(idum(15).ne.0)write(j,*)'P =  -0.1400'
   if(idum(16).ne.0)write(j,*)'S =  -0.1100'
   if(idum(17).ne.0)write(j,*)'CL = -0.0697'
   if(idum(35).ne.0)write(j,*)'BR = -0.0573'
   if(idum(53).ne.0)write(j,*)'I  = -0.0433'
   write(j,*)' }'

   if(rest)then
      write(j,*)'ReadInitialCharges = yes'
   else
      write(j,*)'ReadInitialCharges = no'
   endif
   write(j,*)'MaxSCCIterations = 200'
!     write(j,*)'Mixer=anderson{}'         WORSE
   write(j,*)'MaxAngularMomentum = {'
   if(idum(1).ne.0) write(j,*)'H = "s"'
   if(idum(5).ne.0) write(j,*)'B = "p"'
   if(idum(6).ne.0) write(j,*)'C = "p"'
   if(idum(7).ne.0) write(j,*)'N = "p"'
   if(idum(8).ne.0) write(j,*)'O = "p"'
   if(idum(9).ne.0) write(j,*)'F = "p"'
   if(idum(15).ne.0)write(j,*)'P = "d"'
   if(idum(16).ne.0)write(j,*)'S = "d"'
   if(idum(17).ne.0)write(j,*)'CL= "d"'
   if(idum(35).ne.0)write(j,*)'BR= "d"'
   if(idum(53).ne.0)write(j,*)'I= "d"'
   write(j,*)' }'


   p='"/usr/local/dftb+/'
   write(j,*)'SlaterKosterFiles = {'

   if(hhmod)then
      if(idum(1).ne.0)write(j,*)'h-h =',trim(p),'hh-mod.spl"'
   else
      if(idum(1).ne.0)write(j,*)'h-h =',trim(p),'hh.spl"'
   endif

   k=j
   do i=1,54
      do j=1,i
         if(i.eq.1.and.j.eq.1) cycle
         if(i.eq.j) then
            call splstr(p,i,i,a)
            if(idum(i).gt.0)write(k,*) trim(a)
         else
            call splstr(p,i,j,a)
            if(idum(i).gt.0.and.idum(j).gt.0)write(k,*) trim(a)
            call splstr(p,j,i,a)
            if(idum(i).gt.0.and.idum(j).gt.0)write(k,*) trim(a)
         endif
      enddo
   enddo

   write(k,*)' }'

   if(uhf)then
      write(k,*)'SpinConstants = {'
      if(idum(1).ne.0) write(k,*)'H={ -0.072 }'
      if(idum(6).ne.0) write(k,*)'C={ -0.031 -0.025 -0.025 -0.023 }'
      if(idum(7).ne.0) write(k,*)'N={ -0.033 -0.027 -0.027 -0.026 }'
      if(idum(8).ne.0) write(k,*)'O={ -0.035 -0.030 -0.030 -0.028 }'
!sg extrapolated or taken from CNOF
      if(idum(9).ne.0) write(k,*)'F={ -0.037 -0.033 -0.033 -0.030 }'
      if(idum(15).ne.0)write(k,*)&
      &'P ={ -0.033 -0.027 -0.026 -0.026 0.0 0.0 0.0 0.0 0.0 }'
      if(idum(16).ne.0)write(k,*)&
      &'S ={ -0.035 -0.030 -0.030 -0.030 0.0 0.0 0.0 0.0 0.0 }'
      if(idum(17).ne.0)write(k,*)&
      &'CL={ -0.037 -0.033 -0.033 -0.030 0.0 0.0 0.0 0.0 0.0 }'
      if(idum(35).ne.0)write(k,*)&
      &'BR={ -0.037 -0.033 -0.033 -0.030 0.0 0.0 0.0 0.0 0.0 }'
      if(idum(53).ne.0)write(k,*)&
      &'I={ -0.037 -0.033 -0.033 -0.030 0.0 0.0 0.0 0.0 0.0 }'
      write(k,*)' }'
   endif

   write(ich,*)'}'

   close(ich)
   return
end

subroutine splstr(p,i,j,c)
   use xtb_mctc_accuracy, only: wp

   integer :: i,j

   character(len=* ) :: p
   character(len=80) :: a,b,c
   CHARACTER(len=2 ) :: ELEMNT(107), aa,bb

   DATA ELEMNT/'h ','he',&
   &'li','be','b ','c ','n ','o ','f ','ne',&
   &'na','mg','al','si','p ','s ','cl','ar',&
   &'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
   &'zn','ga','ge','as','se','br','kr',&
   &'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',&
   &'cd','in','sn','sb','te','i ','xe',&
   &'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
   &'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
   &'au','hg','tl','pb','bi','po','at','rn',&
   &'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',&
   &'fm','md','cb','xx','xx','xx','xx','xx'/
   aa=ELEMNT(I)
   bb=ELEMNT(J)

   if(aa(2:2).eq.' '.and.bb(2:2).eq.' ')&
   &write(a,'(a1,a1,''.spl'')')aa,bb
   if(aa(2:2).ne.' '.and.bb(2:2).ne.' ')&
   &write(a,'(a2,a2,''.spl'')')aa,bb
   if(aa(2:2).ne.' '.and.bb(2:2).eq.' ')&
   &write(a,'(a2,a1,''.spl'')')aa,bb
   if(aa(2:2).eq.' '.and.bb(2:2).ne.' ')&
   &write(a,'(a1,a2,''.spl'')')aa,bb

   if(aa(2:2).eq.' '.and.bb(2:2).eq.' ')&
   &write(b,'(a1,''-'',a1)')aa,bb
   if(aa(2:2).ne.' '.and.bb(2:2).ne.' ')&
   &write(b,'(a2,''-'',a2)')aa,bb
   if(aa(2:2).eq.' '.and.bb(2:2).ne.' ')&
   &write(b,'(a1,''-'',a2)')aa,bb
   if(aa(2:2).ne.' '.and.bb(2:2).eq.' ')&
   &write(b,'(a2,''-'',a1)')aa,bb

   write(c,'(a,'' ="'',a,a,''"'')')trim(b),trim(p),trim(a)

end
