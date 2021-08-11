module qcxms_utility
   use common1
   use newcommon
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_constants, only: kB
   use xtb_mctc_convert
   implicit none

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! determine two MOs whose energies add to IEE edum; a third MO (index not used)
   ! represents the virtual MO (vMO) in shake-up to account for right excitation energy

   ! nao = number of AOs (alpha electrons)
   ! emo = orbital energy (eigenvalues)

   subroutine momap(nao,emo,edum,mo1,mo2)
   
      integer  :: nao,mo1,mo2,k,i1,i2 
      integer  :: vmo
   
      real(wp) :: emo(*),edum,dmin,delta,dum
   
      k =0
      i1=0
      i2=0
      dmin=huge(0.0_wp) !1.d+42
   
      do
        mo1=irand(nao)
        mo2=irand(nao)
        vmo=irand(nao/2)+nao
        dum=emo(mo1)

        if ( mo2 > nao / 2 ) then
           mo2 = 0
        else
           dum = dum + emo(mo2)
        endif

        dum = dum + emo(vmo)
        delta = abs(dum - edum)

        if ( delta < dmin ) then
           dmin = delta
           i1 = mo1
           i2 = mo2
        endif

        k = k + 1
        if ( k >= 5000 ) exit 
      enddo
   
      mo1=i1
      mo2=i2
   
   end subroutine momap
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! estimate the electronic temp in Fermi smearing for given IEE and ax
   subroutine setetemp(nfrag,eimp,etemp)
   
      integer  :: nfrag
   
      real(wp) :: eimp,etemp
      real(wp) :: tmp
   
      etemp =  5000. + 20000. * ax
   
      if(eimp.gt.0.and.nfrag.le.1)then
   
         tmp=max(eimp,0.0d0)
   
         etemp = etemp + tmp*ieetemp
   
      endif
   
   end subroutine
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! call the qc code
   subroutine qccall(iprog,fout)
   
      integer :: iprog
   
      character(len=* ) :: fout
      character(len=80) :: atmp
   
      external system
   
      calls = calls + 1
   
   ! DFTB+
      if(iprog.eq.0) write(atmp,'(a,''dftb+ > '',a)') trim(path),trim(fout)
   
   ! TM
      if(iprog.eq.2)then
         if(shell.eq.1) write(atmp,'(''( '',a,''ridft > '',a,'' ) > & /dev/null'')') trim(path),trim(fout)
         if(shell.eq.2) write(atmp,'(''ridft > '',a,'' 2> /dev/null'')')trim(fout)
      endif
   
   ! ORCA
      if(iprog.eq.3) write(atmp,'(''orca ORCA.INPUT > '',a)') trim(fout)
   
   ! MNDO99
      if(iprog.eq.5) write(atmp,'(a,''mndo99 < inp > '',a)') trim(path),trim(fout)
   
      call system(atmp)
   
   ! TM GRAD
      if(iprog.eq.2)then
         if(shell.eq.1) write(atmp,'(''( '',a,''rdgrad >> '',a,'' ) > & /dev/null'')') trim(path),trim(fout)
         if(shell.eq.2) write(atmp,'(''rdgrad >> '',a,'' 2> /dev/null'')')trim(fout)
         call system(atmp)
      endif
   
   end subroutine
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! add a gaussian in imp energy spectrum r() at energy x, width 0.1 eV
   ! just plotting
   subroutine gauss(r,x)
   
      integer :: j
   
      real(wp) :: r(1000),x
      real(wp) :: xmi,xma,st,ee,edif,width
   
      xmi=0
      xma=70
      st=(xma-xmi)/1000
      width=1/1.0**2
      do j=1,1000
         ee=xmi+j*st
         edif=ee-x
         r(j)=r(j)+exp(-width*edif**2)
      enddo
   
   end subroutine
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! create Temporary TMP directories
   subroutine mdtmpdir(it)
      implicit none
   
      integer  :: it
      character(len=80) :: fname
      external :: system
   
      if(it.lt.10000)write(fname,'(''mkdir TMPQCXMS/TMP.'',i4)')it
      if(it.lt.1000) write(fname,'(''mkdir TMPQCXMS/TMP.'',i3)')it
      if(it.lt.100)  write(fname,'(''mkdir TMPQCXMS/TMP.'',i2)')it
      if(it.lt.10)   write(fname,'(''mkdir TMPQCXMS/TMP.'',i1)')it
   
      call system(fname)
   end
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine copytb(it)
      implicit none
   
      integer :: it
      character(len=80) :: fname
      external system
   
      if(it.ge.10000)stop 'error 1 inside copytb'
   
      if(it.ge.1000)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i4)')it
         call system(fname)
   !        write(fname,'(''cp charges.bin TMPQCXMS/TMP.'',i4)')it
   !        call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i4)')it
         call system(fname)
         return
      endif
   
      if(it.ge.100)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i3)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i3)')it
         call system(fname)
         return
      endif
   
      if(it.ge.10)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i2)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i2)')it
         call system(fname)
         return
      endif
   
      if(it.ge.0)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i1)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i1)')it
         call system(fname)
         return
      endif
   
      stop 'error 2 inside copytb'
   
   end subroutine
   
   subroutine copymop(it)
      implicit none
   
      integer :: it
      character(len=80) :: fname
      external system
   
      if(it.ge.10000)stop 'error 1 inside copymop'
   
      if(it.ge.1000)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i4)')it
         call system(fname)
         write(fname,'(''cp inp.den TMPQCXMS/TMP.'',i4)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i4)')it
         call system(fname)
         return
      endif
   
      if(it.ge.100)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i3)')it
         call system(fname)
         write(fname,'(''cp inp.den TMPQCXMS/TMP.'',i3)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i3)')it
         call system(fname)
         return
      endif
   
      if(it.ge.10)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i2)')it
         call system(fname)
         write(fname,'(''cp inp.den TMPQCXMS/TMP.'',i2)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i2)')it
         call system(fname)
         return
      endif
   
      if(it.ge.0)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i1)')it
         call system(fname)
         write(fname,'(''cp inp.den TMPQCXMS/TMP.'',i1)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i1)')it
         call system(fname)
         return
      endif
   
      stop 'error 2 inside copymop'
   
   end subroutine
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine copymsindo(it)
      implicit none
   
      integer :: it
      character(len=80) :: fname
      external system
   
      if(it.ge.10000)stop 'error 1 inside copymsindo'
   
      if(it.ge.1000)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i4)')it
         call system(fname)
         write(fname,'(''cp DENSITY TMPQCXMS/TMP.'',i4)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i4)')it
         call system(fname)
         return
      endif
   
      if(it.ge.100)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i3)')it
         call system(fname)
         write(fname,'(''cp DENSITY TMPQCXMS/TMP.'',i3)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i3)')it
         call system(fname)
         return
      endif
   
      if(it.ge.10)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i2)')it
         call system(fname)
         write(fname,'(''cp DENSITY TMPQCXMS/TMP.'',i2)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i2)')it
         call system(fname)
         return
      endif
   
      if(it.ge.0)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i1)')it
         call system(fname)
         write(fname,'(''cp DENSITY TMPQCXMS/TMP.'',i1)')it
         call system(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i1)')it
         call system(fname)
         return
      endif
   
      stop 'error 2 inside copymsindo'
   
   end subroutine
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! writing routine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine wrstart(it,nat,xyzr,velor,velofr,eimpr,taddr)
      use xtb_mctc_accuracy, only: wp
   ! it is # traj
   ! "nat" is # atoms
   ! xyz coords
   ! velor velocity
   ! velofr is scaling factor
   ! eimpr is energy impact electron ! eimp
   ! taddr is temp add
   
      implicit none
   
      integer  :: it,nat
      integer  :: j
   
      real(wp) :: xyzr (3,nat)
      real(wp) :: velor(3,nat)
      real(wp) :: velofr    (  nat)
      real(wp) :: eimpr,taddr
   
      character(len=80) :: fname
   
      if(it.lt.10000) write(fname,'(''TMPQCXMS/TMP.'',i4,''/qcxms.start'')')it
      if(it.lt.1000)  write(fname,'(''TMPQCXMS/TMP.'',i3,''/qcxms.start'')')it
      if(it.lt.100)   write(fname,'(''TMPQCXMS/TMP.'',i2,''/qcxms.start'')')it
      if(it.lt.10)    write(fname,'(''TMPQCXMS/TMP.'',i1,''/qcxms.start'')')it
   
      open(unit=3,file=fname)
      write(3,'(2i4)') it,nat
      write(3,'(2D22.14)') eimpr,taddr
   
      do j=1,nat
         write(3,'(7D22.14)') xyzr (1:3,j),velor(1:3,j),velofr(j)
      enddo
   
      close(3)
   
   end subroutine
   
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! calculate # of valence electrons in mopac and dftb
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine valel(at,el)
   
      integer  :: at
      real(wp) :: el
   
      if(at.le.2)then
         el=at
      elseif(at.le.10)then
         el=at-2
      elseif(at.le.18)then
         el=at-10
      elseif(at.le.36)then
         el=at-18
         if(at.gt.28) el=at-28
      endif
   
   end subroutine
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! spin value of molecule
   ! returns -1 if no electrons are found (e.g. H+)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine getspin(nat,ic,chrg,isp)
   
      integer nat,i,j,chrg,isp,ic(nat)
   
      j=0
   
      do i = 1, nat
         j = j + ic(i)
      enddo
   
      j = j - chrg
      isp = 1 + mod(j,2)
   
      if ( j < 1 ) isp = -1
   
   end subroutine
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Boltzmann population for temp. t and energies e
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine boltz(units,nfrag,temp,ip,fchrg)
   
      integer  :: nfrag,units,i
   
      real(wp) :: ip(nfrag),fchrg(nfrag)
      real(wp) :: temp,f,esum,const
   
      ! kcal/mol
      if (units==1) const = autokcal
      ! eV
      if (units==2) const = autoev
   
      f = temp * kB * const
   
      esum = 0
   
      do i = 1, nfrag
         esum = esum + exp(-ip(i)/f)
      enddo
   
      do i = 1, nfrag
         fchrg(i) = exp(-ip(i)/f)/esum
      enddo
   
   end subroutine
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! integer random number n=<irand<=1
   function irand(n) result(rnd)
   integer  :: n, rnd
   real(wp) :: x,nx
   
   call random_number(x)
   nx = n * x + 1
   rnd = int(nx)
   if ( rnd > n ) rnd=n
   
   end
   
end module qcxms_utility
