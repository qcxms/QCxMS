module qcxms_utility
   use common1
   use newcommon
   use mctc_io, only : to_symbol, read_structure, write_structure
   use mctc_io_filetype, only : filetype
   use mctc_io_structure, only : structure_type
   use mctc_env_error, only : error_type
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_constants, only: kB
   use xtb_mctc_convert
   implicit none

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! determine two MOs whose energies add to IEE edum; a third MO (index not used)
   ! represents the virtual MO (vMO) in shake-up to account for right excitation energy

   ! ihomo = number (integer) of HOMO
   ! emo = orbital energy (eigenvalues)
   ! edum = total energy (ehomo + exc)

   subroutine momap(ihomo,emo,edum,mo1,mo2)
   
      integer  :: ihomo,mo1,mo2,k,i1,i2 
      integer  :: vmo
   
      real(wp) :: emo(*),edum,dmin,delta,dum
   
      k  = 0
      i1 = 0
      i2 = 0
      dmin = huge(0.0_wp) 
   
      do while (k <= 5000)

        mo1 = irand(ihomo) ! irand -> randomization function (see below) 
        mo2 = irand(ihomo)
        vmo = irand(ihomo/2) + ihomo
        dum = emo(mo1)

        ! if electron is unpaired 
        if ( mo2 > ihomo / 2 ) then
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

      enddo
   
      mo1 = i1
      mo2 = i2
   
   end subroutine momap
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! estimate the electronic temp in Fermi smearing for given IEE and ax
   subroutine setetemp(nfrag,eimp,etemp)
   
      integer  :: nfrag
   
      real(wp) :: eimp,etemp
      real(wp) :: tmp
   
      etemp =  5000. + 20000. * ax
   
      if (eimp > 0 .and. nfrag <= 1) then
   
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
   
      call execute_command_line(atmp)
   
   ! TM GRAD
      if(iprog.eq.2)then
         if(shell.eq.1) write(atmp,'(''( '',a,''rdgrad >> '',a,'' ) > & /dev/null'')') trim(path),trim(fout)
         if(shell.eq.2) write(atmp,'(''rdgrad >> '',a,'' 2> /dev/null'')')trim(fout)
         call execute_command_line(atmp)
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
   
      if (it < 10000) write(fname,'(''mkdir TMPQCXMS/TMP.'',i4)') it
      if (it < 1000)  write(fname,'(''mkdir TMPQCXMS/TMP.'',i3)') it
      if (it < 100)   write(fname,'(''mkdir TMPQCXMS/TMP.'',i2)') it
      if (it < 10)    write(fname,'(''mkdir TMPQCXMS/TMP.'',i1)') it
   
      call execute_command_line(fname)
   end
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine copytb(it)
      implicit none
   
      integer :: it
      character(len=80) :: fname
   
      if(it.ge.10000)stop 'error 1 inside copytb'
   
      if(it.ge.1000)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i4)')it
         call execute_command_line(fname)
   !        write(fname,'(''cp charges.bin TMPQCXMS/TMP.'',i4)')it
   !        call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i4)')it
         call execute_command_line(fname)
         return
      endif
   
      if(it.ge.100)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i3)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i3)')it
         call execute_command_line(fname)
         return
      endif
   
      if(it.ge.10)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i2)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i2)')it
         call execute_command_line(fname)
         return
      endif
   
      if(it.ge.0)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i1)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i1)')it
         call execute_command_line(fname)
         return
      endif
   
      stop 'error 2 inside copytb'
   
   end subroutine
   
   subroutine copymop(it)
      implicit none
   
      integer :: it
      character(len=80) :: fname
   
      if(it.ge.10000)stop 'error 1 inside copymop'
   
      if(it.ge.1000)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i4)')it
         call execute_command_line(fname)
         write(fname,'(''cp inp.den TMPQCXMS/TMP.'',i4)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i4)')it
         call execute_command_line(fname)
         return
      endif
   
      if(it.ge.100)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i3)')it
         call execute_command_line(fname)
         write(fname,'(''cp inp.den TMPQCXMS/TMP.'',i3)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i3)')it
         call execute_command_line(fname)
         return
      endif
   
      if(it.ge.10)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i2)')it
         call execute_command_line(fname)
         write(fname,'(''cp inp.den TMPQCXMS/TMP.'',i2)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i2)')it
         call execute_command_line(fname)
         return
      endif
   
      if(it.ge.0)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i1)')it
         call execute_command_line(fname)
         write(fname,'(''cp inp.den TMPQCXMS/TMP.'',i1)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i1)')it
         call execute_command_line(fname)
         return
      endif
   
      stop 'error 2 inside copymop'
   
   end subroutine
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine copymsindo(it)
      implicit none
   
      integer :: it
      character(len=80) :: fname
   
      if(it.ge.10000)stop 'error 1 inside copymsindo'
   
      if(it.ge.1000)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i4)')it
         call execute_command_line(fname)
         write(fname,'(''cp DENSITY TMPQCXMS/TMP.'',i4)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i4)')it
         call execute_command_line(fname)
         return
      endif
   
      if(it.ge.100)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i3)')it
         call execute_command_line(fname)
         write(fname,'(''cp DENSITY TMPQCXMS/TMP.'',i3)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i3)')it
         call execute_command_line(fname)
         return
      endif
   
      if(it.ge.10)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i2)')it
         call execute_command_line(fname)
         write(fname,'(''cp DENSITY TMPQCXMS/TMP.'',i2)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i2)')it
         call execute_command_line(fname)
         return
      endif
   
      if(it.ge.0)then
         write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i1)')it
         call execute_command_line(fname)
         write(fname,'(''cp DENSITY TMPQCXMS/TMP.'',i1)')it
         call execute_command_line(fname)
         write(fname,'(''cp coord TMPQCXMS/TMP.'',i1)')it
         call execute_command_line(fname)
         return
      endif
   
      stop 'error 2 inside copymsindo'
   
   end subroutine
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! writing routine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine wrstart(traj,mol,nat,xyzr,velor,velofr,eimpr,taddr)

   ! traj is # traj
   ! "nat" is # atoms
   ! xyz coords
   ! velor velocity
   ! velofr is scaling factor
   ! eimpr is energy impact electron ! eimp
   ! taddr is temp add
      integer  :: traj,nat
      integer  :: j
      integer  :: io_wr_xyz, io_wr_info
   
      real(wp) :: xyzr (3,nat)
      real(wp) :: velor(3,nat)
      real(wp) :: velofr    (  nat)
      real(wp) :: eimpr,taddr
   
      character(len=80) :: fname_xyz, fname_info
   
      type(error_type), allocatable :: error
      type(structure_type) :: mol

      mol%xyz = xyzr
      mol%nat = nat

      if(traj < 10000) write(fname_xyz,'(''TMPQCXMS/TMP.'',i4,''/start.xyz'')')traj
      if(traj < 1000)  write(fname_xyz,'(''TMPQCXMS/TMP.'',i3,''/start.xyz'')')traj
      if(traj < 100)   write(fname_xyz,'(''TMPQCXMS/TMP.'',i2,''/start.xyz'')')traj
      if(traj < 10)    write(fname_xyz,'(''TMPQCXMS/TMP.'',i1,''/start.xyz'')')traj
   
      call write_structure(mol, fname_xyz, error, filetype%xyz) 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(traj < 10000) write(fname_info,'(''TMPQCXMS/TMP.'',i4,''/qcxms.start'')')traj
      if(traj < 1000)  write(fname_info,'(''TMPQCXMS/TMP.'',i3,''/qcxms.start'')')traj
      if(traj < 100)   write(fname_info,'(''TMPQCXMS/TMP.'',i2,''/qcxms.start'')')traj
      if(traj < 10)    write(fname_info,'(''TMPQCXMS/TMP.'',i1,''/qcxms.start'')')traj

      open(file=fname_info, newunit= io_wr_info)
      write(io_wr_info,'(i4)') traj
      write(io_wr_info,'(2D22.14)') eimpr,taddr

      do j=1,nat
         write(io_wr_info,'(7D22.14)') velor(1:3,j),velofr(j)
      enddo

      close(io_wr_info)
   
   end subroutine wrstart
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read qcxms.start file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rdstart(itrj,mol,nat,xyz,velo,velof,tadd,eimp)
  
     integer  :: itrj,nat
     integer  :: j,ndum
     integer  :: io_info, io_xyz
     integer  :: ierror
      integer,allocatable  :: iat (:)

     real(wp) :: xyz (3,nat)
     real(wp) :: velo(3,nat)
     real(wp) :: velof   (nat)
     real(wp) :: tadd,eimp
  
     character(len=80) :: fname

     type(error_type), allocatable :: error
     type(structure_type), intent(out) :: mol

     call read_structure(mol, 'start.xyz', error, filetype%xyz)

     ndum = mol%nat
     if (ndum /= nat) stop '- error in rdstart -'
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     fname='qcxms.start'
  
     open(file=fname, newunit=io_info, status='old', &
       action='read', iostat=ierror)

     if(ierror > 0) stop ' - Missing qcxms.start file! -'

     read (io_info,*) itrj
     read (io_info,'(2D22.14)') eimp, tadd

     do j = 1, nat
        read (io_info,'(7D22.14)') velo(1:3,j), velof(j)
     enddo
 
     close(io_info)

     !fname='qcxms.start'
  
     !open(file=fname,newunit=io_start,status='old', &
     !  action='read',iostat=ierror)

     !if(ierror > 0) stop ' - Missing qcxms.start file! -'

     !read (io_start,*) itrj,ndum
     !read (io_start,'(2D22.14)') eimp,tadd
  
     !if (ndum /= nat) stop '- error in rdstart -'
  
     !do j = 1, nat
     !   read (io_start,'(7D22.14)') xyz(1:3,j), velo(1:3,j), velof(j)
     !enddo
  
     !close(io_start)

  end subroutine rdstart
   
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
   
   end subroutine valel

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
   
      j = j - abs(chrg)
      isp = 1 + mod(j,2)
      
      if ( j < 1 ) isp = -1
   
   end subroutine getspin
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Boltzmann population for temp. t and energies e
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine boltz(units,nfrag,mchrg,temp,ip,fragchrg2)
   
      integer  :: nfrag,units,i, j, mchrg
   
      real(wp),intent(in)  :: ip(nfrag,mchrg)
      real(wp),intent(out) :: fragchrg2(nfrag,mchrg)
      real(wp) :: temp,f,esum,const
   
      ! kcal/mol
      if (units==1) const = autokcal
      ! eV
      if (units==2) const = autoev
   
      f = temp * kB * const

      esum = 0
   
      do i = 1, nfrag
        do j = 1, mchrg
          esum = esum + exp(-ip(i,j)/f)
        enddo
      enddo
   
      do i = 1, nfrag
        do j = 1, mchrg
          fragchrg2(i,j) = exp(-ip(i,j)/f) / esum
        enddo
      enddo
   
   end subroutine boltz
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! integer random number n=<irand<=1
   function irand(n) result(rnd)
   integer  :: n, rnd
   real(wp) :: x,nx
   
   call random_number(x)
   nx = n * x + 1
   rnd = int(nx)
   if ( rnd > n ) rnd = n
   
   end function irand
   
end module qcxms_utility
