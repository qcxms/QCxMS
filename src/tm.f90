module qcxms_use_turbomole
  use common1
  use readcommon
  use qcxms_utility, only: getspin
  use xtb_mctc_accuracy, only:wp
  use xtb_mctc_symbols, only: toSymbol 
  implicit none

  contains 

  subroutine rdtmgrad(nat,g,chrg,spin,edum)

     integer  :: nat,i,j
     integer  :: nn,io_grad
     real(wp) :: g(3,nat),xx(10),edum,chrg(nat),spin(nat)
  
     open(file='gradient',newunit=io_grad)
     i=0
     do
       read(io_grad,'(a)',iostat=iocheck)line
       if (iocheck > 0)then 
         write(*,*) 'Something is wrong in tm.f90. Exiting...'
       elseif( iocheck < 0) then !EOF
         exit
       else
  
        if(index(line,'energy').ne.0)then
           call readl(line,xx,nn)
           edum =xx(nn-1)
        endif
        if(index(line,'$end').eq.0)then
           i=i+1
        else
           rewind io_grad
           exit
        endif
       endif
     enddo

  !     write(*,*)'line',i
     outer: do
       do j=1,i-nat
          read(io_grad,'(a)',iostat=iocheck)line
          if (iocheck<0)exit outer !EOF
       end do
       do j=1,nat
          read(io_grad,*)g(1,j),g(2,j),g(3,j)
       end do
     end do  outer
     close(io_grad)
  
     open(file='job.last',newunit=io_grad)
     do 
       read(io_grad,'(a)',iostat=iocheck)line
       if (iocheck > 0) error stop 'Something is wrong in tm.f90. Exiting...'
       if( iocheck < 0) exit
  
       if(index(line,'atomic populations from total density:').ne.0)then
          read(io_grad,'(a)',iostat=iocheck)line
          if( iocheck < 0) exit
          read(io_grad,'(a)',iostat=iocheck)line
          if( iocheck < 0) exit 
          do i=1,nat
             read(io_grad,'(a)',iostat=iocheck)line
             if( iocheck < 0) exit 
             call readl(line,xx,nn)
             chrg(i)=xx(2)
          enddo
       endif
       if(index(line,'Unpaired electrons ').ne.0)then
          read(io_grad,'(a)',iostat=iocheck)line
          if( iocheck < 0) exit 
          read(io_grad,'(a)',iostat=iocheck)line
          if( iocheck < 0) exit 
          spin(1:nat)=0
          do i=1,nat
             read(io_grad,'(a)',iostat=iocheck)line
             if( iocheck < 0) exit 
             call readl(line,xx,nn)
             if(nn.lt.2) exit
             if(nn.gt.1)spin(int(xx(1)))=xx(2)
          enddo
       endif
     enddo
     close(io_grad)
  
     return
  end subroutine rdtmgrad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

  subroutine rdtmenergy(fname,edum)
     integer  :: nn
     integer  :: io_tm 
     real(wp) :: edum,xx(10)
     character(len=*) :: fname
  
     open(file=fname, newunit=io_tm)
     edum=0.0d0
     read(io_tm,'(a)')line
     read(io_tm,'(a)')line
     call readl(line,xx,nn)
     edum=xx(2)
     close(io_tm)
  
  end subroutine rdtmenergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  
  subroutine setfermi(temp)
     integer  :: io_control 
     real(wp) :: temp
     character(len=80) :: atmp
  
  !      write(atmp,'(''$fermi tmstrt='',F8.1,'' tmend='',F8.1)')temp,temp
     write(atmp,'(''$fermi tmstrt='',F8.1,'' hlcrt=-1.0E0  stop=1.E-99 addTS noerf '')')temp
  
     call execute_command_line('kdg fermi')
     call execute_command_line('kdg end')
     open(file='control', newunit=io_control, ACCESS='APPEND')
     write(io_control,'(a)')atmp
     write(io_control,'(''$end'')')
     close(io_control)
  
  end subroutine setfermi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  
  subroutine dscftm

     if(shell.eq.1) &
       call execute_command_line('( /usr/local/bin/ridft >  job.last ) > & /dev/null')
     if(shell.eq.2) &
       call execute_command_line('( ridft >  job.last  2>   /dev/null')
  
  end subroutine dscftm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  
  subroutine gradtm

  !     if(shell.eq.1) call execute_command_line('( grad >> job.last ) > & /dev/null')
  !     if(shell.eq.2) call execute_command_line('  grad >> job.last  2>   /dev/null')
  
     if(shell.eq.1) &
       call execute_command_line('( /usr/local/bin/rdgrad >> job.last ) > & /dev/null')
     if(shell.eq.2) &
       call execute_command_line('  rdgrad >> job.last  2>   /dev/null')
  
  end subroutine gradtm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  
  subroutine initm(n,iat,chrg,spin)
  
     integer  :: n,iat(n),chrg,spin
     integer  :: idum,jdum,io,maxiter,intmem,nopen
     integer  :: scfconv,i

     real(wp) :: thize,thime
  
     character(len=20) :: atmp,gridi
     character(len=30) :: basi
  
     logical :: strange_elem
  
  ! set strange element to false
     strange_elem=.false.
  ! save info
     jdum=spin
     call getspin(n,iat,chrg,idum)
     nopen=idum/2
     spin=idum
  ! check for strange elements (Cu and Pd)
     do i=1,n
        if(iat(i).eq.29.or.iat(i).eq.46) strange_elem=.true.
     enddo
  
     intmem=3000
     maxiter=400
     scfconv=6
     gridi='m4'
     basi='TZVP'

     if(bas.eq.1)basi='SV'
     if(bas.eq.2)basi='SV'
     if(bas.eq.3)basi='SV(P)'
     if(bas.eq.4)basi='SVP'
     if(bas.eq.5)basi='TZVP'
     if(bas.eq.6)basi='ma-def2-svp'
     if(bas.eq.7)basi='ma-def2-tzvp'
     if(bas.eq.8)basi='ma-def2-tzvpp'
     if(bas.eq.9)basi='def2-SV(P)'
     if(bas.eq.10)basi='def2-SVP'
     if(bas.eq.11)basi='def2-TZVP'
     if(bas.eq.12)basi='def2-QZVP'
     if(bas.eq.13)basi='QZVP'
     if(bas.eq.14)basi='ma-def2-QZVP'

     io =11
     thize=1.d-9
     thime=1

     open(unit=io,file='define.inp')
     write(io,*)'    '
     write(io,*)'    '
     write(io,'(''a coord'')')
     write(io,'('' sy c1'')')
     write(io,*)'*'
     write(io,'('' no '')')
     write(io,'('' b all '',a20)') basi
     write(io,*)'*'
     write(io,*)'eht '
     write(io,*)'    '
     if(strange_elem) write(io,*)'    '
     write(io,*)chrg
     if(n.eq.1)write(io,*)'    '
     write(io,*)'n'
     if(nopen.ne.0)then
        write(io,*)'uf ',nopen
     else
        write(io,*)'s'
     endif
     write(io,*)'*'
     write(io,*)'    '
     write(io,*)'    '
     write(io,*)'    '
     write(io,*)'    '
     write(io,*)'    '
     write(io,*)'    '
     write(io,*)'    '
     write(io,*)'dft'
     write(io,*)'on '
     write(io,*)'func'
     write(io,*)'b-p'
     write(io,*)'grid'
     write(io,'(a20)') gridi
     write(io,*)'q'
     write(io,*)'scf'
     write(io,*)'iter'
     write(io,*)maxiter
     write(io,*)'conv'
     write(io,*)scfconv
     write(io,*)'thi'
     write(io,*)thize
     write(io,*)int(thime)
     write(io,*)'ints'
     write(io,*)'y'
     write(io,*)intmem,' twoint'
     write(io,*)'*'
     write(io,*)'    '
     write(io,*)'q'
     close(io)
  
     call execute_command_line('define < define.inp > define.out')
  !     now the setup: hybrid func. b3lyp, grid m4 (savings with m3 are marginal
  !     and m3 grid with fermi smearing produces a noisy gradient!). Real savings
  !     come with rij, extol 2.500, and scfonv 6. CAB 6.10.15
  !     gridsize m3 for GGAs - tests June 2016 CAB
     call execute_command_line('kdg scfdump')
     call execute_command_line('kdg end    ')
     call execute_command_line('kdg dft    ')
     call execute_command_line("echo '$dft' >> control")
  !     PBE0
     if(func.eq.0) then
        call execute_command_line("echo ' functional pbe0'   >> control")
        call execute_command_line("echo ' gridsize m4      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$disp3 -bj        '   >> control")
        call execute_command_line("echo '$extol 2.500      '   >> control")
        call execute_command_line("echo '$end              '   >> control")
  !     PBE12 and PBE38 do not work
     elseif(func.eq.2.or.func.eq.3.or.func.eq.9.or.func.eq.13) then
        error stop 'PBE12/PBE38/B3PW91/REVPBE not implemented with Turbomole.'
  !     M062X without dispersion!
     elseif(func.eq.4) then
        call execute_command_line("echo ' functional m062x'   >> control")
        call execute_command_line("echo ' gridsize m4      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$extol 2.500      '   >> control")
        call execute_command_line("echo '$end              '   >> control")
  !     PBE-D3BJ
     elseif(func.eq.5) then
        call execute_command_line("echo ' functional pbe'   >> control")
        call execute_command_line("echo ' gridsize m3      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$disp3 -bj        '   >> control")
        call execute_command_line("echo '$end              '   >> control")
  !     B97-D3BJ. CAB 27.10.15
     elseif(func.eq.6) then
        call execute_command_line("echo ' functional b97-d'   >> control")
        call execute_command_line("echo ' gridsize m3      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$disp3 -bj        '   >> control")
        call execute_command_line("echo '$end              '   >> control")
  !     B3LYP-D3BJ
     elseif(func.eq.7)  then
        call execute_command_line("echo ' functional b3-lyp'   >> control")
        call execute_command_line("echo ' gridsize m4      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$disp3 -bj        '   >> control")
        call execute_command_line("echo '$extol 2.500      '   >> control")
        call execute_command_line("echo '$end              '   >> control")
  !     PW6B95-D3BJ
     elseif(func.eq.8)  then
        call execute_command_line("echo ' functional pw6b95'   >> control")
        call execute_command_line("echo ' gridsize m4      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$disp3 -bj        '   >> control")
        call execute_command_line("echo '$extol 2.500      '   >> control")
        call execute_command_line("echo '$end              '   >> control")
  !     BLYP-D3BJ
     elseif(func.eq.10)  then
        call execute_command_line("echo ' functional b-lyp'   >> control")
        call execute_command_line("echo ' gridsize m3      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$disp3 -bj        '   >> control")
        call execute_command_line("echo '$end              '   >> control")
  !     BP86-D3BJ
     elseif(func.eq.11)  then
        call execute_command_line("echo ' functional b-p'   >> control")
        call execute_command_line("echo ' gridsize m3      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$disp3 -bj        '   >> control")
        call execute_command_line("echo '$end              '   >> control")
  !     TPSS-D3BJ
     elseif(func.eq.12) then
        call execute_command_line("echo ' functional tpss'   >> control")
        call execute_command_line("echo ' gridsize m3      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$disp3 -bj        '   >> control")
        call execute_command_line("echo '$end              '   >> control")
  !     BH-LYP-D3BJ
     elseif(func.eq.15) then
        call execute_command_line("echo ' functional bh-lyp'   >> control")
        call execute_command_line("echo ' gridsize m4      '   >> control")
        call execute_command_line("echo '$pop              '   >> control")
        call execute_command_line("echo '$rij              '   >> control")
        call execute_command_line("echo '$disp3 -bj        '   >> control")
        call execute_command_line("echo '$end              '   >> control")
     endif
  
  
  end subroutine initm

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     This is to update the coord file with the actual coordinates of the
  !     selected fragment. Necessary for the define input of any subsequent
  !     trajectory.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine wrcoord(n,xyz,iat)

     integer  :: i
     integer  :: iat(n)
     integer  :: n
     integer  :: io_coord

     real(wp) ::  xyz(3,n)
  
     write(*,*)'updating coord file...'
  
     open(file='coord',newunit=io_coord,status='replace')
     write(io_coord,'(a)')'$coord'
     do i=1,n
        write(io_coord,'(3F20.14,6x,a2)')&
        &     xyz(1,i),xyz(2,i),xyz(3,i),toSymbol(iat(i))
     enddo
     write(io_coord,'(a)')'$end'
     close(unit=io_coord)
  
  end subroutine wrcoord


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !copy stuff for turbomole
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine copytm(it)
     integer  :: it
     character(len=80) :: fname
  
     if(it.ge.10000)error stop 'error 1 inside copytm'
  
     call execute_command_line('cp coord coord.original')
  
     if(it.ge.1000)then
        write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i4)')it
        call execute_command_line(fname)
        write(fname,'(''cp alpha   TMPQCXMS/TMP.'',i4)')it
        call execute_command_line(fname)
        write(fname,'(''cp beta    TMPQCXMS/TMP.'',i4)')it
        call execute_command_line(fname)
        write(fname,'(''cp control TMPQCXMS/TMP.'',i4)')it
        call execute_command_line(fname)
        write(fname,'(''cp basis   TMPQCXMS/TMP.'',i4)')it
        call execute_command_line(fname)
        write(fname,'(''cp coord   TMPQCXMS/TMP.'',i4)')it
        call execute_command_line(fname)
        write(fname,'(''mv coord.original TMPQCXMS/TMP.'',i4)')it
        call execute_command_line(fname)
        return
     endif
  
     if(it.ge.100)then
        write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i3)')it
        call execute_command_line(fname)
        write(fname,'(''cp alpha   TMPQCXMS/TMP.'',i3)')it
        call execute_command_line(fname)
        write(fname,'(''cp beta    TMPQCXMS/TMP.'',i3)')it
        call execute_command_line(fname)
        write(fname,'(''cp control TMPQCXMS/TMP.'',i3)')it
        call execute_command_line(fname)
        write(fname,'(''cp basis   TMPQCXMS/TMP.'',i3)')it
        call execute_command_line(fname)
        write(fname,'(''cp coord   TMPQCXMS/TMP.'',i3)')it
        call execute_command_line(fname)
        write(fname,'(''mv coord.original TMPQCXMS/TMP.'',i3)')it
        call execute_command_line(fname)
        return
     endif
  
     if(it.ge.10)then
        write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i2)')it
        call execute_command_line(fname)
        write(fname,'(''cp alpha   TMPQCXMS/TMP.'',i2)')it
        call execute_command_line(fname)
        write(fname,'(''cp beta    TMPQCXMS/TMP.'',i2)')it
        call execute_command_line(fname)
        write(fname,'(''cp control TMPQCXMS/TMP.'',i2)')it
        call execute_command_line(fname)
        write(fname,'(''cp basis   TMPQCXMS/TMP.'',i2)')it
        call execute_command_line(fname)
        write(fname,'(''cp coord   TMPQCXMS/TMP.'',i2)')it
        call execute_command_line(fname)
        write(fname,'(''mv coord.original TMPQCXMS/TMP.'',i2)')it
        call execute_command_line(fname)
        return
     endif
  
     if(it.ge.0)then
        write(fname,'(''cp qcxms.in TMPQCXMS/TMP.'',i1)')it
        call execute_command_line(fname)
        write(fname,'(''cp alpha   TMPQCXMS/TMP.'',i1)')it
        call execute_command_line(fname)
        write(fname,'(''cp beta    TMPQCXMS/TMP.'',i1)')it
        call execute_command_line(fname)
        write(fname,'(''cp control TMPQCXMS/TMP.'',i1)')it
        call execute_command_line(fname)
        write(fname,'(''cp basis   TMPQCXMS/TMP.'',i1)')it
        call execute_command_line(fname)
        write(fname,'(''cp coord   TMPQCXMS/TMP.'',i1)')it
        call execute_command_line(fname)
        write(fname,'(''mv coord.original TMPQCXMS/TMP.'',i1)')it
        call execute_command_line(fname)
        return
     endif
  
     error stop 'error 2 inside copytm'
  
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! for caclulating IPs with TMOL (more copy needed)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine copytm_ip(it)
     implicit none
  
     integer  :: it
     character(len=80) :: fname
  
     if(it.ge.10000)error stop 'error 1 inside copytm'
  
     call execute_command_line('cp coord coord.original')
  
     if(it.ge.1000)then
        write(fname,'(''mv coord.original TMPQCXMS/TMP.'',i4)')it
        call execute_command_line(fname)
        return
     endif
  
     if(it.ge.100)then
        write(fname,'(''mv coord.original TMPQCXMS/TMP.'',i3)')it
        call execute_command_line(fname)
        return
     endif
  
     if(it.ge.10)then
        write(fname,'(''mv coord.original TMPQCXMS/TMP.'',i2)')it
        call execute_command_line(fname)
        return
     endif
  
     if(it.ge.0)then
        write(fname,'(''mv coord.original TMPQCXMS/TMP.'',i1)')it
        call execute_command_line(fname)
        return
     endif
  
     error stop 'error 2 inside copytm'
  
  end subroutine


end module qcxms_use_turbomole
