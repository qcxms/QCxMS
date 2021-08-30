module qcxms_use_mndo
   use common1
   use readcommon
   use qcxms_utility, only: getspin
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   implicit none

   contains 

! read MNDO99 E and gradient
  subroutine mndograd(fname,nat,g,chrg,spin,edum)
  
     integer  :: nat,i,j,nn,idum1,idum2
     integer  :: io_mndo
     real(wp) :: g(3,nat),xx(20),edum,chrg(nat),spin(nat)
  
     character(len=*  ) :: fname
  
     chrg = 1.d-4
     spin = 0
     open(unit=io_mndo,file=fname)
rd1: do
        read(io_mndo,'(a)',iostat=iocheck)line
        if (iocheck > 0) stop 'Failed to read MNDO input file'
        if (iocheck < 0) exit rd1 ! EOF
  
        if (index(line,'SCF TOTAL ENERGY') /= 0) then
           call readl(line,xx,nn)
           edum =xx(1) * evtoau 
        endif

        if(index(line,'NET ATOMIC AND ORBITAL SPIN').ne.0)then
           do i = 1,3
               read(io_mndo,'(a)',iostat=iocheck)line
               if(iocheck < 0) exit rd1
           enddo

           do j=1,nat
              read(io_mndo,'(a)',iostat=iocheck)line
              if(iocheck < 0) exit rd1
              call readl(line,xx,nn)
              spin(j)=xx(2)
           enddo
        endif
  
        if(index(line,'NET ATOMIC CHARGES AND').ne.0)then
           do i = 1,3
               read(io_mndo,'(a)',iostat=iocheck)line
               if(iocheck < 0) exit rd1
           enddo

           do j=1,nat
              read(io_mndo,'(a)',iostat=iocheck)line
              if(iocheck < 0) exit rd1
              call readl(line,xx,nn)
              chrg(j)=xx(2)
           enddo
        endif

     enddo rd1
     close(io_mndo)
  
     open( unit=io_mndo,file='fort.15' )
rd2: do 
        read( io_mndo,'(a)', iostat=iocheck ) line
        if ( iocheck > 0 ) stop 'could not open fort.15!'
        if ( iocheck < 0 ) exit rd2

        if ( index(line,' CARTESIAN GRADIENT') /= 0 ) then
           do j = 1, nat
              read ( io_mndo,*, iostat=iocheck) idum1, idum2, g(1:3,j)
              if ( iocheck < 0 ) exit rd2

              g(1:3,j)=g(1:3,j)/(autokcal * aatoau) !   627.509541d0 /0.529177260d0)
           enddo
        endif
     enddo rd2
     close(io_mndo)

  end subroutine mndograd
  
  
  ! write MNDO99 input
  subroutine mndoout(nat,xyz,ic,chrg,isp,etemp,iconv,idiis)
  
     integer  :: chrg,nat,ic(*),isp,iconv,idiis
     integer  :: i,m,iuhf
     integer  :: io_input
  
     real(wp) :: xyz(3,nat),etemp
  
  
     if(isp.eq.0) call getspin(nat,ic,chrg,isp)
  
     if(isp.eq.1)isp=0
  
     open( unit=io_input, file='inp', status='replace')
  
     write(io_input,'(''nsav15=3 idiis='',i2,'' jop=-2 ipubo=1 nprint=-1 iscf='',i1, &
     &         '' kitscf=200 iplscf='',i1,'' iop='',i2,'' +'')')&
     &    idiis,iconv,iconv,-ihamilt
  
     if ( etemp < 0.1 ) then
        write(io_input,'(''ifast=2 igeom=1 ktrial=11 iuhf=1 kharge='',&
        &       i1,'' iform=1 imult='',i1)')&
        &chrg,isp
  
     else
        if ( etemp >= 10000 ) &
        &  write(io_input,'(''ifast=2 igeom=1 ktrial=11 iuhf=5 kharge='',&
        &       i1,'' iform=1 imult='',i1,'' ifermi='',i5)')&
        &       chrg,isp,int(etemp)
  
        if ( etemp >= 1000 .and. etemp < 10000) &
        &  write(io_input,'(''ifast=2 igeom=1 ktrial=11 iuhf=5 kharge='',&
        &       i1,'' iform=1 imult='',i1,'' ifermi='',i4)')&
        &  chrg,isp,int(etemp)
  
        if ( etemp >= 100 .and. etemp < 1000 ) &
        &  write(io_input,'(''ifast=2 igeom=1 ktrial=11 iuhf=5 kharge='',&
        &       i1,'' iform=1 imult='',i1,'' ifermi='',i3)')&
        &  chrg,isp,int(etemp)
     endif
  
     write(io_input,*)
     write(io_input,*)
  
     m=1
     do i=1, nat
        write(io_input,'(i4,3(F22.12,i3))') &
        &   ic(i),autoaa * xyz(1,i),m,autoaa*xyz(2,i),m,autoaa*xyz(3,i),m
     enddo
  
     write(io_input,*)'0'
  
     close(io_input)
  
  end subroutine mndoout

end module qcxms_use_mndo
