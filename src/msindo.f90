module qcxms_use_msindo
   use newcommon
   use readcommon
   use qcxms_utility, only: getspin
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only: toSymbol 
   implicit none

  contains

  ! get gradient with actual coords
  subroutine getmsindograd(first,nat,ic,chrg,spin,xyz,etemp,nav,&
  &                        g,edum,achrg,aspin)
  
     integer :: nat,i,j,k,ic(nat),nn,cyc
     integer :: chrg,spin,nav
     integer :: io_msindo, io_tmp
  
     real(wp) :: g(3,nat),xx(20),edum,etemp
     real(wp) :: xyz(3,nat)
     real(wp) :: achrg(nat)
     real(wp) :: aspin(nat)
  
     logical :: first
  
     if(spin.eq.0) call getspin(nat,ic,chrg,spin)
  
     cyc=200
     if(nav.gt.10)cyc=999
  
     edum=0
     open(unit=io_msindo,file='inp')
  !     write(io_msindo,'(''NDDO'')')
  !     write(io_msindo,'(''       SHD2(1)           0.60'')')
  !     write(io_msindo,'(''       SHD1(8)           0.50'')')
  !     write(io_msindo,'(''       SHD2(1)           0.60'')')
  !     write(io_msindo,'(''       SHD1(8)           0.55'')')
     write(io_msindo,'(''***'')')
     write(io_msindo,'(''QCxMS generated'')')
     write(io_msindo,'(''UHF NGIV=8 MULTIP='',i1,'' CHARGE='',i1)')spin,chrg
     write(io_msindo,'(''NDDO'')')
     write(io_msindo,'(''MAXCYC='',i3)')cyc
     write(io_msindo,'(''NAV='',i2)')nav
  
     if(etemp.gt.1.)write(io_msindo,'(''FERMI_SMEARING T_SMEAR='',f6.0)')etemp
  
     write(io_msindo,'(''DELEN 1.D-8'')')
     write(io_msindo,'(''NOSYM'')')
  
     if(first)then
        write(io_msindo,'(''IDEN=8'')')
     else
        write(io_msindo,'(''IDEN=10'')')
     endif
     write(io_msindo,'(''OPTDIAMEM'')')
     write(io_msindo,'(''PRINTOPTS=MULAP'')')
     write(io_msindo,'(''CARTOPT GRADONLY'')')
     write(io_msindo,'(''CARTES'')')
     write(io_msindo,'('':End'')')
     do i=1,nat
        write(io_msindo,'(4x,a2,3F20.12)')&
        & toSymbol(ic(i)),&
        & xyz(1,i)*autoaa, &
        & xyz(2,i)*autoaa, &
        & xyz(3,i)*autoaa
     enddo
     write(io_msindo,'('':End'')')
     write(io_msindo,'('':End'')')
     write(io_msindo,'(''END'')')
     close(io_msindo)
  
     calls = calls + 1
     call execute_command_line('/usr/local/bin/msindo < inp > job.last')
  
  
     open(unit=io_tmp,file='tmpcoord')
     write(io_tmp,'(''$coord'')')
     do j=1,nat
        write(io_tmp,'(3F20.12,10x,a2)')&
        &xyz(1,j),xyz(2,j),xyz(3,j),toSymbol(ic(j))
     enddo
     write(io_tmp,'(''$end'')')
     close(io_tmp)
  
     open(unit=io_tmp,file='job.last')
  rd: do
        read(io_tmp,'(a)',iostat=iocheck)line
        if (iocheck < 0 ) exit rd
        if (iocheck > 0 ) stop 'Error in tmp file '
  
        if(index(line,'TOTAL ENERGY = ').ne.0)then
           call readl(line,xx,nn)
           edum =xx(1)
        endif
  
        if(index(line,'GROSS ATOMIC POPULATION').ne.0)then
           do k = 1, 3 
              read(io_tmp,'(a)',iostat=iocheck)line
              if (iocheck < 0 ) exit rd
           enddo
           do j=1,nat
              read(io_tmp,'(a)')line
              call readl(line,xx,nn)
              achrg(j)=xx(2)
           enddo
        endif
  
        if(index(line,'SPIN DENSITY').ne.0)then
           do k = 1, 3
              read(io_tmp,'(a)',iostat=iocheck)line
              if (iocheck < 0 ) exit rd
           enddo
           do j=1,nat
              read(io_tmp,'(a)')line
              call readl(line,xx,nn)
              aspin(j)=xx(2)
           enddo
        endif
  
        if(index(line,'FIRST DERIVATIVES ').ne.0)then
           read(io_tmp,'(a)',iostat=iocheck)line
           if (iocheck < 0 ) exit rd
           do j=1,nat
              read(io_tmp,'(a)')line
              call readl(line,xx,nn)
              g(1,j)=xx(2)
              g(2,j)=xx(3)
              g(3,j)=xx(4)
           enddo
        endif
     enddo rd
     close(io_tmp)
  
  end subroutine getmsindograd

end module qcxms_use_msindo
