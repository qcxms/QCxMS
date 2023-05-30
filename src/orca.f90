module qcxms_use_orca 
   use common1
   use readcommon
   use qcxms_utility, only: getspin
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only: toSymbol 
   implicit none

   contains

     subroutine orcaout(nat,xyz,iat,chrg,spin,etemp,grad,ECP)

        integer ::  nat
        integer ::  chrg
        integer ::  iat(nat)
        integer ::  i,j,spin,idum(100),specbas(4)
        integer ::  io_orca

        real(wp) :: xyz (3,nat)
        real(wp) :: etemp

        character(len=30) :: basis(14)

        logical :: grad,wrb
        logical :: ECP

        data specbas /7,8,9,14/
     
     
        ECP = .False.
     ! stuff for SVx output to ORCA
        wrb=.false.
        idum = 0
        do i=1,nat
           idum(iat(i))=idum(iat(i))+1
        enddo
        do i=1,4
           if(idum(specbas(i)) > 0) wrb=.true.
        enddo
     
        if(spin == 0) call getspin(nat,iat,chrg,spin)
     
        basis(1)  = 'SV'
        basis(2)  = 'SV'
        basis(3)  = 'SV(P)'
        basis(4)  = 'SVP'
        basis(5)  = 'TZVP'
        basis(6)  = 'ma-def2-SVP'
        basis(7)  = 'ma-def2-TZVP'
        basis(8)  = 'ma-def2-TZVPP'
        basis(9)  = 'def2-SV(P)'
        basis(10) = 'def2-SVP'
        basis(11) = 'def2-TZVP'
        basis(12) = 'def2-QZVP'
        basis(13) = 'QZVP'
        basis(14) = 'ma-def2-QZVP'

        ! open and write to ORCA input file
        open(file='ORCA.INPUT', newunit=io_orca)
        
        ! settings for ORCA 4, for ORCA 5 default settings should be fine
        if ( orca_version == 4 ) then 
     ! hybrid vs other funcs.... nat is number of atoms
        if ( func <= 4 .and. nat < 60 .and. .not. noconv ) then
         !if ( No_eTemp ) then
           write(io_orca,'(''! CONV SMALLPRINT NOSOSCF RIJK DEF2/JK'')')
         !else 
         !  write(io_orca,'(''! CONV SMALLPRINT NORI NOSOSCF'')')
         !endif

        elseif ( func  ==  7 .and. nat < 60 .and. .not.noconv ) then
           write(io_orca,'(''! CONV SMALLPRINT NORI NOSOSCF'')')

        elseif ( func  ==  8 .and. nat < 60 .and. .not.noconv ) then
           write(io_orca,'(''! CONV SMALLPRINT NORI NOSOSCF'')')

        elseif ( func  ==  9 .and. nat < 60 .and. .not.noconv ) then
           write(io_orca,'(''! CONV SMALLPRINT NORI NOSOSCF'')')

        elseif ( func <= 4 .and. noconv)then
           write(io_orca,'(''! DIRECT SMALLPRINT NORI NOSOSCF'')')

        elseif ( func  ==  7 .and. noconv) then
           write(io_orca,'(''! DIRECT SMALLPRINT NORI NOSOSCF'')')

        elseif ( func  ==  8 .and. noconv) then
           write(io_orca,'(''! DIRECT SMALLPRINT NORI NOSOSCF'')')

        elseif ( func  ==  9 .and. noconv) then
           write(io_orca,'(''! DIRECT SMALLPRINT NORI NOSOSCF'')')

        ! RI
        else
           write(io_orca,'(''!  DEF2/J SMALLPRINT NOSOSCF'')')
        endif

        else
           write(io_orca,'(''!  DEF2/J SMALLPRINT'')')
        end if
     
        ! Set mayer and finalgrid
        if ( orca_version == 4 ) write(io_orca,'(''! NOFINALGRID NOMAYER'')')
        if ( orca_version == 5 ) then
           write(io_orca,'(''! NOMAYER'')')
       endif
     
        write(io_orca,'(''! UHF'')')

        ! IF PBEh-3c
        if (func /= 14) then
           write(io_orca,'(''! '',a)')trim(basis(bas))
        elseif(func == 14) then
           write(io_orca,'(''! def2-mSVP'')')
        endif

        ! IF ECP true
        if (ecp .and. func /= 14)then
           write(io_orca,'(''! ecp('',a,'')'')')trim(basis(bas))
        elseif(ecp .and. func  ==  14)then
           write(io_orca,'(''! ecp(def2-mSVP)'')')
        endif

        ! Set Grid Mashes
        if ( orca_version == 4 ) write(io_orca,'(''! GRID'',i1)')   grid_orca 
        if ( orca_version == 5 ) write(io_orca,'(''! DEFGRID'',i1)')grid_orca 
        

     ! func 0,1,2 are pbe hybrid family
     ! func 3 is lda (commented out)
        if ( func == 0  ) write(io_orca,'(''! PBE0 D4'')')
        if ( func == 1  ) write(io_orca,'(''! D4'')')
        if ( func == 2  ) write(io_orca,'(''! D4'')')
        !if ( func == 1  ) write(io_orca,'(''! D3BJ'')')
        !if ( func == 2  ) write(io_orca,'(''! D3BJ'')')
     !  if ( func == 3  ) write(io_orca,'(''! D3BJ'')')
        if ( func == 4  ) write(io_orca,'(''! M062X D3zero'')')
        if ( func == 5  ) write(io_orca,'(''! PBE D4'')')
        !if ( func == 5  ) write(io_orca,'(''! PBE D3BJ'')')
        if ( func == 6  ) write(io_orca,'(''! B97-D3 '')')
        if ( func == 7  ) write(io_orca,'(''! B3LYP D4'')')
        !if ( func == 7  ) write(io_orca,'(''! B3LYP D3BJ'')')
        if ( func == 8  ) write(io_orca,'(''! PW6B95 D3BJ'')')
        if ( func == 9  ) write(io_orca,'(''! B3PW91 D3BJ'')')
        if ( func == 10 ) write(io_orca,'(''! BLYP D3BJ '')')
        if ( func == 11 ) write(io_orca,'(''! BP86 D3BJ'')')
        if ( func == 12 ) write(io_orca,'(''! TPSS D3BJ'')')
        if ( func == 13 ) write(io_orca,'(''! REVPBE D3BJ'')')
        if ( func == 14 ) write(io_orca,'(''! PBEH-3c'')')
        if ( func == 15 ) write(io_orca,'(''! BHLYP D3BJ'')')

        write(io_orca,'(''%output'')')
        write(io_orca,'('' Print[ P_AtPopMO_L] 1  end'')')
        write(io_orca,'(''%output'')')
        write(io_orca,'('' Print[ P_AtPopMO_M] 1  end'')')

        ! Set multicores
        if ( nproc_orca > 0 ) then
          write(io_orca,'(''%pal'')')
          write(io_orca,'('' nprocs '',i2)') nproc_orca
          write(io_orca,'(''end'')')
        endif

        write(io_orca,'(''%method'')')

        if(grad)then
           write(io_orca,'('' runtyp gradient'')')
        else
           write(io_orca,'('' runtyp energy  '')')
        endif
     
        if(func == 1.or.func == 2)then
           write(io_orca,'(''   method dfgto'')')
           write(io_orca,'(''   ldaopt C_PWLDA'')')
        endif
     
     !      if(func == 3)then
     !      write(io_orca,'(''   exchange    X_SLATER'')')
     !      write(io_orca,'(''   correlation C_VWN3'')')
     !      write(io_orca,'(''   scalHFX = 0.500'')')
     !      write(io_orca,'(''   scalDFX = 0.500'')')
     !      endif
     
     
     ! pbe12, pbe38,pbe0
        if(func == 1)then
           write(io_orca,'(''   exchange    X_PBE'')')
           write(io_orca,'(''   correlation C_PBE'')')
           write(io_orca,'(''   scalHFX = 0.375'')')
           write(io_orca,'(''   scalDFX = 0.625'')')
        endif
        if(func == 2)then
           write(io_orca,'(''   exchange    X_PBE'')')
           write(io_orca,'(''   correlation C_PBE'')')
           write(io_orca,'(''   scalHFX = 0.500'')')
           write(io_orca,'(''   scalDFX = 0.500'')')
        endif
     !      if(func == 0)then
     !      write(io_orca,'(''   exchange    X_PBE'')')
     !      write(io_orca,'(''   correlation C_PBE'')')
     !      write(io_orca,'(''   scalHFX = 0.250'')')
     !      write(io_orca,'(''   scalDFX = 0.750'')')
     !      endif
     
        write(io_orca,'(''end'')')
     
        write(io_orca,'(''%elprop'')')
        write(io_orca,'('' dipole false'')')
        write(io_orca,'(''end'')')
     
        write(io_orca,'(''%scf'')')
        !> test if fermi-smearing is even important
        if(etemp > 10.0 .and. .not. No_eTemp ) write(io_orca,'('' SmearTemp '',F7.0)') etemp
        write(io_orca,'('' maxcore   '',i6  )') qcmem
        write(io_orca,'('' MaxIter  400''     )')

        if(spin > 1)write(io_orca,'('' LShift  0.30''     )')
        write(io_orca,'(''end'')')
     
        if(bas == 2.and.wrb)then
           write(io_orca,'(''%basis    '')')
           if(idum(9) > 0)then
              write(io_orca,'('' addGTO 9  '')')
              write(io_orca,'('' D 1      '')')
              write(io_orca,'('' 1 1.2 1.0'')')
              write(io_orca,'('' end'')')
           endif
           if(idum(8) > 0)then
              write(io_orca,'('' addGTO 8  '')')
              write(io_orca,'('' D 1      '')')
              write(io_orca,'('' 1 1.0 1.0'')')
              write(io_orca,'('' end'')')
           endif
           if(idum(7) > 0)then
              write(io_orca,'('' addGTO 7  '')')
              write(io_orca,'('' D 1      '')')
              write(io_orca,'('' 1 0.8 1.0'')')
              write(io_orca,'('' end'')')
           endif
           if(idum(14) > 0)then
              write(io_orca,'('' addGTO 14 '')')
              write(io_orca,'('' D 1      '')')
              write(io_orca,'('' 1 0.35 1.0'')')
              write(io_orca,'('' end'')')
           endif
           write(io_orca,'(''end'')')
        endif
     
        write(io_orca,'(''* xyz '',2i3)')chrg,spin
        do i=1,nat
           write(io_orca,300)toSymbol(iat(i)),xyz(1,i)*autoaa,xyz(2,i)*autoaa,xyz(3,i)*autoaa
        enddo
     300 format(a2,3(f22.12))
        write(io_orca,'(''*'')')
     
        close(io_orca)
     
     end subroutine orcaout
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     subroutine rdorcagrad(fname,nat,g,achrg,aspin,edum)

        integer :: nat
        integer :: i,j,nn
        integer :: io_orca

        real(wp) :: g(3,nat),edum
        real(wp) :: achrg(nat),aspin(nat)
        real(wp) :: xx(20)
     
        character(len=*)  :: fname

        logical :: ex
     
        open(file=fname, newunit=io_orca, status='old')
        do
           read ( unit=io_orca, FMT='(a)', iostat=iocheck ) line
           if ( iocheck > 0 ) stop 'error in rdorcagrad!'
           if ( iocheck < 0 ) exit !EOF

           if ( index(line,'FINAL SINGLE POINT ENERGY') /= 0 ) then
              call readl(line,xx,nn)
              edum=xx(1)
           endif

           if ( index(line,'LOEWDIN ATOMIC CHARGES AND SPIN') /= 0 ) then
              read(unit=io_orca, FMT='(a)') line

              do i=1, nat
                 read(unit=io_orca,FMT='(a)') line
                 call readl(line,xx,nn)
                 achrg(i) = xx(2)
                 if ( nn == 3 ) aspin(i) = xx(3)
              enddo

           endif

        enddo 
        close(io_orca)
     
        inquire(file='ORCA.INPUT.engrad',exist=ex)
        if(.not.ex) return
     
        open(file='ORCA.INPUT.engrad', newunit=io_orca, status='old' )
        read(io_orca,'(a)')line !#
        read(io_orca,'(a)')line !#
        read(io_orca,'(a)')line !#
        read(io_orca,*) i       !no. of atms
        read(io_orca,'(a)')line !#
        read(io_orca,'(a)')line !#
        read(io_orca,'(a)')line !#
        read(io_orca,*) edum    !The current total energy in Eh
        read(io_orca,'(a)')line !#
        read(io_orca,'(a)')line !#
        read(io_orca,'(a)')line !#
        do j=1,nat
           read(io_orca,*)g(1,j) ! The current gradient in Eh/bohr
           read(io_orca,*)g(2,j) ! The current gradient in Eh/bohr
           read(io_orca,*)g(3,j) !  The current gradient in Eh/bohr
        enddo
        close(io_orca)
     
     end subroutine rdorcagrad
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     subroutine rdorcaen(fname,energy)

        integer :: nn
        integer :: io_orca

        real(wp) :: xx(20)
        real(wp) :: energy
     
        character(len=*)  :: fname
     
        open(file=fname, newunit=io_orca, status='old')

        do
          read (unit=io_orca, FMT='(a)', iostat=iocheck) line
          if ( iocheck > 0 ) stop 'error in rdorcaen!'
          if ( iocheck < 0 ) exit ! EOF

          if (index(line,'FINAL SINGLE POINT ENERGY') /= 0 ) then
             call readl(line,xx,nn)
             energy=xx(1)
             exit
          endif

        enddo

        close(io_orca)
     
     end subroutine rdorcaen


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine copyorc(line)
   
   integer :: line
   character(len=80) :: command
   
   if(line >= 10000)stop 'Too many folders. Please reduce size'
   
   if(line >= 1000)then
      write(command,'(''cp qcxms.in TMPQCXMS/TMP.'',i4)')line
      call execute_command_line(command)
      return
   endif
   
   if(line >= 100)then
      write(command,'(''cp qcxms.in TMPQCXMS/TMP.'',i3)')line
      call execute_command_line(command)
      return
   endif
   
   if(line >= 10)then
      write(command,'(''cp qcxms.in TMPQCXMS/TMP.'',i2)')line
      call execute_command_line(command)
      return
   endif
   
   if(line >= 0)then
      write(command,'(''cp qcxms.in TMPQCXMS/TMP.'',i1)')line
      call execute_command_line(command)
      return
   endif
   
   stop 'error 2 inside copyorc'
   
   end subroutine


end module qcxms_use_orca
