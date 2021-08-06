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

        integer :: nat
        integer :: chrg
        integer :: iat(nat)
        integer :: i,ich,j,spin,idum(100),specbas(4)

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
           if(idum(specbas(i)).gt.0) wrb=.true.
        enddo
     
        if(spin.eq.0) call getspin(nat,iat,chrg,spin)
     
        basis(1)='SV'
        basis(2)='SV'
        basis(3)='SV(P)'
        basis(4)='SVP'
        basis(5)='TZVP'
        basis(6)='ma-def2-SVP'
        basis(7)='ma-def2-TZVP'
        basis(8)='ma-def2-TZVPP'
        basis(9)='def2-SV(P)'
        basis(10)='def2-SVP'
        basis(11)='def2-TZVP'
        basis(12)='def2-QZVP'
        basis(13)='QZVP'
        basis(14)='ma-def2-QZVP'
        ich=87
        open(unit=ich,file='ORCA.INPUT')
     ! hybrid vs other funcs.... nat is number of atoms
        if(func.le.4.and.nat.lt.60.and.noconv.eq. .false.)then
           write(ich,'(''! CONV SMALLPRINT NORI NOSOSCF'')')
        elseif (func .eq. 7.and.nat.lt.60.and.noconv.eq. .false.) then
           write(ich,'(''! CONV SMALLPRINT NORI NOSOSCF'')')
        elseif (func .eq. 8.and.nat.lt.60.and. noconv .eq. .false.) then
           write(ich,'(''! CONV SMALLPRINT NORI NOSOSCF'')')
        elseif (func .eq. 9.and.nat.lt.60 .and. noconv .eq. .false.) then
           write(ich,'(''! CONV SMALLPRINT NORI NOSOSCF'')')
        elseif(func.le.4.and.noconv.eq. .True.)then
           write(ich,'(''! DIRECT SMALLPRINT NORI NOSOSCF'')')
        elseif (func .eq. 7.and.noconv.eq. .True.) then
           write(ich,'(''! DIRECT SMALLPRINT NORI NOSOSCF'')')
        elseif (func .eq. 8.and.noconv.eq. .True.) then
           write(ich,'(''! DIRECT SMALLPRINT NORI NOSOSCF'')')
        elseif (func .eq. 9.and.noconv.eq. .True.) then
           write(ich,'(''! DIRECT SMALLPRINT NORI NOSOSCF'')')
        else
     ! RI
           write(ich,'(''!  DEF2/J SMALLPRINT NOSOSCF'')')
        endif
     
        write(ich,'(''! NOFINALGRID NOMAYER'')')
     
        write(ich,'(''! UHF'')')
     ! IF PBEH3C
        if(func.ne.14) then
           write(ich,'(''! '',a)')trim(basis(bas))
        elseif(func.eq.14) then
           write(ich,'(''! def2-mSVP'')')
        endif
     ! IF ECP true
        if (ecp .eq. .true. .and. func.ne.14)then
           write(ich,'(''! ecp('',a,'')'')')trim(basis(bas))
        elseif(ecp .eq. .true. .and. func .eq. 14)then
           write(ich,'(''! ecp(def2-mSVP)'')')
        endif
        write(ich,'(''! GRID'',i1)')grid
     ! func 0,1,2 are pbe hybrid family
     ! func 3 is lda (commented out)
        if ( func == 0  ) write(ich,'(''! PBE0 D3'')')
        if ( func == 1  ) write(ich,'(''! D3BJ'')')
        if ( func == 2  ) write(ich,'(''! D3BJ'')')
     !  if ( func == 3  ) write(ich,'(''! D3BJ'')')
        if ( func == 4  ) write(ich,'(''! M062X D3zero'')')
        if ( func == 5  ) write(ich,'(''! PBE D3BJ'')')
        if ( func == 6  ) write(ich,'(''! B97-D3 '')')
        if ( func == 7  ) write(ich,'(''! B3LYP D3BJ'')')
        if ( func == 8  ) write(ich,'(''! PW6B95 D3BJ'')')
        if ( func == 9  ) write(ich,'(''! B3PW91 D3BJ'')')
        if ( func == 10 ) write(ich,'(''! BLYP D3BJ '')')
        if ( func == 11 ) write(ich,'(''! BP86 D3BJ'')')
        if ( func == 12 ) write(ich,'(''! TPSS D3BJ'')')
        if ( func == 13 ) write(ich,'(''! REVPBE D3BJ'')')
        if ( func == 14 ) write(ich,'(''! PBEH-3c'')')
        if ( func == 15 ) write(ich,'(''! BHLYP D3BJ'')')
        write(ich,'(''%output'')')
        write(ich,'('' Print[ P_AtPopMO_L] 1  end'')')
        write(ich,'(''%output'')')
        write(ich,'('' Print[ P_AtPopMO_M] 1  end'')')
        write(ich,'(''%method'')')
        if(grad)then
           write(ich,'('' runtyp gradient'')')
        else
           write(ich,'('' runtyp energy  '')')
        endif
     
        if(func.eq.1.or.func.eq.2)then
           write(ich,'(''   method dfgto'')')
           write(ich,'(''   ldaopt C_PWLDA'')')
        endif
     
     !      if(func.eq.3)then
     !      write(ich,'(''   exchange    X_SLATER'')')
     !      write(ich,'(''   correlation C_VWN3'')')
     !      write(ich,'(''   scalHFX = 0.500'')')
     !      write(ich,'(''   scalDFX = 0.500'')')
     !      endif
     
     
     ! pbe12, pbe38,pbe0
        if(func.eq.1)then
           write(ich,'(''   exchange    X_PBE'')')
           write(ich,'(''   correlation C_PBE'')')
           write(ich,'(''   scalHFX = 0.375'')')
           write(ich,'(''   scalDFX = 0.625'')')
        endif
        if(func.eq.2)then
           write(ich,'(''   exchange    X_PBE'')')
           write(ich,'(''   correlation C_PBE'')')
           write(ich,'(''   scalHFX = 0.500'')')
           write(ich,'(''   scalDFX = 0.500'')')
        endif
     !      if(func.eq.0)then
     !      write(ich,'(''   exchange    X_PBE'')')
     !      write(ich,'(''   correlation C_PBE'')')
     !      write(ich,'(''   scalHFX = 0.250'')')
     !      write(ich,'(''   scalDFX = 0.750'')')
     !      endif
     
        write(ich,'(''end'')')
     
        write(ich,'(''%elprop'')')
        write(ich,'('' dipole false'')')
        write(ich,'(''end'')')
     
        write(ich,'(''%scf'')')
        if(etemp.gt.10.0)&
        &write(ich,'('' SmearTemp '',F7.0)') etemp
        write(ich,'('' maxcore   '',i6  )') qcmem
        write(ich,'('' MaxIter  400''     )')
        if(spin.gt.1)write(ich,'('' LShift  0.30''     )')
        write(ich,'(''end'')')
     
        if(bas.eq.2.and.wrb)then
           write(ich,'(''%basis    '')')
           if(idum(9).gt.0)then
              write(ich,'('' addGTO 9  '')')
              write(ich,'('' D 1      '')')
              write(ich,'('' 1 1.2 1.0'')')
              write(ich,'('' end'')')
           endif
           if(idum(8).gt.0)then
              write(ich,'('' addGTO 8  '')')
              write(ich,'('' D 1      '')')
              write(ich,'('' 1 1.0 1.0'')')
              write(ich,'('' end'')')
           endif
           if(idum(7).gt.0)then
              write(ich,'('' addGTO 7  '')')
              write(ich,'('' D 1      '')')
              write(ich,'('' 1 0.8 1.0'')')
              write(ich,'('' end'')')
           endif
           if(idum(14).gt.0)then
              write(ich,'('' addGTO 14 '')')
              write(ich,'('' D 1      '')')
              write(ich,'('' 1 0.35 1.0'')')
              write(ich,'('' end'')')
           endif
           write(ich,'(''end'')')
        endif
     
        write(ich,'(''* xyz '',2i3)')chrg,spin
        do i=1,nat
           write(ich,300)toSymbol(iat(i)),xyz(1,i)*autoaa,xyz(2,i)*autoaa,xyz(3,i)*autoaa
        enddo
     300 format(a2,3(f22.12))
        write(ich,'(''*'')')
     
        close(ich)
     
     end subroutine orcaout
     
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     
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
        read(io_orca,'(a)')line
        read(io_orca,'(a)')line
        read(io_orca,'(a)')line
        read(io_orca,*) i
        read(io_orca,'(a)')line
        read(io_orca,'(a)')line
        read(io_orca,'(a)')line
        read(io_orca,*) edum
        read(io_orca,'(a)')line
        read(io_orca,'(a)')line
        read(io_orca,'(a)')line
        do j=1,nat
           read(io_orca,*)g(1,j)
           read(io_orca,*)g(2,j)
           read(io_orca,*)g(3,j)
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
   external :: system
   
   if(line.ge.10000)stop 'error 1 inside copyorc'
   
   if(line.ge.1000)then
      write(command,'(''cp qcxms.in TMPQCXMS/TMP.'',i4)')line
      call system(command)
      write(command,'(''cp coord TMPQCXMS/TMP.'',i4)')line
      call system(command)
      return
   endif
   
   if(line.ge.100)then
      write(command,'(''cp qcxms.in TMPQCXMS/TMP.'',i3)')line
      call system(command)
      write(command,'(''cp coord TMPQCXMS/TMP.'',i3)')line
      call system(command)
      return
   endif
   
   if(line.ge.10)then
      write(command,'(''cp qcxms.in TMPQCXMS/TMP.'',i2)')line
      call system(command)
      write(command,'(''cp coord TMPQCXMS/TMP.'',i2)')line
      call system(command)
      return
   endif
   
   if(line.ge.0)then
      write(command,'(''cp qcxms.in TMPQCXMS/TMP.'',i1)')line
      call system(command)
      write(command,'(''cp coord TMPQCXMS/TMP.'',i1)')line
      call system(command)
      return
   endif
   
   stop 'error 2 inside copyorc'
   
   end subroutine


end module qcxms_use_orca
