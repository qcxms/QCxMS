module qcxms_analyse
  use common1
  use qcxms_iniqm, only: eqm
  use qcxms_info, only: qcstring
  use qcxms_utility, only: getspin
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_convert
  use xtb_mctc_symbols, only: toSymbol 
  implicit none

  contains


subroutine analyse(iprog,nuc,iat,iatf,axyz,list,nfrag,etemp,fragip, mchrg, &
      natf,ipok,icoll,isec,nometal,ECP)
    
  integer :: nuc  
  integer :: iprog 
  integer :: list(nuc)
  integer :: iat (nuc)
  integer :: icoll,isec
  integer :: natf(10)
  integer :: i,j,k
  integer :: nfrag
  integer :: progi,itry,useprog(4)  
  integer :: isave,jsave,ksave,gsave
  integer :: iatf(nuc,10)
  !integer :: idum(nuc,10)
  integer :: neutfspin,ionfspin !fragment spin (for metals etc.)
  integer :: spin_neut,spin_ion !save these spins 
  integer :: fiter !Number of spin iterations 
  !integer :: sp(3),sn(3),sn0,sp0
  integer :: nb,nel
  integer :: mchrg, lpchrg, dump_chrg
  integer :: io_xyz 
  
  real(wp) :: axyz(3,nuc)
  real(wp) :: fragip(nfrag,abs(mchrg)),etemp
  real(wp) :: xyzf(3,nuc,10)
  !real(wp) :: dum (3,nuc,10)
  real(wp) :: z,E_neut,E_ion,cema(3,10),rf(10*(10+1)/2)
  real(wp) :: t2,t1,w2,w1
  real(wp) :: gsen(3),gsep(3)
  real(wp) :: dsave
  real(wp) :: lowest_neut, lowest_ion 
  
  character(len=80) :: fname
  character(len=20) :: line, line2
  
  logical :: ipok 
  logical :: nometal, ECP
  logical :: metal = .false. !if fragment has metal
  logical :: ipcalc
  logical :: spec_calc = .false.
  
  ipok = .true.
  ! timings
  t1 = 0.0_wp
  t2 = 0.0_wp
  w1 = 0.0_wp
  w2 = 0.0_wp

  write(*,'('' computing average fragment structures ...'')')
  call avg_frag_struc(nuc,iat,iatf,axyz,list,nfrag, natf, xyzf)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !write fragments with average geometries      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> CID
  if(method == 3) then 
    do i=1,nfrag
      cema(1:3,i) = 0
      z           = 0
  
      if (icoll < 10) write(fname,'(i1,''.'',i1,''.'',i1,''.xyz'')') icoll, isec, i
      if (icoll >= 10) write(fname,'(i2,''.'',i1,''.'',i1,''.xyz'')') icoll, isec, i

      open (file = fname, newunit = io_xyz)
  
      write(io_xyz,*) natf(i)
      write(io_xyz,*)
  
      do j=1,natf(i)
        write(io_xyz,'(a2,5x,3F18.8)') toSymbol(iatf(j,i)), xyzf(1:3,j,i) * autoaa 
        cema(1:3,i) = cema(1:3,i) + xyzf(1:3,j,i) * iatf(j,i)
        z = z + iatf(j,i)
      enddo

      close(io_xyz)
      cema(1:3,i) = cema(1:3,i) / z
     enddo   

  !> EI/DEA
  else
    do i=1,nfrag
      cema(1:3,i) = 0
      z           = 0
  
      write(fname,'(i1,''.'',i1,''.xyz'')') isec, i

      open (file = fname, newunit = io_xyz)

      write(io_xyz,*)natf(i)
      write(io_xyz,*)

      do j=1,natf(i)
         write(io_xyz,'(a2,5x,3F18.8)') toSymbol(iatf(j,i)), xyzf(1:3,j,i) * autoaa 
         cema(1:3,i) = cema(1:3,i) + xyzf(1:3,j,i) * iatf(j,i)
         z = z + iatf(j,i)
      enddo

      close(io_xyz)

     cema(1:3,i) = cema(1:3,i) / z
    enddo 
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  k=0
  do i = 1, nfrag
    do j = 1, i
      k = k + 1
      rf(k) = 0
      if (i /= j) then
         rf(k) = sqrt((cema(1,i)-cema(1,j))**2 + &
                      (cema(2,i)-cema(2,j))**2 + &
                      (cema(3,i)-cema(3,j))**2) * autoaa
      endif
    enddo
  enddo   
  
  if(nfrag > 1) then
    write(*,'(2x,a)') 'inter fragment distances (Angst.)'
    call print_matrix(rf,nfrag,0)
  else
    return
  endif
  
  ! PBE0/SVx for semi, SV(P) is too costly then but negligible in a DFT
  ! save original run parameters
  isave = bas
  jsave = func
  ksave = ihamilt
  gsave = gfnver   !save gfnver

  ipcalc =  .False.

  ! save etemp for XTB !THIS IS NOT ETEMP! BUT AVERAGE TEMP OF FRAGMENT!
  dsave = eTemp

  !> if nothing is defined, calculate the IPs/EAs with SV(P) for + and TZVP for -
  if ( mchrg < 0 ) then
    if ( iprog == 3 ) then      ! ORCA
              bas = 7          !ma-def2-TZVP
      if(ecp) bas = 11         !def2-TZVP

    elseif ( iprog == 2 ) then  ! TM
              bas = 11         !def2-TZVP
      if(ecp) bas = 11         !def2-TZVP
    endif
          
  else
            !bas = 3           !SV(P)
            bas = isave
    if(ecp) bas = 9           !def2-SV(P)
  endif
  
  itry = 1

  ! If IP calculation fails, try other QC codes,
  ! especially XTB (build-in) and ORCA (free)

  ! MOPAC, XTB, XTB2,ORCA
  if (iprog == 1) then
    useprog(1) = iprog
    useprog(2) = 7   
    useprog(3) = 8   
    useprog(4) = 3  

  ! TM, ORCA, XTB, XTB2
  elseif (iprog == 2) then
    useprog(1) = iprog
    useprog(2) = 3  
    useprog(3) = 7   
    useprog(4) = 8   

  ! ORCA, TMOL, XTB, XTB2
  elseif (iprog == 3) then
    useprog(1) = iprog
    useprog(2) = 2   
    useprog(3) = 7   
    useprog(4) = 8  
  
  ! MSINDO, XTB, XTB2, ORCA
  elseif (iprog == 4) then
    useprog(1) = iprog
    useprog(2) = 7   
    useprog(3) = 8   
    useprog(4) = 3   
  
  ! MNDO ,XTB, XTB2, ORCA
  elseif (iprog == 5) then
    useprog(1) = iprog
    useprog(2) = 7   
    useprog(3) = 8   
    useprog(4) = 3   
  
  ! XTB, XTB, XTB2, ORCA
  elseif (iprog == 7) then
    useprog(1) = iprog
    useprog(2) = iprog
    useprog(3) = 8
    useprog(4) = 3
  
  ! XTB2, XTB, ORCA
   elseif (iprog == 8) then
    useprog(1) = iprog
    useprog(2) = iprog
    useprog(3) = 7
    useprog(4) = 3
   endif
  
  call timing(t1,w1)
  
  do 
    progi=useprog(itry)
   
    ! defaults for the IP calc. (PM6 for MOPAC, OM2 for MNDO99 if ORCA fails)      
    if (progi == 1) ihamilt = 4
    if (progi == 5) ihamilt = 6
    
    ! If IP program is XTB
    if (progi == 7) then
      etemp  = 300.0d0
      ipcalc = .True.
      gfnver = 1
  
    elseif (progi == 8) then
      etemp  = 300.0d0
      ipcalc = .True.
      gfnver = 3

    elseif ( progi == 3 .and. method == 3 .and. etemp < 3000.0_wp ) then
      etemp  = 3000.0d0
    endif
    
    call qcstring(progi,line,line2) 
  
    if ( mchrg < 0 ) then !method == 4) then
      write(*,'(/,'' computing EAs with '',(a20),'' at (K) '',f7.0)')trim(line2),dsave
    else
      write(*,'(/,'' computing IPs with '',(a20),'' at (K) '',f7.0)')trim(line2),dsave
    endif
  
    fragip  = 0
    lowest_neut = 0.0_wp
    lowest_ion  = 0.0_wp
    
    !sn = 0
    !sp = 0
  
frg:do i = 1,nfrag
      gsen  = 0.0d0
      gsep  = 0.0d0
      metal = .False.
  
      !> check for metal 
      if ( .not. nometal ) then
        do k = 1, natf(i)
          if ( iatf(k,i) >= 22 .and. iatf(k,i) <= 30 ) metal = .True.
        enddo
      endif
    
      ! find spin for ion and neutral of metal              
      if (metal) then
        fiter = 3
        call getspin(natf(i),iatf(1,i),0,neutfspin)
        call getspin(natf(i),iatf(1,i),mchrg,ionfspin)
   
      ! no metal (That means eqm (iniqm.f), will assign spin by itself)
      else  
        metal = .False.
        fiter = 1
        neutfspin = -1
        ionfspin = -1                  
      endif
    
      ! MOPAC IP is unreliable for H and other atoms           
      if(  progi == 1.and.natf(i) == 1) then
        if ( mchrg < 0 ) stop 'MOPAC CANT BE USED FOR EA!'
        if ( mchrg > 1 ) stop 'MOPAC CANT BE USED FOR MULTIPLE CHARGES!'
        !^ can be fixed, but MOPAC should not be used anyways, so...
        fragip(i,1:abs(mchrg)) =  valip(iatf(1,i))
        E_neut    = 1.d-6
        E_ion    = fragip(i,mchrg) * evtoau 
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> For everything else than MOPAC
      !> Calculate the netr. and ion energies to get IPs/EAs
      if ( progi /= 1 )then

mult_n: do k = 1, fiter         !ITER OVER MULTIPLICITES

          if (metal) neutfspin = neutfspin + 2

          !> 1. Calculate Neutral energy (mcharge=0)
          call eqm(progi,natf(i),xyzf(1:3,1,i),iatf(1,i),0,neutfspin, &
            etemp,.true.,ipok,E_neut,nel,nb,ECP,spec_calc)
    
          if(metal .and. lowest_neut > E_neut)then
            !gsen(k) = E_neut
            !sn(k)   = neutfspin
            spin_neut   = neutfspin
            lowest_neut = E_neut
          endif
          if (metal) then
            write(*,*) 'E-Neut',E_neut 
              write(*,*) 'neut spin', neutfspin
            endif
        enddo mult_n
  
          !> 2. Calculate Ion energies (mcharg= +/- mchrg)
chrg:   do lpchrg = 1, abs(mchrg)
          if (mchrg < 0 ) then
            dump_chrg = -1 * lpchrg
          else
            dump_chrg = lpchrg
          endif

mult_i:   do k=1,fiter         !ITER OVER MULTIPLICITES

            if (metal) ionfspin  = ionfspin  + 2

            call eqm(progi,natf(i),xyzf(1:3,1,i),iatf(1,i),dump_chrg,ionfspin,  &
              etemp,.true.,ipok,E_ion,nel,nb,ECP,spec_calc)

            if (metal .and. lowest_ion > E_ion )then
              spin_ion    = ionfspin
              lowest_ion  = E_ion
            endif
            if (metal) then
              write(*,*) 'E-ION', E_ion 
              write(*,*) 'ION SPin', ionfspin
            endif

          enddo mult_i

    
            ! Select lowest values - calculate vertical IP/EA from groundstate neutral to groundstate ion
            ! in regards to spin multiplicity
            !if(metal) then
            !  E_neut = minval(gsen)
            !  E_ion  = minval(gsep)
            !  ! save neutral-ion (of metal) lowest energy spin
            !  sn0 = 0
            !  sp0 = 0
            !  lowest_neut  = E_neut
            !  lowest_ion   = E_ion
  
            !    if (gsen(j)  ==  E_neut)  sn0 = sn(j)              
            !    if (gsep(j)  ==  E_ion)   sp0 = sp(j)
            !  enddo
            !endif
    
            if (E_ion /= 0 .and. E_neut /= 0) then
              fragip(i,lpchrg) = (E_ion - E_neut) * autoev
            
              if(metal) then
              fragip(i,lpchrg) = (lowest_ion - lowest_neut) * autoev
              endif
  

              !> For 3d metal containing species
              if (metal) then
                write(*,'('' fragment '',i2,'' E(N)='',F12.4,''  E(I)='',F12.4,5x,'' &
                  &       IP/EA(eV)='',F8.2,5x,'' Mult.:'',i2,'' (N) and '',i2,'' (I)'')') &
                  &       i,lowest_neut,lowest_ion,fragip(i,lpchrg),spin_neut,spin_ion
                  !&       i,E_neut,E_ion,fragip(i,lpchrg),sn0,sp0
              else
                write(*,'('' fragment '',i2,'' E(N)='',F12.4,''  E(I)='',F12.4,5x,'' &
                  &       IP/EA(eV)='',F8.2)') i,E_neut,E_ion,fragip(i,lpchrg)
              endif
   
              !> set some boundaries, but not sure how accurate these are
              if ( mchrg < 0 ) then 
                if (fragip(i,lpchrg) > 40.0_wp .or. fragip(i,lpchrg) < -35.0_wp)then 
                  ipok = .false.
                endif
              elseif (mchrg == 1) then
                if (fragip(i,lpchrg) < 0.0_wp  .or. fragip(i,lpchrg) > 50.0_wp) then 
                  ipok = .false.
                endif
              elseif ( mchrg > 1 ) then
                if (fragip(i,lpchrg) < 0.0_wp  .or. fragip(i,lpchrg) > 100.0_wp) then 
                  ipok = .false.
                endif
              elseif ( mchrg > 1 ) then
                if (fragip(i,lpchrg) < 0.0_wp  .or. fragip(i,lpchrg) > 100.0_wp) then 
                  ipok = .false.
                endif
              endif
            endif

          !!!!!
          !enddo mult
        enddo chrg
      !!!!!
      endif 
    enddo frg
   
    ! if failed try another code      
    if ( .not. ipok ) then
      itry = itry + 1
      if (itry <= 3) then
         write(*,*) '* Try: ', itry, ' failed *'
         cycle
      else
      ! total failure, use Mpop in main                 
        fragip(1:nfrag,1:mchrg)=0
        ipok=.false.
        exit
      endif
    else
      exit ! finish all good
    endif   

   enddo
   
   ! restore original settings
   etemp   = dsave
   bas     = isave
   func    = jsave
   ihamilt = ksave
   gfnver  = gsave

   ipcalc  = .False.
   
   call timing(t2,w2)
   if( mchrg < 0 ) then ! method == 4)then
      write(*,'(/,'' wall time for EA (s)'',F10.1,/)')(w2-w1)
   else
      write(*,'(/,'' wall time for IP (s)'',F10.1,/)')(w2-w1)
   endif
    
    
end subroutine analyse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avg_frag_struc(nuc,iat,iatf,axyz,list, nfrag, natf, xyzf)

  integer :: nuc  
  integer :: natf(10)
  integer :: iat (nuc)
  integer :: nfrag
  integer :: list(nuc)
  integer :: iatf(nuc,10)
  integer :: idum(nuc,10)
  integer :: i,j,k

  real(wp) :: axyz(3,nuc)
  real(wp) :: nxyz(3,nfrag) 
  real(wp) :: xyzf(3,nuc,10)
  real(wp) :: dum (3,nuc,10)
    
  
  xyzf = 0
  iatf = 0
  do i=1,nuc
    j=list(i)
    xyzf(1:3,i,j)=axyz(1:3,i)
    iatf(    i,j)=iat(    i)
  enddo   

  dum  = xyzf
  idum = iatf
  
  xyzf = 0
  iatf = 0 
  do i=1,nfrag
    k=0
    do j=1,nuc
      if(idum(j,i) /= 0)then
        k=k+1
        xyzf(1:3,k,i)=dum(1:3,j,i)     
        iatf(    k,i)=idum(   j,i)     
      endif
    enddo   
    natf(i)=k
  enddo    
  
  do i = 1, nfrag
    do j = 1, natf(i)
      nxyz(1:3,i) = xyzf(1:3,j,i)
!      write(*,*)  nxyz(:,i) 
    enddo
 !   write(*,*)  
  enddo

end subroutine avg_frag_struc 
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
! FUNCTION FOR ATOM IP (USED ONLY FOR MOPAC)
function valip(i) result(get_ip)

  integer  :: i

  real(wp) :: ip(94), get_ip
  
  ip(1 )=-13.598
  ip(2 )=-24.587
  ip(3 )=-3.540   
  ip(4 )=-5.600
  ip(5 )=-8.298
  ip(6 )=-11.260
  ip(7 )=-14.534
  ip(8 )=-13.618
  ip(9 )=-17.423
  ip(10)=-21.565
  ip(11)=-3.091
  ip(12)=-4.280
  ip(13)=-5.986
  ip(14)=-8.152
  ip(15)=-10.487
  ip(16)=-10.360
  ip(17)=-12.968
  ip(18)=-15.760
  ip(19)=-2.786
  ip(20)=-3.744
  ip(21)=-9.450
  ip(22)=-10.495
  ip(23)=-10.370
  ip(24)=-10.642
  ip(25)=-13.017
  ip(26)=-14.805
  ip(27)=-14.821
  ip(28)=-13.820
  ip(29)=-14.100
  ip(30)=-4.664
  ip(31)=-5.999
  ip(32)=-7.899
  ip(33)=-9.789
  ip(34)=-9.752
  ip(35)=-11.814
  ip(36)=-14.000
  ip(37)=-2.664
  ip(38)=-3.703
  ip(39)=-7.173
  ip(40)=-8.740
  ip(41)=-10.261
  ip(42)=-11.921
  ip(43)=-12.59
  ip(44)=-14.214
  ip(45)=-15.333
  ip(46)=-8.337
  ip(47)=-17.401
  ip(48)=-4.686
  ip(49)=-5.786
  ip(50)=-7.344
  ip(51)=-8.608
  ip(52)=-9.010
  ip(53)=-10.451
  ip(54)=-12.130
  ip(55)=-2.544
  ip(56)=-3.333
  ip(57)=-7.826
  ip(58)=-7.594
  ip(59)=-4.944
  ip(60)=-4.879
  ip(61)=-4.813
  ip(62)=-4.754
  ip(63)=-4.615
  ip(64)=-7.915
  ip(65)=-4.617
  ip(66)=-4.566
  ip(67)=-4.520
  ip(68)=-4.487
  ip(69)=-4.441
  ip(70)=-4.378
  ip(71)=-5.428
  ip(72)=-8.419
  ip(73)=-10.786
  ip(74)=-12.293
  ip(75)=-13.053
  ip(76)=-15.450
  ip(77)=-17.779
  ip(78)=-19.695
  ip(79)=-21.567
  ip(80)=-5.521
  ip(81)=-6.108
  ip(82)=-7.417
  ip(83)=-7.286
  ip(84)=-8.417
  ip(85)=-10.7
  ip(86)=-10.748
  ip(87)=-2.637
  ip(88)=-3.412
  ip(89)=-6.97
  ip(90)=-9.951
  ip(91)=-8.09
  ip(92)=-9.115
  ip(93)=-9.243
  ip(94)=-6.324
  
  get_ip = -ip(i)

end function valip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module qcxms_analyse
