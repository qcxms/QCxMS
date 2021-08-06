module qcxms_write_fragments
  use common1
  use newcommon
  use qcxms_analyse, only: analyse
  use qcxms_fragments
  use qcxms_utility
  use xtb_mctc_accuracy, only : wp
  implicit none

  contains
    
  subroutine manage_fragments(nuc, iat, xyz, axyz, chrg, spin, mass, imass, &
    &  iprog, aTlast, itrj, icoll, isec, list, chrgcont,  &
    &  tcont, nfrag, metal3d, ECP, btf, maxsec, dtime, asave, io_res)

    integer, intent(in) :: nuc,  iprog
    integer  :: iat(nuc)
    integer  :: i, j, k, l, m
    integer  :: io_res
    integer  :: itrj, icoll,isec
    integer  :: list(nuc)
    integer  :: imass(nuc)
    integer, intent(out)  :: nfrag
    integer  :: fragat(200,10)
    integer  :: natf(10)
    integer  :: num_frags
    integer  :: frag_number
    integer  :: tcont
    integer  :: maxsec
    integer  :: idum1(200)
    integer  :: idum2(200)
    
    real(wp) :: xyz(3, nuc) 
    real(wp) :: axyz(3, nuc)
    real(wp) :: mass(nuc)
    real(wp) :: chrg(nuc)
    real(wp) :: spin(nuc)
    real(wp) :: fragm(10)
    real(wp) :: fragchrg(10), fragchrg2(10)
    real(wp) :: fragip (10)
    real(wp) :: fragspin(10)
    real(wp) :: btf
    real(wp) :: aTlast
    real(wp) :: z
    real(wp) :: largest_chrg, frag_atms
    real(wp) :: chrgcont,chrgcont2
    real(wp) :: dtime

    logical :: ipok
    logical :: metal3d 
    logical :: ECP 

    character(len=:), allocatable :: adum
    character(len=120) :: asave
    character(len=80)  :: fragf(nuc)


    l=0
    num_frags = 0
    !chrgcont2 = 1.0
    frag_number = 0

    ! calculate fragement structure, atoms, masses
    call fragment_structure(nuc,iat,xyz,3.0_wp,1,0,list)
    call fragmass(nuc,iat,list,mass,imass,nfrag,fragm,fragf,fragat)

    write(*,'('' fragment assigment list:'',80i1)')(list(k),k=1,nuc)

    ! compute fragment IP/EA
    if (method == 3 .or. method == 4) then ! fix average geometry for CID (axyz)
       call analyse(iprog,nuc,iat,xyz,list,aTlast,fragip, &
                  natf,ipok,icoll,isec,metal3d,ECP)
    else
       call analyse(iprog,nuc,iat,axyz,list,aTlast,fragip, &
                  natf,ipok,icoll,isec,metal3d,ECP)
    endif

    ! compute charge from IPs at finite temp.
    fragchrg2=1.0_wp

    if(nfrag.gt.1)then
       if (method == 2 .or. method ==4) fragip = -1.0_wp * fragip !neg.ion mode
       call boltz(2,nfrag,aTlast*btf,fragip,fragchrg2)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute charge/spin on fragment(s) according to Mpops
    fragchrg = 0
    fragspin = 0
    do k = 1, nuc
       frag_number = list(k)

       ! dftb
       if ( prog.eq.0 ) then
          call valel(iat(k),z)
          fragchrg(frag_number) = fragchrg(frag_number) + z - chrg(k)

       ! mopac
       elseif (prog.eq.1.or.prog.eq.5)then
          fragchrg(frag_number) = fragchrg(frag_number) + chrg(k)

       ! msindo
       elseif (prog.eq.4)then
          call valel(iat(k),z)
          fragchrg(frag_number) = fragchrg(frag_number) + z - chrg(k)

       ! tm/orca/xtb
       else
          fragchrg(frag_number) = fragchrg(frag_number) + chrg(k)
       endif
          fragspin(frag_number) = fragspin(frag_number) + spin(k)
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! take Mpop if IP calc failed
    if( .not.ipok )fragchrg2 = fragchrg


    largest_chrg = -1
    ! continune trajectory for frag with largest charge
    do i = 1, nfrag
       frag_atms = float(natf(i))
       if ( frag_atms * fragchrg2(i) > largest_chrg) then
          largest_chrg = frag_atms * fragchrg2(i)
          tcont = i ! set the NUMBER of the frag with largest chrg
       endif
    enddo

    ! if not the first run, take charge from earlier run
    chrgcont2 = chrgcont
  
    ! do not continue if there are a) no frags or b) max reached
    if ( nfrag == 1 .or. isec == maxsec + 1 ) then
       tcont = 0
    else 
       chrgcont = fragchrg2(tcont) * chrgcont2
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! save masses of fragments for mspec (->mass.agr)
    ! save all fragments
    write(*,'(''  mass'',16x,''formula'',18x,''q pop'',3x,''spin'', &
      &      4x,''|q IPB|'',2x,''diss time (ps)'')')

loop:do j = 1, nfrag

       adum=''
       if( j == tcont ) adum = ' ~'

       if ( fragm(j) < 10 ) then
          write(*,'('' M='',F4.2,3x,a20,2x,3i4,3F8.3,f9.3,a)')      &
          &    fragm(j),trim(fragf(j)),itrj,icoll,isec,fragchrg(j), &
          &    fragspin(j),fragchrg2(j),dtime,adum
       endif

       if( fragm(j) >= 10 .and. fragm(j) < 100 ) then
          write(*,'('' M='',F5.2,2x,a20,2x,3i4,3F8.3,f9.3,a)')      &
          &    fragm(j),trim(fragf(j)),itrj,icoll,isec,fragchrg(j), &
          &    fragspin(j),fragchrg2(j),dtime,adum
       endif

       if (fragm(j) >= 100 .and. fragm(j) < 1000 ) then
          write(*,'('' M='',F6.2,1x,a20,2x,3i4,3F8.3,f9.3,a)')      &
          &    fragm(j),trim(fragf(j)),itrj,icoll,isec,fragchrg(j), &
          &    fragspin(j),fragchrg2(j),dtime,adum
       endif

       l = 0

       do k = 1, 200

          if(fragat(k,j).ne.0)then
             l = l + 1
             idum1(l) = fragat(k,j)
             idum2(l) = k
          endif

       enddo

       ! write into .res file:
       ! charge(chrg), trajectory(itrj), number of collision event(icoll), 
       ! numbers of fragments (isec,j),
       ! types of atoms in fragment (l), atomic number (idum2), amount of this 
       ! atom typ (idum1) from m = 1 to l

CIDEI: if ( method == 3 .or. method == 4 ) then
          if (tcont > 0) then
             if (tcont == j) then
                write(asave,'(F10.7,2i5,2i2,2x,i3,2x,20(i4,i3))')        &
                &             fragchrg2(j)*chrgcont2,itrj,icoll,isec,j,  &
                &             l,(idum2(m),idum1(m),m=1,l)
             elseif (tcont /= j) then
                write(io_res,'(F10.7,2i5,2i2,2x,i3,2x,20(i4,i3))')       &
                &             fragchrg2(j)*chrgcont2,itrj,icoll,isec,j,  &
                &             l,(idum2(m),idum1(m),m=1,l)
             endif

          elseif (tcont == 0 ) then
            if ( Temprun ) then
                write(io_res,'(F10.7,2i5,2i2,2x,i3,2x,20(i4,i3))')      &
                &             fragchrg2(j)*chrgcont2,itrj,icoll,isec,j, &
                &             l,(idum2(m),idum1(m),m=1,l)
            else  
                write(asave,'(F10.7,2i5,2i2,2x,i3,2x,20(i4,i3))')        &
                &             fragchrg2(j)*chrgcont2,itrj,icoll,isec,j,  &
                &             l,(idum2(m),idum1(m),m=1,l)
            endif
          endif

       else !EI has no icoll

          if (tcont > 0) then
             if (tcont == j) then
                 write(asave,'(F10.7,2i5,2i2,2x,i3,2x,20(i4,i3))')  &
                 &             fragchrg2(j)*chrgcont2,itrj,isec,j,  &
                 &             l,(idum2(m),idum1(m),m=1,l)
             elseif (tcont /= j) then
                write(io_res,'(F10.7,2i5,2i2,2x,i3,2x,20(i4,i3))')  &
                &             fragchrg2(j)*chrgcont2,itrj,isec,j,   &
                &             l,(idum2(m),idum1(m),m=1,l)
             endif

          elseif (tcont == 0) then
               write(io_res,'(F10.7,2i5,2i2,2x,i3,2x,20(i4,i3))')   &
               &             fragchrg2(j)*chrgcont2,itrj,isec,j,    &
               &             l,(idum2(m),idum1(m),m=1,l)
          endif
       endif CIDEI

    enddo loop

    write(*,*)

  end subroutine manage_fragments
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  subroutine reset_structures (nuc, xyz, velo, iat, mass, imass, tcont, list, chrg, spin, grad, axyz)
!
!   integer  :: nuc
!   integer  :: iat(nuc),    store_iat(nuc)
!   integer  :: imass(nuc),  store_imass(nuc)
!!   integer  :: store_iat(nuc)
!!   integer  :: store_imass(nuc)
!   integer  :: list(nuc)
!   integer  :: k, j, i
!   integer, intent(in)  :: tcont
!
!   real(wp) :: xyz(3, nuc), store_xyz(3,nuc)
!   real(wp) :: velo(3,nuc), store_velo(3,nuc)
!   real(wp) :: mass(nuc)  , store_mass(nuc)
!!   real(wp) :: store_xyz(3,nuc)
!!   real(wp) :: store_velo(3,nuc)
!!   real(wp) :: store_mass(nuc)
!   real(wp) :: cema(3)
!   real(wp) :: count_iat
!   real(wp) :: chrg(nuc)
!   real(wp) :: spin(nuc)
!   real(wp) :: grad(3,nuc)
!   real(wp) :: axyz(3,nuc)
!   real(wp) :: velof(nuc)
!   real(wp) :: tadd
!   real(wp) :: eimp
!
!
!   real(wp), allocatable :: xyz0  (:, :)
!   real(wp), allocatable :: axyz0 (:, :)
!   real(wp), allocatable :: grad0 (:, :)
!   real(wp), allocatable :: velo0 (:, :)
!   real(wp), allocatable :: chrg0 (:)
!   real(wp), allocatable :: spin0 (:)
!   real(wp), allocatable :: mass0 (:)  
!   integer, allocatable  :: imass0(:)
!   integer, allocatable  :: iat0  (:)
!   integer, allocatable  :: list0 (:)
!
!   allocate(xyz0 (3,nuc),      &
!   &        axyz0(3,nuc),      &
!   &        grad0(3,nuc),      &
!   &        velo0(3,nuc),      &
!   &        chrg0(nuc),        &
!   &        spin0(nuc),        &
!   &        imass0(nuc),       &
!   &        mass0(nuc),        &
!   &        iat0(nuc),         &
!   &        list0(nuc))        
!
!
!
!   write(*,*) 'XYZ', xyz
!   write(*,*) 'IAT', iat
!
!    k = 0
!    count_iat = 0
!    cema=0
!    if (tcont > 0) then
!       do i=1,nuc
!          j=list(i)
!          if(j == tcont) then
!             k=k+1     
!             store_xyz  (1:3,k) = xyz (1:3,i)
!             store_velo (1:3,k) = velo(1:3,i)
!             store_iat  (    k) = iat (    i)  
!             store_mass (    k) = mass(    i)  
!             store_imass(    k) = imass(   i)  
!             count_iat          = count_iat + store_iat(k)
!             cema(1:3)          = cema(1:3) + store_xyz(1:3,k) * store_iat(k)
!          endif   
!       enddo
!    elseif (tcont == 0) then
!        do i=1,nuc
!           store_xyz (1:3,i) = xyz (1:3,i)
!           store_velo(1:3,i) = velo(1:3,i)
!           store_iat (    i) = iat (    i)  
!           store_mass(    i) = mass(    i)  
!           store_imass(   i) = imass(   i)  
!           count_iat         = count_iat + store_iat(i)
!           cema(1:3)         = cema(1:3) + store_xyz(1:3,i) * store_iat(i)
!        enddo
!    endif
!
!    ! move to Center-of-mass
!    cema(1:3) = cema(1:3) / count_iat
!
!!    ! if fragmented
!    if(tcont > 0)then
!       !if(nuc.le.10) nmax=nmax0/10
!       nuc = k
!    endif
!!
!!    ! do not continue with low masses (user)
!!    if( sum(mass(1:nuc)) * autoamu < minmass )then
!!       littlemass = .true.
!!       exit
!!    endif
!!
!!    ! do not continue with small fragments
!!    if(nuc <= 5)then
!!       small = .true.
!!       exit
!!    endif
!
!
!    deallocate(xyz,axyz,grad,velo,chrg,spin,iat,list,mass,imass)
!
!    allocate(xyz (3,nuc), &
!    &        axyz(3,nuc), &
!    &        grad(3,nuc), &
!    &        velo(3,nuc), &
!    &        chrg(nuc),   &
!    &        spin(nuc),   &
!    &        iat (nuc),   &
!    &        list(nuc),   &
!    &        imass(nuc),  &
!    &        mass(nuc))
!
!    do i=1,3
!       xyz(i,1:nuc) = store_xyz (i,1:nuc) - cema(i)
!    enddo
!
!    velo(1:3,1:nuc) = store_velo(1:3,1:nuc)
!    iat (    1:nuc) = store_iat (    1:nuc)
!    mass(    1:nuc) = store_mass(    1:nuc)
!    imass(   1:nuc) = store_imass(   1:nuc)
!    velof=0
!    tadd=0
!    eimp=0
!
!   write(*,*) 'IAT', iat
!
!  end subroutine reset_structures 


end module qcxms_write_fragments
