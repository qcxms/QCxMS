module qcxms_write_fragments
  use common1
  use cidcommon
  use qcxms_analyse, only: analyse
  use qcxms_fragments
  use qcxms_utility
  use xtb_mctc_accuracy, only : wp
  implicit none

  contains
    
  subroutine manage_fragments(nuc, iat, xyz, axyz, mchrg, qat, spin, mass, &
    &  imass, iprog, aTlast, itrj, icoll, isec, list, chrgcont,  &
    &  tcont, nfrag, metal3d, ECP, btf, maxsec, dtime, asave, io_res)

    integer, intent(in) :: nuc,  iprog
    integer  :: iat(nuc)
    integer  :: i, j, k, l, m, n
    integer  :: io_res
    integer  :: itrj, icoll,isec
    integer  :: list(nuc)
    integer  :: imass(nuc)
    integer, intent(out)  :: nfrag
    integer  :: fragat(200,10)
    integer  :: natf(10)
    integer  :: nearest_chrg
    integer  :: frag_number
    integer  :: tcont
    integer  :: maxsec
    integer  :: idum1(200)
    integer  :: idum2(200)
    integer  :: mchrg, mz_chrg
    integer  :: count_ip
    integer  :: mat_index(2)
    integer  :: MATsize
    integer :: iatf(nuc,10)
    integer, allocatable :: save_fragID(:), save_chrgID(:)
    
    real(wp),intent(inout) :: xyz(3, nuc) 
    real(wp),intent(in)    :: axyz(3, nuc)
    real(wp) :: mass(nuc)
    real(wp) :: qat(nuc)
    real(wp) :: spin(nuc)
    real(wp) :: fragm(10)
    real(wp) :: fragspin(10)
    real(wp) :: btf
    real(wp) :: aTlast
    real(wp) :: z
    real(wp) :: largest_chrg, frag_atms
    real(wp),intent(inout) :: chrgcont
    real(wp) :: chrgcont2
    real(wp) :: dtime
    real(wp) :: fragchrg(10)
    real(wp), allocatable :: norm_chrg(:)
    real(wp), allocatable :: fragip (:,:), fragchrg2(:,:)
    real(wp), allocatable :: ip_diff(:,:),ip_diff2(:,:)
    real(wp), allocatable :: ip_ranking(:)
    real(wp), allocatable :: fragchrg3(:)
    logical :: ipok
    logical :: metal3d 
    logical :: ECP 

    character(len=:), allocatable :: adum
    character(len=120) :: asave
    character(len=80)  :: fragf(nuc)


    l=0
    frag_number = 0

    ! calculate fragement structure, atoms, masses
    call fragment_structure(nuc,iat,xyz,3.0_wp,1,0,list)
    call fragmass(nuc,iat,list,mass,imass,nfrag,fragm,fragf,fragat)

    nfrag=maxval(list)

    allocate ( fragip(nfrag,  abs(mchrg)), fragchrg2(nfrag, abs(mchrg)) )
    allocate ( ip_diff(nfrag, abs(mchrg)) )
    allocate ( ip_diff2(nfrag,abs(mchrg)) )
    allocate ( fragchrg3(nfrag))
    allocate ( norm_chrg(nfrag))

    write(*,'('' fragment assigment list:'',80i1)')(list(k),k=1,nuc)

    !> compute fragment IP/EA per frag. and chrg.
    if ( method == 3 .and. .not. Temprun ) then ! fix average geometry for CID (axyz)
         call analyse(iprog,nuc,iat,iatf,xyz,list,nfrag,aTlast,fragip, mchrg, &
                  natf,ipok,icoll,isec,metal3d,ECP)
    else
       call analyse(iprog,nuc,iat,iatf,axyz,list,nfrag,aTlast,fragip, mchrg, &
                  natf,ipok,icoll,isec,metal3d,ECP)
    endif

    ! compute charge from IPs at finite temp.
    !fragchrg2 = chrgcont !mchrg ! 1.0_wp
    fragchrg2 = mchrg 
    fragchrg3 = 0.0_wp 

    !> get the relative ionization potentials, depending on charges 
    !> i.e. 1->2, 2->3 etc.
    if ( nfrag > 1 ) then

      do i = 1, nfrag
        fragip(i,0) = 0.0_wp
        do j = 1, abs(mchrg)
          ip_diff(i,j) = fragip(i,j) - fragip(i,j-1)

          !> H-Atoms have to be considerd seperately
          !if ( iatf(natf(i),i) == 1 .and. j == 1 ) then
          !  ip_diff(i,j) = fragip(i,j)
          !  write(*,*) 'FRAG IP H-ATOM', ip_diff(i,j)
          !elseif ( iatf(natf(i),i) == 1 .and. j > 1 ) then
          !  ip_diff(i,j) = huge(0.0_wp)
          !  write(*,*) 'FRAG IP 2nd H-ATOM', ip_diff(i,j)
          !endif

        enddo

      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> determine the ip ranking by charge and fragment
      MATsize = size(ip_diff)
      allocate(ip_ranking(MATsize))
      allocate(save_fragID(MATsize))
      allocate(save_chrgID(MATsize))


      ip_ranking = 0.0_wp

      write(*,'(''Pos.'',5x,'' IP: '',3x, '' Frag. Chrg. '', 2x)') 

      !> we want to only boltzmann weigh the 2nd, 3rd, ... charge, because
      !  otherwise the second charge will be double counted.
      !> determine the 1st charge and save, weigh the 2nd charge etc.
      do count_ip = 1, MATsize 

        !> locate the positions in the matrix and mask them for each IP-diff.
        !  this accounts for the correct IP ranking
        if ( mchrg > 0 ) then
          mat_index = minloc(ip_diff, mask = ip_diff > ip_ranking(count_ip-1))
        else
          mat_index = maxloc(ip_diff, mask = ip_diff < ip_ranking(count_ip-1))
        endif

        !> save the indexes for the next run
        ip_ranking(count_ip) = ip_diff(mat_index(1),mat_index(2))

        save_fragID(count_ip) = mat_index(1)
        save_chrgID(count_ip) = mat_index(2)

        if(verbose) write(*,'((i2), 5x, (f5.2), 3x, 2(i2))') &
            & count_ip,  ip_ranking(count_ip), mat_index(1), mat_index(2)

      enddo

      !> set 2nd IP value so it can be manipulated without losing info
      ip_diff2    = ip_diff 

      !> set the value with lowest IP to large number, so boltz will mostly
      !  ignore it -> important for correct values
      do count_ip = 1, abs(mchrg)-1
        ip_diff2(save_fragID(count_ip),save_chrgID(count_ip)) = huge(0.0_wp)
      enddo

      !> for negative ions
      if ( mchrg < 0 ) ip_diff2    = -1.0_wp * ip_diff2 

      !> do boltzmann for the fragment IPs
      call boltz(2,nfrag,abs(mchrg),aTlast*btf,ip_diff2,fragchrg2)

      !> set all charges to = 1 for all strucs that have to be ignored 
      do count_ip = 1, abs(mchrg)-1
         fragchrg2(save_fragID(count_ip),save_chrgID(count_ip)) = 1.0_wp
      enddo


      deallocate(ip_ranking)
      deallocate(save_fragID)
      deallocate(save_chrgID)

    endif


    !> sum the charges to get overall value for the entire fragments
    write(*,*) "fragchrg2"
    if ( nfrag > 1) then
      do i = 1, nfrag
        do j = 1, abs(mchrg)
          write(*,*) i,j,fragchrg2(i,j) !* chrgcont !(chrgcont / abs(mchrg))
          fragchrg3(i) = fragchrg3(i) + fragchrg2(i,j) !* chrgcont 
            ! (chrgcont / abs(mchrg))
        enddo
      enddo
    else
      fragchrg3 = chrgcont
    endif

    write(*,*) 
    write(*,*) 'Summed charges per fragment'
    do i = 1, nfrag
      write(*,*) i, fragchrg3(i)
    enddo


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! continune trajectory for frag with largest positive charge
    !if ( method /= 2 .and. method /= 4 ) then
    if ( mchrg > 0 ) then
      largest_chrg = -1
      do i = 1, nfrag
          frag_atms = float(natf(i))
          if ( frag_atms * fragchrg3(i) > largest_chrg ) then
             largest_chrg = frag_atms * fragchrg3(i)
             tcont = i ! set the NUMBER of the frag with largest charge
          endif
      enddo

    ! continune trajectory for frag with largest negative charge
    else
      largest_chrg = 1
      do i = 1, nfrag
          frag_atms = float(natf(i))
          if ( frag_atms * fragchrg3(i) < largest_chrg ) then
             largest_chrg = frag_atms * fragchrg3(i)
             tcont = i ! set the NUMBER of the frag with largest charge
          endif
      enddo
    endif


    write(*,*) 
    write(*,*) 'Normalized charges per fragment'
    do i = 1, nfrag
      norm_chrg(i) = fragchrg3(i) / fragchrg3(tcont) 
      write(*,*) i, norm_chrg(i)
    enddo

    chrgcont2 = chrgcont   ! <- normalized
    chrgcont =  norm_chrg(tcont)

    ! if not the first run, take charge from earlier run
    !chrgcont2 = chrgcont / abs(mchrg)
    !chrgcont2 = chrgcont   ! <- normalized

  
    ! do not continue if there are a) no frags or b) max reached
    if ( nfrag == 1 .or. isec == maxsec + 1 ) then
       tcont = 0

    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    elseif ( nfrag > 1 .and. mchrg > 1 ) then
      write(*,*) 'Angepasste charges per fragment'
      nearest_chrg = nint(fragchrg3(tcont))
      if ( nearest_chrg == 0 ) nearest_chrg = 1
      do i = 1, nfrag
        fragchrg3(i) = norm_chrg(i) * nearest_chrg
        write(*,*) i, fragchrg3(i)
      enddo
      write(*,*) 

    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    elseif ( nfrag > 1 .and. mchrg == 1 ) then
      do i = 1, nfrag
        fragchrg3(i) = fragchrg3(i) * chrgcont2 !/ abs(mchrg)
      enddo
      chrgcont = fragchrg3(tcont) 
    endif
    write(*,*) 

    if ( nfrag > 1 )then
      mchrg = nint(fragchrg3(tcont))
      if ( mchrg == 1 ) chrgcont2 = fragchrg3(tcont) 
      write(*,*)'NEW MCHRG', mchrg 
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute charge/spin on fragment(s) according to Mpops
    fragspin = 0
    fragchrg = 0
    do k = 1, nuc
      frag_number = list(k)

       ! dftb
       if ( prog == 0 ) then
          call valel(iat(k),z)
          fragchrg(frag_number) = fragchrg(frag_number) + z - qat(k)

       ! mopac
       elseif (prog == 1 .or. prog == 5) then
          fragchrg(frag_number) = fragchrg(frag_number) + qat(k)

       ! msindo
       elseif (prog == 4) then
          call valel(iat(k),z)
          fragchrg(frag_number) = fragchrg(frag_number) + z - qat(k)

       ! tm/orca/xtb
       else
          fragchrg(frag_number) = fragchrg(frag_number) + qat(k)
       endif

       fragspin(frag_number) = fragspin(frag_number) + spin(k)
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! save masses of fragments for mspec (->mass.agr)
    ! save all fragments
    write(*,'(''  mass'',16x,''formula'',22x,''q pop'',3x,''spin'', &
      &      4x,''|q IPB|'',2x,''diss time (ps)'')')

loop:do j = 1, nfrag

       adum=''
       if( j == tcont ) adum = ' ~'

       if ( fragm(j) < 10 ) then
          write(*,'('' M='',F4.2,3x,a20,2x,4i4,3F8.3,f9.3,a)')        &
          &    fragm(j),trim(fragf(j)),itrj,icoll,isec,j,fragchrg(j), &
!          &    fragspin(j),fragchrg3(j)*chrgcont,dtime,adum
!          &    fragspin(j),norm_chrg(j),dtime,adum
!          &    fragspin(j),fragchrg3(j)*chrgcont2,dtime,adum
!          &    fragspin(j),fragchrg3(j),dtime,adum
          &    fragspin(j),fragchrg3(j),dtime,adum
       endif

       if( fragm(j) >= 10 .and. fragm(j) < 100 ) then
          write(*,'('' M='',F5.2,2x,a20,2x,4i4,3F8.3,f9.3,a)')        &
          &    fragm(j),trim(fragf(j)),itrj,icoll,isec,j,fragchrg(j), &
!          &    fragspin(j),fragchrg3(j)*chrgcont,dtime,adum
!          &    fragspin(j),norm_chrg(j),dtime,adum
!         &    fragspin(j),fragchrg3(j)*chrgcont2,dtime,adum
          &    fragspin(j),fragchrg3(j),dtime,adum
       endif

       if (fragm(j) >= 100 .and. fragm(j) < 1000 ) then
          write(*,'('' M='',F6.2,1x,a20,2x,4i4,3F8.3,f9.3,a)')        &
          &    fragm(j),trim(fragf(j)),itrj,icoll,isec,j,fragchrg(j), &
!          &    fragspin(j),fragchrg3(j)*chrgcont,dtime,adum
!          &    fragspin(j),norm_chrg(j),dtime,adum
!          &    fragspin(j),fragchrg3(j)*chrgcont2,dtime,adum
          &    fragspin(j),fragchrg3(j),dtime,adum
       endif

       l = 0

       do k = 1, 200
          if(fragat(k,j) /= 0)then
             l = l + 1
             idum1(l) = fragat(k,j)
             idum2(l) = k
          endif
       enddo

       ! write into .res file:
       ! charge(qat), trajectory(itrj), number of collision event(icoll), 
       ! numbers of fragments (isec,j),
       ! types of atoms in fragment (l), atomic number (idum2), amount of this 
       ! atom typ (idum1) from m = 1 to l

!CIDEI: if ( method == 3 .or. method == 4 ) then
CIDEI: if ( method == 3 ) then !.or. method == 4 ) then
        if ( tcont > 0 ) then
             if ( tcont == j ) then
                write(asave,'(F10.7,i3,2i5,2i2,2x,i3,2x,20(i4,i3))')        &
                !&             fragchrg3(j)*chrgcont2,itrj,icoll,isec,j,  &
                &             fragchrg3(j),mchrg,itrj,icoll,isec,j,  &
                &             l,(idum2(m),idum1(m),m=1,l)



             elseif (tcont /= j) then
               mz_chrg = nint( fragchrg3(j))
               if ( mz_chrg == 0 ) mz_chrg = 1
               write(io_res,'(F10.7,i3,2i5,2i2,2x,i3,2x,20(i4,i3))')       &
                !&             fragchrg3(j)*chrgcont2,itrj,icoll,isec,j,  &
!>mchrg falschmuss nearest  !&             fragchrg3(j),mchrg,itrj,icoll,isec,j,  &
                &             fragchrg3(j),mz_chrg,itrj,icoll,isec,j,  &
                &             l,(idum2(m),idum1(m),m=1,l)
             endif

          elseif (tcont == 0 ) then
            if ( Temprun ) then
                write(io_res,'(F10.7,i3,2i5,2i2,2x,i3,2x,20(i4,i3))')      &
                &             fragchrg3(j),mchrg,itrj,icoll,isec,j, &
                &             l,(idum2(m),idum1(m),m=1,l)
            else  
                write(asave,'(F10.7,i3,2i5,2i2,2x,i3,2x,20(i4,i3))')        &
                &             fragchrg3(j),mchrg,itrj,icoll,isec,j,  &
                &             l,(idum2(m),idum1(m),m=1,l)
            endif
          endif

       else !EI has no icoll

          if ( tcont > 0 ) then
             if ( tcont == j ) then
                 write(asave,'(F10.7,i3,2i5,2i2,2x,i3,2x,20(i4,i3))')  &
                 &             fragchrg3(j),mchrg,itrj,isec,j,  &
                 &             l,(idum2(m),idum1(m),m=1,l)
             elseif (tcont /= j) then
                write(io_res,'(F10.7,i3,2i5,2i2,2x,i3,2x,20(i4,i3))')  &
                &             fragchrg3(j),mchrg,itrj,isec,j,   &
                &             l,(idum2(m),idum1(m),m=1,l)
             endif

          elseif (tcont == 0) then
               write(io_res,'(F10.7,i3,2i5,2i2,2x,i3,2x,20(i4,i3))')   &
               &             fragchrg3(j),mchrg,itrj,isec,j,    &
               &             l,(idum2(m),idum1(m),m=1,l)
          endif
       endif CIDEI

    enddo loop

    deallocate (fragip,     &
                fragchrg2,  &
                ip_diff,    &
                ip_diff2,   &
                fragchrg3 )

    write(*,*)

  end subroutine manage_fragments
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module qcxms_write_fragments
