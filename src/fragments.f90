module qcxms_fragments
  use covalent_radii, only: Rad
  use xtb_mctc_convert
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_symbols, only: toSymbol 
  implicit none

  contains

  subroutine fragmass(nat,iat,list,mass,imass,nfrag,fragx,fragf,fragat)
  
     integer  :: i,j,k,nn
     integer  :: nat,iat(nat),list(nat),nfrag
     integer  :: imass(nat)
     integer  :: fragat(200,*)
     integer  :: fragel(200,10)
  
     real(wp) ::  mass(nat)
     real(wp) ::  fragx(*)
     real(wp) ::  fragm(10) ! max 10 frags
  
  
     character(len=80) :: fragf(nat)
     character(len=5 ) :: aa
     character(len=80) :: aaa
     character(len=2 ) :: a2

     fragm =0
     fragel=0
     do i=1,nat
        fragm(list(i)) = fragm(list(i)) + mass(i)
        j = iat(i)
        ! isotope check
        if ( imass(i).gt.0 ) j = 100 + imass(i)
        fragel(j,list(i)) = fragel(j,list(i)) + 1
     enddo
  
     nfrag=0
     do i=1,10
        if(fragm(i).gt.0)then
           nfrag=nfrag+1
           aaa=''
           k=0
           do j=1,200
              if(fragel(j,i).ne.0)then
                 a2=toSymbol(j)
                 if(fragel(j,i).gt.99 )then
                    if(a2(2:2).eq.' ')then
                       write(aa,'(a1,i3)')toSymbol(j),fragel(j,i)
                       nn=4
                    else
                       write(aa,'(a2,i3)')toSymbol(j),fragel(j,i)
                       nn=5
                    endif
                 endif
                 if(fragel(j,i).lt.100)then
                    if(a2(2:2).eq.' ')then
                       write(aa,'(a1,i2)')toSymbol(j),fragel(j,i)
                       nn=3
                    else
                       write(aa,'(a2,i2)')toSymbol(j),fragel(j,i)
                       nn=4
                    endif
                 endif
                 if(fragel(j,i).lt.10)then
                    if(a2(2:2).eq.' ')then
                       write(aa,'(a1,i1)')toSymbol(j),fragel(j,i)
                       nn=2
                    else
                       write(aa,'(a2,i1)')toSymbol(j),fragel(j,i)
                       nn=3
                    endif
                 endif
                 aaa(k+1:k+nn)=aa(1:nn)
                 k=k+nn
              endif
           enddo
           fragf(nfrag) = aaa
           fragx(nfrag) = fragm(i) * autoamu 
           fragat(1:200,nfrag)=fragel(1:200,i)
        endif
     enddo
  
  end subroutine fragmass
  
  !
  !  Subroutine for definition of two or more fragments
  !  if at1 = 0 :  global separation in (nonbonded parts), beginning with atom at2
  !  if at1 and at2 <> 0 : define fragments only if a at1-at2 bond (and no other) exists
  !  if at1 and at2 = 0 : delete all fragment assignments
  !  no bond if rcut times the cov.dist.
  
  subroutine fragment_structure(nat,oz,xyz,rcut,at1,at2,frag)
  
     integer  :: at1,at2, nat
     integer  :: i,j
     integer  :: attotal,currentfrag
     integer  :: oz(nat),frag(nat)

     real(wp),intent(in) ::  xyz(3,nat)
     real(wp) :: rcov, r
     real(wp) :: rcut
  
     logical  :: finish
     logical  :: connect(nat,nat)
  
  
     connect(1:nat,1:nat)=.false.
  
     do i = 1,nat-1
        do j = i+1, nat
           r = sqrt((xyz(1,i)-xyz(1,j))**2 + (xyz(2,i)-xyz(2,j))**2 &
           & + (xyz(3,i)-xyz(3,j))**2)
           rcov = rcut * 0.5_wp * (Rad(oz(i)) + Rad(oz(j)))
           if(r.lt.rcov)then
              connect(i,j)=.true.
              connect(j,i)=.true.
           endif
        enddo
     enddo
     if ((at1.eq.0).and.(at2.eq.0)) then
        do i = 1,nat
           frag(i) = 1
        end do
        return
     else
  
        do i = 1,nat
           frag(i) = 0
        end do
  
        frag(at1) = 1
        attotal=1
  
        if (at2.ne.0) then
           connect(at1,at2) = .false.
           connect(at2,at1) = .false.
        endif
  
        finish=.false.
        currentfrag=0
  
        do while (attotal.ne.nat)
  
           currentfrag=currentfrag + 1
  
  !  cycle through atoms and find connected ones
  
           do while (.not.(finish))
              finish=.true.
              do i = 1,nat
                 if (frag(i).eq.currentfrag) then
                    do j = 1,nat
                       if (connect(i,j)) then
                          if (frag(j).eq.0) then
                             frag(j) = currentfrag
                             attotal = attotal + 1
                             finish = .false.
                          elseif (frag(j).eq.currentfrag) then
                             cycle
                          endif
                       endif
                    enddo
                 endif
              enddo
           enddo
  
  ! find the first atom in the next fragment
  
           do i = 1,nat
              if (frag(i).eq.0) then
                 frag(i) = currentfrag + 1
                 attotal = attotal + 1
                 exit
              endif
           enddo
           finish=.false.
        enddo
  
     endif
     return
  end subroutine fragment_structure

end module qcxms_fragments
