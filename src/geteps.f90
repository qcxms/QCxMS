module qcxms_eps
  use readcommon
  use qcxms_utility, only: valel
  use xtb_mctc_accuracy, only: wp
  implicit none

  contains

! get occ orbital energies for an estimate of the c impact energy
  subroutine geteps(fname,nuc,iat,chrg,eps,mopop,ncore,ihomo,nb)
  
     integer  :: nuc,iat(nuc),chrg,ihomo,nb
     integer  :: i,nn,ncore,homo,nel,ispin,no,j,mo
     integer  :: k,kk,nv
  
     real(wp) ::  eps(2,*)
     real(wp) :: mopop(5*nuc,nuc)
     real(wp) :: xx(50),vel,dum
  
    !character(len=80) :: line
     character(len=*)  :: fname
  
     mopop=0
     ncore=0
     nel=-chrg
     ispin=0
  ! get # of core levels and val electrons
     do i=1,nuc
        call valel(iat(i),vel)      !get val. e-
        ncore=ncore+iat(i)-int(vel) !sum over all core e-
        nel=nel+int(vel)            !sum over all val. e-
     enddo
     ncore=ncore/2  ! 2 e- per MO
     homo=nel/2     ! 2 e- per MO
     if(2*(nel/2).ne.nel) homo=nel/2+1 ! open-shell
     ihomo=homo        !ihomo= homo, only even e-
     nb=nel/2          !nb= homo, even or uneven e-
     
  ! read only nel/2 lowest virt
     nv=ihomo/2
  
     if(ihomo.gt.10*nuc.or.nb.gt.10*nuc) stop 'weird stuff in geteps'
  
     open(unit=11,file=fname)
     do 
        read(11,'(a)',iostat=iocheck) line
        if (iocheck < 0) exit
  
  ! orca, eps
        if(index(line,'NO   OCC          E').ne.0)then
           ispin=ispin+1
    !       write(*,*) 'SPIN check', ispin 
           if(ispin.eq.1)no=ihomo
           if(ispin.eq.2)no=nb
           do i=1,ncore
              read(11,'(a)') line
           enddo
           j=0
           do i=ncore+1,no+ncore+nv
              j=j+1
              read(11,'(a)') line
              call readl(line,xx,nn)
              if(j.gt.ihomo)  xx(4)=-xx(4)
              eps(ispin,j)=xx(4)
           enddo
        endif
        if(index(line,'LOEWDIN ATOM POPULATIONS PER MO').ne.0)then
           kk=0
          ! new set of mos
          do 
             if(index(line,' -------- ').ne.0)then
                do 
                   read(11,'(a)',iostat=iocheck) line
                   if (iocheck < 0) exit
  
                   call readl(line,xx,nn)
  
                   if(nn.lt.1) then
                      kk=kk+no-1
                      exit
                   else
                      no=nn
                   endif
                   j=int(xx(1))+1
                   do i=2,nn
                      mo=kk+i-1
                      if(mo.gt.ncore.and.mo.le.ihomo+ncore)then
                         mopop(mo-ncore,j)=xx(i)*0.01
                      endif
                   enddo
                enddo
                if (iocheck < 0) exit
                if(kk-ncore.gt.ihomo) exit 
             endif
          enddo
        endif
  
     enddo
     close (11)
  
     if(ispin.eq.1)eps(2,1:ihomo)=eps(1,1:ihomo)
  
     do i=1,ihomo
        dum=sum(mopop(i,1:nuc))
        if(abs(dum-1.0d0).gt.0.05)then
           write(*,*)'Loewdin populations are weird or not computed'
           write(*,*)'taking unity instead'
           mopop = 1
           exit
        endif
     enddo
  
  end subroutine geteps
  
end module qcxms_eps
