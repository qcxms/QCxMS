module qcxms_impact
   use qcxms_mdinit, only: ekinet
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   implicit none

   contains 

! scale factors velof such that the requested kinetic energy
! eimp+e0 is obtained
! velof = velo. scaling-factor, eimp=eimp,ff=fadd*nstep,e0=ekinstart
     subroutine impactscale(nuc,velo,mass,velof,eimp,ff,e0)
     
        integer  :: nuc
        integer  :: k,j
     
        real(wp) :: velo(3,nuc)
        real(wp) :: mass(nuc)
        real(wp) :: velof(nuc),ff,e0,eimp
        real(wp) :: E,T
        real(wp) :: Esoll,scal
        real(wp) :: v(3,nuc)
     
        Esoll=eimp*ff+e0
        scal=0.0
        k=0
     !10 k=k+1
        do
          k=k+1
          do j=1,3
             v(j,1:nuc)=velo(j,1:nuc)*(1.0d0+velof(1:nuc)*scal)
          enddo
          call ekinet(nuc,v,mass,E,T)
          scal=scal+0.0002
          if(Esoll-E.gt.0.001.and.k.lt.20000)then !goto 10
            cycle
          else
            exit
          endif
          
        enddo
     
        do j=1,3
           velo(j,1:nuc)=velo(j,1:nuc)*(1.0d0+velof(1:nuc)*scal)
        enddo
     
        if(k.ge.200000) then
           write(*,*) E,Esoll,k,scal
           stop 'error in impactscale'
        endif
     
     end subroutine impactscale
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
     subroutine calctrelax(n,emo,na,i,trelax,t)
     
        integer :: n,na,i
        integer :: k,j
     
        real(wp) :: emo(10*n),t,trelax
        real(wp) :: te,de,alp
     
        alp=0.5d0*autoev
     
        t=0
        k=i
     ! i=mo1
     ! n=number of atoms
     ! na=number of alpha elecs
     ! trelax => trelax=2000, if not changed by input
        do j=i+1,na
           de=(emo(k) - emo(j))
           te=trelax*exp(alp*de)
           t=t+te
           k=k+1
        enddo
     
     end subroutine calctrelax

end module qcxms_impact
