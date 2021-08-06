
! determine automatically the two parameters in the P(E)
! function for given energy range exc, the number of bonds
! nbnd and the input parameter ieebond (av energy/per bond in eV,
! normally 0.5 eV/bond)
subroutine getieeab(a,b,n,ityp,exc,nbnd,ieebond)
   use xtb_mctc_accuracy, only: wp
   implicit none

   integer  :: ityp,nbnd,k

   real(wp) :: a,b,n,pmax,ieemax,exc,eav,st,ieebond

   st=0.005
   a=0.0
   b=0.0
   k=0

   do
!10 a=a+st
     a=a+st
     b=b+st*7
     a=min(a,0.3)
     k=k+1
     call getmaxiee(a,b,n,ityp,exc,ieemax,pmax,eav)
     if(k.gt.10000) stop 'internal error inside getieeab'
     !if(eav/nbnd.lt.ieebond) goto 10
     if(eav/nbnd.ge.ieebond) exit 
   enddo

end

! get maximum value (vmax) at which energy (ieemax) and
! average energy s for given parameters in the P(E)
! distribution function
subroutine getmaxiee(a,b,n,ityp,exc,ieemax,vmax,s)
   use xtb_mctc_accuracy, only: wp
   implicit none

   integer  :: ityp
   real(wp) :: a,b,n,val,x,vmax,ieemax,exc,s,m

   x=0.001
   vmax=-1
   ieemax=0
   s=0
   m=0

   do
!10 if(ityp.eq.0) call gauss0(a,b,n,x,val)
      if(ityp.eq.0) call gauss0(a,b,n,x,val)
      if(ityp.eq.1) call poiss0(a,b,n,x,val)
      if(val.gt.vmax)then
         vmax=val
         ieemax=x
      endif
      x=x+0.01
      s=s+val*x
      m=m+val
      !if(x.lt.exc) goto 10
      if(x.ge.exc) exit
   enddo

   s=s/m

end


! Gaussian energy distribution, n= # val el, a an b are parameters

subroutine gauss0(a,b,n,x,p)
   use xtb_mctc_accuracy, only: wp
   implicit none

   real(wp) :: a,b,n,x,p

   p = exp( -a * (x-n*b)**2/n )

end


! Poisson energy distribution, n= # val el, a an b are parameters

subroutine poiss0(a,b,n,x,t18)
   use xtb_mctc_accuracy, only: wp
   implicit none

   real(wp) :: a,b
   real(wp) :: k,z,n,x,t18
   real(wp) :: t2,t8,t14,t17
   
   z=b
   k=1./a
   t2 = k/N
   t8 = dlog(z/k*N/x)
   t14 = dexp(t2*x*(0.1D1+t8)-1.D0*z)
   t17 = (t2*x+0.1D1)**(-0.5D0)
   t18 = t14*t17
end

