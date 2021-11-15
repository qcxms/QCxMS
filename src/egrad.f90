!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! important routine:
! get energy, gradient, charge and spin density on atoms from QC
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine egrad(first,nuc,xyz,iat,chrg,spin,etemp,E,grad,qat,aspin,ECP,gradfail)
   use common1
   use qcxms_tblite, only: get_xtb_egrad, gfn1_xtb, gfn2_xtb
   use qcxms_use_mndo, only: mndoout, mndograd 
   use qcxms_use_msindo, only: getmsindograd 
   use qcxms_use_orca, only: orcaout, rdorcagrad 
   use qcxms_use_turbomole, only: initm, rdtmgrad, setfermi
   use qcxms_utility, only: qccall, getspin
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only: toSymbol 
   implicit none

   integer :: nuc
   integer :: chrg
   integer :: spin
   integer :: iat(nuc)
   integer :: i,idum,stat

   real(wp) :: xyz (3,nuc)
   real(wp) :: grad(3,nuc)
   real(wp) :: qat(nuc) ! atomic charge (q-at)
   real(wp) :: aspin(nuc)
   real(wp) :: E,etemp
   real(wp) :: d3a1,d3a2,d3s8,d3atm
   real(wp) :: cn(nuc),dcn(3,nuc,nuc)
   real(wp) :: egcp,ggcp(3,nuc),disp

   character(len=20) :: atmp

   logical :: first
   logical :: ECP
   logical :: ok,gradfail

   disp   = 0
   E    = 0
   grad = 0
   aspin= 0
   qat= 0
   gradfail=.false.

!ccccccccccccccccccccccccccc
! dftb+ (DFTB3-D3(BJ))
!ccccccccccccccccccccccccccc

   if(prog.eq.0)then
      call dftbout(nuc,xyz,iat,chrg,spin,.true.,etemp,1.d-5)
      call qccall(0,'job.last')
      call dftbgrad(nuc,grad,qat,aspin,E)
      call checkqc(nuc,E,grad,qat,ok,method)
      if(.not.ok) then
         call dftbout(nuc,xyz,iat,chrg,spin,.false.,etemp,1.d-4)
         call qccall(0,'job.last2')
         call dftbgrad(nuc,grad,qat,aspin,E)
         call checkqc(nuc,E,grad,qat,ok,method)
         if(ok)write(*,*)'healed!'
      endif
   endif

!ccccccccccccccccccccccccccc
! mopac
!ccccccccccccccccccccccccccc

   if(prog.eq.1)then
      call getmopgrad(nuc,iat,xyz,grad,chrg,etemp,.false.,E,qat,aspin)
      call checkqc(nuc,E,grad,qat,ok,method)
   endif

!ccccccccccccccccccccccccccc
! TM
!ccccccccccccccccccccccccccc

   if(prog.eq.2)then
      call system('rm -rf gradient energy dscf_problem')
      open(unit=87,file='coord')
      write(87,'(a)')'$coord'
      do i=1,nuc
         write(87,'(3F16.10,10x,a2)')&
         &xyz(1,i),xyz(2,i),xyz(3,i),toSymbol(iat(i))
      enddo
      write(87,'(a)')'$end'
      close(87)
      if(etemp.gt.0) call setfermi(etemp)
      call qccall(2,'job.last')
      call rdtmgrad (nuc,grad,qat,aspin,E)
      call checkqc(nuc,E,grad,qat,ok,method)
   endif

!ccccccccccccccccccccccccccc
! ORCA
!ccccccccccccccccccccccccccc

   if(prog.eq.3)then
      call orcaout(nuc,xyz,iat,chrg,spin,etemp,.true.,ECP)
      call qccall(3,'job.last')
      call rdorcagrad('job.last',nuc,grad,qat,aspin,E)
! second attempt
      call checkqc(nuc,E,grad,qat,ok,method)
      if(.not.ok) then
         call qccall(3,'job.last2')
         call rdorcagrad('job.last2',nuc,grad,qat,aspin,E)
         call checkqc(nuc,E,grad,qat,ok,method)
         if(ok)write(*,*)'healed!'
      endif
   endif

!ccccccccccccccccccccccccccc
! MSINDO
!ccccccccccccccccccccccccccc

   if(prog.eq.4)then
      call getmsindograd(first,nuc,iat,chrg,spin,xyz,etemp,4,grad,E,qat,aspin)

      call checkqc2(nuc,E,grad,qat,ok,method)
      if(.not.ok) call getmsindograd(first,nuc,iat,chrg,spin,xyz,etemp,14,grad,E,qat,aspin)

      call checkqc2(nuc,E,grad,qat,ok,method)
      if(.not.ok) call getmsindograd(first,nuc,iat,chrg,spin,xyz,etemp,16,grad,E,qat,aspin)

      call checkqc(nuc,E,grad,qat,ok,method)
   endif

!ccccccccccccccccccccccccccc
! MNDO99
!ccccccccccccccccccccccccccc

   if(prog.eq.5)then
      call mndoout(nuc,xyz,iat,chrg,spin,etemp,5,10)
      call qccall(5,'job.last')
      call mndograd('job.last',nuc,grad,qat,aspin,E)
      call checkqc(nuc,E,grad,qat,ok,method)
      if(.not.ok) then
         call mndoout(nuc,xyz,iat,chrg,spin,etemp,5,20)
         call qccall(5,'job.last2')
         call mndograd('job.last2',nuc,grad,qat,aspin,E)
         call checkqc(nuc,E,grad,qat,ok,method)
         if(ok)write(*,*)'healed!'
      endif
      if(.not.ok) then
         call mndoout(nuc,xyz,iat,chrg,spin,etemp,4,2)
         call qccall(5,'job.last2')
         call mndograd('job.last2',nuc,grad,qat,aspin,E)
         call checkqc(nuc,E,grad,qat,ok,method)
         if(ok)write(*,*)'healed!'
      endif
   endif


!cccccccccccccccccccccccccc
!  XTB
!ccccccccccccccccccccccccc
   if(prog.eq.6)then
! idum = spin
      call getspin(nuc,iat,chrg,idum)
      call callxtb(nuc,xyz,iat,chrg,idum,etemp,E,grad,qat,aspin)
      call checkqc(nuc,E,grad,qat,ok,method)
      if(.not.ok)then
         gradfail=.true.
         write(*,*) 'GRAD failed!'
      endif
   endif

!cccccccccccccccccccccccccc
! GFN1-xTB
!ccccccccccccccccccccccccc
   if(prog.eq.7)then
      call getspin(nuc,iat,chrg,idum)
!      if (with_tblite) then
         call get_xtb_egrad(iat, xyz, chrg, idum, gfn1_xtb, etemp, &
            & "xtb.last", qat, E, grad, stat, .false.)
         ok = stat == 0
!      else
!         call qcxtb2(nuc,xyz,iat,chrg,idum,E,etemp,grad,qat,aspin,.false.)
!      end if
      call checkqc(nuc,E,grad,qat,ok,method)
      if(.not.ok)then
         gradfail=.true.
         write(*,*) 'QC calc failed!'
      endif
   endif

!cccccccccccccccccccccccccc
! GFN2-xTB
!ccccccccccccccccccccccccc
   if(prog.eq.8)then
      call getspin(nuc,iat,chrg,idum)
!      if (with_tblite) then
         call get_xtb_egrad(iat, xyz, chrg, idum, gfn2_xtb, etemp, &
            &               "xtb.last", qat, E, grad, stat, .false.)
         ok = stat == 0
!      else
!         call qcxtb2(nuc,xyz,iat,chrg,idum,E,etemp,grad,qat,aspin,.false.)
!      end if
      call checkqc(nuc,E,grad,qat,ok,method)
      if(.not.ok)then
         gradfail=.true.
         write(*,*) 'QC calc failed!'
      endif
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gcp and D3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(.not.ok) return

   if(gcp) then
      atmp='dft/sv'
      if(bas.eq.2)atmp='dft/svx'
      if(bas.eq.3.or.bas.eq.9)atmp='dft/sv(p)'
      if(bas.eq.4.or.bas.eq.10)atmp='dft/svp'
      call calcgcp(nuc,iat,xyz,.false.,atmp,egcp,ggcp)
      grad(1:3,1:nuc)=grad(1:3,1:nuc)+ggcp(1:3,1:nuc)
      E=E+egcp
   endif

! return with ec=0 if parameters are zero (see input.f)
! make d3 standalone only available for dftb(0), msindo(4) or mndo99(5)
! in TM and ORCA it will be employed intrinsically
   if (prog.eq.0 .or. prog.eq. 4 .or. prog.eq. 5)then
!      call gdispxtb(nuc,xyz,iat,.false.,grad,ec)
      call gdisp(nuc,iat,xyz,d3a1,d3a2,d3s8,d3atm,disp,grad,cn,dcn,.false.)
      E = E + disp
   endif
!ccccccccccccccccccccccccccc
! done
!ccccccccccccccccccccccccccc

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute norm of electronic gradient
function gnorm(nuc,grad)
   use xtb_mctc_accuracy, only: wp
   implicit none
   integer  :: nuc,i
   real(wp) :: gnorm 
   real(wp) :: grad(3,nuc)
   real(wp) :: gn

   gn=0
   do i=1,nuc
      gn=gn+grad(1,i)**2+grad(2,i)**2+grad(2,i)**2
   enddo

   gnorm=sqrt(gn)

end

! is E, grad and charge reasonable ?
subroutine checkqc(nuc,e,grad,qat,ok,method)
   use xtb_mctc_accuracy, only: wp
   implicit none

   integer  :: nuc,method
   real(wp) :: grad(3,nuc),e,qat(nuc)
   real(wp) :: gnorm,gn
   logical  :: ok


   ok=.false.

   if(abs(E).lt.1.d-8) then
      write(*,*)'energy seems to be bogus!'
      return
   endif

   gn=gnorm(nuc,grad)
   if(gn.lt.1.d-8.or.gn.gt.20.0) then
      write(*,*)'grad seems to be bogus!'
      E=0
      return
   endif

   if (method .eq. 0 .OR. method .eq. 1 .OR. method .eq. 3) then
      if(abs(maxval(qat)).lt.1.d-5) then
         write(*,*)'WF seems to be bogus!',sum(qat)
         E=0
         return
      endif
   endif
   ok=.true.

end

! is E, grad and charge reasonable ? same as above without print
subroutine checkqc2(nuc,e,grad,qat,ok,method)
   use xtb_mctc_accuracy, only: wp
   implicit none

   integer  :: nuc,method

   real(wp) :: grad(3,nuc),e,qat(nuc)
   real(wp) :: gnorm,gn

   logical  :: ok

   ok=.false.

   if(abs(E).lt.1.d-8) then
      return
   endif

   gn=gnorm(nuc,grad)
   if(gn.lt.1.d-8.or.gn.gt.20.0) then
      E=0
      return
   endif
   if (method .eq. 0 .OR. method .eq. 1 .OR. method .eq. 3) then
      if(abs(maxval(qat)).lt.1.d-5) then
         E=0
         return
      endif
   endif
   ok=.true.

end
