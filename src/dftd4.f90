module dftd4
   use xtb_mctc_accuracy, only: wp
   implicit none

   real(wp) :: pi,thopi,ootpi
   parameter ( pi = 3.141592653589793238462643383279502884197_wp )
   parameter ( thopi = 3._wp/pi )
   parameter ( ootpi = 0.5_wp/pi )

   real(wp) :: k1,k4,k5,k6
!  global ad hoc parameters
   parameter (k1=16.0)

!  new parameter set for covalent coordination number
   parameter (k4=4.10451_wp)
   parameter (k5=19.08857_wp)
!  the k6 parameter is not included in the original publication,
!  but is strictly necessary to obtain resonable results, beware.
   parameter (k6=2*11.28174_wp**2)

   real(wp),parameter :: zeff(*) = (/ &
   &  1,2,  &
   &  3,4,5,6,7,8,9,10,  &
   &  11,12,13,14,15,16,17,18,  &
   &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  &
   &  9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  &
   &  9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  &
   &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26/)

!  Semiempirical Evaluation of the GlobalHardness of the Atoms of 103
!  Elements of the Periodic Table Using the Most Probable Radii as
!  their Size Descriptors DULAL C. GHOSH, NAZMUL ISLAM 2009 in 
!  Wiley InterScience (www.interscience.wiley.com).
!  DOI 10.1002/qua.22202
!  values in the paper multiplied by two because 
!  (ii:ii)=(IP-EA)=d^2 E/dN^2 but the hardness
!  definition they use is 1/2d^2 E/dN^2 (in Eh)
   real(wp),parameter :: gam(*) = (/ &
  &0.47259288,0.92203391,0.17452888,0.25700733,0.33949086,0.42195412, &
  &0.50438193,0.58691863,0.66931351,0.75191607,0.17964105,0.22157276, &
  &0.26348578,0.30539645,0.34734014,0.38924725,0.43115670,0.47308269, &
  &0.17105469,0.20276244,0.21007322,0.21739647,0.22471039,0.23201501, &
  &0.23933969,0.24665638,0.25398255,0.26128863,0.26859476,0.27592565, &
  &0.30762999,0.33931580,0.37235985,0.40273549,0.43445776,0.46611708, &
  &0.15585079,0.18649324,0.19356210,0.20063311,0.20770522,0.21477254, &
  &0.22184614,0.22891872,0.23598621,0.24305612,0.25013018,0.25719937, &
  &0.28784780,0.31848673,0.34912431,0.37976593,0.41040808,0.44105777, &
  &0.05019332,0.06762570,0.08504445,0.10247736,0.11991105,0.13732772, &
  &0.15476297,0.17218265,0.18961288,0.20704760,0.22446752,0.24189645, &
  &0.25932503,0.27676094,0.29418231,0.31159587,0.32902274,0.34592298, &
  &0.36388048,0.38130586,0.39877476,0.41614298,0.43364510,0.45104014, &
  &0.46848986,0.48584550,0.12526730,0.14268677,0.16011615,0.17755889, &
  &0.19497557,0.21240778,0.07263525,0.09422158,0.09920295,0.10418621, &
  &0.14235633,0.16394294,0.18551941,0.22370139 /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This was added for QCxMS (former aocommon), only needed for subr, 'eself'
   real(wp),parameter :: gam3(*) = (/ &
  &-0.02448,0.178614,0.194034,0.154068,0.173892,0.167160,0.156306, &
  &0.161466,0.163314,0.170862,0.256128,0.189060,0.146310,0.136686, &
  &0.123558,0.122070,0.119424,0.115368,0.175938,0.152214,0.240300, &
  &0.204654,0.186552,0.178824,0.174474,0.205038,0.188772,0.180462, &
  &0.180354,0.176508,0.158250,0.139212,0.131868,0.119712,0.115476, &
  &0.108444,0.091032,0.076980,0.102642,0.095346,0.088266,0.086364, &
  &0.085254,0.088242,0.087774,0.088470,0.091314,0.090372,0.110862, &
  &0.093588,0.079908,0.074082,0.069342,0.064638,0.077826,0.059021, &
  &0.073614,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
  &0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000, &
  &0.000000,0.085236,0.077820,0.074112,0.072060,0.070188,0.069246, &
  &0.072042,0.070380,0.072570,0.108096,0.084150,0.070350,0.069540, &
  &0.064866,0.059687,0.000000,0.000000,0.000000,0.000000,0.000000, &
  &0.000000,0.000000,0.000000 /) 
                                                                   
!  pauling EN's 
   real(wp),parameter :: en(*) = (/ &
  &         2.200,3.000,0.980,1.570,2.040,2.550,3.040,3.440,3.980 &
  &        ,4.500,0.930,1.310,1.610,1.900,2.190,2.580,3.160,3.500 &
  &        ,0.820,1.000,1.360,1.540,1.630,1.660,1.550,1.830,1.880 &
  &        ,1.910,1.900,1.650,1.810,2.010,2.180,2.550,2.960,3.000 &
  &        ,0.820,0.950,1.220,1.330,1.600,2.160,1.900,2.200,2.280 &
  &        ,2.200,1.930,1.690,1.780,1.960,2.050,2.100,2.660,2.600 &
  &,0.79,0.89,1.10,1.12,1.13,1.14,1.15,1.17,1.18,1.20,1.21,1.22 &
  &,1.23,1.24,1.25,1.26,1.27,1.3,1.5,2.36,1.9,2.2,2.20,2.28,2.54 &
  &,2.00,1.62,2.33,2.02,2.0,2.2,2.2,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5/)

!  D3 radii
   real(wp),parameter :: rcov(*) = (/ &
  & 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, &
  & 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
  & 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
  & 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
  & 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
  & 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
  & 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
  & 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, &
  & 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
  & 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
  & 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
  & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
  & 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
  & 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
  & 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
  & 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
  & 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
  & 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
  & 3.82984466, 3.85504098, 3.88023730, 3.90543362 /)

!  r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
!  the large number of digits is just to keep the results consistent
!  with older versions. They should not imply any higher accuracy
!  than the old values
   real(wp),parameter :: r4r2(*) = (/ &
  &   2.00734898,  1.56637132,  5.01986934,  3.85379032,  3.64446594, &
  &   3.10492822,  2.71175247,  2.59361680,  2.38825250,  2.21522516, &
  &   6.58585536,  5.46295967,  5.65216669,  4.88284902,  4.29727576, &
  &   4.04108902,  3.72932356,  3.44677275,  7.97762753,  7.07623947, &
  &   6.60844053,  6.28791364,  6.07728703,  5.54643096,  5.80491167, &
  &   5.58415602,  5.41374528,  5.28497229,  5.22592821,  5.09817141, &
  &   6.12149689,  5.54083734,  5.06696878,  4.87005108,  4.59089647, &
  &   4.31176304,  9.55461698,  8.67396077,  7.97210197,  7.43439917, &
  &   6.58711862,  6.19536215,  6.01517290,  5.81623410,  5.65710424, &
  &   5.52640661,  5.44263305,  5.58285373,  7.02081898,  6.46815523, &
  &   5.98089120,  5.81686657,  5.53321815,  5.25477007, 11.02204549, &
  &  10.15679528,  9.35167836,  9.06926079,  8.97241155,  8.90092807, &
  &   8.85984840,  8.81736827,  8.79317710,  7.89969626,  8.80588454, &
  &   8.42439218,  8.54289262,  8.47583370,  8.45090888,  8.47339339, &
  &   7.83525634,  8.20702843,  7.70559063,  7.32755997,  7.03887381, &
  &   6.68978720,  6.05450052,  5.88752022,  5.70661499,  5.78450695, &
  &   7.79780729,  7.26443867,  6.78151984,  6.67883169,  6.39024318, &
  &   6.09527958, 11.79156076, 11.10997644,  9.51377795,  8.67197068, &
  &   8.77140725,  8.65402716,  8.53923501,  8.85024712 /)

   integer, dimension(86)      :: refn
   integer, dimension(7,86)    :: refc
   real(wp),dimension(7,86)    :: refq
   real(wp),dimension(7,86)    :: refh
   real(wp),dimension(7,86)    :: dftq,pbcq,gffq,solq
   real(wp),dimension(7,86)    :: dfth,pbch,gffh,solh
   real(wp),dimension(7,86)    :: hcount 
   real(wp),dimension(7,86)    :: ascale
   real(wp),dimension(7,86)    :: refcovcn
   real(wp),dimension(7,86)    :: refcn
   integer, dimension(7,86)    :: refsys 
   real(wp),dimension(23,7,86) :: alphaiw
   real(wp),dimension(23,7,86) :: refal
   real(wp),dimension(8)       :: secq
   real(wp),dimension(8)       :: dfts,pbcs,gffs,sols
   real(wp),dimension(8)       :: sscale
   real(wp),dimension(8)       :: seccn
   real(wp),dimension(8)       :: seccnd3
   real(wp),dimension(23,8)    :: secaiw

   include 'param_d4.fh'

contains

subroutine prmolc6(molc6,molc8,molpol,nat,at,  &
           &       cn,covcn,q,qlmom,c6ab,alpha,rvdw,hvol)
   use iso_fortran_env, only : id => output_unit
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_symbols, only: toSymbol 
!   use precision, only : wp => dp
   implicit none
   real(wp),intent(in)  :: molc6
   real(wp),intent(in)  :: molc8
   real(wp),intent(in)  :: molpol
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in),optional :: cn(nat)
   real(wp),intent(in),optional :: covcn(nat)
   real(wp),intent(in),optional :: q(nat)
   real(wp),intent(in),optional :: qlmom(3,nat)
   real(wp),intent(in),optional :: c6ab(nat,nat)
   real(wp),intent(in),optional :: alpha(nat)
   real(wp),intent(in),optional :: rvdw(nat)
   real(wp),intent(in),optional :: hvol(nat)
   real(wp),parameter :: autoaa = 0.52917726_wp
   integer :: i
   if(present(cn).or.present(covcn).or.present(q).or.present(c6ab) &
   &   .or.present(alpha).or.present(rvdw).or.present(hvol)) then
   write(id,'(a)')
   write(id,'(''   #   Z   '')',advance='no')
   if(present(cn))   write(id,'(''        CN'')',advance='no')
   if(present(covcn))write(id,'(''     covCN'')',advance='no')
   if(present(q))    write(id,'(''         q'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(s)'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(p)'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(d)'')',advance='no')
   if(present(c6ab)) write(id,'(''      C6AA'')',advance='no')
   if(present(alpha))write(id,'(''      α(0)'')',advance='no')
   if(present(rvdw)) write(id,'(''    RvdW/Å'')',advance='no')
   if(present(hvol)) write(id,'(''    relVol'')',advance='no')
!   write(*,'(a)')
   do i=1,nat
      write(*,'(i4,x,i3,x,a2)',advance='no') &
      &     i,at(i),toSymbol(at(i))
      if(present(cn))   write(id,'(f10.3)',advance='no')cn(i)
      if(present(covcn))write(id,'(f10.3)',advance='no')covcn(i)
      if(present(q))    write(id,'(f10.3)',advance='no')q(i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(1,i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(2,i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(3,i)
      if(present(c6ab)) write(id,'(f10.3)',advance='no')c6ab(i,i)
      if(present(alpha))write(id,'(f10.3)',advance='no')alpha(i)
      if(present(rvdw)) write(id,'(f10.3)',advance='no')rvdw(i)*autoaa
      if(present(hvol)) write(id,'(f10.3)',advance='no')hvol(i)
      write(*,'(a)')
   enddo
   endif
   write(id,'(/,x''Mol. C6AA /au*bohr^6 :''f18.6,'// &
   &         '/,x''Mol. C8AA /au*bohr^6 :''f18.6,'// &
   &         '/,x''Mol. a(0) /au        :''f18.6,/)') &
   &          molc6,molc8,molpol
end subroutine prmolc6

subroutine mdisp(nat,ndim,at,q,xyz,g_a,g_c, &
           &     gw,c6abns,molc6,molc8,molpol,aout,cout,rout,vout)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : thopi,gam, &
!  &  trapzd,zeta,r4r2, &
!  &  refn,refq,refal
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat) 
   real(wp),intent(in)  :: xyz(3,nat) 
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: molc6
   real(wp),intent(out) :: molc8
   real(wp),intent(out) :: molpol
   real(wp),intent(out),optional :: aout(23,nat)
   real(wp),intent(out),optional :: cout(nat,nat)
   real(wp),intent(out),optional :: rout(nat)
   real(wp),intent(out),optional :: vout(nat)

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: qmod,oth,iz
   real(wp),allocatable :: zetvec(:)
   real(wp),allocatable :: rvdw(:)
   real(wp),allocatable :: phv(:)
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   parameter (oth=1._wp/3._wp)
   
   allocate( zetvec(ndim),rvdw(nat),phv(nat),c6ab(nat,nat),aw(23,nat), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   molc6  = 0._wp
   molc8  = 0._wp
   molpol = 0._wp

   k = 0
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, refn(ia)
         k = k+1
         itbl(ii,i) = k
         zetvec(k) = gw(k) * zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
         aw(:,i) = aw(:,i) + zetvec(k) * refal(:,ii,ia)
      enddo
!     van-der-Waals radius, alpha = 4/3 pi r**3 <=> r = (3/(4pi) alpha)**(1/3)
      rvdw(i) = (0.25_wp*thopi*aw(1,i))**oth
!     pseudo-hirshfeld volume
      phv(i) = aw(1,i)/refal(1,1,ia)
      c6ab(i,i) = thopi * trapzd(aw(:,i)**2)
      molpol = molpol + aw(1,i)
      molc6  = molc6  + c6ab(i,i)
      molc8 = molc8 + 3*r4r2(ia)**2*c6ab(i,i)
      do j = 1, i-1
         ja = at(j)
         c6ab(j,i) = thopi * trapzd(aw(:,i)*aw(:,j))
         c6ab(i,j) = c6ab(j,i)
         molc6 = molc6 + 2*c6ab(j,i)
         molc8 = molc8 + 6*r4r2(ia)*r4r2(ja)*c6ab(j,i)
      enddo
   enddo

   if (present(aout)) aout = aw
   if (present(vout)) vout = phv
   if (present(rout)) rout = rvdw
   if (present(cout)) cout = c6ab

end subroutine mdisp

pure subroutine covncoord(nat,at,xyz,cn,cn_thr)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : k1,k4,k5,k6,rcov,en
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: cn_thr
   real(wp),intent(out) :: cn(nat)

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, rr, r2, xn

   cn = 0._wp

   do i=1,nat
      do j=1,i-1
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle 
         r=sqrt(r2)
!        covalent distance in bohr
         rco=rcov(at(j)) + rcov(at(i))
         rr=rco/r
         den=k4*exp(-(abs((en(at(i))-en(at(j))))+ k5)**2/k6 )
!        counting function exponential has a better long-range 
!        behavior than MHGs inverse damping
         cn(i)=cn(i)+den/(1.d0+exp(-k1*(rr-1.0d0)))
         cn(j)=cn(j)+den/(1.d0+exp(-k1*(rr-1.0d0)))
      enddo
   enddo

end subroutine covncoord

subroutine gncoord(nat,at,xyz,cn_thr,dcn)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : k1,k4,k5,k6,rcov
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: cn_thr
   real(wp),intent(out) :: dcn(3,nat,nat)
   integer  :: i, j
   real(wp) :: r, r2, rij(3)
   real(wp) :: rcovij
   real(wp) :: expterm
   real(wp) :: dcn_

   dcn = 0._wp

!  δCN/δrij = δ/δrij k4·exp(-(ΔEN+k5)²/k6)/(1+exp(-k1((ra+rb)/r-1)))
   do i = 1, nat
      do j = 1, i-1
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=(rcov(at(i))+rcov(at(j)))
         expterm=exp(-k1*(rcovij/r-1.0d0))
         dcn_= (-k1*rcovij*expterm)/ &
         &     (r2*((expterm+1.0d0)**2))
         dcn(:,i,i)= dcn_*rij/r + dcn(:,i,i)
         dcn(:,j,j)=-dcn_*rij/r + dcn(:,j,j)
         dcn(:,i,j)= dcn_*rij/r
         dcn(:,j,i)=-dcn_*rij/r
      enddo
   enddo

end subroutine gncoord

pure subroutine dist_r2(nat,xyz,r2)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
   implicit none
   integer, intent(in)  :: nat
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: r2(nat,nat)

   integer :: i,j

   r2 = 0._wp

   do i=1,nat
      do j=1,i-1
         r2(i,j) = sum( (xyz(:,i)-xyz(:,j))**2 )
         r2(j,i) = r2(i,j)
      enddo
   enddo

end subroutine dist_r2

pure subroutine dist_bj(oor6ab,oor8ab,r2,nat,at,a1,a2)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : r4r2
   implicit none
   integer, intent(in)  :: nat,at(nat)
   real(wp),intent(in)  :: a1,a2,r2(nat,nat)
   real(wp),intent(out) :: oor6ab(nat,nat),oor8ab(nat,nat)

   integer  :: i,j
   real(wp) :: cutoff

   oor6ab = 0._wp
   oor8ab = 0._wp

   do i=1,nat
      do j=1,i-1
         cutoff      = a1*sqrt(3._wp*r4r2(at(i))*r4r2(at(j)))+a2
         oor6ab(i,j) = 1._wp/(r2(i,j)**3 + cutoff**6)
         oor6ab(j,i) = oor6ab(i,j)
         oor8ab(i,j) = 1._wp/(r2(i,j)**4 + cutoff**8)
         oor8ab(j,i) = oor8ab(i,j)
      enddo
   enddo

end subroutine dist_bj

elemental function zeta(a,c,qref,qmod)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
   implicit none
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: zeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1._wp - exp( c * ( 1._wp - qref/qmod ) ) ) )
   endif

end function zeta

elemental function dzeta(a,c,qref,qmod)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : zeta
   implicit none
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: dzeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      dzeta = 0._wp
   else
      dzeta = - a * c * exp( c * ( 1._wp - qref/qmod ) ) &
      &           * zeta(a,c,qref,qmod) * qref / ( qmod**2 )
   endif

end function dzeta

pure function trapzd(pol)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
   implicit none
   real(wp),intent(in) :: pol(23)
   real(wp)            :: trapzd

   real(wp)            :: tmp1, tmp2
   real(wp),parameter  :: freq(23) = (/ &
&   0.000001_wp,0.050000_wp,0.100000_wp, &
&   0.200000_wp,0.300000_wp,0.400000_wp, &
&   0.500000_wp,0.600000_wp,0.700000_wp, &
&   0.800000_wp,0.900000_wp,1.000000_wp, &
&   1.200000_wp,1.400000_wp,1.600000_wp, &
&   1.800000_wp,2.000000_wp,2.500000_wp, &
&   3.000000_wp,4.000000_wp,5.000000_wp, &
&   7.500000_wp,10.00000_wp /)

!  do average between trap(1)-trap(22) .and. trap(2)-trap(23)
   tmp1 = 0.5_wp * ( &
&  ( freq (2) - freq (1) ) * ( pol (2) + pol (1) )+ &
&  ( freq (4) - freq (3) ) * ( pol (4) + pol (3) )+ &
&  ( freq (6) - freq (5) ) * ( pol (6) + pol (5) )+ &
&  ( freq (8) - freq (7) ) * ( pol (8) + pol (7) )+ &
&  ( freq(10) - freq (9) ) * ( pol(10) + pol (9) )+ &
&  ( freq(12) - freq(11) ) * ( pol(12) + pol(11) )+ &
&  ( freq(14) - freq(13) ) * ( pol(14) + pol(13) )+ &
&  ( freq(16) - freq(15) ) * ( pol(16) + pol(15) )+ &
&  ( freq(18) - freq(17) ) * ( pol(18) + pol(17) )+ &
&  ( freq(20) - freq(19) ) * ( pol(20) + pol(19) )+ &
&  ( freq(22) - freq(21) ) * ( pol(22) + pol(21) ))
   tmp2 = 0.5_wp * ( &
&  ( freq (3) - freq (2) ) * ( pol (3) + pol (2) )+ &
&  ( freq (5) - freq (4) ) * ( pol (5) + pol (4) )+ &
&  ( freq (7) - freq (6) ) * ( pol (7) + pol (6) )+ &
&  ( freq (9) - freq (8) ) * ( pol (9) + pol (8) )+ &
&  ( freq(11) - freq(10) ) * ( pol(11) + pol(10) )+ &
&  ( freq(13) - freq(12) ) * ( pol(13) + pol(12) )+ &
&  ( freq(15) - freq(14) ) * ( pol(15) + pol(14) )+ &
&  ( freq(17) - freq(16) ) * ( pol(17) + pol(16) )+ &
&  ( freq(19) - freq(18) ) * ( pol(19) + pol(18) )+ &
&  ( freq(21) - freq(20) ) * ( pol(21) + pol(20) )+ &
&  ( freq(23) - freq(22) ) * ( pol(23) + pol(22) ))

   trapzd = (tmp1+tmp2)

end function trapzd

elemental function cngw(wf,cn,cnref)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
   implicit none
   real(wp),intent(in) :: wf,cn,cnref
   real(wp)            :: cngw ! CN-gaussian-weight

   intrinsic :: exp

   cngw = exp ( -wf * ( cn - cnref )**2 )

end function cngw

elemental function dcngw(wf,cn,cnref)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : cngw
   implicit none
   real(wp),intent(in) :: wf,cn,cnref
   real(wp) :: dcngw

   dcngw = 2*wf*(cnref-cn)*cngw(wf,cn,cnref)

end function dcngw

subroutine d4dim(nat,at,g_a,g_c,mode,ndim)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : gam, &
!  &   zeta, &
!  &   ascale,hcount, &
!  &   alphaiw,refal, &
!  &   refsys, &
!  &   refcovcn,refcn, &
!  &   refq,refh, &
!  &   secq, &
!  &   sscale,secaiw, &
!  &   refn, refc
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: g_a,g_c
   integer, intent(in)  :: mode
   integer, intent(out) :: ndim

   integer  :: i,ia,is,icn,j
   integer  :: cncount(0:15)
   real(wp) :: sec_al(23),iz

   intrinsic :: nint

   select case(mode)
   case(1,2)
!     print'(x''* using PBE0/def2-TZVP Hirshfeld charges'')'
      refq = dftq
      refh = dfth
      secq = dfts
!  case(2)
!     refq = pbcq
!     refh = pbch
!     secq = pbcs
   case(3)
!     print'(x''* using classical Gasteiger charges'')'
      refq = gffq
      refh = gffh
      secq = gffs
   case(4)
!     print'(x''* using GFN2-xTB//GBSA(H2O) charges'')'
      refq = solq
      refh = solh
      secq = sols
   end select

   ndim = 0

!* set up refc und refal, also obtain the dimension of the dispmat
   do i = 1, nat
      cncount = 0
      cncount(0) = 1
      ia = at(i)
      do j = 1, refn(ia)
         is = refsys(j,ia)
         iz = zeff(is)
         sec_al = sscale(is)*secaiw(:,is) &
         &  * zeta(g_a,gam(is)*g_c,secq(is)+iz,refh(j,ia)+iz)
         icn =nint(refcn(j,ia))
         cncount(icn) = cncount(icn) + 1
         refal(:,j,ia) = max(ascale(j,ia)*(alphaiw(:,j,ia)-hcount(j,ia)*sec_al),0.0_wp)
      enddo
      do j = 1, refn(ia)
         icn = cncount(nint(refcn(j,ia)))
         refc(j,ia) = icn*(icn+1)/2
      enddo
      ndim = ndim + refn(ia)
   enddo

end subroutine d4dim

subroutine d4(nat,ndim,at,wf,g_a,g_c,covcn,gw,c6abns)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : thopi, &
!  &   trapzd,cngw, &
!  &   refn, &
!  &   refc,refal, &
!  &   refcovcn
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: wf,g_a,g_c
   real(wp),intent(in)  :: covcn(nat)
   real(wp),intent(out) :: gw(ndim)
   real(wp),intent(out) :: c6abns(ndim,ndim)

   integer  :: i,ia,is,icn,ii,iii,j,jj,ja,k,l
   integer,allocatable :: itbl(:,:)
   real(wp) :: twf,norm

   intrinsic :: maxval

   allocate( itbl(7,nat), source = 0 )

   gw = 0._wp
   c6abns = 0._wp

   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      norm = 0.0_wp
      do ii = 1, refn(ia)
         do iii = 1, refc(ii,ia)
            twf = iii*wf
            norm = norm + cngw(twf,covcn(i),refcovcn(ii,ia))
         enddo
      enddo
      norm = 1._wp / norm
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         do iii = 1, refc(ii,ia)
            twf = iii*wf
            gw(k) = gw(k) + cngw(twf,covcn(i),refcovcn(ii,ia)) * norm
         enddo
!    --- okay, if we run out of numerical precision, gw(k) will be NaN.
!        In case it is NaN, it will not match itself! So we can rescue
!        this exception. This can only happen for very high CNs.
         if (gw(k).ne.gw(k)) then
            if (maxval(refcovcn(:refn(ia),ia)).eq.refcovcn(ii,ia)) then
               gw(k) = 1.0_wp
            else
               gw(k) = 0.0_wp
            endif
         endif
         do j = 1, i-1
            ja = at(j)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6abns(l,k) = thopi * trapzd(refal(:,ii,ia)*refal(:,jj,ja))
               c6abns(k,l) = c6abns(l,k)
            enddo
         enddo
      enddo
   enddo

end subroutine d4

pure subroutine build_dispmat(nat,ndim,at,xyz,s6,s8,a1,a2,c6abns, &
                &             dispmat,r2out,oor6out,oor8out)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : r4r2, &
!  &  dist_r2,dist_bj, &
!  &  refn
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: s6,s8,a1,a2
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: dispmat(ndim,ndim)
   real(wp),intent(out),optional :: r2out(nat,nat)
   real(wp),intent(out),optional :: oor6out(nat,nat)
   real(wp),intent(out),optional :: oor8out(nat,nat)

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: c8abns
   real(wp),allocatable :: r2ab(:,:)
   real(wp),allocatable :: oor6ab(:,:)
   real(wp),allocatable :: oor8ab(:,:)

   intrinsic :: present

   allocate( r2ab(nat,nat), oor6ab(nat,nat), oor8ab(nat,nat), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   call dist_r2(nat,xyz,r2ab)
   call dist_bj(oor6ab,oor8ab,r2ab,nat,at,a1,a2)
 
   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         do j = 1, i-1
            ja = at(j)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c8abns = 3.0_wp * r4r2(ia)*r4r2(ja) * c6abns(k,l)
               dispmat(k,l) = &
               &  - s6 * ( c6abns(k,l) * oor6ab(i,j) ) &
               &  - s8 * ( c8abns      * oor8ab(i,j) )
               dispmat(l,k) = dispmat(k,l)
            enddo
         enddo
      enddo
   enddo

   if (present(r2out))   r2out = r2ab
   if (present(oor6out)) oor6out = oor6ab
   if (present(oor8out)) oor8out = oor8ab

end subroutine build_dispmat

pure subroutine build_wdispmat(nat,ndim,at,xyz,s6,s8,a1,a2,c6abns,gw, &
                &              wdispmat,r2out,oor6out,oor8out)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : r4r2, &
!  &  dist_r2,dist_bj, &
!  &  refn
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: s6,s8,a1,a2
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(out) :: wdispmat(ndim,ndim)
   real(wp),intent(out),optional :: r2out(nat,nat)
   real(wp),intent(out),optional :: oor6out(nat,nat)
   real(wp),intent(out),optional :: oor8out(nat,nat)

   integer :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: c8abns
   real(wp),allocatable :: r2ab(:,:)
   real(wp),allocatable :: oor6ab(:,:)
   real(wp),allocatable :: oor8ab(:,:)

   intrinsic :: present

   allocate( r2ab(nat,nat), oor6ab(nat,nat), oor8ab(nat,nat), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   call dist_r2(nat,xyz,r2ab)
   call dist_bj(oor6ab,oor8ab,r2ab,nat,at,a1,a2)
 
   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         do j = 1, i-1
            ja = at(j)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c8abns = 3.0_wp * r4r2(ia)*r4r2(ja) * c6abns(k,l)
               wdispmat(k,l) = gw(k)*gw(l) * ( &
               &  - s6 * ( c6abns(k,l) * oor6ab(i,j) ) &
               &  - s8 * ( c8abns      * oor8ab(i,j) ) )
               wdispmat(l,k) = wdispmat(k,l)
            enddo
         enddo
      enddo
   enddo

   if (present(r2out))   r2out = r2ab
   if (present(oor6out)) oor6out = oor6ab
   if (present(oor8out)) oor8out = oor8ab

end subroutine build_wdispmat

subroutine disppot(nat,ndim,at,q,g_a,g_c,wdispmat,gw,hdisp)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : gam, &
!  &  zeta,dzeta, &
!  &  refn, &
!  &  refq
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat)
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: wdispmat(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(out) :: hdisp(nat)

   integer  :: i,ii,k,ia
   real(wp) :: qmod,iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),zerovec(ndim),dumvec(ndim), source = 0._wp )

   hdisp   = 0.0_wp

   k = 0
   do i = 1, nat
       ia = at(i)
       iz = zeff(ia)
       do ii = 1, refn(ia)
          k = k + 1
          if (gw(k).lt.gw_cut) cycle
          zerovec(k) = dzeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
          zetavec(k) =  zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
      enddo
   enddo
!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim) 
   call dsymv('U',ndim,1._wp,wdispmat,ndim,zetavec,1,0._wp,dumvec,1)
!  call dgemv('N',ndim,ndim,1._wp,wdispmat,ndim,zetavec, &
!  &     1,0._wp,dumvec,1)
!  get atomic reference contributions
   k = 0
   do i = 1, nat
      ia = at(i)
      hdisp(i) = sum(dumvec(k+1:k+refn(ia))*zerovec(k+1:k+refn(ia)))
      k = k + refn(ia)
   enddo

   deallocate(zetavec,zerovec,dumvec)

end subroutine disppot

function edisp_scc(nat,ndim,at,q,g_a,g_c,wdispmat,gw) result(ed)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : gam, &
!  &  zeta,dzeta, &
!  &  refn, &
!  &  refq
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat)
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: wdispmat(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp) :: ed

   integer  :: i,ii,k,ia
   real(wp) :: qmod,iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),dumvec(ndim), source = 0._wp )

   ed = 0.0_wp

   k = 0
   do i = 1, nat
       ia = at(i)
       iz = zeff(ia)
       do ii = 1, refn(ia)
          k = k + 1
          if (gw(k).lt.gw_cut) cycle
          zetavec(k) =  zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
      enddo
   enddo
!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim) 
   call dsymv('U',ndim,0.5_wp,wdispmat,ndim,zetavec,1,0.0_wp,dumvec,1)
!  call dgemv('N',ndim,ndim,0.5_wp,wdispmat,ndim,zetavec, &
!  &           1,0.0_wp,dumvec,1)
   ed = dot_product(dumvec,zetavec)

   deallocate(zetavec,dumvec)

end function edisp_scc


subroutine edisp(nat,ndim,at,q,xyz,s6,s8,s9,a1,a2,alp,g_a,g_c, &
           &     gw,c6abns,mbd,E,aout)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : gam, &
!  &  zeta,build_dispmat,dispabc, &
!  &  refn,refq,refal
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat) 
   real(wp),intent(in)  :: xyz(3,nat) 
   real(wp),intent(in)  :: s6,s8,s9,a1,a2
   integer, intent(in)  :: alp
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(out) :: E
   real(wp),intent(out),optional :: aout(23,nat)

   integer  :: i,ii,ia,k,ij,l,j,jj,ja
   integer, allocatable :: itbl(:,:)
   real(wp) :: Embd,qmod,c6ij,c6ij_ns,oor6,oor8,r2,cutoff,iz
   real(wp),allocatable :: dispmat(:,:)
   real(wp),allocatable :: zetvec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp),allocatable :: dumvec(:)
   real(wp),allocatable :: c6ab(:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: oor6ab(:,:)

   intrinsic :: present,sqrt,sum
   
   allocate( zetvec(ndim),aw(23,nat),oor6ab(nat,nat), &
   &         zerovec(ndim),c6ab(nat*(nat+1)/2), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   e = 0.0_wp

   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo


   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         zetvec(k) = gw(k) * zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
         zerovec(k) = gw(k) * zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz)
         aw(:,i) = aw(:,i) + zerovec(k) * refal(:,ii,ia)
      enddo
   enddo

!$OMP parallel private(i,j,ia,ja,ij,k,l,r2,oor6,oor8,cutoff,c6ij,c6ij_ns) &
!$omp&         shared(c6ab) reduction(-:E)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
         r2 = sum( (xyz(:,i)-xyz(:,j))**2 )
         cutoff = a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+a2
         oor6 = 1._wp/(r2**3 + cutoff**6)
         oor6ab(i,j) = oor6
         oor6ab(j,i) = oor6
         oor8 = 1._wp/(r2**4 + cutoff**8)
         c6ij_ns = 0.0_wp
         c6ij = 0.0_wp
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij_ns = c6ij_ns + zerovec(k)*zerovec(l)*c6abns(k,l)
               c6ij = c6ij + zetvec(k)*zetvec(l)*c6abns(k,l)
            enddo
         enddo
         c6ab(ij) = c6ij_ns
         E = E - c6ij*(s6*oor6 + s8*3._wp*r4r2(ia)*r4r2(ja)*oor8)
      enddo
   enddo
!$omp enddo
!$omp end parallel

   select case(mbd)
   case(1) ! full RPA-like MBD
!     print'(x,''* MBD effects calculated by RPA like scheme'')'
      call dispmb(Embd,aw,xyz,oor6ab,nat)
      E = E + s9*Embd
   case(2) ! Axilrod-Teller-Muto three-body term
!     print'(x,''* MBD effects calculated by ATM formula'')'
      call dispabc(nat,at,xyz,aw,s9,a1,a2,alp,Embd)
      E = E + Embd
   case(3) ! D3-like approximated ATM term
!     print'(x,''* MBD effects approximated by ATM formula'')'
      call apprabc(nat,at,xyz,c6ab,s9,a1,a2,alp,Embd)
      E = E + Embd
   case default
      Embd = 0.0_wp
   end select

   if (present(aout)) then
      aout = 0._wp
      do i = 1, nat
         ia = at(i)
         do ii = 1, refn(ia)
            aout(:,i) = aout(:,i) + zetvec(k) * refal(:,ii,ia)
         enddo
      enddo
   endif

end subroutine edisp

subroutine gcovncoord(nat,at,xyz,cn_thr,cn,dcn)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : k1,k4,k5,k6,rcov,en
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: cn_thr
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(out) :: dcn(nat,nat)
   integer  :: i,j
   real(wp) :: r, r2, rij(3)
   real(wp) :: den
   real(wp) :: rcovij
   real(wp) :: expterm

   intrinsic :: sum,sqrt,exp,abs

   cn = 0._wp
   dcn = 0._wp

!  ∂CN/∂rij = ∂/∂rij k4·exp(-(∂EN+k5)²/k6)/(1+exp(-k1((ra+rb)/r-1)))
   do i = 1, nat
      do j = 1, i-1
         rij = xyz(:,j) - xyz(:,i)
         r2 = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         den=k4*exp(-(abs((en(at(i))-en(at(j))))+ k5)**2/k6 )
         rcovij=(rcov(at(i))+rcov(at(j)))
         expterm=exp(-k1*(rcovij/r-1._wp))
         cn(i) = cn(i) + den/(1._wp + expterm)
         cn(j) = cn(j) + den/(1._wp + expterm)
         dcn(i,j)=(-k1*rcovij*expterm*den)/ &
         &     (r2*((1._wp + expterm)**2))
         dcn(j,i)=dcn(i,j)
      enddo
   enddo

end subroutine gcovncoord


!* compute D4 gradient
!  ∂E/∂rij = ∂/∂rij (W·D·W)
!          = ∂W/∂rij·D·W  + W·∂D/∂rij·W + W·D·∂W/∂rij
!  ∂W/∂rij = ∂(ζ·w)/∂rij = ζ·∂w/∂rij = ζ·∂w/∂CN·∂CN/∂rij
!  ∂ζ/∂rij = 0
subroutine dispgrad(nat,ndim,at,q,xyz, &
           &        s6,s8,s9,a1,a2,alp,wf,g_a,g_c, &
           &        c6abns,mbd, &
           &        g,eout,aout)
!   use precision, only : wp => dp
   use xtb_mctc_accuracy, only: wp
!  use dftd4, only : &
!  &  r4r2,gam,en,rcov, &
!  &  zeta,cngw,dabcgrad, &
!  &  refn, refc, refal, &
!  &  refq, refcovcn
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat)
!  real(wp),intent(in)  :: cn(nat) ! calculate on-the-fly
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: s6,s8,s9,a1,a2
   integer, intent(in)  :: alp
   real(wp),intent(in)  :: wf,g_a,g_c
!  real(wp),intent(in)  :: gw(ndim) ! calculate on-the-fly
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(inout)        :: g(3,nat)
   real(wp),intent(out),optional :: eout
   real(wp),intent(out),optional :: aout(23,nat)


   integer  :: i,ii,iii,j,jj,k,l,ia,ja,ij
   integer, allocatable :: itbl(:,:)
   real(wp) :: iz
   real(wp) :: qmod,eabc,ed
   real(wp) :: norm,dnorm
   real(wp) :: dexpw,expw
   real(wp) :: twf,tgw,r4r2ij
   real(wp) :: rij(3),r,r2,r4,r6,r8,R0
   real(wp) :: oor6,oor8,door6,door8
   real(wp) :: c8abns,disp,x1,x2,x3
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: rcovij,expterm,den,dcndr
   real(wp) :: drdx(3),dtmp
   real(wp),allocatable :: r2ab(:)
   real(wp),allocatable :: dc6dr(:)
   real(wp),allocatable :: dc6dcn(:)
   real(wp),allocatable :: zvec(:)
   real(wp),allocatable :: dzvec(:)
   real(wp),allocatable :: gw(:)
   real(wp),allocatable :: dgw(:)
   real(wp),allocatable :: cn(:)
   real(wp) :: cn_thr,r_thr,gw_thr
   parameter(cn_thr = 1600.0_wp)
   parameter( r_thr=5000._wp)
   parameter(gw_thr=0.000001_wp)
!  timing
!  real(wp) :: time0,time1
!  real(wp) :: wall0,wall1

   intrinsic :: present,sqrt,sum,maxval,exp,abs

   allocate( dc6dr(nat*(nat+1)/2),dc6dcn(nat),  &
   &         r2ab(nat*(nat+1)/2),cn(nat),  &
   &         zvec(ndim),dzvec(ndim),  &
   &         gw(ndim),dgw(ndim),  &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   ed = 0.0_wp

!  derivatives of the covalent coordination number
!  ∂CN/∂rij = ∂/∂rij k4·exp(-(∂EN+k5)²/(2·k6²))/(1+exp(-k1(-1)))
   call covncoord(nat,at,xyz,cn,cn_thr)

!  precalc
   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

!$OMP parallel private(i,ii,iii,ia,iz,k,norm,dnorm,twf,tgw,dexpw,expw)  &
!$omp&         shared (gw,dgw,zvec,dzvec)
!$omp do
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      norm  = 0.0_wp
      dnorm = 0.0_wp
      do ii=1,refn(ia)
         do iii = 1, refc(ii,ia)
            twf = iii*wf
            tgw = cngw(twf,cn(i),refcovcn(ii,ia))
            norm  =  norm + tgw
            dnorm = dnorm + 2*twf*(refcovcn(ii,ia)-cn(i))*tgw
         enddo
      enddo
      norm = 1._wp/norm
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         dexpw=0.0_wp
         expw=0.0_wp
         do iii = 1, refc(ii,ia)
            twf = wf*iii
            tgw = cngw(twf,cn(i),refcovcn(ii,ia))
            expw  =  expw + tgw
            dexpw = dexpw + 2*twf*(refcovcn(ii,ia)-cn(i))*tgw
         enddo

         ! save
         gw(k) = expw*norm
         if (gw(k).ne.gw(k)) then
            if (maxval(refcovcn(:refn(ia),ia)).eq.refcovcn(ii,ia)) then
               gw(k) = 1.0_wp
            else
               gw(k) = 0.0_wp
            endif
         endif
         zvec(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz) * gw(k)
         ! NEW: q=0 for ATM
         gw(k) =  zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz) * gw(k)

         dgw(k) = dexpw*norm-expw*dnorm*norm**2
         if (dgw(k).ne.dgw(k)) then
            dgw(k) = 0.0_wp
         endif
         dzvec(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz) * dgw(k)
         ! NEW: q=0 for ATM
         dgw(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz) * dgw(k)
      enddo
   enddo
!$omp end do
!$omp end parallel

!$OMP parallel private(i,j,ia,ja,ij,k,l,c6ij,dic6ij,djc6ij,disp,  &
!$omp&                 rij,r2,r,r4r2ij,r0,oor6,oor8,door6,door8)  &
!$omp&         shared(r2ab) reduction(+:dc6dr,dc6dcn) reduction(-:ed)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
         rij = xyz(:,j) - xyz(:,i)
         r2 = sum( rij**2 )
         r2ab(ij) = r2
         if (r2.gt.r_thr) cycle
         ! temps
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*c6abns(k,l)
               dic6ij = dic6ij + dzvec(k)*zvec(l)*c6abns(k,l)
               djc6ij = djc6ij + zvec(k)*dzvec(l)*c6abns(k,l)
            enddo
         enddo

         r = sqrt(r2)

         r4r2ij = 3*r4r2(ia)*r4r2(ja)
         r0 = a1*sqrt(r4r2ij) + a2

         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         door6 = -6*r2**2*r*oor6**2
         door8 = -8*r2**3*r*oor8**2

         disp = s6*oor6 + s8*r4r2ij*oor8

         ! save
         dc6dcn(i) = dc6dcn(i) + dic6ij*disp
         dc6dcn(j) = dc6dcn(j) + djc6ij*disp
         dc6dr(ij) = dc6dr(ij) + c6ij*(s6*door6 + s8*r4r2ij*door8)

         ed = ed - c6ij*disp
      enddo
   enddo
!$omp enddo
!$omp end parallel

   select case(mbd)
!  case(1) ! full RPA-like MBD
!     print'(x,''* MBD effects calculated by RPA like scheme'')'
!     call raise('W','MBD gradient not fully implemented yet')
!     call mbdgrad(nat,xyz,aw,daw,oor6ab,g,embd)
   case(1,2) ! Axilrod-Teller-Muto three-body term
      if(mbd.eq.1) then
!        call raise('W','MBD gradient not fully implemented yet')
!        print'(''MBD gradient not fully implemented yet'')'
!        print'(x,''* calculate MBD effects with ATM formula instead'')'
      else
!        print'(x,''* MBD effects calculated by ATM formula'')'
      endif
!     call dabcgrad(nat,ndim,at,xyz,s9,a1,a2,alp,dcn,gw,dgw,itbl,g,embd)
   case(3) ! D3-like approximated ATM term
!     print'(x,''* MBD effects approximated by ATM formula'')'
      call dabcappr(nat,ndim,at,xyz,s9,a1,a2,alp,  &
           &        r2ab,gw,dgw,c6abns,itbl,dc6dr,dc6dcn,eabc)
   end select

!$OMP parallel private(i,j,ia,ja,ij,rij,r2,r,drdx,den,rcovij,  &
!$omp&                 expterm,dcndr,dtmp) reduction(+:g)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j

         rij = xyz(:,j) - xyz(:,i)
         r2 = sum( rij**2 )
         r = sqrt(r2)
         drdx = rij/r

         if (r2.lt.cn_thr) then
            den=k4*exp(-(abs((en(ia)-en(ja)))+ k5)**2/k6 )
            rcovij=rcov(ia)+rcov(ja)
            expterm=exp(-k1*(rcovij/r-1._wp))
            dcndr=(-k1*rcovij*expterm*den)/(r2*((1._wp + expterm)**2))
         else
            dcndr=0.0_wp
         endif

         dtmp = dc6dr(ij) + (dc6dcn(i)+dc6dcn(j))*dcndr

         g(:,i) = g(:,i) + dtmp * drdx
         g(:,j) = g(:,j) - dtmp * drdx
      enddo
   enddo
!$omp enddo
!$omp end parallel

!  print*,ed,eabc

   if (present(eout)) eout = ed + eabc

   if (present(aout)) then
      aout = 0._wp
      do i = 1, nat
         ia = at(i)
         do ii = 1, refn(ia)
            aout(:,i) = aout(:,i) + zvec(k) * refal(:,ii,ia)
         enddo
      enddo
   endif


end subroutine dispgrad

subroutine apprabc(nat,at,xyz,c6ab,s9,a1,a2,alp,E)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : &
!  &  r4r2, &
!  &  thopi, &
!  &  trapzd
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: c6ab(nat*(nat+1)/2)
   real(wp),intent(in)  :: s9,a1,a2
   integer, intent(in)  :: alp
   real(wp),intent(out) :: E

   integer  :: i,j,k,ia,ja,ka,ij,ik,jk
   real(wp) :: rij(3),rjk(3),rik(3),r2ij,r2jk,r2ik,cij,cjk,cik,cijk
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk,fdmp
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
!  parameter(thf = 3._wp/4._wp)

   intrinsic :: sum,sqrt

   E = 0.0_wp

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
         rij  = xyz(:,j) - xyz(:,i)
         r2ij = sum(rij**2)
         cij  = a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+a2
         do k = 1, j-1
            ka = at(k)
            ik = i*(i-1)/2 + k
            jk = j*(j-1)/2 + k
            rik   = xyz(:,i) - xyz(:,k)
            r2ik  = sum(rik**2)
            cik   = a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+a2
            rjk   = xyz(:,k) - xyz(:,j)
            r2jk  = sum(rjk**2)
            cjk   = a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+a2
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            c9ijk = s9*sqrt(c6ab(ij)*c6ab(jk)*c6ab(ik))
            atm = ( 0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp
            fdmp = one/(one+six*((cijk/rijk)**oth)**alp)
            oor9ijk = atm/rijk**3*fdmp
            E = E + c9ijk * oor9ijk
         enddo
      enddo
   enddo

end subroutine apprabc


subroutine dispabc(nat,at,xyz,aw,s9,a1,a2,alp,E)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : &
!  &  r4r2, &
!  &  thopi, &
!  &  trapzd
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: aw(23,nat)
   real(wp),intent(in)  :: s9,a1,a2
   integer, intent(in)  :: alp
   real(wp),intent(out) :: E

   integer  :: i,j,k,ia,ja,ka
   real(wp) :: rij(3),rjk(3),rik(3),r2ij,r2jk,r2ik,cij,cjk,cik,cijk
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk,fdmp
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
   parameter(thf = 3._wp/4._wp)

   intrinsic :: sum,sqrt

   E = 0.0_wp

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         rij  = xyz(:,j) - xyz(:,i)
         r2ij = sum(rij**2)
         cij  = (a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+a2)
         do k = 1, j-1
            ka = at(k)
            rik   = xyz(:,i) - xyz(:,k)
            r2ik  = sum(rik**2)
            cik   = (a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+a2)
            rjk   = xyz(:,k) - xyz(:,j)
            r2jk  = sum(rjk**2)
            cjk   = (a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+a2)
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            c9ijk = s9*thopi*trapzd( aw(:,i)*aw(:,j)*aw(:,k) )
            atm = ( 0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp
            fdmp = one/(one+six*(thf*(cijk/rijk)**oth)**alp)
            oor9ijk = atm/rijk**3*fdmp
            E = E + c9ijk * oor9ijk
         enddo
      enddo
   enddo

end subroutine dispabc

subroutine abcappr(nat,ndim,at,xyz,g_a,g_c,s9,a1,a2,alp,gw,r2ab,c6abns,eabc)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : r4r2, &
!  &  trapzd,thopi, &
!  &  refn,refal
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: s9,a1,a2
   integer, intent(in)  :: alp
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: r2ab(nat*(nat+1)/2)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: eabc

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   integer  :: ij,jk,ik
   integer, allocatable :: itbl(:,:)
   real(wp),allocatable :: c6ab(:),zvec(:),c(:)
   real(wp) :: r2ij,r2jk,r2ik,iz
   real(wp) :: cij,cjk,cik,cijk
   real(wp) :: fdmp,dtmp,oor9tmp,c9tmp
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk
   real(wp) :: drij(3),drjk(3),drik(3)
   real(wp) :: oorij,oorjk,oorik
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: dic9ijk,djc9ijk,dkc9ijk
   real(wp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(wp) :: dum1,dum2,dum3
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
!  parameter(thf = 3._wp/4._wp)
   real(wp) :: r_thr,gw_thr
   parameter( r_thr=1600._wp)
   parameter(gw_thr=0.0001_wp)

   intrinsic :: sqrt

   allocate( c6ab(nat*(nat+1)/2),zvec(ndim),c(nat*(nat+1)/2),  &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   eabc = 0.0_wp

!  precalc
   k = 0
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, refn(ia)
         k = k+1
         itbl(ii,i) = k
         ! NEW: q=0 for ATM
         zvec(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz) * gw(k)
      enddo
   enddo

!$OMP parallel private(i,ia,j,ja,ij,r2ij,c6ij)  &
!$omp&         shared (c6ab,c)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         c(ij) = a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+a2
         if(r2ij.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*c6abns(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

!$OMP parallel private(i,j,ij,ia,ja,k,ka,ik,jk,atm,fdmp,  &
!$omp&                 r2ij,cij,r2ik,r2jk,cik,cjk,r2ijk,rijk,cijk, &
!$omp&                 c9ijk,oor9ijk) &
!$omp&         reduction(+:eabc)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         if(r2ij.gt.r_thr) cycle
         cij  = c(ij)
         do k = 1, j-1
            ka = at(k)
            ik = i*(i-1)/2 + k
            jk = j*(j-1)/2 + k
            r2ik  = r2ab(ik)
            r2jk  = r2ab(jk)
            if((r2ik.gt.r_thr).or.(r2jk.gt.r_thr)) cycle
            cik   = c(ik)
            cjk   = c(jk)
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik

            atm = ((0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp)/(rijk**3)

            fdmp = one/(one+six*((cijk/rijk)**oth)**alp)

            c9ijk = s9*sqrt(c6ab(ij)*c6ab(jk)*c6ab(ik))

            oor9ijk = atm*fdmp
            eabc = eabc + c9ijk*oor9ijk

         enddo ! k/C
      enddo ! j/B
   enddo ! i/A
!$omp enddo
!$omp end parallel

   deallocate( c6ab,c,zvec )

end subroutine abcappr

subroutine dabcappr(nat,ndim,at,xyz,s9,a1,a2,alp,  &
                &        r2ab,zvec,dzvec,c6abns,itbl,dc6dr,dc6dcn,eout)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : r4r2, &
!  &  trapzd,thopi, &
!  &  refn,refal
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: s9,a1,a2
   integer, intent(in)  :: alp
   real(wp),intent(in)  :: r2ab(nat*(nat+1)/2)
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: dzvec(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: itbl(7,nat)
   real(wp),intent(inout)        :: dc6dr(nat*(nat+1)/2)
   real(wp),intent(inout)        :: dc6dcn(nat)
   real(wp),intent(out),optional :: eout

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   integer  :: ij,jk,ik
   real(wp),allocatable :: c6ab(:),dc6ab(:,:)
   real(wp) :: r2ij,r2jk,r2ik
   real(wp) :: cij,cjk,cik,cijk
   real(wp) :: fdmp,dtmp,oor9tmp,c9tmp
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk
   real(wp) :: drij(3),drjk(3),drik(3)
   real(wp) :: oorij,oorjk,oorik
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: dic9ijk,djc9ijk,dkc9ijk
   real(wp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(wp) :: eabc
   real(wp) :: dum1,dum2,dum3
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
!  parameter(thf = 3._wp/4._wp)
   real(wp) :: r_thr,gw_thr
   parameter( r_thr=1600._wp)
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt

   allocate( c6ab(nat*(nat+1)/2),dc6ab(nat,nat),  &
   &         source = 0.0_wp )

   eabc = 0.0_wp

!$OMP parallel private(i,ia,j,ja,ij,r2ij,c6ij,dic6ij,djc6ij)  &
!$omp&         shared (c6ab,dc6ab)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         if(r2ij.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*c6abns(k,l)
               dic6ij = dic6ij + dzvec(k)*zvec(l)*c6abns(k,l)
               djc6ij = djc6ij + zvec(k)*dzvec(l)*c6abns(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
         dc6ab(i,j) = dic6ij
         dc6ab(j,i) = djc6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

!$OMP parallel private(i,j,ij,ia,ja,k,ka,ik,jk,oorjk,oorik,atm,fdmp,  &
!$omp&                 r2ij,cij,oorij,r2ik,r2jk,cik,cjk,r2ijk,rijk,cijk, &
!$omp&                 dijatm,djkatm,dikatm,dtmp,dijfdmp,djkfdmp,dikfdmp,  &
!$omp&                 c9ijk,oor9ijk,dic9ijk,djc9ijk,dkc9ijk) &
!$omp&         reduction(+:eabc,dc6dr) reduction(-:dc6dcn)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         if(r2ij.gt.r_thr) cycle
         cij  = a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+a2
         oorij = 1._wp/sqrt(r2ij)
         do k = 1, j-1
            ka = at(k)
            ik = i*(i-1)/2 + k
            jk = j*(j-1)/2 + k
            r2ik  = r2ab(ik)
            r2jk  = r2ab(jk)
            if((r2ik.gt.r_thr).or.(r2jk.gt.r_thr)) cycle
            cik   = a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+a2
            cjk   = a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+a2
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            oorjk = 1._wp/sqrt(r2jk)
            oorik = 1._wp/sqrt(r2ik)

            atm = ((0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp)/(rijk**3)
            dijatm=-0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
            &      +r2ij*(3._wp*r2jk**2+2._wp*r2jk*r2ik+3._wp*r2ik**2) &
            &      -5._wp*(r2jk-r2ik)**2*(r2jk+r2ik)) &
            &      /(r2ijk*rijk**3)*oorij
            djkatm=-0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
            &      +r2jk*(3._wp*r2ik**2+2._wp*r2ik*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2ik-r2ij)**2*(r2ik+r2ij)) &
            &      /(r2ijk*rijk**3)*oorjk
            dikatm=-0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
            &      +r2ik*(3._wp*r2jk**2+2._wp*r2jk*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2jk-r2ij)**2*(r2jk+r2ij)) &
            &      /(r2ijk*rijk**3)*oorik

            fdmp = one/(one+six*((cijk/rijk)**oth)**alp)
            dtmp = -(oth*six*alp*((cijk/rijk)**oth)**alp)*fdmp**2
            dijfdmp = dtmp*oorij
            djkfdmp = dtmp*oorjk
            dikfdmp = dtmp*oorik

            c9ijk = s9*sqrt(c6ab(ij)*c6ab(jk)*c6ab(ik))

            oor9ijk = atm*fdmp
            eabc = eabc + c9ijk*oor9ijk

            dc6dr(ij) = dc6dr(ij) + (atm*dijfdmp - dijatm*fdmp)*c9ijk
            dc6dr(ik) = dc6dr(ik) + (atm*dikfdmp - dikatm*fdmp)*c9ijk
            dc6dr(jk) = dc6dr(jk) + (atm*djkfdmp - djkatm*fdmp)*c9ijk
            dic9ijk = dc6ab(i,j)/c6ab(ij) + dc6ab(i,k)/c6ab(ik)
            djc9ijk = dc6ab(j,i)/c6ab(ij) + dc6ab(j,k)/c6ab(jk)
            dkc9ijk = dc6ab(k,j)/c6ab(jk) + dc6ab(k,i)/c6ab(ik)
            dc6dcn(i) = dc6dcn(i) - 0.5_wp*c9ijk*oor9ijk*dic9ijk
            dc6dcn(j) = dc6dcn(j) - 0.5_wp*c9ijk*oor9ijk*djc9ijk
            dc6dcn(k) = dc6dcn(k) - 0.5_wp*c9ijk*oor9ijk*dkc9ijk

         enddo ! k/C
      enddo ! j/B
   enddo ! i/A
!$omp enddo
!$omp end parallel

   if (present(eout)) eout=eabc

end subroutine dabcappr


!* here is the theory for the ATM-gradient (SAW, 180224)
! EABC = WA·WB·WC·DABC
! ∂EABC/∂X = ∂/∂X(WA·WB·WC·DABC)
!          = ∂WA/∂X·WB·WC·DABC + WA·∂WB/∂X·WC·DABC + WA·WB·∂WC/∂X·WC·DABC
!            + WA·WB·WC·∂DABC/∂X
! ∂/∂X =  ∂rAB/∂X·∂/∂rAB +  ∂rBC/∂X·∂/∂rBC +  ∂rCA/∂X·∂/∂rCA
!      = (δAX-δBX)∂/∂rAB + (δBX-δCX)∂/∂rBC + (δCX-δAX)∂/∂rCA
! ∂EABC/∂A = ∑A,ref ∑B,ref ∑C,ref
!            + (∂WA/∂rAB-∂WA/∂rCA)·WB·WC·DABC
!            + WA·∂WB/∂rAB·WC·DABC
!            - WA·WB·∂WC/∂rCA·DABC
!            + WA·WB·WC·(∂DABC/∂rAB-∂DABC/∂rCA)
! ∂EABC/∂B = ∑A,ref ∑B,ref ∑C,ref
!            - ∂WA/∂rAB·WB·WC·DABC
!            + WA·(∂WB/∂rBC-∂WB/∂rAB)·WC·DABC
!            + WA·WB·∂WC/∂rBC·DABC
!            + WA·WB·WC·(∂DABC/∂rBC-∂DABC/∂rAB)
! ∂EABC/∂C = ∑A,ref ∑B,ref ∑C,ref
!            + ∂WA/∂rCA·WB·WC·DABC
!            - WA·∂WB/∂rBC·WC·DABC
!            + WA·WB·(∂WC/∂rCA-∂WC/∂rBC)·DABC
!            + WA·WB·WC·(∂DABC/∂rCA-∂DABC/∂rBC)
! ∂WA/∂rAB = ∂CNA/∂rAB·∂WA/∂CNA w/ ζ=1 and WA=wA
! ATM = 3·cos(α)cos(β)cos(γ)+1
!     = 3/8(r²AB+r²BC-r²CA)(r²AB+r²CA-r²BC)(r²BC+r²CA-r²AB)/(r²BC·r²CA·r²AB)+1
! ∂ATM/∂rAB = 3/4(2r⁶AB-r⁶BC-r⁶CA-r⁴AB·r²BC-r⁴AB·r²CA+r⁴BC·r²CA+r²BC·r⁴CA)
!             /(r³AB·r²BC·r²CA)
! DABC = C9ABCns·f·ATM/(rAB·rBC·rCA)³
! f = 1/(1+6(¾·∛[RAB·RBC·RCA/(rAB·rBC·rCA)])¹⁶)
! ∂(f/(r³AB·r³BC·r³CA)/∂rAB = 
!   ⅓·((6·(16-9)·(¾·∛[RAB·RBC·RCA/(rAB·rBC·rCA)])¹⁶-9)·f²/(r⁴AB·r³BC·r³CA)
subroutine dabcgrad(nat,ndim,at,xyz,s9,a1,a2,alp,dcn,zvec,dzvec,itbl,g,eout)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : r4r2, &
!  &  trapzd,thopi, &
!  &  refn,refal
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: s9,a1,a2
   integer, intent(in)  :: alp
   real(wp),intent(in)  :: dcn(nat,nat)
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: dzvec(ndim)
   integer, intent(in)  :: itbl(7,nat)
   real(wp),intent(inout)        :: g(3,nat)
   real(wp),intent(out),optional :: eout

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   real(wp) :: rij(3),rjk(3),rik(3)
   real(wp) :: r2ij,r2jk,r2ik
   real(wp) :: cij,cjk,cik,cijk
   real(wp) :: fdmp,dtmp,oor9tmp,c9tmp
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk
   real(wp) :: drij(3),drjk(3),drik(3)
   real(wp) :: oorij,oorjk,oorik
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(wp) :: eabc
   real(wp) :: dum1,dum2,dum3
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
   parameter(thf = 3._wp/4._wp)
   real(wp) :: r_thr,gw_thr
   parameter( r_thr=1600._wp)
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt,sum

   eabc = 0._wp

!$omp parallel private(ia,ja,ka,l,m,n, &
!$omp&         rij,rjk,rik,r2ij,r2jk,r2ik, &
!$omp&         cij,cjk,cik,cijk, &
!$omp&         fdmp,dtmp,oor9tmp,c9tmp, &
!$omp&         atm,r2ijk,c9ijk,oor9ijk,rijk, &
!$omp&         drij,drjk,drik,oorij,oorjk,oorik, &
!$omp&         dijfdmp,dikfdmp,djkfdmp, &
!$omp&         dijatm,dikatm,djkatm, &
!$omp&         dijoor9ijk,djkoor9ijk,dikoor9ijk, &
!$omp&         x1,x2,x3,x4,x5,x6,x7,x8,x9) &
!$omp&         reduction(+:g,eabc)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = at(j)
!    --- all distances, cutoff radii ---
         rij  = xyz(:,j) - xyz(:,i)
         r2ij = sum(rij**2)
         if(r2ij.gt.r_thr) cycle
         cij  = (a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+a2)
         oorij = 1._wp/sqrt(r2ij)
         do k = 1, j-1
!           if(k.eq.j) cycle
!           if(i.eq.k) cycle
            ka = at(k)
            rik   = xyz(:,i) - xyz(:,k)
            rjk   = xyz(:,k) - xyz(:,j)
            r2ik  = sum(rik**2)
            r2jk  = sum(rjk**2)
            if((r2ik.gt.r_thr).or.(r2jk.gt.r_thr)) cycle
            cik   = (a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+a2)
            cjk   = (a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+a2)
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            oorjk = 1._wp/sqrt(r2jk)
            oorik = 1._wp/sqrt(r2ik)

            x2 = 0._wp
            x4 = 0._wp
            x6 = 0._wp
            c9ijk = 0._wp

!       --- sum up all references ---
            do ii = 1, refn(ia) ! refs of A
               l = itbl(ii,i)
               do jj = 1, refn(ja) ! refs of B
                  m = itbl(jj,j)
                  do kk = 1, refn(ka) ! refs of C
                     n = itbl(kk,k)
                     if ((zvec(l)*zvec(m)*zvec(n)).lt.gw_thr) cycle
                     c9tmp = s9*thopi*trapzd(refal(:,ii,ia)*refal(:,jj,ja) &
                     &                      *refal(:,kk,ka))

                     c9ijk = c9ijk + c9tmp*zvec(n)*zvec(m)*zvec(l)
!                --- intermediates ---
!                    ∂WA/∂CNA·WB·WC
                     x2 = x2 - dzvec(l)*zvec(m)*zvec(n)*c9tmp
!                    WA·∂WB/∂CNB·WC
                     x4 = x4 - dzvec(m)*zvec(l)*zvec(n)*c9tmp
!                    WA·WB·∂WC/∂CNC
                     x6 = x6 - dzvec(n)*zvec(m)*zvec(l)*c9tmp

                  enddo ! refs of k/C
               enddo ! refs of j/B
            enddo ! refs of i/A

!       --- geometrical term and r⁻³AB·r⁻³BC·r⁻³CA ---
!           ATM = 3·cos(α)cos(β)cos(γ)+1
!               = 3/8(r²AB+r²BC-r²CA)(r²AB+r²CA-r²BC)(r²BC+r²CA-r²AB)
!                 /(r²BC·r²CA·r²AB)+1
            atm = ((0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp)/(rijk**3)
            dijatm=-0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
            &      +r2ij*(3._wp*r2jk**2+2._wp*r2jk*r2ik+3._wp*r2ik**2) &
            &      -5._wp*(r2jk-r2ik)**2*(r2jk+r2ik)) &
            &      /(r2ijk*rijk**3)*oorij
            djkatm=-0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
            &      +r2jk*(3._wp*r2ik**2+2._wp*r2ik*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2ik-r2ij)**2*(r2ik+r2ij)) &
            &      /(r2ijk*rijk**3)*oorjk
            dikatm=-0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
            &      +r2ik*(3._wp*r2jk**2+2._wp*r2jk*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2jk-r2ij)**2*(r2jk+r2ij)) &
            &      /(r2ijk*rijk**3)*oorik

!       --- damping function ---
!           1/(1+6(¾·∛[RAB·RBC·RCA/(rAB·rBC·rCA)])¹⁶)
            fdmp = one/(one+six*(thf*(cijk/rijk)**oth)**alp)
            dtmp = -(oth*six*alp*(thf*(cijk/rijk)**oth)**alp)*fdmp**2
            dijfdmp = dtmp*oorij
            djkfdmp = dtmp*oorjk
            dikfdmp = dtmp*oorik

!       --- intermediates ---
!           ∂WA/∂rAB·WB·WC·DABC = ∂CNA/∂rAB·(∂WA/∂CNA·WB·WC)·DABC
            x1 = x2*dcn(i,j)*( atm*fdmp )
!           ∂WA/∂rCA·WB·WC·DABC = ∂CNA/∂rCA·(∂WA/∂CNA·WB·WC)·DABC
            x2 = x2*dcn(i,k)*( atm*fdmp )
!           WA·∂WB/∂rBC·WC·DABC = ∂CNB/∂rBC·(WA·∂WB/∂rBC·WC)·DABC
            x3 = x4*dcn(j,k)*( atm*fdmp )
!           WA·∂WB/∂rAB·WC·DABC = ∂CNB/∂rAB·(WA·∂WB/∂rAB·WC)·DABC
            x4 = x4*dcn(i,j)*( atm*fdmp )
!           WA·WB·∂WC/∂rCA·DABC = ∂CNC/∂rCA·(WA·WB·∂WC/∂rCA)·DABC
            x5 = x6*dcn(i,k)*( atm*fdmp )
!           WA·WB·∂WC/∂rBC·DABC = ∂CNC/∂rBC·(WA·WB·∂WC/∂rBC)·DABC
            x6 = x6*dcn(j,k)*( atm*fdmp )
!           WA·WB·WC·∂DABC/∂rAB
            x7 = c9ijk*( atm*dijfdmp-dijatm*fdmp )
!           WA·WB·WC·∂DABC/∂rBC
            x8 = c9ijk*( atm*djkfdmp-djkatm*fdmp )
!           WA·WB·WC·∂DABC/∂rCA
            x9 = c9ijk*( atm*dikfdmp-dikatm*fdmp )

!       --- build everything together ---
            eabc = eabc + c9ijk*atm*fdmp

!           ∂rAB/∂A = -∂rAB/∂B
            drij = rij*oorij
!           ∂rBC/∂B = -∂rBC/∂C
            drjk = rjk*oorjk
!           ∂rCA/∂C = -∂rCA/∂A
            drik = rik*oorik

!           ∂EABC/∂A =
!           + (∂WA/∂rAB-∂WA/∂rCA)·WB·WC·DABC
!           + WA·∂WB/∂rAB·WC·DABC
!           - WA·WB·∂WC/∂rCA·DABC
!           + WA·WB·WC·(∂DABC/∂rAB-∂DABC/∂rCA)
            g(:,i) = g(:,i) + ( &
            &        + (x1+x4+x7)*drij &
            &        - (x2+x5+x9)*drik )
!           ∂EABC/∂B =
!           - ∂WA/∂rAB·WB·WC·DABC
!           + WA·(∂WB/∂rBC-∂WB/∂rAB)·WC·DABC
!           + WA·WB·∂WC/∂rBC·DABC
!           + WA·WB·WC·(∂DABC/∂rBC-∂DABC/∂rAB)
            g(:,j) = g(:,j) + ( &
            &        - (x1+x4+x7)*drij &
            &        + (x3+x6+x8)*drjk )
!           ∂EABC/∂C =
!           + ∂WA/∂rCA·WB·WC·DABC
!           - WA·∂WB/∂rBC·WC·DABC
!           + WA·WB·(∂WC/∂rCA-∂WC/∂rBC)·DABC
!           + WA·WB·WC·(∂DABC/∂rCA-∂DABC/∂rBC)
            g(:,k) = g(:,k) + ( &
            &        + (x2+x5+x9)*drik &
            &        - (x3+x6+x8)*drjk )

         enddo ! k/C
      enddo ! j/B
   enddo ! i/A
!$omp enddo
!$omp endparallel

   if (present(eout)) eout=eabc

end subroutine dabcgrad


subroutine d4_numgrad(nat,ndim,at,q,xyz, &
           &          s6,s8,s9,a1,a2,alp,wf,g_a,g_c,mbd,numg)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : covncoord,d4,edisp
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat) 
   real(wp),intent(in)  :: xyz(3,nat) 
   real(wp),intent(in)  :: s6,s8,s9,a1,a2
   integer, intent(in)  :: alp
   real(wp),intent(in)  :: wf,g_a,g_c
   integer, intent(in)  :: mbd
   real(wp),intent(inout) :: numg(3,nat)

   integer  :: k,l
   real(wp) :: er,el,step,cn_thr
   real(wp),allocatable :: xyzl(:,:)
   real(wp),allocatable :: xyzr(:,:)
   real(wp),allocatable :: xyzdup(:,:)
   real(wp),allocatable :: c6abns(:,:)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: gw(:)
   parameter(step=1e-4_wp)
   parameter(cn_thr=1600._wp)

   allocate( xyzr(3,nat),xyzl(3,nat), source=xyz )
   allocate( c6abns(ndim,ndim),cn(nat),gw(ndim), &
   &         source = 0._wp )

   do l = 1, 3
      do k = 1, nat
         er = 0._wp
         el = 0._wp
         cn=0._wp
         gw=0._wp
         c6abns=0._wp
         xyzr = xyz
         xyzr(l,k) = xyz(l,k) + step
         call covncoord(nat,at,xyzr,cn,cn_thr)
         call d4(nat,ndim,at,wf,g_a,g_c,cn,gw,c6abns)
         call edisp(nat,ndim,at,q,xyzr,s6,s8,s9,a1,a2,alp,g_a,g_c, &
         &          gw,c6abns,mbd,er)
         cn=0._wp
         gw=0._wp
         c6abns=0._wp
         xyzl = xyz
         xyzl(l,k) = xyz(l,k) - step
         call covncoord(nat,at,xyzl,cn,cn_thr)
         call d4(nat,ndim,at,wf,g_a,g_c,cn,gw,c6abns)
         call edisp(nat,ndim,at,q,xyzl,s6,s8,s9,a1,a2,alp,g_a,g_c, &
         &          gw,c6abns,mbd,el)
         numg(l,k) = numg(l,k) + (er-el)/(2*step)
      enddo
   enddo

end subroutine d4_numgrad

subroutine dispmb(E,aw,xyz,oor6ab,nat)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : trapzd, ootpi
   implicit none
   integer, intent(in)  :: nat
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: aw(23,nat)
   real(wp),intent(in)  :: oor6ab(nat,nat)
   real(wp),intent(out) :: E

   integer  :: i,j,ii,jj,k
   integer  :: info
   real(wp) :: tau(3,3),spur(23),d_,d2,r(3),r2,alpha
   real(wp) :: two(23),atm(23),d3
   real(wp),allocatable :: T (:,:)
   real(wp),allocatable :: A (:,:)
   real(wp),allocatable :: AT(:,:)
   real(wp),allocatable :: F (:,:)
   real(wp),allocatable :: F_(:,:)
   real(wp),allocatable :: d (:)
   real(wp),allocatable :: w (:)

   intrinsic :: sum,sqrt,minval,log

   allocate( T(3*nat,3*nat),  A(3*nat,3*nat), AT(3*nat,3*nat), &
   &         F(3*nat,3*nat), F_(3*nat,3*nat),  d(3*nat), &
   &         w(12*nat), &
   &         source = 0.0_wp )

   spur = 0.0_wp

   do i = 1, 3*nat
      F(i,i) = 1.0_wp
   enddo

   do i = 1, nat
      do j  = 1, i-1
         r  = xyz(:,j) - xyz(:,i)
         r2 = sum(r**2)
         do ii = 1, 3
            tau(ii,ii) = (3*r(ii)*r(ii)-r2)/r2
            do jj = ii+1, 3
               tau(ii,jj) = (3*r(ii)*r(jj))/r2
               tau(jj,ii) = tau(ii,jj)
            enddo
         enddo
         tau = tau*sqrt(oor6ab(i,j))
         T(3*i-2:3*i,3*j-2:3*j) = tau
         T(3*j-2:3*j,3*i-2:3*i) = tau
      enddo
   enddo

   !call prmat(6,T,3*nat,3*nat,'T')

   do k = 1, 23
      A = 0.0_wp
      do i =  1, nat
         alpha = sqrt(aw(k,i))
         A(3*i-2,3*i-2) = alpha
         A(3*i-1,3*i-1) = alpha
         A(3*i  ,3*i  ) = alpha
      enddo

      AT  = 0.0d0 
      call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,T, &
  &              3*nat,0.0_wp,F_,3*nat)
      call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,F_,3*nat,A, &
  &              3*nat,0.0_wp,AT,3*nat)

      F_ = F - AT

      d = 0.0d0
      call dsyev('N','U',3*nat,F_,3*nat,d,w,12*nat,info)
      if (info.ne.0) then
!        call raise('W','MBD eigenvalue not solvable')
         print'(x''* MBD eigenvalue not solvable'')'
         E = 0.0_wp
         return
      endif
      if (minval(d).le.0.0d0) then
!        call raise('W','Negative MBD eigenvalue occurred')
         print'(x''* Negative MBD eigenvalue occurred'')'
         E = 0.0_wp
         return
      endif

      call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,AT,3*nat,AT, &
  &              3*nat,0.0_wp,F_,3*nat)
!     call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,F_,3*nat,AT, &
! &              3*nat,0.0_wp,A,3*nat)
       
      d_ = 1.0_wp; d2 = 0.0_wp!; d3 = 0.0_wp
      do i = 1, 3*nat
         d_ = d_ * d(i)
         d2 = d2 - F_(i,i)
!        d3 = d3 - A(i,i)
      enddo
      spur(k) = log(d_) - d2*0.5
!     two(k) = d2/2.0_wp
!     atm(k) = d3/3.0_wp
   enddo

   E = trapzd(spur)*ooTPI
   !print*,'     full contribution', trapzd(spur)*ooTPI
   !print*,' manybody contribution', trapzd(spur-two)*ooTPI
   !print*,'  twobody contribution', trapzd(two)*ootpi
   !print*,'threebody contribution', trapzd(atm)*ootpi

   deallocate(T,A,AT,F,F_,d)
end subroutine dispmb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING! !!
!! This implementation of the MBD gradient is incorrect, DO NOT USE! !!
!! !WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING! !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mbdgrad(nat,xyz,aw,daw,oor6ab,g,E)
   use xtb_mctc_accuracy, only: wp
!   use precision, only : wp => dp
!  use dftd4, only : trapzd, ootpi
   implicit none
   integer, intent(in)  :: nat
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: aw(23,nat)
   real(wp),intent(in)  :: daw(23,3,nat)
   real(wp),intent(in)  :: oor6ab(nat,nat)
   real(wp),intent(inout)        :: g(3,nat)
   real(wp),intent(out),optional :: E

   integer  :: i,ii,j,jj,k,kk,l,ll,m,n
   integer  :: info
   real(wp) :: spur(23)
   real(wp) :: fdmp,dfdmp
   real(wp) :: dum,dum2,ddum(3)
   real(wp) :: r(3),r2,drij(3)
   real(wp) :: alpha
   real(wp) :: tau(3,3),dtau(3,3,3)
   real(wp),allocatable :: dspur(:,:,:)
   real(wp),allocatable :: dA(:,:)
   real(wp),allocatable :: dT(:,:)
   real(wp),allocatable :: dAT(:,:)
   real(wp),allocatable :: nT(:,:,:) ! nabla T
   real(wp),allocatable :: T(:,:)
   real(wp),allocatable :: A(:,:)
   real(wp),allocatable :: AT(:,:)
   real(wp),allocatable :: dF(:,:)
   real(wp),allocatable :: F(:,:)
   real(wp),allocatable :: tmp1(:,:)
   real(wp),allocatable :: tmp2(:,:)
   real(wp),allocatable :: d(:)
   real(wp),allocatable :: w(:)

   intrinsic :: sum,sqrt,minval

   ! the MBD calculation in this RPA-like scheme needs quite a lot
   ! of memory, which might be problematic for big systems.
   ! We could use one or maybe to matrices less which would
   ! obfuscate the code
   allocate( T(3*nat,3*nat),A(3*nat,3*nat),AT(3*nat,3*nat), &
   &         F(3*nat,3*nat),dF(3*nat,3*nat), &
   &         tmp1(3*nat,3*nat),tmp2(3*nat,3*nat), &
   &         dT(3*nat,3*nat),dA(3*nat,3*nat), &
   &         dAT(3*nat,3*nat),nT(3*nat,3*nat,3), &
   &         d(3*nat),w(12*nat),dspur(23,3,nat), &
   &         source = 0.0_wp )

   spur = 0.0_wp

   do i = 1, 3*nat
      F(i,i) = 1.0_wp
   enddo

!-----------------------------------------------------------------------
!  interaction tensor setup                                    SAW 1708
!-----------------------------------------------------------------------
   do i = 1, nat
      do j  = 1, i-1
         r  = xyz(:,j) - xyz(:,i)
         r2 = sum(r**2)
         do ii = 1, 3
            tau(ii,ii) = (3*r(ii)**2-r2)/r2
            do jj = ii+1, 3
               tau(ii,jj) = (3*r(ii)*r(jj))/r2
               tau(jj,ii) = tau(ii,jj)
            enddo
         enddo

         fdmp = sqrt(oor6ab(i,j))
         !print*, fdmp
         tau = tau*fdmp
         T(3*i-2:3*i,3*j-2:3*j) = tau
         T(3*j-2:3*j,3*i-2:3*i) = tau
         
!-----------------------------------------------------------------------
!        derivative of interaction tensor                      SAW 1802
!-----------------------------------------------------------------------
         dfdmp = -3*oor6ab(i,j)**2/fdmp*r2**2 ! *sqrt(r2)

         do ii = 1, 3
            dtau(ii,ii,ii) = (15*r(ii)*(0.6_wp*r2-r(ii)**2))/r2**2
            do jj = ii+1, 3
               dtau(jj,ii,ii) = (15*r(jj)*(0.2_wp*r2-r(ii)**2))/r2**2
               dtau(ii,jj,ii) = dtau(jj,ii,ii)
               dtau(ii,ii,jj) = dtau(jj,ii,ii)
               dtau(jj,jj,ii) = (15*r(ii)*(0.2_wp*r2-r(jj)**2))/r2**2
               dtau(jj,ii,jj) = dtau(jj,jj,ii)
               dtau(ii,jj,jj) = dtau(jj,jj,ii)
               do kk = jj+1, 3
                  dtau(ii,jj,kk) = -(15*r(ii)*r(jj)*r(kk))/r2**2
                  dtau(jj,kk,ii) = dtau(ii,jj,kk)
                  dtau(kk,ii,jj) = dtau(ii,jj,kk)
                  dtau(kk,jj,ii) = dtau(ii,jj,kk)
                  dtau(ii,kk,jj) = dtau(ii,jj,kk)
                  dtau(jj,ii,kk) = dtau(ii,jj,kk)
               enddo
            enddo
         enddo

         !drij = r/sqrt(r2)

         dtau(:,:,1) = ( dtau(:,:,1)*fdmp + tau*dfdmp*r(1) )
         dtau(:,:,2) = ( dtau(:,:,2)*fdmp + tau*dfdmp*r(2) )
         dtau(:,:,3) = ( dtau(:,:,3)*fdmp + tau*dfdmp*r(3) )

         !print*,dtau
         nT(3*i-2:3*i,3*j-2:3*j,1:3) = dtau
         nT(3*j-2:3*j,3*i-2:3*i,1:3) = dtau

      enddo
   enddo
        !call prmat(6,T,3*nat,3*nat,'T')

!-----------------------------------------------------------------------
! RPA-like calculation of MBD energy                           SAW 1708
!-----------------------------------------------------------------------
! EMBD = 1/(2π)∫dω Tr{log[1-AT]} = 1/(2π)∫dω log[∏(i)Λii]
! E(2) = 1/(2π)∫dω ½ Tr{(AT)²} ! this two-body energy has to be removed
!-----------------------------------------------------------------------
   do k = 1, 23
      A = 0.0_wp
      do i =  1, nat
         alpha = sqrt(aw(k,i))
         A(3*i-2,3*i-2) = alpha
         A(3*i-1,3*i-1) = alpha
         A(3*i  ,3*i  ) = alpha
      enddo

      AT  = 0.0_wp
      call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,T, &
           &     3*nat,0.0_wp,tmp1,3*nat)
      call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp1,3*nat,A, &
           &     3*nat,0.0_wp,AT,3*nat)
        !call prmat(6,AT,3*nat,3*nat,'AT')

      tmp1 = F - AT
        !call prmat(6,dF,3*nat,3*nat,'1-AT')

      d = 0.0d0
      call dsyev('V','U',3*nat,tmp1,3*nat,d,w,12*nat,info)
        !call prmat(6,tmp1,3*nat,3*nat,'F eigv.')
      if (info.ne.0) then
!        call raise('W','MBD eigenvalue not solvable')
         print'(x''* MBD eigenvalue not solvable'')'
         E = 0.0_wp
         return
      endif
      if (minval(d).le.0.0d0) then
!        call raise('W','Negative MBD eigenvalue occurred')
         print'(x''* Negative MBD eigenvalue occurred'')'
         E = 0.0_wp
         return
      endif

!     two-body contribution to energy
      call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,AT,3*nat,AT, &
           &     3*nat,0.0_wp,tmp2,3*nat)
        !call prmat(6,dF,3*nat,3*nat,'dF')

       
!     get the trace
      dum = 1.0_wp; dum2 = 0.0_wp
      do i = 1, 3*nat
         dum = dum * d(i)
         dum2 = dum2 + tmp2(i,i)
      enddo
      spur(k) = log(dum) + 0.5_wp * dum2
      !print*, log(dum), 0.5_wp*dum2,spur(k)

!-----------------------------------------------------------------------
! MBD gradient calculation                                     SAW 1803
!-----------------------------------------------------------------------
! Some theory:
! EMBD = 1/(2π)∫dω Tr{log[1-AT]} = 1/(2π)∫dω log[∏(i)Λii]
! E(2) = 1/(2π)∫dω ½ Tr{(AT)²}
! ∇EMBD = 1/(2π)∫dω Tr{∇(log[1-AT])} = 1/(2π)∫dω Tr{(1-AT)⁻¹·(∇AT+A∇T)}
!       = 1/(2π)∫dω Tr{(1-A^½TA^½)⁻¹·((∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½))}
! (1-AT)⁻¹ = U·Λ⁻¹·U†  w/ U⁻¹ = U†  if (1-AT) is symmetrized
! Aijkl = δij·δkl·∑(ref)Wi·α
! ∂/∂X Aijkl = δij·δkl·∑(ref)∂Wi/∂CNi·∂CNi/∂rij·∂rij/∂X·α
!-----------------------------------------------------------------------
! eigenvectors are still saved on tmp1, eigenvalues are still on d
! we (over)use tmp2 to hold intermediate results
!-----------------------------------------------------------------------

!     get the inverse of 1-AT by using its eigenvalues
      dF = 0.0_wp
      do i = 1, 3*nat
         dF(i,i) = 1._wp/d(i)
      enddo
        !call prmat(6,dF,3*nat,3*nat,'dF')
!     (1-AT)⁻¹ = U·Λ⁻¹·U†
      call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp1,3*nat,dF, &
           &     3*nat,0.0_wp,tmp2,3*nat)
!     call dgemm('N','T',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,tmp1, &
!          &     3*nat,0.0_wp,AT,3*nat)
        !call prmat(6,dF,3*nat,3*nat,'dF')
!     unfortunately this might not be enough, we have to substract the
!     twobody gradient for the dipole-dipole interaction, which is in
!     fact not easily accessable from here.
!     If this is correct,
!        E(2) = 1/(2π)∫dω ½ Tr{(AT)²}
!       ∇E(2) = 1/(2π)∫dω Tr{(AT)·(∇AT+A∇T)}
!     then the MBD gradient w/o two-body contrib. could be represented by
!     dF = dF - AT
!     or more easier by the use of dgemm's beta by replacing the last
!     dgemm by:
!     call dgemm('N','T',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,tmp1, &
!          &     3*nat,1.0_wp,AT,3*nat)
!     But we would have to consider AT instead of dF in the following code
!     which would be at least confusing

!-----------------------------------------------------------------------
!     Efficient MBD-Gradient calculation, O(N³)                SAW 1803
!-----------------------------------------------------------------------
      do l = 1, 3
         dA = 0.0_wp
         do i = 1, nat
            dA(3*i-2,3*i-2) = 0.5*daw(k,l,i)/sqrt(aw(k,i))
            dA(3*i-1,3*i-1) = 0.5*daw(k,l,i)/sqrt(aw(k,i))
            dA(3*i  ,3*i  ) = 0.5*daw(k,l,i)/sqrt(aw(k,i))
         enddo

!        (∇A^½)TA^½
         call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,dA,3*nat,T, &
              &     3*nat,0.0_wp,tmp2,3*nat)
         call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,A, &
              &     3*nat,0.0_wp,dAT,3*nat)

!        (∇A^½)TA^½+A^½(∇T)A^½ (please note the use of beta=1.0!)
         call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,nT(:,:,l), &
              &     3*nat,0.0_wp,tmp2,3*nat)
         call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,A, &
              &     3*nat,1.0_wp,dAT,3*nat)

!        last term and we have: (∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½)
         call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,T, &
              &     3*nat,0.0_wp,tmp2,3*nat)
         call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,dA, &
              &     3*nat,1.0_wp,dAT,3*nat)

!        by multiplying the prefactor we get to:
!        ((1-A^½TA^½)⁻¹+A^½TA^½)·((∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½))
         call dgemm('N','N',3*nat,3*nat,3*nat,1.0_wp,AT,3*nat,dAT, &
              &     3*nat,0.0_wp,tmp2,3*nat)
!        now the other term of
!        [(1-A^½TA^½)⁻¹+A^½TA^½),((∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½)]
         call dgemm('N','N',3*nat,3*nat,3*nat,-1.0_wp,dAT,3*nat,AT, &
              &     3*nat,1.0_wp,tmp2,3*nat)

         do i = 1, nat
            dspur(k,l,i) = dspur(k,l,i) &
            &              + tmp2(3*i-2,3*i-2) &
            &              + tmp2(3*i-1,3*i-1) &
            &              + tmp2(3*i  ,3*i  )
         enddo ! all atoms
         
      enddo ! cartesian parts

   enddo ! k, imaginary frequencies

   E = trapzd(spur)*ootpi
   !print*, E
   do i = 1, nat
      do j = 1, 3
         g(j,i) = g(j,i) + ootpi*trapzd(dspur(:,j,i))
         !print*,ootpi*trapzd(dspur(:,j,i))
      enddo
   enddo

   deallocate( T,A,AT,F,tmp1,tmp2,dT,dA,dAT,dF,nT,d,w,dspur )

end subroutine mbdgrad

end module dftd4
