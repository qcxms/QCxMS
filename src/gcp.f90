subroutine calcgcp(n,iz,xyz,echo,method,ebsse,g)
!use io_reader
implicit none
integer n,max_elem,max_para
integer nn,nb,i,zz,j
!parameter (maxat =500)
! max elements possible
parameter (max_elem=94)
! actually parameterized
parameter (max_para=36)
real*8 autoang
parameter (autoang =0.52917726d0)
! coordinates in au
real*8 xyz(3,n)
! gradient
real*8 g  (3,n)
! cardinal numbers of elements
integer   iz(n)
integer   iz2(n)
! frozen coordinate
!integer ifrez(maxat)
! #electrons of elements
integer   nel(max_elem)
! # basis functions                   
integer nbas(max_elem)
real*8 emiss(max_elem)
! parameters
real*8 p(4) 
real*8 xva(n),xvb(n)
real*8  xx(10),ebsse,r
real*8 t0,t1
logical ex,echo,grad,DoHess,helpme,parfile,localF,warn,lib,parm,test,verb
character*20 method,ctmp
character*2 esym
character*80 arg(10)
character*80 atmp,ftmp,dtmp
character*80 basname        

!** subroutine ?? **
  lib=.true. 
!*******************

!************
!* Defaults *
!************
warn=.false.
parfile=.false.
parm=.false.
DoHess=.false.
helpme=.false.
grad=.true.
localF=.false.
verb=.false.
test=.false.
basname=''


if(echo)then 
write(*,*)'---------------------------------------------------'        
write(*,*)'                   g C P'        
write(*,*)'---------------------------------------------------'        
endif

if(n.lt.1)     stop 'no atoms' 
!if(n.gt.maxat) stop 'too many atoms' 

if (helpme) call help(.true.)

!************************
!* modify input string  *
!* 1. use lower case    *
!* 2. delete hyphends   *
!************************      
call lower_case(method)
do while (index(method,'-').ne.0) 
  i=scan(method,'-')
  ctmp=trim(method(:(i-1)))//trim(method((i+1):))
  method=trim(ctmp)
enddo

!****************************************************************************
!* read parameter file                                                      *
!* format:                                                                  *
!* basis sigma eta alpha beta                                               *
!* To read in nbas and emiss, generate param file with                      *
!* with '-parfile', edit copy it to .gcppar:                              *
!* A '#' and the beginning will tell the program to read a different format *
!****************************************************************************
if(method.eq.''.or.method.eq.'file')then
  call system('hostname > .tmpx')
  open(unit=43,file='.tmpx')
  read(43,'(a)')ftmp
  close(43,status='delete')  
  write(dtmp,'(''~/.gcppar.'',a)')trim(ftmp)
  if(localF) dtmp='.gcppar'
  inquire(file=dtmp,exist=ex)
  if(ex)then
     if(echo)write(*,*) 'reading ',trim(dtmp)
     open(unit=42,file=dtmp)
         read(42,'(a)',end=9)ftmp
     call charXsplit(ftmp,basname,1)
     call lower_case(basname)
         call readl(ftmp,xx,nn)
         p(1:4)=xx(1:4)
      if(adjustl(trim(basname)).eq.'#') then
       print*,' found extended param file ! '
       parm=.true.    
       read(42,'(x,F8.5)',end=9) p(1:4)
       do i=1,36
          read(42,'(x,I3,x,F8.5)',end=9) nbas(i), emiss(i)
       enddo
      endif
9   close(42)
     if(echo.and..not.parm) then
      write(*,*) 'loading ',trim(basname),' params'
      call setparam(emiss,nbas,p,basname)
     endif
  else
   print*,''
   print*, 'found not param file in ',trim(dtmp)
   stop '* ERROR *'
  endif
method=trim(basname)
!*******************
!* built-in method *
!*******************
else
    call setparam(emiss,nbas,p,method)
endif


!**********************************************************
!* check for missing parameters for current elements      *
!* elements >36 will be treated as their lower homologues *
!* eg: Ag -> Cu, In -> Ga                                 *
!* prepare virtual bf array xva,xvb                       *
!**********************************************************
! backup original
iz2=iz
! re-arrange elements
do i=1,n
zz=iz(i)
select case (zz)
 case(37:54)
  iz(i)=iz(i)-18
  print*, '** NOTE ** -> element ',esym(zz),' will be treated as ', esym(iz(i))
  warn=.true.
 case(55:57)
  iz(i)=iz(i)-18*2
  print*, '** NOTE ** -> element ',esym(zz),' will be treated as ', esym(iz(i))
  warn=.true.
 case(58:71,90:94)
  iz(i)=21
  print*, '** NOTE ** -> element ',esym(zz),' will be treated as ', esym(iz(i))
  warn=.true.
 case(72:89)
  iz(i)=iz(i)-2*18-14
  print*, '** NOTE ** -> element ',esym(zz),' will be treated as ', esym(iz(i))
  warn=.true.
end select 
enddo

! set electrons
do i=1,maxval(iz)
  nel(i)=i
enddo


!* handle ECP basis set
if(index(method,'vmb').ne.0) then
print*, 'VMB basis has ECPs. Adjusting electrons'
do i=5,10
nel(i)=nel(i)-2
enddo
do i=11,18
nel(i)=nel(i)-10
enddo
endif

! count bf
nb=0
do nn=1,n
   nb=nb+nbas(iz(nn))
enddo
if(nb.lt.1) stop 'Nbf setup gone wrong'

! set virtuals and look for missing parameters
xva=0
do i=1,n
  xva(i)=(nbas(iz(i))-0.5d0*nel(iz(i)))  
  if(emiss(iz(i)).gt.01e-6.and.xva(i).lt.0.0d0) then
  print*, 'element/emiss/nvirt/nel/nbas'
  print*,esym(iz(i)),emiss(iz(i)),xva(i),nel(iz(i)),nbas(iz(i))
   stop 'negative number of virtual orbitals. Something is wrong in the parameters !'
  endif
  if( emiss(iz(i)).lt.0.1e-6) then
    print*, '** WARNING ** -> element ',esym(iz(i)),' has no parameters       (no contribution)!'
    xva(i)=0
  cycle
  warn=.true.
  endif
  if(xva(i).lt.0.5d0) then
    print*, '** WARNING ** -> element ',esym(iz(i)),' has no virtual orbitals (no contribution)!'
  warn=.true.
  endif
enddo
xvb=xva


!*********************************
!* Parameter handling/setup done *
!*********************************
  if(echo)then
     write(*,*)''
     write(*,'(2x,''level '',3x,a12)') adjustr(trim(method))
     write(*,'(2x,''Nbf   '',3x,I12)')  nb
     write(*,'(2x,''Atoms '',3x,I12)')  n
     write(*,*)''
     write(*,'(2x,a)') 'Parameters: '
     write(*,'(2x,a,2x,f10.4)') 'sigma ',p(1)
     write(*,'(2x,a,2x,f10.4)') 'eta   ',p(2)
     write(*,'(2x,a,2x,f10.4)') 'alpha ',p(3)
     write(*,'(2x,a,2x,f10.4)') 'beta  ',p(4)
     write(*,*)''
   endif



! **************************
! * print parameter table  *
! **************************
if(echo) then
  write(*,*) ' '
  write(*,*) 'element parameters ',trim(method)
  write(*,'(2x,3(a4,2x,a6,3x,a4,4x))') 'elem','emiss','nbas','elem','emiss','nbas','elem','emiss','nbas'
  do i=1,36,3
      write(*,'(2x,3(A4,2x,F8.5,2x,I3,4x))') esym(i), emiss(i),nbas(i), &
                                           esym(i+1), emiss(i+1),nbas(i+1),&
                                           esym(i+2), emiss(i+2),nbas(i+2)
  enddo
  write(*,*) ' '
endif


! **************************
! * write parameter file   *
! **************************
if(parfile) then
  print*, 'printing extended parameter file ...'
  print*, 'modfiy & copy gcp.param to $HOME/.gcppar.$HOSTNAME'
  open(unit=22,file='gcp.param')
    write(22,*) ' # ',trim(method),' (comment line)'
    write(22,'(x,F8.5)') p(1:4)
    do i=1,36
      write(22,'(x,I3,x,F8.5)') nbas(i), emiss(i)
    enddo
  close(22)
endif


if(test) then
do i=1,n
 do j=1,n
    if(i.eq.j) cycle
r=sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2)*0.52917726d0
 if(r.lt.0.7) then
 write(*,'('' Distance between '',I5,''('',a2,'')'',a,I5,     &
 ''('',a2,'') only about '',F8.2,'' angstrom !'')') i,esym(iz(i)),' <-> ',j,esym(iz(j)),r
endif
enddo
enddo

call done('test run finished',6)
endif


!**************************
!* calc energy & gradient *
!**************************
      g(1:3,1:n)=0.0d0
  call e(n,max_elem,emiss,xyz,iz,p,ebsse,g,grad,echo,xva,xvb)

! restore original cardinal numbers
iz=iz2
! set cart. frozen coords to zero
!if(maxval(ifrez,n).eq.1.and.grad) then 
!if(echo) print*,' setting gradient of frozen cartesians to zero'
!  do i=1,n
!   if(ifrez(i).eq.1) g(1:3,i)=0.0d0 
!  enddo
!endif

!    if(grad.and..not.verb) call wregrad(maxat,n,xyz,iz,ebsse,g,echo)

  if(verb) then
   write(*,'(''gradient: Ggcp'')')
   do i=1,n
    write(*,'(3E22.13)')g(1,i),g(2,i),g(3,i)
   enddo
  endif

    if(echo)then
!      write(*,*) '                '
!      write(*,*) '** gCP correction ** '
!      write(*,'(2x,a7,F18.10,'' / (a.u.) || '',x,F11.4,'' / (kcal/mol)'')') 'Egcp:  ', gcp,gcp*2625.497656d0
!      write(*,*) '                '
       if(grad)then
       write(*,*)'|G|=',sum(abs(g(1:3,1:n)))
       endif
       write(*,*)'---------------------------------------------------'        
    endif
!   open(unit=1,file='.CPC')
!   write(1,*) ebsse
!   close(1)

! calc hessian
!     call hessian(n,max_elem,emiss,xyz,iz,p,ebsse,xva,xvb,DoHess)
 
!   call cpu_time(t1)
!   if(echo) call prtim(6,t1-t0,'t','gCP ')
    if(warn) write(*,*) 'Carefully read the notifications/warnings given at loadup'
!   if(echo) call done('gCP',6)
end 


!***************************************************
!* This subroutine calculates the gcp correction *
!***************************************************
subroutine e(n,max_elem,emiss,xyz,iz,p,ebsse,g,grad,echo,xva,xvb)
implicit none             
integer iz(n),n,max_elem,np
integer iat,jat
real*8  xyz(3,n)
real*8  g  (3,n)
real*8  ebsse,tg,t0,t1
real*8  emiss(max_elem),dum22
real*8  p(4),thrR,thrE
real*8 xva(*),xvb(*)
real*8 va,vb
real*8 r,dx,dy,dz
real*8 damp,tmp,ea,ecp,dum,tmpa,tmpb,dum2
real*8 sab,p1,p2,p3,p4
real*8 gs(3),gab
real*8 za(36),zb(36)
logical echo,grad
  
! Two threshold. thrR: distance cutoff thrE: numerical noise cutoff
thrR=60            ! 60 bohr
thrE=epsilon(1.d0) ! cut below machine precision rounding

if(echo) then
write(*,*) '  '     
write(*,'('' cutoffs: '',F5.1)') 
write(*,'(''   distance [bohr]'',F5.1)') thrR
write(*,'(''   noise    [a.u.]'',Es8.1)') thrE
write(*,*) '  '     
endif

g=0
ecp=0.0d0
dum=0
tg=0d0

p1=abs(p(1))
p2=abs(p(2))
p3=abs(p(3))
p4=abs(p(4))

call setzet(p2,za,zb)
gs=0
if(echo) write(*,'(2x,a5,2x,a5,2x,a5,2x,a7,2x,a14,4x,a15)') &
  '#','ON','sites','Nvirt','Emiss','BSSE [kcal/mol]'

! Loop over all i atoms
      do iat=1,n
       va=xva(iat)
         ea=0.0d0

         np=0
! the BSSE due to atom jat, Loop over all j atoms
         do jat=1,n
            if(iat.eq.jat) cycle
            dx=(xyz(1,iat)-xyz(1,jat))
            dy=(xyz(2,iat)-xyz(2,jat))
            dz=(xyz(3,iat)-xyz(3,jat))
            r=sqrt(dx*dx+dy*dy+dz*dz)

! # of bf that are available from jat
            vb=xvb(jat)
! check virtuals
            if(vb.lt.0.5) cycle
! distance cutoff
            if(r.gt.thrR) cycle
! calulate slater overlap sab
          call ssovl(r,iat,jat,iz,za(iz(iat)),zb(iz(jat)),sab)
! noise cutoff(sab)
           if(abs(sab).lt.thrE) cycle
! evaluate gcp central expression
           tmpa=exp(-p3*r**p4)
           tmpb=sqrt(vb*Sab)
           damp=tmpa/tmpb
! noise cutoff(damp)
           if(abs(damp).lt.thrE) cycle
           ea=ea+emiss(iz(iat))*damp

! sites counter (i.e. # atoms contributing to the 'atomic bsse')
           np=np+1

! gradient for i,j pair
         if(grad)then
           call cpu_time(t0)
          call gsovl(r,iat,jat,iz,za(iz(iat)),zb(iz(jat)),gab)

          gs(1)=gab*dx
          gs(2)=gab*dy
          gs(3)=gab*dz


          dum=exp(-p3*r**p4)*(-1d0/2d0)
          dum2=2d0*p3*p4*r**p4*sab/r
          dum22=r*sab**(3d0/2d0)
          tmpb=dum22*sqrt(vb)

          tmpa=dum2*dx+gs(1)
          tmp=dum*tmpa/tmpb
          g(1,iat)=g(1,iat)+tmp*emiss(iz(iat))

          tmpa=dum2*dy+gs(2)
          tmp=dum*tmpa/tmpb
          g(2,iat)=g(2,iat)+tmp*emiss(iz(iat))

          tmpa=dum2*dz+gs(3)
          tmp=dum*tmpa/tmpb
          g(3,iat)=g(3,iat)+tmp*emiss(iz(iat))

         if(va.lt.0.5) cycle

          tmpb=dum22*sqrt(va)
          
          tmpa=dum2*(-dx)-gs(1)
          tmp=dum*tmpa/tmpb

           g(1,iat)=g(1,iat)-tmp*emiss(iz(jat))

          tmpa=dum2*(-dy)-gs(2)
          tmp=dum*tmpa/tmpb
           g(2,iat)=g(2,iat)-tmp*emiss(iz(jat))

          tmpa=dum2*(-dz)-gs(3)
          tmp=dum*tmpa/tmpb
           g(3,iat)=g(3,iat)-tmp*emiss(iz(jat))

           call cpu_time(t1)
           tg=tg+t1-t0
         endif


! end of j-loop
       enddo

         if(echo) &
       write(*,'(2x,3(I5,2x),2x,F5.1,2x,F14.4,2x,F14.4,2x,a)')  &
                iat,iz(iat),np,va,emiss(iz(iat)), ea*627.5099*p(1)

         ecp=ecp+ea             
! end of i-loop
      enddo
!     if(grad.and.echo)  call prtim(6,tg,'t','gradient')

      ebsse=ecp*p(1)
      g=g*p(1)
      end

!*****************
!* print timings *
!*****************
subroutine prtim(io,tt,is,string)
integer io
real*8 tt,t,sec
integer day,hour,min
character*(*) string
character*1 is

t=tt
day=idint(t/86400)
t=t-day*86400
hour=idint(t/3600)
t=t-hour*3600
min=idint(t/60)
t=t-60*min
sec=t

if(day.ne.0)then
   if(is=='w')write(io,2)trim(string),day,hour,min,sec
   if(is=='t')write(io,22)trim(string),day,hour,min,sec
   return
endif
if(hour.ne.0)then
   if(is=='w')write(io,3)trim(string),hour,min,sec
   if(is=='t')write(io,33)trim(string),hour,min,sec
   return
endif
if(min .ne.0)then
   if(is=='w')write(io,4)trim(string),min,sec
   if(is=='t')write(io,44)trim(string),min,sec
   return
endif
  if(is=='w')write(io,5)trim(string),sec
  if(is=='t')write(io,55)trim(string),sec
return

 1    format('======================================')
 2    format('wall-time for ',a,2x,i3,' d',i3,' h',i3,' m',f5.1,' s')
 3    format('wall-time for ',a,2x,i3,' h',i3,' m',f5.1,' s')
 4    format('wall-time for ',a,2x,i3,' m',f5.1,' s')
 5    format('wall-time for ',a,2x,f5.1,' s')

 22    format('cpu-time for ',a,2x,i3,' d',i3,' h',i3,' m',f5.1,' s')
 33    format('cpu-time for ',a,2x,i3,' h',i3,' m',f5.1,' s')
 44    format('cpu-time for ',a,2x,i3,' m',f5.1,' s')
 55    format('cpu-time for ',a,2x,f5.1,' s')

return
end

!*********************************
!* Load all necessary parameters *
!*********************************
subroutine setparam(emiss,nbas,p,method)
implicit none
integer, parameter :: mpar=86 ! maximal
integer, parameter :: apar=36 ! actual
integer nbas(mpar)
character*(*) method
real(8) emiss(mpar),p(*)
real(8) HFsv(apar),HFminis(apar),HF631gd(apar),HFsvp(apar),HFtz(apar),HFvmb(apar),HFminisd(apar)
integer BASsv(apar),BASminis(apar),BAS631gd(apar),BAStz(apar),BASsvp(apar),BASvmb(apar),BASminisd(apar)

! ***************
! *  Emiss data *
! ***************
data HFsv/ & !H-Kr (no 3d)
0.009037,0.008843,&  ! He,He
0.204189,0.107747,0.049530,0.055482,0.072823,0.100847,0.134029,0.174222,&  ! Li-Ne
0.315616,0.261123,0.168568,0.152287,0.146909,0.168248,0.187882,0.211160,&  !Na -Ar
0.374252,0.460972,&  ! K-Ca
0.444886,0.404993,0.378406,0.373439,0.361245,0.360014,0.362928,0.243801,0.405299,0.396510,&   ! 3d-TM
0.362671,0.360457,0.363355,0.384170,0.399698,0.417307/ !Ga-Kr

data HFminis/ &
0.042400,0.028324,&
0.252661,0.197201,0.224237,0.279950,0.357906,0.479012,0.638518,0.832349, &
1.232920,1.343390,1.448280,1.613360,1.768140,1.992010,2.233110,2.493230, &
3.029640,3.389980,&  ! H-Ca
10*0,6*0/

data HF631GD/ &! H-Ca + Br (no 3d)
0.010083,0.008147,&
0.069260,0.030540,0.032736,0.021407,0.024248,0.036746,0.052733,0.075120,& 
0.104255,0.070691,0.100260,0.072534,0.054099,0.056408,0.056025,0.057578,&
0.079198,0.161462,&
10*0.0, &
0.000000,0.000000,0.000000,0.000000,0.381049,0.000000/ 

data HFsvp /  & ! H-Kr
0.008107,0.008045,&
0.113583,0.028371,0.049369,0.055376,0.072785,0.100310,0.133273,0.173600,&
0.181140,0.125558,0.167188,0.149843,0.145396,0.164308,0.182990,0.205668,&
0.200956,0.299661, &
0.325995,0.305488,0.291723,0.293801,0.29179,0.296729,0.304603,0.242041,0.354186,0.350715,&
0.350021,0.345779,0.349532,0.367305,0.382008,0.399709/

data HFtz /&  ! H-Kr
0.007577,0.003312,&
0.086763,0.009962,0.013964,0.005997,0.004731,0.005687,0.006367,0.007511,&
0.077721,0.050003,0.068317,0.041830,0.025796,0.025512,0.023345,0.022734,&
0.097241,0.099167,&
0.219194,0.189098,0.164378,0.147238,0.137298,0.12751,0.118589,0.0318653,0.120985,0.0568313, &
0.090996,0.071820,0.063562,0.064241,0.061848,0.061021/

data HFvmb/&
0.042400,0.028324,&
0.252661,0.197201,0.156009,0.164586,0.169273,0.214704,0.729138,0.336072,&
0.262329,0.231722,0.251169,0.287795,0.213590,0.250524,0.728579,0.260658, &
2*0,&
16*0/

data HFminisd/& !Al-Ar MINIS + Ahlrichs "P" funktions
0.0,0.0,&
8*0.0,&
2*0.0,1.446950,1.610980,1.766610,1.988230,2.228450,2.487960,&
2*0.0,&
16*0.0/

! *********************
! * nr. of basis fkt  *
! *********************
data BASsv/2*2,2*3,6*9,2*7,6*13,2*11,10*21,6*27/           
data BASminis/2*1,2*2,6*5,2*6,6*9,2*10,16*0/
data BAS631gd/2,5,8*14,8*18,2*22,10*0,4*0,32,0/
data BASsvp/2*5,9,9,6*14,15,18,6*18,24,24,10*31,6*32/  
data BAStz/2*6,14,19,6*31,2*32,6*37,33,36,9*45,48,6*48/  
data BASvmb/2*1,2*2,6*4,2*1,6*4,2*0,16*0/ ! minimal basis set with ECPs 
data BASminisd/2*0,2*0,6*0,2*0,6*14,2*0,16*0/

! **************************************
! * load data into emiss() and nbas()  *
! **************************************

emiss=0d0
nbas =0d0

select case (method)
!*****************
!* Hartree-Fock  *
!*****************
  case ('hf/sv') ! RMS=0.3218975
     emiss(1:apar)=HFsv(1:apar)
     nbas(1:apar)=BASsv(1:apar)
     p(1)=0.1724d0
     p(2)=1.2804d0
     p(3)=0.8568d0
     p(4)=1.2342d0
  case ('hf/svp') ! RMS=0.4065
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     p(1)=0.2054d0
     p(2)=1.3157d0
     p(3)=0.8136d0
     p(4)=1.2572d0
  case ('hf/631gd') ! RMS= 0.40476
     emiss(1:apar)=HF631gd(1:apar)
     nbas(1:apar)=BAS631gd(1:apar)
     p(1)=0.2048d0
     p(2)=1.5652d0
     p(3)=0.9447d0
     p(4)=1.2100d0
  case ('hf/minis2') ! RMS= 0.311279   ! MIDI for O,N
    emiss(1:apar)=HFminis(1:apar)
    nbas(1:apar)=BASminis(1:apar)
    emiss(7)=0.323457 
    emiss(8)=0.449484
    nbas(7)=9
    nbas(8)=9
     p(1)=0.0638d0 
     p(2)=1.6284d0 
     p(3)=1.1567d0 
     p(4)=1.2156d0 
  case ('hf/minis1') ! RMS= 0.287  ! MIDI for O
    emiss(1:apar)=HFminis(1:apar)
    nbas(1:apar)=BASminis(1:apar)
    emiss(8)=0.449484
    nbas(8)=9
     p(1)=0.0794d0 
     p(2)=1.3631d0 
     p(3)=1.1173d0 
     p(4)=1.2051d0 
  case ('hf/minis') ! RMS= 0.3040
    emiss(1:apar)=HFminis(1:apar)
    nbas(1:apar)=BASminis(1:apar)
     p(1)=0.1290d0
     p(2)=1.1526d0            
     p(3)=1.1549d0                 
     p(4)=1.1763d0                 
  case ('hf/minix') 
!H-Mg: MINIS
    emiss(1:12)=HFminis(1:12)
    nbas(1:12)=BASminis(1:12)
!Al-Ar: MINIS+d (stan. TM exp.)
    emiss(13:18)=HFminisd(13:18)
    nbas(13:18)=BASminisd(13:18)
!K-Zn: SV
    emiss(19:30)=HFsv(19:30)
    nbas(19:30)=BASsv(19:30)
!Ga-Kr: SVP 
    emiss(31:36)=HFsvp(31:36)
    nbas(31:36)=BASsvp(31:36)
! Li,Be,Na,Mg MINIS + p-fkt
    emiss(3)=0.177871
    emiss(4)=0.171596
    nbas(3:4)=5
    emiss(11)=1.114110
    emiss(12)=1.271150
    nbas(11:12)=9
! fit param of HF/MINIS
     p(1)=0.1290d0
     p(2)=1.1526d0            
     p(3)=1.1549d0                 
     p(4)=1.1763d0                 
  case ('hf/vmb') ! RMS= 0.305
! ECPs 
    emiss(1:apar)=HFvmb(1:apar)
    nbas(1:apar)=BASvmb(1:apar)
     p(1)=0.1829d0
     p(2)=1.3366d0
     p(3)=0.9911d0
     p(4)=1.2893d0
  case ('hf/vmbd') ! RMS= 0.3017
! plus 1d(0.7)
! ECPs 
    emiss(1:apar)=HFvmb(1:apar)
    nbas(1:apar)=BASvmb(1:apar)
   emiss(8)=0.216427
   emiss(9)=0.729749
   nbas(8)=9
   nbas(9)=9
     p(1)=0.2083d0
     p(2)=1.6204d0
     p(3)=0.9898d0
     p(4)=1.3493d0
  case ('hf/tz') !  RMS= 0.1150
     emiss(1:apar)=HFtz(1:apar)
     nbas(1:apar)=BAStz(1:apar)
     p(1)=0.3127d0
     p(2)=1.9914d0
     p(3)=1.0216d0
     p(4)=1.2833d0
!************************
!* KS-DFT hybrid(B3LYP) *
!************************
  case ('dft/sv','b3lyp/sv') ! RMS= 0.557
     emiss(1:apar)=HFsv(1:apar)
     nbas(1:apar)=BASsv(1:apar)
     p(1)=0.4048d0
     p(2)=1.1626d0
     p(3)=0.8652d0
     p(4)=1.2375d0
  case ('dft/sv(p)','b3lyp/sv(p)') ! RMS= 0.57 ! def2-SV(P)  
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     p(1)=0.2424d0
     p(2)=1.2371d0
     p(3)=0.6076d0
     p(4)=1.4078d0
  case ('dft/svx','b3lyp/svx') ! RMS=  0.56 ! def2-SV(P/h,c)  = SV at h,c
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     emiss(6)=HFsv(6)
     nbas(6)=BASsv(6)
     p(1)=0.1861d0
     p(2)=1.3200d0
     p(3)=0.6171d0
     p(4)=1.4019d0
  case ('dft/svp','b3lyp/svp') ! RMS=0.6498
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     p(1)=0.2990d0
     p(2)=1.2605d0
     p(3)=0.6438d0
     p(4)=1.3694d0
  case ('dft/631gd','b3lyp/631gd') ! RMS=  0.47856
    emiss(1:apar)=HF631gd(1:apar)
    nbas(1:apar)=BAS631gd(1:apar)
     p(1)=0.3405d0
     p(2)=1.6127d0
     p(3)=0.8589d0
     p(4)=1.2830d0
  case ('dft/minis','b3lyp/minis') ! RMS= 0.3400
    emiss(1:apar)=HFminis(1:apar)
    nbas(1:apar)=BASminis(1:apar)
     p(1)=0.2059d0
     p(2)=0.9722d0
     p(3)=1.1961d0
     p(4)=1.1456d0
  case ('dft/tz','b3lyp/tz') ! RMS=0.19648
     emiss(1:apar)=HFtz(1:apar)
     nbas(1:apar)=BAStz(1:apar) 
     p(1)=0.2905d0
     p(2)=2.2495d0
     p(3)=0.8120d0
     p(4)=1.4412d0
!*****************
!* special cases *
!*****************
  case ('blyp/minis','gga/minis') ! RMS= 0.3462
  emiss(1:apar)=HFminis(1:apar)
  nbas(1:apar)=BASminis(1:apar)
  p(1)=0.1566d0
  p(2)=1.0271d0
  p(3)=1.0732d0
  p(4)=1.1968d0
  case ('tpss/minis') ! RMS= 
  emiss(1:apar)=HFminis(1:apar)
  nbas(1:apar)=BASminis(1:apar)
  p(1)=0.22982d0
  p(2)=1.35401d0
  p(3)=1.47633d0
  p(4)=1.11300d0
  case ('tpss/svp') ! RMS=  0.618
  emiss(1:apar)=HFsvp(1:apar)
  nbas(1:apar)=BASsvp(1:apar)
  p(1)=0.6647d0
  p(2)=1.3306d0
  p(3)=1.0792d0
  p(4)=1.1651d0
  case ('gga/svp','blyp/svp') ! RMS=
    emiss(1:apar)=HFsvp(1:apar)
  nbas(1:apar)=BASsvp(1:apar)
  p(1)=0.6823d0
  p(2)=1.2491d0
  p(3)=0.8225d0
  p(4)=1.2811d0
  case ('blyp/sv','gga/sv') ! RMS = 0.6652
  emiss(1:apar)=HFsv(1:apar)
  nbas(1:apar)=BASsv(1:apar)
  p(1)=0.2727d0
  p(2)=1.4022d0
  p(3)=0.8055d0
  p(4)=1.3000d0
  case ('blyp/tz','gga/tz') !RMS = 0.21408
  emiss(1:apar)=HFtz(1:apar)
  nbas(1:apar)=BAStz(1:apar)
  p(1)=0.1182d0
  p(2)=1.0631d0
  p(3)=1.0510d0
  p(4)=1.1287d0
  case ('pw6b95/minis')  ! RMS = 0.3279929
  emiss(1:apar)=HFminis(1:apar)
  nbas(1:apar)=BASminis(1:apar)
  p(1)=0.21054d0
  p(2)=1.25458d0
  p(3)=1.35003d0
  p(4)=1.14061d0
  case ('pw6b95/svp')  ! RMS = 0.58312
  emiss(1:apar)=HFsvp(1:apar)
  nbas(1:apar)=BASsvp(1:apar)
  p(1)=0.3098d0
  p(2)=1.2373d0
  p(3)=0.6896d0
  p(4)=1.3347d0
!***********************
!* load just emiss/nbf *
!**********************
case('sv')
     emiss(1:apar)=HFsv(1:apar)
     nbas(1:apar)=BASsv(1:apar) 
case('sv(p)')
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     emiss(1)=0.009037d0
     nbas(1)=2
case ('svx') ! RMS=  ! def2-SV(P/h,c)  = SV at h,c
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar)
     emiss(1)=HFsv(1)
     nbas(1)=BASsv(1)
     emiss(6)=HFsv(6)
     nbas(6)=BASsv(6)
case('svp')
     emiss(1:apar)=HFsvp(1:apar)
     nbas(1:apar)=BASsvp(1:apar) 
case('minis')
     emiss(1:apar)=HFminis(1:apar)
     nbas(1:apar)=BASminis(1:apar) 
case('631gd')
     emiss(1:apar)=HF631gd(1:apar)
     nbas(1:apar)=BAS631gd(1:apar) 
case('tz')
     emiss(1:apar)=HFtz(1:apar)
     nbas(1:apar)=BAStz(1:apar) 
case('vmb')
     emiss(1:apar)=HFvmb(1:apar)
     nbas(1:apar)=BASvmb(1:apar) 
case('vmbd')
     emiss(1:apar)=HFvmb(1:apar)
     nbas(1:apar)=BASvmb(1:apar) 
   emiss(8)=0.216427
   emiss(9)=0.729749
   nbas(8)=9
   nbas(9)=9
case ('minix') 
!H-Mg: MINIS
    emiss(1:12)=HFminis(1:12)
    nbas(1:12)=BASminis(1:12)
!Al-Ar: MINIS+d (stan. TM exp.)
    emiss(13:18)=HFminisd(13:18)
    nbas(13:18)=BASminisd(13:18)
!K-Zn: SV
    emiss(19:30)=HFsv(19:30)
    nbas(19:30)=BASsv(19:30)
!Ga-Kr: SVP 
    emiss(31:36)=HFsvp(31:36)
    nbas(31:36)=BASsvp(31:36)
! Li,Be,Na,Mg MINIS + p-fkt
    emiss(3)=0.177871
    emiss(4)=0.171596
    nbas(3:4)=5
    emiss(11)=1.114110
    emiss(12)=1.271150
    nbas(11:12)=9
case default
  print*, '** ',trim(method),' **'
  stop 'not implemented'
end select

end



!******************************************************************************
!* calculates the s-type overlap integral over 1s, 2s and 3s slater functions 
!* added support for 3s functions
!* ovl = overlap integral
!* za  = slater exponent atom A
!* zb  = slater exponent atom B
!* R   = distance between atom A and B
!* Inspired by mopac7.0
!******************************************************************************
subroutine ssovl(r,iat,jat,iz,xza,xzb,ovl)
implicit none
integer ii,shell(72)
logical debug
real(8) za,zb,R,ovl,ax,bx,norm,R05
real(8) A0,A1,A2,A4,A3,A5,A6
real(8) B0,B1,B2,B3,B4,B5,B6,bint
integer na,nb
real(8) Bxx0,Bxx1,Bxx2,xx,Bxx4,Bxx6
real(8) Bxx3,Bxx5
data shell/                 &
!          h, he
          1,1               &
!         li-ne
          ,2,2,2,2,2,2,2,2, &
!         na-ar
          3,3,3,3,3,3,3,3,  & 
! 4s,5s will be treated as 3s
!         k-rn , no f-elements
          54*3/
! ...
real(kind=8) xza,xzb
integer iat,jat,iz(*)

       za=xza
       zb=xzb
       na=iz(iat)
       nb=iz(jat)
debug=.false.
!debug=.true.

! ii selects kind of ovl by multiplying the shell
! kind    <1s|1s>  <2s|1s>  <2s|2s>  <1s|3s>  <2s|3s>  <3s|3s>
! case:      1        2        4       3        6         9
!
ii=shell(na)*shell(nb)
if(debug) write(*,*) 'shell', ii

R05=R*0.5
ax=(za+zb)*R05
bx=(zb-za)*R05

! same elements
if(za.eq.zb.OR.abs(za-zb).lt.0.1) then
  select case (ii)
   case (1)
    ovl=0.25d0*sqrt((za*zb*R*R)**3)*(A2(ax)*Bint(bx,0)-Bint(bx,2)*A0(ax))
   case (2)  
    ovl = SQRT(1.D00/3.D00)
    if(shell(na).lt.shell(nb)) then
    ! <1s|2s>
      norm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125D00
      ovl=ovl*norm*(A3(ax)*Bint(bx,0)-Bint(bx,3)*A0(ax)+A2(ax)*Bint(bx,1)-Bint(bx,2)*A1(ax))
     else
    ! switch za/zb to get <2s|1s>
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125D00
      ovl=ovl*norm*(A3(ax)*Bint(bx,0)-Bint(bx,3)*A0(ax)+A2(ax)*Bint(bx,1)-Bint(bx,2)*A1(ax))
    endif
   case (4)
    norm=SQRT((ZA*ZB)**5)*(R**5)*0.0625d0
    ovl=norm* (A4(ax)*Bint(bx,0)+Bint(bx,4)*A0(ax)-2.0d0*A2(ax)*Bint(bx,2))*(1d0/3d0)
   case(3) 
    if(shell(na).lt.shell(nb)) then
      norm=SQRT((ZA**3)*(ZB**7)/7.5d00)*(R**5)*0.0625d00
      ovl=norm*(A4(ax)*Bint(bx,0)-Bint(bx,4)*A0(ax)+2.d0*(A3(ax)*Bint(bx,1)-Bint(bx,3)*A1(ax)))/sqrt(3.d0)
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((ZA**3)*(ZB**7)/7.5d00)*(R**5)*0.0625d00
      ovl=norm*(A4(ax)*Bint(bx,0)-Bint(bx,4)*A0(ax)+2.d0*(A3(ax)*Bint(bx,1)-Bint(bx,3)*A1(ax)))/sqrt(3.d0)
    endif
   case(6) 
    if(shell(na).lt.shell(nb)) then
      norm=SQRT((za**5)*(zb**7)/7.5D00)*(R**6)*0.03125D00
      ovl=norm*(A5(ax)*Bint(bx,0)+A4(ax)*Bint(bx,1)-2d0*(A3(ax)*Bint(bx,2)+A2(ax)*Bint(bx,3))+A1(ax)*Bint(bx,4)+A0(ax)*Bint(bx,5))/3.d0
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((za**5)*(zb**7)/7.5D00)*(R**6)*0.03125D00
      ovl=norm*(A5(ax)*Bint(bx,0)+A4(ax)*Bint(bx,1)-2d0*(A3(ax)*Bint(bx,2)+A2(ax)*Bint(bx,3))+A1(ax)*Bint(bx,4)+A0(ax)*Bint(bx,5))/3.d0
    endif
   case(9)
      norm=sqrt((ZA*ZB*R*R)**7)/480.d0
      ovl=norm*(A6(ax)*Bint(bx,0)-3.d0*(A4(ax)*Bint(bx,2)-A2(ax)*Bint(bx,4))-A0(ax)*Bint(bx,6))/3.D00
   end select
else ! different elements
   select case (ii)
   case (1)
      norm=0.25d0*sqrt((za*zb*R*R)**3)
      ovl=(A2(ax)*B0(bx)-B2(bx)*A0(ax))*norm
   case (2)
      ovl = SQRT(1.D00/3.D00)
    if(shell(na).lt.shell(nb)) then
    ! <1s|2s>
      norm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125D00 
      ovl=ovl*norm*(A3(ax)*B0(bx)-B3(bx)*A0(ax)+A2(ax)*B1(bx)-B2(bx)*A1(ax))
     else
    ! switch za/zb to get <2s|1s>
      xx=za
      za=zb
      zb=xx 
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((ZA**3)*(ZB**5))*(R**4)*0.125D00 
      ovl=ovl*norm*(A3(ax)*B0(bx)-B3(bx)*A0(ax)+A2(ax)*B1(bx)-B2(bx)*A1(ax))
    endif
   case (4) ! <2s|2s>
      norm=SQRT((ZA*ZB)**5)*(R**5)*0.0625D00
      ovl=norm* (A4(ax)*B0(bx)+B4(bx)*A0(ax)-2.0D00*A2(ax)*B2(bx))*(1d0/3d0)
   case(3)  ! <1s|3s> + <3s|1s>
    if(shell(na).lt.shell(nb)) then
      norm=SQRT((ZA**3)*(ZB**7)/7.5d00)*(R**5)*0.0625d00
      ovl=norm*(A4(ax)*B0(bx)-B4(bx)*A0(ax)+2.d0*(A3(ax)*B1(bx)-B3(bx)*A1(ax)))/sqrt(3.d0)
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((ZA**3)*(ZB**7)/7.5d00)*(R**5)*0.0625d00
      ovl=norm*(A4(ax)*B0(bx)-B4(bx)*A0(ax)+2.d0*(A3(ax)*B1(bx)-B3(bx)*A1(ax)))/sqrt(3.d0)
    endif
   case(6)  ! <2s|3s> + <3s|2s>
    if(shell(na).lt.shell(nb)) then
      norm=SQRT((za**5)*(zb**7)/7.5D00)*(R**6)*0.03125D00
      ovl=norm*(A5(ax)*B0(bx)+A4(ax)*B1(bx)-2d0*(A3(ax)*B2(bx)+A2(ax)*B3(bx))+A1(ax)*B4(bx)+A0(ax)*B5(bx))/3.d0
    else
      xx=za
      za=zb
      zb=xx
      ax=(za+zb)*R05
      bx=(zb-za)*R05
      norm=SQRT((za**5)*(zb**7)/7.5D00)*(R**6)*0.03125D00
      ovl=norm*(A5(ax)*B0(bx)+A4(ax)*B1(bx)-2d0*(A3(ax)*B2(bx)+A2(ax)*B3(bx))+A1(ax)*B4(bx)+A0(ax)*B5(bx))/3.d0
    endif
    case(9) ! <3s|3>
      norm=sqrt((ZA*ZB*R*R)**7)/1440.d0
!      ovl=norm*(A6(ax)*B0(bx)-3.d0*(A4(ax)*B2(bx)-A2(ax)*B4(bx))-A0(ax)*Bint(bx,6))
      ovl=norm*(A6(ax)*B0(bx)-3.d0*(A4(ax)*B2(bx)-A2(ax)*B4(bx))-A0(ax)*B6(bx))
   end select
endif
return
end subroutine ssovl


!****************************************
!* A(x) auxiliary integrals             *
!* Quantenchemie - Ein Lehrgang Vol 5   *
!* p. 570  eq. 11.4.14                  *           
!****************************************

real(8) pure function A0(x)
! Hilfsintegral A_0
implicit none
real(8), intent(in) :: x
A0=exp(-x)/x
return 
end function

real(8) pure function A1(x)
! Hilfsintegral A_1
implicit none
real(8), intent(in) :: x
A1=((1+x)*exp(-x))/(x**2)
return
end function


real(8) pure function A2(x)
! Hilfsintegral A_2
implicit none
real(8), intent(in) :: x
A2=((2d0+2d0*x+x**2)*exp(-x))/x**3
return
end function


real(8) pure function A3(x)
! Hilfsintegral A_3
implicit none
real(8), intent(in) :: x
real(8) xx
real(8) x2,x3,x4
x2=x*x
x3=x2*x
x4=x3*x
xx=(6d0+6d0*x+3d0*x2+x3)
A3=(xx*exp(-x))/x4
return
end function


real(8) pure function A4(x)
! Hilfsintegral A_4
implicit none
real(8), intent(in) :: x
real(8) xx
real(8) x2,x3,x4,x5
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
xx=(24d0+24d0*x+12d0*x2+4d0*x3+x4)
A4=(xx*exp(-x))/x5
return
end function

real(8) pure function A5(x)
! Hilfsintegral A_5
implicit none
real(8), intent(in) :: x
real(8) xx
real(8) x2,x3,x4,x5,x6
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
xx=(120d0+120d0*x+60d0*x2+20d0*x3+5d0*x4+x5)
A5=(xx*exp(-x))/x6
return
end function

real(8) pure function A6(x)
! Hilfsintegral A_6
implicit none
real(8), intent(in) :: x
real(8) xx
real(8) x2,x3,x4,x5,x6,x7
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
x7=x6*x
xx=(720d0+720d0*x+360d0*x2+120d0*x3+30d0*x4+6d0*x5+x6)
A6=(xx*exp(-x))/x7
return
end function


!**************************************
!* B(x) auxiliary integrals           *
!* Quantenchemie - Ein Lehrgang Vol 5 *
!* p. 570  eq. 11.4.14b               *
!**************************************


real(8) pure function B0(x)
! Hilfsintegral B_0
implicit none
real(8), intent(in) :: x
B0=(exp(x)-exp(-x))/x
return
end function

real(8) pure function B1(x)
! Hilfsintegral B_1
implicit none
real(8), intent(in) :: x
real(8) x2,x3
x2=x*x
x3=x2*x
B1=((1d0-x)*exp(x)-(1d0+x)*exp(-x))/x2
return
end function

real(8) pure function B2(x)
! Hilfsintegral B_2
implicit none
real(8), intent(in) :: x
real(8) x2,x3
x2=x*x
x3=x2*x
B2=(((2d0-2*x+x2)*exp(x)) - ((2d0+2d0*x+x2)*exp(-x)))/x3
return
end function

real(8) pure function B3(x)
! Hilfsintegral B_3
implicit none
real(8), intent(in) :: x
real(8) xx,yy
real(8) x2,x3,x4
x2=x*x
x3=x2*x
x4=x3*x
xx=(6d0-6d0*x+3d0*x2-x3)*exp(x)/x4
yy=(6d0+6d0*x+3d0*x2+x3)*exp(-x)/x4
B3=xx-yy
return
end function


real(8) pure function B4(x)
! Hilfsintegral B_4
implicit none
real(8), intent(in) :: x
real(8) xx,yy
real(8) x2,x3,x4,x5
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
xx=(24d0-24d0*x+12d0*x2-4d0*x3+x4)*exp(x)/x5
yy=(24d0+24d0*x+12d0*x2+4d0*x3+x4)*exp(-x)/x5
B4=xx-yy
return
end function

real(8) pure function B5(x)
! Hilfsintegral B_5
implicit none
real(8), intent(in) :: x
real(8) xx,yy
real(8) x2,x3,x4,x5,x6
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
xx=(120d0-120*x+60*x2-20*x3+5*x4-x5)*exp(x)/x6
yy=(120d0+120*x+60*x2+20*x3+5*x4+x5)*exp(-x)/x6
B5=xx-yy
return
end function

real(8) function B6(x)
! Hilfsintegral B_6
implicit none
real(8), intent(in) :: x
real(8) x2,x3,x4,x5,x6,x7,yy,xx
x2=x*x
x3=x2*x
x4=x3*x
x5=x4*x
x6=x5*x
x7=x6*x
xx=(720d0 - 720d0*x+ 360d0*x2 - 120d0*x3 + 30d0*x4 - 6d0*x5 + x6)*exp(x)/x7
yy=(720d0 + 720d0*x + 360d0*x2 + 120d0*x3 + 30d0*x4 + 6d0*x5 + x6)*exp(-x)/x7
B6=xx-yy
return
end function


real*8 function bint(x,k)
! calculates B_k(x)
! general summation formula
! 'infinite' sum is numerically unstable. 12 terms seem
! accurate enough
implicit none
real(8), intent(in) :: x
real(8) xx,yy
integer, intent(in) :: k
integer i
integer(8) fact
bint=0

if(abs(x).lt.1e-6) then
do i=0,k
   bint=(1.d0+(-1d0)**i)/(dble(i)+1.d0)
end do
return
endif

do i=0,12
xx=1d0-((-1d0)**(k+i+1))
yy=dble(fact(i))*dble((k+i+1))
bint=bint+xx/yy*(-x)**i
enddo


end function bint


! faculty function
integer(8) function fact(N)
implicit none
integer j,n
fact=1
do j=2,n
  fact=fact*j
enddo
return 

end




subroutine gsovl(r,iat,jat,iz,xza,xzb,g)
! GRADIENT 
! calculates the s-type overlap integral over 1s,2s,3s slater functions                  
! ovl = overlap integral                                                                  
! za  = slater exponent atom A                                                            
! zb  = slater exponent atom B                                                            
! R   = distance between atom A and B                                                     
implicit none                                                                             
integer ii,shell(72)                                                              
logical debug                                                                             
real(8) ax,bx,R05,za,zb,R                                                             
integer na,nb                                                                             
data shell/                 &                                                             
!          h, he                                                                          
          1,1               &                                                             
!         li-ne                                                                           
          ,2,2,2,2,2,2,2,2, &                                                             
!         na-ar                                                                           
          3,3,3,3,3,3,3,3,  &                                                              
! 4s,5s will be treated as 3s
!         k-rn , no f-elements
          54*3/
! ...                                                                                     
real*8 g,Fa,Fb
!--------------------- set exponents ---------------------------------------
real(kind=8) xza,xzb
real(kind=8) xx
integer iat,jat,iz(*)
logical lsame

       za=xza
       zb=xzb
       na=iz(iat)
       nb=iz(jat)
!----------------------------------------------------------------------------

debug=.false.
!debug=.true.

! ii selects kind of ovl by multiplying the shell
! kind    <1s|1s>  <2s|1s>  <2s|2s>  <1s|3s>  <2s|3s>  <3s|3s>
! case:      1        2        4       3        6         9   
!                                                             
ii=shell(na)*shell(nb)                                        
if(debug) write(*,*) 'gshell', ii                              
R05=R*0.5
ax=(za+zb)*R05
Fa=(za+zb)
bx=(zb-za)*R05
Fb=(zb-za)
lsame=.false.
!         
! same elements
if(za.eq.zb.OR.abs(za-zb).lt.0.1) then
lsame=.true.
! call arguments: gtype(exponents,argumentDeriv.,distance,gradient,(Switch shell),sameElement)
  select case (ii)                                      
   case (1)                                             
     call g1s1s(za,zb,Fa,Fb,R,g,lsame)
   case (2)             
    if(shell(na).lt.shell(nb)) then
      call g2s1s(za,zb,Fa,Fb,R,g,.false.,lsame)
     else
      xx=za
      za=zb
      zb=xx
      call g2s1s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case (4)                    
     call g2s2s(za,zb,Fa,Fb,R,g,lsame)
   case(3)                    
    if(shell(na).lt.shell(nb)) then
    call g1s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else
      xx=za
      za=zb
      zb=xx
    call g1s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case(6) 
    if(shell(na).lt.shell(nb)) then
    call g2s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else
      xx=za
      za=zb
      zb=xx
    call g2s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif
   case(9)                                                         
    call g3s3s(za,zb,Fa,Fb,R,g,lsame)
   end select                                                                
else ! different elements
lsame=.false.                                                    
   select case (ii)                                                          
   case (1)                                                                  
     call g1s1s(za,zb,Fa,Fb,R,g,lsame)
   return
   case (2)  ! <1s|2s>                                                                  
    if(shell(na).lt.shell(nb)) then                                          
      call g2s1s(za,zb,Fa,Fb,R,g,.false.,lsame)
     else                                                                    
      xx=za                                                                  
      za=zb                                                                  
      zb=xx                                                                  
      call g2s1s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif                                                                    
   case (4) ! <2s|2s>                                                        
      call g2s2s(za,zb,Fa,Fb,R,g,lsame)
   case(3)  ! <1s|3s> + <3s|1s>                                              
    if(shell(na).lt.shell(nb)) then                                          
    call g1s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else                                                                                  
      xx=za                                                                               
      za=zb                                                                               
      zb=xx                                                                               
    call g1s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif                                                                                 
   case(6)  ! <2s|3s> + <3s|2s>                                                           
    if(shell(na).lt.shell(nb)) then                                                       
    call g2s3s(za,zb,Fa,Fb,R,g,.false.,lsame)
    else                                                                                                       
      xx=za                                                                                                    
      za=zb                                                                                                    
      zb=xx                                                                                                    
    call g2s3s(za,zb,Fa,Fb,R,g,.true.,lsame)
    endif                                                                                                      
    case(9) ! <3s|3>                                                                                           
    call g3s3s(za,zb,Fa,Fb,R,g,lsame)
   end select                                                                                                  
endif                                                                                                          

return
end subroutine gsovl


!-------------------------------------------------------------
! Maple was used to find the analy. derivatives of
! the slater integrals (expressions over A,B aux. integrals) 
! Optimized fortran code by maple with some human-corrections
!-------------------------------------------------------------
subroutine g1s1s(za,zb,Fa,Fb,R,g,sameElement)
! slater overlap derv.
! derivative of explicit integral expression
! using maple 
implicit real(8) (t)
real(8) za,zb,Fa,Fb
real(8) g,R
logical sameElement

if(sameElement) then
  t1 = za ** 2
  t3 = zb ** 2
  t5 = t1 * za * t3 * zb
  t6 = R ** 2
  t7 = t6 ** 2
  t10 = Fa * R
  t14 = exp(-0.5D0 * t10)
  t17 = sqrt(t5 * t7 * t6)
  g = -(1d0/3d0) * t5 * t7 / Fa * (0.2D1 + t10) * t14 / t17
  return
else
  t1 = za ** 2
  t3 = zb ** 2
  t5 = t1 * za * t3 * zb
  t6 = Fb ** 2
  t7 = Fb * R
  t8 = 0.5D0 * t7
  t9 = exp(t8)
  t12 = exp(-t8)
  t15 = t6 * Fa
  t22 = Fa ** 2
  t23 = t22 * t9
  t27 = t22 * t12
  t31 = t6 * Fb
  t32 = R * t31
  t37 = t22 * Fa
  t38 = R * t37
  t43 = R ** 2
  t44 = t43 * t31
  t51 = t43 * t37
  t56 = 0.4D1 * t6 * t9 - 0.4D1 * t6 * t12 + 0.2D1 * t15 * R * t9 -          &
  0.2D1 * t15 * R * t12 - 0.4D1 * t23 + 0.2D1 * t23 * t7 + 0.4D1 * t27       &
  + 0.2D1 * t27 * t7 - 0.2D1 * t32 * t9 - 0.2D1 * t32 * t12 - 0.2D1 * t38 *  &
  t9 + 0.2D1 * t38 * t12 - 0.1D1 * t44 * Fa * t9 - 0.1D1                     &
  * t44 * Fa * t12 + t51 * t9 * Fb + t51 * t12 * Fb
  t61 = exp(-0.5D0 * Fa * R)
  t62 = t43 ** 2
  t65 = sqrt(t5 * t62 * t43)
  g = -0.2D1 * t5 * R * t56 * t61 / t65 / t31 / t37
  return
endif


end subroutine g1s1s


subroutine g2s1s(za,zb,Fa,Fb,R,g,switch,lsame)
! slater overlap derv.
! derivative of explicit integral expression
! using maple 
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical switch
logical lsame
norm=(1d0/24d0)*sqrt(za**3*zb**5*3d0)

if(switch) then
Fb=-Fb
endif

if(lsame) then

      t1 = Fa * R
      t3 = exp(-0.5000000000D0 * t1)
      t6 = Fa ** 2
      t7 = R ** 2
      g = -0.1000000000D-8 * R * t3 * (0.5333333380D10 + 0.2666666670D10 &
      * t1 + 0.1333333333D10 * t6 * t7) / t6
      g=g*norm

else

      t3 = exp(-0.5000000000D0 * Fa * R)
      t4 = Fa ** 2
      t5 = t4 * Fa
      t6 = Fb * R
      t7 = 0.5000000000D0 * t6
      t8 = exp(t7)
      t9 = t5 * t8
      t11 = Fb ** 2
      t12 = t11 * Fa
      t15 = exp(-t7)
      t18 = t4 ** 2
      t19 = R * t18
      t22 = t11 ** 2
      t29 = Fb * t4
      t36 = R ** 2
      t37 = t36 * t18
      t44 = t36 * R
      t48 = -0.12D2 * t9 + 0.4D1 * t12 * t8 - 0.4D1 * t12 * t15 &
      - 0.6D1 * t19 * t8 - 0.6D1 * t22 * t8 * R - 0.6D1 * t22 * t15 &
      * R + 0.4D1 * t29 * t15 - 0.4D1 * t29 * t8 + 0.6D1 * t19 * t15 &
      + 0.2D1 * t37 * t8 * Fb + 0.4D1 * t37 * t15 * Fb + t44 * t18 * t15 * t11
      t49 = t5 * t15
      t51 = t11 * Fb
      t58 = t51 * Fa
      t59 = R * t8
      t76 = t36 * t15
      t79 = t22 * Fa
      t87 = 0.12D2 * t49 - 0.12D2 * t51 * t15 - 0.1D1 * t22 * t4 * t15 * t44 &
      + 0.4D1 * t58 * t59 - 0.8D1 * t58 * R * t15 + 0.4D1 * t9 * t6 + 0.8D1 * &
      t49 * t6 + 0.2D1 * t49 * t11 * t36 + 0.4D1 * t11 * t4 * t59 - 0.2D1 * t51 &
      * t4 * t76 - 0.2D1 * t79 * t36 * t8 - 0.4D1 * t79 * t76 + 0.12D2 * t51 * t8
      g = -0.16D2 * t3 * (t48 + t87) / t36 / t22 / t18
      g=g*norm
endif

return
end subroutine g2s1s

subroutine g2s2s(za,zb,Fa,Fb,R,g,SameElement)
! slater overlap derv.
! derivative of explicit integral expression
! using maple 
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical SameElement

norm=1d0/(16d0*3d0)*SQRT((ZA*ZB)**5)

if(SameElement) then

      t2 = R ** 2
      t5 = Fa ** 2
      t9 = t5 * Fa
      t10 = t2 ** 2
      t16 = exp(-Fa * R / 0.2D1)
      g = (-0.4266666666D2 * R - 0.2133333333D2 * Fa * t2 - 0.2133333333D1 &
          * t5 * t2 * R - 0.1066666666D1 * t9 * t10) * t16 / t9
      g=g*norm


return
else
      t1 = R ** 2
      t3 = 0.3840000000D3 * t1 * Fb
      t4 = t1 * R
      t5 = Fb ** 2
      t7 = 0.6400000000D2 * t4 * t5
      t8 = 0.7680000000D3 * R
      t10 = Fa ** 2
      t11 = t10 ** 2
      t12 = t11 * Fa
      t14 = Fb * R
      t15 = 0.768000000D3 * t14
      t17 = 0.1280000000D3 * t5 * t1
      t21 = 0.256000000D3 * t5 * R
      t22 = t5 * Fb
      t24 = 0.1280000000D3 * t22 * t1
      t26 = t10 * Fa
      t28 = t5 ** 2
      t30 = 0.1280000000D3 * t1 * t28
      t32 = 0.256000000D3 * t22 * R
      t33 = 0.512000000D3 * t5
      t34 = t28 * Fb
      t36 = 0.6400000000D2 * t4 * t34
      t40 = 0.768000000D3 * t28 * R
      t42 = 0.3840000000D3 * t1 * t34
      t45 = 0.1536000000D4 * t28
      t47 = 0.7680000000D3 * t34 * R
      t51 = exp(-0.5D0 * Fa * R)
      t53 = 0.5D0 * t14
      t54 = exp(-t53)
      t68 = exp(t53)
      g = (((t3 + t7 + t8) * t12 + (0.1536000000D4 + t15 + t17) * t11 + &
      (-t21 - t24) * t26 + (t30 - t32 - t33 + t36) * t10 + (t40 + t42) *&
      Fa + t45 + t47) * t51 * t54 + ((t3 - t8 - t7) * t12 + (-0.1536000000D4 &
      + t15 - t17) * t11 + (-t24 + t21) * t26 + (-t30 + t33 - t32 &
      + t36) * t10 + (-t40 + t42) * Fa + t47 - t45) * t51 * t68) / t1 /  &
      t12 / t34


      g=g*norm


return
endif

end subroutine g2s2s






subroutine g1s3s(za,zb,Fa,Fb,R,g,switch,lsame)
! slater overlap derv.
! derivative of explicit integral expression
! using maple
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical switch
logical lsame

if(switch) Fb=-Fb

norm=SQRT((ZA**3)*(ZB**7)/7.5d00)/(16d0*sqrt(3d0))

if(lsame) then

  t1 = Fa * R
  t3 = exp(-0.5000000000D0 * t1)
  t4 = R ** 2
  g = -0.1600000000D1 * t3 * t4 * R * (0.2D1 + t1) / Fa
  g=g*norm
else

      t3 = exp(-0.5000D0 * Fa * R)
      t4 = Fb ** 2
      t5 = t4 ** 2
      t6 = t5 * Fb
      t7 = t6 * Fa
      t8 = R ** 2
      t9 = Fb * R
      t10 = 0.50D0 * t9
      t11 = exp(t10)
      t15 = exp(-t10)
      t16 = t8 * t15
      t19 = Fa ** 2
      t21 = t8 * R
      t22 = t21 * t15
      t25 = t19 * Fa
      t27 = t8 ** 2
      t31 = t19 ** 2
      t32 = t31 * Fa
      t33 = t8 * t32
      t45 = t4 * Fb
      t48 = t31 * t15
      t55 = t4 * t25
      t56 = t11 * R
      t59 = t15 * R
      t62 = t5 * Fa
      t73 = -0.6D1 * t7 * t8 * t11 - 0.18D2 * t7 * t16 - 0.6D1 * t6 * t19 &
      * t22 - 0.1D1 * t6 * t25 * t27 * t15 + 0.6D1 * t33 * t11 * Fb + 0.18D2 &
      * t33 * t15 * Fb + 0.6D1 * t21 * t32 * t15 * t4 + t27 * t32* t15 * t45 &
      + 0.2D1 * t48 * t45 * t21 + 0.12D2 * t48 * t4 * t8 + 0.12D2 * t55 * t56 &
      + 0.12D2 * t55 * t59 + 0.12D2 * t62 * t56 - 0.36D2 * t62 * t59 - 0.12D2 &
      * t5 * t19 * t16 - 0.2D1 * t5 * t25 * t22
      t74 = t31 * t11
      t79 = t45 * t19
      t92 = R * t32
      t95 = t45 * Fa
      t100 = Fb * t25
      t111 = 0.12D2 * t74 * t9 + 0.36D2 * t48 * t9 + 0.12D2 * t79 * t56 - 0.12D2  &
      * t79 * t59 + 0.48D2 * t5 * t11 - 0.24D2 * t6 * t11 * R - 0.24D2 * t6 * t15 &
      * R - 0.24D2 * t92 * t11 + 0.24D2 * t95 * t11 - 0.24D2 * t95 * t15 + 0.24D2 &
      * t100 * t15 + 0.24D2 * t92 * t15 - 0.24D2 * t100 * t11 - 0.48D2 * t5 * t15 &
      - 0.48D2 * t74 + 0.48D2 * t48
      g = -0.32D2 * t3 * (t73 + t111) / t8 / t6 / t32

      g=g*norm
endif

return
end subroutine g1s3s


subroutine g2s3s(za,zb,Fa,Fb,R,g,switch,lsame)
! slater overlap derv.
! derivative of explicit integral expression
! using maple
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical switch
logical lsame
norm=sqrt((za**5)*(zb**7)/7.5D00)/96.d0

if(switch) Fb=-Fb

if(lsame) then
      t1 = Fa * R
      t3 = exp(-0.5000000000D0 * t1)
      t6 = Fa ** 2
      t7 = R ** 2
      t14 = t6 ** 2
      t15 = t7 ** 2
      g = -0.2000000000D-8 * R * t3 * (0.1280000000D12 + 0.6400000000D11 &
      * t1 + 0.1280000000D11 * t6 * t7 + 0.1066666670D10 * t6 * Fa * t7 &
      * R + 0.533333333D9 * t14 * t15) / t14
      g=g*norm
else

      t3 = exp(-0.5D0 * Fa * R)
      t4 = Fb ** 2
      t5 = t4 ** 2
      t6 = Fa ** 2
      t7 = t6 * Fa
      t8 = t5 * t7
      t9 = R ** 2
      t11 = 0.50D0 * Fb * R
      t12 = exp(t11)
      t13 = t9 * t12
      t16 = t6 ** 2
      t17 = t16 * Fa
      t18 = exp(-t11)
      t21 = t5 * Fb
      t28 = t9 * t18
      t32 = t9 * R
      t33 = t32 * t18
      t36 = t5 * t4
      t38 = t9 ** 2
      t39 = t38 * t18
      t41 = t21 * Fa
      t42 = R * t12
      t45 = t16 * t6
      t46 = t4 * Fb
      t49 = t46 * t16
      t52 = -0.6D1 * t8 * t13 + 0.120D3 * t17 * t18 + 0.120D3 * t21 * t18 &
       - 0.120D3 * t17 * t12 - 0.120D3 * t21 * t12 - 0.6D1 * t8 * t28 - 0.2D1 &
      * t5 * t16 * t33 + t36 * t7 * t39 - 0.48D2 * t41 * t42 + t45 * t46 * t39 - 0.6D1 * t49 * t13
      t54 = R * t18
      t60 = t46 * t6
      t63 = Fb * t16
      t66 = t5 * Fa
      t69 = t4 * t7
      t72 = t36 * t6
      t75 = t32 * t12
      t78 = Fb * t9
      t84 = Fb * t17
      t87 = -0.24D2 * t46 * t7 * t54 - 0.24D2 * t5 * t6 * t42 + 0.24D2 *&
       t60 * t12 + 0.24D2 * t63 * t18 - 0.24D2 * t66 * t12 - 0.24D2 * t69 &
      * t18 + 0.9D1 * t72 * t33 + 0.3D1 * t72 * t75 + 0.24D2 * t78 * t45 &
      * t12 - 0.6D1 * t49 * t28 + 0.48D2 * t84 * t42
      t102 = t21 * t6
      t105 = t4 * t17
      t113 = t45 * t4
      t118 = 0.72D2 * t84 * t54 + 0.72D2 * t41 * t54 + 0.36D2 * t78 * t45 &
      * t18 + 0.2D1 * t46 * t17 * t33 + 0.24D2 * t4 * t16 * t42 - 0.6D1   &
      * t102 * t13 - 0.6D1 * t105 * t13 + 0.18D2 * t105 * t28 + 0.2D1 &
      * t21 * t7 * t33 - 0.3D1 * t113 * t75 + 0.9D1 * t113 * t33
      t121 = t36 * Fa
      t130 = R * t45
      t145 = 0.18D2 * t102 * t28 + 0.24D2 * t121 * t13 + 0.36D2 * t121 * &
       t28 - 0.24D2 * t60 * t18 - 0.24D2 * t63 * t12 + 0.60D2 * t130 * t18 &
      + 0.60D2 * t36 * t18 * R + 0.24D2 * t69 * t12 + 0.60D2 * t36 * t12 * &
      R - 0.60D2 * t130 * t12 + 0.24D2 * t66 * t18
      g = 0.128D3 * t3 * (t52 + t87 + t118 + t145) / t9 / t36 / t45

  g=g*norm
endif

return
end subroutine g2s3s


subroutine g3s3s(za,zb,Fa,Fb,R,g,SameElement)
! slater overlap derv.
! derivative of explicit integral expression
! using maple 
implicit double precision (t)
real(8) za,zb,Fa,Fb
real(8) g,R,norm
logical SameElement

norm=sqrt((ZA*ZB)**7)/1440.d0

if(SameElement) then

      t1 = Fa * R
      t3 = exp(-0.5000000000D0 * t1)
      t5 = Fa ** 2
      t6 = t5 ** 2
      t7 = t6 * Fa
      t8 = R ** 2
      t9 = t8 ** 2
      g = -0.2000000000D-8 * t3 * R * (0.457142857D9 * t7 * t9 &
      * R + 0.7680000000D12 * t1 + 0.1536000000D12 * t5 * t8 &
      + 0.1280000000D11 * t5 * Fa * t8 * R + 0.914285715D9 * t6 * t9 + 0.1536000000D13) / t7


      g=g*norm
return
else


      t3 = exp(-0.5000000000D0 * Fa * R)
      t4 = Fa ** 2
      t5 = t4 ** 2
      t6 = t5 * t4
      t7 = Fb * R
      t8 = 0.5000000000D0 * t7
      t9 = exp(-t8)
      t10 = t6 * t9
      t13 = Fb ** 2
      t14 = t13 * Fb
      t15 = t13 ** 2
      t16 = t15 * t14
      t17 = R ** 2
      t18 = t17 * R
      t19 = t16 * t18
      t23 = exp(t8)
      t24 = t6 * t23
      t27 = t5 * Fa
      t28 = t27 * t13
      t29 = R * t23
      t32 = t6 * t13
      t33 = t17 * t9
      t36 = t15 * Fb
      t37 = t4 * t36
      t38 = t9 * R
      t43 = t17 * t23
      t46 = t4 * Fa
      t47 = t5 * t46
      t48 = t47 * t18
      t52 = t47 * t17
      t65 = 0.120D3 * t10 * t7 - 0.12D2 * t19 * t4 * t9 + 0.120D3 &
      * t24 * t7 + 0.24D2 * t28 * t29 + 0.24D2 * t32 * t33 + 0.24D2 * t37 &
      * t38 - 0.24D2 * t28 * t38 - 0.24D2 * t32 * t43 - 0.12D2 * t48 * t13 &
      * t23 + 0.60D2 * t52 * t23 * Fb + 0.12D2 * t48 * t13 * t9 + 0.60D2 &
      * t52 * t9 * Fb - 0.12D2 * t19 * t4 * t23
      t66 = t17 ** 2
      t67 = t16 * t66
      t74 = t27 * t14
      t77 = t6 * t14
      t78 = t18 * t23
      t81 = t46 * t15
      t86 = t27 * t15
      t89 = t5 * t36
      t90 = t18 * t9
      t97 = t46 * t36
      t104 = -0.1D1 * t67 * t46 * t9 - 0.1D1 * t67 * t46 * t23 - 0.12D2 &
      * t74 * t43 + 0.2D1 * t77 * t78 - 0.24D2 * t81 * t29 + 0.24D2 * t81 &
      * t38 + 0.2D1 * t86 * t78 + 0.2D1 * t89 * t90 - 0.2D1 * t86 * t90 &
      + 0.24D2 * t37 * t29 + 0.12D2 * t97 * t33 + 0.2D1 * t89 * t78 - 0.12D2 * t74 * t33
      t108 = t5 * t14
      t111 = t15 * t13
      t112 = t111 * t4
      t117 = t111 * t46
      t122 = t111 * Fa
      t129 = t4 * t15
      t132 = t47 * R
      t139 = 0.2D1 * t77 * t90 - 0.24D2 * t108 * t38 + 0.24D2 * t112 * t43 &
      - 0.24D2 * t112 * t33 + 0.2D1 * t117 * t78 - 0.2D1 * t117 * t90 + 0.120D3 &
      * t122 * t29 - 0.120D3 * t122 * t38 + 0.12D2 * t97 * t43 - 0.48D2 * t129 &
      * t23 + 0.120D3 * t132 * t9 - 0.120D3 * t132 * t23 + 0.240D3 * t111 * t23
      t140 = t47 * t66
      t145 = t16 * R
      t150 = t16 * t17
      t160 = t5 * t13
      t170 = t140 * t14 * t23 + t140 * t14 * t9 - 0.120D3 * t145 * t9 - 0.24D2 &
      * t108 * t29 - 0.60D2 * t150 * Fa * t23 - 0.240D3 * t111 * t9 - 0.240D3 &
      * t24 + 0.240D3 * t10 + 0.48D2 * t129 * t9 - 0.48D2 * t160 * t9 + 0.48D2 &
      * t160 * t23 - 0.120D3 * t145 * t23 - 0.60D2 * t150 * Fa * t9
      g = -0.768D3 * t3 * (t65 + t104 + t139 + t170) / t17 / t47 / t16

      g=g*norm
return
endif

end subroutine g3s3s

subroutine done(aa,io)
! normal program termination
implicit none
integer io
character(*)aa

 write(io,'(a)') ' '
 write(io,'(5x,''normal termination of '',a)') adjustl(trim(aa))
 stop
end subroutine done

!************************************************************
!* reads a turbomole (bohr) or xmol (angst)rom file.        *
!* Tests if xmol starts with "number of atoms + blank" or   *
!* directly with the coordinates.                           *
!************************************************************
subroutine tmolrd(xyz,iat,nat,infile,echo)
!use io_reader
implicit none                        
character*2 cc         
character*80  atmp           
character*(*) infile                 
real(kind=8) xyz(3,10000),xx(5)            
integer iat(10000),nat,i,nn,j    
real(kind=8) bohr                          
logical da,echo             

bohr=0.52917726d0
i=0            
inquire(file=infile,exist=da)
select case (da)             
case (.true.)                
      if(echo) write(*,'('' reading...'',$)')

open(unit=3,file=infile)
! test for tmol or xyz file

 read(3,'(a)') atmp ! $coord
rewind(3)                   
if(index(atmp,'$coord').ne.0) then

 ! count number of atoms
 do while (da)          
  read(3,'(a)',end=100) atmp ! $coord
   if(index(atmp,'$coord').ne.0) cycle
   if(index(atmp,'$').ne.0) exit      
   i=i+1                              
  enddo                               
 nat=i                                
 100 continue                         
 
 rewind(unit=3)
 
 ! read TMOL file
 read(3,*) atmp ! $coord
 do j=1,nat             
   read(3,*) xyz(1,j),xyz(2,j),xyz(3,j),cc
   call elem(cc,iat(j))                   
   xyz(1:3,j)=xyz(1:3,j)             
  enddo                                                 
 if(echo) write(*,*) ' Turbomole file [bohr] :  ', trim(infile)          

 close(3)

 else ! xyz file
       read(3,'(a)',end=101) atmp
! check for first two lines
       call readl(atmp,xx,nn)
        if(nn.gt.1) then   ! more than one argument found, assuming they are coords
           do 
            nat=nat+1
            read(3,'(a)',end=123) atmp
           enddo
          else
            nat=idint(xx(1))
           read(3,'(a)',end=101) atmp  !titel line
        endif
 123   if(nn.gt.1) rewind(3)
       do i=1,nat
            read(3,'(a)') atmp
            call readl(atmp,xx,nn)
            call elem(atmp,iat(i))
            xyz(1:3,i)=xx(1:3)*1d0/bohr
       enddo
 101  close(3)
      if(echo) write(*,'(5x,'' XYZ file [angst]: '',a)')  trim(infile)
      endif

case (.false.)
  write(*,*) ' no input file <',trim(infile) ,'> found !! '
end select                                                 
end subroutine      

!********************************
!* convert a word to lower case *
!********************************
      subroutine lower_case(word)
      character (len=*) , intent(in out) :: word
      integer :: i,ic,nlen
      nlen = len(word)
      do i=1,nlen
      ic = ichar(word(i:i))
      if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
      end do
      end subroutine lower_case


!********************************************
!* split a string s into n separate words w *
!********************************************
subroutine charsplit(s,n,w)
implicit none
integer i,n,k
character*80, intent(in) :: s
character*80, intent(out) :: w(n)
character*80   a,aa

aa=adjustl(s)
do i=1,n
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
enddo
return
end subroutine


!***************************************************
!* split a string s iand the return the x'ths word *
!***************************************************
subroutine charXsplit(s,wx,x)
implicit none
integer i,k,x
character*80, intent(in) :: s
character*80, intent(out) ::wx
character*80   w(20)
character*80   a,aa

aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
wx=w(x)
return
end subroutine



subroutine help(h)
implicit none
logical h

if(.not.h) return

print*,'                                  '
print*,' *    ----   HELP    ----    *    '
print*,'                                   '
print*,' gcp <coordinates> -level <method>   '
print*,'                                   '
print*,' example:                                 '
print*,' -> gcp water.xyz -func dft/631gd     '
print*,'                                  '
print*,'                                  '
print*,' Valid input coordinates are      '
print*,'   - TURBOMOLE files (in bohr)    '
print*,'   - XMOL files (in angstrom)     '
print*,'                                  '
print*,' available parametrisation for <method>        '
print*,'        HF/MINIS    DFT/MINIS     '
print*,'        HF/SV       DFT/SV '
print*,'        HF/SVP      DFT/SVP'
print*,'        HF/631Gd    DFT/631Gd'
print*,'        HF/TZ       DFT/TZ'
print*,'                                  '
print*,' note: <method> does not need to be          '
print*,'       capitelized                           '
print*,'                                  '
print*,' command line options:            '
print*,'    -h                       this print out   '
print*,'    -level <string>    specify <method>   '
print*,'      if <method> = file or left empty a parameter   '
print*,'      file from ~/.gcppar.$HOSTNAME will be read in    '
print*,'    -l                       same as -level   '
print*,'    -grad                    request gradient '
print*,'    -noprint                 surpress output   '
print*,'    -parfile                 print <gcp.param> '
print*,'    -local                   use  .gcppar in work-dir '
print*,'    -hess                    request hessian  '
print*,'    -test                    stop after parameter setup  '
print*,'    -verb                    verbose output: print gradient to stdout  '
print*,'                             instead of gcp_gradient     '
print*,'                                  '
print*,'                                  '
print*,' *    ----   ****    ----    *    '

call done('help',6)
end subroutine help

subroutine setzet(eta,za,zb)    
!*************************
!* set slater exponents  *
!*************************
implicit none
integer i                 
real(8) za(36),zb(36)
real(8) ZS(36),ZP(36),ZD(36),eta
data ZS /1.2000,1.6469,0.6534,1.0365,1.3990,1.7210,2.0348,2.2399,2.5644,2.8812,&
0.8675,1.1935,1.5143,1.7580,1.9860,2.1362,2.3617,2.5796,0.9362,1.2112,&
1.2870,1.3416,1.3570,1.3804,1.4761,1.5465,1.5650,1.5532,1.5781,1.7778,&
2.0675,2.2702,2.4546,2.5680,2.7523,2.9299/
data ZP /0.0000,0.0000,0.5305,0.8994,1.2685,1.6105,1.9398,2.0477,2.4022,2.7421,&
0.6148,0.8809,1.1660,1.4337,1.6755,1.7721,2.0176,2.2501,0.6914,0.9329,&
0.9828,1.0104,0.9947,0.9784,1.0641,1.1114,1.1001,1.0594,1.0527,1.2448,&
1.5073,1.7680,1.9819,2.0548,2.2652,2.4617/
data ZD /0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,&
0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,&
2.4341,2.6439,2.7809,2.9775,3.2208,3.4537,3.6023,3.7017,3.8962,2.0477,&
2.4022,2.7421,0.6148,0.8809,1.1660,1.4337/

  do i=1,36
select case (i)
  case(:2)
    za(i)=ZS(i)
  case(3:20,31:)
    za(i)=( ZS(i)+ZP(i) )/2d0
  case(21:30)
    za(i)=( ZS(i)+ZP(i)+ZD(i) )/3d0
end select 
  enddo


  za=za*eta
  zb=za

 return
end subroutine setzet

!CHARACTER*2 FUNCTION ESYM(I)
!CHARACTER*2 ELEMNT(94)
!DATA ELEMNT/'h ','he',&
!'li','be','b ','c ','n ','o ','f ','ne',&
!'na','mg','al','si','p ','s ','cl','ar',&
!'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
!'zn','ga','ge','as','se','br','kr',&
!'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',&
!'cd','in','sn','sb','te','i ','xe',&
!'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
!'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
!'au','hg','tl','pb','bi','po','at','rn',&
!'fr','ra','ac','th','pa','u ','np','pu'/ 
!ESYM=ELEMNT(I)
!RETURN
!END
