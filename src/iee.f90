module qcxms_iee
  use xtb_mctc_accuracy, only: wp
  implicit none

  contains

! determine automatically the two parameters in the P(E)
! function for given energy range exc, the number of bonds
! nbnd and the input parameter ieeatm (av energy/per bond in eV,
! normally 0.5 eV/bond)

  subroutine getieeab(iee_a,iee_b,ieeel,ityp,exc,nbnd,ieeatm)
  
     integer  :: ityp,nbnd,k
  
     real(wp) :: iee_a,iee_b,ieeel,pmax,ieemax,exc,E_avg,st,ieeatm

  
     st = 0.005_wp
     iee_a = 0.0_wp
     iee_b = 0.0_wp
     k = 0
  
     do
       k = k + 1

       iee_a = iee_a + st
       iee_a = min(iee_a,0.3)
       iee_b = iee_b + st * 7

       call getmaxiee(iee_a,iee_b,ieeel,ityp,exc,ieemax,pmax,E_avg)

       if (k > 10000) stop 'internal error inside getieeab'
       if (E_avg / nbnd >= ieeatm) exit 


     enddo
  
  end subroutine getieeab
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! get maximum value (pmax) at which energy (ieemax) and
  ! average energy E_avg for given parameters in the P(E)
  ! distribution function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in fact, only getting pmax here is the important part
  subroutine getmaxiee(iee_a,iee_b,ieeel,ityp,exc,ieemax,pmax,E_avg)
  
     integer  :: ityp
     real(wp) :: iee_a,iee_b,ieeel,val,x,pmax,ieemax,exc,E_avg,m
  
     x = 0.001_wp
     pmax = -1.0_wp
     ieemax = 0.0_wp
     E_avg = 0.0_wp
     m = 0.0_wp
  
     do
        if ( ityp == 0 ) call gauss0(iee_a,iee_b,ieeel,x,val)
        if ( ityp == 1 ) call poiss0(iee_a,iee_b,ieeel,x,val)

        if ( val > pmax ) then
           pmax   = val
           ieemax = x
        endif

        x = x + 0.01_wp
        E_avg = E_avg + val * x
        m = m + val

        if ( x >= exc ) exit
     enddo
  
     E_avg = E_avg / m
  
  end subroutine getmaxiee
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Gaussian energy distribution, ieeel= # val el, iee_a and iee_b are parameters
  subroutine gauss0(iee_a,iee_b,ieeel,x,p)
  
     real(wp) :: iee_a,iee_b,ieeel,x,p
  
     p = exp( -iee_a * (x - ieeel * iee_b)**2 / ieeel )
  
  end subroutine gauss0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gauss1(height,posi,breite,ergebnis,ntraj, nuc)

    integer :: ntraj
    integer :: i, nuc
  
     real(wp) :: height, posi, breite, wert, x
     
     real(wp) :: ergebnis(ntraj)

     !x = posi / ntraj
     x = posi / (ntraj/2)
     ! write(*,*) x

     do i = 1, ntraj
     !do i = 1, 10

      wert = x * i
      !write(*,*) i
      !write(*,*) 'WERT', i, wert
      !write(*,*) 'position',posi

       !ergebnis(i) = height *  exp(-1* (wert - posi)**2 / (2* breite) )
       !ergebnis(i) = exp(-1* bla**2 )
       ergebnis(i) = (wert-posi) / 8 

      !write(*,*) 'ergebnis', ergebnis(i)
    enddo
  
  end subroutine gauss1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Poisson energy distribution, ieeel= # val el, iee_a and iee_b are parameters
  subroutine poiss0(iee_a,iee_b,ieeel,x,t18)
  
     real(wp) :: iee_a,iee_b
     real(wp) :: k,z,ieeel,x,t18
     real(wp) :: t2,t8,t14,t17
     
     z = iee_b
     k = 1.0_wp / iee_a
     t2 = k / ieeel
     t8 = dlog(z / k * ieeel / x)
     t14 = dexp(t2 * x * (0.1D1 + t8) - 1.D0 * z)
     t17 = (t2 * x + 0.1D1)**(-0.5_wp)
     t18 = t14 * t17

  end subroutine poiss0

end module qcxms_iee
