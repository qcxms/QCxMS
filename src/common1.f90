module common1
  use xtb_mctc_accuracy, only: wp 
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!
  integer  ::  shell
  integer  ::  bas
  integer  ::  func
  integer  ::  method
  integer  ::  gfnver
  integer  ::  prog
  integer  ::  qcmem
  integer  ::  ihamilt
  integer  ::  orca_version
  integer  ::  grid_orca
  !integer  ::  grid_tmol

  real(wp) ::  a1,a2,s8
  real(wp) ::  ax
  real(wp) ::  ieetemp
  real(wp) ::  cab(100,100)
  real(wp) ::  alpc
  real(wp) ::  ipshift,eashift


  logical  ::  gcp
  logical  ::  hhmod
  logical  ::  tbcorr
  logical  ::  tun
  logical  ::  noconv
  logical  ::  slowconv
  logical  ::  XTBMO
  !!!!!!!!!!!!!!!!!!!!!!!!!

end module common1
