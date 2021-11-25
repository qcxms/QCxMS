module common1
  use xtb_mctc_accuracy, only: wp 
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!
  integer  :: shell
  integer  :: bas
  integer  :: func
  integer  :: method
  integer  :: gfnver
  integer  :: prog
  integer  :: qcmem
  integer  :: ihamilt
  integer  :: orca_version
  integer  :: grid_orca
  !integer  ::  grid_tmol
  integer  :: nproc_orca

  real(wp) ::  a1,a2,s8
  real(wp) ::  ax
  real(wp) ::  ieetemp
  real(wp) ::  cab(100,100)
  real(wp) ::  alpc

  logical  ::  gcp
  logical  ::  hhmod
  logical  ::  tbcorr
  logical  ::  tun
  logical  ::  noconv
  logical  ::  slowconv
  logical  ::  XTBMO
  logical  ::  No_eTemp
  !!!!!!!!!!!!!!!!!!!!!!!!!

end module common1
