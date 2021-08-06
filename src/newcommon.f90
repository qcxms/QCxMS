module newcommon
  implicit none

  integer      :: calls

  character(len=80) :: path = '/usr/local/bin/'
  character(len=80) :: xtbpath = '$HOME/.XTBPARAM/'
  character(len=80) :: xtbhome
  character(len=20) :: solvent

  logical  :: TempRun
  logical  :: Verbose

end module newcommon
