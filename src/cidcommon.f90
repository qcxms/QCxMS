module cidcommon
  use xtb_mctc_accuracy, only: wp
  implicit none

  public :: cell_settings

  type :: cell_settings
    real(wp) :: TGas  
    real(wp) :: PGas    
    real(wp) :: lchamb
  end type


  type :: gas_parameters
    integer  :: IndAtom
    integer  :: Iatom
    real(wp) :: mIatom
    real(wp) :: rIatom
  end type


  type(cell_settings) :: cell
  type(gas_parameters) :: gas

end module cidcommon

