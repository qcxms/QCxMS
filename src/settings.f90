module get_settings
  use mctc_env_accuracy, only : wp
  implicit none

  private

  public  :: run_settings, collision_type

  type  :: run_settings


    real(wp) :: tadd  
    real(wp) :: eimp  

    real(wp),allocatable :: xyzr(:,:,:)
    real(wp),allocatable :: velor(:,:,:)
    !real(wp),allocatable :: velof(:)
    real(wp),allocatable :: eimpr(:)
    real(wp),allocatable :: taddr(:)
    real(wp),allocatable :: velofr(:,:)

  end type


  !type  :: charge

  !  integer  :: mchrg = 0
  !  integer  :: mchrg_prod

  !  real(wp) :: chrgcont  

  !end type

  type  :: collision_type

    integer :: set_coll
    integer :: max_coll
  end type


end module get_settings
