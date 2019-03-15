module parameters
!--USE ONLY IN THE MAIN. DO NOT REFER TO THEM IN THE OTHER.

!private
  !---parameters which will be inputed from a file
  integer       , parameter :: nSteps  = 5000
  real(kind = 8), parameter :: initialT = 1.0d0, presetT = 1.5d0
  real(kind = 8), parameter :: dt = 0.001d0
  real(kind = 8), parameter :: rc = 4.0d0 !, rc2 = rc**2
  integer       , parameter :: heatingSteps = 1000
  integer       , parameter :: oInterval = 100

  real(kind = 8), parameter :: margin = 0.5d0
  real(kind = 8), parameter :: rcPlusMargin = rc+margin !or ideally rcPlusMargin^2

  real(kind = 8), allocatable :: cell(:)

  !public :: readParams

contains

!  subroutine readParams
!    integer :: i
!
!    do
!
!    enddo
!  end subroutine

end module
