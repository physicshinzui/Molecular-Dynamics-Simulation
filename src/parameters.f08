module parameters
!--USE ONLY IN THE MAIN. DO NOT REFER TO THEM IN THE OTHER.

!---parameters which will be inputed from a file
  integer       , parameter :: nSteps    = 5000
  real(kind = 8), parameter :: initialT  = 1.0d0, presetT = 1.5d0
  real(kind = 8), parameter :: dt        = 0.001d0
  real(kind = 8), parameter :: rc        = 2.0d0 !, rc2 = rc**2
  integer       , parameter :: oInterval = 100
  logical       , parameter :: isPBC     = .True. !Curently not used

  !--For thermostat
  integer       , parameter :: heatingSteps = 1000
  logical       , parameter :: isThermostat =  .False. !.True.

  !--For book keeping
  real(kind = 8), parameter :: margin = 0.5d0
  real(kind = 8), parameter :: rcPlusMargin = rc+margin !or ideally rcPlusMargin^2

  real(kind = 8), allocatable :: cell(:)

  !---latice system specific parameters
  real(kind = 8), parameter :: latticeWidth = 1.5d0 !(4.0d0 / 0.817657)**(1/3) !1.0d0
  integer       , parameter :: nLattices    = 2

end module
