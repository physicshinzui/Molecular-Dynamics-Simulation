module Thermostats
  use Variables
  implicit none

contains

  !pure function thermostat(velocities, initialT, aimedT, heatingSteps, istep) result(modifiedVelocities)
   function thermostat(sys, initialT, aimedT, heatingSteps, istep) result(modifiedSys)
     !--this is the interface throught which we can access various thermostats
    type(system), intent(in)   :: sys
    real(kind = 8), intent(in) :: initialT, aimedT !velocities(:,:)
    integer, intent(in) :: heatingSteps, istep
    type(system) :: modifiedSys
    !real(kind = 8), allocatable :: modifiedVelocities(:,:)
    real(kind = 8) :: deltaT

    deltaT = (aimedT-initialT)/heatingSteps

    if (istep >= 1 .and. istep <= heatingSteps ) then
      modifiedSys%velocities = scaleVelocities(sys%velocities, sys%masses, initialT + deltaT*istep )

    elseif(istep > heatingSteps) then
      modifiedSys%velocities = scaleVelocities(sys%velocities, sys%masses, aimedT)

    else
      error stop "Thermostat Error: Negative/zero loop numbers are obtained, which cannot be handled."

    endif

    !---this is so stupid expressions. I do not want to do them.
    modifiedSys%positions = sys%positions
    modifiedSys%masses = sys%masses
    modifiedSys%forces = sys%forces
    modifiedSys%potential = sys%potential
    modifiedSys%kinetic = sys%kinetic


  end function

  pure function squreSum(velocities, weight)
    real(kind = 8), intent(in) :: velocities(:,:)
    real(kind = 8), intent(in) :: weight(:)
    real(kind = 8) :: squreSum
    integer :: i, j

    squreSum = 0
    do j = 1, size(velocities, 2)
      do i = 1, size(velocities, 1)
        squreSum = squreSum + weight(i) * velocities(i,j)**2 !+ velocities(i,2)**2 + velocities(i,3)**2)
      enddo
    enddo
  end function

  pure function Tt(velocities, masses)
    !--instantaneus temperature
    real(kind = 8), intent(in) :: velocities(:,:), masses(:)
    real(kind = 8) :: Tt
    integer :: nAtoms, nDimensions
    nAtoms = size(velocities, 1)
    nDimensions = size(velocities, 2)

    Tt = squreSum(velocities, masses) / (nDimensions*(nAtoms - 1))
  end function Tt

  pure function scaleVelocities(velocities, masses, presetT) result(scaledVelocity)
    !
    !Args:
    !    velocities(nAtoms, nDim):
    !    presetT: temperature set beforehand
    !    masses: masses of atoms
    !Returns:
    !    scaledVelocity(nAtoms, nDim): velocities scaled to fit the preset temperature

    real(kind = 8), intent(in)  :: velocities(:,:), masses(:), presetT
    !real(kind = 8), intent(in)  :: presetT
    real(kind = 8), allocatable :: scaledVelocity(:,:)
    !real(kind = 8), allocatable :: mass(:), ones(:)
    real(kind = 8) :: squaredVelo
    real(kind = 8) :: w, instantaneousTt, scaledTt
    integer :: i, nAtoms, nDimensions

    nAtoms = size(velocities, 1)
    nDimensions = size(velocities, 2)
    !allocate(mass(nAtoms), ones(nAtoms))
    !mass(:) = 1.0d0 !???Hard coded. this will be set to each atom in future.

!    instantaneousTt = squreSum(velocities, mass) / (nDimensions*(nAtoms - 1))
    instantaneousTt = Tt(velocities, masses)

      !--Negative tempreature is not allowed.
    if (presetT >= 0.0d0) then

      if (instantaneousTt /= 0.0d0) then
        w = sqrt(presetT / instantaneousTt)
      else
        w = 0.0d0
      endif

      scaledVelocity = w * velocities(:,:)

    else
       error stop "Preset T must be positive." !Implimented from F2015

    endif

  end function


end module Thermostats
