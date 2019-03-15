module integrators
  use systemVariables
  use Potentials
  use BookKeeping
  implicit none

private

  public :: kineticEnergy
  public :: pairList
  public :: updatePairList
  public :: interactionsWithBookKeeping
  public :: interactions
  public :: velocityVerlet

contains

  pure function kineticEnergy(velocities)
    real(kind = 8), intent(in) :: velocities(:,:)
    real(kind = 8) :: kineticEnergy
    integer :: i, j, nAtoms, nDim
    kineticEnergy = 0.0d0
    nAtoms = size(velocities, 1)
    nDim = size(velocities, 2)
    do j = 1, nDim
      do i = 1, nAtoms
        kineticEnergy = kineticEnergy + velocities(i,j)**2
      enddo
    enddo
    kineticEnergy = 0.5*kineticEnergy
  end function kineticEnergy

!add pairlist as a variable
  pure function interactionsWithBookKeeping(pList, velocities, positions, rcPlusMargin, cell, margin, dt) result(fp)
    !{-# Calculate interactions among each particle -#}
    !Args:
    !  positions(No of atoms, no of dimensinos):
    !  rcPlusMargin: cut-off radius
    !  cell(no of dimensions): periodic boundary box size
    !  margin:
    !  dt: time step
    !Returns:
    !  fp: a structure containig forces(no of atoms, no of dimensions)
    !                            and total potential of a sysem
    type(pair)    , intent(in) :: pList(:)
    real(kind = 8), intent(in)  :: velocities(:,:), positions(:,:), rcPlusMargin, cell(:), margin, dt
    real(kind = 8), allocatable :: forces(:,:), dr(:)
    real(kind = 8) :: potentialCounter, r, ecut
    integer :: i, j, ipair, nAtoms, nDimensions, nPairs
    type(forceAndPotential) :: fp
    type(pair), allocatable :: updatedPList(:)

    nAtoms      = size(positions, 1)
    nDimensions = size(positions, 2)
    allocate(forces(nAtoms, nDimensions))

    forces    = 0.0d0
    potentialCounter = 0.0d0
    ecut = lennardJones(rcPlusMargin)

    !pList = pairList(positions, rcPlusMargin, cell)
    !updatedPList = updatePairList(pList, velocities, positions, rcPlusMargin, cell, margin, dt)
    nPairs = size(pList)

    do ipair = 1, nPairs
      i = updatedPList(ipair)%i
      j = updatedPList(ipair)%j

      !---difference of positional vectors
      dr = positions(i, :) - positions(j, :) ! vector like v1 - v2
      !---periodic boundary correction, which realises minimum image convention???
      dr(:) = dr(:) - cell(:) * (nint(dr(:) / cell(:)))
      !---Norm of dr, e.g., sqrt( dr(i, 1)**2 + dr(i, 2)**2 + dr(i, 3)**2 )
      r = norm2(dr) ! This should be squred r to reduce computational cost like that r2 = dot_product(dr,dr)
      if (r < rcPlusMargin) then
        potentialCounter = potentialCounter + lennardJones(r) - ecut !???Why ecut???
                                              !^ should be potential('type', r) which specifies a type of potentials.
        forces(i,:)  = forces(i,:) + ljForce(r) * dr(:)
        forces(j,:)  = forces(j,:) - ljForce(r) * dr(:)
                                                  !^ I wanna include it in ljForce
      endif

    enddo
    !---fp, which is a structure, is the one returned
    fp%totalPotential = potentialCounter
    fp%forces         = forces

  end function


!--Calculation of interactions without book keeping
  function interactions(positions, rc, cell) result(fp)
    !{-# Calculate interactions among each particle -#}
    !Args:
    !  positions(No of atoms, no of dimensinos):
    !  rc: cut-off radius
    !  cell(no of dimensions): periodic boundary box size
    !Returns:
    !  fp: a structure containig forces(no of atoms, no of dimensions)
    !                            and total potential of a sysem
    real(kind = 8), intent(in)  :: positions(:,:), rc, cell(:)
    real(kind = 8), allocatable :: forces(:,:), dr(:), cellHalf(:)
    real(kind = 8) :: potentialCounter, r, ecut !, r2, rc2
    integer :: i, j, nAtoms, nDimensions
    type(forceAndPotential) :: fp

    nAtoms      = size(positions, 1)
    nDimensions = size(positions, 2)
    allocate(forces(nAtoms, nDimensions))

    forces    = 0.0d0
    potentialCounter = 0.0d0
    cellHalf = cell*0.5d0
    ecut = lennardJones(rc)
!    rc2 = rc**2
    do i = 1, nAtoms - 1
      do j = i + 1, nAtoms
        !---difference of positional vectors
        dr = positions(i, :) - positions(j, :) ! vector like v1 - v2

        !---periodic boundary correction, which realises minimum image convention???
        dr(:) = dr(:) - cell(:) * (nint(dr(:) / cell(:)))

        !---Norm of dr, e.g., sqrt( dr(i, 1)**2 + dr(i, 2)**2 + dr(i, 3)**2 )
        r = norm2(dr) ! This should be squred r to reduce computational cost like that r2 = dot_product(dr,dr)

        if (r < rc) then
          potentialCounter = potentialCounter + lennardJones(r) - ecut !???Why ecut???
                                                !^ should be potential('type', r) which specifies a type of potentials.
          forces(i,:)  = forces(i,:) + ljForce(r) * dr(:)
          forces(j,:)  = forces(j,:) - ljForce(r) * dr(:)
                                                    !^ I wanna include it in ljForce

        endif

      enddo
    enddo
    !---fp, which is a structure, is the one returned
    fp%totalPotential = potentialCounter
    fp%forces = forces

  end function interactions

!--I want to a function "integrator" from which I can access any integrator.
  function velocityVerlet(sys, rc, dt) result(updatedSys)
    !{-# impliment velocity veret algorithm #-}
    !Args:
    !  sys: system structure defined in "variable module".
    !       For more detial, please refer to variables.f08
    !  rc: cut-off radius
    !  cell(no of dimensions): periodic boundary cell
    !Returns:
    !  updateSys: this is the same structure as sys, but is time-evolved sys.
    type(system)  , intent(in) :: sys
    real(kind = 8), intent(in) :: rc, dt
    type(system)               :: updatedSys
    type(forceAndPotential)    :: updatedFP

    !---sys's variables are inherited with updatedSys
    updatedSys = sys

    !---update positions: r(t+dt) = r(t) + (1/2)v(t)dt
    updatedSys%positions = sys%positions + sys%velocities * dt + 0.5d0 * sys%forces * dt**2
    !---calculate forces and potential
    !updatedFP = interactionsWithBookKeeping(updatedSys%positions, rc, cell)
    updatedFP = interactions(updatedSys%positions, rc, updatedSys%cell)
    !---update velocities: v(t+dt) = v(t) + (1/2) (F(t) + F(t+dt)) dt
    updatedSys%velocities = sys%velocities + (sys%forces + updatedFP%forces) * dt * 0.5d0

    !--updatedSys is the one returned.
    updatedSys%forces    = updatedFP%forces
    updatedSys%potential = updatedFP%totalPotential
    updatedSys%kinetic   = kineticEnergy(updatedSys%velocities)

  end function velocityVerlet

!  function velocityVerletWithBookKeeping(sys, pList, rcPlusMargin, dt, cell, margin) result(updatedSys)
!    !{-# impliment velocity veret algorithm #-}
!    !Args:
!    !  sys: system structure defined in "variable module".
!    !       For more detial, please refer to variables.f08
!    !  rc: cut-off radius
!    !  cell(no of dimensions): periodic boundary cell
!    !Returns:
!    !  updateSys: this is the same structure as sys, but is time-evolved sys.
!    type(system)  , intent(in) :: sys
!    real(kind = 8), intent(in) :: rcPlusMargin, dt, cell(:)
!    type(pair)    , intent(in) :: pList(:)
!    type(system)               :: updatedSys
!    type(forceAndPotential)    :: updatedFP
!
!    !---update positions: r(t+dt) = r(t) + (1/2)v(t)dt
!    updatedSys%positions = sys%positions + sys%velocities * dt + 0.5d0 * sys%forces * dt**2
!    !---calculate forces and potential
!    updatedFP = interactionsWithBookKeeping(updatedSys%positions, rcPlusMargin, cell)
!    !---update velocities: v(t+dt) = v(t) + (1/2) (F(t) + F(t+dt)) dt
!    updatedSys%velocities = sys%velocities + (sys%forces + updatedFP%forces) * dt * 0.5d0
!
!    !--updatedSys is the one returned.
!    updatedSys%forces    = updatedFP%forces
!    updatedSys%potential = updatedFP%totalPotential
!    updatedSys%kinetic   = kineticEnergy(updatedSys%velocities)
!    updatedSys%masses    = sys%masses
!
!  end function velocityVerletWithBookKeeping



!--I want to make a more small parts of integtors,
!  pure function velocityVerlet(positions, velocities, forces, dt) result(rvf)

!  end function

  pure function leapFrog(sys, rc, cell, dt) result(updatedSys)
    real(kind = 8), intent(in) :: rc, dt, cell(:)
    type(system), intent(in) :: sys
    type(system) :: updatedSys
  end function leapFrog

!--this is positional verlet algorithm, do not use this.
  subroutine verlet(sysCurrent, sysPrevious, forces, dt, kineticE)
    real(kind = 8), intent(in)    :: dt
    real(kind = 8), intent(in)    :: forces(:,:)
    type(system)  , intent(inout) :: sysCurrent, sysPrevious
    real(kind = 8), intent(out)   :: kineticE
    integer :: i, j, nAtoms, nDimensions
    type(system) :: sysTmp

    nAtoms      = size(forces, 1)
    nDimensions = size(forces, 2)

    sysTmp = sysCurrent !@@@@@@@@@@@@@@
    kineticE = 0.0d0
    do i = 1, nAtoms
      !r(t+dt)= 2*r(t) - r(t-dt) + f*dt**2
      sysCurrent%positions(i,:)  = 2.0d0 * sysCurrent%positions(i,:) - sysPrevious%positions(i,:) +  forces(i,:)* dt**2
      !v(t+dt)= ( r(t+dt) - r(t-dt) ) / (2 * dt)
      sysCurrent%velocities(i,:) = (sysCurrent%positions(i,:) - sysPrevious%positions(i,:)) / (2.0d0 * dt)
      do j = 1, nDimensions
        kineticE = kineticE + sysCurrent%velocities(i,j)**2
      enddo
    enddo
    sysPrevious = sysTmp !@@@@@@@@@@@@
    kineticE = 0.5 * kineticE
  end subroutine


end module integrators
