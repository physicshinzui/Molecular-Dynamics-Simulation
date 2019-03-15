function forces(positions, rc, cell) !This function will be moved to a different module.
  real(kind = 8), intent(in)  :: positions(:,:), rc, cell(:)
  real(kind = 8), allocatable :: forces(:,:)
  real(kind = 8), allocatable :: dr(:), drPBCcorrected(:)
  real(kind = 8) :: r, ecut, totPotential
  integer :: i, j, nAtoms, nDimensions, nPairs

  nAtoms      = size(positions, 1)
  nDimensions = size(positions, 2)
  allocate(forces(nAtoms, nDimensions))
  forces      = 0.0d0
  totPotential = 0.0d0
  ecut = lennardJones(rc)
  do i = 1, nAtoms - 1
    do j = i + 1, nAtoms
      dr = positions(i, :) - positions(j, :) ! vector like v1 - v2
      drPBCcorrected = dr(:) - cell(:) * (nint(dr(:) / cell(:)))

!          print*, "dbg:", size(positions, 2) == size(dr)
!          print*, "size of dr:", size(dr)

!          print*, "dr              = ", dr
!          print*, "cell            = ", cell
!          print*, "dr/cell         = ", dr / cell
!          print*, "cell*nint(dr/cell)", cell * (nint(dr(:) / cell(:)))

      r = norm2(dr)
      if (r < rc) then
!            print*, "norm2(dr) is less than rc. dr=", r, "rc=",rc
        totPotential = totPotential + lennardJones(r)
        forces(i,:) = forces(i,:) + ljForce(r) * dr(:) !(dE/dr * dx, dE/dr * dy, dE/dr * dz)
        forces(j,:) = forces(j,:) - ljForce(r) * dr(:)
!???            en=en+4*r6i* (r6i-l) -ecut
      endif

    enddo
  enddo
  print*, "total potential = ", totPotential
end function
