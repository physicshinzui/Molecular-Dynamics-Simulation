module BookKeeping
  implicit none

  type pair
    integer :: i, j
  end type

contains
  pure function pairList(positions, rcPlusMargin, cell)
    real(kind = 8), intent(in) :: positions(:,:), rcPlusMargin, cell(:)
    real(kind = 8), allocatable :: dr(:)
    integer :: i, j, k, nAtoms, nDimensions
    real(kind = 8) :: r
    type(pair), allocatable :: pairList(:)

    nAtoms      = size(positions, 1)
    nDimensions = size(positions, 2)
    allocate(pairList( (nAtoms * (nAtoms - 1)) / 2) )

    k = 1
    do i = 1, nAtoms - 1
      do j = i + 1, nAtoms
        !---difference of positional vectors
        dr = positions(i, :) - positions(j, :) ! vector like v1 - v2
        !---periodic boundary correction, which realises minimum image convention???
        dr(:) = dr(:) - cell(:) * (nint(dr(:) / cell(:)))
        !---Norm of dr, e.g., sqrt( dr(i, 1)**2 + dr(i, 2)**2 + dr(i, 3)**2 )
        r = norm2(dr) ! This should be squred r to reduce computational cost like that r2 = dot_product(dr,dr)
        if (r < rcPlusMargin) then
          pairList(k)%i = i
          pairList(k)%j = j
          k = k + 1
        endif

      enddo
    enddo

  end function pairList

  function updatePairList(pList, velocities, positions, rcPlusMargin, cell, margin, dt) result(updatedPairList)
    type(pair)    , intent(in) :: pList(:)
    real(kind = 8), intent(in) :: velocities(:,:), positions(:,:), rcPlusMargin, cell(:), margin, dt
    type(pair), allocatable :: updatedPairList(:)
    integer :: i, j
    real(kind = 8) :: vmax, vmax2, v2, marginLength

    vmax2 = 0.0d0
    do i = 1, size(velocities, 1)
      v2 = dot_product(velocities(i,:), velocities(i,:))
      if (v2 > vmax2) vmax2 = v2
    enddo

    vmax = sqrt(vmax2)
    marginLength = marginLength - vmax*2.0d0*dt

    if (marginLength < 0.0d0) then
      marginLength = margin
      updatedPairList = pairList(positions, rcPlusMargin, cell)

    else
      updatedPairList = pList

    endif

  end function

end module
