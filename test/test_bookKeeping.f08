module bookKeeping

  type pair
    integer :: i, j
  end type

contains

  function pairList(positions, rc, cell)
    real(kind = 8), intent(in) :: positions(:,:), rc, cell(:)
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
        if (r < rc) then
          pairList(k)%i = i
          pairList(k)%j = j
          k = k + 1
        endif

      enddo
    enddo

  end function pairList

end module

program main
  use bookKeeping
  implicit none
  integer :: i
  type(pair), allocatable :: ppairs(:)
  real(8) :: arr(1:3,1) = reshape( (/ 1.0d0, 2.0d0, 4.0d0 /), (/ 3, 1 /) ), cell(1) = 10.0d0
  ppairs = pairList(arr, 1.0d0, cell)

  do i  = 1, size(ppairs)
    !print*, ppairs(i)%i, ppairs(i)%j
  enddo
end program
