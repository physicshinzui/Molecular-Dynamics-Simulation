module analysis
  implicit none

contains

  pure function centreOfGeometry(positions) result(cog)
    real(kind = 8), intent(in)  :: positions(:,:)
    real(kind = 8), allocatable :: cog(:)
    integer :: i, j, nAtoms, nDimensions
    nAtoms = size(positions, 1)
    nDimensions = size(positions, 2)
    allocate(cog(nDimensions))

    cog(:) = 0.0d0
    do j = 1, nDimensions
      do i = 1, nAtoms
        cog(j) = cog(j) + positions(i,j)
      enddo
    enddo
    cog = cog / nAtoms

  end function centreOfGeometry

  pure function distances(positions)
    real(kind = 8), intent(in)  :: positions(:,:)
    real(kind = 8), allocatable :: distances(:)
    integer :: i, j, ipair, nAtoms, nPairs
    real(kind = 8) :: dr(size(positions, 1))

    nAtoms = size(positions, 2)
    nPairs = (nAtoms * (nAtoms - 1)) / 2.0d0
    allocate(distances(nPairs))
    distances(:) = 0.0d0

    ipair = 1
    do i = 1, nAtoms - 1
      do j = i + 1, nAtoms
        dr(:) = positions(:,i) - positions(:,j)
        distances(ipair) = sqrt(sum(dr**2)) !**0.5
        ipair = ipair + 1
      enddo
    enddo

  end function

end module

!program main
!  use analysis
!  implicit none
!  integer :: i, array(2,2)
!  real(8) :: pos(2,2) = transpose(reshape( (/ 1, 0, 0, 1, 0, 1 /), shape(array))), dists(size(pos, 2) * (size(pos, 2) - 1))
!
!  print*, pos(:,:)
!  print*, distances(pos)
!
!  !!!!FIXME
!
!end program
