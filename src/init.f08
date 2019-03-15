module initialisation
  implicit none
contains

  !---Box mular: http://mathworld.wolfram.com/Box-MullerTransformation.html
  function velocityGenerator(nDimensions, nAtoms) result(velocities)
    !This cannnot be a pure function because random number is generated here.
    integer         , intent(in) :: nAtoms, nDimensions
    real(kind = 8) :: velocities(nAtoms, nDimensions)
    real(kind = 8) :: u0(nDimensions), u1(nDimensions)
    real(kind = 8), parameter :: pi8 = 4 * atan(1.0_8)
    integer :: i

    velocities(:,:) = 0.0d0

    do i = 1, nAtoms
        call random_number(u0)
        call random_number(u1)
        velocities(i,:) = sqrt(-2 * log(u0(:))) * cos(2 * pi8 * u1(:))
    enddo

    !# to cancell translation of momentum of a system
    do i = 1, nDimensions
      velocities(:, i) = velocities(:, i) - (sum(velocities(:, i)) / nAtoms)
    enddo
  end function

end module
