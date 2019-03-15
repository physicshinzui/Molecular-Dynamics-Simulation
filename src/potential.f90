module potential 
    implicit none

contains

    pure function lennardJones(r) result(potential)
      real(kind = 8), intent(in) :: r
      real(kind = 8) :: potential
      integer :: nAtoms, nPairs
      potential = 4 * (r**(-12) - r**(-6))
    end function

end module
