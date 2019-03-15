module systemVariables
  implicit none

  type system
      real(kind = 8), allocatable :: positions(:,:)
      real(kind = 8), allocatable :: velocities(:,:)
      real(kind = 8), allocatable :: masses(:)
      real(kind = 8), allocatable :: forces(:,:)
      real(kind = 8), allocatable :: cell(:)
      real(kind = 8)              :: potential, kinetic

  end type

  type forceAndPotential
    real(kind = 8), allocatable :: forces(:,:)
    real(kind = 8) :: totalPotential
  end type

end module systemVariables
