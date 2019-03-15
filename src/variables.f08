module Variables
  implicit none
  
  type system
      real(kind = 8), allocatable :: positions(:,:)     !px(:), py(:), pz(:)
      real(kind = 8), allocatable :: velocities(:,:)    !vx(:), vy(:), vz(:)
      real(kind = 8), allocatable :: masses(:)
      real(kind = 8) :: potential, kinetic
      real(kind = 8), allocatable :: forces(:,:)
  end type

  type forceAndPotential
    real(kind = 8), allocatable :: forces(:,:)
    real(kind = 8) :: totalPotential
  end type

end module Variables
