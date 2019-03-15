module lattice
  implicit none

  private
  public :: latticeGenetator

contains

  pure function latticeGenetator(nLattices, width, latticeType) result(latticeCoords)
    integer       , intent(in)  :: nLattices
    real(kind = 8), intent(in)  :: width ! width of lattice
    character(len=2), intent(in):: latticeType
    real(kind = 8), allocatable :: latticeCoords(:,:)

    if (latticeType == "1D") then
      latticeCoords = lattice1D(nLattices, width)

    elseif (latticeType == "2D") then
      latticeCoords = lattice2D(nLattices, width)

    elseif (latticeType == "FC") then
      latticeCoords = FaceCenteredCubicLattice(nLattices, width)
    endif

  end function

  pure function FaceCenteredCubicLattice(nLattices, width) result(latticeCoords)
    integer       , intent(in)  :: nLattices
    real(kind = 8), intent(in)  :: width ! width of lattice
    real(kind = 8), allocatable :: latticeCoords(:,:)
    integer :: nParticles
    integer :: x, y, z, iParticle

    nParticles = 4 * nLattices**3 !???? I dont understand

    allocate(latticeCoords(nParticles, 1:3))
    latticeCoords(:,:) = 0.0

    iParticle = 1
    do x = 0, nLattices - 1 ! why -1? because
      do y = 0, nLattices - 1
        do z = 0, nLattices - 1
          latticeCoords(iParticle, 1:3) = (/ x*width, y*width, z*width /)
          iParticle = iParticle + 1
          latticeCoords(iParticle, 1:3)  = (/ (0.5d0+x)*width, (0.5d0+y)*width, z*width /)
          iParticle = iParticle + 1
          latticeCoords(iParticle, 1:3)  = (/ (0.5d0+x)*width, y*width, (0.5d0+z)*width /)
          iParticle = iParticle + 1
          latticeCoords(iParticle, 1:3)  = (/ x*width, (0.5d0+y)*width, (0.5d0+z)*width /)
          iParticle = iParticle + 1
        enddo
      enddo
    enddo
  end function

  pure function lattice2D(nLattices, width) result(latticeCoords)
    integer       , intent(in)  :: nLattices
    real(kind = 8), intent(in)  :: width ! width of lattice
    real(kind = 8), allocatable :: latticeCoords(:,:)
    integer :: x, y, iLattice, i

    allocate(latticeCoords(nLattices**2, 1:2))
    latticeCoords(:,:) = 0.0
    iLattice = 1
    do x = 0, nLattices - 1
      do y = 0, nLattices - 1
        latticeCoords(iLattice, 1:2) = (/ x * width, y * width/)
        iLattice = iLattice + 1
      enddo
    enddo

  end function

  pure function lattice1D(nLattices, width) result(latticeCoords)
    integer       , intent(in)  :: nLattices
    real(kind = 8), intent(in)  :: width ! width of lattice
    real(kind = 8), allocatable :: latticeCoords(:,:)
    integer :: x, iLattice

    allocate(latticeCoords(nLattices, 1))
    latticeCoords(:,:) = 0.0

    do iLattice = 1, nLattices
      latticeCoords(iLattice, 1) = iLattice * width
    enddo
  end function

end module
