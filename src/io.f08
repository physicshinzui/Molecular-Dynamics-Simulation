module IO
  implicit none
  integer, private :: i !, nAtoms

contains

  subroutine energyIO(potential, kinetic, Tt, maxForce)
    real(kind = 8), intent(in) :: potential, kinetic, Tt, maxForce
    open(unit=123, file="energy.dat", action="write")
    write(123,'(5f12.3)') potential, kinetic, potential + kinetic, Tt, maxForce
    !close(unit=11)
  end subroutine energyIO

  subroutine writeLatticeCoordsAsPDB(coordinates)
    real(kind=8), intent(in) :: coordinates(:,:)
    integer :: nAtoms, nDimensions

    nAtoms = size(coordinates, 1)
    open(unit=11, file = "lattice.pdb")
    nDimensions = size(coordinates, 2)

    write(11, '(a)') "MODEL"
    do i  = 1, nAtoms
      if (nDimensions == 1) then
        write(11, '(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,f8.3, 2f8.3)') &
                  "ATOM  ", i, "LT  ", "LTC", " ", i, coordinates(i, :), 0.0d0, 0.0d0

      elseif (nDimensions == 2) then
        write(11, '(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,2f8.3, f8.3)') &
                    "ATOM  ", i, "LT  ", "LTC", " ", i, coordinates(i, :), 0.0d0

      elseif (nDimensions == 3) then
      write(11, '(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)') &
                  "ATOM  ", i, "LT  ", "LTC", " ", i, coordinates(i, :)
      endif

    enddo
    write(11, '(a)') "ENDMDL"
  !  close(11)
  end subroutine

  subroutine writeVelocities(velocities)
    real(kind = 8), intent(in) :: velocities(:,:)
    open(12, file = "initialVelocities.out")
    do i = 1, size(velocities, 1)
      write(12, '(3f10.3)') velocities(i,:)
    enddo
    write(12,*) "#sum of velocities (is it zero??)", sum(velocities(:, 1))
    close(12)
  end subroutine


end module
