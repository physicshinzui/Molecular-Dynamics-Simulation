program main
    use variables
    use parameters !This in future is read from an input file.
    use initialisation
    use analysis
    use IO
    use lattice
    use PBC
    use thermostats
    use integrators
    implicit none
    real(kind = 8) ::  t1, t2
    integer        :: i, j, nAtoms, nDimensions
    type(system)   :: sys
    type(forceAndPotential) :: fp
    real(kind = 8), allocatable :: centreOfg(:) !system's cog

    !---latice system specific parameters
    real(kind = 8), parameter :: latticeWidth = 1.5d0 !(4.0d0 / 0.817657)**(1/3) !1.0d0
    integer       , parameter :: nLattices    = 3

    call cpu_time(t1)

    !---FFC generation for argon. Note: sys%positions are allocated if an array is substituted into it.
    sys%positions = latticeGenetator(nLattices, latticeWidth, "FCC")
    nAtoms        = size(sys%positions, 1)
    nDimensions   = size(sys%positions, 2)
    allocate(sys%masses(nAtoms))
    allocate(cell(nDimensions))
    sys%masses(:) = 1.0d0
    cell(:) = 0.0d0

    !---for FCC
    cell(:) = latticeWidth*nLattices  !@@@@@@@@
    !---for normal lattices
!    cell(:) = latticeWidth*(nLattices - 1) !@@@@@@@@

    centreOfg = centreOfGeometry(sys%positions)
    print*, centreOfg
    call writeLatticeCoordsAsPDB(sys%positions)

    write(*, "(a, i5, /a, 3f10.3)")"No. of atoms : ", nAtoms, "Cell: ", cell

    !---initialise velocities
    sys%velocities = scaleVelocities( velocityGenerator(nDimensions, nAtoms), sys%masses, initialT)
    call writeVelocities(sys%velocities)

    !---initial forces and potential
    fp = interactions(sys%positions, rc, cell)
    sys%forces    = fp%forces
    sys%potential = fp%totalPotential

    !---------------------------------
    !(5) MD loop
    write(*, "(a)") "#---MD is started---"
    do i = 1, nSteps
      !---sys is updated iteratively through md loop.
      !updatedPairList = updatePairList(pList, velocities, positions, rc, cell, margin, dt)
      sys = velocityVerlet(sys, rc, dt, cell)

      sys%positions  = periodic(sys%positions, cell, centreOfg)
      !sys = thermostat(sys, initialT, presetT, heatingSteps, i)
      !sys%velocities = scaleVelocities(sys%velocities, presetT) !simple v-scaling

      if (any(sys%positions > 10**3 )) then
        print*, sys%positions(:,1)
        stop "STOP: positions are far away."
      endif

      if (mod(i,oInterval) == 0) then

        write(*, "(a, i10)") "Loop No.", i
        call writeLatticeCoordsAsPDB(sys%positions)
        write(*,"(a, 3f10.3)") "    Sum of Forces(r) = ", sum(sys%forces)
        print*,""
        !call energyIO(sys%potential, sys%kinetic, Tt(sys%velocities, sys%masses))

      endif
      call energyIO(sys%potential, sys%kinetic, Tt(sys%velocities, sys%masses))

    enddo

    call cpu_time(t2)
    write(*,"(a,f12.3,a)") "Terminated: ", t2 - t1, " sec"

end program
