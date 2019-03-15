program main
  use PBC
  use systemVariables
  implicit none
  type(system) :: sys, sysUpdated, systmp
  integer :: i

  allocate(sys%positions(3,2))
  allocate(sys%velocities(3,2))
  allocate(sys%masses(3))
  allocate(sys%forces(3,2))
  allocate(sys%cell(2))

  sys%positions(1:3,1:2)  = reshape( (/ -1, 2, 3, 4, 5, 6 /), (/3,2/))
  sys%velocities(1:3,1:2) = reshape( (/ -1, 2, 3, 4, 5, 6 /), (/3,2/))
  !Note reshape generates
  ![ [1,2,3],
  !  [4,5,6] ]
  sys%masses(1:3)         = 1.0d0
  sys%forces(1:3,1:2)     = 1.0d0
  sys%potential           = 0.0d0
  sys%kinetic             = 0.0d0
  sys%cell(1:2)           = 1.0d0

  print*, "shape: ", shape(sys%positions)
  sysUpdated = periodic(sys)

  do i = 1, size(sysUpdated%positions,1)
    print*, sysUpdated%positions(i,:)!, sysUpdated%positions(i,:)
  enddo

end program
