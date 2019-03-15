module PBC
  use systemVariables
  implicit none
contains

!---this is currently specialised only for FCC lattice
    function periodic(sys, centerSys) result(pbcCorrectedSys)
      type(system)  , intent(in) :: sys
      real(kind = 8), intent(in) :: centerSys(:)
      type(system) :: pbcCorrectedSys
      real(kind = 8), allocatable :: halfCell(:), centerOfcell(:)
      integer :: iatom, idim

      !---Make pbcCorrectedSys remember variables of sys which are not updated here.
      pbcCorrectedSys = sys

      pbcCorrectedSys%positions(:,:) = 0.0d0
      halfCell     = sys%cell * 0.5d0
      !centerOfcell =

      do idim = 1, size(sys%positions, 2)
        do iatom = 1, size(sys%positions, 1)

          if (sys%positions(iatom, idim) < centerSys(idim) - halfCell(idim) ) then
!          if (sys%positions(iatom, idim) < 0.0d0 ) then
            pbcCorrectedSys%positions(iatom, idim) = sys%positions(iatom, idim) + sys%cell(idim)
            !print*, "1"

!          elseif (sys%positions(iatom, idim) > sys%cell(idim)) then
          elseif (sys%positions(iatom, idim) > centerSys(idim) + halfCell(idim)  ) then
            pbcCorrectedSys%positions(iatom, idim) = sys%positions(iatom, idim) - sys%cell(idim)
            !print*, "2"

          else
            pbcCorrectedSys%positions(iatom, idim) = sys%positions(iatom, idim)
            !print*, "3"

          endif

!---I wanna generalise this function in future
!          if (coords(iatom, idim) > centerSys(idim) + halfCell(idim) ) then
!            pbcCorrectedSys%positions(iatom, idim) = sys%positions(iatom, idim) - sys%cell(idim)
!
!          elseif (coords(iatom, idim) < centerSys(idim) - halfCell(idim)) then
!            pbcCorrectedSys%positions(iatom, idim) = sys%positions(iatom, idim) + sys%cell(idim)
!
!          else
!            pbcCorrectedSys%positions(iatom, idim) = sys%positions(iatom, idim)

!          endif

        enddo
      enddo
    end function periodic

end module
