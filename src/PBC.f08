module PBC
  use systemVariables
  implicit none
contains

    pure function periodic(sys, centerSys) result(pbcCorrectedSys)
      !
      !Args:
      !  sys: defined in "systemVariables" module.
      !  centerSys: Center of cell.
      !Return:
      !  pbcCorrectedSys: positions are corrected if they are over the cell.
      !
      type(system)  , intent(in)  :: sys
      real(kind = 8), intent(in)  :: centerSys(:)
      type(system)                :: pbcCorrectedSys
      real(kind = 8), allocatable :: halfCell(:), centerOfcell(:)
      integer :: iatom, idim

      !---Make pbcCorrectedSys remember variables of sys which are not updated here.
      pbcCorrectedSys = sys

      pbcCorrectedSys%positions(:,:) = 0.0d0
      halfCell     = sys%cell * 0.5d0

      do idim = 1, size(sys%positions, 2)
        do iatom = 1, size(sys%positions, 1)

          if (sys%positions(iatom, idim) < centerSys(idim) - halfCell(idim) ) then
            pbcCorrectedSys%positions(iatom, idim) = sys%positions(iatom, idim) + sys%cell(idim)

          elseif (sys%positions(iatom, idim) > centerSys(idim) + halfCell(idim)  ) then
            pbcCorrectedSys%positions(iatom, idim) = sys%positions(iatom, idim) - sys%cell(idim)

          else
            pbcCorrectedSys%positions(iatom, idim) = sys%positions(iatom, idim)

          endif

        enddo
      enddo
    end function periodic

end module
