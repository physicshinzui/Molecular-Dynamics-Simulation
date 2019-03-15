module PBC
  implicit none
contains
!---this is currently specialised only for FCC lattice
!    what I want to do are
!        (1) calculate center of geormetry (cog) of a system
!        (2) calculate center of geormetry (cog) of a system
!    (This would be a code by python)
    function periodic(coords, cell, centreOfGeometry) result(coordsPbc)
      real(kind = 8), intent(in) :: coords(:,:), cell(:), centreOfGeometry(:)
      real(kind = 8), allocatable :: coordsPbc(:,:), halfCell(:)
      integer :: iatom, idim

      coordsPbc = coords
      coordsPbc(:,:) = 0.0d0
      halfCell = cell * 0.5d0

      do idim = 1, size(coords, 2)
        do iatom = 1, size(coords, 1)

          if (coords(iatom, idim) < 0.0d0 ) then
            coordsPbc(iatom, idim) = coords(iatom, idim) + cell(idim)

          elseif (coords(iatom, idim) > cell(idim)) then
            coordsPbc(iatom, idim) = coords(iatom, idim) - cell(idim)

          else
            coordsPbc(iatom, idim) = coords(iatom, idim)

          endif

!---I wanna generalise this function in future
!          if (coords(iatom, idim) > centreOfGeometry(idim) + halfCell(idim) ) then
!            coordsPbc(iatom, idim) = coords(iatom, idim) - cell(idim)
!
!          elseif (coords(iatom, idim) < - centreOfGeometry(idim) + halfCell(idim)) then
!            coordsPbc(iatom, idim) = coords(iatom, idim) + cell(idim)
!
!          else
!            coordsPbc(iatom, idim) = coords(iatom, idim)
!
!          endif

        enddo
      enddo
    end function periodic

end module 
