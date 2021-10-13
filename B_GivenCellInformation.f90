subroutine B_GivenCellInformation(i, j, ni, nj, CellCenter, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector, &
								  numCell, cenNeighCell, cenFaceCell, vecFaceCell)
implicit none
integer :: i, j, ni, nj
integer :: numCell(4, 2)
real 	:: CellCenter(0:ni, 0:nj, 2), IFaceCenter(ni, nj-1, 2), IFaceVector(ni, nj-1, 2), JFaceCenter(ni-1, nj, 2), &
		   JFaceVector(ni-1, nj, 2)
real, dimension(4,2) :: cenNeighCell, cenFaceCell, vecFaceCell

! Numbers of neighbour cells
numCell(1, :) = [i-1, j]
numCell(2, :) = [i, j+1]
numCell(3, :) = [i+1, j]
numCell(4, :) = [i, j-1]

! Centers of neighbour cells
cenNeighCell(1, :) = CellCenter(i-1, j, :)
cenNeighCell(2, :) = CellCenter(i, j+1, :)
cenNeighCell(3, :) = CellCenter(i+1, j, :)
cenNeighCell(4, :) = CellCenter(i, j-1, :)

! Centers of faces of a given cell
cenFaceCell(1, :) = IFaceCenter(i, j, :)
cenFaceCell(2, :) = JFaceCenter(i, j+1, :)
cenFaceCell(3, :) = IFaceCenter(i+1, j, :)
cenFaceCell(4, :) = JFaceCenter(i, j, :)

! Vectors of faces of a given cell
vecFaceCell(1, :) = - IFaceVector(i, j, :)
vecFaceCell(2, :) = JFaceVector(i, j+1, :)
vecFaceCell(3, :) = IFaceVector(i+1, j, :)
vecFaceCell(4, :) = - JFaceVector(i, j, :)

end subroutine