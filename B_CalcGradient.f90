subroutine B_CalcGradient(NI, NJ, P, GradP, CellVolume, CellCenter, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
implicit none
integer :: i, j, NI, NJ, IFace, iCell, jCell
integer :: numberCell(4, 2)
real 	:: P(0:NI, 0:NJ), GradP(0:NI, 0:NJ, 2), CellVolume(NI-1, NJ-1), CellCenter(0:NI, 0:NJ, 2), &
		   IFaceCenter(NI, NJ-1, 2), IFaceVector(NI, NJ-1, 2), JFaceCenter(NI-1, NJ, 2), JFaceVector(NI-1, NJ, 2)
real, dimension(4,2) :: centerNeighbourCell, centerFaceCell, vectorFaceCell
real	:: d1, d2, phiFace, LinearInterpol

do j = 1, NJ - 1
	do i = 1, NI - 1
	
		call B_GivenCellInformation(i, j, ni, nj, CellCenter, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector, numberCell, &
									centerNeighbourCell, centerFaceCell, vectorFaceCell)
									
		do iFace = 1, 4
			iCell = numberCell(iFace, 1)
			jCell = numberCell(iFace, 2)
			d1 = Norm2(centerFaceCell(iFace, :) - CellCenter(i, j, :))
			d2 = Norm2(centerFaceCell(iFace, :) - centerNeighbourCell(iFace, :))
			phiFace = LinearInterpol(d1, d2, P(i, j), P(iCell, jCell))
			GradP(i, j, :) = GradP(i, j, :) + phiFace * vectorFaceCell(iFace, :)
		end do
		GradP(i, j, :) = GradP(i, j, :) / CellVolume(i, j)
		
		! GradP(i, j, :) = (- (P(i-1, j) + P(i, j)) / 2.0 * IFaceVector(i, j, :) + &
						 ! (P(i+1, j) + P(i, j)) / 2.0 * IFaceVector(i+1, j, :) - &
						 ! (P(i, j-1) + P(i, j)) / 2.0 * JFaceVector(i, j, :) + &
						 ! (P(i, j+1) + P(i, j)) / 2.0 * JFaceVector(i, j+1, :)) &
						 ! / CellVolume(i, j)
	end do
end do

end Subroutine