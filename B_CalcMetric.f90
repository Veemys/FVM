subroutine calcmetric(NI, NJ, x, y, cellcenter, cellvolume, iface_center, iface_vector, jface_center, jface_vector)
implicit none

integer :: i, j, NI, NJ, nbound, ibound, iout, jbound, jout
real(8) :: x(NI, NJ), y(NI, NJ), &									! input: nodes coordinates
		   cellcenter(0:NI, 0:NJ, 2), cellvolume(NI-1, NJ-1), &		! output: cell centers and volumes
		   iface_center(NI, NJ-1, 2), iface_vector(NI, NJ-1, 2), &	! face centers and vectors for I-faces
		   jface_center(NI-1, NJ, 2), jface_vector(NI-1, NJ, 2)		! face centers and vectors for J-faces
real(8) :: r(2)

!=== FACE CENTERS AND FACE VECTORS ===
! I-DIRECTION
do j = 1, NJ - 1
	do i = 1, NI
	
		r(1) = x(i, j+1) - x(i, j)  	! r = vector from one node to another
		r(2) = y(i, j+1) - y(i, j)
		iface_vector(i, j, 1) = r(2) 	! IFaceVector = r rotated on 90 degree
		iface_vector(i, j, 2) = - r(1) 	! IFaceVector directed to increasing I-index
		iface_center(i, j, 1) = 0.5 * (x(i, j) + x(i, j+1))
		iface_center(i, j, 2) = 0.5 * (y(i, j) + y(i, j+1))
		
	end do
end do

! J-DIRECTION
do j = 1, NJ
    do i = 1, NI - 1
	
		r(1) = x(i+1, j) - x(i, j)  	! r = vector from one node to another
		r(2) = y(i+1, j) - y(i, j)
		jface_vector(i, j, 1) = - r(2) 	! JFaceVector = r rotated on -90 degree
		jface_vector(i, j, 2) = r(1) 	! JFaceVector directed to increasing J-index 
		jface_center(i, j, 1) = 0.5 * (x(i, j) + x(i+1, j))
		jface_center(i, j, 2) = 0.5 * (y(i, j) + y(i+1, j))
		
	end do
end do


!=== CELL VOLUMES ===
do j = 1, NJ - 1
	do i = 1, NI - 1
	
		r(1) = x(i+1, j+1) - X(i, j)
		r(2) = Y(i+1, j+1) - Y(i, j)
		cellvolume(i, j) = 0.5 * dot_product(iface_vector(i, j, :), r) & ! sum surfaces of two triangles
						   + 0.5 * dot_product(jface_vector(i, j, :), r)
	
	end do
end do


!=== CELL CENTERS ===
! FOR INNER CELLS: CENTER OF CONTOUR (sum of FaceCenter*FaceLength/Perimeter)
do j = 1, NJ - 1
	do i = 1, NI - 1
	
		cellcenter(i, j, :) = (iface_center(i, j, :) * Norm2(iface_vector(i, j, :)) + &
							  iface_center(i+1, j, :) * Norm2(iface_vector(i+1, j, :)) + &
                              jface_center(i, j, :) * Norm2(jface_vector(i, j, :)) + &
                              jface_center(i, j+1, :) * Norm2(jface_vector(i, j+1, :))) &
							  / (Norm2(iface_vector(i, j, :)) + Norm2(iface_vector(i+1, j, :)) + &
							  Norm2(jface_vector(i, j, :)) + Norm2(jface_vector(i, j+1, :)))
		
	end do
end do

! FOR DUMMY CELLS ON BOUNDARIES: CELL CENTER = FACE CENTER
! I-BOUNDARIES -----------------------------------------------------
do nbound = 1, 2

	if (nbound == 1) then
		ibound = 1; iout = 0
    else
		ibound = NI; iout = NI
    end if
	
    do j = 1, NJ - 1
		cellcenter(iout, j, :) = iface_center(ibound, j, :)
    end do
	
end do

! J-BOUNDARIES -----------------------------------------------------
do nbound = 1, 2

	if (nbound == 1) then
		jbound = 1; jout = 0
    else
		jbound = NJ; jout = NJ
    end if
	
    do i = 1, NI - 1
		cellcenter(i, jout, :) = jface_center(i, jbound, :)
    end do
	
end do

end subroutine