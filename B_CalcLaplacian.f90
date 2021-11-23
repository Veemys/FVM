subroutine calclaplacian(NI, NJ, p, grad, laplacian, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
implicit none

integer						:: i, j, NI, NJ, iface, icell, jcell
real(8) 					:: p(0:NI, 0:NJ), grad(0:NI, 0:NJ, 2), laplacian(0:NI, 0:NJ), cellvolume(NI-1, NJ-1), cellcenter(0:NI, 0:NJ, 2), &
							   iface_center(NI, NJ-1, 2), iface_vector(NI, NJ-1, 2), jface_center(NI-1, NJ, 2), jface_vector(NI-1, NJ, 2), &
							   d1, d2, linear_interpolation, dp_dn_face
integer, dimension(4, 2) 	:: ncell(4, 2)
real(8), dimension(4, 2) 	:: centerneigh_cell, centerface_cell, vectorface_cell
real(8), dimension(2) 		:: gradinterpol, vectorcell_cell, ivec_face, ivec_cell

laplacian = 0.0

do j = 1, NJ - 1
	do i = 1, NI - 1
	
		call cell_information(i, j, NI, NJ, cellcenter, iface_center, iface_vector, jface_center, jface_vector, ncell, &
							  centerneigh_cell, centerface_cell, vectorface_cell)
				
		do iface = 1, 4
		
			icell = ncell(iface, 1)
			jcell = ncell(iface, 2)
			d1 = Norm2(centerface_cell(iface, :) - cellcenter(i, j, :))
			d2 = Norm2(centerface_cell(iface, :) - cellcenter(icell, jcell, :))
			
			vectorcell_cell(:) = cellcenter(icell, jcell, :) - cellcenter(i, j, :)
			ivec_cell(:) = vectorcell_cell(:) / Norm2(vectorcell_cell(:))
			ivec_face(:) = vectorface_cell(iface, :) / Norm2(vectorface_cell(iface, :))
			
			if (d2 < 1e-7) then
				gradinterpol(:) = grad(icell, jcell, :)
				dp_dn_face = 5.0 / 3.0 * (p(icell, jcell) - p(i, j)) / Norm2(vectorcell_cell(:)) - &
							 2.0 / 3.0 * dot_product(grad(i, j, :), ivec_face(:)) + &
							 dot_product(ivec_face(:) - ivec_cell(:), gradinterpol(:))
			else
				gradinterpol(1) = linear_interpolation(d1, d2, grad(i, j, 1), grad(icell, jcell, 1))
				gradinterpol(2) = linear_interpolation(d1, d2, grad(i, j, 2), grad(icell, jcell, 2))
				dp_dn_face = (p(icell, jcell) - p(i, j)) / Norm2(vectorcell_cell(:)) + & 
							 dot_product(ivec_face(:) - ivec_cell(:), gradinterpol(:))
			end if
			
			laplacian(i, j) = laplacian(i, j) + dp_dn_face * Norm2(vectorface_cell(iface, :))
			
		end do
		
		laplacian(i, j) = laplacian(i, j) / cellvolume(i, j)
				
	end do
end do

end subroutine