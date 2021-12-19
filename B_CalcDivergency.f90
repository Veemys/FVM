subroutine calcdiv(NI, NJ, v, p, div, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, mode)
implicit none

integer						:: i, j, NI, NJ, iface, icell, jcell, mode
real(8) 					:: v(0:NI, 0:NJ, 2), p(0:NI, 0:NJ), div(0:NI, 0:NJ), grad(0:NI, 0:NJ, 2), & 
							   cellvolume(NI-1, NJ-1), cellcenter(0:NI, 0:NJ, 2), &
							   iface_center(NI, NJ-1, 2), iface_vector(NI, NJ-1, 2), jface_center(NI-1, NJ, 2), jface_vector(NI-1, NJ, 2), &
							   d1, d2, linear_interpolation, p_central_scheme, p_1stupwind_scheme, p_2ndupwind_scheme, Gface, pface
integer, dimension(4, 2) 	:: ncell(4, 2)
real(8), dimension(4, 2) 	:: centerneigh_cell, centerface_cell, vectorface_cell
real(8), dimension(2) 		:: vface, d1_vec, d2_vec

div = 0.0

do j = 1, NJ - 1
	do i = 1, NI - 1
	
		call cell_information(i, j, NI, NJ, cellcenter, iface_center, iface_vector, jface_center, jface_vector, ncell, &
							  centerneigh_cell, centerface_cell, vectorface_cell)
				
		do iface = 1, 4
		
			icell = ncell(iface, 1)
			jcell = ncell(iface, 2)
			d1 = Norm2(centerface_cell(iface, :) - cellcenter(i, j, :))
			d2 = Norm2(centerface_cell(iface, :) - cellcenter(icell, jcell, :))
			vface(1) = linear_interpolation(d1, d2, v(i, j, 1), v(icell, jcell, 1))
			vface(2) = linear_interpolation(d1, d2, v(i, j, 2), v(icell, jcell, 2))
			
			select case (mode)
			
				! velocity divergency
				case (0)
				
					div(i, j) = div(i, j) + dot_product(vface(:), vectorface_cell(iface, :))
					
				! Central scheme
				case(1)
				
					pface = linear_interpolation(d1, d2, p(i, j), p(icell, jcell))
					div(i, j) = div(i, j) + dot_product(pface * vface(:), vectorface_cell(iface, :))
					
				! 1st upwind
				case(2)
				
					Gface = dot_product(vface(:), vectorface_cell(iface, :))
					if (Gface >= 0.0) then
						pface = p(i, j)
					else
						pface = p(icell, jcell)
						if (d2 < 1e-6) pface = 2.0 * p(icell, jcell) - p(i, j)
					end if
					div(i, j) = div(i, j) + dot_product(pface * vface(:), vectorface_cell(iface, :))
					
				! 2st upwind
				case(3)
				
					Gface = dot_product(vface(:), vectorface_cell(iface, :))
					d1_vec(:) = centerface_cell(iface, :) - cellcenter(i, j, :)
					d2_vec(:) = centerface_cell(iface, :) - cellcenter(icell, jcell, :)
					if (Gface >= 0.0) then
						pface = p(i, j) + dot_product(d1_vec(:), grad(i, j, :))
					else
						pface = p(icell, jcell) + dot_product(d2_vec(:), grad(icell, jcell, :))
						if (d2 < 1e-6) pface = 2.0 * p(icell, jcell) - p(i, j) + 4.0 * (p(i, j) - p(icell, jcell)) &
											   - 3.0 * dot_product(grad(i, j, :), -d1_vec(:))
					end if
					div(i, j) = div(i, j) + dot_product(pface * vface(:), vectorface_cell(iface, :))
					
			end select
			
		end do
		
		div(i, j) = div(i, j) / cellvolume(i, j)
		
	end do
end do

end subroutine