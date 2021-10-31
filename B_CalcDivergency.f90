subroutine calcdiv(NI, NJ, v, div, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
implicit none

integer						:: i, j, NI, NJ, iface, icell, jcell
real(8) 					:: v(0:NI, 0:NJ, 2), div(0:NI, 0:NJ), cellvolume(NI-1, NJ-1), cellcenter(0:NI, 0:NJ, 2), &
							   iface_center(NI, NJ-1, 2), iface_vector(NI, NJ-1, 2), jface_center(NI-1, NJ, 2), jface_vector(NI-1, NJ, 2), &
							   d1, d2, linear_interpolation
integer, dimension(4, 2) 	:: ncell(4, 2)
real(8), dimension(4, 2) 	:: centerneigh_cell, centerface_cell, vectorface_cell
real(8), dimension(2) 		:: vface

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
			
			div(i, j) = div(i, j) + dot_product(vface(:), vectorface_cell(iface, :))
			
		end do
		
		div(i, j) = div(i, j) / cellvolume(i, j)
		
	end do
end do

end Subroutine