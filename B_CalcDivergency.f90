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
				case (0) 																				! Velocity divergency
					div(i, j) = div(i, j) + dot_product(vface(:), vectorface_cell(iface, :))
				case(1)																					! Central scheme
					pface = p_central_scheme(p(i, j), p(icell, jcell), d1, d2)
					div(i, j) = div(i, j) + dot_product(pface * vface(:), vectorface_cell(iface, :))
				case(2)																					! 1st upwind
					Gface = dot_product(vface(:), vectorface_cell(iface, :))
					pface = p_1stupwind_scheme(p(i, j), p(icell, jcell), Gface)
					div(i, j) = div(i, j) + dot_product(pface * vface(:), vectorface_cell(iface, :))
				case(3)																					! 2st upwind
					Gface = dot_product(vface(:), vectorface_cell(iface, :))
					d1_vec = centerface_cell(iface, :) - cellcenter(i, j, :)
					d2_vec = centerface_cell(iface, :) - cellcenter(icell, jcell, :)
					pface = p_2ndupwind_scheme(p(i, j), p(icell, jcell), d1_vec, d2_vec, grad(i, j, :), grad(icell, jcell, :), Gface)
					div(i, j) = div(i, j) + dot_product(pface * vface(:), vectorface_cell(iface, :))
			end select
			
		end do
		
		div(i, j) = div(i, j) / cellvolume(i, j)
		
	end do
end do

end subroutine

! Central scheme
real(8) function p_central_scheme(p1, p2, d1, d2)
implicit none

real(8) :: d1, d2, p1, p2, linear_interpolation

p_central_scheme = linear_interpolation(d1, d2, p1, p2)

end function

! 1st upwind scheme
real(8) function p_1stupwind_scheme(p1, p2, Gface)
implicit none

real(8) :: p1, p2, Gface

if (Gface >= 0.0) then
	p_1stupwind_scheme = p1
else
	p_1stupwind_scheme = p2
end if

end function

! 2nd upsind scheme
real(8) function p_2ndupwind_scheme(p1, p2, d1, d2, grad1, grad2, Gface)
implicit none

real(8) 				:: p1, p2, Gface
real(8), dimension(2) 	:: d1, d2, grad1, grad2

if (Gface >= 0.0) then
	p_2ndupwind_scheme = p1 + dot_product(d1, grad1)
else
	p_2ndupwind_scheme = p2 + dot_product(d2, grad2)
end if
	
end function