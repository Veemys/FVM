! Function of pressure
real(8) function pressure(x, y)
implicit none

real(8) :: x, y

pressure = x**3.0 + y**3.0

end function

! Function of velocity
real(8) function velocity(x, y, coord)
implicit none

integer :: coord
real(8) :: x, y

if (coord == 1) then
	velocity = 1 + x * y**2.0
else
	velocity = 1 + y * x
end if

end function 

! Calculation of exact values of gradient
real(8) function calcgrad_exact(x, y, coord)
implicit none

integer :: coord
real(8) :: x, y

if (coord == 1) then
	calcgrad_exact = 3.0 * x**2.0
else
	calcgrad_exact = 3.0 * y**2.0
end if

end function

! Calculation of exact values of divergency
real(8) function calcdiv_exact(x, y)
implicit none

real(8) :: x, y

calcdiv_exact = y**2.0 + x

end function

! Calculation of exact values of rotor
real(8) function calcrot_exact(x, y)
implicit none

real(8) :: x, y

calcrot_exact = y - 2.0 * x * y ! dVy/dx - dVx/dy

end function

! Calculation of error of gradient
subroutine calcgrad_error(NI, NJ, grad, gradexact, graderror)
implicit none

integer 							:: i, j, NI, NJ
real(8), dimension(0:NI, 0:NJ, 2)	:: grad, gradexact, graderror

do j = 0, NJ
	do i = 0, NI
		graderror(i, j, :) = abs(grad(i, j, :) - gradexact(i, j, :)) / maxval(gradexact(i, j, :))
	end do
end do

print*, 'Max Error of Gradient =', maxval(graderror(1:NI-1, 1:NJ-1, :))

end subroutine

! Calculation of error of divergency
subroutine calcdiv_error(NI, NJ, div, divexact, diverror)
implicit none

integer 						:: i, j, NI, NJ
real(8), dimension(0:NI, 0:NJ) 	:: div, divexact, diverror

do j = 0, NJ
	do i = 0, NI
		diverror(i, j) = abs(div(i, j) - divexact(i, j)) / divexact(i, j)
	end do
end do

print*, 'Max Error of Divergency =', maxval(diverror(1:NI-1, 1:NJ-1))

end subroutine

! Calculation of error of rotor
subroutine calcrot_error(NI, NJ, rot, rotexact, roterror)
implicit none

integer 						:: i, j, NI, NJ
real(8), dimension(0:NI, 0:NJ) 	:: rot, rotexact, roterror

do j = 0, NJ
	do i = 0, NI
		roterror(i, j) = abs(rot(i, j) - rotexact(i, j)) / rotexact(i, j)
	end do
end do

print*, 'Max Error of Rotor =', maxval(roterror(1:NI-1, 1:NJ-1))

end subroutine

! Linear Interpolation
real(8) function linear_interpolation(d1, d2, x1, x2)
implicit none

real(8) :: x1, x2, d1, d2

linear_interpolation = (x1 * d2 + x2 * d1) / (d1 + d2)

end function

! Information about given cell
subroutine cell_information(i, j, NI, NJ, cellcenter, iface_center, iface_vector, jface_center, jface_vector, &
							ncell, centerneigh_cell, centerface_cell, vectorface_cell)
implicit none

integer 					:: i, j, NI, NJ
real(8) 					:: cellcenter(0:NI, 0:NJ, 2), iface_center(NI, NJ-1, 2), iface_vector(NI, NJ-1, 2), & 
							   jface_center(NI-1, NJ, 2), jface_vector(ni-1, nj, 2)
integer, dimension(4, 2) 	:: ncell
real(8), dimension(4,2) 	:: centerneigh_cell, centerface_cell, vectorface_cell

! Numbers of neighbour cells
ncell(1, :) = [i-1, j]
ncell(2, :) = [i, j+1]
ncell(3, :) = [i+1, j]
ncell(4, :) = [i, j-1]

! Centers of neighbour cells
centerneigh_cell(1, :) = cellcenter(i-1, j, :)
centerneigh_cell(2, :) = cellcenter(i, j+1, :)
centerneigh_cell(3, :) = cellcenter(i+1, j, :)
centerneigh_cell(4, :) = cellcenter(i, j-1, :)

! Centers of faces of a given cell
centerface_cell(1, :) = iface_center(i, j, :)
centerface_cell(2, :) = jface_center(i, j+1, :)
centerface_cell(3, :) = iface_center(i+1, j, :)
centerface_cell(4, :) = jface_center(i, j, :)

! Vectors of faces of a given cell
vectorface_cell(1, :) = - iface_vector(i, j, :)
vectorface_cell(2, :) = jface_vector(i, j+1, :)
vectorface_cell(3, :) = iface_vector(i+1, j, :)
vectorface_cell(4, :) = - jface_vector(i, j, :)

end subroutine