! calculation of the gradient field via Green-Gauss method
subroutine calcgrad_greengauss(NI, NJ, p, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
implicit none

integer 					:: i, j, NI, NJ, iface, icell, jcell
real(8) 					:: p(0:NI, 0:NJ), grad(0:NI, 0:NJ, 2), cellvolume(NI-1, NJ-1), cellcenter(0:NI, 0:NJ, 2), &
							   iface_center(NI, NJ-1, 2), iface_vector(NI, NJ-1, 2), jface_center(NI-1, NJ, 2), jface_vector(NI-1, NJ, 2), &
							   d1, d2, pface, linear_interpolation, pface_correct
integer, dimension(4, 2)	:: ncell(4, 2)
real(8), dimension(4, 2) 	:: centerneigh_cell, centerface_cell, vectorface_cell
real(8), dimension(2) 		:: errorcoord, gradinterpol, gradcorrect

do j = 1, NJ - 1
	do i = 1, NI - 1
	
		call cell_information(i, j, ni, nj, cellcenter, iface_center, iface_vector, jface_center, jface_vector, ncell, &
							  centerneigh_cell, centerface_cell, vectorface_cell)
		
		gradcorrect(:) = 0.0
		
		do iface = 1, 4
		
			icell = ncell(iface, 1)
			jcell = ncell(iface, 2)
			d1 = Norm2(centerface_cell(iface, :) - cellcenter(i, j, :))
			d2 = Norm2(centerface_cell(iface, :) - cellcenter(icell, jcell, :))
			pface = linear_interpolation(d1, d2, p(i, j), p(icell, jcell))
			
			errorcoord(1) = linear_interpolation(d1, d2, cellcenter(i, j, 1), cellcenter(icell, jcell, 1))
			errorcoord(2) = linear_interpolation(d1, d2, cellcenter(i, j, 2), cellcenter(icell, jcell, 2))
			gradinterpol(1) = linear_interpolation(d1, d2, grad(i, j, 1), grad(icell, jcell, 1))
			gradinterpol(2) = linear_interpolation(d1, d2, grad(i, j, 2), grad(icell, jcell, 2))
			pface_correct = pface + dot_product(centerface_cell(iface, :) - errorcoord(:), gradinterpol(:))
			
			gradcorrect(:) = gradcorrect(:) + pface_correct * vectorface_cell(iface, :)
			
		end do
		
		grad(i, j, :) = gradcorrect(:) / cellvolume(i, j)
		
	end do
end do

end subroutine

! calculation of the gradient field via least square method
subroutine calcgrad_leastsquare(NI, NJ, p, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
implicit none

integer 					:: i, j, NI, NJ, iface, icell, jcell, m, n
real(8), parameter			:: eps = 1e-16
real(8) 					:: p(0:NI, 0:NJ), grad(0:NI, 0:NJ, 2), cellvolume(NI-1, NJ-1), cellcenter(0:NI, 0:NJ, 2), &
							   iface_center(NI, NJ-1, 2), iface_vector(NI, NJ-1, 2), jface_center(NI-1, NJ, 2), jface_vector(NI-1, NJ, 2), &
							   d1, d2, pface, linear_interpolation, pface_correct
integer, dimension(4, 2)	:: ncell(4, 2)
real(8), dimension(4, 2) 	:: centerneigh_cell, centerface_cell, vectorface_cell
real(8), dimension(2, 2) 	:: A
real(8), dimension(2) 		:: r, b

do j = 1, NJ - 1
	do i = 1, NI - 1
		
		call cell_information(i, j, ni, nj, cellcenter, iface_center, iface_vector, jface_center, jface_vector, ncell, &
							  centerneigh_cell, centerface_cell, vectorface_cell)
		
		A = 0.0
		b = 0.0
		
		do iface = 1, 4
			
			icell = ncell(iface, 1)
			jcell = ncell(iface, 2)
			r(:) = cellcenter(icell, jcell, :) - cellcenter(i, j, :)
			do n = 1, 2
				do m = 1, 2
					A(n, m) = A(n, m) + r(n) * r(m) / Norm2(r)
				end do
			end do
			b(:) = b(:) + r(:) * (p(icell, jcell) - p(i, j)) / Norm2(r)
			
		end do
		
		call gaussseidel_method(2, A, b, eps, grad(i, j, :))
		
	end do
end do

end subroutine

! Gauss-Seidel method
subroutine gaussseidel_method(n, A, b, eps, x)
implicit none

integer 					:: i, j, n
real(8) 					:: sum1, sum2, eps
real(8), dimension(n)		:: b, x, x_new
real(8), dimension(n, n)	:: A

x = 0.0
x_new = 0.1

do while (Norm2(x_new - x) >= eps)
	
	x = x_new

	do i = 1, n
	
		sum1 = 0.0
		do j = 1, i - 1
			sum1 = sum1 + A(i, j) * x_new(j)
		end do
		
		sum2 = 0.0
		do j = i + 1, n
			sum2 = sum2 + A(i, j) * x(j)
		end do
		
		x_new(i) = (b(i) - sum1 - sum2) / A(i, i)
		
	end do
	
end do

end subroutine

! choice of Grad Scheme
subroutine calcgrad_choice_scheme(ni, nj, gradient_scheme, maxiter_correct, p, grad, cellvolume, cellcenter, &
								  iface_center, iface_vector, jface_center, jface_vector)

integer :: ni, nj, gradient_scheme, maxiter_correct
real(8) :: p(0:NI, 0:NJ), grad(0:NI, 0:NJ, 2), cellvolume(NI-1, NJ-1), cellcenter(0:NI, 0:NJ, 2), &
			 iface_center(NI, NJ-1, 2), iface_vector(NI, NJ-1, 2), jface_center(NI-1, NJ, 2), jface_vector(NI-1, NJ, 2)

grad = 0.0
select case (gradient_scheme)
	case (1)
		do k = 1, maxiter_correct
			call calcgrad_greengauss(ni, nj, p, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
		end do
	case (2)
		call calcgrad_leastsquare(ni, nj, p, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
end select

end subroutine