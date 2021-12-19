subroutine Solver(ni, nj, v, t, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, &
				  Re, Pr, eps, maxiter, mode, gradient_scheme, maxiter_correct)
implicit none

integer		:: i, j, ni, nj, iter, maxiter, mode, gradient_scheme, maxiter_correct
real(8)		:: D, Re, Pr, eps, dt, R_m_p
real(8) 	:: v(0:ni, 0:nj, 2), t(0:ni, 0:nj), R_m(0:ni, 0:nj), div(0:ni, 0:nj), grad(0:ni, 0:nj, 2), laplacian(0:NI, 0:NJ), & 
			   cellvolume(ni-1, nj-1), cellcenter(0:ni, 0:nj, 2), &
			   iface_center(ni, nj-1, 2), iface_vector(ni, nj-1, 2), jface_center(ni-1, nj, 2), jface_vector(ni-1, nj, 2)

D = 1.0 / Re / Pr
dt = 0.002
t(0:ni, 0:nj) = 0.5
R_m(0:ni, 0:nj) = 1.0
iter = 0
call bound_cond(ni, nj, t)
call calcgrad_choice_scheme(ni, nj, gradient_scheme, maxiter_correct, t, grad, cellvolume, cellcenter, &
							iface_center, iface_vector, jface_center, jface_vector)
call calcdiv(ni, nj, v, t, div, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, mode)
call calclaplacian(ni, nj, t, grad, laplacian, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)

do while (maxval(abs(R_m)) > eps .and. iter < maxiter)
	
	call bound_cond(ni, nj, t)
	call calcgrad_leastsquare(NI, NJ, t, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
	call calcdiv(ni, nj, v, t, div, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, mode)
	call calclaplacian(ni, nj, t, grad, laplacian, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
	
	R_m = div - D * laplacian

	do j = 1, nj - 1
		do i = 1, ni - 1
		
			! call choice_time_step(i, j, ni, nj, v, D, cellvolume, cellcenter, iface_center, iface_vector, jface_center, jface_vector, dt)
			
			t(i, j) = t(i, j) - dt * R_m(i, j)
			
		end do
	end do
	
	iter = iter + 1
	if (mod(iter, 1) == 0) write(*,*) 'Iter number =', iter, 'Res =', maxval(abs(R_m))
	
end do

end subroutine

! boundary conditions
subroutine bound_cond(ni, nj, t)
implicit none

integer 						:: ni, nj
real(8), dimension(0:ni, 0:nj) 	:: t

! bottom and top wall: dt_dn = 0.0
t(:, 0) = t(:, 1)
t(:, nj) = t(:, nj - 1)

! right wall: t = 1.0
t(ni, :) = 1.0

! left wall: t = 0.0
t(0, :) = 0.0

end subroutine

! choice of time step
subroutine choice_time_step(i, j, ni, nj, v, D, cellvolume, cellcenter, iface_center, iface_vector, jface_center, jface_vector, dt)
implicit none

integer 			:: i, j, ni, nj, iface, icell, jcell
integer				:: ncell(4, 2)
real(8), parameter 	:: CFL = 0.9, VNM = 0.4
real(8)				:: dt, sum1, sum2, d1, d2, D, dt_conv, dt_diff, linear_interpolation
real(8) 			:: v(0:ni, 0:nj, 2), cellvolume(ni-1, nj-1), cellcenter(0:ni, 0:nj, 2), &
					   iface_center(ni, nj-1, 2), iface_vector(ni, nj-1, 2), jface_center(ni-1, nj, 2), jface_vector(ni-1, nj, 2), &
					   centerneigh_cell(4, 2), centerface_cell(4, 2), vectorface_cell(4, 2), v_face(2)

call cell_information(i, j, ni, nj, cellcenter, iface_center, iface_vector, jface_center, jface_vector, &
					  ncell, centerneigh_cell, centerface_cell, vectorface_cell)

sum1 = 0.0; sum2 = 0.0
do iface = 1, 4
	
	icell = ncell(iface, 1)
	jcell = ncell(iface, 2)
	d1 = Norm2(centerface_cell(iface, :) - cellcenter(i, j, :))
	d2 = Norm2(centerface_cell(iface, :) - cellcenter(icell, jcell, :))
	v_face(1) = linear_interpolation(d1, d2, v(i, j, 1), v(icell, jcell, 1))
	v_face(2) = linear_interpolation(d1, d2, v(i, j, 2), v(icell, jcell, 2))
	
	sum1 = sum1 + dot_product(abs(v_face(:)), abs(vectorface_cell(iface, :)))
	sum2 = sum2 + D / Norm2(vectorface_cell(iface, :))
	
end do

dt_conv = CFL / (sum1 / cellvolume(i, j))
dt_diff = VNM / (2.0 * sum2)

dt = 1.0 / (1.0 / dt_conv + 1.0 / dt_diff)

end subroutine