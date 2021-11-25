subroutine Solver(ni, nj, v, t, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, &
				  Re, Pr, eps, maxiter, mode)
implicit none

integer						:: i, j, ni, nj, iter, maxiter, mode
real(8)						:: D, Re, Pr, eps, dt
real(8) 					:: v(0:ni, 0:nj, 2), t(0:ni, 0:nj), R_m(ni-1, nj-1), div(0:ni, 0:nj), grad(0:ni, 0:nj, 2), laplacian(0:NI, 0:NJ), & 
							   cellvolume(ni-1, nj-1), cellcenter(0:ni, 0:nj, 2), &
							   iface_center(ni, nj-1, 2), iface_vector(ni, nj-1, 2), jface_center(ni-1, nj, 2), jface_vector(ni-1, nj, 2)

D = 1.0 / Re / Pr
dt = 1.0e-3
t(0:ni, 0:nj) = 0.5
R_m(ni-1, nj-1) = 1.0
iter = 0.0
call calcgrad_leastsquare(NI, NJ, t, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
call calcdiv(ni, nj, v, t, div, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, mode)
call calclaplacian(ni, nj, t, grad, laplacian, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)

do while (maxval(R_m) > eps .and. iter < maxiter)
	
	call bound_cond(ni, nj, t)
	call calcgrad_leastsquare(NI, NJ, t, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
	call calcdiv(ni, nj, v, t, div, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, mode)
	call calclaplacian(ni, nj, t, grad, laplacian, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
	print*, iter, maxval(grad), maxval(div), maxval(laplacian)
	
	do j = 1, nj - 1
		do i = 1, ni - 1
			R_m(i, j) = div(i, j) - D * laplacian(i, j)
			t(i, j) = t(i, j) - dt * R_m(i, j)
		end do
	end do
	
	! if (mod(iter, 10) == 0) then
		! print*, iter, maxval(t)
	! end if
	iter = iter + 1
	!write(*,*) 'Iter number =', iter, 'Res =', maxval(R_m)
	
end do

end subroutine

! biundary cond
subroutine bound_cond(ni, nj, t)
implicit none

integer 						:: ni, nj
real(8), dimension(0:ni, 0:nj) 	:: t

! bottom wall: t = 1.0
t(:, 0) = 1.0

! top wall: t = 0.0
t(:, nj) = 0.0

! side wall: dt_dn = 0.0
t(0, :) = t(1, :)
t(ni, :) = t(ni - 1, :)

end subroutine