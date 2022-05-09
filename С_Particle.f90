!*****************************************************************************************************************************
subroutine solver_part(rho_part, d_part, x_0_part, y_0_part, u_0_part, v_0_part, omega_0_part, rho_envi, x, y, v, &
	cellvolume, iface_vector, jface_vector, dt, nt, ni, nj, mu, &
	x_part, y_part, u_part, v_part, omega_part, i_part, j_part, &
	rotV, gradP)
	implicit none

	integer, parameter :: io_trac = 69, io_forces = 228
	integer :: iter, nt, ni, nj, i_part_old, j_part_old, i_part, j_part, ip1, jp1, St
	real(8), parameter :: pi = 4.0 * atan(1.0)
	real(8) :: x_part, y_part, u_part, v_part, omega_part, x_0_part, y_0_part, u_0_part, v_0_part, omega_0_part, F_x, F_y, M_z, &
			v_ref, l_ref, Sk, rho_part, d_part, rho_envi, mu, u_envi, v_envi, dt
	real(8) :: iface_vector(ni, nj-1, 2), jface_vector(ni-1, nj, 2)
	real(8) :: x_part_new, y_part_new, u_part_new, v_part_new, omega_part_new
	real(8) :: v(0:ni, 0:nj, 2), x(ni, nj), y(ni, nj), cellvolume(ni-1, nj-1), v_rel(0:nt, 2)
	real(8) :: x_cross, y_cross, x_part_reflect, y_part_reflect
	real(8) :: rotV(0:ni, 0:nj), gradP(0:ni, 0:nj, 2)
	real(8) :: F_Stokes(2), F_Buoyancy(2), F_Magnus(2), F_Saffman(2), F_Basset(2), F_Add_mass(2)
	real(8) :: u_part_old, v_part_old, u_envi_old, v_envi_old
	real(8) :: calc_abs_vector				! functions

	v_ref = 0.05
	l_ref = 0.05
	Sk = rho_part * d_part**2 * v_ref / (mu * l_ref)
	print*, 'Stokes number = ', Sk

	open(io_trac, file = 'part_trac.plt')
	write(io_trac,*) 'VARIABLES = "X", "Y"'
	open(io_forces, file = 'forces.plt')
	write(io_forces,*) 'VARIABLES = "t", "F_Stokes", "F_Buoyancy", "F_Magnus", "F_Saffman", "F_Basset", "F_Add_mass", "F_SUMM"'

	x_part = x_0_part
	y_part = y_0_part
	u_part = u_0_part
	v_part = v_0_part
	omega_part = omega_0_part

	F_Stokes = 0.0
	F_Buoyancy = 0.0
	F_Magnus = 0.0
	F_Saffman = 0.0
	F_Basset = 0.0
	F_Add_mass = 0.0

	u_part_old = u_part
	v_part_old = v_part
	u_envi_old = u_envi
	v_envi_old = v_envi

	F_x = 0.0
	F_y = 0.0
	M_z = 0.0
	! integral_basset_force = 0.0

	call find_part_location(x_part, y_part, ni, nj, x, y, cellvolume, i_part, j_part)
	i_part_old = -1
	j_part_old = -1

	v_rel = 0.0;
	u_envi = v(i_part, j_part, 1)
	v_envi = v(i_part, j_part, 2)
	v_rel(0, 1) = u_envi - u_part
	v_rel(0, 2) = v_envi - v_part

	call write_part_trac(x_part, y_part, io_trac)
	call write_forces(dt, F_Stokes, F_Buoyancy, F_Magnus, F_Saffman, F_Basset, F_Add_mass, F_x, F_y, io_forces, &
		i_part_old, j_part_old, i_part, j_part)

	do iter = 1, nt
		
		! calc index of cell with particle
		i_part_old = i_part
		j_part_old = j_part
		i_part = - 1
		j_part = - 1
		call find_part_location(x_part, y_part, ni, nj, x, y, cellvolume, i_part, j_part)
		
		u_envi = v(i_part, j_part, 1)
		v_envi = v(i_part, j_part, 2)
		v_rel(iter, 1) = u_envi - u_part
		v_rel(iter, 2) = v_envi - v_part
		
		call calc_force(nt, dt, iter, ni, nj, i_part, j_part, rho_envi, rho_part, mu, d_part, u_part, v_part, omega_part, &
			u_envi, v_envi, v_rel, rotV, gradP, F_x, F_y, M_z, &
			u_part_old, v_part_old, u_envi_old, v_envi_old, &
			F_Stokes, F_Buoyancy, F_Magnus, F_Saffman, F_Basset, F_Add_mass)
		
		x_part_new = x_part + dt * u_part
		y_part_new = y_part + dt * v_part
		u_part_new = u_part + dt * F_x
		v_part_new = v_part + dt * F_y
		omega_part_new = omega_part + dt * M_z

		! Check boundary
		i_part_old = i_part
		j_part_old = j_part
		call find_part_location(x_part_new, y_part_new, ni, nj, x, y, cellvolume, i_part, j_part)
		
		if (i_part == - 1 .and. j_part == - 1) then
			!write(*,*) 'u_old = ', u_part_new, 'v_old = ', v_part_new
			call C_Boundary(x, y, iface_vector, jface_vector, ni, nj, x_part, y_part, u_part, v_part, omega_part, d_part,&
				x_part_new, y_part_new, u_part_new, v_part_new, omega_part_new, x_cross, y_cross, ip1, jp1, St, &
				x_part_reflect, y_part_reflect)
			!write(*,*) 'u_new = ', u_part_new, 'v_new = ', v_part_new
			call write_part_trac(x_cross, y_cross, io_trac)
			!call write_forces(iter * dt, F_Stokes, F_Buoyancy, F_Magnus, F_Saffman, F_Basset, F_Add_mass, F_x, F_y, io_forces, &
			!	i_part_old, j_part_old, i_part, j_part)

			if (St == 1 .or. St == 2) then
				! call write_part_info(x_part, y_part, iter * dt, F_Stokes, F_Buoyancy, F_Magnus, F_Saffman, F_Basset, F_Add_mass, F_x, F_y, 69)
				x_part_new = x_part_reflect
				y_part_new = y_part_reflect
			end if

			if (St == 3 .or. St == 4) then
				! call write_part_info(x_part, y_part, iter * dt, F_Stokes, F_Buoyancy, F_Magnus, F_Saffman, F_Basset, F_Add_mass, F_x, F_y, 69)
				return
			end if

		end if
		
		u_part_old = u_part
		v_part_old = v_part
		u_envi_old = u_envi
		v_envi_old = v_envi

		x_part = x_part_new
		y_part = y_part_new
		u_part = u_part_new
		v_part = v_part_new
		omega_part = omega_part_new

		call write_part_trac(x_part, y_part, io_trac)
		call write_forces(iter * dt, F_Stokes, F_Buoyancy, F_Magnus, F_Saffman, F_Basset, F_Add_mass, F_x, F_y, io_forces, &
			i_part_old, j_part_old, i_part, j_part)
		
		if (mod(iter, 100) == 0) then
			print*, 'time =', iter * dt, 'i =', i_part, 'j =', j_part, 'x =', x_part, 'y =', y_part
			!print*, F_x, F_y, M_z
		end if

	end do
	close(io_trac)
	close(io_forces)

end subroutine

!*****************************************************************************************************************************
subroutine find_part_location(x_part, y_part, ni, nj, x, y, cellvolume, i_part, j_part)
	implicit none

	integer :: ni, nj, i_part, j_part, i, j
	real(8) :: x_part, y_part, x(ni, nj), y(ni, nj), cellvolume(ni-1, nj-1), s, eps, triang_square

	i_part = - 1
	j_part = - 1
	eps = 1e-13

	do j = 1, nj - 1
		do i = 1, ni - 1
			
			s = triang_square(x(i, j), y(i, j), x(i, j+1), y(i, j+1), x_part, y_part) + &
			triang_square(x(i+1, j+1), y(i+1, j+1), x(i, j+1), y(i, j+1), x_part, y_part) + &
			triang_square(x(i, j), y(i, j), x(i+1, j), y(i+1, j), x_part, y_part) + &
			triang_square(x(i+1, j+1), y(i+1, j+1), x(i+1, j), y(i+1, j), x_part, y_part)
				
			if (abs(s - cellvolume(i, j)) <= eps) then
				i_part = i
				j_part = j
				return
			end if
			
		end do
	end do
		
end subroutine

!*****************************************************************************************************************************
subroutine C_Boundary(x, y, iface_vector, jface_vector, ni, nj, x_part, y_part, u_part, v_part, omega_part, d_part, &
						x_part_new, y_part_new, u_part_new, v_part_new, omega_part_new, x_cross, y_cross, ip1, jp1, St, &
						x_part_reflect, y_part_reflect)
	implicit none

	integer :: ni, nj, ip1, jp1, St, i, j, icr ! icr - index + cross
	real(8) :: x(ni, nj), y(ni, nj)
	real(8) :: iface_vector(ni, nj-1, 2), jface_vector(ni-1, nj, 2)
	real(8) :: d_part
	real(8) :: x_part_new, y_part_new, u_part_new, v_part_new, omega_part_new, x_part, y_part, u_part, v_part, omega_part
	real(8) :: x_cross, y_cross, x_part_reflect, y_part_reflect

	St = 0
	icr = 0

	do i = 1, ni - 1

		call cross_edges(x_part, y_part, x_part_new, y_part_new, x(i, 1), y(i, 1), x(i+1, 1), y(i+1, 1), x_cross, y_cross, icr)
		if (icr == 1) then
			St = 1
			call reflect(x_part, y_part, x_part_new, y_part_new, x_part_reflect, y_part_reflect, x(i, 1), y(i, 1), x(i+1, 1), y(i+1, 1))
			write(*,*) 'Particle reflected from boundary', St, 'and now part. pos. is', x_part_reflect, y_part_reflect
			call calc_velo_after_reflect(x(i, 1), y(i, 1), x(i+1, 1), y(i+1, 1), u_part, v_part, omega_part, u_part_new, v_part_new, &
				omega_part_new, d_part)
			return
		end if
		
		call cross_edges(x_part, y_part, x_part_new, y_part_new, x(i, nj), y(i, nj), x(i+1, nj), y(i+1, nj), x_cross, y_cross, icr)
		if (icr == 1) then
			St = 2
			call reflect(x_part, y_part, x_part_new, y_part_new, x_part_reflect, y_part_reflect, x(i, nj), y(i, nj), x(i+1, nj), y(i+1, nj))
			write(*,*) 'Particle reflected from boundary', St, 'and now part. pos. is', x_part_reflect, y_part_reflect
			call calc_velo_after_reflect(x(i, nj), y(i, nj), x(i+1, nj), y(i+1, nj), u_part, v_part, omega_part, u_part_new, v_part_new, &
				omega_part_new, d_part)
			return
		end if

	end do

	do j = 1, nj - 1

		call cross_edges(x_part, y_part, x_part_new, y_part_new, x(1, j), y(1, j), x(1, j+1), y(1, j+1), x_cross, y_cross, icr)
		if (icr == 1) then
			St = 3
			write(*,*) 'Particle cross boundaty', St, 'in corrdinate', x_cross, y_cross
			return
		end if
		
		call cross_edges(x_part, y_part, x_part_new, y_part_new, x(ni, j), y(ni, j), x(ni, j+1), y(ni, j+1), x_cross, y_cross, icr)
		if (icr == 1) then
			St = 4
			write(*,*) 'Particle cross boundaty', St, 'in corrdinate', x_cross, y_cross
			return
		end if

	end do

	end subroutine

!*****************************************************************************************************************************
subroutine cross_edges(x_part, y_part, x_part_new, y_part_new, x1, y1, x2, y2, x_cross, y_cross, icr)
	implicit none

	integer :: icr
	real(8) :: x_part, y_part, x_part_new, y_part_new
	real(8) :: x1, y1, x2, y2
	real(8) :: x_cross, y_cross
	real(8) :: u_a, u_b

	icr = 0

	if (abs((y_part_new - y_part)  * (x2 - x1) - (x_part_new - x_part)  * (y2 - y1)) < 1e-10) return

	u_a = ((x_part_new - x_part)  * (y1 - y_part) - (y_part_new - y_part)  * (x1 - x_part)) / &
		((y_part_new - y_part)  * (x2 - x1) - (x_part_new - x_part)  * (y2 - y1))

	u_b = ((x2 - x1)  * (y1 - y_part) - (y2 - y1)  * (x1 - x_part)) / &
		((y_part_new - y_part)  * (x2 - x1) - (x_part_new - x_part)  * (y2 - y1))

	if ((u_a >= 0 .and. u_a <= 1) .and. (u_b >= 0 .and. u_b <= 1)) then
		x_cross = x1 + u_a * (x2 - x1)
		y_cross = y1 + u_a * (y2 - y1)
		icr = 1
	end if

end subroutine

!*****************************************************************************************************************************
subroutine calc_force(nt, dt, iter, ni, nj, i_part, j_part, rho_envi, rho_part, mu, d_part, u_part, v_part, omega_part, &
	u_envi, v_envi, v_rel, rotV, gradP, F_x, F_y, M_z, &
	u_part_old, v_part_old, u_envi_old, v_envi_old, &
	F_Stokes, F_Buoyancy, F_Magnus, F_Saffman, F_Basset, F_Add_mass)
	implicit none

	real(8), parameter :: pi = 4.0 * atan(1.0)

	integer :: nt, iter, ni, nj, i_part, j_part
	real(8) :: dt, m_part, inert_moment_part, d_part, rho_part, u_part, v_part, omega_part, Re_part, Re_omega, C_L
	real(8) :: u_envi, v_envi, rho_envi, mu, v_rel(0:nt, 2), v_r_abs, omega_rel
	real(8) :: u_part_old, v_part_old, u_envi_old, v_envi_old
	real(8) :: rotV(0:ni, 0:nj), gradP(0:ni, 0:nj, 2)
	real(8) :: f, integral_basset_force(2)
	real(8) :: F_Stokes(2), F_Buoyancy(2), F_Magnus(2), F_Saffman(2), m_add, F_Basset(2), F_Add_mass(2), F_Gravity(2)
	real(8) :: F_x, F_y, M_z

	m_part = 4.0 / 3.0 * pi * (0.5 * d_part)**3.0 * rho_part
	m_add = 2. / 3. * pi * (0.5 * d_part)**3 * rho_envi ! Add mass force
	inert_moment_part = 1. / 60. * rho_part

	v_rel(iter, 1) = u_envi - u_part ! relative velocity
	v_rel(iter, 2) = v_envi - v_part
	v_r_abs = sqrt((v_rel(iter, 1))**2 + (v_rel(iter, 2))**2) ! abs of relative velocity
	omega_rel = 1.0 / 2.0 * rotV(i_part, j_part) - omega_part

	Re_part = d_part * rho_envi * v_r_abs / mu
	Re_omega = d_part**2 * abs(omega_rel) * rho_envi / mu
	C_L = 64 / Re_omega

	! Stokes force
	f = 1 + 0.179 * sqrt(Re_part) + 0.013 * Re_part
	! F_Stokes(1) = f * mu * 18 / d_part / d_part / rho_part * v_rel(1) ! сразу деленное на массу
	! F_Stokes(2) = f * mu * 18 / d_part / d_part / rho_part * v_rel(2)
	F_Stokes(1) = 3 * f * pi * d_part * mu * v_rel(iter, 1) / (m_part + m_add)
	F_Stokes(2) = 3 * f * pi * d_part * mu * v_rel(iter, 2) / (m_part + m_add)

	! Buoyancy force
	F_Buoyancy(1) = - 4.0 / 3.0 * pi * (d_part / 2.0)**3 * gradP(i_part, j_part, 1) / (m_part + m_add)
	F_Buoyancy(2) = - 4.0 / 3.0 * pi * (d_part / 2.0)**3 * gradP(i_part, j_part, 2) / (m_part + m_add)

	! Magnus force
	F_Magnus(1) = - pi * d_part**3 / 8.0 * rho_envi * v_rel(iter, 2) * omega_rel / (m_part + m_add)
	F_Magnus(2) = pi * d_part**3 / 8.0 * rho_envi * v_rel(iter, 1) * omega_rel / (m_part + m_add)

	! Saffman force
	F_Saffman(1) = sqrt(2.0) / 4.0 * 6.46 * d_part**2 * (mu * rho_envi)**0.5 * v_rel(iter, 2) * 1. / 2. * & 
		rotV(i_part, j_part) / sqrt(1. / 2. * abs(rotV(i_part, j_part))) / (m_part + m_add)
	F_Saffman(2) = - sqrt(2.0) / 4.0  * 6.46 * d_part**2 * (mu * rho_envi)**0.5 * v_rel(iter, 1) * 1. / 2. * &
		rotV(i_part, j_part) / sqrt(1. / 2. * abs(rotV(i_part, j_part))) / (m_part + m_add)

	! Basset force
	integral_basset_force = 0.0
	call calc_integral_basset_force(iter, nt, dt, v_rel, integral_basset_force)
	F_Basset(:) = 3. / 2. * d_part**2 * sqrt(pi * rho_envi * mu) * integral_basset_force(:) / (m_part + m_add)
	! F_Basset(1) = 3. / 2. * d_part**2 * sqrt(pi * rho_envi * mu) * integral_basset_force(1) / (m_part + m_add)
	! F_Basset(2) = 3. / 2. * d_part**2 * sqrt(pi * rho_envi * mu) * integral_basset_force(2) / (m_part + m_add)

	! Add mass force
	F_Add_mass(1) = - m_add * (u_part_old - u_part) / dt
	F_Add_mass(2) = - m_add * (v_rel(iter, 2) - v_rel(iter - 1, 2)) / dt

	! Gravity
	F_Gravity(1) = 0.0
	F_Gravity(2) = 0.0

	! F_Stokes(1) = 0.0
	! F_Stokes(2) = 0.0
	! F_Buoyancy(1) = 0.0
	! F_Buoyancy(2) = 0.0
	! F_Basset(1) = 0.0
	! F_Basset(2) = 0.0
	! F_Saffman(1) = 0.0
	! F_Saffman(2) = 0.0
	! F_Magnus(1) = 0.0
	! F_Magnus(2) = 0.0
	! m_add = 0.0

	! Final force
	F_x = (F_Stokes(1) + F_Buoyancy(1) + F_Magnus(1) + F_Saffman(1) + F_Basset(1) + F_Gravity(1)) ! / (m_part + m_add)
	F_y = (F_Stokes(2) + F_Buoyancy(2) + F_Magnus(2) + F_Saffman(2) + F_Basset(2) + F_Gravity(2)) ! / (m_part + m_add)
	M_z = 1. / 64. * C_L * abs(omega_rel) * omega_rel / inert_moment_part

end subroutine

!*****************************************************************************************************************************
subroutine reflect(x_part, y_part, x_part_new, y_part_new, x_part_reflect, y_part_reflect, x1, y1, x2, y2)
	implicit none

	real(8) :: x1, x2, y1, y2, x_part, y_part, x_part_new, y_part_new, x_part_reflect, y_part_reflect
	real(8) :: a, b, c

	a = y2 - y1
	b = x1 - x2
	c = - (a * x1 + b * y1)

	x_part_reflect = ((b**2 - a**2) * x_part_new - 2 * a * b * y_part_new - 2 * c * a) / (a**2 + b**2)
	y_part_reflect = (- 2 * a * b * x_part_new + (a**2 - b**2) * y_part_new - 2 * c * b) / (a**2 + b**2)

end subroutine

!*****************************************************************************************************************************
subroutine calc_velo_after_reflect(x1, y1, x2, y2, u_old, v_old, omega_old, u_new, v_new, omega_new, d_part)
	implicit none

	real(8), parameter :: e = 0.7
	real(8)			   :: x1, y1, x2, y2
	real(8) 		   :: u_old, v_old, omega_old, u_new, v_new, omega_new
	real(8)			   :: n_x, n_y, tau_x, tau_y, d
	real(8) 		   :: d_part
	real(8)			   :: U_ct_x, U_ct_y, U_c_x, U_c_y, U_x, U_y

	d = sqrt((y2 - y1)**2 + (x2 - x1)**2)

	n_x = -(y2 - y1) / d
	n_y = (x2 - x1) / d
	tau_x = (x2 - x1) / d
	tau_y = (y2 - y1) / d

! ********** simple model **********

	! u_new = (u_old * tau_x + v_old * tau_y) * tau_x - (u_old * n_x + v_old * n_y) * n_x
	! v_new = (u_old * tau_x + v_old * tau_y) * tau_y - (u_old * n_x + v_old * n_y) * n_y
	! omega_new = 0.0

! ********** improved model **********

	U_x = u_old
	U_y = v_old

	U_c_x = U_x - d_part / 2.0 * omega_old * n_y
	U_c_y = U_y + d_part / 2.0 *  omega_old * n_x

	U_ct_x = U_c_x - (U_c_x * n_x + U_c_y * n_y) * n_x
	U_ct_y = U_c_y - (U_c_x * n_x + U_c_y * n_y) * n_y

	u_new = u_old - ((1 + e) * (U_x * n_x + U_y * n_y) * n_x + 2.0 / 7.0 * sqrt(U_ct_x**2 + U_ct_y**2) * tau_x)
	v_new = v_old - ((1 + e) * (U_x * n_x + U_y * n_y) * n_y + 2.0 / 7.0 * sqrt(U_ct_x**2 + U_ct_y**2) * tau_y)
	omega_new = omega_old - 5.0 / (7.0 * d_part / 2.0) * sqrt(U_ct_x**2 + U_ct_y**2) * (n_x * tau_y - n_y * tau_x)

end subroutine

!*******************************************************************************************************************************
subroutine write_part_trac(x_part, y_part, io)
	implicit none

	integer :: io
	real(8) :: x_part, y_part

	write(io,*) x_part, y_part

end subroutine

!*******************************************************************************************************************************
subroutine write_forces(t, F_Stokes, F_Buoyancy, F_Magnus, F_Saffman, F_Basset, F_Add_mass, F_x, F_y, io, &
	i_part_old, j_part_old, i_part, j_part)
	implicit none

	integer :: io, i_part_old, j_part_old, i_part, j_part
	real(8) :: t
	real(8) :: F_Stokes(2), F_Buoyancy(2), F_Magnus(2), F_Saffman(2), F_Basset(2), F_Add_mass(2), F_x, F_y
	real(8) :: calc_abs_vector

	if (i_part == i_part_old .and. j_part == j_part_old) return

	write(io,*) t, calc_abs_vector(F_Stokes), calc_abs_vector(F_Buoyancy), calc_abs_vector(F_Magnus), &
		calc_abs_vector(F_Saffman), calc_abs_vector(F_Basset), calc_abs_vector(F_Add_mass), sqrt(F_x**2 + F_y**2)

end subroutine

!*******************************************************************************************************************************
subroutine calc_integral_basset_force(iter, nt, dt, v_rel, integral_basset_force)
	implicit none

	integer :: iter, nt, i 						! iter - number of time step
	real(8) :: dt, integral_basset_force(2)
	real(8) :: v_rel(0:nt, 2)

	integral_basset_force(1) = 0.0
	integral_basset_force(2) = 0.0
	do i = 1, iter
		integral_basset_force(1) = integral_basset_force(1) + (v_rel(i, 1) - v_rel(i - 1, 1)) / sqrt(dt * (iter - (i - 1)))
		integral_basset_force(2) = integral_basset_force(2) + (v_rel(i, 2) - v_rel(i - 1, 2)) / sqrt(dt * (iter - (i - 1)))
	end do

end subroutine

!*****************************************************************************************************************************
real(8) function triang_square(x1, y1, x2, y2, x3, y3)
	implicit none

	real(8) :: x1, y1, x2, y2, x3, y3, a, b, c, p

	a = sqrt((x1 - x2)**2 + (y1 - y2)**2)
	b = sqrt((x3 - x2)**2 + (y3 - y2)**2)
	c = sqrt((x1 - x3)**2 + (y1 - y3)**2)

	p = (a + b + c) / 2.0

	triang_square = sqrt((p - a) * (p - b) * (p - c) * p)

end function