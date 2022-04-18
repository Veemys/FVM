!!*****************************************************************************************************************************
subroutine solver_part(rho_part, d_part, x_0_part, y_0_part, u_0_part, v_0_part, omega_0_part, rho_envi, x, y, v, &
                       cellvolume, iface_vector, jface_vector, dt, nt, ni, nj, mu, &
					   x_part, y_part, u_part, v_part, omega_part, i_part, j_part)
implicit none

integer :: iter, nt, ni, nj, i_part, j_part, ip1, jp1, St
real(8) :: x_part, y_part, u_part, v_part, omega_part, x_0_part, y_0_part, u_0_part, v_0_part, omega_0_part, F_x, F_y, M_z, &
           v_ref, l_ref, Sk, rho_part, d_part, rho_envi, mu, u_envi, v_envi, dt
real(8) :: iface_vector(ni, nj-1, 2), jface_vector(ni-1, nj, 2)
real(8) :: x_part_new, y_part_new, u_part_new, v_part_new, omega_part_new
real(8) :: v(0:ni, 0:nj, 2), x(ni, nj), y(ni, nj), cellvolume(ni-1, nj-1)
real(8) :: x_cross, y_cross, x_part_reflect, y_part_reflect

v_ref = 1.0
l_ref = 1.0

Sk = rho_part * d_part**2 * v_ref / (mu * l_ref)
print*, 'Stokes number = ', Sk

open(69, file = 'part_trac.plt')
write(69,*) 'VARIABLES = "X", "Y"'

x_part = x_0_part
y_part = y_0_part
u_part = u_0_part
v_part = v_0_part
omega_part = omega_0_part
write(69,*) x_part, y_part

! v = v_ref

do iter = 1, nt
	
	! calc index of cell with particle
	i_part = - 1
	j_part = - 1
	call find_part_location(x_part, y_part, ni, nj, x, y, cellvolume, i_part, j_part)
	
	u_envi = v(i_part, j_part, 1)
	v_envi = v(i_part, j_part, 2)
	
	! F_x = 0.0
	! F_y = 0.0
	! M_z = 0.0
	
	call calc_force(rho_envi, rho_part, mu, d_part, u_part, v_part, omega_part, u_envi, v_envi, F_x, F_y, M_z)
	
	x_part_new = x_part + dt * u_part
	y_part_new = y_part + dt * v_part
	u_part_new = u_part + dt * F_x
	v_part_new = v_part + dt * F_y
	omega_part_new = omega_part + dt * M_z

	! Check boundary
	
	call find_part_location(x_part_new, y_part_new, ni, nj, x, y, cellvolume, i_part, j_part)
	
	if (i_part == - 1 .and. j_part == - 1) then
		call C_Boundary(x, y, iface_vector, jface_vector, ni, nj, x_part, y_part, u_part, v_part, omega_part, &
						x_part_new, y_part_new, u_part_new, v_part_new, omega_part_new, x_cross, y_cross, ip1, jp1, St, &
						x_part_reflect, y_part_reflect)

		if (St == 1 .or. St == 2) then
			write(69,*) x_cross, y_cross
			!write(69,*) x_part_new, y_part_new
			!write(69,*) x_part_reflect, y_part_reflect
			x_part_new = x_part_reflect
			y_part_new = y_part_reflect
		end if

		if (St == 3 .or. St == 4) then
			write(69,*) x_cross, y_cross
			exit
		end if

	end if
	
	x_part = x_part_new
	y_part = y_part_new
	u_part = u_part_new
	v_part = v_part_new
	omega_part = omega_part_new
	write(69,*) x_part, y_part
	
	print*, 'time =', iter * dt, 'i =', i_part, 'j =', j_part, 'x =', x_part, 'y =', y_part
	
end do
close(69)

end subroutine

!!*****************************************************************************************************************************
subroutine find_part_location(x_part, y_part, ni, nj, x, y, cellvolume, i_part, j_part)
implicit none

integer :: ni, nj, i_part, j_part, i, j
real(8) :: x_part, y_part, x(ni, nj), y(ni, nj), cellvolume(ni-1, nj-1), s, eps, triang_square

i_part = - 1
j_part = - 1
eps = 1e-10

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

!!*****************************************************************************************************************************
real(8) function triang_square(x1, y1, x2, y2, x3, y3)
implicit none

real(8) :: x1, y1, x2, y2, x3, y3, a, b, c, p

a = sqrt((x1 - x2)**2 + (y1 - y2)**2)
b = sqrt((x3 - x2)**2 + (y3 - y2)**2)
c = sqrt((x1 - x3)**2 + (y1 - y3)**2)

p = (a + b + c) / 2.0

triang_square = sqrt((p - a) * (p - b) * (p - c) * p)

end function

!!*****************************************************************************************************************************
subroutine C_Boundary(x, y, iface_vector, jface_vector, ni, nj, x_part, y_part, u_part, v_part, omega_part, &
						x_part_new, y_part_new, u_part_new, v_part_new, omega_part_new, x_cross, y_cross, ip1, jp1, St, &
						x_part_reflect, y_part_reflect)
implicit none

integer :: ni, nj, ip1, jp1, St, i, j, icr ! icr - index + cross
real(8) :: x(ni, nj), y(ni, nj)
real(8) :: iface_vector(ni, nj-1, 2), jface_vector(ni-1, nj, 2)
real(8) :: x_part_new, y_part_new, u_part_new, v_part_new, omega_part_new, x_part, y_part, u_part, v_part, omega_part
real(8) :: x_cross, y_cross, x_part_reflect, y_part_reflect

St = 0
icr = 0

do i = 1, ni - 1

	call cross_edges(x_part, y_part, x_part_new, y_part_new, x(i, 1), y(i, 1), x(i+1, 1), y(i+1, 1), x_cross, y_cross, icr)
	if (icr == 1) then
		St = 1
		call reflect(x_part, y_part, x_part_new, y_part_new, x_part_reflect, y_part_reflect, x(i, 1), y(i, 1), x(i+1, 1), y(i+1, 1))
		!write(*,*) 'Particle cross boundaty', St, 'in corrdinate', x_cross, y_cross
		write(*,*) 'Particle reflected from boundary', St, 'and now part. pos. is', x_part_reflect, y_part_reflect
		write(*,*) u_part_new, v_part_new
		call calc_velo_after_reflect(x(i, 1), y(i, 1), x(i+1, 1), y(i+1, 1), u_part, v_part, omega_part, u_part_new, v_part_new, &
			omega_part_new)
		write(*,*) u_part_new, v_part_new
		!x_cross = x_part_reflect
		!y_cross = y_part_reflect
		return
	end if
	
	call cross_edges(x_part, y_part, x_part_new, y_part_new, x(i, nj), y(i, nj), x(i+1, nj), y(i+1, nj), x_cross, y_cross, icr)
	if (icr == 1) then
		St = 2
		call reflect(x_part, y_part, x_part_new, y_part_new, x_part_reflect, y_part_reflect, x(i, nj), y(i, nj), x(i+1, nj), y(i+1, nj))
		! write(*,*) 'Particle cross boundaty', St, 'in corrdinate', x_cross, y_cross
		write(*,*) 'Particle reflected from boundary', St, 'and now part. pos. is', x_part_reflect, y_part_reflect
		call calc_velo_after_reflect(x(i, nj), y(i, nj), x(i+1, nj), y(i+1, nj), u_part, v_part, omega_part, u_part_new, v_part_new, &
			omega_part_new)
		!x_cross = x_part_reflect
		!y_cross = y_part_reflect
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

!!*****************************************************************************************************************************
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

! write(*,*) u_a, u_b

if ((u_a >= 0 .and. u_a <= 1) .and. (u_b >= 0 .and. u_b <= 1)) then
	x_cross = x1 + u_a * (x2 - x1)
	y_cross = y1 + u_a * (y2 - y1)
	icr = 1
end if

end subroutine

!!*****************************************************************************************************************************
subroutine calc_force(rho_envi, rho_part, mu, d_part, u_part, v_part, omega_part, u_envi, v_envi, F_x, F_y, M_z)
implicit none

real(8), parameter :: pi = 4.0 * atan(1.0)
real(8) :: m_part, d_part, rho_part, u_part, v_part, omega_part, Re_part
real(8) :: u_envi, v_envi, rho_envi, mu, v_r
real(8) :: f, F_Stokes(2), F_Gravity(2), F_x, F_y, M_z

m_part = 4.0 / 3.0 * pi * (0.5 * d_part)**3.0 * rho_part

v_r = sqrt((u_part - u_envi)**2 + (v_part - v_envi)**2) ! abs of relative vrlocity
Re_part = d_part * rho_envi * v_r / mu

! Stokes force
f = 1 + 0.179 * sqrt(Re_part) + 0.013 * Re_part
F_Stokes(1) = f * mu * 18 / d_part / d_part / rho_part * (u_envi - u_part)
F_Stokes(2) = f * mu * 18 / d_part / d_part / rho_part * (v_envi - v_part)

! Gravity
F_Gravity(1) = 0.0
F_Gravity(2) = 0.0

! Final force
F_x = F_Stokes(1) + F_Gravity(1)
F_y = F_Stokes(2) + F_Gravity(2)
M_z = 0.0

end subroutine

!!*****************************************************************************************************************************
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

!!*****************************************************************************************************************************
subroutine calc_velo_after_reflect(x1, y1, x2, y2, u_old, v_old, omega_old, u_new, v_new, omega_new)
implicit none

real(8) :: x1, y1, x2, y2
real(8) :: u_old, v_old, omega_old, u_new, v_new, omega_new
real(8) :: n_x, n_y, tau_x, tau_y, d

d = sqrt((y2 - y1)**2 + (x2 - x1)**2)

n_x = -(y2 - y1) / d
n_y = (x2 - x1) / d
tau_x = (x2 - x1) / d
tau_y = (y2 - y1) / d

u_new = (u_old * tau_x + v_old * tau_y) * tau_x + (u_old * n_x + v_old * n_y) * n_x
v_new = - ((u_old * tau_x + v_old * tau_y) * tau_y + (u_old * n_x + v_old * n_y) * n_y)
omega_new = 0.0

end subroutine