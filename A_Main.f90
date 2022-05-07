program main
implicit none

! variables declarations
integer 									:: i, j, k, ni, nj, gradCorrectIter, maxiter_correct, gradient_scheme, div_mode, mode, maxiter
integer, parameter 							:: io = 12
character(*), parameter 					:: input_file = 'input.txt', input_particle_file = 'input_particle.txt', output_file = 'data.plt', &
											output_file_temp = 'Tfield.plt'
character 									:: mesh_file * 30, inputdata * 30, ctmp
real(8) 									:: pressure, velocity, calcgrad_exact, calcdiv_exact, calcrot_exact, calclaplacian_exact, rtmp, Re, Pr, eps
real(8), allocatable, dimension(:, :) 		:: x, y, p, t, cellvolume, div, divexact, diverror, rot, rotexact, roterror, &
										  	laplacian, laplacianexact, laplacianerror
real(8), allocatable, dimension(:, :, :) 	:: v, grad, gradexact, graderror
real(8), allocatable, dimension(:, :, :) 	:: cellcenter, iface_center, iface_vector, jface_center, jface_vector

! particles information
integer :: nt, time_iter, i_part, j_part
real(8) :: rho_envi, rho_part, mu, d_part, x_0_part, y_0_part, u_0_part, v_0_part, omega_0_part, dt, &
		x_part, y_part, u_part, v_part, omega_part, F_x, F_y, M_z, Sk, u_envi, v_envi
real(8) :: v_ref, l_ref
		
! read input file
write(*, *) 'Read input file: ', input_file
open(io, file = input_file)
read(io, *) mesh_file  		! read name of file with computational mesh
read(io, *) inputdata		! read name of file with FloS-input fields
read(io, *) maxiter_correct ! number of iterations for gradient correction
read(io, *) gradient_scheme	! scheme for gradient computation
read(io, *) div_mode		! scheme for divergency computation
read(io, *) mode			! way of calc. of oprators
read(io, *) Re				! Reynolds number
read(io, *) Pr				! Prandtl number
read(io, *) maxiter			! maxiter for solving conv-dif eq
read(io, *) eps				! accurance for solving conv-dif eq
close(io)

! read nodes number (ni, nj) from file with mesh
write(*, *) 'Read nodes number from file: ', mesh_file
open(io, file = mesh_file)
read(io, *) ni, nj
write(*, *) 'NI, NJ = ', ni, nj

! allocate all arrays
write(*, *) 'Allocate arrays...'       
allocate(x(NI, NJ)) 					! mesh nodes X-coordinates
allocate(y(NI, NJ))						! mesh nodes Y-coordinates
allocate(p(0:NI, 0:NJ))   				! Pressure
allocate(t(0:NI, 0:NJ))					! Temperature
allocate(v(0:NI, 0:NJ, 2))				! Velocity
allocate(div(0:NI, 0:NJ))				! div(V)
allocate(divexact(0:NI, 0:NJ))			! Exact value of divergency
allocate(diverror(0:NI, 0:NJ))			! Values of the divergency error
allocate(rot(0:NI, 0:NJ))				! rot(V)
allocate(rotexact(0:NI, 0:NJ))			! Exact value of rotor
allocate(roterror(0:NI, 0:NJ))			! Values of the rotor error
allocate(laplacian(0:NI, 0:NJ))			! div(grad(p))
allocate(laplacianexact(0:NI, 0:NJ))	! Exact value of laplacian
allocate(laplacianerror(0:NI, 0:NJ))	! Values of the laplacian error
allocate(grad(0:NI, 0:NJ, 2)) 			! grad(p)
allocate(gradexact(0:NI, 0:NJ, 2))		! Exact value of gradient 
allocate(graderror(0:NI, 0:NJ, 2))		! Values of the gradient error
allocate(cellvolume(NI-1, NJ-1))   		! Cell Volumes    
allocate(cellcenter(0:NI, 0:NJ, 2)) 	! Cell Centers
allocate(iface_center(NI, NJ-1, 2))		! Face Centers for I-faces
allocate(iface_vector(NI, NJ-1, 2)) 	! Face Vectors for I-faces
allocate(jface_center(NI-1, NJ, 2)) 	! Face Centers for J-faces
allocate(jface_vector(NI-1, NJ, 2)) 	! Face Vectors for J-faces

! read grid, calculate metric, initiate fields
write(*, *) 'Read mesh from file: ', mesh_file

select case(mode)

	case(0)

		read(io, *) ((x(i, j), y(i, j), i = 1, ni), j = 1, nj)
		close(io)

		write(*, *) 'Calculate metric...'       
		call calc_metric(ni, nj, x, y, cellcenter, cellvolume, iface_center, iface_vector, jface_center, jface_vector)
  
		write(*, *) 'Initiate fields...'
		do j = 0, nj
			do  i = 0, ni
	
				p(i, j) = pressure(cellcenter(I, J, 1), cellcenter(I, J, 2))
				v(i, j, 1) = velocity(cellcenter(I, J, 1), cellcenter(I, J, 2), 1)
				v(i, j, 2) = velocity(cellcenter(I, J, 1), cellcenter(I, J, 2), 2)
				gradexact(i, j, 1) = calcgrad_exact(cellcenter(I, J, 1), cellcenter(I, J, 2), 1)
				gradexact(i, j, 2) = calcgrad_exact(cellcenter(I, J, 1), cellcenter(I, J, 2), 2)
				divexact(i, j) = calcdiv_exact(cellcenter(i, j, 1), cellcenter(i, j, 2), p(i, j), v(i, j, :), gradexact(i, j, :), div_mode)
				rotexact(i, j) = calcrot_exact(cellcenter(i, j, 1), cellcenter(i, j, 2))
				laplacianexact(i, j) = calclaplacian_exact(cellcenter(i, j, 1), cellcenter(i, j, 2))
		
			end do
		end do
	
	case(1)
		
		read(io, *) ((x(i, j), y(i, j), rtmp, i = 1, ni), j = 1, nj)
		close(io)
		
		write(*, *) 'Calculate metric...'       
		call calc_metric(ni, nj, x, y, cellcenter, cellvolume, iface_center, iface_vector, jface_center, jface_vector)
		
		write(*, *) 'Initiate fields...'
		open(io, file = inputdata)
		read(io, *) ctmp
		read(io, *) ctmp
		read(io, *) ((rtmp, rtmp, v(i, j, 1), v(i, j, 2), rtmp, p(i, j), rtmp, rtmp, i = 0, ni), j = 0, nj)
		close(io)

end select

! reading particle's information
open(io, file = input_particle_file)
read(io, *) rho_envi, rho_part 	! density enviroment and particle
read(io, *) mu 					! enviroment viscosity
read(io, *) d_part 				! particle's diameter
read(io, *) x_0_part, y_0_part 	! initial position of particle
read(io, *) u_0_part, v_0_part 	! initial velocity of particle
read(io, *) omega_0_part 		! initioal circular velocity of particle
read(io, *) dt 					! time step
read(io, *) nt 					! number of time step
close(io)


! calculate gradient
write(*, *) 'Calculate gradient...'
call calcgrad_choice_scheme(ni, nj, gradient_scheme, maxiter_correct, p, grad, cellvolume, cellcenter, &
							iface_center, iface_vector, jface_center, jface_vector)
write(*, *) 'COMPLETE!'

! calculate divergency
write(*, *) 'Calculate divergency...'
call calcdiv(NI, NJ, v, p, div, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, div_mode)
write(*, *) 'COMPLETE!'

! calculate rotor
write(*, *) 'Calculate rotor...'
call calcrot(NI, NJ, v, rot, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
write(*, *) 'COMPLETE!'

! calculate laplacian
write(*, *) 'Calculate laplacian...'
call calclaplacian(NI, NJ, p, grad, laplacian, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
write(*, *) 'COMPLETE!'


write(*, *) 'Finding particle traectory...'
call solver_part(rho_part, d_part, x_0_part, y_0_part, u_0_part, v_0_part, omega_0_part, rho_envi, x, y, v, &
				cellvolume, iface_vector, jface_vector, dt, nt, ni, nj, mu, &
				x_part, y_part, u_part, v_part, omega_part, i_part, j_part, &
				rot, grad)

! calculate errors
! call calcgrad_error(NI, NJ, grad, gradexact, graderror)
! call calcdiv_error(NI, NJ, div, divexact, diverror)
! call calcrot_error(NI, NJ, rot, rotexact, roterror)
! call calclaplacian_error(NI, NJ, laplacian, laplacianexact, laplacianerror)

! ! output field
! write(*, *) 'Output fields to file: ', output_file
! open(io, file = output_file)
! call output_fields(io, NI, NJ, x, y, p, grad, graderror, v, div, diverror, rot, roterror, laplacian, laplacianerror)
! close(io)

! solving equation
! if (mode == 1) then
	
! 	write(*, *) 'Calculate temperature field'
! 	call Solver(ni, nj, v, t, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, &
! 				Re, Pr, eps, maxiter, div_mode, gradient_scheme, maxiter_correct)

! 	write(*, *) 'Output field of temperature to file: ', output_file_temp
! 	open(io, file = output_file_temp)
! 	call output_temperature(io, ni, nj, x, y, t)
! 	close(io)
	
! end if

end program main