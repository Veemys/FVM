program main
implicit none

integer 								:: i, j, k, NI, NJ, dradCorrectIter, maxiter_correct, gradient_scheme, div_mode, mode, maxiter
integer, parameter 						:: io = 12 												! input-output unit
character(*), parameter 				:: inputfile = 'input.txt', outputfile = 'data.plt', outputfile_temp = 'Tfield.plt'		! names of input and output files
character 								:: meshfile * 30, inputdata * 30, ctmp   				! name of file with computational mesh
real(8) 								:: pressure, velocity, calcgrad_exact, calcdiv_exact, calcrot_exact, calclaplacian_exact, rtmp, Re, Pr, eps
real(8), allocatable, dimension(:,:) 	:: x, y, p, t, cellvolume, div, divexact, diverror, rot, rotexact, roterror, &
										   laplacian, laplacianexact, laplacianerror
real(8), allocatable, dimension(:,:,:) 	:: v, grad, gradexact, graderror
real(8), allocatable, dimension(:,:,:) 	:: cellcenter, iface_center, iface_vector, jface_center, jface_vector

!===  READ INPUT FILE ===
write(*,*) 'Read input file: ', inputfile
open(io, file = inputfile)
read(io, *) meshfile  		! read name of file with computational mesh
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

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
write(*,*) 'Read nodes number from file: ', meshfile
open(io, file = meshfile)
read(io, *) ni, nj
write(*,*) 'NI, NJ = ', ni, nj

!=== ALLOCATE ALL ARRAYS ===
write(*,*) 'Allocate arrays'       
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

!===  READ GRID, CALCULATE METRIC, INITIATE FIELDS ===
write(*,*) 'Read mesh from file: ', meshfile

select case(mode)

	case(0)

		read(io, *) ((x(i, j), y(i, j), i = 1, ni), j = 1, nj)
		close(io)

		write(*,*) 'Calculate metric'       
		call calcmetric(ni, nj, x, y, cellcenter, cellvolume, iface_center, iface_vector, jface_center, jface_vector)
  
		write(*,*) 'Initiate fields'
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
		
		write(*,*) 'Calculate metric'       
		call calcmetric(ni, nj, x, y, cellcenter, cellvolume, iface_center, iface_vector, jface_center, jface_vector)
		
		write(*,*) 'Initiate fields'
		open(io, file = inputdata)
		read(io, *) ctmp
		read(io, *) ctmp
		read(io, *) ((rtmp, rtmp, v(i, j, 1), v(i, j, 2), rtmp, p(i, j), rtmp, rtmp, i = 0, ni), j = 0, nj)
		close(io)

end select

!=== CALCULATE GRADIENT ===
write(*,*) 'Calculate gradient'
call calcgrad_choice_scheme(ni, nj, gradient_scheme, maxiter_correct, p, grad, cellvolume, cellcenter, &
							iface_center, iface_vector, jface_center, jface_vector)

!=== Calculate divergency ===
write(*,*) 'Calculate divergency'
call calcdiv(NI, NJ, v, p, div, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, div_mode)

!=== Calculate rotor ===
write(*,*) 'Calculate rotor'
call calcrot(NI, NJ, v, rot, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)

write(*,*) 'Calculate laplacian'
call calclaplacian(NI, NJ, p, grad, laplacian, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)

!=== Calculate errors ===
call calcgrad_error(NI, NJ, grad, gradexact, graderror)
call calcdiv_error(NI, NJ, div, divexact, diverror)
call calcrot_error(NI, NJ, rot, rotexact, roterror)
call calclaplacian_error(NI, NJ, laplacian, laplacianexact, laplacianerror)

!=== OUTPUT FIELDS ===
write(*,*) 'Output fields to file: ', outputfile
open(io, file = outputfile)
call output_fields(io, NI, NJ, x, y, p, grad, graderror, v, div, diverror, rot, roterror, laplacian, laplacianerror)
close(io)

!=== SOLVING EQUATION ===
if (mode == 1) then
	
	write(*,*) 'Calculate temperature field'
	call Solver(ni, nj, v, t, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector, &
				Re, Pr, eps, maxiter, mode, gradient_scheme, maxiter_correct)

	write(*,*) 'Output field of temperature to file: ', outputfile_temp
	open(io, file = outputfile_temp)
	call output_temperature(io, ni, nj, x, y, t)
	close(io)
	
end if

end program main