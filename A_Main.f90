program main
implicit none

integer 								:: i, j, k, NI, NJ, dradCorrectIter, maxiter_correct, gradient_scheme
integer, parameter 						:: io = 12 												! input-output unit
character(*), parameter 				:: inputfile = 'input.txt', outputfile = 'data.plt' 	! names of input and output files
character 								:: meshfile * 30   										! name of file with computational mesh
real(8) 								:: pressure, velocity, calcgrad_exact, calcdiv_exact, calcrot_exact
real(8), allocatable, dimension(:,:) 	:: x, y, p, cellvolume, div, divexact, diverror, rot, rotexact, roterror
real(8), allocatable, dimension(:,:,:) 	:: grad, v, gradexact, graderror
real(8), allocatable, dimension(:,:,:) 	:: cellcenter, iface_center, iface_vector, jface_center, jface_vector

!===  READ INPUT FILE ===
write(*,*) 'Read input file: ', inputfile
open(io, file = inputfile)
read(io, *) meshfile  		! read name of file with computational mesh
read(io, *) maxiter_correct ! number of iterations for gradient correction
read(io, *) gradient_scheme	! scheme for gradient computation
close(io)

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
write(*,*) 'Read nodes number from file: ', meshfile
open(io, file = meshfile)
read(io, *) NI, NJ
write(*,*) 'NI, NJ = ', NI, NJ

!=== ALLOCATE ALL ARRAYS ===
write(*,*) 'Allocate arrays'       
allocate(x(NI, NJ)) 				! mesh nodes X-coordinates
allocate(y(NI, NJ))					! mesh nodes Y-coordinates
allocate(p(0:NI, 0:NJ))   			! Pressure
allocate(v(0:NI, 0:NJ, 2))			! Velocity
allocate(div(0:NI, 0:NJ))			! div(V)
allocate(divexact(0:NI, 0:NJ))		! Exact value of divergency
allocate(diverror(0:NI, 0:NJ))		! Values of the divergency error
allocate(rot(0:NI, 0:NJ))			! rot(V)
allocate(rotexact(0:NI, 0:NJ))		! Exact value of rotor
allocate(roterror(0:NI, 0:NJ))		! Values of the rotor error
allocate(grad(0:NI, 0:NJ, 2)) 		! grad(p)
allocate(gradexact(0:NI, 0:NJ, 2))	! Exact value of gradient 
allocate(graderror(0:NI, 0:NJ, 2))	! Values of the gradient error
allocate(cellvolume(NI-1, NJ-1))   	! Cell Volumes    
allocate(cellcenter(0:NI, 0:NJ, 2)) ! Cell Centers
allocate(iface_center(NI, NJ-1, 2))	! Face Centers for I-faces
allocate(iface_vector(NI, NJ-1, 2)) ! Face Vectors for I-faces
allocate(jface_center(NI-1, NJ, 2)) ! Face Centers for J-faces
allocate(jface_vector(NI-1, NJ, 2)) ! Face Vectors for J-faces

!===  READ GRID ===
write(*,*) 'Read mesh from file: ', meshfile
read(io, *) ((x(i, j), y(i, j), i = 1, ni), j = 1, nj)
close(io)

!=== CALCULATE METRIC ===
write(*,*) 'Calculate metric'       
call calcmetric(NI, NJ, x, y, cellcenter, cellvolume, iface_center, iface_vector, jface_center, jface_vector) 
  
!=== INITIATE FIELDS ===
write(*,*) 'Initiate fields'       
do j = 0, NJ
	do  i = 0, NI
	
		p(i, j) = pressure(cellcenter(I, J, 1), cellcenter(I, J, 2))
		v(i, j, 1) = velocity(cellcenter(I, J, 1), cellcenter(I, J, 2), 1)
		v(i, j, 2) = velocity(cellcenter(I, J, 1), cellcenter(I, J, 2), 2)
		gradexact(i, j, 1) = calcgrad_exact(cellcenter(I, J, 1), cellcenter(I, J, 2), 1)
		gradexact(i, j, 2) = calcgrad_exact(cellcenter(I, J, 1), cellcenter(I, J, 2), 2)
		divexact(i, j) = calcdiv_exact(cellcenter(i, j, 1), cellcenter(i, j, 2))
		rotexact(i, j) = calcrot_exact(cellcenter(i, j, 1), cellcenter(i, j, 2))
		
	end do
end do

!=== CALCULATE GRADIENT ===
write(*,*) 'Calculate gradient'
if (gradient_scheme == 1) then
	do k = 1, maxiter_correct
		call calcgrad_greengauss(NI, NJ, p, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
	end do
else
	call calcgrad_leastsquare(NI, NJ, p, grad, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)
end if

!=== Calculate divergency ===
write(*,*) 'Calculate divergency'
call calcdiv(NI, NJ, v, div, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)

!=== Calculate rotor ===
write(*,*) 'Calculate rotor'
call calcrot(NI, NJ, v, rot, cellvolume, cellcenter, iface_center, jface_center, iface_vector, jface_vector)

!=== Calculate errors ===
call calcgrad_error(NI, NJ, grad, gradexact, graderror)
call calcdiv_error(NI, NJ, div, divexact, diverror)
call calcrot_error(NI, NJ, rot, rotexact, roterror)

!=== OUTPUT FIELDS ===
write(*,*) 'Output fields to file: ', outputfile
open(io, file = outputfile)
call output_fields(io, NI, NJ, x, y, p, grad, graderror, v, div, diverror, rot, roterror)
close(io)

end program main