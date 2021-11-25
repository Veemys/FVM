subroutine output_fields(io, NI, NJ, x, y, p, grad, graderror, v, div, diverror, rot, roterror, laplacian, laplacianerror)
implicit none

integer 							:: io, NI, NJ
real(8), dimension(NI, NJ) 			:: x, y
real(8), dimension(0:NI, 0:NJ) 		:: p, div, diverror, rot, roterror, laplacian, laplacianerror
real(8), dimension(0:NI, 0:NJ, 2) 	:: grad, graderror, v

write(io, *) 'VARIABLES = "X", "Y", "P", "Vx", "Vy", "GradPx", "GradPy", "DivV", "RotV", "LaplacianP", &
						  "DivV_err", "RotV_err", "LaplacianP_error", "GradPx_err", "GradPy_err"'
write(io, *) 'ZONE I=', NI, ', J=', NJ, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
write(io, '(100F14.7)') x(1:NI, 1:NJ) 
write(io, '(100F14.7)') y(1:NI, 1:NJ)
write(io, '(100F14.7)') p(1:NI-1, 1:NJ-1)
write(io, '(100F14.7)') v(1:NI-1, 1:NJ-1, 1)
write(io, '(100F14.7)') v(1:NI-1, 1:NJ-1, 2)
write(io, '(100F25.17)') grad(1:NI-1, 1:NJ-1, 1)
write(io, '(100F25.17)') grad(1:NI-1, 1:NJ-1, 2)
write(io, '(100F25.17)') div(1:NI-1, 1:NJ-1)
write(io, '(100F25.17)') rot(1:NI-1, 1:NJ-1)
write(io, '(100F25.17)') laplacian(1:NI-1, 1:NJ-1)
write(io, '(100F25.17)') diverror(1:NI-1, 1:NJ-1)
write(io, '(100F25.17)') roterror(1:NI-1, 1:NJ-1)
write(io, '(100F25.17)') laplacianerror(1:NI-1, 1:NJ-1)
write(io, '(100F25.17)') graderror(1:NI-1, 1:NJ-1, 1)
write(io, '(100F25.17)') graderror(1:NI-1, 1:NJ-1, 2)

end subroutine

! out temp field
subroutine output_temperature(io, ni, nj, x, y, t)
implicit none

integer								:: io, i, j, ni, nj
real(8), dimension(0:ni, 0:nj) 		:: x, y, t

write(io, *) 'VARIABLES = "X", "Y", "T"'
write(io, *) 'ZONE I=', ni, ', J=', nj, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
write(io, '(100F14.7)') x(1:ni, 1:nj) 
write(io, '(100F14.7)') y(1:ni, 1:nj)
write(io, '(100F14.7)') t(1:ni-1, 1:nj-1)

end subroutine