subroutine output_fields(io, NI, NJ, x, y, p, grad, graderror, v, div, diverror, rot, roterror)
implicit none

integer 							:: io, NI, NJ
real(8), dimension(NI, NJ) 			:: x, y
real(8), dimension(0:NI, 0:NJ) 		:: p, div, diverror, rot, roterror
real(8), dimension(0:NI, 0:NJ, 2) 	:: grad, graderror, v

write(io, *) 'VARIABLES = "X", "Y", "P", "Vx", "Vy", "GradPx", "GradPy", "DivV", "RotV", &
						  "DivV_err", "RotV_err", "GradPx_err", "GradPy_err"'
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
write(io, '(100F25.17)') diverror(1:NI-1, 1:NJ-1)
write(io, '(100F25.17)') roterror(1:NI-1, 1:NJ-1)
write(io, '(100F25.17)') graderror(1:NI-1, 1:NJ-1, 1)
write(io, '(100F25.17)') graderror(1:NI-1, 1:NJ-1, 2)

end Subroutine