subroutine B_OutputFields(IO, NI, NJ, X, Y, P, GradP, gradPError)

real, dimension(NI, NJ) :: X, Y
real, dimension(0:NI, 0:NJ) :: P
real, dimension(0:NI, 0:NJ, 2) :: GradP, gradPError

write(IO, *) 'VARIABLES = "X", "Y", "P", "GradPx", "GradPy", "GradPx_err", "GradPy_err"'
write(IO, *) 'ZONE I=',NI, ', J=',NJ, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
write(IO, '(100F14.7)') X(1:NI, 1:NJ) 
write(IO, '(100F14.7)') Y(1:NI, 1:NJ)
write(IO, '(100F14.7)') P(1:NI-1, 1:NJ-1)
write(IO, '(100F25.17)') GradP(1:NI-1, 1:NJ-1, 1)
write(IO, '(100F25.17)') GradP(1:NI-1, 1:NJ-1, 2)
write(IO, '(100F25.17)') gradPError(1:NI-1, 1:NJ-1, 1)
write(IO, '(100F25.17)') gradPError(1:NI-1, 1:NJ-1, 2)

end Subroutine 
