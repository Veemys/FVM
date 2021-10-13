real function Pressure(X, Y)
implicit none
real :: x, y

Pressure = x**2 + y

end function

real function LinearInterpol(d1, d2, x1, x2)
implicit none
real :: x1, x2, d1, d2

LinearInterpol = (x1 * d2 + x2 * d1) / (d1 + d2)

end function

subroutine CalcGradPExact(ni, nj, x, y, gradPExact)
implicit none
integer :: i, j, ni, nj
real :: x, y
real :: gradPExact(0:ni, 0:nj, 2)

gradPExact(i, j, 1) = 2.0 * x
gradPExact(i, j, 2) = 1.0

end subroutine