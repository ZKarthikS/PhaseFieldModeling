module testfun 
    implicit none
contains

real function f(x)
    real, intent(in) :: x
    f = x**2
end function f

subroutine calc_(x, y, f_)
    real :: x, f_(1:10), y(1:10)
    integer :: i
    do i = 1,10
        f_(i) = y(i)/2
    enddo
    return
end subroutine calc_


end module testfun