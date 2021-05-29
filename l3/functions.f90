real(kind=8) function burger (t, x ,nu) result(burg)
    implicit none
    
    real(kind=8), intent(in) :: t, x, nu
    real(kind=8) :: pi = 4 * atan(1.0_8)
    
    burg = -2 * nu * (-(-8 * t + 2 * x) * exp(-(-4 * t + x)** 2 / &
		(4 * nu * (t + 1))) / (4 * nu * (t + 1)) - &
		(-8 * t + 2 * x - 4 * pi) * exp(-(-4 * t + x - 2 * pi)** 2 / &
		(4 * nu * (t + 1))) / (4 * nu * (t + 1))) / &
			(exp(-(-4 * t + x - 2 * pi)** 2 / (4 * nu * (t + 1))) + &
				exp(-(-4 * t + x)** 2 / (4 * nu * (t + 1)))) + 4

end function burger

subroutine spaceevenly(start, finish, number, arrayout)
    implicit none
    
    integer, intent(in) :: start, finish
    real(kind=8), intent(in) :: number
    real(kind=8) :: div
    real(kind=8), dimension(0:finish), intent(out) :: arrayout
    integer :: i
    
    div = number / ((finish-1) - (start-1))
    
    do i=start, finish
        arrayout(i) = div * i
    end do

end subroutine spaceevenly
    
!subroutine copy_real_array(source, destination, size)
    !implicit none
    