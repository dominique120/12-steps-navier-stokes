module functions
    contains
    
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
    
subroutine diffuse(timesteps, u, nu, dt, dx, dy, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny, timesteps
    real(kind=8), intent(in) :: nu, dt, dx, dy
    real(kind=8), intent(out), dimension(nx,ny) :: u
    real(kind=8), dimension(nx,ny) :: un
    integer :: i, j, time
    real(kind=8) :: Ydisip, Xdisip, Yderiv, Xderiv
    time = timesteps

    do time = 1, time
            un = u
        
            do j=2, nx-1
                do i=2, ny-1
                    Xderiv = (nu * dt) / (dx**2)
                    Yderiv = (nu * dt) / (dy**2)
                    Xdisip = un(j,i+1) - (2 * un(j,i)) + un(j,i-1)
                    Ydisip = un(j+1,i) - (2 * un(j,i)) + un(j-1,i)
                    
                    u(j,i) = un(j,i) + (Xderiv * Xdisip) + (Yderiv * Ydisip)

                    u(1,j) = 1
                    u(nx,j) = 1
                    u(i,1) = 1
                    u(i,ny) = 1
        
                    u(1,1) = 1
                    u(1,nx) = 1
                    u(ny,1) = 1
                    u(ny,nx) = 1
                end do
            end do
        end do
end subroutine 

end module functions