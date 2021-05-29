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


subroutine laplace2d(p, y, dx, dy, l1nom_target, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(kind=8), intent(in) :: dx, dy, l1nom_target
    real(kind=8), intent(in), dimension(nx) :: y
    real(kind=8), intent(out), dimension(nx, ny) :: p
    real(kind=8), dimension(nx, ny) :: pn
    integer :: i, j, l1norm = 1, whilecounter =1
    real(kind=8) :: derivYpow2, derivXpow2, F1term, S2term, numerator, denominator
    real(kind=8) :: lnorm_last, psum, pnsum
    
    derivYpow2 = dy**2
	derivXpow2 = dx**2
    
    ! Counter still not working with the l1norm condition
    do while (whilecounter < 708)! (l1norm > l1nom_target)
        pn = p
        do i = 2, nx-1
            do j = 2, ny-1
                F1term = p(j, i + 1) + p(j, i - 1)
                S2term = p(j + 1, i) + p(j - 1, i)
                numerator = (derivYpow2 * F1term) + (derivXpow2 * S2term)
                denominator = 2 * (derivXpow2 + derivYpow2)
                    
                p(j,i) = numerator / denominator
                
                !print*, "in main while loop -> iteration:)", i, " - ", j
            end do       
        end do
   
        ! Reset Boundary conditions
        call boundary_conditions(p, y, nx, ny)
        
        ! calculate l1norm
        do i = 1, nx
            do j = 1, ny
                psum = psum + abs(p(j,i))
				pnsum = pnsum + abs(pn(j,i))
            end do
        end do
        lnorm_last = l1norm
		l1norm = (psum - pnsum) / pnsum
    ! print*, "loop counter: ", whilecounter  
    whilecounter = whilecounter +1
    end do
    
end subroutine

subroutine boundary_conditions(p, y, nx,ny)
    implicit none
    
    real(kind=8), intent(out), dimension(nx, ny) :: p
    real(kind=8), intent(in), dimension(nx) :: y
    integer, intent(in) :: nx, ny
    integer :: i
    
    do i=1, ny
        p(i,1) = 0
    end do
    
    do i=1, ny
        p(i,ny) = y(i)
    end do
    
    do i=1, nx
        p(1,i) = p(2,i)
    end do
    
    do i=1, nx
        p(nx,i) = p(nx-1,i)
    end do

end subroutine

    
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

subroutine write_csv (unit_number, filename, p, nx, ny) 
    implicit none
    integer,  intent(in) :: unit_number
    character(len = *),  intent(in) :: filename
    integer, intent(in) :: nx, ny
    real(kind=8), intent(in), dimension(nx,ny) :: p
    integer :: i, j

10  format (i1)
15  format (f10.5,",")
20  format ("")
    
    open(unit=unit_number, file=filename)
    do i=1, nx
        do j=1, ny
            write(unit=unit_number, fmt=15, advance="no"), p(i,j)
        end do
        write(unit=unit_number,fmt=20)
    end do
    
    close(unit=unit_number)

end subroutine

end module functions