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
    
program l10
    use pyplot_module, only : pyplot
 
    implicit none
    
    type(pyplot) :: plt

    integer, parameter :: nx = 50, ny = 50, nt = 100
    integer, parameter :: xmin = 0, xmax = 2, ymin = 0, ymax = 1
    real(kind=8) :: dx = real(xmax - xmin) / real(nx - 1)
    real(kind=8) :: dy = real(ymax - ymin) / real(ny - 1)
    integer :: i, j, time, istat
    real(kind=8),dimension(nx) :: y, x
    real(kind=8), dimension(nx, ny) :: p, pd, b
    real(kind=8) :: F1term, S2term, numerator, denominator
    real(kind=8) :: source_term, derivYpow2, derivXpow2
    character:: input
    
    ! Initialize arrays
    do i=1,nx
        do j=1, ny
            p(i,j) = 0
            pd(i,j) = 0
            b(i,j) = 0
        end do
    end do 
    
    do i=1, nx
        y(i) = 0
        x(i) = 0
    end do
    
    
    !!!!! NOT WORKING !!
    !!!!! Inserts NaN into arrays and cant plot porperly.
    !!!!! Values correct, verified via csv plot
    ! Generate values for y
    call spaceevenly(1, ny, ymax, y)
    ! Generate values for x
    call spaceevenly(1, nx, xmax , x)
    
    b(int(ny/4),int(nx/4)) = 100
    b(int(3*ny/4),int(3*nx/4)) = -100
    
    derivYpow2 = dy**2
	derivXpow2 = dx**2

    time = nt
    do time = 1, time
            pd = p
            do j=2, nx-1
                do i=2, ny-1

                    F1term = p(j, i + 1) + p(j, i - 1)
                    S2term = p(j + 1, i) + p(j - 1, i)
                    source_term = b(j,i) * derivXpow2 * derivYpow2
                
                    numerator = (derivYpow2 * F1term) + (derivXpow2 * S2term) - source_term
                    denominator = 2 * (derivXpow2 + derivYpow2)
                    
                    p(j,i) = numerator / denominator
                    
                    p(1,j) = 0
                    p(nx,j) = 0
                    p(i,1) = 0
                    p(i,ny) = 0
        
                    p(1,1) = 0
                    p(1,nx) = 0
                    p(ny,1) = 0
                    p(ny,nx) = 0
                    
                end do
            end do
    end do
    
        ! write csv
    call write_csv(1, "l10.csv", p, nx, ny) 
    
    
    ! Plotting 
    call plt%initialize(grid=.true.,&
                        title='Plot of p',legend=.false.,mplot3d=.true.,use_colormap=.true.)
    
    call plt%plot_surface(y,x,p,label='value of p', istat=istat, &
                            cmap="cm.cividis")
    
    call plt%savefig('l10.png', &
                      pyfile='l10.py', &
                      dpi='200',istat=istat)
    
    print*, "Plotting done "
    
    print*, "ready"
    read*, input
    
end program l10

