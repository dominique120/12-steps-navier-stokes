program l7
    use pyplot_module, only : pyplot
    use functions
    
    implicit none
    
    type(pyplot) :: plt

    ! Variables
    integer, parameter :: nx = 31, ny = 31, nt = 50, c = 1
    real(kind=8), parameter :: dx = real(2) / (nx -1 )
    real(kind=8), parameter :: dy = real(2) / (ny -1 )
    real(kind=8), parameter :: sigma = 0.25, nu = 0.05
    real(kind=8), parameter :: dt = sigma * dx * dy / nu
    real(kind=8), dimension(nx,ny) :: u, un
    real(kind=8), dimension(nx) :: x
    real(kind=8), dimension(ny) :: y
    integer :: i, j, n, time, istat
    character :: input
    real(kind=8) :: plotX = 2, plotY = 2
    
    real(kind=8) :: Ydisip, Xdisip, Yderiv, Xderiv
    
10  format (i1)  
15  format (f10.5,",")   
20  format ("")  
    
    !Generate values for y
    call spaceevenly(1, ny, plotY, y)
    !Generate values for x
    call spaceevenly(1, nx, plotX , x)
    
    ! Generate ones for entire array "u" un"
    do i=1,nx
        do j=1, ny
            u(i,j) = 1
            un(i,j) = 1
        end do
    end do 
    
    ! Set boundary conditions for "u" 
    do i=1,nx
        do j=1, ny
            if (i >= 0.5 / dx  .and. i <= 1 / dx + 1 .and. &
                j >= 0.5 / dy .and. j <= 1 / dy + 1) then
                u(i,j) = 2
            end if
        end do
    end do
    

    ! Main loop
    time = 50
    call diffuse(time, u, nu, dt, dx, dy, nx, ny)
    
    ! Write to csv
    open(unit=1, file="l7.csv")
    do i=1, nx
        do j=1, ny
            write(unit=1, fmt=15, advance="no"), u(i,j)
        end do
        write(unit=1,fmt=20)
    end do
    
    print*, "File written"
    
    
    call plt%initialize(grid=.true.,&
                        title='Plot of u',legend=.false.,mplot3d=.true.,use_colormap=.true.)
    
    call plt%plot_surface(x,y,u,label='value of u', istat=istat, &
                            cmap="cm.cividis")
    
    call plt%savefig('l7.png', &
                      pyfile='l7.py', &
                      dpi='200',istat=istat)
    
    print*, "Plotting done "
    close(unit=1)
    !read*, input

end program l7

