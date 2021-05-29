program l9
    use pyplot_module, only : pyplot
    use functions
    
    implicit none
    
    type(pyplot) :: plt

    
    integer, parameter :: nx = 31, ny = 31, c = 1
    real(kind=8), parameter :: dx = real(2) / (nx -1)
    real(kind=8), parameter :: dy = real(2) / (ny -1)
    real(kind=8) :: pi = 4 * atan(1.0)
    character :: input
    real(kind=8), dimension(nx) :: x
    real(kind=8), dimension(ny) :: y
    integer :: i, j, istat
    real(kind=8) :: plotY = 1, plotX = 2, l1norm
    real(kind=8), dimension(nx, ny) :: p
    
    
10  format (i1)
15  format (f10.5,",")
20  format ("")
    
    ! Initialize Arrays
    do i=1,nx
        do j=1, ny
            p(i,j) = 0
        end do
    end do 
    
    do i=1, nx
        y(i) = 0
    end do
    
    ! Generate values for y
    call spaceevenly(1, ny, plotY, y)
    !Generate values for x
    call spaceevenly(1, nx, plotX , x)
    
    ! Asign values to p   
    call boundary_conditions(p, y, nx,ny)
    
    ! Laplace2d
    l1norm = 0.0004
    call laplace2d(p, y, dx, dy, l1norm, nx, ny)
    
    ! write csv
    call write_csv(1, "l9.csv", p, nx, ny) 
    
    
    ! Plotting 
    call plt%initialize(grid=.true.,&
                        title='Plot of p',legend=.false.,mplot3d=.true.,use_colormap=.true.)
    
    call plt%plot_surface(y,x,p,label='value of p', istat=istat, &
                            cmap="cm.cividis")
    
    call plt%savefig('lp.png', &
                      pyfile='l9.py', &
                      dpi='200',istat=istat)
    
    print*, "Plotting done "
    
    print*, "ready"
    read*, input
end program l9

