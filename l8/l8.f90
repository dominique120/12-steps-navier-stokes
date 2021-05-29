program l8
    use pyplot_module, only : pyplot
    use functions
    
    implicit none
    
    type(pyplot) :: plt

    ! Variables
    integer, parameter :: nx = 41, ny = 41, nt = 120, c = 1
    real(kind=8), parameter :: dx = real(2) / (nx -1 )
    real(kind=8), parameter :: dy = real(2) / (ny -1 )
    real(kind=8), parameter :: sigma = 0.0009, nu = 0.01
    real(kind=8), parameter :: dt = sigma * dx * dy / nu
    real(kind=8), dimension(nx,ny) :: u, un, v, vn
    real(kind=8), dimension(nx) :: x
    real(kind=8), dimension(ny) :: y
    integer :: i, j, n, time, istat
    character :: input
    real(kind=8) :: plotX = 2, plotY = 2
    
    real(kind=8) :: dtOdx, dtOdy, nudtOdx2, nudtOdy2, F1termU, S2termU, T3termU, F4termU, F1termV, S2termV, T3termV, F4termV


    ! Body of l8
10  format (i1)  
15  format (f10.5,",")   
20  format ("")  
    
    !Generate values for y
    call spaceevenly(1, ny, plotY, y)
    !Generate values for x
    call spaceevenly(1, nx, plotX , x)
    
    ! Generate ones for entire array "u" "un" "v" "vn"
    do i=1,nx
        do j=1, ny
            u(i,j) = 1
            un(i,j) = 1
            v(i,j) = 1
            vn(i,j) = 1
        end do
    end do 
    
    ! Set boundary conditions for "u" "v"
    do i=1,nx
        do j=1, ny
            if (i >= 0.5 / dx  .and. i <= 1 / dx + 1 .and. &
                j >= 0.5 / dy .and. j <= 1 / dy + 1) then
                u(i,j) = 2
                v(i,j) = 2
            end if
        end do
    end do
    
    
    time = nt
    do time = 1, time
            un = u
            do j=2, nx-1
                do i=2, ny-1
                    nudtOdx2 = (nu * dt) / (dx**2)
                    nudtOdy2 = (nu * dt) / (dy**2)
                    dtOdx = dt / dx
                    dtOdy = dt / dy
                    
                    F1termU = u(j,i) * (u(j,i) - u(j-1,i))
                    S2termU = v(j,i) * (u(j,i) - u(j,i-1))
                    T3termU = u(j+1,i) - (2 * u(j,i)) + u(j-1,i)
                    F4termU = u(j,i+1) - (2 * u(j,i)) + u(j,i-1)
                    
                    F1termV = v(j,i) * (v(j,i) - v(j-1,i))
                    S2termV = v(j,i) * (v(j,i) - v(j,i-1))
                    T3termV = v(j+1,i) - (2 * v(j,i)) + v(j-1,i)
                    F4termV = v(j,i+1) - (2 * v(j,i)) + v(j,i-1)

                    u(j,i) = un(j,i) - (dtOdx * F1termU) - &
                        (dtOdy * S2termU) + &
                        (nudtOdx2 * T3termU) + &
                        (nudtOdy2 * F4termU)
                    
                    v(j,i) = vn(j,i) - (dtOdx * F1termV) - &
                        (dtOdy * S2termV) + &
                        (nudtOdx2 * T3termV) + &
                        (nudtOdy2 * F4termV)

                    u(1,j) = 1
                    u(nx,j) = 1
                    u(i,1) = 1
                    u(i,ny) = 1
        
                    u(1,1) = 1
                    u(1,nx) = 1
                    u(ny,1) = 1
                    u(ny,nx) = 1
                    
                    v(1,j) = 1
                    v(nx,j) = 1
                    v(i,1) = 1
                    v(i,ny) = 1
        
                    v(1,1) = 1
                    v(1,nx) = 1
                    v(ny,1) = 1
                    v(ny,nx) = 1
                end do
            end do
        end do
    

        ! Write to csv
    open(unit=1, file="l8.csv")
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
    
    call plt%savefig('l8.png', &
                      pyfile='l88.py', &
                      dpi='200',istat=istat)
    
    print*, "Plotting done "
    close(unit=1)
    read*, input

end program l8

