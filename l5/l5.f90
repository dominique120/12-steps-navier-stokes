program l5
    use pyplot_module, only : pyplot
    
    implicit none
    
    type(pyplot) :: plt

    ! Variables
    integer, parameter :: nx = 101, ny = 101, nt = 80, c = 1
    real(kind=8), parameter :: dx = real(2) / (nx -1 )
    real(kind=8), parameter :: dy = real(2) / (ny -1 )
    real(kind=8), parameter :: sigma = 0.2
    real(kind=8), parameter :: dt = sigma * dx
    real(kind=8), dimension(nx,ny) :: u, un, v, vn
    real(kind=8), dimension(nx) :: x
    real(kind=8), dimension(ny) :: y
    integer :: i, j, n, time, istat
    character :: input
    real(kind=8) :: plotX = 2, plotY = 2
    
10  format (i1)  
15  format (f10.5,",")   
20  format ("")    
    
    !Generate values for y
    call spaceevenly(1, ny, plotY, y)
    !Generate values for x
    call spaceevenly(1, nx, plotX , x)
    
    ! Generate ones for entire array "u" "v" "vn" un"
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
  
    ! Print array
    do i=1, nx
        do j=1, ny
            write(*, 10, advance="no"), int(u(i,j))
        end do
        print*,""
    end do
      
    open(unit=1, file="l5.csv")
    
    ! Main loop
    do time = 1, nt + 1
        un = u
        vn = v
        
        do j=1, nx
            do i=1, ny
                !u(j,i) = (un(j,i) - (c * dt / dx * (un(j,i) - un(j,merge(1,i-1,i.eq.1)))) - &
                !    (c * dt / dy * (un(j,i) - un(merge(1,j-1,j.eq.1),i))));

                 u(j,i) = (un(j,i) - (un(j,i) * c * dt / dx * (un(j,i) - un(j,merge(1,i-1,i.eq.1)))) - &
                    vn(j,i) * c * dt / dy * (un(j,i) - un(merge(1,j-1,j.eq.1),i)))

                v(j,i) = (vn(j,i) - (un(j,i) * c * dt / dx * (vn(j,i) - vn(j,merge(1,i-1,i.eq.1)))) - &
                    vn(j,i) * c * dt / dy * (vn(j,i) - vn(merge(1,j-1,j.eq.1),i)))

                
                u(1,j) = 1;
                u(nx-1,j) = 1;
                u(i,1) = 1;
                u(i,ny-1) = 1;
                
                v(1,j) = 1;
                v(nx-1,j) = 1;
                v(i,1) = 1;
                v(i,ny-1) = 1;
            end do
        end do
    end do
    
    ! Write to csv
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
    
    call plt%savefig('l5.png', &
                      pyfile='l5.py', &
                      dpi='200',istat=istat)
    close(unit=1)
    read*, input
    
end program l5

