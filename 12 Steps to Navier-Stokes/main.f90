program ns
    use, intrinsic :: iso_fortran_env, only : wp => real64
    use pyplot_module, only : pyplot
    
    implicit none
    
    type(pyplot) :: plt
    
    real(kind=8), parameter :: nx = 41 ! Number of grid points
    real(kind=8) :: dx, nt, dt, c, i, n, j, nu, sigma
    real(kind=8), dimension(nx) :: u, un, grid ! 
    integer :: istat
    character :: input 
    
    !nx is declared and assigned above        
    dx = 2 / (nx-1) ! Distance between points
    nt = 20         ! Number of timesteps
    nu = 0.3        ! Viscosity
    sigma = 0.2
    dt = sigma * dx**2 / nu
    
    ! Initialize u to match our initial condition of 
    ! u = 2 in the following range: 0.5 > x > 1 
    ! and 1 everywhere else
    do i=1, nx
        u(i) = 1
        if (i .ge. (.5 / dx) .and. i .le. (1 / dx + 1)) then 
            u(i) = 2
        end if 
    end do
    
    ! Use this array to represent the gridpoints for plotting
    do i=1, nx
        grid(i) = i
    end do
    
    ! Copy values from u to un
    do i=1, nx
        un(i) = u(i)
    end do
    
    ! Solve convection
    do n=1, nt
        do j =1, nx
            un(j)=u(j)
        end do
        
        do i=2, nx-1
            u(i) = un(i) + nu * dt / dx**2 * (un(i+1) - 2 * un(i) + un(i-1))
        end do
    end do
    
        
    
    ! Print the values contained in u
    do i=1, nx
        print*, u(i)
    end do
    
    print*, "plotting"
    call plt%initialize(grid=.true.,xlabel='timestep',&
                        title='Plot of u over time',legend=.true.)
    call plt%add_plot(grid,u,label='value of u',linestyle='b-o',markersize=3,linewidth=1,istat=istat)
    call plt%savefig('plot.png', pyfile='plot.py', dpi='320', istat=istat)
    
    read*, input  
end program ns
 