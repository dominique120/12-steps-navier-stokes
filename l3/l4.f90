program ns
    use, intrinsic :: iso_fortran_env, only : wp => real64
    use pyplot_module, only : pyplot
    
    implicit none
    
    type(pyplot) :: plt
    
    integer, parameter :: nx = 101 ! Number of grid points
    integer, parameter :: nt = 100 ! Number of timesteps
    real(kind=8) :: dx, dt, c, nu
    real(kind=8) :: pi, t, x1, x2, x3, sigma
    real(kind=8) :: burger
    real(kind=8), dimension(nx) :: u, un, x, u_analytical
    integer :: istat, finish = nx
    integer :: n, i, j
    character :: input
    character(len=3) :: index_str
    pi = 4 * atan(1.0)
    
    !nx is declared and assigned above        
    dx = 2 * pi / (nx-1) ! Distance between points     
    nu = 0.07        ! Viscosity
    sigma = 0.2
    dt = dx * nu
    t = 0
    
    ! Populate "x" with evenly spaced numbers -> numpy.linspace
    ! Initialize "x"
    ! Done
    call spaceevenly(1, finish,  2 * pi, x)
    
    ! Print x     
    do i = 1, nx
        print*, x(i)
    end do
    
    ! Initialize u   
    ! Done
    do i = 1, nx
        u(i) = burger(t, x(i), nu)
    end do
    
    ! Initialize un
    ! Done
    do i = 1, nx
        un(i) = 1
    end do
    
    !main loop    
    do n=1, nt
        un = u
        do i=2, (nx - 1)
            u(i) = un(i) - un(i) * dt / dx * (un(i) - un(i-1)) + nu * dt / dx**2 * &
                (un(i+1) - 2 * un(i) + un(i-1))
            
            u(1) = un(1) - un(1) * dt / dx * (un(1) - un(nx-2)) + nu * dt / dx**2 * &
                (un(2) - 2 * un(1) + un(nx-2))
            
            u(nx) = u(1)
        end do
        write(index_str, '(I0)') n
        !call plt%initialize(grid=.true.,xlabel='timestep',&
                    !title='Plot of u over time',legend=.true.)
        !call plt%add_plot(x,u,label='value of u',linestyle='b-o',markersize=3,linewidth=1,istat=istat)

        !call plt%savefig('plot'//index_str//'.png', pyfile='plot'//index_str//'.py', dpi='220', istat=istat)
        print*, "fig ", n, " saved ot of ", nt
   end do
    
    !initialize Analitical
    do i = 1, nx
        u_analytical(i) = burger(nt*dt, x(i), nu)
    end do
    
    print*, "un"
    do i=1, nx
        print*, un(i)
    end do
    
    print*, "u"
    do i=1, nx
        print*, u(i)
    end do
    
    print*, "u_analytical"
    do i=1, nx
        print*, u_analytical(i)
    end do
    
    print*, "x"
    do i=1, nx
        print*, x(i)
    end do
    

!    print*, "plotting"
!    call plt%initialize(grid=.true.,xlabel='timestep',&
!                        title='Plot of u over time',legend=.true.)
!    call plt%add_plot(x,u,label='value of u',linestyle='b-o',markersize=3,linewidth=1,istat=istat)
!    call plt%add_plot(x,un,label='value of un',linestyle='b-o',markersize=3,linewidth=1,istat=istat)
!    call plt%add_plot(x,u_analytical,label='value of u_analytical',linestyle='b-o',markersize=3,linewidth=1,istat=istat)
!    call plt%savefig('plot.png', pyfile='plot.py', dpi='320', istat=istat)

    read*, INPUT
    
end program ns