! ----------------
! ---
! Semi Implicit Spectral -- Phase Field Code
! For Solving Allen-Cahn Equation
! ---
! ----------------

program fft_ca_v1

    use init_grain_micro_
    use prepare_fft_
    use free_energ_fft_ca_v1_
    use write_vtk_grid_values_

    implicit none

    include '/usr/local/include/fftw3.f'    

    integer :: NxNy, iflag, isolve, ngrain
    
    ! Simulation Cell Parameters, Material Parameters and Time Integration Parameters
    integer, parameter   :: Nx = 64, Ny = 64, nstep =5000, nprint = 500
    real(8), parameter :: dx = 0.5, dy = 0.5, dtime = 0.005, coefA = 1.0, mobil = 5.0, grcoef = 0.1
    
    integer(8) :: forward, backward, forw2
    real(8) :: ttime, time0, time1, grain_sum, ncount, dummy!, numer, denom
    real(8), allocatable :: etas(:,:,:), etas_(:,:), glist(:), eta(:,:), dfdeta(:,:)
    real(8), allocatable :: kx(:), ky(:), k2(:,:), k4(:,:), eta2(:,:)!, k2_shift(:,:)
    integer :: istep, igrain, i, j
    integer :: iret
    double complex, allocatable :: etak(:,:), dfdetak(:,:), eta_interim(:,:)
    double complex :: numer, denom

    ! Get initial time

    call cpu_time(time0)

    open(unit = 2, file = 'area_frac.out')

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(1) ! Change the number inside the brackets to change the number of threads

    call dfftw_planner_nthreads(iret)

    NxNy = Nx*Ny
    
    ttime = 0.0


    iflag = 1
    isolve = 1


    ! if( .not. allocated(glist) ) allocate(glist(2))
    
    if(isolve==2) then
        if( .not. allocated(etas_) ) allocate(etas_(Nx*Ny, ngrain))
    else
        if( .not. allocated(etas) ) allocate(etas(Nx, Ny, ngrain))
    endif

    if(.not. allocated(eta) ) allocate(eta(Nx,Ny))
    if(.not. allocated(dfdeta) ) allocate(dfdeta(Nx,Ny))

    if( .not. allocated(kx) ) allocate(kx(Nx+1))
    if( .not. allocated(ky) ) allocate(ky(Ny+1))
    if( .not. allocated(k2) ) allocate(k2(Nx,Ny))
    ! if( .not. allocated(k2_shift) ) allocate(k2_shift(Nx,Ny))
    if( .not. allocated(k4) ) allocate(k4(Nx,Ny))


    etas = 0.0
    eta = 0.0
    dfdeta = 0.0
    
    ! --
    ! --- Generate initial grain_structure
    ! -- 

    call init_grain_micro(Nx, Ny, dx, dy, iflag, isolve, etas, ngrain, glist)

    ! --
    ! --- Prepare the fft
    ! --

    call prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)

    ! --
    ! --- Evolve
    ! --

    if(.not. allocated(etak)) allocate(etak(Nx,Ny))
    if(.not. allocated(dfdetak)) allocate(dfdetak(Nx,Ny))
    if(.not. allocated(eta_interim)) allocate(eta_interim(Nx,Ny))
    if(.not. allocated(eta2)) allocate(eta2(Nx,Ny))

    etak = (0.0, 0.0)
    dfdetak = (0.0, 0.0)
    eta2 = 0.0

    call dfftw_plan_dft_2d(forward, Nx, Ny, eta, etak, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_2d(forw2, Nx, Ny, dfdeta, dfdetak, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_2d(backward, Nx, Ny, etak, eta_interim, FFTW_BACKWARD, FFTW_ESTIMATE)

    do istep = 1,nstep
        
        ttime = ttime + dtime

        do igrain = 1,ngrain
            
            if (glist(igrain)==1) then

                do i=1,Nx
                    do j=1,Ny

                        eta(i,j) = etas(i,j,igrain)

                        ! 
                        ! -- Derivative of free energy
                        !

                        call free_energ_fft_ca_v1(i, j, ngrain, etas, eta, igrain, dummy, Nx, Ny)

                        dfdeta(i,j) = dummy

                    enddo
                enddo

                etak = eta
                dfdetak = dfdeta
 
                call dfftw_execute_dft(forward, etak, etak)

                call dfftw_execute_dft(forw2, dfdetak, dfdetak)

                !--
                !---- Time integration
                !-- 

                do i=1,Nx
                    do j=1,Ny

                        numer = dtime*mobil*dfdetak(i,j)

                        denom = 1.0 + dtime*coefA*mobil*grcoef*k2(i,j)

                        etak(i,j) = (etak(i,j) - numer)/denom

                    enddo
                enddo

                call dfftw_execute_dft(backward, etak, etak)
                
                eta = REALPART(etak)

                ! --
                ! --

                eta = eta/(1.0*NxNy)
                grain_sum = 0.0

                do i=1,Nx
                    do j=1,Ny

                        !-- For small deviations

                        if(eta(i,j) >= 0.9999) eta(i,j) = 0.9999

                        if(eta(i,j) < 0.0001) eta(i,j) = 0.0001

                        grain_sum = grain_sum + eta(i,j)

                        etas(i,j,igrain) = eta(i,j)

                    enddo
                enddo

                ! -- Check Volume fraction of current grain :

                grain_sum = grain_sum
                if(grain_sum <= 0.001) then
                    glist(igrain) = 0
                    write(*, *) 'grain: ', igrain, ' is eliminated'
                endif

            endif
        
        enddo

        if(mod(istep,nprint) == 0) then
            write(*, *) 'done step: ', istep

            write(2, '(F14.6)') ttime

            eta2 = 0.0

            do igrain = 1,ngrain
                ncount = 0.0
                do i=1,Nx
                    do j=1,Ny

                        eta2(i,j) = eta2(i,j) + etas(i,j,igrain)**2

                        if(etas(i,j,igrain) >= 0.5) then
                            ncount = ncount+1
                        endif

                    enddo
                enddo
                ncount = ncount/(1.0*NxNy)
                write(2,'(F14.6)') REAL(ncount)
            enddo

            call write_vtk_grid_values(Nx, Ny, dx, dy, istep, eta2)

        endif
        
    enddo
    

    close(2)

    ! if(allocated(etas_)) deallocate(etas_)
    ! if(allocated(eta)) deallocate(eta)
    ! if(allocated(etas)) deallocate(etas)
    ! if(allocated(glist)) deallocate(glist)
    ! if(allocated(kx)) deallocate(kx)
    ! if(allocated(ky)) deallocate(ky)
    ! if(allocated(k2)) deallocate(k2)
    ! if(allocated(k4)) deallocate(k4)
    ! if(allocated(dfdeta)) deallocate(dfdeta)

    call dfftw_destroy_plan(forward)
    call dfftw_destroy_plan(forw2)
    call dfftw_destroy_plan(backward)

    call dfftw_cleanup_threads();

    call cpu_time(time1)

    write(*, *) "Time Taken : ", time1-time0
        


end program fft_ca_v1
