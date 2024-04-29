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
    use fft_shift_

    implicit none

    include '/usr/local/include/fftw3.f'
    ! include './fftw3.f'
    

    integer :: NxNy, iflag, isolve, ngrain
    integer, parameter   :: Nx = 64, Ny = 64, nstep =5000, nprint = 500
    integer(8) :: forward, backward, forw2
    real(8), parameter :: dx = 0.5, dy = 0.5, dtime = 0.005, coefA = 1.0, mobil = 5.0, grcoef = 0.1
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
    open(unit = 3, file = 'eta.txt')

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(1)

    call dfftw_planner_nthreads(iret)


    ! -- Simulation Cell parameters

    ! integer, parameter :: Nx = 64
    ! integer, parameter :: Ny = 64
    NxNy = Nx*Ny

    ! real, parameter :: dx = 0.5
    ! real, parameter :: dy = 0.5

    ! -- Time integration parameters

    ! integer, parameter :: nstep = 5000
    ! integer, parameter :: nprint = 50
    
    ! real, parameter :: dtime = 0.005
    ! real, parameter :: coefA = 1.0
    ttime = 0.0


    ! -- Material Parameters

    ! real, parameter :: mobil = 5.0
    ! real, parameter :: grcoef = 0.1



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

    ! write(*, *) 'NGrain, glist : ', ngrain, glist

    ! write(*,*) 'etas : '

    ! do igrain=1,ngrain
    !     do i=1,Nx
    !         write(*,*) etas(i,:,igrain)
    !     enddo
    ! enddo

    ! write (*, *) 'etas(blah) : ', etas(62,62,1), etas(62,62,2)
    ! write (*, *) 'glist : ', glist
    ! write (*, *) 'ngrain : ', ngrain


    ! --
    ! --- Prepare the fft
    ! --

    call prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)

    ! write(*,*) 'k2 : '
    ! do i=1,Nx 
    !     write(*,*) k2(i,:)
    ! enddo

    ! call fft_shift(k2, k2_shift, Nx)

    ! write(*,*) k4(2,1), k4(2,2), k4(64,64)


    ! --
    ! --- Evolve
    ! --

    ! eta(2,2) = etas(2,2,1)

    ! call free_energ_fft_ca_v1(2,2,ngrain, etas, eta, 2, dfdeta, Nx, Ny)

    ! write(*,*) dfdeta(2,2)

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

                ! write(*,*) 'dfdeta : '
                ! do i=1,Nx
                !     write(*,*) dfdeta(i,:)
                ! enddo

                ! write(*,*) 'eta : '
                ! do i=1,Nx
                !     write(*,*) eta(i,:)
                ! enddo

                etak = eta
                dfdetak = dfdeta
 
                call dfftw_execute_dft(forward, etak, etak)

                call dfftw_execute_dft(forw2, dfdetak, dfdetak)

                ! write(*,*) 'etak : '
                ! do i=1,Nx
                !     write(*,*) etak(i,:)
                ! enddo

                ! write(*,*) 'dfdetak : '
                ! do i=1,Nx
                !     write(*,*) dfdetak(i,:)
                ! enddo

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

                ! write(*,*) 'eta : '
                ! do i=1,Nx
                !     write(*,*) eta(i,:)
                ! enddo

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
                        
                        ! write(3, *) eta(i,j)

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

                        ! write(3, *) etas(i,j,igrain)

                    enddo
                enddo
                ! write(2,'(F14.6)') REAL(ncount), numer, denom
                ! write(3, '(F14.6)') eta2
                ! write(3,*) '\n\n\n\n'
                ncount = ncount/(1.0*NxNy)
                write(2,'(F14.6)') REAL(ncount)
            enddo

            call write_vtk_grid_values(Nx, Ny, dx, dy, istep, eta2)

        endif

        ! write(3, *) eta(16,16)!, numer, denom
        ! write(3, *) numer, denom
        ! write(3, *) dfdeta
        ! do i=1,Nx
        !     write(3,*) eta(i,:)
        ! enddo
        
    enddo

    ! write(3, *) etas(Nx,Ny,1), etas(Nx,Ny,2)
    

    close(2)
    close(3)

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

    write(*, *) time1-time0
        


end program fft_ca_v1