! ----------------
! ---
! Semi Implicit Spectral -- Phase Field Code
! For Solving Allen-Cahn Equation
! Parallelized using MPI
! ---
! ----------------

module nodeinfo
    implicit none
    integer :: rank
    integer :: nprocs
end module nodeinfo


program fft_ca_v1

    use init_grain_micro_
    use prepare_fft_
    use free_energ_fft_ca_v1_
    use write_vtk_grid_values_

    use nodeinfo
    use mpi

    implicit none

    include '/usr/local/include/fftw3.f'    

    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: status1

    integer :: NxNy, iflag, isolve, ngrain
    
    ! Simulation Cell Parameters, Material Parameters and Time Integration Parameters
    integer, parameter   :: Nx = 64, Ny = 64, nstep =5000, nprint = 500
    real(8), parameter :: dx = 0.5, dy = 0.5, dtime = 0.005, coefA = 1.0, mobil = 5.0, grcoef = 0.1
    
    integer(8) :: forward, backward, forw2
    real(8) :: ttime, time0, time1, grain_sum, ncount, dummy
    real(8), allocatable :: etas(:,:,:), etas_(:,:), glist(:), eta(:,:), dfdeta(:,:), eta_temp(:,:), dfdeta_temp(:,:)
    real(8), allocatable :: kx(:), ky(:), k2(:,:), k4(:,:), eta2(:,:)!, k2_shift(:,:)
    integer :: istep, igrain, i, j, k
    integer :: istart, iend, istart2, iend2
    double complex, allocatable :: etak(:,:), dfdetak(:,:), etak_temp(:,:)
    double complex :: numer, denom


    ! Get initial time

    if(rank==0) then
        call cpu_time(time0)
    endif

    call MPI_INIT(ierr)
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    open(unit = 2, file = 'area_frac.out')
    
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

    if(.not. allocated(eta_temp) ) allocate(eta_temp(Nx,Ny))
    if(.not. allocated(dfdeta_temp) ) allocate(dfdeta_temp(Nx,Ny))

    if( .not. allocated(kx) ) allocate(kx(Nx+1))
    if( .not. allocated(ky) ) allocate(ky(Ny+1))
    if( .not. allocated(k2) ) allocate(k2(Nx,Ny))
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
    if(.not. allocated(etak_temp)) allocate(etak_temp(Nx,Ny))
    if(.not. allocated(eta2)) allocate(eta2(Nx,Ny))

    etak = (0.0, 0.0)
    dfdetak = (0.0, 0.0)
    eta2 = 0.0

    call dfftw_plan_dft_2d(forward, Nx, Ny, eta, etak, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_2d(forw2, Nx, Ny, dfdeta, dfdetak, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_2d(backward, Nx, Ny, etak, etak_temp, FFTW_BACKWARD, FFTW_ESTIMATE)

    call para_range(1,Ny, rank, istart, iend)

    do istep = 1,nstep
        
        ttime = ttime + dtime

        do igrain = 1,ngrain
            
            if (glist(igrain)==1) then

                call MPI_BCAST(etas, Nx*Ny*ngrain, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                do i=1,Nx
                    do j=istart, iend

                        eta(i,j) = etas(i,j,igrain)

                        ! 
                        ! -- Derivative of free energy
                        !

                        call free_energ_fft_ca_v1(i, j, ngrain, etas, eta, igrain, dummy, Nx, Ny)

                        dfdeta(i,j) = dummy

                    enddo
                enddo

                if(rank /= 0) then
                    call MPI_SEND(istart, 1, MPI_INT, 0, rank+nprocs, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(iend, 1, MPI_INT, 0, rank+nprocs*2, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(eta, Nx*Ny, MPI_REAL8, 0, rank, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(dfdeta, Nx*Ny, MPI_REAL8, 0, rank+nprocs*3, MPI_COMM_WORLD, ierr)
                endif

                if(rank==0) then
                    do i=1,nprocs-1
                        call MPI_RECV(istart2, 1, MPI_INT, i, i+nprocs, MPI_COMM_WORLD, status1, ierr)
                        call MPI_RECV(iend2, 1, MPI_INT, i, i+nprocs*2, MPI_COMM_WORLD, status1, ierr)
                        call MPI_RECV(eta_temp, Nx*Ny, MPI_REAL8, i, i, MPI_COMM_WORLD, status1, ierr)
                        call MPI_RECV(dfdeta_temp, Nx*Ny, MPI_REAL8, i, i+nprocs*3, MPI_COMM_WORLD, status1, ierr)

                        do j = istart2, iend2
                            eta(:,j) = eta_temp(:,j)
                            dfdeta(:,j) = dfdeta_temp(:,j)
                        enddo
                    enddo
                endif


                if(rank==0) then
                    
                    etak = eta
                    dfdetak = dfdeta
 
                    call dfftw_execute_dft(forward, etak, etak)

                    call dfftw_execute_dft(forw2, dfdetak, dfdetak)
                
                endif

                call MPI_BCAST(etak, Nx*Ny, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
                call MPI_BCAST(dfdetak, Nx*Ny, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                !--
                !---- Time integration
                !-- 

                do i=1,Nx
                    do j=istart, iend

                        numer = dtime*mobil*dfdetak(i,j)

                        denom = 1.0 + dtime*coefA*mobil*grcoef*k2(i,j)

                        etak(i,j) = (etak(i,j) - numer)/denom

                    enddo
                enddo

                if(rank /= 0) then
                    call MPI_SEND(istart, 1, MPI_INT, 0, rank+nprocs, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(iend, 1, MPI_INT, 0, rank+nprocs*2, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(etak, Nx*Ny, MPI_DOUBLE_COMPLEX, 0, rank, MPI_COMM_WORLD, ierr)
                endif

                if(rank==0) then
                    do i=1,nprocs-1
                        call MPI_RECV(istart2, 1, MPI_INT, i, i+nprocs, MPI_COMM_WORLD, status1, ierr)
                        call MPI_RECV(iend2, 1, MPI_INT, i, i+nprocs*2, MPI_COMM_WORLD, status1, ierr)
                        call MPI_RECV(etak_temp, Nx*Ny, MPI_DOUBLE_COMPLEX, i, i, MPI_COMM_WORLD, status1, ierr)

                        do j = istart2, iend2
                            etak(:,j) = etak_temp(:,j)
                        enddo
                    enddo
                endif

                if(rank==0) then
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

                endif

                if(rank==0) then

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

            endif
        
        enddo

        if(rank==0) then


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
        endif
        
    enddo    

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

    if(rank==0) then
        call cpu_time(time1)
        
        write(*, *) time1-time0
    endif

    call MPI_FINALIZE(ierr)
        
end program fft_ca_v1


subroutine para_range(n1, n2, irank, istart, iend)
    use nodeinfo
    
    implicit none
    integer :: n1, n2, irank, istart, iend, iwork1, iwork2
    iwork1 = (n2-n1+1)/nprocs
    iwork2 = MOD(n2-n1+1, nprocs)

    istart = n1 + irank*iwork1 + MIN(irank, iwork2)
    iend = istart + iwork1 - 1

    if(iwork2 > irank) iend = iend+1
    return

 end subroutine para_range
