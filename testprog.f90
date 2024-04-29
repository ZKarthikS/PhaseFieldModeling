program check
    ! use, intrinsic :: iso_c_binding 
    use testfun
    use prepare_fft_
    use fft_shift_

    ! use mpi

    implicit none

    include '/usr/local/include/fftw3.f'
    ! include '/home/zkarthik/fftw-3.3.10/mpi/fftw3-mpi.f03'

    real :: x, xbegin, xend, f_(1:10)
    real(8) :: y(1:10), eta(3,3), in_(7), const
    real(8) :: eta_(6,6), time0, time1
    real(8), allocatable :: kx(:), ky(:), k2(:,:), k4(:,:), k4_shift(:,:)
    double complex :: etak(3,3), out_(7)
    double complex :: etak_(6,6), eta_new(3,3), etak_new(3,3)
    complex, parameter :: IJ = (0,1)
    integer :: i, steps, iret
    integer(8) :: forward, backward, forward1, backward1, forward2, backward2
    character(len=10) :: a,b,fmt,x1
    character(len=20) :: fname

    real, allocatable :: z(:), z_(:,:)

    call cpu_time(time0)

    call dfftw_init_threads(iret)
    call dfftw_plan_with_nthreads(4)

    ! write(*,*) 'Please enter value of x'
    ! read(*,*) x

    ! eta = transpose(reshape((/ 1, 1, 3, 6, 5, 6, 8, 8, 9 /), shape(eta)))
    eta_ = transpose(reshape((/1, 1, 3, 2, 4, 6, 7, 9, 9, 10, 11, 12, 12, 13, 14, &
            16, 19, 16, 18, 20, 21, 22, 23, 25, 24, &
            26, 27, 28, 31, 30, 29, 35, 32, 34, 33, 36/), shape(eta_)))

    ! eta_new = (reshape((/(45,0), (-13.5,7.7942), (-4.5,2.5981), (0,0), &
    !         (-4.5,-2.5981), (0,0), (0,0), (0,0), (0,0)/), shape(eta_new)))

    eta_new = transpose((reshape((/(1,0), (2,0), (3,0), (4,0), &
            (5,0), (6,0), (7,0), (8,0), (9,0)/), shape(eta_new))))

    ! if(mod(x,2.0)==0) then
    !     real, allocatable :: z(:)
    ! else
    !     real, allocatable :: z(:,:)
    ! endif

    ! if(mod(x,2.0)==0) then
    !     if(.not. allocated(z)) allocate(z(5))
    !     write(*,*) 'z : ', z
    ! else
    !     if(.not. allocated(z_)) allocate(z_(2,2))
    !     write(*,*) 'z : ', z_
    ! endif

    ! write(*,*) 'Obtained f(x) : ', f(x)

    ! a = 'time_'
    ! b = '.vtk'

    ! fmt = '(1x,i5)'

    ! write(x1, fmt) i

    ! x1 = adjustl(x1)

    ! fname = 'time_'//trim(x1)//'.vtk'

    ! write(*,*) 'fname : ' , fname

    ! write(a,'1x,i5', b) "time_", i, ".vtk"
    ! write('time_','(1x,i5)','.vtk') i

    ! i = i+1

    !fname = write('time_','1x,i5','.vtk') i
    const = 1.0
    call prepare_fft(kx, ky, k2, k4, 3, 3, const, const)
    call fft_shift(k4, k4_shift, 3)

    call dfftw_plan_dft_r2c_2d(forward, 3, 3, eta, etak, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(backward, 3, 3, etak, eta, FFTW_ESTIMATE)

    call dfftw_plan_dft_2d(forward2, 3, 3, eta_new, etak_new, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_2d(backward2, 3, 3, etak_new, eta_new, FFTW_BACKWARD, FFTW_ESTIMATE)

    call dfftw_plan_dft_r2c_2d(forward1, 6, 6, eta_, etak_, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(backward1, 6, 6, etak_, eta_, FFTW_ESTIMATE)

    ! call dfftw_plan_dft_r2c_1d(forward1, 7, in_, out_, FFTW_ESTIMATE)
    ! call dfftw_plan_dft_c2r_1d(backward1, 7, out_, in_, FFTW_ESTIMATE)
    
    ! do i = 1,10
    !     y(i) = i
    ! enddo

    ! call calc_(x, y, f_)

    ! write(*,*) 'New Obtained f = x/2 ', f_

    ! write(*,*) 'Eta : '
    ! do i=1,6
    !     write(*,*)  eta_(i, :)
    ! enddo

    write (*,*) 'ETA:'
    do i=1,3
        write(*,*)  eta_new(i, :)
    enddo

    ! k4 = transpose(reshape((/ 1, 2, 3, 2, 3, 4, 3, 4, 5 /), shape(eta)))

    call dfftw_execute_dft_r2c(forward1, eta_, etak_)
    ! write(*,*) 'Initial ETAK : '
    ! do i=1,6
    !     write(*,*)  etak_(i, :)
    ! enddo 
    ! etak_ = 2.0*IJ*etak_
    ! write(*,*) 'Final ETAK : '
    ! do i=1,6
    !     write(*,*)  etak_(i, :)
    ! enddo 
    
    ! call dfftw_execute_dft_c2r(backward1, etak_, eta_)

    ! eta_new = 2*IJ*eta_new

    ! call dfftw_execute_dft_c2r(backward, eta_new, eta)

    call dfftw_execute_dft(forward2, eta_new, etak_new)

    write (*,*) 'ETAK:'
    do i=1,3
        write(*,*)  etak_new(i, :)
    enddo

    etak_new = IJ*etak_new*k4

    write (*,*) 'ETAK Final:'
    do i=1,3
        write(*,*)  etak_new(i, :)
    enddo    

    call dfftw_execute_dft(backward2, etak_new, eta_new)

    
    ! write (*,*) 'ETA:'
    ! do i=1,6
    !     write(*,*)  eta_(i, :)/(36.0)
    ! enddo

    write (*,*) 'ETA:'
    do i=1,3
        write(*,*)  eta_new(i, :)/(9.0)
    enddo

    eta = REALPART(eta_new)/9.0

    write (*,*) 'ETA:'
    do i=1,3
        write(*,*)  eta(i, :)
    enddo

    ! write(*,*) 'k4'
    ! write(*,*) k4

    ! write(*,*) 'k2'
    ! write(*,*) k2

    ! write(*,*) 'kx'
    ! write(*,*) kx

    ! in_ = (/1, 2, 3, 4, 5, 6, 7/)

    ! ! write(*,*) in_
    ! ! write(*,*) out_
    
    ! call dfftw_execute_dft_r2c(forward1, in_, out_)

    ! write(*,*) out_
    ! out_ = out_*1.0

    ! call dfftw_execute_dft_c2r(backward1, out_, in_)

    ! write(*,*) in_/(7.0)

    ! do i=1,6
    !     write(*,*)  k4(i, :)
    ! enddo

    ! write(*,*) 'Well WEll, lets take a break'

    ! do i=1,6
    !     write(*,*)  k4_shift(i, :)
    ! enddo



    call dfftw_destroy_plan(forward, eta, etak)
    call dfftw_destroy_plan(backward, etak, eta)
    ! call dfftw_destroy_plan(forward1, in_, out_)
    ! call dfftw_destroy_plan(backward1, out_, in_)

    call dfftw_cleanup_threads();

    call cpu_time(time1)
    write(*,*) time1-time0

end program check