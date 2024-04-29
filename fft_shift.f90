module fft_shift_
    implicit none

contains

subroutine fft_shift(k4, k4_shift, n)

    ! Nx, Ny, dx, dy are inputs

    real(8), allocatable :: k4(:,:), k4_shift(:,:)
    integer :: n
    integer :: i,j, l1, l2, m1, m2, N_

    if( .not. allocated(k4) ) allocate(k4(n,n))
    if( .not. allocated(k4_shift) ) allocate(k4_shift(n,n))

    ! do i=1,n
    ! enddo

    N_ = int(n/2) + 1
    do i=1, N_*n
        
        l1 = mod(i,N_)
        l2 = i/N_ + 1

        if(l1==0) then
            l1 = N_
            l2 = l2-1
        endif

        m1 = mod(i,n)
        m2 = i/n + 1

        if(m1==0) then
            m1 = n
            m2 = m2-1
        endif

        k4_shift(m1, m2) = k4(l1,l2)

    enddo



    ! if(allocated(kappax)) deallocate(kappax)
    ! if(allocated(kappay)) deallocate(kappay)

    return
end subroutine fft_shift


end module fft_shift_
