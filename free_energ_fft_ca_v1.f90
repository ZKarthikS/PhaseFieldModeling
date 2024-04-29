module free_energ_fft_ca_v1_
    implicit none

contains

subroutine free_energ_fft_ca_v1(i, j, ngrain, etas, eta, igrain, dummy, Nx, Ny)

    ! i, j, ngrain, etas, eta, igrain - Inputs
    ! Nx, Ny - Also inputs

    integer :: i, j, igrain, jgrain, ngrain, Nx, Ny
    real(8) :: A, B, sum, dummy
    real(8), allocatable :: etas(:,:,:), eta(:,:)!, dfdeta(:,:)

    if( .not. allocated(etas) ) allocate(etas(Nx, Ny, ngrain))
    if( .not. allocated(eta) ) allocate(eta(Nx, Ny))
    ! if( .not. allocated(dfdeta) ) allocate(dfdeta(Nx,Ny))

    A = 1.0
    B = 1.0
    sum = 0.0

    do jgrain = 1, ngrain

        if(jgrain /= igrain) then
            sum = sum + etas(i,j,jgrain)**2
        endif

    enddo

    ! dfdeta(i,j) = A*(2.0*B*eta(i,j)*sum + eta(i,j)**3 - eta(i,j))
    dummy  = A*(2.0*B*eta(i,j)*sum + eta(i,j)**3 - eta(i,j))

    return 
end subroutine free_energ_fft_ca_v1


end module free_energ_fft_ca_v1_

