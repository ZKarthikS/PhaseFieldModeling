module prepare_fft_
    implicit none

contains

subroutine prepare_fft(kx, ky, k2, k4, Nx, Ny, dx, dy)

    ! Nx, Ny, dx, dy are inputs

    real(8), allocatable :: kx(:), ky(:), k2(:,:), k4(:,:), kappax(:), kappay(:)
    real(8), parameter :: pi = 4*atan(1.d0)
    real(8) :: dx, dy, delkx, delky, fk1, fk2
    integer :: Nx, Ny, Nx21, Ny21, Nx2, Ny2
    integer :: i,j

    if( .not. allocated(kx) ) allocate(kx(Nx+1))
    if( .not. allocated(ky) ) allocate(ky(Ny+1))
    if( .not. allocated(kappax) ) allocate(kappax(Nx))
    if( .not. allocated(kappay) ) allocate(kappay(Ny))
    if( .not. allocated(k2) ) allocate(k2(Nx,Ny))
    if( .not. allocated(k4) ) allocate(k4(Nx,Ny))

    ! do i=1,n
    ! enddo

    ! kappax(1) = -Nx/2.0
    ! kappay(1) = -Ny/2.0
    
    ! do i=2,(Nx/2+1)
    !     kappax(i) = (i-1)-Nx/2.0
    !     kappax(Nx+2-i) = Nx/2.0-(i-1)
    ! enddo

    ! do i=2,(Ny/2+1)
    !     kappay(i) = (i-1)-Ny/2.0
    !     kappay(Ny+2-i) = Ny/2.0-(i-1)
    ! enddo

    ! kx(1:Nx/2) = kappax(CEILING(Nx/2.0)+1:Nx)*2.0*pi/(Nx*dx)
    ! kx(Nx/2+1:Nx) = kappax(1:CEILING(Nx/2.0))*2.0*pi/(Nx*dx)

    ! ky(1:Ny/2) = kappay(CEILING(Ny/2.0)+1:Ny)*2.0*pi/(Ny*dy)
    ! ky(Ny/2+1:Ny) = kappay(1:CEILING(Ny/2.0))*2.0*pi/(Ny*dy)

    ! ky(1:Ny/2) = kappay(CEILING(Ny/2.0)+1:Ny)*2.0*pi/(Ny*dy)
    ! ky(Ny/2+1:Ny) = kappay(1:CEILING(Ny/2.0)+1)*2.0*pi/(Ny*dy)

    Nx21 = Nx/2 + 1
    Ny21 = Ny/2 + 1

    Nx2 = Nx+2
    Ny2 = Ny+2

    delkx = (2.0*pi)/(Nx*dx)
    delky = (2.0*pi)/(Ny*dy)

    do i=1,Nx21
        fk1 = (i-1)*delkx
        kx(i) = fk1;
        kx(Nx2-i) = -fk1;
    enddo

    do j=1,Ny21
        fk2 = (j-1)*delky
        ky(j) = fk2
        ky(Ny2-j) = -fk2
    enddo 

    do i=1,Nx
        do j=1,Ny

            k2(i,j) = kx(i)**2 + Ky(j)**2

        enddo
    enddo
    
    ! print *, kx
    ! print *, ky

    k4 = k2**2

    ! if(allocated(kappax)) deallocate(kappax)
    ! if(allocated(kappay)) deallocate(kappay)

    return
end subroutine prepare_fft


end module prepare_fft_
