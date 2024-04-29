module write_vtk_grid_values_
    implicit none

contains

subroutine write_vtk_grid_values(Nx, Ny, dx, dy, istep, data1)

    integer :: Nx, Ny, istep, Nz, npoin, i, j, ii
    real(8) :: dx, dy, x, y, z
    real(8), allocatable :: data1(:,:)
    character(len=20) :: fmt, x1, fname

    if( .not. allocated(data1) ) allocate(data1(Nx, Ny))


    fmt = '(i5)'
    write(x1, fmt) istep
    x1 = adjustl(x1)

    fname = 'time_'//trim(x1)//'.vtk'

    open (unit = 1, file = fname)

    Nz = 1

    npoin = Nx*Ny*Nz

    ! Start writing ASCII VTK file


    ! header of VTK file

    write(1, '(A)') '# vtk DataFile Version 2.0'
    write(1, '(A)') 'time_10.vtk'
    write(1, '(A)') 'ASCII'
    write(1, '(A)') 'DATASET STRUCTURED_GRID'

    !----- Coords of Grid points

    fmt = '(1x,i7)'
    write(1, '(A, i7, i7, i7)') 'DIMENSIONS ', Nx, Ny, Nz
    write(1, '(A, i7, A)') 'POINTS ', npoin, ' float'

    do i=1,Nx
        do j=1,Ny

            x = (i-1)*dx
            y = (j-1)*dy
            z = 0.0

            write(1, '(F14.6, F14.6, F14.6)') x, y, z

        enddo
    enddo

    !------ Write Grid point values

    write(1, '(A, i5)') 'POINT_DATA ', npoin
    write(1, '(A)') 'SCALARS CON float 1'
    write(1, '(A)') 'LOOKUP_TABLE default'

    do i=1,Nx
        do j=1,Ny
            ii=(i-1)*nx+j

            write(1, '(F14.6)') data1(i,j)
        enddo
    enddo

    close(1)

end subroutine write_vtk_grid_values


end module write_vtk_grid_values_

