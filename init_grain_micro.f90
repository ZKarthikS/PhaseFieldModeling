module init_grain_micro_
    implicit none

contains

subroutine init_grain_micro(Nx, Ny, dx, dy, iflag, isolve, etas, ngrain, glist)

    ! Nx, Ny, dx, dy, iflag, isolve - Inputs

    integer :: Nx, Ny, iflag, isolve, ngrain
    real(8) :: dx, dy
    real(8), allocatable :: etas(:, :, :), glist(:), etas_(:,:)

    integer :: x0, y0, i, j, ii, ndime
    integer :: nvpoin, nvnode, nvelem, ipoin, jpoin
    integer :: ielem, jelem, idime, inode, igrain, jnode
    real(8) :: radius, xlength, knode, theta, xv1, yv1
    real(8) :: xv2, yv2, p1x, p1y, p2x, p2y, x1, x2, tx1, mnode
    real(8), parameter :: twopi = 8*atan(1.d0)
    real(8), parameter :: epsilon = 1.0E-4
    real(8), allocatable :: dummy(:), vcord(:,:), vlnods(:,:), nnode2(:), gx(:), gy(:)

    
    if(isolve==2) then
        if( .not. allocated(etas_) ) allocate(etas_(Nx*Ny, ngrain))
    else
        if( .not. allocated(etas) ) allocate(etas(Nx, Ny, ngrain))
    endif


    !--------------------------------------------
    ! Generate two grains
    !--------------------------------------------

    if(iflag==1) then

        ngrain = 2

        if( .not. allocated(glist) ) allocate(glist(ngrain))

        ! etas(,1) => first grain
        ! etas(,2) => second grain

        x0 = Nx/2
        y0 = Ny/2

        radius = 14.0       ! Radius of second grain

        do i=1,Nx
            do j=1,Ny

                ii = (i-1)*Nx+j

                if(isolve==2) then
                    etas_(ii,1)=1.0
                    etas_(ii,2)=0.0
                else
                    etas(i,j,1)=1.0
                    etas(i,j,2)=0.0
                endif

                xlength = SQRT(REAL((i-x0)**2 + (j-y0)**2))

                if(xlength<=radius) then

                    if(isolve==2) then
                        etas_(ii,1) = 0.0
                        etas_(ii,2) = 1.0
                    else
                        etas(i,j,1) = 0.0
                        etas(i,j,2) = 1.0
                    endif
                endif
                
            enddo
        enddo

    endif


    !-------------------------------------------------
    ! Generate Polycrystal microstructure
    !-------------------------------------------------

    if(iflag==2) then

        !--------
        ! Read data generated from voroni_1.m
        !--------

        open (unit = 1, file = './grain_25.inp')

        ndime = 2

        read(1,*) nvpoin
        read(1,*) nvnode
        read(1,*) nvelem
        read(1,*) ngrain

        if(.not. allocated(dummy)) allocate(dummy(2))
        if(.not. allocated(vcord)) allocate(vcord(nvelem, ndime))
        if(allocated(glist)) deallocate(glist)
        if( .not. allocated(glist) ) allocate(glist(ngrain))

        do ipoin=1,nvpoin
            read(1,*) jpoin

            read(1,*) dummy(1), dummy(2)
            
            do idime = 1,ndime
                vcord(jpoin, idime) = dummy(idime)
            enddo
        enddo

        if(allocated(dummy)) deallocate(dummy)
        if(.not. allocated(dummy)) allocate(dummy(nvnode+1))
        if(.not. allocated(vlnods)) allocate(vlnods(nvelem, nvnode+1))

        do ielem=1,nvelem
            read(1,*) jelem

            do i=1,nvnode+1
                read(1,*) dummy(i)
            enddo

            do inode=1,nvnode+1
                vlnods(ielem, inode) = dummy(inode)
            enddo
        enddo

        !-------

        if(.not. allocated(nnode2)) allocate(nnode2(nvelem))

        do ielem = 1,nvelem

            jnode = 0

            do inode = 1,nvnode
                knode = vlnods(ielem, inode)
                if(knode /= 0) jnode  = jnode + 1
            enddo

            nnode2(ielem) = jnode
        enddo


        !----------
        ! Form the grid
        !----------

        if(.not. allocated(gx)) allocate(gx(Nx))
        if(.not. allocated(gy)) allocate(gy(Ny))

        do i=1,Nx
            gx(i) = i*dx
        enddo

        do j=1,Ny
            gy(j) = j*dy
        enddo


        !---------------
        ! Initialize order parameters
        !---------------

        do i=1,Nx
            do j=1,Ny
                do igrain=1,ngrain

                    if(isolve==2) then
                        ii = (i-1)*Nx+j
                        etas_(ii,igrain) = 0.0
                    else
                        etas(i,j,igrain) = 0.0
                    endif

                enddo
            enddo
        enddo

        !--------
        !--------

        do i=1,Nx
            do j=1,Ny

                ii = (i-1)*Nx+j

                do ielem=1,nvelem

                    igrain = vlnods(ielem, nvnode+1)
                    theta = 0.0
                    mnode = nnode2(ielem)

                    do inode=1,INT(mnode)

                        knode = vlnods(ielem, inode)

                        xv1 = vcord(INT(knode),1)
                        yv1 = vcord(INT(knode),2)

                        jnode = vlnods(ielem, inode+1)
                        if(inode==mnode) jnode = vlnods(ielem,1)

                        xv2 = vcord(jnode,1)
                        yv2 = vcord(jnode,2)

                        p1x = xv1 - gx(i)
                        p1y = yv1 - gy(j)

                        p2x = xv2 - gx(i)
                        p2y = yv2 - gy(j)

                        x1 = SQRT(p1x**2 + p1y**2)
                        x2 = SQRT(p2x**2 + p2y**2)

                        if(x1*x2<=epsilon) then 
                            theta = twopi
                        else
                            tx1 = ((p1x*p2x + p1y*p2y)/(x1*x2))
                        endif

                        if(abs(tx1)>=1.0) tx1 = 0.9999999999

                        theta = theta + acos(tx1)

                    enddo

                    if (abs(theta-twopi) <= epsilon) then
                        
                        if (isolve==2) then
                            etas_(ii, igrain) = 1.0
                        else
                            etas(i,j,igrain) = 1.0
                        endif

                    endif

                enddo
            
            enddo
        enddo


        !------
        close(1)

    endif

    ! --- Initialize glist

    do igrain=1,ngrain
        glist(igrain) = 1.0
    enddo


    ! if (allocated(etas)) deallocate(etas)
    ! if (allocated(glist)) deallocate(glist)
    ! if (allocated(etas_)) deallocate(etas_)
    ! if (allocated(dummy)) deallocate(dummy)
    ! if (allocated(vcord)) deallocate(vcord)
    ! if (allocated(vlnods)) deallocate(vlnods)
    ! if (allocated(nnode2)) deallocate(nnode2)
    ! if (allocated(gx)) deallocate(gx)
    ! if (allocated(gy)) deallocate(gy)

end subroutine init_grain_micro


end module init_grain_micro_

