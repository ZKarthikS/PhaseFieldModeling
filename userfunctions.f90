module userfunctions
    implicit none
contains

real function init_grain_micro(Nx,Ny,dx,dy,iflag,isolve)
    real, intent(in) :: Nx, Ny, dx, dy, iflag, isolve
    real, intent(inout) :: etas, ngrain, glist   ! etas = real(Nx,Ny,2), ngrain = int?,   glist =  real(ngrain)
    


! %-------------
! % generate two grains
! %-------------

! if(iflag == 1)

!     ngrain =2;

!     %-etas(, 1) first grain
!     %-etas(, 2) second grain

!     x0 = Nx/2;
!     y0 = Ny/2;

!     radius = 14.0; % radius of second grain

!     for i=1:Nx
!         for j=1:Ny

!             ii=(i-1)*Nx+j;

!             if(isolve == 2)
!                 etas(ii,1)=1.0;
!                 etas(ii,2)=0.0;
!             else
!                 etas(i,j,1) =1.0;
!                 etas(i,j,2) =0.0;
!             end

!             xlength =sqrt((i-x0)^2+(j-y0)^2);
!             if(xlength <= radius)

!                 if(isolve == 2)
!                     etas(ii,1)=0.0;
!                     etas(ii,2)=1.0;
!                 else
!                     etas(i,j,1)=0.0;
!                     etas(i,j,2)=1.0;
!                 end

!             end %if
!         end %j
!     end %i

! end %iflag





end function init_grain_micro

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

real function prepare_fft(Nx,Ny,dx,dy)
    real, intent(in) :: Nx, Ny, dx, dy
    real, intent(inout) :: kx,ky,k2,k4

    ! kx = real(Nx + 2)
    ! ky = real(Ny + 2)
    ! k2, k4 = real(Nx, Ny) 


    ! function [kx,ky,k2,k4] = prepare_fft(Nx,Ny,dx,dy)

    !     format long;
        
    !     Nx21 = Nx/2 + 1;
    !     Ny21 = Ny/2 + 1;
        
    !     Nx2=Nx+2;
    !     Ny2=Ny+2;
        
    !     %--
        
    !     delkx=(2.0*pi)/(Nx*dx);
    !     delky=(2.0*pi)/(Ny*dy);
        
    !     %--
        
    !     for i=1:Nx21
    !         fk1=(i-1)*delkx;
    !         kx(i)=fk1;
    !         kx(Nx2-i)=-fk1;
    !     end
        
    !     for j=1:Ny21
    !         fk2=(j-1)*delky;
    !         ky(j)=fk2;
    !         ky(Ny2-j)=-fk2;
    !     end
        
    !     %---
        
    !     for i=1:Nx
    !         for j=1:Ny
        
    !             k2(i,j)=kx(i)^2+ky(j)^2;
        
    !         end
    !     end
        
    !     %%--
        
    !     k4 = k2.^2;
        
    !     end %endfunction


end function prepare_fft


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


real function free_energ_fft_ca_v1(i,j,ngrain,etas,eta,igrain)


    real, intent(in) :: i, j, ngrain, etas, eta, igrain
    real, intent(inout) :: dfdeta


    ! i,j int
    ! ngrain int
    ! etas real(Nx,Ny,ngrain)
    ! eta, dfdeta = real(Nx,Ny)





! function [dfdeta] = free_energ_fft_ca_v1(i,j,ngrain,etas,eta,igrain)

!     format long;
    
!     A=1.0;
!     B=1.0;
    
!     sum=0.0;
    
!     for jgrain=1:ngrain
    
!         if(jgrain ~= igrain)
!             sum = sum +etas(i,j,jgrain).^2;
!         end
    
!     end
    
!     dfdeta = A*(2.0*B* eta(i,j)*sum + eta(i,j)^3 - eta(i,j));
    
!     end %endfunction



end module userfunctions