!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for velocity and temperature             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CreateInitialConditions
      use param
      use local_arrays, only: vy,vx,temp,vz,xi
      use decomp_2d, only: xstart,xend
      use mpih
      implicit none
      integer :: j,k,i
      real :: xxx,yyy,eps,rdm

      !eps=0.01d0
      eps = 1.0d0
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=1,nxm
           vz(k,j,i)=0.0d0
           yyy=xm(k) 
           xxx=yc(j)            
           vy(k,j,i)=(2.0d0*yyy-6.0d0*yyy**2+4.0d0*yyy**3)*sin(3*xxx)*eps
           vz(k,j,i)=yyy*(1.0d0-yyy)**3*cos(4.3*zm(i))*eps

           yyy=xc(k)          
           xxx=ym(j)
           vx(k,j,i)=-yyy**2*(1.0d0-yyy)**2*cos(3.1*xxx)*eps

         enddo
        enddo
      enddo

      !assign linear temperature profile in the nodes k=2 to k=nxm
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=2,nxm
      temp(k,j,i)=tempbp(j,i)-(tempbp(j,i)-temptp(j,i))*xc(k)/xc(nx)
      enddo
      enddo
      enddo

      !temp perturbation
      !eps = 1.0d-5
      !do i=xstart(3),xend(3)
      !do j=xstart(2),xend(2)
      !do k=1,nxm
      !   call random_number(rdm)
      !   temp(k,j,i)=temp(k,j,i)+eps*rdm
      !end do
      !end do
      !end do

      !assign the boundary conditions at k=1 and k=nx
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      temp(1 ,j,i) = tempbp(j,i)
      temp(nx,j,i) = temptp(j,i)
      end do
      end do

      !assign linear temperature profile in the nodes k=2 to k=nxm
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=2,nxm
      xi(k,j,i)=xibp(j,i)-(xibp(j,i)-xitp(j,i))*xc(k)/xc(nx)
      enddo
      enddo
      enddo

      !assign the boundary conditions at k=1 and k=nx
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      xi(1 ,j,i) = xibp(j,i)
      xi(nx,j,i) = xitp(j,i)
      end do
      end do

      return                                                            
      end                                                               
