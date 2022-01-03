module UpdateVel

   use param
   use local_arrays, only : vx,rhs,rux,qcap,pr,vy,ruy,dph, &
       &                    vz,dq,ruz
   use decomp_2d, only: xstart,xend

   implicit none

   private

   public :: ImplicitAndUpdateVX, ImplicitAndUpdateVY, ImplicitAndUpdateVZ

contains

   subroutine ImplicitAndUpdateVX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the X (vertical) direction and call !
!     the implicit solver. After this routine, the        !
!     vertical velocity has been updated to the new       !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: jc,kc
      integer :: km,kp,ic
      real    :: alre,udx3
      real    :: amm,acc,app,dxp,dxxvx

      alre=al

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vx,pr) &
!$OMP   SHARED(kmv,kpv,am3ck,ac3ck,ap3ck) &
!$OMP   SHARED(al,ga,ro,alre,dt,qcap) &
!$OMP   SHARED(udx3c,rhs,rux) &
!$OMP   PRIVATE(ic,jc,kc,km,kp) &
!$OMP   PRIVATE(amm,acc,app,udx3) &
!$OMP   PRIVATE(dxxvx,dxp)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      do kc=2,nxm
      km=kc-1
      kp=kc+1
      udx3 = al*udx3c(kc)
      amm=am3ck(kc)
      acc=ac3ck(kc)
      app=ap3ck(kc)

!   Second derivative in x-direction of vx
!
            dxxvx=vx(kp,jc,ic)*app+vx(kc,jc,ic)*acc+vx(km,jc,ic)*amm

!  component of grad(pr) along x direction
!
            dxp=(pr(kc,jc,ic)-pr(km,jc,ic))*udx3

!    Calculate right hand side of Eq. 5 (VO96)
!
            rhs(kc,jc,ic)=(ga*qcap(kc,jc,ic)+ro*rux(kc,jc,ic) &
            &             +alre*dxxvx-dxp)*dt 

!    Store the non-linear terms for the calculation of 
!    the next timestep

            rux(kc,jc,ic)=qcap(kc,jc,ic)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

!  Solve equation and update velocity

      call SolveImpEqnUpdate_X()

!  Set boundary conditions on the vertical velocity at top
!  and bottom plates. This seems necessary.

      vx(1,:,:)=0.0d0
      vx(nx,:,:)=0.0d0

   end subroutine ImplicitAndUpdateVX
!-------------------------------------------------------------------------------
   subroutine ImplicitAndUpdateVY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the y (horizontal) dimension        !
!     and call the implicit solver                        !
!     After this routine, the velocity field in y has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: kc,jmm,jc,ic
      integer :: kpp,kmm
      real    :: alre,udy
      real    :: amm,acc,app
      real    :: dyp,dxxvy

      alre=al
      udy=dy*al

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vy,pr) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
!$OMP   SHARED(dy,al,ga,ro,alre,dt,dph) &
!$OMP   SHARED(udy,udx3m,rhs,ruy) &
!$OMP   PRIVATE(ic,jc,kc,kmm,kpp,jmm) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(dyp,dxxvy)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      jmm=jc-1
      do kc=1,nxm
      kmm=kmv(kc)
      kpp=kpv(kc)
      amm=am3sk(kc)
      acc=ac3sk(kc)
      app=ap3sk(kc)


!   Second derivative in x-direction of vy
!
!
            dxxvy=vy(kpp,jc,ic)*app+vy(kc,jc,ic)*acc+vy(kmm,jc,ic)*amm

!   component of grad(pr) along y direction
!
            dyp=(pr(kc,jc,ic)-pr(kc,jmm,ic))*udy

!    Calculate right hand side of Eq. 5 (VO96)
!
            rhs(kc,jc,ic)=(ga*dph(kc,jc,ic)+ro*ruy(kc,jc,ic) &
            &             +alre*dxxvy-dyp)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

            ruy(kc,jc,ic)=dph(kc,jc,ic)
      enddo
      enddo
      enddo

!  Solve equation and update velocity

      call SolveImpEqnUpdate_YZ(vy,rhs)
      
   end subroutine ImplicitAndUpdateVY
!-------------------------------------------------------------------------------
   subroutine ImplicitAndUpdateVZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the z (horizontal) dimension        !
!     and call the implicit solver.                       !
!     After this routine, the velocity field in z has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: kc,jc,ic,imm
      integer :: kmm,kpp
      real    :: alre,amm,acc,app,udz
      real    :: dxxvz,dzp

      alre=al
      udz=dz*al

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vz,pr) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
!$OMP   SHARED(dz,al,ga,ro,alre,dt,dq) &
!$OMP   SHARED(udx3m,rhs,ruz) &
!$OMP   PRIVATE(ic,jc,kc,imm,kmm,kpp) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(dxxvz,dzp)
      do ic=xstart(3),xend(3)
      imm=ic-1
      do jc=xstart(2),xend(2)
      do kc=1,nxm
      kmm=kmv(kc)
      kpp=kpv(kc)
      amm=am3sk(kc)
      acc=ac3sk(kc)
      app=ap3sk(kc)

!   Second derivative in x-direction of vz
!
            dxxvz=vz(kpp,jc,ic)*app+vz(kc,jc,ic)*acc+vz(kmm,jc,ic)*amm
      
!   component of grad(pr) along z direction
!
            dzp=(pr(kc,jc,ic)-pr(kc,jc,imm))*dz*al

!    Calculate right hand side of Eq. 5 (VO96)
!
            rhs(kc,jc,ic)=(ga*dq(kc,jc,ic)+ro*ruz(kc,jc,ic) &
                          +alre*dxxvz-dzp)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

            ruz(kc,jc,ic)=dq(kc,jc,ic)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO


      call SolveImpEqnUpdate_YZ(vz,rhs)

   end subroutine ImplicitAndUpdateVZ
!-------------------------------------------------------------------------------
   subroutine SolveImpEqnUpdate_X
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any the vertical direction, and updates it to    !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none
      real, dimension(nx) :: amkl,apkl,ackl
      real :: amkT(nx-1),apkT(nx-1)
      real :: appk(nx-2)
      real :: ackT(nx)
      integer :: jc,kc,info,ic,nrhs
      integer :: ipkv(nx)
      real :: betadx,ackl_b

      betadx=beta*al

      amkl(1)=0.d0
      apkl(1)=0.d0
      ackl(1)=1.d0
      do kc=2,nxm
        ackl_b=1.0d0/(1.0d0-ac3ck(kc)*betadx)
        amkl(kc)=-am3ssk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3ssk(kc)*betadx*ackl_b
      enddo
      amkl(nx)=0.d0
      apkl(nx)=0.d0
      ackl(nx)=1.d0

      amkT=amkl(2:nx)
      apkT=apkl(1:(nx-1))
      ackT=ackl(1:nx)

      call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)

      nrhs=(xend(3)-xstart(3)+1)*(xend(2)-xstart(2)+1)
      do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=2,nxm
              ackl_b=1.0/(1.0-ac3ck(kc)*betadx)
              rhs(kc,jc,ic)=rhs(kc,jc,ic)*ackl_b
            end do
         end do
      end do

      call dgttrs('N',nx,nrhs,amkT,ackT,apkT,appk,ipkv,rhs,nx,info)

       do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=2,nxm
              vx(kc,jc,ic)=vx(kc,jc,ic) + rhs(kc,jc,ic)
             end do
          end do
      end do

   end subroutine SolveImpEqnUpdate_X
!-------------------------------------------------------------------------------
   subroutine SolveImpEqnUpdate_YZ(q,rhs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any of the horizontal directions, and updates    !
!     it to time t+dt                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      real,intent(inout) :: q(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo, &
     &      xstart(3)-lvlhalo:xend(3)+lvlhalo)
      real,intent(inout) :: rhs(1:nx,xstart(2):xend(2), &
     &      xstart(3):xend(3))
      real, dimension(nx) :: amkl,apkl,ackl
      integer :: jc,kc,info,ipkv(nxm),ic,nrhs
      real :: betadx,ackl_b
      real :: amkT(nxm-1),ackT(nxm),apkT(nxm-1),appk(nx-3)

      betadx=beta*al

      do kc=1,nxm
        ackl_b=1.0d0/(1.0d0-ac3sk(kc)*betadx)
        amkl(kc)=-am3sk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3sk(kc)*betadx*ackl_b
      enddo

      amkT=amkl(2:nxm)
      apkT=apkl(1:(nxm-1))
      ackT=ackl(1:nxm)

      call dgttrf(nxm,amkT,ackT,apkT,appk,ipkv,info)

      nrhs=(xend(3)-xstart(3)+1)*(xend(2)-xstart(2)+1)
       do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=1,nxm
               ackl_b=1.0/(1.0-ac3sk(kc)*betadx)
               rhs(kc,jc,ic)=rhs(kc,jc,ic)*ackl_b
             end do
          end do
      end do

      call dgttrs('N',nxm,nrhs,amkT,ackT,apkT,appk,ipkv,rhs,nx,info)

       do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=1,nxm
              q(kc,jc,ic)=q(kc,jc,ic) + rhs(kc,jc,ic)
             end do
          end do
      end do

   end subroutine SolveImpEqnUpdate_YZ
!-------------------------------------------------------------------------------
end module UpdateVel
