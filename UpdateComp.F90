module UpdateComp

   use param
   use local_arrays, only: xi,xiro,ruxi,rhs
   use decomp_2d, only: xstart,xend

   implicit none

   private

   public :: ImplicitAndUpdateComp

contains

   subroutine ImplicitAndUpdateComp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the composition and call the implicit solver.       !
!     After this routine, the composition has been        !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: jc,kc,ic
      integer :: km,kp
      real    :: alpec,dxxt
      real    :: app,acc,amm

      alpec=al*osc

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,xi) &
!$OMP   SHARED(kmv,kpv,am3ck,ac3ck,ap3ck) &
!$OMP   SHARED(ga,ro,alpec,dt) &
!$OMP   SHARED(rhs,ruxi,xiro) &
!$OMP   PRIVATE(ic,jc,kc,km,kp) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(dxxt)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      do kc=2,nxm

!   Calculate second derivative of composition in the x-direction.
!   This is the only term calculated implicitly for composition.

               dxxt= xi(kc+1,jc,ic)*ap3ck(kc) &
               &    +xi(kc  ,jc,ic)*ac3ck(kc) &
               &    +xi(kc-1,jc,ic)*am3ck(kc)


!    Calculate right hand side of Eq. 5 (VO96)

            rhs(kc,jc,ic)=(ga*xiro(kc,jc,ic)+ro*ruxi(kc,jc,ic)+alpec*dxxt)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

            ruxi(kc,jc,ic)=xiro(kc,jc,ic)

        enddo
       enddo
      enddo
!$OMP END PARALLEL DO


!  Solve equation and update composition

      call SolveImpEqnUpdate_Comp()

!  Set boundary conditions on the composition field at top
!  and bottom plates. This seems necessary.

       xi(1,xstart(2):xend(2),xstart(3):xend(3)) &
       &  = xibp(xstart(2):xend(2),xstart(3):xend(3))

       xi(nx,xstart(2):xend(2),xstart(3):xend(3)) &
       &  = xitp(xstart(2):xend(2),xstart(3):xend(3))

   end subroutine ImplicitAndUpdateComp
!-----------------------------------------------------------------------------
   subroutine SolveImpEqnUpdate_Comp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     composition, and updates it to time t+dt            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension(nx) :: amkl,apkl,ackl
      integer :: jc,kc,info,ipkv(nx),ic,nrhs
      real :: betadx,ackl_b
      real :: amkT(nx-1),ackT(nx),apkT(nx-1),appk(nx-2)

!     Calculate the coefficients of the tridiagonal matrix
!     The coefficients are normalized to prevent floating
!     point errors.

      betadx=0.5d0*al*dt*osc

      amkl(1)=0.d0
      apkl(1)=0.d0
      ackl(1)=1.d0
      do kc=2,nxm
        ackl_b=1.0d0/(1.-ac3ssk(kc)*betadx)
        amkl(kc)=-am3ssk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3ssk(kc)*betadx*ackl_b
      enddo
      amkl(nx)=0.d0
      apkl(nx)=0.d0
      ackl(nx)=1.d0

      amkT=amkl(2:nx)
      apkT=apkl(1:nxm)
      ackT=ackl(1:nx)

!     Call to LAPACK library to factor tridiagonal matrix.
!     No solving is done in this call.

      call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)
     
      nrhs=(xend(3)-xstart(3)+1)*(xend(2)-xstart(2)+1)
      do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
           do kc=2,nxm
              ackl_b=1.0/(1.0-ac3ssk(kc)*betadx)
              rhs(kc,jc,ic)=rhs(kc,jc,ic)*ackl_b
           end do
        end do
      end do
      
      call dgttrs('N',nx,nrhs,amkT,ackT,apkT,appk,ipkv,rhs,nx,info)

       do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=2,nxm
              xi(kc,jc,ic)=xi(kc,jc,ic) + rhs(kc,jc,ic)
             end do
          end do
      end do

   end subroutine SolveImpEqnUpdate_Comp
!-----------------------------------------------------------------------------
end module UpdateComp
