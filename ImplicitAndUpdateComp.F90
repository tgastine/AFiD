!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateComp.F90                      !
!    CONTAINS: subroutine ImplicitAndUpdateComp           !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the composition and call the implicit solver.       !
!     After this routine, the composition has been        !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateComp
      use param
      use local_arrays, only: xi,xiro,ruxi,rhs
      use decomp_2d, only: xstart,xend
      implicit none
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
                    +xi(kc  ,jc,ic)*ac3ck(kc) &
                    +xi(kc-1,jc,ic)*am3ck(kc)


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

      call SolveImpEqnUpdate_Comp

!  Set boundary conditions on the composition field at top
!  and bottom plates. This seems necessary.

       xi(1,xstart(2):xend(2),xstart(3):xend(3)) &
          = xibp(xstart(2):xend(2),xstart(3):xend(3))

       xi(nx,xstart(2):xend(2),xstart(3):xend(3)) &
          = xitp(xstart(2):xend(2),xstart(3):xend(3))


      return
      end
