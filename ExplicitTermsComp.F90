!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsComp.F90                          !
!    CONTAINS: subroutine ExplicitTermsTemp               !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the chemical composition.                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ExplicitTermsComp
      use param
      use local_arrays, only: vy,vx,xi,vz,xiro
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc,ic
      integer :: km,kp,jm,jp,im,ip
      real    :: htx,hty,htz,udy,udz
      real    :: udzq,udyq
      real    :: dyyt,dzzt

      udz=dz*0.25
      udy=dy*0.25
      udzq=dzq*osc
      udyq=dyq*osc

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,vz,vy,vx,nxm) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
!$OMP   SHARED(udy,udzq,udyq,udx3c,xi,xiro) &
!$OMP   PRIVATE(ic,jc,kc,im,ip,km,kp,jm,jp) &
!$OMP   PRIVATE(htx,hty,htz,dyyt,dzzt)
      do ic=xstart(3),xend(3)
       im=ic-1
       ip=ic+1
       do jc=xstart(2),xend(2)
        jm=jc-1
        jp=jc+1
        do kc=2,nxm
         km=kc-1
         kp=kc+1
!
!
!    rho vz term
!
!
!                d  rho q_z
!             -----------
!                d   z      
!
      htz=((vz(km,jc,ip)+vz(kc,jc,ip))*(xi(kc,jc,ip)+xi(kc,jc,ic))- &
           (vz(km,jc,ic)+vz(kc,jc,ic))*(xi(kc,jc,ic)+xi(kc,jc,im)) &
          )*udz
!
!
!    rho vy term
!
!
!                d  rho q_y 
!             -----------
!                d   y      
!
      hty=((vy(kc,jp,ic)+vy(km,jp,ic))*(xi(kc,jp,ic)+xi(kc,jc,ic))- &
           (vy(kc,jc,ic)+vy(km,jc,ic))*(xi(kc,jc,ic)+xi(kc,jm,ic)) &
          )*udy
!
!    rho vx term
!
!
!                 d  rho q_x 
!                -----------
!                 d   x      
!
      htx=((vx(kp,jc,ic)+vx(kc,jc,ic))*(xi(kp,jc,ic)+xi(kc,jc,ic))- &
           (vx(kc,jc,ic)+vx(km,jc,ic))*(xi(kc,jc,ic)+xi(km,jc,ic)) &
          )*udx3c(kc)*0.25d0
!
!
!   zz second derivatives of xi
!
            dzzt=(xi(kc,jc,ip) &
             -2.0*xi(kc,jc,ic) &
                 +xi(kc,jc,im))*udzq
      
!
!   yy second derivatives of xi
!
            dyyt=(xi(kc,jp,ic) &
             -2.0*xi(kc,jc,ic) &
                 +xi(kc,jm,ic))*udyq
!
            xiro(kc,jc,ic) = -(htx+hty+htz)+dyyt+dzzt
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
