module explicit

   use param, only: dz, dy, dyq, dzq, BuoFac, ChemFac, opr, osc, &
       &            udx3c, nxm, kpv, kmv, udx3m
   use local_arrays, only: vy,vx,xi,vz,xiro,temp,hro,qcap,dph,dq
   use decomp_2d, only: xstart,xend

   implicit none

   private

   public :: ExplicitTermsComp, ExplicitTermsTemp, ExplicitTermsVX, &
   &         ExplicitTermsVY, ExplicitTermsVZ

contains

   subroutine ExplicitTermsComp

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
            dzzt=(xi(kc,jc,ip)-2.0d0*xi(kc,jc,ic)+xi(kc,jc,im))*udzq
      
!
!   yy second derivatives of xi
!
            dyyt=(xi(kc,jp,ic)-2.0d0*xi(kc,jc,ic)+xi(kc,jm,ic))*udyq
!
            xiro(kc,jc,ic) = -(htx+hty+htz)+dyyt+dzzt
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

   end subroutine ExplicitTermsComp
!------------------------------------------------------------------------------
   subroutine ExplicitTermsTemp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the temperature.                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: jc,kc,ic
      integer :: km,kp,jm,jp,im,ip
      real    :: htx,hty,htz,udy,udz
      real    :: udzq,udyq
      real    :: dyyt,dzzt

      udz=dz*0.25d0
      udy=dy*0.25d0
      udzq=dzq*opr
      udyq=dyq*opr

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,vz,vy,vx,nxm) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
!$OMP   SHARED(udy,udzq,udyq,udx3c,temp,hro) &
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
      htz=((vz(km,jc,ip)+vz(kc,jc,ip))*(temp(kc,jc,ip)+temp(kc,jc,ic))- &
           (vz(km,jc,ic)+vz(kc,jc,ic))*(temp(kc,jc,ic)+temp(kc,jc,im)) &
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
      hty=((vy(kc,jp,ic)+vy(km,jp,ic))*(temp(kc,jp,ic)+temp(kc,jc,ic))- &
           (vy(kc,jc,ic)+vy(km,jc,ic))*(temp(kc,jc,ic)+temp(kc,jm,ic)) &
          )*udy
!
!    rho vx term
!
!
!                 d  rho q_x 
!                -----------
!                 d   x      
!
      htx=((vx(kp,jc,ic)+vx(kc,jc,ic))*(temp(kp,jc,ic)+temp(kc,jc,ic))- &
           (vx(kc,jc,ic)+vx(km,jc,ic))*(temp(kc,jc,ic)+temp(km,jc,ic)) &
          )*udx3c(kc)*0.25d0
!
!
!   zz second derivatives of temp
!
            dzzt=(temp(kc,jc,ip)-2.0d0*temp(kc,jc,ic)+temp(kc,jc,im))*udzq
      
!
!   yy second derivatives of temp
!
            dyyt=(temp(kc,jp,ic)-2.0*temp(kc,jc,ic)+temp(kc,jm,ic))*udyq
!
            hro(kc,jc,ic) = -(htx+hty+htz)+dyyt+dzzt
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

   end subroutine ExplicitTermsTemp
!------------------------------------------------------------------------------
   subroutine ExplicitTermsVX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the x (vertical) dimension          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: jc,kc
      integer :: km,kp,jmm,jpp,ic,imm,ipp
      real    :: hxx,hxy,hxz
      real    :: udz,udy,tempit
      real    :: udzq,udyq
      real    :: dzzvx,dyyvx

      udy=dy*0.25d0
      udz=dz*0.25d0

      udyq=dyq
      udzq=dzq

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vz,vy,vx,dz,dy) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz,ChemFac) &
!$OMP   SHARED(udy,udzq,udyq,udx3c,qcap,temp,xi,BuoFac) &
!$OMP   PRIVATE(ic,jc,kc,imm,ipp,km,kp) &
!$OMP   PRIVATE(jmm,jpp,tempit) &
!$OMP   PRIVATE(hxx,hxy,hxz,dzzvx,dyyvx)
      do ic=xstart(3),xend(3)
       imm=ic-1
       ipp=ic+1
       do jc=xstart(2),xend(2)
        jmm=jc-1
        jpp=jc+1
        do kc=2,nxm
         km=kc-1
         kp=kc+1
!
!    vx vz term
!
!
!                d  q_x q_t 
!             -----------
!                d   t      
!
!
      hxz=(((vz(kc,jc,ipp)+vz(km,jc,ipp)) &
           *(vx(kc,jc,ipp)+vx(kc,jc,ic))) &
          -((vz(kc,jc,ic)+vz(km,jc,ic)) &
           *(vx(kc,jc,ic)+vx(kc,jc,imm))))*udz
!
!    vx vy term
!
!
!                d  q_x q_r 
!             -----------
!                d   r      
!
      hxy=(((vy(kc,jpp,ic)+vy(km,jpp,ic)) &
           *(vx(kc,jpp,ic)+vx(kc,jc,ic))) &
          -((vy(kc,jc,ic)+vy(km,jc,ic)) &
           *(vx(kc,jc,ic)+vx(kc,jmm,ic))))*udy
!
!    vx vx term
!
!
!                 d  q_x q_x 
!                -----------
!                 d   x      
!
      hxx=((vx(kp,jc,ic)+vx(kc,jc,ic))*(vx(kp,jc,ic)+vx(kc,jc,ic)) &
          -(vx(kc,jc,ic)+vx(km,jc,ic))*(vx(kc,jc,ic)+vx(km,jc,ic)) &
          )*udx3c(kc)*0.25d0
!
!  add the buoyancy term
!
          tempit=BuoFac*temp(kc,jc,ic)+ChemFac*xi(kc,jc,ic)

!
!   z second derivatives of vx
!
            dzzvx=(vx(kc,jc,imm)-2.0d0*vx(kc,jc,ic)+vx(kc,jc,ipp))*udzq
!
!   y second derivatives of vx
!
            dyyvx=(vx(kc,jmm,ic)-2.0d0*vx(kc,jc,ic)+vx(kc,jpp,ic))*udyq


          qcap(kc,jc,ic) =-(hxx+hxy+hxz)+dyyvx+dzzvx+tempit
            
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

   end subroutine ExplicitTermsVX
!------------------------------------------------------------------------------
   subroutine ExplicitTermsVY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the y (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
      integer :: kpp,kmm
      real    :: udzq,udyq
      real    :: udy,udz,hyx,hyy,hyz 
      real    :: dyyvy, dzzvy

      udyq=dyq
      udzq=dzq

      udy=dy*0.25d0
      udz=dz*0.25d0

!$OMP  PARALLEL DO &
!$OMP  DEFAULT(none) &
!$OMP  SHARED(xstart,xend,vz,vy,vx,dz,dy) &
!$OMP  SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
!$OMP  SHARED(udy,udzq,udyq,udx3m,dph,nxm) &
!$OMP  PRIVATE(ic,jc,kc,imm,ipp,kmm,kp,kpp) &
!$OMP  PRIVATE(jmm,jpp) &
!$OMP  PRIVATE(hyx,hyy,hyz,dyyvy,dzzvy)
      do ic=xstart(3),xend(3)
       imm=ic-1
       ipp=ic+1
       do jc=xstart(2),xend(2)
        jmm=jc-1
        jpp=jc+1
        do kc=1,nxm
         kmm=kmv(kc)
         kpp=kpv(kc)
         kp=kc+1
!
!     vy vx term
!
!                 d  q_x q_r 
!                -----------
!                 d   x      
!
      hyx=((vx(kp,jc,ic)+vx(kp,jmm,ic))*(vy(kpp,jc,ic)+vy(kc,jc,ic)) &
          -(vx(kc,jc,ic)+vx(kc,jmm,ic))*(vy(kc,jc,ic)+vy(kmm,jc,ic)) &
          )*udx3m(kc)*0.25d0

!     
!     vy vy term
!
!                 d  q_r q_r 
!                ------------
!                 d   r      
!
      hyy=( (vy(kc,jpp,ic)+vy(kc,jc,ic)) &
           *(vy(kc,jpp,ic)+vy(kc,jc,ic)) &
           -(vy(kc,jmm,ic)+vy(kc,jc,ic)) &
           *(vy(kc,jmm,ic)+vy(kc,jc,ic)) &
          )*udy
!
!     vz vy term
!
!                 d  q_t q_r 
!                ------------
!                 d   t      
!
      hyz=( (vy(kc,jc,ipp)+vy(kc,jc,ic)) &
           *(vz(kc,jc,ipp)+vz(kc,jmm,ipp)) &
           -(vy(kc,jc,ic)+vy(kc,jc,imm)) &
           *(vz(kc,jc,ic)+vz(kc,jmm,ic)) &
          )*udz

!
!   yy second derivative of vy
!
            dyyvy=(vy(kc,jpp,ic)-2.0d0*vy(kc,jc,ic)+vy(kc,jmm,ic))*udyq

!
!   zz second derivative of vy
!
            dzzvy=(vy(kc,jc,ipp)-2.0d0*vy(kc,jc,ic)+vy(kc,jc,imm))*udzq


            dph(kc,jc,ic)=-(hyx+hyy+hyz)+dyyvy+dzzvy
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

   end subroutine ExplicitTermsVY
!------------------------------------------------------------------------------
   subroutine ExplicitTermsVZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the z (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
      integer :: kmm,kpp
      real    :: hzx,hzy,hzz,udy,udz
      real    :: udyq,udzq
      real    :: dzzvz,dyyvz

!
      udyq=dyq
      udzq=dzq

      udy=dy*0.25d0
      udz=dz*0.25d0

!$OMP  PARALLEL DO &
!$OMP  DEFAULT(none) &
!$OMP  SHARED(xstart,xend,vz,vy,vx,dz,dy,udx3m) &
!$OMP  SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
!$OMP  SHARED(udy,udzq,udyq,dq,nxm) &
!$OMP  PRIVATE(ic,jc,kc,imm,ipp,kmm,kp,kpp) &
!$OMP  PRIVATE(jmm,jpp) &
!$OMP  PRIVATE(hzz,hzy,hzx,dzzvz,dyyvz)

      do ic=xstart(3),xend(3)
       imm=ic-1
       ipp=ic+1
       do jc=xstart(2),xend(2)
        jmm=jc-1
        jpp=jc+1
        do kc=1,nxm
         kmm=kmv(kc)
         kpp=kpv(kc)
         kp=kc+1
      
!     vz vz term
!
!
!                 d  q_t q_t 
!                ------------
!                 d   t      
!
      hzz=( (vz(kc,jc,ipp)+vz(kc,jc,ic)) &
           *(vz(kc,jc,ipp)+vz(kc,jc,ic)) &
           -(vz(kc,jc,imm)+vz(kc,jc,ic)) &
           *(vz(kc,jc,imm)+vz(kc,jc,ic)) &
          )*udz

!     vz vy term
!
!
!                 d  q_t q_r 
!                ------------
!                 d   r      
!
      hzy=( (vy(kc,jpp,ic)+vy(kc,jpp,imm)) &
           *(vz(kc,jpp,ic)+vz(kc,jc,ic)) &
           -(vy(kc,jc,ic)+vy(kc,jc,imm)) &
           *(vz(kc,jc,ic)+vz(kc,jmm,ic)) &
          )*udy
!
!     vz vx term
!
!
!                 d  q_t q_x 
!                -----------
!                 d   x      
!
      hzx=((vx(kp,jc,ic)+vx(kp,jc,imm))*(vz(kpp,jc,ic)+vz(kc,jc,ic)) &
          -(vx(kc,jc,ic)+vx(kc,jc,imm))*(vz(kc,jc,ic)+vz(kmm,jc,ic)) &
          )*udx3m(kc)*0.25d0
!
!
!
!   11 second derivative of vz
!
            dzzvz=(vz(kc,jc,ipp)-2.0d0*vz(kc,jc,ic)+vz(kc,jc,imm))*udzq
!
!   22 second derivative of vz
!
            dyyvz=(vz(kc,jpp,ic)-2.0d0*vz(kc,jc,ic)+vz(kc,jmm,ic))*udyq

!
        dq(kc,jc,ic)=-(hzx+hzy+hzz)+dyyvz+dzzvz
!
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

   end subroutine ExplicitTermsVZ
!------------------------------------------------------------------------------
end module explicit
