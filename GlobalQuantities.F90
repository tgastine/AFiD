!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: GlobalQuantities.F90                           !
!    CONTAINS: subroutine GlobalQuantities                !
!                                                         ! 
!    PURPOSE: Calculate maximum velocity and temperature. !
!     volume averaged Nusselt number and Reynolds number  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine GlobalQuantities
      use param
      use local_arrays, only: vy,vx,vz,temp,xi
      use decomp_2d, only: xstart,xend
      use mpih
      implicit none
      integer :: jc,kc,kp,ic
      real :: anusin,vol,vxcen,fac2,tempcen,xicen,xinusin
      real :: vy_rms_vol,vz_rms_vol
      real :: vx_rms_vol,vzvyvx_rms_vol

!EP   Initialize
      vmax(1)=-huge(0.0d0)
      vmax(2)=-huge(0.0d0)
      vmax(3)=-huge(0.0d0)
      tempmax=-huge(0.0d0)
      tempmin=huge(0.0d0)
      tempm=0.0d0
      ximax=-huge(0.0d0)
      ximin=huge(0.0d0)
      xim=0.0d0
      anusin=0.d0 
      vx_rms_vol = 0.0d0
      vy_rms_vol = 0.0d0
      vz_rms_vol = 0.0d0
      vzvyvx_rms_vol = 0.0d0
      vmax = 0.0d0
      vol = 1.d0/(alx3*dx*real(nzm)*real(nym))


!EP   Loop over volume
      do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
          do kc=1,nxm
          kp = kc + 1
          fac2 = g3rc(kc)
          vmax(1) = max(vmax(1),abs(vz(kc,jc,ic)))
          vmax(2) = max(vmax(2),abs(vy(kc,jc,ic)))
          vmax(3) = max(vmax(3),abs(vx(kc,jc,ic)))
          tempmax = max(tempmax,temp(kc,jc,ic))
          tempmin = min(tempmin,temp(kc,jc,ic))
          ximax = max(ximax,xi(kc,jc,ic))
          ximin = min(ximin,xi(kc,jc,ic))
          vxcen = (vx(kc,jc,ic)+vx(kp,jc,ic))*0.5d0
          tempcen = (temp(kc,jc,ic)+temp(kp,jc,ic))*0.5d0
          xicen = (xi(kc,jc,ic)+xi(kp,jc,ic))*0.5d0
          anusin=anusin+tempcen*vxcen*fac2
          xinusin=anusin+xicen*vxcen*fac2
          tempm=tempm+tempcen*fac2
          xim=xim+xicen*fac2
          vx_rms_vol = vx_rms_vol + fac2*vx(kc,jc,ic)**2
          vy_rms_vol = vy_rms_vol + fac2*vy(kc,jc,ic)**2
          vz_rms_vol = vz_rms_vol + fac2*vz(kc,jc,ic)**2
          vzvyvx_rms_vol = vzvyvx_rms_vol + fac2* &
     &    (vz(kc,jc,ic)**2+vy(kc,jc,ic)**2+vx(kc,jc,ic)**2)  
          enddo
        enddo
      enddo

!EP   Reduce

      call MpiSumRealScalar(tempm)
      call MpiSumRealScalar(xim)
      call MpiSumRealScalar(anusin)
      call MpiSumRealScalar(xinusin)
      call MpiSumRealScalar(vx_rms_vol)
      call MpiSumRealScalar(vy_rms_vol)
      call MpiSumRealScalar(vz_rms_vol)
      call MpiSumRealScalar(vzvyvx_rms_vol)
      call MpiMinRealScalar(tempmin)
      call MpiMaxRealScalar(tempmax)
      call MpiMinRealScalar(ximin)
      call MpiMaxRealScalar(ximax)
      call MpiMaxRealScalar(vmax(1))
      call MpiMaxRealScalar(vmax(2))
      call MpiMaxRealScalar(vmax(3))
       
!EP   Write logs
      if(ismaster) then

      anusin=1.d0 + pra*anusin*vol
      xinusin=1.d0 + sc*xinusin*vol

      open(95,file='nu_vol.out',status='unknown',access='sequential', &
        position='append')
      write(95,'(es20.12,2es16.8)') time, anusin, xinusin
      close(95)

      tempm=tempm*vol
      xim=xim*vol
      vx_rms_vol=dsqrt(vx_rms_vol*vol)
      vy_rms_vol=dsqrt(vy_rms_vol*vol)
      vz_rms_vol=dsqrt(vz_rms_vol*vol)
      vzvyvx_rms_vol=dsqrt(vzvyvx_rms_vol*vol)

       open(94,file='rms_vel.out',status='unknown',position='append', &
        access='sequential')
       write(94,'(es20.12,4es16.8)') time,vz_rms_vol,vy_rms_vol,vx_rms_vol,vzvyvx_rms_vol
       close(94)

      endif

      return   
      end     
