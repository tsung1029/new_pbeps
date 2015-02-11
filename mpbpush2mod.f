!-----------------------------------------------------------------------
!
      module mpbpush2d
!
! Fortran90 interface to parallel 2d PIC Fortran77 library mpbpush2lib.f
! mpbpush2mod.f contains multi-tasking interface procedures to process
!              particles with magnetic fields:
!              defines module mpbpush2d
! djpost => impgjpost2 deposits current density, with various
!           interpolations and optimizations.
!           calls MPGJPOST2, MPGSJPOST2, MPGSJOST2X, MPGJPOST2L,
!           MPGSJPOST2L, MPGSJOST2XL, MPGJPOST22, MPGSJPOST22,
!           MPGSJOST22X, MPGJPOST22L, MPGSJPOST22L, or MPGSJOST22XL
! push => impgbpush2 push particles with magnetic field and 2 component
!         electric field, with various interpolations and optimizations.
!         calls MPGBPUSH2, MPGSBPUSH2, MPGBPUSH2L, or MPGSBPUSH2L
! push3 => impgbpush23 push particles with magnetic field and 3
!          component electric field, with various interpolations and
!          optimizations
!          calls MPGBPUSH23, MPGSBPUSH23, MPGBPUSH23L, MPGSBPUSH23L,
!          MPGBPUSH22, MPGSBPUSH22, MPGBPUSH22L, or MPGSBPUSH22L
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 12, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use pbpush2d, only: wtimer, retard
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: djpost, push, push3, retard
!
! define interface to Fortran77 procedures
!
      interface
         subroutine MPJDOST2(part,cux,cuy,cuz,npp,nps,noff,qm,dt,nx,idim&
     &p,npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: cux, cuy, cuz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGJPOST2(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJPOST2(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,np&
     &max,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSJOST2X(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,idim&
     &p,npmax,nblok,nxv,nxvyp,npd,n27,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nxvyp, npd, n27, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(n27,npd,nblok,nmt+1) :: nn
         real, dimension(n27,npd,nblok,nmt+1) :: amxy
         real, dimension(3*nxvyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJOST2X(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,ny,&
     &idimp,npmax,nblok,nxv,nxvyp,npd,n27,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n27
         integer :: ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(n27,npd,nblok,nmt+1) :: nn
         real, dimension(n27,npd,nblok,nmt+1) :: amxy
         real, dimension(3*nxvyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPJDOST2L(part,cux,cuy,cuz,npp,nps,noff,qm,dt,nx,idi&
     &mp,npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: cux, cuy, cuz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGJPOST2L(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,np&
     &max,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJPOST2L(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,n&
     &pmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSJOST2XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,idi&
     &mp,npmax,nblok,nxv,nxvyp,npd,n12,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nxvyp, npd, n12, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(n12,npd,nblok,nmt+1) :: nn
         real, dimension(n12,npd,nblok,nmt+1) :: amxy
         real, dimension(3*nxvyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJOST2XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,ny&
     &,idimp,npmax,nblok,nxv,nxvyp,npd,n12,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n12
         integer :: ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(n12,npd,nblok,nmt+1) :: nn
         real, dimension(n12,npd,nblok,nmt+1) :: amxy
         real, dimension(3*nxvyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGJPOST22(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,np&
     &max,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJPOST22(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,n&
     &pmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJOST22X(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,ny&
     &,idimp,npmax,nblok,nxv,nxvyp,npd,n18,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n18
         integer :: ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(n18,npd,nblok,nmt+1) :: nn
         real, dimension(n18,npd,nblok,nmt+1) :: amxy
         real, dimension(2*nxvyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGJPOST22L(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,n&
     &pmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJPOST22L(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,&
     &npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSJOST22XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,n&
     &y,idimp,npmax,nblok,nxv,nxvyp,npd,n8,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n8
         integer :: ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(n8,npd,nblok,nmt+1) :: nn
         real, dimension(n8,npd,nblok,nmt+1) :: amxy
         real, dimension(2*nxvyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPBPUSH2(part,fx,fy,bx,by,bz,npp,nps,noff,qbm,dt,ek,&
     &nx,idimp,npmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy, bx, by, bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPBPUSH2L(part,fx,fy,bx,by,bz,npp,nps,noff,qbm,dt,ek&
     &,nx,idimp,npmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy, bx, by, bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,&
     &nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPBPUSH2CQ(part,fx,fy,bx,by,bz,npp,nps,noff,qbm,dt,e&
     &k,nx,idimp,npmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy, bx, by, bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH2CQ(part,fxy,bxy,npp,nps,noff,qbm,dt,ek,nx,n&
     &y,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPBPUSH2CL(part,fx,fy,bx,by,bz,npp,nps,noff,qbm,dt,e&
     &k,nx,idimp,npmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy, bx, by, bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH2CL(part,fxy,bxy,npp,nps,noff,qbm,dt,ek,nx,n&
     &y,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,&
     &nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,&
     &nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek&
     &,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH22(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSBPUSH22(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGBPUSH22L(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSBPUSH22L(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ek,&
     &nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface djpost
         module procedure impgjpost2
      end interface
!
      interface push
         module procedure impgbpush2
!        module procedure impgbpush2c
      end interface
!
      interface push3
         module procedure impgbpush23
      end interface
!
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impgjpost2(part,cu,npp,noff,qm,dt,tdjpost,nx,ny,ipbc&
     &,inorder,djopt)
! multi-tasking current deposit, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, tdjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: npp, noff
! local data
         integer:: idimp, npmax, nblok, nxv, nypmx, nxyp, nnxyp
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128
         integer, parameter :: n8 = 8, n12 = 12, n18 = 18, n27 = 27
         integer, dimension(ntasks) :: idtask
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),size(cu,4),nta&
!    &sks) :: cup
         real, dimension(:,:,:,:,:), allocatable, save :: cup
         integer, save :: szbuf = 0
         integer, dimension(n27,npd,size(part,3),ntasks+1) :: nn
         real, dimension(n27,npd,size(part,3),ntasks+1) :: amxy
         real :: tj
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(cu,2); nypmx = size(cu,3); nxyp = nxv*nypmx
         nnxyp = size(cu,1)*nxyp
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! check if size of buffer has changed
         if (szbuf < nnxyp*nblok) then
            if (szbuf /= 0) deallocate(cup)
! allocate buffer
            allocate(cup(size(cu,1),nxv,nypmx,nblok,ntasks))
            szbuf = nnxyp*nblok
         endif
! initialize timer
         call wtimer(tj,dtime,-1)
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSJPOST22L(part,cu,npp,nps,noff,qm,dt,nx,ny,idi&
     &mp,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MPGSJOST22XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,npd,n8,ipbc,cup,idtask,nmt,ierr)
               else
                  call MPGJPOST22L(part,cu,npp,nps,noff,qm,dt,nx,ny,idim&
     &p,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSJPOST22(part,cu,npp,nps,noff,qm,dt,nx,ny,idim&
     &p,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MPGSJOST22X(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx&
     &,ny,idimp,npmax,nblok,nxv,nxyp,npd,n18,ipbc,cup,idtask,nmt,ierr)
               else
                  call MPGJPOST22(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp&
     &,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSJPOST2L(part,cu,npp,nps,noff,qm,dt,nx,ny,idim&
     &p,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MPGSJOST2XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx&
     &,ny,idimp,npmax,nblok,nxv,nxyp,npd,n12,ipbc,cup,idtask,nmt,ierr)
               else
                  call MPGJPOST2L(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp&
     &,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSJPOST2(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp&
     &,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MPGSJOST2X(part,cu,npp,nps,noff,nn,amxy,qm,dt,nx,&
     &ny,idimp,npmax,nblok,nxv,nxyp,npd,n27,ipbc,cup,idtask,nmt,ierr)
               else
                  call MPGJPOST2(part,cu,npp,nps,noff,qm,dt,nx,ny,idimp,&
     &npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tj,dtime)
         tdjpost = tdjpost + tj
         end subroutine impgjpost2
!
         subroutine impgbpush2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,tpush&
     &,nx,ny,ipbc,inorder,popt)
! multi-tasking particle push with 2d electromagnetic fields,
! 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
         real, dimension(ntasks) :: ekp
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MPGSBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,&
     &nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
            else
               call MPGBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
            else
               call MPGBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine impgbpush2
!
         subroutine impgbpush2cq(part,fxy,bxy,npp,noff,qbm,dt,ek,tpush,n&
     &x,ny,ipbc,inorder)
! multi-tasking particle push with 2d electromagnetic fields,
! 1d partition, with correction to Boris Mover
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm, dt, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, nmt, order
         integer :: ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
         real, dimension(ntasks) :: ekp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tpush,dtime,-1)
         if (order==LINEAR) then
            call MPGBPUSH2CL(part,fxy,bxy,npp,nps,noff,qbm,dt,ek,nx,ny,i&
     &dimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         else
            call MPGBPUSH2CQ(part,fxy,bxy,npp,nps,noff,qbm,dt,ek,nx,ny,i&
     &dimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tpush,dtime)
         end subroutine impgbpush2cq
!
         subroutine impgbpush23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,tpus&
     &h,nx,ny,ipbc,inorder,popt)
! multi-tasking particle push with 2-1/2d electromagnetic fields,
! 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
         real, dimension(ntasks) :: ekp
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSBPUSH22L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc&
     &,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MPGBPUSH22L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSBPUSH22(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MPGBPUSH22(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,e&
     &k,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc&
     &,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MPGBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MPGBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,e&
     &k,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine impgbpush23
!
      end module mpbpush2d
