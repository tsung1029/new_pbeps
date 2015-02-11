!-----------------------------------------------------------------------
!
      module mprbpush2d
!
! Fortran90 interface to parallel 2d PIC Fortran77 library
! mprbpush2mod.f contains multi-tasking interface procedures to process
!                relativistic particles with magnetic fields:
!                defines module mprbpush2d
! rdjpost => impgrjpost2 deposits relativistic current density, with
!            various interpolations and optimizations.
!            calls MPGRJPOST2, MPGSRJPOST2, MPGSRJOST2X, MPGRJPOST2L,
!            MPGSRJPOST2L, MPGSRJOST2XL, MPGRJPOST22, MPGSRJPOST22,
!            MPGSRJOST22X, MPGRJPOST22L, MPGSRJPOST22L, or
!            MPGSRJOST22XL
! rpush => impgrpush2 push relativistic particles with 2 component
!          electric field, with various interpolations and optimizations
!          calls MPGRPUSH2, MPGSRPUSH2, MPGRPUSH2L, or MPGSRPUSH2L
! rpush => impgrbpush2 push relativistic particles with magnetic field
!          and 2 component electric field, with various interpolations
!          and optimizations.
!          calls MPGRBPUSH2, MPGSRBPUSH2, MPGRBPUSH2L, or MPGSRBPUSH2L
! rpush3 => impgrbpush23 push relativistic particles with magnetic field
!           and 3 component electric field, with various interpolations
!           and optimizations.
!           calls MPGRBPUSH23, MPGSRBPUSH23, MPGRBPUSH23L,
!           MPGSRBPUSH23L, MPGRBPUSH22, MPGSRBPUSH22, MPGRBPUSH22L,
!           or MPGSRBPUSH22L
! rpushzf => imprpush2zf push 2d relativistic particles with no forces.
!            calls MPRPUSH2ZF
! rpush3zf => imprpush23zf push 2-1/2d relativistic particles with no
!             forces.
!             calls MPRPUSH23ZF, or MPRPUSH2ZF
! rgcjpost => imprpush23zf deposits time-centered relativistic current
!             density with 2d electrostatic fields.
!             calls MPGRCJPOST2, or MPGRCJPOST2L
! mprbpush2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 12, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use prbpush2d, only: wtimer, retard, ipcptov2
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: rdjpost, rpush, rpush3, retard, ipcptov2
      public :: rpushzf, rpush3zf, rgcjpost
!
! buffer data for current deposit
      real, dimension(:,:,:,:,:), allocatable :: cup
      integer :: szbuf = 0
      save
!
! define interface to Fortran77 procedures
!
      interface
         subroutine MPGRJPOST2(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idimp&
     &,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJPOST2(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idim&
     &p,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJOST2X(part,cu,npp,nps,noff,nn,amxy,qm,dt,ci,nx&
     &,ny,idimp,npmax,nblok,nxv,nxvyp,npd,n27,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n27
         integer :: ipbc, nmt, ierr
         real :: qm, dt, ci
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
         subroutine MPGRJPOST2L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idim&
     &p,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJPOST2L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idi&
     &mp,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJOST2XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,ci,n&
     &x,ny,idimp,npmax,nblok,nxv,nxvyp,npd,n12,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n12
         integer :: ipbc, nmt, ierr
         real :: qm, dt, ci
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
         subroutine MPGRJPOST22(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idim&
     &p,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJPOST22(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idi&
     &mp,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJOST22X(part,cu,npp,nps,noff,nn,amxy,qm,dt,ci,n&
     &x,ny,idimp,npmax,nblok,nxv,nxvyp,npd,n18,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n18
         integer :: ipbc, nmt, ierr
         real :: qm, dt, ci
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
         subroutine MPGRJPOST22L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,idi&
     &mp,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJPOST22L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,id&
     &imp,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRJOST22XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,ci,&
     &nx,ny,idimp,npmax,nblok,nxv,nxvyp,npd,n8,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n8
         integer :: ipbc, nmt, ierr
         real :: qm, dt, ci
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
         subroutine MPGRPUSH2(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,i&
     &dimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRPUSH2(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny&
     &,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,e&
     &k,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci&
     &,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci&
     &,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci&
     &,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,c&
     &i,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRBPUSH22(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ci,e&
     &k,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRBPUSH22(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ci,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRBPUSH22L(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ci,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRBPUSH22L(part,fxy,bz,npp,nps,noff,qbm,dt,dtc,ci&
     &,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPRPUSH2ZF(part,npp,nps,dt,ci,ek,idimp,npmax,nblok,n&
     &x,ny,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nx, ny, ipbc, nmt, ierr
         real :: dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPRPUSH23ZF(part,npp,nps,dt,ci,ek,idimp,npmax,nblok,&
     &nx,ny,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nx, ny, ipbc, nmt, ierr
         real :: dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRCJPOST2(part,fxy,npp,nps,noff,cu,qm,qbm,dt,ci,id&
     &imp,npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRCJPOST2L(part,fxy,npp,nps,noff,cu,qm,qbm,dt,ci,i&
     &dimp,npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdjpost
         module procedure impgrjpost2
      end interface
!
      interface rpush
         module procedure impgrpush2
         module procedure impgrbpush2
      end interface
!
      interface rpush3
         module procedure impgrbpush23
      end interface
!
      interface rpushzf
         module procedure imprpush2zf
      end interface
!
      interface rpush3zf
         module procedure imprpush23zf
      end interface
!
      interface rgcjpost
         module procedure imprpush23zf
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impgrpush2(part,fxy,npp,noff,qbm,dt,ci,ek,tpush,nx,n&
     &y,ipbc,inorder,popt)
! multi-tasking relativistic particle push with 2d electrostatic fields,
! 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp
         integer :: nmt, order, opt, ierr
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
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tpush,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MPGSRPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny&
     &,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
            else
               call MPGRPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSRPUSH2(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
            else
               call MPGRPUSH2(part,fxy,npp,nps,noff,qbm,dt,ci,ek,nx,ny,i&
     &dimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tpush,dtime)
         end subroutine impgrpush2
!
         subroutine impgrjpost2(part,cu,npp,noff,qm,dt,ci,tdjpost,nx,ny,&
     &ipbc,inorder,djopt)
! multi-tasking, relativisitc current deposit, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, ci, tdjpost
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
         integer, dimension(n27,npd,size(part,3),ntasks+1) :: nn
         real, dimension(n27,npd,size(part,3),ntasks+1) :: amxy
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
         call wtimer(tdjpost,dtime,-1)
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSRJPOST22L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny&
     &,idimp,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MPGSRJOST22XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,&
     &ci,nx,ny,idimp,npmax,nblok,nxv,nxyp,npd,n8,ipbc,cup,idtask,nmt,ier&
     &r)
               else
                  call MPGRJPOST22L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,&
     &idimp,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSRJPOST22(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,&
     &idimp,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MPGSRJOST22X(part,cu,npp,nps,noff,nn,amxy,qm,dt,c&
     &i,nx,ny,idimp,npmax,nblok,nxv,nxyp,npd,n18,ipbc,cup,idtask,nmt,ier&
     &r)
                else
                  call MPGRJPOST22(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,i&
     &dimp,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSRJPOST2L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,&
     &idimp,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MPGSRJOST2XL(part,cu,npp,nps,noff,nn,amxy,qm,dt,c&
     &i,nx,ny,idimp,npmax,nblok,nxv,nxyp,npd,n12,ipbc,cup,idtask,nmt,ier&
     &r)
               else
                  call MPGRJPOST2L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,i&
     &dimp,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSRJPOST2(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,i&
     &dimp,npmax,nblok,nxv,nxyp,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MPGSRJOST2X(part,cu,npp,nps,noff,nn,amxy,qm,dt,ci&
     &,nx,ny,idimp,npmax,nblok,nxv,nxyp,npd,n27,ipbc,cup,idtask,nmt,ierr&
     &)
                else
                  call MPGRJPOST2(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,id&
     &imp,npmax,nblok,nxv,nypmx,ipbc,cup,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tdjpost,dtime)
         end subroutine impgrjpost2
!
         subroutine impgrbpush2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,t&
     &push,nx,ny,ipbc,inorder,popt)
! multi-tasking, relativistic particle push with 2d electromagnetic
! fields,1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp
         integer :: nmt, order, opt, ierr
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
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tpush,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MPGSRBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci&
     &,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
            else
               call MPGRBPUSH2L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSRBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
            else
               call MPGRBPUSH2(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,ci,e&
     &k,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tpush,dtime)
         end subroutine impgrbpush2
!
         subroutine impgrbpush23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,&
     &tpush,nx,ny,ipbc,inorder,popt)
! multi-tasking, relativistic particle push with 2-1/2d electromagnetic
! fields, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp
         integer :: nmt, order, opt, ierr
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
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tpush,dtime,-1)
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSRBPUSH22L(part,fxy,bxy,npp,nps,noff,qbm,dt,dt&
     &c,ci,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MPGRBPUSH22L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc&
     &,ci,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSRBPUSH22(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc&
     &,ci,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MPGRBPUSH22(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,&
     &ci,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSRBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dt&
     &c,ci,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MPGRBPUSH23L(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc&
     &,ci,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSRBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc&
     &,ci,ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MPGRBPUSH23(part,fxy,bxy,npp,nps,noff,qbm,dt,dtc,&
     &ci,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tpush,dtime)
         end subroutine impgrbpush23
!
         subroutine imprpush2zf(part,npp,dt,ci,ek,tpush,nx,ny,ipbc)
! multi-tasking relativistic particle push with no forces, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         real :: dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, nmt, ierr
         real :: tp
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
         real, dimension(ntasks) :: ekp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nmt = ntasks
! initialize timer
         call wtimer(tp,dtime,-1)
         call MPRPUSH2ZF(part,npp,nps,dt,ci,ek,idimp,npmax,nblok,nx,ny,i&
     &pbc,ekp,idtask,nmt,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine imprpush2zf
!
         subroutine imprpush23zf(part,npp,dt,ci,ek,tpush,nx,ny,ipbc,ndim&
     &)
! multi-tasking 2d relativistic particle push with no forces,
! 1d partition
         implicit none
         integer :: nx, ny, ipbc, ndim
         real :: dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, nmt, ierr
         real :: tp
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
         real, dimension(ntasks) :: ekp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nmt = ntasks
! initialize timer
         call wtimer(tp,dtime,-1)
         select case(ndim)
         case (2)
            call MPRPUSH2ZF(part,npp,nps,dt,ci,ek,idimp,npmax,nblok,nx,n&
     &y,ipbc,ekp,idtask,nmt,ierr)
         case (3)
            call MPRPUSH23ZF(part,npp,nps,dt,ci,ek,idimp,npmax,nblok,nx,&
     &ny,ipbc,ekp,idtask,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine imprpush23zf
!
         subroutine impgrcjpost2(part,fxy,npp,noff,cu,qm,qbm,dt,ci,tdcjp&
     &ost,inorder)
! multi-tasking relativitic current density deposit
! with 2d electrostatic fields, 1d partition
         implicit none
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, cu
         integer, dimension(:), pointer :: npp, noff
         integer, optional :: inorder
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, order, nmt, ierr
         integer :: nnxyp
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),size(cu,4),nta&
!    &sks) :: cup
         real :: tdc
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3)
         nnxyp = size(cu,1)*nxv*nypmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! check if size of buffer has changed
         if (szbuf < nnxyp*nblok) then
            if (szbuf /= 0) deallocate(cup)
! allocate buffer
            allocate(cup(size(cu,1),nxv,nypmx,nblok,ntasks))
            szbuf = nnxyp*nblok
         endif
! initialize timer
         call wtimer(tdc,dtime,-1)
         if (order==LINEAR) then
            call MPGRCJPOST2L(part,fxy,npp,nps,noff,cu,qm,qbm,dt,ci,idim&
     &p,npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         else
            call MPGRCJPOST2(part,fxy,npp,nps,noff,cu,qm,qbm,dt,ci,idimp&
     &,npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine impgrcjpost2
!

      end module mprbpush2d
