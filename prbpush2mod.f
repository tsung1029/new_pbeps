!-----------------------------------------------------------------------
!
      module prbpush2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library prbpush2lib.f
! prbpush2mod.f contains interface procedures to process relativistic
!               particles with magnetic fields:
!               defines module prbpush2d
! rdjpost => ipgrjpost2 deposits relativistic current density, with
!            various interpolations and optimizations.
!            calls PGRJPOST2, PGSRJPOST2, PGSRJOST2X, PGRJPOST2L,
!            PGSRJPOST2L, PGSRJOST2XL, PGRJPOST22, PGSRJPOST22,
!            PGSRJOST22X, PGRJPOST22L, PGSRJPOST22L, or PGSRJOST22XL
! rpush => ipgrpush2 push relativistic particles with 2 component
!          electric field, with various interpolations and optimizations
!          calls PGRPUSH2, PGSRPUSH2, PGRPUSH2L, or PGSRPUSH2L
! rpush => ipgrbpush2 push relativistic particles with magnetic field
!          and 2 component electric field, with various interpolations
!          and optimizations.
!          calls PGRBPUSH2, PGSRBPUSH2, PGRBPUSH2L, or PGSRBPUSH2L
! rpush3 => ipgrbpush23 push relativistic particles with magnetic field
!           and 3 component electric field, with various interpolations
!           and optimizations.
!           calls PGRBPUSH23, PGSRBPUSH23, PGRBPUSH23L, PGSRBPUSH23L,
!           PGRBPUSH22, PGSRBPUSH22, PGRBPUSH22L, or PGSRBPUSH22L
! rpushzf => iprpush2zf push 2d relativistic particles with no forces.
!            calls PRPUSH2ZF
! rpush3zf => iprpush23zf push 2-1/2d relativistic particles with no
!             forces.
!             calls PRPUSH23ZF, or PRPUSH2ZF
! retard => iprretard2 retard relativistic particle position a half
!           time-step.
!           calls PRRETARD2, or PRRETARD22
! icptov2 converts momentum to velocity for relativistic particles.
!         calls PCPTOV2, or PCPTOV22
! rgcjpost => ipgrcjpost2 deposits time-centered relativistic current
!             density with 2d electrostatic fields.
!             calls PGRCJPOST2, or PGRCJPOST2L
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: september 24, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer
      public :: rdjpost, rpush, rpush3, retard, ipcptov2
      public :: rpushzf, rpush3zf, rgcjpost
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PGRJPOST2(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,npma&
     &x,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRJPOST2(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,npm&
     &ax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRJOST2X(part,cu,npp,noff,nn,amxy,qm,dt,ci,nx,ny,i&
     &dimp,npmax,nblok,nxv,nxvyp,npd,n27,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n27
         integer :: ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n27,npd,nblok) :: nn
         real, dimension(n27,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGRJPOST2L(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRJPOST2L(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,np&
     &max,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRJOST2XL(part,cu,npp,noff,nn,amxy,qm,dt,ci,nx,ny,&
     &idimp,npmax,nblok,nxv,nxvyp,npd,n12,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n12
         integer :: ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n12,npd,nblok) :: nn
         real, dimension(n12,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGRJPOST22(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRJPOST22(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,np&
     &max,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRJOST22X(part,cu,npp,noff,nn,amxy,qm,dt,ci,nx,ny,&
     &idimp,npmax,nblok,nxv,nxvyp,npd,n18,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n18
         integer :: ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n18,npd,nblok) :: nn
         real, dimension(n18,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGRJPOST22L(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,np&
     &max,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRJPOST22L(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,n&
     &pmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRJOST22XL(part,cu,npp,noff,nn,amxy,qm,dt,ci,nx,ny&
     &,idimp,npmax,nblok,nxv,nxvyp,npd,n8,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n8
         integer :: ipbc
         real :: qm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n8,npd,nblok) :: nn
         real, dimension(n8,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGRPUSH2(part,fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idimp,&
     &npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRPUSH2(part,fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idimp&
     &,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRPUSH2L(part,fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idimp&
     &,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRPUSH2L(part,fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idim&
     &p,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx,&
     &ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,&
     &nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBPUSH22(part,fxy,bz,npp,noff,qbm,dt,dtc,ci,ek,nx,&
     &ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRBPUSH22(part,fxy,bz,npp,noff,qbm,dt,dtc,ci,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBPUSH22L(part,fxy,bz,npp,noff,qbm,dt,dtc,ci,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRBPUSH22L(part,fxy,bz,npp,noff,qbm,dt,dtc,ci,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PRRETARD2(part,npp,dtc,ci,nx,ny,idimp,npmax,nblok,ip&
     &bc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, ipbc
         real :: dtc, ci
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PCPTOV2(part,npp,ci,idimp,npmax,nblok)
         integer :: idimp, npmax, nblok
         real :: ci
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PCPTOV22(part,npp,ci,idimp,npmax,nblok)
         integer :: idimp, npmax, nblok
         real :: ci
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PRPUSH2ZF(part,npp,dt,ci,ek,idimp,npmax,nblok,nx,ny,&
     &ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, ipbc
         real :: dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PRPUSH23ZF(part,npp,dt,ci,ek,idimp,npmax,nblok,nx,ny&
     &,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, ipbc
         real :: dt, ci, ek
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PGRCJPOST2(part,fxy,npp,noff,cu,qm,qbm,dt,ci,idimp,n&
     &pmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRCJPOST2L(part,fxy,npp,noff,cu,qm,qbm,dt,ci,idimp,&
     &npmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdjpost
         module procedure ipgrjpost2
      end interface
!
      interface rpush
         module procedure ipgrpush2
         module procedure ipgrbpush2
      end interface
!
      interface rpush3
         module procedure ipgrbpush23
      end interface
!
      interface rpushzf
         module procedure iprpush2zf
      end interface
!
      interface rpush3zf
         module procedure iprpush23zf
      end interface
!
      interface retard
         module procedure iprretard2
      end interface
!
      interface rgcjpost
         module procedure ipgrcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!
         subroutine ipgrpush2(part,fxy,npp,noff,qbm,dt,ci,ek,tpush,nx,ny&
     &,ipbc,inorder,popt)
! push relativistic particles with 2d electrostatic fields, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, order, opt
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSRPUSH2L(part,fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idim&
     &p,npmax,nblok,nxv,nxyp,ipbc)
            else
               call PGRPUSH2L(part,fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idimp&
     &,npmax,nblok,nxv,nypmx,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSRPUSH2(part,fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idimp&
     &,npmax,nblok,nxv,nxyp,ipbc)
            else
               call PGRPUSH2(part,fxy,npp,noff,qbm,dt,ci,ek,nx,ny,idimp,&
     &npmax,nblok,nxv,nypmx,ipbc)
            endif
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgrpush2
!
         subroutine ipgrjpost2(part,cu,npp,noff,qm,dt,ci,tdjpost,nx,ny,i&
     &pbc,inorder,djopt)
! deposit relativisitc current, 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, dt, ci, tdjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, nxyp, ipbc
         integer :: order, opt
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128
         integer, parameter :: n8 = 8, n12 = 12, n18 = 18, n27 = 27
         integer, dimension(n27,npd,size(part,3)) :: nn
         real, dimension(n27,npd,size(part,3)) :: amxy
         real :: tj
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(cu,2); nypmx = size(cu,3); nxyp = nxv*nypmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tj,dtime,-1)
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSRJPOST22L(part,cu,npp,noff,qm,dt,ci,nx,ny,idim&
     &p,npmax,nblok,nxv,nxyp,ipbc)
               else if (opt==VECTOR) then
                  call PGSRJOST22XL(part,cu,npp,noff,nn,amxy,qm,dt,ci,nx&
     &,ny,idimp,npmax,nblok,nxv,nxyp,npd,n8,ipbc)
               else
                  call PGRJPOST22L(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp&
     &,npmax,nblok,nxv,nypmx,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSRJPOST22(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp&
     &,npmax,nblok,nxv,nxyp,ipbc)
               else if (opt==VECTOR) then
                  call PGSRJOST22X(part,cu,npp,noff,nn,amxy,qm,dt,ci,nx,&
     &ny,idimp,npmax,nblok,nxv,nxyp,npd,n18,ipbc)
               else
                  call PGRJPOST22(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,&
     &npmax,nblok,nxv,nypmx,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSRJPOST2L(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp&
     &,npmax,nblok,nxv,nxyp,ipbc)
               else if (opt==VECTOR) then
                  call PGSRJOST2XL(part,cu,npp,noff,nn,amxy,qm,dt,ci,nx,&
     &ny,idimp,npmax,nblok,nxv,nxyp,npd,n12,ipbc)
               else
                  call PGRJPOST2L(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,&
     &npmax,nblok,nxv,nypmx,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSRJPOST2(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,&
     &npmax,nblok,nxv,nxyp,ipbc)
               else if (opt==VECTOR) then
                  call PGSRJOST2X(part,cu,npp,noff,nn,amxy,qm,dt,ci,nx,n&
     &y,idimp,npmax,nblok,nxv,nxyp,npd,n27,ipbc)
               else
                  call PGRJPOST2(part,cu,npp,noff,qm,dt,ci,nx,ny,idimp,n&
     &pmax,nblok,nxv,nypmx,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tj,dtime)
         tdjpost = tdjpost + tj
         end subroutine ipgrjpost2
!
         subroutine ipgrbpush2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,tp&
     &ush,nx,ny,ipbc,inorder,popt)
! push relativistic particles with 2d electromagnetic fields,
! 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, order, opt
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSRBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
            else
               call PGRBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSRBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
            else
               call PGRBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,nx,&
     &ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
            endif
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgrbpush2
!
         subroutine ipgrbpush23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,t&
     &push,nx,ny,ipbc,inorder,popt)
! push relativistic particles with 2-1/2d electromagnetic fields,
! 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, order, opt
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
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
                  call PGSRBPUSH22L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
               else
                  call PGRBPUSH22L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,e&
     &k,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSRBPUSH22(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,e&
     &k,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
               else
                  call PGRBPUSH22(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek&
     &,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSRBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,&
     &ek,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
               else
                  call PGRBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,e&
     &k,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSRBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,e&
     &k,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
               else
                  call PGRBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek&
     &,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgrbpush23
!
         subroutine iprretard2(part,npp,dtc,ci,nx,ny,ipbc,ndim)
! retards relativistic particle positions half time-step
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: ndim
         real :: dtc, ci
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, nd
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nd = 3
         if (present(ndim)) nd = ndim
         select case(nd)
         case (2)
            call PRRETARD22(part,npp,dtc,ci,nx,ny,idimp,npmax,nblok,ipbc&
     &)
         case (3)
            call PRRETARD2(part,npp,dtc,ci,nx,ny,idimp,npmax,nblok,ipbc)
         end select
         end subroutine iprretard2
!
         subroutine ipcptov2(part,npp,ci,ndim)
! convert momentum to velocity for relativistic particles, 1d partition
         implicit none
         real :: ci
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
         integer, optional :: ndim
! local data
         integer :: idimp, npmax, nblok, nd
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nd = 3
         if (present(ndim)) nd = ndim
         select case(nd)
         case (2)
            call PCPTOV22(part,npp,ci,idimp,npmax,nblok)
         case (3)
            call PCPTOV2(part,npp,ci,idimp,npmax,nblok)
         end select
         end subroutine ipcptov2
!
         subroutine iprpush2zf(part,npp,dt,ci,ek,tpush,nx,ny,ipbc)
! push relativistic particles with no forces, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         real :: dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
! initialize timer
         call wtimer(tp,dtime,-1)
         call PRPUSH2ZF(part,npp,dt,ci,ek,idimp,npmax,nblok,nx,ny,ipbc)
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine iprpush2zf
!
         subroutine iprpush23zf(part,npp,dt,ci,ek,tpush,nx,ny,ipbc,ndim)
! push 2d relativistic particles with no forces, 1d partition
         implicit none
         integer :: nx, ny, ipbc, ndim
         real :: dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
! initialize timer
         call wtimer(tp,dtime,-1)
         select case(ndim)
         case (2)
            call PRPUSH2ZF(part,npp,dt,ci,ek,idimp,npmax,nblok,nx,ny,ipb&
     &c)
         case (3)
            call PRPUSH23ZF(part,npp,dt,ci,ek,idimp,npmax,nblok,nx,ny,ip&
     &bc)
         end select
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine iprpush23zf
!
         subroutine ipgrcjpost2(part,fxy,npp,noff,cu,qm,qbm,dt,ci,tdcjpo&
     &st,inorder)
! deposit relativitic current density with 2d electrostatic fields,
! 1d partition
         implicit none
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, cu
         integer, dimension(:), pointer :: npp, noff
         integer, optional :: inorder
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, order
         real :: tdc
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tdc,dtime,-1)
         if (order==LINEAR) then
            call PGRCJPOST2L(part,fxy,npp,noff,cu,qm,qbm,dt,ci,idimp,npm&
     &ax,nblok,nxv,nypmx)
         else
            call PGRCJPOST2(part,fxy,npp,noff,cu,qm,qbm,dt,ci,idimp,npma&
     &x,nblok,nxv,nypmx)
         endif
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine ipgrcjpost2
!
      end module prbpush2d
 