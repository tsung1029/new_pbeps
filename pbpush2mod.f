!-----------------------------------------------------------------------
!
      module pbpush2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pbpush2lib.f
! pbpush2mod.f contains interface procedures to process particles with
!              magnetic fields:
!              defines module pbpush2d
! djpost => ipgjpost2 deposits current density, with various
!           interpolations and optimizations.
!           calls PGJPOST2, PGSJPOST2, PGSJOST2X, PGJPOST2L, PGSJPOST2L,
!           PGSJOST2XL, PGJPOST22, PGSJPOST22, PGSJOST22X, PGJPOST22L,
!           PGSJPOST22L, or PGSJOST22XL
! push => ipgbpush2 push particles with magnetic field and 2 component
!         electric field, with various interpolations and optimizations.
!         calls PGBPUSH2, PGSBPUSH2, PGBPUSH22L, or PGSBPUSH2L
! push3 => ipgbpush23 push particles with magnetic field and 3 component
!          electric field, with various interpolations and optimizations
!          calls PGBPUSH23, PGSBPUSH23, PGBPUSH23L, PGSBPUSH23L,
!          PGBPUSH22, PGSBPUSH22, GBPUSH22L, or PGSBPUSH22L
! retard => ipretard2 retard particle position a half time-step.
!           calls PRETARD2
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 5, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer
      public :: djpost, push, push3, retard
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PJDOST2(part,cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,npm&
     &ax,nblok,nxv,nypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: cux, cuy, cuz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGJPOST2(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nb&
     &lok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSJPOST2(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,n&
     &blok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PSJOST2X(part,cu,npp,noff,nn,amxy,qm,dt,nx,idimp,npm&
     &ax,nblok,nxv,nxvyp,npd,n27)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nxvyp, npd, n27
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n27,npd,nblok) :: nn
         real, dimension(n27,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGSJOST2X(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp&
     &,npmax,nblok,nxv,nxvyp,npd,n27,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n27
         integer :: ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n27,npd,nblok) :: nn
         real, dimension(n27,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PJDOST2L(part,cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,np&
     &max,nblok,nxv,nypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: cux, cuy, cuz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGJPOST2L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,n&
     &blok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSJPOST2L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,&
     &nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PSJOST2XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,idimp,np&
     &max,nblok,nxv,nxvyp,npd,n12)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nxvyp, npd, n12
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n12,npd,nblok) :: nn
         real, dimension(n12,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGSJOST2XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idim&
     &p,npmax,nblok,nxv,nxvyp,npd,n12,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n12
         integer :: ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n12,npd,nblok) :: nn
         real, dimension(n12,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGJPOST22(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,n&
     &blok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSJPOST22(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,&
     &nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSJOST22X(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idim&
     &p,npmax,nblok,nxv,nxvyp,npd,n18,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n18
         integer :: ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n18,npd,nblok) :: nn
         real, dimension(n18,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGJPOST22L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,&
     &nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSJPOST22L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax&
     &,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSJOST22XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idi&
     &mp,npmax,nblok,nxv,nxvyp,npd,n8,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxvyp, npd, n8
         integer :: ipbc
         real :: qm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2*nxvyp,nblok) :: cu
         integer, dimension(nblok) :: npp, noff
         integer, dimension(n8,npd,nblok) :: nn
         real, dimension(n8,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PBPUSH2(part,fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,id&
     &imp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy, bx, by, bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,i&
     &dimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PBPUSH2L(part,fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,i&
     &dimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy, bx, by, bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny&
     &,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
       interface
         subroutine PBPUSH2CQ(part,fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,&
     &idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy, bx, by, bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH2CQ(part,fxy,bxy,npp,noff,qbm,dt,ek,nx,ny,idi&
     &mp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PBPUSH2CL(part,fx,fy,bx,by,bz,npp,noff,qbm,dt,ek,nx,&
     &idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy, bx, by, bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH2CL(part,fxy,bxy,npp,noff,qbm,dt,ek,nx,ny,idi&
     &mp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny&
     &,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny&
     &,idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy
         real, dimension(3,nxv,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,n&
     &y,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy
         real, dimension(3,nxyp,nblok) :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH22(part,fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,i&
     &dimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSBPUSH22(part,fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBPUSH22L(part,fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSBPUSH22L(part,fxy,bz,npp,noff,qbm,dt,dtc,ek,nx,ny&
     &,idimp,npmax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PRETARD2(part,npp,dtc,nx,ny,idimp,npmax,nblok,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, ipbc
         real :: dtc
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         end subroutine
      end interface    
!
! define generic interface to Fortran90 library
!
      interface djpost
         module procedure ipgjpost2
      end interface
!
      interface push
         module procedure ipgbpush2
!        module procedure ipgbpush2c
      end interface
!
      interface push3
         module procedure ipgbpush23
      end interface
!
      interface retard
         module procedure ipretard2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!
         subroutine ipgjpost2(part,cu,npp,noff,qm,dt,tdjpost,nx,ny,ipbc,&
     &inorder,djopt)
! deposit current, 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, dt, tdjpost
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
                  call PGSJPOST22L(part,cu,npp,noff,qm,dt,nx,ny,idimp,np&
     &max,nblok,nxv,nxyp,ipbc)
               else if (opt==VECTOR) then
                  call PGSJOST22XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,&
     &idimp,npmax,nblok,nxv,nxyp,npd,n8,ipbc)
               else
                  call PGJPOST22L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSJPOST22(part,cu,npp,noff,qm,dt,nx,ny,idimp,npm&
     &ax,nblok,nxv,nxyp,ipbc)

               else if (opt==VECTOR) then
                  call PGSJOST22X(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,i&
     &dimp,npmax,nblok,nxv,nxyp,npd,n18,ipbc)
               else
                  call PGJPOST22(part,cu,npp,noff,qm,dt,nx,ny,idimp,npma&
     &x,nblok,nxv,nypmx,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSJPOST2L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npm&
     &ax,nblok,nxv,nxyp,ipbc)
               else if (opt==VECTOR) then
                  call PGSJOST2XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,i&
     &dimp,npmax,nblok,nxv,nxyp,npd,n12,ipbc)
               else
                  call PGJPOST2L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npma&
     &x,nblok,nxv,nypmx,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSJPOST2(part,cu,npp,noff,qm,dt,nx,ny,idimp,npma&
     &x,nblok,nxv,nxyp,ipbc)
               else if (opt==VECTOR) then
                  call PGSJOST2X(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,id&
     &imp,npmax,nblok,nxv,nxyp,npd,n27,ipbc)
               else
                  call PGJPOST2(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax&
     &,nblok,nxv,nypmx,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tj,dtime)
         tdjpost = tdjpost + tj
         end subroutine ipgjpost2
!
         subroutine ipgbpush2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,tpush,&
     &nx,ny,ipbc,inorder,popt)
! push particles with 2d electromagnetic fields, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
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
               call PGSBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny&
     &,idimp,npmax,nblok,nxv,nxyp,ipbc)
            else
               call PGBPUSH2L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nypmx,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,&
     &idimp,npmax,nblok,nxv,nxyp,ipbc)
            else
               call PGBPUSH2(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,ny,i&
     &dimp,npmax,nblok,nxv,nypmx,ipbc)
            endif
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgbpush2
!
         subroutine ipgbpush2cq(part,fxy,bxy,npp,noff,qbm,dt,ek,tpush,nx&
     &,ny,ipbc,inorder)
! push particles with 2d electromagnetic fields, 1d partition
! with correction to Boris Mover
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm, dt, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, order
         real :: tp
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tp,dtime,-1)
         if (order==LINEAR) then
            call PGBPUSH2CL(part,fxy,bxy,npp,noff,qbm,dt,ek,nx,ny,idimp,&
     &npmax,nblok,nxv,nypmx,ipbc)
         else
            call PGBPUSH2CQ(part,fxy,bxy,npp,noff,qbm,dt,ek,nx,ny,idimp,&
     &npmax,nblok,nxv,nypmx,ipbc)
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgbpush2cq
!
         subroutine ipgbpush23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,tpush&
     &,nx,ny,ipbc,inorder,popt)
! push particles with 2-1/2d electromagnetic fields, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
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
                  call PGSBPUSH22L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
               else
                  call PGBPUSH22L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSBPUSH22(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
               else
                  call PGBPUSH22(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,&
     &ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,n&
     &x,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
               else
                  call PGBPUSH23L(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx&
     &,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                 call PGSBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,&
     &ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
               else
                  call PGBPUSH23(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,nx,&
     &ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgbpush23
!
         subroutine ipretard2(part,npp,dtc,nx,ny,ipbc,ndim)
! retards particle positions half time-step
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: ndim
         real :: dtc
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PRETARD2(part,npp,dtc,nx,ny,idimp,npmax,nblok,ipbc)
         end subroutine ipretard2
!
      end module pbpush2d
 