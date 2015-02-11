!-----------------------------------------------------------------------
!
      module prdpush2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library prdpush2lib.f
! prdpush2mod.f contains interface procedures to process relativistic
!               particles with darwin electric and magnetic fields:
!               defines module prdpush2d
! rdmjpost => ipgrmjpost2 deposits relativistic momentum flux, with
!             various interpolations and optimizations.
!             calls PGRMJPOST2, PGSRMJPOST2, PGRMJPOST2L, PGSRMJPOST2L,
!             PGRMJPOST22, PGSRMJPOST22, PGRMJPOST22L, or PGSRMJPOST22L
! rdcjpost => ipgrdcjpost2 deposits relativistic momentum flux,
!             acceleration density, and current density, with various
!             interpolations and optimizations.
!             calls PGRDCJPOST2, PGSRDCJPOST2, PGRDCJPOST2L,
!             PGSRDCJPOST2L, PGRDCJPOST22, PGSRDCJPOST22, PGRDCJPOST22L,
!             or PGSRDCJPOST22L
! written by viktor k. decyk, ucla
! copyright 2006, regents of the university of california
! update: october 22, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use p0d, only: wtimer
      use prbpush2d, only: rdjpost, rpush3, ipcptov2
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, rdmjpost, rdcjpost
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PGRMJPOST2(part,amu,npp,noff,qm,ci,idimp,npmax,nblok&
     &,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRMJPOST2(part,amu,npp,noff,qm,ci,idimp,npmax,nblo&
     &k,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,&
     &dt,ci,idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm&
     &,dt,ci,idimp,npmax,nblok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRMJPOST22(part,amu,npp,noff,qm,ci,idimp,npmax,nblo&
     &k,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRMJPOST22(part,amu,npp,noff,qm,ci,idimp,npmax,nbl&
     &ok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRDCJPOST22(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,&
     &dt,ci,idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRDCJPOST22(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm&
     &,dt,ci,idimp,npmax,nblok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRMJPOST2L(part,amu,npp,noff,qm,ci,idimp,npmax,nblo&
     &k,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRMJPOST2L(part,amu,npp,noff,qm,ci,idimp,npmax,nbl&
     &ok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm&
     &,dt,ci,idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qb&
     &m,dt,ci,idimp,npmax,nblok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRMJPOST22L(part,amu,npp,noff,qm,ci,idimp,npmax,nbl&
     &ok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRMJPOST22L(part,amu,npp,noff,qm,ci,idimp,npmax,nb&
     &lok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRDCJPOST22L(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm&
     &,dt,ci,idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSRDCJPOST22L(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qb&
     &m,dt,ci,idimp,npmax,nblok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdmjpost
         module procedure ipgrmjpost2
      end interface
!
      interface rdcjpost
         module procedure ipgrdcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipgrmjpost2(part,amu,npp,noff,qm,ci,tdcjpost,inorder&
     &,djopt)
! deposit momentum flux with 2-1/2d electromagnetic fields,
! with relativistic particles, and 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, ci, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: amu
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp
         integer :: order, opt
         real :: tdc
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(amu,2); nypmx = size(amu,3); nxyp = nxv*nypmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tdc,dtime,-1)
         select case(size(amu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSRMJPOST22L(part,amu,npp,noff,qm,ci,idimp,npmax&
     &,nblok,nxv,nxyp)
               else
                  call PGRMJPOST22L(part,amu,npp,noff,qm,ci,idimp,npmax,&
     &nblok,nxv,nypmx)

               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSRMJPOST22(part,amu,npp,noff,qm,ci,idimp,npmax,&
     &nblok,nxv,nxyp)
               else
                  call PGRMJPOST22(part,amu,npp,noff,qm,ci,idimp,npmax,n&
     &blok,nxv,nypmx)
               endif
            endif
         case (4)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSRMJPOST2L(part,amu,npp,noff,qm,ci,idimp,npmax,&
     &nblok,nxv,nxyp)
               else
                  call PGRMJPOST2L(part,amu,npp,noff,qm,ci,idimp,npmax,n&
     &blok,nxv,nypmx)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSRMJPOST2(part,amu,npp,noff,qm,ci,idimp,npmax,n&
     &blok,nxv,nxyp)
               else
                  call PGRMJPOST2(part,amu,npp,noff,qm,ci,idimp,npmax,nb&
     &lok,nxv,nypmx)
               endif
            endif
         end select
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine ipgrmjpost2
!
         subroutine ipgrdcjpost2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm&
     &,dt,ci,tdcjpost,inorder,djopt)
! deposit momentum flux, acceleration density, and current density
! with 2-1/2d electromagnetic fields, with relativistic particles, and
! 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy, cu, dcu, amu
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, order, opt
         real :: tdc
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tdc,dtime,-1)
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSRDCJPOST22L(part,fxy,bxy,npp,noff,cu,dcu,amu,q&
     &m,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp)
               else
                  call PGRDCJPOST22L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm&
     &,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSRDCJPOST22(part,fxy,bxy,npp,noff,cu,dcu,amu,qm&
     &,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp)
               else
                  call PGRDCJPOST22(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,&
     &qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSRDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm&
     &,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp)
               else
                  call PGRDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,&
     &qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSRDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,&
     &qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp)
               else
                  call PGRDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,q&
     &bm,dt,ci,idimp,npmax,nblok,nxv,nypmx)
               endif
            endif
         end select
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine ipgrdcjpost2
!
      end module prdpush2d
