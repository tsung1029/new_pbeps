!-----------------------------------------------------------------------
!
      module pdpush2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pdpush2lib.f
! pdpush2mod.f contains interface procedures to process particles with
!              darwin electric and magnetic fields:
!              defines module pdpush2d
! dmjpost => ipgmjpost2 deposits momentum flux, with various
!            interpolations and optimizations.
!            calls PGMJPOST2, PGSMJPOST2, PGMJPOST2L, PGSMJPOST2L,
!            PGMJPOST22, PGSMJPOST22, PGMJPOST22L, or PGSMJPOST22L
! dcjpost => ipgdcjpost2 deposits momentum flux, acceleration density,
!            and current density, with various interpolations and
!            optimizations.
!            calls PGDCJPOST2, PGSDCJPOST2, PGDCJPOST2L, PGSDCJPOST2L,
!            PGDCJPOST22, PGSDCJPOST22, PGDCJPOST22L, or PGSDCJPOST22L
! written by viktor k. decyk, ucla
! copyright 2006, regents of the university of california
! update: december 12, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, dmjpost, dcjpost
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PGMJPOST2(part,amu,npp,noff,qm,idimp,npmax,nblok,nxv&
     &,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSMJPOST2(part,amu,npp,noff,qm,idimp,npmax,nblok,nx&
     &v,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,d&
     &t,idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,&
     &dt,idimp,npmax,nblok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGMJPOST22(part,amu,npp,noff,qm,idimp,npmax,nblok,nx&
     &v,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSMJPOST22(part,amu,npp,noff,qm,idimp,npmax,nblok,n&
     &xv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGDCJPOST22(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,d&
     &t,idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSDCJPOST22(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,&
     &dt,idimp,npmax,nblok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGMJPOST2L(part,amu,npp,noff,qm,idimp,npmax,nblok,nx&
     &v,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSMJPOST2L(part,amu,npp,noff,qm,idimp,npmax,nblok,n&
     &xv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,&
     &dt,idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm&
     &,dt,idimp,npmax,nblok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGMJPOST22L(part,amu,npp,noff,qm,idimp,npmax,nblok,n&
     &xv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSMJPOST22L(part,amu,npp,noff,qm,idimp,npmax,nblok,&
     &nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGDCJPOST22L(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,&
     &dt,idimp,npmax,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSDCJPOST22L(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm&
     &,dt,idimp,npmax,nblok,nxv,nxyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dmjpost
         module procedure ipgmjpost2
      end interface
!
      interface dcjpost
         module procedure ipgdcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipgmjpost2(part,amu,npp,noff,qm,tdcjpost,inorder,djo&
     &pt)
! deposit momentum flux with 2-1/2d electromagnetic fields, 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, tdcjpost
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
                  call PGSMJPOST22L(part,amu,npp,noff,qm,idimp,npmax,nbl&
     &ok,nxv,nxyp)
               else
                  call PGMJPOST22L(part,amu,npp,noff,qm,idimp,npmax,nblo&
     &k,nxv,nypmx)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSMJPOST22(part,amu,npp,noff,qm,idimp,npmax,nblo&
     &k,nxv,nxyp)
               else
                  call PGMJPOST22(part,amu,npp,noff,qm,idimp,npmax,nblok&
     &,nxv,nypmx)
               endif
            endif
         case (4)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSMJPOST2L(part,amu,npp,noff,qm,idimp,npmax,nblo&
     &k,nxv,nxyp)
               else
                  call PGMJPOST2L(part,amu,npp,noff,qm,idimp,npmax,nblok&
     &,nxv,nypmx)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSMJPOST2(part,amu,npp,noff,qm,idimp,npmax,nblok&
     &,nxv,nxyp)
               else
                  call PGMJPOST2(part,amu,npp,noff,qm,idimp,npmax,nblok,&
     &nxv,nypmx)
               endif
            endif
         end select
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine ipgmjpost2
!
         subroutine ipgdcjpost2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,&
     &dt,tdcjpost,inorder,djopt)
! deposit momentum flux, acceleration density, and current density
! with 2-1/2d electromagnetic fields, 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, tdcjpost
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
                  call PGSDCJPOST22L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm&
     &,qbm,dt,idimp,npmax,nblok,nxv,nxyp)
               else
                  call PGDCJPOST22L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,&
     &qbm,dt,idimp,npmax,nblok,nxv,nypmx)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSDCJPOST22(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,&
     &qbm,dt,idimp,npmax,nblok,nxv,nxyp)
               else
                  call PGDCJPOST22(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,q&
     &bm,dt,idimp,npmax,nblok,nxv,nypmx)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call PGSDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,&
     &qbm,dt,idimp,npmax,nblok,nxv,nxyp)
               else
                  call PGDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,q&
     &bm,dt,idimp,npmax,nblok,nxv,nypmx)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call PGSDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,q&
     &bm,dt,idimp,npmax,nblok,nxv,nxyp)
               else
                   call PGDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,q&
     &bm,dt,idimp,npmax,nblok,nxv,nypmx)
               endif
            endif
         end select
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine ipgdcjpost2
!
      end module pdpush2d
