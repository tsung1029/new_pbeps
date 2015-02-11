!-----------------------------------------------------------------------
!
      module mprdpush2d
!
! Fortran90 interface to parallel 2d PIC Fortran77 library
! mprdpush2lib.f
! mprdpush2mod.f contains multi-tasking interface procedures to process
!                relativistic particles with darwin electric and
!                magnetic fields:
!                defines module mprdpush2d
! rdmjpost => impgrmjpost2 deposits relativistic momentum flux, with
!             various interpolations and optimizations.
!             calls MPGRMJPOST2, MPGSRMJPOST2, MPGRMJPOST2L,
!             MPGSRMJPOST2L, MPGRMJPOST22, MPGSRMJPOST22, MPGRMJPOST22L,
!             or MPGSRMJPOST22L
! rdcjpost => impgrdcjpost2 deposits relativistic momentum flux,
!             acceleration density, and current density, with various
!             interpolations and optimizations.
!             calls MPGRDCJPOST2, MPGSRDCJPOST2, MPGRDCJPOST2L,
!             MPGSRDCJPOST2L, MPGRDCJPOST22, MPGSRDCJPOST22,
!             MPGRDCJPOST22L, or MPGSRDCJPOST22L
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 12, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use prdpush2d, only: wtimer
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, rdmjpost, rdcjpost
!
! buffer data for  momentum flux
      real, dimension(:,:,:,:,:), allocatable :: amup
      integer :: szbufa = 0
      save
!
! define interface to Fortran77 procedures
!
      interface
         subroutine MPGRMJPOST2(part,amu,npp,nps,noff,qm,ci,idimp,npmax,&
     &nblok,nxv,nypmx,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(4,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRMJPOST2(part,amu,npp,nps,noff,qm,ci,idimp,npmax&
     &,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(4,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRMJPOST2L(part,amu,npp,nps,noff,qm,ci,idimp,npmax&
     &,nblok,nxv,nypmx,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(4,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRMJPOST2L(part,amu,npp,nps,noff,qm,ci,idimp,npma&
     &x,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(4,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRMJPOST22(part,amu,npp,nps,noff,qm,ci,idimp,npmax&
     &,nblok,nxv,nypmx,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRMJPOST22(part,amu,npp,nps,noff,qm,ci,idimp,npma&
     &x,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRMJPOST22L(part,amu,npp,nps,noff,qm,ci,idimp,npma&
     &x,nblok,nxv,nypmx,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRMJPOST22L(part,amu,npp,nps,noff,qm,ci,idimp,npm&
     &ax,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm&
     &,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ie&
     &rr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxv,nypmx,nblok,nmt) :: cup, dcup
         real, dimension(4,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,q&
     &m,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ie&
     &rr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxyp,nblok,nmt) :: cup, dcup
         real, dimension(4,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,q&
     &m,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,i&
     &err)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxv,nypmx,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxv,nypmx,nblok,nmt) :: cup, dcup
         real, dimension(4,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,&
     &qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,i&
     &err)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(3,nxyp,nblok) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(3,nxyp,nblok,nmt) :: cup, dcup
         real, dimension(4,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRDCJPOST22(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm&
     &,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ie&
     &rr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRDCJPOST22(part,fxy,bz,npp,nps,noff,cu,dcu,amu,q&
     &m,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ie&
     &rr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGRDCJPOST22L(part,fxy,bz,npp,nps,noff,cu,dcu,amu,q&
     &m,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,i&
     &err)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSRDCJPOST22L(part,fxy,bz,npp,nps,noff,cu,dcu,amu,&
     &qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,i&
     &err)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdmjpost
         module procedure impgrmjpost2
      end interface
!
      interface rdcjpost
         module procedure impgrdcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impgrmjpost2(part,amu,npp,noff,qm,ci,tdcjpost,inorde&
     &r,djopt)
! multi-tasking momentum flux deposit, 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, ci, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: amu
         integer, dimension(:), pointer :: npp, noff
! local data
         integer:: idimp, npmax, nblok, nxv, nypmx, nxyp, nnxyp
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
!        real, dimension(size(amu,1),size(amu,2),size(amu,3),size(amu,4)&
!    &,ntasks) :: amup
         real :: tdc
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(amu,2); nypmx = size(amu,3); nxyp = nxv*nypmx
         nnxyp = size(amu,1)*nxyp
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! check if size of momentum flux buffer has changed
         if (szbufa < nnxyp*nblok) then
            if (szbufa /= 0) deallocate(amup)
! allocate buffer
            allocate(amup(size(amu,1),nxv,nypmx,nblok,ntasks))
            szbufa = nnxyp*nblok
         endif
! initialize timer
         call wtimer(tdc,dtime,-1)
         select case(size(amu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSRMJPOST22L(part,amu,npp,nps,noff,qm,ci,idimp,&
     &npmax,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
               else
                  call MPGRMJPOST22L(part,amu,npp,nps,noff,qm,ci,idimp,n&
     &pmax,nblok,nxv,nypmx,amup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSRMJPOST22(part,amu,npp,nps,noff,qm,ci,idimp,n&
     &pmax,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
               else
                  call MPGRMJPOST22(part,amu,npp,nps,noff,qm,ci,idimp,np&
     &max,nblok,nxv,nypmx,amup,idtask,nmt,ierr)
               endif
            endif
         case (4)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSRMJPOST2L(part,amu,npp,nps,noff,qm,ci,idimp,n&
     &pmax,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
               else
                  call MPGRMJPOST2L(part,amu,npp,nps,noff,qm,ci,idimp,np&
     &max,nblok,nxv,nypmx,amup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSRMJPOST2(part,amu,npp,nps,noff,qm,ci,idimp,np&
     &max,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
               else
                  call MPGRMJPOST2(part,amu,npp,nps,noff,qm,ci,idimp,npm&
     &ax,nblok,nxv,nypmx,amup,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine impgrmjpost2
!
         subroutine impgrdcjpost2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qb&
     &m,dt,ci,tdcjpost,inorder,djopt)
! multi-tasking deposit momentum flux, acceleration density, and current
! density with 2-1/2d electromagnetic fields, 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy, cu, dcu, amu
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, nnxyp
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(ntasks) :: idtask
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),size(cu,4),nta&
!    &sks) :: cup
!        real, dimension(size(dcu,1),size(dcu,2),size(dcu,3),size(dcu,4)&
!    &,ntasks) :: dcup
!        real, dimension(size(amu,1),size(amu,2),size(amu,3),size(amu,4)&
!    &,ntasks) :: amup
         real, dimension(:,:,:,:,:), allocatable, save :: cup, dcup
         integer, save :: szbufc = 0, szbufd = 0
         real :: tdc
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! check if size of current buffer has changed
         nnxyp = size(cu,1)*nxyp
         if (szbufc < nnxyp*nblok) then
            if (szbufc /= 0) deallocate(cup)
            allocate(cup(size(cu,1),nxv,nypmx,nblok,ntasks))
            szbufc = nnxyp*nblok
         endif
! check if size of acceleration buffer has changed
         nnxyp = size(dcu,1)*nxyp
         if (szbufd < nnxyp*nblok) then
            if (szbufd /= 0) deallocate(dcup)
            allocate(dcup(size(dcu,1),nxv,nypmx,nblok,ntasks))
            szbufd = nnxyp*nblok
         endif
! check if size of momentum flux buffer has changed
         nnxyp = size(amu,1)*nxyp
         if (szbufa < nnxyp*nblok) then
            if (szbufa /= 0) deallocate(amup)
            allocate(amup(size(amu,1),nxv,nypmx,nblok,ntasks))
            szbufa = nnxyp*nblok
         endif
! initialize timer
         call wtimer(tdc,dtime,-1)
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSRDCJPOST22L(part,fxy,bxy,npp,nps,noff,cu,dcu,&
     &amu,qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,n&
     &mt,ierr)
               else
                  call MPGRDCJPOST22L(part,fxy,bxy,npp,nps,noff,cu,dcu,a&
     &mu,qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,n&
     &mt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSRDCJPOST22(part,fxy,bxy,npp,nps,noff,cu,dcu,a&
     &mu,qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nm&
     &t,ierr)
               else
                  call MPGRDCJPOST22(part,fxy,bxy,npp,nps,noff,cu,dcu,am&
     &u,qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nm&
     &t,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSRDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,a&
     &mu,qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nm&
     &t,ierr)
               else
                  call MPGRDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,am&
     &u,qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nm&
     &t,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSRDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,am&
     &u,qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt&
     &,ierr)
               else
                  call MPGRDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu&
     &,qm,qbm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt&
     &,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine impgrdcjpost2
!
      end module mprdpush2d
