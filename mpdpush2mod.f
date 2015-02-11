!-----------------------------------------------------------------------
!
      module mpdpush2d
!
! Fortran90 interface to parallel 2d PIC Fortran77 library mpdpush2lib.f
! mpdpush2mod.f contains multi-tasking interface procedures to process
!               particles with darwin electric and magnetic fields:
!               defines module mpdpush2d
! dmjpost => impgmjpost2 deposits momentum flux, with various
!            interpolations and optimizations.
!            calls MPGMJPOST2, MPGSMJPOST2, MPGMJPOST2L, MPGSMJPOST2L,
!            MPGMJPOST22, MPGSMJPOST22, MPGMJPOST22L, or MPGSMJPOST22L
! dcjpost => impgdcjpost2 deposits momentum flux, acceleration density,
!            and current density, with various interpolations and
!            optimizations.
!            calls MPGDCJPOST2, MPGSDCJPOST2, MPGDCJPOST2L,
!            MPGSDCJPOST2L, MPGDCJPOST22, MPGSDCJPOST22, MPGDCJPOST22L,
!            or MPGSDCJPOST22L
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 12, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use pdpush2d, only: wtimer
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, dmjpost, dcjpost
!
! buffer data for  momentum flux
      real, dimension(:,:,:,:,:), allocatable :: amup
      integer :: szbufa = 0
      save
!
! define interface to Fortran77 procedures
!
      interface
         subroutine MPGMJPOST2(part,amu,npp,nps,noff,qm,idimp,npmax,nblo&
     &k,nxv,nypmx,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(4,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSMJPOST2(part,amu,npp,nps,noff,qm,idimp,npmax,nbl&
     &ok,nxv,nxyp,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(4,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGMJPOST2L(part,amu,npp,nps,noff,qm,idimp,npmax,nbl&
     &ok,nxv,nypmx,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(4,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSMJPOST2L(part,amu,npp,nps,noff,qm,idimp,npmax,nb&
     &lok,nxv,nxyp,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(4,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(4,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGMJPOST22(part,amu,npp,nps,noff,qm,idimp,npmax,nbl&
     &ok,nxv,nypmx,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSMJPOST22(part,amu,npp,nps,noff,qm,idimp,npmax,nb&
     &lok,nxv,nxyp,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGMJPOST22L(part,amu,npp,nps,noff,qm,idimp,npmax,nb&
     &lok,nxv,nypmx,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSMJPOST22L(part,amu,npp,nps,noff,qm,idimp,npmax,n&
     &blok,nxv,nxyp,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: amu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm,&
     &qbm,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt
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
         subroutine MPGSDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm&
     &,qbm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, qbm, dt
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
         subroutine MPGDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm&
     &,qbm,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt
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
         subroutine MPGSDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,q&
     &m,qbm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, qbm, dt
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
         subroutine MPGDCJPOST22(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm,&
     &qbm,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSDCJPOST22(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm&
     &,qbm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxyp,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxyp,nblok,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGDCJPOST22L(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm&
     &,qbm,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, dcu, amu
         real, dimension(nxv,nypmx,nblok) :: bz
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSDCJPOST22L(part,fxy,bz,npp,nps,noff,cu,dcu,amu,q&
     &m,qbm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm, qbm, dt
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
      interface dmjpost
         module procedure impgmjpost2
      end interface
!
      interface dcjpost
         module procedure impgdcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impgmjpost2(part,amu,npp,noff,qm,tdcjpost,inorder,dj&
     &opt)
! multi-tasking momentum flux deposit, 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, tdcjpost
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
                  call MPGSMJPOST22L(part,amu,npp,nps,noff,qm,idimp,npma&
     &x,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
               else
                  call MPGMJPOST22L(part,amu,npp,nps,noff,qm,idimp,npmax&
     &,nblok,nxv,nypmx,amup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSMJPOST22(part,amu,npp,nps,noff,qm,idimp,npmax&
     &,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
               else
                  call MPGMJPOST22(part,amu,npp,nps,noff,qm,idimp,npmax,&
     &nblok,nxv,nypmx,amup,idtask,nmt,ierr)
               endif
            endif
         case (4)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSMJPOST2L(part,amu,npp,nps,noff,qm,idimp,npmax&
     &,nblok,nxv,nxyp,amup,idtask,nmt,ierr)
               else
                  call MPGMJPOST2L(part,amu,npp,nps,noff,qm,idimp,npmax,&
     &nblok,nxv,nypmx,amup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSMJPOST2(part,amu,npp,nps,noff,qm,idimp,npmax,&
     &nblok,nxv,nxyp,amup,idtask,nmt,ierr)
               else
                  call MPGMJPOST2(part,amu,npp,nps,noff,qm,idimp,npmax,n&
     &blok,nxv,nypmx,amup,idtask,nmt,ierr)
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
         end subroutine impgmjpost2
!
         subroutine impgdcjpost2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm&
     &,dt,tdcjpost,inorder,djopt)
! multi-tasking deposit momentum flux, acceleration density, and current
! density with 2-1/2d electromagnetic fields, 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, tdcjpost
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
                  call MPGSDCJPOST22L(part,fxy,bxy,npp,nps,noff,cu,dcu,a&
     &mu,qm,qbm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,i&
     &err)
               else
                  call MPGDCJPOST22L(part,fxy,bxy,npp,nps,noff,cu,dcu,am&
     &u,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,i&
     &err)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSDCJPOST22(part,fxy,bxy,npp,nps,noff,cu,dcu,am&
     &u,qm,qbm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ie&
     &rr)
               else
                  call MPGDCJPOST22(part,fxy,bxy,npp,nps,noff,cu,dcu,amu&
     &,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ie&
     &rr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MPGSDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,am&
     &u,qm,qbm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ie&
     &rr)
               else
                  call MPGDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,amu&
     &,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ie&
     &rr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MPGSDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu&
     &,qm,qbm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ier&
     &r)
               else
                  call MPGDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,&
     &qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ier&
     &r)
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
         end subroutine impgdcjpost2
!
      end module mpdpush2d
