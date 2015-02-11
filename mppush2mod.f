!-----------------------------------------------------------------------
!
      module mppush2d
!
! Fortran90 interface to parallel 2d PIC Fortran77 library mppush2lib.f
! mppush2mod.f contains multi-tasking interface procedures to process
!              particles:
!              defines module mppush2d
! dpost => impgpost2 deposits charge density, with various
!          interpolations and optimizations.
!          calls MPGPOST2, MPGSPOST2, MPGSOST2X, MPGPOST2L, MPGSPOST2L,
!          or MPGSOST2XL
! push => impgpush2 push particles, with various interpolations and
!         optimizations.
!         calls MPGPUSH2, MPGSPUSH2, MPGPUSH2L, or MPGSPUSH2L
! pushzf => imppush2zf, push particles with no forces.
!           calls MPPUSH2ZF
! sortp => impsortp2y sorts particles by y grid using memory-conserving
!          bin sort, with various interpolations.
!          calls MPSORTP2Y, or MPSORTP2YL
! sortp => imdpsortp2y sorts particles by y grid using optimized bin
!          sort with various interpolations.
!          calls MPDSORTP2Y, or MPDSORTP2YL
! gcjpost => impgcjpost2 deposits time-centered current density with
!            2d electrostatic fields.
!            calls MPGCJPOST2, or MPGCJPOST2L
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 12, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use ppush2d, only: wtimer, countp, prmove, initmomt2, premoment2, &
     &primoment2
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: dpost, push, sortp, countp, prmove, pushzf, gcjpost
      public :: initmomt2, premoment2, primoment2
!
! define interface to Fortran77 procedures
      interface
         subroutine MPDOST2(part,q,npp,nps,noff,qm,nx,idimp,npmax,nblok,&
     &nxv,nypmx,qp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nxv,nypmx,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGPOST2(part,q,npp,nps,noff,qm,idimp,npmax,nblok,nx&
     &v,nypmx,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nxv,nypmx,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSPOST2(part,q,npp,nps,noff,qm,idimp,npmax,nblok,n&
     &xv,nxyp,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxyp,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nxyp,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSOST2X(part,q,npp,nps,noff,nn,amxy,qm,nx,idimp,npm&
     &ax,nblok,nxv,nxvyp,npd,nine,qp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nxvyp, npd, nine, nmt
         integer :: ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxvyp,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(nine,npd,nblok,nmt+1) :: nn
         real, dimension(nine,npd,nblok,nmt+1) :: amxy
         real, dimension(nxvyp,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSOST2X(part,q,npp,nps,noff,nn,amxy,qm,idimp,npmax&
     &,nblok,nxv,nxvyp,npd,nine,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxvyp, npd, nine, nmt
         integer :: ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxvyp,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(nine,npd,nblok,nmt+1) :: nn
         real, dimension(nine,npd,nblok,nmt+1) :: amxy
         real, dimension(nxvyp,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPDOST2L(part,q,npp,nps,noff,qm,nx,idimp,npmax,nblok&
     &,nxv,nypmx,qp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nxv,nypmx,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGPOST2L(part,q,npp,nps,noff,qm,idimp,npmax,nblok,n&
     &xv,nypmx,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nxv,nypmx,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSPOST2L(part,q,npp,nps,noff,qm,idimp,npmax,nblok,&
     &nxv,nxyp,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp, nmt, ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxyp,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nxyp,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSOST2XL(part,q,npp,nps,noff,nn,amxy,qm,nx,idimp,np&
     &max,nblok,nxv,nxvyp,npd,ifour,qp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nxvyp, npd, ifour, nmt
         integer :: ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxvyp,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(ifour,npd,nblok,nmt+1) :: nn
         real, dimension(ifour,npd,nblok,nmt+1) :: amxy
         real, dimension(nxvyp,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSOST2XL(part,q,npp,nps,noff,nn,amxy,qm,idimp,npma&
     &x,nblok,nxv,nxvyp,npd,ifour,qp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxvyp, npd, ifour, nmt
         integer :: ierr
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxvyp,nblok) :: q
         integer, dimension(nblok) :: npp, nps, noff
         integer, dimension(ifour,npd,nblok,nmt+1) :: nn
         real, dimension(ifour,npd,nblok,nmt+1) :: amxy
         real, dimension(nxvyp,nblok,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPPUSH2(part,fx,fy,npp,nps,noff,qbm,dt,ek,nx,idimp,n&
     &pmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGPUSH2(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idimp&
     &,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSPUSH2(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idim&
     &p,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPPUSH2L(part,fx,fy,npp,nps,noff,qbm,dt,ek,nx,idimp,&
     &npmax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idim&
     &p,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGSPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idi&
     &mp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc, nmt
         integer :: ierr
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSORTP2Y(part,pt,ip,npic,npp,nps,noff,nyp,idimp,npm&
     &ax,nblok,nypm1,npicp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nypm1, nmt, ierr
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(npmax,nblok) :: pt
         integer, dimension(npmax,nblok) :: ip
         integer, dimension(nypm1,nblok) :: npic
         integer, dimension(nblok) :: npp, nps, noff, nyp
         integer, dimension(nypm1,nblok,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPSORTP2YL(part,pt,ip,npic,npp,nps,noff,nyp,idimp,np&
     &max,nblok,nypm1,npicp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nypm1, nmt, ierr
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(npmax,nblok) :: pt
         integer, dimension(npmax,nblok) :: ip
         integer, dimension(nypm1,nblok) :: npic
         integer, dimension(nblok) :: npp, nps, noff, nyp
         integer, dimension(nypm1,nblok,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPDSORTP2Y(parta,partb,npic,npp,nps,noff,nyp,idimp,n&
     &pmax,nblok,nypm1,npicp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nypm1, nmt, ierr
         real, dimension(idimp,npmax,nblok) :: parta, partb
         integer, dimension(nypm1,nblok) :: npic
         integer, dimension(nblok) :: npp, nps, noff, nyp
         integer, dimension(nypm1,nblok,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPDSORTP2YL(parta,partb,npic,npp,nps,noff,nyp,idimp,&
     &npmax,nblok,nypm1,npicp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nypm1, nmt, ierr
         real, dimension(idimp,npmax,nblok) :: parta, partb
         integer, dimension(nypm1,nblok) :: npic
         integer, dimension(nblok) :: npp, nps, noff, nyp
         integer, dimension(nypm1,nblok,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPPUSH2ZF(part,npp,nps,dt,ek,idimp,npmax,nblok,nx,ny&
     &,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nx, ny, ipbc, nmt, ierr
         real :: dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGCJPOST2(part,fxy,npp,nps,noff,cu,qm,qbm,dt,idimp,&
     &npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu
         integer, dimension(nblok) :: npp, nps, noff
         real, dimension(2,nxv,nypmx,nblok,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPGCJPOST2L(part,fxy,npp,nps,noff,cu,qm,qbm,dt,idimp&
     &,npmax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx, nmt, ierr
         real :: qm, qbm, dt
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
      interface dpost
         module procedure impgpost2
      end interface
!
      interface push
         module procedure impgpush2
      end interface
!
      interface pushzf
         module procedure imppush2zf
      end interface
!
      interface sortp
         module procedure impsortp2y
         module procedure imdpsortp2y
      end interface
!
      interface gcjpost
         module procedure impgcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impgpost2(part,q,qm,npp,noff,tdpost,inorder,dopt)
! multi-tasking charge deposit, 1d partition
         implicit none
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: npp, noff
! local data
         integer:: idimp, npmax, nblok, nxv, nypmx, nxyp
         integer :: nmt, order, opt, ierr
         integer, dimension(size(npp)) :: nps
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128, ifour = 4, nine = 9
         integer, dimension(ntasks) :: idtask
!        real, dimension(size(q,1),size(q,2),size(q,3),ntasks) :: qp
         real, dimension(:,:,:,:), allocatable, save :: qp
         integer, save :: szbuf = 0
         integer, dimension(nine,npd,size(part,3),ntasks+1) :: nn
         real, dimension(nine,npd,size(part,3),ntasks+1) :: amxy
         real :: td
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(q,1); nypmx = size(q,2); nxyp = nxv*nypmx
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(dopt)) opt = dopt
! check if size of buffer has changed
         if (szbuf < nxyp*nblok) then
            if (szbuf /= 0) deallocate(qp)
! allocate buffer
            allocate(qp(nxv,nypmx,nblok,ntasks))
            szbuf = nxyp*nblok
         endif
! initialize timer
         call wtimer(td,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MPGSPOST2L(part,q,npp,nps,noff,qm,idimp,npmax,nblok,&
     &nxv,nxyp,qp,idtask,nmt,ierr)
            else if (opt==VECTOR) then
               call MPGSOST2XL(part,q,npp,nps,noff,nn,amxy,qm,idimp,npma&
     &x,nblok,nxv,nxyp,npd,ifour,qp,idtask,nmt,ierr)
            else
               call MPGPOST2L(part,q,npp,nps,noff,qm,idimp,npmax,nblok,n&
     &xv,nypmx,qp,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSPOST2(part,q,npp,nps,noff,qm,idimp,npmax,nblok,n&
     &xv,nxyp,qp,idtask,nmt,ierr)
            else if (opt==VECTOR) then
               call MPGSOST2X(part,q,npp,nps,noff,nn,amxy,qm,idimp,npmax&
     &,nblok,nxv,nxyp,npd,nine,qp,idtask,nmt,ierr)
            else
               call MPGPOST2(part,q,npp,nps,noff,qm,idimp,npmax,nblok,nx&
     &v,nypmx,qp,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(td,dtime)
         tdpost = tdpost + td
         end subroutine impgpost2
!
         subroutine impgpush2(part,fxy,npp,noff,qbm,dt,ek,tpush,nx,ny,ip&
     &bc,inorder,popt)
! multi-tasking particle push with 2d electrostatic fields, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy
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
               call MPGSPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idi&
     &mp,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
            else
               call MPGPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idim&
     &p,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
            endif
         else
            if (opt==LOOKAHEAD) then
               call MPGSPUSH2(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idim&
     &p,npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
            else
               call MPGPUSH2(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idimp&
     &,npmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
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
         end subroutine impgpush2
!
         subroutine impsortp2y(part,pt,ip,npp,noff,nyp,npic,tsort,inorde&
     &r)
! multi-tasking particle sort by y grid using memory-conserving bin sort
         implicit none
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: pt
         integer, dimension(:,:), pointer :: ip, npic
         integer, dimension(:), pointer :: npp, noff, nyp
! local data
         integer :: idimp, npmax, nblok, nypm1, nmt, order, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(size(npic,1),size(npic,2),ntasks) :: npicp
         integer, dimension(ntasks) :: idtask
         real :: ts
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); nypm1 = size(npic,1)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,dtime,-1)
         if (order==LINEAR) then
            call MPSORTP2YL(part,pt,ip,npic,npp,nps,noff,nyp,idimp,npmax&
     &,nblok,nypm1,npicp,idtask,nmt,ierr)
         else
            call MPSORTP2Y(part,pt,ip,npic,npp,nps,noff,nyp,idimp,npmax,&
     &nblok,nypm1,npicp,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(ts,dtime)
         tsort = tsort + ts
         end subroutine impsortp2y
!
         subroutine imdpsortp2y(parta,partb,npp,noff,nyp,npic,tsort,inor&
     &der)
! multi-tasking particle sort by y grid using optimized bin sort
         implicit none
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:,:), pointer :: parta, partb
         integer, dimension(:), pointer :: npp, noff, nyp
         integer, dimension(:,:), pointer :: npic
! local data
         integer :: idimp, npmax, nblok, nypm1, nmt, order, ierr
         integer, dimension(size(npp)) :: nps
         integer, dimension(size(npic,1),size(npic,2),ntasks) :: npicp
         integer, dimension(ntasks) :: idtask
         real, dimension(:,:,:), pointer :: part
         real :: ts
         double precision :: dtime
         idimp = size(parta,1); npmax = size(parta,2)
         nblok = size(parta,3); nypm1 = size(npic,1)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,dtime,-1)
         if (order==LINEAR) then
            call MPDSORTP2YL(parta,partb,npic,npp,nps,noff,nyp,idimp,npm&
     &ax,nblok,nypm1,npicp,idtask,nmt,ierr)
         else
            call MPDSORTP2Y(parta,partb,npic,npp,nps,noff,nyp,idimp,npma&
     &x,nblok,nypm1,npicp,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         part => parta
         parta => partb
         partb => part
! record time
         call wtimer(ts,dtime)
         tsort = tsort + ts 
         end subroutine imdpsortp2y
!
         subroutine imppush2zf(part,npp,dt,ek,tpush,nx,ny,ipbc)
! multi-tasking particle push with no forces, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         real :: dt, ek, tpush
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
         call MPPUSH2ZF(part,npp,nps,dt,ek,idimp,npmax,nblok,nx,ny,ipbc,&
     &ekp,idtask,nmt,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine imppush2zf
!
         subroutine impgcjpost2(part,fxy,npp,noff,cu,qm,qbm,dt,tdcjpost,&
     &inorder)
! multi-tasking current density deposit with 2d electrostatic fields,
! 1d partition
         implicit none
         real :: qm, qbm, dt, tdcjpost
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
         real, dimension(:,:,:,:,:), allocatable, save :: cup
         integer, save :: szbuf = 0
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
            call MPGCJPOST2L(part,fxy,npp,nps,noff,cu,qm,qbm,dt,idimp,np&
     &max,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         else
            call MPGCJPOST2(part,fxy,npp,nps,noff,cu,qm,qbm,dt,idimp,npm&
     &ax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine impgcjpost2
!
      end module mppush2d
