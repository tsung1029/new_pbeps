!-----------------------------------------------------------------------
!
      module ppush2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library ppush2lib.f
! ppush2mod.f contains interface procedures to process particles:
!             defines module ppush2d
! dpost => ipgpost2 deposits charge density, with various interpolations
!          and optimizations.
!          calls PGPOST2, PGSPOST2, PGSOST2X, PGPOST2L, PGSPOST2L, or
!          PGSOST2XL
! push => ipgpush2 push particles, with various interpolations and
!         optimizations.
!         calls PGPUSH2, PGSPUSH2, PGPUSH2L, or PGSPUSH2L
! pushzf => ippush2zf, push particles with no forces.
!           calls PPUSH2ZF
! sortp => ipsortp2y sorts particles by y grid using memory-conserving
!          bin sort, with various interpolations.
!          calls PSORTP2Y, or PSORTP2YL
! sortp => ipdsortp2y sorts particles by y grid using optimized bin sort
!          with various interpolations.
!          calls PDSORTP2Y, or PDSORTP2YL
! countp => ipcount2y counts number of particles per cell in y grid.
!           calls PCOUNT2YL
! prmove => iprmove2, removes particles which would normally be
!           reflected.
!           calls PRMOVE2
! gcjpost => ipgcjpost2 deposits time-centered current density with
!            2d electrostatic fields.
!            calls PGCJPOST2, or PGCJPOST2L
! initmomt2 calculates local initial momentum for each processor, for 2
!           or 2-1/2d code.
! premoment2 prints out electron and field momentum, calculates total
!            momentum, for 2 or 2-1/2d code.  assumes values are local
!            to this processor.
! primoment2 prints out ion momentum, adds total momentum, for 2 or
!            2-1/2d code.
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 10, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer
      public :: dpost, push, sortp, countp, prmove, pushzf, gcjpost
      public :: initmomt2, premoment2, primoment2
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PDOST2(part,q,npp,noff,qm,nx,idimp,npmax,nblok,nxv,n&
     &ypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGPOST2(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,nyp&
     &mx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSPOST2(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,nx&
     &yp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxyp,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PSOST2X(part,q,npp,noff,nn,amxy,qm,nx,idimp,npmax,nb&
     &lok,nxv,nxvyp,npd,nine)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nxvyp, npd, nine
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxvyp,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         integer, dimension(nine,npd,nblok) :: nn
         real, dimension(nine,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGSOST2X(part,q,npp,noff,nn,amxy,qm,idimp,npmax,nblo&
     &k,nxv,nxvyp,npd,nine)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxvyp, npd, nine
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxvyp,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         integer, dimension(nine,npd,nblok) :: nn
         real, dimension(nine,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PDOST2L(part,q,npp,noff,qm,nx,idimp,npmax,nblok,nxv,&
     &nypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGPOST2L(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,ny&
     &pmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSPOST2L(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,n&
     &xyp)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxyp
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxyp,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PSOST2XL(part,q,npp,noff,nn,amxy,qm,nx,idimp,npmax,n&
     &blok,nxv,nxvyp,npd,ifour)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nxvyp, npd, ifour
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxvyp,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         integer, dimension(ifour,npd,nblok) :: nn
         real, dimension(ifour,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PGSOST2XL(part,q,npp,noff,nn,amxy,qm,idimp,npmax,nbl&
     &ok,nxv,nxvyp,npd,ifour)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nxvyp, npd, ifour
         real :: qm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxvyp,nblok) :: q
         integer, dimension(nblok) :: npp, noff
         integer, dimension(ifour,npd,nblok) :: nn
         real, dimension(ifour,npd,nblok) :: amxy
         end subroutine
      end interface
      interface
         subroutine PPUSH2(part,fx,fy,npp,noff,qbm,dt,ek,nx,idimp,npmax,&
     &nblok,nxv,nypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGPUSH2(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npma&
     &x,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSPUSH2(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npm&
     &ax,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PPUSH2L(part,fx,fy,npp,noff,qbm,dt,ek,nx,idimp,npmax&
     &,nblok,nxv,nypmx)
         implicit none
         integer :: nx, idimp, npmax, nblok, nxv, nypmx
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok) :: fx, fy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGPUSH2L(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGSPUSH2L(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,np&
     &max,nblok,nxv,nxyp,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxyp,nblok) :: fxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PSORTP2Y(part,pt,ip,npic,npp,noff,nyp,idimp,npmax,nb&
     &lok,nypm1)
         implicit none
         integer :: idimp, npmax, nblok, nypm1
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(npmax,nblok) :: pt
         integer, dimension(npmax,nblok) :: ip
         integer, dimension(nypm1,nblok) :: npic
         integer, dimension(nblok) :: npp, noff, nyp
         end subroutine
      end interface
      interface
         subroutine PSORTP2YL(part,pt,ip,npic,npp,noff,nyp,idimp,npmax,n&
     &blok,nypm1)
         implicit none
         integer :: idimp, npmax, nblok, nypm1
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(npmax,nblok) :: pt
         integer, dimension(npmax,nblok) :: ip
         integer, dimension(nypm1,nblok) :: npic
         integer, dimension(nblok) :: npp, noff, nyp
         end subroutine
      end interface
      interface
         subroutine PDSORTP2Y(parta,partb,npic,npp,noff,nyp,idimp,npmax,&
     &nblok,nypm1)
         implicit none
         integer :: idimp, npmax, nblok, nypm1
         real, dimension(idimp,npmax,nblok) :: parta, partb
         integer, dimension(nypm1,nblok) :: npic
         integer, dimension(nblok) :: npp, noff, nyp
         end subroutine
      end interface
      interface
         subroutine PDSORTP2YL(parta,partb,npic,npp,noff,nyp,idimp,npmax&
     &,nblok,nypm1)
         implicit none
         integer :: idimp, npmax, nblok, nypm1
         real, dimension(idimp,npmax,nblok) :: parta, partb
         integer, dimension(nypm1,nblok) :: npic
         integer, dimension(nblok) :: npp, noff, nyp
         end subroutine
      end interface
      interface
         subroutine PCOUNT2YL(part,isign,npic,npp,noff,nyp,idimp,npmax,n&
     &blok,nypm1)
         implicit none
         integer :: isign, idimp, npmax, nblok, nypm1
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nypm1,nblok) :: npic
         integer, dimension(nblok) :: npp, noff, nyp
         end subroutine
      end interface
      interface
         subroutine PRMOVE2(part,npp,ihole,jss,nx,ny,idimp,npmax,nblok,i&
     &dps,ntmax,ipbc)
         implicit none
         integer nx, ny, idimp, npmax, nblok, idps, ntmax, ipbc
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(idps,nblok) :: jss
         end subroutine
      end interface
      interface
         subroutine PPUSH2ZF(part,npp,dt,ek,idimp,npmax,nblok,nx,ny,ipbc&
     &)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, ipbc
         real :: dt, ek
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PGCJPOST2(part,fxy,npp,noff,cu,qm,qbm,dt,idimp,npmax&
     &,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGCJPOST2L(part,fxy,npp,noff,cu,qm,qbm,dt,idimp,npma&
     &x,nblok,nxv,nypmx)
         implicit none
         integer :: idimp, npmax, nblok, nxv, nypmx
         real :: qm, qbm, dt
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2,nxv,nypmx,nblok) :: fxy, cu
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dpost
         module procedure ipgpost2
      end interface
!
      interface push
         module procedure ipgpush2
      end interface
!
      interface pushzf
         module procedure ippush2zf
      end interface
!
      interface sortp
         module procedure ipsortp2y
         module procedure ipdsortp2y
      end interface
!
      interface countp
         module procedure ipcount2y
      end interface
!
      interface prmove
         module procedure iprmove2
      end interface
!
      interface gcjpost
         module procedure ipgcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!
         subroutine ipgpost2(part,q,qm,npp,noff,tdpost,inorder,dopt)
! deposit charge, 1d partition
         implicit none
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, order, opt
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128, ifour = 4, nine = 9
         integer, dimension(nine,npd,size(part,3)) :: nn
         real, dimension(nine,npd,size(part,3)) :: amxy
         real :: td
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(q,1); nypmx = size(q,2); nxyp = nxv*nypmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(dopt)) opt = dopt
! initialize timer
         call wtimer(td,dtime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call PGSPOST2L(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,n&
     &xyp)
            else if (opt==VECTOR) then
               call PGSOST2XL(part,q,npp,noff,nn,amxy,qm,idimp,npmax,nbl&
     &ok,nxv,nxyp,npd,ifour)
            else
               call PGPOST2L(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,ny&
     &pmx)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSPOST2(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,nx&
     &yp)
            else if (opt==VECTOR) then
               call PGSOST2X(part,q,npp,noff,nn,amxy,qm,idimp,npmax,nblo&
     &k,nxv,nxyp,npd,nine)
            else
               call PGPOST2(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,nyp&
     &mx)
            endif
         endif
! record time
         call wtimer(td,dtime)
         tdpost = tdpost + td
         end subroutine ipgpost2
!
         subroutine ipgpush2(part,fxy,npp,noff,qbm,dt,ek,tpush,nx,ny,ipb&
     &c,inorder,popt)
! push particles with 2d electrostatic fields, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
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
               call PGSPUSH2L(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,np&
     &max,nblok,nxv,nxyp,ipbc)
            else
               call PGPUSH2L(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc)
            endif
         else
            if (opt==LOOKAHEAD) then
               call PGSPUSH2(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npm&
     &ax,nblok,nxv,nxyp,ipbc)
            else
               call PGPUSH2(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npma&
     &x,nblok,nxv,nypmx,ipbc)
            endif
         endif
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ipgpush2
!
         subroutine ipsortp2y(part,pt,ip,npp,noff,nyp,npic,tsort,inorder&
     &)
! sort particles by y grid using memory-conserving bin sort
         implicit none
         real :: tsort
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: pt
         integer, dimension(:,:), pointer :: ip, npic
         integer, dimension(:), pointer :: npp, noff, nyp
         integer, optional :: inorder
! local data
         integer :: idimp, npmax, nblok, nypm1, order
         real :: ts
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); nypm1 = size(npic,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,dtime,-1)
         if (order==LINEAR) then
            call PSORTP2YL(part,pt,ip,npic,npp,noff,nyp,idimp,npmax,nblo&
     &k,nypm1)
         else
            call PSORTP2Y(part,pt,ip,npic,npp,noff,nyp,idimp,npmax,nblok&
     &,nypm1)
         endif
! record time
         call wtimer(ts,dtime)
         tsort = tsort + ts
         end subroutine ipsortp2y
!
         subroutine ipdsortp2y(parta,partb,npp,noff,nyp,npic,tsort,inord&
     &er)
! sort particles by y grid using optimized bin sort
         implicit none
         real :: tsort
         real, dimension(:,:,:), pointer :: parta, partb
         integer, dimension(:), pointer :: npp, noff, nyp
         integer, dimension(:,:), pointer :: npic
         integer, optional :: inorder
! local data
         real, dimension(:,:,:), pointer :: part
         integer :: idimp, npmax, nblok, nypm1, order
         real :: ts
         double precision :: dtime
         idimp = size(parta,1); npmax = size(parta,2)
         nblok = size(parta,3); nypm1 = size(npic,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,dtime,-1)
         if (order==LINEAR) then
            call PDSORTP2YL(parta,partb,npic,npp,noff,nyp,idimp,npmax,nb&
     &lok,nypm1)
         else
            call PDSORTP2Y(parta,partb,npic,npp,noff,nyp,idimp,npmax,nbl&
     &ok,nypm1)
         endif
         part => parta
         parta => partb
         partb => part
! record time
         call wtimer(ts,dtime)
         tsort = tsort + ts 
         end subroutine ipdsortp2y
!
         subroutine ipcount2y(part,npic,npp,noff,nyp)
! counts number of particles per cell in y grid
         implicit none
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:,:), pointer :: npic
         integer, dimension(:), pointer :: npp, noff, nyp
! local data
         integer :: isign = 1
         integer :: idimp, npmax, nblok, nypm1
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); nypm1 = size(npic,1)
         call PCOUNT2YL(part,isign,npic,npp,noff,nyp,idimp,npmax,nblok,n&
     &ypm1)
         end subroutine ipcount2y
!
         subroutine iprmove2(part,npp,nx,ny,nbmax,idps,ipbc)
! removes particles which would normally be reflected
         implicit none
         integer :: nx, ny, nbmax, idps, ipbc
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(2*nbmax,size(part,3)) :: ihole
         integer, dimension(idps,size(part,3)) :: jss
         integer :: idimp, npmax, nblok, ntmax
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         ntmax = 2*nbmax
         call PRMOVE2(part,npp,ihole,jss,nx,ny,idimp,npmax,nblok,idps,nt&
     &max,ipbc)
         end subroutine iprmove2
!
         subroutine ippush2zf(part,npp,dt,ek,tpush,nx,ny,ipbc)
! push particles with no forces, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         real :: dt, ek, tpush
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
         call PPUSH2ZF(part,npp,dt,ek,idimp,npmax,nblok,nx,ny,ipbc)
! record time
         call wtimer(tp,dtime)
         tpush = tpush + tp
         end subroutine ippush2zf
!
         subroutine ipgcjpost2(part,fxy,npp,noff,cu,qm,qbm,dt,tdcjpost,i&
     &norder)
! deposit current density with 2d electrostatic fields, 1d partition
         implicit none
         real :: qm, qbm, dt, tdcjpost
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
            call PGCJPOST2L(part,fxy,npp,noff,cu,qm,qbm,dt,idimp,npmax,n&
     &blok,nxv,nypmx)
         else
            call PGCJPOST2(part,fxy,npp,noff,cu,qm,qbm,dt,idimp,npmax,nb&
     &lok,nxv,nypmx)
         endif
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine ipgcjpost2
!
         subroutine initmomt2(part,npp,px,py,pz,ndim)
! calculate local initial momentum for each processor,
! for 2 or 2-1/2d code
         real :: px, py, pz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:) :: npp
         integer, optional :: ndim
! local data
         integer :: j, l, nd
         double precision :: sum1, sum2, sum3
         nd = 3
         if (present(ndim)) nd = ndim
! calculate momentum at t=t-dt/2
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum3 = 0.0d0
         select case(nd)
         case (2)
            do l = 1, size(part,3)
            do j = 1, npp(l)
            sum1 = sum1 + part(3,j,l)
            sum2 = sum2 + part(4,j,l)
            enddo
            enddo
            px = sum1
            py = sum2
            pz = 0.0
         case (3)
            do l = 1, size(part,3)
            do j = 1, npp(l)
            sum1 = sum1 + part(3,j,l)
            sum2 = sum2 + part(4,j,l)
            sum3 = sum3 + part(5,j,l)
            enddo
            enddo
            px = sum1
            py = sum2
            pz = sum3
         end select
         end subroutine initmomt2
!
         subroutine premoment2(part,itime,npp,msg,id0,iunit,px,py,pz,sx,&
     &sy,sz,wx,wy,wz,ndim,nprint)
! print out electron and field momentum, calculate total momentum.
! assumes values of px, py, pz, sx, sy, sz are local to this processor
! for 2 or 2-1/2d code
         integer :: itime, id0, iunit
         real :: px, py, pz, sx, sy, sz, wx, wy, wz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
         double precision, dimension(:) :: msg
         integer, optional :: ndim, nprint
! local data
         integer :: j, l, nd, np
         real :: tx, ty, tz
         double precision :: sum1, sum2, sum3
         double precision, dimension(6) :: sum6, work6
  991    format (' T = ',i7)
  994    format (' electron momentum = ',3e14.7)
  996    format (' field momentum = ',3e14.7)
         nd = 3
         if (present(ndim)) nd = ndim
         np = 1
         if (present(nprint)) np = nprint
         if (np < 0) return
         if (np==0) then
            px = 0.0; py = 0.0; pz = 0.0
         else if (np==1) then
            if (id0==0) write (iunit,991) itime
         endif
! calculate and print electron momentum
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum3 = 0.0d0
         select case(nd)
         case (2)
            do l = 1, size(part,3)
            do j = 1, npp(l)
            sum1 = sum1 + part(3,j,l)
            sum2 = sum2 + part(4,j,l)
            enddo
            enddo
            sum6(1) = 0.5*(px + sum1)
            sum6(2) = 0.5*(py + sum2)
            sum6(3) = 0.0
         case (3)
            do l = 1, size(part,3)
            do j = 1, npp(l)
            sum1 = sum1 + part(3,j,l)
            sum2 = sum2 + part(4,j,l)
            sum3 = sum3 + part(5,j,l)
            enddo
            enddo
            sum6(1) = 0.5*(px + sum1)
            sum6(2) = 0.5*(py + sum2)
            sum6(3) = 0.5*(pz + sum3)
         end select
         sum6(4) = sx
         sum6(5) = sy
         sum6(6) = sz
         call PDSUM(sum6,work6,6,1)
         px = sum6(1)
         py = sum6(2)
         pz = sum6(3)
         tx = sum6(4)
         ty = sum6(5)
         tz = sum6(6)
         if (np==1) then
! save momentum values
            msg(1:6) = sum6
            if (id0==0) then
! print electron momentum
               write (iunit,994) px, py, pz
! print field momentum
               write (iunit,996) tx, ty, tz
            endif
! calculate total momentum
            wx = px + tx
            wy = py + ty
            wz = pz + tz
         endif
         px = sum1
         py = sum2
         pz = sum3
         end subroutine premoment2
!
         subroutine primoment2(parti,nppi,msg,id0,iunit,rmass,px,py,pz,w&
     &x,wy,wz,ndim,nprint)
! print out ion momentum, adds total momentum, for 2 or 2-1/2d code
         integer :: id0, iunit
         real :: rmass, px, py, pz, wx, wy, wz
         real, dimension(:,:,:), pointer :: parti
         integer, dimension(:), pointer :: nppi
         double precision, dimension(:) :: msg
         integer, optional :: ndim, nprint
! local data
         integer :: j, l, nd, np
         real :: at1
         double precision :: sum0, sum1, sum2
         double precision, dimension(3) :: sum3, work3
  995    format (' ion momentum = ',3e14.7)
         nd = 3
         if (present(ndim)) nd = ndim
         np = 1
         if (present(nprint)) np = nprint
         if (np < 0) return
         at1 = 0.5*rmass
         if (np==0) then
            px = 0.0; py = 0.0; pz = 0.0
         endif
! calculate and print ion momentum
         sum0 = 0.0d0
         sum1 = 0.0d0
         sum2 = 0.0d0
         select case(nd)
         case (2)
            do l = 1, size(parti,3)
            do j = 1, nppi(l)
            sum0 = sum0 + parti(3,j,l)
            sum1 = sum1 + parti(4,j,l)
            enddo
            enddo
            sum3(1) = at1*(px + sum0)
            sum3(2) = at1*(py + sum1)
            sum3(3) = 0.0
         case (3)
            do l = 1, size(parti,3)
            do j = 1, nppi(l)
            sum0 = sum0 + parti(3,j,l)
            sum1 = sum1 + parti(4,j,l)
            sum2 = sum2 + parti(5,j,l)
            enddo
            enddo
            sum3(1) = at1*(px + sum0)
            sum3(2) = at1*(py + sum1)
            sum3(3) = at1*(pz + sum2)
         end select
         call PDSUM(sum3,work3,3,1)
         px = sum3(1)
         py = sum3(2)
         pz = sum3(3)
         if (np==1) then
! print ion momentum
            if (id0==0) write (iunit,995) px, py, pz
! add to total momentum
            wx = wx + px
            wy = wy + py
            wz = wz + pz
! save momentum values
            msg(1:3) = sum3
            msg(4:6) = (/wx,wy,wz/)
         endif
         px = sum0
         py = sum1
         pz = sum2
         end subroutine primoment2
!
      end module ppush2d
