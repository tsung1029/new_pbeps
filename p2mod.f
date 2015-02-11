!-----------------------------------------------------------------------
!
      module p2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library p2lib.f
! p2mod.f contains interface procedures for communications with 1d
!         partitions:
!         defines module p2d
! dcomp => idcomp2 finds uniform 1d partition boundaries in 2d code.
!          calls DCOMP2, or DCOMP2L
! pmove => iwpmove2 moves particles to appropriate processor non-uniform
!          1d partition boundaries in 2d code.
!          calls WPMOVE2, or WPXMOV2
! pmove => iwpdmove2 moves particles to appropriate processor
!          non-uniform 1d partition boundaries in 2d code, and returns
!          load imbalance.
!          calls WPMOVE2, or WPXMOV2
! pmoves => iwpmoves2 moves particles to appropriate processor
!           non-uniform 1d partition boundaries in 2d code.  Optimized
!           requiring that maximum number of particle passes must be
!           known.
!           calls WPMOVES2, or WPXMOVS2
! pcguard => ipcguard2 copies guard cells in y for uniform, periodic 2d
!            vector data.
!            calls PCGUARD2, or PCGUARD2L
! pcguard => ipdguard2 copies guard cells in y for uniform, periodic 2d
!            scalar data.
!            calls PCGUARD2, or PCGUARD2L
! pcguardp => ipcguard2p copies guard cells in y for uniform,
!             non-periodic 2d vector data.
!             calls PCGUARD2, PCGUARD2L, PLBGUARD2, or PLCGUARD2
! pcguardp => ipdguard2p copies guard cells in y for uniform.
!             non-periodic 2d scalar data.
!             calls PCGUARD2, PCGUARD2L, PLDGUARD2, or PLCGUARD2L
! pncguardp => ipncguard2p copies guard cells in y for non-uniform,
!              non-periodic 2d vector data.
!              calls PNCGUARD2, PNCGUARD2L, PNLBGUARD2, PNLCGUARD2, or
!              PNLCGUARD2L
! pncguardp => ipndguard2p copies guard cells in y for non-uniform,
!              non-periodic 2d scalar data.
!              calls PNCGUARD2, PNCGUARD2L, PNLDGUARD2, or PNLCGUARD2L
! paguard => ipaguard2 add2 guard cells in y for uniform, periodic 2d
!            scalar data.
!            calls PAGUARD2, or PAGUARD2L
! paguard => ipacguard2 adds guard cells in y for uniform, periodic 2d
!            vector data.
!            calls PACGUARD2, PACGUARD22, PACGUARD2L, or PACGUARD22L
! pamcguard => ipamcguard2 adds guard cells in y for uniform, periodic
!              2d tensor data.
!              calls PAMCGUARD2, or PAMCGUARD2L
! paguardp => ipaguard2p adds guard cells in y for uniform, non-periodic
!             2d scalar data.
!             calls PAGUARD2, PAGUARD2L, PLAGUARD2, PLAGUARD2, or
!             PLAGUARD2L
! paguardp => ipacguard2p adds guard cells in y for uniform,
!             non-periodic 2d vector data.
!             calls PACGUARD2, PACGUARD22, PACGUARD2L, PACGUARD22L,
!             PLACGUARD2, PLACGUARDS2, PLACGUARD22, PLACGUARDS22,
!             PLACGUARD2L, or PLACGUARD22L
! pnaguardp => ipnaguard2p adds guard cells in y for non-uniform,
!              non-periodic 2d scalar data.
!              calls PNAGUARD2, PNAGUARD2L, PNLAGUARD2, PNLAGUARDS2, or
!              PNLAGUARD2L
! pnaguardp => ipnacguard2p adds guard cells in y for non-uniform,
!              non-periodic 2d vector data.
!              calls PNACGUARD2, PNACGUARD22, PNACGUARD2L, PNACGUARD22L,
!              PNLACGUARD2, PNLACGUARDS2, PNLACGUARD22, PNLACGUARDS22,
!              PNLACGUARD2L, or PNLACGUARD22L
! pfmove => ipfmove2 moves 2d scalar data between uniform and
!           non-uniform 1d partitions.
!           calls PFMOVE2
! pfmove => ipfcmove2 moves 2d vector data between uniform and
!           non-uniform 1d partitions.
!           calls PFMOVE2
! pfmove => ipnfmove2 moves 2d scalar data between two different
!           non-uniform 1d partitions.
!           calls PFMOVE2
! repart => irepartd2 finds new 1d partitions from old partition and
!           particle information.
!           calls REPARTD2
! fnoff => ifnoff2 finds new 1d partitions arrays from edges.
!          calls FNOFF2
! dblsin => ipdblsin2c doubles array in each dimension for 2d vector
!           data for dirichlet boundary conditions.
!           calls PDBLSIN2B, or PDBLSIN2C
! dblsin => ipdblsin2d doubles array in each dimension for 2d scalar
!           data for dirichlet boundary conditions.
!           calls PDBLSIN2D
! dblcos => ipdblcos2c doubles array in each dimension for 2d vector
!           data for neumann boundary conditions.
!           calls PDBLCOS2B, or PDBLCOS2C
! dblcos => ipdblcos2d doubles array in each dimension for 2d scalar
!           data for neumann boundary conditions.
!           calls PDBLCOS2D
! hafdbl => iphafdbl2c copies 2d vector data from double to normal array
!           in each dimension.
!           calls PHAFDBL2D
! hafdbl => iphafdbl2d copies 2d scalar data from double to normal array
!           in each dimension.
!           calls PHAFDBL2B, or PHAFDBL2C
! zdbl => izpdbl2c doubles array in each dimension for 2d vector data,
!         zeroing copies for open boundary conditions.
!         calls PZDBL2B, or PZDBL2C
! zdbl => izpdbl2d doubles array in each dimension for 2d scalar data,
!         zeroing copies for open boundary conditions.
!         calls PZDBL2D
! plsum => ipsum2 perform2 global sum of 2d real array.
!          calls PSUM
! plbcast => ipbcast2 broadcasts 2d real array.
!            calls PBCAST
! writebf => ipwrite2 collects a subset of a distributed real 2d scalar
!            array and writes it to a direct access binary file.
!            calls PWRITE2
! writebf => ipcwrite2 collects a subset of a distributed complex 2d
!            scalar array and writes it to a direct access binary file.
!            calls PCWRITE2
! writebf => ipvcwrite2 collects a subset of a distributed complex 2d
!            vector array and writes it to a direct access binary file.
!            calls PCWRITE2
! readbf => ipread2 reads a subset of a distributed real 2d scalar array
!           from a direct access binary file and distributes it.
!           calls PREAD2
! readbf => ipcread2 reads a subset of a distributed complex 2d scalar
!           array from a direct access binary file and distributes it.
!           calls PCREAD2
! readbf => ipvcread2 reads a subset of a distributed complex 2d vector
!           array from a direct access binary file and distributes it.
!           calls PCREAD2
! wrdata => ipwrdata2 collects distributed real 2d scalar data and
!           writes it to a fortran unformatted sequential file.
!           calls PWRDATA
! wrdata => ipwrrdata2 collects distributed real 2d vector data and
!           writes it to a fortran unformatted sequential file.
!           calls PWRDATA
! wrdata => ipwrcdata2 collects distributed complex 2d vector data and
!           writes it to a fortran unformatted sequential file.
!           calls PWRDATA
! rddata => iprddata2 reads real 2d scalar data from a fortran
!           unformatted sequential file and distributes it, for uniform
!           1d partitions
!           calls PRDDATA
! rddata => iprdrdata2 reads real 2d vector data from a fortran
!           unformatted sequential file and distributes it.
!           calls PWRDATA
! rddata => iprdcdata2 reads complex 2d vector data from a fortran
!           unformatted sequential file, and distributes it.
!           calls PRDDATA
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 5, 2009
!
      use p0d
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PPINIT, PPID, PPEXIT, HARTBEAT
      public :: get_funit, pwtimer, plsum, plmax, plscan, plbcast
      public :: writebf, readbf, wtimer, wrdata, rddata
      public :: dcomp, pmove, pcguard, pcguardp
      public :: paguard, pamcguard, paguardp, pncguardp, pnaguardp
      public :: pfmove, repart
      public :: fnoff, dblsin, dblcos, hafdbl, zdbl, pmoves
!
! buffer data for particle managers
      real, dimension(:,:,:), allocatable :: sbufl, sbufr, rbufl, rbufr
      integer, dimension(:,:), allocatable :: ihole
      integer :: szbuf = 0
      save
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine DCOMP2(edges,nyp,noff,ny,kstrt,nvp,idps,nblok)
         implicit none
         integer :: ny, kstrt, nvp, idps, nblok
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: nyp, noff
         end subroutine
      end interface
      interface
         subroutine DCOMP2L(edges,nyp,noff,ny,kstrt,nvp,idps,nblok)
         implicit none
         integer :: ny, kstrt, nvp, idps, nblok
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: nyp, noff
         end subroutine
      end interface
      interface
         subroutine PCGUARD2(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nxv, nypmx, kyp, kblok
         real, dimension(nxv,nypmx,kblok) :: f
         end subroutine
      end interface
      interface
         subroutine PNCGUARD2(f,scs,nyp,kstrt,nvp,nxv,nypmx,nblok,mter)
         implicit none
         integer :: kstrt, nvp, nxv, nypmx, nblok, mter
         real, dimension(nxv,nypmx,nblok) :: f
         real, dimension(nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nxv, nypmx, kyp, kblok
         real, dimension(nxv,nypmx,kblok) :: f
         end subroutine
      end interface
      interface
         subroutine PNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
         implicit none
         integer :: kstrt, nvp, nxv, nypmx, nblok
         real, dimension(nxv,nypmx,nblok) :: f
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngd&
     &s)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
         real, dimension(3,nxv,nypmx,kblok) :: f
         real, dimension(3,nxv,ngds,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PACGUARD22(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ng&
     &ds)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
         real, dimension(2,nxv,nypmx,kblok) :: f
         real, dimension(2,nxv,ngds,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PAGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds&
     &)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
         real, dimension(nxv,nypmx,kblok) :: f
         real, dimension(nxv,ngds,kblok) :: scr
         end subroutine
      end interface
      interface
        subroutine PAMCGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngd&
     &s,ndim)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds, ndim
         real, dimension(ndim,nxv,nypmx,kblok) :: f
         real, dimension(ndim,nxv,ngds,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNACGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblo&
     &k,ngds,mter)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
         real, dimension(3,nxv,nypmx,nblok) :: f
         real, dimension(3,nxv,ngds,nblok) :: scr
         real, dimension(3,nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNACGUARD22(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nbl&
     &ok,ngds,mter)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
         real, dimension(2,nxv,nypmx,nblok) :: f
         real, dimension(2,nxv,ngds,nblok) :: scr
         real, dimension(2,nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNAGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok&
     &,ngds,mter)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
         real, dimension(nxv,nypmx,nblok) :: f
         real, dimension(nxv,ngds,nblok) :: scr
         real, dimension(nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNAMCGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nbl&
     &ok,ngds,ndim,mter)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, ndim, mter
         real, dimension(ndim,nxv,nypmx,nblok) :: f
         real, dimension(ndim,nxv,ngds,nblok) :: scr
         real, dimension(ndim,nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(3,nxv,nypmx,kblok) :: f
         real, dimension(3,nxv,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PACGUARD22L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(2,nxv,nypmx,kblok) :: f
         real, dimension(2,nxv,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PAGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(nxv,nypmx,kblok) :: f
         real, dimension(nxv,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PAMCGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,n&
     &dim)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ndim
         real, dimension(ndim,nxv,nypmx,kblok) :: f
         real, dimension(ndim,nxv,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNACGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, nblok
         real, dimension(3,nxv,nypmx,nblok) :: f
         real, dimension(3,nxv,nblok) :: scr
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNACGUARD22L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, nblok
         real, dimension(2,nxv,nypmx,nblok) :: f
         real, dimension(2,nxv,nblok) :: scr
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNAGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, nblok
         real, dimension(nxv,nypmx,nblok) :: f
         real, dimension(nxv,nblok) :: scr
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNAMCGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,&
     &ndim)
         implicit none
         integer kstrt, nvp, nx, nxv, nypmx, nblok, ndim
         real, dimension(ndim,nxv,nypmx,nblok) :: f
         real, dimension(ndim,nxv,nblok) :: scr
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PDBLSIN2C(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok&
     &,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: cu
         real, dimension(2,2*nxv,kyp2,k2blok) :: cu2
         end subroutine
      end interface
      interface
         subroutine PDBLSIN2D(q,q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k&
     &2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: q
         real, dimension(2*nxv,kyp2,k2blok) :: q2
         end subroutine
      end interface
      interface
         subroutine PDBLSIN2B(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok&
     &,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: cu
         real, dimension(3,2*nxv,kyp2,k2blok) :: cu2
         end subroutine
      end interface
      interface
         subroutine PDBLCOS2C(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok&
     &,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: cu
         real, dimension(2,2*nxv,kyp2,k2blok) :: cu2
         end subroutine
      end interface
      interface
         subroutine PDBLCOS2D(q,q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k&
     &2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: q
         real, dimension(2*nxv,kyp2,k2blok) :: q2
         end subroutine
      end interface
      interface
         subroutine PDBLCOS2B(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok&
     &,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: cu
         real, dimension(3,2*nxv,kyp2,k2blok) :: cu2
         end subroutine
      end interface
      interface
         subroutine PHAFDBL2C(fxy,fxy2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kbl&
     &ok,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: fxy
         real, dimension(2,2*nxv,kyp2,k2blok) :: fxy2
         end subroutine
      end interface
      interface
         subroutine PHAFDBL2D(q,q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k&
     &2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: q
         real, dimension(2*nxv,kyp2,k2blok) :: q2
         end subroutine
      end interface
      interface
         subroutine PHAFDBL2B(fxy,fxy2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kbl&
     &ok,k2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: fxy
         real, dimension(3,2*nxv,kyp2,k2blok) :: fxy2
         end subroutine
      end interface
      interface
         subroutine PLCGUARD2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(2,nxv,nypmx,kblok) :: f
         end subroutine
      end interface
      interface
         subroutine PLDGUARD2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(nxv,nypmx,kblok) :: f
         end subroutine
      end interface
      interface
         subroutine PLBGUARD2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(3,nxv,nypmx,kblok) :: f
         end subroutine
      end interface
      interface
         subroutine PNLCGUARD2(f,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mt&
     &er)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok, mter
         real, dimension(2,nxv,nypmx,nblok) :: f
         real, dimension(2,nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNLDGUARD2(f,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mt&
     &er)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok, mter
         real, dimension(nxv,nypmx,nblok) :: f
         real, dimension(nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNLBGUARD2(f,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mt&
     &er)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok, mter
         real, dimension(3,nxv,nypmx,nblok) :: f
         real, dimension(3,nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nxv, nypmx, kyp, kblok
         real, dimension(nxv,nypmx,kblok) :: f
         end subroutine
      end interface
      interface
         subroutine PNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
         implicit none
         integer :: kstrt, nvp, nxv, nypmx, nblok
         real, dimension(nxv,nypmx,nblok) :: f
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLACGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ng&
     &ds)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
         real, dimension(3,nxv,nypmx,kblok) :: f
         real, dimension(3,nxv,ngds,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PLACGUARD22(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,n&
     &gds)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
         real, dimension(2,nxv,nypmx,kblok) :: f
         real, dimension(2,nxv,ngds,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PLAGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngd&
     &s)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
         real, dimension(nxv,nypmx,kblok) :: f
         real, dimension(nxv,ngds,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNLACGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nbl&
     &ok,ngds,mter)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
         real, dimension(3,nxv,nypmx,nblok) :: f
         real, dimension(3,nxv,ngds,nblok) :: scr
         real, dimension(3,nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNLACGUARD22(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nb&
     &lok,ngds,mter)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
         real, dimension(2,nxv,nypmx,nblok) :: f
         real, dimension(2,nxv,ngds,nblok) :: scr
         real, dimension(2,nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNLAGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblo&
     &k,ngds,mter)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
         real, dimension(3,nxv,nypmx,nblok) :: f
         real, dimension(3,nxv,ngds,nblok) :: scr
         real, dimension(3,nxv,nblok) :: scs
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLACGUARDS2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(3,nxv,nypmx,kblok) :: f
         end subroutine
      end interface
      interface
         subroutine PLACGUARDS22(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(2,nxv,nypmx,kblok) :: f
         end subroutine
      end interface
      interface
         subroutine PLAGUARDS2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(nxv,nypmx,kblok) :: f
         end subroutine
      end interface
      interface
         subroutine PNLACGUARDS2(f,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mter&
     &)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok, mter
         real, dimension(3,nxv,nypmx,nblok) :: f
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNLACGUARDS22(f,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mte&
     &r)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok, mter
         real, dimension(2,nxv,nypmx,nblok) :: f
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNLAGUARDS2(f,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mter)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok, mter
         real, dimension(nxv,nypmx,nblok) :: f
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLACGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(3,nxv,nypmx,kblok) :: f
         real, dimension(3,nxv,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PLACGUARD22L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(2,nxv,nypmx,kblok) :: f
         real, dimension(2,nxv,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PLAGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, kyp, kblok
         real, dimension(nxv,nypmx,kblok) :: f
         real, dimension(nxv,kblok) :: scr
         end subroutine
      end interface
      interface
         subroutine PNLACGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok
         real, dimension(3,nxv,nypmx,nblok) :: f
         real, dimension(3,nxv,nblok) :: scr
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNLACGUARD22L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok&
     &)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok
         real, dimension(2,nxv,nypmx,nblok) :: f
         real, dimension(2,nxv,nblok) :: scr
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PNLAGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
         implicit none
         integer :: kstrt, nvp, nx, nxv, nypmx, nblok
         real, dimension(nxv,nypmx,nblok) :: f
         real, dimension(nxv,nblok) :: scr
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PZDBL2C(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k&
     &2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: cu
         real, dimension(2,2*nxv,kyp2,k2blok) :: cu2
         end subroutine
      end interface
      interface
         subroutine PZDBL2D(q,q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2b&
     &lok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: q
         real, dimension(2*nxv,kyp2,k2blok) :: q2
         end subroutine
      end interface
      interface
         subroutine PZDBL2B(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k&
     &2blok)
         implicit none
         integer :: nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
         real :: cu
         real, dimension(3,2*nxv,kyp2,k2blok) :: cu2
         end subroutine
      end interface
      interface
         subroutine PMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,&
     &jsr,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         end subroutine
      end interface
      interface
         subroutine PXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,&
     &jsr,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,maskp,&
     &info)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(npmax,nblok) :: maskp
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         end subroutine
      end interface
      interface
         subroutine WPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole&
     &,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,in&
     &fo)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax
         real :: th
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         end subroutine
      end interface
      interface
         subroutine WPXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole&
     &,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,ma&
     &skp,info)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax
         real :: th
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(npmax,nblok) :: maskp
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         end subroutine
      end interface
      interface
         subroutine WPMOVES2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihol&
     &e,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,i&
     &nfo)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax
         real :: th
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         end subroutine
      end interface
      interface
         subroutine WPXMOVS2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihol&
     &e,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,m&
     &askp,info)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax
         real :: th
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(npmax,nblok) :: maskp
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         end subroutine
      end interface
      interface
         subroutine PFMOVE2(f,g,noff,nyp,noffs,nyps,noffd,nypd,jsr,jsl,i&
     &sign,kyp,kstrt,nvp,nxv,nypmx,nblok,idps,mter,ierr)
         implicit none
         integer :: isign, kyp, kstrt, nvp, nxv, nypmx, nblok, idps
         integer :: mter, ierr
!        real, dimension(*) :: f
         real :: f
         real, dimension(nxv,nypmx,nblok) :: g
         integer, dimension(nblok) :: noff, nyp
         integer, dimension(nblok) :: noffs, nyps, noffd, nypd
         integer, dimension(idps,nblok) :: jsl, jsr
         end subroutine
      end interface
      interface
         subroutine REPARTD2(edges,edg,eds,eg,es,et2,npic,noff,nyp,anpav&
     &,nypmin,nypmax,kstrt,nvp,nblok,idps,nypm)
         implicit none
         integer :: nypmin, nypmax, kstrt, nvp, nblok, idps, nypm
         real :: anpav
         real, dimension(idps,nblok) :: edges
         real, dimension(nypm,nblok) :: edg, eds
         real, dimension(idps,nblok) :: eg, es
         real, dimension(2*idps,nblok) :: et2
         integer, dimension(nypm,nblok) :: npic
         integer, dimension(nblok) :: noff, nyp
         end subroutine
      end interface
      interface
         subroutine FNOFF2(edges,noff,nyp,nypmin,nypmax,nblok,idps)
         implicit none
         integer :: nypmin, nypmax, nblok, idps
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: noff, nyp
         end subroutine
      end interface
      interface
         subroutine PTPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd&
     &,jblok,kblok)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
         integer :: jblok, kblok
         complex, dimension(nxv,kypd,kblok) :: f
         complex, dimension(nyv,kxpd,jblok) :: g
         complex, dimension(kxp,kyp,kblok) :: s
         complex, dimension(kxp,kyp,jblok) :: t
         end subroutine
      end interface
      interface
         subroutine P2TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kyp&
     &d,jblok,kblok)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
         integer :: jblok, kblok
         complex, dimension(2,nxv,kypd,kblok) :: f
         complex, dimension(2,nyv,kxpd,jblok) :: g
         complex, dimension(2,kxp,kyp,kblok) :: s
         complex, dimension(2,kxp,kyp,jblok) :: t
         end subroutine
      end interface
      interface
         subroutine P3TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kyp&
     &d,jblok,kblok)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
         integer :: jblok, kblok
         complex, dimension(3,nxv,kypd,kblok) :: f
         complex, dimension(3,nyv,kxpd,jblok) :: g
         complex, dimension(3,kxp,kyp,kblok) :: s
         complex, dimension(3,kxp,kyp,jblok) :: t
         end subroutine
      end interface
      interface
         subroutine PTPOSEX(f,g,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,jb&
     &lok,kblok)
         implicit none
         integer  :: nx, ny, kstrt, nxv, nyv, kxp, kyp
         integer :: kxpd, kypd, jblok, kblok
         complex, dimension(nxv*kypd*kblok) :: f
         complex, dimension(nyv*kxpd*jblok) :: g
         end subroutine
      end interface
      interface
         subroutine P2TPOSEX(f,g,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j&
     &blok,kblok)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp
         integer :: kxpd, kypd, jblok, kblok
         complex, dimension(2*nxv*kypd*kblok) :: f
         complex, dimension(2*nyv*kxpd*jblok) :: g
         end subroutine
      end interface
      interface
         subroutine P3TPOSEX(f,g,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j&
     &blok,kblok)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp
         integer :: kxpd, kypd, jblok, kblok
         complex, dimension(3*nxv*kypd*kblok) :: f
         complex, dimension(3*nyv*kxpd*jblok) :: g
         end subroutine
      end interface
      interface
         subroutine PNTPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kyp&
     &d,jblok,kblok,ndim)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
         integer :: jblok, kblok, ndim
         complex, dimension(3,nxv,kypd,kblok) :: f
         complex, dimension(3,nyv,kxpd,jblok) :: g
         complex, dimension(3,kxp,kyp,kblok) :: s
         complex, dimension(3,kxp,kyp,jblok) :: t
         end subroutine
      end interface
      interface
         subroutine PNTPOSEX(f,g,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j&
     &blok,kblok,ndim)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp
         integer :: kxpd, kypd, jblok, kblok, ndim
         complex, dimension(3*nxv*kypd*kblok) :: f
         complex, dimension(3*nyv*kxpd*jblok) :: g
         end subroutine
      end interface
      interface
         subroutine PRTPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kyp&
     &d,jblok,kblok)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
         integer :: jblok, kblok
         real, dimension(nxv,kypd,kblok) :: f
         real, dimension(nyv,kxpd,jblok) :: g
         real, dimension(kxp+1,kyp+1,kblok) :: s
         real, dimension(kxp+1,kyp+1,jblok) :: t
         end subroutine
      end interface
      interface
         subroutine PR2TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,ky&
     &pd,jblok,kblok)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
         integer :: jblok, kblok
         real, dimension(2,nxv,kypd,kblok) :: f
         real, dimension(2,nyv,kxpd,jblok) :: g
         real, dimension(2,kxp+1,kyp+1,kblok) :: s
         real, dimension(2,kxp+1,kyp+1,jblok) :: t
         end subroutine
      end interface
      interface
         subroutine PR3TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,ky&
     &pd,jblok,kblok)
         implicit none
         integer :: nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
         integer :: jblok, kblok
         real, dimension(3,nxv,kypd,kblok) :: f
         real, dimension(3,nyv,kxpd,jblok) :: g
         real, dimension(3,kxp+1,kyp+1,kblok) :: s
         real, dimension(3,kxp+1,kyp+1,jblok) :: t
         end subroutine
      end interface
      interface
         subroutine PWRITE2(f,nx,kyp,nxv,kypmx,nblok,iunit,nrec,lrec,nam&
     &e)
         implicit none
         integer :: nx, kyp, nxv, kypmx, nblok, iunit, nrec, lrec
         character(len=*) :: name
!        real, dimension(*) :: f
         real :: f
         end subroutine
      end interface
      interface
         subroutine PREAD2(f,nx,kyp,nxv,kypmx,nblok,iunit,nrec,lrec,name&
     &,ierror)
         implicit none
         integer :: nx, kyp, nxv, kypmx, nblok, iunit, nrec, lrec
         integer :: ierror
         character(len=*) :: name
!        real, dimension(*) :: f
         real :: f
         end subroutine
      end interface
      interface
         subroutine PCWRITE2(f,nx,ny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,&
     &name)
         implicit none
         integer :: nx, ny, kxp, nyv, kxpd, jblok, iunit, nrec, lrec
         character(len=*) :: name
         complex, dimension(nyv,kxpd,jblok) :: f
         end subroutine
      end interface
      interface
         subroutine PCREAD2(f,nx,ny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,n&
     &ame,ierror)
         implicit none
         integer :: nx, ny, kxp, nyv, kxpd, jblok, iunit, nrec, lrec
         integer :: ierror
         character(len=*) :: name
         complex, dimension(nyv,kxpd,jblok) :: f
         end subroutine
      end interface
      interface
         subroutine PCDIFF2(f,g,kstrt,nvp,nx,nyp,nxv,nypd,nyp1,nblok)
         implicit none
         integer :: kstrt, nvp, nx, nyp, nxv, nypd, nyp1, nblok
         real, dimension(nxv,nypd,nblok) :: f
         real, dimension(nxv,nyp1,nblok) :: g
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dcomp
         module procedure idcomp2
      end interface
!
      interface pmove
!        module procedure ipmove2
!        module procedure ipdmove2
         module procedure iwpmove2
         module procedure iwpdmove2
      end interface
!
      interface pmoves
         module procedure iwpmoves2
      end interface
!
      interface pcguard
         module procedure ipcguard2
         module procedure ipdguard2
      end interface
!
      interface pcguardp
         module procedure ipcguard2p
         module procedure ipdguard2p
      end interface
!
      interface pncguardp
         module procedure ipncguard2p
         module procedure ipndguard2p
      end interface
!
      interface paguard
         module procedure ipaguard2
         module procedure ipacguard2
      end interface
!
      interface pamcguard
         module procedure ipamcguard2
      end interface
!
      interface paguardp
         module procedure ipaguard2p
         module procedure ipacguard2p
      end interface
!
      interface pnaguardp
         module procedure ipnaguard2p
         module procedure ipnacguard2p
      end interface
!
      interface pfmove
         module procedure ipfmove2
         module procedure ipfcmove2
         module procedure ipnfmove2
      end interface
!
      interface repart
         module procedure irepartd2
      end interface
!
      interface fnoff
         module procedure ifnoff2
      end interface
!
      interface dblsin
         module procedure ipdblsin2c
         module procedure ipdblsin2d
      end interface
!
      interface dblcos
         module procedure ipdblcos2c
         module procedure ipdblcos2d
      end interface
!
      interface hafdbl
         module procedure iphafdbl2c
         module procedure iphafdbl2d
      end interface
!
      interface zdbl
         module procedure izpdbl2c
         module procedure izpdbl2d
      end interface
!
      interface plsum
         module procedure ipsum2
      end interface
!
      interface plbcast
         module procedure ipbcast2
      end interface
!
      interface writebf
         module procedure ipwrite2
         module procedure ipcwrite2
         module procedure ipvcwrite2
      end interface
!
      interface readbf
         module procedure ipread2
         module procedure ipcread2
         module procedure ipvcread2
      end interface
!
      interface wrdata
         module procedure ipwrdata2
         module procedure ipwrrdata2
         module procedure ipwrcdata2
      end interface
!
      interface rddata
         module procedure iprddata2
         module procedure iprdrdata2
         module procedure iprdcdata2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!
         subroutine idcomp2(edges,nyp,noff,ny,kstrt,nvp,inorder)
! find uniform 1d partition boundaries in 2d code
         implicit none
         integer :: ny, kstrt, nvp
         integer, optional :: inorder
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: nyp, noff
! local data
         integer :: idps, nblok, order
         idps = size(edges,1); nblok = size(edges,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call DCOMP2L(edges,nyp,noff,ny,kstrt,nvp,idps,nblok)
         else
            call DCOMP2(edges,nyp,noff,ny,kstrt,nvp,idps,nblok)
         endif
         end subroutine idcomp2
!
         subroutine ipcguard2(f,kstrt,nvp,kyp,inorder)
! copy guard cells in y for uniform, periodic 2d vector data
         implicit none
         integer :: kstrt, nvp, kyp
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kblok, order
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); kblok = size(f,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
         else
            call PCGUARD2(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
         endif
         end subroutine ipcguard2
!
         subroutine ipdguard2(f,kstrt,nvp,kyp,inorder)
! copy guard cells in y for uniform, periodic 2d scalar data
         implicit none
         integer :: kstrt, nvp, kyp
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kblok, order
         nxv = size(f,1); nypmx = size(f,2); kblok = size(f,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
         else
            call PCGUARD2(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
         endif
         end subroutine ipdguard2
!
         subroutine ipcguard2p(f,kstrt,nvp,nx,kyp,ipbc,inorder)
! copy guard cells in y for uniform 2d vector data
         implicit none
         integer :: kstrt, nvp, nx, kyp, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kblok, order
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); kblok = size(f,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            else
               call PCGUARD2(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               call PLCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            else
               nxv = size(f,2)
               if (size(f,1)==2) then
                  call PLCGUARD2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
               else if (size(f,1)==3) then
                  call PLBGUARD2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
               endif
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call PCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            else
               call PCGUARD2(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            endif
         endif
         end subroutine ipcguard2p
!
         subroutine ipdguard2p(f,kstrt,nvp,nx,kyp,ipbc,inorder)
! copy guard cells in y for uniform 2d scalar data
         implicit none
         integer :: kstrt, nvp, nx, kyp, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kblok, order
         nxv = size(f,1); nypmx = size(f,2); kblok = size(f,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            else
               call PCGUARD2(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               call PLCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            else
               call PLDGUARD2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call PCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            else
               call PCGUARD2(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
            endif
         endif
         end subroutine ipdguard2p
!
         subroutine ipncguard2p(f,nyp,kstrt,nvp,nx,mter,ipbc,inorder)
! copy guard cells in y for non-uniform 2d vector data
         implicit none
         integer :: kstrt, nvp, nx, mter, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxv, nypmx, nblok, order
         real, dimension(size(f,1)*size(f,2),size(f,4)) :: scs
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); nblok = size(f,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
            else
               call PNCGUARD2(f,scs,nyp,kstrt,nvp,nxv,nypmx,nblok,mter)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               call PNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
            else
               nxv = size(f,2)
               if (size(f,1)==2) then
                  call PNLCGUARD2(f,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok&
     &,mter)
               else if (size(f,1)==3) then
                  call PNLBGUARD2(f,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok&
     &,mter)
               endif
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call PNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
            else
               call PNCGUARD2(f,scs,nyp,kstrt,nvp,nxv,nypmx,nblok,mter)
            endif
         endif
         end subroutine ipncguard2p
!
         subroutine ipndguard2p(f,nyp,kstrt,nvp,nx,mter,ipbc,inorder)
! copy guard cells in y for non-uniform 2d scalar data
         implicit none
         integer :: kstrt, nvp, nx, mter, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxv, nypmx, nblok, order
         real, dimension(size(f,1),size(f,3)) :: scs
         nxv = size(f,1); nypmx = size(f,2); nblok = size(f,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
            else
               call PNCGUARD2(f,scs,nyp,kstrt,nvp,nxv,nypmx,nblok,mter)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               call PNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
            else
               call PNLDGUARD2(f,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mt&
     &er)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call PNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
            else
               call PNCGUARD2(f,scs,nyp,kstrt,nvp,nxv,nypmx,nblok,mter)
            endif
         endif
         end subroutine ipndguard2p
!
         subroutine ipacguard2(f,kstrt,nvp,nx,kyp,ngds)
! add guard cells in y for uniform, periodic 2d vector data
         implicit none
         integer :: kstrt, nvp, nx, kyp, ngds
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kblok
         real, dimension(size(f,1)*size(f,2),ngds,size(f,4)) :: scr
         nxv = size(f,2); nypmx = size(f,3); kblok = size(f,4)
         if (ngds.eq.1) then
            if (size(f,1)==2) then
               call PACGUARD22L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
            else if (size(f,1)==3) then
               call PACGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
            endif
         else
            if (size(f,1)==2) then
               call PACGUARD22(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ng&
     &ds)
            else if (size(f,1)==3) then
               call PACGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngd&
     &s)
            endif
         endif
         end subroutine ipacguard2
!
         subroutine ipaguard2(f,kstrt,nvp,nx,kyp,ngds)
! add guard cells in y for uniform, periodic 2d scalar data
         implicit none
         integer :: kstrt, nvp, nx, kyp, ngds
         real, dimension(:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kblok
         real, dimension(size(f,1),ngds,size(f,3)) :: scr
         nxv = size(f,1); nypmx = size(f,2); kblok = size(f,3)
         if (ngds.eq.1) then
            call PAGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
         else
            call PAGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds)
         endif
         end subroutine ipaguard2
!
         subroutine ipamcguard2(f,kstrt,nvp,nx,kyp,ngds)
! add guard cells in y for uniform, periodic 2d tensor data
         implicit none
         integer :: kstrt, nvp, nx, kyp, ngds
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kblok, ndim
         real, dimension(size(f,1)*size(f,2),ngds,size(f,4)) :: scr
         nxv = size(f,2); nypmx = size(f,3); kblok = size(f,4)
         ndim = size(f,1)
         if (ngds.eq.1) then
            call PAMCGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ndim&
     &)
         else
            call PAMCGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds,&
     &ndim)
         endif
         end subroutine ipamcguard2
!
         subroutine ipacguard2p(f,kstrt,nvp,nx,kyp,ngds,ipbc)
! add guard cells in y for uniform 2d vector data
         implicit none
         integer :: kstrt, nvp, nx, kyp, ngds, ipbc
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kblok
         real, dimension(size(f,1),size(f,2),ngds,size(f,4)) :: scr
         nxv = size(f,2); nypmx = size(f,3); kblok = size(f,4)
         if (ipbc==1) then
            if (ngds.eq.1) then
               if (size(f,1)==2) then
                  call PACGUARD22L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblo&
     &k)
               else if (size(f,1)==3) then
                  call PACGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok&
     &)
               endif
            else
               if (size(f,1)==2) then
                  call PACGUARD22(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok&
     &,ngds)
               else if (size(f,1)==3) then
                  call PACGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,&
     &ngds)
               endif
            endif
         else if (ipbc==2) then
            if (ngds.eq.1) then
               if (size(f,1)==2) then
                  call PLACGUARD22L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kbl&
     &ok)
               else if (size(f,1)==3) then
                  call PLACGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblo&
     &k)
               endif
            else
               if (size(f,1)==2) then
                  call PLACGUARD22(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblo&
     &k,ngds)
                  call PLACGUARDS22(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
               else if (size(f,1)==3) then
                  call PLACGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok&
     &,ngds)
                  call PLACGUARDS2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
               endif
            endif
         else if (ipbc==3) then
            if (ngds.eq.1) then
               if (size(f,1)==2) then
                  call PACGUARD22L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblo&
     &k)
               else if (size(f,1)==3) then
                  call PACGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok&
     &)
               endif
            else
               if (size(f,1)==2) then
                  call PACGUARD22(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok&
     &,ngds)
               else if (size(f,1)==3) then
                  call PACGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,&
     &ngds)
               endif
            endif
         endif
         end subroutine ipacguard2p
!
         subroutine ipaguard2p(f,kstrt,nvp,nx,kyp,ngds,ipbc)
! add guard cells in y for uniform 2d scalar data
         implicit none
         integer :: kstrt, nvp, nx, kyp, ngds, ipbc
         real, dimension(:,:,:), pointer :: f
! local data
         integer :: nxv, nypmx, kblok
         real, dimension(size(f,1),ngds,size(f,3)) :: scr
         nxv = size(f,1); nypmx = size(f,2); kblok = size(f,3)
         if (ipbc==1) then
            if (ngds.eq.1) then
               call PAGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
            else
               call PAGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds&
     &)
            endif
         else if (ipbc==2) then
            if (ngds.eq.1) then
               call PLAGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
            else
               call PLAGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngd&
     &s)
               call PLAGUARDS2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
            endif
         else if (ipbc==3) then
            if (ngds.eq.1) then
               call PAGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
            else
               call PAGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds&
     &)
            endif
         endif
         end subroutine ipaguard2p
!
         subroutine ipnacguard2p(f,nyp,kstrt,nvp,nx,mter,ngds,ipbc)
! add guard cells in y for non-uniform 2d vector data
         implicit none
         integer :: kstrt, nvp, nx, mter, ngds, ipbc
         real, dimension(:,:,:,:), pointer :: f
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxv, nypmx, nblok
         real, dimension(size(f,1),size(f,2),ngds,size(f,4)) :: scr
         real, dimension(size(f,1),size(f,2),size(f,4)) :: scs
         nxv = size(f,2); nypmx = size(f,3); nblok = size(f,4)
         if (ipbc==1) then
            if (ngds.eq.1) then
               if (size(f,1)==2) then
                  call PNACGUARD22L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nbl&
     &ok)
               else if (size(f,1)==3) then
                  call PNACGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblo&
     &k)
               endif
            else
               if (size(f,1)==2) then
                  call PNACGUARD22(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,&
     &nblok,ngds,mter)
               else if (size(f,1)==3) then
                  call PNACGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,n&
     &blok,ngds,mter)
               endif
            endif
         else if (ipbc==2) then
            if (ngds.eq.1) then
               if (size(f,1)==2) then
                  call PNLACGUARD22L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nb&
     &lok)
               else if (size(f,1)==3) then
                  call PNLACGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nbl&
     &ok)
               endif
            else
               if (size(f,1)==2) then
                  call PNLACGUARD22(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx&
     &,nblok,ngds,mter)
                  call PNLACGUARDS22(f,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,&
     &mter)
               else if (size(f,1)==3) then
                  call PNLACGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,&
     &nblok,ngds,mter)
                  call PNLACGUARDS2(f,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,m&
     &ter)
               endif
            endif
         else if (ipbc==3) then
            if (ngds.eq.1) then
               if (size(f,1)==2) then
                  call PNACGUARD22L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nbl&
     &ok)
               else if (size(f,1)==3) then
                  call PNACGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblo&
     &k)
               endif
            else
               if (size(f,1)==2) then
                  call PNACGUARD22(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,&
     &nblok,ngds,mter)
               else if (size(f,1)==3) then
                  call PNACGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,n&
     &blok,ngds,mter)
               endif
            endif
         endif
         end subroutine ipnacguard2p
!
         subroutine ipnaguard2p(f,nyp,kstrt,nvp,nx,mter,ngds,ipbc)
! add guard cells in y for non-uniform 2d scalar data
         implicit none
         integer :: kstrt, nvp, nx, mter, ngds, ipbc
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxv, nypmx, nblok
         real, dimension(size(f,1),ngds,size(f,3)) :: scr
         real, dimension(size(f,1),size(f,3)) :: scs
         nxv = size(f,1); nypmx = size(f,2); nblok = size(f,3)
         if (ipbc==1) then
            if (ngds.eq.1) then
               call PNAGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
            else
               call PNAGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok&
     &,ngds,mter)
            endif
         else if (ipbc==2) then
            if (ngds.eq.1) then
               call PNLAGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
            else
               call PNLAGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblo&
     &k,ngds,mter)
               call PNLAGUARDS2(f,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mter)
            endif
         else if (ipbc==3) then
            if (ngds.eq.1) then
               call PNAGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
            else
               call PNAGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok&
     &,ngds,mter)
            endif
         endif
         end subroutine ipnaguard2p
!
         subroutine ipdblsin2c(cu,cu2,nx,ny,kstrt,kyp,kyp2,inorder)
! double array in each dimension for 2d vector data
! for dirichlet boundary conditions
         implicit none
         integer :: nx, ny, kstrt, kyp, kyp2
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: cu
         real, dimension(:,:,:,:), pointer :: cu2
! local data
         integer :: nxv, kypd, kblok, k2blok, order
         nxv = size(cu,2);  kypd = size(cu,3); kblok = size(cu,4)
         k2blok = size(cu2,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(cu,1)==3) then
            if (order==LINEAR) then
               call PDBLSIN2B(cu(1,1,1,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,k&
     &yp2,kblok,k2blok)
            else
               call PDBLSIN2B(cu(1,2,2,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,k&
     &yp2,kblok,k2blok)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call PDBLSIN2C(cu(1,1,1,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,k&
     &yp2,kblok,k2blok)
            else
               call PDBLSIN2C(cu(1,2,2,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,k&
     &yp2,kblok,k2blok)
            endif
         endif

         end subroutine ipdblsin2c
!
         subroutine ipdblsin2d(q,q2,nx,ny,kstrt,kyp,kyp2,inorder)
! double array in each dimension for 2d scalar data
! for dirichlet boundary conditions
         implicit none
         integer :: nx, ny, kstrt, kyp, kyp2
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         real, dimension(:,:,:), pointer :: q2
! local data
         integer :: nxv, kypd, kblok, k2blok, order
         nxv = size(q,1);  kypd = size(q,2); kblok = size(q,3)
         k2blok = size(q2,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PDBLSIN2D(q(1,1,1),q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kbl&
     &ok,k2blok)
         else
            call PDBLSIN2D(q(2,2,1),q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kbl&
     &ok,k2blok)
         endif
         end subroutine ipdblsin2d
!
         subroutine ipdblcos2c(cu,cu2,nx,ny,kstrt,kyp,kyp2,inorder)
! double array in each dimension for 2d vector data
! for neumann boundary conditions
         implicit none
         integer :: nx, ny, kstrt, kyp, kyp2
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: cu
         real, dimension(:,:,:,:), pointer :: cu2
! local data
         integer :: nxv, kypd, kblok, k2blok, order
         nxv = size(cu,2);  kypd = size(cu,3); kblok = size(cu,4)
         k2blok = size(cu2,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(cu,1)==3) then
            if (order==LINEAR) then
               call PDBLCOS2B(cu(1,1,1,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,k&
     &yp2,kblok,k2blok)
            else
               call PDBLCOS2B(cu(1,2,2,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,k&
     &yp2,kblok,k2blok)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call PDBLCOS2C(cu(1,1,1,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,k&
     &yp2,kblok,k2blok)
            else
               call PDBLCOS2C(cu(1,2,2,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,k&
     &yp2,kblok,k2blok)
            endif
         endif

         end subroutine ipdblcos2c
!
         subroutine ipdblcos2d(q,q2,nx,ny,kstrt,kyp,kyp2,inorder)
! double array in each dimension for 2d scalar data
! for neumann boundary conditions
         implicit none
         integer :: nx, ny, kstrt, kyp, kyp2
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         real, dimension(:,:,:), pointer :: q2
! local data
         integer :: nxv, kypd, kblok, k2blok, order
         nxv = size(q,1);  kypd = size(q,2); kblok = size(q,3)
         k2blok = size(q2,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PDBLCOS2D(q(1,1,1),q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kbl&
     &ok,k2blok)
         else
            call PDBLCOS2D(q(2,2,1),q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kbl&
     &ok,k2blok)
         endif
         end subroutine ipdblcos2d
!
         subroutine iphafdbl2c(fxy,fxy2,nx,ny,kstrt,kyp,kyp2,inorder)
! copy from double to normal array in each dimension for 2d vector data
         implicit none
         integer :: nx, ny, kstrt, kyp, kyp2
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:,:), pointer :: fxy2
! local data
         integer :: nxv, kypd, kblok, k2blok, order
         nxv = size(fxy,2);  kypd = size(fxy,3); kblok = size(fxy,4)
         k2blok = size(fxy2,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(fxy,1)==2) then
            if (order==LINEAR) then
               call PHAFDBL2C(fxy(1,1,1,1),fxy2,nx,ny,kstrt,nxv,kyp,kypd&
     &,kyp2,kblok,k2blok)
            else
               call PHAFDBL2C(fxy(1,2,2,1),fxy2,nx,ny,kstrt,nxv,kyp,kypd&
     &,kyp2,kblok,k2blok)
            endif
         else if (size(fxy,1)==3) then
            if (order==LINEAR) then
               call PHAFDBL2B(fxy(1,1,1,1),fxy2,nx,ny,kstrt,nxv,kyp,kypd&
     &,kyp2,kblok,k2blok)
            else
               call PHAFDBL2B(fxy(1,2,2,1),fxy2,nx,ny,kstrt,nxv,kyp,kypd&
     &,kyp2,kblok,k2blok)
            endif
         endif
         end subroutine iphafdbl2c
!
         subroutine iphafdbl2d(q,q2,nx,ny,kstrt,kyp,kyp2,inorder)
! copy from double to normal array in each dimension for 2d scalar data
         implicit none
         integer :: nx, ny, kstrt, kyp, kyp2
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         real, dimension(:,:,:), pointer :: q2
! local data
         integer :: nxv, kypd, kblok, k2blok, order
         nxv = size(q,1);  kypd = size(q,2); kblok = size(q,3)
         k2blok = size(q2,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PHAFDBL2D(q(1,1,1),q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kbl&
     &ok,k2blok)
         else
            call PHAFDBL2D(q(2,2,1),q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kbl&
     &ok,k2blok)
         endif
         end subroutine iphafdbl2d
!
         subroutine izpdbl2c(cu,cu2,nx,ny,kstrt,kyp,kyp2,inorder)
! double array in each dimension for 2d vector data, zeroing copies
! for open boundary conditions
         implicit none
         integer :: nx, ny, kstrt, kyp, kyp2
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: cu
         real, dimension(:,:,:,:), pointer :: cu2
! local data
         integer :: nxv, kypd, kblok, k2blok, order
         nxv = size(cu,2);  kypd = size(cu,3); kblok = size(cu,4)
         k2blok = size(cu2,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(cu,1)==3) then
            if (order==LINEAR) then
               call PZDBL2B(cu(1,1,1,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp&
     &2,kblok,k2blok)
            else
               call PZDBL2B(cu(1,2,2,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp&
     &2,kblok,k2blok)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call PZDBL2C(cu(1,1,1,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp&
     &2,kblok,k2blok)
            else
               call PZDBL2C(cu(1,2,2,1),cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp&
     &2,kblok,k2blok)
            endif
         endif
         end subroutine izpdbl2c
!
         subroutine izpdbl2d(q,q2,nx,ny,kstrt,kyp,kyp2,inorder)
! double array in each dimension for 2d scalar data, zeroing copies
! for open boundary conditions
         implicit none
         integer :: nx, ny, kstrt, kyp, kyp2
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         real, dimension(:,:,:), pointer :: q2
! local data
         integer :: nxv, kypd, kblok, k2blok, order
         nxv = size(q,1);  kypd = size(q,2); kblok = size(q,3)
         k2blok = size(q2,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PZDBL2D(q(1,1,1),q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok&
     &,k2blok)
         else
            call PZDBL2D(q(2,2,1),q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok&
     &,k2blok)
         endif
         end subroutine izpdbl2d
!
         subroutine ipmove2(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vt,i&
     &err)
! particle manager: moves particles to appropriate processor
! non-uniform 1d partition boundaries in 2d code
         implicit none
         integer :: ny, kstrt, nvp, nbmax, vt, ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(7) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
! check if size of buffers has changed
         if (szbuf /= nbmax*size(part,3)) then
            if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call PXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr&
     &,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,maskp,inf&
     &o)
         else
            call PMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr&
     &,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
         endif
         if (info(1) /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         ierr = info(1)
         end subroutine ipmove2
!
         subroutine ipdmove2(part,edges,npp,anpav,pibal,tmove,ny,kstrt,n&
     &vp,nbmax,vt,ierr)
! particle manager: moves particles to appropriate processor
! non-uniform 1d partition boundaries in 2d code
! returns load imbalance
         implicit none
         integer :: ny, kstrt, nvp, nbmax, vt, ierr
         real :: anpav, pibal
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(7) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
! check if size of buffers has changed
         if (szbuf /= nbmax*size(part,3)) then
            if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call PXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr&
     &,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,maskp,inf&
     &o)
         else
            call PMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr&
     &,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
         endif
         if (info(1) /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
! calculate percent imbalance
         anpav = real(info(7))/real(nvp)
         if (anpav > 0.0) then
            pibal = max(real(info(2))-anpav,anpav-real(info(3)))/anpav
         endif
         ierr = info(1)
         end subroutine ipdmove2
!
         subroutine iwpmove2(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vt,&
     &ierr)
! particle manager: moves particles to appropriate processor
! non-uniform 1d partition boundaries in 2d code
         implicit none
         integer :: ny, kstrt, nvp, nbmax, vt, ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(7) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
         th = 0.0
! check if size of buffers has changed
         if (szbuf /= nbmax*size(part,3)) then
            if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call WPXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
     &r,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,maskp&
     &,info)
         else
            call WPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
     &r,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
         endif
         if (info(1) /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
         ierr = info(1)
         end subroutine iwpmove2
!
         subroutine iwpdmove2(part,edges,npp,anpav,pibal,tmove,ny,kstrt,&
     &nvp,nbmax,vt,ierr)
! particle manager: moves particles to appropriate processor
! non-uniform 1d partition boundaries in 2d code
! returns load imbalance
         implicit none
         integer :: ny, kstrt, nvp, nbmax, vt, ierr
         real :: anpav, pibal
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(7) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
         th = 0.0
! check if size of buffers has changed
         if (szbuf /= nbmax*size(part,3)) then
            if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call WPXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
     &r,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,maskp&
     &,info)
         else
            call WPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
     &r,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
         endif
         if (info(1) /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
! calculate percent imbalance
         anpav = real(info(7))/real(nvp)
         if (anpav > 0.0) then
            pibal = max(real(info(2))-anpav,anpav-real(info(3)))/anpav
         endif
         ierr = info(1)
         end subroutine iwpdmove2
!
         subroutine iwpmoves2(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vt&
     &,ierr)
! particle manager: moves particles to appropriate processor
! non-uniform 1d partition boundaries in 2d code
         implicit none
         integer :: ny, kstrt, nvp, nbmax, vt, ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(7) :: info
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
! set maximum number of passes
         info(5) = 2
         th = 0.0
! check if size of buffers has changed
         if (szbuf /= nbmax*size(part,3)) then
            if (szbuf /= 0) deallocate(sbufl,sbufr,rbufl,rbufr,ihole)
! allocate buffers
            allocate(sbufl(idimp,nbmax,size(part,3)))
            allocate(sbufr(idimp,nbmax,size(part,3)))
            allocate(rbufl(idimp,nbmax,size(part,3)))
            allocate(rbufr(idimp,nbmax,size(part,3)))
            allocate(ihole(ntmax,size(part,3)))
            szbuf = nbmax*size(part,3)
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call WPXMOVS2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j&
     &sr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,mask&
     &p,info)
         else
            call WPMOVES2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j&
     &sr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info&
     &)
         endif
         if (info(1) /= 0) then
            call MP_END
            call PPABORT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
         ierr = info(1)
         end subroutine iwpmoves2
!
         subroutine ipfcmove2(f,noff,nyp,isign,tfmove,kyp,kstrt,nvp,idps&
     &,mter,ierr,inorder)
! partition manager: moves 2d vector data between uniform
! and non-uniform 1d partitions
         implicit none
         integer :: isign, kyp, kstrt, nvp, idps, mter, ierr
         real :: tfmove
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: f
         integer, dimension(:), pointer :: noff, nyp
! local data
         real, dimension(size(f,1),size(f,2),size(f,3),size(f,4)) :: g
         integer, dimension(size(f,4)) :: noffs, nyps, noffd, nypd
         integer, dimension(idps,size(f,4)) :: jsl, jsr
         integer :: nxv, nypmx, nblok, order
         real :: tfm
         double precision :: dtime
         nxv = size(f,1)*size(f,2); nypmx = size(f,3); nblok = size(f,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tfm,dtime,-1)
         if (order==LINEAR) then
            call PFMOVE2(f(1,1,1,1),g,noff,nyp,noffs,nyps,noffd,nypd,jsr&
     &,jsl,isign,kyp,kstrt,nvp,nxv,nypmx,nblok,idps,mter,ierr)
         else
            call PFMOVE2(f(1,1,2,1),g,noff,nyp,noffs,nyps,noffd,nypd,jsr&
     &,jsl,isign,kyp,kstrt,nvp,nxv,nypmx,nblok,idps,mter,ierr)
         endif
! record time
         call wtimer(tfm,dtime)
         tfmove = tfmove + tfm
         end subroutine ipfcmove2
!
         subroutine ipfmove2(f,noff,nyp,isign,tfmove,kyp,kstrt,nvp,idps,&
     &mter,ierr,inorder)
! partition manager: moves 2d scalar data between uniform
! and non-uniform 1d partitions
         implicit none
         integer :: isign, kyp, kstrt, nvp, idps, mter, ierr
         real :: tfmove
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: noff, nyp
! local data
         real, dimension(size(f,1),size(f,2),size(f,3)) :: g
         integer, dimension(size(f,3)) :: noffs, nyps, noffd, nypd
         integer, dimension(idps,size(f,3)) :: jsl, jsr
         integer :: nxv, nypmx, nblok, order
         real :: tfm
         double precision :: dtime
         nxv = size(f,1); nypmx = size(f,2); nblok = size(f,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tfm,dtime,-1)
         if (order==LINEAR) then
            call PFMOVE2(f(1,1,1),g,noff,nyp,noffs,nyps,noffd,nypd,jsr,j&
     &sl,isign,kyp,kstrt,nvp,nxv,nypmx,nblok,idps,mter,ierr)
         else
            call PFMOVE2(f(1,2,1),g,noff,nyp,noffs,nyps,noffd,nypd,jsr,j&
     &sl,isign,kyp,kstrt,nvp,nxv,nypmx,nblok,idps,mter,ierr)
         endif
! record time
         call wtimer(tfm,dtime)
         tfmove = tfmove + tfm
         end subroutine ipfmove2
!
         subroutine ipnfmove2(f,noff,nyp,noffs,nyps,tfmove,kstrt,nvp,idp&
     &s,ierr,inorder)
! partition manager: moves 2d scalar data between two different
! non-uniform 1d partitions
! noffs and nyps are modified by this call
         implicit none
         integer :: kstrt, nvp, idps, ierr
         real :: tfmove
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: noff, nyp, noffs, nyps
! local data
         real, dimension(size(f,1),size(f,2),size(f,3)) :: g
         integer, dimension(size(f,3)) :: noffd, nypd
         integer, dimension(idps,size(f,3)) :: jsl, jsr
         integer :: isign = 0, mter = 0, kyp = 1, nxv, nypmx, nblok
         integer :: order
         real :: tfm
         double precision :: dtime
         nxv = size(f,1); nypmx = size(f,2); nblok = size(f,3)
         noffd = noff; nypd = nyp
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tfm,dtime,-1)
         if (order==LINEAR) then
            call PFMOVE2(f(1,1,1),g,noff,nyp,noffs,nyps,noffd,nypd,jsr,j&
     &sl,isign,kyp,kstrt,nvp,nxv,nypmx,nblok,idps,mter,ierr)
         else
            call PFMOVE2(f(1,2,1),g,noff,nyp,noffs,nyps,noffd,nypd,jsr,j&
     &sl,isign,kyp,kstrt,nvp,nxv,nypmx,nblok,idps,mter,ierr)
         endif
! record time
         call wtimer(tfm,dtime)
         tfmove = tfmove + tfm
         end subroutine ipnfmove2
!
         subroutine irepartd2(edges,npic,noff,nyp,anpav,kstrt,nvp,nypmx,&
     &nterg,ierr,inorder)
! finds new 1d partitions from old partition and particle information
         implicit none
         integer :: kstrt, nvp, nypmx, nterg, ierr
         real :: anpav
         integer, optional :: inorder
         real, dimension(:,:), pointer :: edges
         integer, dimension(:,:), pointer :: npic
         integer, dimension(:), pointer :: noff, nyp
! local data
         real, dimension(size(npic,1),size(edges,2)) :: edg, eds
         real, dimension(size(edges,1),size(edges,2)) :: eg, es
         real, dimension(2*size(edges,1),size(edges,2)) :: et2
         integer :: nypmin, nypmax, nblok, idps, nypm, order
         idps = size(edges,1); nblok = size(edges,2)
         nypm = size(npic,1)
         ierr = 0
         order = QUADRATIC
         if (present(inorder)) order = inorder
         call REPARTD2(edges,edg,eds,eg,es,et2,npic,noff,nyp,anpav,nypmi&
     &n,nypmax,kstrt,nvp,nblok,idps,nypm)
         if (order==LINEAR) then
            nypmax = nypmax + 1
         else
            nypmax = nypmax + 3
         endif
         if ((nypmin.lt.1).or.(nypmax.gt.nypmx)) then
            write (2,*) 'Field size error: nypmin,nypmax=',nypmin,nypmax
            ierr = 1
         endif
         nterg = nypmin - 1
         end subroutine irepartd2
!
         subroutine ifnoff2(edges,noff,nyp,nypmx,nterg,ierr,inorder)
! finds new 1d partitions arrays from edges
         implicit none
         integer :: nypmx, nterg, ierr
         integer, optional :: inorder
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: noff, nyp
! local data
         integer :: nypmin, nypmax, nblok, idps, order
         idps = size(edges,1); nblok = size(edges,2)
         ierr = 0
         order = QUADRATIC
         if (present(inorder)) order = inorder
         call FNOFF2(edges,noff,nyp,nypmin,nypmax,nblok,idps)
         if (order==LINEAR) then
            nypmax = nypmax + 1
         else
            nypmax = nypmax + 3
         endif
         if ((nypmin.lt.1).or.(nypmax.gt.nypmx)) then
            write (2,*) 'Field size error: nypmin,nypmax=',nypmin,nypmax
            ierr = 1
         endif
         nterg = nypmin - 1
         end subroutine ifnoff2
!
         subroutine ipsum2(f)
! perform global sum of 2d real array
         implicit none
         real, dimension(:,:) :: f
! local data
         integer :: nxyp, nblok
         real, dimension(size(f,1),size(f,2)) :: g
         nxyp = size(f,1)*size(f,2); nblok = 1
         call PSUM(f,g,nxyp,nblok)
         end subroutine ipsum2
!
         subroutine ipbcast2(f)
! broadcast 2d real array
         implicit none
         real, dimension(:,:) :: f
! local data
         integer :: nxp
         nxp = size(f)
         call PBCAST(f,nxp)
         end subroutine ipbcast2
!
         subroutine ipwrite2(f,nx,kyp,iunit,nrec,name,order)
! collects a subset of a distributed real 2d scalar array and writes it
! to a direct access binary file, for uniform 1d partitions
         implicit none
         integer :: nx, kyp, iunit, nrec
         real, dimension(:,:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, kypmx, nblok, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1)
            lrec = nx*kyp*lrec
         endif
         nxv = size(f,1); kypmx = size(f,2); nblok = size(f,3)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call PWRITE2(f(1,1,1),nx,kyp,nxv,kypmx,nblok,iunit,nrec,l&
     &rec,name)
            else
               call PWRITE2(f(2,2,1),nx,kyp,nxv,kypmx,nblok,iunit,nrec,l&
     &rec,name)
            endif
         else
            if (inorder==LINEAR) then
               call PWRITE2(f(1,1,1),nx,kyp,nxv,kypmx,nblok,iunit,nrec,l&
     &rec,noname)
            else
               call PWRITE2(f(2,2,1),nx,kyp,nxv,kypmx,nblok,iunit,nrec,l&
     &rec,noname)
            endif
         endif
         end subroutine ipwrite2
!
         subroutine ipread2(f,nx,kyp,iunit,nrec,ierr,name,order)
! reads a subset of a distributed real 2d scalar array from a direct
! access binary file and distributes it, for uniform 1d partitions
         implicit none
         integer :: nx, kyp, iunit, nrec, ierr
         real, dimension(:,:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, kypmx, nblok, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1)
            lrec = nx*kyp*lrec
         endif
         nxv = size(f,1); kypmx = size(f,2); nblok = size(f,3)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call PREAD2(f(1,1,1),nx,kyp,nxv,kypmx,nblok,iunit,nrec,lr&
     &ec,name,ierr)
            else
               call PREAD2(f(2,2,1),nx,kyp,nxv,kypmx,nblok,iunit,nrec,lr&
     &ec,name,ierr)
            endif
         else
            if (inorder==LINEAR) then
               call PREAD2(f(1,1,1),nx,kyp,nxv,kypmx,nblok,iunit,nrec,lr&
     &ec,noname,ierr)
            else
               call PREAD2(f(2,2,1),nx,kyp,nxv,kypmx,nblok,iunit,nrec,lr&
     &ec,noname,ierr)
            endif
         endif
         end subroutine ipread2
!
         subroutine ipcwrite2(f,nx,ny,kxp,iunit,nrec,name)
! collects a subset of a distributed complex 2d scalar array and writes
! it to a direct access binary file, for uniform 1d partitions
         implicit none
         integer :: nx, ny, kxp, iunit, nrec
         complex, dimension(:,:,:), pointer :: f
         character(len=*), optional :: name
! local data
         integer :: nyv, kxpd, jblok, lrec
         character(len=1) :: noname = ' '
         nyv = size(f,1); kxpd = size(f,2); jblok = size(f,3)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1)
            lrec = ny*lrec
         endif
         if (present(name)) then
            call PCWRITE2(f,nx,ny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,nam&
     &e)
         else
            call PCWRITE2(f,nx,ny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,non&
     &ame)
         endif
         end subroutine ipcwrite2
!
         subroutine ipcread2(f,nx,ny,kxp,iunit,nrec,ierr,name)
! reads a subset of a distributed complex 2d scalar array from a direct
! access binary file and distributes it, for uniform 1d partitions
         implicit none
         integer :: nx, ny, kxp, iunit, nrec, ierr
         complex, dimension(:,:,:), pointer :: f
         character(len=*), optional :: name
! local data
         integer :: nyv, kxpd, jblok, lrec
         character(len=1) :: noname = ' '
         nyv = size(f,1); kxpd = size(f,2); jblok = size(f,3)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1)
            lrec = ny*lrec
         endif
         if (present(name)) then
            call PCREAD2(f,nx,ny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,name&
     &,ierr)
         else
            call PCREAD2(f,nx,ny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,nona&
     &me,ierr)
         endif
         end subroutine ipcread2
!
         subroutine ipvcwrite2(f,nx,ny,kxp,iunit,nrec,name)
! collects a subset of a distributed complex 2d vector array and writes
! it to a direct access binary file, for uniform 1d partitions
         implicit none
         integer :: nx, ny, kxp, iunit, nrec
         complex, dimension(:,:,:,:), pointer :: f
         character(len=*), optional :: name
! local data
         integer :: nny, nyv, kxpd, jblok, lrec
         character(len=1) :: noname = ' '
         nny = size(f,1)*ny
         nyv = size(f,1)*size(f,2)
         kxpd = size(f,3); jblok = size(f,4)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1,1)
            lrec = nny*lrec
         endif
         if (present(name)) then
            call PCWRITE2(f,nx,nny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,na&
     &me)
         else
            call PCWRITE2(f,nx,nny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,no&
     &name)
         endif
         end subroutine ipvcwrite2
!
         subroutine ipvcread2(f,nx,ny,kxp,iunit,nrec,ierr,name)
! reads a subset of a distributed complex 2d vector array from a direct
! access binary file and distributes it, for uniform 1d partitions
         implicit none
         integer :: nx, ny, kxp, iunit, nrec, ierr
         complex, dimension(:,:,:,:), pointer :: f
         character(len=*), optional :: name
! local data
         integer :: nny, nyv, kxpd, jblok, lrec
         character(len=1) :: noname = ' '
         nny = size(f,1)*ny
         nyv = size(f,1)*size(f,2)
         kxpd = size(f,3); jblok = size(f,4)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1,1)
            lrec = nny*lrec
         endif
         if (present(name)) then
            call PCREAD2(f,nx,nny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,nam&
     &e,ierr)
         else
            call PCREAD2(f,nx,nny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,non&
     &ame,ierr)
         endif
         end subroutine ipvcread2
!
         subroutine ipwrdata2(f,nvp,iunit)
! collects distributed real 2d scalar data and writes it to a
! fortran unformatted sequential file, for uniform 1d partitions
         implicit none
         integer :: nvp, iunit
         real, dimension(:,:,:), pointer :: f
! local data
         integer :: nxv, kyp, nblok
         nxv = size(f,1); kyp = size(f,2)
         nblok = size(f,3)
         call PWRDATA(f,nvp,nxv*kyp,nblok,iunit)
         end subroutine ipwrdata2
!
         subroutine iprddata2(f,nvp,iunit,ierror)
! reads real 2d scalar data from a fortran unformatted sequential file
! and distributes it, for uniform 1d partitions
         implicit none
         integer :: nvp, iunit, ierror
         real, dimension(:,:,:), pointer :: f
! local data
         integer :: nxv, kyp, nblok
         nxv = size(f,1); kyp = size(f,2)
         nblok = size(f,3)
         call PRDDATA(f,nvp,nxv*kyp,nblok,iunit,ierror)
         end subroutine iprddata2
!
         subroutine ipwrrdata2(f,nvp,iunit)
! collects distributed real 2d vector data and writes it to a
! fortran unformatted sequential file, for uniform 1d partitions
         implicit none
         integer :: nvp, iunit
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, kyp, nblok
         nxv = size(f,1)*size(f,2); kyp = size(f,3)
         nblok = size(f,4)
         call PWRDATA(f,nvp,nxv*kyp,nblok,iunit)
         end subroutine ipwrrdata2
!
         subroutine iprdrdata2(f,nvp,iunit,ierror)
! reads real 2d vector data from a fortran unformatted sequential file
! and distributes it, for uniform 1d partitions
         implicit none
         integer :: nvp, iunit, ierror
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, kyp, nblok
         nxv = size(f,1)*size(f,2); kyp = size(f,3)
         nblok = size(f,4)
         call PRDDATA(f,nvp,nxv*kyp,nblok,iunit,ierror)
         end subroutine iprdrdata2
!
         subroutine ipwrcdata2(f,nvp,iunit)
! collects distributed complex 2d vector data and writes it to a
! fortran unformatted sequential file, for uniform 1d partitions
         implicit none
         integer :: nvp, iunit
         complex, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, kyp, nblok
         nxv = size(f,1)*size(f,2); kyp = size(f,3)
         nblok = size(f,4)
         call PWRDATA(f,nvp,2*nxv*kyp,nblok,iunit)
         end subroutine ipwrcdata2
!
         subroutine iprdcdata2(f,nvp,iunit,ierror)
! reads complex 2d vector data from a fortran unformatted sequential file
! and distributes it, for uniform 1d partitions
         implicit none
         integer :: nvp, iunit, ierror
         complex, dimension(:,:,:,:), pointer :: f
! local data
         integer :: nxv, kyp, nblok
         nxv = size(f,1)*size(f,2); kyp = size(f,3)
         nblok = size(f,4)
         call PRDDATA(f,nvp,2*nxv*kyp,nblok,iunit,ierror)
         end subroutine iprdcdata2
!
      end module p2d
