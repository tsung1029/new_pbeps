!-----------------------------------------------------------------------
!
      module pinit2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pinit2lib.f
! pinit2mod.f defines namelists containing input and output variables.
!             contains interface procedures to initialize particle
!             co-ordinates:
!             defines module pinit2d
! distr => ipistr2 initializes x, y and vx, vy co-ordinates for 2d code,
!          with uniform density and maxwellian velocity with drift.
!          calls PISTR2
! distr => ipistrh2 initializes x, y and vx, vy, vz co-ordinates for
!          magnetized 2-1/2d codes, with uniform density and maxwellian
!          velocity with drift.
!          calls PISTR2H
! distr => ipbpistr2 calculates guiding centers for magnetized
!          2-1/2d codes.
!          calls PGBDISTR2L, or PGBZDISTR2L
! distr => iprbpistr2 calculates guiding centers for relativistic,
!          magnetized 2-1/2d codes.
!          calls PGRBDISTR2L, or PGRBZDISTR2L
! pvdistr => ipvistr2 initializes x, y and vx, vy co-ordinates for 2d
!            code, with uniform density and maxwellian velocity with
!            drift, using parallel random number generator.
!            calls PVISTR2
! pvdistr => ipvistrh2 initializes x, y and vx, vy, vz co-ordinates for
!            magnetized 2-1/2d codes, with uniform density and
!            maxwellian velocity with drift, using parallel random
!            number generator.
!            calls PVISTR2H
! fdistr => ipfdistr2 initializes x, y co-ordinates for 2d code, with
!           various density profiles.
!           calls PFDISTR2, and PRDISTR2
! vdistr => ipvdistr2 initializes vx, vy co-ordinates for 2d code, with
!           maxwellian velocity distribution with drift.
!           calls PVDISTR2
! vdistr => ipvdistrg2 initializes vx, vy co-ordinates for 2-1/2d code,
!           with maxwellian velocity distribution with drift.
!           calls PVDISTR2, or PVDISTR2H
! vfdistr => ipvfdistr2 initializes x, y co-ordinates for 2d code, with
!            various density profiles, using parallel random number
!            generator.
!            calls PFDISTR2, and PVRDISTR2
! vvdistr => ipvvdistr2 initializes vx, vy co-ordinates for 2d code,
!            with maxwellian velocity distribution with drift, using
!            parallel random number generator.
!            calls PVVISTR2
! vvdistr => ipvvdistrg2 initializes vx, vy co-ordinates for 2-1/2d
!            code, with maxwellian velocity distribution with drift,
!            using parallel random number generator.
!            calls PVVISTR2, or PVVISTR2H
! fedges => ifedges2 finds new 1d partitions from initial analytic
!           distribution function.
!           calls FEDGES2
! sendnml =>  sendnml2 packs 2d namelist variables into a double
!             precision buffer and broadcasts them to other nodes.
!             calls PBDCAST
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: september 15, 2011
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR,&
     & PERIODIC_2D, DIRICHLET_2D, DIRICHLET_PERIODIC_2D, VACUUM_2D,     &
     &VACUUM_3D, NEUMANN_2D
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PERIODIC_2D, DIRICHLET_2D, DIRICHLET_PERIODIC_2D
      public :: VACUUM_2D, VACUUM_3D, NEUMANN_2D
      public :: idrun, idrun0, idcode, indx, indy, npx, npy, npxb, npyb
      public :: inorder, popt, dopt, djopt, nustrt, ntr
      public :: ntw, ntp, ntd, nta, ntv, nts, ntm, ntj, nte
      public :: ndw, ndp, ndd, nda, ndv, nds, ndm, ndj, nde
      public :: tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0 
      public :: vdx, vdy, vdz, vtdx, vtdy, vtdz
      public :: psolve, relativity, omx, omy, omz, ci, ax, ay
      public :: ndim, ndc, movion, npxi, npyi, npxbi, npybi
      public :: qmi, rmass, rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0
      public :: vdxi, vdyi, vdzi, rtempdxi, rtempdyi, rtempdzi, v0, w0
      public :: sortime, sortimi, nplot, idpal, ndstyle, sntasks
      public :: itpon, ionoff, nsrand, ndprof, nsrandi, ndprofi
      public :: ampdx, scaledx, shiftdx, ampdy, scaledy, shiftdy
      public :: ampdxi, scaledxi, shiftdxi, ampdyi, scaledyi, shiftdyi
      public :: modesxd, modesyd, modesxp, modesyp, modesxa, modesya
      public :: modesxj, modesyj, modesxe, modesye
      public :: imbalance, monitor
      public :: pinput2, sendnml, distr, pvdistr, ldistr, fdistr, vdistr
      public :: vfdistr, vvdistr, fedges
      public :: t0, ceng, indian, rlprec, inprec, pden2d, ndrec, fdname
      public :: ppot2d, nprec, fpname, pvpot2d, narec, faname, pvcur2d
      public :: njrec, fjname, pem2d, nerec, fename
!
! Namelist Input
      save
! idrun/idrun0 = run identifier for current/old run
! idcode = code identifier
      integer :: idrun = 0, idrun0 = 0, idcode = 0
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! npx/npy = initial number of particles distributed in x/y direction
!     integer :: indx =   5, indy =   6, npx =      96, npy =     192
      integer :: indx =   6, indy =   7, npx =     384, npy =     768
!     integer :: indx =   7, indy =   8, npx =    1280, npy =    2560
! npxb/npyb = initial number of particles in beam in x/y direction
!     integer :: npxb =   0, npyb =   0
      integer :: npxb =  32, npyb =  64
!     integer :: npxb = 128, npyb = 256
!     integer :: npxb = 384, npyb = 768
! inorder = interpolation order
! popt = particle optimization scheme
! dopt = charge deposit optimization scheme
! djopt = current deposit optimization scheme
      integer :: inorder = LINEAR, popt = STANDARD, dopt = LOOKAHEAD
      integer :: djopt = STANDARD
! nustrt = (0,1,2) = this is an (old start,new start,restart)
! ntr = number of time steps between restart routine
      integer :: nustrt = 1, ntr = 0
! ntw, ndw = number of time steps between energy diagnostic
! ntp, ndp = number of time steps between potential diagnostic
! ntd, ndd = number of time steps between ion density diagnostic
! nta, nda = number of time steps between vector potential diagnostic
! ntv, ndv = number of time steps between velocity-space diagnostic
! nts, nds = number of time steps between phase space diagnostic
      integer :: ntw = 1, ntp = 0, ntd = 0, nta = 0, ntv = 0, nts = 0
      integer :: ndw = 1, ndp = 0, ndd = 0, nda = 0, ndv = 0, nds = 0
! ntm, ndm = number of time steps between momentum diagnostic
! ntj, ndj = number of time steps between ion current diagnostic
! nte, nde = number of time steps between electromagnetic diagnostic
      integer :: ntm = 0, ntj = 0, nte = 0
      integer :: ndm = 0, ndj = 0, nde = 0
! tend = time at end of simulation, in units of plasma frequency
! dt = time interval between successive calculations
      real :: tend =  65.000, dt = 0.2000000e+00
! qme = charge on electron, in units of e
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
      real :: vdx = 0.0, vdy = 0.0, vdz = 0.0
! vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
      real :: vtdx = 1.0, vtdy = 1.0, vtdz = 1.0
! psolve = type of poisson solver = (1,2,3)
      integer :: psolve = PERIODIC_2D
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: relativity = 0
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
! ci = reciprical of velocity of light
      real :: ci = 1.0
! ax/ay = half-width of particle in x/y direction
!     real :: ax = .816497, ay = .816497
!     real :: ax = .866025, ay = .866025
      real :: ax = .912871, ay = .912871
! ndim = number of velocity dimensions = 2 or 3
      integer :: ndim = 3
! ndc = number of corrections in darwin iteration
      integer :: ndc = 2
! movion = (0,1) = (no,yes) move the ions
! npxi/npyi = initial number of ions distributed in x/y/z direction
      integer :: movion = 0, npxi =  384, npyi =  768
! npxbi/npybi = initial number of ions in beam in x/y/z direction
      integer :: npxbi =   32, npybi =   64
! qmi = charge on ion, in units of 3
! rmass = ion/electron mass ratio
      real :: qmi = 1.0, rmass = 100.0
! rtempxi/rtempyi/rtempzi = electron/ion temperature ratio of background
! ions in x/y/z direction
      real :: rtempxi = 1.0, rtempyi = 1.0, rtempzi = 1.0
! vxi0/vyi0/vzi0 = drift velocity of ions in x/y/z direction
      real :: vxi0 = 0.0, vyi0 = 0.0, vzi0 = 0.0
! vdxi/vdyi/vdzi = drift velocity of beam ions in x/y/z direction
      real :: vdxi = 0.0, vdyi = 0.0, vdzi = 0.0
! rtempdxi/rtempdyi/rtempdzi = electron/ion temperature ratio of beam
! ions in x/y/z direction
      real :: rtempdxi = 1.0, rtempdyi = 1.0, rtempdzi = 1.0
! v0 = external pump strength, in units vos/vthermal
! w0 = external pump frequency, in units of wpe
      real :: v0 = 0.0, w0 = 0.0
! sortime = number of time steps between electron sorting
! sortimi = number of time steps between ion sorting
      integer :: sortime = 50, sortimi = 250
! nplot = maximum number of plots per page
! idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
! ndstyle = (1,2,3) = display (color map,contour plot,both)
      integer :: nplot = 4, idpal = 1, ndstyle = 1
! sntasks = (-1,n) = set maximum number of tasks (-1 = number of cpus-1)
      integer :: sntasks = -1
! itpon = time when external pump is turned on (-1=never)
! ionoff = time when ions are frozen and their charge saved (-1=never)
      integer :: itpon = -1, ionoff = -1
! nsrand = (0,1) = (no,yes) randomize spatially positions locally
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4)
! nsrandi = (0,1) = (no,yes) randomize spatially ion positions locally
! ndprofi = ion profile (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4)
      integer :: nsrand = 0, ndprof = 0, nsrandi = 0, ndprofi = 0
! ampdx/ampdx = amplitude of density compared to uniform in x/y
! scaledx/scaledx = scale length for spatial coordinate in x/y
! shiftdx/shiftdx = shift of spatial coordinate in x/y
      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0
      real :: ampdy = 0.0, scaledy = 0.0, shiftdy = 0.0
! ampdxi/ampdxi = amplitude of ion density compared to uniform in x/y
! scaledxi/scaledxi = scale length for spatial ion coordinate in x/y
! shiftdxi/shiftdxi = shift of spatial ion coordinate in x/y
      real :: ampdxi = 0.0, scaledxi = 0.0, shiftdxi = 0.0
      real :: ampdyi = 0.0, scaledyi = 0.0, shiftdyi = 0.0
! modesxd/modesyd = number of modes in x/y to keep for ion density 
! diagnostic
      integer :: modesxd = 11, modesyd = 11
! modesxp/modesyp = number of modes in x/y to keep for potential 
! diagnostic
      integer :: modesxp = 11, modesyp = 11
! modesxa/modesya = number of modes in x/y to keep for vector potential
! diagnostic
      integer :: modesxa = 11, modesya = 11
! modesxj/modesyj = number of modes in x/y to keep for ion current
! diagnostic
      integer :: modesxj = 11, modesyj = 11
! modesxe/modesye = number of modes in x/y to keep for electromagnetic
! diagnostic
      integer :: modesxe = 11, modesye = 11
!
! imbalance = load imbalance fraction repartition trigger
! (< 0.  to suppress repartition)
      real :: imbalance = .08
! monitor = (0,1,2) = (disable,normal,extended) error processing
      integer :: monitor = 0
! define namelist
      namelist /pinput2/ idrun, idrun0, idcode, indx, indy, npx, npy,   &
     &npxb, npyb, inorder, popt, dopt, djopt, nustrt, ntr, ntw, ntp,    &
     &ntd, nta, ntv, nts, ntm, ntj, nte, ndw, ndp, ndd, nda, ndv, nds,  &
     &ndm, ndj, nde, tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0, vdx,  &
     &vdy, vdz, vtdx, vtdy, vtdz, psolve, relativity, omx, omy, omz, ci,&
     &ax, ay, ndim, ndc, movion, npxi, npyi, npxbi, npybi, qmi, rmass,  &
     &rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0, vdxi, vdyi, vdzi,    &
     &rtempdxi, rtempdyi, rtempdzi, v0, w0, sortime, sortimi, nplot,    &
     &idpal, ndstyle, sntasks, itpon, ionoff, nsrand, ndprof, ampdx,    &
     &scaledx, shiftdx, ampdy, scaledy, shiftdy, nsrandi, ndprofi,      &
     &ampdxi, scaledxi, shiftdxi, ampdyi, scaledyi, shiftdyi, modesxd,  &
     &modesyd, modesxp, modesyp, modesxa, modesya, modesxj, modesyj,    &
     &modesxe, modesye, imbalance, monitor
!
! t0 = initial time value
! ceng = energy normalization
      real :: t0 = 0.0, ceng = 1.0
! indian = (0,1) = architecture is (little-endian,big-endian)
! rlprec = (0,1) = default reals are (normal,double-precision)
! inprec = (0,1) = default integers are (normal,double-precision)
      integer :: indian = 1, rlprec = 1, inprec = 0
!
! Namelist output for ion density diagnostic
! ndrec = current record number for ion density writes
      integer :: ndrec = 0
! fdname = file name for potential diagnostic
      character(len=32) :: fdname
! define namelist
      namelist /pden2d/ idrun, indx, indy, ntd, modesxd, modesyd, psolve&
     &, ndrec, t0, tend, dt, ceng, indian, rlprec, inprec, fdname 
!
! Namelist output for potential diagnostic
! nprec = current record number for potential writes
      integer :: nprec = 0
! fpname = file name for potential diagnostic
      character(len=32) :: fpname
! define namelist
      namelist /ppot2d/ idrun, indx, indy, ntp, modesxp, modesyp, psolve&
     &, omx, omy, omz, nprec, t0, tend, dt, ceng, indian, rlprec, inprec&
     &, fpname
!
! Namelist output for vector potential diagnostic
! narec = current record number for vector potential writes
      integer :: narec = 0
! faname = file name for potential diagnostic
      character(len=32) :: faname
! define namelist
      namelist /pvpot2d/ idrun, indx, indy, nta, modesxa, modesya,      &
     &psolve, omx, omy, omz, ci, narec, t0, tend, dt, ceng, indian,     &
     &rlprec, inprec, faname
!
! Namelist output for ion current diagnostic
! njrec = current record number for ion current writes
      integer :: njrec = 0
! fjname = file name for ion current diagnostic
      character(len=32) :: fjname
! define namelist
      namelist /pvcur2d/ idrun, indx, indy, ntj, modesxj, modesyj,      &
     &psolve, ndim, omx, omy, omz, ci, njrec, t0, tend, dt, ceng,       &
     &indian, rlprec, inprec, fjname 
!
! Namelist output for electromagnetic diagnostic
! nerec = current record number for electromagnetic writes
      integer :: nerec = 0
! fename = file name for potential diagnostic
      character(len=32) :: fename
! define namelist
      namelist /pem2d/ idrun, indx, indy, nte, modesxe, modesye, psolve,&
     &omx, omy, omz, ci, nerec, t0, tend, dt, ceng, indian, rlprec,     &
     &inprec, fename
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx&
     &,ny,idimp,npmax,nblok,idps,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, idps, ipbc
         integer :: ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,i&
     &dimp,npmax,nblok,ipbc,vranx,vrany,kstrt,nvp,ndv,nvrp,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, ipbc, kstrt
         integer :: nvp, ndv, nvrp, ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany
         end subroutine
      end interface
      interface
         subroutine PISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,n&
     &px,npy,nx,ny,idimp,npmax,nblok,idps,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, idps, ipbc
         integer :: ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,np&
     &y,nx,ny,idimp,npmax,nblok,ipbc,vranx,vrany,vranz,kstrt,nvp,ndv,nvr&
     &p,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, ipbc, kstrt
         integer :: nvp, ndv, nvrp, ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany, vranz
         end subroutine
      end interface
      interface
         subroutine PLDISTR2(part,nps,anlx,anly,npx,npy,nx,ny,idimp,npma&
     &x,nblok,kstrt,nvp,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, kstrt, nvp
         integer :: ipbc, ierr
         real :: anlx, anly
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: nps
         end subroutine
      end interface
      interface
         subroutine PFDISTR2(part,nps,fnx,argx1,argx2,argx3,fny,argy1,ar&
     &gy2,argy3,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, kstrt, nvp
         integer :: ipbc, ierr
         real :: argx1, argx2, argx3, argy1, argy2, argy3
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: nps
         real, external :: fnx, fny
         end subroutine
      end interface
      interface
         subroutine PRDISTR2(part,nps,fnx,argx1,argx2,argx3,fny,argy1,ar&
     &gy2,argy3,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, kstrt, nvp
         integer :: ipbc
         real :: argx1, argx2, argx3, argy1, argy2, argy3
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: nps
         real, external :: fnx, fny
         end subroutine
      end interface
      interface
         subroutine PVRDISTR2(part,nps,fnx,argx1,argx2,argx3,fny,argy1,a&
     &rgy2,argy3,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt,nvp,n&
     &dv,nvrp,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, kstrt, nvp
         integer :: ndv, nvrp, ipbc, ierr
         real :: argx1, argx2, argx3, argy1, argy2, argy3
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany
         real, external :: fnx, fny
         end subroutine
      end interface
      interface
         subroutine PVDISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,&
     &npmax,nblok,kstrt,nvp,ierr)
         implicit none
         integer :: npx, npy, idimp, npmax, nblok, kstrt, nvp, ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVVISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,&
     &npmax,nblok,vranx,vrany,kstrt,nvp,ndv,nvrp,ierr)
         implicit none
         integer :: npx, npy, idimp, npmax, nblok, kstrt, nvp, ndv, nvrp
         integer :: ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany
         end subroutine
      end interface
      interface
         subroutine PVDISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,n&
     &py,idimp,npmax,nblok,kstrt,nvp,ierr)
         implicit none
         integer :: npx, npy, idimp, npmax, nblok, kstrt, nvp, ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVVISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,n&
     &py,idimp,npmax,nblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ierr)
         implicit none
         integer :: npx, npy, idimp, npmax, nblok, kstrt, nvp, ndv, nvrp
         integer :: ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany, vranz
         end subroutine
      end interface
      interface
         subroutine PBDISTR2L(part,bx,by,bz,npp,noff,qbm,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx
         real :: qbm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok):: bx, by, bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBDISTR2L(part,bxy,npp,noff,qbm,nx,ny,idimp,npmax,n&
     &blok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm
         real, dimension(idimp,npmax,nblok) :: part
         real :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBZDISTR2L(part,bz,npp,noff,qbm,nx,ny,idimp,npmax,n&
     &blok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm
         real, dimension(idimp,npmax,nblok) :: part
         real :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBDISTR2L(part,bxy,npp,noff,qbm,ci,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBZDISTR2L(part,bz,npp,noff,qbm,ci,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine FEDGES2(edges,noff,nyp,fny,arg1,arg2,arg3,ny,nypmin,&
     &nypmax,kstrt,nvp,nblok,idps,ipbc)
         implicit none
         integer :: ny, nypmin, nypmax, kstrt, nvp, nblok, idps, ipbc
         real :: arg1, arg2, arg3
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: noff, nyp
         real, external :: fny
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface distr
         module procedure ipistr2
         module procedure ipistrh2
         module procedure ipbpistr2
         module procedure iprbpistr2
      end interface
!
      interface pvdistr
         module procedure ipvistr2
         module procedure ipvistrh2
      end interface
!
      interface ldistr
         module procedure ipldistr2
      end interface
!
      interface fdistr
         module procedure ipfdistr2
      end interface
!
      interface vdistr
         module procedure ipvdistr2
!        module procedure ipvdistrh2
         module procedure ipvdistrg2
      end interface
!
      interface vfdistr
         module procedure ipvfdistr2
      end interface
!
      interface vvdistr
         module procedure ipvvdistr2
!        module procedure ipvvdistrh2
         module procedure ipvvdistrg2
      end interface
!
      interface fedges
         module procedure ifedges2
      end interface
!
      interface sendnml
         module procedure sendnml2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine sendnml2()
! this subroutine packs 2d namelist variables into a double precision
! buffer and broadcasts them to other nodes
         integer, parameter :: lenml = 115
         double precision, dimension(lenml) :: ddata
! pack data
         ddata(1) = idrun; ddata(2) = idrun0; ddata(3) = idcode
         ddata(4) = indx; ddata(5) = indy
         ddata(6) = npx; ddata(7) = npy
         ddata(8) = npxb; ddata(9) = npyb
         ddata(10) = inorder; ddata(11) = popt; ddata(12) = dopt
         ddata(13) = djopt;; ddata(14) = nustrt; ddata(15) = ntr
         ddata(16) = ntw; ddata(17) = ntp; ddata(18) = ntd
         ddata(19) = nta; ddata(20) = ntv; ddata(21) = nts
         ddata(22) = ntm; ddata(23) = ntj; ddata(24) = nte
         ddata(25) = ndw; ddata(26) = ndp; ddata(27) = ndd
         ddata(28) = nda; ddata(29) = ndv; ddata(30) = nds
         ddata(31) = ndm; ddata(32) = ndj; ddata(33) = nde
         ddata(34) = tend; ddata(35) = dt; ddata(36) = qme
         ddata(37) = vtx; ddata(38) = vty; ddata(39) = vtz
         ddata(40) = vx0; ddata(41) = vy0; ddata(42) = vz0
         ddata(43) = vdx; ddata(44) = vdy; ddata(45) = vdz
         ddata(46) = vtdx; ddata(47) = vtdy; ddata(48) = vtdz
         ddata(49) = psolve; ddata(50) = relativity
         ddata(51) = omx; ddata(52) = omy; ddata(53) = omz
         ddata(54) = ci; ddata(55) = ax; ddata(56) = ay
         ddata(57) = ndim; ddata(58) = ndc; ddata(59) = movion
         ddata(60) = npxi; ddata(61) = npyi
         ddata(62) = npxbi; ddata(63) = npybi
         ddata(64) = qmi; ddata(65) = rmass
         ddata(66) = rtempxi; ddata(67) = rtempyi; ddata(68) = rtempzi
         ddata(69) = vxi0; ddata(70) = vyi0; ddata(71) = vzi0
         ddata(72) = vdxi; ddata(73) = vdyi; ddata(74) = vdzi
         ddata(75) = rtempdxi; ddata(76) = rtempdyi
         ddata(77) = rtempdzi; ddata(78) = v0; ddata(79) = w0
         ddata(80) = sortime; ddata(81) = sortimi
         ddata(82) = nplot; ddata(83) = idpal; ddata(84) = ndstyle
         ddata(85) = sntasks; ddata(86) = itpon; ddata(87) = ionoff
         ddata(88) = nsrand; ddata(89) = ndprof
         ddata(90) = ampdx; ddata(91) = scaledx; ddata(92) = shiftdx
         ddata(93) = ampdy; ddata(94) = scaledy; ddata(95) = shiftdy
         ddata(96) = nsrandi; ddata(97) = ndprofi
         ddata(98) = ampdxi; ddata(99) = scaledxi; ddata(100) = shiftdxi
         ddata(101) = ampdyi; ddata(102) = scaledyi
         ddata(103) = shiftdyi
         ddata(104) = modesxd; ddata(105) = modesyd
         ddata(106) = modesxp; ddata(107) = modesyp
         ddata(108) = modesxa; ddata(109) = modesya
         ddata(110) = modesxj; ddata(111) = modesyj
         ddata(112) = modesxe; ddata(113) = modesye
         ddata(114) = imbalance; ddata(115) = monitor
! broadcast data
         call PBDCAST(ddata,lenml)
! unpack data
         idrun = ddata(1); idrun0 = ddata(2); idcode = ddata(3)
         indx = ddata(4); indy = ddata(5)
         npx = ddata(6); npy = ddata(7)
         npxb = ddata(8); npyb = ddata(9)
         inorder = ddata(10); popt = ddata(11); dopt = ddata(12)
         djopt = ddata(13); nustrt = ddata(14); ntr = ddata(15)
         ntw = ddata(16); ntp = ddata(17); ntd = ddata(18)
         nta = ddata(19); ntv = ddata(20); nts = ddata(21)
         ntm = ddata(22); ntj = ddata(23); nte = ddata(24)
         ndw = ddata(25); ndp = ddata(26); ndd = ddata(27)
         nda = ddata(28); ndv = ddata(29); nds = ddata(30)
         ndm = ddata(31); ndj = ddata(32); nde = ddata(33)
         tend = ddata(34); dt = ddata(35); qme = ddata(36)
         vtx = ddata(37); vty = ddata(38); vtz = ddata(39)
         vx0 = ddata(40); vy0 = ddata(41); vz0 = ddata(42)
         vdx = ddata(43); vdy = ddata(44); vdz = ddata(45)
         vtdx = ddata(46); vtdy = ddata(47); vtdz = ddata(48)
         psolve = ddata(49); relativity = ddata(50)
         omx = ddata(51); omy = ddata(52); omz = ddata(53)
         ci = ddata(54); ax = ddata(55); ay = ddata(56)
         ndim = ddata(57); ndc = ddata(58); movion = ddata(59)
         npxi = ddata(60); npyi = ddata(61)
         npxbi = ddata(62); npybi = ddata(63)
         qmi = ddata(64); rmass = ddata(65)
         rtempxi = ddata(66); rtempyi = ddata(67); rtempzi = ddata(68)
         vxi0 = ddata(69); vyi0 = ddata(70); vzi0 = ddata(71)
         vdxi = ddata(72); vdyi = ddata(73); vdzi = ddata(74)
         rtempdxi = ddata(75); rtempdyi = ddata(76)
         rtempdzi = ddata(77); v0 = ddata(78); w0 = ddata(79)
         sortime = ddata(80); sortimi = ddata(81)
         nplot = ddata(82); idpal = ddata(83); ndstyle = ddata(84)
         sntasks = ddata(85); itpon = ddata(86); ionoff = ddata(87)
         nsrand = ddata(88); ndprof = ddata(89)
         ampdx = ddata(90); scaledx = ddata(91); shiftdx = ddata(92)
         ampdy = ddata(93); scaledy = ddata(94); shiftdy = ddata(95)
         nsrandi = ddata(96); ndprofi = ddata(97)
         ampdxi = ddata(98); scaledxi = ddata(99); shiftdxi = ddata(100)
         ampdyi = ddata(101); scaledyi = ddata(102)
         shiftdyi = ddata(103)
         modesxd = ddata(104); modesyd = ddata(105)
         modesxp = ddata(106); modesyp = ddata(107)
         modesxa = ddata(108); modesya = ddata(109)
         modesxj = ddata(110); modesyj = ddata(111)
         modesxe = ddata(112); modesye = ddata(113)
         imbalance = ddata(114); monitor = ddata(115)
         end subroutine sendnml2
!
         subroutine ipistr2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,n&
     &x,ny,ipbc)
! calculates initial particle co-ordinates and velocities in 2d
! with uniform density and maxwellian velocity with drift
         implicit none
         integer :: npx, npy, nx, ny, ipbc
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, idps, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); idps = size(edges,1)
         call PISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,id&
     &imp,npmax,nblok,idps,ipbc,ierr)
         end subroutine ipistr2
!
         subroutine ipistrh2(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
     &npx,npy,nx,ny,ipbc)
! calculates initial particle co-ordinates and velocities in 2-1/2d
! with uniform density and maxwellian velocity with drift
         implicit none
         integer :: npx, npy, nx, ny, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, idps, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); idps = size(edges,1)
         call PISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy&
     &,nx,ny,idimp,npmax,nblok,idps,ipbc,ierr)
         end subroutine ipistrh2
!
         subroutine ipvistr2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,&
     &ipbc,kstrt,nvp)
! calculates initial particle co-ordinates and velocities in 2d
! with uniform density and maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, nx, ny, ipbc, kstrt, nvp
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         call PVISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,idimp,n&
     &pmax,nblok,ipbc,vranx,vrany,kstrt,nvp,ndv,nvrp,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvistr2
!
         subroutine ipvistrh2(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,n&
     &py,nx,ny,ipbc,kstrt,nvp)
! calculates initial particle co-ordinates and velocities in 2-1/2d
! with uniform density and maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, nx, ny, ipbc, kstrt, nvp
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany, vranz
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         call PVISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,nx,n&
     &y,idimp,npmax,nblok,ipbc,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ierr&
     &)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvistrh2
!
         subroutine ipldistr2(part,nps,anlx,anly,npx,npy,nx,ny,kstrt,nvp&
     &,ipbc)
! calculates initial particle co-ordinates in 2d
! with bi-linear density profile
         implicit none
         integer :: npx, npy, nx, ny, kstrt, nvp, ipbc
         real :: anlx, anly
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: nps
! local data
         integer :: idimp, npmax, nblok, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PLDISTR2(part,nps,anlx,anly,npx,npy,nx,ny,idimp,npmax,nblo&
     &k,kstrt,nvp,ipbc,ierr)
         end subroutine ipldistr2
!
         subroutine ipfdistr2(part,nps,ampx,scalex,shiftx,ampy,scaley,sh&
     &ifty,npx,npy,nx,ny,kstrt,nvp,ipbc,ndpro,nsran)
! calculates initial particle co-ordinates in 2d
! with various density profiles
         implicit none
         integer :: npx, npy, nx, ny, kstrt, nvp, ipbc, ndpro, nsran
         real :: ampx, scalex, shiftx, ampy, scaley, shifty
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: nps
! local data
         integer :: idimp, npmax, nblok, ierr
         real :: sxi, syi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         sxi = 0.
         if (scalex /= 0.) sxi = 1.0/scalex
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         zero = 0.0
! uniform density
         if (ndpro==0) then
            call PFDISTR2(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,zero&
     &,zero,zero,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,z&
     &ero,zero,zero,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
! linear density
         else if (ndpro==1) then
            call PFDISTR2(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1,&
     &ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
! sinusoidal density
         else if (ndpro==2) then
            call PFDISTR2(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1,&
     &ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
! gaussian density
         else if (ndpro==3) then
            call PFDISTR2(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1,&
     &ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
! hyperbolic secant squared density
         else if (ndpro==4) then
            call PFDISTR2(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1,&
     &ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
         endif
         end subroutine ipfdistr2
!
         subroutine ipvdistr2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,kstrt&
     &,nvp)
! calculates initial particle velocities in 2d
! with maxwellian velocity with drift
         implicit none
         integer :: npx, npy, kstrt, nvp
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PVDISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,npmax,&
     &nblok,kstrt,nvp,ierr)
         end subroutine ipvdistr2
!
         subroutine ipvdistrh2(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,&
     &npy,kstrt,nvp)
! calculates initial particle velocities in 2-1/2d
! with maxwellian velocity with drift
         implicit none
         integer :: npx, npy, kstrt, nvp
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PVDISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idi&
     &mp,npmax,nblok,kstrt,nvp,ierr)
         end subroutine ipvdistrh2
!
         subroutine ipvdistrg2(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,&
     &npy,kstrt,nvp,ndim)
! calculates initial particle velocities in 2d or 2-1/2d
! with maxwellian velocity with drift
         implicit none
         integer :: npx, npy, kstrt, nvp
         integer, optional :: ndim
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, nd, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); nd = 3
         if (present(ndim)) nd = ndim
         select case(nd)
         case (2)
            call PVDISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,npm&
     &ax,nblok,kstrt,nvp,ierr)
         case (3)
            call PVDISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,&
     &idimp,npmax,nblok,kstrt,nvp,ierr)
         end select
         end subroutine ipvdistrg2
!
         subroutine ipvfdistr2(part,nps,ampx,scalex,shiftx,ampy,scaley,s&
     &hifty,npx,npy,nx,ny,kstrt,nvp,ipbc,ndpro,nsran)
! calculates initial particle co-ordinates in 2d
! with various density profiles
! using parallel random number generator
         implicit none
         integer :: npx, npy, nx, ny, kstrt, nvp, ipbc, ndpro
         integer :: nsran
         real :: ampx, scalex, shiftx, ampy, scaley, shifty
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         real :: sxi, syi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         sxi = 0.
         if (scalex /= 0.) sxi = 1.0/scalex
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         zero = 0.0
! uniform density
         if (ndpro==0) then
            call PFDISTR2(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,zero&
     &,zero,zero,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,&
     &zero,zero,zero,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt,n&
     &vp,ndv,nvrp,ipbc,ierr)
            endif
! linear density
         else if (ndpro==1) then
            call PFDISTR2(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1&
     &,ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt&
     &,nvp,ndv,nvrp,ipbc,ierr)
            endif
! sinusoidal density
         else if (ndpro==2) then
            call PFDISTR2(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1&
     &,ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt&
     &,nvp,ndv,nvrp,ipbc,ierr)
            endif
! gaussian density
         else if (ndpro==3) then
            call PFDISTR2(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1&
     &,ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt&
     &,nvp,ndv,nvrp,ipbc,ierr)
            endif
! hyperbolic secant squared density
         else if (ndpro==4) then
            call PFDISTR2(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1&
     &,ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt&
     &,nvp,ndv,nvrp,ipbc,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvfdistr2
!
         subroutine ipvvdistr2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,kstr&
     &t,nvp)
! calculates initial particle velocities in 2d
! with maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, kstrt, nvp
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         call PVVISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,npmax,&
     &nblok,vranx,vrany,kstrt,nvp,ndv,nvrp,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvvdistr2
!
         subroutine ipvvdistrh2(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx&
     &,npy,kstrt,nvp)
! calculates initial particle velocities in 2-1/2d
! with maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, kstrt, nvp
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany, vranz
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         call PVVISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idi&
     &mp,npmax,nblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvvdistrh2
!
         subroutine ipvvdistrg2(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx&
     &,npy,kstrt,nvp,ndim)
! calculates initial particle velocities in 2d or 2-1/2d
! with maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, kstrt, nvp
         integer, optional :: ndim
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, nd, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany, vranz
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); nd = 3
         nvrp = (ndv - 1)/nvp + 1
         if (present(ndim)) nd = ndim
         select case(nd)
         case (2)
            call PVVISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,npm&
     &ax,nblok,vranx,vrany,kstrt,nvp,ndv,nvrp,ierr)
         case (3)
            call PVVISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,&
     &idimp,npmax,nblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvvdistrg2
!
         subroutine ipbpistr2(part,bxy,npp,noff,qbm,nx,ny,ipbc,inorder)
! reinterprets curent particle positions as positions of guiding centers
! and calculates the actual particle positions for 2d
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, order
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(bxy,2); nypmx = size(bxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               call PGBZDISTR2L(part,bxy(1,1,1,1),npp,noff,qbm,nx,ny,idi&
     &mp,npmax,nblok,nxv,nypmx,ipbc)
            else
               call PGBZDISTR2L(part,bxy(1,2,2,1),npp,noff,qbm,nx,ny,idi&
     &mp,npmax,nblok,nxv,nypmx,ipbc)
            endif
         case (3)
            if (order==LINEAR) then
               call PGBDISTR2L(part,bxy(1,1,1,1),npp,noff,qbm,nx,ny,idim&
     &p,npmax,nblok,nxv,nypmx,ipbc)
            else
               call PGBDISTR2L(part,bxy(1,2,2,1),npp,noff,qbm,nx,ny,idim&
     &p,npmax,nblok,nxv,nypmx,ipbc)
            endif
         end select
         end subroutine ipbpistr2
!
         subroutine iprbpistr2(part,bxy,npp,noff,qbm,ci,nx,ny,ipbc,inord&
     &er)
! reinterprets curent particle positions as positions of guiding centers
! and calculates the actual particle positions for relativistic 2d
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm, ci
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, order
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(bxy,2); nypmx = size(bxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               call PGRBZDISTR2L(part,bxy(1,1,1,1),npp,noff,qbm,ci,nx,ny&
     &,idimp,npmax,nblok,nxv,nypmx,ipbc)
            else
               call PGRBZDISTR2L(part,bxy(1,2,2,1),npp,noff,qbm,ci,nx,ny&
     &,idimp,npmax,nblok,nxv,nypmx,ipbc)
            endif
         case (3)
            if (order==LINEAR) then
               call PGRBDISTR2L(part,bxy(1,1,1,1),npp,noff,qbm,ci,nx,ny,&
     &idimp,npmax,nblok,nxv,nypmx,ipbc)
            else
               call PGRBDISTR2L(part,bxy(1,2,2,1),npp,noff,qbm,ci,nx,ny,&
     &idimp,npmax,nblok,nxv,nypmx,ipbc)
            endif
         end select
         end subroutine iprbpistr2
!
         subroutine ifedges2(edges,noff,nyp,ampy,scaley,shifty,ny,kstrt,&
     &nvp,nypmx,ipbc,ndpro,nterg,ierr,inorder)
! finds new 1d partitions from initial analytic distribution function
         implicit none
         integer :: ny, kstrt, nvp, nypmx, ipbc, ndpro, nterg, ierr
         real :: ampy, scaley, shifty
         integer, optional :: inorder
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: noff, nyp
! local data
         integer :: idps, nblok, nypmin, nypmax, order
         real :: syi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         idps = size(edges,1); nblok = size(edges,2)
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         zero = 0.0
         ierr = 0
         order = QUADRATIC
         if (present(inorder)) order = inorder
! uniform density
         if (ndpro==0) then
            call FEDGES2(edges,noff,nyp,FLDISTR1,zero,zero,zero,ny,nypmi&
     &n,nypmax,kstrt,nvp,nblok,idps,ipbc)
! linear density
         else if (ndpro==1) then
            call FEDGES2(edges,noff,nyp,FLDISTR1,ampy,syi,shifty,ny,nypm&
     &in,nypmax,kstrt,nvp,nblok,idps,ipbc)
! sinusoidal density
         else if (ndpro==2) then
            call FEDGES2(edges,noff,nyp,FSDISTR1,ampy,syi,shifty,ny,nypm&
     &in,nypmax,kstrt,nvp,nblok,idps,ipbc)
! gaussian density
         else if (ndpro==3) then
            call FEDGES2(edges,noff,nyp,FGDISTR1,ampy,syi,shifty,ny,nypm&
     &in,nypmax,kstrt,nvp,nblok,idps,ipbc)
! hyperbolic secant squared density
         else if (ndpro==4) then
            call FEDGES2(edges,noff,nyp,FHDISTR1,ampy,syi,shifty,ny,nypm&
     &in,nypmax,kstrt,nvp,nblok,idps,ipbc)
         endif
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
         end subroutine ifedges2
!
      end module pinit2d
