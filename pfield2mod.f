!-----------------------------------------------------------------------
!
      module pfield2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pfield2lib.f
! pfield2mod.f contains procedures to manage guard cells and solve
!              fields in fourier space:
!              defines module pfield2d
! cguard => ipcguard2x copies guard cells in x for non-uniform,
!           periodic 2d vector data.
!           calls PCGUARD2X, PBGUARD2X, PDGUARD2X, PCGUARD2XL,
!           PBGUARD2XL, or PDGUARD2XL
! cguard => ipdguard2x copies guard cells in x for non-uniform,
!           periodic 2d scalar data.
!           calls PDGUARD2X,  or PDGUARD2XL
! sguard => ipsguard2 initialize field for scalar array with various
!           interpolations.
!           calls PSGUARD2, or PSGUARD2L
! sguard => ipscguard2 initialize field for 2 or 3 component vector
!           array with various interpolations.
!           calls PSCGUARD2, PSCGUARD22, PSCGUARD2L, or PSCGUARD22L
! sguard => ipsfguard2 initialize 2 or 3 component field with scaled
!           vector array with various interpolations.
!           calls PSCFGUARD2, PSCFGUARD22, PSCFGUARD2L, or PSCFGUARD22L
! sguard => ipsmcguard2 initialize 2d tensor field  with various
!           interpolations.
!           calls PSMCGUARD2, PSMCGUARD22, PSMCGUARD2L, or PSMCGUARD22L
! aguard => ipaguard2x adds guard cells in x for non-uniform, periodic
!           2d scalar data, with various interpolations.
!           calls PAGUARD2X,  or PAGUARD2XL
! aguard => ipacguard2x adds guard cells in x for non-uniform, periodic
!           2d vector data, with various interpolations.
!           calls PACGUARD2X, PACGUARD22X, PACGUARD2XL, or PACGUARD22XL
! amcguard => ipamcguard2x adds guard cells in x for for non-uniform,
!             periodic 2d tensor field, with various interpolations.
!             calls PAMCGUARD2X, or PAMCGUARD2XL
! zguard => ipzguard2 zeros out guard cells in periodic 2d scalar field.
!           calls PZGUARD2, or PZGUARD2L
! pois_init => ippois22init initializes tables for field solvers.
!              calls PPOIS22
! pois => ippois2 solves poisson equation for electric force, potential,
!         or smoothing.
!         calls PPOISP2
! pois => ippois22 solves 2d poisson equation for electric force.
!         calls PPOIS22
! spois => ipspois2 smoother for 2d periodic scalar field.
!          calls PPOISP2
! pois3 => ippois23 solves 2-1/2d poisson equation for electric force.
!          calls PPOIS22, or PPOIS23
! cuperp => ipcuperp2 calculates the transverse part of periodic 2-1/2d
!           vector field.
!           calls PCUPERP2, or PCUPERP22
! bpois => jpbpois23 solves 2-1/2d vector poisson equation for magnetic
!          force.
!          calls PBPOISP23, or PBPOISP22
! sbpois => ipsbpois23 smoother for 2-1/2d periodic vector field.
!           calls PBPOISP23, or PBPOISP22
! apois => ipapois23 solves 2-1/2d vector poisson equation for vector
!          potential.
!          calls PBPOISP23, or PBPOISP22
! ibpois => iipbpois23 solves vector poisson equation for magnetic field
!           calls IPBPOISP23
! maxwel => ipmaxwel2 solves maxwell equation for electric and magnetic
!           fields.
!           calls PMAXWEL2
! emfield => ipemfield2 calculates periodic electric and magnetic forces
!            from fields given by maxwell and poisson equations.
!            calls PEMFIELD2
! emfieldr => ipemfieldr2 calculates electric and magnetic forces from
!             fields given by maxwell and poisson equations for
!             sine-cosine transforms.
!             calls PEMFIELDR2
! avpot => ipavpot23 calculates vector potential from magnetic field.
!          calls PAVPOT23
! avrpot => ipavrpot23 calculates radiative vector potential from
!           magnetic field and current.
!           calls PAVRPOT23
! gtmodes => ipgtmodes2 extracts selected fourier components from
!            potential array.
!            calls PGTMODES2
! gtmodes => ipgtvmodes2 extracts selected fourier components from
!            vector potential array.
!            calls PGTVMODES2
! ptmodes => ipptmodes2 places selected fourier components into potential
!            array.
!            calls PPTMODES2
! ptmodes => ipptvmodes2 places selected fourier components into vector
!            potential array.
!            calls PPTVMODES2
! poynt => ippoynt2 calculates the momentum in the electromagnetic field
!          calls PPOYNT2
! poynt => ipdpoynt2 calculates the momentum in the darwin field.
!          calls PDPOYNT2, or PDPOYNT22
! dcuperp => ipdcuperp23 calculate transverse derivative of 2-1/2d
!            current density from momentum flux.
!            calls PDCUPERP23, or PDCUPERP22
! adcuperp => ipadcuperp23 calculate transverse derivative of 2-1/2d
!             current density from momentum flux and acceleration
!             density.
!             calls PADCUPERP23, or PADCUPERP22
! epois_init => ipepois23init initializes tables for darwin field solver
!               calls PEPOISP23
! epois => ipepois23 solves 2-1/2d vector vector poisson equation for
!          transverse electric force.
!          calls PEPOISP23, or PEPOISP22
! iepois => iipepois23 solves 2-1/2d vector vector poisson equation for
!           transverse electric field.
!           calls PEPOISP23, or PEPOISP22
! addqei => ipaddqei2 adds electron and ion densities.
!           calls PADDQEI2
! addqei => ipaddqei2x adds electron and ion densities, and calculates
!           maximum and minimum plasma frequency.
!           calls PADDQEI2X
! baddext => ipbaddext2 adds constant to magnetic field in real space for
!            2-1/2d code.
!            calls PBADDEXT2, or PBADDEXT22
! imoment => ipimoment2 calculates ion momentum from integral of qi*fxy.
!            calls PIMOMENT2, or PIMOMENT22
! addfields => ipaddvrfield2 calculates a = b + c for distributed real
!              vector fields.
!              calls PADDVRFIELD2
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 15, 2009
!
      use globals, only: LINEAR, QUADRATIC
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: PSCGUARD2, PSGUARD2, PSCGUARD2L, PSGUARD2L
      public :: PACGUARD2X, PAGUARD2X, PACGUARD2XL, PAGUARD2XL
      public :: cguard, bguard, sguard, aguard, zguard
      public :: pois_init, pois, spois, pois3, cuperp, bpois, sbpois
      public :: ibpois, maxwel, emfield, emfieldr, apois, avpot, avrpot
      public :: gtmodes, ptmodes, poynt
      public :: amcguard, dcuperp, adcuperp, epois_init, epois, iepois
      public :: addqei, baddext
      public :: imoment, ipdivf2, ipgradf2, ipcurlf2, ipcurlf22
      public :: addfields
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PBGUARD2X(bxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD2(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD22(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSGUARD2(q,nyp,qi0,nx,nxe,nypmx,nblok)
         implicit none
         real :: qi0
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD2X(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD22X(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PAGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PZGUARD2(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PCGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PDGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PBGUARD2XL(bxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD2L(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD22L(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSGUARD2L(q,nyp,qi0,nx,nxe,nypmx,nblok)
         implicit none
         real :: qi0
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD2XL(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD22XL(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PAGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PZGUARD2L(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PPOISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(nyv,kxp,jblok) :: q, fx, fy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,ny&
     &v,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(nyv,kxp,jblok) :: q
         complex, dimension(2,nyv,kxp,jblok) :: fxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,ny&
     &v,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(nyv,kxp,jblok) :: q
         complex, dimension(2,nyv,kxp,jblok) :: fxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, ndim, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: f
         complex, dimension(nyv,kxp,jblok) :: df
         end subroutine
      end interface
      interface
         subroutine PGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, ndim, nyv, kxp, jblok
         complex, dimension(nyv,kxp,jblok) :: df
         complex, dimension(3,nyv,kxp,jblok) :: f
         end subroutine
      end interface
      interface
         subroutine PCURLF2(f,g,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PCURLF22(f,g,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(2,nyv,kxp,jblok) :: f
         complex, dimension(nyv,kxp,jblok) :: g
         end subroutine
      end interface
      interface
         subroutine PCUPERP2(cu,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERP22(cu,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(2,nyv,kxp,jblok) :: cu
         end subroutine
      end interface
      interface
         subroutine PBPOISP23(cu,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,ks&
     &trt,nyv,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(3,nyv,kxp,jblok) :: cu
         complex, dimension(3,nyv,kxp,jblok) :: bxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PBPOISP22(cu,bxy,bz,isign,ffc,ax,ay,affp,ci,wm,nx,ny&
     &,kstrt,nyv,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(2,nyv,kxp,jblok) :: cu
         complex, dimension(2,nyv,kxp,jblok) :: bxy
         complex, dimension(nyv,kxp,jblok) :: bz
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine IPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,jblo&
     &k,nyhd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(3,nyv,kxp,jblok) :: cu
         complex, dimension(3,nyv,kxp,jblok) :: bxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt&
     &,nyv,kxp,jblok,nyhd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(3,nyv,kxp,jblok) :: exy, bxy, cu
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,jblo&
     &k,nyhd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(3,nyv,kxp,jblok) :: fxy, exy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PEMFIELDR2(fxy,exy,ffd,isign,nx,ny,kstrt,nyv,kxp2,j2&
     &blok,nyd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(3,nyv,kxp2+1,j2blok) :: fxy, exy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: bxy, axy
         end subroutine
      end interface
      interface
         subroutine PAVRPOT23(axy,bxy,ffc,affp,ci,nx,ny,kstrt,nyv,kxp,jb&
     &lok,nyhd)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok, nyhd
         real :: affp, ci
         complex, dimension(3,nyv,kxp,jblok) :: axy, bxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PGTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,&
     &kxp,jblok,nt,modesxpd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, kstrt, nyv, kxp, jblok
         integer :: nt, modesxpd, modesyd
         complex, dimension(nyv,kxp,jblok) :: pot
         complex, dimension(nt,modesyd,modesxpd,jblok) :: pott
         end subroutine
      end interface
      interface
         subroutine PPTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,&
     &kxp,jblok,nt,modesxpd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, kstrt, nyv, kxp, jblok
         integer :: nt, modesxpd, modesyd
         complex, dimension(nyv,kxp,jblok) :: pot
         complex, dimension(nt,modesyd,modesxpd,jblok) :: pott
         end subroutine
      end interface
      interface
         subroutine PGTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,ks&
     &trt,nyv,kxp,jblok,nt,modesxpd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, ndim, kstrt, nyv, kxp
         integer :: jblok, nt, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp,jblok) :: vpot
         complex, dimension(nt,ndim,modesyd,modesxpd,jblok) :: vpott
         end subroutine
      end interface
      interface
         subroutine PPTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,ks&
     &trt,nyv,kxp,jblok,nt,modesxpd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, ndim, kstrt, nyv, kxp
         integer :: jblok, nt, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp,jblok) :: vpot
         complex, dimension(nt,ndim,modesyd,modesxpd,jblok) :: vpott
         end subroutine
      end interface
      interface
         subroutine PPOYNT2(q,exy,bxy,ffc,affp,sx,sy,sz,nx,ny,kstrt,nyv,&
     &kxp,jblok,nyhd)
         implicit none
         real :: affp, sx, sy, sz
         integer :: nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(nyv,kxp,jblok) :: q
         complex, dimension(3,nyv,kxp,jblok) :: exy, bxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PDPOYNT2(q,cu,ffc,affp,ci,sx,sy,sz,nx,ny,kstrt,nyv,k&
     &xp,jblok,nyhd)
         implicit none
         real :: affp, ci, sx, sy, sz
         integer :: nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(nyv,kxp,jblok) :: q
         complex, dimension(3,nyv,kxp,jblok) :: cu
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PDPOYNT22(q,cu,ffc,affp,ci,sx,sy,nx,ny,kstrt,nyv,kxp&
     &,jblok,nyhd)
         implicit none
         real :: affp, ci, sx, sy
         integer :: nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(nyv,kxp,jblok) :: q
         complex, dimension(2,nyv,kxp,jblok) :: cu
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PSCFGUARD2(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
         implicit none
         real :: q2m0
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cus, cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCFGUARD22(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
         implicit none
         real :: q2m0
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cus, cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSMCGUARD2(amu,nyp,x2y2m0,xym0,zxm0,zym0,nx,nxe,nypm&
     &x,nblok)
         implicit none
         real :: x2y2m0, xym0, zxm0, zym0
         integer nx, nxe, nypmx, nblok
         real, dimension(4,nxe,nypmx,nblok) :: amu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSMCGUARD22(amu,nyp,x2y2m0,xym0,nx,nxe,nypmx,nblok)
         implicit none
         real :: x2y2m0, xym0
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: amu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PAMCGUARD2X(amu,nyp,nx,nxe,nypmx,nblok,ndim)
         implicit none
         integer nx, nxe, nypmx, nblok, ndim
         real, dimension(ndim,nxe,nypmx,nblok) :: amu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCFGUARD2L(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
         implicit none
         real :: q2m0
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cus, cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCFGUARD22L(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
         implicit none
         real :: q2m0
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cus, cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSMCGUARD2L(amu,nyp,x2y2m0,xym0,zxm0,zym0,nx,nxe,nyp&
     &mx,nblok)
         implicit none
         real :: x2y2m0, xym0, zxm0, zym0
         integer nx, nxe, nypmx, nblok
         real, dimension(4,nxe,nypmx,nblok) :: amu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSMCGUARD22L(amu,nyp,x2y2m0,xym0,nx,nxe,nypmx,nblok)
         implicit none
         real :: x2y2m0, xym0
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: amu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PAMCGUARD2XL(amu,nyp,nx,nxe,nypmx,nblok,ndim)
         implicit none
         integer nx, nxe, nypmx, nblok, ndim
         real, dimension(ndim,nxe,nypmx,nblok) :: amu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PDCUPERP22(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(2,nyv,kxp,jblok) :: dcu
         complex, dimension(2,nyv,kxp,jblok) :: amu
         end subroutine
      end interface
      interface
         subroutine PADCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: dcu
         complex, dimension(4,nyv,kxp,jblok) :: amu
         end subroutine
      end interface
      interface
         subroutine PADCUPERP22(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(2,nyv,kxp,jblok) :: dcu
         complex, dimension(2,nyv,kxp,jblok) :: amu
         end subroutine
      end interface
      interface
         subroutine PDCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: dcu
         complex, dimension(4,nyv,kxp,jblok) :: amu
         end subroutine
      end interface
      interface
         subroutine PEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,&
     &ny,kstrt,nyv,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, wp0, ci, wf
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(3,nyv,kxp,jblok) :: dcu, exy
         complex, dimension(nyhd,kxp,jblok) :: ffe
         end subroutine
      end interface
      interface
         subroutine PEPOISP22(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,&
     &ny,kstrt,nyv,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, wp0, ci, wf
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(2,nyv,kxp,jblok) :: dcu, exy
         complex, dimension(nyhd,kxp,jblok) :: ffe
         end subroutine
      end interface
      interface
         subroutine PADDQEI2(qe,qi,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
!        real, dimension(*) :: qe, qi
         real :: qe, qi
         integer,  dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PADDQEI2X(qe,qi,nyp,qbme,qbmi,wpmax,wpmin,nx,nxe,nyp&
     &mx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: qbme, qbmi ,wpmax, wpmin
!        real, dimension(*) :: qe, qi
         real :: qe, qi
         integer,  dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PBADDEXT2(bxy,nyp,omx,omy,omz,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: omx, omy, omz
!        real, dimension(*) :: bxy
         real :: bxy
         integer,  dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PBADDEXT22(bz,nyp,omz,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: omz
!        real, dimension(*) :: bz
         real :: bz
         integer,  dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PIMOMENT2(qi,fxy,nyp,pxi,pyi,pzi,dt,nx,nxe,nypmx,nbl&
     &ok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: pxi, pyi, pzi, dt
!        real, dimension(*) :: qi
!        real, dimension(*) :: fxy
         real :: qi, fxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PADDVRFIELD2(a,b,c,ndim,nxe,nypmx,nblok)
         implicit none
         integer :: ndim, nxe, nypmx, nblok
         real, dimension(ndim,nxe,nypmx,nblok) :: a, b, c
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface cguard
         module procedure ipcguard2x
         module procedure ipdguard2x
      end interface
!
      interface bguard
         module procedure ipbguard2x
      end interface
!
      interface sguard
         module procedure ipscguard2
         module procedure ipsguard2
         module procedure ipsfguard2
         module procedure ipsmcguard2
      end interface
!      
      interface aguard
         module procedure ipacguard2x
         module procedure ipaguard2x
      end interface
!
      interface zguard
         module procedure ipzguard2
      end interface
!
       interface pois_init
         module procedure ippois22init
      end interface
! 
      interface pois
         module procedure ippois2
         module procedure ippois22
      end interface
!
      interface spois
         module procedure ipspois2
      end interface
!
      interface pois3
         module procedure ippois23
      end interface
!
      interface cuperp
         module procedure ipcuperp2
      end interface
!
      interface bpois
         module procedure jpbpois23
      end interface
!
      interface sbpois
         module procedure ipsbpois23
      end interface
!
      interface apois
         module procedure ipapois23
      end interface
!
      interface ibpois
         module procedure iipbpois23
      end interface
!
      interface maxwel
         module procedure ipmaxwel2
      end interface
!
      interface emfield
         module procedure ipemfield2
      end interface
!
      interface emfieldr
         module procedure ipemfieldr2
      end interface
!
      interface avpot
         module procedure ipavpot23
      end interface
!
      interface avrpot
         module procedure ipavrpot23
      end interface
!
      interface gtmodes
         module procedure ipgtmodes2
         module procedure ipgtvmodes2
      end interface
!
      interface ptmodes
         module procedure ipptmodes2
         module procedure ipptvmodes2
      end interface
!
      interface poynt
         module procedure ippoynt2
         module procedure ipdpoynt2
      end interface
!
      interface amcguard
         module procedure ipamcguard2x
      end interface
!
      interface dcuperp
         module procedure ipdcuperp23
      end interface
!
      interface adcuperp
         module procedure ipadcuperp23
      end interface
!
       interface epois_init
         module procedure ipepois23init
      end interface
!
      interface epois
         module procedure ipepois23
      end interface
!
      interface iepois
         module procedure iipepois23
      end interface
!
      interface addqei
         module procedure ipaddqei2
         module procedure ipaddqei2x
      end interface
!
      interface baddext
         module procedure ipbaddext2
      end interface
!
      interface imoment
         module procedure ipimoment2
      end interface
!
      interface addfields
         module procedure ipaddvrfield2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipcguard2x(fxy,nyp,nx,inorder)
! copy guard cells in x for non-uniform, periodic 2d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: fxy
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(fxy,2); nypmx = size(fxy,3); nblok = size(fxy,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(fxy,1))
         case (1)
            if (order==LINEAR) then
               call PDGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
            else
               call PDGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
            endif
         case (2)
            if (order==LINEAR) then
               call PCGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
            else
               call PCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
            endif
         case (3)
            if (order==LINEAR) then
               call PBGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
            else
               call PBGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
            endif
         end select
         end subroutine ipcguard2x
!
         subroutine ipdguard2x(q,nyp,nx,inorder)
! copy guard cells in x for non-uniform, periodic 2d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PDGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
         else
            call PDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipdguard2x
!
         subroutine ipbguard2x(bxy,nyp,nx,inorder)
! copy guard cells in x for non-uniform, periodic 2d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: bxy
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(bxy,2); nypmx = size(bxy,3); nblok = size(bxy,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PBGUARD2XL(bxy,nyp,nx,nxe,nypmx,nblok)
         else
            call PBGUARD2X(bxy,nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipbguard2x
!
         subroutine ipscguard2(cu,nyp,xj0,yj0,zj0,nx,inorder)
! initialize periodic 2d vector field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: xj0, yj0, zj0
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(cu,2); nypmx = size(cu,3); nblok = size(cu,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call PSCGUARD22L(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
            else
               call PSCGUARD22(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
            endif
         case (3)
            if (order==LINEAR) then
               call PSCGUARD2L(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
            else
               call PSCGUARD2(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
            endif
         end select
         end subroutine ipscguard2
!
         subroutine ipsguard2(q,nyp,qi0,nx,inorder)
! initialize periodic 2d scalar field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: qi0
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PSGUARD2L(q,nyp,qi0,nx,nxe,nypmx,nblok)
         else
            call PSGUARD2(q,nyp,qi0,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipsguard2
!
         subroutine ipacguard2x(cu,nyp,nx,inorder)
! add guard cells in x for non-uniform, periodic 2d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(cu,2); nypmx = size(cu,3); nblok = size(cu,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call PACGUARD22XL(cu,nyp,nx,nxe,nypmx,nblok)
            else
               call PACGUARD22X(cu,nyp,nx,nxe,nypmx,nblok)
            endif
         case (3)
            if (order==LINEAR) then
               call PACGUARD2XL(cu,nyp,nx,nxe,nypmx,nblok)
            else
               call PACGUARD2X(cu,nyp,nx,nxe,nypmx,nblok)
            endif
         end select
         end subroutine ipacguard2x
!
         subroutine ipaguard2x(q,nyp,nx,inorder)
! add guard cells in x for non-uniform, periodic 2d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PAGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
         else
            call PAGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipaguard2x
!
         subroutine ipzguard2(q,nyp,nx,inorder)
! zero out guard cells in periodic 2d scalar field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PZGUARD2L(q,nyp,nx,nxe,nypmx,nblok)
         else
            call PZGUARD2(q,nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipzguard2
!
         subroutine ippois2init(ffc,ax,ay,affp,nx,ny,kstrt)
! initialize 2d periodic poisson solver
         implicit none
         integer :: nx, ny, kstrt
         real :: ax, ay, affp
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = 0, nyv, kxp, jblok, nyhd
         real :: we
         complex, dimension(1,1,1) :: q, fx, fy
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kx&
     &p,jblok,nyhd)
         end subroutine ippois2init
!
         subroutine ippois2(q,fx,ffc,we,nx,ny,kstrt)
! poisson solver for periodic 2d potential
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q, fx, ffc
! local data
         integer :: isign = 1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp
         complex, dimension(1,1,1) :: fy
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kx&
     &p,jblok,nyhd)
         end subroutine ippois2
!
         subroutine ipspois2(q,fy,ffc,nx,ny,kstrt)
! smoother for periodic 2d scalar field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: q, fy, ffc
! local data
         integer :: isign = 2, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp, we
         complex, dimension(1,1,1) :: fx
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kx&
     &p,jblok,nyhd)
         end subroutine ipspois2
!
         subroutine ippois22init(ffc,ax,ay,affp,nx,ny,kstrt)
! initialize 2d periodic electric field solver
         implicit none
         integer :: nx, ny, kstrt
         real :: ax, ay, affp
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = 0, nyv, kxp, jblok, nyhd
         real :: we
         complex, dimension(1,1,1) :: q
         complex, dimension(2,1,1,1) :: fxy
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp,&
     &jblok,nyhd)
         end subroutine ippois22init
!
         subroutine ippois22(q,fxy,ffc,we,tfield,nx,ny,kstrt)
! poisson solver for periodic 2d electric field
         implicit none
         integer :: nx, ny, kstrt
         real :: we, tfield
         complex, dimension(:,:,:), pointer :: q, ffc
         complex, dimension(:,:,:,:), pointer :: fxy
! local data
         integer :: isign = -1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp, tf
         double precision :: dtime
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
! initialize timer
         call wtimer(tf,dtime,-1)
         call PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp,&
     &jblok,nyhd)
! record time
         call wtimer(tf,dtime)
         tfield = tfield + tf
         end subroutine ippois22
!
         subroutine ippois23(q,fxy,ffc,we,tfield,nx,ny,kstrt)
! poisson solver for periodic 2-1/2d electric field
         implicit none
         integer :: nx, ny, kstrt
         real :: we, tfield
         complex, dimension(:,:,:), pointer :: q, ffc
         complex, dimension(:,:,:,:), pointer :: fxy
! local data
         integer :: isign = -1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp, tf
         double precision :: dtime
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
! initialize timer
         call wtimer(tf,dtime,-1)
         select case(size(fxy,1))
         case (2)
            call PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,k&
     &xp,jblok,nyhd)
         case (3)
            call PPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,k&
     &xp,jblok,nyhd)
         end select
! record time
         call wtimer(tf,dtime)
         tfield = tfield + tf
         end subroutine ippois23
!
         subroutine ipdivf2(f,df,nx,ny,kstrt)
! calculates the divergence of periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: df
! local data
         integer :: ndim, nyv, kxp, jblok
         ndim = size(f,1)
         nyv = size(f,2); kxp = size(f,3); jblok = size(f,4)
         call PDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp,jblok)
         end subroutine ipdivf2
!
         subroutine ipgradf2(df,f,nx,ny,kstrt)
! calculates the gradient of periodic 2d scalar field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: df
         complex, dimension(:,:,:,:), pointer :: f
! local data
         integer :: ndim, nyv, kxp, jblok
         ndim = size(f,1)
         nyv = size(df,1); kxp = size(df,2); jblok = size(df,3)
         call PGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp,jblok)
         end subroutine ipgradf2
!
         subroutine ipcurlf2(f,g,nx,ny,kstrt)
! calculates the curl of periodic 2-1/2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: f, g
! local data
         integer :: nyv, kxp, jblok
         nyv = size(f,2); kxp = size(f,3); jblok = size(f,4)
         call PCURLF2(f,g,nx,ny,kstrt,nyv,kxp,jblok)
         end subroutine ipcurlf2
!
         subroutine ipcurlf22(f,g,nx,ny,kstrt)
! calculates the curl of periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
! local data
         integer :: nyv, kxp, jblok
         nyv = size(f,2); kxp = size(f,3); jblok = size(f,4)
         call PCURLF22(f,g,nx,ny,kstrt,nyv,kxp,jblok)
         end subroutine ipcurlf22
!
         subroutine ipcuperp2(cu,nx,ny,kstrt)
! calculates the transverse part of periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: cu
! local data
         integer :: nyv, kxp, jblok
         nyv = size(cu,2); kxp = size(cu,3); jblok = size(cu,4)
         select case(size(cu,1))
         case (2)
            call PCUPERP22(cu,nx,ny,kstrt,nyv,kxp,jblok)
         case (3)
            call PCUPERP2(cu,nx,ny,kstrt,nyv,kxp,jblok)
         end select
         end subroutine ipcuperp2
!
         subroutine jpbpois23(cu,bxy,ffc,ci,wm,tfield,nx,ny,kstrt)
! calculates static magnetic field for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm, tfield
         complex, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = -1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp, tf
         double precision :: dtime
         complex, dimension(2,1,1,1) :: bxy0
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
! initialize timer
         call wtimer(tf,dtime,-1)
         select case(size(cu,1))
         case (2)
            call PBPOISP22(cu,bxy0,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,&
     &kstrt,nyv,kxp,jblok,nyhd)
         case (3)
            call PBPOISP23(cu,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,kstrt&
     &,nyv,kxp,jblok,nyhd)
         end select
! record time
         call wtimer(tf,dtime)
         tfield = tfield + tf
         end subroutine jpbpois23
!
         subroutine ipsbpois23(cu,bxy,ffc,nx,ny,kstrt)
! smoother for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = 2, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp, ci, wm
         complex, dimension(1,1,1) :: bz0
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         select case(size(cu,1))
         case (2)
            call PBPOISP22(cu,bxy,bz0,isign,ffc,ax,ay,affp,ci,wm,nx,ny,k&
     &strt,nyv,kxp,jblok,nyhd)
         case (3)
            call PBPOISP23(cu,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,kstrt&
     &,nyv,kxp,jblok,nyhd)
         end select
         end subroutine ipsbpois23
!
         subroutine ipapois23(cu,axy,ffc,ci,wm,nx,ny,kstrt)
! calculates static vector potential for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: cu, axy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = 1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp
         complex, dimension(1,1,1) :: bz0
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         select case(size(cu,1))
         case (2)
            call PBPOISP22(cu,axy,bz0,isign,ffc,ax,ay,affp,ci,wm,nx,ny,k&
     &strt,nyv,kxp,jblok,nyhd)
         case (3)
            call PBPOISP23(cu,axy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,kstrt&
     &,nyv,kxp,jblok,nyhd)
         end select
         end subroutine ipapois23
!
         subroutine iipbpois23(cu,bxy,ffc,ci,wm,nx,ny,kstrt)
! calculates static magnetic field for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         complex, dimension(:,:,:), pointer :: ffc
         complex, dimension(:,:,:,:), pointer :: cu, bxy
! local data
         integer :: nyv, kxp, jblok, nyhd
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call IPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,jblok,nyhd&
     &)
         end subroutine iipbpois23
!
         subroutine ipmaxwel2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,tfield,nx,&
     &ny,kstrt)
! calculates maxwell's equation for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         real :: affp, ci, dt, wf, wm, tfield
         complex, dimension(:,:,:), pointer :: ffc
         complex, dimension(:,:,:,:), pointer :: exy, bxy, cu
! local data
         integer :: nyv, kxp, jblok, nyhd
         real :: tf
         double precision :: dtime
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
! initialize timer
         call wtimer(tf,dtime,-1)
         call PMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt,nyv,k&
     &xp,jblok,nyhd)
! record time
         call wtimer(tf,dtime)
         tfield = tfield + tf
         end subroutine ipmaxwel2
!
         subroutine ipemfield2(fxy,exy,ffc,isign,nx,ny,kstrt)
! combines and smooths periodic 2d vector fields
         implicit none
         integer :: isign, nx, ny, kstrt, nyhd
         complex, dimension(:,:,:,:), pointer :: fxy, exy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: nyv, kxp, jblok
         nyv = size(fxy,2); kxp = size(fxy,3); jblok = size(fxy,4)
         nyhd = size(ffc,1) 
         call PEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd&
     &)
         end subroutine ipemfield2 
!
         subroutine ipemfieldr2(fxy,exy,ffd,isign,nx,ny,kstrt)
! combines and smooths real 2d vector fields
         implicit none
         integer :: isign, nx, ny, kstrt, nyd
         real, dimension(:,:,:,:), pointer :: fxy, exy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: nyv, kxp2, j2blok
         nyv = size(fxy,2); j2blok = size(fxy,4)
         nyd = size(ffd,1); kxp2 = size(ffd,2); 
         call PEMFIELDR2(fxy,exy,ffd,isign,nx,ny,kstrt,nyv,kxp2,j2blok,n&
     &yd)
         end subroutine ipemfieldr2 
!
         subroutine ipavpot23(bxy,axy,nx,ny,kstrt)
! calculates periodic 2d vector potential from magnetic field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: bxy, axy
! local data
         integer :: nyv, kxp, jblok
         nyv = size(bxy,2); kxp = size(bxy,3); jblok = size(bxy,4)
         call PAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp,jblok)
         end subroutine ipavpot23
!
         subroutine ipavrpot23(axy,bxy,ffc,affp,ci,nx,ny,kstrt)
! calculates periodic 2d radiative vector potential from magnetic field
! and current
         implicit none
         integer :: nx, ny, kstrt
         real :: affp, ci
         complex, dimension(:,:,:,:), pointer :: axy, bxy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: nyv, kxp, jblok, nyhd
         nyv = size(axy,2); kxp = size(axy,3); jblok = size(axy,4)
         nyhd = size(ffc,1)
         call PAVRPOT23(axy,bxy,ffc,affp,ci,nx,ny,kstrt,nyv,kxp,jblok,ny&
     &hd)
         end subroutine ipavrpot23
!
         subroutine ipgtmodes2(pot,pott,nx,ny,modesx,modesy,kstrt)
! extracts lowest order modes from periodic 2d scalar field
         implicit none
         integer :: nx, ny, modesx, modesy, kstrt
         complex, dimension(:,:,:), pointer :: pot
         complex, dimension(:,:,:), pointer :: pott
! local data
         integer :: it, nyv, kxp, jblok, nt, modesxpd, modesyd
         it = 1; nt = 1
         nyv = size(pot,1); kxp = size(pot,2); jblok = size(pot,3)
         modesyd = size(pott,1); modesxpd = size(pott,2)
         call PGTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,kxp,jb&
     &lok,nt,modesxpd,modesyd)
         end subroutine ipgtmodes2
!
         subroutine ipptmodes2(pot,pott,nx,ny,modesx,modesy,kstrt)
! extracts lowest order modes to periodic 2d scalar field
         implicit none
         integer :: nx, ny, modesx, modesy, kstrt
         complex, dimension(:,:,:), pointer :: pot
         complex, dimension(:,:,:), pointer :: pott
! local data
         integer :: it, nyv, kxp, jblok, nt, modesxpd, modesyd
         it = 1; nt = 1
         nyv = size(pot,1); kxp = size(pot,2); jblok = size(pot,3)
         modesyd = size(pott,1); modesxpd = size(pott,2)
         call PPTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,kxp,jb&
     &lok,nt,modesxpd,modesyd)
         end subroutine ipptmodes2
!
         subroutine ipgtvmodes2(vpot,vpott,nx,ny,modesx,modesy,kstrt)
! extracts lowest order modes from periodic 2d vector field
         implicit none
         integer :: nx, ny, modesx, modesy, kstrt
         complex, dimension(:,:,:,:), pointer :: vpot
         complex, dimension(:,:,:,:), pointer :: vpott
! local data
         integer :: ndim, it, nyv, kxp, jblok, nt, modesxpd, modesyd
         it = 1; nt = 1
         ndim = size(vpot,1)
         nyv = size(vpot,2); kxp = size(vpot,3); jblok = size(vpot,4)
         modesyd = size(vpott,2); modesxpd = size(vpott,3)
         call PGTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,kstrt,ny&
     &v,kxp,jblok,nt,modesxpd,modesyd)
         end subroutine ipgtvmodes2
!
         subroutine ipptvmodes2(vpot,vpott,nx,ny,modesx,modesy,kstrt)
! extracts lowest order modes to periodic 2d vector field
         implicit none
         integer :: nx, ny, modesx, modesy, kstrt
         complex, dimension(:,:,:,:), pointer :: vpot
         complex, dimension(:,:,:,:), pointer :: vpott
! local data
         integer :: ndim, it, nyv, kxp, jblok, nt, modesxpd, modesyd
         it = 1; nt = 1
         ndim = size(vpot,1)
         nyv = size(vpot,2); kxp = size(vpot,3); jblok = size(vpot,4)
         modesyd = size(vpott,2); modesxpd = size(vpott,3)
         call PPTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,kstrt,ny&
     &v,kxp,jblok,nt,modesxpd,modesyd)
         end subroutine ipptvmodes2
!
         subroutine ippoynt2(q,exy,bxy,ffc,affp,sx,sy,sz,nx,ny,kstrt)
! calculates the momentum in the electromagnetic field
         implicit none
         integer :: nx, ny, kstrt
         real :: affp, sx, sy, sz
         complex, dimension(:,:,:), pointer :: q
         complex, dimension(:,:,:,:), pointer :: exy, bxy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: nyv, kxp, jblok, nyhd
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOYNT2(q,exy,bxy,ffc,affp,sx,sy,sz,nx,ny,kstrt,nyv,kxp,jb&
     &lok,nyhd)
         end subroutine ippoynt2
!
         subroutine ipdpoynt2(q,cu,ffc,affp,ci,sx,sy,sz,nx,ny,kstrt)
! calculates the momentum in the darwin field
         implicit none
         integer :: nx, ny, kstrt
         real :: affp, ci, sx, sy, sz
         complex, dimension(:,:,:), pointer :: q
         complex, dimension(:,:,:,:), pointer :: cu
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: nyv, kxp, jblok, nyhd
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         select case(size(cu,1))
         case (2)
            call PDPOYNT22(q,cu,ffc,affp,ci,sx,sy,nx,ny,kstrt,nyv,kxp,jb&
     &lok,nyhd)
            sz = 0.0
         case (3)
            call PDPOYNT2(q,cu,ffc,affp,ci,sx,sy,sz,nx,ny,kstrt,nyv,kxp,&
     &jblok,nyhd)
         end select
         end subroutine ipdpoynt2
!
         subroutine ipsfguard2(cus,cu,nyp,q2m0,nx,inorder)
! initialize periodic 2d vector field with scaled field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: q2m0
         real, dimension(:,:,:,:), pointer :: cu, cus
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(cus,2); nypmx = size(cus,3); nblok = size(cus,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cus,1))
         case (2)
            if (order==LINEAR) then
               call PSCFGUARD22L(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
            else
               call PSCFGUARD22(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
            endif
         case (3)
            if (order==LINEAR) then
               call PSCFGUARD2L(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
            else
               call PSCFGUARD2(cus,cu,nyp,q2m0,nx,nxe,nypmx,nblok)
            endif
         end select
         end subroutine ipsfguard2
!
         subroutine ipsmcguard2(amu,nyp,x2y2m0,xym0,zxm0,zym0,nx,inorder&
     &)
! initialize periodic 2d tensor field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: x2y2m0, xym0, zxm0, zym0
         real, dimension(:,:,:,:), pointer :: amu
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(amu,2); nypmx = size(amu,3); nblok = size(amu,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(amu,1))
         case (2)
            if (order==LINEAR) then
               call PSMCGUARD22L(amu,nyp,x2y2m0,xym0,nx,nxe,nypmx,nblok)
            else
               call PSMCGUARD22(amu,nyp,x2y2m0,xym0,nx,nxe,nypmx,nblok)
            endif
         case (4)
            if (order==LINEAR) then
               call PSMCGUARD2L(amu,nyp,x2y2m0,xym0,zxm0,zym0,nx,nxe,nyp&
     &mx,nblok)
            else
               call PSMCGUARD2(amu,nyp,x2y2m0,xym0,zxm0,zym0,nx,nxe,nypm&
     &x,nblok)
            endif
         end select
         end subroutine ipsmcguard2
!
         subroutine ipamcguard2x(amu,nyp,nx,inorder)
! add guard cells in x for non-uniform, periodic 2d tensor data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: amu
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, ndim, order
         nxe = size(amu,2); nypmx = size(amu,3); nblok = size(amu,4)
         ndim = size(amu,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PAMCGUARD2XL(amu,nyp,nx,nxe,nypmx,nblok,ndim)
         else
            call PAMCGUARD2X(amu,nyp,nx,nxe,nypmx,nblok,ndim)
         endif
         end subroutine ipamcguard2x
!
         subroutine ipdcuperp23(dcu,amu,nx,ny,kstrt)
! calculates the transverse part of periodic 2d vector field
! from momentum flux tensor
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: dcu, amu
! local data
         integer :: nyv, kxp, jblok
         nyv = size(dcu,2); kxp = size(dcu,3); jblok = size(dcu,4)
         select case(size(dcu,1))
         case (2)
            call PDCUPERP22(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
         case (3)
            call PDCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
         end select
         end subroutine ipdcuperp23
!
         subroutine ipadcuperp23(dcu,amu,nx,ny,kstrt)
! calculates the transverse part of periodic 2d vector field
! from acceleration vector and momentum flux tensor
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: dcu, amu
! local data
         integer :: nyv, kxp, jblok
         nyv = size(dcu,2); kxp = size(dcu,3); jblok = size(dcu,4)
         select case(size(dcu,1))
         case (2)
            call PADCUPERP22(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
         case (3)
            call PADCUPERP23(dcu,amu,nx,ny,kstrt,nyv,kxp,jblok)
         end select
         end subroutine ipadcuperp23
!
         subroutine ipepois23init(ffe,ax,ay,affp,wp0,ci,nx,ny,kstrt)
! initialize 2d periodic transverse electric field solver
         implicit none
         integer :: nx, ny, kstrt
         real :: ax, ay, affp, wp0, ci
         complex, dimension(:,:,:), pointer :: ffe
! local data
         integer :: isign = 0, nyv, kxp, jblok, nyhd
         real :: wf
         complex, dimension(2,1,1,1) :: dcu
         complex, dimension(2,1,1,1) :: exy
         nyv = size(dcu,2)
         nyhd = size(ffe,1); kxp = size(ffe,2); jblok = size(ffe,3)
         call PEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,kst&
     &rt,nyv,kxp,jblok,nyhd)
         end subroutine ipepois23init
!
         subroutine ipepois23(dcu,exy,ffe,affp,ci,wf,tfield,nx,ny,kstrt)
! calculates transverse electric field for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         real :: affp, ci, wf, tfield
         complex, dimension(:,:,:,:), pointer :: dcu, exy
         complex, dimension(:,:,:), pointer :: ffe
! local data
         integer :: isign = -1, nyv, kxp, jblok, nyhd
         real :: ax, ay, wp0, tf
         double precision :: dtime
         nyv = size(dcu,2)
         nyhd = size(ffe,1); kxp = size(ffe,2); jblok = size(ffe,3)
! initialize timer
         call wtimer(tf,dtime,-1)
         select case(size(dcu,1))
         case (2)
            call PEPOISP22(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,&
     &kstrt,nyv,kxp,jblok,nyhd)
         case (3)
            call PEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,&
     &kstrt,nyv,kxp,jblok,nyhd)
         end select
! record time
         call wtimer(tf,dtime)
         tfield = tfield + tf
         end subroutine ipepois23
!
         subroutine iipepois23(dcu,exy,ffe,affp,ci,wf,tfield,nx,ny,kstrt&
     &)
! calculates transverse electric field for periodic 2d vector field
! without smoothing
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wf, tfield
         complex, dimension(:,:,:,:), pointer :: dcu, exy
         complex, dimension(:,:,:), pointer :: ffe
! local data
         integer :: isign = 1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp, wp0, tf
         double precision :: dtime
         nyv = size(dcu,2)
         nyhd = size(ffe,1); kxp = size(ffe,2); jblok = size(ffe,3)
! initialize timer
         call wtimer(tf,dtime,-1)
         select case(size(dcu,1))
         case (2)
            call PEPOISP22(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,&
     &kstrt,nyv,kxp,jblok,nyhd)
         case (3)
            call PEPOISP23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,&
     &kstrt,nyv,kxp,jblok,nyhd)
         end select
! record time
         call wtimer(tf,dtime)
         tfield = tfield + tf
         end subroutine iipepois23
!
         subroutine ipaddqei2(qe,qi,nyp,nx,inorder)
! adds electron and ion densities
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: qe, qi
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(qe,1); nypmx = size(qe,2); nblok = size(qe,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PADDQEI2(qe(1,1,1),qi(1,1,1),nyp,nx,nxe,nypmx,nblok)
         else
            call PADDQEI2(qe(2,2,1),qi(2,2,1),nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipaddqei2
!
         subroutine ipaddqei2x(qe,qi,qbme,qbmi,wpmax,wpmin,nyp,nx,inorde&
     &r)
! adds electron and ion densities, and calculates maximum and minimum
! plasma frequency
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: qbme, qbmi, wpmax, wpmin
         real, dimension(:,:,:), pointer :: qe, qi
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         double precision, dimension(2) :: sum2, work2
         nxe = size(qe,1); nypmx = size(qe,2); nblok = size(qe,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PADDQEI2X(qe(1,1,1),qi(1,1,1),nyp,qbme,qbmi,wpmax,wpmin&
     &,nx,nxe,nypmx,nblok)
         else
            call PADDQEI2X(qe(2,2,1),qi(2,2,1),nyp,qbme,qbmi,wpmax,wpmin&
     &,nx,nxe,nypmx,nblok)
         endif
! maximum/minimum over the y direction
         sum2(1) = wpmax
         sum2(2) = -wpmin
         call PDMAX(sum2,work2,2,1)
         wpmax = sum2(1)
         wpmin = -sum2(2)
         end subroutine ipaddqei2x
!
         subroutine ipbaddext2(bxy,nyp,omx,omy,omz,nx,inorder)
! adds constant to magnetic field
         implicit none
         integer :: nx
         real :: omx, omy, omz
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: bxy
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(bxy,2); nypmx = size(bxy,3); nblok = size(bxy,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               call PBADDEXT22(bxy(1,1,1,1),nyp,omz,nx,nxe,nypmx,nblok)
            else
               call PBADDEXT22(bxy(1,2,2,1),nyp,omz,nx,nxe,nypmx,nblok)
            endif
         case (3)
            if (order==LINEAR) then
               call PBADDEXT2(bxy(1,1,1,1),nyp,omx,omy,omz,nx,nxe,nypmx,&
     &nblok)
            else
               call PBADDEXT2(bxy(1,2,2,1),nyp,omx,omy,omz,nx,nxe,nypmx,&
     &nblok)
            endif
         end select
         end subroutine ipbaddext2
!
         subroutine ipimoment2(qi,fxy,nyp,msg,id0,iunit,px,py,pz,dt,wx,w&
     &y,wz,nx,inorder)
! calculate ion momentum from integral of qi*fxy,
! and prints it, and adds it to total momentum, for 2 or 2-1/2d code
         implicit none
         integer :: nx, id0, iunit
         real :: px, py, pz, dt, wx, wy, wz
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: qi
         real, dimension(:,:,:,:), pointer :: fxy
         integer, dimension(:), pointer :: nyp
         double precision, dimension(:) :: msg
! local data
         integer :: nxe, nypmx, nblok, order
         real :: sx, sy, sz
         double precision, dimension(3) :: sum3, work3
  995    format (' ion momentum = ',3e14.7)
         nxe = size(fxy,2); nypmx = size(fxy,3); nblok = size(fxy,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! calculate and print ion momentum
         select case(size(fxy,1))
         case (2)
            if (order==LINEAR) then
               call PIMOMENT22(qi(1,1,1),fxy(1,1,1,1),nyp,sx,sy,dt,nx,nx&
     &e,nypmx,nblok)
            else
               call PIMOMENT22(qi(2,2,1),fxy(1,2,2,1),nyp,sx,sy,dt,nx,nx&
     &e,nypmx,nblok)
            endif
            sz = 0.0
         case (3)
            if (order==LINEAR) then
               call PIMOMENT2(qi(1,1,1),fxy(1,1,1,1),nyp,sx,sy,sz,dt,nx,&
     &nxe,nypmx,nblok)
            else
               call PIMOMENT2(qi(2,2,1),fxy(1,2,2,1),nyp,sx,sy,sz,dt,nx,&
     &nxe,nypmx,nblok)
            endif
         end select
! sum over the x, y and z directions
         sum3(1) = sx
         sum3(2) = sy
         sum3(3) = sz
         call PDSUM(sum3,work3,3,1)
         px = px + sum3(1)
         py = py + sum3(2)
         pz = pz + sum3(3)
         if (id0==0) write (iunit,995) px, py, pz
! add to total momentum
         wx = wx + px
         wy = wy + py
         wz = wz + pz
! save momentum values
         msg(1:3) = (/px,py,pz/)
         msg(4:6) = (/wx,wy,wz/)
         end subroutine ipimoment2
!
         subroutine ipaddvrfield2(a,b,c)
! calculate a = b + c for distributed real vector fields
         implicit none
         real, dimension(:,:,:,:), pointer :: a, b, c
! local data
         integer :: ndim, nxe, nypmx, nblok
         ndim = size(a,1); nxe = size(a,2)
         nypmx = size(a,3); nblok = size(a,4)
         call PADDVRFIELD2(a,b,c,ndim,nxe,nypmx,nblok)
         end subroutine ipaddvrfield2
!
      end module pfield2d

