!-----------------------------------------------------------------------
!
      module pdfield2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pdfield2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: february 25, 2006
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: PLCGUARD2X, PLDGUARD2X, PLBGUARD2X, PLSCGUARD2X
      public :: PLSCGUARD22X, PLSGUARD2X, PLACGUARD2X, PLACGUARD22X
      public :: PLAGUARD2X, PLSCGUARD2XL, PLSCGUARD22XL, PLSGUARD2XL
      public :: laguard, lcguard
      public :: poisd_init, poisdx, poisd, cuperpd, bpoisdx, bpoisd
      public :: ibpoisd, maxweld, cmfieldd, emfieldd, cpfieldd, avpotd
      public :: ipdivfd2, ipgradfd2, ipcurlfd2, ipcurlfd22
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PLCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLBGUARD2X(bxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD2X(cu,kstrt,nvp,noff,nyp,xj0,yj0,zj0,nx,ny,&
     &ngx,ngy,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: kstrt, nvp, nx, ny, ngx, ngy, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: noff, nyp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD22X(cu,kstrt,nvp,noff,nyp,xj0,yj0,nx,ny,ngx&
     &,ngy,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer :: kstrt, nvp, nx, ny, ngx, ngy, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: noff, nyp
         end subroutine
      end interface
      interface
         subroutine PLSGUARD2X(q,kstrt,nvp,noff,nyp,qi0,nx,ny,ngx,ngy,nx&
     &e,nypmx,nblok)
         implicit none
         real :: qi0
         integer :: kstrt, nvp, nx, ny, ngx, ngy, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: noff, nyp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD2(cu,kstrt,nvp,nyp,xj0,yj0,zj0,nx,ngx,ngy,n&
     &xe,nypmx,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: kstrt, nvp, nx, ngx, ngy, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD22(cu,kstrt,nvp,nyp,xj0,yj0,nx,ngx,ngy,nxe,&
     &nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer :: kstrt, nvp, nx, ngx, ngy, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLSGUARD2(q,kstrt,nvp,nyp,qi0,nx,ngx,ngy,nxe,nypmx,n&
     &blok)
         implicit none
         real :: qi0
         integer :: kstrt, nvp, nx, ngx, ngy, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLACGUARD2X(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLACGUARD22X(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLAGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD2XL(cu,kstrt,nvp,noff,nyp,xj0,yj0,zj0,nx,ny&
     &,ngx,ngy,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: kstrt, nvp, nx, ny, ngx, ngy, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: noff, nyp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD22XL(cu,kstrt,nvp,noff,nyp,xj0,yj0,nx,ny,ng&
     &x,ngy,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer :: kstrt, nvp, nx, ny, ngx, ngy, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: noff, nyp
         end subroutine
      end interface
      interface
         subroutine PLSGUARD2XL(q,kstrt,nvp,noff,nyp,qi0,nx,ny,ngx,ngy,n&
     &xe,nypmx,nblok)
         implicit none
         real :: qi0
         integer :: kstrt, nvp, nx, ny, ngx, ngy, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: noff, nyp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD2L(cu,kstrt,nvp,nyp,xj0,yj0,zj0,nx,ngx,ngy,&
     &nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer :: kstrt, nvp, nx, ngx, ngy, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLSCGUARD22L(cu,kstrt,nvp,nyp,xj0,yj0,nx,ngx,ngy,nxe&
     &,nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer :: kstrt, nvp, nx, ngx, ngy, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PLSGUARD2L(q,kstrt,nvp,nyp,qi0,nx,ngx,ngy,nxe,nypmx,&
     &nblok)
         implicit none
         real :: qi0
         integer :: kstrt, nvp, nx, ngx, ngy, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PPOISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt&
     &,ny2d,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(ny2d,kxp2,j2blok) :: q, fx, fy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPOISD2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(nyv,kxp2+1,j2blok) :: q, fx, fy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPOISDX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,&
     &ny2d,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(ny2d,kxp2,j2blok) :: q
         complex, dimension(2,ny2d,kxp2,j2blok) :: fxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,n&
     &yv,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(nyv,kxp2+1,j2blok) :: q
         real, dimension(2,nyv,kxp2+1,j2blok) :: fxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPOISDX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,&
     &ny2d,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(ny2d,kxp2,j2blok) :: q
         complex, dimension(3,ny2d,kxp2,j2blok) :: fxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPOISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,n&
     &yv,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(nyv,kxp2+1,j2blok) :: q
         real, dimension(3,nyv,kxp2+1,j2blok) :: fxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ndim, nyv, kxp2, j2blok
         real, dimension(3,nyv,kxp2+1,j2blok) :: f
         real, dimension(nyv,kxp2+1,j2blok) :: df
         end subroutine
      end interface
      interface
         subroutine PGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ndim, nyv, kxp2, j2blok
         real, dimension(nyv,kxp2+1,j2blok) :: df
         real, dimension(3,nyv,kxp2+1,j2blok) :: f
         end subroutine
      end interface
      interface
         subroutine PCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(3,nyv,kxp2+1,j2blok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PCURLFD22(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(2,nyv,kxp2+1,j2blok) :: f
         real, dimension(nyv,kxp2+1,j2blok) :: g
         end subroutine
      end interface
      interface
         subroutine PCUPERPDX2(cu,nx,ny,kstrt,ny2d,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok
         complex, dimension(3,ny2d,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERPD2(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(3,nyv,kxp2+1,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERPDX22(cu,nx,ny,kstrt,ny2d,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok
         complex, dimension(2,ny2d,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERPD22(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(2,nyv,kxp2+1,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PBPOISDX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,k&
     &strt,ny2d,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(3,ny2d,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PBPOISD23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,ks&
     &trt,nyv,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(3,nyv,kxp2+1,j2blok) :: cu, bxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PBPOISDX22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,n&
     &y,kstrt,ny2d,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(2,ny2d,kxp2,j2blok) :: cu, bxy
         complex, dimension(ny2d,kxp2,j2blok) :: bz
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PBPOISD22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny&
     &,kstrt,nyv,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(2,nyv,kxp2+1,j2blok) :: cu, bxy
         real, dimension(nyv,kxp2+1,j2blok) :: bz
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine IPBPOISDX23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,ny2d,kxp2,j&
     &2blok,nyd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(3,ny2d,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine IPBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,j2b&
     &lok,nyd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(3,nyv,kxp2+1,j2blok) :: cu, bxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PMAXWELDX2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kst&
     &rt,ny2d,kxp2,j2blok,nyd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(3,ny2d,kxp2,j2blok) :: exy, bxy, cu
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PMAXWELD2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstr&
     &t,nyv,kxp2,j2blok,nyd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(3,nyv,kxp2+1,j2blok) :: exy, bxy, cu
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PDMFIELDD2(q2,q,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
         complex, dimension(ny2d,kxp2,j2blok) :: q2
         real, dimension(nyv,kxp2+1,j2blok) :: q
         end subroutine
      end interface
      interface
         subroutine PCMFIELDD2(cu2,cu,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
         complex, dimension(3,ny2d,kxp2,j2blok) :: cu2
         real, dimension(3,nyv,kxp2+1,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PEMFIELDD2(fxy,exy,ffd,isign,nx,ny,kstrt,ny2d,nyv,kx&
     &p2,j2blok,nyd)
         implicit none
         integer :: isign, nx, ny, kstrt, ny2d, nyv, kxp2, j2blok, nyd
         complex, dimension(3,ny2d,kxp2,j2blok) :: fxy
         real, dimension(3,nyv,kxp2+1,j2blok) :: exy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPMFIELDD2(pot2,pot,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok&
     &)
         implicit none
         integer :: nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
         complex, dimension(ny2d,kxp2,j2blok) :: pot2
         real, dimension(nyv,kxp2+1,j2blok) :: pot
         end subroutine
      end interface
      interface
         subroutine PCPFIELDD2(fxy,exy,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
         complex, dimension(3,ny2d,kxp2,j2blok) :: fxy
         real, dimension(3,nyv,kxp2+1,j2blok) :: exy
         end subroutine
      end interface
      interface
         subroutine PAVPOTDX23(bxy,axy,nx,ny,kstrt,ny2d,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok
         complex, dimension(3,ny2d,kxp2,j2blok) :: bxy, axy
         end subroutine
      end interface
      interface
         subroutine PAVPOTD23(bxy,axy,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(3,nyv,kxp2+1,j2blok) :: bxy, axy
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!      
      interface laguard
         module procedure iplaguard2x
      end interface
!
      interface lcguard
         module procedure iplcguard2x
         module procedure ipldguard2x
      end interface
!
      interface poisd_init
         module procedure ippoisd22init
      end interface
!
      interface poisdx
         module procedure ippoisdx2
         module procedure ipspoisdx2
         module procedure ippoisdx22
      end interface
!
      interface poisd
         module procedure ippoisd2
         module procedure ipspoisd2
         module procedure ippoisd23
      end interface
!
      interface cuperpd
         module procedure ipcuperpd2
      end interface
!
      interface bpoisdx
         module procedure jpbpoisdx23
      end interface
!
      interface bpoisd
         module procedure jpbpoisd23
      end interface
!
      interface ibpoisd
         module procedure jipbpoisd23
      end interface
!
      interface maxweld
         module procedure ipmaxweld2
      end interface
!
      interface cmfieldd
         module procedure ipcmfieldd2
         module procedure ipdmfieldd2
      end interface
!
      interface emfieldd
         module procedure ipemfieldd2
      end interface
!
      interface cpfieldd
         module procedure ipcpfieldd2
      end interface
!
      interface avpotd
         module procedure ipavpotd23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine iplaguard2x(q,nyp,nx,inorder)
! add guard cells in x for non-uniform 2d scalar data
! disable quadratic interpolation at edge
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
         if (order==QUADRATIC) then
            call PLAGUARD2X(q(2,1,1),nyp,nx-2,nxe,nypmx,nblok)
         endif
         end subroutine iplaguard2x
!
         subroutine iplcguard2x(fxy,nyp,nx,inorder)
! copy guard cells in x for non-uniform 2d vector data
! disable quadratic interpolation at edge
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
         if (size(fxy,1)==2) then
            if (order==QUADRATIC) then
               call PLCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
            endif
         else if (size(fxy,1)==3) then
            if (order==QUADRATIC) then
               call PLBGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
            endif
         endif
         end subroutine iplcguard2x
!
         subroutine ipldguard2x(q,nyp,nx,inorder)
! copy guard cells in x for non-uniform 2d scalar data
! disable quadratic interpolation at edge
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
         if (order==QUADRATIC) then
            call PLDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipldguard2x
!
         subroutine ippoisdx2(q,fx,ffd,we,nx,ny,kstrt)
! poisson solver for 2d potential, conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q, fx, ffd
! local data
         integer :: isign = 1, ny2d, kxp2, j2blok, nyd
         real :: ax, ay, affp
         complex, dimension(1,1,1) :: fy
         ny2d = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         nyd = size(ffd,1)
         call PPOISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2d,&
     &kxp2,j2blok,nyd)
         end subroutine ippoisdx2
!
         subroutine ipspoisdx2(q,fy,ffd,nx,ny,kstrt)
! smoother for 2d scalar field, conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: q, fy, ffd
! local data
         integer :: isign = 2, ny2d, kxp2, j2blok, nyd
         real :: ax, ay, affp, we
         complex, dimension(1,1,1) :: fx
         ny2d = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         nyd = size(ffd,1)
         call PPOISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2d,&
     &kxp2,j2blok,nyd)
         end subroutine ipspoisdx2
!
         subroutine ippoisd22init(ffd,ax,ay,affp,nx,ny,kstrt)
! initialize 2d electric field solver, conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: ax, ay, affp
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: isign = 0, ny2d, kxp2, j2blok, nyd
         real :: we
         complex, dimension(1,1,1) :: q
         complex, dimension(2,1,1,1) :: fxy
         ny2d = size(q,1)
         nyd = size(ffd,1); kxp2 = size(ffd,2); j2blok = size(ffd,3)
         call PPOISDX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2d,k&
     &xp2,j2blok,nyd)
         end subroutine ippoisd22init
!
         subroutine ippoisdx22(q,fxy,ffd,we,nx,ny,kstrt)
! poisson solver for 2d electric field, conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q, ffd
         complex, dimension(:,:,:,:), pointer :: fxy
! local data
         integer :: isign = -1, ny2d, kxp2, j2blok, nyd
         real :: ax, ay, affp
         ny2d = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         nyd = size(ffd,1)
         if (size(fxy,1)==2) then
            call PPOISDX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2&
     &d,kxp2,j2blok,nyd)
         else if (size(fxy,1)==3) then
            call PPOISDX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2&
     &d,kxp2,j2blok,nyd)
         endif
         end subroutine ippoisdx22
!
         subroutine ippoisd2(q,fx,ffd,we,nx,ny,kstrt)
! poisson solver for 2d potential, conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         real, dimension(:,:,:), pointer :: q, fx
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: isign = 1, nyv, kxp2, j2blok, nyd
         real :: ax, ay, affp
         real, dimension(1,1,1) :: fy
         nyv = size(q,1); j2blok = size(q,3)
         nyd = size(ffd,1); kxp2 = size(ffd,2)
         call PPOISD2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,kx&
     &p2,j2blok,nyd)
         end subroutine ippoisd2
!
         subroutine ipspoisd2(q,fy,ffd,nx,ny,kstrt)
! smoother for 2d scalar field, conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real, dimension(:,:,:), pointer :: q, fy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: isign = 2, nyv, kxp2, j2blok, nyd
         real :: ax, ay, affp, we
         real, dimension(1,1,1) :: fx
         nyv = size(q,1); j2blok = size(q,3)
         nyd = size(ffd,1); kxp2 = size(ffd,2)
         call PPOISD2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,kx&
     &p2,j2blok,nyd)
         end subroutine ipspoisd2
!
         subroutine ippoisd23(q,fxy,ffd,we,nx,ny,kstrt)
! poisson solver for 2-1/2d electric field, conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         real, dimension(:,:,:), pointer :: q
         real, dimension(:,:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: isign = -1, nyv, kxp2, j2blok, nyd
         real :: ax, ay, affp
         nyv = size(q,1); j2blok = size(q,3)
         nyd = size(ffd,1); kxp2 = size(ffd,2)
         if (size(fxy,1)==2) then
            call PPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,&
     &kxp2,j2blok,nyd)
         else if (size(fxy,1)==3) then
            call PPOISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,&
     &kxp2,j2blok,nyd)
         endif
         end subroutine ippoisd23
!
         subroutine ipdivfd2(f,df,nx,ny,kstrt,kxp2)
! calculates the divergence of periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt, kxp2
         real, dimension(:,:,:,:), pointer :: f
         real, dimension(:,:,:), pointer :: df
! local data
         integer :: ndim, nyv, j2blok
         ndim = size(f,1); nyv = size(f,2); j2blok = size(f,4)
         call PDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
         end subroutine ipdivfd2
!
         subroutine ipgradfd2(df,f,nx,ny,kstrt,kxp2)
! calculates the gradient of periodic 2d scalar field
         implicit none
         integer :: nx, ny, kstrt, kxp2
         real, dimension(:,:,:), pointer :: df
         real, dimension(:,:,:,:), pointer :: f
! local data
         integer :: ndim, nyv, j2blok
         ndim = size(f,1); nyv = size(df,1); j2blok = size(df,3)
         call PGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
         end subroutine ipgradfd2
!
         subroutine ipcurlfd2(f,g,nx,ny,kstrt,kxp2)
! calculates the curl of periodic 2-1/2d vector field
         implicit none
         integer :: nx, ny, kstrt, kxp2
         real, dimension(:,:,:,:), pointer :: f, g
! local data
         integer :: nyv, j2blok
         nyv = size(f,2); j2blok = size(f,4)
         call PCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
         end subroutine ipcurlfd2
!
         subroutine ipcurlfd22(f,g,nx,ny,kstrt,kxp2)
! calculates the curl of periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt, kxp2
         real, dimension(:,:,:,:), pointer :: f
         real, dimension(:,:,:), pointer :: g
! local data
         integer :: nyv, j2blok
         nyv = size(f,2); j2blok = size(f,4)
         call PCURLFD22(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
         end subroutine ipcurlfd22
!
         subroutine ipcuperpd2(cu,nx,ny,kstrt,kxp2)
! calculates transverse part of 2d vector field, conducting boundaries
         implicit none
         integer :: nx, ny, kstrt, kxp2
         real, dimension(:,:,:,:), pointer :: cu
! local data
         integer :: nyv, j2blok
         nyv = size(cu,2); j2blok = size(cu,4)
         if (size(cu,1)==2) then
            call PCUPERPD22(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
         else if (size(cu,1)==3) then
            call PCUPERPD2(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
         endif
         end subroutine ipcuperpd2
!
         subroutine jpbpoisdx23(cu,bxy,ffd,ci,wm,nx,ny,kstrt)
! calculates static vector potential for 2d vector field,
! conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: isign = 1, ny2d, kxp2, j2blok, nyd
         real :: ax, ay, affp
         ny2d = size(cu,2); kxp2 = size(cu,3); j2blok = size(cu,4)
         nyd = size(ffd,1)
         if (size(cu,1)==3) then
            call PBPOISDX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,kstr&
     &t,ny2d,kxp2,j2blok,nyd)
         endif
         end subroutine jpbpoisdx23
!
         subroutine jpbpoisd23(cu,bxy,ffd,ci,wm,nx,ny,kstrt)
! calculates static vector potential for 2d vector field,
! conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         real, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: isign = 1, nyv, kxp2, j2blok, nyd
         real :: ax, ay, affp
         nyv = size(cu,2); j2blok = size(cu,4)
         nyd = size(ffd,1); kxp2 = size(ffd,2)
         if (size(cu,1)==3) then
            call PBPOISD23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,kstrt&
     &,nyv,kxp2,j2blok,nyd)
         endif
         end subroutine jpbpoisd23
!
         subroutine jipbpoisd23(cu,bxy,ffd,ci,wm,nx,ny,kstrt)
! calculates static magnetic field for periodic 2d vector field
! conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         real, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: nyv, kxp2, j2blok, nyd
         nyv= size(cu,2); j2blok = size(cu,4)
         nyd = size(ffd,1); kxp2 = size(ffd,2)
         call IPBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,j2blok,ny&
     &d)
         end subroutine jipbpoisd23
!
         subroutine ipmaxweld2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kst&
     &rt)
! calculates maxwell's equation for 2d vector field,
! conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: affp, ci, dt, wf, wm
         real, dimension(:,:,:,:), pointer :: exy, bxy, cu
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: nyv, kxp2, j2blok, nyd
         nyv = size(cu,2); j2blok = size(cu,4)
         nyd = size(ffd,1); kxp2 = size(ffd,2)
         call PMAXWELD2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstrt,nyv,&
     &kxp2,j2blok,nyd)
         end subroutine ipmaxweld2
!
         subroutine ipcmfieldd2(cu2,cu,nx,ny,kstrt)
! copies from double to normal array in y dimension for 2d vector data
! conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: cu2
         real, dimension(:,:,:,:), pointer :: cu
! local data
         integer :: ny2d, kxp2, j2blok, nyv
         ny2d = size(cu2,2); kxp2 = size(cu2,3); j2blok = size(cu2,4)
         nyv = size(cu,2)
         call PCMFIELDD2(cu2,cu,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         end subroutine ipcmfieldd2
!
         subroutine ipdmfieldd2(q2,q,nx,ny,kstrt)
! copies from double to normal array in y dimension for 2d scalar data
! conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: q2
         real, dimension(:,:,:), pointer :: q
! local data
         integer :: ny2d, kxp2, j2blok, nyv
         ny2d = size(q2,1); kxp2 = size(q2,2); j2blok = size(q2,3)
         nyv = size(q,1)
         call PDMFIELDD2(q2,q,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         end subroutine ipdmfieldd2
!
         subroutine ipemfieldd2(fxy,exy,ffd,isign,nx,ny,kstrt)
! combines and smooths 2d vector fields, conducting boundaries
         implicit none
         integer :: isign, nx, ny, kstrt, nyd
         complex, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:,:), pointer :: exy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: ny2d, kxp2, j2blok, nyv
         ny2d = size(fxy,2); kxp2 = size(fxy,3); j2blok = size(fxy,4)
         nyv = size(exy,2); nyd = size(ffd,1)
         call PEMFIELDD2(fxy,exy,ffd,isign,nx,ny,kstrt,ny2d,nyv,kxp2,j2b&
     &lok,nyd)
         end subroutine ipemfieldd2
!
         subroutine ipcpfieldd2(fxy,exy,nx,ny,kstrt)
! combines 2d electric fields, conducting boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:,:), pointer :: exy
! local data
         integer :: ny2d, kxp2, j2blok, nyv
         ny2d = size(fxy,2); kxp2 = size(fxy,3); j2blok = size(fxy,4)
         nyv = size(exy,2)
         call PCPFIELDD2(fxy,exy,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         end subroutine ipcpfieldd2
!
         subroutine ipavpotd23(bxy,axy,nx,ny,kstrt,kxp2)
! calculates 2d vector potential from magnetic field
! conducting boundaries
         implicit none
         integer :: nx, ny, kstrt, kxp2
         real, dimension(:,:,:,:), pointer :: bxy, axy
! local data
         integer :: nyv, j2blok
         nyv = size(bxy,2); j2blok = size(bxy,4)
         call PAVPOTD23(bxy,axy,nx,ny,kstrt,nyv,kxp2,j2blok)
         end subroutine ipavpotd23
!
      end module pdfield2d
