!-----------------------------------------------------------------------
!
      module pnfield2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pnfield2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: february 17, 2006
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: poisn_init, poisn, cuperpn, bpoisn
      public :: ibpoisn, maxweln, cmfieldn, emfieldn, cpfieldn, avpotn
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PPOISNX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt&
     &,ny2d,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(ny2d,kxp2,j2blok) :: q, fx, fy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPOISNX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,&
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
         subroutine PPOISNX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,&
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
         subroutine PPOISN23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,n&
     &yv,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(nyv,kxp2,j2blok) :: q
         real, dimension(3,nyv,kxp2,j2blok) :: fxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PCUPERPNX2(cu,nx,ny,kstrt,ny2d,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok
         complex, dimension(3,ny2d,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERPN2(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(3,nyv,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERPNX22(cu,nx,ny,kstrt,ny2d,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok
         complex, dimension(2,ny2d,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PBPOISNX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,k&
     &strt,ny2d,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(3,ny2d,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PBPOISN23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,ks&
     &trt,nyv,kxp2,j2blok,nyd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(3,nyv,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PBPOISNX22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,n&
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
         subroutine IPBPOISNX23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,ny2d,kxp2,j&
     &2blok,nyd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(3,ny2d,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine IPBPOISN23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,j2b&
     &lok,nyd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(3,nyv,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PMAXWELNX2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kst&
     &rt,ny2d,kxp2,j2blok,nyd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok, nyd
         complex, dimension(3,ny2d,kxp2,j2blok) :: exy, bxy, cu
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PMAXWELN2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstr&
     &t,nyv,kxp2,j2blok,nyd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(3,nyv,kxp2,j2blok) :: exy, bxy, cu
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PDMFIELDN2(q2,q,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
         complex, dimension(ny2d,kxp2,j2blok) :: q2
         real, dimension(nyv,kxp2,j2blok) :: q
         end subroutine
      end interface
      interface
         subroutine PCMFIELDN2(cu2,cu,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
         complex, dimension(3,ny2d,kxp2,j2blok) :: cu2
         real, dimension(3,nyv,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PEMFIELDN2(fxy,exy,ffd,isign,nx,ny,kstrt,ny2d,nyv,kx&
     &p2,j2blok,nyd)
         implicit none
         integer :: isign, nx, ny, kstrt, ny2d, nyv, kxp2, j2blok, nyd
         complex, dimension(3,ny2d,kxp2,j2blok) :: fxy
         real, dimension(3,nyv,kxp2,j2blok) :: exy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PPMFIELDN2(pot2,pot,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok&
     &)
         implicit none
         integer :: nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
         complex, dimension(ny2d,kxp2,j2blok) :: pot2
         real, dimension(nyv,kxp2,j2blok) :: pot
         end subroutine
      end interface
      interface
         subroutine PCPFIELDN2(fxy,exy,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
         complex, dimension(3,ny2d,kxp2,j2blok) :: fxy
         real, dimension(3,nyv,kxp2,j2blok) :: exy
         end subroutine
      end interface
      interface
         subroutine PAVPOTNX23(bxy,axy,nx,ny,kstrt,ny2d,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok
         complex, dimension(3,ny2d,kxp2,j2blok) :: bxy, axy
         end subroutine
      end interface
      interface
         subroutine PAVPOTN23(bxy,axy,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(3,nyv,kxp2,j2blok) :: bxy, axy
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface poisn_init
         module procedure ippoisn22init
      end interface
!
      interface poisn
         module procedure ippoisn2
         module procedure ipspoisn2
         module procedure ippoisn22
      end interface
!
      interface cuperpn
         module procedure ipcuperpn2
      end interface
!
      interface bpoisn
         module procedure jpbpoisn23
      end interface
!
      interface ibpoisn
         module procedure jipbpoisn23
      end interface
!
      interface maxweln
         module procedure ipmaxweln2
      end interface
!
      interface cmfieldn
         module procedure ipcmfieldn2
         module procedure ipdmfieldn2
      end interface
!
      interface emfieldn
         module procedure ipemfieldn2
      end interface
!
      interface cpfieldn
         module procedure ipcpfieldn2
      end interface
!
      interface avpotn
         module procedure ipavpotn23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ippoisn2(q,fx,ffd,we,nx,ny,kstrt)
! poisson solver for 2d potential, neumann boundaries
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
         call PPOISNX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2d,&
     &kxp2,j2blok,nyd)
         end subroutine ippoisn2
!
         subroutine ipspoisn2(q,fy,ffd,nx,ny,kstrt)
! smoother for 2d scalar field, neumann boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: q, fy, ffd
! local data
         integer :: isign = 2, ny2d, kxp2, j2blok, nyd
         real :: ax, ay, affp, we
         complex, dimension(1,1,1) :: fx
         ny2d = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         nyd = size(ffd,1)
         call PPOISNX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2d,&
     &kxp2,j2blok,nyd)
         end subroutine ipspoisn2
!
         subroutine ippoisn22init(ffd,ax,ay,affp,nx,ny,kstrt)
! initialize 2d electric field solver, neumann boundaries
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
         call PPOISNX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2d,k&
     &xp2,j2blok,nyd)
         end subroutine ippoisn22init
!
         subroutine ippoisn22(q,fxy,ffd,we,nx,ny,kstrt)
! poisson solver for 2d electric field, neumann boundaries
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
            call PPOISNX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2&
     &d,kxp2,j2blok,nyd)
         else if (size(fxy,1)==3) then
            call PPOISNX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2&
     &d,kxp2,j2blok,nyd)
         endif
         end subroutine ippoisn22
!
         subroutine ipcuperpn2(cu,nx,ny,kstrt,kxp2)
! calculates transverse part of 2d vector field, neumann boundaries
         implicit none
         integer :: nx, ny, kstrt, kxp2
         real, dimension(:,:,:,:), pointer :: cu
! local data
         integer :: nyv, j2blok
         nyv = size(cu,2); j2blok = size(cu,4)
         call PCUPERPN2(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
         end subroutine ipcuperpn2
!
         subroutine jpbpoisn23(cu,bxy,ffd,ci,wm,nx,ny,kstrt)
! caculates static vector potential for 2d vector field,
! neumann boundaries
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
         call PBPOISNX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,kstrt,n&
     &y2d,kxp2,j2blok,nyd)
         end subroutine jpbpoisn23
!
         subroutine jipbpoisn23(cu,bxy,ffd,ci,wm,nx,ny,kstrt)
! calculates static magnetic field for periodic 2d vector field
! neumann boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         real, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: nyv, kxp2, j2blok, nyd
         nyv= size(cu,2); j2blok = size(cu,4)
         nyd = size(ffd,1); kxp2 = size(ffd,2)
         call IPBPOISN23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,j2blok,ny&
     &d)
         end subroutine jipbpoisn23
!
         subroutine ipmaxweln2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kst&
     &rt)
! calculates maxwell's equation for 2d vector field,
! neumann boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: affp, ci, dt, wf, wm
         real, dimension(:,:,:,:), pointer :: exy, bxy, cu
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: nyv, kxp2, j2blok, nyd
         nyv = size(cu,2); j2blok = size(cu,4)
         nyd = size(ffd,1); kxp2 = size(ffd,2)
         call PMAXWELN2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstrt,nyv,&
     &kxp2,j2blok,nyd)
         end subroutine ipmaxweln2
!
         subroutine ipcmfieldn2(cu2,cu,nx,ny,kstrt)
! copies from double to normal array in y dimension for 2d vector data
! neumann boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: cu2
         real, dimension(:,:,:,:), pointer :: cu
! local data
         integer :: ny2d, kxp2, j2blok, nyv
         ny2d = size(cu2,2); kxp2 = size(cu2,3); j2blok = size(cu2,4)
         nyv = size(cu,2)
         call PCMFIELDN2(cu2,cu,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         end subroutine ipcmfieldn2
!
         subroutine ipdmfieldn2(q2,q,nx,ny,kstrt)
! copies from double to normal array in y dimension for 2d scalar data
! neumann boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: q2
         real, dimension(:,:,:), pointer :: q
! local data
         integer :: ny2d, kxp2, j2blok, nyv
         ny2d = size(q2,1); kxp2 = size(q2,2); j2blok = size(q2,3)
         nyv = size(q,1)
         call PDMFIELDN2(q2,q,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         end subroutine ipdmfieldn2
!
         subroutine ipemfieldn2(fxy,exy,ffd,isign,nx,ny,kstrt)
! combines and smooths 2d vector fields, neumann boundaries
         implicit none
         integer :: isign, nx, ny, kstrt, nyd
         complex, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:,:), pointer :: exy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: ny2d, kxp2, j2blok, nyv
         ny2d = size(fxy,2); kxp2 = size(fxy,3); j2blok = size(fxy,4)
         nyv = size(exy,2); nyd = size(ffd,1)
         call PEMFIELDN2(fxy,exy,ffd,isign,nx,ny,kstrt,ny2d,nyv,kxp2,j2b&
     &lok,nyd)
         end subroutine ipemfieldn2
!
         subroutine ipcpfieldn2(fxy,exy,nx,ny,kstrt)
! combines 2d electric fields, neumann boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:,:), pointer :: exy
! local data
         integer :: ny2d, kxp2, j2blok, nyv
         ny2d = size(fxy,2); kxp2 = size(fxy,3); j2blok = size(fxy,4)
         nyv = size(exy,2)
         call PCPFIELDN2(fxy,exy,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
         end subroutine ipcpfieldn2
!
         subroutine ipavpotn23(bxy,axy,nx,ny,kstrt,kxp2)
! calculates 2d vector potential from magnetic field
! neumann boundaries
         implicit none
         integer :: nx, ny, kstrt, kxp2
         real, dimension(:,:,:,:), pointer :: bxy, axy
! local data
         integer :: nyv, j2blok
         nyv = size(bxy,2); j2blok = size(bxy,4)
         call PAVPOTN23(bxy,axy,nx,ny,kstrt,nyv,kxp2,j2blok)
         end subroutine ipavpotn23
!
      end module pnfield2d
