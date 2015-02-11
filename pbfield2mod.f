!-----------------------------------------------------------------------
!
      module pbfield2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pbfield2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: january 17, 2006
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: PMSCGUARD2, PMSCGUARD22, PMSGUARD2
      public :: PMSCGUARD2L, PMSCGUARD22L, PMSGUARD2L
      public :: sglsin, hafsgl, poism_init, poism, poism3, cuperpm
      public :: bpoism, ibpoism, maxwelm, cmfieldm, emfieldm, cpfieldm
      public :: avpotm
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PMSCGUARD2(cu,nyp,xj0,yj0,zj0,nx,ngx,nxe,nypmx,nblok&
     &)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, ngx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PMSCGUARD22(cu,nyp,xj0,yj0,nx,ngx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer :: nx, ngx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PMSGUARD2(q,nyp,qi0,nx,ngx,nxe,nypmx,nblok)
         implicit none
         real :: qi0
         integer :: nx, ngx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PMSCGUARD2L(cu,nyp,xj0,yj0,zj0,nx,ngx,nxe,nypmx,nblo&
     &k)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, ngx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PMSCGUARD22L(cu,nyp,xj0,yj0,nx,ngx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer :: nx, ngx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PMSGUARD2L(q,nyp,qi0,nx,ngx,nxe,nypmx,nblok)
         implicit none
         real :: qi0
         integer :: nx, ngx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSGLSIN2C(cu,cu1,nx,kyp,nxv,nypmx,kypd,kblok)
         implicit none
         integer :: nx, kyp, nxv, nypmx, kypd, kblok
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(2,2*nxv,kypd,kblok) :: cu1
         end subroutine
      end interface
      interface
         subroutine PSGLSIN2D(q,q1,nx,kyp,nxv,nypmx,kypd,kblok)
         implicit none
         integer :: nx, kyp, nxv, nypmx, kypd, kblok
!        real, dimension(*) :: q
         real :: q
         real, dimension(2*nxv,kypd,kblok) :: q1
         end subroutine
      end interface
      interface
         subroutine PSGLSIN2B(cu,cu1,nx,kyp,nxv,nypmx,kypd,kblok)
         implicit none
         integer :: nx, kyp, nxv, nypmx, kypd, kblok
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,2*nxv,kypd,kblok) :: cu1
         end subroutine
      end interface
      interface
         subroutine PHAFSGL2C(fxy,fxy1,nx,kyp,nxv,nypmx,kypd,kblok)
         implicit none
         integer :: nx, kyp, nxv, nypmx, kypd, kblok
!        real, dimension(*) :: fxy
         real :: fxy
         real, dimension(2,2*nxv,kypd,kblok) :: fxy1
         end subroutine
      end interface
      interface
         subroutine PHAFSGL2D(q,q1,nx,kyp,nxv,nypmx,kypd,kblok)
         implicit none
         integer :: nx, kyp, nxv, nypmx, kypd, kblok
!        real, dimension(*) :: q
         real :: q
         real, dimension(2*nxv,kypd,kblok) :: q1
         end subroutine
      end interface
      interface
         subroutine PHAFSGL2B(fxy,fxy1,nx,kyp,nxv,nypmx,kypd,kblok)
         implicit none
         integer :: nx, kyp, nxv, nypmx, kypd, kblok
!        real, dimension(*) :: fxy
         real :: fxy
         real, dimension(3,2*nxv,kypd,kblok) :: fxy1
         end subroutine
      end interface
      interface
         subroutine PPOISMX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt&
     &,nyv,kxp2,j2blok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyhd
         complex, dimension(nyv,kxp2,j2blok) :: q, fx, fy
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PPOISMX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp2,j2blok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyhd
         complex, dimension(nyv,kxp2,j2blok) :: q
         complex, dimension(2,nyv,kxp2,j2blok) :: fxy
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PPOISMX23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp2,j2blok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyhd
         complex, dimension(nyv,kxp2,j2blok) :: q
         complex, dimension(3,nyv,kxp2,j2blok) :: fxy
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PPOISM23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,n&
     &yvh,kxp2,j2blok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyvh, kxp2, j2blok, nyhd
         real, dimension(2*nyvh,kxp2,j2blok) :: q
         real, dimension(3,2*nyvh,kxp2,j2blok) :: fxy
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PCUPERPMX2(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         complex, dimension(3,nyv,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERPM2(cu,nx,ny,kstrt,nyvh,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyvh, kxp2, j2blok
         real, dimension(3,2*nyvh,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERPMX22(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         real, dimension(2,2*nyv,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PBPOISMX23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,k&
     &strt,nyv,kxp2,j2blok,nyhd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyhd
         complex, dimension(3,nyv,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PBPOISM23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,ks&
     &trt,nyvh,kxp2,j2blok,nyhd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyvh, kxp2, j2blok, nyhd
         real, dimension(3,2*nyvh,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PBPOISMX22(cu,bxy,bz,isign,ffb,ax,ay,affp,ci,wm,nx,n&
     &y,kstrt,nyv,kxp2,j2blok,nyhd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyhd
         complex, dimension(2,nyv,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyv,kxp2,j2blok) :: bz
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine IPBPOISMX23(cu,bxy,ffb,ci,wm,nx,ny,kstrt,nyv,kxp2,j2&
     &blok,nyhd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok, nyhd
         complex, dimension(3,nyv,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine IPBPOISM23(cu,bxy,ffb,ci,wm,nx,ny,kstrt,nyvh,kxp2,j2&
     &blok,nyhd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, kstrt, nyvh, kxp2, j2blok, nyhd
         real, dimension(3,2*nyvh,kxp2,j2blok) :: cu, bxy
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PMAXWELMX2(exy,bxy,cu,ffb,affp,ci,dt,wf,wm,nx,ny,kst&
     &rt,nyv,kxp2,j2blok,nyhd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok, nyhd
         complex, dimension(3,nyv,kxp2,j2blok) :: exy, bxy, cu
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PMAXWELM2(exy,bxy,cu,ffb,affp,ci,dt,wf,wm,nx,ny,kstr&
     &t,nyvh,kxp2,j2blok,nyhd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, kstrt, nyvh, kxp2, j2blok, nyhd
         real, dimension(3,2*nyvh,kxp2,j2blok) :: exy, bxy, cu
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PDMFIELDM2(q1,q,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, nyvh, kxp2, j2blok
         complex, dimension(nyv,kxp2,j2blok) :: q1
         real, dimension(2*nyvh,kxp2,j2blok) :: q
         end subroutine
      end interface
      interface
         subroutine PCMFIELDM2(cu1,cu,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, nyvh, kxp2, j2blok
         complex, dimension(3,nyv,kxp2,j2blok) :: cu1
         real, dimension(3,2*nyvh,kxp2,j2blok) :: cu
         end subroutine
      end interface
      interface
         subroutine PEMFIELDM2(fxy,exy,ffb,isign,nx,ny,kstrt,nyv,nyvh,kx&
     &p2,j2blok,nyhd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, nyvh, kxp2, j2blok, nyhd
         complex, dimension(3,nyv,kxp2,j2blok) :: fxy
         real, dimension(3,2*nyvh,kxp2,j2blok) :: exy
         complex, dimension(nyhd,kxp2,j2blok) :: ffb
         end subroutine
      end interface
      interface
         subroutine PPMFIELDM2(pot1,pot,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok&
     &)
         implicit none
         integer :: nx, ny, kstrt, nyv, nyvh, kxp2, j2blok
         complex, dimension(nyv,kxp2,j2blok) :: pot1
         real, dimension(2*nyvh,kxp2,j2blok) :: pot
         end subroutine
      end interface
      interface
         subroutine PCPFIELDM2(fxy,exy,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, nyvh, kxp2, j2blok
         complex, dimension(3,nyv,kxp2,j2blok) :: fxy
         real, dimension(3,2*nyvh,kxp2,j2blok) :: exy
         end subroutine
      end interface
      interface
         subroutine PAVPOTMX23(bxy,axy,nx,ny,kstrt,nyv,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp2, j2blok
         complex, dimension(3,nyv,kxp2,j2blok) :: bxy, axy
         end subroutine
      end interface
      interface
         subroutine PAVPOTM23(bxy,axy,nx,ny,kstrt,nyvh,kxp2,j2blok)
         implicit none
         integer :: nx, ny, kstrt, nyvh, kxp2, j2blok
         real, dimension(3,2*nyvh,kxp2,j2blok) :: bxy, axy
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface sglsin
         module procedure ipsglsin2c
         module procedure ipsglsin2d
      end interface
!
      interface hafsgl
         module procedure iphafsgl2c
         module procedure iphafsgl2d
      end interface
!
      interface poism_init
         module procedure ippoism22init
      end interface
!
      interface poism
         module procedure ippoism2
         module procedure ipspoism2
         module procedure ippoism22
      end interface
!
      interface poism3
         module procedure ippoism23
      end interface
!
      interface cuperpm
         module procedure ipcuperpm2
      end interface
!
      interface bpoism
         module procedure jpbpoism23
      end interface
!
      interface ibpoism
         module procedure jipbpoism23
      end interface
!
      interface maxwelm
         module procedure ipmaxwelm2
      end interface
!
      interface cmfieldm
         module procedure ipcmfieldm2
         module procedure ipdmfieldm2
      end interface
!
      interface emfieldm
         module procedure ipemfieldm2
      end interface
!
      interface cpfieldm
         module procedure ipcpfieldm2
      end interface
!
      interface avpotm
         module procedure ipavpotm23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipsglsin2c(cu,cu1,nx,kyp,inorder)
! double array in x dimension for 2d vector data
! for mixed periodic/dirichlet boundary conditions
         implicit none
         integer :: nx, kyp
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: cu
         real, dimension(:,:,:,:), pointer :: cu1
! local data
         integer :: nxv, nypmx, kypd, kblok, order
         nxv = size(cu,2);  nypmx = size(cu,3); kblok = size(cu,4)
         kypd = size(cu1,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(cu,1)==3) then
            if (order==LINEAR) then
               call PSGLSIN2B(cu(1,1,1,1),cu1,nx,kyp,nxv,nypmx,kypd,kblo&
     &k)
            else
               call PSGLSIN2B(cu(1,2,2,1),cu1,nx,kyp,nxv,nypmx,kypd,kblo&
     &k)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call PSGLSIN2C(cu(1,1,1,1),cu1,nx,kyp,nxv,nypmx,kypd,kblo&
     &k)
            else
               call PSGLSIN2C(cu(1,2,2,1),cu1,nx,kyp,nxv,nypmx,kypd,kblo&
     &k)
            endif
         endif
         end subroutine ipsglsin2c
!
         subroutine ipsglsin2d(q,q1,nx,kyp,inorder)
! double array in x dimension for 2d scalar data
! for mixed periodic/dirichlet boundary conditions
         implicit none
         integer :: nx, kyp
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         real, dimension(:,:,:), pointer :: q1
! local data
         integer :: nxv, nypmx, kypd, kblok, order
         nxv = size(q,1);  nypmx = size(q,2); kblok = size(q,3)
         kypd = size(q1,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PSGLSIN2D(q(1,1,1),q1,nx,kyp,nxv,nypmx,kypd,kblok)
         else
            call PSGLSIN2D(q(2,2,1),q1,nx,kyp,nxv,nypmx,kypd,kblok)
         endif
         end subroutine ipsglsin2d
!
         subroutine iphafsgl2c(fxy,fxy1,nx,kyp,inorder)
! copy from double to normal array in x dimension for 2d vector data
         implicit none
         integer :: nx, kyp
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:,:), pointer :: fxy1
! local data
         integer :: nxv, nypmx, kypd, kblok, order
         nxv = size(fxy,2);  nypmx = size(fxy,3); kblok = size(fxy,4)
         kypd = size(fxy1,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(fxy,1)==2) then
            if (order==LINEAR) then
               call PHAFSGL2C(fxy(1,1,1,1),fxy1,nx,kyp,nxv,nypmx,kypd,kb&
     &lok)
            else
               call PHAFSGL2C(fxy(1,2,2,1),fxy1,nx,kyp,nxv,nypmx,kypd,kb&
     &lok)
            endif
         else if (size(fxy,1)==3) then
            if (order==LINEAR) then
               call PHAFSGL2B(fxy(1,1,1,1),fxy1,nx,kyp,nxv,nypmx,kypd,kb&
     &lok)
            else
               call PHAFSGL2B(fxy(1,2,2,1),fxy1,nx,kyp,nxv,nypmx,kypd,kb&
     &lok)
            endif
         endif
         end subroutine iphafsgl2c
!
         subroutine iphafsgl2d(q,q1,nx,kyp,inorder)
! copy from double to normal array in x dimension for 2d scalar data
         implicit none
         integer :: nx, kyp
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         real, dimension(:,:,:), pointer :: q1
! local data
         integer :: nxv, nypmx, kypd, kblok, order
         nxv = size(q,1);  nypmx = size(q,2); kblok = size(q,3)
         kypd = size(q1,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PHAFSGL2D(q(1,1,1),q1,nx,kyp,nxv,nypmx,kypd,kblok)
         else
            call PHAFSGL2D(q(2,2,1),q1,nx,kyp,nxv,nypmx,kypd,kblok)
         endif
         end subroutine iphafsgl2d
!
         subroutine ippoism2(q,fx,ffb,we,nx,ny,kstrt)
! poisson solver for 2d potential, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q, fx, ffb
! local data
         integer :: isign = 1, nyv, kxp2, j2blok, nyhd
         real :: ax, ay, affp
         complex, dimension(1,1,1) :: fy
         nyv = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         nyhd = size(ffb,1)
         call PPOISMX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,nyv,k&
     &xp2,j2blok,nyhd)
         end subroutine ippoism2
!
         subroutine ipspoism2(q,fy,ffb,nx,ny,kstrt)
! smoother for 2d scalar field, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: q, fy, ffb
! local data
         integer :: isign = 2, nyv, kxp2, j2blok, nyhd
         real :: ax, ay, affp, we
         complex, dimension(1,1,1) :: fx
         nyv = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         nyhd = size(ffb,1)
         call PPOISMX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,nyv,k&
     &xp2,j2blok,nyhd)
         end subroutine ipspoism2
!
         subroutine ippoism22init(ffb,ax,ay,affp,nx,ny,kstrt)
! initialize 2d electric field solver,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: ax, ay, affp
         complex, dimension(:,:,:), pointer :: ffb
! local data
         integer :: isign = 0, nyv, kxp2, j2blok, nyhd
         real :: we
         complex, dimension(1,1,1) :: q
         complex, dimension(2,1,1,1) :: fxy
         nyv = size(q,1)
         nyhd = size(ffb,1); kxp2 = size(ffb,2); j2blok = size(ffb,3)
         call PPOISMX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,nyv,kx&
     &p2,j2blok,nyhd)
         end subroutine ippoism22init
!
         subroutine ippoism22(q,fxy,ffb,we,nx,ny,kstrt)
! poisson solver for 2d electric field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q, ffb
         complex, dimension(:,:,:,:), pointer :: fxy
! local data
         integer :: isign = -1, nyv, kxp2, j2blok, nyhd
         real :: ax, ay, affp
         nyv = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         nyhd = size(ffb,1)
         if (size(fxy,1)==2) then
            call PPOISMX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,nyv&
     &,kxp2,j2blok,nyhd)
         else if (size(fxy,1)==3) then
            call PPOISMX23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,nyv&
     &,kxp2,j2blok,nyhd)
         endif
         end subroutine ippoism22
!
         subroutine ippoism23(q,fxy,ffb,we,nx,ny,kstrt)
! poisson solver for 2-1/2d electric field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         real, dimension(:,:,:), pointer :: q
         real, dimension(:,:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: ffb
! local data
         integer :: isign = -1, nyvh, kxp2, j2blok, nyhd
         real :: ax, ay, affp
         nyvh = size(q,1)/2; kxp2 = size(q,2); j2blok = size(q,3)
         nyhd = size(ffb,1)
         call PPOISM23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,nyvh,kx&
     &p2,j2blok,nyhd)
         end subroutine ippoism23
!
         subroutine ipcuperpm2(cu,nx,ny,kstrt)
! calculates transverse part of 2d vector field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         real, dimension(:,:,:,:), pointer :: cu
! local data
         integer :: nyvh, kxp2, j2blok
         nyvh = size(cu,2)/2; kxp2 = size(cu,3); j2blok = size(cu,4)
         call PCUPERPM2(cu,nx,ny,kstrt,nyvh,kxp2,j2blok)
         end subroutine ipcuperpm2
!
         subroutine jpbpoism23(cu,bxy,ffb,ci,wm,nx,ny,kstrt)
! caculates static vector potential for 2d vector field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffb
! local data
         integer :: isign = 1, nyv, kxp2, j2blok, nyhd
         real :: ax, ay, affp
         nyv = size(cu,2); kxp2 = size(cu,3); j2blok = size(cu,4)
         nyhd = size(ffb,1)
         call PBPOISMX23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,kstrt,n&
     &yv,kxp2,j2blok,nyhd)
         end subroutine jpbpoism23
!
         subroutine jipbpoism23(cu,bxy,ffb,ci,wm,nx,ny,kstrt)
! calculates static magnetic field for periodic 2d vector field
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         real, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffb
! local data
         integer :: nyvh, kxp2, j2blok, nyhd
         nyvh = size(cu,2)/2; kxp2 = size(cu,3); j2blok = size(cu,4)
         nyhd = size(ffb,1)
         call IPBPOISM23(cu,bxy,ffb,ci,wm,nx,ny,kstrt,nyvh,kxp2,j2blok,n&
     &yhd)
         end subroutine jipbpoism23
!
         subroutine ipmaxwelm2(exy,bxy,cu,ffb,affp,ci,dt,wf,wm,nx,ny,kst&
     &rt)
! calculates maxwell's equation for 2d vector field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: affp, ci, dt, wf, wm
         real, dimension(:,:,:,:), pointer :: exy, bxy, cu
         complex, dimension(:,:,:), pointer :: ffb
! local data
         integer :: nyvh, kxp2, j2blok, nyhd
         nyvh = size(cu,2)/2; kxp2 = size(cu,3); j2blok = size(cu,4)
         nyhd = size(ffb,1)
         call PMAXWELM2(exy,bxy,cu,ffb,affp,ci,dt,wf,wm,nx,ny,kstrt,nyvh&
     &,kxp2,j2blok,nyhd)
         end subroutine ipmaxwelm2
!
         subroutine ipcmfieldm2(cu1,cu,nx,ny,kstrt)
! copies from double to normal array in y dimension for 2d vector data
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: cu1
         real, dimension(:,:,:,:), pointer :: cu
! local data
         integer :: nyv, kxp2, j2blok, nyvh
         nyv = size(cu1,2); kxp2 = size(cu1,3); j2blok = size(cu1,4)
         nyvh = size(cu,2)/2
         call PCMFIELDM2(cu1,cu,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
         end subroutine ipcmfieldm2
!
         subroutine ipdmfieldm2(q1,q,nx,ny,kstrt)
! copies from double to normal array in y dimension for 2d scalar data
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: q1
         real, dimension(:,:,:), pointer :: q
! local data
         integer :: nyv, kxp2, j2blok, nyvh
         nyv = size(q1,1); kxp2 = size(q1,2); j2blok = size(q1,3)
         nyvh = size(q,1)/2
         call PDMFIELDM2(q1,q,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
         end subroutine ipdmfieldm2
!
         subroutine ipemfieldm2(fxy,exy,ffb,isign,nx,ny,kstrt)
! combines and smooths 2d vector fields,
! mixed conducting/periodic boundaries
         implicit none
         integer :: isign, nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:,:), pointer :: exy
         complex, dimension(:,:,:), pointer :: ffb
! local data
         integer :: nyv, kxp2, j2blok, nyvh, nyhd
         nyv = size(fxy,2); kxp2 = size(fxy,3); j2blok = size(fxy,4)
         nyvh = size(exy,2)/2; nyhd = size(ffb,1)
         call PEMFIELDM2(fxy,exy,ffb,isign,nx,ny,kstrt,nyv,nyvh,kxp2,j2b&
     &lok,nyhd)
         end subroutine ipemfieldm2
!
         subroutine ipcpfieldm2(fxy,exy,nx,ny,kstrt)
! combines 2d electric fields, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:,:), pointer :: exy
! local data
         integer :: nyv, kxp2, j2blok, nyvh
         nyv = size(fxy,2); kxp2 = size(fxy,3); j2blok = size(fxy,4)
         nyvh = size(exy,2)/2
         call PCPFIELDM2(fxy,exy,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
         end subroutine ipcpfieldm2
!
         subroutine ipavpotm23(bxy,axy,nx,ny,kstrt)
! calculates 2d vector potential from magnetic field
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny, kstrt
         real, dimension(:,:,:,:), pointer :: bxy, axy
! local data
         integer :: nyvh, kxp2, j2blok
         nyvh = size(bxy,2)/2; kxp2 = size(bxy,3); j2blok = size(bxy,4)
         call PAVPOTM23(bxy,axy,nx,ny,kstrt,nyvh,kxp2,j2blok)
         end subroutine ipavpotm23
!
      end module pbfield2d
