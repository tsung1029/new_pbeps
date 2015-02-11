!-----------------------------------------------------------------------
!
      module pcfield2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pcfield2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: march 15, 2006
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: poisc2_init, poisc3_init, poisc
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PFORMC2(ffg,f,ft,bs,br,fpotc,mixup2,sct2,affp,ar,ind&
     &x1,indy1,kstrt,nxv,ny2d,kxp2,kyp2,j2blok,k2blok,kxp2d,ny1d,nxhy2,n&
     &xyh2)
         implicit none
         integer :: indx1, indy1, kstrt, nxv, ny2d, kxp2, kyp2
         integer :: j2blok, k2blok, kxp2d, ny1d, nxhy2, nxyh2
         real :: ar, affp
         real, dimension(4,ny1d,kxp2d,j2blok) :: ffg
         real, dimension(2*nxv,kyp2,k2blok) :: f
         complex, dimension(ny2d,kxp2,j2blok) :: ft
         complex, dimension(kxp2,kyp2,k2blok) :: bs
         complex, dimension(kxp2,kyp2,j2blok) :: br
         integer, dimension(nxhy2) :: mixup2
         complex, dimension(nxyh2) :: sct2
         real, external :: fpotc
         end subroutine
      end interface
      interface
         subroutine PPOISC2(q,fx,fy,isign,ffg,we,nx,ny,kstrt,ny2d,kxp2,j&
     &2blok,ny1d,kxp2d)
         implicit none
         real :: we
         integer :: isign, nx, ny, kstrt, ny2d, kxp2, j2blok
         integer :: ny1d, kxp2d 
         complex, dimension(ny2d,kxp2,j2blok) :: q, fx, fy
         real, dimension(4,ny1d,kxp2d,j2blok) :: ffg
         end subroutine
      end interface
      interface
         subroutine PPOISC22(q,fxy,ffg,we,nx,ny,kstrt,ny2d,kxp2,j2blok,n&
     &y1d,kxp2d)
         implicit none
         real :: we
         integer :: nx, ny, kstrt, ny2d, kxp2, j2blok, ny1d, kxp2d 
         complex, dimension(ny2d,kxp2,j2blok) :: q
         complex, dimension(2,ny2d,kxp2,j2blok) :: fxy
         real, dimension(4,ny1d,kxp2d,j2blok) :: ffg
         end subroutine
      end interface
      interface
         function POTC3(r,affp,ari,ifun)
         implicit none
         integer :: ifun
         real :: POTC3, r, affp, ari
         end function
      end interface
      interface
         function POTC2(r,affp,ari,ifun)
         implicit none
         integer :: ifun
         real :: POTC2, r, affp, ari
         end function
      end interface 
!
! define generic interfaces to Fortran90 library
!
      interface poisc2_init
         module procedure ippoisc2init
      end interface
!
      interface poisc3_init
         module procedure ippoisc3init
      end interface
!
      interface poisc
         module procedure ippoisc2
         module procedure ipspoisc2
         module procedure ippoisc22
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ippoisc2(q,fx,ffg,we,nx,ny,kstrt)
! poisson solver for 2d potential, open boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q, fx
         real, dimension(:,:,:,:), pointer :: ffg
! local data
         integer :: isign = 1, ny2d, kxp2, j2blok, ny1d, kxp2d
         complex, dimension(1,1,1) :: fy
         ny2d = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         ny1d = size(ffg,2); kxp2d = size(ffg,3)
         call PPOISC2(q,fx,fy,isign,ffg,we,nx,ny,kstrt,ny2d,kxp2,j2blok,&
     &ny1d,kxp2d)
         end subroutine ippoisc2
!
         subroutine ipspoisc2(q,fy,ffg,nx,ny,kstrt)
! smoother for 2d scalar field, open boundaries
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: q, fy
         real, dimension(:,:,:,:), pointer :: ffg
! local data
         integer :: isign = 2, ny2d, kxp2, j2blok, ny1d, kxp2d
         real :: we
         complex, dimension(1,1,1) :: fx
         ny2d = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         ny1d = size(ffg,2); kxp2d = size(ffg,3)
         call PPOISC2(q,fx,fy,isign,ffg,we,nx,ny,kstrt,ny2d,kxp2,j2blok,&
     &ny1d,kxp2d)
         end subroutine ipspoisc2
!
         subroutine ippoisc2init(ffg,f,ft,mixup2,sct2,ar,affp,indx,indy,&
     &kstrt)
! initialize 2d poisson solver, open boundary conditions
         implicit none
         integer :: indx, indy, kstrt
         real :: ar, affp
         real, dimension(:,:,:,:), pointer :: ffg
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: ft
         integer, dimension(:), pointer :: mixup2
         complex, dimension(:), pointer :: sct2
! local data
         integer :: indx1, indy1, nxv, ny2d, kxp2, kyp2, j2blok, k2blok
         integer :: kxp2d, ny1d, nxhy2, nxyh2
         complex, dimension(size(ft,2),size(f,2),size(f,3)) :: bs
         complex, dimension(size(ft,2),size(f,2),size(ft,3)) :: br
         real, external :: POTC2
         indx1 = indx + 1; indy1 = indy + 1
         ny1d = size(ffg,2); kxp2d = size(ffg,3)
         nxv = size(f,1)/2; kyp2 = size(f,2); k2blok = size(f,3)
         ny2d = size(ft,1); kxp2 = size(ft,2); j2blok = size(ft,3)
         nxhy2 = size(mixup2); nxyh2 = size(sct2)
         call PFORMC2(ffg,f,ft,bs,br,POTC2,mixup2,sct2,affp,ar,indx1,ind&
     &y1,kstrt,nxv,ny2d,kxp2,kyp2,j2blok,k2blok,kxp2d,ny1d,nxhy2,nxyh2)
         end subroutine ippoisc2init
!
         subroutine ippoisc3init(ffg,f,ft,mixup2,sct2,ar,affp,indx,indy,&
     &kstrt)
! initialize 2d poisson solver with 3d fields, open boundary conditions
         implicit none
         integer :: indx, indy, kstrt
         real :: ar, affp
         real, dimension(:,:,:,:), pointer :: ffg
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: ft
         integer, dimension(:), pointer :: mixup2
         complex, dimension(:), pointer :: sct2
! local data
         integer :: indx1, indy1, nxv, ny2d, kxp2, kyp2, j2blok, k2blok
         integer :: kxp2d, ny1d, nxhy2, nxyh2
         complex, dimension(size(ft,2),size(f,2),size(f,3)) :: bs
         complex, dimension(size(ft,2),size(f,2),size(ft,3)) :: br
         real, external :: POTC3
         indx1 = indx + 1; indy1 = indy + 1
         ny1d = size(ffg,2); kxp2d = size(ffg,3)
         nxv = size(f,1)/2; kyp2 = size(f,2); k2blok = size(f,3)
         ny2d = size(ft,1); kxp2 = size(ft,2); j2blok = size(ft,3)
         nxhy2 = size(mixup2); nxyh2 = size(sct2)
         call PFORMC2(ffg,f,ft,bs,br,POTC3,mixup2,sct2,affp,ar,indx1,ind&
     &y1,kstrt,nxv,ny2d,kxp2,kyp2,j2blok,k2blok,kxp2d,ny1d,nxhy2,nxyh2)
         end subroutine ippoisc3init
!
         subroutine ippoisc22(q,fxy,ffg,we,nx,ny,kstrt)
! poisson solver for 2d electric field, open boundaries
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q
         complex, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:,:), pointer :: ffg
! local data
         integer :: ny2d, kxp2, j2blok, ny1d, kxp2d
         ny2d = size(q,1); kxp2 = size(q,2); j2blok = size(q,3)
         ny1d = size(ffg,2); kxp2d = size(ffg,3)
         if (size(fxy,1)==2) then
            call PPOISC22(q,fxy,ffg,we,nx,ny,kstrt,ny2d,kxp2,j2blok,ny1d&
     &,kxp2d)
!        else if (size(fxy,1)==3) then
!           call PPOISC23(q,fxy,ffg,we,nx,ny,kstrt,ny2d,kxp2,j2blok,ny1d&
!    &,kxp2d)
         endif
         end subroutine ippoisc22
!
      end module pcfield2d
