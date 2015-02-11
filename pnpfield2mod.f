!-----------------------------------------------------------------------
!
      module pnpfield2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 libraries
! pfield2lib.f, pdfield2lib.f, pbfield2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 15, 2009
!
      use globals, only: LINEAR, QUADRATIC
      use pfield2d
      use pdfield2d
      use pbfield2d
      use pcfield2d
      use pnfield2d
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: cguard, bguard, sguard, aguard, zguard
      public :: pois_init, pois, spois, pois3, cuperp, bpois, sbpois
      public :: ibpois, maxwel, emfield, emfieldr, apois, avpot, avrpot
      public :: gtmodes, ptmodes, poynt
      public :: amcguard, dcuperp, adcuperp, epois_init, epois, iepois
      public :: addqei, baddext
      public :: imoment, ipdivf2, ipgradf2, ipcurlf2, ipcurlf22
      public :: poisd_init, poisdx, poisd, cuperpd, bpoisdx, bpoisd
      public :: ibpoisd, maxweld, cmfieldd, emfieldd, cpfieldd, avpotd
      public :: ipdivfd2, ipgradfd2, ipcurlfd2, ipcurlfd22
      public :: sglsin, hafsgl, poism_init, poism, poism3, cuperpm
      public :: bpoism, ibpoism, maxwelm, cmfieldm, emfieldm, cpfieldm
      public :: avpotm
      public :: sguardp, aguardp, cguardp
      public :: poisc2_init, poisc3_init, poisc
      public :: poisn_init, poisn, cuperpn, bpoisn
      public :: ibpoisn, maxweln, cmfieldn, emfieldn, cpfieldn, avpotn
      public :: addfields
!
! define generic interfaces to Fortran90 library
!
      interface sguardp
         module procedure ipsguard2xp
         module procedure ipscguard2xp
      end interface
! 
      interface aguardp
         module procedure ipaguard2xp
         module procedure ipacguard2xp
      end interface
!
      interface cguardp
         module procedure ipcguard2xp
         module procedure ipdguard2xp
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipscguard2xp(cu,kstrt,nvp,noff,nyp,xj0,yj0,zj0,nx,ny&
     &,ipbc,inorder)
! initialize non-uniform 2d vector field
         implicit none
         integer :: kstrt, nvp, nx, ny, ipbc
         integer, optional :: inorder
         real :: xj0, yj0, zj0
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: noff, nyp
! local data
         integer :: ngx = 1, ngy = 1, nxe, nypmx, nblok, order
         nxe = size(cu,2); nypmx = size(cu,3); nblok = size(cu,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(cu,1)==3) then
            if (ipbc==1) then
               if (order==LINEAR) then
                  call PSCGUARD2L(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
               else
                  call PSCGUARD2(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
               endif
            else if (ipbc==2) then
               if (order==LINEAR) then
                  call PLSCGUARD2XL(cu,kstrt,nvp,noff,nyp,xj0,yj0,zj0,nx&
     &,ny,ngx,ngy,nxe,nypmx,nblok)
               else
                  call PLSCGUARD2X(cu,kstrt,nvp,noff,nyp,xj0,yj0,zj0,nx,&
     &ny,ngx,ngy,nxe,nypmx,nblok)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call PMSCGUARD2L(cu,nyp,xj0,yj0,zj0,nx,ngx,nxe,nypmx,n&
     &blok)
               else
                  call PMSCGUARD2(cu,nyp,xj0,yj0,zj0,nx,ngx,nxe,nypmx,nb&
     &lok)
               endif
            endif
         else if (size(cu,1)==2) then
            if (ipbc==1) then
               if (order==LINEAR) then
                  call PSCGUARD22L(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
               else
                  call PSCGUARD22(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
               endif
            else if (ipbc==2) then
               if (order==LINEAR) then
                  call PLSCGUARD22XL(cu,kstrt,nvp,noff,nyp,xj0,yj0,nx,ny&
     &,ngx,ngy,nxe,nypmx,nblok)
               else
                  call PLSCGUARD22X(cu,kstrt,nvp,noff,nyp,xj0,yj0,nx,ny,&
     &ngx,ngy,nxe,nypmx,nblok)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call PMSCGUARD22L(cu,nyp,xj0,yj0,nx,ngx,nxe,nypmx,nblo&
     &k)
               else
                  call PMSCGUARD22(cu,nyp,xj0,yj0,nx,ngx,nxe,nypmx,nblok&
     &)
               endif
            endif
         endif
         end subroutine ipscguard2xp
!
         subroutine ipsguard2xp(q,kstrt,nvp,noff,nyp,qi0,nx,ny,ipbc,inor&
     &der)
! initialize non-uniform 2d scalar field
         implicit none
         integer :: kstrt, nvp, nx, ny, ipbc
         integer, optional :: inorder
         real :: qi0
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: noff, nyp
! local data
         integer :: ngx = 1, ngy = 1, nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PSGUARD2L(q,nyp,qi0,nx,nxe,nypmx,nblok)
            else
               call PSGUARD2(q,nyp,qi0,nx,nxe,nypmx,nblok)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               call PLSGUARD2XL(q,kstrt,nvp,noff,nyp,qi0,nx,ny,ngx,ngy,n&
     &xe,nypmx,nblok)
            else
               call PLSGUARD2X(q,kstrt,nvp,noff,nyp,qi0,nx,ny,ngx,ngy,nx&
     &e,nypmx,nblok)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call PMSGUARD2L(q,nyp,qi0,nx,ngx,nxe,nypmx,nblok)
            else
               call PMSGUARD2(q,nyp,qi0,nx,ngx,nxe,nypmx,nblok)
            endif
         endif
         end subroutine ipsguard2xp
!
         subroutine ipacguard2xp(cu,nyp,nx,ipbc,inorder)
! add guard cells in x for non-uniform 2d vector data
         implicit none
         integer :: nx, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(cu,2); nypmx = size(cu,3); nblok = size(cu,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(cu,1)==3) then
            if (ipbc==1) then
               if (order==LINEAR) then
                  call PACGUARD2XL(cu,nyp,nx,nxe,nypmx,nblok)
               else
                  call PACGUARD2X(cu,nyp,nx,nxe,nypmx,nblok)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call PLACGUARD2X(cu(1,2,1,1),nyp,nx-2,nxe,nypmx,nblok)
               endif
            else if (ipbc==3) then
               if (order==QUADRATIC) then
                  call PLACGUARD2X(cu(1,2,1,1),nyp,nx-2,nxe,nypmx,nblok)
               endif
            endif
         else if (size(cu,1)==2) then
            if (ipbc==1) then
               if (order==LINEAR) then
                  call PACGUARD22XL(cu,nyp,nx,nxe,nypmx,nblok)
               else
                  call PACGUARD22X(cu,nyp,nx,nxe,nypmx,nblok)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call PLACGUARD22X(cu(1,2,1,1),nyp,nx-2,nxe,nypmx,nblok&
     &)
               endif
            else if (ipbc==3) then
               if (order==QUADRATIC) then
                  call PLACGUARD22X(cu(1,2,1,1),nyp,nx-2,nxe,nypmx,nblok&
     &)
               endif
            endif
         endif
         end subroutine ipacguard2xp
!
         subroutine ipaguard2xp(q,nyp,nx,ipbc,inorder)
! add guard cells in x for non-uniform 2d scalar data
         implicit none
         integer :: nx, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PAGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
            else
               call PAGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
            endif
         else if (ipbc==2) then
            if (order==QUADRATIC) then
               call PLAGUARD2X(q(2,1,1),nyp,nx-2,nxe,nypmx,nblok)
            endif
         else if (ipbc==3) then
            if (order==QUADRATIC) then
               call PLAGUARD2X(q(2,1,1),nyp,nx-2,nxe,nypmx,nblok)
            endif
         endif
         end subroutine ipaguard2xp
!
         subroutine ipcguard2xp(fxy,nyp,nx,ipbc,inorder)
! copy guard cells in x for non-uniform 2d vector data
         implicit none
         integer :: nx, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: fxy
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(fxy,2); nypmx = size(fxy,3); nblok = size(fxy,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(fxy,1)==2) then
            if (ipbc==1) then
               if (order==LINEAR) then
                  call PCGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
               else
                  call PCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call PLCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
               endif
            else if (ipbc==3) then
               if (order==QUADRATIC) then
                  call PLCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
               endif
            endif
         else if (size(fxy,1)==3) then
            if (ipbc==1) then
               if (order==LINEAR) then
                  call PBGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
               else
                  call PBGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call PLBGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
               endif
            else if (ipbc==3) then
               if (order==QUADRATIC) then
                  call PLBGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
               endif
            endif
         endif
         end subroutine ipcguard2xp
!
         subroutine ipdguard2xp(q,nyp,nx,ipbc,inorder)
! copy guard cells in x for non-uniform 2d scalar data
         implicit none
         integer :: nx, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call PDGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
            else
               call PDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
            endif
         else if (ipbc==2) then
            if (order==QUADRATIC) then
               call PLDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
            endif
         else if (ipbc==3) then
            if (order==QUADRATIC) then
               call PLDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
            endif
         endif
         end subroutine ipdguard2xp
!
      end module pnpfield2d
