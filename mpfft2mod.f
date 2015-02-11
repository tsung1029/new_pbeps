!-----------------------------------------------------------------------
!
      module mpfft2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library mpfft2lib.f
! mpfft2mod.f contains multi-tasking interface procedures to perform
!             ffts:
!             defines module mpfft2d
! fft => impfft2r performs 2d scalar real to complex fft and its inverse
!        calls MPFFT2R
! fft => impfft2r2 performs 2d vector complex to real fft for 2
!        component vectors.
!        calls MPFFT2R2
! fft => impfft2r3 performs 2d vector real to complex fft for 1, 2, or 3
!        component vectors and their inverses.
!        calls MPFFT2R, MPFFT2R2, or MPFFT2R3
! fftn => impfft2rn perform 2d vector real to complex fft for n
!         component vectors and their inverses.
!         calls MPFFT2RN
! fftn => imp1fft2rn performs 2d real to complex fft for a scalar and n
!         component vector and their inverses.
!         calls MP2FFT2RN
! fftn => imp2fft2rn perform 2d real to complex fft for two n component
!         vectors and their inverses.
!         calls MP2FFT2RN
! fsst => impfsst2r performs 2d scalar real sine-sine transform.
!         calls MPFSST2R
! fcct => impfcct2r performs 2d scalar real cosine-cosine transform.
!         calls MPFCCT2R
! fcst => impfcst2r3 performs 2d vector real cosine-sine transforms
!         transforms for 2 or 3 component vector arrays.
!         calls  MPFCST2R2, or MPFCST2R3
! fsct => impfsct2r3 perform 2d vector real sine-cosine transforms
!         for 2 or 3 component vector arrays.
!         calls MPFSCT2R2, or MPFSCT2R3
! fsft => impfsft2r performs 2d scalar real sine/periodic transform.
!         calls MPFSFT2R
! fcft => impfcft2r performs 2d scalar real cosine/periodic transform.
!         calls MPFCFT2R
! fcsft => impfcsft2r3 performs 2d vector real cosine-sine/periodic
!          transforms for 2 or 3 component vector arrays.
!          calls  MPFCSFT2R2, or MPFCSFT2R3
! fscft => impfscft2r3 performs 2d vector real sine-cosine/periodic
!          transforms for 2 or 3 component vector arrays.
!          calls  MPFSCFT2R2, or MPFSCFT2R3
! fdsft => impfdsft2rx performs 2d scalar real half sine/periodic
!          transform.
!          calls MPFDSFT2RX
! fdcft => impfdcft2rx performs 2d scalar real half cosine/periodic
!          transform.
!          calls MPFDCFT2RX
! fdcsft => impfdcsft2r3 performs 2d vector real half
!           cosine-sine/periodic transforms for 2 or 3 component vector
!           arrays.
!           calls MPFDCSFT2R2, or MPFDCSFT2R3
! fdscft => impfdscft2r3 performs 2d vector real half
!           sine-cosine/periodic transforms for 2 or 3 component vector
!           arrays.
!           calls MPFDSCFT2R2, or MPFDSCFT2R3
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 12, 2009
!
      use globals, only: LINEAR, QUADRATIC
      use pfft2d, only: wtimer, fft_init, fftc_init, fst_init, fdt_init
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: fft_init, fft, fftn, fftc_init
      public :: fst_init, fsst, fcct, fcst, fsct
      public :: fsft, fcft, fcsft, fscft
      public :: fdt_init, fdsft, fdcft, fdcsft, fdscft
!
! buffer data for transpose
      complex, dimension(:,:,:,:), allocatable :: bs, br
      integer :: szbuf = 0
      save
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MPFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,in&
     &dy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftas&
     &k,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(nyv,kxp,jblok) :: g
         complex, dimension(kxp,kyp,kblok) :: bs
         complex, dimension(kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT2R2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,i&
     &ndy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,ifta&
     &sk,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,nyv,kxp,jblok) :: g
         complex, dimension(2,kxp,kyp,kblok) :: bs
         complex, dimension(2,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT2R3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,i&
     &ndy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,ifta&
     &sk,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,nyv,kxp,jblok) :: g
         complex, dimension(3,kxp,kyp,kblok) :: bs
         complex, dimension(3,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT2RX(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,ks&
     &trt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask,nmt&
     &,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT2RX2(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,k&
     &strt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask,nm&
     &t,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT2RX3(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,k&
     &strt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask,nm&
     &t,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT2RN(f,g,bs,br,ss,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd,kx&
     &yip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, ndim, nxhyd, nxyhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(ndim,nyv,kxp,jblok) :: g
         complex, dimension(ndim,kxp,kyp,kblok) :: bs
         complex, dimension(ndim,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         complex, dimension(ndim,nxvh,nmt+1) :: ss
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFFT2RXN(f,g,ss,isign,ntpose,mixup,sct,ttp,indx,ind&
     &y,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd,kxyip,i&
     &ftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, ndim, nxhyd, nxyhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(ndim,nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         complex, dimension(ndim,nxvh,nmt+1) :: ss
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MP2FFT2RN(f1,f2,g1,g2,bs,br,ss,isign,ntpose,mixup,sc&
     &t,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim1,ndim&
     &2,nxhyd,nxyhd,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, ndim1, ndim2
         integer :: nxhyd, nxyhd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f1, f2
         real :: f1, f2
         complex, dimension(ndim1,nyv,kxp,jblok) :: g1
         complex, dimension(ndim2,nyv,kxp,jblok) :: g2
         complex, dimension(ndim1+ndim2,kxp,kyp,kblok) :: bs
         complex, dimension(ndim1+ndim2,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         complex, dimension(ndim1+ndim2,nxvh,nmt+1) :: ss
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFSST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxy&
     &ip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         real, dimension(nyv,kxp2d,jblok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFSCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxy&
     &ip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         real, dimension(nyv,kxp2d,jblok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFCST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxy&
     &ip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         real, dimension(nyv,kxp2d,jblok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFCCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxy&
     &ip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         real, dimension(nyv,kxp2d,jblok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFCST2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kx&
     &yip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         real, dimension(2,nyv,kxp2d,jblok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFSCT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kx&
     &yip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         real, dimension(2,nyv,kxp2d,jblok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFCST2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kx&
     &yip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         real, dimension(3,nyv,kxp2d,jblok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFSCT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kx&
     &yip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         real, dimension(3,nyv,kxp2d,jblok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFSFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,k&
     &xyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(nyvh,kxp2d,j2blok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFCFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,k&
     &xyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(nyvh,kxp2d,j2blok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFCSFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,ind&
     &x,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd&
     &,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,nyvh,kxp2d,j2blok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFSCFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,ind&
     &x,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd&
     &,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,nyvh,kxp2d,j2blok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFCSFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,ind&
     &x,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd&
     &,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,nyvh,kxp2d,j2blok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFSCFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,ind&
     &x,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd&
     &,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,nyvh,kxp2d,j2blok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFDSFT2RX(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,t&
     &tp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(nyvh,kxp2d,j2blok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFDCFT2RX(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,t&
     &tp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(nyvh,kxp2d,j2blok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFDCSFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,&
     &ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxh&
     &yd,nxyd,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,nyvh,kxp2d,j2blok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFDSCFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,&
     &ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxh&
     &yd,nxyd,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,nyvh,kxp2d,j2blok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFDCSFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,&
     &ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxh&
     &yd,nxyd,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,nyvh,kxp2d,j2blok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MPFDSCFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,&
     &ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxh&
     &yd,nxyd,kxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         integer :: nmt, ierr
         real :: ttp
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,nyvh,kxp2d,j2blok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         integer, dimension(nmt) :: kxyip, iftask
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface fft
         module procedure impfft2r
         module procedure impfft2r2
         module procedure impfft2r3
!        module procedure impfft2rx
!        module procedure impfft2rx2
!        module procedure impfft2rx3
      end interface
!
      interface fftn
         module procedure impfft2rn
!        module procedure impfft2rxn
         module procedure imp1fft2rn
         module procedure imp2fft2rn
      end interface
!
      interface fsst
         module procedure impfsst2r
      end interface
!
      interface fcct
         module procedure impfcct2r
      end interface
!
      interface fcst
         module procedure impfcst2r3
      end interface
!
      interface fsct
         module procedure impfsct2r3
      end interface
!
      interface fsft
         module procedure impfsft2r
      end interface
!
      interface fcft
         module procedure impfcft2r
      end interface
!
      interface fcsft
         module procedure impfcsft2r3
      end interface
!
      interface fscft
         module procedure impfscft2r3
      end interface
!
      interface fdsft
         module procedure impfdsft2rx
      end interface
!
      interface fdcft
         module procedure impfdcft2rx
      end interface
!
      interface fdcsft
         module procedure impfdcsft2r3
      end interface
!
      interface fdscft
         module procedure impfdscft2r3
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine impfft2r(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,ky&
     &p,inorder)
! perform multi-tasking 2d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv, nmt
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order, ierr
!        complex, dimension(size(g,2),kyp,size(f,3)) :: bs
!        complex, dimension(size(g,2),kyp,size(g,3)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! check if size of buffers has changed
         if (szbuf < kxp*kyp) then
            if (szbuf /= 0) deallocate(bs,br)
! allocate buffers
            allocate(bs(1,kxp,kyp,kblok),br(1,kxp,kyp,jblok))
            szbuf = kxp*kyp
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,i&
     &ftask,nmt,ierr)
         else
            call MPFFT2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,i&
     &ftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft2r
!
         subroutine impfft2r2(f,g,mixup,sct,tfft,indx,indy,kstrt,kyp,ino&
     &rder)
! perform multi-tasking 2d vector real to complex fft
! for 2 component vectors
         implicit none
         integer :: indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 1, ntpose = 1, nxvh, nyv, nmt
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order, ierr
!        complex, dimension(2,size(g,3),kyp,size(f,4)) :: bs
!        complex, dimension(2,size(g,3),kyp,size(g,4)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! check if size of buffers has changed
         if (szbuf < 2*kxp*kyp) then
            if (szbuf /= 0) deallocate(bs,br)
! allocate buffers
            allocate(bs(2,kxp,kyp,kblok),br(2,kxp,kyp,jblok))
            szbuf = 2*kxp*kyp
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyi&
     &p,iftask,nmt,ierr)
         else
            call MPFFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyi&
     &p,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft2r2
!
         subroutine impfft2r3(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,k&
     &yp,inorder)
! perform multi-tasking 2d vector real to complex fft
! for 3 component vectors
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv, nmt
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order, ierr
!        complex, dimension(3,size(g,3),kyp,size(f,4)) :: bs
!        complex, dimension(3,size(g,3),kyp,size(g,4)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! check if size of buffers has changed
         if (szbuf < 3*kxp*kyp) then
            if (szbuf /= 0) deallocate(bs,br)
! allocate buffers
            allocate(bs(3,kxp,kyp,kblok),br(3,kxp,kyp,jblok))
            szbuf = 3*kxp*kyp
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         select case(size(f,1))
         case (1)
            if (order==LINEAR) then
               call MPFFT2R(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,tt&
     &p,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kx&
     &yip,iftask,nmt,ierr)
            else
               call MPFFT2R(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,tt&
     &p,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kx&
     &yip,iftask,nmt,ierr)
            endif
         case (2)
            if (order==LINEAR) then
               call MPFFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,k&
     &xyip,iftask,nmt,ierr)
            else
               call MPFFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,k&
     &xyip,iftask,nmt,ierr)
            endif
         case (3)
            if (order==LINEAR) then
               call MPFFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,k&
     &xyip,iftask,nmt,ierr)
            else
               call MPFFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,k&
     &xyip,iftask,nmt,ierr)
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft2r3
!
         subroutine impfft2rx(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,k&
     &yp,inorder)
! perform multi-tasking optimized 2d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv, kxp, kypd, nmt
         integer :: jblok, kblok, nxhyd, nxyhd, order, ierr
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT2RX(f(1,1,1),g,isign,ntpose,mixup,sct,ttp,indx,ind&
     &y,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask&
     &,nmt,ierr)
         else
            call MPFFT2RX(f(2,2,1),g,isign,ntpose,mixup,sct,ttp,indx,ind&
     &y,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask&
     &,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft2rx
!
         subroutine impfft2rx2(f,g,mixup,sct,tfft,indx,indy,kstrt,kyp,in&
     &order)
! perform multi-tasking optimized 2d vector real to complex fft
! for 2 component vector
         implicit none
         integer :: indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 1, ntpose = 1, nxvh, nyv, kxp, nmt
         integer :: kypd, jblok, kblok, nxhyd, nxyhd, order, ierr
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT2RX2(f(1,1,1,1),g,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,ift&
     &ask,nmt,ierr)
         else
            call MPFFT2RX2(f(1,2,2,1),g,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,ift&
     &ask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft2rx2
!
         subroutine impfft2rx3(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,&
     &kyp,inorder)
! perform multi-tasking optimized 2d vector real to complex fft
! for 3 component vector
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv, kxp, nmt
         integer :: kypd, jblok, kblok, nxhyd, nxyhd, order, ierr
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call MPFFT2RX2(f(1,1,1,1),g,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,&
     &iftask,nmt,ierr)
            else
               call MPFFT2RX2(f(1,2,2,1),g,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,&
     &iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call MPFFT2RX3(f(1,1,1,1),g,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,&
     &iftask,nmt,ierr)
            else
               call MPFFT2RX3(f(1,2,2,1),g,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,&
     &iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft2rx3
!
         subroutine impfft2rn(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,k&
     &yp,inorder)
! perform multi-tasking 2d vector real to complex fft
! for n component vectors
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, ndim, nxvh, nyv, nmt
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order, ierr
!        complex, dimension(size(f,1),size(g,3),kyp,size(f,4)) :: bs
!        complex, dimension(size(f,1),size(g,3),kyp,size(g,4)) :: br
         complex, dimension(size(f,1),size(f,2)/2,ntasks+1) :: ss
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! check if size of buffers has changed
         if (szbuf < ndim*kxp*kyp) then
            if (szbuf /= 0) deallocate(bs,br)
! allocate buffers
            allocate(bs(ndim,kxp,kyp,kblok),br(ndim,kxp,kyp,jblok))
            szbuf = ndim*kxp*kyp
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT2RN(f(1,1,1,1),g,bs,br,ss,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nx&
     &yhd,kxyip,iftask,nmt,ierr)
         else
            call MPFFT2RN(f(1,2,2,1),g,bs,br,ss,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nx&
     &yhd,kxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft2rn
!
         subroutine impfft2rxn(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,&
     &kyp,inorder)
! perform multi-tasking optimized 2d vector real to complex fft
! for n component vector
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, ndim, nxvh, nyv, kxp, nmt
         integer :: kypd, jblok, kblok, nxhyd, nxyhd, order, ierr
         complex, dimension(size(f,1),size(f,2)/2,ntasks+1) :: ss
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFFT2RXN(f(1,1,1,1),g,ss,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd,k&
     &xyip,iftask,nmt,ierr)
         else
            call MPFFT2RXN(f(1,2,2,1),g,ss,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd,k&
     &xyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfft2rxn
!
         subroutine imp1fft2rn(f1,f2,g1,g2,isign,mixup,sct,tfft,indx,ind&
     &y,kstrt,kyp,inorder)
! perform multi-tasking 2d real to complex fft for a scalar and n
! component vector
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f1
         real, dimension(:,:,:,:), pointer :: f2
         complex, dimension(:,:,:), pointer :: g1
         complex, dimension(:,:,:,:), pointer :: g2
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, ndim1, ndim2, ndim, nxvh, nyv, nmt
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order, ierr
!        complex, dimension(1+size(f2,1),size(g1,2),kyp,size(f1,3)) ::  &
!    &bs
!        complex, dimension(1+size(f2,1),size(g1,2),kyp,size(g1,3)) ::  &
!    &br
         complex, dimension(1+size(f2,1),size(f1,1)/2,ntasks+1):: ss
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         ndim1 = 1; ndim2 = size(f2,1); ndim = ndim1 + ndim2
         nxvh = size(f1,1)/2; kypd = size(f1,2); kblok = size(f1,3)
         nyv = size(g1,1); kxp = size(g1,2); jblok = size(g1,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! check if size of buffers has changed
         if (szbuf < ndim*kxp*kyp) then
            if (szbuf /= 0) deallocate(bs,br)
! allocate buffers
            allocate(bs(ndim,kxp,kyp,kblok),br(ndim,kxp,kyp,jblok))
            szbuf = ndim*kxp*kyp
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MP2FFT2RN(f1(1,1,1),f2(1,1,1,1),g1,g2,bs,br,ss,isign,nt&
     &pose,mixup,sct,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kbl&
     &ok,ndim1,ndim2,nxhyd,nxyhd,kxyip,iftask,nmt,ierr)
         else
            call MP2FFT2RN(f1(2,2,1),f2(1,2,2,1),g1,g2,bs,br,ss,isign,nt&
     &pose,mixup,sct,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kbl&
     &ok,ndim1,ndim2,nxhyd,nxyhd,kxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine imp1fft2rn
!
         subroutine imp2fft2rn(f1,f2,g1,g2,isign,mixup,sct,tfft,indx,ind&
     &y,kstrt,kyp,inorder)
! perform multi-tasking 2d real to complex fft for two n component
! vectors
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f1, f2
         complex, dimension(:,:,:,:), pointer :: g1, g2
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, ndim1, ndim2, ndim, nxvh, nyv, nmt
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order, ierr
!        complex, dimension(size(f1,1)+size(f2,1),size(g1,3),kyp,size(f1&
!    &,4)) :: bs
!        complex, dimension(size(f1,1)+size(f2,1),size(g1,3),kyp,size(g1&
!    &,4)) :: br
         complex, dimension(size(f1,1)+size(f2,1),size(f1,2)/2,ntasks+1)&
     &:: ss
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         ndim1 = size(f1,1); ndim2 = size(f2,1); ndim = ndim1 + ndim2
         nxvh = size(f1,2)/2; kypd = size(f1,3); kblok = size(f1,4)
         nyv = size(g1,2); kxp = size(g1,3); jblok = size(g1,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! check if size of buffers has changed
         if (szbuf < ndim*kxp*kyp) then
            if (szbuf /= 0) deallocate(bs,br)
! allocate buffers
            allocate(bs(ndim,kxp,kyp,kblok),br(ndim,kxp,kyp,jblok))
            szbuf = ndim*kxp*kyp
         endif
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MP2FFT2RN(f1(1,1,1,1),f2(1,1,1,1),g1,g2,bs,br,ss,isign,&
     &ntpose,mixup,sct,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,k&
     &blok,ndim1,ndim2,nxhyd,nxyhd,kxyip,iftask,nmt,ierr)
         else
            call MP2FFT2RN(f1(1,2,2,1),f2(1,2,2,1),g1,g2,bs,br,ss,isign,&
     &ntpose,mixup,sct,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,k&
     &blok,ndim1,ndim2,nxhyd,nxyhd,kxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine imp2fft2rn
!
         subroutine impfsst2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,&
     &kxp2,kyp,inorder)
! perform multi-tasking 2d scalar real sine-sine transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         real, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyv, kypd, kxp2d, nmt
         integer :: jblok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp2d = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFSST2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd&
     &,kxyip,iftask,nmt,ierr)
         else
            call MPFSST2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd&
     &,kxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfsst2r
!
         subroutine impfcct2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,&
     &kxp2,kyp,inorder)
! perform multi-tasking 2d scalar real sine-sine transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         real, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyv, kypd, kxp2d, nmt
         integer :: jblok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp2d = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFCCT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd&
     &,kxyip,iftask,nmt,ierr)
         else
            call MPFCCT2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd&
     &,kxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfcct2r
!
         subroutine impfcst2r3(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt&
     &,kxp2,kyp,inorder)
! perform multi-tasking 2d vector real cosine-sine transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         real, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyv, kypd, kxp2d, nmt
         integer :: jblok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp2d = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call MPFCST2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFCST2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call MPFCST2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFCST2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfcst2r3
!
         subroutine impfsct2r3(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt&
     &,kxp2,kyp,inorder)
! perform multi-tasking 2d vector real sine-cosine transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         real, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyv, kypd, kxp2d, nmt
         integer :: jblok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp2d = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call MPFSCT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFSCT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call MPFSCT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFSCT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd,kxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfsct2r3
!
         subroutine impfsft2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,&
     &kxp2,kyp,inorder)
! perform multi-tasking 2d scalar real sine/periodic transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d, nmt
         integer :: j2blok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyvh = size(g,1); kxp2d = size(g,2); j2blok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFSFT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nx&
     &yd,kxyip,iftask,nmt,ierr)
         else
            call MPFSFT2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nx&
     &yd,kxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfsft2r
!
         subroutine impfcft2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,&
     &kxp2,kyp,inorder)
! perform multi-tasking 2d scalar real cosine/periodic transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d, nmt
         integer :: j2blok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyvh = size(g,1); kxp2d = size(g,2); j2blok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFCFT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nx&
     &yd,kxyip,iftask,nmt,ierr)
         else
            call MPFCFT2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nx&
     &yd,kxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfcft2r
!
         subroutine impfcsft2r3(f,g,isign,mixup,sctd,tfft,indx,indy,kstr&
     &t,kxp2,kyp,inorder)
! perform multi-tasking 2d vector real cosine-sine/periodic transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d, nmt
         integer :: j2blok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyvh = size(g,2); kxp2d = size(g,3); j2blok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call MPFCSFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFCSFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd,kxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call MPFCSFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFCSFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd,kxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfcsft2r3
!
         subroutine impfscft2r3(f,g,isign,mixup,sctd,tfft,indx,indy,kstr&
     &t,kxp2,kyp,inorder)
! perform multi-tasking 2d vector realosine-cosine/periodic transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d, nmt
         integer :: j2blok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyvh = size(g,2); kxp2d = size(g,3); j2blok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call MPFSCFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFSCFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd,kxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call MPFSCFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFSCFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd,kxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfscft2r3
!
         subroutine impfdsft2rx(f,g,isign,mixup,sctd,sctdx,tfft,indx,ind&
     &y,kstrt,kxp2,kyp,inorder)
! perform multi-tasking 2d scalar real half sine/periodic transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d, nmt
         integer :: j2blok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyvh = size(g,1); kxp2d = size(g,2); j2blok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFDSFT2RX(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,sct&
     &dx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,&
     &nxhyd,nxyd,kxyip,iftask,nmt,ierr)
         else
            call MPFDSFT2RX(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,sct&
     &dx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,&
     &nxhyd,nxyd,kxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfdsft2rx
!
         subroutine impfdcft2rx(f,g,isign,mixup,sctd,sctdx,tfft,indx,ind&
     &y,kstrt,kxp2,kyp,inorder)
! perform multi-tasking 2d scalar real half cosine/periodic transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d, nmt
         integer :: j2blok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyvh = size(g,1); kxp2d = size(g,2); j2blok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call MPFDCFT2RX(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,sct&
     &dx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,&
     &nxhyd,nxyd,kxyip,iftask,nmt,ierr)
         else
            call MPFDCFT2RX(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,sct&
     &dx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,&
     &nxhyd,nxyd,kxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfdcft2rx
!
         subroutine impfdcsft2r3(f,g,isign,mixup,sctd,sctdx,tfft,indx,in&
     &dy,kstrt,kxp2,kyp,inorder)
! perform multi-tasking 2d vector real half cosine-sine/periodic
! transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d, nmt
         integer :: j2blok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyvh = size(g,2); kxp2d = size(g,3); j2blok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call MPFDCSFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFDCSFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd,kxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call MPFDCSFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFDCSFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd,kxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfdcsft2r3
!
         subroutine impfdscft2r3(f,g,isign,mixup,sctd,sctdx,tfft,indx,in&
     &dy,kstrt,kxp2,kyp,inorder)
! perform multi-tasking 2d vector real half sine-cosine/periodic
! transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d, nmt
         integer :: j2blok, kblok, nxhyd, nxyd, order, ierr
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         integer, dimension(ntasks) :: kxyip, iftask
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyvh = size(g,2); kxp2d = size(g,3); j2blok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call MPFDSCFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFDSCFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd,kxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call MPFDSCFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd,kxyip,iftask,nmt,ierr)
            else
               call MPFDSCFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd,kxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine impfdscft2r3
!
      end module mpfft2d

