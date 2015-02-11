!-----------------------------------------------------------------------
!
      module pfft2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pfft2lib.f
! pfft2mod.f contains interface procedures to perform ffts:
!            defines module pfft2d
! fft_init => iwpfft2rinit initializes 2d real to complex fft tables.
!             calls WPFFT2RINIT
! fft => iwpfft2r performs 2d scalar real to complex fft and its inverse
!        calls WPFFT2R
! fft => iwpfft2r2 performs 2d vector complex to real fft for 2
!        component vectors.
!        calls WPFFT2R2
! fft => iwpfft2r3 performs 2d vector real to complex fft for 1, 2, or 3
!        component vectors and their inverses.
!        calls WPFFT2R, WPFFT2R2, or WPFFT2R3
! fft => ipfft2c performs 2d scalar complex to complex fft and its
!        inverse.
!        calls PFFT2C
! fftn => iwpfft2rn perform 2d vector real to complex fft for n
!         component vectors and their inverses.
!         calls WPFFT2RN
! fftn => iwp1fft2rn performs 2d real to complex fft for a scalar and n
!         component vector and their inverses.
!         calls WP2FFT2RN
! fftn => iwp2fft2rn perform 2d real to complex fft for two n component
!         vectors and their inverses.
!         calls WP2FFT2RN
! fftc_init => ipfft2cinit initializes 2d complex to complex fft tables.
!              calls PFFT2C
! fst_init => iwfst2rinit initializes 2d real sine-cosine transform.
!             tables.
!             calls WPFST2RINIT
! fsst => iwpfsst2r performs 2d scalar real sine-sine transform.
!         calls WPFSST2R
! fcct => iwpfcct2r performs 2d scalar real cosine-cosine transform.
!         calls WPFCCT2R
! fcst => iwpfcst2r3 performs 2d vector real cosine-sine transforms
!         transforms for 2 or 3 component vector arrays.
!         calls  WPFCST2R2, or WPFCST2R3
! fsct => iwpfsct2r3 perform 2d vector real sine-cosine transforms
!         for 2 or 3 component vector arrays.
!         calls WPFSCT2R2, or WPFSCT2R3
! fsft => iwpfsft2r performs 2d scalar real sine/periodic transform.
!         calls WPFSFT2R
! fcft => iwpfcft2r performs 2d scalar real cosine/periodic transform.
!         calls WPFCFT2R
! fcsft => iwpfcsft2r3 performs 2d vector real cosine-sine/periodic
!          transforms for 2 or 3 component vector arrays.
!          calls  WPFCSFT2R2, or WPFCSFT2R3
! fscft => iwpfscft2r3 performs 2d vector real sine-cosine/periodic
!          transforms for 2 or 3 component vector arrays.
!          calls  WPFSCFT2R2, or WPFSCFT2R3
! fdt_init => iwpfdt2rinit initializes 2d real half sine-cosine/periodic
!             transform tables.
!             calls WFDT2RINIT
! fdsft => iwpfdsft2rx performs 2d scalar real half sine/periodic
!          transform.
!          calls WPFDSFT2RX
! fdcft => iwpfdcft2rx performs 2d scalar real half cosine/periodic
!          transform.
!          calls WPFDCFT2RX
! fdcsft => iwpfdcsft2r3 performs 2d vector real half
!           cosine-sine/periodic transforms for 2 or 3 component vector
!           arrays.
!           calls WPFDCSFT2R2, or WPFDCSFT2R3
! fdscft => iwpfdscft2r3 performs 2d vector real half
!           sine-cosine/periodic transforms for 2 or 3 component vector
!           arrays.
!           calls WPFDSCFT2R2, or WPFDSCFT2R3
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 12, 2009
!
      use globals, only: LINEAR, QUADRATIC
      use p0d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: wtimer
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
         subroutine PFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,ks&
     &trt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: f
         complex, dimension(nyv,kxp,jblok) :: g
         complex, dimension(kxp,kyp,kblok) :: bs
         complex, dimension(kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT2R2(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,k&
     &strt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: f
         complex, dimension(2,nyv,kxp,jblok) :: g
         complex, dimension(2,kxp,kyp,kblok) :: bs
         complex, dimension(2,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT2R3(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,k&
     &strt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: f
         complex, dimension(3,nyv,kxp,jblok) :: g
         complex, dimension(3,kxp,kyp,kblok) :: bs
         complex, dimension(3,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT2RX(f,g,isign,ntpose,mixup,sct,indx,indy,kstrt,n&
     &xvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: f
         complex, dimension(nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT2RX2(f,g,isign,ntpose,mixup,sct,indx,indy,kstrt,&
     &nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: f
         complex, dimension(2,nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT2RX3(f,g,isign,ntpose,mixup,sct,indx,indy,kstrt,&
     &nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: f
         complex, dimension(3,nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine PFFT2C(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,ks&
     &trt,nxv,nyv,kxp,kyp,kypd,jblok,kblok,nxyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxv, nyv, kxp, kyp
         integer :: kypd, jblok, kblok, nxyd, nxyhd
!        complex, dimension(*) :: f
         complex :: f
         complex, dimension(nyv,kxp,jblok) :: g
         complex, dimension (kxp,kyp,kblok) :: bs
         complex, dimension (kxp,kyp,jblok) :: br
         integer, dimension(nxyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer :: indx, indy
         integer :: nxhyd, nxyhd
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,in&
     &dy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: ttp
         real :: f
         complex, dimension(nyv,kxp,jblok) :: g
         complex, dimension(kxp,kyp,kblok) :: bs
         complex, dimension(kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT2RX(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,ks&
     &trt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: ttp
         real :: f
         complex, dimension(nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT2R2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,i&
     &ndy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: ttp
         real :: f
         complex, dimension(2,nyv,kxp,jblok) :: g
         complex, dimension(2,kxp,kyp,kblok) :: bs
         complex, dimension(2,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT2R3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,i&
     &ndy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: ttp
         real :: f
         complex, dimension(3,nyv,kxp,jblok) :: g
         complex, dimension(3,kxp,kyp,kblok) :: bs
         complex, dimension(3,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT2RX2(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,k&
     &strt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: ttp
         real :: f
         complex, dimension(2,nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT2RX3(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,k&
     &strt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
         real :: ttp
         real :: f
         complex, dimension(3,nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WPFFT2RN(f,g,bs,br,ss,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, ndim, nxhyd, nxyhd
         real :: ttp
         real :: f
         complex, dimension(ndim,nyv,kxp,jblok) :: g
         complex, dimension(ndim,kxp,kyp,kblok) :: bs
         complex, dimension(ndim,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         complex, dimension(ndim,nxvh) :: ss
         end subroutine
      end interface
      interface
         subroutine WPFFT2RXN(f,g,ss,isign,ntpose,mixup,sct,ttp,indx,ind&
     &y,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, ndim, nxhyd, nxyhd
         real :: ttp
         real :: f
         complex, dimension(ndim,nyv,kxp,jblok) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         complex, dimension(ndim,nxvh) :: ss
         end subroutine
      end interface
      interface
         subroutine WP2FFT2RN(f1,f2,g1,g2,bs,br,ss,isign,ntpose,mixup,sc&
     &t,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim1,ndim&
     &2,nxhyd,nxyhd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp, kyp, kypd, jblok, kblok, ndim1, ndim2
         integer :: nxhyd, nxyhd
         real :: ttp
         real :: f1, f2
         complex, dimension(ndim1,nyv,kxp,jblok) :: g1
         complex, dimension(ndim2,nyv,kxp,jblok) :: g2
         complex, dimension(ndim1+ndim2,kxp,kyp,kblok) :: bs
         complex, dimension(ndim1+ndim2,kxp,kyp,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         complex, dimension(ndim1+ndim2,nxvh) :: ss
         end subroutine
      end interface
      interface
         subroutine WPFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
         implicit none
         integer :: indx, indy
         integer :: nxhyd, nxyd
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFSST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         real, dimension(nyv,kxp2d,jblok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFSCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         real, dimension(nyv,kxp2d,jblok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFCST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         real, dimension(nyv,kxp2d,jblok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFCCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         real, dimension(nyv,kxp2d,jblok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFCST2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         real, dimension(2,nyv,kxp2d,jblok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFSCT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         real, dimension(2,nyv,kxp2d,jblok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFCST2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         real, dimension(3,nyv,kxp2d,jblok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFSCT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx&
     &,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyv
         integer :: kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         real, dimension(3,nyv,kxp2d,jblok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,jblok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFSFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(nyvh,kxp2d,j2blok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFCFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,&
     &indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(nyvh,kxp2d,j2blok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFCSFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,ind&
     &x,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd&
     &)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(2,nyvh,kxp2d,j2blok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFSCFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,ind&
     &x,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd&
     &)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(2,nyvh,kxp2d,j2blok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFCSFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,ind&
     &x,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd&
     &)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(3,nyvh,kxp2d,j2blok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFSCFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,ind&
     &x,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd&
     &)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(3,nyvh,kxp2d,j2blok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WPFDT2RINIT(sctdx,indx,nxd)
         implicit none
         integer :: indx, nxd
         complex, dimension(nxd) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WPFDSFT2RX(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,t&
     &tp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhy&
     &d,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(nyvh,kxp2d,j2blok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WPFDCFT2RX(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,t&
     &tp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhy&
     &d,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(nyvh,kxp2d,j2blok) :: g
         real, dimension(kxp2+1,kyp+1,kblok) :: bs
         real, dimension(kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WPFDCSFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,&
     &ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxh&
     &yd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(2,nyvh,kxp2d,j2blok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WPFDSCFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,&
     &ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxh&
     &yd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(2,nyvh,kxp2d,j2blok) :: g
         real, dimension(2,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(2,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WPFDCSFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,&
     &ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxh&
     &yd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(3,nyvh,kxp2d,j2blok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WPFDSCFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,&
     &ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxh&
     &yd,nxyd)
         implicit none
         integer :: isign, ntpose, indx, indy, kstrt, nxvh, nyvh
         integer :: kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
         real :: ttp
         real :: f
         complex, dimension(3,nyvh,kxp2d,j2blok) :: g
         real, dimension(3,kxp2+1,kyp+1,kblok) :: bs
         real, dimension(3,kxp2+1,kyp+1,j2blok) :: br
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxvh) :: sctdx
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface fft_init
!        module procedure ipfft2rinit
!        module procedure ipfft2rxinit
         module procedure iwpfft2rinit
      end interface
!
      interface fft
!        module procedure ipfft2r
!        module procedure ipfft2r2
!        module procedure ipfft2r3
!        module procedure ipfft2rx
!        module procedure ipfft2rx2
!        module procedure ipfft2rx3
         module procedure iwpfft2r
         module procedure iwpfft2r2
         module procedure iwpfft2r3
!        module procedure iwpfft2rx
!        module procedure iwpfft2rx2
!        module procedure iwpfft2rx3
         module procedure ipfft2c
      end interface
!
      interface fftn
         module procedure iwpfft2rn
!        module procedure iwpfft2rxn
         module procedure iwp1fft2rn
         module procedure iwp2fft2rn
      end interface
!
      interface fftc_init
         module procedure ipfft2cinit
      end interface
!
      interface fst_init
         module procedure iwpfst2rinit
      end interface
!
      interface fsst
         module procedure iwpfsst2r
      end interface
!
      interface fcct
         module procedure iwpfcct2r
      end interface
!
      interface fcst
         module procedure iwpfcst2r3
      end interface
!
      interface fsct
         module procedure iwpfsct2r3
      end interface
!
      interface fsft
         module procedure iwpfsft2r
      end interface
!
      interface fcft
         module procedure iwpfcft2r
      end interface
!
      interface fcsft
         module procedure iwpfcsft2r3
      end interface
!
      interface fscft
         module procedure iwpfscft2r3
      end interface
!
      interface fdt_init
         module procedure iwpfdt2rinit
      end interface
!
      interface fdsft
         module procedure iwpfdsft2rx
      end interface
!
      interface fdcft
         module procedure iwpfdcft2rx
      end interface
!
      interface fdcsft
         module procedure iwpfdcsft2r3
      end interface
!
      interface fdscft
         module procedure iwpfdscft2r3
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipfft2rinit(mixup,sct,indx,indy)
! initialize 2d real to complex fft
         implicit none
         integer :: indx, indy
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 0, ntpose = 1, kstrt = 1, nxvh = 1, nyv = 1
         integer :: kxp = 1, kyp = 1, kypd = 1, jblok = 1, kblok = 1
         integer :: nxhyd, nxyhd
         real :: f
         complex, dimension(1,1,1) :: g, bs, br
         nxhyd = size(mixup); nxyhd = size(sct)
         call PFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,kstrt,nx&
     &vh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         end subroutine ipfft2rinit
!
         subroutine ipfft2r(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,kyp&
     &,inorder)
! perform 2d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order
!        complex, dimension(size(g,2),kyp,size(f,3)) :: bs
!        complex, dimension(size(g,2),kyp,size(g,3)) :: br
         real :: tf
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyhd = size(sct)
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
            call PFFT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sct,indx,ind&
     &y,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         else
            call PFFT2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sct,indx,ind&
     &y,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft2r
!
         subroutine ipfft2r2(f,g,mixup,sct,tfft,indx,indy,kstrt,kyp,inor&
     &der)
! perform 2d vector real to complex fft for 2 component vectors
         implicit none
         integer :: indx, indy, kstrt, kyp
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 1, ntpose = 1, nxvh, nyv
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order
!        complex, dimension(2,size(g,3),kyp,size(f,4)) :: bs
!        complex, dimension(2,size(g,3),kyp,size(g,4)) :: br
         real :: tf
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
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
            call PFFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,indx,&
     &indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         else
            call PFFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,indx,&
     &indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft2r2
!
         subroutine ipfft2r3(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,ky&
     &p,inorder)
! perform 2d vector real to complex fft for 3 component vectors
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order
!        complex, dimension(3,size(g,3),kyp,size(f,4)) :: bs
!        complex, dimension(3,size(g,3),kyp,size(g,4)) :: br
         real :: tf
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
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
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call PFFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            else
               call PFFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call PFFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            else
               call PFFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            endif
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft2r3
!
         subroutine ipfft2rxinit(mixup,sct,indx,indy)
! initialize optimized 2d real to complex fft
         implicit none
         integer :: indx, indy
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 0, ntpose = 1, kstrt = 1, nxvh = 1, nyv = 1
         integer :: kxp = 1, kyp = 1, kypd = 1
         integer :: jblok = 1, kblok = 1
         integer :: nxhyd, nxyhd
         real :: f
         complex, dimension(1,1,1) :: g
         nxhyd = size(mixup); nxyhd = size(sct)
         call PFFT2RX(f,g,isign,ntpose,mixup,sct,indx,indy,kstrt,nxvh,ny&
     &v,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         end subroutine ipfft2rxinit
!
         subroutine ipfft2rx(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,ky&
     &p,inorder)
! perform optimized 2d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv, kxp, kypd
         integer :: jblok, kblok, nxhyd, nxyhd, order
         real :: tf
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call PFFT2RX(f(1,1,1),g,isign,ntpose,mixup,sct,indx,indy,kst&
     &rt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         else
            call PFFT2RX(f(2,2,1),g,isign,ntpose,mixup,sct,indx,indy,kst&
     &rt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft2rx
!
         subroutine ipfft2rx2(f,g,mixup,sct,tfft,indx,indy,kstrt,kyp,ino&
     &rder)
! perform optimized 2d vector real to complex fft for 2 component vector
         implicit none
         integer :: indx, indy, kstrt, kyp
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 1, ntpose = 1, nxvh, nyv, kxp
         integer :: kypd, jblok, kblok, nxhyd, nxyhd, order
         real :: tf
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call PFFT2RX2(f(1,1,1,1),g,isign,ntpose,mixup,sct,indx,indy,&
     &kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         else
            call PFFT2RX2(f(1,2,2,1),g,isign,ntpose,mixup,sct,indx,indy,&
     &kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft2rx2
!
         subroutine ipfft2rx3(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,k&
     &yp,inorder)
! perform optimized 2d vector real to complex fft for 3 component vector
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv, kxp
         integer :: kypd, jblok, kblok, nxhyd, nxyhd, order
         real :: tf
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call PFFT2RX2(f(1,1,1,1),g,isign,ntpose,mixup,sct,indx,in&
     &dy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            else
               call PFFT2RX2(f(1,2,2,1),g,isign,ntpose,mixup,sct,indx,in&
     &dy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call PFFT2RX3(f(1,1,1,1),g,isign,ntpose,mixup,sct,indx,in&
     &dy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            else
               call PFFT2RX3(f(1,2,2,1),g,isign,ntpose,mixup,sct,indx,in&
     &dy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            endif
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft2rx3
!
         subroutine iwpfft2rinit(mixup,sct,indx,indy)
! initialize 2d real to complex fft
         implicit none
         integer :: indx, indy
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhyd, nxyhd
         nxhyd = size(mixup); nxyhd = size(sct)
         call WPFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         end subroutine iwpfft2rinit
!
         subroutine iwpfft2r(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,ky&
     &p,inorder)
! perform 2d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order
!        complex, dimension(size(g,2),kyp,size(f,3)) :: bs
!        complex, dimension(size(g,2),kyp,size(g,3)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyhd = size(sct)
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
            call WPFFT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         else
            call WPFFT2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sct,ttp,ind&
     &x,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft2r
!
         subroutine iwpfft2r2(f,g,mixup,sct,tfft,indx,indy,kstrt,kyp,ino&
     &rder)
! perform 2d vector real to complex fft for 2 component vectors
         implicit none
         integer :: indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 1, ntpose = 1, nxvh, nyv
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order
!        complex, dimension(2,size(g,3),kyp,size(f,4)) :: bs
!        complex, dimension(2,size(g,3),kyp,size(g,4)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
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
            call WPFFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         else
            call WPFFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,ttp,&
     &indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft2r2
!
         subroutine iwpfft2r3(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,k&
     &yp,inorder)
! perform 2d vector real to complex fft for 3 component vectors
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order
!        complex, dimension(3,size(g,3),kyp,size(f,4)) :: bs
!        complex, dimension(3,size(g,3),kyp,size(g,4)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
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
               call WPFFT2R(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,tt&
     &p,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            else
               call WPFFT2R(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,tt&
     &p,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            endif
         case (2)
            if (order==LINEAR) then
               call WPFFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            else
               call WPFFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            endif
         case (3)
            if (order==LINEAR) then
               call WPFFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            else
               call WPFFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            endif
         end select
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft2r3
!
         subroutine iwpfft2rx(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,k&
     &yp,inorder)
! perform optimized 2d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv, kxp, kypd
         integer :: jblok, kblok, nxhyd, nxyhd, order
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFFT2RX(f(1,1,1),g,isign,ntpose,mixup,sct,ttp,indx,ind&
     &y,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         else
            call WPFFT2RX(f(2,2,1),g,isign,ntpose,mixup,sct,ttp,indx,ind&
     &y,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft2rx
!
         subroutine iwpfft2rx2(f,g,mixup,sct,tfft,indx,indy,kstrt,kyp,in&
     &order)
! perform optimized 2d vector real to complex fft for 2 component vector
         implicit none
         integer :: indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 1, ntpose = 1, nxvh, nyv, kxp
         integer :: kypd, jblok, kblok, nxhyd, nxyhd, order
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFFT2RX2(f(1,1,1,1),g,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         else
            call WPFFT2RX2(f(1,2,2,1),g,isign,ntpose,mixup,sct,ttp,indx,&
     &indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft2rx2
!
         subroutine iwpfft2rx3(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,&
     &kyp,inorder)
! perform optimized 2d vector real to complex fft for 3 component vector
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxvh, nyv, kxp
         integer :: kypd, jblok, kblok, nxhyd, nxyhd, order
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call WPFFT2RX2(f(1,1,1,1),g,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            else
               call WPFFT2RX2(f(1,2,2,1),g,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call WPFFT2RX3(f(1,1,1,1),g,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            else
               call WPFFT2RX3(f(1,2,2,1),g,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd)
            endif
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft2rx3
!
         subroutine iwpfft2rn(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,k&
     &yp,inorder)
! perform 2d vector real to complex fft for n component vectors
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, ndim, nxvh, nyv
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order
!        complex, dimension(size(f,1),size(g,3),kyp,size(f,4)) :: bs
!        complex, dimension(size(f,1),size(g,3),kyp,size(g,4)) :: br
         complex, dimension(size(f,1),size(f,2)/2) :: ss
         real :: tf, ttp
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
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
            call WPFFT2RN(f(1,1,1,1),g,bs,br,ss,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nx&
     &yhd)
         else
            call WPFFT2RN(f(1,2,2,1),g,bs,br,ss,isign,ntpose,mixup,sct,t&
     &tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nx&
     &yhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft2rn
!
         subroutine iwpfft2rxn(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,&
     &kyp,inorder)
! perform optimized 2d vector real to complex fft for n component vector
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         complex, dimension(size(f,1),size(f,2)/2) :: ss
         integer :: ntpose = 1, ndim, nxvh, nyv, kxp
         integer :: kypd, jblok, kblok, nxhyd, nxyhd, order
         real :: tf, ttp
         double precision :: dtime
         ndim = size(f,1); nxvh = size(f,2)/2; kypd = size(f,3)
         kblok = size(f,4)
         nyv = size(g,2); kxp = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyhd = size(sct)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFFT2RXN(f(1,1,1,1),g,ss,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd)
         else
            call WPFFT2RXN(f(1,2,2,1),g,ss,isign,ntpose,mixup,sct,ttp,in&
     &dx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfft2rxn
!
         subroutine iwp1fft2rn(f1,f2,g1,g2,isign,mixup,sct,tfft,indx,ind&
     &y,kstrt,kyp,inorder)
! perform 2d real to complex fft for a scalar and n component vector
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
         integer :: ntpose = 1, ndim1, ndim2, ndim, nxvh, nyv
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order
!        complex, dimension(1+size(f2,1),size(g1,2),kyp,size(f1,3)) ::  &
!    &bs
!        complex, dimension(1+size(f2,1),size(g1,2),kyp,size(g1,3)) ::  &
!    &br
         complex, dimension(1+size(f2,1),size(f1,1)/2) :: ss
         real :: tf, ttp
         double precision :: dtime
         ndim1 = 1; ndim2 = size(f2,1); ndim = ndim1 + ndim2
         nxvh = size(f1,1)/2; kypd = size(f1,2); kblok = size(f1,3)
         nyv = size(g1,1); kxp = size(g1,2); jblok = size(g1,3)
         nxhyd = size(mixup); nxyhd = size(sct)
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
            call WP2FFT2RN(f1(1,1,1),f2(1,1,1,1),g1,g2,bs,br,ss,isign,nt&
     &pose,mixup,sct,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kbl&
     &ok,ndim1,ndim2,nxhyd,nxyhd)
         else
            call WP2FFT2RN(f1(2,2,1),f2(1,2,2,1),g1,g2,bs,br,ss,isign,nt&
     &pose,mixup,sct,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kbl&
     &ok,ndim1,ndim2,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwp1fft2rn
!
         subroutine iwp2fft2rn(f1,f2,g1,g2,isign,mixup,sct,tfft,indx,ind&
     &y,kstrt,kyp,inorder)
! perform 2d real to complex fft for two n component vectors
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f1, f2
         complex, dimension(:,:,:,:), pointer :: g1, g2
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, ndim1, ndim2, ndim, nxvh, nyv
         integer :: kxp, kypd, jblok, kblok, nxhyd, nxyhd, order
!        complex, dimension(size(f1,1)+size(f2,1),size(g1,3),kyp,size(f1&
!    &,4)) :: bs
!        complex, dimension(size(f1,1)+size(f2,1),size(g1,3),kyp,size(g1&
!    &,4)) :: br
         complex, dimension(size(f1,1)+size(f2,1),size(f1,2)/2) :: ss
         real :: tf, ttp
         double precision :: dtime
         ndim1 = size(f1,1); ndim2 = size(f2,1); ndim = ndim1 + ndim2
         nxvh = size(f1,2)/2; kypd = size(f1,3); kblok = size(f1,4)
         nyv = size(g1,2); kxp = size(g1,3); jblok = size(g1,4)
         nxhyd = size(mixup); nxyhd = size(sct)
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
            call WP2FFT2RN(f1(1,1,1,1),f2(1,1,1,1),g1,g2,bs,br,ss,isign,&
     &ntpose,mixup,sct,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,k&
     &blok,ndim1,ndim2,nxhyd,nxyhd)
         else
            call WP2FFT2RN(f1(1,2,2,1),f2(1,2,2,1),g1,g2,bs,br,ss,isign,&
     &ntpose,mixup,sct,ttp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,k&
     &blok,ndim1,ndim2,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwp2fft2rn
!
         subroutine ipfft2c(f,g,isign,mixup,sct,tfft,indx,indy,kstrt,kyp&
     &,inorder)
! perform 2d scalar complex to complex fft
         implicit none
         integer :: isign, indx, indy, kstrt, kyp
         integer, optional :: inorder
         real :: tfft
         complex, dimension(:,:,:), pointer :: f, g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ntpose = 1, nxv, nyv
         integer :: kxp, kypd, jblok, kblok, nxyd, nxyhd, order
!        complex, dimension(size(g,2),kyp,size(f,3)) :: bs
!        complex, dimension(size(g,2),kyp,size(g,3)) :: br
         real :: tf
         double precision :: dtime
         nxv = size(f,1); kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp = size(g,2); jblok = size(g,3)
         nxyd = size(mixup); nxyhd = size(sct)
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
            call PFFT2C(f(1,1,1),g,bs,br,isign,ntpose,mixup,sct,indx,ind&
     &y,kstrt,nxv,nyv,kxp,kyp,kypd,jblok,kblok,nxyd,nxyhd)
         else
            call PFFT2C(f(2,2,1),g,bs,br,isign,ntpose,mixup,sct,indx,ind&
     &y,kstrt,nxv,nyv,kxp,kyp,kypd,jblok,kblok,nxyd,nxyhd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft = tfft + tf
         end subroutine ipfft2c
!
         subroutine ipfft2cinit(mixup,sct,indx,indy)
! initialize 2d complex to complex fft
         implicit none
         integer :: indx, indy
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 0, ntpose = 1, kstrt = 1, nxv = 1, nyv = 1
         integer :: kxp = 1, kyp = 1, kypd = 1
         integer :: jblok = 1, kblok = 1
         integer :: nxyd, nxyhd
         complex :: f
         complex, dimension(1,1,1) :: g, bs, br
         nxyd = size(mixup); nxyhd = size(sct)
         call PFFT2C(f,g,bs,br,isign,ntpose,mixup,sct,indx,indy,kstrt,nx&
     &v,nyv,kxp,kyp,kypd,jblok,kblok,nxyd,nxyhd)
         end subroutine ipfft2cinit
!
         subroutine iwpfst2rinit(mixup,sctd,indx,indy)
! initialize 2d real sine-cosine transforms
         implicit none
         integer :: indx, indy
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhyd, nxyd
         nxhyd = size(mixup); nxyd = size(sctd)
         call WPFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
         end subroutine iwpfst2rinit
!
         subroutine iwpfsst2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,&
     &kxp2,kyp,inorder)
! perform 2d scalar real sine-sine transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f, g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyv, kypd, kxp2d
         integer :: jblok, kblok, nxhyd, nxyd, order
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp2d = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFSST2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd&
     &)
         else
            call WPFSST2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd&
     &)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfsst2r
!
         subroutine iwpfcct2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,&
     &kxp2,kyp,inorder)
! perform 2d scalar real cosine-cosine transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f, g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyv, kypd, kxp2d
         integer :: jblok, kblok, nxhyd, nxyd, order
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyv = size(g,1); kxp2d = size(g,2); jblok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFCCT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd&
     &)
         else
            call WPFCCT2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd&
     &)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfcct2r
!
         subroutine iwpfcst2r3(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt&
     &,kxp2,kyp,inorder)
! perform 2d vector real cosine-sine transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f, g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyv, kypd, kxp2d
         integer :: jblok, kblok, nxhyd, nxyd, order
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp2d = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call WPFCST2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd)
            else
               call WPFCST2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call WPFCST2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd)
            else
               call WPFCST2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd)
            endif
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfcst2r3
!
         subroutine iwpfsct2r3(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt&
     &,kxp2,kyp,inorder)
! perform 2d vector real sine-cosine transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f, g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyv, kypd, kxp2d
         integer :: jblok, kblok, nxhyd, nxyd, order
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyv = size(g,2); kxp2d = size(g,3); jblok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call WPFSCT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd)
            else
               call WPFSCT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call WPFSCT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd)
            else
               call WPFSCT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sctd&
     &,ttp,indx,indy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhy&
     &d,nxyd)
            endif
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfsct2r3
!
         subroutine iwpfsft2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,&
     &kxp2,kyp,inorder)
! perform 2d scalar real sine/periodic transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d
         integer :: j2blok, kblok, nxhyd, nxyd, order
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyvh = size(g,1); kxp2d = size(g,2); j2blok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFSFT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nx&
     &yd)
         else
            call WPFSFT2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nx&
     &yd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfsft2r
!
         subroutine iwpfcft2r(f,g,isign,mixup,sctd,tfft,indx,indy,kstrt,&
     &kxp2,kyp,inorder)
! perform 2d scalar real cosine/periodic transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d
         integer :: j2blok, kblok, nxhyd, nxyd, order
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyvh = size(g,1); kxp2d = size(g,2); j2blok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFCFT2R(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nx&
     &yd)
         else
            call WPFCFT2R(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,ttp,i&
     &ndx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nx&
     &yd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfcft2r
!
         subroutine iwpfcsft2r3(f,g,isign,mixup,sctd,tfft,indx,indy,kstr&
     &t,kxp2,kyp,inorder)
! perform 2d vector real cosine-sine/periodic transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d
         integer :: j2blok, kblok, nxhyd, nxyd, order
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyvh = size(g,2); kxp2d = size(g,3); j2blok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call WPFCSFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd)
            else
               call WPFCSFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call WPFCSFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd)
            else
               call WPFCSFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfcsft2r3
!
         subroutine iwpfscft2r3(f,g,isign,mixup,sctd,tfft,indx,indy,kstr&
     &t,kxp2,kyp,inorder)
! perform 2d vector real sine-cosine/periodic transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d
         integer :: j2blok, kblok, nxhyd, nxyd, order
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyvh = size(g,2); kxp2d = size(g,3); j2blok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call WPFSCFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd)
            else
               call WPFSCFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call WPFSCFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd)
            else
               call WPFSCFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sct&
     &d,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,n&
     &xhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfscft2r3
!
         subroutine iwpfdt2rinit(sctdx,indx)
! initialize 2d real sine-cosine/periodic transforms
         implicit none
         integer :: indx
         complex, dimension(:), pointer :: sctdx
! local data
         integer :: nxd
         nxd = size(sctdx)
         call WPFDT2RINIT(sctdx,indx,nxd)
         end subroutine iwpfdt2rinit
!
         subroutine iwpfdsft2rx(f,g,isign,mixup,sctd,sctdx,tfft,indx,ind&
     &y,kstrt,kxp2,kyp,inorder)
! perform 2d scalar real half sine/periodic transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d
         integer :: j2blok, kblok, nxhyd, nxyd, order
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyvh = size(g,1); kxp2d = size(g,2); j2blok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFDSFT2RX(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,sct&
     &dx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,&
     &nxhyd,nxyd)
         else
            call WPFDSFT2RX(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,sct&
     &dx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,&
     &nxhyd,nxyd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfdsft2rx
!
         subroutine iwpfdcft2rx(f,g,isign,mixup,sctd,sctdx,tfft,indx,ind&
     &y,kstrt,kxp2,kyp,inorder)
! perform 2d scalar real half cosine/periodic transform
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d
         integer :: j2blok, kblok, nxhyd, nxyd, order
         real, dimension(kxp2+1,kyp+1,size(f,3)) :: bs
         real, dimension(kxp2+1,kyp+1,size(g,3)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,1)/2; kypd = size(f,2); kblok = size(f,3)
         nyvh = size(g,1); kxp2d = size(g,2); j2blok = size(g,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (order==LINEAR) then
            call WPFDCFT2RX(f(1,1,1),g,bs,br,isign,ntpose,mixup,sctd,sct&
     &dx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,&
     &nxhyd,nxyd)
         else
            call WPFDCFT2RX(f(2,2,1),g,bs,br,isign,ntpose,mixup,sctd,sct&
     &dx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,&
     &nxhyd,nxyd)
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfdcft2rx
!
         subroutine iwpfdcsft2r3(f,g,isign,mixup,sctd,sctdx,tfft,indx,in&
     &dy,kstrt,kxp2,kyp,inorder)
! perform 2d vector real half cosine-sine/periodic transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d
         integer :: j2blok, kblok, nxhyd, nxyd, order
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyvh = size(g,2); kxp2d = size(g,3); j2blok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call WPFDCSFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd)
            else
               call WPFDCSFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call WPFDCSFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd)
            else
               call WPFDCSFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfdcsft2r3
!
         subroutine iwpfdscft2r3(f,g,isign,mixup,sctd,sctdx,tfft,indx,in&
     &dy,kstrt,kxp2,kyp,inorder)
! perform 2d vector real half sine-cosine/periodic transforms
         implicit none
         integer :: isign, indx, indy, kstrt, kxp2, kyp
         integer, optional :: inorder
         real, dimension(2) :: tfft
         real, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:,:), pointer :: g
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         integer :: ntpose = 1, nxvh, nyvh, kypd, kxp2d
         integer :: j2blok, kblok, nxhyd, nxyd, order
         real, dimension(size(f,1),kxp2+1,kyp+1,size(f,4)) :: bs
         real, dimension(size(f,1),kxp2+1,kyp+1,size(g,4)) :: br
         real :: tf, ttp
         double precision :: dtime
         nxvh = size(f,2)/2; kypd = size(f,3); kblok = size(f,4)
         nyvh = size(g,2); kxp2d = size(g,3); j2blok = size(g,4)
         nxhyd = size(mixup); nxyd = size(sctd)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,dtime,-1)
         if (size(f,1)==2) then
            if (order==LINEAR) then
               call WPFDSCFT2R2(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd)
            else
               call WPFDSCFT2R2(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (order==LINEAR) then
               call WPFDSCFT2R3(f(1,1,1,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd)
            else
               call WPFDSCFT2R3(f(1,2,2,1),g,bs,br,isign,ntpose,mixup,sc&
     &td,sctdx,ttp,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,&
     &kblok,nxhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,dtime)
         tfft(1) = tfft(1) + tf
         tfft(2) = tfft(2) + ttp
         end subroutine iwpfdscft2r3
!
      end module pfft2d

