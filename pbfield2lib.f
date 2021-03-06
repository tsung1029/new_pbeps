c 2d parallel PIC library for solving field equations with mixed
c dirichlet/periodic boundary conditions
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: may 18, 2004
c-----------------------------------------------------------------------
      subroutine PMSCGUARD2(cu,nyp,xj0,yj0,zj0,nx,ngx,nxe,nypmx,nblok)
c initialize extended non-periodic field in x and periodic in y
c nyp(l) = number of primary gridpoints in particle partition l
c ngx = (0,1) = number of grid cells away from edge
c quadratic interpolation, for distributed data
      implicit none
      real cu, xj0, yj0, zj0
      integer nyp, nx, ngx, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
      integer i, j, k, l, nyp3, nxg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      do 60 l = 1, nblok
      nyp3 = nyp(l) + 3
c handle interior grid points in y
      do 30 k = 1, nyp(l)
      do 10 j = 2, nxg
      cu(1,j+ngx+1,k+1,l) = xj0
      cu(2,j+ngx+1,k+1,l) = yj0
      cu(3,j+ngx+1,k+1,l) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k+1,l) = 0.
      cu(i,2,k+1,l) = 0.
      cu(i,nx+2,k+1,l) = 0.
      cu(i,nx+3,k+1,l) = 0.
   20 continue
      cu(1,ngx+2,k+1,l) = .5*xj0
      cu(2,ngx+2,k+1,l) = .5*yj0
      cu(3,ngx+2,k+1,l) = .5*zj0
      cu(1,nx-ngx+2,k+1,l) = .5*xj0
      cu(2,nx-ngx+2,k+1,l) = .5*yj0
      cu(3,nx-ngx+2,k+1,l) = .5*zj0
   30 continue
c zero out guard cells
      do 50 j = 1, nx3
      do 40 i = 1, 3
      cu(i,j,1,l) = 0.
      cu(i,j,nyp3-1,l) = 0.
      cu(i,j,nyp3,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMSCGUARD22(cu,nyp,xj0,yj0,nx,ngx,nxe,nypmx,nblok)
c initialize extended non-periodic field in x and periodic in y
c nyp(l) = number of primary gridpoints in particle partition l
c ngx = (0,1) = number of grid cells away from edge
c quadratic interpolation, for distributed data
      implicit none
      real cu, xj0, yj0
      integer nyp, nx, ngx, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
      integer i, j, k, l, nyp3, nxg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      do 60 l = 1, nblok
      nyp3 = nyp(l) + 3
c handle interior grid points in y
      do 30 k = 1, nyp(l)
      do 10 j = 2, nxg
      cu(1,j+ngx+1,k+1,l) = xj0
      cu(2,j+ngx+1,k+1,l) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k+1,l) = 0.
      cu(i,2,k+1,l) = 0.
      cu(i,nx+2,k+1,l) = 0.
      cu(i,nx+3,k+1,l) = 0.
   20 continue
      cu(1,ngx+2,k+1,l) = .5*xj0
      cu(2,ngx+2,k+1,l) = .5*yj0
      cu(1,nx-ngx+2,k+1,l) = .5*xj0
      cu(2,nx-ngx+2,k+1,l) = .5*yj0
   30 continue
c zero out guard cells
      do 50 j = 1, nx3
      do 40 i = 1, 2
      cu(i,j,1,l) = 0.
      cu(i,j,nyp3-1,l) = 0.
      cu(i,j,nyp3,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMSGUARD2(q,nyp,qi0,nx,ngx,nxe,nypmx,nblok)
c initialize extended non-periodic scalar field in x and periodic in y
c nyp(l) = number of primary gridpoints in particle partition l
c ngx = (0,1) = number of grid cells away from edge
c quadratic interpolation, for distributed data
      implicit none
      real q, qi0
      integer nyp, nx, ngx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
      integer j, k, l, nyp3, nxg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      do 40 l = 1, nblok
      nyp3 = nyp(l) + 3
c handle interior grid points in y
      do 20 k = 1, nyp(l)
      do 10 j = 2, nxg
      q(j+ngx+1,k+1,l) = qi0
   10 continue
      q(1,k+1,l) = 0.
      q(2,k+1,l) = 0.
      q(nx+2,k+1,l) = 0.
      q(nx+3,k+1,l) = 0.
      q(ngx+2,k+1,l) = .5*qi0
      q(nx-ngx+2,k+1,l) = .5*qi0
   20 continue
c zero out guard cells
      do 30 j = 1, nx3
      q(j,1,l) = 0.
      q(j,nyp3-1,l) = 0.
      q(j,nyp3,l) = 0.
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMSCGUARD2L(cu,nyp,xj0,yj0,zj0,nx,ngx,nxe,nypmx,nblok)
c initialize extended non-periodic field in x and periodic in y
c nyp(l) = number of primary gridpoints in particle partition l
c ngx = (0,1) = number of grid cells away from edge
c linear interpolation, for distributed data
      implicit none
      real cu, xj0, yj0, zj0
      integer nyp, nx, ngx, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
      integer i, j, k, l, nyp1, nxg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      do 60 l = 1, nblok
      nyp1 = nyp(l) + 1
c handle interior grid points in y
      do 30 k = 1, nyp(l)
      do 10 j = 2, nxg
      cu(1,j+ngx,k,l) = xj0
      cu(2,j+ngx,k,l) = yj0
      cu(3,j+ngx,k,l) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k,l) = 0.
      cu(i,nx+1,k,l) = 0.
   20 continue
      cu(1,ngx+1,k,l) = .5*xj0
      cu(2,ngx+1,k,l) = .5*yj0
      cu(3,ngx+1,k,l) = .5*zj0
      cu(1,nx-ngx+1,k,l) = .5*xj0
      cu(2,nx-ngx+1,k,l) = .5*yj0
      cu(3,nx-ngx+1,k,l) = .5*zj0
   30 continue
c zero out guard cells
      do 50 j = 1, nx1
      do 40 i = 1, 3
      cu(i,j,nyp1,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMSCGUARD22L(cu,nyp,xj0,yj0,nx,ngx,nxe,nypmx,nblok)
c initialize extended non-periodic field in x and periodic in y
c nyp(l) = number of primary gridpoints in particle partition l
c ngx = (0,1) = number of grid cells away from edge
c linear interpolation, for distributed data
      implicit none
      real cu, xj0, yj0
      integer nyp, nx, ngx, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
      integer i, j, k, l, nyp1, nxg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      do 60 l = 1, nblok
      nyp1 = nyp(l) + 1
c handle interior grid points in y
      do 30 k = 1, nyp(l)
      do 10 j = 2, nxg
      cu(1,j+ngx,k,l) = xj0
      cu(2,j+ngx,k,l) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k,l) = 0.
      cu(i,nx+1,k,l) = 0.
   20 continue
      cu(1,ngx+1,k,l) = .5*xj0
      cu(2,ngx+1,k,l) = .5*yj0
      cu(1,nx-ngx+1,k,l) = .5*xj0
      cu(2,nx-ngx+1,k,l) = .5*yj0
   30 continue
c zero out guard cells
      do 50 j = 1, nx1
      do 40 i = 1, 2
      cu(i,j,nyp1,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMSGUARD2L(q,nyp,qi0,nx,ngx,nxe,nypmx,nblok)
c initialize extended non-periodic scalar field in x and periodic in y
c nyp(l) = number of primary gridpoints in particle partition l
c ngx = (0,1) = number of grid cells away from edge
c linear interpolation, for distributed data
      implicit none
      real q, qi0
      integer nyp, nx, ngx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
      integer j, k, l, nyp1, nxg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      do 40 l = 1, nblok
      nyp1 = nyp(l) + 1
c handle interior grid points in y
      do 20 k = 1, nyp(l)
      do 10 j = 2, nxg
      q(j+ngx,k,l) = qi0
   10 continue
      q(1,k,l) = 0.
      q(nx+1,k,l) = 0.
      q(ngx+1,k,l) = .5*qi0
      q(nx-ngx+1,k,l) = .5*qi0
   20 continue
c zero out guard cells
      do 30 j = 1, nx1
      q(j,nyp1,l) = 0.
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSGLSIN2C(cu,cu1,nx,kyp,nxv,nypmx,kypd,kblok)
c this subroutine creates a doubled vector array cu1 from a vector array
c cu, so that various 1d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in y,
c and y component is an odd function in x.
c Asummes vector cu vanishes at end points
c linear interpolation for distributed data
c nx = system length in x direction
c kyp = number of data values per block in y
c nxv = second dimension of input array cu, must be >= nx
c nypmx = third dimension of input array cu, must be >= kyp
c kzpd = third dimension of output array cu1, must be >= kyp
c kblok = number of data blocks in y
      implicit none
      real cu, cu1
      integer nx, kyp, nxv, nypmx, kypd, kblok
      dimension cu(2,nxv,nypmx,kblok), cu1(2,2*nxv,kypd,kblok)
c local data
      integer j, k, l, nxs
c copy to double array
      nxs = nx - 1
      do 30 l = 1, kblok
      do 20 k = 1, kyp
      do 10 j = 1, nxs
      cu1(1,j+1,k,l) = cu(1,j+1,k,l)
      cu1(2,j+1,k,l) = cu(2,j+1,k,l)
      cu1(1,nx+j+1,k,l) = cu(1,nx-j+1,k,l)
      cu1(2,nx+j+1,k,l) = -cu(2,nx-j+1,k,l)
   10 continue
      cu1(1,1,k,l) = 0.
      cu1(2,1,k,l) = 0.
      cu1(1,nx+1,k,l) = 0.
      cu1(2,nx+1,k,l) = 0.
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSGLSIN2D(q,q1,nx,kyp,nxv,nypmx,kypd,kblok)
c this subroutine creates a doubled array q1 from an array q, so that
c a 1d sine transform can be performed with a 2d real to complex fft.
c linear interpolation for distributed data
c nx = system length in x direction
c kyp = number of data values per block in y
c nxv = second dimension of input array cu, must be >= nx
c nypmx = third dimension of input array cu, must be >= kyp
c kzpd = third dimension of output array cu1, must be >= kyp
c kblok = number of data blocks in y
      implicit none
      real q, q1
      integer nx, kyp, nxv, nypmx, kypd, kblok
      dimension q(nxv,nypmx,kblok), q1(2*nxv,kypd,kblok)
c local data
      integer j, k, l, nxs
c copy to double array
      nxs = nx - 1
      do 30 l = 1, kblok
      do 20 k = 1, kyp
      do 10 j = 1, nxs
      q1(j+1,k,l) = q(j+1,k,l)
      q1(nx+j+1,k,l) = -q(nx-j+1,k,l)
   10 continue
      q1(1,k,l) = 0.
      q1(nx+1,k,l) = 0.
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSGLSIN2B(cu,cu1,nx,kyp,nxv,nypmx,kypd,kblok)
c this subroutine creates a doubled vector array cu1 from a vector array
c cu, so that various 1d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in y,
c and y component is an odd function in x.
c asummes vector cu vanishes at end points
c linear interpolation for distributed data
c nx = system length in x direction
c kyp = number of data values per block in y
c nxv = second dimension of input array cu, must be >= nx
c nypmx = third dimension of input array cu, must be >= kyp
c kzpd = third dimension of output array cu1, must be >= kyp
c kblok = number of data blocks in y
      implicit none
      real cu, cu1
      integer nx, kyp, nxv, nypmx, kypd, kblok
      dimension cu(3,nxv,nypmx,kblok), cu1(3,2*nxv,kypd,kblok)
c local data
      integer j, k, l, nxs
c copy to double array
      nxs = nx - 1
      do 30 l = 1, kblok
      do 20 k = 1, kyp
      do 10 j = 1, nxs
      cu1(1,j+1,k,l) = cu(1,j+1,k,l)
      cu1(2,j+1,k,l) = cu(2,j+1,k,l)
      cu1(3,j+1,k,l) = cu(3,j+1,k,l)
      cu1(1,nx+j+1,k,l) = cu(1,nx-j+1,k,l)
      cu1(2,nx+j+1,k,l) = -cu(2,nx-j+1,k,l)
      cu1(3,nx+j+1,k,l) = -cu(3,nx-j+1,k,l)
   10 continue
      cu1(1,1,k,l) = 0.
      cu1(2,1,k,l) = 0.
      cu1(3,1,k,l) = 0.
      cu1(1,nx+1,k,l) = 0.
      cu1(2,nx+1,k,l) = 0.
      cu1(3,nx+1,k,l) = 0.
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFSGL2C(fxy,fxy1,nx,kyp,nxv,nypmx,kypd,kblok)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c nx = system length in x direction
c nxv = second dimension of input array cu, must be >= nx
c nypmx = third dimension of input array cu, must be >= kyp
c kzpd = third dimension of output array cu1, must be >= kyp
c kblok = number of data blocks in y
      implicit none
      real fxy, fxy1
      integer nx, kyp, nxv, nypmx, kypd, kblok
      dimension fxy(2,nxv,nypmx,kblok), fxy1(2,2*nxv,kypd,kblok)
c local data
      integer j, k, l, nx1
      nx1 = nx + 1
      do 30 l = 1, kblok
      do 20 k = 1, kyp
      do 10 j = 1, nx1
      fxy(1,j,k,l) = fxy1(1,j,k,l)
      fxy(2,j,k,l) = fxy1(2,j,k,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFSGL2D(q,q1,nx,kyp,nxv,nypmx,kypd,kblok)
c this subroutine copies data from a double array to regular array
c with guard cells for scalar field and linear interpolation
c nx = system length in x direction
c nxv = second dimension of input array cu, must be >= nx
c nypmx = third dimension of input array cu, must be >= kyp
c kzpd = third dimension of output array cu1, must be >= kyp
c kblok = number of data blocks in y
      implicit none
      real q, q1
      integer nx, kyp, nxv, nypmx, kypd, kblok
      dimension q(nxv,nypmx,kblok), q1(2*nxv,kypd,kblok)
c local data
      integer j, k, l, nx1
      nx1 = nx + 1
      do 30 l = 1, kblok
      do 20 k = 1, kyp
      do 10 j = 1, nx1
      q(j,k,l) = q1(j,k,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFSGL2B(fxy,fxy1,nx,kyp,nxv,nypmx,kypd,kblok)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c nx = system length in x direction
c nxv = second dimension of input array cu, must be >= nx
c nypmx = third dimension of input array cu, must be >= kyp
c kzpd = third dimension of output array cu1, must be >= kyp
c kblok = number of data blocks in y
      implicit none
      real fxy, fxy1
      integer nx, kyp, nxv, nypmx, kypd, kblok
      dimension fxy(3,nxv,nypmx,kblok), fxy1(3,2*nxv,kypd,kblok)
c local data
      integer j, k, l, nx1
      nx1 = nx + 1
      do 30 l = 1, kblok
      do 20 k = 1, kyp
      do 10 j = 1, nx1
      fxy(1,j,k,l) = fxy1(1,j,k,l)
      fxy(2,j,k,l) = fxy1(2,j,k,l)
      fxy(3,j,k,l) = fxy1(3,j,k,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISMX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,ny
     1v,kxp2,j2blok,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function,
c with mixed dirichlet/periodic boundary conditions for distributed data
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd
c output: ffb
c for isign = -1, input: q,ffb,isign,nx,ny,kstrt,nyv,kxp2,j2blok,nyhd,
c output: fx,fy,we
c approximate flop count is: 24*nxc*nyc + 6*(nxc + nyc)
c for isign = 1, input: q,ffb,isign,nx,ny,kstrt,nyv,kxp2,j2blok,nyhd,
c output: fx,we
c approximate flop count is: 14*nxc*nyc + 5*(nxc + nyc)
c for isign = 2, input: q,ffb,isign,nx,ny,kstrt,nyv,kxp2,j2blok,nyhd,
c output: fy
c approximate flop count is: 3*nxc*nyc + 1*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky)*s(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky)*s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fx(k,j,l) = x component of complex force/charge,
c fy(k,j,l) = y component of complex force/charge,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c real(ffb(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fx, fy, ffb, zero, zt1, zt2
      dimension q(nyv,kxp2,j2blok)
      dimension fx(nyv,kxp2,j2blok), fy(nyv,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffb(k,j,l) = cmplx(affp,1.)
      else
         ffb(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffb(k,j,l))*aimag(ffb(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         zt1 = cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         zt2 = at2*zt1
         fx(k,j,l) = zt2
         fx(k1,j,l) = conjg(zt2)
         zt2 = at3*zt1
         fy(k,j,l) = zt2
         fy(k1,j,l) = -conjg(zt2)
         wp = wp + 2.0*at1*q(k,j,l)*conjg(q(k,j,l))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffb(1,j,l))*aimag(ffb(1,j,l))
         fx(1,j,l) = cmplx(dkx*at1*aimag(q(1,j,l)),0.)
         fx(k1,j,l) = zero
         fy(1,j,l) = zero
         fy(k1,j,l) = zero
         wp = wp + at1*aimag(q(1,j,l))**2
      endif
   60 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 70 k = 1, ny
         fx(k,1,l) = zero
         fy(k,1,l) = zero
   70    continue
      endif
   80 continue
   90 continue
      we = float(nx*ny)*wp
      return
c calculate potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, nyh
         k1 = ny2 - k
         at2 = real(ffb(k,j,l))
         at1 = at2*aimag(ffb(k,j,l))
         zt1 = at2*q(k,j,l)
         fx(k,j,l) = zt1
         fx(k1,j,l) = -conjg(zt1)
         wp = wp + 2.0*at1*q(k,j,l)*conjg(q(k,j,l))
  110    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = real(ffb(1,j,l))
         at1 = at2*aimag(ffb(1,j,l))
         fx(1,j,l) = cmplx(0.,at2*aimag(q(1,j,l)))
         fx(k1,j,l) = zero
         wp = wp + at1*aimag(q(1,j,l))**2
      endif
  120 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 130 k = 1, ny
         fx(k,1,l) = zero
  130    continue
      endif
  140 continue
  150 continue
      we = float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffb(k,j,l))
         zt1 = at1*q(k,j,l)
         fy(k,j,l) = zt1
         fy(k1,j,l) = -conjg(zt1)
  170    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffb(1,j,l))
         fy(1,j,l) = cmplx(0.,at1*aimag(q(1,j,l)))
         fy(k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 190 k = 1, ny
         fy(k,1,l) = zero
  190    continue
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISMX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,nyv
     1,kxp2,j2blok,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with mixed dirichlet/periodic boundary conditions for distributed data
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd,
c output: ffb
c for isign /= 0, input: q,ffb,isign,nx,ny,kstrt,nyv,kxp2,j2blok,nyhd,
c output: fxy,we
c approximate flop count is: 24*nxc*nyc + 6*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c fxy(3,k,j,l) = zero,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated.
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c real(ffb(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffb, zero, zt1, zt2
      dimension q(nyv,kxp2,j2blok), fxy(2,nyv,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffb(k,j,l) = cmplx(affp,1.)
      else
         ffb(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffb(k,j,l))*aimag(ffb(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         zt2 = cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         zt1 = at2*zt2
         zt2 = at3*zt2
         fxy(1,k,j,l) = zt1
         fxy(2,k,j,l) = zt2
         fxy(1,k1,j,l) = conjg(zt1)
         fxy(2,k1,j,l) = -conjg(zt2)
         wp = wp + 2.0*at1*q(k,j,l)*conjg(q(k,j,l))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffb(1,j,l))*aimag(ffb(1,j,l))
         fxy(1,1,j,l) = cmplx(dkx*at1*aimag(q(1,j,l)),0.)
         fxy(2,1,j,l) = zero
         fxy(1,k1,j,l) = zero
         fxy(2,k1,j,l) = zero
         wp = wp + at1*aimag(q(1,j,l))**2
      endif
   60 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 70 k = 1, ny
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = zero
   70    continue
      endif
   80 continue
   90 continue
      we = float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISMX23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,nyv
     1,kxp2,j2blok,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with mixed dirichlet/periodic boundary conditions for distributed data
c Zeros out z component.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd,
c output: ffb
c for isign /= 0, input: q,ffb,isign,nx,ny,kstrt,nyv,kxp2,j2blok,nyhd,
c output: fxy,we
c approximate flop count is: 24*nxc*nyc + 6*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c fxy(3,k,j,l) = zero,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated.
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c real(ffb(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffb, zero, zt1, zt2
      dimension q(nyv,kxp2,j2blok), fxy(3,nyv,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffb(k,j,l) = cmplx(affp,1.)
      else
         ffb(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         at1 = real(ffb(k,j,l))*aimag(ffb(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         zt2 = cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         zt1 = at2*zt2
         zt2 = at3*zt2
         fxy(1,k,j,l) = zt1
         fxy(2,k,j,l) = zt2
         fxy(3,k,j,l) = zero
         fxy(1,k1,j,l) = conjg(zt1)
         fxy(2,k1,j,l) = -conjg(zt2)
         fxy(3,k1,j,l) = zero
         wp = wp + 2.0*at1*q(k,j,l)*conjg(q(k,j,l))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = real(ffb(1,j,l))*aimag(ffb(1,j,l))
         fxy(1,1,j,l) = cmplx(dkx*at1*aimag(q(1,j,l)),0.)
         fxy(2,1,j,l) = zero
         fxy(3,1,j,l) = zero
         fxy(1,k1,j,l) = zero
         fxy(2,k1,j,l) = zero
         fxy(3,k1,j,l) = zero
         wp = wp + at1*aimag(q(1,j,l))**2
      endif
   60 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 70 k = 1, ny
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = zero
         fxy(3,k,1,l) = zero
   70    continue
      endif
   80 continue
   90 continue
      we = float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISM23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,kstrt,nyvh
     1,kxp2,j2blok,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with mixed dirichlet/periodic boundary conditions for distributed data
c Zeros out z component.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd,
c output: ffb
c for isign /= 0, input: q,ffb,isign,nx,ny,kstrt,nyv,kxp2,j2blok,nyhd,
c output: fxy,we
c approximate flop count is: 20*nxc*nyc + 6*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c fxy(3,k,j,l) = zero,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated.
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c real(ffb(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyvh = first dimension of field arrays, must be >= nyh
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffb, zero, zt1
      dimension q(nyvh,kxp2,j2blok), fxy(3,nyvh,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      nyh = ny/2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffb(k,j,l) = cmplx(affp,1.)
      else
         ffb(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         at1 = real(ffb(k,j,l))*aimag(ffb(k,j,l))
         at2 = dkx*at1
         at3 = dny*float(k - 1)*at1
         zt1 = cmplx(aimag(q(k,j,l)),-real(q(k,j,l)))
         fxy(1,k,j,l) = at2*zt1
         fxy(2,k,j,l) = at3*zt1
         fxy(3,k,j,l) = zero
         wp = wp + 2.0*at1*q(k,j,l)*conjg(q(k,j,l))
   50    continue
c mode numbers ky = 0, ny/2
         at1 = real(ffb(1,j,l))*aimag(ffb(1,j,l))
         fxy(1,1,j,l) = cmplx(dkx*at1*aimag(q(1,j,l)),0.)
         fxy(2,1,j,l) = zero
         fxy(3,1,j,l) = zero
         wp = wp + at1*aimag(q(1,j,l))**2
      endif
   60 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 70 k = 1, nyh
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = zero
         fxy(3,k,1,l) = zero
   70    continue
      endif
   80 continue
   90 continue
      we = float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPMX2(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with mixed dirichlet/periodic boundary conditions for distributed data
c input: all, output: cu
c approximate flop count is: 29*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0.
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp2 = number of data values per block
c j2blok = number of data blocks
      complex cu, zero, zt1, zt2
      dimension cu(3,nyv,kxp2,j2blok)
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         zt1 = at1*(dkx*cu(1,k,j,l) + dky*cu(2,k,j,l))
         cu(1,k,j,l) = cu(1,k,j,l) - dkx*zt1
         cu(2,k,j,l) = cu(2,k,j,l) - dky*zt1
         zt2 = conjg(zt1)
         cu(1,k1,j,l) = cu(1,k1,j,l) - dkx*zt2
         cu(2,k1,j,l) = cu(2,k1,j,l) + dky*zt2
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         cu(1,1,j,l) = zero
         cu(1,k1,j,l) = zero
         cu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         cu(2,k,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         cu(1,1,1,l) = zero
         cu(2,1,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPM2(cu,nx,ny,kstrt,nyvh,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with mixed dirichlet/periodic boundary conditions for distributed data
c input: all, output: cu
c approximate flop count is: 20*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0.
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyvh = first dimension of field arrays, must be >= nyh
c kxp2 = number of data values per block
c j2blok = number of data blocks
      complex cu, zero, zt1
      dimension cu(3,nyvh,kxp2,j2blok)
      nyh = ny/2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         zt1 = at1*(dkx*cu(1,k,j,l) + dky*cu(2,k,j,l))
         cu(1,k,j,l) = cu(1,k,j,l) - dkx*zt1
         cu(2,k,j,l) = cu(2,k,j,l) - dky*zt1
   10    continue
c mode numbers ky = 0, ny/2
         cu(1,1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         cu(2,k,1,l) = zero
   30    continue
         cu(1,1,1,l) = zero
         cu(2,1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPMX22(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with mixed dirichlet/periodic boundary conditions for distributed data
c input: all, output: cu
c approximate flop count is: 29*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0.
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c kxp2 = number of data values per block
c j2blok = number of data blocks
      complex cu, zero, zt1, zt2
      dimension cu(2,nyv,kxp2,j2blok)
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         zt1 = at1*(dkx*cu(1,k,j,l) + dky*cu(2,k,j,l))
         cu(1,k,j,l) = cu(1,k,j,l) - dkx*zt1
         cu(2,k,j,l) = cu(2,k,j,l) - dky*zt1
         zt2 = conjg(zt1)
         cu(1,k1,j,l) = cu(1,k1,j,l) - dkx*zt2
         cu(2,k1,j,l) = cu(2,k1,j,l) + dky*zt2
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         cu(1,1,j,l) = zero
         cu(1,k1,j,l) = zero
         cu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         cu(2,k,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         cu(1,1,1,l) = zero
         cu(2,1,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISMX23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,kstr
     1t,nyv,kxp2,j2blok,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed dirichlet/periodic boundary conditions for distributed data
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd
c output: ffb
c for isign = -1, input: cu,ffb,isign,ci,nx,ny,kstrt,nyv,kxp2,j2blok,
c                        nyhd
c output: bxy,wm
c approximate flop count is: 52*nxc*nyc + 11*(nxc + nyc)
c for isign = 1, input: cu,ffb,isign,ci,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd
c output: bxy,wm
c approximate flop count is: 43*nxc*nyc + 10*(nxc + nyc)
c for isign = 2, input: cu,ffb,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: bxy
c approximate flop count is: 13*nxc*nyc + 3*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0.
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c bz(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c bz(kx,ky) = cuz(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c real(ffb(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt2, zt3
      dimension cu(3,nyv,kxp2,j2blok)
      dimension bxy(3,nyv,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffb(k,j,l) = cmplx(affp,1.)
      else
         ffb(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,j,l))*aimag(ffb(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(cu(3,k,j,l)),real(cu(3,k,j,l)))
         zt2 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt3 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         zt3 = at3*zt2 - at2*zt3
         zt2 = -at3*zt1
         zt1 = at2*zt1
         bxy(1,k,j,l) = zt1
         bxy(2,k,j,l) = zt2
         bxy(3,k,j,l) = zt3
         bxy(1,k1,j,l) = -conjg(zt1)
         bxy(2,k1,j,l) = conjg(zt2)
         bxy(3,k1,j,l) = conjg(zt3)
         wp = wp + 2.0*at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)
     1*conjg(cu(2,k,j,l)) + cu(3,k,j,l)*conjg(cu(3,k,j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = ci2*real(ffb(1,j,l))*aimag(ffb(1,j,l))
         at2 = dkx*at1
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(at2*aimag(cu(3,1,j,l)),0.)
         bxy(3,1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + aimag(cu(3,1,j,l))**2)
      endif
   60 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,1,l))*aimag(ffb(k,1,l))
         at2 = dky*at1
         zt2 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         bxy(1,k,1,l) = zero
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = -at2*zt2
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*cu(1,k,1,l)*conjg(cu(1,k,1,l))
   70    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
   80 continue
   90 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffb(k,j,l))
         at1 = at2*aimag(ffb(k,j,l))
         zt1 = at2*cu(1,k,j,l)
         zt2 = at2*cu(2,k,j,l)
         zt3 = at2*cu(3,k,j,l)
         bxy(1,k,j,l) = zt1
         bxy(2,k,j,l) = zt2
         bxy(3,k,j,l) = zt3
         bxy(1,k1,j,l) = conjg(zt1)
         bxy(2,k1,j,l) = -conjg(zt2)
         bxy(3,k1,j,l) = -conjg(zt3)
         wp = wp + 2.0*at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)
     1*conjg(cu(2,k,j,l)) + cu(3,k,j,l)*conjg(cu(3,k,j,l)))
  110    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = ci2*real(ffb(1,j,l))
         at1 = at2*aimag(ffb(1,j,l))
         bxy(1,1,j,l) = cmplx(at2*real(cu(1,1,j,l)),0.)
         bxy(2,1,j,l) = cmplx(0.,at2*aimag(cu(2,1,j,l)))
         bxy(3,1,j,l) = cmplx(0.,at2*aimag(cu(3,1,j,l)))
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + aimag(cu(3,1,j,l))**2)
      endif
  120 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 130 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffb(k,1,l))
         at1 = at2*aimag(ffb(k,1,l))
         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*cu(1,k,1,l)*conjg(cu(1,k,1,l))
  130    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
  140 continue
  150 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffb(k,j,l))
         zt1 = at1*cu(1,k,j,l)
         zt2 = at1*cu(2,k,j,l)
         zt3 = at1*cu(3,k,j,l)
         bxy(1,k,j,l) = zt1
         bxy(2,k,j,l) = zt2
         bxy(3,k,j,l) = zt3
         bxy(1,k1,j,l) = conjg(zt1)
         bxy(2,k1,j,l) = -conjg(zt2)
         bxy(3,k1,j,l) = -conjg(zt3)
  170    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffb(1,j,l))
         bxy(1,1,j,l) = cmplx(at1*real(cu(1,1,j,l)),0.)
         bxy(2,1,j,l) = cmplx(0.,at1*aimag(cu(2,1,j,l)))
         bxy(3,1,j,l) = cmplx(0.,at1*aimag(cu(3,1,j,l)))
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 190 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffb(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
  190    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISM23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,kstrt
     1,nyvh,kxp2,j2blok,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed dirichlet/periodic boundary conditions for distributed data
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd
c output: ffb
c for isign = -1, input: cu,ffb,isign,ci,nx,ny,kstrt,nyv,kxp2,j2blok,
c                        nyhd
c output: bxy,wm
c approximate flop count is: 49*nxc*nyc + 11*(nxc + nyc)
c for isign = 1, input: cu,ffb,isign,ci,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd
c output: bxy,wm
c approximate flop count is: 36*nxc*nyc + 10*(nxc + nyc)
c for isign = 2, input: cu,ffb,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: bxy
c approximate flop count is: 6*nxc*nyc + 3*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0.
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c bz(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c bz(kx,ky) = cuz(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c real(ffb(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyvh = first dimension of field arrays, must be >= nyh
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt2, zt3
      dimension cu(3,nyvh,kxp2,j2blok)
      dimension bxy(3,nyvh,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      nyh = ny/2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffb(k,j,l) = cmplx(affp,1.)
      else
         ffb(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,j,l))*aimag(ffb(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(cu(3,k,j,l)),real(cu(3,k,j,l)))
         zt2 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt3 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         bxy(1,k,j,l) = at2*zt1
         bxy(2,k,j,l) = -at3*zt1
         bxy(3,k,j,l) = at3*zt2 - at2*zt3
         wp = wp + 2.0*at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)
     1*conjg(cu(2,k,j,l)) + cu(3,k,j,l)*conjg(cu(3,k,j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         at1 = ci2*real(ffb(1,j,l))*aimag(ffb(1,j,l))
         at2 = dkx*at1
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(at2*aimag(cu(3,1,j,l)),0.)
         bxy(3,1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + aimag(cu(3,1,j,l))**2)
      endif
   60 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,1,l))*aimag(ffb(k,1,l))
         at2 = dky*at1
         zt2 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         bxy(1,k,1,l) = zero
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = -at2*zt2
         wp = wp + at1*cu(1,k,1,l)*conjg(cu(1,k,1,l))
   70    continue
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
      endif
   80 continue
   90 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, nyh
         at2 = ci2*real(ffb(k,j,l))
         at1 = at2*aimag(ffb(k,j,l))
         bxy(1,k,j,l) = at2*cu(1,k,j,l)
         bxy(2,k,j,l) = at2*cu(2,k,j,l)
         bxy(3,k,j,l) = at2*cu(3,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)
     1*conjg(cu(2,k,j,l)) + cu(3,k,j,l)*conjg(cu(3,k,j,l)))
  110    continue
c mode numbers ky = 0, ny/2
         at2 = ci2*real(ffb(1,j,l))
         at1 = at2*aimag(ffb(1,j,l))
         bxy(1,1,j,l) = cmplx(at2*real(cu(1,1,j,l)),0.)
         bxy(2,1,j,l) = cmplx(0.,at2*aimag(cu(2,1,j,l)))
         bxy(3,1,j,l) = cmplx(0.,at2*aimag(cu(3,1,j,l)))
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + aimag(cu(3,1,j,l))**2)
      endif
  120 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 130 k = 2, nyh
         at2 = ci2*real(ffb(k,1,l))
         at1 = at2*aimag(ffb(k,1,l))
         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = zero
         wp = wp + at1*cu(1,k,1,l)*conjg(cu(1,k,1,l))
  130    continue
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
      endif
  140 continue
  150 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, nyh
         at1 = aimag(ffb(k,j,l))
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
         bxy(3,k,j,l) = at1*cu(3,k,j,l)
  170    continue
c mode numbers ky = 0, ny/2
         at1 = aimag(ffb(1,j,l))
         bxy(1,1,j,l) = cmplx(at1*real(cu(1,1,j,l)),0.)
         bxy(2,1,j,l) = cmplx(0.,at1*aimag(cu(2,1,j,l)))
         bxy(3,1,j,l) = cmplx(0.,at1*aimag(cu(3,1,j,l)))
      endif
  180 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 190 k = 2, nyh
         at1 = aimag(ffb(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = zero
  190    continue
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISMX22(cu,bxy,bz,isign,ffb,ax,ay,affp,ci,wm,nx,ny,k
     1strt,nyv,kxp2,j2blok,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed dirichlet/periodic boundary conditions for distributed data
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd
c output: ffb
c for isign = -1, input: cu,ffb,isign,ci,nx,ny,kstrt,nyv,kxp2,j2blok,
c                        nyhd
c output: bz,wm
c approximate flop count is: 33*nxc*nyc + 8*(nxc + nyc)
c for isign = 1, input: cu,ffb,isign,ci,nx,ny,kstrt,nyv,kxp2,j2blok,
c                       nyhd
c output: bxy,wm
c approximate flop count is: 28*nxc*nyc + 7*(nxc + nyc)
c for isign = 2, input: cu,ffb,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd,
c output: bxy
c approximate flop count is: 8*nxc*nyc + 2*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bz(kx=pi) = 0, bz(ky=pi) = 0.
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c real(ffb(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, bz, ffb, zero, zt1, zt2
      dimension cu(2,nyv,kxp2,j2blok)
      dimension bxy(2,nyv,kxp2,j2blok), bz(nyv,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffb(k,j,l) = cmplx(affp,1.)
      else
         ffb(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 100
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,j,l))*aimag(ffb(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt2 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         zt1 = at3*zt1 - at2*zt2
         bz(k,j,l) = zt1
         bz(k1,j,l) = conjg(zt1)
         wp = wp + 2.0*at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)
     1*conjg(cu(2,k,j,l)))
   50    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = ci2*real(ffb(1,j,l))*aimag(ffb(1,j,l))
         at2 = dkx*at1
         bz(1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         bz(k1,j,l) = zero
         wp = wp + at1*aimag(cu(2,1,j,l))**2
      endif
   60 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 70 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,1,l))*aimag(ffb(k,1,l))
         at2 = dky*at1
         zt2 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         bz(k,1,l) = -at2*zt2
         bz(k1,1,l) = zero
         wp = wp + at1*cu(1,k,1,l)*conjg(cu(1,k,1,l))
   70    continue
         k1 = nyh + 1
         bz(1,1,l) = zero
         bz(k1,1,l) = zero
      endif
   80 continue
   90 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffb(k,j,l))
         at1 = at2*aimag(ffb(k,j,l))
         zt1 = at2*cu(1,k,j,l)
         zt2 = at2*cu(2,k,j,l)
         bxy(1,k,j,l) = zt1
         bxy(2,k,j,l) = zt2
         bxy(1,k1,j,l) = conjg(zt1)
         bxy(2,k1,j,l) = -conjg(zt2)
         wp = wp + 2.0*at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)
     1*conjg(cu(2,k,j,l)))
  110    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = ci2*real(ffb(1,j,l))
         at1 = at2*aimag(ffb(1,j,l))
         bxy(1,1,j,l) = cmplx(at2*real(cu(1,1,j,l)),0.)
         bxy(2,1,j,l) = cmplx(0.,at2*aimag(cu(2,1,j,l)))
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         wp = wp + at1*aimag(cu(2,1,j,l))**2
      endif
  120 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 130 k = 2, nyh
         k1 = ny2 - k
         at2 = ci2*real(ffb(k,1,l))
         at1 = at2*aimag(ffb(k,1,l))
         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(2,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         wp = wp + at1*cu(1,k,1,l)*conjg(cu(1,k,1,l))
  130    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
      endif
  140 continue
  150 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffb(k,j,l))
         zt1 = at1*cu(1,k,j,l)
         zt2 = at1*cu(2,k,j,l)
         bxy(1,k,j,l) = zt1
         bxy(2,k,j,l) = zt2
         bxy(1,k1,j,l) = conjg(zt1)
         bxy(2,k1,j,l) = -conjg(zt2)
  170    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = aimag(ffb(1,j,l))
         bxy(1,1,j,l) = cmplx(at1*real(cu(1,1,j,l)),0.)
         bxy(2,1,j,l) = cmplx(0.,at1*aimag(cu(2,1,j,l)))
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 190 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffb(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
  190    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine IPBPOISMX23(cu,bxy,ffb,ci,wm,nx,ny,kstrt,nyv,kxp2,j2blo
     1k,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field with mixed dirichlet/periodic boundary conditions
c for distributed data
c input: cu,ffb,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd, output: bxy,wm
c approximate flop count is: 52*nxc*nyc + 11*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0.
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c real(ffb(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt2, zt3
      dimension cu(3,nyv,kxp2,j2blok)
      dimension bxy(3,nyv,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffb(k,j,l))
         zt1 = cmplx(-aimag(cu(3,k,j,l)),real(cu(3,k,j,l)))
         zt2 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt3 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         zt3 = at3*zt2 - at2*zt3
         zt2 = -at3*zt1
         zt1 = at2*zt1
         bxy(1,k,j,l) = zt1
         bxy(2,k,j,l) = zt2
         bxy(3,k,j,l) = zt3
         bxy(1,k1,j,l) = -conjg(zt1)
         bxy(2,k1,j,l) = conjg(zt2)
         bxy(3,k1,j,l) = conjg(zt3)
         wp = wp + 2.0*at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)
     1*conjg(cu(2,k,j,l)) + cu(3,k,j,l)*conjg(cu(3,k,j,l)))
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at1 = ci2*real(ffb(1,j,l))
         at2 = dkx*at1
         at1 = at1*aimag(ffb(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(at2*aimag(cu(3,1,j,l)),0.)
         bxy(3,1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + aimag(cu(3,1,j,l))**2)
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,1,l))
         at2 = dky*at1
         at1 = at1*aimag(ffb(k,1,l))
         zt2 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         bxy(1,k,1,l) = zero
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = -at2*zt2
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*cu(1,k,1,l)*conjg(cu(1,k,1,l))
   30    continue
         k1 = nyh + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      wm = float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine IPBPOISM23(cu,bxy,ffb,ci,wm,nx,ny,kstrt,nyvh,kxp2,j2blo
     1k,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field with mixed dirichlet/periodic boundary conditions
c for distributed data
c input: cu,ffb,ci,nx,ny,kstrt,nyv,kxp,jblok,nyhd, output: bxy,wm
c approximate flop count is: 49*nxc*nyc + 11*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bz(ky=pi) = 0.
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c real(ffb(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyvh = first dimension of field arrays, must be >= nyh
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt2, zt3
      dimension cu(3,nyvh,kxp2,j2blok)
      dimension bxy(3,nyvh,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      nyh = ny/2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffb(k,j,l))
         zt1 = cmplx(-aimag(cu(3,k,j,l)),real(cu(3,k,j,l)))
         zt2 = cmplx(-aimag(cu(2,k,j,l)),real(cu(2,k,j,l)))
         zt3 = cmplx(-aimag(cu(1,k,j,l)),real(cu(1,k,j,l)))
         bxy(1,k,j,l) = at2*zt1
         bxy(2,k,j,l) = -at3*zt1
         bxy(3,k,j,l) = at3*zt2 - at2*zt3
         wp = wp + 2.0*at1*(cu(1,k,j,l)*conjg(cu(1,k,j,l)) + cu(2,k,j,l)
     1*conjg(cu(2,k,j,l)) + cu(3,k,j,l)*conjg(cu(3,k,j,l)))
   10    continue
c mode numbers ky = 0, ny/2
         at1 = ci2*real(ffb(1,j,l))
         at2 = dkx*at1
         at1 = at1*aimag(ffb(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(at2*aimag(cu(3,1,j,l)),0.)
         bxy(3,1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + aimag(cu(3,1,j,l))**2)
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         dky = dny*float(k - 1)
         at1 = ci2*real(ffb(k,1,l))
         at2 = dky*at1
         at1 = at1*aimag(ffb(k,1,l))
         zt2 = cmplx(-aimag(cu(1,k,1,l)),real(cu(1,k,1,l)))
         bxy(1,k,1,l) = zero
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = -at2*zt2
         wp = wp + at1*cu(1,k,1,l)*conjg(cu(1,k,1,l))
   30    continue
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
      endif
   40 continue
   50 continue
      wm = float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWELMX2(exy,bxy,cu,ffb,affp,ci,dt,wf,wm,nx,ny,kstrt,
     1nyv,kxp2,j2blok,nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields
c with mixed dirichlet/periodic boundary conditions for distributed data
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 158*nxc*nyc + 32*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0.
c and similarly for bx, by, bz.
c cu(i,k,j,l) = i-th component of complex current density and
c exy(i,k,j,l) = i-th component of complex electric field,
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny = system length in x/y direction
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffb
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      dimension exy(3,nyv,kxp2,j2blok), bxy(3,nyv,kxp2,j2blok)
      dimension cu(3,nyv,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      if (ci.le.0.) return
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
      if (kstrt.gt.nx) go to 50
c calculate the electromagnetic fields
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffb(k,j,l))
c update magnetic field half time step, ky > 0
         zt1 = cmplx(-aimag(exy(3,k,j,l)),real(exy(3,k,j,l)))
         zt2 = cmplx(-aimag(exy(2,k,j,l)),real(exy(2,k,j,l)))
         zt3 = cmplx(-aimag(exy(1,k,j,l)),real(exy(1,k,j,l)))
         zt4 = bxy(1,k,j,l) - dth*(dky*zt1)
         zt5 = bxy(2,k,j,l) + dth*(dkx*zt1)
         zt6 = bxy(3,k,j,l) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k,j,l) + cdt*(dky*zt1) - afdt*cu(1,k,j,l)
         zt8 = exy(2,k,j,l) - cdt*(dkx*zt1) - afdt*cu(2,k,j,l)
         zt9 = exy(3,k,j,l) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,k,j,l)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         zt4 = zt4 - dth*(dky*zt1)
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
         ws = ws + 2.0*anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conj
     1g(zt9))
         wp = wp + 2.0*anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conj
     1g(zt6))
         exy(1,k,j,l) = zt7
         exy(2,k,j,l) = zt8
         exy(3,k,j,l) = zt9
         bxy(1,k,j,l) = zt4
         bxy(2,k,j,l) = zt5
         bxy(3,k,j,l) = zt6
c update electric and  magnetic fields, ky < 0
         exy(1,k1,j,l) = conjg(zt7)
         exy(2,k1,j,l) = -conjg(zt8)
         exy(3,k1,j,l) = -conjg(zt9)
         bxy(1,k1,j,l) = -conjg(zt4)
         bxy(2,k1,j,l) = conjg(zt5)
         bxy(3,k1,j,l) = conjg(zt6)
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         afdt = adt*aimag(ffb(1,j,l))
c update magnetic field half time step
         at1 = aimag(exy(3,1,j,l))
         at2 = aimag(exy(2,1,j,l))
         at5 = real(bxy(2,1,j,l)) - dth*(dkx*at1)
         at6 = real(bxy(3,1,j,l)) + dth*(dkx*at2)
c update electric field whole time step
         at8 = at2 - cdt*(dkx*at6) - afdt*aimag(cu(2,1,j,l))
         at9 = at1 + cdt*(dkx*at5) - afdt*aimag(cu(3,1,j,l))
c update magnetic field half time step and store electric field
         at5 = at5 - dth*(dkx*at9)
         at6 = at6 + dth*(dkx*at8)
         ws = ws + anorm*(at8*at8 + at9*at9)
         wp = wp + anorm*(at5*at5 + at6*at6)
         exy(1,1,j,l) = zero
         exy(2,1,j,l) = cmplx(0.,at8)
         exy(3,1,j,l) = cmplx(0.,at9)
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(at5,0.)
         bxy(3,1,j,l) = cmplx(at6,0.)
         exy(1,k1,j,l) = zero
         exy(2,k1,j,l) = zero
         exy(3,k1,j,l) = zero
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffb(k,1,l))
c update magnetic field half time step
         zt3 = cmplx(-aimag(exy(1,k,1,l)),real(exy(1,k,1,l)))
         zt6 = bxy(3,k,1,l) + dth*(dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt7 = exy(1,k,1,l) + cdt*(dky*zt1) - afdt*cu(1,k,1,l)
c update magnetic field half time step and store electric field
         zt3 = cmplx(-aimag(zt7),real(zt7))
         zt6 = zt6 + dth*(dky*zt3)
         ws = ws + anorm*zt7*conjg(zt7)
         wp = wp + anorm*zt6*conjg(zt6)
         exy(1,k,1,l) = zt7
         exy(2,k,1,l) = zero
         exy(3,k,1,l) = zero
         bxy(1,k,1,l) = zero
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = zt6
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         exy(1,1,1,l) = zero
         exy(2,1,1,l) = zero
         exy(3,1,1,l) = zero
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWELM2(exy,bxy,cu,ffb,affp,ci,dt,wf,wm,nx,ny,kstrt,n
     1yvh,kxp2,j2blok,nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields
c with mixed dirichlet/periodic boundary conditions for distributed data
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 146*nxc*nyc + 32*(nxc + nyc)
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ez(ky=pi) = 0.
c and similarly for bx, by, bz.
c cu(i,k,j,l) = i-th component of complex current density and
c exy(i,k,j,l) = i-th component of complex electric field,
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c aimag(ffb(k,j,l)) = finite-size particle shape factor s
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny = system length in x/y direction
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c nyvh = first dimension of field arrays, must be >= nyh
c nyhd = first dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffb
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      dimension exy(3,nyvh,kxp2,j2blok), bxy(3,nyvh,kxp2,j2blok)
      dimension cu(3,nyvh,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      if (ci.le.0.) return
      nyh = ny/2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
      if (kstrt.gt.nx) go to 50
c calculate the electromagnetic fields
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffb(k,j,l))
c update magnetic field half time step, ky > 0
         zt1 = cmplx(-aimag(exy(3,k,j,l)),real(exy(3,k,j,l)))
         zt2 = cmplx(-aimag(exy(2,k,j,l)),real(exy(2,k,j,l)))
         zt3 = cmplx(-aimag(exy(1,k,j,l)),real(exy(1,k,j,l)))
         zt4 = bxy(1,k,j,l) - dth*(dky*zt1)
         zt5 = bxy(2,k,j,l) + dth*(dkx*zt1)
         zt6 = bxy(3,k,j,l) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt2 = cmplx(-aimag(zt5),real(zt5))
         zt3 = cmplx(-aimag(zt4),real(zt4))
         zt7 = exy(1,k,j,l) + cdt*(dky*zt1) - afdt*cu(1,k,j,l)
         zt8 = exy(2,k,j,l) - cdt*(dkx*zt1) - afdt*cu(2,k,j,l)
         zt9 = exy(3,k,j,l) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,k,j,l)
c update magnetic field half time step and store electric field
         zt1 = cmplx(-aimag(zt9),real(zt9))
         zt2 = cmplx(-aimag(zt8),real(zt8))
         zt3 = cmplx(-aimag(zt7),real(zt7))
         zt4 = zt4 - dth*(dky*zt1)
         zt5 = zt5 + dth*(dkx*zt1)
         zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
         ws = ws + 2.0*anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conj
     1g(zt9))
         wp = wp + 2.0*anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conj
     1g(zt6))
         exy(1,k,j,l) = zt7
         exy(2,k,j,l) = zt8
         exy(3,k,j,l) = zt9
         bxy(1,k,j,l) = zt4
         bxy(2,k,j,l) = zt5
         bxy(3,k,j,l) = zt6
   10    continue
c mode numbers ky = 0, ny/2
         afdt = adt*aimag(ffb(1,j,l))
c update magnetic field half time step
         at1 = aimag(exy(3,1,j,l))
         at2 = aimag(exy(2,1,j,l))
         at5 = real(bxy(2,1,j,l)) - dth*(dkx*at1)
         at6 = real(bxy(3,1,j,l)) + dth*(dkx*at2)
c update electric field whole time step
         at8 = at2 - cdt*(dkx*at6) - afdt*aimag(cu(2,1,j,l))
         at9 = at1 + cdt*(dkx*at5) - afdt*aimag(cu(3,1,j,l))
c update magnetic field half time step and store electric field
         at5 = at5 - dth*(dkx*at9)
         at6 = at6 + dth*(dkx*at8)
         ws = ws + anorm*(at8*at8 + at9*at9)
         wp = wp + anorm*(at5*at5 + at6*at6)
         exy(1,1,j,l) = zero
         exy(2,1,j,l) = cmplx(0.,at8)
         exy(3,1,j,l) = cmplx(0.,at9)
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(at5,0.)
         bxy(3,1,j,l) = cmplx(at6,0.)
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffb(k,1,l))
c update magnetic field half time step
         zt3 = cmplx(-aimag(exy(1,k,1,l)),real(exy(1,k,1,l)))
         zt6 = bxy(3,k,1,l) + dth*(dky*zt3)
c update electric field whole time step
         zt1 = cmplx(-aimag(zt6),real(zt6))
         zt7 = exy(1,k,1,l) + cdt*(dky*zt1) - afdt*cu(1,k,1,l)
c update magnetic field half time step and store electric field
         zt3 = cmplx(-aimag(zt7),real(zt7))
         zt6 = zt6 + dth*(dky*zt3)
         ws = ws + anorm*zt7*conjg(zt7)
         wp = wp + anorm*zt6*conjg(zt6)
         exy(1,k,1,l) = zt7
         exy(2,k,1,l) = zero
         exy(3,k,1,l) = zero
         bxy(1,k,1,l) = zero
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = zt6
   30    continue
         exy(1,1,1,l) = zero
         exy(2,1,1,l) = zero
         exy(3,1,1,l) = zero
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
      endif
   40 continue
   50 continue
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PDMFIELDM2(q1,q,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
c this subroutine copies the charge density into a smaller array
      implicit none
      integer nx, ny, kstrt, nyv, nyvh, kxp2, j2blok
      complex q1, q
      dimension q1(nyv,kxp2,j2blok), q(nyvh,kxp2,j2blok)
      integer j, k, l, nyh
      if (kstrt.gt.nx) return
      nyh = ny/2
      do 30 l = 1, j2blok
      do 20 j = 1, kxp2
      do 10 k = 1, nyh
      q(k,j,l) = q1(k,j,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCMFIELDM2(cu1,cu,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
c this subroutine copies the current into a smaller array
      implicit none
      integer nx, ny, kstrt, nyv, nyvh, kxp2, j2blok
      complex cu1, cu
      dimension cu1(3,nyv,kxp2,j2blok), cu(3,nyvh,kxp2,j2blok)
      integer j, k, l, nyh
      if (kstrt.gt.nx) return
      nyh = ny/2
      do 30 l = 1, j2blok
      do 20 j = 1, kxp2
      do 10 k = 1, nyh
      cu(1,k,j,l) = cu1(1,k,j,l)
      cu(2,k,j,l) = cu1(2,k,j,l)
      cu(3,k,j,l) = cu1(3,k,j,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PEMFIELDM2(fxy,exy,ffb,isign,nx,ny,kstrt,nyv,nyvh,kxp2,
     1j2blok,nyhd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign <= 0
c adds image charges appropriate for electric field if isign >= 0
c or appropriate for magnetic field if isign < 0
c includes additional smoothing for isign /= 0
      implicit none
      integer isign, nx, ny, kstrt, nyv, nyvh, kxp2, j2blok, nyhd
      complex fxy, exy, ffb
      dimension fxy(3,nyv,kxp2,j2blok), exy(3,nyvh,kxp2,j2blok)
      dimension ffb(nyhd,kxp2,j2blok)
      complex zero, zt1, zt2, zt3
      integer j, k, l, nyh, ny2, ks, joff, k1
      real at1
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      zero = cmplx(0.,0.)
      if (kstrt.gt.nx) return
c add the fields
      if (isign.gt.0) then
         do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
         joff = kxp2*(l + ks) - 1
         do 20 j = 1, kxp2
         if ((j+joff).gt.0) then
            do 10 k = 2, nyh
            k1 = ny2 - k
            at1 = aimag(ffb(k,j,l))
            zt1 = exy(1,k,j,l)*at1
            zt2 = exy(2,k,j,l)*at1
            zt3 = exy(3,k,j,l)*at1
            fxy(1,k,j,l) = fxy(1,k,j,l) + zt1
            fxy(2,k,j,l) = fxy(2,k,j,l) + zt2
            fxy(3,k,j,l) = fxy(3,k,j,l) + zt3
            fxy(1,k1,j,l) = fxy(1,k1,j,l) + conjg(zt1)
            fxy(2,k1,j,l) = fxy(2,k1,j,l) - conjg(zt2)
            fxy(3,k1,j,l) = fxy(3,k1,j,l) - conjg(zt3)
   10       continue
c mode numbers ky = 0, ny/2
            at1 = aimag(ffb(1,j,l))
            fxy(1,1,j,l) = fxy(1,1,j,l) + exy(1,1,j,l)*at1
            fxy(2,1,j,l) = fxy(2,1,j,l) + exy(2,1,j,l)*at1
            fxy(3,1,j,l) = fxy(3,1,j,l) + exy(3,1,j,l)*at1
         endif
   20    continue
c mode numbers kx = 0, nx
         if ((l+ks).eq.0) then
            do 30 k = 1, nyh
            at1 = aimag(ffb(k,1,l))
            fxy(1,k,1,l) = fxy(1,k,1,l) + exy(1,k,1,l)*at1
            fxy(2,k,1,l) = fxy(2,k,1,l) + exy(2,k,1,l)*at1
            fxy(3,k,1,l) = fxy(3,k,1,l) + exy(3,k,1,l)*at1
   30       continue
         endif
   40    continue
c copy the magnetic fields
      else if (isign.lt.0) then
         do 80 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
         joff = kxp2*(l + ks) - 1
         do 60 j = 1, kxp2
         if ((j+joff).gt.0) then
            do 50 k = 2, nyh
            k1 = ny2 - k
            at1 = aimag(ffb(k,j,l))
            zt1 = exy(1,k,j,l)*at1
            zt2 = exy(2,k,j,l)*at1
            zt3 = exy(3,k,j,l)*at1
            fxy(1,k,j,l) = zt1
            fxy(2,k,j,l) = zt2
            fxy(3,k,j,l) = zt3
            fxy(1,k1,j,l) = -conjg(zt1)
            fxy(2,k1,j,l) = conjg(zt2)
            fxy(3,k1,j,l) = conjg(zt3)
   50       continue
c mode numbers ky = 0, ny/2
            k1 = nyh + 1
            at1 = aimag(ffb(1,j,l))
            fxy(1,1,j,l) = exy(1,1,j,l)*at1
            fxy(2,1,j,l) = exy(2,1,j,l)*at1
            fxy(3,1,j,l) = exy(3,1,j,l)*at1
            fxy(1,k1,j,l) = zero
            fxy(2,k1,j,l) = zero
            fxy(3,k1,j,l) = zero
         endif
   60    continue
c mode numbers kx = 0, nx
         if ((l+ks).eq.0) then
            do 70 k = 2, nyh
            k1 = ny2 - k
            at1 = aimag(ffb(k,1,l))
            fxy(1,k,1,l) = exy(1,k,1,l)*at1
            fxy(2,k,1,l) = exy(2,k,1,l)*at1
            fxy(3,k,1,l) = exy(3,k,1,l)*at1
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
   70       continue
            k1 = nyh + 1
            at1 = aimag(ffb(1,1,l))
            fxy(1,1,1,l) = exy(1,1,1,l)*at1
            fxy(2,1,1,l) = exy(2,1,1,l)*at1
            fxy(3,1,1,l) = exy(3,1,1,l)*at1
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
         endif
   80    continue
c copy the electric fields
      else
         do 120 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
         joff = kxp2*(l + ks) - 1
         do 100 j = 1, kxp2
         if ((j+joff).gt.0) then
            do 90 k = 2, nyh
            k1 = ny2 - k
            zt1 = exy(1,k,j,l)
            zt2 = exy(2,k,j,l)
            zt3 = exy(3,k,j,l)
            fxy(1,k,j,l) = zt1
            fxy(2,k,j,l) = zt2
            fxy(3,k,j,l) = zt3
            fxy(1,k1,j,l) = conjg(zt1)
            fxy(2,k1,j,l) = -conjg(zt2)
            fxy(3,k1,j,l) = -conjg(zt3)
   90       continue
c mode numbers ky = 0, ny/2
            k1 = nyh + 1
            fxy(1,1,j,l) = exy(1,1,j,l)
            fxy(2,1,j,l) = exy(2,1,j,l)
            fxy(3,1,j,l) = exy(3,1,j,l)
            fxy(1,k1,j,l) = zero
            fxy(2,k1,j,l) = zero
            fxy(3,k1,j,l) = zero
         endif
  100    continue
c mode numbers kx = 0, nx
         if ((l+ks).eq.0) then
            do 110 k = 2, nyh
            k1 = ny2 - k
            fxy(1,k,1,l) = exy(1,k,1,l)
            fxy(2,k,1,l) = exy(2,k,1,l)
            fxy(3,k,1,l) = exy(3,k,1,l)
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
  110       continue
            k1 = nyh + 1
            fxy(1,1,1,l) = exy(1,1,1,l)
            fxy(2,1,1,l) = exy(2,1,1,l)
            fxy(3,1,1,l) = exy(3,1,1,l)
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
         endif
  120    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMFIELDM2(pot1,pot,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
c copies image charges appropriate for potential
      implicit none
      integer nx, ny, kstrt, nyv, nyvh, kxp2, j2blok
      complex pot1, pot
      dimension pot1(nyv,kxp2,j2blok), pot(nyvh,kxp2,j2blok)
      complex zero, zt1
      integer j, k, l, nyh, ny2, ks, joff, k1
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      zero = cmplx(0.,0.)
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         zt1 = pot(k,j,l)
         pot1(k,j,l) = zt1
         pot1(k1,j,l) = -conjg(zt1)
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         pot1(1,j,l) = pot(1,j,l)
         pot1(k1,j,l) = zero
         endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         pot1(k,1,l) = pot(k,1,l)
         pot1(k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         pot1(1,1,l) = pot(1,1,l)
         pot1(k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCPFIELDM2(fxy,exy,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok)
c this subroutine copies complex vector fields and adds image charges
c appropriate for electric field
      implicit none
      integer nx, ny, kstrt, nyv, nyvh, kxp2, j2blok
      complex fxy, exy
      dimension fxy(3,nyv,kxp2,j2blok), exy(3,nyvh,kxp2,j2blok)
c local data
      integer isign, nyhd
      complex ffb
      dimension ffb(1,1,1)
      isign = 0
      nyhd = 1
      call PEMFIELDM2(fxy,exy,ffb,isign,nx,ny,kstrt,nyv,nyvh,kxp2,j2blok
     1,nyhd)
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVPOTMX23(bxy,axy,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with mixed dirichlet/periodic boundary conditions
c for distributed data.
c input: bxy,nx,ny,kstrt,nyv,kxp2,j2blok, output: axy
c approximate flop count is: 28*nxc*nyc + 3*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0.
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c axy(i,k,j,l) = i-th component of complex vector potential,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny
      complex bxy, axy, zero, zt1, zt2, zt3
      dimension bxy(3,nyv,kxp2,j2blok), axy(3,nyv,kxp2,j2blok)
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate vector potential
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(bxy(3,k,j,l)),real(bxy(3,k,j,l)))
         zt2 = cmplx(-aimag(bxy(2,k,j,l)),real(bxy(2,k,j,l)))
         zt3 = cmplx(-aimag(bxy(1,k,j,l)),real(bxy(1,k,j,l)))
         zt3 = at3*zt2 - at2*zt3
         zt2 = -at3*zt1
         zt1 = at2*zt1
         axy(1,k,j,l) = zt1
         axy(2,k,j,l) = zt2
         axy(3,k,j,l) = zt3
         axy(1,k1,j,l) = conjg(zt1)
         axy(2,k1,j,l) = -conjg(zt2)
         axy(3,k1,j,l) = -conjg(zt3)
   10    continue
c mode numbers ky = 0, ny/2
         k1 = nyh + 1
         at2 = 1.0/dkx
         axy(1,1,j,l) = zero
         axy(2,1,j,l) = cmplx(0.,-at2*real(bxy(3,1,j,l)))
         axy(3,1,j,l) = cmplx(0.,at2*real(bxy(2,1,j,l)))
         axy(1,k1,j,l) = zero
         axy(2,k1,j,l) = zero
         axy(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at2 = 1.0/dky
         zt1 = cmplx(-aimag(bxy(3,k,1,l)),real(bxy(3,k,1,l)))
         axy(1,k,1,l) = at2*zt1
         axy(2,k,1,l) = zero
         axy(3,k,1,l) = zero
         axy(1,k1,1,l) = zero
         axy(2,k1,1,l) = zero
         axy(3,k1,1,l) = zero
   30    continue
         k1 = nyh + 1
         axy(1,1,1,l) = zero
         axy(2,1,1,l) = zero
         axy(3,1,1,l) = zero
         axy(1,k1,1,l) = zero
         axy(2,k1,1,l) = zero
         axy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVPOTM23(bxy,axy,nx,ny,kstrt,nyvh,kxp2,j2blok)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with mixed dirichlet/periodic boundary conditions
c for distributed data.
c input: bxy,nx,ny,kstrt,nyv,kxp2,j2blok, output: axy
c approximate flop count is: 21*nxc*nyc + 3*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0.
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c axy(i,k,j,l) = i-th component of complex vector potential,
c for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c nx/ny = system length in x/y direction
c nyvh = first dimension of field arrays, must be >= nyh
      complex bxy, axy, zero, zt1, zt2, zt3
      dimension bxy(3,nyvh,kxp2,j2blok), axy(3,nyvh,kxp2,j2blok)
      nyh = ny/2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate vector potential
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, nyh
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = dky*at1
         at3 = dkx*at1
         zt1 = cmplx(-aimag(bxy(3,k,j,l)),real(bxy(3,k,j,l)))
         zt2 = cmplx(-aimag(bxy(2,k,j,l)),real(bxy(2,k,j,l)))
         zt3 = cmplx(-aimag(bxy(1,k,j,l)),real(bxy(1,k,j,l)))
         axy(1,k,j,l) = at2*zt1
         axy(2,k,j,l) = -at3*zt1
         axy(3,k,j,l) = at3*zt2 - at2*zt3
   10    continue
c mode numbers ky = 0, ny/2
         at2 = 1.0/dkx
         axy(1,1,j,l) = zero
         axy(2,1,j,l) = cmplx(0.,-at2*real(bxy(3,1,j,l)))
         axy(3,1,j,l) = cmplx(0.,at2*real(bxy(2,1,j,l)))
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, nyh
         dky = dny*float(k - 1)
         at2 = 1.0/dky
         zt1 = cmplx(-aimag(bxy(3,k,1,l)),real(bxy(3,k,1,l)))
         axy(1,k,1,l) = at2*zt1
         axy(2,k,1,l) = zero
         axy(3,k,1,l) = zero
   30    continue
         axy(1,1,1,l) = zero
         axy(2,1,1,l) = zero
         axy(3,1,1,l) = zero
      endif
   40 continue
   50 continue
      return
      end
