c 2d parallel PIC library for solving field equations with non-periodic
c boundary conditions
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: may 14, 2004
c-----------------------------------------------------------------------
      subroutine PBNDRYV2(q,ffc,bv,nx,ny,kstrt,nyv,kxp,jblok,nyhd)
c this subroutine calculates the boundary values of electric field of
c the periodic solution of poisson's equation in Fourier space from the
c charge density.  The results are used in calculating the solution of a
c laplacian in order to satisfy non-periodic boundary conditions.
c algorithm used in described in V. K. Decyk and J. M. Dawson,
c Journal of Computational Physics 30, 407 (1979).
c input: q, ffp, nx, ny, nxvh, nxhd, nyhd, output: bv
c approximate flop count = 16*nxc*nyc + 7*nyc + 7*nxc
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c q(j,k) = input complex charge density,
c real(ffc(j,k)) = finite-size particle shape factor s,
c aimag(ffc(j,k)) = potential green's function g,
c all for for fourier mode (j-1,k-1)
c bv = boundary fields, bv(k,3) = KmPm and bv(k,4) = PIm, except
c imag(PI0) = net charge density rho
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nxhd = must be >= nx/2
c nyhd = must be >= ny/2
      implicit none
      complex q, ffc, bv
      integer nx, ny, kstrt, nyv, kxp, jblok, nyhd
      dimension q(nyv,kxp,jblok), ffc(nyhd,kxp,jblok), bv(nyhd,4)
c local data
      integer nxh, nyh, ny2, j, k, k1, l, ks, joff
      real dnx, dny, dky, sum1, sum2, sum3, sum4, at1, at2
      nxh = nx/2
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      if (kstrt.gt.nxh) return
      do 40 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
c calculate KmPm and PIm
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
      do 10 j = 1, kxp
      if ((j+joff).gt.0) then
         at2 = real(ffc(k,j,l))
         at1 = dky*at2
         at2 = dnx*float(j + joff)*at2
         sum1 = sum1 + at1*real(q(k,j,l) + q(k1,j,l))
         sum2 = sum2 + at1*aimag(q(k,j,l) - q(k1,j,l))
         sum3 = sum3 + at2*aimag(q(k,j,l) + q(k1,j,l))
         sum4 = sum4 + at2*real(q(k1,j,l) - q(k,j,l))
      endif
   10 continue
      bv(k,3) = cmplx(sum1,sum2)
      bv(k,4) = cmplx(sum3,sum4)
      if ((l+ks).eq.0) bv(k,3) = bv(k,3) + (dky*real(ffc(k,1,l)))*q(k,1,
     1l)
   20 continue
c calculate P0 and PI0
      sum1 = 0.
      sum2 = 0.
      do 30 j = 1, kxp
      if ((j+joff).gt.0) then
         at1 = real(ffc(1,j,l))
         sum1 = sum1 + at1*real(q(1,j,l))
         sum2 = sum2 + (dnx*float(j + joff)*at1)*aimag(q(1,j,l))
      endif
   30 continue
      bv(1,3) = cmplx(sum1 + sum1,0.)
c imaginary part of bv(1,4) contains net charge rho00
      bv(1,4) = cmplx(sum2 + sum2,0.)
      if ((l+ks).eq.0) bv(1,4) = bv(1,4) + cmplx(0.,real(ffc(1,1,l))*rea
     1l(q(1,1,l)))
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISB2(fx,fy,isign,ffb,bv,bcd,mixup,sct,t,we,affp,indx
     1,ny,kstrt,nyv,kxp,jblok,nxhd,nyhd)
c not yet complete !!!!
c this subroutine finds corrections to 2d poisson's equation for 
c force/charge or potential with vacuum boundary conditions and 
c external surface charge, for distributed data.
c average potential across system is zero.  a periodic solution is
c assumed to have been found first with ppoisp2, and boundary values
c with bndryv2
c algorithm used in described in V. K. Decyk and J. M. Dawson,
c Journal of Computational Physics 30, 407 (1979).
c for isign = 0, input: isign,indx,nxh,ny,nyh,nxvh
c                output: ffb,bcd
c                scratch: mixup,sct,t
c for isign = -1, input:  isign,fx,fy,ffb,bv,bcd,affp,indx,ny
c                         nxvh,nxhd,nyhd
c                 output: fx, fy, bv, we
c approximate flop count is: 54*nxc*nyc + + 115*nyc + 3*nxc
c for isign = 1, input:  isign,fx,ffb,bv,bcd,affp,indx,ny,nxvh,nxhd,nyhd
c                output: fx, bv, we
c approximate flop count is: 27*nxc*nyc + + 74*nyc + 6*nxc
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign < 0, the force/charge correction is calculated:
c exc(km,x) = -Am*exp(-km*(Lx-x)) + Bm*exp(-km*x)
c eyc(km,x) = -sqrt(-1)*(Am*exp(-km*(Lx-x)) + Bm*exp(-km*x))
c exc(k0,x) = -4*pi*rho00*(Lx/2-x) - A0
c eyc(k0,x) = 0.
c where Am = .5*(4*pi*sigma(x=Lx,k) + PIm - km*Pm),
c       Bm = .5*(4*pi*sigma(x=0,k) - PIm - km*Pm),
c and   A0 = 2*pi*sigma(x=Lx) - 2*pi*sigma(x=0) + PI0
c where PIm and Pm are the periodic ex and phi at the boundaries
c the calculations are done in fourier space and are added to the
c periodic forces already in fx, fy
c on output, bv = value of electric fields on right boundary:
c bv(k,5) = ex(x=Lx), bv(k,6) = ey(x=Lx)
c if isign = 1, potential correction is calculated:
c potc(km,x) = (Am*exp(-km*(Lx-x)) + Bm*exp(-km*x)/km
c potc(k0,x) = 2*pi*rho00*x*(Lx-x) - A0*(Lx/2-x) - P0
c the calculation is done in fourier space and is added to the
c periodic potential already in fx.
c on output, bv = value of potential on right boundary:
c bv(k,6) = phi(x=Lx)
c if isign = 0, form factor arrays ffb and bcd are prepared
c on input, fx and/or fy contain periodic part of solution
c on output, fx and/or fy contain total solution
c ffb(j,k) = (1/nx)*inverse fft(exp(-dky*float(nx + 1 - j))))
c real(ffb(j,1)) = (1/nx)*inverse fft((j - 1)*(nx + 1 - j))))
c aimag(ffb(j,1)) = (1/nx)*inverse fft((nx/2 + 1 - j)))
c on input, bv = input surface charge and boundary values
c for fourier mode k-1:
c bv(k,1) = 4*pi*sigma(x=0), bv(k,2) = 4*pi*sigma(x=Lx)
c bv(k,3) = KmPm, bv(k,4) = PIm
c both are normalized in the same way as the electric field.
c bcd(k) = exp(-ky*Lx)
c mixup = array of bit reversed addresses for fft
c sct = sine/cosine table for fft
c t = complex scratch array, used during initialiation of fft tables
c we = bounded corrections to periodic electric field energy
c affp = normalization constant for poisson's equation
c indx = exponent which determines length in x direction, where nx=2**indx
c ny = system length in y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nxhd = must be >= nx/2
c nyhd = must be >= ny/2
      implicit none
      complex fx, fy, ffb, bv, sct, t
      integer isign, mixup, indx, ny, kstrt, nyv, kxp, jblok, nxhd, nyhd
      real bcd, we, affp
      dimension fx(nyv,kxp,jblok), fy(nyv,kxp,jblok)
      dimension ffb(nyhd,kxp,jblok), bv(nyhd,6), bcd(nyhd)
      dimension mixup(nxhd), sct(nxhd), t(nxhd)
c local data
      double precision wp, wb
      complex zc, zd, zt1, zt2, zt3, zt4
      integer nx, nxh, nx1, nyh, ny2
      integer is, j, j1, j2, j3, k, k1, ks, l, joff
      real dny, anx, anxi, dky, at1, at2, at3, at4, rho, rholx, dkyi
      real sum1, sum2, sum3, sum4
      nx = 2**indx
      nxh = nx/2
      nx1 = nx + 1
      nyh = ny/2
      ny2 = ny + 2
      ks = kstrt - 2
      dny = 6.28318530717959/float(ny)
      anx = float(nx)
c initialization
      if (isign.ne.0) go to 50
      if (kstrt.gt.nxh) return
c prepare fft tables
      is = 0
c     call FFT1RX(ffb,t,is,mixup,sct,indx,nx,nxh)
      is = -1
c prepare form factor array
      do 10 j = 1, nxh
      j1 = j - 1
      j2 = j1 + j1
      j3 = j2 + 1
c     ffb(j,1) = cmplx(float(j2*(nx - j2)),float(j3*(nx - j3)))
c     ffb(j,2) = cmplx(float(nxh - j2),float(nxh - j3))
   10 continue
c     call FFT1RX(ffb(1,1),t,is,mixup,sct,indx,nx,nxh)
c     call FFT1RX(ffb(1,2),t,is,mixup,sct,indx,nx,nxh)
      do 20 j = 1, nxh
c     ffb(j,1) = cmplx(real(ffb(j,1)),aimag(ffb(j,2)))
   20 continue
      do 40 k = 2, nyh
      dky = dny*float(k - 1)
      do 30 j = 1, nxh
      j2 = j + j
      j1 = j2 - 1
c     ffb(j,k) = cmplx(exp(-amin1(50.,dky*float(nx1 - j1))),exp(-amin1(5
c    10.,dky*float(nx1 - j2))))
   30 continue
      bcd(k) = exp(-amin1(50.,dky*anx))
c     call FFT1RX(ffb(1,k),t,is,mixup,sct,indx,nx,nxh)
   40 continue
      return
   50 if (isign.gt.0) go to 100
c calculate force/charge and sum field energy
      anxi = 1./anx
      wp = 0.0d0
      wb = 0.0d0
      do 90 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
c find constants for solution of homogeneous equation
      zc = .5*(bv(k,1) - bv(k,3) - bv(k,4))
      zd = .5*(bv(k,2) - bv(k,3) + bv(k,4))
      zt1 = zc - zd
      zt2 = cmplx(0.,1.)*conjg(zc + zd)
c boundary fields
      zt3 = zc + zd*bcd(k)
      zt4 = zc*bcd(k) + zd
c calculate internal and boundary energy corrections
      at2 = anxi/dky
      wp = wp + (aimag(zt2*bv(k,3)) - real(conjg(zt1)*bv(k,4)))*at2*(1. 
     1- bcd(k))
      wb = wb + (conjg(bv(k,1))*(bv(k,3) + zt3) + conjg(bv(k,2))*(bv(k,3
     1) + zt4))*at2
c homogenous electric field in x direction at x = Lx
      bv(k,5) = zc*bcd(k) - zd
c homogenous electric field in y direction at x = Lx
      bv(k,6) = -cmplx(0.,1.)*zt4
c homogenous electric field in x direction at x = 0
c     bv(k,7) = zc - zd*bcd(k)
c homogenous electric field in y direction at x = 0
c     bv(k,8) = -cmplx(0.,1.)*zt3
c calculate extra term in homogeneous solution
      zc = zc*(1. - bcd(k))*anxi
      zd = -cmplx(0.,1.)*zc
c find boundary values of periodic solution
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
c add solutions of homogeneous equation to periodic solution
      do 60 j = 2, nxh
      sum1 = sum1 + real(fx(k,j,l) + fx(k1,j,l))
      sum2 = sum2 + aimag(fx(k,j,l) - fx(k1,j,l))
      sum3 = sum3 + real(fy(k,j,l) + fy(k1,j,l))
      sum4 = sum4 + aimag(fy(k,j,l) - fy(k1,j,l))
      fx(k,j,l) = fx(k,j,l) + zt1*real(ffb(k,j,l)) + conjg(zt2)*aimag(ff
     1b(k,j,l)) + zc
      fx(k1,j,l) = fx(k1,j,l) + conjg(zt1)*real(ffb(k,j,l)) - zt2*aimag(
     1ffb(k,j,l)) + conjg(zc)
      fy(k,j,l) = fy(k,j,l) + conjg(zt2)*real(ffb(k,j,l)) - zt1*aimag(ff
     1b(k,j,l)) + zd
      fy(k1,j,l) = fy(k1,j,l) + zt2*real(ffb(k,j,l)) + conjg(zt1)*aimag(
     1ffb(k,j,l)) + conjg(zd)
   60 continue
c modes with n = 0, nx/2 are special
      sum1 = sum1 + real(fx(k,1,l))
      sum2 = sum2 + aimag(fx(k,1,l))
      sum3 = sum3 + real(fy(k,1,l))
      sum4 = sum4 + aimag(fy(k,1,l))
      fx(k,1,l) = zt1*real(ffb(k,1,l)) + zc
      fx(k1,1,l) = conjg(zt1)*aimag(ffb(k,1,l)) + conjg(zc)
      fy(k,1,l) = fy(k,1,l) + conjg(zt2)*real(ffb(k,1,l)) + zd
      fy(k1,1,l) = zt2*aimag(ffb(k,1,l)) + conjg(zd)
c electric field in x direction at x = Lx
      bv(k,5) = bv(k,5) + cmplx(sum1,sum2)
c electric field in y direction at x = Lx
      bv(k,6) = bv(k,6) + cmplx(sum3,sum4)
c electric field in x direction at x = 0
c     bv(k,7) = bv(k,7) + cmplx(sum1,sum2)
c electric field in y direction at x = 0
c     bv(k,8) = bv(k,8) + cmplx(sum3,sum4)
   70 continue
c find constants for solution of homogeneous equation
      rho = aimag(bv(1,4))
      rholx = .5*rho*anx
c find constants for solution of homogeneous equation
      at1 = rho*aimag(ffb(1,1,l))
      at2 = -(.5*(bv(1,2) - bv(1,1)) + bv(1,4))
      at3 = .5*at2*anx
c calculate energies
      wp = wp - .5*(at2*bv(1,4) - rho*(rholx*anx/6. - 2.*real(bv(1,3))))
      wb = wb - .5*(bv(1,2) - bv(1,1))*at2
      we = anx*float(ny)*(wp + wb)/affp
c find boundary values of periodic solution
      sum1 = 0.
      sum2 = 0.
c add solution of homogeneous equation to periodic solution
      do 80 j = 2, nxh
      sum1 = sum1 + real(fx(1,j,l))
      sum2 = sum2 + real(fy(1,j,l))
      if (rho.ne.0.) fx(1,j,l) = fx(1,j,l) - cmplx(at1,rho*aimag(ffb(1,j
     1,l)))
   80 continue
      sum1 = sum1 + sum1
      sum2 = sum2 + sum2
      fx(1,1,l) = cmplx(at2 - at1,-at1)
c electric field in x direction at x = Lx
      bv(1,5) = cmplx(at2+rholx,0.) + cmplx(sum1,0.)
c electric field in y direction at x = Lx
      bv(1,6) = cmplx(sum2,0.)
c electric field in x direction at x = 0
c     bv(1,7) = cmplx(at2-rholx,0.) + cmplx(sum1,0.)
c electric field in y direction at x = 0
c     bv(1,8) = cmplx(sum2,0.)
   90 continue
      return
c calculate potential and sum field energy
  100 anxi = 1./anx
      wp = 0.0d0
      wb = 0.0d0
      do 140 l = 1, jblok
      joff = kxp*(l + ks) - 1
      do 120 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dkyi = .5/dky
c find constants for solution of homogeneous equation
      zc = dkyi*(bv(k,1) - bv(k,3) - bv(k,4))
      zd = dkyi*(bv(k,2) - bv(k,3) + bv(k,4))
      zt1 = zc + zd
      zt2 = cmplx(0.,1.)*conjg(zc - zd)
c boundary potentials
      zt3 = zc + zd*bcd(k)
      zt4 = zc*bcd(k) + zd
c calculate internal and boundary energy corrections
      at2 = anxi/dky
      wp = wp + (real(conjg(zt1)*bv(k,3)) - aimag(zt2*bv(k,4)))*anxi*(1.
     1- bcd(k))
      wb = wb + (conjg(bv(k,1))*(bv(k,3) + zt3*dky) + conjg(bv(k,2))*(bv
     1(k,3) + zt4*dky))*at2
c calculate extra term in homogeneous solution
      zc = zc*(1. - bcd(k))*anxi
c find boundary values of periodic solution
      sum1 = 0.
      sum2 = 0.
c add solutions of homogeneous equation to periodic solution
      do 110 j = 2, nxh
      sum1 = sum1 + real(fx(k,j,l) + fx(k1,j,l))
      sum2 = sum2 + aimag(fx(k,j,l) - fx(k1,j,l))
      fx(k,j,l) = fx(k,j,l) + zt1*real(ffb(k,j,l)) + conjg(zt2)*aimag(ff
     1b(k,j,l)) + zc
      fx(k1,j,l) = fx(k1,j,l) + conjg(zt1)*real(ffb(k,j,l)) - zt2*aimag(
     1ffb(k,j,l)) + conjg(zc)
  110 continue
c modes with n = 0, nx/2 are special
      sum1 = sum1 + real(fx(k,1,l))
      sum2 = sum2 + aimag(fx(k,1,l))
      fx(k,1,l) = fx(k,1,l) + zt1*real(ffb(k,1,l)) + zc
      fx(k1,1,l) = conjg(zt1)*aimag(ffb(k,1,l)) + conjg(zc)
c potential at x = Lx
      bv(k,6) = zt4 + cmplx(sum1,sum2)
c potential at x = 0
c     bv(k,5) = zt3 + cmplx(sum1,sum2)
  120 continue
c find constants for solution of homogeneous equation
      rho = aimag(bv(1,4))
      rholx = .5*rho*anx
c find constants for solution of homogeneous equation
      at2 = -(.5*(bv(1,2) - bv(1,1)) + bv(1,4))
      at3 = .5*at2*anx
      at1 = at2*aimag(ffb(1,1,l))
      at4 = -bv(1,3)
c calculate energies
      wp = wp - .5*(at2*bv(1,4) - rho*(rholx*anx/6. + 2.*at4))
      wb = wb - .5*(bv(1,2) - bv(1,1))*at2
      we = anx*float(ny)*(wp + wb)/affp
c homogenous potential at x = Lx
      bv(1,6) = cmplx(-at3,0.)
c homogenous potential at x = 0
c     bv(1,5) = cmplx(at3,0.)
      at3 = .5*rho
c find boundary values of periodic solution
      sum1 = 0.
c add solution of homogeneous equation to periodic solution
      do 130 j = 2, nxh
      sum1 = sum1 + real(fx(1,j,l))
      fx(1,j,l) = fx(1,j,l) + cmplx(at1 + at3*real(ffb(1,j,l)),at2*aimag
     1(ffb(1,j,l)))
  130 continue
      sum1 = sum1 + sum1
      fx(1,1,l) = at3*conjg(ffb(1,1,l)) + cmplx(at1 + at4,at1)
c potential at x = Lx
      bv(1,6) = bv(1,6) + cmplx(sum1+at4,0.)
c potential at x = 0
c     bv(1,5) = bv(1,5) + cmplx(sum1+at4,0.)
  140 continue
      return
      end
