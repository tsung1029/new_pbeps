c 2d parallel PIC library for solving field equations with open (vacuum)
c boundary conditions
c written by viktor k. decyk, ucla
c copyright 1991, regents of the university of california
c update: march 15, 2006
c-----------------------------------------------------------------------
      subroutine PFORMC2(ffg,f,ft,bs,br,fpotc,mixup2,sct2,affp,ar,indx1,
     1indy1,kstrt,nxv,ny2d,kxp2,kyp2,j2blok,k2blok,kxp2d,ny1d,nxhy2,nxyh
     22)
c this subroutine calculates the form factor array ffg needed by field
c solvers with open (vacuum) boundary conditions using hockney's method.
c the four green's functions calculated are:
c g(kx,ky) = affp*inverse FFT of potr
c s(kx,ky) = inverse FFT of the density of a finite-sized particle
c gx(kx,ky) = affp*s(kx,ky)*inverse FFT of (x/r)*Er
c gy(kx,ky) = affp*s(kx,ky)*inverse FFT of (y/r)*Er
c where the fields due to the finite-sized particles are given by fpotc
c input: fpotc,mixup2,sct2,affp,ar,indx1,indy1,kstrt,nxv,ny2d,kxp2,kyp2,
c        j2blok,k2blok,kxp2d,ny1d,nxhy2,nxyh2)
c output: ffg, f
c ffg(1,k,j,l) = potential green's function g
c ffg(2,k,j,l) = finite-size particle shape factor s
c ffg(3,k,j,l) = x component of electric field green's function gx
c ffg(4,k,j,l) = y component of electric field green's function gy
c on processor 0, ffg(i,k,kxp2+1,l) = ffg(i,k,NX+1,l)
c on other processors, ffg(i,k,kxp2+1,l) = ffg(i,k,1,l) on processor 0
c f, ft = scratch arrays used by FFT
c fpotc = a function which calculates green's function
c mixup2/sct2 = bit-reverse and sine-cosine table used by FFT
c affp = normalization constant = nx*ny/np, where np=number of particles
c ar = half-width of particle in r direction
c indx1/indy1 = exponent which determines FFT length in x/y direction,
c where 2*nx=2**indx1, 2*ny=2**indy1
c kstrt = starting data block number
c nxv = half of first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
c kxp2/kyp2 = number of data values per block in x/y
c j2blok/k2blok = number of data blocks in x/y
c kxp2d = third dimension of ffg arrays, must be >= nx+1
c ny1d = second dimension of field arrays, must be >= ny+1
c nxhy2 = maximum of (nx,2*ny)
c nxyh2 = maximum of (nx,ny)
      implicit none
      real ffg, f
      complex ft, bs, br
      integer mixup2
      complex sct2
      real affp, ar
      integer indx1, indy1, kstrt, kxp2d, ny1d, nxv, ny2d, kxp2, kyp2
      integer j2blok, k2blok, nxhy2, nxyh2
      dimension ffg(4,ny1d,kxp2d,j2blok)
      dimension f(2*nxv,kyp2,k2blok), ft(ny2d,kxp2,j2blok)
      dimension bs(kxp2,kyp2,k2blok), br(kxp2,kyp2,j2blok)
      dimension mixup2(nxhy2), sct2(nxyh2)
      real fpotc
      external fpotc
c local data
      integer ntpose, nx, ny, ny1, nx2, ny2, isign, j, k, l, j1, k1, ks
      integer joff, koff, ifun
      real an, ari, at1, x, y, r, ttp
      real POTC2
      external POTC2
      data ntpose /1/
      nx2 = 2**(indx1)
      ny2 = 2**(indy1)
      nx = nx2/2
      ny = ny2/2
      ny1 = ny + 1
      ks = kstrt - 2
      ari = 0.0
      if (ar.gt.0.) ari = 1.0/ar
      an = float(nx2*ny2)
c calculate potential green's function
      ifun = 1
      if (kstrt.gt.ny2) go to 40
      do 30 l = 1, k2blok
      koff = kyp2*(l + ks) - 1
      do 20 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at1 = float(k1)**2
      do 10 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at1 + float(j1)**2)
      f(j,k,l) = fpotc(r,affp,ari,1)
   10 continue
   20 continue
   30 continue
      isign = -1
      call WPFFT2R(f,ft,bs,br,isign,ntpose,mixup2,sct2,ttp,indx1,indy1,k
     1strt,nxv,ny2d,kxp2,kyp2,kyp2,j2blok,k2blok,nxhy2,nxyh2)
   40 if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
      do 60 j = 1, kxp2
      do 50 k = 1, ny1
      ffg(ifun,k,j,l) = an*real(ft(k,j,l))
   50 continue
   60 continue
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 + 2 - k
         ffg(ifun,k,kxp2+1,l) = an*real(ft(k1,1,l))
   70    continue
         ffg(ifun,1,kxp2+1,l) = an*aimag(ft(1,1,l))
         ffg(ifun,ny1,kxp2+1,l) = an*aimag(ft(ny1,1,l))
      else
         do 80 k = 1, ny1
         ffg(ifun,k,kxp2+1,l) = 0.0
   80    continue
      endif
   90 continue
c calculate particle smoothing function
  100 ifun = ifun + 1
      if (kstrt.gt.ny2) go to 140
      do 130 l = 1, k2blok
      koff = kyp2*(l + ks) - 1
      do 120 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at1 = float(k1)**2
      do 110 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at1 + float(j1)**2)
      f(j,k,l) = POTC2(r,affp,ari,2)
  110 continue
  120 continue
  130 continue
      isign = -1
      call WPFFT2R(f,ft,bs,br,isign,ntpose,mixup2,sct2,ttp,indx1,indy1,k
     1strt,nxv,ny2d,kxp2,kyp2,kyp2,j2blok,k2blok,nxhy2,nxyh2)
  140 if (kstrt.gt.nx) go to 200
      do 190 l = 1, j2blok
      do 160 j = 1, kxp2
      do 150 k = 1, ny1
      ffg(ifun,k,j,l) = an*real(ft(k,j,l))
  150 continue
  160 continue
      if ((l+ks).eq.0) then
         do 170 k = 2, ny
         k1 = ny2 + 2 - k
         ffg(ifun,k,kxp2+1,l) = an*real(ft(k1,1,l))
  170    continue
         ffg(ifun,1,kxp2+1,l) = an*aimag(ft(1,1,l))
         ffg(ifun,ny1,kxp2+1,l) = an*aimag(ft(ny1,1,l))
      else
         do 180 k = 1, ny1
         ffg(ifun,k,kxp2+1,l) = 0.0
  180    continue
      endif
  190 continue
c calculate green's function for x component of electric field
  200 ifun = ifun + 1
      if (kstrt.gt.ny2) go to 240
      do 230 l = 1, k2blok
      koff = kyp2*(l + ks) - 1
      do 220 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at1 = float(k1)**2
      do 210 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      x = float(j1)
      r = sqrt(at1 + x*x)
      f(j,k,l) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k,l) = f(j,k,l)*(x/r)
  210 continue
  220 continue
  230 continue
      isign = -1
      call WPFFT2R(f,ft,bs,br,isign,ntpose,mixup2,sct2,ttp,indx1,indy1,k
     1strt,nxv,ny2d,kxp2,kyp2,kyp2,j2blok,k2blok,nxhy2,nxyh2)
  240 if (kstrt.gt.nx) go to 300
      do 290 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 260 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 250 k = 1, ny1
         ffg(ifun,k,j,l) = an*aimag(ft(k,j,l))
  250    continue
      endif
  260 continue
      if ((l+ks).eq.0) then
         do 270 k = 2, ny
         k1 = ny2 + 2 - k
         ffg(ifun,k,1,l) = an*real(ft(k,1,l))
         ffg(ifun,k,kxp2+1,l) = an*real(ft(k1,1,l))
  270    continue
         ffg(ifun,1,1,l) = an*real(ft(1,1,l))
         ffg(ifun,1,kxp2+1,l) = an*aimag(ft(1,1,l))
         ffg(ifun,ny1,1,l) = an*real(ft(ny1,1,l))
         ffg(ifun,ny1,kxp2+1,l) = an*aimag(ft(ny1,1,l))
      else
         do 280 k = 1, ny1
         ffg(ifun,k,kxp2+1,l) = 0.0
  280    continue
      endif
  290 continue
c calculate green's function for y component of electric field
  300 ifun = ifun + 1
      if (kstrt.gt.ny2) go to 340
      do 330 l = 1, k2blok
      koff = kyp2*(l + ks) - 1
      do 320 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      y = float(k1)
      at1 = y*y
      do 310 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at1 + float(j1)**2)
      f(j,k,l) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k,l) = f(j,k,l)*(y/r)
  310 continue
  320 continue
  330 continue
      isign = -1
      call WPFFT2R(f,ft,bs,br,isign,ntpose,mixup2,sct2,ttp,indx1,indy1,k
     1strt,nxv,ny2d,kxp2,kyp2,kyp2,j2blok,k2blok,nxhy2,nxyh2)
  340 if (kstrt.gt.nx) go to 400
      do 390 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 360 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 350 k = 2, ny
         ffg(ifun,k,j,l) = an*aimag(ft(k,j,l))
  350    continue
         ffg(ifun,1,j,l) = an*real(ft(1,j,l))
         ffg(ifun,ny1,j,l) = an*real(ft(ny1,j,l))
      endif
  360 continue
      if ((l+ks).eq.0) then
         do 370 k = 2, ny
         k1 = ny2 + 2 - k
         ffg(ifun,k,1,l) = an*aimag(ft(k,1,l))
         ffg(ifun,k,kxp2+1,l) = an*aimag(ft(k1,1,l))
  370    continue
         ffg(ifun,1,1,l) = an*real(ft(1,1,l))
         ffg(ifun,1,kxp2+1,l) = an*aimag(ft(1,1,l))
         ffg(ifun,ny1,1,l) = an*real(ft(ny1,1,l))
         ffg(ifun,ny1,kxp2+1,l) = an*aimag(ft(ny1,1,l))
      else
         do 380 k = 1, ny1
         ffg(ifun,k,kxp2+1,l) = 0.0
  380    continue
      endif
  390 continue
c copy ffg(i,k,1,l) on node 0 to ffg(i,k,kxp2+1,l) on other nodes
  400 do 410 l = 1, j2blok
      call P0COPY(ffg(1,1,1,l),ffg(1,1,kxp2+1,l),4*ny1d)
  410 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISC2(q,fx,fy,isign,ffg,we,nx,ny,kstrt,ny2d,kxp2,j2bl
     1ok,ny1d,kxp2d)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function, with open (vacuum)
c boundary conditions using hockney's method, for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate convolution
c for isign = -1, 
c input: q,ffg,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,ny1d,kxp2d
c output: fx,fy,we
c approximate flop count is: 44*nx*ny + 36*(nx + ny)
c for isign = 1,
c input: q,ffg,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,ny1d,kxp2d
c output: fx,we
c approximate flop count is: 22*nx*ny + 24(nx + ny)
c for isign = 2,
c input: q,ffg,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,ny1d,kxp2d
c output: fy
c approximate flop count is: 4*nx*ny + 2*(nx + ny)
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = gx(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = gy(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c gx(kx,ky) = s(kx,ky)*inverse FFT of (x/r)*Er
c gy(kx,ky) = s(kx,ky)*inverse FFT of (y/r)*Er
c where Er is the electric field of a single finite-sized particle
c s(kx,ky) = inverse FFT of the density of a finite-sized particle
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)
c where g(kx,ky) = affp*inverse FFT of potr
c where potr is the potential of a single finite-sized particle
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fx(k,j,l) = x component of complex force/charge,
c fy(k,j,l) = y component of complex force/charge,
c ffg(1,k,j,l) = potential green's function g
c ffg(2,k,j,l) = finite-size particle shape factor s
c ffg(3,k,j,l) = x component of electric field green's function gx
c ffg(4,k,j,l) = y component of electric field green's function gy
c all for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c the ffg array is calculated by the subroutine PFORMC2
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c ny2d = first dimension of field arrays, must be >= 2*ny
c kxp2 = number of data values per block
c j2blok = number of data blocks
c electric field energy is also calculated and returned in we
c ny1d = second dimension of ffg array, must be >= ny+1
c kxp2d = third dimension of ffg array, must be >= kxp2+1
      implicit none
      real ffg
      complex q, fx, fy
      integer isign, nx, ny, kstrt, ny2d, kxp2, j2blok, ny1d, kxp2d
      real we
      dimension q(ny2d,kxp2,j2blok)
      dimension fx(ny2d,kxp2,j2blok), fy(ny2d,kxp2,j2blok)
      dimension ffg(4,ny1d,kxp2d,j2blok)
c local data
      double precision wp
      integer j, k, l, k1, ny22, ks, kx1, joff, ny1
      real at1, at2, at3, at4
      complex zt1, zt2
      if (isign.eq.0) return
      ny1 = ny + 1
      ny22 = ny + ny + 2
      ks = kstrt - 2
      kx1 = 1
      if (isign.gt.0) go to 70
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 60
      do 50 l = 1, j2blok
      if ((l+ks).gt.0) kx1 = kxp2 + 1
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      at3 = -1.0
      do 20 j = 1, kxp2
      at3 = -at3
      if ((j+joff).gt.0) then
         at2 = ffg(4,1,j,l)
         do 10 k = 2, ny
         k1 = ny22 - k
         at1 = at3*ffg(3,k,kx1,l)
         at2 = -at2
         zt1 = cmplx(at1,ffg(3,k,j,l))
         zt2 = cmplx(at2,ffg(4,k,j,l))
         fx(k,j,l) = zt1*q(k,j,l)
         fx(k1,j,l) = zt1*q(k1,j,l)
         fy(k,j,l) = zt2*q(k,j,l)
         fy(k1,j,l) = conjg(zt2)*q(k1,j,l)
         wp = wp + ffg(1,k,j,l)*(q(k,j,l)*conjg(q(k,j,l)) + q(k1,j,l)*co
     1njg(q(k1,j,l)))
   10    continue
c mode number ky = 0
         at1 = at3*ffg(3,1,kx1,l)
         zt1 = cmplx(at1,ffg(3,1,j,l))
         fx(1,j,l) = zt1*q(1,j,l)
         fy(1,j,l) = ffg(4,1,j,l)*q(1,j,l)
         wp = wp + ffg(1,1,j,l)*(q(1,j,l)*conjg(q(1,j,l)))
c mode number ky = ny
         at1 = at3*ffg(3,ny1,kx1,l)
         zt1 = cmplx(at1,ffg(3,ny1,j,l))
         fx(ny1,j,l) = zt1*q(ny1,j,l)
         fy(ny1,j,l) = ffg(4,ny1,j,l)*q(ny1,j,l)
         wp = wp + ffg(1,ny1,j,l)*(q(ny1,j,l)*conjg(q(ny1,j,l)))
      endif
   20 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
c mode number kx = 0
         at3 = ffg(4,1,1,l)
         do 30 k = 2, ny
         at3 = -at3
         zt1 = cmplx(at3,ffg(4,k,1,l))
         fx(k,1,l) = ffg(3,k,1,l)*q(k,1,l)
         fy(k,1,l) = zt1*q(k,1,l)
         wp = wp + ffg(1,k,1,l)*(q(k,1,l)*conjg(q(k,1,l)))
   30    continue
c mode number kx = nx/2
         at3 = ffg(4,1,kxp2+1,l)
         do 40 k = 2, ny
         k1 = ny22 - k
         at3 = -at3
         zt1 = cmplx(at3,ffg(4,k,kxp2+1,l))
         fx(k1,1,l) = ffg(3,k,kxp2+1,l)*q(k1,1,l)
         fy(k1,1,l) = zt1*q(k1,1,l)
         wp = wp + ffg(1,k,kxp2+1,l)*(q(k1,1,l)*conjg(q(k1,1,l)))
   40    continue
c mode numbers ky = 0, kx = 0, nx/2
         fx(1,1,l) = cmplx(ffg(3,1,1,l)*real(q(1,1,l)),ffg(3,1,kxp2+1,l)
     1*aimag(q(1,1,l)))
         fy(1,1,l) = cmplx(ffg(4,1,1,l)*real(q(1,1,l)),ffg(4,1,kxp2+1,l)
     1*aimag(q(1,1,l)))
         wp = wp + .5*(ffg(1,1,1,l)*real(q(1,1,l))**2 + ffg(1,1,kxp2+1,l
     1)*aimag(q(1,1,l))**2)
c mode numbers ky = ny/2, kx = 0, nx/2
         fx(ny1,1,l) = cmplx(ffg(3,ny1,1,l)*real(q(ny1,1,l)),ffg(3,ny1,k
     1xp2+1,l)*aimag(q(ny1,1,l)))
         fy(ny1,1,l) = cmplx(ffg(4,ny1,1,l)*real(q(ny1,1,l)),ffg(4,ny1,k
     1xp2+1,l)*aimag(q(ny1,1,l)))
         wp = wp + .5*(ffg(1,ny1,1,l)*real(q(ny1,1,l))**2 + ffg(1,ny1,kx
     1p2+1,l)*aimag(q(ny1,1,l))**2)
      endif
   50 continue
   60 continue
      we = 4.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
   70 if (isign.gt.1) go to 140
      wp = 0.0d0
      if (kstrt.gt.nx) go to 130
      do 120 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 90 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 80 k = 2, ny
         k1 = ny22 - k
         at2 = ffg(1,k,j,l)
         at1 = at2*ffg(2,k,j,l)
c        at1 = at2
         fx(k,j,l) = at2*q(k,j,l)
         fx(k1,j,l) = at2*q(k1,j,l)
         wp = wp + at1*(q(k,j,l)*conjg(q(k,j,l)) + q(k1,j,l)*conjg(q(k1,
     1j,l)))
   80    continue
c mode number ky = 0
         at2 = ffg(1,1,j,l)
         at1 = at2*ffg(2,1,j,l)
c        at1 = at2
         fx(1,j,l) = at2*q(1,j,l)
         wp = wp + at1*(q(1,j,l)*conjg(q(1,j,l)))
c mode number ky = ny
         at2 = ffg(1,ny1,j,l)
         at1 = at2*ffg(2,ny1,j,l)
c        at1 = at2
         fx(ny1,j,l) = at2*q(ny1,j,l)
         wp = wp + at1*(q(ny1,j,l)*conjg(q(ny1,j,l)))
      endif
   90 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
c mode number kx = 0
         do 100 k = 2, ny
         at2 = ffg(1,k,1,l)
         at1 = at2*ffg(2,k,1,l)
c        at1 = at2
         fx(k,1,l) = at2*q(k,1,l)
         wp = wp + at1*(q(k,1,l)*conjg(q(k,1,l)))
  100    continue
c mode number kx = nx/2
         do 110 k = 2, ny
         k1 = ny22 - k
         at2 = ffg(1,k,kxp2+1,l)
         at1 = at2*ffg(2,k,kxp2+1,l)
c        at1 = at2
         fx(k1,1,l) = at2*q(k1,1,l)
         wp = wp + at1*(q(k1,1,l)*conjg(q(k1,1,l)))
  110    continue
c mode numbers ky = 0, kx = 0, nx/2
         at2 = ffg(1,1,1,l)
         at1 = at2*ffg(2,1,1,l)
c        at1 = at2
         at4 = ffg(1,1,kxp2+1,l)
         at3 = at4*ffg(2,1,kxp2+1,l)
c        at3 = at4
         fx(1,1,l) = cmplx(at2*real(q(1,1,l)),at4*aimag(q(1,1,l)))
         wp = wp + .5*(at1*real(q(1,1,l))**2 + at3*aimag(q(1,1,l))**2)
c mode numbers ky = ny/2, kx = 0, nx/2
         at2 = ffg(1,ny1,1,l)
         at1 = at2*ffg(2,ny1,1,l)
c        at1 = at2
         at4 = ffg(1,ny1,kxp2+1,l)
         at3 = at4*ffg(2,ny1,kxp2+1,l)
c        at3 = at4
         fx(ny1,1,l) = cmplx(at2*real(q(ny1,1,l)),at4*aimag(q(ny1,1,l)))
         wp = wp + .5*(at1*real(q(ny1,1,l))**2 + at3*aimag(q(ny1,1,l))**
     12)
      endif
  120 continue
  130 continue
      we = 4.0*float(nx*ny)*wp
      return
c calculate smoothing
  140 if (kstrt.gt.nx) go to 200
      do 190 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 160 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 150 k = 2, ny
         k1 = ny22 - k
         at1 = ffg(2,k,j,l)
         fy(k,j,l) = at1*q(k,j,l)
         fy(k1,j,l) = at1*q(k1,j,l)
  150    continue
c mode number ky = 0
         fy(1,j,l) = ffg(2,1,j,l)*q(1,j,l)
c mode number ky = ny
         fy(ny1,j,l) = ffg(2,ny1,j,l)*q(ny1,j,l)
      endif
  160 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
c mode number kx = 0
         do 170 k = 2, ny
         fy(k,1,l) = ffg(2,k,1,l)*q(k,1,l)
  170    continue
c mode number kx = nx/2
         do 180 k = 2, ny
         k1 = ny22 - k
         fy(k1,1,l) = ffg(2,k,kxp2+1,l)*q(k1,1,l)
  180    continue
c mode numbers ky = 0, kx = 0, nx/2
         at1 = ffg(2,1,1,l)
         at3 = ffg(2,1,kxp2+1,l)
         fy(1,1,l) = cmplx(at1*real(q(1,1,l)),at3*aimag(q(1,1,l)))
c mode numbers ky = ny/2, kx = 0, nx/2
         at1 = ffg(2,ny1,1,l)
         at3 = ffg(2,ny1,kxp2+1,l)
         fy(ny1,1,l) = cmplx(at1*real(q(ny1,1,l)),at3*aimag(q(ny1,1,l)))
      endif
  190 continue
  200 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISC22(q,fxy,ffg,we,nx,ny,kstrt,ny2d,kxp2,j2blok,ny1d
     1,kxp2d)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with open (vacuum) boundary conditions using hockney's method,
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate convolution
c input: q,ffg,nx,ny,kstrt,ny2d,kxp2,j2blok,ny1d,kxp2d
c output: fxy,we
c approximate flop count is: 44*nx*ny + 36*(nx + ny)
c force/charge is calculated using the equations:
c fx(kx,ky) = gx(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = gy(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c gx(kx,ky) = s(kx,ky)*inverse FFT of (x/r)*Er
c gy(kx,ky) = s(kx,ky)*inverse FFT of (y/r)*Er
c where Er is the electric field of a single finite-sized particle
c s(kx,ky) = inverse FFT of the density of a finite-sized particle
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c ffg(1,k,j,l) = potential green's function g
c ffg(2,k,j,l) = finite-size particle shape factor s
c ffg(3,k,j,l) = x component of electric field green's function gx
c ffg(4,k,j,l) = y component of electric field green's function gy
c all for fourier mode (jj-1,k-1), where jj = j + kxp2*(l - 1)
c the ffg array is calculated by the subroutine PFORMC2
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c ny2d = first dimension of field arrays, must be >= 2*ny
c kxp2 = number of data values per block
c j2blok = number of data blocks
c electric field energy is also calculated and returned in we
c ny1d = second dimension of ffg array, must be >= ny+1
c kxp2d = third dimension of ffg array, must be >= kxp2+1
      implicit none
      real ffg
      complex q, fxy
      integer nx, ny, kstrt, ny2d, kxp2, j2blok, ny1d, kxp2d
      real we
      dimension q(ny2d,kxp2,j2blok), fxy(2,ny2d,kxp2,j2blok)
      dimension ffg(4,ny1d,kxp2d,j2blok)
c local data
      double precision wp
      integer j, k, l, k1, ny22, ks, kx1, joff, ny1
      real at1, at2, at3
      complex zt1, zt2
      ny1 = ny + 1
      ny22 = ny + ny + 2
      ks = kstrt - 2
      kx1 = 1
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 60
      do 50 l = 1, j2blok
      if ((l+ks).gt.0) kx1 = kxp2 + 1
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      at3 = -1.0
      do 20 j = 1, kxp2
      at3 = -at3
      if ((j+joff).gt.0) then
         at2 = ffg(4,1,j,l)
         do 10 k = 2, ny
         k1 = ny22 - k
         at1 = at3*ffg(3,k,kx1,l)
         at2 = -at2
         zt1 = cmplx(at1,ffg(3,k,j,l))
         zt2 = cmplx(at2,ffg(4,k,j,l))
         fxy(1,k,j,l) = zt1*q(k,j,l)
         fxy(1,k1,j,l) = zt1*q(k1,j,l)
         fxy(2,k,j,l) = zt2*q(k,j,l)
         fxy(2,k1,j,l) = conjg(zt2)*q(k1,j,l)
         wp = wp + ffg(1,k,j,l)*(q(k,j,l)*conjg(q(k,j,l)) + q(k1,j,l)*co
     1njg(q(k1,j,l)))
   10    continue
c mode number ky = 0
         at1 = at3*ffg(3,1,kx1,l)
         zt1 = cmplx(at1,ffg(3,1,j,l))
         fxy(1,1,j,l) = zt1*q(1,j,l)
         fxy(2,1,j,l) = ffg(4,1,j,l)*q(1,j,l)
         wp = wp + ffg(1,1,j,l)*(q(1,j,l)*conjg(q(1,j,l)))
c mode number ky = ny
         at1 = at3*ffg(3,ny1,kx1,l)
         zt1 = cmplx(at1,ffg(3,ny1,j,l))
         fxy(1,ny1,j,l) = zt1*q(ny1,j,l)
         fxy(2,ny1,j,l) = ffg(4,ny1,j,l)*q(ny1,j,l)
         wp = wp + ffg(1,ny1,j,l)*(q(ny1,j,l)*conjg(q(ny1,j,l)))
      endif
   20 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
c mode number kx = 0
         at3 = ffg(4,1,1,l)
         do 30 k = 2, ny
         at3 = -at3
         zt1 = cmplx(at3,ffg(4,k,1,l))
         fxy(1,k,1,l) = ffg(3,k,1,l)*q(k,1,l)
         fxy(2,k,1,l) = zt1*q(k,1,l)
         wp = wp + ffg(1,k,1,l)*(q(k,1,l)*conjg(q(k,1,l)))
   30    continue
c mode number kx = nx/2
         at3 = ffg(4,1,kxp2+1,l)
         do 40 k = 2, ny
         k1 = ny22 - k
         at3 = -at3
         zt1 = cmplx(at3,ffg(4,k,kxp2+1,l))
         fxy(1,k1,1,l) = ffg(3,k,kxp2+1,l)*q(k1,1,l)
         fxy(2,k1,1,l) = zt1*q(k1,1,l)
         wp = wp + ffg(1,k,kxp2+1,l)*(q(k1,1,l)*conjg(q(k1,1,l)))
   40    continue
c mode numbers ky = 0, kx = 0, nx/2-
         fxy(1,1,1,l) = cmplx(ffg(3,1,1,l)*real(q(1,1,l)),ffg(3,1,kxp2+1
     1,l)*aimag(q(1,1,l)))
         fxy(2,1,1,l) = cmplx(ffg(4,1,1,l)*real(q(1,1,l)),ffg(4,1,kxp2+1
     1,l)*aimag(q(1,1,l)))
         wp = wp + .5*(ffg(1,1,1,l)*real(q(1,1,l))**2 + ffg(1,1,kxp2+1,l
     1)*aimag(q(1,1,l))**2)
c mode numbers ky = ny/2, kx = 0, nx/2
         fxy(1,ny1,1,l) = cmplx(ffg(3,ny1,1,l)*real(q(ny1,1,l)),ffg(3,ny
     11,kxp2+1,l)*aimag(q(ny1,1,l)))
         fxy(2,ny1,1,l) = cmplx(ffg(4,ny1,1,l)*real(q(ny1,1,l)),ffg(4,ny
     11,kxp2+1,l)*aimag(q(ny1,1,l)))
         wp = wp + .5*(ffg(1,ny1,1,l)*real(q(ny1,1,l))**2 + ffg(1,ny1,kx
     1p2+1,l)*aimag(q(ny1,1,l))**2)
      endif
   50 continue
   60 continue
      we = 4.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      function POTC3(r,affp,ari,ifun)
c this function calculates the fields for finite-size gaussian particles
c in 3D:
c if ifun = 1, calculate potential function
c POTC3 = (affp/(4*pi))*erfn(r/(ar*sqrt(2.)))/r, for r > 0.
c POTC3 = (affp/(4*pi))*sqrt(2./3.14159265358979)/ar, for r = 0.
c if ifun = 2, calculate particle shape function
c POTC3 = exp(-(r/(sqrt(2.)*ar))**2)/(sqrt(2.*pi)*ar)**3, for r > 0.
c POTC3 = 1./(sqrt(2.*pi)*ar)**3, for r = 0.
c if ifun = 3, calculate radial electric field
c POTC3 = (affp/(4*pi))*(1/r)*(erf(r/(sqrt(2.)*ar))/r -
c exp(-(r/(sqrt(2.)*ar))**2)*sqrt(2./3.14159265358979)/ar, for r > 0.
c POTC3 = 0.0, for r = 0.
c where erfn is the error function
c and where the finite-size particle density is given by:
c rho(r) = exp(-(r/sqrt(2)*ar)**2)/(sqrt(2*pi)*ar)**3
c affp = 4*pi*e**2/(me*(omega0**2)*delta**3) = 1/(n0*delta**3)
c where n0*delta**3 = number density per grid
c r = radial coordinate
c affp = normalization constant
c ari = 1/ar = inverse of particle size function
c (ari = 0., means use point particle result)
c ifun = (1,2,3) = calculate (potential,shape,electric field)
      implicit none
      real r, affp, ari
      integer ifun
c local data
c pi4i = 1/4*pi, sqt2i = 1./sqrt(2.), sqt2pi = sqrt(2./pi)
      real pi4i, sqt2i, sqt2pi
      parameter(pi4i=0.5/6.28318530717959)
      parameter(sqt2i=0.707106781186548,sqt2pi=0.797884560802865)
      real POTC3, erfn
      external erfn
      real anorm, at1, ri
      anorm = affp*pi4i
c calculate potential function
      if (ifun.eq.1) then
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC3 = anorm*sqt2pi*ari
            else
               POTC3 = anorm*erfn(r*sqt2i*ari)/r
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC3 = 0.0
            else
               POTC3 = anorm/r
            endif
         endif
c calculate particle shape function
      else if (ifun.eq.2) then
         anorm = affp*(.5*sqt2pi*ari)**3
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC3 = anorm
            else
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC3 = anorm*exp(-(at1*at1))
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC3 = affp
            else
               POTC3 = 0.0
            endif
         endif
c calculate radial electric field
      else if (ifun.eq.3) then
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC3 = 0.0
            else
               ri = 1.0/r
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC3 = anorm*ri*(erfn(at1)*ri - sqt2pi*ari*exp(-(at1*at1
     1)))
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC3 = 0.0
            else
               POTC3 = anorm/(r*r)
            endif
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      function POTC2(r,affp,ari,ifun)
c this function calculates the fields for finite-size gaussian particles
c in 2D:
c if ifun = 1, calculate potential function
c POTC2 = -(affp/(4*pi))*(e1(r**2/(2*ar**2)) + ln(r**2)), for r > 0.
c POTC2 = -(affp/(4*pi))*(ln(2) - gamma + 2*ln(ar), for r = 0.
c if ifun = 2, calculate particle shape function
c POTC2 = exp(-(r/(sqrt(2.)*ar))**2)/(sqrt(2.*pi)*ar)**2, for r > 0.
c POTC2 = 1./(sqrt(2.*pi)*ar)**2, for r = 0.
c if ifun = 3, calculate radial electric field
c POTC2 = 2*(1 - exp(-(r/(sqrt(2.)*ar))**2)/r, for r > 0.
c POTC2 = 0.0, for r = 0.
c where e1 is the exponential integral
c and where the finite-size particle density is given by:
c rho(r) = exp(-(r/sqrt(2)*ar)**2)/(2*pi*ar**2), qm = q/e
c affp = 4*pi*e**2/(me*(omega0**2)*delta**2) = 1/(n0*delta**2)
c where n0*delta**2 = number density per grid
c r = radial coordinate
c affp = normalization constant
c ari = 1/ar = inverse of particle size function
c (ari = 0., means use point particle result)
c ifun = 1 = calculate (potential)
      implicit none
      real r, affp, ari
      integer ifun
c local data
c pi4i = 1/4*pi, sqt2i = 1./sqrt(2.), sqt2pi = sqrt(2./pi)
      real pi4i, sqt2i, sqt2pi
      parameter(pi4i=0.5/6.28318530717959)
      parameter(sqt2i=0.707106781186548,sqt2pi=0.797884560802865)
      real POTC2, e1ln
      external e1ln
      real anorm, at1
c calculate potential function
      if (ifun.eq.1) then
         anorm = -affp*pi4i
c finite-size particles
         if (ari.gt.0.) then
            POTC2 = anorm*(e1ln((r*sqt2i*ari)**2) - 2.0*alog(sqt2i*ari))
c point particles
         else
            if (r.eq.0.) then
               POTC2 = 0.0
            else
               POTC2 = 2.0*anorm*alog(r)
            endif
         endif
c calculate particle shape function
      else if (ifun.eq.2) then
         anorm = affp*(.5*sqt2pi*ari)**2
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC2 = anorm
            else
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC2 = anorm*exp(-(at1*at1))
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC2 = affp
            else
               POTC2 = 0.0
            endif
         endif
c calculate radial electric field
      else if (ifun.eq.3) then
         anorm = 2.*affp*pi4i
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC2 = 0.0
            else
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC2 = anorm*(1.0 - exp(-(at1*at1)))/r
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC2 = 0.0
            else
               POTC2 = anorm/r
            endif
         endif
      endif
      return
      end

