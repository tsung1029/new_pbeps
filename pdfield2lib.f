c 2d parallel PIC library for solving field equations with dirichlet
c boundary conditions
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: july 21, 2006
c-----------------------------------------------------------------------
      subroutine PLCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
c this subroutine replicates field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c for distributed data
      implicit none
      real fxy
      integer nyp, nx, nxe, nypmx, nblok
      dimension fxy(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp3
      do 20 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 10 k = 1, nyp3
      fxy(1,1,k,l) = 2.*fxy(1,2,k,l) - fxy(1,3,k,l)
      fxy(2,1,k,l) = 2.*fxy(2,2,k,l) - fxy(2,3,k,l)
      fxy(1,nx+3,k,l) = 2.*fxy(1,nx+2,k,l) - fxy(1,nx+1,k,l)
      fxy(2,nx+3,k,l) = 2.*fxy(2,nx+2,k,l) - fxy(2,nx+1,k,l)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
c this subroutine replicates scalar field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c for distributed data
      implicit none
      real q
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp3
      do 20 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 10 k = 1, nyp3
      q(1,k,l) = 2.*q(2,k,l) - q(3,k,l)
      q(nx+3,k,l) = 2.*q(nx+2,k,l) - q(nx+1,k,l)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLBGUARD2X(bxy,nyp,nx,nxe,nypmx,nblok)
c this subroutine replicates field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c for distributed data
      implicit none
      real bxy
      integer nyp, nx, nxe, nypmx, nblok
      dimension bxy(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp3
      do 20 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 10 k = 1, nyp3
      bxy(1,1,k,l) = 2.*bxy(1,2,k,l) - bxy(1,3,k,l)
      bxy(2,1,k,l) = 2.*bxy(2,2,k,l) - bxy(2,3,k,l)
      bxy(3,1,k,l) = 2.*bxy(3,2,k,l) - bxy(3,3,k,l)
      bxy(1,nx+3,k,l) = 2.*bxy(1,nx+2,k,l) - bxy(1,nx+1,k,l)
      bxy(2,nx+3,k,l) = 2.*bxy(2,nx+2,k,l) - bxy(2,nx+1,k,l)
      bxy(3,nx+3,k,l) = 2.*bxy(3,nx+2,k,l) - bxy(3,nx+1,k,l)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD2X(cu,kstrt,nvp,noff,nyp,xj0,yj0,zj0,nx,ny,ngx
     1,ngy,nxe,nypmx,nblok)
c initialize extended non-periodic field
c cu(i,j+1,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c kstrt = starting data block number
c nvp = number of real or virtual processors
c noff(l) = lowermost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c xj0/yj0/zj0 = initialization constants in x/y/z direction
c nx/ny = system length in x/y direction
c ngx/ngy = (0,1) number of grid cells away from edge
c nxe = first dimension of charge array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells
c nblok = number of particle partitions
c quadratic interpolation, for distributed data
      implicit none
      real cu, xj0, yj0, zj0
      integer kstrt, nvp, noff, nyp, nx, ny, ngx, ngy, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok), noff(nblok)
      integer i, j, k, l, ks, kk, nyp3, nxg, nx3
      real chx, chy, chz
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      ks = kstrt - 2
      chx = .5*xj0
      chy = .5*yj0
      chz = .5*zj0
      do 120 l = 1, nblok
      nyp3 = nyp(l) + 3
c handle first grid point in y
      if ((l+ks).eq.0) then
         do 20 j = 1, nx3
         do 10 i = 1, 3
         cu(i,j,2,l) = 0.
   10    continue
   20    continue
      endif
c handle interior grid points in y
      do 70 k = 1, nyp(l)
      kk = k + noff(l)
      if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
         do 30 j = 2, nxg
         cu(1,j+ngx+1,k+1,l) = xj0
         cu(2,j+ngx+1,k+1,l) = yj0
         cu(3,j+ngx+1,k+1,l) = zj0
   30    continue
         do 40 i = 1, 3
         cu(i,1,k+1,l) = 0.
         cu(i,2,k+1,l) = 0.
         cu(i,nx+2,k+1,l) = 0.
         cu(i,nx+3,k+1,l) = 0.
   40    continue
         cu(1,ngx+2,k+1,l) = chx
         cu(2,ngx+2,k+1,l) = chy
         cu(3,ngx+2,k+1,l) = chz
         cu(1,nx-ngx+2,k+1,l) = chx
         cu(2,nx-ngx+2,k+1,l) = chy
         cu(3,nx-ngx+2,k+1,l) = chz
      else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
         do 50 j = 2, nxg
         cu(1,j+ngx+1,k+1,l) = chx
         cu(2,j+ngx+1,k+1,l) = chy
         cu(3,j+ngx+1,k+1,l) = chz
   50    continue
         do 60 i = 1, 3
         cu(i,1,k+1,l) = 0.
         cu(i,2,k+1,l) = 0.
         cu(i,nx+2,k+1,l) = 0.
         cu(i,nx+3,k+1,l) = 0.
   60    continue
         cu(1,ngx+2,k+1,l) = .5*chx
         cu(2,ngx+2,k+1,l) = .5*chy
         cu(3,ngx+2,k+1,l) = .5*chz
         cu(1,nx-ngx+2,k+1,l) = .5*chx
         cu(2,nx-ngx+2,k+1,l) = .5*chy
         cu(3,nx-ngx+2,k+1,l) = .5*chz
      endif
   70 continue
c guard cells in y
      do 90 j = 1, nx3
      do 80 i = 1, 3
      cu(i,j,1,l) = 0.
      cu(i,j,nyp3-1,l) = 0.
      cu(i,j,nyp3,l) = 0.
   80 continue
   90 continue
c handle last grid point in y
      if (((nyp(l)+noff(l)).lt.(ny-ngy+1)).and.((l+ks).eq.(nvp-1))) then
         do 100 j = 2, nxg
         cu(1,j+ngx+1,nyp3-ngy-1,l) = chx
         cu(2,j+ngx+1,nyp3-ngy-1,l) = chy
         cu(3,j+ngx+1,nyp3-ngy-1,l) = chz
  100    continue
         do 110 i = 1, 3
         cu(i,1,nyp3-ngy-1,l) = 0.
         cu(i,2,nyp3-ngy-1,l) = 0.
         cu(i,nx+2,nyp3-ngy-1,l) = 0.
         cu(i,nx+3,nyp3-ngy-1,l) = 0.
  110    continue
         cu(1,ngx+2,nyp3-ngy-1,l) = .5*chx
         cu(2,ngx+2,nyp3-ngy-1,l) = .5*chy
         cu(3,ngx+2,nyp3-ngy-1,l) = .5*chz
         cu(1,nx-ngx+2,nyp3-ngy-1,l) = .5*chx
         cu(2,nx-ngx+2,nyp3-ngy-1,l) = .5*chy
         cu(3,nx-ngx+2,nyp3-ngy-1,l) = .5*chz
      endif
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD22X(cu,kstrt,nvp,noff,nyp,xj0,yj0,nx,ny,ngx,ng
     1y,nxe,nypmx,nblok)
c initialize extended non-periodic field
c cu(i,j+1,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l) - 1
c kstrt = starting data block number
c nvp = number of real or virtual processors
c noff(l) = lowermost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c xj0/yj0 = initialization constants in x/y direction
c nx/ny = system length in x/y direction
c ngx/ngy = (0,1) number of grid cells away from edge
c nxe = first dimension of charge array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells
c nblok = number of particle partitions
c quadratic interpolation, for distributed data
      implicit none
      real cu, xj0, yj0
      integer kstrt, nvp, noff, nyp, nx, ny, ngx, ngy, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok), noff(nblok)
      integer i, j, k, l, ks, kk, nyp3, nxg, nx3
      real chx, chy
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      ks = kstrt - 2
      chx = .5*xj0
      chy = .5*yj0
      do 120 l = 1, nblok
      nyp3 = nyp(l) + 3
c handle first grid point in y
      if ((l+ks).eq.0) then
         do 20 j = 1, nx3
         do 10 i = 1, 2
         cu(i,j,2,l) = 0.
   10    continue
   20    continue
      endif
c handle interior grid points in y
      do 70 k = 1, nyp(l)
      kk = k + noff(l)
      if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
         do 30 j = 2, nxg
         cu(1,j+ngx+1,k+1,l) = xj0
         cu(2,j+ngx+1,k+1,l) = yj0
   30    continue
         do 40 i = 1, 2
         cu(i,1,k+1,l) = 0.
         cu(i,2,k+1,l) = 0.
         cu(i,nx+2,k+1,l) = 0.
         cu(i,nx+3,k+1,l) = 0.
   40    continue
         cu(1,ngx+2,k+1,l) = chx
         cu(2,ngx+2,k+1,l) = chy
         cu(1,nx-ngx+2,k+1,l) = chx
         cu(2,nx-ngx+2,k+1,l) = chy
      else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
         do 50 j = 2, nxg
         cu(1,j+ngx+1,k+1,l) = chx
         cu(2,j+ngx+1,k+1,l) = chy
   50    continue
         do 60 i = 1, 2
         cu(i,1,k+1,l) = 0.
         cu(i,2,k+1,l) = 0.
         cu(i,nx+2,k+1,l) = 0.
         cu(i,nx+3,k+1,l) = 0.
   60    continue
         cu(1,ngx+2,k+1,l) = .5*chx
         cu(2,ngx+2,k+1,l) = .5*chy
         cu(1,nx-ngx+2,k+1,l) = .5*chx
         cu(2,nx-ngx+2,k+1,l) = .5*chy
      endif
   70 continue
c guard cells in y
      do 90 j = 1, nx3
      do 80 i = 1, 2
      cu(i,j,1,l) = 0.
      cu(i,j,nyp3-1,l) = 0.
      cu(i,j,nyp3,l) = 0.
   80 continue
   90 continue
c handle last grid point in y
      if (((nyp(l)+noff(l)).lt.(ny-ngy+1)).and.((l+ks).eq.(nvp-1))) then
         do 100 j = 2, nxg
         cu(1,j+ngx+1,nyp3-ngy-1,l) = chx
         cu(2,j+ngx+1,nyp3-ngy-1,l) = chy
  100    continue
         do 110 i = 1, 2
         cu(i,1,nyp3-ngy-1,l) = 0.
         cu(i,2,nyp3-ngy-1,l) = 0.
         cu(i,nx+2,nyp3-ngy-1,l) = 0.
         cu(i,nx+3,nyp3-ngy-1,l) = 0.
  110    continue
         cu(1,ngx+2,nyp3-ngy-1,l) = .5*chx
         cu(2,ngx+2,nyp3-ngy-1,l) = .5*chy
         cu(1,nx-ngx+2,nyp3-ngy-1,l) = .5*chx
         cu(2,nx-ngx+2,nyp3-ngy-1,l) = .5*chy
      endif
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSGUARD2X(q,kstrt,nvp,noff,nyp,qi0,nx,ny,ngx,ngy,nxe,n
     1ypmx,nblok)
c initialize extended non-periodic scalar field
c q(j+1,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c kstrt = starting data block number
c nvp = number of real or virtual processors
c noff(l) = lowermost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c qi0 = initialization constant
c nx/ny = system length in x/y direction
c ngx/ngy = (0,1) number of grid cells away from edge
c nxe = first dimension of charge array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells
c nblok = number of particle partitions
c quadratic interpolation, for distributed data
      implicit none
      real q, qi0
      integer kstrt, nvp, noff, nyp, nx, ny, ngx, ngy, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok), noff(nblok)
      integer j, k, l, ks, kk, nyp3, nxg, nx3
      real qh
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      ks = kstrt - 2
      qh = .5*qi0
      do 70 l = 1, nblok
      nyp3 = nyp(l) + 3
c handle first grid point in y
      if ((l+ks).eq.0) then
         do 10 j = 1, nx3
         q(j,2,l) = 0.
   10    continue
      endif
c handle interior grid points in y
      do 40 k = 1, nyp(l)
      kk = k + noff(l)
      if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
         do 20 j = 2, nxg
         q(j+ngx+1,k+1,l) = qi0
   20    continue
         q(1,k+1,l) = 0.
         q(2,k+1,l) = 0.
         q(nx+2,k+1,l) = 0.
         q(nx+3,k+1,l) = 0.
         q(ngx+2,k+1,l) = qh
         q(nx-ngx+2,k+1,l) = qh
      else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
         do 30 j = 2, nxg
         q(j+ngx+1,k+1,l) = qh
   30    continue
         q(1,k+1,l) = 0.
         q(2,k+1,l) = 0.
         q(nx+2,k+1,l) = 0.
         q(nx+3,k+1,l) = 0.
         q(ngx+2,k+1,l) = .5*qh
         q(nx-ngx+2,k+1,l) = .5*qh
      endif
   40 continue
c guard cells in y
      do 50 j = 1, nx3
      q(j,1,l) = 0.
      q(j,nyp3-1,l) = 0.
      q(j,nyp3,l) = 0.
   50 continue
c handle last grid point in y
      if (((nyp(l)+noff(l)).lt.(ny-ngy+1)).and.((l+ks).eq.(nvp-1))) then
         do 60 j = 2, nxg
         q(j+ngx+1,nyp3-ngy-1,l) = qh
   60    continue
         q(1,nyp3-ngy-1,l) = 0.
         q(2,nyp3-ngy-1,l) = 0.
         q(nx+2,nyp3-ngy-1,l) = 0.
         q(nx+3,nyp3-ngy-1,l) = 0.
         q(ngx+2,nyp3-ngy-1,l) = .5*qh
         q(nx-ngx+2,nyp3-ngy-1,l) = .5*qh
      endif
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD2(cu,kstrt,nvp,nyp,xj0,yj0,zj0,nx,ngx,ngy,nxe,
     1nypmx,nblok)
c initialize extended non-periodic field
c ngx/ngy = (0,1) number of grid cells away from edge
c quadratic interpolation, for distributed data
      implicit none
      real cu, xj0, yj0, zj0
      integer kstrt, nvp, nyp, nx, ngx, ngy, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
      integer i, j, k, l, ks, nyp3, nxg, nx3
      real chx, chy, chz
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      ks = kstrt - 2
      chx = .5*xj0
      chy = .5*yj0
      chz = .5*zj0
      do 120 l = 1, nblok
      nyp3 = nyp(l) + 3
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
      cu(1,ngx+2,k+1,l) = chx
      cu(2,ngx+2,k+1,l) = chy
      cu(3,ngx+2,k+1,l) = chz
      cu(1,nx-ngx+2,k+1,l) = chx
      cu(2,nx-ngx+2,k+1,l) = chy
      cu(3,nx-ngx+2,k+1,l) = chz
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 3
      cu(i,j,1,l) = 0.
      cu(i,j,nyp3-1,l) = 0.
      cu(i,j,nyp3,l) = 0.
   40 continue
   50 continue
      if ((l+ks).eq.0) then
c fix previous deposit when grid is only one wide.
         if ((ngy+4).gt.nyp3) then
            chx = -chx
            chy = -chy
            chz = -chz
         endif
c handle first grid point in y
         do 70 j = 1, nx3
         do 60 i = 1, 3
         cu(i,j,2,l) = 0.
   60    continue
   70    continue
         do 80 j = 2, nxg
         cu(1,j+ngx+1,ngy+2,l) = chx
         cu(2,j+ngx+1,ngy+2,l) = chy
         cu(3,j+ngx+1,ngy+2,l) = chz
   80    continue
         do 90 i = 1, 3
         cu(i,1,ngy+2,l) = 0.
         cu(i,2,ngy+2,l) = 0.
         cu(i,nx+2,ngy+2,l) = 0.
         cu(i,nx+3,ngy+2,l) = 0.
   90    continue
         cu(1,ngx+2,ngy+2,l) = .5*chx
         cu(2,ngx+2,ngy+2,l) = .5*chy
         cu(3,ngx+2,ngy+2,l) = .5*chz
         cu(1,nx-ngx+2,ngy+2,l) = .5*chx
         cu(2,nx-ngx+2,ngy+2,l) = .5*chy
         cu(3,nx-ngx+2,ngy+2,l) = .5*chz
         chx = .5*xj0
         chy = .5*yj0
         chz = .5*zj0
      endif
      if ((l+ks).eq.(nvp-1)) then
c handle last grid point in y
         do 100 j = 2, nxg
         cu(1,j+ngx+1,nyp3-ngy-1,l) = chx
         cu(2,j+ngx+1,nyp3-ngy-1,l) = chy
         cu(3,j+ngx+1,nyp3-ngy-1,l) = chz
  100    continue
         do 110 i = 1, 3
         cu(i,1,nyp3-ngy-1,l) = 0.
         cu(i,2,nyp3-ngy-1,l) = 0.
         cu(i,nx+2,nyp3-ngy-1,l) = 0.
         cu(i,nx+3,nyp3-ngy-1,l) = 0.
  110    continue
         cu(1,ngx+2,nyp3-ngy-1,l) = .5*chx
         cu(2,ngx+2,nyp3-ngy-1,l) = .5*chy
         cu(3,ngx+2,nyp3-ngy-1,l) = .5*chz
         cu(1,nx-ngx+2,nyp3-ngy-1,l) = .5*chx
         cu(2,nx-ngx+2,nyp3-ngy-1,l) = .5*chy
         cu(3,nx-ngx+2,nyp3-ngy-1,l) = .5*chz
      endif
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD22(cu,kstrt,nvp,nyp,xj0,yj0,nx,ngx,ngy,nxe,nyp
     1mx,nblok)
c initialize extended non-periodic field
c ngx/ngy = (0,1) number of grid cells away from edge
c quadratic interpolation, for distributed data
      implicit none
      real cu, xj0, yj0
      integer kstrt, nvp, nyp, nx, ngx, ngy, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
      integer i, j, k, l, ks, nyp3, nxg, nx3
      real chx, chy
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      ks = kstrt - 2
      chx = .5*xj0
      chy = .5*yj0
      do 120 l = 1, nblok
      nyp3 = nyp(l) + 3
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
      cu(1,ngx+2,k+1,l) = chx
      cu(2,ngx+2,k+1,l) = chy
      cu(1,nx-ngx+2,k+1,l) = chx
      cu(2,nx-ngx+2,k+1,l) = chy
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 2
      cu(i,j,1,l) = 0.
      cu(i,j,nyp3-1,l) = 0.
      cu(i,j,nyp3,l) = 0.
   40 continue
   50 continue
      if ((l+ks).eq.0) then
c fix previous deposit when grid is only one wide.
         if ((ngy+4).gt.nyp3) then
            chx = -chx
            chy = -chy
         endif
c handle first grid point in y
         do 70 j = 1, nx3
         do 60 i = 1, 2
         cu(i,j,2,l) = 0.
   60    continue
   70    continue
         do 80 j = 2, nxg
         cu(1,j+ngx+1,ngy+2,l) = chx
         cu(2,j+ngx+1,ngy+2,l) = chy
   80    continue
         do 90 i = 1, 2
         cu(i,1,ngy+2,l) = 0.
         cu(i,2,ngy+2,l) = 0.
         cu(i,nx+2,ngy+2,l) = 0.
         cu(i,nx+3,ngy+2,l) = 0.
   90    continue
         cu(1,ngx+2,ngy+2,l) = .5*chx
         cu(2,ngx+2,ngy+2,l) = .5*chy
         cu(1,nx-ngx+2,ngy+2,l) = .5*chx
         cu(2,nx-ngx+2,ngy+2,l) = .5*chy
         chx = .5*xj0
         chy = .5*yj0
      endif
      if ((l+ks).eq.(nvp-1)) then
c handle last grid point in y
         do 100 j = 2, nxg
         cu(1,j+ngx+1,nyp3-ngy-1,l) = chx
         cu(2,j+ngx+1,nyp3-ngy-1,l) = chy
  100    continue
         do 110 i = 1, 2
         cu(i,1,nyp3-ngy-1,l) = 0.
         cu(i,2,nyp3-ngy-1,l) = 0.
         cu(i,nx+2,nyp3-ngy-1,l) = 0.
         cu(i,nx+3,nyp3-ngy-1,l) = 0.
  110    continue
         cu(1,ngx+2,nyp3-ngy-1,l) = .5*chx
         cu(2,ngx+2,nyp3-ngy-1,l) = .5*chy
         cu(1,nx-ngx+2,nyp3-ngy-1,l) = .5*chx
         cu(2,nx-ngx+2,nyp3-ngy-1,l) = .5*chy
      endif
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSGUARD2(q,kstrt,nvp,nyp,qi0,nx,ngx,ngy,nxe,nypmx,nblo
     1k)
c initialize extended non-periodic scalar field
c ngx/ngy = (0,1) number of grid cells away from edge
c quadratic interpolation, for distributed data
      implicit none
      real q, qi0
      integer kstrt, nvp, nyp, nx, ngx, ngy, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
      integer j, k, l, ks, nyp3, nxg, nx3
      real qh
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      ks = kstrt - 2
      qh = .5*qi0
      do 70 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 20 k = 1, nyp(l)
      do 10 j = 2, nxg
      q(j+ngx+1,k+1,l) = qi0
   10 continue
      q(1,k+1,l) = 0.
      q(2,k+1,l) = 0.
      q(nx+2,k+1,l) = 0.
      q(nx+3,k+1,l) = 0.
      q(ngx+2,k+1,l) = qh
      q(nx-ngx+2,k+1,l) = qh
   20 continue
      do 30 j = 1, nx3
      q(j,1,l) = 0.
      q(j,nyp3-1,l) = 0.
      q(j,nyp3,l) = 0.
   30 continue
      if ((l+ks).eq.0) then
c fix previous deposit when grid is only one wide.
         if ((ngy+4).gt.nyp3) qh = -qh
c handle first grid point in y
         do 40 j = 1, nx3
         q(j,2,l) = 0.
   40    continue
         do 50 j = 2, nxg
         q(j+ngx+1,ngy+2,l) = qh
   50    continue
         q(1,ngy+2,l) = 0.
         q(2,ngy+2,l) = 0.
         q(nx+2,ngy+2,l) = 0.
         q(nx+3,ngy+2,l) = 0.
         q(ngx+2,ngy+2,l) = .5*qh
         q(nx-ngx+2,ngy+2,l) = .5*qh
         qh = .5*qi0
      endif
      if ((l+ks).eq.(nvp-1)) then
c handle last grid point in y
         do 60 j = 2, nxg
         q(j+ngx+1,nyp3-ngy-1,l) = qh
   60    continue
         q(1,nyp3-ngy-1,l) = 0.
         q(2,nyp3-ngy-1,l) = 0.
         q(nx+2,nyp3-ngy-1,l) = 0.
         q(nx+3,nyp3-ngy-1,l) = 0.
         q(ngx+2,nyp3-ngy-1,l) = .5*qh
         q(nx-ngx+2,nyp3-ngy-1,l) = .5*qh
      endif
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARD2X(cu,nyp,nx,nxe,nypmx,nblok)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp3
c add up guard cells
      do 30 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 20 k = 1, nyp3
      do 10 i = 1, 3
      cu(i,2,k,l) = cu(i,2,k,l) + 2.*cu(i,1,k,l)
      cu(i,3,k,l) = cu(i,3,k,l) - cu(i,1,k,l)
      cu(i,nx+1,k,l) = cu(i,nx+1,k,l) - cu(i,nx+3,k,l)
      cu(i,nx+2,k,l) = cu(i,nx+2,k,l) + 2.*cu(i,nx+3,k,l)
      cu(i,1,k,l) = 0.
      cu(i,nx+3,k,l) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARD22X(cu,nyp,nx,nxe,nypmx,nblok)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c for distributed data
      implicit none
      real cu
      integer nyp, nx, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
c local data
      integer i, k, l, nyp3
c add up guard cells
      do 30 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 20 k = 1, nyp3
      do 10 i = 1, 2
      cu(i,2,k,l) = cu(i,2,k,l) + 2.*cu(i,1,k,l)
      cu(i,3,k,l) = cu(i,3,k,l) - cu(i,1,k,l)
      cu(i,nx+1,k,l) = cu(i,nx+1,k,l) - cu(i,nx+3,k,l)
      cu(i,nx+2,k,l) = cu(i,nx+2,k,l) + 2.*cu(i,nx+3,k,l)
      cu(i,1,k,l) = 0.
      cu(i,nx+3,k,l) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLAGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c for distributed scalar data
      implicit none
      real q
      integer nyp, nx, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
c local data
      integer k, l, nyp3
c add up guard cells
      do 20 l = 1, nblok
      nyp3 = nyp(l) + 3
      do 10 k = 1, nyp3
      q(2,k,l) = q(2,k,l) + 2.*q(1,k,l)
      q(3,k,l) = q(3,k,l) - q(1,k,l)
      q(nx+1,k,l) = q(nx+1,k,l) - q(nx+3,k,l)
      q(nx+2,k,l) = q(nx+2,k,l) + 2.*q(nx+3,k,l)
      q(1,k,l) = 0.
      q(nx+3,k,l) = 0.
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD2XL(cu,kstrt,nvp,noff,nyp,xj0,yj0,zj0,nx,ny,ng
     1x,ngy,nxe,nypmx,nblok)
c initialize extended non-periodic field
c cu(i,j,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l)
c kstrt = starting data block number
c nvp = number of real or virtual processors
c noff(l) = lowermost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c xj0/yj0/zj0 = initialization constants in x/y/z direction
c nx/ny = system length in x/y direction
c ngx/ngy = (0,1) number of grid cells away from edge
c nxe = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells
c nblok = number of particle partitions
c linear interpolation, for distributed data
      implicit none
      real cu, xj0, yj0, zj0
      integer kstrt, nvp, noff, nyp, nx, ny, ngx, ngy, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok), noff(nblok)
      integer i, j, k, l, ks, kk, nyp1, nxg, nx1
      real chx, chy, chz
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      ks = kstrt - 2
      chx = .5*xj0
      chy = .5*yj0
      chz = .5*zj0
      do 120 l = 1, nblok
      nyp1 = nyp(l) + 1
c handle first grid point in y
      if ((l+ks).eq.0) then
         do 20 j = 1, nx1
         do 10 i = 1, 3
         cu(i,j,1,l) = 0.
   10    continue
   20    continue
      endif
c handle interior grid points in y
      do 70 k = 1, nyp(l)
      kk = k + noff(l)
      if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
         do 30 j = 2, nxg
         cu(1,j+ngx,k,l) = xj0
         cu(2,j+ngx,k,l) = yj0
         cu(3,j+ngx,k,l) = zj0
   30    continue
         do 40 i = 1, 3
         cu(i,1,k,l) = 0.
         cu(i,nx+1,k,l) = 0.
   40    continue
         cu(1,ngx+1,k,l) = chx
         cu(2,ngx+1,k,l) = chy
         cu(3,ngx+1,k,l) = chz
         cu(1,nx-ngx+1,k,l) = chx
         cu(2,nx-ngx+1,k,l) = chy
         cu(3,nx-ngx+1,k,l) = chz
      else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
         do 50 j = 2, nxg
         cu(1,j+ngx,k,l) = chx
         cu(2,j+ngx,k,l) = chy
         cu(3,j+ngx,k,l) = chz
   50    continue
         do 60 i = 1, 3
         cu(i,1,k,l) = 0.
         cu(i,nx+1,k,l) = 0.
   60    continue
         cu(1,ngx+1,k,l) = .5*chx
         cu(2,ngx+1,k,l) = .5*chy
         cu(3,ngx+1,k,l) = .5*chz
         cu(1,nx-ngx+1,k,l) = .5*chx
         cu(2,nx-ngx+1,k,l) = .5*chy
         cu(3,nx-ngx+1,k,l) = .5*chz
      endif
   70 continue
c guard cells in y
      do 90 j = 1, nx1
      do 80 i = 1, 3
      cu(i,j,nyp1,l) = 0.
   80 continue
   90 continue
c handle last grid point in y
      if (((nyp(l)+noff(l)).lt.(ny-ngy+1)).and.((l+ks).eq.(nvp-1))) then
         do 100 j = 2, nxg
         cu(1,j+ngx,nyp1-ngy,l) = chx
         cu(2,j+ngx,nyp1-ngy,l) = chy
         cu(3,j+ngx,nyp1-ngy,l) = chz
  100    continue
         do 110 i = 1, 3
         cu(i,1,nyp1-ngy,l) = 0.
         cu(i,nx+1,nyp1-ngy,l) = 0.
  110    continue
         cu(1,ngx+1,nyp1-ngy,l) = .5*chx
         cu(2,ngx+1,nyp1-ngy,l) = .5*chy
         cu(3,ngx+1,nyp1-ngy,l) = .5*chz
         cu(1,nx-ngx+1,nyp1-ngy,l) = .5*chx
         cu(2,nx-ngx+1,nyp1-ngy,l) = .5*chy
         cu(3,nx-ngx+1,nyp1-ngy,l) = .5*chz
      endif
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD22XL(cu,kstrt,nvp,noff,nyp,xj0,yj0,nx,ny,ngx,n
     1gy,nxe,nypmx,nblok)
c initialize extended non-periodic field
c cu(i,j,k,l) = ith component of current density at grid point (j,kk),
c where kk = k + noff(l)
c kstrt = starting data block number
c nvp = number of real or virtual processors
c noff(l) = lowermost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c xj0/yj0 = initialization constants in x/y/ direction
c nx/ny = system length in x/y direction
c ngx/ngy = (0,1) number of grid cells away from edge
c nxe = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells
c nblok = number of particle partitions
c linear interpolation, for distributed data
      implicit none
      real cu, xj0, yj0
      integer kstrt, nvp, noff, nyp, nx, ny, ngx, ngy, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok), noff(nblok)
      integer i, j, k, l, ks, kk, nyp1, nxg, nx1
      real chx, chy
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      ks = kstrt - 2
      chx = .5*xj0
      chy = .5*yj0
      do 120 l = 1, nblok
      nyp1 = nyp(l) + 1
c handle first grid point in y
      if ((l+ks).eq.0) then
         do 20 j = 1, nx1
         do 10 i = 1, 2
         cu(i,j,1,l) = 0.
   10    continue
   20    continue
      endif
c handle interior grid points in y
      do 70 k = 1, nyp(l)
      kk = k + noff(l)
      if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
         do 30 j = 2, nxg
         cu(1,j+ngx,k,l) = xj0
         cu(2,j+ngx,k,l) = yj0
   30    continue
         do 40 i = 1, 2
         cu(i,1,k,l) = 0.
         cu(i,nx+1,k,l) = 0.
   40    continue
         cu(1,ngx+1,k,l) = chx
         cu(2,ngx+1,k,l) = chy
         cu(1,nx-ngx+1,k,l) = chx
         cu(2,nx-ngx+1,k,l) = chy
      else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
         do 50 j = 2, nxg
         cu(1,j+ngx,k,l) = chx
         cu(2,j+ngx,k,l) = chy
   50    continue
         do 60 i = 1, 2
         cu(i,1,k,l) = 0.
         cu(i,nx+1,k,l) = 0.
   60    continue
         cu(1,ngx+1,k,l) = .5*chx
         cu(2,ngx+1,k,l) = .5*chy
         cu(1,nx-ngx+1,k,l) = .5*chx
         cu(2,nx-ngx+1,k,l) = .5*chy
      endif
   70 continue
c guard cells in y
      do 90 j = 1, nx1
      do 80 i = 1, 2
      cu(i,j,nyp1,l) = 0.
   80 continue
   90 continue
c handle last grid point in y
      if (((nyp(l)+noff(l)).lt.(ny-ngy+1)).and.((l+ks).eq.(nvp-1))) then
         do 100 j = 2, nxg
         cu(1,j+ngx,nyp1-ngy,l) = chx
         cu(2,j+ngx,nyp1-ngy,l) = chy
  100    continue
         do 110 i = 1, 2
         cu(i,1,nyp1-ngy,l) = 0.
         cu(i,nx+1,nyp1-ngy,l) = 0.
  110    continue
         cu(1,ngx+1,nyp1-ngy,l) = .5*chx
         cu(2,ngx+1,nyp1-ngy,l) = .5*chy
         cu(1,nx-ngx+1,nyp1-ngy,l) = .5*chx
         cu(2,nx-ngx+1,nyp1-ngy,l) = .5*chy
      endif
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSGUARD2XL(q,kstrt,nvp,noff,nyp,qi0,nx,ny,ngx,ngy,nxe,
     1nypmx,nblok)
c initialize extended non-periodic scalar field
c q(j,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l)
c kstrt = starting data block number
c nvp = number of real or virtual processors
c noff(l) = lowermost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c qi0 = initialization constant
c nx/ny = system length in x/y directionn
c ngx/ngy = (0,1) number of grid cells away from edge
c nxe = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells
c nblok = number of particle partitions
c linear interpolation, for distributed data
      implicit none
      real q, qi0
      integer kstrt, nvp, noff, nyp, nx, ny, ngx, ngy, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok), noff(nblok)
      integer j, k, l, ks, kk, nyp1, nxg, nx1
      real qh
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      ks = kstrt - 2
      qh = .5*qi0
      do 70 l = 1, nblok
      nyp1 = nyp(l) + 1
c handle first grid point in y
      if ((l+ks).eq.0) then
         do 10 j = 1, nx1
         q(j,1,l) = 0.
   10    continue
      endif
c handle interior grid points in y
      do 40 k = 1, nyp(l)
      kk = k + noff(l)
      if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
         do 20 j = 2, nxg
         q(j+ngx,k,l) = qi0
   20    continue
         q(1,k,l) = 0.
         q(nx+1,k,l) = 0.
         q(ngx+1,k,l) = qh
         q(nx-ngx+1,k,l) = qh
      else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
         do 30 j = 2, nxg
         q(j+ngx,k,l) = qh
   30    continue
         q(1,k,l) = 0.
         q(nx+1,k,l) = 0.
         q(ngx+1,k,l) = .5*qh
         q(nx-ngx+1,k,l) = .5*qh
      endif
   40 continue
c guard cells in y
      do 50 j = 1, nx1
      q(j,nyp1,l) = 0.
   50 continue
c handle last grid point in y
      if (((nyp(l)+noff(l)).lt.(ny-ngy+1)).and.((l+ks).eq.(nvp-1))) then
         do 60 j = 2, nxg
         q(j+ngx,nyp1-ngy,l) = qh
   60    continue
         q(1,nyp1-ngy,l) = 0.
         q(nx+1,nyp1-ngy,l) = 0.
         q(ngx+1,nyp1-ngy,l) = .5*qh
         q(nx-ngx+1,nyp1-ngy,l) = .5*qh
      endif
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD2L(cu,kstrt,nvp,nyp,xj0,yj0,zj0,nx,ngx,ngy,nxe
     1,nypmx,nblok)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c linear interpolation, for distributed data
      implicit none
      real cu, xj0, yj0, zj0
      integer kstrt, nvp, nyp, nx, ngx, ngy, nxe, nypmx, nblok
      dimension cu(3,nxe,nypmx,nblok), nyp(nblok)
      integer i, j, k, l, ks, nyp1, nxg, nx1
      real chx, chy, chz
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      ks = kstrt - 2
      chx = .5*xj0
      chy = .5*yj0
      chz = .5*zj0
      do 120 l = 1, nblok
      nyp1 = nyp(l) + 1
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
      cu(1,ngx+1,k,l) = chx
      cu(2,ngx+1,k,l) = chy
      cu(3,ngx+1,k,l) = chz
      cu(1,nx-ngx+1,k,l) = chx
      cu(2,nx-ngx+1,k,l) = chy
      cu(3,nx-ngx+1,k,l) = chz
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 3
      cu(i,j,nyp1,l) = 0.
   40 continue
   50 continue
      if ((l+ks).eq.0) then
c fix previous deposit when grid is only one wide.
         if ((ngy+2).gt.nyp1) then
            chx = -chx
            chy = -chy
            chz = -chz
         endif
c handle first grid point in y
         do 70 j = 1, nx1
         do 60 i = 1, 3
         cu(i,j,1,l) = 0.
   60    continue
   70    continue
         do 80 j = 2, nxg
         cu(1,j+ngx,ngy+1,l) = chx
         cu(2,j+ngx,ngy+1,l) = chy
         cu(3,j+ngx,ngy+1,l) = chz
   80    continue
         do 90 i = 1, 3
         cu(i,1,ngy+1,l) = 0.
         cu(i,nx+1,ngy+1,l) = 0.
   90    continue
         cu(1,ngx+1,ngy+1,l) = .5*chx
         cu(2,ngx+1,ngy+1,l) = .5*chy
         cu(3,ngx+1,ngy+1,l) = .5*chz
         cu(1,nx-ngx+1,ngy+1,l) = .5*chx
         cu(2,nx-ngx+1,ngy+1,l) = .5*chy
         cu(3,nx-ngx+1,ngy+1,l) = .5*chz
         chx = .5*xj0
         chy = .5*yj0
         chz = .5*zj0
      endif
      if ((l+ks).eq.(nvp-1)) then
c handle last grid point in y
         do 100 j = 2, nxg
         cu(1,j+ngx,nyp1-ngy,l) = chx
         cu(2,j+ngx,nyp1-ngy,l) = chy
         cu(3,j+ngx,nyp1-ngy,l) = chz
  100    continue
         do 110 i = 1, 3
         cu(i,1,nyp1-ngy,l) = 0.
         cu(i,nx+1,nyp1-ngy,l) = 0.
  110    continue
         cu(1,ngx+1,nyp1-ngy,l) = .5*chx
         cu(2,ngx+1,nyp1-ngy,l) = .5*chy
         cu(3,ngx+1,nyp1-ngy,l) = .5*chz
         cu(1,nx-ngx+1,nyp1-ngy,l) = .5*chx
         cu(2,nx-ngx+1,nyp1-ngy,l) = .5*chy
         cu(3,nx-ngx+1,nyp1-ngy,l) = .5*chz
      endif
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD22L(cu,kstrt,nvp,nyp,xj0,yj0,nx,ngx,ngy,nxe,ny
     1pmx,nblok)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c linear interpolation, for distributed data
      implicit none
      real cu, xj0, yj0
      integer kstrt, nvp, nyp, nx, ngx, ngy, nxe, nypmx, nblok
      dimension cu(2,nxe,nypmx,nblok), nyp(nblok)
      integer i, j, k, l, ks, nyp1, nxg, nx1
      real chx, chy
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      ks = kstrt - 2
      chx = .5*xj0
      chy = .5*yj0
      do 120 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 30 k = 1, nyp(l)
      do 10 j = 2, nxg
      cu(1,j+ngx,k,l) = xj0
      cu(2,j+ngx,k,l) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k,l) = 0.
      cu(i,nx+1,k,l) = 0.
   20 continue
      cu(1,ngx+1,k,l) = chx
      cu(2,ngx+1,k,l) = chy
      cu(1,nx-ngx+1,k,l) = chx
      cu(2,nx-ngx+1,k,l) = chy
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 2
      cu(i,j,nyp1,l) = 0.
   40 continue
   50 continue
      if ((l+ks).eq.0) then
c fix previous deposit when grid is only one wide.
         if ((ngy+2).gt.nyp1) then
            chx = -chx
            chy = -chy
         endif
c handle first grid point in y
         do 70 j = 1, nx1
         do 60 i = 1, 2
         cu(i,j,1,l) = 0.
   60    continue
   70    continue
         do 80 j = 2, nxg
         cu(1,j+ngx,ngy+1,l) = chx
         cu(2,j+ngx,ngy+1,l) = chy
   80    continue
         do 90 i = 1, 2
         cu(i,1,ngy+1,l) = 0.
         cu(i,nx+1,ngy+1,l) = 0.
   90    continue
         cu(1,ngx+1,ngy+1,l) = .5*chx
         cu(2,ngx+1,ngy+1,l) = .5*chy
         cu(1,nx-ngx+1,ngy+1,l) = .5*chx
         cu(2,nx-ngx+1,ngy+1,l) = .5*chy
         chx = .5*xj0
         chy = .5*yj0
      endif
      if ((l+ks).eq.(nvp-1)) then
c handle last grid point in y
         do 100 j = 2, nxg
         cu(1,j+ngx,nyp1-ngy,l) = chx
         cu(2,j+ngx,nyp1-ngy,l) = chy
  100    continue
         do 110 i = 1, 2
         cu(i,1,nyp1-ngy,l) = 0.
         cu(i,nx+1,nyp1-ngy,l) = 0.
  110    continue
         cu(1,ngx+1,nyp1-ngy,l) = .5*chx
         cu(2,ngx+1,nyp1-ngy,l) = .5*chy
         cu(1,nx-ngx+1,nyp1-ngy,l) = .5*chx
         cu(2,nx-ngx+1,nyp1-ngy,l) = .5*chy
      endif
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSGUARD2L(q,kstrt,nvp,nyp,qi0,nx,ngx,ngy,nxe,nypmx,nbl
     1ok)
c initialize extended non-periodic scalar field
c ngx/ngy = (0,1) = number of grid cells away from edge
c linear interpolation, for distributed data
      implicit none
      real q, qi0
      integer kstrt, nvp, nyp, nx, ngx, ngy, nxe, nypmx, nblok
      dimension q(nxe,nypmx,nblok), nyp(nblok)
      integer j, k, l, ks, nyp1, nxg, nx1
      real qh
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      ks = kstrt - 2
      qh = .5*qi0
      do 70 l = 1, nblok
      nyp1 = nyp(l) + 1
      do 20 k = 1, nyp(l)
      do 10 j = 2, nxg
      q(j+ngx,k,l) = qi0
   10 continue
      q(1,k,l) = 0.
      q(nx+1,k,l) = 0.
      q(ngx+1,k,l) = qh
      q(nx-ngx+1,k,l) = qh
   20 continue
      do 30 j = 1, nx1
      q(j,nyp1,l) = 0.
   30 continue
      if ((l+ks).eq.0) then
c fix previous deposit when grid is only one wide.
         if ((ngy+2).gt.nyp1) qh = -qh
c handle first grid point in y
         do 40 j = 1, nx1
         q(j,1,l) = 0.
   40    continue
         do 50 j = 2, nxg
         q(j+ngx,ngy+1,l) = qh
   50    continue
         q(1,ngy+1,l) = 0.
         q(nx+1,ngy+1,l) = 0.
         q(ngx+1,ngy+1,l) = .5*qh
         q(nx-ngx+1,ngy+1,l) = .5*qh
         qh = .5*qi0
      endif
      if ((l+ks).eq.(nvp-1)) then
c handle last grid point in y
         do 60 j = 2, nxg
         q(j+ngx,nyp1-ngy,l) = qh
   60    continue
         q(1,nyp1-ngy,l) = 0.
         q(nx+1,nyp1-ngy,l) = 0.
         q(ngx+1,nyp1-ngy,l) = .5*qh
         q(nx-ngx+1,nyp1-ngy,l) = .5*qh
      endif
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny
     12d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function, with dirichlet
c boundary conditions (zero potential), for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-sin, sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign = -1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,fy,we
c approximate flop count is: 11*nx*ny
c for isign = 1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,we
c approximate flop count is: 6*nx*ny
c for isign = 2, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fy
c approximate flop count is: 2*nx*ny
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fx(k,j,l) = x component of complex force/charge,
c fy(k,j,l) = y component of complex force/charge,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex q, fx, fy, ffd, zero
      dimension q(ny2d,kxp2,j2blok)
      dimension fx(ny2d,kxp2,j2blok), fy(ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
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
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
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
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*real(q(k,j,l))
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fx(k,j,l) = cmplx(0.,at2)
         fx(k1,j,l) = cmplx(0.,-at2)
         fy(k,j,l) = cmplx(0.,at3)
         fy(k1,j,l) = cmplx(0.,at3)
         wp = wp + at1*real(q(k,j,l))**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fx(1,j,l) = zero
      fx(ny+1,j,l) = zero
      fy(1,j,l) = zero
      fy(ny+1,j,l) = zero
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         fx(k,1,l) = zero
         fx(k1,1,l) = zero
         fy(k,1,l) = zero
         fy(k1,1,l) = zero
   70    continue
      endif
   80 continue
   90 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
  100 if (isign.gt.1) go to 160
      wp = 0.0d0
      if (kstrt.gt.nx) go to 150
      do 140 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, ny
         k1 = ny2 - k
         at2 = real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         at3 = at2*real(q(k,j,l))
         fx(k,j,l) = cmplx(at3,0.)
         fx(k1,j,l) = cmplx(-at3,0.)
         wp = wp + at1*real(q(k,j,l))**2
  110    continue
      endif
c mode numbers ky = 0, ny
      fx(1,j,l) = zero
      fx(ny+1,j,l) = zero
  120 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 130 k = 2, ny
         k1 = ny2 - k
         fx(k,1,l) = zero
         fx(k1,1,l) = zero
  130    continue
      endif
  140 continue
  150 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  160 if (kstrt.gt.nx) go to 210
      do 200 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,j,l))
         at2 = at1*real(q(k,j,l))
         fy(k,j,l) = cmplx(at2,0.)
         fy(k1,j,l) = cmplx(-at2,0.)
  170    continue
      endif
c mode numbers ky = 0, ny
      fy(1,j,l) = zero
      fy(ny+1,j,l) = zero
  180 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 190 k = 2, ny
         k1 = ny2 - k
         fy(k,1,l) = zero
         fy(k1,1,l) = zero
  190    continue
      endif
  200 continue
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISD2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny
     1v,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function, with dirichlet
c boundary conditions (zero potential), for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-sin, sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign = -1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,fy,we
c approximate flop count is: 10*nx*ny
c for isign = 1, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fx,we
c approximate flop count is: 5*nx*ny
c for isign = 2, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fy
c approximate flop count is: 1*nx*ny
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fx(k,j,l) = x component of complex force/charge,
c fy(k,j,l) = y component of complex force/charge,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension q(nyv,kxp2+1,j2blok)
      dimension fx(nyv,kxp2+1,j2blok), fy(nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 110
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*q(k,j,l)
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fx(k,j,l) = at2
         fy(k,j,l) = at3
         wp = wp + at1*q(k,j,l)**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fx(1,j,l) = 0.
      fx(ny+1,j,l) = 0.
      fy(1,j,l) = 0.
      fy(ny+1,j,l) = 0.
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         fx(k,1,l) = 0.
         fy(k,1,l) = 0.
   70    continue
      endif
      do 80 k = 1, ny1
      fx(k,kxp2+1,l) = 0.
      fy(k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
  110 if (isign.gt.1) go to 180
      wp = 0.0d0
      if (kstrt.gt.nx) go to 170
      do 160 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 130 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 120 k = 2, ny
         at2 = real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         at3 = at2*q(k,j,l)
         fx(k,j,l) = at3
         wp = wp + at1*q(k,j,l)**2
  120    continue
      endif
c mode numbers ky = 0, ny
      fx(1,j,l) = 0.
      fx(ny+1,j,l) = 0.
  130 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 140 k = 2, ny
         fx(k,1,l) = 0.
  140    continue
      endif
      do 150 k = 1, ny1
      fx(k,kxp2+1,l) = 0.
  150 continue
  160 continue
  170 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  180 if (kstrt.gt.nx) go to 240
      do 230 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 200 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 190 k = 2, ny
         at1 = aimag(ffd(k,j,l))
         at2 = at1*q(k,j,l)
         fy(k,j,l) = at2
  190    continue
      endif
c mode numbers ky = 0, ny
      fy(1,j,l) = 0.
      fy(ny+1,j,l) = 0.
  200 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 210 k = 2, ny
         fy(k,1,l) = 0.
  210    continue
      endif
      do 220 k = 1, ny1
      fy(k,kxp2+1,l) = 0.
  220 continue
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISDX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2
     1d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fxy,we
c approximate flop count is: 11*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex q, fxy, ffd, zero
      dimension q(ny2d,kxp2,j2blok), fxy(2,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
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
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*real(q(k,j,l))
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fxy(1,k,j,l) = cmplx(0.,at2)
         fxy(2,k,j,l) = cmplx(0.,at3)
         fxy(1,k1,j,l) = cmplx(0.,-at2)
         fxy(2,k1,j,l) = cmplx(0.,at3)
         wp = wp + at1*real(q(k,j,l))**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fxy(1,1,j,l) = zero
      fxy(2,1,j,l) = zero
      fxy(1,ny+1,j,l) = zero
      fxy(2,ny+1,j,l) = zero
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = zero
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
   70    continue
      endif
   80 continue
   90 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,
     1kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fxy,we
c approximate flop count is: 10*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension q(nyv,kxp2+1,j2blok), fxy(2,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*q(k,j,l)
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fxy(1,k,j,l) = at2
         fxy(2,k,j,l) = at3
         wp = wp + at1*q(k,j,l)**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fxy(1,1,j,l) = 0.
      fxy(2,1,j,l) = 0.
      fxy(1,ny+1,j,l) = 0.
      fxy(2,ny+1,j,l) = 0.
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         fxy(1,k,1,l) = 0.
         fxy(2,k,1,l) = 0.
   70    continue
      endif
      do 80 k = 1, ny1
      fxy(1,k,kxp2+1,l) = 0.
      fxy(2,k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISDX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,ny2
     1d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential),
c for distributed data.  Zeros out z component
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fxy,we
c approximate flop count is: 11*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c fxy(3,k,j,l) = zero,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex q, fxy, ffd, zero
      dimension q(ny2d,kxp2,j2blok), fxy(3,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
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
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 90
      do 80 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*real(q(k,j,l))
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fxy(1,k,j,l) = cmplx(0.,at2)
         fxy(2,k,j,l) = cmplx(0.,at3)
         fxy(3,k,j,l) = zero
         fxy(1,k1,j,l) = cmplx(0.,-at2)
         fxy(2,k1,j,l) = cmplx(0.,at3)
         fxy(3,k1,j,l) = zero
         wp = wp + at1*real(q(k,j,l))**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fxy(1,1,j,l) = zero
      fxy(2,1,j,l) = zero
      fxy(3,1,j,l) = zero
      fxy(1,ny+1,j,l) = zero
      fxy(2,ny+1,j,l) = zero
      fxy(3,ny+1,j,l) = zero
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         fxy(1,k,1,l) = zero
         fxy(2,k,1,l) = zero
         fxy(3,k,1,l) = zero
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
         fxy(3,k1,1,l) = zero
   70    continue
      endif
   80 continue
   90 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,kstrt,nyv,
     1kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential),
c for distributed data.  Zeros out z component
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0,input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c               output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c                output: fxy,we
c approximate flop count is: 10*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c q(k,j,l) = complex charge density for fourier mode (jj-1,k-1)
c fxy(1,k,j,l) = x component of complex force/charge,
c fxy(2,k,j,l) = y component of complex force/charge,
c fxy(3,k,j,l) = zero,
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c if isign = 0, form factor array is prepared
c ffd(k,2*j,l) = finite-size particle shape factor s
c ffd(k,2*j-1,l) = potential green's function g
c all for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension q(nyv,kxp2+1,j2blok), fxy(3,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 40
      if (kstrt.gt.nx) return
c prepare form factor array
      do 30 l = 1, j2blok
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
c calculate force/charge and sum field energy
   40 wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         at1 = real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at3 = -at1*q(k,j,l)
         at2 = dkx*at3
         at3 = dny*float(k - 1)*at3
         fxy(1,k,j,l) = at2
         fxy(2,k,j,l) = at3
         fxy(3,k,j,l) = 0.
         wp = wp + at1*q(k,j,l)**2
   50    continue
      endif
c mode numbers ky = 0, ny
      fxy(1,1,j,l) = 0.
      fxy(2,1,j,l) = 0.
      fxy(3,1,j,l) = 0.
      fxy(1,ny+1,j,l) = 0.
      fxy(2,ny+1,j,l) = 0.
      fxy(3,ny+1,j,l) = 0.
   60 continue
c mode number kx = 0
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         fxy(1,k,1,l) = 0.
         fxy(2,k,1,l) = 0.
         fxy(3,k,1,l) = 0.
   70    continue
      endif
      do 80 k = 1, ny1
      fxy(1,k,kxp2+1,l) = 0.
      fxy(2,k,kxp2+1,l) = 0.
      fxy(3,k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PDIVFD2(f,df,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
c this subroutine calculates the divergence in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the charge density from the electric field
c input: all except df, output: df
c approximate flop count is: 6*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny - 1, and nvp = number of procs
c the divergence is calculated using the equation:
c df(kx,ky) = sqrt(-1)*(kx*fx(kx,ky)+ky*fy(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for df(kx=pi) = df(ky=pi) = df(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
      real f, df
      dimension f(ndim,nyv,kxp2+1,j2blok), df(nyv,kxp2+1,j2blok)
      if (ndim.lt.2) return
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the divergence
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         df(k,j,l) = -(dkx*f(1,k,j,l) + dky*f(2,k,j,l))
   10    continue
      endif
c mode numbers ky = 0, ny
      df(1,j,l) = 0.
      df(ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         df(k,1,l) = 0.
   30    continue
      endif
      do 40 k = 1, ny1
      df(k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRADFD2(df,f,nx,ny,kstrt,ndim,nyv,kxp2,j2blok)
c this subroutine calculates the gradient in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the electric field from the potential
c input: all except f, output: f
c approximate flop count is: 4*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny - 1, and nvp = number of procs
c the gradient is calculated using the equations:
c fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
c fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for fx(kx=pi) = fy(kx=pi) = 0, fx(ky=pi) = fy(ky=pi) = 0,
c and fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
      real df, f
      dimension df(nyv,kxp2+1,j2blok), f(ndim,nyv,kxp2+1,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the gradient
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         f(1,k,j,l) = dkx*df(k,j,l)
         f(2,k,j,l) = dky*df(k,j,l)
   10    continue
      endif
c mode numbers ky = 0, ny
      f(1,1,j,l) = 0.
      f(2,1,j,l) = 0.
      f(1,ny+1,j,l) = 0.
      f(2,ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         f(1,k,1,l) = 0.
         f(2,k,1,l) = 0.
   30    continue
      endif
      do 40 k = 1, ny1
      f(1,k,kxp2+1,l) = 0.
      f(2,k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCURLFD2(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates the curl in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the magnetic field from the vector potential
c input: all except g, output: g
c approximate flop count is: 8*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny - 1, and nvp = number of procs
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
      real f, g
      dimension f(3,nyv,kxp2+1,j2blok), g(3,nyv,kxp2+1,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the curl
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         g(1,k,j,l) = dky*f(3,k,j,l)
         g(2,k,j,l) = -dkx*f(3,k,j,l)
         g(3,k,j,l) = dkx*f(2,k,j,l) - dky*f(1,k,j,l)
   10    continue
c mode numbers ky = 0, ny
         g(1,1,j,l) = 0.
         g(2,1,j,l) = 0.
         g(3,1,j,l) = dkx*f(2,1,j,l)
      endif
      g(1,ny+1,j,l) = 0.
      g(2,ny+1,j,l) = 0.
      g(3,ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         g(1,k,1,l) = 0.
         g(2,k,1,l) = 0.
         g(3,k,1,l) = -dky*f(1,k,1,l)
   30    continue
         g(1,1,1,l) = 0.
         g(2,1,1,l) = 0.
         g(3,1,1,l) = 0.
      endif
      do 40 k = 1, ny1
      g(1,k,kxp2+1,l) = 0.
      g(2,k,kxp2+1,l) = 0.
      g(3,k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCURLFD22(f,g,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates the curl in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the magnetic field from the vector potential
c input: all except g, output: g
c approximate flop count is: 32*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny - 1, and nvp = number of procs
c the curl is calculated using the equations:
c g(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c kxp2 = number of data values per block
c j2blok = number of data blocks
      real f, g
      dimension f(2,nyv,kxp2+1,j2blok), g(nyv,kxp2+1,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the curl
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         g(k,j,l) = dkx*f(2,k,j,l) - dky*f(1,k,j,l)
   10    continue
c mode numbers ky = 0, ny
         g(1,j,l) = dkx*f(2,1,j,l)
      endif
      g(ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         g(k,1,l) = -dky*f(1,k,1,l)
   30    continue
         g(1,1,l) = 0.
      endif
      do 40 k = 1, ny1
      g(k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPDX2(cu,nx,ny,kstrt,ny2d,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c ny2d = first dimension of field arrays, must be >= 2*ny
c j2blok = number of data blocks
c kxp2 = number of data values per block
      complex cu, zero
      dimension cu(3,ny2d,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = at1*(dkx*aimag(cu(1,k,j,l)) + dky*aimag(cu(2,k,j,l)))
         at3 = aimag(cu(1,k,j,l)) - dkx*at2
         at4 = aimag(cu(2,k,j,l)) - dky*at2
         cu(1,k,j,l) = cmplx(0.,at3)
         cu(2,k,j,l) = cmplx(0.,at4)
         cu(1,k1,j,l) = cmplx(0.,-at3)
         cu(2,k1,j,l) = cmplx(0.,at4)
   10    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         cu(1,1,j,l) = zero
         cu(1,k1,j,l) = zero
         cu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         cu(2,k,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
   30    continue
         k1 = ny + 1
         cu(1,1,1,l) = zero
         cu(2,1,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPD2(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c j2blok = number of data blocks
c kxp2 = number of data values per block
      dimension cu(3,nyv,kxp2+1,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = at1*(dkx*cu(1,k,j,l) + dky*cu(2,k,j,l))
         at3 = cu(1,k,j,l) - dkx*at2
         at4 = cu(2,k,j,l) - dky*at2
         cu(1,k,j,l) = at3
         cu(2,k,j,l) = at4
   10    continue
c mode numbers ky = 0, ny
         cu(1,1,j,l) = 0.
         cu(1,ny+1,j,l) = 0.
         cu(2,ny+1,j,l) = 0.
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         cu(2,k,1,l) = 0.
   30    continue
         cu(1,1,1,l) = 0.
         cu(2,1,1,l) = 0.
      endif
      do 40 k = 1, ny1
      cu(1,k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPDX22(cu,nx,ny,kstrt,ny2d,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c ny2d = first dimension of field arrays, must be >= 2*ny
c j2blok = number of data blocks
c kxp2 = number of data values per block
      complex cu, zero
      dimension cu(2,ny2d,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = at1*(dkx*aimag(cu(1,k,j,l)) + dky*aimag(cu(2,k,j,l)))
         at3 = aimag(cu(1,k,j,l)) - dkx*at2
         at4 = aimag(cu(2,k,j,l)) - dky*at2
         cu(1,k,j,l) = cmplx(0.,at3)
         cu(2,k,j,l) = cmplx(0.,at4)
         cu(1,k1,j,l) = cmplx(0.,-at3)
         cu(2,k1,j,l) = cmplx(0.,at4)
   10    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         cu(1,1,j,l) = zero
         cu(1,k1,j,l) = zero
         cu(2,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         cu(2,k,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
   30    continue
         k1 = ny + 1
         cu(1,1,1,l) = zero
         cu(2,1,1,l) = zero
         cu(1,k1,1,l) = zero
         cu(2,k1,1,l) = zero
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPD22(cu,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,k,j,l) = i-th component of complex current density and
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
c j2blok = number of data blocks
c kxp2 = number of data values per block
      dimension cu(2,nyv,kxp2+1,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate transverse part of current
      if (kstrt.gt.nx) return
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = at1*(dkx*cu(1,k,j,l) + dky*cu(2,k,j,l))
         at3 = cu(1,k,j,l) - dkx*at2
         at4 = cu(2,k,j,l) - dky*at2
         cu(1,k,j,l) = at3
         cu(2,k,j,l) = at4
   10    continue
c mode numbers ky = 0, ny
         cu(1,1,j,l) = 0.
         cu(1,ny+1,j,l) = 0.
         cu(2,ny+1,j,l) = 0.
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         cu(2,k,1,l) = 0.
   30    continue
         cu(1,1,1,l) = 0.
         cu(2,1,1,l) = 0.
      endif
      do 40 k = 1, ny1
      cu(1,k,kxp2+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISDX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,kstr
     1t,ny2d,kxp2,j2blok,nyd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 21*nxc*nyc + 11*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 15*nxc*nyc + 8*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 8*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
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
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex cu, bxy, ffd, zero
      dimension cu(3,ny2d,kxp2,j2blok), bxy(3,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
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
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
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
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bxy(1,k,j,l) = cmplx(0.,at2*real(cu(3,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,-at3*real(cu(3,k,j,l)))
         bxy(3,k,j,l) = cmplx(at2*aimag(cu(1,k,j,l))-at3*aimag(cu(2,k,j,
     1l)),0.)
         bxy(1,k1,j,l) = bxy(1,k,j,l)
         bxy(2,k1,j,l) = -bxy(2,k,j,l)
         bxy(3,k1,j,l) = bxy(3,k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12 + real(cu(3,k,j,l))**2)
   50    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = zero
         bxy(3,1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + real(cu(3,1,j,l))**2)
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bxy(1,k,1,l) = zero
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = cmplx(at2*aimag(cu(1,k,1,l)),0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2 + real(cu(3,k,1,l))**2)
   70    continue
         k1 = ny + 1
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
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, ny
         k1 = ny2 - k
         at2 = ci2*real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at2*aimag(cu(1,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,at2*aimag(cu(2,k,j,l)))
         bxy(3,k,j,l) = cmplx(at2*real(cu(3,k,j,l)),0.)
         bxy(1,k1,j,l) = -bxy(1,k,j,l)
         bxy(2,k1,j,l) = bxy(2,k,j,l)
         bxy(3,k1,j,l) = -bxy(3,k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12 + real(cu(3,k,j,l))**2)
  110    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at2 = ci2*real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at2*aimag(cu(2,1,j,l)))
         bxy(3,1,j,l) = zero
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + real(cu(3,1,j,l))**2)
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, ny
         k1 = ny2 - k
         at2 = ci2*real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = cmplx(0.,at2*aimag(cu(1,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2 + real(cu(3,k,1,l))**2)
  130    continue
         k1 = ny + 1
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
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at1*aimag(cu(1,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,at1*aimag(cu(2,k,j,l)))
         bxy(3,k,j,l) = cmplx(at1*real(cu(3,k,j,l)),0.)
         bxy(1,k1,j,l) = -bxy(1,k,j,l)
         bxy(2,k1,j,l) = bxy(2,k,j,l)
         bxy(3,k1,j,l) = -bxy(3,k,j,l)
  170    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at1*aimag(cu(2,1,j,l)))
         bxy(3,1,j,l) = zero
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 190 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = cmplx(0.,at1*aimag(cu(1,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
  190    continue
         k1 = ny + 1
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
      subroutine PBPOISD23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,kstrt
     1,nyv,kxp2,j2blok,nyd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 20*nxc*nyc + 8*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 13*nxc*nyc + 8*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 3*nxc*nyc + 1*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
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
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension cu(3,nyv,kxp2+1,j2blok), bxy(3,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
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
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 110
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bxy(1,k,j,l) = at2*cu(3,k,j,l)
         bxy(2,k,j,l) = -at3*cu(3,k,j,l)
         bxy(3,k,j,l) = at3*cu(2,k,j,l) - at2*cu(1,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2 + cu(3,k,j,l
     1)**2)
   50    continue
c mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = 0.
         bxy(3,1,j,l) = at2*cu(2,1,j,l)
         wp = wp + at1*(cu(2,1,j,l)**2 + cu(3,1,j,l)**2)
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
      bxy(3,ny+1,j,l) = 0.
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bxy(1,k,1,l) = 0.
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = -at2*cu(1,k,1,l)
         wp = wp + at1*(cu(1,k,1,l)**2 + cu(3,k,1,l)**2)
   70    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
      do 80 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
      bxy(3,k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  110 if (isign.gt.1) go to 180
      wp = 0.0d0
      if (kstrt.gt.nx) go to 170
      do 160 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 130 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 120 k = 2, ny
         at2 = ci2*real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at2*cu(1,k,j,l)
         bxy(2,k,j,l) = at2*cu(2,k,j,l)
         bxy(3,k,j,l) = at2*cu(3,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2 + cu(3,k,j,l
     1)**2)
  120    continue
c mode numbers ky = 0, ny
         at2 = ci2*real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at2*cu(2,1,j,l)
         bxy(3,1,j,l) = 0.
         wp = wp + at1*(cu(2,1,j,l)**2 + cu(3,1,j,l)**2)
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
      bxy(3,ny+1,j,l) = 0.
  130 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 140 k = 2, ny
         at2 = ci2*real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = 0.
         wp = wp + at1*(cu(1,k,1,l)**2 + cu(3,k,1,l)**2)
  140    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
      do 150 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
      bxy(3,k,kxp2+1,l) = 0.
  150 continue
  160 continue
  170 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  180 if (kstrt.gt.nx) go to 240
      do 230 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 200 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 190 k = 2, ny
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
         bxy(3,k,j,l) = at1*cu(3,k,j,l)
  190    continue
c mode numbers ky = 0, ny
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
         bxy(3,1,j,l) = 0.
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
      bxy(3,ny+1,j,l) = 0.
  200 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 210 k = 2, ny
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = 0.
  210    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
      do 220 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
      bxy(3,k,kxp2+1,l) = 0.
  220 continue
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISDX22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny,k
     1strt,ny2d,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bz,wm
c approximate flop count is: 15*nxc*nyc + 8*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 11*nxc*nyc + 6*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 6*nxc*nyc
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex vector potential,
c bz(k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex cu, bxy, bz, ffd, zero
      dimension cu(2,ny2d,kxp2,j2blok), bxy(2,ny2d,kxp2,j2blok)
      dimension bz(ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
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
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
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
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bz(k,j,l) = cmplx(at2*aimag(cu(1,k,j,l))-at3*aimag(cu(2,k,j,l))
     1,0.)
         bz(k1,j,l) = bz(k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12)
   50    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bz(1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         bz(k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2)
      endif
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bz(k,1,l) = cmplx(at2*aimag(cu(1,k,1,l)),0.)
         bz(k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2)
   70    continue
         k1 = ny + 1
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
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 120 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 110 k = 2, ny
         k1 = ny2 - k
         at2 = ci2*real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at2*aimag(cu(1,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,at2*aimag(cu(2,k,j,l)))
         bxy(1,k1,j,l) = -bxy(1,k,j,l)
         bxy(2,k1,j,l) = bxy(2,k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12)
  110    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at2 = ci2*real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at2*aimag(cu(2,1,j,l)))
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2)
      endif
  120 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 130 k = 2, ny
         k1 = ny2 - k
         at2 = ci2*real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = cmplx(0.,at2*aimag(cu(1,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2)
  130    continue
         k1 = ny + 1
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
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 180 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 170 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at1*aimag(cu(1,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,at1*aimag(cu(2,k,j,l)))
         bxy(1,k1,j,l) = -bxy(1,k,j,l)
         bxy(2,k1,j,l) = bxy(2,k,j,l)
  170    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at1*aimag(cu(2,1,j,l)))
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
      endif
  180 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 190 k = 2, ny
         k1 = ny2 - k
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = cmplx(0.,at1*aimag(cu(1,k,1,l)))
         bxy(2,k,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
  190    continue
         k1 = ny + 1
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
      subroutine PBPOISD22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny,ks
     1trt,nyv,kxp2,j2blok,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,kstrt,ny2d,kxp,j2blok,nyd
c output: ffd
c for isign = -1, input:cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bz,wm
c approximate flop count is: 15*nxc*nyc + 7*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy,wm
c approximate flop count is: 10*nxc*nyc + 6*(nxc + nyc)
c for isign = 2, input: cu,ffd,isign,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd
c output: bxy
c approximate flop count is: 2*nxc*nyc + 1*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c if isign = 0, form factor array is prepared
c if isign < 0, magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c             s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex vector potential,
c bz(k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension cu(2,nyv,kxp2+1,j2blok), bxy(2,nyv,kxp2+1,j2blok)
      dimension bz(nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
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
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffd(k,j,l) = cmplx(affp,1.)
      else
         ffd(k,j,l) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
   30 continue
      return
   40 if (isign.gt.0) go to 110
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 100
      do 90 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 60 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 50 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))*aimag(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         bz(k,j,l) = at3*cu(2,k,j,l) - at2*cu(1,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2)
   50    continue
c mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j,l))*aimag(ffd(1,j,l))
         at2 = dkx*at1
         bz(1,j,l) = at2*cu(2,1,j,l)
         wp = wp + at1*cu(2,1,j,l)**2
      endif
      bz(ny+1,j,l) = 0.
   60 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 70 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))*aimag(ffd(k,1,l))
         at2 = dky*at1
         bz(k,1,l) = -at2*cu(1,k,1,l)
         wp = wp + at1*cu(1,k,1,l)**2
   70    continue
         bz(1,1,l) = 0.
      endif
      do 80 k = 1, ny1
      bz(k,kxp2+1,l) = 0.
   80 continue
   90 continue
  100 continue
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
  110 if (isign.gt.1) go to 180
      wp = 0.0d0
      if (kstrt.gt.nx) go to 170
      do 160 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 130 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 120 k = 2, ny
         at2 = ci2*real(ffd(k,j,l))
         at1 = at2*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at2*cu(1,k,j,l)
         bxy(2,k,j,l) = at2*cu(2,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2)
  120    continue
c mode numbers ky = 0, ny
         at2 = ci2*real(ffd(1,j,l))
         at1 = at2*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at2*cu(2,1,j,l)
         wp = wp + at1*cu(2,1,j,l)**2
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
  130 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 140 k = 2, ny
         at2 = ci2*real(ffd(k,1,l))
         at1 = at2*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at2*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
         wp = wp + at1*cu(1,k,1,l)**2
  140    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
      endif
      do 150 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
  150 continue
  160 continue
  170 continue
      wm = float(nx*ny)*wp
      return
c calculate smoothing
  180 if (kstrt.gt.nx) go to 240
      do 230 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 200 j = 1, kxp2
      if ((j+joff).gt.0) then
         do 190 k = 2, ny
         at1 = aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at1*cu(1,k,j,l)
         bxy(2,k,j,l) = at1*cu(2,k,j,l)
  190    continue
c mode numbers ky = 0, ny
         at1 = aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at1*cu(2,1,j,l)
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
  200 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 210 k = 2, ny
         at1 = aimag(ffd(k,1,l))
         bxy(1,k,1,l) = at1*cu(1,k,1,l)
         bxy(2,k,1,l) = 0.
  210    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
      endif
      do 220 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
  220 continue
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine IPBPOISDX23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,ny2d,kxp2,j2bl
     1ok,nyd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c input: cu,ffd,ci,nx,ny,kstrt,ny2d,kxp2,j2blok,nyd, output: bxy,wm
c approximate flop count is: 21*nxc*nyc + 11*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*,
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex cu, bxy, ffd, zero
      dimension cu(3,ny2d,kxp2,j2blok), bxy(3,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = cmplx(0.,at2*real(cu(3,k,j,l)))
         bxy(2,k,j,l) = cmplx(0.,-at3*real(cu(3,k,j,l)))
         bxy(3,k,j,l) = cmplx(at2*aimag(cu(1,k,j,l))-at3*aimag(cu(2,k,j,
     1l)),0.)
         bxy(1,k1,j,l) = bxy(1,k,j,l)
         bxy(2,k1,j,l) = -bxy(2,k,j,l)
         bxy(3,k1,j,l) = bxy(3,k,j,l)
         wp = wp + 2.0*at1*(aimag(cu(1,k,j,l))**2 + aimag(cu(2,k,j,l))**
     12 + real(cu(3,k,j,l))**2)
   10    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         at1 = ci2*real(ffd(1,j,l))
         at2 = dkx*at1
         at1 = at1*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = zero
         bxy(3,1,j,l) = cmplx(-at2*aimag(cu(2,1,j,l)),0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         wp = wp + at1*(aimag(cu(2,1,j,l))**2 + real(cu(3,1,j,l))**2)
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))
         at2 = dky*at1
         at1 = at1*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = zero
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = cmplx(at2*aimag(cu(1,k,1,l)),0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         wp = wp + at1*(aimag(cu(1,k,1,l))**2 + real(cu(3,k,1,l))**2)
   30    continue
         k1 = ny + 1
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
      subroutine IPBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,kstrt,nyv,kxp2,j2blok
     1,nyd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field with dirichlet boundary conditions (zero potential),
c for distributed data.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c input: cu,ffd,ci,nx,ny,kstrt,nyv,kxp2,j2blok,nyd, output: bxy,wm
c approximate flop count is: 20*nxc*nyc + 9*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c cu(i,k,j,l) = i-th component of complex current density and
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
c real(ffd(k,j,l)) = potential green's function g
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp
      complex ffd
      dimension cu(3,nyv,kxp2+1,j2blok), bxy(3,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.nx) go to 60
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,j,l))
         at2 = dky*at1
         at3 = dkx*at1
         at1 = at1*aimag(ffd(k,j,l))
         bxy(1,k,j,l) = at2*cu(3,k,j,l)
         bxy(2,k,j,l) = -at3*cu(3,k,j,l)
         bxy(3,k,j,l) = at3*cu(2,k,j,l) - at2*cu(1,k,j,l)
         wp = wp + 2.0*at1*(cu(1,k,j,l)**2 + cu(2,k,j,l)**2 + cu(3,k,j,l
     1)**2)
   10    continue
c mode numbers ky = 0, ny
         at1 = ci2*real(ffd(1,j,l))
         at2 = dkx*at1
         at1 = at1*aimag(ffd(1,j,l))
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = 0.
         bxy(3,1,j,l) = at2*cu(2,1,j,l)
         wp = wp + at1*(cu(2,1,j,l)**2 + cu(3,1,j,l)**2)
      endif
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
      bxy(3,ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         at1 = ci2*real(ffd(k,1,l))
         at2 = dky*at1
         at1 = at1*aimag(ffd(k,1,l))
         bxy(1,k,1,l) = 0.
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = -at2*cu(1,k,1,l)
         wp = wp + at1*(cu(1,k,1,l)**2 + cu(3,k,1,l)**2)
   30    continue
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
      do 40 k = 1, ny1
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
      bxy(3,k,kxp2+1,l) = 0.
   40 continue
   50 continue
   60 continue
      wm = float(nx*ny)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWELDX2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstrt,
     1ny2d,kxp2,j2blok,nyd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with dirichlet boundary
c conditions (zero potential).
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 62*nxc*nyc + 36*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
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
c where kx = pi*j/nx, ky = pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers
c and similarly for bx, by, bz.
c cu(i,k,j,l) = i-th component of complex current density and
c exy(i,k,j,l) = i-th component of complex electric field,
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
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
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      double precision wp, ws
      complex exy, bxy, cu, ffd
      complex zero
      dimension exy(3,ny2d,kxp2,j2blok), bxy(3,ny2d,kxp2,j2blok)
      dimension cu(3,ny2d,kxp2,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      if (ci.le.0.) return
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
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
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffd(k,j,l))
         at7 = aimag(exy(1,k,j,l))
         at8 = aimag(exy(2,k,j,l))
         at9 = real(exy(3,k,j,l))
c update magnetic field half time step, ky > 0
         at4 = aimag(bxy(1,k,j,l)) - dth*(dky*at9)
         at5 = aimag(bxy(2,k,j,l)) + dth*(dkx*at9)
         at6 = real(bxy(3,k,j,l)) + dth*(dkx*at8 - dky*at7)
c update electric field whole time step
         at7 = at7 + cdt*(dky*at6) - afdt*aimag(cu(1,k,j,l))
         at8 = at8 - cdt*(dkx*at6) - afdt*aimag(cu(2,k,j,l))
         at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*real(cu(3,k,j,l))
c update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 + dth*(dkx*at8 - dky*at7)
         ws = ws + 2.0*anorm*(at7*at7 + at8*at8 + at9*at9)
         wp = wp + 2.0*anorm*(at4*at4 + at5*at5 + at6*at6)
         exy(1,k,j,l) = cmplx(0.,at7)
         exy(2,k,j,l) = cmplx(0.,at8)
         exy(3,k,j,l) = cmplx(at9,0.)
         bxy(1,k,j,l) = cmplx(0.,at4)
         bxy(2,k,j,l) = cmplx(0.,at5)
         bxy(3,k,j,l) = cmplx(at6,0.)
c update electric and magnetic fields, ky < 0
         exy(1,k1,j,l) = cmplx(0.,-at7)
         exy(2,k1,j,l) = cmplx(0.,at8)
         exy(3,k1,j,l) = cmplx(-at9,0.)
         bxy(1,k1,j,l) = cmplx(0.,at4)
         bxy(2,k1,j,l) = cmplx(0.,-at5)
         bxy(3,k1,j,l) = cmplx(at6,0.)
   10    continue
c mode numbers ky = 0, ny
         k1 = ny + 1
         afdt = adt*aimag(ffd(1,j,l))
         at8 = aimag(exy(2,1,j,l))
         at9 = real(exy(3,1,j,l))
c update magnetic field half time step, ky > 0
         at5 = aimag(bxy(2,1,j,l)) + dth*(dkx*at9)
         at6 = real(bxy(3,1,j,l)) + dth*(dkx*at8)
c update electric field whole time step
         at8 = at8 - cdt*(dkx*at6) - afdt*aimag(cu(2,1,j,l))
         at9 = at9 - cdt*(dkx*at5) - afdt*real(cu(3,1,j,l))
c update magnetic field half time step and store electric field
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 + dth*(dkx*at8)
         ws = ws + anorm*(at8*at8 + at9*at9)
         wp = wp + anorm*(at5*at5 + at6*at6)
         exy(1,1,j,l) = zero
         exy(2,1,j,l) = cmplx(0.,at8)
         exy(3,1,j,l) = cmplx(at9,0.)
         bxy(1,1,j,l) = zero
         bxy(2,1,j,l) = cmplx(0.,at5)
         bxy(3,1,j,l) = cmplx(at6,0.)
         bxy(1,k1,j,l) = zero
         bxy(2,k1,j,l) = zero
         bxy(3,k1,j,l) = zero
         exy(1,k1,j,l) = zero
         exy(2,k1,j,l) = zero
         exy(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffd(k,1,l))
         at7 = aimag(exy(1,k,1,l))
         at9 = real(exy(3,k,1,l))
c update magnetic field half time step, ky > 0
         at4 = aimag(bxy(1,k,1,l)) - dth*(dky*at9)
         at6 = real(bxy(3,k,1,l)) - dth*(dky*at7)
c update electric field whole time step
         at7 = at7 + cdt*(dky*at6) - afdt*aimag(cu(1,k,1,l))
         at9 = at9 + cdt*(dky*at4) - afdt*real(cu(3,k,1,l))
c update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at6 = at6 - dth*(dky*at7)
         ws = ws + anorm*(at7*at7 + at9*at9)
         wp = wp + anorm*(at4*at4 + at6*at6)
         exy(1,k,1,l) = cmplx(0.,at7)
         exy(2,k,1,l) = zero
         exy(3,k,1,l) = cmplx(at9,0.)
         bxy(1,k,1,l) = cmplx(0.,at4)
         bxy(2,k,1,l) = zero
         bxy(3,k,1,l) = cmplx(at6,0.)
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
   30    continue
         k1 = ny + 1
         bxy(1,1,1,l) = zero
         bxy(2,1,1,l) = zero
         bxy(3,1,1,l) = zero
         exy(1,1,1,l) = zero
         exy(2,1,1,l) = zero
         exy(3,1,1,l) = zero
         bxy(1,k1,1,l) = zero
         bxy(2,k1,1,l) = zero
         bxy(3,k1,1,l) = zero
         exy(1,k1,1,l) = zero
         exy(2,k1,1,l) = zero
         exy(3,k1,1,l) = zero
      endif
   40 continue
   50 continue
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWELD2(exy,bxy,cu,ffd,affp,ci,dt,wf,wm,nx,ny,kstrt,n
     1yv,kxp2,j2blok,nyd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with dirichlet boundary
c conditions (zero potential).
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 62*nxc*nyc + 36*(nxc + nyc)
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
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
c where kx = pi*j/nx, ky = pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers
c and similarly for bx, by, bz.
c cu(i,k,j,l) = i-th component of complex current density and
c exy(i,k,j,l) = i-th component of complex electric field,
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c aimag(ffd(k,j,l)) = finite-size particle shape factor s
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
c nyv = first dimension of field arrays, must be >= ny+1
c nyd = first dimension of form factor array, must be >= ny
      double precision wp, ws
      complex ffd
      dimension exy(3,nyv,kxp2+1,j2blok), bxy(3,nyv,kxp2+2,j2blok)
      dimension cu(3,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      if (ci.le.0.) return
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
      if (kstrt.gt.nx) go to 60
c calculate the electromagnetic fields
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffd(k,j,l))
         at7 = exy(1,k,j,l)
         at8 = exy(2,k,j,l)
         at9 = exy(3,k,j,l)
c update magnetic field half time step, ky > 0
         at4 = bxy(1,k,j,l) - dth*(dky*at9)
         at5 = bxy(2,k,j,l) + dth*(dkx*at9)
         at6 = bxy(3,k,j,l) + dth*(dky*at7 - dkx*at8)
c update electric field whole time step
         at7 = at7 - cdt*(dky*at6) - afdt*cu(1,k,j,l)
         at8 = at8 + cdt*(dkx*at6) - afdt*cu(2,k,j,l)
         at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*cu(3,k,j,l)
c update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 + dth*(dky*at7 - dkx*at8)
         ws = ws + 2.0*anorm*(at7*at7 + at8*at8 + at9*at9)
         wp = wp + 2.0*anorm*(at4*at4 + at5*at5 + at6*at6)
         exy(1,k,j,l) = at7
         exy(2,k,j,l) = at8
         exy(3,k,j,l) = at9
         bxy(1,k,j,l) = at4
         bxy(2,k,j,l) = at5
         bxy(3,k,j,l) = at6
   10    continue
c mode numbers ky = 0, ny
         afdt = adt*aimag(ffd(1,j,l))
         at8 = exy(2,1,j,l)
         at9 = exy(3,1,j,l)
c update magnetic field half time step, ky > 0
         at5 = bxy(2,1,j,l) + dth*(dkx*at9)
         at6 = bxy(3,1,j,l) - dth*(dkx*at8)
c update electric field whole time step
         at8 = at8 + cdt*(dkx*at6) - afdt*cu(2,1,j,l)
         at9 = at9 - cdt*(dkx*at5) - afdt*cu(3,1,j,l)
c update magnetic field half time step and store electric field
         at5 = at5 + dth*(dkx*at9)
         at6 = at6 - dth*(dkx*at8)
         ws = ws + anorm*(at8*at8 + at9*at9)
         wp = wp + anorm*(at5*at5 + at6*at6)
         exy(1,1,j,l) = 0.
         exy(2,1,j,l) = at8
         exy(3,1,j,l) = at9
         bxy(1,1,j,l) = 0.
         bxy(2,1,j,l) = at5
         bxy(3,1,j,l) = at6
      endif
      exy(1,ny+1,j,l) = 0.
      exy(2,ny+1,j,l) = 0.
      exy(3,ny+1,j,l) = 0.
      bxy(1,ny+1,j,l) = 0.
      bxy(2,ny+1,j,l) = 0.
      bxy(3,ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         afdt = adt*aimag(ffd(k,1,l))
         at7 = exy(1,k,1,l)
         at9 = exy(3,k,1,l)
c update magnetic field half time step, ky > 0
         at4 = bxy(1,k,1,l) - dth*(dky*at9)
         at6 = bxy(3,k,1,l) + dth*(dky*at7)
c update electric field whole time step
         at7 = at7 - cdt*(dky*at6) - afdt*cu(1,k,1,l)
         at9 = at9 + cdt*(dky*at4) - afdt*cu(3,k,1,l)
c update magnetic field half time step and store electric field
         at4 = at4 - dth*(dky*at9)
         at6 = at6 + dth*(dky*at7)
         ws = ws + anorm*(at7*at7 + at9*at9)
         wp = wp + anorm*(at4*at4 + at6*at6)
         exy(1,k,1,l) = at7
         exy(2,k,1,l) = 0.
         exy(3,k,1,l) = at9
         bxy(1,k,1,l) = at4
         bxy(2,k,1,l) = 0.
         bxy(3,k,1,l) = at6
   30    continue
         exy(1,1,1,l) = 0.
         exy(2,1,1,l) = 0.
         exy(3,1,1,l) = 0.
         bxy(1,1,1,l) = 0.
         bxy(2,1,1,l) = 0.
         bxy(3,1,1,l) = 0.
      endif
      do 40 k = 1, ny1
      exy(1,k,kxp2+1,l) = 0.
      exy(2,k,kxp2+1,l) = 0.
      exy(3,k,kxp2+1,l) = 0.
      bxy(1,k,kxp2+1,l) = 0.
      bxy(2,k,kxp2+1,l) = 0.
      bxy(3,k,kxp2+1,l) = 0.
   40 continue
   50 continue
   60 continue
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PDMFIELDD2(q2,q,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c this subroutine copies the charge density into a smaller array
c which would have been created by a fast sine transform in x and y
      implicit none
      integer nx, ny, kstrt, nyv, ny2d, kxp2, j2blok
      complex q2
      real q
      dimension q2(ny2d,kxp2,j2blok), q(nyv,kxp2+1,j2blok)
      integer j, k, l, ny1
      if (kstrt.gt.nx) return
      ny1 = ny + 1
      do 40 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      do 20 j = 1, kxp2
      do 10 k = 1, ny1
      q(k,j,l) = -real(q2(k,j,l))
   10 continue
   20 continue
c mode number kx = nx/2 is zero.
      do 30 k = 1, ny1
      q(k,kxp2+1,l) = 0.
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFMFIELDD2(q2,q,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c this subroutine copies the charge density into a smaller array
c which would have been created by a fast cosine transform in x and y
      implicit none
      integer nx, ny, kstrt, nyv, ny2d, kxp2, j2blok
      complex q2
      real q
      dimension q2(ny2d,kxp2,j2blok), q(nyv,kxp2+1,j2blok)
      integer j, k, l, ny1
      if (kstrt.gt.nx) return
      ny1 = ny + 1
      do 40 l = 1, j2blok
c mode numbers kx > 0 and 0 < ky < ny
      do 20 j = 1, kxp2
      do 10 k = 1, ny1
      q(k,j,l) = real(q2(k,j,l))
   10 continue
   20 continue
c mode number kx = nx/2 is assumed to be zero.
c if not zero, it is actually on the wrong processor
      do 30 k = 1, ny1
      q(k,kxp2+1,l) = 0.
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCMFIELDD2(cu2,cu,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c this subroutine copies the current into a smaller array
      implicit none
      integer nx, ny, kstrt, nyv, ny2d, kxp2, j2blok
      complex cu2
      real cu
      dimension cu2(3,ny2d,kxp2,j2blok), cu(3,nyv,kxp2+1,j2blok)
      integer j, k, l, ny1
      if (kstrt.gt.nx) return
      ny1 = ny + 1
      do 40 l = 1, j2blok
      do 20 j = 1, kxp2
      do 10 k = 1, ny1
      cu(1,k,j,l) = -aimag(cu2(1,k,j,l))
      cu(2,k,j,l) = -aimag(cu2(2,k,j,l))
      cu(3,k,j,l) = -real(cu2(3,k,j,l))
   10 continue
   20 continue
c mode number kx = nx/2 is assumed to be zero for cu(1)
c if not zero, it is actually on the wrong processor
      do 30 k = 1, ny1
      cu(1,k,kxp2+1,l) = 0.
      cu(2,k,kxp2+1,l) = 0.
      cu(3,k,kxp2+1,l) = 0.
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PEMFIELDD2(fxy,exy,ffd,isign,nx,ny,kstrt,ny2d,nyv,kxp2,
     1j2blok,nyd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign <= 0
c adds image charges appropriate for electric field if isign >= 0
c or appropriate for magnetic field if isign < 0
c includes additional smoothing for isign /= 0
      implicit none
      integer isign, nx, ny, kstrt, ny2d, nyv, kxp2, j2blok, nyd
      complex fxy, ffd
      real exy
      dimension fxy(3,ny2d,kxp2,j2blok), exy(3,nyv,kxp2+1,j2blok)
      dimension ffd(nyd,kxp2,j2blok)
      integer j, k, l, ny2, ks, joff, k1
      real at1
      complex zero
      ny2 = 2*ny + 2
      ks = kstrt - 2
      zero = cmplx(0.,0.)
      if (kstrt.gt.nx) return
c add the fields
      if (isign.gt.0) then
         do 40 l = 1, j2blok
         do 20 j = 1, kxp2
         joff = kxp2*(l + ks) - 1
         if ((j+joff).gt.0) then
            do 10 k = 2, ny
            k1 = ny2 - k
            at1 = aimag(ffd(k,j,l))
            fxy(1,k,j,l) = fxy(1,k,j,l) + cmplx(0.,-exy(1,k,j,l)*at1)
            fxy(2,k,j,l) = fxy(2,k,j,l) + cmplx(0.,-exy(2,k,j,l)*at1)
            fxy(3,k,j,l) = fxy(3,k,j,l) + cmplx(-exy(3,k,j,l)*at1,0.)
            fxy(1,k1,j,l) = fxy(1,k1,j,l) + cmplx(0.,exy(1,k,j,l)*at1)
            fxy(2,k1,j,l) = fxy(2,k1,j,l) + cmplx(0.,-exy(2,k,j,l)*at1)
            fxy(3,k1,j,l) = fxy(3,k1,j,l) + cmplx(exy(3,k,j,l)*at1,0.)
   10       continue
            at1 = aimag(ffd(1,j,l))
            fxy(1,1,j,l) = fxy(1,1,j,l) + cmplx(0.,-exy(1,1,j,l)*at1)
            fxy(2,1,j,l) = fxy(2,1,j,l) + cmplx(0.,-exy(2,1,j,l)*at1)
            fxy(3,1,j,l) = fxy(3,1,j,l) + cmplx(-exy(3,1,j,l)*at1,0.)
         endif
   20    continue
         if ((l+ks).eq.0) then
            do 30 k = 1, ny
            at1 = aimag(ffd(k,1,l))
            fxy(1,k,1,l) = fxy(1,k,1,l) + cmplx(0.,-exy(1,k,1,l)*at1)
            fxy(2,k,1,l) = fxy(2,k,1,l) + cmplx(0.,-exy(2,k,1,l)*at1)
            fxy(3,k,1,l) = fxy(3,k,1,l) + cmplx(-exy(3,k,1,l)*at1,0.)
   30       continue
         endif
   40    continue
c copy the magnetic fields
      else if (isign.lt.0) then
         do 80 l = 1, j2blok
         do 60 j = 1, kxp2
         joff = kxp2*(l + ks) - 1
         if ((j+joff).gt.0) then
            do 50 k = 2, ny
            k1 = ny2 - k
            at1 = aimag(ffd(k,j,l))
            fxy(1,k,j,l) = cmplx(0.,-exy(1,k,j,l)*at1)
            fxy(2,k,j,l) = cmplx(0.,-exy(2,k,j,l)*at1)
            fxy(3,k,j,l) = cmplx(exy(3,k,j,l)*at1,0.)
            fxy(1,k1,j,l) = cmplx(0.,-exy(1,k,j,l)*at1)
            fxy(2,k1,j,l) = cmplx(0.,exy(2,k,j,l)*at1)
            fxy(3,k1,j,l) = cmplx(exy(3,k,j,l)*at1,0.)
   50       continue
            k1 = ny + 1
            at1 = aimag(ffd(1,j,l))
            fxy(1,1,j,l) = cmplx(0.,-exy(1,1,j,l)*at1)
            fxy(2,1,j,l) = cmplx(0.,-exy(2,1,j,l)*at1)
            fxy(3,1,j,l) = cmplx(exy(3,1,j,l)*at1,0.)
            fxy(1,k1,j,l) = zero
            fxy(2,k1,j,l) = zero
            fxy(3,k1,j,l) = zero
         endif
   60    continue
         if ((l+ks).eq.0) then
            do 70 k = 2, ny
            k1 = ny2 - k
            at1 = aimag(ffd(k,1,l))
            fxy(1,k,1,l) = cmplx(0.,-exy(1,k,1,l)*at1)
            fxy(2,k,1,l) = cmplx(0.,-exy(2,k,1,l)*at1)
            fxy(3,k,1,l) = cmplx(exy(3,k,1,l)*at1,0.)
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
   70       continue
            k1 = ny + 1
            at1 = aimag(ffd(1,1,l))
            fxy(1,1,1,l) = cmplx(-exy(1,1,1,l)*at1,0.)
            fxy(2,1,1,l) = cmplx(-exy(2,1,1,l)*at1,0.)
            fxy(3,1,1,l) = cmplx(exy(3,1,1,l)*at1,0.)
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
         endif
   80    continue
c copy the electric fields
      else
         do 120 l = 1, j2blok
         do 100 j = 1, kxp2
         joff = kxp2*(l + ks) - 1
         if ((j+joff).gt.0) then
            do 90 k = 2, ny
            k1 = ny2 - k
            fxy(1,k,j,l) = cmplx(0.,-exy(1,k,j,l))
            fxy(2,k,j,l) = cmplx(0.,-exy(2,k,j,l))
            fxy(3,k,j,l) = cmplx(-exy(3,k,j,l),0.)
            fxy(1,k1,j,l) = cmplx(0.,exy(1,k,j,l))
            fxy(2,k1,j,l) = cmplx(0.,-exy(2,k,j,l))
            fxy(3,k1,j,l) = cmplx(exy(3,k,j,l),0.)
   90       continue
            k1 = ny + 1
            fxy(1,1,j,l) = cmplx(0.,-exy(1,1,j,l))
            fxy(2,1,j,l) = cmplx(0.,-exy(2,1,j,l))
            fxy(3,1,j,l) = cmplx(-exy(3,1,j,l),0.)
            fxy(1,k1,j,l) = cmplx(0.,-exy(1,k1,j,l))
            fxy(2,k1,j,l) = cmplx(0.,-exy(2,k1,j,l))
            fxy(3,k1,j,l) = cmplx(-exy(3,k1,j,l),0.)
         endif
  100    continue
         if ((l+ks).eq.0) then
            do 110 k = 2, ny
            k1 = ny2 - k
            fxy(1,k,1,l) = cmplx(0.,-exy(1,k,1,l))
            fxy(2,k,1,l) = cmplx(0.,-exy(2,k,1,l))
            fxy(3,k,1,l) = cmplx(-exy(3,k,1,l),0.)
c mode number kx = nx/2 is assumed to be zero for fxy(1)
c if not zero, it is actually on the wrong processor
            fxy(1,k1,1,l) = zero
            fxy(2,k1,1,l) = zero
            fxy(3,k1,1,l) = zero
  110       continue
            k1 = ny + 1
            fxy(1,1,1,l) = cmplx(-exy(1,1,1,l),0.)
            fxy(2,1,1,l) = cmplx(-exy(2,1,1,l),0.)
            fxy(3,1,1,l) = cmplx(-exy(3,1,1,l),0.)
            fxy(1,k1,1,l) = cmplx(-exy(1,k1,1,l),0.)
            fxy(2,k1,1,l) = cmplx(-exy(2,k1,1,l),0.)
            fxy(3,k1,1,l) = cmplx(-exy(3,k1,1,l),0.)
         endif
  120    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMFIELDD2(pot2,pot,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c copies image charges appropriate for potential
      implicit none
      integer nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
      complex pot2
      real pot
      dimension pot2(ny2d,kxp2,j2blok), pot(nyv,kxp2+1,j2blok)
      integer j, k, l, ny2, ks, joff, k1
      complex zero
      ny2 = 2*ny + 2
      ks = kstrt - 2
      zero = cmplx(0.,0.)
      if (kstrt.gt.nx) return
      do 40 l = 1, j2blok
      do 20 j = 1, kxp2
      joff = kxp2*(l + ks) - 1
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         pot2(k,j,l) = cmplx(-pot(k,j,l),0.)
         pot2(k1,j,l) = cmplx(pot(k,j,l),0.)
   10    continue
         k1 = ny + 1
         pot2(1,j,l) = cmplx(-pot(1,j,l),0.)
         pot2(k1,j,l) = cmplx(-pot(k1,j,l),0.)
      endif
   20 continue
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         pot2(k,1,l) = cmplx(-pot(k,1,l),0.)
         pot2(k1,1,l) = zero
   30    continue
         k1 = ny + 1
         pot2(1,1,l) = cmplx(-pot(1,1,l),0.)
         pot2(k1,1,l) = cmplx(-pot(k1,1,l),0.)
         endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBMFIELDD2(fxy,bxy,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c copies image charges appropriate for magnetic field
      implicit none
      integer nx, ny, kstrt, ny2d, nyv, kxp2, j2blok
      complex fxy
      real bxy
      dimension fxy(3,ny2d,kxp2,j2blok), bxy(3,nyv,kxp2+1,j2blok)
      integer j, k, l, ny2, ks, joff, k1
      complex zero
      ny2 = 2*ny + 2
      ks = kstrt - 2
      zero = cmplx(0.,0.)
      if (kstrt.gt.nx) return
c copy the magnetic fields
         do 40 l = 1, j2blok
         do 20 j = 1, kxp2
         joff = kxp2*(l + ks) - 1
         if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         fxy(1,k,j,l) = cmplx(0.,-bxy(1,k,j,l))
         fxy(2,k,j,l) = cmplx(0.,-bxy(2,k,j,l))
         fxy(3,k,j,l) = cmplx(bxy(3,k,j,l),0.)
         fxy(1,k1,j,l) = cmplx(0.,-bxy(1,k,j,l))
         fxy(2,k1,j,l) = cmplx(0.,bxy(2,k,j,l))
         fxy(3,k1,j,l) = cmplx(bxy(3,k,j,l),0.)
   10    continue
         k1 = ny + 1
         fxy(1,1,j,l) = cmplx(0.,-bxy(1,1,j,l))
         fxy(2,1,j,l) = cmplx(0.,-bxy(2,1,j,l))
         fxy(3,1,j,l) = cmplx(bxy(3,1,j,l),0.)
         fxy(1,k1,j,l) = cmplx(0.,-bxy(1,k1,j,l))
         fxy(2,k1,j,l) = cmplx(0.,-bxy(2,k1,j,l))
         fxy(3,k1,j,l) = cmplx(bxy(3,k1,j,l),0.)
      endif
   20 continue
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         fxy(1,k,1,l) = cmplx(0.,-bxy(1,k,1,l))
         fxy(2,k,1,l) = cmplx(0.,-bxy(2,k,1,l))
         fxy(3,k,1,l) = cmplx(bxy(3,k,1,l),0.)
c mode number kx = nx/2 is assumed to be zero for fxy(2,3)
c if not zero, it is actually on the wrong processor
         fxy(1,k1,1,l) = zero
         fxy(2,k1,1,l) = zero
         fxy(3,k1,1,l) = zero
   30    continue
         k1 = ny + 1
         fxy(1,1,1,l) = cmplx(-bxy(1,1,1,l),0.)
         fxy(2,1,1,l) = cmplx(-bxy(2,1,1,l),0.)
         fxy(3,1,1,l) = cmplx(bxy(3,1,1,l),0.)
         fxy(1,k1,1,l) = cmplx(-bxy(1,k1,1,l),0.)
         fxy(2,k1,1,l) = cmplx(-bxy(2,k1,1,l),0.)
         fxy(3,k1,1,l) = cmplx(bxy(3,k1,1,l),0.)
      endif
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCPFIELDD2(fxy,exy,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok)
c this subroutine copies complex vector fields and adds image charges
c appropriate for electric field
      implicit none
      integer nx, ny, kstrt, ny2d, nyv, kxp2, j2blok, nyd
      complex fxy
      real exy
      dimension fxy(3,ny2d,kxp2,j2blok), exy(3,nyv,kxp2+1,j2blok)
c local data
      integer isign
      complex ffd
      dimension ffd(1,1,1)
      isign = 0
      call PEMFIELDD2(fxy,exy,ffd,isign,nx,ny,kstrt,ny2d,nyv,kxp2,j2blok
     1,nyd)
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVPOTDX23(bxy,axy,nx,ny,kstrt,ny2d,kxp2,j2blok)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with dirichlet boundary conditions (zero potential).
c for distributed data.
c input: bxy,nx,ny,kstrt,ny2d,kxp2,j2blok, output: axy
c approximate flop count is: 14*nxc*nyc + 4*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c axy(i,k,j,l) = i-th component of complex vector potential,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c ny2d = first dimension of field arrays, must be >= 2*ny
      complex bxy, axy, zero
      dimension bxy(3,ny2d,kxp2,j2blok), axy(3,ny2d,kxp2,j2blok)
      ny2 = 2*ny + 2
      ks = kstrt - 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c calculate vector potential
      if (kstrt.gt.nx) go to 50
      do 40 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = dky*at1
         at3 = dkx*at1
         at4 = real(bxy(3,k,j,l))
         at5 = aimag(bxy(2,k,j,l))
         at6 = aimag(bxy(1,k,j,l))
         at6 = at2*at6 - at3*at5
         at5 = -at3*at4
         at4 = at2*at4
         axy(1,k,j,l) = cmplx(0.,at4)
         axy(2,k,j,l) = cmplx(0.,at5)
         axy(3,k,j,l) = cmplx(at6,0.)
         axy(1,k1,j,l) = cmplx(0.,-at4)
         axy(2,k1,j,l) = cmplx(0.,at5)
         axy(3,k1,j,l) = cmplx(-at6,0.)
   10    continue
c mode numbers ky = 0, ny/2
         k1 = ny + 1
         at2 = 1.0/dkx
         at4 = real(bxy(3,1,j,l))
         at5 = aimag(bxy(2,1,j,l))
         axy(1,1,j,l) = zero
         axy(2,1,j,l) = cmplx(0.,-at2*at4)
         axy(3,1,j,l) = cmplx(-at2*at5,0.)
         axy(1,k1,j,l) = zero
         axy(2,k1,j,l) = zero
         axy(3,k1,j,l) = zero
      endif
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         k1 = ny2 - k
         dky = dny*float(k - 1)
         at2 = 1.0/dky
         at4 = real(bxy(3,k,1,l))
         at6 = aimag(bxy(1,k,1,l))
         axy(1,k,1,l) = cmplx(0.,at2*at4)
         axy(2,k,1,l) = zero
         axy(3,k,1,l) = cmplx(at2*at6,0.)
         axy(1,k1,1,l) = zero
         axy(2,k1,1,l) = zero
         axy(3,k1,1,l) = zero
   30    continue
         k1 = ny + 1
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
      subroutine PAVPOTD23(bxy,axy,nx,ny,kstrt,nyv,kxp2,j2blok)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with dirichlet boundary conditions (zero potential).
c for distributed data.
c input: bxy,nx,ny,kstrt,ny2d,kxp2,j2blok, output: axy
c approximate flop count is: 12*nxc*nyc + 4*(nxc + nyc)
c and nxc*nyc divides
c where nxc = (nx/2-1)/nvp, nyc = ny/2 - 1, and nvp = number of procs
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c bxy(i,k,j,l) = i-th component of complex magnetic field,
c axy(i,k,j,l) = i-th component of complex vector potential,
c for fourier mode (jj-1,k-1), where jj = j + kxp*(l - 1)
c nx/ny = system length in x/y direction
c j2blok = number of data blocks
c kxp2 = number of data values per block
c kstrt = starting data block number
c nyv = first dimension of field arrays, must be >= ny+1
      dimension bxy(3,nyv,kxp2+1,j2blok), axy(3,nyv,kxp2+1,j2blok)
      ks = kstrt - 2
      ny1 = ny + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate vector potential
      if (kstrt.gt.nx) go to 60
      do 50 l = 1, j2blok
c mode numbers 0 < kx < nx and 0 < ky < ny
      joff = kxp2*(l + ks) - 1
      do 20 j = 1, kxp2
      dkx = dnx*float(j + joff)
      dkx2 = dkx*dkx
      if ((j+joff).gt.0) then
         do 10 k = 2, ny
         dky = dny*float(k - 1)
         at1 = 1./(dky*dky + dkx2)
         at2 = dky*at1
         at3 = dkx*at1
         at4 = bxy(3,k,j,l)
         at5 = bxy(2,k,j,l)
         at6 = bxy(1,k,j,l)
         axy(1,k,j,l) = -at2*at4
         axy(2,k,j,l) = at3*at4
         axy(3,k,j,l) = at2*at6 - at3*at5
   10    continue
c mode numbers ky = 0, ny/2
         at2 = 1.0/dkx
         at4 = bxy(3,1,j,l)
         at5 = bxy(2,1,j,l)
         axy(1,1,j,l) = 0.
         axy(2,1,j,l) = at2*at4
         axy(3,1,j,l) = -at2*at5
      endif
      axy(1,ny+1,j,l) = 0.
      axy(2,ny+1,j,l) = 0.
      axy(3,ny+1,j,l) = 0.
   20 continue
c mode numbers kx = 0, nx/2
      if ((l+ks).eq.0) then
         do 30 k = 2, ny
         dky = dny*float(k - 1)
         at2 = 1.0/dky
         at4 = bxy(3,k,1,l)
         at6 = bxy(1,k,1,l)
         axy(1,k,1,l) = -at2*at4
         axy(2,k,1,l) = 0.
         axy(3,k,1,l) = at2*at6
   30    continue
         axy(1,1,1,l) = 0.
         axy(2,1,1,l) = 0.
         axy(3,1,1,l) = 0.
      endif
      do 40 k = 1, ny1
      axy(1,k,kxp2+1,l) = 0.
      axy(2,k,kxp2+1,l) = 0.
      axy(3,k,kxp2+1,l) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
