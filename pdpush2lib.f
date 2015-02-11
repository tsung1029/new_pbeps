c-----------------------------------------------------------------------
c 2d parallel PIC library for pushing particles with darwin electric and
c magnetic fields and depositing current and derivative of current
c pdpush2lib.f contains procedures to process particles with darwin
c              electric and magnetic fields:
c PGMJPOST2 deposits momentum flux for 2-1/2d code, quadratic
c           interpolation, STANDARD optimization, and distributed data.
c PGSMJPOST2 deposits momentum flux for 2-1/2d code, quadratic
c            interpolation, LOOKAHEAD optimization, and distributed data
c PGDCJPOST2 deposits momentum flux, acceleration density and current
c            density for 2-1/2d code, quadratic interpolation, STANDARD
c            optimization, and distributed data.
c PGSDCJPOST2 deposits momentum flux, acceleration density and current
c             density for 2-1/2d code, quadratic interpolation,
c             LOOKAHEAD optimization, and distributed data.
c PGMJPOST22 deposits momentum flux for 2d code, quadratic
c            interpolation, STANDARD optimization, and distributed data.
c PGSMJPOST22 deposits momentum flux for 2d code, quadratic
c             interpolation, LOOKAHEAD optimization, and distributed
c             data.
c PGDCJPOST22 deposits momentum flux, acceleration density and current
c             density for 2d code, quadratic interpolation, STANDARD
c             optimization, and distributed data.
c PGSDCJPOST22 deposits momentum flux, acceleration density and current
c              density for 2d code, quadratic interpolation, LOOKAHEAD
c              optimization, and distributed data.
c PGMJPOST2L deposits momentum flux for 2-1/2d code, linear
c            interpolation, STANDARD optimization, and distributed data.
c PGSMJPOST2L deposits momentum flux for 2-1/2d code, linear
c             interpolation, LOOKAHEAD optimization, and distributed
c             data.
c PGDCJPOST2L deposits momentum flux, acceleration density and current
c             density for 2-1/2d code, linear interpolation, STANDARD
c             optimization, and distributed data.
c PGSDCJPOST2L deposits momentum flux, acceleration density and current
c              density for 2-1/2d code, linear interpolation, LOOKAHEAD
c              optimization, and distributed data.
c PGMJPOST22L deposits momentum flux for 2d code, linear interpolation,
c             STANDARD optimization, and distributed data.
c PGSMJPOST22LL deposits momentum flux for 2d code, linear interpolation,
c              LOOKAHEAD optimization, and distributed data.
c PGDCJPOST22L deposits momentum flux, acceleration density and current
c              density for 2d code, linear interpolation, STANDARD
c              optimization, and distributed data.
c PGSDCJPOST22L deposits momentum flux, acceleration density and current
c               density for 2d code, linear interpolation, LOOKAHEAD
c               optimization, and distributed data.
c written by viktor k. decyk, ucla
c copyright 2006, regents of the university of california
c update: august 29, 2009
c-----------------------------------------------------------------------
      subroutine PGMJPOST2(part,amu,npp,noff,qm,idimp,npmax,nblok,nxv,ny
     1pmx)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation, and distributed data.
c scalar version using guard cells, for distributed data
c 112 flops/particle, 41 loads, 36 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c part(5,n,l) = z velocity of particle n at t - dt/2 in partition l
c amu(i,j+1,k,l) = ith component of momentum flux at grid point (j,kk)
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of flux array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, amu, qm
      dimension part(idimp,npmax,nblok), amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, nl, np, ml, mp
      real qmh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy, vz, v1, v2, v3, v4
      qmh = .5*qm
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      amx = qm*(.75 - dxp*dxp)
      ml = mm - mnoff
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = qmh*(.5 - dxp)**2
      np = nl + 2
      dxp = qmh*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      amu(1,nl,mm,l) = amu(1,nl,mm,l) + v1*dx
      amu(2,nl,mm,l) = amu(2,nl,mm,l) + v2*dx
      amu(3,nl,mm,l) = amu(3,nl,mm,l) + v3*dx
      amu(4,nl,mm,l) = amu(4,nl,mm,l) + v4*dx
      dx = dxl*dyl
      amu(1,nn,mm,l) = amu(1,nn,mm,l) + v1*dy
      amu(2,nn,mm,l) = amu(2,nn,mm,l) + v2*dy
      amu(3,nn,mm,l) = amu(3,nn,mm,l) + v3*dy
      amu(4,nn,mm,l) = amu(4,nn,mm,l) + v4*dy
      dy = amx*dyl
      amu(1,np,mm,l) = amu(1,np,mm,l) + v1*dz
      amu(2,np,mm,l) = amu(2,np,mm,l) + v2*dz
      amu(3,np,mm,l) = amu(3,np,mm,l) + v3*dz
      amu(4,np,mm,l) = amu(4,np,mm,l) + v4*dz
      dz = dxp*dyl
      amu(1,nl,ml,l) = amu(1,nl,ml,l) + v1*dx
      amu(2,nl,ml,l) = amu(2,nl,ml,l) + v2*dx
      amu(3,nl,ml,l) = amu(3,nl,ml,l) + v3*dx
      amu(4,nl,ml,l) = amu(4,nl,ml,l) + v4*dx
      dx = dxl*dyp
      amu(1,nn,ml,l) = amu(1,nn,ml,l) + v1*dy
      amu(2,nn,ml,l) = amu(2,nn,ml,l) + v2*dy
      amu(3,nn,ml,l) = amu(3,nn,ml,l) + v3*dy
      amu(4,nn,ml,l) = amu(4,nn,ml,l) + v4*dy
      dy = amx*dyp
      amu(1,np,ml,l) = amu(1,np,ml,l) + v1*dz
      amu(2,np,ml,l) = amu(2,np,ml,l) + v2*dz
      amu(3,np,ml,l) = amu(3,np,ml,l) + v3*dz
      amu(4,np,ml,l) = amu(4,np,ml,l) + v4*dz
      dz = dxp*dyp
      amu(1,nl,mp,l) = amu(1,nl,mp,l) + v1*dx
      amu(2,nl,mp,l) = amu(2,nl,mp,l) + v2*dx
      amu(3,nl,mp,l) = amu(3,nl,mp,l) + v3*dx
      amu(4,nl,mp,l) = amu(4,nl,mp,l) + v4*dx
      amu(1,nn,mp,l) = amu(1,nn,mp,l) + v1*dy
      amu(2,nn,mp,l) = amu(2,nn,mp,l) + v2*dy
      amu(3,nn,mp,l) = amu(3,nn,mp,l) + v3*dy
      amu(4,nn,mp,l) = amu(4,nn,mp,l) + v4*dy
      amu(1,np,mp,l) = amu(1,np,mp,l) + v1*dz
      amu(2,np,mp,l) = amu(2,np,mp,l) + v2*dz
      amu(3,np,mp,l) = amu(3,np,mp,l) + v3*dz
      amu(4,np,mp,l) = amu(4,np,mp,l) + v4*dz
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSMJPOST2(part,amu,npp,noff,qm,idimp,npmax,nblok,nxv,n
     1xyp)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation, and distributed data.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 112 flops/particle, 41 loads, 36 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c part(5,n,l) = z velocity of particle n at t - dt/2 in partition l
c amu(i,n,l) = ith component of momentum flux at grid point (j,kk),
c where n = j + 1 + nxv*k, and kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyp = first actual dimension of current array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nxyp
      real part, amu, qm
      dimension part(idimp,npmax,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      integer nop, nnn, mmn, j, l, mnoff, nn, mm, ml, mn, mp
      real dxn, dyn, qmh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy, vz, v1, v2, v3, v4
      real dx1, dy1, dx2, dy2, dx3
      qmh = .5*qm
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l) + .5
      mmn = part(2,1,l) + .5
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp(l)
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j,l) + .5
      mmn = part(2,j,l) + .5
      dxp = dxn
      dyp = dyn
      dxn = part(1,j,l) - float(nnn)
      dyn = part(2,j,l) - float(mmn)
      ml = mm + nn
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      mn = ml + nxv
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mmn = mmn - mnoff
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j-1,l)
      vy = part(4,j-1,l)
      vz = part(5,j-1,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      dx1 = amu(1,mn,l) + v1*dx
      dy1 = amu(2,mn,l) + v2*dx
      amy = amu(3,mn,l) + v3*dx
      vx = amu(4,mn,l) + v4*dx
      dx2 = amu(1,mn+1,l) + v1*dy
      dy2 = amu(2,mn+1,l) + v2*dy
      dx3 = amu(3,mn+1,l) + v3*dy
      vy = amu(4,mn+1,l) + v4*dy
      dx = amu(1,mn+2,l) + v1*dz
      dy = amu(2,mn+2,l) + v2*dz
      vz = amu(3,mn+2,l) + v3*dz
      dz = amu(4,mn+2,l) + v4*dz
      amu(1,mn,l) = dx1
      amu(2,mn,l) = dy1
      amu(3,mn,l) = amy
      amu(4,mn,l) = vx
      amu(1,mn+1,l) = dx2
      amu(2,mn+1,l) = dy2
      amu(3,mn+1,l) = dx3
      amu(4,mn+1,l) = vy
      amu(1,mn+2,l) = dx
      amu(2,mn+2,l) = dy
      amu(3,mn+2,l) = vz
      amu(4,mn+2,l) = dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = amu(1,ml,l) + v1*dx
      dy1 = amu(2,ml,l) + v2*dx
      amy = amu(3,ml,l) + v3*dx
      vx = amu(4,ml,l) + v4*dx
      dx2 = amu(1,ml+1,l) + v1*dy
      dy2 = amu(2,ml+1,l) + v2*dy
      dyl = amu(3,ml+1,l) + v3*dy
      vy = amu(4,ml+1,l) + v4*dy
      dx = amu(1,ml+2,l) + v1*dz
      dy = amu(2,ml+2,l) + v2*dz
      vz = amu(3,ml+2,l) + v3*dz
      dz = amu(4,ml+2,l) + v4*dz
      amu(1,ml,l) = dx1
      amu(2,ml,l) = dy1
      amu(3,ml,l) = amy
      amu(4,ml,l) = vx
      amu(1,ml+1,l) = dx2
      amu(2,ml+1,l) = dy2
      amu(3,ml+1,l) = dyl
      amu(4,ml+1,l) = vy
      amu(1,ml+2,l) = dx
      amu(2,ml+2,l) = dy
      amu(3,ml+2,l) = vz
      amu(4,ml+2,l) = dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = amu(1,mp,l) + v1*dx
      dy1 = amu(2,mp,l) + v2*dx
      amy = amu(3,mp,l) + v3*dx
      vx = amu(4,mp,l) + v4*dx
      dxl = amu(1,mp+1,l) + v1*dy
      amx = amu(2,mp+1,l) + v2*dy
      dxp = amu(3,mp+1,l) + v3*dy
      vy = amu(4,mp+1,l) + v4*dy
      dx = amu(1,mp+2,l) + v1*dz
      dy = amu(2,mp+2,l) + v2*dz
      vz = amu(3,mp+2,l) + v3*dz
      dz = amu(4,mp+2,l) + v4*dz
      amu(1,mp,l) = dx1
      amu(2,mp,l) = dy1
      amu(3,mp,l) = amy
      amu(4,mp,l) = vx
      amu(1,mp+1,l) = dxl
      amu(2,mp+1,l) = amx
      amu(3,mp+1,l) = dxp
      amu(4,mp+1,l) = vy
      amu(1,mp+2,l) = dx
      amu(2,mp+2,l) = dy
      amu(3,mp+2,l) = vz
      amu(4,mp+2,l) = dz
   10 continue
      nop = npp(l)
c deposit momentum flux for last particle
      nn = nnn + 1
      mm = nxv*mmn
      ml = mm + nn
      amx = qm*(.75 - dxn*dxn)
      amy = .75 - dyn*dyn
      mn = ml + nxv
      dxl = qmh*(.5 - dxn)**2
      dxp = qmh*(.5 + dxn)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      vz = part(5,nop,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      amu(1,mn,l) = amu(1,mn,l) + v1*dx
      amu(2,mn,l) = amu(2,mn,l) + v2*dx
      amu(3,mn,l) = amu(3,mn,l) + v3*dx
      amu(4,mn,l) = amu(4,mn,l) + v4*dx
      amu(1,mn+1,l) = amu(1,mn+1,l) + v1*dy
      amu(2,mn+1,l) = amu(2,mn+1,l) + v2*dy
      amu(3,mn+1,l) = amu(3,mn+1,l) + v3*dy
      amu(4,mn+1,l) = amu(4,mn+1,l) + v4*dy
      amu(1,mn+2,l) = amu(1,mn+2,l) + v1*dz
      amu(2,mn+2,l) = amu(2,mn+2,l) + v2*dz
      amu(3,mn+2,l) = amu(3,mn+2,l) + v3*dz
      amu(4,mn+2,l) = amu(4,mn+2,l) + v4*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      amu(1,ml,l) = amu(1,ml,l) + v1*dx
      amu(2,ml,l) = amu(2,ml,l) + v2*dx
      amu(3,ml,l) = amu(3,ml,l) + v3*dx
      amu(4,ml,l) = amu(4,ml,l) + v4*dx
      amu(1,ml+1,l) = amu(1,ml+1,l) + v1*dy
      amu(2,ml+1,l) = amu(2,ml+1,l) + v2*dy
      amu(3,ml+1,l) = amu(3,ml+1,l) + v3*dy
      amu(4,ml+1,l) = amu(4,ml+1,l) + v4*dy
      amu(1,ml+2,l) = amu(1,ml+2,l) + v1*dz
      amu(2,ml+2,l) = amu(2,ml+2,l) + v2*dz
      amu(3,ml+2,l) = amu(3,ml+2,l) + v3*dz
      amu(4,ml+2,l) = amu(4,ml+2,l) + v4*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      amu(1,mp,l) = amu(1,mp,l) + v1*dx
      amu(2,mp,l) = amu(2,mp,l) + v2*dx
      amu(3,mp,l) = amu(3,mp,l) + v3*dx
      amu(4,mp,l) = amu(4,mp,l) + v4*dx
      amu(1,mp+1,l) = amu(1,mp+1,l) + v1*dy
      amu(2,mp+1,l) = amu(2,mp+1,l) + v2*dy
      amu(3,mp+1,l) = amu(3,mp+1,l) + v3*dy
      amu(4,mp+1,l) = amu(4,mp+1,l) + v4*dy
      amu(1,mp+2,l) = amu(1,mp+2,l) + v1*dz
      amu(2,mp+2,l) = amu(2,mp+2,l) + v2*dz
      amu(3,mp+2,l) = amu(3,mp+2,l) + v3*dz
      amu(4,mp+2,l) = amu(4,mp+2,l) + v4*dz
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,dt,i
     1dimp,npmax,nblok,nxv,nypmx)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation.
c scalar version using guard cells, for distributed data
c 407 flops/particle, 1 divide, 150 loads, 80 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c dcu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c dcu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c dcu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c dcu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c dcu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c dcu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c dcu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c dcu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c part(5,n,l) = z velocity of particle n at t - dt/2 in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c fxy(3,j+1,k,l) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j+1,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j+1,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j+1,k,l) = ith component of current density
c at grid point j,kk for i = 1, 3
c dcu(i,j+1,k,l) = ith component of acceleration density
c at grid point j,kk for i = 1, 3
c amu(i,j+1,k,l) = ith component of momentum flux
c at grid point j,kk for i = 1, 4
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of flux array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension cu(3,nxv,nypmx,nblok), dcu(3,nxv,nypmx,nblok)
      dimension amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, nl, np, ml, mp
      real qtmh, dti, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      ml = mm - mnoff
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = .5*(.5 - dxp)**2
      np = nl + 2
      dxp = .5*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))
      dz = amy*(dxl*fxy(3,nl,mm,l) + amx*fxy(3,nn,mm,l) + dxp*fxy(3,np,m
     1m,l)) + dyl*(dxl*fxy(3,nl,ml,l) + amx*fxy(3,nn,ml,l) + dxp*fxy(3,n
     2p,ml,l)) + dyp*(dxl*fxy(3,nl,mp,l) + amx*fxy(3,nn,mp,l) + dxp*fxy(
     33,np,mp,l))
c find magnetic field
      ox = amy*(dxl*bxy(1,nl,mm,l) + amx*bxy(1,nn,mm,l) + dxp*bxy(1,np,m
     1m,l)) + dyl*(dxl*bxy(1,nl,ml,l) + amx*bxy(1,nn,ml,l) + dxp*bxy(1,n
     2p,ml,l)) + dyp*(dxl*bxy(1,nl,mp,l) + amx*bxy(1,nn,mp,l) + dxp*bxy(
     31,np,mp,l))
      oy = amy*(dxl*bxy(2,nl,mm,l) + amx*bxy(2,nn,mm,l) + dxp*bxy(2,np,m
     1m,l)) + dyl*(dxl*bxy(2,nl,ml,l) + amx*bxy(2,nn,ml,l) + dxp*bxy(2,n
     2p,ml,l)) + dyp*(dxl*bxy(2,nl,mp,l) + amx*bxy(2,nn,mp,l) + dxp*bxy(
     32,np,mp,l))
      oz = amy*(dxl*bxy(3,nl,mm,l) + amx*bxy(3,nn,mm,l) + dxp*bxy(3,np,m
     1m,l)) + dyl*(dxl*bxy(3,nl,ml,l) + amx*bxy(3,nn,ml,l) + dxp*bxy(3,n
     2p,ml,l)) + dyp*(dxl*bxy(3,nl,mp,l) + amx*bxy(3,nn,mp,l) + dxp*bxy(
     33,np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      amu(1,nl,mm,l) = amu(1,nl,mm,l) + v1*dx
      amu(2,nl,mm,l) = amu(2,nl,mm,l) + v2*dx
      amu(3,nl,mm,l) = amu(3,nl,mm,l) + v3*dx
      amu(4,nl,mm,l) = amu(4,nl,mm,l) + v4*dx
      dcu(1,nl,mm,l) = dcu(1,nl,mm,l) + vx*dx
      dcu(2,nl,mm,l) = dcu(2,nl,mm,l) + vy*dx
      dcu(3,nl,mm,l) = dcu(3,nl,mm,l) + vz*dx
      cu(1,nl,mm,l) = cu(1,nl,mm,l) + ox*dx
      cu(2,nl,mm,l) = cu(2,nl,mm,l) + oy*dx
      cu(3,nl,mm,l) = cu(3,nl,mm,l) + oz*dx
      dx = dxl*dyl
      amu(1,nn,mm,l) = amu(1,nn,mm,l) + v1*dy
      amu(2,nn,mm,l) = amu(2,nn,mm,l) + v2*dy
      amu(3,nn,mm,l) = amu(3,nn,mm,l) + v3*dy
      amu(4,nn,mm,l) = amu(4,nn,mm,l) + v4*dy
      dcu(1,nn,mm,l) = dcu(1,nn,mm,l) + vx*dy
      dcu(2,nn,mm,l) = dcu(2,nn,mm,l) + vy*dy
      dcu(3,nn,mm,l) = dcu(3,nn,mm,l) + vz*dy
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + ox*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + oy*dy
      cu(3,nn,mm,l) = cu(3,nn,mm,l) + oz*dy
      dy = amx*dyl
      amu(1,np,mm,l) = amu(1,np,mm,l) + v1*dz
      amu(2,np,mm,l) = amu(2,np,mm,l) + v2*dz
      amu(3,np,mm,l) = amu(3,np,mm,l) + v3*dz
      amu(4,np,mm,l) = amu(4,np,mm,l) + v4*dz
      dcu(1,np,mm,l) = dcu(1,np,mm,l) + vx*dz
      dcu(2,np,mm,l) = dcu(2,np,mm,l) + vy*dz
      dcu(3,np,mm,l) = dcu(3,np,mm,l) + vz*dz
      cu(1,np,mm,l) = cu(1,np,mm,l) + ox*dz
      cu(2,np,mm,l) = cu(2,np,mm,l) + oy*dz
      cu(3,np,mm,l) = cu(3,np,mm,l) + oz*dz
      dz = dxp*dyl
      amu(1,nl,ml,l) = amu(1,nl,ml,l) + v1*dx
      amu(2,nl,ml,l) = amu(2,nl,ml,l) + v2*dx
      amu(3,nl,ml,l) = amu(3,nl,ml,l) + v3*dx
      amu(4,nl,ml,l) = amu(4,nl,ml,l) + v4*dx
      dcu(1,nl,ml,l) = dcu(1,nl,ml,l) + vx*dx
      dcu(2,nl,ml,l) = dcu(2,nl,ml,l) + vy*dx
      dcu(3,nl,ml,l) = dcu(3,nl,ml,l) + vz*dx
      cu(1,nl,ml,l) = cu(1,nl,ml,l) + ox*dx
      cu(2,nl,ml,l) = cu(2,nl,ml,l) + oy*dx
      cu(3,nl,ml,l) = cu(3,nl,ml,l) + oz*dx
      dx = dxl*dyp
      amu(1,nn,ml,l) = amu(1,nn,ml,l) + v1*dy
      amu(2,nn,ml,l) = amu(2,nn,ml,l) + v2*dy
      amu(3,nn,ml,l) = amu(3,nn,ml,l) + v3*dy
      amu(4,nn,ml,l) = amu(4,nn,ml,l) + v4*dy
      dcu(1,nn,ml,l) = dcu(1,nn,ml,l) + vx*dy
      dcu(2,nn,ml,l) = dcu(2,nn,ml,l) + vy*dy
      dcu(3,nn,ml,l) = dcu(3,nn,ml,l) + vz*dy
      cu(1,nn,ml,l) = cu(1,nn,ml,l) + ox*dy
      cu(2,nn,ml,l) = cu(2,nn,ml,l) + oy*dy
      cu(3,nn,ml,l) = cu(3,nn,ml,l) + oz*dy
      dy = amx*dyp
      amu(1,np,ml,l) = amu(1,np,ml,l) + v1*dz
      amu(2,np,ml,l) = amu(2,np,ml,l) + v2*dz
      amu(3,np,ml,l) = amu(3,np,ml,l) + v3*dz
      amu(4,np,ml,l) = amu(4,np,ml,l) + v4*dz
      dcu(1,np,ml,l) = dcu(1,np,ml,l) + vx*dz
      dcu(2,np,ml,l) = dcu(2,np,ml,l) + vy*dz
      dcu(3,np,ml,l) = dcu(3,np,ml,l) + vz*dz
      cu(1,np,ml,l) = cu(1,np,ml,l) + ox*dz
      cu(2,np,ml,l) = cu(2,np,ml,l) + oy*dz
      cu(3,np,ml,l) = cu(3,np,ml,l) + oz*dz
      dz = dxp*dyp
      amu(1,nl,mp,l) = amu(1,nl,mp,l) + v1*dx
      amu(2,nl,mp,l) = amu(2,nl,mp,l) + v2*dx
      amu(3,nl,mp,l) = amu(3,nl,mp,l) + v3*dx
      amu(4,nl,mp,l) = amu(4,nl,mp,l) + v4*dx
      dcu(1,nl,mp,l) = dcu(1,nl,mp,l) + vx*dx
      dcu(2,nl,mp,l) = dcu(2,nl,mp,l) + vy*dx
      dcu(3,nl,mp,l) = dcu(3,nl,mp,l) + vz*dx
      cu(1,nl,mp,l) = cu(1,nl,mp,l) + ox*dx
      cu(2,nl,mp,l) = cu(2,nl,mp,l) + oy*dx
      cu(3,nl,mp,l) = cu(3,nl,mp,l) + oz*dx
      amu(1,nn,mp,l) = amu(1,nn,mp,l) + v1*dy
      amu(2,nn,mp,l) = amu(2,nn,mp,l) + v2*dy
      amu(3,nn,mp,l) = amu(3,nn,mp,l) + v3*dy
      amu(4,nn,mp,l) = amu(4,nn,mp,l) + v4*dy
      dcu(1,nn,mp,l) = dcu(1,nn,mp,l) + vx*dy
      dcu(2,nn,mp,l) = dcu(2,nn,mp,l) + vy*dy
      dcu(3,nn,mp,l) = dcu(3,nn,mp,l) + vz*dy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + ox*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + oy*dy
      cu(3,nn,mp,l) = cu(3,nn,mp,l) + oz*dy
      amu(1,np,mp,l) = amu(1,np,mp,l) + v1*dz
      amu(2,np,mp,l) = amu(2,np,mp,l) + v2*dz
      amu(3,np,mp,l) = amu(3,np,mp,l) + v3*dz
      amu(4,np,mp,l) = amu(4,np,mp,l) + v4*dz
      dcu(1,np,mp,l) = dcu(1,np,mp,l) + vx*dz
      dcu(2,np,mp,l) = dcu(2,np,mp,l) + vy*dz
      dcu(3,np,mp,l) = dcu(3,np,mp,l) + vz*dz
      cu(1,np,mp,l) = cu(1,np,mp,l) + ox*dz
      cu(2,np,mp,l) = cu(2,np,mp,l) + oy*dz
      cu(3,np,mp,l) = cu(3,np,mp,l) + oz*dz
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSDCJPOST2(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,dt,
     1idimp,npmax,nblok,nxv,nxyp)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 407 flops/particle, 1 divide, 150 loads, 80 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c dcu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c dcu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c dcu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c dcu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c dcu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c dcu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c dcu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c dcu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = velocity vx of particle n at t - dt/2 in partition l
c part(4,n,l) = velocity vy of particle n at t - dt/2 in partition l
c part(5,n,l) = velocity vz of particle n at t - dt/2 in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c fxy(3,j+1,k,l) = z component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bxy(1,j+1,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j+1,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,n,l) = ith component of current density at grid point j,kk
c where n = j + 1 + nxv*k, for i = 1, 3
c dcu(i,n,l) = ith component of acceleration density at grid point j,kk
c where n = j + 1 + nxv*k, for i = 1, 3
c amu(i,n,l) = ith component of momentum flux at grid point j,kk
c where n = j + 1 + nxv*k, for i = 1, 4
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+3
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nxyp
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension cu(3,nxyp,nblok), dcu(3,nxyp,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      integer nop1, j, l, mnoff, nop, nnn, mmn, nn, mm, ml, mn, mp
      real qtmh, dti, dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4
      real dx1, dy1, dx2, dy2, dx3, dy3, dx4, dy4, dx5, dy5, dx6
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l) + .5
      mmn = part(2,1,l) + .5
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l) + .5
      mmn = part(2,j+1,l) + .5
      dx = dxn
      dy = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      ml = mm + nn
      amx = .75 - dx*dx
      amy = .75 - dy*dy
      mn = ml + nxv
      dxl = .5*(.5 - dx)**2
      dxp = .5*(.5 + dx)**2
      mp = mn + nxv
      dyl = .5*(.5 - dy)**2
      dyp = .5*(.5 + dy)**2
      mmn = mmn - mnoff
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l))
      dz = amy*(dxl*fxy(3,mn,l) + amx*fxy(3,mn+1,l) + dxp*fxy(3,mn+2,l))
     1 + dyl*(dxl*fxy(3,ml,l) + amx*fxy(3,ml+1,l) + dxp*fxy(3,ml+2,l)) +
     2 dyp*(dxl*fxy(3,mp,l) + amx*fxy(3,mp+1,l) + dxp*fxy(3,mp+2,l))
c find magnetic field
      ox = amy*(dxl*bxy(1,mn,l) + amx*bxy(1,mn+1,l) + dxp*bxy(1,mn+2,l))
     1 + dyl*(dxl*bxy(1,ml,l) + amx*bxy(1,ml+1,l) + dxp*bxy(1,ml+2,l)) +
     2 dyp*(dxl*bxy(1,mp,l) + amx*bxy(1,mp+1,l) + dxp*bxy(1,mp+2,l)) 
      oy = amy*(dxl*bxy(2,mn,l) + amx*bxy(2,mn+1,l) + dxp*bxy(2,mn+2,l))
     1 + dyl*(dxl*bxy(2,ml,l) + amx*bxy(2,ml+1,l) + dxp*bxy(2,ml+2,l)) +
     2 dyp*(dxl*bxy(2,mp,l) + amx*bxy(2,mp+1,l) + dxp*bxy(2,mp+2,l)) 
      oz = amy*(dxl*bxy(3,mn,l) + amx*bxy(3,mn+1,l) + dxp*bxy(3,mn+2,l))
     1 + dyl*(dxl*bxy(3,ml,l) + amx*bxy(3,ml+1,l) + dxp*bxy(3,ml+2,l)) +
     2 dyp*(dxl*bxy(3,mp,l) + amx*bxy(3,mp+1,l) + dxp*bxy(3,mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      dx1 = amu(1,mn,l) + v1*dx
      dy1 = amu(2,mn,l) + v2*dx
      amy = amu(3,mn,l) + v3*dx
      dx2 = amu(4,mn,l) + v4*dx
      dy2 = amu(1,mn+1,l) + v1*dy
      dx3 = amu(2,mn+1,l) + v2*dy
      dy3 = amu(3,mn+1,l) + v3*dy
      dx4 = amu(4,mn+1,l) + v4*dy
      dy4 = amu(1,mn+2,l) + v1*dz
      dx5 = amu(2,mn+2,l) + v2*dz
      dy5 = amu(3,mn+2,l) + v3*dz
      dx6 = amu(4,mn+2,l) + v4*dz
      amu(1,mn,l) = dx1
      amu(2,mn,l) = dy1
      amu(3,mn,l) = amy
      amu(4,mn,l) = dx2
      amu(1,mn+1,l) = dy2
      amu(2,mn+1,l) = dx3
      amu(3,mn+1,l) = dy3
      amu(4,mn+1,l) = dx4
      amu(1,mn+2,l) = dy4
      amu(2,mn+2,l) = dx5
      amu(3,mn+2,l) = dy5
      amu(4,mn+2,l) = dx6
      dx1 = dcu(1,mn,l) + vx*dx
      dy1 = dcu(2,mn,l) + vy*dx
      amy = dcu(3,mn,l) + vz*dx
      dx2 = dcu(1,mn+1,l) + vx*dy
      dy2 = dcu(2,mn+1,l) + vy*dy
      dx3 = dcu(3,mn+1,l) + vz*dy
      dy3 = dcu(1,mn+2,l) + vx*dz
      dx4 = dcu(2,mn+2,l) + vy*dz
      dy4 = dcu(3,mn+2,l) + vz*dz
      dcu(1,mn,l) = dx1
      dcu(2,mn,l) = dy1
      dcu(3,mn,l) = amy
      dcu(1,mn+1,l) = dx2
      dcu(2,mn+1,l) = dy2
      dcu(3,mn+1,l) = dx3
      dcu(1,mn+2,l) = dy3
      dcu(2,mn+2,l) = dx4
      dcu(3,mn+2,l) = dy4
      dx1 = cu(1,mn,l) + ox*dx
      dy1 = cu(2,mn,l) + oy*dx
      amy = cu(3,mn,l) + oz*dx
      dx2 = cu(1,mn+1,l) + ox*dy
      dy2 = cu(2,mn+1,l) + oy*dy
      dx3 = cu(3,mn+1,l) + oz*dy
      dy3 = cu(1,mn+2,l) + ox*dz
      dx4 = cu(2,mn+2,l) + oy*dz
      dy4 = cu(3,mn+2,l) + oz*dz
      cu(1,mn,l) = dx1
      cu(2,mn,l) = dy1
      cu(3,mn,l) = amy
      cu(1,mn+1,l) = dx2
      cu(2,mn+1,l) = dy2
      cu(3,mn+1,l) = dx3
      cu(1,mn+2,l) = dy3
      cu(2,mn+2,l) = dx4
      cu(3,mn+2,l) = dy4
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = amu(1,ml,l) + v1*dx
      dy1 = amu(2,ml,l) + v2*dx
      amy = amu(3,ml,l) + v3*dx
      dx2 = amu(4,ml,l) + v4*dx
      dy2 = amu(1,ml+1,l) + v1*dy
      dx3 = amu(2,ml+1,l) + v2*dy
      dy3 = amu(3,ml+1,l) + v3*dy
      dx4 = amu(4,ml+1,l) + v4*dy
      dy4 = amu(1,ml+2,l) + v1*dz
      dx5 = amu(2,ml+2,l) + v2*dz
      dy5 = amu(3,ml+2,l) + v3*dz
      dx6 = amu(4,ml+2,l) + v4*dz
      amu(1,ml,l) = dx1
      amu(2,ml,l) = dy1
      amu(3,ml,l) = amy
      amu(4,ml,l) = dx2
      amu(1,ml+1,l) = dy2
      amu(2,ml+1,l) = dx3
      amu(3,ml+1,l) = dy3
      amu(4,ml+1,l) = dx4
      amu(1,ml+2,l) = dy4
      amu(2,ml+2,l) = dx5
      amu(3,ml+2,l) = dy5
      amu(4,ml+2,l) = dx6
      dx1 = dcu(1,ml,l) + vx*dx
      dy1 = dcu(2,ml,l) + vy*dx
      amy = dcu(3,ml,l) + vz*dx
      dx2 = dcu(1,ml+1,l) + vx*dy
      dy2 = dcu(2,ml+1,l) + vy*dy
      dx3 = dcu(3,ml+1,l) + vz*dy
      dy3 = dcu(1,ml+2,l) + vx*dz
      dx4 = dcu(2,ml+2,l) + vy*dz
      dy4 = dcu(3,ml+2,l) + vz*dz
      dcu(1,ml,l) = dx1
      dcu(2,ml,l) = dy1
      dcu(3,ml,l) = amy
      dcu(1,ml+1,l) = dx2
      dcu(2,ml+1,l) = dy2
      dcu(3,ml+1,l) = dx3
      dcu(1,ml+2,l) = dy3
      dcu(2,ml+2,l) = dx4
      dcu(3,ml+2,l) = dy4
      dx1 = cu(1,ml,l) + ox*dx
      dy1 = cu(2,ml,l) + oy*dx
      amy = cu(3,ml,l) + oz*dx
      dx2 = cu(1,ml+1,l) + ox*dy
      dy2 = cu(2,ml+1,l) + oy*dy
      dx3 = cu(3,ml+1,l) + oz*dy
      dy3 = cu(1,ml+2,l) + ox*dz
      dx4 = cu(2,ml+2,l) + oy*dz
      dy4 = cu(3,ml+2,l) + oz*dz
      cu(1,ml,l) = dx1
      cu(2,ml,l) = dy1
      cu(3,ml,l) = amy
      cu(1,ml+1,l) = dx2
      cu(2,ml+1,l) = dy2
      cu(3,ml+1,l) = dx3
      cu(1,ml+2,l) = dy3
      cu(2,ml+2,l) = dx4
      cu(3,ml+2,l) = dy4
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = amu(1,mp,l) + v1*dx
      dy1 = amu(2,mp,l) + v2*dx
      amy = amu(3,mp,l) + v3*dx
      dx2 = amu(4,mp,l) + v4*dx
      dy2 = amu(1,mp+1,l) + v1*dy
      dx3 = amu(2,mp+1,l) + v2*dy
      dy3 = amu(3,mp+1,l) + v3*dy
      dx4 = amu(4,mp+1,l) + v4*dy
      dy4 = amu(1,mp+2,l) + v1*dz
      dx5 = amu(2,mp+2,l) + v2*dz
      dy5 = amu(3,mp+2,l) + v3*dz
      dx6 = amu(4,mp+2,l) + v4*dz
      amu(1,mp,l) = dx1
      amu(2,mp,l) = dy1
      amu(3,mp,l) = amy
      amu(4,mp,l) = dx2
      amu(1,mp+1,l) = dy2
      amu(2,mp+1,l) = dx3
      amu(3,mp+1,l) = dy3
      amu(4,mp+1,l) = dx4
      amu(1,mp+2,l) = dy4
      amu(2,mp+2,l) = dx5
      amu(3,mp+2,l) = dy5
      amu(4,mp+2,l) = dx6
      dx1 = dcu(1,mp,l) + vx*dx
      dy1 = dcu(2,mp,l) + vy*dx
      amy = dcu(3,mp,l) + vz*dx
      dx2 = dcu(1,mp+1,l) + vx*dy
      dy2 = dcu(2,mp+1,l) + vy*dy
      dx3 = dcu(3,mp+1,l) + vz*dy
      dy3 = dcu(1,mp+2,l) + vx*dz
      dx4 = dcu(2,mp+2,l) + vy*dz
      dy4 = dcu(3,mp+2,l) + vz*dz
      dcu(1,mp,l) = dx1
      dcu(2,mp,l) = dy1
      dcu(3,mp,l) = amy
      dcu(1,mp+1,l) = dx2
      dcu(2,mp+1,l) = dy2
      dcu(3,mp+1,l) = dx3
      dcu(1,mp+2,l) = dy3
      dcu(2,mp+2,l) = dx4
      dcu(3,mp+2,l) = dy4
      dx1 = cu(1,mp,l) + ox*dx
      dy1 = cu(2,mp,l) + oy*dx
      amy = cu(3,mp,l) + oz*dx
      dx2 = cu(1,mp+1,l) + ox*dy
      dy2 = cu(2,mp+1,l) + oy*dy
      dx3 = cu(3,mp+1,l) + oz*dy
      dy3 = cu(1,mp+2,l) + ox*dz
      dx4 = cu(2,mp+2,l) + oy*dz
      dy4 = cu(3,mp+2,l) + oz*dz
      cu(1,mp,l) = dx1
      cu(2,mp,l) = dy1
      cu(3,mp,l) = amy
      cu(1,mp+1,l) = dx2
      cu(2,mp+1,l) = dy2
      cu(3,mp+1,l) = dx3
      cu(1,mp+2,l) = dy3
      cu(2,mp+2,l) = dx4
      cu(3,mp+2,l) = dy4
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      ml = mm + nn
      amx = .75 - dxn*dxn
      amy = .75 - dyn*dyn
      mn = ml + nxv
      dxl = .5*(.5 - dxn)**2
      dxp = .5*(.5 + dxn)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l))
      dz = amy*(dxl*fxy(3,mn,l) + amx*fxy(3,mn+1,l) + dxp*fxy(3,mn+2,l))
     1 + dyl*(dxl*fxy(3,ml,l) + amx*fxy(3,ml+1,l) + dxp*fxy(3,ml+2,l)) +
     2 dyp*(dxl*fxy(3,mp,l) + amx*fxy(3,mp+1,l) + dxp*fxy(3,mp+2,l)) 
c find magnetic field
      ox = amy*(dxl*bxy(1,mn,l) + amx*bxy(1,mn+1,l) + dxp*bxy(1,mn+2,l))
     1 + dyl*(dxl*bxy(1,ml,l) + amx*bxy(1,ml+1,l) + dxp*bxy(1,ml+2,l)) +
     2 dyp*(dxl*bxy(1,mp,l) + amx*bxy(1,mp+1,l) + dxp*bxy(1,mp+2,l)) 
      oy = amy*(dxl*bxy(2,mn,l) + amx*bxy(2,mn+1,l) + dxp*bxy(2,mn+2,l))
     1 + dyl*(dxl*bxy(2,ml,l) + amx*bxy(2,ml+1,l) + dxp*bxy(2,ml+2,l)) +
     2 dyp*(dxl*bxy(2,mp,l) + amx*bxy(2,mp+1,l) + dxp*bxy(2,mp+2,l)) 
      oz = amy*(dxl*bxy(3,mn,l) + amx*bxy(3,mn+1,l) + dxp*bxy(3,mn+2,l))
     1 + dyl*(dxl*bxy(3,ml,l) + amx*bxy(3,ml+1,l) + dxp*bxy(3,ml+2,l)) +
     2 dyp*(dxl*bxy(3,mp,l) + amx*bxy(3,mp+1,l) + dxp*bxy(3,mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      vz = part(5,nop,l)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      amu(1,mn,l) = amu(1,mn,l) + v1*dx
      amu(2,mn,l) = amu(2,mn,l) + v2*dx
      amu(3,mn,l) = amu(3,mn,l) + v3*dx
      amu(4,mn,l) = amu(4,mn,l) + v4*dx
      amu(1,mn+1,l) = amu(1,mn+1,l) + v1*dy
      amu(2,mn+1,l) = amu(2,mn+1,l) + v2*dy
      amu(3,mn+1,l) = amu(3,mn+1,l) + v3*dy
      amu(4,mn+1,l) = amu(4,mn+1,l) + v4*dy
      amu(1,mn+2,l) = amu(1,mn+2,l) + v1*dz
      amu(2,mn+2,l) = amu(2,mn+2,l) + v2*dz
      amu(3,mn+2,l) = amu(3,mn+2,l) + v3*dz
      amu(4,mn+2,l) = amu(4,mn+2,l) + v4*dz
      dcu(1,mn,l) = dcu(1,mn,l) + vx*dx
      dcu(2,mn,l) = dcu(2,mn,l) + vy*dx
      dcu(3,mn,l) = dcu(3,mn,l) + vz*dx
      dcu(1,mn+1,l) = dcu(1,mn+1,l) + vx*dy
      dcu(2,mn+1,l) = dcu(2,mn+1,l) + vy*dy
      dcu(3,mn+1,l) = dcu(3,mn+1,l) + vz*dy
      dcu(1,mn+2,l) = dcu(1,mn+2,l) + vx*dz
      dcu(2,mn+2,l) = dcu(2,mn+2,l) + vy*dz
      dcu(3,mn+2,l) = dcu(3,mn+2,l) + vz*dz
      cu(1,mn,l) = cu(1,mn,l) + ox*dx
      cu(2,mn,l) = cu(2,mn,l) + oy*dx
      cu(3,mn,l) = cu(3,mn,l) + oz*dx
      cu(1,mn+1,l) = cu(1,mn+1,l) + ox*dy
      cu(2,mn+1,l) = cu(2,mn+1,l) + oy*dy
      cu(3,mn+1,l) = cu(3,mn+1,l) + oz*dy
      cu(1,mn+2,l) = cu(1,mn+2,l) + ox*dz
      cu(2,mn+2,l) = cu(2,mn+2,l) + oy*dz
      cu(3,mn+2,l) = cu(3,mn+2,l) + oz*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      amu(1,ml,l) = amu(1,ml,l) + v1*dx
      amu(2,ml,l) = amu(2,ml,l) + v2*dx
      amu(3,ml,l) = amu(3,ml,l) + v3*dx
      amu(4,ml,l) = amu(4,ml,l) + v4*dx
      amu(1,ml+1,l) = amu(1,ml+1,l) + v1*dy
      amu(2,ml+1,l) = amu(2,ml+1,l) + v2*dy
      amu(3,ml+1,l) = amu(3,ml+1,l) + v3*dy
      amu(4,ml+1,l) = amu(4,ml+1,l) + v4*dy
      amu(1,ml+2,l) = amu(1,ml+2,l) + v1*dz
      amu(2,ml+2,l) = amu(2,ml+2,l) + v2*dz
      amu(3,ml+2,l) = amu(3,ml+2,l) + v3*dz
      amu(4,ml+2,l) = amu(4,ml+2,l) + v4*dz
      dcu(1,ml,l) = dcu(1,ml,l) + vx*dx
      dcu(2,ml,l) = dcu(2,ml,l) + vy*dx
      dcu(3,ml,l) = dcu(3,ml,l) + vz*dx
      dcu(1,ml+1,l) = dcu(1,ml+1,l) + vx*dy
      dcu(2,ml+1,l) = dcu(2,ml+1,l) + vy*dy
      dcu(3,ml+1,l) = dcu(3,ml+1,l) + vz*dy
      dcu(1,ml+2,l) = dcu(1,ml+2,l) + vx*dz
      dcu(2,ml+2,l) = dcu(2,ml+2,l) + vy*dz
      dcu(3,ml+2,l) = dcu(3,ml+2,l) + vz*dz
      cu(1,ml,l) = cu(1,ml,l) + ox*dx
      cu(2,ml,l) = cu(2,ml,l) + oy*dx
      cu(3,ml,l) = cu(3,ml,l) + oz*dx
      cu(1,ml+1,l) = cu(1,ml+1,l) + ox*dy
      cu(2,ml+1,l) = cu(2,ml+1,l) + oy*dy
      cu(3,ml+1,l) = cu(3,ml+1,l) + oz*dy
      cu(1,ml+2,l) = cu(1,ml+2,l) + ox*dz
      cu(2,ml+2,l) = cu(2,ml+2,l) + oy*dz
      cu(3,ml+2,l) = cu(3,ml+2,l) + oz*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      amu(1,mp,l) = amu(1,mp,l) + v1*dx
      amu(2,mp,l) = amu(2,mp,l) + v2*dx
      amu(3,mp,l) = amu(3,mp,l) + v3*dx
      amu(4,mp,l) = amu(4,mp,l) + v4*dx
      amu(1,mp+1,l) = amu(1,mp+1,l) + v1*dy
      amu(2,mp+1,l) = amu(2,mp+1,l) + v2*dy
      amu(3,mp+1,l) = amu(3,mp+1,l) + v3*dy
      amu(4,mp+1,l) = amu(4,mp+1,l) + v4*dy
      amu(1,mp+2,l) = amu(1,mp+2,l) + v1*dz
      amu(2,mp+2,l) = amu(2,mp+2,l) + v2*dz
      amu(3,mp+2,l) = amu(3,mp+2,l) + v3*dz
      amu(4,mp+2,l) = amu(4,mp+2,l) + v4*dz
      dcu(1,mp,l) = dcu(1,mp,l) + vx*dx
      dcu(2,mp,l) = dcu(2,mp,l) + vy*dx
      dcu(3,mp,l) = dcu(3,mp,l) + vz*dx
      dcu(1,mp+1,l) = dcu(1,mp+1,l) + vx*dy
      dcu(2,mp+1,l) = dcu(2,mp+1,l) + vy*dy
      dcu(3,mp+1,l) = dcu(3,mp+1,l) + vz*dy
      dcu(1,mp+2,l) = dcu(1,mp+2,l) + vx*dz
      dcu(2,mp+2,l) = dcu(2,mp+2,l) + vy*dz
      dcu(3,mp+2,l) = dcu(3,mp+2,l) + vz*dz
      cu(1,mp,l) = cu(1,mp,l) + ox*dx
      cu(2,mp,l) = cu(2,mp,l) + oy*dx
      cu(3,mp,l) = cu(3,mp,l) + oz*dx
      cu(1,mp+1,l) = cu(1,mp+1,l) + ox*dy
      cu(2,mp+1,l) = cu(2,mp+1,l) + oy*dy
      cu(3,mp+1,l) = cu(3,mp+1,l) + oz*dy
      cu(1,mp+2,l) = cu(1,mp+2,l) + ox*dz
      cu(2,mp+2,l) = cu(2,mp+2,l) + oy*dz
      cu(3,mp+2,l) = cu(3,mp+2,l) + oz*dz
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGMJPOST22(part,amu,npp,noff,qm,idimp,npmax,nblok,nxv,n
     1ypmx)
c for 2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation, for distributed data
c scalar version using guard cells
c 74 flops/particle, 22 loads, 18 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vj*vk, where jk = xx-yy,xy, for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c amu(i,j+1,k,l) = ith component of momentum flux at grid point (j,kk)
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of flux array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, amu, qm
      dimension part(idimp,npmax,nblok), amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, nl, np, ml, mp
      real qmh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy, v1, v2
      qmh = .5*qm
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      amx = qm*(.75 - dxp*dxp)
      ml = mm - mnoff
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = qmh*(.5 - dxp)**2
      np = nl + 2
      dxp = qmh*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j,l)
      vy = part(4,j,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      amu(1,nl,mm,l) = amu(1,nl,mm,l) + v1*dx
      amu(2,nl,mm,l) = amu(2,nl,mm,l) + v2*dx
      dx = dxl*dyl
      amu(1,nn,mm,l) = amu(1,nn,mm,l) + v1*dy
      amu(2,nn,mm,l) = amu(2,nn,mm,l) + v2*dy
      dy = amx*dyl
      amu(1,np,mm,l) = amu(1,np,mm,l) + v1*dz
      amu(2,np,mm,l) = amu(2,np,mm,l) + v2*dz
      dz = dxp*dyl
      amu(1,nl,ml,l) = amu(1,nl,ml,l) + v1*dx
      amu(2,nl,ml,l) = amu(2,nl,ml,l) + v2*dx
      dx = dxl*dyp
      amu(1,nn,ml,l) = amu(1,nn,ml,l) + v1*dy
      amu(2,nn,ml,l) = amu(2,nn,ml,l) + v2*dy
      dy = amx*dyp
      amu(1,np,ml,l) = amu(1,np,ml,l) + v1*dz
      amu(2,np,ml,l) = amu(2,np,ml,l) + v2*dz
      dz = dxp*dyp
      amu(1,nl,mp,l) = amu(1,nl,mp,l) + v1*dx
      amu(2,nl,mp,l) = amu(2,nl,mp,l) + v2*dx
      amu(1,nn,mp,l) = amu(1,nn,mp,l) + v1*dy
      amu(2,nn,mp,l) = amu(2,nn,mp,l) + v2*dy
      amu(1,np,mp,l) = amu(1,np,mp,l) + v1*dz
      amu(2,np,mp,l) = amu(2,np,mp,l) + v2*dz
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSMJPOST22(part,amu,npp,noff,qm,idimp,npmax,nblok,nxv,
     1nxyp)
c for 2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation, and distributed data.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 74 flops/particle, 22 loads, 18 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vj*vk, where jk = xx-yy,xy, for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c amu(i,n,l) = ith component of momentum flux at grid point (j,kk),
c where n = j + 1 + nxv*k, and kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyp = first actual dimension of current array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nxyp
      real part, amu, qm
      dimension part(idimp,npmax,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      integer nop, nnn, mmn, j, l, mnoff, nn, mm, ml, mn, mp
      real dxn, dyn, qmh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy, v1, v2
      real dx1, dy1, dx2, dy2
      qmh = .5*qm
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l) + .5
      mmn = part(2,1,l) + .5
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp(l)
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j,l) + .5
      mmn = part(2,j,l) + .5
      dxp = dxn
      dyp = dyn
      dxn = part(1,j,l) - float(nnn)
      dyn = part(2,j,l) - float(mmn)
      ml = mm + nn
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      mn = ml + nxv
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mmn = mmn - mnoff
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j-1,l)
      vy = part(4,j-1,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      dx1 = amu(1,mn,l) + v1*dx
      dy1 = amu(2,mn,l) + v2*dx
      dx2 = amu(1,mn+1,l) + v1*dy
      dy2 = amu(2,mn+1,l) + v2*dy
      dx = amu(1,mn+2,l) + v1*dz
      dy = amu(2,mn+2,l) + v2*dz
      amu(1,mn,l) = dx1
      amu(2,mn,l) = dy1
      amu(1,mn+1,l) = dx2
      amu(2,mn+1,l) = dy2
      amu(1,mn+2,l) = dx
      amu(2,mn+2,l) = dy
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = amu(1,ml,l) + v1*dx
      dy1 = amu(2,ml,l) + v2*dx
      dx2 = amu(1,ml+1,l) + v1*dy
      dy2 = amu(2,ml+1,l) + v2*dy
      dx = amu(1,ml+2,l) + v1*dz
      dy = amu(2,ml+2,l) + v2*dz
      amu(1,ml,l) = dx1
      amu(2,ml,l) = dy1
      amu(1,ml+1,l) = dx2
      amu(2,ml+1,l) = dy2
      amu(1,ml+2,l) = dx
      amu(2,ml+2,l) = dy
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = amu(1,mp,l) + v1*dx
      dy1 = amu(2,mp,l) + v2*dx
      dxl = amu(1,mp+1,l) + v1*dy
      amx = amu(2,mp+1,l) + v2*dy
      dx = amu(1,mp+2,l) + v1*dz
      dy = amu(2,mp+2,l) + v2*dz
      amu(1,mp,l) = dx1
      amu(2,mp,l) = dy1
      amu(1,mp+1,l) = dxl
      amu(2,mp+1,l) = amx
      amu(1,mp+2,l) = dx
      amu(2,mp+2,l) = dy
   10 continue
      nop = npp(l)
c deposit momentum flux for last particle
      nn = nnn + 1
      mm = nxv*mmn
      amx = qm*(.75 - dxn*dxn)
      amy = .75 - dyn*dyn
      ml = mm + nn
      dxl = qmh*(.5 - dxn)**2
      dxp = qmh*(.5 + dxn)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
      mp = mn + nxv
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      amu(1,mn,l) = amu(1,mn,l) + v1*dx
      amu(2,mn,l) = amu(2,mn,l) + v2*dx
      amu(1,mn+1,l) = amu(1,mn+1,l) + v1*dy
      amu(2,mn+1,l) = amu(2,mn+1,l) + v2*dy
      amu(1,mn+2,l) = amu(1,mn+2,l) + v1*dz
      amu(2,mn+2,l) = amu(2,mn+2,l) + v2*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      amu(1,ml,l) = amu(1,ml,l) + v1*dx
      amu(2,ml,l) = amu(2,ml,l) + v2*dx
      amu(1,ml+1,l) = amu(1,ml+1,l) + v1*dy
      amu(2,ml+1,l) = amu(2,ml+1,l) + v2*dy
      amu(1,ml+2,l) = amu(1,ml+2,l) + v1*dz
      amu(2,ml+2,l) = amu(2,ml+2,l) + v2*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      amu(1,mp,l) = amu(1,mp,l) + v1*dx
      amu(2,mp,l) = amu(2,mp,l) + v2*dx
      amu(1,mp+1,l) = amu(1,mp+1,l) + v1*dy
      amu(2,mp+1,l) = amu(2,mp+1,l) + v2*dy
      amu(1,mp+2,l) = amu(1,mp+2,l) + v1*dz
      amu(2,mp+2,l) = amu(2,mp+2,l) + v2*dz
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGDCJPOST22(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,dt,i
     1dimp,npmax,nblok,nxv,nypmx)
c for 2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation.
c scalar version using guard cells, for distributed data
c 235 flops/particle, 1 divide, 58 loads, 54 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*vj, where j = x,y, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c dcu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c dcu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c dcu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c dcu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c dcu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c dcu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c dcu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c dcu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*dvj/dt, where j = x,y, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*vj*vk, where jk = xx-yy,xy, for i = 1, 2
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = velocity vx of particle n at t - dt/2 in partition l
c part(4,n,l) = velocity vy of particle n at t - dt/2 in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bz(j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j+1,k,l) = ith component of current density
c at grid point j,kk for i = 1, 2
c dcu(i,j+1,k,l) = ith component of acceleration density
c at grid point j,kk for i = 1, 2
c amu(i,j+1,k,l) = ith component of momentum flux
c at grid point j,kk for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension cu(2,nxv,nypmx,nblok), dcu(2,nxv,nypmx,nblok)
      dimension amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, nl, np, ml, mp
      real qtmh, dti, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, ox, oy, oz
      real acx, acy, omzt, omt, anorm, rot1, rot2
      real vx, vy, v1, v2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn + 1
      ml = mm - mnoff
      amx = .75 - dxp*dxp
      mm = mm + 2
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = .5*(.5 - dxp)**2
      np = nl + 2
      dxp = .5*(.5 + dxp)**2
      mm = ml + 1
      dyl = .5*(.5 - dyp)**2
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))
c find magnetic field
      oz = amy*(dxl*bz(nl,mm,l) + amx*bz(nn,mm,l) + dxp*bz(np,mm,l)) + d
     1yl*(dxl*bz(nl,ml,l) + amx*bz(nn,ml,l) + dxp*bz(np,ml,l)) + dyp*(dx
     2l*bz(nl,mp,l) + amx*bz(nn,mp,l) + dxp*bz(np,mp,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,j,l)
      vy = part(4,j,l)
      acx = vx + dx
      acy = vy + dy
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      amu(1,nl,mm,l) = amu(1,nl,mm,l) + v1*dx
      amu(2,nl,mm,l) = amu(2,nl,mm,l) + v2*dx
      dcu(1,nl,mm,l) = dcu(1,nl,mm,l) + vx*dx
      dcu(2,nl,mm,l) = dcu(2,nl,mm,l) + vy*dx
      cu(1,nl,mm,l) = cu(1,nl,mm,l) + ox*dx
      cu(2,nl,mm,l) = cu(2,nl,mm,l) + oy*dx
      dx = dxl*dyl
      amu(1,nn,mm,l) = amu(1,nn,mm,l) + v1*dy
      amu(2,nn,mm,l) = amu(2,nn,mm,l) + v2*dy
      dcu(1,nn,mm,l) = dcu(1,nn,mm,l) + vx*dy
      dcu(2,nn,mm,l) = dcu(2,nn,mm,l) + vy*dy
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + ox*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + oy*dy
      dy = amx*dyl
      amu(1,np,mm,l) = amu(1,np,mm,l) + v1*dz
      amu(2,np,mm,l) = amu(2,np,mm,l) + v2*dz
      dcu(1,np,mm,l) = dcu(1,np,mm,l) + vx*dz
      dcu(2,np,mm,l) = dcu(2,np,mm,l) + vy*dz
      cu(1,np,mm,l) = cu(1,np,mm,l) + ox*dz
      cu(2,np,mm,l) = cu(2,np,mm,l) + oy*dz
      dz = dxp*dyl
      amu(1,nl,ml,l) = amu(1,nl,ml,l) + v1*dx
      amu(2,nl,ml,l) = amu(2,nl,ml,l) + v2*dx
      dcu(1,nl,ml,l) = dcu(1,nl,ml,l) + vx*dx
      dcu(2,nl,ml,l) = dcu(2,nl,ml,l) + vy*dx
      cu(1,nl,ml,l) = cu(1,nl,ml,l) + ox*dx
      cu(2,nl,ml,l) = cu(2,nl,ml,l) + oy*dx
      dx = dxl*dyp
      amu(1,nn,ml,l) = amu(1,nn,ml,l) + v1*dy
      amu(2,nn,ml,l) = amu(2,nn,ml,l) + v2*dy
      dcu(1,nn,ml,l) = dcu(1,nn,ml,l) + vx*dy
      dcu(2,nn,ml,l) = dcu(2,nn,ml,l) + vy*dy
      cu(1,nn,ml,l) = cu(1,nn,ml,l) + ox*dy
      cu(2,nn,ml,l) = cu(2,nn,ml,l) + oy*dy
      dy = amx*dyp
      amu(1,np,ml,l) = amu(1,np,ml,l) + v1*dz
      amu(2,np,ml,l) = amu(2,np,ml,l) + v2*dz
      dcu(1,np,ml,l) = dcu(1,np,ml,l) + vx*dz
      dcu(2,np,ml,l) = dcu(2,np,ml,l) + vy*dz
      cu(1,np,ml,l) = cu(1,np,ml,l) + ox*dz
      cu(2,np,ml,l) = cu(2,np,ml,l) + oy*dz
      dz = dxp*dyp
      amu(1,nl,mp,l) = amu(1,nl,mp,l) + v1*dx
      amu(2,nl,mp,l) = amu(2,nl,mp,l) + v2*dx
      dcu(1,nl,mp,l) = dcu(1,nl,mp,l) + vx*dx
      dcu(2,nl,mp,l) = dcu(2,nl,mp,l) + vy*dx
      cu(1,nl,mp,l) = cu(1,nl,mp,l) + ox*dx
      cu(2,nl,mp,l) = cu(2,nl,mp,l) + oy*dx
      amu(1,nn,mp,l) = amu(1,nn,mp,l) + v1*dy
      amu(2,nn,mp,l) = amu(2,nn,mp,l) + v2*dy
      dcu(1,nn,mp,l) = dcu(1,nn,mp,l) + vx*dy
      dcu(2,nn,mp,l) = dcu(2,nn,mp,l) + vy*dy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + ox*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + oy*dy
      amu(1,np,mp,l) = amu(1,np,mp,l) + v1*dz
      amu(2,np,mp,l) = amu(2,np,mp,l) + v2*dz
      dcu(1,np,mp,l) = dcu(1,np,mp,l) + vx*dz
      dcu(2,np,mp,l) = dcu(2,np,mp,l) + vy*dz
      cu(1,np,mp,l) = cu(1,np,mp,l) + ox*dz
      cu(2,np,mp,l) = cu(2,np,mp,l) + oy*dz
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSDCJPOST22(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,dt,
     1idimp,npmax,nblok,nxv,nxyp)
c for 2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 235 flops/particle, 1 divide, 58 loads, 54 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*vj, where j = x,y, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c dcu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c dcu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c dcu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c dcu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c dcu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c dcu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c dcu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c dcu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*dvj/dt, where j = x,y, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*vj*vk, where jk = xx-yy,xy, for i = 1, 2
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = velocity vx of particle n at t - dt/2 in partition l
c part(4,n,l) = velocity vy of particle n at t - dt/2 in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bz(j+1,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j+1,k,l) = ith component of current density
c at grid point j,kk for i = 1, 2
c dcu(i,j+1,k,l) = ith component of acceleration density
c at grid point j,kk for i = 1, 2
c amu(i,j+1,k,l) = ith component of momentum flux
c at grid point j,kk for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+3
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nxyp
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension cu(2,nxyp,nblok), dcu(2,nxyp,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      integer nop1, j, l, mnoff, nop, nnn, mmn, nn, mm, ml, mn, mp
      real qtmh, dti, dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, ox, oy, oz
      real acx, acy, omzt, omt, anorm, rot1, rot2
      real vx, vy, v1, v2
      real dx1, dy1, dx2, dy2, dx3, dy3, dx4, dy4, dx5
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l) + .5
      mmn = part(2,1,l) + .5
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l) + .5
      mmn = part(2,j+1,l) + .5
      dx = dxn
      dy = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      ml = mm + nn
      amx = .75 - dx*dx
      amy = .75 - dy*dy
      mn = ml + nxv
      dxl = .5*(.5 - dx)**2
      dxp = .5*(.5 + dx)**2
      mp = mn + nxv
      dyl = .5*(.5 - dy)**2
      dyp = .5*(.5 + dy)**2
      mmn = mmn - mnoff
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l)) 
c find magnetic field
      oz = amy*(dxl*bz(mn,l) + amx*bz(mn+1,l) + dxp*bz(mn+2,l)) + dyl*(d
     1xl*bz(ml,l) + amx*bz(ml+1,l) + dxp*bz(ml+2,l)) + dyp*(dxl*bz(mp,l)
     2 + amx*bz(mp+1,l) + dxp*bz(mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,j,l)
      vy = part(4,j,l)
      acx = vx + dx
      acy = vy + dy
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      dx1 = amu(1,mn,l) + v1*dx
      dy1 = amu(2,mn,l) + v2*dx
      dy2 = amu(1,mn+1,l) + v1*dy
      dx3 = amu(2,mn+1,l) + v2*dy
      dy4 = amu(1,mn+2,l) + v1*dz
      dx5 = amu(2,mn+2,l) + v2*dz
      amu(1,mn,l) = dx1
      amu(2,mn,l) = dy1
      amu(1,mn+1,l) = dy2
      amu(2,mn+1,l) = dx3
      amu(1,mn+2,l) = dy4
      amu(2,mn+2,l) = dx5
      dx1 = dcu(1,mn,l) + vx*dx
      dy1 = dcu(2,mn,l) + vy*dx
      dx2 = dcu(1,mn+1,l) + vx*dy
      dy2 = dcu(2,mn+1,l) + vy*dy
      dy3 = dcu(1,mn+2,l) + vx*dz
      dx4 = dcu(2,mn+2,l) + vy*dz
      dcu(1,mn,l) = dx1
      dcu(2,mn,l) = dy1
      dcu(1,mn+1,l) = dx2
      dcu(2,mn+1,l) = dy2
      dcu(1,mn+2,l) = dy3
      dcu(2,mn+2,l) = dx4
      dx1 = cu(1,mn,l) + ox*dx
      dy1 = cu(2,mn,l) + oy*dx
      dx2 = cu(1,mn+1,l) + ox*dy
      dy2 = cu(2,mn+1,l) + oy*dy
      dy3 = cu(1,mn+2,l) + ox*dz
      dx4 = cu(2,mn+2,l) + oy*dz
      cu(1,mn,l) = dx1
      cu(2,mn,l) = dy1
      cu(1,mn+1,l) = dx2
      cu(2,mn+1,l) = dy2
      cu(1,mn+2,l) = dy3
      cu(2,mn+2,l) = dx4
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = amu(1,ml,l) + v1*dx
      dy1 = amu(2,ml,l) + v2*dx
      dy2 = amu(1,ml+1,l) + v1*dy
      dx3 = amu(2,ml+1,l) + v2*dy
      dy4 = amu(1,ml+2,l) + v1*dz
      dx5 = amu(2,ml+2,l) + v2*dz
      amu(1,ml,l) = dx1
      amu(2,ml,l) = dy1
      amu(1,ml+1,l) = dy2
      amu(2,ml+1,l) = dx3
      amu(1,ml+2,l) = dy4
      amu(2,ml+2,l) = dx5
      dx1 = dcu(1,ml,l) + vx*dx
      dy1 = dcu(2,ml,l) + vy*dx
      dx2 = dcu(1,ml+1,l) + vx*dy
      dy2 = dcu(2,ml+1,l) + vy*dy
      dy3 = dcu(1,ml+2,l) + vx*dz
      dx4 = dcu(2,ml+2,l) + vy*dz
      dcu(1,ml,l) = dx1
      dcu(2,ml,l) = dy1
      dcu(1,ml+1,l) = dx2
      dcu(2,ml+1,l) = dy2
      dcu(1,ml+2,l) = dy3
      dcu(2,ml+2,l) = dx4
      dx1 = cu(1,ml,l) + ox*dx
      dy1 = cu(2,ml,l) + oy*dx
      dx2 = cu(1,ml+1,l) + ox*dy
      dy2 = cu(2,ml+1,l) + oy*dy
      dy3 = cu(1,ml+2,l) + ox*dz
      dx4 = cu(2,ml+2,l) + oy*dz
      cu(1,ml,l) = dx1
      cu(2,ml,l) = dy1
      cu(1,ml+1,l) = dx2
      cu(2,ml+1,l) = dy2
      cu(1,ml+2,l) = dy3
      cu(2,ml+2,l) = dx4
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = amu(1,mp,l) + v1*dx
      dy1 = amu(2,mp,l) + v2*dx
      dy2 = amu(1,mp+1,l) + v1*dy
      dx3 = amu(2,mp+1,l) + v2*dy
      dy4 = amu(1,mp+2,l) + v1*dz
      dx5 = amu(2,mp+2,l) + v2*dz
      amu(1,mp,l) = dx1
      amu(2,mp,l) = dy1
      amu(1,mp+1,l) = dy2
      amu(2,mp+1,l) = dx3
      amu(1,mp+2,l) = dy4
      amu(2,mp+2,l) = dx5
      dx1 = dcu(1,mp,l) + vx*dx
      dy1 = dcu(2,mp,l) + vy*dx
      dx2 = dcu(1,mp+1,l) + vx*dy
      dy2 = dcu(2,mp+1,l) + vy*dy
      dy3 = dcu(1,mp+2,l) + vx*dz
      dx4 = dcu(2,mp+2,l) + vy*dz
      dcu(1,mp,l) = dx1
      dcu(2,mp,l) = dy1
      dcu(1,mp+1,l) = dx2
      dcu(2,mp+1,l) = dy2
      dcu(1,mp+2,l) = dy3
      dcu(2,mp+2,l) = dx4
      dx1 = cu(1,mp,l) + ox*dx
      dy1 = cu(2,mp,l) + oy*dx
      dx2 = cu(1,mp+1,l) + ox*dy
      dy2 = cu(2,mp+1,l) + oy*dy
      dy3 = cu(1,mp+2,l) + ox*dz
      dx4 = cu(2,mp+2,l) + oy*dz
      cu(1,mp,l) = dx1
      cu(2,mp,l) = dy1
      cu(1,mp+1,l) = dx2
      cu(2,mp+1,l) = dy2
      cu(1,mp+2,l) = dy3
      cu(2,mp+2,l) = dx4
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      amx = .75 - dxn*dxn
      amy = .75 - dyn*dyn
      ml = mm + nn
      dxl = .5*(.5 - dxn)**2
      dxp = .5*(.5 + dxn)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
      mp = mn + nxv
c find electric field
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1 + dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) +
     2 dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l)) 
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1 + dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) +
     2 dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l)) 
c find magnetic field
      oz = amy*(dxl*bz(mn,l) + amx*bz(mn+1,l) + dxp*bz(mn+2,l)) + dyl*(d
     1xl*bz(ml,l) + amx*bz(ml+1,l) + dxp*bz(ml+2,l)) + dyp*(dxl*bz(mp,l)
     2 + amx*bz(mp+1,l) + dxp*bz(mp+2,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      acx = vx + dx
      acy = vy + dy
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      amu(1,mn,l) = amu(1,mn,l) + v1*dx
      amu(2,mn,l) = amu(2,mn,l) + v2*dx
      amu(1,mn+1,l) = amu(1,mn+1,l) + v1*dy
      amu(2,mn+1,l) = amu(2,mn+1,l) + v2*dy
      amu(1,mn+2,l) = amu(1,mn+2,l) + v1*dz
      amu(2,mn+2,l) = amu(2,mn+2,l) + v2*dz
      dcu(1,mn,l) = dcu(1,mn,l) + vx*dx
      dcu(2,mn,l) = dcu(2,mn,l) + vy*dx
      dcu(1,mn+1,l) = dcu(1,mn+1,l) + vx*dy
      dcu(2,mn+1,l) = dcu(2,mn+1,l) + vy*dy
      dcu(1,mn+2,l) = dcu(1,mn+2,l) + vx*dz
      dcu(2,mn+2,l) = dcu(2,mn+2,l) + vy*dz
      cu(1,mn,l) = cu(1,mn,l) + ox*dx
      cu(2,mn,l) = cu(2,mn,l) + oy*dx
      cu(1,mn+1,l) = cu(1,mn+1,l) + ox*dy
      cu(2,mn+1,l) = cu(2,mn+1,l) + oy*dy
      cu(1,mn+2,l) = cu(1,mn+2,l) + ox*dz
      cu(2,mn+2,l) = cu(2,mn+2,l) + oy*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      amu(1,ml,l) = amu(1,ml,l) + v1*dx
      amu(2,ml,l) = amu(2,ml,l) + v2*dx
      amu(1,ml+1,l) = amu(1,ml+1,l) + v1*dy
      amu(2,ml+1,l) = amu(2,ml+1,l) + v2*dy
      amu(1,ml+2,l) = amu(1,ml+2,l) + v1*dz
      amu(2,ml+2,l) = amu(2,ml+2,l) + v2*dz
      dcu(1,ml,l) = dcu(1,ml,l) + vx*dx
      dcu(2,ml,l) = dcu(2,ml,l) + vy*dx
      dcu(1,ml+1,l) = dcu(1,ml+1,l) + vx*dy
      dcu(2,ml+1,l) = dcu(2,ml+1,l) + vy*dy
      dcu(1,ml+2,l) = dcu(1,ml+2,l) + vx*dz
      dcu(2,ml+2,l) = dcu(2,ml+2,l) + vy*dz
      cu(1,ml,l) = cu(1,ml,l) + ox*dx
      cu(2,ml,l) = cu(2,ml,l) + oy*dx
      cu(1,ml+1,l) = cu(1,ml+1,l) + ox*dy
      cu(2,ml+1,l) = cu(2,ml+1,l) + oy*dy
      cu(1,ml+2,l) = cu(1,ml+2,l) + ox*dz
      cu(2,ml+2,l) = cu(2,ml+2,l) + oy*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      amu(1,mp,l) = amu(1,mp,l) + v1*dx
      amu(2,mp,l) = amu(2,mp,l) + v2*dx
      amu(1,mp+1,l) = amu(1,mp+1,l) + v1*dy
      amu(2,mp+1,l) = amu(2,mp+1,l) + v2*dy
      amu(1,mp+2,l) = amu(1,mp+2,l) + v1*dz
      amu(2,mp+2,l) = amu(2,mp+2,l) + v2*dz
      dcu(1,mp,l) = dcu(1,mp,l) + vx*dx
      dcu(2,mp,l) = dcu(2,mp,l) + vy*dx
      dcu(1,mp+1,l) = dcu(1,mp+1,l) + vx*dy
      dcu(2,mp+1,l) = dcu(2,mp+1,l) + vy*dy
      dcu(1,mp+2,l) = dcu(1,mp+2,l) + vx*dz
      dcu(2,mp+2,l) = dcu(2,mp+2,l) + vy*dz
      cu(1,mp,l) = cu(1,mp,l) + ox*dx
      cu(2,mp,l) = cu(2,mp,l) + oy*dx
      cu(1,mp+1,l) = cu(1,mp+1,l) + ox*dy
      cu(2,mp+1,l) = cu(2,mp+1,l) + oy*dy
      cu(1,mp+2,l) = cu(1,mp+2,l) + ox*dz
      cu(2,mp+2,l) = cu(2,mp+2,l) + oy*dz
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGMJPOST2L(part,amu,npp,noff,qm,idimp,npmax,nblok,nxv,n
     1ypmx)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using first-order spline interpolation
c scalar version using guard cells, for distributed data
c 51 flops/particle, 21 loads, 16 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c part(5,n,l) = z velocity of particle n at t - dt/2 in partition l
c amu(i,j,k,l) = ith component of momentum flux at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of flux array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, amu, qm
      dimension part(idimp,npmax,nblok), amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, np, mp
      real dxp, dyp, amx, amy
      real dx, dy, vx, vy, vz, v1, v2, v3, v4
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = qm*(part(1,j,l) - float(nn))
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit momentum flux
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      amu(1,np,mp,l) = amu(1,np,mp,l) + v1*dx
      amu(2,np,mp,l) = amu(2,np,mp,l) + v2*dx
      amu(3,np,mp,l) = amu(3,np,mp,l) + v3*dx
      amu(4,np,mp,l) = amu(4,np,mp,l) + v4*dx
      dx = dxp*amy
      amu(1,nn,mp,l) = amu(1,nn,mp,l) + v1*dy
      amu(2,nn,mp,l) = amu(2,nn,mp,l) + v2*dy
      amu(3,nn,mp,l) = amu(3,nn,mp,l) + v3*dy
      amu(4,nn,mp,l) = amu(4,nn,mp,l) + v4*dy
      dy = amx*amy
      amu(1,np,mm,l) = amu(1,np,mm,l) + v1*dx
      amu(2,np,mm,l) = amu(2,np,mm,l) + v2*dx
      amu(3,np,mm,l) = amu(3,np,mm,l) + v3*dx
      amu(4,np,mm,l) = amu(4,np,mm,l) + v4*dx
      amu(1,nn,mm,l) = amu(1,nn,mm,l) + v1*dy
      amu(2,nn,mm,l) = amu(2,nn,mm,l) + v2*dy
      amu(3,nn,mm,l) = amu(3,nn,mm,l) + v3*dy
      amu(4,nn,mm,l) = amu(4,nn,mm,l) + v4*dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSMJPOST2L(part,amu,npp,noff,qm,idimp,npmax,nblok,nxv,
     1nxyp)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using first-order linear interpolation
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 51 flops/particle, 21 loads, 16 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in in partition l partition l
c part(5,n,l) = z velocity of particle n at t - dt/2
c amu(i,n,l) = ith component of momentum flux at grid point j,kk
c where n = j + nxv*(k-1) and kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyp = actual first dimension of current array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nxyp
      real part, amu, qm
      dimension part(idimp,npmax,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      integer nnn, mmn, nn, mm, mp, j, l, mnoff, nop
      real dxn, dyn, dxp, dyp, amx, amy, dx, dy, dz, vx, vy, vz
      real v1, v2, v3, v4, dx1, dy1
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l)
      mmn = part(2,1,l)
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp(l)
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j,l)
      mmn = part(2,j,l)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j,l) - float(nnn)
      dyn = part(2,j,l) - float(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyp
      mmn = mmn - mnoff
c deposit momentum flux
      dx = dxp*dyp
      dz = amx*dyp
      vx = part(3,j-1,l)
      vy = part(4,j-1,l)
      vz = part(5,j-1,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      dx1 = amu(1,mp+1,l) + v1*dx
      dy1 = amu(2,mp+1,l) + v2*dx
      dyp = amu(3,mp+1,l) + v3*dx
      vx = amu(4,mp+1,l) + v4*dx
      dx = amu(1,mp,l) + v1*dz
      dy = amu(2,mp,l) + v2*dz
      vy = amu(3,mp,l) + v3*dz
      dz = amu(4,mp,l) + v4*dz
      amu(1,mp+1,l) = dx1
      amu(2,mp+1,l) = dy1
      amu(3,mp+1,l) = dyp
      amu(4,mp+1,l) = vx
      amu(1,mp,l) = dx
      amu(2,mp,l) = dy
      amu(3,mp,l) = vy
      amu(4,mp,l) = dz
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1,l) + v1*dx
      amx = amu(2,mm+1,l) + v2*dx
      dyp = amu(3,mm+1,l) + v3*dx
      vx = amu(4,mm+1,l) + v4*dx
      dx = amu(1,mm,l) + v1*dz
      dy = amu(2,mm,l) + v2*dz
      vy = amu(3,mm,l) + v3*dz
      dz = amu(4,mm,l) + v4*dz
      amu(1,mm+1,l) = dxp
      amu(2,mm+1,l) = amx
      amu(3,mm+1,l) = dyp
      amu(4,mm+1,l) = vx
      amu(1,mm,l) = dx
      amu(2,mm,l) = dy
      amu(3,mm,l) = vy
      amu(4,mm,l) = dz
   10 continue
      nop = npp(l)
c deposit momentum flux for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit momentum flux
      dx = dxp*dyn
      dy = amx*dyn
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      vz = part(5,nop,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      v3 = vz*vx
      v4 = vz*vy
      amu(1,mp+1,l) = amu(1,mp+1,l) + v1*dx
      amu(2,mp+1,l) = amu(2,mp+1,l) + v2*dx
      amu(3,mp+1,l) = amu(3,mp+1,l) + v3*dx
      amu(4,mp+1,l) = amu(4,mp+1,l) + v4*dx
      amu(1,mp,l) = amu(1,mp,l) + v1*dy
      amu(2,mp,l) = amu(2,mp,l) + v2*dy
      amu(3,mp,l) = amu(3,mp,l) + v3*dy
      amu(4,mp,l) = amu(4,mp,l) + v4*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1,l) = amu(1,mm+1,l) + v1*dx
      amu(2,mm+1,l) = amu(2,mm+1,l) + v2*dx
      amu(3,mm+1,l) = amu(3,mm+1,l) + v3*dx
      amu(4,mm+1,l) = amu(4,mm+1,l) + v4*dx
      amu(1,mm,l) = amu(1,mm,l) + v1*dy
      amu(2,mm,l) = amu(2,mm,l) + v2*dy
      amu(3,mm,l) = amu(3,mm,l) + v3*dy
      amu(4,mm,l) = amu(4,mm,l) + v4*dy
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,dt,
     1idimp,npmax,nblok,nxv,nypmx)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation.
c scalar version using guard cells, for distributed data
c 218 flops/particle, 1 divide, 69 loads, 40 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = velocity vx of particle n at t - dt/2 in partition l
c part(4,n,l) = velocity vy of particle n at t - dt/2 in partition l
c part(5,n,l) = velocity vz of particle n at t - dt/2 in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c fxy(3,j,k,l) = z component of force/charge at grid (j,kk)
c where kk = k + noff(l) - 1
c that is, convolution of electric field over particle shape
c bxy(1,j,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j,k,l) = ith component of current density
c at grid point j,kk for i = 1, 3
c dcu(i,j,k,l) = ith component of acceleration density
c at grid point j,kk for i = 1, 3
c amu(i,j,k,l) = ith component of momentum flux
c at grid point j,kk for i = 1, 4
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension cu(3,nxv,nypmx,nblok), dcu(3,nxv,nypmx,nblok)
      dimension amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, np, mp
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp,l) + amx*fxy(1,nn,mp,l)) + amy*(dxp*fxy(
     11,np,mm,l) + amx*fxy(1,nn,mm,l))
      dy = dyp*(dxp*fxy(2,np,mp,l) + amx*fxy(2,nn,mp,l)) + amy*(dxp*fxy(
     12,np,mm,l) + amx*fxy(2,nn,mm,l))
      dz = dyp*(dxp*fxy(3,np,mp,l) + amx*fxy(3,nn,mp,l)) + amy*(dxp*fxy(
     13,np,mm,l) + amx*fxy(3,nn,mm,l))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp,l) + amx*bxy(1,nn,mp,l)) + amy*(dxp*bxy(
     11,np,mm,l) + amx*bxy(1,nn,mm,l))
      oy = dyp*(dxp*bxy(2,np,mp,l) + amx*bxy(2,nn,mp,l)) + amy*(dxp*bxy(
     12,np,mm,l) + amx*bxy(2,nn,mm,l))
      oz = dyp*(dxp*bxy(3,np,mp,l) + amx*bxy(3,nn,mp,l)) + amy*(dxp*bxy(
     13,np,mm,l) + amx*bxy(3,nn,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = dxp*dyp
      dy = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      amu(1,np,mp,l) = amu(1,np,mp,l) + v1*dx
      amu(2,np,mp,l) = amu(2,np,mp,l) + v2*dx
      amu(3,np,mp,l) = amu(3,np,mp,l) + v3*dx
      amu(4,np,mp,l) = amu(4,np,mp,l) + v4*dx
      dcu(1,np,mp,l) = dcu(1,np,mp,l) + vx*dx
      dcu(2,np,mp,l) = dcu(2,np,mp,l) + vy*dx
      dcu(3,np,mp,l) = dcu(3,np,mp,l) + vz*dx
      cu(1,np,mp,l) = cu(1,np,mp,l) + ox*dx
      cu(2,np,mp,l) = cu(2,np,mp,l) + oy*dx
      cu(3,np,mp,l) = cu(3,np,mp,l) + oz*dx
      dx = dxp*amy
      amu(1,nn,mp,l) = amu(1,nn,mp,l) + v1*dy
      amu(2,nn,mp,l) = amu(2,nn,mp,l) + v2*dy
      amu(3,nn,mp,l) = amu(3,nn,mp,l) + v3*dy
      amu(4,nn,mp,l) = amu(4,nn,mp,l) + v4*dy
      dcu(1,nn,mp,l) = dcu(1,nn,mp,l) + vx*dy
      dcu(2,nn,mp,l) = dcu(2,nn,mp,l) + vy*dy
      dcu(3,nn,mp,l) = dcu(3,nn,mp,l) + vz*dy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + ox*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + oy*dy
      cu(3,nn,mp,l) = cu(3,nn,mp,l) + oz*dy
      dy = amx*amy
      amu(1,np,mm,l) = amu(1,np,mm,l) + v1*dx
      amu(2,np,mm,l) = amu(2,np,mm,l) + v2*dx
      amu(3,np,mm,l) = amu(3,np,mm,l) + v3*dx
      amu(4,np,mm,l) = amu(4,np,mm,l) + v4*dx
      dcu(1,np,mm,l) = dcu(1,np,mm,l) + vx*dx
      dcu(2,np,mm,l) = dcu(2,np,mm,l) + vy*dx
      dcu(3,np,mm,l) = dcu(3,np,mm,l) + vz*dx
      cu(1,np,mm,l) = cu(1,np,mm,l) + ox*dx
      cu(2,np,mm,l) = cu(2,np,mm,l) + oy*dx
      cu(3,np,mm,l) = cu(3,np,mm,l) + oz*dx
      amu(1,nn,mm,l) = amu(1,nn,mm,l) + v1*dy
      amu(2,nn,mm,l) = amu(2,nn,mm,l) + v2*dy
      amu(3,nn,mm,l) = amu(3,nn,mm,l) + v3*dy
      amu(4,nn,mm,l) = amu(4,nn,mm,l) + v4*dy
      dcu(1,nn,mm,l) = dcu(1,nn,mm,l) + vx*dy
      dcu(2,nn,mm,l) = dcu(2,nn,mm,l) + vy*dy
      dcu(3,nn,mm,l) = dcu(3,nn,mm,l) + vz*dy
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + ox*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + oy*dy
      cu(3,nn,mm,l) = cu(3,nn,mm,l) + oz*dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSDCJPOST2L(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,dt
     1,idimp,npmax,nblok,nxv,nxyp)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using first-order spline
c interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 218 flops/particle, 1 divide, 69 loads, 40 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj/dt, where j = x,y,z, for i = 1, 3
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj*vk, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = velocity vx of particle n at t - dt/2 in partition l
c part(4,n,l) = velocity vy of particle n at t - dt/2 in partition l
c part(5,n,l) = velocity vz of particle n at t - dt/2 in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l,l) = y component of force/charge at grid (j,kk)
c fxy(3,j,k) = z component of force/charge at grid (j,kk)
c where kk = k + noff(l) - 1
c that is, convolution of electric field over particle shape
c bxy(1,j,k,l) = x component of magnetic field at grid (j,kk)
c bxy(2,j,k,l) = y component of magnetic field at grid (j,kk)
c bxy(3,j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,n,l) = ith component of current density at grid point j,kk
c where n = j + nxv*(k-1)
c dcu(i,n,l) = ith component of acceleration density at grid point j,kk
c where n = j + nxv*(k-1)
c amu(i,n,l) = ith component of momentum at grid point j,kk
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nxyp
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension cu(3,nxyp,nblok), dcu(3,nxyp,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      integer nnn, mmn, nop1, j, l, mnoff, nop, nn, mm, mp
      real qtmh, dti, dxn, dyn, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4, dx1, dy1, dx2, dy2, dx3, dy3
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l)
      mmn = part(2,1,l)
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l)
      mmn = part(2,j+1,l)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      mm = mm + nn
      amx = 1. - dxp
      mp = mm + nxv
      amy = 1. - dyp
      mmn = mmn - mnoff
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxp*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyp*(dxp*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxp*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
      dz = dyp*(dxp*fxy(3,mp+1,l) + amx*fxy(3,mp,l)) + amy*(dxp*fxy(3,mm
     1+1,l) + amx*fxy(3,mm,l))
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1,l) + amx*bxy(1,mp,l)) + amy*(dxp*bxy(1,mm
     1+1,l) + amx*bxy(1,mm,l))
      oy = dyp*(dxp*bxy(2,mp+1,l) + amx*bxy(2,mp,l)) + amy*(dxp*bxy(2,mm
     1+1,l) + amx*bxy(2,mm,l))
      oz = dyp*(dxp*bxy(3,mp+1,l) + amx*bxy(3,mp,l)) + amy*(dxp*bxy(3,mm
     1+1,l) + amx*bxy(3,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j,l)
      vy = part(4,j,l)
      vz = part(5,j,l)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = dxp*dyp
      dz = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      dx1 = amu(1,mp+1,l) + v1*dx
      dy1 = amu(2,mp+1,l) + v2*dx
      dyp = amu(3,mp+1,l) + v3*dx
      dx2 = amu(4,mp+1,l) + v4*dx
      dy2 = amu(1,mp,l) + v1*dz
      dx3 = amu(2,mp,l) + v2*dz
      dy3 = amu(3,mp,l) + v3*dz
      dy = amu(4,mp,l) + v4*dz
      amu(1,mp+1,l) = dx1
      amu(2,mp+1,l) = dy1
      amu(3,mp+1,l) = dyp
      amu(4,mp+1,l) = dx2
      amu(1,mp,l) = dy2
      amu(2,mp,l) = dx3
      amu(3,mp,l) = dy3
      amu(4,mp,l) = dy
      dx1 = dcu(1,mp+1,l) + vx*dx
      dy1 = dcu(2,mp+1,l) + vy*dx
      dyp = dcu(3,mp+1,l) + vz*dx
      dx2 = dcu(1,mp,l) + vx*dz
      dy2 = dcu(2,mp,l) + vy*dz
      dy = dcu(3,mp,l) + vz*dz
      dcu(1,mp+1,l) = dx1
      dcu(2,mp+1,l) = dy1
      dcu(3,mp+1,l) = dyp
      dcu(1,mp,l) = dx2
      dcu(2,mp,l) = dy2
      dcu(3,mp,l) = dy
      dx1 = cu(1,mp+1,l) + ox*dx
      dy1 = cu(2,mp+1,l) + oy*dx
      dyp = cu(3,mp+1,l) + oz*dx
      dx2 = cu(1,mp,l) + ox*dz
      dy2 = cu(2,mp,l) + oy*dz
      dy = cu(3,mp,l) + oz*dz
      cu(1,mp+1,l) = dx1
      cu(2,mp+1,l) = dy1
      cu(3,mp+1,l) = dyp
      cu(1,mp,l) = dx2
      cu(2,mp,l) = dy2
      cu(3,mp,l) = dy
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1,l) + v1*dx
      amx = amu(2,mm+1,l) + v2*dx
      dyp = amu(3,mm+1,l) + v3*dx
      dx1 = amu(4,mm+1,l) + v4*dx
      dy1 = amu(1,mm,l) + v1*dz
      dx2 = amu(2,mm,l) + v2*dz
      dy2 = amu(3,mm,l) + v3*dz
      dy = amu(4,mm,l) + v4*dz
      amu(1,mm+1,l) = dxp
      amu(2,mm+1,l) = amx
      amu(3,mm+1,l) = dyp
      amu(4,mm+1,l) = dx1
      amu(1,mm,l) = dy1
      amu(2,mm,l) = dx2
      amu(3,mm,l) = dy2
      amu(4,mm,l) = dy
      dxp = dcu(1,mm+1,l) + vx*dx
      amx = dcu(2,mm+1,l) + vy*dx
      dyp = dcu(3,mm+1,l) + vz*dx
      dx1 = dcu(1,mm,l) + vx*dz
      dy1 = dcu(2,mm,l) + vy*dz
      dy = dcu(3,mm,l) + vz*dz
      dcu(1,mm+1,l) = dxp
      dcu(2,mm+1,l) = amx
      dcu(3,mm+1,l) = dyp
      dcu(1,mm,l) = dx1
      dcu(2,mm,l) = dy1
      dcu(3,mm,l) = dy
      dxp = cu(1,mm+1,l) + ox*dx
      amx = cu(2,mm+1,l) + oy*dx
      dyp = cu(3,mm+1,l) + oz*dx
      dx1 = cu(1,mm,l) + ox*dz
      dy1 = cu(2,mm,l) + oy*dz
      dy = cu(3,mm,l) + oz*dz
      cu(1,mm+1,l) = dxp
      cu(2,mm+1,l) = amx
      cu(3,mm+1,l) = dyp
      cu(1,mm,l) = dx1
      cu(2,mm,l) = dy1
      cu(3,mm,l) = dy
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxn*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyn*(dxn*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxn*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
      dz = dyn*(dxn*fxy(3,mp+1,l) + amx*fxy(3,mp,l)) + amy*(dxn*fxy(3,mm
     1+1,l) + amx*fxy(3,mm,l))
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1,l) + amx*bxy(1,mp,l)) + amy*(dxn*bxy(1,mm
     1+1,l) + amx*bxy(1,mm,l))
      oy = dyn*(dxn*bxy(2,mp+1,l) + amx*bxy(2,mp,l)) + amy*(dxn*bxy(2,mm
     1+1,l) + amx*bxy(2,mm,l))
      oz = dyn*(dxn*bxy(3,mp+1,l) + amx*bxy(3,mp,l)) + amy*(dxn*bxy(3,mm
     1+1,l) + amx*bxy(3,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      vz = part(5,nop,l)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxn
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      oz = 0.5*(dz + vz)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      vz = dti*(dz - vz)
      dx = dxp*dyn
      dy = amx*dyn
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      amu(1,mp+1,l) = amu(1,mp+1,l) + v1*dx
      amu(2,mp+1,l) = amu(2,mp+1,l) + v2*dx
      amu(3,mp+1,l) = amu(3,mp+1,l) + v3*dx
      amu(4,mp+1,l) = amu(4,mp+1,l) + v4*dx
      amu(1,mp,l) = amu(1,mp,l) + v1*dy
      amu(2,mp,l) = amu(2,mp,l) + v2*dy
      amu(3,mp,l) = amu(3,mp,l) + v3*dy
      amu(4,mp,l) = amu(4,mp,l) + v4*dy
      dcu(1,mp+1,l) = dcu(1,mp+1,l) + vx*dx
      dcu(2,mp+1,l) = dcu(2,mp+1,l) + vy*dx
      dcu(3,mp+1,l) = dcu(3,mp+1,l) + vz*dx
      dcu(1,mp,l) = dcu(1,mp,l) + vx*dy
      dcu(2,mp,l) = dcu(2,mp,l) + vy*dy
      dcu(3,mp,l) = dcu(3,mp,l) + vz*dy
      cu(1,mp+1,l) = cu(1,mp+1,l) + ox*dx
      cu(2,mp+1,l) = cu(2,mp+1,l) + oy*dx
      cu(3,mp+1,l) = cu(3,mp+1,l) + oz*dx
      cu(1,mp,l) = cu(1,mp,l) + ox*dy
      cu(2,mp,l) = cu(2,mp,l) + oy*dy
      cu(3,mp,l) = cu(3,mp,l) + oz*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1,l) = amu(1,mm+1,l) + v1*dx
      amu(2,mm+1,l) = amu(2,mm+1,l) + v2*dx
      amu(3,mm+1,l) = amu(3,mm+1,l) + v3*dx
      amu(4,mm+1,l) = amu(4,mm+1,l) + v4*dx
      amu(1,mm,l) = amu(1,mm,l) + v1*dy
      amu(2,mm,l) = amu(2,mm,l) + v2*dy
      amu(3,mm,l) = amu(3,mm,l) + v3*dy
      amu(4,mm,l) = amu(4,mm,l) + v4*dy
      dcu(1,mm+1,l) = dcu(1,mm+1,l) + vx*dx
      dcu(2,mm+1,l) = dcu(2,mm+1,l) + vy*dx
      dcu(3,mm+1,l) = dcu(3,mm+1,l) + vz*dx
      dcu(1,mm,l) = dcu(1,mm,l) + vx*dy
      dcu(2,mm,l) = dcu(2,mm,l) + vy*dy
      dcu(3,mm,l) = dcu(3,mm,l) + vz*dy
      cu(1,mm+1,l) = cu(1,mm+1,l) + ox*dx
      cu(2,mm+1,l) = cu(2,mm+1,l) + oy*dx
      cu(3,mm+1,l) = cu(3,mm+1,l) + oz*dx
      cu(1,mm,l) = cu(1,mm,l) + ox*dy
      cu(2,mm,l) = cu(2,mm,l) + oy*dy
      cu(3,mm,l) = cu(3,mm,l) + oz*dy
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGMJPOST22L(part,amu,npp,noff,qm,idimp,npmax,nblok,nxv,
     1nypmx)
c for 2d code, this subroutine calculates particle momentum flux
c using first-order spline interpolation
c scalar version using guard cells, for distributed data
c 33 flops/particle, 12 loads, 8 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vj*vk, where jk = xx-yy,xy, for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c amu(i,j,k,l) = ith component of momentum flux at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of flux array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, amu, qm
      dimension part(idimp,npmax,nblok), amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, np, mp
      real dxp, dyp, amx, amy
      real dx, dy, vx, vy, v1, v2
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = qm*(part(1,j,l) - float(nn))
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit momentum flux
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j,l)
      vy = part(4,j,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      amu(1,np,mp,l) = amu(1,np,mp,l) + v1*dx
      amu(2,np,mp,l) = amu(2,np,mp,l) + v2*dx
      dx = dxp*amy
      amu(1,nn,mp,l) = amu(1,nn,mp,l) + v1*dy
      amu(2,nn,mp,l) = amu(2,nn,mp,l) + v2*dy
      dy = amx*amy
      amu(1,np,mm,l) = amu(1,np,mm,l) + v1*dx
      amu(2,np,mm,l) = amu(2,np,mm,l) + v2*dx
      amu(1,nn,mm,l) = amu(1,nn,mm,l) + v1*dy
      amu(2,nn,mm,l) = amu(2,nn,mm,l) + v2*dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSMJPOST22L(part,amu,npp,noff,qm,idimp,npmax,nblok,nxv
     1,nxyp)
c for 2d code, this subroutine calculates particle momentum flux
c using first-order linear interpolation
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 33 flops/particle, 12 loads, 8 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vj*vk, where jk = xx-yy,xy, for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c amu(i,n,l) = ith component of momentum flux at grid point j,kk
c where n = j + nxv*(k-1) and kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyp = actual first dimension of current array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nxyp
      real part, amu, qm
      dimension part(idimp,npmax,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      integer nnn, mmn, nn, mm, mp, j, l, mnoff, nop
      real dxn, dyn, dxp, dyp, amx, amy, dx, dy, dz, vx, vy
      real v1, v2, dx1, dy1
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l)
      mmn = part(2,1,l)
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      do 10 j = 2, npp(l)
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j,l)
      mmn = part(2,j,l)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j,l) - float(nnn)
      dyn = part(2,j,l) - float(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyp
      mmn = mmn - mnoff
c deposit momentum flux
      dx = dxp*dyp
      dz = amx*dyp
      vx = part(3,j-1,l)
      vy = part(4,j-1,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      dx1 = amu(1,mp+1,l) + v1*dx
      dy1 = amu(2,mp+1,l) + v2*dx
      dx = amu(1,mp,l) + v1*dz
      dy = amu(2,mp,l) + v2*dz
      amu(1,mp+1,l) = dx1
      amu(2,mp+1,l) = dy1
      amu(1,mp,l) = dx
      amu(2,mp,l) = dy
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1,l) + v1*dx
      amx = amu(2,mm+1,l) + v2*dx
      dx = amu(1,mm,l) + v1*dz
      dy = amu(2,mm,l) + v2*dz
      amu(1,mm+1,l) = dxp
      amu(2,mm+1,l) = amx
      amu(1,mm,l) = dx
      amu(2,mm,l) = dy
   10 continue
      nop = npp(l)
c deposit momentum flux for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit momentum flux
      dx = dxp*dyn
      dy = amx*dyn
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      v1 = vx*vx - vy*vy
      v2 = vx*vy
      amu(1,mp+1,l) = amu(1,mp+1,l) + v1*dx
      amu(2,mp+1,l) = amu(2,mp+1,l) + v2*dx
      amu(1,mp,l) = amu(1,mp,l) + v1*dy
      amu(2,mp,l) = amu(2,mp,l) + v2*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1,l) = amu(1,mm+1,l) + v1*dx
      amu(2,mm+1,l) = amu(2,mm+1,l) + v2*dx
      amu(1,mm,l) = amu(1,mm,l) + v1*dy
      amu(2,mm,l) = amu(2,mm,l) + v2*dy
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGDCJPOST22L(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,dt,
     1idimp,npmax,nblok,nxv,nypmx)
c for 2d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation.
c scalar version using guard cells, for distributed data
c 120 flops/particle, 1 divide, 28 loads, 24 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj/dt, where j = x,y, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj*vk, where jk = xx-yy,xy, for i = 1, 2
c where vj = 0.5*(vj(t+dt/2)+vj(t-dt/2),
c and vk = 0.5*(vk(t+dt/2)+vk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = velocity vx of particle n at t - dt/2 in partition l
c part(4,n,l) = velocity vy of particle n at t - dt/2 in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k),l = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bz(j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j,k,l) = ith component of current density
c at grid point j,kk for i = 1, 2
c dcu(i,j,k,l) = ith component of acceleration density
c at grid point j,kk for i = 1, 2
c amu(i,j,k,l) = ith component of momentum flux
c at grid point j,kk for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension cu(2,nxv,nypmx,nblok), dcu(2,nxv,nypmx,nblok)
      dimension amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, np, mp
      real qtmh, dti, dxp, dyp, amx, amy, dx, dy, ox, oy, oz
      real acx, acy, omzt, omt, anorm, rot1, rot2
      real vx, vy, v1, v2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm - mnoff
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp,l) + amx*fxy(1,nn,mp,l)) + amy*(dxp*fxy(
     11,np,mm,l) + amx*fxy(1,nn,mm,l))
      dy = dyp*(dxp*fxy(2,np,mp,l) + amx*fxy(2,nn,mp,l)) + amy*(dxp*fxy(
     12,np,mm,l) + amx*fxy(2,nn,mm,l))
c find magnetic field
      oz = dyp*(dxp*bz(np,mp,l) + amx*bz(nn,mp,l)) + amy*(dxp*bz(np,mm,l
     1) + amx*bz(nn,mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,j,l)
      vy = part(4,j,l)
      acx = vx + dx
      acy = vy + dy
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      dx = dxp*dyp
      dy = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      amu(1,np,mp,l) = amu(1,np,mp,l) + v1*dx
      amu(2,np,mp,l) = amu(2,np,mp,l) + v2*dx
      dcu(1,np,mp,l) = dcu(1,np,mp,l) + vx*dx
      dcu(2,np,mp,l) = dcu(2,np,mp,l) + vy*dx
      cu(1,np,mp,l) = cu(1,np,mp,l) + ox*dx
      cu(2,np,mp,l) = cu(2,np,mp,l) + oy*dx
      dx = dxp*amy
      amu(1,nn,mp,l) = amu(1,nn,mp,l) + v1*dy
      amu(2,nn,mp,l) = amu(2,nn,mp,l) + v2*dy
      dcu(1,nn,mp,l) = dcu(1,nn,mp,l) + vx*dy
      dcu(2,nn,mp,l) = dcu(2,nn,mp,l) + vy*dy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + ox*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + oy*dy
      dy = amx*amy
      amu(1,np,mm,l) = amu(1,np,mm,l) + v1*dx
      amu(2,np,mm,l) = amu(2,np,mm,l) + v2*dx
      dcu(1,np,mm,l) = dcu(1,np,mm,l) + vx*dx
      dcu(2,np,mm,l) = dcu(2,np,mm,l) + vy*dx
      cu(1,np,mm,l) = cu(1,np,mm,l) + ox*dx
      cu(2,np,mm,l) = cu(2,np,mm,l) + oy*dx
      amu(1,nn,mm,l) = amu(1,nn,mm,l) + v1*dy
      amu(2,nn,mm,l) = amu(2,nn,mm,l) + v2*dy
      dcu(1,nn,mm,l) = dcu(1,nn,mm,l) + vx*dy
      dcu(2,nn,mm,l) = dcu(2,nn,mm,l) + vy*dy
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + ox*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + oy*dy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSDCJPOST22L(part,fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,dt
     1,idimp,npmax,nblok,nxv,nxyp)
c for 2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using first-order spline
c interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 120 flops/particle, 1 divide, 28 loads, 24 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c acceleration density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj/dt, where j = x,y, for i = 1, 2
c where dvj = (vj(t+dt/2)-vj(t-dt/2))/dt
c momentum flux is approximated by values at the grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj*vk, where jk = xx-yy,xy, for i = 1, 2
c where vj = vj(t-dt/2) and vk = vk(t-dt/2)
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = velocity vx of particle n at t - dt/2 in partition l
c part(4,n,l) = velocity vy of particle n at t - dt/2 in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c bz(j,k,l) = z component of magnetic field at grid (j,kk)
c that is, the convolution of magnetic field over particle shape
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,n,l) = ith component of current density at grid point j,kk
c where n = j + nxv*(k-1)
c dcu(i,n,l) = ith component of acceleration density at grid point j,kk
c where n = j + nxv*(k-1)
c amu(i,n,l) = ith component of momentum at grid point j,kk
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nxyp
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension cu(2,nxyp,nblok), dcu(2,nxyp,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      integer nnn, mmn, nop1, j, l, mnoff, nop, nn, mm, mp
      real qtmh, dti, dxn, dyn, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, omzt, omt, anorm, rot1, rot2
      real vx, vy, v1, v2, dx1, dy1, dx2, dy2, dx3
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      do 20 l = 1, nblok
      if (npp(l).lt.1) go to 20
      mnoff = noff(l)
c begin first particle
      nnn = part(1,1,l)
      mmn = part(2,1,l)
      dxn = part(1,1,l) - float(nnn)
      dyn = part(2,1,l) - float(mmn)
      mmn = mmn - mnoff
      nop1 = npp(l) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1,l)
      mmn = part(2,j+1,l)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      mm = mm + nn
      amx = 1. - dxp
      mp = mm + nxv
      amy = 1. - dyp
      mmn = mmn - mnoff
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxp*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyp*(dxp*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxp*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
c find magnetic field
      oz = dyp*(dxp*bz(mp+1,l) + amx*bz(mp,l)) + amy*(dxp*bz(mm+1,l) + a
     1mx*bz(mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,j,l)
      vy = part(4,j,l)
      acx = vx + dx
      acy = vy + dy
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      dx = dxp*dyp
      dz = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      dx1 = amu(1,mp+1,l) + v1*dx
      dy1 = amu(2,mp+1,l) + v2*dx
      dy2 = amu(1,mp,l) + v1*dz
      dx3 = amu(2,mp,l) + v2*dz
      amu(1,mp+1,l) = dx1
      amu(2,mp+1,l) = dy1
      amu(1,mp,l) = dy2
      amu(2,mp,l) = dx3
      dx1 = dcu(1,mp+1,l) + vx*dx
      dy1 = dcu(2,mp+1,l) + vy*dx
      dx2 = dcu(1,mp,l) + vx*dz
      dy2 = dcu(2,mp,l) + vy*dz
      dcu(1,mp+1,l) = dx1
      dcu(2,mp+1,l) = dy1
      dcu(1,mp,l) = dx2
      dcu(2,mp,l) = dy2
      dx1 = cu(1,mp+1,l) + ox*dx
      dy1 = cu(2,mp+1,l) + oy*dx
      dx2 = cu(1,mp,l) + ox*dz
      dy2 = cu(2,mp,l) + oy*dz
      cu(1,mp+1,l) = dx1
      cu(2,mp+1,l) = dy1
      cu(1,mp,l) = dx2
      cu(2,mp,l) = dy2
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1,l) + v1*dx
      amx = amu(2,mm+1,l) + v2*dx
      dy1 = amu(1,mm,l) + v1*dz
      dx2 = amu(2,mm,l) + v2*dz
      amu(1,mm+1,l) = dxp
      amu(2,mm+1,l) = amx
      amu(1,mm,l) = dy1
      amu(2,mm,l) = dx2
      dxp = dcu(1,mm+1,l) + vx*dx
      amx = dcu(2,mm+1,l) + vy*dx
      dx1 = dcu(1,mm,l) + vx*dz
      dy1 = dcu(2,mm,l) + vy*dz
      dcu(1,mm+1,l) = dxp
      dcu(2,mm+1,l) = amx
      dcu(1,mm,l) = dx1
      dcu(2,mm,l) = dy1
      dxp = cu(1,mm+1,l) + ox*dx
      amx = cu(2,mm+1,l) + oy*dx
      dx1 = cu(1,mm,l) + ox*dz
      dy1 = cu(2,mm,l) + oy*dz
      cu(1,mm+1,l) = dxp
      cu(2,mm+1,l) = amx
      cu(1,mm,l) = dx1
      cu(2,mm,l) = dy1
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxn*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyn*(dxn*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxn*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
c find magnetic field
      oz = dyn*(dxn*bz(mp+1,l) + amx*bz(mp,l)) + amy*(dxn*bz(mm+1,l) + a
     1mx*bz(mm,l))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,nop,l)
      vy = part(4,nop,l)
      acx = vx + dx
      acy = vy + dy
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxn
      ox = 0.5*(dx + vx)
      oy = 0.5*(dy + vy)
      vx = dti*(dx - vx)
      vy = dti*(dy - vy)
      dx = dxp*dyn
      dy = amx*dyn
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      amu(1,mp+1,l) = amu(1,mp+1,l) + v1*dx
      amu(2,mp+1,l) = amu(2,mp+1,l) + v2*dx
      amu(1,mp,l) = amu(1,mp,l) + v1*dy
      amu(2,mp,l) = amu(2,mp,l) + v2*dy
      dcu(1,mp+1,l) = dcu(1,mp+1,l) + vx*dx
      dcu(2,mp+1,l) = dcu(2,mp+1,l) + vy*dx
      dcu(1,mp,l) = dcu(1,mp,l) + vx*dy
      dcu(2,mp,l) = dcu(2,mp,l) + vy*dy
      cu(1,mp+1,l) = cu(1,mp+1,l) + ox*dx
      cu(2,mp+1,l) = cu(2,mp+1,l) + oy*dx
      cu(1,mp,l) = cu(1,mp,l) + ox*dy
      cu(2,mp,l) = cu(2,mp,l) + oy*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1,l) = amu(1,mm+1,l) + v1*dx
      amu(2,mm+1,l) = amu(2,mm+1,l) + v2*dx
      amu(1,mm,l) = amu(1,mm,l) + v1*dy
      amu(2,mm,l) = amu(2,mm,l) + v2*dy
      dcu(1,mm+1,l) = dcu(1,mm+1,l) + vx*dx
      dcu(2,mm+1,l) = dcu(2,mm+1,l) + vy*dx
      dcu(1,mm,l) = dcu(1,mm,l) + vx*dy
      dcu(2,mm,l) = dcu(2,mm,l) + vy*dy
      cu(1,mm+1,l) = cu(1,mm+1,l) + ox*dx
      cu(2,mm+1,l) = cu(2,mm+1,l) + oy*dx
      cu(1,mm,l) = cu(1,mm,l) + ox*dy
      cu(2,mm,l) = cu(2,mm,l) + oy*dy
   20 continue
      return
      end
