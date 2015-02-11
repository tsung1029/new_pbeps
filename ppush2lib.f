c-----------------------------------------------------------------------
c 2d parallel PIC library for pushing particles and depositing charge
c ppush2lib.f contains procedures to process particles:
c PGPOST2 deposits charge density, quadratic interpolation, STANDARD
c         optimization, and distributed data.
c PGSPOST2 deposits charge density, quadratic interpolation, LOOKAHEAD
c          optimization, and distributed data.
c PGSOST2X deposits charge density, quadratic interpolation, VECTOR
c          optimization, and distributed data.
c PGPOST2L deposits charge density, linear interpolation, STANDARD
c          optimization, and distributed data.
c PGSPOST2L deposits charge density, linear interpolation, LOOKAHEAD
c           optimization, and distributed data.
c PGSOST2XL deposits charge density, linear interpolation, VECTOR
c           optimization, and distributed data.
c PGPUSH2 push particles, quadratic interpolation, STANDARD
c         optimization, and distributed data.
c PGSPUSH2 push particles, quadratic interpolation, LOOKAHEAD
c          optimization, and distributed data.
c PGPUSH2L push particles, linear interpolation, STANDARD optimization,
c          and distributed data.
c PGSPUSH2L push particles, linear interpolation, LOOKAHEAD
c           optimization, and distributed data.
c PSORTP2Y sort particles by y grid, quadratic interpolation, memory
c          conserving algorithm, and distributed data.
c PSORTP2YL sort particles by y grid, linear interpolation, memory
c           conserving algorithm, and distributed data.
c PDSORTP2Y sort particles by y grid, quadratic interpolation, high
c           performance algorithm, and distributed data.
c PDSORTP2YL sort particles by y grid, linear interpolation, high
c            performance algorithm, and distributed data.
c PCOUNT2YL counts particles by y grid, accumulating into npic array,
c           with distributed data.
c PRMOVE2 remove particles instead of reflecting at boundary, with
c         distributed data.
c PPUSH2ZF update particle co-ordinates for particles with fixed
c          velocities, and distributed data.
c PGCJPOST2 deposits time-centered particle current density, quadratic
c           interpolation, and distributed data.
c PGCJPOST2L deposits time-centered particle current density, linear
c            interpolation, and distributed data.
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: august 22, 2009
c-----------------------------------------------------------------------
      subroutine PDOST2(part,q,npp,noff,qm,nx,idimp,npmax,nblok,nxv,nypm
     1x)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries
c and distributed data.
c baseline scalar distributed version
c 43 flops/particle, 11 loads, 9 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c nx = system length in x direction
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
      dimension part(idimp,npmax,nblok), q(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qmh = .5*qm
      do 20 l = 1, nblok
      mnoff = noff(l) - 2
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn
      if (nl.lt.1) nl = nl + nx
      amx = qm*(.75 - dxp*dxp)
      mm = mm - mnoff
      amy = .75 - dyp*dyp
      nn = nn + 1
      if (nn.gt.nx) nn = nn - nx
      dxl = qmh*(.5 - dxp)**2
      np = nn + 1
      if (np.gt.nx) np = np - nx
      dxp = qmh*(.5 + dxp)**2
      ml = mm - 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c deposit charge
      q(nl,mm,l) = q(nl,mm,l) + dxl*amy
      q(nn,mm,l) = q(nn,mm,l) + amx*amy
      q(np,mm,l) = q(np,mm,l) + dxp*amy
      q(nl,ml,l) = q(nl,ml,l) + dxl*dyl
      q(nn,ml,l) = q(nn,ml,l) + amx*dyl
      q(np,ml,l) = q(np,ml,l) + dxp*dyl
      q(nl,mp,l) = q(nl,mp,l) + dxl*dyp
      q(nn,mp,l) = q(nn,mp,l) + amx*dyp
      q(np,mp,l) = q(np,mp,l) + dxp*dyp
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGPOST2(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,nypmx)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries
c and distributed data.
c scalar version using guard cells, for distributed data
c 43 flops/particle, 11 loads, 9 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j+1,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
      dimension part(idimp,npmax,nblok), q(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
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
c deposit charge
      q(nl,mm,l) = q(nl,mm,l) + dxl*amy
      q(nn,mm,l) = q(nn,mm,l) + amx*amy
      q(np,mm,l) = q(np,mm,l) + dxp*amy
      q(nl,ml,l) = q(nl,ml,l) + dxl*dyl
      q(nn,ml,l) = q(nn,ml,l) + amx*dyl
      q(np,ml,l) = q(np,ml,l) + dxp*dyl
      q(nl,mp,l) = q(nl,mp,l) + dxl*dyp
      q(nn,mp,l) = q(nn,mp,l) + amx*dyp
      q(np,mp,l) = q(np,mp,l) + dxp*dyp
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSPOST2(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,nxyp)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries
c and distributed data.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 43 flops/particle, 11 loads, 9 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j+1,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx+3
c nxyp = first actual dimension of charge array, must be >= nxv*nypmx
      dimension part(idimp,npmax,nblok), q(nxyp,nblok)
      dimension npp(nblok), noff(nblok)
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
c deposit charge
      dx = q(mn,l) + dxl*amy
      dy = q(mn+1,l) + amx*amy
      amy = q(mn+2,l) + dxp*amy
      dx1 = q(ml,l) + dxl*dyl
      dy1 = q(ml+1,l) + amx*dyl
      dyl = q(ml+2,l) + dxp*dyl
      dxl = q(mp,l) + dxl*dyp
      amx = q(mp+1,l) + amx*dyp
      dyp = q(mp+2,l) + dxp*dyp
      q(mn,l) = dx
      q(mn+1,l) = dy
      q(mn+2,l) = amy
      q(ml,l) = dx1
      q(ml+1,l) = dy1
      q(ml+2,l) = dyl
      q(mp,l) = dxl
      q(mp+1,l) = amx
      q(mp+2,l) = dyp
   10 continue
c deposit charge for last particle
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
c deposit charge
      q(mn,l) = q(mn,l) + dxl*amy
      q(mn+1,l) = q(mn+1,l) + amx*amy
      q(mn+2,l) = q(mn+2,l) + dxp*amy
      q(ml,l) = q(ml,l) + dxl*dyl
      q(ml+1,l) = q(ml+1,l) + amx*dyl
      q(ml+2,l) = q(ml+2,l) + dxp*dyl
      q(mp,l) = q(mp,l) + dxl*dyp
      q(mp+1,l) = q(mp+1,l) + amx*dyp
      q(mp+2,l) = q(mp+2,l) + dxp*dyp
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSOST2X(part,q,npp,noff,nn,amxy,qm,nx,idimp,npmax,nblok
     1,nxv,nxvyp,npd,nine)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized distributed version with 1d addressing
c 43 flops/particle, 29 loads, 27 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized charge deposition
c amxy = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c nx = system length in x direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx
c nxvyp = nxv*nypmx, first actual dimension of charge array
c npd = size of scratch buffers for vectorized charge deposition
c nine = number of independent weights
      dimension part(idimp,npmax,nblok), q(nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(nine,npd,nblok), amxy(nine,npd,nblok)
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l) - 1
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l) + .5
      dxl = part(1,i+jb,l) - float(n)
      n = n + 1
      if (n.gt.nx) n = n - nx
      amx = qm*(.75 - dxl*dxl)
      np = n + 1
      if (np.gt.nx) np = np - nx
      dxp = .5*qm*(.5 + dxl)**2
      nl = n - 1
      if (nl.lt.1) nl = nl + nx
      dxl = .5*qm*(.5 - dxl)**2
      m = part(2,i+jb,l) + .5
      dyl = part(2,i+jb,l) - float(m)
      m = nxv*(m - mnoff)
      amy = .75 - dyl*dyl
      mp = m + nxv
      dyp = .5*(.5 + dyl)**2
      ml = m - nxv
      dyl = .5*(.5 - dyl)**2
      nn(1,i,l) = n + m
      nn(2,i,l) = np + m
      nn(3,i,l) = nl + m
      nn(4,i,l) = n + mp
      nn(5,i,l) = np + mp
      nn(6,i,l) = nl + mp
      nn(7,i,l) = n + ml
      nn(8,i,l) = np + ml
      nn(9,i,l) = nl + ml
      amxy(1,i,l) = amx*amy
      amxy(2,i,l) = dxp*amy
      amxy(3,i,l) = dxl*amy
      amxy(4,i,l) = amx*dyp
      amxy(5,i,l) = dxp*dyp
      amxy(6,i,l) = dxl*dyp
      amxy(7,i,l) = amx*dyl
      amxy(8,i,l) = dxp*dyl
      amxy(9,i,l) = dxl*dyl
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 9
      q(nn(k,i,l),l) = q(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSOST2X(part,q,npp,noff,nn,amxy,qm,idimp,npmax,nblok,n
     1xv,nxvyp,npd,nine)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized version with guard cells and 1d addressing
c 43 flops/particle, 29 loads, 27 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j+1,k+1,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized charge deposition
c amxy = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx
c nxvyp = nxv*nypmx, first actual dimension of charge array
c npd = size of scratch buffers for vectorized charge deposition
c nine = number of independent weights
      dimension part(idimp,npmax,nblok), q(nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(nine,npd,nblok), amxy(nine,npd,nblok)
      qmh = .5*qm
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l)
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l) + .5
      m = part(2,i+jb,l) + .5
      dxp = part(1,i+jb,l) - float(n)
      dyp = part(2,i+jb,l) - float(m)
      n = n + 1
      m = nxv*(m - mnoff)
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      ml = m + n
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mp = mn + nxv
      nn(1,i,l) = mn
      nn(2,i,l) = mn + 1
      nn(3,i,l) = mn + 2
      nn(4,i,l) = ml
      nn(5,i,l) = ml + 1
      nn(6,i,l) = ml + 2
      nn(7,i,l) = mp
      nn(8,i,l) = mp + 1
      nn(9,i,l) = mp + 2
      amxy(1,i,l) = dxl*amy
      amxy(2,i,l) = amx*amy
      amxy(3,i,l) = dxp*amy
      amxy(4,i,l) = dxl*dyl
      amxy(5,i,l) = amx*dyl
      amxy(6,i,l) = dxp*dyl
      amxy(7,i,l) = dxl*dyp
      amxy(8,i,l) = amx*dyp
      amxy(9,i,l) = dxp*dyp
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 9
      q(nn(k,i,l),l) = q(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDOST2L(part,q,npp,noff,qm,nx,idimp,npmax,nblok,nxv,nyp
     1mx)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c and distributed data.
c baseline scalar distributed version
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c nx = system length in x direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
      dimension part(idimp,npmax,nblok), q(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
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
      np = nn + 1
      if (np.gt.nx) np = np - nx
      amy = 1. - dyp
      mp = mm + 1
c deposit charge
      q(nn,mm,l) = q(nn,mm,l) + amx*amy
      q(np,mm,l) = q(np,mm,l) + dxp*amy
      q(nn,mp,l) = q(nn,mp,l) + amx*dyp
      q(np,mp,l) = q(np,mp,l) + dxp*dyp
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGPOST2L(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,nypmx
     1)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c and distributed data.
c scalar version using guard cells, for distributed data
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
      dimension part(idimp,npmax,nblok), q(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
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
c deposit charge
      q(np,mp,l) = q(np,mp,l) + dxp*dyp
      q(nn,mp,l) = q(nn,mp,l) + amx*dyp
      q(np,mm,l) = q(np,mm,l) + dxp*amy
      q(nn,mm,l) = q(nn,mm,l) + amx*amy
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSPOST2L(part,q,npp,noff,qm,idimp,npmax,nblok,nxv,nxyp
     1)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c and distributed data.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx+1
c nxyp = actual first dimension of charge array, must be >= nxv*nypmx
      dimension part(idimp,npmax,nblok), q(nxyp,nblok)
      dimension npp(nblok), noff(nblok)
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
c deposit charge
      dx1 = q(mp+1,l) + dxp*dyp
      dyp = q(mp,l) + amx*dyp
      dxp = q(mm+1,l) + dxp*amy
      amy = q(mm,l) + amx*amy
      q(mp+1,l) = dx1
      q(mp,l) = dyp
      q(mm+1,l) = dxp
      q(mm,l) = amy
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit charge
      q(mp+1,l) = q(mp+1,l) + dxp*dyn
      q(mp,l) = q(mp,l) + amx*dyn
      q(mm+1,l) = q(mm+1,l) + dxp*amy
      q(mm,l) = q(mm,l) + amx*amy
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSOST2XL(part,q,npp,noff,nn,amxy,qm,nx,idimp,npmax,nblo
     1k,nxv,nxvyp,npd,ifour)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized distributed version with 1d addressing
c 17 flops/particle, 14 loads, 12 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized charge deposition
c amxy = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c nx = system length in x direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx
c nxvyp = nxv*nypmx, first actual dimension of charge array
c npd = size of scratch buffers for vectorized charge deposition
c ifour = number of independent weights
      dimension part(idimp,npmax,nblok), q(nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(ifour,npd,nblok), amxy(ifour,npd,nblok)
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l)
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l)
      dxp = qm*(part(1,i+jb,l) - float(n))
      n = n + 1
      amx = qm - dxp
      np = n + 1
      if (np.gt.nx) np = np - nx
      m = part(2,i+jb,l)
      dyp = part(2,i+jb,l) - float(m)
      m = nxv*(m - mnoff)
      amy = 1. - dyp
      mp = m + nxv
      nn(1,i,l) = n + m
      nn(2,i,l) = np + m
      nn(3,i,l) = n + mp
      nn(4,i,l) = np + mp
      amxy(1,i,l) = amx*amy
      amxy(2,i,l) = dxp*amy
      amxy(3,i,l) = amx*dyp
      amxy(4,i,l) = dxp*dyp
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 4
      q(nn(k,i,l),l) = q(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSOST2XL(part,q,npp,noff,nn,amxy,qm,idimp,npmax,nblok,
     1nxv,nxvyp,npd,ifour)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized distributed version with guard cells and 1d addressing
c 17 flops/particle, 14 loads, 12 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c q(j,k,l) = charge density at grid point (j,kk),
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c nn = scratch address array for vectorized charge deposition
c amxy = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx+1
c nxvyp = nxv*nypmx, first actual dimension of charge array
c npd = size of scratch buffers for vectorized charge deposition
c ifour = number of independent weights
      dimension part(idimp,npmax,nblok), q(nxvyp,nblok)
      dimension npp(nblok), noff(nblok)
      dimension nn(ifour,npd,nblok), amxy(ifour,npd,nblok)
c parallel loop
      do 50 l = 1, nblok
      mnoff = noff(l)
      npb = npd
      if (npp(l).gt.npd) then
         ipp = float(npp(l) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(l) - (ipp - 1)*npd
      do 10 i = 1, npb
c find interpolation weights
      n = part(1,i+jb,l)
      m = part(2,i+jb,l)
      dxp = qm*(part(1,i+jb,l) - float(n))
      dyp = part(2,i+jb,l) - float(m)
      n = n + 1
      mm = nxv*(m - mnoff)
      amx = qm - dxp
      mm = mm + n
      amy = 1. - dyp
      mp = mm + nxv
      nn(4,i,l) = mm
      nn(3,i,l) = mm + 1
      nn(2,i,l) = mp
      nn(1,i,l) = mp + 1
      amxy(1,i,l) = dxp*dyp
      amxy(2,i,l) = amx*dyp
      amxy(3,i,l) = dxp*amy
      amxy(4,i,l) = amx*amy
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 4
      q(nn(k,i,l),l) = q(nn(k,i,l),l) + amxy(k,i,l)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPUSH2(part,fx,fy,npp,noff,qbm,dt,ek,nx,idimp,npmax,nbl
     1ok,nxv,nypmx)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with periodic boundary conditions,
c for distributed data.
c baseline scalar distributed version 
c 80 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fx(j,k,l) = x component of force/charge at grid (j,kk)
c fy(j,k,l) = y component of force/charge at grid (j,kk)
c in other words, fx/fy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx = system length in x direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      zero = 0.
      anx = float(nx)
      qtm = qbm*dt
      sum1 = 0.0d0
      do 20 l = 1, nblok
      mnoff = noff(l) - 2
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l) + .5
      mm = part(2,j,l) + .5
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nl = nn
      if (nl.lt.1) nl = nl + nx
      mm = mm - mnoff
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      nn = nn + 1
      if (nn.gt.nx) nn = nn - nx
      dxl = .5*(.5 - dxp)**2
      np = nn + 1
      if (np.gt.nx) np = np - nx
      dxp = .5*(.5 + dxp)**2
      ml = mm - 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c find acceleration
      dx = amy*(dxl*fx(nl,mm,l) + amx*fx(nn,mm,l) + dxp*fx(np,mm,l)) + d
     1yl*(dxl*fx(nl,ml,l) + amx*fx(nn,ml,l) + dxp*fx(np,ml,l)) + dyp*(dx
     2l*fx(nl,mp,l) + amx*fx(nn,mp,l) + dxp*fx(np,mp,l))
      dy = amy*(dxl*fy(nl,mm,l) + amx*fy(nn,mm,l) + dxp*fy(np,mm,l)) + d
     1yl*(dxl*fy(nl,ml,l) + amx*fy(nn,ml,l) + dxp*fy(np,ml,l)) + dyp*(dx
     2l*fy(nl,mp,l) + amx*fy(nn,mp,l) + dxp*fy(np,mp,l))
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGPUSH2(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,n
     1blok,nxv,nypmx,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with various boundary conditions,
c for distributed data.
c scalar version using guard cells, for distributed data
c 80 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
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
      mp = ml + 2
      dyp = .5*(.5 + dyp)**2
c find acceleration
      dx = amy*(dxl*fxy(1,nl,mm,l) + amx*fxy(1,nn,mm,l) + dxp*fxy(1,np,m
     1m,l)) + dyl*(dxl*fxy(1,nl,ml,l) + amx*fxy(1,nn,ml,l) + dxp*fxy(1,n
     2p,ml,l)) + dyp*(dxl*fxy(1,nl,mp,l) + amx*fxy(1,nn,mp,l) + dxp*fxy(
     31,np,mp,l))
      dy = amy*(dxl*fxy(2,nl,mm,l) + amx*fxy(2,nn,mm,l) + dxp*fxy(2,np,m
     1m,l)) + dyl*(dxl*fxy(2,nl,ml,l) + amx*fxy(2,nn,ml,l) + dxp*fxy(2,n
     2p,ml,l)) + dyp*(dxl*fxy(2,nl,mp,l) + amx*fxy(2,nn,mp,l) + dxp*fxy(
     32,np,mp,l))
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSPUSH2(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,
     1nblok,nxv,nxyp,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with various boundary conditions,
c for distributed data.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 80 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = second virtual dimension of field array, must be >= nx+3
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
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
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1,l) - float(nnn)
      dyn = part(2,j+1,l) - float(mmn)
      ml = mm + nn
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      mn = ml + nxv
      dxl = .5*(.5 - dxp)**2
      dxp = .5*(.5 + dxp)**2
      mp = mn + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mmn = mmn - mnoff
c find acceleration
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1+ dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) + 
     2dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l))
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1+ dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) + 
     2dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l))
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
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
c find acceleration
      dx = amy*(dxl*fxy(1,mn,l) + amx*fxy(1,mn+1,l) + dxp*fxy(1,mn+2,l))
     1+ dyl*(dxl*fxy(1,ml,l) + amx*fxy(1,ml+1,l) + dxp*fxy(1,ml+2,l)) + 
     2dyp*(dxl*fxy(1,mp,l) + amx*fxy(1,mp+1,l) + dxp*fxy(1,mp+2,l))
      dy = amy*(dxl*fxy(2,mn,l) + amx*fxy(2,mn+1,l) + dxp*fxy(2,mn+2,l))
     1+ dyl*(dxl*fxy(2,ml,l) + amx*fxy(2,ml+1,l) + dxp*fxy(2,ml+2,l)) + 
     2dyp*(dxl*fxy(2,mp,l) + amx*fxy(2,mp+1,l) + dxp*fxy(2,mp+2,l))
c new velocity
      dx = part(3,nop,l) + qtm*dx
      dy = part(4,nop,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,nop,l))**2 + (dy + part(4,nop,l))**2
      part(3,nop,l) = dx
      part(4,nop,l) = dy
c new position
      dx = part(1,nop,l) + dx*dt
      dy = part(2,nop,l) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PPUSH2L(part,fx,fy,npp,noff,qbm,dt,ek,nx,idimp,npmax,nb
     1lok,nxv,nypmx)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with periodic boundary conditions,
c for distributed data.
c baseline scalar distributed version
c 42 flops/particle, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fx(j,k,l) = x component of force/charge at grid (j,kk)
c fy(j,k,l) = y component of force/charge at grid (j,kk)
c in other words, fx/fy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx = system length in x direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      zero = 0.
      anx = float(nx)
      qtm = qbm*dt
      sum1 = 0.0d0
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
      np = nn + 1
      if (np.gt.nx) np = np - nx
      amy = 1. - dyp
      mp = mm + 1
c find acceleration
      dx = amy*(amx*fx(nn,mm,l) + dxp*fx(np,mm,l)) + dyp*(amx*fx(nn,mp,l
     1) + dxp*fx(np,mp,l))
      dy = amy*(amx*fy(nn,mm,l) + dxp*fy(np,mm,l)) + dyp*(amx*fy(nn,mp,l
     1) + dxp*fy(np,mp,l))
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGPUSH2L(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,
     1nblok,nxv,nypmx,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions,
c for distributed data.
c scalar version using guard cells, for distributed data
c 42 flops/particle, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
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
c find acceleration
      dx = dyp*(dxp*fxy(1,np,mp,l) + amx*fxy(1,nn,mp,l)) + amy*(dxp*fxy(
     11,np,mm,l) + amx*fxy(1,nn,mm,l))
      dy = dyp*(dxp*fxy(2,np,mp,l) + amx*fxy(2,nn,mp,l)) + amy*(dxp*fxy(
     12,np,mm,l) + amx*fxy(2,nn,mm,l))
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSPUSH2L(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax
     1,nblok,nxv,nxyp,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions,
c for distributed data.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996)
c 42 flops/particle, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c in other words, fxy are the convolutions of the electric field
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = second virtual dimension of field array, must be >= nx+1
c nxyp = second actual dimension of field array, must be >= nxv*nypmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok)
      dimension npp(nblok), noff(nblok)
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
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
c find acceleration
      dx = dyp*(dxp*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxp*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyp*(dxp*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxp*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
c new velocity
      dx = part(3,j,l) + qtm*dx
      dy = part(4,j,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j,l))**2 + (dy + part(4,j,l))**2
      part(3,j,l) = dx
      part(4,j,l) = dy
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
      nop = npp(l)
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find acceleration
      dx = dyn*(dxn*fxy(1,mp+1,l) + amx*fxy(1,mp,l)) + amy*(dxn*fxy(1,mm
     1+1,l) + amx*fxy(1,mm,l))
      dy = dyn*(dxn*fxy(2,mp+1,l) + amx*fxy(2,mp,l)) + amy*(dxn*fxy(2,mm
     1+1,l) + amx*fxy(2,mm,l))
c new velocity
      dx = part(3,nop,l) + qtm*dx
      dy = part(4,nop,l) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,nop,l))**2 + (dy + part(4,nop,l))**2
      part(3,nop,l) = dx
      part(4,nop,l) = dy
c new position
      dx = part(1,nop,l) + dx*dt
      dy = part(2,nop,l) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,l)
            part(4,nop,l) = -part(4,nop,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,l)
            part(3,nop,l) = -part(3,nop,l)
         endif
      endif
c set new position
      part(1,nop,l) = dx
      part(2,nop,l) = dy
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PSORTP2Y(part,pt,ip,npic,npp,noff,nyp,idimp,npmax,nblok
     1,nypm1)
c this subroutine sorts particles by y grid
c quadratic interpolation, spatial decomposition in y direction
c part = particle array
c part(2,n,m) = position y of particle n in partition m
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c npp(m) = number of particles in partition m
c noff(m) = backmost global gridpoint in particle partition m
c nyp(m) = number of primary gridpoints in particle partition m
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nypm1 = maximum size of particle partition plus one
      dimension part(idimp,npmax,nblok), pt(npmax,nblok)
      dimension ip(npmax,nblok), npic(nypm1,nblok)
      dimension npp(nblok), noff(nblok), nyp(nblok)
      do 80 l = 1, nblok
      mnoff = noff(l) - 1
      nyp1 = nyp(l) + 1
c clear counter array
      do 10 k = 1, nyp1
      npic(k,l) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp(l)
      m = part(2,j,l) + .5
      m = m - mnoff
      npic(m,l) = npic(m,l) + 1
      ip(j,l) = m
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyp1
      ist = npic(k,l)
      npic(k,l) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, npp(l)
      m = ip(j,l)
      npic(m,l) = npic(m,l) + 1
      ip(j,l) = npic(m,l)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, npp(l)
      pt(ip(j,l),l) = part(i,j,l)
   50 continue
      do 60 j = 1, npp(l)
      part(i,j,l) = pt(j,l)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSORTP2YL(part,pt,ip,npic,npp,noff,nyp,idimp,npmax,nblo
     1k,nypm1)
c this subroutine sorts particles by y grid
c linear interpolation, spatial decomposition in y direction
c part = particle array
c part(2,n,m) = position y of particle n in partition m
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c npp(m) = number of particles in partition m
c noff(m) = backmost global gridpoint in particle partition m
c nyp(m) = number of primary gridpoints in particle partition m
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nypm1 = maximum size of particle partition plus one
      dimension part(idimp,npmax,nblok), pt(npmax,nblok)
      dimension ip(npmax,nblok), npic(nypm1,nblok)
      dimension npp(nblok), noff(nblok), nyp(nblok)
      do 80 l = 1, nblok
      mnoff = noff(l) - 1
      nyp1 = nyp(l) + 1
c clear counter array
      do 10 k = 1, nyp1
      npic(k,l) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp(l)
      m = part(2,j,l)
      m = m - mnoff
      npic(m,l) = npic(m,l) + 1
      ip(j,l) = m
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyp1
      ist = npic(k,l)
      npic(k,l) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, npp(l)
      m = ip(j,l)
      npic(m,l) = npic(m,l) + 1
      ip(j,l) = npic(m,l)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, npp(l)
      pt(ip(j,l),l) = part(i,j,l)
   50 continue
      do 60 j = 1, npp(l)
      part(i,j,l) = pt(j,l)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDSORTP2Y(parta,partb,npic,npp,noff,nyp,idimp,npmax,nbl
     1ok,nypm1)
c this subroutine sorts particles by y grid
c quadratic interpolation, spatial decomposition in y direction
c parta/partb = input/output particle array
c part(2,n,m) = position y of particle n in partition m
c npic = address offset for reordering particles
c npp(m) = number of particles in partition m
c noff(m) = backmost global gridpoint in particle partition m
c nyp(m) = number of primary gridpoints in particle partition m
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nypm1 = maximum size of particle partition plus one
      implicit none
      integer idimp, npmax, nblok, nypm1
      integer npic, npp, noff, nyp
      real parta, partb
      dimension parta(idimp,npmax,nblok), partb(idimp,npmax,nblok)
      dimension npic(nypm1,nblok)
      dimension npp(nblok), noff(nblok), nyp(nblok)
c local data
      integer i, j, k, l, m, mnoff, nyp1, isum, ist, ip
      do 60 l = 1, nblok
      mnoff = noff(l) - 1
      nyp1 = nyp(l) + 1
c clear counter array
      do 10 k = 1, nyp1
      npic(k,l) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp(l)
      m = parta(2,j,l) + .5
      m = m - mnoff
      npic(m,l) = npic(m,l) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyp1
      ist = npic(k,l)
      npic(k,l) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, npp(l)
      m = parta(2,j,l) + .5
      m = m - mnoff
      ip = npic(m,l) + 1
      do 40 i = 1, idimp
      partb(i,ip,l) = parta(i,j,l)
   40 continue
      npic(m,l) = ip
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDSORTP2YL(parta,partb,npic,npp,noff,nyp,idimp,npmax,nb
     1lok,nypm1)
c this subroutine sorts particles by y grid
c linear interpolation, spatial decomposition in y direction
c parta/partb = input/outpu particle array
c part(2,n,m) = position y of particle n in partition m
c npic = address offset for reordering particles
c npp(m) = number of particles in partition m
c noff(m) = backmost global gridpoint in particle partition m
c nyp(m) = number of primary gridpoints in particle partition m
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nypm1 = maximum size of particle partition plus one
      implicit none
      integer idimp, npmax, nblok, nypm1
      integer npic, npp, noff, nyp
      real parta, partb
      dimension parta(idimp,npmax,nblok), partb(idimp,npmax,nblok)
      dimension npic(nypm1,nblok)
      dimension npp(nblok), noff(nblok), nyp(nblok)
c local data
      integer i, j, k, l, m, mnoff, nyp1, isum, ist, ip
      do 60 l = 1, nblok
      mnoff = noff(l) - 1
      nyp1 = nyp(l) + 1
c clear counter array
      do 10 k = 1, nyp1
      npic(k,l) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp(l)
      m = parta(2,j,l)
      m = m - mnoff
      npic(m,l) = npic(m,l) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyp1
      ist = npic(k,l)
      npic(k,l) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, npp(l)
      m = parta(2,j,l)
      m = m - mnoff
      ip = npic(m,l) + 1
      do 40 i = 1, idimp
      partb(i,ip,l) = parta(i,j,l)
   40 continue
      npic(m,l) = ip
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCOUNT2YL(part,isign,npic,npp,noff,nyp,idimp,npmax,nblo
     1k,nypm1)
c this subroutine counts particles by y grid, accumulating into npic
c npic array is initialized to zero, unless isign = 0
c part = particle array
c isign = (0,1) = (no,yes) initialize npic array to zero
c part(2,j,l) = position y of particle j in partition l
c npic = number of particles per grid
c npp(l) = number of particles in partition l
c noff(l) = backmost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nypm1 = maximum size of particle partition plus one
      implicit none
      real part
      integer npic, npp, noff, nyp
      integer isign, idimp, npmax, nblok, nypm1
      dimension part(idimp,npmax,nblok)
      dimension npic(nypm1,nblok), npp(nblok), noff(nblok), nyp(nblok)
      integer j, l, m, mnoff, nyp1
      do 30 l = 1, nblok
      mnoff = noff(l) - 1
      nyp1 = nyp(l) + 1
      if (isign.ne.0) then
         do 10 j = 1, nyp1
         npic(j,l) = 0
   10    continue
      endif
c find how many particles in each grid
      do 20 j = 1, npp(l)
      m = part(2,j,l)
      m = m - mnoff
      npic(m,l) = npic(m,l) + 1
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PRMOVE2(part,npp,ihole,jss,nx,ny,idimp,npmax,nblok,idps
     1,ntmax,ipbc)
c this subroutine removes particles which would normally be reflected
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c npp(l) = number of particles in partition l
c ihole = location of holes left in particle arrays
c jss(idps,l) = scratch array for particle partition l
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,-1,-2,-3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      real part
      integer npp, ihole, jss
      integer nx, ny, idimp, npmax, nblok, idps, ntmax, ipbc
      dimension part(idimp,npmax,nblok)
      dimension npp(nblok), ihole(ntmax,nblok), jss(idps,nblok)
c local data
      integer i, j, j1, j2, l, nter
      real edgelx, edgely, edgerx, edgery, dx, dy
c set boundary values
      if (ipbc.eq.(-1)) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.(-2)) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.(-3)) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      nter = 0
c buffer outgoing particles
      do 60 l = 1, nblok
   10 jss(1,l) = 0
      jss(2,l) = 0
      do 20 j = 1, npp(l)
      dx = part(1,j,l)
      dy = part(2,j,l)
c periodic boundary conditions
      if (ipbc.eq.(-1)) then
         if (dx.lt.edgelx) part(1,j,l) = dx + edgerx
         if (dx.ge.edgerx) part(1,j,l) = part(1,j,l) - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.(-2)) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            if (jss(1,l).lt.ntmax) then
               jss(1,l) = jss(1,l) + 1
               ihole(jss(1,l),l) = j
            else
               jss(2,l) = 1
               go to 30
            endif
         else if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            if (jss(1,l).lt.ntmax) then
               jss(1,l) = jss(1,l) + 1
               ihole(jss(1,l),l) = j
            else
               jss(2,l) = 1
               go to 30
            endif
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.(-3)) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            if (jss(1,l).lt.ntmax) then
               jss(1,l) = jss(1,l) + 1
               ihole(jss(1,l),l) = j
            else
               jss(2,l) = 1
               go to 30
            endif
         endif
      endif
   20 continue
   30 continue
c fill up holes in particle array with particles from bottom
      do 50 j = 1, jss(1,l)
      j1 = npp(l) - j + 1
      j2 = jss(1,l) - j + 1
      if (j1.gt.ihole(j2,l)) then
c move particle only if it is below current hole
         do 40 i = 1, idimp
         part(i,ihole(j2,l),l) = part(i,j1,l)
   40    continue
      endif
   50 continue
      npp(l) = npp(l) - jss(1,l)
c check if buffer overflowed and more particles remain to be checked
      if (jss(2,l).gt.0) then
         nter = nter + 1
         go to 10
      endif
      jss(2,l) = nter
   60 continue
c information
      nter = 0
      do 70 l = 1, nblok 
      nter = max0(nter,jss(2,l))
   70 continue
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, ntmax=', ntmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPUSH2ZF(part,npp,dt,ek,idimp,npmax,nblok,nx,ny,ipbc)
c for 2d code, this subroutine updates particle co-ordinates for
c particles with fixed velocities, with various boundary conditions,
c for distributed data
c 7 flops/particle, 4 loads, 2 stores
c input: all, output: part, ek
c equations used are:
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c npp(l) = number of particles in partition l
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t+dt/2)**2+vy(t+dt/2)**2)
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npp, idimp, npmax, nblok, nx, ny, ipbc
      real part, dt, ek
      dimension part(idimp,npmax,nblok)
      dimension npp(nblok)
      integer j, l
      real edgelx, edgely, edgerx, edgery, dx, dy
      double precision sum1
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgerx = float(nx-1)
      endif
      do 20 l = 1, nblok
      do 10 j = 1, npp(l)
c velocity
      dx = part(3,j,l)
      dy = part(4,j,l)
c average kinetic energy
      sum1 = sum1 + dx**2 + dy**2
c new position
      dx = part(1,j,l) + dx*dt
      dy = part(2,j,l) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,l)
            part(4,j,l) = -part(4,j,l)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,l)
            part(3,j,l) = -part(3,j,l)
         endif
      endif
c set new position
      part(1,j,l) = dx
      part(2,j,l) = dy
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGCJPOST2(part,fxy,npp,noff,cu,qm,qbm,dt,idimp,npmax,nb
     1lok,nxv,nypmx)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation.
c scalar version using guard cells, for distributed data
c 96 flops/particle, 40 loads, 18 stores
c input: all, output: cu
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
c and qci = qm*vj, where j = x,y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = x velocity of particle n at t - dt/2 in partition l
c part(4,n,l) = y velocity of particle n at t - dt/2 in partition l
c fxy(1,j+1,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j+1,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c over the particle shape, where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j+1,k,l) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current density array, must be >= nx+3
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp, noff, idimp, npmax, nblok, nxv, nypmx
      real part, fxy, cu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, nl, np, ml, mp
      real qtmh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy
      qtmh = 0.5*qbm*dt
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
c new velocity
      vx = part(3,j,l) + qtmh*dx
      vy = part(4,j,l) + qtmh*dy
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      cu(1,nl,mm,l) = cu(1,nl,mm,l) + vx*dx
      cu(2,nl,mm,l) = cu(2,nl,mm,l) + vy*dx
      dx = dxl*dyl
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + vx*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + vy*dy
      dy = amx*dyl
      cu(1,np,mm,l) = cu(1,np,mm,l) + vx*dz
      cu(2,np,mm,l) = cu(2,np,mm,l) + vy*dz
      dz = dxp*dyl
      cu(1,nl,ml,l) = cu(1,nl,ml,l) + vx*dx
      cu(2,nl,ml,l) = cu(2,nl,ml,l) + vy*dx
      dx = dxl*dyp
      cu(1,nn,ml,l) = cu(1,nn,ml,l) + vx*dy
      cu(2,nn,ml,l) = cu(2,nn,ml,l) + vy*dy
      dy = amx*dyp
      cu(1,np,ml,l) = cu(1,np,ml,l) + vx*dz
      cu(2,np,ml,l) = cu(2,np,ml,l) + vy*dz
      dz = dxp*dyp
      cu(1,nl,mp,l) = cu(1,nl,mp,l) + vx*dx
      cu(2,nl,mp,l) = cu(2,nl,mp,l) + vy*dx
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + vx*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + vy*dy
      cu(1,np,mp,l) = cu(1,np,mp,l) + vx*dz
      cu(2,np,mp,l) = cu(2,np,mp,l) + vy*dz
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGCJPOST2L(part,fxy,npp,noff,cu,qm,qbm,dt,idimp,npmax,n
     1blok,nxv,nypmx)
c for 2d code, this subroutine calculates particle current density
c using first-order spline interpolation.
c scalar version using guard cells, for distributed data
c 52 flops/particle, 20 loads, 0 stores
c input: all, output: cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1.m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n,l) = position x of particle n at t in partition l
c part(2,n,l) = position y of particle n at t in partition l
c part(3,n,l) = velocity vx of particle n at t - dt/2 in partition l
c part(4,n,l) = velocity vy of particle n at t - dt/2 in partition l
c fxy(1,j,k,l) = x component of force/charge at grid (j,kk)
c fxy(2,j,k,l) = y component of force/charge at grid (j,kk)
c that is, convolution of electric field over particle shape
c where kk = k + noff(l) - 1
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l.
c cu(i,j,k,l) = ith component of current density
c at grid point j,kk for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c nxv = first dimension of current density array, must be >= nx+1
c nypmx = maximum size of particle partition, including guard cells.
      implicit none
      integer npp,noff, idimp, npmax, nblok, nxv, nypmx
      real part, fxy, cu, qm, qbm, dt
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), noff(nblok)
      integer j, l, mnoff, nn, mm, np, mp
      real qtmh, dxp, dyp, amx, amy, dx, dy, vx, vy
      qtmh = 0.5*qbm*dt
      do 20 l = 1, nblok
      mnoff = noff(l) - 1
      do 10 j = 1, npp(l)
c find interpolation weights
      nn = part(1,j,l)
      mm = part(2,j,l)
      dxp = part(1,j,l) - float(nn)
      dyp = part(2,j,l) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp,l) + amx*fxy(1,nn,mp,l)) + amy*(dxp*fxy(
     11,np,mm,l) + amx*fxy(1,nn,mm,l))
      dy = dyp*(dxp*fxy(2,np,mp,l) + amx*fxy(2,nn,mp,l)) + amy*(dxp*fxy(
     12,np,mm,l) + amx*fxy(2,nn,mm,l))
c new velocity
      vx = part(3,j,l) + qtmh*dx
      vy = part(4,j,l) + qtmh*dy
c deposit current density
      amx = qm*amx
      dxp = qm*dxp
      dx = dxp*dyp
      dy = amx*dyp
      cu(1,np,mp,l) = cu(1,np,mp,l) + vx*dx
      cu(2,np,mp,l) = cu(2,np,mp,l) + vy*dx
      dx = dxp*amy
      cu(1,nn,mp,l) = cu(1,nn,mp,l) + vx*dy
      cu(2,nn,mp,l) = cu(2,nn,mp,l) + vy*dy
      dy = amx*amy
      cu(1,np,mm,l) = cu(1,np,mm,l) + vx*dx
      cu(2,np,mm,l) = cu(2,np,mm,l) + vy*dx
      cu(1,nn,mm,l) = cu(1,nn,mm,l) + vx*dy
      cu(2,nn,mm,l) = cu(2,nn,mm,l) + vy*dy
   10 continue
   20 continue
      return
      end
