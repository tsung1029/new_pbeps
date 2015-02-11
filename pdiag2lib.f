c-----------------------------------------------------------------------
c 2d parallel PIC library for diagnostics
c pdiag2lib.f contains diagnostic procedures:
c PVDIST2 calculates 2 component velocity distribution and velocity
c         moments, for distributed data.
c PVDIST23 calculates 3 component velocity distribution and velocity
c         moments, for distributed data.
c PSDIST2 calculates entropy from 2 component velocity distribution,
c         for distributed data.
c PSDIST23 calculates entropy from 3 component velocity distribution,
c          for distributed data.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: november 6, 2009
c-----------------------------------------------------------------------
      subroutine PVDIST2(part,fv,fvm,npp,idimp,npmax,nblok,nmv,nmvf)
c for 2d code, this subroutine calculates 2d velocity distribution
c and velocity moments
c input: all except fvm, output: fv, fvm
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c fv = distribution function, number of particles in each velocity range
c maximum velocity (used for scaling) is contained in first element fv.
c vdrift for i-th dimension is contained in fvm(1,i)
c vth for i-th dimension is contained in fvm(2,i)
c entropy for i-th dimension is contained in fvm(3,i), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
c is uniform in space and distributions in each dimension are
c independent.
c npp(l) = number of particles in partition l
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, npmax, nblok, nmv, nmvf
      integer npp
      real part, fv, fvm
      dimension part(idimp,npmax,nblok)
      dimension fv(nmvf,2,nblok), fvm(3,2,nblok)
      dimension npp(nblok)
c local data
      integer j, l, nvx, nvy
      real anmv, svx, svy
      double precision sumvx, sumvy, sumvx2, sumvy2, anp
      double precision sum5, work5
      dimension sum5(5), work5(5)
      anmv = real(nmv)
      do 40 l = 1, nblok
      svx = anmv/fv(1,1,l)
      svy = anmv/fv(1,2,l)
c zero out distribution
      do 10 j = 2, nmvf
      fv(j,1,l) = 0.
      fv(j,2,l) = 0.
   10 continue
c count particles in each velocity region
      anp = 0.0d0
      anmv = anmv + 2.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      do 20 j = 1, npp(l)
      anp = anp + 1.
      nvx = part(3,j,l)*svx + anmv
      sumvx = sumvx + part(3,j,l)
      sumvx2 = sumvx2 + part(3,j,l)**2
      nvy = part(4,j,l)*svy + anmv
      sumvy = sumvy + part(4,j,l)
      sumvy2 = sumvy2 + part(4,j,l)**2
      if ((nvx.ge.2).and.(nvx.le.nmvf)) fv(nvx,1,l) = fv(nvx,1,l) + 1.
      if ((nvy.ge.2).and.(nvy.le.nmvf)) fv(nvy,2,l) = fv(nvy,2,l) + 1.
   20 continue
      sum5(1) = sumvx
      sum5(2) = sumvy
      sum5(3) = sumvx2
      sum5(4) = sumvy2
      sum5(5) = anp
      call PDSUM(sum5,work5,5,1)
      sumvx = sum5(1)
      sumvy = sum5(2)
      sumvx2 = sum5(3)
      sumvy2 = sum5(4)
      anp = sum5(5)
c calculate velocity moments
      if (anp.ne.0.0d0) anp = 1.0d0/anp
      sumvx = sumvx*anp
      fvm(1,1,l) = sumvx
      fvm(2,1,l) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(1,2,l) = sumvy
      fvm(2,2,l) = dsqrt(sumvy2*anp - sumvy**2)
c calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      do 30 j = 2, nmvf
      if (fv(j,1,l).gt.0.) then
         sumvx = sumvx + fv(j,1,l)
         sumvx2 = sumvx2 + fv(j,1,l)*dlog(dble(fv(j,1,l)*svx))
      endif
      if (fv(j,2,l).gt.0.) then
         sumvy = sumvy + fv(j,2,l)
         sumvy2 = sumvy2 + fv(j,2,l)*dlog(dble(fv(j,2,l)*svy))
      endif
   30 continue
      sum5(1) = sumvx
      sum5(2) = sumvy
      sum5(3) = sumvx2
      sum5(4) = sumvy2
      call PDSUM(sum5,work5,4,1)
      sumvx = sum5(1)
      sumvy = sum5(2)
      sumvx2 = sum5(3)
      sumvy2 = sum5(4)
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      fvm(3,1,l) = sumvx
      fvm(3,2,l) = sumvy
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PVDIST23(part,fv,fvm,npp,idimp,npmax,nblok,nmv,nmvf)
c for 2-1/2d code, this subroutine calculates 3d velocity distribution
c and velocity moments
c input: all except fvm, output: fv, fvm
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c fv = distribution function, number of particles in each velocity range
c maximum velocity (used for scaling) is contained in first element fv.
c vdrift for i-th dimension is contained in fvm(1,i)
c vth for i-th dimension is contained in fvm(2,i)
c entropy for i-th dimension is contained in fvm(3,i), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
c is uniform in space and distributions in each dimension are
c independent.
c npp(l) = number of particles in partition l
c idimp = size of phase space = 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, npmax, nblok, nmv, nmvf
      integer npp
      real part, fv, fvm
      dimension part(idimp,npmax,nblok)
      dimension fv(nmvf,3,nblok), fvm(3,3,nblok)
      dimension npp(nblok)
c local data
      integer j, l, nvx, nvy, nvz
      real anmv, svx, svy, svz
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2, anp
      double precision sum7, work7
      dimension sum7(7), work7(7)
      anmv = real(nmv)
      do 40 l = 1, nblok
      svx = anmv/fv(1,1,l)
      svy = anmv/fv(1,2,l)
      svz = anmv/fv(1,3,l)
c zero out distribution
      do 10 j = 2, nmvf
      fv(j,1,l) = 0.
      fv(j,2,l) = 0.
      fv(j,3,l) = 0.
   10 continue
c count particles in each velocity region
      anp = 0.0d0
      anmv = anmv + 2.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 20 j = 1, npp(l)
      anp = anp + 1.
      nvx = part(3,j,l)*svx + anmv
      sumvx = sumvx + part(3,j,l)
      sumvx2 = sumvx2 + part(3,j,l)**2
      nvy = part(4,j,l)*svy + anmv
      sumvy = sumvy + part(4,j,l)
      sumvy2 = sumvy2 + part(4,j,l)**2
      nvz = part(5,j,l)*svz + anmv
      sumvz = sumvz + part(5,j,l)
      sumvz2 = sumvz2 + part(5,j,l)**2
      if ((nvx.ge.2).and.(nvx.le.nmvf)) fv(nvx,1,l) = fv(nvx,1,l) + 1.
      if ((nvy.ge.2).and.(nvy.le.nmvf)) fv(nvy,2,l) = fv(nvy,2,l) + 1.
      if ((nvz.ge.2).and.(nvz.le.nmvf)) fv(nvz,3,l) = fv(nvz,3,l) + 1.
   20 continue
      sum7(1) = sumvx
      sum7(2) = sumvy
      sum7(3) = sumvz
      sum7(4) = sumvx2
      sum7(5) = sumvy2
      sum7(6) = sumvz2
      sum7(7) = anp
      call PDSUM(sum7,work7,7,1)
      sumvx = sum7(1)
      sumvy = sum7(2)
      sumvz = sum7(3)
      sumvx2 = sum7(4)
      sumvy2 = sum7(5)
      sumvz2 = sum7(6)
      anp = sum7(7)
c calculate velocity moments
      if (anp.ne.0.0d0) anp = 1.0d0/anp
      sumvx = sumvx*anp
      fvm(1,1,l) = sumvx
      fvm(2,1,l) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(1,2,l) = sumvy
      fvm(2,2,l) = dsqrt(sumvy2*anp - sumvy**2)
      sumvz = sumvz*anp
      fvm(1,3,l) = sumvz
      fvm(2,3,l) = dsqrt(sumvz2*anp - sumvz**2)
c calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 30 j = 2, nmvf
      if (fv(j,1,l).gt.0.) then
         sumvx = sumvx + fv(j,1,l)
         sumvx2 = sumvx2 + fv(j,1,l)*dlog(dble(fv(j,1,l)*svx))
      endif
      if (fv(j,2,l).gt.0.) then
         sumvy = sumvy + fv(j,2,l)
         sumvy2 = sumvy2 + fv(j,2,l)*dlog(dble(fv(j,2,l)*svy))
      endif
      if (fv(j,3,l).gt.0.) then
         sumvz = sumvz + fv(j,3,l)
         sumvz2 = sumvz2 + fv(j,3,l)*dlog(dble(fv(j,3,l)*svz))
      endif
   30 continue
      sum7(1) = sumvx
      sum7(2) = sumvy
      sum7(3) = sumvz
      sum7(4) = sumvx2
      sum7(5) = sumvy2
      sum7(6) = sumvz2
      call PDSUM(sum7,work7,6,1)
      sumvx = sum7(1)
      sumvy = sum7(2)
      sumvz = sum7(3)
      sumvx2 = sum7(4)
      sumvy2 = sum7(5)
      sumvz2 = sum7(6)
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fvm(3,1,l) = sumvx
      fvm(3,2,l) = sumvy
      fvm(3,3,l) = sumvz
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSDIST2(fv,fvm,nblok,nmv,nmvf)
c for 2d code, this subroutine calculates entropy from 2d velocity
c distribution
c input: all except fvm, output: fvm
c fv = distribution function, number of particles in each velocity range
c maximum velocity (used for scaling) is contained in first element fv.
c entropy for i-th dimension is contained in fvm(3,i), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
c is uniform in space and distributions in each dimension are
c independent.
c nblok = number of particle partitions
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer nblok, nmv, nmvf
      real fv, fvm
      dimension fv(nmvf,2,nblok), fvm(3,2,nblok)
c local data
      integer j, l
      real anmv, svx, svy
      double precision sumvx, sumvy, sumvx2, sumvy2
      anmv = real(nmv)
      do 20 l = 1, nblok
      svx = anmv/fv(1,1,l)
      svy = anmv/fv(1,2,l)
c calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      do 10 j = 2, nmvf
      if (fv(j,1,l).gt.0.) then
         sumvx = sumvx + fv(j,1,l)
         sumvx2 = sumvx2 + fv(j,1,l)*dlog(dble(fv(j,1,l)*svx))
      endif
      if (fv(j,2,l).gt.0.) then
         sumvy = sumvy + fv(j,2,l)
         sumvy2 = sumvy2 + fv(j,2,l)*dlog(dble(fv(j,2,l)*svy))
      endif
   10 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      fvm(3,1,l) = sumvx
      fvm(3,2,l) = sumvy
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSDIST23(fv,fvm,nblok,nmv,nmvf)
c for 2-1/2d code, this subroutine calculates entropy from 3d velocity
c distribution
c input: all except fvm, output: fvm
c fv = distribution function, number of particles in each velocity range
c maximum velocity (used for scaling) is contained in first element fv.
c entropy for i-th dimension is contained in fvm(3,i), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
c is uniform in space and distributions in each dimension are
c independent.
c nblok = number of particle partitions
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer nblok, nmv, nmvf
      real fv, fvm
      dimension fv(nmvf,3,nblok), fvm(3,3,nblok)
c local data
      integer j, l
      real anmv, svx, svy, svz
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2
      anmv = real(nmv)
      do 20 l = 1, nblok
      svx = anmv/fv(1,1,l)
      svy = anmv/fv(1,2,l)
      svz = anmv/fv(1,3,l)
c calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 10 j = 2, nmvf
      if (fv(j,1,l).gt.0.) then
         sumvx = sumvx + fv(j,1,l)
         sumvx2 = sumvx2 + fv(j,1,l)*dlog(dble(fv(j,1,l)*svx))
      endif
      if (fv(j,2,l).gt.0.) then
         sumvy = sumvy + fv(j,2,l)
         sumvy2 = sumvy2 + fv(j,2,l)*dlog(dble(fv(j,2,l)*svy))
      endif
      if (fv(j,3,l).gt.0.) then
         sumvz = sumvz + fv(j,3,l)
         sumvz2 = sumvz2 + fv(j,3,l)*dlog(dble(fv(j,3,l)*svz))
      endif
   10 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fvm(3,1,l) = sumvx
      fvm(3,2,l) = sumvy
      fvm(3,3,l) = sumvz
   20 continue
      return
      end
