c-----------------------------------------------------------------------
c 2d parallel PIC library for MPI communications
c p2lib.f contains basic communications procedures for 2d code with
c         1d partitions:
c DCOMP2 determines 1d domain decomposition parameters, quadratic
c        interpolation.
c DCOMP2L determines 1d domain decomposition parameters, linear
c         interpolation.
c PCGUARD2 copy guard cells in y for 2d scalar array, quadratic
c          interpolation, and distributed data with uniform partition.
c PNCGUARD2 copy guard cells in y for 2d scalar array, quadratic
c           interpolation, and distributed data with non-uniform
c           partition.
c PCGUARD2L copy guard cells in y for 2d scalar array, linear
c           interpolation, and distributed data with uniform partition.
c PNCGUARD2L copy guard cells in y for 2d scalar array, linear
c            interpolation, and distributed data with non-uniform
c            partition.
c PACGUARD2 add guard cells in y for 3 component vector array,
c           quadratic interpolation, and distributed data with uniform
c           partition.
c PACGUARD22 add guard cells in y for 2 component vector array,
c            quadratic interpolation, and distributed data with uniform
c            partition.
c PAGUARD2 add guard cells in y for scalar array, quadratic
c          interpolation, and distributed data with uniform partition.
c PAMCGUARD2 add guard cells in y for n component tensor array,
c            quadratic interpolation, and distributed data with uniform
c            partition.
c PNACGUARD2 add guard cells in y for 3 component vector array,
c            quadratic interpolation, and distributed data with
c            non-uniform partition.
c PNACGUARD22 add guard cells in y for 2 component vector array,
c             quadratic interpolation, and distributed data with
c             non-uniform partition.
c PNAGUARD2 add guard cells in y for scalar array, quadratic
c           interpolation, and distributed data with non-uniform
c           partition.
c PNAMCGUARD2 add guard cells in y for n component vector array,
c             quadratic interpolation, and distributed data with
c             non-uniform partition.
c PACGUARD2L add guard cells in y for 3 component vector array,
c            linear interpolation, and distributed data with uniform
c            partition.
c PACGUARD22L add guard cells in y for 2 component vector array,
c             quadratic interpolation, and distributed data with uniform
c             partition.
c PAGUARD2L add guard cells in y for scalar array, linear interpolation,
c           and distributed data with uniform partition.
c PAMCGUARD2L add guard cells in y for n component tensor array,
c             linear interpolation, and distributed data with uniform
c             partition.
c PNACGUARD2L add guard cells in y for 3 component vector array,
c             linear interpolation, and distributed data with
c             non-uniform partition.
c PNACGUARD22L add guard cells in y for 2 component vector array,
c              linear interpolation, and distributed data with
c              non-uniform partition.
c PNAGUARD2L add guard cells in y for scalar array, linear
c            interpolation, and distributed data with non-uniform
c            partition.
c PNAMCGUARD2L add guard cells in y for n component vector array,
c              linear interpolation, and distributed data with
c              non-uniform partition.
c PDBLSIN2C creates a doubled 2 component vector array to perform 2d
c           sine/cosine transforms with real to complex ffts, for
c           distributed data.
c PDBLSIN2D creates a doubled scalar array to perform 2d sine transform
c           with real to complex fft, for distributed data.
c PDBLSIN2B creates a doubled 3 component vector array to perform 2d
c           sine/cosine transforms with real to complex ffts, for
c           distributed data.
c PDBLCOS2C creates a doubled 2 component vector array to perform 2d
c           cosine/sine transforms with real to complex ffts, for
c           distributed data.
c PDBLCOS2D creates a doubled scalar array to perform 2d cosine
c           transform with real to complex fft, for distributed data.
c PDBLCOS2B creates a doubled 3 component vector array to perform 2d
c           cosine/sine transforms with real to complex ffts, for
c           distributed data.
c PHAFDBL2C copies data from a double array to regular array with guard
c           cells for 2 component vector field and linear interpolation,
c           for distributed data.
c PHAFDBL2D copies data from a double array to regular array with guard
c           cells for scalar field and linear interpolation, for
c           distributed data.
c PHAFDBL2B copies data from a double array to regular array with guard
c           cells for 3 component vector field and linear interpolation,
c           for distributed data.
c PLCGUARD2 copy guard cells in y for 2 component vector array,
c           quadratic interpolation, and distributed non-periodic data
c           with uniform partition, replicating data to reduce
c           interpolation to linear at edges.
c PLDGUARD2 copy guard cells in y for scalar array, quadratic
c           interpolation, and distributed non-periodic data with
c           uniform partition, replicating data to reduce interpolation
c           to linear at edges.
c PLBGUARD2 copy guard cells in y for 3 component vector array,
c           quadratic interpolation, and distributed non-periodic data
c           with uniform partition, replicating data to reduce
c           interpolation to linear at edges.
c PNLCGUARD2 copy guard cells in y for 2 component vector array,
c            quadratic interpolation, and distributed non-periodic data
c            with non-uniform partition, replicating data to reduce
c            interpolation to linear at edges.
c PNLDGUARD2 copy guard cells in y for scalar array, quadratic
c            interpolation, and distributed non-periodic data with
c            non-uniform partition, replicating data to reduce
c            interpolation to linear at edges.
c PNLBGUARD2 copy guard cells in y for 3 component vector array,
c            quadratic interpolation, and distributed non-periodic data
c            with non-uniform partition, replicating data to reduce
c            interpolation to linear at edges.
c PLCGUARD2L copy guard cells in y for scalar array, linear
c            interpolation, and distributed non-periodic data with
c            uniform partition.
c PNLCGUARD2L copy guard cells in y for scalar array, linear
c             interpolation, and distributed non-periodic data with
c             non-uniform partition.
c PLACGUARD2 add guard cells in y for 3 component vector array,
c            quadratic interpolation, and distributed non-periodic data
c            with uniform partition, reducing interpolation to linear at
c            edges.
c PLACGUARD22 add guard cells in y for 2 component vector array,
c             quadratic interpolation, and distributed non-periodic data
c             with uniform partition, reducing interpolation to linear
c             at edges.
c PLAGUARD2 add guard cells in y for scalar array, quadratic
c           interpolation, and distributed non-periodic data with
c           uniform partition, reducing interpolation to linear at edges.
c PNLACGUARD2 add guard cells in y for 3 component vector array,
c             quadratic interpolation, and distributed non-periodic data
c             with non-uniform partition, reducing interpolation to
c             linear at edges.
c PNLACGUARD22 add guard cells in y for 2 component vector array,
c              quadratic interpolation, and distributed non-periodic
c              data with non-uniform partition, reducing interpolation
c              to linear at edges.
c PNLAGUARD2 add guard cells in y for scalar array, quadratic
c            interpolation, and distributed non-periodic data with
c            uniform partition, reducing interpolation to linear at
c            edges.
c PLACGUARDS2 corrects current density in y for 3 component vector
c             array for particle boundary conditions which keep
c             particles one grid away from the edges, quadratic
c             interpolation, and distributed non-periodic data with
c             uniform partition.
c PLACGUARDS22 corrects current density in y for 2 component vector
c              array for particle boundary conditions which keep
c              particles one grid away from the edges, quadratic
c              interpolation, and distributed non-periodic data with
c              uniform partition.
c PLAGUARDS2 corrects current density in y for scalar array for particle
c            boundary conditions which keep particles one grid away from
c            the edges, quadratic interpolation, and distributed
c            non-periodic data with uniform partition.
c PNLACGUARDS2 corrects current density in y for 3 component vector
c              array for particle boundary conditions which keep
c              particles one grid away from the edges, quadratic
c              interpolation, and distributed non-periodic data with
c              non-uniform partition.
c PNLACGUARDS22 corrects current density in y for 2 component vector
c               array for particle boundary conditions which keep
c               particles one grid away from the edges, quadratic
c               interpolation, and distributed non-periodic data with
c               non-uniform partition.
c PNLAGUARDS2 corrects current density in y for scalar array for particle
c             boundary conditions which keep particles one grid away from
c             the edges, quadratic interpolation, and distributed
c             non-periodic data with non-uniform partition.
c PLACGUARD2L add guard cells in y for 3 component vector array,
c             linear interpolation, and distributed non-periodic data
c             with uniform partition.
c PLACGUARD22L add guard cells in y for 2 component vector array,
c              linear interpolation, and distributed non-periodic data
c              with uniform partition.
c PLAGUARD2L add guard cells in y for scalar array, linear
c            interpolation, and distributed non-periodic data with
c            uniform partition.
c PNLACGUARD2L add guard cells in y for 3 component vector array,
c              linear interpolation, and distributed non-periodic data
c              with non-uniform partition.
c PNLACGUARD22L add guard cells in y for 2 component vector array,
c               linear interpolation, and distributed non-periodic data
c               with non-uniform partition.
c PNLAGUARD2L add guard cells in y for scalar array, linear
c             interpolation, and distributed non-periodic data with
c             non-uniform partition.
c PZDBL2C creates a doubled 2 component vector array to perform 2d
c         convolutions with real to complex ffts, for distributed data.
c PZDBL2D creates a doubled scalar array to perform 2d convolutions with
c         real to complex ffts, for distributed data.
c PZDBL2B creates a doubled 3 component vector array to perform 2d
c         convolutions with real to complex ffts, for distributed data.
c PMOVE2 moves particles into appropriate spatial regions with periodic
c        boundary conditions.
c PXMOV2 moves particles into appropriate spatial regions with periodic
c        boundary conditions, for vector processors.
c WPMOVE2 wrapper function for particle manager, which moves particles
c         into appropriate spatial regions with periodic boundary
c         conditions.
c WPXMOV2 wrapper function for particle manager, which moves particles
c         into appropriate spatial regions with periodic boundary
c         conditions, for vector processors.
c WPMOVES2 wrapper function for particle manager, which moves particles
c          into appropriate spatial regions with periodic boundary
c          conditions.  Optimized for fixed number of particle passes.
c WPXMOVS2 wrapper function for particle manager, which moves particles
c          into appropriate spatial regions with periodic boundary
c          conditions, for vector processors.  Optimized for fixed
c          number of particle passes.
c PMOVEH2 creates ihole list of particles which are leaving this
c         processor.
c PMOVEHX2 creates ihole list  of particles which are leaving this
c          processor, for vector processors.
c PMOVES2 moves particles into appropriate spatial regions with periodic
c         boundary conditions.  Assumes ihole list has been found.
c PMOVESS2 moves particles into appropriate spatial regions with
c          periodic boundary conditions.  Assumes ihole list has been
c          found.  Optimized for fixed number of particle passes.
c PFMOVE2 moves fields into appropriate spatial regions, between
c         non-uniform and uniform partitions.
c REPART2 finds new partitions boundaries (edges, noff, nyp) from old
c         partition information.
c REPARTD2 finds new partitions boundaries (edges, noff, nyp) from old
c          partition information.
c FNOFF2 finds new partitions arrays (noff,nyp) from edges.
c PTPOSE performs a transpose of a complex scalar array, distributed
c        in y, to a complex scalar array, distributed in x.
c P2TPOSE performs a transpose of a 2 component complex vector array,
c         distributed in y, to a 2 component complex vector array,
c         distributed in x.
c P3TPOSE performs a transpose of a 3 component complex vector array,
c         distributed in y, to a 3 component complex vector array,
c         distributed in x.
c PTPOSEX performs a transpose of a complex scalar array, distributed
c         in y, to a complex scalar array, distributed in x, with
c         bufferless algorithm.
c P2TPOSEX performs a transpose of a 2 component complex vector array,
c          distributed in y, to a 2 component complex vector array,
c          distributed in x, with bufferless algorithm.
c P3TPOSEX performs a transpose of a 3 component complex vector array,
c          distributed in y, to a 3 component complex vector array,
c          distributed in x, with bufferless algorithm.
c PNTPOSE performs a transpose of an n component complex vector array,
c         distributed in y, to an n component complex vector array,
c         distributed in x.
c PNTPOSEX performs a transpose of an n component complex vector array,
c          distributed in y, to an n component complex vector array,
c          distributed in x, with bufferless algorithm.
c PN2TPOSE performs a transpose of two n component complex vector
c          arrays, distributed in y, to two n component complex vector
c          array, distributed in x.
c PRTPOSE performs a transpose of a real scalar array, distributed in y,
c         to a complex real array, distributed in x.
c PR2TPOSE performs a transpose of a 2 component real vector array,
c          distributed in y, to a 2 component real vector array,
c          distributed in x.
c PR3TPOSE performs a transpose of a 3 component real vector array,
c          distributed in y, to a 3 component real vector array,
c          distributed in x.
c PWRITE2 collects distributed real 2d data f and writes to a direct
c         access binary file.
c PREAD2 reads real 2d data from a direct access binary file and
c        distributes it.
c PCWRITE2 collects distributed complex 2d data f and writes to a direct
c          access binary file.
c PCREAD2 reads complex 2d data from a direct access binary file and
c         distributes it.
c PCDIFF2 performs a centered finite difference calculation in the
c         second index of a two dimensional array whose second index is
c         distributed.
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: october 4, 2010
c-----------------------------------------------------------------------
      subroutine DCOMP2(edges,nyp,noff,ny,kstrt,nvp,idps,nblok)
c this subroutine determines spatial boundaries for particle
c decomposition, calculates number of grid points in each spatial
c region, and the offset of these grid points from the global address
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c nyp(l) = number of primary gridpoints in particle partition l.
c noff(l) = lowermost global gridpoint in particle partition l.
c ny = system length in y direction
c kstrt = starting data block number
c nvp = number of real or virtual processors
c idps = number of partition boundaries
c nblok = number of particle partitions.
      implicit none
      real edges
      integer nyp, noff, ny, kstrt, nvp, idps, nblok
      dimension edges(idps,nblok)
      dimension nyp(nblok), noff(nblok)
c local data
      integer ks, kb, kr, l
      real at1
      ks = kstrt - 2
      at1 = float(ny)/float(nvp)
      do 10 l = 1, nblok
      kb = l + ks
      edges(1,l) = at1*float(kb)
      noff(l) = edges(1,l) + .5
      edges(2,l) = at1*float(kb + 1)
      kr = edges(2,l) + .5
      nyp(l) = kr - noff(l)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DCOMP2L(edges,nyp,noff,ny,kstrt,nvp,idps,nblok)
c this subroutine determines spatial boundaries for particle
c decomposition, calculates number of grid points in each spatial
c region, and the offset of these grid points from the global address
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c nyp(l) = number of primary gridpoints in particle partition l.
c noff(l) = lowermost global gridpoint in particle partition l.
c ny = system length in y direction
c kstrt = starting data block number
c nvp = number of real or virtual processors
c idps = number of partition boundaries
c nblok = number of particle partitions.
      dimension edges(idps,nblok)
      dimension nyp(nblok), noff(nblok)
      ks = kstrt - 2
      at1 = float(ny)/float(nvp)
      do 10 l = 1, nblok
      kb = l + ks
      edges(1,l) = at1*float(kb)
      noff(l) = edges(1,l)
      edges(2,l) = at1*float(kb + 1)
      kr = edges(2,l)
      nyp(l) = kr - noff(l)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCGUARD2(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes three extra guard cells.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer kstrt, nvp, nxv, nypmx, kyp, kblok
      dimension f(nxv,nypmx,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ks, moff, kr, krr, kl, kll, ngc, l
c     integer j
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 30 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         if (krr.gt.nvp) krr = krr - nvp
         kll = kll - 1
         if (kll.lt.1) kll = kll + nvp
         ngc = 1
      endif
c this segment is used for shared memory computers
c     do 10 j = 1, nxv
c     f(j,1,l) = f(j,kyp+1,kl)
c     f(j,kyp+2,l) = f(j,2,kr)
c     f(j,kyp+3,l) = f(j,ngc+1,krr)
c  10 continue
c this segment is used for mpi computers
      call MPI_IRECV(f(1,1,l),nxv,mreal,kl-1,moff+3,lgrp,msid,ierr)
      call MPI_SEND(f(1,kyp+1,l),nxv,mreal,kr-1,moff+3,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(f(1,kyp+2,l),ngc*nxv,mreal,kr-1,moff+4,lgrp,msid,ie
     1rr)
      call MPI_SEND(f(1,2,l),ngc*nxv,mreal,kl-1,moff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         call MPI_IRECV(f(1,kyp+3,l),ngc*nxv,mreal,krr-1,moff+6,lgrp,msi
     1d,ierr)
         call MPI_SEND(f(1,2,l),ngc*nxv,mreal,kll-1,moff+6,lgrp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNCGUARD2(f,scs,nyp,kstrt,nvp,nxv,nypmx,nblok,mter)
c this subroutine copies data to guard cells in non-uniform partitions
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c scs(j,l) = scratch array for field partition l
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scs
      integer nyp
      integer kstrt, nvp, nxv, nypmx, nblok, mter
      dimension f(nxv,nypmx,nblok), scs(nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ks, moff, kr, krr, kl, kll, ngc, nps, l
c     integer j
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 30 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr + 1
      if (krr.gt.nvp) krr = krr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl - 1
      if (kll.lt.1) kll = kll + nvp
      ngc = 0
c special case of only one grid per processor
      if (nyp(l).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyp(kr).eq.1) then
c        do 10 j = 1, nxv
c        f(j,1,l) = f(j,nyp(kl)+1,kl)
c        f(j,nyp(l)+2,l) = f(j,2,kr)
c        f(j,nyp(l)+3,l) = f(j,2,krr)
c  10    continue
c     else
c        do 20 j = 1, nxv
c        f(j,1,l) = f(j,nyp(kl)+1,kl)
c        f(j,nyp(l)+2,l) = f(j,2,kr)
c        f(j,nyp(l)+3,l) = f(j,3,kr)
c  20    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(f(1,1,l),nxv,mreal,kl-1,moff+3,lgrp,msid,ierr)
      call MPI_SEND(f(1,nyp(l)+1,l),nxv,mreal,kr-1,moff+3,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(f(1,nyp(l)+2,l),2*nxv,mreal,kr-1,moff+4,lgrp,msid,i
     1err)
      call MPI_SEND(f(1,2,l),(2-ngc)*nxv,mreal,kl-1,moff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor
      if (mter.ge.1) go to 30
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (nps.eq.nxv) then
         call MPI_IRECV(f(1,nyp(l)+3,l),nxv,mreal,krr-1,moff+6,lgrp,msid
     1,ierr)
      else
         call MPI_IRECV(scs,nxv,mreal,krr-1,moff+6,lgrp,msid,ierr)
      endif
      call MPI_SEND(f(1,2,l),nxv,mreal,kll-1,moff+6,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes one extra guard cell.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f
      integer kstrt, nvp, nxv, nypmx, kyp, kblok
      dimension f(nxv,nypmx,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ks, moff, kl, kr, l
c     integer j
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 20 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) then
         kr = kr - nvp
      endif
      kl = l + ks
      if (kl.lt.1) then
         kl = kl + nvp
      endif
c this loop is used for shared memory computers
c     do 10 j = 1, nxv
c     f(j,kyp+1,l) = f(j,1,kr)
c  10 continue
c this segment is used for mpi computers
      call MPI_IRECV(f(1,kyp+1,l),nxv,mreal,kr-1,moff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,l),nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
c this subroutine copies data to guard cells in non-uniform partitions
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes one extra guard cell.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cell.
c nblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f
      integer nyp
      integer kstrt, nvp, nxv, nypmx, nblok
      dimension f(nxv,nypmx,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ks, moff, kl, kr, l
c     integer j
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 20 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1)  kl = kl + nvp
c this loop is used for shared memory computers
c     do 10 j = 1, nxv
c     f(j,nyp(l)+1,l) = f(j,1,kr)
c  10 continue
c this segment is used for mpi computers
      call MPI_IRECV(f(1,nyp(l)+1,l),nxv,mreal,kr-1,moff+2,lgrp,msid,ier
     1r)
      call MPI_SEND(f(1,1,l),nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(3,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes three extra guard cells.
c scr(3,j,ngds,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = second dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
      dimension f(3,nxv,nypmx,kblok), scr(3,nxv,ngds,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, kblok
         do 20 j = 1, nx3
         do 10 m = 1, 3
         f(m,j,2,l) = f(m,j,2,l) + f(m,j,kyp+2,l)
         f(m,j,3,l) = f(m,j,3,l) + f(m,j,kyp+3,l)
         f(m,j,kyp+1,l) = f(m,j,kyp+1,l) + f(m,j,1,l)
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         if (krr.gt.nvp) krr = krr - nvp
         kll = kll - 1
         if (kll.lt.1) kll = kll + nvp
         ngc = 1
      endif
c this segment is used for shared memory computers
c     do 50 j = 1, nx3
c     do 40 m = 1, 3
c     scr(m,j,1,l) = f(m,j,kyp+2,kl)
c     scr(m,j,2,l) = f(m,j,kyp+3,kll)
c     scr(m,j,3,l) = f(m,j,1,kr)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*ngc*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kyp+2,l),3*ngc*nxv,mreal,kr-1,moff+1,lgrp,ierr
     1)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,3,l),3*nxv,mreal,kr-1,moff+2,lgrp,msid,ierr
     1)
      call MPI_SEND(f(1,1,1,l),3*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         call MPI_IRECV(scr(1,1,2,l),3*ngc*nxv,mreal,kll-1,moff+5,lgrp,m
     1sid,ierr)
         call MPI_SEND(f(1,1,kyp+3,l),3*ngc*nxv,mreal,krr-1,moff+5,lgrp,
     1ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
c add up the guard cells
      do 70 j = 1, nx3
      do 60 m = 1, 3
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,ngc+1,l) = f(m,j,ngc+1,l) + scr(m,j,2,l)
      f(m,j,kyp+1,l) = f(m,j,kyp+1,l) + scr(m,j,3,l)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD22(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(2,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes three extra guard cells.
c scr(2,j,ngds,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = second dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
      dimension f(2,nxv,nypmx,kblok), scr(2,nxv,ngds,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, kblok
         do 20 j = 1, nx3
         do 10 m = 1, 2
         f(m,j,2,l) = f(m,j,2,l) + f(m,j,kyp+2,l)
         f(m,j,3,l) = f(m,j,3,l) + f(m,j,kyp+3,l)
         f(m,j,kyp+1,l) = f(m,j,kyp+1,l) + f(m,j,1,l)
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         if (krr.gt.nvp) krr = krr - nvp
         kll = kll - 1
         if (kll.lt.1) kll = kll + nvp
         ngc = 1
      endif
c this segment is used for shared memory computers
c     do 50 j = 1, nx3
c     do 40 m = 1, 2
c     scr(m,j,1,l) = f(m,j,kyp+2,kl)
c     scr(m,j,2,l) = f(m,j,kyp+3,kll)
c     scr(m,j,3,l) = f(m,j,1,kr)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,2*ngc*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kyp+2,l),2*ngc*nxv,mreal,kr-1,moff+1,lgrp,ierr
     1)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,3,l),2*nxv,mreal,kr-1,moff+2,lgrp,msid,ierr
     1)
      call MPI_SEND(f(1,1,1,l),2*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         call MPI_IRECV(scr(1,1,2,l),2*ngc*nxv,mreal,kll-1,moff+5,lgrp,m
     1sid,ierr)
         call MPI_SEND(f(1,1,kyp+3,l),2*ngc*nxv,mreal,krr-1,moff+5,lgrp,
     1ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
c add up the guard cells
      do 70 j = 1, nx3
      do 60 m = 1, 2
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,ngc+1,l) = f(m,j,ngc+1,l) + scr(m,j,2,l)
      f(m,j,kyp+1,l) = f(m,j,kyp+1,l) + scr(m,j,3,l)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes three extra guard cells.
c scr(j,ngds,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
      dimension f(nxv,nypmx,kblok), scr(nxv,ngds,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l
      dimension istatus(lstat)
      nx3 = nx + 3
c special case for one processor
      if (nvp.eq.1) then
         do 20 l = 1, kblok
         do 10 j = 1, nx3
         f(j,2,l) = f(j,2,l) + f(j,kyp+2,l)
         f(j,3,l) = f(j,3,l) + f(j,kyp+3,l)
         f(j,kyp+1,l) = f(j,kyp+1,l) + f(j,1,l)
   10    continue
   20    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 50 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         if (krr.gt.nvp) krr = krr - nvp
         kll = kll - 1
         if (kll.lt.1) kll = kll + nvp
         ngc = 1
      endif
c this segment is used for shared memory computers
c     do 30 j = 1, nx3
c     scr(j,1,l) = f(j,kyp+2,kl)
c     scr(j,2,l) = f(j,kyp+3,kll)
c     scr(j,3,l) = f(j,1,kr)
c  30 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,ngc*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,kyp+2,l),ngc*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,3,l),nxv,mreal,kr-1,moff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,l),nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         call MPI_IRECV(scr(1,2,l),ngc*nxv,mreal,kll-1,moff+5,lgrp,msid,
     1ierr)
         call MPI_SEND(f(1,kyp+3,l),ngc*nxv,mreal,krr-1,moff+5,lgrp,ierr
     1)
         call MPI_WAIT(msid,istatus,ierr)
      endif
c add up the guard cells
      do 40 j = 1, nx3
      f(j,2,l) = f(j,2,l) + scr(j,1,l)
      f(j,ngc+1,l) = f(j,ngc+1,l) + scr(j,2,l)
      f(j,kyp+1,l) = f(j,kyp+1,l) + scr(j,3,l)
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAMCGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds,
     1ndim)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(ndim,j,k,l) = real data for grid j,k in particle partition l. number
c of grids per partition is uniform and includes three extra guard cells
c scr(ndim,j,ngds,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = second dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c ngds = number of guard cells
c ndim = first dimension of f
c quadratic interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds, ndim
      dimension f(ndim,nxv,nypmx,kblok), scr(ndim,nxv,ngds,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, kblok
         do 20 j = 1, nx3
         do 10 m = 1, ndim
         f(m,j,2,l) = f(m,j,2,l) + f(m,j,kyp+2,l)
         f(m,j,3,l) = f(m,j,3,l) + f(m,j,kyp+3,l)
         f(m,j,kyp+1,l) = f(m,j,kyp+1,l) + f(m,j,1,l)
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         if (krr.gt.nvp) krr = krr - nvp
         kll = kll - 1
         if (kll.lt.1) kll = kll + nvp
         ngc = 1
      endif
c this segment is used for shared memory computers
c     do 50 j = 1, nx3
c     do 40 m = 1, ndim
c     scr(m,j,1,l) = f(m,j,kyp+2,kl)
c     scr(m,j,2,l) = f(m,j,kyp+3,kll)
c     scr(m,j,3,l) = f(m,j,1,kr)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,ndim*ngc*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kyp+2,l),ndim*ngc*nxv,mreal,kr-1,moff+1,lgrp,i
     1err)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,3,l),ndim*nxv,mreal,kr-1,moff+2,lgrp,msid,i
     1err)
      call MPI_SEND(f(1,1,1,l),ndim*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         call MPI_IRECV(scr(1,1,2,l),ndim*ngc*nxv,mreal,kll-1,moff+5,lgr
     1p,msid,ierr)
         call MPI_SEND(f(1,1,kyp+3,l),ndim*ngc*nxv,mreal,krr-1,moff+5,lg
     1rp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
c add up the guard cells
      do 70 j = 1, nx3
      do 60 m = 1, ndim
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,ngc+1,l) = f(m,j,ngc+1,l) + scr(m,j,2,l)
      f(m,j,kyp+1,l) = f(m,j,kyp+1,l) + scr(m,j,3,l)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNACGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,n
     1gds,mter)
c this subroutine adds data from guard cells in non-uniform partitions
c f(3,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c scr(3,j,ngds,l) = scratch array for field partition l
c scs(3,j,l) = scratch array for field partition l
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c ngds = number of guard cells
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scr, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
      dimension f(3,nxv,nypmx,nblok), scr(3,nxv,ngds,nblok)
      dimension scs(3,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, nblok
         do 20 j = 1, nx3
         do 10 m = 1, 3
         f(m,j,2,l) = f(m,j,2,l) + f(m,j,nyp(l)+2,l)
         f(m,j,3,l) = f(m,j,3,l) + f(m,j,nyp(l)+3,l)
         f(m,j,nyp(l)+1,l) = f(m,j,nyp(l)+1,l) + f(m,j,1,l)
         f(m,j,1,l) = 0.
         f(m,j,nyp(l)+2,l) = 0.
         f(m,j,nyp(l)+3,l) = 0.
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 110 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr + 1
      if (krr.gt.nvp) krr = krr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl - 1
      if (kll.lt.1) kll = kll + nvp
c this segment is used for shared memory computers
c     if (nyp(kl).eq.1) then
c        do 50 j = 1, nx3
c        do 40 m = 1, 3
c        scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl) + f(m,j,nyp(kll)+3,kll)
c        scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c        scr(m,j,3,l) = f(m,j,1,kr)
c  40    continue
c  50    continue
c     else
c        do 70 j = 1, nx3
c        do 60 m = 1, 3
c        scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl)
c        scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c        scr(m,j,3,l) = f(m,j,1,kr)
c  60    continue
c  70    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr,6*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nyp(l)+2,l),6*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,3,l),3*nxv,mreal,kr-1,moff+2,lgrp,msid,ierr
     1)
      call MPI_SEND(f(1,1,1,l),3*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor
      if (mter.ge.1) go to 80
      call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      call MPI_IRECV(scs,3*nxv,mreal,kll-1,moff+5,lgrp,nsid,ierr)
      call MPI_SEND(nyp(l),1,mint,kr-1,moff+3,lgrp,ierr)
      call MPI_SEND(f(1,1,nyp(l)+3,l),3*nxv,mreal,krr-1,moff+5,lgrp,ierr
     1)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 50 j = 1, nx3
         do 40 m = 1, 3
         scr(m,j,1,l) = scr(m,j,1,l) + scs(m,j,l)
   40    continue
   50    continue
      endif
c add up the guard cells
   80 do 100 j = 1, nx3
      do 90 m = 1, 3
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,3,l) = f(m,j,3,l) + scr(m,j,2,l)
      f(m,j,nyp(l)+1,l) = f(m,j,nyp(l)+1,l) + scr(m,j,3,l)
      f(m,j,1,l) = 0.
      f(m,j,nyp(l)+2,l) = 0.
      f(m,j,nyp(l)+3,l) = 0.
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNACGUARD22(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,
     1ngds,mter)
c this subroutine adds data from guard cells in non-uniform partitions
c f(2,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c scr(2,j,ngds,l) = scratch array for field partition l
c scs(2,j,l) = scratch array for field partition l
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c ngds = number of guard cells
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scr, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
      dimension f(2,nxv,nypmx,nblok), scr(2,nxv,ngds,nblok)
      dimension scs(2,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, nblok
         do 20 j = 1, nx3
         do 10 m = 1, 2
         f(m,j,2,l) = f(m,j,2,l) + f(m,j,nyp(l)+2,l)
         f(m,j,3,l) = f(m,j,3,l) + f(m,j,nyp(l)+3,l)
         f(m,j,nyp(l)+1,l) = f(m,j,nyp(l)+1,l) + f(m,j,1,l)
         f(m,j,1,l) = 0.
         f(m,j,nyp(l)+2,l) = 0.
         f(m,j,nyp(l)+3,l) = 0.
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 110 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr + 1
      if (krr.gt.nvp) krr = krr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl - 1
      if (kll.lt.1) kll = kll + nvp
c this segment is used for shared memory computers
c     if (nyp(kl).eq.1) then
c        do 50 j = 1, nx3
c        do 40 m = 1, 2
c        scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl) + f(m,j,nyp(kll)+3,kll)
c        scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c        scr(m,j,3,l) = f(m,j,1,kr)
c  40    continue
c  50    continue
c     else
c        do 70 j = 1, nx3
c        do 60 m = 1, 2
c        scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl)
c        scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c        scr(m,j,3,l) = f(m,j,1,kr)
c  60    continue
c  70    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr,4*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nyp(l)+2,l),4*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,3,l),2*nxv,mreal,kr-1,moff+2,lgrp,msid,ierr
     1)
      call MPI_SEND(f(1,1,1,l),2*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor
      if (mter.ge.1) go to 80
      call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      call MPI_IRECV(scs,2*nxv,mreal,kll-1,moff+5,lgrp,nsid,ierr)
      call MPI_SEND(nyp(l),1,mint,kr-1,moff+3,lgrp,ierr)
      call MPI_SEND(f(1,1,nyp(l)+3,l),2*nxv,mreal,krr-1,moff+5,lgrp,ierr
     1)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 50 j = 1, nx3
         do 40 m = 1, 2
         scr(m,j,1,l) = scr(m,j,1,l) + scs(m,j,l)
   40    continue
   50    continue
      endif
c add up the guard cells
   80 do 100 j = 1, nx3
      do 90 m = 1, 2
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,3,l) = f(m,j,3,l) + scr(m,j,2,l)
      f(m,j,nyp(l)+1,l) = f(m,j,nyp(l)+1,l) + scr(m,j,3,l)
      f(m,j,1,l) = 0.
      f(m,j,nyp(l)+2,l) = 0.
      f(m,j,nyp(l)+3,l) = 0.
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNAGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,ng
     1ds,mter)
c this subroutine adds data from guard cells in non-uniform partitions
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c scr(j,ngds,l) = scratch array for field partition l
c scs(j,l) = scratch array for field partition l
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c ngds = number of guard cells
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scr, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
      dimension f(nxv,nypmx,nblok), scr(nxv,ngds,nblok)
      dimension scs(nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l
      dimension istatus(lstat)
      nx3 = nx + 3
c special case for one processor
      if (nvp.eq.1) then
         do 20 l = 1, nblok
         do 10 j = 1, nx3
         f(j,2,l) = f(j,2,l) + f(j,nyp(l)+2,l)
         f(j,3,l) = f(j,3,l) + f(j,nyp(l)+3,l)
         f(j,nyp(l)+1,l) = f(j,nyp(l)+1,l) + f(j,1,l)
         f(j,1,l) = 0.
         f(j,nyp(l)+2,l) = 0.
         f(j,nyp(l)+3,l) = 0.
   10    continue
   20    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 70 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr + 1
      if (krr.gt.nvp) krr = krr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl - 1
      if (kll.lt.1) kll = kll + nvp
c this segment is used for shared memory computers
c     if (nyp(kl).eq.1) then
c        do 30 j = 1, nx3
c        scr(j,1,l) = f(j,nyp(kl)+2,kl) + f(j,nyp(kll)+3,kll)
c        scr(j,2,l) = f(j,nyp(kl)+3,kl)
c        scr(j,3,l) = f(j,1,kr)
c  30    continue
c     else
c        do 40 j = 1, nx3
c        scr(j,1,l) = f(j,nyp(kl)+2,kl)
c        scr(j,2,l) = f(j,nyp(kl)+3,kl)
c        scr(j,3,l) = f(j,1,kr)
c  40    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr,2*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,nyp(l)+2,l),2*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,3,l),nxv,mreal,kr-1,moff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,l),nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor
      if (mter.ge.1) go to 50
      call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      call MPI_IRECV(scs,nxv,mreal,kll-1,moff+5,lgrp,nsid,ierr)
      call MPI_SEND(nyp(l),1,mint,kr-1,moff+3,lgrp,ierr)
      call MPI_SEND(f(1,nyp(l)+3,l),nxv,mreal,krr-1,moff+5,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 30 j = 1, nx3
         scr(j,1,l) = scr(j,1,l) + scs(j,l)
   30    continue
      endif
c add up the guard cells
   50 do 60 j = 1, nx3
      f(j,2,l) = f(j,2,l) + scr(j,1,l)
      f(j,3,l) = f(j,3,l) + scr(j,2,l)
      f(j,nyp(l)+1,l) = f(j,nyp(l)+1,l) + scr(j,3,l)
      f(j,1,l) = 0.
      f(j,nyp(l)+2,l) = 0.
      f(j,nyp(l)+3,l) = 0.
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNAMCGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,
     1ngds,ndim,mter)
c this subroutine adds data from guard cells in non-uniform partitions
c f(ndim,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c scr(ndim,j,ngds,l) = scratch array for field partition l
c scs(ndim,j,l) = scratch array for field partition l
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c ngds = number of guard cells
c ndim = first dimension of f
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scr, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, ndim, mter
      dimension f(ndim,nxv,nypmx,nblok), scr(ndim,nxv,ngds,nblok)
      dimension scs(ndim,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, nblok
         do 20 j = 1, nx3
         do 10 m = 1, ndim
         f(m,j,2,l) = f(m,j,2,l) + f(m,j,nyp(l)+2,l)
         f(m,j,3,l) = f(m,j,3,l) + f(m,j,nyp(l)+3,l)
         f(m,j,nyp(l)+1,l) = f(m,j,nyp(l)+1,l) + f(m,j,1,l)
         f(m,j,1,l) = 0.
         f(m,j,nyp(l)+2,l) = 0.
         f(m,j,nyp(l)+3,l) = 0.
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 110 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      krr = kr + 1
      if (krr.gt.nvp) krr = krr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
      kll = kl - 1
      if (kll.lt.1) kll = kll + nvp
c this segment is used for shared memory computers
c     if (nyp(kl).eq.1) then
c        do 50 j = 1, nx3
c        do 40 m = 1, ndim
c        scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl) + f(m,j,nyp(kll)+3,kll)
c        scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c        scr(m,j,3,l) = f(m,j,1,kr)
c  40    continue
c  50    continue
c     else
c        do 70 j = 1, nx3
c        do 60 m = 1, ndim
c        scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl)
c        scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c        scr(m,j,3,l) = f(m,j,1,kr)
c  60    continue
c  70    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr,2*ndim*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nyp(l)+2,l),2*ndim*nxv,mreal,kr-1,moff+1,lgrp,
     1ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,3,l),ndim*nxv,mreal,kr-1,moff+2,lgrp,msid,i
     1err)
      call MPI_SEND(f(1,1,1,l),ndim*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor
      if (mter.ge.1) go to 80
      call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      call MPI_IRECV(scs,ndim*nxv,mreal,kll-1,moff+5,lgrp,nsid,ierr)
      call MPI_SEND(nyp(l),1,mint,kr-1,moff+3,lgrp,ierr)
      call MPI_SEND(f(1,1,nyp(l)+3,l),ndim*nxv,mreal,krr-1,moff+5,lgrp,i
     1err)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 50 j = 1, nx3
         do 40 m = 1, ndim
         scr(m,j,1,l) = scr(m,j,1,l) + scs(m,j,l)
   40    continue
   50    continue
      endif
c add up the guard cells
   80 do 100 j = 1, nx3
      do 90 m = 1, ndim
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,3,l) = f(m,j,3,l) + scr(m,j,2,l)
      f(m,j,nyp(l)+1,l) = f(m,j,nyp(l)+1,l) + scr(m,j,3,l)
      f(m,j,1,l) = 0.
      f(m,j,nyp(l)+2,l) = 0.
      f(m,j,nyp(l)+3,l) = 0.
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(3,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes one extra guard cell.
c scr(3,j,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(3,nxv,nypmx,kblok), scr(3,nxv,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, kblok
         do 20 j = 1, nx1
         do 10 m = 1, 3
         f(m,j,1,l) = f(m,j,1,l) + f(m,j,kyp+1,l)
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) then
         kr = kr - nvp
      endif
      kl = l + ks
      if (kl.lt.1) then
         kl = kl + nvp
      endif
c this segment is used for shared memory computers
c     do 50 j = 1, nx1
c     do 40 m = 1, 3
c     scr(m,j,l) = f(m,j,kyp+1,kl)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kyp+1,l),3*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 70 j = 1, nx1
      do 60 m = 1, 3
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD22L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(2,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes one extra guard cell.
c scr(2,j,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(2,nxv,nypmx,kblok), scr(2,nxv,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, kblok
         do 20 j = 1, nx1
         do 10 m = 1, 2
         f(m,j,1,l) = f(m,j,1,l) + f(m,j,kyp+1,l)
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) then
         kr = kr - nvp
      endif
      kl = l + ks
      if (kl.lt.1) then
         kl = kl + nvp
      endif
c this segment is used for shared memory computers
c     do 50 j = 1, nx1
c     do 40 m = 1, 2
c     scr(m,j,l) = f(m,j,kyp+1,kl)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,2*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kyp+1,l),2*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 70 j = 1, nx1
      do 60 m = 1, 2
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes one extra guard cell.
c scr(j,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(nxv,nypmx,kblok), scr(nxv,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 20 l = 1, kblok
         do 10 j = 1, nx1
         f(j,1,l) = f(j,1,l) + f(j,kyp+1,l)
   10    continue
   20    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 50 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) then
         kr = kr - nvp
      endif
      kl = l + ks
      if (kl.lt.1) then
         kl = kl + nvp
      endif
c this segment is used for shared memory computers
c     do 30 j = 1, nx1
c     scr(j,l) = f(j,kyp+1,kl)
c  30 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,kyp+1,l),nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 40 j = 1, nx1
      f(j,1,l) = f(j,1,l) + scr(j,l)
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAMCGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ndim
     1)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(ndim,j,k,l) = real data for grid j,k in particle partition l. number
c of grids per partition is uniform and includes one extra guard cell.
c scr(ndim,j,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = second dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c ndim = first dimension of f
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ndim
      dimension f(ndim,nxv,nypmx,kblok), scr(ndim,nxv,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, kblok
         do 20 j = 1, nx1
         do 10 m = 1, ndim
         f(m,j,1,l) = f(m,j,1,l) + f(m,j,kyp+1,l)
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, kblok
      kr = l + ks + 2
      if (kr.gt.nvp) then
         kr = kr - nvp
      endif
      kl = l + ks
      if (kl.lt.1) then
         kl = kl + nvp
      endif
c this segment is used for shared memory computers
c     do 50 j = 1, nx1
c     do 40 m = 1, ndim
c     scr(m,j,l) = f(m,j,kyp+1,kl)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,ndim*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kyp+1,l),ndim*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 70 j = 1, nx1
      do 60 m = 1, ndim
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNACGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
c this subroutine adds data from guard cells in non-uniform partitions
c f(3,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes one extra guard cell.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c scr(3,j,l) = scratch array for field partition l
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok
      dimension f(3,nxv,nypmx,nblok), scr(3,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, nblok
         do 20 j = 1, nx1
         do 10 m = 1, 3
         f(m,j,1,l) = f(m,j,1,l) + f(m,j,nyp(l)+1,l)
         f(m,j,nyp(l)+1,l) = 0.
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
c this segment is used for shared memory computers
c     do 50 j = 1, nx1
c     do 40 m = 1, 3
c     scr(m,j,l) = f(m,j,nyp(kl)+1,kl)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nyp(l)+1,l),3*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 70 j = 1, nx1
      do 60 m = 1, 3
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
      f(m,j,nyp(l)+1,l) = 0.
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNACGUARD22L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
c this subroutine adds data from guard cells in non-uniform partitions
c f(2,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes one extra guard cell.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c scr(2,j,l) = scratch array for field partition l
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok
      dimension f(2,nxv,nypmx,nblok), scr(2,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, nblok
         do 20 j = 1, nx1
         do 10 m = 1, 2
         f(m,j,1,l) = f(m,j,1,l) + f(m,j,nyp(l)+1,l)
         f(m,j,nyp(l)+1,l) = 0.
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
c this segment is used for shared memory computers
c     do 50 j = 1, nx1
c     do 40 m = 1, 2
c     scr(m,j,l) = f(m,j,nyp(kl)+1,kl)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,2*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nyp(l)+1,l),2*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 70 j = 1, nx1
      do 60 m = 1, 2
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
      f(m,j,nyp(l)+1,l) = 0.
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNAGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
c this subroutine adds data from guard cells in non-uniform partitions
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes one extra guard cell.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c scr(j,l) = scratch array for field partition l
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok
      dimension f(nxv,nypmx,nblok), scr(nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 20 l = 1, nblok
         do 10 j = 1, nx1
         f(j,1,l) = f(j,1,l) + f(j,nyp(l)+1,l)
         f(j,nyp(l)+1,l) = 0.
   10    continue
   20    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 50 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
c this segment is used for shared memory computers
c     do 30 j = 1, nx1
c     scr(j,l) = f(j,nyp(kl)+1,kl)
c  30 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,nyp(l)+1,l),nxv,mreal,kr-1,moff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 40 j = 1, nx1
      f(j,1,l) = f(j,1,l) + scr(j,l)
      f(j,nyp(l)+1,l) = 0.
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNAMCGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,ndi
     1m)
c this subroutine adds data from guard cells in non-uniform partitions
c f(ndim,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes one extra guard cell.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c scr(ndim,j,l) = scratch array for field partition l
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c ndim = first dimension of f
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, ndim
      dimension f(ndim,nxv,nypmx,nblok), scr(ndim,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
c special case for one processor
      if (nvp.eq.1) then
         do 30 l = 1, nblok
         do 20 j = 1, nx1
         do 10 m = 1, ndim
         f(m,j,1,l) = f(m,j,1,l) + f(m,j,nyp(l)+1,l)
         f(m,j,nyp(l)+1,l) = 0.
   10    continue
   20    continue
   30    continue
         return
      endif
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, nblok
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
c this segment is used for shared memory computers
c     do 50 j = 1, nx1
c     do 40 m = 1, ndim
c     scr(m,j,l) = f(m,j,nyp(kl)+1,kl)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,ndim*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nyp(l)+1,l),ndim*nxv,mreal,kr-1,moff+1,lgrp,ie
     1rr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 70 j = 1, nx1
      do 60 m = 1, ndim
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
      f(m,j,nyp(l)+1,l) = 0.
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLSIN2C(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2
     1blok)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in y,
c and y component is an odd function in x.
c linear interpolation for distributed data
c cu2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = second dimension of input array cu, must be >= nx+1
c kyp = number of data values per block in y
c kypd = third dimension of input array cu, must be >= kyp+1
c kyp2 = third dimension of output array cu2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real cu, cu2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension cu(2,nxv,kypd,kblok), cu2(2,2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, nxs, nys, ny2, kyb, kyb2, ks, koff, moff
      integer kk, ll, lm, k1, k2, joff
      real at1, at2
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      ks = kstrt - 2
      moff = kypd + kyb
c copy to double array in x direction
      do 80 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = koff/kyp + 1
      koff = kyp*(l + ks)
      lm = koff/kyp2 + 1
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, nys
         do 10 j = 1, nxs
         cu2(1,j+1,k+1,l) = cu(1,j+1,k+1,l)
         cu2(2,j+1,k+1,l) = cu(2,j+1,k+1,l)
         cu2(1,nx+j+1,k+1,l) = cu(1,nx-j+1,k+1,l)
         cu2(2,nx+j+1,k+1,l) = -cu(2,nx-j+1,k+1,l)
         cu2(1,j+1,ny+k+1,l) = -cu(1,j+1,ny-k+1,l)
         cu2(2,j+1,ny+k+1,l) = cu(2,j+1,ny-k+1,l)
         cu2(1,nx+j+1,ny+k+1,l) = -cu(1,nx-j+1,ny-k+1,l)
         cu2(2,nx+j+1,ny+k+1,l) = -cu(2,nx-j+1,ny-k+1,l)
   10    continue
         cu2(1,1,k+1,l) = cu(1,1,k+1,l)
         cu2(2,1,k+1,l) = 0.
         cu2(1,nx+1,k+1,l) = cu(1,nx+1,k+1,l)
         cu2(2,nx+1,k+1,l) = 0.
         cu2(1,1,k+ny+1,l) = -cu(1,1,ny-k+1,l)
         cu2(2,1,k+ny+1,l) = 0.
         cu2(1,nx+1,k+ny+1,l) = -cu(1,nx+1,ny-k+1,l)
         cu2(2,nx+1,k+ny+1,l) = 0.
   20    continue
         do 30 j = 1, nxs
         cu2(1,j+1,1,l) = 0.
         cu2(2,j+1,1,l) = cu(2,j+1,1,l)
         cu2(1,j+nx+1,1,l) = 0.
         cu2(2,j+nx+1,1,l) = -cu(2,nx-j+1,1,l)
         cu2(1,j+1,ny+1,l) = 0.
         cu2(2,j+1,ny+1,l) = cu(2,j+1,ny+1,l)
         cu2(1,j+nx+1,ny+1,l) = 0.
         cu2(2,j+nx+1,ny+1,l) = -cu(2,nx-j+1,ny+1,l)
   30    continue
         cu2(1,1,1,l) = 0.
         cu2(2,1,1,l) = 0.
         cu2(1,nx+1,1,l) = 0.
         cu2(2,nx+1,1,l) = 0.
         cu2(1,1,ny+1,l) = 0.
         cu2(2,1,ny+1,l) = 0.
         cu2(1,nx+1,ny+1,l) = 0.
         cu2(2,nx+1,ny+1,l) = 0.
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 50 k = 1, kyp
c        do 40 j = 1, nx
c        do 35 i = 1, 2
c        cu2(i,j,k,l) = cu(i,j,k,ll)
c  35    continue
c  40    continue
c  50    continue
c        if (kyp.lt.kyp2) then
c           do 70 k = 1, kyp
c           do 60 j = 1, nx
c           do 55 i = 1, 2
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll+1)
c  55       continue
c  60       continue
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(cu2(1,1,1,l),2*kyp*nxv,mreal,ll-1,moff+1,lgrp,ms
     1id,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),2*kyp*nxv,mreal,ll,moff+1,lg
     1rp,nsid,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(cu(1,1,1,l),2*kyp*nxv,mreal,lm-1,moff+1,lgrp,ierr
     1)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 50 k = 2, kyp
         k1 = kyp - k + 2
         k2 = (k1 - 1)/2 + 1
         joff = nxv*(k1 - 2*k2 + 1)
         do 40 j = 1, nxv
         do 35 i = 1, 2
         cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
   35    continue
   40    continue
   50    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 70 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 60 j = 1, nxv
            do 55 i = 1, 2
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   55       continue
   60       continue
   70       continue
         endif
      endif
   80 continue
c copy to double array in y direction
      do 140 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(l + ks)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 90 j = 1, nx
c           do 85 i = 1, 2
c           cu2(i,j,1,l) = cu(i,j,1,ll+1)
c  85       continue
c  90       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 110 k = 1, kyp
c           do 100 j = 1, nx
c           do 95 i = 1, 2
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll)
c  95       continue
c 100       continue
c 110       continue
c        endif
c        if (kyp.gt.1) then
c           do 130 k = 2, kyp
c           do 120 j = 1, nx
c           do 115 i = 1, 2
c           cu2(i,j,k,l) = cu(i,j,k,ll-1)
c 115       continue
c 120       continue
c 130       continue
c        endif
c     endif
c ny+1 point is special
c     if (kyp2*(l+ks).eq.ny) then
c        do 136 j = 1, nx+1
c        do 135 i = 1, 2
c        cu2(i,j,1,l) = cu(i,j,kyp+1,kyb)
c 135    continue
c 136    continue
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(cu2(1,1,1,l),2*nxv,mreal,ll,moff+2,lgrp,lsid,
     1ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),2*kyp*nxv,mreal,ll-1,moff+2,
     1lgrp,msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(cu2(1,1,2,l),2*(kyp-1)*nxv,mreal,ll-2,moff+2,
     1lgrp,nsid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               call MPI_SEND(cu(1,1,1,l),2*nxv,mreal,lm,moff+2,lgrp,ierr
     1)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(cu(1,1,2,l),2*(kyp-1)*nxv,mreal,lm-1,moff+2
     1,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,l),2*kyp*nxv,mreal,lm-1,moff+2,lgrp,i
     1err)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 100 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 90 j = 1, nxv
            do 85 i = 1, 2
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   85       continue
   90       continue
  100       continue
         endif
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 k = 3, kyp
            k1 = kyp - k + 3
            k2 = k1/2 + 1
            joff = nxv*(k1 - 2*k2 + 2)
            do 110 j = 1, nxv
            do 105 i = 1, 2
            cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
  105       continue
  110       continue
  120       continue
         endif
      endif
c ny+1 point is special
      if (kyp2*(l+ks).eq.ny) then
         call MPI_IRECV(cu2(1,1,1,l),2*nxv,mreal,kyb-1,moff+3,lgrp,msid,
     1ierr)
      endif
      if (kyp*(l+ks+1).eq.ny) then
         kk = ny/kyp2
         call MPI_SEND(cu(1,1,kyp+1,l),2*nxv,mreal,kk,moff+3,lgrp,ierr)
      endif
      if (kyp2*(l+ks).eq.ny) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
  140 continue
c create odd array
      do 200 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 190 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 150 j = 1, nxs
         cu2(1,j+1,k,l) = 0.
         cu2(2,j+1,k,l) = cu2(2,j+1,k,l)
         cu2(1,j+nx+1,k,l) = 0.
         cu2(2,j+nx+1,k,l) = -cu2(2,nx-j+1,k,l)
  150    continue
         cu2(1,1,k,l) = 0.
         cu2(2,1,k,l) = 0.
         cu2(1,nx+1,k,l) = 0.
         cu2(2,nx+1,k,l) = 0.
      else if (kk.le.ny) then
         do 160 j = 1, nxs
         cu2(1,nx+j+1,k,l) = cu2(1,nx-j+1,k,l)
         cu2(2,nx+j+1,k,l) = -cu2(2,nx-j+1,k,l)
  160    continue
         cu2(1,1,k,l) = cu2(1,1,k,l)
         cu2(2,1,k,l) = 0.
         cu2(1,nx+1,k,l) = cu2(1,nx+1,k,l)
         cu2(2,nx+1,k,l) = 0.
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 170 j = 1, nxs
            cu2(1,nx+j+1,k,l) = -cu2(1,nx-j+1,k,l)
            cu2(2,nx+j+1,k,l) = -cu2(2,nx-j+1,k,l)
  170       continue
            cu2(1,1,k,l) = -cu2(1,1,k,l)
            cu2(2,1,k,l) = 0.
            cu2(1,nx+1,k,l) = -cu2(1,nx+1,k,l)
            cu2(2,nx+1,k,l) = 0.
         else
            do 180 j = 1, nxs
            cu2(1,nx+j+1,kyp2-k+2,l) = -cu2(1,nx-j+1,k,l)
            cu2(2,nx+j+1,kyp2-k+2,l) = -cu2(2,nx-j+1,k,l)
  180       continue
            if (k.le.(kyp2/2+1)) then
               at1 = -cu2(1,1,kyp2-k+2,l)
               at2 = -cu2(1,nx+1,kyp2-k+2,l)
               cu2(1,1,kyp2-k+2,l) = -cu2(1,1,k,l)
               cu2(1,nx+1,kyp2-k+2,l) = -cu2(1,nx+1,k,l)
               cu2(1,1,k,l) = at1
               cu2(1,nx+1,k,l) = at2
            endif
            cu2(2,1,kyp2-k+2,l) = 0.
            cu2(2,nx+1,kyp2-k+2,l) = 0.
         endif
      endif
  190 continue
  200 continue
c finish odd array
      do 230 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 220 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 210 j = 1, nxs
         cu2(1,nx-j+1,k,l) = cu2(1,nx+j+1,k,l)
         cu2(2,nx-j+1,k,l) = -cu2(2,nx+j+1,k,l)
  210    continue
      endif
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLSIN2D(q,q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2bl
     1ok)
c this subroutine creates an odd array q2 from an array q, so that
c a 2d sine transform can be performed with a 2d real to complex fft.
c linear interpolation for distributed data
c q2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = first dimension of input array q, must be >= nx+1
c kyp = number of data values per block in y
c kypd = second dimension of input array q, must be >= kyp+1
c kyp2 = second dimension of output array q2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real q, q2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension q(nxv,kypd,kblok), q2(2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer j, k, l, nxs, nys, ny2, kyb, kyb2, ks, koff, moff
      integer kk, ll, lm, k1, k2, joff
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      ks = kstrt - 2
      moff = kypd + kyb
c copy to double array in x direction
      do 80 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = koff/kyp + 1
      koff = kyp*(l + ks)
      lm = koff/kyp2 + 1
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, nys
         do 10 j = 1, nxs
         q2(j+1,k+1,l) = q(j+1,k+1,l)
         q2(nx+j+1,k+1,l) = -q(nx-j+1,k+1,l)
         q2(j+1,ny+k+1,l) = -q(j+1,ny-k+1,l)
         q2(nx+j+1,ny+k+1,l) = q(nx-j+1,ny-k+1,l)
   10    continue
         q2(1,k+1,l) = 0.
         q2(nx+1,k+1,l) = 0.
         q2(1,k+ny+1,l) = 0.
         q2(nx+1,k+ny+1,l) = 0.
   20    continue
         do 30 j = 1, nx
         q2(j,1,l) = 0.
         q2(j+nx,1,l) = 0.
         q2(j,ny+1,l) = 0.
         q2(j+nx,ny+1,l) = 0.
   30    continue
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 50 k = 1, kyp
c        do 40 j = 1, nx
c        q2(j,k,l) = q(j,k,ll)
c  40    continue
c  50    continue
c        if (kyp.lt.kyp2) then
c           do 70 k = 1, kyp
c           do 60 j = 1, nx
c           q2(j,k+kyp,l) = q(j,k,ll+1)
c  60       continue
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(q2(1,1,l),kyp*nxv,mreal,ll-1,moff+1,lgrp,msid,ie
     1rr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(q2(1,kyp+1,l),kyp*nxv,mreal,ll,moff+1,lgrp,ns
     1id,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(q(1,1,l),kyp*nxv,mreal,lm-1,moff+1,lgrp,ierr)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 50 k = 2, kyp
         k1 = kyp - k + 2
         k2 = (k1 - 1)/2 + 1
         joff = nxv*(k1 - 2*k2 + 1)
         do 40 j = 1, nxv
         q2(j,k1,l) = q2(j+joff,k2,l)
   40    continue
   50    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 70 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 60 j = 1, nxv
            q2(j,k1+kyp,l) = q2(j+joff,k2+kyp,l)
   60       continue
   70       continue
         endif
      endif
   80 continue
c copy to double array in y direction
      do 140 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(l + ks)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 90 j = 1, nx
c           q2(j,1,l) = q(j,1,ll+1)
c  90       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 110 k = 1, kyp
c           do 100 j = 1, nx
c           q2(j,k+kyp,l) = q(j,k,ll)
c 100       continue
c 110       continue
c        endif
c        if (kyp.gt.1) then
c           do 130 k = 2, kyp
c           do 120 j = 1, nx
c           q2(j,k,l) = q(j,k,ll-1)
c 120       continue
c 130       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(q2(1,1,l),nxv,mreal,ll,moff+2,lgrp,lsid,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(q2(1,kyp+1,l),kyp*nxv,mreal,ll-1,moff+2,lgrp,
     1msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(q2(1,2,l),(kyp-1)*nxv,mreal,ll-2,moff+2,lgrp,
     1nsid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               call MPI_SEND(q(1,1,l),nxv,mreal,lm,moff+2,lgrp,ierr)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(q(1,2,l),(kyp-1)*nxv,mreal,lm-1,moff+2,lgrp
     1,ierr)
            endif
         else
            call MPI_SEND(q(1,1,l),kyp*nxv,mreal,lm-1,moff+2,lgrp,ierr)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 100 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 90 j = 1, nxv
            q2(j,k1+kyp,l) = q2(j+joff,k2+kyp,l)
   90       continue
  100       continue
         endif
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 k = 3, kyp
            k1 = kyp - k + 3
            k2 = k1/2 + 1
            joff = nxv*(k1 - 2*k2 + 2)
            do 110 j = 1, nxv
            q2(j,k1,l) = q2(j+joff,k2,l)
  110       continue
  120       continue
         endif
      endif
  140 continue
c create odd array
      do 200 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 190 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 150 j = 1, nx
         q2(j,k,l) = 0.
         q2(j+nx,k,l) = 0.
  150    continue
      else if (kk.le.ny) then
         do 160 j = 1, nxs
         q2(nx+j+1,k,l) = -q2(nx-j+1,k,l)
  160    continue
         q2(1,k,l) = 0.
         q2(nx+1,k,l) = 0.
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 170 j = 1, nxs
            q2(nx+j+1,k,l) = q2(nx-j+1,k,l)
  170       continue
         else
            do 180 j = 1, nxs
            q2(nx+j+1,kyp2-k+2,l) = q2(nx-j+1,k,l)
  180       continue
         endif
         q2(1,k,l) = 0.
         q2(nx+1,k,l) = 0.
      endif
  190 continue
  200 continue
c finish odd array
      do 230 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 220 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 210 j = 1, nxs
         q2(nx-j+1,k,l) = -q2(nx+j+1,k,l)
  210    continue
         q2(nx+1,k,l) = -q2(nx+1,k,l)
      endif
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLSIN2B(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2
     1blok)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in y,
c y component is an odd function in x, and the z component is an odd
c function in both x and y.
c linear interpolation for distributed data
c cu2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = first dimension of input array q, must be >= nx+1
c kyp = number of data values per block in y
c kypd = second dimension of input array q, must be >= kyp+1
c kyp2 = second dimension of output array q2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real cu, cu2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension cu(3,nxv,kypd,kblok), cu2(3,2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, nxs, nys, ny2, kyb, kyb2, ks, koff, moff
      integer kk, ll, lm, k1, k2, joff
      real at1, at2
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      ks = kstrt - 2
      moff = kypd + kyb
c copy to double array in x direction
      do 80 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = koff/kyp + 1
      koff = kyp*(l + ks)
      lm = koff/kyp2 + 1
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, nys
         do 10 j = 1, nxs
         cu2(1,j+1,k+1,l) = cu(1,j+1,k+1,l)
         cu2(2,j+1,k+1,l) = cu(2,j+1,k+1,l)
         cu2(3,j+1,k+1,l) = cu(3,j+1,k+1,l)
         cu2(1,nx+j+1,k+1,l) = cu(1,nx-j+1,k+1,l)
         cu2(2,nx+j+1,k+1,l) = -cu(2,nx-j+1,k+1,l)
         cu2(3,nx+j+1,k+1,l) = -cu(3,nx-j+1,k+1,l)
         cu2(1,j+1,ny+k+1,l) = -cu(1,j+1,ny-k+1,l)
         cu2(2,j+1,ny+k+1,l) = cu(2,j+1,ny-k+1,l)
         cu2(3,j+1,ny+k+1,l) = -cu(3,j+1,ny-k+1,l)
         cu2(1,nx+j+1,ny+k+1,l) = -cu(1,nx-j+1,ny-k+1,l)
         cu2(2,nx+j+1,ny+k+1,l) = -cu(2,nx-j+1,ny-k+1,l)
         cu2(3,nx+j+1,ny+k+1,l) = cu(3,nx-j+1,ny-k+1,l)
   10    continue
         cu2(1,1,k+1,l) = cu(1,1,k+1,l)
         cu2(2,1,k+1,l) = 0.
         cu2(3,1,k+1,l) = 0.
         cu2(1,nx+1,k+1,l) = cu(1,nx+1,k+1,l)
         cu2(2,nx+1,k+1,l) = 0.
         cu2(3,nx+1,k+1,l) = 0.
         cu2(1,1,k+ny+1,l) = -cu(1,1,ny-k+1,l)
         cu2(2,1,k+ny+1,l) = 0.
         cu2(3,1,k+ny+1,l) = 0.
         cu2(1,nx+1,k+ny+1,l) = -cu(1,nx+1,ny-k+1,l)
         cu2(2,nx+1,k+ny+1,l) = 0.
         cu2(3,nx+1,k+ny+1,l) = 0.
   20    continue
         do 30 j = 1, nxs
         cu2(1,j+1,1,l) = 0.
         cu2(2,j+1,1,l) = cu(2,j+1,1,l)
         cu2(3,j+1,1,l) = 0.
         cu2(1,j+nx+1,1,l) = 0.
         cu2(2,j+nx+1,1,l) = -cu(2,nx-j+1,1,l)
         cu2(3,j+nx+1,1,l) = 0.
         cu2(1,j+1,ny+1,l) = 0.
         cu2(2,j+1,ny+1,l) = cu(2,j+1,ny+1,l)
         cu2(3,j+1,ny+1,l) = 0.
         cu2(1,j+nx+1,ny+1,l) = 0.
         cu2(2,j+nx+1,ny+1,l) = -cu(2,nx-j+1,ny+1,l)
         cu2(3,j+nx+1,ny+1,l) = 0.
   30    continue
         cu2(1,1,1,l) = 0.
         cu2(2,1,1,l) = 0.
         cu2(3,1,1,l) = 0.
         cu2(1,nx+1,1,l) = 0.
         cu2(2,nx+1,1,l) = 0.
         cu2(3,nx+1,1,l) = 0.
         cu2(1,1,ny+1,l) = 0.
         cu2(2,1,ny+1,l) = 0.
         cu2(3,1,ny+1,l) = 0.
         cu2(1,nx+1,ny+1,l) = 0.
         cu2(2,nx+1,ny+1,l) = 0.
         cu2(3,nx+1,ny+1,l) = 0.
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 50 k = 1, kyp
c        do 40 j = 1, nx+1
c        do 35 i = 1, 3
c        cu2(i,j,k,l) = cu(i,j,k,ll)
c  35    continue
c  40    continue
c  50    continue
c        if (kyp.lt.kyp2) then
c           do 70 k = 1, kyp
c           do 60 j = 1, nx+1
c           do 55 i = 1, 3
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll+1)
c  55       continue
c  60       continue
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(cu2(1,1,1,l),3*kyp*nxv,mreal,ll-1,moff+1,lgrp,ms
     1id,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),3*kyp*nxv,mreal,ll,moff+1,lg
     1rp,nsid,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(cu(1,1,1,l),3*kyp*nxv,mreal,lm-1,moff+1,lgrp,ierr
     1)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 50 k = 2, kyp
         k1 = kyp - k + 2
         k2 = (k1 - 1)/2 + 1
         joff = nxv*(k1 - 2*k2 + 1)
         do 40 j = 1, nxv
         do 35 i = 1, 3
         cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
   35    continue
   40    continue
   50    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 70 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 60 j = 1, nxv
            do 55 i = 1, 3
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   55       continue
   60       continue
   70       continue
         endif
      endif
   80 continue
c copy to double array in y direction
      do 140 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(l + ks)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 90 j = 1, nx+1
c           do 85 i = 1, 3
c           cu2(i,j,1,l) = cu(i,j,1,ll+1)
c  85       continue
c  90       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 110 k = 1, kyp
c           do 100 j = 1, nx+1
c           do 95 i = 1, 3
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll)
c  95       continue
c 100       continue
c 110       continue
c        endif
c        if (kyp.gt.1) then
c           do 130 k = 2, kyp
c           do 120 j = 1, nx+1
c           do 115 i = 1, 3
c           cu2(i,j,k,l) = cu(i,j,k,ll-1)
c 115       continue
c 120       continue
c 130       continue
c        endif
c     endif
c ny+1 point is special
c     if (kyp2*(l+ks).eq.ny) then
c        do 136 j = 1, nx+1
c        do 135 i = 1, 3
c        cu2(i,j,1,l) = cu(i,j,kyp+1,kyb)
c 135    continue
c 136    continue
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(cu2(1,1,1,l),3*nxv,mreal,ll,moff+2,lgrp,lsid,
     1ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),3*kyp*nxv,mreal,ll-1,moff+2,
     1lgrp,msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(cu2(1,1,2,l),3*(kyp-1)*nxv,mreal,ll-2,moff+2,
     1lgrp,nsid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               call MPI_SEND(cu(1,1,1,l),3*nxv,mreal,lm,moff+2,lgrp,ierr
     1)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(cu(1,1,2,l),3*(kyp-1)*nxv,mreal,lm-1,moff+2
     1,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,l),3*kyp*nxv,mreal,lm-1,moff+2,lgrp,i
     1err)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 100 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 90 j = 1, nxv
            do 85 i = 1, 3
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   85       continue
   90       continue
  100       continue
         endif
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 k = 3, kyp
            k1 = kyp - k + 3
            k2 = k1/2 + 1
            joff = nxv*(k1 - 2*k2 + 2)
            do 110 j = 1, nxv
            do 105 i = 1, 3
            cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
  105       continue
  110       continue
  120       continue
         endif
      endif
c ny+1 point is special
      if (kyp2*(l+ks).eq.ny) then
         call MPI_IRECV(cu2(1,1,1,l),3*nxv,mreal,kyb-1,moff+3,lgrp,msid,
     1ierr)
      endif
      if (kyp*(l+ks+1).eq.ny) then
         kk = ny/kyp2
         call MPI_SEND(cu(1,1,kyp+1,l),3*nxv,mreal,kk,moff+3,lgrp,ierr)
      endif
      if (kyp2*(l+ks).eq.ny) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
  140 continue
c create odd array
      do 200 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 190 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 150 j = 1, nxs
         cu2(1,j+1,k,l) = 0.
         cu2(2,j+1,k,l) = cu2(2,j+1,k,l)
         cu2(3,j+1,k,l) = 0.
         cu2(1,j+nx+1,k,l) = 0.
         cu2(2,j+nx+1,k,l) = -cu2(2,nx-j+1,k,l)
         cu2(3,j+nx+1,k,l) = 0.
  150    continue
         cu2(1,1,k,l) = 0.
         cu2(2,1,k,l) = 0.
         cu2(3,1,k,l) = 0.
         cu2(1,nx+1,k,l) = 0.
         cu2(2,nx+1,k,l) = 0.
         cu2(3,nx+1,k,l) = 0.
      else if (kk.le.ny) then
         do 160 j = 1, nxs
         cu2(1,nx+j+1,k,l) = cu2(1,nx-j+1,k,l)
         cu2(2,nx+j+1,k,l) = -cu2(2,nx-j+1,k,l)
         cu2(3,nx+j+1,k,l) = -cu2(3,nx-j+1,k,l)
  160    continue
         cu2(1,1,k,l) = cu2(1,1,k,l)
         cu2(2,1,k,l) = 0.
         cu2(3,1,k,l) = 0.
         cu2(1,nx+1,k,l) = cu2(1,nx+1,k,l)
         cu2(2,nx+1,k,l) = 0.
         cu2(3,nx+1,k,l) = 0.
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 170 j = 1, nxs
            cu2(1,nx+j+1,k,l) = -cu2(1,nx-j+1,k,l)
            cu2(2,nx+j+1,k,l) = -cu2(2,nx-j+1,k,l)
            cu2(3,nx+j+1,k,l) = cu2(3,nx-j+1,k,l)
  170       continue
            cu2(1,1,k,l) = -cu2(1,1,k,l)
            cu2(2,1,k,l) = 0.
            cu2(3,1,k,l) = 0.
            cu2(1,nx+1,k,l) = -cu2(1,nx+1,k,l)
            cu2(2,nx+1,k,l) = 0.
            cu2(3,nx+1,k,l) = 0.
         else
            do 180 j = 1, nxs
            cu2(1,nx+j+1,kyp2-k+2,l) = -cu2(1,nx-j+1,k,l)
            cu2(2,nx+j+1,kyp2-k+2,l) = -cu2(2,nx-j+1,k,l)
            cu2(3,nx+j+1,kyp2-k+2,l) = cu2(3,nx-j+1,k,l)
  180       continue
            if (k.le.(kyp2/2+1)) then
               at1 = -cu2(1,1,kyp2-k+2,l)
               at2 = -cu2(1,nx+1,kyp2-k+2,l)
               cu2(1,1,kyp2-k+2,l) = -cu2(1,1,k,l)
               cu2(1,nx+1,kyp2-k+2,l) = -cu2(1,nx+1,k,l)
               cu2(1,1,k,l) = at1
               cu2(1,nx+1,k,l) = at2
            endif
            cu2(2,1,kyp2-k+2,l) = 0.
            cu2(3,1,kyp2-k+2,l) = 0.
            cu2(2,nx+1,kyp2-k+2,l) = 0.
            cu2(3,nx+1,kyp2-k+2,l) = 0.
         endif
      endif
  190 continue
  200 continue
c finish odd array
      do 230 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 220 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 210 j = 1, nxs
         cu2(1,nx-j+1,k,l) = cu2(1,nx+j+1,k,l)
         cu2(2,nx-j+1,k,l) = -cu2(2,nx+j+1,k,l)
         cu2(3,nx-j+1,k,l) = -cu2(3,nx+j+1,k,l)
  210    continue
      endif
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLCOS2C(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2
     1blok)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d cosine/sine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in y,
c and y component is an odd function in x.
c linear interpolation for distributed data
c cu2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = second dimension of input array cu, must be >= nx+1
c kyp = number of data values per block in y
c kypd = third dimension of input array cu, must be >= kyp+1
c kyp2 = third dimension of output array cu2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real cu, cu2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension cu(2,nxv,kypd,kblok), cu2(2,2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, nxs, nys, ny2, kyb, kyb2, ks, koff, moff
      integer kk, ll, lm, k1, k2, joff
      real at1, at2
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      ks = kstrt - 2
      moff = kypd + kyb
c copy to double array in x direction
      do 80 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = koff/kyp + 1
      koff = kyp*(l + ks)
      lm = koff/kyp2 + 1
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, nys
         do 10 j = 1, nxs
         cu2(1,j+1,k+1,l) = cu(1,j+1,k+1,l)
         cu2(2,j+1,k+1,l) = cu(2,j+1,k+1,l)
         cu2(1,nx+j+1,k+1,l) = -cu(1,nx-j+1,k+1,l)
         cu2(2,nx+j+1,k+1,l) = cu(2,nx-j+1,k+1,l)
         cu2(1,j+1,ny+k+1,l) = cu(1,j+1,ny-k+1,l)
         cu2(2,j+1,ny+k+1,l) = -cu(2,j+1,ny-k+1,l)
         cu2(1,nx+j+1,ny+k+1,l) = -cu(1,nx-j+1,ny-k+1,l)
         cu2(2,nx+j+1,ny+k+1,l) = -cu(2,nx-j+1,ny-k+1,l)
   10    continue
         cu2(1,1,k+1,l) = 0.
         cu2(2,1,k+1,l) = cu(2,1,k+1,l)
         cu2(1,nx+1,k+1,l) = 0.
         cu2(2,nx+1,k+1,l) = cu(2,nx+1,k+1,l)
         cu2(1,1,k+ny+1,l) = 0.
         cu2(2,1,k+ny+1,l) = -cu(2,1,ny-k+1,l)
         cu2(1,nx+1,k+ny+1,l) = 0.
         cu2(2,nx+1,k+ny+1,l) = -cu(2,nx+1,ny-k+1,l)
   20    continue
         do 30 j = 1, nxs
         cu2(1,j+1,1,l) = cu(1,j+1,1,l)
         cu2(2,j+1,1,l) = 0.
         cu2(1,j+nx+1,1,l) = -cu(1,nx-j+1,1,l)
         cu2(2,j+nx+1,1,l) = 0.
         cu2(1,j+1,ny+1,l) = cu(1,j+1,ny+1,l)
         cu2(2,j+1,ny+1,l) = 0.
         cu2(1,j+nx+1,ny+1,l) = -cu(1,nx-j+1,ny+1,l)
         cu2(2,j+nx+1,ny+1,l) = 0.
   30    continue
         cu2(1,1,1,l) = 0.
         cu2(2,1,1,l) = 0.
         cu2(1,nx+1,1,l) = 0.
         cu2(2,nx+1,1,l) = 0.
         cu2(1,1,ny+1,l) = 0.
         cu2(2,1,ny+1,l) = 0.
         cu2(1,nx+1,ny+1,l) = 0.
         cu2(2,nx+1,ny+1,l) = 0.
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 50 k = 1, kyp
c        do 40 j = 1, nx
c        do 35 i = 1, 2
c        cu2(i,j,k,l) = cu(i,j,k,ll)
c  35    continue
c  40    continue
c  50    continue
c        if (kyp.lt.kyp2) then
c           do 70 k = 1, kyp
c           do 60 j = 1, nx
c           do 55 i = 1, 2
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll+1)
c  55       continue
c  60       continue
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(cu2(1,1,1,l),2*kyp*nxv,mreal,ll-1,moff+1,lgrp,ms
     1id,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),2*kyp*nxv,mreal,ll,moff+1,lg
     1rp,nsid,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(cu(1,1,1,l),2*kyp*nxv,mreal,lm-1,moff+1,lgrp,ierr
     1)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 50 k = 2, kyp
         k1 = kyp - k + 2
         k2 = (k1 - 1)/2 + 1
         joff = nxv*(k1 - 2*k2 + 1)
         do 40 j = 1, nxv
         do 35 i = 1, 2
         cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
   35    continue
   40    continue
   50    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 70 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 60 j = 1, nxv
            do 55 i = 1, 2
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   55       continue
   60       continue
   70       continue
         endif
      endif
   80 continue
c copy to double array in y direction
      do 140 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(l + ks)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 90 j = 1, nx
c           do 85 i = 1, 2
c           cu2(i,j,1,l) = cu(i,j,1,ll+1)
c  85       continue
c  90       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 110 k = 1, kyp
c           do 100 j = 1, nx
c           do 95 i = 1, 2
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll)
c  95       continue
c 100       continue
c 110       continue
c        endif
c        if (kyp.gt.1) then
c           do 130 k = 2, kyp
c           do 120 j = 1, nx
c           do 115 i = 1, 2
c           cu2(i,j,k,l) = cu(i,j,k,ll-1)
c 115       continue
c 120       continue
c 130       continue
c        endif
c     endif
c ny+1 point is special
c     if (kyp2*(l+ks).eq.ny) then
c        do 136 j = 1, nx+1
c        do 135 i = 1, 2
c        cu2(i,j,1,l) = cu(i,j,kyp+1,kyb)
c 135    continue
c 136    continue
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(cu2(1,1,1,l),2*nxv,mreal,ll,moff+2,lgrp,lsid,
     1ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),2*kyp*nxv,mreal,ll-1,moff+2,
     1lgrp,msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(cu2(1,1,2,l),2*(kyp-1)*nxv,mreal,ll-2,moff+2,
     1lgrp,nsid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               call MPI_SEND(cu(1,1,1,l),2*nxv,mreal,lm,moff+2,lgrp,ierr
     1)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(cu(1,1,2,l),2*(kyp-1)*nxv,mreal,lm-1,moff+2
     1,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,l),2*kyp*nxv,mreal,lm-1,moff+2,lgrp,i
     1err)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 100 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 90 j = 1, nxv
            do 85 i = 1, 2
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   85       continue
   90       continue
  100       continue
         endif
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 k = 3, kyp
            k1 = kyp - k + 3
            k2 = k1/2 + 1
            joff = nxv*(k1 - 2*k2 + 2)
            do 110 j = 1, nxv
            do 105 i = 1, 2
            cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
  105       continue
  110       continue
  120       continue
         endif
      endif
c ny+1 point is special
      if (kyp2*(l+ks).eq.ny) then
         call MPI_IRECV(cu2(1,1,1,l),2*nxv,mreal,kyb-1,moff+3,lgrp,msid,
     1ierr)
      endif
      if (kyp*(l+ks+1).eq.ny) then
         kk = ny/kyp2
         call MPI_SEND(cu(1,1,kyp+1,l),2*nxv,mreal,kk,moff+3,lgrp,ierr)
      endif
      if (kyp2*(l+ks).eq.ny) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
  140 continue
c create even array
      do 200 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 190 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 150 j = 1, nxs
         cu2(1,j+1,k,l) = cu2(1,j+1,k,l)
         cu2(2,j+1,k,l) = 0.
         cu2(1,j+nx+1,k,l) = -cu2(1,nx-j+1,k,l)
         cu2(2,j+nx+1,k,l) = 0.
  150    continue
         cu2(1,1,k,l) = 0.
         cu2(2,1,k,l) = 0.
         cu2(1,nx+1,k,l) = 0.
         cu2(2,nx+1,k,l) = 0.
      else if (kk.le.ny) then
         do 160 j = 1, nxs
         cu2(1,nx+j+1,k,l) = -cu2(1,nx-j+1,k,l)
         cu2(2,nx+j+1,k,l) = cu2(2,nx-j+1,k,l)
  160    continue
         cu2(1,1,k,l) = 0.
         cu2(2,1,k,l) = cu2(2,1,k,l)
         cu2(1,nx+1,k,l) = 0.
         cu2(2,nx+1,k,l) = cu2(2,nx+1,k,l)
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 170 j = 1, nxs
            cu2(1,nx+j+1,k,l) = -cu2(1,nx-j+1,k,l)
            cu2(2,nx+j+1,k,l) = -cu2(2,nx-j+1,k,l)
  170       continue
            cu2(1,1,k,l) = 0.
            cu2(2,1,k,l) = -cu2(2,1,k,l)
            cu2(1,nx+1,k,l) = 0.
            cu2(2,nx+1,k,l) = -cu2(2,nx+1,k,l)
         else
            do 180 j = 1, nxs
            cu2(1,nx+j+1,kyp2-k+2,l) = -cu2(1,nx-j+1,k,l)
            cu2(2,nx+j+1,kyp2-k+2,l) = -cu2(2,nx-j+1,k,l)
  180       continue
            if (k.le.(kyp2/2+1)) then
               at1 = -cu2(2,1,kyp2-k+2,l)
               at2 = -cu2(2,nx+1,kyp2-k+2,l)
               cu2(2,1,kyp2-k+2,l) = -cu2(2,1,k,l)
               cu2(2,nx+1,kyp2-k+2,l) = -cu2(2,nx+1,k,l)
               cu2(2,1,k,l) = at1
               cu2(2,nx+1,k,l) = at2
            endif
            cu2(1,1,kyp2-k+2,l) = 0.
            cu2(1,nx+1,kyp2-k+2,l) = 0.
         endif
      endif
  190 continue
  200 continue
c finish even array
      do 230 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 220 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 210 j = 1, nxs
         cu2(1,nx-j+1,k,l) = -cu2(1,nx+j+1,k,l)
         cu2(2,nx-j+1,k,l) = cu2(2,nx+j+1,k,l)
  210    continue
      endif
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLCOS2D(q,q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2bl
     1ok)
c this subroutine creates an even array q2 from an array q, so that
c a 2d cosine transform can be performed with a 2d real to complex fft.
c linear interpolation for distributed data
c q2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = first dimension of input array q, must be >= nx+1
c kyp = number of data values per block in y
c kypd = second dimension of input array q, must be >= kyp+1
c kyp2 = second dimension of output array q2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real q, q2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension q(nxv,kypd,kblok), q2(2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer j, k, l, nxs, nys, ny2, kyb, kyb2, ks, koff, moff
      integer kk, ll, lm, k1, k2, joff
      real at1, at2
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      ks = kstrt - 2
      moff = kypd + kyb
c copy to double array in x direction
      do 80 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = koff/kyp + 1
      koff = kyp*(l + ks)
      lm = koff/kyp2 + 1
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, nys
         do 10 j = 1, nxs
         q2(j+1,k+1,l) = q(j+1,k+1,l)
         q2(nx+j+1,k+1,l) = q(nx-j+1,k+1,l)
         q2(j+1,ny+k+1,l) = q(j+1,ny-k+1,l)
         q2(nx+j+1,ny+k+1,l) = q(nx-j+1,ny-k+1,l)
   10    continue
         q2(1,k+1,l) = q(1,k+1,l)
         q2(nx+1,k+1,l) = q(nx+1,k+1,l)
         q2(1,k+ny+1,l) = q(1,ny-k+1,l)
         q2(nx+1,k+ny+1,l) = q(nx+1,ny-k+1,l)
   20    continue
         do 30 j = 1, nxs
         q2(j+1,1,l) = q(j+1,1,l)
         q2(j+nx+1,1,l) = q(nx-j+1,1,l)
         q2(j+1,ny+1,l) = q(j+1,ny+1,l)
         q2(j+nx+1,ny+1,l) = q(nx-j+1,ny+1,l)
   30    continue
         q2(1,1,l) = q(1,1,l)
         q2(nx+1,1,l) = q(nx+1,1,l)
         q2(1,ny+1,l) = q(1,ny+1,l)
         q2(nx+1,ny+1,l) = q(nx+1,ny+1,l)
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 50 k = 1, kyp
c        do 40 j = 1, nx
c        q2(j,k,l) = q(j,k,ll)
c  40    continue
c  50    continue
c        if (kyp.lt.kyp2) then
c           do 70 k = 1, kyp
c           do 60 j = 1, nx
c           q2(j,k+kyp,l) = q(j,k,ll+1)
c  60       continue
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(q2(1,1,l),kyp*nxv,mreal,ll-1,moff+1,lgrp,msid,ie
     1rr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(q2(1,kyp+1,l),kyp*nxv,mreal,ll,moff+1,lgrp,ns
     1id,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(q(1,1,l),kyp*nxv,mreal,lm-1,moff+1,lgrp,ierr)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 50 k = 2, kyp
         k1 = kyp - k + 2
         k2 = (k1 - 1)/2 + 1
         joff = nxv*(k1 - 2*k2 + 1)
         do 40 j = 1, nxv
         q2(j,k1,l) = q2(j+joff,k2,l)
   40    continue
   50    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 70 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 60 j = 1, nxv
            q2(j,k1+kyp,l) = q2(j+joff,k2+kyp,l)
   60       continue
   70       continue
         endif
      endif
   80 continue
c copy to double array in y direction
      do 140 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(l + ks)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 90 j = 1, nx
c           q2(j,1,l) = q(j,1,ll+1)
c  90       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 110 k = 1, kyp
c           do 100 j = 1, nx
c           q2(j,k+kyp,l) = q(j,k,ll)
c 100       continue
c 110       continue
c        endif
c        if (kyp.gt.1) then
c           do 130 k = 2, kyp
c           do 120 j = 1, nx
c           q2(j,k,l) = q(j,k,ll-1)
c 120       continue
c 130       continue
c        endif
c     endif
c ny+1 point is special
c     if (kyp2*(l+ks).eq.ny) then
c        do 135 j = 1, nx+1
c        q2(j,1,l) = q(j,kyp+1,kyb)
c 135    continue
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(q2(1,1,l),nxv,mreal,ll,moff+2,lgrp,lsid,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(q2(1,kyp+1,l),kyp*nxv,mreal,ll-1,moff+2,lgrp,
     1msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(q2(1,2,l),(kyp-1)*nxv,mreal,ll-2,moff+2,lgrp,
     1nsid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               call MPI_SEND(q(1,1,l),nxv,mreal,lm,moff+2,lgrp,ierr)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(q(1,2,l),(kyp-1)*nxv,mreal,lm-1,moff+2,lgrp
     1,ierr)
            endif
         else
            call MPI_SEND(q(1,1,l),kyp*nxv,mreal,lm-1,moff+2,lgrp,ierr)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 100 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 90 j = 1, nxv
            q2(j,k1+kyp,l) = q2(j+joff,k2+kyp,l)
   90       continue
  100       continue
         endif
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 k = 3, kyp
            k1 = kyp - k + 3
            k2 = k1/2 + 1
            joff = nxv*(k1 - 2*k2 + 2)
            do 110 j = 1, nxv
            q2(j,k1,l) = q2(j+joff,k2,l)
  110       continue
  120       continue
         endif
      endif
c ny+1 point is special
      if (kyp2*(l+ks).eq.ny) then
         call MPI_IRECV(q2(1,1,l),nxv,mreal,kyb-1,moff+3,lgrp,msid,ierr)
      endif
      if (kyp*(l+ks+1).eq.ny) then
         kk = ny/kyp2
         call MPI_SEND(q(1,kyp+1,l),nxv,mreal,kk,moff+3,lgrp,ierr)
      endif
      if (kyp2*(l+ks).eq.ny) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
  140 continue
c create even array
      do 200 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 190 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 150 j = 1, nxs
         q2(j+1,k,l) = q2(j+1,k,l)
         q2(j+nx+1,k,l) = q2(nx-j+1,k,l)
  150    continue
         q2(1,k,l) = q2(1,k,l)
         q2(nx+1,k,l) = q2(nx+1,k,l)
      else if (kk.le.ny) then
         do 160 j = 1, nxs
         q2(nx+j+1,k,l) = q2(nx-j+1,k,l)
  160    continue
         q2(1,k,l) = q2(1,k,l)
         q2(nx+1,k,l) = q2(nx+1,k,l)
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 170 j = 1, nxs
            q2(nx+j+1,k,l) = q2(nx-j+1,k,l)
  170       continue
            q2(1,k,l) = q2(1,k,l)
            q2(nx+1,k,l) = q2(nx+1,k,l)
         else
            do 180 j = 1, nxs
            q2(nx+j+1,kyp2-k+2,l) = q2(nx-j+1,k,l)
  180       continue
            if (k.le.(kyp2/2+1)) then
               at1 = q2(1,kyp2-k+2,l)
               at2 = q2(nx+1,kyp2-k+2,l)
               q2(1,kyp2-k+2,l) = q2(1,k,l)
               q2(nx+1,kyp2-k+2,l) = q2(nx+1,k,l)
               q2(1,k,l) = at1
               q2(nx+1,k,l) = at2
            endif
         endif
      endif
  190 continue
  200 continue
c finish even array
      do 230 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 220 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 210 j = 1, nxs
         q2(nx-j+1,k,l) = q2(nx+j+1,k,l)
  210    continue
      endif
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLCOS2B(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2
     1blok)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d cosine/sine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in x,
c y component is an odd function in x, and the z component is an even
c function in both x and y.
c linear interpolation for distributed data
c cu2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = first dimension of input array q, must be >= nx+1
c kyp = number of data values per block in y
c kypd = second dimension of input array q, must be >= kyp+1
c kyp2 = second dimension of output array q2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real cu, cu2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension cu(3,nxv,kypd,kblok), cu2(3,2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, nxs, nys, ny2, kyb, kyb2, ks, koff, moff
      integer kk, ll, lm, k1, k2, joff
      real at1, at2
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      ks = kstrt - 2
      moff = kypd + kyb
c copy to double array in x direction
      do 80 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = koff/kyp + 1
      koff = kyp*(l + ks)
      lm = koff/kyp2 + 1
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, nys
         do 10 j = 1, nxs
         cu2(1,j+1,k+1,l) = cu(1,j+1,k+1,l)
         cu2(2,j+1,k+1,l) = cu(2,j+1,k+1,l)
         cu2(3,j+1,k+1,l) = cu(3,j+1,k+1,l)
         cu2(1,nx+j+1,k+1,l) = -cu(1,nx-j+1,k+1,l)
         cu2(2,nx+j+1,k+1,l) = cu(2,nx-j+1,k+1,l)
         cu2(3,nx+j+1,k+1,l) = cu(3,nx-j+1,k+1,l)
         cu2(1,j+1,ny+k+1,l) = cu(1,j+1,ny-k+1,l)
         cu2(2,j+1,ny+k+1,l) = -cu(2,j+1,ny-k+1,l)
         cu2(3,j+1,ny+k+1,l) = cu(3,j+1,ny-k+1,l)
         cu2(1,nx+j+1,ny+k+1,l) = -cu(1,nx-j+1,ny-k+1,l)
         cu2(2,nx+j+1,ny+k+1,l) = -cu(2,nx-j+1,ny-k+1,l)
         cu2(3,nx+j+1,ny+k+1,l) = cu(3,nx-j+1,ny-k+1,l)
   10    continue
         cu2(1,1,k+1,l) = 0.
         cu2(2,1,k+1,l) = cu(2,1,k+1,l)
         cu2(3,1,k+1,l) = cu(3,1,k+1,l)
         cu2(1,nx+1,k+1,l) = 0.
         cu2(2,nx+1,k+1,l) = cu(2,nx+1,k+1,l)
         cu2(3,nx+1,k+1,l) = cu(3,nx+1,k+1,l)
         cu2(1,1,k+ny+1,l) = 0.
         cu2(2,1,k+ny+1,l) = -cu(2,1,ny-k+1,l)
         cu2(3,1,k+ny+1,l) = cu(3,1,ny-k+1,l)
         cu2(1,nx+1,k+ny+1,l) = 0.
         cu2(2,nx+1,k+ny+1,l) = -cu(2,nx+1,ny-k+1,l)
         cu2(3,nx+1,k+ny+1,l) = cu(3,nx+1,ny-k+1,l)
   20    continue
         do 30 j = 1, nxs
         cu2(1,j+1,1,l) = cu(1,j+1,1,l)
         cu2(2,j+1,1,l) = 0.
         cu2(3,j+1,1,l) = cu(3,j+1,1,l)
         cu2(1,j+nx+1,1,l) = -cu(1,nx-j+1,1,l)
         cu2(2,j+nx+1,1,l) = 0.
         cu2(3,j+nx+1,1,l) = cu(3,nx-j+1,1,l)
         cu2(1,j+1,ny+1,l) = cu(1,j+1,ny+1,l)
         cu2(2,j+1,ny+1,l) = 0.
         cu2(3,j+1,ny+1,l) = cu(3,j+1,ny+1,l)
         cu2(1,j+nx+1,ny+1,l) = -cu(1,nx-j+1,ny+1,l)
         cu2(2,j+nx+1,ny+1,l) = 0.
         cu2(3,j+nx+1,ny+1,l) = cu(3,nx-j+1,ny+1,l)
   30    continue
         cu2(1,1,1,l) = 0.
         cu2(2,1,1,l) = 0.
         cu2(3,1,1,l) = cu(3,1,1,l)
         cu2(1,nx+1,1,l) = 0.
         cu2(2,nx+1,1,l) = 0.
         cu2(3,nx+1,1,l) = cu(3,nx+1,1,l)
         cu2(1,1,ny+1,l) = 0.
         cu2(2,1,ny+1,l) = 0.
         cu2(3,1,ny+1,l) = cu(3,1,ny+1,l)
         cu2(1,nx+1,ny+1,l) = 0.
         cu2(2,nx+1,ny+1,l) = 0.
         cu2(3,nx+1,ny+1,l) = cu(3,nx+1,ny+1,l)
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 50 k = 1, kyp
c        do 40 j = 1, nx
c        do 35 i = 1, 3
c        cu2(i,j,k,l) = cu(i,j,k,ll)
c  35    continue
c  40    continue
c  50    continue
c        if (kyp.lt.kyp2) then
c           do 70 k = 1, kyp
c           do 60 j = 1, nx
c           do 55 i = 1, 3
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll+1)
c  55       continue
c  60       continue
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(cu2(1,1,1,l),3*kyp*nxv,mreal,ll-1,moff+1,lgrp,ms
     1id,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),3*kyp*nxv,mreal,ll,moff+1,lg
     1rp,nsid,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(cu(1,1,1,l),3*kyp*nxv,mreal,lm-1,moff+1,lgrp,ierr
     1)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 50 k = 2, kyp
         k1 = kyp - k + 2
         k2 = (k1 - 1)/2 + 1
         joff = nxv*(k1 - 2*k2 + 1)
         do 40 j = 1, nxv
         do 35 i = 1, 3
         cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
   35    continue
   40    continue
   50    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 70 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 60 j = 1, nxv
            do 55 i = 1, 3
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   55       continue
   60       continue
   70       continue
         endif
      endif
   80 continue
c copy to double array in y direction
      do 140 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(l + ks)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 90 j = 1, nx
c           do 85 i = 1, 3
c           cu2(i,j,1,l) = cu(i,j,1,ll+1)
c  85       continue
c  90       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 110 k = 1, kyp
c           do 100 j = 1, nx
c           do 95 i = 1, 3
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll)
c  95       continue
c 100       continue
c 110       continue
c        endif
c        if (kyp.gt.1) then
c           do 130 k = 2, kyp
c           do 120 j = 1, nx
c           do 115 i = 1, 3
c           cu2(i,j,k,l) = cu(i,j,k,ll-1)
c 115       continue
c 120       continue
c 130       continue
c        endif
c     endif
c ny+1 point is special
c     if (kyp2*(l+ks).eq.ny) then
c        do 136 j = 1, nx+1
c        do 135 i = 1, 3
c        cu2(i,j,1,l) = cu(i,j,kyp+1,kyb)
c 135    continue
c 136    continue
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(cu2(1,1,1,l),3*nxv,mreal,ll,moff+2,lgrp,lsid,
     1ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),3*kyp*nxv,mreal,ll-1,moff+2,
     1lgrp,msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(cu2(1,1,2,l),3*(kyp-1)*nxv,mreal,ll-2,moff+2,
     1lgrp,nsid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               call MPI_SEND(cu(1,1,1,l),3*nxv,mreal,lm,moff+2,lgrp,ierr
     1)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(cu(1,1,2,l),3*(kyp-1)*nxv,mreal,lm-1,moff+2
     1,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,l),3*kyp*nxv,mreal,lm-1,moff+2,lgrp,i
     1err)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 100 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 90 j = 1, nxv
            do 85 i = 1, 3
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   85       continue
   90       continue
  100       continue
         endif
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 k = 3, kyp
            k1 = kyp - k + 3
            k2 = k1/2 + 1
            joff = nxv*(k1 - 2*k2 + 2)
            do 110 j = 1, nxv
            do 105 i = 1, 3
            cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
  105       continue
  110       continue
  120       continue
         endif
      endif
c ny+1 point is special
      if (kyp2*(l+ks).eq.ny) then
         call MPI_IRECV(cu2(1,1,1,l),3*nxv,mreal,kyb-1,moff+3,lgrp,msid,
     1ierr)
      endif
      if (kyp*(l+ks+1).eq.ny) then
         kk = ny/kyp2
         call MPI_SEND(cu(1,1,kyp+1,l),3*nxv,mreal,kk,moff+3,lgrp,ierr)
      endif
      if (kyp2*(l+ks).eq.ny) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
  140 continue
c create even array
      do 200 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 190 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 150 j = 1, nxs
         cu2(1,j+1,k,l) = cu2(1,j+1,k,l)
         cu2(2,j+1,k,l) = 0.
         cu2(3,j+1,k,l) = cu2(3,j+1,k,l)
         cu2(1,j+nx+1,k,l) = -cu2(1,nx-j+1,k,l)
         cu2(2,j+nx+1,k,l) = 0.
         cu2(3,j+nx+1,k,l) = cu2(3,nx-j+1,k,l)
  150    continue
         cu2(1,1,k,l) = 0.
         cu2(2,1,k,l) = 0.
         cu2(3,1,k,l) = cu2(3,1,k,l)
         cu2(1,nx+1,k,l) = 0.
         cu2(2,nx+1,k,l) = 0.
         cu2(3,nx+1,k,l) = cu2(3,nx+1,k,l)
      else if (kk.le.ny) then
         do 160 j = 1, nxs
         cu2(1,nx+j+1,k,l) = -cu2(1,nx-j+1,k,l)
         cu2(2,nx+j+1,k,l) = cu2(2,nx-j+1,k,l)
         cu2(3,nx+j+1,k,l) = cu2(3,nx-j+1,k,l)
  160    continue
         cu2(1,1,k,l) = 0.
         cu2(2,1,k,l) = cu2(2,1,k,l)
         cu2(3,1,k,l) = cu2(3,1,k,l)
         cu2(1,nx+1,k,l) = 0.
         cu2(2,nx+1,k,l) = cu2(2,nx+1,k,l)
         cu2(3,nx+1,k,l) = cu2(3,nx+1,k,l)
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 170 j = 1, nxs
            cu2(1,nx+j+1,k,l) = -cu2(1,nx-j+1,k,l)
            cu2(2,nx+j+1,k,l) = -cu2(2,nx-j+1,k,l)
            cu2(3,nx+j+1,k,l) = cu2(3,nx-j+1,k,l)
  170       continue
            cu2(1,1,k,l) = 0.
            cu2(2,1,k,l) = -cu2(2,1,k,l)
            cu2(3,1,k,l) = cu2(3,1,k,l)
            cu2(1,nx+1,k,l) = 0.
            cu2(2,nx+1,k,l) = -cu2(2,nx+1,k,l)
            cu2(3,nx+1,k,l) = cu2(3,nx+1,k,l)
         else
            do 180 j = 1, nxs
            cu2(1,nx+j+1,kyp2-k+2,l) = -cu2(1,nx-j+1,k,l)
            cu2(2,nx+j+1,kyp2-k+2,l) = -cu2(2,nx-j+1,k,l)
            cu2(3,nx+j+1,kyp2-k+2,l) = cu2(3,nx-j+1,k,l)
  180       continue
            if (k.le.(kyp2/2+1)) then
               at1 = -cu2(2,1,kyp2-k+2,l)
               at2 = -cu2(2,nx+1,kyp2-k+2,l)
               cu2(2,1,kyp2-k+2,l) = -cu2(2,1,k,l)
               cu2(2,nx+1,kyp2-k+2,l) = -cu2(2,nx+1,k,l)
               cu2(2,1,k,l) = at1
               cu2(2,nx+1,k,l) = at2
               at1 = cu2(3,1,kyp2-k+2,l)
               at2 = cu2(3,nx+1,kyp2-k+2,l)
               cu2(3,1,kyp2-k+2,l) = cu2(3,1,k,l)
               cu2(3,nx+1,kyp2-k+2,l) = cu2(3,nx+1,k,l)
               cu2(3,1,k,l) = at1
               cu2(3,nx+1,k,l) = at2
            endif
            cu2(1,1,kyp2-k+2,l) = 0.
            cu2(1,nx+1,kyp2-k+2,l) = 0.
         endif
      endif
  190 continue
  200 continue
c finish even array
      do 230 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 220 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 210 j = 1, nxs
         cu2(1,nx-j+1,k,l) = -cu2(1,nx+j+1,k,l)
         cu2(2,nx-j+1,k,l) = cu2(2,nx+j+1,k,l)
         cu2(3,nx-j+1,k,l) = cu2(3,nx+j+1,k,l)
  210    continue
      endif
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFDBL2C(fxy,fxy2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,
     1k2blok)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c for distributed data
c fxy array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = second dimension of output array fxy, must be >= nx
c kyp = number of data values per block in y
c kypd = third dimension of output array fxy, must be >= kyp+1
c kyp2 = third dimension of output array fxy2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real fxy, fxy2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension fxy(2,nxv,kypd,kblok), fxy2(2,2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer j, k, l, nx1, ny1, kyb, kyb2, kyp1, ks, joff, koff, moff
      integer kk, ll, lm
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      kyp1 = kyp + 1
      ks = kstrt - 2
      moff = kypd + kyb
      do 90 l = 1, k2blok
      koff = kyp2*(l + ks)
      lm = koff/kyp + 1
      koff = kyp*(l + ks)
      ll = koff/kyp2 + 1
      koff = koff - kyp2*(ll - 1)
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, ny1
         do 10 j = 1, nx1
         fxy(1,j,k,l) = fxy2(1,j,k,l)
         fxy(2,j,k,l) = fxy2(2,j,k,l)
   10    continue
   20    continue
         go to 90
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 40 k = 1, kyp1
c           do 30 j = 1, nx1
c           fxy(1,j,k,l) = fxy2(1,j,k,ll)
c           fxy(2,j,k,l) = fxy2(2,j,k,ll)
c  30       continue
c  40       continue
c        else
c           do 60 k = 1, kyp
c           do 50 j = 1, nx1
c           fxy(1,j,k,l) = fxy2(1,j,k+koff,ll)
c           fxy(2,j,k,l) = fxy2(2,j,k+koff,ll)
c  50       continue
c  60       continue
c           do 70 j = 1, nx1
c           fxy(1,j,kyp+1,l) = fxy2(1,j,1,ll+1)
c           fxy(2,j,kyp+1,l) = fxy2(2,j,1,ll+1)
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(fxy(1,1,1,l),2*nxv*kyp1,mreal,ll-1,moff+3,lgr
     1p,msid,ierr)
         else
            call MPI_IRECV(fxy(1,1,1,l),2*nxv*kyp,mreal,ll-1,moff+3,lgrp
     1,msid,ierr)
            call MPI_IRECV(fxy(1,1,kyp+1,l),2*nxv,mreal,ll,moff+3,lgrp,n
     1sid,ierr)
         endif
      endif
c pack data and send it
      if (lm.le.kyb) then
         if (kyp.lt.kyp2) then
            do 40 k = 2, kyp1
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do 30 j = 1, nxv
            fxy2(1,j+joff,kk,l) = fxy2(1,j,k,l)
            fxy2(2,j+joff,kk,l) = fxy2(2,j,k,l)
   30       continue
   40       continue
            call MPI_SEND(fxy2(1,1,1,l),2*nxv*kyp1,mreal,lm-1,moff+3,lgr
     1p,ierr)
            do 60 k = 2, kyp
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do 50 j = 1, nxv
            fxy2(1,j+joff,kk+kyp,l) = fxy2(1,j,k+kyp,l)
            fxy2(2,j+joff,kk+kyp,l) = fxy2(2,j,k+kyp,l)
   50       continue
   60       continue
            call MPI_SEND(fxy2(1,1,kyp+1,l),2*nxv*kyp,mreal,lm,moff+3,lg
     1rp,ierr)
         else
            do 80 k = 2, kyp
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do 70 j = 1, nxv
            fxy2(1,j+joff,kk,l) = fxy2(1,j,k,l)
            fxy2(2,j+joff,kk,l) = fxy2(2,j,k,l)
   70       continue
   80       continue
            call MPI_SEND(fxy2(1,1,1,l),2*nxv*kyp,mreal,lm-1,moff+3,lgrp
     1,ierr)
         endif
         if (lm.gt.1) then
            call MPI_SEND(fxy2(1,1,1,l),2*nxv,mreal,lm-2,moff+3,lgrp,ier
     1r)
         endif
      else if (lm.eq.(kyb+1)) then
         call MPI_SEND(fxy2(1,1,1,l),2*nxv,mreal,lm-2,moff+3,lgrp,ierr)
      endif
c wait for data
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFDBL2D(q,q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2bl
     1ok)
c this subroutine copies data from a double array to regular array
c with guard cells for scalar field and linear interpolation
c for distributed data
c q2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = second dimension of output array fxy, must be >= nx
c kyp = number of data values per block in y
c kypd = third dimension of output array fxy, must be >= kyp+1
c kyp2 = third dimension of output array fxy2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real q, q2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension q(nxv,kypd,kblok), q2(2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer j, k, l, nx1, ny1, kyb, kyb2, kyp1, ks, joff, koff, moff
      integer  kk, ll, lm
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      kyp1 = kyp + 1
      ks = kstrt - 2
      moff = kypd + kyb
      do 90 l = 1, k2blok
      koff = kyp2*(l + ks)
      lm = koff/kyp + 1
      koff = kyp*(l + ks)
      ll = koff/kyp2 + 1
      koff = koff - kyp2*(ll - 1)
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, ny1
         do 10 j = 1, nx1
         q(j,k,l) = q2(j,k,l)
   10    continue
   20    continue
         go to 90
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 40 k = 1, kyp1
c           do 30 j = 1, nx1
c           q(j,k,l) = q2(j,k,ll)
c  30       continue
c  40       continue
c        else
c           do 60 k = 1, kyp
c           do 50 j = 1, nx1
c           q(j,k,l) = q2(j,k+koff,ll)
c  50       continue
c  60       continue
c           do 70 j = 1, nx1
c           q(j,kyp+1,l) = q2(j,1,ll+1)
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(q(1,1,l),nxv*kyp1,mreal,ll-1,moff+4,lgrp,msid
     1,ierr)
         else
            call MPI_IRECV(q(1,1,l),nxv*kyp,mreal,ll-1,moff+4,lgrp,msid,
     1ierr)
            call MPI_IRECV(q(1,kyp+1,l),nxv,mreal,ll,moff+4,lgrp,nsid,ie
     1rr)
         endif
      endif
c pack data and send it
      if (lm.le.kyb) then
         if (kyp.lt.kyp2) then
            do 40 k = 2, kyp1
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do 30 j = 1, nxv
            q2(j+joff,kk,l) = q2(j,k,l)
   30       continue
   40       continue
            call MPI_SEND(q2(1,1,l),nxv*kyp1,mreal,lm-1,moff+4,lgrp,ierr
     1)
            do 60 k = 2, kyp
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do 50 j = 1, nxv
            q2(j+joff,kk+kyp,l) = q2(j,k+kyp,l)
   50       continue
   60       continue
            call MPI_SEND(q2(1,kyp+1,l),nxv*kyp,mreal,lm,moff+4,lgrp,ier
     1r)
         else
            do 80 k = 2, kyp
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do 70 j = 1, nxv
            q2(j+joff,kk,l) = q2(j,k,l)
   70       continue
   80       continue
            call MPI_SEND(q2(1,1,l),nxv*kyp,mreal,lm-1,moff+4,lgrp,ierr)
         endif
         if (lm.gt.1) then
            call MPI_SEND(q2(1,1,l),nxv,mreal,lm-2,moff+4,lgrp,ierr)
         endif
      else if (lm.eq.(kyb+1)) then
         call MPI_SEND(q2(1,1,l),nxv,mreal,lm-2,moff+4,lgrp,ierr)
      endif
c wait for data
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFDBL2B(fxy,fxy2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,
     1k2blok)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c for distributed data
c fxy array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = second dimension of output array fxy, must be >= nx
c kyp = number of data values per block in y
c kypd = third dimension of output array fxy, must be >= kyp+1
c kyp2 = third dimension of output array fxy2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real fxy, fxy2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension fxy(3,nxv,kypd,kblok), fxy2(3,2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer j, k, l, nx1, ny1, kyb, kyb2, kyp1, ks, joff, koff, moff
      integer kk, ll, lm
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      kyp1 = kyp + 1
      ks = kstrt - 2
      moff = kypd + kyb
      do 90 l = 1, k2blok
      koff = kyp2*(l + ks)
      lm = koff/kyp + 1
      koff = kyp*(l + ks)
      ll = koff/kyp2 + 1
      koff = koff - kyp2*(ll - 1)
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, ny1
         do 10 j = 1, nx1
         fxy(1,j,k,l) = fxy2(1,j,k,l)
         fxy(2,j,k,l) = fxy2(2,j,k,l)
         fxy(3,j,k,l) = fxy2(3,j,k,l)
   10    continue
   20    continue
         go to 90
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 40 k = 1, kyp1
c           do 30 j = 1, nx1
c           fxy(1,j,k,l) = fxy2(1,j,k,ll)
c           fxy(2,j,k,l) = fxy2(2,j,k,ll)
c           fxy(3,j,k,l) = fxy2(3,j,k,ll)
c  30       continue
c  40       continue
c        else
c           do 60 k = 1, kyp
c           do 50 j = 1, nx1
c           fxy(1,j,k,l) = fxy2(1,j,k+koff,ll)
c           fxy(2,j,k,l) = fxy2(2,j,k+koff,ll)
c           fxy(3,j,k,l) = fxy2(3,j,k+koff,ll)
c  50       continue
c  60       continue
c           do 70 j = 1, nx1
c           fxy(1,j,kyp+1,l) = fxy2(1,j,1,ll+1)
c           fxy(2,j,kyp+1,l) = fxy2(2,j,1,ll+1)
c           fxy(3,j,kyp+1,l) = fxy2(3,j,1,ll+1)
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(fxy(1,1,1,l),3*nxv*kyp1,mreal,ll-1,moff+3,lgr
     1p,msid,ierr)
         else
            call MPI_IRECV(fxy(1,1,1,l),3*nxv*kyp,mreal,ll-1,moff+3,lgrp
     1,msid,ierr)
            call MPI_IRECV(fxy(1,1,kyp+1,l),3*nxv,mreal,ll,moff+3,lgrp,n
     1sid,ierr)
         endif
      endif
c pack data and send it
      if (lm.le.kyb) then
         if (kyp.lt.kyp2) then
            do 40 k = 2, kyp1
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do 30 j = 1, nxv
            fxy2(1,j+joff,kk,l) = fxy2(1,j,k,l)
            fxy2(2,j+joff,kk,l) = fxy2(2,j,k,l)
            fxy2(3,j+joff,kk,l) = fxy2(3,j,k,l)
   30       continue
   40       continue
            call MPI_SEND(fxy2(1,1,1,l),3*nxv*kyp1,mreal,lm-1,moff+3,lgr
     1p,ierr)
            do 60 k = 2, kyp
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do 50 j = 1, nxv
            fxy2(1,j+joff,kk+kyp,l) = fxy2(1,j,k+kyp,l)
            fxy2(2,j+joff,kk+kyp,l) = fxy2(2,j,k+kyp,l)
            fxy2(3,j+joff,kk+kyp,l) = fxy2(3,j,k+kyp,l)
   50       continue
   60       continue
            call MPI_SEND(fxy2(1,1,kyp+1,l),3*nxv*kyp,mreal,lm,moff+3,lg
     1rp,ierr)
         else
            do 80 k = 2, kyp
            kk = (k - 1)/2 + 1
            joff = nxv*(k - 2*kk + 1)
            do 70 j = 1, nxv
            fxy2(1,j+joff,kk,l) = fxy2(1,j,k,l)
            fxy2(2,j+joff,kk,l) = fxy2(2,j,k,l)
            fxy2(3,j+joff,kk,l) = fxy2(3,j,k,l)
   70       continue
   80       continue
            call MPI_SEND(fxy2(1,1,1,l),3*nxv*kyp,mreal,lm-1,moff+3,lgrp
     1,ierr)
         endif
         if (lm.gt.1) then
            call MPI_SEND(fxy2(1,1,1,l),3*nxv,mreal,lm-2,moff+3,lgrp,ier
     1r)
         endif
      else if (lm.eq.(kyb+1)) then
         call MPI_SEND(fxy2(1,1,1,l),3*nxv,mreal,lm-2,moff+3,lgrp,ierr)
      endif
c wait for data
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLCGUARD2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.  for vector data
c the field is replicated so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes three extra guard cells.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(2,nxv,nypmx,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx2, ks, moff, kr, krr, kl, kll, ngc, j, l
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 70 l = 1, kblok
      kr = l + ks + 2
      krr = kr
      kl = l + ks
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         kll = kll - 1
         ngc = 1
      endif
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nxv
c        f(1,j,1,l) = f(1,j,kyp+1,kl)
c        f(2,j,1,l) = f(2,j,kyp+1,kl)
c  10    continue
c     else
c        do 20 j = 1, nxv
c        f(1,j,1,l) = f(1,j,3,l)
c        f(2,j,1,l) = f(2,j,3,l)
c  20    continue
c     endif
c     if (kr.le.nvp) then
c        do 30 j = 1, nxv
c        f(1,j,kyp+2,l) = f(1,j,2,kr)
c        f(2,j,kyp+2,l) = f(2,j,2,kr)
c        f(1,j,kyp+3,l) = f(1,j,ngc+1,krr)
c        f(2,j,kyp+3,l) = f(2,j,ngc+1,krr)
c  30    continue
c     else
c        do 40 j = 2, nx2
c        f(1,j,kyp+3,l) = 2.*f(1,j,kyp+2,l) - f(1,j,kyp+1,l)
c        f(2,j,kyp+3,l) = 2.*f(2,j,kyp+2,l) - f(2,j,kyp+1,l)
c  40    continue
c     endif
c     if (kyp.eq.1) then
c        if ((kl.eq.0).and.(kr.le.nvp)) then
c           do 50 j = 1, nxv
c           f(1,j,1,l) = f(1,j,2,kr)
c           f(2,j,1,l) = f(2,j,2,kr)
c  50       continue
c        endif
c        if (kr.eq.nvp) then
c           do 60 j = 1, nxv
c           f(1,j,kyp+3,l) = f(1,j,kyp+2,kr)
c           f(2,j,kyp+3,l) = f(2,j,kyp+2,kr)
c  60       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kl-1,moff+3,lgrp,msid,ier
     1r)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,kyp+1,l),2*nxv,mreal,kr-1,moff+3,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nxv
         f(1,j,1,l) = f(1,j,3,l)
         f(2,j,1,l) = f(2,j,3,l)
   20    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(f(1,1,kyp+2,l),2*ngc*nxv,mreal,kr-1,moff+4,lgrp,
     1msid,ierr)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,2,l),2*ngc*nxv,mreal,kl-1,moff+4,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 40 j = 2, nx2
         f(1,j,kyp+3,l) = 2.*f(1,j,kyp+2,l) - f(1,j,kyp+1,l)
         f(2,j,kyp+3,l) = 2.*f(2,j,kyp+2,l) - f(2,j,kyp+1,l)
   40    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (krr.le.nvp) then
            call MPI_IRECV(f(1,1,kyp+3,l),2*nxv,mreal,krr-1,moff+6,lgrp,
     1msid,ierr)
         else if (kr.le.nvp) then
            call MPI_IRECV(f(1,1,kyp+3,l),2*nxv,mreal,kr-1,moff+6,lgrp,m
     1sid,ierr)
         endif
         if ((kl.eq.0).and.(kr.le.nvp)) then
            call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kr-1,moff+6,lgrp,nsid,
     1ierr)
         endif
         if (kll.ge.1) then
            call MPI_SEND(f(1,1,2,l),2*nxv,mreal,kll-1,moff+6,lgrp,ierr)
         else if (kl.eq.1) then
            call MPI_SEND(f(1,1,2,l),2*nxv,mreal,kl-1,moff+6,lgrp,ierr)
         endif
         if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
            call MPI_SEND(f(1,1,kyp+2,l),2*nxv,mreal,kl-1,moff+6,lgrp,ie
     1rr)
         endif
         if (kr.le.nvp) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
         if (kl.eq.0) then
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
   70 continue
c fix left edge
      do 100 l = 1, kblok
      kl = l + ks
      if (kl.eq.0) then
         do 90 j = 2, nx2
         f(1,j,1,l) = 2.*f(1,j,2,l) - f(1,j,1,l)
         f(2,j,1,l) = 2.*f(2,j,2,l) - f(2,j,1,l)
   90    continue
      endif
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLDGUARD2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is replicated so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes three extra guard cells.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(nxv,nypmx,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx2, ks, moff, kr, krr, kl, kll, ngc, j, l
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 70 l = 1, kblok
      kr = l + ks + 2
      krr = kr
      kl = l + ks
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         kll = kll - 1
         ngc = 1
      endif
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nxv
c        f(j,1,l) = f(j,kyp+1,kl)
c  10    continue
c     else
c        do 20 j = 1, nxv
c        f(j,1,l) = f(j,3,l)
c  20    continue
c     endif
c     if (kr.le.nvp) then
c        do 30 j = 1, nxv
c        f(j,kyp+2,l) = f(j,2,kr)
c        f(j,kyp+3,l) = f(j,ngc+1,krr)
c  30    continue
c     else
c        do 40 j = 2, nx2
c        f(j,kyp+3,l) = 2.*f(j,kyp+2,l) - f(j,kyp+1,l)
c  40    continue
c     endif
c     if (kyp.eq.1) then
c        if ((kl.eq.0).and.(kr.le.nvp)) then
c           do 50 j = 1, nxv
c           f(j,1,l) = f(j,2,kr)
c  50       continue
c        endif
c        if (kr.eq.nvp) then
c           do 60 j = 1, nxv
c           f(j,kyp+3,l) = f(j,kyp+2,kr)
c  60       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(f(1,1,l),nxv,mreal,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,kyp+1,l),nxv,mreal,kr-1,moff+3,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nxv
         f(j,1,l) = f(j,3,l)
   20    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(f(1,kyp+2,l),ngc*nxv,mreal,kr-1,moff+4,lgrp,msid
     1,ierr)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,2,l),ngc*nxv,mreal,kl-1,moff+4,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 40 j = 2, nx2
         f(j,kyp+3,l) = 2.*f(j,kyp+2,l) - f(j,kyp+1,l)
   40    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (krr.le.nvp) then
            call MPI_IRECV(f(1,kyp+3,l),nxv,mreal,krr-1,moff+6,lgrp,msid
     1,ierr)
         else if (kr.le.nvp) then
            call MPI_IRECV(f(1,kyp+3,l),nxv,mreal,kr-1,moff+6,lgrp,msid,
     1ierr)
         endif
         if ((kl.eq.0).and.(kr.le.nvp)) then
            call MPI_IRECV(f(1,1,l),nxv,mreal,kr-1,moff+6,lgrp,nsid,ierr
     1)
         endif
         if (kll.ge.1) then
            call MPI_SEND(f(1,2,l),nxv,mreal,kll-1,moff+6,lgrp,ierr)
         else if (kl.eq.1) then
            call MPI_SEND(f(1,2,l),nxv,mreal,kl-1,moff+6,lgrp,ierr)
         endif
         if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
            call MPI_SEND(f(1,kyp+2,l),nxv,mreal,kl-1,moff+6,lgrp,ierr)
         endif
         if (kr.le.nvp) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
         if (kl.eq.0) then
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
   70 continue
c fix left edge
      do 100 l = 1, kblok
      kl = l + ks
      if (kl.eq.0) then
         do 90 j = 2, nx2
         f(j,1,l) = 2.*f(j,2,l) - f(j,1,l)
   90    continue
      endif
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLBGUARD2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.  for vector data
c the field is replicated so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes three extra guard cells.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(3,nxv,nypmx,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx2, ks, moff, kr, krr, kl, kll, ngc, j, l
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 70 l = 1, kblok
      kr = l + ks + 2
      krr = kr
      kl = l + ks
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         kll = kll - 1
         ngc = 1
      endif
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nxv
c        f(1,j,1,l) = f(1,j,kyp+1,kl)
c        f(2,j,1,l) = f(2,j,kyp+1,kl)
c        f(3,j,1,l) = f(3,j,kyp+1,kl)
c  10    continue
c     else
c        do 20 j = 1, nxv
c        f(1,j,1,l) = f(1,j,3,l)
c        f(2,j,1,l) = f(2,j,3,l)
c        f(3,j,1,l) = f(3,j,3,l)
c  20    continue
c     endif
c     if (kr.le.nvp) then
c        do 30 j = 1, nxv
c        f(1,j,kyp+2,l) = f(1,j,2,kr)
c        f(2,j,kyp+2,l) = f(2,j,2,kr)
c        f(3,j,kyp+2,l) = f(3,j,2,kr)
c        f(1,j,kyp+3,l) = f(1,j,ngc+1,krr)
c        f(2,j,kyp+3,l) = f(2,j,ngc+1,krr)
c        f(3,j,kyp+3,l) = f(3,j,ngc+1,krr)
c  30    continue
c     else
c        do 40 j = 2, nx2
c        f(1,j,kyp+3,l) = 2.*f(1,j,kyp+2,l) - f(1,j,kyp+1,l)
c        f(2,j,kyp+3,l) = 2.*f(2,j,kyp+2,l) - f(2,j,kyp+1,l)
c        f(3,j,kyp+3,l) = 2.*f(3,j,kyp+2,l) - f(3,j,kyp+1,l)
c  40    continue
c     endif
c     if (kyp.eq.1) then
c        if ((kl.eq.0).and.(kr.le.nvp)) then
c           do 50 j = 1, nxv
c           f(1,j,1,l) = f(1,j,2,kr)
c           f(2,j,1,l) = f(2,j,2,kr)
c           f(3,j,1,l) = f(3,j,2,kr)
c  50       continue
c        endif
c        if (kr.eq.nvp) then
c           do 60 j = 1, nxv
c           f(1,j,kyp+3,l) = f(1,j,kyp+2,kr)
c           f(2,j,kyp+3,l) = f(2,j,kyp+2,kr)
c           f(3,j,kyp+3,l) = f(3,j,kyp+2,kr)
c  60       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kl-1,moff+3,lgrp,msid,ier
     1r)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,kyp+1,l),3*nxv,mreal,kr-1,moff+3,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nxv
         f(1,j,1,l) = f(1,j,3,l)
         f(2,j,1,l) = f(2,j,3,l)
         f(3,j,1,l) = f(3,j,3,l)
   20    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(f(1,1,kyp+2,l),3*ngc*nxv,mreal,kr-1,moff+4,lgrp,
     1msid,ierr)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,2,l),3*ngc*nxv,mreal,kl-1,moff+4,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 40 j = 2, nx2
         f(1,j,kyp+3,l) = 2.*f(1,j,kyp+2,l) - f(1,j,kyp+1,l)
         f(2,j,kyp+3,l) = 2.*f(2,j,kyp+2,l) - f(2,j,kyp+1,l)
         f(3,j,kyp+3,l) = 2.*f(3,j,kyp+2,l) - f(3,j,kyp+1,l)
   40    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (krr.le.nvp) then
            call MPI_IRECV(f(1,1,kyp+3,l),3*nxv,mreal,krr-1,moff+6,lgrp,
     1msid,ierr)
         else if (kr.le.nvp) then
            call MPI_IRECV(f(1,1,kyp+3,l),3*nxv,mreal,kr-1,moff+6,lgrp,m
     1sid,ierr)
         endif
         if ((kl.eq.0).and.(kr.le.nvp)) then
            call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kr-1,moff+6,lgrp,nsid,
     1ierr)
         endif
         if (kll.ge.1) then
            call MPI_SEND(f(1,1,2,l),3*nxv,mreal,kll-1,moff+6,lgrp,ierr)
         else if (kl.eq.1) then
            call MPI_SEND(f(1,1,2,l),3*nxv,mreal,kl-1,moff+6,lgrp,ierr)
         endif
         if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
            call MPI_SEND(f(1,1,kyp+2,l),3*nxv,mreal,kl-1,moff+6,lgrp,ie
     1rr)
         endif
         if (kr.le.nvp) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
         if (kl.eq.0) then
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
   70 continue
c fix left edge
      do 100 l = 1, kblok
      kl = l + ks
      if (kl.eq.0) then
         do 90 j = 2, nx2
         f(1,j,1,l) = 2.*f(1,j,2,l) - f(1,j,1,l)
         f(2,j,1,l) = 2.*f(2,j,2,l) - f(2,j,1,l)
         f(3,j,1,l) = 2.*f(3,j,2,l) - f(3,j,1,l)
   90    continue
      endif
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLCGUARD2(f,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mter)
c this subroutine copies data to guard cells in non-uniform partitions
c for vector data the field is replicated so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction.
c f(2,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c scs(2,j,l) = scratch array for field partition l
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, mter
      dimension f(2,nxv,nypmx,nblok), scs(2,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx2, ks, moff, kr, krr, kl, kll, ngc, nps, j, l
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 80 l = 1, nblok
      kr = l + ks + 2
      krr = kr + 1
      kl = l + ks
      kll = kl - 1
      ngc = 0
c special case of only one grid per processor
      if (nyp(l).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nxv
c        f(1,j,1,l) = f(1,j,nyp(kl)+1,kl)
c        f(2,j,1,l) = f(2,j,nyp(kl)+1,kl)
c  10    continue
c     else
c        do 20 j = 1, nxv
c        f(1,j,1,l) = f(1,j,3,l)
c        f(2,j,1,l) = f(2,j,3,l)
c  20    continue
c     endif
c     if (kr.le.nvp) then
c        if (nyp(kr).eq.1) then
c           do 30 j = 1, nxv
c           f(1,j,nyp(l)+2,l) = f(1,j,2,kr)
c           f(2,j,nyp(l)+2,l) = f(2,j,2,kr)
c           f(1,j,nyp(l)+3,l) = f(1,j,2,krr)
c           f(2,j,nyp(l)+3,l) = f(2,j,2,krr)
c  30       continue
c        else
c           do 40 j = 1, nxv
c           f(1,j,nyp(l)+2,l) = f(1,j,2,kr)
c           f(2,j,nyp(l)+2,l) = f(2,j,2,kr)
c           f(1,j,nyp(l)+3,l) = f(1,j,3,kr)
c           f(2,j,nyp(l)+3,l) = f(2,j,3,kr)
c  40       continue
c        endif
c     else
c        do 50 j = 2, nx2
c        f(1,j,nyp(l)+3,l) = 2.*f(1,j,nyp(l)+2,l) - f(1,j,nyp(l)+1,l)
c        f(2,j,nyp(l)+3,l) = 2.*f(2,j,nyp(l)+2,l) - f(2,j,nyp(l)+1,l)
c  50    continue
c     endif
c     if (nyp(l).eq.1) then
c        if ((kl.eq.0).and.(kr.le.nvp)) then
c           do 60 j = 1, nxv
c           f(1,j,1,l) = f(1,j,2,kr)
c           f(2,j,1,l) = f(2,j,2,kr)
c  60       continue
c        endif
c        if (kr.eq.nvp) then
c           do 70 j = 1, nxv
c           f(1,j,nyp(l)+3,l) = f(1,j,nyp(kr)+2,kr)
c           f(2,j,nyp(l)+3,l) = f(2,j,nyp(kr)+2,kr)
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kl-1,moff+3,lgrp,msid,ier
     1r)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,nyp(l)+1,l),2*nxv,mreal,kr-1,moff+3,lgrp,ie
     1rr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nxv
         f(1,j,1,l) = f(1,j,3,l)
         f(2,j,1,l) = f(2,j,3,l)
   20    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(f(1,1,nyp(l)+2,l),4*nxv,mreal,kr-1,moff+4,lgrp,m
     1sid,ierr)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,2,l),2*(2-ngc)*nxv,mreal,kl-1,moff+4,lgrp,i
     1err)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 j = 2, nx2
         f(1,j,nyp(l)+3,l) = 2.*f(1,j,nyp(l)+2,l) - f(1,j,nyp(l)+1,l)
         f(2,j,nyp(l)+3,l) = 2.*f(2,j,nyp(l)+2,l) - f(2,j,nyp(l)+1,l)
   50    continue
      endif
c special case of only one grid per processor
      if (mter.ge.1) go to 80
      if (kr.le.nvp) call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (krr.le.nvp) then
         if (nps.eq.(2*nxv)) then
            call MPI_IRECV(f(1,1,nyp(l)+3,l),2*nxv,mreal,krr-1,moff+6,lg
     1rp,msid,ierr)
         else
            call MPI_IRECV(scs,2*nxv,mreal,krr-1,moff+6,lgrp,msid,ierr)
         endif
      else if (kr.le.nvp) then
         if (nps.eq.(2*nxv)) then
            call MPI_IRECV(f(1,1,nyp(l)+3,l),2*nxv,mreal,kr-1,moff+6,lgr
     1p,msid,ierr)
         else
            call MPI_IRECV(scs,2*nxv,mreal,kr-1,moff+6,lgrp,msid,ierr)
         endif
      endif
      if ((kl.eq.0).and.(kr.le.nvp)) then
         if (ngc.eq.1) then
            call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kr-1,moff+6,lgrp,nsid,
     1ierr)
         else
            call MPI_IRECV(scs,2*nxv,mreal,kr-1,moff+6,lgrp,nsid,ierr)
         endif
      endif
      if (kll.ge.1) then
         call MPI_SEND(f(1,1,2,l),2*nxv,mreal,kll-1,moff+6,lgrp,ierr)
      else if (kl.eq.1) then
         call MPI_SEND(f(1,1,2,l),2*nxv,mreal,kl-1,moff+6,lgrp,ierr)
      endif
      if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
         call MPI_SEND(f(1,1,nyp(l)+2,l),2*nxv,mreal,kl-1,moff+6,lgrp,ie
     1rr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      if (kl.eq.0) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
   80 continue
c fix left edge
      do 110 l = 1, nblok
      kl = l + ks
      if (kl.eq.0) then
         do 100 j = 2, nx2
         f(1,j,1,l) = 2.*f(1,j,2,l) - f(1,j,1,l)
         f(2,j,1,l) = 2.*f(2,j,2,l) - f(2,j,1,l)
  100    continue
      endif
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLDGUARD2(f,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mter)
c this subroutine copies data to guard cells in non-uniform partitions
c for scalar data the field is replicated so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction.
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c scs(j,l) = scratch array for field partition l
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, mter
      dimension f(nxv,nypmx,nblok), scs(nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx2, ks, moff, kr, krr, kl, kll, ngc, nps, j, l
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 80 l = 1, nblok
      kr = l + ks + 2
      krr = kr + 1
      kl = l + ks
      kll = kl - 1
      ngc = 0
c special case of only one grid per processor
      if (nyp(l).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nxv
c        f(j,1,l) = f(j,nyp(kl)+1,kl)
c  10    continue
c     else
c        do 20 j = 1, nxv
c        f(j,1,l) = f(j,3,l)
c  20    continue
c     endif
c     if (kr.le.nvp) then
c        if (nyp(kr).eq.1) then
c           do 30 j = 1, nxv
c           f(j,nyp(l)+2,l) = f(j,2,kr)
c           f(j,nyp(l)+3,l) = f(j,2,krr)
c  30       continue
c        else
c           do 40 j = 1, nxv
c           f(j,nyp(l)+2,l) = f(j,2,kr)
c           f(j,nyp(l)+3,l) = f(j,3,kr)
c  40       continue
c        endif
c     else
c        do 50 j = 2, nx2
c        f(j,nyp(l)+3,l) = 2.*f(j,nyp(l)+2,l) - f(j,nyp(l)+1,l)
c  50    continue
c     endif
c     if (nyp(l).eq.1) then
c        if ((kl.eq.0).and.(kr.le.nvp)) then
c           do 60 j = 1, nxv
c           f(j,1,l) = f(j,2,kr)
c  60       continue
c        endif
c        if (kr.eq.nvp) then
c           do 70 j = 1, nxv
c           f(j,nyp(l)+3,l) = f(j,nyp(kr)+2,kr)
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(f(1,1,l),nxv,mreal,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,nyp(l)+1,l),nxv,mreal,kr-1,moff+3,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nxv
         f(j,1,l) = f(j,3,l)
   20    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(f(1,nyp(l)+2,l),2*nxv,mreal,kr-1,moff+4,lgrp,msi
     1d,ierr)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,2,l),(2-ngc)*nxv,mreal,kl-1,moff+4,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 j = 2, nx2
         f(j,nyp(l)+3,l) = 2.*f(j,nyp(l)+2,l) - f(j,nyp(l)+1,l)
   50    continue
      endif
c special case of only one grid per processor
      if (mter.ge.1) go to 80
      if (kr.le.nvp) call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (krr.le.nvp) then
         if (nps.eq.nxv) then
            call MPI_IRECV(f(1,nyp(l)+3,l),nxv,mreal,krr-1,moff+6,lgrp,m
     1sid,ierr)
         else
            call MPI_IRECV(scs,nxv,mreal,krr-1,moff+6,lgrp,msid,ierr)
         endif
      else if (kr.le.nvp) then
         if (nps.eq.nxv) then
            call MPI_IRECV(f(1,nyp(l)+3,l),nxv,mreal,kr-1,moff+6,lgrp,ms
     1id,ierr)
         else
            call MPI_IRECV(scs,nxv,mreal,kr-1,moff+6,lgrp,msid,ierr)
         endif
      endif
      if ((kl.eq.0).and.(kr.le.nvp)) then
         if (ngc.eq.1) then
            call MPI_IRECV(f(1,1,l),nxv,mreal,kr-1,moff+6,lgrp,nsid,ierr
     1)
         else
            call MPI_IRECV(scs,nxv,mreal,kr-1,moff+6,lgrp,nsid,ierr)
         endif
      endif
      if (kll.ge.1) then
         call MPI_SEND(f(1,2,l),nxv,mreal,kll-1,moff+6,lgrp,ierr)
      else if (kl.eq.1) then
         call MPI_SEND(f(1,2,l),nxv,mreal,kl-1,moff+6,lgrp,ierr)
      endif
      if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
         call MPI_SEND(f(1,nyp(l)+2,l),nxv,mreal,kl-1,moff+6,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      if (kl.eq.0) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
   80 continue
c fix left edge
      do 110 l = 1, nblok
      kl = l + ks
      if (kl.eq.0) then
         do 100 j = 2, nx2
         f(j,1,l) = 2.*f(j,2,l) - f(j,1,l)
  100    continue
      endif
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLBGUARD2(f,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mter)
c this subroutine copies data to guard cells in non-uniform partitions
c for vector data the field is replicated so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction.
c f(3,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c scs(3,j,l) = scratch array for field partition l
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, mter
      dimension f(3,nxv,nypmx,nblok), scs(3,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx2, ks, moff, kr, krr, kl, kll, ngc, nps, j, l
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 80 l = 1, nblok
      kr = l + ks + 2
      krr = kr + 1
      kl = l + ks
      kll = kl - 1
      ngc = 0
c special case of only one grid per processor
      if (nyp(l).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nxv
c        f(1,j,1,l) = f(1,j,nyp(kl)+1,kl)
c        f(2,j,1,l) = f(2,j,nyp(kl)+1,kl)
c        f(3,j,1,l) = f(3,j,nyp(kl)+1,kl)
c  10    continue
c     else
c        do 20 j = 1, nxv
c        f(1,j,1,l) = f(1,j,3,l)
c        f(2,j,1,l) = f(2,j,3,l)
c        f(3,j,1,l) = f(3,j,3,l)
c  20    continue
c     endif
c     if (kr.le.nvp) then
c        if (nyp(kr).eq.1) then
c           do 30 j = 1, nxv
c           f(1,j,nyp(l)+2,l) = f(1,j,2,kr)
c           f(2,j,nyp(l)+2,l) = f(2,j,2,kr)
c           f(3,j,nyp(l)+2,l) = f(3,j,2,kr)
c           f(1,j,nyp(l)+3,l) = f(1,j,2,krr)
c           f(2,j,nyp(l)+3,l) = f(2,j,2,krr)
c           f(3,j,nyp(l)+3,l) = f(3,j,2,krr)
c  30       continue
c        else
c           do 40 j = 1, nxv
c           f(1,j,nyp(l)+2,l) = f(1,j,2,kr)
c           f(2,j,nyp(l)+2,l) = f(2,j,2,kr)
c           f(3,j,nyp(l)+2,l) = f(3,j,2,kr)
c           f(1,j,nyp(l)+3,l) = f(1,j,3,kr)
c           f(2,j,nyp(l)+3,l) = f(2,j,3,kr)
c           f(3,j,nyp(l)+3,l) = f(3,j,3,kr)
c  40       continue
c        endif
c     else
c        do 50 j = 2, nx2
c        f(1,j,nyp(l)+3,l) = 2.*f(1,j,nyp(l)+2,l) - f(1,j,nyp(l)+1,l)
c        f(2,j,nyp(l)+3,l) = 2.*f(2,j,nyp(l)+2,l) - f(2,j,nyp(l)+1,l)
c        f(3,j,nyp(l)+3,l) = 2.*f(3,j,nyp(l)+2,l) - f(3,j,nyp(l)+1,l)
c  50    continue
c     endif
c     if (nyp(l).eq.1) then
c        if ((kl.eq.0).and.(kr.le.nvp)) then
c           do 60 j = 1, nxv
c           f(1,j,1,l) = f(1,j,2,kr)
c           f(2,j,1,l) = f(2,j,2,kr)
c           f(3,j,1,l) = f(3,j,2,kr)
c  60       continue
c        endif
c        if (kr.eq.nvp) then
c           do 70 j = 1, nxv
c           f(1,j,nyp(l)+3,l) = f(1,j,nyp(kr)+2,kr)
c           f(2,j,nyp(l)+3,l) = f(2,j,nyp(kr)+2,kr)
c           f(3,j,nyp(l)+3,l) = f(3,j,nyp(kr)+2,kr)
c  70       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kl-1,moff+3,lgrp,msid,ier
     1r)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,nyp(l)+1,l),3*nxv,mreal,kr-1,moff+3,lgrp,ie
     1rr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nxv
         f(1,j,1,l) = f(1,j,3,l)
         f(2,j,1,l) = f(2,j,3,l)
         f(3,j,1,l) = f(3,j,3,l)
   20    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(f(1,1,nyp(l)+2,l),6*nxv,mreal,kr-1,moff+4,lgrp,m
     1sid,ierr)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,2,l),3*(2-ngc)*nxv,mreal,kl-1,moff+4,lgrp,i
     1err)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 j = 2, nx2
         f(1,j,nyp(l)+3,l) = 2.*f(1,j,nyp(l)+2,l) - f(1,j,nyp(l)+1,l)
         f(2,j,nyp(l)+3,l) = 2.*f(2,j,nyp(l)+2,l) - f(2,j,nyp(l)+1,l)
         f(3,j,nyp(l)+3,l) = 2.*f(3,j,nyp(l)+2,l) - f(3,j,nyp(l)+1,l)
   50    continue
      endif
c special case of only one grid per processor
      if (mter.ge.1) go to 80
      if (kr.le.nvp) call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (krr.le.nvp) then
         if (nps.eq.(3*nxv)) then
            call MPI_IRECV(f(1,1,nyp(l)+3,l),3*nxv,mreal,krr-1,moff+6,lg
     1rp,msid,ierr)
         else
            call MPI_IRECV(scs,3*nxv,mreal,krr-1,moff+6,lgrp,msid,ierr)
         endif
      else if (kr.le.nvp) then
         if (nps.eq.(3*nxv)) then
            call MPI_IRECV(f(1,1,nyp(l)+3,l),3*nxv,mreal,kr-1,moff+6,lgr
     1p,msid,ierr)
         else
            call MPI_IRECV(scs,3*nxv,mreal,kr-1,moff+6,lgrp,msid,ierr)
         endif
      endif
      if ((kl.eq.0).and.(kr.le.nvp)) then
         if (ngc.eq.1) then
            call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kr-1,moff+6,lgrp,nsid,
     1ierr)
         else
            call MPI_IRECV(scs,3*nxv,mreal,kr-1,moff+6,lgrp,nsid,ierr)
         endif
      endif
      if (kll.ge.1) then
         call MPI_SEND(f(1,1,2,l),3*nxv,mreal,kll-1,moff+6,lgrp,ierr)
      else if (kl.eq.1) then
         call MPI_SEND(f(1,1,2,l),3*nxv,mreal,kl-1,moff+6,lgrp,ierr)
      endif
      if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
         call MPI_SEND(f(1,1,nyp(l)+2,l),3*nxv,mreal,kl-1,moff+6,lgrp,ie
     1rr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      if (kl.eq.0) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
   80 continue
c fix left edge
      do 110 l = 1, nblok
      kl = l + ks
      if (kl.eq.0) then
         do 100 j = 2, nx2
         f(1,j,1,l) = 2.*f(1,j,2,l) - f(1,j,1,l)
         f(2,j,1,l) = 2.*f(2,j,2,l) - f(2,j,1,l)
         f(3,j,1,l) = 2.*f(3,j,2,l) - f(3,j,1,l)
  100    continue
      endif
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLCGUARD2L(f,kstrt,nvp,nxv,nypmx,kyp,kblok)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes one extra guard cell.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f
      integer kstrt, nvp, nxv, nypmx, kyp, kblok
      dimension f(nxv,nypmx,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ks, moff, kl, kr, l
c     integer j
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 20 l = 1, kblok
      kr = l + ks + 2
      kl = l + ks
c this loop is used for shared memory computers
c     if (kr.le.nvp) then
c        do 10 j = 1, nxv
c        f(j,kyp+1,l) = f(j,1,kr)
c  10    continue
c     endif
c this segment is used for mpi computers
      if (kr.le.nvp) then
         call MPI_IRECV(f(1,kyp+1,l),nxv,mreal,kr-1,moff+2,lgrp,msid,ier
     1r)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,l),nxv,mreal,kl-1,moff+2,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLCGUARD2L(f,nyp,kstrt,nvp,nxv,nypmx,nblok)
c this subroutine copies data to guard cells in non-uniform partitions
c guard cell on last processor is presumed already set.
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes one extra guard cell.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cell.
c nblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f
      integer nyp
      integer kstrt, nvp, nxv, nypmx, nblok
      dimension f(nxv,nypmx,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ks, moff, kl, kr, l
c     integer j
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c copy to guard cells
      do 20 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
c this loop is used for shared memory computers
c     if (kr.le.nvp) then
c        do 10 j = 1, nxv
c        f(j,nyp(l)+1,l) = f(j,1,kr)
c  10    continue
c     endif
c this segment is used for mpi computers
      if (kr.le.nvp) then
         call MPI_IRECV(f(1,nyp(l)+1,l),nxv,mreal,kr-1,moff+2,lgrp,msid,
     1ierr)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,l),nxv,mreal,kl-1,moff+2,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(3,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes three extra guard cells.
c scr(3,j,ngds,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
      dimension f(3,nxv,nypmx,kblok), scr(3,nxv,ngds,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 170 l = 1, kblok
      kr = l + ks + 2
      krr = kr
      kl = l + ks
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         kll = kll - 1
         ngc = 1
      endif
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 20 j = 1, nx3
c        do 10 m = 1, 3
c        scr(m,j,1,l) = f(m,j,kyp+2,kl)
c        scr(m,j,2,l) = f(m,j,kyp+3,kll)
c  10    continue
c  20    continue
c     else
c        do 40 j = 1, nx3
c        do 30 m = 1, 3
c        scr(m,j,1,l) = 2.*f(m,j,1,l)
c        scr(m,j,2,l) = -f(m,j,1,l)
c  30    continue
c  40    continue
c     endif
c     if (kr.le.nvp) then
c        do 60 j = 1, nx3
c        do 50 m = 1, 3
c        scr(m,j,3,l) = f(m,j,1,kr)
c  50    continue
c  60    continue
c     else
c        do 80 j = 1, nx3
c        do 70 m = 1, 3
c        scr(m,j,3,l) = -f(m,j,kyp+3,l)
c        f(m,j,kyp+2,l) = f(m,j,kyp+2,l) + 2.*f(m,j,kyp+3,l)
c        f(m,j,kyp+3,l) = 0.
c  70    continue
c  80    continue
c     endif
c     if (kyp.eq.1) then
c        if (kl.eq.1) then
c           do 100 j = 1, nx3
c           do 90 m = 1, 3
c           scr(m,j,1,l) = f(m,j,kyp+2,kl)
c           scr(m,j,2,l) = -f(m,j,1,kl)
c  90       continue
c 100       continue
c        else if (kl.eq.0) then
c           do 120 j = 1, nx3
c           do 110 m = 1, 3
c           scr(m,j,2,l) = 0.
c 110       continue
c 120       continue
c        endif
c last point is special with only one grid
c        if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
c           do 140 j = 1, nx3
c           do 130 m = 1, 3
c           f(m,j,kyp+2,l) = f(m,j,kyp+2,l) + f(m,j,kyp+3,kl)
c 130    continue
c 140    continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,3*ngc*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,kyp+2,l),3*ngc*nxv,mreal,kr-1,moff+1,lgrp,i
     1err)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nx3
         do 10 m = 1, 3
         scr(m,j,1,l) = 2.*f(m,j,1,l)
         scr(m,j,2,l) = -f(m,j,1,l)
   10    continue
   20    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(scr(1,1,3,l),3*nxv,mreal,kr-1,moff+2,lgrp,msid,i
     1err)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,1,l),3*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 40 j = 1, nx3
         do 30 m = 1, 3
         scr(m,j,3,l) = -f(m,j,kyp+3,l)
         f(m,j,kyp+2,l) = f(m,j,kyp+2,l) + 2.*f(m,j,kyp+3,l)
         f(m,j,kyp+3,l) = 0.
   30    continue
   40    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (kll.ge.1) then
            call MPI_IRECV(scr(1,1,2,l),3*nxv,mreal,kll-1,moff+5,lgrp,ms
     1id,ierr)
         else if (kl.eq.1) then
            call MPI_IRECV(scr(1,1,2,l),3*nxv,mreal,kl-1,moff+5,lgrp,msi
     1d,ierr)
         endif
         if (krr.le.nvp) then
            call MPI_SEND(f(1,1,kyp+3,l),3*nxv,mreal,krr-1,moff+5,lgrp,i
     1err)
         endif
         if ((kl.eq.0).and.(kr.le.nvp)) then
            call MPI_SEND(f(1,1,1,l),3*nxv,mreal,kr-1,moff+5,lgrp,ierr)
         endif
         if (kl.ge.1) then
            call MPI_WAIT(msid,istatus,ierr)
            if (kl.eq.1) then
               do 60 j = 1, nx3
               do 50 m = 1, 3
               scr(m,j,2,l) = -scr(m,j,2,l)
   50          continue
   60          continue
            endif
         else
            do 80 j = 1, nx3
            do 70 m = 1, 3
            scr(m,j,2,l) = 0.
   70       continue
   80       continue
         endif
c last point is special with only one grid
         if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
            call MPI_IRECV(f(1,1,kyp+3,l),3*nxv,mreal,kl-1,moff+6,lgrp,m
     1sid,ierr)
         endif
         if (kr.eq.nvp) then
            call MPI_SEND(f(1,1,kyp+3,l),3*nxv,mreal,kr-1,moff+6,lgrp,ie
     1rr)
         endif
         if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 140 j = 1, nx3
            do 130 m = 1, 3
            f(m,j,kyp+2,l) = f(m,j,kyp+2,l) + f(m,j,kyp+3,l)
            f(m,j,kyp+3,l) = 0.
  130       continue
  140       continue
         endif
      endif
c add up the guard cells
      do 160 j = 1, nx3
      do 150 m = 1, 3
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,ngc+1,l) = f(m,j,ngc+1,l) + scr(m,j,2,l)
      f(m,j,kyp+1,l) = f(m,j,kyp+1,l) + scr(m,j,3,l)
  150 continue
  160 continue
  170 continue
c zero out the left edge
      do 200 l = 1, kblok
      kl = l + ks
      if (kl.eq.0) then
         do 190 j = 1, nx3
         do 180 m = 1, 3
         f(m,j,1,l) = 0.
  180    continue
  190    continue
      endif
  200 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARD22(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds
     1)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(2,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes three extra guard cells.
c scr(2,j,ngds,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
      dimension f(2,nxv,nypmx,kblok), scr(2,nxv,ngds,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 170 l = 1, kblok
      kr = l + ks + 2
      krr = kr
      kl = l + ks
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         kll = kll - 1
         ngc = 1
      endif
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nx3
c        do 10 m = 1, 2
c        scr(m,j,1,l) = f(m,j,kyp+2,kl)
c        scr(m,j,2,l) = f(m,j,kyp+3,kll)
c  10    continue
c  20    continue
c     else
c        do 40 j = 1, nx3
c        do 30 m = 1, 2
c        scr(m,j,1,l) = 2.*f(m,j,1,l)
c        scr(m,j,2,l) = -f(m,j,1,l)
c  30    continue
c  40    continue
c     endif
c     if (kr.le.nvp) then
c        do 60 j = 1, nx3
c        do 50 m = 1, 2
c        scr(m,j,3,l) = f(m,j,1,kr)
c  50    continue
c  60    continue
c     else
c        do 80 j = 1, nx3
c        do 70 m = 1, 2
c        scr(m,j,3,l) = -f(m,j,kyp+3,l)
c        f(m,j,kyp+2,l) = f(m,j,kyp+2,l) + 2.*f(m,j,kyp+3,l)
c        f(m,j,kyp+3,l) = 0.
c  70    continue
c  80    continue
c     endif
c     if (kyp.eq.1) then
c        if (kl.eq.1) then
c           do 100 j = 1, nx3
c           do 90 m = 1, 2
c           scr(m,j,1,l) = f(m,j,kyp+2,kl)
c           scr(m,j,2,l) = -f(m,j,1,kl)
c  90       continue
c 100       continue
c        else if (kl.eq.0) then
c           do 120 j = 1, nx3
c           do 110 m = 1, 2
c           scr(m,j,2,l) = 0.
c 110       continue
c 120       continue
c        endif
c last point is special with only one grid
c        if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
c           do 140 j = 1, nx3
c           do 130 m = 1, 2
c           f(m,j,kyp+2,l) = f(m,j,kyp+2,l) + f(m,j,kyp+3,kl)
c 130    continue
c 140    continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,2*ngc*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,kyp+2,l),2*ngc*nxv,mreal,kr-1,moff+1,lgrp,i
     1err)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nx3
         do 10 m = 1, 2
         scr(m,j,1,l) = 2.*f(m,j,1,l)
         scr(m,j,2,l) = -f(m,j,1,l)
   10    continue
   20    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(scr(1,1,3,l),2*nxv,mreal,kr-1,moff+2,lgrp,msid,i
     1err)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,1,l),2*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 40 j = 1, nx3
         do 30 m = 1, 2
         scr(m,j,3,l) = -f(m,j,kyp+3,l)
         f(m,j,kyp+2,l) = f(m,j,kyp+2,l) + 2.*f(m,j,kyp+3,l)
         f(m,j,kyp+3,l) = 0.
   30    continue
   40    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (kll.ge.1) then
            call MPI_IRECV(scr(1,1,2,l),2*nxv,mreal,kll-1,moff+5,lgrp,ms
     1id,ierr)
         else if (kl.eq.1) then
            call MPI_IRECV(scr(1,1,2,l),2*nxv,mreal,kl-1,moff+5,lgrp,msi
     1d,ierr)
         endif
         if (krr.le.nvp) then
            call MPI_SEND(f(1,1,kyp+3,l),2*nxv,mreal,krr-1,moff+5,lgrp,i
     1err)
         endif
         if ((kl.eq.0).and.(kr.le.nvp)) then
            call MPI_SEND(f(1,1,1,l),2*nxv,mreal,kr-1,moff+5,lgrp,ierr)
         endif
         if (kl.ge.1) then
            call MPI_WAIT(msid,istatus,ierr)
            if (kl.eq.1) then
               do 60 j = 1, nx3
               do 50 m = 1, 2
               scr(m,j,2,l) = -scr(m,j,2,l)
   50          continue
   60          continue
            endif
         else
            do 80 j = 1, nx3
            do 70 m = 1, 2
            scr(m,j,2,l) = 0.
   70       continue
   80       continue
         endif
c last point is special with only one grid
         if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
            call MPI_IRECV(f(1,1,kyp+3,l),2*nxv,mreal,kl-1,moff+6,lgrp,m
     1sid,ierr)
         endif
         if (kr.eq.nvp) then
            call MPI_SEND(f(1,1,kyp+3,l),2*nxv,mreal,kr-1,moff+6,lgrp,ie
     1rr)
         endif
         if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 140 j = 1, nx3
            do 130 m = 1, 2
            f(m,j,kyp+2,l) = f(m,j,kyp+2,l) + f(m,j,kyp+3,l)
            f(m,j,kyp+3,l) = 0.
  130       continue
  140       continue
         endif
      endif
c add up the guard cells
      do 160 j = 1, nx3
      do 150 m = 1, 2
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,ngc+1,l) = f(m,j,ngc+1,l) + scr(m,j,2,l)
      f(m,j,kyp+1,l) = f(m,j,kyp+1,l) + scr(m,j,3,l)
  150 continue
  160 continue
  170 continue
c zero out the left edge
      do 200 l = 1, kblok
      kl = l + ks
      if (kl.eq.0) then
         do 190 j = 1, nx3
         do 180 m = 1, 2
         f(m,j,1,l) = 0.
  180    continue
  190    continue
      endif
  200 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLAGUARD2(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes three extra guard cells.
c scr(j,ngds,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok, ngds
      dimension f(nxv,nypmx,kblok), scr(nxv,ngds,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 90 l = 1, kblok
      kr = l + ks + 2
      krr = kr
      kl = l + ks
      kll = kl
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         kll = kll - 1
         ngc = 1
      endif
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nx3
c        scr(j,1,l) = f(j,kyp+2,kl)
c        scr(j,2,l) = f(j,kyp+3,kll)
c  10    continue
c     else
c        do 20 j = 1, nx3
c        scr(j,1,l) = 2.*f(j,1,l)
c        scr(j,2,l) = -f(j,1,l)
c  20    continue
c     endif
c     if (kr.le.nvp) then
c        do 30 j = 1, nx3
c        scr(j,3,l) = f(j,1,kr)
c  30    continue
c     else
c        do 40 j = 1, nx3
c        scr(j,3,l) = -f(j,kyp+3,l)
c        f(j,kyp+2,l) = f(j,kyp+2,l) + 2.*f(j,kyp+3,l)
c        f(j,kyp+3,l) = 0.
c  40    continue
c     endif
c     if (kyp.eq.1) then
c        if (kl.eq.1) then
c           do 50 j = 1, nx3
c           scr(j,1,l) = f(j,kyp+2,kl)
c           scr(j,2,l) = -f(j,1,kl)
c  50       continue
c        else if (kl.eq.0) then
c           do 60 j = 1, nx3
c           scr(j,2,l) = 0.
c  60       continue
c        endif
c last point is special with only one grid
c        if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
c           do 70 j = 1, nx3
c           f(j,kyp+2,l) = f(j,kyp+2,l) + f(j,kyp+3,kl)
c  70    continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,ngc*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,kyp+2,l),ngc*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 10 j = 1, nx3
         scr(j,1,l) = 2.*f(j,1,l)
         scr(j,2,l) = -f(j,1,l)
   10    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(scr(1,3,l),nxv,mreal,kr-1,moff+2,lgrp,msid,ierr)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,l),nxv,mreal,kl-1,moff+2,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nx3
         scr(j,3,l) = -f(j,kyp+3,l)
         f(j,kyp+2,l) = f(j,kyp+2,l) + 2.*f(j,kyp+3,l)
         f(j,kyp+3,l) = 0.
   20    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (kll.ge.1) then
            call MPI_IRECV(scr(1,2,l),nxv,mreal,kll-1,moff+5,lgrp,msid,i
     1err)
         else if (kl.eq.1) then
            call MPI_IRECV(scr(1,2,l),nxv,mreal,kl-1,moff+5,lgrp,msid,ie
     1rr)
         endif
         if (krr.le.nvp) then
            call MPI_SEND(f(1,kyp+3,l),nxv,mreal,krr-1,moff+5,lgrp,ierr)
         endif
         if ((kl.eq.0).and.(kr.le.nvp)) then
            call MPI_SEND(f(1,1,l),nxv,mreal,kr-1,moff+5,lgrp,ierr)
         endif
         if (kl.ge.1) then
            call MPI_WAIT(msid,istatus,ierr)
            if (kl.eq.1) then
               do 30 j = 1, nx3
               scr(j,2,l) = -scr(j,2,l)
   30          continue
            endif
         else
            do 40 j = 1, nx3
            scr(j,2,l) = 0.
   40       continue
         endif
c last point is special with only one grid
         if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
            call MPI_IRECV(f(1,kyp+3,l),nxv,mreal,kl-1,moff+6,lgrp,msid,
     1ierr)
         endif
         if (kr.eq.nvp) then
            call MPI_SEND(f(1,kyp+3,l),nxv,mreal,kr-1,moff+6,lgrp,ierr)
         endif
         if ((kl.eq.(nvp-1)).and.(kl.ge.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 70 j = 1, nx3
            f(j,kyp+2,l) = f(j,kyp+2,l) + f(j,kyp+3,l)
            f(j,kyp+3,l) = 0.
   70       continue
         endif
      endif
c add up the guard cells
      do 80 j = 1, nx3
      f(j,2,l) = f(j,2,l) + scr(j,1,l)
      f(j,ngc+1,l) = f(j,ngc+1,l) + scr(j,2,l)
      f(j,kyp+1,l) = f(j,kyp+1,l) + scr(j,3,l)
   80 continue
   90 continue
c zero out the left edge
      do 110 l = 1, kblok
      kl = l + ks
      if (kl.eq.0) then
         do 100 j = 1, nx3
         f(j,1,l) = 0.
  100    continue
      endif
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLACGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,
     1ngds,mter)
c this subroutine adds data from guard cells in non-uniform partitions
c for vector data, the field is added up so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction.
c f(3,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c scr(3,j,ngds,l) = scratch array for particle partition l
c scs(3,j,l) = scratch array for field partition l
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c ngds = number of guard cells
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scr, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
      dimension f(3,nxv,nypmx,nblok), scr(3,nxv,ngds,nblok)
      dimension scs(3,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 220 l = 1, nblok
      kr = l + ks + 2
      krr = kr + 1
      kl = l + ks
      kll = kl - 1
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        if (nyp(kl).eq.1) then
c           do 20 j = 1, nx3
c           do 10 m = 1, 3
c           scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl) + f(m,j,nyp(kll)+2,kll)
c           scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c  10       continue
c  20       continue
c        else
c           do 40 j = 1, nx3
c           do 30 m = 1, 3
c           scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl)
c           scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c  30       continue
c  40       continue
c        endif
c     else
c        do 60 j = 1, nx3
c        do 50 m = 1, 3
c        scr(m,j,1,l) = 2.*f(m,j,1,l)
c        scr(m,j,2,l) = -f(m,j,1,l)
c  50    continue
c  60    continue
c     endif
c     if (kr.le.nvp) then
c        do 80 j = 1, nx3
c        do 70 m = 1, 3
c        scr(m,j,3,l) = f(m,j,1,kr)
c  70    continue
c  80    continue
c     else
c        do 100 j = 1, nx3
c        do 90 m = 1, 3
c        scr(m,j,3,l) = -f(m,j,nyp(l)+3,l)
c        f(m,j,nyp(l)+2,l) = f(m,j,nyp(l)+2,l) + 2.*f(m,j,nyp(l)+3,l)
c  90    continue
c 100    continue
c     endif
c     if (kl.ge.1) then
c        if ((nyp(kl).eq.1).and.(kl.eq.1)) then
c           do 120 j = 1, nx3
c           do 110 m = 1, 3
c           scr(m,j,1,l) = scr(m,j,1,l) - f(m,j,1,kl)
c 110       continue
c 120       continue
c        endif
c     else if (nyp(l).eq.1) then 
c        do 160 j = 1, nx3
c        do 150 m = 1, 3
c        scr(m,j,2,l) = 0.
c 150    continue
c 160    continue
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,6*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,nyp(l)+2,l),6*nxv,mreal,kr-1,moff+1,lgrp,ie
     1rr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 60 j = 1, nx3
         do 50 m = 1, 3
         scr(m,j,1,l) = 2.*f(m,j,1,l)
         scr(m,j,2,l) = -f(m,j,1,l)
   50    continue
   60    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(scr(1,1,3,l),3*nxv,mreal,kr-1,moff+2,lgrp,msid,i
     1err)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,1,l),3*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 j = 1, nx3
         do 90 m = 1, 3
         scr(m,j,3,l) = -f(m,j,nyp(l)+3,l)
         f(m,j,nyp(l)+2,l) = f(m,j,nyp(l)+2,l) + 2.*f(m,j,nyp(l)+3,l)
   90    continue
  100    continue
      endif
c special case of only one grid per processor
      if (mter.ge.1) go to 170
      if (kl.ge.1) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (kll.ge.1) then
         call MPI_IRECV(scs,3*nxv,mreal,kll-1,moff+5,lgrp,nsid,ierr)
      else if (kl.eq.1) then
         call MPI_IRECV(scs,3*nxv,mreal,kl-1,moff+5,lgrp,nsid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(nyp(l),1,mint,kr-1,moff+3,lgrp,ierr)
      endif
      if (krr.le.nvp) then
         call MPI_SEND(f(1,1,nyp(l)+3,l),3*nxv,mreal,krr-1,moff+5,lgrp,i
     1err)
      endif
      if ((kl.eq.0).and.(kr.le.nvp)) then
         call MPI_SEND(f(1,1,1,l),3*nxv,mreal,kr-1,moff+5,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         if (ngc.eq.1) then
            if (kl.eq.1) then
               do 120 j = 1, nx3
               do 110 m = 1, 3
               scr(m,j,1,l) = scr(m,j,1,l) - scs(m,j,l)
  110          continue
  120          continue
            else
               do 140 j = 1, nx3
               do 130 m = 1, 3
               scr(m,j,1,l) = scr(m,j,1,l) + scs(m,j,l)
  130          continue
  140          continue
            endif
         endif
      else if (nyp(l).eq.1) then
         do 160 j = 1, nx3
         do 150 m = 1, 3
         scr(m,j,2,l) = 0.
  150    continue
  160    continue
      endif
c add up the guard cells
  170 do 190 j = 1, nx3
      do 180 m = 1, 3
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,3,l) = f(m,j,3,l) + scr(m,j,2,l)
      f(m,j,nyp(l)+1,l) = f(m,j,nyp(l)+1,l) + scr(m,j,3,l)
      f(m,j,1,l) = 0.
      f(m,j,nyp(l)+3,l) = 0.
  180 continue
  190 continue
      if (kr.le.nvp) then
         do 210 j = 1, nx3
         do 200 m = 1, 3
         f(m,j,nyp(l)+2,l) = 0.
  200    continue
  210    continue
      endif
  220 continue
c zero out the left edge
      do 250 l = 1, nblok
      kl = l + ks
      if (kl.eq.0) then
         do 240 j = 1, nx3
         do 230 m = 1, 3
         f(m,j,1,l) = 0.
  230    continue
  240    continue
      endif
  250 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLACGUARD22(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok
     1,ngds,mter)
c this subroutine adds data from guard cells in non-uniform partitions
c for vector data, the field is added up so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction.
c f(2,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c scr(2,j,ngds,l) = scratch array for particle partition l
c scs(2,j,l) = scratch array for field partition l
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c ngds = number of guard cells
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scr, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
      dimension f(2,nxv,nypmx,nblok), scr(2,nxv,ngds,nblok)
      dimension scs(2,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l, m
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 220 l = 1, nblok
      kr = l + ks + 2
      krr = kr + 1
      kl = l + ks
      kll = kl - 1
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        if (nyp(kl).eq.1) then
c           do 20 j = 1, nx3
c           do 10 m = 1, 2
c           scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl) + f(m,j,nyp(kll)+2,kll)
c           scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c  10       continue
c  20       continue
c        else
c           do 40 j = 1, nx3
c           do 30 m = 1, 2
c           scr(m,j,1,l) = f(m,j,nyp(kl)+2,kl)
c           scr(m,j,2,l) = f(m,j,nyp(kl)+3,kl)
c  30       continue
c  40       continue
c        endif
c     else
c        do 60 j = 1, nx3
c        do 50 m = 1, 2
c        scr(m,j,1,l) = 2.*f(m,j,1,l)
c        scr(m,j,2,l) = -f(m,j,1,l)
c  50    continue
c  60    continue
c     endif
c     if (kr.le.nvp) then
c        do 80 j = 1, nx3
c        do 70 m = 1, 2
c        scr(m,j,3,l) = f(m,j,1,kr)
c  70    continue
c  80    continue
c     else
c        do 100 j = 1, nx3
c        do 90 m = 1, 2
c        scr(m,j,3,l) = -f(m,j,nyp(l)+3,l)
c        f(m,j,nyp(l)+2,l) = f(m,j,nyp(l)+2,l) + 2.*f(m,j,nyp(l)+3,l)
c  90    continue
c 100    continue
c     endif
c     if (kl.ge.1) then
c        if ((nyp(kl).eq.1).and.(kl.eq.1)) then
c           do 120 j = 1, nx3
c           do 110 m = 1, 2
c           scr(m,j,1,l) = scr(m,j,1,l) - f(m,j,1,kl)
c 110       continue
c 120       continue
c        endif
c     else if (nyp(l).eq.1) then 
c        do 160 j = 1, nx3
c        do 150 m = 1, 2
c        scr(m,j,2,l) = 0.
c 150    continue
c 160    continue
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,4*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,nyp(l)+2,l),4*nxv,mreal,kr-1,moff+1,lgrp,ie
     1rr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 60 j = 1, nx3
         do 50 m = 1, 2
         scr(m,j,1,l) = 2.*f(m,j,1,l)
         scr(m,j,2,l) = -f(m,j,1,l)
   50    continue
   60    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(scr(1,1,3,l),2*nxv,mreal,kr-1,moff+2,lgrp,msid,i
     1err)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,1,l),2*nxv,mreal,kl-1,moff+2,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 j = 1, nx3
         do 90 m = 1, 2
         scr(m,j,3,l) = -f(m,j,nyp(l)+3,l)
         f(m,j,nyp(l)+2,l) = f(m,j,nyp(l)+2,l) + 2.*f(m,j,nyp(l)+3,l)
   90    continue
  100    continue
      endif
c special case of only one grid per processor
      if (mter.ge.1) go to 170
      if (kl.ge.1) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (kll.ge.1) then
         call MPI_IRECV(scs,2*nxv,mreal,kll-1,moff+5,lgrp,nsid,ierr)
      else if (kl.eq.1) then
         call MPI_IRECV(scs,2*nxv,mreal,kl-1,moff+5,lgrp,nsid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(nyp(l),1,mint,kr-1,moff+3,lgrp,ierr)
      endif
      if (krr.le.nvp) then
         call MPI_SEND(f(1,1,nyp(l)+3,l),2*nxv,mreal,krr-1,moff+5,lgrp,i
     1err)
      endif
      if ((kl.eq.0).and.(kr.le.nvp)) then
         call MPI_SEND(f(1,1,1,l),2*nxv,mreal,kr-1,moff+5,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         if (ngc.eq.1) then
            if (kl.eq.1) then
               do 120 j = 1, nx3
               do 110 m = 1, 2
               scr(m,j,1,l) = scr(m,j,1,l) - scs(m,j,l)
  110          continue
  120          continue
            else
               do 140 j = 1, nx3
               do 130 m = 1, 2
               scr(m,j,1,l) = scr(m,j,1,l) + scs(m,j,l)
  130          continue
  140          continue
            endif
         endif
      else if (nyp(l).eq.1) then
         do 160 j = 1, nx3
         do 150 m = 1, 2
         scr(m,j,2,l) = 0.
  150    continue
  160    continue
      endif
c add up the guard cells
  170 do 190 j = 1, nx3
      do 180 m = 1, 2
      f(m,j,2,l) = f(m,j,2,l) + scr(m,j,1,l)
      f(m,j,3,l) = f(m,j,3,l) + scr(m,j,2,l)
      f(m,j,nyp(l)+1,l) = f(m,j,nyp(l)+1,l) + scr(m,j,3,l)
      f(m,j,1,l) = 0.
      f(m,j,nyp(l)+3,l) = 0.
  180 continue
  190 continue
      if (kr.le.nvp) then
         do 210 j = 1, nx3
         do 200 m = 1, 3
         f(m,j,nyp(l)+2,l) = 0.
  200    continue
  210    continue
      endif
  220 continue
c zero out the left edge
      do 250 l = 1, nblok
      kl = l + ks
      if (kl.eq.0) then
         do 240 j = 1, nx3
         do 230 m = 1, 2
         f(m,j,1,l) = 0.
  230    continue
  240    continue
      endif
  250 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLAGUARD2(f,scr,scs,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,n
     1gds,mter)
c this subroutine adds data from guard cells in non-uniform partitions
c for scalar data, the field is added up so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction.
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c scr(j,ngds,l) = scratch array for particle partition l
c scs(j,l) = scratch array for field partition l
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c ngds = number of guard cells
c mter = (0,1) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f, scr, scs
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, ngds, mter
      dimension f(nxv,nypmx,nblok), scr(nxv,ngds,nblok)
      dimension scs(nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer nx3, ks, moff, kr, krr, kl, kll, ngc, j, l
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 120 l = 1, nblok
      kr = l + ks + 2
      krr = kr + 1
      kl = l + ks
      kll = kl - 1
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        if (nyp(kl).eq.1) then
c           do 10 j = 1, nx3
c           scr(j,1,l) = f(j,nyp(kl)+2,kl) + f(j,nyp(kll)+2,kll)
c           scr(j,2,l) = f(j,nyp(kl)+3,kl)
c  10       continue
c        else
c           do 20 j = 1, nx3
c           scr(j,1,l) = f(j,nyp(kl)+2,kl)
c           scr(j,2,l) = f(j,nyp(kl)+3,kl)
c  20       continue
c     else
c        do 30 j = 1, nx3
c        scr(j,1,l) = 2.*f(j,1,l)
c        scr(j,2,l) = -f(j,1,l)
c  30    continue
c     endif
c     if (kr.le.nvp) then
c        do 40 j = 1, nx3
c        scr(j,3,l) = f(j,1,kr)
c  40    continue
c     else
c        do 50 j = 1, nx3
c        scr(j,3,l) = -f(j,nyp(l)+3,l)
c        f(j,nyp(l)+2,l) = f(j,nyp(l)+2,l) + 2.*f(j,nyp(l)+3,l)
c  50    continue
c     endif
c     if (kl.ge.1) then
c        if ((nyp(kl).eq.1).and.(kl.eq.1)) then
c           do 60 j = 1, nx3
c           scr(j,1,l) = scr(j,1,l) - f(j,1,kl)
c  60       continue
c        endif
c     else if (nyp(l).eq.1) then
c        do 80 j = 1, nx3
c        scr(j,2,l) = 0.
c  80    continue
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,2*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,nyp(l)+2,l),2*nxv,mreal,kr-1,moff+1,lgrp,ierr
     1)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 30 j = 1, nx3
         scr(j,1,l) = 2.*f(j,1,l)
         scr(j,2,l) = -f(j,1,l)
   30    continue
      endif
      if (kr.le.nvp) then
         call MPI_IRECV(scr(1,3,l),nxv,mreal,kr-1,moff+2,lgrp,msid,ierr)
      endif
      if (kl.ge.1) then
         call MPI_SEND(f(1,1,l),nxv,mreal,kl-1,moff+2,lgrp,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 j = 1, nx3
         scr(j,3,l) = -f(j,nyp(l)+3,l)
         f(j,nyp(l)+2,l) = f(j,nyp(l)+2,l) + 2.*f(j,nyp(l)+3,l)
   50    continue
      endif
c special case of only one grid per processor
      if (mter.ge.1) go to 90
      if (kl.ge.1) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (kll.ge.1) then
         call MPI_IRECV(scs,nxv,mreal,kll-1,moff+5,lgrp,nsid,ierr)
      else if (kl.eq.1) then
         call MPI_IRECV(scs,nxv,mreal,kl-1,moff+5,lgrp,nsid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(nyp(l),1,mint,kr-1,moff+3,lgrp,ierr)
      endif
      if (krr.le.nvp) then
         call MPI_SEND(f(1,nyp(l)+3,l),nxv,mreal,krr-1,moff+5,lgrp,ierr)
      endif
      if ((kl.eq.0).and.(kr.le.nvp)) then
         call MPI_SEND(f(1,1,l),nxv,mreal,kr-1,moff+5,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         if (ngc.eq.1) then
            if (kl.eq.1) then
               do 60 j = 1, nx3
               scr(j,1,l) = scr(j,1,l) - scs(j,l)
   60          continue
            else
               do 70 j = 1, nx3
               scr(j,1,l) = scr(j,1,l) + scs(j,l)
   70          continue
            endif
         endif
      else if (nyp(l).eq.1) then
         do 80 j = 1, nx3
         scr(j,2,l) = 0.
   80    continue
      endif
c add up the guard cells
   90 do 100 j = 1, nx3
      f(j,2,l) = f(j,2,l) + scr(j,1,l)
      f(j,3,l) = f(j,3,l) + scr(j,2,l)
      f(j,nyp(l)+1,l) = f(j,nyp(l)+1,l) + scr(j,3,l)
      f(j,1,l) = 0.
      f(j,nyp(l)+3,l) = 0.
  100 continue
      if (kr.le.nvp) then
         do 110 j = 1, nx3
         f(j,nyp(l)+2,l) = 0.
  110    continue
      endif
  120 continue
c zero out the left edge
      do 140 l = 1, nblok
      kl = l + ks
      if (kl.eq.0) then
         do 130 j = 1, nx3
         f(j,1,l) = 0.
  130    continue
      endif
  140 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARDS2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine corrects current density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(3,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes three extra guard cells.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(3,nxv,nypmx,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ks, moff, kr, krr, kl, kll, j, l, m
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 210 l = 1, kblok
      kr = l + ks + 2
      krr = kr
      kl = l + ks
      kll = kl
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         kll = kll - 1
      endif
c fix edges if all points are on the same processor
      if (kl.eq.0) then
         if (kyp.gt.2) then
            do 20 j = 2, nx
            do 10 m = 1, 3
            f(m,j+1,3,l) = f(m,j+1,3,l) + 2.*f(m,j+1,2,l)
            f(m,j+1,4,l) = f(m,j+1,4,l) - f(m,j+1,2,l)
            f(m,j+1,2,l) = 0.
   10       continue
   20       continue
         else if (kyp.eq.2) then
            do 40 j = 2, nx
            do 30 m = 1, 3
            f(m,j+1,3,l) = f(m,j+1,3,l) + 2.*f(m,j+1,2,l)
   30       continue
   40       continue
         endif
      endif
      if (kr.eq.(nvp+1)) then
         if (kyp.gt.1) then
            do 60 j = 2, nx
            do 50 m = 1, 3
            f(m,j+1,kyp,l) = f(m,j+1,kyp,l) - f(m,j+1,kyp+2,l)
            f(m,j+1,kyp+1,l) = f(m,j+1,kyp+1,l) + 2.*f(m,j+1,kyp+2,l)
            f(m,j+1,kyp+2,l) = 0.
   50       continue
   60       continue
         else if (kyp.eq.1) then
            do 80 j = 2, nx
            do 70 m = 1, 3
            f(m,j+1,kyp+1,l) = f(m,j+1,kyp+1,l) + 2.*f(m,j+1,kyp+2,l)
   70       continue
   80       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (kyp.eq.2) then
c        if (kl.eq.1) then
c           do 120 j = 2, nx
c           do 110 m = 1, 3
c           f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,2,kl)
c 110       continue
c 120       continue
c        endif
c     else if (kyp.eq.1) then
c        if (kl.eq.1) then
c           do 140 j = 2, nx
c           do 130 m = 1, 3
c           f(m,j+1,2,l) = f(m,j+1,2,l) + 2.*f(m,j+1,2,kl)
c 130       continue
c 140       continue
c        endif
c        if (kll.eq.1) then
c           do 160 j = 2, nx
c           do 150 m = 1, 3
c           f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,2,kll)
c 150       continue
c 160       continue
c        endif
c        if (kr.eq.nvp) then
c           do 180 j = 2, nx
c           do 170 m = 1, 3
c           f(m,j+1,kyp+1,l) = f(m,j+1,kyp+1,l) - f(m,j+1,kyp+2,kr)
c 170       continue
c 180       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kyp.eq.2) then
         if (kl.eq.1) then
            call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kl-1,moff+1,lgrp,msid,
     1ierr)
         endif
         if (kl.eq.0) then
            call MPI_SEND(f(1,1,2,l),3*nxv,mreal,kr-1,moff+1,lgrp,ierr)
            do 100 j = 2, nx
            do 90 m = 1, 3
            f(m,j+1,2,l) = 0.
   90       continue
  100       continue
         endif
         if (kl.eq.1) then
            call MPI_WAIT(msid,istatus,ierr)
            do 120 j = 2, nx
            do 110 m = 1, 3
            f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  110       continue
  120       continue
         endif
      else if (kyp.eq.1) then
         if (kl.eq.1) then
            call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kl-1,moff+1,lgrp,msid,
     1ierr)
         endif
         if (kll.eq.1) then
            call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kll-1,moff+1,lgrp,nsid
     1,ierr)
         endif
         if (kl.eq.0) then
            call MPI_SEND(f(1,1,2,l),3*nxv,mreal,kr-1,moff+1,lgrp,ierr)
            call MPI_SEND(f(1,1,2,l),3*nxv,mreal,krr-1,moff+1,lgrp,ierr)
            do 140 j = 2, nx
            do 130 m = 1, 3
            f(m,j+1,2,l) = 0.
  130       continue
  140       continue
         endif
         if (kl.eq.1) then
            call MPI_WAIT(msid,istatus,ierr)
            do 160 j = 2, nx
            do 150 m = 1, 3
            f(m,j+1,2,l) = f(m,j+1,2,l) + 2.*f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  150       continue
  160       continue
         endif
         if (kll.eq.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 180 j = 2, nx
            do 170 m = 1, 3
            f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  170       continue
  180       continue
         endif
         if (kr.eq.nvp) then
            call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kr-1,moff+2,lgrp,msid,
     1ierr)
         endif
         if (kr.eq.(nvp+1)) then
            call MPI_SEND(f(1,1,kyp+2,l),3*nxv,mreal,kl-1,moff+2,lgrp,ie
     1rr)
         endif
         if (kr.eq.nvp) then
            call MPI_WAIT(msid,istatus,ierr)
            do 200 j = 2, nx
            do 190 m = 1, 3
            f(m,j+1,kyp+1,l) = f(m,j+1,kyp+1,l) - f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  190       continue
  200       continue
         endif
      endif
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARDS22(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine corrects current density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(2,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes three extra guard cells.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(2,nxv,nypmx,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ks, moff, kr, krr, kl, kll, j, l, m
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 210 l = 1, kblok
      kr = l + ks + 2
      krr = kr
      kl = l + ks
      kll = kl
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         kll = kll - 1
      endif
c fix edges if all points are on the same processor
      if (kl.eq.0) then
         if (kyp.gt.2) then
            do 20 j = 2, nx
            do 10 m = 1, 2
            f(m,j+1,3,l) = f(m,j+1,3,l) + 2.*f(m,j+1,2,l)
            f(m,j+1,4,l) = f(m,j+1,4,l) - f(m,j+1,2,l)
            f(m,j+1,2,l) = 0.
   10       continue
   20       continue
         else if (kyp.eq.2) then
            do 40 j = 2, nx
            do 30 m = 1, 2
            f(m,j+1,3,l) = f(m,j+1,3,l) + 2.*f(m,j+1,2,l)
   30       continue
   40       continue
         endif
      endif
      if (kr.eq.(nvp+1)) then
         if (kyp.gt.1) then
            do 60 j = 2, nx
            do 50 m = 1, 2
            f(m,j+1,kyp,l) = f(m,j+1,kyp,l) - f(m,j+1,kyp+2,l)
            f(m,j+1,kyp+1,l) = f(m,j+1,kyp+1,l) + 2.*f(m,j+1,kyp+2,l)
            f(m,j+1,kyp+2,l) = 0.
   50       continue
   60       continue
         else if (kyp.eq.1) then
            do 80 j = 2, nx
            do 70 m = 1, 2
            f(m,j+1,kyp+1,l) = f(m,j+1,kyp+1,l) + 2.*f(m,j+1,kyp+2,l)
   70       continue
   80       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (kyp.eq.2) then
c        if (kl.eq.1) then
c           do 120 j = 2, nx
c           do 110 m = 1, 2
c           f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,2,kl)
c 110       continue
c 120       continue
c        endif
c     else if (kyp.eq.1) then
c        if (kl.eq.1) then
c           do 140 j = 2, nx
c           do 130 m = 1, 2
c           f(m,j+1,2,l) = f(m,j+1,2,l) + 2.*f(m,j+1,2,kl)
c 130       continue
c 140       continue
c        endif
c        if (kll.eq.1) then
c           do 160 j = 2, nx
c           do 150 m = 1, 2
c           f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,2,kll)
c 150       continue
c 160       continue
c        endif
c        if (kr.eq.nvp) then
c           do 180 j = 2, nx
c           do 170 m = 1, 2
c           f(m,j+1,kyp+1,l) = f(m,j+1,kyp+1,l) - f(m,j+1,kyp+2,kr)
c 170       continue
c 180       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kyp.eq.2) then
         if (kl.eq.1) then
            call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kl-1,moff+1,lgrp,msid,
     1ierr)
         endif
         if (kl.eq.0) then
            call MPI_SEND(f(1,1,2,l),2*nxv,mreal,kr-1,moff+1,lgrp,ierr)
            do 100 j = 2, nx
            do 90 m = 1, 2
            f(m,j+1,2,l) = 0.
   90       continue
  100       continue
         endif
         if (kl.eq.1) then
            call MPI_WAIT(msid,istatus,ierr)
            do 120 j = 2, nx
            do 110 m = 1, 2
            f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  110       continue
  120       continue
         endif
      else if (kyp.eq.1) then
         if (kl.eq.1) then
            call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kl-1,moff+1,lgrp,msid,
     1ierr)
         endif
         if (kll.eq.1) then
            call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kll-1,moff+1,lgrp,nsid
     1,ierr)
         endif
         if (kl.eq.0) then
            call MPI_SEND(f(1,1,2,l),2*nxv,mreal,kr-1,moff+1,lgrp,ierr)
            call MPI_SEND(f(1,1,2,l),2*nxv,mreal,krr-1,moff+1,lgrp,ierr)
            do 140 j = 2, nx
            do 130 m = 1, 2
            f(m,j+1,2,l) = 0.
  130       continue
  140       continue
         endif
         if (kl.eq.1) then
            call MPI_WAIT(msid,istatus,ierr)
            do 160 j = 2, nx
            do 150 m = 1, 2
            f(m,j+1,2,l) = f(m,j+1,2,l) + 2.*f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  150       continue
  160       continue
         endif
         if (kll.eq.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 180 j = 2, nx
            do 170 m = 1, 2
            f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  170       continue
  180       continue
         endif
         if (kr.eq.nvp) then
            call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kr-1,moff+2,lgrp,msid,
     1ierr)
         endif
         if (kr.eq.(nvp+1)) then
            call MPI_SEND(f(1,1,kyp+2,l),2*nxv,mreal,kl-1,moff+2,lgrp,ie
     1rr)
         endif
         if (kr.eq.nvp) then
            call MPI_WAIT(msid,istatus,ierr)
            do 200 j = 2, nx
            do 190 m = 1, 2
            f(m,j+1,kyp+1,l) = f(m,j+1,kyp+1,l) - f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  190       continue
  200       continue
         endif
      endif
  210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLAGUARDS2(f,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine corrects the charge density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes three extra guard cells.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(nxv,nypmx,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ks, moff, kr, krr, kl, kll, j, l
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 110 l = 1, kblok
      kr = l + ks + 2
      krr = kr
      kl = l + ks
      kll = kl
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = krr + 1
         kll = kll - 1
      endif
c fix edges if all points are on the same processor
      if (kl.eq.0) then
         if (kyp.gt.2) then
            do 10 j = 2, nx
            f(j+1,3,l) = f(j+1,3,l) + 2.*f(j+1,2,l)
            f(j+1,4,l) = f(j+1,4,l) - f(j+1,2,l)
            f(j+1,2,l) = 0.
   10       continue
         else if (kyp.eq.2) then
            do 20 j = 2, nx
            f(j+1,3,l) = f(j+1,3,l) + 2.*f(j+1,2,l)
   20       continue
         endif
      endif
      if (kr.eq.(nvp+1)) then
         if (kyp.gt.1) then
            do 30 j = 2, nx
            f(j+1,kyp,l) = f(j+1,kyp,l) - f(j+1,kyp+2,l)
            f(j+1,kyp+1,l) = f(j+1,kyp+1,l) + 2.*f(j+1,kyp+2,l)
            f(j+1,kyp+2,l) = 0.
   30       continue
         else if (kyp.eq.1) then
            do 40 j = 2, nx
            f(j+1,kyp+1,l) = f(j+1,kyp+1,l) + 2.*f(j+1,kyp+2,l)
   40       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (kyp.eq.2) then
c        if (kl.eq.1) then
c           do 60 j = 2, nx
c           f(j+1,2,l) = f(j+1,2,l) - f(j+1,2,kl)
c  60       continue
c        endif
c     else if (kyp.eq.1) then
c        if (kl.eq.1) then
c           do 70 j = 2, nx
c           f(j+1,2,l) = f(j+1,2,l) + 2.*f(j+1,2,kl)
c  70       continue
c        endif
c        if (kll.eq.1) then
c           do 80 j = 2, nx
c           f(j+1,2,l) = f(j+1,2,l) - f(j+1,2,kll)
c  80       continue
c        endif
c        if (kr.eq.nvp) then
c           do 90 j = 2, nx
c           f(j+1,kyp+1,l) = f(j+1,kyp+1,l) - f(j+1,kyp+2,kr)
c  90       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kyp.eq.2) then
         if (kl.eq.1) then
            call MPI_IRECV(f(1,1,l),nxv,mreal,kl-1,moff+1,lgrp,msid,ierr
     1)
         endif
         if (kl.eq.0) then
            call MPI_SEND(f(1,2,l),nxv,mreal,kr-1,moff+1,lgrp,ierr)
            do 50 j = 2, nx
            f(j+1,2,l) = 0.
   50       continue
         endif
         if (kl.eq.1) then
            call MPI_WAIT(msid,istatus,ierr)
            do 60 j = 2, nx
            f(j+1,2,l) = f(j+1,2,l) - f(j+1,1,l)
            f(j+1,1,l) = 0.
   60       continue
         endif
      else if (kyp.eq.1) then
         if (kl.eq.1) then
            call MPI_IRECV(f(1,1,l),nxv,mreal,kl-1,moff+1,lgrp,msid,ierr
     1)
         endif
         if (kll.eq.1) then
            call MPI_IRECV(f(1,1,l),nxv,mreal,kll-1,moff+1,lgrp,nsid,ier
     1r)
         endif
         if (kl.eq.0) then
            call MPI_SEND(f(1,2,l),nxv,mreal,kr-1,moff+1,lgrp,ierr)
            call MPI_SEND(f(1,2,l),nxv,mreal,krr-1,moff+1,lgrp,ierr)
            do 70 j = 2, nx
            f(j+1,2,l) = 0.
   70       continue
         endif
         if (kl.eq.1) then
            call MPI_WAIT(msid,istatus,ierr)
            do 80 j = 2, nx
            f(j+1,2,l) = f(j+1,2,l) + 2.*f(j+1,1,l)
            f(j+1,1,l) = 0.
   80       continue
         endif
         if (kll.eq.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 90 j = 2, nx
            f(j+1,2,l) = f(j+1,2,l) - f(j+1,1,l)
            f(j+1,1,l) = 0.
   90       continue
         endif
         if (kr.eq.nvp) then
            call MPI_IRECV(f(1,1,l),nxv,mreal,kr-1,moff+2,lgrp,msid,ierr
     1)
         endif
         if (kr.eq.(nvp+1)) then
            call MPI_SEND(f(1,kyp+2,l),nxv,mreal,kl-1,moff+2,lgrp,ierr)
         endif
         if (kr.eq.nvp) then
            call MPI_WAIT(msid,istatus,ierr)
            do 100 j = 2, nx
            f(j+1,kyp+1,l) = f(j+1,kyp+1,l) - f(j+1,1,l)
            f(j+1,1,l) = 0.
  100       continue
         endif
      endif
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLACGUARDS2(f,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mter)
c this subroutine corrects current density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(3,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c mter = (0,1,2) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, mter
      dimension f(3,nxv,nypmx,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ks, moff, kr, krr, kl, kll, ngc, nps, j, l, m
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 240 l = 1, nblok
      kr = l + ks + 2
      krr = kr + 1
      kl = l + ks
      kll = kl - 1
c fix edges if all points are on the same processor
      if (kl.eq.0) then
         if (nyp(l).gt.2) then
            do 20 j = 2, nx
            do 10 m = 1, 3
            f(m,j+1,3,l) = f(m,j+1,3,l) + 2.*f(m,j+1,2,l)
            f(m,j+1,4,l) = f(m,j+1,4,l) - f(m,j+1,2,l)
            f(m,j+1,2,l) = 0.
   10       continue
   20       continue
         else if (nyp(l).eq.2) then
            do 40 j = 2, nx
            do 30 m = 1, 3
            f(m,j+1,3,l) = f(m,j+1,3,l) + 2.*f(m,j+1,2,l)
   30       continue
   40       continue
         endif
      endif
      if (kr.eq.(nvp+1)) then
         if (nyp(l).gt.1) then
            do 60 j = 2, nx
            do 50 m = 1, 3
            f(m,j+1,nyp(l),l) = f(m,j+1,nyp(l),l) - f(m,j+1,nyp(l)+2,l)
            f(m,j+1,nyp(l)+1,l) = f(m,j+1,nyp(l)+1,l) + 2.*f(m,j+1,nyp(l
     1)+2,l)
            f(m,j+1,nyp(l)+2,l) = 0.
   50       continue
   60       continue
         else if (nyp(l).eq.1) then
            do 80 j = 2, nx
            do 70 m = 1, 3
            f(m,j+1,nyp(l)+1,l) = f(m,j+1,nyp(l)+1,l) + 2.*f(m,j+1,nyp(l
     1)+2,l)
   70       continue
   80       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (mter.ge.2) go to 240
c     if (mter.eq.1) go to 150
c     if (kll.eq.1) then
c        if (nyp(kll).eq.1).and.(nyp(kl).eq.1)) then
c           do 100 j = 2, nx
c           do 90 m = 1, 3
c           f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,2,kll)
c  90       continue
c 100       continue
c        endif
c     endif
c     if (kr.eq.nvp) then
c        if (nyp(kr).eq.1) then
c           do 140 j = 2, nx
c           do 130 m = 1, 3
c           f(m,j+1,nyp(l)+1,l) = f(m,j+1,nyp(l)+1,l) - f(m,j+1,nyp(kr)+
c    12,kr)
c 130       continue
c 140       continue
c        endif
c     endif
c 150 if (kl.eq.1) then
c        if (nyp(kl).le.2) then
c           if (nyp(kl).eq.1) then
c              if (nyp(l).gt.1) then
c                 do 190 j = 2, nx
c                 do 180 m = 1, 3
c                 f(m,j+1,3,l) = f(m,j+1,3,l) - f(m,j+1,2,kl)
c 180             continue
c 190             continue
c              endif
c              do 210 j = 2, nx
c              do 200 m = 1, 3
c              f(m,j+1,2,l) = f(m,j+1,2,l) + 2.*f(m,j+1,2,kl)
c 200          continue
c 210          continue
c           else
c              do 230 j = 2, nx
c              do 220 m = 1, 3
c              f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,2,kl)
c 220          continue
c 230          continue
c           endif
c        endif
c     endif
c this segment is used for mpi computers
      if (mter.ge.2) go to 240
      if (mter.eq.1) go to 150
      if (kll.eq.1) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+4,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kll-1,moff+1,lgrp,nsid,ie
     1rr)
      endif
      if (kl.eq.1) then
         call MPI_SEND(nyp(l),1,mint,kr-1,moff+4,lgrp,ierr)
      endif
      if (kl.eq.0) then
         if (nyp(l).eq.1) then
            call MPI_SEND(f(1,1,2,l),3*nxv,mreal,krr-1,moff+1,lgrp,ierr)
         else
            nps = 0
            call MPI_SEND(f(1,1,2,l),nps,mreal,krr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (kll.eq.1) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               do 100 j = 2, nx
               do 90 m = 1, 3
               f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,1,l)
   90          continue
  100          continue
            endif
            do 120 j = 2, nx
            do 110 m = 1, 3
            f(m,j+1,1,l) = 0.
  110       continue
  120       continue
         endif
      endif
      if (kr.eq.nvp) then
         call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kr-1,moff+2,lgrp,msid,ier
     1r)
      endif
      if (kr.eq.(nvp+1)) then
         if (nyp(l).eq.1) then
            call MPI_SEND(f(1,1,nyp(l)+2,l),3*nxv,mreal,kl-1,moff+2,lgrp
     1,ierr)
         else
            nps = 0
            call MPI_SEND(f(1,1,nyp(l)+2,l),nps,mreal,kl-1,moff+2,lgrp,i
     1err)
         endif
      endif
      if (kr.eq.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            do 140 j = 2, nx
            do 130 m = 1, 3
            f(m,j+1,nyp(l)+1,l) = f(m,j+1,nyp(l)+1,l) - f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  130       continue
  140       continue
         endif
      endif
  150 if (kl.eq.1) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+5,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,1,l),3*nxv,mreal,kl-1,moff+1,lgrp,nsid,ier
     1r)
      endif
      if ((kl.eq.0).and.(kr.le.nvp)) then
         call MPI_SEND(nyp(l),1,mint,kr-1,moff+5,lgrp,ierr)
         if (nyp(l).le.2) then
            call MPI_SEND(f(1,1,2,l),3*nxv,mreal,kr-1,moff+1,lgrp,ierr)
            do 170 j = 2, nx
            do 160 m = 1, 3
            f(m,j+1,2,l) = 0.
  160       continue
  170       continue
         else
            nps = 0
            call MPI_SEND(f(1,1,2,l),nps,mreal,kr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (kl.eq.1) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               if (nyp(l).gt.1) then
                  do 190 j = 2, nx
                  do 180 m = 1, 3
                  f(m,j+1,3,l) = f(m,j+1,3,l) - f(m,j+1,1,l)
  180             continue
  190             continue
               endif
               do 210 j = 2, nx
               do 200 m = 1, 3
               f(m,j+1,2,l) = f(m,j+1,2,l) + 2.*f(m,j+1,1,l)
               f(m,j+1,1,l) = 0.
  200          continue
  210          continue
           else
               do 230 j = 2, nx
               do 220 m = 1, 3
               f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,1,l)
               f(m,j+1,1,l) = 0.
  220          continue
  230          continue
            endif
         endif
      endif
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLACGUARDS22(f,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mter)
c this subroutine corrects current density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(2,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c mter = (0,1,2) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, mter
      dimension f(2,nxv,nypmx,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ks, moff, kr, krr, kl, kll, ngc, nps, j, l, m
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 240 l = 1, nblok
      kr = l + ks + 2
      krr = kr + 1
      kl = l + ks
      kll = kl - 1
c fix edges if all points are on the same processor
      if (kl.eq.0) then
         if (nyp(l).gt.2) then
            do 20 j = 2, nx
            do 10 m = 1, 2
            f(m,j+1,3,l) = f(m,j+1,3,l) + 2.*f(m,j+1,2,l)
            f(m,j+1,4,l) = f(m,j+1,4,l) - f(m,j+1,2,l)
            f(m,j+1,2,l) = 0.
   10       continue
   20       continue
         else if (nyp(l).eq.2) then
            do 40 j = 2, nx
            do 30 m = 1, 2
            f(m,j+1,3,l) = f(m,j+1,3,l) + 2.*f(m,j+1,2,l)
   30       continue
   40       continue
         endif
      endif
      if (kr.eq.(nvp+1)) then
         if (nyp(l).gt.1) then
            do 60 j = 2, nx
            do 50 m = 1, 2
            f(m,j+1,nyp(l),l) = f(m,j+1,nyp(l),l) - f(m,j+1,nyp(l)+2,l)
            f(m,j+1,nyp(l)+1,l) = f(m,j+1,nyp(l)+1,l) + 2.*f(m,j+1,nyp(l
     1)+2,l)
            f(m,j+1,nyp(l)+2,l) = 0.
   50       continue
   60       continue
         else if (nyp(l).eq.1) then
            do 80 j = 2, nx
            do 70 m = 1, 2
            f(m,j+1,nyp(l)+1,l) = f(m,j+1,nyp(l)+1,l) + 2.*f(m,j+1,nyp(l
     1)+2,l)
   70       continue
   80       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (mter.ge.2) go to 240
c     if (mter.eq.1) go to 150
c     if (kll.eq.1) then
c        if (nyp(kll).eq.1).and.(nyp(kl).eq.1)) then
c           do 100 j = 2, nx
c           do 90 m = 1, 2
c           f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,2,kll)
c  90       continue
c 100       continue
c        endif
c     endif
c     if (kr.eq.nvp) then
c        if (nyp(kr).eq.1) then
c           do 140 j = 2, nx
c           do 130 m = 1, 2
c           f(m,j+1,nyp(l)+1,l) = f(m,j+1,nyp(l)+1,l) - f(m,j+1,nyp(kr)+
c    12,kr)
c 130       continue
c 140       continue
c        endif
c     endif
c 150 if (kl.eq.1) then
c        if (nyp(kl).le.2) then
c           if (nyp(kl).eq.1) then
c              if (nyp(l).gt.1) then
c                 do 190 j = 2, nx
c                 do 180 m = 1, 2
c                 f(m,j+1,3,l) = f(m,j+1,3,l) - f(m,j+1,2,kl)
c 180             continue
c 190             continue
c              endif
c              do 210 j = 2, nx
c              do 200 m = 1, 2
c              f(m,j+1,2,l) = f(m,j+1,2,l) + 2.*f(m,j+1,2,kl)
c 200          continue
c 210          continue
c           else
c              do 230 j = 2, nx
c              do 220 m = 1, 2
c              f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,2,kl)
c 220          continue
c 230          continue
c           endif
c        endif
c     endif
c this segment is used for mpi computers
      if (mter.ge.2) go to 240
      if (mter.eq.1) go to 150
      if (kll.eq.1) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+4,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kll-1,moff+1,lgrp,nsid
     1,ierr)
      endif
      if (kl.eq.1) then
         call MPI_SEND(nyp(l),1,mint,kr-1,moff+4,lgrp,ierr)
      endif
      if (kl.eq.0) then
         if (nyp(l).eq.1) then
            call MPI_SEND(f(1,1,2,l),2*nxv,mreal,krr-1,moff+1,lgrp,ierr)
         else
            nps = 0
            call MPI_SEND(f(1,1,2,l),nps,mreal,krr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (kll.eq.1) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               do 100 j = 2, nx
               do 90 m = 1, 2
               f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,1,l)
   90          continue
  100          continue
            endif
            do 120 j = 2, nx
            do 110 m = 1, 2
            f(m,j+1,1,l) = 0.
  110       continue
  120       continue
         endif
      endif
      if (kr.eq.nvp) then
         call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kr-1,moff+2,lgrp,msid,ier
     1r)
      endif
      if (kr.eq.(nvp+1)) then
         if (nyp(l).eq.1) then
            call MPI_SEND(f(1,1,nyp(l)+2,l),2*nxv,mreal,kl-1,moff+2,lgrp
     1,ierr)
         else
            nps = 0
            call MPI_SEND(f(1,1,nyp(l)+2,l),nps,mreal,kl-1,moff+2,lgrp,i
     1err)
         endif
      endif
      if (kr.eq.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            do 140 j = 2, nx
            do 130 m = 1, 2
            f(m,j+1,nyp(l)+1,l) = f(m,j+1,nyp(l)+1,l) - f(m,j+1,1,l)
            f(m,j+1,1,l) = 0.
  130       continue
  140       continue
         endif
      endif
  150 if (kl.eq.1) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+5,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,1,l),2*nxv,mreal,kl-1,moff+1,lgrp,nsid,ier
     1r)
      endif
      if ((kl.eq.0).and.(kr.le.nvp)) then
         call MPI_SEND(nyp(l),1,mint,kr-1,moff+5,lgrp,ierr)
         if (nyp(l).le.2) then
            call MPI_SEND(f(1,1,2,l),2*nxv,mreal,kr-1,moff+1,lgrp,ierr)
            do 170 j = 2, nx
            do 160 m = 1, 2
            f(m,j+1,2,l) = 0.
  160       continue
  170       continue
         else
            nps = 0
            call MPI_SEND(f(1,1,2,l),nps,mreal,kr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (kl.eq.1) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               if (nyp(l).gt.1) then
                  do 190 j = 2, nx
                  do 180 m = 1, 2
                  f(m,j+1,3,l) = f(m,j+1,3,l) - f(m,j+1,1,l)
  180             continue
  190             continue
               endif
               do 210 j = 2, nx
               do 200 m = 1, 2
               f(m,j+1,2,l) = f(m,j+1,2,l) + 2.*f(m,j+1,1,l)
               f(m,j+1,1,l) = 0.
  200          continue
  210          continue
           else
               do 230 j = 2, nx
               do 220 m = 1, 2
               f(m,j+1,2,l) = f(m,j+1,2,l) - f(m,j+1,1,l)
               f(m,j+1,1,l) = 0.
  220          continue
  230          continue
            endif
         endif
      endif
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLAGUARDS2(f,nyp,kstrt,nvp,nx,nxv,nypmx,nblok,mter)
c this subroutine corrects charge density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes three extra guard cells.
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c mter = (0,1,2) = (no,yes) pass data to next processor only
c quadratic interpolation, for distributed data
      implicit none
      real f
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok, mter
      dimension f(nxv,nypmx,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ks, moff, kr, krr, kl, kll, ngc, nps, j, l
      dimension istatus(lstat)
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 130 l = 1, nblok
      kr = l + ks + 2
      krr = kr + 1
      kl = l + ks
      kll = kl - 1
c fix edges if all points are on the same processor
      if (kl.eq.0) then
         if (nyp(l).gt.2) then
            do 10 j = 2, nx
            f(j+1,3,l) = f(j+1,3,l) + 2.*f(j+1,2,l)
            f(j+1,4,l) = f(j+1,4,l) - f(j+1,2,l)
            f(j+1,2,l) = 0.
   10       continue
         else if (nyp(l).eq.2) then
            do 20 j = 2, nx
            f(j+1,3,l) = f(j+1,3,l) + 2.*f(j+1,2,l)
   20       continue
         endif
      endif
      if (kr.eq.(nvp+1)) then
         if (nyp(l).gt.1) then
            do 30 j = 2, nx
            f(j+1,nyp(l),l) = f(j+1,nyp(l),l) - f(j+1,nyp(l)+2,l)
            f(j+1,nyp(l)+1,l) = f(j+1,nyp(l)+1,l) + 2.*f(j+1,nyp(l)+2,l)
            f(j+1,nyp(l)+2,l) = 0.
   30       continue
         else if (nyp(l).eq.1) then
            do 40 j = 2, nx
            f(j+1,nyp(l)+1,l) = f(j+1,nyp(l)+1,l) + 2.*f(j+1,nyp(l)+2,l)
   40       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (mter.ge.2) go to 130
c     if (mter.eq.1) go to 80
c     if (kll.eq.1) then
c        if (nyp(kll).eq.1).and.(nyp(kl).eq.1)) then
c           do 50 j = 2, nx
c           f(j+1,2,l) = f(j+1,2,l) - f(j+1,2,kll)
c  50       continue
c        endif
c     endif
c     if (kr.eq.nvp) then
c        if (nyp(kr).eq.1) then
c           do 70 j = 2, nx
c           f(j+1,nyp(l)+1,l) = f(j+1,nyp(l)+1,l) - f(j+1,nyp(kr)+2,kr)
c  70       continue
c        endif
c     endif
c  80 if (kl.eq.1) then
c        if (nyp(kl).le.2) then
c           if (nyp(kl).eq.1) then
c              if (nyp(l).gt.1) then
c                 do 100 j = 2, nx
c                 f(j+1,3,l) = f(j+1,3,l) - f(j+1,2,kl)
c 100             continue
c              endif
c              do 110 j = 2, nx
c              f(j+1,2,l) = f(j+1,2,l) + 2.*f(j+1,2,kl)
c 110          continue
c           else
c              do 120 j = 2, nx
c              f(j+1,2,l) = f(j+1,2,l) - f(j+1,2,kl)
c 120          continue
c           endif
c        endif
c     endif
c this segment is used for mpi computers
      if (mter.ge.2) go to 130
      if (mter.eq.1) go to 80
      if (kll.eq.1) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+4,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,l),nxv,mreal,kll-1,moff+1,lgrp,nsid,ierr)
      endif
      if (kl.eq.1) then
         call MPI_SEND(nyp(l),1,mint,kr-1,moff+4,lgrp,ierr)
      endif
      if (kl.eq.0) then
         if (nyp(l).eq.1) then
            call MPI_SEND(f(1,2,l),nxv,mreal,krr-1,moff+1,lgrp,ierr)
         else
            nps = 0
            call MPI_SEND(f(1,2,l),nps,mreal,krr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (kll.eq.1) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               do 50 j = 2, nx
               f(j+1,2,l) = f(j+1,2,l) - f(j+1,1,l)
   50          continue
            endif
            do 60 j = 2, nx
            f(j+1,1,l) = 0.
   60       continue
         endif
      endif
      if (kr.eq.nvp) then
         call MPI_IRECV(f(1,1,l),nxv,mreal,kr-1,moff+2,lgrp,msid,ierr)
      endif
      if (kr.eq.(nvp+1)) then
         if (nyp(l).eq.1) then
            call MPI_SEND(f(1,nyp(l)+2,l),nxv,mreal,kl-1,moff+2,lgrp,ier
     1r)
         else
            nps = 0
            call MPI_SEND(f(1,nyp(l)+2,l),nps,mreal,kl-1,moff+2,lgrp,ier
     1r)
         endif
      endif
      if (kr.eq.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            do 70 j = 2, nx
            f(j+1,nyp(l)+1,l) = f(j+1,nyp(l)+1,l) - f(j+1,1,l)
            f(j+1,1,l) = 0.
   70       continue
         endif
      endif
   80 if (kl.eq.1) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+5,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,l),nxv,mreal,kl-1,moff+1,lgrp,nsid,ierr)
      endif
      if ((kl.eq.0).and.(kr.le.nvp)) then
         call MPI_SEND(nyp(l),1,mint,kr-1,moff+5,lgrp,ierr)
         if (nyp(l).le.2) then
            call MPI_SEND(f(1,2,l),nxv,mreal,kr-1,moff+1,lgrp,ierr)
            do 90 j = 2, nx
            f(j+1,2,l) = 0.
   90       continue
         else
            nps = 0
            call MPI_SEND(f(1,2,l),nps,mreal,kr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (kl.eq.1) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               if (nyp(l).gt.1) then
                  do 100 j = 2, nx
                  f(j+1,3,l) = f(j+1,3,l) - f(j+1,1,l)
  100             continue
               endif
               do 110 j = 2, nx
               f(j+1,2,l) = f(j+1,2,l) + 2.*f(j+1,1,l)
               f(j+1,1,l) = 0.
  110          continue
            else
               do 120 j = 2, nx
               f(j+1,2,l) = f(j+1,2,l) - f(j+1,1,l)
               f(j+1,1,l) = 0.
  120          continue
            endif
         endif
      endif
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c no copying is done at the boundary edges.
c f(3,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes one extra guard cell.
c scr(3,j,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(3,nxv,nypmx,kblok), scr(3,nxv,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, kblok
      kr = l + ks + 2
      kl = l + ks
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nx1
c        do 10 m = 1, 3
c        scr(m,j,l) = f(m,j,kyp+1,kl)
c  10    continue
c  20    continue
c     else
c        do 40 j = 1, nx1
c        do 30 m = 1, 3
c        scr(m,j,l) = 0.
c  30    continue
c  40    continue
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,3*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,kyp+1,l),3*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 40 j = 1, nx1
         do 30 m = 1, 3
         scr(m,j,l) = 0.
   30    continue
   40    continue
      endif
c add up the guard cells
      do 60 j = 1, nx1
      do 50 m = 1, 3
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
   50 continue
   60 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARD22L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c no copying is done at the boundary edges.
c f(2,j,k,l) = real data for grid j,k in particle partition l. number of
c grids per partition is uniform and includes one extra guard cell.
c scr(2,j,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(2,nxv,nypmx,kblok), scr(2,nxv,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 80 l = 1, kblok
      kr = l + ks + 2
      kl = l + ks
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nx1
c        do 10 m = 1, 2
c        scr(m,j,l) = f(m,j,kyp+1,kl)
c  10    continue
c  20    continue
c     else
c        do 40 j = 1, nx1
c        do 30 m = 1, 2
c        scr(m,j,l) = 0.
c  30    continue
c  40    continue
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,2*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,kyp+1,l),2*nxv,mreal,kr-1,moff+1,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 40 j = 1, nx1
         do 30 m = 1, 2
         scr(m,j,l) = 0.
   30    continue
   40    continue
      endif
c add up the guard cells
      do 60 j = 1, nx1
      do 50 m = 1, 2
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
   50 continue
   60 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLAGUARD2L(f,scr,kstrt,nvp,nx,nxv,nypmx,kyp,kblok)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c no copying is done at the boundary edges.
c f(j,k,l) = real data for grid j,k in particle partition l.  the number
c grids per partition is uniform and includes one extra guard cell.
c scr(j,k) = scratch array for particle partition k
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells.
c kyp = number of complex grids in each field partition.
c kblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer kstrt, nvp, nx, nxv, nypmx, kyp, kblok
      dimension f(nxv,nypmx,kblok), scr(nxv,kblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 40 l = 1, kblok
      kr = l + ks + 2
      kl = l + ks
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nx1
c        scr(j,l) = f(j,kyp+1,kl)
c  10    continue
c     else
c        do 20 j = 1, nx1
c        scr(j,l) = 0.
c  20    continue
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,kyp+1,l),nxv,mreal,kr-1,moff+1,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nx1
         scr(j,l) = 0.
   20    continue
      endif
c add up the guard cells
      do 30 j = 1, nx1
      f(j,1,l) = f(j,1,l) + scr(j,l)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLACGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
c this subroutine adds data from guard cells in non-uniform partitions
c for vector data.  no copying is done at the boundary edges.
c f(3,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes one extra guard cell.
c scr(3,j,l) = scratch array for particle partition l
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok
      dimension f(3,nxv,nypmx,nblok), scr(3,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 90 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nx1
c        do 10 m = 1, 3
c        scr(m,j,l) = f(m,j,nyp(kl)+1,kl)
c  10    continue
c  20    continue
c     else
c        do 40 j = 1, nx1
c        do 30 m = 1, 3
c        scr(m,j,l) = 0.
c  30    continue
c  40    continue
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,3*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,nyp(l)+1,l),3*nxv,mreal,kr-1,moff+1,lgrp,ie
     1rr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 40 j = 1, nx1
         do 30 m = 1, 3
         scr(m,j,l) = 0.
   30    continue
   40    continue
      endif
c add up the guard cells
      do 60 j = 1, nx1
      do 50 m = 1, 3
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
   50 continue
   60 continue
      if (kr.le.nvp) then
         do 80 j = 1, nx1
         do 70 m = 1, 3
         f(m,j,nyp(l)+1,l) = 0.
   70    continue
   80    continue
      endif
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLACGUARD22L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
c this subroutine adds data from guard cells in non-uniform partitions
c for vector data.  no copying is done at the boundary edges.
c f(2,j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes one extra guard cell.
c scr(2,j,l) = scratch array for particle partition l
c nyp(l) = number of primary gridpoints in field partition l
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok
      dimension f(2,nxv,nypmx,nblok), scr(2,nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l, m
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 90 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nx1
c        do 10 m = 1, 2
c        scr(m,j,l) = f(m,j,nyp(kl)+1,kl)
c  10    continue
c  20    continue
c     else
c        do 40 j = 1, nx1
c        do 30 m = 1, 2
c        scr(m,j,l) = 0.
c  30    continue
c  40    continue
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,2*nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,1,nyp(l)+1,l),2*nxv,mreal,kr-1,moff+1,lgrp,ie
     1rr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 40 j = 1, nx1
         do 30 m = 1, 2
         scr(m,j,l) = 0.
   30    continue
   40    continue
      endif
c add up the guard cells
      do 60 j = 1, nx1
      do 50 m = 1, 2
      f(m,j,1,l) = f(m,j,1,l) + scr(m,j,l)
   50 continue
   60 continue
      if (kr.le.nvp) then
         do 80 j = 1, nx1
         do 70 m = 1, 2
         f(m,j,nyp(l)+1,l) = 0.
   70    continue
   80    continue
      endif
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLAGUARD2L(f,scr,nyp,kstrt,nvp,nx,nxv,nypmx,nblok)
c this subroutine adds data from guard cells in non-uniform partitions
c for scalar data.  no copying is done at the boundary edges.
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes one extra guard cell.
c scr(j,l) = scratch array for particle partition l
c nyp(l) = number of primary gridpoints in field partition l
c it is assumed the nyp(l) > 0.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, including guard cells.
c nblok = number of field partitions.
c linear interpolation, for distributed data
      implicit none
      real f, scr
      integer nyp
      integer kstrt, nvp, nx, nxv, nypmx, nblok
      dimension f(nxv,nypmx,nblok), scr(nxv,nblok)
      dimension nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer nx1, ks, moff, kl, kr, j, l
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = kstrt - 2
      moff = nypmx*nvp
c add guard cells
      do 50 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
c this segment is used for shared memory computers
c     if (kl.ge.1) then
c        do 10 j = 1, nx1
c        scr(j,l) = f(j,nyp(kl)+1,kl)
c  10    continue
c     else
c        do 20 j = 1, nx1
c        scr(j,l) = 0.
c  20    continue
c     endif
c this segment is used for mpi computers
      if (kl.ge.1) then
         call MPI_IRECV(scr,nxv,mreal,kl-1,moff+1,lgrp,msid,ierr)
      endif
      if (kr.le.nvp) then
         call MPI_SEND(f(1,nyp(l)+1,l),nxv,mreal,kr-1,moff+1,lgrp,ierr)
      endif
      if (kl.ge.1) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 20 j = 1, nx1
         scr(j,l) = 0.
   20    continue
      endif
c add up the guard cells
      do 30 j = 1, nx1
      f(j,1,l) = f(j,1,l) + scr(j,l)
   30 continue
      if (kr.le.nvp) then
         do 40 j = 1, nx1
         f(j,nyp(l)+1,l) = 0.
   40    continue
      endif
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZDBL2C(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2bl
     1ok)
c this subroutine creates an zeroed array cu2 from an array cu, so that
c a 2d real to complex transform will perform a correct convolution.
c linear interpolation for distributed data
c cu2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = second dimension of input array cu, must be >= nx
c kyp = number of data values per block in y
c kypd = third dimension of input array q, must be >= kyp
c kyp2 = third dimension of output array q2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real cu, cu2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension cu(2,nxv,kypd,kblok), cu2(2,2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer i, j, k, l, ny2, kyb, kyb2, ks, koff, moff
      integer ll, lm, k1, k2, joff
      dimension istatus(lstat)
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      ks = kstrt - 2
      moff = kypd + kyb
c copy to double array in x direction
      do 100 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = koff/kyp + 1
      koff = kyp*(l + ks)
      lm = koff/kyp2 + 1
c special case for one processor
      if (kyb2.eq.1) then
         do 30 k = 1, ny
         do 20 j = 1, nx
         do 10 i = 1, 2
         cu2(i,j,k,l) = cu(i,j,k,l)
         cu2(i,nx+j,k,l) = 0.
         cu2(i,j,ny+k,l) =  0.
         cu2(i,nx+j,ny+k,l) =  0.
   10    continue
   20    continue
   30    continue
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 60 k = 1, kyp
c        do 50 j = 1, nx
c        do 40 i = 1, 2
c        cu2(i,j,k,l) = cu(i,j,k,ll)
c  40    continue
c  50    continue
c  60    continue
c        if (kyp.lt.kyp2) then
c           do 90 k = 1, kyp
c           do 80 j = 1, nx
c           do 70 i = 1, 2
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll+1)
c  70       continue
c  80       continue
c  90       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(cu2(1,1,1,l),2*kyp*nxv,mreal,ll-1,moff+1,lgrp,ms
     1id,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),2*kyp*nxv,mreal,ll,moff+1,lg
     1rp,nsid,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(cu(1,1,1,l),2*kyp*nxv,mreal,lm-1,moff+1,lgrp,ierr
     1)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 k = 2, kyp
         k1 = kyp - k + 2
         k2 = (k1 - 1)/2 + 1
         joff = nxv*(k1 - 2*k2 + 1)
         do 50 j = 1, nxv
         do 40 i = 1, 2
         cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
   40    continue
   50    continue
   60 continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 90 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 80 j = 1, nxv
            do 70 i = 1, 2
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   70       continue
   80       continue
   90       continue
         endif
      endif
  100 continue
c zero out remainder of array
      do 160 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 150 k = 1, kyp2
      k1 = k + koff
      if (k1.le.ny) then
         do 120 j = 1, nx
         do 110 i = 1, 2
         cu2(i,j+nx,k,l) = 0.
  110    continue
  120    continue
      else
         do 140 j = 1, nx
         do 130 i = 1, 2
         cu2(i,j,k,l) = 0.
         cu2(i,j+nx,k,l) = 0.
  130    continue
  140    continue
      endif
  150 continue
  160 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZDBL2D(q,q2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2blok
     1)
c this subroutine creates an zeroed array q2 from an array q, so that
c a 2d real to complex transform will perform a correct convolution.
c linear interpolation for distributed data
c q2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = first dimension of input array q, must be >= nx
c kyp = number of data values per block in y
c kypd = second dimension of input array q, must be >= kyp
c kyp2 = second dimension of output array q2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real q, q2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension q(nxv,kypd,kblok), q2(2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer j, k, l, ny2, kyb, kyb2, ks, koff, moff
      integer ll, lm, k1, k2, joff
      dimension istatus(lstat)
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      ks = kstrt - 2
      moff = kypd + kyb
c copy to double array in x direction
      do 70 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = koff/kyp + 1
      koff = kyp*(l + ks)
      lm = koff/kyp2 + 1
c special case for one processor
      if (kyb2.eq.1) then
         do 20 k = 1, ny
         do 10 j = 1, nx
         q2(j,k,l) = q(j,k,l)
         q2(nx+j,k,l) = 0.
         q2(j,ny+k,l) =  0.
         q2(nx+j,ny+k,l) =  0.
   10    continue
   20    continue
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 40 k = 1, kyp
c        do 30 j = 1, nx
c        q2(j,k,l) = q(j,k,ll)
c  30    continue
c  40    continue
c        if (kyp.lt.kyp2) then
c           do 60 k = 1, kyp
c           do 50 j = 1, nx
c           q2(j,k+kyp,l) = q(j,k,ll+1)
c  50       continue
c  60       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(q2(1,1,l),kyp*nxv,mreal,ll-1,moff+1,lgrp,msid,ie
     1rr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(q2(1,kyp+1,l),kyp*nxv,mreal,ll,moff+1,lgrp,ns
     1id,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(q(1,1,l),kyp*nxv,mreal,lm-1,moff+1,lgrp,ierr)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 2, kyp
         k1 = kyp - k + 2
         k2 = (k1 - 1)/2 + 1
         joff = nxv*(k1 - 2*k2 + 1)
         do 30 j = 1, nxv
         q2(j,k1,l) = q2(j+joff,k2,l)
   30    continue
   40    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 60 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 50 j = 1, nxv
            q2(j,k1+kyp,l) = q2(j+joff,k2+kyp,l)
   50       continue
   60       continue
         endif
      endif
   70 continue
c zero out remainder of array
      do 110 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 100 k = 1, kyp2
      k1 = k + koff
      if (k1.le.ny) then
         do 80 j = 1, nx
         q2(j+nx,k,l) = 0.
   80    continue
      else
         do 90 j = 1, nx
         q2(j,k,l) = 0.
         q2(j+nx,k,l) = 0.
   90    continue
      endif
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZDBL2B(cu,cu2,nx,ny,kstrt,nxv,kyp,kypd,kyp2,kblok,k2bl
     1ok)
c this subroutine creates an zeroed array cu2 from an array cu, so that
c a 2d real to complex transform will perform a correct convolution.
c linear interpolation for distributed data
c cu2 array may be modified
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nxv = second dimension of input array cu, must be >= nx
c kyp = number of data values per block in y
c kypd = third dimension of input array q, must be >= kyp
c kyp2 = third dimension of output array q2, must be >= kyp2
c kblok = number of data blocks in y
c k2blok = number of data blocks in y for doubled data
      implicit none
      real cu, cu2
      integer nx, ny, kstrt, nxv, kyp, kypd, kyp2, kblok, k2blok
      dimension cu(3,nxv,kypd,kblok), cu2(3,2*nxv,kyp2,k2blok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer i, j, k, l, ny2, kyb, kyb2, ks, koff, moff
      integer ll, lm, k1, k2, joff
      dimension istatus(lstat)
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      ks = kstrt - 2
      moff = kypd + kyb
c copy to double array in x direction
      do 100 l = 1, k2blok
      koff = kyp2*(l + ks)
      ll = koff/kyp + 1
      koff = kyp*(l + ks)
      lm = koff/kyp2 + 1
c special case for one processor
      if (kyb2.eq.1) then
         do 30 k = 1, ny
         do 20 j = 1, nx
         do 10 i = 1, 3
         cu2(i,j,k,l) = cu(i,j,k,l)
         cu2(i,nx+j,k,l) = 0.
         cu2(i,j,ny+k,l) =  0.
         cu2(i,nx+j,ny+k,l) =  0.
   10    continue
   20    continue
   30    continue
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 60 k = 1, kyp
c        do 50 j = 1, nx
c        do 40 i = 1, 3
c        cu2(i,j,k,l) = cu(i,j,k,ll)
c  40    continue
c  50    continue
c  60    continue
c        if (kyp.lt.kyp2) then
c           do 90 k = 1, kyp
c           do 80 j = 1, nx
c           do 70 i = 1, 3
c           cu2(i,j,k+kyp,l) = cu(i,j,k,ll+1)
c  70       continue
c  80       continue
c  90       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(cu2(1,1,1,l),3*kyp*nxv,mreal,ll-1,moff+1,lgrp,ms
     1id,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(cu2(1,1,kyp+1,l),3*kyp*nxv,mreal,ll,moff+1,lg
     1rp,nsid,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(cu(1,1,1,l),3*kyp*nxv,mreal,lm-1,moff+1,lgrp,ierr
     1)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 k = 2, kyp
         k1 = kyp - k + 2
         k2 = (k1 - 1)/2 + 1
         joff = nxv*(k1 - 2*k2 + 1)
         do 50 j = 1, nxv
         do 40 i = 1, 3
         cu2(i,j,k1,l) = cu2(i,j+joff,k2,l)
   40    continue
   50    continue
   60 continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 90 k = 2, kyp
            k1 = kyp - k + 2
            k2 = (k1 - 1)/2 + 1
            joff = nxv*(k1 - 2*k2 + 1)
            do 80 j = 1, nxv
            do 70 i = 1, 3
            cu2(i,j,k1+kyp,l) = cu2(i,j+joff,k2+kyp,l)
   70       continue
   80       continue
   90       continue
         endif
      endif
  100 continue
c zero out remainder of array
      do 160 l = 1, k2blok
      koff = kyp2*(l + ks)
      do 150 k = 1, kyp2
      k1 = k + koff
      if (k1.le.ny) then
         do 120 j = 1, nx
         do 110 i = 1, 3
         cu2(i,j+nx,k,l) = 0.
  110    continue
  120    continue
      else
         do 140 j = 1, nx
         do 130 i = 1, 3
         cu2(i,j,k,l) = 0.
         cu2(i,j+nx,k,l) = 0.
  130    continue
  140    continue
      endif
  150 continue
  160 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr
     1,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
c this subroutine moves particles into appropriate spatial regions
c periodic boundary conditions
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c npp(l) = number of particles in partition l
c sbufl = buffer for particles being sent to lower processor
c sbufr = buffer for particles being sent to upper processor
c rbufl = buffer for particles being received from lower processor
c rbufr = buffer for particles being received from upper processor
c ihole = location of holes left in particle arrays
c jsl(idps,l) = number of particles going down in particle partition l
c jsr(idps,l) = number of particles going up in particle partition l
c jss(idps,l) = scratch array for particle partition l
c ny = system length in y direction
c kstrt = starting data block number
c nvp = number of real or virtual processors
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c idps = number of partition boundaries
c nbmax =  size of buffers for passing particles between processors
c ntmax =  size of hole array for particles leaving processors
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows
c info(5) = maximum number of particle passes required
c info(6) = total number of particles on entry
c info(7) = total number of particles on exit
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer iy
      parameter(iy=2)
      integer ierr, l, ks, iter, npr, nps, npt, kb, kl, kr, j, j1, j2
      integer i, nbsize, nter, mter, itermax
      integer msid, istatus
      integer ibflg, iwork
      real any, yt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
      any = float(ny)
      ks = kstrt - 2
      nbsize = idimp*nbmax
      iter = 2
      nter = 0
      do 5 j = 1, 7
      info(j) = 0
    5 continue
      itermax = 2000
c debugging section: count total number of particles before move
      npr = 0
      do 10 l = 1, nblok
      npr = npr + npp(l)
   10 continue
c buffer outgoing particles
   20 mter = 0
      do 50 l = 1, nblok
      kb = l + ks
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 30 j = 1, npp(l)
      yt = part(iy,j,l)
c particles going down
      if (yt.lt.edges(1,l)) then
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            if (kb.eq.0) yt = yt + any
            do 23 i = 1, idimp
            sbufl(i,jsl(1,l),l) = part(i,j,l)
   23       continue
            sbufl(iy,jsl(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = j
         else
            jss(2,l) = 1
            go to 40
         endif
c particles going up
      else if (yt.ge.edges(2,l)) then
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            if ((kb+1).eq.nvp) yt = yt - any
            do 27 i = 1, idimp
            sbufr(i,jsr(1,l),l) = part(i,j,l)
   27       continue
            sbufr(iy,jsr(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = j
         else
            jss(2,l) = 1
            go to 40
         endif
      endif
   30 continue
   40 jss(1,l) = jsl(1,l) + jsr(1,l)
   50 continue
c check for full buffer condition
      nps = 0
      do 90 l = 1, nblok
      nps = max0(nps,jss(2,l))
   90 continue
      ibflg(3) = nps
c copy particle buffers
  100 iter = iter + 2
      mter = mter + 1
      do 130 l = 1, nblok
c get particles from below and above
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
c this segment is used for shared memory computers
c     jsl(2,l) = jsr(1,kl)
c     do 110 j = 1, jsl(2,l)
c     do 105 i = 1, idimp
c     rbufl(i,j,l) = sbufr(i,j,kl)
c 105 continue
c 110 continue
c     jsr(2,l) = jsl(1,kr)
c     do 120 j = 1, jsr(2,l)
c     do 115 i = 1, idimp
c     rbufr(i,j,l) = sbufl(i,j,kr)
c 115 continue
c 120 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,l),mreal,kr-1,iter-1,lgrp,msid(3)
     1,ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,l),mreal,kl-1,iter,lgrp,msid(4),i
     1err)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,l) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,l) = nps/idimp
  130 continue
c check if particles must be passed further
      nps = 0
      do 160 l = 1, nblok
c check if any particles coming from above belong here
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 140 j = 1, jsr(2,l)
      if (rbufr(iy,j,l).lt.edges(1,l)) jsl(1,l) = jsl(1,l) + 1
      if (rbufr(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
  140 continue
      if (jsr(1,l).ne.0) write (2,*) 'Info: particles returning up'
c check if any particles coming from below belong here
      do 150 j = 1, jsl(2,l)
      if (rbufl(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
      if (rbufl(iy,j,l).lt.edges(1,l)) jss(2,l) = jss(2,l) + 1
  150 continue
      if (jss(2,l).ne.0) write (2,*) 'Info: particles returning down'
      jsl(1,l) = jsl(1,l) + jss(2,l)
      nps = max0(nps,jsl(1,l)+jsr(1,l))
  160 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
      if (nps.eq.0) go to 210
c remove particles which do not belong here
      do 200 l = 1, nblok
      kb = l + ks
c first check particles coming from above
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 170 j = 1, jsr(2,l)
      yt = rbufr(iy,j,l)
c particles going down
      if (yt.lt.edges(1,l)) then
         jsl(1,l) = jsl(1,l) + 1
         if (kb.eq.0) yt = yt + any
         rbufr(iy,j,l) = yt
         do 163 i = 1, idimp
         sbufl(i,jsl(1,l),l) = rbufr(i,j,l)
  163    continue
c particles going up, should not happen
      elseif (yt.ge.edges(2,l)) then
         jsr(1,l) = jsr(1,l) + 1
         if ((kb+1).eq.nvp) yt = yt - any
         rbufr(iy,j,l) = yt
         do 165 i = 1, idimp
         sbufr(i,jsr(1,l),l) = rbufr(i,j,l)
  165    continue
c particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 167 i = 1, idimp
         rbufr(i,jss(2,l),l) = rbufr(i,j,l)
  167    continue
      endif
  170 continue
      jsr(2,l) = jss(2,l)
c next check particles coming from below
      jss(2,l) = 0
      do 180 j = 1, jsl(2,l)
      yt = rbufl(iy,j,l)
c particles going up
      if (yt.ge.edges(2,l)) then
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            if ((kb+1).eq.nvp) yt = yt - any
            rbufl(iy,j,l) = yt
            do 173 i = 1, idimp
            sbufr(i,jsr(1,l),l) = rbufl(i,j,l)
  173       continue    
         else
            jss(2,l) = 2*npmax
            go to 190
         endif
c particles going down, should not happen
      elseif (yt.lt.edges(1,l)) then
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            if (kb.eq.0) yt = yt + any
            rbufl(iy,j,l) = yt
            do 175 i = 1, idimp
            sbufl(i,jsl(1,l),l) = rbufl(i,j,l)
  175       continue
         else
            jss(2,l) = 2*npmax
            go to 190
         endif
c particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 177 i = 1, idimp
         rbufl(i,jss(2,l),l) = rbufl(i,j,l)
  177    continue
      endif
  180 continue
  190 jsl(2,l) = jss(2,l)
  200 continue
c check if move would overflow particle array
  210 nps = 0
      npt = npmax
      do 220 l = 1, nblok
      jss(2,l) = npp(l) + jsl(2,l) + jsr(2,l) - jss(1,l)  
      nps = max0(nps,jss(2,l))
      npt = min0(npt,jss(2,l))
  220 continue
      ibflg(1) = nps
      ibflg(4) = -npt
      call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
      do 260 l = 1, nblok
c distribute incoming particles from buffers
c distribute particles coming from below into holes
      jss(2,l) = min0(jss(1,l),jsl(2,l))
      do 230 j = 1, jss(2,l)
      do 225 i = 1, idimp
      part(i,ihole(j,l),l) = rbufl(i,j,l)
  225 continue
  230 continue
      if (jss(1,l).gt.jsl(2,l)) then
         jss(2,l) = min0(jss(1,l)-jsl(2,l),jsr(2,l))
      else
         jss(2,l) = jsl(2,l) - jss(1,l)
      endif
      do 240 j = 1, jss(2,l)
c no more particles coming from below
c distribute particles coming from above into holes
      if (jss(1,l).gt.jsl(2,l)) then
         do 233 i = 1, idimp
         part(i,ihole(j+jsl(2,l),l),l) = rbufr(i,j,l)
  233    continue
      else
c no more holes
c distribute remaining particles from below into bottom
         do 237 i = 1, idimp
         part(i,j+npp(l),l) = rbufl(i,j+jss(1,l),l)
  237    continue
      endif
  240 continue
      if (jss(1,l).le.jsl(2,l)) then
         npp(l) = npp(l) + (jsl(2,l) - jss(1,l))
         jss(1,l) = jsl(2,l)
      endif
      jss(2,l) = jss(1,l) - (jsl(2,l) + jsr(2,l))
      if (jss(2,l).gt.0) then
         jss(1,l) = (jsl(2,l) + jsr(2,l))
         jsr(2,l) = jss(2,l)
      else
         jss(1,l) = jss(1,l) - jsl(2,l)
         jsr(2,l) = -jss(2,l)
      endif
      do 250 j = 1, jsr(2,l)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (jss(2,l).gt.0) then
         j1 = npp(l) - j + 1
         j2 = jss(1,l) + jss(2,l) - j + 1
         if (j1.gt.ihole(j2,l)) then
c move particle only if it is below current hole
            do 243 i = 1, idimp
            part(i,ihole(j2,l),l) = part(i,j1,l)
  243       continue
         endif
      else
c no more holes
c distribute remaining particles from above into bottom
         do 247 i = 1, idimp
         part(i,j+npp(l),l) = rbufr(i,j+jss(1,l),l)
  247    continue
      endif
  250 continue
      if (jss(2,l).gt.0) then
         npp(l) = npp(l) - jsr(2,l)
      else
         npp(l) = npp(l) + jsr(2,l)
      endif
      jss(1,l) = 0
  260 continue
c check if any particles have to be passed further
      info(5) = max0(info(5),mter)
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 100
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         go to 280
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(4) = nter
         go to 20
      endif
c debugging section: count total number of particles after move
      nps = 0
      do 270 l = 1, nblok
      nps = nps + npp(l)
  270 continue
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      info(6) = ibflg(1)
      info(7) = ibflg(2)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         info(1) = 1
      endif
c information
  280 if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr
     1,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,maskp,inf
     2o)
c this subroutine moves particles into appropriate spatial regions
c periodic boundary conditions
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c npp(l) = number of particles in partition l
c sbufl = buffer for particles being sent to lower processor
c sbufr = buffer for particles being sent to upper processor
c rbufl = buffer for particles being received from lower processor
c rbufr = buffer for particles being received from upper processor
c ihole = location of holes left in particle arrays
c jsl(idps,l) = number of particles going down in particle partition l
c jsr(idps,l) = number of particles going up in particle partition l
c jss(idps,l) = scratch array for particle partition l
c ny = system length in y direction
c kstrt = starting data block number
c nvp = number of real or virtual processors
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c idps = number of partition boundaries
c nbmax =  size of buffers for passing particles between processors
c ntmax =  size of hole array for particles leaving processors
c maskp = scratch array for particle addresses
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows
c info(5) = maximum number of particle passes required
c info(6) = total number of particles on entry
c info(7) = total number of particles on exit
c optimized for vector processor
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, maskp, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok), maskp(npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer iy
      parameter(iy=2)
      integer ierr, l, ks, iter, npr, nps, npt, kb, kl, kr, j, j1, j2
      integer i, nbsize, nter, mter, itermax
      integer msid, istatus
      integer ibflg, iwork
      real any, yt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
      any = float(ny)
      ks = kstrt - 2
      nbsize = idimp*nbmax
      iter = 2
      nter = 0
      do 5 j = 1, 7
      info(j) = 0
    5 continue
      itermax = 2000
c debugging section: count total number of particles before move
      npr = 0
      do 10 l = 1, nblok
      npr = npr + npp(l)
   10 continue
c buffer outgoing particles
   20 mter = 0
      do 80 l = 1, nblok
      jss(1,l) = 0
      jss(2,l) = 0
c find mask function for particles out of bounds
      do 30 j = 1, npp(l)
      yt = part(iy,j,l)
      if ((yt.ge.edges(2,l)).or.(yt.lt.edges(1,l))) then
         jss(1,l) = jss(1,l) + 1
         maskp(j,l) = 1
      else
         maskp(j,l) = 0
      endif
   30 continue
c set flag if hole buffer would overflow
      if (jss(1,l).gt.ntmax) then
         jss(1,l) = ntmax
         jss(2,l) = 1
      endif
c accumulate location of holes
      do 40 j = 2, npp(l)
      maskp(j,l) = maskp(j,l) + maskp(j-1,l)
   40 continue
c store addresses of particles out of bounds
      do 50 j = 2, npp(l)
      if ((maskp(j,l).gt.maskp(j-1,l)).and.(maskp(j,l).le.ntmax)) then
         ihole(maskp(j,l),l) = j
      endif
   50 continue
      if (maskp(1,l).gt.0) ihole(1,l) = 1
      kb = l + ks
      jsl(1,l) = 0
      jsr(1,l) = 0
c load particle buffers
      do 60 j = 1, jss(1,l)
      yt = part(iy,ihole(j,l),l)
c particles going down
      if (yt.lt.edges(1,l)) then
         if (kb.eq.0) yt = yt + any
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            do 53 i = 1, idimp
            sbufl(i,jsl(1,l),l) = part(i,ihole(j,l),l)
   53       continue
            sbufl(iy,jsl(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = ihole(j,l)
         else
            jss(2,l) = 1
c           go to 70
         endif
c particles going up
      else
         if ((kb+1).eq.nvp) yt = yt - any
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            do 57 i = 1, idimp
            sbufr(i,jsr(1,l),l) = part(i,ihole(j,l),l)
   57       continue
            sbufr(iy,jsr(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = ihole(j,l)
         else
            jss(2,l) = 1
c           go to 70
         endif
      endif
   60 continue
   70 jss(1,l) = jsl(1,l) + jsr(1,l)
   80 continue
c check for full buffer condition
      nps = 0
      do 90 l = 1, nblok
      nps = max0(nps,jss(2,l))
   90 continue
      ibflg(3) = nps
c copy particle buffers
  100 iter = iter + 2
      mter = mter + 1
      do 130 l = 1, nblok
c get particles from below and above
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
c this segment is used for shared memory computers
c     jsl(2,l) = jsr(1,kl)
c     do 110 j = 1, jsl(2,l)
c     do 105 i = 1, idimp
c     rbufl(i,j,l) = sbufr(i,j,kl)
c 105 continue
c 110 continue
c     jsr(2,l) = jsl(1,kr)
c     do 120 j = 1, jsr(2,l)
c     do 115 i = 1, idimp
c     rbufr(i,j,l) = sbufl(i,j,kr)
c 115 continue
c 120 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,l),mreal,kr-1,iter-1,lgrp,msid(3)
     1,ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,l),mreal,kl-1,iter,lgrp,msid(4),i
     1err)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,l) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,l) = nps/idimp
  130 continue
c check if particles must be passed further
      nps = 0
      do 160 l = 1, nblok
c check if any particles coming from above belong here
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 140 j = 1, jsr(2,l)
      if (rbufr(iy,j,l).lt.edges(1,l)) jsl(1,l) = jsl(1,l) + 1
      if (rbufr(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
  140 continue
      if (jsr(1,l).ne.0) write (2,*) 'Info: particles returning up'
c check if any particles coming from below belong here
      do 150 j = 1, jsl(2,l)
      if (rbufl(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
      if (rbufl(iy,j,l).lt.edges(1,l)) jss(2,l) = jss(2,l) + 1
  150 continue
      if (jss(2,l).ne.0) write (2,*) 'Info: particles returning down'
      jsl(1,l) = jsl(1,l) + jss(2,l)
      nps = max0(nps,jsl(1,l)+jsr(1,l))
  160 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
      if (nps.eq.0) go to 210
c remove particles which do not belong here
      do 200 l = 1, nblok
      kb = l + ks
c first check particles coming from above
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 170 j = 1, jsr(2,l)
      yt = rbufr(iy,j,l)
c particles going down
      if (yt.lt.edges(1,l)) then
         jsl(1,l) = jsl(1,l) + 1
         if (kb.eq.0) yt = yt + any
         rbufr(iy,j,l) = yt
         do 163 i = 1, idimp
         sbufl(i,jsl(1,l),l) = rbufr(i,j,l)
  163    continue
c particles going up, should not happen
      elseif (yt.ge.edges(2,l)) then
         jsr(1,l) = jsr(1,l) + 1
         if ((kb+1).eq.nvp) yt = yt - any
         rbufr(iy,j,l) = yt
         do 165 i = 1, idimp
         sbufr(i,jsr(1,l),l) = rbufr(i,j,l)
  165    continue
c particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 167 i = 1, idimp
         rbufr(i,jss(2,l),l) = rbufr(i,j,l)
  167    continue
      endif
  170 continue
      jsr(2,l) = jss(2,l)
c next check particles coming from below
      jss(2,l) = 0
      do 180 j = 1, jsl(2,l)
      yt = rbufl(iy,j,l)
c particles going up
      if (yt.ge.edges(2,l)) then
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            if ((kb+1).eq.nvp) yt = yt - any
            rbufl(iy,j,l) = yt
            do 173 i = 1, idimp
            sbufr(i,jsr(1,l),l) = rbufl(i,j,l)
  173       continue
         else
            jss(2,l) = 2*npmax
            go to 190
         endif
c particles going down, should not happen
      elseif (yt.lt.edges(1,l)) then
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            if (kb.eq.0) yt = yt + any
            rbufl(iy,j,l) = yt
            do 175 i = 1, idimp
            sbufl(i,jsl(1,l),l) = rbufl(i,j,l)
  175       continue
         else
            jss(2,l) = 2*npmax
            go to 190
         endif
c particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 177 i = 1, idimp
         rbufl(i,jss(2,l),l) = rbufl(i,j,l)
  177    continue
      endif
  180 continue
  190 jsl(2,l) = jss(2,l)
  200 continue
c check if move would overflow particle array
  210 nps = 0
      npt = npmax
      do 220 l = 1, nblok
      jss(2,l) = npp(l) + jsl(2,l) + jsr(2,l) - jss(1,l)  
      nps = max0(nps,jss(2,l))
      npt = min0(npt,jss(2,l))
  220 continue
      ibflg(1) = nps
      ibflg(4) = -npt
      call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
      do 260 l = 1, nblok
c distribute incoming particles from buffers
c distribute particles coming from below into holes
      jss(2,l) = min0(jss(1,l),jsl(2,l))
      do 230 j = 1, jss(2,l)
      do 225 i = 1, idimp
      part(i,ihole(j,l),l) = rbufl(i,j,l)
  225 continue
  230 continue
      if (jss(1,l).gt.jsl(2,l)) then
         jss(2,l) = min0(jss(1,l)-jsl(2,l),jsr(2,l))
      else
         jss(2,l) = jsl(2,l) - jss(1,l)
      endif
      do 240 j = 1, jss(2,l)
c no more particles coming from below
c distribute particles coming from above into holes
      if (jss(1,l).gt.jsl(2,l)) then
         do 233 i = 1, idimp
         part(i,ihole(j+jsl(2,l),l),l) = rbufr(i,j,l)
  233    continue
      else
c no more holes
c distribute remaining particles from below into bottom
         do 237 i = 1, idimp
         part(i,j+npp(l),l) = rbufl(i,j+jss(1,l),l)
  237    continue
      endif
  240 continue
      if (jss(1,l).le.jsl(2,l)) then
         npp(l) = npp(l) + (jsl(2,l) - jss(1,l))
         jss(1,l) = jsl(2,l)
      endif
      jss(2,l) = jss(1,l) - (jsl(2,l) + jsr(2,l))
      if (jss(2,l).gt.0) then
         jss(1,l) = (jsl(2,l) + jsr(2,l))
         jsr(2,l) = jss(2,l)
      else
         jss(1,l) = jss(1,l) - jsl(2,l)
         jsr(2,l) = -jss(2,l)
      endif
      do 250 j = 1, jsr(2,l)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (jss(2,l).gt.0) then
         j1 = npp(l) - j + 1
         j2 = jss(1,l) + jss(2,l) - j + 1
         if (j1.gt.ihole(j2,l)) then
c move particle only if it is below current hole
            do 243 i = 1, idimp
            part(i,ihole(j2,l),l) = part(i,j1,l)
  243       continue
         endif
      else
c no more holes
c distribute remaining particles from above into bottom
         do 247 i = 1, idimp
         part(i,j+npp(l),l) = rbufr(i,j+jss(1,l),l)
  247    continue
      endif
  250 continue
      if (jss(2,l).gt.0) then
         npp(l) = npp(l) - jsr(2,l)
      else
         npp(l) = npp(l) + jsr(2,l)
      endif
      jss(1,l) = 0
  260 continue
c check if any particles have to be passed further
      info(5) = max0(info(5),mter)
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 100
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         go to 280
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(4) = nter
         go to 20
      endif
c debugging section: count total number of particles after move
      nps = 0
      do 270 l = 1, nblok
      nps = nps + npp(l)
  270 continue
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      info(6) = ibflg(1)
      info(7) = ibflg(2)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         info(1) = 1
      endif
c information
  280 if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js
     1r,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
c wrapper function for particle manager
c info(6) = total number of particles on entry
c info(7) = total number of particles on exit
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
c local data
      integer j, l, npr, nps, nter
      real tf
      double precision dtime
      integer ibflg, iwork
      dimension ibflg(2), iwork(2)
      do 10 j = 1, 7
      info(j) = 0
   10 continue
      th = 0.0
c debugging section: count total number of particles before move
      npr = 0
      do 20 l = 1, nblok
      npr = npr + npp(l)
   20 continue
c find outgoing particles
   30 nter = info(4)
      call PWTIMERA(-1,tf,dtime)
      call PMOVEH2(part,edges,npp,ihole,jss,idimp,npmax,nblok,idps,ntmax
     1)
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVES2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl,
     1jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
c particle overflow error
      if (info(1).gt.0) return
c buffer overflowed and more particles remain to be checked
      if (info(4).gt.nter) go to 30
c iteration overflow
      if (info(1).lt.0) go to 50
c debugging section: count total number of particles after move
      nps = 0
      do 40 l = 1, nblok
      nps = nps + npp(l)
   40 continue
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      info(6) = ibflg(1)
      info(7) = ibflg(2)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         info(1) = 1
      endif
c information
   50 nter = info(4)
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WPXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js
     1r,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,maskp
     2,info)
c wrapper function for particle manager
c info(6) = total number of particles on entry
c info(7) = total number of particles on exit
c optimized for vector processor
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, ihole, jsr, jsl, jss, maskp, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok), maskp(npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
c local data
      integer j, l, npr, nps, nter
      real tf
      double precision dtime
      integer ibflg, iwork
      dimension ibflg(2), iwork(2)
      do 10 j = 1, 7
      info(j) = 0
   10 continue
      th = 0.0
c debugging section: count total number of particles before move
      npr = 0
      do 20 l = 1, nblok
      npr = npr + npp(l)
   20 continue
c find outgoing particles
   30 nter = info(4)
      call PWTIMERA(-1,tf,dtime)
      call PMOVEHX2(part,edges,npp,ihole,jss,idimp,npmax,nblok,idps,ntma
     1x,maskp)
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVES2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl,
     1jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
c particle overflow error
      if (info(1).gt.0) return
c buffer overflowed and more particles remain to be checked
      if (info(4).gt.nter) go to 30
c iteration overflow
      if (info(1).lt.0) go to 50
c debugging section: count total number of particles after move
      nps = 0
      do 40 l = 1, nblok
      nps = nps + npp(l)
   40 continue
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      info(6) = ibflg(1)
      info(7) = ibflg(2)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         info(1) = 1
      endif
c information
   50 nter = info(4)
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WPMOVES2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j
     1sr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info
     2)
c wrapper function for particle manager
c info(5) = maximum number of particle passes required
c info(5) must be set on entry
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
c local data
      integer j
      real tf
      double precision dtime
      do 10 j = 1, 4
      info(j) = 0
   10 continue
      th = 0.0
c find outgoing particles
      call PWTIMERA(-1,tf,dtime)
      call PMOVEH2(part,edges,npp,ihole,jss,idimp,npmax,nblok,idps,ntmax
     1)
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVESS2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl
     1,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
      return
      end
c-----------------------------------------------------------------------
      subroutine WPXMOVS2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j
     1sr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,mask
     2p,info)
c wrapper function for particle manager
c info(5) = maximum number of particle passes required
c info(5) must be set on entry
c optimized for vector processor
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, ihole, jsr, jsl, jss, maskp, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok), maskp(npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
c local data
      integer j
      real tf
      double precision dtime
      do 10 j = 1, 4
      info(j) = 0
   10 continue
      th = 0.0
c find outgoing particles
      call PWTIMERA(-1,tf,dtime)
      call PMOVEHX2(part,edges,npp,ihole,jss,idimp,npmax,nblok,idps,ntma
     1x,maskp)
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVESS2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl
     1,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVEH2(part,edges,npp,ihole,jss,idimp,npmax,nblok,idps
     1,ntmax)
c this subroutine determines list of particles which are leaving this
c processor
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c npp(l) = number of particles in partition l
c ihole = location of holes left in particle arrays
c jss(1,l) = number of particles leaving, for particle partition l
c jss(2,l) = (0,1) = (no,yes) ihole overflowed, for particle partition l
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
      implicit none
      real part, edges
      integer npp, ihole, jss
      integer idimp, npmax, nblok, idps, ntmax
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
c local data
c iy = partitioned co-ordinate
      integer iy
      parameter(iy=2)
      integer j, l
      real yt
c find particles out of bounds
      do 20 l = 1, nblok
      jss(1,l) = 0
      jss(2,l) = 0
      do 10 j = 1, npp(l)
      yt = part(iy,j,l)
      if ((yt.ge.edges(2,l)).or.(yt.lt.edges(1,l))) then
         if (jss(1,l).lt.ntmax) then
            jss(1,l) = jss(1,l) + 1
            ihole(jss(1,l),l) = j
         else
            jss(2,l) = 1
            go to 20
         endif
      endif
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVEHX2(part,edges,npp,ihole,jss,idimp,npmax,nblok,idp
     1s,ntmax,maskp)
c this subroutine determines list of particles which are leaving this
c processor
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c npp(l) = number of particles in partition l
c ihole = location of holes left in particle arrays
c jss(1,l) = number of particles leaving, for particle partition l
c jss(2,l) = (0,1) = (no,yes) ihole overflowed, for particle partition l
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c idps = number of partition boundaries
c ntmax =  size of hole array for particles leaving processors
c maskp = scratch array for particle addresses
c optimized for vector processor
      implicit none
      real part, edges
      integer npp, ihole, jss, maskp
      integer idimp, npmax, nblok, idps, ntmax
      dimension part(idimp,npmax,nblok), maskp(npmax,nblok)
      dimension edges(idps,nblok), npp(nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
c local data
c iy = partitioned co-ordinate
      integer iy
      parameter(iy=2)
      integer j, l
      real yt
c find particles out of bounds
      do 40 l = 1, nblok
      jss(1,l) = 0
      jss(2,l) = 0
c find mask function for particles out of bounds
      do 10 j = 1, npp(l)
      yt = part(iy,j,l)
      if ((yt.ge.edges(2,l)).or.(yt.lt.edges(1,l))) then
         jss(1,l) = jss(1,l) + 1
         maskp(j,l) = 1
      else
         maskp(j,l) = 0
      endif
   10 continue
c set flag if hole buffer would overflow
      if (jss(1,l).gt.ntmax) then
         jss(1,l) = ntmax
         jss(2,l) = 1
      endif
c accumulate location of holes
      do 20 j = 2, npp(l)
      maskp(j,l) = maskp(j,l) + maskp(j-1,l)
   20 continue
c store addresses of particles out of bounds
      do 30 j = 2, npp(l)
      if ((maskp(j,l).gt.maskp(j-1,l)).and.(maskp(j,l).le.ntmax)) then
         ihole(maskp(j,l),l) = j
      endif
   30 continue
      if (maskp(1,l).gt.0) ihole(1,l) = 1
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVES2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js
     1r,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
c this subroutine moves particles into appropriate spatial regions
c periodic boundary conditions
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c npp(l) = number of particles in partition l
c sbufl = buffer for particles being sent to lower processor
c sbufr = buffer for particles being sent to upper processor
c rbufl = buffer for particles being received from lower processor
c rbufr = buffer for particles being received from upper processor
c ihole = location of holes left in particle arrays
c jsl(idps,l) = number of particles going down in particle partition l
c jsr(idps,l) = number of particles going up in particle partition l
c jss(idps,l) = scratch array for particle partition l
c ny = system length in y direction
c kstrt = starting data block number
c nvp = number of real or virtual processors
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c idps = number of partition boundaries
c nbmax =  size of buffers for passing particles between processors
c ntmax =  size of hole array for particles leaving processors
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows
c info(5) = maximum number of particle passes required
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
c iy = partitioned co-ordinate
      integer iy
      parameter(iy=2)
      integer ierr, l, ks, iter, nps, npt, kb, kl, kr, j, j1, j2, i
      integer nbsize, nter, mter, itermax
      integer msid, istatus
      integer ibflg, iwork
      real any, yt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
      any = float(ny)
      ks = kstrt - 2
      nbsize = idimp*nbmax
      iter = 2
      nter = info(4)
      itermax = 2000
      mter = 0
c buffer outgoing particles
      do 50 l = 1, nblok
      kb = l + ks
      jsl(1,l) = 0
      jsr(1,l) = 0
c load particle buffers
      do 30 j = 1, jss(1,l)
      yt = part(iy,ihole(j,l),l)
c particles going down
      if (yt.lt.edges(1,l)) then
         if (kb.eq.0) yt = yt + any
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            do 10 i = 1, idimp
            sbufl(i,jsl(1,l),l) = part(i,ihole(j,l),l)
   10       continue
            sbufl(iy,jsl(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = ihole(j,l)
         else
            jss(2,l) = 1
c           go to 40
         endif
c particles going up
      else
         if ((kb+1).eq.nvp) yt = yt - any
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            do 20 i = 1, idimp
            sbufr(i,jsr(1,l),l) = part(i,ihole(j,l),l)
   20       continue
            sbufr(iy,jsr(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = ihole(j,l)
         else
            jss(2,l) = 1
c           go to 40
         endif
      endif
   30 continue
   40 jss(1,l) = jsl(1,l) + jsr(1,l)
   50 continue
c check for full buffer condition
      nps = 0
      do 60 l = 1, nblok
      nps = max0(nps,jss(2,l))
   60 continue
      ibflg(3) = nps
c copy particle buffers
   70 iter = iter + 2
      mter = mter + 1
      do 120 l = 1, nblok
c get particles from below and above
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
c this segment is used for shared memory computers
c     jsl(2,l) = jsr(1,kl)
c     do 90 j = 1, jsl(2,l)
c     do 80 i = 1, idimp
c     rbufl(i,j,l) = sbufr(i,j,kl)
c  80 continue
c  90 continue
c     jsr(2,l) = jsl(1,kr)
c     do 110 j = 1, jsr(2,l)
c     do 100 i = 1, idimp
c     rbufr(i,j,l) = sbufl(i,j,kr)
c 100 continue
c 110 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,l),mreal,kr-1,iter-1,lgrp,msid(3)
     1,ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,l),mreal,kl-1,iter,lgrp,msid(4),i
     1err)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,l) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,l) = nps/idimp
  120 continue
c check if particles must be passed further
      nps = 0
      do 150 l = 1, nblok
c check if any particles coming from above belong here
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 130 j = 1, jsr(2,l)
      if (rbufr(iy,j,l).lt.edges(1,l)) jsl(1,l) = jsl(1,l) + 1
      if (rbufr(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
  130 continue
      if (jsr(1,l).ne.0) write (2,*) 'Info: particles returning up'
c check if any particles coming from below belong here
      do 140 j = 1, jsl(2,l)
      if (rbufl(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
      if (rbufl(iy,j,l).lt.edges(1,l)) jss(2,l) = jss(2,l) + 1
  140 continue
      if (jss(2,l).ne.0) write (2,*) 'Info: particles returning down'
      jsl(1,l) = jsl(1,l) + jss(2,l)
      nps = max0(nps,jsl(1,l)+jsr(1,l))
  150 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
      if (nps.eq.0) go to 260
c remove particles which do not belong here
      do 250 l = 1, nblok
      kb = l + ks
c first check particles coming from above
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 190 j = 1, jsr(2,l)
      yt = rbufr(iy,j,l)
c particles going down
      if (yt.lt.edges(1,l)) then
         jsl(1,l) = jsl(1,l) + 1
         if (kb.eq.0) yt = yt + any
         rbufr(iy,j,l) = yt
         do 160 i = 1, idimp
         sbufl(i,jsl(1,l),l) = rbufr(i,j,l)
  160    continue
c particles going up, should not happen
      elseif (yt.ge.edges(2,l)) then
         jsr(1,l) = jsr(1,l) + 1
         if ((kb+1).eq.nvp) yt = yt - any
         rbufr(iy,j,l) = yt
         do 170 i = 1, idimp
         sbufr(i,jsr(1,l),l) = rbufr(i,j,l)
  170    continue
c particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 180 i = 1, idimp
         rbufr(i,jss(2,l),l) = rbufr(i,j,l)
  180    continue
      endif
  190 continue
      jsr(2,l) = jss(2,l)
c next check particles coming from below
      jss(2,l) = 0
      do 230 j = 1, jsl(2,l)
      yt = rbufl(iy,j,l)
c particles going up
      if (yt.ge.edges(2,l)) then
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            if ((kb+1).eq.nvp) yt = yt - any
            rbufl(iy,j,l) = yt
            do 200 i = 1, idimp
            sbufr(i,jsr(1,l),l) = rbufl(i,j,l)
  200       continue
         else
            jss(2,l) = 2*npmax
            go to 240
         endif
c particles going down, should not happen
      elseif (yt.lt.edges(1,l)) then
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            if (kb.eq.0) yt = yt + any
            rbufl(iy,j,l) = yt
            do 210 i = 1, idimp
            sbufl(i,jsl(1,l),l) = rbufl(i,j,l)
  210       continue
         else
            jss(2,l) = 2*npmax
            go to 240
         endif
c particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 220 i = 1, idimp
         rbufl(i,jss(2,l),l) = rbufl(i,j,l)
  220    continue
      endif
  230 continue
  240 jsl(2,l) = jss(2,l)
  250 continue
c check if move would overflow particle array
  260 nps = 0
      npt = npmax
      do 270 l = 1, nblok
      jss(2,l) = npp(l) + jsl(2,l) + jsr(2,l) - jss(1,l)  
      nps = max0(nps,jss(2,l))
      npt = min0(npt,jss(2,l))
  270 continue
      ibflg(1) = nps
      ibflg(4) = -npt
      call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
      do 360 l = 1, nblok
c distribute incoming particles from buffers
c distribute particles coming from below into holes
      jss(2,l) = min0(jss(1,l),jsl(2,l))
      do 290 j = 1, jss(2,l)
      do 280 i = 1, idimp
      part(i,ihole(j,l),l) = rbufl(i,j,l)
  280 continue
  290 continue
      if (jss(1,l).gt.jsl(2,l)) then
         jss(2,l) = min0(jss(1,l)-jsl(2,l),jsr(2,l))
      else
         jss(2,l) = jsl(2,l) - jss(1,l)
      endif
      do 320 j = 1, jss(2,l)
c no more particles coming from below
c distribute particles coming from above into holes
      if (jss(1,l).gt.jsl(2,l)) then
         do 300 i = 1, idimp
         part(i,ihole(j+jsl(2,l),l),l) = rbufr(i,j,l)
  300    continue
      else
c no more holes
c distribute remaining particles from below into bottom
         do 310 i = 1, idimp
         part(i,j+npp(l),l) = rbufl(i,j+jss(1,l),l)
  310    continue
      endif
  320 continue
      if (jss(1,l).le.jsl(2,l)) then
         npp(l) = npp(l) + (jsl(2,l) - jss(1,l))
         jss(1,l) = jsl(2,l)
      endif
      jss(2,l) = jss(1,l) - (jsl(2,l) + jsr(2,l))
      if (jss(2,l).gt.0) then
         jss(1,l) = (jsl(2,l) + jsr(2,l))
         jsr(2,l) = jss(2,l)
      else
         jss(1,l) = jss(1,l) - jsl(2,l)
         jsr(2,l) = -jss(2,l)
      endif
      do 350 j = 1, jsr(2,l)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (jss(2,l).gt.0) then
         j1 = npp(l) - j + 1
         j2 = jss(1,l) + jss(2,l) - j + 1
         if (j1.gt.ihole(j2,l)) then
c move particle only if it is below current hole
            do 330 i = 1, idimp
            part(i,ihole(j2,l),l) = part(i,j1,l)
  330       continue
         endif
      else
c no more holes
c distribute remaining particles from above into bottom
         do 340 i = 1, idimp
         part(i,j+npp(l),l) = rbufr(i,j+jss(1,l),l)
  340    continue
      endif
  350 continue
      if (jss(2,l).gt.0) then
         npp(l) = npp(l) - jsr(2,l)
      else
         npp(l) = npp(l) + jsr(2,l)
      endif
      jss(1,l) = 0
  360 continue
c check if any particles have to be passed further
      info(5) = max0(info(5),mter)
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 70
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(4) = nter
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVESS2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j
     1sr,jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
c this subroutine moves particles into appropriate spatial regions
c periodic boundary conditions
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c npp(l) = number of particles in partition l
c sbufl = buffer for particles being sent to lower processor
c sbufr = buffer for particles being sent to upper processor
c rbufl = buffer for particles being received from lower processor
c rbufr = buffer for particles being received from upper processor
c ihole = location of holes left in particle arrays
c jsl(idps,l) = number of particles going down in particle partition l
c jsr(idps,l) = number of particles going up in particle partition l
c jss(idps,l) = scratch array for particle partition l
c ny = system length in y direction
c kstrt = starting data block number
c nvp = number of real or virtual processors
c idimp = size of phase space = 4
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions.
c idps = number of partition boundaries
c nbmax =  size of buffers for passing particles between processors
c ntmax =  size of hole array for particles leaving processors
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows
c info(5) = maximum number of particle passes required
c info(5) must be set on entry
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
c iy = partitioned co-ordinate
      integer iy
      parameter(iy=2)
      integer ierr, l, ks, iter, nps, npt, kb, kl, kr, j, j1, j2, i
      integer nbsize, nter, mter, itermax
      integer msid, istatus
      integer ibflg
      real any, yt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4)
      any = float(ny)
      ks = kstrt - 2
      nbsize = idimp*nbmax
      iter = 2
c     nter = info(4)
      nter = 0
      itermax = 2000
      mter = 0
c buffer outgoing particles
      do 50 l = 1, nblok
      kb = l + ks
      jsl(1,l) = 0
      jsr(1,l) = 0
c load particle buffers
      do 30 j = 1, jss(1,l)
      yt = part(iy,ihole(j,l),l)
c particles going down
      if (yt.lt.edges(1,l)) then
         if (kb.eq.0) yt = yt + any
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            do 10 i = 1, idimp
            sbufl(i,jsl(1,l),l) = part(i,ihole(j,l),l)
   10       continue
            sbufl(iy,jsl(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = ihole(j,l)
         else
            jss(2,l) = 1
c           go to 40
         endif
c particles going up
      else
         if ((kb+1).eq.nvp) yt = yt - any
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            do 20 i = 1, idimp
            sbufr(i,jsr(1,l),l) = part(i,ihole(j,l),l)
   20       continue
            sbufr(iy,jsr(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = ihole(j,l)
         else
            jss(2,l) = 1
c           go to 40
         endif
      endif
   30 continue
   40 jss(1,l) = jsl(1,l) + jsr(1,l)
   50 continue
c check for full buffer condition
      nps = 0
      do 60 l = 1, nblok
      nps = max0(nps,jss(2,l))
   60 continue
      ibflg(3) = nps
c copy particle buffers
   70 iter = iter + 2
      mter = mter + 1
      do 120 l = 1, nblok
c get particles from below and above
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
c this segment is used for shared memory computers
c     jsl(2,l) = jsr(1,kl)
c     do 90 j = 1, jsl(2,l)
c     do 80 i = 1, idimp
c     rbufl(i,j,l) = sbufr(i,j,kl)
c  80 continue
c  90 continue
c     jsr(2,l) = jsl(1,kr)
c     do 110 j = 1, jsr(2,l)
c     do 100 i = 1, idimp
c     rbufr(i,j,l) = sbufl(i,j,kr)
c 100 continue
c 110 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,l),mreal,kr-1,iter-1,lgrp,msid(3)
     1,ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,l),mreal,kl-1,iter,lgrp,msid(4),i
     1err)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,l) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,l) = nps/idimp
  120 continue
c check if particles must be passed further
      nps = 0
      do 150 l = 1, nblok
c check if any particles coming from above belong here
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 130 j = 1, jsr(2,l)
      if (rbufr(iy,j,l).lt.edges(1,l)) jsl(1,l) = jsl(1,l) + 1
      if (rbufr(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
  130 continue
      if (jsr(1,l).ne.0) write (2,*) 'Info: particles returning up'
c check if any particles coming from below belong here
      do 140 j = 1, jsl(2,l)
      if (rbufl(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
      if (rbufl(iy,j,l).lt.edges(1,l)) jss(2,l) = jss(2,l) + 1
  140 continue
      if (jss(2,l).ne.0) write (2,*) 'Info: particles returning down'
      jsl(1,l) = jsl(1,l) + jss(2,l)
      nps = max0(nps,jsl(1,l)+jsr(1,l))
  150 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
      if ((nps.eq.0).or.(mter.eq.info(5))) go to 260
c remove particles which do not belong here
      do 250 l = 1, nblok
      kb = l + ks
c first check particles coming from above
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 190 j = 1, jsr(2,l)
      yt = rbufr(iy,j,l)
c particles going down
      if (yt.lt.edges(1,l)) then
         jsl(1,l) = jsl(1,l) + 1
         if (kb.eq.0) yt = yt + any
         rbufr(iy,j,l) = yt
         do 160 i = 1, idimp
         sbufl(i,jsl(1,l),l) = rbufr(i,j,l)
  160    continue
c particles going up, should not happen
      elseif (yt.ge.edges(2,l)) then
         jsr(1,l) = jsr(1,l) + 1
         if ((kb+1).eq.nvp) yt = yt - any
         rbufr(iy,j,l) = yt
         do 170 i = 1, idimp
         sbufr(i,jsr(1,l),l) = rbufr(i,j,l)
  170    continue
c particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 180 i = 1, idimp
         rbufr(i,jss(2,l),l) = rbufr(i,j,l)
  180    continue
      endif
  190 continue
      jsr(2,l) = jss(2,l)
c next check particles coming from below
      jss(2,l) = 0
      do 230 j = 1, jsl(2,l)
      yt = rbufl(iy,j,l)
c particles going up
      if (yt.ge.edges(2,l)) then
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            if ((kb+1).eq.nvp) yt = yt - any
            rbufl(iy,j,l) = yt
            do 200 i = 1, idimp
            sbufr(i,jsr(1,l),l) = rbufl(i,j,l)
  200       continue
         else
            jss(2,l) = 2*npmax
            go to 240
         endif
c particles going down, should not happen
      elseif (yt.lt.edges(1,l)) then
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            if (kb.eq.0) yt = yt + any
            rbufl(iy,j,l) = yt
            do 210 i = 1, idimp
            sbufl(i,jsl(1,l),l) = rbufl(i,j,l)
  210       continue
         else
            jss(2,l) = 2*npmax
            go to 240
         endif
c particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 220 i = 1, idimp
         rbufl(i,jss(2,l),l) = rbufl(i,j,l)
  220    continue
      endif
  230 continue
  240 jsl(2,l) = jss(2,l)
  250 continue
c check if move would overflow particle array
  260 nps = 0
      npt = npmax
      do 270 l = 1, nblok
      jss(2,l) = npp(l) + jsl(2,l) + jsr(2,l) - jss(1,l)  
      nps = max0(nps,jss(2,l))
      npt = min0(npt,jss(2,l))
  270 continue
      ibflg(1) = nps
      ibflg(4) = -npt
c     call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
      do 360 l = 1, nblok
c distribute incoming particles from buffers
c distribute particles coming from below into holes
      jss(2,l) = min0(jss(1,l),jsl(2,l))
      do 290 j = 1, jss(2,l)
      do 280 i = 1, idimp
      part(i,ihole(j,l),l) = rbufl(i,j,l)
  280 continue
  290 continue
      if (jss(1,l).gt.jsl(2,l)) then
         jss(2,l) = min0(jss(1,l)-jsl(2,l),jsr(2,l))
      else
         jss(2,l) = jsl(2,l) - jss(1,l)
      endif
      do 320 j = 1, jss(2,l)
c no more particles coming from below
c distribute particles coming from above into holes
      if (jss(1,l).gt.jsl(2,l)) then
         do 300 i = 1, idimp
         part(i,ihole(j+jsl(2,l),l),l) = rbufr(i,j,l)
  300    continue
      else
c no more holes
c distribute remaining particles from below into bottom
         do 310 i = 1, idimp
         part(i,j+npp(l),l) = rbufl(i,j+jss(1,l),l)
  310    continue
      endif
  320 continue
      if (jss(1,l).le.jsl(2,l)) then
         npp(l) = npp(l) + (jsl(2,l) - jss(1,l))
         jss(1,l) = jsl(2,l)
      endif
      jss(2,l) = jss(1,l) - (jsl(2,l) + jsr(2,l))
      if (jss(2,l).gt.0) then
         jss(1,l) = (jsl(2,l) + jsr(2,l))
         jsr(2,l) = jss(2,l)
      else
         jss(1,l) = jss(1,l) - jsl(2,l)
         jsr(2,l) = -jss(2,l)
      endif
      do 350 j = 1, jsr(2,l)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (jss(2,l).gt.0) then
         j1 = npp(l) - j + 1
         j2 = jss(1,l) + jss(2,l) - j + 1
         if (j1.gt.ihole(j2,l)) then
c move particle only if it is below current hole
            do 330 i = 1, idimp
            part(i,ihole(j2,l),l) = part(i,j1,l)
  330       continue
         endif
      else
c no more holes
c distribute remaining particles from above into bottom
         do 340 i = 1, idimp
         part(i,j+npp(l),l) = rbufr(i,j+jss(1,l),l)
  340    continue
      endif
  350 continue
      if (jss(2,l).gt.0) then
         npp(l) = npp(l) - jsr(2,l)
      else
         npp(l) = npp(l) + jsr(2,l)
      endif
      jss(1,l) = 0
  360 continue
c check if any particles have to be passed further
c     info(5) = max0(info(5),mter)
      if (mter.eq.info(5)) then
         if (ibflg(2).gt.0) then
            write (2,*) 'Exceeded maximum number of passes = ', mter
            info(1) = -1
            return
         endif
c        write (2,*) 'Info: particles being passed further = ', ibflg(2)
      else
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 70
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(4) = nter
         write (2,*) 'Buffer overflow = ', nter
         info(1) = -2
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PFMOVE2(f,g,noff,nyp,noffs,nyps,noffd,nypd,jsr,jsl,isig
     1n,kyp,kstrt,nvp,nxv,nypmx,nblok,idps,mter,ierr)
c this subroutine moves fields into appropriate spatial regions,
c between non-uniform and uniform partitions
c f(j,k,l) = real data for grid j,k in field partition l.
c the grid is non-uniform and includes extra guard cells.
c g(j,k,l) = scratch data for grid j,k in field partition l.
c noff(l) = lowermost global gridpoint in field partition l
c nyp(l) = number of primary gridpoints in field partition l
c noffs(l)/nyps(l) = source or scratch arrays for field partition l
c noffd(l)/nypd(l) = destination or scratch arrays for field partition l
c jsl(idps,l) = number of field elements going down in field partition l
c jsr(idps,l) = number of field elements going up in field partition l
c isign = -1, move from non-uniform (noff/nyp) to uniform (kyp) fields
c isign = 1, move from uniform (kyp) to non-uniform (noff/nyp) fields
c if isign = 0, the noffs/nyps contains the source partition, noffd/nypd
c    contains the destination partition, and noff/nyp, kyp are not used.
c    the source partition noffs/nyps is modified.
c kyp = number of complex grids in each uniform field partition.
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of field partition, must be >= kyp+1
c nblok = number of field partitions.
c idps = number of partition boundaries
c mter = number of shifts required
c if mter = 0, then number of shifts is determined and returned
c ierr = (0,1) = (no,yes) error condition exists
      implicit none
      real f, g
      integer noff, nyp, noffs, nyps, noffd, nypd, jsr, jsl
      integer isign, kyp, kstrt, nvp, nxv, nypmx, nblok, idps, mter
      integer ierr
      dimension f(nxv,nypmx,nblok), g(nxv,nypmx,nblok)
      dimension noff(nblok), nyp(nblok)
      dimension noffs(nblok), nyps(nblok), noffd(nblok), nypd(nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer j, k, l
      integer nbsize, ks, iter, npr, nps, nter, koff, kl, kr, kk, nypmn
      integer msid, istatus
      integer ibflg, iwork
      dimension istatus(lstat)
      dimension ibflg(2), iwork(2)
c exit if certain flags are set
      if (mter.lt.0) return
      ks = kstrt - 2
      nbsize = nxv*nypmx
      nypmn = nypmx
      iter = 2
      ierr = 0
c move from non-uniform to uniform fields
      if (isign.lt.0) then
c copy non-uniform partition parameters
         do 10 l = 1, nblok
         noffs(l) = noff(l)
         nyps(l) = nyp(l)
         koff = kyp*(l + ks)
         noffd(l) = koff
         nypd(l) = kyp
c extend partition to include ny+1 grid
         if ((l+ks).eq.(nvp-1)) then
            nyps(l) = nyps(l) + 1
            nypd(l) = nypd(l) + 1
         endif
   10    continue
c move from uniform to non-uniform fields
      else if (isign.gt.0) then
c set uniform partition parameters
         do 20 l = 1, nblok
         koff = kyp*(l + ks)
         noffs(l) = koff
         nyps(l) = kyp
         noffd(l) = noff(l)
         nypd(l) = nyp(l)
c extend partition to include ny+1 grid
         if ((l+ks).eq.(nvp-1)) then
            nyps(l) = nyps(l) + 1
            nypd(l) = nypd(l) + 1
         endif
   20    continue
      endif
c determine number of outgoing grids
   30 do 50 l = 1, nblok
      kl = noffd(l)
      kr = kl + nypd(l)
      jsl(1,l) = 0
      jsr(1,l) = 0
      do 40 k = 1, nyps(l)
      kk = k + noffs(l)
c fields going up
      if (kk.gt.kr) then
         jsr(1,l) = jsr(1,l) + 1
c fields going down
      else if (kk.le.kl) then
         jsl(1,l) = jsl(1,l) + 1
      endif
   40 continue
   50 continue
c copy fields
      iter = iter + 2
      npr = 0
      nter = 0
c get fields from below
      do 80 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
      jsl(2,l) = 0
      jsr(2,l) = 0
c this segment is used for shared memory computers
c     if (noffs(l).gt.noffd(l)) then     
c        jsl(2,l) = jsr(1,kl)
c        do 70 k = 1, jsl(2,l)
c        do 60 j = 1, nxv
c        g(j,k,l) = f(j,k+nyps(kl)-jsr(1,kl),kl)
c  60    continue
c  70    continue
c     endif
c this segment is used for mpi computers
c post receive from left
      if (noffs(l).gt.noffd(l)) then      
         call MPI_IRECV(g,nbsize,mreal,kl-1,iter-1,lgrp,msid,ierr)
      endif
c send fields to right
      if (jsr(1,l).gt.0) then
         call MPI_SEND(f(1,nyps(l)-jsr(1,l)+1,l),nxv*jsr(1,l),mreal,kr-1
     1,iter-1,lgrp,ierr)
      endif
c wait for fields to arrive
      if (noffs(l).gt.noffd(l)) then 
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsl(2,l) = nps/nxv
      endif
   80 continue
c adjust field
      do 150 l = 1, nblok
c adjust field size
      nyps(l) = nyps(l) - jsr(1,l)
c do not allow move to overflow field array
      jsr(1,l) = max0((nyps(l)+jsl(2,l)-nypmn),0)
      nyps(l) = nyps(l) - jsr(1,l)
      if (jsr(1,l).gt.0) then
         npr = max0(npr,jsr(1,l))
c save whatever is possible into end of g
         kk = min0(jsr(1,l),nypmn-jsl(2,l))
         do 100 k = 1, kk
         do  90 j = 1, nxv
         g(j,nypmn-kk+k,l) = f(j,nyps(l)+k,l)
   90    continue
  100    continue
      endif
c shift data which is staying, if necessary
      if ((nyps(l).gt.0).and.(jsl(2,l).gt.0)) then
         do 120 k = 1, nyps(l)
         kk = nyps(l) - k + 1
         do 110 j = 1, nxv
         f(j,kk+jsl(2,l),l) = f(j,kk,l)
  110    continue
  120    continue
      endif
c insert data coming from left
      do 140 k = 1, jsl(2,l)
      do 130 j = 1, nxv
      f(j,k,l) = g(j,k,l)
  130 continue
  140 continue
c adjust field size and offset
      nyps(l) = nyps(l) + jsl(2,l)
      noffs(l) = noffs(l) - jsl(2,l)
  150 continue
c get fields from above
      do 180 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
c this segment is used for shared memory computers
c     if ((noffs(l)+nyps(l)).lt.(noffd(l)+nypd(l))) then 
c        jsr(2,l) = jsl(1,kr)
c        do 170 k = 1, jsr(2,l)
c        do 160 j = 1, nxv
c        g(j,k,l) =  f(j,k,kr)
c 160    continue
c 170    continue
c     endif
c this segment is used for mpi computers
c post receive from right
      if ((noffs(l)+nyps(l)).lt.(noffd(l)+nypd(l))) then     
         call MPI_IRECV(g,nbsize,mreal,kr-1,iter,lgrp,msid,ierr)
      endif
c send fields to left
      if (jsl(1,l).gt.0) then
         call MPI_SEND(f(1,1,l),nxv*jsl(1,l),mreal,kl-1,iter,lgrp,ierr)
      endif
c wait for fields to arrive
      if ((noffs(l)+nyps(l)).lt.(noffd(l)+nypd(l))) then  
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsr(2,l) = nps/nxv
      endif
  180 continue
c adjust field
      do 240 l = 1, nblok
c adjust field size
      nyps(l) = nyps(l) - jsl(1,l)
      noffs(l) = noffs(l) + jsl(1,l)
c shift data which is staying, if necessary
      if ((nyps(l).gt.0).and.(jsl(1,l).gt.0)) then
         do 200 k = 1, nyps(l)
         do 190 j = 1, nxv
         f(j,k,l) = f(j,k+jsl(1,l),l)
  190    continue
  200    continue
      endif
c do not allow move to overflow field array
      jsl(1,l) = max0((nyps(l)+jsr(2,l)-nypmn),0)
      if (jsl(1,l).gt.0) then
         npr = max0(npr,jsl(1,l))
         jsr(2,l) = jsr(2,l) - jsl(1,l)
c do not process if prior error
      else if (jsr(1,l).gt.0) then
         go to 230
      endif
c insert data coming from right
      do 220 k = 1, jsr(2,l)
      do 210 j = 1, nxv
      f(j,k+nyps(l),l) = g(j,k,l)
  210 continue
  220 continue
c adjust field size and offset
      nyps(l) = nyps(l) + jsr(2,l)
c check if new partition is uniform
  230 nter = nter + abs(nyps(l)-nypd(l)) + abs(noffs(l)-noffd(l))
  240 continue
c calculate number of iterations
      nps = iter/2 - 1
      if (nps.le.mter) then
c process errors
         if (npr.ne.0) then
            ierr = npr
            write (2,*) 'local field overflow error, ierr = ', ierr
            go to 250
         endif
         if (nps.lt.mter) go to 30
         go to 250
      endif
c process errors
      ibflg(1) = npr
      ibflg(2) = nter
      call PIMAX(ibflg,iwork,2,1)
c field overflow error
      if (ibflg(1).ne.0) then
         ierr = ibflg(1)
         write (2,*) 'global field overflow error, ierr = ', ierr
         go to 250
      endif
c check if any fields have to be passed further
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: fields being passed further = ', ibflg(2)
         go to 30
      endif
      mter = nps
c restore partition to normal
  250 do 260 l = 1, nblok
      if ((l+ks).eq.(nvp-1)) then
         nyps(l) = nyps(l) - 1
         nypd(l) = nypd(l) - 1
      endif
  260 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine REPART2(edges,eg,es,et3,npp,noff,nyp,npav,nypmin,nypmax
     1,kstrt,nvp,nblok,idps,inorder)
c this subroutines finds new partitions boundaries (edges, noff, nyp)
c from old partition information (npp, nyp).
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c eg/es/et3 = scratch arrays
c npp(l) = number of particles in partition l
c noff(l) = lowermost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c npav = average number of particles per partition desired
c nypmin = minimum value of nyp in new partition
c nypmax = maximum value of nyp plus guard cells in new partition
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nblok = number of field partitions.
c idps = number of partition boundaries
c inorder = (1,2) (linear,quadratic) interpolation used
      implicit none
      real edges, eg, es, et3
      integer npp, noff, nyp
      integer npav, nypmin, nypmax, kstrt, nvp, nblok, idps, inorder
      dimension edges(idps,nblok)
      dimension eg(idps,nblok), es(idps,nblok), et3(3*idps,nblok)
      dimension npp(nblok), noff(nblok), nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, iter, nter, kl, kr, l, ierr
      real anpav, apav, anpl, anpr
      integer msid, istatus
      integer ibflg, iwork
      dimension istatus(lstat)
      dimension ibflg(2), iwork(2)
c exit if flag is set
      ks = kstrt - 2
      iter = 2
      anpav = real(npav)
c copy number of particles and grid in current partition
      do 10 l = 1, nblok
      edges(1,l) = npp(l)
      edges(2,l) = nyp(l)
      et3(1,l) = edges(1,l)
      et3(2,l) = edges(2,l)
   10 continue
c perform running sum
      call PSCAN(edges,eg,es,idps,nblok)
      do 20 l = 1, nblok
      es(1,l) = et3(1,l)
      es(2,l) = et3(2,l)
      et3(1,l) = edges(1,l)
      et3(2,l) = edges(2,l)
      et3(3,l) = et3(1,l)
      et3(4,l) = et3(2,l)
      et3(5,l) = es(1,l)
      et3(6,l) = es(2,l)
      edges(2,l) = 1.0
   20 continue
c move partitions
   30 iter = iter + 2
c get partition from left
      do 40 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
c apav = desired number of particles on processor to left
      apav = real(kl)*anpav
c anpl = deficit of particles on processor to left
      anpl = apav - et3(1,l) + es(1,l)
c anpr = excess of particles on current processor
      anpr = et3(1,l) - apav - anpav
c this segment is used for shared memory computers
c     if (anpl.lt.0.) then
c        eg(1,l) = es(1,kl)
c        eg(2,l) = es(2,kl)
c     endif
c this segment is used for mpi computers
c post receive from left
      if (anpl.lt.0.) then
         call MPI_IRECV(eg,idps,mreal,kl-1,iter-1,lgrp,msid,ierr)
      endif
c send partition to right
      if (anpr.gt.0.) then
         call MPI_SEND(es,idps,mreal,kr-1,iter-1,lgrp,ierr)
      endif
c wait for partition to arrive
      if (anpl.lt.0.) call MPI_WAIT(msid,istatus,ierr)
   40 continue
c find new partitions
      nter = 0
      do 50 l = 1, nblok
      kl = l + ks
      apav = real(kl)*anpav
      anpl = apav - et3(1,l) + es(1,l)
      anpr = et3(1,l) - apav - anpav
c left boundary is on the left
      if (anpl.lt.0.) then
         if ((anpl+eg(1,l)).ge.0.) then
            edges(1,l) = (et3(2,l) - es(2,l)) + anpl*eg(2,l)/eg(1,l)
c left boundary is even further to left
         else
            nter = nter + 1
         endif
c left boundary is inside
      else if (et3(1,l).ge.apav) then
         edges(1,l) = (et3(2,l) - es(2,l)) + anpl*es(2,l)/es(1,l)
      endif
c right going data will need to be sent
      if (anpr.gt.es(1,l)) nter = nter + 1
      if (kl.gt.0) then
         et3(1,l) = et3(1,l) - es(1,l)
         et3(2,l) = et3(2,l) - es(2,l)
         es(1,l) = eg(1,l)
         es(2,l) = eg(2,l)
      endif
   50 continue
c get more data from left
      if (nter.gt.0) go to 30
      iter = nvp + 2
c restore partition data
      do 60 l = 1, nblok
      et3(1,l) = et3(3,l)
      et3(2,l) = et3(4,l)
      es(1,l) = et3(5,l)
      es(2,l) = et3(6,l)
   60 continue
c continue moving partitions
   70 iter = iter + 2
c get partition from right
      do 80 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
      apav = real(kl)*anpav
      anpl = apav - et3(1,l) + es(1,l)
c this segment is used for shared memory computers
c     if (et3(1,l).lt.apav) then
c        eg(1,l) = es(1,kr)
c        eg(2,l) = es(2,kr)
c     endif
c this segment is used for mpi computers
c post receive from right
      if (et3(1,l).lt.apav) then
         call MPI_IRECV(eg,idps,mreal,kr-1,iter,lgrp,msid,ierr)
      endif
c send partition to left
      if (anpl.gt.anpav) then
         call MPI_SEND(es,idps,mreal,kl-1,iter,lgrp,ierr)
      endif
c wait for partition to arrive
      if (et3(1,l).lt.apav) call MPI_WAIT(msid,istatus,ierr)
   80 continue
c find new partitions
      nter = 0
      do 90 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
      apav = real(kl)*anpav
      anpl = apav - et3(1,l) + es(1,l)
      anpr = et3(1,l) - apav - anpav
c left boundary is on the right
      if (et3(1,l).lt.apav) then
         if ((et3(1,l)+eg(1,l)).ge.apav) then
            edges(1,l) = et3(2,l) - (anpr + anpav)*eg(2,l)/eg(1,l)
c left boundary is even further to right
         else
            nter = nter + 1
         endif
      endif
c left going data will need to be sent
      if ((anpl-es(1,l)).gt.anpav) nter = nter + 1
      if (kr.le.nvp) then
         et3(1,l) = et3(1,l) + eg(1,l)
         et3(2,l) = et3(2,l) + eg(2,l)
         es(1,l) = eg(1,l)
         es(2,l) = eg(2,l)
      endif
   90 continue
c get more data from right
      if (nter.gt.0) go to 70
c send left edge to processor on right
      iter = 2
      do 100 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
c this segment is used for shared memory computers
c     if (kr.le.nvp) then
c        edges(2,l) = edges(1,kr)
c     else
c        edges(2,l) = et3(4,l)
c     endif
c this segment is used for mpi computers
c post receive from right
      if (kr.le.nvp) then
         call MPI_IRECV(edges(2,l),1,mreal,kr-1,iter,lgrp,msid,ierr)
      endif
c send left edge to left
      if (kl.gt.0) then
         call MPI_SEND(edges(1,l),1,mreal,kl-1,iter,lgrp,ierr)
      endif
c wait for edge to arrive
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         edges(2,l) = et3(4,l)
      endif
  100 continue
c calculate number of grids and offsets in new partitions
      do 110 l = 1, nblok
      kl = edges(1,l) + .5
      noff(l) = kl
      kr = edges(2,l) + .5
      nyp(l) = kr - kl
      if (inorder.eq.1) then
         edges(1,l) = real(kl)
         edges(2,l) = real(kr)
      endif
  110 continue
c find minimum and maximum partition size
      nypmin = nyp(1)
      nypmax = nyp(1)
      do 120 l = 1, nblok
      nypmin = min0(nypmin,nyp(l))
      nypmax = max0(nypmax,nyp(l))
  120 continue
      ibflg(1) = -nypmin
      ibflg(2) = nypmax
      call PIMAX(ibflg,iwork,2,1)
      nypmin = -ibflg(1)
      nypmax = ibflg(2) + 1
      if (inorder.eq.2) nypmax = nypmax + 2
      return
      end
c-----------------------------------------------------------------------
      subroutine REPARTD2(edges,edg,eds,eg,es,et2,npic,noff,nyp,anpav,ny
     1pmin,nypmax,kstrt,nvp,nblok,idps,nypm)
c this subroutines finds new partitions boundaries (edges,noff,nyp)
c from old partition information (npic,nyp).
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c edg/eds/eg/es/et2 = scratch arrays
c npic(l) = number of particles per grid in partition l
c noff(l) = lowermost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c anpav = average number of particles per partition desired
c nypmin/nypmax = minimum/maximum value of nyp in new partition
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nblok = number of field partitions.
c idps = number of partition boundaries
c nypm = maximum size of particle partition
      implicit none
      real edges, edg, eds, eg, es, et2
      integer npic, noff, nyp
      real anpav
      integer nypmin, nypmax, kstrt, nvp, nblok, idps, nypm
      dimension edges(idps,nblok)
      dimension edg(nypm,nblok), eds(nypm,nblok)
      dimension eg(idps,nblok), es(idps,nblok), et2(2*idps,nblok)
      dimension npic(nypm,nblok), noff(nblok), nyp(nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, iter, nter, nyp1, npav, k1, kl, kr, k, l, ierr
      real sum1, at1, at2, apav, anpl, anpr
      integer msid, istatus
      integer ibflg, iwork
      dimension istatus(lstat)
      dimension ibflg(2), iwork(2)
c exit if flag is set
      ks = kstrt - 2
      iter = 2
c copy number of particles and grid in current partition
      do 20 l = 1, nblok
      sum1 = 0.
      do 10 k = 1, nyp(l)
      at1 = npic(k,l)
      sum1 = sum1 + at1
      eds(k,l) = at1
   10 continue
      edges(1,l) = sum1
      edges(2,l) = nyp(l)
      et2(1,l) = edges(1,l)
      et2(2,l) = edges(2,l)
   20 continue
c perform running sum
      call PSCAN(edges,eg,es,idps,nblok)
      do 30 l = 1, nblok
      es(1,l) = et2(1,l)
      es(2,l) = et2(2,l)
      et2(1,l) = edges(1,l)
      et2(2,l) = edges(2,l)
      et2(3,l) = et2(1,l)
      et2(4,l) = et2(2,l)
      eg(1,l) = 0.
      eg(2,l) = 0.
      edges(2,l) = 1.0
   30 continue
c move partitions
   40 iter = iter + 2
c get partition from left
      do 60 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
c apav = desired number of particles on processor to left
      npav = real(kl)*anpav + .5
      apav = real(npav)
c anpl = deficit of particles on processor to left
      anpl = apav - et2(1,l) + es(1,l)
c anpr = excess of particles on current processor
      npav = real(kl+1)*anpav + .5
      anpr = et2(1,l) - real(npav)
c this segment is used for shared memory computers
c     if (anpl.lt.0.) then
c        nyp1 = es(2,kl)
c        do 50 k = 1, nyp1
c        edg(k,l) = eds(k,kl)
c  50    continue
c        eg(1,l) = es(1,kl)
c        eg(2,l) = es(2,kl)
c     endif
c this segment is used for mpi computers
c post receive from left
      if (anpl.lt.0.) then
         call MPI_IRECV(edg,nypm,mreal,kl-1,iter-1,lgrp,msid,ierr)
      endif
c send partition to right
      if (anpr.gt.0.) then
         nyp1 = es(2,l)
         call MPI_SEND(eds,nyp1,mreal,kr-1,iter-1,lgrp,ierr)
      endif
c wait for partition to arrive
      if (anpl.lt.0.) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nyp1,ierr)
         eg(2,l) = nyp1
         sum1 = 0.
         do 50 k = 1, nyp1
         sum1 = sum1 + edg(k,l)
   50    continue
         eg(1,l) = sum1
      endif
   60 continue
c find new partitions
      nter = 0
      do 100 l = 1, nblok
      kl = l + ks
      npav = real(kl)*anpav + .5
      apav = real(npav)
      anpl = apav - et2(1,l) + es(1,l)
      npav = real(kl+1)*anpav + .5
      anpr = et2(1,l) - real(npav)
c left boundary is on the left
      if (anpl.lt.0.) then
         if ((anpl+eg(1,l)).ge.0.) then
            nyp1 = eg(2,l)
            k1 = nyp1
            sum1 = 0.
   70       at1 = sum1
            sum1 = sum1 - edg(k1,l)
            k1 = k1 - 1
            if ((sum1.gt.anpl).and.(k1.gt.0)) go to 70
            at1 = real(nyp1 - k1 - 1) + (anpl - at1)/(sum1 - at1)
            edges(1,l) = (et2(2,l) - es(2,l)) - at1
c left boundary is even further to left
         else
            nter = nter + 1
         endif
c left boundary is inside
      else if (et2(1,l).ge.apav) then
         nyp1 = es(2,l)
         k1 = 1
         sum1 = 0.
   80    at1 = sum1
         sum1 = sum1 + eds(k1,l)
         k1 = k1 + 1
         if ((sum1.lt.anpl).and.(k1.le.nyp1)) go to 80
         at2 = real(k1 - 2)
         if (sum1.gt.at1) at2 = at2 + (anpl - at1)/(sum1 - at1)
         edges(1,l) = (et2(2,l) - es(2,l)) + at2
      endif
c right going data will need to be sent
      if (anpr.gt.es(1,l)) nter = nter + 1
      if (kl.gt.0) then
         nyp1 = eg(2,l)
         do 90 k = 1, nyp1
         eds(k,l) = edg(k,l)
   90    continue
         et2(1,l) = et2(1,l) - es(1,l)
         et2(2,l) = et2(2,l) - es(2,l)
         es(1,l) = eg(1,l)
         es(2,l) = eg(2,l)
      endif
  100 continue
c get more data from left
      if (nter.gt.0) go to 40
      iter = nvp + 2
c restore partition data
      do 120 l = 1, nblok
      sum1 = 0.
      do 110 k = 1, nyp(l)
      at1 = npic(k,l)
      sum1 = sum1 + at1
      eds(k,l) = at1
  110 continue
      et2(1,l) = et2(3,l)
      et2(2,l) = et2(4,l)
      es(1,l) = sum1
      es(2,l) = nyp(l)
      eg(1,l) = 0.
      eg(2,l) = 0.
  120 continue
c continue moving partitions
  130 iter = iter + 2
c get partition from right
      do 150 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
      npav = real(kl)*anpav + .5
      apav = real(npav)
      npav = real(kl-1)*anpav + .5
      anpl = real(npav) - et2(1,l) + es(1,l)
c this segment is used for shared memory computers
c     if (et2(1,l).lt.apav) then
c        nyp1 = es(2,kr)
c        do 140 k = 1, nyp1
c        edg(k,l) = eds(k,kr)
c 140    continue
c        eg(1,l) = es(1,kr)
c        eg(2,l) = es(2,kr)
c     endif
c this segment is used for mpi computers
c post receive from right
      if (et2(1,l).lt.apav) then
         call MPI_IRECV(edg,nypm,mreal,kr-1,iter,lgrp,msid,ierr)
      endif
c send partition to left
      if (anpl.gt.0.) then
         nyp1 = es(2,l)
         call MPI_SEND(eds,nyp1,mreal,kl-1,iter,lgrp,ierr)
      endif
c wait for partition to arrive
      if (et2(1,l).lt.apav) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nyp1,ierr)
         eg(2,l) = nyp1
         sum1 = 0.
         do 140 k = 1, nyp1
         sum1 = sum1 + edg(k,l)
  140    continue
         eg(1,l) = sum1
      endif
  150 continue
c find new partitions
      nter = 0
      do 180 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
      npav = real(kl)*anpav + .5
      apav = real(npav)
      npav = real(kl-1)*anpav + .5
      anpl = real(npav) - et2(1,l) + es(1,l)
      anpr = et2(1,l) - apav
c left boundary is on the right
      if (et2(1,l).lt.apav) then
         if ((et2(1,l)+eg(1,l)).ge.apav) then
            nyp1 = eg(2,l)
            k1 = 1
            sum1 = 0.
            at2 = -anpr
  160       at1 = sum1
            sum1 = sum1 + edg(k1,l)
            k1 = k1 + 1
            if ((sum1.lt.at2).and.(k1.le.nyp1)) go to 160
            at1 = real(k1 - 2) + (at2 - at1)/(sum1 - at1)
            edges(1,l) = et2(2,l) + at1
c left boundary is even further to right
         else
            nter = nter + 1
         endif
      endif
c left going data will need to be sent
      if (anpl.gt.es(1,l)) nter = nter + 1
      if (kr.le.nvp) then
         nyp1 = eg(2,l)
         do 170 k = 1, nyp1
         eds(k,l) = edg(k,l)
  170    continue
         et2(1,l) = et2(1,l) + eg(1,l)
         et2(2,l) = et2(2,l) + eg(2,l)
         es(1,l) = eg(1,l)
         es(2,l) = eg(2,l)
      endif
  180 continue
c get more data from right
      if (nter.gt.0) go to 130
c send left edge to processor on right
      iter = 2
      do 190 l = 1, nblok
      kr = l + ks + 2
      kl = l + ks
c this segment is used for shared memory computers
c     if (kr.le.nvp) then
c        edges(2,l) = edges(1,kr)
c     else
c        edges(2,l) = et2(4,l)
c     endif
c this segment is used for mpi computers
c post receive from right
      if (kr.le.nvp) then
         call MPI_IRECV(edges(2,l),1,mreal,kr-1,iter,lgrp,msid,ierr)
      endif
c send left edge to left
      if (kl.gt.0) then
         call MPI_SEND(edges(1,l),1,mreal,kl-1,iter,lgrp,ierr)
      endif
c wait for edge to arrive
      if (kr.le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         edges(2,l) = et2(4,l)
      endif
  190 continue
c calculate number of grids and offsets in new partitions
      do 200 l = 1, nblok
      kl = edges(1,l) + .5
      noff(l) = kl
      kr = edges(2,l) + .5
      nyp(l) = kr - kl
      edges(1,l) = real(kl)
      edges(2,l) = real(kr)
  200 continue
c find minimum and maximum partition size
      nypmin = nyp(1)
      nypmax = nyp(1)
      do 210 l = 1, nblok
      nypmin = min0(nypmin,nyp(l))
      nypmax = max0(nypmax,nyp(l))
  210 continue
      ibflg(1) = -nypmin
      ibflg(2) = nypmax
      call PIMAX(ibflg,iwork,2,1)
      nypmin = -ibflg(1)
      nypmax = ibflg(2)
      return
      end
c-----------------------------------------------------------------------
      subroutine FNOFF2(edges,noff,nyp,nypmin,nypmax,nblok,idps)
c this subroutines finds new partitions arrays (noff,nyp) from edges
c edges(1,l) = lower boundary of particle partition l
c edges(2,l) = upper boundary of particle partition l
c noff(l) = lowermost global gridpoint in particle partition l
c nyp(l) = number of primary gridpoints in particle partition l
c nypmin/nypmax = minimum/maximum value of nyp in new partition
c nblok = number of field partitions.
c idps = number of partition boundaries
      implicit none
      real edges
      integer noff, nyp, nypmin, nypmax, nblok, idps
      dimension edges(idps,nblok)
      dimension noff(nblok), nyp(nblok)
c local data
      integer kl, kr, l
      integer ibflg, iwork
      dimension ibflg(2), iwork(2)
c calculate number of grids and offsets in new partitions
      do 10 l = 1, nblok
      kl = edges(1,l) + .5
      noff(l) = kl
      kr = edges(2,l) + .5
      nyp(l) = kr - kl
   10 continue
c find minimum and maximum partition size
      nypmin = nyp(1)
      nypmax = nyp(1)
      do 20 l = 1, nblok
      nypmin = min0(nypmin,nyp(l))
      nypmax = max0(nypmax,nyp(l))
   20 continue
      ibflg(1) = -nypmin
      ibflg(2) = nypmax
      call PIMAX(ibflg,iwork,2,1)
      nypmin = -ibflg(1)
      nypmax = ibflg(2)
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,jb
     1lok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxpd = second dimension of f/g
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      complex f, g, s, t
      dimension f(nxv,kypd,kblok), g(nyv,kxpd,jblok)
      dimension s(kxp,kyp,kblok), t(kxp,kyp,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 20 k = 1, kyp
c     do 10 j = 1, kxp
c     g(k+koff,j,l) = f(j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         call MPI_IRECV(t(1,1,l),kxp*kyp,mcplx,ir-1,ir+kxym+1,lgrp,msid,
     1ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp
         do 10 j = 1, kxp
         s(j,k,l) = f(j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,l),kxp*kyp,mcplx,is-1,l+ks+kxym+2,lgrp,ierr
     1)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kyp
         do 30 j = 1, kxp
         g(k+koff,j,l) = t(j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P2TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j
     1blok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:2,k+kyp*(m-1),j,l) = f(1:2,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxpd = second dimension of f/g
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      complex f, g, s, t
      dimension f(2,nxv,kypd,kblok), g(2,nyv,kxpd,jblok)
      dimension s(2,kxp,kyp,kblok), t(2,kxp,kyp,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 20 k = 1, kyp
c     do 10 j = 1, kxp
c     g(1,k+koff,j,l) = f(1,j+joff,k,i)
c     g(2,k+koff,j,l) = f(2,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         call MPI_IRECV(t(1,1,1,l),2*kxp*kyp,mcplx,ir-1,ir+kxym+1,lgrp,m
     1sid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp
         do 10 j = 1, kxp
         s(1,j,k,l) = f(1,j+joff,k,l)
         s(2,j,k,l) = f(2,j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,1,l),2*kxp*kyp,mcplx,is-1,l+ks+kxym+2,lgrp,
     1ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kyp
         do 30 j = 1, kxp
         g(1,k+koff,j,l) = t(1,j,k,l)
         g(2,k+koff,j,l) = t(2,j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j
     1blok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:3,k+kyp*(m-1),j,l) = f(1:3,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxpd = second dimension of f/g
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      complex f, g, s, t
      dimension f(3,nxv,kypd,kblok), g(3,nyv,kxpd,jblok)
      dimension s(3,kxp,kyp,kblok), t(3,kxp,kyp,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 20 k = 1, kyp
c     do 10 j = 1, kxp
c     g(1,k+koff,j,l) = f(1,j+joff,k,i)
c     g(2,k+koff,j,l) = f(2,j+joff,k,i)
c     g(3,k+koff,j,l) = f(3,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         call MPI_IRECV(t(1,1,1,l),3*kxp*kyp,mcplx,ir-1,ir+kxym+1,lgrp,m
     1sid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp
         do 10 j = 1, kxp
         s(1,j,k,l) = f(1,j+joff,k,l)
         s(2,j,k,l) = f(2,j+joff,k,l)
         s(3,j,k,l) = f(3,j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,1,l),3*kxp*kyp,mcplx,is-1,l+ks+kxym+2,lgrp,
     1ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kyp
         do 30 j = 1, kxp
         g(1,k+koff,j,l) = t(1,j,k,l)
         g(2,k+koff,j,l) = t(2,j,k,l)
         g(3,k+koff,j,l) = t(3,j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOSEX(f,g,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,jblok
     1,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages.
c f = complex input array
c g = complex output array
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kxp/kyp = number of data values per block in x/y
c kypd/kxpd = second dimension of f/g
c jblok/kblok = number of data blocks in x/y
c optimized version
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp
      integer kxpd, kypd, jblok, kblok
      complex f, g
      dimension f(nxv*kypd*kblok), g(nyv*kxpd*jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, l, i, joff, koff, k, j
      integer jkblok, kxym, mtr, ntr, mntr, msid
      integer ir0, is0, ii, ir, is, ioff, ierr, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 20 k = 1, kyp
c     do 10 j = 1, kxp
c     g(k+koff+nyv*(j-1+kxpd*(l-1))) = f(j+joff+nxv*(k-1+kypd*(i-1)))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
c transpose local data
      do 50 l = 1, jkblok
      ioff = kxb*(l - 1) - 1
      koff = kypd*(l - 1) - 1
      do 40 i = 1, kxym
      is0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      do 30 ii = 1, ntr
      if (kstrt.le.ny) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         is = kyp*(is + ioff) - 1
         do 20 k = 1, kyp
         do 10 j = 1, kxp
         g(j+kxp*(k+is)) = f(j+joff+nxv*(k+koff))
   10    continue
   20    continue
      endif
   30 continue
   40 continue
   50 continue
c exchange data
      do 80 l = 1, jkblok
      ioff = kxb*(l - 1) - 1
      koff = kyb*(l - 1) - 1
      do 70 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 60 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         call MPI_IRECV(f(1+kxp*kyp*(ir+koff)),kxp*kyp,mcplx,ir-1,ir+kxy
     1m+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         call MPI_SEND(g(1+kxp*kyp*(is+ioff)),kxp*kyp,mcplx,is-1,l+ks+kx
     1ym+2,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   60 continue
   70 continue
   80 continue
c transpose local data
      do 130 l = 1, jkblok
      ioff = kyb*(l - 1) - 1
      joff = kxpd*(l - 1) - 1
      do 120 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      do 110 ii = 1, mtr
      if (kstrt.le.nx) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*(ir - 1)
         ir = kyp*(ir + ioff) - 1
         do 100 k = 1, kyp
         do 90 j = 1, kxp
         g(k+koff+nyv*(j+joff)) = f(j+kxp*(k+ir))
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P2TPOSEX(f,g,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,jblo
     1k,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:2,k+kyp*(m-1),j,l) = f(1:2,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages.
c f = complex input array
c g = complex output array
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kxp/kyp = number of data values per block in x/y
c kypd/kxpd = second dimension of f/g
c jblok/kblok = number of data blocks in x/y
c optimized version
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp
      integer kxpd, kypd, jblok, kblok
      complex f, g
      dimension f(2*nxv*kypd*kblok), g(2*nyv*kxpd*jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, l, i, joff, koff, k, j
      integer jkblok, kxym, mtr, ntr, mntr, msid
      integer ir0, is0, ii, ir, is, ioff, ierr, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks) - 1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1) - 1
c     do 20 k = 1, kyp
c     do 10 j = 1, kxp
c     g(1+2*(k+koff+nyv*(j-1+kxpd*(l-1)))) = f(1+2*(j+joff+nxv*(k-1+kypd
c    1*(i-1))))
c     g(2+2*(k+koff+nyv*(j-1+kxpd*(l-1)))) = f(2+2*(j+joff+nxv*(k-1+kypd
c    1*(i-1))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
c transpose local data
      do 50 l = 1, jkblok
      ioff = kxb*(l - 1) - 1
      koff = kypd*(l - 1) - 1
      do 40 i = 1, kxym
      is0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      do 30 ii = 1, ntr
      if (kstrt.le.ny) then
         is = is0 + kxym*(ii - 1)
         joff = 2*kxp*(is - 1)
         is = kyp*(is + ioff) - 1
         do 20 k = 1, kyp
         do 10 j = 1, 2*kxp
         g(j+2*kxp*(k+is)) = f(j+joff+2*nxv*(k+koff))
   10    continue
   20    continue
      endif
   30 continue
   40 continue
   50 continue
      do 80 l = 1, jkblok
      ioff = kxb*(l - 1) - 1
      koff = kyb*(l - 1) - 1
      do 70 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 60 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         call MPI_IRECV(f(1+2*kxp*kyp*(ir+koff)),2*kxp*kyp,mcplx,ir-1,ir
     1+kxym+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         call MPI_SEND(g(1+2*kxp*kyp*(is+ioff)),2*kxp*kyp,mcplx,is-1,l+k
     1s+kxym+2,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   60 continue
   70 continue
   80 continue
c transpose local data
      do 130 l = 1, jkblok
      ioff = kyb*(l - 1) - 1
      joff = kxpd*(l - 1) - 1
      do 120 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      do 110 ii = 1, mtr
      if (kstrt.le.nx) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*(ir - 1)
         ir = kyp*(ir + ioff) - 1
         do 100 k = 1, kyp
         do 90 j = 1, kxp
         g(2*(k+koff+nyv*(j+joff))-1) = f(2*(j+kxp*(k+ir))-1)
         g(2*(k+koff+nyv*(j+joff))) = f(2*(j+kxp*(k+ir)))
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOSEX(f,g,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,jblo
     1k,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:3,k+kyp*(m-1),j,l) = f(1:3,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages.
c f = complex input array
c g = complex output array
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kxp/kyp = number of data values per block in x/y
c kypd/kxpd = second dimension of f/g
c jblok/kblok = number of data blocks in x/y
c optimized version
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp
      integer kxpd, kypd, jblok, kblok
      complex f, g
      dimension f(3*nxv*kypd*kblok), g(3*nyv*kxpd*jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, l, i, joff, koff, k, j
      integer jkblok, kxym, mtr, ntr, mntr, msid
      integer ir0, is0, ii, ir, is, ioff, ierr, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks) - 1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1) - 1
c     do 20 k = 1, kyp
c     do 10 j = 1, kxp
c     g(1+3*(k+koff+nyv*(j-1+kxpd*(l-1)))) = f(1+3*(j+joff+nxv*(k-1+kypd
c    1*(i-1))))
c     g(2+3*(k+koff+nyv*(j-1+kxpd*(l-1)))) = f(2+3*(j+joff+nxv*(k-1+kypd
c    1*(i-1))))
c     g(3+3*(k+koff+nyv*(j-1+kxpd*(l-1)))) = f(3+3*(j+joff+nxv*(k-1+kypd
c    1*(i-1))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
c transpose local data
      do 50 l = 1, jkblok
      ioff = kxb*(l - 1) - 1
      koff = kypd*(l - 1) - 1
      do 40 i = 1, kxym
      is0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      do 30 ii = 1, ntr
      if (kstrt.le.ny) then
         is = is0 + kxym*(ii - 1)
         joff = 3*kxp*(is - 1)
         is = kyp*(is + ioff) - 1
         do 20 k = 1, kyp
         do 10 j = 1, 3*kxp
         g(j+3*kxp*(k+is)) = f(j+joff+3*nxv*(k+koff))
   10    continue
   20    continue
      endif
   30 continue
   40 continue
   50 continue
      do 80 l = 1, jkblok
      ioff = kxb*(l - 1) - 1
      koff = kyb*(l - 1) - 1
      do 70 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 60 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         call MPI_IRECV(f(1+3*kxp*kyp*(ir+koff)),3*kxp*kyp,mcplx,ir-1,ir
     1+kxym+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         call MPI_SEND(g(1+3*kxp*kyp*(is+ioff)),3*kxp*kyp,mcplx,is-1,l+k
     1s+kxym+2,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   60 continue
   70 continue
   80 continue
c transpose local data
      do 130 l = 1, jkblok
      ioff = kyb*(l - 1) - 1
      joff = kxpd*(l - 1) - 1
      do 120 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      do 110 ii = 1, mtr
      if (kstrt.le.nx) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*(ir - 1)
         ir = kyp*(ir + ioff) - 1
         do 100 k = 1, kyp
         do 90 j = 1, kxp
         g(3*(k+koff+nyv*(j+joff))-2) = f(3*(j+kxp*(k+ir))-2)
         g(3*(k+koff+nyv*(j+joff))-1) = f(3*(j+kxp*(k+ir))-1)
         g(3*(k+koff+nyv*(j+joff))) = f(3*(j+kxp*(k+ir)))
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNTPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j
     1blok,kblok,ndim)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxpd = second dimension of f/g
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
c ndim = leading dimension of arrays f and g
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok, ndim
      complex f, g, s, t
      dimension f(ndim,nxv,kypd,kblok), g(ndim,nyv,kxpd,jblok)
      dimension s(ndim,kxp,kyp,kblok), t(ndim,kxp,kyp,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j, n
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 50 l = 1, jblok
c     joff = kxp*(l + ks)
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 30 k = 1, kyp
c     do 20 j = 1, kxp
c     do 10 n = 1, ndim
c     g(n,k+koff,j,l) = f(n,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 90 l = 1, jkblok
      do 80 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 70 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         call MPI_IRECV(t(1,1,1,l),ndim*kxp*kyp,mcplx,ir-1,ir+kxym+1,lgr
     1p,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 30 k = 1, kyp
         do 20 j = 1, kxp
         do 10 n = 1, ndim
         s(n,j,k,l) = f(n,j+joff,k,l)
   10    continue
   20    continue
   30    continue
         call MPI_SEND(s(1,1,1,l),ndim*kxp*kyp,mcplx,is-1,l+ks+kxym+2,lg
     1rp,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 60 k = 1, kyp
         do 50 j = 1, kxp
         do 40 n = 1, ndim
         g(n,k+koff,j,l) = t(n,j,k,l)
   40    continue
   50    continue
   60    continue
      endif
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNTPOSEX(f,g,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,jblo
     1k,kblok,ndim)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages, but
c waits after each pair.
c f = complex input array
c g = complex output array
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kxp/kyp = number of data values per block in x/y
c kypd/kxpd = second dimension of f/g
c jblok/kblok = number of data blocks in x/y
c ndim = leading dimension of arrays f and g
c optimized version
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp
      integer kxpd, kypd, jblok, kblok, ndim
      complex f, g
      dimension f(ndim*nxv*kypd*kblok), g(ndim*nyv*kxpd*jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, l, i, joff, koff, k, j, n
      integer jkblok, kxym, mtr, ntr, mntr, msid
      integer ir0, is0, ii, ir, is, ioff, ierr, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 50 l = 1, jblok
c     joff = kxp*(l + ks) - 1
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1) - 1
c     do 30 k = 1, kyp
c     do 20 j = 1, kxp
c     do 10 n = 1, ndim
c     g(n+ndim*(k+koff+nyv*(j-1+kxpd*(l-1)))) = f(n+ndim*(j+joff+nxv*(k-
c    11+kypd*(i-1))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
c transpose local data
      do 50 l = 1, jkblok
      ioff = kxb*(l - 1) - 1
      koff = kypd*(l - 1) - 1
      do 40 i = 1, kxym
      is0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      do 30 ii = 1, ntr
      if (kstrt.le.ny) then
         is = is0 + kxym*(ii - 1)
         joff = ndim*kxp*(is - 1)
         is = kyp*(is + ioff) - 1
         do 20 k = 1, kyp
         do 10 j = 1, ndim*kxp
         g(j+ndim*kxp*(k+is)) = f(j+joff+ndim*nxv*(k+koff))
   10    continue
   20    continue
      endif
   30 continue
   40 continue
   50 continue
      do 80 l = 1, jkblok
      ioff = kxb*(l - 1) - 1
      koff = kyb*(l - 1) - 1
      do 70 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 60 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         call MPI_IRECV(f(1+ndim*kxp*kyp*(ir+koff)),ndim*kxp*kyp,mcplx,i
     1r-1,ir+kxym+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         call MPI_SEND(g(1+ndim*kxp*kyp*(is+ioff)),ndim*kxp*kyp,mcplx,is
     1-1,l+ks+kxym+2,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   60 continue
   70 continue
   80 continue
c transpose local data
      do 140 l = 1, jkblok
      ioff = kyb*(l - 1) - 1
      joff = kxpd*(l - 1) - 1
      do 130 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      do 120 ii = 1, mtr
      if (kstrt.le.nx) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*(ir - 1) - 1
         ir = kyp*(ir + ioff) - 1
         do 110 k = 1, kyp
         do 100 j = 1, kxp
         do 90 n = 1, ndim
         g(ndim*(k+koff+nyv*(j+joff))+n) = f(ndim*(j+kxp*(k+ir)-1)+n)
   90    continue
  100    continue
  110    continue
      endif
  120 continue
  130 continue
  140 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PN2TPOSE(f1,f2,g1,g2,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kx
     1pd,kypd,jblok,kblok,ndim1,ndim2)
c this subroutine performs a transpose of two matrices f1 and f2,
c distributed in y, to two matrices g1 and g2, distributed in x,
c that is,
c g1(1:ndim1,k+kyp*(m-1),j,l) = f1(1:ndim1,j+kxp*(l-1),k,m), and
c g2(1:ndim2,k+kyp*(m-1),j,l) = f2(1:ndim2,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f1, f2 = complex input arrays
c g1, g2 = complex output arrays
c s, t = complex scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxpd = second dimension of f/g
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
c ndim1 = leading dimension of arrays f1 and g1
c ndim2 = leading dimension of arrays f2 and g2
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok, ndim1, ndim2
      complex f1, f2, g1, g2, s, t
      dimension f1(ndim1,nxv,kypd,kblok), f2(ndim2,nxv,kypd,kblok)
      dimension g1(ndim1,nyv,kxpd,jblok), g2(ndim2,nyv,kxpd,jblok)
      dimension s(ndim1+ndim2,kxp,kyp,kblok)
      dimension t(ndim1+ndim2,kxp,kyp,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, ndim
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j, n
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
      ndim = ndim1 + ndim2
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     do 90 l = 1, jblok
c     joff = kxp*(l + ks)

c     do 40 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 30 k = 1, kyp
c     do 20 j = 1, kxp
c     do 10 n = 1, ndim1
c     g1(n,k+koff,j,l) = f1(n,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c     do 80 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 70 k = 1, kyp
c     do 60 j = 1, kxp
c     do 50 n = 1, ndim2
c     g2(n,k+koff,j,l) = f2(n,j+joff,k,i)
c  50 continue
c  60 continue
c  70 continue
c  80 continue
c  90 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 130 l = 1, jkblok
      do 120 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 110 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         call MPI_IRECV(t(1,1,1,l),ndim*kxp*kyp,mcplx,ir-1,ir+kxym+1,lgr
     1p,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 30 k = 1, kyp
         do 20 j = 1, kxp
         do 10 n = 1, ndim1
         s(n,j,k,l) = f1(n,j+joff,k,l)
   10    continue
   20    continue
   30    continue
         do 60 k = 1, kyp
         do 50 j = 1, kxp
         do 40 n = 1, ndim2
         s(n+ndim1,j,k,l) = f2(n,j+joff,k,l)
   40    continue
   50    continue
   60    continue
         call MPI_SEND(s(1,1,1,l),ndim*kxp*kyp,mcplx,is-1,l+ks+kxym+2,lg
     1rp,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 100 k = 1, kyp
         do 90 j = 1, kxp
         do 70 n = 1, ndim1
         g1(n,k+koff,j,l) = t(n,j,k,l)
   70    continue
         do 80 n = 1, ndim2
         g2(n,k+koff,j,l) = t(n+ndim1,j,k,l)
   80    continue
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PRTPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j
     1blok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(k+kyp*(m-1),j,l) = f(j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c includes an extra guard cell for last row and column
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = real input array
c g = real output array
c s, t = real scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv = first dimension of f >= nx+1
c nyv = first dimension of g >= ny+1
c kypd = second dimension of f >= kyp+1
c kxpd = second dimension of g >= kxp+1
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      real f, g, s, t
      dimension f(nxv,kypd,kblok), g(nyv,kxpd,jblok)
      dimension s(kxp+1,kyp+1,kblok), t(kxp+1,kyp+1,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, kxp1, kyp1, kxpt, kypt
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c set constants to receive extra guard cells
      kxp1 = kxp + 1
      kyp1 = kyp + 1
      kxpt = kxp
      if (kstrt.eq.kxb) kxpt = kxp1
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     kypt = kyp
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     if ((l+ks).eq.(kxb-1)) kxpt = kxp1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     if (i.eq.kyb) kypt = kyp1
c     do 20 k = 1, kypt
c     do 10 j = 1, kxpt
c     g(k+koff,j,l) = f(j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         kypt = kyp
         if (ir.eq.kyb) kypt = kyp1
         call MPI_IRECV(t(1,1,l),kxp1*kyp1,mreal,ir-1,ir+kxym+1,lgrp,msi
     1d,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp1
         do 10 j = 1, kxp1
         s(j,k,l) = f(j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,l),kxp1*kyp1,mreal,is-1,l+ks+kxym+2,lgrp,ie
     1rr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kypt
         do 30 j = 1, kxpt
         g(k+koff,j,l) = t(j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PR2TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,
     1jblok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:2,k+kyp*(m-1),j,l) = f(1:2,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c includes an extra guard cell for last row and column
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = real input array
c g = real output array
c s, t = real scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv = first dimension of f >= nx+1
c nyv = first dimension of g >= ny+1
c kypd = second dimension of f >= kyp+1
c kxpd = second dimension of g >= kxp+1
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      real f, g, s, t
      dimension f(2,nxv,kypd,kblok), g(2,nyv,kxpd,jblok)
      dimension s(2,kxp+1,kyp+1,kblok), t(2,kxp+1,kyp+1,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, kxp1, kyp1, kxpt, kypt
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c set constants to receive extra guard cells
      kxp1 = kxp + 1
      kyp1 = kyp + 1
      kxpt = kxp
      if (kstrt.eq.kxb) kxpt = kxp1
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     kypt = kyp
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     if ((l+ks).eq.(kxb-1)) kxpt = kxp1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     if (i.eq.kyb) kypt = kyp1
c     do 20 k = 1, kypt
c     do 10 j = 1, kxpt
c     g(1,k+koff,j,l) = f(1,j+joff,k,i)
c     g(2,k+koff,j,l) = f(2,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         kypt = kyp
         if (ir.eq.kyb) kypt = kyp1
         call MPI_IRECV(t(1,1,1,l),2*kxp1*kyp1,mreal,ir-1,ir+kxym+1,lgrp
     1,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp1
         do 10 j = 1, kxp1
         s(1,j,k,l) = f(1,j+joff,k,l)
         s(2,j,k,l) = f(2,j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,1,l),2*kxp1*kyp1,mreal,is-1,l+ks+kxym+2,lgr
     1p,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kypt
         do 30 j = 1, kxpt
         g(1,k+koff,j,l) = t(1,j,k,l)
         g(2,k+koff,j,l) = t(2,j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PR3TPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,
     1jblok,kblok)
c this subroutine performs a transpose of a matrix f, distributed in y,
c to a matrix g, distributed in x, that is,
c g(1:3,k+kyp*(m-1),j,l) = f(1:3,j+kxp*(l-1),k,m), where
c 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kyp
c and where indices l and m can be distributed across processors.
c includes an extra guard cell for last row and column
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = real input array
c g = real output array
c s, t = real scratch arrays
c nx/ny = number of points in x/y
c kstrt = starting data block number
c nxv = first dimension of f >= nx+1
c nyv = first dimension of g >= ny+1
c kypd = second dimension of f >= kyp+1
c kxpd = second dimension of g >= kxp+1
c kxp/kyp = number of data values per block in x/y
c jblok/kblok = number of data blocks in x/y
      implicit none
      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd
      integer jblok, kblok
      real f, g, s, t
      dimension f(3,nxv,kypd,kblok), g(3,nyv,kxpd,jblok)
      dimension s(3,kxp+1,kyp+1,kblok), t(3,kxp+1,kyp+1,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ks, kxb, kyb, kxp1, kyp1, kxpt, kypt
      integer jkblok, kxym, mtr, ntr, mntr
      integer l, i, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      ks = kstrt - 2
      kxb = nx/kxp
      kyb = ny/kyp
c set constants to receive extra guard cells
      kxp1 = kxp + 1
      kyp1 = kyp + 1
      kxpt = kxp
      if (kstrt.eq.kxb) kxpt = kxp1
c this segment is used for shared memory computers
c     if (kstrt.gt.nx) return
c     kypt = kyp
c     do 40 l = 1, jblok
c     joff = kxp*(l + ks)
c     if ((l+ks).eq.(kxb-1)) kxpt = kxp1
c     do 30 i = 1, kyb
c     koff = kyp*(i - 1)
c     if (i.eq.kyb) kypt = kyp1
c     do 20 k = 1, kypt
c     do 10 j = 1, kxpt
c     g(1,k+koff,j,l) = f(1,j+joff,k,i)
c     g(2,k+koff,j,l) = f(2,j+joff,k,i)
c     g(3,k+koff,j,l) = f(3,j+joff,k,i)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 70 l = 1, jkblok
      do 60 i = 1, kxym
      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1
      is0 = ir0
      do 50 ii = 1, mntr
c post receive
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         kypt = kyp
         if (ir.eq.kyb) kypt = kyp1
         call MPI_IRECV(t(1,1,1,l),3*kxp1*kyp1,mreal,ir-1,ir+kxym+1,lgrp
     1,msid,ierr)
      endif
c send data
      if ((kstrt.le.ny).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxp*(is - 1)
         do 20 k = 1, kyp1
         do 10 j = 1, kxp1
         s(1,j,k,l) = f(1,j+joff,k,l)
         s(2,j,k,l) = f(2,j+joff,k,l)
         s(3,j,k,l) = f(3,j+joff,k,l)
   10    continue
   20    continue
         call MPI_SEND(s(1,1,1,l),3*kxp1*kyp1,mreal,is-1,l+ks+kxym+2,lgr
     1p,ierr)
      endif
c receive data
      if ((kstrt.le.nx).and.(ii.le.mtr)) then
         koff = kyp*(ir - 1)
         call MPI_WAIT(msid,istatus,ierr)
         do 40 k = 1, kypt
         do 30 j = 1, kxpt
         g(1,k+koff,j,l) = t(1,j,k,l)
         g(2,k+koff,j,l) = t(2,j,k,l)
         g(3,k+koff,j,l) = t(3,j,k,l)
   30    continue
   40    continue
      endif
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PWRITE2(f,nx,kyp,nxv,kypmx,nblok,iunit,nrec,lrec,name)
c this subroutine collects distributed real 2d data f and writes to a
c direct access binary file
c f = input data to be written, modified on node 0
c nx/kyp = length of data f in x/y on each processor to write
c nxv = first dimension of data array f, must be >= nx
c kypmx = second dimension of data array f, must be >= kyp
c nblok = number of parallel partitions
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open new file and write first record
c if nrec = 0, open old file, do not write
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: f, nx, kyp, nxv, kypmx, nblok, iunit, nrec, lrec, fname
c output: nrec
      implicit none
      integer nx, kyp, nxv, kypmx, nblok, iunit, nrec, lrec
      real f
      character*(*) name
      dimension f(nxv,kypmx,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, nvp, idproc, np, ioff, id, nrec0, i, j, k, l
      integer ierr
      dimension istatus(lstat)
c this segment is used for shared memory computers
c     if (nrec.lt.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='replace')
c        nrec = 1
c open old file
c     else if (nrec.eq.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='old')
c     endif
c     write (unit=iunit,rec=nrec) (((f(j,k,l),j=1,nx),k=1,kyp),l=1,nblok
c    1)
c     nrec = nrec + 1
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 receives messages from other nodes
      if (idproc.eq.0) then
         if (nrec.lt.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='replace')
            nrec = 1
c open old file
         else if (nrec.eq.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
         endif
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
            id = 1
            call MPI_RECV(f,nxv*kyp,mreal,id,99,lworld,istatus,ierr)
         endif
c first write data for node 0
         nrec0 = nrec
         write (unit=iunit,rec=nrec) (((f(j,k,l),j=1,nx),k=1,kyp),l=1,nb
     1lok)
         nrec = nrec + 1
c then write data from remaining nodes
         do 10 i = 2, np
            id = i - ioff
            call MPI_RECV(f,nxv*kyp,mreal,id,99,lworld,istatus,ierr)
            write (unit=iunit,rec=nrec) (((f(j,k,l),j=1,nx),k=1,kyp),l=1
     1,nblok)
            nrec = nrec + 1
   10    continue
c read data back for node 0
         read (unit=iunit,rec=nrec0) (((f(j,k,l),j=1,nx),k=1,kyp),l=1,nb
     1lok)
c other nodes send data to node 0
      elseif (idproc.le.(nproc+1)) then
         call MPI_SEND(f,nxv*kyp,mreal,0,99,lworld,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PREAD2(f,nx,kyp,nxv,kypmx,nblok,iunit,nrec,lrec,name,ie
     1rror)
c this subroutine reads real 2d data f from a direct access binary file
c and distributes it
c f = output data to be read
c nx/kyp = length of data f in x/y on each processor to read
c nxv = first dimension of data array f, must be >= nx
c kypmx = second dimension of data array f, must be >= kyp
c nblok = number of parallel partitions
c iunit = fortran unit number
c nrec = current record number for read, if nrec > 0
c if nrec < 0, open old file and read first record
c if nrec = 0, open old file, do not read
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c ierror = error indicator
c input: nx, kyp, nxv, kypmx, nblok, iunit, nrec, lrec, fname
c output: f, nrec, ierror
      implicit none
      integer nx, kyp, nxv, kypmx, nblok, iunit, nrec, lrec, ierror
      real f
      character*(*) name
      dimension f(nxv,kypmx,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, nvp, idproc, np, ioff, id, nrec0, i, j, k, l
      integer ierr
      dimension istatus(lstat)
      ierror = 0
c this segment is used for shared memory computers
c     if (nrec.le.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='old')
c        if (nrec.eq.0) return
c        nrec = 1
c     endif
c     read (unit=iunit,rec=nrec,err=10) (((f(j,k,l),j=1,nx),k=1,kyp),l=1
c    1,nblok)
c     nrec = nrec + 1
c     return
c 10 ierror = 1
c     return
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
         if (nrec.le.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
            if (nrec.eq.0) return
            nrec = 1
         endif
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
         endif
c first read data for remaining nodes
         nrec0 = nrec
         nrec = nrec + 1
         do 30 i = 2, np
            read (unit=iunit,rec=nrec,err=10) (((f(j,k,l),j=1,nx),k=1,ky
     1p),l=1,nblok)
            go to 20
   10       ierror = ierror + 1
   20       nrec = nrec + 1
            id = i - ioff
            call MPI_SEND(f,nxv*kyp,mreal,id,98,lworld,ierr)
   30    continue
c then read data from node 0
         read (unit=iunit,rec=nrec0,err=40) (((f(j,k,l),j=1,nx),k=1,kyp)
     1,l=1,nblok)
         go to 50
   40    ierror = ierror + 1
   50    if (ioff.eq.0) then
            id = 1
            call MPI_SEND(f,nxv*kyp,mreal,id,98,lworld,ierr)
         endif
c other nodes receive data from node 0
      else if (idproc.le.(nproc+1)) then
         if (nrec.eq.0) return
         call MPI_RECV(f,nxv*kyp,mreal,0,98,lworld,istatus,ierr)
      endif
c check for error condition
      call PIMAX(ierror,ierr,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine PCWRITE2(f,nx,ny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,nam
     1e)
c this subroutine collects distributed complex 2d data f and writes to a
c direct access binary file
c f = input data to be written, modified on node 0
c nx/ny = total length of data f in x/y to write
c kxp = maximum length of data f in x on each processor to write
c nyv = first dimension of data array f, must be >= ny
c kxpd = second dimension of data array f, must be >= kxp
c jblok = number of parallel partitions
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open new file and write first record
c if nrec = 0, open old file, do not write
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: f, nx, ny, kxp, nyv, kxpd, jblok, iunit, nrec, lrec, name
c output: nrec
      implicit none
      integer nx, ny, kxp, nyv, kxpd, jblok, iunit, nrec, lrec
      complex f
      character*(*) name
      dimension f(nyv,kxpd,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus
      integer ks, kxpp, nvp, idproc, np, ioff, id, nrec0, i, j, k, l
      integer ierr
      dimension istatus(lstat)
      do 50 l = 1, jblok
c this segment is used for shared memory computers
c     if (nrec.lt.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='replace')
c        nrec = 1
c open old file
c     else if (nrec.eq.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='old')
c     endif
c     do 10 j = 1, nx
c     write (unit=iunit,rec=nrec) (f(k,j,l),k=1,ny)
c     nrec = nrec + 1
c  10 continue
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 receives messages from other nodes
      if (idproc.eq.0) then
         if (nrec.lt.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='replace')
            nrec = 1
c open old file
         else if (nrec.eq.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
         return
         endif
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
            kxpp = min(nx,kxp)
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
            id = 1
            call MPI_RECV(f,nyv*kxp,mcplx,id,99,lworld,istatus,ierr)
            call MPI_GET_COUNT(istatus,mcplx,kxpp,ierr)
            kxpp = kxpp/nyv
         endif
c first write data for node 0
         nrec0 = nrec
         do 10 j = 1, kxpp
            write (unit=iunit,rec=nrec) (f(k,j,l),k=1,ny)
            nrec = nrec + 1
   10    continue
c then write data from remaining nodes
         do 30 i = 2, np
            id = i - ioff
            call MPI_RECV(f,nyv*kxp,mcplx,id,99,lworld,istatus,ierr)
            call MPI_GET_COUNT(istatus,mcplx,kxpp,ierr)
            kxpp = kxpp/nyv
            do 20 j = 1, kxpp
               write (unit=iunit,rec=nrec) (f(k,j,l),k=1,ny)
               nrec = nrec + 1
   20       continue
   30    continue
c read data back for node 0
         do 40 j = 1, kxpp
            read (unit=iunit,rec=nrec0) (f(k,j,l),k=1,ny)
            nrec0 = nrec0 + 1
   40    continue 
c other nodes send data to node 0
      elseif (idproc.le.(nproc+1)) then
         if (nrec.eq.0) return
c find amount of data to write
         ks = idproc - 2
         if (nvp.eq.nproc) ks = ks + 1
         kxpp = nx - kxp*(l + ks)
         if (kxpp.gt.kxp) then
            kxpp = kxp
         else if (kxpp.le.0) then
            kxpp = 0
         endif
c send the data
         call MPI_SEND(f,nyv*kxpp,mcplx,0,99,lworld,ierr)
      endif
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCREAD2(f,nx,ny,kxp,nyv,kxpd,jblok,iunit,nrec,lrec,name
     1,ierror)
c this subroutine reads complex 2d data f from a direct access binary
c file and distributes it
c f = output data to be read
c nx/ny = total length of data f in x/y to read
c kxp = maximum length of data f in x on each processor to read
c nyv = first dimension of data array f, must be >= ny
c kxpd = second dimension of data array f, must be >= kxp
c jblok = number of parallel partitions
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open old file and read first record
c if nrec = 0, open old file, do not read
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: nx, ny, kxp, nyv, kxpd, jblok, iunit, nrec, lrec, name
c output: f, nrec, ierror
      implicit none
      integer nx, ny, kxp, nyv, kxpd, jblok, iunit, nrec, lrec
      integer ierror
      complex f
      character*(*) name
      dimension f(nyv,kxpd,jblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus
      integer kxpp, nvp, idproc, np, ioff, id, nrec0, i, j, k, l, ierr
      dimension istatus(lstat)
      ierror = 0
      do 80 l = 1, jblok
c this segment is used for shared memory computers
c     if (nrec.le.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='old')
c        if (nrec.eq.0) return
c        nrec = 1
c     endif
c     do 10 j = 1, nx
c     read (unit=iunit,rec=nrec,err=20) (f(k,j,l),k=1,ny)
c     nrec = nrec + 1
c  10 continue
c     return
c  20 ierror = 1
c     return
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
         if (nrec.le.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
            if (nrec.eq.0) return
            nrec = 1
         endif
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
         endif
c first read data for remaining nodes
         nrec0 = nrec
         nrec = nrec + min(nx,kxp)
         do 40 i = 2, np
            kxpp = nx - kxp*(l + i - 2)
            if (kxpp.gt.kxp) then
               kxpp = kxp
            else if (kxpp.lt.0) then
               kxpp = 0
            endif
            do 10 j = 1, kxpp
               read (unit=iunit,rec=nrec,err=20) (f(k,j,l),k=1,ny)
               nrec = nrec + 1
   10       continue
            go to 30
   20       ierror = ierror + 1
            nrec = nrec - 1
   30       id = i - ioff
            call MPI_SEND(f,nyv*kxpp,mcplx,id,98,lworld,ierr)
   40    continue
c then read data from node 0
         kxpp = min(nx,kxp)
         do 50 j = 1, kxpp
            read (unit=iunit,rec=nrec0,err=60) (f(k,j,l),k=1,ny)
            nrec0 = nrec0 + 1
   50    continue
         go to 70
   60    ierror = ierror + 1
   70    if (ioff.eq.0) then
            id = 1
            call MPI_SEND(f,nyv*kxpp,mcplx,id,98,lworld,ierr)
         endif
c other nodes receive data from node 0
      else if (idproc.le.(nproc+1)) then
         if (nrec.eq.0) return
         call MPI_RECV(f,nyv*kxp,mcplx,0,98,lworld,istatus,ierr)
      endif
   80 continue
c check for error condition
      call PIMAX(ierror,ierr,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine PCDIFF2(f,g,kstrt,nvp,nx,nyp,nxv,nypd,nyp1,nblok)
c this subroutine performs a centered finite difference calculation in
c the second index of a two dimensional array whose second index is
c distributed.  specifically, for 1 < j < nx:
c if 2 < k < nyp, g(j,k,l) = f(j,k+1,l) - f(j,k-1,l)
c if k = 1        g(j,k,l) = f(j,2,l) - f(j,nyp,l-1)
c if k = nyp      g(j,k,l) = f(j,1,l+1) - f(j,nyp-1,l)
c period boundary conditions are assumed for the l subscript
c f = input array
c g = output array
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nx = extent of first index
c nyp = extent of second index
c nxv = first dimension of arrays f and g
c nypd = second dimension of f, >= nyp
c npy1 = second dimension of g, minimum(2,nyp)
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      real f, g
      integer kstrt, nvp, nx, nyp, nxv, nypd, nyp1, nblok
      dimension f(nxv,nypd,nblok), g(nxv,nyp1,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer nyps, ks, kl, kr, msid, j, k, l, ierr, istatus
      dimension istatus(lstat)
      nyps = nyp - 1
      ks = kstrt - 2
      do 20 l = 1, nblok
c find left and right neighbors
      kl = l + ks
      if (kl.lt.1) then
         kl = kl + nvp
      endif
      kr = l + ks + 2
      if (kr.gt.nvp) then
         kr = kr - nvp
      endif
c copy edge k values from left and right
c this segment is for shared memory computers
c     do 10 j = 1, nx
c     g(j,1,l) = f(j,nyp,kl)
c     g(j,nyp1,l) = f(j,1,kr)
c  10 continue
c this segment is used for mpi computers
      call MPI_IRECV(g(1,1,l),nx,mreal,kl-1,l+nyp,lgrp,msid,ierr)
      call MPI_SEND(f(1,nyp,l),nx,mreal,kr-1,l+nyp,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(g(1,nyp1,l),nx,mreal,kr-1,l+nyp+1,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,l),nx,mreal,kl-1,l+nyp+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
   20 continue
c perform centered finite difference
      do 70 l = 1, nblok
c subtract edge k values
      if (nyp.gt.1) then
         do 30 j = 1, nx
         g(j,1,l) = f(j,2,l) - g(j,1,l)
         g(j,nyp,l) = g(j,nyp1,l) - f(j,nyps,l)
   30    continue
c special case of only one k value per processor
      else
         do 40 j = 1, nx
         g(j,1,l) = g(j,nyp1,l) - g(j,1,l)
   40    continue
      endif
c subtract interior k values
      do 60 k = 2, nyps
      do 50 j = 1, nx
      g(j,k,l) = f(j,k+1,l) - f(j,k-1,l)
   50 continue
   60 continue
   70 continue
      return
      end
