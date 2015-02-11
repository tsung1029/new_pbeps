c-----------------------------------------------------------------------      program PNTPOSE2      implicit none      integer indx, indy, indnvp, mshare, ndim      integer nx, ny, nxh, nvp, kyp, kxp, kxb, kyb, kblok, jblok      integer nxv, nyv, nxvhc indnvp = exponent which determines number of virtual processorsc mshare = (0,1) = (no,yes) architecture is shared memory      parameter( indx =   7, indy =   8, indnvp =   2, mshare =   0)      parameter( ndim = 2)      parameter(nx=2**indx,ny=2**indy,nxh=nx/2)      parameter(nvp=2**indnvp,kyp=(ny-1)/nvp+1,kxp=(nxh-1)/nvp+1)      parameter(kxb=nxh/kxp,kyb=ny/kyp)      parameter(kblok=1+mshare*(ny/kyp-1),jblok=1+mshare*(nxh/kxp-1))      parameter(nxv=nx+2,nyv=ny+2,nxvh=nxv/2)      integer i, j, k, l      integer idproc, kstrt, ks, kk, k1, kxp2, joff, j1, jj, koff      integer msid, mrid      dimension msid(kxb), mrid(kyb)      real time, epsmax, eps      real f, t, g, h      dimension f(ndim,nxv,kyp,kblok), t(ndim,2*nyv,kxp,jblok)      dimension g(ndim,nxv,ny), h(ndim,nyv,nx)      double precision ranorm      complex sbuff, bbuff      dimension sbuff(ndim,kxp,kyp,kblok), bbuff(ndim,kxp,kyp,jblok)c initialize for parallel processing      call PPINIT(idproc,nvp)      kstrt = idproc + 1      ks = kstrt - 2c create test function      do 40 k = 1, ny      kk = (k - 1)/kyp      k1 = k - kyp*kk      do 30 j = 1, nx      do 20 i = 1, ndim      g(i,j,k) = ranorm()      h(i,k,j) = g(i,j,k)      do 10 l = 1, kblok      if (kk.eq.(l+ks)) f(i,j,k1,l) = g(i,j,k)   10 continue   20 continue   30 continue   40 continue      call TIMERA(-1,'total   ',time)c start special test casec     call PNTPOSE(f,t,sbuff,bbuff,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypc    1,jblok,kblok,ndim)c     call PNTPOSE(t,f,bbuff,sbuff,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kyp,kxpc    1,kblok,jblok,ndim)cc     call PNTPOSEX(f,t,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kyp,jblok,kbloc    1k,ndim)c     call PNTPOSEX(t,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kyp,kxp,kblok,jbloc    1k,ndim)c      call PNTPOSEY(f,t,msid,mrid,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxb,kyb,     1kxp,kyp,jblok,kblok,ndim)      call PNTPOSEY(t,f,mrid,msid,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kyb,kxb,     1kyp,kxp,kblok,jblok,ndim)cc end special test casec     do 50 i = 1, 100c transform to fourier spacec     isign = -1c     call pfft2r(f,t,sbuff,bbuff,isign,ntpose,mixup,sct,indx,indy,kstrtc    1,nxvh,nyv,kxp,kyp,jblok,kblok,nxhy,nxyh)c     call fft2rx(g,isign,mmixup,ssct,indx,indy,nxvh,nyv,nxhy,nxyh)c transform to real spacec     isign = 1c     call pfft2r(f,t,sbuff,bbuff,isign,ntpose,mixup,sct,indx,indy,kstrtc    1,nxvh,nyv,kxp,kyp,jblok,kblok,nxhy,nxyh)c     call fft2rx(g,isign,mmixup,ssct,indx,indy,nxvh,nyv,nxhy,nxyh)c  50 continue      call TIMERA(1,'total   ',time)c      kxp2 = kxp + kxp      epsmax = 0.      if (kstrt.gt.nxh) go to 100      do 90 l = 1, jblok      joff = kxp2*(l + ks)      do 80 k = 1, ny      do 70 j = 1, kxp2      j1 = j + joff      jj = (j - 1)/2 + 1      if ((j-1).eq.(2*(jj-1))) then         kk = 2*k-1      else         kk = 2*k      endif      do 60 i = 1, ndim      eps = abs(t(i,kk,jj,l) - h(i,k,j1))      if (eps.gt.epsmax) then         write (71,*) i,j1,k,t(i,kk,jj,l),h(i,k,j1),eps         epsmax = eps      endif   60 continue   70 continue   80 continue   90 continue  100 continue      write (71,*) 'first local epsmax=',epsmax      call PSUM(epsmax,eps,1,1)      write (71,*) 'first global epsmax=',epsmaxc      epsmax = 0.      if (kstrt.gt.ny) go to 150      do 140 l = 1, kblok      koff = kyp*(l + ks)      do 130 k = 1, kyp      k1 = k + koff      do 120 j = 1, nx      do 110 i = 1, ndim      eps = abs(f(i,j,k,l) - g(i,j,k1))      if (eps.gt.epsmax) then         write (71,*) i,j,k1,f(i,j,k,l),g(i,j,k1),eps         epsmax = eps      endif  110 continue  120 continue  130 continue  140 continue  150 continue      write (71,*) 'second local epsmax=',epsmax      call PSUM(epsmax,eps,1,1)      write (71,*) 'second global epsmax=',epsmax      call PPEXIT      stop      endc-----------------------------------------------------------------------      subroutine PPINIT(idproc,nvp)c this subroutine initializes parallel processingc input: nvp, output: idprocc idproc = processor idc nvp = number of real or virtual processors requested      implicit none      integer idproc, nvpc get definition of MPI constants      include 'mpif.h'c common block for parallel processing      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworldc lstat = length of status array      parameter(lstat=10)c nproc = number of real or virtual processors obtainedc lgrp = current communicatorc mreal = default datatype for realsc mint = default datatype for integersc mcplx = default datatype for complex type      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworldc local data      integer ierror, ndprec      save /PPARMS/c ndprec = (0,1) = (no,yes) use (normal,autodouble) precision      data ndprec /1/c this segment is used for shared memory computersc     nproc = nvpc     idproc = 0c this segment is used for mpi computers      if (MPI_STATUS_SIZE.gt.lstat) then         write (2,*) ' status size too small, actual/required = ', lstat     1, MPI_STATUS_SIZE         stop      endifc initialize the MPI execution environment      call MPI_INIT(ierror)      if (ierror.ne.0) stop      lgrp = MPI_COMM_WORLDc determine the rank of the calling process in the communicator      call MPI_COMM_RANK(lgrp,idproc,ierror)c determine the size of the group associated with a communicator      call MPI_COMM_SIZE(lgrp,nproc,ierror)c set default datatypes         mint = MPI_INTEGERc single precision      if (ndprec.eq.0) then         mreal = MPI_REAL         mcplx = MPI_COMPLEXc double precision      else         mreal = MPI_DOUBLE_PRECISION         mcplx = MPI_DOUBLE_COMPLEX      endifc requested number of processors not obtained      if (nproc.ne.nvp) then         write (2,*) ' processor number error: nvp, nproc=', nvp, nproc         call PPEXIT         stop      endif      return      endc-----------------------------------------------------------------------      subroutine PPEXITc this subroutine terminates parallel processing      implicit nonec common block for parallel processing      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworldc lgrp = current communicator      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld      integer ierrorc synchronize processes      call MPI_BARRIER(lgrp,ierror)c terminate MPI execution environment      call MPI_FINALIZE(ierror)      return      endc-----------------------------------------------------------------------      subroutine PNTPOSE(f,g,s,t,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,j     1blok,kblok,ndim)c this subroutine performs a transpose of a matrix f, distributed in y,c to a matrix g, distributed in x, that is,c g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), wherec 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kypc and where indices l and m can be distributed across processors.c this subroutine sends and receives one message at a time, eitherc synchronously or asynchronously. it uses a minimum of system resourcesc f = complex input arrayc g = complex output arrayc s, t = complex scratch arraysc nx/ny = number of points in x/yc kstrt = starting data block numberc nxv/nyv = first dimension of f/gc kypd/kxpd = second dimension of f/gc kxp/kyp = number of data values per block in x/yc jblok/kblok = number of data blocks in x/yc ndim = leading dimension of arrays f and g      implicit none      integer nx, ny, kstrt, nxv, nyv, kxp, kyp, kxpd, kypd      integer jblok, kblok, ndim      complex f, g, s, t      dimension f(ndim,nxv,kypd,kblok), g(ndim,nyv,kxpd,jblok)      dimension s(ndim,kxp,kyp,kblok), t(ndim,kxp,kyp,jblok)c common block for parallel processing      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworldc lstat = length of status array      parameter(lstat=10)c lgrp = current communicatorc mcplx = default datatype for complex      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworldc local data      integer ks, kxb, kyb      integer jkblok, kxym, mtr, ntr, mntr      integer l, i, joff, koff, k, j, n      integer ir0, is0, ii, ir, is, ierr, msid, istatus      dimension istatus(lstat)      ks = kstrt - 2      kxb = nx/kxp      kyb = ny/kypc this segment is used for shared memory computersc     if (kstrt.gt.nx) returnc     do 50 l = 1, jblokc     joff = kxp*(l + ks)c     do 40 i = 1, kybc     koff = kyp*(i - 1)c     do 30 k = 1, kypc     do 20 j = 1, kxpc     do 10 n = 1, ndimc     g(n,k+koff,j,l) = f(n,j+joff,k,i)c  10 continuec  20 continuec  30 continuec  40 continuec  50 continuec this segment is used for mpi computers      jkblok = max0(jblok,kblok)      kxym = min0(kxb,kyb)      mtr = kyb/kxym      ntr = kxb/kxym      mntr = max0(mtr,ntr)      do 90 l = 1, jkblok      do 80 i = 1, kxym      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1      is0 = ir0      do 70 ii = 1, mntrc post receive      if ((kstrt.le.nx).and.(ii.le.mtr)) then         ir = ir0 + kxym*(ii - 1)         call MPI_IRECV(t(1,1,1,l),ndim*kxp*kyp,mcplx,ir-1,ir+kxym+1,lgr     1p,msid,ierr)      endifc send data      if ((kstrt.le.ny).and.(ii.le.ntr)) then         is = is0 + kxym*(ii - 1)         joff = kxp*(is - 1)         do 30 k = 1, kyp         do 20 j = 1, kxp         do 10 n = 1, ndim         s(n,j,k,l) = f(n,j+joff,k,l)   10    continue   20    continue   30    continue         call MPI_SEND(s(1,1,1,l),ndim*kxp*kyp,mcplx,is-1,l+ks+kxym+2,lg     1rp,ierr)      endifc receive data      if ((kstrt.le.nx).and.(ii.le.mtr)) then         koff = kyp*(ir - 1)         call MPI_WAIT(msid,istatus,ierr)         do 60 k = 1, kyp         do 50 j = 1, kxp         do 40 n = 1, ndim         g(n,k+koff,j,l) = t(n,j,k,l)   40    continue   50    continue   60    continue      endif   70 continue   80 continue   90 continue      return      endc-----------------------------------------------------------------------      subroutine PNTPOSEX(f,g,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxpd,kypd,jblo     1k,kblok,ndim)c this subroutine performs a transpose of a matrix f, distributed in y,c to a matrix g, distributed in x, that is,c g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), wherec 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kypc and where indices l and m can be distributed across processors.c this subroutine sends and receives multiple asynchronous messages, butc waits after each pair.c f = complex input arrayc g = complex output arrayc nx/ny = number of points in x/yc kstrt = starting data block numberc nxv/nyv = first dimension of f/gc kxp/kyp = number of data values per block in x/yc kypd/kxpd = second dimension of f/gc jblok/kblok = number of data blocks in x/yc ndim = leading dimension of arrays f and gc optimized version      implicit none      integer nx, ny, kstrt, nxv, nyv, kxp, kyp      integer kxpd, kypd, jblok, kblok, ndim      complex f, g      dimension f(ndim*nxv*kypd*kblok), g(ndim*nyv*kxpd*jblok)c common block for parallel processing      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworldc lstat = length of status array      parameter(lstat=10)c lgrp = current communicatorc mcplx = default datatype for complex      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworldc local data      integer ks, kxb, kyb, l, i, joff, koff, k, j, n      integer jkblok, kxym, mtr, ntr, mntr, msid      integer ir0, is0, ii, ir, is, ioff, ierr, istatus      dimension istatus(lstat)      ks = kstrt - 2      kxb = nx/kxp      kyb = ny/kypc this segment is used for shared memory computersc     if (kstrt.gt.nx) returnc     do 50 l = 1, jblokc     joff = kxp*(l + ks) - 1c     do 40 i = 1, kybc     koff = kyp*(i - 1) - 1c     do 30 k = 1, kypc     do 20 j = 1, kxpc     do 10 n = 1, ndimc     g(n+ndim*(k+koff+nyv*(j-1+kxpd*(l-1)))) = f(n+ndim*(j+joff+nxv*(k-c    11+kypd*(i-1))))c  10 continuec  20 continuec  30 continuec  40 continuec  50 continuec this segment is used for mpi computers      jkblok = max0(jblok,kblok)      kxym = min0(kxb,kyb)      mtr = kyb/kxym      ntr = kxb/kxym      mntr = max0(mtr,ntr)c transpose local data      do 50 l = 1, jkblok      ioff = kxb*(l - 1) - 1      koff = kypd*(l - 1) - 1      do 40 i = 1, kxym      is0 = iand(kxym-1,ieor(l+ks,i-1)) + 1      do 30 ii = 1, ntr      if (kstrt.le.ny) then         is = is0 + kxym*(ii - 1)         joff = ndim*kxp*(is - 1)         is = kyp*(is + ioff) - 1         do 20 k = 1, kyp         do 10 j = 1, ndim*kxp         g(j+ndim*kxp*(k+is)) = f(j+joff+ndim*nxv*(k+koff))   10    continue   20    continue      endif   30 continue   40 continue   50 continue      do 80 l = 1, jkblok      ioff = kxb*(l - 1) - 1      koff = kyb*(l - 1) - 1      do 70 i = 1, kxym      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1      is0 = ir0      do 60 ii = 1, mntrc post receive      if ((kstrt.le.nx).and.(ii.le.mtr)) then         ir = ir0 + kxym*(ii - 1)         call MPI_IRECV(f(1+ndim*kxp*kyp*(ir+koff)),ndim*kxp*kyp,mcplx,i     1r-1,ir+kxym+1,lgrp,msid,ierr)      endifc send data      if ((kstrt.le.ny).and.(ii.le.ntr)) then         is = is0 + kxym*(ii - 1)         call MPI_SEND(g(1+ndim*kxp*kyp*(is+ioff)),ndim*kxp*kyp,mcplx,is     1-1,l+ks+kxym+2,lgrp,ierr)      endifc receive data      if ((kstrt.le.nx).and.(ii.le.mtr)) then         call MPI_WAIT(msid,istatus,ierr)      endif   60 continue   70 continue   80 continuec transpose local data      do 140 l = 1, jkblok      ioff = kyb*(l - 1) - 1      joff = kxpd*(l - 1) - 1      do 130 i = 1, kxym      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1      do 120 ii = 1, mtr      if (kstrt.le.nx) then         ir = ir0 + kxym*(ii - 1)         koff = kyp*(ir - 1) - 1         ir = kyp*(ir + ioff) - 1         do 110 k = 1, kyp         do 100 j = 1, kxp         do 90 n = 1, ndim         g(ndim*(k+koff+nyv*(j+joff))+n) = f(ndim*(j+kxp*(k+ir)-1)+n)   90    continue  100    continue  110    continue      endif  120 continue  130 continue  140 continue      return      endc-----------------------------------------------------------------------      subroutine PNTPOSEY(f,g,msid,mrid,nx,ny,kstrt,nxv,nyv,kxp,kyp,kxb,     1kyb,kxpd,kypd,jblok,kblok,ndim)c this subroutine performs a transpose of a matrix f, distributed in y,c to a matrix g, distributed in x, that is,c g(1:ndim,k+kyp*(m-1),j,l) = f(1:ndim,j+kxp*(l-1),k,m), wherec 1 <= j <= kxp, 1 <= k <= kyp, 1 <= l <= nx/kxp, 1 <= m <= ny/kypc and where indices l and m can be distributed across processors,c msid, mrid = scratch arrays for identifying asynchronous messagesc this subroutine sends and receives multiple asynchronous messagesc simultaneously.c f = complex input arrayc g = complex output arrayc nx/ny = number of points in x/yc kstrt = starting data block numberc nxv/nyv = first dimension of f/gc kxp/kyp = number of data values per block in x/yc kxb/kyb = number of processors in x/yc kypd/kxpd = second dimension of f/gc jblok/kblok = number of data blocks in x/yc ndim = leading dimension of arrays f and gc optimized version      implicit none      integer nx, ny, kstrt, nxv, nyv, kxp, kyp      integer kxb, kyb, kxpd, kypd, jblok, kblok, ndim      integer msid, mrid      complex f, g      dimension f(ndim*nxv*kypd*kblok), g(ndim*nyv*kxpd*jblok)      dimension msid(kxb), mrid(kyb)c common block for parallel processing      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworldc lstat = length of status array      parameter(lstat=10)c lgrp = current communicatorc mcplx = default datatype for complex      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworldc local data      integer ks, l, i, joff, koff, k, j, n      integer jkblok, kxym, mtr, ntr, mntr      integer ir0, is0, ii, ir, is, ioff, ierr, istatus      dimension istatus(lstat)      ks = kstrt - 2c this segment is used for shared memory computersc     if (kstrt.gt.nx) returnc     do 50 l = 1, jblokc     joff = kxp*(l + ks) - 1c     do 40 i = 1, kybc     koff = kyp*(i - 1) - 1c     do 30 k = 1, kypc     do 20 j = 1, kxpc     do 10 n = 1, ndimc     g(n+ndim*(k+koff+nyv*(j-1+kxpd*(l-1)))) = f(n+ndim*(j+joff+nxv*(k-c    11+kypd*(i-1))))c  10 continuec  20 continuec  30 continuec  40 continuec  50 continuec this segment is used for mpi computers      jkblok = max0(jblok,kblok)      kxym = min0(kxb,kyb)      mtr = kyb/kxym      ntr = kxb/kxym      mntr = max0(mtr,ntr)c transpose local data      do 50 l = 1, jkblok      ioff = kxb*(l - 1) - 1      koff = kypd*(l - 1) - 1      do 40 i = 1, kxym      is0 = iand(kxym-1,ieor(l+ks,i-1)) + 1      do 30 ii = 1, ntr      if (kstrt.le.ny) then         is = is0 + kxym*(ii - 1)         joff = ndim*kxp*(is - 1)         is = kyp*(is + ioff) - 1         do 20 k = 1, kyp         do 10 j = 1, ndim*kxp         g(j+ndim*kxp*(k+is)) = f(j+joff+ndim*nxv*(k+koff))   10    continue   20    continue      endif   30 continue   40 continue   50 continuec post all receives      do 80 l = 1, jkblok      koff = kyb*(l - 1) - 1      do 70 i = 1, kxym      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1      do 60 ii = 1, mntr      if ((kstrt.le.nx).and.(ii.le.mtr)) then         ir = ir0 + kxym*(ii - 1)         call MPI_IRECV(f(1+ndim*kxp*kyp*(ir+koff)),ndim*kxp*kyp,mcplx,i     1r-1,ir+kxym+1,lgrp,mrid(ir),ierr)      endif   60 continue   70 continue   80 continuec post all sends      do 110 l = 1, jkblok      ioff = kxb*(l - 1) - 1      do 100 i = 1, kxym      is0 = iand(kxym-1,ieor(l+ks,i-1)) + 1      do 90 ii = 1, mntr      if ((kstrt.le.ny).and.(ii.le.ntr)) then         is = is0 + kxym*(ii - 1)         call MPI_ISEND(g(1+ndim*kxp*kyp*(is+ioff)),ndim*kxp*kyp,mcplx,i     1s-1,l+ks+kxym+2,lgrp,msid(is),ierr)      endif   90 continue  100 continue  110 continuec wait for all messages      do 140 l = 1, jkblok      do 130 i = 1, kxym      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1      is0 = ir0      do 120 ii = 1, mntrc wait for data to arrive      if ((kstrt.le.nx).and.(ii.le.mtr)) then         ir = ir0 + kxym*(ii - 1)         call MPI_WAIT(mrid(ir),istatus,ierr)      endifc make sure data was successfully sent      if ((kstrt.le.ny).and.(ii.le.ntr)) then         is = is0 + kxym*(ii - 1)         call MPI_WAIT(msid(is),istatus,ierr)      endif  120 continue  130 continue  140 continuec transpose local data      do 200 l = 1, jkblok      ioff = kyb*(l - 1) - 1      joff = kxpd*(l - 1) - 1      do 190 i = 1, kxym      ir0 = iand(kxym-1,ieor(l+ks,i-1)) + 1      do 180 ii = 1, mtr      if (kstrt.le.nx) then         ir = ir0 + kxym*(ii - 1)         koff = kyp*(ir - 1) - 1         ir = kyp*(ir + ioff) - 1         do 170 k = 1, kyp         do 160 j = 1, kxp         do 150 n = 1, ndim         g(ndim*(k+koff+nyv*(j+joff))+n) = f(ndim*(j+kxp*(k+ir)-1)+n)  150    continue  160    continue  170    continue      endif  180 continue  190 continue  200 continue      return      endc-----------------------------------------------------------------------      function ranorm()c this program calculates a random number y from a gaussian distributionc with zero mean and unit variance, according to the method ofc mueller and box:c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),c where x is a random number uniformly distributed on (0,1).c written for the ibm by viktor k. decyk, ucla      integer r1,r2,r4,r5      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/      data iflg,r0 /0,0.0d0/      if (iflg.eq.0) go to 10      ranorm = r0      r0 = 0.0d0      iflg = 0      return   10 isc = 65536      asc = dble(isc)      bsc = asc*asc      i1 = r1 - (r1/isc)*isc      r3 = h1l*dble(r1) + asc*h1u*dble(i1)      i1 = r3/bsc      r3 = r3 - dble(i1)*bsc      bsc = 0.5d0*bsc      i1 = r2/isc      isc = r2 - i1*isc      r0 = h1l*dble(r2) + asc*h1u*dble(isc)      asc = 1.0d0/bsc      isc = r0*asc      r2 = r0 - dble(isc)*bsc      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))      isc = r3*asc      r1 = r3 - dble(isc)*bsc      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))      isc = 65536      asc = dble(isc)      bsc = asc*asc      i1 = r4 - (r4/isc)*isc      r3 = h2l*dble(r4) + asc*h1u*dble(i1)      i1 = r3/bsc      r3 = r3 - dble(i1)*bsc      bsc = 0.5d0*bsc      i1 = r5/isc      isc = r5 - i1*isc      r0 = h2l*dble(r5) + asc*h1u*dble(isc)      asc = 1.0d0/bsc      isc = r0*asc      r5 = r0 - dble(isc)*bsc      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))      isc = r3*asc      r4 = r3 - dble(isc)*bsc      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)      ranorm = temp*dsin(r0)      r0 = temp*dcos(r0)      iflg = 1      return      endc-----------------------------------------------------------------------      subroutine TIMERA(icntrl,chr,time)c this subroutine performs timingc input: icntrl, chrc icntrl = (-1,0,1) = (initialize,ignore,read) clockc clock should be initialized before it is read!c chr = character variable for labeling timingsc time = elapsed time in secondsc written for mpi      implicit none      integer icntrl      character*8 chr      real timec get definition of MPI constants      include 'mpif.h'c common block for parallel processing      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworldc lgrp = current communicatorc mreal = default datatype for reals      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworldc local data      integer idproc, ierr      real nclock, mclock      double precision jclock      save jclock   91 format (1x,a8,1x,'max/min real time = ',e14.7,1x,e14.7,1x,'sec')      data jclock /0.0d0/      if (icntrl.eq.0) return      if (icntrl.eq.1) go to 10c initialize clock      call MPI_BARRIER(lgrp,ierr)      jclock = MPI_WTIME()      returnc read clock and write time difference from last clock initialization   10 nclock = real(MPI_WTIME() - jclock)      call MPI_ALLREDUCE(nclock,time,1,mreal,MPI_MIN,lgrp,ierr)      mclock = time      call MPI_ALLREDUCE(nclock,time,1,mreal,MPI_MAX,lgrp,ierr)      call MPI_COMM_RANK(lgrp,idproc,ierr)      if (idproc.eq.0) write (6,91) chr, time, mclock      return      endc-----------------------------------------------------------------------      subroutine PSUM(f,g,nxp,nblok)c this subroutine performs a parallel sum of a vector, that is:c f(j,k) = sum over k of f(j,k)c assumes the number of processors nproc is a power of two.c the algorithm performs partial sums in binary pairs, as follows:c first, adjacent processors exchange vectors and sum them.  next,c processors separated by 2 exchange the new vectors and sum them, thenc those separated by 4, up to processors separated by nproc/2.  at thec end, all processors contain the same summation.c f = input and output datac g = scratch arrayc nxp = number of data values in vectorc nblok = number of data blocksc written by viktor k. decyk, ucla      implicit none      real f, g      integer nxp, nblok      dimension f(nxp,nblok), g(nxp,nblok)c common block for parallel processing      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworldc lstat = length of status array      parameter(lstat=10)c nproc = number of real or virtual processors obtainedc lgrp = current communicatorc mreal = default datatype for reals      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworldc local data      integer istatus, ierr, msid      integer idproc, kstrt, ks, l, kxs, k, kb, lb, j      dimension istatus(lstat)c find processor idc this line is used for shared memory computersc     idproc = 0c this line is used for mpi computers      call MPI_COMM_RANK(lgrp,idproc,ierr)      kstrt = idproc + 1      if (kstrt.gt.nproc) return      ks = kstrt - 2      l = 1      kxs = 1c main iteration loop   10 if (kxs.ge.nproc) go to 60c shift data      do 30 k = 1, nblok      kb = k + ks      lb = kb/kxs      kb = kb + 1      lb = lb - 2*(lb/2)c this loop is used for shared memory computersc     do 20 j = 1, nxpc     if (lb.eq.0) thenc        g(j,k) = f(j,kb+kxs)c     elsec        g(j,k) = f(j,kb-kxs)c     endifc  20 continuec this segment is used for mpi computers      if (lb.eq.0) then         call MPI_IRECV(g,nxp,mreal,kb+kxs-1,l+nxp,lgrp,msid,ierr)         call MPI_SEND(f,nxp,mreal,kb+kxs-1,l+nxp,lgrp,ierr)      else         call MPI_IRECV(g,nxp,mreal,kb-kxs-1,l+nxp,lgrp,msid,ierr)         call MPI_SEND(f,nxp,mreal,kb-kxs-1,l+nxp,lgrp,ierr)      endif      call MPI_WAIT(msid,istatus,ierr)   30 continuec perform sum      do 50 k = 1, nblok      do 40 j = 1, nxp      f(j,k) = f(j,k) + g(j,k)   40 continue   50 continue      l = l + 1      kxs = kxs + kxs      go to 10   60 return      end