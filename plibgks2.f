c general 2d parallel gks graphics library
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: april 30, 2013
c-----------------------------------------------------------------------
      subroutine PGRCLOSE
c this subroutine deactivates workstation and closes gks
c idwk = workstation identifier
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
      parameter(lstat=10)
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c pause if plots are still pending
      if (((iplot.ne.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c this segment is used for shared memory computers
c        idproc = 0
c this segment is used for mpi computers
         call MPI_COMM_RANK(lworld,idproc,ierr)
c read code from input device, if present
         if (idproc.eq.0) call readrc(irc)
c this segment is used for mpi computers
         call PBICAST(irc,1)
      endif
c deactivate workstation
      call gdawk(idwk)
c close workstation
      call gclwk(idwk)
c close gks
      call gclks
      return
      end
c-----------------------------------------------------------------------
      subroutine PCARPET(f,g,nvp,label,isc,ist,nx,ny,nxv,nypmx,nblok,chr
     1,ntc,irc)
c this subroutine displays an array f as a color raster image, for
c distributed data
c a 256 color palette must have been defined prior to this call. 
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = distributed field array to be plotted
c g = scratch array for receiving messages
c nvp = number of real or virtual processors requested
c label = long character string label for plot
c isc = power of 2 scale of range of values of f
c ist = flag for choosing positive and/or negative values
c the range of values of f are given by fmax and fmin.
c if ist = 0, then fmax = 2**isc and fmin = -2**isc.
c if ist = 1, then fmax = 2**isc and fmin = 0.
c if ist = -1, then fmax = 0 and fmin = -2**isc.
c if ist = 2, then fmax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells
c nblok = number of particle partitions
c chr = additional long character string comment for plot
c ntc = number of valid colors, should be power of 2, <= 256
c irc = return code (0 = normal return)
      implicit none
      real f, g
      integer nvp, isc, ist, nx, ny, nxv, nypmx, nblok, ntc, irc
      character*(*) label, chr
      dimension f(nxv,nypmx,nblok), g(nxv,nypmx)
c data in common blocks
      integer idwk, ncols, iplot, nplot, iclr, iupd, idstr, idloc, nclsp
      integer ifrg, isx, isy, kprime
      real rx, ry
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
c (default=when iplot=nplot-1)
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
      integer npald, lxm, lym, lupt, ipal, img8
c npald = number of palette entries
      parameter(npald=256)
c lupt = (0,1) = (no,yes) pixel lookup table needed
c ipal = integer pixel lookup table
      dimension ipal(npald)
c lxm, lym = maximum number of pixels in x, y
      parameter(lxm=720,lym=540)
c img8 = integer image array
      dimension img8(lxm*lym)
      common /movicm/ lupt,ipal,img8
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c mreal = default datatype for reals
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msg
      integer ntx, nty, istyle, nyp1, is, j, k, l, it, ir, npl1, npl
      integer iy, ix, id, lxs, lys, lxp, lyp, ncv, ic, i1
      integer joff, idproc, i, lvp, ierr
      real gmax, hmax
      real dv, smin, tmin, smax, tmax, csize, algdvi, fmax, fmin
      real rmax, rmin, xmin, xmax, ymin, ymax, apl, sx, sy, orx, ory
      real smn, smx, tmn, tmx, chh, xmn, ymn, ac
      double precision range
      dimension istatus(lstat)
      dimension msg(2)
      dimension gmax(2), hmax(2)
      dimension range(2)
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /11,11/
c csize = vertical size of characters
      data csize /0.034/
c istyle = (0,1) = color map (fills area,preserves aspect ratio)
      data istyle /1/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
      nyp1 = (ny - 1)/nvp + 2
c this segment is used for shared memory computers
c     idproc = 0
c     lvp = nvp
c this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
c find scales for plot
      is = isc
c nodes with data find range
      if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
         fmax = f(1,1,1)
         fmin = fmax
         do 30 l = 1, nblok
         do 20 k = 1, nyp1
         do 10 j = 1, nx
         fmax = amax1(fmax,f(j,k,l))
         fmin = amin1(fmin,f(j,k,l))
   10    continue
   20    continue
   30    continue
c this segment is used for shared memory computers
         hmax(1) = -fmin
         hmax(2) = fmax
c this line is used for mpi computers
         call PMAX(hmax,gmax,2,1)
         fmin = -hmax(1)
         fmax = hmax(2)
         if (fmax.eq.0.) fmax = 1.0e-35
         rmax = fmax - fmin
         if (rmax.eq.0.) rmax = 1.0e-35
         rmin = fmin
         ymax = abs(fmax)
         is = alog(ymax)*algdvi
         if (ymax.ge.1.) is = is + 1
         if (ymax.le.dv**(is-1)) is = is - 1
         ymin = abs(fmin)  
         if (ymin.gt.0.) then
            it = alog(ymin)*algdvi
            if (ymin.ge.1.) it = it + 1
            if (ymin.le.dv**(it-1)) it = it - 1
         endif
         if (fmax.gt.0.) then
            if (fmin.gt.0.) then
               fmin = dv**(it - 1)
            else if (fmin.lt.0.) then
               fmin = -dv**it
            endif
            fmax = dv**is
         else
            fmax = -dv**(is - 1)
            fmin = -dv**it
         endif
         if (ist.eq.0) then
            if (ymin.gt.ymax) then
               fmax = dv**it
            else
               fmin = -dv**is
            endif
         else if (ist.eq.2) then
            ir = alog(rmax)*algdvi
            if (rmax.ge.1.) ir = ir + 1
            if (rmax.le.dv**(ir-1)) ir = ir - 1
            fmin = rmin
            fmax = rmin + dv**ir
         endif
      else
         fmax = dv**is
         fmin = -fmax
      endif
c broadcast range to diagnostic node
      range(1) = fmin
      range(2) = fmax
      call HARTBEAT(range,2)
      fmin = range(1)
      fmax = range(2)
c clip range if necessary
      if (ist.eq.1) then
         fmin = 0.
      else if (ist.eq.(-1)) then
         fmax = 0.  
      endif
c parameters for plots
      xmin = 0.
      xmax = float(nx - 1)
      ymin = 0.
      ymax = float(ny)
c find location for plot
      npl1 = sqrt(float(nplot-1)) + 0.0001
      npl = npl1 + 1
      apl = 1./float(npl)
      iy = iplot/npl
      ix = iplot - iy*npl
      sx = apl*rx
      sy = apl*ry
      orx = sx*float(ix)
      ory = sy*float(npl1 - iy)
      smn = orx + sx*smin
      smx = orx + sx*smax
      tmn = ory + sy*tmin
      tmx = ory + sy*tmax
      chh = sy*csize
c fill area
      xmn = smn
      ymn = tmn
c preserve aspect ratio
      if (istyle.eq.1) then
         if (nx.gt.ny) ymn = tmx - (tmx - tmn)*float(ny)/float(nx)
         if (ny.gt.nx) xmn = smx - (smx - smn)*float(nx)/float(ny)
      endif
c set size of raster image
      lxs = lxm
      lys = lym
      lxp = min(2,lxs)
      lyp = lys/nvp
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c create sample palette
      ac = float(ntc - 1)/float(lys - 1)
c rescale factor ncv
      ncv = 256/ntc
      do 50 k = 1, lys
      joff = lxp*(k - 1)
      ic = ac*float(lys - k) + 0.999999
c rescale index
      ic = ic*ncv
c lookup table required
      if (lupt.eq.1) ic = ipal(ic+1)
      do 40 j = 1, lxp
      img8(j+joff) = ic
   40 continue
   50 continue
c draw grid and labels, call identity transformation
      call tickz(xmin,xmax,ymin,ymax,orx,ory,xmn,smx,ymn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample palette
      call dsmpal(img8,fmin,fmax,orx,smn,tmn,tmx,lxp,lys,chh)
c map parts of f to color raster image
      ic = lyp
      do 80 l = 1, nblok
c this segment is used for shared memory computers
c     do 70 k = 1, nyp1
c     do 60 j = 1, nx
c     g(j,k) = f(j,k,l)
c  60 continue
c  70 continue
c map f to color raster image
c     i1 = lxs*(lys - lyp*l) + 1
c     if (l.eq.nblok) ic = lyp
c     call mraster(g,img8(i1),fmin,fmax,nx,nyp1,nxv,lxs,ic,ntc)
c this segment is used for mpi computers
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
      if (idproc.eq.0) then
c special diagnostic node present
         if (joff.eq.0) then
            id = 1
            call MPI_RECV(f,nxv*nyp1,mreal,id,101,lworld,istatus,ierr)
         endif
c processor 0 maps his (or her) own data
         i1 = lxs*(lys - lyp) + 1
         call mraster(f,img8(i1),fmin,fmax,nx,nyp1,nxv,lxs,ic,ntc)
c then collects data from remaining nodes
         do 60 i = 2, nvp
         id = i - joff
         call MPI_RECV(g,nxv*nyp1,mreal,id,101,lworld,istatus,ierr)
         i1 = lxs*(lys - lyp*i) + 1
c then map the remote data
         if (i.eq.nvp) ic = lyp
         call mraster(g,img8(i1),fmin,fmax,nx,nyp1,nxv,lxs,ic,ntc)
   60    continue
c other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         id = idproc + joff
         call MPI_SEND(f,nxv*nyp1,mreal,0,101,lworld,ierr)
         i1 = lxs*(lys - lyp*id) + 1
c then map the local data
         if (id.eq.nproc) ic = lyp
         call mraster(f,img8(i1),fmin,fmax,nx,nyp1,nxv,lxs,ic,ntc)
      endif
   80 continue
c copy with lookup table
      if (lupt.eq.1) then
         do 100 k = 1, lys
         joff = lxs*(k - 1)
         do 90 j = 1, lxs
         img8(j+joff) = ipal(img8(j+joff)+1)
   90    continue
  100    continue
      endif
c cell array
c special case for rs/6000 with graPHIGS gks
c     call gca(xmn,tmx,smx,ymn,lys,lxs,1,1,lys,lxs,img8)
      call gca(xmn,tmx,smx,ymn,lxs,lys,1,1,lxs,lys,img8)
c update workstation, perform
      call guwk(idwk,1)
c update plot number
      iplot = iplot + 1
      if (iplot.eq.nplot) iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         if (idproc.eq.0) call readrc(irc)
c this segment is used for mpi computers
         msg(1) = irc
         msg(2) = nplot
         call PBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
c reset defaults
      iclr = 0
      iupd = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine PCONTUR(f,g,lf,nvp,label,isc,ist,nx,ny,nxv,nypmx,nblok,
     1chr,nc,irc)
c this subroutine displays an array f as a contour plot.
c a maximum of ncols colors are used, used in order from lowest to
c highest contour: blue, green, cyan, foreground, yellow, magenta, red
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = field array to be plotted
c lf = scratch field array
c g = scratch array for receiving messages
c nvp = number of real or virtual processors requested
c label = long character string label for plot
c isc = power of 2 scale of range of values of f
c ist = flag for choosing positive and/or negative values
c the range of values of f are given by fmax and fmin.
c if ist = 0, then fmax = 2**isc and fmin = -2**isc.
c if ist = 1, then fmax = 2**isc and fmin = 0.
c if ist = -1, then fmax = 0 and fmin = -2**isc.
c if ist = 2, then fmax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c nx/ny = length of field f in x/y direction
c nxv = first dimension of field array f, must be >= nx
c nypmx = maximum size of particle partition, including guard cells
c nblok = number of particle partitions
c chr = additional long character string comment for plot
c nc = number of contour lines
c irc = return code (0 = normal return)
      implicit none
      real f, g
      integer lf
      integer nvp, isc, ist, nx, ny, nxv, nypmx, nblok, nc, irc
      character*(*) label, chr
      dimension f(nxv,nypmx,nblok), lf(nxv,nypmx,nblok), g(nxv,nypmx)
c data in common blocks
      integer idwk, ncols, iplot, nplot, iclr, iupd, idstr, idloc, nclsp
      integer ifrg, isx, isy, kprime
      real rx, ry
c idwk = workstation identifier
c ncols = number of foreground colors available for line plotting
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
c (default=when iplot=nplot-1)
c kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c mreal = default datatype for reals
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, icv, icolor, msg
      integer ntx, nty, istyle, nyp1, is, j, k, l, it, ir, npl1, npl
      integer iy, ix, id, ntc, joff, idproc, i, lvp, ierr
      real gmax, hmax
      real dv, smin, tmin, smax, tmax, csize, algdvi, fmax, fmin
      real rmax, rmin, xmin, xmax, ymin, ymax, apl, sx, sy, orx, ory
      real smn, smx, tmx, tmn, chh, xmn, ymn, dyp, ymnp, tmxp
      double precision range
      dimension istatus(lstat)
c icolor = color index, used in order from lowest to highest contour
      dimension icv(8), icolor(8), msg(2)
      dimension gmax(2), hmax(2)
      dimension range(2)
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /11,11/
c csize = vertical size of characters
      data csize /0.034/
c icv = location in kprime array of color indices needed for icolor
      data icv /1,3,8,6,2,5,7,4/
c istyle = (0,1) = contour plot (fills area,preserves aspect ratio)
      data istyle /1/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
      nyp1 = (ny - 1)/nvp + 2
c this segment is used for shared memory computers
c     idproc = 0
c     lvp = nvp
c this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
c find scales for plot
      is = isc
c nodes with data find range
      if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
         fmax = f(1,1,1)
         fmin = fmax
         do 30 l = 1, nblok
         do 20 k = 1, nyp1
         do 10 j = 1, nx
         fmax = amax1(fmax,f(j,k,l))
         fmin = amin1(fmin,f(j,k,l))
   10    continue
   20    continue
   30    continue
c this segment is used for shared memory computers
         hmax(1) = -fmin
         hmax(2) = fmax
c this line is used for mpi computers
         call PMAX(hmax,gmax,2,1)
         fmin = -hmax(1)
         fmax = hmax(2)
         if (fmax.eq.0.) fmax = 1.0e-35
         rmax = fmax - fmin
         if (rmax.eq.0.) rmax = 1.0e-35
         rmin = fmin
         ymax = abs(fmax)
         is = alog(ymax)*algdvi
         if (ymax.ge.1.) is = is + 1
         if (ymax.le.dv**(is-1)) is = is - 1
         ymin = abs(fmin)  
         if (ymin.gt.0.) then
            it = alog(ymin)*algdvi
            if (ymin.ge.1.) it = it + 1
            if (ymin.le.dv**(it-1)) it = it - 1
         endif
         if (fmax.gt.0.) then
            if (fmin.gt.0.) then
               fmin = dv**(it - 1)
            else if (fmin.lt.0.) then
               fmin = -dv**it
            endif
            fmax = dv**is
         else
            fmax = -dv**(is - 1)
            fmin = -dv**it
         endif
         if (ist.eq.0) then
            if (ymin.gt.ymax) then
               fmax = dv**it
            else
               fmin = -dv**is
            endif
         else if (ist.eq.2) then
            ir = alog(rmax)*algdvi
            if (rmax.ge.1.) ir = ir + 1
            if (rmax.le.dv**(ir-1)) ir = ir - 1
            fmin = rmin
            fmax = rmin + dv**ir
         endif
      else
         fmax = dv**is
         fmin = -fmax
      endif
c broadcast range to diagnostic node
      range(1) = fmin
      range(2) = fmax
      call HARTBEAT(range,2)
      fmin = range(1)
      fmax = range(2)
c clip range if necessary
      if (ist.eq.1) then
         fmin = 0.
      else if (ist.eq.(-1)) then
         fmax = 0.  
      endif
c parameters for plots
      xmin = 0.
      xmax = float(nx - 1)
      ymin = 0.
      ymax = float(ny)
c find location for plot
      npl1 = sqrt(float(nplot-1)) + 0.0001
      npl = npl1 + 1
      apl = 1./float(npl)
      iy = iplot/npl
      ix = iplot - iy*npl
      sx = apl*rx
      sy = apl*ry
      orx = sx*float(ix)
      ory = sy*float(npl1 - iy)
      smn = orx + sx*smin
      smx = orx + sx*smax
      tmn = ory + sy*tmin
      tmx = ory + sy*tmax
      chh = sy*csize
c fill area
      xmn = smn
      ymn = tmn
c preserve aspect ratio
      if (istyle.eq.1) then
         if (nx.gt.ny) ymn = tmx - (tmx - tmn)*float(ny)/float(nx)
         if (ny.gt.nx) xmn = smx - (smx - smn)*float(nx)/float(ny)
      endif
      dyp = (tmx - ymn)/float(nvp)
c calculate color indices
      ntc = ncols + 1
      do 40 i = 1, ntc
      icolor(i) =  kprime(icv(i))
   40 continue
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c draw grid and labels, call identity transformation
      call tickz(xmin,xmax,ymin,ymax,orx,ory,xmn,smx,ymn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample contours
      call dsmpcn(icolor,fmin,fmax,orx,smn,tmn,tmx,chh,ntc,nc)
c display parts of f as contour plot
      do 70 l = 1, nblok
c this segment is used for shared memory computers
c     do 60 k = 1, nyp1
c     do 50 j = 1, nx
c     g(j,k) = f(j,k,l)
c  50 continue
c  60 continue
c     ymnp = ymn + dyp*float(l - 1)
c     tmxp = ymnp + dyp
c draw contour plot
c     call dcontr(g,lf,icolor,fmin,fmax,xmn,smx,ymnp,tmxp,nx,nyp1,nxv,nt
c    1c,nc)
c  70 continue
c this segment is used for mpi computers
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
      if (idproc.eq.0) then
c special diagnostic node present
         if (joff.eq.0) then
            id = 1
            call MPI_RECV(f,nxv*nyp1,mreal,id,101,lworld,istatus,ierr)
         endif
c processor 0 maps his (or her) own data
         tmxp = ymn + dyp
c draw contour map
         call dcontr(f,lf,icolor,fmin,fmax,xmn,smx,ymn,tmxp,nx,nyp1,nxv,
     1ntc,nc)
c then collects data from remaining nodes
         do 50 i = 2, nvp
         id = i - joff
         call MPI_RECV(g,nxv*nyp1,mreal,id,101,lworld,istatus,ierr)
         ymnp = ymn + dyp*float(i - 1)
         tmxp = ymnp + dyp
c then draw contour map of remote data
         call dcontr(g,lf,icolor,fmin,fmax,xmn,smx,ymnp,tmxp,nx,nyp1,nxv
     1,ntc,nc)
   50    continue
c other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         id = idproc + joff - 1
         call MPI_SEND(f,nxv*nyp1,mreal,0,101,lworld,ierr)
         ymnp = ymn + dyp*float(id)
         tmxp = ymnp + dyp
c then draw contour map of local data
         call dcontr(f,lf,icolor,fmin,fmax,xmn,smx,ymnp,tmxp,nx,nyp1,nxv
     1,ntc,nc)
      endif
   70 continue
c update workstation, perform
      call guwk(idwk,1)
c update plot number
      iplot = iplot + 1
      if (iplot.eq.nplot) iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         if (idproc.eq.0) call readrc(irc)
c this segment is used for mpi computers
         msg(1) = irc
         msg(2) = nplot
         call PBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
c reset defaults
      iclr = 0
      iupd = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRASP23(part,f,npp,label,itime,isc,nx,ny,iyp,ixp,idimp
     1,npmax,nblok,irc)
c for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
c part(1,n,l) = position x of particle n in partition l
c part(2,n,l) = position y of particle n in partition l
c part(3,n,l) = velocity vx of particle n in partition l
c part(4,n,l) = velocity vy of particle n in partition l
c part(5,n,l) = velocity vz of particle n in partition l
c f = scratch array for receiving messages
c npp(l) = number of particles in partition l
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c nx/ny = system length in x/y direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 4 or 5
c npmax = maximum number of particles in each partition
c nblok = number of particle partitions
c irc = return code (0 = normal return)
      implicit none
      real part, f
      integer npp
      integer itime, isc, nx, ny, iyp, ixp, idimp, npmax, nblok, irc
      character*(*) label
      dimension part(idimp,npmax,nblok), f(2*npmax), npp(nblok)
c data in common blocks
      integer idwk, ncols, iplot, nplot, iclr, iupd, idstr, idloc, nclsp
      integer ifrg, isx, isy, kprime
      real rx, ry
c idwk = workstation identifier
c ncols = number of foreground colors available for line plotting
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
c (default=when iplot=nplot-1)
c kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c mint = default datatype for integers
c mreal = default datatype for reals
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msg
      integer ntx, nty, is, npl1, npl, iy, ix, mks, ngs, nrt
      integer i, j, l, nd, koff, icol, joff, id, idproc, nvp, lvp, ierr
      real dv, smin, tmin, smax, tmax, zero, csize, algdvi, apl
      real fmax, ymin, ymax, xmin, xmax, aplx, aply, orx, ory, smn, smx
      real tmn, tmx, chh, dd
      double precision range
      dimension istatus(lstat)
      dimension msg(2)
      character*26 lbl
      character*2 lblsp(5)
      character*10 chrs(2)
      dimension range(4)
      save lblsp
   91 format (1x,a2,' VERSUS ',a2,', T = ',i7)
      data lblsp /' X',' Y','VX','VY','VZ'/
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /9,9/
c csize = vertical size of characters
      data zero,csize /0.,0.034/
c sample labels
      data chrs /'BACKGROUND','   BEAM   '/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
c this segment is used for shared memory computers
c     idproc = 0
c     lvp = nproc
c this segment is used for mpi computers
      call MPI_COMM_RANK(lworld,idproc,ierr)
      call MPI_COMM_SIZE(lworld,lvp,ierr)
c find y scale for plot
      if (iyp.le.1) then
         if (iyp.lt.1) iyp = 1
         ymin = zero
         ymax = float(nx)
      elseif (iyp.eq.2) then
         ymin = zero
         ymax = float(ny)
      elseif (iyp.gt.2) then
         if (iyp.gt.idimp) iyp = idimp
         is = isc
         if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
            fmax = abs(part(iyp,1,1))
            do 20 l = 1, nblok
            do 10 j = 1, npp(l)
            fmax = amax1(fmax,abs(part(iyp,j,l)))
   10       continue
   20       continue
c this line is used for mpi computers
            call PMAX(fmax,ymax,1,1)
            ymax = fmax
            if (ymax.eq.0.) ymax = 1.0e-35
            is = alog(ymax)*algdvi
            if (ymax.ge.1.) is = is + 1
            if (ymax.le.dv**(is-1)) is = is - 1
         endif
         ymax = dv**is
         ymin = -ymax
      endif
c find x scale for plot
      if (ixp.le.1) then
         if (ixp.lt.1) ixp = 1
         xmin = zero
         xmax = float(nx)
      elseif (ixp.eq.2) then
         xmin = zero
         xmax = float(ny)
      elseif (ixp.gt.2) then
         if (ixp.gt.idimp) ixp = idimp
         is = isc
         if ((abs(is).gt.116).and.(lvp.le.(2*nproc))) then
            fmax = abs(part(ixp,1,1))
            do 40 l = 1, nblok
            do 30 j = 1, npp(l)
            fmax = amax1(fmax,abs(part(ixp,j,l)))
   30       continue
   40       continue
c this segment is used for mpi computers
            call PMAX(fmax,xmax,1,1)
            xmax = fmax
            if (xmax.eq.0.) xmax = 1.0e-35
            is = alog(xmax)*algdvi
            if (xmax.ge.1.) is = is + 1
            if (xmax.le.dv**(is-1)) is = is - 1
         endif
         xmax = dv**is
         xmin = -xmax
      endif
c broadcast range to diagnostic node
      range(1) = ymax
      range(2) = ymin
      range(3) = xmax
      range(4) = xmin
      call HARTBEAT(range,2)
      ymax = range(1)
      ymin = range(2)
      xmax = range(3)
      xmin = range(4)
c find location for plot
      npl1 = sqrt(float(nplot-1)) + 0.0001
      npl = npl1 + 1
      apl = 1./float(npl)
      iy = iplot/npl
      ix = iplot - iy*npl
      aplx = apl*rx
      aply = apl*ry
      orx = aplx*float(ix)
      ory = aply*float(npl1 - iy)
      smn = orx + aplx*smin
      smx = orx + aplx*smax
      tmn = ory + aply*tmin
      tmx = ory + aply*tmax
      chh = aply*csize
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c write labels
      write (lbl,91) lblsp(iyp), lblsp(ixp), itime
c select point as marker
      mks = 1
c draw grid and labels, call identity transformation
      call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,lab
     1el,lbl,chh)
c display sample markers with dots of minimum visible size
      dd = (smax - smin)*float(npl)*2.0e-3
      ngs = 1
c     call dsmpln(orx,smn,tmn,tmx,ngs,1,1,mks,chrs,chh)
      call dsdmpln(orx,smn,tmn,tmx,dd,ngs,1,1,mks,chrs,chh)
c define transformation number 2
      nrt = 1
c set window
      call gswn(nrt,xmin,xmax,ymin,ymax)
c set viewport
      call gsvp(nrt,smn,smx,tmn,tmx)
c select normalization transformation
      call gselnt(nrt)
c set marker size scale factor, 1.0 = nominal
      call gsmksc(1.0)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
      call gsmk(mks)
c set clipping indicator, 1 = on
      call gsclip(1)
c plot particles
      do 90 l = 1, nblok
c this segment is used for mpi computers
      nvp = nproc
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
      if (idproc.eq.0) then
c no special diagnostic node
         if (joff.eq.1) then
            nd = npp(l)
            do 50 j = 1, nd
            f(j) = part(ixp,j,l)
            f(nd+j) = part(iyp,j,l)
   50       continue
c special diagnostic node present
         else
            nvp = lvp - nproc
            id = 1
            call MPI_RECV(f,2*npmax,mreal,id,103,lworld,istatus,ierr)
c determine how many particles to plot
            call MPI_GET_COUNT(istatus,mreal,nd,ierr)
            nd = nd/2
         endif
         do 60 i = 1, nvp
c processor 0 plots his (or her) own data
         if (i.gt.1) then
c collects data from remaining nodes
            id = i - joff
            call MPI_RECV(f,2*npmax,mreal,id,103,lworld,istatus,ierr)
c determine how many particles to plot
            call MPI_GET_COUNT(istatus,mreal,nd,ierr)
            nd = nd/2
         endif
         if (nd.eq.0) go to 60
         koff = 0
         icol = i + 1 - ncols*(i/ncols)
         icol = kprime(icol+1)
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gspmci(icol)
c treat dots by drawing a line to itself
c        call spdots(f(1),f(nd+1),nd,icol,nd)
c treat dots by drawing a line to itself with non-zero width in x
         dd = (xmax - xmin)*float(npl)*2.0e-3
         call spddots(f(1),f(nd+1),dd,nd,icol,nd)
c update workstation, perform
         call guwk(idwk,1)
   60    continue
c other nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         id = idproc + joff - 1
         nd = npp(l)
         do 70 j = 1, nd
         f(j) = part(ixp,j,l)
         f(nd+j) = part(iyp,j,l)
   70    continue
         call MPI_SEND(f,2*nd,mreal,0,103,lworld,ierr)
c plot particles
         if (nd.eq.0) go to 90
         koff = 0
         icol = id + 2 - ncols*((id+1)/ncols)
         icol = kprime(icol+1)
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gspmci(icol)
c treat dots by drawing a line to itself
c        call spdots(f(1),f(nd+1),nd,icol,nd)
c treat dots by drawing a line to itself with non-zero width in x
         dd = (xmax - xmin)*float(npl)*2.0e-3
         call spddots(f(1),f(nd+1),dd,nd,icol,nd)
c update workstation, perform
         call guwk(idwk,1)
      endif
   90 continue
c update plot number
      iplot = iplot + 1
      if (iplot.eq.nplot) iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         if (idproc.eq.0) call readrc(irc)
c this segment is used for mpi computers
         msg(1) = irc
         msg(2) = nplot
         call PBICAST(msg,2)
         irc = msg(1)
         nplot = msg(2)
      endif
c reset defaults
      iclr = 0
      iupd = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine dsdmpln(orgx,smin,tmin,tmax,dd,ngs,nlts,nmks,mks,chrs,c
     1hh)
c this subroutine displays line or marker samples, with short character
c labels placed underneath
c dots are displayed with with non-zero width in x
c orgx = x origin of window
c smin = minimum x value of plotting window
c tmin/tmax = range of y values of plotting window
c dd = smallest visible size, such as (xmax - xmin)*4.0e-3
c ngs = number of subarrays being plotted
c nlts = number of line types available
c nmks = number of markers available
c mks = flag to determine whether lines or markers are used
c chrs = array of ngs short character labels for line or marker samples,
c each should be less than or equal to 10 characters in length
c chh = character height
      character*(*) chrs(ngs)
c rx, ry = ndc coordinates of upper-right corner of workstation window
c ncols = number of foreground colors available for line plotting
c kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c x,y = scratch variables for plotting
      dimension x(4), y(2)
c omit samples if there is only one curve
      if (ngs.le.1) return
      mkr = abs(mks)
c set marker size scale factor, 1.0 = nominal
      if (mkr.ne.0) call gsmksc(1.0)
c draw line or marker samples
      stx = .4*chh*(rx/ry)
      x(1) = orgx + 4.*stx
      x(2) = smin - 4.*stx
      dx = (x(2) - x(1))/3.
      x(3) = x(1) + dx
      x(4) = x(2) - dx
      ax = orgx + 2.*stx
      dy = (tmax - tmin - 4.*chh)/real(ngs + 1)
      ay = tmax - chh
c draw samples, first cycle through colors, then line or marker types
      do 30 k = 1, ngs
      icol = k + 1 - ncols*(k/ncols)
      icol = kprime(icol+1)
      y(1) = ay - dy*real(k)
      y(2) = y(1)
      if ((mks.eq.0).or.((mks.lt.0).and.(k.eq.1))) then
c set polyline color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gsplci(icol)
         ltype = 1 + (k - 1)/ncols - nlts*((k - 1)/(nlts*ncols))
c set linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
         call gsln(ltype)
c draw polyline
         call gpl(2,x,y)
      else
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gspmci(icol)
         mtype = mkr + (k - 1)/ncols - nmks*((mkr - 1 + (k - 1)/ncols)/n
     1mks)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mtype)
c dots
         if (mtype.eq.1) then
c treat dots by drawing a line to itself
c           call spdots(x(3),y,2,icol,3)
c treat dots by drawing a line to itself with non-zero width in x
            call spddots(x(3),y,dd,2,icol,3)
         else
c draw polymarker
            call gpm(2,x(3),y)
         endif
      endif
c draw label underneath sample line or marker
      y(1) = y(1) - 2.*chh
c draw text
      call gtx(ax,y(1),chrs(k))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine spddots(x,y,dd,npt,icol,nxbs)
c this subroutine draws dot markers by drawing a line to itself
c with non-zero width in x
c x, y = arrays to be plotted
c dd = smallest visible size, such as (xmax - xmin)*4.0e-3
c npt = number of points to be plotted
c icol = color index
c nxbs = dimension of x, y arrays
      dimension x(nxbs), y(nxbs)
c xs, ys = scratch arrays for plotting
      dimension xs(2), ys(2)
c set polyline color index
      call gsplci(icol)
      do 10 j = 1, npt
      xs(1) = x(j)
      ys(1) = y(j)
      xs(2) = xs(1) + dd
      ys(2) = ys(1)
c draw polyline
      call gpl(2,xs,ys)
   10 continue
      return
      end
