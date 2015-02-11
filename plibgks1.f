c general 1d parallel gks graphics library
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: october 30, 2008
c-----------------------------------------------------------------------
      subroutine PGRCLOSE1
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
      subroutine PDISPR(f,g,nvp,label,xmin,xmax,isc,ist,mks,nx,nxv,nblok
     1,ngs,chr,chrs,irc)
c this subroutine plots ngs subarrays of the array f, on a common graph,
c each plot with nx points, versus a linear function in x,
c where xmin < x < xmax, for distributed data.
c depending on the number of colors in the display device, each subarray
c is plotted with a different color, given in order by:
c blue, red, yellow, cyan, magenta, green, and foreground.
c after all the colors are cycled through, then different line styles
c are cycled through if mks=0, in order: solid, dash, dot, dash-dot,
c or different marker types if mks>0: dot, plus, star, circle, cross.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = distributed field array to be plotted
c g = scratch array for receiving messages
c nvp = number of real or virtual processors requested
c label = long character string label for plot
c xmin/xmax = range of x values in plot
c isc = power of 2 scale of y coordinate for plot
c ist = flag for choosing positive and/or negative values
c the plots have a common scale in y given by ymax and ymin.
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist = 1, then ymax = 2**isc and ymin = 0.
c if ist = -1, then ymax = 0 and ymin = -2**isc.
c if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plots, determined by the absolute value of f
c mks = flag to determine whether lines or markers are used,
c mks=0 means cycle through lines styles, mks > 0 means cycle through
c marker styles, using the value of mks as the initial marker,
c mks < 0 means draw the first subarray with a line, then subsequent
c subarrays with markers, using the value of abs(mks) as the initial
c marker.
c nx = number of points plotted in each subarray
c nxv = first dimension of field array f, must be >= nx/nvp
c nblok = number of particle partitions
c ngs = second dimension of array f, number of subarrays to be plotted
c chr = additional long character string comment for plot
c chrs = array of ngs short character labels used by subroutine tickd
c to label individual line or marker samples
c irc = return code (0 = normal return)
      implicit none
      real f, g
      integer nvp, isc, ist, mks, nx, nxv, nblok, ngs, irc
      real xmin, xmax
      character*(*) label, chr
      dimension f(nxv,ngs,nblok), g(nxv,ngs)
      character*(*) chrs(ngs)
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
      integer istatus, msg
      integer nlts, nmks, ntx, nty, nxp, nxp1, i, j, k, l, is, it, ir
      integer npl1, npl, iy, ix, nrt, idproc, joff, id, lvp, nl, ierr
      real gmax, hmax
      real dv, smin, tmin, smax, tmax, csize, algdvi, fmax, fmin
      real rmax, rmin, ymin, ymax, dx, apl, aplx, aply, orx, ory
      real smn, smx, tmn, tmx, chh, xmn
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
c nlts = number of line types available
c nmks = number of markers available
c ntx/nty = number of ticks in grid in x/y direction
      data nlts,nmks,ntx,nty /4,5,11,9/
c csize = vertical size of characters
      data csize /0.034/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
      nxp = (nx - 1)/nvp + 1
      nxp1 = nxp + 1
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
         do 20 k = 1, ngs
         do 10 j = 1, nxp
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
      ymin = range(1)
      ymax = range(2)
c clip range if necessary
      if (ist.eq.1) then
         ymin = 0.
      else if (ist.eq.(-1)) then
         ymax = 0.  
      endif
c parameters for plots
      dx = xmax - xmin
      if (nx.gt.1) dx = dx/float(nx)
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
c draw grid and labels, call identity transformation
      call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample lines or markers
      call dsmpln(orx,smn,tmn,tmx,ngs,nlts,nmks,mks,chrs,chh)
c define transformation number 2
      nrt = 1
c set window
      call gswn(nrt,xmin,xmax,ymin,ymax)
c set viewport
      call gsvp(nrt,smn,smx,tmn,tmx)
c select normalization transformation
      call gselnt(nrt)
c set linewidth scale factor, 1.0 = nominal
      call gslwsc(1.0)
c use markers
      if (mks.ne.0) then
c set marker size scale factor, 1.0 = nominal
         call gsmksc(1.0)
      endif
c set clipping indicator, 1 = on
      call gsclip(1)
      do 50 l = 1, nblok
c this segment is used for shared memory computers
c     do 70 k = 1, ngs
c     do 60 j = 1, nxv
c     g(j,k) = f(j,k,l)
c  60 continue
c  70 continue
c plot curves
c     call mplotit(g,xmin,dx,nxv,nblok,ngs,nxp1,mks,nlts,nmks)
c this segment is used for mpi computers
      if (lvp.eq.nproc) then
         joff = 1
      else
         joff = 0
      endif
      nl = nxv*(ngs - 1) + nxp1
      if (idproc.eq.0) then
c special diagnostic node present
         if (joff.eq.0) then
            id = 1
            call MPI_RECV(f,nl,mreal,id,102,lworld,istatus,ierr)
         endif
c processor 0 plots his (or her) own data
         call mplotit(f,xmin,dx,nxv,nblok,ngs,nxp1,mks,nlts,nmks)
c then collects data from remaining nodes
         do 40 i = 2, nvp
         id = i - joff
         call MPI_RECV(g,nl,mreal,id,102,lworld,istatus,ierr)
c plot curves
         xmn = xmin + dx*float(nxp*(i - 1))
         call mplotit(g,xmn,dx,nxv,nblok,ngs,nxp1,mks,nlts,nmks)
   40    continue
c some nodes send data to node 0
      else if (idproc.le.(nvp-joff)) then
         call MPI_SEND(f,nl,mreal,0,102,lworld,ierr)
c and plot the local data
         id = idproc + joff - 1
         xmn = xmin + dx*float(nxp*id)
         call mplotit(f,xmn,dx,nxv,nblok,ngs,nxp1,mks,nlts,nmks)
c other nodes just plot the local data
      else if (idproc.le.nproc) then
         id = min(idproc+joff,nvp) - 1
         xmn = xmin + dx*float(nxp*id)
         call mplotit(f,xmn,dx,nxv,nblok,ngs,nxp1,mks,nlts,nmks)
      endif
   50 continue
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
      subroutine PDISPS(f,g,nvp,label,xmin,xmax,isc,ist,nx,nxv,nblok,chr
     1,irc)
c this subroutine plots an array f versus a linear function in x,
c where xmin < x < xmax.  It is plotted in solid line style, in blue
c if color is available, for distributed data.
c multiple plots per page can be displayed by dividing the screen into
c n x n subregions, where n*n is the next largest integer >= nplot
c the location (ix,iy) of a plot in the subregions is determined by
c the parameter iplot = ix + iy*n
c f = distributed field array to be plotted
c g = scratch array for receiving messages
c nvp = number of real or virtual processors requested
c label = long character string label for plot
c xmin/xmax = range of x values in plot
c isc = power of 2 scale of y coordinate for plot
c ist = flag for choosing positive and/or negative values
c the plot has a scale in y given by ymax and ymin.
c if ist = 0, then ymax = 2**isc and ymin = -2**isc.
c if ist = 1, then ymax = 2**isc and ymin = 0.
c if ist = -1, then ymax = 0 and ymin = -2**isc.
c if ist = 2, then ymin = fmin, ymax = fmin + 2**ir,
c where fmin/fmax are the function minimum/maximum, 
c and ir = power of 2 scale for (fmax - fmin)
c if abs(isc) < 116, then the isc value passed is used for scale.
c if abs(isc) > 116, then the program finds the minimum value of isc
c which will contain the plot, determined by the absolute value of f
c nx = number of points plotted
c nxv = first dimension of field array f, must be >= nx/nvp
c nblok = number of particle partitions
c chr = additional long character string comment for plot
c irc = return code (0 = normal return)
      implicit none
      real f, g
      integer nvp, isc, ist, nx, nxv, nblok, irc
      real xmin, xmax
      character*(*) label, chr
      dimension f(nxv,nblok), g(nxv)
      character*1 chs(1)
      data chs /' '/
      call PDISPR(f,g,nvp,label,xmin,xmax,isc,ist,0,nx,nxv,nblok,1,chr,c
     1hs,irc)
      return
      end
c-----------------------------------------------------------------------
      subroutine mplotit(f,xmin,dx,nxv,nblok,ngs,nxp1,mks,nlts,nmks)
c this subroutine plots curves, first cycle through colors, then line or
c marker types
c f = distributed field array to be plotted
c xmin/dx = range/increment of x values in plot
c ngs = second dimension of array f, number of subarrays to be plotted
c nxp1 = number of points per processor plotted in each subarray
c mks = flag to determine whether lines or markers are used,
c mks=0 means cycle through lines styles, mks > 0 means cycle through
c marker styles, using the value of mks as the initial marker,
c mks < 0 means draw the first subarray with a line, then subsequent
c subarrays with markers, using the value of abs(mks) as the initial
c marker.
c nlts = number of line types available
c nmks = number of markers available
      implicit none
      integer nxv, nblok, ngs, nxp1, mks, nlts, nmks
      real xmin, dx
      real f
      dimension f(nxv,ngs,nblok)
c data in common blocks
      integer idwk, ncols, iplot, nplot, iclr, iupd, idstr, idloc, nclsp
      integer ifrg, isx, isy, kprime
      real rx, ry
c ncols = number of foreground colors available for line plotting
c kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c local data
      integer i, j, k, l, nxs, mkr, icol, nxb, npts, nptl, ltype, mtype
      integer npt, js
c nxbs = length of scratch variable for plotting
      integer nxbs
      real x, y
      parameter(nxbs=65)
      dimension x(nxbs), y(nxbs)
      nxs = nxbs - 1
      mkr = abs(mks)
c plot curves, first cycle through colors, then line or marker types
      do 70 k = 1, ngs
      icol = k + 1 - ncols*(k/ncols)
      icol = kprime(icol+1)
c use line types
      if ((mks.eq.0).or.((mks.lt.0).and.(k.eq.1))) then
c blocking parameters for plots
         nxb = (nxp1 - 2)/nxs + 1
         npts = nxbs
c length of last block
         nptl = nxp1 - nxs*(nxb - 1)
c set polyline color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gsplci(icol)
         ltype = 1 + (k - 1)/ncols - nlts*((k - 1)/(nlts*ncols))
c set linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
         call gsln(ltype)
c plot curve
         npt = npts
         do 30 l = 1, nblok
c loop over number of blocks
         do 20 j = 1, nxb
         js = nxs*(j - 1)
         if (j.eq.nxb) npt = nptl
c calculate x,y axes for block
         do 10 i = 1, npt
         x(i) = xmin + dx*float(i + js - 1)
         y(i) = f(i+js,k,l)
   10    continue
c draw polyline
         call gpl(npt,x,y)
   20    continue
   30    continue
c use markers
      else
c blocking parameters for plots
         nxb = (nxp1 - 1)/nxs + 1
         npts = nxs
c length of last block
         nptl = nxp1 - nxs*(nxb - 1)
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
         call gspmci(icol)
         mtype = mkr + (k - 1)/ncols - nmks*((mkr - 1 + (k - 1)/ncols)/n
     1mks)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
         call gsmk(mtype)
c plot polymarkers
         npt = npts
         do 60 l = 1, nblok
c loop over number of blocks
         do 50 j = 1, nxb
         js = nxs*(j - 1)
         if (j.eq.nxb) npt = nptl
c calculate x,y axes for block
         do 40 i = 1, npt
         x(i) = xmin + dx*float(i + js - 1)
         y(i) = f(i+js,k,l)
   40    continue
c dots
         if (mtype.eq.1) then
c treat dots by drawing a line to itself
            call spdots(x,y,npt,icol,nxbs)
         else
c draw polymarker
            call gpm(npt,x,y)
         endif
   50    continue
   60    continue
      endif
   70 continue
      return
      end
