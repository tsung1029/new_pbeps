c Null general 2d parallel gks graphics library
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: april 9, 2003
c-----------------------------------------------------------------------
      subroutine PGRCLOSE1
c this subroutine deactivates workstation and closes gks
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
      irc = 0
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
      irc = 0
      return
      end
