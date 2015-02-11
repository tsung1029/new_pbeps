      subroutine gitwks(idcon,iwtype)
c this is a site-dependent subroutine which returns
c connection identifier and workstation type
c version for null library
c iwtype = workstation type
c 0 = null device
      iwtype = 0
c idcon = connection identifier, 1 seems to work
      idcon = 1
      return
      end
c gks device driver for null output
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: november 3, 1997
      subroutine gqops(istat)
c inquire operating state value
c input arguments: none
c istat = operating state (0=gks closed,1=gks open,2=workstation open,
c 3=workstation active,4=segment open)
      istat = 0
      return
      end
      subroutine gopks(nerrfl,meml)
c open gks
c input arguments: all
c nerrfl = error file unit number, 6 for terminal
c meml = storage limit for one segment, in bytes
      return
      end
      subroutine gopwk(idwk,idcon,iwtype)
c open workstation
c input arguments: all
c idwk = workstation identifier
c idcon = connection identifier
c iwtype = workstation type
      return
      end
      subroutine gqopwk(n,ierr,now,idwk)
c inquire set number of open workstations
c input arguments: n
c n = set member requested
c ierr = error indicator
c now = number of open workstations
c idwk = workstation identifier of the nth member
      now = 0
      idwk = 0
      ierr = 0
      return
      end
      subroutine gqwkca(iwtype,ierr,iwkca)
c inquire workstation category
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator
c iwkca = workstation category
c (0 = output, 1 = input, 2 = outin, 3 = wiss, 4 = mo, 5 = mi)
      iwkca = 0
      ierr = 0
      return
      end
      subroutine gscnid(idcon,connam)
c set connection identifier
c input arguments: all
c idcon = connection identifier
c connam = connection name
      character*8 connam
      return
      end
      subroutine gqwkc(idwk,ierr,idcon,iwtype)
c inquire workstation connection and type
c input arguments: idwk
c idwk = workstation identifier
c ierr = error indicator
c idcon = connection identifier
c iwtype = workstation type
      idcon = 0
      iwtype = 0
      ierr = 0
      return
      end
      subroutine gacwk(idwk)
c activate workstation
c input arguments: all
c idwk = workstation identifier
      return
      end
      subroutine gqacwk(n,ierr,naw,idwk)
c inquire set of active workstations
c input arguments: n
c n = set member requested
c ierr = error indicator
c naw = number of active workstations
c idwk = workstation identifier of the nth member
      naw = 0
      idwk = 0
      ierr = 0
      return
      end
      subroutine gqwks(idwk,ierr,istat)
c inquire workstation state
c input arguments: idwk
c idwk = workstation identifier
c ierr = error indicator
c istat = workstation state (0 = inactive, 1 = active)
      istat = 0
      ierr = 0
      return
      end
      subroutine gqcf(iwtype,ierr,ncoli,iscol,npci)
c inquire color facilities
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator (0=inquiry successful)
c ncoli = number of colors available
c iscol = color availability indicator, 0 = monochrome, 1 = color
c npci = number of predefined color indices
      ncoli = 2
      iscol = 0
      npci = 2
      ierr = 0
      return
      end
      subroutine gscr(idwk,ic,cr,cg,cb)
c set color representation
c input arguments: all
c idwk = workstation identifier
c ic = color index
c cr/cg/cb = red/green/blue component (0 < cr,cg,cb < 1)
      return
      end
      subroutine gqeci(idwk,n,ierr,ncoli,icol)
c inquire list element of color indices
c input arguments: idwk, n
c idwk = workstation identifier
c n = requested list element
c ierr = error indicator
c ncoli = number of indices currently defined
c icol = color index of requested element
      ncoli = 2
      icol = 0
      ierr = 0
      return
      end
      subroutine gqdsp(iwtype,ierr,idcun,dcx,dcy,lx,ly)
c inquire maximum display surface size
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator (0=inquiry successful)
c idcun = device coordinate units, (0=meters,1=other)
c dcx/dcy = width/height in device coordinate units
c lx/ly = width/height in device (raster) units
      idcun = 1
      dcx = 1.0
      dcy = 1.0
      lx = 1
      ly = 1
      ierr = 0
      return
      end
      subroutine gswkwn(idwk,xmin,xmax,ymin,ymax)
c set workstation window
c input arguments: all
c idwk = workstation identifier
c xmin/xmax = window x coordinates in ndc
c ymin/ymax = window y coordinates in ndc
      return
      end
      subroutine gswkvp(idwk,xmin,xmax,ymin,ymax)
c set workstation viewport
c input arguments: all
c idwk = workstation identifier
c xmin/xmax = viewport x coordinates in device coordinates
c ymin/ymax = viewport y coordinates in device coordinates
      return
      end
      subroutine gqsts(idwk,idstr,mldr,ierr,mode,iesw,lstr,str,ipet,eare
     1a,lenb,ipos,ldr,datar)
c inquire string device state
c input arguments: idwk, idstr, mldr
c idwk = workstation identifier
c idstr = string device number
c mldr = maximum dimension of data record array
c ierr = error indicator (0=inquiry successful)
c mode = operating mode (0=request,1=sample,2=event)
c iesw = echo switch (0=no echo,1=echo)
c lstr = number of characters in initial string
c str = initial string
c ipet = prompt/echo type (1=normal)
c earea(1)/earea(2) = echo area x coordinates in device coordinates
c earea(3)/earea(4) = echo area y coordinates in device coordinates
c lenb = input buffer size
c ipos = inital edit position
c ldr = length of data record array
c datar = data record array
      character*(*) str
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 0
      lstr = 1
      ipet = 1
      earea(1) = 0.
      earea(2) = 1.
      earea(3) = 0.
      earea(4) = 1.
      lenb = 1
      ipos = 1
      ldr = 1
      ierr = 140
      return
      end
      subroutine gqlcs(idwk,idloc,itype,mldr,ierr,mode,iesw,nrt,pxi,pyi,
     1ipet,earea,ldr,datar)
c inquire locator device state
c input arguments: idwk, idloc, itype, mldr
c idwk = workstation identifier
c idloc = locator device number
c itype = type of returned value (0=set,1=realized)
c mldr = maximum dimension of data record array
c ierr = error indicator (0=inquiry successful)
c mode = operating mode (0=request,1=sample,2=event)
c iesw = echo switch (0=no echo,1=echo)
c nrt = normalization transformation number for initial position
c pxi/pyi = initial position, in world coordinates
c ipet = prompt/echo type (1=normal)
c earea(1)/earea(2) = echo area x coordinates in device coordinates
c earea(3)/earea(4) = echo area y coordinates in device coordinates
c ldr = length of data record array
c datar = data record array
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 0
      nrt = 0
      pxi = .5
      pyi = .5
      ipet = 1
      earea(1) = 0.
      earea(2) = 1.
      earea(3) = 0.
      earea(4) = 1.
      ldr = 1
      ierr = 140
      return
      end
      subroutine ginst(idwk,idstr,lstr,str,ipet,xmin,xmax,ymin,ymax,lenb
     1,ipos,ldr,datar)
c initialize string device
c input arguments: all
c idwk = workstation identifier
c idstr = string device number
c lstr = number of characters in initial string
c str = initial string
c ipet = prompt/echo type (1=normal)
c xmin/xmax = echo area x coordinates in device coordinates
c ymin/ymax = echo area y coordinates in device coordinates
c lenb = input buffer size
c ipos = inital edit position
c ldr = length of data record array
c datar = data record array
      return
      end
      subroutine ginlc(idwk,idloc,nrt,px,py,ipet,xmin,xmax,ymin,ymax,ldr
     1,datar)
c initialize locator device
c input arguments: all
c idwk = workstation identifier
c idloc = locator device number
c nrt = initial normalization transformation number
c px/py = initial locator position, in world coordinates
c ipet = prompt/echo type (1=normal)
c xmin/xmax = echo area x coordinates in device coordinates
c ymin/ymax = echo area y coordinates in device coordinates
c ldr = length of data record array
c datar = data record array
      return
      end
      subroutine gdawk(idwk)
c deactivate workstation
c input arguments: all
c idwk = workstation identifier
      return
      end
      subroutine gclwk(idwk)
c close workstation
c input arguments: all
c idwk = workstation identifier
      return
      end
      subroutine gclks
c close gks
      return
      end
      subroutine gclrwk(idwk,icofl)
c clear workstation
c input arguments: all
c idwk = workstation identifier
c icofl = control flag (0=conditionally,1=always)
      return
      end
      subroutine gqcntn(ierr,nrt)
c inquire current normalization transformation number
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c nrt = transformation number (0 <= nrt <= 25)
      nrt = 0
      ierr = 0
      return
      end
      subroutine gqnt(nrt,ierr,window,viewpt)
c inquire normalization transformation
c input arguments: nrt
c nrt = transformation number (0 <= nrt <= 25)
c ierr = error indicator (0=inquiry successful)
c window(1)/window(2) = window x coordinates in world coordinates
c window(3)/window(4) = window y coordinates in world coordinates
c viewpt(1)/viewpt(2) = viewport x coordinates in ndc
c viewpt(3)/viewpt(4) = viewport y coordinates in ndc
      dimension window(4), viewpt(4)
      window(1) = 0.
      window(2) = 1.
      window(3) = 0.
      window(4) = 1.
      viewpt(1) = 0.
      viewpt(2) = 1.
      viewpt(3) = 0.
      viewpt(4) = 1.
      ierr = 0
      return
      end
      subroutine gswn(nrt,xmin,xmax,ymin,ymax)
c set window
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
c xmin/xmax = window x coordinates in world coordinates
c ymin/ymax = window y coordinates in world coordinates
      return
      end
      subroutine gsvp(nrt,xmin,xmax,ymin,ymax)
c set viewport
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
c xmin/xmax = viewport x coordinates in ndc
c ymin/ymax = viewport y coordinates in ndc
      return
      end
      subroutine gselnt(nrt)
c select normalization transformation
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
      return
      end
      subroutine gqln(ierr,ltype)
c inquire linetype
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      ltype = 1
      ierr = 0
      return
      end
      subroutine gslwsc(alwsc)
c set linewidth scale factor
c input arguments: all
c alwsc = linewidth scale factor, (alwsc > 1.0)
      return
      end
      subroutine gsplci(icol)
c set polyline color index
c input arguments: all
c icol = color index
      return
      end
      subroutine gsln(ltype)
c set linetype
c input arguments: all
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      return
      end
      subroutine gpl(n,px,py)
c draw polyline
c input arguments: all
c n = number of points to be connected by a polyline
c px/py = x/y coordinates of points in world coordinates
      dimension px(n), py(n)
      return
      end
      subroutine gqmk(ierr,mtype)
c inquire marker type
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c mtype = marker type
c 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross, 6 = square,
c 7 = square with cross, 8 = diamond, 9 = diamond with cross,
c 10 = filled circle
      mtype = 1
      ierr = 0
      return
      end
      subroutine gsmksc(amksc)
c set marker size scale factor
c input arguments: all
c amksc = linewidth scale factor, (amksc > 1.0)
      return
      end
      subroutine gspmci(imcol)
c set polymarker color index
c input arguments: all
c imcol = polymarker color index
      return
      end
      subroutine gsmk(mtype)
c set marker type
c input arguments: all
c mtype = marker type
c 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross, 6 = square,
c 7 = square with cross, 8 = diamond, 9 = diamond with cross,
c 10 = filled circle
      return
      end
      subroutine gpm(n,px,py)
c draw polymarker
c input arguments: all
c n = number of points to be connected by a polymarker
c px/py = x/y coordinates of points in world coordinates
      dimension px(n), py(n)
      return
      end
      subroutine gqtxal(ierr,itxalh,itxalv)
c inquire text alignment
c input arguments: none
c ierr = error indicator
c itxalh = horizontal component
c 0 = normal, 1 = left, 2 = center, 3 = right
c itxalv = vertical component:
c = 0 normal, 1 = top, 2 = cap, 3 = half, 4 = base, 5 = bottom
      itxalh = 0
      itxalv = 0
      ierr = 0
      return
      end
      subroutine gqchh(ierr,chh)
c inquire character height
c input arguments: none
c ierr = error indicator
c chh = character height, in world coordinates
      chh = .01
      ierr = 0
      return
      end
      subroutine gqtxp(ierr,itxp)
c inquire text path
c input arguments: none
c ierr = error indicator
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      itxp = 0
      ierr = 0
      return
      end
      subroutine gqchup(ierr,chux,chuy)
c inquire character up vector
c input arguments: none
c ierr = error indicator
c chux/chuy = up vector, in world coordinates
      chuy = 0.
      chux = 1.0
      ierr = 0
      return
      end
      subroutine gstxal(itxalh,itxalv)
c set text alignment
c input arguments: all
c itxalh = horizontal component:
c 0 = normal, 1 = left, 2 = center, 3 = right
c itxalv = vertical component:
c = 0 normal, 1 = top, 2 = cap, 3 = half, 4 = base, 5 = bottom
      return
      end
      subroutine gstxfp(nfont,iprec)
c set text font
c input arguments: all
c nfont = character font number
c iprec = text precision (0=string,1=character,2=stroke)
      return
      end
      subroutine gstxp(itxp)
c set text path
c input arguments: all
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      return
      end
      subroutine gstxci(itcol)
c set text color index
c input arguments: all
c itcol = text color index
      return
      end
      subroutine gschh(chh)
c set character height
c input arguments: all
c chh = character height, in world coordinates
      return
      end
      subroutine gschup(chux,chuy)
c set character up vector
c input arguments: all
c chux/chuy = up vector, in world coordinates
      return
      end
      subroutine gschxp(chxp)
c set character expansion factor
c input arguments: all
c chxp = character expansion factor (>0.)
      return
      end
      subroutine gtx(px,py,chars)
c display text
c input arguments: all
c px/py = starting x/y position of text, in world coordinates
c chars = test string to be displayed
      character*(*) chars
      return
      end
      subroutine gsfaci(ifcol)
c set fill area color index
c input arguments: all
c ifcol = color index
      return
      end
      subroutine gsfais(ints)
c set fill area interior style
c input arguments: all
c ints = desired interior style:
c 0 = hollow (default), 1 = solid, 2 = pattern, 3 = hatch
      return
      end
      subroutine gfa(n,px,py)
c fill area
c input arguments: all
c n = number of points in fill area
c px,py = arrays of points, in world coordinates
      dimension px(n), py(n)
      return
      end
      subroutine gca(px,py,qx,qy,icxd,icyd,ncs,nrs,idx,idy,icola)
c cell array
c input arguments: all
c px,py = lower-left cell corner, in world coordinates
c qx,qy = upper-right cell corner, in world coordinates
c icxd,icyd = color index array dimensions
c ncs,nrs = starting column and row in the color index array
c idx,idy = number of columns and rows in the cell array
c icola = color index array
      dimension icola(icxd,icyd)
      return
      end
      subroutine gqclip(ierr,indcl,clrect)
c inquire clipping indicator
c input arguments: none
c ierr = error indicator
c indcl = clipping indicator (0=no clip, 1=clip)
c clrect = clipping rectangle, in ndc
      dimension clrect(4)
      indcl = 1
      clrect(1) = 0.
      clrect(2) = 1.
      clrect(3) = 0.
      clrect(4) = 1.
      ierr = 0
      return
      end
      subroutine gsclip(iclsw)
c set clipping indicator
c input arguments: all
c iclsw = clipping switch (0=no clip,1=clip)
      return
      end
      subroutine guwk(idwk,iregfl)
c update workstation
c input arguments: all
c idwk = workstation identifier
c iregfl = regeneration flag (0=postponed,1=perform)
      return
      end
      subroutine grqst(idwk,idstr,istat,lostr,str)
c request string
c input arguments: idwk, idstr
c idwk = workstation identifier
c idstr = string device number
c istat = return status (0=none,1=ok)
c lostr = number of characters in string
c str = returned string
      character*(*) str
      istat = 0
      lostr = 0
      return
      end
      subroutine grqlc(idwk,idloc,istat,nrt,px,py)
c request locator
c input arguments: idwk, idloc
c idwk = workstation identifier
c idloc = locator device number
c istat = return status (0=none,1=ok)
c nrt = normalization transformation number
c px/py = position, in world coordinates
      istat = 0
      nrt = 0
      return
      end
      subroutine gsstm(idwk,idstr,mode,iesw)
c set string mode
c input arguments: all
c idwk = workstation identifier
c idstr = string device number
c mode  =  mode of operation (0 = request,1 = sample,2 = event)
c iesw  =  echo switch (0 = no echo,1 = echo)
      return
      end
      subroutine gslcm(idwk,idloc,mode,iesw)
c set locator mode
c input arguments: all
c idwk = workstation identifier
c idloc = locator device number
c mode  =  mode of operation (0 = request,1 = sample,2 = event)
c iesw  =  echo switch (0 = no echo,1 = echo)
      return
      end
      subroutine gwait(tout,idwk,icl,idnr)
c await event
c input arguments: tout
c tout = timeout period, seconds
c idwk = workstation identifier
c icl = input class:
c (0=no class,1=locator,2=stroke,3=valuator,4=choice,5=pick,6=string)
c idnr = device number
      idwk = 0
      icl = 0
      idnr = 0
      return
      end
      subroutine ggtst(lostr,str)
c get string
c input arguments: none
c lostr = number of characters in string
c str = input string
      character*(*) str
      lostr = 0
      return
      end
      subroutine ggtlc(nrt,px,py)
c get locator
c input arguments: none
c nrt = normalization transformation number
c px/py = position, in world coordinates
      nrt = 0
      return
      end
      subroutine gesc(idfct,ldi,datai,mldr,ldr,datar)
c escape function
c input arguments: idfct, ldi, datai, mldr
c idfct = escape function identifier
c ldi = length of input data record
c datai = input data record
c mldr = maximum length of output data record
c ldr = length of output data record
c datar = output data record
      character*80 datai(ldi), datar(mldr)
      ldr = 0
c escape functions not supported
      ierr = 180
      return
      end
      integer function kywait()
c special function to request keystroke
c returns: ascii code for keystroke
      kywait = 0
      return
      end
