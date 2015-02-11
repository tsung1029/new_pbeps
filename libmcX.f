      subroutine gitwks(idcon,iwtype)
c this is a site-dependent subroutine which returns
c connection identifier and workstation type
c version for macintosh graphics
c iwtype = workstation type
c 1 = mac console
      iwtype = 1
c idcon = connection identifier, 1 seems to work
      idcon = 1
      return
      end
c gks device driver for macintosh graphics using CarbonLib
c written by viktor k. decyk, ucla
c copyright 1997, regents of the university of california
c version for apple macintosh PowerPC with Absoft compiler
c update: february 21, 2002
      block data
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      common /macevnts/ sevent, mevent, strevt
      integer*2 sevent(8), mevent(8)
      character*8 strevt
c sevent = last string event
c mevent = last mouse event
c strevt = last string emulation event
      save /gksmc/, /macevnts/
      data kosv /0/
      data sevent, mevent /8*0, 8*0/
      data strevt /'        '/
      end
      subroutine gqops(istat)
c inquire operating state value
c input arguments: none
c istat = operating state (0=gks closed,1=gks open,2=workstation open,
c 3=workstation active,4=segment open)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      istat = kosv
      return
      end
      subroutine gopks(nerrfl,meml)
c open gks
c input arguments: all
c nerrfl = error file unit number, 6 for terminal
c meml = storage limit for one segment, in bytes
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      kosv = 1
      return
      end
      subroutine gopwk(idwk,idcon,iwtype)
c open workstation
c input arguments: all
c idwk = workstation identifier
c idcon = connection identifier
c iwtype = workstation type
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
c idc = current connection identifier
c myls = y screen offset, used because origin is at top of screen
c minx/miny = location of upper left hand corner of window
      minx = 6
      miny = 40
c llx/lly = current window width/height
      llx = 720
      lly = 540
c define a safe window size
      call setpsz(minx,miny,llx,lly)
c return size of a window
      call gitsiz(lx,ly)
      myls = ly - 1
      kosv = 2
      iwk = idwk
      idc = idcon
      return
      end
      subroutine gqopwk(n,ierr,now,idwk)
c inquire set number of open workstations
c input arguments: n
c n = set member requested
c ierr = error indicator
c now = number of open workstations
c idwk = workstation identifier of the nth member
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
      now = 0
      ierr = 0
      if (kosv.ge.2) now = 1
      if ((n.lt.0).or.(n.gt.1)) ierr = -12
      if (n.eq.1) idwk = iwk
      return
      end
      subroutine gqwkca(iwtype,ierr,iwkca)
c inquire workstation category
c input arguments: iwtype
c iwtype = workstation type
c ierr = error indicator
c iwkca = workstation category
c (0 = output, 1 = input, 2 = outin, 3 = wiss, 4 = mo, 5 = mi)
      ierr = 0
      if (iwtype.eq.1) then
         iwkca = 2
      else
         ierr = 22
      endif
      return
      end
      subroutine gscnid(idcon,connam)
c set connection identifier
c input arguments: all
c idcon = connection identifier
c connam = connection name
      character*8 connam
c     open(unit=idcon,file=connam,form='formatted',status='unknown')
      return
      end
      subroutine gqwkc(idwk,ierr,idcon,iwtype)
c inquire workstation connection and type
c input arguments: idwk
c idwk = workstation identifier
c ierr = error indicator
c idcon = connection identifier
c iwtype = workstation type
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwk = current workstation identifier
c idc = current connection identifier
      ierr = 0
      idcon = idc
      if (idwk.eq.iwk) then
         iwtype = 1
      else
         ierr = 20
      endif
      return
      end
      subroutine gacwk(idwk)
c activate workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c create a window
      call cr8win('gkslib')
c create mac defaults
      call mcdflts
      kosv = 3
      return
      end
      subroutine gqacwk(n,ierr,naw,idwk)
c inquire set of active workstations
c input arguments: n
c n = set member requested
c ierr = error indicator
c naw = number of active workstations
c idwk = workstation identifier of the nth member
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
      naw = 0
      ierr = 0
      if (kosv.ge.3) naw = 1
      if ((n.lt.0).or.(n.gt.1)) ierr = -12
      if (n.eq.1) idwk = iwk
      return
      end
      subroutine gqwks(idwk,ierr,istat)
c inquire workstation state
c input arguments: idwk
c idwk = workstation identifier
c ierr = error indicator
c istat = workstation state (0 = inactive, 1 = active)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
      ierr = 0
      istat = 0
      if (kosv.lt.2) ierr = 25
      if (kosv.ge.3) istat = 1
      if (idwk.ne.iwk) then
         istat = 0
         ierr = 20
      endif
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
      ierr = 0
c get available colors
      call gitnc(ncoli)
      iscol = 0
      if (ncoli.gt.2) iscol = 1
      npci = ncoli
      return
      end
      subroutine gscr(idwk,ic,cr,cg,cb)
c set color representation
c input arguments: all
c idwk = workstation identifier
c ic = color index
c cr/cg/cb = red/green/blue component (0 < cr,cg,cb < 1)
c convert to integer*2
      icr = 65535.*cr
      if (icr.gt.32767) icr = icr - 65536
      icg = 65535.*cg
      if (icg.gt.32767) icg = icg - 65536
      icb = 65535.*cb
      if (icb.gt.32767) icb = icb - 65536
c set a color map entry to a specified RGB value
      call crctab(ic,icr,icg,icb)
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwk = current workstation identifier
      ierr = 0
c get available colors
      call gitnc(ncoli)
      if (idwk.ne.iwk) ierr = -12
      if ((n.lt.0).or.(n.gt.ncoli)) ierr = -12
      if ((n.gt.0).and.(n.le.ncoli)) icol = n - 1
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
c return size of a window
      call gitsiz(llx,lly)
c llx/lly = current window width/height
      lx = llx
      ly = lly
      dcx = float(lx - 1)
      dcy = float(ly - 1)
      ierr = 0
      return
      end
      subroutine gswkwn(idwk,xmin,xmax,ymin,ymax)
c set workstation window
c input arguments: all
c idwk = workstation identifier
c xmin/xmax = window x coordinates in ndc
c ymin/ymax = window y coordinates in ndc
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c wwnvp = workstation window and viewport
c trans = normalization transformations
c clip to maximum window size
      wwnvp(1) = amax1(xmin,0.)
      wwnvp(2) = amin1(xmax,1.)
      wwnvp(3) = amax1(ymin,0.)
      wwnvp(4) = amin1(ymax,1.)
c redefine viewport of transformation 0
      trans(1,1) = wwnvp(1)
      trans(2,1) = wwnvp(2)
      trans(3,1) = wwnvp(3)
      trans(4,1) = wwnvp(4)
      return
      end
      subroutine gswkvp(idwk,xmin,xmax,ymin,ymax)
c set workstation viewport
c input arguments: all
c idwk = workstation identifier
c xmin/xmax = viewport x coordinates in device coordinates
c ymin/ymax = viewport y coordinates in device coordinates
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c wwnvp = workstation window and viewport
c return size of a window
      call gitsiz(llx,lly)
c llx/lly = current window width/height
      dcx = float(llx - 1)
      dcy = float(lly - 1)
c clip to maximum screen size
      wwnvp(5) = amax1(xmin,0.)
      wwnvp(7) = amax1(ymin,0.)
      wwnvp(6) = amin1(xmax,dcx)
      wwnvp(8) = amin1(ymax,dcy)
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c cea = character echo area
      character*(*) str
      character*80 datar(mldr)
      dimension earea(4)
      mode = 0
      iesw = 1
      lstr = 1
      str(1:1) = ' '
      ipet = 1
      do 10 j = 1, 4
      earea(j) = cea(j)
   10 continue
      lenb = 80
      ipos = 1
      ldr = 1
      ierr = 0
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
      iesw = 1
      nrt = 0
      pxi = .5
      pyi = .5
      ipet = 1
c return size of a window
      call gitsiz(llx,lly)
c llx/lly = current window width/height
      earea(1) = 0.
      earea(2) = float(llx - 1)
      earea(3) = 0.
      earea(4) = float(lly - 1)
      ldr = 1
      ierr = 0
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c idst = current string device number
c cea = character echo area
      character*80 datar(ldr)
      character*(*) str
c save string device number
      idst = idstr
c save character echo area
      cea(1) = xmin
      cea(2) = xmax
      cea(3) = ymin
      cea(4) = ymax
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c idlc = current locator device number
      character*80 datar(ldr)
      integer GetCursor
      external GetCursor
      integer chandle, cursptr
c save locator device number
      idlc = idloc
c cross cursor
      if (ipet.eq.2) then
c get a handle to a specified 'CURS' resource
         chandle = GetCursor(val(2))
         cursptr = long(chandle)
c change the shape of the mouse cursor
         call SetCursor(val(cursptr))
c beam cursor
      elseif (ipet.eq.3) then
c get a handle to a specified 'CURS' resource
         chandle = GetCursor(val(1))
         cursptr = long(chandle)
c change the shape of the mouse cursor
         call SetCursor(val(cursptr))
c normal cursor
      else
c initialize a cursor to the standard arrow
         call InitCursor()
      endif
      return
      end
      subroutine gdawk(idwk)
c deactivate workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c close window
      call delwind
      kosv = 2
      return
      end
      subroutine gclwk(idwk)
c close workstation
c input arguments: all
c idwk = workstation identifier
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      kosv = 1
      return
      end
      subroutine gclks
c close gks
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
      kosv = 0
      return
      end
      subroutine gclrwk(idwk,icofl)
c clear workstation
c input arguments: all
c idwk = workstation identifier
c icofl = control flag (0=conditionally,1=always)
c erase screen
      call eraser
      return
      end
      subroutine gqcntn(ierr,nrt)
c inquire current normalization transformation number
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c nrt = transformation number (0 <= nrt <= 25)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
      ierr = 0
      nrt = nrtn
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c trans = normalization transformations
      dimension window(4), viewpt(4)
      ierr = 0
      if ((nrt.lt.0).or.(nrt.ge.maxt)) then
         nrt = 0
         ierr = 1
      endif
      n = nrt + 1
      do 10 j = 1, 4
      window(j) = trans(j+4,n)
      viewpt(j) = trans(j,n)
   10 continue
      return
      end
      subroutine gswn(nrt,xmin,xmax,ymin,ymax)
c set window
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
c xmin/xmax = window x coordinates in world coordinates
c ymin/ymax = window y coordinates in world coordinates
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c trans = normalization transformations
      if ((nrt.lt.0).or.(nrt.ge.maxt)) nrt = 0
      n = nrt + 1
c store transformation
      trans(5,n) = xmin
      trans(6,n) = xmax
      trans(7,n) = ymin
      trans(8,n) = ymax
      return
      end
      subroutine gsvp(nrt,xmin,xmax,ymin,ymax)
c set viewport
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
c xmin/xmax = viewport x coordinates in ndc
c ymin/ymax = viewport y coordinates in ndc
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c trans = normalization transformations
      if ((nrt.lt.0).or.(nrt.ge.maxt)) nrt = 0
      n = nrt + 1
c store transformation
      trans(1,n) = xmin
      trans(2,n) = xmax
      trans(3,n) = ymin
      trans(4,n) = ymax
      return
      end
      subroutine gselnt(nrt)
c select normalization transformation
c input arguments: all
c nrt = transformation number (1 <= nrt <= 25)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c kcf = clipping flag
c idvt = device transformation
c wwnvp = workstation window and viewport
c trans = normalization transformations
      integer*2 rect(4)
      if ((nrt.lt.0).or.(nrt.ge.maxt)) nrt = 0
      nrtn = nrt
      n = nrt + 1
c clip to workstation window
      xmin = amax1(wwnvp(1),trans(1,n))
      xmax = amin1(wwnvp(2),trans(2,n))
      ymin = amax1(wwnvp(3),trans(3,n))
      ymax = amin1(wwnvp(4),trans(4,n))
c convert from viewport to screen coordinates
      scx = (wwnvp(6) - wwnvp(5))/(wwnvp(2) - wwnvp(1))
      idvt(1) = (xmin - wwnvp(1))*scx + wwnvp(5) + .5
      idvt(2) = (xmax - wwnvp(1))*scx + wwnvp(5) + .5
      scy = (wwnvp(8) - wwnvp(7))/(wwnvp(4) - wwnvp(3))
      idvt(3) = (ymin - wwnvp(3))*scy + wwnvp(7) + .5
      idvt(4) = (ymax - wwnvp(3))*scy + wwnvp(7) + .5
c return size of a window
      call gitsiz(lx,ly)
c clipping is off
      if (kcf.eq.0) then
         rect(1) = ly - wwnvp(8) - 1
         rect(2) = wwnvp(5)
         rect(3) = ly - wwnvp(7) - 1
         rect(4) = wwnvp(6) 
c clipping is on
      else
         rect(1) = ly - idvt(4) - 1
         rect(2) = idvt(1)
         rect(3) = ly - idvt(3) - 1
         rect(4) = idvt(2)
      endif
c set and save clipping region
      call stclip(rect)
      return
      end
      subroutine gqln(ierr,ltype)
c inquire linetype
c input arguments: none
c ierr = error indicator (0=inquiry successful)
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c ltc = current line type
      ltype = ltc
      ierr = 0
      return
      end
      subroutine gslwsc(alwsc)
c set linewidth scale factor
c input arguments: all
c alwsc = linewidth scale factor, (alwsc > 1.0)
      integer*2 ix, iy
      lws = alwsc
      if (lws.lt.1) lws = 1
      if (lws.gt.3) lws = 3
      ix = lws
      iy = lws
c set dimensions of pen for current GrafPort
      call PenSize(val2(ix),val2(iy))
      return
      end
      subroutine gsplci(icol)
c set polyline color index
c input arguments: all
c icol = color index
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcl = current line color
      lcl = icol
      return
      end
      subroutine gsln(ltype)
c set linetype
c input arguments: all
c ltype = linetype, 1 = solid, 2 = dash, 3 = dot, 4 = dash-dot
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c ltc = current line type
      ltc = 1
      if ((ltype.ge.1).and.(ltype.le.5)) ltc = ltype
      return
      end
      subroutine gpl(n,px,py)
c draw polyline
c input arguments: all
c n = number of points to be connected by a polyline
c px/py = x/y coordinates of points in world coordinates
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c lcl = current line color
c myls = y screen offset, used because origin is at top of screen
c ltc = current line type
c idvt = device transformation
c trans = normalization transformations
      dimension px(n), py(n)
      integer*2 ix, iy, my
c select current color
      call scolor(lcl)
      my = myls
c calculate transformation factors
      m = nrtn + 1
      scx = (idvt(2) - idvt(1))/(trans(6,m) - trans(5,m))
      aminx = idvt(1) - trans(5,m)*scx + .5
      scy = (idvt(4) - idvt(3))/(trans(8,m) - trans(7,m))
      aminy = idvt(3) - trans(7,m)*scy + .5
c convert to screen coordinates
      ix = px(1)*scx + aminx
      iy = py(1)*scy + aminy
c set pen location without drawing
      call MoveTo(val2(ix),val2(my - iy))
      do 10 j = 2, n
c convert to screen coordinates
      ix = px(j)*scx + aminx
      iy = py(j)*scy + aminy
      call dashln(ix,my-iy,ltc)
   10 continue
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c mtc = current marker type
      mtype = mtc
      ierr = 0
      return
      end
      subroutine gsmksc(amksc)
c set marker size scale factor
c input arguments: all
c amksc = linewidth scale factor, (amksc > 1.0)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c ams = current marker size scale factor
      ams = amksc
      return
      end
      subroutine gspmci(imcol)
c set polymarker color index
c input arguments: all
c imcol = polymarker color index
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcm = current marker color
      lcm = imcol
      return
      end
      subroutine gsmk(mtype)
c set marker type
c input arguments: all
c mtype = marker type
c 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross, 6 = square,
c 7 = square with cross, 8 = diamond, 9 = diamond with cross,
c 10 = filled circle
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c mtc = current marker type
      if ((mtype.ge.1).and.(mtype.le.10)) then
c set marker symbol
         mtc = mtype
      else
c set marker symbol to star
         mtc = 3
      endif
      return
      end
      subroutine gpm(n,px,py)
c draw polymarker
c input arguments: all
c n = number of points to be connected by a polymarker
c px/py = x/y coordinates of points in world coordinates
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c lcm = current marker color
c myls = y screen offset, used because origin is at top of screen
c idvt = device transformation
c trans = normalization transformations
      dimension px(n), py(n)
      character*1 amks(10)
      integer*2 CharWidth
      external CharWidth     
      integer*2 ix, iy, my, if(4)
      save amks
      data amks /'.','+','*','o','x','#','H','V','W','@'/
c select current color
      call scolor(lcm)
      my = myls
c calculate transformation factors
      m = nrtn + 1
      scx = (idvt(2) - idvt(1))/(trans(6,m) - trans(5,m))
      aminx = idvt(1) - trans(5,m)*scx + .5
      scy = (idvt(4) - idvt(3))/(trans(8,m) - trans(7,m))
      aminy = idvt(3) - trans(7,m)*scy + .5
c special case of dot markers
      if (mtc.eq.1) then
         do 10 j = 1, n
c convert to screen coordinates
         ix = px(j)*scx + aminx
         iy = py(j)*scy + aminy
c set pen location without drawing
         call MoveTo(val2(ix),val2(my - iy))
c draw a line to specified coordinates
         call LineTo(val2(ix),val2(my - iy))
   10    continue
c other markers
      else
c set character height
         is = 10.*ams
c set point size for subsequent text drawing
         call TextSize(val2(is))
c obtain font sizing information
         call GetFontInfo(if)
         ich = if(1)
c get width of one character
         icw = CharWidth(val2(ichar(amks(mtc))))
c calculate shift
         dx = .5*float(icw)
         dy = .5*float(ich)
c shift origin
         aminx = aminx - dx
         aminy = aminy - dy
         do 20 j = 1, n
c convert to screen coordinates
         ix = px(j)*scx + aminx
         iy = py(j)*scy + aminy
c set pen location without drawing
         call MoveTo(val2(ix),val2(my - iy))
c draw a character at current pen location
         call DrawChar(val2(ichar(amks(mtc))))
   20    continue
      endif
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itx = current horizontal text alignment
c ity = current vertical text alignment
      ierr = 0
      itxalh = itx
      itxalv = ity
      return
      end
      subroutine gqchh(ierr,chh)
c inquire character height
c input arguments: none
c ierr = error indicator
c chh = character height, in world coordinates
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c idvt = device transformation
c trans = normalization transformations
      integer*2 if(4)
c obtain font sizing information
      call GetFontInfo(if)
      ich = if(1)
c convert to world coordinates
      m = nrtn + 1
      scy = (trans(8,m) - trans(7,m))/float(idvt(4) - idvt(3))
      chh = float(ich)*scy
      ierr = 0
      return
      end
      subroutine gqtxp(ierr,itxp)
c inquire text path
c input arguments: none
c ierr = error indicator
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itxtp = current text path
      itxp = itxtp
      ierr = 0
      return
      end
      subroutine gqchup(ierr,chux,chuy)
c inquire character up vector
c input arguments: none
c ierr = error indicator
c chux/chuy = up vector, in world coordinates
      chux = 0.
      chuy = 1.
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itx = current horizontal text alignment
c ity = current vertical text alignment
      if ((itxalh.ge.0).and.(itxalh.le.3)) itx = itxalh
      if ((itxalv.ge.0).and.(itxalh.le.3)) ity = itxalv
      return
      end
      subroutine gstxfp(nfont,iprec)
c set text font
c input arguments: all
c nfont = character font number
c iprec = text precision (0=string,1=character,2=stroke)
c Monaco
      if (nfont.eq.1) then
         nft = 4
c systemFont, Chicago
      elseif (nfont.eq.2) then
         nft = 0
c Symbol
      elseif (nfont.eq.3) then
         nft = 23
c applFont, Geneva
      else
         nft = 1
      endif
c select font for subsequent text drawing
      call TextFont(val2(nft))
      return
      end
      subroutine gstxp(itxp)
c set text path
c input arguments: all
c itxp = text path (0=right(default), 1=left, 2=up, 3=down)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c itxtp = current text path
      itxtp = itxp
      return
      end
      subroutine gstxci(itcol)
c set text color index
c input arguments: all
c itcol = text color index
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lct = current text color
      lct = itcol
      return
      end
      subroutine gschh(chh)
c set character height
c input arguments: all
c chh = character height, in world coordinates
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c cch = current character height
      cch = chh
      return
      end
      subroutine gschup(chux,chuy)
c set character up vector
c input arguments: all
c chux/chuy = up vector, in world coordinates
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c chuv = character up vector
      chuv(1) = chux
      chuv(2) = chuy
      return
      end
      subroutine gschxp(chxp)
c set character expansion factor
c input arguments: all
c chxp = character expansion factor (>0.)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c chx = character width expansion factor
      if (chxp.gt.0.) then
         chx = chxp
      else
         chx = 1.0
      endif
      return
      end
      subroutine gtx(px,py,chars)
c display text
c input arguments: all
c px/py = starting x/y position of text, in world coordinates
c chars = test string to be displayed
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c itx = current horizontal text alignment
c ity = current vertical text alignment
c lct = current text color
c myls = y screen offset, used because origin is at top of screen
c idvt = device transformation
c cch = current character height
c trans = normalization transformations
      character*(*) chars
      integer*2 TextWidth
      external TextWidth
      integer*2 ix, iy, my, if(4)
c select current
      call scolor(lct)
      n = len(chars)
      my = myls
c calculate transformation factors
      m = nrtn + 1
      scx = (idvt(2) - idvt(1))/(trans(6,m) - trans(5,m))
      aminx = idvt(1) - trans(5,m)*scx + .5
      scy = (idvt(4) - idvt(3))/(trans(8,m) - trans(7,m))
      aminy = idvt(3) - trans(7,m)*scy + .5
c convert to screen coordinates
      ix = px*scx + aminx
      iy = py*scy + aminy
      is = cch*scy
c set point size for subsequent text drawing
      call TextSize(val2(is))
c obtain font sizing information
      call GetFontInfo(if)
      ich = if(1)
c get width of unformatted text
      icw = TextWidth(chars,val2(0),val2(n))
c determine horizontal offset
      if (itx.eq.2) then
         ix = ix - .5*float(icw)
      elseif (itx.eq.3) then
         ix = ix - icw
      endif
c determine vertical offset
      if (ity.eq.3) then
         iy = iy - .5*float(ich)
      elseif ((ity.eq.1).or.(ity.eq.2)) then
         iy = iy - ich
      endif
c set pen location without drawing
      call MoveTo(val2(ix),val2(my - iy))
c draw text from any arbitrary buffer
      call DrawText(chars,val2(0),val2(n))
      return
      end
      subroutine gsfaci(ifcol)
c set fill area color index
c input arguments: all
c ifcol = color index
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lcf = current fill color
      lcf = ifcol
      return
      end
      subroutine gsfais(ints)
c set fill area interior style
c input arguments: all
c ints = desired interior style:
c 0 = hollow (default), 1 = solid, 2 = pattern, 3 = hatch
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c infs = interior fill style
      infs = ints
      return
      end
      subroutine gfa(n,px,py)
c fill area
c input arguments: all
c n = number of points in fill area
c px,py = arrays of points, in world coordinates
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c lcf = current fill color
c myls = y screen offset, used because origin is at top of screen
c infs = interior fill style
c idvt = device transformation
c trans = normalization transformations
      dimension px(n), py(n)
      integer OpenPoly
      external OpenPoly
      integer*2 ix, iy, my, jx, jy, pattern(4)
c select current color
      call scolor(lcf)
      my = myls
c calculate transformation factors
      m = nrtn + 1
      scx = (idvt(2) - idvt(1))/(trans(6,m) - trans(5,m))
      aminx = idvt(1) - trans(5,m)*scx + .5
      scy = (idvt(4) - idvt(3))/(trans(8,m) - trans(7,m))
      aminy = idvt(3) - trans(7,m)*scy + .5
c start recording of a polygon definition
      mpoly = OpenPoly()
c balance a previous HidePen; make pen visible
      if (infs.eq.0) call ShowPen
c convert to screen coordinates
      jx = px(1)*scx + aminx
      jy = py(1)*scy + aminy
c set pen location without drawing
      call MoveTo(val2(jx),val2(my - jy))
      do 10 j = 2, n
c convert to screen coordinates
      ix = px(j)*scx + aminx
      iy = py(j)*scy + aminy
c draw a line to specified coordinates
      call LineTo(val2(ix),val2(my-iy))
   10 continue
c draw a line to specified coordinates
      call LineTo(val2(jx),val2(my-jy))
c make subsequent pen motion invisible
      if (infs.eq.0) call HidePen
c stop recording polygon vertices
      call ClosePoly
c solid fill
      if (infs.eq.1) then
c get pattern from an indexed 'PAT#' resource
         call GetIndPattern(pattern,val2(0),val2(1))
c Fill interior of a polygon with specified pattern
         call FillPoly(val(mpoly),pattern)
c pattern fill
      elseif (infs.eq.2) then
c get pattern from an indexed 'PAT#' resource
         call GetIndPattern(pattern,val2(0),val2(12))
c Fill interior of a polygon with specified pattern
         call FillPoly(val(mpoly),pattern)
c hatch fill
      elseif (infs.eq.3) then
c get pattern from an indexed 'PAT#' resource
         call GetIndPattern(pattern,val2(0),val2(16))
c Fill interior of a polygon with specified pattern
         call FillPoly(val(mpoly),pattern)
      endif
c deallocate all storage for a polygon
      call KillPoly(val(mpoly))
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
      parameter(lxm=1024,lym=1024)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c myls = y screen offset, used because origin is at top of screen
c idvt = device transformation
c wwnvp = workstation window and viewport
c trans = normalization transformations
      dimension icola(icxd,icyd)
      character*1 image(lxm*lym)
      integer*2 rect(4)
c lxn, lyn = maximum address of pixels in x, y
      lxn = wwnvp(6)
      lyn = wwnvp(8)
c find location of upper left and lower right hand corner of image
      xu = px
      xl = qx
      yu = amax1(py,qy)
      yl = amin1(py,qy)
c calculate transformation factors
      m = nrtn + 1
      scx = float(idvt(2) - idvt(1))/(trans(6,m) - trans(5,m))
      aminx = float(idvt(1)) - trans(5,m)*scx
      scy = float(idvt(4) - idvt(3))/(trans(8,m) - trans(7,m))
      aminy = float(idvt(3)) - trans(7,m)*scy
c convert to screen coordinates
      ax = xu*scx + aminx
      ay = yu*scy + aminy
      bx = xl*scx + aminx
      by = yl*scy + aminy
c clip to workstation viewport
      ix0 = amin1(amax1(ax,wwnvp(5)),wwnvp(6)) + .5
      iy0 = amin1(amax1(ay,wwnvp(7)),wwnvp(8)) + .5
      ix1 = amin1(amax1(bx,wwnvp(5)),wwnvp(6)) + .5
      iy1 = amin1(amax1(by,wwnvp(7)),wwnvp(8)) + .5
c calculate size of rescaled raster image
      rect(1) = myls - iy0
      rect(2) = ix0
      rect(3) = myls - iy1 + 1
      rect(4) = ix1 + 1
c comment: image dimension must be even, preferrably a multiple of 4
      ldx = 2*(idx/2)
      ldy = min0(idy,lxm*lym/ldx)
c outer loop over rows
      do 20 k = 1, ldy
      kk = nrs + k - 1
      joff = ldx*(k - 1)
c inner loop over bytes in row
      do 10 j = 1, ldx
      image(j+joff) = char(icola(j+ncs-1,kk))
   10 continue
   20 continue
c copy offscreen image to screen
      call copyimg(image,rect,8,ldx,ldy)
      return
      end
      subroutine gqclip(ierr,indcl,clrect)
c inquire clipping indicator
c input arguments: none
c ierr = error indicator
c indcl = clipping indicator (0=no clip, 1=clip)
c clrect = clipping rectangle, in ndc
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c nrtn = current normalization transformation number
c kcf = clipping flag
c trans = normalization transformations
      dimension clrect(4)
      indcl = kcf
c find clipping rectangle
      n = nrtn + 1
      do 10 j = 1, 4
      clrect(j) = trans(j,n)
   10 continue
      ierr = 0
      return
      end
      subroutine gsclip(iclsw)
c set clipping indicator
c input arguments: all
c iclsw = clipping switch (0=no clip,1=clip)
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kcf = clipping flag
c idvt = device transformation
c wwnvp = workstation window and viewport
      integer*2 rect(4)
c return size of a window
      call gitsiz(lx,ly)
c clipping is on, but was off
      if ((kcf.eq.0).and.(iclsw.eq.1)) then
         rect(1) = ly - idvt(4) - 1
         rect(2) = idvt(1)
         rect(3) = ly - idvt(3) - 1
         rect(4) = idvt(2)
c set and save clipping region
         call stclip(rect)
c clipping is off, but was on
      elseif ((kcf.eq.1).and.(iclsw.eq.0)) then
         rect(1) = ly - wwnvp(8) - 1
         rect(2) = wwnvp(5)
         rect(3) = ly - wwnvp(7) - 1
         rect(4) = wwnvp(6)
c set and save clipping region
         call stclip(rect)
      endif
      kcf = iclsw
      return
      end
      subroutine guwk(idwk,iregfl)
c update workstation
c input arguments: all
c idwk = workstation identifier
c iregfl = regeneration flag (0=postponed,1=perform)
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c wptr = pointer to cwindow structure
      integer GetWindowPort
      external GetWindowPort
c flush window buffer to screen
      call QDFlushPortBuffer(val4(GetWindowPort(val4(wptr))),val4(0))
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lct = current text color
c myls = y screen offset, used because origin is at top of screen
c cea = character echo area
      character*(*) str
      integer*2 myEventMask, mouseDown, keyDown, autoKey, updateEvt
      integer*2 osEvt, resus
c myEventMask looks for mouse, keyboard, update, system and quit events
c     parameter(myEventMask=1150-32768,mouseDown=1,keyDown=3,autoKey=5)
c myEventMask looks for mouse, keyboard, system and quit events
      parameter(myEventMask=1086-32768,mouseDown=1,keyDown=3,autoKey=5)
      parameter(updateEvt=6,osEvt=15,resus=256)
      integer*2 ix, iy, my, color(3), event(8)
      lostr = 0
      n = len(str)
      my = myls
      ix = cea(1)
      iy = cea(3)
c return current foreground color
      call GetForeColor(color)
c select current color to current text color
      call scolor(lct)
c multifinder-aware way to obtain events
   10 call WaitNextEvent(val2(myEventMask),event,val4(-1),val4(0))
      if ((event(1).eq.osEvt).and.(iand(event(2),resus).eq.resus)) then
c resume event
         if ((event(3)-2*(event(3)/2)).eq.1) then
            call savscrn(-1,ipch)
c suspend event
         else
            call savscrn(1,ipch)
         endif
c try again
         go to 10
      endif
      if (event(1).eq.mouseDown) then
c check alternative mouse events
         call auxmaus(event,ierr,str)
c get another event for window or processed file menu events
c this alternative will always wait for a keyboard event
c        if (ierr.le.0) then
c get another event for processed file menu or drag events only
c this alternative treats window events as a carriage return
         if (ierr.lt.0) then
            go to 10
c process keyboard alternative
         elseif (ierr.gt.0) then
c numbers
            if (ierr.gt.128) then
               lostr = ierr - 128
c single character code
            else
               lostr = 1
            endif
         endif
         go to 20
      endif
c try again if not keyboard
      if ((event(1).ne.keyDown).and.(event(1).ne.autoKey)) go to 10
      is = event(3) - 256*(event(3)/256)
c quit if control character or too many characters
      if ((is.lt.32).or.(lostr.ge.n)) go to 20
      lostr = lostr + 1
      str(lostr:lostr) = char(is)
c echo character
c set pen location without drawing
      call MoveTo(val2(ix),val2(my - iy))
c draw text from any arbitrary buffer
      call DrawText(str(1:lostr),val2(0),val2(lostr))
c look for more characters
      go to 10
c restore previous color
c set foreground color to best match for current device
   20 call RGBForeColor(color)
      istat = 1
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c myls = y screen offset, used because origin is at top of screen
c wwnvp = workstation window and viewport
c trans = normalization transformations
      integer*2 myEventMask, mouseDown, updateEvt, osEvt, resus
c myEventMask looks for mouse, update, system and quit events
c     parameter(myEventMask=1094-32768,mouseDown=1,updateEvt=6)
c myEventMask looks for mouse, system and quit events
      parameter(myEventMask=1030-32768,mouseDown=1,updateEvt=6)
      parameter(osEvt=15,resus=256)
      integer*2 event(8)
      character*1 str
c multifinder-aware way to obtain events
   10 call WaitNextEvent(val2(myEventMask),event,val4(-1),val4(0))
c check for system events
      if ((event(1).eq.osEvt).and.(iand(event(2),resus).eq.resus)) then
c resume event
         if ((event(3)-2*(event(3)/2)).eq.1) then
            call savscrn(-1,ipch)
c suspend event
         else
            call savscrn(1,ipch)
         endif
c try again
         go to 10
      endif
c try again if not mouse
      if (event(1).ne.mouseDown) go to 10
c check alternative mouse events
      call auxmaus(event,ierr,str)
c get another event for all non-window events
      if (ierr.ne.0) go to 10
c obtain local coordinates of global point
      call GlobalToLocal(event(6))
      qx = event(7)
      qy = myls - event(6)
c calculate transformation factors
      scx = (trans(6,1) - trans(5,1))/(wwnvp(6) - wwnvp(5))
      scy = (trans(8,1) - trans(7,1))/(wwnvp(8) - wwnvp(7))
c convert from device to world coordinates
      px = (qx - wwnvp(5))*scx + trans(5,1)
      py = (qy - wwnvp(7))*scy + trans(7,1)
      nrt = 0
      istat = 1
      return
      end
      subroutine gsstm(idwk,idstr,mode,iesw)
c set string mode
c input arguments: all
c idwk = workstation identifier
c idstr = string device number
c mode  =  mode of operation (0 = request,1 = sample,2 = event)
c iesw  =  echo switch (0 = no echo,1 = echo)
      integer GetMenuHandle
      external GetMenuHandle
      character*256 str
c given a menu ID, obtain a handle to the menu
      menuh = GetMenuHandle(val2(3))
c get the text of a menu item
      call GetMenuItemText(val4(menuh),val2(3),str)
      if ((str(2:8).eq.'Animate').and.(mode.eq.2)) then
c change the text of a menu item
         call SetMenuItemText(val4(menuh),val2(3),char(5)//'Pause')
      elseif ((str(2:6).eq.'Pause').and.(mode.eq.0)) then
c change the text of a menu item
         call SetMenuItemText(val4(menuh),val2(3),char(7)//'Animate')
      endif
      return
      end
      subroutine gslcm(idwk,idloc,mode,iesw)
c set locator mode
c input arguments: all
c idwk = workstation identifier
c idloc = locator device number
c mode  =  mode of operation (0 = request,1 = sample,2 = event)
c iesw  =  echo switch (0 = no echo,1 = echo)
      integer GetMenuHandle
      external GetMenuHandle
      character*256 str
c given a menu ID, obtain a handle to the menu
      menuh = GetMenuHandle(val2(3))
c get the text of a menu item
      call GetMenuItemText(val4(menuh),val2(3),str)
      if ((str(2:8).eq.'Animate').and.(mode.eq.2)) then
c change the text of a menu item
         call SetMenuItemText(val4(menuh),val2(3),char(5)//'Pause')
      elseif ((str(2:6).eq.'Pause').and.(mode.eq.0)) then
c change the text of a menu item
         call SetMenuItemText(val4(menuh),val2(3),char(7)//'Animate')
      endif
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
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c iwk = current workstation identifier
c idst = current string device number
c idlc = current locator device number
      common /macevnts/ sevent, mevent, strevt
      integer*2 sevent(8), mevent(8)
      character*8 strevt
c sevent = last string event
c mevent = last mouse event
c strevt = last string emulation event
      logical*1 EventAvail, GetNextEvent
      integer*2 FindWindow
      integer TickCount
      external TickCount, EventAvail, GetNextEvent, FindWindow
      integer*2 myEventMask, mouseDown, keyDown, autoKey, updateEvt
      integer*2 osEvt, resus
c myEventMask looks for mouse, keyboard, update, system and quit events
c     parameter(myEventMask=1150-32768,mouseDown=1,keyDown=3,autoKey=5)
c myEventMask looks for mouse, keyboard, system and quit events
      parameter(myEventMask=1086-32768,mouseDown=1,keyDown=3,autoKey=5)
      parameter(updateEvt=6,osEvt=15,resus=256)
      integer*2 event(8), intw(2), mouseloc
      integer longw, wptr, loc
      equivalence(longw,intw)
      idwk = iwk
c give DAs a chance to perform periodic actions
c     call SystemTask()
c get current system tick count
      jclock = TickCount() + 60.*tout
c get an event without removing it from the queue
   10 if (EventAvail(val2(myEventMask),event)) then
c obtain next available event of specified type(s)
         call GetNextEvent(val2(myEventMask),event)
c keyboard event
         if ((event(1).eq.keyDown).or.(event(1).eq.autoKey)) then
            do 20 i = 1, 8
            sevent(i) = event(i)
   20       continue
            icl = 6
            idnr = idst
c mouse event
         elseif (event(1).eq.mouseDown) then
            intw(1) = event(6)
            intw(2) = event(7)
c loc = position in global coordinates where mouse event took place
            loc = longw
c see which window part, including menu bar, is at a point
            mouseloc = FindWindow(val4(loc),wptr)
c mouse event in content region
            if (mouseloc.eq.3) then
               do 30 i = 1, 8
               mevent(i) = event(i)
   30          continue
               icl = 1
               idnr = idlc
c possible keyboard simulation event
            else
c check alternative mouse events
               call auxmaus(event,ierr,strevt)
               if (ierr.ge.0) then
                  sevent(1) = mouseDown
                  sevent(2) = ierr
                  icl = 6
                  idnr = idst
               else
                  icl = -1
                  idnr = 0
               endif
            endif
c reject any other event
         else
c check for system events
            if ((event(1).eq.osEvt).and.(iand(event(2),resus).eq.resus))
     1then
c resume event
               if ((event(3)-2*(event(3)/2)).eq.1) then
                  call savscrn(-1,ipch)
c suspend event
               else
                  call savscrn(1,ipch)
               endif
            endif
            icl = -1
            idnr = 0
         endif
c no event
      else
         icl = 0
         idnr = 0
      endif
c check for other events
      if (icl.lt.0) go to 10
c wait for time out period if no event
      if ((icl.eq.0).and.(TickCount().lt.jclock)) go to 10
      return
      end
      subroutine ggtst(lostr,str)
c get string
c input arguments: none
c lostr = number of characters in string
c str = input string
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c lct = current text color
c myls = y screen offset, used because origin is at top of screen
c cea = character echo area
      common /macevnts/ sevent, mevent, strevt
      integer*2 sevent(8), mevent(8)
      character*8 strevt
c sevent = last string event
c mevent = last mouse event
c strevt = last string emulation event
      character*(*) str
      integer mouseDown, keyDown, autoKey
      parameter(mouseDown=1,keyDown=3,autoKey=5)
      integer*2 ix, iy, my, color(3)
      lostr = 1
      my = myls
      ix = cea(1)
      iy = cea(3)
c return current foreground color
      call GetForeColor(color)
c select current color to current text color
      call scolor(lct)
      if (sevent(1).ne.0) then
c check if keyboard event
         if ((sevent(1).eq.keyDown).or.(sevent(1).eq.autoKey)) then
            is = sevent(3) - 256*(sevent(3)/256)
            if (is.lt.32) then
               lostr = 0
            else
               str(1:1) = char(is)
            endif
c check alternative mouse events
         elseif (sevent(1).eq.mouseDown) then
            ierr = sevent(2)
            str = strevt
            strevt = '        '
c this alternative treats window events as a carriage return
            if (ierr.le.0) then
               lostr = 0
c process numbers keyboard alternative
            elseif (ierr.gt.128) then
               lostr = ierr - 128
            endif
         endif
c echo character
c set pen location without drawing
         call MoveTo(val2(ix),val2(my - iy))
c draw text from any arbitrary buffer
         call DrawText(str(1:lostr),val2(0),val2(lostr))
c clear last string event
         do 10 i = 1, 8
         sevent(i) = 0
   10    continue
      endif
      return
      end
      subroutine ggtlc(nrt,px,py)
c get locator
c input arguments: none
c nrt = normalization transformation number
c px/py = position, in world coordinates
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c myls = y screen offset, used because origin is at top of screen
c wwnvp = workstation window and viewport
c trans = normalization transformations
      common /macevnts/ sevent, mevent, strevt
      integer*2 sevent(8), mevent(8)
      character*8 strevt
c sevent = last string event
c mevent = last mouse event
      integer*2 mouseDown
      parameter(mouseDown=1)
      if (mevent(1).ne.mouseDown) return
c obtain local coordinates of global point
      call GlobalToLocal(mevent(6))
      qx = mevent(7)
      qy = myls - mevent(6)
c calculate transformation factors
      scx = (trans(6,1) - trans(5,1))/(wwnvp(6) - wwnvp(5))
      scy = (trans(8,1) - trans(7,1))/(wwnvp(8) - wwnvp(7))
c convert from device to world coordinates
      px = (qx - wwnvp(5))*scx + trans(5,1)
      py = (qy - wwnvp(7))*scy + trans(7,1)
      nrt = 0
c clear last mouse event
      do 10 i = 1, 8
      mevent(i) = 0
   10 continue
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
      subroutine mcdflts
c this subroutines creates default tables for macintosh driver for gks
      parameter(maxt=3)
      common /gksmc/ kosv,iwk,idc,nrtn,itx,ity,kcf,lcl,lcm,lct,lcf,myls,
     1ltc,mtc,infs,itxtp,idst,idlc,idvt,ams,chx,cch,wwnvp,cea,trans,chuv
      dimension idvt(4), wwnvp(8), cea(4), trans(8,maxt), chuv(2)
c kosv = operating state value
c iwk = current workstation identifier
c idc = current connection identifier
c nrtn = current normalization transformation number
c itx = current horizontal text alignment
c ity = current vertical text alignment
c kcf = clipping flag
c lcl = current line color
c lcm = current marker color
c lct = current text color
c lcf = current fill color
c myls = y screen offset, used because origin is at top of screen
c ltc = current line type
c mtc = current marker type
c infs = interior fill style
c itxtp = current text path
c idst = current string device number
c idlc = current locator device number
c idvt = device transformation
c ams = current marker size scale factor
c chx = character width expansion factor
c cch = current character height
c wwnvp = workstation window and viewport
c cea = character echo area
c trans = normalization transformations
c chuv = character up vector
      integer*2 rect(4), color(3), ix, iy
c set default workstation window to square
      wwnvp(1) = 0.
      wwnvp(2) = 1.
      wwnvp(3) = 0.
      wwnvp(4) = 1.
c set default workstation viewport to square
c return size of a window
      call gitsiz(lx,ly)
c lx/ly = current window width/height
      if (lx.gt.ly) then
         wwnvp(5) = .5*float(lx - ly)
         wwnvp(8) = float(ly - 1)
         wwnvp(6) = wwnvp(5) + wwnvp(8)
         wwnvp(7) = 0.
      else
         wwnvp(6) = float(lx - 1)
         wwnvp(7) = .5*float(ly - lx)
         wwnvp(5) = 0.
         wwnvp(8) = wwnvp(6) + wwnvp(7)
      endif
c make default character echo area equal to workstation viewport
      do 10 j = 1, 4
      cea(j) = wwnvp(j+4)
   10 continue
c default transformation
      nrtn = 0
      trans(1,1) = 0.
      trans(2,1) = 1.
      trans(3,1) = 0.
      trans(4,1) = 1.
      trans(5,1) = 0.
      trans(6,1) = 1.
      trans(7,1) = 0.
      trans(8,1) = 1.
      do 30 k = 2, maxt
      do 20 j = 1, 8
      trans(j,k) = trans(j,1)
   20 continue
   30 continue
c convert from viewport to screen coordinates
      idvt(1) = wwnvp(5) + .5
      idvt(2) = wwnvp(6) + .5
      idvt(3) = wwnvp(7) + .5
      idvt(4) = wwnvp(8) + .5
c set default clipping state to on
      kcf = 1
c default colors
c set background color to black
      icr = 0
      icg = 0
      icb = 0
c set a color map entry to a specified RGB value
      call crctab(0,icr,icg,icb)
      color(1) = icr
      color(2) = icg
      color(3) = icb
c set background color to best match for current device
      call RGBBackColor(color)
c set foreground color to white
      icr = -1
      icg = -1
      icb = -1
c set a color map entry to a specified RGB value
      call crctab(1,icr,icg,icb)
c select current color
      call scolor(1)
      lcl = 1
      lcm = lcl
      lct = lcl
      lcf = lcl
      rect(1) = 0
      rect(2) = 0
      rect(3) = ly
      rect(4) = lx
c fill rectangle with background pattern
      call EraseRect(rect)
      rect(1) = ly - idvt(4) - 1
      rect(2) = idvt(1)
      rect(3) = ly - idvt(3) - 1
      rect(4) = idvt(2)
c set and save clipping region
      call stclip(rect)
c default line style
c set line style to solid
      ltc = 1
c default marker
c set marker symbol to star
      mtc = 3
c default text alignment
      itx = 0
      ity = 0
c default line width
      ix = 1
      iy = 1
c set dimensions of pen for current GrafPort
      call PenSize(val2(ix),val2(iy))
c default marker size scale factor
      ams = 1.0
c itxtp = current text path
      itxtp = 0
c default character width expansion factor
      chx = 1.0
c default character height 
      cch = 0.01
c set default character up vector
      chuv(1) = 0.
      chuv(2) = 1.
c set default fill area interior style to hollow
      infs = 0
c set default string and locator device numbers
      idst = 0
      idlc = 0
      return
      end
c internal library for macintosh images
      subroutine setpsz(mx,my,lx,ly)
c this subroutine defines a safe window size
c window placement is requested and a safe window placement is stored
c in wrect
c mx, my = location of upper left corner
c lx, ly = preferred size of window
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c wrect = window rectangle
      integer GetMainDevice, GetDeviceList, GetNextDevice
      external GetMainDevice, GetDeviceList, GetNextDevice
      integer handle, cptr
      integer*2 trect(4)
      save /gmcdev/
c get handle to main graphics device that carries a menu bar
      handle = GetMainDevice()
      cptr = long(handle)
c get size of screen
      trect(3) = word(cptr+38)
      trect(4) = word(cptr+40)
c set top, left, bottom, and right sides of main window rect
      trect(1) = min0(trect(3),my)
      trect(2) = min0(trect(4),mx)
      trect(3) = min0(trect(3),ly+trect(1))
      trect(4) = min0(trect(4),lx+trect(2))
c provide a handle to the first gDevice in the DeviceList
      handle = GetDeviceList()
c clear screen size
      wrect(1) = 0
      wrect(2) = 0
      wrect(3) = 0
      wrect(4) = 0
      nscrns = 0
   10 nscrns = nscrns + 1
      cptr = long(handle)
c find smallest rectangle enclosing two rectangles
      call UnionRect(val4(cptr+34),wrect,wrect)
c allocate new gDevice structure
      handle = GetNextDevice(cptr)
      if (handle.ne.0) go to 10
c set top, left, bottom, and right sides of multiple window rects
      if (nscrns.eq.1) then
         wrect(1) = trect(1)
         wrect(2) = trect(2)
         wrect(3) = trect(3)
         wrect(4) = trect(4)
      else
         wrect(1) = wrect(1) + 40
         wrect(2) = wrect(2) + 40
         wrect(3) = wrect(3) - 40
         wrect(4) = wrect(4) - 40
      endif
      return
      end
      subroutine cr8win(name)
c this subroutine creates a window from data stored in wrect
c in addition, a default clipping area is defined and menus are added
c input arguments: name
c name = label for window
      parameter(nclmax=256,lctab=4*(nclmax+1))
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c wrect = window rectangle
c wptr = pointer to cwindow structure
c crect = current clipping region
      character*(*) name
      integer NewCWindow, newMenu, GetWindowPort
      external NewCWindow, newMenu, GetWindowPort
      integer*2 trect(4)
      character*1 str
      l = len(name)
c create a new color window
      wptr = NewCWindow(val4(0),wrect,char(l)//name,val1(.true.),val2(4)
     1,val4(-1),val1(.false.),val4(0))
c activate a GrafPort
      call SetPort(val4(GetWindowPort(val4(wptr))))
c keep a rectangular area from being updated
      call GetPortBounds(GetWindowPort(val4(wptr)),trect)
      call ValidWindowRect(val4(wptr),trect)
c define default current clipping region
      crect(1) = 0
      crect(2) = 0
      crect(3) = wrect(3) - wrect(1) + 1
      crect(4) = wrect(4) - wrect(2) + 1
c menu 1 = apple menu
c get a handle to an empty menu
      menuh = NewMenu(val2(1),char(1)//char(20))
c append names of selected resource type to menu
c    call AddResMenu(val4(menuh),val4('DRVR'))
c add a menu to the menu list
      call InsertMenu(val4(menuh),val2(0))
c menu 2 = file menu
c get a handle to an empty menu
      menuh = NewMenu(val2(2),char(4)//'File')
c add one or more items to a menu
      call AppendMenu(val4(menuh),char(14)//'Save Window...')
c add one or more items to a menu
      call AppendMenu(val4(menuh),char(12)//'Print Window')
c add one or more items to a menu
      call AppendMenu(val4(menuh),char(20)//'Print Inverse Window')
c add one or more items to a menu
      call AppendMenu(val4(menuh),char(4)//'Quit')
c add a menu to the menu list
      call InsertMenu(val4(menuh),val2(0))
c display the titles of all the menus in the menu list
      call DrawMenuBar
c create extra menus
      call mcextra(0,ierr,str)
c allocate relocatable block from current heap zone for color table
      cthandle = NewHandle(val4(2*lctab))
c lock a handle's data area (keep it from moving)
      call HLock(val4(cthandle))
      return
      end
      subroutine gitsiz(lx,ly)
c this subroutine returns information about screen size
c input arguments: none
c lx, ly = the size of the screen
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c wrect = window rectangle
      lx = wrect(4) - wrect(2) + 1
      ly = wrect(3) - wrect(1) + 1
      return
      end
      subroutine gitnc(nc)
c this subroutine returns information about available colors
c input arguments: none
c nc = number of colors
c ncolors = number of available colors
      integer GetMainDevice
      external GetMainDevice
      integer handle, cptr
c get handle to main graphics device that carries a menu bar
      handle = GetMainDevice()
      cptr = long(handle)
c find number of available colors
      handle = long(cptr+22)
      cptr = long(handle)
      ncolors = min(word(cptr+32),16)
      nc = 2**ncolors
      return
      end
      subroutine crctab(ic,icr,icg,icb)
c set a color map entry to a specified RGB value
c input arguments: all
c ic = color index
c icr/icg/icb = rgb color values, -32768 < color value < 32767
      parameter(nclmax=256,lctab=4*(nclmax+1))
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c cthandle = handle to color table structure
      integer GetCTSeed
      integer ctseed
      integer*2 cts2(2)
      equivalence(ctseed,cts2)
      save istart
      data istart /0/
c initialize default table on first entry
      if (istart.eq.0) then
         longw = long(cthandle)
c get unique seed value for color table
         ctseed = GetCTSeed()
         word(longw) = cts2(1)
         word(longw+2) = cts2(2)
c set ctFlags for PixMap
         word(longw+4) = 0
         word(longw+6) = nclmax - 1
c set all colors to white
         do 10 i = 1, nclmax
            icc = longw + 8*i
            word(icc) = i - 1
            word(icc+2) = -1
            word(icc+4) = -1
            word(icc+6) = -1
   10    continue
         istart = 1
      endif
c ignore indices out of range
      if ((ic.lt.0).or.(ic.gt.255)) return
      icc = long(cthandle) + 8*ic + 8
      word(icc) = ic
      word(icc+2) = icr
      word(icc+4) = icg
      word(icc+6) = icb
      return
      end
      subroutine scolor(ic)
c select current color
c input arguments: all
c ic = color index
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c cthandle = handle to color table structur
      data lastic /-1/
      save lastic
      if ((ic.lt.0).or.(ic.gt.255)) return
      if (ic.ne.lastic) then
c set foreground color to best match for current device
         call RGBForeColor(val4(long(cthandle)+8*ic+10))
c        lastic = ic
      endif
      return
      end
      subroutine stclip(rect)
c this subroutine sets and saves clipping rectangle
c input arguments: all
c rect = clipping rectangle
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c crect = current clipping region
      integer*2 rect(4)
      crect(1) = rect(1)
      crect(2) = rect(2)
      crect(3) = rect(3) + 1
      crect(4) = rect(4) + 1
c set clipping region to a rectangle
      call ClipRect(crect)
      return
      end
      subroutine eraser
c erase screen
c this subroutine sets and stores clipping rectangle
c input arguments: none
c rect = clipping rectangle
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
      integer*2 rect(4)
c wrect = window rectangle
c crect = current clipping region
      rect(1) = 0
      rect(2) = 0
      rect(3) = wrect(3) - wrect(1) + 1
      rect(4) = wrect(4) - wrect(2) + 1
c set clipping region to a rectangle
      call ClipRect(rect)
c fill rectangle with background pattern
      call EraseRect(rect)
c set clipping region to a rectangle
      call ClipRect(crect)
      return
      end
      subroutine copyimg(image,rect,nbit,idx,idy)
c this subroutine copies offscreen image to screen
c input arguments: all
c image = input image data
c rect = destination image size
c nbit = bits per pixel (always a power of 2)
c idx/idy = image dimensions
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c wptr = pointer to cwindow structure
c cthandle = handle to color table structure
      integer NewGWorld, GetGWorldPixMap, GetPortBitMapForCopyBits
      integer GetWindowPort
      external NewGWorld, GetGWorldPixMap, GetPortBitMapForCopyBits
      external GetWindowPort
      character*1 image(idx*idy)
      integer*2 rect(4)
      integer*2 nbit2, ierr, rbytes
      integer*2 prect(4), color(3)
      integer gwptr, pxmptr, pxmhdl, baddr, input, output
c set source size
      prect(1) = 0
      prect(2) = 0
      prect(3) = idy
      prect(4) = idx
      nbit2 = nbit
      ierr = NewGWorld(gwptr,val2(nbit2),prect,val4(cthandle),val4(0),va
     1l4(0))
c get a handle to the pixel map for an offscreen graphics world
      pxmhdl = GetGWorldPixMap(val4(gwptr))
c lock the offscreen buffer in memory for duration of a draw
      call LockPixels(val4(pxmhdl))
      pxmptr = long(pxmhdl)
      baddr = long(pxmptr)
      rbytes = iand(word(pxmptr+4),16383)
c copybits needs background white and foreground black
      call srcolor(1,color)
c copy image to pixel map
      do 112 k = 1, idy
      joff = baddr + rbytes*(k - 1) - 1
      ioff = idx*(k - 1)
      do 111 j = 1, idx
      byte(j+joff) = ichar(image(j+ioff))
  111 continue
  112 continue
      input = GetPortBitMapForCopyBits(val4(gwptr))
      output = GetPortBitMapForCopyBits(val4(GetWindowPort(val4(wptr))))
c copy bitMap, with optional scaling, clipping, etc.
      call CopyBits(val4(input),val4(output),prect,rect,val2(64),val(0))
c restore background to black and foreground to previous color
      call srcolor(-1,color)
c unlock the buffer used by an offscreen graphics world
      call UnlockPixels(val4(pxmhdl))
c dispose of a GWorld structure and and substructures
      call DisposeGWorld(val4(gwptr))
      return
      end
      subroutine srcolor(isr,color)
c this subroutine saves or restores current foreground color to make
c copybits work properly.
c if isr = 1, current color is stored in array color, background is
c set to white, and foreground to black
c if isr = 0, current color is stored in array color, background is
c set to black, and foreground to white
c if isr = -1, background color is set to black, and current color is
c set to array color
c input arguments: isr, color (only if isr = -1)
c output arguments: color (only if isr = 0,1)
c isr = (-1,0,1) = (restore,invert-save,save) color
c color = rgb values of color
      integer*2 color(3)
      integer*2 black(3), white(3)
      data black /0,0,0/
      data white /-1,-1,-1/
      save black, white
c normal copy mode
      if (isr.eq.1) then
c return current foreground color
         call GetForeColor(color)
c set background color to best match for current device
         call RGBBackColor(white)
c set foreground color to best match for current device
         call RGBForeColor(black)
c invert copy mode
      elseif (isr.eq.0) then
c return current foreground color
         call GetForeColor(color)
c set background color to best match for current device
         call RGBBackColor(black)
c set foreground color to best match for current device
         call RGBForeColor(white)
c restore colors
      elseif (isr.eq.-1) then
c set background color to best match for current device
         call RGBBackColor(black)
c set foreground color to best match for current device
         call RGBForeColor(color)
      endif
      return
      end
      subroutine auxmaus(event,ierr,estr)
c check alternative mouse menu events and replace with keyboard
c equivalents if appropriate
c input arguments: event
c event = event record
c ierr = return code for extras menu selection (-1 for non-menu event)
c estr = returned string for keyboard equivalent of menu selection
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c wrect = window rectangle
      integer*2 event(8)
      character*(*) estr
      integer*2 FindWindow
      integer MenuSelect, GetMenuHandle
      external FindWindow, MenuSelect, GetMenuHandle
      integer*2 intw(2), rect(4), mouseloc
      integer longw, loc, cptr, menuinfo, menuh
      character*256 str
      equivalence(longw,intw)
      intw(1) = event(6)
      intw(2) = event(7)
c loc = position in global coordinates where mouse event took place
      loc = longw
      ierr = 0
c see which window part, including menu bar, is at a point
      mouseloc = FindWindow(val4(loc),cptr)
c check if mouse is in drag region
      if (mouseloc.eq.4) then
         rect(1) = 0
         rect(2) = 0
         rect(3) = 2*(wrect(3) - wrect(1))
         rect(4) = 2*(wrect(4) - wrect(2))
c track the mouse and move a window
         call DragWindow(val4(cptr),val4(loc),rect)
         ierr = -1
c check if mouse is in menu bar
      elseif (mouseloc.eq.1) then
c initiate user selection of a menu item
         menuinfo = MenuSelect(val4(loc))
         longw = menuinfo
c given a menu ID, obtain a handle to the menu
         menuh = GetMenuHandle(val2(intw(1)))
         ierr = -1
         if (menuh.eq.0) return
c get the text of a menu item
         call GetMenuItemText(val4(menuh),val2(intw(2)),str)
c apple menu
         if (intw(1).eq.1) then
            continue
c file menu
         elseif (intw(1).eq.2) then
c quit
            if (str(2:5).eq.'Quit') then
               ierr = 1
               estr(1:1) = 'q'
c              call delwind
c              stop
c print a inverse copy of the current window
            elseif (str(2:21).eq.'Print Inverse Window') then
               call printit(0)
c print a copy of the current window
            elseif (str(2:13).eq.'Print Window') then
               call printit(1)
c print a copy of the current window
            elseif (str(2:15).eq.'Save Window...') then
               call saveit
            endif
c extras menu
         elseif (intw(1).eq.3) then
c obtain menu item
            menun = intw(2)
c check if pause was selected
            if (str(2:6).eq.'Pause') then
               ierr = 0
            elseif (menun.gt.0) then
               call mcextra(menun,ierr,estr)
            endif
         endif
c highlight or unhighlight a menu bar
      call HiliteMenu(val2(0))
      endif
      return
      end
      subroutine savscrn(isr,ipch)
c this subroutine saves a copy of the window to a Picture
c if isr = 1, current screen is stored into picture
c if isr = 0, current screen is stored with inverted colors
c if isr = -1, current screen is set to previous stored picture
c input argument: isr
c ipch = current picture structure handle
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c wrect = window rectangle
c wptr = pointer to cwindow structure
c crect = current clipping region
      logical*1 EmptyRect
      integer OpenCPicture, GetPortBitMapForCopyBits
      external OpenCPicture, EmptyRect, GetPortBitMapForCopyBits
      external GetWindowPort
      integer*2 picparam(12), srect(4), color(3)
      integer pich, out
      save pich, picparam
c pich = handle to picture structure
      data pich /0/
c set default OpenCPicParams structure
      data picparam /0,0,0,0,72,0,72,0,-2,0,0,0/
c source srect is full screen
      srect(1) = 0
      srect(2) = 0
      srect(3) = wrect(3) - wrect(1)
      srect(4) = wrect(4) - wrect(2)
c save screen
      if ((isr.eq.0).or.(isr.eq.1)) then
c release memory used by a picture definition
         if (pich.ne.0)  call KillPicture(val(pich))
c set OpenCPicParams structure
         do 10 i = 1, 4
            picparam(i) = srect(i)
   10    continue
c set clipping region to a rectangle
         call ClipRect(srect)
c copybits needs background white and foreground black
         call srcolor(1,color)
c Invert all pixels enclosed by a rectangle
         if (isr.eq.0) call InvertRect(srect)
c begin recording a picture definition
         pich = OpenCPicture(picparam)
         out = GetPortBitMapForCopyBits(val4(GetWindowPort(val4(wptr))))
c copy bitMap, with optional scaling, clipping, etc.
         call CopyBits(val4(out),val4(out),srect,srect,val2(64),val(0))
c stop recording picture information
         call ClosePicture
c lock a handle's data area (keep it from moving)
         call HLock(val4(pich))
         iptr = long(pich)
c determine if a rectangle is empty
         if (EmptyRect(val4(iptr+2))) then
c play a system alert sound
            call SysBeep(val2(6))
         endif
c restore background to black and foreground to previous color
         call srcolor(-1,color)
c set clipping region to a rectangle
         call ClipRect(crect)
c restore screen
      elseif (isr.eq.-1) then
c quit if no previous picture
         if (pich.eq.0) return
c activate a GrafPort
         call SetPort(val4(GetWindowPort(val4(wptr))))
c set clipping region to a rectangle
         call ClipRect(srect)
c Draw a pre-defined, scaled to desired size
         call DrawPicture(val(pich),srect)
c unlock a handle's data (allowing it to be moved)
         call HUnlock(val4(pich))
c release memory used by a picture definition
         call KillPicture(val(pich))
         pich = 0
c set clipping region to a rectangle
         call ClipRect(crect)
c keep a rectangular area from being updated
         call GetPortBounds(GetWindowPort(val4(wptr)),srect)
         call ValidWindowRect(val4(wptr),srect)
      endif
c return picture structure handle
      ipch = pich
      return
      end
      subroutine printit(isr)
c this subroutine prints a copy of the current window
c if isr = 1, current screen is stored into picture
c if isr = 0, current screen is stored with inverted colors
c external functions
      integer PMCreateSession, PMCreatePageFormat
      integer PMSessionDefaultPageFormat, PMSessionPageSetupDialog
      integer PMCreatePrintSettings, PMSessionDefaultPrintSettings
      integer PMSessionPrintDialog, PMSessionBeginDocument
      integer PMSessionError
      integer PMSessionBeginPage, PMSessionGetGraphicsContext
      integer PMSessionEndPage, PMSessionEndDocument
      external NewHandle, PMCreateSession, PMCreatePageFormat
      external PMSessionDefaultPageFormat, PMSessionPageSetupDialog
      external PMCreatePrintSettings, PMSessionDefaultPrintSettings
      external PMSessionPrintDialog, PMSessionBeginDocument
      external PMSessionError
      external PMSessionBeginPage, PMSessionGetGraphicsContext
      external PMSessionEndPage, PMSessionEndDocument
c local data
      integer*1 result
      integer*2 drect(4)
      integer kPMGraphicsContextQuickdraw  
      integer session, pformat, settings, oss, oldgp, grafptr
      save drect
c destination drect for 8.5 x 11 postscript printer
      data drect /0,0,405,540/
      data session, pformat, settings /0,0,0/
c allocate a printing session object and initialize it
      oss = PMCreateSession(session)
      if (oss.ne.0) return
c allocate a new page format object
      oss = PMCreatePageFormat(pformat)
      if (oss.ne.0) go to 50
c save screen to picture with inverted colors
      call savscrn(isr,ipch)
c assign default parameter values to an existing page format object
      oss = PMSessionDefaultPageFormat(val4(session),val4(pformat))
      if (oss.ne.0) go to 40
c display the page setup dialog and record the user's selections
c     oss = PMSessionPageSetupDialog(val4(session),val4(pformat),result)
c     if ((oss.ne.0).or.(result.eq.0)) go to 40
c allocate a new print settings object
      oss = PMCreatePrintSettings(settings)
      if (oss.ne.0) go to 40
c assign default parameter values to an existing print settings object
      oss = PMSessionDefaultPrintSettings(val4(session),val4(settings))
      if (oss.ne.0) go to 30
c display the print dialog and record the user's selections
      oss = PMSessionPrintDialog(val4(session),val4(settings),val4(pform
     1at),result)
      if ((oss.ne.0).or.(result.eq.0)) go to 30
c begin a print job
      oss = PMSessionBeginDocument(val4(session),val4(settings),val4(pfo
     1rmat))
      if (oss.ne.0) go to 30
c inform the printing system that the drawing is part of new page
      oss = PMSessionBeginPage(val4(session),val4(pformat),val4(0))
      if (oss.ne.0) go to 20
c find which grafPort is currently active
      call GetPort(oldgp)
      kPMGraphicsContextQuickdraw = CFStringCreateWithCharacters(val4(0)
     1,'com.apple.graphicscontext.quickdraw',val4(35))
      oss = PMSessionGetGraphicsContext(val4(session),val4(kPMGraphicsCo
     1ntextQuickdraw),grafptr)
      if (oss.ne.0) go to 10
c activate a grafPort
      call SetPort(val4(grafptr))
c Draw a pre-defined, scaled to desired size
      call DrawPicture(val(ipch),drect)
c activate a grafPort
      call SetPort(val4(oldgp))
c complete drawing of current page
   10 oss =  PMSessionEndPage(val4(session))
c end the print job
   20 oss = PMSessionEndDocument(val4(session))
c decrement the reference count for a printing object
   30 call PMRelease(val4(settings))
c decrement the reference count for a printing object
   40 call PMRelease(val4(pformat))
      if (isr.eq.0) then
c restore screen from picture (with inverted colors)
         call savscrn(-1,ipch)
c save screen to picture with inverted colors
         call savscrn(0,ipch)
      endif
c restore screen from picture (with normal colors)
      call savscrn(-1,ipch)
c decrement the reference count for a printing object
   50 call PMRelease(val4(session))
      return
      end
      subroutine saveit
c this subroutine saves a copy of the current window as PICT
      implicit none
c external functions
      integer*2 NavGetDefaultDialogOptions, NavPutFile, AEGetNthDesc
      integer*2 AEGetDescData, FSpDelete, FSpCreate, FSpOpenDF
      integer*2 FSWrite, FlushVol, FSClose, AEDisposeDesc
      integer*2 NavCompleteSave, NavDisposeReply
      integer GetHandleSize
      external NavGetDefaultDialogOptions,  NavPutFile, AEGetNthDesc
      external AEGetDescData, FSpDelete, FSpCreate, FSpOpenDF
      external FSWrite, FlushVol, FSClose, AEDisposeDesc
      external NavCompleteSave, NavDisposeReply, GetHandleSize
c local data
      integer kNavDontAutoTranslate, kNavDontAddTranslateItems
      parameter(kNavDontAutoTranslate=2,kNavDontAddTranslateItems=4)
      integer kNavTranslateInPlace, typeFSS
      parameter(kNavTranslateInPlace=0,typeFSS=1718842144)
      integer*2 ierr, frefn
      integer*2 doptions(1024), reply(128), filespec(35)
      integer i, ipch, ftype, fcreator, nerr, keyword, isize, iptr
      integer descr(2), header(128)
      character*16 prompt
      character*8 default
      data prompt, default /'Save document as','gks.pict'/
      data ftype, fcreator /'PICT','ttxt'/
c determine default attributes or behavior for dialog boxes
      ierr =  NavGetDefaultDialogOptions(doptions)
      if (ierr.ne.0) return
c save screen to picture
      call savscrn(1,ipch)
c set specific options for dialog box
      doptions(3) = kNavDontAutoTranslate + kNavDontAddTranslateItems
      doptions(646) = 256*16 + ichar(prompt(1:1))
      do 10 i = 1, 7
      doptions(646+i) = 256*ichar(prompt(2*i:2*i)) + ichar(prompt(2*i+1:
     12*i+1)) 
   10 continue
      doptions(654) = 256*ichar(prompt(16:16)) 
      doptions(518) = 256*8 + ichar(default(1:1))
      do 20 i = 1, 3
      doptions(518+i) = 256*ichar(default(2*i:2*i)) + ichar(default(2*i+
     11:2*i+1)) 
   20 continue
      doptions(522) = 256*ichar(default(8:8)) 
c display a save dialog
      ierr = NavPutFile(val4(0),reply,doptions,val4(0),val4(ftype),val4(
     1fcreator),val4(0))
      if (ierr.ne.0) go to 60
      nerr = reply(2)/256
c check if ValidRecord is false
      if (nerr.eq.0) go to 50
c copy a descriptor record from a descriptor list
      ierr = AEGetNthDesc(reply(4),val4(1),val4(typeFSS),keyword,descr)
      if (ierr.ne.0) go to 50
c get data from specified descriptor record
      ierr = AEGetDescData(descr,filespec,val4(70))
      if (ierr.ne.0) go to 40
c check if replacing flag is true
      nerr = reply(2) - 256*nerr
      if (nerr.eq.1) then
c remove a closed file
         ierr = FSpDelete(filespec)
      endif
c create a new file and set the type and creator
      ierr = FSpCreate(filespec,val4(fcreator),val4(ftype),val4(-1))
      if (ierr.ne.0) go to 40
c create an access path to the data fork of a file
      ierr = FSpOpenDF(filespec,val1(3),frefn)
      if (ierr.ne.0) go to 40
      do 30 i = 1, 128
      header(i) = 0
   30 continue
      isize = 512
c write data from memory to a file
      ierr = FSWrite(val2(frefn),isize,header)
      if (ierr.ne.0) go to 40
c get the size of a handle's data area
      isize = GetHandleSize(val4(ipch))
      iptr = long(ipch)
c write data from memory to a file
      ierr = FSWrite(val2(frefn),isize,val4(iptr))
      if (ierr.ne.0) go to 40
c update disk with any unwritten data
      ierr = FlushVol(val4(0),val2(filespec(1)))
      if (ierr.ne.0) go to 40
c close a file
      ierr = FSClose(val2(frefn))
      if (ierr.ne.0) go to 40
c complete a save operation and perform any needed translation
      ierr = NavCompleteSave(reply,val4(kNavTranslateInPlace))
c deallocate memory used by a descriptor record
   40 ierr = AEDisposeDesc(descr)
c release memory allocated for a NavReplyRecord structure
   50 ierr = NavDisposeReply(reply)
c restore screen from picture
   60 call savscrn(-1,ipch)
      return
      end
      subroutine delwind
c this subroutine closes a window
c input arguments: none
      integer wptr, cthandle
      integer*2 wrect(4), crect(4)
      common /gmcdev/ wrect, crect, wptr, cthandle
c wptr = pointer to cwindow structure
      integer GetMenuHandle
      external GetMenuHandle
      integer*2 menuid
      character*1 str
c unlock a handle's data area (allowing it to be moved)
      call HUnlock(val4(cthandle))
c free allocation created via NewHandle of color table
      call DisposeHandle(val4(cthandle))
c destroy extra menus
      call mcextra(-1,ierr,str)
c remove menus
      do 10 i = 1, 2
      menuid = i
c remove a menu from the menulist
      call DeleteMenu(val2(menuid))
c given a menu ID, obtain a handle to the menu
      menuh = GetMenuHandle(val2(menuid))
c release memory menu created via NewMenu
      call DisposeMenu(val4(menuh))
   10 continue
c remove window from screen; keep WindowRecord
      call DisposeWindow(val4(wptr))
      return
      end
      subroutine dashln(i,j,lnt)
c this subroutine draws a dashed line from (icx,icy) to (i,j)
c input arguments: all
c i, j = location of point
c lnt = line type
c 1 = solid, 2 = short-dash, 3 = dot, 4 = dash-dot, 5 = long-dash
      dimension lt(4,4)
      integer*2 i, j
      integer*2 point(2), icx, icy, ix, iy
      save ns,nl,ltype
c ltype = current line type
      data ltype /0/
c lt = software dash type patterns
      data lt /9,6,9,6,5,5,5,5,14,6,4,6,23,7,23,7/
      l = lnt - 1
c use current style if line style is known
      if ((l.lt.0).or.(l.gt.7)) l = ltype
c reset line style, if necessary
      if (l.ne.ltype) then
         ltype = l
c reset software line style
         if ((ltype.ge.1).and.(ltype.le.4)) then
         ns = 0
         nl = lt(ns+1,ltype)
         endif
      endif
c use solid style if requested or unknown
      if ((ltype.eq.0).or.(ltype.gt.4)) then
c draw a line to specified coordinates
         call LineTo(val2(i),val2(j))
         return
      endif
c obtain current pen position
      call GetPen(point)
      icx = point(2)
      icy = point(1)
c software dashed line
      cost = float(i - icx)
      sint = float(j - icy)
      alen = sqrt(cost*cost + sint*sint)
      len = alen + .5
c current pattern segment is longer than line length
      if (nl.ge.len) go to 20
c find starting location and direction cosines
      ax0 = float(icx) + .5
      ay0 = float(icy) + .5
      cost = cost/alen
      sint = sint/alen
c iterate pattern segments
   10 anl = float(nl)
c find end coordinate of next segment
      ix = ax0 + anl*cost
      iy = ay0 + anl*sint
c dark or bright vector flag
      it = ns - (ns/2)*2
c draw line
      if (it.eq.0) then
c draw a line to specified coordinates
         call LineTo(val2(ix),val2(iy))
c set pen location without drawing
      else
         call MoveTo(val2(ix),val2(iy))
      endif
c increment pattern segment index
      ns = ns + 1
      if (ns.eq.4) ns = 0
c add length of next pattern segment
      nl = nl + lt(ns+1,ltype)
c do next segment
      if (nl.lt.len) go to 10
c finish up last segment, which may be incomplete
   20 it = ns - (ns/2)*2
c draw line
      if (it.eq.0) then
c draw a line to specified coordinates
         call LineTo(val2(i),val2(j))
c set pen location without drawing
      else
         call MoveTo(val2(i),val2(j))
      endif
c adjust length of next pattern segment
      nl = nl - len
c if segment complete, reset to next pattern segment
      if (nl.eq.0) then
         ns = ns + 1
         if (ns.eq.4) ns = 0
         nl = lt(ns+1,ltype)
      endif
      return
      end
      subroutine dimagx(image,lx,ly,lz,lenb,npix,nf)
c this subroutine displays raster image stored in character array
c if necessary, data is copied to a character array
c an identity transformation has been assumed
c image = uncompressed single image
c lx, ly = the size of the image, in pixels
c lz = width of picture, in bytes
c lenb = size of picture, in bytes
c npix = number of pixels per byte
c nf = current frame number being processed
c optimized for macintosh graphics
c npald = number of palette entries
      parameter(npald=256)
c lxm, lym = maximum number of pixels in x, y
      parameter(lxm=720,lym=540)
c idwk = workstation identifier
c rx, ry = ndc coordinates of upper-right corner of workstation window
c ifrg = index of foreground color
c isx, isy = display width, height, in raster units
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c lupt = (0,1) = (no,yes) pixel lookup table needed
      common /movicm/ lupt,ipal,img8
      character*1 image(lenb)
c ipal = integer pixel lookup table
      dimension ipal(npald)
c img8 = character*1 image array
      character*1 img8(lxm*lym)
      character*8 lbl
      integer*2 rect(4), ix, iy, my, it
      save nfl,n,csize,alx,lbl,img8
c nfl = previous frame number
      data nfl /1/
c csize = vertical size of characters
c alx = x coordinate for numerical labels
      data csize,alx /.02,.86/
c lbl = character variable for numerical label
      data lbl /'#       '/
c nbit = the number of colors, pixel depth
      nbit = 8/npix
c calculate maximum size of image
      lxs = lx
      if (lxs.gt.isx) lxs = isx
      lys = ly
      if (lys.gt.isy) lys = isy
      if (lxs.gt.lxm) lxs = lxm
      if (lys.gt.lym) lys = lym
c normalize characters and label location
      ix = alx*float(isx)
      iy = 0
      it = csize*float(isy)
c eight bit color
      if (nbit.eq.8) then
c display raster without modification
         if (lupt.eq.0) then
            do 20 k = 1, lys
            ioff = lz*(k - 1)
            joff = lxs*(k - 1)
            do 10 j = 1, lxs
            img8(j+joff) = image(j+ioff)
   10       continue
   20       continue
c copy to character array with lookup table
         else
            do 40 k = 1, lys
            ioff = lz*(k - 1)
            joff = lxs*(k - 1)
            do 30 j = 1, lxs
            img8(j+joff) = char(ipal(ichar(image(j+ioff))+1))
   30       continue
   40       continue
         endif
c nbits per pixel
      else
c maximum width
         lzs = lxs/npix
c convert from nbits per pixel to 8 bits per pixel
         ntc = 2**nbit
         npixm = npix - 1
c copy to character array without lookup table
         if (lupt.eq.0) then
            do 70 k = 1, lys
c loop over bytes in image
            j1 = lxs*(k - 1) + npix + 1
            ioff = lz*(k - 1)
            do 60 j = 1, lzs
            itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
            do 50 i = 1, npixm
            it1 = itc/ntc
            img8(j1-i) = char(itc-it1*ntc)
            itc = it1
   50       continue
            img8(j1-npix) = char(itc)
            j1 = j1 + npix
   60       continue
   70       continue
c copy to character array with lookup table
         else
            do 100 k = 1, lys
c loop over bytes in image
            j1 = lxs*(k - 1) + npix + 1
            ioff = lz*(k - 1)
            do 90 j = 1, lzs
            itc = ichar(image(j+ioff))
c innermost loop over pixels in byte
            do 80 i = 1, npixm
            it1 = itc/ntc
            img8(j1-i) = char(ipal(itc-it1*ntc+1))
            itc = it1
   80       continue
            img8(j1-npix) = char(ipal(itc+1))
            j1 = j1 + npix
   90       continue
  100       continue
         endif
      endif
c my = y screen offset, used because origin is at top of screen
      my = isy - 1
      rect(1) = isy - lys
      rect(2) = 0
      rect(3) = isy
      rect(4) = lxs
c clear workstation, if frame number is not advanced
      if (nf.le.nfl) then
         call gclrwk(idwk,1)
c otherwise, erase label
      else
c set point size for subsequent text drawing
         call TextSize(val2(it))
c select current color
         call scolor(0)
c set pen location without drawing
         call MoveTo(val2(ix),val2(my - iy))
c draw text from any arbitrary buffer
         call DrawText(lbl,val2(0),val2(8))
      endif
c copy offscreen image to screen
      call copyimg(img8,rect,8,lxs,lys)
c first find how many digits in nf
      id = 0
      n = nf
  110 id = id + 1
      n = n/10
      if (n.gt.0) go to 110
c create label template
      lbl = '#       '
c create left justfied label
      is = ichar('0')
      if (id.gt.7) id = 7
      ls = 10**(id - 1)
      nt = nf
      do 120 i = 1, id
      i1 = i + 1
      n = nt/ls
      lbl(i1:i1) = char(n+is)
      nt = nt - n*ls
      ls = ls/10
  120 continue
c save current frame number
      nfl = nf
c set point size for subsequent text drawing
      call TextSize(val2(it))
c select current color
      call scolor(ifrg)
c set length of string
      n = id + 1
c set pen location without drawing
      call MoveTo(val2(ix),val2(my - iy))
c draw text from any arbitrary buffer
      call DrawText(lbl,val2(0),val2(n))
      return
      end
      integer function kywait()
c special function to request keystroke
c returns: ascii code for keystroke
      integer*2 myEventMask, mouseDown, keyDown, autoKey, osEvt, resus
c myEventMask looks for mouse, keyboard, system and quit events
      parameter(myEventMask=1086-32768,mouseDown=1,keyDown=3,autoKey=5)
      parameter(osEvt=15,resus=256)
      integer*2 event(8)
      character*1 str
c multifinder-aware way to obtain events
   10 call WaitNextEvent(val2(myEventMask),event,val4(-1),val4(0))
      if ((event(1).eq.osEvt).and.(iand(event(2),resus).eq.resus)) then
c resume event
         if ((event(3)-2*(event(3)/2)).eq.1) then
            call savscrn(-1,ipch)
c suspend event
         else
            call savscrn(1,ipch)
         endif
c try again
         go to 10
      endif
      if (event(1).eq.mouseDown) then
c check alternative mouse events
         call auxmaus(event,ierr,str)
c get another event
         go to 10
      endif
c try again if not keyboard
      if ((event(1).ne.keyDown).and.(event(1).ne.autoKey)) go to 10
      kywait = event(3) - 256*(event(3)/256)
      return
      end
c internal mactintosh library for enhancing libgks1 images
      subroutine mcextra(kreate,ierr,str)
c this subroutine creates, destroys, or processes extras menus used by
c libgks1 images
c input arguments: kcreate
c kreate = (-1,0,1,2,...) = (destroy,create,process) extras menus
c ierr = return code for extras menu selection (-1 for non-menu event)
c str = returned string for keyboard equivalent of menu selection
      character*(*) str
      integer NewMenu, GetMenuHandle
      external NewMenu, GetMenuHandle
      integer*2 menuid
      save menuh
      data menuh /0/
c menu 3 = extras menu
      menuid = 3
      ierr = -1
c create extra menus
      if (kreate.eq.0) then
c get a handle to an empty menu
         menuh = NewMenu(val2(menuid),char(6)//'Extras')
c add one or more items to a menu
         call AppendMenu(val4(menuh),char(11)//'Save & Quit')
c add one or more items to a menu
         call AppendMenu(val4(menuh),char(6)//'Modify')
c add one or more items to a menu
         call AppendMenu(val4(menuh),char(7)//'Animate')
c add one or more items to a menu
         call AppendMenu(val4(menuh),char(5)//'Reset')
c add one or more items to a menu
         call AppendMenu(val4(menuh),char(12)//'Plotparms...')
c add one or more items to a menu
         call AppendMenu(val4(menuh),char(12)//'Set Frame...')
c add a menu to the menu list
         call InsertMenu(val4(menuh),val2(0))
c display the titles of all the menus in the menu list
         call DrawMenuBar
c delete extra menus
      elseif (kreate.eq.-1) then
         if (menuh.eq.0) return
c remove a menu from the menulist
         call DeleteMenu(val2(menuid))
c given a menu ID, obtain a handle to the menu
         menuh = GetMenuHandle(val2(menuid))
c release memory menu created via NewMenu
         call DisposeMenu(val4(menuh))
         menuh = 0
c process save selection
      elseif (kreate.eq.1) then
         ierr = 2
         str(1:1) = 's'
c process modify selection
      elseif (kreate.eq.2) then
         ierr = 3
         str(1:1) = 'm'
c process animation selection
      elseif (kreate.eq.3) then
         ierr = 4
         str(1:1) = 'a'
c process reset selection
      elseif (kreate.eq.4) then
         ierr = 6
         str(1:1) = 'r'
c process plotparms selection
      elseif (kreate.eq.5) then
c        ierr = 7
c        str(1:1) = 'p'
         call plparm
         ierr = -1
c process frame number selection
      elseif (kreate.eq.6) then
c get user input
         call idialog('Enter Frame Number: ',str,ierr)
         if (ierr.gt.0) ierr = 128 + ierr
      endif
      return
      end
      subroutine idialog(chars,str,ierr)
c this subroutine creates a dialog box and asks for integer input
c input argument: chars
c chars = static prompt text for dialog box
c str = output string representing number
c ierr = error code (-1,string length) for (invalid,valid) string
      character*(*) chars, str
      integer*2 TextWidth
      integer NewHandle, NewDialog
      external TextWidth, NewHandle, NewDialog
      integer ctrlItem, statText, editText, itemDisable, ok, cancel
      parameter(ctrlItem=4,statText=8,editText=16,itemDisable=128)
      parameter(ok=1,cancel=2)
      character*12 cnum
      integer*2 item, itemt, if(4), rect(4)
      integer dptr, dlisth
c horizontal spacing between items
      data lsp /16/
c set size of data items (number of 2 byte items)
      data lok,lcancel /1,3/
c save screen to picture
      call savscrn(1,ipch)
      nh = len(chars)/2
c make sure n is even, truncate if necessary
      n = 2*nh
c save current font size
c obtain font sizing information
      call GetFontInfo(if)
c set point size for subsequent text drawing
      ich = 12
      call TextSize(val2(ich))
c get width of unformatted text
      icw = TextWidth(chars,val2(0),val2(n))
c get width of unformatted text
      icn = TextWidth('0123456789',val2(0),val2(10)) + 12
c get width of unformatted text
      icc = TextWidth('Cancel',val2(0),val2(6)) + 12
      icnw = max(icn,icw)
      nitems = 4
c allocate relocatable block from current heap zone
      dlisth = NewHandle(val(2+14*nitems+n+2*(lok+lcancel)))
      dptr = long(dlisth)
c set size and location for OK ctrlItem
      rect(1) = 2.5*ich
      rect(2) = 2*lsp + icnw
      rect(3) = rect(1) + 1.25*ich
      rect(4) = rect(2) + icc
      ldptr = dptr + 2
c add OK ctrlItem to item list
      call citlist(val4(ldptr),rect,ctrlItem,'OK',lok)
c set size and location for Cancel ctrlItem
      rect(1) = .5*ich
      rect(3) = rect(1) + 1.25*ich
      ldptr = ldptr + 14 + 2*lok
c add Cancel ctrlItem to item list
      call citlist(val4(ldptr),rect,ctrlItem,'Cancel',lcancel)
c set size for editText item
      rect(1) = 2.5*ich
      rect(2) = lsp + .5*(icnw - icn)
      rect(3) = rect(1) + 1.25*ich
      rect(4) = rect(2) + icn
      ldptr = ldptr + 14 + 2*lcancel
c add editText to item list
      call citlist(val4(ldptr),rect,editText+itemDisable,'00',0)
c set size for statText item
      rect(1) = .5*ich
      rect(2) = lsp + .5*(icnw - icw)
      rect(3) = rect(1) + ich
      rect(4) = rect(2) + icw
      ldptr = ldptr + 14
c add disabled statText to item list
      call citlist(val4(ldptr),rect,statText+itemDisable,chars,nh)
c write header for item list
      word(dptr) = nitems - 1
c place dialog box in upper left corner
      rect(1) = 10
      rect(2) = 10
      rect(3) = rect(1) + 4.5*ich
      rect(4) = rect(2) + (icnw + icc + 3*lsp)
c obtain global (screen) value of local point
      call LocalToGlobal(rect)
      call LocalToGlobal(rect(3))
c create a new dialog box
      dptr = NewDialog(val4(0),rect,char(0),val1(.true.),val2(1),val4(-1
     1),val1(.false.),val4(0),val4(dlisth))
c begin user interaction in a modal dialog
      call ModalDialog(val4(0),item)
c ok
      if (item.eq.1) then
c obtain dialog item type, handle, and rectangle
         call GetDialogItem(val4(dptr),val2(3),itemt,ith,irp)
c obtain a copy of the text of an editText item
         call GetDialogItemText(val4(ith),cnum)
c check for invalid characters
         call ckchar(cnum,str,ierr)
         if (ierr.eq.0) ierr = -1
c cancel
      elseif (item.eq.2) then
         ierr = -1
      endif
c free allocation created via NewHandle
      call DisposeHandle(dlisth)
c close dialog and release all related memory
      call DisposeDialog(val4(dptr))
c set point size for subsequent text drawing
      call TextSize(val2(if(1)))
c restore screen from picture
      call savscrn(-1,ipch)
      return
      end
      subroutine plparm
c this subroutine allows one to make changes to certain of the plotting
c parameters stored in the common blocks plotcm and pextra.
c currently, only the variables nplot and ndi can be changed.
c input arguments: none
      integer*2 TextWidth
      integer NewHandle, NewDialog
      external TextWidth, NewHandle, NewDialog
      integer ctrlItem, statText, editText, itemDisable, ok, cancel
      parameter(ctrlItem=4,statText=8,editText=16,itemDisable=128)
      parameter(ok=1,cancel=2)
c nvars = number of variables to be modified
      parameter (nvars=2)
c nplot = number of plots per page
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c ndi = increment between frames
      common /pextra/ ndi
      dimension nvalue(nvars)
      character*8 code(nvars)
      character*12 output, cnum
      character*32 prompt, helpv(nvars)
      integer*2 item, itemt, if(4), rect(4)
      integer dptr, dlisth
      save prompt,code,helpv
c horizontal spaceing between items
      data lsp /16/
c set size of data items (number of 2 byte items)
      data lok,lcancel /1,3/
      data prompt /'select variable name for help   '/
      data code /'NPLOT = ','NDI   = '/
      data helpv  /'NPLOT = number of plots per page','NDI = increment b
     1etween frames  '/   
c save screen to picture
      call savscrn(1,ipch)
c first copy variables
      nvalue(1) = nplot
      nvalue(2) = ndi
c save current font size
c obtain font sizing information
      call GetFontInfo(if)
c set point size for subsequent text drawing
      ich = 12
      call TextSize(val2(ich))
c get width of unformatted text
      icw = TextWidth('NPLOTX= ',val2(0),val2(8)) + 12
c get width of unformatted text
      icc = TextWidth('Cancel',val2(0),val2(6)) + 12
      lvars = lsp + 2*icw
      nitems = 2*nvars + 3
c allocate relocatable block from current heap zone
      dlisth = NewHandle(val(2+14*nitems+16*nvars+2*(lok+lcancel)+32))
      dptr = long(dlisth)
c set size and location for OK ctrlItem
      rect(1) = 4.5*ich
      rect(2) = nvars*lvars - icc
      rect(3) = rect(1) + 1.25*ich
      rect(4) = rect(2) + icc
      ldptr = dptr + 2
c add OK ctrlItem to item list
      call citlist(val4(ldptr),rect,ctrlItem,'OK',lok)
      ldptr = ldptr + 14 + 2*lok
c set size and location for Cancel ctrlItem
      rect(2) = (nvars-1)*lvars + 4
      rect(4) = rect(2) + icc
c add Cancel ctrlItem to item list
      call citlist(val4(ldptr),rect,ctrlItem,'Cancel',lcancel)
      ldptr = ldptr + 14 + 2*lcancel
c display variables
      rect(1) = 0.5*ich
      rect(2) = 0
      rect(3) = rect(1) + 1.25*ich
      rect(4) = 0
      do 10 j = 1, nvars
c set size for editText item
      rect(2) = rect(4) + lsp
      rect(4) = rect(2) + icw
c add editText to item list
      call citlist(val4(ldptr),rect,statText,code(j),4)
      ldptr = ldptr + 22
c set size for editText item
      rect(2) = rect(4)
      rect(4) = rect(2) + icw
c convert 32-bit integer to string of decimal digits
      call NumToString(val4(nvalue(j)),output)
c pad with blanks
      ls = ichar(output(1:1)) + 2
      if (ls.le.9) output(ls:9) = '        '
c add editText to item list
      call citlist(val4(ldptr),rect,editText+itemDisable,output(2:9),4)
      ldptr = ldptr + 22
   10 continue
c leave room for help
c set size for statText item
      rect(1) = 2.5*ich
      rect(2) = lsp
      rect(3) = rect(1) + ich
      rect(4) = rect(2) + (2*lvars - lsp)
c add disabled statText to item list
      call citlist(val4(ldptr),rect,statText+itemDisable,prompt,16)
c write header for item list
      word(dptr) = nitems - 1
c place dialog box in upper left corner
      rect(1) = 10
      rect(2) = 10
      rect(3) = rect(1) + 6.5*ich
      rect(4) = rect(2) + (2*lvars + lsp)
c obtain global (screen) value of local point
      call LocalToGlobal(rect)
      call LocalToGlobal(rect(3))
c create a new dialog box
      dptr = NewDialog(val4(0),rect,char(0),val1(.true.),val2(1),val4(-1
     1),val1(.false.),val4(0),val4(dlisth))
c draw the contents of a dialog box
c     call DrawDialog(val4(dptr))
c     call kywait()
c begin user interaction in a modal dialog
   20 call ModalDialog(val4(0),item)
c ok
      if (item.eq.1) then
         nerrs = 0
         do 30 j = 1, nvars
c obtain dialog item type, handle, and rectangle
         call GetDialogItem(val4(dptr),val2(2*(j+1)),itemt,ith,irp)
c obtain a copy of the text of an editText item
         call GetDialogItemText(val4(ith),cnum)
c check for invalid characters
         call ckchar(cnum,output,ierr)
c valid characters found
         if (ierr.gt.0) then
c convert string of decimal digits to binary number
            call StringToNum(cnum,nvalue(j))
c count number of invalid inputs
         else
            nvalue(j) = 0
            nerrs = nerrs + 1
         endif
   30    continue
c invalid input
         if (nerrs.gt.0) then
c play a system alert sound
            call SysBeep(val2(6))
            go to 20
         else
            nplot = nvalue(1)
            ndi = nvalue(2)
         endif
      endif
c check if help selected
      do 40 j = 1, nvars
      if (item.eq.(2*j+1)) then
c obtain dialog item type, handle, and rectangle
         call GetDialogItem(val4(dptr),val2(nitems),itemt,ith,irp)
c specify the text of an editText item and draw it
         call SetDialogItemText(val4(ith),char(32)//helpv(j))
         go to 20
      endif
   40 continue
c free allocation created via NewHandle
      call DisposeHandle(dlisth)
c close dialog and release all related memory
      call DisposeDialog(val4(dptr))
c set point size for subsequent text drawing
      call TextSize(val2(if(1)))
c restore screen from picture
      call savscrn(-1,ipch)
      return
      end
      subroutine citlist(itlist,rect,itype,idata,nd)
c add item to item list for dialogs, described in Inside Macintosh,
c vol. I, page 427.
c input arguments: rect, itype, data
c itlist = output item list
c rect = rect for displaying item in dialog box
c itype = item type
c idata = item data
c nd = length of item data
      integer*2 itlist(nd+7), rect(4), idata(*)
c handle placeholder
      itlist(1) = 0
      itlist(2) = 0
c item rect
      do 10 i = 1, 4
      itlist(2+i) = rect(i)
   10 continue
c set type
      if (itype.ge.128) then
         itlist(7) = 256*(itype - 128) - 32768
      else
         itlist(7) = 256*itype
      endif
c set length
      itlist(7) = itlist(7) + (nd + nd)
c copy item data
      do 20 i = 1, nd
      itlist(7+i) = idata(i)
   20 continue   
      return
      end
      subroutine ckchar(cnum,str,lstr)
c this subroutine checks for non-numeric characters in pascal string cnum
c the integer part is returned in a fortran string whose length is lstr
c input cnum is also modified if partially valid string is found
c input argument: cnum
c cnum = pascal style string containing integers and possible errors
c str = fortran style with trailing non-integers removed
c lstr = number of digits in integer number (0 if no valid digits)
      character*(*) cnum, str
c get length of string
      lstr = ichar(cnum(1:1))
c check if characters represent a number
      if (lstr.gt.0) then
         i = 1
         num = 0
   10    iv = ichar(cnum(i+1:i+1)) - ichar('0')
c abort if not a number
         if ((iv.lt.0).or.(iv.gt.9)) go to 20
         num = iv + 10*num
         i = i + 1
         if (i.le.lstr) go to 10
c invalid string
   20    if (i.eq.1) then
            lstr = 0
c string truncated due to invalid character
         elseif ((i.gt.1).and.(i.le.lstr)) then
c convert 32-bit integer to string of decimal digits
            call NumToString(val4(num),cnum)
c get new length of string
            lstr = ichar(cnum(1:1))
c convert from pascal string to fortran string
            str(1:lstr) = cnum(2:lstr+1)
c no errors
         else
            str(1:lstr) = cnum(2:lstr+1)
         endif
c null string
      else
         lstr = 0
      endif
      return
      end
