c-----------------------------------------------------------------------
c 2d parallel PIC multi-tasking library for fast fourier transforms
c mpfft2lib.f contains multi-tasking procedures to perform ffts:
c MPFFT2R multi-tasking wrapper for WPFFT2R
c MPFFT2RX multi-tasking wrapper for WPFFT2RX
c MPFFT2R2 multi-tasking wrapper for WPFFT2R2
c MPFFT2R3 multi-tasking wrapper for WPFFT2R3
c MPFFT2RX2 multi-tasking wrapper for WPFFT2RX2
c MPFFT2RX3 multi-tasking wrapper for WPFFT2RX3
c MPFFT2RN multi-tasking wrapper for WPFFT2RN
c MPFFT2RXN multi-tasking wrapper for WPFFT2RXN
c MP2FFT2RN multi-tasking wrapper for WP2FFT2RN
c MPFSST2R multi-tasking wrapper for WPFSST2R
c MPFSCT2R multi-tasking wrapper for WPFSCT2R
c MPFCST2R multi-tasking wrapper for WPFCST2R
c MPFCCT2R multi-tasking wrapper for WPFCCT2R
c MPFCST2R2 multi-tasking wrapper for WPFCST2R2
c MPFSCT2R2 multi-tasking wrapper for WPFSCT2R2
c MPFCST2R3 multi-tasking wrapper for WPFCST2R3
c MPFSCT2R3 multi-tasking wrapper for WPFSCT2R3
c MPFSFT2R multi-tasking wrapper for WPFSFT2R
c MPFCFT2R multi-tasking wrapper for WPFCFT2R
c MPFCSFT2R2 multi-tasking wrapper for WPFCSFT2R2
c MPFSCFT2R2 multi-tasking wrapper for WPFSCFT2R2
c MPFCSFT2R3 multi-tasking wrapper for WPFCSFT2R3
c MPFSCFT2R3 multi-tasking wrapper for WPFSCFT2R3
c MPFDSFT2RX multi-tasking wrapper for WPFDSFT2RX
c MPFDCFT2RX multi-tasking wrapper for WPFDCFT2RX
c MPFDCSFT2R2 multi-tasking wrapper for WPFDCSFT2R2
c MPFDSCFT2R2 multi-tasking wrapper for WPFDSCFT2R2
c MPFDCSFT2R3 multi-tasking wrapper for WPFDCSFT2R3
c MPFDSCFT2R3 multi-tasking wrapper for WPFDSCFT2R3
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: december 10, 2008
c-----------------------------------------------------------------------
      subroutine MPFFT2R(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,
     1kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask,n
     2mt,ierr)
c multi-tasking real to complex fft
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      integer kxyip, iftask, nmt, ierr
      real ttp
      complex f, g, bs, br, sct
      dimension f(nxvh,kypd,kblok), g(nyv,kxp,jblok)
      dimension bs(kxp,kyp,kblok), br(kxp,kyp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFFT2RXX, PFFT2RXY
      data nargs /14/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      kxpp = kxp/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RXX,nargs,f,isign,mixup,sct,in
     1dx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jb
     1lok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RXY,nargs,g,isign,mixup,sct,in
     1dx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp
     1,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd
     1,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RXY,nargs,g,isign,mixup,sct,in
     1dx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kb
     1lok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RXX,nargs,f,isign,mixup,sct,in
     1dx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call PFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT2RX(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,kstrt
     1,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask,nmt,ie
     2rr)
c multi-tasking real to complex fft
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh
      integer nyv, kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
      integer kxyip, iftask, nmt, ierr
      real ttp
      complex f, g, sct
      dimension f(nxvh,kypd,kblok), g(nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFFT2RXX, PFFT2RXY
      data nargs /14/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      kxpp = kxp/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RXX,nargs,f,isign,mixup,sct,in
     1dx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,k
     1blok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RXY,nargs,g,isign,mixup,sct,in
     1dx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblo
     1k,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblo
     1k,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RXY,nargs,g,isign,mixup,sct,in
     1dx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFFT2RXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,k
     1xp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,j
     1blok)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RXX,nargs,f,isign,mixup,sct,in
     1dx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call PFFT2RXX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh,
     1kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT2R2(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy
     1,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask,
     2nmt,ierr)
c multi-tasking real to complex fft
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      integer kxyip, iftask, nmt, ierr
      real ttp
      complex f, g, bs, br, sct
      dimension f(2,nxvh,kypd,kblok), g(2,nyv,kxp,jblok)
      dimension bs(2,kxp,kyp,kblok), br(2,kxp,kyp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFFT2R2XX, PFFT2R2XY
      data nargs /14/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      kxpp = kxp/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R2XX,nargs,f,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh
     1,kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P2TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,j
     1blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R2XY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P2TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kx
     1p,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P2TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kyp
     1d,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R2XY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P2TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,k
     1blok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R2XX,nargs,f,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call PFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh
     1,kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT2R3(f,g,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy
     1,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask,
     2nmt,ierr)
c multi-tasking real to complex fft
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      integer kxyip, iftask, nmt, ierr
      real ttp
      complex f, g, bs, br, sct
      dimension f(3,nxvh,kypd,kblok), g(3,nyv,kxp,jblok)
      dimension bs(3,kxp,kyp,kblok), br(3,kxp,kyp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFFT2R3XX, PFFT2R3XY
      data nargs /14/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      kxpp = kxp/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R3XX,nargs,f,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh
     1,kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,j
     1blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R3XY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kx
     1p,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kyp
     1d,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R3XY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,k
     1blok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R3XX,nargs,f,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call PFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh
     1,kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT2RX2(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,kstr
     1t,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask,nmt,i
     2err)
c multi-tasking real to complex fft
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh
      integer nyv, kxp, kyp, kypd, jblok, kblok, nxhyd, nxyhd
      integer kxyip, iftask, nmt, ierr
      real ttp
      complex f, g, sct
      dimension f(2,nxvh,kypd,kblok), g(2,nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFFT2R2XX, PFFT2R2XY
      data nargs /14/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      kxpp = kxp/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R2XX,nargs,f,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh
     1,kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P2TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,
     1kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R2XY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P2TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kbl
     1ok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P2TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jbl
     1ok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R2XY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFFT2R2XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P2TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,
     1jblok)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R2XX,nargs,f,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call PFFT2R2XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh
     1,kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT2RX3(f,g,isign,ntpose,mixup,sct,ttp,indx,indy,kstr
     1t,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,nxhyd,nxyhd,kxyip,iftask,nmt,i
     2err)
c multi-tasking real to complex fft
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, nxhyd, nxyhd, mixup
      integer kxyip, iftask, nmt, ierr
      real ttp
      complex f, g, sct
      dimension f(3,nxvh,kypd,kblok), g(3,nyv,kxp,jblok)
      dimension mixup(nxhyd), sct(nxyhd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFFT2R3XX, PFFT2R3XY
      data nargs /14/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      kxpp = kxp/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R3XX,nargs,f,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh
     1,kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,
     1kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R3XY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kbl
     1ok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jbl
     1ok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R3XY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFFT2R3XY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,
     1jblok)
         call PWTIMERA(-1,ttp,dtime)
c start x fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2R3XX,nargs,f,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call PFFT2R3XX(f,isign,mixup,sct,indx,indy,kstrt,kypi,kypl,nxvh
     1,kypd,kblok,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT2RN(f,g,bs,br,ss,isign,ntpose,mixup,sct,ttp,indx,i
     1ndy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd,kxyip
     2,iftask,nmt,ierr)
c multi-tasking real to complex fft
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, ndim, nxhyd, nxyhd, mixup
      integer kxyip, iftask, nmt, ierr
      real ttp
      complex f, g, bs, br, ss, sct
      dimension f(3,nxvh,kypd,kblok), g(3,nyv,kxp,jblok)
      dimension bs(3,kxp,kyp,kblok), br(3,kxp,kyp,jblok)
      dimension ss(ndim,nxvh,nmt+1)
      dimension mixup(nxhyd), sct(nxyhd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer i, nargs, margs, nx, ny, nxh, nmtt, kxpi, kxpp, kxpl
      integer kypi, kypp, kypl
      real tf
      double precision dtime
      external PFFT2RNXX, PFFT2RNXY
      data nargs, margs /15,16/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      nmtt = nmt + 1
      kxpp = kxp/nmtt
      kypp = kyp/nmtt
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXX,margs,f,ss(1,1,i),isign,m
     1ixup,sct,indx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,ndim,nxhyd,
     2nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT2RNXX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,kstrt,k
     1ypi,kypl,nxvh,kypd,kblok,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,j
     1blok,kblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,ndim,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT2RNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kx
     1p,kblok,jblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOSE(f,g,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kyp
     1d,jblok,kblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,ndim,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFFT2RNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOSE(g,f,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,k
     1blok,jblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXX,margs,f,ss(1,1,i),isign,m
     1ixup,sct,indx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,ndim,nxhyd,
     2nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call PFFT2RNXX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,kstrt,k
     1ypi,kypl,nxvh,kypd,kblok,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT2RXN(f,g,ss,isign,ntpose,mixup,sct,ttp,indx,indy,k
     1strt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim,nxhyd,nxyhd,kxyip,ifta
     2sk,nmt,ierr)
c multi-tasking real to complex fft
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, ndim, nxhyd, nxyhd, mixup
      integer kxyip, iftask, nmt, ierr
      real ttp
      complex f, g, ss, sct
      dimension f(3,nxvh,kypd,kblok), g(3,nyv,kxp,jblok)
      dimension ss(ndim,nxvh,nmt+1)
      dimension mixup(nxhyd), sct(nxyhd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer i, nargs, margs, nx, ny, nxh, nmtt, kxpi, kxpp, kxpl
      integer kypi, kypp, kypl
      real tf
      double precision dtime
      external PFFT2RNXX, PFFT2RNXY
      data nargs, margs /15,16/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      nmtt = nmt + 1
      kxpp = kxp/nmtt
      kypp = kyp/nmtt
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXX,margs,f,ss(1,1,i),isign,m
     1ixup,sct,indx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,ndim,nxhyd,
     2nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT2RNXX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,kstrt,k
     1ypi,kypl,nxvh,kypd,kblok,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jblok,
     1kblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,ndim,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT2RNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kbl
     1ok,jblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOSEX(f,g,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,kxp,kypd,jbl
     1ok,kblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXY,nargs,g,isign,mixup,sct,i
     1ndx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,ndim,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFFT2RNXY(g,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv,
     1kxp,jblok,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOSEX(g,f,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,kypd,kxp,kblok,
     1jblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXX,margs,f,ss(1,1,i),isign,m
     1ixup,sct,indx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,ndim,nxhyd,
     2nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call PFFT2RNXX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,kstrt,k
     1ypi,kypl,nxvh,kypd,kblok,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MP2FFT2RN(f1,f2,g1,g2,bs,br,ss,isign,ntpose,mixup,sct,t
     1tp,indx,indy,kstrt,nxvh,nyv,kxp,kyp,kypd,jblok,kblok,ndim1,ndim2,n
     2xhyd,nxyhd,kxyip,iftask,nmt,ierr)
c multi-tasking two real to complex ffts
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, kstrt, nxvh, nyv, kxp, kyp
      integer kypd, jblok, kblok, ndim1, ndim2, nxhyd, nxyhd, mixup
      integer kxyip, iftask, nmt, ierr
      real ttp
      complex f1, f2, g1, g2, bs, br, ss, sct
      dimension f1(ndim1,nxvh,kypd,kblok), f2(ndim2,nxvh,kypd,kblok)
      dimension g1(ndim1,nyv,kxp,jblok), g2(ndim2,nyv,kxp,jblok)
      dimension bs(ndim1+ndim2,kxp,kyp,kblok)
      dimension br(ndim1+ndim2,kxp,kyp,jblok)
      dimension ss(ndim1+ndim2,nxvh,nmt+1)
      dimension mixup(nxhyd), sct(nxyhd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer i, nargs, margs, nx, ny, nxh, nmtt, kxpi, kxpp, kxpl
      integer kypi, kypp, kypl
      real tf
      double precision dtime
      external PFFT2RNXX, PFFT2RNXY
      data nargs, margs /15,16/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      nmtt = nmt + 1
      kxpp = kxp/nmtt
      kypp = kyp/nmtt
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start first x fft tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXX,margs,f1,ss(1,1,i),isign,
     1mixup,sct,indx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,ndim1,nxhy
     2d,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish first x fft
         call PFFT2RNXX(f1,ss(1,1,nmtt),isign,mixup,sct,indx,indy,kstrt,
     1kypi,kypl,nxvh,kypd,kblok,ndim1,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start second x fft tasks
         do 30 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXX,margs,f2,ss(1,1,i),isign,
     1mixup,sct,indx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,ndim2,nxhy
     2d,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish second x fft
         call PFFT2RNXX(f2,ss(1,1,nmtt),isign,mixup,sct,indx,indy,kstrt,
     1kypi,kypl,nxvh,kypd,kblok,ndim2,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   40    continue
c transpose f arrays to g
         call PWTIMERA(-1,ttp,dtime)
         call PN2TPOSE(f1,f2,g1,g2,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,kyp,k
     1xp,kypd,jblok,kblok,ndim1,ndim2)
         call PWTIMERA(1,ttp,dtime)
c start first y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXY,nargs,g1,isign,mixup,sct,
     1indx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,ndim1,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish first y fft
         call PFFT2RNXY(g1,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv
     1,kxp,jblok,ndim1,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c start second y fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXY,nargs,g2,isign,mixup,sct,
     1indx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,ndim2,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish second y fft
         call PFFT2RNXY(g2,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv
     1,kxp,jblok,ndim2,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose g arrays to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PN2TPOSE(g1,g2,f1,f2,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kx
     1p,kypd,kxp,kblok,jblok,ndim1,ndim2)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f arrays to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PN2TPOSE(f1,f2,g1,g2,bs,br,nxh,ny,kstrt,nxvh,nyv,kxp,ky
     1p,kxp,kypd,jblok,kblok,ndim1,ndim2)
            call PWTIMERA(1,tf,dtime)
         endif
c start first y fft tasks
         do 90 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXY,nargs,g1,isign,mixup,sct,
     1indx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,ndim1,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   90    continue
c finish first y fft
         call PFFT2RNXY(g1,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv
     1,kxp,jblok,ndim1,nxhyd,nxyhd)
c wait for tasks to complete
         do 100 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  100    continue
c start second y fft tasks
         do 110 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXY,nargs,g2,isign,mixup,sct,
     1indx,indy,kstrt,kxyip(i),kxpp,nyv,kxp,jblok,ndim2,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish second y fft
         call PFFT2RNXY(g2,isign,mixup,sct,indx,indy,kstrt,kxpi,kxpl,nyv
     1,kxp,jblok,ndim2,nxhyd,nxyhd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  120    continue
c transpose g arrays to f
         call PWTIMERA(-1,ttp,dtime)
         call PN2TPOSE(g1,g2,f1,f2,br,bs,ny,nxh,kstrt,nyv,nxvh,kyp,kxp,k
     1ypd,kxp,kblok,jblok,ndim1,ndim2)
         call PWTIMERA(1,ttp,dtime)
c start first x fft tasks
         do 130 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXX,margs,f1,ss(1,1,i),isign,
     1mixup,sct,indx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,ndim1,nxhy
     2d,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish first x fft
         call PFFT2RNXX(f1,ss(1,1,nmtt),isign,mixup,sct,indx,indy,kstrt,
     1kypi,kypl,nxvh,kypd,kblok,ndim1,nxhyd,nxyhd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  140    continue
c start second x fft tasks
         do 150 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT2RNXX,margs,f2,ss(1,1,i),isign,
     1mixup,sct,indx,indy,kstrt,kxyip(i),kypp,nxvh,kypd,kblok,ndim2,nxhy
     2d,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  150    continue
c finish second x fft
         call PFFT2RNXX(f2,ss(1,1,nmtt),isign,mixup,sct,indx,indy,kstrt,
     1kypi,kypl,nxvh,kypd,kblok,ndim2,nxhyd,nxyhd)
c wait for tasks to complete
         do 160 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  160    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFSST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxyip,
     2iftask,nmt,ierr)
c multi-tasking real sine/sine transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFST2RXX, PFST2RXY
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y sine tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXY,nargs,g,isign,mixup,sctd,i
     1ndx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxpl
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y sine tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXY,nargs,g,isign,mixup,sctd,i
     1ndx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxpl
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x sine tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFSCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxyip,
     2iftask,nmt,ierr)
c multi-tasking real sine/cosine transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFST2RXX, PFCT2RXY
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y cosine tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXY,nargs,g,isign,mixup,sctd,i
     1ndx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxpl
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y cosine tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXY,nargs,g,isign,mixup,sctd,i
     1ndx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxpl
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x sine tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFCST2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxyip,
     2iftask,nmt,ierr)
c multi-tasking real cosine/sine transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFCT2RXX, PFST2RXY
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y sine tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXY,nargs,g,isign,mixup,sctd,i
     1ndx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxpl
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y sine tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXY,nargs,g,isign,mixup,sctd,i
     1ndx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine transform
         call PFST2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxpl
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x cosine tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFCCT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxyip,
     2iftask,nmt,ierr)
c multi-tasking real cosine/cosine transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2*nxvh,kypd,kblok), g(nyv,kxp2d,jblok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFCT2RXX, PFCT2RXY
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,ky
     1pd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y cosine tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXY,nargs,g,isign,mixup,sctd,i
     1ndx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxpl
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d
     1,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y cosine tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXY,nargs,g,isign,mixup,sctd,i
     1ndx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y cosine transform
         call PFCT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxpl
     1,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kxp
     12d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x cosine tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x scoine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFCST2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxyip
     2,iftask,nmt,ierr)
c multi-tasking 2 real sine-cosine transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok), g(2,nyv,kxp2d,jblok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFCST2R2X, PFSCT2R2Y
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCST2R2X,nargs,f,isign,mixup,sctd,
     1indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transform
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y sine-cosine tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine-cosine transform
         call PFSCT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y sine-cosine tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine-cosine transforms
         call PFSCT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x sine-cosine tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCST2R2X,nargs,f,isign,mixup,sctd,
     1indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFSCT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxyip
     2,iftask,nmt,ierr)
c multi-tasking 2 real sine-cosine transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real f, g, bs, br, ttp
      complex sctd
      dimension f(2,2*nxvh,kypd,kblok), g(2,nyv,kxp2d,jblok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFSCT2R2X, PFCST2R2Y
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCT2R2X,nargs,f,isign,mixup,sctd,
     1indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transform
         call PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y sine-cosine tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCST2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine-cosine transform
         call PFCST2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y sine-cosine tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCST2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine-cosine transforms
         call PFCST2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x sine-cosine tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCT2R2X,nargs,f,isign,mixup,sctd,
     1indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFCST2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxyip
     2,iftask,nmt,ierr)
c multi-tasking 3 real sine-cosine transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok), g(3,nyv,kxp2d,jblok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFCSST2R3X, PFSCST2R3Y
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCSST2R3X,nargs,f,isign,mixup,sctd
     1,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transform
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1l,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y sine-cosine tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCST2R3Y,nargs,g,isign,mixup,sctd
     1,indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine-cosine transform
         call PFSCST2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1pl,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y sine-cosine tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCST2R3Y,nargs,g,isign,mixup,sctd
     1,indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine-cosine transforms
         call PFSCST2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1pl,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x sine-cosine tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCSST2R3X,nargs,f,isign,mixup,sctd
     1,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1l,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFSCT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,in
     1dy,kstrt,nxvh,nyv,kxp2,kyp,kypd,kxp2d,jblok,kblok,nxhyd,nxyd,kxyip
     2,iftask,nmt,ierr)
c multi-tasking 3 real sine-cosine transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyv
      integer kxp2, kyp, kypd, kxp2d, jblok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real f, g, bs, br, ttp
      complex sctd
      dimension f(3,2*nxvh,kypd,kblok), g(3,nyv,kxp2d,jblok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,jblok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i
      real tf
      double precision dtime
      external PFSCCT2R3X, PFCSCT2R3Y
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCCT2R3X,nargs,f,isign,mixup,sctd
     1,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transform
         call PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1l,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2d,k
     1ypd,jblok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y sine-cosine tasks
         do 30 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCSCT2R3Y,nargs,g,isign,mixup,sctd
     1,indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine-cosine transform
         call PFCSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1pl,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,jblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,nyv,kxp2,kyp,kxp2
     1d,kypd,jblok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y sine-cosine tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCSCT2R3Y,nargs,g,isign,mixup,sctd
     1,indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyv,kxp2d,jblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine-cosine transforms
         call PFCSCT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kx
     1pl,nyv,kxp2d,jblok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,nyv,2*nxvh,kyp,kxp2,kypd,kx
     1p2d,kblok,jblok)
         call PWTIMERA(1,ttp,dtime)
c start x sine-cosine tasks
         do 70 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCCT2R3X,nargs,f,isign,mixup,sctd
     1,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1l,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   80    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFSFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,kxyi
     2p,iftask,nmt,ierr)
c multi-tasking real sine/periodic transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh1, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i, j, l
      real tf
      double precision dtime
      external PFST2RXX, PFDFT2RXY
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c start x sine tasks
         do 30 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   40    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,ky
     1pd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kx
     1p2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 100 l = 1, j2blok
         do 90 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   90    continue
  100    continue
c start x sine tasks
         do 110 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFST2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x sine transform
         call PFST2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFCFT2R(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,ind
     1y,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,kxyi
     2p,iftask,nmt,ierr)
c multi-tasking real cosine/periodic transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh1, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i, j, l
      real tf
      double precision dtime
      external PFCT2RXX, PFDFT2RXY
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c start x cosine tasks
         do 30 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   40    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,ky
     1pd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kx
     1p2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 100 l = 1, j2blok
         do 90 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   90    continue
  100    continue
c start x cosine tasks
         do 110 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCT2RXX,nargs,f,isign,mixup,sctd,i
     1ndx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x cosine transform
         call PFCT2RXX(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl,
     1nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFCSFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,i
     1ndy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,kx
     2yip,iftask,nmt,ierr)
c multi-tasking 2 real cosine-sine/periodic transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd
      dimension f(2,nxvh,kypd,kblok), g(2,nyvh,kxp2d,j2blok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh1, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i, j, l
      real tf
      double precision dtime
      external PFCST2R2X, PFDFT2R2Y
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c start x cosine-sine tasks
         do 40 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCST2R2X,nargs,f,isign,mixup,sctd,
     1indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   40    continue
c finish x cosine-sine transform
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 50 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   50    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 60 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   60    continue
c finish y fft
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 70 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   70    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 80 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   80    continue
c finish y fft
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 90 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   90    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 120 l = 1, j2blok
         do 110 j = 1, nxh1
         do 100 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
  100    continue
  110    continue
  120    continue
c start x cosine-sine tasks
         do 130 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCST2R2X,nargs,f,isign,mixup,sctd,
     1indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish x cosine-sine transforms
         call PFCST2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  140    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFSCFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,i
     1ndy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,kx
     2yip,iftask,nmt,ierr)
c multi-tasking 2 real sine-cosine/periodic transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd
      dimension f(2,nxvh,kypd,kblok), g(2,nyvh,kxp2d,j2blok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh1, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i, j, l
      real tf
      double precision dtime
      external PFSCT2R2X, PFDFT2R2Y
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c start x sine-cosine tasks
         do 40 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCT2R2X,nargs,f,isign,mixup,sctd,
     1indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   40    continue
c finish x sine-cosine transform
         call PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 50 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   50    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 60 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   60    continue
c finish y fft
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 70 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   70    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 80 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   80    continue
c finish y fft
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 90 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   90    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 120 l = 1, j2blok
         do 110 j = 1, nxh1
         do 100 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
  100    continue
  110    continue
  120    continue
c start x sine-cosine tasks
         do 130 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCT2R2X,nargs,f,isign,mixup,sctd,
     1indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish x sine-cosine transforms
         call PFSCT2R2X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kypl
     1,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  140    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFCSFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,i
     1ndy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,kx
     2yip,iftask,nmt,ierr)
c multi-tasking 3 real cosine-sine/periodic transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd
      dimension f(3,nxvh,kypd,kblok), g(3,nyvh,kxp2d,j2blok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh1, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i, j, l
      real tf
      double precision dtime
      external PFCSST2R3X, PFDFT2R3Y
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c start x cosine-sine tasks
         do 40 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCSST2R3X,nargs,f,isign,mixup,sctd
     1,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   40    continue
c finish x cosine-sine transform
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1l,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 50 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   50    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 60 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R3Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   60    continue
c finish y fft
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 70 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   70    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 80 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R3Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   80    continue
c finish y fft
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 90 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   90    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 120 l = 1, j2blok
         do 110 j = 1, nxh1
         do 100 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
  100    continue
  110    continue
  120    continue
c start x cosine-sine tasks
         do 130 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFCSST2R3X,nargs,f,isign,mixup,sctd
     1,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish x cosine-sine transforms
         call PFCSST2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1l,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  140    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFSCFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,ttp,indx,i
     1ndy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,nxyd,kx
     2yip,iftask,nmt,ierr)
c multi-tasking 3 real sine-cosine/periodic transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd
      dimension f(3,nxvh,kypd,kblok), g(3,nyvh,kxp2d,j2blok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nxh1, kxpi, kxpp, kxpl, kypi, kypp, kypl
      integer i, j, l
      real tf
      double precision dtime
      external PFSCCT2R3X, PFDFT2R3Y
      data nargs /15/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c start x sine-scoine tasks
         do 40 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCCT2R3X,nargs,f,isign,mixup,sctd
     1,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   40    continue
c finish x sine-cosine transform
         call PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1l,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 50 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   50    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 60 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R3Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   60    continue
c finish y fft
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 70 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   70    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 80 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R3Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   80    continue
c finish y fft
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 90 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   90    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 120 l = 1, j2blok
         do 110 j = 1, nxh1
         do 100 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
  100    continue
  110    continue
  120    continue
c start x sine-cosine tasks
         do 130 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFSCCT2R3X,nargs,f,isign,mixup,sctd
     1,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish x sine-cosine transforms
         call PFSCCT2R3X(f,isign,mixup,sctd,indx,indy,kstrt,kyp,kypi,kyp
     1l,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  140    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFDSFT2RX(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp,
     1indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,n
     2xyd,kxyip,iftask,nmt,ierr)
c multi-tasking real sine/periodic transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nx, ny, nxh1, kxpi, kxpp, kxpl
      integer kypi, kypp, kypl, i, j, l
      real tf
      double precision dtime
      external PFDST2RXX, PFDFT2RXY
      data nargs, margs /15,16/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c start x sine tasks
         do 30 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDST2RXX,margs,f,isign,mixup,sctd,
     1sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd
     2)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c perform x sine transform
         call PFDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,kyp
     1i,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   40    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,ky
     1pd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kx
     1p2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 100 l = 1, j2blok
         do 90 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   90    continue
  100    continue
c start x sine tasks
         do 110 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDST2RXX,margs,f,isign,mixup,sctd,
     1sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd
     2)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x sine transform
         call PFDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,kyp
     1i,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFDCFT2RX(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp,
     1indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,n
     2xyd,kxyip,iftask,nmt,ierr)
c multi-tasking real cosine/periodic transform
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(nxvh,kypd,kblok), g(nyvh,kxp2d,j2blok)
      dimension bs(kxp2+1,kyp+1,kblok), br(kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nx, ny, nxh1, kxpi, kxpp, kxpl
      integer kypi, kypp, kypl, i, j, l
      real tf
      double precision dtime
      external PFDCT2RXX, PFDFT2RXY
      data nargs, margs /15,16/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 20 l = 1, j2blok
         do 10 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
c start x cosine tasks
         do 30 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDCT2RXX,margs,f,isign,mixup,sctd,
     1sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd
     2)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish x cosine transform
         call PFDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,kyp
     1i,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   40    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2d
     1,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 50 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,ky
     1pd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PRTPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kx
     1p2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 70 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2RXY,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish y fft
         call PFDFT2RXY(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PRTPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd,
     1kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 100 l = 1, j2blok
         do 90 j = 1, nxh1
         f(j,kyp+1,l) = cmplx(0.,0.)
   90    continue
  100    continue
c start x cosine tasks
         do 110 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDCT2RXX,margs,f,isign,mixup,sctd,
     1sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxyd
     2)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x cosine transform
         call PFDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,kyp
     1i,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFDCSFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp
     1,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,
     2nxyd,kxyip,iftask,nmt,ierr)
c multi-tasking 2 real cosine-sine/periodic transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(2,nxvh,kypd,kblok), g(2,nyvh,kxp2d,j2blok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nx, ny, nxh1, kxpi, kxpp, kxpl
      integer kypi, kypp, kypl, i, j, l
      real tf
      double precision dtime
      external PFDCST2R2X, PFDFT2R2Y
      data nargs, margs /15,16/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c start x cosine-sine tasks
         do 40 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDCST2R2X,margs,f,isign,mixup,sctd
     1,sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxy
     2d)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   40    continue
c finish x cosine-sine transform
         call PFDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,ky
     1pi,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 50 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   50    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 60 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   60    continue
c finish y fft
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 70 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   70    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 80 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   80    continue
c finish y fft
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 90 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   90    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 120 l = 1, j2blok
         do 110 j = 1, nxh1
         do 100 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
  100    continue
  110    continue
  120    continue
c start x cosine-sine tasks
         do 130 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDCST2R2X,margs,f,isign,mixup,sctd
     1,sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxy
     2d)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish x cosine-sine transforms
         call PFDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,ky
     1pi,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  140    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFDSCFT2R2(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp
     1,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,
     2nxyd,kxyip,iftask,nmt,ierr)
c multi-tasking 2 real sine-cosine/periodic transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(2,nxvh,kypd,kblok), g(2,nyvh,kxp2d,j2blok)
      dimension bs(2,kxp2+1,kyp+1,kblok), br(2,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nx, ny, nxh1, kxpi, kxpp, kxpl
      integer kypi, kypp, kypl, i, j, l
      real tf
      double precision dtime
      external PFDSCT2R2X, PFDFT2R2Y
      data nargs, margs /15,16/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c start x sine-cosine tasks
         do 40 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDSCT2R2X,margs,f,isign,mixup,sctd
     1,sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxy
     2d)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   40    continue
c finish x sine-cosine transform
         call PFDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,ky
     1pi,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 50 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   50    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 60 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   60    continue
c finish y fft
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 70 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   70    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR2TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 80 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R2Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   80    continue
c finish y fft
         call PFDFT2R2Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 90 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   90    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR2TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 120 l = 1, j2blok
         do 110 j = 1, nxh1
         do 100 i = 1, 2
         f(i,j,kyp+1,l) = cmplx(0.,0.)
  100    continue
  110    continue
  120    continue
c start x sine-cosine tasks
         do 130 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDSCT2R2X,margs,f,isign,mixup,sctd
     1,sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nxy
     2d)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish x sine-cosine transforms
         call PFDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,ky
     1pi,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  140    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFDCSFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp
     1,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,
     2nxyd,kxyip,iftask,nmt,ierr)
c multi-tasking 3 real cosine-sine/periodic transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(3,nxvh,kypd,kblok), g(3,nyvh,kxp2d,j2blok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nx, ny, nxh1, kxpi, kxpp, kxpl
      integer kypi, kypp, kypl, i, j, l
      real tf
      double precision dtime
      external PFDCSST2R3X, PFDFT2R3Y
      data nargs, margs /15,16/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c start x cosine-sine tasks
         do 40 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDCSST2R3X,margs,f,isign,mixup,sct
     1d,sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nx
     2yd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   40    continue
c finish x cosine-sine transform
         call PFDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,k
     1ypi,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 50 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   50    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 60 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R3Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   60    continue
c finish y fft
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 70 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   70    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 80 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R3Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   80    continue
c finish y fft
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 90 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   90    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 120 l = 1, j2blok
         do 110 j = 1, nxh1
         do 100 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
  100    continue
  110    continue
  120    continue
c start x cosine-sine tasks
         do 130 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDCSST2R3X,margs,f,isign,mixup,sct
     1d,sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nx
     2yd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish x cosine-sine transforms
         call PFDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,k
     1ypi,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  140    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFDSCFT2R3(f,g,bs,br,isign,ntpose,mixup,sctd,sctdx,ttp
     1,indx,indy,kstrt,nxvh,nyvh,kxp2,kyp,kypd,kxp2d,j2blok,kblok,nxhyd,
     2nxyd,kxyip,iftask,nmt,ierr)
c multi-tasking 3 real sine-cosine/periodic transforms
c kxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, mixup, indx, indy, kstrt, nxvh, nyvh
      integer kxp2, kyp, kypd, kxp2d, j2blok, kblok, nxhyd, nxyd
      integer kxyip, iftask, nmt, ierr
      real bs, br, ttp
      complex f, g, sctd, sctdx
      dimension f(3,nxvh,kypd,kblok), g(3,nyvh,kxp2d,j2blok)
      dimension bs(3,kxp2+1,kyp+1,kblok), br(3,kxp2+1,kyp+1,j2blok)
      dimension mixup(nxhyd), sctd(nxyd), sctdx(2*nxvh)
      dimension kxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nx, ny, nxh1, kxpi, kxpp, kxpl
      integer kypi, kypp, kypl, i, j, l
      real tf
      double precision dtime
      external PFDSCCT2R3X, PFDFT2R3Y
      data nargs, margs /15,16/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh1 = nx/2 + 1
      kxpp = kxp2/(nmt + 1)
      kypp = kyp/(nmt + 1)
      kxpi = kxpp*nmt
      kypi = kypp*nmt
      kxpl = kxp2 - kxpi
      kypl = kyp - kypi
      kxpi = kxpi + 1
      kypi = kypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c zero out guard cells
         do 30 l = 1, j2blok
         do 20 j = 1, nxh1
         do 10 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
   10    continue
   20    continue
   30    continue
c start x sine-scoine tasks
         do 40 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDSCCT2R3X,margs,f,isign,mixup,sct
     1d,sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nx
     2yd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   40    continue
c finish x sine-cosine transform
         call PFDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,k
     1ypi,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 50 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   50    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,kxp2
     1d,kypd,j2blok,kblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 60 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R3Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   60    continue
c finish y fft
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 70 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   70    continue
c transpose g array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,k
     1ypd,kxp2d,kblok,j2blok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to g
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PR3TPOSE(f,g,bs,br,nx,ny,kstrt,2*nxvh,2*nyvh,kxp2,kyp,k
     1xp2d,kypd,j2blok,kblok)
            call PWTIMERA(1,tf,dtime)
         endif
c start y fft tasks
         do 80 i = 1, nmt
         kxyip(i) = kxpp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDFT2R3Y,nargs,g,isign,mixup,sctd,
     1indx,indy,kstrt,kxp2,kxyip(i),kxpp,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   80    continue
c finish y fft
         call PFDFT2R3Y(g,isign,mixup,sctd,indx,indy,kstrt,kxp2,kxpi,kxp
     1l,nyvh,kxp2d,j2blok,nxhyd,nxyd)
c wait for tasks to complete
         do 90 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   90    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PR3TPOSE(g,f,br,bs,ny,nx,kstrt,2*nyvh,2*nxvh,kyp,kxp2,kypd
     1,kxp2d,kblok,j2blok)
         call PWTIMERA(1,ttp,dtime)
c zero out guard cells
         do 120 l = 1, j2blok
         do 110 j = 1, nxh1
         do 100 i = 1, 3
         f(i,j,kyp+1,l) = cmplx(0.,0.)
  100    continue
  110    continue
  120    continue
c start x sine-cosine tasks
         do 130 i = 1, nmt
         kxyip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFDSCCT2R3X,margs,f,isign,mixup,sct
     1d,sctdx,indx,indy,kstrt,kyp,kxyip(i),kypp,nxvh,kypd,kblok,nxhyd,nx
     2yd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish x sine-cosine transforms
         call PFDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,kstrt,kyp,k
     1ypi,kypl,nxvh,kypd,kblok,nxhyd,nxyd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  140    continue
      endif
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
