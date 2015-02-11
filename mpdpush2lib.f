c-----------------------------------------------------------------------
c 2d PIC parallel multi-tasking library for pushing particles
c with darwin electric magnetic fields and depositing current
c and derivative of current
c mpdpush2lib.f contains multi-tasking procedures to process particles
c              with darwin electric and magnetic fields:
c MPGMJPOST2 multi-tasking wrapper for PGMJPOST2
c MPGSMJPOST2 multi-tasking wrapper for PGSMJPOST2
c MPGMJPOST2L multi-tasking wrapper for PGMJPOST2L
c MPGSMJPOST2L multi-tasking wrapper for PGSMJPOST2L
c MPGMJPOST22 multi-tasking wrapper for PGMJPOST22
c MPGSMJPOST22 multi-tasking wrapper for PGSMJPOST22
c MPGMJPOST22L multi-tasking wrapper for PGMJPOST22L
c MPGSMJPOST22L multi-tasking wrapper for PGSMJPOST22L
c MPGDCJPOST2 multi-tasking wrapper for PGDCJPOST2
c MPGSDCJPOST2 multi-tasking wrapper for PGSDCJPOST2
c MPGDCJPOST2L multi-tasking wrapper for PGDCJPOST2L
c MPGSDCJPOST2L multi-tasking wrapper for PGSDCJPOST2L
c MPGDCJPOST22 multi-tasking wrapper for PGDCJPOST22
c MPGSDCJPOST22 multi-tasking wrapper for PGSDCJPOST22
c MPGDCJPOST22L multi-tasking wrapper for PGDCJPOST22L
c MPGSDCJPOST22L multi-tasking wrapper for PGSDCJPOST22L
c written by viktor k. decyk, ucla
c copyright 2006, regents of the university of california
c update: december 11, 2009
c-----------------------------------------------------------------------
      subroutine MPGMJPOST2(part,amu,npp,nps,noff,qm,idimp,npmax,nblok,n
     1xv,nypmx,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(4,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGMJPOST2
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear momentum flux arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 4
      amup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGMJPOST2,nargs,part(1,npo,1),amup(1,1
     1,1,1,i),nps,noff,qm,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining momentum flux
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGMJPOST2(part(1,npo,1),amu,npp,noff,qm,idimp,npmax,nblok,nxv
     1,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 4
      amu(m,j,k,l) = amu(m,j,k,l) + amup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSMJPOST2(part,amu,npp,nps,noff,qm,idimp,npmax,nblok,
     1nxv,nxyp,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(4,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSMJPOST2
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear momentum flux arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 4
      amup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSMJPOST2,nargs,part(1,npo,1),amup(1,
     11,1,i),nps,noff,qm,idimp,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining momentum flux
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSMJPOST2(part(1,npo,1),amu,npp,noff,qm,idimp,npmax,nblok,nx
     1v,nxyp)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 4
      amu(m,j,l) = amu(m,j,l) + amup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGMJPOST2L(part,amu,npp,nps,noff,qm,idimp,npmax,nblok,
     1nxv,nypmx,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(4,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGMJPOST2L
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear momentum flux arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 4
      amup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGMJPOST2L,nargs,part(1,npo,1),amup(1,
     11,1,1,i),nps,noff,qm,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining momentum flux
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGMJPOST2L(part(1,npo,1),amu,npp,noff,qm,idimp,npmax,nblok,nx
     1v,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 4
      amu(m,j,k,l) = amu(m,j,k,l) + amup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSMJPOST2L(part,amu,npp,nps,noff,qm,idimp,npmax,nblok
     1,nxv,nxyp,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(4,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSMJPOST2L
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear momentum flux arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 4
      amup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSMJPOST2L,nargs,part(1,npo,1),amup(1
     1,1,1,i),nps,noff,qm,idimp,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining momentum flux
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSMJPOST2L(part(1,npo,1),amu,npp,noff,qm,idimp,npmax,nblok,n
     1xv,nxyp)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 4
      amu(m,j,l) = amu(m,j,l) + amup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGMJPOST22(part,amu,npp,nps,noff,qm,idimp,npmax,nblok,
     1nxv,nypmx,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGMJPOST22
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear momentum flux arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 2
      amup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGMJPOST22,nargs,part(1,npo,1),amup(1,
     11,1,1,i),nps,noff,qm,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining momentum flux
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGMJPOST22(part(1,npo,1),amu,npp,noff,qm,idimp,npmax,nblok,nx
     1v,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 2
      amu(m,j,k,l) = amu(m,j,k,l) + amup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSMJPOST22(part,amu,npp,nps,noff,qm,idimp,npmax,nblok
     1,nxv,nxyp,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSMJPOST22
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear momentum flux arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 2
      amup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSMJPOST22,nargs,part(1,npo,1),amup(1
     1,1,1,i),nps,noff,qm,idimp,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining momentum flux
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSMJPOST22(part(1,npo,1),amu,npp,noff,qm,idimp,npmax,nblok,n
     1xv,nxyp)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 2
      amu(m,j,l) = amu(m,j,l) + amup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGMJPOST22L(part,amu,npp,nps,noff,qm,idimp,npmax,nblok
     1,nxv,nypmx,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGMJPOST22L
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear momentum flux arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 2
      amup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGMJPOST22L,nargs,part(1,npo,1),amup(1
     1,1,1,1,i),nps,noff,qm,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining momentum flux
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGMJPOST22L(part(1,npo,1),amu,npp,noff,qm,idimp,npmax,nblok,n
     1xv,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 2
      amu(m,j,k,l) = amu(m,j,k,l) + amup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSMJPOST22L(part,amu,npp,nps,noff,qm,idimp,npmax,nblo
     1k,nxv,nxyp,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSMJPOST22L
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear momentum flux arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 2
      amup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSMJPOST22L,nargs,part(1,npo,1),amup(
     11,1,1,i),nps,noff,qm,idimp,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining momentum flux
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSMJPOST22L(part(1,npo,1),amu,npp,noff,qm,idimp,npmax,nblok,
     1nxv,nxyp)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 2
      amu(m,j,l) = amu(m,j,l) + amup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm,qbm
     1,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension cu(3,nxv,nypmx,nblok), dcu(3,nxv,nypmx,nblok)
      dimension amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxv,nypmx,nblok,nmt), dcup(3,nxv,nypmx,nblok,nmt)
      dimension amup(4,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGDCJPOST2
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 3
      cup(m,j,k,l,i) = 0.
      dcup(m,j,k,l,i) = 0.
      amup(m,j,k,l,i) = 0.
   20 continue
      amup(4,j,k,l,i) = 0.
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGDCJPOST2,nargs,part(1,npo,1),fxy,bxy
     1,nps,noff,cup(1,1,1,1,i),dcup(1,1,1,1,i),amup(1,1,1,1,i),qm,qbm,dt
     2,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGDCJPOST2(part(1,npo,1),fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,d
     1t,idimp,npmax,nblok,nxv,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 3
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
      dcu(m,j,k,l) = dcu(m,j,k,l) + dcup(m,j,k,l,i)
      amu(m,j,k,l) = amu(m,j,k,l) + amup(m,j,k,l,i)
   90 continue
      amu(4,j,k,l) = amu(4,j,k,l) + amup(4,j,k,l,i)
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm,qb
     1m,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension cu(3,nxyp,nblok), dcu(3,nxyp,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxyp,nblok,nmt), dcup(3,nxyp,nblok,nmt)
      dimension amup(4,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSDCJPOST2
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 3
      cup(m,j,l,i) = 0.
      dcup(m,j,l,i) = 0.
      amup(m,j,l,i) = 0.
   20 continue
      amup(4,j,l,i) = 0.
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSDCJPOST2,nargs,part(1,npo,1),fxy,bx
     1y,nps,noff,cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idim
     2p,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSDCJPOST2(part(1,npo,1),fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,
     1dt,idimp,npmax,nblok,nxv,nxyp)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 3
      cu(m,j,l) = cu(m,j,l) + cup(m,j,l,i)
      dcu(m,j,l) = dcu(m,j,l) + dcup(m,j,l,i)
      amu(m,j,l) = amu(m,j,l) + amup(m,j,l,i)
   80 continue
      amu(4,j,l) = amu(4,j,l) + amup(4,j,l,i)
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm,qb
     1m,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxv,nypmx,nblok), bxy(3,nxv,nypmx,nblok)
      dimension cu(3,nxv,nypmx,nblok), dcu(3,nxv,nypmx,nblok)
      dimension amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxv,nypmx,nblok,nmt), dcup(3,nxv,nypmx,nblok,nmt)
      dimension amup(4,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGDCJPOST2L
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 3
      cup(m,j,k,l,i) = 0.
      dcup(m,j,k,l,i) = 0.
      amup(m,j,k,l,i) = 0.
   20 continue
      amup(4,j,k,l,i) = 0.
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGDCJPOST2L,nargs,part(1,npo,1),fxy,bx
     1y,nps,noff,cup(1,1,1,1,i),dcup(1,1,1,1,i),amup(1,1,1,1,i),qm,qbm,d
     2t,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGDCJPOST2L(part(1,npo,1),fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,
     1dt,idimp,npmax,nblok,nxv,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 3
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
      dcu(m,j,k,l) = dcu(m,j,k,l) + dcup(m,j,k,l,i)
      amu(m,j,k,l) = amu(m,j,k,l) + amup(m,j,k,l,i)
   90 continue
      amu(4,j,k,l) = amu(4,j,k,l) + amup(4,j,k,l,i)
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm,q
     1bm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(3,nxyp,nblok), bxy(3,nxyp,nblok)
      dimension cu(3,nxyp,nblok), dcu(3,nxyp,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(3,nxyp,nblok,nmt), dcup(3,nxyp,nblok,nmt)
      dimension amup(4,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSDCJPOST2L
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 3
      cup(m,j,l,i) = 0.
      dcup(m,j,l,i) = 0.
      amup(m,j,l,i) = 0.
   20 continue
      amup(4,j,l,i) = 0.
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSDCJPOST2L,nargs,part(1,npo,1),fxy,b
     1xy,nps,noff,cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idi
     2mp,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSDCJPOST2L(part(1,npo,1),fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm
     1,dt,idimp,npmax,nblok,nxv,nxyp)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 3
      cu(m,j,l) = cu(m,j,l) + cup(m,j,l,i)
      dcu(m,j,l) = dcu(m,j,l) + dcup(m,j,l,i)
      amu(m,j,l) = amu(m,j,l) + amup(m,j,l,i)
   80 continue
      amu(4,j,l) = amu(4,j,l) + amup(4,j,l,i)
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGDCJPOST22(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm,qbm
     1,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension cu(2,nxv,nypmx,nblok), dcu(2,nxv,nypmx,nblok)
      dimension amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), dcup(2,nxv,nypmx,nblok,nmt)
      dimension amup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGDCJPOST22
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 2
      cup(m,j,k,l,i) = 0.
      dcup(m,j,k,l,i) = 0.
      amup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGDCJPOST22,nargs,part(1,npo,1),fxy,bz
     1,nps,noff,cup(1,1,1,1,i),dcup(1,1,1,1,i),amup(1,1,1,1,i),qm,qbm,dt
     2,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGDCJPOST22(part(1,npo,1),fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,d
     1t,idimp,npmax,nblok,nxv,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 2
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
      dcu(m,j,k,l) = dcu(m,j,k,l) + dcup(m,j,k,l,i)
      amu(m,j,k,l) = amu(m,j,k,l) + amup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSDCJPOST22(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm,qb
     1m,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension cu(2,nxyp,nblok), dcu(2,nxyp,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxyp,nblok,nmt), dcup(2,nxyp,nblok,nmt)
      dimension amup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSDCJPOST22
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 2
      cup(m,j,l,i) = 0.
      dcup(m,j,l,i) = 0.
      amup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSDCJPOST22,nargs,part(1,npo,1),fxy,b
     1z,nps,noff,cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idim
     2p,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSDCJPOST22(part(1,npo,1),fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,
     1dt,idimp,npmax,nblok,nxv,nxyp)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 2
      cu(m,j,l) = cu(m,j,l) + cup(m,j,l,i)
      dcu(m,j,l) = dcu(m,j,l) + dcup(m,j,l,i)
      amu(m,j,l) = amu(m,j,l) + amup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGDCJPOST22L(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm,qb
     1m,dt,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), bz(nxv,nypmx,nblok)
      dimension cu(2,nxv,nypmx,nblok), dcu(2,nxv,nypmx,nblok)
      dimension amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), dcup(2,nxv,nypmx,nblok,nmt)
      dimension amup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGDCJPOST22L
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 50 l = 1, nblok
      nps(l) = npt
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 m = 1, 2
      cup(m,j,k,l,i) = 0.
      dcup(m,j,k,l,i) = 0.
      amup(m,j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGDCJPOST22L,nargs,part(1,npo,1),fxy,b
     1z,nps,noff,cup(1,1,1,1,i),dcup(1,1,1,1,i),amup(1,1,1,1,i),qm,qbm,d
     2t,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining current
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      call PGDCJPOST22L(part(1,npo,1),fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,
     1dt,idimp,npmax,nblok,nxv,nypmx)
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 120 l = 1, nblok
      do 110 k = 1, nypmx
      do 100 j = 1, nxv
      do 90 m = 1, 2
      cu(m,j,k,l) = cu(m,j,k,l) + cup(m,j,k,l,i)
      dcu(m,j,k,l) = dcu(m,j,k,l) + dcup(m,j,k,l,i)
      amu(m,j,k,l) = amu(m,j,k,l) + amup(m,j,k,l,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSDCJPOST22L(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm,q
     1bm,dt,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok), bz(nxyp,nblok)
      dimension cu(2,nxyp,nblok), dcu(2,nxyp,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxyp,nblok,nmt), dcup(2,nxyp,nblok,nmt)
      dimension amup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSDCJPOST22L
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 j = 1, nxyp
      do 20 m = 1, 2
      cup(m,j,l,i) = 0.
      dcup(m,j,l,i) = 0.
      amup(m,j,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSDCJPOST22L,nargs,part(1,npo,1),fxy,
     1bz,nps,noff,cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idi
     2mp,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGSDCJPOST22L(part(1,npo,1),fxy,bz,npp,noff,cu,dcu,amu,qm,qbm
     1,dt,idimp,npmax,nblok,nxv,nxyp)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 l = 1, nblok
      do 90 j = 1, nxyp
      do 80 m = 1, 2
      cu(m,j,l) = cu(m,j,l) + cup(m,j,l,i)
      dcu(m,j,l) = dcu(m,j,l) + dcup(m,j,l,i)
      amu(m,j,l) = amu(m,j,l) + amup(m,j,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
