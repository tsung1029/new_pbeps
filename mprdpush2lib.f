c-----------------------------------------------------------------------
c 2d PIC parallel multi-tasking library for pushing relativistic
c particles with darwin electric magnetic fields and depositing current
c and derivative of current
c mprdpush2lib.f contains multi-tasking procedures to process
c                relativistic particles with darwin electric and
c                magnetic fields:
c MPGRMJPOST2 multi-tasking wrapper for PGRMJPOST2
c MPGSRMJPOST2 multi-tasking wrapper for PGSRMJPOST2
c MPGRMJPOST2L multi-tasking wrapper for PGRMJPOST2L
c MPGSRMJPOST2L multi-tasking wrapper for PGSRMJPOST2L
c MPGRMJPOST22 multi-tasking wrapper for PGRMJPOST22
c MPGSRMJPOST22 multi-tasking wrapper for PGSRMJPOST22
c MPGRMJPOST22L multi-tasking wrapper for PGRMJPOST22L
c MPGSRMJPOST22L multi-tasking wrapper for PGSRMJPOST22L
c MPGRDCJPOST2 multi-tasking wrapper for PGRDCJPOST2
c MPGSRDCJPOST2 multi-tasking wrapper for PGSRDCJPOST2
c MPGRDCJPOST2L multi-tasking wrapper for PGRDCJPOST2L
c MPGSRDCJPOST2L multi-tasking wrapper for PGSRDCJPOST2L
c MPGRDCJPOST22 multi-tasking wrapper for PGRDCJPOST22
c MPGSRDCJPOST22 multi-tasking wrapper for PGSRDCJPOST22
c MPGRDCJPOST22L multi-tasking wrapper for PGRDCJPOST22L
c MPGSRDCJPOST22L multi-tasking wrapper for PGSRDCJPOST22L
c written by viktor k. decyk, ucla
c copyright 2006, regents of the university of california
c update: december 11, 2009
c-----------------------------------------------------------------------
      subroutine MPGRMJPOST2(part,amu,npp,nps,noff,qm,ci,idimp,npmax,nbl
     1ok,nxv,nypmx,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(4,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRMJPOST2
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),PGRMJPOST2,nargs,part(1,npo,1),amup(1,
     11,1,1,i),nps,noff,qm,ci,idimp,npmax,nblok,nxv,nypmx)
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
      call PGRMJPOST2(part(1,npo,1),amu,npp,noff,qm,ci,idimp,npmax,nblok
     1,nxv,nypmx)
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
      subroutine MPGSRMJPOST2(part,amu,npp,nps,noff,qm,ci,idimp,npmax,nb
     1lok,nxv,nxyp,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(4,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSRMJPOST2
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),PGSRMJPOST2,nargs,part(1,npo,1),amup(1
     1,1,1,i),nps,noff,qm,ci,idimp,npmax,nblok,nxv,nxyp)
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
      call PGSRMJPOST2(part(1,npo,1),amu,npp,noff,qm,ci,idimp,npmax,nblo
     1k,nxv,nxyp)
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
      subroutine MPGRMJPOST2L(part,amu,npp,nps,noff,qm,ci,idimp,npmax,nb
     1lok,nxv,nypmx,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(4,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(4,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRMJPOST2L
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),PGRMJPOST2L,nargs,part(1,npo,1),amup(1
     1,1,1,1,i),nps,noff,qm,ci,idimp,npmax,nblok,nxv,nypmx)
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
      call PGRMJPOST2L(part(1,npo,1),amu,npp,noff,qm,ci,idimp,npmax,nblo
     1k,nxv,nypmx)
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
      subroutine MPGSRMJPOST2L(part,amu,npp,nps,noff,qm,ci,idimp,npmax,n
     1blok,nxv,nxyp,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(4,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(4,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSRMJPOST2L
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),PGSRMJPOST2L,nargs,part(1,npo,1),amup(
     11,1,1,i),nps,noff,qm,ci,idimp,npmax,nblok,nxv,nxyp)
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
      call PGSRMJPOST2L(part(1,npo,1),amu,npp,noff,qm,ci,idimp,npmax,nbl
     1ok,nxv,nxyp)
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
      subroutine MPGRMJPOST22(part,amu,npp,nps,noff,qm,ci,idimp,npmax,nb
     1lok,nxv,nypmx,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRMJPOST22
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),PGRMJPOST22,nargs,part(1,npo,1),amup(1
     1,1,1,1,i),nps,noff,qm,ci,idimp,npmax,nblok,nxv,nypmx)
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
      call PGRMJPOST22(part(1,npo,1),amu,npp,noff,qm,ci,idimp,npmax,nblo
     1k,nxv,nypmx)
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
      subroutine MPGSRMJPOST22(part,amu,npp,nps,noff,qm,ci,idimp,npmax,n
     1blok,nxv,nxyp,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSRMJPOST22
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),PGSRMJPOST22,nargs,part(1,npo,1),amup(
     11,1,1,i),nps,noff,qm,ci,idimp,npmax,nblok,nxv,nxyp)
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
      call PGSRMJPOST22(part(1,npo,1),amu,npp,noff,qm,ci,idimp,npmax,nbl
     1ok,nxv,nxyp)
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
      subroutine MPGRMJPOST22L(part,amu,npp,nps,noff,qm,ci,idimp,npmax,n
     1blok,nxv,nypmx,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGRMJPOST22L
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),PGRMJPOST22L,nargs,part(1,npo,1),amup(
     11,1,1,1,i),nps,noff,qm,ci,idimp,npmax,nblok,nxv,nypmx)
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
      call PGRMJPOST22L(part(1,npo,1),amu,npp,noff,qm,ci,idimp,npmax,nbl
     1ok,nxv,nypmx)
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
      subroutine MPGSRMJPOST22L(part,amu,npp,nps,noff,qm,ci,idimp,npmax,
     1nblok,nxv,nxyp,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), amu(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension amup(2,nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l, m
      external PGSRMJPOST22L
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),PGSRMJPOST22L,nargs,part(1,npo,1),amup
     1(1,1,1,i),nps,noff,qm,ci,idimp,npmax,nblok,nxv,nxyp)
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
      call PGSRMJPOST22L(part(1,npo,1),amu,npp,noff,qm,ci,idimp,npmax,nb
     1lok,nxv,nxyp)
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
      subroutine MPGRDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm,qb
     1m,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci
      real cup, dcup, amup
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
      external PGRDCJPOST2
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGRDCJPOST2,nargs,part(1,npo,1),fxy,bx
     1y,nps,noff,cup(1,1,1,1,i),dcup(1,1,1,1,i),amup(1,1,1,1,i),qm,qbm,d
     2t,ci,idimp,npmax,nblok,nxv,nypmx)
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
      call PGRDCJPOST2(part(1,npo,1),fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,
     1dt,ci,idimp,npmax,nblok,nxv,nypmx)
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
      subroutine MPGSRDCJPOST2(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm,q
     1bm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci
      real cup, dcup, amup
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
      external PGSRDCJPOST2
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSRDCJPOST2,nargs,part(1,npo,1),fxy,b
     1xy,nps,noff,cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci,
     2idimp,npmax,nblok,nxv,nxyp)
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
      call PGSRDCJPOST2(part(1,npo,1),fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm
     1,dt,ci,idimp,npmax,nblok,nxv,nxyp)
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
      subroutine MPGRDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm,q
     1bm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr
     2)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci
      real cup, dcup, amup
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
      external PGRDCJPOST2L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGRDCJPOST2L,nargs,part(1,npo,1),fxy,b
     1xy,nps,noff,cup(1,1,1,1,i),dcup(1,1,1,1,i),amup(1,1,1,1,i),qm,qbm,
     2dt,ci,idimp,npmax,nblok,nxv,nypmx)
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
      call PGRDCJPOST2L(part(1,npo,1),fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm
     1,dt,ci,idimp,npmax,nblok,nxv,nypmx)
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
      subroutine MPGSRDCJPOST2L(part,fxy,bxy,npp,nps,noff,cu,dcu,amu,qm,
     1qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr
     2)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci
      real cup, dcup, amup
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
      external PGSRDCJPOST2L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSRDCJPOST2L,nargs,part(1,npo,1),fxy,
     1bxy,nps,noff,cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci
     2,idimp,npmax,nblok,nxv,nxyp)
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
      call PGSRDCJPOST2L(part(1,npo,1),fxy,bxy,npp,noff,cu,dcu,amu,qm,qb
     1m,dt,ci,idimp,npmax,nblok,nxv,nxyp)
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
      subroutine MPGRDCJPOST22(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm,qb
     1m,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci
      real cup, dcup, amup
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
      external PGRDCJPOST22
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGRDCJPOST22,nargs,part(1,npo,1),fxy,b
     1z,nps,noff,cup(1,1,1,1,i),dcup(1,1,1,1,i),amup(1,1,1,1,i),qm,qbm,d
     2t,ci,idimp,npmax,nblok,nxv,nypmx)
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
      call PGRDCJPOST22(part(1,npo,1),fxy,bz,npp,noff,cu,dcu,amu,qm,qbm,
     1dt,ci,idimp,npmax,nblok,nxv,nypmx)
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
      subroutine MPGSRDCJPOST22(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm,q
     1bm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci
      real cup, dcup, amup
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
      external PGSRDCJPOST22
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSRDCJPOST22,nargs,part(1,npo,1),fxy,
     1bz,nps,noff,cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci,
     2idimp,npmax,nblok,nxv,nxyp)
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
      call PGSRDCJPOST22(part(1,npo,1),fxy,bz,npp,noff,cu,dcu,amu,qm,qbm
     1,dt,ci,idimp,npmax,nblok,nxv,nxyp)
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
      subroutine MPGRDCJPOST22L(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm,q
     1bm,dt,ci,idimp,npmax,nblok,nxv,nypmx,cup,dcup,amup,idtask,nmt,ierr
     2)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci
      real cup, dcup, amup
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
      external PGRDCJPOST22L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGRDCJPOST22L,nargs,part(1,npo,1),fxy,
     1bz,nps,noff,cup(1,1,1,1,i),dcup(1,1,1,1,i),amup(1,1,1,1,i),qm,qbm,
     2dt,ci,idimp,npmax,nblok,nxv,nypmx)
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
      call PGRDCJPOST22L(part(1,npo,1),fxy,bz,npp,noff,cu,dcu,amu,qm,qbm
     1,dt,ci,idimp,npmax,nblok,nxv,nypmx)
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
      subroutine MPGSRDCJPOST22L(part,fxy,bz,npp,nps,noff,cu,dcu,amu,qm,
     1qbm,dt,ci,idimp,npmax,nblok,nxv,nxyp,cup,dcup,amup,idtask,nmt,ierr
     2)
c parallel multitasking momentum flux, acceleration, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci
      real cup, dcup, amup
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
      external PGSRDCJPOST22L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PGSRDCJPOST22L,nargs,part(1,npo,1),fxy
     1,bz,nps,noff,cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci
     2,idimp,npmax,nblok,nxv,nxyp)
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
      call PGSRDCJPOST22L(part(1,npo,1),fxy,bz,npp,noff,cu,dcu,amu,qm,qb
     1m,dt,ci,idimp,npmax,nblok,nxv,nxyp)
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
