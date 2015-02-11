c 2d PIC parallel multi-tasking library for pushing particles
c and depositing charge
c mppush2lib.f contains multi-tasking procedures to process particles:
c MPGPOST2 multi-tasking wrapper for PGPOST2
c MPGSPOST2 multi-tasking wrapper for PGSPOST2
c MPGSOST2X multi-tasking wrapper for PGSOST2X
c MPGPOST2L multi-tasking wrapper for PGPOST2L
c MPGSPOST2L multi-tasking wrapper for PGSPOST2L
c MPGSOST2XL multi-tasking wrapper for PGSOST2XL
c MPGPUSH2 multi-tasking wrapper for PGPUSH2
c MPGSPUSH2 multi-tasking wrapper for PGSPUSH2
c MPGPUSH2L multi-tasking wrapper for PGPUSH2L
c MPGSPUSH2L multi-tasking wrapper for PGSPUSH2L
c MPSORTP2Y multi-tasking wrapper for PSORTP2Y
c MPSORTP2YL multi-tasking wrapper for PSORTP2YL
c MPDSORTP2Y multi-tasking wrapper for PDSORTP2Y
c MPDSORTP2YL multi-tasking wrapper for PDSORTP2YL
c MPPUSH2ZF multi-tasking wrapper for PPUSH2ZF
c MPGCJPOST2 multi-tasking wrapper for PGCJPOST2
c MPGCJPOST2L multi-tasking wrapper for PGCJPOST2L
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: december 5, 2009
c-----------------------------------------------------------------------
      subroutine MPDOST2(part,q,npp,nps,noff,qm,nx,idimp,npmax,nblok,nxv
     1,nypmx,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension qp(nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l
      external PDOST2
      data nargs /11/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 k = 1, nypmx
      do 20 j = 1, nx
      qp(j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PDOST2,nargs,part(1,npo,1),qp(1,1,1,i)
     1,nps,noff,qm,nx,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining charge
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PDOST2(part(1,npo,1),q,npp,noff,qm,nx,idimp,npmax,nblok,nxv,n
     1ypmx)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 100 l = 1, nblok
      do 90 k = 1, nypmx
      do 80 j = 1, nx
      q(j,k,l) = q(j,k,l) + qp(j,k,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGPOST2(part,q,npp,nps,noff,qm,idimp,npmax,nblok,nxv,n
     1ypmx,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension qp(nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l
      external PGPOST2
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 k = 1, nypmx
      do 20 j = 1, nxv
      qp(j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGPOST2,nargs,part(1,npo,1),qp(1,1,1,i
     1),nps,noff,qm,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining charge
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGPOST2(part(1,npo,1),q,npp,noff,qm,idimp,npmax,nblok,nxv,nyp
     1mx)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 100 l = 1, nblok
      do 90 k = 1, nypmx
      do 80 j = 1, nxv
      q(j,k,l) = q(j,k,l) + qp(j,k,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSPOST2(part,q,npp,nps,noff,qm,idimp,npmax,nblok,nxv,
     1nxyp,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension qp(nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSPOST2
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, nxyp
      qp(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSPOST2,nargs,part(1,npo,1),qp(1,1,i)
     1,nps,noff,qm,idimp,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PGSPOST2(part(1,npo,1),q,npp,noff,qm,idimp,npmax,nblok,nxv,nx
     1yp)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 l = 1, nblok
      do 70 j = 1, nxyp
      q(j,l) = q(j,l) + qp(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPSOST2X(part,q,npp,nps,noff,nn,amxy,qm,nx,idimp,npmax,
     1nblok,nxv,nxvyp,npd,nine,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, amxy, qm, qp
      integer npp, nps, noff, nn
      integer nx, idimp, npmax, nblok, nxv, nxvyp, npd, nine
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(nine,npd,nblok,nmt+1), amxy(nine,npd,nblok,nmt+1)
      dimension qp(nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PSOST2X
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, nxvyp
      qp(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PSOST2X,nargs,part(1,npo,1),qp(1,1,i),
     1nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,nx,idimp,npmax,nblok,nxv
     2,nxvyp,npd,nine)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PSOST2X(part(1,npo,1),q,npp,noff,nn,amxy,qm,nx,idimp,npmax,nb
     1lok,nxv,nxvyp,npd,nine)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 l = 1, nblok
      do 70 j = 1, nxvyp
      q(j,l) = q(j,l) + qp(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSOST2X(part,q,npp,nps,noff,nn,amxy,qm,idimp,npmax,nb
     1lok,nxv,nxvyp,npd,nine,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, amxy, qm, qp
      integer npp, nps, noff, nn
      integer idimp, npmax, nblok, nxv, nxvyp, npd, nine
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(nine,npd,nblok,nmt+1), amxy(nine,npd,nblok,nmt+1)
      dimension qp(nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSOST2X
      data nargs /14/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, nxvyp
      qp(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSOST2X,nargs,part(1,npo,1),qp(1,1,i)
     1,nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,idimp,npmax,nblok,nxv,n
     2xvyp,npd,nine)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PGSOST2X(part(1,npo,1),q,npp,noff,nn,amxy,qm,idimp,npmax,nblo
     1k,nxv,nxvyp,npd,nine)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 l = 1, nblok
      do 70 j = 1, nxvyp
      q(j,l) = q(j,l) + qp(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPDOST2L(part,q,npp,nps,noff,qm,nx,idimp,npmax,nblok,nx
     1v,nypmx,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension qp(nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l
      external PDOST2L
      data nargs /11/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 k = 1, nypmx
      do 20 j = 1, nx
      qp(j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PDOST2L,nargs,part(1,npo,1),qp(1,1,1,i
     1),nps,noff,qm,nx,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining charge
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PDOST2L(part(1,npo,1),q,npp,noff,qm,nx,idimp,npmax,nblok,nxv,
     1nypmx)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 100 l = 1, nblok
      do 90 k = 1, nypmx
      do 80 j = 1, nx
      q(j,k,l) = q(j,k,l) + qp(j,k,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGPOST2L(part,q,npp,nps,noff,qm,idimp,npmax,nblok,nxv,
     1nypmx,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension qp(nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l
      external PGPOST2L
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 40 l = 1, nblok
      nps(l) = npt
      do 30 k = 1, nypmx
      do 20 j = 1, nxv
      qp(j,k,l,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGPOST2L,nargs,part(1,npo,1),qp(1,1,1,
     1i),nps,noff,qm,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining charge
      npo = npl + 1
      do 60 l = 1, nblok
      npp(l) = npp(l) - npl
   60 continue
      call PGPOST2L(part(1,npo,1),q,npp,noff,qm,idimp,npmax,nblok,nxv,ny
     1pmx)
      do 70 l = 1, nblok
      npp(l) = npp(l) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 100 l = 1, nblok
      do 90 k = 1, nypmx
      do 80 j = 1, nxv
      q(j,k,l) = q(j,k,l) + qp(j,k,l,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSPOST2L(part,q,npp,nps,noff,qm,idimp,npmax,nblok,nxv
     1,nxyp,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nxyp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension qp(nxyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSPOST2L
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, nxyp
      qp(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSPOST2L,nargs,part(1,npo,1),qp(1,1,i
     1),nps,noff,qm,idimp,npmax,nblok,nxv,nxyp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PGSPOST2L(part(1,npo,1),q,npp,noff,qm,idimp,npmax,nblok,nxv,n
     1xyp)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 l = 1, nblok
      do 70 j = 1, nxyp
      q(j,l) = q(j,l) + qp(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPSOST2XL(part,q,npp,nps,noff,nn,amxy,qm,nx,idimp,npmax
     1,nblok,nxv,nxvyp,npd,ifour,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, amxy, qm, qp
      integer npp, nps, noff, nn
      integer nx, idimp, npmax, nblok, nxv, nxvyp, npd, ifour
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(ifour,npd,nblok,nmt+1), amxy(ifour,npd,nblok,nmt+1)
      dimension qp(nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PSOST2XL
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, nxvyp
      qp(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PSOST2XL,nargs,part(1,npo,1),qp(1,1,i)
     1,nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,nx,idimp,npmax,nblok,nx
     2v,nxvyp,npd,ifour)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PSOST2XL(part(1,npo,1),q,npp,noff,nn,amxy,qm,nx,idimp,npmax,n
     1blok,nxv,nxvyp,npd,ifour)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 l = 1, nblok
      do 70 j = 1, nxvyp
      q(j,l) = q(j,l) + qp(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSOST2XL(part,q,npp,nps,noff,nn,amxy,qm,idimp,npmax,n
     1blok,nxv,nxvyp,npd,ifour,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, amxy, qm, qp
      integer npp, nps, noff, nn
      integer idimp, npmax, nblok, nxv, nxvyp, npd, ifour
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), q(nxvyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension nn(ifour,npd,nblok,nmt+1), amxy(ifour,npd,nblok,nmt+1)
      dimension qp(nxvyp,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, l
      external PGSOST2XL
      data nargs /14/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 l = 1, nblok
      nps(l) = npt
      do 20 j = 1, nxvyp
      qp(j,l,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSOST2XL,nargs,part(1,npo,1),qp(1,1,i
     1),nps,noff,nn(1,1,1,i+1),amxy(1,1,1,i+1),qm,idimp,npmax,nblok,nxv,
     2nxvyp,npd,ifour)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      call PGSOST2XL(part(1,npo,1),q,npp,noff,nn,amxy,qm,idimp,npmax,nbl
     1ok,nxv,nxvyp,npd,ifour)
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 l = 1, nblok
      do 70 j = 1, nxvyp
      q(j,l) = q(j,l) + qp(j,l,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPUSH2(part,fx,fy,npp,nps,noff,qbm,dt,ek,nx,idimp,npma
     1x,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PPUSH2
      data nargs /14/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PPUSH2,nargs,part(1,npo,1),fx,fy,nps,n
     1off,qbm,dt,ekp(i),nx,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PPUSH2(part(1,npo,1),fx,fy,npp,noff,qbm,dt,ek,nx,idimp,npmax,
     1nblok,nxv,nypmx)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGPUSH2(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idimp,np
     1max,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGPUSH2
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGPUSH2,nargs,part(1,npo,1),fxy,nps,no
     1ff,qbm,dt,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGPUSH2(part(1,npo,1),fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npma
     1x,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSPUSH2(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idimp,n
     1pmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSPUSH2
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSPUSH2,nargs,part(1,npo,1),fxy,nps,n
     1off,qbm,dt,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSPUSH2(part(1,npo,1),fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npm
     1ax,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPUSH2L(part,fx,fy,npp,nps,noff,qbm,dt,ek,nx,idimp,npm
     1ax,nblok,nxv,nypmx,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, nblok, nxv, nypmx, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fx(nxv,nypmx,nblok), fy(nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PPUSH2L
      data nargs /14/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PPUSH2L,nargs,part(1,npo,1),fx,fy,nps,
     1noff,qbm,dt,ekp(i),nx,idimp,npmax,nblok,nxv,nypmx)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PPUSH2L(part(1,npo,1),fx,fy,npp,noff,qbm,dt,ek,nx,idimp,npmax
     1,nblok,nxv,nypmx)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idimp,n
     1pmax,nblok,nxv,nypmx,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGPUSH2L
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGPUSH2L,nargs,part(1,npo,1),fxy,nps,n
     1off,qbm,dt,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGPUSH2L(part(1,npo,1),fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npm
     1ax,nblok,nxv,nypmx,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSPUSH2L(part,fxy,npp,nps,noff,qbm,dt,ek,nx,ny,idimp,
     1npmax,nblok,nxv,nxyp,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxyp,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PGSPUSH2L
      data nargs /15/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSPUSH2L,nargs,part(1,npo,1),fxy,nps,
     1noff,qbm,dt,ekp(i),nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PGSPUSH2L(part(1,npo,1),fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,np
     1max,nblok,nxv,nxyp,ipbc)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPSORTP2Y(part,pt,ip,npic,npp,nps,noff,nyp,idimp,npmax,
     1nblok,nypm1,npicp,idtask,nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, pt
      integer ip, npic, npp, nps, noff, nyp
      integer idimp, npmax, nblok, nypm1
      integer npicp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), pt(npmax,nblok)
      dimension ip(npmax,nblok), npic(nypm1,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok), nyp(nblok)
      dimension npicp(nypm1,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PSORTP2Y
      data nargs /11/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle sorting tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      call MP_TASKSTART(idtask(i),PSORTP2Y,nargs,part(1,npo,1),pt(npo,1)
     1,ip(npo,1),npicp(1,1,i),nps,noff,nyp,idimp,npmax,nblok,nypm1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c sort remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PSORTP2Y(part(1,npo,1),pt(npo,1),ip(npo,1),npic,npp,noff,nyp,
     1idimp,npmax,nblok,nypm1)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPSORTP2YL(part,pt,ip,npic,npp,nps,noff,nyp,idimp,npmax
     1,nblok,nypm1,npicp,idtask,nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, pt
      integer ip, npic, npp, nps, noff, nyp
      integer idimp, npmax, nblok, nypm1
      integer npicp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), pt(npmax,nblok)
      dimension ip(npmax,nblok), npic(nypm1,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok), nyp(nblok)
      dimension npicp(nypm1,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PSORTP2YL
      data nargs /11/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle sorting tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      call MP_TASKSTART(idtask(i),PSORTP2YL,nargs,part(1,npo,1),pt(npo,1
     1),ip(npo,1),npicp(1,1,i),nps,noff,nyp,idimp,npmax,nblok,nypm1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c sort remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PSORTP2YL(part(1,npo,1),pt(npo,1),ip(npo,1),npic,npp,noff,nyp
     1,idimp,npmax,nblok,nypm1)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPDSORTP2Y(parta,partb,npic,npp,nps,noff,nyp,idimp,npma
     1x,nblok,nypm1,npicp,idtask,nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real parta, partb
      integer npic, npp, nps, noff, nyp
      integer idimp, npmax, nblok, nypm1
      integer npicp, idtask, nmt, ierr
      dimension parta(idimp,npmax,nblok), partb(idimp,npmax,nblok)
      dimension npic(nypm1,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok), nyp(nblok)
      dimension npicp(nypm1,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PDSORTP2Y
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle sorting tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      call MP_TASKSTART(idtask(i),PDSORTP2Y,nargs,parta(1,npo,1),partb(1
     1,npo,1),npicp(1,1,i),nps,noff,nyp,idimp,npmax,nblok,nypm1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c sort remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PDSORTP2Y(parta(1,npo,1),partb(1,npo,1),npic,npp,noff,nyp,idi
     1mp,npmax,nblok,nypm1)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPDSORTP2YL(parta,partb,npic,npp,nps,noff,nyp,idimp,npm
     1ax,nblok,nypm1,npicp,idtask,nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real parta, partb
      integer npic, npp, nps, noff, nyp
      integer idimp, npmax, nblok, nypm1
      integer npicp, idtask, nmt, ierr
      dimension parta(idimp,npmax,nblok), partb(idimp,npmax,nblok)
      dimension npic(nypm1,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok), nyp(nblok)
      dimension npicp(nypm1,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PDSORTP2YL
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle sorting tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      call MP_TASKSTART(idtask(i),PDSORTP2YL,nargs,parta(1,npo,1),partb(
     11,npo,1),npicp(1,1,i),nps,noff,nyp,idimp,npmax,nblok,nypm1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c sort remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PDSORTP2YL(parta(1,npo,1),partb(1,npo,1),npic,npp,noff,nyp,id
     1imp,npmax,nblok,nypm1)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPUSH2ZF(part,npp,nps,dt,ek,idimp,npmax,nblok,nx,ny,ip
     1bc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, dt, ek, ekp
      integer npp, nps
      integer idimp, npmax, nblok, nx, ny, ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension npp(nblok), nps(nblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, l
      external PPUSH2ZF
      data nargs /10/
c find minimum of npp
      npo = npp(1)
      do 10 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 l = 1, nblok
      nps(l) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PPUSH2ZF,nargs,part(1,npo,1),nps,dt,ek
     1p(i),idimp,npmax,nblok,nx,ny,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 l = 1, nblok
      npp(l) = npp(l) - npl
   40 continue
      call PPUSH2ZF(part(1,npo,1),npp,dt,ek,idimp,npmax,nblok,nx,ny,ipbc
     1)
      do 50 l = 1, nblok
      npp(l) = npp(l) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGCJPOST2(part,fxy,npp,nps,noff,cu,qm,qbm,dt,idimp,npm
     1ax,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
c parallel multitasking current density deposit
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, cup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGCJPOST2
      data nargs /13/
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
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGCJPOST2,nargs,part(1,npo,1),fxy,nps,
     1noff,cup(1,1,1,1,i),qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx)
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
      call PGCJPOST2(part(1,npo,1),fxy,npp,noff,cu,qm,qbm,dt,idimp,npmax&
     &,nblok,nxv,nypmx)
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
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGCJPOST2L(part,fxy,npp,nps,noff,cu,qm,qbm,dt,idimp,np
     1max,nblok,nxv,nypmx,cup,idtask,nmt,ierr)
c parallel multitasking current density deposit
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, cup
      integer npp, nps, noff
      integer idimp, npmax, nblok, nxv, nypmx
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension fxy(2,nxv,nypmx,nblok), cu(2,nxv,nypmx,nblok)
      dimension npp(nblok), nps(nblok), noff(nblok)
      dimension cup(2,nxv,nypmx,nblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGCJPOST2L
      data nargs /13/
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
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGCJPOST2L,nargs,part(1,npo,1),fxy,nps
     1,noff,cup(1,1,1,1,i),qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx)
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
      call PGCJPOST2L(part(1,npo,1),fxy,npp,noff,cu,qm,qbm,dt,idimp,npma
     1x,nblok,nxv,nypmx)
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
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
