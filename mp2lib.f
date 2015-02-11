c-----------------------------------------------------------------------
c 2d parallel PIC multi-tasking library for MPI communications
c mp2lib.f contains basic communications multi-tasking procedures for 2d
c          code with 1d partitions:
c MPMOVE2 multi-tasking wrapper for PMOVEH2 plus PMOVES2
c MPXMOV2 multi-tasking wrapper for PMOVEHX2 plus PMOVES2
c MPMOVES2 multi-tasking wrapper for PMOVEH2 plus PMOVESS2
c MPXMOVS2 multi-tasking wrapper for PMOVEHX2 plus PMOVESS2
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: december 11, 2009
c-----------------------------------------------------------------------
      subroutine MPMOVE2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ihol
     1e,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,i
     2nfo,jssp,idtask,nmt,ierr)
c multi-tasking particle manager
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, npq, ihole, jsr, jsl, jss, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      integer jssp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok), npq(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
      dimension jssp(idps,nblok,nmt), idtask(nmt)
c local data
      integer i, j, l, npr, nps, nter, nargs, npt, npl, npo
      integer mpo, mpt, mpl, mps
      integer ibflg, iwork
      real tf
      double precision dtime
      dimension ibflg(2), iwork(2)
      external PMOVEH2
      data nargs /10/
      nter = 0
      do 10 j = 1, 7
      info(j) = 0
   10 continue
      th = 0.0
c debugging section: count total number of particles before move
      npr = 0
      do 20 l = 1, nblok
      npr = npr + npp(l)
   20 continue
c find minimum of npp
      npo = npp(1)
      do 30 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   30 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
c find ihole partitions
      mpt = ntmax/(nmt + 1)
      mpl = mpt*nmt
      ierr = 0
c find outgoing particles
   40 nter = info(4)
      call PWTIMERA(-1,tf,dtime)
c start tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
      mpo = mpt*(i - 1) + 1
      do 50 l = 1, nblok
      npq(l) = npt
   50 continue
      call MP_TASKSTART(idtask(i),PMOVEH2,nargs,part(1,npo,1),edges,npq,
     1ihole(mpo,1),jssp(1,1,i),idimp,npt,nblok,idps,mpt)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c finish finding outgoing particles
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      npl = npmax - npl
      mpo = mpl + 1
      mps = ntmax - mpl
      call PMOVEH2(part(1,npo,1),edges,npp,ihole(mpo,1),jss,idimp,npl,nb
     1lok,idps,mps)
      npl = npmax - npl
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c compress ihole array
      if (i.eq.1) go to 110
      npo = npt*(i - 1)
      mpo = mpt*(i - 1)
      do 100 l = 1, nblok
      do 90 j = 1, jssp(1,l,i)
      ihole(jssp(1,l,i-1)+j,l) = ihole(j+mpo,l) + npo
   90 continue
      jssp(1,l,i) = jssp(1,l,i) + jssp(1,l,i-1)
      jssp(2,l,i) = max0(jssp(2,l,i),jssp(2,l,i-1))
  100 continue
  110 continue
      if (nmt.gt.0) then
         do 130 l = 1, nblok
         do 120 j = 1, jss(1,l)
         ihole(jssp(1,l,nmt)+j,l) = ihole(j+mpl,l) + npl
  120    continue
         jss(1,l) = jss(1,l) + jssp(1,l,nmt)
         jss(2,l) = max0(jss(2,l),jssp(2,l,nmt))
  130    continue
      endif
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVES2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl,
     1jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
c particle overflow error
      if (info(1).gt.0) return
c buffer overflowed and more particles remain to be checked
      if (info(4).gt.nter) go to 40
c iteration overflow
      if (info(1).lt.0) go to 150
c debugging section: count total number of particles after move
      nps = 0
      do 140 l = 1, nblok
      nps = nps + npp(l)
  140 continue
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      info(6) = ibflg(1)
      info(7) = ibflg(2)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         info(1) = 1
      endif
c information
  150 nter = info(4)
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPXMOV2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ihol
     1e,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,m
     2askp,info,jssp,idtask,nmt,ierr)
c multi-tasking particle manager
c optimized for vector processor
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, npq, ihole, jsr, jsl, jss, maskp, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      integer jssp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), maskp(npmax,nblok)
      dimension edges(idps,nblok), npp(nblok), npq(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
      dimension jssp(idps,nblok,nmt), idtask(nmt)
c local data
      integer i, j, l, npr, nps, nter, nargs, npt, npl, npo
      integer mpo, mpt, mpl, mps
      integer ibflg, iwork
      real tf
      double precision dtime
      dimension ibflg(2), iwork(2)
      external PMOVEHX2
      data nargs /11/
      nter = 0
      do 10 j = 1, 7
      info(j) = 0
   10 continue
      th = 0.0
c debugging section: count total number of particles before move
      npr = 0
      do 20 l = 1, nblok
      npr = npr + npp(l)
   20 continue
c find minimum of npp
      npo = npp(1)
      do 30 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   30 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
c find ihole partitions
      mpt = ntmax/(nmt + 1)
      mpl = mpt*nmt
      ierr = 0
c find outgoing particles
   40 nter = info(4)
      call PWTIMERA(-1,tf,dtime)
c start tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
      mpo = mpt*(i - 1) + 1
      do 50 l = 1, nblok
      npq(l) = npt
   50 continue
      call MP_TASKSTART(idtask(i),PMOVEHX2,nargs,part(1,npo,1),edges,npq
     1,ihole(mpo,1),jssp(1,1,i),idimp,npt,nblok,idps,mpt,maskp(npo,1))
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c finish finding outgoing particles
      npo = npl + 1
      do 70 l = 1, nblok
      npp(l) = npp(l) - npl
   70 continue
      npl = npmax - npl
      mpo = mpl + 1
      mps = ntmax - mpl
      call PMOVEHX2(part(1,npo,1),edges,npp,ihole(mpo,1),jss,idimp,npl,n
     1blok,idps,mps,maskp(npo,1))
      npl = npmax - npl
      do 80 l = 1, nblok
      npp(l) = npp(l) + npl
   80 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c compress ihole array
      if (i.eq.1) go to 110
      npo = npt*(i - 1)
      mpo = mpt*(i - 1)
      do 100 l = 1, nblok
      do 90 j = 1, jssp(1,l,i)
      ihole(jssp(1,l,i-1)+j,l) = ihole(j+mpo,l) + npo
   90 continue
      jssp(1,l,i) = jssp(1,l,i) + jssp(1,l,i-1)
      jssp(2,l,i) = max0(jssp(2,l,i),jssp(2,l,i-1))
  100 continue
  110 continue
      if (nmt.gt.0) then
         do 130 l = 1, nblok
         do 120 j = 1, jss(1,l)
         ihole(jssp(1,l,nmt)+j,l) = ihole(j+mpl,l) + npl
  120    continue
         jss(1,l) = jss(1,l) + jssp(1,l,nmt)
         jss(2,l) = max0(jss(2,l),jssp(2,l,nmt))
  130    continue
      endif
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVES2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl,
     1jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
c particle overflow error
      if (info(1).gt.0) return
c buffer overflowed and more particles remain to be checked
      if (info(4).gt.nter) go to 40
c iteration overflow
      if (info(1).lt.0) go to 150
c debugging section: count total number of particles after move
      nps = 0
      do 140 l = 1, nblok
      nps = nps + npp(l)
  140 continue
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      info(6) = ibflg(1)
      info(7) = ibflg(2)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         info(1) = 1
      endif
c information
  150 nter = info(4)
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPMOVES2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho
     1le,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,
     2info,jssp,idtask,nmt,ierr)
c multi-tasking particle manager
c info(5) = maximum number of particle passes required
c info(5) must be set on entry
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, npq, ihole, jsr, jsl, jss, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      integer jssp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok), npq(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
      dimension jssp(idps,nblok,nmt), idtask(nmt)
c local data
      integer i, j, l, nargs, npt, npl, npo
      integer mpo, mpt, mpl, mps
      real tf
      double precision dtime
      external PMOVEH2
      data nargs /10/
      do 10 j = 1, 4
      info(j) = 0
   10 continue
      th = 0.0
c find minimum of npp
      npo = npp(1)
      do 20 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   20 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
c find ihole partitions
      mpt = ntmax/(nmt + 1)
      mpl = mpt*nmt
      ierr = 0
c find outgoing particles
      call PWTIMERA(-1,tf,dtime)
c start tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
      mpo = mpt*(i - 1) + 1
      do 30 l = 1, nblok
      npq(l) = npt
   30 continue
      call MP_TASKSTART(idtask(i),PMOVEH2,nargs,part(1,npo,1),edges,npq,
     1ihole(mpo,1),jssp(1,1,i),idimp,npt,nblok,idps,mpt)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c finish finding outgoing particles
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      npl = npmax - npl
      mpo = mpl + 1
      mps = ntmax - mpl
      call PMOVEH2(part(1,npo,1),edges,npp,ihole(mpo,1),jss,idimp,npl,nb
     1lok,idps,mps)
      npl = npmax - npl
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c compress ihole array
      if (i.eq.1) go to 90
      npo = npt*(i - 1)
      mpo = mpt*(i - 1)
      do 80 l = 1, nblok
      do 70 j = 1, jssp(1,l,i)
      ihole(jssp(1,l,i-1)+j,l) = ihole(j+mpo,l) + npo
   70 continue
      jssp(1,l,i) = jssp(1,l,i) + jssp(1,l,i-1)
      jssp(2,l,i) = max0(jssp(2,l,i),jssp(2,l,i-1))
   80 continue
   90 continue
      if (nmt.gt.0) then
         do 110 l = 1, nblok
         do 100 j = 1, jss(1,l)
         ihole(jssp(1,l,nmt)+j,l) = ihole(j+mpl,l) + npl
  100    continue
         jss(1,l) = jss(1,l) + jssp(1,l,nmt)
         jss(2,l) = max0(jss(2,l),jssp(2,l,nmt))
  110    continue
      endif
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVESS2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl
     1,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPXMOVS2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho
     1le,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,
     2maskp,info,jssp,idtask,nmt,ierr)
c multi-tasking particle manager
c info(5) = maximum number of particle passes required
c info(5) must be set on entry
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, npq, ihole, jsr, jsl, jss, maskp, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      integer jssp, idtask, nmt, ierr
      dimension part(idimp,npmax,nblok), maskp(npmax,nblok)
      dimension edges(idps,nblok), npp(nblok), npq(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
      dimension jssp(idps,nblok,nmt), idtask(nmt)
c local data
      integer i, j, l, nargs, npt, npl, npo
      integer mpo, mpt, mpl, mps
      real tf
      double precision dtime
      external PMOVEHX2
      data nargs /11/
      do 10 j = 1, 4
      info(j) = 0
   10 continue
      th = 0.0
c find minimum of npp
      npo = npp(1)
      do 20 l = 1, nblok
      if (npo.gt.npp(l)) npo = npp(l)
   20 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
c find ihole partitions
      mpt = ntmax/(nmt + 1)
      mpl = mpt*nmt
      ierr = 0
c find outgoing particles
      call PWTIMERA(-1,tf,dtime)
c start tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
      mpo = mpt*(i - 1) + 1
      do 30 l = 1, nblok
      npq(l) = npt
   30 continue
      call MP_TASKSTART(idtask(i),PMOVEHX2,nargs,part(1,npo,1),edges,npq
     1,ihole(mpo,1),jssp(1,1,i),idimp,npt,nblok,idps,mpt,maskp(npo,1))
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c finish finding outgoing particles
      npo = npl + 1
      do 50 l = 1, nblok
      npp(l) = npp(l) - npl
   50 continue
      npl = npmax - npl
      mpo = mpl + 1
      mps = ntmax - mpl
      call PMOVEHX2(part(1,npo,1),edges,npp,ihole(mpo,1),jss,idimp,npl,n
     1blok,idps,mps,maskp(npo,1))
      npl = npmax - npl
      do 60 l = 1, nblok
      npp(l) = npp(l) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c compress ihole array
      if (i.eq.1) go to 90
      npo = npt*(i - 1)
      mpo = mpt*(i - 1)
      do 80 l = 1, nblok
      do 70 j = 1, jssp(1,l,i)
      ihole(jssp(1,l,i-1)+j,l) = ihole(j+mpo,l) + npo
   70 continue
      jssp(1,l,i) = jssp(1,l,i) + jssp(1,l,i-1)
      jssp(2,l,i) = max0(jssp(2,l,i),jssp(2,l,i-1))
   80 continue
   90 continue
      if (nmt.gt.0) then
         do 110 l = 1, nblok
         do 100 j = 1, jss(1,l)
         ihole(jssp(1,l,nmt)+j,l) = ihole(j+mpl,l) + npl
  100    continue
         jss(1,l) = jss(1,l) + jssp(1,l,nmt)
         jss(2,l) = max0(jss(2,l),jssp(2,l,nmt))
  110    continue
      endif
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVESS2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl
     1,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
      return
      end
