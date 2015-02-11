!-----------------------------------------------------------------------
!
      module mp2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library mp2lib.f
! mp2mod.f contains multi-tasking interface procedures for
!          communications with 1d partitions:
!         defines module mp2d
! pmove => impmove2 moves particles to appropriate processor non-uniform
!          1d partition boundaries in 2d code.
!          calls MPMOVE2, or MPXMOV2
! pmove => impdmove2 moves particles to appropriate processor
!          non-uniform 1d partition boundaries in 2d code, and returns
!          load imbalance.
!          calls MPMOVE2, or MPXMOV2
! pmoves => impmoves2 moves particles to appropriate processor
!           non-uniform 1d partition boundaries in 2d code.  Optimized
!           requiring that maximum number of particle passes must be
!           known.
!           calls MPMOVES2, or MPXMOVS2
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 12, 2009
!
      use p0d
      use mp0d, only: ntasks
      use p2d, only: dcomp, pcguard, pcguardp, paguard, pamcguard,      &
     &paguardp, pncguardp, pnaguardp, pfmove, repart, fnoff, dblsin,    &
     &dblcos, hafdbl, zdbl, plsum, plbcast, writebf, readbf, wrdata,    &
     &rddata
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PPINIT, PPID, PPEXIT, HARTBEAT
      public :: get_funit, pwtimer, plsum, plmax, plscan, plbcast
      public :: writebf, readbf, wtimer, wrdata, rddata
      public :: dcomp, pmove, pcguard, pcguardp
      public :: paguard, pamcguard, paguardp, pncguardp, pnaguardp
      public :: pfmove, repart
      public :: fnoff, dblsin, dblcos, hafdbl, zdbl, pmoves
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MPMOVE2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,i&
     &hole,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntma&
     &x,info,jssp,idtask,nmt,ierr)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax, nmt, ierr
         real :: th
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp, npq
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         integer, dimension(idps,nblok,nmt) :: jssp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPXMOV2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,i&
     &hole,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntma&
     &x,maskp,info,jssp,idtask,nmt,ierr)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax, nmt, ierr
         real :: th
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(npmax,nblok) :: maskp
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp, npq
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         integer, dimension(idps,nblok,nmt) :: jssp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPMOVES2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,&
     &ihole,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntm&
     &ax,info,jssp,idtask,nmt,ierr)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax, nmt, ierr
         real :: th
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp, npq
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         integer, dimension(idps,nblok,nmt) :: jssp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPXMOVS2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,&
     &ihole,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntm&
     &ax,maskp,info,jssp,idtask,nmt,ierr)
         implicit none
         integer :: ny, kstrt, nvp, idimp, npmax, nblok, idps
         integer :: nbmax, ntmax, nmt, ierr
         real :: th
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(npmax,nblok) :: maskp
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp, npq
         real, dimension(idimp,nbmax,nblok) :: sbufl, sbufr
         real, dimension(idimp,nbmax,nblok) :: rbufl, rbufr
         integer, dimension(idps,nblok) :: jsl, jsr, jss
         integer, dimension(ntmax,nblok) :: ihole
         integer, dimension(7) :: info
         integer, dimension(idps,nblok,nmt) :: jssp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface pmove
         module procedure impmove2
         module procedure impdmove2
      end interface
!
      interface pmoves
         module procedure impmoves2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains  
!
         subroutine impmove2(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vt,&
     &ierr)
! multi-tasking particle manager
! non-uniform 1d partition boundaries in 2d code
         implicit none
         integer :: ny, kstrt, nvp, nbmax, vt, ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, idps, ntmax, nmt
         integer, dimension(ntasks) :: idtask
         integer, dimension(size(edges,1),size(edges,2),ntasks) :: jssp
         integer, dimension(size(npp)) :: npq
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         real, dimension(size(part,1),nbmax,size(part,3)) :: sbufl
         real, dimension(size(part,1),nbmax,size(part,3)) :: sbufr
         real, dimension(size(part,1),nbmax,size(part,3)) :: rbufl
         real, dimension(size(part,1),nbmax,size(part,3)) :: rbufr
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(2*nbmax,size(part,3)) :: ihole
         integer, dimension(7) :: info
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
         nmt = ntasks
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call MPXMOV2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ihol&
     &e,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,m&
     &askp,info,jssp,idtask,nmt,ierr)
         else
            call MPMOVE2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ihol&
     &e,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,i&
     &nfo,jssp,idtask,nmt,ierr)
         endif
         if ((ierr /= 0) .or. (info(1) /= 0)) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
         ierr = info(1)
         end subroutine impmove2
!
         subroutine impdmove2(part,edges,npp,anpav,pibal,tmove,ny,kstrt,&
     &nvp,nbmax,vt,ierr)
! multi-tasking particle manager
! non-uniform 1d partition boundaries in 2d code
! returns load imbalance
         implicit none
         integer :: ny, kstrt, nvp, nbmax, vt, ierr
         real :: anpav, pibal
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, idps, ntmax, nmt
         integer, dimension(ntasks) :: idtask
         integer, dimension(size(edges,1),size(edges,2),ntasks) :: jssp
         integer, dimension(size(npp)) :: npq
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         real, dimension(size(part,1),nbmax,size(part,3)) :: sbufl
         real, dimension(size(part,1),nbmax,size(part,3)) :: sbufr
         real, dimension(size(part,1),nbmax,size(part,3)) :: rbufl
         real, dimension(size(part,1),nbmax,size(part,3)) :: rbufr
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(2*nbmax,size(part,3)) :: ihole
         integer, dimension(7) :: info
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
         nmt = ntasks
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call MPXMOV2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ihol&
     &e,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,m&
     &askp,info,jssp,idtask,nmt,ierr)
         else
            call MPMOVE2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ihol&
     &e,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,i&
     &nfo,jssp,idtask,nmt,ierr)
         endif
         if ((ierr /= 0) .or. (info(1) /= 0)) then
            call MP_END
            call PPEXIT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
! calculate percent imbalance
         anpav = real(info(7))/real(nvp)
         if (anpav > 0.0) then
            pibal = max(real(info(2))-anpav,anpav-real(info(3)))/anpav
         endif
         ierr = info(1)
         end subroutine impdmove2
!
         subroutine impmoves2(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vt&
     &,ierr)
! multi-tasking particle manager
! non-uniform 1d partition boundaries in 2d code
         implicit none
         integer :: ny, kstrt, nvp, nbmax, vt, ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, idps, ntmax, nmt
         integer, dimension(ntasks) :: idtask
         integer, dimension(size(edges,1),size(edges,2),ntasks) :: jssp
         integer, dimension(size(npp)) :: npq
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
         real, dimension(size(part,1),nbmax,size(part,3)) :: sbufl
         real, dimension(size(part,1),nbmax,size(part,3)) :: sbufr
         real, dimension(size(part,1),nbmax,size(part,3)) :: rbufl
         real, dimension(size(part,1),nbmax,size(part,3)) :: rbufr
         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
         integer, dimension(size(edges,1),size(edges,2)) :: jss
         integer, dimension(2*nbmax,size(part,3)) :: ihole
         integer, dimension(7) :: info
         real :: tf, th
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
         nmt = ntasks
! set maximum number of passes
         info(5) = 1
! initialize timer
         call wtimer(tf,dtime,-1)
         if (vt.eq.1) then
            call MPXMOVS2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,&
     &maskp,info,jssp,idtask,nmt,ierr)
         else
            call MPMOVES2(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho&
     &le,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,&
     &info,jssp,idtask,nmt,ierr)
         endif
         if ((ierr /= 0) .or. (info(1) /= 0)) then
            call MP_END
            call PPABORT
            stop
         endif
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
         ierr = info(1)
         end subroutine impmoves2
!
      end module mp2d
