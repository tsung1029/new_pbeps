!-----------------------------------------------------------------------
! * * * periodic 2d darwin particle simulation kernel code * * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves electrons and ions, with periodic
! electromagnetic forces obtained by solving non-radiative form of
! maxwell's equation with fast fourier transforms
! using algorithm similar to that described in
! J. Busnardo-Neto, P. L. Pritchett, A. T. Lin, and J. M. Dawson,
! J. Computational Phys. 23, 300 (1977).
! portable gcpic kernel code, using algorithm described in:
! p. c. liewer and v. k. decyk, j. computational phys. 85, 302 (1989).
! written by viktor k. decyk, ucla
! for mpi distributed memory multiprocessing macintosh computers
! copyright 2009, regents of the university of california
! update: May 1, 2013
      program pdbeps2
      use pinit2d
      use pempush2d
      use pfield2d
      use pdiag2d
      use pemsimul2d
      use mp0d, only: mpinit, ncpus
      implicit none
! idps = number of partition boundaries
! mshare = (0,1) = (no,yes) architecture is shared memory
      integer :: idps =    2, mshare =   0
! nmv = number of segments in v for velocity distribution
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      integer :: nmv = 40, vect = 0, ipbc = 1
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10
      integer :: iuvi = 20, iup = 11, iua = 15, iuj = 25, ium = 21
      integer :: iuer = 2
!     integer :: npxy, npxyb, np, npxyi, npxybi, npi
      double precision :: npxy, npxyb, np, npxyi, npxybi, npi
      integer :: nx, ny, nxh, nyh, nyv, nxe, nxeh, idimp
      integer :: nloop, nvp, nblok, npmax, npimax = 0, kyp, kxp, nypmx
      integer :: kyb, kxb, kxyb, kbmin, kblok, jbmin, jblok
      integer :: ngds, nxyh, nxhy, nx1, nypm1, nbmax
      integer :: idproc, id0, kstrt, itime, itime0, ntime, isign
      integer :: k, it, itw, ierr, iur1, iur2, irc
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0, ltime = 0.0, tloop = 0.0, ts = 0.0
      real :: we = 0.0, wf = 0.0, wm = 0.0, wke = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tdjpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tdjposti = 0.0, tsorti = 0.0
      real :: tfield = 0.0, totpush = 0.0, totpushi = 0.0
      real :: tdcjpost = 0.0, tdcjposti = 0.0
      real :: qbme, qbmi, affp, qi0, omt, q2m0, wp0, wpmax, wpmin, etx
      real :: sx = 0.0, sy = 0.0, sz = 0.0
      real :: pxe = 0.0, pye = 0.0, pze = 0.0
      real :: pxi = 0.0, pyi = 0.0, pzi = 0.0
      real :: wx, wy, wz
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      double precision :: dtime, ldtime
      real, dimension(:,:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:,:), pointer :: qe, qi
      real, dimension(:,:,:,:), pointer :: fxyze, exyze, bxyze
      real, dimension(:,:,:,:), pointer :: cu, cus, amu
      complex, dimension(:,:,:), pointer :: qt
      complex, dimension(:,:,:,:), pointer :: cut, cur, amut, fxyt, bxyt
      complex, dimension(:,:), pointer :: q2m
      complex, dimension(:,:,:), pointer :: ffc, ffe
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:,:), pointer  :: edges
      integer, dimension(:), pointer :: nyp, noff
      integer, dimension(:), pointer :: npp, nppi, nps
      real, dimension(:,:), pointer :: pt
      integer, dimension(:,:), pointer :: ip, npic
      real, dimension(:,:,:), pointer :: sfield
      complex, dimension(:,:,:), pointer :: sfieldt, dent, pott
      real, dimension(:,:,:,:), pointer :: vfield
      complex, dimension(:,:,:,:), pointer :: vfieldt, vpott, vcurt
      real, dimension(:,:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
! wtot = total energy
      real, dimension(7) :: wtot
! time = timing array
      real, dimension(2) :: tfft = 0.0, time = 0.0
      real, dimension(2) :: tmove = 0.0, tmovi = 0.0
! msg = heartbeat array
      double precision, dimension(13) :: msg
      character(len=10) :: cdrun
      character(len=32) :: fname
      character(len=12) :: label
      integer, external :: NDIAN, NDPREC, IDPREC
  991 format (' T = ',i7)
! get unit number for error file
      iuer = get_funit(iuer)
! nvp = number of real or virtual processors
! initialize for parallel processing
      call PPINIT(idproc,id0,nvp)
      kstrt = idproc + 1
! read namelist
      if (id0==0) then
         iuin = get_funit(iuin)
         open(unit=iuin,file='pinput2',form='formatted',status='old')
         read (iuin,pinput2)
      endif
! override input data
      idcode = 3
      psolve = 1
! broadcast namelist to other nodes
      call sendnml()
! set monitor flag
      call SET_MON(monitor)
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
! text output file
      if (id0==0) then
         iuot = get_funit(iuot)
         fname = 'poutput2.'//cdrun
         open(unit=iuot,file=trim(fname),form='formatted',status=       &
     &'replace')
      endif
! np = total number of electrons in simulation
!     npxy = npx*npy; npxyb = npxb*npyb; np = npxy + npxyb
      npxy = dble(npx)*dble(npy); npxyb = dble(npxb)*dble(npyb)
      np = npxy + npxyb
! npi = total number of ions in simulation
!     npxyi = npxi*npyi; npxybi = npxbi*npybi; npi = npxyi + npxybi
      npxyi = dble(npxi)*dble(npyi); npxybi = dble(npxbi)*dble(npybi)
      npi = npxyi + npxybi
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = max(1,ny/2)
      nyv = ny + 2; nxe = nx + 4
! kyp = number of complex grids in each field partition in y direction
! nypmx = maximum size of particle partition, including guard cells.
      kyp = (ny - 1)/nvp + 1; nypmx = kyp + 3
! ngds = number of guard cells
      ngds = 3*((idps - 1)/2 + 1)
!     ax = .866025; ay = .866025; az = .866025
      if (inorder==LINEAR) then
         ax = .912871; ay = .912871
         nxe = nx + 2; ; nypmx = kyp + 1
         ngds = (idps - 1)/2 + 1
      endif
      nxeh = nxe/2
! check if too many processors
      if (nvp > ny) then
         write (2,*) 'Too many processors requested, ny, nvp=', ny, nvp
         call PPEXIT
         stop
      endif
! initialize for multiprocessing
      ntasks = mpinit(sntasks)
      call MP_SETSTACK(262144)
      if (dopt==VECTOR) vect = 1
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! nblok = number of particle partitions
      nblok = 1 + mshare*(nvp - 1)
! npmax = maximum number of electrons in each partition
      npmax = (np/nvp)*1.25
      if (movion==1) npimax = (npi/nvp)*1.25
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! kyb = number of processors in y
! kxb = number of processors in x
      kyb = ny/kyp; kxb = nxh/kxp
! kxyb = maximum(kxb,kyb)
      kxyb = max(kxb,kyb)
! kblok = number of field partitions in y direction
      kbmin = 1 + (1 - mshare)*(kxyb/kxb - 1)
      kblok = 1 + mshare*(ny/kyp - 1)
! jblok = number of field partitions in x direction
      jbmin = 1 + (1 - mshare)*(kxyb/kyb - 1)
      jblok = 1 + mshare*(nxh/kxp - 1)
! nxyh = maximum(nx,ny)/2
      nxyh = max(nx,ny)/2
! nxhy = maximum(nx/2,ny)
      nxhy = max(nxh,ny)
! dimensions for index and sorting arrays
      nx1 = nx + 1; nypm1 = kyp + 1
! nbmax = size of buffer for passing particles between processors
      nbmax = 1 + (2*(npxy*vty + npxyb*vtdy) + 1.4*npxyb*abs(vdy))*dt/ny
!     nbmax = 2*nbmax
! idimp = dimension of phase space = 4 or 5
      idimp = 2 + ndim
      if (movion==1) then
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
      endif
! initialize time constants
      itime0 = 0
      itime = itime0
      ntime = itime + itime0
! diagnostic information needed by diagnostic nodes
! set default diagnostic file names
      if (ntd > 0) fdname = 'pdenk2.'//cdrun
      if (ntp > 0) fpname = 'ppotk2.'//cdrun
      if (nta > 0) faname = 'pvpotk2.'//cdrun
      if (ntj > 0) fjname = 'pvcurk2.'//cdrun
! energy time history
      if (ndw > 0) then
         allocate(wt((nloop-1)/ndw-(itime0/ndw)+1,7))
         itw = 0
      endif
! open restart files
      if (id0==0) then
         if (nustrt==0) then
            call restart_open(nustrt,ntr,idrun0,iur1,iur2,iuer)
         else
            call restart_open(nustrt,ntr,idrun,iur1,iur2,iuer)
         endif
      endif
! open graphics device
      call GROPEN
      call SETNPLT(nplot,irc)
      call STPALIT(idpal)
!
! diagnostic nodes have special processing
      if (idproc < 0) call diag2nodes
!
! part(1,n,l) = position x of particle n in partition l
! part(2,n,l) = position y of particle n in partition l
! part(3,n,l) = velocity vx of particle n in partition l
! part(4,n,l) = velocity vy of particle n in partition l
! part(5,n,l) = velocity vz of particle n in partition l
      allocate(part(idimp,npmax,nblok))
! maskp = scratch array for particle addresses
!     allocate(maskp(npmax,inblok))
! in real space, qe(j+1,k,l) = charge density at grid point (j,kk)
! in real space, qi(j+1,k,l) = ion charge density at grid point (j,kk)
! at grid point (j,kk), where kk = k + noff(l) - 1
      allocate(qe(nxe,nypmx*kbmin,kblok),qi(nxe,nypmx*kbmin,kblok))
! fxyze(i,j+1,k,l) = i component of longitudinal electric force/charge
! exyze(i,j+1,k,l) = i component of transverse electric force/charge
! at grid point (j,kk), where kk = k + noff(l) - 1
      allocate(fxyze(ndim,nxe,nypmx*kbmin,kblok))
      allocate(exyze(ndim,nxe,nypmx*kbmin,kblok))
! bxyze(i,j+1,k,l) = i component of magnetic force at grid point (j,kk)
      allocate(bxyze(2*ndim-3,nxe,nypmx*kbmin,kblok))
! in real space, cu(ndim,j+1,k,l) = i component of current density
! in real space, cus(ndim,j+1,k,l) = i component of acceleration density
! in real space, amu(2*ndim-2,j+1,k,l) = i component of momentum flux
! at grid point (j,kk), where kk = k + noff(l) - 1
      allocate(cu(ndim,nxe,nypmx*kbmin,kblok))
      allocate(cus(ndim,nxe,nypmx*kbmin,kblok))
      allocate(amu(2*ndim-2,nxe,nypmx*kbmin,kblok))
! qt(k,j,l) = complex charge density for fourier mode jj-1,k-1
! cut(k,j,l) = complex current density for fourier mode jj-1,k-1
! cur(k,j,l) = complex acceleration density for fourier mode jj-1,k-1
! amut(k,j,l) = complex momentum flux for fourier mode jj-1,k-1
! where jj = j + kxp*(l - 1)
      allocate(qt(nyv,kxp,jblok),cut(ndim,nyv,kxp,jblok))
      allocate(cur(ndim,nyv,kxp,jblok),amut(2*ndim-2,nyv,kxp,jblok))
! fxyt(i,k,j,l) = i component of force/charge for fourier mode jj-1,k-1
! bxyt(i,k,j,l) = i component of magnetic field for mode jj-1,k-1
! where jj = j + kxp*(l - 1)
      allocate(fxyt(ndim,nyv,kxp,jblok),bxyt(2*ndim-3,nyv,kxp,jblok))
! shift constants
      allocate(q2m(1,1))
! ffc = form factor arrays for poisson solvers
      allocate(ffc(nyh,kxp,jblok),ffe(nyh,kxp,jblok))
! mixup, sct = arrays for fft
      allocate(mixup(nxhy),sct(nxyh))
! edges(1,l) = lower boundary of particle partition l
! edges(2,l) = upper boundary of particle partition l
      allocate(edges(idps,nblok))
! nyp(l) = number of primary gridpoints in particle partition l.
! noff(l) = lowermost global gridpoint in particle partition l.
      allocate(nyp(nblok),noff(nblok))
! npp(l) = number of particles in partition l
! nps(l) = starting address of particles in partition l
      allocate(npp(nblok),nps(nblok))
! initialize parallel timer
      call pwtimer(time,dtime,-1)
! initialize constants
      qbme = qme
      affp = dble(nx)*dble(ny)/np
      qbmi = zero
      omt = sqrt(omx*omx + omy*omy + omz*omz)
      q2m0 = qbme*qme*np/(dble(nx)*dble(ny))
      if (movion==1) then
         qbmi = qmi/rmass
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
         vtdzi = vtz/sqrt(rmass*rtempdzi)
         q2m0 = q2m0 + qbmi*qmi*npi/(dble(nx)*dble(ny))
      endif
      wp0 = q2m0*affp
      q2m = cmplx(q2m0,wp0)
! set initial time
      t0 = dt*real(itime0)
! determine number format and default precisions
      indian = NDIAN()
      rlprec = NDPREC()
      inprec = IDPREC()
! calculate partition variables
      call dcomp(edges,nyp,noff,ny,kstrt,nvp,inorder)
! initialize external magnetic field
      bxyze = 0.0
      if (omt > 0) then
         call baddext(bxyze,nyp,omx,omy,omz,nx,inorder)
         call pcguard(bxyze,kstrt,nvp,kyp,inorder)
         call cguard(bxyze,nyp,nx,inorder)
      endif
      fxyze = 0.0; cu = 0.0
! prepare fft tables
      call fft_init(mixup,sct,indx,indy)
! calculate form factors
      call pois_init(ffc,ax,ay,affp,nx,ny,kstrt)
      call epois_init(ffe,ax,ay,affp,wp0,ci,nx,ny,kstrt)
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npimax,nblok),nppi(nblok))
         allocate(parti2(0,0,0))
      endif
! new start
      if (nustrt==1) then
! initialize electrons
         nps = 1
         npp = 0
! background electrons
         if (npxy > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,   &
     &shiftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            call vdistr(part,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,   &
     &kstrt,nvp,ndim)
!           call vfdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,  &
!    &shiftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
!           call vvdistr(part,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,  &
!    &kstrt,nvp,ndim)
         endif
! beam electrons
         nps = npp + 1
         if (npxyb > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,   &
     &shiftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            call vdistr(part,npp,nps,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,   &
     &npyb,kstrt,nvp,ndim)
!           call vfdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,  &
!    &shiftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
!           call vvdistr(part,npp,nps,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,  &
!    &npyb,kstrt,nvp,ndim)
         endif
! initialize background charge density
         if (movion==0) then
            qi0 = -qme/affp
            call dpostg(part,qi,npp,noff,nyp,kstrt,nvp,nx,kyp,ngds,-qme,&
     &tdpost,inorder,dopt)
! debug
!           call sguard(qi,nyp,qi0,nx,inorder)
         endif
! fix guiding centers for electrons
         if (omt > 0) then
            if (relativity==1) then
               call distr(part,bxyze,npp,noff,qbme,ci,nx,ny,ipbc,inorder&
     &)
            else
               call distr(part,bxyze,npp,noff,qbme,nx,ny,ipbc,inorder)
            endif
         endif
! move electrons into appropriate spatial regions
         call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
! calculate initial electron momentum
         if (ntm > 0) call initmomt2(part,npp,pxe,pye,pze,ndim)
! initialize ions
         if (movion==1) then
            nps = 1
            nppi = 0
! background ions
            if (npxyi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,   &
     &scaledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
               call vdistr(parti,nppi,nps,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0,&
     &npxi,npyi,kstrt,nvp,ndim)
!              call vfdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,  &
!    &scaledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
!              call vvdistr(parti,nppi,nps,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0&
!    &,npxi,npyi,kstrt,nvp,ndim)
            endif
! beam ions
            nps = nppi + 1
            if (npxybi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,   &
     &scaledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi&
     &)
               call vdistr(parti,nppi,nps,vtdxi,vtdyi,vtdzi,vdxi,vdyi,  &
     &vdzi,npxbi,npybi,kstrt,nvp,ndim)
!              call vfdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,  &
!    &scaledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi&
!    &)
!              call vvdistr(parti,nppi,nps,vtdxi,vtdyi,vtdzi,vdxi,vdyi, &
!    &vdzi,npxbi,npybi,kstrt,nvp,ndim)
            endif
! fix guiding centers for ions
            if (omt > 0) then
               if (relativity==1) then
                  call distr(parti,bxyze,nppi,noff,qbmi,ci,nx,ny,ipbc,  &
     &inorder)
               else
                  call distr(parti,bxyze,nppi,noff,qbmi,nx,ny,ipbc,     &
     &inorder)
               endif
            endif
! move ions into appropriate spatial regions
            call pmove(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,  &
     &ierr)
! calculate initial ion momentum
            if (ntm > 0) call initmomt2(parti,nppi,pxi,pyi,pzi,ndim)
         endif
! freeze the ions now
         if ((movion==1).and.(ntime==ionoff)) then
! initialize ion charge density
            call dpostg(parti,qi,nppi,noff,nyp,kstrt,nvp,nx,kyp,ngds,qmi&
     &,tdposti,inorder,dopt)
! delete ions
            deallocate(parti,nppi)
            movion = 0
         endif
! restart
      else
! read restart files
         it = 0
         call restart_bread(iur1,iur2,id0,it,itime,itime0,nvp,npp,part, &
     &movion,nppi,parti,qi,irc,iuer,q2m)
         if (irc /= 0) go to 400
! extend run
         if (nustrt==0) then
            itime0 = itime + itime0
            t0 = dt*real(itime0)
            itime = 0
            ntime = itime + itime0
            if (id0==0) then
               if (iur1 >= 0) close (unit=iur1)
               if (iur2 >= 0) close (unit=iur2)
               call restart_open(1,ntr,idrun,iur1,iur2,iuer)
            endif
            go to 490
         endif
! initiate momentum diagnostic
         if (ntm > 0) then
            call initmomt2(part,npp,pxe,pye,pze,ndim)
            if (movion==1) call initmomt2(parti,nppi,pxi,pyi,pzi,ndim)
         endif
! read diagnostics
         call restart_dread(it,id0,itime,itw,wt,iud,ndrec,fdname,iup,   &
     &nprec,fpname,irc,iuer,iua,narec,faname,iuj,njrec,fjname)
         if (irc /= 0) go to 400
         ntime = itime + itime0
         t0 = dt*real(itime0)
         if (id0==0) rewind it
         go to 490
! handle error
  400    if (id0==0) then
            write (iuer,*) 'Restart Error, irc = ', irc
         endif
         go to 3000
      endif
!
! sorting arrays
  490 allocate(pt(max(npmax,npimax),nblok))
      allocate(ip(max(npmax,npimax),nblok),npic(nypm1,nblok))
      if (sortime > 0) then
         allocate(part2(idimp,npmax,nblok))
      else
         allocate(part2(0,0,0))
      endif
! reduce size of particle manager buffers
!     nbmax = nbmax/2
! initialize diagnostics
! open initial diagnostic metafile
      if (id0==0) then
         iudm = get_funit(iudm)
         fname = 'pdiag2.init.'//cdrun
         open(unit=iudm,file=trim(fname),form='formatted',status=       &
     &'replace')
      endif
! ion density or potential diagnostics
      if ((ntp > 0) .or. (ndp > 0) .or. (ntd > 0) .or. (ndd > 0)) then
         allocate(sfield(nxe,nypmx*kbmin,kblok))
         allocate(sfieldt(nyv,kxp,jblok))
      endif
! ion density diagnostics
      call initmodediag(dent,ntd,id0,nxh,nyh,kxp,modesxd,modesyd,jblok, &
     &iud,ndrec,fdname)
      if (ntd > 0) then
         if (id0==0) then
            ceng = zero
            write (iudm,pden2d,iostat=irc)
         endif
      endif
! velocity diagnostics
      fname = 'fv2.'//cdrun
      call initveldiag(fv,fvm,vtx,vty,vtz,ntv,ndv,id0,nmv,ndim,nblok,iuv&
     &,fname)
      if (movion==1) then
         fname = 'fvi2.'//cdrun
         call initveldiag(fvi,fvmi,vtxi,vtyi,vtzi,ntv,ndv,id0,nmv,ndim, &
     &nblok,iuvi,fname)
      endif
! potential diagnostics
      call initmodediag(pott,ntp,id0,nxh,nyh,kxp,modesxp,modesyp,jblok, &
     &iup,nprec,fpname)
      if (ntp > 0) then
         if (id0==0) then
            ceng = affp
            write (iudm,ppot2d,iostat=irc)
         endif
      endif
! vector potential or ion current diagnostics
      if ((nta > 0) .or. (nda > 0) .or. (ntj > 0) .or. (ndj > 0)) then
         allocate(vfield(ndim,nxe,nypmx*kbmin,kblok))
         allocate(vfieldt(ndim,nyv,kxp,jblok))
         vfield = 0.0
      endif
! vector potential diagnostics
      call initvmodediag(vpott,nta,id0,nxh,nyh,kxp,ndim,modesxa,modesya,&
     &jblok,iua,narec,faname)
      if (nta > 0) then
         if (id0==0) then
            ceng = affp
            write (iudm,pvpot2d,iostat=irc)
         endif
      endif
! ion current diagnostics
      call initvmodediag(vcurt,ntj,id0,nxh,nyh,kxp,ndim,modesxj,modesyj,&
     &jblok,iuj,njrec,fjname)
      if (ntj > 0) then
         if (id0==0) then
            ceng = zero
            write (iudm,pvcur2d,iostat=irc)
         endif
      endif
! momentum diagnostics
      fname = 'pmomentum2.'//cdrun
      if (ntm > 0) then
         if (id0==0) then
            ium = get_funit(ium)
            open(unit=ium,file=trim(fname),form='formatted',status=     &
     &'unknown')
         endif
      endif
! write out and close input file
      if (id0==0) then
         write (iudm,pinput2,iostat=irc)
         close(unit=iudm)
      endif
! record time
      call pwtimer(time,dtime)
! send initial CPU Time to diagnostic nodes
      msg(1) = time(1); msg(2) = time(2)
      call HARTBEAT(msg,2)
      call wtimer(tloop,ldtime,-1)
      if (id0==0) then
         write (iuot,*) 'init max/min real time=',time(1),time(2), 'sec'
      endif
      write (iuer,*) 'Info: Initialization Complete'
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
! send time step to diagnostic nodes
      msg(1) = ntime
      call HARTBEAT(msg,1)
      if (id0==0) write (iuot,991) ntime
      write (label,991) ntime
      call LOGNAME(label)
! initialize current density and momenum flux to background
      call sguard(cu,nyp,zero,zero,zero,nx,inorder)
      call sguard(amu,nyp,zero,zero,zero,zero,nx,inorder)
! deposit current and momentum flux for electrons
      call djpostg(part,cu,npp,noff,qme,zero,ci,tdjpost,nx,ny,ipbc,     &
     &relativity,inorder,djopt)
      call dmjpostg(part,amu,npp,noff,qme,ci,tdcjpost,relativity,inorder&
     &,djopt)
! deposit electron charge density
      call dpostg(part,qe,npp,noff,nyp,kstrt,nvp,nx,kyp,ngds,qme,tdpost,&
     &inorder,dopt)
! save electron current for ion current diagnostic
      if (ndc==0) call vcurdiagprep(cu,vfield,ntj,ndj,ntime)
! deposit current, momentum flux, and charge density for ions
      if (movion==1) then
         call djpostg(parti,cu,nppi,noff,qmi,zero,ci,tdjposti,nx,ny,ipbc&
     &,relativity,inorder,djopt)
         call dmjpostg(parti,amu,nppi,noff,qmi,ci,tdcjposti,relativity, &
     &inorder,djopt)
         call dpostg(parti,qi,nppi,noff,nyp,kstrt,nvp,nx,kyp,ngds,qmi,  &
     &tdposti,inorder,dopt)
! freeze the ions
         if (ntime==ionoff) then
            deallocate(parti,nppi)
            movion = 0
         endif
      endif
! ion current diagnostic
      if (ndc==0) then
         call vcurdiag(cu,cut,vfield,vcurt,vfieldt,ffc,nyp,mixup,sct,   &
     &tfft,ntj,ndj,nx,ny,modesxj,modesyj,iuj,njrec,indx,indy,ntime,nvp, &
     &kstrt,kxp,kyp,ngds,ndstyle,irc,inorder)
         if (irc==1) go to 2000
      endif
! add guard cells for current and momentum flux in x direction
      call aguard(cu,nyp,nx,inorder)
      call amcguard(amu,nyp,nx,inorder)
! add guard cells for current and momentum flux in y direction
      call paguard(cu,kstrt,nvp,nx,kyp,ngds)
      call pamcguard(amu,kstrt,nvp,nx,kyp,ngds)
! add electron and ion densities
      call addqei(qe,qi,qbme,qbmi,wpmax,wpmin,nyp,nx,inorder)
      wp0 = 0.5*(wpmax + wpmin)
! recalculate form factors
      if ((wp0 > 1.15*q2m0) .or. (wp0 < 0.85*q2m0)) then
         q2m0 = wp0
         wp0 = affp*wp0
         q2m = cmplx(q2m0,wp0)
         call epois_init(ffe,ax,ay,affp,wp0,ci,nx,ny,kstrt)
         if (id0==0) then
            write (iuer,*) ntime, 'new shift constants,q2m0,wp0=', q2m
         endif
      endif
! ion density diagnostic
      call dendiag(qt,qi,sfield,dent,sfieldt,ffc,nyp,mixup,sct,tfft,ntd,&
     &ndd,nx,ny,modesxd,modesyd,iud,ndrec,indx,indy,ntime,nvp,kstrt,kxp,&
     &kyp,ndstyle,irc,inorder)
      if (irc==1) go to 2000
! velocity diagnostic
      call veldiag(part,fv,fvm,npp,msg,ntv,ndv,id0,nmv,iuv,ntime,       &
     &' ELECTRON',irc)
      if (irc==1) go to 2000
      if (movion==1) then
         call veldiag(parti,fvi,fvmi,nppi,msg,ntv,ndv,id0,nmv,iuvi,ntime&
     &,' ION',irc)
         if (irc==1) go to 2000
      endif
! phase space diagnostic
      fname = ' ELECTRON PHASE SPACE'
      call phasediag(part,npp,nts,nds,nx,ny,ntime,fname,irc)
      if (irc==1) go to 2000
      if (movion==1) then
         fname = ' ION PHASE SPACE'
         call phasediag(parti,nppi,nts,nds,nx,ny,ntime,fname,irc)
         if (irc==1) go to 2000
      endif
! transform charge to fourier space
      isign = -1
      call fft(qe,qt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! potential diagnostic
      call potdiag(qt,sfield,pott,sfieldt,ffc,nyp,mixup,sct,tfft,ntp,ndp&
     &,nx,ny,modesxp,modesyp,iup,nprec,indx,indy,ntime,nvp,kstrt,kxp,kyp&
     &,ndstyle,irc,inorder)
      if (irc==1) go to 2000
! calculate longitudinal electric force in fourier space
      call pois3(qt,fxyt,ffc,we,tfield,nx,ny,kstrt)
! transform longitudinal electric force to real space
      isign = 1
      call fft(fxyze,fxyt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,     &
     &inorder)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(fxyze,kstrt,nvp,kyp,inorder)
      call cguard(fxyze,nyp,nx,inorder)
! transform current to fourier space
      isign = -1
      call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! take transverse part of current
      call cuperp(cut,nx,ny,kstrt)
! calculate magnetic field in fourier space
      call bpois(cut,bxyt,ffc,ci,wm,tfield,nx,ny,kstrt)
! calculate the momentum in the darwin field
      if (ntm > 0) then
         it = ntime/ntm
         if (ntime==ntm*it) then
            call poynt(qt,cut,ffc,affp,ci,sx,sy,sz,nx,ny,kstrt)
         endif
      endif
! save current for vector potential diagnostic
      if (ndc==0) call vpotdiagprep(cut,vfieldt,nta,nda,ntime)
! transform magnetic field to real space
      isign = 1
      call fft(bxyze,bxyt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,     &
     &inorder)
      if (omt > 0) call baddext(bxyze,nyp,omx,omy,omz,nx,inorder)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(bxyze,kstrt,nvp,kyp,inorder)
      call cguard(bxyze,nyp,nx,inorder)
! transform momentum flux to fourier space
      isign = -1
      call fftn(amu,amut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,      &
     &inorder)
! take transverse part of time derivative of current
      call dcuperp(cur,amut,nx,ny,kstrt)
! calculate convective part of transverse electric field
      call epois(cur,cut,ffe,affp,ci,wf,tfield,nx,ny,kstrt)
! transform transverse electric field to real space
      isign = 1
      call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
!     call fftn(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(cu,kstrt,nvp,kyp,inorder)
      call cguard(cu,nyp,nx,inorder)
! add longitudinal and transverse fields
!     exyze = cu + fxyze
      call addfields(exyze,cu,fxyze)
!
! inner iteration loop
      do k = 1, ndc
! initialize current, acceleration density and momentum flux
      call sguard(cus,cu,nyp,q2m0,nx,inorder)
      call sguard(cu,nyp,zero,zero,zero,nx,inorder)
      call sguard(amu,nyp,zero,zero,zero,zero,nx,inorder)
! deposit electron current and acceleration density and momentum flux
      call dcjpostg(part,exyze,bxyze,npp,noff,cu,cus,amu,qme,qbme,dt,ci,&
     &tdcjpost,relativity,inorder,djopt)
! save electron current for ion current diagnostic
      if (k==ndc) call vcurdiagprep(cu,vfield,ntj,ndj,ntime)
! deposit ion current and acceleration density and momentum flux
      if (movion==1) then
         call dcjpostg(parti,exyze,bxyze,nppi,noff,cu,cus,amu,qmi,qbmi, &
     &dt,ci,tdcjposti,relativity,inorder,djopt)
      endif
! ion current diagnostic
      if (k==ndc) then
         call vcurdiag(cu,cut,vfield,vcurt,vfieldt,ffc,nyp,mixup,sct,   &
     &tfft,ntj,ndj,nx,ny,modesxj,modesyj,iuj,njrec,indx,indy,ntime,nvp, &
     &kstrt,kxp,kyp,ngds,ndstyle,irc,inorder)
         if (irc==1) go to 2000
      endif
! add guard cells for current, acceleration density, and momentum flux
      call aguard(cu,nyp,nx,inorder)
      call aguard(cus,nyp,nx,inorder)
      call amcguard(amu,nyp,nx,inorder)
      call paguard(cu,kstrt,nvp,nx,kyp,ngds)
      call paguard(cus,kstrt,nvp,nx,kyp,ngds)
      call pamcguard(amu,kstrt,nvp,nx,kyp,ngds)
! transform current to fourier space
      isign = -1
      call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! take transverse part of current
      call cuperp(cut,nx,ny,kstrt)
! calculate magnetic field in fourier space
      call bpois(cut,bxyt,ffc,ci,wm,tfield,nx,ny,kstrt)
! calculate the momentum in the darwin field
      if (ntm > 0) then
         it = ntime/ntm
         if (ntime==ntm*it) then
            call poynt(qt,cut,ffc,affp,ci,sx,sy,sz,nx,ny,kstrt)
         endif
      endif
! save current for vector potential diagnostic
      if (k==ndc) call vpotdiagprep(cut,vfieldt,nta,nda,ntime)
! transform magnetic field to real space
      isign = 1
      call fft(bxyze,bxyt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,     &
     &inorder)
      if (omt > 0) call baddext(bxyze,nyp,omx,omy,omz,nx,inorder)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(bxyze,kstrt,nvp,kyp,inorder)
      call cguard(bxyze,nyp,nx,inorder)
! transform acceleration density and momentum flux to fourier space
      isign = -1
      call fft(cus,cur,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
      call fftn(amu,amut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,      &
     &inorder)
! take transverse part of time derivative of current
      call adcuperp(cur,amut,nx,ny,kstrt)
! calculate convective part of transverse electric field
      call epois(cur,cut,ffe,affp,ci,wf,tfield,nx,ny,kstrt)
! transform transverse electric field to real space
      isign = 1
      call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(cu,kstrt,nvp,kyp,inorder)
      call cguard(cu,nyp,nx,inorder)
! add longitudinal and transverse fields
!     exyze = cu + fxyze
      call addfields(exyze,cu,fxyze)
      enddo
!
! vector potential diagnostic
      call vpotdiag(vfieldt,vfield,vpott,cur,ffc,nyp,mixup,sct,ci,tfft, &
     &nta,nda,nx,ny,modesxa,modesya,iua,narec,indx,indy,ntime,nvp,kstrt,&
     &kxp,kyp,ndstyle,irc,inorder)
      if (irc==1) go to 2000
! external pump
!     if ((itpon > 0).and.(ntime >= itpon)) then
!        etx = (v0*vtx)*w0*cos(w0*dt*(ntime - itpon))
!        fxyze(1,:,:,:) = fxyze(1,:,:,:) + etx
!     endif
! particle push and charge density update
      wke = 0.
! push electrons
      call push3g(part,exyze,bxyze,npp,noff,qbme,dt,dt,ci,wke,tpush,nx, &
     &ny,ipbc,relativity,inorder,popt)
! move electrons into appropriate spatial regions
!     call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
      call pmoves(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
! push ions
      if (movion==1) then
         wki = 0.
         call push3g(parti,exyze,bxyze,nppi,noff,qbmi,dt,dt,ci,wki,     &
     &tpushi,nx,ny,ipbc,relativity,inorder,popt)
         wki = wki*rmass
! move ions into appropriate spatial regions
!        call pmove(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ierr)
         call pmoves(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ierr&
     &)
      endif
! calculate electron momentum
      call emomtdiag(part,npp,msg,pxe,pye,pze,sx,sy,sz,wx,wy,wz,ntm,id0,&
     &ium,ntime,ndim)
! calculate ion and total momentum
      call imomtdiag(parti,qi,exyze,nppi,nyp,msg,dt,rmass,pxi,pyi,pzi,wx&
     &,wy,wz,ntm,movion,id0,ium,nx,ntime,ndim,inorder)
! sort electrons
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
!           call sortp(part,pt,ip,npp,noff,nyp,npic,tsort,inorder)
            call sortp(part,part2,npp,noff,nyp,npic,tsort,inorder)
         endif
      endif
! sort ions
      if ((movion==1) .and. (sortimi > 0)) then
         if (mod(ntime,sortimi)==0) then
            call sortp(parti,pt,ip,nppi,noff,nyp,npic,tsorti,inorder)
!           call sortp(parti,parti2,nppi,noff,nyp,npic,tsorti,inorder)
         endif
      endif
! energy diagnostic
      call dmenergy(wt,wtot,msg,we,wf,wm,wke,wki,ntw,ndw,id0,itw,iuot,  &
     &ntime)
      itime = itime + 1
      ntime = itime + itime0
! restart file
      if (ntr > 0) then
         it = ntime/ntr
         if (ntime==ntr*it) then
            it = iur1 + mod(it-1,2)*(iur2 - iur1)
! write file
            call restart_bwrite(it,id0,itime,itime0,nvp,npp,part,movion,&
     &nppi,parti,qi,q2m)
            call restart_dwrite(it,id0,itime,itw,wt,ndrec,fdname,nprec, &
     &fpname,narec,faname,njrec,fjname)
         endif
      endif
      call wtimer(tloop,ldtime)
      ltime = ltime + tloop
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
! send QUIT message to diagnostic nodes
      msg = -1.
      call HARTBEAT(msg,1)
! energy diagnostic
      if (ndw > 0) then
         ts = t0 + dt*real(ndw)*(itime0-(itime0/ndw)*ndw)
         call displayw(wt,ts,dt*real(ndw),itw,irc)
! check error return code
         if (irc==1) go to 3000
      endif
      call pwtimer(time,dtime)
! send main CPU Time to diagnostic nodes
      msg(1:2) = time; msg(3) = tpush; msg(4) = tdpost; msg(5) = tdjpost
      msg(6) = tdcjpost; msg(7) = tsort; msg(8:9) = tmove
      msg(10:11) = tfft; msg(12) = tfield; msg(13) = ltime
      call HARTBEAT(msg,13)
      if (id0==0) then
         write (iuot,*) 'processor partition used: nvp = ', nvp
         write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
         write (iuot,*) 'main max/min real time=',time(1),time(2), 'sec'
         totpush = tpush + tdpost + tdjpost + tdcjpost
         write (iuot,*) 'electron push time=', tpush, 'sec'
         write (iuot,*) 'electron charge deposit time=', tdpost, 'sec'
         write (iuot,*) 'electron current deposit time = ', tdjpost,    &
     &'sec'
         write (iuot,*) 'electron current derivative deposit time = ',  &
     &tdcjpost, 'sec'
         write (iuot,*) 'total electron push time=', totpush, 'sec'
         write (iuot,*) 'electron sort time=', tsort, 'sec'
         write (iuot,*) 'electron move time=', tmove, 'sec'
         totpush = totpush + tsort + tmove(1)
         write (iuot,*) 'total electron time=', totpush, 'sec'
      endif
      if (movion==1) then
         msg(1) = tpushi; msg(2) = tdposti; msg(3) = tdjposti
         msg(4) = tdcjposti; msg(5) = tsorti; msg(6:7) = tmovi
         call HARTBEAT(msg,7)
         if (id0==0) then
            totpushi = tpushi + tdposti + tdjposti + tdcjposti
            write (iuot,*) 'ion push time=', tpushi, 'sec'
            write (iuot,*) 'ion charge deposit time=', tdposti, 'sec'
            write (iuot,*) 'ion current deposit time = ',tdjposti, 'sec'
            write (iuot,*) 'ion current derivative deposit time = ',    &
     &tdcjposti, 'sec'
            write (iuot,*) 'total ion push time=', totpushi, 'sec'
            write (iuot,*) 'ion sort time=', tsorti
            write (iuot,*) 'ion move time=', tmovi
            totpushi = totpushi + tsorti + tmovi(1)
            write (iuot,*) 'total ion time=', totpushi, 'sec'
         endif
      endif
      if (id0==0) then
         write (iuot,*) 'total fft time=', tfft, 'sec'
         write (iuot,*) 'field time=', tfield
         time(1) = time(1) - (totpush + totpushi + tfft(1) + tfield)
         write (iuot,*) 'other time=', time(1), ltime, 'sec'
! write final diagnostic metafile
         fname = 'pdiag2.'//cdrun
         open(unit=iudm,file=trim(fname),form='formatted',status=       &
     &'replace')
! ion density diagnostics
         if (ntd > 0) then
            ndrec = ndrec - 1
            ceng = zero
            write (iudm,pden2d,iostat=irc)
            if (irc /= 0) write (iuer,*) 'pden2d namelist not written'
         endif
! potential diagnostics
         if (ntp > 0) then
            nprec = nprec - 1
            ceng = affp
            write (iudm,ppot2d,iostat=irc)
            if (irc /= 0) write (iuer,*) 'ppot2d namelist not written'
         endif
! vector potential diagnostics
         if (nta > 0) then
            narec = narec - 1
            ceng = affp
            write (iudm,pvpot2d,iostat=irc)
            if (irc /= 0) write (iuer,*) 'pvpot2d namelist not written'
         endif
! ion current diagnostics
         if (ntj > 0) then
            njrec = njrec - 1
            ceng = zero
            write (iudm,pvcur2d,iostat=irc)
            if (irc /= 0) write (iuer,*) 'pvcur2d namelist not written'
         endif
! write out input file
         write (iudm,pinput2,iostat=irc)
         if (irc /= 0) write (iuer,*) 'pinput2 namelist not written'
! done
         write (iuot,*) '* * * q.e.d. * * *'
         close(unit=iudm)
         close(unit=iuot)
         close(unit=iuer)
      endif
! close graphics device
 3000 call PGRCLOSE
      call MP_END
      call PPEXIT
      stop
!
      contains
!
         subroutine diag2nodes
         implicit none
         integer :: j, jt, nt
! diagnostic nodes have special processing
  991    format (' T = ',i7)
  992    format (' field, kinetic, total energies = ',3e14.7)
  993    format (' electric(l,t), magnetic energies = ',3e14.7)
  994    format (' electron momentum = ',3e14.7)
  995    format (' ion momentum = ',3e14.7)
  996    format (' field momentum = ',3e14.7)
  997    format (' total momentum = ',3e14.7)
! allocate data for restart and/or phase space diagnostic
         if ((nts > 0) .or. (nustrt /= 1) .or. (ntr > 0)) then
            allocate(part(idimp,max(npmax,npimax),nblok))
            allocate(npp(nblok))
            parti => part
            nppi => npp
         endif
         if (movion==0) then
            allocate(qi(nxe,nypmx*kbmin,kblok))
         else
            allocate(qi(0,0,0))
         endif
         if ((nustrt /= 1).or.(ntr > 0)) allocate(q2m(1,1))
         if ((ndp > 0) .or. (ndd > 0)) then
            allocate(sfield(nxe,nypmx*kbmin,kblok))
         endif
         if ((nda > 0) .or. (ndj > 0)) then
            allocate(vfield(ndim,nxe,nypmx*kbmin,kblok))
            vfield = 0.0
         endif
! restart
         if (nustrt /= 1) then
            it = 0
            call restart_bread(iur1,iur2,id0,it,itime,itime0,nvp,npp,   &
     &part,movion,nppi,parti,qi,irc,iuer,q2m)
            if (irc /= 0) go to 30
! extend run
            if (nustrt==0) then
               itime0 = itime + itime0
               t0 = dt*real(itime0)
               itime = 0
               if (id0==0) then
                  if (iur1 >= 0) close (unit=iur1)
                  if (iur2 >= 0) close (unit=iur2)
                  call restart_open(1,ntr,idrun,iur1,iur2,iuer)
               endif
               go to 40
            endif
! read diagnostics
            call restart_dread(it,id0,itime,itw,wt,iud,ndrec,fdname,iup,&
     &nprec,fpname,irc,iuer,iua,narec,faname,iuj,njrec,fjname)
            if (irc /= 0) go to 30
            t0 = dt*real(itime0)
            if (id0==0) rewind it
            go to 40
! handle error'
   30       if (id0==0) then
               write (iuer,*) 'Restart Error, irc = ', irc
            endif
            call PGRCLOSE
            call MP_END
            call PPEXIT
            stop
         endif
! initialize diagnostics
   40    if (id0==0) then
            iudm = get_funit(iudm)
            fname = 'pdiag2.init.'//cdrun
            open(unit=iudm,file=trim(fname),form='formatted',status=    &
     &'replace')
         endif
! ion density diagnostics
         call initmodediag(dent,ntd,id0,nxh,nyh,kxp,modesxd,modesyd,    &
     &jblok,iud,ndrec,fdname)
         if (ntd > 0) then
            if (id0==0) then
               ceng = zero
               write (iudm,pden2d,iostat=irc)
            endif
         endif
! velocity diagnostics
         fname = 'fv2.'//cdrun
         call initveldiag(fv,fvm,vtx,vty,vtz,ntv,ndv,id0,nmv,ndim,nblok,&
     &iuv,fname)
         if (movion==1) then
            fname = 'fvi2.'//cdrun
            call initveldiag(fvi,fvmi,vtxi,vtyi,vtzi,ntv,ndv,id0,nmv,   &
     &ndim,nblok,iuvi,fname)
         endif
! potential diagnostics
         call initmodediag(pott,ntp,id0,nxh,nyh,kxp,modesxp,modesyp,    &
     &jblok,iup,nprec,fpname)
         if (ntp > 0) then
            if (id0==0) then
               ceng = affp
               write (iudm,ppot2d,iostat=irc)
            endif
         endif
! vector potential diagnostics
      call initvmodediag(vpott,nta,id0,nxh,nyh,kxp,ndim,modesxa,modesya,&
     &jblok,iua,narec,faname)
      if (nta > 0) then
         if (id0==0) then
            ceng = affp
            write (iudm,pvpot2d,iostat=irc)
         endif
      endif
! ion current diagnostics
      call initvmodediag(vcurt,ntj,id0,nxh,nyh,kxp,ndim,modesxj,modesyj,&
     &jblok,iuj,njrec,fjname)
      if (ntj > 0) then
         if (id0==0) then
            ceng = zero
            write (iudm,pvcur2d,iostat=irc)
         endif
      endif
! momentum diagnostics
         fname = 'pmomentum2.'//cdrun
         if (ntm > 0) then
            if (id0==0) then
               ium = get_funit(ium)
               open(unit=ium,file=trim(fname),form='formatted',status=  &
     &'unknown')
            endif
         endif
! write out and close input file
         if (id0==0) then
            write (iudm,pinput2,iostat=irc)
            close(unit=iudm)
         endif
! get initial CPU Time
         call HARTBEAT(msg,2)
         time(1) = msg(1); time(2) = msg(2)
         call wtimer(tloop,ldtime,-1)
         if (id0==0) then
            write (iuot,*) 'init max/min real time=', time(1), time(2), &
     &'sec'
         endif
! get time step
   10    call HARTBEAT(msg,1)
         it = msg(1)
! end main interation loop
         if (it < 0) then
! energy diagnostic
            if (ndw > 0) then
               ts = t0 + dt*real(ndw)*(itime0-(itime0/ndw)*ndw)
               call displayw(wt,ts,dt*real(ndw),itw,irc)
! check error return code
               if (irc==1) go to 20
            endif
! get main CPU Time
            call HARTBEAT(msg,13)
            time = msg(1:2); tpush = msg(3); tdpost = msg(4)
            tdjpost = msg(5); tdcjpost = msg(6); tsort = msg(7)
            tmove = msg(8:9); tfft = msg(10:11); tfield = msg(12)
            ltime = msg(13)
            if (id0==0) then
               write (iuot,*) 'processor partition used: nvp = ', nvp
               write (iuot,*) ncpus, ' processors found, ', ntasks+1,   &
     &' used'
               write (iuot,*) 'main max/min real time=',time(1),time(2),&
     &'sec'
               totpush = tpush + tdpost + tdjpost + tdcjpost
               write (iuot,*) 'electron push time=', tpush, 'sec'
               write (iuot,*) 'electron charge deposit time=', tdpost,  &
     &'sec'
               write (iuot,*) 'electron current deposit time = ',       &
     &tdjpost, 'sec'
               write (iuot,*) 'electron current derivative deposit time &
     &= ', tdcjpost, 'sec'
               write (iuot,*) 'total electron push time=',totpush, 'sec'
               write (iuot,*) 'electron sort time=', tsort, 'sec'
               write (iuot,*) 'electron move time=', tmove, 'sec'
               totpush = totpush + tsort + tmove(1)
               write (iuot,*) 'total electron time=', totpush, 'sec'
            endif
            if (movion==1) then
               call HARTBEAT(msg,7)
               tpushi = msg(1); tdposti = msg(2); tdjposti = msg(3)
               tdcjposti = msg(4); tsorti = msg(5); tmovi = msg(6:7)
               if (id0==0) then
                  totpushi = tpushi + tdposti + tdjposti + tdcjposti
                  write (iuot,*) 'ion push time=', tpushi, 'sec'
                  write (iuot,*) 'ion charge deposit time=', tdposti,   &
     &'sec'
                  write (iuot,*) 'ion current deposit time = ', tdjposti&
     &, 'sec'
                  write (iuot,*) 'ion current derivative deposit time = &
     &', tdcjposti, 'sec'
                  write (iuot,*) 'total ion push time=', totpushi, 'sec'
                  write (iuot,*) 'ion sort time=', tsorti
                  write (iuot,*) 'ion move time=', tmovi
                  totpushi = totpushi + tsorti + tmovi(1)
                  write (iuot,*) 'total ion time=', totpushi, 'sec'
               endif
            endif
            if (id0==0) then
               write (iuot,*) 'total fft time=', tfft, 'sec'
               write (iuot,*) 'field time=', tfield
               time(1) = time(1) - (totpush + totpushi + tfft(1) +      &
     &tfield)
               write (iuot,*) 'other time=', time(1), ltime, 'sec'
! write final diagnostic metafile
               fname = 'pdiag2.'//cdrun
               open(unit=iudm,file=trim(fname),form='formatted',status= &
     &'replace')
! ion density diagnostics
               if (ntd > 0) then
                  ndrec = ndrec - 1
                  ceng = zero
                  write (iudm,pden2d,iostat=irc)
               if (irc/=0) write (iuer,*) 'pden2d namelist not written'
               endif
! potential diagnostics
               if (ntp > 0) then
                  nprec = nprec - 1
                  ceng = affp
                  write (iudm,ppot2d,iostat=irc)
               if (irc/=0) write (iuer,*) 'ppot2d namelist not written'
               endif
! vector potential diagnostics
               if (nta > 0) then
                  narec = narec - 1
                  ceng = affp
                  write (iudm,pvpot2d,iostat=irc)
               if (irc/=0) write (iuer,*) 'pvpot2d namelist not written'
               endif
! ion current diagnostics
               if (ntj > 0) then
                  njrec = njrec - 1
                  ceng = zero
                  write (iudm,pvcur2d,iostat=irc)
               if (irc/=0) write (iuer,*) 'pvcur2d namelist not written'
               endif
! write out input file
               write (iudm,pinput2,iostat=irc)
               if (irc/=0) write (iuer,*) 'pinput2 namelist not written'
! done
               write (iuot,*) '* * * q.e.d. * * *'
               close(unit=iudm)
               close(unit=iuot)
               close(unit=iuer)
            endif
   20       call PGRCLOSE
            call MP_END
            call PPEXIT
            stop
         else
            ntime = it
         endif
! start main iteration loop
         if (id0==0) write (iuot,991) ntime
         write (label,991) ntime
         call LOGNAME(label)
! ion current diagnostic
         if (ndc==0) then
            if ((ntj > 0) .or. (ndj > 0)) then
               it = -1; if (ntj > 0) it = ntime - ntj*(ntime/ntj)
               jt = -1; if (ndj > 0) jt = ntime - ndj*(ntime/ndj)
               nt = 2*modesyj - 1
               if (it==0) call writebf(vcurt,modesxj,nt,kxp,iuj,njrec)
! display absolute value of ion current
               if (jt==0) then
                  call displayv(vfield,nvp,' ION CURRENT',ntime,999,1,  &
     &ndstyle,nx,ny,irc,inorder)
                  if (irc==1) go to 10
               endif
            endif
         endif
! ion density diagnostic
         if ((ntd > 0) .or. (ndd > 0)) then
            it = -1; if (ntd > 0) it = ntime - ntd*(ntime/ntd)
            jt = -1; if (ndd > 0) jt = ntime - ndd*(ntime/ndd)
            nt = 2*modesyd - 1
            if (it==0) call writebf(dent,modesxd,nt,kxp,iud,ndrec)
! display ion density
            if (jt==0) then
               call displays(sfield,nvp,' ION DENSITY',ntime,999,2,     &
     &ndstyle,nx,ny,irc,inorder)
               if (irc==1) go to 10
            endif
         endif
! velocity diagnostic
         if ((ntv > 0) .or. (ndv > 0)) then
            it = -1; if (ntv > 0) it = ntime - ntv*(ntime/ntv)
            jt = -1; if (ndv > 0) jt = ntime - ndv*(ntime/ndv)
            nt = size(fvm,2)
            if ((it==0) .or. (jt==0)) then
               nt = size(fvm,2)
! get electron velocity moments
               call HARTBEAT(msg,3*nt)
               do k = 1, nt
               do j = 1, 3
               fvm(j,k,1) = msg(j+3*(k-1))
               enddo
               enddo
! print out electron velocity moments
               if (it==0) then
                  nt = ntime/ntv
                  if (id0==0) then
                     write (iuv,*) nt, fvm(1,:,1), fvm(2,:,1),          &
     &sum(fvm(3,:,1))
                  endif
               endif
! display electron velocity distributions
               if (jt==0) then
                  call displayfv(fv,fvm,' ELECTRON',ntime,nmv,2,ndim,irc&
     &)
                  if (irc==1) go to 10
               endif
               if (movion==1) then
                  nt = size(fvmi,2)
! get ion velocity moments
                  call HARTBEAT(msg,3*nt)
                  do k = 1, nt
                  do j = 1, 3
                  fvmi(j,k,1) = msg(j+3*(k-1))
                  enddo
                  enddo
! print out ion velocity moments
                  if (it==0) then
                     nt = ntime/ntv
                     if (id0==0) then
                        write (iuvi,*) nt, fvmi(1,:,1), fvmi(2,:,1),    &
     &sum(fvmi(3,:,1))
                     endif
                  endif
! display ion velocity distributions
                  if (jt==0) then
                     call displayfv(fvi,fvmi,' ION',ntime,nmv,2,ndim,irc&
     &)
                     if (irc==1) go to 10
                  endif
               endif
            endif
         endif
! phase space diagnostic
         fname = ' ELECTRON PHASE SPACE'
         call phasediag(part,npp,nts,nds,nx,ny,ntime,fname,irc)
         if (irc==1) go to 10
         if (movion==1) then
            fname = ' ION PHASE SPACE'
            call phasediag(parti,nppi,nts,nds,nx,ny,ntime,fname,irc)
            if (irc==1) go to 10
         endif
! potential diagnostic
         if ((ntp > 0) .or. (ndp > 0)) then
            it = -1; if (ntp > 0) it = ntime - ntp*(ntime/ntp)
            jt = -1; if (ndp > 0) jt = ntime - ndp*(ntime/ndp)
            nt = 2*modesyp - 1
            if (it==0) call writebf(pott,modesxp,nt,kxp,iup,nprec)
! display potential
            if (jt==0) then
               call displays(sfield,nvp,' POTENTIAL',ntime,999,0,ndstyle&
     &,nx,ny,irc,inorder)
               if (irc==1) go to 10
            endif
         endif
! ion current diagnostic
         if (ndc > 0) then
            if ((ntj > 0) .or. (ndj > 0)) then
               it = -1; if (ntj > 0) it = ntime - ntj*(ntime/ntj)
               jt = -1; if (ndj > 0) jt = ntime - ndj*(ntime/ndj)
               nt = 2*modesyj - 1
               if (it==0) call writebf(vcurt,modesxj,nt,kxp,iuj,njrec)
! display absolute value of ion current
               if (jt==0) then
                  call displayv(vfield,nvp,' ION CURRENT',ntime,999,1,  &
     &ndstyle,nx,ny,irc,inorder)
                  if (irc==1) go to 10
               endif
            endif
         endif
! vector potential diagnostic
         if ((nta > 0) .or. (nda > 0)) then
            it = -1; if (nta > 0) it = ntime - nta*(ntime/nta)
            jt = -1; if (nda > 0) jt = ntime - nda*(ntime/nda)
            nt = 2*modesya - 1
            if (it==0) call writebf(vpott,modesxa,nt,kxp,iua,narec)
! display absolute value of vector potential
            if (jt==0) then
               call displayv(vfield,nvp,' VECTOR POTENTIAL',ntime,999,1,&
     &ndstyle,nx,ny,irc,inorder)
               if (irc==1) go to 10
            endif
         endif
! calculate momentum
         if (ntm > 0) then
            it = ntime/ntm
            it = ntime - ntm*it + 1
            if (it > 1) it = it - ntm
            if (it==1) then
               if (id0==0) write (ium,991) ntime
! get electron momentum values
               call HARTBEAT(msg,6)
               wtot(1:6) = msg(1:6)
! print electron and field momentum
               if (id0==0) then
                  write (ium,994) wtot(1), wtot(2), wtot(3)
                  write (ium,996) wtot(4), wtot(5), wtot(6)
               endif
! get ion momentum values
               call HARTBEAT(msg,6)
               wtot(1:6) = msg(1:6)
! print ion and total momentum
               if (id0==0) then
                  write (ium,995) wtot(1), wtot(2), wtot(3)
                  write (ium,997) wtot(4), wtot(5), wtot(6)
               endif
            endif
         endif
! energy diagnostic
         if ((ntw > 0) .or. (ndw > 0)) then
            it = -1; if (ntw > 0) it = ntime - ntw*((ntime+1)/ntw) + 1
            jt = -1; if (ndw > 0) jt = ntime - ndw*((ntime+1)/ndw) + 1
            if ((it==0) .or. (jt==0)) then
! get energy values
               call HARTBEAT(msg,7)
               wtot = msg(1:7)
               if (it==0) then
                  if (id0==0) then
                     write (iuot,992) wtot(1), wtot(2), wtot(4)
                     write (iuot,993) wtot(5), wtot(6), wtot(7)
                  endif
               endif
               if (jt==0) then
                  itw = itw + 1
                  wt(itw,:) = wtot
               endif
            endif
         endif
         itime = itime + 1
         ntime = itime + itime0
! restart file
         if (ntr > 0) then
            it = ntime/ntr
            if (ntime==ntr*it) then
               it = iur1 + mod(it-1,2)*(iur2 - iur1)
! write file
               call restart_bwrite(it,id0,itime,itime0,nvp,npp,part,    &
     &movion,nppi,parti,qi,q2m)
               call restart_dwrite(it,id0,itime,itw,wt,ndrec,fdname,    &
     &nprec,fpname,narec,faname,njrec,fjname)
            endif
         endif
         call wtimer(tloop,ldtime)
         ltime = ltime + tloop
         go to 10
         end subroutine
!
      end program pdbeps2
