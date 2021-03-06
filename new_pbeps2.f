!-----------------------------------------------------------------------
! * * * periodic 2d electrostatic particle simulation kernel code * * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge, advancing particles, and solving the
! field.  the code moves electrons and ions, with periodic electrostatic
! forces obtained by solving poisson's equation with fast fourier
! transforms.
! portable gcpic kernel code, using algorithm described in:
! p. c. liewer and v. k. decyk, j. computational phys. 85, 302 (1989).
! written by viktor k. decyk, ucla
! for mpi distributed memory computers
! copyright 2000, regents of the university of california
! update: may 1, 2013
      program pbeps2
      use pinit2d
      use pespush2d
      use pfield2d
      use pdiag2d
      use psimul2d
      use mp0d, only: mpinit, ncpus
      implicit none
! idps = number of partition boundaries
! idimp = dimension of phase space = 4
! mshare = (0,1) = (no,yes) architecture is shared memory
      integer :: idps =    2, idimp =   4, mshare =   0
! nmv = number of segments in v for velocity distribution
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      integer :: nmv = 40, vect = 0, ipbc = 1
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10
      integer :: iuvi = 20, iup = 11, ium = 21, iuer = 2
!     integer :: npxy, npxyb, np, npxyi, npxybi, npi
      double precision :: npxy, npxyb, np, npxyi, npxybi, npi
      integer :: nx, ny, nxh, nyh, nyv, nxe, nxeh
      integer :: nloop, nvp, nblok, npmax, npimax = 0, kyp, kxp, nypmx
      integer :: kyb, kxb, kxyb, kbmin, kblok, jbmin, jblok
      integer :: ngds, nxyh, nxhy, nx1, nypm1, nbmax
      integer :: idproc, id0, kstrt, itime, itime0, ntime, isign
      integer :: it, itw, ierr, iur1, iur2, irc
      integer :: ntasks
      real :: zero = 0.0, ltime = 0.0, tloop = 0.0, ts = 0.0
      real :: we = 0.0, wke = 0.0, wki = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tsorti = 0.0
      real :: tfield = 0.0, totpush = 0.0, totpushi = 0.0
      real :: qbme, qbmi, affp, qi0, omt, etx
      real :: sx = 0.0, sy = 0.0, sz = 0.0
      real :: pxe = 0.0, pye = 0.0, pze = 0.0
      real :: pxi = 0.0, pyi = 0.0, pzi = 0.0
      real :: wx, wy, wz
      real :: vtxi, vtyi, vtdxi, vtdyi
      double precision :: dtime, ldtime
      real, dimension(:,:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:,:), pointer :: qe, qi
      real, dimension(:,:,:,:), pointer :: fxye, bxyze
      complex, dimension(:,:,:), pointer :: qt
      complex, dimension(:,:,:,:), pointer :: fxyt
      complex, dimension(:,:,:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:,:), pointer  :: edges
      integer, dimension(:), pointer :: nyp, noff
      integer, dimension(:), pointer :: npp, nppi, nps
      real, dimension(:,:), pointer :: pt
      integer, dimension(:,:), pointer :: ip, npic
      real, dimension(:,:,:), pointer :: sfield
      complex, dimension(:,:,:), pointer :: sfieldt, dent, pott
      real, dimension(:,:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
! wtot = total energy
      real, dimension(4) :: wtot
! time = timing array
      real, dimension(2) :: tfft = 0.0, time = 0.0
      real, dimension(2) :: tmove = 0.0, tmovi = 0.0
! msg = heartbeat array
      double precision, dimension(11) :: msg
      character(len=10) :: cdrun
      character(len=32) :: fname
      character(len=12) :: label
      integer, external :: NDIAN, NDPREC, IDPREC
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
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
      idcode = 1
      psolve = 1
      ndim = 2
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
!     call MP_SETSTACK(262144)
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
!     idimp = 2 + ndim
      if (movion==1) then
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
      endif
! initialize time constants
      itime0 = 0
      itime = itime0
      ntime = itime + itime0
! diagnostic information needed by diagnostic nodes
! set default diagnostic file names
      if (ntd > 0) fdname = 'pdenk2.'//cdrun
      if (ntp > 0) fpname = 'ppotk2.'//cdrun
! energy time history
      if (ndw > 0) then
         allocate(wt((nloop-1)/ndw-(itime0/ndw)+1,4))
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
      allocate(part(idimp,npmax,nblok))
! maskp = scratch array for particle addresses
!     allocate(maskp(npmax,inblok))
! in real space, qe(j+1,k,l) = charge density at grid point (j,kk)
! in real space, qi(j+1,k,l) = ion charge density at grid point (j,kk)
! at grid point (j,kk), where kk = k + noff(l) - 1
      allocate(qe(nxe,nypmx*kbmin,kblok),qi(nxe,nypmx*kbmin,kblok))
! in real space, fxye(i,j+1,k,l) = i component of force/charge at 
! grid point (j,kk)
! in other words, fxye are the convolutions of the electric field
! over the particle shape, where kk = k + noff(l) - 1
      allocate(fxye(ndim,nxe,nypmx*kbmin,kblok))
! qt(k,j,l) = complex charge density for fourier mode jj-1,k-1
! fxyt(1,k,j,l) = x component of force/charge for fourier mode jj-1,k-1
! fxyt(2,k,j,l) = y component of force/charge for fourier mode jj-1,k-1
! where jj = j + kxp*(l - 1)
      allocate(qt(nyv,kxp,jblok),fxyt(ndim,nyv,kxp,jblok))
! ffc = form factor array for poisson solver
      allocate(ffc(nyh,kxp,jblok))
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
!     omt = sqrt(omx*omx + omy*omy + omz*omz)
! debug
      omt = abs(omz)
! end debug
      if (movion==1) then
         qbmi = qmi/rmass
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
      endif
! set initial time
      t0 = dt*real(itime0)
! determine number format and default precisions
      indian = NDIAN()
      rlprec = NDPREC()
      inprec = IDPREC()
! calculate partition variables
      call dcomp(edges,nyp,noff,ny,kstrt,nvp,inorder)
! initialize external magnetic field
      if (omt > 0) then
!        allocate(bxyze(3,nxe,nypmx*kbmin,kblok))
         allocate(bxyze(1,nxe,nypmx*kbmin,kblok))
         bxyze = 0.0
         call baddext(bxyze,nyp,omx,omy,omz,nx,inorder)
         call pcguard(bxyze,kstrt,nvp,kyp,inorder)
         call cguard(bxyze,nyp,nx,inorder)
      endif
! prepare fft tables
      call fft_init(mixup,sct,indx,indy)
! calculate form factors
      call pois_init(ffc,ax,ay,affp,nx,ny,kstrt)
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
!        if (npxy > 0) call distr(part,edges,npp,nps,vtx,vty,vx0,vy0,npx&
!    &,npy,nx,ny,ipbc)
         if (npxy > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,   &
     &shiftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            call vdistr(part,npp,nps,vtx,vty,vx0,vy0,npx,npy,kstrt,nvp)
!           call vfdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,  &
!    &shiftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
!           call vvdistr(part,npp,nps,vtx,vty,vx0,vy0,npx,npy,kstrt,nvp)
         endif
! beam electrons
         nps = npp + 1
!        if (npxyb > 0) call distr(part,edges,npp,nps,vtdx,vtdy,vdx,vdy,&
!    &npxb,npyb,nx,ny,ipbc)
         if (npxyb > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,   &
     &shiftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            call vdistr(part,npp,nps,vtdx,vtdy,vdx,vdy,npxb,npyb,kstrt, &
     &nvp)
!           call vfdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,  &
!    &shiftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
!           call vvdistr(part,npp,nps,vtdx,vtdy,vdx,vdy,npxb,npyb,kstrt,&
!    &nvp)
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
!           if (npxyi > 0) call distr(parti,edges,nppi,nps,vtxi,vtyi,   &
!    &vxi0,vyi0,npxi,npyi,nx,ny,ipbc)
            if (npxyi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,   &
     &scaledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
               call vdistr(parti,nppi,nps,vtxi,vtyi,vxi0,vyi0,npxi,npyi,&
     &kstrt,nvp)
!              call vfdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,  &
!    &scaledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
!              call vvdistr(parti,nppi,nps,vtxi,vtyi,vxi0,vyi0,npxi,npyi&
!    &,kstrt,nvp)
            endif
! beam ions
            nps = nppi + 1
!           if (npxybi > 0) call distr(parti,edges,nppi,nps,vtdxi,vtdyi,&
!    &vdxi,vdyi,npxbi,npybi,nx,ny,ipbc)
            if (npxybi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,   &
     &scaledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi&
     &)
               call vdistr(parti,nppi,nps,vtdxi,vtdyi,vdxi,vdyi,npxbi,  &
     &npybi,kstrt,nvp)
!              call vfdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,  &
!    &scaledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi&
!    &)
!              call vvdistr(parti,nppi,nps,vtdxi,vtdyi,vdxi,vdyi,npxbi, &
!    &npybi,kstrt,nvp)
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
     &movion,nppi,parti,qi,irc,iuer)
         if (irc /= 0) go to 400
! initiate momentum diagnostic
         if (ntm > 0) then
            call initmomt2(part,npp,pxe,pye,pze,ndim)
            if (movion==1) call initmomt2(parti,nppi,pxi,pyi,pzi,ndim)
         endif
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
! read diagnostics
         call restart_dread(it,id0,itime,itw,wt,iud,ndrec,fdname,iup,   &
     &nprec,fpname,irc,iuer)
         if (irc /= 0) go to 400
         ntime = itime + itime0
         t0 = dt*real(itime0)
         if (id0==0) rewind it
         go to 490
! handle error'
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
      call initveldiag(fv,fvm,vtx,vty,zero,ntv,ndv,id0,nmv,ndim,nblok,  &
     &iuv,fname)
      if (movion==1) then
         fname = 'fvi2.'//cdrun
         call initveldiag(fvi,fvmi,vtxi,vtyi,zero,ntv,ndv,id0,nmv,ndim, &
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
! deposit electron charge density
      call dpostg(part,qe,npp,noff,nyp,kstrt,nvp,nx,kyp,ngds,qme,tdpost,&
     &inorder,dopt)
! deposit ion charge density
      if (movion==1) then
         call dpostg(parti,qi,nppi,noff,nyp,kstrt,nvp,nx,kyp,ngds,qmi,  &
     &tdposti,inorder,dopt)
! freeze the ions
         if (ntime==ionoff) then
            deallocate(parti,nppi)
            movion = 0
         endif
      endif
! add electron and ion densities
      call addqei(qe,qi,nyp,nx,inorder)
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
! calculate force/charge in fourier space
      call pois(qt,fxyt,ffc,we,tfield,nx,ny,kstrt)
! transform force/charge to real space
      isign = 1
      call fft(fxye,fxyt,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(fxye,kstrt,nvp,kyp,inorder)
      call cguard(fxye,nyp,nx,inorder)
! external pump
!     if ((itpon > 0).and.(ntime >= itpon)) then
!        etx = (v0*vtx)*w0*cos(w0*dt*(ntime - itpon))
!        fxye(1,:,:,:) = fxye(1,:,:,:) + etx
!     endif
! particle push and charge density update
      wke = 0.
! push electrons
      if (omt <= 0) then
         call pushg(part,fxye,npp,noff,qbme,dt,ci,wke,tpush,nx,ny,ipbc, &
     &relativity,inorder,popt)
      else
         call bpushg(part,fxye,bxyze,npp,noff,qbme,dt,ci,wke,tpush,nx,ny&
     &,ipbc,relativity,inorder,popt)
      endif
! move electrons into appropriate spatial regions
!     call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
      call pmoves(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
! push ions
      if (movion==1) then
         wki = 0.
         if (omt <= 0) then
            call pushg(parti,fxye,nppi,noff,qbmi,dt,ci,wki,tpushi,nx,ny,&
     &ipbc,relativity,inorder,popt)
         else
            call bpushg(parti,fxye,bxyze,nppi,noff,qbmi,dt,ci,wki,tpushi&
     &,nx,ny,ipbc,relativity,inorder,popt)
         endif
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
      call imomtdiag(parti,qi,fxye,nppi,nyp,msg,dt,rmass,pxi,pyi,pzi,wx,&
     &wy,wz,ntm,movion,id0,ium,nx,ntime,ndim,inorder)
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
      call esenergy(wt,wtot,msg,we,wke,wki,ntw,ndw,id0,itw,iuot,ntime)
      itime = itime + 1
      ntime = itime + itime0
! restart file
      if (ntr > 0) then
         it = ntime/ntr
         if (ntime==ntr*it) then
            it = iur1 + mod(it-1,2)*(iur2 - iur1)
! write file
            call restart_bwrite(it,id0,itime,itime0,nvp,npp,part,movion,&
     &nppi,parti,qi)
            call restart_dwrite(it,id0,itime,itw,wt,ndrec,fdname,nprec, &
     &fpname)
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
      msg(1:2) = time; msg(3) = tpush; msg(4) = tdpost; msg(5) = tsort
      msg(6:7) = tmove; msg(8:9) = tfft; msg(10) = tfield
      msg(11) = ltime
      call HARTBEAT(msg,11)
      if (id0==0) then
         write (iuot,*) 'processor partition used: nvp = ', nvp
         write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
         write (iuot,*) 'main max/min real time=',time(1),time(2), 'sec'
         totpush = tpush + tdpost
         write (iuot,*) 'electron push time=', tpush, 'sec'
         write (iuot,*) 'electron charge deposit time=', tdpost, 'sec'
         write (iuot,*) 'total electron push time=', totpush, 'sec'
         write (iuot,*) 'electron sort time=', tsort, 'sec'
         write (iuot,*) 'electron move time=', tmove, 'sec'
         totpush = totpush + tsort + tmove(1)
         write (iuot,*) 'total electron time=', totpush, 'sec'
      endif
      if (movion==1) then
         msg(1) = tpushi; msg(2) = tdposti; msg(3) = tsorti
         msg(4:5) = tmovi
         call HARTBEAT(msg,5)
         if (id0==0) then
            totpushi = tpushi + tdposti
            write (iuot,*) 'ion push time=', tpushi, 'sec'
            write (iuot,*) 'ion charge deposit time=', tdposti, 'sec'
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
! write out input file
         write (iudm,pinput2,iostat=irc)
         if (irc /= 0) write (iuot,*) 'pinput2 namelist not written'
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
         integer :: j, k, jt, nt
         real, dimension(6) :: mtot
! diagnostic nodes have special processing
  991    format (' T = ',i7)
  992    format (' field, kinetic, total energies = ',3e14.7)
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
         if ((ndp > 0) .or. (ndd > 0)) then
            allocate(sfield(nxe,nypmx*kbmin,kblok))
         endif
! restart
         if (nustrt /= 1) then
            it = 0
            call restart_bread(iur1,iur2,id0,it,itime,itime0,nvp,npp,   &
     &part,movion,nppi,parti,qi,irc,iuer)
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
     &nprec,fpname,irc,iuer)
            if (irc /= 0) go to 30
            t0 = dt*real(itime0)
            if (id0==0) rewind it
            go to 40
! handle error
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
         call initveldiag(fv,fvm,vtx,vty,zero,ntv,ndv,id0,nmv,ndim,nblok&
     &,iuv,fname)
         if (movion==1) then
            fname = 'fvi2.'//cdrun
            call initveldiag(fvi,fvmi,vtxi,vtyi,zero,ntv,ndv,id0,nmv,   &
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
            call HARTBEAT(msg,11)
            time = msg(1:2); tpush = msg(3); tdpost = msg(4)
            tsort = msg(5); tmove = msg(6:7); tfft = msg(8:9)
            tfield = msg(10); ltime = msg(11)
            if (id0==0) then
               write (iuot,*) 'processor partition used: nvp = ', nvp
               write (iuot,*) ncpus, ' processors found, ', ntasks+1,   &
     &'used'
               write (iuot,*) 'main max/min real time=',time(1),time(2),&
     &'sec'
               totpush = tpush + tdpost
               write (iuot,*) 'electron push time=', tpush, 'sec'
               write (iuot,*) 'electron charge deposit time=', tdpost,  &
     &'sec'
               write (iuot,*) 'total electron push time=',totpush, 'sec'
               write (iuot,*) 'electron sort time=', tsort, 'sec'
               write (iuot,*) 'electron move time=', tmove, 'sec'
               totpush = totpush + tsort + tmove(1)
               write (iuot,*) 'total electron time=', totpush, 'sec'
            endif
            if (movion==1) then
               call HARTBEAT(msg,5)
               tpushi = msg(1); tdposti = msg(2); tsorti = msg(3)
               tmovi = msg(4:5)
               if (id0==0) then
                  totpushi = tpushi + tdposti
                  write (iuot,*) 'ion push time=', tpushi, 'sec'
                  write (iuot,*) 'ion charge deposit time=', tdposti,   &
     &'sec'
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
! calculate momentum
         if (ntm > 0) then
            it = ntime/ntm
            it = ntime - ntm*it + 1
            if (it > 1) it = it - ntm
            if (it==1) then
               if (id0==0) write (ium,991) ntime
! get electron momentum values
               call HARTBEAT(msg,6)
               mtot = msg(1:6)
! print electron and field momentum
               if (id0==0) then
                  write (ium,994) mtot(1), mtot(2), mtot(3)
                  write (ium,996) mtot(4), mtot(5), mtot(6)
               endif
! get ion momentum values
               call HARTBEAT(msg,6)
               mtot = msg(1:6)
! print ion and total momentum
               if (id0==0) then
                  write (ium,995) mtot(1), mtot(2), mtot(3)
                  write (ium,997) mtot(4), mtot(5), mtot(6)
               endif
            endif
         endif
! energy diagnostic
         if ((ntw > 0) .or. (ndw > 0)) then
            it = -1; if (ntw > 0) it = ntime - ntw*((ntime+1)/ntw) + 1
            jt = -1; if (ndw > 0) jt = ntime - ndw*((ntime+1)/ndw) + 1
            if ((it==0) .or. (jt==0)) then
! get energy values
               call HARTBEAT(msg,4)
               wtot = msg(1:4)
               if (it==0) then
                  if (id0==0) write (iuot,992) wtot(1), wtot(2), wtot(4)
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
     &movion,nppi,parti,qi)
               call restart_dwrite(it,id0,itime,itw,wt,ndrec,fdname,    &
     &nprec,fpname)
            endif
         endif
         call wtimer(tloop,ldtime)
         ltime = ltime + tloop
         go to 10
         end subroutine
!
      end program pbeps2
