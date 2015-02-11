!-----------------------------------------------------------------------
! * * * periodic 2d electromagnetic particle simulation kernel code * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves electrons and ions, with periodic
! electromagnetic forces obtained by solving maxwell's equation with
! fast fourier transforms.
! portable gcpic kernel code, using algorithm described in:
! p. c. liewer and v. k. decyk, j. computational phys. 85, 302 (1989).
! written by viktor k. decyk, ucla
! for mpi distributed memory multiprocessing macintosh g4 computers
! copyright 2000, regents of the university of california
! update: november 4, 2009
      program pbbeps2
      use pinit2d
      use pempush2d
      use pfield2d
      use pdiag2d
!     use psimul2d
      use mp0d, only: mpinit, ncpus
      implicit none
! idps = number of partition boundaries
! idimp = dimension of phase space = 5
! mshare = (0,1) = (no,yes) architecture is shared memory
      integer :: idps =    2, idimp =   5, mshare =   0
! nmv = number of segments in v for velocity distribution
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      integer :: nmv = 40, vect = 0, ipbc = 1
      integer :: npxy, npxyb, np, npxyi, npxybi, npi
      integer :: nx, ny, nxh, nyh, nyv, nxe, nxeh
      integer :: nloop, nvp, nblok, npmax, npimax = 0, kyp, kxp, nypmx
      integer :: kyb, kxb, kxyb, kbmin, kblok, jbmin, jblok
      integer :: ngds, nxyh, nxhy, nx1, nypm1, nbmax
      integer :: idproc, id0, kstrt, itime, isign, isc, irc, ierr
      integer :: it, modesy2p
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0, ltime = 0.0, tloop = 0.0
      real :: we = 0.0, wf = 0.0, wm = 0.0, wef = 0.0,  wke = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tdjpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tdjposti = 0.0, tsorti = 0.0
      real :: totpush = 0.0, totpushi = 0.0
      real :: qbme, qbmi, affp, dth, qi0, etx
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      double precision :: dtime
      real, dimension(:,:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:,:), pointer :: qe, qi
      real, dimension(:,:,:,:), pointer :: cu, fxyze, bxyze
      complex, dimension(:,:,:,:), pointer :: exyz, bxyz
      complex, dimension(:,:,:), pointer :: qt
      complex, dimension(:,:,:,:), pointer :: cut, fxyt, bxyt
      complex, dimension(:,:,:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:,:), pointer  :: edges
      integer, dimension(:), pointer :: nyp, noff
      integer, dimension(:), pointer :: npp, nppi, nps
      real, dimension(:,:), pointer :: pt
      integer, dimension(:,:), pointer :: ip, npic
      real, dimension(:,:,:), pointer :: sfield
      complex, dimension(:,:,:), pointer :: sfieldt, pott
      real, dimension(:,:,:,:), pointer :: vfield
      real, dimension(:,:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
      integer, dimension(2) :: ktime
! wtot = total energy
      real, dimension(7) :: wtot
! time = timing array
      real, dimension(2) :: tfft = 0.0, time = 0.0
      real, dimension(2) :: tmove = 0.0, tmovi = 0.0
      real, dimension(2) :: tjmove = 0.0, tjmovi = 0.0
! msg = heartbeat array
      double precision, dimension(12) :: msg
      character(len=10) :: cdrun
      character(len=32) :: fname
      character(len=12) :: label
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
  993 format (' electric(l,t), magnetic energies = ',3e14.7)
! nvp = number of real or virtual processors
! initialize for parallel processing
      call PPINIT(idproc,id0,nvp)
      kstrt = idproc + 1
! read namelist
      if (id0==0) then
         open(unit=8,file='pinput2',form='formatted',status='old')
         read (8,pinput2)
! create string from idrun
         write (cdrun,'(i10)') idrun
         cdrun = adjustl(cdrun)
! text output file
         fname = 'poutput2.'//cdrun
         open(unit=18,file=trim(fname),form='formatted',status='replace'&
     &)
! open initial diagnostic metafile
         fname = 'pdiag2.init.'//cdrun
         open(unit=19,file=trim(fname),form='formatted',status='replace'&
     &)
      endif
! broadcast namelist to other nodes
      call sendnml()
! np = total number of electrons in simulation
      npxy = npx*npy; npxyb = npxb*npyb; np = npxy + npxyb
! npi = total number of ions in simulation
      npxyi = npxi*npyi; npxybi = npxbi*npybi; npi = npxyi + npxybi
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = ny/2
      nyv = ny + 2; nxe = nx + 4
! kyp = number of complex grids in each field partition in y direction
! nypmx = maximum size of particle partition, including guard cells.
      kyp = (ny - 1)/nvp + 1; nypmx = kyp + 3
! ngds = number of guard cells
      ngds = 3*((idps - 1)/2 + 1)
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
      npmax = (np/nvp)*1.05 + 7000
      if (movion==1) npimax = (npi/nvp)*1.05 + 7000
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
      if (movion==1) then
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
      endif
! diagnostic information needed by diagnostic nodes
! velocity diagnostics
      if (ntv > 0) then
         allocate(fv(2*nmv+2,3,nblok),fvm(3,3,nblok))
         if (id0==0) then
            fname = 'fv2.'//cdrun
            open(unit=10,file=trim(fname),form='formatted',status='unkno&
     &wn')
! write captions
            write (10,*) 'it vdx vdy vdz vtx vty vtz'
         endif
         if (movion==1) then
            allocate(fvi(2*nmv+2,3,nblok),fvmi(3,3,nblok))
            if (id0==0) then
               fname = 'fvi2.'//cdrun
               open(unit=20,file=trim(fname),form='formatted',status='un&
     &known')
! write captions
               write (20,*) 'it vdxi vdyi vdzi vtxi vtyi vtzi'
            endif
         endif
      endif
! density or potential diagnostics
      if ((ntp > 0) .or. (ntd > 0)) then
         allocate(sfield(nxe,nypmx*kbmin,kblok))
      endif
! potential diagnostics
      if (ntp > 0) then
         if (modesxp > nxh) modesxp = nxh
         if (modesyp > nyh) modesyp = nyh
         modesy2p = 2*modesyp - 1
         allocate(pott(modesy2p,min(modesxp,kxp),jblok))
         if (id0==0) then
            write (19,ppot2d,iostat=irc)
         endif
      endif
! write out and close input file
      if (id0==0) then
         write (19,pinput2,iostat=irc)
         close(unit=19)
      endif
! vector potential diagnostics
      if (nta > 0) then
         allocate(vfield(3,nxe,nypmx*kbmin,kblok))
      endif
! energy diagnostics
      if (ntw > 0) allocate(wt((nloop-1)/ntw+1,7))
! open restart files
      if (nustrt /= 1) then
         if (id0==0) then
            fname = 'rstrt1.'//cdrun
            open(unit=16,file=trim(fname),form='unformatted',status='old&
     &')
            fname = 'rstrt2.'//cdrun
            open(unit=17,file=trim(fname),form='unformatted',status='old&
     &')
         endif
      else if (ntr > 0) then
         if (id0==0) then
            fname = 'rstrt1.'//cdrun
            open(unit=16,file=trim(fname),form='unformatted',status='unk&
     &nown')
            fname = 'rstrt2.'//cdrun
            open(unit=17,file=trim(fname),form='unformatted',status='unk&
     &nown')
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
! in real space, qe(j+1,k,l) = charge density at grid point (j,kk)
! in real space, cu(3,j+1,k,l) = current density at grid point (j,kk)
      allocate(qe(nxe,nypmx*kbmin,kblok),cu(3,nxe,nypmx*kbmin,kblok))
! in real space, fxyze(i,j+1,k,l) = i component of force/charge at 
! grid point (j,kk)
! in other words, fxyze are the convolutions of the electric field
! over the particle shape, where kk = k + noff(l) - 1
      allocate(fxyze(3,nxe,nypmx*kbmin,kblok))
! bxyze(i,j+1,k,l) = i component of magnetic field at grid point (j,kk)
      allocate(bxyze(3,nxe,nypmx*kbmin,kblok))
! in fourier space, exyz = transverse electric field
! in fourier space, bxyz = magnetic field    
      allocate(exyz(3,nyv,kxp,jblok),bxyz(3,nyv,kxp,jblok))
! qt(k,j,l) = complex charge density for fourier mode jj-1,k-1
! cut(k,j,l) = complex current density for fourier mode jj-1,k-1
      allocate(qt(nyv,kxp,jblok),cut(3,nyv,kxp,jblok))
! fxyt(i,k,j,l) = i component of force/charge for fourier mode jj-1,k-1
! bxyt(i,k,j,l) = i component of magnetic field for mode jj-1,k-1
! where jj = j + kxp*(l - 1)
      allocate(fxyt(3,nyv,kxp,jblok),bxyt(3,nyv,kxp,jblok))
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
! sorting arrays
      allocate(pt(max(npmax,npimax),nblok))
      allocate(ip(max(npmax,npimax),nblok),npic(nypm1,nblok))
      allocate(part2(idimp,npmax,nblok))
! initialize parallel timer
      call pwtimer(time,dtime,-1)
! initialize constants
      itime = 0
      qbme = qme
      affp = float(nx*ny)/float(np)
      dth = .5*dt
      if (movion==1) then
         qbmi = qmi/rmass
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
         vtdzi = vtz/sqrt(rmass*rtempdzi)
      endif
! diagnostics
! velocity diagnostics
      if (ntv > 0) then
         fv(1,:,:) = 8.*max(vtx,vty,vtz)
         if (movion==1) fvi(1,:,:) = 8.*max(vtxi,vtyi,vtzi)
      endif
! density or potential diagnostics
      if ((ntp > 0) .or. (ntd > 0)) then
         allocate(sfieldt(nyv,kxp,jblok))
      endif
! calculate partition variables
      call dcomp(edges,nyp,noff,ny,kstrt,nvp,inorder)
! initialize electromagnetic fields
      bxyze = 0.
      bxyz = cmplx(0.,0.)
      exyz = cmplx(0.,0.)
      call bguard(bxyze,nyp,nx,inorder)
! prepare fft tables
      call fft_init(mixup,sct,indx,indy)
! calculate form factors
      call pois_init(ffc,ax,ay,affp,nx,ny,kstrt)
! allocate background charge density
      if (movion==0) allocate(qi(nxe,nypmx*kbmin,kblok))
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npimax,nblok),nppi(nblok))
         nullify(parti2)
      endif
! new start
      if (nustrt==1) then
! initialize electrons
         nps = 1
         npp = 0
! background electrons
!        if (npxy > 0) call distr(part,edges,npp,nps,vtx,vty,vx0,vy0,npx&
!    &,npy,nx,ny,ipbc)
!        if (npxy > 0) call distr(part,edges,npp,nps,vtx,vty,vtz,vx0,vy0&
!    &,vz0,npx,npy,nx,ny,ipbc)
         if (npxy > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            call vdistr(part,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,kst&
     &rt,nvp)
         endif
! beam electrons
         nps = npp + 1
!        if (npxyb > 0) call distr(part,edges,npp,nps,vtdx,vtdy,vdx,vdy,&
!    &npxb,npyb,nx,ny,ipbc)
!        if (npxyb > 0) call distr(part,edges,npp,nps,vtdx,vtdy,vtdz,vdx&
!    &,vdy,vdz,npxb,npyb,nx,ny,ipbc)
         if (npxyb > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            call vdistr(part,npp,nps,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,npy&
     &b,kstrt,nvp)
         endif
! fix guiding centers
!        if (relativity==1) then
!           call distr(part,bxyze,npp,noff,qbme,ci,nx,ny,ipbc,inorder)
!        else
!           call distr(part,bxyze,npp,noff,qbme,nx,ny,ipbc,inorder)
!        endif
! move electrons into appropriate spatial regions
         call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! initialize ions
         if (movion==1) then
            nps = 1
            nppi = 0
! background ions
!           if (npxyi > 0) call distr(parti,edges,nppi,nps,vtxi,vtyi,vxi&
!    &0,vyi0,npxi,npyi,nx,ny,ipbc)
!           if (npxyi > 0) call distr(parti,edges,nppi,nps,vtxi,vtyi,vtz&
!    &i,vxi0,vyi0,vzi0,npxi,npyi,nx,ny,ipbc)
            if (npxyi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sca&
     &ledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
               call vdistr(parti,nppi,nps,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0,&
     &npxi,npyi,kstrt,nvp)
            endif
! beam ions
            nps = nppi + 1
!           if (npxybi > 0) call distr(parti,edges,nppi,nps,vtdxi,vtdyi,&
!    &vdxi,vdyi,npxbi,npybi,nx,ny,ipbc)
!           if (npxybi > 0) call distr(parti,edges,nppi,nps,vtdxi,vtdyi,&
!    &vtdzi,vdxi,vdyi,vdzi,npxbi,npybi,nx,ny,ipbc)
            if (npxybi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sca&
     &ledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
               call vdistr(parti,nppi,nps,vtdxi,vtdyi,vtdzi,vdxi,vdyi,vd&
     &zi,npxbi,npybi,kstrt,nvp)
            endif
! fix guiding centers
!           if (relativity==1) then
!              call distr(parti,bxyze,nppi,noff,qbmi,ci,nx,ny,ipbc,inord&
!    &er)
!           else
!              call distr(parti,bxyze,nppi,noff,qbmi,nx,ny,ipbc,inorder)
!           endif
! move ions into appropriate spatial regions
            call pmove(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ie&
     &rr)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
         endif
! initialize background charge density
         if (movion==0) then
            qi0 = -qme/affp
            call sguard(qi,nyp,zero,nx,inorder)
            call dpost(part,qi,-qme,npp,noff,tdpost,inorder,dopt)
! debug
!           call sguard(qi,nyp,qi0,nx,inorder)
! freeze the ions now
         else if ((movion==1).and.(itime==ionoff)) then
            allocate(qi(nxe,nypmx*kbmin,kblok))
! initialize ion charge density to zero
            call sguard(qi,nyp,zero,nx,inorder)
! deposit ion charge
            call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! delete ions
            deallocate(parti,nppi)
            movion = 0
         endif
! retard electron velocities to deposit current
         if (relativity==1) then
            call retard(part,npp,dth,ci,nx,ny,ipbc)
         else
            call retard(part,npp,dth,nx,ny,ipbc)
         endif
! move electrons into appropriate spatial regions
         call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! retard ion velocities to deposit current
         if (movion==1) then
            if (relativity==1) then
               call retard(parti,nppi,dth,ci,nx,ny,ipbc)
            else
               call retard(parti,nppi,dth,nx,ny,ipbc)
            endif
! move ions into appropriate spatial regions
            call pmove(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ie&
     &rr)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
         endif
! restart
      else
! determine most recent restart file
         if (id0==0) then
            read (16,iostat=ierr) ktime(1)
            if (ierr /= 0) ktime(1) = -1
            read (17,iostat=ierr) ktime(2)
            if (ierr /= 0) ktime(2) = -1
            if (ktime(1) > ktime(2)) then
               ktime(2) = 16
            else
               ktime(1) = ktime(2)
               ktime(2) = 17
            endif
         endif
         call plbcast(ktime)
         itime = ktime(1)
         if (itime < 0) go to 400
! read restart file
         it = ktime(2)
         call rddata(part,npp,it,ierr)
         if (ierr /= 0) go to 400
         call rddata(exyz,nvp,it,ierr)
         if (ierr /= 0) go to 400
         call rddata(bxyz,nvp,it,ierr)
         if (ierr /= 0) go to 400
         if (movion==1) then
            call rddata(parti,nppi,it,ierr)
            if (ierr /= 0) go to 400
         endif
         if (movion==0) then
            call rddata(qi,nvp,it,ierr)
            if (ierr /= 0) go to 400
         endif
         if (ntw > 0) then
            call rddata(wt,1,it,ierr)
            if (ierr /= 0) go to 400
            call plbcast(wt)
         endif
         if (ntp > 0) then
            if (id0==0) then
               read (it,iostat=ierr) ktime(1)
               if (ierr /= 0) ktime(1) = -1
               irc = 0
               fname = 'ppotk2.'//cdrun
               call writebf(pott,modesxp,modesy2p,kxp,11,irc,trim(fname)&
     &)
            endif
            call plbcast(ktime)
            nprec = ktime(1)
            if (nprec< 0) go to 400
         endif
         if (id0==0) then
            read (it,iostat=ierr) ktime(1)
            if (ierr /= 0) ktime(1) = -1
            rewind it
         endif
         call plbcast(ktime)
         irc = ktime(1)
         if (irc==itime) go to 490
! handle error
  400    if (id0==0) write (18,*) 'Restart Error'
         go to 3000
      endif
! record time
  490 call pwtimer(time,dtime)
! send initial CPU Time to diagnostic nodes
      msg(1) = time(1); msg(2) = time(2)
      call HARTBEAT(msg,2)
      if (id0==0) then
         write (18,*) 'init max/min real time=', time(1), time(2), 'sec'
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= itime) go to 2000
! send time step to diagnostic nodes
      msg(1) = itime
      call HARTBEAT(msg,1)
      if (id0==0) write (18,991) itime
      write (label,991) itime
      call LOGNAME(label)
! initialize current density to background
      call sguard(cu,nyp,zero,zero,zero,nx,inorder)
! deposit current
      if (relativity==1) then
         call rdjpost(part,cu,npp,noff,qme,dth,ci,tdjpost,nx,ny,ipbc,ino&
     &rder,djopt)
      else
         call djpost(part,cu,npp,noff,qme,dth,tdjpost,nx,ny,ipbc,inorder&
     &,djopt)
      endif
! move electrons into appropriate spatial regions
      call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
      if (ierr /= 0) then
         call MP_END
         call PPEXIT
         stop
      endif
! initialize charge density to background
      call sguard(qe,nyp,zero,nx,inorder)
! deposit charge
      call dpost(part,qe,qme,npp,noff,tdpost,inorder,dopt)
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
            sfield = -qe
! add guard cells for density in x
            call aguard(sfield,nyp,nx,inorder)
! add guard cells for density in y
            call paguard(sfield,kstrt,nvp,nx,kyp,ngds)
! transform density to fourier space
            isign = -1
            call fft(sfield,qt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,&
     &inorder)
! calculate smoothing in fourier space
            call spois(qt,sfieldt,ffc,nx,ny,kstrt)
! transform density to real space
            isign = 1
            call fft(sfield,sfieldt,isign,mixup,sct,tfft,indx,indy,kstrt&
     &,kyp,inorder)
! copy to guard cells
            call pcguard(sfield,kstrt,nvp,kyp,inorder)
            call cguard(sfield,nyp,nx,inorder)
! display density
            call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! deposit ion current
      if (movion==1) then
         if (relativity==1) then
            call rdjpost(parti,cu,nppi,noff,qmi,dth,ci,tdjposti,nx,ny,ip&
     &bc,inorder,djopt)
         else
            call djpost(parti,cu,nppi,noff,qmi,dth,tdjposti,nx,ny,ipbc,i&
     &norder,djopt)
         endif
! move ions into appropriate spatial regions
         call pmove(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! add ion density
         call dpost(parti,qe,qmi,nppi,noff,tdposti,inorder,dopt)
      else
         qe = qe + qi
      endif
! add guard cells for current in x direction
      call aguard(cu,nyp,nx,inorder)
! add guard cells for current in y direction
      call paguard(cu,kstrt,nvp,nx,kyp,ngds)
! add guard cells for density in x
      call aguard(qe,nyp,nx,inorder)
! add guard cells for density in y
      call paguard(qe,kstrt,nvp,nx,kyp,ngds)
! freeze the ions
      if ((movion==1).and.(itime==ionoff)) then
         allocate(qi(nxe,nypmx*kbmin,kblok))
! initialize ion charge density to zero
         call sguard(qi,nyp,zero,nx,inorder)
! deposit ion charge
         call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! delete ions
         deallocate(parti,nppi)
         movion = 0
      endif
! velocity diagnostic
      if (ntv > 0) then
         it = itime/ntv
         if (itime==ntv*it) then
! calculate electron distribution function and moments
            call vdist(part,fv,fvm,npp,nmv,3)
            call plsum(fv(:,:,1))
            fv(1,:,:) = 8.*max(vtx,vty,vtz)
! display velocity distributions
            call displayfv(fv,fvm,' ELECTRON',itime,nmv,2,3,irc)
            if (irc==1) go to 2000
! print out velocity moments
            if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)
            if (movion==1) then
! calculate ion distribution function and moments
               call vdist(parti,fvi,fvmi,nppi,nmv,3)
               call plsum(fvi(:,:,1))
               fvi(1,:,:) = 8.*max(vtxi,vtyi,vtzi)
! display velocity distributions
               call displayfv(fvi,fvmi,' ION',itime,nmv,2,3,irc)
               if (irc==1) go to 2000
! print out velocity moments
               if (id0==0) write (20,*) it, fvmi(1,:,:), fvmi(2,:,:)
            endif
         endif
      endif
! phase space diagnostic
      if (nts > 0) then
         it = itime/nts
         if (itime==nts*it) then
            isc = 999
! plot electrons vx versus x
            call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,ny,&
     &3,1,irc)
            if (irc==1) go to 2000
! plot electrons vy versus y
            call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,ny,&
     &4,2,irc)
            if (irc==1) go to 2000
            if (movion==1) then
! plot ions vx versus x
               call grasp(parti,nppi,' ION PHASE SPACE',itime,isc,nx,ny,&
     &3,1,irc)
               if (irc==1) go to 2000
! plot ions vy versus y
               call grasp(parti,nppi,' ION PHASE SPACE',itime,isc,nx,ny,&
     &4,2,irc)
               if (irc==1) go to 2000
            endif
         endif
      endif
! transform charge to fourier space
      isign = -1
      call fft(qe,qt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! potential diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! calculate potential in fourier space
            call pois(qt,sfieldt,ffc,we,nx,ny,kstrt)
! store selected fourier modes
            call gtmodes(sfieldt,pott,nx,ny,modesxp,modesyp,kstrt)
! transform potential to real space
            isign = 1
            call fft(sfield,sfieldt,isign,mixup,sct,tfft,indx,indy,kstrt&
     &,kyp,inorder)
! copy to guard cells
            call pcguard(sfield,kstrt,nvp,kyp,inorder)
            call cguard(sfield,nyp,nx,inorder)
! display potential
            call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
!           g(:,1) = sfield(nxh,:,1)
!           call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorder)
!           if (irc==1) go to 2000
! write diagnostic output
            if (nprec==0) then
               nprec = -1
               fname = 'ppotk2.'//cdrun
               call writebf(pott,modesxp,modesy2p,kxp,11,nprec,trim(fnam&
     &e))
            else
               call writebf(pott,modesxp,modesy2p,kxp,11,nprec)
            endif
         endif
      endif
! transform current to fourier space
      isign = -1
      call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
!     call fftn(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! take transverse part of current
      call cuperp(cut,nx,ny,kstrt)
! calculate magnetic field in fourier space
      if (itime==0) then
         call ibpois(cut,bxyz,ffc,ci,wm,nx,ny,kstrt)
         wf = 0.
      else
         call maxwel(exyz,bxyz,cut,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt)
      endif
! vector potential diagnostic
      if (nta > 0) then
         it = itime/nta
         if (itime==nta*it) then
! calculate vector potential in fourier space
!           call bpois(cut,bxyt,ffc,ci,wm,nx,ny,kstrt)
            call avpot(bxyz,bxyt,nx,ny,kstrt)
! transform vector potential to real space
            isign = 1
            call fft(vfield,bxyt,isign,mixup,sct,tfft,indx,indy,kstrt,ky&
     &p,inorder)
            call pcguard(vfield,kstrt,nvp,kyp,inorder)
            call bguard(vfield,nyp,nx,inorder)
! display vector potential
            call displayv(vfield,nvp,' VECTOR POTENTIAL',itime,999,1,nds&
     &tyle,nx,ny,irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! calculate force/charge in fourier space
      call pois3(qt,fxyt,ffc,we,nx,ny,kstrt)
! add longitudinal and transverse electric fields
      isign = 1
      call emfield(fxyt,exyz,ffc,isign,nx,ny,kstrt)
! copy magnetic field
      isign = -1
      call emfield(bxyt,bxyz,ffc,isign,nx,ny,kstrt)
! transform force/charge to real space
      isign = 1
      call fft(fxyze,fxyt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inord&
     &er)
!     call fftn(fxyze,fxyt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inor&
!    &der)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(fxyze,kstrt,nvp,kyp,inorder)
      call bguard(fxyze,nyp,nx,inorder)
! transform magnetic field to real space
      isign = 1
      call fft(bxyze,bxyt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inord&
     &er)
!     call fftn(bxyze,bxyt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inor&
!    &der)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(bxyze,kstrt,nvp,kyp,inorder)
      call bguard(bxyze,nyp,nx,inorder)
! external pump
      if ((itpon > 0).and.(itime >= itpon)) then
         etx = (v0*vtx)*w0*cos(w0*dt*(itime - itpon))
         fxyze(1,:,:,:) = fxyze(1,:,:,:) + etx
      endif
! particle push and charge density update
      wke = 0.
! push electrons
      if (relativity==1) then
!        call rpush(part,fxyze,npp,noff,qbme,dt,ci,wke,tpush,nx,ny,ipbc,&
!    &inorder,popt)
!        call rpush(part,fxyze,bxyze,npp,noff,qbme,dt,dth,ci,wke,tpush,n&
!    &x,ny,ipbc,inorder,popt)
         call rpush3(part,fxyze,bxyze,npp,noff,qbme,dt,dth,ci,wke,tpush,&
     &nx,ny,ipbc,inorder,popt)
      else
!        call push(part,fxyze,npp,noff,qbme,dt,wke,tpush,nx,ny,ipbc,inor&
!    &der,popt)
!        call push(part,fxyze,bxyze,npp,noff,qbme,dt,dth,wke,tpush,nx,ny&
!    &,ipbc,inorder,popt)
         call push3(part,fxyze,bxyze,npp,noff,qbme,dt,dth,wke,tpush,nx,n&
     &y,ipbc,inorder,popt)
      endif
! move electrons into appropriate spatial regions
      call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
      if (ierr /= 0) then
         call MP_END
         call PPEXIT
         stop
      endif
! push ions
      if (movion==1) then
         wki = 0.
         if (relativity==1) then
!           call rpush(parti,fxyze,nppi,noff,qbmi,dt,ci,wki,tpushi,nx,ny&
!    &,ipbc,inorder,popt)
!           call rpush(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dth,ci,wki,tp&
!    &ushi,nx,ny,ipbc,inorder,popt)
            call rpush3(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dth,ci,wki,t&
     &pushi,nx,ny,ipbc,inorder,popt)
         else
!           call push(parti,fxyze,nppi,noff,qbmi,dt,wki,tpushi,nx,ny,ipb&
!    &c,inorder,popt)
!           call push(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dth,wki,tpushi&
!    &,nx,ny,ipbc,inorder,popt)
            call push3(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dth,wki,tpush&
     &i,nx,ny,ipbc,inorder,popt)
         endif
         wki = wki*rmass
! move ions into appropriate spatial regions
         call pmove(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
      endif
! sort electrons
      if (sortime > 0) then
         if (mod(itime,sortime)==0) then
!           call sortp(part,pt,ip,npp,noff,nyp,npic,tsort,inorder)
            call sortp(part,part2,npp,noff,nyp,npic,tsort,inorder)
         endif
      endif
! sort ions
      if ((movion==1) .and. (sortimi > 0)) then
         if (mod(itime,sortimi)==0) then
            call sortp(parti,pt,ip,nppi,noff,nyp,npic,tsorti,inorder)
!           call sortp(parti,parti2,nppi,noff,nyp,npic,tsorti,inorder)
         endif
      endif
! energy diagnostic
      if (ntw > 0) then
         it = itime/ntw
         if (itime==ntw*it) then
            wef = we + wf + wm
            wtot(1) = wef
            wtot(2) = wke
            wtot(3) = wki
            wtot(4) = wef + wke + wki
            wtot(5) = we
            wtot(6) = wf
            wtot(7) = wm
            call plsum(wtot)
! send energy values to diagnostic node
            msg(1:7) = wtot
            call HARTBEAT(msg,7)
            if (id0==0) then
               write (18,992) wtot(1), wtot(2), wtot(4)
               write (18,993) wtot(5), wtot(6), wtot(7)
            endif
            wt(it+1,:) = wtot
         endif
      endif
      itime = itime + 1
! restart file
      if (ntr > 0) then
         it = itime/ntr
         if (itime==ntr*it) then
            it = 16 + mod(it-1,2)
            if (id0==0) write (it) itime
            call wrdata(part,npp,it)
            call wrdata(exyz,nvp,it)
            call wrdata(bxyz,nvp,it)
            if (movion==1) call wrdata(parti,nppi,it)
            if (movion==0) call wrdata(qi,nvp,it)
            if (ntw > 0) call wrdata(wt,1,it)
            if ((ntp > 0) .and. (id0==0)) write (it) nprec
            if (id0==0) then
               write (it) itime
               end file it
               rewind it
            endif
         endif
      endif
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
! send QUIT message to diagnostic nodes
      msg = -1.
      call HARTBEAT(msg,1)
! energy diagnostic
      if (ntw > 0) then
         it = (itime - 1)/ntw + 1
         call displayw(wt,dt*real(ntw),it,irc)
! check error return code
         if (irc==1) go to 3000
      endif
      call pwtimer(time,dtime)
! send main CPU Time to diagnostic nodes
      msg(1:2) = time; msg(3) = tpush; msg(4) = tdpost; msg(5) = tdjpost
      msg(6) = tsort; msg(7:8) = tmove; msg(9:10) = tjmove
      msg(11:12) = tfft
      call HARTBEAT(msg,12)
      if (id0==0) then
         write (18,*) 'processor partition used: nvp = ', nvp
         write (18,*) ncpus, ' processors found, ', ntasks+1, ' used'
         write (18,*) 'main max/min real time=', time(1), time(2), 'sec'
         totpush = tpush + tdpost + tdjpost
         write (18,*) 'electron push time=', tpush, 'sec'
         write (18,*) 'electron charge deposit time=', tdpost, 'sec'
         write (18,*) 'electron current deposit time = ', tdjpost, 'sec'
         write (18,*) 'total electron push time=', totpush, 'sec'
         write (18,*) 'electron sort time=', tsort, 'sec'
         tmove = tmove + tjmove
         write (18,*) 'electron move time=', tmove, 'sec'
         totpush = totpush + tsort + tmove(1)
         write (18,*) 'total electron time=', totpush, 'sec'
      endif
      if (movion==1) then
         msg(1) = tpushi; msg(2) = tdposti; msg(3) = tsorti
         msg(4:5) = tmovi; msg(6:7) = tjmovi
         call HARTBEAT(msg,7)
         if (id0==0) then
            totpushi = tpushi + tdposti + tdjposti
            write (18,*) 'ion push time=', tpushi, 'sec'
            write (18,*) 'ion charge deposit time=', tdposti, 'sec'
            write (18,*) 'ion current deposit time = ', tdjposti, 'sec'
            write (18,*) 'total ion push time=', totpushi, 'sec'
            write (18,*) 'ion sort time=', tsorti
            tmovi = tmovi + tjmovi
            write (18,*) 'ion move time=', tmovi
            totpushi = totpushi + tsorti + tmovi(1)
            write (18,*) 'total ion time=', totpushi, 'sec'
         endif
      endif
      if (id0==0) then
         write (18,*) 'total fft time=', tfft, 'sec'
         time(1) = time(1) - (totpush + totpushi + tfft(1))
         write (18,*) 'other time=', time(1), 'sec'
! write final diagnostic metafile
         fname = 'pdiag2.'//cdrun
         open(unit=19,file=trim(fname),form='formatted',status='replace'&
     &)
! potential diagnostics
         if (ntp > 0) then
            nprec = nprec - 1
            write (19,ppot2d,iostat=irc)
         endif
         write (19,pinput2,iostat=irc)
         if (irc /= 0) write (18,*) 'pinput2 namelist not written'
! done
         write (18,*) '* * * q.e.d. * * *'
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
! diagnostic nodes have special processing
  991    format (' T = ',i7)
  992    format (' field, kinetic, total energies = ',3e14.7)
  993    format (' electric(l,t), magnetic energies = ',3e14.7)
! allocate data for restart and/or phase space diagnostic
         if ((nts > 0) .or. (nustrt /= 1) .or. (ntr > 0)) then
            allocate(part(idimp,max(npmax,npimax),nblok))
            allocate(npp(nblok))
            if ((nustrt /= 1).or.(ntr > 0)) then
               allocate(exyz(3,nyv,kxp,jblok),bxyz(3,nyv,kxp,jblok))
            endif
         endif
         if (movion==0) allocate(qi(nxe,nypmx*kbmin,kblok))
! restart
         if (nustrt /= 1) then
! determine most recent restart file
            if (id0==0) then
               read (16,iostat=ierr) ktime(1)
               if (ierr /= 0) ktime(1) = -1
               read (17,iostat=ierr) ktime(2)
               if (ierr /= 0) ktime(2) = -1
               if (ktime(1) > ktime(2)) then
                  ktime(2) = 16
               else
                  ktime(1) = ktime(2)
                  ktime(2) = 17
               endif
            endif
            call plbcast(ktime)
            itime = ktime(1)
            if (itime < 0) go to 30
! read restart file
            it = ktime(2)
            call rddata(part,npp,it,ierr)
            if (ierr /= 0) go to 30
            call rddata(exyz,nvp,it,ierr)
            if (ierr /= 0) go to 30
            call rddata(bxyz,nvp,it,ierr)
            if (ierr /= 0) go to 30
            if (movion==1) then
               call rddata(part,npp,it,ierr)
               if (ierr /= 0) go to 30
            endif
            if (movion==0) then
               call rddata(qi,nvp,it,ierr)
               if (ierr /= 0) go to 30
            endif
            if (ntw > 0) then
               call rddata(wt,1,it,ierr)
               if (ierr /= 0) go to 30
               call plbcast(wt)
            endif
            if (ntp > 0) then
               if (id0==0) then
                  read (it,iostat=ierr) ktime(1)
                  if (ierr /= 0) ktime(1) = -1
                  irc = 0
                  fname = 'ppotk2.'//cdrun
                  call writebf(pott,modesxp,modesy2p,kxp,11,irc,trim(fna&
     &me))
               endif
               call plbcast(ktime)
               nprec = ktime(1)
               if (nprec< 0) go to 30
            endif
            if (id0==0) then
               read (it,iostat=ierr) ktime(1)
            if (ierr /= 0) ktime(1) = -1
               rewind it
            endif
            call plbcast(ktime)
            irc = ktime(1)
            if (irc==itime) go to 40
! handle error
   30       if (id0==0) write (18,*) 'Restart Error'
            call PGRCLOSE
            call MP_END
            call PPEXIT
            stop
         endif
! get initial CPU Time
   40    call HARTBEAT(msg,2)
         time(1) = msg(1); time(2) = msg(2)
         if (id0==0) then
         write (18,*) 'init max/min real time=', time(1), time(2), 'sec'
         endif
! get time step
   10    call HARTBEAT(msg,1)
         it = msg(1)
         if (it < 0) then
! energy diagnostic
            if (ntw > 0) then
               it = (itime - 1)/ntw + 1
               call displayw(wt,dt*real(ntw),it,irc)
! check error return code
               if (irc==1) go to 20
            endif
! get main CPU Time
            call HARTBEAT(msg,12)
            time = msg(1:2); tpush = msg(3); tdpost = msg(4)
            tdjpost = msg(5); tsort = msg(6); tmove = msg(7:8)
            tjmove = msg(9:10); tfft = msg(11:12)
            if (id0==0) then
               write (18,*) 'processor partition used: nvp = ', nvp
               write (18,*) ncpus, ' processors found, ', ntasks+1, ' us&
     &ed'
               write (18,*) 'main max/min real time=', time(1), time(2),&
     &'sec'
               totpush = tpush + tdpost + tdjpost
               write (18,*) 'electron push time=', tpush, 'sec'
               write (18,*) 'electron charge deposit time=', tdpost, 'se&
     &c'
               write (18,*) 'electron current deposit time = ', tdjpost,&
     &'sec'
               write (18,*) 'total electron push time=', totpush, 'sec'
               write (18,*) 'electron sort time=', tsort, 'sec'
               tmove = tmove + tjmove
               write (18,*) 'electron move time=', tmove, 'sec'
               totpush = totpush + tsort + tmove(1)
               write (18,*) 'total electron time=', totpush, 'sec'
            endif
            if (movion==1) then
               call HARTBEAT(msg,7)
               tpushi = msg(1); tdposti = msg(2); tsorti = msg(3)
               tmovi = msg(4:5); tjmovi = msg(6:7)
               if (id0==0) then
                  totpushi = tpushi + tdposti + tdjposti
                  write (18,*) 'ion push time=', tpushi, 'sec'
                  write (18,*) 'ion charge deposit time=', tdposti, 'sec&
     &'
                  write (18,*) 'ion current deposit time = ', tdjposti, &
     &'sec'
                  write (18,*) 'total ion push time=', totpushi, 'sec'
                  write (18,*) 'ion sort time=', tsorti
                  tmovi = tmovi + tjmovi
                  write (18,*) 'ion move time=', tmovi
                  totpushi = totpushi + tsorti + tmovi(1)
                  write (18,*) 'total ion time=', totpushi, 'sec'
               endif
            endif
            if (id0==0) then
               write (18,*) 'total fft time=', tfft, 'sec'
               time(1) = time(1) - (totpush + totpushi + tfft(1))
               write (18,*) 'other time=', time(1), 'sec'
! write final diagnostic metafile
               fname = 'pdiag2.'//cdrun
               open(unit=19,file=trim(fname),form='formatted',status='re&
     &place')
! potential diagnostics
               if (ntp > 0) then
                  nprec = nprec - 1
                  write (19,ppot2d,iostat=irc)
               endif
               write (19,pinput2,iostat=irc)
               if (irc /= 0) write (18,*) 'pinput2 namelist not written'
! done
               write (18,*) '* * * q.e.d. * * *'
            endif
   20       call PGRCLOSE
            call MP_END
            call PPEXIT
            stop
         else
            itime = it
         endif
         if (id0==0) write (18,991) itime
         write (label,991) itime
         call LOGNAME(label)
! density diagnostic
         if (ntd > 0) then
            it = itime/ntd
            if (itime==ntd*it) then
! display density
               call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle&
     &,nx,ny,irc,inorder)
               if (irc==1) go to 10
            endif
         endif
! velocity diagnostic
         if (ntv > 0) then
            it = itime/ntv
            if (itime==ntv*it) then
! display velocity distributions
               call displayfv(fv,fvm,' ELECTRON',itime,nmv,2,3,irc)
               if (irc==1) go to 10
! print out velocity moments
               if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)
               if (movion==1) then
! display velocity distributions
                  call displayfv(fvi,fvmi,' ION',itime,nmv,2,3,irc)
                  if (irc==1) go to 10
! print out velocity moments
                  if (id0==0) write (20,*) it, fvmi(1,:,:), fvmi(2,:,:)
              endif
            endif
         endif
! phase space diagnostic
         if (nts > 0) then
            it = itime/nts
            if (itime==nts*it) then
               isc = 999
! plot electrons vx versus x
               call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,&
     &ny,3,1,irc)
               if (irc==1) go to 10
! plot electrons vy versus y
               call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,&
     &ny,4,2,irc)
               if (irc==1) go to 10
               if (movion==1) then
! plot ions vx versus x
                  call grasp(part,npp,' ION PHASE SPACE',itime,isc,nx,ny&
     &,3,1,irc)
                  if (irc==1) go to 10
! plot ions vy versus y
                  call grasp(part,npp,' ION PHASE SPACE',itime,isc,nx,ny&
     &,4,2,irc)
                  if (irc==1) go to 10
               endif
            endif
         endif
! potential diagnostic
         if (ntp > 0) then
            it = itime/ntp
            if (itime==ntp*it) then
! display potential
               call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle&
     &,nx,ny,irc,inorder)
               if (irc==1) go to 10
!              call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorde&
!    &r)
!              if (irc==1) go to 10
! write diagnostic output
               if (nprec==0) then
                  nprec = -1
                  fname = 'ppotk2.'//cdrun
                  call writebf(pott,modesxp,modesy2p,kxp,11,nprec,trim(f&
     &name))
               else
                  call writebf(pott,modesxp,modesy2p,kxp,11,nprec)
               endif
            endif
         endif
! vector potential diagnostic
         if (nta > 0) then
            it = itime/nta
            if (itime==nta*it) then
! display vector potential
               call displayv(vfield,nvp,' VECTOR POTENTIAL',itime,999,1,&
     &ndstyle,nx,ny,irc,inorder)
               if (irc==1) go to 10
            endif
         endif
! energy diagnostic
         if (ntw > 0) then
            it = itime/ntw
            if (itime==ntw*it) then
! get energy values
               call HARTBEAT(msg,7)
               wtot = msg(1:7)
               if (id0==0) then
                  write (18,992) wtot(1), wtot(2), wtot(4)
                  write (18,993) wtot(5), wtot(6), wtot(7)
               endif
            endif
         endif
! restart file
         itime = itime + 1
         if (ntr > 0) then
            it = itime/ntr
            if (itime==ntr*it) then
               it = 16 + mod(it-1,2)
               if (id0==0) write (it) itime
               call wrdata(part,npp,it)
               call wrdata(exyz,nvp,it)
               call wrdata(bxyz,nvp,it)
               if (movion==1) call wrdata(part,npp,it)
               if (movion==0) call wrdata(qi,nvp,it)
               if (ntw > 0) call wrdata(wt,1,it)
               if ((ntp > 0) .and. (id0==0)) write (it) nprec
               if (id0==0) then
                  write (it) itime
                  end file it
                  rewind it
               endif
            endif
         endif
         go to 10
         end subroutine
!
      end program pbbeps2
