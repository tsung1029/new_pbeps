!-----------------------------------------------------------------------
! * * * periodic 2d electrostatic particle simulation kernel code * * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge, advancing particles, and solving the
! field.  the code moves electrons and ions, with periodic electrostatic
! forces obtained by solving poisson's equation with fast fourier
! transforms.  the only diagnostic is particle and field energy.
! portable gcpic kernel code, using algorithm described in:
! p. c. liewer and v. k. decyk, j. computational phys. 85, 302 (1989).
! written by viktor k. decyk, ucla
! for mpi distributed memory computers
! update: november 14, 2009
      program pbeps2
      use pinit2d
      use pespush2d
      use pnpfield2d
      use pdiag2d
!     use psimul2d
      use mp0d, only: mpinit, ncpus
      implicit none
! idps = number of partition boundaries
! idimp = dimension of phase space = 4
! mshare = (0,1) = (no,yes) architecture is shared memory
      integer :: idps =    2, idimp =   4, mshare =   0
! nmv = number of segments in v for velocity distribution
      integer :: nmv = 40, vect = 0
      integer :: npxy, npxyb, np, npxyi, npxybi, npi
      integer :: nx, ny, nxh, nyh, nyv, nxe, nxeh
      integer :: ipbc, nloop, nvp, nblok, npav, npmax, npimax = 0, kyp
      integer :: kxp, nypmx, kyb, kxb, kxyb, kbmin, kblok, jbmin, jblok
      integer :: ngds, nxyh, nxhy, nx1, nypm1, nbmax
      integer :: idproc, id0, kstrt, itime, isign, isc, irc, ierr
      integer :: nterf, nterg, it, modesy2p
      integer :: ntasks
      real :: zero = 0.0, ltime = 0.0, tloop = 0.0, ts = 0.0
      real :: we = 0.0, wke = 0.0, wki = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tsorti = 0.0
      real :: tfield = 0.0, totpush = 0.0, totpushi = 0.0
      real :: trepart = 0.0, tfmove = 0.0
      real :: qbme, qbmi, affp, qi0, etx, anpav = 0.0, pibal = 0.0
      real :: vtxi, vtyi, vtdxi, vtdyi
      double precision :: dtime, etime
      real, dimension(:,:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:,:), pointer :: qe, qi
      real, dimension(:,:,:,:), pointer :: fxye
      complex, dimension(:,:,:), pointer :: qt
      complex, dimension(:,:,:,:), pointer :: fxyt
      complex, dimension(:,:,:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:,:), pointer  :: edges
      integer, dimension(:), pointer :: nyp, noff, nypu, noffu
      integer, dimension(:), pointer :: npp, nppi, nps
      real, dimension(:,:), pointer :: pt
      integer, dimension(:,:), pointer :: ip, npic
      integer, dimension(:), pointer :: nyps, noffs
      real, dimension(:,:,:), pointer :: sfield
      complex, dimension(:,:,:), pointer :: sfieldt, pott
      real, dimension(:,:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
! semi-periodic boundary conditions
      integer :: indx1, nx2, nx2e
      integer :: kxp2, j2blok, kxb2, kxyb1, kbmin1
      real, dimension(:,:,:), pointer :: q1, sfield1
      real, dimension(:,:,:,:), pointer :: fxy1
      complex, dimension(:,:,:,:), pointer :: fxyt1
      complex, dimension(:,:,:), pointer :: qt1, sfieldt1
      complex, dimension(:,:,:), pointer :: ffb
      integer, dimension(:), pointer :: mixup1
      complex, dimension(:), pointer :: sct1
! dirichlet boundary conditions
      real, dimension(:,:,:), pointer :: qd, sfieldd
      real, dimension(:,:,:,:), pointer :: fxyd
! dirichlet or neumann boundary conditions
      integer :: indy1, ny2, kyp2, kyb2, kxyb2
      integer :: kbmin2, k2blok, jbmin2
      real, dimension(:,:,:), pointer :: q2, sfield2
      real, dimension(:,:,:,:), pointer :: fxy2
      complex, dimension(:,:,:,:), pointer :: fxyt2
      complex, dimension(:,:,:), pointer :: qt2, sfieldt2
      complex, dimension(:,:,:), pointer :: ffd
      integer, dimension(:), pointer :: mixup2
      complex, dimension(:), pointer :: sct2
! vacuum boundary conditions
      real, dimension(:,:,:,:), pointer :: ffg
!
      integer, dimension(2) :: ktime
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
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
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
      kyp = (ny - 1)/nvp + 1
! nypm1 = dimension for index and sorting arrays
      if (imbalance >= 0.0) then
         nypm1 = 3*kyp + 1
      else
         nypm1 = kyp + 1
      endif
! nypmx = maximum size of particle partition, including guard cells.
      nypmx = nypm1 + 2
! ngds = number of guard cells
      ngds = 3*((idps - 1)/2 + 1)
      if (inorder==LINEAR) then
         ax = .912871; ay = .912871
         nxe = nx + 2; ; nypmx = nypm1
         ngds = (idps - 1)/2 + 1
      endif
      nxeh = nxe/2
! check if too many processors
      if (nvp > ny) then
         write (2,*) 'Too many processors requested, ny, nvp=', ny, nvp
         call PPEXIT
         stop
      endif
! boundary conditions
      ipbc = psolve
      if ((psolve==VACUUM_2D).or.(psolve==VACUUM_3D)) ipbc = 2
      if (psolve==NEUMANN_2D) ipbc = 2
! initialize for multiprocessing
      ntasks = mpinit()
      if (dopt==VECTOR) vect = 1
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! nblok = number of particle partitions
      nblok = 1 + mshare*(nvp - 1)
! npav = average number of particles per processor
! npmax = maximum number of electrons in each partition
      npav = np/nvp; npmax = npav*1.2 + 7000
      if (movion==1) npimax = (npi/nvp)*1.2 + 7000
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
! dimension for file output
      nx1 = nx + 1
! nbmax = size of buffer for passing particles between processors
      nbmax = 1 + (2*(npxy*vty + npxyb*vtdy) + 1.4*npxyb*abs(vdy))*dt/ny
      if (movion==1) then
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
      endif
! diagnostic information needed by diagnostic nodes
! velocity diagnostics
      if (ntv > 0) then
         allocate(fv(2*nmv+2,2,nblok),fvm(3,2,nblok))
         if (id0==0) then
            fname = 'fv2.'//cdrun
            open(unit=10,file=trim(fname),form='formatted',status='unkno&
     &wn')
! write captions
            write (10,*) 'it vdx vdy vtx vty'
         endif
         if (movion==1) then
            allocate(fvi(2*nmv+2,2,nblok),fvmi(3,2,nblok))
            if (id0==0) then
               fname = 'fvi2.'//cdrun
               open(unit=20,file=trim(fname),form='formatted',status='un&
     &known')
! write captions
             write (20,*) 'it vdxi vdyi vtxi vtyi'
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
! energy diagnostics
      if (ntw > 0) allocate(wt((nloop-1)/ntw+1,4))
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
      allocate(part(idimp,npmax,nblok))
! in real space, qe(j+1,k,l) = charge density at grid point (j,kk)
! in real space, fxye(i,j+1,k,l) = i component of force/charge at 
! grid point (j,kk)
! in other words, fxye are the convolutions of the electric field
! over the particle shape, where kk = k + noff(l) - 1
      allocate(qe(nxe,nypmx*kbmin,kblok),fxye(2,nxe,nypmx*kbmin,kblok))
! qt(k,j,l) = complex charge density for fourier mode jj-1,k-1
! fxyt(1,k,j,l) = x component of force/charge for fourier mode jj-1,k-1
! fxyt(2,k,j,l) = y component of force/charge for fourier mode jj-1,k-1
! where jj = j + kxp*(l - 1)
      allocate(qt(nyv,kxp,jblok),fxyt(2,nyv,kxp,jblok))
! ffc = form factor array for poisson solver
      allocate(ffc(nyh,kxp,jblok))
! mixup, sct = arrays for fft
      allocate(mixup(nxhy),sct(nxyh))
! edges(1,l) = lower boundary of particle partition l
! edges(2,l) = upper boundary of particle partition l
      allocate(edges(idps,nblok))
! nyp(l) = number of primary gridpoints in particle partition l.
! noff(l) = lowermost global gridpoint in particle partition l.
! nypu(l) = number of primary gridpoints in uniform partition l.
! noffu(l) = lowermost global gridpoint in uniform partition l.
      allocate(nyp(nblok),noff(nblok),nypu(nblok),noffu(nblok))
! npp(l) = number of particles in partition l
! nps(l) = starting address of particles in partition l
      allocate(npp(nblok),nps(nblok))
! sorting arrays
      allocate(pt(max(npmax,npimax),nblok))
      allocate(ip(max(npmax,npimax),nblok),npic(nypm1,nblok))
      allocate(part2(idimp,npmax,nblok))
! data for moving field partitions
      allocate(nyps(nblok),noffs(nblok))
! non-periodic boundary conditions
      indx1 = indx+1; indy1 = indy+1;
      nx2 = 2*nx; ny2 = 2*ny; nx2e = 2*nxe
      kxp2 = (nx - 1)/nvp + 1; kxb2 = nx/kxp2
      j2blok = 1 + mshare*(nx/kxp2 - 1)
! dirichlet boundary conditions
      if (psolve==DIRICHLET_2D) then
         allocate(qd(nyv,kxp2+1,j2blok))
         allocate(fxyd(2,nyv,kxp2+1,j2blok))
         allocate(sfieldd(nyv,kxp2+1,j2blok))
!
         kyp2 = (ny2-1)/nvp + 1; kyb2 = ny2/kyp2
         kxyb2 = max(kxb2,kyb2)
         kbmin2 = 1 + (1 - mshare)*(kxyb2/kxb2 - 1)
         k2blok = 1 + mshare*(ny2/kyp2 - 1)
         jbmin2 = 1 + (1 - mshare)*(kxyb2/kyb2 - 1)
         allocate(q2(nx2e,kyp2*kbmin2,k2blok))
         allocate(fxy2(2,nx2e,kyp2*kbmin2,k2blok))
         allocate(sfield2(nx2e,kyp2*kbmin2,k2blok))
         allocate(qt2(ny2,kxp2,j2blok),fxyt2(2,ny2,kxp2,j2blok))
         allocate(sfieldt2(ny2,kxp2,j2blok))
         allocate(ffd(ny,kxp2,j2blok))
         allocate(mixup2(2*nxhy), sct2(2*nxyh))
! semi-periodic boundary conditions
      else if (psolve==DIRICHLET_PERIODIC_2D) then
         kxyb1 = max(kxb2,kyb)
         kbmin1 = 1 + (1 - mshare)*(kxyb1/kxb2 - 1)
         allocate(q1(nx2e,kyp*kbmin1,kblok))
         allocate(fxy1(2,nx2e,kyp*kbmin1,kblok))
         allocate(sfield1(nx2e,kyp*kbmin1,kblok))
         allocate(qt1(nyv,kxp2,j2blok),fxyt1(2,nyv,kxp2,j2blok))
         allocate(sfieldt1(nyv,kxp2,j2blok))
         allocate(ffb(nyh,kxp2,j2blok))
         allocate(mixup1(max(nx,ny)),sct1(max(nx,nyh)))
! vacuum boundary conditions
      else if ((psolve==VACUUM_2D).or.(psolve==VACUUM_3D)) then
         kyp2 = (ny2-1)/nvp + 1; kyb2 = ny2/kyp2
         kxyb2 = max(kxb2,kyb2)
         kbmin2 = 1 + (1 - mshare)*(kxyb2/kxb2 - 1)
         k2blok = 1 + mshare*(ny2/kyp2 - 1)
         jbmin2 = 1 + (1 - mshare)*(kxyb2/kyb2 - 1)
         allocate(q2(nx2e,kyp2*kbmin2,k2blok))
         allocate(fxy2(2,nx2e,kyp2*kbmin2,k2blok))
         allocate(sfield2(nx2e,kyp2*kbmin2,k2blok))
         allocate(qt2(ny2,kxp2,j2blok),fxyt2(2,ny2,kxp2,j2blok))
         allocate(sfieldt2(ny2,kxp2,j2blok))
         allocate(ffg(4,ny+1,kxp2+1,j2blok))
         allocate(mixup2(2*nxhy), sct2(2*nxyh))
! neumann boundary conditions
      else if (psolve==NEUMANN_2D) then
         kyp2 = (ny2-1)/nvp + 1; kyb2 = ny2/kyp2
         kxyb2 = max(kxb2,kyb2)
         kbmin2 = 1 + (1 - mshare)*(kxyb2/kxb2 - 1)
         k2blok = 1 + mshare*(ny2/kyp2 - 1)
         jbmin2 = 1 + (1 - mshare)*(kxyb2/kyb2 - 1)
         allocate(q2(nx2e,kyp2*kbmin2,k2blok))
         allocate(fxy2(2,nx2e,kyp2*kbmin2,k2blok))
         allocate(sfield2(nx2e,kyp2*kbmin2,k2blok))
         allocate(qt2(ny2,kxp2,j2blok),fxyt2(2,ny2,kxp2,j2blok))
         allocate(sfieldt2(ny2,kxp2,j2blok))
         allocate(ffd(ny,kxp2,j2blok))
         allocate(mixup2(2*nxhy), sct2(2*nxyh))
      endif
!
! initialize parallel timer
      call pwtimer(time,dtime,-1)
! initialize constants
      itime = 0
      nterf = 0
      nterg = kyp - 1
      qbme = qme
      if (ipbc==1) then
         affp = float(nx*ny)/float(np)
      else if (ipbc==2) then
         affp = float((nx-2)*(ny-2))/float(np)
      else if (ipbc.eq.3) then
         affp = float((nx-2)*ny)/float(np)
      endif
      if (movion==1) then
         qbmi = qmi/rmass
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
      endif
! diagnostics
! velocity diagnostics
      if (ntv > 0) then
         fv(1,:,:) = 8.*max(vtx,vty)
         if (movion==1) fvi(1,:,:) = 8.*max(vtxi,vtyi)
      endif
! density or potential diagnostics
      if ((ntp > 0) .or. (ntd > 0)) then
         allocate(sfieldt(nyv,kxp,jblok))
      endif
! calculate uniform partition variables
      call dcomp(edges,nypu,noffu,ny,kstrt,nvp,inorder)
! prepare fft tables
      call fft_init(mixup,sct,indx,indy)
! calculate form factors
      call pois_init(ffc,ax,ay,affp,nx,ny,kstrt)
! dirichlet boundary conditions
      if (psolve==DIRICHLET_2D) then
         call fst_init(mixup,sct2,indx,indy)
         call fft_init(mixup2,sct2,indx1,indy1)
         call poisd_init(ffd,ax,ay,affp,nx,ny,kstrt)
! semi-periodic boundary conditions
      else if (psolve==DIRICHLET_PERIODIC_2D) then
         call fft_init(mixup1,sct1,indx1,indy)
         call poism_init(ffb,ax,ay,affp,nx,ny,kstrt)
! vacuum boundary conditions
      else if (psolve==VACUUM_3D) then
         call fft_init(mixup2,sct2,indx1,indy1)
         call poisc3_init(ffg,q2,qt2,mixup2,sct2,ax,affp,indx,indy,kstrt&
     &)
! neumann boundary conditions
      else if (psolve==NEUMANN_2D) then
         call fft_init(mixup2,sct2,indx1,indy1)
         call poisn_init(ffd,ax,ay,affp,nx,ny,kstrt)
      endif
!
! allocate background charge density
      if (movion==0) allocate(qi(nxe,nypmx*kbmin,kblok))
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npimax,nblok),nppi(nblok))
         nullify(parti2)
      endif
! debug
      if (ipbc==2) vdy = 0.
! new start
      if (nustrt==1) then
! initialize electrons
         nps = 1
         npp = 0
! background electrons
!        if (npxy > 0) call distr(part,edges,npp,nps,vtx,vty,vx0,vy0,npx&
!    &,npy,nx,ny,ipbc)
         if (npxy > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            call vdistr(part,npp,nps,vtx,vty,vx0,vy0,npx,npy,kstrt,nvp)
         endif
! beam electrons
         nps = npp + 1
!        if (npxyb > 0) call distr(part,edges,npp,nps,vtdx,vtdy,vdx,vdy,&
!    &npxb,npyb,nx,ny,ipbc)
         if (npxyb > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            call vdistr(part,npp,nps,vtdx,vtdy,vdx,vdy,npxb,npyb,kstrt,n&
     &vp)
         endif
! find new partition analytically
         if (imbalance >= 0.0) then
            call fedges(edges,noff,nyp,ampdy,scaledy,shiftdy,ny,kstrt,nv&
     &p,nypmx,ipbc,ndprof,nterg,ierr,inorder)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
! use uniform partition
         else
            noff = noffu; nyp = nypu
         endif
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
            if (npxyi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sca&
     &ledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
               call vdistr(parti,nppi,nps,vtxi,vtyi,vxi0,vyi0,npxi,npyi,&
     &kstrt,nvp)
            endif
! beam ions
            nps = nppi + 1
!           if (npxybi > 0) call distr(parti,edges,nppi,nps,vtdxi,vtdyi,&
!    &vdxi,vdyi,npxbi,npybi,nx,ny,ipbc)
            if (npxybi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sca&
     &ledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
               call vdistr(parti,nppi,nps,vtdxi,vtdyi,vdxi,vdyi,npxbi,np&
     &ybi,kstrt,nvp)
            endif
! use electron partition for ions
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
            call sguardp(qi,kstrt,nvp,noff,nyp,zero,nx,ny,ipbc,inorder)
            call dpost(part,qi,-qme,npp,noff,tdpost,inorder,dopt)
! debug
!           call sguardp(qi,kstrt,nvp,noff,nyp,qi0,nx,ny,ipbc,inorder)
! freeze the ions now
         else if ((movion==1).and.(itime==ionoff)) then
            allocate(qi(nxe,nypmx*kbmin,kblok))
! initialize ion charge density to zero
            call sguardp(qi,kstrt,nvp,noff,nyp,zero,nx,ny,ipbc,inorder)
! deposit ion charge
            call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! delete ions
            deallocate(parti,nppi)
            movion = 0
         endif
! add guard cells for ion density in x
         if (movion==0) then
            call aguardp(qi,nyp,nx,ipbc,inorder)
! add guard cells for ion density in y
            call pnaguardp(qi,nyp,kstrt,nvp,nx,nterg,ngds,ipbc)
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
         call rddata(edges,nvp,it,ierr)
         if (ierr /= 0) go to 400
         call fnoff(edges,noff,nyp,nypmx,nterg,ierr,inorder)
         if (ierr /= 0) go to 400
         if (movion==1) then
            call rddata(parti,nppi,it,ierr)
            if (ierr /= 0) go to 400
         endif
         if (movion==0) then
            call rddata(qi,nvp,it,ierr)
            if (ierr /= 0) go to 400
            isign = 1
            call pfmove(qi,noff,nyp,isign,tfmove,kyp,kstrt,nvp,idps,nter&
     &f,ierr,inorder)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
            call zguard(qi,nyp,nx,inorder)
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
! initialize charge density to background
      call sguardp(qe,kstrt,nvp,noff,nyp,zero,nx,ny,ipbc,inorder)
! deposit charge
      call dpost(part,qe,qme,npp,noff,tdpost,inorder,dopt)
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
            sfield = -qe
! add guard cells for density in x
            call aguardp(sfield,nyp,nx,ipbc,inorder)
! add guard cells for density in y
            call pnaguardp(sfield,nyp,kstrt,nvp,nx,nterg,ngds,ipbc)
! move density to uniform field partition
            isign = -1
            call pfmove(sfield,noff,nyp,isign,tfmove,kyp,kstrt,nvp,idps,&
     &nterf,ierr,inorder)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
         endif
      endif
! add ion density
      if (movion==1) then
         call dpost(parti,qe,qmi,nppi,noff,tdposti,inorder,dopt)
      else
         qe = qe + qi
      endif
! add guard cells for density in x
      call aguardp(qe,nyp,nx,ipbc,inorder)
! add guard cells for density in y
      call pnaguardp(qe,nyp,kstrt,nvp,nx,nterg,ngds,ipbc)
! freeze the ions
      if ((movion==1).and.(itime==ionoff)) then
         allocate(qi(nxe,nypmx*kbmin,kblok))
! initialize ion charge density to zero
         call sguardp(qi,kstrt,nvp,noff,nyp,zero,nx,ny,ipbc,inorder)
! deposit ion charge
         call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! add guard cells for ion density in x
         call aguardp(qi,nyp,nx,ipbc,inorder)
! add guard cells for ion density in y
         call pnaguardp(qi,nyp,kstrt,nvp,nx,nterg,ngds,ipbc)
! delete ions
         deallocate(parti,nppi)
         movion = 0
      endif
! velocity diagnostic
      if (ntv > 0) then
         it = itime/ntv
         if (itime==ntv*it) then
! calculate electron distribution function and moments
            call vdist(part,fv,fvm,npp,nmv,2)
            call plsum(fv(:,:,1))
            fv(1,:,:) = 8.*max(vtx,vty)
! display velocity distributions
            call displayfv(fv,fvm,' ELECTRON',itime,nmv,2,2,irc)
            if (irc==1) go to 2000
! print out velocity moments
            if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)
            if (movion==1) then
! calculate ion distribution function and moments
               call vdist(parti,fvi,fvmi,nppi,nmv,2)
               call plsum(fvi(:,:,1))
               fvi(1,:,:) = 8.*max(vtxi,vtyi)
! display velocity distributions
               call displayfv(fvi,fvmi,' ION',itime,nmv,2,2,irc)
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
! move charge density to uniform field partition
      isign = -1
      call pfmove(qe,noff,nyp,isign,tfmove,kyp,kstrt,nvp,idps,nterf,ierr&
     &,inorder)
      if (ierr /= 0) then
         call MP_END
         call PPEXIT
         stop
      endif
!
! dirichlet boundary conditions
!
      if (psolve==DIRICHLET_2D) then
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
! transform electron density to fourier space
!           call dblsin(sfield,q2,nx,ny,kstrt,kyp,kyp2,inorder)
            isign = -1
!           call fft(q2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp&
!    &2,LINEAR)
            call fsst(sfield,qd,isign,mixup,sct2,tfft,indx,indy,kstrt,kx&
     &p2,kyp,inorder)
! calculate smoothing in fourier space
!           call poisdx(qt2,sfieldt2,ffd,nx,ny,kstrt)
            call poisd(qd,sfieldd,ffd,nx,ny,kstrt)
! transform electron density to real space
            isign = 1
!           call fft(sfield2,sfieldt2,isign,mixup2,sct2,tfft,indx1,indy1&
!    &,kstrt,kyp2,LINEAR)
!           call hafdbl(sfield,sfield2,nx,ny,kstrt,kyp,kyp2,inorder)
            call fsst(sfield,sfieldd,isign,mixup,sct2,tfft,indx,indy,kst&
     &rt,kxp2,kyp,inorder)
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
! display electron density
            call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
!     call dblsin(qe,q2,nx,ny,kstrt,kyp,kyp2,inorder)
      isign = -1
!     call fft(q2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2,LINE&
!    &AR)
      call fsst(qe,qd,isign,mixup,sct2,tfft,indx,indy,kstrt,kxp2,kyp,ino&
     &rder)
! potential diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! solve for potential
!           call poisdx(qt2,sfieldt2,ffd,we,nx,ny,kstrt)
            call poisd(qd,sfieldd,ffd,nx,ny,kstrt)
            isign = 1
!           call fft(sfield2,sfieldt2,isign,mixup2,sct2,tfft,indx1,indy1&
!    &,kstrt,kyp2,LINEAR)
!           call hafdbl(sfield,sfield2,nx,ny,kstrt,kyp,kyp2,inorder)
            call fsst(sfield,sfieldd,isign,mixup,sct2,tfft,indx,indy,kst&
     &rt,kxp2,kyp,inorder)
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
! display potential
            call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
!           g(:,1) = sfield(nxh,:,1)
!           call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorder)
!           if (irc==1) go to 2000
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,kyp,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
!     call poisdx(qt2,fxyt2,ffd,we,nx,ny,kstrt)
      call poisd(qd,fxyd,ffd,we,nx,ny,kstrt)
      isign = 1
!     call fft(fxy2,fxyt2,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2,LINEAR&
!    &)
!     call hafdbl(fxye,fxy2,nx,ny,kstrt,kyp,kyp2,inorder)
      call fcst(fxye,fxyd,isign,mixup,sct2,tfft,indx,indy,kstrt,kxp2,kyp&
     &,inorder)
!
! semi-periodic boundary conditions
!
      else if (psolve==DIRICHLET_PERIODIC_2D) then
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
! transform electron density to fourier space
            call sglsin(sfield,q1,nx,kyp,inorder)
            isign = -1
            call fft(q1,qt1,isign,mixup1,sct1,tfft,indx1,indy,kstrt,kyp,&
     &LINEAR)
! calculate smoothing in fourier space
            call poism(qt1,sfieldt1,ffb,nx,ny,kstrt)
! transform electron density to real space
            isign = 1
            call fft(sfield1,sfieldt1,isign,mixup1,sct1,tfft,indx1,indy,&
     &kstrt,kyp,LINEAR)
            call hafsgl(sfield,sfield1,nx,kyp,inorder)
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
! display electron density
            call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
      call sglsin(qe,q1,nx,kyp,inorder)
      isign = -1
      call fft(q1,qt1,isign,mixup1,sct1,tfft,indx1,indy,kstrt,kyp,LINEAR&
     &)
! potential diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! solve for potential
            call poism(qt1,sfieldt1,ffb,we,nx,ny,kstrt)
            isign = 1
            call fft(sfield1,sfieldt1,isign,mixup1,sct1,tfft,indx1,indy,&
     &kstrt,kyp,LINEAR)
            call hafsgl(sfield,sfield1,nx,kyp,inorder)
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
! display potential
            call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
!           g(:,1) = sfield(nxh,:,1)
!           call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorder)
!           if (irc==1) go to 2000
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,kyp,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
      call poism(qt1,fxyt1,ffb,we,nx,ny,kstrt)
      call fft(fxy1,fxyt1,mixup1,sct1,tfft,indx1,indy,kstrt,kyp,LINEAR)
      call hafsgl(fxye,fxy1,nx,kyp,inorder)
!
! vacuum boundary conditions
!
      else if (psolve==VACUUM_3D) then
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
! transform electron density to fourier space
            call zdbl(sfield,q2,nx,ny,kstrt,kyp,kyp2,inorder)
            isign = -1
            call fft(q2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp&
     &2,LINEAR)
! calculate smoothing in fourier space
            call poisc(qt2,sfieldt2,ffg,nx,ny,kstrt)
! transform electron density to real space
            isign = 1
            call fft(sfield2,sfieldt2,isign,mixup2,sct2,tfft,indx1,indy1&
     &,kstrt,kyp2,LINEAR)
            call hafdbl(sfield,sfield2,nx,ny,kstrt,kyp,kyp2,inorder)
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
! display electron density
            call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
      call zdbl(qe,q2,nx,ny,kstrt,kyp,kyp2,inorder)
      isign = -1
      call fft(q2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2,LINE&
     &AR)
! potential diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! solve for potential
            call poisc(qt2,sfieldt2,ffg,we,nx,ny,kstrt)
            isign = 1
            call fft(sfield2,sfieldt2,isign,mixup2,sct2,tfft,indx1,indy1&
     &,kstrt,kyp2,LINEAR)
            call hafdbl(sfield,sfield2,nx,ny,kstrt,kyp,kyp2,inorder)
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
! display potential
            call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
!           g(:,1) = sfield(nxh,:,1)
!           call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorder)
!           if (irc==1) go to 2000
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,kyp,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
      call poisc(qt2,fxyt2,ffg,we,nx,ny,kstrt)
      call fft(fxy2,fxyt2,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2,LINEAR&
     &)
      call hafdbl(fxye,fxy2,nx,ny,kstrt,kyp,kyp2,inorder)
!
! neumann boundary conditions
!
      else if (psolve==NEUMANN_2D) then
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
! transform electron density to fourier space
            call dblcos(sfield,q2,nx,ny,kstrt,kyp,kyp2,inorder)
            isign = -1
            call fft(q2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp&
     &2,LINEAR)
! calculate smoothing in fourier space
            call poisn(qt2,sfieldt2,ffd,nx,ny,kstrt)
! transform electron density to real space
            isign = 1
            call fft(sfield2,sfieldt2,isign,mixup2,sct2,tfft,indx1,indy1&
     &,kstrt,kyp2,LINEAR)
            call hafdbl(sfield,sfield2,nx,ny,kstrt,kyp,kyp2,inorder)
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
! display electron density
            call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
      call dblcos(qe,q2,nx,ny,kstrt,kyp,kyp2,inorder)
      isign = -1
      call fft(q2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2,LINE&
     &AR)
! potential diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! solve for potential
            call poisn(qt2,sfieldt2,ffd,we,nx,ny,kstrt)
            isign = 1
            call fft(sfield2,sfieldt2,isign,mixup2,sct2,tfft,indx1,indy1&
     &,kstrt,kyp2,LINEAR)
            call hafdbl(sfield,sfield2,nx,ny,kstrt,kyp,kyp2,inorder)
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
! display potential
            call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
!           g(:,1) = sfield(nxh,:,1)
!           call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorder)
!           if (irc==1) go to 2000
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,kyp,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
      call poisn(qt2,fxyt2,ffd,we,nx,ny,kstrt)
      call fft(fxy2,fxyt2,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2,LINEAR&
     &)
      call hafdbl(fxye,fxy2,nx,ny,kstrt,kyp,kyp2,inorder)
!
! periodic boundary conditions
!
      else if (psolve==PERIODIC_2D) then
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
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
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
! display density
            call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
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
            call pcguardp(sfield,kstrt,nvp,nx,kyp,ipbc,inorder)
            call cguardp(sfield,nypu,nx,ipbc,inorder)
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
! calculate force/charge in fourier space
      call pois(qt,fxyt,ffc,we,tfield,nx,ny,kstrt)
! transform force/charge to real space
      isign = 1
      call fft(fxye,fxyt,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
!
      endif
!
! move force/charge to non-uniform field partition
      isign = 1
      call pfmove(fxye,noff,nyp,isign,tfmove,kyp,kstrt,nvp,idps,nterf,ie&
     &rr,inorder)
      if (ierr /= 0) then
         call MP_END
         call PPEXIT
         stop
      endif
! copy to guard cells
      call pncguardp(fxye,nyp,kstrt,nvp,nx,nterg,ipbc,inorder)
      call cguardp(fxye,nyp,nx,ipbc,inorder)
! external pump
      if ((itpon > 0).and.(itime >= itpon)) then
         etx = (v0*vtx)*w0*cos(w0*dt*(itime - itpon))
         fxye(1,:,:,:) = fxye(1,:,:,:) + etx
      endif
! particle push and charge density update
      wke = 0.
! push electrons
      if (relativity==1) then
         call rpush(part,fxye,npp,noff,qbme,dt,ci,wke,tpush,nx,ny,ipbc,i&
     &norder,popt)
      else
         call push(part,fxye,npp,noff,qbme,dt,wke,tpush,nx,ny,ipbc,inord&
     &er,popt)
      endif
! move electrons into appropriate spatial regions
      call pmove(part,edges,npp,anpav,pibal,tmove,ny,kstrt,nvp,nbmax,vec&
     &t,ierr)
      if (ierr /= 0) then
         call MP_END
         call PPEXIT
         stop
      endif
! push ions
      if (movion==1) then
         wki = 0.
         if (relativity==1) then
            call rpush(parti,fxye,nppi,noff,qbmi,dt,ci,wki,tpushi,nx,ny,&
     &ipbc,inorder,popt)
         else
            call push(parti,fxye,nppi,noff,qbmi,dt,wki,tpushi,nx,ny,ipbc&
     &,inorder,popt)
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
!
! begin repartitioning
      if (imbalance >= 0.0) then
      if (pibal > imbalance) then
! initialize timer
         call wtimer(ts,etime,-1)
! count the number of electrons per cell
         call countp(part,npic,npp,noff,nyp)
! save old repartitioning boundaries
         noffs = noff
         nyps = nyp
! determine new repartitioning boundaries
         call repart(edges,npic,noff,nyp,anpav,kstrt,nvp,nypmx,nterg,ier&
     &r,inorder)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         nterf = 0
         if (movion==0) then
! move background ion density to new field partition.
            call pfmove(qi,noff,nyp,noffs,nyps,tfmove,kstrt,nvp,idps,ier&
     &r,inorder)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
! zero out guard cells
            call zguard(qi,nyp,nx,inorder)
         endif
! move electrons into new spatial regions
         call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! use electron partition for ions
         if (movion==1) then
! move ions into new spatial regions
            call pmove(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ie&
     &rr)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
         endif
! record time
         call wtimer(ts,etime)
         trepart = trepart + ts
         if (id0==0) then
            write (18,*) 'repartitioning complete, imbalance = ', pibal
         endif
      endif
      endif
! end repartitioning
!
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
            wtot(1) = we
            wtot(2) = wke
            wtot(3) = wki
            wtot(4) = we + wke + wki
            call plsum(wtot)
! send energy values to diagnostic node
            msg(1:4) = wtot
            call HARTBEAT(msg,4)
            if (id0==0) write (18,992) wtot(1), wtot(2), wtot(4)
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
            call wrdata(edges,nvp,it)
            if (movion==1) call wrdata(parti,nppi,it)
            if (movion==0) then
               isign = -1
               call pfmove(qi,noff,nyp,isign,tfmove,kyp,kstrt,nvp,idps,n&
     &terf,ierr,inorder)
               if (ierr /= 0) then
                  call MP_END
                  call PPEXIT
                  stop
               endif
               call wrdata(qi,nvp,it)
               isign = 1
               call pfmove(qi,noff,nyp,isign,tfmove,kyp,kstrt,nvp,idps,n&
     &terf,ierr,inorder)
               if (ierr /= 0) then
                  call MP_END
                  call PPEXIT
                  stop
               endif
               call zguard(qi,nyp,nx,inorder)
            endif
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
         call displayw(wt,zero,dt*real(ntw),it,irc)
! check error return code
         if (irc==1) go to 3000
      endif
      call pwtimer(time,dtime)
! send main CPU Time to diagnostic nodes
      msg(1:2) = time; msg(3) = tpush; msg(4) = tdpost; msg(5) = tsort
      msg(6:7) = tmove; msg(8:9) = tfft; msg(10) = tfmove
      msg(11) = trepart
      call HARTBEAT(msg,11)
      if (id0==0) then
         write (18,*) 'processor partition used: nvp = ', nvp
         write (18,*) ncpus, ' processors found, ', ntasks+1, ' used'
         write (18,*) 'main max/min real time=', time(1), time(2), 'sec'
         totpush = tpush + tdpost
         write (18,*) 'electron push time=', tpush, 'sec'
         write (18,*) 'electron charge deposit time=', tdpost, 'sec'
         write (18,*) 'total electron push time=', totpush, 'sec'
         write (18,*) 'electron sort time=', tsort, 'sec'
         write (18,*) 'electron move time=', tmove, 'sec'
         totpush = totpush + tsort + tmove(1)
         write (18,*) 'total electron time=', totpush, 'sec'
      endif
      if (movion==1) then
         msg(1) = tpushi; msg(2) = tdposti; msg(3) = tsorti
         msg(4:5) = tmovi
         call HARTBEAT(msg,5)
         if (id0==0) then
            totpushi = tpushi + tdposti
            write (18,*) 'ion push time=', tpushi, 'sec'
            write (18,*) 'ion charge deposit time=', tdposti, 'sec'
            write (18,*) 'total ion push time=', totpushi, 'sec'
            write (18,*) 'ion sort time=', tsorti
            write (18,*) 'ion move time=', tmovi
            totpushi = totpushi + tsorti + tmovi(1)
            write (18,*) 'total ion time=', totpushi, 'sec'
         endif
      endif
      if (id0==0) then
         write (18,*) 'total fft time=', tfft, 'sec'
         write (18,*) 'total partition time=', tfmove, 'sec'
         write (18,*) 'total repartition time=', trepart, 'sec'
         tfmove = tfmove + trepart
         time(1) = time(1) - (totpush + totpushi + tfft(1) + tfmove)
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
         if (irc /= 0) write (18,*) 'pinput2 namelist not written'! done
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
! allocate data for restart and/or phase space diagnostic
         if ((nts > 0) .or. (nustrt /= 1) .or. (ntr > 0)) then
            allocate(part(idimp,max(npmax,npimax),nblok))
            allocate(npp(nblok))
            allocate(edges(idps,nblok))
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
            call rddata(edges,nvp,it,ierr)
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
               call displayw(wt,zero,dt*real(ntw),it,irc)
! check error return code
               if (irc==1) go to 20
            endif
! get main CPU Time
            call HARTBEAT(msg,11)
            time = msg(1:2); tpush = msg(3); tdpost = msg(4)
            tsort = msg(5); tmove = msg(6:7); tfft = msg(8:9)
            tfmove = msg(10); trepart = msg(11)
            if (id0==0) then
               write (18,*) 'processor partition used: nvp = ', nvp
               write (18,*) ncpus, ' processors found, ', ntasks+1, ' us&
     &ed'
               write (18,*) 'main max/min real time=', time(1), time(2),&
     &'sec'
               totpush = tpush + tdpost
               write (18,*) 'electron push time=', tpush, 'sec'
               write (18,*) 'electron charge deposit time=', tdpost, 'se&
     &c'
               write (18,*) 'total electron push time=', totpush, 'sec'
               write (18,*) 'electron sort time=', tsort, 'sec'
               write (18,*) 'electron move time=', tmove, 'sec'
               totpush = totpush + tsort + tmove(1)
               write (18,*) 'total electron time=', totpush, 'sec'
            endif
            if (movion==1) then
               call HARTBEAT(msg,5)
               tpushi = msg(1); tdposti = msg(2); tsorti = msg(3)
               tmovi = msg(4:5)
               if (id0==0) then
                  totpushi = tpushi + tdposti
                  write (18,*) 'ion push time=', tpushi, 'sec'
                  write (18,*) 'ion charge deposit time=', tdposti, 'sec&
     &'
                  write (18,*) 'total ion push time=', totpushi, 'sec'
                  write (18,*) 'ion sort time=', tsorti
                  write (18,*) 'ion move time=', tmovi
                  totpushi = totpushi + tsorti + tmovi(1)
                  write (18,*) 'total ion time=', totpushi, 'sec'
               endif
            endif
            if (id0==0) then
               write (18,*) 'total fft time=', tfft, 'sec'
               write (18,*) 'total partition time=', tfmove, 'sec'
               write (18,*) 'total repartition time=', trepart, 'sec'
               tfmove = tfmove + trepart
               time(1) = time(1) - (totpush + totpushi + tfft(1) + tfmov&
     &e)
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
! velocity diagnostic
         if (ntv > 0) then
            it = itime/ntv
            if (itime==ntv*it) then
! display velocity distributions
               call displayfv(fv,fvm,' ELECTRON',itime,nmv,2,2,irc)
               if (irc==1) go to 10
! print out velocity moments
               if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)
               if (movion==1) then
! display velocity distributions
                  call displayfv(fvi,fvmi,' ION',itime,nmv,2,2,irc)
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
! debug: only periodic writes currently supported
               if (ipbc.eq.1) then
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
         endif
! energy diagnostic
         if (ntw > 0) then
            it = itime/ntw
            if (itime==ntw*it) then
! get energy values
               call HARTBEAT(msg,4)
               wtot = msg(1:4)
               if (id0==0) write (18,992) wtot(1), wtot(2), wtot(4)
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
               call wrdata(edges,nvp,it)
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
      end program pbeps2
