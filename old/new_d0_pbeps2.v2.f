!-----------------------------------------------------------------------! * * * periodic 2d electrostatic particle simulation kernel code * * *! this is a simple 2d skeleton particle-in-cell code designed for! exploring new computer architectures.  it contains the critical pieces! needed for depositing charge, advancing particles, and solving the! field.  the code moves only electrons, with periodic electrostatic! forces obtained by solving poisson's equation with fast fourier! transforms.  the only diagnostic is particle and field energy.! portable gcpic kernel code, using algorithm described in:! p. c. liewer and v. k. decyk, j. computational phys. 85, 302 (1989).! written by viktor k. decyk, ucla! for mpi distributed memory computers! update: january 19, 2004      program pbeps2      use pinit2d!     use mprbpush2d!     use mppush2d      use prbpush2d      use ppush2d      use pdfield2d!     use mpfield2d      use pdiag2d      use p2d      implicit none! idps = number of partition boundaries! idimp = dimension of phase space = 4! mshare = (0,1) = (no,yes) architecture is shared memory! nmv = number of segments in v for velocity distribution      integer :: idps =    2, idimp =   4, mshare =   0, nmv = 40      integer :: vect = 0      integer :: npxy, npxyb, np, npxyi, npxybi, npi      integer :: nx, ny, nxh, nyh, nyv, nxe, nxeh      integer :: ipbc, nloop, nvp, nblok, npmax, nypmx, npimax = 0      integer :: kyp, kxp, kyb, kxb, kxyb, kbmin, kblok, jbmin, jblok      integer :: ngds, nxyh, nxhy, nx1, nbmax      integer :: idproc, id0, kstrt, itime, isign, isc, irc, ierr      integer :: it, lprec, nrec      integer :: ntasks      real :: zero = 0.0, wki = 0.      real :: qbme, qbmi, affp, qi0, etx, we, wke      real :: vtxi, vtyi, vtdxi, vtdyi      double precision :: dtime      real, dimension(:,:,:), allocatable :: part, parti      real, dimension(:,:,:), allocatable :: qe, qi      real, dimension(:,:,:,:), allocatable :: fxye      complex, dimension(:,:,:), allocatable :: qt      complex, dimension(:,:,:,:), allocatable :: fxyt      complex, dimension(:,:,:), allocatable :: ffc      integer, dimension(:), allocatable :: mixup      complex, dimension(:), allocatable :: sct      real, dimension(:,:), allocatable  :: edges      integer, dimension(:), allocatable :: nyp, noff      integer, dimension(:), allocatable :: npp, nppi, nps      real, dimension(:,:,:), allocatable :: sfield      complex, dimension(:,:,:), allocatable :: sfieldt      real, dimension(:,:,:), allocatable :: fv, fvm, fvi, fvmi      real, dimension(:,:), allocatable :: wt! dirichlet boundary conditions      integer :: indx1, indy1, nx2, ny2, nx2e, nxhy2, nxyh2      integer :: kyp2, kxp2, kyb2, kxb2, kxyb2      integer :: kbmin2, k2blok, jbmin2, j2blok      real, dimension(:,:,:), allocatable :: q2, sfield2      real, dimension(:,:,:,:), allocatable :: fxy2      complex, dimension(:,:,:,:), allocatable :: fxyt2      complex, dimension(:,:,:), allocatable :: qt2, sfieldt2      complex, dimension(:,:,:), allocatable :: ffd      integer, dimension(:), allocatable :: mixup2      complex, dimension(:), allocatable :: sct2!      integer, dimension(2) :: ktime! wtot = total energy      real, dimension(4) :: wtot! time = timing array      real, dimension(2) :: time! msg = heartbeat array      double precision, dimension(4) :: msg      character(len=16) :: fname      character(len=12) :: label, potname      character(len=4) :: nlabel  991 format (' T = ',i7)  992 format (' field, kinetic, total energies = ',3e14.7)! nvp = number of real or virtual processors! initialize for parallel processing      call PPINIT(idproc,id0,nvp)      kstrt = idproc + 1! read namelist      if (id0==0) then         open(unit=8,file='pinput2',form='formatted',status='old')         read (8,pinput2)         open(unit=18,file='poutput2',form='formatted',status='unknown')         write (18,pinput2)      endif! broadcast namelist to other nodes      call sendnml()! np = total number of electrons in simulation      npxy = npx*npy; npxyb = npxb*npyb; np = npxy + npxyb! npi = total number of ions in simulation      npxyi = npxi*npyi; npxybi = npxbi*npybi; npi = npxyi + npxybi      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = ny/2      nxe = nx + 4; nyv = ny + 2! kyp = number of complex grids in each field partition in y direction! nypmx = maximum size of particle partition, including guard cells.      kyp = (ny - 1)/nvp + 1; nypmx = kyp + 3! ngds = number of guard cells      ngds = 3*((idps - 1)/2 + 1)      if (inorder==LINEAR) then         ax = .912871; ay = .912871         nxe = nx + 2; ; nypmx = kyp + 1         ngds = (idps - 1)/2 + 1      endif! check if too many processors      if (nvp > ny) then         write (2,*) 'Too many processors requested, ny, nvp=', ny, nvp         call PPEXIT         stop      endif! boundary conditions      ipbc = psolve! initialize for multiprocessing!     call mpinit(ntasks)!     write (18,*) ntasks+1, ' processors found and used'      if (popt==VECTOR) vect = 1      nxeh = nxe/2! nloop = number of time steps in simulation      nloop = tend/dt + .0001! nblok = number of particle partitions      nblok = 1 + mshare*(nvp - 1)! npmax = maximum number of electrons in each partition      npmax = (np/nvp)*1.01 + 7000      if (movion==1) npimax = (npi/nvp)*1.01 + 7000! kxp = number of complex grids in each field partition in x direction      kxp = (nxh - 1)/nvp + 1! kyb = number of processors in y! kxb = number of processors in x      kyb = ny/kyp; kxb = nxh/kxp! kxyb = maximum(kxb,kyb)      kxyb = max(kxb,kyb)! kblok = number of field partitions in y direction      kbmin = 1 + (1 - mshare)*(kxyb/kxb - 1)      kblok = 1 + mshare*(ny/kyp - 1)! jblok = number of field partitions in x direction      jbmin = 1 + (1 - mshare)*(kxyb/kyb - 1)      jblok = 1 + mshare*(nxh/kxp - 1)! nxyh = maximum(nx,ny)/2      nxyh = max(nx,ny)/2! nxhy = maximum(nx/2,ny)      nxhy = max(nxh,ny)! dimensions for index and sorting arrays      nx1 = nx + 1! nbmax = size of buffer for passing particles between processors      nbmax = 1 + (2*(npxy*vty + npxyb*vtdy) + 1.4*npxyb*abs(vdy))*dt/ny      if (movion==1) then         vtxi = vtx/sqrt(rmass*rtempxi)         vtyi = vty/sqrt(rmass*rtempyi)      endif! diagnostic information needed by diagnostic nodes! velocity diagnostics      if (ntv > 0) then         allocate(fv(2*nmv+2,2,nblok),fvm(2,2,nblok))         if (id0==0) then            open(unit=10,file='fv2',form='formatted',status='unknown')! write captions            write (10,*) 'it vdx vdy vtx vty'         endif         if (movion==1) then            allocate(fvi(2*nmv+2,2,nblok),fvmi(2,2,nblok))            if (id0==0) then             open(unit=20,file='fvi2',form='formatted',status='unknown')! write captions             write (20,*) 'it vdxi vdyi vtxi vtyi'            endif         endif      endif! density or potential diagnostics      if ((ntp > 0) .or. (ntd > 0)) then         allocate(sfield(nxe,nypmx*kbmin,kblok))         inquire(iolength=lprec) sfield(1,1,1)         lprec = nx1*kyp*lprec      endif! potential diagnostics      if (ntp > 0) then         potname = 'pot2'         if (id0==0) call PHEAD2(11,nx1,kyp,nvp,trim(potname))      endif! energy diagnostics      if (ntw > 0) allocate(wt((nloop-1)/ntw+1,4))! open restart files      if ((nustrt /= 1) .or. (ntr > 0)) then         if (id0==0) then            open(unit=16,file='rstrt1',form='unformatted',status='unknow&     &n')            open(unit=17,file='rstrt2',form='unformatted',status='unknow&     &n')         endif      endif! open graphics device      call GROPEN      call SETNPLT(nplot,irc)      call STPALIT(idpal)!! diagnostic nodes have special processing      if (idproc < 0) call diag2nodes!! part(1,n,l) = position x of particle n in partition l! part(2,n,l) = position y of particle n in partition l! part(3,n,l) = velocity vx of particle n in partition l! part(4,n,l) = velocity vy of particle n in partition l      allocate(part(idimp,npmax,nblok))! in real space, qe(j+1,k,l) = charge density at grid point (j,kk)! in real space, fxye(i,j+1,k,l) = i component of force/charge at ! grid point (j,kk)! in other words, fxye are the convolutions of the electric field! over the particle shape, where kk = k + noff(l) - 1      allocate(qe(nxe,nypmx*kbmin,kblok),fxye(2,nxe,nypmx*kbmin,kblok))! qt(k,j,l) = complex charge density for fourier mode jj-1,k-1! fxyt(1,k,j,l) = x component of force/charge for fourier mode jj-1,k-1! fxyt(2,k,j,l) = y component of force/charge for fourier mode jj-1,k-1! where jj = j + kxp*(l - 1)      allocate(qt(nyv,kxp,jblok),fxyt(2,nyv,kxp,jblok))! ffc = form factor array for poisson solver      allocate(ffc(nyh,kxp,jblok))! mixup, sct = arrays for fft      allocate(mixup(nxhy),sct(nxyh))! edges(1,l) = lower boundary of particle partition l! edges(2,l) = upper boundary of particle partition l      allocate(edges(idps,nblok))! nyp(l) = number of primary gridpoints in particle partition l.! noff(l) = lowermost global gridpoint in particle partition l.      allocate(nyp(nblok),noff(nblok))! npp(l) = number of particles in partition l! nps(l) = starting address of particles in partition l      allocate(npp(nblok),nps(nblok))! dirichlet boundary conditions      if (ipbc==2) then         indx1 = indx+1; indy1 = indy+1; nx2 = 2*nx; ny2 = 2*ny         nx2e = 2*nxe; nxhy2 = 2*nxhy; nxyh2 = 2*nxyh         kyp2 = (ny2-1)/nvp + 1; kxp2 = (nx-1)/nvp + 1         kyb2 = ny2/kyp2; kxb2 = nx/kxp2         kxyb2 = max(kxb2,kyb2)         kbmin2 = 1 + (1 - mshare)*(kxyb2/kxb2 - 1)         k2blok = 1 + mshare*(ny2/kyp2 - 1)         jbmin2 = 1 + (1 - mshare)*(kxyb2/kyb2 - 1)         j2blok = 1 + mshare*(nx/kxp2 - 1)         allocate(q2(nx2e,kyp2*kbmin2,k2blok))         allocate(fxy2(2,nx2e,kyp2*kbmin2,k2blok))         allocate(sfield2(nx2e,kyp2*kbmin2,k2blok))         allocate(qt2(ny2,kxp2,j2blok),fxyt2(2,ny2,kxp2,j2blok))         allocate(sfieldt2(ny2,kxp2,j2blok))         allocate(ffd(ny,kxp2,j2blok))         allocate(mixup2(nxhy2),sct2(nxyh2))      endif!! initialize parallel timer      call pwtimer(time,dtime,-1)! initialize constants      itime = 0      qbme = qme      if (ipbc==1) then         affp = float(nx*ny)/float(np)      else if (ipbc==2) then         affp = float((nx-2)*(ny-2))/float(np)      endif      if (movion==1) then         qbmi = qmi/rmass         vtdxi = vtx/sqrt(rmass*rtempdxi)         vtdyi = vty/sqrt(rmass*rtempdyi)      endif! diagnostics! velocity diagnostics      if (ntv > 0) then         fv(1,:,:) = 8.*max(vtx,vty)         if (movion==1) fvi(1,:,:) = 8.*max(vtxi,vtyi)      endif! density or potential diagnostics      if ((ntp > 0) .or. (ntd > 0)) then         allocate(sfieldt(nyv,kxp,jblok))      endif! calculate partition variables      call dcomp(edges,nyp,noff,ny,kstrt,nvp,inorder)! prepare fft tables      call fft_init(mixup,sct,indx,indy)! calculate form factors      call pois_init(ffc,ax,ay,affp,nx,ny,kstrt)! dirichlet boundary conditions      if (ipbc==2) then         call fft_init(mixup2,sct2,indx1,indy1)         call poisd_init(ffd,ax,ay,affp,nx,ny,kstrt)      endif!! allocate background charge density      if (movion==0) allocate(qi(nxe,nypmx*kbmin,kblok))! allocate ion data      if (movion==1) allocate(parti(idimp,npimax,nblok),nppi(nblok))! debug      if (ipbc==2) vdy = 0.! new start      if (nustrt==1) then! initialize ions         if (movion==1) then            nps = 1            nppi = 0! background ions            if (npxyi > 0) call distr(parti,edges,nppi,nps,vtxi,vtyi,vxi&     &0,vyi0,npxi,npyi,nx,ny,ipbc)!           if (npxyi > 0) then!              call sdistr(parti,nps,zero,zero,zero,zero,zero,zero,npxi,&!    &npyi,nx,ny,kstrt,nvp,ipbc)!              call vdistr(parti,nppi,nps,vtxi,vtyi,vxi0,vyi0,npxi,npyi,&!    &kstrt,nvp)!           endif! beam ions            nps = nppi + 1            if (npxybi > 0) call distr(parti,edges,nppi,nps,vtdxi,vtdyi,&     &vdxi,vdyi,npxbi,npybi,nx,ny,ipbc)!           if (npxybi > 0) then!              call sdistr(parti,nps,zero,zero,zero,zero,zero,zero,npxbi&!    &,npybi,nx,ny,kstrt,nvp,ipbc)!              call vdistr(parti,nppi,nps,vtdxi,vtdyi,vdxi,vdyi,npxbi,np&!    &ybi,kstrt,nvp)!           endif! move ions into appropriate spatial regions            call pmove(parti,edges,nppi,ny,kstrt,nvp,nbmax,vect,ierr)            if (ierr /= 0) then!              call MP_END               call PPEXIT               stop            endif         endif! initialize electrons         nps = 1         npp = 0! background electrons         if (npxy > 0) call distr(part,edges,npp,nps,vtx,vty,vx0,vy0,npx&     &,npy,nx,ny,ipbc)!        if (npxy > 0) then!           call sdistr(part,nps,zero,zero,zero,zero,zero,zero,npx,npy,n&!    &x,ny,kstrt,nvp,ipbc)!           call vdistr(part,npp,nps,vtx,vty,vx0,vy0,npx,npy,kstrt,nvp)!        endif! beam electrons         nps = npp + 1         if (npxyb > 0) call distr(part,edges,npp,nps,vtdx,vtdy,vdx,vdy,&     &npxb,npyb,nx,ny,ipbc)!        if (npxyb > 0) then!           call sdistr(part,nps,zero,zero,zero,zero,zero,zero,npxb,npyb&!    &,nx,ny,kstrt,nvp,ipbc)!           call vdistr(part,npp,nps,vtdx,vtdy,vdx,vdy,npxb,npyb,kstrt,n&!    &vp)!        endif! move electrons into appropriate spatial regions         call pmove(part,edges,npp,ny,kstrt,nvp,nbmax,vect,ierr)         if (ierr /= 0) then!           call MP_END            call PPEXIT            stop         endif! initialize background charge density         if (movion==0) then            qi0 = -qme/affp            call sguardp(qi,kstrt,nvp,noff,nyp,qi0,nx,ny,ipbc,inorder)!           call sguardp(qi,kstrt,nvp,nyp,qi0,nx,ipbc,inorder)!           call sguardp(qi,kstrt,nvp,nyp,zero,nx,ipbc,inorder)!           call dpost(part,qi,-qme,npp,noff,inorder,dopt)! freeze the ions now         else if ((movion==1).and.(itime==ionoff)) then            allocate(qi(nxe,nypmx*kbmin,kblok))! initialize ion charge density to zero            call sguardp(qi,kstrt,nvp,noff,nyp,zero,nx,ny,ipbc,inorder)!           call sguardp(qi,kstrt,nvp,nyp,zero,nx,ipbc,inorder)! deposit ion charge            call dpost(parti,qi,qmi,nppi,noff,inorder,dopt)!           call dpost(parti,qi,qmi,nppi,nps,noff,ntasks,inorder,dopt)! delete ions            deallocate(parti,nppi)            movion = 0         endif! restart      else! determine most recent restart file         if (id0==0) then            read (16,iostat=ierr) ktime(1)            if (ierr /= 0) ktime(1) = -1            read (17,iostat=ierr) ktime(2)            if (ierr /= 0) ktime(2) = -1            if (ktime(1) > ktime(2)) then               ktime(2) = 16            else               ktime(1) = ktime(2)               ktime(2) = 17            endif         endif         call plbcast(ktime)         itime = ktime(1)         if (itime < 0) go to 400! read restart file         it = ktime(2)         call rddata(part,npp,it,ierr)         if (ierr /= 0) go to 400         if (movion==1) then            call rddata(parti,nppi,it,ierr)            if (ierr /= 0) go to 400         endif         if (movion==0) then            call rddata(qi,nvp,it,ierr)            if (ierr /= 0) go to 400         endif         if (ntw > 0) then            call rddata(wt,1,it,ierr)            if (ierr /= 0) go to 400            call plbcast(wt)         endif         if (id0==0) then            read (it,iostat=ierr) ktime(1)            if (ierr /= 0) ktime(1) = -1            rewind it         endif         call plbcast(ktime)         irc = ktime(1)         if (irc==itime) go to 490! handle error  400    write (18,*) 'Restart Error'         go to 3000      endif! record time  490 call pwtimer(time,dtime)! send initial CPU Time to diagnostic nodes      msg(1) = time(1); msg(2) = time(2)      call HARTBEAT(msg,2)      if (id0==0) then         write (18,*) 'init max/min real time=', time(1), time(2), 'sec'      endif!! * * * start main iteration loop * * *!  500 if (nloop <= itime) go to 2000! send time step to diagnostic nodes      msg(1) = itime      call HARTBEAT(msg,1)      if (id0==0) write (18,991) itime      write (label,991) itime      call LOGNAME(label)! initialize charge density to background      call sguardp(qe,kstrt,nvp,noff,nyp,zero,nx,ny,ipbc,inorder)!     call sguardp(qe,kstrt,nvp,nyp,zero,nx,ipbc,inorder)! deposit charge      call dpost(part,qe,qme,npp,noff,inorder,dopt)!     call dpost(part,qe,qme,npp,nps,noff,ntasks,inorder,dopt)! density diagnostic      if (ntd > 0) then         it = itime/ntd         if (itime==ntd*it) then            sfield = -qe! add guard cells for density in x            call aguardp(sfield,nyp,nx,ipbc,inorder)! add guard cells for density in y            call paguardp(sfield,kstrt,nvp,nx,kyp,ngds,ipbc)         endif      endif      if (movion==1) then         call dpost(parti,qe,qmi,nppi,noff,inorder,dopt)!        call dpost(parti,qe,qmi,nppi,nps,noff,ntasks,inorder,dopt)      else         qe = qe + qi      endif! add guard cells for density in x      call aguardp(qe,nyp,nx,ipbc,inorder)! add guard cells for density in y      call paguardp(qe,kstrt,nvp,nx,kyp,ngds,ipbc)! freeze the ions      if ((movion==1).and.(itime==ionoff)) then         allocate(qi(nxe,nypmx*kbmin,kblok))! initialize ion charge density to zero         call sguardp(qi,kstrt,nvp,noff,nyp,zero,nx,ny,ipbc,inorder)!        call sguardp(qi,kstrt,nvp,nyp,zero,nx,ipbc,inorder)! deposit ion charge         call dpost(parti,qi,qmi,nppi,noff,inorder,dopt)!        call dpost(parti,qi,qmi,nppi,nps,noff,ntasks,inorder,dopt)! delete ions         deallocate(parti,nppi)         movion = 0      endif! velocity diagnostic      if (ntv > 0) then         it = itime/ntv         if (itime==ntv*it) then! calculate electron distribution function and moments            call vdist(part,fv,fvm,npp,nmv,2)            call plsum(fv(:,:,1))            fv(1,:,:) = 8.*max(vtx,vty)! display velocity distributions            call displayfv(fv,fvm,' ELECTRON',itime,nmv,2,2,irc)            if (irc==1) go to 2000! print out velocity moments            if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)            if (movion==1) then! calculate ion distribution function and moments               call vdist(parti,fvi,fvmi,nppi,nmv,2)               call plsum(fvi(:,:,1))               fvi(1,:,:) = 8.*max(vtxi,vtyi)! display velocity distributions               call displayfv(fvi,fvmi,' ION',itime,nmv,2,2,irc)               if (irc==1) go to 2000! print out velocity moments               if (id0==0) write (20,*) it, fvmi(1,:,:), fvmi(2,:,:)            endif         endif      endif! phase space diagnostic      if (nts > 0) then         it = itime/nts         if (itime==nts*it) then            isc = 999! plot electrons vx versus x            call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,ny,&     &3,1,irc)            if (irc==1) go to 2000! plot electrons vy versus y            call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,ny,&     &4,2,irc)            if (irc==1) go to 2000            if (movion==1) then! plot ions vx versus x               call grasp(parti,nppi,' ION PHASE SPACE',itime,isc,nx,ny,&     &3,1,irc)               if (irc==1) go to 2000! plot ions vy versus y               call grasp(parti,nppi,' ION PHASE SPACE',itime,isc,nx,ny,&     &4,2,irc)               if (irc==1) go to 2000            endif         endif      endif!! dirichlet boundary conditions!      if (ipbc==2) then! density diagnostic      if (ntd > 0) then         it = itime/ntd         if (itime==ntd*it) then! transform electron density to fourier space            call dblsin(sfield,q2,nx,ny,kstrt,kyp,kyp2,inorder)            isign = -1            call fft(q2,qt2,isign,mixup2,sct2,indx1,indy1,kstrt,kyp2,LIN&     &EAR)!           call fft(q2,qt2,isign,mixup2,sct2,indx1,indy1,kstrt,kyp2,ntas&!    &ks,LINEAR)! calculate smoothing in fourier space            call poisd(qt2,sfieldt2,ffd,nx,ny,kstrt)! transform electron density to real space            isign = 1            call fft(sfield2,sfieldt2,isign,mixup2,sct2,indx1,indy1,kstr&     &t,kyp2,LINEAR)!           call fft(sfield2,sfieldt2,isign,mixup2,sct2,indx1,indy1,kstr&!    &t,kyp2,ntasks,LINEAR)            call hafdbl(sfield,sfield2,nx,ny,kstrt,kyp,kyp2,inorder)            call plcguard(sfield,kstrt,nvp,nx,kyp,inorder)            call lcguard(sfield,nyp,nx,inorder)! display electron density            call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle,nx&     &,ny,irc,inorder)            if (irc==1) go to 2000         endif      endif! transform charge to fourier space      call dblsin(qe,q2,nx,ny,kstrt,kyp,kyp2,inorder)      isign = -1      call fft(q2,qt2,isign,mixup2,sct2,indx1,indy1,kstrt,kyp2,LINEAR)!     call fft(q2,qt2,isign,mixup2,sct2,indx1,indy1,kstrt,kyp2,ntasks,LI&!    &NEAR)! potential diagnostic      if (ntp > 0) then         it = itime/ntp         if (itime==ntp*it) then! solve for potential            call poisd(qt2,sfieldt2,ffd,we,nx,ny,kstrt)            isign = 1            call fft(sfield2,sfieldt2,isign,mixup2,sct2,indx1,indy1,kstr&     &t,kyp2,LINEAR)!           call fft(sfield2,sfieldt2,isign,mixup2,sct2,indx1,indy1,kstr&!    &t,kyp2,ntasks,LINEAR)            call hafdbl(sfield,sfield2,nx,ny,kstrt,kyp,kyp2,inorder)            call plcguard(sfield,kstrt,nvp,nx,kyp,inorder)            call lcguard(sfield,nyp,nx,inorder)! display potential            call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle,nx&     &,ny,irc,inorder)            if (irc==1) go to 2000!           g(:,1) = sfield(nxh,:,1)!           call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorder)!           if (irc==1) go to 2000! write diagnostic output            write (nlabel,'(i4)') it            fname = trim(potname)//'_'//trim(adjustl(nlabel))            nrec = -lprec!           call writef(sfield,nx1,kyp,11,nrec,trim(fname),inorder)         endif      endif! calculate force/charge in fourier space      call poisd(qt2,fxyt2,ffd,we,nx,ny,kstrt)      call fft(fxy2,fxyt2,mixup2,sct2,indx1,indy1,kstrt,kyp2,LINEAR)!     call fft(fxy2,fxyt2,mixup2,sct2,indx1,indy1,kstrt,kyp2,ntasks,LINE&!    &AR)      call hafdbl(fxye,fxy2,nx,ny,kstrt,kyp,kyp2,inorder)      call plcguard(fxye,kstrt,nvp,nx,kyp,inorder)      call lcguard(fxye,nyp,nx,inorder)!! periodic boundary conditions!      else if (ipbc==1) then! density diagnostic      if (ntd > 0) then         it = itime/ntd         if (itime==ntd*it) then! transform density to fourier space            isign = -1            call fft(sfield,qt,isign,mixup,sct,indx,indy,kstrt,kyp,inord&     &er)!           call fft(sfield,qt,isign,mixup,sct,indx,indy,kstrt,kyp,ntask&!    &s,inorder)! calculate smoothing in fourier space            call pois(qt,sfieldt,ffc,nx,ny,kstrt)! transform density to real space            isign = 1            call fft(sfield,sfieldt,isign,mixup,sct,indx,indy,kstrt,kyp,&     &inorder)!           call fft(sfield,sfieldt,isign,mixup,sct,indx,indy,kstrt,kyp,&!    &ntasks,inorder)! copy to guard cells            call pcguard(sfield,kstrt,nvp,kyp,inorder)            call cguard(sfield,nyp,nx,inorder)! display density            call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle,nx&     &,ny,irc,inorder)            if (irc==1) go to 2000         endif      endif! transform charge to fourier space      isign = -1      call fft(qe,qt,isign,mixup,sct,indx,indy,kstrt,kyp,inorder)!     call fft(qe,qt,isign,mixup,sct,indx,indy,kstrt,kyp,ntasks,inorder)! potential diagnostic      if (ntp > 0) then         it = itime/ntp         if (itime==ntp*it) then! calculate potential in fourier space            call pois(qt,sfieldt,ffc,we,nx,ny,kstrt)! transform potential to real space            isign = 1            call fft(sfield,sfieldt,isign,mixup,sct,indx,indy,kstrt,kyp,&     &inorder)!           call fft(sfield,sfieldt,isign,mixup,sct,indx,indy,kstrt,kyp,&!    &ntasks,inorder)! copy to guard cells            call pcguard(sfield,kstrt,nvp,kyp,inorder)            call cguard(sfield,nyp,nx,inorder)! display potential            call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle,nx&     &,ny,irc,inorder)            if (irc==1) go to 2000!           g(:,1) = sfield(nxh,:,1)!           call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorder)!           if (irc==1) go to 2000! write diagnostic output            write (nlabel,'(i4)') it            fname = trim(potname)//'_'//trim(adjustl(nlabel))            nrec = -lprec!           call writef(sfield,nx1,kyp,11,nrec,trim(fname),inorder)         endif      endif! calculate force/charge in fourier space      call pois(qt,fxyt,ffc,we,nx,ny,kstrt)! transform force/charge to real space      call fft(fxye,fxyt,mixup,sct,indx,indy,kstrt,kyp,inorder)!     call fft(fxye,fxyt,mixup,sct,indx,indy,kstrt,kyp,ntasks,inorder)! copy data from field to particle partition, and copy to guard cells      call pcguard(fxye,kstrt,nvp,kyp,inorder)      call cguard(fxye,nyp,nx,inorder)!      endif!! external pump      if ((itpon > 0).and.(itime >= itpon)) then         etx = (v0*vtx)*w0*cos(w0*dt*(itime - itpon))         fxye(1,:,:,:) = fxye(1,:,:,:) + etx      endif! particle push and charge density update      wke = 0.! push electrons      if (relativity==1) then         call push(part,fxye,npp,noff,qbme,dt,ci,wke,nx,ny,ipbc,inorder,&     &popt)!        call push(part,fxye,npp,nps,noff,qbme,dt,ci,wke,nx,ny,ipbc,ntas&!    &ks,inorder,popt)      else         call push(part,fxye,npp,noff,qbme,dt,wke,nx,ny,ipbc,inorder,pop&     &t)!        call push(part,fxye,npp,nps,noff,qbme,dt,wke,nx,ny,ipbc,ntasks,&!    &inorder,popt)      endif! move electrons into appropriate spatial regions      call pmove(part,edges,npp,ny,kstrt,nvp,nbmax,vect,ierr)      if (ierr /= 0) then!        call MP_END         call PPEXIT         stop      endif! push ions      if (movion==1) then         wki = 0.         if (relativity==1) then            call push(parti,fxye,nppi,noff,qbmi,dt,ci,wki,nx,ny,ipbc,ino&     &rder,popt)!           call push(parti,fxye,nppi,nps,noff,qbmi,dt,ci,wki,nx,ny,ipbc&!    &,ntasks,inorder,popt)         else            call push(parti,fxye,nppi,noff,qbmi,dt,wki,nx,ny,ipbc,inorde&     &r,popt)!           call push(parti,fxye,nppi,nps,noff,qbmi,dt,wki,nx,ny,ipbc,nt&!    &asks,inorder,popt)         endif         wki = wki*rmass! move ions into appropriate spatial regions         call pmove(parti,edges,nppi,ny,kstrt,nvp,nbmax,vect,ierr)         if (ierr /= 0) then!           call MP_END            call PPEXIT            stop         endif      endif! sort electrons      if (sortime > 0) then         if (mod(itime,sortime)==0) then            call sortp(part,npp,noff,kyp,inorder)         endif      endif! sort ions      if ((movion==1).and.(sortimi > 0)) then         if (mod(itime,sortimi)==0) then            call sortp(parti,nppi,noff,kyp,inorder)         endif      endif! energy diagnostic      if (ntw > 0) then         it = itime/ntw         if (itime==ntw*it) then            wtot(1) = we            wtot(2) = wke            wtot(3) = wki            wtot(4) = we + wke + wki            call plsum(wtot)! send energy values to diagnostic node            msg = wtot            call HARTBEAT(msg,4)            if (id0==0) write (18,992) wtot(1), wtot(2), wtot(4)            wt(it+1,:) = wtot         endif      endif      itime = itime + 1! restart file      if (ntr > 0) then         it = itime/ntr         if (itime==ntr*it) then            it = 16 + mod(it-1,2)            if (id0==0) write (it) itime            call wrdata(part,npp,it)            if (movion==1) call wrdata(parti,nppi,it)            if (movion==0) call wrdata(qi,nvp,it)            if (ntw > 0) call wrdata(wt,1,it)            if (id0==0) then               write (it) itime               end file it               rewind it            endif         endif      endif      go to 500 2000 continue!! * * * end main iteration loop * * *!! send QUIT message to diagnostic nodes      msg = -1.      call HARTBEAT(msg,1)! energy diagnostic      if (ntw > 0) then         it = (itime - 1)/ntw + 1         call displayw(wt,dt*real(ntw),it,irc)! check error return code         if (irc==1) go to 3000      endif      call pwtimer(time,dtime)! send main CPU Time to diagnostic nodes      msg(1) = time(1); msg(2) = time(2)      call HARTBEAT(msg,2)      if (id0==0) then         write (18,*) 'main max/min real time=', time(1), time(2), 'sec'         write (18,*) '* * * q.e.d. * * *'      endif! close graphics device 3000 call PGRCLOSE!     call MP_END      call PPEXIT      stop!      contains!         subroutine diag2nodes         implicit none! diagnostic nodes have special processing  991    format (' T = ',i7)  992    format (' field, kinetic, total energies = ',3e14.7)! allocate data for restart and/or phase space diagnostic         if ((nts > 0) .or. (nustrt /= 1) .or. (ntr > 0)) then            allocate(part(idimp,max(npmax,npimax),nblok))            allocate(npp(nblok))         endif         if (movion==0) allocate(qi(nxe,nypmx*kbmin,kblok))! restart         if (nustrt /= 1) then! determine most recent restart file            if (id0==0) then               read (16,iostat=ierr) ktime(1)               if (ierr /= 0) ktime(1) = -1               read (17,iostat=ierr) ktime(2)               if (ierr /= 0) ktime(2) = -1               if (ktime(1) > ktime(2)) then                  ktime(2) = 16               else                  ktime(1) = ktime(2)                  ktime(2) = 17               endif            endif            call plbcast(ktime)            itime = ktime(1)            if (itime < 0) go to 30! read restart file            it = ktime(2)            call rddata(part,npp,it,ierr)            if (ierr /= 0) go to 30            if (movion==1) then               call rddata(part,npp,it,ierr)               if (ierr /= 0) go to 30            endif            if (movion==0) then               call rddata(qi,nvp,it,ierr)               if (ierr /= 0) go to 30            endif            if (ntw > 0) then               call rddata(wt,1,it,ierr)               if (ierr /= 0) go to 30               call plbcast(wt)            endif            if (id0==0) then               read (it,iostat=ierr) ktime(1)            if (ierr /= 0) ktime(1) = -1               rewind it            endif            call plbcast(ktime)            irc = ktime(1)            if (irc==itime) go to 40! handle error   30       if (id0==0) write (18,*) 'Restart Error'            call PGRCLOSE!           call MP_END            call PPEXIT            stop         endif! get initial CPU Time   40    call HARTBEAT(msg,2)         time(1) = msg(1); time(2) = msg(2)         if (id0==0) then         write (18,*) 'init max/min real time=', time(1), time(2), 'sec'         endif! get time step   10    call HARTBEAT(msg,1)         it = msg(1)         if (it < 0) then! energy diagnostic            if (ntw > 0) then               it = (itime - 1)/ntw + 1               call displayw(wt,dt*real(ntw),it,irc)! check error return code               if (irc==1) go to 20            endif! get main CPU Time            call HARTBEAT(msg,2)            time(1) = msg(1); time(2) = msg(2)            if (id0==0) then               write (18,*) 'main max/min real time=', time(1), time(2),&     &'sec'               write (18,*) '* * * q.e.d. * * *'            endif   20       call PGRCLOSE!           call MP_END            call PPEXIT            stop         else            itime = it         endif         if (id0==0) write (18,991) itime         write (label,991) itime         call LOGNAME(label)! velocity diagnostic         if (ntv > 0) then            it = itime/ntv            if (itime==ntv*it) then! display velocity distributions               call displayfv(fv,fvm,' ELECTRON',itime,nmv,2,2,irc)               if (irc==1) go to 10! print out velocity moments               if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)               if (movion==1) then! display velocity distributions                  call displayfv(fvi,fvmi,' ION',itime,nmv,2,2,irc)                  if (irc==1) go to 10! print out velocity moments                  if (id0==0) write (20,*) it, fvmi(1,:,:), fvmi(2,:,:)              endif            endif         endif! phase space diagnostic         if (nts > 0) then            it = itime/nts            if (itime==nts*it) then               isc = 999! plot electrons vx versus x               call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,&     &ny,3,1,irc)               if (irc==1) go to 10! plot electrons vy versus y               call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,&     &ny,4,2,irc)               if (irc==1) go to 10               if (movion==1) then! plot ions vx versus x                  call grasp(part,npp,' ION PHASE SPACE',itime,isc,nx,ny&     &,3,1,irc)                  if (irc==1) go to 10! plot ions vy versus y                  call grasp(part,npp,' ION PHASE SPACE',itime,isc,nx,ny&     &,4,2,irc)                  if (irc==1) go to 10               endif            endif         endif! density diagnostic         if (ntd > 0) then            it = itime/ntd            if (itime==ntd*it) then! display density               call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle&     &,nx,ny,irc,inorder)               if (irc==1) go to 10            endif         endif! potential diagnostic         if (ntp > 0) then            it = itime/ntp            if (itime==ntp*it) then! display potential               call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle&     &,nx,ny,irc,inorder)               if (irc==1) go to 10!              call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorde&!    &r)!              if (irc==1) go to 10               write (nlabel,'(i4)') it               fname = trim(potname)//'_'//trim(adjustl(nlabel))               nrec = -lprec!              call writef(sfield,nx1,kyp,11,nrec,trim(fname),inorder)            endif         endif! energy diagnostic         if (ntw > 0) then            it = itime/ntw            if (itime==ntw*it) then! get energy values               call HARTBEAT(msg,4)               wtot = msg               if (id0==0) write (18,992) wtot(1), wtot(2), wtot(4)            endif         endif! restart file         itime = itime + 1         if (ntr > 0) then            it = itime/ntr            if (itime==ntr*it) then               it = 16 + mod(it-1,2)               if (id0==0) write (it) itime               call wrdata(part,npp,it)               if (movion==1) call wrdata(part,npp,it)               if (movion==0) call wrdata(qi,nvp,it)               if (ntw > 0) call wrdata(wt,1,it)               if (id0==0) then                  write (it) itime                  end file it                  rewind it               endif            endif         endif         go to 10         end subroutine!      end program pbeps2