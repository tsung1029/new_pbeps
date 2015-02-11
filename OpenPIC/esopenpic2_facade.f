!-----------------------------------------------------------------------! Facade for 2D Electrostatic PIC code with Vacuum Poisson solver! written by viktor k. decyk, ucla! update: april 25, 2008!      module esopenpic2_facade!      use globals, only: LINEAR, QUADRATIC      use pinit2d!     use mppush2d      use ppush2d!     use mprbpush2d      use prbpush2d!     use mpfft2d      use pfft2d      use pdfield2d      use pcfield2d      use pfield2d      use p2d      implicit none      private      public :: init_parallel, allocate_fields, allocate_particles      public :: initialize_particles, initial_partition      public :: particle_manager, deposit_charge, push_particles      public :: repartition, sort_particles      public :: initialize_csolver, add_guardcells, cdensity      public :: cefieldpot, cpotential, cefield! idps = number of partition boundaries! idimp = dimension of phase space = 4      integer :: idps =    2, idimp =   4! ipbc = particle boundary condition = (0,1,2,3) =! (none,2d periodic,2d reflecting,mixed reflecting in x/periodic in y)      integer :: ipbc = 2! lcsolver = (1,2,3) = solve for (potential,efield,both)      integer :: lcsolver = -1      integer :: nx, ny, lvp, kyp      integer :: nbmax, nterg, ngds, nterf      real :: anpav!     real, dimension(:,:), pointer :: part2      real, dimension(:), pointer :: pt      integer, dimension(:), pointer :: ip, npic! vacuum boundary conditions      integer :: indx1, indy1, kxp2, kyp2      real, dimension(:,:), pointer :: q2      real, dimension(:,:,:), pointer :: fxy2      complex, dimension(:,:), pointer :: qt2, pott2      complex, dimension(:,:,:), pointer :: fxyt2      real, dimension(:,:,:), pointer :: ffg      integer, dimension(:), pointer :: mixup2      complex, dimension(:), pointer :: sct2      save!      contains!         subroutine init_parallel(edges,nypu,noffu,indx,indy,kstrt,nvp, &     &inorder)! initializes parallel processing with uniform partitions! with 1D domain decomposition! input: indx, indy, inorder! allocates and returns values for edges, nypu, noffu! returns values for nx, ny, kstrt, nvp!        use mp0d, only: mpinit         implicit none         integer :: indx, indy, kstrt, nvp, inorder!        integer :: ntasks         real, dimension(:), pointer :: edges         integer, dimension(:), pointer :: nypu, noffu! local data         integer :: idproc, id0, lkyp! nvp = number of real or virtual processors! initialize for parallel processing         call PPINIT1(idproc,id0,nvp)         kstrt = idproc + 1! calculate global constants         nx = 2**indx; ny = 2**indy! kyp = number of complex grids in each field partition in y direction         lkyp = (ny - 1)/nvp + 1! edges(1) = lower boundary of particle partition l! edges(2) = upper boundary of particle partition l         allocate(edges(idps))! nyp(1) = number of primary gridpoints in particle partition l.! noff(1) = lowermost global gridpoint in particle partition l.         allocate(nypu(1),noffu(1))! calculate uniform partition variables         call dcomp(edges,nypu,noffu,ny,kstrt,nvp,inorder)! initialize for multiprocessing!        ntasks = mpinit()!        if (id0==0) write (18,*) ntasks+1, ' processors found and used'! save global copy         lvp = nvp         end subroutine init_parallel!         subroutine allocate_fields(qe,qi,fxye,imbalance,inorder)! allocate real space fields with non-uniform 1D domain decomposition! input: inorder! allocates qe, qi, fxye, npic! defines certain global constants:! kyp, nterg, ndgs, nterf         implicit none         integer :: inorder         real :: imbalance         real, dimension(:,:), pointer :: qe, qi         real, dimension(:,:,:), pointer :: fxye! local data         integer :: nxh, kxp, kyb, kxb, kbmin         integer :: nypm1, nxe, nypmx! calculate global constants         nxh = nx/2! kyp = number of complex grids in each field partition in y direction         kyp = (ny - 1)/lvp + 1! kxp = number of complex grids in each field partition in x direction         kxp = (nxh - 1)/lvp + 1! kyb = number of processors in y! kxb = number of processors in x         kyb = ny/kyp; kxb = nxh/kxp         kbmin = max(kxb,kyb)/kxb! size of spatial regions including extra cells for repartitioning         if (imbalance >= 0.0) then            nypm1 = 2*kyp + 1         else            nypm1 = kyp + 1         endif! ngds = number of guard cells this array has! nxe = size of x coordinate, including guard cells (must be even)! nypmx = maximum size of particle partition, including guard cells.! ar = half-width of particle in r direction         if (inorder==LINEAR) then            ngds = (idps - 1)/2 + 1            nxe = nx + 2; nypmx = nypm1         else            ngds = 3*((idps - 1)/2 + 1)            nxe = nx + 4; nypmx = nypm1 + 2         endif! qe = real space charge density         allocate(qe(nxe,nypmx*kbmin))! qi = real space ion background charge density         allocate(qi(nxe,nypmx*kbmin))! fxye = real space force/charge, or convolution of the electric field! over the particle shape         allocate(fxye(2,nxe,nypmx*kbmin))! sorting array         allocate(npic(nypm1))         nterf = 0         nterg = kyp - 1         end subroutine allocate_fields!         subroutine allocate_particles(part,npp,np)! allocate particle data part, npp, pt, ip! with 1D domain decomposition         implicit none         integer :: np         real, dimension(:,:), pointer :: part         integer, dimension(:), pointer :: npp!        real, dimension(:,:), pointer :: part2! local data         integer :: nvp, npmax         nvp = lvp! npmax = maximum number of particles in each partition         npmax = (np/nvp)*1.2! part(1,n) = position x of particle n in partition l! part(2,n) = position y of particle n in partition l! part(3,n) = velocity vx of particle n in partition l! part(4,n) = velocity vy of particle n in partition l         allocate(part(idimp,npmax))! npp(1) = number of particles in partition l         allocate(npp(1))! sorting arrays         allocate(pt(npmax),ip(npmax))!        allocate(part2(idimp,npmax))         end subroutine allocate_particles!         subroutine initialize_particles(part,npp,ampdx,scaledx,shiftdx,&     &ampdy,scaledy,shiftdy,vtx,vty,vx0,vy0,npx,npy,kstrt,ndprof,nsrand,&     &dt)! initialize particles with 1D domain decomposition! defines global constant nbmax         implicit none         integer :: npx, npy, kstrt, ndprof, nsrand         real :: ampdx, scaledx, shiftdx, ampdy, scaledy, shiftdy         real :: vtx, vty, vx0, vy0, dt         real, dimension(:,:), pointer :: part         integer, dimension(:), pointer :: npp! local data         integer :: np, nvp         integer, dimension(size(npp)), target :: lnps         integer, dimension(:), pointer :: nps! initialize particles         np = npx*npy         nvp = lvp         nps => lnps         nps = 1         npp = 0         if (np > 0) then            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&     &ftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)            call vdistr(part,npp,nps,vtx,vty,vx0,vy0,npx,npy,kstrt,nvp)         endif! nbmax = size of buffer for passing particles between processors         nbmax = 1 + (2*np*vty + 1.4*np*abs(vy0))*dt/ny         end subroutine initialize_particles!         subroutine initial_partition(edges,noff,nyp,noffu,nypu,fxye,amp&     &dy,scaledy,shiftdy,kstrt,ndprof,imbalance,inorder)! finds new non-uniform partition boundaries (edges,noff,nyp)! analytically from initial distribution function         implicit none         integer :: kstrt, ndprof, inorder         real :: ampdy, scaledy, shiftdy, imbalance         real, dimension(:), pointer :: edges         integer, dimension(:), pointer :: nyp, noff, nypu, noffu         real, dimension(:,:,:), pointer :: fxye! local data         integer :: nypmx, ierr         nypmx = size(fxye,3)! find non-uniform partition         allocate(noff(size(noffu)),nyp(size(nypu)))         if (imbalance >= 0.0) then            call fedges(edges,noff,nyp,ampdy,scaledy,shiftdy,ny,kstrt,lv&     &p,nypmx,ipbc,ndprof,nterg,ierr,inorder)            if (ierr /= 0) then!              call MP_END               call PPEXIT               stop            endif! use uniform partition         else            noff = noffu; nyp = nypu         endif         end subroutine initial_partition!         subroutine particle_manager(part,edges,npp,kstrt,tmove,pibal)! move particles into appropriate spatial regions! with 1D decomposition         implicit none         integer :: kstrt         real :: tmove         real, dimension(:,:), pointer :: part         real, dimension(:), pointer :: edges         integer, dimension(:), pointer :: npp         real, optional :: pibal! local data         integer :: ierr         real :: lpibal         call pmove(part,edges,npp,anpav,lpibal,tmove,ny,kstrt,lvp,nbmax&     &,ierr)         if (ierr /= 0) then!           call MP_END            call PPEXIT            stop         endif         if (present(pibal)) pibal = lpibal         end subroutine particle_manager!         subroutine deposit_charge(part,qe,qme,tdpost,npp,noff,inorder, &     &dopt)! deposit charge         implicit none         integer :: inorder, dopt         real :: qme, tdpost         real, dimension(:,:), pointer :: part         real, dimension(:,:), pointer :: qe         integer, dimension(:), pointer :: npp         integer, dimension(:), pointer :: noff         call dpost(part,qe,qme,npp,noff,tdpost,inorder,dopt)         end subroutine deposit_charge!         subroutine push_particles(part,fxye,npp,noff,qbme,dt,ci,wke,   &     &tpush,relativity,inorder,popt)! push particles         implicit none         integer :: relativity         real :: qbme, dt, ci, wke, tpush         real, dimension(:,:), pointer :: part         real, dimension(:,:,:), pointer :: fxye         integer, dimension(:), pointer :: npp         integer, dimension(:), pointer :: noff         integer :: inorder, popt         if (relativity==1) then            call rpush(part,fxye,npp,noff,qbme,dt,ci,wke,tpush,nx,ny,ipb&     &c,inorder,popt)         else            call push(part,fxye,npp,noff,qbme,dt,wke,tpush,nx,ny,ipbc,in&     &order,popt)         endif         end subroutine push_particles!         subroutine repartition(part,qi,npp,edges,noff,nyp,kstrt,tmove, &     &tfmove,inorder)! find new non-uniform partitions         implicit none         integer :: kstrt, inorder         real :: tmove, tfmove         real, dimension(:,:), pointer :: part         real, dimension(:,:), pointer :: qi         integer, dimension(:), pointer :: npp         real, dimension(:), pointer :: edges         integer, dimension(:), pointer :: nyp, noff! local data         integer :: nypmx, ierr         real :: pibal         integer, dimension(size(nyp)), target :: lnyps         integer, dimension(size(noff)), target :: lnoffs         integer, dimension(:), pointer :: nyps, noffs         nypmx = size(qi,2)         nyps => lnyps; noffs => lnoffs! count the number of particles per cell         call countp(part,npic,npp,noff,nyp)! save old repartitioning boundaries         noffs = noff         nyps = nyp! determine new repartitioning boundaries         call repart(edges,npic,noff,nyp,anpav,kstrt,lvp,nypmx,nterg,ier&     &r,inorder)         if (ierr /= 0) then!           call MP_END            call PPEXIT            stop         endif         nterf = 0! move background ion density to new field partition.         call pfmove(qi,noff,nyp,noffs,nyps,tfmove,kstrt,lvp,idps,ierr,i&     &norder)         if (ierr /= 0) then!           call MP_END            call PPEXIT            stop         endif! zero out guard cells         call zguard(qi,nyp,nx,inorder)! move particles into new spatial regions         call pmove(part,edges,npp,anpav,pibal,tmove,ny,kstrt,lvp,nbmax,&     &ierr)         if (ierr /= 0) then!           call MP_END            call PPEXIT            stop         endif         end subroutine repartition!         subroutine sort_particles(part,npp,noff,nyp,tsort,inorder)! sort particles         implicit none         integer :: inorder         real :: tsort         real, dimension(:,:), pointer :: part         integer, dimension(:), pointer :: npp         integer, dimension(:), pointer :: nyp, noff         call sortp(part,pt,ip,npp,noff,nyp,npic,tsort,inorder)!        call sortp(part,part2,npp,noff,nyp,npic,tsort,inorder)         end subroutine sort_particles!         subroutine initialize_csolver(qe,pot,ar,affp,indx,indy,csolver,&     &kstrt,nvp,inorder)! initialize 2D Vacuum Poisson solver with 1D domain decomposition! csolver = (1,2,3) = solve for (potential,efield,both)! input: ar, affp, indx, indy, csolver, kstrt, nvp, inorder! the size of qe is used to allocate internal arrays! allocates pot, if requested         implicit none         integer :: indx, indy, csolver, kstrt, nvp, inorder         real :: ar, affp         real, dimension(:,:), pointer :: qe, pot! local data         integer :: nxyh, nxhy, nxe, nypmx! vacuum boundary conditions         integer :: nx2, ny2, nx2e         integer :: kxb2, kyb2, kxyb2, kbmin2! set size of fft tables         nxyh = max(nx,ny)/2; nxhy = max(nx/2,ny)! kyp = number of complex grids in each field partition in y direction         kyp = (ny - 1)/nvp + 1! nypmx = maximum size of particle partition in y, including guard cells         if (inorder==LINEAR) then            nypmx = kyp + 1         else            nypmx = kyp + 3         endif         nxe = size(qe,1)! pot = real space potential         if ((csolver==1) .or. (csolver==3)) then            allocate(pot(nxe,nypmx))         else            nullify(pot)         endif         nypmx = size(qe,2)! vacuum boundary conditions         indx1 = indx + 1; indy1 = indy + 1         nx2 = 2*nx; ny2 = 2*ny; nx2e = 2*nxe         kyp2 = (ny2-1)/nvp + 1; kyb2 = ny2/kyp2         kxp2 = (nx - 1)/nvp + 1; kxb2 = nx/kxp2         kxyb2 = max(kxb2,kyb2)         kbmin2 = kxyb2/kxb2         allocate(q2(nx2e,kyp2*kbmin2))         allocate(qt2(ny2,kxp2))         allocate(ffg(4,ny+1,kxp2+1))         allocate(mixup2(2*nxhy), sct2(2*nxyh))         nullify(pott2,fxy2,fxyt2)         if ((csolver==1) .or. (csolver==3)) then            allocate(pott2(ny2,kxp2))         else            nullify(pott2)         endif         if ((csolver==2) .or. (csolver==3)) then            allocate(fxy2(2,nx2e,kyp2*kbmin2))            allocate(fxyt2(2,ny2,kxp2))         else            nullify(fxy2,fxyt2)         endif! vacuum boundary conditions         call fft_init(mixup2,sct2,indx1,indy1)!        call poisc2_init(ffg,q2,qt2,mixup2,sct2,ar,affp,indx,indy,kstrt&!    &)         call poisc3_init(ffg,q2,qt2,mixup2,sct2,ar,affp,indx,indy,kstrt&     &)! save global copy         lcsolver = csolver         end subroutine initialize_csolver!         subroutine add_guardcells(qi,nyp,kstrt,inorder)! add guard cells with 1D domain decomposition! input: qi, nyp, kstrt, inorder         implicit none         integer :: kstrt, inorder         real, dimension(:,:), pointer :: qi         integer, dimension(:), pointer :: nyp! add guard cells for density in x, only for quadratic interpolation         call laguard(qi,nyp,nx,inorder)! add guard cells for density in y         call pnlaguard(qi,nyp,kstrt,lvp,nx,nterg,ngds)         end subroutine add_guardcells!         subroutine cdensity(qe,noff,nyp,kstrt,tfmove,tfft,inorder)! add guard cells, move to uniform partition and transform density! to fourier space for vacuum boundary conditions! with 1D domain decomposition! input: qe, noff, nyp, kstrt, inorder         implicit none         integer :: kstrt, inorder         real :: tfmove         real, dimension(2) :: tfft         real, dimension(:,:), pointer :: qe         integer, dimension(:), pointer :: noff, nyp! local data         integer :: isign, ierr! return if error         if ((lcsolver < 1) .or. (lcsolver > 3)) return! add guard cells for density in x, only for quadratic interpolation         call laguard(qe,nyp,nx,inorder)! add guard cells for density in y and z         call pnlaguard(qe,nyp,kstrt,lvp,nx,nterg,ngds)! move charge density to uniform field partition         isign = -1         call pfmove(qe,noff,nyp,isign,tfmove,kyp,kstrt,lvp,idps,nterf, &     &ierr,inorder)         if (ierr /= 0) then!           call MP_END            call PPEXIT            stop         endif! copy charge density array qe to triple size array q3         call zdbl(qe,q2,nx,ny,kstrt,kyp,kyp2,inorder)! transform charge to fourier space         isign = -1         call fft(q2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2, &     &LINEAR)         end subroutine cdensity!         subroutine cefieldpot(pot,fxye,noff,nyp,nypu,we,kstrt,tfmove,  &     &tfft,inorder)! calculates potential and electric field! with vacuum boundary conditions and 1D domain decomposition! assumes qt2 has been calculated by cdensity function! input: noff, nyp, nypu, kstrt, inorder! output: we, pot, fxye, tfmove, tfft         implicit none         integer :: kstrt, inorder         real :: we, tfmove         real, dimension(2) :: tfft         real, dimension(:,:), pointer :: pot         real, dimension(:,:,:), pointer :: fxye         integer, dimension(:), pointer :: noff, nyp, nypu! local data         integer :: isign, ierr! return if error         if (lcsolver /= 3) return! solve for potential         call poisc(qt2,pott2,ffg,we,nx,ny,kstrt)         isign = 1         call fft(q2,pott2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2&     &,LINEAR)         call hafdbl(pot,q2,nx,ny,kstrt,kyp,kyp2,inorder)! copy to guard cells in y         call pnlcguard(pot,nypu,kstrt,lvp,nx,nterg,inorder)! copy to guard cells in x, only for quadratic interpolation         call lcguard(pot,nypu,nx,inorder)! solve for force/charge in fourier space         call poisc(qt2,fxyt2,ffg,we,nx,ny,kstrt)         call fft(fxy2,fxyt2,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2,   &     &LINEAR)         call hafdbl(fxye,fxy2,nx,ny,kstrt,kyp,kyp2,inorder)! move force/charge to non-uniform field partition         isign = 1         call pfmove(fxye,noff,nyp,isign,tfmove,kyp,kstrt,lvp,idps,nterf&     &,ierr,inorder)         if (ierr /= 0) then!           call MP_END            call PPEXIT            stop         endif! copy to guard cells in y         call pnlcguard(fxye,nyp,kstrt,lvp,nx,nterg,inorder)! copy to guard cells in x, only for quadratic interpolation         call lcguard(fxye,nyp,nx,inorder)         end subroutine cefieldpot!         subroutine cpotential(pot,nypu,we,kstrt,tfft,inorder)! calculates potential with vacuum boundary conditions! and 1D domain decomposition! assumes qt2 has been calculated by cdensity function! input: nypu, kstrt, inorder! output: we, pot, tfft         implicit none         integer :: kstrt, inorder         real :: we         real, dimension(2) :: tfft         real, dimension(:,:), pointer :: pot         integer, dimension(:), pointer :: nypu! local data         integer :: isign! return if error         if ((lcsolver /= 1) .and. (lcsolver /= 3)) return! solve for potential         call poisc(qt2,pott2,ffg,we,nx,ny,kstrt)         isign = 1         call fft(q2,pott2,isign,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2&     &,LINEAR)         call hafdbl(pot,q2,nx,ny,kstrt,kyp,kyp2,inorder)! copy to guard cells in y         call pnlcguard(pot,nypu,kstrt,lvp,nx,nterg,inorder)! copy to guard cells in x, only for quadratic interpolation         call lcguard(pot,nypu,nx,inorder)         end subroutine cpotential!         subroutine cefield(fxye,noff,nyp,we,kstrt,tfmove,tfft,inorder)! calculates electric field with vacuum boundary conditions! and 1D domain decomposition! assumes qt2 has been calculated by cdensity function! input: noff, nyp, kstrt, inorder! output: we, fxye, tfmove, tfft         implicit none         integer :: kstrt, inorder         real :: we, tfmove         real, dimension(2) :: tfft         real, dimension(:,:,:), pointer :: fxye         integer, dimension(:), pointer :: noff, nyp! local data         integer :: isign, ierr! return if error         if ((lcsolver /= 2) .and. (lcsolver /= 3)) return! solve for force/charge in fourier space         call poisc(qt2,fxyt2,ffg,we,nx,ny,kstrt)         call fft(fxy2,fxyt2,mixup2,sct2,tfft,indx1,indy1,kstrt,kyp2,   &     &LINEAR)         call hafdbl(fxye,fxy2,nx,ny,kstrt,kyp,kyp2,inorder)! move force/charge to non-uniform field partition         isign = 1         call pfmove(fxye,noff,nyp,isign,tfmove,kyp,kstrt,lvp,idps,nterf&     &,ierr,inorder)         if (ierr /= 0) then!           call MP_END            call PPEXIT            stop         endif! copy to guard cells in y         call pnlcguard(fxye,nyp,kstrt,lvp,nx,nterg,inorder)! copy to guard cells in x, only for quadratic interpolation         call lcguard(fxye,nyp,nx,inorder)         end subroutine cefield!      end module esopenpic2_facade