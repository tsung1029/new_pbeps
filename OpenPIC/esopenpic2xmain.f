!-----------------------------------------------------------------------! 2D Electrostatic PIC code with Vacuum Poisson solver! with 1D domain decomposition! version compatible with both shared and distributed memory programming! written by viktor k. decyk, ucla! update: april 25, 2008!      program esopenpic2x_main      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR      use esopenpic2x_facade      use p2d, only: wtimer, plsum      implicit none! indx/indy = exponent which determines length in x/y direction,! where nx=2**indx, ny=2**indy! npx/npy = initial number of particles distributed in x/y direction      integer :: indx =   6, indy =   7, npx =     384, npy =     768! inorder = interpolation order (LINEAR = 1, QUADRATIC = 2)      integer :: inorder = QUADRATIC! popt = particle optimization scheme! dopt = charge deposit optimization scheme! djopt = current deposit optimization scheme      integer :: popt = STANDARD, dopt = LOOKAHEAD! dt = time interval between successive calculations      real :: dt = 0.2000000e+00! qme = charge on electron, in units of e! vtx/vty = thermal velocity of electrons in x/y direction! vx0/vy0 = drift velocity of electrons in x/y direction      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0! relativity = (no,yes) = (0,1) = relativity is used      integer :: relativity = 0! ci = reciprical of velocity of light      real :: ci = 1.0! sortime = number of time steps between electron sorting      integer :: sortime = 50! nsrand = (0,1) = (no,yes) randomize spatially positions locally! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,!                        hyperbolic secant squared=4)      integer :: nsrand = 0, ndprof = 0! ampdx/ampdx = amplitude of density compared to uniform in x/y! scaledx/scaledx = scale length for spatial coordinate in x/y! shiftdx/shiftdx = shift of spatial coordinate in x/y      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0      real :: ampdy = 0.0, scaledy = 0.0, shiftdy = 0.0! imbalance = load imbalance fraction repartition trigger! (< 0.  to suppress repartion)      real :: imbalance = .08! csolver = (1,2,3) = solve for (potential,efield,both)      integer :: csolver = 2      integer :: kstrt, nvp, np, itime      real :: we = 0.0, wke = 0.0, pibal = 0.0      real :: tpush = 0.0, tdpost = 0.0, tmove = 0.0, tfmove = 0.0      real :: tsort = 0.0, time = 0.0      real :: qbme, ar, affp      double precision :: dtime      real, dimension(:,:,:), pointer :: part!     real, dimension(:,:,:), pointer :: part2      real, dimension(:,:,:), pointer :: qe, qi, pot      real, dimension(:,:,:,:), pointer :: fxye      real, dimension(:,:), pointer :: edges      integer, dimension(:), pointer :: nyp, noff, nypu, noffu      integer, dimension(:), pointer :: npp!      real, dimension(2) :: tfft = 0.0      real, dimension(3) :: wtot!      call init_parallel(edges,nypu,noffu,indx,indy,kstrt,nvp,inorder)      call allocate_fields(qe,qi,fxye,imbalance,inorder)      np = npx*npy      call allocate_particles(part,npp,np)!! initialize timer      call wtimer(time,dtime,-1)! initialize constants      qbme = qme! ar = half-width of particle in r direction      if (inorder==LINEAR) then         ar = .912871      else         ar = .866667      endif! affp = normalization constant for poisson solver      affp = real((2**indx-2)*(2**indy-2))/real(np)!! initialize particles      call initialize_particles(part,npp,ampdx,scaledx,shiftdx,ampdy,   &     &scaledy,shiftdy,vtx,vty,vx0,vy0,npx,npy,kstrt,ndprof,nsrand,dt)! initialize solver for vacuum boundary conditions      call initialize_csolver(qe,pot,ar,affp,indx,indy,csolver,kstrt,nvp&     &,inorder)! find initial non-uniform partition analytically      call initial_partition(edges,noff,nyp,noffu,nypu,fxye,ampdy,      &     &scaledy,shiftdy,kstrt,ndprof,imbalance,inorder)! move particles into appropriate spatial regions      call particle_manager(part,edges,npp,kstrt,tmove)! initialize background charge density      qi = 0.0      call deposit_charge(part,qi,-qme,tdpost,npp,noff,inorder,dopt)      call add_guardcells(qi,nyp,kstrt,inorder)! record time      call wtimer(time,dtime)      if (kstrt==1) write (*,*) 'initialization time=', time      call wtimer(time,dtime,-1)!! * * * start main iteration loop * * *!      do itime = 0, 324!! initialize charge density      qe = 0.0! deposit charge      call deposit_charge(part,qe,qme,tdpost,npp,noff,inorder,dopt)! add ion charge density      qe = qe + qi!! vacuum boundary conditions!! add guard cells, move to uniform partition, and transform density      call cdensity(qe,noff,nyp,kstrt,tfmove,tfft,inorder)! solve for electric field      call cefield(fxye,noff,nyp,we,kstrt,tfmove,tfft,inorder)!! push particles      wke = 0.      call push_particles(part,fxye,npp,noff,qbme,dt,ci,wke,tpush,      &     &relativity,inorder,popt)! move particles into appropriate spatial regions      call particle_manager(part,edges,npp,kstrt,tmove,pibal)!! repartition      if (imbalance >= 0.0) then         if (pibal > imbalance) then            call repartition(part,qi,npp,edges,noff,nyp,kstrt,tmove,    &     &tfmove,inorder)            if (kstrt==1) then               write (*,*) 'repartitioning complete, imbalance = ',pibal            endif         endif      endif!! sort particles      if (sortime > 0) then         if (mod(itime,sortime)==0) then            call sort_particles(part,npp,noff,nyp,tsort,inorder)            if (kstrt==1) write (*,*) 'particles sorted'         endif      endif! energy diagnostic      wtot(1) = we      wtot(2) = wke      wtot(3) = we + wke      call plsum(wtot)      if (kstrt==1) write (*,*) 'time, energies=', itime, wtot!      enddo!! * * * end main iteration loop * * *!! record timings      call wtimer(time,dtime)      if (kstrt==1) then         write (*,*) 'partition used: nvp = ', nvp         write (*,*) 'total main time=', time, 'sec'         time = tpush + tdpost         write (*,*) 'push and deposit time = ', time, 'sec'         write (*,*) 'sort time = ', tsort         write (*,*) 'particle manager time = ', tmove, 'sec'         time = time + tsort + tmove         write (*,*) 'total particle time = ', time, 'sec'         write (*,*) 'fft, transpose time = ', tfft, 'sec'         write (*,*) 'moving partition time = ', tfmove, 'sec'      endif!     call MP_END      call PPEXIT      stop!      end program