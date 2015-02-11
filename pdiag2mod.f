!-----------------------------------------------------------------------
!
      module pdiag2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pdiag2lib.f
! pdiag2mod.f contains diagnostic procedures:
!             defines module pdiag2d
! bfopen => bfopen2 opens binary file for real 2d scalar data, with 1d
!           domain decomposition.
! bfopen => bfcopen2 opens binary file for complex 2d scalar data, with
!           1d domain decomposition.
! bfopen => bfvcopen2 opens binary file for complex 2d vector data, with
!           1d domain decomposition.
! vdist => ipvdist2 calculates 2 or 3 component velocity distribution
!          and velocity moments.
!          calls PVDIST2 or PVDIST23
! sdist => ipsdist2 calculates entropy from 2 or 3 component velocity
!          distribution.
!          calls PSDIST2 or PSDIST23
! grasp => ipgrasp23 displays particle in phase space.
!          calls PGRASP23
! displayfv => ipdisplayv2 displays velocity distribution functions.
!              calls PDISPR
! displays => ipdscaler2 displays 2d scalar field in real space.
!             calls PCARPET, or PCONTUR
! displays => ipdscaler1 displays 1d scalar field in real space.
!             calls PDISPS
! displayv => ipdvector2 displays amplitude of a 2d vector field in real
!             space.
!             calls PCARPET, or PCONTUR
! displayw => ipdisplayw displays time history of electric field,
!             kinetic and total energies.
!             calls PDISPR
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 14, 2009
!
      use globals, only: LINEAR, QUADRATIC
      use pinit2d, only: idrun, indx, indy, ntp, ntd, nta, ntj, nte,    &
     &psolve, tend, dt, ndim, omx, omy, omz, ci, t0, ceng, indian,      &
     &rlprec, inprec, pden2d, modesxd, modesyd, ndrec, fdname, ppot2d,  &
     &modesxp, modesyp, nprec, fpname, pvpot2d, modesxa, modesya, narec,&
     &faname, pvcur2d, modesxj, modesyj, njrec, fjname, pem2d, modesxe, &
     &modesye, nerec, fename
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: GROPEN, SETNPLT, STPALIT, GRCLOSE, PGRCLOSE
      public :: idrun, indx, indy, ntp, ntd, nta, nte, psolve
      public :: tend, dt, omx, omy, omz, ci, t0, ceng, indian, rlprec
      public :: inprec
      public :: pden2d, modesxd, modesyd, ndrec, fdname
      public :: ppot2d, modesxp, modesyp, nprec, fpname
      public :: pvpot2d, modesxa, modesya, narec, faname
      public :: pem2d, modesxe, modesye, nerec, fename
      public :: grasp, vdist, sdist, displayfv, displays, displayv
      public :: displayw, bfopen
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GROPEN
         implicit none
         end subroutine
      end interface
      interface
         subroutine SETNPLT(nplt,irc)
         implicit none
         integer :: nplt, irc
         end subroutine
      end interface
      interface
         subroutine STPALIT(idpal)
         implicit none
         integer :: idpal
         end subroutine
      end interface
      interface
         subroutine GRCLOSE
         implicit none
         end subroutine
      end interface
      interface
         subroutine PGRCLOSE
         implicit none
         end subroutine
      end interface
      interface
         subroutine PGRASP23(part,f,npp,label,itime,isc,nx,ny,iyp,ixp,id&
     &imp,npmax,nblok,irc)
         implicit none
         integer :: itime, isc, nx, ny, iyp, ixp, idimp, npmax, nblok
         integer :: irc
         character(len=*) :: label
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(2*npmax) :: f
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PVDIST2(part,fv,fvm,npp,idimp,npmax,nblok,nmv,nmvf)
         implicit none
         integer :: idimp, npmax, nblok, nmv, nmvf
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nmvf,2,nblok) :: fv
         real, dimension(2,2,nblok) :: fvm
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PVDIST23(part,fv,fvm,npp,idimp,npmax,nblok,nmv,nmvf)
         implicit none
         integer :: idimp, npmax, nblok, nmv, nmvf
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nmvf,3,nblok) :: fv
         real, dimension(2,3,nblok) :: fvm
         integer, dimension(nblok) :: npp
         end subroutine
      end interface
      interface
         subroutine PSDIST2(fv,fvm,nblok,nmv,nmvf)
         implicit none
         integer :: nblok, nmv, nmvf
         real, dimension(nmvf,2,nblok) :: fv
         real, dimension(2,2,nblok) :: fvm
         end subroutine
      end interface
      interface
         subroutine PSDIST23(fv,fvm,nblok,nmv,nmvf)
         implicit none
         integer :: nblok, nmv, nmvf
         real, dimension(nmvf,3,nblok) :: fv
         real, dimension(2,3,nblok) :: fvm
         end subroutine
      end interface
      interface
         subroutine PDISPS(f,g,nvp,label,xmin,xmax,isc,ist,nx,nxv,nblok,&
     &chr,irc)
         implicit none
         integer :: nvp, isc, ist, nx, nxv, nblok, irc
         real :: xmin, xmax
         character(len=*) :: label, chr
!        real, dimension(*) :: f
         real :: f
         real, dimension(nxv) :: g
         end subroutine
      end interface
      interface
         subroutine DISPR(f,label,xmin,xmax,isc,ist,mks,nx,nxv,ngs,chr,c&
     &hrs,irc)
         implicit none
         integer :: isc, ist, mks, nx, nxv, ngs, irc
         real :: xmin, xmax
         character(len=*) :: label, chr
         character(len=*), dimension(ngs) :: chrs
!        real, dimension(*) :: f
         real :: f
         end subroutine
      end interface
      interface
         subroutine PCARPET(f,g,nvp,label,isc,ist,nx,ny,nxv,nypmx,nblok,&
     &chr,ntc,irc)
         implicit none
         integer :: nvp, isc, ist, nx, ny, nxv, nypmx, nblok, ntc, irc
         character(len=*) :: label, chr
!        real, dimension(*) :: f
         real :: f
         real, dimension(nxv,nypmx) :: g
         end subroutine
      end interface
      interface
         subroutine PCONTUR(f,g,lf,nvp,label,isc,ist,nx,ny,nxv,nypmx,nbl&
     &ok,chr,nc,irc)
         implicit none
         integer :: nvp, isc, ist, nx, ny, nxv, nypmx, nblok, nc, irc
         character(len=*) label, chr
!        real, dimension(*) :: f
         real :: f
         real, dimension(nxv,nypmx) :: g
         integer, dimension(nxv,nypmx,nblok) :: lf
         end subroutine
      end interface
      interface
         subroutine PDISPR(f,g,nvp,label,xmin,xmax,isc,ist,mks,nx,nxv,nb&
     &lok,ngs,chr,chrs,irc)
         integer :: nvp, isc, ist, mks, nx, nxv, nblok, ngs, irc
         real :: xmin, xmax
         character(len=*) :: label, chr
!        real, dimension(*) :: f
         real :: f
         real, dimension(nxv,ngs) :: g
         character(len=*), dimension(ngs) :: chrs
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface grasp
         module procedure ipgrasp23
      end interface
!
      interface vdist
         module procedure ipvdist2
      end interface
!
      interface sdist
         module procedure ipsdist2
      end interface
!
      interface displayfv
         module procedure ipdisplayv2
      end interface
!
      interface displays
         module procedure ipdscaler2
         module procedure ipdscaler1
      end interface
!
      interface displayv
         module procedure ipdvector2
      end interface
!
      interface displayw
         module procedure ipdisplayw
      end interface
!
      interface bfopen
         module procedure bfopen2
         module procedure bfcopen2
         module procedure bfvcopen2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipgrasp23(part,npp,label,itime,isc,nx,ny,iyp,ixp,irc&
     &)
! displays phase space
         implicit none
         integer :: itime, isc, nx, ny, iyp, ixp, irc
         character(len=*) :: label
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         real, dimension(2*size(part,2)) :: f
         integer :: idimp, npmax, nblok
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PGRASP23(part,f,npp,label,itime,isc,nx,ny,iyp,ixp,idimp,np&
     &max,nblok,irc)
         end subroutine ipgrasp23
!
         subroutine ipvdist2(part,fv,fvm,npp,nmv,idimv)
! calculates 2d velocity distribution and velocity moments
! idimv = (2,3) = calculate distribution for (vx,vy) or (vx,vy,vz)
         integer :: nmv, idimv
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fv, fvm
         integer, dimension(:), pointer :: npp
! local data
         integer :: idimp, npmax, nblok, nmvf
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); nmvf = size(fv,1)
         if (idimv==2) then
            call PVDIST2(part,fv,fvm,npp,idimp,npmax,nblok,nmv,nmvf)
         else if (idimv==3) then
            call PVDIST23(part,fv,fvm,npp,idimp,npmax,nblok,nmv,nmvf)
         endif
         end subroutine ipvdist2
!
         subroutine ipsdist2(fv,fvm,nmv,idimv)
! calculates entropy from 2d velocity distribution
! idimv = (2,3) = calculate distribution for (vx,vy) or (vx,vy,vz)
         integer :: nmv, idimv
         real, dimension(:,:,:), pointer :: fv, fvm
! local data
         integer :: nblok, nmvf
         nblok = size(fv,3); nmvf = size(fv,1)
         if (idimv==2) then
            call PSDIST2(fv,fvm,nblok,nmv,nmvf)
         else if (idimv==3) then
            call PSDIST23(fv,fvm,nblok,nmv,nmvf)
         endif
         end subroutine ipsdist2
!
         subroutine ipdisplayv2(fv,fvm,label,itime,nmv,idt,idimv,irc)
! displays velocity distribution functions
! fv = velocity distribution
! fvm = velocity moments
! itime = current time step
! nmv = number of velocity intervals
! idt = (1,2,3) = display (individual,composite,both) functions
! idimv = (2,3) = display (vx,vy) or (vx,vy,vz)
! irc = return code (0 = normal return)
         implicit none
         integer :: itime, nmv, idt, idimv, irc
         real, dimension(:,:,:), pointer :: fv, fvm
         character(len=*) :: label
! isc = 999 = use default display scale
! ist = 1 = display positive values only
! mks = 0 = cycle through line styles
         integer :: isc = 999, ist = 1, mks = 0
         integer :: i, nmvf, nmv2, nblok
         real :: vmax, vmin
         real, dimension(size(fv,1),size(fv,2)) :: g
         character(len=12) :: c
         character(len=2) :: cs
         character(len=54) :: lbl
         character(len=45) :: chr
         character(len=10), dimension(3) :: chrs
   91    format(', T =',i7)
   92    format(' VD =',f9.6,' VTH =',f9.5)
   93    format(' VTX =',f9.5,' VTY =',f9.5)
   94    format(' VTX =',f9.5,' VTY =',f9.5,' VTZ =',f9.5)
! chrs = short array of characters to label individual line samples
         data chrs /'    VX    ','    VY    ','    VZ    '/
         if ((idimv /= 2) .and. (idimv /= 3)) return
         nmvf = size(fv,1); nblok = size(fv,3)
         nmv2 = 2*nmv
         write (c,91) itime
! each velocity distributions on its own plot
         if (idt /= 2) then
            do i = 1, idimv
            cs = trim(adjustl(chrs(i)))
            lbl = trim(label)//' VELOCITY DISTRIBUTION VERSUS '//cs//c
            write (chr,92) fvm(1,i,1), fvm(2,i,1)
            vmax = fv(1,i,1)
            vmin = -vmax
            call PDISPR(fv(2,1,1),g,1,lbl,vmin,vmax,isc,ist,mks,nmv2,nmv&
     &f,nblok,1,chr,chrs(i),irc)
            if (irc==1) return
            enddo
         endif
! all velocity distributions on common plot
         if (idt /= 1) then
            lbl = trim(label)//' VELOCITY DISTRIBUTIONS VERSUS '//'V'//c
            if (idimv==2) then
               write (chr,93) fvm(2,1,1), fvm(2,2,1)
               vmax = max(fv(1,1,1),fv(1,2,1))
               vmin = -vmax
            else if (idimv==3) then
               write (chr,94) fvm(2,1,1), fvm(2,2,1), fvm(2,3,1)
               vmax = max(fv(1,1,1),fv(1,2,1),fv(1,3,1))
               vmin = -vmax
            endif
            call PDISPR(fv(2,1,1),g,1,lbl,vmin,vmax,isc,ist,mks,nmv2,nmv&
     &f,nblok,idimv,chr,chrs,irc)
         endif
         end subroutine ipdisplayv2
!
         subroutine ipdscaler1(f,nvp,label,itime,isc,ist,nx,irc,inorder)
! displays 1d scalar field in real space
! f = 1d scalar field in real space
! nvp = number of real or virtual processors requested
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
! ist = flag for choosing positive and/or negative values
! nx = system length in x direction
! irc = return code (0 = normal return)
         implicit none
         integer :: nvp, itime, isc, ist, nx, irc
         integer, optional :: inorder
         character(len=*) label
         real, dimension(:,:), pointer :: f
! local data
         real, dimension(size(f,1)) :: g
         integer :: nxv, nblok, order
         real :: xmin, xmax
         character(len=12) :: lbl
   91    format(' T = ',i7)
         nxv = size(f,1); nblok = size(f,2)
         xmin = 0.0; xmax = real(nx)
         write (lbl,91) itime
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PDISPS(f(1,1),g,nvp,label,xmin,xmax,isc,ist,nx,nxv,nblo&
     &k,lbl,irc)
         else
            nxv = nxv - 1
            call PDISPS(f(2,1),g,nvp,label,xmin,xmax,isc,ist,nx,nxv,nblo&
     &k,lbl,irc)
         endif
         end subroutine ipdscaler1
!
         subroutine ipdscaler2(pot,nvp,label,itime,isc,ist,idt,nx,ny,irc&
     &,inorder)
! displays 2d scalar field in real space
! pot = 2d scalar field in real space
! nvp = number of real or virtual processors requested
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of pot
! ist = flag for choosing positive and/or negative values
! idt = (1,2,3) = display (color map,contour plot,both)
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
         implicit none
         integer :: nvp, itime, isc, ist, idt, nx, ny, irc
         integer, optional :: inorder
         character(len=*) label
         real, dimension(:,:,:), target :: pot
! local data
         real, dimension(size(pot,1),size(pot,2)) :: g
         integer, dimension(size(pot,1),size(pot,2),size(pot,3)) :: lf
! ntc = number of valid colors, should be power of 2, <= 256
! nc = number of contour lines
         real, dimension(:,:,:), pointer :: plt
         integer :: ntc = 16, nc = 16
         integer :: nxv, nypmx, nblok, mx, order
         character(len=12) :: lbl
   91    format(' T = ',i7)
         plt => pot
         nxv = size(pot,1); nypmx = size(pot,2); nblok = size(pot,3)
         mx = nx
! plot guard cells if present
         if ((mx+1) <= nxv) mx = mx + 1
         write (lbl,91) itime
         order = QUADRATIC
         if (present(inorder)) order = inorder
! color map plot for all values
         if (idt /= 2) then
            if (order==LINEAR) then
               call PCARPET(plt(1,1,1),g,nvp,label,isc,ist,mx,ny,nxv,nyp&
     &mx,nblok,lbl,ntc,irc)
            else
               call PCARPET(plt(2,2,1),g,nvp,label,isc,ist,mx,ny,nxv,nyp&
     &mx,nblok,lbl,ntc,irc)
            endif
         endif
! contour map for all values
         if (idt /= 1) then
            if (order==LINEAR) then
               call PCONTUR(plt(1,1,1),g,lf,nvp,label,isc,ist,mx,ny,nxv,&
     &nypmx,nblok,lbl,nc,irc)
            else
               call PCONTUR(plt(2,2,1),g,lf,nvp,label,isc,ist,mx,ny,nxv,&
     &nypmx,nblok,lbl,nc,irc)
            endif
         endif
         end subroutine ipdscaler2
!
         subroutine ipdvector2(vpot,nvp,label,itime,isc,ist,idt,nx,ny,ir&
     &c,inorder)
! displays 2d vector field in real space
! vpot = 2d vector field in real space
! nvp = number of real or virtual processors requested
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of pot
! ist = flag for choosing positive and/or negative values
! idt = (1,2,3) = display (color map,contour plot,both)
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
         implicit none
         integer :: nvp, itime, isc, ist, idt, nx, ny, irc
         integer, optional :: inorder
         character(len=*) label
         real, dimension(:,:,:,:), pointer :: vpot
! local data
         integer :: i, j, k, l, mx, order
         real :: sum1
         real, dimension(size(vpot,2),size(vpot,3),size(vpot,4)) :: pot
! calculate array size with guard cells
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            mx = nx + 1
         else
            mx = nx + 3
         endif
         mx = min(mx,size(vpot,2))
! find absolute value of vector field
         do l = 1, size(vpot,4)
         do k = 1, size(vpot,3)
         do j = 1, mx
            sum1 = 0.0
            do i = 1, size(vpot,1)
            sum1 = sum1 + vpot(i,j,k,l)**2
            enddo
            pot(j,k,l) = sqrt(sum1)
         enddo
         enddo
         enddo
! display amplitude
         call ipdscaler2(pot,nvp,label,itime,isc,ist,idt,nx,ny,irc,inord&
     &er)
         end subroutine ipdvector2
!
         subroutine ipdisplayw(wt,t0,dtw,nt,irc)
! displays time history of electric field, kinetic, and total energies
! wt = time history array for energies
! t0 = initial energy time
! dtw = time between energy values
! nt = number of energy values to be displayed
! irc = return code (0 = normal return)
         integer :: nt, irc
         real :: t0, dtw
         real, dimension(:,:), pointer :: wt
! isc = 999 = use default display scale
! ist = 2 = display minimum range
! mks = 0 = cycle through line styles
         integer :: isc = 999, ist = 2, mks = 0
         integer :: i, ntwd, ns, nd, nblok = 1
         real :: tmin, tmax
         real, dimension(size(wt,1),size(wt,2)) :: g
         character(len=36) :: lbl
         character(len=20), dimension(7) :: cs 
         character(len=10), dimension(7) :: chrs
! chrs = short array of characters to label individual line samples
         data cs /' TOTAL FIELD        ',' ELECTRON KINETIC   ',' ION KI&
     &NETIC     ',' TOTAL              ',' ES FIELD           ',' ET FIE&
     &LD        ',' MAGNETIC FIELD     '/
         data chrs /'TOT FIELD ','ELECT ENRG',' ION ENRG ','TOTAL ENRG',&
     &' EL FIELD ',' ET FIELD ',' B FIELD  '/
! quit if array is empty or incorrect
         if (nt <= 0) return
         ntwd = size(wt,1)
         ns = min(size(wt,2),7)
         nd = nt - 1
! tmin/tmax = range of time values in plot
         tmin = t0
         tmax = t0 + dtw*(nt - 1)
! display individual energies
         do i = 1, ns
         lbl = trim(cs(i))//' ENERGY VERSUS TIME'
         call PDISPR(wt(1,i),g,1,lbl,tmin,tmax,isc,ist,mks,nd,ntwd,nblok&
     &,1,' ',chrs(i),irc)
         if (irc==1) return
         enddo
! all energies on common plot
         lbl = ' ENERGIES VERSUS TIME'
         call PDISPR(wt(1,1),g,1,lbl,tmin,tmax,isc,ist,mks,nd,ntwd,nblok&
     &,ns,' ',chrs,irc)
         end subroutine ipdisplayw
!
         subroutine bfopen2(f,nx,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for real 2d scalar data with 1d domain decomposition
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nx, iunit, nrec
         real, dimension(:,:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, ierr
         if (nrec > 0) return
         inquire(iolength=lrec) f(1,1,1)
         lrec = nx*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfopen2
!
         subroutine bfcopen2(f,nz,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 2d scalar data with 1d domain decomposition
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nz, iunit, nrec
         complex, dimension(:,:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, ierr
         if (nrec > 0) return
         inquire(iolength=lrec) f(1,1,1)
         lrec = nz*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfcopen2
!
         subroutine bfvcopen2(f,nz,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 2d vector data with 1d domain decomposition
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nz, iunit, nrec
         complex, dimension(:,:,:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, nnz, ierr
         if (nrec > 0) return
         nnz = size(f,1)*nz
         inquire(iolength=lrec) f(1,1,1,1)
         lrec = nnz*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfvcopen2
!
      end module pdiag2d
