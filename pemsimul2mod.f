!-----------------------------------------------------------------------
!
      module pemsimul2d
! Higher level subroutines for electromagnetics
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: november 17, 2009
      use globals, only: LINEAR, QUADRATIC
      use pdiag2d, only: displayv, bfopen
      use psimul2d
      use pempush2d, only: retard, djpost, rdjpost, push3, rpush3,      &
     &pushzf, rpush3zf, dmjpost, rdmjpost, dcjpost, rdcjpost,           &
     &premoment2, fft, get_funit, plsum, writebf, paguard, pcguard
      use pfield2d, only: aguard, cguard, cuperp, sbpois, apois, avpot, &
     &avrpot, gtmodes, poynt
      implicit none
      private
      public :: restart_open, restart_bwrite, restart_dwrite
      public :: restart_bread, restart_dread
      public :: dpostg, pushg, bpushg, pushzfg
      public :: initmodediag, initveldiag
      public :: dendiag, phasediag, veldiag, potdiag
      public :: emomtdiag, imomtdiag, esenergy
      public :: retardg, djpostg, push3g, push3zfg, dmjpostg, dcjpostg
      public :: initvmodediag, vpotdiagprep, vcurdiagprep, vpotrdiagprep
      public :: vpotdiag, vcurdiag, avpotdiag, avcurdiag, vpotrdiag
      public :: fmomtdiag, dmenergy, emenergy
!
      contains
!
         subroutine retardg(part,npp,dtc,ci,nx,ny,ipbc,relativity,ndim)
! retards particle positions half time-step to deposit current
         implicit none
         integer :: nx, ny, ipbc, relativity
         integer, optional :: ndim
         real :: dtc, ci
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
         if (dtc==0.0) return
         if (relativity==1) then
            call retard(part,npp,dtc,ci,nx,ny,ipbc,ndim)
         else
            call retard(part,npp,dtc,nx,ny,ipbc,ndim)
         endif
         end subroutine retardg
!
         subroutine djpostg(part,cu,npp,noff,qm,dt,ci,tdjpost,nx,ny,ipbc&
     &,relativity,inorder,djopt)
! deposit current
         implicit none
         integer :: nx, ny, ipbc, relativity
         integer, optional :: inorder, djopt
         real :: qm, dt, ci, tdjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: npp, noff
! deposit current
         if (relativity==1) then
            call rdjpost(part,cu,npp,noff,qm,dt,ci,tdjpost,nx,ny,ipbc,in&
     &order,djopt)
         else
            call djpost(part,cu,npp,noff,qm,dt,tdjpost,nx,ny,ipbc,inorde&
     &r,djopt)
         endif
         end subroutine djpostg
!
         subroutine push3g(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,tpush,&
     &nx,ny,ipbc,relativity,inorder,popt)
! push particles with 2-1/2d electromagnetic fields
         implicit none
         integer :: nx, ny, ipbc, relativity
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
         if (relativity==1) then
            call rpush3(part,fxy,bxy,npp,noff,qbm,dt,dtc,ci,ek,tpush,nx,&
     &ny,ipbc,inorder,popt)
         else
            call push3(part,fxy,bxy,npp,noff,qbm,dt,dtc,ek,tpush,nx,ny,i&
     &pbc,inorder,popt)
         endif
         end subroutine push3g
!
         subroutine push3zfg(part,npp,dt,ci,ek,tpush,nx,ny,ipbc,ndim,rel&
     &ativity)
! push particles with no forces
         implicit none
         integer :: nx, ny, ipbc, ndim, relativity
         real :: dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! push particles
         if (relativity==1) then
            call rpush3zf(part,npp,dt,ci,ek,tpush,nx,ny,ipbc,ndim)
         else
            call pushzf(part,npp,dt,ek,tpush,nx,ny,ipbc)
         endif
         end subroutine push3zfg
!
         subroutine dmjpostg(part,amu,npp,noff,qm,ci,tdcjpost,relativity&
     &,inorder,djopt)
! deposit momentum flux with 2-1/2d electromagnetic fields
         implicit none
         integer :: relativity
         integer, optional :: inorder, djopt
         real :: qm, ci, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: amu
         integer, dimension(:), pointer :: npp, noff
! deposit momentum flux
         if (relativity==1) then
            call rdmjpost(part,amu,npp,noff,qm,ci,tdcjpost,inorder,djopt&
     &)
         else
            call dmjpost(part,amu,npp,noff,qm,tdcjpost,inorder,djopt)
         endif
         end subroutine dmjpostg
!
         subroutine dcjpostg(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,dt,&
     &ci,tdcjpost,relativity,inorder,djopt)
! deposit momentum flux, acceleration density, and current density
! with 2-1/2d electromagnetic fields
         implicit none
         integer :: relativity
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy, cu, dcu, amu
         integer, dimension(:), pointer :: npp, noff
         if (relativity==1) then
            call rdcjpost(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,dt,ci,&
     &tdcjpost,inorder,djopt)
         else
            call dcjpost(part,fxy,bxy,npp,noff,cu,dcu,amu,qm,qbm,dt,tdcj&
     &post,inorder,djopt)
         endif
         end subroutine dcjpostg
!
         subroutine initvmodediag(vcurt,ntj,id0,nxh,nyh,kxp,ndim,modesxj&
     &,modesyj,jblok,iuj,njrec,fjname)
! initialize vector mode diagnostic
         implicit none
         integer :: ntj, id0, nxh, nyh, kxp, ndim
         integer :: modesxj, modesyj, jblok, iuj, njrec
         character(len=*) :: fjname
         complex, dimension(:,:,:,:), pointer :: vcurt
! local data
         integer :: modesy2j
         if (ntj <= 0) return
         if (modesxj > nxh) modesxj = nxh
         if (modesyj > nyh) modesyj = nyh
         modesy2j = 2*modesyj - 1
         allocate(vcurt(ndim,modesy2j,min(modesxj,kxp),jblok))
! open output file
         if (id0==0) then
            if (njrec==0) then
               iuj = get_funit(iuj); njrec = -1
               call bfopen(vcurt,modesy2j,iuj,njrec,trim(fjname))
            endif
         else
            if (njrec==0) njrec = 1
         endif
         end subroutine initvmodediag
!
         subroutine vpotdiagprep(cut,vfieldt,nta,nda,ntime)
! save data for vector potential diagnostic
         implicit none
         integer :: nta, nda, ntime
         complex, dimension(:,:,:,:), pointer :: cut, vfieldt
! local data
         integer :: it, jt
         if ((nta > 0) .or. (nda > 0)) then
            it = -1; if (nta > 0) it = ntime - nta*(ntime/nta)
            jt = -1; if (nda > 0) jt = ntime - nda*(ntime/nda)
! save old current
            if ((it==0) .or. (jt==0)) then
               vfieldt = cut
            endif
         endif
         end subroutine vpotdiagprep
!
         subroutine vcurdiagprep(cu,vfield,ntj,ndj,ntime)
! save data for static (darwin) ion current diagnostic
         implicit none
         integer :: ntj, ndj, ntime
         real, dimension(:,:,:,:), pointer :: cu, vfield
! local data
         integer :: it, jt
         if ((ntj > 0) .or. (ndj > 0)) then
            it = -1; if (ntj > 0) it = ntime - ntj*(ntime/ntj)
            jt = -1; if (ndj > 0) jt = ntime - ndj*(ntime/ndj)
! save old current
            if ((it==0) .or. (jt==0)) then
               vfield = cu
            endif
         endif
         end subroutine vcurdiagprep
!
         subroutine vpotrdiagprep(cu,vfield,nta,nda,ntime,diff)
! save data for electromagnetic diagnostic
         implicit none
         integer :: nta, nda, ntime
         real, dimension(:,:,:,:), pointer :: cu, vfield
         logical, optional :: diff
! local data
         integer :: is, js, it, jt
         logical :: ldiff
         if ((nta > 0) .or. (nda > 0)) then
            ldiff = .false.
            if (present(diff)) ldiff = diff
            is = -1; if (nta > 0) is = ntime - nta*(ntime/nta)
            js = -1; if (nda > 0) js = ntime - nda*(ntime/nda)
            it = -1; if (nta > 1) it = ntime - nta*((ntime-1)/nta) - 1
            jt = -1; if (nda > 1) jt = ntime - nda*((ntime-1)/nda) - 1
! save current
            if ((is==0) .or. (js==0) .or. (it==0) .or. (jt==0)) then
               if (ldiff) then
                  vfield = cu - vfield
               else
                  vfield = cu
               endif
            endif
         endif
         end subroutine vpotrdiagprep
!
         subroutine vpotdiag(cut,vfield,vpott,vfieldt,ffc,nyp,mixup,sct,&
     &ci,tfft,nta,nda,nx,ny,modesxa,modesya,iua,narec,indx,indy,ntime,nv&
     &p,kstrt,kxp,kyp,ndstyle,irc,inorder)
! static (darwin) vector potential diagnostic
         implicit none
         integer :: nta, nda, nx, ny, modesxa, modesya, iua, narec
         integer :: indx, indy, ntime, nvp, kstrt, kxp, kyp
         integer :: ndstyle, irc
         real :: ci
         real, dimension(2) :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: vfield
         complex, dimension(:,:,:,:), pointer :: cut, vpott, vfieldt
         complex, dimension(:,:,:), pointer :: ffc
         integer, dimension(:), pointer :: nyp
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2a, isign
         real :: wm
         irc = 0
         if ((nta > 0) .or. (nda > 0)) then
            it = -1; if (nta > 0) it = ntime - nta*(ntime/nta)
            jt = -1; if (nda > 0) jt = ntime - nda*(ntime/nda)
            if ((it==0) .or. (jt==0)) then
! calculate vector potential in fourier space
               call apois(cut,vfieldt,ffc,ci,wm,nx,ny,kstrt)
! store selected fourier modes
               if (it==0) then
                  modesy2a = 2*modesya - 1
                  call gtmodes(vfieldt,vpott,nx,ny,modesxa,modesya,kstrt&
     &)
! write diagnostic output
                  call writebf(vpott,modesxa,modesy2a,kxp,iua,narec)
               endif
! transform vector potential to real space
               if (jt==0) then
                  isign = 1
                  call fft(vfield,vfieldt,isign,mixup,sct,tfft,indx,indy&
     &,kstrt,kyp,inorder)
                  call pcguard(vfield,kstrt,nvp,kyp,inorder)
                  call cguard(vfield,nyp,nx,inorder)
! display absolute value of vector potential
                  call displayv(vfield,nvp,' VECTOR POTENTIAL',ntime,999&
     &,1,ndstyle,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine vpotdiag
!
         subroutine vcurdiag(cu,cut,vfield,vcurt,vfieldt,ffc,nyp,mixup,s&
     &ct,tfft,ntj,ndj,nx,ny,modesxj,modesyj,iuj,njrec,indx,indy,ntime,nv&
     &p,kstrt,kxp,kyp,ngds,ndstyle,irc,inorder)
! static (darwin) ion current diagnostic
         implicit none
         integer :: ntj, ndj, nx, ny, modesxj, modesyj, iuj, njrec
         integer :: indx, indy, ntime, nvp, kstrt, kxp, kyp, ngds
         integer :: ndstyle, irc
         real, dimension(2) :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: cu, vfield
         complex, dimension(:,:,:,:), pointer :: cut, vcurt, vfieldt
         complex, dimension(:,:,:), pointer :: ffc
         integer, dimension(:), pointer :: nyp
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2j, isign
         irc = 0
         if ((ntj > 0) .or. (ndj > 0)) then
            it = -1; if (ntj > 0) it = ntime - ntj*(ntime/ntj)
            jt = -1; if (ndj > 0) jt = ntime - ndj*(ntime/ndj)
            if ((it==0) .or. (jt==0)) then
               vfield = cu - vfield
! add guard cells for current
               call aguard(vfield,nyp,nx,inorder)
               call paguard(vfield,kstrt,nvp,nx,kyp,ngds)
! transform ion current to fourier space
               call fft(vfield,cut,isign,mixup,sct,tfft,indx,indy,kstrt,&
     &kyp,inorder)
! take transverse part of current
               call cuperp(cut,nx,ny,kstrt) 
! calculate smoothing in fourier space
               call sbpois(cut,vfieldt,ffc,nx,ny,kstrt)
! store selected fourier modes
               if (it==0) then
                  modesy2j = 2*modesyj - 1
                  call gtmodes(vfieldt,vcurt,nx,ny,modesxj,modesyj,kstrt&
     &)
! write diagnostic output
                  call writebf(vcurt,modesxj,modesy2j,kxp,iuj,njrec)
               endif
! transform ion current to real space
               if (jt==0) then
                  isign = 1
                  call fft(vfield,vfieldt,isign,mixup,sct,tfft,indx,indy&
     &,kstrt,kyp,inorder)
                  call pcguard(vfield,kstrt,nvp,kyp,inorder)
                  call cguard(vfield,nyp,nx,inorder)
! display absolute value of ion current
                  call displayv(vfield,nvp,' ION CURRENT',ntime,999,1,  &
     &ndstyle,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine vcurdiag
!
         subroutine avpotdiag(bxyz,vfield,vpott,vfieldt,nyp,mixup,sct,tf&
     &ft,nta,nda,nx,ny,modesxa,modesya,iua,narec,indx,indy,ntime,nvp,kst&
     &rt,kxp,kyp,ndstyle,irc,inorder)
! vector potential diagnostic
         implicit none
         integer :: nta, nda, nx, ny, modesxa, modesya, iua, narec
         integer :: indx, indy, ntime, nvp, kstrt, kxp, kyp
         integer :: ndstyle, irc
         real, dimension(2) :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: vfield
         complex, dimension(:,:,:,:), pointer :: bxyz, vpott, vfieldt
         integer, dimension(:), pointer :: nyp
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2a, isign
         irc = 0
         if ((nta > 0) .or. (nda > 0)) then
            it = -1; if (nta > 0) it = ntime - nta*(ntime/nta)
            jt = -1; if (nda > 0) jt = ntime - nda*(ntime/nda)
            if ((it==0) .or. (jt==0)) then
! calculate vector potential in fourier space
               call avpot(bxyz,vfieldt,nx,ny,kstrt)
! store selected fourier modes
               if (it==0) then
                  modesy2a = 2*modesya - 1
                  call gtmodes(vfieldt,vpott,nx,ny,modesxa,modesya,kstrt&
     &)
! write diagnostic output
                  call writebf(vpott,modesxa,modesy2a,kxp,iua,narec)
               endif
! transform vector potential to real space
               if (jt==0) then
                  isign = 1
                  call fft(vfield,vfieldt,isign,mixup,sct,tfft,indx,indy&
     &,kstrt,kyp,inorder)
                  call pcguard(vfield,kstrt,nvp,kyp,inorder)
                  call cguard(vfield,nyp,nx,inorder)
! display absolute value of vector potential
                  call displayv(vfield,nvp,' VECTOR POTENTIAL',ntime,999&
     &,1,ndstyle,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine avpotdiag
!
         subroutine avcurdiag(cui,cu,cut,vfield,vcurt,vfieldt,ffc,nyp,mi&
     &xup,sct,tfft,ntj,ndj,nx,ny,modesxj,modesyj,iuj,njrec,indx,indy,nti&
     &me,nvp,kstrt,kxp,kyp,ngds,ndstyle,irc,inorder)
! ion current diagnostic for electromagnetic code
         implicit none
         integer :: ntj, ndj, nx, ny, modesxj, modesyj, iuj, njrec
         integer :: indx, indy, ntime, nvp, kstrt, kxp, kyp, ngds
         integer :: ndstyle, irc
         real, dimension(2) :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: cui, cu, vfield
         complex, dimension(:,:,:,:), pointer :: cut, vcurt, vfieldt
         complex, dimension(:,:,:), pointer :: ffc
         integer, dimension(:), pointer :: nyp
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: is, js, it, jt, modesy2j, isign
         irc = 0
         if ((ntj > 0) .or. (ndj > 0)) then
            is = -1; if (ntj > 0) is = ntime - ntj*((ntime)/ntj)
            js = -1; if (ndj > 0) js = ntime - ndj*((ntime)/ndj)
            it = -1; if (ntj > 0) it = ntime - ntj*((ntime-1)/ntj) - 1
            jt = -1; if (ndj > 0) jt = ntime - ndj*((ntime-1)/ndj) - 1
! save current if needed next time and not saved below
            if (((is==0) .or. (js==0)).and.((it/=0) .and. (jt/=0))) then
               cu = cui
            endif
            if ((it==0) .or. (jt==0)) then
! calculate averaged ion current
               vfield = 0.5*(cui + cu)
! save current if needed next time and not saved above
               if ((is==0) .or. (js==0)) then
                  cu = cui
               endif
! add guard cells for current
               call aguard(vfield,nyp,nx,inorder)
               call paguard(vfield,kstrt,nvp,nx,kyp,ngds)
! transform ion current to fourier space
               isign = -1
               call fft(vfield,cut,isign,mixup,sct,tfft,indx,indy,kstrt,&
     &kyp,inorder)
! take transverse part of current
               call cuperp(cut,nx,ny,kstrt) 
! calculate smoothing in fourier space
               call sbpois(cut,vfieldt,ffc,nx,ny,kstrt)
! store selected fourier modes
               if (it==0) then
                  modesy2j = 2*modesyj - 1
                  call gtmodes(vfieldt,vcurt,nx,ny,modesxj,modesyj,kstrt&
     &)
! write diagnostic output
                  call writebf(vcurt,modesxj,modesy2j,kxp,iuj,njrec)
               endif
! transform ion current to real space
               if (jt==0) then
                  isign = 1
                  call fft(vfield,vfieldt,isign,mixup,sct,tfft,indx,indy&
     &,kstrt,kyp,inorder)
                  call pcguard(vfield,kstrt,nvp,kyp,inorder)
                  call cguard(vfield,nyp,nx,inorder)
! display absolute value of ion current
                  call displayv(vfield,nvp,' ION CURRENT',ntime,999,1,  &
     &ndstyle,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine avcurdiag
!
         subroutine vpotrdiag(bxyz,cut,vfield,vpotr,vfieldt,ffc,nyp,mixu&
     &p,sct,tfft,affp,ci,nte,nde,nx,ny,modesxe,modesye,iue,nerec,indx,in&
     &dy,ntime,nvp,kstrt,kxp,kyp,ndstyle,irc,inorder)
! electromagnetic diagnostic
         implicit none
         integer :: nte, nde, nx, ny, modesxe, modesye, iue, nerec
         integer :: indx, indy, ntime, nvp, kstrt, kxp, kyp
         integer :: ndstyle, irc
         real :: affp, ci
         real, dimension(2) :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: vfield
         complex, dimension(:,:,:,:), pointer :: bxyz, cut, vpotr
         complex, dimension(:,:,:,:), pointer :: vfieldt
         complex, dimension(:,:,:), pointer :: ffc
         integer, dimension(:), pointer :: nyp
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2e, isign
         irc = 0
         if ((nte > 0) .or. (nde > 0)) then
            it = -1; if (nte > 0) it = ntime - nte*((ntime-1)/nte) - 1
            jt = -1; if (nde > 0) jt = ntime - nde*((ntime-1)/nde) - 1
            if ((it==0) .or. (jt==0)) then
! calculate averaged radiative vector potential
               vfieldt = 0.5*(vfieldt + cut)
               call avrpot(vfieldt,bxyz,ffc,affp,ci,nx,ny,kstrt)
! store selected fourier modes
               if (it==0) then
                  modesy2e = 2*modesye - 1
                  call gtmodes(vfieldt,vpotr,nx,ny,modesxe,modesye,kstrt&
     &)
! write diagnostic output
                  call writebf(vpotr,modesxe,modesy2e,kxp,iue,nerec)
               endif
! transform radiative vector potential to real space
               if (jt==0) then
                  isign = 1
                  call fft(vfield,vfieldt,isign,mixup,sct,tfft,indx,indy&
     &,kstrt,kyp,inorder)
                  call pcguard(vfield,kstrt,nvp,kyp,inorder)
                  call cguard(vfield,nyp,nx,inorder)
! display absolute value of radiative vector potential
                  call displayv(vfield,nvp,' RADIATIVE VPOTENTIAL',ntime&
     &,999,1,ndstyle,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine vpotrdiag
!
         subroutine fmomtdiag(part,q,ffc,exyz,bxyz,npp,msg,affp,pxe,pye,&
     &pze,sx,sy,sz,wx,wy,wz,ntm,id0,ium,kstrt,nx,ny,ntime)
! calculate electron and field momentum
         implicit none
         integer :: ntm, id0, ium, kstrt, nx, ny, ntime
         real :: affp, pxe, pye, pze, sx, sy, sz, wx, wy, wz
         real, dimension(:,:,:), pointer :: part
         complex, dimension(:,:,:), pointer :: q
         complex, dimension(:,:,:,:), pointer :: exyz, bxyz
         complex, dimension(:,:,:), pointer :: ffc
         integer, dimension(:), pointer :: npp
         double precision, dimension(:) :: msg
! local data
         integer :: it
         if (ntm > 0) then
            it = ntime/ntm
! calculate the momentum in the electromagnetic field
            if (ntime==ntm*it) then
               call poynt(q,exyz,bxyz,ffc,affp,sx,sy,sz,nx,ny,kstrt)
            endif
! calculate the momentum in the electrons
            it = ntime - ntm*it + 1
            if (it > 1) it = it - ntm
            if (it >= 0) then
               call premoment2(part,ntime,npp,msg,id0,ium,pxe,pye,pze,sx&
     &,sy,sz,wx,wy,wz,nprint=it)
            endif
! send momentum values to diagnostic node
            if (it==1) call HARTBEAT(msg,6)
         endif
         end subroutine fmomtdiag
!
         subroutine dmenergy(wt,wtot,msg,we,wf,wm,wke,wki,ntw,ndw,id0,it&
     &w,iuot,ntime)
! darwin electromagnetic energy diagnostic
         implicit none
         integer :: ntw, ndw, id0, itw, iuot, ntime
         real :: we, wf, wm, wke, wki
         real, dimension(:,:), pointer :: wt
         real, dimension(7) :: wtot
         double precision, dimension(:) :: msg
! local data
         integer :: it, jt
         real :: wef
  992    format (' field, kinetic, total energies = ',3e14.7)
  993    format (' electric(l,t), magnetic energies = ',3e14.7)
         if ((ntw > 0) .or. (ndw > 0)) then
            it = -1; if (ntw > 0) it = ntime - ntw*((ntime+1)/ntw) + 1
            jt = -1; if (ndw > 0) jt = ntime - ndw*((ntime+1)/ndw) + 1
            if ((it==0) .or. (jt==0)) then
               wef = we + wm
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
         end subroutine dmenergy
!
         subroutine emenergy(wt,wtot,msg,we,wf,wm,wke,wki,ntw,ndw,id0,it&
     &w,iuot,ntime)
! electromagnetic energy diagnostic
         implicit none
         integer :: ntw, ndw, id0, itw, iuot, ntime
         real :: we, wf, wm, wke, wki
         real, dimension(:,:), pointer :: wt
         real, dimension(7) :: wtot
         double precision, dimension(:) :: msg
! local data
         integer :: it, jt
         real :: wef
  992    format (' field, kinetic, total energies = ',3e14.7)
  993    format (' electric(l,t), magnetic energies = ',3e14.7)
         if ((ntw > 0) .or. (ndw > 0)) then
            it = -1; if (ntw > 0) it = ntime - ntw*((ntime+1)/ntw) + 1
            jt = -1; if (ndw > 0) jt = ntime - ndw*((ntime+1)/ndw) + 1
            if ((it==0) .or. (jt==0)) then
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
         end subroutine emenergy
!
      end module pemsimul2d