!-----------------------------------------------------------------------
!
      module psimul2d
! Higher level subroutines for electrostatics
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: november 13, 2009
!
      use globals, only: LINEAR, QUADRATIC
      use pdiag2d, only: vdist, sdist, displayfv, displays, grasp,      &
     &bfopen
      use pespush2d, only: dpost, push, rpush, pushzf, rpushzf,         &
     &premoment2, primoment2, fft, get_funit, plbcast, plsum, wrdata,   &
     &rddata, writebf, paguard, pcguard
      use pfield2d, only: sguard, aguard, cguard, spois, pois, gtmodes, &
     &imoment
      implicit none
      private
      public :: restart_open, restart_bwrite, restart_dwrite
      public :: restart_bread, restart_dread
      public :: dpostg, pushg, bpushg, pushzfg
      public :: initmodediag, initveldiag
      public :: dendiag, phasediag, veldiag, potdiag
      public :: emomtdiag, imomtdiag, esenergy
!
      contains
!
         subroutine restart_open(nustrt,ntr,idrun,iur1,iur2,iuer)
! open restart files
         implicit none
         integer, intent(in) :: nustrt, ntr, idrun
         integer :: iur1, iur2, iuer
! local data
         integer :: ierr
         character(len=10) :: cdrun
         character(len=32) :: fname
! create string from idrun
         write (cdrun,'(i10)') idrun
         cdrun = adjustl(cdrun)
         iur1 = 0; iur2 = 0
! open old restart files
         if ((nustrt==2).or.(nustrt==0)) then
            iur1 = get_funit(16)
            fname = 'rstrt1.'//cdrun
            open(unit=iur1,file=trim(fname),form='unformatted',status='o&
     &ld',iostat=ierr)
            if (ierr /= 0)  then
               iur1 = -1
               write (iuer,*) 'Cannot open restart file1=',fname
            endif
            iur2 = get_funit(17)
            fname = 'rstrt2.'//cdrun
            open(unit=iur2,file=trim(fname),form='unformatted',status='o&
     &ld',iostat=ierr)
            if (ierr /= 0)  then
               iur2 = -1
               write (iuer,*) 'Cannot open restart file2=',fname
            endif
! open new restart files
         else if (ntr > 0) then
            iur1 = get_funit(16)
            fname = 'rstrt1.'//cdrun
            open(unit=iur1,file=trim(fname),form='unformatted',status='u&
     &nknown')
            iur2 = get_funit(17)
            fname = 'rstrt2.'//cdrun
            open(unit=iur2,file=trim(fname),form='unformatted',status='u&
     &nknown')
         endif
         end subroutine restart_open
!
         subroutine restart_bwrite(kunit,id0,itime,itime0,nvp,npp,part,m&
     &ovion,nppi,parti,qi,q2m,ef,bf)
! write file for basic restart or continuation
         implicit none
         integer, intent(in) :: kunit, id0, itime, itime0, nvp, movion
         real, dimension(:,:,:), pointer :: part, parti
         integer, dimension(:), pointer :: npp, nppi
         real, dimension(:,:,:), pointer :: qi
         complex, dimension(:,:), pointer, optional :: q2m
         complex, dimension(:,:,:,:), pointer, optional :: ef, bf
! local data
         integer :: j, it
         character(len=1), dimension(8), save :: arch = ' '
         integer, external :: NDIAN, NDPREC, IDPREC
! determine architecture information
         if (arch(1)==' ') then
! determine if number format is big or little endian
            it = NDIAN()
            if (it==0) then
               arch(1) = 'L'
            else if (it==1) then
               arch(1) = 'B'
            endif
! determine if default reals are double precision
            it = NDPREC()
            if (it==0) then
               arch(2) = 'S'
            else if (it==1) then
               arch(2) = 'D'
            endif
! determine if default integers are double precision
            it = IDPREC()
            if (it==0) then
               arch(3) = 'S'
            else if (it==1) then
               arch(3) = 'D'
            endif
         endif
! write out architecture information
         if (id0==0) write (kunit) (arch(j),j=1,8)
! write out current and initial time
         if (id0==0) write (kunit) itime, itime0
! write out number of processors
         if (id0==0) write (kunit) nvp
! write out size of particle array
         if (id0==0) write (kunit) size(part,1)
! write out electrons, if non-zero
         call wrdata(part,npp,kunit)
! write out if ions are moving
         if (id0==0) write (kunit) movion
         if (movion > 0) then
! write out size of ion array
            if (id0==0) write (kunit) size(parti,1)
! write out number of ions, size of particle array
            call wrdata(parti,nppi,kunit)
         else if (movion==0) then
! write out ion density, if ions are not moving
            if (id0==0) write (kunit) size(qi)
            call wrdata(qi,nvp,kunit)
         endif
! write out shift constants
         it = 0
         if (present(q2m)) then
            if (id0==0) write (kunit) size(q2m)
            call wrdata(q2m,1,kunit)
         else
            if (id0==0) write (kunit) it
         endif
! write out electromagnetic fields, if present
         it = 0
         if (present(ef)) then
            if (id0==0) write (kunit) size(ef)
            call wrdata(ef,nvp,kunit)
         else
            if (id0==0) write (kunit) it
         endif
         it = 0
         if (present(bf)) then
            if (id0==0) write (kunit) size(bf)
            call wrdata(bf,nvp,kunit)
         else
            if (id0==0) write (kunit) it
         endif
! write out current time for later confirmation
         if (id0==0) write (kunit) itime
         end subroutine restart_bwrite
!
         subroutine restart_dwrite(kunit,id0,itime,itw,wt,ndrec,fdname,n&
     &prec,fpname,narec,faname,njrec,fjname,nerec,fename)
! write diagnostic file for restart or continuation
         implicit none
         integer, intent(in) :: kunit, id0, itime
         integer :: itw, ndrec, nprec
         real, dimension(:,:), pointer :: wt
         character(len=*), intent(in) :: fdname, fpname
         integer, optional :: narec, njrec, nerec
         character(len=*), intent(in), optional :: faname, fjname
         character(len=*), intent(in), optional :: fename
! local data
         integer :: it, na, nj, ne, irc
         character(len=32) :: fname
         na = 0; nj = 0; ne = 0
         if ((present(narec)).and.(present(faname))) na = 1
         if ((present(njrec)).and.(present(fjname))) nj = 1
         if ((present(nerec)).and.(present(fename))) ne = 1
! write out number of diagnostics in file
         it = 6
         if (id0==0) write (kunit) it
! write out current energy time step
         if (id0==0) write (kunit) itw
         if (itw > 0) then
! write out size of energy array
            if (id0==0) write (kunit) size(wt,2)
! write out energy values
            call wrdata(wt,1,kunit)
         endif
! write out ion density diagnostic write location
         if (id0==0) write (kunit) ndrec
! write out record length
         if (ndrec > 0) then
            if (id0==0) then
               inquire(file=fdname,recl=it,iostat=irc)
               if (irc /= 0) it = 0
               write (kunit) it
               if (it > 0) then
                  fname = fdname
                  write (kunit) fname
               endif
            endif
         endif
! write out potential diagnostic write location
         if (id0==0) write (kunit) nprec
! write out record length
         if (nprec > 0) then
            if (id0==0) then
               inquire(file=fpname,recl=it,iostat=irc)
               if (irc /= 0) it = 0
               write (kunit) it
               if (it > 0) then
                  fname = fpname
                  write (kunit) fname
               endif
            endif
         endif
! write out vector potential diagnostic write location
         if (na==1) then
            if (id0==0) write (kunit) narec
! write out record length
            if (narec > 0) then
               if (id0==0) then
                  inquire(file=faname,recl=it,iostat=irc)
                  if (irc /= 0) it = 0
                  write (kunit) it
                  if (it > 0) then
                     fname = faname
                     write (kunit) fname
                  endif
               endif
            endif
         else
            if (id0==0) write (kunit) na
         endif
! write out ion current diagnostic write location
         if (nj==1) then
            if (id0==0) write (kunit) njrec
! write out record length
            if (njrec > 0) then
               if (id0==0) then
                  inquire(file=fjname,recl=it,iostat=irc)
                  if (irc /= 0) it = 0
                  write (kunit) it
                  if (it > 0) then
                     fname = fjname
                     write (kunit) fname
                  endif
               endif
            endif
         else
            if (id0==0) write (kunit) nj
         endif
! write out electromagnetic diagnostic write location
         if (ne==1) then
            if (id0==0) write (kunit) nerec
! write out record length
            if (nerec > 0) then
               if (id0==0) then
                  inquire(file=fename,recl=it,iostat=irc)
                  if (irc /= 0) it = 0
                  write (kunit) it
                  if (it > 0) then
                     fname = fename
                  write (kunit) fname
                  endif
               endif
            endif
         else
            if (id0==0) write (kunit) ne
         endif
! write current time for later confirmation
         if (id0==0) then
            write (kunit) itime
            end file kunit
            rewind kunit
         endif
         end subroutine restart_dwrite
!
         subroutine restart_bread(iur1,iur2,id0,kunit,itime,itime0,nvp,n&
     &pp,part,movion,nppi,parti,qi,irc,iuer,q2m,ef,bf)
! read file for basic restart or continuation
         implicit none
         integer, intent(in) :: iur1, iur2, id0, nvp, movion
         integer :: iuer, kunit, itime, itime0, irc
         real, dimension(:,:,:), pointer :: part, parti
         integer, dimension(:), pointer :: npp, nppi
         real, dimension(:,:,:), pointer :: qi
         complex, dimension(:,:), pointer, optional :: q2m
         complex, dimension(:,:,:,:), pointer, optional :: ef, bf
! local data
         integer :: j, it
         integer, dimension(2) :: ktime
         character(len=1), dimension(8), save :: arch = ' '
         integer, external :: NDIAN, NDPREC, IDPREC
         if (id0==0) then
! determine most recent restart file
            if (kunit <= 0) then
               ktime(1) = -1
               read (iur1,iostat=irc) (arch(j),j=1,8)
               if (irc==0) then
                  read (iur1,iostat=irc) ktime(1)
                  if (irc /= 0) ktime(1) = -1
                  rewind iur1
               endif
               ktime(2) = -1
               read (iur2,iostat=irc) (arch(j),j=1,8)
               if (irc==0) then
                  read (iur2,iostat=irc) ktime(2)
                  if (irc /= 0) ktime(2) = -1
                  rewind iur2
               endif
               if (ktime(1) > ktime(2)) then
                  ktime(2) = iur1
               else if (ktime(2) >= 0) then
                  ktime(1) = ktime(2)
                  ktime(2) = iur2
               endif
! switch restart files
            else
               if (kunit==iur1) then
                  ktime(1) = 0
                  ktime(2) = iur2
               else if (kunit==iur2) then
                  ktime(1) = 0
                  ktime(2) = iur1
               else
                  ktime(1) = 0
                  ktime(2) = -1
               endif
            endif
         endif
         call plbcast(ktime)
! check if unit number is valid
         irc = 99; if (ktime(2) < 0) return
         kunit = ktime(2)
! read in architecture information
         if (id0==0)  then
            read (kunit,err=10) (arch(j),j=1,8)
! convert restart architecture information to integer
            ktime(1) = -1
! determine if number format is big or little endian
            if (arch(1)=='L') then
              ktime(1) = 0
            else if (arch(1)=='B') then
              ktime(1) = 1
            endif
! determine if default reals are double precision
            if (arch(2)=='S') then
               ktime(1) = 2*ktime(1)
            else if (arch(2)=='D') then
               ktime(1) = 2*ktime(1) + 1
            endif
! determine if default integers are double precision
            if (arch(3)=='S') then
               ktime(1) = 2*ktime(1)
            else if (arch(3)=='D') then
               ktime(1) = 2*ktime(1) + 1
            endif
! convert current architecture information to integer
            it = 4*NDIAN() + 2*NDPREC() + IDPREC()
            ktime(1) = abs(ktime(1)-it)
         endif
         call plbcast(ktime)
         irc = 1; if (ktime(1) > 0) go to 10
! read in current and initial time
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1), ktime(2)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 2; if (ktime(1) < 0) go to 10
         itime = ktime(1); itime0 = ktime(2)
! read in number of processors
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 3; if (ktime(1) /= nvp) go to 10
! read in size of particle array
         if (id0==0) then
           read (kunit,iostat=irc) ktime(1)
           if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 4; if (ktime(1) /= size(part,1)) go to 10
! read in electrons, if non-zero
         call rddata(part,npp,kunit,it)
         irc = 5; if (it /= 0) go to 10
! read in if ions are moving
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 6; if (ktime(1) /= movion) go to 10
         if (movion > 0) then
! read in size of ion array
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = 7; if (ktime(1) /= size(parti,1)) go to 10
! read in ions, if non-zero
            call rddata(parti,nppi,kunit,it)
            irc = 8; if (it /= 0) go to 10
! read in ion density, if ions are not moving
         else if (movion==0) then
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = 9; if (ktime(1) > size(qi)) go to 10
            call rddata(qi,nvp,kunit,it)
            irc = 10; if (it /= 0) go to 10
         endif
! read in shift constants, if present
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 11; if (ktime(1) < 0) go to 10
         if (ktime(1) > 0) then
            if (.not.present(q2m)) go to 10
            if (ktime(1) > size(q2m)) go to 10
            irc = 12; call rddata(q2m,1,kunit,it)
            if (it /= 0) go to 10
         endif
! read in first electromagnetic field, if present
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 13; if (ktime(1) < 0) go to 10
         if (ktime(1) > 0) then
            if (.not.present(ef)) go to 10
            if (ktime(1) > size(ef)) go to 10
            irc = 14; call rddata(ef,nvp,kunit,it)
            if (it /= 0) go to 10
         endif
! read in second electromagnetic field, if present
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 15; if (ktime(1) < 0) go to 10
         if (ktime(1) > 0) then
            if (.not.present(bf)) go to 10
            if (ktime(1) > size(bf)) go to 10
            irc = 16; call rddata(bf,nvp,kunit,it)
            if (it /= 0) go to 10
         endif
! read in current time for confirmation
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 17; if (ktime(1) /= itime) go to 10
         irc = 0
         return
! write out errors
   10    if (id0==0) then
            write (iuer,*) 'Basic Restart Error, irc = ', irc
            if (irc==1) then
               write (iuer,*) 'Architecture=', (arch(j),j=1,8)
            else if (irc==3) then
               write (iuer,*) 'nvp=', ktime(1)
            else if (irc==4) then
               write (iuer,*) 'size(part,1)=', ktime(1)
            else if (irc==6) then
               write (iuer,*) 'movion=', ktime(1)
            else if (irc==7) then
               write (iuer,*) 'size(parti,1)=', ktime(1)
            else if (irc==9) then
               write (iuer,*) 'size(qi)=', ktime(1)
            else if (irc==11) then
               write (iuer,*) 'size(q2m)=', ktime(1)
            else if (irc==13) then
               write (iuer,*) 'size(ef)=', ktime(1)
            else if (irc==15) then
               write (iuer,*) 'size(bf)=', ktime(1)
            else if (irc==17) then
               write (iuer,*) 'confirmation itime,it=',itime,ktime(1)
            endif
         endif
         end subroutine restart_bread
!
         subroutine restart_dread(kunit,id0,itime,itw,wt,iud,ndrec,fdnam&
     &e,iup,nprec,fpname,irc,iuer,iua,narec,faname,iuj,njrec,fjname,iue,&
     &nerec,fename)
! read diagnostic file for restart or continuation
         implicit none
         integer, intent(in) :: kunit, id0, itime
         integer :: itw, iud, ndrec, iup, nprec, irc, iuer
         real, dimension(:,:), pointer :: wt
         character(len=*) :: fdname, fpname
         integer, optional :: iua, narec, iuj, njrec, iue, nerec
         character(len=*), optional :: faname, fjname, fename
! local data
         integer :: it, nt, na, nj, ne, ierr
         integer, dimension(1) :: ktime
         character(len=32) :: fname
         na = 0; nj = 0; ne = 0
         if ((present(iua)).and.(present(narec)).and.(present(faname))) &
     &na = 1
         if ((present(iuj)).and.(present(njrec)).and.(present(fjname))) &
     &nj = 1
         if ((present(iue)).and.(present(nerec)).and.(present(fename))) &
     &ne = 1
! check if unit number is valid
         irc = -99; if (kunit < 0) return
! read in number of diagnostics in file
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -98; if (ktime(1) < 0) go to 20
         nt = ktime(1)
! read in current energy time step
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -1; if (ktime(1) < 0) go to 20
         itw = ktime(1)
         if (itw > 0) then
! read in size of energy array
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -2
            if ((itw > size(wt,1)).or.(ktime(1) > size(wt,2))) go to 20
! read in energy values
            irc = -3; call rddata(wt,1,kunit,it)
            if (it /= 0) go to 20
         endif
         if (nt==1) go to 10
! read in ion density diagnostic write location
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -4; if (ktime(1) < 0) go to 20
         ndrec = ktime(1)
! read in record length and open file
         if (ndrec > 0) then
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -5; if (ktime(1) < 1) go to 20
            it = ktime(1)
            if (id0==0) then
               ktime(1) = 0
               read (kunit,iostat=irc) fname
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -6; if (ktime(1) < 0) go to 20
            if (id0==0) then
               fdname = fname
               iud = get_funit(iud)
               open(unit=iud,file=fdname,form='unformatted',access='dire&
     &ct',recl=it,status='old',iostat=ktime(1))
            endif
            call plbcast(ktime)
            irc = -7; if (ktime(1) /= 0) go to 20
         endif
         if (nt==2) go to 10
! read in potential diagnostic write location
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -8; if (ktime(1) < 0) go to 20
         nprec = ktime(1)
! read in record length and open file
         if (nprec > 0) then
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -9; if (ktime(1) < 1) go to 20
            it = ktime(1)
            if (id0==0) then
               ktime(1) = 0
               read (kunit,iostat=irc) fname
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -10; if (ktime(1) < 0) go to 20
            if (id0==0) then
               fpname = fname
               iup = get_funit(iup)
               open(unit=iup,file=fpname,form='unformatted',access='dire&
     &ct',recl=it,status='old',iostat=ktime(1))
            endif
            call plbcast(ktime)
            irc = -11; if (ktime(1) /= 0) go to 20
         endif
         if (nt==3) go to 10
! read in vector potential diagnostic write location
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -12; if (ktime(1) < 0) go to 20
         it = ktime(1)
         if (na==1) then
            narec = it
! read in record length and open file
            if (narec > 0) then
               if (id0==0) then
                  read (kunit,iostat=irc) ktime(1)
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -13; if (ktime(1) < 1) go to 20
               it = ktime(1)
               if (id0==0) then
                  ktime(1) = 0
                  read (kunit,iostat=irc) fname
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -14; if (ktime(1) < 0) go to 20
               if (id0==0) then
                  faname = fname
                  iua = get_funit(iua)
                  open(unit=iua,file=faname,form='unformatted',access='d&
     &irect',recl=it,status='old',iostat=ktime(1))
               endif
               call plbcast(ktime)
               irc = -15; if (ktime(1) /= 0) go to 20
            endif
         else
            if (it > 0) go to 20
         endif
         if (nt==4) go to 10
! read in ion current diagnostic diagnostic write location
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -16; if (ktime(1) < 0) go to 20
         it = ktime(1)
         if (nj==1) then
            njrec = it
! read in record length and open file
            if (njrec > 0) then
               if (id0==0) then
                  read (kunit,iostat=irc) ktime(1)
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -17; if (ktime(1) < 1) go to 20
               it = ktime(1)
               if (id0==0) then
                  ktime(1) = 0
                  read (kunit,iostat=irc) fname
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -18; if (ktime(1) < 0) go to 20
               if (id0==0) then
                  fjname = fname
                  iuj = get_funit(iuj)
                  open(unit=iuj,file=fjname,form='unformatted',access='d&
     &irect',recl=it,status='old',iostat=ktime(1))
               endif
               call plbcast(ktime)
               irc = -19; if (ktime(1) /= 0) go to 20
            endif
         else
            if (it > 0) go to 20
         endif
         if (nt==5) go to 10
! read in electromagnetic diagnostic write location
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -20; if (ktime(1) < 0) go to 20
         it = ktime(1)
         if (ne==1) then
            nerec = it
! read in record length and open file
            if (nerec > 0) then
               if (id0==0) then
                  read (kunit,iostat=irc) ktime(1)
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -21; if (ktime(1) < 1) go to 20
               it = ktime(1)
               if (id0==0) then
                  ktime(1) = 0
                  read (kunit,iostat=irc) fname
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -22; if (ktime(1) < 0) go to 20
               if (id0==0) then
                  fename = fname
                  iue = get_funit(iue)
                  open(unit=iue,file=fename,form='unformatted',access='d&
     &irect',recl=it,status='old',iostat=ierr)
               endif
               call plbcast(ktime)
               irc = -23; if (ktime(1) /= 0) go to 20
            endif
         else
            if (it > 0) go to 20
         endif
! read in current time for confirmation
   10    if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -20; if (ktime(1) /= itime) go to 20
         irc = 0
         if (id0==0) rewind kunit
         return
! write out errors
   20    if (id0==0) then
            write (iuer,*) 'Diagnostic Restart Error, irc = ', irc
            if (irc==(-1)) then
               write (iuer,*) 'itw=', itw
            else if (irc==(-2)) then
               write (iuer,*) 'itw,size(wt,2)=', itw, it
            else if ((irc==(-5)).or.(irc==(-9)).or.(irc==(-13)).or.(irc=&
     &=(-17)).or.(irc==(-21))) then
               write (iuer,*) 'recl=', it
            else if ((irc==(-7)).or.(irc==(-11)).or.(irc==(-15)).or.(irc&
     &==(-19)).or.(irc==(-23))) then
               write (iuer,*) 'fname=', fname
            else if (irc==(-20)) then
               write (iuer,*) 'confirmation itime,it=', itime, it
            endif
         endif
         end subroutine restart_dread
!
         subroutine dpostg(part,q,npp,noff,nyp,kstrt,nvp,nx,kyp,ngds,qm,&
     &tdpost,inorder,dopt)
! deposit charge
         implicit none
         integer :: kstrt, nvp, nx, kyp, ngds
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: npp, noff, nyp
! local data
         real, save :: zero = 0.0
! initialize charge density to zero
         call sguard(q,nyp,zero,nx,inorder)
! deposit charge density
         call dpost(part,q,qm,npp,noff,tdpost,inorder,dopt)
! add guard cells for charge density
         call aguard(q,nyp,nx,inorder)
         call paguard(q,kstrt,nvp,nx,kyp,ngds)
         end subroutine dpostg
!
         subroutine pushg(part,fxy,npp,noff,qbm,dt,ci,ek,tpush,nx,ny,ipb&
     &c,relativity,inorder,popt)
! push particles with 2d electrostatic fields
         implicit none
         integer :: nx, ny, ipbc, relativity
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy
         integer, dimension(:), pointer :: npp, noff
! push particles
         if (relativity==1) then
            call rpush(part,fxy,npp,noff,qbm,dt,ci,ek,tpush,nx,ny,ipbc,i&
     &norder,popt)
         else
            call push(part,fxy,npp,noff,qbm,dt,ek,tpush,nx,ny,ipbc,inord&
     &er,popt)
         endif
         end subroutine pushg
!
         subroutine bpushg(part,fxy,bxy,npp,noff,qbm,dt,ci,ek,tpush,nx,n&
     &y,ipbc,relativity,inorder,popt)
! push particles with 2d electromagnetic fields
         implicit none
         integer :: nx, ny, ipbc, relativity
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, bxy
         integer, dimension(:), pointer :: npp, noff
! push particles
         if (relativity==1) then
            call rpush(part,fxy,bxy,npp,noff,qbm,dt,dt,ci,ek,tpush,nx,ny&
     &,ipbc,inorder,popt)
         else
            call push(part,fxy,bxy,npp,noff,qbm,dt,dt,ek,tpush,nx,ny,ipb&
     &c,inorder,popt)
         endif
         end subroutine bpushg
!
         subroutine pushzfg(part,npp,dt,ci,ek,tpush,nx,ny,ipbc,relativit&
     &y)
! push particles with no forces
         implicit none
         integer :: nx, ny, ipbc, relativity
         real :: dt, ci, ek, tpush
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! push particles
         if (relativity==1) then
            call rpushzf(part,npp,dt,ci,ek,tpush,nx,ny,ipbc)
         else
            call pushzf(part,npp,dt,ek,tpush,nx,ny,ipbc)
         endif
         end subroutine pushzfg
!
         subroutine initmodediag(dent,ntd,id0,nxh,nyh,kxp,modesxd,modesy&
     &d,jblok,iud,ndrec,fdname)
! initialize mode diagnostic
         implicit none
         integer :: ntd, id0, nxh, nyh, kxp
         integer :: modesxd, modesyd, jblok, iud, ndrec
         character(len=*) :: fdname
         complex, dimension(:,:,:), pointer :: dent
! local data
         integer :: modesy2d
         if (ntd <= 0) return
         if (modesxd > nxh) modesxd = nxh
         if (modesyd > nyh) modesyd = nyh
         modesy2d = 2*modesyd - 1
         allocate(dent(modesy2d,min(modesxd,kxp),jblok))
! open output file
         if (id0==0) then
            if (ndrec==0) then
               iud = get_funit(iud); ndrec = -1
               call bfopen(dent,modesy2d,iud,ndrec,trim(fdname))
            endif
         else
            if (ndrec==0) ndrec = 1
         endif
         end subroutine initmodediag
!
         subroutine initveldiag(fv,fvm,vtx,vty,vtz,ntv,ndv,id0,nmv,ndim,&
     &nblok,iuv,fvname)
! initialize velocity diagnostic
         implicit none
         integer :: ntv, ndv, id0, nmv, ndim, nblok, iuv
         real :: vtx, vty, vtz
         real, dimension(:,:,:), pointer :: fv, fvm
         character(len=*) :: fvname
         if ((ntv <= 0) .and. (ndv <= 0)) return
         allocate(fv(2*nmv+2,ndim,nblok),fvm(3,ndim,nblok))
! fix velocity range
         if (ndim==2) then
            fv(1,:,:) = 8.*max(vtx,vty)
         elseif (ndim==3) then
            fv(1,:,:) = 8.*max(vtx,vty,vtz)
         endif
         if (ntv > 0) then
            if (id0==0) then
               iuv = get_funit(iuv)
               open(unit=iuv,file=trim(fvname),form='formatted',status='&
     &unknown')
! write captions
               if (ndim==2) then
                  write (iuv,*) 'it vdx vdy vtx vty sk'
               else if (ndim==3) then
                  write (iuv,*) 'it vdx vdy vdz vtx vty vtz sk'
               endif
            endif
         endif
         end subroutine initveldiag
!
         subroutine dendiag(qt,qi,sfield,dent,sfieldt,ffc,nyp,mixup,sct,&
     &tfft,ntd,ndd,nx,ny,modesxd,modesyd,iud,ndrec,indx,indy,ntime,nvp,k&
     &strt,kxp,kyp,ndstyle,irc,inorder)
! ion density diagnostic
         implicit none
         integer :: ntd, ndd, nx, ny, modesxd, modesyd, iud, ndrec
         integer :: indx, indy, ntime, nvp, kstrt, kxp, kyp
         integer :: ndstyle, irc
         real, dimension(2) :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: qi, sfield
         complex, dimension(:,:,:), pointer :: qt, sfieldt, ffc, dent
         integer, dimension(:), pointer :: nyp
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2d, isign
         irc = 0
         if ((ntd > 0) .or. (ndd > 0)) then
            it = -1; if (ntd > 0) it = ntime - ntd*(ntime/ntd)
            jt = -1; if (ndd > 0) jt = ntime - ndd*(ntime/ndd)
            if ((it==0) .or. (jt==0)) then
               sfield = qi
! transform ion density to fourier space
               isign = -1
               call fft(sfield,qt,isign,mixup,sct,tfft,indx,indy,kstrt,k&
     &yp,inorder)
! calculate smoothing in fourier space
               call spois(qt,sfieldt,ffc,nx,ny,kstrt)
! store selected fourier modes
               if (it==0) then
                  modesy2d = 2*modesyd - 1
                  call gtmodes(sfieldt,dent,nx,ny,modesxd,modesyd,kstrt)
! write diagnostic output
                  call writebf(dent,modesxd,modesy2d,kxp,iud,ndrec)
               endif
! transform ion density to real space
               if (jt==0) then
                  isign = 1
                  call fft(sfield,sfieldt,isign,mixup,sct,tfft,indx,indy&
     &,kstrt,kyp,inorder)
! copy to guard cells
                  call pcguard(sfield,kstrt,nvp,kyp,inorder)
                  call cguard(sfield,nyp,nx,inorder)
! display ion density
                  call displays(sfield,nvp,' ION DENSITY',ntime,999,2,nd&
     &style,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine dendiag
!
         subroutine veldiag(part,fv,fvm,npp,msg,ntv,ndv,id0,nmv,iuv,ntim&
     &e,label,irc)
! velocity diagnostic
         implicit none
         integer :: ntv, ndv, id0, nmv, iuv, ntime, irc
         character(len=*) :: label
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fv, fvm
         integer, dimension(:), pointer :: npp
         double precision, dimension(:) :: msg
! local data
         integer :: j, k, it, jt, nt, idimv
         real, dimension(size(fvm,2)) :: scale
         irc = 0
         idimv = size(fvm,2)
         if ((ntv > 0) .or. (ndv > 0)) then
            it = -1; if (ntv > 0) it = ntime - ntv*(ntime/ntv)
            jt = -1; if (ndv > 0) jt = ntime - ndv*(ntime/ndv)
            if ((it==0) .or. (jt==0)) then
! calculate particle distribution function and moments
               call vdist(part,fv,fvm,npp,nmv,idimv)
! sum distribution over processors, but save velocity scale
               scale = fv(1,:,1)
               call plsum(fv(:,:,1))
               fv(1,:,1) = scale
! calculate entropy
               call sdist(fv,fvm,nmv,idimv)
! send moments to diagnostic node
               do k = 1, idimv
               do j = 1, 3
               msg(j+3*(k-1)) = fvm(j,k,1)
               enddo
               enddo
               call HARTBEAT(msg,3*idimv)
! print out velocity moments
               if (it==0) then
                  if (id0==0) then
                     nt = ntime/ntv
                     write (iuv,*) nt, fvm(1,:,1), fvm(2,:,1), sum(fvm(3&
     &,:,1))
                  endif
               endif
! display velocity distributions
               if (jt==0) then
                  call displayfv(fv,fvm,label,ntime,nmv,2,idimv,irc)
               endif
            endif
         endif
         end subroutine veldiag
!
         subroutine phasediag(part,npp,nts,nds,nx,ny,ntime,label,irc)
! phase space diagnostic
         implicit none
         integer :: nts, nds, nx, ny, ntime, irc
         character(len=*) :: label
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
! local data
         integer :: it, jt, isc
         irc = 0
         if ((nts > 0) .or. (nds > 0)) then
            it = -1; if (nts > 0) it = ntime - nts*(ntime/nts)
            jt = -1; if (nds > 0) jt = ntime - nds*(ntime/nds)
            if ((it==0) .or. (jt==0)) then
               isc = 999
! plot electrons vx versus x
               call grasp(part,npp,label,ntime,isc,nx,ny,3,1,irc)
               if (irc==1) return
! plot electrons vy versus y
               call grasp(part,npp,label,ntime,isc,nx,ny,4,2,irc)
            endif
         endif
         end subroutine phasediag
!
         subroutine potdiag(qt,sfield,pott,sfieldt,ffc,nyp,mixup,sct,tff&
     &t,ntp,ndp,nx,ny,modesxp,modesyp,iup,nprec,indx,indy,ntime,nvp,kstr&
     &t,kxp,kyp,ndstyle,irc,inorder)
! potential diagnostic
         implicit none
         integer :: ntp, ndp, nx, ny, modesxp, modesyp, iup, nprec
         integer :: indx, indy, ntime, nvp, kstrt, kxp, kyp
         integer :: ndstyle, irc
         real, dimension(2) :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: sfield
         complex, dimension(:,:,:), pointer :: qt, pott, sfieldt
         complex, dimension(:,:,:), pointer :: ffc
         integer, dimension(:), pointer :: nyp
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2p, isign
         real :: we
         irc = 0
         if ((ntp > 0) .or. (ndp > 0)) then
            it = -1; if (ntp > 0) it = ntime - ntp*(ntime/ntp)
            jt = -1; if (ndp > 0) jt = ntime - ndp*(ntime/ndp)
            if ((it==0) .or. (jt==0)) then
! calculate potential in fourier space
               call pois(qt,sfieldt,ffc,we,nx,ny,kstrt)
! store selected fourier modes
               if (it==0) then
                  modesy2p = 2*modesyp - 1
                  call gtmodes(sfieldt,pott,nx,ny,modesxp,modesyp,kstrt)
! write diagnostic output
                  call writebf(pott,modesxp,modesy2p,kxp,iup,nprec)
               endif
! transform potential to real space
               if (jt==0) then
                  isign = 1
                  call fft(sfield,sfieldt,isign,mixup,sct,tfft,indx,indy&
     &,kstrt,kyp,inorder)
! copy to guard cells
                  call pcguard(sfield,kstrt,nvp,kyp,inorder)
                  call cguard(sfield,nyp,nx,inorder)
! display potential
                  call displays(sfield,nvp,' POTENTIAL',ntime,999,0,ndst&
     &yle,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine potdiag
!
         subroutine emomtdiag(part,npp,msg,pxe,pye,pze,sx,sy,sz,wx,wy,wz&
     &,ntm,id0,ium,ntime,ndim)
! calculate electron momentum
         implicit none
         integer :: ntm, id0, ium, ntime, ndim
         real :: pxe, pye, pze, sx, sy, sz, wx, wy, wz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp
         double precision, dimension(:) :: msg
! local data
         integer :: it
         if (ntm > 0) then
            it = ntime/ntm
! calculate the momentum in the electrons
            it = ntime - ntm*it + 1
            if (it > 1) it = it - ntm
            if (it >= 0) then
               call premoment2(part,ntime,npp,msg,id0,ium,pxe,pye,pze,sx&
     &,sy,sz,wx,wy,wz,ndim,nprint=it)
            endif
! send momentum values to diagnostic node
            if (it==1) call HARTBEAT(msg,6)
         endif
         end subroutine emomtdiag
!
         subroutine imomtdiag(parti,qi,fxy,nppi,nyp,msg,dt,rmass,pxi,pyi&
     &,pzi,wx,wy,wz,ntm,movion,id0,ium,nx,ntime,ndim,inorder)
! calculate ion and field momentum
         implicit none
         integer :: ntm, movion, id0, ium, nx, ntime, ndim
         real :: dt, rmass, pxi, pyi, pzi, wx, wy, wz
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: parti
         integer, dimension(:), pointer :: nppi, nyp
         real, dimension(:,:,:), pointer :: qi
         real, dimension(:,:,:,:), pointer :: fxy
         double precision, dimension(:) :: msg
! local data
         integer :: it
  997    format (' total momentum = ',3e14.7)
         if (ntm > 0) then
            it = ntime/ntm
            it = ntime - ntm*it + 1
            if (it > 1) it = it - ntm
            if (it >= 0) then
! calculate ion momentum
               if (movion==0) then
                  if (it==1) then
                     call imoment(qi,fxy,nyp,msg,id0,ium,pxi,pyi,pzi,dt,&
     &wx,wy,wz,nx,inorder)
                  endif
               else if (movion==1) then
                  call primoment2(parti,nppi,msg,id0,ium,rmass,pxi,pyi,p&
     &zi,wx,wy,wz,ndim,nprint=it)
               endif
! print total momentum
               if (it==1) then
! send momentum values to diagnostic node
                  call HARTBEAT(msg,6)
                  if (id0==0) then
                     write (ium,997) wx, wy, wz
                  endif
               endif
            endif
         endif
         end subroutine imomtdiag
!
         subroutine esenergy(wt,wtot,msg,we,wke,wki,ntw,ndw,id0,itw,iuot&
     &,ntime)
! electrostatic energy diagnostic
         implicit none
         integer :: ntw, ndw, id0, itw, iuot, ntime
         real :: we, wke, wki
         real, dimension(:,:), pointer :: wt
         real, dimension(4) :: wtot
         double precision, dimension(:) :: msg
! local data
         integer :: it, jt
  992    format (' field, kinetic, total energies = ',3e14.7)
         if ((ntw > 0) .or. (ndw > 0)) then
            it = -1; if (ntw > 0) it = ntime - ntw*((ntime+1)/ntw) + 1
            jt = -1; if (ndw > 0) jt = ntime - ndw*((ntime+1)/ndw) + 1
            if ((it==0) .or. (jt==0)) then
               wtot(1) = we
               wtot(2) = wke
               wtot(3) = wki
               wtot(4) = we + wke + wki
               call plsum(wtot)
! send energy values to diagnostic node
               msg(1:4) = wtot
               call HARTBEAT(msg,4)
               if (it==0) then
                  if (id0==0) write (iuot,992) wtot(1), wtot(2), wtot(4)
               endif
               if (jt==0) then
                  itw = itw + 1
                  wt(itw,:) = wtot
               endif
            endif
         endif
         end subroutine esenergy
!
      end module psimul2d