#Makefile for 2D parallel PIC codes in new_pbeps2.source

GOBJS = nullpgks2.o nullpgks1.o
CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

# Makefile Absoft compiler with MacOS X

# MacMPI
MPIFC = f90
MPIFC77 = f77
MPIOBJS = MacMPIf77.o LnxMPI_S.o

FC90 = f90
FC77 = f77
CC = gcc

OPTS90 = -O3 -N113 -YEXT_NAMES=LCS -YEXT_SFX=_
OPTS77 = -O -N113 -f -N15
CCOPTS = -O
MOPTS = -s
MBOPTS = -s -N11
LOPTS = -plainappl
LEGACY =

MPOBJS = MacMPf77.o LnxMP.o

# Mac graphics
#GOBJS = plibgks2.o plibgks1.o libgks2.o libgks1.o libmcX.o
#LIBS = $(CARBON)
# No graphics
LIBS = $(CARBON)

# Makefile Nag compiler with MacOS X

# MacMPI
#MPIFC = f95
#MPIFC77 = f95
#MPIOBJS = MacMPIf77.o LnxMPI_S.o

#FC90 = f95
#FC77 = f95
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS = -framework carbon 
#LEGACY = -dusty

#MPOBJS = MacMPf77.o LnxMP.o

# No graphics
#LIBS =

# Makefile IBM compiler with MacOS X

# MacMPI
#MPIFC = xlf90
#MPIFC77 = xlf
#MPIOBJS = MacMPIf77.o LnxMPI_S.o

#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qarch=g5 -qextname
#OPTS77 = -O3 -qautodbl=dbl4 -qarch=g5 -qextname -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o

# No graphics
#LIBS = $(CARBON)

# Makefile IBM compiler with LAM MPI and MacOS X

# LAM MPI
#MPIFC = mpif77
#MPIFC77 = mpif77
#MPIOBJS = nullLOG.o

#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qarch=g5
#OPTS77 = -O3 -qautodbl=dbl4 -qarch=g5 -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

#MPOBJS = MacMPxlf.o LnxMP.o

# No graphics
#LIBS = $(CARBON)

# Makefile Intel compiler with OpenMPI and MacOS X

# MacMPI
#MPIFC = ifort
#MPIFC77 = ifort
#MPIOBJS = MacMPIf77.o LnxMPI_S.o
# OpenMPI
#MPIFC = mpif90
#MPIFC77 = mpif77
#MPIOBJS = nullLOG.o

#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = plibgks2.o plibgks1.o libgks2.o libgks1.o libygl.o
#LIBS = $(CARBON) -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# No graphics
#LIBS = $(CARBON)

# Makefile g95 compiler with MacOS X

# MacMPI
#MPIFC = g95
#MPIFC77 = g95
#MPIOBJS = MacMPIf77.o LnxMPI_S.o
# OpenMPI
#MPIFC = openmpif90
#MPIFC77 = openmpif77
#MPIOBJS = nullLOG.o

#FC90 = g95
#FC77 = g95
#CC = gcc

#OPTS90 = -O3 -r8 -fno-second-underscore
#OPTS77 = -O3 -r8 -fno-second-underscore
#CCOPTS = -O
#MOPTS = -fstatic
#MBOPTS = -fstatic
#LOPTS =
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = plibgks2.o plibgks1.o libgks2.o libgks1.o libygl.o
#LIBS = $(CARBON) -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# No graphics
#LIBS = $(CARBON) -lSystemStubs

# Makefile gfortran compiler with MacOS X

# MacMPI
#MPIFC = gfortran
#MPIFC77 = gfortran
#MPIOBJS = MacMPIf77.o LnxMPI_S.o
# OpenMPI
#MPIFC = mpif90
#MPIFC77 = mpif77
#MPIOBJS = nullLOG.o

#FC90 = gfortran
#FC77 = gfortran
#CC = gcc

#OPTS90 = -O3 -fdefault-real-8
#OPTS77 = -O3 -fdefault-real-8
#CCOPTS = -O
#MOPTS = -fno-automatic
#MBOPTS = -fno-automatic
#LOPTS =
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = plibgks2.o plibgks1.o libgks2.o libgks1.o libygl.o
#LIBS = $(CARBON) -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# No graphics
#LIBS = $(CARBON) -lSystemStubs

# Makefile Intel compiler with Linux

# MacMPI
#MPIFC = ifc
#MPIFC77 = ifc
#MPIOBJS = MacMPIf77.o LnxMPI_S.o

#FC90 = ifc
#FC77 = ifc
#CC = gcc

#OPTS90 = -O3 -tpp7 -xW -r8
#OPTS77 = -O3 -tpp7 -xW -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

# No graphics
#LIBS =

# Makefile Intel compiler with MPI and Linux

#MPIFC = mpiifort
#MPIFC77 = mpiifort
#MPIOBJS = nullLOG.o

#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -tpp7 -xW -r8
#OPTS77 = -O3 -tpp7 -xW -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

# No graphics
#LIBS =

# Makefile gfortran compiler with Linux

# MacMPI
#MPIFC = gfortran
#MPIFC77 = gfortran
#MPIOBJS = MacMPIf77.o LnxMPI_S.o
# OpenMPI
#MPIFC = mpif90
#MPIFC77 = mpif77
#MPIOBJS = nullLOG.o

#FC90 = gfortran
#FC77 = gfortran
#CC = gcc

#OPTS90 = -O3 -fdefault-real-8
#OPTS77 = -O3 -fdefault-real-8
#CCOPTS = -O
#MOPTS = -fno-automatic
#MBOPTS = -fno-automatic
#LOPTS = -lpthread
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = plibgks2.o plibgks1.o libgks2.o libgks1.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
#LIBS = -L/usr/lib64 -lYgl -L/usr/X11R6/lib -lX11
# No graphics
#LIBS =

# Makefile Pathscale compiler with AMD Opteron and Linux

#MPIFC = mpipathf90
#MPIFC77 = mpipathf90
#MPIOBJS = nullLOG.o

#FC90 = pathf90
#FC77 = pathf90
#CC = gcc

#OPTS90 = -O3 -r8 
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -static-data
#MBOPTS = -static-data
#LOPTS = -lpthread
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

# No graphics
#LIBS =

# Makefile PGI compiler with AMD Opteron and Linux

#MPIFC = mpipgf90
#MPIFC77 = mpipgf90
#MPIOBJS = nullLOG.o

#FC90 = pgf90
#FC77 = pgf90
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -Msave
#MBOPTS = -Msave
#LOPTS = -lpthread
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

#LIBS =

#

GESOBJS = pespush2mod.o psimul2mod.o

MGESOBJS = mpespush2mod.o mpsimul2mod.o

ESOBJS = globals.o pinit2mod.o pdiag2mod.o prbpush2mod.o pbpush2mod.o \
ppush2mod.o pfft2mod.o pfield2mod.o p2mod.o p0mod.o mp0mod.o pinit2lib.o \
prbpush2lib.o pbpush2lib.o ppush2lib.o pfft2lib.o pfield2lib.o pdiag2lib.o \
p2lib.o p0lib.o

MESOBJS = mprbpush2mod.o mpbpush2mod.o mppush2mod.o mpfft2mod.o mp2mod.o \
mprbpush2lib.o mpbpush2lib.o mppush2lib.o mpfft2lib.o mp2lib.o 

GEMOBJS = pempush2mod.o pemsimul2mod.o psimul2mod.o

MGEMOBJS = mpempush2mod.o mpemsimul2mod.o mpsimul2mod.o

EMOBJS = prdpush2mod.o pdpush2mod.o prdpush2lib.o pdpush2lib.o

MEMOBJS = mprdpush2mod.o mpdpush2mod.o mprdpush2lib.o mpdpush2lib.o

DESOBJS = pdfield2mod.o pbfield2mod.o pcfield2mod.o pnfield2mod.o \
pnpfield2mod.o pdfield2lib.o pbfield2lib.o pcfield2lib.o pnfield2lib.o 

NPOBJS = nullMP.o

# Linkage rules

all : new_pbeps2.out new_d0_pbeps2.out new_pbbeps2.out new_d0_pbbeps2.out \
      new_pdbeps2.out

threaded: new_mpbeps2.out new_d0_mpbeps2.out new_mpbbeps2.out new_d0_mpbbeps2.out \
          new_mpdbeps2.out 

new_pbeps2.out : new_pbeps2.o $(GESOBJS) $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_pbeps2.out \
	new_pbeps2.o $(GESOBJS) $(ESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_d0_pbeps2.out : new_d0_pbeps2.o $(GESOBJS) $(ESOBJS) $(DESOBJS) $(MPIOBJS) \
                    $(NPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_pbeps2.out \
	new_d0_pbeps2.o $(GESOBJS) $(ESOBJS) $(DESOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) \
        $(LIBS)

new_mpbeps2.out : new_mpbeps2.o $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPIOBJS) $(MPOBJS) \
                  $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_mpbeps2.out \
        new_mpbeps2.o  $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPIOBJS) $(MPOBJS) $(GOBJS) \
        $(LIBS)

new_d0_mpbeps2.out : new_d0_mpbeps2.o $(MGESOBJS) $(ESOBJS) $(DESOBJS) $(MESOBJS) \
                     $(MPIOBJS) $(MPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_mpbeps2.out \
	new_d0_mpbeps2.o $(MGESOBJS) $(ESOBJS) $(DESOBJS) $(MESOBJS) $(MPIOBJS) \
        $(MPOBJS) $(GOBJS) $(LIBS)

new_pbbeps2.out : new_pbbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(MPIOBJS) $(NPOBJS) \
                  $(GOBJS)
	$(MPIFC) $(OPTS90) $(OPTS90) $(LOPTS) -o new_pbbeps2.out \
	new_pbbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) \
        $(LIBS)

new_d0_pbbeps2.out : new_d0_pbbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(DESOBJS) \
                     $(MPIOBJS) $(NPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_pbbeps2.out \
	new_d0_pbbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(DESOBJS) $(MPIOBJS) \
        $(NPOBJS) $(GOBJS) $(LIBS)

new_mpbbeps2.out : new_mpbbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                   $(MEMOBJS) $(MPIOBJS) $(MPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_mpbbeps2.out \
        new_mpbbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
        $(MPIOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

new_d0_mpbbeps2.out : new_d0_mpbbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                     $(MEMOBJS) $(DESOBJS) $(MPIOBJS) $(MPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_d0_mpbbeps2.out \
        new_d0_mpbbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
        $(DESOBJS) $(MPIOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

new_pdbeps2.out : new_pdbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(MPIOBJS) $(NPOBJS) \
                  $(GOBJS)
	$(MPIFC) $(OPTS90) $(OPTS90) $(LOPTS) -o new_pdbeps2.out \
	new_pdbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(MPIOBJS) $(NPOBJS) $(GOBJS) \
        $(LIBS)

new_mpdbeps2.out : new_mpdbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                   $(MEMOBJS) $(MPIOBJS) $(MPOBJS) $(GOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o new_mpdbeps2.out \
        new_mpdbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
        $(MPIOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

# Compilation rules

#MacMPIcf.o : MacMPIcf.c
#	$(CC) $(CCOPTS) -c MacMPIcf.c

MacMPIf77.o : MacMPIf77.c
	$(CC) $(CCOPTS) -c MacMPIf77.c

#MacMPIxlf.o : MacMPIxlf.c
#	$(CC) $(CCOPTS) -c MacMPIxlf.c

#MacMPI_S.o : MacMPI_S.c
#	$(CC) $(CCOPTS) -c -I /Developer/Headers/FlatCarbon MacMPI_S.c

LnxMPI_S.o : LnxMPI_S.c
	$(CC) $(CCOPTS) -c LnxMPI_S.c

#MacMPcf.o : MacMPcf.c
#	$(CC) $(CCOPTS) -c MacMPcf.c

MacMPf77.o : MacMPf77.c
	$(CC) $(CCOPTS) -c MacMPf77.c

#MacMPxlf.o : MacMPxlf.c
#	$(CC) $(CCOPTS) -c MacMPxlf.c

#MacMP.o : MacMP.c
#	$(CC) $(CCOPTS) -c -I /Developer/Headers/FlatCarbon MacMP.c

LnxMP.o : LnxMP.c
	$(CC) $(CCOPTS) -c LnxMP.c

LPProcessors.o : LPProcessors.c
	$(CC) $(CCOPTS) -c LPProcessors.c

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

libmcX.o : libmcX.f
	$(FC77) $(OPTS77) -c libmcX.f

libygl.o : libygl.f
	$(FC77) $(OPTS77) -c libygl.f

libgks1.o : libgks1.f
	$(FC77) $(OPTS77) -c libgks1.f

libgks2.o : libgks2.f
	$(FC77) $(OPTS77) -c libgks2.f

plibgks1.o : plibgks1.f
	$(FC77) $(OPTS77) -c plibgks1.f

plibgks2.o : plibgks2.f
	$(FC77) $(OPTS77) -c plibgks2.f

nullpgks1.o : nullpgks1.f
	$(FC77) $(OPTS77) -c nullpgks1.f

nullpgks2.o : nullpgks2.f
	$(FC77) $(OPTS77) -c nullpgks2.f

nullLOG.o : nullLOG.f
	$(FC77) $(OPTS77) -c nullLOG.f

nullMP.o : nullMP.f
	$(FC77) $(OPTS77) -c nullMP.f

p0lib.o : p0lib.f
	$(MPIFC77) $(OPTS77) $(LEGACY) -c p0lib.f

mp2lib.o : mp2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mp2lib.f

p2lib.o : p2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c p2lib.f

pinit2lib.o : pinit2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c pinit2lib.f

mppush2lib.o : mppush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mppush2lib.f

ppush2lib.o : ppush2lib.f
	$(FC77) $(OPTS77) -c ppush2lib.f

mpbpush2lib.o : mpbpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mpbpush2lib.f

pbpush2lib.o : pbpush2lib.f
	$(FC77) $(OPTS77) -c pbpush2lib.f

mpdpush2lib.o : mpdpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mpdpush2lib.f

pdpush2lib.o : pdpush2lib.f
	$(FC77) $(OPTS77) -c pdpush2lib.f

mprbpush2lib.o : mprbpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mprbpush2lib.f

prbpush2lib.o : prbpush2lib.f
	$(FC77) $(OPTS77) -c prbpush2lib.f

mprdpush2lib.o : mprdpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mprdpush2lib.f

prdpush2lib.o : prdpush2lib.f
	$(FC77) $(OPTS77) -c prdpush2lib.f

mpfft2lib.o : mpfft2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c mpfft2lib.f

pfft2lib.o : pfft2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c pfft2lib.f

pfield2lib.o : pfield2lib.f
	$(FC77) $(OPTS77) -c pfield2lib.f

pdfield2lib.o : pdfield2lib.f
	$(FC77) $(OPTS77) -c pdfield2lib.f

pbfield2lib.o : pbfield2lib.f
	$(FC77) $(OPTS77) -c pbfield2lib.f

pcfield2lib.o : pcfield2lib.f
	$(FC77) $(OPTS77) -c pcfield2lib.f

pnfield2lib.o : pnfield2lib.f
	$(FC77) $(OPTS77) -c pnfield2lib.f

pdiag2lib.o : pdiag2lib.f
	$(FC77) $(OPTS77) -c pdiag2lib.f

globals.o : globals.f
	$(FC90) $(OPTS90) -c globals.f

pinit2mod.o : pinit2mod.f globals.o
	$(FC90) $(OPTS90) -c pinit2mod.f

mppush2mod.o : mppush2mod.f ppush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mppush2mod.f

ppush2mod.o : ppush2mod.f p0mod.o
	$(FC90) $(OPTS90) -c ppush2mod.f

mpbpush2mod.o : mpbpush2mod.f pbpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mpbpush2mod.f

pbpush2mod.o : pbpush2mod.f p0mod.o
	$(FC90) $(OPTS90) -c pbpush2mod.f

mpdpush2mod.o : mpdpush2mod.f pdpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mpdpush2mod.f

pdpush2mod.o : pdpush2mod.f p0mod.o
	$(FC90) $(OPTS90) -c pdpush2mod.f

mprbpush2mod.o : mprbpush2mod.f prbpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mprbpush2mod.f

prbpush2mod.o : prbpush2mod.f p0mod.o
	$(FC90) $(OPTS90) -c prbpush2mod.f

mprdpush2mod.o : mprdpush2mod.f prdpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mprdpush2mod.f

prdpush2mod.o : prdpush2mod.f p0mod.o
	$(FC90) $(OPTS90) -c prdpush2mod.f

mpfft2mod.o : mpfft2mod.f pfft2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mpfft2mod.f

pfft2mod.o : pfft2mod.f p0mod.o
	$(FC90) $(OPTS90) -c pfft2mod.f

pfield2mod.o : pfield2mod.f p0mod.o
	$(FC90) $(OPTS90) -c pfield2mod.f

pdfield2mod.o : pdfield2mod.f globals.o
	$(FC90) $(OPTS90) -c pdfield2mod.f

pbfield2mod.o : pbfield2mod.f globals.o
	$(FC90) $(OPTS90) -c pbfield2mod.f

pcfield2mod.o : pcfield2mod.f globals.o
	$(FC90) $(OPTS90) -c pcfield2mod.f

pnfield2mod.o : pnfield2mod.f globals.o
	$(FC90) $(OPTS90) -c pnfield2mod.f

pnpfield2mod.o : pnpfield2mod.f pfield2mod.o pdfield2mod.o pbfield2mod.o \
                 pcfield2mod.o pnfield2mod.o
	$(FC90) $(OPTS90) -c pnpfield2mod.f

pdiag2mod.o : pdiag2mod.f pinit2mod.o
	$(FC90) $(OPTS90) -c pdiag2mod.f

p0mod.o : p0mod.f globals.o
	$(FC90) $(OPTS90) $(LEGACY) -c p0mod.f

mp2mod.o : mp2mod.f p2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mp2mod.f

p2mod.o : p2mod.f p0mod.o
	$(FC90) $(OPTS90) $(LEGACY) -c p2mod.f

mpespush2mod.o : mpespush2mod.f mprbpush2mod.o mpbpush2mod.o mppush2mod.o \
                 mpfft2mod.o mp2mod.o
	$(FC90) $(OPTS90) -c mpespush2mod.f

pespush2mod.o : pespush2mod.f prbpush2mod.o pbpush2mod.o ppush2mod.o pfft2mod.o \
                p2mod.o
	$(FC90) $(OPTS90) -c pespush2mod.f

mpempush2mod.o : mpempush2mod.f mprbpush2mod.o mprdpush2mod.o mpbpush2mod.o \
                 mpdpush2mod.o mppush2mod.o mpfft2mod.o mp2mod.o
	$(FC90) $(OPTS90) -c mpempush2mod.f

pempush2mod.o : pempush2mod.f prbpush2mod.o prdpush2mod.o pbpush2mod.o \
                pdpush2mod.o ppush2mod.o pfft2mod.o p2mod.o
	$(FC90) $(OPTS90) -c pempush2mod.f

psimul2mod.o : psimul2mod.f pdiag2mod.o pespush2mod.o pfield2mod.o
	$(FC90) $(OPTS90) -c psimul2mod.f

mpsimul2mod.o : psimul2mod.f pdiag2mod.o mpespush2mod.o pfield2mod.o
	$(FC90) $(OPTS90) -o mpsimul2mod.o -c psimul2mod.f

pemsimul2mod.o : pemsimul2mod.f psimul2mod.o pempush2mod.o
	$(FC90) $(OPTS90) -c pemsimul2mod.f

mpemsimul2mod.o : pemsimul2mod.f mpsimul2mod.o mpempush2mod.o
	$(FC90) $(OPTS90) -o mpemsimul2mod.o -c pemsimul2mod.f

mp0mod.o : mp0mod.f
	$(FC90) $(OPTS90) -c mp0mod.f

new_pbeps2.o : new_pbeps2.f psimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_pbeps2.f

new_mpbeps2.o : new_pbeps2.f mpsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mpbeps2.o -c new_pbeps2.f

new_pdbeps2.o : new_pdbeps2.f pemsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_pdbeps2.f

new_mpdbeps2.o : new_pdbeps2.f mpemsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mpdbeps2.o -c new_pdbeps2.f

new_pbbeps2.o : new_pbbeps2.f pemsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_pbbeps2.f

new_mpbbeps2.o : new_pbbeps2.f mpemsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mpbbeps2.o -c new_pbbeps2.f

new_d0_pbeps2.o : new_d0_pbeps2.f pespush2mod.o pnpfield2mod.o pdiag2mod.o \
                  mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_d0_pbeps2.f

new_d0_mpbeps2.o : new_d0_pbeps2.f mpespush2mod.o pnpfield2mod.o pdiag2mod.o \
                   mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_d0_mpbeps2.o -c new_d0_pbeps2.f

new_d0_pbbeps2.o : new_d0_pbbeps2.f pempush2mod.o pnpfield2mod.o pdiag2mod.o \
                   mp0mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -c new_d0_pbbeps2.f

new_d0_mpbbeps2.o : new_d0_pbbeps2.f mpempush2mod.o pnpfield2mod.o pdiag2mod.o \
                    mp0mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -o new_d0_mpbbeps2.o -c new_d0_pbbeps2.f

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f *.out
