# Makefile for new_d0_mpbbeps2 with MacOS 9OPTS90 = -O -N113OPTS77 = -O -N113 -N2OBJS = new_d0_pbbeps2.f.o �globals.f.o pinit2mod.f.o mpempush2mod.f.o mprbpush2mod.f.o mpbpush2mod.f.o mppush2mod.f.o �prbpush2mod.f.o pbpush2mod.f.o ppush2mod.f.o mpfft2mod.f.o pfft2mod.f.o pfield2mod.f.o �pdfield2mod.f.o pbfield2mod.f.o pcfield2mod.f.o pnfield2mod.f.o pnpfield2mod.f.o �pdiag2mod.f.o mp2mod.f.o mp0mod.f.o p2mod.f.o p0mod.f.o pinit2lib.f.o mprbpush2lib.f.o �mpbpush2lib.f.o mppush2lib.f.o prbpush2lib.f.o pbpush2lib.f.o ppush2lib.f.o mpfft2lib.f.o �pfft2lib.f.o pfield2lib.f.o pdfield2lib.f.o pbfield2lib.f.o pcfield2lib.f.o pnfield2lib.f.o �pdiag2lib.f.o mp2lib.f.o p2lib.f.o p0lib.f.o MacMPI_X.f.o MacMP.f.o �plibgks2.f.o plibgks1.f.o libgks2.f.o libgks1.f.o libmcX.f.o# Linkage rulenew_d0_mpbbeps2 �� {OBJS}	f90 {OPTS90} -carbonappl -o new_d0_mpbbeps2 {OBJS} �	"{SharedLibraries}"MPLibrary# Compilation rulesMacMPI_X.f.o � MacMPI_X.f	f77 -O -c MacMPI_X.fMacMP.f.o � MacMP.f	f77 -O -c MacMP.flibmcX.f.o � libmcX.f	f77 {OPTS77} -c libmcX.flibgks1.f.o � libgks1.f	f77 {OPTS77} -c libgks1.flibgks2.f.o � libgks2.f	f77 {OPTS77} -c libgks2.fplibgks1.f.o � plibgks1.f	f77 {OPTS77} -c plibgks1.fplibgks2.f.o � plibgks2.f	f77 {OPTS77} -c plibgks2.fp0lib.f.o � p0lib.f	f90 {OPTS90} -c p0lib.fmp2lib.f.o � mp2lib.f	f90 {OPTS90} -c mp2lib.fp2lib.f.o � p2lib.f	f90 {OPTS90} -c p2lib.fpinit2lib.f.o � pinit2lib.f	f90 {OPTS90} -c pinit2lib.fmprbpush2lib.f.o � mprbpush2lib.f	f90 {OPTS90} -c mprbpush2lib.fmpbpush2lib.f.o � mpbpush2lib.f	f90 {OPTS90} -c mpbpush2lib.fmppush2lib.f.o � mppush2lib.f	f90 {OPTS90} -c mppush2lib.fprbpush2lib.f.o � prbpush2lib.f	f90 {OPTS90} -c prbpush2lib.fpbpush2lib.f.o � pbpush2lib.f	f90 {OPTS90} -c pbpush2lib.fppush2lib.f.o � ppush2lib.f	f90 {OPTS90} -c ppush2lib.fmpfft2lib.f.o � mpfft2lib.f	f90 {OPTS90} -c mpfft2lib.fpfft2lib.f.o � pfft2lib.f	f90 {OPTS90} -c pfft2lib.fpfield2lib.f.o � pfield2lib.f	f90 {OPTS90} -c pfield2lib.fpdfield2lib.f.o � pdfield2lib.f	f90 {OPTS90} -c pdfield2lib.fpbfield2lib.f.o � pbfield2lib.f	f90 {OPTS90} -c pbfield2lib.fpcfield2lib.f.o � pcfield2lib.f	f90 {OPTS90} -c pcfield2lib.fpnfield2lib.f.o � pnfield2lib.f	f90 {OPTS90} -c pnfield2lib.fpdiag2lib.f.o � pdiag2lib.f	f90 {OPTS90} -c pdiag2lib.fglobals.f.o � globals.f	f90 {OPTS90} -c globals.fpinit2mod.f.o � pinit2mod.f globals.f.o	f90 {OPTS90} -c pinit2mod.fppush2mod.f.o � ppush2mod.f p0mod.f.o	f90 {OPTS90} -c ppush2mod.fpbpush2mod.f.o � pbpush2mod.f p0mod.f.o	f90 {OPTS90} -c pbpush2mod.fprbpush2mod.f.o � prbpush2mod.f p0mod.f.o	f90 {OPTS90} -c prbpush2mod.fmppush2mod.f.o � mppush2mod.f ppush2mod.f.o mp0mod.f.o	f90 {OPTS90} -c mppush2mod.fmpbpush2mod.f.o � mpbpush2mod.f pbpush2mod.f.o mp0mod.f.o	f90 {OPTS90} -c mpbpush2mod.fmprbpush2mod.f.o � mprbpush2mod.f prbpush2mod.f.o mp0mod.f.o	f90 {OPTS90} -c mprbpush2mod.fpfft2mod.f.o � pfft2mod.f p0mod.f.o	f90 {OPTS90} -c pfft2mod.fpfield2mod.f.o � pfield2mod.f globals.f.o	f90 {OPTS90} -c pfield2mod.fpdfield2mod.f.o � pdfield2mod.f globals.f.o	f90 {OPTS90} -c pdfield2mod.fpbfield2mod.f.o � pbfield2mod.f globals.f.o	f90 {OPTS90} -c pbfield2mod.fpcfield2mod.f.o � pcfield2mod.f globals.f.o	f90 {OPTS90} -c pcfield2mod.fpnfield2mod.f.o � pnfield2mod.f globals.f.o	f90 {OPTS90} -c pnfield2mod.fpnpfield2mod.f.o � pnpfield2mod.f pfield2mod.f.o pdfield2mod.f.o pbfield2mod.f.o �                   pcfield2mod.f.o pnfield2mod.f.o	f90 {OPTS90} -c pnpfield2mod.fmpfft2mod.f.o � mpfft2mod.f pfft2mod.f.o mp0mod.f.o	f90 {OPTS90} -c mpfft2mod.fpdiag2mod.f.o � pdiag2mod.f pinit2mod.f.o	f90 {OPTS90} -c pdiag2mod.fp0mod.f.o � p0mod.f globals.f.o	f90 {OPTS90} -c p0mod.fp2mod.f.o � p2mod.f p0mod.f.o	f90 {OPTS90} -c p2mod.fmp0mod.f.o � mp0mod.f	f90 {OPTS90} -c mp0mod.fmp2mod.f.o � mp2mod.f p2mod.f.o mp0mod.f.o 	f90 {OPTS90} -c mp2mod.fmpempush2mod.f.o � mpempush2mod.f mprbpush2mod.f.o mpbpush2mod.f.o �                   mppush2mod.f.o mpfft2mod.f.o mp2mod.f.o	f90 {OPTS90} -c mpempush2mod.fnew_d0_pbbeps2.f.o � new_d0_pbbeps2.f mpempush2mod.f.o pnpfield2mod.f.o �                     pdiag2mod.f.o mp2mod.f.o	f90 {OPTS90} -N11 -s -c new_d0_pbbeps2.fclean � Delete -i new_d0_mpbbeps2 {OBJS}