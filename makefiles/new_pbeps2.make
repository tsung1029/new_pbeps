# Makefile for new_pbeps2 with MacOS 9OPTS90 = -O -N113OPTS77 = -O -N113 -N2OBJS = new_pbeps2.f.o �globals.f.o pinit2mod.f.o pespush2mod.f.o prbpush2mod.f.o ppush2mod.f.o �pfft2mod.f.o pfield2mod.f.o pdiag2mod.f.o p2mod.f.o p0mod.f.o mp0mod.f.o �pinit2lib.f.o prbpush2lib.f.o ppush2lib.f.o pfft2lib.f.o pfield2lib.f.o �pdiag2lib.f.o p2lib.f.o p0lib.f.o MacMPI_X.f.o nullMP.f.o �plibgks2.f.o plibgks1.f.o libgks2.f.o libgks1.f.o libmcX.f.o# Linkage rulenew_pbeps2 �� {OBJS}	f90 {OPTS90} -carbonappl -o new_pbeps2 {OBJS}# Compilation rulesMacMPI_X.f.o � MacMPI_X.f	f77 -O -c MacMPI_X.fnullMP.f.o � nullMP.f	f77 -O -c nullMP.flibmcX.f.o � libmcX.f	f77 {OPTS77} -c libmcX.flibgks1.f.o � libgks1.f	f77 {OPTS77} -c libgks1.flibgks2.f.o � libgks2.f	f77 {OPTS77} -c libgks2.fplibgks1.f.o � plibgks1.f	f77 {OPTS77} -c plibgks1.fplibgks2.f.o � plibgks2.f	f77 {OPTS77} -c plibgks2.fp0lib.f.o � p0lib.f	f90 {OPTS90} -c p0lib.fp2lib.f.o � p2lib.f	f90 {OPTS90} -c p2lib.fpinit2lib.f.o � pinit2lib.f	f90 {OPTS90} -c pinit2lib.fppush2lib.f.o � ppush2lib.f	f90 {OPTS90} -c ppush2lib.fprbpush2lib.f.o � prbpush2lib.f	f90 {OPTS90} -c prbpush2lib.fpfft2lib.f.o � pfft2lib.f	f90 {OPTS90} -c pfft2lib.fpfield2lib.f.o � pfield2lib.f	f90 {OPTS90} -c pfield2lib.fpdiag2lib.f.o � pdiag2lib.f	f90 {OPTS90} -c pdiag2lib.fglobals.f.o � globals.f	f90 {OPTS90} -c globals.fpinit2mod.f.o � pinit2mod.f globals.f.o	f90 {OPTS90} -c pinit2mod.fppush2mod.f.o � ppush2mod.f p0mod.f.o	f90 {OPTS90} -c ppush2mod.fprbpush2mod.f.o � prbpush2mod.f p0mod.f.o	f90 {OPTS90} -c prbpush2mod.fpfft2mod.f.o � pfft2mod.f p0mod.f.o	f90 {OPTS90} -c pfft2mod.fpfield2mod.f.o � pfield2mod.f globals.f.o	f90 {OPTS90} -c pfield2mod.fpdiag2mod.f.o � pdiag2mod.f pinit2mod.f.o	f90 {OPTS90} -c pdiag2mod.fp0mod.f.o � p0mod.f globals.f.o	f90 {OPTS90} -c p0mod.fp2mod.f.o � p2mod.f p0mod.f.o	f90 {OPTS90} -c p2mod.fpespush2mod.f.o � pespush2mod.f prbpush2mod.f.o ppush2mod.f.o pfft2mod.f.o �                  p2mod.f.o  	f90 {OPTS90} -c pespush2mod.fmp0mod.f.o � mp0mod.f	f90 {OPTS90} -c mp0mod.fnew_pbeps2.f.o � new_pbeps2.f pespush2mod.f.o pfield2mod.f.o pdiag2mod.f.o �                 mp0mod.f.o	f90 {OPTS90} -s -c new_pbeps2.fclean � Delete -i new_pbeps2 {OBJS}