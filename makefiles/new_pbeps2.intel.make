# Makefile for new_pbeps2 using Intel compilerFC = ifcOPTS90 = -O3 -tpp7 -xW -r8OBJS = new_pbeps2.o \globals.o pinit2mod.o ppush2mod.o prbpush2mod.o pfield2mod.o \pdiag2mod.o p2mod.o p0mod.o pinit2lib.o prbpush2lib.o \ppush2lib.o pfield2lib.o pdiag2lib.o p2lib.o p0lib.o \MacMPIf77.o LnxMPI_S.o nullpgks2.o nullpgks1.o# Linkage rulenew_pbeps2.out : $(OBJS)	$(FC) $(OPTS90) -lm -o new_pbeps2.out $(OBJS)# Compilation rulesLnxMPI_S.o : LnxMPI_S.c	gcc -O -c LnxMPI_S.cMacMPIf77.o : MacMPIf77.c	gcc -O -c MacMPIf77.cnullpgks1.o : nullpgks1.f	$(FC) $(OPTS90) -c nullpgks1.fnullpgks2.o : nullpgks2.f	$(FC) $(OPTS90) -c nullpgks2.fp0lib.o : p0lib.f	$(FC) $(OPTS90) -c p0lib.fp2lib.o : p2lib.f	$(FC) $(OPTS90) -c p2lib.fpinit2lib.o : pinit2lib.f	$(FC) $(OPTS90) -c pinit2lib.fppush2lib.o : ppush2lib.f	$(FC) $(OPTS90) -c ppush2lib.fprbpush2lib.o : prbpush2lib.f	$(FC) $(OPTS90) -c prbpush2lib.fpfield2lib.o : pfield2lib.f	$(FC) $(OPTS90) -c pfield2lib.fpdiag2lib.o : pdiag2lib.f	$(FC) $(OPTS90) -c pdiag2lib.fglobals.o : globals.f	$(FC) $(OPTS90) -c globals.fpinit2mod.o : pinit2mod.f globals.o	$(FC) $(OPTS90) -c pinit2mod.fppush2mod.o : ppush2mod.f p0mod.o	$(FC) $(OPTS90) -c ppush2mod.fprbpush2mod.o : prbpush2mod.f p0mod.o	$(FC) $(OPTS90) -c prbpush2mod.fpfield2mod.o : pfield2mod.f globals.o	$(FC) $(OPTS90) -c pfield2mod.fpdiag2mod.o : pdiag2mod.f pinit2mod.o	$(FC) $(OPTS90) -c pdiag2mod.fp0mod.o : p0mod.f globals.o	$(FC) $(OPTS90) -c p0mod.fp2mod.o : p2mod.f p0mod.o	$(FC) $(OPTS90) -c p2mod.fnew_pbeps2.o : new_pbeps2.f ppush2mod.o prbpush2mod.o pfield2mod.o \               pdiag2mod.o p2mod.o	$(FC) $(OPTS90) -c new_pbeps2.fclean :# Delete -i new_pbeps2.out $(OBJS)