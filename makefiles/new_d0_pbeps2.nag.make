# Makefile for new_d0_pbeps2 with MacOS X using Nag compilerFC = f95OPTS90 = -O3 -r8OBJS = new_d0_pbeps2.o \globals.o pinit2mod.o prbpush2mod.o ppush2mod.o pfield2mod.o \pdfield2mod.o pbfield2mod.o pcfield2mod.o pnpfield2mod.o \pdiag2mod.o p2mod.o p0mod.o pinit2lib.o prbpush2lib.o ppush2lib.o \pfield2lib.o pdfield2lib.o pbfield2lib.o pcfield2lib.o pdiag2lib.o \p2lib.o p0lib.o MacMPIf77.o MacMPI_S.o nullpgks2.o nullpgks1.o# Linkage rulenew_d0_pbeps2.out : $(OBJS)	$(FC) $(OPTS90) -framework carbon -o new_d0_pbeps2.out $(OBJS)# Compilation rulesMacMPI_S.o : MacMPI_S.c	gcc -O -c -I /Developer/Headers/FlatCarbon MacMPI_S.cMacMPIf77.o : MacMPIf77.c	gcc -O -c MacMPIf77.cnullpgks1.o : nullpgks1.f	$(FC) $(OPTS90) -c nullpgks1.fnullpgks2.o : nullpgks2.f	$(FC) $(OPTS90) -c nullpgks2.fp0lib.o : p0lib.f	$(FC) $(OPTS90) -dusty -c p0lib.fp2lib.o : p2lib.f	$(FC) $(OPTS90) -dusty -c p2lib.fpinit2lib.o : pinit2lib.f	$(FC) $(OPTS90) -dusty -c pinit2lib.fppush2lib.o : ppush2lib.f	$(FC) $(OPTS90) -c ppush2lib.fprbpush2lib.o : prbpush2lib.f	$(FC) $(OPTS90) -c prbpush2lib.fpfield2lib.o : pfield2lib.f	$(FC) $(OPTS90) -dusty -c pfield2lib.fpdfield2lib.o : pdfield2lib.f	$(FC) $(OPTS90) -c pdfield2lib.fpbfield2lib.o : pbfield2lib.f	$(FC) $(OPTS90) -c pbfield2lib.fpcfield2lib.o : pcfield2lib.f	$(FC) $(OPTS90) -c pcfield2lib.fpdiag2lib.o : pdiag2lib.f	$(FC) $(OPTS90) -c pdiag2lib.fglobals.o : globals.f	$(FC) $(OPTS90) -c globals.fpinit2mod.o : pinit2mod.f globals.o	$(FC) $(OPTS90) -c pinit2mod.fppush2mod.o : ppush2mod.f p0mod.o	$(FC) $(OPTS90) -c ppush2mod.fprbpush2mod.o : prbpush2mod.f p0mod.o	$(FC) $(OPTS90) -c prbpush2mod.fpfield2mod.o : pfield2mod.f globals.o	$(FC) $(OPTS90) -c pfield2mod.fpdfield2mod.o : pdfield2mod.f globals.o	$(FC) $(OPTS90) -c pdfield2mod.fpbfield2mod.o : pbfield2mod.f globals.o	$(FC) $(OPTS90) -c pbfield2mod.fpcfield2mod.o : pcfield2mod.f globals.o	$(FC) $(OPTS90) -c pcfield2mod.fpnpfield2mod.o : pnpfield2mod.f pfield2mod.o pdfield2mod.o pbfield2mod.o \                 pcfield2mod.o	$(FC) $(OPTS90) -c pnpfield2mod.fpdiag2mod.o : pdiag2mod.f pinit2mod.o	$(FC) $(OPTS90) -c pdiag2mod.fp0mod.o : p0mod.f globals.o	$(FC) $(OPTS90) -c p0mod.fp2mod.o : p2mod.f p0mod.o	$(FC) $(OPTS90) -dusty -c p2mod.fnew_d0_pbeps2.o : new_d0_pbeps2.f ppush2mod.o prbpush2mod.o \                  pnpfield2mod.o pdiag2mod.o p2mod.o p0mod.o	$(FC) $(OPTS90) -c new_d0_pbeps2.fclean :# Delete -i new_d0_pbeps2.out $(OBJS)