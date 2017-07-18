include ../Makefile.inc

$(EXE): $(OBJS)
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS) module.o
	mv $(EXE) ../$(EXE)

#$(OBJS) : base.mod

#base.mod: module.F90
#	$(FC) -c $(CPP) $(OPT) module.F90

main.F90 : module.o cg.o

bo.o : module.o
fileio.o : module.o bo.o
init.o : param.o module.o
pot.o : param.o module.o bo.o
qeq.o : module.o
stress.o : module.o
cg.o : module.o
comm.o : module.o

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h ../rxmd
