include ../Makefile.inc

$(EXE): $(OBJS)
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS) module.o
	mv $(EXE) ../$(EXE)

#$(OBJS) : base.mod

#base.mod: module.F90
#	$(FC) -c $(CPP) $(OPT) module.F90

##FIXME##
#main.F90 : module.o cg.o
main.F90 : module.o bo.o pot.o

bo.o : module.o
fileio.o : module.o bo.o qeq.o
init.o : param.o module.o pot.o
pot.o : param.o module.o bo.o memory.o
qeq.o : module.o
stress.o : module.o
#cg.o : module.o pot.o
comm.o : module.o qeq.o
module.o : memory.o

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h ../rxmd
