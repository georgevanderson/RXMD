include ../Makefile.inc

$(EXE): $(OBJS)
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS) module.o
	mv $(EXE) ../$(EXE)

$(OBJS) : base.mod

base.mod: module.F90
	$(FC) -c $(CPP) $(OPT) module.F90

main.F90 : cg.o

bo.o : param.o
fileio.o : param.o
init.o : param.o
pot.o : param.o
qeq.o : param.o
stress.o : param.o

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h ../rxmd
