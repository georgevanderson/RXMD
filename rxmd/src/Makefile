include ../Makefile.inc

$(EXE): $(OBJS)
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS) module.o
	mv $(EXE) ../$(EXE)

##FIXME##
#main.F90 : module.o cg.o
main.F90 : module.o bo.o pot.o support.o init.o

bo.o : module.o
fileio.o : module.o bo.o qeq_vars.o support.o
init.o : module.o param.o pot.o qeq_vars.o qeq.o
pot.o : module.o param.o qeq_vars.o qeq.o bo.o memory.o support.o list.o comm.o
qeq.o : module.o list.o comm.o qeq_vars.o
stress.o : module.o comm.o
#cg.o : module.o pot.o  #FIXME#
comm.o : module.o support.o qeq_vars.o
module.o : memory.o
list.o : module.o param.o support.o
support.o : module.o param.o 


clean:
	rm -f PI* *.o *.mod *.MOD mpif.h ../rxmd
