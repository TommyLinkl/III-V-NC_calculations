OBJS = main.o fit.o init.o read.o pot.o states.o hamiltonian.o diag.o write.o vector.o rand.o complex.o

# Tells the program to use the intel c compiler - faster than the gnu compiler 
CC = icx

MKLROOT = /opt/intel/oneapi/mkl/latest
# Gives the path to the intel compilers needed in this program
LIB =  -lpthread -lfftw3 ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
CFLAGS = -fast -Wall -Wextra -std=c99 -DMKL_ILP64 -I"${MKLROOT}/include" -qopenmp 
# CFLAGS = -O0 -g -std=c99 -DMKL_ILP64 -I"${MKLROOT}/include"

fit: $(OBJS)
	$(CC) $(CFLAGS) -o soFit.x $(OBJS) $(LIB)
clean:
	/bin/rm *.o *.x

