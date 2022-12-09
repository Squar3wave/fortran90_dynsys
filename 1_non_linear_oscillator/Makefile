
FC = gfortran
OMP = -fopenmp

FFLAGS = -Wall -fcheck=all -O3 $(OMP) -fmax-errors=3
TARGET = osc_non_lin

SOURCES = osc_non_lin.f90 runge_kutta_4.f90 carepackage.f90

OBJS = $(SOURCES:.f90=.o)

#LAPACKDIR = ../lapack-3.9.0
#LAPACK = -L$(LAPACKDIR) -llapack -lrefblas

MKLDIR = 
MKLLIB = -L$(MKLDIR) -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5

#LIBS = $(LAPACK)


.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o: 
	$(FC) $(FFLAGS) -c  $<

all: $(OBJS)
	$(FC) $(OMP) -o $(TARGET) $(OBJS) #$(LIBS) 

clean:
	rm -f *.o *.mod *.dat $(TARGET)
	find -name '*.dat' -delete 
	rm -f *.png

distclean: clean
	rm $(TARGET)


osc_non_lin.o: carepackage.o runge_kutta_4.o
runge_kutta_4.o: carepackage.o
