
FC = gfortran
OMP = -fopenmp

FFLAGS = -O4 $(OMP) -fmax-errors=3 -ffree-line-length-none -fcheck=all -Wall
TARGET = multifract_albero

SOURCES = albero.f90 carepackage.f90

OBJS = $(SOURCES:.f90=.o)


.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o: 
	$(FC) $(FFLAGS) -c  $<

all: $(OBJS)
	$(FC) $(OMP) -o $(TARGET) $(OBJS)

clean:
	rm -f *.o *.mod $(TARGET)
	find -name '*.dat' -delete
	rm -f fit.log
	rm -f *.png

distclean: clean
	rm $(TARGET)


albero.o: carepackage.o
