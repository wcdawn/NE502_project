# Makefile

FC = gfortran
FFLAGS = -O2 -fbacktrace -fbounds-check

.SUFFIXES :
.SUFFIXES : .f95 .o

default: uncertainty.exe

.f95.o :
	$(FC) $(FFLAGS) -c $*.f95

uncertainty.exe: tableio.o calculations.o uncertainty.o
	$(FC) $(FFLAGS) tableio.o calculations.o uncertainty.o -o $@

clean:
	rm -f core *.o *.mod *.exe

query:
	cvs -n update -P

update:
	cvs update -P

