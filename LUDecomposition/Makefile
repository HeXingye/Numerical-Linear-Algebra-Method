## makefile for problem1
## PROGRAM TO SOLVE AX=B using Gaussian elimination with partial pivoting
OBJECTS = LinAl.o main.o

FSOURCE    = LinAl.f90 main.f90 

MODULES = linal.mod

.PHONY: clean

main.exe: $(MODULES) $(OBJECTS)
	gfortran LinAl.f90 main.f90 -o main
	./main Amat.dat	Bmat.dat
       
%.o: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8  $<

%.mod: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8  $<

clean:
	rm -f ${ALLOBJ} ${EXECBIN} *.o *.exe



