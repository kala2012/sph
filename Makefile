CC=icc
CFLAGS=-g3 -D_USE_MATH_DEFINES `pkg-config --cflags --libs mkl` -O2
OBJ= $(SRC:.c=.o)
OBJ_FORTRAN= $(SRC_FORTRAN:.f=.of)
F77=gfortran

EXEC=sph
SRC=derivatives.c inoutflow.c  inputbp.c  input.c  output.c  parameterfile.c  sph.c System.c
SRC_FORTRAN=derivatives.f inoutflow.f  inputbp.f  input.f  output.f  parameterfile.f  sph.f 

all: $(EXEC)

sph: $(OBJ)
	${CC} -Wall -std=c99 -g -o $@ $^ $(CFLAGS) -lm -limf 

sphf90: $(OBJ_FORTRAN)
	$(F77)  -g -o $@ $^ $(CFLAGS) -lm -limf 
%.o: %.c
	${CC} -Wall -g -std=c99 -o $@ -c $< $(CFLAGS) 

%.of: %.f
	${F77} -Wall -g -o $@ -c $< ${FLAGS}

.PHONY: clean mrproper

clean: 
	rm -rf *.o

mrproper: clean
	rm -Rf $(EXEC)
