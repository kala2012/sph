CC=icc
CFLAGS=-g -D_USE_MATH_DEFINES
OBJ= $(SRC:.c=.o)
EXEC=sph
SRC=derivatives.c inoutflow.c  inputbp.c  input.c  output.c  parameterfile.c  sph.c System.c

all: $(EXEC)

sph: $(OBJ)
	${CC} -Wall -std=c99 -g  -o $@ $^ $(CFLAGS) -lm

%.o: %.c
	${CC} -Wall -g -std=c99 -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean: 
	rm -rf *.o

mrproper: clean
	rm -Rf $(EXEC)
