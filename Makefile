CC=gcc
CFLAGS=-g -O2 
OBJ= $(SRC:.c=.o)
EXEC=sph
SRC=derivatives.c inoutflow.c  inputbp.c  input.c  output.c  parameterfile.c  sph.c System.c

all: $(EXEC)

sph: $(OBJ)
	gcc -std=c99 -g -o $@ $^ $(CFLAGS) -lm

%.o: %.c
	gcc -g -std=c99 -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean: 
	rm -rf *.o

mrproper: clean
	rm -Rf $(EXEC)
