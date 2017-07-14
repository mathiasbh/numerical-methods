CFLAGS  = -std=c99 -Wall -O
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
MAINS = mainA mainB mainC
OBJECTS = mainA.o mainB.o mainC.o qrgs.o qrgivens.o

all: output.dat

output.dat : $(MAINS) Makefile
	./$(word 1,$^) > $@
	./$(word 2,$^) >> $@
	./$(word 3,$^) >> $@

mainA : mainA.o qrgs.o print_functions.o
mainB : mainB.o qrgivens.o print_functions.o
mainC : mainC.o qrLU.o print_functions.o


clean :
	rm -rf $(OBJECTS) $(MAINS) *~ *.o output.dat
