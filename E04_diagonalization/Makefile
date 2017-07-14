CFLAGS  = -std=c99 -Wall -O3
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
MAINS   = main
OBJECTS = main.o jacobi.o print_functions.o
HEADER  = jacobi.h

all: Output.txt

Output.txt : $(MAINS) Makefile
	./$< > $@


main : $(OBJECTS)

main.o : $(HEADER)

jacobi.o : $(HEADER)


clean :
	$(RM) $(OBJECTS) $(MAINS) *~ *.o Output.txt
