CFLAGS  = -std=c99 -Wall -O3
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
MAINS   = main
OBJECTS = main.o newton.o qrgivens.o newtonquadratic.o testfunc.o
HEADER  = newton.h

all: Output.ABC.txt

Output.ABC.txt : $(MAINS) Makefile
	./$< > $@
	cat $@


main : $(OBJECTS)

main.o : $(HEADER)


clean :
	$(RM) $(OBJECTS) $(MAINS) *~ *.o Output.ABC.txt
