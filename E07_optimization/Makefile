CFLAGS  = -std=c99 -Wall -O3
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
MAINS   = main
OBJECTS = main.o newtonMinimize.o qrgivens.o testfunc.o
HEADER  = header.h
OUTPUT  = OutputAB.txt

all: $(OUTPUT)

$(OUTPUT) : $(MAINS) Makefile
	./$< > $@
	cat $@


main : $(OBJECTS)

main.o : $(HEADER)

clean :
	$(RM) $(OBJECTS) $(MAINS) *~ *.o $(OUTPUT)
