CC = gcc
CFLAGS  = -std=c99 -Wall -O3
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
OUTPUTA = OutputABC.txt
OBJECTS = integrator.o rkstep5.o driver.o
HEADER  = header.h header1.h
MAINS   = main
OUT     = $(OUTPUTA)

all: $(OUT)
main : main.o $(OBJECTS) $(HEADER) 

$(OUTPUTA) : main Makefile ; ./$< > $@

clean: ; $(RM) $(MAINS) $(OUT) *.o *~

