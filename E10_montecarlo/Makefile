CC = gcc
CFLAGS  = -std=c99 -Wall -O3
#CFLAGS += $$(gsl-config --cflags)
#LDLIBS += $$(gsl-config --libs)
LDLIBS  = -lm
OUTPUTA = OutputErrConv.pdf
OBJECTS = mcplain.o integrator.o mctol.o plainrand.o
HEADER  = header.h
MAINS   = mainA mainB
OUT     = OutputA.txt OutputA.log.txt OutputB.log.txt $(OUTPUTA)

all: $(OUT)

mainA : mainA.o $(OBJECTS) $(HEADER) 
mainB : mainB.o $(OBJECTS) $(HEADER)

OutputA.txt : mainA Makefile ; ./$< 1> $@ 2> OutputA.log.txt
OutputB.log.txt : mainB Makefile ; ./$< 2> $@

$(OUTPUTA) : ploterr.pyxplot Makefile ; cat $< | pyxplot

ploterr.pyxplot : OutputA.txt Makefile
	echo 'set terminal pdf' > $@
	echo 'set output "$(OUTPUTA)"' >> $@
	echo 'set xlabel "$$N$$"' >> $@
	echo 'set ylabel "error"' >> $@
	echo 'plot \ ' >> $@
	echo '[][] "$<" using 1:2 with lines notitle \ ' >> $@
	echo ', "" using 1:3.3/(sqrt($$1)) with lines title "$$1/\sqrt{N}$$" \ ' >> $@


clean: ; $(RM) $(MAINS) *.o *~ $(OUT) *.pyxplot

