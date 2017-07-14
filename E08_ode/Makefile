CC = gcc
CFLAGS  = -std=c99 -Wall -O3
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
OUTPUTDA = OutputA.txt
OUTPUTDB = OutputB.txt
OUTPUTPB1 = OutputB1.pdf
OUTPUTPB2 = OutputB2.pdf
OUTPUTDC = OutputC.txt
OUTPUTPC = OutputC.pdf
OBJECTS = rkstep12.o rkstep45.o rkstep2.o rkstep5.o driver.o
HEADER  = header.h
MAINS   = mainA mainB mainC
OUT     = $(OUTPUTDA) $(OUTPUTDB) $(OUTPUTDC) $(OUTPUTPB1) $(OUTPUTPB2) $(OUTPUTPC)

all: $(OUT)
mainA : mainA.o $(OBJECTS) $(HEADER) 
mainB : mainB.o $(OBJECTS) $(HEADER)
mainC : mainC.o $(OBJECTS) $(HEADER)
$(OUTPUTDA) : mainA Makefile ; ./$< > $@
$(OUTPUTDB) : mainB Makefile ; ./$< > $@
$(OUTPUTDC) : mainC Makefile ; ./$< > $@ 2> OutputC_comparison.txt


$(OUTPUTPB1) : plotB1.pyxplot Makefile ; cat $< | pyxplot
$(OUTPUTPB2) : plotB2.pyxplot Makefile ; cat $< | pyxplot
$(OUTPUTPC) : plotC.pyxplot Makefile ; cat $< | pyxplot

plotB1.pyxplot : $(OUTPUTDB) Makefile
	echo 'set terminal pdf' > $@
	echo 'set output "$(OUTPUTPB1)"' >> $@
	echo 'set xlabel "$$x$$"' >> $@
	echo 'set ylabel "$$y$$"' >> $@
	echo 'plot \ ' >> $@
	echo '[][] "$<" index 0 using 1:2 with lines notitle \ ' >> $@
	echo ',      "" index 0 using 1:3 with lines notitle \ ' >> $@


plotB2.pyxplot : $(OUTPUTDB) Makefile
	echo 'set terminal pdf' > $@
	echo 'set output "$(OUTPUTPB2)"' >> $@
	echo 'set xlabel "$$x$$"' >> $@
	echo 'set ylabel "$$y$$"' >> $@
	echo 'plot \ ' >> $@
	echo '[][] "$<" index 1 using 1:2 with lines notitle \ ' >> $@
	echo ',      "" index 1 using 1:3 with lines notitle \ ' >> $@


plotC.pyxplot : $(OUTPUTDC) Makefile
	echo 'set terminal pdf' > $@
	echo 'set output "$(OUTPUTPC)"' >> $@
	echo 'set xlabel "$$x$$"' >> $@
	echo 'set ylabel "$$y$$"' >> $@
	echo 'plot \ ' >> $@
	echo '[][] "$<" index 0 using 1:2 with lines title "RK12" \ ' >> $@
	echo ',      "" index 1 using 1:2 with lines title "RK45" \ ' >> $@
	echo ',      "" index 2 using 1:2 with lines title "RK2" \ ' >> $@
	echo ',      "" index 3 using 1:2 with lines title "RK5" \ ' >> $@

clean :
	$(RM) *.o $(MAINS) *~ $(OUT) *.pyxplot OutputC_comparison.txt

