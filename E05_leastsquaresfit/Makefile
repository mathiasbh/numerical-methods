CFLAGS  = -std=c99 -Wall -O
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
OBJECTS = main.o qrgivens.o lsfit.o jacobi.o lssvd.o
HEADERS = header.h 
DATA    = Output.ABC.txt

all: $(DATA) PlotA.pdf PlotB.pdf PlotC.pdf

$(DATA) : main Makefile
	./$< > $@ 2> Output.Console.txt

main : $(OBJECTS)

PlotA.pdf : plotA.pyxplot Makefile ; cat $< | pyxplot
PlotB.pdf : plotB.pyxplot Makefile ; cat $< | pyxplot
PlotC.pdf : plotC.pyxplot Makefile ; cat $< | pyxplot

plotA.pyxplot : $(DATA) Makefile
	echo 'set terminal pdf monochrome' > $@
	echo 'set output "PlotA.pdf"' >> $@
	echo 'set xlabel "$$x$$"' >> $@
	echo 'set ylabel "$$y$$"' >> $@
	echo 'plot \ ' >> $@
	echo '[0:3][4:18] "$<" index 0 with yerrorbars title "Data" \ ' >> $@
	echo ',      "$<" index 1 using 1:2 with lines linetype 1 title "Fit A" \ ' >> $@
	echo ',      "$<" index 1 using 1:3 with lines linetype 2 notitle \ ' >> $@
	echo ',      "$<" index 1 using 1:4 with lines linetype 2 notitle \ ' >> $@

plotB.pyxplot : $(DATA) Makefile
	echo 'set terminal pdf monochrome' > $@
	echo 'set output "PlotB.pdf"' >> $@
	echo 'set xlabel "$$x$$"' >> $@
	echo 'set ylabel "$$y$$"' >> $@
	echo 'set key left' >> $@
	echo 'plot \ ' >> $@
	echo '[][] "$<" index 2 with yerrorbars title "Data" \ ' >> $@
	echo ',      "$<" index 3 using 1:2 with lines linetype 1 title "Fit B" \ ' >> $@
	echo ',      "$<" index 3 using 1:3 with lines linetype 2 notitle \ ' >> $@
	echo ',      "$<" index 3 using 1:4 with lines linetype 2 notitle \ ' >> $@


plotC.pyxplot : $(DATA) Makefile
	echo 'set terminal pdf monochrome' > $@
	echo 'set output "PlotC.pdf"' >> $@
	echo 'set xlabel "$$x$$"' >> $@
	echo 'set ylabel "$$y$$"' >> $@
	echo 'plot \ ' >> $@
	echo '[0:3][4:18] "$<" index 0 with yerrorbars title "Data" \ ' >> $@
	echo ',      "$<" index 4 using 1:2 with lines linetype 1 title "Fit C" \ ' >> $@
	echo ',      "$<" index 4 using 1:3 with lines linetype 2 notitle \ ' >> $@
	echo ',      "$<" index 4 using 1:4 with lines linetype 2 notitle \ ' >> $@

clean : 
	$(RM) $(DATA) main $(OBJECTS) *~ *.pyxplot PlotA.pdf PlotB.pdf PlotC.pdf

