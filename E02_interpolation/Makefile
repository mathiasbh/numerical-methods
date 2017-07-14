CFLAGS = -std=c99 -Wall -O
LDLIBS = -lm
PLOT_INTERPOLATE = Plot_interpolate.pdf
PLOT_DERIVATIVE = Plot_derivative.pdf 
PLOT_INTEGRAL = Plot_integral.pdf
DATA = data.txt

all : $(DATA) $(PLOT_INTERPOLATE) $(PLOT_DERIVATIVE) $(PLOT_INTEGRAL)
$(PLOT_INTERPOLATE) : PyxInterpolation.pyxplot Makefile ; cat $< | pyxplot
$(PLOT_DERIVATIVE) : PyxDerivative.pyxplot Makefile ; cat $< | pyxplot
$(PLOT_INTEGRAL) : PyxIntegral.pyxplot Makefile ; cat $< | pyxplot

data.txt : main ; ./$< > $@
main : linterp.o qspline.o cspline.o bsearch.o

# Plot of data using pyxplot
PyxInterpolation.pyxplot : $(DATA) Makefile
	echo 'set terminal pdf' > $@
	echo 'set output "$(PLOT_INTERPOLATE)"' >> $@
	echo 'set xlabel "$$x$$"' >> $@
	echo 'set ylabel "$$y$$"' >> $@
	echo 'plot \ ' >> $@
	echo '[][] "$(word 1,$^)" index 0 with points title "Data" \ ' >> $@
	echo ',    "$(word 1,$^)" index 1 using 1:2 with lines color black title "$$e^{-x^2}$$" \ ' >> $@
	echo ',    "$(word 1,$^)" index 2 with lines title "Linear" \ ' >> $@
	echo ',    "$(word 1,$^)" index 3 using 1:2 with lines title "Quadratic" \ ' >> $@
	echo ',    "$(word 1,$^)" index 4 using 1:2 with lines title "Cubic" \ ' >> $@


# Plot of data using pyxplot
PyxDerivative.pyxplot : $(DATA) Makefile
	echo 'set terminal pdf' > $@
	echo 'set output "$(PLOT_DERIVATIVE)"' >> $@
	echo 'set xlabel "$$x$$"' >> $@
	echo 'set ylabel "$$y$$"' >> $@
	echo 'plot \ ' >> $@
	echo '[][] "$(word 1,$^)" index 1 using 1:3 with lines title "$$-2xe^{-x^2}$$" \ ' >> $@
	echo ',    "$(word 1,$^)" index 3 using 1:3 with lines title "Quadratic" \ ' >> $@
	echo ',    "$(word 1,$^)" index 4 using 1:3 with lines title "Cubic" \ ' >> $@


# Plot of data using pyxplot
PyxIntegral.pyxplot : $(DATA) Makefile
	echo 'set terminal pdf' > $@
	echo 'set output "$(PLOT_INTEGRAL)"' >> $@
	echo 'set xlabel "$$x$$"' >> $@
	echo 'set ylabel "$$y$$"' >> $@
	echo 'set key bottom right' >> $@
	echo 'plot \ ' >> $@
	echo '[][] "$(word 1,$^)" index 1 using 1:4 with lines title "$$\\frac{\\sqrt{\\pi}}{2} \\textup{erf}(x)$$" \ ' >> $@
	echo ',    "$(word 1,$^)" index 3 using 1:4 with lines title "Quadratic" \ ' >> $@
	echo ',    "$(word 1,$^)" index 4 using 1:4 with lines title "Cubic" \ ' >> $@


clean :
	rm -rf main $(PLOT_INTERPOLATE) $(PLOT_DERIVATIVE) $(PLOT_INTEGRAL) $(DATA) *~ *.pyxplot cspline.o qspline.o bsearch.o linterp.o 
