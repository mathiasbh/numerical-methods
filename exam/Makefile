CC = gcc
CFLAGS  = -std=c99 -Wall -O3
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
OBJECTS = qrgivens.o print_functions.o inviter.o
HEADER  = header.h
MAINS   = main
PLOT    = OutputC_steps.pdf
OUT     = OutputAB.txt OutputC.txt $(PLOT)
ARCHIVE = mbh-exam.tgz

all: $(OUT) $(ARCHIVE)

main : main.o $(OBJECTS) $(HEADER)

OutputAB.txt : main Makefile ; ./$< 1> $@ 2> OutputC.txt

$(PLOT) : ploterr.pyxplot Makefile ; cat $< | pyxplot

ploterr.pyxplot : OutputC.txt Makefile
	echo 'set terminal pdf' > $@
	echo 'set output "$(PLOT)"' >> $@
	echo 'set ylabel "steps"' >> $@
	echo 'set xlabel "shift away from true eigenvalue"' >> $@
	echo 'plot \ ' >> $@
	echo '[][] "$<" using 1:2 with lines title "Inverse iteration" \ ' >> $@
	echo ',"" using 1:3 with lines title "Rayleigh iteration" \ ' >> $@


$(ARCHIVE) : main main.c main.o Makefile header.h $(OBJECTS) $(OUT) ploterr.pyxplot README.txt
	tar --file $@ --create --exclude=*.tgz $^
	tar --file $@ --list


.PHONY : clean

clean: ; $(RM) $(MAINS) *~ $(OUT) *.pyxplot inviter.o main.o

