CFLAGS = -Wall -O

all : A B

A : out.A.txt
# target = out.A.txt depends on "mainA"
# the command is ./mainA > out.A.txt
# meaning, run program mainA and send 
# output to out.A.txt
out.A.txt : mainA ; ./$< > $@
mainA : mainA.c


B : out.B.txt
out.B.txt : mainB ; ./$< > $@
mainB : mainB.c hello.o user.o

# Exercise C: Check if outputs are correct
check : checkA checkB

checkA : out.A.txt
	@echo "check of output A:"
	@echo "hello, mathias" > correct.txt
	@diff --brief correct.txt $< # brief: report when different
	@echo "Output A is OK."
	@rm -f correct.txt

checkB : out.B.txt
	@echo "check of output B:"
	@echo "hello, mathias" > correct.txt
	@diff --brief correct.txt $< # brief: report when different
	@echo "Output B is OK."
	@rm -f correct.txt

clean : cleanA cleanB
cleanA : ; rm -f mainA out.A.txt
cleanB : ; rm -f mainB out.B.txt hello.o user.o 
