CC=g++
# CCOPTS=-g -Wall -DDEBUG
CCOPTS=-O3 -g -Wall 
#CCOPTS=-O3 -g -Wall -D_GLIBCXX_DEBUG

COGtriangles: cogmaker.o cogtriangles.o
	$(CC) $(CCOPTS) -o COGtriangles cogmaker.o cogtriangles.o

cogtriangles.o: cogtriangles.cpp os.h cogmaker.h
	$(CC) $(CCOPTS) -c cogtriangles.cpp

cogmaker.o: cogmaker.cpp cogmaker.h
	$(CC) $(CCOPTS) -c cogmaker.cpp

clean:
	rm -f *.o COGtriangles
