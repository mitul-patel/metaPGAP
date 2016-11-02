CC=g++
#CCOPTS=-g -Wall -DDEBUG
CCOPTS=-O2

COGlse: mklse.o reader.o
	$(CC) $(CCOPTS) -o COGlse mklse.o reader.o

mklse.o: mklse.cpp graph.h logger.h 
	$(CC) $(CCOPTS) -c mklse.cpp

reader.o: reader.cpp reader.h
	$(CC) $(CCOPTS) -c reader.cpp

clean:
	rm *.o COGlse
