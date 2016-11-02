CC=g++
#CCOPTS=-g -Wall -DDEBUG
CCOPTS=-O2

COGmakehash: ereader.o COGmakehash.o
	$(CC) $(CCOPTS) -o COGmakehash ereader.o COGmakehash.o

COGmakehash.o: COGmakehash.cpp ereader.h
	$(CC) $(CCOPTS) -c COGmakehash.cpp

ereader.o: ereader.cpp ereader.h
	$(CC) $(CCOPTS) -c ereader.cpp

clean:
	rm *.o COGmakehash
