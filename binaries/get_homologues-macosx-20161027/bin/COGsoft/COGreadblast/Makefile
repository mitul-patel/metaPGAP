CC=g++
# CCOPTS=-g -Wall -DDEBUG
CCOPTS=-O2

COGreadblast: reader.o blastconvglob.o bc.o blastconv.o cogreadblast.o
	$(CC) $(CCOPTS) -o COGreadblast reader.o blastconvglob.o bc.o blastconv.o cogreadblast.o

cogreadblast.o: cogreadblast.cpp reader.h blastconv.h
	$(CC) $(CCOPTS) -c cogreadblast.cpp

blastconv.o: blastconv.cpp blastconv.h blastconvglob.h 
	$(CC) $(CCOPTS) -c blastconv.cpp

reader.o: reader.cpp reader.h
	$(CC) $(CCOPTS) -c reader.cpp

bc.o: bc.cpp bc.h
	$(CC) $(CCOPTS) -c bc.cpp

blastconvglob.o: blastconvglob.cpp blastconvglob.h
	$(CC) $(CCOPTS) -c blastconvglob.cpp

clean:
	rm *.o COGreadblast
