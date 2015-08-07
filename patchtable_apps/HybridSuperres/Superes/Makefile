CFLAGS=-fmax-errors=3 -std=c++0x -DBUILD_DEBUG=0 -fopenmp -g -I../../../patchtable -I../../../solver
LFLAGS=-lopencv_core -lopencv_imgproc -lopencv_highgui
CC=g++
DEPS=superes.o pfm.o patchtable.o timer.o util.o

superes:$(DEPS)
	$(CC) $(DEPS) -o superes $(LFLAGS) $(CFLAGS)

%.o:%.cpp
	$(CC) -c $< $(CFLAGS)

%.o:../../../solver/%.cpp
	$(CC) -c $< $(CFLAGS)

%.o:../../../patchtable/%.cpp
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f *.o
	rm -f superes
    
