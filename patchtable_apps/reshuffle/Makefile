CFLAGS=-fmax-errors=3 -std=c++0x -DBUILD_DEBUG=0 -fopenmp -g
LFLAGS=-lopencv_core -lopencv_imgproc -lopencv_highgui
CC=g++
DEPS=main.o reshuffle.o pfm.o patchtable.o timer.o util.o

reshuffle:$(DEPS)
	$(CC) $(DEPS) -o reshuffle $(LFLAGS) $(CFLAGS)

%.o:%.cpp
	$(CC) -c $< $(CFLAGS)

%.o:../../solver/%.cpp
	$(CC) -c $< $(CFLAGS)

%.o:../../patchtable/%.cpp
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f *.o
	rm -f reshuffle
    
