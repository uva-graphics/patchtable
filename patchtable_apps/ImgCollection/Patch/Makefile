CFLAGS_BASE=-fmax-errors=3 -std=c++11 -fopenmp -g -I../../../patchtable -I../../../solver
CFLAGS_RELEASE=$(CFLAGS_BASE) -O3 -DBUILD_DEBUG=0
CFLAGS_DEBUG=$(CFLAGS_BASE) -g -DBUILD_DEBUG=0
LFLAGS=-lpng -lopencv_core -lopencv_imgproc -lopencv_highgui
CC=g++
DEPS=main.o PatchMatch.o pfm.o patchtable.o timer.o util.o

CFLAGS=$(CFLAGS_RELEASE)
debug: CFLAGS = $(CFLAGS_DEBUG)

main:$(DEPS)
	$(CC) $(DEPS) -o main $(LFLAGS) $(CFLAGS)

%.o:%.cpp
	$(CC) -c $< $(CFLAGS)

%.o:../../../solver/%.cpp
	$(CC) -c $< $(CFLAGS)

%.o:../../../patchtable/%.cpp
	$(CC) -c $< $(CFLAGS)

debug: $(DEPS)
	$(CC) $(DEPS) -o main $(LFLAGS) $(CFLAGS)
    
clean:
	rm -f *.o
	rm -f main
    
