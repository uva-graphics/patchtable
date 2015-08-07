CFLAGS_BASE=-fopenmp -funroll-loops -fmax-errors=3 -std=c++0x -fabi-version=6 -fomit-frame-pointer -msse4.2
CFLAGS_RELEASE=$(CFLAGS_BASE) -O3 -DBUILD_DEBUG=0
CFLAGS_DEBUG=$(CFLAGS_BASE) -g -DBUILD_DEBUG=1
LFLAGS=-lpng -ldl -lopencv_core -lopencv_imgproc -lopencv_flann -lopencv_video -lANN
CC=g++
DEPS=patchtable.o patchtable_main.o array.o timer.o util.o pfm.o analogy.o
H_DEPS=gck.h patchtable.h analogy.h

CFLAGS=$(CFLAGS_RELEASE)
debug: CFLAGS = $(CFLAGS_DEBUG)

patchtable: $(DEPS) $(H_DEPS)
	$(CC) -o patchtable $(CFLAGS) $(DEPS) $(LFLAGS)

%.o:%.cpp
	$(CC) -c $< $(CFLAGS)

%.o:../solver/%.cpp
	$(CC) -c $< $(CFLAGS)

gck.h: gck.py
	python gck.py

debug: $(DEPS) $(H_DEPS)
	python gck.py
	$(CC) -o patchtable $(CFLAGS) $(DEPS) $(LFLAGS)

clean:
	rm -f *.o
	rm -f patchtable

