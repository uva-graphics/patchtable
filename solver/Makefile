CXX=g++
CXXFLAGS=-fmax-errors=3

all: CXXFLAGS += -O3 -funroll-loops -DBUILD_DEBUG=0 -I/opt/X11/include -I ../extern/jsoncpp-src-0.5.0/bin/include -L ../extern/jsoncpp-src-0.5.0/bin -I ../../Halide/include -L ../../Halide/bin
debug: CXXFLAGS += -g -DBUILD_DEBUG=1 -I/opt/X11/include -I ../../Halide/include -L ../../Halide/bin
all: LDFLAGS += -Isrc -L/opt/X11/lib -ladept -lX11 -lpthread -lpng -lHalide -ldl -ljson -llbfgs
debug: LDFLAGS += -Isrc -ladept -lX11 -L/opt/X11/lib -lpthread -lpng -lHalide -ldl -ljson -llbfgs

all: experiment_time solve_kernel test
debug: experiment_time solve_kernel test
	
experiment_time: experiment_time.cpp timer.h timer.cpp pfm.h pfm.cpp
	$(CXX) -o experiment_time $(CXXFLAGS) experiment_time.cpp timer.cpp $(LDFLAGS)

solve_kernel: solve_kernel.cpp array.h array.cpp timer.h timer.cpp util.cpp filter.h util.h pareto.h params.h params.cpp solver.cpp solver.h simplify.h mutate.h test_solver.h test_solver.cpp filter.cpp filter_json.cpp filter_json.h mutate.cpp pfm.h pfm.cpp
	$(CXX) -o solve_kernel $(CXXFLAGS) solve_kernel.cpp array.cpp params.cpp timer.cpp util.cpp solver.cpp test_solver.cpp filter.cpp filter_json.cpp mutate.cpp -std=c++0x $(LDFLAGS) -Wall

test: test.cpp array.h timer.h timer.cpp util.cpp filter.h util.h pareto.h pfm.h pfm.cpp
	$(CXX) -o test $(CXXFLAGS) test.cpp -std=c++0x -Wall $(LDFLAGS)

clean:
	rm -f experiment_time solve_kernel test

