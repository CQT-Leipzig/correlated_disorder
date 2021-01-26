CC = g++
CCFLAGS = -O3 -std=c++11 -Wall -Wpedantic

LDFLAGS = 
LDLIBS = -lfftw3

all: example
example: example_generate_configuration example_1D_correlation

example_generate_configuration:	example_generate_configuration.cpp
	$(CC) $(CCFLAGS) example_generate_configuration.cpp -o example $(LDFLAGS) $(LDLIBS)

example_1D_correlation:	example_1D_correlation.cpp
	$(CC) $(CCFLAGS) example_1D_correlation.cpp -o example_1D $(LDFLAGS) $(LDLIBS)

clean:
	rm -f *.o example
