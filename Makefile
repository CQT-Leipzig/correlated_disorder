CC = g++
CCFLAGS = -O3 -std=c++11 -Wall -Wpedantic

LDFLAGS = 
LDLIBS = -lfftw3

all: example
example: example_generate_configuration

example_generate_configuration:	example_generate_configuration.cpp
	$(CC) $(CCFLAGS) example_generate_configuration.cpp -o example $(LDFLAGS) $(LDLIBS)

clean:
	rm -f *.o example
