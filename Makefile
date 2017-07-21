CC=g++ -O3 -std=c++11 -Wall -Wpedantic

all: example
example: example_generate_configuration

example_generate_configuration:	example_generate_configuration.cpp
	$(CC) example_generate_configuration.cpp -o example
