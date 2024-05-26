PHONY: main

main: main.cpp
	g++ -o main -g -pg -std=c++17 -Wall -ggdb3 main.cpp 

clean:
	rm -rf main build