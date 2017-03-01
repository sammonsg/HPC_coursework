all: task1

task1: main.cpp
	@g++ main.cpp -o code.o
	@./code.o 10000 4 12000 172800000 210000

clean:
	@rm -rf *.o
