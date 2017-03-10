all:
	compile
	task1
	task2
	task3
	task4
	task5

task1: compile
	@./a.out 10000 4 12000 172800000 210000 1
	# Do things

task2: compile
	@./a.out 10000 4 12000 172800000 210000 2 1 10 7850
	# Do things

task3:
	compile
	@./a.out 10000 4 12000 172800000 210000 3 1 10 7850
	# Do things

task4:
	compile
	@./a.out 10000 4 12000 172800000 210000 4
	# Do things

task5:
	compile
	@./a.out 10000 4 12000 172800000 210000 5
	# Do things

clean:
	@rm -rf *.out


compile:
	@g++ -std=c++11 main.cpp mat_builder.cpp helpers.cpp solve.cpp \
	-L/usr/local/lib -llapack -lblas
