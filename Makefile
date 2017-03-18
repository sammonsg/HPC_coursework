compile:
	@g++ -std=c++11 main.cpp mat_builder.cpp helpers.cpp solve.cpp solve_parallel.cpp \
	-L/usr/local/lib -llapack -lblas

compile_mpi:
	@mpicxx main.cpp mat_builder.cpp helpers.cpp solve.cpp solve_parallel.cpp \
	-L/usr/local/lib -llapack -lblas

run1:
	@./a.out 10000 24 12000 14400000 210000 1
	# Do things

run2:
	@./a.out 10000 24 12000 14400000 210000 2 1 100000 7850
	# Do things

run3:
	@./a.out 10000 24 12000 14400000 210000 3 1 100000 7850
	# Do things

run4:
	@mpiexec -np 4 ./a.out 10000 4 12000 14400000 210000 4 1 100000 7850
	# Do things

run5:
	@./a.out 10000 4 12000 14400000 210000 5
	# Do things

clean:
	@rm -rf *.out

task1: compile run1 clean
task2: compile run2 clean
task3: compile run3 clean
task4: compile run4 clean
task5: compile run5 clean

all: compile run1 run2 run3 compile_mpi run4 run5
