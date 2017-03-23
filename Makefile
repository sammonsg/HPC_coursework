compile:
	@g++ -std=c++11 main.cpp mat_builder.cpp helpers.cpp solve.cpp \
	-L/usr/local/lib -llapack -lblas

compile_mpi:
	@mpicxx -std=c++11 main.cpp mat_builder.cpp helpers.cpp solve.cpp solve_parallel.cpp \
	-L/usr/local/lib -llapack -lblas

run1:
	@./a.out 10000 24 12000 14400000 210000 1
	# Do things

run2:
	@./a.out 10000 24 12000 14400000 210000 2 1 20000 7850
	# Do things

run3:
	@./a.out 10000 24 12000 14400000 210000 3 1 20000 7850
	# Do things

run4:
	@mpiexec -np 2 ./a.out 10000 24 12000 14400000 210000 4 1 20000 7850
	# Do things

run5:
	@mpiexec -np 4 ./a.out 10000 8 12000 14400000 210000 5 1 20000 7850
	# Do things

clean:
	@rm -rf *.out

task1: compile_mpi run1 clean
task2: compile_mpi run2 clean
task3: compile_mpi run3 clean
task4: compile_mpi run4 clean
task5: compile_mpi run5 clean

all: compile_mpi run1 run2 run3 run4 run5
