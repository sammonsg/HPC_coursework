all: task1

task1: main.cpp mat_builder.cpp helpers.cpp solve.cpp
	@g++ main.cpp mat_builder.cpp helpers.cpp solve.cpp \
	-L/usr/local/lib -lclapack -lcblas
	@./a.out 10000 4 12000 172800000 210e9
	#python makemycharts.py

clean:
	@rm -rf *.out

sayhello:
	cd ../../
