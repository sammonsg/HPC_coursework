all: task1

task1: main.cpp mat_builder.cpp helpers.cpp
	@g++ main.cpp mat_builder.cpp helpers.cpp
	@./a.out 10000 6 12000 172800000 210e9
	#python makemycharts.py

clean:
	@rm -rf *.out

sayhello:
	cd ../../
