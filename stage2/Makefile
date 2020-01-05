all: display measure

display: display.cu
	nvcc -ccbin g++ -std=c++11 -m64 -lGL -lGLU -lglut -o display.out display.cu

measure: measure.cu
	nvcc -ccbin g++ -std=c++11 -m64 -o measure.out measure.cu
