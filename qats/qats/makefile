#GCC = g++-4.2 -O3 -fopenmp
GCC = g++ -O3

all: call_qats call_qats_indices

call_qats: call_qats.cpp qats.o
	$(GCC) -o call_qats call_qats.cpp qats.o

call_qats_indices: call_qats_indices.cpp qats.o
	$(GCC) -o call_qats_indices call_qats_indices.cpp qats.o

qats.o: qats.cpp qats.h
	$(GCC) -c qats.cpp