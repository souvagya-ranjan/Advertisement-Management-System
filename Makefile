CC=mpic++ 
CFLAGS=-std=c++11 -g -fopenmp -lmpi

a3: main.cpp
	$(CC) $(CFLAGS) main.cpp -fopenmp -o a3

run: main.cpp
	mpirun -np 8 ./a3 --taskid=2 --inputpath=A3_test/test1/test-input-1.gra --headerpath=A3_test/test1/test-header-1.dat --outputpath=Output/ayushT2_output1.txt --verbose=1 --startk=1 --endk=3 --p=10

clean:
	rm -rf *.o a3