CC = g++
MPI = mpic++
CFLAGS = -Wall -Wextra -Wpedantic -O3 --std=c++11

main.x : main.o random.o
	$(MPI) random.o main.o -o main.x
main.o : main.cpp
	$(MPI) -c main.cpp -o main.o $(CFLAGS)
random.o : rng/random.cpp rng/random.h
	$(MPI) -c rng/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.x *.out
