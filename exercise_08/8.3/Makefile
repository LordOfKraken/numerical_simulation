CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : main.o random.o
	$(CC) random.o main.o -o main.x
main.o : qmc1d.cpp
	$(CC) -c qmc1d.cpp -o main.o $(CFLAGS)
random.o : rng/random.cpp rng/random.h
	$(CC) -c rng/random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o main.x potential.dat kinetic.dat probability.dat
