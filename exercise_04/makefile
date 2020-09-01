CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.x : main.o
	$(CC) main.o -o main.x
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.x
	rm data/*.out
	rm -rf frames/*.xyz
