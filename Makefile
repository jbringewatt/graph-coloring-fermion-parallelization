CC=g++ -O2
CFLAGS=-std=c++11 
INC=-I /usr/include/boost
L=-lboost_graph 
all: main clean

main: main.o
	$(CC) $(CFLAGS) -o main main.o $(INC) $(L) 

main.o: main.cpp
	$(CC) -o main.o main.cpp $(INC) $(L) -c 

clean:
	rm -f *.o