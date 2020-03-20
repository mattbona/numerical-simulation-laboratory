CC = g++
CFLAGS = -Wall -O3 --std=c++11

LSN_Exercise_011.exe : LSN_Exercise_011.o random.o
	$(CC) random.o LSN_Exercise_011.o -o LSN_Exercise_011.exe
LSN_Exercise_011.o : LSN_Exercise_011.cpp
	$(CC) -c LSN_Exercise_011.cpp -o LSN_Exercise_011.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o LSN_Exercise_011.exe seed.out
