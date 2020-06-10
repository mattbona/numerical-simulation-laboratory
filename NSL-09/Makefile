CC = g++
CFLAGS = -Wall -O3 --std=c++11

LSN_Exercise_091.exe : LSN_Exercise_091.o random.o
	$(CC) random.o LSN_Exercise_091.o -o LSN_Exercise_091.exe
LSN_Exercise_091.o : LSN_Exercise_091.cpp
	$(CC) -c LSN_Exercise_091.cpp -o LSN_Exercise_091.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o LSN_Exercise_091.exe
