CC = g++
CFLAGS = -Wall -O3 --std=c++11

LSN_Exercise_101.exe : LSN_Exercise_101.o random.o chromosome.o
	$(CC) random.o chromosome.o LSN_Exercise_101.o -o LSN_Exercise_101.exe
LSN_Exercise_101.o :  LSN_Exercise_101.cpp LSN_Exercise_101.h
	$(CC) -c LSN_Exercise_101.cpp -o LSN_Exercise_101.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
chromosome.o : chromosome.cpp chromosome.h
	$(CC) -c chromosome.cpp -o chromosome.o $(CFLAGS)

clean :
	rm *.o LSN_Exercise_101.exe
