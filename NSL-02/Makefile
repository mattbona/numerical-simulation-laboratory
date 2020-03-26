## Define Variables
CC = g++
CFLAGS = -Wall -O3 --std=c++11

## Define path to subdirectories
OBJ_DIR = OBJ_DIR

## Define the object files needed by the excersises
OBJ02_1 = $(OBJ_DIR)/LSN_Exercise_02_1.o $(OBJ_DIR)/random.o $(OBJ_DIR)/statistical_functions.o

EX02_1.exe : $(OBJ02_1)
	$(CC) $(OBJ02_1) -o EX02_1.exe
$(OBJ_DIR)/LSN_Exercise_02_1.o : LSN_Exercise_02_1.cpp
	$(CC) -c LSN_Exercise_02_1.cpp -o $(OBJ_DIR)/LSN_Exercise_02_1.o $(CFLAGS)

$(OBJ_DIR)/random.o : random.cpp random.h
	$(CC) -c random.cpp -o $(OBJ_DIR)/random.o $(CFLAGS)
$(OBJ_DIR)/statistical_functions.o : statistical_functions.cpp statistical_functions.h
	$(CC) -c statistical_functions.cpp -o $(OBJ_DIR)/statistical_functions.o $(CFLAGS)


clean :
	rm $(OBJ_DIR)/*.o *.exe
