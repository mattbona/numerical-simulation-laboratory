## Define Variables
CC = g++
CFLAGS = -Wall -O3 --std=c++11

## Define path to subdirectories
OBJ_DIR = OBJ_DIR

## Define the object files needed by the excersises
OBJ011 = $(OBJ_DIR)/LSN_Exercise_011.o $(OBJ_DIR)/random.o $(OBJ_DIR)/statistical_functions.o
OBJ012 = $(OBJ_DIR)/LSN_Exercise_012.o $(OBJ_DIR)/random.o $(OBJ_DIR)/statistical_functions.o $(OBJ_DIR)/central_limit.o


EX011.exe : $(OBJ011)
	$(CC) $(OBJ011) -o EX011.exe
$(OBJ_DIR)/LSN_Exercise_011.o : LSN_Exercise_011.cpp
	$(CC) -c LSN_Exercise_011.cpp -o $(OBJ_DIR)/LSN_Exercise_011.o $(CFLAGS)

EX012.exe : $(OBJ012)
	$(CC) $(OBJ012) -o EX012.exe
$(OBJ_DIR)/LSN_Exercise_012.o : LSN_Exercise_012.cpp
	$(CC) -c LSN_Exercise_012.cpp -o $(OBJ_DIR)/LSN_Exercise_012.o $(CFLAGS)
$(OBJ_DIR)/central_limit.o : central_limit.cpp central_limit.h
	$(CC) -c central_limit.cpp -o $(OBJ_DIR)/central_limit.o $(CFLAGS)

$(OBJ_DIR)/random.o : random.cpp random.h
	$(CC) -c random.cpp -o $(OBJ_DIR)/random.o $(CFLAGS)
$(OBJ_DIR)/statistical_functions.o : statistical_functions.cpp statistical_functions.h
	$(CC) -c statistical_functions.cpp -o $(OBJ_DIR)/statistical_functions.o $(CFLAGS)


clean :
	rm $(OBJ_DIR)/*.o *.exe
