## Define Variables
CC = g++
CFLAGS = -Wall -O3 --std=c++11

## Define path to subdirectories
OBJ_DIR = object_files_dir

## Define the object files needed by the excersises
OBJ04_0 = $(OBJ_DIR)/MolDyn_NVE.o  $(OBJ_DIR)/statistical_functions.o

MolDyn_NVE.exe : $(OBJ04_0)
	$(CC) $(OBJ04_0) -o MolDyn_NVE.exe
$(OBJ_DIR)/MolDyn_NVE.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o $(OBJ_DIR)/MolDyn_NVE.o $(CFLAGS)

$(OBJ_DIR)/statistical_functions.o : statistical_functions.cpp statistical_functions.h
	$(CC) -c statistical_functions.cpp -o $(OBJ_DIR)/statistical_functions.o $(CFLAGS)

clean :
	rm $(OBJ_DIR)/*.o *.exe
