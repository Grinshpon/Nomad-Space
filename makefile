#SHELL := pwsh -NoProfile

NAME=test

ifeq ($(OS), Windows_NT)
  EXE=$(NAME).exe
  REMOVE=powershell rm
else
  EXE=$(NAME)
  REMOVE=rm
endif

#ASSET_DIR=assets
#TEXTURES=$(wildcard $(ASSET_DIR)/*.png)
SRC_DIR=src
INC_DIR=include
MOD_DIR=mod
ASSET_DIR=assets
MAIN=$(SRC_DIR)/main.cpp
SRC_FILES=$(wildcard $(SRC_DIR)/*.cpp)
MOD_FILES=$(patsubst $(SRC_DIR)/%.cpp, $(MOD_DIR)/%.pcm, $(filter-out $(MAIN), $(SRC_FILES)))

OLC_PGE_HPP=$(INC_DIR)/olcPixelGameEngine.hpp
OLC_PGE_CPP=$(INC_DIR)/olcPixelGameEngine.cpp

CC=clang++
#CC=cl

ifeq ($(CC), cl)
  STD= -std:c++latest /experimental:module
  MFLAGS=
else
  STD= -std=c++2a
  MFLAGS= -Xclang -emit-module-interface
endif

#extra flags
EFLAGS= -I$(INC_DIR) -I$(ASSET_DIR) #-O3
#-stdlib=libc++
#pixelgameengine compilation flags
GFLAGS= -luser32 -lgdi32 -lopengl32 -lgdiplus -lShlwapi -ldwmapi #-lstdc++fs #-lc++fs#-lstdc++fs #-llibc++fs

$(MOD_DIR)/%.pcm: $(SRC_DIR)/%.cpp
	$(CC) $(STD) $(EFLAGS) -c $< $(MFLAGS) -o $@

$(MAIN): $(SRC_FILES)

$(EXE): $(MOD_FILES) $(MAIN)
	$(CC) $(STD) $(EFLAGS) $(GFLAGS) -fprebuilt-module-path=$(MOD_DIR) $(SRC_FILES) $(OLC_PGE_CPP) -o $@

build: $(EXE)

run: build
	./$(EXE)

.PHONY: clean
clean:
	$(REMOVE) $(MOD_FILES)
