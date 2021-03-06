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
MOD_FILES=$(wildcard $(SRC_DIR)/*.ixx)
PCM_FILES=$(patsubst $(SRC_DIR)/%.ixx, $(MOD_DIR)/%.pcm, $(MOD_FILES))

OLC_PGE_HPP=$(INC_DIR)/olcPixelGameEngine.hpp
OLC_PGE_CPP=$(INC_DIR)/olcPixelGameEngine.cpp

CC=clang++
#CC=cl

ifeq ($(CC), cl)
  STD= -std:c++latest /experimental:module
  MFLAGS=
else
  STD= -std=c++2a
  MFLAGS= -Xclang -fimplicit-modules -fimplicit-module-maps -emit-module-interface
endif

#extra flags
EFLAGS= -I$(INC_DIR) -I$(ASSET_DIR) #-O3
#-stdlib=libc++
#pixelgameengine compilation flags
GFLAGS= -luser32 -lgdi32 -lopengl32 -lgdiplus -lShlwapi -ldwmapi #-lstdc++fs #-lc++fs#-lstdc++fs #-llibc++fs

$(MOD_DIR)/%.pcm: $(SRC_DIR)/%.ixx
	$(CC) $(STD) $(EFLAGS) -c $< $(MFLAGS) -emit-module-interface -o $@

$(MAIN): $(SRC_FILES)

$(EXE): $(PCM_FILES) $(MAIN)
	$(CC) $(STD) $(EFLAGS) $(GFLAGS) -fprebuilt-module-path=$(MOD_DIR) $(SRC_FILES) $(OLC_PGE_CPP) -o $@

build: $(EXE)

run: build
	./$(EXE)

.PHONY: clean
clean:
	$(REMOVE) $(PCM_FILES)
