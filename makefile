CC := g++
SRC := src
INC := lib
PROJECT_NAME := main

SRC_FILES := $(shell find $(SRC) -name '*.cpp')
LIB_FILES := $(shell find $(INC) -name '*.h')

inc := -I inc -lgsl -lm -Ilib

CFLAGS := -Wall
OPTIMLEVEL := -O3


all:
	$(CC) $(CFLAGS) $(inc) $(SRC_FILES) $(OPTIMLEVEL) -o  $(PROJECT_NAME) 

.PHONY: clean
clean:
	rm -rf $(PROJECT_NAME)