CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: BMSAMS



BMSAMS: BMSAMS.c
	$(CC) $(CFLAGS) alignmentFunctions.c -lm commonFunctions.c -lm BMSAMS.c -o $(BIN)/BMSAMS
	

clean:
	rm -rf $(BIN)/BMSAMS
