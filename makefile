all: FMAlign


CC = g++


CFLAGS	= -Wall -O3 -m64 -static --debug -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11
CLINK	= -lm -O3 --debug -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11


BIN_DIR = bin
ROOT_DIR = $(shell pwd)
FMTREE_LIB_DIR = FMtree


MAIN_OBJS = load_seqs.o chain.o common_seed.o subseqs.o command_line.o my_main.o


FMTREE_OBJS = \
$(FMTREE_LIB_DIR)/bwt.o\
$(FMTREE_LIB_DIR)/saca-k.o



HEADER = Auxiliary.h


FMAlign: $(FMTREE_OBJS) $(MAIN_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

my_main.o : my_main.cpp $(HEADER)
	$(CC) $(CFLAGS) -c my_main.cpp

command_line.o : command_line.cpp $(HEADER)
	$(CC) $(CFLAGS) -c command_line.cpp

subseqs.o : subseqs.cpp $(HEADER)
	$(CC) $(CFLAGS) -c subseqs.cpp

load_seqs.o : load_seqs.cpp $(HEADER)
	$(CC) $(CFLAGS) -c load_seqs.cpp

chain.o : chain.cpp 
	$(CC) $(CFLAGS) -c chain.cpp

common_seed.o : common_seed.cpp 
	$(CC) $(CFLAGS) -c common_seed.cpp
	
bwt.o : $(FMTREE_LIB_DIR)/bwt.cpp
	$(CC) $(CFLAGS) -c $(FMTREE_LIB_DIR)/bwt.cpp

saca-k.o : $(FMTREE_LIB_DIR)/saca-k.cpp $(HEADER)
	$(CC) $(CFLAGS) -c $(FMTREE_LIB_DIR)/saca-k.cpp
	
clean:
	-rm -f $(ROOT_DIR)/*.o
	-rm -f $(FMTREE_LIB_DIR)/*.o
