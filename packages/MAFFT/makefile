all: mafftWrapper


CC = g++

CFLAGS	= -Wall -O3 -m64 -static --debug -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11



MAIN_OBJS = progress.o mafftWrapper.o


mafftWrapper: $(MAIN_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

mafftWrapper.o : mafftWrapper.cpp
	$(CC) $(CFLAGS) -c mafftWrapper.cpp
	
lprogress.o : progress.cpp
	$(CC) $(CFLAGS) -c progress.cpp

clean:
	-rm *.o