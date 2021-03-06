CC = g++ 
CFLAGS = -Wall -g -std=c++11 -pthread
LINKFLAGS = -lrt
# if mac use this LINFLAGS
# LINKFLAGS =

all : mpi user client server
clean : 
	rm /dev/shm/*
	pkill -f user
	rm server
	rm client
simpi : simpi.cpp simpi.h
	$(CC) $(CFLAGS) -c simpi.cpp -o simpi $(LINKFLAGS)
mpi : mpi.cpp user  simpi simpi.h
	$(CC) $(CFLAGS) mpi.cpp -o mpi simpi $(LINKFLAGS)
user : user.cpp simpi simpi.h
	$(CC) $(CFLAGS) user.cpp -o user simpi $(LINKFLAGS)
client : simpclient.cpp 
	$(CC) $(CFLAGS) simpclient.cpp -o client $(LINKFLAGS)
server : simpserver.cpp 
	$(CC) $(CFLAGS) simpserver.cpp -o server $(LINKFLAGS)
ipreader : ipreader.cpp 
	$(CC) $(CFLAGS) ipreader.cpp -o ipreader $(LINKFLAGS)
