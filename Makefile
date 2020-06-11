CC = g++ 
CFLAGS = -Wall -g -std=c++11
LINKFLAGS = -lrt
# if mac use this LINFLAGS
# LINKFLAGS =

all : mpi user
clean : 
	rm -f user mpi /dev/shm/simpi_shared_mem
simpi : simpi.cpp simpi.h
	$(CC) $(CFLAGS) -c simpi.cpp -o simpi $(LINKFLAGS)
mpi : mpi.cpp user  simpi simpi.h
	$(CC) $(CFLAGS) mpi.cpp -o mpi simpi $(LINKFLAGS)
user : user.cpp simpi simpi.h
	$(CC) $(CFLAGS) user.cpp -o user simpi $(LINKFLAGS)
