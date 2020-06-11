#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include <cstring>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "simpi.h"

int main(int argc, char* argv[])
{
  if (argc != 3) {
    printf("Usage: ./prog2 prog1 <num workers>");
    exit(2);
  }

  char progname[100];
  int fd;
  // strcpy(progname, "./");
  strcpy(progname, argv[1]);

  int numWorkers = atoi(argv[2]);

  // create shared mem for workers
  size_t synchObjectSize =
      sizeof(synch_object) + sizeof(int) * (numWorkers + 1);
  fd = shm_open(SYNCH_OBJECT_MEM_NAME, O_RDWR | O_CREAT, 0777);
  if (fd == -1) {
    perror("Unable to create synch info: ");
    exit(1);
  }
  ftruncate(fd, synchObjectSize);
  synch_object* shared_mem = (synch_object*)mmap(
      NULL, synchObjectSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

  if (shared_mem == MAP_FAILED) {
    perror("Unable to mmap shared memory: ");
    exit(1);
  }
  // initialize ready to zero
  for (int i = 0; i <= numWorkers; i++) {
    shared_mem->ready[i] = 0;
  }
  shared_mem->par_count = numWorkers;
  for (int i = 0; i < numWorkers; i++) {
    std::string worker_count_str = std::to_string(numWorkers);
    std::string worker_id_str = std::to_string(i);
    char* args[] = {progname, const_cast<char*>(worker_id_str.c_str()),
                    const_cast<char*>(worker_count_str.c_str()), NULL};
    if (fork() == 0) {
      execv(progname, args);
    }
  }
  exit(0);
}
