#include <signal.h>
#include <string.h>

#include "simpi.h"
#define MATRIX_DIMENSION_X 10
#define MATRIX_DIMENSION_Y 10

int par_id;

void segfault_printer(int dummy)
{
  char buf[20];
  sprintf(buf, "%d: segfaulted\n", par_id);
  write(STDOUT_FILENO, buf, strlen(buf));
  exit(1);
}

int main(int argc, char* argv[])
{
  printf("user starting...\n");
  signal(SIGSEGV, segfault_printer);
  par_id = atoi(argv[1]);
  int num_workers = atoi(argv[2]);
  int num_workstations = atoi(argv[3]);
  SIMPI_INIT(par_id, num_workers, num_workstations);
  matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  matrix C(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  vector D(10);
  SIMPI_SYNCH();

  for (int y = 0; y < MATRIX_DIMENSION_Y; y++) {
    for (int x = 0; x < MATRIX_DIMENSION_X; x++) {
      A.get(x, y) = rand()%10 + 1;
    }
  }
  SIMPI_SYNCH();

  A.inverse(&C);
  std::cout << A;
  
  std::cout << C;

  //C= A.inverse();

  //std::cout << C; 
  //SIMPI_FINALIZE();
}
