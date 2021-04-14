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
  int workstationid = atoi(argv[4]);
  SIMPI_INIT(par_id, num_workers, num_workstations, workstationid);
  //user code goes here 
  matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  matrix B(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  matrix C(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  matrix D(10, 10);
  matrix E(10, 10);
  matrix F(10,10);
  matrix G(10, 10);
  matrix out(10, 10);
  //vector D(10);
  SIMPI_SYNCH();

  for (int y = 0; y < MATRIX_DIMENSION_Y; y++) {
    for (int x = 0; x < MATRIX_DIMENSION_X; x++) {
      A.get(x, y) = rand()%10 + 1;
    }
  }
  for (int y = 0; y < MATRIX_DIMENSION_Y; y++) {
    for (int x = 0; x < MATRIX_DIMENSION_X; x++) {
      B.get(x, y) = rand()%10 + 1;
    }
  }
  for (int y = 0; y < 10; y++) {
    for (int x = 0; x < 10; x++) {
      D.get(x, y) = rand()%10 + 1;
    }
  }
  for (int y = 0; y < 10; y++) {
    for (int x = 0; x < 10; x++) {
      E.get(x, y) = rand()%10 + 1;
    }
  }

  SIMPI_SYNCH();

  
  std::cout << A;
  
  std::cout << B;
  C= A * B;
  SIMPI_SYNCH();
  F = D * E;
  //std::cout << C;
  SIMPI_DISTRIBUTE(C, G);
  SIMPI_SYNCH();
  std::cout << G;
  SIMPI_SYNCH();
  SIMPI_DISTRIBUTE(F, out);
  SIMPI_SYNCH();
  std::cout << out; 
  SIMPI_FINALIZE();
  exit(0);
}
