#include <signal.h>
#include <string.h>
#include <math.h>

#include "simpi.h"
#define MATRIX_DIMENSION_X 20
#define MATRIX_DIMENSION_Y 20

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
  matrix E(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  matrix F(10,10);
  matrix G(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);

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
      //lower
      B.get(x,y) = x==y;
      //upper
      C.get(x,y) = A.get(x,y);
    }
  }
  // for (int y = 0; y < 10; y++) {
  //   for (int x = 0; x < 10; x++) {
  //     E.get(x, y) = rand()%10 + 1;
  //   }
  // }

  SIMPI_SYNCH();

  
  
  // std::cout << A;
  
  // std::cout << B;
  // C= A * B;
  // SIMPI_SYNCH();
  // //F = D + E;
  // //std::cout << C;
  // SIMPI_DISTRIBUTE(C, G);
  SIMPI_SYNCH();
  // std::cout << G;
  //time_t start,end;
  //start=clock();
  std::cout << A;
  std::cout << "starting LU decomp\n";
  A.newluDecomposition(B, C);

  if (par_id == 0) {
    FILE *file = fopen("file.txt", "w");
    for (int y = 0; y < MATRIX_DIMENSION_Y; y++) {
      for (int x = 0; x < MATRIX_DIMENSION_X; x++) {
        fprintf(file, "%d, ", A.get(x, y));
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }
  // E = B*C;
  // SIMPI_DISTRIBUTE(E,E,1);

  // int num_errors = 0;

  // for (int i = 0; i < MATRIX_DIMENSION_X; i++) {
  //   for (int j = 0; j < MATRIX_DIMENSION_Y; j++) {
  //     if (abs(A.get(i, j) - E.get(i, j)) > 0.001) {
  //       num_errors++;
  //     }
  //   }
  // }
  // printf("errors: %d\n", num_errors);

  // std::cout << A;
  // std::cout << E;
  // std::cout << B;
  // std::cout << C;
  //SIMPI_SYNCH();
  //end=clock();
  // std::cout << C;
  //if (par_id == 1) { 
   // float t= ((float) end-start)/CLOCKS_PER_SEC;
    //std::cout << '\n' << t;
  //}
  //F = D + E;
  //std::cout << C;
  //SIMPI_DISTRIBUTE(C, G);
  //SIMPI_SYNCH();
  
 //SIMPI_SYNCH();
  //SIMPI_DISTRIBUTE(F, out);
  //SIMPI_SYNCH();
  //std::cout << out; 
  //SIMPI_FINALIZE();
  exit(0);
}
