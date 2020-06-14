#include <signal.h>
#include <string.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <unistd.h>
// #include "simpi.h"
#include "simpi.cpp"
#include "eigenvalue.cpp"
#include "eigenvector.cpp"

int par_id;

void segfault_printer(int dummy)
{
  char buf[20];
  sprintf(buf, "%d: segfaulted\n", par_id);
  write(STDOUT_FILENO, buf, strlen(buf));
  exit(1);
}
void test_eigenvalue(matrix A,int par_id, int par_count){
   int xy=A.get_x();
   struct timeval stop,start;
   gettimeofday(&start,NULL);
   double *evals = eigenvalue(A,par_id,par_count);
   gettimeofday(&stop,NULL);
   if (par_id==0){
      printf("Runtime: %d microseconds.\n",(stop.tv_usec-start.tv_usec));
      printf("\n---pid:%d--------\n",getpid());
      for (int i=0;i<xy;i++){
         std::cout<<"value"<<i<<":"<<evals[i]<<" ";
      }
      std::cout<<std::endl;
      printf("-----------\n");
   }
}
void test_eigenvector(matrix A,int par_id,int par_count){
   int xy=A.get_x();
   // double V[xy];
   double *V = (double*)mmap(NULL,sizeof(double)*xy,PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
   // std::cout<<"V address1: "<<V<<std::endl;
   // std::cout<<par_id<<"finding vectors"<<std::endl;
   struct timeval stop,start;
   gettimeofday(&start,NULL);
   eigenvector(A,V,par_id,par_count);
   if (par_id==0){
      gettimeofday(&stop,NULL);
      printf("Runtime: %d microseconds.\n",(stop.tv_usec-start.tv_usec));
     for (int i=0;i<xy;i++){
        std::cout<<"V["<<i<<"] ";
        for (int j=0;j<xy;j++){
           std::cout<<V[j+i*xy]<<" ";
        }
        std::cout<<std::endl;
     }
   }
   else{
      printf("non parent\n");
   }
}
int main(int argc, char* argv[])
{
  signal(SIGSEGV, segfault_printer);
  int pid = getpid();
  printf("pid: %d\n",pid);
  par_id = atoi(argv[1]);
  // std::cout<<"par_id "<<par_id<<std::endl;
  unsigned long num_workers = atoi(argv[2]);
  SIMPI_INIT(par_id, num_workers);
  int *pids;
  int pidsfd;
  if (par_id==0){
     pidsfd = shm_open("pids", O_RDWR|O_CREAT, 0777);
     if (pidsfd==-1){
        printf("pids shared memory failed\n");
     }
     ftruncate(pidsfd,sizeof(int)*num_workers);
     pids = (int*)mmap(NULL,sizeof(int)*num_workers,PROT_READ|PROT_WRITE,MAP_SHARED,pidsfd,0);
     pids[num_workers-1] = getpid();
     if (fork()==0){
        sleep(15);
        printf("\n");
        for (int i=0;i<num_workers;i++){
           try{
             kill(pids[i],SIGKILL);
             printf("CLEANUP: %d killed\n",pids[i]);
          }
          catch(int e){
             //pass
          }
       }
       shm_unlink("pids");
       munmap(pids,sizeof(int)*num_workers);
       exit(0);
     }
  }
  main_simpi->synch();
  if (par_id!=0){
     // sleep(2);
     pidsfd = shm_open("pids", O_RDWR, 0777);
     pids = (int*)mmap(NULL,sizeof(int)*num_workers,PROT_READ|PROT_WRITE,MAP_SHARED,pidsfd,0);
     pids[par_id-1] = getpid();
  }
  // std::cout<<"pids-address: "<<pids<<std::endl;
  // SIMPI_SYNCH();
  // main_simpi->synch();



  // vector test:
  int dim = 3;
  matrix B = matrix(dim,dim);
  B.arr[0] = 1;
  B.arr[1] = 2;
  B.arr[2] = 1;
  B.arr[3] = 6;
  B.arr[4] = -1;
  B.arr[5] = 0;
  B.arr[6] = -1;
  B.arr[7] = -2;
  B.arr[8] = -1;





 //  int l = 1;
 //  for (int i=0;i<dim*dim;i++){
 //     // A.arr[i] = 1;
 //     if (i%(dim+1)==0){
 //        B.arr[i] = l;
 //        l++;
 //     }
 //     else{
 //        B.arr[i] = 0;
 //     }
 //  }
 //  int xy = 6;
 //  int n = 5;
 //  int m = 1;
 //  matrix A = matrix(xy,xy);
 //  for (int i=0; i<xy;i++){
 //     for (int j=0;j<xy;j++){
 //        if (j>i){
 //           A.arr[j+xy*i] = n;
 //           n++;
 //        }
 //        else if (j==i){
 //           A.arr[j+xy*i] = m;
 //           m++;
 //        }
 //        else{
 //           A.arr[j+xy*i] = 0;
 //        }
 //     }
 // }

 // A.arr[0] = 0.45;
 // A.arr[5] = 0.1;
 // if (par_id==0){
 // printf("A1\n");
 //    for (int i=0; i<xy*xy;i++){
 //       printf("%.2f ",A.arr[i]);
 //       if ((i+1)%xy==0 ){
 //         printf("\n");
 //      }
 //    }
 // }
  // test_eigenvalue(A,par_id,num_workers);

  test_eigenvector(B,par_id,num_workers);
  // std::cout<<"done\n";
  // SIMPI_SYNCH(); // causing segfault
  // SIMPI_FINALIZE(); // also causing segfault
}
