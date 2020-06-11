#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#define SYNCH_OBJECT_MEM_NAME "/simpi_shared_mem"
#define UNIQUE_ID_SIZE 23
typedef struct matrix_metadata {
  char unique_id[UNIQUE_ID_SIZE];
  int file_descriptor;
  size_t size;
  double* matrix_data;
} matrix_metadata;

typedef struct synch_object {
  int par_count;
  char last_matrix_id[UNIQUE_ID_SIZE];
  int ready[];
} synch_object;

// static methods
void SIMPI_INIT(int par_id, size_t synch_size);
void SIMPI_SYNCH();
void SIMPI_FINALIZE();

class simpi {
 public:
  simpi(int _id, int _num_workers);
  ~simpi();
  int get_id() { return id; }
  int get_num_workers() { return num_workers; }
  synch_object* get_synch_info() { return synch_info; }

  std::pair<std::string, double*> create_matrix(int x, int y);

  void free_matrix(std::string unique_id);
  void synch();

 private:
  int id;
  int num_workers;
  int shm_fd;
  synch_object* synch_info;
  std::map<std::string, matrix_metadata> matrix_info;
  std::string sync_shared_mem_name;
  std::string get_shared_mem_name();
};

class vector  // similar stuff for vector
{
 private:
  int dim;
  simpi* mysimpi = NULL;  // for later reference
  std::string unique_id;

 public:
  double* arr;
  int get_size() { return dim; }
  void set(int pos, double val) { arr[pos] = val; }
  double& get(int pos) { return arr[pos]; }
  void print();
  vector(int a);
  ~vector();
};

class matrix  // similar stuff for vector
{
 private:
  int xdim, ydim;
  std::string unique_id;
  simpi* mysimpi = NULL;  // for later reference
 public:
  double* arr;
  matrix(int x, int y);
  ~matrix();
  void print();
  int get_x() { return xdim; }
  int get_y() { return ydim; }
  int determinant(double* A, int n, int order);
  void adjoint(double* A, double* adj, int order, int par_id, int par_count);
  void getCofactor(double* A, double* temp, int p, int q, int n, int order);
  double get_algbera(int pos) { return arr[pos]; }
  void set(int pos, int val) { arr[pos] = val; }
  void inverse_old(matrix *inverse);
  void luDecomposition(matrix* lower, matrix* upper);
  void inverse(matrix* inverse);
  void backward_substitution(float* b, float* x);
  void forward_substitution(float *b, float* x); 
  void solveSystem(vector* constants, vector* solution);
  void jacobi(vector* constants, vector* solution);
  void failSafe(vector* constants, vector* solution);
  bool isDiagonallyDominant();
  friend std::ostream& operator<<(std::ostream& out, const matrix& m);
  double& get(int x, int y) { return arr[x + y * xdim]; }
};
