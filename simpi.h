#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <arpa/inet.h> 
#include <thread>
#include <exception>      // std::exception_ptr, std::current_esxception, std::rethrow_exception
#include <stdexcept>
#include <array>
#include <sys/ioctl.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sys/socket.h> 
#include <netinet/in.h> 
#include <stdlib.h> 
#include <stdio.h>
#include <netdb.h>
#include <signal.h>
#include <vector>
#include <chrono>
#define PORT 8080 
#define SYNCH_OBJECT_MEM_NAME "/simpi_shared_mem"
#define UNIQUE_ID_SIZE 23


struct data_info
{
    int start;
    int end;
    double *arr;
};

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
void SIMPI_INIT(int par_id, size_t synch_size, int num_workstations, int _workstation_id);
void SIMPI_SYNCH();
void SIMPI_FINALIZE();


class client{
    public:
        int port = 8080;
        int sock; // = 0;
    int setup_client(){
        int valread; 
        struct sockaddr_in serv_addr;  
        if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 
            { 
                printf("\n Socket creation error \n"); 
                return 0; 
            }          
        serv_addr.sin_family = AF_INET; 
        serv_addr.sin_port = htons(port);
        //This is where you put the address of the main server(thats running the client code because in TCP terms we are using it oposite)
        if(inet_pton(AF_INET, "129.65.221.26", &serv_addr.sin_addr)<=0)  
        { 
            printf("\nInvalid address/ Address not supported \n"); 
            return 0; 
        }
        if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
            { 
                printf("\nConnection Failed \n"); 
                return 0; 
            } 
    }
};



class simpi {
 public:
  simpi(int _id, int _num_workers, int _num_workstatons, int _workstation_id);
  ~simpi();
  int get_id() { return id; }
  int get_num_workers() { return num_workers; }
  int get_num_workstations() {return num_workstations;}
  int get_workstation_id() {return workstationid;}
  int get_start() {return start;}
  int get_end() {return end;}
  int get_socket() {return simpi_socket;}
  void set_start(int start_) {start = start_;}
  void set_end(int end_) {end = end_;}
  //client get_client() {return c;}
  synch_object* get_synch_info() { return synch_info; }

  std::pair<std::string, double*> create_matrix(int x, int y);

  void free_matrix(std::string unique_id);
  void synch();

 private:
  int id;
  int num_workers;
  int shm_fd;
  int num_workstations;
  int workstationid;
  int start;
  int end;
  int simpi_socket;
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
  friend int operator*(vector &lhs, vector &rhs);
  simpi *get_simpi() { return mysimpi; }
   vector &scalar_vector_mult(int other);
  vector &add( vector other);
  vector &subtract(vector other);
  int multiply(vector other);
  bool vector_is_equal(vector other);
  friend std::ostream &operator<<(std::ostream &out, const vector &m);
  friend int operator*(vector &lhs, vector &rhs);
  friend vector &operator*(int lhs, vector &rhs);
	friend vector &operator*(vector &lhs, int rhs);
	friend vector &operator+(vector &lhs, vector &rhs);
	friend void operator+=(vector &lhs, vector &rhs);
	friend void operator*=(vector &lhs, int rhs);
	friend void operator-=(vector &lhs, vector &rhs);
	friend vector &operator-(vector &lhs, vector &rhs);
	friend bool operator==(vector &lhs, vector &rhs);
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
  friend std::ostream &operator<<(std::ostream &out, const matrix &m);
	friend matrix &operator*(matrix &lhs, matrix &rhs);
	friend matrix &operator*(int lhs, matrix &rhs);
	friend matrix &operator*(matrix &lhs, int rhs);
	friend void operator*=(matrix &lhs, matrix &rhs);
	friend void operator*=(matrix &lhs, int rhs);
	friend matrix &operator+(matrix &lhs, matrix &rhs);
	friend void operator+=(matrix &lhs, matrix &rhs);
	friend void operator-=(matrix &lhs, matrix &rhs);
	friend matrix &operator-(matrix &lhs, matrix &rhs);
	friend bool operator==(matrix &lhs, matrix &rhs);
  matrix &multiply(matrix other);
  matrix &add(matrix other);
  matrix &scalar_matrix_mult( int other);
  matrix &subtract(matrix other);
  matrix &transpose();
  bool matrix_is_equal(matrix other);
  double *eigenvalue(matrix Anaught,int par_id, int par_count);
  void eigenvector(matrix A, double *V, int par_id, int par_count);
  matrix &SIMPI_DISTRIBUTE();
};

void SIMPI_DISTRIBUTE(matrix m, const matrix &m1);
