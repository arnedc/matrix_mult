long seed;

struct cell {
  int colnumber;
  double value;
  struct cell *next;
};

struct cellrow {
  struct cell *head;
};

struct row {
  int *colnumber;
  double *value;
};

struct sparse_CRS {
  int *colnumber;
  double *value;
  int *rowptr;
};

typedef struct cellrow *sparse2_matrix;
typedef double *dense_matrix;
typedef struct row *sparse_matrix;

void create_matrix(double *,long,int, double, int);
void destroy_dense(dense_matrix);
void multiply_matrix(int n, dense_matrix matA, dense_matrix matB, int lld, dense_matrix matprod);
void multiply_matrix_blocked(int n, dense_matrix matA, dense_matrix matB, int blocksize, int lld, dense_matrix matprod);
void multiply_blas_matrix(int, double *, double *, double *);

sparse_matrix create_sparse_matrix(long, int, double);
void destroy_sparse(int,sparse_matrix);
sparse_matrix multiply_sparse_matrix(int, struct row *, struct row *);
struct sparse_CRS *create_sparse_CRS(long, int, double);
void destroy_sparse_CRS(struct sparse_CRS *);

sparse2_matrix create_sparse2_matrix(long, int, double);
void destroy_sparse2(int,sparse2_matrix);
sparse2_matrix multiply_sparse2_matrix(int , struct cellrow *, struct cellrow *);
sparse2_matrix multiply_sparse2B_matrix(int, sparse2_matrix, sparse2_matrix);

dense_matrix sparse2dense(int,sparse2_matrix);
sparse_matrix sparse2sparse(int,sparse2_matrix);

void printdense(int, dense_matrix, char*);
void printsparse(int, sparse_matrix, char *);
void printsparse2(int, sparse2_matrix, char *);
void printsparse_CRS(int, struct sparse_CRS *, char *);

  


