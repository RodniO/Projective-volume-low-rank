#pragma once

struct Matrix {
    void *ptr;
};

struct IntVector {
    void *ptr;
};

enum class AccType : int {
    Half = 0,
    Single = 1,
    Double = 2,
};

extern "C" {

void pvlr_matrix_create(Matrix *matrix, int rows, int cols);

void pvlr_matrix_destroy(Matrix *matrix);

void pvlr_matrix_shape(Matrix const *vector, int *rows, int *cols);

void pvlr_matrix_from_buffer(Matrix *matrix, double const *buffer);

void pvlr_matrix_to_buffer(Matrix const *matrix, double *buffer);

void pvlr_int_vector_create(IntVector *vector, int size);

void pvlr_int_vector_destroy(IntVector *vector);

void pvlr_int_vector_shape(IntVector const *vector, int *size);

void pvlr_int_vector_from_buffer(IntVector *vector, int const *buffer);

void pvlr_int_vector_to_buffer(IntVector const *vector, int *buffer);

void pvlr_aca(Matrix const *A, int rows, int cols, int max_rank, int jpmax,
              Matrix *U, Matrix *V, IntVector const *perm, double rho,
              double rel_err, double rel_err_fro, double abs_err, IntVector *ix,
              IntVector *jx, int max_steps);

// Simplest CUR approximation with $U = \hat A^{-1}$.
//
// Find dominant k-by-k submatrix $\hat A$ and construct $C \hat A^{-1} R$
// approximation.
void pvlr_maxvol(Matrix const *A, int M, int N, int rank, IntVector const *per1,
                 IntVector const *per2, Matrix *C, Matrix *UR, int maxsteps,
                 int maxswaps, Matrix *CA, bool premaxvol);

void pvlr_maxvol2(Matrix const *A, int M, int N, int rank, IntVector *per1,
                  IntVector *per2, Matrix *C, Matrix *UR, int maxsteps,
                  int maxswaps, Matrix *CA, bool premaxvol);

void pvlr_maxvol_proj(Matrix const *A, int M, int N, int rank, int k, int l,
                      IntVector const *per1, IntVector const *per2, Matrix *C,
                      Matrix *UR, int maxsteps, int maxswaps,
                      AccType acc_type = AccType::Double);

} // extern "C"
