#include <stdio.h>
#include <math.h>

typedef struct {
    double  real[2][2]; 
    double  imag[2][2];
}  GL2CMatrix;

extern void copy_GL2C(GL2CMatrix* A, GL2CMatrix* B);
extern void zero_out_GL2C(GL2CMatrix* A);
extern void identity_GL2C(GL2CMatrix* A);
extern void multiply_GL2C(GL2CMatrix* A, GL2CMatrix* B, GL2CMatrix* C); 
extern void trace_GL2C(GL2CMatrix* A, double* real, double* imag);
extern void hash_GL2C(GL2CMatrix* A, long* hash, int bits);
extern void inverse_SL2C(GL2CMatrix* A, GL2CMatrix* B);
extern double norm_GL2(GL2CMatrix* A);
extern int is_one(GL2CMatrix* A, int bits);
