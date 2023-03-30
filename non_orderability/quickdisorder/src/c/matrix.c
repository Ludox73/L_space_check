#include "matrix.h"

void copy_GL2C(GL2CMatrix* A, GL2CMatrix* B){
    int i, j;
    for (i=0; i < 2; i++){
	for (j=0; j < 2; j++){
	    B->real[i][j] = A->real[i][j];
	    B->imag[i][j] = A->imag[i][j];
	}
    }
}

void zero_out_GL2C(GL2CMatrix* A){
    int i, j;
    for (i=0; i < 2; i++){
	for (j=0; j < 2; j++){
	    A->real[i][j] = 0.0;
	    A->imag[i][j] = 0.0;
	}
    }
}

void identity_GL2C(GL2CMatrix* A){
    int i, j;
    for (i=0; i < 2; i++){
	for (j=0; j < 2; j++){
	    A->real[i][j] = (i == j) ? 1.0 : 0.0;
	    A->imag[i][j] = 0.0;
	}
    }
}

void multiply_GL2C(GL2CMatrix* A, GL2CMatrix* B, GL2CMatrix* C){
    int i, j, k;
    zero_out_GL2C(C);
    for (i=0; i < 2; i++){
	for (j=0; j < 2; j++){
	    for (k=0; k < 2; k++){
		C->real[i][j] += (A->real[i][k])*(B->real[k][j]) - (A->imag[i][k])*(B->imag[k][j]);
		C->imag[i][j] +=  (A->real[i][k])*(B->imag[k][j]) + (A->imag[i][k])*(B->real[k][j]);
	    }
	}
    }
}

void trace_GL2C(GL2CMatrix* A, double* real, double* imag){
    int i;
    double a, b; 
    a = 0;
    b = 0;
    for (i=0; i < 2; i++){
	a += A->real[i][i];
	b += A->imag[i][i];
    }
    *real = a;
    *imag = b;
}

void hash_GL2C(GL2CMatrix* A, long* hash, int bits){
    int i, j, N; 
    N = 1 << bits; 
    for (i=0; i < 2; i++){
	for (j=0; j < 2; j++){
	    hash[4*i + 2*j] = N*A->real[i][j];
	    hash[4*i + 2*j + 1] = N*A->imag[i][j];
	}
    }
}
			
void inverse_SL2C(GL2CMatrix* A, GL2CMatrix* B){
    B->real[0][0] = A-> real[1][1];
    B->real[1][1] = A-> real[0][0];
    B->real[0][1] = -A->real[0][1];
    B->real[1][0] = -A->real[1][0];
    B->imag[0][0] = A->imag[1][1];
    B->imag[1][1] = A->imag[0][0];
    B->imag[0][1] = -A->imag[0][1];
    B->imag[1][0] = -A->imag[1][0];
}

double norm_GL2(GL2CMatrix* A){
    int i, j;
    double ans, entry_size;
    ans = 0.0;
    for (i=0; i < 2; i++){
	for (j=0; j < 2; j++){
	    entry_size = sqrt( A->real[i][j]*A->real[i][j] +
			       A->imag[i][j]*A->imag[i][j] );
	    ans = fmax(ans, entry_size);
	}
    }
    return log2(ans);
}

int is_one(GL2CMatrix* A, int bits){
    GL2CMatrix B;
    long hash[8];
    int i;
    copy_GL2C(A, &B);
    B.real[0][0] += -1;
    B.real[1][1] += -1;
    hash_GL2C(&B, hash, bits);
    for (i=0; i < 8; i++){
	if(hash[i] != 0){
	    return 0;
	}
    }
    return 1;
}
		
	
