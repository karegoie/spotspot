#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

typedef struct {
	int rows;
	int cols;
	double** data;
} Matrix;

Matrix* create_matrix(int rows, int cols);
void free_matrix(Matrix* m);
void print_matrix(const Matrix* m);
Matrix* matrix_multiply(const Matrix* A, const Matrix* B);

#endif // MATRIX_H
