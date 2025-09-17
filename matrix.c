#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

Matrix* create_matrix(int rows, int cols) {
	if (rows <= 0 || cols <= 0) {
		return NULL;
	}

	Matrix* m = (Matrix*)malloc(sizeof(Matrix));
	if (!m) {
		return NULL;
	}
	m->rows = rows;
	m->cols = cols;
	m->data = (double**)malloc((size_t)rows * sizeof(double*));
	if (!m->data) {
		free(m);
		return NULL;
	}

	for (int i = 0; i < rows; ++i) {
		m->data[i] = (double*)malloc((size_t)cols * sizeof(double));
		if (!m->data[i]) {
			// free previously allocated rows
			for (int k = 0; k < i; ++k) {
				free(m->data[k]);
			}
			free(m->data);
			free(m);
			return NULL;
		}
		// initialize to 0.0
		for (int j = 0; j < cols; ++j) {
			m->data[i][j] = 0.0;
		}
	}

	return m;
}

void free_matrix(Matrix* m) {
	if (!m) return;
	if (m->data) {
		for (int i = 0; i < m->rows; ++i) {
			free(m->data[i]);
		}
		free(m->data);
	}
	free(m);
}

void print_matrix(const Matrix* m) {
	if (!m || !m->data) {
		printf("(null)\n");
		return;
	}
	for (int i = 0; i < m->rows; ++i) {
		for (int j = 0; j < m->cols; ++j) {
			printf("%g%s", m->data[i][j], (j + 1 == m->cols) ? "" : " ");
		}
		printf("\n");
	}
}

Matrix* matrix_multiply(const Matrix* A, const Matrix* B) {
	if (!A || !B || !A->data || !B->data) {
		return NULL;
	}
	if (A->cols != B->rows) {
		return NULL;
	}

	Matrix* C = create_matrix(A->rows, B->cols);
	if (!C) {
		return NULL;
	}

	for (int i = 0; i < A->rows; ++i) {
		for (int k = 0; k < A->cols; ++k) {
			double aik = A->data[i][k];
			if (aik == 0.0) continue; // small optimization
			for (int j = 0; j < B->cols; ++j) {
				C->data[i][j] += aik * B->data[k][j];
			}
		}
	}

	return C;
}

