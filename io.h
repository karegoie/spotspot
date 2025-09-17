#ifndef IO_H
#define IO_H

#include "matrix.h"

/* Load CSV (comma-separated) into a Matrix (rows x cols). Returns NULL on error. */
Matrix* load_csv_matrix(const char* filename);

/* Save Matrix to CSV. Returns 0 on success, non-zero on error. */
int save_csv_matrix(const Matrix* m, const char* filename);

#endif /* IO_H */