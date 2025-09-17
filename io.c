#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Matrix* load_csv_matrix(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) return NULL;

    /* First pass: count rows and max columns using fgets */
    int rows = 0;
    int cols = 0;
    char linebuf[65536];
    while (fgets(linebuf, sizeof(linebuf), f)) {
        if (linebuf[0] == '\0' || linebuf[0] == '\n') continue;
        rows++;
        int ccount = 1;
        for (char* p = linebuf; *p; ++p) if (*p == ',') ccount++;
        if (ccount > cols) cols = ccount;
    }
    if (rows <= 0 || cols <= 0) { fclose(f); return NULL; }
    rewind(f);

    Matrix* m = create_matrix(rows, cols);
    if (!m) { fclose(f); return NULL; }

    /* Second pass: parse values */
    char buf[65536];
    int r = 0;
    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '\0' || buf[0] == '\n') continue;
        int c = 0;
        char* p = buf;
        while (c < cols) {
            char* end = NULL;
            double v = strtod(p, &end);
            if (p == end) v = 0.0; /* default if parse fails */
            m->data[r][c++] = v;
            if (!end || *end == '\n' || *end == '\0') break;
            if (*end == ',') end++;
            p = end;
        }
        for (; c < cols; ++c) m->data[r][c] = 0.0;
        r++;
        if (r >= rows) break;
    }
    fclose(f);
    return m;
}

int save_csv_matrix(const Matrix* m, const char* filename) {
    if (!m || !filename) return -1;
    FILE* f = fopen(filename, "w");
    if (!f) return -2;
    for (int i = 0; i < m->rows; ++i) {
        for (int j = 0; j < m->cols; ++j) {
            if (j) fputc(',', f);
            fprintf(f, "%.10g", m->data[i][j]);
        }
        fputc('\n', f);
    }
    fclose(f);
    return 0;
}
