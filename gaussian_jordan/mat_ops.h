#ifndef MAT_OPS
#define MAT_OPS
#include <inttypes.h>
#include <stdbool.h> 
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "signal.h"

int easy_rand(int min, int max);

void s_xorshift96(uint32_t usr_def_seed);

uint32_t xorshift96(int min, int max);

void free_ptr(int r, double **p);

double** mat_zeros(int row, int col);

double** mat_ones(int row, int col);

double** mat_rand(int row, int col, double min, double max);

double** mat_mul(int m, int p, int n, double** mat_l, double** mat_r);

double** mat_trans(int row, int col, double** mat_old);

double** mat_diag(int col, double** arr);

double** mat_add(int m, int n, double** mat1, double** mat2);

void print_mat(int r, int c, double** mat);

void copy_mat(int r, int c, double** new, double** old);

#endif
