#pragma once

double *vecAdd(double *v1, double *v2);
double *vecSubtract(double *v1, double *v2);
double *vecConstMult(double c, double *v);
double dot(double v1[], double v2[], int length);
double *norm(double v[]);
double *inv(double *A, int size);
double *cofactor(double *A, int size);
double *minor(double *A, int i, int j, int size);
double determinant(double *A, int size);
double *transpose(double *A, int num_rows, int num_columns);
double *matMult(double *A1, double *A2, int rows1, int cols1, int rows2, int cols2);
double *constMult(double *A, double c, int num_rows, int num_cols);
void printMatrix(double *A, int num_rows, int num_cols);
