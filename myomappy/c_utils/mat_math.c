#pragma once
#include "mat_math.h"
#include "array_accessing.c"
#include <math.h>


/*
compute the dot product of two vectors
inputs:
- v1: double array, the first vector
- v2: double array, the second vector
- length: int, length of v1 and v2 arrays
output:
- returns double, the dot product of the two vectors
*/
double dot(double v1[], double v2[], int length) {
    double ret_sum = 0;
    for (int i = 0; i < length; i++) {
        ret_sum += v1[i] * v2[i];
    }
    return(ret_sum);
}


/*
calculate and return the inverse of a square matrix. assumes the matrix is invertible (det(A) != 0)
inputs:
- *A: pointer to 2D double array, the matrix
- size: int, the size of the matrix (e.g., 2 for a 2x2 matrix)
output:
- returns a pointer to a newly allocated 2D double array, the inverse matrix of A
*/
double *inv(double *A, int size) {
    double *adjunct = transpose(cofactor(A, size), size, size);
    double det = determinant(A, size);
    double *A_inv = constMult(adjunct, 1 / det, size, size);
    return(A_inv);
}


/*
get the cofactor of a square matrix
inputs:
- *A: pointer to 2D double array, the matrix
- size: int, the size of the square matrix (e.g. 2 for a 2x2 matrix)
output:
- returns a pointer to a newly allocated 2D double array, the cofactor matrix of A
*/
double *cofactor(double *A, int size) {
    double *C = malloc(size * size * sizeof(double));  // cofactor matrix to return
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double curr_element = pow(-1.0, i + j + 2);  // uses 'math' indexing for i + j calculation, so offset by 2. shouldnt actually matter
            curr_element *= determinant(minor(A, i + 1, j + 1, size), size - 1);
            // printf("%f\n", curr_element);
            setElementDouble2D(C, curr_element, i, j, size, size);
        }
    }
    return(C);
}


/*
get the minor of a square matrix at a given row and colum
inputs:
- *A: pointer to 2D double array, the matrix
- i: int, the row to calculate the minor at, with 'math' indexing (starting from 1)
- j: int, the column to calculate the minor at, with 'math' indexing (starting from 1)
- size: int, the size of the matrix
output:
- returns a double pointer to a newly allocated 2D array, the resultant minor matrix
*/
double *minor(double *A, int i, int j, int size) {
    double *M = malloc((size - 1) * (size - 1) * sizeof(double));  // minor matrix to return
    int ind = 0;
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            if (row != (i-1) && col != (j-1)) {
                *(M + ind) = *(A + (row * size) + col);
                ind++;
            }
        }
    }
    return(M);
}


/*
calculate the determinant of a square matrix
inputs:
- *A: pointer 2D double array, the matrix
- size: int, the size of the square matrix (e.g., 3)
output:
- returns a double, the determinant of the given matrix at *A
*/
double determinant(double *A, int size) {
    if (size == 2) {  // catch case
        double a11 = getElementDouble2D(A, 0, 0, 2, 2);
        double a12 = getElementDouble2D(A, 0, 1, 2, 2);
        double a21 = getElementDouble2D(A, 1, 0, 2, 2);
        double a22 = getElementDouble2D(A, 1, 1, 2, 2);
        return((a11*a22) - (a12*a21));
    } else {
        double ret_determinant = 0;
        for (int i = 0; i < size; i++) {
            double *curr_minor = minor(A, 1, i+1, size);
            if (i % 2 == 0) {  // even, so add
                ret_determinant += getElementDouble2D(A, 0, i, size, size) * determinant(curr_minor, size - 1);
            } else {
                ret_determinant -= getElementDouble2D(A, 0, i, size, size) * determinant(curr_minor, size - 1);
            }
        }
        return(ret_determinant);
    }
}


/*
get the transpose of a matrix
inputs:
- *A: pointer to 2D double array, the matrix
- num_rows: int, the number of rows in (height of) the matrix
- num_cols: int, the number of columns in (width of) the matrix
output:
- returns a pointer to a newly allocated 2D double array T, the transpose of A. If A is
    shape (m, n), T is shape (n, m)
*/
double *transpose(double *A, int num_rows, int num_cols) {
    double *T = malloc(num_rows * num_cols * sizeof(double));
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            double curr_element = getElementDouble2D(A, i, j, num_rows, num_cols);
            setElementDouble2D(T, curr_element, j, i, num_cols, num_rows);
        }
    }
    return(T);
}


/*
multiply two matrices together via A1 @ A2
inputs:
- *A1: pointer to 2D double array, the first matrix
- *A2: pointer to 2D double array, the second matrix
- rows1: int, the number of rows in (height of) A1
- cols1: int, the number of columns in (width of) A1 (must equal rows2)
- rows2: int, the number of rows in (height of) A2 (must equal cols1)
- cols2: int, the number of columns in (width of) A2
*/
double *matMult(double *A1, double *A2, int rows1, int cols1, int rows2, int cols2) {
    double *ret_mat = malloc(rows1 * cols2 * sizeof(double));
    for (int i = 0; i < rows1; i++) {  // for every row of first matrix
        for (int j = 0; j < cols2; j++) {  // for every column of second matrix
            double curr_dot = 0;
            for (int ind = 0; ind < cols1; ind++) {
                curr_dot += getElementDouble2D(A1, i, ind, rows1, cols1) * getElementDouble2D(A2, ind, j, rows2, cols2);
            }
            setElementDouble2D(ret_mat, curr_dot, i, j, rows1, cols2);
        }
    }
    return(ret_mat);
}


/*
multiply a matrix by a constant element-wise, returning the new matrix
inputs:
- *A: pointer to 2D double array, the matrix
- c: double, the constant to multiply the matrix by
- num_rows: int, the number of rows in (height of) the matrix
- num_cols: int, the number of columns in (width of) the matrix
output:
- returns a pointer to a newly allocated 2D double array, the output matrix of the multiplication
*/
double *constMult(double *A, double c, int num_rows, int num_cols) {
    double *ret_matrix = malloc(num_rows * num_cols * sizeof(double));  // allocate new matrix to return
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            double curr_element = getElementDouble2D(A, i, j, num_rows, num_cols);
            setElementDouble2D(ret_matrix, curr_element * c, i, j, num_rows, num_cols);
        }
    }
    return(ret_matrix);
}


/*
print a given double matrix to the console
inputs:
- *A: pointer to 2D double array (matrix to print)
- num_rows: int, the number of rows in (height of) the matrix
- num_cols: int, the number of columns in (width of) the matrix
output:
- void, prints matrix to the console row-by-row
*/
void printMatrix(double *A, int num_rows, int num_cols) {
    printf("----- matrix -----\n");
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            printf("%f ", *(A + (i * num_cols) + j));
        }
        printf("\n");
    }
    printf("------------------\n");
}
