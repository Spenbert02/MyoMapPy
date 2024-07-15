#pragma once


/*
get/set element from 2D double array
- **arr: pointer to 2D double array
- val (setElementDouble2D only): double, value to set at arr[i, j]
- i: int, first dimension index to access
- j: int, second dimension index to access
- dim1: int, size of first dimension of array
- dim2: int, size of second dimension of array
*/
double getElementDouble2D(double *arr, int i, int j, int dim1, int dim2) {
    unsigned long long offset = (i * dim2) + j;
    double retVal = *(arr + offset);
    return(retVal);
}
void setElementDouble2D(double *arr, double val, int i, int j, int dim1, int dim2) {
    unsigned long long offset = (i * dim2) + j;
    *(arr + offset) = val;
}


/*
get/set element from 3D int array
inputs:
- ***arr: pointer to 3D int array
- val (setElementInt3D only): int, value to set at arr[i, j, k]
- i: int, first dimension index to access
- j: int, second dimension index to access
- k: int, third dimension index to access
- dim1: int, size of first dimension of array
- dim2: int, size of second dimension of array
- dim3: int, size of third dimension of array
*/
int getElementInt3D(int *arr, int i, int j, int k, int dim1, int dim2, int dim3) {
    int offset = (i * dim2 * dim3) + (j * dim3) + k;
    int retVal = *(arr + offset);
    return(retVal);
}
void setElementInt3D(int *arr, int val, int i, int j, int k, int dim1, int dim2, int dim3) {
    unsigned long long offset = (i * dim2 * dim3) + (j * dim3) + k;
    *(arr + offset) = val;
}


/*
get/set element from 3D double array
inputs:
- ***arr: pointer to 3D double array
- val (setElementDouble3D only): double, value to set at arr[i, j, k]
- i: int, first dimension index to access
- j: int, second dimension index to access
- k: int, third dimension index to access
- dim1: int, size of first dimension of arr
- dim2: int, size of second dimension of arr
- dim3: int, size of third dimension of arr
output:
- returns a double, the value at arr[i, j, k]
*/
double getElementDouble3D(int *arr, int i, int j, int k, int dim1, int dim2, int dim3) {
    unsigned long long offset = (i * dim2 * dim3) + (j * dim3) + k;
    double retVal = *(arr + offset);
    return(retVal);
}
void setElementDouble3D(int *arr, double val, int i, int j, int k, int dim1, int dim2, int dim3) {
    unsigned long long offset = (i * dim2 * dim3) + (j * dim3) + k;
    *(arr + offset) = val;
}


/*
get/set element from 4D double array
inputs:
- ****arr: pointer to 4D double array
- val (setElementDouble4D only): double, the value to set at arr[i, j, k, l]
- i: int, first dimension index to access
- j: int, second dimension index to access
- k: int, third dimension index to access
- l: int, fourth dimension index to access
- dim1: int, size of first dimension of arr
- dim2: int, size of second dimension of arr
- dim3: int, size of third dimension of arr
- dim4: int, size of fourth dimension of arr
output:
- returns a double, the value at arr[i, j, k, l]
*/
double getElementDouble4D(double *arr, int i, int j, int k, int l, int dim1, int dim2, int dim3, int dim4) {
    unsigned long long offset = (i * dim2 * dim3 * dim4) + (j * dim3 * dim4) + (k * dim4) + l;
    double retVal = *(arr + offset);
    return(retVal);
}
void setElementDouble4D(double *arr, double val, int i, int j, int k, int l, int dim1, int dim2, int dim3, int dim4) {
    unsigned long long offset = (i * dim2 * dim3 * dim4) + (j * dim3 * dim4) + (k * dim4) + l;
    *(arr + offset) = val;
}
