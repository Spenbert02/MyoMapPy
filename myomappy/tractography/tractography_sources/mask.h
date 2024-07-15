#pragma once

int pointInMask(int *mask, double *affine, double *p, int dim1, int dim2, int dim3);
int pointInVoxel(double *affine, int x_vox, int y_vox, int z_vox, double *p);
int *getTrackingMask(int *mask, double *affine, int dim1, int dim2, int dim3);
double *getTrackingAffine(double *affine);
int isContinuousVolume(int *mask, int dim1, int dim2, int dim3);
int getMaskSum(int *mask, int dim1, int dim2, int dim3);
