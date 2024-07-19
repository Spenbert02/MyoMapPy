#pragma once

double *generateSeeds(int *mask, double *affine, int seed_density, int dim1, int dim2, int dim3);
double voxelSideLength(double *affine);
int *getVoxelIndices(double *v_to_w_affine, double *x0, double *out_remainders);
double *voxelToWorld(double *affine, double *x1);
int pointInMask(int *mask, double *affine, double *p, int dim1, int dim2, int dim3);
int *getTrackingMask(int *mask, double *affine, int dim1, int dim2, int dim3);
double *getTrackingAffine(double *affine);
int isContinuousVolume(int *mask, int dim1, int dim2, int dim3);
int getMaskSum(int *mask, int dim1, int dim2, int dim3);
