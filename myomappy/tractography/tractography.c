#include <stdio.h>
#include "Streamline.c"
#include "array_accessing.c"
#include "mask.c"
#include "mat_math.c"


// hard-coded maximum number of points per streamline
int MAX_SL_POINTS = 1000;

// tractography method enumerations
const int EULER = 0;
const int RUNGE_KUTTA_4 = 1;


/*
generate a set of streamlines from input data, mask, affine, and tractography parameters
inputs:
- dirs_v: pointer to 4D float numpy array, the normalized direction vector in each voxel
- mask_v: pointer to 3D int numpy array, the binary mask for the scan
- affine_v: pointer to shape (4, 4) float numpy array, the voxel-to-world affine matrix for the scan
- scan_dim1: int, the first dimension (# x-voxels) of the scan, dirs.shape[0]=mask.shape[0]
- scan_dim2: int, the second dimension (# y-voxels) of the scan, dirs.shape[1]=mask.shape[1]
- scan_dim3: int, the third dimension (# z-voxels) of the scan, dirs.shape[2]=mask.shape[2]
- seed_density: int, the seed density per voxel (e.g., seed_density=2 yields 2x2x2=8 seeds per voxel)
- step_size: double, the step size for tractography in the same units as affine
- max_length: double, the maximum streamline length for tractography in the same units as affine
    if -1.0, no max length is used during tracking
- max_prop_angle: double, the maximum allowable propogation angle during tractography in degrees
    if -1.0, no max prop angle is used during tracking
- tract_method: int (enumeration), the tractography method to use
    0: euler
    1: 4th order runge-kutta
- out_streams: pointer to 3D float numpy array, where the output streamlines are stored. Should be
    of shape (n, myomappy.tractography.MAX_SL_POINTS, 3), where n is twice the total number of seed
    points
*/
void trackStreams(double *dirs, int *mask, double *affine,
                  size_t scan_dim1, size_t scan_dim2, size_t scan_dim3,
                  int seed_density, double step_size, double max_length, double max_prop_angle, int tract_method,
                  double *out_streams, int *out_stream_lengths) {  // testing

    // number of streamlines that will be tracked
    int NUM_SLS = 2 * (seed_density * seed_density * seed_density) * getMaskSum(mask, scan_dim1, scan_dim2, scan_dim3);


    // TESTING
    // printf("%d\n", NUM_SLS);
    // printf("%f\n", getElementDouble4D(dirs, 10, 10, 10, 2, scan_dim1, scan_dim2, scan_dim3, 3));
    // printf("%d\n", dirs[(10 * scan_dim2 * scan_dim3 * 3) + (10 * scan_dim3 * 3) + (10 * 3) + 2]);
    // printf("%d\n", getElementInt3D(mask, 10, 11, 13, scan_dim1, scan_dim2, scan_dim3));
    // printf("%d\n", getMaskSum(mask, scan_dim1, scan_dim2, scan_dim3));
    // double v1[4] = {1, 2, 3, 4};
    // double v2[4] = {5, 6, 7, 8};
    // printf("%f\n", dot(v1, v2, 4));
    // double A[3][3] = {{1, 4, 3}, {5, 5, 6}, {10, 8, 9}};
    // double *M = minor(A, 2, 3, 3);
    // printMatrix(A, 3, 3);
    // printf("%f\n", determinant(A, 3));
    // printMatrix(M, 2, 2);
    // printf("%f\n", determinant(M, 2));
    // double B[4][4] = {{0.8, 0, 0, 25}, {0, 0.8, 0, 50}, {0, 0, 0.8, 75}, {0, 0, 0, 1}};
    // printMatrix(B, 4, 4);
    // printf("%f\n", determinant(B, 4));
    // double *C = cofactor(B, 4);
    // double *multed = constMult(C, 2, 4, 4);
    // printMatrix(C, 4, 4);
    // printMatrix(multed, 4, 4);
    // double *B_inv = inv(B, 4);
    // printMatrix(B, 4, 4);
    // printMatrix(B_inv, 4, 4);
    // for (int i = -1; i < 2; i += 2) {
    //     printf("%d\n", i);
    // }
    // setElementInt3D(mask, 1, 0, 0, 0, scan_dim1, scan_dim2, scan_dim3);
    // printf("%d\n", isContinuousVolume(mask, scan_dim1, scan_dim2, scan_dim3));
    // double A[2][3] = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    // double B[3][2] = {{7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}};
    // printMatrix(matMult(A, B, 2, 3, 3, 2), 2, 2);
    // int *temp_tracking_mask = getTrackingMask(mask, affine, scan_dim1, scan_dim2, scan_dim3);
    // double *temp_tracking_affine = getTrackingAffine(affine);
    // printf("%d\n", isContinuousVolume(temp_tracking_mask, scan_dim1-1, scan_dim2-1, scan_dim3-1));
    // for (int i = 0; i < scan_dim1 - 1; i++) {
    //     for (int j = 0; j < scan_dim2 - 1; j++) {
    //         for (int k = 0; k < scan_dim3 - 1; k++) {
    //             setElementInt3D(tracking_mask, getElementInt3D(temp_tracking_mask, i, j, k, scan_dim1-1, scan_dim2-1, scan_dim3-1), i, j, k, scan_dim1-1, scan_dim2-1, scan_dim3-1);
    //         }
    //     }
    // }
    // for (int i = 0; i < 4; i++) {
    //     for (int j = 0; j < 4; j++) {
    //         setElementDouble2D(tracking_affine, getElementDouble2D(temp_tracking_affine, i, j, 4, 4), i, j, 4, 4);
    //     }
    // }
    double p[3] = {0, 0, 0};
    int idk = pointInMask(mask, affine, p, scan_dim1, scan_dim2, scan_dim3);
}
