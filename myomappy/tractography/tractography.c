#include <stdio.h>
#include "Streamline.c"
#include "array_accessing.c"
#include "mask.c"
#include "mat_math.c"
#include "VecField.c"
#include "propagation.c"
// #include "cmemcounter.h"  // for debugging, from: https://github.com/lemire/CMemoryUsage


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
- out_stream_lengths: pointer to 1D int numpy array, where the lengths of the output streamlines are stored
*/
void trackStreams(double *dirs, int *mask, double *affine,
                  int scan_dim1, int scan_dim2, int scan_dim3,
                  int seed_density, double step_size, double max_length, double max_prop_angle, int tract_method,
                  double *out_streams, int *out_stream_lengths) {

    // number of streamlines that will be tracked
    int NUM_SLS = 2 * (seed_density * seed_density * seed_density) * getMaskSum(mask, scan_dim1, scan_dim2, scan_dim3);

    // array of seed points - every seed must correspond to two streamlines. i.e., seed_points is a (NUM_SLS/2) by 3 array
    double *seed_points = generateSeeds(mask, affine, seed_density, scan_dim1, scan_dim2, scan_dim3);

    // creating vector field
    struct VecField vec_field = VecField.newVecField(dirs, mask, affine, scan_dim1, scan_dim2, scan_dim3);

    // tracking bidirectionally from each seed point and saving to streamline array
    int curr_sl_ind = 0;
    for (int i = 0; i < NUM_SLS/2; i++) {  // for each seed point
        // printf("\t\t\t\t\t\t\t%d\n", malloced_memory_usage);  // TESTING
        for (int sign = -1; sign < 2; sign += 2) {  // (-1, 1)
            struct Streamline streamline = Streamline.newSL();  // create streamline object

            double s0[3] = {  // raw seed point, not offset
                getElementDouble2D(seed_points, i, 0, NUM_SLS/2, 3),
                getElementDouble2D(seed_points, i, 1, NUM_SLS/2, 3),
                getElementDouble2D(seed_points, i, 2, NUM_SLS/2, 3)
            };

            // printf("%d / %d\n", i, NUM_SLS/2);  // TESTING

            if (pointInMask(vec_field.tracking_mask, vec_field.tracking_affine, s0, scan_dim1-1, scan_dim2-1, scan_dim3-1) == 0) {  // seed point not in tracking mask M2, cant propogate
                // printf("not in M2 %d / %d\n", i, NUM_SLS/2);  // TESTING
                streamline.append(&streamline, s0);
            } else {  // seed point is in tracking mask, so we can get an initial direction and offset the seed points
                // printf("in M2 %d / %d\n", i, NUM_SLS/2);  // TESTING
                double v0[3];  // vector field at s_0
                vec_field.getValue(&vec_field, s0, v0);
                double *s_offset = vecAdd(s0, vecConstMult(sign*pow(10, -6), v0));  // s_offset = s0 + (sign * 10^-6 * v0)
                streamline.append(&streamline, s_offset);  // append offset seed point
                double *sense_vec = vecConstMult(sign, v0);
                int flag = IN_TRACKING_REGION;  // have ensure that seed point is in M2
                // printf(" %d-%d", i, sign);
                while (flag == IN_TRACKING_REGION) {
                    // printf("%d idk1 %d\n", i, streamline.length);  // TESTING

                    // propagateStream internally appends the next point if appropriate, and returns the appropriate flag
                    // not that sense_vec is only needed for the first tracking step (afterwards, the trajectory of the
                    // streamline is used for sense checking)
                    flag = propagateStream(&vec_field, &streamline, tract_method, step_size, max_length, max_prop_angle, sense_vec);

                    // printf("-p%d", flag);  // TESTING
                    // printf("idk2\n");  // TESTING
                }

                free(s_offset); free(sense_vec);  // free pointers
            }

            // copy streamline into out_streams array (and out_lengths)
            for (int i = 0; i < streamline.length; i++) {
                double curr_p[3];
                streamline.getPoint(&streamline, i, curr_p);
                for (int j = 0; j < 3; j++) {
                    setElementDouble3D(out_streams, curr_p[j], curr_sl_ind, i, j, NUM_SLS, MAX_SL_POINTS, 3);
                }
            }
            out_stream_lengths[curr_sl_ind] = streamline.length;  // store length in lengths array
            streamline.freeSL(&streamline);  // free streamline memory
            curr_sl_ind++;  // update for next streamline
        }
    }


    // // TESTING
    // double *tracking_affine = getTrackingAffine(affine);
    // int *tracking_mask = getTrackingMask(mask, affine, scan_dim1, scan_dim2, scan_dim3);
    // int in_m1 = 0;
    // int in_m2 = 0;
    // int total = 0;
    // for (int i = 0; i < NUM_SLS/2; i++) {
    //     double x0[3] = {
    //         getElementDouble2D(seed_points, i, 0, NUM_SLS/2, 3),
    //         getElementDouble2D(seed_points, i, 1, NUM_SLS/2, 3),
    //         getElementDouble2D(seed_points, i, 2, NUM_SLS/2, 3)
    //     };
    //     in_m1 += pointInMask(mask, affine, x0, scan_dim1, scan_dim2, scan_dim3);
    //     in_m2 += pointInMask(tracking_mask, tracking_affine, x0, scan_dim1-1, scan_dim2-1, scan_dim3-1);
    //     total++;
    // }
    // printf("total: %d\nin m1: %d\nin m2: %d\n", total, in_m1, in_m2);
    // printf("---\n%d\n", getMaskSum(mask, scan_dim1, scan_dim2, scan_dim3));

    // // TESTING
    // for (int i = 0; i < NUM_SLS/2; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         double val = getElementDouble2D(seed_points, i, j, NUM_SLS/2, 3);
    //         setElementDouble2D(test_seeds, val, i, j, NUM_SLS/2, 3);
    //     }
    // }
    // // printMatrix(seed_points, NUM_SLS/2, 3);

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
    
    // // point in mask testing
    // int *test_mask = malloc(scan_dim1 * scan_dim2 * scan_dim3 * sizeof(int));
    // for (int i = 0; i < scan_dim1; i++) {
    //     for (int j = 0; j < scan_dim2; j++) {
    //         for (int k = 0; k < scan_dim3; k++) {
    //             if (!(i == 18 && j == 24 && k == 9)) {
    //                 setElementInt3D(test_mask, 0, i, j, k, scan_dim1, scan_dim2, scan_dim3);
    //             } else {
    //                 setElementInt3D(test_mask, 1, i, j, k, scan_dim1, scan_dim2, scan_dim3);
    //             }
    //         }
    //     }
    // }
    // // double v_1[4][1] = {{18}, {24}, {9}, {1}};
    // // double v_2[4][1] = {{17.5}, {23.5}, {8.5}, {1}};
    // // double v_3[4][1] = {{18.5}, {24.5}, {9.5}, {1}};
    // // double *vp_1 = matMult(affine, v_1, 4, 4, 4, 1);
    // // double *vp_2 = matMult(affine, v_2, 4, 4, 4, 1);
    // // double *vp_3 = matMult(affine, v_3, 4, 4, 4, 1);
    // // printMatrix(vp_1, 4, 1);
    // // printMatrix(vp_2, 4, 1);
    // // printMatrix(vp_3, 4, 1);
    // double p_a[3] = {26.5, 27.6, -158.3};
    // double p_b[3] = {26.6, 27.5, -158.4};
    // double p_c[3] = {31.5, 22.6, -163.3};
    // double p_d[3] = {31.6, 22.5, -163.4};
    // double p_center[3] = {29.009, 25.042, -160.831};
    // int idk_a = pointInMask(test_mask, affine, p_a, scan_dim1, scan_dim2, scan_dim3);
    // int idk_b = pointInMask(test_mask, affine, p_b, scan_dim1, scan_dim2, scan_dim3);
    // int idk_c = pointInMask(test_mask, affine, p_c, scan_dim1, scan_dim2, scan_dim3);
    // int idk_d = pointInMask(test_mask, affine, p_d, scan_dim1, scan_dim2, scan_dim3);
    // int idk_center = pointInMask(test_mask, affine, p_center, scan_dim1, scan_dim2, scan_dim3);
    // printf("%d\n%d\n%d\n%d\n%d\n", idk_a, idk_b, idk_c, idk_d, idk_center);

    // // vec field testing
    // // int *tracking_mask = getTrackingMask(mask, affine, scan_dim1, scan_dim2, scan_dim3);
    // // printf("%d\n", getElementInt3D(tracking_mask, 19, 25, 10, scan_dim1-1, scan_dim2-1, scan_dim3-1));
    // struct VecField vec_field = VecField.newVecField(dirs, mask, affine, scan_dim1, scan_dim2, scan_dim3);
    // // printing dir field values at voxel 19, 25, 10 in tracking mask
    // double *dir_vals = malloc(8 * 3 * sizeof(double));
    // for (int i = 0; i < 2; i++) {
    //     for (int j = 0; j < 2; j++) {
    //         for (int k = 0; k < 2; k++) {
    //             for (int a = 0; a < 3; a++) {
    //                 dir_vals[3*(4*i + 2*j + k) + a] = getElementDouble4D(dirs, 19+i, 25+j, 10+k, a, scan_dim1, scan_dim2, scan_dim3, 3);
    //             }
    //         }
    //     }
    // }
    // printMatrix(dir_vals, 8, 3);
    // // double x0[3] = { 34.00900269 , 20.042499269999993 , -165.83100128 };
    // double x0[3] = { 39.00900269 , 15.042499269999993 , -170.83100128 };
    // double out[3];
    // vec_field.getValue(&vec_field, x0, out);
    // printMatrix(out, 3, 1);
}
