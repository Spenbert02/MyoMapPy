#pragma once
#include <math.h>
#include "mask.h"
#include "array_accessing.c"
#include "mat_math.c"


/*
generate uniformly distributed seeds throughout a mask
inputs:
- *mask: pointer to 3D int array, the binary mask
- *affine: pointer to 4x4 double array, the voxel-to-world affine matrix
- seed_density: int, the number of seeds per voxel per axis (i.e., seed_density=2 yields 2x2x2=8 seeds per voxel)
- dim1: int, size of the first dimension of the scan (mask)
- dim2: int, size of the second dimension of the scan (mask)
- dim3: int, size of the third dimension of the scan (mask)
output:
- returns pointer to a newly allocated 2D double array, the list of seed points. the array has
    shape (seed_density^3 * <# voxels in mask>, 3)
*/
double *generateSeeds(int *mask, double *affine, int seed_density, int dim1, int dim2, int dim3) {
    int num_seeds = seed_density * seed_density * seed_density * getMaskSum(mask, dim1, dim2, dim3);  // number of seed points (and streamlines)
    double *ret_seeds = malloc(3 * num_seeds * sizeof(double));  // return matrix for seed points (num_seeds x 3)

    int curr_seed = 0;
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            for (int k = 0; k < dim3; k++) {  // for every voxel
                if (getElementInt3D(mask, i, j, k, dim1, dim2, dim3) == 1) {  // if voxel in mask, create seeds
                    double x1_min[3] = {i - 0.5 + (1/(2*seed_density)), j - 0.5 + (1/(2*seed_density)), k - 0.5 + (1/(2*seed_density))};  // 'minimum' seed point (nearest to origin in voxel space)
                    for (int si = 0; si < seed_density; si++) {
                        for (int sj = 0; sj < seed_density; sj++) {
                            for (int sk = 0; sk < seed_density; sk++) {  // for every seed point we need to create (si is 'seed index' i)
                                // printf("%f, %f, %f\n", curr_x1[0], curr_x1[1], curr_x1[2]);  // TESTING
                                double curr_x1[3] = {  // get current seed point
                                    x1_min[0] + (si*(1.0/seed_density)),
                                    x1_min[1] + (sj*(1.0/seed_density)),
                                    x1_min[2] + (sk*(1.0/seed_density))
                                };
                                double *curr_x0 = voxelToWorld(affine, curr_x1);  // convert seed point in voxel coords to world coords
                                for (int a = 0; a < 3; a++) {  // copy current seed into ret_seed array
                                    setElementDouble2D(ret_seeds, curr_x0[a], curr_seed, a, num_seeds, 3);
                                }
                                curr_seed += 1;  // increment current seed
                            }
                        }
                    }
                }
            }
        }
    }

    return(ret_seeds);
}


/*
get the side length of a voxel. IMPORTANT: assumes voxels are isotropic
inputs:
- *affine: pointer to 4x4 double array, the voxel-to-world affine transformation matrix
output:
- double, the side length of the voxel
*/
double voxelSideLength(double *affine) {
    double x1_ax1[4][1] = {  // (1, 0, 0) in voxel coordinates (x1)
        {1},
        {0},
        {0},
        {1}
    };
    double x1_ax0[4][1] = {  // (0, 0, 0) in voxel coordinates (x1)
        {0},
        {0},
        {0},
        {0}
    };
    double *x0_ax1 = matMult(affine, x1_ax1, 4, 4, 4, 1);  // x1=(1, 0, 0) converted to world coordinates (x0)
    double *x0_ax0 = matMult(affine, x1_ax0, 4, 4, 4, 1);  // x0=(0, 0, 0) converted to world coordinates (x1)

    double side_vec[3] = {
        x0_ax1[0] - x0_ax0[0],
        x0_ax1[1] - x0_ax0[1],
        x0_ax1[2] - x0_ax0[2]
    };

    double ret_val = sqrt(pow(side_vec[0], 2) + pow(side_vec[1], 2) + pow(side_vec[2], 2));
    free(x0_ax0); free(x0_ax1);  // free pointers
    return(ret_val);
}


/*
ALWAYS call this to convert world coordinates to voxel coordinates. this ensures
    that rounding errors are handled consistently
inputs:
- *v_to_w_affine: pointer to 4x4 double array, the voxel-to-world (A_0_2) affine matrix. internally,
    the inverse of this is taken
- *x0: pointer to length-3 double array, the point of interest in world coordinates (x0)
- *out_remainders: pointer to length-3 double array, the remainder of each voxel index after accounting
    for floating point rounding. necessary during trilinear interpolation
output:
- pointer to length-3 int array, the voxel indices (x2) of the point (x0). remainder values
    are stored in out_remainders
*/
int *getVoxelIndices(double *v_to_w_affine, double *x0, double *out_remainders) {
    double *A_2_0 = inv(v_to_w_affine, 4);  // get world-to-voxel matrix (inverse of A_0_2)
    double x_0_aff[4][1] = {
        {x0[0]},
        {x0[1]},
        {x0[2]},
        {1.0}
    };
    double *x_2_aff = matMult(A_2_0, x_0_aff, 4, 4, 4, 1);
    double x_2[3];  // voxel coordinates in double form
    x_2[0] = x_2_aff[0];
    x_2[1] = x_2_aff[1];
    x_2[2] = x_2_aff[2];
    
    // converting double voxel indices to integers
    int *ret_indices = malloc(3 * sizeof(int));
    for (int i = 0; i < 3; i++) {  // for each axis
        double voxel = x_2[i] + 0.5;
        if (fabs(voxel - floor(voxel)) < pow(10, -6)) {  // is marginally greater than integer, round down
            ret_indices[i] = (int) (voxel + 0.1);  // give 0.1 increment so necessarily rounds down (i.e., 9.0000001 = floor(9.0000001 + 0.1) = 9)
            out_remainders[i] = 0.0;
        } else if (fabs(voxel - floor(voxel) - 1.0) < pow(10, -6)) {  // is marginally less than integer, round up to voxel
            ret_indices[i] = ((int) (voxel - 0.1)) + 1;  // decrease by 0.1 so necessarily rounds up (i.e., 9.999999 = floor(9.999999 - 0.1) + 1 = 10)
            out_remainders[i] = 0.0;
        } else {  // 'voxel' is true decimal, not rounded at all
            ret_indices[i] = (int) voxel;  // effectively take floor of value
            out_remainders[i] = voxel - floor(voxel);
        }
    }
    // printf("\t\t\t\t\t%d, %d, %d\n", ret_indices[0], ret_indices[1], ret_indices[2]);
    free(A_2_0); free(x_2_aff);  // free pointers
    return(ret_indices);
}


/*
convert a point in voxel coordinates (x1) to world coordinates (x0), given the
    voxel-to-world affine matrix (A_0_1)
inputs:
- *affine: pointer to 4x4 double array, the voxel-to-world affine transformation matrix
- *x1: pointer to length-3 double array, the point in voxel coordinates
output:
- returns pointer to length-3 double array, the point in world coordinates (x0)
*/
double *voxelToWorld(double *affine, double *x1) {
    double x1_aff[4][1] = {
        {x1[0]},
        {x1[1]},
        {x1[2]},
        {1.0}
    };
    double *x0_aff = matMult(affine, x1_aff, 4, 4, 4, 1);
    double *ret_x0 = malloc(3 * sizeof(double));
    ret_x0[0] = x0_aff[0];
    ret_x0[1] = x0_aff[1];
    ret_x0[2] = x0_aff[2];

    free(x0_aff);  // free pointer
    return(ret_x0);
}


/*
determine whether a point is inside a mask
inputs:
- *mask: pointer to 3D int array, the mask
- *affine: pointer to 4x4 double array, the voxel-to-world affine matrix
- *p: pointer to length-3 double array, a point in world coordinates
- dim1: int, the size of the first dimension of the mask
- dim2: int, the size of the second dimension of the mask
- dim3: int, the size of the third dimension of the mask
output:
- int, 1 if the point is in the mask, otherwise 0
*/
int pointInMask(int *mask, double *affine, double *p, int dim1, int dim2, int dim3) {
    // printf("\taaaaaaah1\n");  // TESTING
    double *trash = malloc(3 * sizeof(int));  // remainder of world-to-voxel conversion, dont care about but function requires
    // printf("\taaaaaaah2\n");  // TESTING
    int *x2 = getVoxelIndices(affine, p, trash);
    // printf("\taaaaaaah3 (%d, %d, %d) (%d, %d, %d)\n", x2[0], x2[1], x2[2], dim1, dim2, dim3);  // TESTING
    if (getElementInt3D(mask, x2[0], x2[1], x2[2], dim1, dim2, dim3) == 1) {
        // printf("\taaaaah4\n");
        free(trash); free(x2);  // free pointers
        return(1);
    } else {
        // printf("\taaaaah4\n");
        free(trash); free(x2);  // free pointers
        return(0);
    }
}


/*
get the tracking ('inset') mask for a given mask. Every voxel in the tracking mask contains all 8 vertices
    inside the original mask, such that trilinear interpolation is possible
inputs:
- *mask: pointer to 3D int array, the original mask
- *affine: pointer to 4x4 double array, the original affine
- dim1: int, size of the first dimension of the mask
- dim2: int, size of the second dimension of the mask
- dim3: int, size of the third dimension of the mask
output:
- returns pointer to a newly allocated 3D double array, the tracking mask. note that the tracking mask has
    shape: tracking_mask.shape = mask.shape - (1, 1, 1)
*/
int *getTrackingMask(int *mask, double *affine, int dim1, int dim2, int dim3) {
    double *tracking_affine = getTrackingAffine(affine);  // getting offset tracking affine
    int *tracking_mask = malloc((dim1 - 1) * (dim2 - 1) * (dim3 - 1) * sizeof(int));
    for (int i = 0; i < (dim1 - 1); i++) {
        for (int j = 0; j < (dim2 - 1); j++) {
            for (int k = 0; k < (dim3 - 1); k++) {  // for every voxel in new tracking mask
                int offsets[8][3] = {  // corner offsets of first voxel in new mask (offset centers of original mask 8 voxels)
                    {0, 0, 0},
                    {0, 0, 1},
                    {0, 1, 0},
                    {1, 0, 0},
                    {0, 1, 1},
                    {1, 0, 1},
                    {1, 1, 0},
                    {1, 1, 1}
                };
                int curr_voxel_in_mask = 1;  // value of current voxel in new mask
                for (int a = 0; a < 8; a++) {  // for each corner (center of old voxels)
                    if (getElementInt3D(mask, i + offsets[a][0], j + offsets[a][1], k + offsets[a][2], dim1, dim2, dim3) == 0) {
                        curr_voxel_in_mask = 0;  // if any corners not in old mask, dont include this voxel in new mask
                    }
                }
                setElementInt3D(tracking_mask, curr_voxel_in_mask, i, j, k, dim1 - 1, dim2 - 1, dim3 - 1);
            }
        }
    }
    free(tracking_affine);  // free pointer
    return(tracking_mask);
}


/*
convert original affine to "tracking" affine, which is inset by half a voxel. If A_0_1
    is the original affine matrix, an offset transform A_1_2 is defined to yield the
    new affine matrix, A_0_2 = A_0_1 @ A_1_2. Note that with this transform, each dimension
    in the mask decreases by 1, and each voxel in the tracking mask is included only if
    the vertex of each voxel is in the original mask.
inputs:
- *affine: pointer to 2D double array, the original affine matrix
output:
- returns a pointer to a newly allocated 4x4 double array, the "tracking" affine matrix
*/
double *getTrackingAffine(double *affine) {
    double transform_affine[4][4] = {  // A_1_2, gives new_affine = old_affine @ A_1_2
        {1, 0, 0, 0.5},
        {0, 1, 0, 0.5},
        {0, 0, 1, 0.5},
        {0, 0, 0, 1}
    };
    double *new_affine = matMult(affine, transform_affine, 4, 4, 4, 4);
    return(new_affine);
}


/*
determine if a mask is continuous. A mask is defined as continuous if, for every voxel
    in the mask, there is at least one face-adjacent mask that is also in the mask
inputs:
- *mask: pointer to 3D int array, the mask
- dim1: int, the size of the first dimension of the mask
- dim2: int, the size of the second dimension of the mask
- dim3: int, the size of the third dimension of the mask
output:
- returns an int, 1 if the mask is continuous, otherwise 0
*/
int isContinuousVolume(int *mask, int dim1, int dim2, int dim3) {
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            for (int k = 0; k < dim3; k++) {
                if (getElementInt3D(mask, i, j, k, dim1, dim2, dim3) == 1) {  // current voxel is in mask
                    int curr_voxel_connected = 0;
                    int offsets[6][3] = {
                        {-1, 0, 0},
                        {1, 0, 0},
                        {0, -1, 0},
                        {0, 1, 0},
                        {0, 0, -1},
                        {0, 0, 1}
                    };
                    for (int a = 0; a < 6; a++) {  // check all face adjacent voxels
                        int curr_i = i + offsets[a][0];
                        int curr_j = j + offsets[a][1];
                        int curr_k = k + offsets[a][2];
                        if (getElementInt3D(mask, curr_i, curr_j, curr_k, dim1, dim2, dim3) == 1) {  // at least one face adjacent voxel in mask
                            curr_voxel_connected = 1;
                        }
                    }
                    if (curr_voxel_connected == 0) {  // no face adjacent voxels in mask
                        return(0);
                    }
                }
            }
        }
    }
    return(1);  // hasn't returned yet, so all voxels are connected (face-adjacent to at least one other voxel in mask)
}


/*
get the total sum of a mask (number of voxels assuming binary 0/1 mask)
inputs:
- *mask: pointer to 3D int array, the mask
- dim1: int, size of first dimension of mask
- dim2: int, size of second dimension of mask
- dim3: int, size of third dimension of mask
output:
- returns an int, the sum of the mask
*/
int getMaskSum(int *mask, int dim1, int dim2, int dim3) {
    int ret_sum = 0;
    for (int i = 0; i < dim1; i++) {  // across x voxels
        for (int j = 0; j < dim2; j++) {  // across y voxels
            for (int k = 0; k < dim3; k++) {  // across z voxels
                ret_sum += getElementInt3D(mask, i, j, k, dim1, dim2, dim3);
            }
        }
    }
    return(ret_sum);
}