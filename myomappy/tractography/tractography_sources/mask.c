#pragma once
#include <math.h>
#include "mask.h"
#include "array_accessing.c"
#include "mat_math.c"


/*
determine wheter a point is inside a mask
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
    // for each voxel in mask, check if point is within. if not within any, return false
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            for (int k = 0; k < dim3; k++) {
                if (getElementInt3D(mask, i, j, k, dim1, dim2, dim3) == 1) {  // if current voxel in mask
                    if (pointInVoxel(affine, i, j, k, p)) {  // if point in current voxel, return true (1)
                        return(1);
                    }
                }
            }
        }
    }
    return(0);  // hasn't returned, so point not in any voxels
}


/*
determine whether a point is in a specified voxel of a mask. points on the voxel boundary
    (on one of the faces) are considered inside the voxel
inputs:
- *affine, pointer to 2D double array, the voxel-to-world affine transform
- x_vox: int, the index of the x-voxel of interest
- y_vox: int, the index of the y-voxel of interest
- z_vox: int, the index of the z-voxel of interest
- *p: pointer to length-3 double array, the point in world coordinates to check
output:
- int, 1 if the point is in the voxel, otherwise 0
*/
int pointInVoxel(double *affine, int x_vox, int y_vox, int z_vox, double *p) {
    double *A_1_0 = inv(affine, 4);  // inverse of affine
    double p_0_aff[4][1] = {{p[0]}, {p[1]}, {p[2]}, {1.0}};
    double *p_1_aff = matMult(A_1_0, p_0_aff, 4, 4, 4, 1);
    double p_1[3] = {  // point converted to voxel coordinates
        getElementDouble2D(p_1_aff, 0, 0, 4, 1),
        getElementDouble2D(p_1_aff, 1, 0, 4, 1),
        getElementDouble2D(p_1_aff, 2, 0, 4, 1),
    }; 

    int vox_borders[3][2] = {  // voxel borders
        {x_vox - 0.5, x_vox + 0.5},
        {y_vox - 0.5, y_vox + 0.5},
        {z_vox - 0.5, z_vox + 0.5}
    };
    
    // test for "on" voxel wall case (10^-6 approximation). if on voxel wall, return True
    for (int i = 0; i < 3; i++) {  // for each axis
        for (int j = 0; j < 2; j++) {  // for each border (min, max)
            if (fabs(p_1[i] - vox_borders[i][j]) < pow(10, -6)) {  // if within 10^-6, consider on voxel face and return true (1)
                return(1);
            }
        }
    }

    // check if inside voxel
    for (int i = 0; i < 3; i++) {  // for each axis
        if (!(p[i] > vox_borders[i][0] && p[i] < vox_borders[i][1])) {  // if not inside voxel
            return(0);
        }
    }
    return(1);  // if hasn't returned, point is within each axis voxel range
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