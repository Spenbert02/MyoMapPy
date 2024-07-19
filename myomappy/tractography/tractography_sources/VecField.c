#pragma once
#include <math.h>
#include "VecField.h"
#include "mask.c"
#include "mat_math.c"
#include "array_accessing.c"

int num_seeds_out_voxel = 0;

/*
get value of the vector field at a point using trilinear interpolation
inputs:
- *this: pointer to this VecField object
- *x0: pointer to length-3 double array, the point of interest in world coordinates
- *out: pointer to length-3 double array, to hold the output value of the vector field at x0
output:
- out is set to the value of the vector field at x0
*/
void getValue(struct VecField *this, double *x0, double *out) {
    // getting base dirs[] indices, the tracking voxel corresponding to x0
    double v_frac[3];
    int *dirs_base_indices = getVoxelIndices(this->tracking_affine, x0, v_frac);

    // TESTING
    if (dirs_base_indices[0] == -1) {
        printf("negative\n");
    }
    if (getElementInt3D(this->tracking_mask, dirs_base_indices[0], dirs_base_indices[1], dirs_base_indices[2],
        this->scan_dim1, this->scan_dim2, this->scan_dim3) == 0) {
            num_seeds_out_voxel += 1;
        // printf("%d, %d, %d\n", dirs_base_indices[0], dirs_base_indices[1], dirs_base_indices[2]);
    }
    // printf("%d\n", num_seeds_out_voxel);

    // defining corner values of voxel
    double *corner_vecs = malloc(2 * 2 * 2 * 3 * sizeof(double));  // 2x2x2x3 array, axes are [i, j, k, 3] (vectors at each corner of tracking voxel)
    for (int di = 0; di < 2; di++) {  // dirs x-index offset
        for (int dj = 0; dj < 2; dj++) {  // dirs y-index offset
            for (int dk = 0; dk < 2; dk++) {  // dirs z-index offset
                for (int a = 0; a < 3; a++) {  // vector index (0, 1, 2)
                    double curr_val = getElementDouble4D(  // getting dirs[_+di, _+dj, _+dk, a]
                        this->dirs,
                        dirs_base_indices[0] + di,
                        dirs_base_indices[1] + dj,
                        dirs_base_indices[2] + dk,
                        a,
                        this->scan_dim1,
                        this->scan_dim2,
                        this->scan_dim3,
                        3
                    );
                    setElementDouble4D(corner_vecs, curr_val, di, dj, dk, a, 2, 2, 2, 3);  // setting corner vec value
                }
            }
        }
    }

    // interpolating and setting out value
    double *ret_vec = trilinearInterpolate(corner_vecs, v_frac);
    for (int i = 0; i < 3; i++) {
        out[i] = ret_vec[i];
    }

    free(dirs_base_indices); free(corner_vecs); free(ret_vec);  // free pointers
}


/*
given the corner vectors of a voxel and the fractional position inside the voxel,
    get the trilinearly interpolated vector at that position
inputs:
- *corner_vals: pointer to shape (2, 2, 2, 3) double array, the voxels at each corner.
    corner_vals[i][j][k][] corresponds to the vector value c_ijk (refer to notes)
- *v_frac: pointer to length-3 double array, the fractional position (x, y, z) in the voxel.
    example: c_000 corresponds to v_frac=[0, 0, 0], and c_111 corresponds to v_frac=[1, 1, 1]
output:
- returns a pointer to a length-3 double array, the vector field evaluated at the provided point,
    using trilinear interpolation
*/
double *trilinearInterpolate(double *corner_vals, double *v_frac) {
    double *ret_vec = malloc(3 * sizeof(double));
    for (int a = 0; a < 3; a++) {  // perform interpolation for each dimension/vector index (x=0, y=1, z=2)
        // getting corner values
        double c_000 = getElementDouble4D(corner_vals, 0, 0, 0, a, 2, 2, 2, 3);
        double c_001 = getElementDouble4D(corner_vals, 0, 0, 1, a, 2, 2, 2, 3);
        double c_010 = getElementDouble4D(corner_vals, 0, 1, 0, a, 2, 2, 2, 3);
        double c_100 = getElementDouble4D(corner_vals, 1, 0, 0, a, 2, 2, 2, 3);
        double c_011 = getElementDouble4D(corner_vals, 0, 1, 1, a, 2, 2, 2, 3);
        double c_110 = getElementDouble4D(corner_vals, 1, 1, 0, a, 2, 2, 2, 3);
        double c_101 = getElementDouble4D(corner_vals, 1, 0, 1, a, 2, 2, 2, 3);
        double c_111 = getElementDouble4D(corner_vals, 1, 1, 1, a, 2, 2, 2, 3);

        // interpolating along x
        double c_00 = (c_000 * (1-v_frac[0])) + (c_100 * v_frac[0]);
        double c_01 = (c_001 * (1-v_frac[0])) + (c_101 * v_frac[0]);
        double c_10 = (c_010 * (1-v_frac[0])) + (c_110 * v_frac[0]);
        double c_11 = (c_011 * (1-v_frac[0])) + (c_111 * v_frac[0]);

        // interpolating along y
        double c_0 = (c_00 * (1-v_frac[1])) + (c_10 * v_frac[1]);
        double c_1 = (c_01 * (1-v_frac[1])) + (c_11 * v_frac[1]);

        // finally interpolating along z
        double c = (c_0 * (1-v_frac[2])) + (c_1 * v_frac[2]);
        ret_vec[a] = c;
    }
    double *normed_ret_vec = norm(ret_vec);
    free(ret_vec);  // free pointer
    return(normed_ret_vec);  // return normalized vector
}


/*
Constructor for VecField class
*/
static struct VecField newVecField(double *dirs, int *mask, double *affine, int scan_dim1, int scan_dim2, int scan_dim3) {
    int *tracking_mask = getTrackingMask(mask, affine, scan_dim1, scan_dim2, scan_dim3);
    double *tracking_affine = getTrackingAffine(affine);
    return (struct VecField){
        .dirs=dirs,
        .mask=mask,
        .affine=affine,
        .scan_dim1=scan_dim1,
        .scan_dim2=scan_dim2,
        .scan_dim3=scan_dim3,
        .tracking_mask=tracking_mask,
        .tracking_affine=tracking_affine,
        .getValue=&getValue
    };
}


const struct VecFieldClass VecField={.newVecField=&newVecField};

