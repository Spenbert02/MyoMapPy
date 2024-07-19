#pragma once

double *trilinearInterpolate(double *corner_vals, double *v_frac);

struct VecField {
    /* member variables */
    // pointer to 4D double array of shape (scan_dim1, scan_dim2, scan_dim3, 3), the vector field. dirs[i, j, k] is associated
    // with the center of voxel [i, j, k] of the original affine/mask
    double *dirs;

    // pointer to 3D int array of shape (scan_dim1, scan_dim2, scan_dim3), the binary mask.
    int *mask;

    // pointer to 4x4 double array, the voxel-to-world (A_0_1) affine transform matrix for the given mask and dirs
    double *affine;

    // int, the size of the corresponding dimension of the original scan
    int scan_dim1;
    int scan_dim2;
    int scan_dim3;

    // pointer to 3D double array of shape (scan_dim1 - 1, scan_dim2 - 1, scan_dim3 - 1), the slightly "inset" tracking mask.
    // internally set within the constructor. only points in this mask can be further propogated. see notes for more logic
    int *tracking_mask;

    // pointer to 4x4 double array, the updated affine for the tracking mask. internally set within the constructor. again,
    // see notes for more logic
    double *tracking_affine;

    /* member functions */
    // get value of the vector field at a given point in world coordinates, x0.
    void (*getValue)(struct VecField *this, double *x0, double *out);
};


extern const struct VecFieldClass {
    struct VecField (*newVecField)();
} VecField;
