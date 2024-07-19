#pragma once
#include <math.h>
#include "propagation.h"
#include "mat_math.c"
#include "mask.c"


/*
propagate a streamline, appending the propagated point (if appropriate) and returning the appropriate flag
inputs:
- *vec_field: pointer to VecField object, used to get directions and mask information
- *streamline: pointer to Streamline object, to propagate. length(streamline) must be >= 1. also,
    streamline[-1] must be in the tracking mask.
- tracking_method: int, enumeration for the tracking method to use (enumerations provided in propagation.h)
- step_size: double, the step size for tracking
- max_length: double, the maximum path length for tracking
- max_prop_angle: double, the maximum propagation angle in degrees
- *sense_vec: pointer to length-3 double, the sense vector at the last point of the streamline, to guide
    propagation. this is only used if length(streamline) = 1, otherwise the last segment of the streamline
    is used for sense checking
output:
- returns an int, the tracking flag (enumerated in propagation.h). If appropriate, the propagated point is
    appended to the streamline
*/
int propagateStream(struct VecField *vec_field, struct Streamline *streamline, int tracking_method,
                    double step_size, double max_length, double max_prop_angle, double *sense_vec) {

    // check if streamline is at max size
    if (streamline->length >= MAX_SL_POINTS) {
        return(TERMINATION_CRITERION_REACHED);
    }

    // printf("groovy\n"); // TESTING

    // get sense vector
    double v_sense[3];
    if (streamline->length == 1) {  // need to use the provided sense vector
        v_sense[0] = sense_vec[0];
        v_sense[1] = sense_vec[1];
        v_sense[2] = sense_vec[2];
    } else {  // use the last segment of the streamline as the sense vector
        double p_neg1[3], p_neg2[3];
        // printf("\t\tgroovy1\n");  // TESTING
        streamline->getPoint(streamline, streamline->length - 1, p_neg1);
        // printf("\t\tgroovy2\n");  // TESTING
        streamline->getPoint(streamline, streamline->length - 2, p_neg2);
        // printf("\t\tgroovy3\n");  // TESTING
        double *seg = vecSubtract(p_neg1, p_neg2);
        v_sense[0] = seg[0];
        v_sense[1] = seg[1];
        v_sense[2] = seg[2];
        free(seg);  // free seg pointer
    }

    // get the tracking direction at the current point, and the corresponding tracking flag
    double dir_vec[3];  // vector at streamline[-1] to use for propagation, accounting for sense
    double p[3];  // streamline [-1]
    int dir_get_flag;
    streamline->getPoint(streamline, streamline->length - 1, p);
    if (tracking_method == RUNGE_KUTTA_4) {
        dir_get_flag = rk4Direction(vec_field, p, v_sense, step_size, dir_vec);
    } else if (tracking_method == EULER) {
        dir_get_flag = eulerDirection(vec_field, p, v_sense, dir_vec);
    }

    if (dir_get_flag == OUT_OF_MASK) {  // intermediate point (i.e., midpoint in RK tracking) was out of tracking mask. cant get direction, need to terminate
        return(OUT_OF_MASK);
    }

    // propagate point
    double *v_step = vecConstMult(step_size, norm(dir_vec));  // complete segment step (double ensuring dir vec is normalized)
    double *p_next = vecAdd(p, v_step);

    // check if appending would reach max length termination criterion (max_length = -1.0 means no max length criteria)
    // length of streamline after appending = length(stream) + step_size must be <= max_length, rearrange to get
    // criteria that length(stream) <= max_length - step_size to continue
    // printf("groovy\n");  // TESTING
    // printf("%d%d ", fabs(max_length + 1.0) < pow(10.0, -6.0), streamline->euclideanLength(streamline) > (max_length - step_size));
    if (!(fabs(max_length + 1.0) < pow(10.0, -6.0)) && (streamline->euclideanLength(streamline) > (max_length - step_size))) {
        return(TERMINATION_CRITERION_REACHED);
    }

    // check if appending would reach propagation angle termination criterion (max_prop_angle = -1.0 means no max prop angle criteria)
    if (streamline->length > 1) {  // need streamline to have one segment
        double *seg_neg1 = norm(v_step);  // normalized segment vector for new segment
        double p_neg1[3];  // last point in streamline before appending new point
        double p_neg2[3];  // second to last point in streamline before appending new point
        // printf("\t\tgroovy1\n");  // TESTING
        // printf("\t\tSL length %d\n", streamline->length);  // TESTING
        streamline->getPoint(streamline, streamline->length - 1, p_neg1);
        // printf("\t\tgroovy2\n");  // TESTING
        streamline->getPoint(streamline, streamline->length - 2, p_neg2);
        // printf("\t\tgroovy3\n");  // TESTING
        double *seg_neg2 = norm(vecSubtract(p_neg1, p_neg2));  // normalized segment vector for last segment of streamline before appending new point
        // calculate prop angle (in degrees), where theta = acos(a.b/(||a||*||b||)), and ||a|| = ||b|| = 1, such that theta = acos(a.b)
        double prop_angle_deg = acos(dot(seg_neg1, seg_neg2, 3)) * (180.0 / M_PI);
        if (!(fabs(max_prop_angle +1.0) < pow(10.0, -6.0)) && (prop_angle_deg > max_prop_angle)) {  // if considering prop angle and greater than max
            return(TERMINATION_CRITERION_REACHED);
        }
        
        free(seg_neg1); free(seg_neg2);  // freeing pointers
    }

    // printf("groovy ahhh\n");  // TESTING

    // termination criterion not reached, so determine where propagated point is, return appropriate flag and append point
    // to streamline if appropriate
    if (pointInMask(vec_field->tracking_mask, vec_field->tracking_affine, p_next,
                    vec_field->scan_dim1-1, vec_field->scan_dim2-1, vec_field->scan_dim3-1) == 1) {  // point is inside tracking mask (M2). append and continue tracking
        // printf("\tgroovy1\n");  // TESTING
        streamline->append(streamline, p_next);
        // printf("\tgroovy1\n");  // TESTING

        free(v_step); free(p_next);// free pointers
        return(IN_TRACKING_REGION);
    } else if (pointInMask(vec_field->mask, vec_field->affine, p_next,
               vec_field->scan_dim1, vec_field->scan_dim2, vec_field->scan_dim3) == 1) {  // point is outside M2 but inside M1. append but stop tracking
        // printf("\tgroovy2\n");  // TESTING
        streamline->append(streamline, p_next);
        free(v_step); free(p_next);// free pointers
        return(IN_BORDER_REGION);
    } else {  // point is outside both M2 and M1. dont append and stop tracking
        // printf("\tgroovy3\n");  // TESTING
        free(v_step); free(p_next);// free pointers
        return(OUT_OF_MASK);
    }
}


/*
get the direction for rk4 tracking at the current point. because the Runge-Kutta method used
    iterative propagation, it is possible for the internally propagated midpoints to extend beyond
    the tracking mask. in such cases where one of the propagated midpoints extends outside the tracking mask,
    the OUT_OF_MASK flag is returned and tracking should be stopped (OUT_OF_MASK returned by propagateStream
    function). if all propagated midpoints remain inside the tracking mask, the IN_TRACKING_REGION flag is
    returned, but it is not gauranteed that the propagated point itself will be in the tracking mask, just
    that the midpoints were.
*/
int rk4Direction(struct VecField *vec_field, double *p, double *sense_vec, double step_size, double *out) {
    return(-99);
}


/*
get the direction for euler tracking at the current point. the given point should ALWAYS be in the tracking
    mask, so the flag returned is always IN_TRACKING_REGION (no intermediate steps). A flag is only returned
    for consistency between tracking method direction "getters"
inputs:
- *vec_field: pointer to VecField object to use for direction getting
- *p: pointer to length-3 double array, the point at which to get the direction
- *sense_vec: pointer to length-3 double array, the sense vector to ensure same sense tracking
- *out: pointer to length-3 double array, the output direction vector
output:
- normalized direction vector is stored in *out
*/
int eulerDirection(struct VecField *vec_field, double *p, double *sense_vec, double *out) {
    vec_field->getValue(vec_field, p, out);  // out is a normalized vector
    if (dot(out, sense_vec, 3) < 0) {  // point obtained is antiparallel to sense vec, flip 'out' vector
        out[0] = -1.0 * out[0];
        out[1] = -1.0 * out[1];
        out[2] = -1.0 * out[2];
    }
    return(IN_TRACKING_REGION);
}