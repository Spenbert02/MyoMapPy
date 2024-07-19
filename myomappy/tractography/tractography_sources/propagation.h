#pragma once
#include "VecField.h"
#include "Streamline.h"

// tracking flags
int IN_TRACKING_REGION = 0;
int IN_BORDER_REGION = 1;
int OUT_OF_MASK = 2;
int TERMINATION_CRITERION_REACHED = 3;

// tractography method enumerations
int EULER = 0;
int RUNGE_KUTTA_4 = 1;

int propagateStream(struct VecField *vec_field, struct Streamline *streamline, int tracking_method, double step_size, double max_length, double max_prop_angle, double *sense_vec);
int rk4Direction(struct VecField *vec_field, double *p, double *sense_vec, double step_size, double *out);
int eulerDirection(struct VecField *vec_field, double *p, double *sense_vec, double *out);
