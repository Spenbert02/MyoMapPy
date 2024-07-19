#pragma once

// hard-coded maximum number of points per streamline
int MAX_SL_POINTS = 1000;

struct Streamline {
    int length;  // number of points in streamline
    double *elements;  // pointer to 2D array storing points
    double (*euclideanLength)(struct Streamline *this);  // get the euclidean length of a streamline
    void (*append)(struct Streamline *this, double *p);  // append point to streamline
    void (*getPoint)(struct Streamline *this, int index, double out[3]);  // get point from streamline
    void (*print)(struct Streamline *this);  // print the values of this streamline (may be very long if streamline is long)
    void (*freeSL) (struct Streamline *this);  // free the memory stored in streamline
};

extern const struct StreamlineClass {
    struct Streamline (*newSL)();
} Streamline;
