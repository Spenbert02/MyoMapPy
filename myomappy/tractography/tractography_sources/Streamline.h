#pragma once


struct Streamline {
    int length;  // number of points in streamline
    int arr_length;  // first dimension of 'elements' array (increases as elements are appended)
    double *elements;  // pointer to 2D array storing points
    void (*append)(struct Streamline *this, double *p);  // append point to streamline
    void (*getPoint)(struct Streamline *this, int index, double out[3]);  // get point from streamline
    void (*print)(struct Streamline *this);  // print the values of this streamline (may be very long if streamline is long)
};


extern const struct StreamlineClass {
    struct Streamline (*newSL)();
} Streamline;
