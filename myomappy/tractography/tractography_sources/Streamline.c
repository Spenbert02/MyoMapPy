#pragma once
#include <stdlib.h>
#include <math.h>
#include "Streamline.h"
#include "mat_math.c"


/*
get the euclidean distance of a streamline. the length of each segment is summed and returned
inputs:
- *this: pointer to Streamline object, to get the euclidean distance of
output:
- returns double, the euclidean length of the streamline
*/
double euclideanLength(struct Streamline *this) {
    double ret_length = 0;  // total length to return

    if (this->length <= 1) {  // if 0 or 1 points, has length 0
        return(0.0);
    }

    for (int i = 0; i < this->length - 1; i++) {  // for every segment in the streamline
        // printf("euc loop %d\n", i);  // TESTING
        double p0[3], p1[3];
        this->getPoint(this, i, p0);
        this->getPoint(this, i+1, p1);
        double *segment = vecSubtract(p0, p1);
        ret_length += sqrt(pow(segment[0], 2) + pow(segment[1], 2) + pow(segment[2], 2));
        free(segment);  // free pointer
    }
    return(ret_length);
}


/*
append a point to the streamline.
inputs:
- *this: pointer to Streamline object to append to
- p: length-3 double array, the point to append to the streamline
output:
- point p is appended to the end of the Streamline
*/
void append(struct Streamline *this, double p[3]) {
    for (int i = 0; i < 3; i++) {  // append array values to new point in array
        int offset = (this->length * 3) + i;
        *((double *)this->elements + offset) = p[i];
    }
    this->length = this->length + 1;  // increment length
}


/*
get point at a given index in the array
inputs:
- *this: pointer to Streamline object to get point from
- index: int, the index from which to get the point
- out: length-3 double array, the array where the point at the given index will be stored
output:
- the values of the point at the given index are stored in out. if the index is out of range,
    out will be set to {-99999.0, -99999.0, -99999.0}
*/
void getPoint(struct Streamline *this, int index, double out[3]) {
    if (index >= this->length) {  // out of bounds access
        out[0] = -99999.0;
        out[1] = -99999.0;
        out[2] = -99999.0;
        exit(1);
    } else {
        for (int i = 0; i < 3; i++) {  // if valid index, access elements
            int offset = (index * 3) + i;
            out[i] = *((double *)this->elements + offset);
        }
    }
}


/*
print the streamline values to the terminal
*/
void print(struct Streamline *this){
    printf("--- streamline ---\n");
    for (int i = 0; i < this->length; i++) {
        printf("%f, %f, %f\n", *((double *)this->elements + (i * 3) + 0), *((double *)this->elements + (i * 3) + 1), *((double *)this->elements + (i * 3) + 2));
    }
    printf("------------------\n");
}


/*
pseudo-destructor, to manually free the memory taken up by a streamline
inputs:
- *this: pointer to Streamline object to free
output:
- all pointers in the Streamline object are freed
*/
void freeSL(struct Streamline *this) {
    free(this->elements);
    // free(this->freeSL);
}


/*
Constructor for Streamline class
*/
static struct Streamline newSL() {
    double *elements = malloc(MAX_SL_POINTS * 3 * sizeof(double));  // allocate max streamline length memory
    return (struct Streamline){.length=0, .elements=elements, .append=&append, .getPoint=&getPoint, .print=&print, .euclideanLength=&euclideanLength, .freeSL=&freeSL};
}


const struct StreamlineClass Streamline={.newSL=&newSL};

