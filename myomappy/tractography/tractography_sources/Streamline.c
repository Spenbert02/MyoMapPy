#pragma once
#include <stdlib.h>
#include "Streamline.h"


/*
append a point to the streamline.
inputs:
- *this: pointer to Streamline object to append to
- p: length-3 double array, the point to append to the streamline
output:
- point p is appended to the end of the Streamline
*/
void append(struct Streamline *this, double p[3]) {
    if (this->length == this->arr_length) {  // once elements array is maxed out, increase size of this->elements
        this->arr_length = this->arr_length * 2;  // double array size

        double *temp_elements = malloc(this->arr_length * 3 * sizeof(double));  // copy old elements into new array
        for (int i = 0; i < this->length; i++) {
            for (int j = 0; j < 3; j++) {
                int offset = (i * 3) + j;
                *((double *)temp_elements + offset) = *((double *)this->elements + offset);
            }
        }

        this->elements = temp_elements;
    }

    // append array values to new point in array
    for (int i = 0; i < 3; i++) {
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
Constructor for Streamline class
*/
static struct Streamline newSL() {
    double *elements = malloc(4 * 3 * sizeof(double));  // size 4 array by default
    return (struct Streamline){.length=0, .arr_length=4, .elements=elements, .append=&append, .getPoint=&getPoint, .print=&print};
}


const struct StreamlineClass Streamline={.newSL=&newSL};

