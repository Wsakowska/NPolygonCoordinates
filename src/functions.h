#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>

#define PI 3.14159265358979323846

typedef struct {
    double x;
    double y;
} Vector;

Vector rotate(Vector v, double theta);
int compareAsc(const void *a, const void *b);
int compareDesc(const void *a, const void *b);

/* Heap implementation (H4) */
void swap_double(double *a, double *b);
void heapify_up(double heap[], int index);
void heapify_down(double heap[], int heap_size, int index);
void heap_insert(double heap[], int *heap_size, double value);
double heap_extract_min(double heap[], int *heap_size);
double pairwise_heap_sum(double arr[], int n);

#endif
