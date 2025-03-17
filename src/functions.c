#include "functions.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

Vector rotate(Vector v, double theta) {
    Vector result;
    result.x = v.x * cos(theta) - v.y * sin(theta);
    result.y = v.x * sin(theta) + v.y * cos(theta);
    return result;
}

int compareAsc(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    else if (da > db) return 1;
    else return 0;
}

int compareDesc(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da > db) return -1;
    else if (da < db) return 1;
    else return 0;
}

void swap_double(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

void heapify_up(double heap[], int index) {
    while (index > 0) {
        int parent = (index - 1) / 2;
        if (fabs(heap[index]) < fabs(heap[parent])) {
            swap_double(&heap[index], &heap[parent]);
            index = parent;
        } else {
            break;
        }
    }
}

void heapify_down(double heap[], int heap_size, int index) {
    while (2 * index + 1 < heap_size) {
        int left = 2 * index + 1;
        int right = left + 1;
        int smallest = index;
        if (left < heap_size && fabs(heap[left]) < fabs(heap[smallest]))
            smallest = left;
        if (right < heap_size && fabs(heap[right]) < fabs(heap[smallest]))
            smallest = right;
        if (smallest != index) {
            swap_double(&heap[index], &heap[smallest]);
            index = smallest;
        } else {
            break;
        }
    }
}

void heap_insert(double heap[], int *heap_size, double value) {
    heap[*heap_size] = value;
    (*heap_size)++;
    heapify_up(heap, *heap_size - 1);
}

double heap_extract_min(double heap[], int *heap_size) {
    double min = heap[0];
    heap[0] = heap[*heap_size - 1];
    (*heap_size)--;
    heapify_down(heap, *heap_size, 0);
    return min;
}

double pairwise_heap_sum(double arr[], int n) {
    double *heap = malloc(n * sizeof(double));
    if (!heap) {
    fprintf(stderr, "Memory allocation error in pairwise_heap_sum\n");
        exit(EXIT_FAILURE);
    }
    int heap_size = 0;
    for (int i = 0; i < n; i++) {
        heap[heap_size] = arr[i];
        heap_size++;
        heapify_up(heap, heap_size - 1);
    }
    while (heap_size > 1) {
        double a = heap_extract_min(heap, &heap_size);
        double b = heap_extract_min(heap, &heap_size);
        double s = a + b;
        heap_insert(heap, &heap_size, s);
    }
    double result = heap[0];
    free(heap);
    return result;
}
