#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "functions.h"

int main() {
    int n_min, n_max, step;
    printf("Enter the minimum number of vertices (n_min): ");
    scanf("%d", &n_min);
    printf("Enter the maximum number of vertices (n_max): ");
    scanf("%d", &n_max);
    printf("Enter the iteration step: ");
    scanf("%d", &step);

    // Opening global output files
    FILE *global_errors = fopen("doc/global_errors.dat", "w");
    FILE *global_sorted_x_pos = fopen("doc/global_sorted_x_pos.dat", "w");
    FILE *global_sorted_x_neg = fopen("doc/global_sorted_x_neg.dat", "w");
    // We also keep files for y-values, which might be useful for further analysis
    FILE *global_sorted_y_pos = fopen("doc/global_sorted_y_pos.dat", "w");
    FILE *global_sorted_y_neg = fopen("doc/global_sorted_y_neg.dat", "w");

    if (global_errors == NULL || global_sorted_x_pos == NULL || global_sorted_x_neg == NULL ||
        global_sorted_y_pos == NULL || global_sorted_y_neg == NULL) {
        perror("Error opening global files");
        exit(EXIT_FAILURE);
    }

    // File headers
    fprintf(global_errors, "# n errorH2 errorH3 errorH4\n");
    fprintf(global_sorted_x_pos, "# n index x_value\n");
    fprintf(global_sorted_x_neg, "# n index x_value\n");
    fprintf(global_sorted_y_pos, "# n index y_value\n");
    fprintf(global_sorted_y_neg, "# n index y_value\n");

    // Iterative loop over n values
    for (int n = n_min; n <= n_max; n += step) {
        printf("\nProcessing n = %d\n", n);

        double theta = 2 * PI / n;
        Vector *w = malloc(n * sizeof(Vector));     // differential vectors
        Vector *v = malloc((n + 1) * sizeof(Vector)); // coordinates
        if (w == NULL || v == NULL) {
            fprintf(stderr, "Memory allocation error for n = %d.\n", n);
            exit(EXIT_FAILURE);
        }

        // H1: Calculate vertices of the n-gon
        v[0].x = 1.0;
        v[0].y = 0.0;

        w[0].x = cos(theta) - 1.0; 
        w[0].y = sin(theta);

        // Something like the shift from (1,0) to (1,0) + w[0]
        v[1].x = v[0].x + w[0].x;
        v[1].y = v[0].y + w[0].y;

        for (int i = 1; i < n; i++) {
            w[i] = rotate(w[i - 1], theta);
            v[i + 1].x = v[i].x + w[i].x;
            v[i + 1].y = v[i].y + w[i].y;
        }

        // H2: Summation of vectors w
        double sum_wx = 0.0, sum_wy = 0.0;
        for (int i = 0; i < n; i++) {
            sum_wx += w[i].x;
            sum_wy += w[i].y;
        }

        // H3: Ordered summation — sorting x and y coordinates
        int countPosX = 0, countNegX = 0;
        for (int i = 0; i < n; i++) {
            if (w[i].x >= 0)
                countPosX++;
            else
                countNegX++;
        }
        double *Xpos = malloc(countPosX * sizeof(double));
        double *Xneg = malloc(countNegX * sizeof(double));
        if (Xpos == NULL || Xneg == NULL) {
            fprintf(stderr, "Memory allocation error for Xpos/Xneg.\n");
            exit(EXIT_FAILURE);
        }
        int posIndex = 0, negIndex = 0;
        for (int i = 0; i < n; i++) {
            if (w[i].x >= 0)
                Xpos[posIndex++] = w[i].x;
            else
                Xneg[negIndex++] = w[i].x;
        }
        qsort(Xpos, countPosX, sizeof(double), compareAsc);
        qsort(Xneg, countNegX, sizeof(double), compareDesc);

        double SxPos = 0.0, SxNeg = 0.0;
        for (int i = 0; i < countPosX; i++) {
            SxPos += Xpos[i];
        }
        for (int i = 0; i < countNegX; i++) {
            SxNeg += Xneg[i];
        }
        double Sx = SxPos + SxNeg;

        int countPosY = 0, countNegY = 0;
        for (int i = 0; i < n; i++) {
            if (w[i].y >= 0)
                countPosY++;
            else
                countNegY++;
        }
        double *Ypos = malloc(countPosY * sizeof(double));
        double *Yneg = malloc(countNegY * sizeof(double));
        if (Ypos == NULL || Yneg == NULL) {
            fprintf(stderr, "Memory allocation error for Ypos/Yneg.\n");
            exit(EXIT_FAILURE);
        }
        posIndex = 0;
        negIndex = 0;
        for (int i = 0; i < n; i++) {
            if (w[i].y >= 0)
                Ypos[posIndex++] = w[i].y;
            else
                Yneg[negIndex++] = w[i].y;
        }
        qsort(Ypos, countPosY, sizeof(double), compareAsc);
        qsort(Yneg, countNegY, sizeof(double), compareDesc);

        double SyPos = 0.0, SyNeg = 0.0;
        for (int i = 0; i < countPosY; i++) {
            SyPos += Ypos[i];
        }
        for (int i = 0; i < countNegY; i++) {
            SyNeg += Yneg[i];
        }
        double Sy = SyPos + SyNeg;

        // H4: Pairwise heap summation
        double H4_SxPos = (countPosX > 0) ? pairwise_heap_sum(Xpos, countPosX) : 0.0;
        double H4_SxNeg = (countNegX > 0) ? pairwise_heap_sum(Xneg, countNegX) : 0.0;
        double H4_Sx = H4_SxPos + H4_SxNeg;

        double H4_SyPos = (countPosY > 0) ? pairwise_heap_sum(Ypos, countPosY) : 0.0;
        double H4_SyNeg = (countNegY > 0) ? pairwise_heap_sum(Yneg, countNegY) : 0.0;
        double H4_Sy = H4_SyPos + H4_SyNeg;

        // Calculate errors (norms)
        double errorH2 = sqrt(sum_wx * sum_wx + sum_wy * sum_wy);
        double errorH3 = sqrt(Sx * Sx + Sy * Sy);
        double errorH4 = sqrt(H4_Sx * H4_Sx + H4_Sy * H4_Sy);

        printf("For n = %d:\n  H2 error: %.15f\n  H3 error: %.15f\n  H4 error: %.15f\n",
               n, errorH2, errorH3, errorH4);

        // Write the results to global_errors.dat
        fprintf(global_errors, "%d %.15f %.15f %.15f\n", n, errorH2, errorH3, errorH4);

        // Write the sorted x values to global files
        for (int i = 0; i < countPosX; i++) {
            fprintf(global_sorted_x_pos, "%d %d %.15f\n", n, i, Xpos[i]);
        }
        for (int i = 0; i < countNegX; i++) {
            fprintf(global_sorted_x_neg, "%d %d %.15f\n", n, i, Xneg[i]);
        }

        // Write the sorted y values to global files
        for (int i = 0; i < countPosY; i++) {
            fprintf(global_sorted_y_pos, "%d %d %.15f\n", n, i, Ypos[i]);
        }
        for (int i = 0; i < countNegY; i++) {
            fprintf(global_sorted_y_neg, "%d %d %.15f\n", n, i, Yneg[i]);
        }

        // Free memory for the current iteration
        free(w);
        free(v);
        free(Xpos);
        free(Xneg);
        free(Ypos);
        free(Yneg);
    }

    // Close global output files
    fclose(global_errors);
    fclose(global_sorted_x_pos);
    fclose(global_sorted_x_neg);
    fclose(global_sorted_y_pos);
    fclose(global_sorted_y_neg);

    // ------------------- GENERATING PLOTS USING GNUPlot -------------------

    // Plot of global errors (H2, H3, H4) in logarithmic scale
    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp == NULL) {
        fprintf(stderr, "Error: Could not open a pipe to GNUPlot.\n");
    } else {
        fprintf(gp, "set terminal pngcairo enhanced size 800,600\n");
        fprintf(gp, "set output 'doc/global_errors.png'\n");
        fprintf(gp, "set title 'Summation errors for different n values (logarithmic scale)'\n");
        fprintf(gp, "set xlabel 'n'\n");
        fprintf(gp, "set ylabel 'Error (norm)'\n");
        fprintf(gp, "set grid\n");
        fprintf(gp, "set logscale y\n");
        fprintf(gp, "plot 'doc/global_errors.dat' using 1:2 with linespoints title 'H2', \\\n");
        fprintf(gp, "     'doc/global_errors.dat' using 1:3 with linespoints title 'H3', \\\n");
        fprintf(gp, "     'doc/global_errors.dat' using 1:4 with linespoints title 'H4'\n");
        fflush(gp);
        pclose(gp);
        printf("'doc/global_errors.png' chart generated.\n");
    }

    // Plot of globally sorted x values (linear scale)
    gp = popen("gnuplot -persistent", "w");
    if (gp == NULL) {
        fprintf(stderr, "Error: Could not open a pipe to GNUPlot.\n");
    } else {
        fprintf(gp, "set terminal pngcairo enhanced size 800,600\n");
        fprintf(gp, "set output 'doc/global_sorted_x.png'\n");
        fprintf(gp, "set title 'Sorted x values for different n'\n");
        fprintf(gp, "set xlabel 'Index'\n");
        fprintf(gp, "set ylabel 'x'\n");
        fprintf(gp, "set grid\n");
        fprintf(gp, "plot 'doc/global_sorted_x_pos.dat' using 2:3 with linespoints title 'positive x', \\\n");
        fprintf(gp, "     'doc/global_sorted_x_neg.dat' using 2:3 with linespoints title 'negative x'\n");
        fflush(gp);
        pclose(gp);
        printf("'doc/global_sorted_x.png' chart generated.\n");
    }

    // Plot of error ratios: H3/H2 and H4/H2 vs n (logarithmic scale)
    gp = popen("gnuplot -persistent", "w");
    if (gp == NULL) {
        fprintf(stderr, "Error: Could not open a pipe to GNUPlot.\n");
    } else {
        fprintf(gp, "set terminal pngcairo enhanced size 800,600\n");
        fprintf(gp, "set output 'doc/global_error_ratios.png'\n");
        fprintf(gp, "set title 'Error ratio: H3/H2 and H4/H2 vs n (logarithmic scale)'\n");
        fprintf(gp, "set xlabel 'n'\n");
        fprintf(gp, "set ylabel 'Error ratio'\n");
        fprintf(gp, "set grid\n");
        fprintf(gp, "set logscale y\n");
        fprintf(gp, "plot 'doc/global_errors.dat' using 1:($3/$2) with linespoints title 'H3/H2', \\\n");
        fprintf(gp, "     'doc/global_errors.dat' using 1:($4/$2) with linespoints title 'H4/H2'\n");
        fflush(gp);
        pclose(gp);
        printf("'doc/global_error_ratios.png' chart generated.\n");
    }

    // Plot of error differences: (H2 - H3) and (H2 - H4) vs n (logarithmic scale)
    gp = popen("gnuplot -persistent", "w");
    if (gp == NULL) {
        fprintf(stderr, "Error: Could not open a pipe to GNUPlot.\n");
    } else {
        fprintf(gp, "set terminal pngcairo enhanced size 800,600\n");
        fprintf(gp, "set output 'doc/global_error_differences.png'\n");
        fprintf(gp, "set title 'Error differences: H2 - H3 and H2 - H4 vs n (logarithmic scale)'\n");
        fprintf(gp, "set xlabel 'n'\n");
        fprintf(gp, "set ylabel 'Error difference'\n");
        fprintf(gp, "set grid\n");
        fprintf(gp, "set logscale y\n");
        fprintf(gp, "plot 'doc/global_errors.dat' using 1:($2-$3) with linespoints title 'H2 - H3', \\\n");
        fprintf(gp, "     'doc/global_errors.dat' using 1:($2-$4) with linespoints title 'H2 - H4'\n");
        fflush(gp);
        pclose(gp);
        printf("'doc/global_error_differences.png' chart generated.\n");
    }

    printf("\nProcess completed.\n");
    return 0;
}