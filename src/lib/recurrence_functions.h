#ifndef _RECURRENCE_FUNCTIONS_H
#define _RECURRENCE_FUNCTIONS_H

#include <stdbool.h>

// -----------------------------------------------------------------------------
// Distance‐matrix & embedding routines
// -----------------------------------------------------------------------------

typedef double (*distance_function)(double[], double[], int);

void embed_time_series(double **embedded_time_series,
                       double *time_series,
                       int n,
                       int dim,
                       int tau);

void manhattan_distance(double **embedded_time_series,
                        int n,
                        int dim,
                        double **D);

void euclidean_distance(double **embedded_time_series,
                        int n,
                        int dim,
                        double **D);

void supremum_distance(double **embedded_time_series,
                       int n,
                       int dim,
                       double **D);

void modulated_distance_2pi(double **embedded_time_series,
                            int n,
                            int dim,
                            double **D);

void modulated_distance_pi(double **embedded_time_series,
                           int n,
                           int dim,
                           double **D);

void calculate_distance_matrix(double **embedded_time_series,
                               int n,
                               int dim,
                               double **D,
                               const char *norm);

// -----------------------------------------------------------------------------
// Recurrence‐matrix construction & threshold adjustment
// -----------------------------------------------------------------------------

void compute_recurrence_matrix(double **D,
                               int n,
                               int **R,
                               double threshold);

double adjust_threshold_via_recurrence_rate(double **D,
                                            int n,
                                            int **R,
                                            double desired_percentage,
                                            int *status,
                                            double *actual_percentage);

// -----------------------------------------------------------------------------
// Basic RQA quantifiers
// -----------------------------------------------------------------------------

double compute_recurrence_rate(int **R, int n);
double compute_determinism(int **R, int n);
double compute_laminarity(int **R, int n);
double compute_trapping_time(int **R, int n);

// -----------------------------------------------------------------------------
// Diagonal‐length distributions
// -----------------------------------------------------------------------------

int* calculate_diagonal_frequency_distribution(int **R, int n);
double compute_average_diagonal_length(int **R, int n);
double compute_longest_diagonal_length(int **R, int n);

// -----------------------------------------------------------------------------
// Vertical‐length distributions
// -----------------------------------------------------------------------------

int* calculate_vertical_frequency_distribution(int **R, int n);
double compute_average_vertical_length(int **R, int n);
double compute_longest_vertical_length(int **R, int n);

// -----------------------------------------------------------------------------
// Entropy measures
// -----------------------------------------------------------------------------

double compute_entropy_diagonal(int **R, int n);
double compute_entropy_vertical(int **R, int n);
double compute_new_entropy_diagonal(int **R, int n);
double compute_new_entropy_vertical(int **R, int n);

// -----------------------------------------------------------------------------
// Divergence / inverse longest‐line
// -----------------------------------------------------------------------------

double compute_divergence(int **R, int n);
double wrap_longest_diagonal_length(int **matrix, int size);
double wrap_longest_vertical_length(int **matrix, int size);

// -----------------------------------------------------------------------------
// Skeletonization (Kraemer & Marwan 2019)
// -----------------------------------------------------------------------------

// Represents a vertical run (length, start_row, col) in the “close‐returns” map.
typedef struct {
    int length;
    int start_row;
    int col;
} Line;

// Compare two Line structs by descending length (for qsort).
int compare_lines_desc(const void *pA, const void *pB);

// Build a “close‐returns” map CR[N][N+1] from RM[N][N]’s lower triangle.
void build_close_returns_map(int **RM,
                             int N,
                             int **CR);

// Scan CR (size N×M) column‐by‐column to extract all vertical runs.
// Allocates *out_lines (caller must free); sets *out_count to the number of runs.
void build_line_list(int **CR,
                     int N,
                     int M,
                     Line **out_lines,
                     int *out_count);

// From the full list of runs “lines” (size = linecount), select a final, non‐overlapping
// subset by sorting by descending length, accepting runs that are still present in CR_copy,
// zeroing out their pixels, and recursively deleting any adjacent runs that touch them.
// Allocates *out_final_lines (caller must free); sets *out_final_count.
void get_final_line_list(Line *lines,
                         int linecount,
                         int **CR_copy,
                         int N,
                         int M,
                         Line **out_final_lines,
                         int *out_final_count);

// Recursively delete any run in CR_copy that touches pixel (r,c).
void delete_line_recursive(int **CR,
                           int N,
                           int M,
                           int r,
                           int c);

// From the selected “final” runs (size = final_count), reconstruct SK[N][N] so that each
// run maps back to “thickness‐1” diagonals in the original RP. Ensures SK[i][i] = 1.
void build_skeletonized_RP_from_runs(Line *final,
                                     int final_count,
                                     int N,
                                     int **SK);

// Top‐level thinning: produce SK[N][N] from RM[N][N] by applying the full
// Kraemer & Marwan (2019) skeletonization pipeline.
void skeletonize_recurrence_matrix(int **RM, int N, int **SK);

#endif  // _RECURRENCE_FUNCTIONS_H
