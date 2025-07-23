# include "common.h"
# include "recurrence_functions.h"
# include "../lib/auxiliary_functions.h" 

void embed_time_series(double **embedded_time_series, double *time_series, int n, int dim, int tau) 
{
    for (int i = 0; i < n - (dim - 1) * tau; i++) {
        for (int j = 0; j < dim; j++) {
            embedded_time_series[i][j] = time_series[i + j * tau];
        }
    }
}

void manhattan_distance(double **embedded_time_series, int n, int dim, double **D) 
{
    double temp;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            temp = 0.0;
            for (int k = 0; k < dim; k++) {
                temp += fabs(embedded_time_series[i][k] - embedded_time_series[j][k]);
            }
            D[i][j] = temp;
            D[j][i] = temp;
        }
    }
}

void euclidean_distance(double **embedded_time_series, int n, int dim, double **D) 
{
    double temp;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            temp = 0.0;
            for (int k = 0; k < dim; k++) {
                temp += pow(embedded_time_series[i][k] - embedded_time_series[j][k], 2);
            }
            D[i][j] = sqrt(temp);
            D[j][i] = sqrt(temp);
        }
    }
}

void supremum_distance(double **embedded_time_series, int n, int dim, double **D) 
{
    double max, temp;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            max = fabs(embedded_time_series[i][0] - embedded_time_series[j][0]);
            for (int k = 1; k < dim; k++) {
                temp = fabs(embedded_time_series[i][k] - embedded_time_series[j][k]);
                if (temp > max) {
                    max = temp;
                }
            }
            D[i][j] = max;
            D[j][i] = max;
        }
    }
}

void modulated_distance_2pi(double **embedded_time_series, int n, int dim, double **D) 
{
    double diff;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            diff = 0.0;
            for (int k = 0; k < dim; k++) {
                double temp = fmod(fabs(embedded_time_series[i][k] - embedded_time_series[j][k]), 2 * M_PI);
                if (temp > M_PI) {
                    temp = 2 * M_PI - temp;
                }
                diff += temp * temp;
            }
            D[i][j] = sqrt(diff);
            D[j][i] = sqrt(diff);
        }
    }
}

void modulated_distance_pi(double **embedded_time_series, int n, int dim, double **D) 
{
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double diff = 0.0;

            for (int k = 0; k < dim; k++) {
                double delta = embedded_time_series[i][k] - embedded_time_series[j][k];

                // Wrap difference to [-π, π)
                while (delta <= -M_PI) delta += 2 * M_PI;
                while (delta >   M_PI) delta -= 2 * M_PI;

                diff += delta * delta;
            }

            D[i][j] = sqrt(diff);
            D[j][i] = D[i][j];  // Symmetric
        }
    }
}

void calculate_distance_matrix(double **embedded_time_series, int n, int dim, double **D, const char *norm) 
{
    if (strcmp(norm, "manhattan") == 0) manhattan_distance(embedded_time_series, n, dim, D);
    else if (strcmp(norm, "euclidean") == 0) euclidean_distance(embedded_time_series, n, dim, D);
    else if (strcmp(norm, "supremum") == 0) supremum_distance(embedded_time_series, n, dim, D);
    else if (strcmp(norm, "modulated_2pi") == 0) modulated_distance_2pi(embedded_time_series, n, dim, D);
    else if (strcmp(norm, "modulated_pi") == 0) modulated_distance_pi(embedded_time_series, n, dim, D);
    else {
        fprintf(stderr, "⚠ WARNING: Unknown norm type '%s'. Using Euclidean as default.\n", norm);
        euclidean_distance(embedded_time_series, n, dim, D);
    }
}

void compute_recurrence_matrix(double **D, int n, int **R, double threshold) 
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            R[i][j] = (D[i][j] < threshold) ? 1 : 0;
        }
    }
}

double compute_recurrence_rate(int **R, int n) 
{
    int count_recurrent_points = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j && R[i][j] == 1) {
                count_recurrent_points++;
            }
        }
    }

    return (double) count_recurrent_points / (double) (n * n);
}

double adjust_threshold_via_recurrence_rate(double **D, int n, int **R, double desired_percentage, int *status, double *actual_percentage) 
{
    if (n <= 0) {
        *status = -1;            // Invalid size
        if (actual_percentage)  *actual_percentage = -1.0;
        return -1.0;
    }

    // Find the maximum value in the upper triangle of D to set high.
    double low = 0.0, high = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (D[i][j] > high) {
                high = D[i][j];
            }
        }
    }

    double mid, recurrence_rate;
    double best_mid = 0.0;
    double best_recurrence_rate = 0.0;      // Will store the rate for best_mid
    double closest_diff = DBL_MAX;

    *status = 1;  // Default: no exact match, but we will return the closest.

    // Binary‐search loop until high‐low is very small.
    while ((high - low) > 1e-10) {
        mid = 0.5 * (low + high);

        compute_recurrence_matrix(D, n, R, mid);
        recurrence_rate = compute_recurrence_rate(R, n) * 100.0;  // in percent

        double diff = fabs(recurrence_rate - desired_percentage);
        if (diff < closest_diff) {
            closest_diff = diff;
            best_mid = mid;
            best_recurrence_rate = recurrence_rate;

            // If we’re within 1e-5 percentage points, consider it “exact.”
            if (diff < 1e-5) {
                *status = 0;
                if (actual_percentage) {
                    *actual_percentage = recurrence_rate;
                }
                return mid;
            }
        }

        if (recurrence_rate < desired_percentage) {
            low = mid;
        } else {
            high = mid;
        }
    }

    // If we exit the loop without hitting the 1e-5 tolerance,
    // best_mid holds the threshold whose recurrence_rate is closest.
    if (actual_percentage) {
        *actual_percentage = best_recurrence_rate;
    }
    return best_mid;
}

int* calculate_diagonal_frequency_distribution(int** R, int n) 
{
    int i, j, k, diagonal_line_length;
    int *diagonal_frequency_distribution = (int*)malloc((n + 1) * sizeof(int));
    
    if (diagonal_frequency_distribution == NULL) {
        // Handle memory allocation failure
        return NULL;
    }

    for(i = 0; i < n + 1; i++)
        diagonal_frequency_distribution[i] = 0;

    for(i = n - 1; i > -1; i--) {
        diagonal_line_length = 0;
        for (j = 0; j < n - i; j++) {
            if (R[i + j][j] == 1) {
                diagonal_line_length += 1;
                if (j == (n - i - 1))
                    diagonal_frequency_distribution[diagonal_line_length] += 1;
            } else {
                if (diagonal_line_length != 0) {
                    diagonal_frequency_distribution[diagonal_line_length] += 1;
                    diagonal_line_length = 0;
                }
            }
        }
    }

    for (k = 1; k < n; k++) {
        diagonal_line_length = 0;
        for (i = 0; i < n - k; i++) {
            j = i + k;
            if (R[i][j] == 1) {
                diagonal_line_length += 1;
                if (j == (n - 1))
                    diagonal_frequency_distribution[diagonal_line_length] += 1;
            } else {
                if (diagonal_line_length != 0) {
                    diagonal_frequency_distribution[diagonal_line_length] += 1;
                    diagonal_line_length = 0;
                }
            }
        }
    }

    return diagonal_frequency_distribution;
}

int* calculate_vertical_frequency_distribution(int** R, int n)
{
    int i, j;
    int *vertical_frequency_distribution = (int*)malloc((n + 1) * sizeof(int));
    
    if (vertical_frequency_distribution == NULL) {
        // Handle memory allocation failure
        return NULL;
    }

    for (i = 0; i < n + 1; i++)
        vertical_frequency_distribution[i] = 0;

    for (j = 0; j < n; j++) {
        int vertical_line_length = 0;
        for (i = 0; i < n; i++) {
            if (R[i][j] == 1) {
                vertical_line_length += 1;
                if (i == (n - 1))
                    vertical_frequency_distribution[vertical_line_length] += 1;
            } else {
                if (vertical_line_length != 0) {
                    vertical_frequency_distribution[vertical_line_length] += 1;
                    vertical_line_length = 0;
                }
            }
        }
    }

    return vertical_frequency_distribution;
}

double compute_average_diagonal_length(int** R, int n) 
{
    double numerator = 0.0;
    double denominator = 0.0;
    int l_min = 2;

    // Calculate diagonal frequency distribution using the helper function
    int *diagonal_frequency_distribution = calculate_diagonal_frequency_distribution(R, n);
    if (diagonal_frequency_distribution == NULL) {
        // Handle error if allocation failed
        return 0.0;
    }

    // Calculate average diagonal length
    for (int l = l_min; l < n; l++) {
        numerator += l * diagonal_frequency_distribution[l];
        denominator += diagonal_frequency_distribution[l];
    }

    // Free the dynamically allocated memory
    free(diagonal_frequency_distribution);

    // Avoid division by zero
    return (denominator == 0) ? 0.0 : (double)numerator / denominator;
}

double compute_longest_diagonal_length(int **R, int n) 
{
    int *diagonal_frequency_distribution = calculate_diagonal_frequency_distribution(R, n);
    if (diagonal_frequency_distribution == NULL) {
        return 0;
    }

    int longest_diagonal_line_length = 0; // Initialize with 0
    for (int l = n - 1; l > 0; l--) {
        if (diagonal_frequency_distribution[l] != 0) {
            longest_diagonal_line_length = l;
            break;
        }
    }

    free(diagonal_frequency_distribution);
    
    return longest_diagonal_line_length;
}

double compute_average_vertical_length(int** R, int n) 
{
    double numerator   = 0.0;
    double denominator = 0.0;
    int v_min = 2;

    int *vertical_frequency_distribution = calculate_vertical_frequency_distribution(R, n);
    if (vertical_frequency_distribution == NULL) {
        return 0.0;
    }

    for (int v = v_min; v < n + 1; v++) {
        numerator += v * vertical_frequency_distribution[v];
        denominator += vertical_frequency_distribution[v];
    }

    // Free the dynamically allocated memory
    free(vertical_frequency_distribution);

    // Avoid division by zero
    return (denominator == 0) ? 0.0 : (double)numerator / denominator;
}

double compute_longest_vertical_length(int **R, int n) 
{
    int *vertical_frequency_distribution = calculate_vertical_frequency_distribution(R, n);
    if (vertical_frequency_distribution == NULL) {
        // Handle error if allocation failed
        return 0.0;
    }

    int longest_vertical_line_length;
    for (int v = n; v > 0; v--) {
        if (vertical_frequency_distribution[v] != 0) {
            longest_vertical_line_length = v;
            break;
        }
    }

    free(vertical_frequency_distribution);
    
    return longest_vertical_line_length;
}

double compute_determinism(int** R, int n) 
{
    int *diagonal_frequency_distribution = calculate_diagonal_frequency_distribution(R, n);
    int l_min = 2;

    if (diagonal_frequency_distribution == NULL) {
        // Handle error if allocation failed
        fprintf(stderr, "Memory allocation for diagonal_frequency_distribution failed.\n");
        return -1.0; // Return an error value
    }

    double numerator = 0.0;
    double denominator = 0.0;

    for (int l = l_min; l < n; l++) {
        numerator += l * diagonal_frequency_distribution[l];
    }
    
    for (int l = 1; l < n; l++) {
        denominator += l * diagonal_frequency_distribution[l];
    }

    // Free the dynamically allocated memory
    free(diagonal_frequency_distribution);

    // Avoid division by zero
    if (denominator == 0) return 0.0;

    double determinism = numerator / denominator;

    // Ensure determinism is bounded by 1
    if (determinism > 1.0) determinism = 1.0;

    return determinism;
}

double compute_laminarity(int **R, int n) 
{
    int *vertical_frequency_distribution = calculate_vertical_frequency_distribution(R, n);
    int v_min = 2;

    if (vertical_frequency_distribution == NULL) {
        // Handle error if allocation failed
        fprintf(stderr, "Memory allocation for vertical_frequency_distribution failed.\n");
        return -1.0; // Return an error value
    }

    double numerator = 0.0;
    double denominator = 0.0;

    for (int v = v_min; v < n + 1; v++) {
        numerator += v * vertical_frequency_distribution[v];
    }

    for (int v = 1; v < n + 1; v++) {
        denominator += v * vertical_frequency_distribution[v];
    }

    free(vertical_frequency_distribution);

    // Avoid division by zero
    if (denominator == 0) return 0.0;

    double laminarity = numerator / denominator;

    // Ensure laminarity is bounded by 1
    if (laminarity > 1.0) laminarity = 1.0;

    return laminarity;
}

double compute_divergence(int **R, int n) 
{
    return 1.0 / compute_longest_diagonal_length(R, n);
}

double compute_trapping_time(int **R, int n) 
{
    int *vertical_frequency_distribution = calculate_vertical_frequency_distribution(R, n);
    int v_min = 2;

    if (vertical_frequency_distribution == NULL) {
        // Handle error if allocation failed
        return 0.0;
    }

    double numerator = 0.0;
    double denominator = 0.0;

    for (int v = v_min; v <= n; v++) {
        numerator += v * vertical_frequency_distribution[v];
        denominator += vertical_frequency_distribution[v];
    }

    free(vertical_frequency_distribution);

    return (denominator == 0) ? 0.0 : (double)numerator / denominator;
}

double compute_entropy_diagonal(int **R, int n) 
{
    int *diagonal_frequency_distribution = calculate_diagonal_frequency_distribution(R, n);
    int l_min = 2;

    if (diagonal_frequency_distribution == NULL) {
        // Handle error if allocation failed
        return 0.0;
    }

    double sum_diagonal_frequency_distribution = 0.0;
    for(int l = l_min; l < n; l++)
    {
        sum_diagonal_frequency_distribution += diagonal_frequency_distribution[l];
    }
        
    double entropy_diagonal_lines = 0.0;

    if (sum_diagonal_frequency_distribution > 0.0)
    {
        for (int l = l_min; l < n; l++)
        {
            if (diagonal_frequency_distribution[l] != 0)
            {
                entropy_diagonal_lines += (diagonal_frequency_distribution[l]/sum_diagonal_frequency_distribution) * log(diagonal_frequency_distribution[l]/sum_diagonal_frequency_distribution);
            }
                
        }
        
        entropy_diagonal_lines *= -1.0;    
    }

    free(diagonal_frequency_distribution);

    return entropy_diagonal_lines;
}

double compute_entropy_vertical(int **R, int n)
{
    int *vertical_frequency_distribution = calculate_vertical_frequency_distribution(R, n);
    int v_min = 2;

    if (vertical_frequency_distribution == NULL) {
        // Handle error if allocation failed
        return 0.0;
    }

    double sum_vertical_frequency_distribution = 0.0;
    
    for (int v = v_min; v < n + 1; v++)
    {
        sum_vertical_frequency_distribution += vertical_frequency_distribution[v];
    }

    double entropy_vertical_lines = 0.0;

    if(sum_vertical_frequency_distribution > 0.0)
    {
        for (int v = v_min; v < n + 1; v++)
        {
            if (vertical_frequency_distribution[v] != 0)
            {
                entropy_vertical_lines += (vertical_frequency_distribution[v]/sum_vertical_frequency_distribution) * log(vertical_frequency_distribution[v]/sum_vertical_frequency_distribution);
            }
                
        }
        entropy_vertical_lines *= -1.0;    
    }

    free(vertical_frequency_distribution);

    return entropy_vertical_lines;
}

double compute_new_entropy_diagonal(int **R, int n) 
{
    // Expect calculate_diagonal_frequency_distribution to return an array of size (n+1)
    int *P = calculate_diagonal_frequency_distribution(R, n);
    const int l_min = 2;

    if (!P) {
        // Allocation or calculation failed
        return 0.0;
    }

    // Compute total count for l >= l_min
    double total = 0.0;
    for (int l = l_min; l <= n; ++l) {
        total += (double)P[l];
    }

    double H = 0.0;
    if (total > 0.0) {
        for (int l = l_min; l <= n; ++l) {
            int count = P[l];
            if (count > 0) {
                double p = (double)count / total;
                H -= p * log(p);
            }
        }
    }

    free(P);
    return H;
}

double compute_new_entropy_vertical(int **R, int n)
{
    // Expect calculate_vertical_frequency_distribution to return an array of size (n+1)
    int *P = calculate_vertical_frequency_distribution(R, n);
    const int v_min = 2;

    if (!P) {
        // Allocation or calculation failed
        return 0.0;
    }

    // Compute total count for v >= v_min
    double total = 0.0;
    for (int v = v_min; v <= n; ++v) {
        total += (double)P[v];
    }

    double H = 0.0;
    if (total > 0.0) {
        for (int v = v_min; v <= n; ++v) {
            int count = P[v];
            if (count > 0) {
                double p = (double)count / total;
                H -= p * log(p);
            }
        }
    }

    free(P);
    return H;
}

double wrap_longest_diagonal_length(int **matrix, int size) 
{
    return (double)compute_longest_diagonal_length(matrix, size);
}

double wrap_longest_vertical_length(int **matrix, int size) 
{
    return (double)compute_longest_vertical_length(matrix, size);
}

int compare_lines_desc(const void *pA, const void *pB) 
{
    const Line *a = (const Line *)pA;
    const Line *b = (const Line *)pB;
    return (b->length - a->length);
}

void delete_line_recursive(int **CR, int N, int M, int r, int c) 
{
    if (r < 0 || r >= N || c < 0 || c >= M) return;
    if (CR[r][c] == 0) return;

    // 1) Find run_start going upward until CR[run_start-1][c] == 0 or run_start == 0
    int run_start = r;
    while (run_start > 0 && CR[run_start - 1][c] == 1) {
        run_start--;
    }

    // 2) Measure run length len
    int len = 0;
    while (run_start + len < N && CR[run_start + len][c] == 1) {
        len++;
    }

    // 3) Zero out those len pixels in column c
    for (int k = 0; k < len; k++) {
        CR[run_start + k][c] = 0;
    }

    // 4) Recurse on neighbors (left/right) for each pixel of that run
    for (int k = 0; k < len; k++) {
        int rr = run_start + k;
        if (c - 1 >= 0)
            delete_line_recursive(CR, N, M, rr, c - 1);
        if (c + 1 < M)
            delete_line_recursive(CR, N, M, rr, c + 1);
    }
}

void build_close_returns_map(int **RM, int N, int **CR) 
{
    int M = N + 1;
    // Zero‐out CR (caller may already have done this)
    for (int i = 0; i < N; i++) {
        for (int c = 0; c < M; c++) {
            CR[i][c] = 0;
        }
    }

    // Fill only lower‐triangle entries from RM
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            if (RM[i][j] == 1) {
                int c = N + (j - i);
                CR[i][c] = 1;
            }
        }
    }
}

void build_line_list(int **CR, int N, int M, Line **out_lines, int *out_count) 
{
    int capacity = 256;
    int count = 0;
    Line *lines = (Line *)malloc(capacity * sizeof(Line));
    if (!lines) {
        fprintf(stderr, "Error: malloc failed in build_line_list.\n");
        exit(EXIT_FAILURE);
    }

    for (int c = 0; c < M; c++) {
        int i = 0;
        while (i < N) {
            if (CR[i][c] == 1) {
                int start_i = i;
                int len = 0;
                while (i < N && CR[i][c] == 1) {
                    len++;
                    i++;
                }
                if (count >= capacity) {
                    capacity *= 2;
                    lines = (Line *)realloc(lines, capacity * sizeof(Line));
                    if (!lines) {
                        fprintf(stderr, "Error: realloc failed in build_line_list.\n");
                        exit(EXIT_FAILURE);
                    }
                }
                lines[count].length    = len;
                lines[count].start_row = start_i;
                lines[count].col       = c;
                count++;
            } else {
                i++;
            }
        }
    }

    *out_lines = lines;
    *out_count = count;
}

void get_final_line_list(Line *lines, int linecount, int **CR_copy, int N, int M, Line **out_final, int *out_final_count)
{
    // 1) Sort by descending length
    qsort(lines, linecount, sizeof(Line), compare_lines_desc);

    int capacity = 256;
    int final_count = 0;
    Line *final = (Line *)malloc(capacity * sizeof(Line));
    if (!final) {
        fprintf(stderr, "Error: malloc failed in get_final_line_list.\n");
        exit(EXIT_FAILURE);
    }

    // 2) Iterate through sorted runs
    for (int idx = 0; idx < linecount; idx++) {
        int len = lines[idx].length;
        int r0  = lines[idx].start_row;
        int c   = lines[idx].col;

        // Check if this entire run is still present:
        int still_present = 1;
        for (int k = 0; k < len; k++) {
            int rr = r0 + k;
            if (rr < 0 || rr >= N || CR_copy[rr][c] == 0) {
                still_present = 0;
                break;
            }
        }
        if (!still_present) continue;

        // Accept the run:
        if (final_count >= capacity) {
            capacity *= 2;
            final = (Line *)realloc(final, capacity * sizeof(Line));
            if (!final) {
                fprintf(stderr, "Error: realloc failed in get_final_line_list.\n");
                exit(EXIT_FAILURE);
            }
        }
        final[final_count++] = lines[idx];

        // Zero out its pixels in CR_copy
        for (int k = 0; k < len; k++) {
            CR_copy[r0 + k][c] = 0;
        }

        // Recursively delete adjacent runs that touch any of its pixels:
        for (int k = 0; k < len; k++) {
            int rr = r0 + k;
            if (c - 1 >= 0) delete_line_recursive(CR_copy, N, M, rr, c - 1);
            if (c + 1 < M)  delete_line_recursive(CR_copy, N, M, rr, c + 1);
        }
    }

    *out_final       = final;
    *out_final_count = final_count;
}

void build_skeletonized_RP_from_runs(Line *final, int final_count, int N, int **SK) 
{
    // Zero‐out SK (caller may have already done this)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            SK[i][j] = 0;
        }
    }

    for (int idx = 0; idx < final_count; idx++) {
        int len = final[idx].length;
        int r0  = final[idx].start_row;
        int c   = final[idx].col;
        for (int k = 0; k < len; k++) {
            int rr = r0 + k;
            int cc = rr + (c - N);
            if (rr < 0 || rr >= N || cc < 0 || cc >= N) continue;
            SK[rr][cc] = 1;
            SK[cc][rr] = 1;
        }
    }

    // Ensure line‐of‐identity
    for (int i = 0; i < N; i++) {
        SK[i][i] = 1;
    }
}

void skeletonize_recurrence_matrix(int **RM, int N, int **SK)
{
    int M = N + 1;

    // 1) Allocate and build CR[N][M]
    int **CR = NULL;
    alloc_2d_int(&CR, N, M);
    build_close_returns_map(RM, N, CR);

    // 2) Allocate CR_copy[N][M] and copy CR into it
    int **CR_copy = NULL;
    alloc_2d_int(&CR_copy, N, M);
    for (int i = 0; i < N; i++) {
        memcpy(CR_copy[i], CR[i], M * sizeof(int));
    }

    // 3) Extract all vertical runs into “lines”
    Line *lines = NULL;
    int linecount = 0;
    build_line_list(CR, N, M, &lines, &linecount);

    // 4) Select final runs (deleting overlaps) → “final_lines”
    Line *final_lines = NULL;
    int final_count = 0;
    get_final_line_list(lines, linecount, CR_copy, N, M, &final_lines, &final_count);

    // 5) Rebuild SK[N][N] from final_lines
    build_skeletonized_RP_from_runs(final_lines, final_count, N, SK);

    // 6) Free temporaries
    free(lines);
    free(final_lines);
    dealloc_2d_int(&CR, N);
    dealloc_2d_int(&CR_copy, N);
}