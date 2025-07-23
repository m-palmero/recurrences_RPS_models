#include "lib/common.h"
#include "RQA_in_RPS_systems.h"
#include "lib/recurrence_functions.h"
#include "lib/auxiliary_functions.h"
#include "lib/gnuplot_functions.h"

//================================================================================
// Global Configuration
//================================================================================

// Define base paths for output files
const char *results_path = "results/";
const char *figs_path = "figs/";

// Global buffers for constructing filenames
char data[100]; 
char filename1[500], filename2[500], filename3[500], filename4[500], filename5[500]; 
char filename6[500], filename7[500], filename8[500], filename9[500], filename10[500];
char filename11[500];

// --- RQA & Time Series Parameters ---
int N = 2000;                   // Length of the time series to be analyzed
int emb_dim = 2;                // Embedding dimension for phase space reconstruction
int tau = 1;                      // Time delay for phase space reconstruction
int species = 1;                // Data column to analyze (0=empty, 1=1st species, 2=2nd, 3=3rd)
int size;                       // (Currently unused) Placeholder for varying size analysis
double eps = 0.001;             // Fixed recurrence threshold (used if define_threshold_via_percentage is false)
int max_length = 2000;          // (Currently unused)
int number_sub_series = 100;    // (Currently unused)
double mobility;                // Mobility parameter (set within analysis functions)
double percentage = 10;         // Target recurrence rate (RR) in percent (used if define_threshold_via_percentage is true)

// --- Analysis Switches ---
bool define_threshold_via_percentage = true; // If true, RR is fixed; if false, eps is fixed.
bool enable_skeletonization = false;         // If true, RQA is performed on the skeletonized recurrence matrix.


//================================================================================
// Helper Functions
//================================================================================

/**
 * @brief Reads a time series from a specific data file.
 * * @param time_series Pointer to the array where the time series will be stored.
 * @param mobility The mobility parameter, used to construct the filename.
 * @param initial The simulation run ID (initial condition), used for the filename.
 * @return int 0 on success, 1 on failure.
 */
int read_file_for_time_series(double *time_series, double mobility, int initial)
{
    // Construct filename based on mobility and run ID
    sprintf(data, "data/mobility/%1.2f-%d.dat", mobility, initial);
    FILE *f = fopen(data, "r");

    int number_of_columns = 4;
    int target_column = species; // 0 = 1st column; 1 = 2nd column; ...
    
    if (f == NULL) 
    {
        fprintf(stderr, "Error opening file %s\n", data);
        return 1;
    }

    // Read the file line by line
    for (int i = 0; i < N; ++i) 
    {
        double values[number_of_columns];
        
        // Read all columns from the current line
        for (int j = 0; j < number_of_columns; ++j) 
        {
            if (fscanf(f, "%lf", &values[j]) != 1) 
            {
                fprintf(stderr, "Error reading values from file\n");
                fclose(f);
                return 1;
            }
        }

        // Store only the value from the target column
        time_series[i] = values[target_column];
    }
    fclose(f);

    return 0;
}

/**
 * @brief Writes a time series to a file named "ts.dat". Used for debugging.
 * * @param time_series The time series data to write.
 * @return int 0 on success.
 */
int give_time_series(double *time_series)
{
    FILE *x = fopen("ts.dat", "w");
    for (int i = 0; i < N; ++i) 
    {
        fprintf(x,"%f\n", time_series[i]);
    }
    fclose(x);
    
    return 0;
}

/**
 * @brief Generates output files for a recurrence plot.
 * It creates "rm.dat" (the full recurrence matrix) and "rp.dat" (the coordinates of recurrent points).
 * * @param R The recurrence matrix.
 * @param n The size of the matrix.
 * @return int 0 on success, 1 on failure.
 */
int give_recurrence_plot(int **R, int n)
{
    int x_rp[n], y_rp[n];

    // File for the full matrix
    sprintf(filename1, "rm.dat");
    FILE *rm = fopen(filename1, "w");
    // File for recurrent point coordinates (for plotting)
    sprintf(filename2, "rp.dat");
    FILE *rp = fopen(filename2, "w");

    if (rm == NULL || rp == NULL) 
    {
        fprintf(stderr, "Error opening files for writing.\n");
        return 1;
    }

    // Iterate through the matrix
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fprintf(rm, "%d ", R[i][j]);
            if (R[i][j] == 1)
            {
                // If it's a recurrent point, save its coordinates
                x_rp[i]=i;
                y_rp[j]=j;
                fprintf(rp,"%d %d\n", x_rp[i], y_rp[j]);
            }
        }
        fprintf(rm, "\n");
    }

    fclose(rm);
    fclose(rp);
    return 0;
}
    
//================================================================================
// Analysis Functions
//================================================================================

/**
 * @brief Performs RQA on an ensemble of simulations, treating each simulation run ('initial' ID) as an independent sample.
 * NOTE: As currently written, this function uses the single global 'mobility' variable to read all data files.
 * It is designed to analyze the statistical distribution of RQA measures for a set of simulations
 * that were all generated with the same mobility parameter.
 * The output is a separate file for each RQA measure, where each row corresponds to a simulation ID.
 */
int ensemble_recurrence_analysis()
{
    int size = 0;
    int ensemble_size = 10000;
    // These flags seem to be for features that are not fully implemented or are deprecated.
    bool species_analysis = false;
    bool varying_sizes = false;

    // Allocate arrays to store RQA measures for each run in the ensemble
    double *e_RR, *e_DET, *e_LAM, *e_TT, *e_ENTR_L, *e_ENTR_V, *e_L_MED, *e_V_MED, *e_DIV, *e_L_MAX, *e_V_MAX;
    alloc_1d_double(&e_RR, ensemble_size);
    alloc_1d_double(&e_DET, ensemble_size);
    alloc_1d_double(&e_LAM, ensemble_size);
    alloc_1d_double(&e_TT, ensemble_size);
    alloc_1d_double(&e_ENTR_L, ensemble_size);
    alloc_1d_double(&e_ENTR_V, ensemble_size);
    alloc_1d_double(&e_L_MED, ensemble_size);
    alloc_1d_double(&e_V_MED, ensemble_size);
    alloc_1d_double(&e_DIV, ensemble_size);
    alloc_1d_double(&e_L_MAX, ensemble_size);
    alloc_1d_double(&e_V_MAX, ensemble_size);

    // Species-specific filename generation
    if (species_analysis)
    {
        switch (species)
        {
        case 1:
            if (varying_sizes)
            {
                sprintf(filename1, "1st_spec_RR_ensemble_%d.dat", size);
                sprintf(filename2, "1st_spec_DET_ensemble_%d.dat", size);
                sprintf(filename3, "1st_spec_LAM_ensemble_%d.dat", size);
                sprintf(filename4, "1st_spec_TT_ensemble_%d.dat", size);
                sprintf(filename5, "1st_spec_ENTR_L_ensemble_%d.dat", size);
            }
            else
            {
                sprintf(filename1, "1st_spec_RR_ensemble.dat");
                sprintf(filename2, "1st_spec_DET_ensemble.dat");
                sprintf(filename3, "1st_spec_LAM_ensemble.dat");
                sprintf(filename4, "1st_spec_TT_ensemble.dat");
                sprintf(filename5, "1st_spec_ENTR_L_ensemble.dat");
            }

            break;
        case 2:
            if (varying_sizes)
            { 
                sprintf(filename1, "2nd_spec_RR_ensemble_%d.dat", size);
                sprintf(filename2, "2nd_spec_DET_ensemble_%d.dat", size);
                sprintf(filename3, "2nd_spec_LAM_ensemble_%d.dat", size);
                sprintf(filename4, "2nd_spec_TT_ensemble_%d.dat", size);
                sprintf(filename5, "2nd_spec_ENTR_L_ensemble_%d.dat", size);
            }
            else
            {
                sprintf(filename1, "2nd_spec_RR_ensemble.dat");
                sprintf(filename2, "2nd_spec_DET_ensemble.dat");
                sprintf(filename3, "2nd_spec_LAM_ensemble.dat");
                sprintf(filename4, "2nd_spec_TT_ensemble.dat");
                sprintf(filename5, "2nd_spec_ENTR_L_ensemble.dat");
            }
            break;
        case 3:
            if (varying_sizes)
            { 
                sprintf(filename1, "3rd_spec_RR_ensemble_%d.dat", size);
                sprintf(filename2, "3rd_spec_DET_ensemble_%d.dat", size);
                sprintf(filename3, "3rd_spec_LAM_ensemble_%d.dat", size);
                sprintf(filename4, "3rd_spec_TT_ensemble_%d.dat", size);
                sprintf(filename5, "3rd_spec_ENTR_L_ensemble_%d.dat", size);
            }
            else
            {
                sprintf(filename1, "3rd_spec_RR_ensemble.dat");
                sprintf(filename2, "3rd_spec_DET_ensemble.dat");
                sprintf(filename3, "3rd_spec_LAM_ensemble.dat");
                sprintf(filename4, "3rd_spec_TT_ensemble.dat");
                sprintf(filename5, "3rd_spec_ENTR_L_ensemble.dat");
            }
            break;
        }
    }
    
    // Create the output directory based on the recurrence rate percentage
    size_t folder_name_size = strlen(results_path) + 256;
    char *folder_name = (char *)malloc(folder_name_size * sizeof(char));
    snprintf(folder_name, folder_name_size, "%sRR=%.2f/", results_path, percentage);
    create_folders_in_path_if_not_exists(folder_name);

    // Prepare filenames for each RQA measure
    sprintf(filename1, "%sRR_N=%d.dat", folder_name, N);
    sprintf(filename2, "%sDET_N=%d.dat", folder_name, N);
    sprintf(filename3, "%sLAM_N=%d.dat", folder_name, N);
    sprintf(filename4, "%sTT_N=%d.dat", folder_name, N);
    sprintf(filename5, "%sENTR_L_N=%d.dat", folder_name, N);
    sprintf(filename6, "%sENTR_V_N=%d.dat", folder_name, N);
    sprintf(filename7, "%sL_MED_N=%d.dat", folder_name, N);
    sprintf(filename8, "%sV_MED_N=%d.dat", folder_name, N);
    sprintf(filename9, "%sDIV_N=%d.dat", folder_name, N);
    sprintf(filename10, "%sL_MAX_N=%d.dat", folder_name, N);
    sprintf(filename11, "%sV_MAX_N=%d.dat", folder_name, N);
    
    // Open all output files for writing
    FILE *out1 = fopen(filename1, "w");
    FILE *out2 = fopen(filename2, "w");
    FILE *out3 = fopen(filename3, "w");
    FILE *out4 = fopen(filename4, "w");
    FILE *out5 = fopen(filename5, "w");
    FILE *out6 = fopen(filename6, "w");
    FILE *out7 = fopen(filename7, "w");
    FILE *out8 = fopen(filename8, "w");
    FILE *out9 = fopen(filename9, "w");
    FILE *out10 = fopen(filename10, "w");
    FILE *out11 = fopen(filename11, "w");

    // Loop through each simulation file in the ensemble
    for (int initial = 1; initial <= ensemble_size; initial++)
    {
        print_prog((double) (initial+1) / (double) ensemble_size);

        // --- RQA Pipeline for a single time series ---
        double *time_series;
        alloc_1d_double(&time_series, N);
        read_file_for_time_series(time_series, mobility, initial);
        give_time_series(time_series); // For debugging

        double **embedded_time_series;
        alloc_2d_double(&embedded_time_series, N, emb_dim);
        embed_time_series(embedded_time_series, time_series, N, emb_dim, tau);

        double **DM;
        alloc_2d_double(&DM, N, N);
        calculate_distance_matrix(embedded_time_series, N, emb_dim, DM, "euclidean");

        int **RM;
        alloc_2d_int(&RM, N, N);
        
        int status = 0;
        double actual_pct = 0.0;
        double threshold;

        if (define_threshold_via_percentage) {
            threshold = adjust_threshold_via_recurrence_rate(DM, N, RM, percentage, &status, &actual_pct);
        } else {
            threshold = eps;
            compute_recurrence_matrix(DM, N, RM, threshold);
        }

        // --- Compute all RQA measures ---
        e_RR[initial] = compute_recurrence_rate(RM, N);
        e_DET[initial] = compute_determinism(RM, N);
        e_LAM[initial] = compute_laminarity(RM, N);
        e_TT[initial] = compute_trapping_time(RM, N);
        e_ENTR_L[initial] = compute_entropy_diagonal(RM, N);
        e_ENTR_V[initial] = compute_entropy_vertical(RM, N);
        e_L_MED[initial] = compute_average_diagonal_length(RM, N);
        e_V_MED[initial] = compute_average_vertical_length(RM, N);
        e_DIV[initial] = compute_divergence(RM, N);
        e_L_MAX[initial] = compute_longest_diagonal_length(RM, N);
        e_V_MAX[initial] = compute_longest_vertical_length(RM, N);
        
        // Write the calculated measure to its corresponding file
        fprintf(out1,"%d %f\n", initial, e_RR[initial]);
        fprintf(out2,"%d %f\n", initial, e_DET[initial]);
        fprintf(out3,"%d %f\n", initial, e_LAM[initial]);
        fprintf(out4,"%d %f\n", initial, e_TT[initial]);
        fprintf(out5,"%d %f\n", initial, e_ENTR_L[initial]);
        fprintf(out6,"%d %f\n", initial, e_ENTR_V[initial]);
        fprintf(out7,"%d %f\n", initial, e_L_MED[initial]);
        fprintf(out8,"%d %f\n", initial, e_V_MED[initial]);
        fprintf(out9,"%d %f\n", initial, e_DIV[initial]);
        fprintf(out10,"%d %f\n", initial, e_L_MAX[initial]);
        fprintf(out11,"%d %f\n", initial, e_V_MAX[initial]);
        
        // Free memory for the current time series
        dealloc_2d_double(&DM, N);
        dealloc_2d_int(&RM, N);
        dealloc_2d_double(&embedded_time_series, N);
        dealloc_1d_double(&time_series);
    }

    // Close all output files
    fclose(out1);
    fclose(out2);
    fclose(out3);
    fclose(out4);
    fclose(out5);
    fclose(out6);
    fclose(out7);
    fclose(out8);
    fclose(out9);
    fclose(out10);
    fclose(out11);

    // Free the memory allocated for the ensemble statistics
    dealloc_1d_double(&e_RR);
    dealloc_1d_double(&e_DET);
    dealloc_1d_double(&e_LAM);
    dealloc_1d_double(&e_TT);
    dealloc_1d_double(&e_ENTR_L);
    dealloc_1d_double(&e_ENTR_V);
    dealloc_1d_double(&e_L_MED);
    dealloc_1d_double(&e_V_MED);
    dealloc_1d_double(&e_DIV);
    dealloc_1d_double(&e_L_MAX);
    dealloc_1d_double(&e_V_MAX);

    return 0;
}

/**
 * @brief **PRIMARY ANALYSIS FUNCTION.**
 * Iterates over a range of mobility values. For each mobility, it performs RQA on an
 * entire ensemble of simulations, then calculates the mean and standard error of each
 * RQA measure. The final output shows how RQA measures change with mobility.
 */
int ensemble_analysis_varying_mobility()
{
    // --- Configuration for this analysis ---
    int ensemble_size = 25; // Number of simulations per mobility value
    double mobility_min  = 0.0;
    double mobility_max  = 0.9;
    double mobility_step = 0.03;
    int num_steps = (int)((mobility_max - mobility_min) / mobility_step) + 1;

    // Arrays to hold the final results: mean and standard error for each mobility step
    double avg_RR[num_steps],    se_RR[num_steps];
    double avg_DET[num_steps],   se_DET[num_steps];
    double avg_LAM[num_steps],   se_LAM[num_steps];
    double avg_TT[num_steps],    se_TT[num_steps];
    double avg_ENTR_L[num_steps], se_ENTR_L[num_steps];
    double avg_ENTR_V[num_steps], se_ENTR_V[num_steps];
    double avg_L_MED[num_steps], se_L_MED[num_steps];
    double avg_V_MED[num_steps], se_V_MED[num_steps];
    double avg_DIV[num_steps],   se_DIV[num_steps];
    double avg_L_MAX[num_steps], se_L_MAX[num_steps];
    double avg_V_MAX[num_steps], se_V_MAX[num_steps];

    // Initialize all statistics arrays to zero
    memset(avg_RR,     0, sizeof(avg_RR));     memset(se_RR,      0, sizeof(se_RR));
    memset(avg_DET,    0, sizeof(avg_DET));    memset(se_DET,     0, sizeof(se_DET));
    memset(avg_LAM,    0, sizeof(avg_LAM));    memset(se_LAM,     0, sizeof(se_LAM));
    memset(avg_TT,     0, sizeof(avg_TT));     memset(se_TT,      0, sizeof(se_TT));
    memset(avg_ENTR_L, 0, sizeof(avg_ENTR_L)); memset(se_ENTR_L,  0, sizeof(se_ENTR_L));
    memset(avg_ENTR_V, 0, sizeof(avg_ENTR_V)); memset(se_ENTR_V,  0, sizeof(se_ENTR_V));
    memset(avg_L_MED,  0, sizeof(avg_L_MED));  memset(se_L_MED,   0, sizeof(se_L_MED));
    memset(avg_V_MED,  0, sizeof(avg_V_MED));  memset(se_V_MED,   0, sizeof(se_V_MED));
    memset(avg_DIV,    0, sizeof(avg_DIV));    memset(se_DIV,     0, sizeof(se_DIV));
    memset(avg_L_MAX,  0, sizeof(avg_L_MAX));  memset(se_L_MAX,   0, sizeof(se_L_MAX));
    memset(avg_V_MAX,  0, sizeof(avg_V_MAX));  memset(se_V_MAX,   0, sizeof(se_V_MAX));

    // --- Main Loop: Iterate over each mobility value ---
    int m = 0; // Index for the current mobility step
    for (double mobility = mobility_min; mobility < mobility_max + 1e-6; mobility += mobility_step)
    {
        // Buffers to store RQA values for the current ensemble (at a fixed mobility)
        double *e_RR, *e_DET, *e_LAM, *e_TT;
        double *e_ENTR_L, *e_ENTR_V, *e_L_MED, *e_V_MED;
        double *e_DIV, *e_L_MAX, *e_V_MAX;
        alloc_1d_double(&e_RR,     ensemble_size); alloc_1d_double(&e_DET,    ensemble_size);
        alloc_1d_double(&e_LAM,    ensemble_size); alloc_1d_double(&e_TT,     ensemble_size);
        alloc_1d_double(&e_ENTR_L, ensemble_size); alloc_1d_double(&e_ENTR_V, ensemble_size);
        alloc_1d_double(&e_L_MED,  ensemble_size); alloc_1d_double(&e_V_MED,  ensemble_size);
        alloc_1d_double(&e_DIV,    ensemble_size); alloc_1d_double(&e_L_MAX,  ensemble_size);
        alloc_1d_double(&e_V_MAX,  ensemble_size);

        printf("Calculating RQA for mobility = %1.2f\n", mobility);

        // --- Inner Loop: Iterate over the ensemble for the current mobility ---
        for (int initial = 0; initial < ensemble_size; initial++)
        {
            print_prog((double)initial / (double)ensemble_size);

            // 1) Read and embed the time series
            double *time_series;
            alloc_1d_double(&time_series, N);
            read_file_for_time_series(time_series, mobility, initial);

            double **embedded_time_series;
            alloc_2d_double(&embedded_time_series, N, emb_dim);
            embed_time_series(embedded_time_series, time_series, N, emb_dim, tau);

            // 2) Build distance matrix
            double **DM;
            alloc_2d_double(&DM, N, N);
            calculate_distance_matrix(embedded_time_series, N, emb_dim, DM, "euclidean");

            // 3) Allocate the full recurrence matrix
            int **RM;
            alloc_2d_int(&RM, N, N);

            // 4) Determine threshold (either fixed RR or fixed epsilon)
            int status = 0;
            double actual_pct = 0.0;
            double threshold;
            if (define_threshold_via_percentage) {
                threshold = adjust_threshold_via_recurrence_rate(DM, N, RM, percentage, &status, &actual_pct);
            } else {
                threshold = eps;
            }

            // 5) Compute the full recurrence matrix using the determined threshold
            compute_recurrence_matrix(DM, N, RM, threshold);

            // 6) Optionally skeletonize the matrix before analysis
            int **RM_to_use = RM; // By default, use the original RM
            int **SK = NULL;
            if (enable_skeletonization) {
                alloc_2d_int(&SK, N, N);
                skeletonize_recurrence_matrix(RM, N, SK);
                RM_to_use = SK; // Use the skeletonized matrix for RQA
            }

            // 7) Compute all RQA measures on the chosen matrix (RM or SK)
            e_RR[initial]     = compute_recurrence_rate(RM_to_use, N);
            e_DET[initial]    = compute_determinism(RM_to_use, N);
            e_LAM[initial]    = compute_laminarity(RM_to_use, N);
            e_TT[initial]     = compute_trapping_time(RM_to_use, N);
            e_ENTR_L[initial] = compute_new_entropy_diagonal(RM_to_use, N);
            e_ENTR_V[initial] = compute_new_entropy_vertical(RM_to_use, N);
            e_L_MED[initial]  = compute_average_diagonal_length(RM_to_use, N);
            e_V_MED[initial]  = compute_average_vertical_length(RM_to_use, N);
            e_DIV[initial]    = compute_divergence(RM_to_use, N);
            e_L_MAX[initial]  = compute_longest_diagonal_length(RM_to_use, N);
            e_V_MAX[initial]  = compute_longest_vertical_length(RM_to_use, N);

            // 8) Clean up skeleton matrix if it was used
            if (SK != NULL) {
                dealloc_2d_int(&SK, N);
            }

            // 9) Accumulate sums to calculate the mean later
            avg_RR[m]     += e_RR[initial];     avg_DET[m]    += e_DET[initial];
            avg_LAM[m]    += e_LAM[initial];    avg_TT[m]     += e_TT[initial];
            avg_ENTR_L[m] += e_ENTR_L[initial]; avg_ENTR_V[m] += e_ENTR_V[initial];
            avg_L_MED[m]  += e_L_MED[initial];  avg_V_MED[m]  += e_V_MED[initial];
            avg_DIV[m]    += e_DIV[initial];    avg_L_MAX[m]  += e_L_MAX[initial];
            avg_V_MAX[m]  += e_V_MAX[initial];

            // 10) Free memory for the current trajectory
            dealloc_2d_double(&DM, N);
            dealloc_2d_int(&RM, N);
            dealloc_2d_double(&embedded_time_series, N);
            dealloc_1d_double(&time_series);
        }
        printf("\n");

        // 11) Compute ensemble means for the current mobility
        avg_RR[m]     /= ensemble_size;     avg_DET[m]    /= ensemble_size;
        avg_LAM[m]    /= ensemble_size;     avg_TT[m]     /= ensemble_size;
        avg_ENTR_L[m] /= ensemble_size;     avg_ENTR_V[m] /= ensemble_size;
        avg_L_MED[m]  /= ensemble_size;     avg_V_MED[m]  /= ensemble_size;
        avg_DIV[m]    /= ensemble_size;     avg_L_MAX[m]  /= ensemble_size;
        avg_V_MAX[m]  /= ensemble_size;

        // 12) Compute ensemble standard errors (sd / sqrt(n)) using GSL
        se_RR[m]     = gsl_stats_sd(e_RR,     1, ensemble_size) / sqrt(ensemble_size);
        se_DET[m]    = gsl_stats_sd(e_DET,    1, ensemble_size) / sqrt(ensemble_size);
        se_LAM[m]    = gsl_stats_sd(e_LAM,    1, ensemble_size) / sqrt(ensemble_size);
        se_TT[m]     = gsl_stats_sd(e_TT,     1, ensemble_size) / sqrt(ensemble_size);
        se_ENTR_L[m] = gsl_stats_sd(e_ENTR_L, 1, ensemble_size) / sqrt(ensemble_size);
        se_ENTR_V[m] = gsl_stats_sd(e_ENTR_V, 1, ensemble_size) / sqrt(ensemble_size);
        se_L_MED[m]  = gsl_stats_sd(e_L_MED,  1, ensemble_size) / sqrt(ensemble_size);
        se_V_MED[m]  = gsl_stats_sd(e_V_MED,  1, ensemble_size) / sqrt(ensemble_size);
        se_DIV[m]    = gsl_stats_sd(e_DIV,    1, ensemble_size) / sqrt(ensemble_size);
        se_L_MAX[m]  = gsl_stats_sd(e_L_MAX,  1, ensemble_size) / sqrt(ensemble_size);
        se_V_MAX[m]  = gsl_stats_sd(e_V_MAX,  1, ensemble_size) / sqrt(ensemble_size);

        // 13) Free per-ensemble buffers
        dealloc_1d_double(&e_RR);     dealloc_1d_double(&e_DET);
        dealloc_1d_double(&e_LAM);    dealloc_1d_double(&e_TT);
        dealloc_1d_double(&e_ENTR_L); dealloc_1d_double(&e_ENTR_V);
        dealloc_1d_double(&e_L_MED);  dealloc_1d_double(&e_V_MED);
        dealloc_1d_double(&e_DIV);    dealloc_1d_double(&e_L_MAX);
        dealloc_1d_double(&e_V_MAX);

        m++; // Move to the next mobility index
    }

    // 14) Prepare output folder and filenames based on thresholding method
    char folder_name[100];
    if (define_threshold_via_percentage) {
        // Create the specific sub-folder name (e.g., "results/RR=10.00")
        snprintf(folder_name, sizeof(folder_name), "%sRR=%.2f", results_path, percentage);
        create_folder_if_not_exists(folder_name);

        // Safely create the full path for each output file
        snprintf(filename1, sizeof(filename1), "%s/RR.dat",    folder_name);
        snprintf(filename2, sizeof(filename2), "%s/DET.dat",   folder_name);
        snprintf(filename3, sizeof(filename3), "%s/LAM.dat",   folder_name);
        snprintf(filename4, sizeof(filename4), "%s/TT.dat",    folder_name);
        snprintf(filename5, sizeof(filename5), "%s/ENTR_L.dat",folder_name);
        snprintf(filename6, sizeof(filename6), "%s/ENTR_V.dat",folder_name);
        snprintf(filename7, sizeof(filename7), "%s/L_MED.dat", folder_name);
        snprintf(filename8, sizeof(filename8), "%s/V_MED.dat", folder_name);
        snprintf(filename9, sizeof(filename9), "%s/DIV.dat",   folder_name);
        snprintf(filename10,sizeof(filename10),"%s/L_MAX.dat", folder_name);
        snprintf(filename11,sizeof(filename11),"%s/V_MAX.dat", folder_name);
    }
    else {
        // Create the specific sub-folder name (e.g., "results/eps=0.0010")
        snprintf(folder_name, sizeof(folder_name), "%seps=%.4f", results_path, eps);
        create_folder_if_not_exists(folder_name);

        // Safely create the full path for each output file
        snprintf(filename1, sizeof(filename1), "%s/RR.dat",    folder_name);
        snprintf(filename2, sizeof(filename2), "%s/DET.dat",   folder_name);
        snprintf(filename3, sizeof(filename3), "%s/LAM.dat",   folder_name);
        snprintf(filename4, sizeof(filename4), "%s/TT.dat",    folder_name);
        snprintf(filename5, sizeof(filename5), "%s/ENTR_L.dat",folder_name);
        snprintf(filename6, sizeof(filename6), "%s/ENTR_V.dat",folder_name);
        snprintf(filename7, sizeof(filename7), "%s/L_MED.dat", folder_name);
        snprintf(filename8, sizeof(filename8), "%s/V_MED.dat", folder_name);
        snprintf(filename9, sizeof(filename9), "%s/DIV.dat",   folder_name);
        snprintf(filename10,sizeof(filename10),"%s/L_MAX.dat", folder_name);
        snprintf(filename11,sizeof(filename11),"%s/V_MAX.dat", folder_name);
    }

    // 15) Open and write the final results to output files
    FILE *out1  = fopen(filename1,  "w"); FILE *out2  = fopen(filename2,  "w");
    FILE *out3  = fopen(filename3,  "w"); FILE *out4  = fopen(filename4,  "w");
    FILE *out5  = fopen(filename5,  "w"); FILE *out6  = fopen(filename6,  "w");
    FILE *out7  = fopen(filename7,  "w"); FILE *out8  = fopen(filename8,  "w");
    FILE *out9  = fopen(filename9,  "w"); FILE *out10 = fopen(filename10, "w");
    FILE *out11 = fopen(filename11, "w");

    if (!out1 || !out2 || !out3 || !out4 || !out5 || !out6 || !out7 || !out8 || !out9 || !out10 || !out11) {
        fprintf(stderr, "Error opening output files.\n");
        return -1;
    }

    for (int idx = 0; idx < num_steps; idx++) {
        double current_mobility = mobility_min + idx * mobility_step;
        // Write data in the format: mobility, mean, standard_error
        fprintf(out1,  "%1.2f %f %f\n", current_mobility, avg_RR[idx],     se_RR[idx]);
        fprintf(out2,  "%1.2f %f %f\n", current_mobility, avg_DET[idx],    se_DET[idx]);
        fprintf(out3,  "%1.2f %f %f\n", current_mobility, avg_LAM[idx],    se_LAM[idx]);
        fprintf(out4,  "%1.2f %f %f\n", current_mobility, avg_TT[idx],     se_TT[idx]);
        fprintf(out5,  "%1.2f %f %f\n", current_mobility, avg_ENTR_L[idx], se_ENTR_L[idx]);
        fprintf(out6,  "%1.2f %f %f\n", current_mobility, avg_ENTR_V[idx], se_ENTR_V[idx]);
        fprintf(out7,  "%1.2f %f %f\n", current_mobility, avg_L_MED[idx],  se_L_MED[idx]);
        fprintf(out8,  "%1.2f %f %f\n", current_mobility, avg_V_MED[idx],  se_V_MED[idx]);
        fprintf(out9,  "%1.2f %f %f\n", current_mobility, avg_DIV[idx],    se_DIV[idx]);
        fprintf(out10, "%1.2f %f %f\n", current_mobility, avg_L_MAX[idx],  se_L_MAX[idx]);
        fprintf(out11, "%1.2f %f %f\n", current_mobility, avg_V_MAX[idx],  se_V_MAX[idx]);
    }

    fclose(out1);  fclose(out2);  fclose(out3);
    fclose(out4);  fclose(out5);  fclose(out6);
    fclose(out7);  fclose(out8);  fclose(out9);
    fclose(out10); fclose(out11);

    return 0;
}

/**
 * @brief Performs a full RQA on a single time series.
 * This is useful for debugging, generating individual recurrence plots, and
 * getting a detailed look at a specific simulation run.
 */
int single_recurrence_analysis()
{
    // --- Configuration for this analysis ---
    int initial = 1;         // The run ID to analyze
    mobility = 0.60;         // The mobility to analyze
    bool single_RQA_values = false; // If true, print RQA values to console

    // --- RQA Pipeline ---
    double *time_series;
    alloc_1d_double(&time_series, N);
    read_file_for_time_series(time_series, mobility, initial);
    
    // Output the time series and generate a plot
    give_time_series(time_series);
    sprintf(filename3, "time_series.png");
    plot_gnuplot_time_series_no_title(filename3, (double)initial, N);

    double **embedded_time_series;
    alloc_2d_double(&embedded_time_series, N, emb_dim);
    embed_time_series(embedded_time_series, time_series, N, emb_dim, tau);

    double **DM;
    alloc_2d_double(&DM, N, N);
    calculate_distance_matrix(embedded_time_series, N, emb_dim, DM, "euclidean");

    int **RM;
    alloc_2d_int(&RM, N, N);
    
    int status = 0;
    double actual_pct = 0.0;
    double threshold;
    if (define_threshold_via_percentage) {
        threshold = adjust_threshold_via_recurrence_rate(DM, N, RM, percentage, &status, &actual_pct);
    } else {
        threshold = eps;
        compute_recurrence_matrix(DM, N, RM, threshold);
    }

    // Optionally skeletonize
    int **RM_to_use = RM;
    int **SK = NULL;
    if (enable_skeletonization) {
        alloc_2d_int(&SK, N, N);
        skeletonize_recurrence_matrix(RM, N, SK);
        RM_to_use = SK;
    }

    // Generate output files and plot for the recurrence plot
    give_recurrence_plot(RM_to_use, N);
    sprintf(filename4, "recurrence_plot.png");
    plot_gnuplot_RPs_no_title(filename4, N);

    // Optionally print all RQA values to the console
    if (single_RQA_values)
    {
        printf("RR = %f\n", compute_recurrence_rate(RM, N));
        printf("DET = %f\n", compute_determinism(RM, N));
        printf("LAM = %f\n", compute_laminarity(RM, N));
        printf("DIV = %f\n", compute_divergence(RM, N));
        printf("TT = %f\n", compute_trapping_time(RM, N));
        printf("ENTR_L = %f\n", compute_entropy_diagonal(RM, N));
        printf("ENTR_V = %f\n", compute_entropy_vertical(RM, N));
        printf("L_MED = %f\n", compute_average_diagonal_length(RM, N));
        printf("L_MAX = %f\n", compute_longest_diagonal_length(RM, N));
        printf("V_MED = %f\n", compute_average_vertical_length(RM, N));
        printf("V_MAX = %f\n", compute_longest_vertical_length(RM, N));
    }
    
    // --- Cleanup ---
    dealloc_2d_double(&DM, N);
    dealloc_2d_int(&RM, N);
    dealloc_1d_double(&time_series);
    dealloc_2d_double(&embedded_time_series,N);
    return 0;
}

/**
 * @brief Computes basic scalar measures for a given time series.
 * Calculates variance, volatility, and dominant frequency (via FFT).
 * * @param ts The input time series.
 * @param index The run ID, used for writing to output files.
 */
void compute_all_ts_scalars(double *ts, int index)
{
    // === Basic statistics: Mean, Variance, Volatility ===
    double mean = 0.0, variance = 0.0, autocorr = 0.0, volatility = 0.0;
    for (int i = 0; i < N; i++) mean += ts[i];
    mean /= N;

    for (int i = 0; i < N - 1; i++) {
        double diff = ts[i + 1] - ts[i];
        volatility += fabs(diff);
        variance += (ts[i] - mean) * (ts[i] - mean);
        autocorr += (ts[i] - mean) * (ts[i + 1] - mean); // Autocorrelation lag-1
    }
    variance += (ts[N - 1] - mean) * (ts[N - 1] - mean);
    volatility /= (N - 1);
    variance /= N;

    // === Dominant Frequency using GSL's Fast Fourier Transform ===
    gsl_fft_real_wavetable *real;
    gsl_fft_real_workspace *work;
    double *ts_copy;
    alloc_1d_double(&ts_copy, N);
    memcpy(ts_copy, ts, N * sizeof(double));

    real = gsl_fft_real_wavetable_alloc(N);
    work = gsl_fft_real_workspace_alloc(N);
    gsl_fft_real_transform(ts_copy, 1, N, real, work);

    // Find the frequency with the highest magnitude (power)
    double max_mag = 0.0;
    int max_idx = 0;
    for (int i = 1; i < N / 2; i++) { // Ignore DC component (i=0) and Nyquist frequency
        double mag = fabs(ts_copy[i]);
        if (mag > max_mag) {
            max_mag = mag;
            max_idx = i;
        }
    }
    double dominant_freq = (double)max_idx / N;

    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
    dealloc_1d_double(&ts_copy);

    // === Write results to files (appending) ===
    FILE *var_out = fopen(filename1, "a");
    FILE *vol_out = fopen(filename3, "a");
    FILE *df_out = fopen(filename5, "a");

    if (!var_out || !vol_out || !df_out) {
        fprintf(stderr, "Error writing scalar output files.\n");
        return;
    }

    fprintf(var_out, "%d %f\n", index, variance);
    fprintf(vol_out, "%d %f\n", index, volatility);
    fprintf(df_out, "%d %f\n", index, dominant_freq);

    fclose(var_out);
    fclose(vol_out);
    fclose(df_out);
}

/**
 * @brief Runs a basic time series analysis over an ensemble of simulations.
 * NOTE: As currently written, this function uses the single global 'mobility' variable to read all data files.
 * It is designed to analyze the statistical distribution of scalar measures (e.g., variance) for a set of simulations
 * that were all generated with the same mobility parameter.
 */
int ensemble_time_series_analysis()
{
    int ensemble_size = 10000;

    // Prepare output folder and filenames
    size_t folder_name_size = strlen(results_path) + 256;
    char *folder_name = (char *)malloc(folder_name_size * sizeof(char));
    snprintf(folder_name, folder_name_size, "%sTS_Analysis/", results_path);
    create_folders_in_path_if_not_exists(folder_name);

    sprintf(filename1, "%svariance_N=%d.dat", folder_name, N);
    sprintf(filename3, "%svolatility_N=%d.dat", folder_name, N);
    sprintf(filename5, "%sfreq_N=%d.dat", folder_name, N);

    // Clear the files once before appending in the loop
    FILE *f;
    f = fopen(filename1, "w"); fclose(f);
    f = fopen(filename3, "w"); fclose(f);
    f = fopen(filename5, "w"); fclose(f);

    // Loop over the ensemble
    for (int initial = 1; initial <= ensemble_size; initial++)
    {
        print_prog((double)(initial + 1) / (double)ensemble_size);

        double *ts;
        alloc_1d_double(&ts, N);

        if (read_file_for_time_series(ts, mobility, initial) != 0) {
            fprintf(stderr, "Skipping index %d due to read error.\n", initial);
            dealloc_1d_double(&ts);
            continue;
        }

        compute_all_ts_scalars(ts, initial);
        dealloc_1d_double(&ts);
    }

    free(folder_name);
    return 0;
}

//================================================================================
// Main Function
//================================================================================
int main() 
{
    // For timing the execution
    struct timeval start, end;
    double time_elapsed;
    srand(time(NULL));
    gettimeofday(&start, NULL);

    // --- Select the analysis to run by uncommenting ONE of the following lines ---

    // ensemble_analysis_varying_mobility(); // Primary analysis: RQA vs. mobility
    // ensemble_recurrence_analysis();       // RQA stats for an ensemble of simulations
    single_recurrence_analysis();      // RQA on a single file for debugging/plotting
    // ensemble_time_series_analysis();      // Basic stats (variance, etc.) across an ensemble of simulations

    // Stop timer and print duration
    gettimeofday(&end, NULL);
    time_elapsed = get_time_elapsed_minutes(start, end);
    printf("Run duration: %1.1f minutes\n", time_elapsed);

    return 0;
}
