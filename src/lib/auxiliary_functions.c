# include "../lib/common.h"
# include "auxiliary_functions.h"

#define PBSTR "=================================================="
#define PBWIDTH 50

#define T 100 //Maximum characteres for file's names.

void print_prog(double percentage)
{
	int val = (int)(percentage * 100);
	int lpad = (int)(percentage * PBWIDTH);
	int rpad = PBWIDTH - lpad;
	printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
    printf("\n");
}

void copy(double *x, double *y, int dim)
{
  for (int i=0; i<dim; i++)	x[i] = y[i];
}

void change(double *x, double a, double b)
{
	x[1,0] = a;
	x[0,1] = b;
}

int alloc_1d_double(double **x, int n)
{
	*x = (double *)malloc(n * sizeof(double));
	return 0;
}

int alloc_2d_double(double ***x, int n, int m)
{
	*x = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (double *)malloc(m * sizeof(double));
	}
	return 0;
}

int alloc_2d_int(int ***x, int n, int m)
{
	*x = (int **)malloc(n * sizeof(int *));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (int *)malloc(m * sizeof(int));
	}
	return 0;
}

int alloc_1d_int(int **x, int n)
{
	*x = (int *)malloc(n * sizeof(int));
	return 0;
}


int dealloc_1d_double(double **x)
{
	free(*x);
	return 0;
}

int dealloc_1d_int(int **x)
{
	free(*x);
	return 0;
}

int dealloc_2d_double(double ***x, int n)
{
	for (int i = 0; i < n; i++)
	{
		free((*x)[i]);
	}
	free(*x);
	return 0;
}

int dealloc_2d_int(int ***x, int n)
{
	for (int i = 0; i < n; i++)
	{
		free((*x)[i]);
	}
	free(*x);
	return 0;
}

double constrainAngle(double *x)
{
    *x = fmod(*x, 2 * M_PI);
    if (*x < 0)
        *x += 2 * M_PI;
    return *x;
}

double moving_window(double *x, int n, int m)
{
	double aux;

	for (int i = 0; i < n - m; i++)
	{
    	aux = 0;
    	for (int j = 0; j < m; j++)
		{
        	aux += x[i+j];
    	}
    	x[i] = aux * (1.0 / m);
	}

	return *x;
}

void create_folder_if_not_exists(const char* folder_name) 
{
    // Check if the folder exists
    struct stat st;
    if (stat(folder_name, &st) == -1) {
        // Folder does not exist, create it
        mkdir(folder_name, 0700); // 0700 gives full permissions, you can adjust it accordingly
        printf("Folder '%s' created!\n", folder_name);
    }
}

int create_folders_in_path_if_not_exists(const char* path) 
{
    char temp_path[1024];
    char* p = NULL;
    struct stat st;
    size_t len;
    int status = 0;

    // Copy the path to a temporary variable
    snprintf(temp_path, sizeof(temp_path), "%s", path);
    len = strlen(temp_path);

    // Iterate through the path, creating each part
    for (p = temp_path + 1; p < temp_path + len; p++) {
        if (*p == '/') {
            *p = '\0';  // Temporarily end the string here

            // Create the directory if it doesn't exist
            if (stat(temp_path, &st) == -1) {
                if (mkdir(temp_path, 0700) == 0) {
                    printf("Directory created: %s\n", temp_path);
                } else {
                    perror("Error creating folder");
                    return -1;  // Error creating folder
                }
            }

            *p = '/';  // Restore the original character
        }
    }

    // Create the final directory (if not created in the loop)
    if (stat(temp_path, &st) == -1) {
        if (mkdir(temp_path, 0700) == 0) {
            printf("Full path: %s\n", temp_path);
        } else {
            perror("Error creating folder");
            return -1;  // Error creating folder
        }
    }

    return status;
}

void remove_dot(char *result, size_t resultSize, double number)
{
    char buffer[50];

    // Convert double to string with high precision
    snprintf(buffer, sizeof(buffer), "%.16f", number);

    // Initialize result to empty string
    result[0] = '\0';

    // Copy characters from buffer to result, skipping the dot
    int j = 0;
    for (int i = 0; i < strlen(buffer) && j < resultSize - 1; i++) {
        if (buffer[i] != '.') {
            result[j++] = buffer[i];
        }
    }
    
    // Ensure null termination
    result[j] = '\0';
}

double get_time_elapsed_minutes(struct timeval start, struct timeval end) 
{
    long seconds, useconds;
    double timeElapsed;

    // Calculate the time elapsed in microseconds
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    timeElapsed = seconds + useconds / 1000000.0;

    // Convert time elapsed to minutes
    return timeElapsed / 60.0;
}

int contains_nan(double **array, int rows, int cols) 
{
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (isnan(array[i][j])) {
                return 1; // Found NaN
            }
        }
    }
    return 0; // No NaNs
}

int check_timeout(clock_t start_time, int time_out_limit) 
{
    clock_t current_time = clock();
    double elapsed_time = (double)(current_time - start_time) / CLOCKS_PER_SEC;
    if (elapsed_time > time_out_limit) {
        fprintf(stderr, "Timeout reached after %f seconds\n", elapsed_time);
        return 1; // Timeout reached
    }
    return 0; // No timeout
}

void build_data_filename(char *buffer, size_t buffer_size, const char *base_path, const char *prefix, int index, const char *extension)
{
    snprintf(buffer, buffer_size, "%s%s%d%s", base_path, prefix, index, extension);
}

// Compare two strings case-insensitively, ignoring leading/trailing whitespace
int strcasecmp_trim(const char *a, const char *b) {
    while (isspace(*a)) a++;
    while (isspace(*b)) b++;
    size_t la = strlen(a);
    size_t lb = strlen(b);
    while (la > 0 && isspace(a[la-1])) la--;
    while (lb > 0 && isspace(b[lb-1])) lb--;
    if (la != lb) return 1;
    for (size_t i = 0; i < la; i++)
        if (tolower(a[i]) != tolower(b[i])) return 1;
    return 0;
}
    // Returns 1 if the strings are the same (case-insensitive, whitespace-insensitive)
int field_matches(const char *a, const char *b) {
    while (isspace(*a)) a++;
    while (isspace(*b)) b++;
    size_t la = strlen(a);
    size_t lb = strlen(b);
    while (la > 0 && isspace(a[la-1])) la--;
    while (lb > 0 && isspace(b[lb-1])) lb--;
    if (la != lb) return 0;
    for (size_t i = 0; i < la; i++)
        if (tolower(a[i]) != tolower(b[i])) return 0;
    return 1;
}