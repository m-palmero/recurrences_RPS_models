# ifndef _auxiliary_functions_h
# define _auxiliary_functions_h

void print_prog(double percentage);
void copy (double *x, double *y, int dim);
int alloc_1d_double(double **x, int n);
int alloc_2d_double(double ***x, int n, int m);
int alloc_2d_int(int ***x, int n, int m);
int alloc_1d_int(int **x, int n);
int dealloc_1d_double(double **x);
int dealloc_1d_int(int **x);
int dealloc_2d_double(double ***x, int n);
int dealloc_2d_int(int ***x, int n);
double constrainAngle(double *x);
double moving_window(double *x, int n, int m);
void create_folder_if_not_exists(const char* folderName);
int create_folders_in_path_if_not_exists(const char* path);
void remove_dot(char *result, size_t resultSize, double number);
double get_time_elapsed_minutes(struct timeval start, struct timeval end);
int contains_nan(double **array, int rows, int cols);
int check_timeout(clock_t start_time, int time_out_limit);
void build_data_filename(char *buffer, size_t buffer_size, const char *base_path, const char *prefix, int index, const char *extension);
int strcasecmp_trim(const char *a, const char *b);
int field_matches(const char *a, const char *b);
# endif