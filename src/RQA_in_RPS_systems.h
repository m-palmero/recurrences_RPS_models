#ifndef RQA_in_RPS_systems_H
#define RQA_in_RPS_systems_H

int read_file_for_time_series(double *time_series, double mobility, int initial);
int give_time_series(double *time_series);
int give_recurrence_plot(int **R, int n);
int ensemble_recurrence_analysis();
int ensemble_analysis_varying_mobility();
int single_recurrence_analysis();
void compute_all_ts_scalars(double *ts, int initial);
int ensemble_time_series_analysis();

# endif
