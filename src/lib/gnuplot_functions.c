# include "../lib/common.h"
# include "gnuplot_functions.h"

#define GNUPLOT "gnuplot --persist"

void plot_gnuplot_time_series_no_title(char png_name[], double parameter, int N)
{	
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	fprintf(gp, "set terminal png enhanced size 2500,1500 font \"CMR10\" 60 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");

	// fprintf(gp, "unset key \n");
	// fprintf(gp, "unset cbtics \n");
	// fprintf(gp, "unset xtics \n");
	// fprintf(gp, "unset ytics \n");
	// fprintf(gp, "unset colorbox \n");
	// fprintf(gp, "unset border \n");
	// fprintf(gp, "set lmargin 0 \n");
	// fprintf(gp, "set rmargin 0 \n");
	// fprintf(gp, "set tmargin 0 \n");
	// fprintf(gp, "set bmargin 0 \n");

	//fprintf(gp, "set title \"Time series test\\n{/*0.7 Control Parameter=%.5f; Max Iterations =%d;}\" \n", parameter, N); 

	//fprintf(gp, "set format y '' \n");
	//fprintf(gp, "set cbtics 0,500,1000 \n");
	//fprintf(gp, "set palette maxcolors 4\n");
	fprintf(gp, "set xrange[0:%d] \n", N);
	fprintf(gp, "set output '%s'\n", png_name);

	fprintf(gp, "plot 'ts.dat' u 0:1 w l lw 4 lc 'black' noti\n");
	//fprintf(gp, "pause -1");
	fclose(gp);
}

void plot_gnuplot_RPs_no_title(char *filename, int N)
{
	/* Start gnuplot */
	FILE *gp = popen("gnuplot", "w");
	if (!gp) {
	fprintf(stderr, "Error: could not open gnuplot pipe: %s\n", strerror(errno));
	}

	int resolution = 1000;
	/* Terminal & output */
	fprintf(gp,
	//"set terminal png enhanced crop size 2500,2500 font 'CMR10' 60 \n"
	"set terminal pngcairo enhanced size %d,%d font 'CMR10,20'\n"
	"set output '%s'\n",
	resolution, resolution, filename);

	/* Layout */
	fprintf(gp,
	"set size square\n"
	"unset key\n"
	"unset colorbox\n"
	"set xrange [0:%d]\n"
	"set yrange [0:%d]\n",
	N, N);

	// fprintf(gp, "unset key \n");
	// fprintf(gp, "unset xtics \n");
	// fprintf(gp, "unset ytics \n");
	// fprintf(gp, "unset cbtics \n");
	// fprintf(gp, "unset colorbox \n");
	// fprintf(gp, "unset border \n");
	// fprintf(gp, "set lmargin 0 \n");
	// fprintf(gp, "set rmargin 0 \n");
	// fprintf(gp, "set tmargin 0 \n");
	// fprintf(gp, "set bmargin 0 \n");


	/* Compute dynamic pointsize */
	double base_ps = (double) resolution / 1000;
	double ps = base_ps * sqrt(((double) resolution/100) / (double)N);

	/* Clip big symbols at the border */
    //fprintf(gp, "set clip points \n");

	/* Plot with dynamic PS */
	fprintf(gp,
	"plot 'rp.dat' using 1:2 "
	"with points pointtype 5 pointsize %.3f linecolor rgb 'black' notitle\n",
	ps);

	/* Close gnuplot */
	fclose(gp);
}