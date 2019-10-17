#ifndef director__

#define director__
#include  <stdio.h>
#include <complex.h>
struct lc_cell
{

  double k, alpha, tau;
  double tau_k[2], tau_d[2];
  double ti, tf, dt;
  double rho0;
  double sigma0[2];
  double dz;
  int nz;
  char output_file_name[200];
  char initial_conditions[200];
  char ic_file_name[200];
  int ic_file_flag;

};

int RhsFunction (double t,
		  const double rho[],
		  double Rhs[],
		  void  * params);

double calculate_total_particle_quantity ( const double rho[],
					   const void  *params);


int jacobian(double t,
	     const double rho[],
	     double * dRhsdrho,
	     double dRhsdt[],
	     void * params);



int print_snapshot_to_file(const double *,
			   const double  ,
			   const double  ,
			   const int     ,
                           const char *  ,
			   int   );

void print_log_file(const struct lc_cell,
		    const double ,
		    const double ,		    
		    const char []);


void print_sigma_time(const struct lc_cell lc,
		    const double * rho,
		    const double  time,
                      FILE * time_file);


double calculate_average_rho ( const double rho[],
                               const void  * params);


double calculate_average_rho_z_1 ( const double rho[],
                                   const void  * params);

double calculate_average_rho_z_2 ( const double rho[],
                                   const void  * params);

#endif
