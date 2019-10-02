#ifndef director__

#define director__
#include  <stdio.h>
#include <complex.h>
struct lc_cell
{

  double k11,k22,k33;
  double n0, ne;
  double ti, tf, dt;
  double viscosity, surf_viscosity[2];
  double pretwist[2], pretilt[2];
  double wa[2], omega_d[2];
  double cell_length;
  double dz;
  double q;
  int nz;
  char output_file_name[200];
  char initial_conditions[200];
  char ic_file_name[200];
  int ic_file_flag;

};


struct optical_setup
{
  double complex Pol[2][2];
  double complex Anal[2][2];
  double complex Ei[2][2];
  double lambda;
};

int frank_energy (double t,
		  const double phi[],
		  double f[],
		  void  * params);

int jacobian(double t,
	     const double phi[],
	     double * dfdphi,
	     double dfdt[],
	     void * params);


int print_phi_time( const double *,
		    const double  ,
		    const double  ,
		    const int );

int print_snapshot_to_file(const double *,
			   const double  ,
			   const double  ,
			   const int     ,
			   FILE   *);

void print_log_file(const struct lc_cell,
		    const struct optical_setup,
		    const double ,
		    const double ,		    
		    const char []);

double optical_transmitance (const double *phi,
			     const double * theta,
			     const int ,
			     const struct lc_cell * lc,
			     struct optical_setup * opt);


#endif
