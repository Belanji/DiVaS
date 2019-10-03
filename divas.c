#include <stdio.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "divas.h"
#include "parser.h"

const static double pi=3.141592653589793;
const lz=2;


int main (int argc, char * argv[]) {

  double  * rho;
  struct lc_cell lc_environment;
  double  tf=50.0;
  double time, dz,  dt=1e-3;
  double timeprint=0.2;
  FILE * time_file, * snapshot_file;
  const char * initial_conditions="standard";
  const char * time_file_name="sigma_time.dat";
  const char * output_file_name="rho_time";
  int timesteper_kind_flag=0;
  int nz;
  int  snapshot_number=0;
  double total_particles;
  
  //Standard values:
  strcpy(lc_environment.initial_conditions,initial_conditions);
  strcpy(lc_environment.output_file_name,output_file_name);


  lc_environment.k=1.0;
  lc_environment.alpha=1.0;
  lc_environment.D_c=1.0;
  lc_environment.cell_length=1.0;

  lc_environment.tau[0]=1.0;
  lc_environment.tau[1]=1.0;

  lc_environment.kappa[0]=1.0;
  lc_environment.kappa[1]=1.0;


  lc_environment.ti=0.;
  lc_environment.tf=50.;
  lc_environment.dt=0.2;

  lc_environment.rho0=1;
  lc_environment.sigma0[0]=0;
  lc_environment.sigma0[1]=0;
  

  
  //Read the parameter values form the input file:
  parse_input_file(  & lc_environment,  & tf, & timeprint , & dt );
  print_log_file( lc_environment, tf, dt, "log.file");


  nz=lc_environment.nz;
  dz=lz/(nz-1);
  lc_environment.dz=dz;
  time=lc_environment.ti;
  
  //Starting the PDE solver:
  gsl_odeiv2_system sys = {RhsFunction, jacobian, nz+2, &lc_environment};


  //Choose the integrator:
  gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-9, 0.0);
  //gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf, 1e-8, 1e-8, 0.0);


  gsl_odeiv2_driver_set_hmax (pde_driver , dt );  
  rho= (double *) malloc( (nz+2)*sizeof(double) );


  if ( strcmp(lc_environment.initial_conditions,"standard") == 0 )
    {


      rho[0]=lc_environment.sigma0[0];
      for (int ii=1; ii<nz+1;ii++)
	{

	  rho[ii]=lc_environment.rho0;
	  
	}
      rho[nz+1]=lc_environment.sigma0[1];
    }
  else if ( strcmp(lc_environment.initial_conditions,"read_from_file") == 0 || strcmp(lc_environment.initial_conditions,"ic_file") == 0)
    {



      int i,j,k,ii,jj,kk;
      double trash_double;	
      FILE * ic_file;
      char string_placeholder[400];
      int read_status;
      int reading_line=1;
  
      ic_file=fopen(lc_environment.ic_file_name,"r");
      if (ic_file== NULL)
	{
	  printf("Unable to find the file \"%s\".\nPlease check your initial condition file name.\n\nAborting the program.\n\n",lc_environment.ic_file_name);
	  exit(0);
	}

      //get the file header:
  
      printf("\nReading initial conditions from \"%s\".\n",lc_environment.ic_file_name);
  
      //removing the header:
      fgets(string_placeholder,400,ic_file);
      reading_line++;


      //Let's work:

      for(k= 0; k< nz; k++)
	{
	  fgets(string_placeholder,400,ic_file);
	  read_status=sscanf(string_placeholder,"%lf %lf\n",&trash_double,&rho[k]);
	  //read_check(read_status,reading_line);


      	      
	  reading_line++;
	}
      
  
    }
  else 
    {

      printf("No initial condition named %s is defined.\nAborting the program.\n\n",lc_environment);
      exit(0);
  
    };

  printf("\n\nStarting calculations\n\n");


  total_particles=calculate_total_particle_quantity ( rho,
						      & lc_environment);
  
  time_file=fopen(time_file_name,"w");
  fprintf(time_file,"#time   sigma_b   sigma_t total_particles\n");
  fprintf(time_file,"%e  %e  %e  %e\n",time, rho[0],rho[nz+1],total_particles);
  
  print_snapshot_to_file(rho,time,lz,dz,nz,output_file_name,snapshot_number);
  printf("snapshot %d: %lf\n",snapshot_number,time);
  snapshot_number++;

  while(time <tf)
    {

      int status = gsl_odeiv2_driver_apply (pde_driver, &time, time+timeprint, rho);      


      if (status != GSL_SUCCESS)
	{

	  printf ("error, return value=%d\n", status);
   
	};

      printf("snapshot %d: %lf\n",snapshot_number,time);
      print_snapshot_to_file(rho,time,lz,dz,nz,output_file_name,snapshot_number);
      snapshot_number++;

      total_particles=calculate_total_particle_quantity ( rho,
						      & lc_environment);
      fprintf(time_file,"%e  %e  %e  %e \n",time, rho[0],rho[nz+1], total_particles);
      fflush(time_file);
	
    };
  
  fprintf(time_file,"\n");
  gsl_odeiv2_driver_free (pde_driver);
  free(rho);
  fclose(time_file);
  return 0;


};
     

    
int RhsFunction (double t, const double rho[], double Rhs[], void * params)
{ 
  struct lc_cell mu = *(struct lc_cell *)params;
  int nz=mu.nz;
  double dz = mu.cell_length/(nz-1);
  double k= mu.k;
  double alpha=mu.alpha;
  double D_c=mu.D_c;
  double tau[2], kappa[2];
  double drho, d2rho, dsigma;
  double GhostRho;
  double z_position;
  
  tau[0]=mu.tau[0];
  tau[1]=mu.tau[1];
  
  
  kappa[0]=mu.kappa[0];
  kappa[1]=mu.kappa[1];
  
  
    /*bottom boundary equations */

  z_position=-lz/2;
  dsigma=kappa[0]*rho[1]-rho[0]/tau[0];
  
  //Extrapolate ghost point:
  GhostRho=rho[2]-2*dz*dsigma/D_c*(1.0+alpha*cos(k*z_position));
  

  drho=(rho[2]-GhostRho)/(2*dz);
  d2rho=(rho[2]+GhostRho-2.0*rho[1])/(dz*dz);
  
  Rhs[0]=dsigma;
  Rhs[1]=D_c*(1.0+alpha*cos(z_position))*d2rho-D_c*alpha*k*sin(k*z_position)*drho;


  /*Bulk equations */
  
  for(int ii=2; ii<nz+1; ii++)
    {

      z_position=-lz/2.+k*dz*(ii-1);
      d2rho=(rho[ii+1]+rho[ii-1]-2.0*rho[ii])/(dz*dz);
      drho=(rho[ii+1]-rho[ii-1])/(2*dz);

      Rhs[ii]= D_c*(1.0+alpha*cos(z_position))*d2rho-D_c*alpha*k*sin(k*z_position)*drho;
          

    };

  
  /* Top boundary equations*/

  dsigma=kappa[1]*rho[nz]-rho[nz+1]/tau[1];
  GhostRho=rho[nz-1]-2*dz*dsigma/D_c*(1.0+alpha);
    
  
    drho=(GhostRho-rho[nz-1])/(2*dz);
    d2rho=(GhostRho+rho[nz-1]-2.0*rho[nz])/(dz*dz);
  
  Rhs[nz]=D_c*(1.0+alpha*cos(k*lz/2))*d2rho-D_c*alpha*k*sin(k*lz/2)*drho;
          
  Rhs[nz+1]=dsigma;

      

      return GSL_SUCCESS;

    };


int jacobian(double t, const double rho[], double * dRhsdrho, double dRhsdt[], void * params)
{
struct lc_cell mu = *(struct lc_cell *)params;
  int nz=mu.nz;
  double dz = mu.cell_length/(nz-1);
  double k= mu.k;
  double alpha=mu.alpha;
  double D_c=mu.D_c;
  double surf_viscosity[2];
  double tau[2], kappa[2];
  double drho, d2rho;
  gsl_matrix_view dRhsdrho_mat= gsl_matrix_view_array (dRhsdrho, nz+2, nz+2);
  
  //tau[0]=mu.tau[0];
  //tau[1]=mu.tau[1];
  //
  //
  //kappa[0]=mu.kappa[0];
  //kappa[1]=mu.kappa[1];
  //
  //
  //
  //
  //
  //
  //
  //gsl_matrix_set_zero( &dRhsdrho_mat.matrix );
  //
  //for(int ii=0; ii<nz+2;ii++)
  //  {
  //
  //    dRhsdt[ii]=0;
  //    
  //  };
  //
  //
  //for(int i=1;i<nz+2;i++){
  //
  //  gsl_matrix_set ( &dRhsdrho_mat.matrix,i,i-1,k/(dz*dz) );
  //  gsl_matrix_set ( &dRhsdrho_mat.matrix,i,i  ,-2.0*k/(dz*dz));
  //  gsl_matrix_set ( &dRhsdrho_mat.matrix,i,i+1,k/(dz*dz) );
  //
  //};
  //
  //
  //gsl_matrix_set ( &dRhsdrho_mat.matrix,0,0,-(k/dz)) ;
  //gsl_matrix_set ( &dRhsdrho_mat.matrix,0,1,k/(dz));
  //
  //gsl_matrix_set( &dRhsdrho_mat.matrix,nz-1,nz-2,k/(dz) );		  
  //gsl_matrix_set( &dRhsdrho_mat.matrix,nz-1,nz-1,(-(k/dz) ));
    
  return GSL_SUCCESS;
  
};


int print_snapshot_to_file(const double * rho,
			   const double time,
                           const double lz,
			   const double dz,
			   const int nz,
                           const char * output_file_prefix,
			   int  snapshot_number)
{

  FILE * snapshot_file;
  char output_file_name[200];

  
  sprintf(output_file_name,"%s_%d.dat",output_file_prefix,snapshot_number);

  snapshot_file=fopen(output_file_name,"w");
  fprintf(snapshot_file,"#z  rho(z)\n");

  
  for(int ii=1; ii<nz+1;ii++)
    {
	  
      fprintf(snapshot_file,"%e  %e\n",(ii-1)*dz-lz/2,rho[ii]);
      

    };
  fprintf(snapshot_file,"\n");
  
  fclose(snapshot_file);
};



void print_log_file(const struct lc_cell lc,
		    const double  tf,
		    const double  dt,
		    const char something[])
{

  printf("\n\nParameters values used:\n\n");
  printf( "cell length:      %lf   \n",lc.cell_length);
  printf( "Number of Layers(Nz):       %d  \n", lc.nz);
  printf( "K, alpha, D_c:  %lf  %lf  %lf\n",lc.k,lc.alpha,lc.D_c );
  printf("rho_0:  %lf\n", lc.rho0);

  printf("\nBoundary conditions:\n\n");
  
  printf( "sigma0_b, tau_b, kappa_b:  %lf  %lf  %lf\n",lc.sigma0[0],lc.tau[0],lc.kappa[0] );
  printf( "sigma0_t, tau_t, kappa_t:  %lf  %lf  %lf\n",lc.sigma0[1],lc.tau[1],lc.kappa[1] );


  printf("\nTime parameters:\n\n");
  printf( "maximum timestep (dt):      %lf \n",dt);
  printf( "Simulation time:            %lf  \n\n",tf);
    
};



double calculate_total_particle_quantity ( const double rho[],
					   const void  * params)
{
  struct lc_cell mu = *(struct lc_cell *)params;
  const int nz=mu.nz;
  const double lz=mu.cell_length;
  const double dz = mu.cell_length/(nz-1);
  double total_particle_quantity;


  total_particle_quantity=rho[0]+rho[nz+1];

  total_particle_quantity+=rho[1]*dz/2.;
  for(int ii=2; ii<nz;ii++)
    {

      total_particle_quantity+=dz*rho[ii];

    }
  total_particle_quantity+=rho[nz]*dz/2.;
  
  return total_particle_quantity;
}
