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
const double lz=2.;


int main (int argc, char * argv[]) {

  double  * rho;
  struct lc_cell lc_environment;
  double  tf=50.0;
  double time, dz,  dt=1e-3;
  double timeprint=0.2;
  FILE * time_file, * snapshot_file;
  const char * initial_conditions="delta";
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

  lc_environment.tau=1.0;

  lc_environment.tau_k[0]=1.0;
  lc_environment.tau_k[1]=1.0;

  lc_environment.tau_a[0]=1.0;
  lc_environment.tau_a[1]=1.0;

  
  lc_environment.tau_d[0]=1.0;
  lc_environment.tau_d[1]=1.0;

  lc_environment.sigma0[0]=0;
  lc_environment.sigma0[1]=0;

  
  lc_environment.ti=0.;
  lc_environment.tf=50.;
  lc_environment.dt=2e-3;
  lc_environment.rho0=1;
  

  
  //Read the parameter values form the input file:
  parse_input_file(  & lc_environment,  & tf, & timeprint , & dt );
  print_log_file( lc_environment, tf, dt, "log.file");


  nz=lc_environment.nz;
  dz=lz/(nz-1);
  lc_environment.dz=dz;
  time=lc_environment.ti;
  
  //Starting the PDE solver:
  gsl_odeiv2_system sys = {RhsFunction, jacobian, nz+4, &lc_environment};


  //Choose the integrator:
  //gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-9, 0.0);
  gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf, 1e-8, 1e-8, 0.0);


  gsl_odeiv2_driver_set_hmax (pde_driver , dt );  




  rho= (double *) malloc( (nz+4)*sizeof(double) );
  /*Important note: simga0 are placed at rho[0] and rho[nz+3]
   while dsigma_dt are placed at rho[1] and rho[nz+2].
  The rho values are placed between 2 and nz+1.*/
  

  if (    strcasecmp(lc_environment.initial_conditions,"standard") == 0
       || strcasecmp(lc_environment.initial_conditions,"homogeneous") == 0)
    {


      rho[0]=lc_environment.sigma0[0];
      rho[1]=0;
            
      for (int ii=2; ii<nz+2;ii++)
	{

	  rho[ii]=lc_environment.rho0;
	  
	}
      rho[nz+2]=0;
      rho[nz+3]=lc_environment.sigma0[1];
    }
  else if (    strcasecmp(lc_environment.initial_conditions,"delta") == 0)
    {


      rho[0]=lc_environment.sigma0[0];
      rho[1]=0;
            
      for (int ii=2; ii<nz+2;ii++) rho[ii]=0;

      rho[nz/2 +2]=2*lc_environment.rho0/dz;
	
	  	  	
      rho[nz+2]=0;
      rho[nz+3]=lc_environment.sigma0[1];
    }
  else if ( strcasecmp(lc_environment.initial_conditions,"read_from_file") == 0 || strcmp(lc_environment.initial_conditions,"ic_file") == 0)
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
  fprintf(time_file,"%e  %e  %e  %e\n",time, rho[0],rho[nz+3],total_particles);
  
  print_snapshot_to_file(rho,time,dz,nz,output_file_name,snapshot_number);
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
      print_snapshot_to_file(rho,time,dz,nz,output_file_name,snapshot_number);
      snapshot_number++;

      total_particles=calculate_total_particle_quantity ( rho,
						      & lc_environment);
      fprintf(time_file,"%e  %e  %e  %e \n",time, rho[0],rho[nz+3], total_particles);
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
  double dz = lz/(nz-1);
  double k=pi*mu.k;
  double alpha=mu.alpha;
  const double tau=mu.tau;
  double tau_d[2], tau_k[2], tau_a[2];
  double drho, d2rho, dsigma, sigma, drho_dt;
  double GhostRho;
  double z_position;
  
  tau_d[0]=mu.tau_d[0];
  tau_d[1]=mu.tau_d[1];
  
  tau_k[0]=mu.tau_k[0];
  tau_k[1]=mu.tau_k[1];

  tau_a[0]=mu.tau_a[0];
  tau_a[1]=mu.tau_a[1];


  /*bottom boundary equations */

  z_position=-lz/2;

  sigma=rho[0];
  dsigma=rho[1];
  
  
  //Extrapolate ghost point:
  GhostRho=rho[3]-2*dz*dsigma/(1.0+alpha*cos(k*z_position));
  

  drho=(rho[3]-GhostRho)/(2*dz);
  d2rho=(rho[3]+GhostRho-2.0*rho[2])/(dz*dz);
//
  drho_dt=(1.0+alpha*cos(k*z_position))*d2rho-alpha*k*sin(k*z_position)*drho;

  
  Rhs[0]=dsigma;
  Rhs[1]=(tau_d[0]*tau_d[0]/16.)*(4*drho_dt/(tau_d[0]*tau_k[0])
                                  -4*dsigma/(tau_d[0]*tau_a[0])
                                  +rho[2]/(tau_a[0]*tau_k[0])
				  -exp(-t*tau_d[0]/(tau_a[0]*4.))*sigma/(tau*tau_a[0]));

  Rhs[2]=drho_dt;

  /*Bulk equations */
  
  for(int ii=3; ii<nz+1; ii++)
    {

      z_position=-lz/2.+dz*(ii-2);
      d2rho=(rho[ii+1]+rho[ii-1]-2.0*rho[ii])/(dz*dz);
      drho=(rho[ii+1]-rho[ii-1])/(2*dz);

      Rhs[ii]=(1.0+alpha*cos(k*z_position))*d2rho-alpha*k*sin(k*z_position)*drho;
          

    };

  
  /* Top boundary equations*/

  z_position=lz/2;
  dsigma=rho[nz+2];
  GhostRho=rho[nz]-2*dz*dsigma/(1.0+alpha*cos(k*z_position));
    
  
  drho=(GhostRho-rho[nz])/(2*dz);
  d2rho=(GhostRho+rho[nz]-2.0*rho[nz+1])/(dz*dz);

  drho_dt=(1.0+alpha*cos(k*z_position))*d2rho-alpha*k*sin(k*z_position)*drho;
  
  Rhs[nz+1]=drho_dt;
  Rhs[nz+2]=(tau_d[1]*tau_d[1]/16.)*(4*drho_dt/(tau_d[1]*tau_k[1])-4*dsigma/(tau_d[1]*tau_a[1])+rho[nz+1]/(tau_a[1]*tau_k[1])
				     -exp(-t*tau_d[1]/(4*tau_a[1]))*rho[nz+3]/(tau*tau_a[1]));
  
  Rhs[nz+3]=dsigma;


  Rhs[nz+1]=0;
  Rhs[nz+2]=0;
  Rhs[nz+3]=0;
  
  return GSL_SUCCESS;
      
    };


int jacobian(double t, const double rho[], double * dRhsdrho, double dRhsdt[], void * params)
{


  struct lc_cell mu = *(struct lc_cell *)params;
  int nz=mu.nz;
  double dz = lz/(nz-1);
  double k=pi*mu.k;
  double alpha=mu.alpha;
  const double tau=mu.tau;
  double tau_d[2], tau_k[2], tau_a[2];
  double drho, d2rho, dsigma, drho_dt;
  double GhostRho;
  double z_position;
  double dz_2=1/(dz*dz);
  double dz_1=1/dz;
  
  
  gsl_matrix_view dRhsdrho_mat= gsl_matrix_view_array (dRhsdrho, nz+4, nz+4);
  
  tau_d[0]=mu.tau_d[0];
  tau_d[1]=mu.tau_d[1];
    
  tau_k[0]=mu.tau_k[0];
  tau_k[1]=mu.tau_k[1];

  tau_a[0]=mu.tau_a[0];
  tau_a[1]=mu.tau_a[1];

  
  

  gsl_matrix_set_zero( &dRhsdrho_mat.matrix );

  double tau_d2=tau_d[0]*tau_d[0];
  double tau_d3=tau_d[0]*tau_d[0]*tau_d[0];
  double tau_a2=tau_a[0]*tau_a[0];
  
  dRhsdt[0]=tau_d3*exp(-t*tau_d[0]*0.25/tau_a[0])*rho[0]/(64*tau*tau_a2); 
  for(int ii=1; ii<nz+4;ii++)
    {
  
      dRhsdt[ii]=0;
      
    };

  z_position=-lz/2;

  //Boundary consitions:
  gsl_matrix_set ( &dRhsdrho_mat.matrix,0,0,-tau_d2*exp(-t*tau_d[0]*0.25/tau_a[0])/(16.*tau*tau_a[0] );


  gsl_matrix_set ( &dRhsdrho_mat.matrix,1,1  ,1);
  gsl_matrix_set ( &dRhsdrho_mat.matrix,2,2  ,1);


  sigma=rho[0];
  dsigma=rho[1];
  
  
  //Extrapolate ghost point:
  GhostRho=rho[3]-2*dz*rho[1]/(1.0+alpha*cos(k*z_position));
  

  drho=(rho[3]-GhostRho)/(2*dz);
  d2rho=(rho[3]+GhostRho-2.0*rho[2])/(dz*dz);
//
  drho_dt=(1.0+alpha*cos(k*z_position))*d2rho-alpha*k*sin(k*z_position)*drho;

  
  Rhs[0]=dsigma;
  Rhs[1]=(tau_d[0]*tau_d[0]/16.)*(4*drho_dt/(tau_d[0]*tau_k[0])
                                  -4*dsigma/(tau_d[0]*tau_a[0])
                                  +rho[2]/(tau_a[0]*tau_k[0])
				  -exp(-t*tau_d[0]/(tau_a[0]*4.))*sigma/(tau*tau_a[0]));

  Rhs[2]=drho_dt;

  
  for(int i=3;i<nz+1;i++)
    {

      gsl_matrix_set(  &dRhsdrho_mat.matrix,i,i-1, (1.0+alpha*cos(k*z_position))*dz_2 - alpha*k*sin(k*z_position)*(-dz_1*0.5) );
      gsl_matrix_set(  &dRhsdrho_mat.matrix,i,i  , (1.0+alpha*cos(k*z_position))*(-2*dz_2)  );
      gsl_matrix_set(  &dRhsdrho_mat.matrix,i,i+1, (1.0+alpha*cos(k*z_position))*dz_2 - alpha*k*sin(k*z_position)*(-dz_1*0.5) );

  
    };

  gsl_matrix_set ( &dRhsdrho_mat.matrix, nz+1, nz+1, 1);
  gsl_matrix_set ( &dRhsdrho_mat.matrix, nz+2, nz+2, 1);
  gsl_matrix_set ( &dRhsdrho_mat.matrix, nz+3, nz+3, 1);


  //gsl_matrix_set ( &dRhsdrho_mat.matrix,0,0,-(k/dz)) ;
  //gsl_matrix_set ( &dRhsdrho_mat.matrix,0,1,k/(dz));
  //
  //gsl_matrix_set( &dRhsdrho_mat.matrix,nz-1,nz-2,k/(dz) );		  
  //gsl_matrix_set( &dRhsdrho_mat.matrix,nz-1,nz-1,(-(k/dz) ));
    
  return GSL_SUCCESS;
  
};


int print_snapshot_to_file(const double * rho,
			   const double time,
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

  
  for(int ii=2; ii<nz+2;ii++)
    {
	  
      fprintf(snapshot_file,"%e  %e\n",(ii-2)*dz-lz/2,rho[ii]);
      

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

  printf( "Number of Layers(Nz):       %d  \n", lc.nz);
  printf( "k(in Pi units):  %lf \n",lc.k);
  printf( "alpha:  %lf \n",lc.alpha);
  printf( "tau:  %e  \n",lc.tau);

  
  printf("\nBoundary conditions:\n\n");

  printf(  "tau_a:  %e  \n",lc.tau_a[0] );
  printf(  "tau_d:  %e  \n",lc.tau_d[0] );
  printf(  "tau_k:  %e  \n",lc.tau_k[0] );




  printf("\nTime parameters:\n\n");
  printf( "maximum timestep (dt):      %e \n",dt);
  printf( "Simulation time:            %lf  \n\n",tf);
    
};
//
//
//
double calculate_total_particle_quantity ( const double rho[],
					   const void  * params)
{
  struct lc_cell mu = *(struct lc_cell *)params;
  const int nz=mu.nz;
  const double dz = lz/(nz-1);
  double total_particle_quantity;


  total_particle_quantity=rho[0]+rho[nz+3];

  total_particle_quantity+=rho[2]*dz/2.;
  for(int ii=3; ii<nz+1;ii++)
    {

      total_particle_quantity+=dz*rho[ii];

    }
  total_particle_quantity+=rho[nz+1]*dz/2.;
  
  return total_particle_quantity;
}
