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



int main (int argc, char * argv[]) {

  const double pi=3.141592653589793;
  double  * rho;
  struct lc_cell lc_environment;
  double  tf=50.0;
  double time, dz,  dt=1e-3;
  double timeprint=0.2;
  FILE * time_file, * snapshot_file;
  const char * initial_conditions="standard";
  const char * time_file_name="rho_bottom_middle_top.dat";
  const char * output_file_name="rho_time.dat";
  int timesteper_kind_flag=0;
  int nz;


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
  
  //Read the parameter values form the input file:
  parse_input_file(  & lc_environment,  & tf, & timeprint , & dt );
  print_log_file( lc_environment, tf, dt, "log.file");


  nz=lc_environment.nz;
  dz=lc_environment.cell_length/(nz-1);
  lc_environment.dz=dz;
  time=lc_environment.ti;
  
  //Starting the PDE solver:
  gsl_odeiv2_system sys = {frank_energy, jacobian, nz+2, &lc_environment};


  //Choose the integrator:
  gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-9, 0.0);
  //gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf, 1e-8, 1e-8, 0.0);


  gsl_odeiv2_driver_set_hmax (pde_driver , dt );


  time_file=fopen(time_file_name,"w");
  snapshot_file=fopen(lc_environment.output_file_name,"w");

  
  rho= (double *) malloc( (nz+2)*sizeof(double) );


  if ( strcmp(lc_environment.initial_conditions,"standard") == 0 )
    {
      for (int ii=1; ii<=nz;ii++)
	{

	  rho[ii]=lc_environment.rho0;
	  
	};
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

    fprintf(time_file,"#time rho[bottom]  rho[middle] rho[top] \n");
  fprintf(time_file,"%f  %f  %f  %f\n",time, rho[0],rho[nz/2], rho[nz]);

  
  print_snapshot_to_file(rho,time,dz,nz,snapshot_file);


  printf("time=%lf/%lf\n",time,tf);	
  while(time <tf)
    {

      int status = gsl_odeiv2_driver_apply (pde_driver, &time, time+timeprint, rho);      


      if (status != GSL_SUCCESS)
	{

	  printf ("error, return value=%d\n", status);
   
	};

      printf("time=%lf/%lf\n",time,tf);
      print_snapshot_to_file(rho,time,dz,nz,snapshot_file);
      fprintf(time_file,"%f  %f  %f  %f\n",time, rho[0],rho[(nz)/2], rho[nz+1]);
      
	
    };
  

  gsl_odeiv2_driver_free (pde_driver);
  free(rho);
  fclose(time_file);
  fclose(snapshot_file);
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
  double surf_viscosity[2];
  double tau[2], kappa[2];
  double drho, d2rho;

  
  tau[0]=mu.tau[0];
  tau[1]=mu.tau[1];
  
  
  kappa[0]=mu.kappa[0];
  kappa[1]=mu.kappa[1];
  
  

  /*Bulk equations */  
  for(int ii=1; ii<=nz; ii++)
    {

      d2rho=(rho[ii+1]+rho[ii-1]-2.0*rho[ii])/(dz*dz);
      drho=(rho[ii+1]-rho[ii-1])/(2*dz);

      Rhs[ii]= D_c*d2rho;
          

    };

  Rhs[0]=0;


  Rhs[nz+1]=0;

      

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
  
  tau[0]=mu.tau[0];
  tau[1]=mu.tau[1];
  
  
  kappa[0]=mu.kappa[0];
  kappa[1]=mu.kappa[1];
  



  
  

  gsl_matrix_set_zero( &dRhsdrho_mat.matrix );
  
  for(int ii=0; ii<nz+2;ii++)
    {

      dRhsdt[ii]=0;
      
    };
  
  
  for(int i=1;i<nz+2;i++){

    gsl_matrix_set ( &dRhsdrho_mat.matrix,i,i-1,k/(dz*dz) );
    gsl_matrix_set ( &dRhsdrho_mat.matrix,i,i  ,-2.0*k/(dz*dz));
    gsl_matrix_set ( &dRhsdrho_mat.matrix,i,i+1,k/(dz*dz) );

  };


  gsl_matrix_set ( &dRhsdrho_mat.matrix,0,0,-(k/dz)) ;
  gsl_matrix_set ( &dRhsdrho_mat.matrix,0,1,k/(dz));

  gsl_matrix_set( &dRhsdrho_mat.matrix,nz-1,nz-2,k/(dz) );		  
  gsl_matrix_set( &dRhsdrho_mat.matrix,nz-1,nz-1,(-(k/dz) ));
    
  return GSL_SUCCESS;
  
};


int print_snapshot_to_file(const double * rho,
			   const double time,
			   const double dz,
			   const int nz,
			   FILE * snapshot_file)
{

      fprintf(snapshot_file,"#time=%f\n",time);

      for(int ii=1;ii<nz+2;ii++)
	{
	  
	  fprintf(snapshot_file,"%f  %f\n",ii*dz,rho[ii]);
      

	};
      fprintf(snapshot_file,"\n\n");

};

int print_rho_time( const double * rho,
		    const double time,
		    const double dz,
		    const int nz)
			       
{
  FILE * snapshot_file;
  char snap_name[50];

  sprintf(snap_name,"time=%f.dat",time);
  snapshot_file=fopen(snap_name,"w");

  fprintf(snapshot_file,"#i  nx        ny         nz         rho   time=%f, dz=%f\n",time,dz);

  for(int ii=0;ii<nz;ii++)
    {

      fprintf(snapshot_file,"%i  %f  %f  %f  %f\n",ii,cos(rho[ii]),sin(rho[ii]), 0.0, rho[ii]);
      

    };

  fclose(snapshot_file);
  
  return 0;
};


void print_log_file(const struct lc_cell lc,
		    const double  tf,
		    const double  dt,
		    const char something[])
{

  printf("\n\n Parameters values used:\n\n");
  printf( "K, alpha, D_c:                        %lf  %lf  %lf\n",lc.k,lc.alpha,lc.D_c );
  printf( "maximum timestep (dt):      %lf \n",dt);
  printf( "Number of Layers(Nz):       %d  \n", lc.nz);
  printf( "cell length:                %lf \n",lc.cell_length);  
  printf( "Simulation time:            %lf  \n\n",tf);
    
};



 
