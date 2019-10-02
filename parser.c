#include "divas.h"
#include "parser.h"
#include <string.h>
#include <stdlib.h>
const static double pi=3.141592653589793;

void error_check(int error_handler,
		  char parser[])

{

  	  if (error_handler <= 0 )
	    {
	    printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	  printf("Please review your input file.\n Aborting the program\n");
	  exit(0);
	    }
	  
}




void parse_input_file(struct lc_cell  * lc,
		      struct optical_setup  * opt,
		      double * tf,
		      double * timeprint,
		      double * dt )
{

  char parser[80];
  char garbage[400];
  int error_handler;

  while (   scanf("%79s",parser) != EOF )
    {


      if ( strcmp(parser,"k11") == 0 || strcmp(parser,"K11") == 0 )
	{

	  error_handler=scanf("%lf",&lc->k11);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
		
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcmp(parser,"k22")== 0 || strcmp(parser,"K22")== 0 )
	{
	

	  error_handler=scanf("%lf",&lc->k22);


	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
		  
	  
	  fgets(garbage,400,stdin);
	  
	  
      
	}
      else if ( strcmp(parser,"k33") == 0 || strcmp(parser,"K33")== 0 )
	{

	  error_handler=scanf("%lf",&lc->k33);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcmp(parser,"n0") == 0 || strcmp(parser,"N0")== 0 )
	{

	  error_handler=scanf("%lf",&lc->n0);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);


	}
            else if ( strcmp(parser,"lambda") == 0 || strcmp(parser,"wavelength")== 0 || strcmp(parser,"wave_length")== 0)
	{

	  error_handler=scanf("%lf",&opt->lambda);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);


	}
            else if ( strcmp(parser,"ne") == 0 || strcmp(parser,"Ne")== 0 )
	{

	  error_handler=scanf("%lf",&lc->ne);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);


	}
      else if (strcmp(parser,"ti") == 0 || strcmp(parser,"start_time") ==0 )
	{


	  error_handler=scanf("%lf",&lc->ti);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if (strcmp(parser,"tf") == 0 || strcmp(parser,"run_time") ==0 )
	{


	  error_handler=scanf("%lf",tf);


	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);

	    };


	  fgets(garbage,400,stdin);


	}      
      else if (strcmp(parser,"dt") == 0 || strcmp(parser,"maximum_timestep") ==0 )
	{


	  error_handler=scanf("%lf",dt);


	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);

	    };


	  fgets(garbage,400,stdin);


	}
      else if ( strcmp(parser,"timeprint") == 0 || strcmp(parser,"print_every") ==0 )
	{


	  error_handler=scanf("%lf", timeprint);


	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };


	  fgets(garbage,400,stdin);


	}
      else if (strcmp(parser,"bulk_viscosity") == 0  || strcmp(parser,"bviscosity") == 0 || strcmp(parser,"bulk_visc") == 0 || strcmp(parser,"gamma1") == 0 || strcmp(parser,"gamma_1") == 0 )
	{


	  error_handler=scanf("%lf",&lc->viscosity);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };

	  
	  fgets(garbage,400,stdin);
	  
	}
      else if (strcmp(parser,"bottom_surface_viscosity") == 0  || strcmp(parser,"bottom_sviscosity") == 0 ||  strcmp(parser,"bsurf_visc") == 0 || strcmp(parser,"b_gammas") == 0 )
	{

	  
	  error_handler=scanf("%lf",&lc->surf_viscosity[0]);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"top_surface_viscosity") == 0 || strcmp(parser,"top_sviscosity") == 0 ||  strcmp(parser,"tsurf_visc") == 0 || strcmp(parser,"t_gammas") == 0 )
	{

	  
	  error_handler=scanf("%lf",&lc->surf_viscosity[1]);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"surface_viscosity") == 0 || strcmp(parser,"sviscosity") == 0 ||  strcmp(parser,"surf_visc") == 0 || strcmp(parser,"gammas") == 0 || strcmp(parser,"gamma_s") == 0 )
	{

	  
	  error_handler=scanf("%lf",&lc->surf_viscosity[0]);
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);

	    };
	  


	  error_handler=scanf("%lf",&lc->surf_viscosity[1]);
	  if (error_handler <= 0 )
	    {

	      printf("The keyword '%s' must be followed by 2 numeric arguments.\n",parser);
	      printf("The first is the bottom %s, while the second is the top %s.\n",parser,parser);	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);

	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"cell_length") == 0  || strcmp(parser,"cell_height") == 0 ||  strcmp(parser,"lenght_z") == 0 || strcmp(parser,"d")==0  )
	{

	  
	  error_handler=scanf("%lf",&lc->cell_length);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
            else if (strcmp(parser,"chiral_power") == 0  || strcmp(parser,"q0") == 0 || strcmp(parser,"q") == 0 )
	{

	  
	  error_handler=scanf("%lf",&lc->q);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"pitch") == 0 || strcmp(parser,"Pitch") == 0  || strcmp(parser,"P0") == 0 || strcmp(parser,"p0") == 0)
	{

	  double lc_pitch;
	  error_handler=scanf("%lf",&lc_pitch);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };

	  lc->q=2.0*pi/lc_pitch;
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"azimutal_ancoring") == 0  || strcmp(parser,"wa" ) == 0  )
	{

	  
	  error_handler=scanf("%lf",&lc->wa[0]);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };

	  error_handler=scanf("%lf",&lc->wa[1]);


	  if (error_handler <= 0 )
	    {

	      printf("The keyword '%s' must be followed by 2 numeric arguments.\n",parser);
	      printf("The first is the bottom %s, while the second is the top %s.\n",parser,parser);	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  
	}
      else if (strcmp(parser,"omega_d") == 0  || strcmp(parser,"omega" ) == 0 || strcmp(parser,"omega" ) == 0 || strcmp(parser,"omegad" ) == 0  )
	{
	  
	  error_handler=scanf("%lf",&lc->omega_d[0]);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };

	  error_handler=scanf("%lf",&lc->omega_d[1]);

	  if (error_handler <= 0 )
	    {

	      printf("The keyword '%s' must be followed by 2 numeric arguments.\n",parser);
	      printf("The first is the bottom %s, while the second is the top %s.\n",parser,parser);	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  
	}
      else if (strcmp(parser,"nz") == 0  || strcmp(parser,"Nz" ) == 0 || strcmp(parser,"NZ" ) == 0 || strcmp(parser,"nZ" ) == 0  )
	{
	  
	  error_handler=scanf("%d",&lc->nz);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  fgets(garbage,400,stdin);
	  
	}
      else if ( strcmp(parser,"initial_conditions") == 0 )
	{

	  error_handler=scanf("%200s",&lc->initial_conditions);
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);
	  

	}
      else if ( strcmp(parser,"ic_file") == 0 || strcmp(parser,"initial_conditions_file") == 0 || strcmp(parser,"input_initial_conditions") == 0 || strcmp(parser,"input_initial_conditions_file") == 0)
	{

	  error_handler=scanf("%200s",&lc->ic_file_name);
	  lc->ic_file_flag=1;
	  error_check(error_handler,parser);
	  fgets(garbage,400,stdin);


	}
      else if ( strcmp(parser,"output_file_name") == 0 ||
		strcmp(parser,"output_file") == 0 )
	{

	  error_handler=scanf("%s",&(lc->output_file_name));
	  error_check(error_handler,parser);
		
	  
	  fgets(garbage,400,stdin);


	}
      else if (parser[0]=='#')
	{

	  fgets(garbage,400,stdin);

	}	  
      else
	{

	  printf("The parser did not recognize the option '%s'. Please review your input file\n", parser);
	  printf("Aborting the program\n");
	  exit(0);
	};
      
    };
};
