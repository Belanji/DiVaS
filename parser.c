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
		      double * tf,
		      double * timeprint,
		      double * dt )
{

  char parser[120];
  char garbage[400];
  int error_handler;

  while (   scanf("%119s",parser) != EOF )
    {


      if ( strcasecmp(parser,"k") == 0 )
	{

	  error_handler=scanf("%lf",&lc->k);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
		
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"alpha")== 0  )
	{
	

	  error_handler=scanf("%lf",&lc->alpha);


	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
		  
	  
	  fgets(garbage,400,stdin);
	  
	  
      
	}
      else if ( strcasecmp(parser,"tau") == 0 )
	{

	  error_handler=scanf("%lf",&(lc->tau));
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"tau_d") == 0 )
	{

	  error_handler=scanf("%lf",&(lc->tau_d[0]));
	  lc->tau_d[1]=lc->tau_d[0];
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"tau_k") == 0 )
	{

	  error_handler=scanf("%lf",&(lc->tau_k[0]));
	  lc->tau_k[1]=lc->tau_k[0];
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcasecmp(parser,"sigma0") == 0 )
	{

	  error_handler=scanf("%lf",&(lc->sigma0[0]));
	  lc->sigma0[0]=lc->sigma0[1];
	  
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
