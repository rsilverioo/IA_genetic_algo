/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
|
| Genetic Algorithm Test Program 
============================================================================*/

#include "ga.h"

int next_fit();    /*--- Forward declaration ---*/

#define MAXPKGS 100  /*these variables can be used in any part of this program*/

float Pkgs[MAXPKGS],     /*packages*/
      Sum_Pkgs;          /*sum of packages weight*/
 

int Num_Pkgs;   /*Actual number of packages*/

/*----------------------------------------------------------------------------
| main()
----------------------------------------------------------------------------*/
main(argc, argv)
   int  argc;
   char *argv[];
{
   GA_Info_Ptr ga_info;

   /*--- Initialize the genetic algorithm ---*/
   printf("Reading GA config\n");
   ga_info = GA_config("bp.cfg", next_fit);
   printf("GA config read successfully\n");

   read_packages(ga_info->user_data);
   
   ga_info->chrom_len = Num_Pkgs;

if(argc > 1)
 {X_select(ga_info, argv[1]);
 };
   /*--- Run the GA ---*/
   printf("Running GA\n");
   GA_run(ga_info);
   printf("GA run completed\n");
   printf("Sum of packages weight = %f\n\n", Sum_Pkgs);
}

/*----------------------------------------------------------------------------
| obj_fun() - user specified objective function
----------------------------------------------------------------------------*/
int next_fit(chrom) 
   Chrom_Ptr chrom;
{
   int i, num_bins;
   float pkg_weight, weight;

   if(chrom->length < 1)
   {chrom->fitness = 0;
    return;
   }
   num_bins = 1;
   weight = 0;

   /*--- Fitness is number of genes out of place for sorted order ---*/
   for(i = 0; i < chrom->length; i++) 
{  pkg_weight = Pkgs[(int)chrom->gene[i]-1];
   
   if(weight + pkg_weight > 1.0)
      { 
       weight = pkg_weight;
        num_bins++;
      } else
	{weight += pkg_weight;
	}
   }

   chrom->fitness = num_bins;


}
/*----------------------------------------------------------------------------------------------------
|read data from file
---------------------------------------------------------------------------------------------------*/

read_packages(filename)
 char *filename;
{ 
  FILE *fid;
  int i;
 
  /*--open datafile----*/
  if((fid = fopen(filename,"r")) == NULL) 
  {
   printf("Errors in opening the packages datafile <%s>\n", filename);
   exit(1);
  }

  /*--get number of packages----*/
  fscanf(fid,"%d", &Num_Pkgs);
  
  if(Num_Pkgs < 1 || Num_Pkgs > MAXPKGS)
  {
    printf("Number of packages ,%d, out of bounds [1...%d]\n", Num_Pkgs, MAXPKGS);
    exit(1);
  }

  /*----get the packages weight and sum -----*/
  Sum_Pkgs = 0;
  for(i=0; i < Num_Pkgs; i++)
  {
     fscanf(fid,"%f", &Pkgs[i]);
     Sum_Pkgs += Pkgs[i];
  }

/*---Close the file--*/
fclose(fid);
}

