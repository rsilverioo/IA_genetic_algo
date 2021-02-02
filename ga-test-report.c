/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
| Genetic Algorithm Test Program 
============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ga.h"

/* Global Variables*/
char graph[500][500]; // could be optimised using malloc...
int nnodes, nedges, prev_best;
int read_once = 1;
float best_x = 0, best_mu = 0;

/* File variables */
FILE *fp;
Pool_Ptr rep_pool;

/* Function prototypes */
int obj_fun(Chrom_Ptr);
int read_instance();

/*----------------------------------------------------------------------------
| main()
----------------------------------------------------------------------------*/
int main()
{
  GA_Info_Ptr ga_info;
  int i, j, k, count;

  /*--- x_rate variation ---*/
  prev_best = 1000;

  // Create .csv file
  fp = fopen("x_rate.csv", "w");
  fprintf(fp, "x_rate, Best, Min, Max, Ave, Tot, Var, SD");
  fclose(fp);

  printf("Starting cross rate change simulations \n");
  for (k = 0; k <= 20; k++)
  {
    /*--- Initialize the genetic algorithm ---*/
    ga_info = GA_config("GAconfig", obj_fun);

    // Fill the matrix "graph" with all the info
    // Also initialize nnodes and nedges
    read_instance(ga_info->user_data);

    // Changing chromosome length
    ga_info->chrom_len = nnodes;
    ga_info->mu_rate = 0.6;

    /*--- Change the config ---*/
    ga_info->rand_seed = rand() % 22000 + 10000;
    ga_info->x_rate = k / 20.0f;

    /*--- Run the GA ---*/
    GA_run(ga_info);

    count = 0;
    for (i = 0; i < ga_info->chrom_len; i++)
    {
      if (ga_info->best->gene[i])
        count++;
    }

    printf("  - x_rate = %f | Media nodos: %d (fitness: %g)\n", ga_info->x_rate, count, ga_info->best->fitness);

    // Save best config for report
    if (prev_best < count)
    {
      prev_best = count;
      best_x = ga_info->x_rate;
    }

    // Add stats to csv file
    fp = fopen("mu_rate.csv", "a");
    rep_pool = ga_info->old_pool;
    fprintf(fp, "\n%f, %d, %G, %G, %.2G, %G, %.2G, %.2G",
            ga_info->x_rate, count, rep_pool->min, rep_pool->max, rep_pool->ave, rep_pool->total_fitness, rep_pool->var, rep_pool->dev);
    fclose(fp);
  }

  /*--- mu_rate variation ---*/
  read_once = 1;
  prev_best = 1000;

  // Create .csv file
  fp = fopen("mu_rate.csv", "w");
  fprintf(fp, "mu_rate, Best, Min, Max, Ave, Tot, Var, SD");
  fclose(fp);

  printf("\nStarting mutation rate change simulations \n");
  for (k = 0; k <= 20; k++)
  {
    /*--- Initialize the genetic algorithm ---*/
    ga_info = GA_config("GAconfig", obj_fun);

    // Fill the matrix "graph" with all the info
    // Also initialize nnodes and nedges
    read_instance(ga_info->user_data);

    // Changing chromosome length
    ga_info->chrom_len = nnodes;
    ga_info->x_rate = 0.6;

    /*--- Change the config ---*/
    ga_info->rand_seed = rand() % 22000 + 10000;
    ga_info->mu_rate = k / 20.0f;

    /*--- Run the GA ---*/
    GA_run(ga_info);

    count = 0;
    for (i = 0; i < ga_info->chrom_len; i++)
    {
      if (ga_info->best->gene[i])
        count++;
    }

    printf("  - mu_rate = %f | Media nodos: %d (fitness: %g)\n", ga_info->mu_rate, count, ga_info->best->fitness);

    // Saving best config for report
    if (prev_best < count)
    {
      prev_best = count;
      best_mu = ga_info->mu_rate;
    }

    // Add stats to csv file
    fp = fopen("mu_rate.csv", "a");
    rep_pool = ga_info->old_pool;
    fprintf(fp, "\n%f, %d, %G, %G, %.2G, %G, %.2G, %.2G",
            ga_info->mu_rate, count, rep_pool->min, rep_pool->max, rep_pool->ave, rep_pool->total_fitness, rep_pool->var, rep_pool->dev);
    fclose(fp);
  }

  printf("\nFinished!\n");
  printf("Best config: x_rate = %f, mu_rate = %f \n\n", best_x, best_mu);

  printf("Press ENTER to close the window.");
  getchar();
}

/*----------------------------------------------------------------------------
| obj_fun() - user specified objective function
----------------------------------------------------------------------------*/
int obj_fun(Chrom_Ptr chrom)
{
  int i, j;
  double val = 0.0;
  int v = 0;
  int a = 0;

  for (i = 0; i < chrom->length; i++)
  {
    v += chrom->gene[i];
  }

  for (i = 0; i < chrom->length; i++)
  {
    if (chrom->gene[i] == 1)
    {
      for (j = i + 1; j < chrom->length; j++)
      {
        if (chrom->gene[j] == 1)
        {
          a += graph[i][j];
        }
      }
    }
  }

  // Target functions
  int alpha = 1;

  // Function 1
  //chrom -> fitness = v - alpha*(v*(v-1)/2-a);

  // Function 2
  //chrom -> fitness = v - alpha*sqrt(pow(v*(v-1)/2, 2)-pow(a, 2));

  // Function 3
  //chrom -> fitness = v + 1/(v*(v-1)/2 - a + 0.01);

  // Function 4
  //chrom -> fitness = (v*(v-1)/2 - a);

  // Function 5
  chrom->fitness = (v * (v - 1) / 2 - a) + 1 / (pow(v, 2) + 0.01);

  return 0;
}

/*----------------------------------------------------------------------------
| read_instance() - read DIMACS format
----------------------------------------------------------------------------*/
int read_instance(char *filename)
{
  char dummy1;
  char dummy2[100];
  int dummy3;
  int n1, n2;
  FILE *inputf;
  int i, j;

  nnodes = 0;
  nedges = 0;

  if ((inputf = fopen(filename, "rt")) == NULL)
  {
    printf("Cannot open file %s\n", filename);
    exit(-1);
  }

  // Read header
  fscanf(inputf, "%c %s %d %d\n", &dummy1, dummy2, &nnodes, &nedges);

  if (read_once)
  {
    read_once = 0;
    printf("Opening %s (%d nodes, %d edges)\n", filename, nnodes, nedges);
  }

  for (i = 0; i < nnodes; i++)
    for (j = 0; j < nnodes; j++)
      graph[i][j] = 0;

  // Skip node list
  for (i = 0; i < nnodes; i++)
    fscanf(inputf, "%c  %d %d\n", &dummy1, &dummy3, &dummy3);

  // Read all edges
  for (i = 0; i < nedges; i++)
  {
    fscanf(inputf, "%c %d %d\n", &dummy1, &n1, &n2);
    graph[n1 - 1][n2 - 1] = 1;
    graph[n2 - 1][n1 - 1] = 1;
  }

  fclose(inputf);
}