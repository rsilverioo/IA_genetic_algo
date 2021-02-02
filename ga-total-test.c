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
int nnodes, nedges;
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
int main(int argc, char *argv[])
{
  GA_Info_Ptr ga_info;
  int i, j, count;
  int node_array[50];

  char *a = argv[1];
  int i_index = atoi(a) - 1;

  // Creamos y/o abrimos el .csv
  if ((fp = fopen("madness.csv", "r")) == NULL)
  {
    fp = fopen("madness.csv", "w");
    fprintf(fp, "Index, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00");
  }
  else
  {
    fclose(fp);
    fp = fopen("madness.csv", "a");
  }

  for (j = 0; j <= 20; j++)
  {
    /*--- Initialize the genetic algorithm ---*/
    ga_info = GA_config("GAconfig", obj_fun);

    // Fill the matrix "graph" with all the info
    // Also initialize nnodes and nedges
    read_instance(ga_info->user_data);

    // Changing chromosome length
    ga_info->chrom_len = nnodes;

    /*--- Change the config ---*/
    ga_info->rand_seed = rand() % 22000 + 10000;
    ga_info->x_rate = i_index / 20.0f;
    ga_info->mu_rate = j / 20.0f;

    printf("\nmu_rate = %f | x_rate = %f", ga_info->mu_rate, ga_info->x_rate);

    /*--- Run the GA ---*/
    GA_run(ga_info);

    count = 0;
    for (i = 0; i < ga_info->chrom_len; i++)
    {
      if (ga_info->best->gene[i])
        count++;
    }

    node_array[j] = count;
  }

  // Add stats to csv file
  fprintf(fp, "\n%f, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d",
          ga_info->x_rate, node_array[0], node_array[1], node_array[2], node_array[3], node_array[4], node_array[5], node_array[6], node_array[7], node_array[8], node_array[9],
          node_array[10], node_array[11], node_array[12], node_array[13], node_array[14], node_array[15], node_array[16], node_array[17], node_array[18], node_array[19], node_array[20]);
  fclose(fp);

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