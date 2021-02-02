/*============================================================================ 
| (c) Resources Allocation in Virtualized Cluster
============================================================================*/ 

#include "ga.h" 

int eval_vms();    /*--- Forward declaration ---*/ 
                     /*these variables can be used in any part of this program*/ 

#define MAXVMS 100    /* The maximum number of tasks */ 
struct{ 
	float time, 
	      mem,
               cpu; 
} VM[MAXVMS]; 

int Num_VMs; 
int Num_Nodes;
float Node_CPU;
float Node_MEM;
float no_of_nodes;
float tot_time;

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
   ga_info = GA_config("vmtest.cfg", eval_vms); 
   printf("GA config read successfully\n"); 

   read_vms(ga_info->user_data); 
   
   ga_info->chrom_len = Num_VMs; 

if(argc > 1) 
 {X_select(ga_info, argv[1]); 
 }; 
 
 /*--- Run the GA ---*/ 
   printf("Running GA\n"); 
   GA_run(ga_info); 
   printf("GA run completed\n"); 
printf("Total used nodes=%G\n",no_of_nodes);
   printf("Total Time is=%G\n",tot_time);
} 

/*---------------------------------------------------------------------------- 
| obj_fun() - user specified objective function 
----------------------------------------------------------------------------*/ 
int eval_vms(chrom) 
   Chrom_Ptr chrom; 
{ 
   int i; 
float fit_fun,max_time;
float tot_mem, vm_time, vm_mem, vm_cpu, tot_cpu; 

/* Trivial case no VMs */ 
   if(chrom->length < 1)  { 
	chrom->fitness = 0; 
	    return; 
   } 

/* Initialization*/ 

tot_time = 0.0; 
max_time = 0.0; 
tot_mem = 0.0; 
tot_cpu = 0.0; 
no_of_nodes =1.0;

/*--Place each Task using next fit Its 2-D Multi Capacity bin Packing---- */ 
for(i = 0; i < chrom->length; i++)  {
	no_of_nodes +=1; 
	vm_time = VM[(int)chrom->gene[i]-1].time; 
   	vm_mem = VM[(int)chrom->gene[i]-1].mem; 
            vm_cpu = VM[(int)chrom->gene[i]-1].cpu;
  
 /* ---Place VM on Bins----*/ 
   
   if(vm_mem + tot_mem > Node_MEM ||  vm_cpu + tot_cpu > Node_CPU ) {/*--Too much memory or CPU---*/ 
       tot_mem = vm_mem;
       tot_cpu   = vm_cpu;
       tot_time += vm_time; 
       max_time = 0;	 
      } 
	else { 
	   tot_time += vm_time;
               tot_mem += vm_mem;
               tot_cpu += vm_cpu;
	} 

/*__Find the longest time at the node--*/ 

if(vm_time > max_time) { 
max_time = vm_time; 
}} 
fit_fun = no_of_nodes*tot_time;
chrom->fitness = fit_fun; 
} 

/*--------------------------------------------------------------------------------------------- 
|read data from file 
---------------------------------------------------------------------------------------------------*/ 

read_vms(filename) 
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

 /*--get number of Nodes----*/ 
fscanf(fid,"%d",&Num_Nodes);
printf("Number of free nodes in the server are:,%d\n", Num_Nodes);

  /*--get number of VMs----*/ 
fscanf(fid,"%d", &Num_VMs); 
  
  if(Num_VMs < 1 || Num_VMs > MAXVMS) 
  { 
    printf("Number of VMS ,%d, out of bounds [1...%d]\n", Num_VMs, MAXVMS); 
    exit(1); 
  } 
printf("Number of VMS to deploy are,%d\n", Num_VMs);

/*get the CPUand Memory constraints for any node*/

fscanf(fid,"%f%f", &Node_CPU, &Node_MEM);


/*----get the VMs time and the memory requirement -----*/ 

  for(i=0; i < Num_VMs; i++) 
   
     fscanf(fid,"%f%f%f", &VM[i].time, &VM[i].cpu, &VM[i].mem); 
     
  

/*---Close the file--*/ 
fclose(fid); 
} 
