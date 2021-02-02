

/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| Chromosome management
|
| Functions:
|    CH_alloc()  - allocate a chrom
|    CH_resize() - resize a chrom
|    CH_free()   - deallocate a chrom
|    CH_valid()  - is a chrom valid?
|    CH_reset()  - reset a chrom
|    CH_copy()   - copy a chrom over another
|    CH_cmp()    - compare two chromosomes
|    CH_print()  - print a chrom
|    CH_verify() - ensure chrom makes sense
============================================================================*/


#define MAXTOK 10  /* Maximum number of tokens on a line */
#define STRLEN 80  /* Length of an input line */


/* Number of chromosome pointers to alloc at a time */
#define PL_ALLOC_SIZE 10 


#include "ga.h"



int GA_generational(), GA_steady_state();

static Chrom_Ptr child1, child2;


int X_simple(), X_uniform(), X_order1(), X_order2(), X_pos(), X_cycle(), 
    X_pmx(), X_uox(), X_rox(), X_asex();


int MU_simple_invert(), MU_simple_random(), MU_swap();
 /* rnd float in [0..1] -- introduced by claudio 10/02/2004 */
int MU_float_random(), MU_float_rnd_pert(), MU_float_LS(), MU_float_gauss_pert();

double gaussian_random(void);

int obj_fun(   Chrom_Ptr chrom);

int RE_append(), RE_by_rank(), RE_first_weaker(), RE_weakest();

int SE_uniform_random(), SE_roulette(), SE_rank_biased();




/*----------------------------------------------------------------------------
| Allocate a chromosome
----------------------------------------------------------------------------*/
Chrom_Ptr CH_alloc(
   int length)
{
   Chrom_Ptr chrom;

   /*--- Error check ---*/
   if(length <= 0) UT_error("CH_alloc: invalid length");

   /*--- Allocate memory for chromosome ---*/
   chrom = (Chrom_Ptr)calloc(1, sizeof(Chrom_Type));
   if(chrom == NULL) UT_error("CH_alloc: chrom alloc failed");
   chrom->length = length;

   /*--- Allocate memory for genes ---*/
   chrom->gene = (Gene_Ptr)calloc(length, sizeof(Gene_Type));
   if(chrom->gene == NULL) UT_error("CH_alloc: gene alloc failed");

   /*--- Put in magic cookie ---*/
   chrom->magic_cookie = CH_cookie;

   /*--- Reset the chromosome ---*/
   CH_reset(chrom);

   return chrom;
}

/*----------------------------------------------------------------------------
| Resize a chromosome
|
| NOTE: after the resize, the chromosome may no longer be valid, 
|       therefore it is reset
----------------------------------------------------------------------------*/
CH_resize(
   Chrom_Ptr chrom,
   int       length)
{
   /*--- Error check ---*/
   if(!CH_valid(chrom)) UT_error("CH_resize: invalid chrom");
   if(length <= 0) UT_error("CH_resize: invalid length");

   /*--- Reallocate memory for genes ---*/
   chrom->gene = (Gene_Ptr)realloc(chrom->gene, length * sizeof(Gene_Type));
   if(chrom->gene == NULL) UT_error("CH_resize: gene realloc failed");
   chrom->length = length;

   /*--- Reset the chromosome ---*/
   CH_reset(chrom);
}

/*----------------------------------------------------------------------------
| De-Allocate a chromosome
----------------------------------------------------------------------------*/
CH_free(
   Chrom_Ptr chrom)
{
   /*--- Error check ---*/
  if(!CH_valid(chrom)) return 0;  //CC

   /*--- Free memory for genes ---*/
   if(chrom->gene != NULL) {
      free(chrom->gene);
      chrom->gene = NULL;
   }

   /*--- Put in NULL magic cookie ---*/
   chrom->magic_cookie = NL_cookie;

   /*--- Free memory for chromosome ---*/
   free(chrom);
}

/*----------------------------------------------------------------------------
| Is a chromosome valid, i.e., has it been allocated by CH_alloc()?
----------------------------------------------------------------------------*/
CH_valid(
   Chrom_Ptr chrom)
{
   /*--- Check for NULL pointers ---*/
   if(chrom == NULL) return FALSE;
   if(chrom->gene == NULL) return FALSE;

   /*--- Check for magic cookie ---*/
   if(chrom->magic_cookie != CH_cookie) return FALSE;

   /*--- Otherwise valid ---*/
   return TRUE;
}

/*----------------------------------------------------------------------------
| Reset a chromosome
----------------------------------------------------------------------------*/
CH_reset(
   Chrom_Ptr chrom)
{
   int i;

   /*--- Error check ---*/
   if(!CH_valid(chrom)) UT_error("CH_reset: invalid chrom");

   /*--- Initialize genes ---*/
   for(i=0; i<chrom->length; i++)
      chrom->gene[i] = (Gene_Type)0;

   /*--- Initialize chromosome ---*/
   chrom->fitness  = 0.0;
   chrom->ptf      = 0.0;
   chrom->index    = -1;
   chrom->idx_min  = 0;
   chrom->idx_max  = chrom->length;
   chrom->parent_1 = -1;
   chrom->parent_2 = -1;
   chrom->xp1      = -1;
   chrom->xp2      = -1;
}

/*----------------------------------------------------------------------------
| Copy a chromosome
----------------------------------------------------------------------------*/
CH_copy(
   Chrom_Ptr src, Chrom_Ptr dst)
{
   Gene_Ptr gene;

   /*--- Error check ---*/
   if(!CH_valid(src)) UT_error("CH_copy: invalid src");
   if(!CH_valid(dst)) UT_error("CH_copy: invalid dst");

   /*--- Resize if necessary ---*/
   if(dst->length != src->length) CH_resize(dst, src->length);

   /*--- Save memory pointed to by gene ---*/
   gene = dst->gene;

   /*--- Copy chrom ---*/
   memcpy(dst, src, sizeof(Chrom_Type));

   /*--- Restore memory pointed to by gene ---*/
   dst->gene = gene;

   /*--- Copy gene ---*/
   memcpy(dst->gene, src->gene, src->length * sizeof(Gene_Type));
}

/*----------------------------------------------------------------------------
| Compare chromosomes A and B:
|    -1  = A is better
|     1  = B is better
|     0  = Same
----------------------------------------------------------------------------*/
CH_cmp(
   GA_Info_Ptr ga_info,
   Chrom_Ptr   a, Chrom_Ptr b)
{
   if(ga_info->minimize)
      if(a->fitness < b->fitness)
         return -1;
      else if(a->fitness > b->fitness)
         return 1;
      else
         return 0;
   else
      if(a->fitness > b->fitness)
         return -1;
      else if(a->fitness < b->fitness)
         return 1;
      else
         return 0;
}

/*----------------------------------------------------------------------------
| Print a chromosome
----------------------------------------------------------------------------*/
CH_print(
   Chrom_Ptr chrom)
{
   int i;

   /*--- Error check ---*/
   if(!CH_valid(chrom)) UT_error("CH_print: invalid chrom");

   printf("==============================================================\n");
   printf("\nChrom: \n");
   for(i=0; i<chrom->length; i++)
      printf("%G ", chrom->gene[i]);
   printf("\n\n");
   printf("fitness = %G, ptf = %G, index = %d, idx_min = %d, idx_max = %d\n", 
      chrom->fitness, chrom->ptf, chrom->index, chrom->idx_min, chrom->idx_max);
   printf("parent_1 = %d, parent_2 = %d, xp1 = %d, xp2 = %d\n", 
      chrom->parent_1, chrom->parent_2, chrom->xp1, chrom->xp2);
   printf("==============================================================\n");
}

/*----------------------------------------------------------------------------
| Verify a chromosome
----------------------------------------------------------------------------*/
CH_verify(
   GA_Info_Ptr ga_info,
   Chrom_Ptr chrom)
{
   int i;
   char *allele_count, err_str[80];

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("CH_verify: invalid ga_info");
   if(!CH_valid(chrom)) UT_error("CH_verify: invalid chrom");

   /*--- Check for invalid length ---*/
   if(chrom->length <= 0) {
      CH_print(chrom);
      UT_error("CH_verify: bad length");
   }

   /*--- Check for invalid idx_min ---*/
   if(chrom->idx_min < 0 || chrom->idx_min > chrom->length) {
      CH_print(chrom);
      UT_error("CH_verify: idx_min out of bounds");
   }

   /*--- Check for invalid idx_max ---*/
   if(chrom->idx_max < 0 || chrom->idx_max > chrom->length) {
      CH_print(chrom);
      UT_error("CH_verify: idx_max out of bounds");
   }

   /*--- Check for invalid permutation ---*/
   if(ga_info->datatype == DT_INT_PERM) {

      /*--- Allocate allele_count vector ---*/
      allele_count = calloc(chrom->length, sizeof(char));
      if(allele_count == NULL) 
         UT_error("CH_verify: cannot alloc allele_count");
   
      /*--- Check each gene in the chromosome ---*/
      for(i=0; i<chrom->length; i++) {

         /*--- Check for allele out of bounds ---*/
         if(chrom->gene[i] < 1 || (int)chrom->gene[i] > chrom->length) {
            CH_print(chrom);
            sprintf(err_str,"CH_verify: gene[%d] = %G is out of bounds", 
                    i, chrom->gene[i]);
            UT_error(err_str);

         /*--- Check for duplicate alleles ---*/
         } else if(++(allele_count[((int)chrom->gene[i])-1]) > 1) {
            CH_print(chrom);
            sprintf(err_str,"CH_verify: gene[%d] = %G is a duplicate", 
                    i, chrom->gene[i]);
            UT_error(err_str);
         }
      }
      free(allele_count);
   }
}
/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| Configuration file and ga_info management
|
| Functions:
|    CF_alloc()    - allocate a ga_info
|    CF_free()     - deallocate a ga_info
|    CF_valid()    - is ga_info valid?
|    CF_reset()    - reset config
|    CF_report()   - print out current config
|    CF_read()     - read config file
|    CF_tokenize() - convert input line to tokens
|    CF_verify()   - ensure ga_info makes sense
============================================================================*/

/*----------------------------------------------------------------------------
| Allocate a GA_Info structure
----------------------------------------------------------------------------*/
GA_Info_Ptr CF_alloc() 
{
   GA_Info_Ptr ga_info;

   /*--- Allocate memory for a ga_info ---*/
   ga_info = (GA_Info_Ptr)calloc(1, sizeof(GA_Info_Type));
   if(ga_info == NULL) UT_error("CF_alloc: alloc failed");

   /*--- Ensure these are NULL ---*/
   ga_info->old_pool = NULL;
   ga_info->new_pool = NULL;
   ga_info->best     = NULL;

   /*--- Put in a magic cookie ---*/
   ga_info->magic_cookie = CF_cookie;

   /*--- Reset ga_info ---*/
   CF_reset(ga_info);

   return ga_info;
}

/*----------------------------------------------------------------------------
| De-Allocate a GA_Info structure
----------------------------------------------------------------------------*/
CF_free(
   GA_Info_Ptr ga_info)
{
   /*--- Error check ---*/
  if(!CF_valid(ga_info)) return GA_ERROR;  //CCC

   /*--- Free pools ---*/
   if(ga_info->old_pool != NULL) PL_free(ga_info->old_pool);
   if(ga_info->new_pool != NULL) PL_free(ga_info->new_pool);
   ga_info->old_pool = ga_info->new_pool = NULL;

   /*--- Free best chrom ---*/
   if(ga_info->best != NULL) CH_free(ga_info->best);
   ga_info->best = NULL;

   /*--- Put in a NULL cookie ---*/
   ga_info->magic_cookie = NL_cookie;

   /*--- Free ga_info ---*/
   free(ga_info);
}

/*----------------------------------------------------------------------------
| Is a ga_info valid, i.e., has it been allocated by CF_alloc()?
----------------------------------------------------------------------------*/
CF_valid(
   GA_Info_Ptr ga_info)
{
   /*--- Check for NULL pointers ---*/
   if(ga_info == NULL) return FALSE;

   /*--- Check for magic cookie ---*/
   if(ga_info->magic_cookie != CF_cookie) return FALSE;

   /*--- Otherwise valid ---*/
   return TRUE;
}

/*----------------------------------------------------------------------------
| Set all ga_info to defaults
----------------------------------------------------------------------------*/
CF_reset(
   GA_Info_Ptr ga_info)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("CF_reset: invalid ga_info");

   /*--- Default basic parameters ---*/
   ga_info->user_data[0]    = '\0';
   ga_info->function_index  = 1;
   ga_info->rand_seed       = 1;
   ga_info->datatype        = DT_INT_PERM;
   ga_info->ip_flag         = IP_RANDOM;
   ga_info->ip_data[0]      = '\0';
   ga_info->chrom_len       = 10;
   ga_info->pool_size       = 100;
   ga_info->iter            = -1;
   ga_info->max_iter        = -1;
   ga_info->bias            = 1.8;
   ga_info->gap             = 0.0;
   ga_info->x_rate          = 1.0;
   ga_info->mu_rate         = 0.0;
   ga_info->scale_factor    = 0.0;
   ga_info->minimize        = TRUE;
   ga_info->elitist         = TRUE;
   ga_info->converged       = FALSE;
   ga_info->use_convergence = TRUE;

   /*--- Default operators ---*/
   SE_select(ga_info, "roulette");
    X_select(ga_info, "order1");
   MU_select(ga_info, "swap");
   RE_select(ga_info, "append");
   GA_select(ga_info, "generational");
   ga_info->EV_fun = NULL;

   /*--- Default report parameters ---*/
   ga_info->rp_type      = RP_SHORT;
   ga_info->rp_interval  = 1;
   ga_info->rp_fid       = stdout;
   ga_info->rp_file[0]   = '\0';

   /*--- Reset pools ---*/
   if(PL_valid(ga_info->old_pool)) PL_reset(ga_info->old_pool);
   if(PL_valid(ga_info->new_pool)) PL_reset(ga_info->new_pool);

   /*--- Reset best ---*/
   if(CH_valid(ga_info->best)) CH_reset(ga_info->best);
}

/*----------------------------------------------------------------------------
| Print out all of the config information
----------------------------------------------------------------------------*/
CF_report(
   GA_Info_Ptr    ga_info)
{
   FILE *fid;
   char *sptr;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("CF_report: invalid ga_info");
   if(ga_info->rp_fid == NULL) UT_error("CF_report: invalid ga_info->rp_fid");

   fid = ga_info->rp_fid;

   /*--- Header ---*/
   fprintf(fid,"\nLibGA Version %s\n%s\n%s\n\n", VERSION, COPYRIGHT, COPYRIGHT2);
   fprintf(fid,"GA Configuration Information:\n");
   fprintf(fid,"-----------------------------\n");

   /*--- Basic info --- */
   fprintf(fid,"Basic Info\n");
   if(ga_info->user_data[0] != '\0')
      fprintf(fid,"   User Data         : %s\n", ga_info->user_data);
      fprintf(fid,"   Fuction Index     : %d\n", ga_info->function_index);
   fprintf(fid,"   Random Seed       : %d\n", ga_info->rand_seed);
   fprintf(fid,"   Data Type         : ");
   switch(ga_info->datatype) {
      case DT_BIT:      fprintf(fid,"Bit\n"); break;
      case DT_INT:      fprintf(fid,"Integer\n"); break;
      case DT_INT_PERM: fprintf(fid,"Integer Permutation\n"); break;
      case DT_REAL:     fprintf(fid,"Real\n"); break;
      default:          fprintf(fid,"Unspecified\n"); break;
   }
   fprintf(fid,"   Init Pool Entered : ");
   switch(ga_info->ip_flag) {
      case IP_RANDOM     : fprintf(fid,"Randomly     \n"); break;
      case IP_FROM_FILE  : fprintf(fid,"From File    \n"); break;
      case IP_INTERACTIVE: fprintf(fid,"Interactively\n"); break;
      default            : fprintf(fid,"Unspecified  \n"); break;
   }
   if(ga_info->ip_flag == IP_FROM_FILE) {
      fprintf(fid,"   Initial Pool File : ");
      if(ga_info->ip_data[0] == 0)
         fprintf(fid,"None\n");
      else if(!strcmp(ga_info->ip_data, "UNSPECIFIED"))
         fprintf(fid,"Unspecified\n");
      else
         fprintf(fid,"%s\n", ga_info->ip_data);
   }
   fprintf(fid,"   Chromosome Length : %d\n", ga_info->chrom_len);
   fprintf(fid,"   Pool Size         : %d\n", ga_info->pool_size);
   fprintf(fid,"   Number of Trials  : ");
   if(ga_info->max_iter < 0)
      fprintf(fid,"Run until convergence\n");
   else
      fprintf(fid,"%d iterations, %s\n", 
              ga_info->max_iter,
              ga_info->use_convergence ? "or until convergence" 
                                       : "ignore convergence"
      );
   fprintf(fid,"   Minimize          : %s\n", 
      ga_info->minimize ? "Yes" : "No");
   fprintf(fid,"   Elitism           : %s\n", 
      ga_info->elitist ? "Yes" : "No");
   fprintf(fid,"   Scale Factor      : %G\n", ga_info->scale_factor);

   /*--- Functions ---*/
   fprintf(fid,"\n");
   fprintf(fid,"Functions\n");
   fprintf(fid,"   GA          : %s (Gap = %G)\n", 
      GA_name(ga_info), ga_info->gap);
   fprintf(fid,"   Selection   : %s ", sptr = SE_name(ga_info));
   if(!strcmp(sptr,"rank_biased")) fprintf(fid,"(Bias = %G)", ga_info->bias);
   fprintf(fid,"\n");
   fprintf(fid,"   Crossover   : %s (Rate = %G)\n", 
      X_name(ga_info), ga_info->x_rate);
   if(ga_info->mu_rate > 0.0)
      fprintf(fid,"   Mutation    : %s (Rate = %G)\n", 
         MU_name(ga_info), ga_info->mu_rate);
   fprintf(fid,"   Replacement : %s\n", RE_name(ga_info));

   /*--- Reports ---*/
   if(ga_info->rp_type != RP_NONE) {
      fprintf(fid,"\n");
      fprintf(fid,"Reports\n"); 
      fprintf(fid,"   Type     : %s\n", 
         ga_info->rp_type == RP_MINIMAL ? "Minimal" :
         ga_info->rp_type == RP_SHORT   ? "Short"   :
         ga_info->rp_type == RP_LONG    ? "Long"    :
         "Unknown");
      fprintf(fid,"   Interval : %d\n", ga_info->rp_interval);
      if(ga_info->rp_file[0] != 0) {
         fprintf(fid,"   File  : ");
         if(!strcmp(ga_info->rp_file, "UNSPECIFIED"))
            fprintf(fid,"Unspecified\n");
         else
            fprintf(fid,"%s\n", ga_info->rp_file);
      }
   }

   /*--- Footer ---*/
   fprintf(fid,"-----------------------------\n");

   /*--- Make sure it is printed immediately ---*/
   fflush(fid);
}

/*----------------------------------------------------------------------------
| Read config file
----------------------------------------------------------------------------*/
CF_read(
   GA_Info_Ptr ga_info,
   char        *cfg_name)
{
   char str[STRLEN], token[MAXTOK][STRLEN];
   int  numtok, CF_tokenize();
   FILE *fid;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("CF_read: invalid ga_info");

   /*--- Try to open config file ---*/
   if((fid = fopen(cfg_name, "r")) == NULL) {
      UT_warn("CF_read: Error opening config file.");
      return GA_ERROR;
   }

   /*--- Read each line from config file ---*/
   while(fgets(str, STRLEN, fid) != NULL) {

      /*--- Convert to tokens ---*/
      if((numtok=CF_tokenize(str, token)) <= 0) continue;

      /*--- Which command? ---*/
      switch(token[0][0]) {

      case 'b': 
         if(!strcmp(token[0], "bias")) {
            if(numtok >= 2 && sscanf(token[1], "%f", &ga_info->bias) == 1)
               ;
            else
               UT_warn("CF_read: Invalid bias response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'c': 
         if(!strcmp(token[0], "chrom_len")) {
            if(numtok >= 2 && sscanf(token[1], "%d", &ga_info->chrom_len) == 1)
               ;
            else
               UT_warn("CF_read: Invalid chrom_len response");
         } else if(!strcmp(token[0], "crossover")) {
            if(numtok >= 2)
               X_select(ga_info, token[1]);
            else
               UT_warn("CF_read: Invalid crossover response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'd': 
         if(!strcmp(token[0], "datatype")) {
            if(numtok >= 2 && !strcmp(token[1], "bit")) 
               ga_info->datatype = DT_BIT;
            else if(numtok >= 2 && !strcmp(token[1], "int")) 
               ga_info->datatype = DT_INT;
            else if(numtok >= 2 && !strcmp(token[1], "int_perm")) 
               ga_info->datatype = DT_INT_PERM;
            else if(numtok >= 2 && !strcmp(token[1], "real")) 
               ga_info->datatype = DT_REAL;
            else
               UT_warn("CF_read: Invalid datatype response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'e': 
         if(!strcmp(token[0], "elitism")) {
            if(numtok >= 2 && !strcmp(token[1], "true"))
               ga_info->elitist = TRUE;
            else if(numtok >= 2 && !strcmp(token[1], "false"))
               ga_info->elitist = FALSE;
            else
               UT_warn("CF_read: Invalid elitism response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'f': 
         if(!strcmp(token[0], "function_index")) {
            if(numtok >= 2)
               sscanf(token[1], "%d", &ga_info->function_index);
            else
               UT_warn("CF_read: Invalid function_index response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'g': 
         if(!strcmp(token[0], "gap")) {
            if(numtok >= 2 && sscanf(token[1], "%f", &ga_info->gap) == 1)
               ;
            else
               UT_warn("CF_read: Invalid gap response");
         } else if(!strcmp(token[0], "ga")) {
            if(numtok >= 2) {
               GA_select(ga_info, token[1]);
               if(!strcmp(token[1], "generational")) {
                  SE_select(ga_info, "roulette");
                  RE_select(ga_info, "append");
                  ga_info->rp_interval = 1;
               } else if(!strcmp(token[1], "steady_state")) {
                  SE_select(ga_info, "rank_biased");
                  RE_select(ga_info, "by_rank");
                  ga_info->rp_interval = 100;
               }
            } else
               UT_warn("CF_read: Invalid ga response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'i': 
         if(!strcmp(token[0], "initpool")) {
            if(numtok >= 2 && !strcmp(token[1], "random")) 
               ga_info->ip_flag = IP_RANDOM;
            if(numtok >= 2 && !strcmp(token[1], "random01")) 
	      ga_info->ip_flag = IP_RANDOM01;
            else if(numtok >= 2 && !strcmp(token[1], "from_file")) {
               ga_info->ip_flag = IP_FROM_FILE;
               if(numtok >= 3) strcpy(ga_info->ip_data, token[2]);
            } else if(numtok >= 2 && !strcmp(token[1], "interactive")) 
               ga_info->ip_flag = IP_INTERACTIVE;
            else
               UT_warn("CF_read: Invalid initpool response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'm': 
         if(!strcmp(token[0], "mutation")) {
            if(numtok >= 2) 
               MU_select(ga_info, token[1]);
            else
               UT_warn("CF_read: Invalid mutation response");
         } else if(!strcmp(token[0], "mu_rate")) {
            if(numtok >= 2)
               sscanf(token[1], "%f", &ga_info->mu_rate);
            else
               UT_warn("CF_read: Invalid mu_rate response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'o': 
         if(!strcmp(token[0], "objective")) {
            if(numtok >= 2 && !strcmp(token[1], "minimize"))
               ga_info->minimize = TRUE;
            else if(numtok >= 2 && !strcmp(token[1], "maximize"))
               ga_info->minimize = FALSE;
            else
               UT_warn("CF_read: Invalid objective response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'p': 
         if(!strcmp(token[0], "pool_size")) {
            if(numtok >= 2)
               sscanf(token[1], "%d", &ga_info->pool_size);
            else
               UT_warn("CF_read: Invalid pool_size response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'r': 
         if(!strcmp(token[0], "replacement")) {
            if(numtok >= 2)
               RE_select(ga_info, token[1]);
            else
               UT_warn("CF_read: Invalid replacement response");
         } else if(!strcmp(token[0], "rp_interval")) {
            if(numtok >= 2)
               sscanf(token[1], "%d", &ga_info->rp_interval);
            else
               UT_warn("CF_read: Invalid rp_interval response");
         } else if(!strcmp(token[0], "rp_type")) {
            if(numtok >= 2 && !strcmp(token[1], "minimal"))
               ga_info->rp_type = RP_MINIMAL;
            else if(numtok >= 2 && !strcmp(token[1], "short"))
               ga_info->rp_type = RP_SHORT;
            else if(numtok >= 2 && !strcmp(token[1], "long"))
               ga_info->rp_type = RP_LONG;
            else if(numtok >= 2 && !strcmp(token[1], "none"))
               ga_info->rp_type = RP_NONE;
            else
               UT_warn("CF_read: Invalid rp_type response");
         } else if(!strcmp(token[0], "rp_file")) {
            if(numtok >= 2) {
               char *file_mode = "a";

               /*--- Save file name ---*/
               strcpy(ga_info->rp_file, token[1]);

               /*--- Get file mode if provided ---*/
               if(numtok >= 3)
                  file_mode = token[2];

               /*--- Open report file ---*/
               ga_info->rp_fid = fopen(ga_info->rp_file, file_mode);
               if(ga_info->rp_fid == NULL)
                  UT_error("CF_read: error opening report file");
            } else
               UT_warn("CF_read: Invalid rp_file response");
         } else if(!strcmp(token[0], "rand_seed")) {
            if(numtok >= 2 && !strcmp(token[1], "my_pid"))
               ga_info->rand_seed = getpid();
            else if(numtok >= 2)
               sscanf(token[1], "%d", &ga_info->rand_seed);
            else
               UT_warn("CF_read: Invalid rand_seed response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 's': 
         if(!strcmp(token[0], "selection")) {
            if(numtok >= 2)
               SE_select(ga_info, token[1]);
            else
               UT_warn("CF_read: Invalid selection response");
         } else if(!strcmp(token[0], "stop_after")) {
            if(numtok == 2 && !strcmp(token[1], "convergence")) {
               ga_info->use_convergence = TRUE;
               ga_info->max_iter = -1;
            } else if(numtok >= 2) {
               sscanf(token[1], "%d", &ga_info->max_iter);
               if(ga_info->max_iter < 1) 
                  UT_warn("CF_read: Invalid number for stop_after");
               ga_info->use_convergence = TRUE;
               if(numtok > 2 && !strcmp(token[2], "ignore_convergence"))
                  ga_info->use_convergence = FALSE;
            } else {
               UT_warn("CF_read: Invalid stop_after response");
            }
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'u': 
         if(!strcmp(token[0], "user_data")) {
            if(numtok >= 2)
               strcpy(ga_info->user_data, token[1]);
            else
               UT_warn("CF_read: Invalid user_data response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      case 'x': 
         if(!strcmp(token[0], "x_rate")) {
            if(numtok >= 2 && sscanf(token[1], "%f", &ga_info->x_rate) == 1)
               ;
            else
               UT_warn("CF_read: Invalid x_rate response");
         } else
            UT_warn("CF_read: Unknown config command");
         break;

      default:
            UT_warn("CF_read: Unknown config command");
      }
   }
   /*--- PATCH 1 BEGIN ---*/
   /* Many thanks to Paul-Erik Raue (peraue@cs.vu.nl) 
    * for finding this bug. 
    * 
    * Close the configuration file 
    */
   fclose(fid);
   /*--- PATCH 1 END ---*/
}

/*----------------------------------------------------------------------------
| Convert input line to tokens
----------------------------------------------------------------------------*/
CF_tokenize(
   char *line, char token[MAXTOK][STRLEN])
{
   int numtok, i, j, len;

   numtok = 0;
   len = strlen(line);
   for(i = 0; i < len; ) {

      /*--- Find token ---*/
      while(isspace(line[i]) && i < len) i++;

      /*--- Skip comments & blank lines ---*/
      if(i >= len || line[i] == '#' || line[i] == '\n') break;

      /*--- Save token ---*/
      for(j = 0; !isspace(line[i]) && i < len; i++, j++) 
         token[numtok][j] = line[i];
      token[numtok++][j] = 0;
   }

   return numtok;
}

/*----------------------------------------------------------------------------
| Verify configuration
----------------------------------------------------------------------------*/
CF_verify(
   GA_Info_Ptr ga_info)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("CF_verify: invalid ga_info");

   switch(ga_info->datatype) {
      case DT_BIT:
      case DT_INT:
      case DT_INT_PERM:
      case DT_REAL:
         break;
      default: UT_error("CF_verify: Invalid datatype");
   }

   switch(ga_info->ip_flag) {
      case IP_FROM_FILE:
         if(ga_info->ip_data[0] == '\0')
            UT_error("CF_verify: no file specified for initpool");
         break;
      case IP_INTERACTIVE:
      case IP_RANDOM:
      case IP_NONE:
	  case IP_RANDOM01:
      
         break;
      default: UT_error("CF_verify: Invalid ip_flag");
   }

   if(ga_info->chrom_len <= 0)
      UT_error("CF_verify: invalid chromosome length");

   if(ga_info->pool_size <= 0)
      UT_error("CF_verify: invalid pool size");

   if(ga_info->SE_fun == NULL)
      UT_error("CF_verify: no selection function specified");

   if(ga_info->X_fun == NULL)
      UT_error("CF_verify: no crossover function specified");

   if(ga_info->x_rate < 0.0 || ga_info->x_rate > 1.0)
      UT_error("CF_verify: invalid crossover rate");

   if(ga_info->mu_rate < 0.0)
      UT_error("CF_verify: invalid mutation rate");

   if(ga_info->mu_rate > 0.0 && ga_info->MU_fun == NULL)
      UT_error("CF_verify: no mutation function specified");

   if(ga_info->RE_fun == NULL)
      UT_error("CF_verify: no replacement function specified");

   if(ga_info->EV_fun == NULL)
      UT_error("CF_verify: no evaluation function specified");

   if(ga_info->GA_fun == NULL)
      UT_error("CF_verify: no ga function specified");

   if(ga_info->gap < 0.0 || ga_info->gap > 1.0)
      UT_error("CF_verify: invalid generation gap");

   if(ga_info->minimize != TRUE && ga_info->minimize != FALSE)
      UT_error("CF_verify: illegal value for minimize");

   if(ga_info->elitist != TRUE && ga_info->elitist != FALSE)
      UT_error("CF_verify: illegal value for elitism");

   switch(ga_info->rp_type) {
      case RP_NONE:
      case RP_MINIMAL:
      case RP_SHORT:
      case RP_LONG:
         break;
      default: UT_error("CF_verify: Invalid report type");
   }

   if(ga_info->rp_interval <= 0)
      UT_error("CF_verify: invalid report interval");
}
/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| Crossover operators 
|
| Bit String Representations
|    X_simple()   - simple crossover
|    X_uniform()  - uniform crossover
|
| Order-Based Integer Representations
|    X_order1()     - order1   (Starkweather, et. al., 1991 GA Conf.)
|    X_order2()     - order2   (Starkweather, et. al., 1991 GA Conf.)
|    X_pos()        - position (Starkweather, et. al., 1991 GA Conf.)
|    X_cycle()      - cycle    (Starkweather, et. al., 1991 GA Conf.)
|    X_pmx()        - PMX      (Starkweather, et. al., 1991 GA Conf.)
|    X_uox()        - uniform order crossover
|    X_rox()        - relative order crossover
|    X_asex()       - asexual crossover
|       X_do_asex() - helper for X_asex()
|
| Interface
|    X_table[]    - used in selection of crossover method
|    X_set_fun()  - set and select user defined crossover function
|    X_select()   - select crossover method by name
|    X_name()     - get name of current crossover function
|    X_fun()      - setup and perform current crossover operator
|
| Utility
|    X_gen_xp()    - generate a random crossover point
|    X_gen_2_xp()  - generate two sorted, random crossover points
|    X_gen_4_xp()  - generate four sorted, random crossover points
|    X_init_kids() - reset children for crossover
|    X_map()       - find allele in a chromosome
|
| NOTE: Crossover points should always be thought of as "inclusive"
============================================================================*/

/*============================================================================
|                           Crossover interface
============================================================================*/
/*----------------------------------------------------------------------------
| Table for selection of crossover method
----------------------------------------------------------------------------*/
FN_Table_Type X_table[] = {
   {  NULL,       NULL      }, /* user defined function */
   { "simple",    X_simple  },
   { "uniform",   X_uniform },
   { "order1",    X_order1  },
   { "order2",    X_order2  },
   { "position",  X_pos,    },
   { "cycle",     X_cycle,  },
   { "pmx",       X_pmx,    },
   { "uox",       X_uox,    },
   { "rox",       X_rox,    },
   { "asexual",   X_asex,   },
   { NULL,        NULL,     }
};

/*----------------------------------------------------------------------------
| Select user crossover method
----------------------------------------------------------------------------*/
X_set_fun(
   GA_Info_Ptr  ga_info,
   char         *fn_name,
   FN_Ptr       fn_ptr)
{
   return FN_set_fun(ga_info, X_table, fn_name, fn_ptr, &ga_info->X_fun);
}

/*----------------------------------------------------------------------------
| Select crossover method
----------------------------------------------------------------------------*/
X_select(
   GA_Info_Ptr  ga_info,
   char         *fn_name)
{
   return FN_select(ga_info, X_table, fn_name, &ga_info->X_fun);
}

/*----------------------------------------------------------------------------
| Crossover name
----------------------------------------------------------------------------*/
char *X_name(
   GA_Info_Ptr  ga_info)
{
   return FN_name(ga_info, X_table, ga_info->X_fun);
}

/*----------------------------------------------------------------------------
| Crossover interface
----------------------------------------------------------------------------*/
X_fun(
   GA_Info_Ptr ga_info,
   Chrom_Ptr   parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr   child_1,Chrom_Ptr child_2)
{
   /*--- Init children ---*/
   X_init_kids(parent_1, parent_2, child_1, child_2);

   /*--- Clone instead of crossover ---*/
   if(ga_info->x_rate < 1.0 && RAND_FRAC() > ga_info->x_rate) {
      CH_copy(parent_1, child_1);
      CH_copy(parent_2, child_2);
      child_1->parent_1 = parent_1->index;
      child_1->parent_2 = parent_2->index;
      child_2->parent_1 = parent_1->index;
      child_2->parent_2 = parent_2->index;
      return OK;
   }

   /*--- No crossover function ---*/
   if(ga_info->X_fun == NULL) 
      return GA_ERROR;

   /*--- Crossover ---*/
   ga_info->X_fun(ga_info, parent_1, parent_2, child_1, child_2);
}

/*============================================================================
|                           Crossover operators
============================================================================*/
/*----------------------------------------------------------------------------
| Simple crossover
|
| A single crossover point is selected at random.  Alleles up to and including
| the crossover point are copied to the respective child.  The remaining
| alleles are copied to the alternate child.
----------------------------------------------------------------------------*/
X_simple(
   GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   unsigned i, xp;
   Gene_Type tmp;

   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype == DT_INT_PERM)
      UT_error("X_simple: bad data type");

   /*--- Cannot yet deal with heterozygous parents ---*/
   if(parent_1->length != parent_2->length)
      UT_error("crossover: heterozygous parents");

   /*--- Random crossover point ---*/
   X_gen_xp(0, parent_1->length-1, &xp);
   child_1->xp1 = xp;
   child_2->xp1 = xp;

   /*--- Half is same as parent ---*/
   for(i = 0; i <= xp; i++) {
      child_1->gene[i] = parent_1->gene[i];
      child_2->gene[i] = parent_2->gene[i];
   }

   /*--- Other half is swapped ---*/
   for(i = xp+1; i < parent_1->length; i++) {
      child_1->gene[i] = parent_2->gene[i];
      child_2->gene[i] = parent_1->gene[i];
   }

   return OK;
}

/*----------------------------------------------------------------------------
| Uniform crossover
|
| Each allele is copied from a parent based on a random flip of a fair coin
----------------------------------------------------------------------------*/
X_uniform(GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   unsigned i;

   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype == DT_INT_PERM)
      UT_error("X_uniform: bad data type");

   /*--- Cannot yet deal with heterozygous parents ---*/
   if(parent_1->length != parent_2->length)
      UT_error("crossover: heterozygous parents");

   for(i = 0; i < parent_1->length; i++) {
      if(RAND_BIT()) {
         child_1->gene[i] = parent_1->gene[i];
         child_2->gene[i] = parent_2->gene[i];
      } else {
         child_1->gene[i] = parent_2->gene[i];
         child_2->gene[i] = parent_1->gene[i];
      }
   }

   return OK;
}

/*----------------------------------------------------------------------------
| Order1 (Starkweather, et. al., 1991 GA Conf.)
|
| "The offspring inherits the elements between the two crossover points,
| inclusive, from the selected parent in the same order and position as they
| appeared in that parent.  The remaining elements are inherited from the
| alternate parent in the order in which they appear in that parent, beginning
| with the first position following the second crossover point and skipping
| over all elements already present in the offspring."
|
| L. Davis. "Applying Adaptive Algorithms to Epistatic Domains". In Proc.
|    International Joint Conference on Artificial Intelligence (1985).
----------------------------------------------------------------------------*/
X_order1(GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   int xp1, xp2;
   int i, p1, p2, c;                         

   
   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype != DT_INT_PERM)
      UT_error("X_order1: bad data type");

   /*--- Cannot yet deal with heterozygous parents ---*/
   if(parent_1->length != parent_2->length)
      UT_error("crossover: heterozygous parents");

   /*--- Select two sorted crossover points ---*/
   X_gen_2_xp(FALSE, 0, parent_1->length, &xp1, &xp2);
   child_1->xp1 = child_2->xp1 = xp1;
   child_1->xp2 = child_2->xp2 = xp2;

   /*--- Info between xp is same as parent ---*/
   for(i = xp1; i <= xp2; i++) {
      child_1->gene[i] = parent_1->gene[i];
      child_2->gene[i] = parent_2->gene[i];
   }

   /*--- Inherit remainder from other parent ---*/
   for (i=0, p1=p2=xp2; i < (parent_1->length - (xp2 - xp1 + 1)); i++) {

      /*--- Index to fill in children ---*/
      c = (xp2 + 1 + i) % parent_1->length;

      /*--- Child 1 gets next unused element in parent 2 ---*/
      while(TRUE) {
         p2 = (p2 + 1) % parent_1->length;
         if(X_map(&parent_2->gene[p2], parent_1, xp1, xp2) < 0) break;
      }

      /*--- Child 2 gets next unused element in parent 1 ---*/
      while(TRUE) {
         p1 = (p1 + 1) % parent_2->length;
         if(X_map(&parent_1->gene[p1], parent_2, xp1, xp2) < 0) break;
      }

      /*--- Transfer to children ---*/
      child_1->gene[c] = parent_2->gene[p2];
      child_2->gene[c] = parent_1->gene[p1];
   }

   return OK; 
}

/*----------------------------------------------------------------------------
| Order2 (Starkweather, et. al., 1991 GA Conf.)
|
| "...several key positions are chosen randomly and the order in which these
| elements appear in one parent is imposed on the other parent to produce
| two offspring..."
|
| G. Syswerda.  "Schedule Optimization Using Genetic Algorithms".  In Handbook
|    of Genetic Algorithms.  L. Davis ed.  Van Nostrand Reinhold: NY (1990).
----------------------------------------------------------------------------*/
X_order2(GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   int xp1, xp2, xp3, xp4;
   int i, j1, j2;                         
   int xidx_1[4], xidx_2[4];

   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype != DT_INT_PERM)
      UT_error("X_order2: bad data type");

   /*--- Cannot yet deal with heterozygous parents ---*/
   if(parent_1->length != parent_2->length)
      UT_error("crossover: heterozygous parents");

   /*--- Select four sorted crossover points ---*/
   X_gen_4_xp(TRUE, 0, parent_1->length, &xp1, &xp2, &xp3, &xp4);
   child_1->xp1 = xp1; child_1->xp2 = xp2;
   child_2->xp1 = xp3; child_2->xp2 = xp4;

   /*--- Children look like parents ---*/
   for (i=0; i < parent_1->length; i++) {
      child_1->gene[i] = parent_1->gene[i];
      child_2->gene[i] = parent_2->gene[i];
   }

   /*--- Map order of xp's in other parent ---*/
   for (i=0, j1=j2=0; i < parent_1->length; i++) {

      /*--- Child_1 uses order in parent_2 ---*/
      if( (int)parent_2->gene[i] == (int)parent_1->gene[xp1] || 
          (int)parent_2->gene[i] == (int)parent_1->gene[xp2] ||
          (int)parent_2->gene[i] == (int)parent_1->gene[xp3] ||
          (int)parent_2->gene[i] == (int)parent_1->gene[xp4]     ) {
         xidx_1[j1++] = i;
      }

      /*--- Child_2 uses order in parent_1 ---*/
      if( (int)parent_1->gene[i] == (int)parent_2->gene[xp1] || 
          (int)parent_1->gene[i] == (int)parent_2->gene[xp2] ||
          (int)parent_1->gene[i] == (int)parent_2->gene[xp3] ||
          (int)parent_1->gene[i] == (int)parent_2->gene[xp4]     ) {
         xidx_2[j2++] = i;
      }
   }

   /*--- Impose ordering of xp's from other parent ---*/
   child_1->gene[xp1] = parent_2->gene[xidx_1[0]];
   child_1->gene[xp2] = parent_2->gene[xidx_1[1]];
   child_1->gene[xp3] = parent_2->gene[xidx_1[2]];
   child_1->gene[xp4] = parent_2->gene[xidx_1[3]];
   child_2->gene[xp1] = parent_1->gene[xidx_2[0]];
   child_2->gene[xp2] = parent_1->gene[xidx_2[1]];
   child_2->gene[xp3] = parent_1->gene[xidx_2[2]];
   child_2->gene[xp4] = parent_1->gene[xidx_2[3]];

   return OK; 
}

/*----------------------------------------------------------------------------
| Position (Starkweather, et. al., 1991 GA Conf.)
|
| "Several random locations in the sequence are selected along with one 
| parent; the elements in those positions are inherited from that parent.
| The remaining elements are inherited in the order in which they appear in
| the alternate parent, skipping over all elements which have already been
| included in the offspring."
|
| G. Syswerda.  "Schedule Optimization Using Genetic Algorithms".  In Handbook
|    of Genetic Algorithms.  L. Davis ed.  Van Nostrand Reinhold: NY (1990).
----------------------------------------------------------------------------*/
X_pos(GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   int xp1, xp2, xp3, xp4;
   int i, j1, j2;                         

   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype != DT_INT_PERM)
      UT_error("X_pos: bad data type");

   /*--- Cannot yet deal with heterozygous parents ---*/
   if(parent_1->length != parent_2->length)
      UT_error("crossover: heterozygous parents");

   /*--- Select four sorted crossover points ---*/
   X_gen_4_xp(FALSE, 0, parent_1->length, &xp1, &xp2, &xp3, &xp4);
   child_1->xp1 = xp1; child_1->xp2 = xp2;
   child_2->xp1 = xp3; child_2->xp2 = xp4;

   /*--- Children get parent's xp values ---*/
   child_1->gene[xp1] = parent_1->gene[xp1];
   child_1->gene[xp2] = parent_1->gene[xp2];
   child_1->gene[xp3] = parent_1->gene[xp3];
   child_1->gene[xp4] = parent_1->gene[xp4];
   child_2->gene[xp1] = parent_2->gene[xp1];
   child_2->gene[xp2] = parent_2->gene[xp2];
   child_2->gene[xp3] = parent_2->gene[xp3];
   child_2->gene[xp4] = parent_2->gene[xp4];

   /*--- Inherit rest using order from other parent ---*/
   for (i=0, j1=j2=0; i < parent_1->length; i++) {

      /*--- Transfer if not a crossover point (child_1) ---*/
      if( (int)parent_2->gene[i] != (int)parent_1->gene[xp1] && 
          (int)parent_2->gene[i] != (int)parent_1->gene[xp2] &&
          (int)parent_2->gene[i] != (int)parent_1->gene[xp3] && 
          (int)parent_2->gene[i] != (int)parent_1->gene[xp4]   ) {

         /*--- Make sure j1 is not a crossover point ---*/
         while(j1 == xp1 || j1 == xp2 || j1 == xp3 || j1 == xp4) j1++;

         child_1->gene[j1++] = parent_2->gene[i];
      }

      /*--- Transfer if not a crossover point (child_2) ---*/
      if( (int)parent_1->gene[i] != (int)parent_2->gene[xp1] && 
          (int)parent_1->gene[i] != (int)parent_2->gene[xp2] &&
          (int)parent_1->gene[i] != (int)parent_2->gene[xp3] && 
          (int)parent_1->gene[i] != (int)parent_2->gene[xp4]   ) {

         /*--- Make sure j2 is not a crossover point ---*/
         while(j2 == xp1 || j2 == xp2 || j2 == xp3 || j2 == xp4) j2++;

         child_2->gene[j2++] = parent_1->gene[i];
      }
   }

   return OK; 
}

/*----------------------------------------------------------------------------
| Cycle (Starkweather, et. al., 1991 GA Conf.)
|
| "A parent sequence and a cycle starting point are randomly selected.  The 
| element at the cycle starting point of the selected parent is inherited by
| the child.  The element which is in the same position in the other parent
| cannot then be placed in this position so its position is found in the
| selected parent and is inherited from that position by the child.  This
| continues until the cycle is completed by encountering the initial item in
| the unselected parent.  Any elements which are not yet present in the
| offspring are inherited from the unselected parent."
|
| I. Oliver, D. Smith, and J. Holland.  "A Study of Permutation Crossover
|   Operators on the Traveling Salesman Problem."  In Proc. Second Interna-
|   tional Conference on Genetic Algorithms and their Applications (1987).
----------------------------------------------------------------------------*/
X_cycle(GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   int xp, i;

   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype != DT_INT_PERM)
      UT_error("X_cycle: bad data type");

   /*--- Cannot yet deal with heterozygous parents ---*/
   if(parent_1->length != parent_2->length)
      UT_error("crossover: heterozygous parents");

   /*--- Select crossover point ---*/
   X_gen_xp(0, parent_1->length, &xp);
   child_1->xp1 = xp;
   child_2->xp1 = xp;

   /*--- Transfer material to children ---*/
   for(i = 0; i < parent_1->length; i++) {
      child_1->gene[i] = parent_2->gene[i];
      child_2->gene[i] = parent_1->gene[i];
   }

   /*--- Crossover (child 1) ---*/
   for (i=xp; ; ) {
      child_1->gene[i] = parent_1->gene[i];
      i = X_map(&parent_2->gene[i], parent_1, 0, parent_1->length - 1);
      if(i == xp) break;
   }

   /*--- Crossover (child 2) ---*/
   for (i=xp; ; ) {
      child_2->gene[i] = parent_2->gene[i];
      i = X_map(&parent_1->gene[i], parent_2, 0, parent_2->length - 1);
      if(i == xp) break;
   }

   return OK; 
}

/*----------------------------------------------------------------------------
| PMX (Starkweather, et. al., 1991 GA Conf.)
|
| "A parent and two crossover sites are selected randomly and the elements
| between the two starting positions in one of the parents are directly 
| inherited by the offspring.  Each element between the two crossover points
| in the alternate parent are mapped to the position held by this element in
| the first parent.  Then the remaining elements are inherited from the
| alternate parent.
|
| D. Goldberg and R. Lingle.  "Alleles, loci, and the Traveling Salesman
| Problem".  In Proc. International Conference on Genetic Algorithms and their
| Applications (1985).
----------------------------------------------------------------------------*/
X_pmx(GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   int xp1, xp2;
   int i, j, k;

   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype != DT_INT_PERM)
      UT_error("X_pmx: bad data type");

   /*--- Cannot yet deal with heterozygous parents ---*/
   if(parent_1->length != parent_2->length)
      UT_error("crossover: heterozygous parents");

   /*--- Select two sorted crossover points ---*/
   X_gen_2_xp(FALSE, 0, parent_1->length, &xp1, &xp2);
   child_1->xp1 = child_2->xp1 = xp1; 
   child_1->xp2 = child_2->xp2 = xp2;

   /*--- Copy info to children ---*/
   for(i = 0; i < parent_1->length; i++) {
      if(i < xp1 || i > xp2) {
         child_1->gene[i] = parent_1->gene[i];
         child_2->gene[i] = parent_2->gene[i];
      } else {
         child_1->gene[i] = parent_2->gene[i];
         child_2->gene[i] = parent_1->gene[i];
      }
   }

   /*--- Fixup mapped elements ---*/
   for(i = 0; i < parent_1->length; i++) {

      /*--- Skip if between xp's ---*/
      if(i >= xp1 && i <= xp2) continue;

      /*--- A mapped element (child_1) ---*/
      if((j = X_map(&child_1->gene[i], child_1, xp1, xp2)) >= 0) {
         while(TRUE) {
            child_1->gene[i] = parent_1->gene[j];
            if((j = X_map(&child_1->gene[i], child_1, xp1, xp2)) < 0) 
               break;
         }
      }

      /*--- A mapped element (child_2) ---*/
      if((j = X_map(&child_2->gene[i], child_2, xp1, xp2)) >= 0) {
         while(TRUE) {
            child_2->gene[i] = parent_2->gene[j];
            if((j = X_map(&child_2->gene[i], child_2, xp1, xp2)) < 0) 
               break;
         }
      }
   }

   return OK; 
}

/*----------------------------------------------------------------------------
| Uniform order crossover
|
| Analogous to uniform crossover, but used for order-based representations.
----------------------------------------------------------------------------*/
X_uox(GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   unsigned i, j1, j2;
   int length;
   static int  m1_length, m2_length;
   static char *m1 = NULL, *m2 = NULL;

   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype != DT_INT_PERM)
      UT_error("X_uox: bad data type");

   /*--- Cannot yet deal with heterozygous parents ---*/
   if(parent_1->length != parent_2->length)
      UT_error("crossover: heterozygous parents");

   /*--- Make room for m1 ---*/
   if(m1 == NULL) {
      m1 = (char *)calloc(parent_1->length, sizeof(char));
      if(m1 == NULL) UT_error("X_uox: m1 alloc failed");
      m1_length = parent_1->length;
   } else {
      if(m1_length != parent_1->length) {
         m1 = (char *)realloc(m1, parent_1->length * sizeof(char));
         if(m1 == NULL) UT_error("X_uox: m1 realloc failed");
         m1_length = parent_1->length;
      }
   }

   /*--- Make room for m2 ---*/
   if(m2 == NULL) {
      m2 = (char *)calloc(parent_2->length, sizeof(char));
      if(m2 == NULL) UT_error("X_uox: m2 alloc failed");
      m2_length = parent_2->length;
   } else {
      if(m2_length != parent_2->length) {
         m2 = (char *)realloc(m2, parent_2->length * sizeof(char));
         if(m2 == NULL) UT_error("X_uox: m2 realloc failed");
         m2_length = parent_2->length;
      }
   }

   /*--- Random mask ---*/
   for(i = 0; i < parent_1->length; i++) {
      m1[i] = m2[i] = (RAND_BIT() ? 1 : 0);
   }

   /*--- Place alleles from mask ---*/
   for(i = 0; i < parent_1->length; i++) {
      if(m1[i]) child_1->gene[i] = parent_1->gene[i];
      else      child_1->gene[i] = -1;
   }
   for(i = 0; i < parent_2->length; i++) {
      if(m2[i]) child_2->gene[i] = parent_2->gene[i];
      else      child_2->gene[i] = -1;
   }

   /*--- Place remaining alleles ---*/
   j1 = 0;
   for(i = 0; i < parent_1->length; i++) {
      if((int)child_1->gene[i] == -1) {
         while(X_map(&parent_2->gene[j1], child_1, 0, child_1->length-1) != -1)
            if(j1 < parent_2->length)
               j1++;
            else
               UT_error("X_uox: invalid j1");
         child_1->gene[i] = parent_2->gene[j1];
      }
   }
   j2 = 0;
   for(i = 0; i < parent_2->length; i++) {
      if((int)child_2->gene[i] == -1) {
         while(X_map(&parent_1->gene[j2], child_2, 0, child_2->length-1) != -1)
            if(j2 < parent_1->length)
               j2++;
            else
               UT_error("X_uox: invalid j2");
         child_2->gene[i] = parent_1->gene[j2];
      }
   }
}

/*----------------------------------------------------------------------------
| Relative order crossover
|
| Not yet implemented
----------------------------------------------------------------------------*/
X_rox(GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   UT_error("X_rox: not yet implemented");

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("X_rox: invalid ga_info");
   if(!CH_valid(parent_1)) UT_error("X_rox: invalid parent_1");
   if(!CH_valid(parent_2)) UT_error("X_rox: invalid parent_2");
   if(!CH_valid(child_1)) UT_error("X_rox: invalid child_1");
   if(!CH_valid(child_2)) UT_error("X_rox: invalid child_2");

   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype != DT_INT_PERM)
      UT_error("X_rox: bad data type");

   /*--- Cannot yet deal with heterozygous parents ---*/
   if(parent_1->length != parent_2->length)
      UT_error("crossover: heterozygous parents");
}

/*----------------------------------------------------------------------------
| Asexual crossover
|
| Two crossover points are selected at random and the elements at those points
| are swapped resulting in a single child.  To remain compatible with the
| other crossover methods which generate two children from two parents, two
| independent crossovers are performed to generate the two children.
|
| This method is the same as a "two-opt" in simulated annealing.
----------------------------------------------------------------------------*/
X_asex(GA_Info_Ptr ga_info,
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   /*--- Make sure datatype is compatible ---*/
   if(ga_info->datatype != DT_INT_PERM)
      UT_error("X_asex: bad data type");

   /*--- Perform asexual crossover ---*/
   X_do_asex(parent_1, child_1);
   X_do_asex(parent_2, child_2);

   return OK;
}

/*----------------------------------------------------------------------------
| Asexual crossover
|
| Perform asexual crossover for X_asex()
----------------------------------------------------------------------------*/
X_do_asex(
   Chrom_Ptr  parent,
   Chrom_Ptr  child)
{
   int i, xp1, xp2;

   /*--- Copy info to child ---*/
   for(i = 0; i < parent->length; i++) {
      child->gene[i] = parent->gene[i];
   }
   child->idx_min = parent->idx_min;

   /*--- Prevent infinite loop ---*/
   if(parent->idx_min >= parent->length-1) return OK;

   /*--- Select two sorted crossover points ---*/
   X_gen_2_xp(TRUE, parent->idx_min, parent->length, &xp1, &xp2);
   child->xp1 = xp1; 
   child->xp2 = xp2;

   /*--- Crossover just swaps xp's ---*/
   child->gene[xp1] = parent->gene[xp2];
   child->gene[xp2] = parent->gene[xp1];

   return OK;
}

/*============================================================================
|                             Utility functions
============================================================================*/
/*----------------------------------------------------------------------------
| Generate crossover point in [idx_min..idx_max-1]
----------------------------------------------------------------------------*/
X_gen_xp(
   int idx_min,int idx_max,int *xp)
{
   *xp = RAND_DOM(idx_min, idx_max - 1);
}

/*----------------------------------------------------------------------------
| Generate two sorted crossover points
----------------------------------------------------------------------------*/
X_gen_2_xp(int unique,int  idx_min,int  idx_max,int  *xp1,int  *xp2)  
{
   /*--- Generate two points ---*/
   X_gen_xp(idx_min, idx_max, xp1);
   X_gen_xp(idx_min, idx_max, xp2);

   /*--- Make sure unique if specified ---*/
   if(unique) {
      while(*xp2 == *xp1) 
         X_gen_xp(idx_min, idx_max, xp2);
   }

   /*--- Make sure they are sorted ---*/
   if(*xp1 > *xp2) 
      UT_iswap(xp1, xp2);
}

/*----------------------------------------------------------------------------
| Generate four sorted crossover points
----------------------------------------------------------------------------*/
X_gen_4_xp(int unique,int idx_min,int idx_max,int *xp1,int *xp2,int *xp3,int *xp4)   
{
   /*--- Generate four points ---*/
   X_gen_xp(idx_min, idx_max, xp1);
   X_gen_xp(idx_min, idx_max, xp2);
   X_gen_xp(idx_min, idx_max, xp3);
   X_gen_xp(idx_min, idx_max, xp4);

   /*--- Make sure unique if specified ---*/
   if(unique) {
      while(*xp2 == *xp1) 
         X_gen_xp(idx_min, idx_max, xp2);

      while(*xp3 == *xp1 || *xp3 == *xp2) 
         X_gen_xp(idx_min, idx_max, xp3);

      while(*xp4 == *xp1 || *xp4 == *xp2 || *xp4 == *xp3) 
         X_gen_xp(idx_min, idx_max, xp4);
   }

   /*--- Make sure they are sorted (use "sorting network") ---*/
   if (*xp1 > *xp2) UT_iswap(xp1, xp2);
   if (*xp3 > *xp4) UT_iswap(xp3, xp4);
   if (*xp2 > *xp3) UT_iswap(xp2, xp3);
   if (*xp1 > *xp2) UT_iswap(xp1, xp2);
   if (*xp3 > *xp4) UT_iswap(xp3, xp4);
}

/*----------------------------------------------------------------------------
| Initialize child chromosomes for crossover
----------------------------------------------------------------------------*/
X_init_kids(
   Chrom_Ptr  parent_1,Chrom_Ptr parent_2,
   Chrom_Ptr  child_1,Chrom_Ptr child_2)
{
   /*--- Assume for now that parents are homozygous ---*/
   if(parent_1->length <= 0) UT_error("crossover: parent_1->length");
   if(parent_2->length <= 0) UT_error("crossover: parent_2->length");
   if(child_1 == NULL) UT_error("X_init_kids: null child_1");
   if(child_2 == NULL) UT_error("X_init_kids: null child_2");

   /*--- Initialize the children ---*/
   CH_reset(child_1);
   CH_reset(child_2);
   child_1->parent_1 = parent_1->index;
   child_1->parent_2 = parent_2->index;
   child_2->parent_1 = parent_1->index;
   child_2->parent_2 = parent_2->index;
}

/*----------------------------------------------------------------------------
| Allele map in Chrom->gene[lo..hi]
----------------------------------------------------------------------------*/
X_map(
   Gene_Type  *allele,
   Chrom_Ptr  chrom,
   int lo, int hi)
{
   int i;

   /*--- Error check ---*/
   if(lo < 0 || lo > hi || hi >= chrom->length) 
      UT_error("X_map: bad range");

   /*--- Find allele in range of genes ---*/
   for(i = lo; i <= hi; i++) 
      if((int)*allele == (int)chrom->gene[i]) 
         return i;

   /*--- Not found ---*/
   return -1;
}
/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| Function table management
|
| Functions:
|    FN_set_fun() - Set user function
|    FN_select()  - Select function by name
|    FN_name()    - Get function name from function pointer
============================================================================*/


/*----------------------------------------------------------------------------
| Set and select user defined function
----------------------------------------------------------------------------*/
FN_set_fun(
   GA_Info_Ptr  ga_info,
   FN_Table_Ptr fn_table,
   char         *fn_name,
   FN_Ptr       fn_ptr,FN_Ptr *rtn_fun)
{
   int len;

   /*--- Check for invalid ga_info ---*/
   if(!CF_valid(ga_info)) UT_error("FN_set_fun: invalid ga_info");
   if(rtn_fun == NULL) UT_error("FN_set_fun: invalid rtn_fun");

   /*--- Free current function name ---*/
   if(fn_table[0].name != NULL) {
      free(fn_table[0].name);
      fn_table[0].name = NULL;
   }

   /*--- Set function name if provided ---*/
   if(fn_name != NULL && (len = strlen(fn_name) + 1) > 1) {

      /*--- Allocate memory for function name ---*/
      fn_table[0].name = (char *)calloc(len, sizeof(char));
      if(fn_table[0].name == NULL) UT_error("FN_set_fun: alloc failed");

      /*--- Copy the function name ---*/
      strcpy(fn_table[0].name, fn_name);
   }

   /*--- Set user function ---*/
   fn_table[0].fun = fn_ptr;

   /*--- Set return function ---*/
   *rtn_fun = fn_ptr;

}

/*----------------------------------------------------------------------------
| Select function by name
----------------------------------------------------------------------------*/
FN_select(
   GA_Info_Ptr  ga_info,
   FN_Table_Ptr fn_table,
   char         *fn_name,
   FN_Ptr       *rtn_fun)
{
   int i;

   /*--- Check for invalid ga_info ---*/
   if(!CF_valid(ga_info)) UT_error("FN_select: invalid ga_info");
   if(rtn_fun == NULL) UT_error("FN_select: invalid rtn_fun");

   /*--- User defined crossover? ---*/
   if(fn_table[0].name != NULL &&
      !strncmp(fn_name, fn_table[0].name, strlen(fn_table[0].name))) 
   {
      /*--- Null user function ---*/
      if(fn_table[0].fun == NULL) 
         UT_warn("FN_select: User function is NULL");

      /*--- Return pointer to user function ---*/
      *rtn_fun = fn_table[0].fun; 
      return OK;  //CCC 
   }

   /*--- Search fn_table for matching function name ---*/
   for(i = 1; fn_table[i].fun != NULL; i++) {

      /*--- Does name match? ---*/
      if(!strncmp(fn_name, 
                  fn_table[i].name, 
                  MIN(strlen(fn_name),strlen(fn_table[i].name)))) 
      {
         /*--- Return pointer to function ---*/
         *rtn_fun = fn_table[i].fun; 
         return OK;  //CCC 
      }
   }

   /*--- Invalid selection ---*/
   UT_error("FN_select: Invalid selection");
}

/*----------------------------------------------------------------------------
| Function name
----------------------------------------------------------------------------*/
char *FN_name(
   GA_Info_Ptr  ga_info,
   FN_Table_Ptr fn_table,
   FN_Ptr       fn_ptr)
{
   int i;

   /*--- Check for invalid ga_info ---*/
   if(!CF_valid(ga_info)) UT_error("FN_name: invalid ga_info");

   /*--- Search for current function in fn_table ---*/
   for(i = 0; i == 0 || fn_table[i].fun != NULL; i++) {

      /*--- Does this function match? ---*/
      if(fn_ptr == fn_table[i].fun) {

         /*--- Function match, but null name ---*/
         if(fn_table[i].name == NULL) 
            return "Unspecified";

         /*--- Function match, valid name ---*/
         else    
            return fn_table[i].name;
      }
   }

   /*--- Function not found ---*/
   return "Unknown";
}

/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| Genetic algorithm 
|
| Operators
|    GA_generational()  - generational GA
|       GA_gen_init()   - initialize generational GA
|       GA_init_trial() - initialize inner loop for generational GA
|    GA_steady_state()  - steady state GA
|       GA_ss_init()    - initialize steady state GA
|    
| Interface
|    GA_table[]   - used in selection of GA method
|    GA_set_fun() - set and select user defined GA function
|    GA_select()  - select GA function by name
|    GA_name()    - get name of current GA function
|    GA_config()  - configure the GA (do only once for each ga_info)
|    GA_reset()   - reset the GA
|    GA_run()     - setup and perform current GA
|                   (GA_fun would be more consistent but less intuitive)
|    
| Utility
|    GA_trial()      - a single iteration of the inner loop
|    GA_cum()        - see if children are the cumulative/historical best
|    GA_gap()        - handle generation gap
============================================================================*/

/*============================================================================
|                                  Interface
============================================================================*/
/*----------------------------------------------------------------------------
| GA table
----------------------------------------------------------------------------*/
FN_Table_Type GA_table[] = {
   { NULL,           NULL            },  /* user defined function */
   { "generational", GA_generational },
   { "steady_state", GA_steady_state },
   { NULL,           NULL            }
};

/*----------------------------------------------------------------------------
| Select user GA method
----------------------------------------------------------------------------*/
GA_set_fun(
   GA_Info_Ptr  ga_info,
   char         *fn_name,
   FN_Ptr       fn_ptr)
{
  printf("AA\n");
   return FN_set_fun(ga_info, GA_table, fn_name, fn_ptr, &ga_info->GA_fun);
}

/*----------------------------------------------------------------------------
| Select GA method
----------------------------------------------------------------------------*/
GA_select(
   GA_Info_Ptr  ga_info,
   char         *fn_name)
{
   return FN_select(ga_info, GA_table, fn_name, &ga_info->GA_fun);
}

/*----------------------------------------------------------------------------
| GA name
----------------------------------------------------------------------------*/
char *GA_name(
   GA_Info_Ptr  ga_info)
{
   return FN_name(ga_info, GA_table, ga_info->GA_fun);
}

/*----------------------------------------------------------------------------
| Configure the genetic algorithm
----------------------------------------------------------------------------*/
GA_Info_Ptr GA_config(char *cfg_name, int  (*EV_fun)(Chrom_Ptr chrom))
{
   GA_Info_Ptr ga_info;

   /*--- Get memory for ga_info ---*/
   ga_info = CF_alloc();

   /*--- Register user's evaluation function if provided ---*/
   if(EV_fun != NULL)
      ga_info->EV_fun = EV_fun;

   /*--- Read config file if provided ---*/
   if(cfg_name != NULL && cfg_name[0] != '\0' && cfg_name[0] != '\n')
      CF_read(ga_info, cfg_name);

   return ga_info;
}

/*----------------------------------------------------------------------------
| Reset the genetic algorithm
----------------------------------------------------------------------------*/
GA_reset(
   GA_Info_Ptr ga_info,
   char *cfg_name)
{
   int  (*EV_fun)();

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("GA_reset: invalid ga_info");

   /*--- Save EV_fun() ---*/
   EV_fun = ga_info->EV_fun;

   /*--- Reset ga_info ---*/
   CF_reset(ga_info);

   /*--- Restore EV_fun() ---*/
   ga_info->EV_fun = EV_fun;

   /*--- Read config file if provided ---*/
   if(cfg_name != NULL && cfg_name[0] != '\0' && cfg_name[0] != '\n')
      CF_read(ga_info, cfg_name);
}

/*----------------------------------------------------------------------------
| Run the GA
----------------------------------------------------------------------------*/
GA_run(
   GA_Info_Ptr ga_info)
{
   /*--- Ensure valid config information ---*/
   CF_verify(ga_info);

   /*--- Print out config information ---*/
   RP_config(ga_info);

   //printf("seed: %d",ga_info->rand_seed);
   /*--- Seed random number generator ---*/
   SEED_RAND(ga_info->rand_seed);
   
   /*--- Run the GA ---*/
   ga_info->GA_fun(ga_info);

}

/*============================================================================
|                                Generational GA
============================================================================*/
/*----------------------------------------------------------------------------
| Generational GA
----------------------------------------------------------------------------*/
GA_generational(
   GA_Info_Ptr ga_info)
{
   Pool_Ptr      tmp_pool;

   /*--- Initialize ---*/
   GA_gen_init(ga_info);

   /*--- Outer loop is for each generation ---*/
   for(ga_info->iter = 0; 
       ga_info->max_iter < 0 || ga_info->iter < ga_info->max_iter; 
       ga_info->iter++) {

      /*--- Check for convergence ---*/
      if(ga_info->use_convergence && ga_info->converged) break;

      /*--- Setup for new set of trials ---*/
      GA_init_trial(ga_info);

      /*--- Handle generation gap ---*/
      GA_gap(ga_info);

      /*--- Inner loop is for each reproduction ---*/
      for( ; ga_info->new_pool->size < ga_info->old_pool->size; ) {
         GA_trial(ga_info);
      }

      /*--- Print report if appropriate ---*/
      RP_report(ga_info, ga_info->new_pool);

      /*--- Swap old and new pools ---*/
      tmp_pool          = ga_info->old_pool;
      ga_info->old_pool = ga_info->new_pool;
      ga_info->new_pool = tmp_pool;
   }

   /*--- Final report ---*/
   RP_final(ga_info);

   /*--- Free genes for children ---*/
   CH_free(child1);
   CH_free(child2);

   return OK;
}
 
/*----------------------------------------------------------------------------
| Initialize the GA
----------------------------------------------------------------------------*/
void GA_gen_init(
   GA_Info_Ptr ga_info)
{
   Pool_Ptr     old_pool, new_pool;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("GA_generational: invalid ga_info");

   /*--- Make sure pool allocation is ok ---*/
   if(!PL_valid(ga_info->old_pool))
      ga_info->old_pool = PL_alloc(ga_info->pool_size);
   if(!PL_valid(ga_info->new_pool) || ga_info->new_pool == ga_info->old_pool)
      ga_info->new_pool = PL_alloc(ga_info->pool_size);
   old_pool = ga_info->old_pool;
   new_pool = ga_info->new_pool;

   /*--- Minimize or maximize? ---*/
   old_pool->minimize = new_pool->minimize = ga_info->minimize;
 
   /*--- Initialize or reset pools ---*/
   PL_generate(ga_info, old_pool);

   /*--- Don't initialize pool second time ---*/
   if(ga_info->ip_flag != IP_NONE) ga_info->ip_flag = IP_NONE;
 
   /*--- Pool must have even number of chromosomes ---*/
   if((old_pool->size % 2) != 0) {
      /*--- Put in a good chromosome ---*/
      if(ga_info->minimize)
         PL_append(old_pool, old_pool->chrom[old_pool->min_index], TRUE);
      else
         PL_append(old_pool, old_pool->chrom[old_pool->max_index], TRUE);
   }
   new_pool->size = 0;

   /*--- Make sure best is allocated ---*/
   if(!CH_valid(ga_info->best)) 
      ga_info->best = CH_alloc(ga_info->chrom_len);
    
   /*--- Save best member of initial pool ---*/
   if(ga_info->minimize)
      CH_copy(old_pool->chrom[old_pool->min_index], ga_info->best);
   else
      CH_copy(old_pool->chrom[old_pool->max_index], ga_info->best);
 
   /*--- No mutations yet ---*/
   ga_info->num_mut = 0;
   ga_info->tot_mut = 0;
 
   /*--- Initial pool report ---*/
   ga_info->iter = -1;
   RP_report(ga_info, old_pool);

   /*--- Allocate genes for children ---*/
   child1 = CH_alloc(ga_info->chrom_len);
   child2 = CH_alloc(ga_info->chrom_len);
}
 
/*----------------------------------------------------------------------------
| Setup for a new set of trials (Generational GA only)
----------------------------------------------------------------------------*/
GA_init_trial(GA_Info_Ptr ga_info)
{


   /*--- Cleanup the new pool ---*/
   ga_info->new_pool->size = 0;
 
   /*--- Reset number of mutations ---*/
   ga_info->num_mut = 0;
 
   /*--- Not elitist ---*/
   if(!ga_info->elitist) return OK;
 
   /*--- Save best members if Elitist ---*/
   if(ga_info->minimize) {
      PL_append(ga_info->new_pool,
                ga_info->old_pool->chrom[ga_info->old_pool->min_index],
                TRUE);
      PL_append(ga_info->new_pool,
                ga_info->old_pool->chrom[ga_info->old_pool->min_index],
                TRUE);
   } else {
      PL_append(ga_info->new_pool,
                ga_info->old_pool->chrom[ga_info->old_pool->max_index],
                TRUE);
      PL_append(ga_info->new_pool,
                ga_info->old_pool->chrom[ga_info->old_pool->max_index],
                TRUE);
   }

   return OK;
}

/*============================================================================
|                                Steady State GA
============================================================================*/
/*----------------------------------------------------------------------------
| Steady-state GA 
----------------------------------------------------------------------------*/
GA_steady_state(
   GA_Info_Ptr ga_info)
{
   /*--- Initialize ---*/
   GA_ss_init(ga_info);
 
   /*--- Outer loop is for each generation ---*/
   for(ga_info->iter = 0; 
       ga_info->max_iter < 0 || ga_info->iter < ga_info->max_iter; 
       ga_info->iter++) {
 
      /*--- Check convergence (only if no mutation) ---*/
      if(ga_info->use_convergence && ga_info->converged) break;
 
      /*--- "Inner loop" is a single reproduction ---*/
      GA_trial(ga_info);
 
      /*--- Print report if appropriate ---*/
      RP_report(ga_info, ga_info->new_pool);
   }
 
   /*--- Final report ---*/
   RP_final(ga_info);

   /*--- Free genes for children ---*/
   CH_free(child1);
   CH_free(child2);
 
   return OK;
}

/*----------------------------------------------------------------------------
| Initialize GA
----------------------------------------------------------------------------*/
GA_ss_init(
   GA_Info_Ptr ga_info)
{
   Pool_Ptr pool;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("GA_steady_state: invalid ga_info");
 
   /*--- Make sure pool allocation is ok ---*/
   if(PL_valid(ga_info->old_pool)) {
      if(ga_info->new_pool != ga_info->old_pool && PL_valid(ga_info->new_pool))
         PL_free(ga_info->new_pool);
      ga_info->new_pool = ga_info->old_pool;
   } else {
      if(PL_valid(ga_info->new_pool)) 
         ga_info->old_pool = ga_info->new_pool;
      else
         ga_info->old_pool = ga_info->new_pool = PL_alloc(ga_info->pool_size);
   }
   pool = ga_info->old_pool;

   /*--- Minimize or maximize ---*/
   pool->minimize = ga_info->minimize;

   /*--- Initialize or reset pool ---*/
   PL_generate(ga_info, pool);

   /*--- Don't initialize pool second time ---*/
   if(ga_info->ip_flag != IP_NONE) ga_info->ip_flag = IP_NONE;

   /*--- Make sure best is allocated ---*/
   if(!CH_valid(ga_info->best)) 
      ga_info->best = CH_alloc(ga_info->chrom_len);

   /*--- Save min and max members of initial pool ---*/
   if(ga_info->minimize) 
      CH_copy(pool->chrom[pool->min_index], ga_info->best);
   else
      CH_copy(pool->chrom[pool->max_index], ga_info->best);

   /*--- No mutations yet ---*/
   ga_info->num_mut = 0;
   ga_info->tot_mut = 0;

   /*--- Initial pool report ---*/
   ga_info->iter = -1;
   RP_report(ga_info, pool);

   /*--- Allocate genes for children ---*/
   child1 = CH_alloc(ga_info->chrom_len);
   child2 = CH_alloc(ga_info->chrom_len);
}

/*============================================================================
|                               GA Inner Loop
============================================================================*/
/*----------------------------------------------------------------------------
| A single trial
----------------------------------------------------------------------------*/
void GA_trial(
   GA_Info_Ptr ga_info)
{
   Chrom_Ptr parent1, parent2;


   /*--- Selection ---*/
   parent1 = SE_fun(ga_info, ga_info->old_pool);
   parent2 = SE_fun(ga_info, ga_info->old_pool);

   /*--- Validate parents ---*/
   CH_verify(ga_info, parent1);
   CH_verify(ga_info, parent2);
   
   /*--- Crossover ---*/
   X_fun(ga_info, parent1, parent2, child1, child2);

   /*--- Mutation ---*/
   MU_fun(ga_info, child1);
   MU_fun(ga_info, child2);
   
   /*--- Evaluate children ---*/
   ga_info->EV_fun(child1);
   ga_info->EV_fun(child2);

   /*--- Validate children ---*/
   CH_verify(ga_info, child1);
   CH_verify(ga_info, child2);

   /*--- Replacement ---*/
   RE_fun(ga_info, ga_info->new_pool, parent1, parent2, child1, child2);
   
   /*--- Best So Far? ---*/
   GA_cum(ga_info, child1, child2);
   
   /*--- Update GA system statistics ---*/
   PL_stats(ga_info, ga_info->new_pool);
}

/*============================================================================
|                               Utility
============================================================================*/
/*----------------------------------------------------------------------------
| See if children are best so far
----------------------------------------------------------------------------*/
GA_cum(
   GA_Info_Ptr ga_info,
   Chrom_Ptr   c1,Chrom_Ptr c2)
{
   if(ga_info->minimize) {
      if(c1->fitness < ga_info->best->fitness)
         CH_copy(c1, ga_info->best);
      if(c2->fitness < ga_info->best->fitness)
         CH_copy(c2, ga_info->best);
   } else {
      if(c1->fitness > ga_info->best->fitness)
         CH_copy(c1, ga_info->best);
      if(c2->fitness > ga_info->best->fitness)
         CH_copy(c2, ga_info->best);
   }
}

/*----------------------------------------------------------------------------
| Handle generation gap
----------------------------------------------------------------------------*/
GA_gap(
   GA_Info_Ptr ga_info)
{
   int i, num_clones;
   Chrom_Ptr parent;

   /*--- Disabled ---*/
   if(ga_info->gap <= 0.0) return OK;

   /*--- How many to copy over ---*/
   num_clones = (int)(ga_info->pool_size * ga_info->gap);

   /*--- Append to new pool ---*/
   for(i=0; 
       i < num_clones && ga_info->new_pool->size < ga_info->old_pool->size;
       i++) 
   {
      /*--- Select a chromosome ---*/
      parent = SE_fun(ga_info, ga_info->old_pool);

      /*--- Put it in the new pool ---*/
      PL_append(ga_info->new_pool, parent, TRUE);

   }

   /*--- Make sure stats are updated ---*/
   if(ga_info->new_pool->size == ga_info->old_pool->size)
      PL_stats(ga_info, ga_info->new_pool);

   return OK;
}
/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| Mutation operators 
|
| Bit String Representations
|    MU_simple_invert() - random bit inversion based on mutation rate
|    MU_simple_random() - random bit value based on mutation rate
|
| Any Representation
|    MU_swap()          - random element swap based on mutation rate
|
| Interface
|    MU_table[]   - used in selection of mutation method
|    MU_set_fun() - set and select user defined mutation function
|    MU_select()  - select mutation function by name
|    MU_name()    - get name of current mutation function
|    MU_fun()     - setup and perform current mutation operator
============================================================================*/

/*============================================================================
|                     Mutation interface
============================================================================*/
/*----------------------------------------------------------------------------
| Mutation table
----------------------------------------------------------------------------*/
FN_Table_Type MU_table[] = {
   {  NULL,           NULL             }, /* user defined function */
   { "simple_invert", MU_simple_invert },
   { "simple_random", MU_simple_random },
   { "swap",          MU_swap          },
   { "float_random",  MU_float_random  },
   { "float_rnd_pert",MU_float_rnd_pert},
   { "float_LS",      MU_float_LS      },
   { "float_gauss_pert",MU_float_gauss_pert},
   { NULL,            NULL             }
};

/*----------------------------------------------------------------------------
| Select user mutation method
----------------------------------------------------------------------------*/
MU_set_fun(
   GA_Info_Ptr  ga_info,
   char         *fn_name,
   FN_Ptr       fn_ptr)
{
   return FN_set_fun(ga_info, MU_table, fn_name, fn_ptr, &ga_info->MU_fun);
}

/*----------------------------------------------------------------------------
| Select mutation method
----------------------------------------------------------------------------*/
MU_select(
   GA_Info_Ptr  ga_info,
   char         *fn_name)
{
   return FN_select(ga_info, MU_table, fn_name, &ga_info->MU_fun);
}

/*----------------------------------------------------------------------------
| Mutation name
----------------------------------------------------------------------------*/
char *MU_name(
   GA_Info_Ptr  ga_info)
{
   return FN_name(ga_info, MU_table, ga_info->MU_fun);
}

/*----------------------------------------------------------------------------
| Mutation interface
----------------------------------------------------------------------------*/
MU_fun(
   GA_Info_Ptr ga_info,
   Chrom_Ptr   chrom)
{
   /*--- Random chance to mutate ---*/
   if(RAND_FRAC() <= ga_info->mu_rate && ga_info->MU_fun != NULL) {
      ga_info->MU_fun(ga_info, chrom);
      ga_info->num_mut++;
      ga_info->tot_mut++;
   }
}

/*============================================================================
|                             Mutation operators
============================================================================*/
/*----------------------------------------------------------------------------
| Simple mutation, invert bit
----------------------------------------------------------------------------*/
MU_simple_invert(
   GA_Info_Ptr ga_info,
   Chrom_Ptr chrom)
{
   int idx;

   /*--- Select bit at random ---*/
   idx = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Invert selected bit ---*/
   chrom->gene[idx] = chrom->gene[idx] ? 0 : 1;
}

/*----------------------------------------------------------------------------
| Simple mutation, random bit
----------------------------------------------------------------------------*/
MU_simple_random(GA_Info_Ptr ga_info,
   Chrom_Ptr chrom)
{
   int idx;

   /*--- Select bit at random ---*/
   idx = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Assign random value to bit ---*/
   chrom->gene[idx] = RAND_BIT();
}

/*----------------------------------------------------------------------------
| Swap random elements
----------------------------------------------------------------------------*/
MU_swap(GA_Info_Ptr ga_info,
   Chrom_Ptr chrom)
{
   Gene_Type tmp;
   int       i, j;

   /*--- Select two bits at random (can be same) ---*/
   i = RAND_DOM(chrom->idx_min, chrom->length-1);
   j = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Swap the elements ---*/
   tmp            = chrom->gene[i];
   chrom->gene[i] = chrom->gene[j];
   chrom->gene[j] = tmp;
}


/*----------------------------------------------------------------------------
|   rnd float perturbation in [0..1] -- introduced by claudio 10/02/2004 
----------------------------------------------------------------------------*/
MU_float_rnd_pert(GA_Info_Ptr ga_info,
   Chrom_Ptr chrom)
{
   int       i;

   /*--- Select one element at random ---*/
   i = RAND_DOM(chrom->idx_min, chrom->length);

   //patch
   if(i!=chrom->length)
   {
       /*--- Generate randomly perturbed element ---*/
       chrom->gene[i] += ga_info->pert_range*(1.0 - 2.0*RAND_FRAC()) ;
	   // + 	 ga_info->mut_bias[i]; 
   
   if( chrom->gene[i]>1)
     chrom->gene[i]=1;
   if( chrom->gene[i]<0)
     chrom->gene[i]=0;
   
   }
}


/*----------------------------------------------------------------------------
|   rnd float in [0..1] -- introduced by claudio 10/02/2004 
----------------------------------------------------------------------------*/
MU_float_random(GA_Info_Ptr ga_info,
   Chrom_Ptr chrom)
{
   int       i;

   /*--- Select one element at random ---*/
   i = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Generate random element ---*/
   chrom->gene[i] = RAND_FRAC();
   if( chrom->gene[i]>1)
     chrom->gene[i]=1;
   if( chrom->gene[i]<0)
     chrom->gene[i]=0;

}

/*----------------------------------------------------------------------------
|   Local Search ***experimental*** -- introduced by claudio 10/02/2004 
----------------------------------------------------------------------------*/
MU_float_LS(
     GA_Info_Ptr ga_info,
     Chrom_Ptr chrom)
{

  int   i;
  float old_fit, new_fit, prev_fit, prev_val,tmp;
  

  while(1)
    {
      old_fit=chrom->fitness;
      
      /*--- loop thru genes ---*/
      for(i=chrom->idx_min;i<chrom->length;i++)
	{
	  prev_fit=chrom->fitness;
	  prev_val= chrom->gene[i];
	  
	  chrom->gene[i] += 0.1*(1.0 - 2.0*RAND_FRAC());
	  if( chrom->gene[i]>1)
	    chrom->gene[i]=1;
	  if( chrom->gene[i]<0)
	    chrom->gene[i]=0;
	  
	  //obj_fun(chrom);
	  new_fit=chrom->fitness;
	  
	  if(new_fit>prev_fit)
	    {
	      chrom->gene[i]=prev_val;
	      chrom->fitness=prev_fit;
	    }
	}
      if(chrom->fitness>=old_fit)
	break;
      

    }  

}



/*----------------------------------------------------------------------------
|   gaussian perturbation in [0..1] -- introduced by claudio 28/10/2004 
----------------------------------------------------------------------------*/
MU_float_gauss_pert(GA_Info_Ptr ga_info,
     Chrom_Ptr chrom)
   
{
   int       i;

   //int q=15;
   int c1,c2;
   float c3, pert;

   
   c1 = 32767;  // c1=(1 << q)-1;
   c2 = 10922;  // c2=(c1 / 3);
   c3 = 3.05185e-05;  // c3=1.0/c1;
   

   // NB pert in [-1,1], gaussian: mean=0, var=1
   pert=(2*(RAND_DOM(0,c2)+RAND_DOM(0,c2)+RAND_DOM(0,c2))-3*(c2))*c3;

   //pert=gaussian_random();

   /*--- Select one element at random ---*/
   i = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Generate randomly perturbed element ---*/
   chrom->gene[i] += ga_info->pert_range*pert ;//+ ga_info->mut_bias[i];

   //   printf("gene %d, bias %g\n",i,ga_info->mut_bias[i]);
   
   if( chrom->gene[i]>1)
     chrom->gene[i]=1;
   if( chrom->gene[i]<0)
     chrom->gene[i]=0;
}



double gaussian_random(void)
{
  static int next_gaussian = 0;
  static double saved_gaussian_value;

  double fac, rsq, v1, v2;

  if (next_gaussian == 0) {
    do {
      v1 = 2.0*RAND_FRAC()-1.0;
      v2 = 2.0*RAND_FRAC()-1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    saved_gaussian_value=v1*fac;
    next_gaussian=1;
    return v2*fac;
  } else {
    next_gaussian=0;
    return saved_gaussian_value;
  }
}





/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| Pool management
|
| Functions:
|    PL_alloc()    - allocate a pool
|    PL_resize()   - resize a pool
|    PL_free()     - deallocate a pool
|    PL_valid()    - is a pool valid?
|    PL_reset()    - reset a pool
|    PL_eval()     - evaluate a pool
|    PL_get_num()  - get a number
|    PL_generate() - generate a pool
|    PL_read()     - read a pool from a file
|    PL_rand()     - generate a random pool
|    PL_stats()    - calculate pool statistics
|    PL_index()    - index a pool
|    PL_update_ptf() - update the percent of total fitness in a pool
|    PL_clean()    - empty a pool
|    PL_append()   - append a chrom to a pool
|    PL_insert()   - insert a chrom in a pool at a spec loc
|    PL_remove()   - remove a chrom from a pool
|    PL_move()     - move a chrom in a pool
|    PL_swap()     - swap two chroms in a pool
|    PL_sort()     - sort a pool
============================================================================*/

/*----------------------------------------------------------------------------
| Allocate a pool
----------------------------------------------------------------------------*/
Pool_Ptr PL_alloc(   int max_size)
{
   Pool_Ptr pool;

   /*--- Error check ---*/
   if(max_size <= 0) UT_error("PL_alloc: invalid max_size");

   /*--- Allocate memory for pool ---*/
   pool = (Pool_Ptr)calloc(1, sizeof(Pool_Type));
   if(pool == NULL) UT_error("PL_alloc: pool alloc failed");
   pool->max_size = max_size;

   /*--- Allocate memory for chromosome pointers ---*/
   pool->chrom = (Chrom_Ptr *)calloc(max_size, sizeof(Chrom_Ptr));
   if(pool->chrom == NULL) UT_error("PL_alloc: chrom alloc failed");

   /*--- Put in magic cookie ---*/
   pool->magic_cookie = PL_cookie;

   /*--- Reset pool ---*/
   PL_reset(pool);

   return pool;
}

/*----------------------------------------------------------------------------
| Resize pool
|
| NOTE: chromosomes still in pool after resize are still valid
----------------------------------------------------------------------------*/
PL_resize(
   Pool_Ptr  pool,
   int       new_size)
{
   int old_size;

   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_resize: invalid pool");
   if(new_size <= 0) UT_error("PL_resize: invalid new_size");

   /*--- Free any excess chromosomes ---*/
   if(new_size < pool->max_size) PL_clean(pool, new_size, pool->max_size);

   /*--- Reallocate memory for chromosome pointers ---*/
   pool->chrom = (Chrom_Ptr *)realloc(pool->chrom, new_size*sizeof(Chrom_Ptr));
   if(pool->chrom == NULL) UT_error("PL_resize: chrom realloc failed");

   /*--- Update pool size ---*/
   old_size       = pool->max_size;
   pool->max_size = new_size;

   /*--- Make any new chromosomes NULL ---*/
   if(new_size > old_size) PL_clean(pool, old_size, new_size);
}

/*----------------------------------------------------------------------------
| De-Allocate pool
----------------------------------------------------------------------------*/
PL_free(   Pool_Ptr pool)
{
   /*--- Error check ---*/
  if(!PL_valid(pool)) return GA_ERROR;  // CCC

   /*--- Release memory for chromosome pointers ---*/
   if(pool->chrom != NULL) {

      /*--- Release memory for chromosomes ---*/
      PL_clean(pool, 0, pool->max_size);

      /*--- Release memory for chromosome pointers ---*/
      free(pool->chrom);
      pool->chrom = NULL;
   }

   /*--- Put in a NULL magic cookie ---*/
   pool->magic_cookie = NL_cookie;

   /*--- Release memory for pool ---*/
   free(pool);
}

/*----------------------------------------------------------------------------
| Is pool valid, i.e., has it been allocated by PL_alloc()?
----------------------------------------------------------------------------*/
PL_valid(  Pool_Ptr pool)
{
   /*--- Check for NULL pointers ---*/
   if(pool == NULL) return FALSE;
   if(pool->chrom == NULL) return FALSE;

   /*--- Check for magic cookie ---*/
   if(pool->magic_cookie != PL_cookie) return FALSE;

   /*--- Otherwise valid ---*/
   return TRUE;
}

/*----------------------------------------------------------------------------
| Reset pool
----------------------------------------------------------------------------*/
PL_reset(  Pool_Ptr pool)
{
   int i;

   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_reset: invalid pool");

   /*--- Reset chromosomes ---*/
   for(i=0; i<pool->max_size; i++) {
      if(CH_valid(pool->chrom[i]))
         CH_reset(pool->chrom[i]);
      else
         pool->chrom[i] = NULL;
   }

   /*--- Reset pool ---*/
   pool->size = 0;
   pool->total_fitness = 0.0;
   pool->min = 0.0;
   pool->max = 0.0;
   pool->ave = 0.0;
   pool->var = 0.0;
   pool->dev = 0.0;
   pool->min_index = -1;
   pool->max_index = -1;
   pool->minimize = TRUE;
   pool->sorted   = FALSE;
}

/*----------------------------------------------------------------------------
| Evaluate pool
----------------------------------------------------------------------------*/
PL_eval(
   GA_Info_Ptr ga_info,
   Pool_Ptr pool)
{
   int i;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("PL_eval: invalid ga_info");
   if(!PL_valid(pool)) UT_error("PL_eval: invalid pool");

   /*--- Evaluate each chromosome ---*/
   for(i = 0; i < pool->size; i++) {
      ga_info->EV_fun(pool->chrom[i]);
	 // printf("%d -> %g\n",i,pool->chrom[i]->fitness);	
   }
}

/*============================================================================
|                            Generate pool
============================================================================*/
/*----------------------------------------------------------------------------
| Read a number for PL_generate()
----------------------------------------------------------------------------*/
char *PL_get_num(   FILE   *fid)
{
   static char str[80]; /* This MUST be static */
   int  ch, len = 0;

   /*--- Search for a digit ---*/
   while(TRUE) {

      /*--- Get a character ---*/
      ch = fgetc(fid);

      /*--- Error or Quit ---*/
      if(ch == EOF || ch == 'q' || ch == 'Q') return NULL;

      /*--- Digit ---*/
      if(isdigit(ch)) break;

      /*--- Comment: skip to end of line ---*/
      if(ch == '#') 
         while((ch = fgetc(fid)) != '\n') 
            ;
   }

   /*--- Make sure digit found above ---*/
   if(!isdigit(ch)) UT_error("PL_get_num: bad digit");

   /*--- Digit found, now put into str ---*/
   while(len < 80 && ch != EOF && !isspace(ch) && ch != '#') {
      str[len++] = ch;
      ch = fgetc(fid);
   }
   str[len] = 0;

   /*--- Return number in string ---*/
   return str;
}

/*----------------------------------------------------------------------------
| Generate pool according to ga_info->ip_flag:
|
|   IP_INTERACTIVE: interactively read pool 
|   IP_FROM_FILE:   read pool from file given in file_name argument  
|   IP_RANDOM:      random pool limited by pool_size & chrom_len
|   IP_RANDOM01:    random pool in [0,1] limited by pool_size & chrom_len
|   IP_NONE:        do nothing
----------------------------------------------------------------------------*/
PL_generate(
   GA_Info_Ptr ga_info,
   Pool_Ptr pool)
{
   FILE *fid;
   long chrom_len;
   char *sptr;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("PL_generate: invalid ga_info");
   if(!PL_valid(pool)) UT_error("PL_generate: invalid pool");

  
   /*--- Get initial pool data ---*/
   switch(ga_info->ip_flag) {

      case IP_INTERACTIVE:

         /*--- Get chrom_len ---*/
         puts("\nEnter chromosome length:");
         if((sptr = PL_get_num(stdin)) == NULL) {
            UT_error("PL_generate: No chrom_len was read");
         }
         if(sscanf(sptr,"%ld",&chrom_len) != 1) 
            UT_error("PL_generate: error reading chrom_len");
         if(chrom_len <= 0) 
            UT_error("PL_generate: invalid chrom_len");
         ga_info->chrom_len = chrom_len;

         /*--- Get the pool ---*/
         puts("\nEnter initial pool (`q' to quit):");
         PL_read(pool, chrom_len, stdin);

         break;

      case IP_FROM_FILE: 

         /*--- Open the file ---*/
         if((fid = fopen(ga_info->ip_data, "r")) == NULL)
            UT_error("PL_generate: Invalid data file");

         /*--- Get chrom_len ---*/
         if((sptr = PL_get_num(fid)) == NULL) {
            UT_error("PL_generate: No chrom_len was read");
         }
         if(sscanf(sptr,"%ld",&chrom_len) != 1) 
            UT_error("PL_generate: error reading chrom_len");
         if(chrom_len <= 0) 
            UT_error("PL_generate: invalid chrom_len");
         ga_info->chrom_len = chrom_len;

         /*--- Get the pool ---*/
         PL_read(pool, chrom_len, fid);

         /*--- Close the file ---*/
         fclose(fid);
         break;

      case IP_RANDOM:
         PL_rand(pool, ga_info->pool_size, ga_info->chrom_len, 
                    ga_info->datatype);
         break;
      case IP_RANDOM01:
         PL_rand01(pool, ga_info->pool_size, ga_info->chrom_len, 
                    ga_info->datatype);
         break;

      case IP_NONE:
         break;

      default:
         UT_error("PL_generate: invalid ip_flag");
   }


   /*--- Evaluate the pool ---*/
   PL_eval(ga_info, pool);


   /*--- Compute pool statistics ---*/
   PL_stats(ga_info, pool);


}

/*----------------------------------------------------------------------------
| Initialize pool from file.
----------------------------------------------------------------------------*/
PL_read(
   Pool_Ptr pool,
   long     chrom_len,
   FILE     *fid)
{
   Chrom_Ptr  chrom;
   double     gene;
   long       i;
   char       *sptr;

   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_read: invalid pool");
   if(chrom_len <= 0) UT_error("PL_read: invalid chrom_len");
   if(fid == NULL) UT_error("PL_read: invalid fid");

   /*--- While there is a chromosome to read ---*/
   while(TRUE) {

      /*--- Allocate or reuse a chromosome ---*/
      if(CH_valid(pool->chrom[pool->size])) {
         chrom = pool->chrom[pool->size];
         CH_resize(chrom, chrom_len);
         pool->chrom[pool->size] = NULL;
      } else {
         chrom = CH_alloc(chrom_len);
      }

      /*--- Read genes ---*/
      for(i=0; i < chrom_len; i++) {

         /*--- Read a gene ---*/
         sptr = PL_get_num(fid);

         /*--- Get number from sptr if valid ---*/
         if(sptr == NULL || sscanf(sptr,"%lf",&gene) != 1) {

            /*--- Free chromosome ---*/
            CH_free(chrom);

            /*--- EOF in middle of chromosome ---*/
            if(i != 0 && (sptr == NULL || *sptr != 'q')) {
               UT_warn("PL_read: premature eof reading chromosome");
            }

            /*--- EOF: no more to read ---*/
            return OK;  // CCC
         }

         /*--- Valid gene was read ---*/
         chrom->gene[i] = (Gene_Type)gene;
      }

      /*--- Put the chromosome into the pool ---*/
      PL_append(pool, chrom, FALSE);
   }
}

/*----------------------------------------------------------------------------
| Initialize random pool.
|
| NOTE: Alleles for datatypes DT_INT and DT_REAL are generated from an
|       arbitrarily chosen domain.  There really needs to be a way for
|       the user to indicate a domain for EACH gene.
----------------------------------------------------------------------------*/
PL_rand(
   Pool_Ptr pool,
   int      pool_size,int   chrom_len, int  datatype)
{
   int       i, j, idx;
   Chrom_Ptr chrom;

   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_rand: invalid pool");
   if(pool_size < 0) UT_error("PL_rand: negative pool_size");
   if(chrom_len < 0) UT_error("PL_rand: negative chrom_len");
   switch(datatype) {
      case DT_BIT:
      case DT_INT:
      case DT_INT_PERM:
      case DT_REAL:
         break;
      default: UT_error("PL_rand: Invalid datatype");
   }

 
   /*--- Make sure there is enough space for pool ---*/
   if(pool_size > pool->max_size) PL_resize(pool, pool_size);

   
   /*--- For each chromosome to be generated ---*/
   for(i = 0; i < pool_size; i++) {

      /*--- Allocate or reuse a chromosome ---*/
      if(CH_valid(pool->chrom[pool->size])) {
         chrom = pool->chrom[pool->size];
         CH_resize(chrom, chrom_len);
         pool->chrom[pool->size] = NULL;
      } else {
         chrom = CH_alloc(chrom_len);
      }

      /*--- Generate random genes ---*/
      switch(datatype) {

         case DT_BIT:
            /*--- Random bit ---*/
            for(j = 0; j < chrom_len; j++)
               chrom->gene[j] = (Gene_Type)RAND_BIT();
            chrom->length = chrom_len;
            break;

         case DT_INT:
            /*--- Random integers from an arbitrary domain ---*/
            for(j = 0; j < chrom_len; j++)
               chrom->gene[j] = (Gene_Type)RAND_DOM(0,chrom_len);
            chrom->length = chrom_len;
            break;

         case DT_INT_PERM:
	   
   
            /*--- Random permutations of integers ---*/
            for(j = 0; j < chrom_len; j++) 
               chrom->gene[j] = (Gene_Type)-1;
	    
	    

            for(j = 0; j < chrom_len; j++) {
               idx = RAND_DOM(0,chrom_len-1);
	 
	       while(chrom->gene[idx] != -1) 
                  idx = RAND_DOM(0,chrom_len-1);

               chrom->gene[idx] = (Gene_Type)(j + 1);

	 
            }
            chrom->length = chrom_len;
	
            break;

         case DT_REAL:
            /*--- Random reals from an arbitrary domain ---*/
            for(j = 0; j < chrom_len; j++)
               chrom->gene[j] = 
                  (double)RAND_DOM(0,chrom_len-1) + (double)RAND_FRAC();
            chrom->length = chrom_len;
            break;

         /*--- Should never reach here ---*/
         default: UT_error("PL_rand: Invalid datatype");
      }

      /*--- Put the chromosome into the pool ---*/
      PL_append(pool, chrom, FALSE);
   }
}
/*----------------------------------------------------------------------------
| Initialize random pool in [0,1].
----------------------------------------------------------------------------*/
PL_rand01(
   Pool_Ptr pool,
   int      pool_size, int  chrom_len, int  datatype)
{
   int       i, j, idx;
   Chrom_Ptr chrom;

   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_rand: invalid pool");
   if(pool_size < 0) UT_error("PL_rand: negative pool_size");
   if(chrom_len < 0) UT_error("PL_rand: negative chrom_len");
   switch(datatype) {
      case DT_BIT:
      case DT_INT:
      case DT_INT_PERM:
      case DT_REAL:
         break;
      default: UT_error("PL_rand: Invalid datatype");
   }

   /*--- Make sure there is enough space for pool ---*/
   if(pool_size > pool->max_size) PL_resize(pool, pool_size);

   /*--- For each chromosome to be generated ---*/
   for(i = 0; i < pool_size; i++) {

      /*--- Allocate or reuse a chromosome ---*/
      if(CH_valid(pool->chrom[pool->size])) {
         chrom = pool->chrom[pool->size];
         CH_resize(chrom, chrom_len);
         pool->chrom[pool->size] = NULL;
      } else {
         chrom = CH_alloc(chrom_len);
      }

      /*--- Generate random genes ---*/
      switch(datatype) {

         case DT_BIT:
            /*--- Random bit ---*/
            for(j = 0; j < chrom_len; j++)
               chrom->gene[j] = (Gene_Type)RAND_BIT();
            chrom->length = chrom_len;
            break;

         case DT_INT:
            /*--- Random integers from an arbitrary domain ---*/
            for(j = 0; j < chrom_len; j++)
               chrom->gene[j] = (Gene_Type)RAND_DOM(0,chrom_len);
            chrom->length = chrom_len;
            break;

         case DT_INT_PERM:
            /*--- Random permutations of integers ---*/
            for(j = 0; j < chrom_len; j++) 
               chrom->gene[j] = (Gene_Type)-1;
            for(j = 0; j < chrom_len; j++) {
               idx = RAND_DOM(0,chrom_len-1);
               while(chrom->gene[idx] != -1) 
                  idx = RAND_DOM(0,chrom_len-1);
               chrom->gene[idx] = (Gene_Type)(j + 1);
            }
            chrom->length = chrom_len;
            break;

         case DT_REAL:
            /*--- Random reals from an arbitrary domain ---*/
            for(j = 0; j < chrom_len; j++)
               chrom->gene[j] =  (double)RAND_FRAC();
            chrom->length = chrom_len;
            break;

         /*--- Should never reach here ---*/
         default: UT_error("PL_rand: Invalid datatype");
      }

      /*--- Put the chromosome into the pool ---*/
      PL_append(pool, chrom, FALSE);
   }
}

/*============================================================================
|                            Pool statistics
============================================================================*/
/*----------------------------------------------------------------------------
| Compute statistics for a pool
----------------------------------------------------------------------------*/
PL_stats(
   GA_Info_Ptr ga_info,
   Pool_Ptr    pool)
{
   unsigned i, min_index, max_index;
   double   min, max, total, var;
   int      no_variance, sorted;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("PL_stats: invalid ga_info");
   if(!PL_valid(pool)) UT_error("PL_stats: invalid pool");


   /*--- Trivial cases ---*/
   if(pool->size <= 0) {
      pool->min           = 0.0;
      pool->max           = 0.0;
      pool->ave           = 0.0;
      pool->var           = 0.0;
      pool->dev           = 0.0;
      pool->total_fitness = 0.0;
      pool->min_index     = -1;
      pool->max_index     = -1;
      pool->best_index    = -1;
      return OK;
   } else if(pool->size == 1) {
      pool->min           = pool->chrom[0]->fitness;
      pool->max           = pool->chrom[0]->fitness;
      pool->ave           = pool->chrom[0]->fitness;
      pool->var           = 0.0;
      pool->dev           = 0.0;
      pool->total_fitness = pool->chrom[0]->fitness;
      pool->min_index     = 0;
      pool->max_index     = 0;
      pool->best_index    = 0;
      return OK;
   }

   /*--- Initialize stats ---*/
   min         = pool->chrom[0]->fitness;
   max         = pool->chrom[0]->fitness;
   min_index   = 0;
   max_index   = 0;
   var         = 0.0;
   total       = 0.0;
   no_variance = TRUE;
   sorted      = TRUE;

   /*--- Collect statistics and find total fitness ---*/
   for(i = 0; i < pool->size; i++) {

      /*--- Error check ---*/
      if(!CH_valid(pool->chrom[i])) UT_error("PL_stats: invalid chrom");

      /*--- Comparative statistics ---*/
      if(i != 0) {

         /*--- Check for variance in fitness ---*/
         if(no_variance && 
            pool->chrom[i]->fitness != pool->chrom[i-1]->fitness)
         {
            no_variance = FALSE;
         }

         /*--- Check for not sorted ---*/
         if(sorted && pool->minimize) {
            if(pool->chrom[i]->fitness < pool->chrom[i-1]->fitness)
               sorted = FALSE;
         } else if(sorted) {
            if(pool->chrom[i]->fitness > pool->chrom[i-1]->fitness)
               sorted = FALSE;
         }
      }

      /*--- Less than min? ---*/
      if(pool->chrom[i]->fitness < min) {
         min = pool->chrom[i]->fitness;
         min_index = i;
      }

      /*--- Greater than max? ---*/
      if(pool->chrom[i]->fitness > max) {
         max = pool->chrom[i]->fitness;
         max_index = i;
      }

      /*--- Update cum. values ---*/
      total += pool->chrom[i]->fitness;
      var   += (float)pool->chrom[i]->fitness * 
               (float)pool->chrom[i]->fitness;

      /*--- Make sure index is set ---*/
      pool->chrom[i]->index = i;
   }

   /*--- Update pool statistics ---*/
   pool->min = min;
   pool->max = max;
   pool->ave = total / pool->size;
   pool->min_index = min_index;
   pool->max_index = max_index;
   if(pool->minimize) pool->best_index = min_index;
   else               pool->best_index = max_index;
   pool->total_fitness = total;

   /*--- Variance and standard deviation ---*/
   var = (var - (pool->ave * total)) / (pool->size - 1);
   if(no_variance || var <= 0.0) {
      pool->var = 0.0;
      pool->dev = 0.0;
      ga_info->converged = TRUE;
   } else {
      pool->var = var;
      pool->dev = sqrt(var);
      ga_info->converged = FALSE;
   }
}

/*----------------------------------------------------------------------------
| Update indices in the pool
----------------------------------------------------------------------------*/
PL_index(   Pool_Ptr pool)
{
   int i;

   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_index: invalid pool");

   /*--- Set index for each chromosome ---*/
   for(i=0; i<pool->size; i++)
      pool->chrom[i]->index = i;
}

/*----------------------------------------------------------------------------
| Update ptf (percent of total fitness) for each chromosome
----------------------------------------------------------------------------*/
PL_update_ptf(
   GA_Info_Ptr ga_info,
   Pool_Ptr    pool)
{
   int i, sf_changed, all_positive;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("PL_update_ptf: invalid ga_info");
   if(!PL_valid(pool)) UT_error("PL_update_ptf: invalid pool");

   /*--- Compute scale factor (to ensure positive fitness) ---*/
   sf_changed = FALSE;
   all_positive = TRUE;
   for(i = 0; i < pool->size; i++) {

      /*--- If scaled fitness not positive ---*/
      if((pool->chrom[i]->fitness + ga_info->scale_factor) <= 0) {

         /*--- Adjust scale_factor so scaled fitness is 1.0 ---*/
         ga_info->scale_factor += 
            1.0 - (pool->chrom[i]->fitness + ga_info->scale_factor);

         /*--- The scale_factor has been changed ---*/
         sf_changed = TRUE;
      } 

      /*--- Make sure scale factor is still needed ---*/
      if(pool->chrom[i]->fitness <= 0) all_positive = FALSE;
   }

   /*--- Scale factor no longer needed ---*/
   if(all_positive) {
      if(ga_info->scale_factor > 0.0) sf_changed = TRUE;
      ga_info->scale_factor = 0.0;
   }

   /*--- Report change in scale factor ---*/
   if(sf_changed && 
      ga_info->rp_type != RP_NONE && 
      ga_info->rp_type != RP_MINIMAL) 
   {
      fprintf(ga_info->rp_fid,
              "New scale factor = %G\n", 
               ga_info->scale_factor);
   }

   /*--- Find total fitness ---*/
   pool->total_fitness = 0;
   for(i = 0; i < pool->size; i++) {
      pool->total_fitness += pool->chrom[i]->fitness + ga_info->scale_factor;
   }

   /*--- Update ptf for each chromosome (minimize) ---*/
   if(ga_info->minimize) {
      double new_total_fitness;

      /*--- Compute new fitness and total_fitness ---*/
      new_total_fitness = 0.0;
      for(i = 0; i < pool->size; i++) {

         /*--- Failed scaling leads to divide by zero ---*/
         if((pool->chrom[i]->fitness + ga_info->scale_factor) <= 0.0)
            UT_error("PL_update_ptf: fitness + scale <= 0.0");

         /*--- Save new fitness in ptf ---*/
         pool->chrom[i]->ptf = 
            pool->total_fitness / 
            (pool->chrom[i]->fitness + ga_info->scale_factor);

         /*--- New total fitness based on new ptf ---*/
         new_total_fitness += pool->chrom[i]->ptf;
      }

      /*--- Compute new ptf ---*/
      for(i = 0; i < pool->size; i++) {

         /*--- Failed scaling leads to divide by zero ---*/
         if(new_total_fitness <= 0.0) 
            UT_error("PL_update_ptf: new_total_fitness <= 0.0");

         /*--- New ptf ---*/
         pool->chrom[i]->ptf *= 100.0 / new_total_fitness; 
      }

   /*--- Update ptf for each chromosome (maximize) ---*/
   } else {

      /*--- Compute new ptf ---*/
      for(i = 0; i < pool->size; i++) {

         /*--- Failed scaling leads to divide by zero ---*/
         if(pool->total_fitness <= 0.0)
            UT_error("PL_update_ptf: pool->total_fitness <= 0.0");

         /*--- New ptf ---*/
         pool->chrom[i]->ptf = 
            ((pool->chrom[i]->fitness + ga_info->scale_factor) / 
            pool->total_fitness) * 100.0;
      }
   }

   return OK;
}

/*============================================================================
|                            Pool manipulation
============================================================================*/
/*----------------------------------------------------------------------------
| Clean out the pool
----------------------------------------------------------------------------*/
PL_clean(
   Pool_Ptr pool,
   int      idx_min, int idx_max)
{
   int i;

   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_clean: invalid pool");
   if(idx_min < 0 || idx_min > pool->max_size) 
      UT_error("PL_clean: invalid idx_min");
   if(idx_max < 0 || idx_max > pool->max_size) 
      UT_error("PL_clean: invalid idx_max");
   if(idx_min > idx_max) UT_error("PL_clean: idx_min > idx_max");

   for(i=idx_min; i<idx_max; i++)
      PL_remove(pool, i);
}
 
/*----------------------------------------------------------------------------
| Append a chromosome to the pool
----------------------------------------------------------------------------*/
PL_append(
   Pool_Ptr  pool,
   Chrom_Ptr chrom,
   int       make_copy)
{
   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_append: invalid pool");
   if(!CH_valid(chrom)) UT_error("PL_append: invalid chrom");

   /*--- Append = insert at end of pool ---*/
   PL_insert(pool, (int)pool->size, chrom, make_copy);
   ++(pool->size); 
}

/*----------------------------------------------------------------------------
| Replace a chromosome in the pool
----------------------------------------------------------------------------*/
PL_insert(
   Pool_Ptr  pool,
   int       index,
   Chrom_Ptr chrom,
   int       make_copy)
{
   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_insert: invalid pool");
   if(!CH_valid(chrom)) UT_error("PL_insert: invalid chrom");
   if(index < 0 || index > pool->max_size) 
      UT_error("PL_insert: invalid index");

   /*--- Realloc for more space ---*/
   if(index == pool->max_size) 
      PL_resize(pool, pool->max_size + PL_ALLOC_SIZE);
 
   /*--- Insert the chromosome ---*/
   if(make_copy) {
      if(!CH_valid(pool->chrom[index])) 
         pool->chrom[index] = CH_alloc(chrom->length);
      CH_copy(chrom, pool->chrom[index]);
   } else {
      if(CH_valid(pool->chrom[index])) 
         PL_remove(pool, index);
      pool->chrom[index] = chrom;
   }
}
 
/*----------------------------------------------------------------------------
| Remove a chromosome from the pool
----------------------------------------------------------------------------*/
PL_remove(
   Pool_Ptr pool,
   int      index)
{
   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_remove: invalid pool");
   if(index < 0 || index >= pool->max_size) 
      UT_error("PL_remove: invalid index");

   if(CH_valid(pool->chrom[index])) CH_free(pool->chrom[index]);
   pool->chrom[index] = NULL;
}

/*----------------------------------------------------------------------------
| Move a chromosome in the pool
----------------------------------------------------------------------------*/
PL_move(
   Pool_Ptr pool,
   int      idx_src, int idx_dst)
{
   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_move: invalid pool");
   if(idx_src < 0 || idx_src >= pool->max_size) 
      UT_error("PL_move: invalid idx_src");
   if(idx_dst < 0 || idx_dst >= pool->max_size) 
      UT_error("PL_move: invalid idx_dst");
 
   if(CH_valid(pool->chrom[idx_dst])) PL_remove(pool, idx_dst);
   pool->chrom[idx_dst] = pool->chrom[idx_src];
   pool->chrom[idx_src] = NULL;
}

/*----------------------------------------------------------------------------
| Swap chromosomes in the pool
----------------------------------------------------------------------------*/
PL_swap(
   Pool_Ptr pool,
   int      idx1, int idx2)
{
   Chrom_Ptr tmp;

   /*--- Error check ---*/
   if(!PL_valid(pool)) UT_error("PL_swap: invalid pool");
   if(idx1 < 0 || idx1 >= pool->max_size) 
      UT_error("PL_swap: invalid idx1");
   if(idx2 < 0 || idx2 >= pool->max_size) 
      UT_error("PL_swap: invalid idx2");
 
   tmp               = pool->chrom[idx1];
   pool->chrom[idx1] = pool->chrom[idx2];
   pool->chrom[idx2] = tmp;
}

/*----------------------------------------------------------------------------
| Sort comparison function for minimizing GA (ascending fitness)
----------------------------------------------------------------------------*/
static PL_cmp_min(   Chrom_Ptr *a,Chrom_Ptr *b)
{
   if((*a)->fitness < (*b)->fitness)
      return -1;
   else if((*a)->fitness > (*b)->fitness)
      return 1;
   else
      return 0;
}

/*----------------------------------------------------------------------------
| Sort comparison function for maximizing GA (descending fitness)
----------------------------------------------------------------------------*/
static PL_cmp_max(Chrom_Ptr *a,Chrom_Ptr *b)
{
   if((*a)->fitness < (*b)->fitness)
      return 1;
   else if((*a)->fitness > (*b)->fitness)
      return -1;
   else
      return 0;
}

/*----------------------------------------------------------------------------
| Sort the pool
----------------------------------------------------------------------------*/
PL_sort(
   GA_Info_Ptr ga_info,
   Pool_Ptr    pool)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("PL_sort: invalid ga_info");
   if(!PL_valid(pool)) UT_error("PL_sort: invalid pool");

   /*--- Sort based on objective ---*/
   if(ga_info->minimize)
      qsort(pool->chrom, pool->size, sizeof(Chrom_Ptr), PL_cmp_min);
   else
      qsort(pool->chrom, pool->size, sizeof(Chrom_Ptr), PL_cmp_max);

   /*--- Reindex ---*/
   PL_index(pool);

   /*--- Pool is now sorted ---*/
   pool->sorted = TRUE;
}
/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| Replacement operators 
|
| Operators
|    RE_append()        - simply append to pool
|    RE_by_rank()       - insert into pool by rank
|       RE_do_by_rank() - helper for RE_by_rank()
|    RE_first_weaker()  - replace first weaker member of pool
|    RE_weakest()       - replace weakest member of pool
|
| Interface
|    RE_table[]   - used in selection of replacement method
|    RE_set_fun() - set and select user defined replacement function
|    RE_select()  - select replacement method by name
|    RE_name()    - get name of current crossover function
|    RE_fun()     - setup and perform current crossover operator
|    
| Utility
|    RE_pick_best() - pick the best two out of four chromosomes
============================================================================*/

/*============================================================================
|                     Replacement interface
============================================================================*/
/*----------------------------------------------------------------------------
| Replacement table
----------------------------------------------------------------------------*/
FN_Table_Type RE_table[] = {
   { NULL,           NULL            },  /* user defined function */
   { "append",       RE_append       },
   { "by_rank",      RE_by_rank      },
   { "first_weaker", RE_first_weaker },
   { "weakest",      RE_weakest      },
   { NULL,           NULL            }
};

/*----------------------------------------------------------------------------
| Select user replacement method
----------------------------------------------------------------------------*/
RE_set_fun(
   GA_Info_Ptr  ga_info,
   char         *fn_name,
   FN_Ptr       fn_ptr)
{
   return FN_set_fun(ga_info, RE_table, fn_name, fn_ptr, &ga_info->RE_fun);
}

/*----------------------------------------------------------------------------
| Select replacement method
----------------------------------------------------------------------------*/
RE_select(
   GA_Info_Ptr  ga_info,
   char         *fn_name)
{
   return FN_select(ga_info, RE_table, fn_name, &ga_info->RE_fun);
}

/*----------------------------------------------------------------------------
| Replacement name
----------------------------------------------------------------------------*/
char *RE_name(   GA_Info_Ptr  ga_info)
{
   return FN_name(ga_info, RE_table, ga_info->RE_fun);
}

/*----------------------------------------------------------------------------
| Replacement interface
----------------------------------------------------------------------------*/
RE_fun(
   GA_Info_Ptr    ga_info,
   Pool_Ptr       pool,
   Chrom_Ptr      p1, Chrom_Ptr p2,Chrom_Ptr c1,Chrom_Ptr c2)
{
   /*--- Error checking ---*/
   if(pool == NULL || pool->size < 0) UT_error("RE_fun: invalid pool");
   if(ga_info == NULL) UT_error("RE_fun: invalid ga_info");
   if(p1 == NULL) UT_error("RE_fun: invalid p1");
   if(p2 == NULL) UT_error("RE_fun: invalid p2");
   if(c1 == NULL) UT_error("RE_fun: invalid c1");
   if(c2 == NULL) UT_error("RE_fun: invalid c2");
   if(ga_info->RE_fun == NULL) 
      UT_error("RE_fun: null replacement function");

   if(ga_info->elitist)
      RE_pick_best(ga_info, p1, p2, c1, c2);

   ga_info->RE_fun(ga_info, pool, p1, p2, c1, c2);
}

/*============================================================================
|                             Replacement Methods
============================================================================*/
/*----------------------------------------------------------------------------
| Append children to pool
----------------------------------------------------------------------------*/
RE_append(
   GA_Info_Ptr    ga_info,
   Pool_Ptr       pool,
   Chrom_Ptr      p1,Chrom_Ptr p2,Chrom_Ptr c1,Chrom_Ptr c2)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RE_append: invalid ga_info");
   if(!PL_valid(pool)) UT_error("RE_append: invalid pool");
   if(!CH_valid(p1)) UT_error("RE_append: invalid p1");
   if(!CH_valid(p2)) UT_error("RE_append: invalid p2");
   if(!CH_valid(c1)) UT_error("RE_append: invalid c1");
   if(!CH_valid(c2)) UT_error("RE_append: invalid c2");

   /*--- Error conditions ---*/
   if(pool == NULL) UT_error("RE_append: null pool");
   if(pool->size < 0) UT_error("RE_append: invalid pool");

   PL_append(pool, c1, TRUE);
   PL_append(pool, c2, TRUE);
}

/*----------------------------------------------------------------------------
| Insert in ranked order
----------------------------------------------------------------------------*/
RE_by_rank(GA_Info_Ptr    ga_info,
   Pool_Ptr       pool,
   Chrom_Ptr      p1, Chrom_Ptr p2,Chrom_Ptr c1,Chrom_Ptr c2)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RE_by_rank: invalid ga_info");
   if(!PL_valid(pool)) UT_error("RE_by_rank: invalid pool");
   if(!CH_valid(p1)) UT_error("RE_by_rank: invalid p1");
   if(!CH_valid(p2)) UT_error("RE_by_rank: invalid p2");
   if(!CH_valid(c1)) UT_error("RE_by_rank: invalid c1");
   if(!CH_valid(c2)) UT_error("RE_by_rank: invalid c2");
   /*--- PATCH 1 BEGIN ---*/
   /* Many thanks to Paul-Erik Raue (peraue@cs.vu.nl) 
    * for finding this bug. 
    * 
    * This operation is invalid under the generational model 
    */
   if(!strcmp(GA_name(ga_info),"generational"))
      UT_error("RE_by_rank: invalid under generational model");
   /*--- PATCH 1 END ---*/

   RE_do_by_rank(ga_info, pool, c1);
   RE_do_by_rank(ga_info, pool, c2);
}

/*----------------------------------------------------------------------------
| Replace first weaker
----------------------------------------------------------------------------*/
RE_first_weaker(GA_Info_Ptr    ga_info,
   Pool_Ptr       pool,
   Chrom_Ptr      p1, Chrom_Ptr p2,Chrom_Ptr c1,Chrom_Ptr c2)
{
   int       i;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RE_first_weaker: invalid ga_info");
   if(!PL_valid(pool)) UT_error("RE_first_weaker: invalid pool");
   if(!CH_valid(p1)) UT_error("RE_first_weaker: invalid p1");
   if(!CH_valid(p2)) UT_error("RE_first_weaker: invalid p2");
   if(!CH_valid(c1)) UT_error("RE_first_weaker: invalid c1");
   if(!CH_valid(c2)) UT_error("RE_first_weaker: invalid c2");
   /*--- PATCH 1 BEGIN ---*/
   /* Many thanks to Paul-Erik Raue (peraue@cs.vu.nl) 
    * for finding this bug. 
    * 
    * This operation is invalid under the generational model 
    */
   if(!strcmp(GA_name(ga_info),"generational"))
      UT_error("RE_first_weaker: invalid under generational model");
   /*--- PATCH 1 END ---*/

   /*--- Insert c1 ---*/
   for(i=0; i<pool->size; i++) 
      if(CH_cmp(ga_info, pool->chrom[i], c1) > 0) {
         PL_insert(pool, i, c1, TRUE);
         break;
      }

   /*--- Insert c2 ---*/
   for(i=0; i<pool->size; i++) 
      if(CH_cmp(ga_info, pool->chrom[i], c2) > 0) {
         PL_insert(pool, i, c2, TRUE);
         break;
      }
}

/*----------------------------------------------------------------------------
| Replace weakest
----------------------------------------------------------------------------*/
RE_weakest(GA_Info_Ptr    ga_info,
   Pool_Ptr       pool,
   Chrom_Ptr      p1, Chrom_Ptr p2,Chrom_Ptr c1,Chrom_Ptr c2)
{
   Chrom_Ptr weakest;
   int       i, index;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RE_weakest: invalid ga_info");
   if(!PL_valid(pool)) UT_error("RE_weakest: invalid pool");
   if(!CH_valid(p1)) UT_error("RE_weakest: invalid p1");
   if(!CH_valid(p2)) UT_error("RE_weakest: invalid p2");
   if(!CH_valid(c1)) UT_error("RE_weakest: invalid c1");
   if(!CH_valid(c2)) UT_error("RE_weakest: invalid c2");
   /*--- PATCH 1 BEGIN ---*/
   /* Many thanks to Paul-Erik Raue (peraue@cs.vu.nl) 
    * for finding this bug. 
    * 
    * This operation is invalid under the generational model 
    */
   if(!strcmp(GA_name(ga_info),"generational"))
      UT_error("RE_weakest: invalid under generational model");
   /*--- PATCH 1 END ---*/

   /*--- Insert c1 ---*/
   weakest = pool->chrom[0];
   index = 0;
   for(i=0; i<pool->size; i++) {
      if(CH_cmp(ga_info, pool->chrom[i], weakest) >= 0) {
         weakest = pool->chrom[i];
         index = i;
      }
   }
   if(CH_cmp(ga_info, weakest, c1) >= 0)
      PL_insert(pool, index, c1, TRUE);

   /*--- Insert c2 ---*/
   weakest = pool->chrom[0];
   index = 0;
   for(i=0; i<pool->size; i++) {
      if(CH_cmp(ga_info, pool->chrom[i], weakest) >= 0) {
         weakest = pool->chrom[i];
         index = i;
      }
   }
   if(CH_cmp(ga_info, weakest, c2) >= 0)
      PL_insert(pool, index, c2, TRUE);
}

/*============================================================================
|                             Utility functions
============================================================================*/
/*----------------------------------------------------------------------------
| Insert a chromosome by rank
----------------------------------------------------------------------------*/
RE_do_by_rank(
   GA_Info_Ptr    ga_info,
   Pool_Ptr       pool,
   Chrom_Ptr      chrom)
{
   int i, cmp;
   int min, max, med;

   /*--- Failure ---*/
   if(CH_cmp(ga_info, pool->chrom[pool->size-1], chrom) <= 0) return OK;

   /*--- Take place of last chrom ---*/
   PL_insert(pool, pool->size-1, chrom, TRUE);

   /*--- Bubble new chrom into position ---*/
   for(i=pool->size-1; i > 0; i--) {
      cmp = CH_cmp(ga_info, pool->chrom[i-1], pool->chrom[i]);
      if(cmp <= 0) break;
      PL_swap(pool, i-1, i);
   }

   return OK;
}

/*----------------------------------------------------------------------------
| Pick best 2 of 4 chromosomes
----------------------------------------------------------------------------*/
RE_pick_best(GA_Info_Ptr    ga_info,
   Chrom_Ptr      p1, Chrom_Ptr p2,Chrom_Ptr c1,Chrom_Ptr c2)
{
   int xp1, xp2;

   /*--- Replace worst child with p1 if p1 better ---*/
   if(CH_cmp(ga_info, c1, c2) > 0) {
      if(CH_cmp(ga_info, c1, p1) > 0) {
         xp1 = c1->xp1; xp2 = c1->xp2;
         CH_copy(p1, c1);
         c1->xp1 = xp1; c1->xp2 = xp2;
      }
   } else {
      if(CH_cmp(ga_info, c2, p1) > 0) {
         xp1 = c2->xp1; xp2 = c2->xp2;
         CH_copy(p1, c2);
         c2->xp1 = xp1; c2->xp2 = xp2;
      }
   }

   /*--- Replace worst child with p2 if p2 better ---*/
   if(CH_cmp(ga_info, c1, c2) > 0) {
      if(CH_cmp(ga_info, c1, p2) > 0) {
         xp1 = c1->xp1; xp2 = c1->xp2;
         CH_copy(p2, c1);
         c1->xp1 = xp1; c1->xp2 = xp2;
      }
   } else {
      if(CH_cmp(ga_info, c2, p2) > 0) {
         xp1 = c2->xp1; xp2 = c2->xp2;
         CH_copy(p2, c2);
         c2->xp1 = xp1; c2->xp2 = xp2;
      }
   }

   /*--- Update indices ---*/
   c1->parent_1 = p1->index;
   c1->parent_2 = p2->index;
   c2->parent_1 = p1->index;
   c2->parent_2 = p2->index;
}
/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| GA reports
|
| Functions:
|    RP_report() - periodic report function
|    RP_config() - configuration report function
|    RP_final()  - final report function
|    RP_time()   - time for a report?
|    RP_short()  - short report format
|    RP_long()   - long report format
============================================================================*/

/*----------------------------------------------------------------------------
| Print a report
----------------------------------------------------------------------------*/
RP_report(
   GA_Info_Ptr    ga_info,
   Pool_Ptr       pool)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RP_report: invalid ga_info");
   if(!PL_valid(pool)) UT_error("RP_report: invalid pool");

   /*--- Is it report time? ---*/
   if(!RP_time(ga_info, pool)) return GA_ERROR;  // CCC

   /*--- Yes, do appropriate report ---*/
   switch(ga_info->rp_type) {
      case RP_NONE:    break;
      case RP_MINIMAL: break;
      case RP_SHORT:   RP_short(ga_info, pool);
                       break;
      case RP_LONG:    RP_long(ga_info, pool);
                       break;
   }
}

/*----------------------------------------------------------------------------
| Print a configuration report
----------------------------------------------------------------------------*/
RP_config(   GA_Info_Ptr    ga_info)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RP_config: invalid ga_info");

   switch(ga_info->rp_type) {

      /*--- These do not have a config report ---*/
      case RP_NONE:    
                       break;

      /*--- All others have a config report ---*/
      default:
                       CF_report(ga_info);
                       break;
   }
}

/*----------------------------------------------------------------------------
| Print a final report
----------------------------------------------------------------------------*/
RP_final( GA_Info_Ptr    ga_info)
{
   int i;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RP_final: invalid ga_info");
   if(ga_info->rp_fid == NULL) UT_error("RP_final: invalid ga_info->rp_fid");

   /*--- Skip final report if RP_NONE ---*/
   if(ga_info->rp_type == RP_NONE) return OK; // CCC

   /*--- Reason for stopping ---*/
   if(ga_info->use_convergence && ga_info->converged) {
      fprintf(ga_info->rp_fid,
              "\nThe GA has converged after %d iterations.\n",
              ga_info->iter);
   } else {
      fprintf(ga_info->rp_fid,
              "\nThe specified number of iterations has been reached.\n");
   }

   /*--- Print best ---*/
   fprintf(ga_info->rp_fid,"\nBest: ");
   for(i = 0; i < ga_info->best->length; i++) {
      fprintf(ga_info->rp_fid,"%G ", ga_info->best->gene[i]);
      if(i % 20 == 19 && i+1 < ga_info->best->length) 
         fprintf(ga_info->rp_fid,"\n      ");
   }
   fprintf(ga_info->rp_fid," (%g)\n\n", ga_info->best->fitness);
}

/*----------------------------------------------------------------------------
| See if time to print a report
----------------------------------------------------------------------------*/
RP_time(
   GA_Info_Ptr    ga_info,
   Pool_Ptr       pool)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RP_time: invalid ga_info");
   if(!PL_valid(pool)) UT_error("RP_time: invalid pool");

   /*--- First generation ---*/
   if(ga_info->iter == 0) return TRUE; 

   /*--- Report interval reached ---*/
   if(!((ga_info->iter+1) % ga_info->rp_interval)) return TRUE;

   /*--- Last generation reached ---*/
   if(ga_info->iter+1 == ga_info->max_iter) return TRUE;

   /*--- GA has converged ---*/
   if(ga_info->use_convergence && ga_info->converged) return TRUE;

   /*--- Otherwise, not time for report ---*/
   return FALSE;
}

/*----------------------------------------------------------------------------
| Print a short report
----------------------------------------------------------------------------*/
RP_short(GA_Info_Ptr    ga_info,
   Pool_Ptr       pool)
{
   int i, j;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RP_short: invalid ga_info");
   if(!PL_valid(pool)) UT_error("RP_short: invalid pool");
   if(ga_info->rp_fid == NULL) UT_error("RP_short: invalid ga_info->rp_fid");

   /*--- Print header first time only ---*/
   if(ga_info->iter < 0) {
      fprintf(ga_info->rp_fid,"\n%s%s\n%s%s\n",
         "Gener    Min      Max      Ave    Variance  ",
         "Std Dev  Tot Fit    Best ",
         "-----  -------  -------  -------  --------  ",
         "-------  -------  -------"
      );
   }

   /*--- Print a line for current iteration ---*/
   fprintf(ga_info->rp_fid,
      "%5d  %7.6G  %7.6G  %7.3G  %8.3G  %7.3G  %7.6G  %7.6G\n", 
      ga_info->iter+1, pool->min, pool->max, pool->ave, pool->var, pool->dev,
      pool->total_fitness, ga_info->best->fitness);
   fflush(ga_info->rp_fid);

   if(ga_info->iter+1 == ga_info->max_iter) {
      /*--- Statistics ---*/
      fprintf(ga_info->rp_fid,
         "\nMin= %G   Max= %G   Ave= %.2G   Tot= %G   Var= %.2G   SD= %.2G\n", 
         pool->min, pool->max, pool->ave, pool->total_fitness, 
         pool->var, pool->dev);

      /*--- Print best ---*/
      fprintf(ga_info->rp_fid,"\nBest: ");
      for(i = 0; i < ga_info->best->length; i++) {
         fprintf(ga_info->rp_fid,"%G ", ga_info->best->gene[i]);
         if(i % 20 == 19 && i+1 < ga_info->best->length) 
            fprintf(ga_info->rp_fid,"\n      ");
      }
      fprintf(ga_info->rp_fid,"(%G)\n", ga_info->best->fitness);

      /*--- Report footer line ---*/
      fprintf(ga_info->rp_fid,"=======================================");
      fprintf(ga_info->rp_fid,"=======================================");
      fprintf(ga_info->rp_fid,"\n");
      fflush(ga_info->rp_fid);
   }
}

/*----------------------------------------------------------------------------
| print a long report
----------------------------------------------------------------------------*/
RP_long(GA_Info_Ptr    ga_info,
   Pool_Ptr       pool)
{
   int i, j;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("RP_long: invalid ga_info");
   if(!PL_valid(pool)) UT_error("RP_long: invalid pool");
   if(ga_info->rp_fid == NULL) UT_error("RP_long: invalid ga_info->rp_fid");

   /*--- Report header line ---*/
   fprintf(ga_info->rp_fid,"\n");
   fprintf(ga_info->rp_fid,"=======================================");
   fprintf(ga_info->rp_fid,"=======================================");
   fprintf(ga_info->rp_fid,"\n");

   /*--- Generation and # mutations ---*/
   fprintf(ga_info->rp_fid,"Generation %d: Mutations = %d (%d total)\n\n", 
      ga_info->iter+1, ga_info->num_mut, ga_info->tot_mut);

   /*--- Pool header ---*/
   fprintf(ga_info->rp_fid," # Parents  XP   Fitness  String\n");
   fprintf(ga_info->rp_fid,"-- ------- ----- -------  ------\n");

   /*--- Print pool ---*/
   for(i = 0; i < pool->size; i++) {
     fprintf(ga_info->rp_fid,"%2d (%2d,%2d) %2d %2d %7G  ",
        i+1, pool->chrom[i]->parent_1 + 1, pool->chrom[i]->parent_2 + 1, 
        pool->chrom[i]->xp1 + 1, pool->chrom[i]->xp2 + 1, 
        pool->chrom[i]->fitness);
     for(j = 0; j < pool->chrom[i]->length; j++) {
        fprintf(ga_info->rp_fid,"%G ", pool->chrom[i]->gene[j]);
        if(j % 15 == 14 && j+1 < pool->chrom[i]->length) 
           fprintf(ga_info->rp_fid,"\n                                  ");
     }
     fprintf(ga_info->rp_fid,"\n");
   }

   /*--- Statistics ---*/
   fprintf(ga_info->rp_fid,
      "\nMin= %G   Max= %G   Ave= %.2G   Tot= %G   Var= %.2G   SD= %.2G\n", 
      pool->min, pool->max, pool->ave, pool->total_fitness, 
      pool->var, pool->dev);

   /*--- Print best ---*/
   fprintf(ga_info->rp_fid,"\nBest: ");
   for(i = 0; i < ga_info->best->length; i++) {
      fprintf(ga_info->rp_fid,"%G ", ga_info->best->gene[i]);
      if(i % 20 == 19 && i+1 < ga_info->best->length) 
         fprintf(ga_info->rp_fid,"\n      ");
   }
   fprintf(ga_info->rp_fid,"(%G)\n", ga_info->best->fitness);

   /*--- Report footer line ---*/
   fprintf(ga_info->rp_fid,"=======================================");
   fprintf(ga_info->rp_fid,"=======================================");
   fprintf(ga_info->rp_fid,"\n");
   fflush(ga_info->rp_fid);

   return(OK);
}
/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
| (c) Copyright IA UPM - Group 5, 2020.  All rights reserved.
|
| Selection operators 
|
| Operators
|    SE_uniform_random()  - just pick one
|    SE_roulette()        - standard roulette
|       SE_max_roulette() - helper for roulette (maximizing)
|       SE_min_roulette() - helper for roulette (minimizing)
|    SE_rank_biased()     - standard linear bias 
|
| Interface
|    SE_table[]   - used in selection of selection method
|    SE_set_fun() - set and select user defined selection function
|    SE_select()  - select selection method by name
|    SE_name()    - get name of current selection method
|    SE_fun()     - setup and perform selection operator
============================================================================*/

/*============================================================================
|                     Selection interface
============================================================================*/
/*----------------------------------------------------------------------------
| Selection table
----------------------------------------------------------------------------*/
FN_Table_Type SE_table[] = {
   { NULL,             NULL              },  /* user defined function */
   { "uniform_random", SE_uniform_random },
   { "roulette",       SE_roulette       },
   { "rank_biased",    SE_rank_biased    },
   { NULL,             NULL              }
};

/*----------------------------------------------------------------------------
| Select user selection method
----------------------------------------------------------------------------*/
SE_set_fun(
   GA_Info_Ptr  ga_info,
   char         *fn_name,
   FN_Ptr       fn_ptr)
{
   return FN_set_fun(ga_info, SE_table, fn_name, fn_ptr, &ga_info->SE_fun);
}

/*----------------------------------------------------------------------------
| Select selection method
----------------------------------------------------------------------------*/
SE_select(
   GA_Info_Ptr  ga_info,
   char         *fn_name)
{
   return FN_select(ga_info, SE_table, fn_name, &ga_info->SE_fun);
}

/*----------------------------------------------------------------------------
| Selection name
----------------------------------------------------------------------------*/
char *SE_name(   GA_Info_Ptr  ga_info)
{
   return FN_name(ga_info, SE_table, ga_info->SE_fun);
}

/*----------------------------------------------------------------------------
| Selection interface
----------------------------------------------------------------------------*/
Chrom_Ptr SE_fun(
   GA_Info_Ptr    ga_info,
   Pool_Ptr       pool)
{
   int       idx;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("SE_fun: invalid ga_info");

   if(ga_info->SE_fun == NULL) UT_error("SE_fun: null SE_fun");

   /*--- Select a chromosome ---*/
   idx = ga_info->SE_fun(ga_info, pool);
   if(idx < 0 || idx >= pool->size) UT_error("SE_fun: invalid idx");
   if(pool->chrom[idx] == NULL) UT_error("SE_fun: null pool->chrom[idx]");

   return pool->chrom[idx];
}

/*============================================================================
|                             Selection Methods
============================================================================*/
/*----------------------------------------------------------------------------
| Uniform random
----------------------------------------------------------------------------*/
SE_uniform_random(GA_Info_Ptr    ga_info,
   Pool_Ptr       pool)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("SE_uniform_random: invalid ga_info");

   return RAND_DOM(0, pool->size-1);
}

/*----------------------------------------------------------------------------
| Roulette
----------------------------------------------------------------------------*/
SE_roulette(GA_Info_Ptr    ga_info,
   Pool_Ptr       pool)
{
   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("SE_roulette: invalid ga_info");

   /*--- Find PTF for each chromosome ---*/                         
   PL_update_ptf(ga_info, pool);

   if(ga_info->minimize)         
      return SE_min_roulette(ga_info, pool);
   else
      return SE_max_roulette(ga_info, pool);
}      

/*----------------------------------------------------------------------------
| Roulette helper (maximizing)
----------------------------------------------------------------------------*/
SE_max_roulette(GA_Info_Ptr    ga_info,   Pool_Ptr       pool)
{
   int   i = 0;
   float val = 0.0, spin_val;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("SE_max_roulette: invalid ga_info");

   /*--- Spin the wheel ---*/
   spin_val = RAND_FRAC() * pool->total_fitness;

   /*--- Find corresponding chromosome ---*/
   while(val < spin_val && i < pool->size)
      val += pool->chrom[i++]->fitness;
   if(i > 0) i--;

   return i;
}

/*----------------------------------------------------------------------------
| Roulette helper (minimizing)
----------------------------------------------------------------------------*/
SE_min_roulette(GA_Info_Ptr    ga_info,
   Pool_Ptr       pool)
{
   int   i = 0;
   float val = 0.0, spin_val;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("SE_min_roulette: invalid ga_info");

   /*--- Spin the wheel (value between 0.0 and 100.0) ---*/
   spin_val = RAND_FRAC() * 100.0;

   /*--- Find corresponding chromosome ---*/
   while(val < spin_val && i < pool->size)
      val += pool->chrom[i++]->ptf;
   if(i > 0) i--;

   return i;
}

/*----------------------------------------------------------------------------
| Rank biased 
----------------------------------------------------------------------------*/
SE_rank_biased(GA_Info_Ptr    ga_info,
   Pool_Ptr       pool)
{
   static int ranked = FALSE;

   /*--- Error check ---*/
   if(!CF_valid(ga_info)) UT_error("SE_rank_biased: invalid ga_info");

   /*--- Rank pool ---*/
   if(!ranked) {
      PL_sort(ga_info, pool);

      /*--- Only rank once if replacement is by_rank ---*/
      if(!strcmp(RE_name(ga_info), "by_rank")) ranked = TRUE;
   }

   /*--- Linear biased selection ---*/
   return pool->size * (ga_info->bias - sqrt(ga_info->bias * ga_info->bias
          - 4.0 * (ga_info->bias-1) * RAND_FRAC())) / 2.0 / (ga_info->bias-1);
}


///  added by claudio 11/5/2005

//*--- Re-evaluate best chromosome ---*
void re_evaluate_pop(GA_Info_Ptr ga_info)
{
	float ZZZ=9999999.9;
	int i;

	for(i=0;i<ga_info->pool_size;i++)
	{
		ga_info->EV_fun(ga_info->old_pool->chrom[i]);
    //obj_fun(ga_info->old_pool->chrom[i]);
	    if(ga_info->old_pool->chrom[i]->fitness<ZZZ)
		{
		  ZZZ=ga_info->old_pool->chrom[i]->fitness;
		  ga_info->best=ga_info->old_pool->chrom[i];
		}
    }
    
}
