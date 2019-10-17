/*!
  \file caa_estimate.c
  \brief Containing the main routines for fitting the age, lga and wgl model.
  \author Geir Storvik and Hanne Rognebakke

  Fitting the age and lga models are performed in the ::main_model1 routine.

  Fitting the wgl model is performed in the ::main_model2 routine.
  
*/
#include "caa.h"
#include "caa_read_write.h"
#include "caa_mcmc.h"
#include "caa_routines.h"
#include "caa_cell_constr.h"
#include "caa_init.h"
#include "caa_evaluate.h" 
#include "caa_sample_gauss.h"
#include "caa_input.h"
#include "caa_estimate.h"

#ifdef LOG_FILE
extern FILE     *g_caa_log; 
#endif

extern FILE     *g_caa_mcmc1; 
extern FILE     *g_caa_mcmc2; 
extern FILE     *g_caa_mcmc_hsz; 
extern FILE     *g_caa_mcmc_hsz_eff;



/*!
  \brief Fitting age and lga model through MCMC simulations.
  \author Geir Storvik and Hanne Rognebakke

  The routine starts to convert input data into 
  approperiate c-structures (as defined in caa.h), 
  and does the necessary initialization.

  Thereafter call the main routine ::MCMC_model1 for performing the MCMC simulations.
*/
int main_model1(Input_common *i_inCommon, Input_age *i_inAge, Input_lga *i_inLga, 
		Input_wgl *i_inHsz, Input_prior *i_inPrior, Data_orig *i_D_orig, 
		Data_CC *i_D_CC)
{
  /* Data structures */
  Age_struct       *age;             /* Current simulations of age-parameters */ 
  Age_struct       *age_mean;        /* Mean of simulations of age-parameters */ 
  LW_struct        *length;          /* Current simulations of lga-parameters */ 
  LW_struct        *length_mean;     /* Mean of simulations of lga-parameters */ 
  Data_age         *D_age;           /* Age data */
  Data_lin         *D_lga;           /* Lga data */
  Data_g_a         *D_g_a;           /* g_a parameters and simulations */
  LW_struct        *length_CC=NULL;  /* Current simulations of lga-parameters for coastal cod */ 
  LW_struct        *length_CC_mean=NULL; /* Mean of simulations of lga-parameters for coastal cod */ 
  Data_lin         *D_lga_CC=NULL;   /* Lga data for coastal cod */
  Data_g_a         *D_g_a_CC=NULL;   /* g_a parameters and simulations for coastal cod */

  int       err=0;
  int       i, l_int, max_boat;
  double    lobs;
  char      buffer[MAX_STR];
  FILE     *fp;
  double   *mod1_mean_inv_lik; 
  double   *age_mean_inv_lik;
  double   *lga_mean_inv_lik;

  mod1_mean_inv_lik = CALLOC(i_D_orig->nHaul,double); 
  age_mean_inv_lik = CALLOC(i_D_orig->nHaul,double);
  lga_mean_inv_lik = CALLOC(i_D_orig->nHaul,double);
 
  #ifdef LOG_FILE
  long      time_now, time_start;
  time_start = clock();
  time_now = clock();
  fprintf(g_caa_log," CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"\nStart main_model1\n");
  #endif

  #ifdef DEBUG_INPUT
  sprintf(buffer,"%s/caa_input_model1.txt",i_inCommon->inputfolder);
  fprintf(stderr,"Write input model1 to file: %s\n",buffer);
  fp = fopen(buffer,"w");   
  err = write_input_model1(fp,i_D_orig,i_inCommon,i_inAge,i_inLga,i_inPrior,i_D_CC);
  fclose(fp);
  if(err)
    {
      write_warning("main_model1:Error calling write_input_model1\n");
      return(err);
    }
  #endif

  if(i_inCommon->inc_hsz)
    {
      #ifdef DEBUG_PROG
      fprintf(stderr,"Haulsize included in model, do haulsize regression\n");
      #endif
      err = main_model_hsz(i_inCommon,i_inHsz,i_inPrior,i_D_orig);
    }

  /* Make age data */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Make age data: makedata_age1\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log," CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"\nStart makedata_age1\n");
  #endif
  err = makedata_age1(i_D_orig->nHaul,i_inAge->nAges,i_inAge->a_vec,
		      i_inAge->cov,i_inAge->num_adj_area,i_inAge->adj_area,
		      &D_age,1);
  if(err)
    {
      write_warning("main_model1:Error calling makedata_age1\n");
      return(err);
    } 
  /* Put boat covariate in D_orig*/ // moved from add_object_info_age_lga because iboat not in input anymore - must be initialized in makedata_age1
  max_boat = 1;
  for(i=0;i<i_D_orig->nHaul;i++)
    {
      i_D_orig->boat[i] = 1;
      if(D_age->glm->xcov[0]->iboat > -1)
	{
 	  i_D_orig->boat[i] = D_age->glm->xcov[0]->c_cov[i][D_age->glm->xcov[0]->iboat];
	  max_boat = MAX(max_boat,i_D_orig->boat[i]);
	}
    }
  i_D_orig->nBoat = max_boat;

  #ifdef DEBUG_PROG
  fprintf(stderr,"Make age data: makedata_age2\n");
  #endif
  err = makedata_age2(i_D_orig,D_age,i_D_CC->class_error);
  if(err)
    {
      write_warning("main_model1:Error calling makedata_age2\n");
      return(err);
    }

  // Construct length interval limits - maybe these should be a part of the original input from R?
  i_D_orig->lstart = CALLOC(i_D_orig->nFish,double); // FREE in end of main_model1
  i_D_orig->lend = CALLOC(i_D_orig->nFish,double);
  for(i=0;i<i_D_orig->nFish;i++)
    {
      lobs = i_D_orig->totlength[i];
      l_int = 0;
      while(lobs > i_D_orig->int_len_lim[l_int])
	l_int++;
      i_D_orig->lstart[i] = i_D_orig->int_len_lim[l_int-1];      
      i_D_orig->lend[i] = i_D_orig->int_len_lim[l_int];      
    }


  /* Make g_a data */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Make g_a data: makedata_g_a\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log," CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"\nStart makedata_g_a\n");
  #endif
  err = makedata_g_a(D_age->glm->ncat,i_inLga->g_a_ncat,i_inLga->g_a_nSeason,
		     D_age->a_vec,i_D_orig->coastal_cod,
		     i_inLga->g_a_model,i_inLga->g_a_par_init,i_inLga->g_a_sample_c,
		     i_inLga->g_a_sample_theta,i_inLga->g_a_sample_gamma,
		     &D_g_a);
  if(err)
    {
      write_warning("main_model1:Error calling makedata_g_a\n");
      return(err);
    }     
  if(i_D_orig->coastal_cod)
    {
      /* Make extra g_a data if coastal cod */
      err = makedata_g_a(D_age->glm->ncat,i_inLga->g_a_ncat,i_inLga->g_a_nSeason,
			 D_age->a_vec,i_D_orig->coastal_cod,
			 i_inLga->g_a_model,i_inLga->g_a_par_init,i_inLga->g_a_sample_c,
			 i_inLga->g_a_sample_theta,i_inLga->g_a_sample_gamma,
			 &D_g_a_CC);
      if(err)
	{
	  write_warning("main_model1:Error calling makedata_g_a\n");
	  return(err);
	}
    }
  
  /* Make lga data */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Make lga data: makedata_lin1\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log," CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"\nStart makedata_lin1\n");
  #endif
  err = makedata_lin1(i_D_orig->nHaul,i_D_orig->haulweight,
		      i_inLga->int_cov,i_inLga->slp_cov,i_inLga->num_adj_area,i_inLga->adj_area,
		      &D_lga,1,0);
  if(err)
    {
      write_warning("main_model1:Error calling makedata_lin1\n");
      return(err);
    }

  /* Make extra lga data if coastal cod */
  if(i_D_orig->coastal_cod)
    {
      #ifdef LOG_FILE
      time_now = clock();
      fprintf(g_caa_log," CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
      fprintf(g_caa_log,"\nStart makedata_lin1 coastal cod\n");
      #endif
      err = makedata_lin1(i_D_orig->nHaul,i_D_orig->haulweight,
			  i_inLga->int_cov,i_inLga->slp_cov,
			  i_inLga->num_adj_area,i_inLga->adj_area,
			  &D_lga_CC,1,0);
      if(err)
	{
	  write_warning("main_model1:Error calling makedata_lin1\n");
	  return(err);
	}
    }
  #ifdef DEBUG_PROG
  fprintf(stderr,"Make lga data: makedata_lga_suff\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"Start makedata_lga_suff\n");
  #endif
  if(i_D_orig->coastal_cod==0)
    {
      err = makedata_lga_suff(D_lga,i_D_orig,D_g_a);
      if(err)
	{
	  write_warning("main_model1:Error calling makedata_lga_suff\n");
	  return(err);
	}
    }
  else
    {
      err = makedata_lga_suff_CC(D_lga,D_lga_CC,i_D_orig,D_g_a,D_g_a_CC,i_D_CC->class_error);
      if(err)
	{
	  write_warning("main_model1:Error calling makedata_lga_suff_CC\n");
	  return(err);
	}
    }
  
  #ifdef DEBUG_PROG
  fprintf(stderr,"Initialize model1\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"Start initialize_model1\n");
  #endif
  err = initialize_model1(i_inCommon->seed,i_inCommon->num_par1,
			  D_age,&age,&age_mean,
			  D_lga,&length,&length_mean,
			  D_lga_CC,&length_CC,&length_CC_mean,
			  D_g_a,D_g_a_CC,
			  i_inAge->errors,i_inAge->A2A,i_inAge->delta_age,
			  i_D_orig->coastal_cod,
			  i_inCommon->inc_hsz,i_inCommon->sim_ar,
			  i_inPrior->age_eff_mean,i_inPrior->age_eff_prec,
			  i_inPrior->age_prec_par,i_inPrior->age_ar,
			  i_inPrior->lga_eff_mean,i_inPrior->lga_eff_prec,
			  i_inPrior->lga_prec_par,i_inPrior->lga_ar,
			  i_inLga->fixed_model,i_inLga->fixed_int,
			  i_inLga->fixed_slp,i_inLga->fixed_tau,
			  i_inLga->fixed_g_a_c,i_inLga->fixed_g_a_theta,
			  i_inLga->fixed_g_a_gamma,i_inLga->fixed_int_CC,
			  i_inLga->fixed_slp_CC,i_inLga->fixed_tau_CC,
			  i_inLga->fixed_g_a_c_CC,i_inLga->fixed_g_a_theta_CC,
			  i_inLga->fixed_g_a_gamma_CC, 
			  i_inLga->cens_model,i_inLga->cens_par);
  if(err)
    {
      write_warning("main_model1:Error calling initialize_model1\n");
      return(err);
    }


  /* Find node numbers in fixed effect graph for age */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Find node numbers: find_node_effect\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"Start find_node_effect\n");
  #endif
  D_age->glm->xcov[0]->n_cov--; //not include haul effect here
  err = find_node_effect(D_age->glm->xcov,D_age->glm->nxcov,age->gr_str_f->in_gr,
			 &(age->gr_str_f->node));
  D_age->glm->xcov[0]->n_cov++;
  if(err)
    {
      write_warning("main_model1:Error calling find_node_effect\n");
      return(err);
    }

  /* Find node numbers in graph for length */
  err = find_node_effect(D_lga->glm->xcov,D_lga->glm->nxcov,length->gr_str->in_gr,
			 &(length->gr_str->node));
  if(err)
    {
      write_warning("main_model1:Error calling find_node_effect\n");
      return(err);
    }
  if(i_D_orig->coastal_cod)
    {
      err = find_node_effect(D_lga_CC->glm->xcov,D_lga_CC->glm->nxcov,length_CC->gr_str->in_gr,
			     &(length_CC->gr_str->node));
      if(err)
	{
	  write_warning("main_model1:Error calling find_node_effect\n");
	  return(err);
	}
    }

  /* Initialize GMRFLib-type graphs for age and lga models */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Initialize GMRFLib-type graphs: init_graph_model1\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"Start init_graph_model1\n");
  #endif
  err = init_graph_model1(D_age,age,D_lga,length,D_lga_CC,length_CC,
			  i_D_orig->coastal_cod,i_inCommon->constr);
  if(err)
    {
      write_warning("main_model1:Error calling init_graph_model1\n");
      return(err);
    }

  /* Starting values of MCMC simulations */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Starting values of MCMC simulations: MCMC_model1_init\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"Start MCMC_model1_init\n");
  #endif

  #ifdef LOG_FILE
  fprintf(g_caa_log,"use MCMC_model1_init\n");
  #endif
  err = MCMC_model1_init(i_D_orig,D_age,age,D_lga,length,D_lga_CC,length_CC,D_g_a,D_g_a_CC,
			 i_D_CC,NULL,i_inCommon->use_debug);
  if(err)
    {
      write_warning("main_model1:Error calling MCMC_model1_init\n");
      return(err);
    }

  /* Open file for printing mcmc parameters in binary format */
  if(i_inCommon->print_format==0)
    {
      fprintf(stderr,"Print parameters to new file  %s  in binary format\n",i_inCommon->filename_mcmc1);
      g_caa_mcmc1 = fopen(i_inCommon->filename_mcmc1, "wb");
    }
  else
    {
      fprintf(stderr,"Print parameters to new file  %s  in ascii format\n",i_inCommon->filename_mcmc1);
      g_caa_mcmc1 = fopen(i_inCommon->filename_mcmc1, "w");
    }
  if(!(g_caa_mcmc1))
    {
      sprintf(buffer,"main_model1: Couldn't open file for writing: %s",i_inCommon->filename_mcmc1);
      write_warning(buffer);
      return(1);
    }
  #ifdef LOG_FILE
  fprintf(g_caa_log,"print model parameters to file %s\n",i_inCommon->filename_mcmc1);
  #endif
  err = write_mcmc1(D_age, age, D_lga, length, D_lga_CC, length_CC,
		    D_g_a, D_g_a_CC, i_inCommon->num_it_outer, 
		    i_inCommon->num_par1, i_D_orig->coastal_cod,
		    i_inCommon->print_boat, age->delta_age, i_inCommon->print_format);
  if(err)
    {
      write_warning("main_model1:Error calling write_mcmc1\n");
      return(err);
    }


  /* Start MCMC simulations */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Start MCMC simulations: MCMC_model1\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"Start MCMC_model1\n");
  #endif
  err = MCMC_model1(i_inCommon,i_D_orig,D_age,age,age_mean,
		    D_lga,length,length_mean,D_lga_CC,length_CC,length_CC_mean,
		    D_g_a,D_g_a_CC,i_D_CC,mod1_mean_inv_lik,age_mean_inv_lik,lga_mean_inv_lik,
		    NULL);
  if(err)
    {
      write_warning("main_model1:Error calling MCMC_model1\n");
      return(err);
    }
  
  /* Close file with mcmc samples */    
  err = fclose(g_caa_mcmc1);
  printf("close mcmc1 file: %d\n",err);
  
  /* Clean up */
  #ifdef DEBUG_PROG
  fprintf(stderr,"\nClean up\n");
  #endif
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"Start cleaning up\n");
  #endif
if(0){  
  
  FREE(mod1_mean_inv_lik);
  FREE(age_mean_inv_lik);
  FREE(lga_mean_inv_lik);
 
  err = re_init_graph_model1(D_age,age,length,length_CC,i_D_orig->coastal_cod);
  if(err)
    {
      write_warning("main_model1:Error calling re_init_graph_model1\n");
      return(err);
    }

  D_age->glm->xcov[0]->n_cov--;
  err = re_find_node_effect(D_age->glm->xcov,D_age->glm->nxcov,age->gr_str_f->in_gr,
			    &(age->gr_str_f->node));
  D_age->glm->xcov[0]->n_cov++;
  if(err)
    {
      write_warning("main_model1:Error calling re_find_node_effect\n");
      return(err);
    }

  err = re_find_node_effect(D_lga->glm->xcov,D_lga->glm->nxcov,length->gr_str->in_gr,
			    &(length->gr_str->node));
  if(i_D_orig->coastal_cod)
    err = re_find_node_effect(D_lga_CC->glm->xcov,D_lga_CC->glm->nxcov,length_CC->gr_str->in_gr,
			      &(length_CC->gr_str->node));
  if(err)
    {
      write_warning("main_model1:Error calling re_find_node_effect\n");
      return(err);
    }

  err = re_initialize_model1(D_age,age,age_mean,D_lga,length,length_mean,
			     D_lga_CC,length_CC,length_CC_mean,D_g_a,
			     i_D_orig->coastal_cod,i_inAge->errors,i_inAge->A2A,i_inCommon->inc_hsz);
  if(err)
    {
      write_warning("main_model1:Error calling re_initialize\n");
      return(err);
    }

  if(i_D_orig->coastal_cod==0)
    {
      err = re_makedata_lga_suff(D_lga);
      if(err)
	{
	  write_warning("main_model1:Error calling re_makedata_lga_suff\n");
	  return(err);
	}
    }
  else
    {
      err = re_makedata_lga_suff_CC(D_lga,D_lga_CC);
      if(err)
	{
	  write_warning("main_model1:Error calling re_makedata_lga_suff_CC\n");
	  return(err);
	}
    }

  err = re_makedata_lin1(i_inLga->num_adj_area,i_inLga->adj_area,
			 &D_lga,1);
  if(err)
    {
      write_warning("main_model1:Error calling re_makedata_lin1\n");
      return(err);
    }
 if(i_D_orig->coastal_cod)
    {
      err = re_makedata_lin1(i_inLga->num_adj_area,i_inLga->adj_area,
			     &D_lga_CC,1);
      if(err)
	{
	  write_warning("main_model1:Error calling re_makedata_lin1\n");
	  return(err);
	}
    }

  err = re_makedata_g_a(&D_g_a);
  if(err)
    {
      write_warning("main_model1:Error calling re_makedata_g_a\n");
      return(err);
    }
  if(i_D_orig->coastal_cod)
    {
      err = re_makedata_g_a(&D_g_a_CC);
      if(err)
	{
	  write_warning("main_model1:Error calling re_makedata_g_a\n");
	  return(err);
	}
    } 
 
  err = re_makedata_age2(D_age);
  if(err)
    {
      write_warning("main_model1:Error calling re_makedata_age2\n");
      return(err);
    }

  err = re_makedata_age1(i_inAge->num_adj_area,i_inAge->adj_area,
			 &D_age,1);
  if(err)
    {
      write_warning("main_model1:Error calling re_makedata_age1\n");
      return(err);
    }
  
  FREE(i_D_orig->lstart);
  FREE(i_D_orig->lend);

}
  
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  #endif


  return(0);
}		/* end of main_model1 */




/*!
  \brief Fitting wgl model through MCMC simulations.
  \author Geir Storvik and Hanne Rognebakke

  The routine starts to convert input data into 
  approperiate c-structures (as defined in caa.h), 
  and does the necessary initialization.

  Thereafter call the main routine ::MCMC_model2 for performing the MCMC simulations.
*/
int main_model2(Input_common *i_inCommon, Input_wgl *i_inWgl,
		Input_prior *i_inPrior, Data_orig *i_D_orig, int i_hsz)
{
  /* Data structures */
  LW_struct      *weight;          /*!< Current simulations of wgl-parameters */ 
  LW_struct      *weight_mean;     /*!< Mean of simulations of wgl-parameters */ 
  Data_lin       *D_wgl;           /*!< Wgl data */
  LW_struct      *weight_CC=NULL;  /*!< Current simulations of wgl-parameters for coastal cod */ 
  LW_struct      *weight_CC_mean=NULL; /*!< Mean of simulations of wgl-parameters for coastal cod */ 
  Data_lin       *D_wgl_CC=NULL;   /*!< Wgl data for coastal cod */

  int       err=0;
  char      buffer[MAX_STR];
  FILE     *fp;

  double   *wgl_mean_inv_lik;

  wgl_mean_inv_lik = CALLOC(i_D_orig->nHaul,double);

  #ifdef LOG_FILE
  long      time_now, time_start;
  time_start = clock();
  time_now = clock();
  fprintf(g_caa_log," CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"\nStart main_model2\n");
  #endif

  if(i_hsz==0)
    i_inCommon->inc_hsz = 0;

  #ifdef DEBUG_INPUT
  if(i_inCommon->inc_hsz) 
    {
      sprintf(buffer,"%s/caa_input_model2_hsz.txt",i_inCommon->inputfolder);
      fprintf(stderr,"Write input model2 to file: %s\n",buffer);
      fp = fopen(buffer,"w");
    } 
  else 
    {
      sprintf(buffer,"%s/caa_input_model2.txt",i_inCommon->inputfolder);
      fprintf(stderr,"Write input model2 to file: %s\n",buffer);
      fp = fopen(buffer,"w");
    }
  err = write_input_model2(fp,i_D_orig,i_inCommon,i_inWgl,i_inPrior);
  fclose(fp);
  if(err)
    {
      write_warning("main_model2:Error calling write_input_model2\n");
      return(err);
    } 
  #endif 

  /* Make wgl data */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Make wgl data: makedata_lin1\n");
  #endif
  err = makedata_lin1(i_D_orig->nHaul,i_D_orig->haulweight,
		      i_inWgl->int_cov,i_inWgl->slp_cov,
		      i_inWgl->num_adj_area,i_inWgl->adj_area,
		      &D_wgl,1,i_hsz);
  if(err)
    {
      write_warning("main_model2:Error calling makedata_lin1\n");
      return(err);
    }

  /* Make extra wgl data if coastal cod */
  if(i_D_orig->coastal_cod)
    {
      err = makedata_lin1(i_D_orig->nHaul,i_D_orig->haulweight,
			  i_inWgl->int_cov,i_inWgl->slp_cov,
			  i_inWgl->num_adj_area,i_inWgl->adj_area,
			  &D_wgl_CC,1,0);
      if(err)
	{
	  write_warning("main_model2:Error calling makedata_lin1\n");
	  return(err);
	}
    }

  #ifdef DEBUG_PROG
  fprintf(stderr,"Make wgl data: sufficient statistics\n");
  #endif
  if(i_D_orig->coastal_cod==0)
    {
      err = makedata_wgl_suff(D_wgl,i_D_orig);
      if(err)
	{
	  write_warning("main_model2:Error calling makedata_wgl_suff\n");
	  return(err);
	}
    }
  else
    {
      err = makedata_wgl_suff_CC(D_wgl,D_wgl_CC,i_D_orig);
      if(err)
	{
	  write_warning("main_model2:Error calling makedata_wgl_suff_CC\n");
	  return(err);
	}
    }

  #ifdef DEBUG_PROG
  fprintf(stderr,"Initialize model2\n");
  #endif
  err = initialize_model2(i_inCommon->seed,i_inCommon->num_par2,
			  D_wgl,&weight,&weight_mean,D_wgl_CC,&weight_CC,&weight_CC_mean,
			  i_inCommon->sim_ar,
			  i_inPrior->wgl_eff_mean,i_inPrior->wgl_eff_prec,
			  i_inPrior->wgl_prec_par,i_inPrior->wgl_ar,
			  i_inWgl->fixed_model,i_inWgl->fixed_int,
			  i_inWgl->fixed_slp,i_inWgl->fixed_tau,
			  i_inWgl->fixed_int_CC,i_inWgl->fixed_slp_CC,
			  i_inWgl->fixed_tau_CC,i_D_orig->coastal_cod);
  if(err)
    {
      write_warning("main_model2:Error calling initialize_model2\n");
      return(err);
    }

  /* Find node numbers in graph for weight */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Find node numbers in graph for weight\n");
  #endif
  err = find_node_effect(D_wgl->glm->xcov,D_wgl->glm->nxcov,weight->gr_str->in_gr,
			 &(weight->gr_str->node));
  if(err)
    {
      write_warning("main_model2:Error calling find_node_effect\n");
      return(err);
    }
  if(i_D_orig->coastal_cod)
    {
      err = find_node_effect(D_wgl_CC->glm->xcov,D_wgl_CC->glm->nxcov,weight_CC->gr_str->in_gr,
			     &(weight_CC->gr_str->node));
      if(err)
	{
	  write_warning("main_model2:Error calling find_node_effect\n");
	  return(err);
	}
    }
  /* Initialize GMRFLib-type graphs for wgl models */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Initialize GMRFLib-type graphs for wgl model\n");
  #endif
  err = make_graph_gauss(D_wgl->glm,weight->gr_str,weight->par,0,i_inCommon->constr);
  if(err)
    {
      write_warning("main_model2:Error calling make_graph_gauss\n");
      return(err);
    }
  if(i_D_orig->coastal_cod)
    {
      err = make_graph_gauss(D_wgl_CC->glm,weight_CC->gr_str,weight_CC->par,0,i_inCommon->constr);
      if(err)
	{
	  write_warning("main_model2:Error calling make_graph_gauss\n");
	  return(err);
	}
    }

  /* Open file for printing mcmc parameters in binary format */
  if(i_inCommon->print_format==0)
    {
      fprintf(stderr,"print parameters to new file  %s  in binary format\n",i_inCommon->filename_mcmc2);
      g_caa_mcmc2 = fopen(i_inCommon->filename_mcmc2, "wb");
    }
  else 
    {
      fprintf(stderr,"print parameters to new file  %s  in ascii format\n",i_inCommon->filename_mcmc2);
      g_caa_mcmc2 = fopen(i_inCommon->filename_mcmc2, "w");
    }
  if(!(g_caa_mcmc2))
    printError("main_model2: Couldn't open file for writing",i_inCommon->filename_mcmc2);
  err = write_mcmc2(D_wgl, weight, D_wgl_CC, weight_CC, i_inCommon->num_it_outer, 
		    i_inCommon->num_par2, i_D_orig->coastal_cod, i_inCommon->print_boat, i_inCommon->print_format);
  if(err)
    {
      write_warning("main_model2:Error calling write_mcmc2\n");
      return(err);
    }

  if(i_hsz)
    {
      /* Open file for printing mcmc parameters in binary format */
      if(!(g_caa_mcmc_hsz = fopen(i_inCommon->filename_hsz_it, "wb")))
	printError("main_model2: Couldn't open file for writing",i_inCommon->filename_hsz_it);
      fprintf(stderr,"print specific haulsize parameters to file %s\n",i_inCommon->filename_hsz_it);
      if(!(g_caa_mcmc_hsz_eff = fopen(i_inCommon->filename_hsz_hauleff, "wb")))
	printError("main_model2: Couldn't open file for writing",i_inCommon->filename_hsz_hauleff);
      fprintf(stderr,"print haulsize effects to file %s\n",i_inCommon->filename_hsz_hauleff);
    }

  /* Run MCMC */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Run MCMC\n");
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Run MCMC\n");
  #endif
  err = MCMC_model2(i_inCommon, 
		    D_wgl,weight,weight_mean,D_wgl_CC,weight_CC,weight_CC_mean,
		    i_D_orig->coastal_cod,wgl_mean_inv_lik);
  if(err)
    {
      write_warning("main_model2:Error calling MCMC_model2\n");
      return(err);
    }

 
  /* Close file with mcmc samples */    
  err = fclose(g_caa_mcmc2);
  printf("close mcmc2 file: %d\n",err);
  if(i_hsz)
    {
      fclose(g_caa_mcmc_hsz);
      fclose(g_caa_mcmc_hsz_eff);
    }



  
  /* Clean up */
  #ifdef DEBUG_PROG
  fprintf(stderr,"Clean up\n");
  #endif
  FREE(wgl_mean_inv_lik);

  err = re_make_graph_gauss(weight->gr_str);
  if(err)
    {
      write_warning("main_model2:Error calling re_make_graph_gauss\n");
      return(err);
    }

  if(i_D_orig->coastal_cod)
    {
      err = re_make_graph_gauss(weight_CC->gr_str);
      if(err)
	{
	  write_warning("main_model2:Error calling re_make_graph_gauss\n");
	  return(err);
	}
    }

  err = re_find_node_effect(D_wgl->glm->xcov,D_wgl->glm->nxcov,weight->gr_str->in_gr,
			    &(weight->gr_str->node));
  if(err)
    {
      write_warning("main_model2:Error calling re_find_node_effect\n");
      return(err);
    }
  if(i_D_orig->coastal_cod)
    {
      err = re_find_node_effect(D_wgl_CC->glm->xcov,D_wgl_CC->glm->nxcov,weight_CC->gr_str->in_gr,
				&(weight_CC->gr_str->node));
      if(err)
	{
	  write_warning("main_model2:Error calling re_find_node_effect\n");
	  return(err);
	}
    }
  
  err = re_initialize_model2(D_wgl,weight,weight_mean,D_wgl_CC,weight_CC,weight_CC_mean,
			     i_D_orig->coastal_cod);
  if(err)
    {
      write_warning("main_model2:Error calling re_initialize_model2\n");
      return(err);
    }

  if(i_D_orig->coastal_cod==0)
    {
      err = re_makedata_wgl_suff(D_wgl);
      if(err)
	{
	  write_warning("main_model2:Error calling re_makedata_wgl_suff\n");
	  return(err);
	}
    }
  else
    {
      err = re_makedata_wgl_suff_CC(D_wgl,D_wgl_CC);
      if(err)
	{
	  write_warning("main_model2:Error calling re_makedata_wgl_suff_CC\n");
	  return(err);
	}
    }

  err = re_makedata_lin1(i_inWgl->num_adj_area,i_inWgl->adj_area,
			 &D_wgl,1);
  if(err)
    {
      write_warning("main_model2:Error calling re_makedata_lin1\n");
      return(err);
    }
  if(i_D_orig->coastal_cod)
    {
      err = re_makedata_lin1(i_inWgl->num_adj_area,i_inWgl->adj_area,
			     &D_wgl_CC,1);
      if(err)
	{
	  write_warning("main_model2:Error calling re_makedata_lin1\n");
	  return(err);
	}
    }

  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"\n CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  #endif


  return(0);
}		/* end of main_model2 */



/*!
  \brief Do haulsize regression 
  \author Hanne Rognebakke

  The routine puts the input data in a format so mail_model2 can be used for doing the haulsize regression.
  It prints the results to file.
*/
int main_model_hsz(Input_common *i_inCommon,Input_wgl *i_inHsz,Input_prior *i_inPrior,Data_orig *i_D_orig)
{
  Input_common *inCommon_hsz;
  Data_orig    *D_orig_hsz;
  char          buffer[MAX_STR];
  int           i,itmp,err=0;

  i_inPrior->wgl_eff_mean = i_inPrior->lga_eff_mean;
  i_inPrior->wgl_eff_prec = i_inPrior->lga_eff_prec;
  i_inPrior->wgl_prec_par = i_inPrior->lga_prec_par;
  i_inPrior->wgl_ar = i_inPrior->lga_ar;

  inCommon_hsz = CALLOC(1,Input_common);
  inCommon_hsz->seed = i_inCommon->seed;
  inCommon_hsz->burn_in = 500;
  inCommon_hsz->num_it_inner = 1;
  inCommon_hsz->num_it_outer = i_inCommon->burn_in+i_inCommon->num_it_inner*i_inCommon->num_it_outer;
  inCommon_hsz->constr = i_inCommon->constr;
  inCommon_hsz->sim_ar = i_inCommon->sim_ar;
  inCommon_hsz->use_debug = i_inCommon->use_debug;
  inCommon_hsz->print_boat = 1;
  inCommon_hsz->num_par2 = CALLOC(2,int);
  itmp = 3; // const. slope + tau_obs + loglik
  for(i=0;i<i_inHsz->int_cov->n_cov;i++)
    {
      itmp += i_inHsz->int_cov->n_lev[i];
      if(i_inHsz->int_cov->random[i])
	itmp += 1;
      if(i_inHsz->int_cov->spatial[i])
	itmp += 1;
    }
  inCommon_hsz->num_par2[0] = itmp;
  inCommon_hsz->num_par2[1] = 0;
  inCommon_hsz->inputfolder = i_inCommon->inputfolder;
  inCommon_hsz->filename_mcmc2 = i_inCommon->filename_hsz_mcmc2;
  inCommon_hsz->filename_hsz_it = i_inCommon->filename_hsz_it;
  inCommon_hsz->filename_hsz_hauleff = i_inCommon->filename_hsz_hauleff;
  inCommon_hsz->inc_hsz = 1;
  inCommon_hsz->print_format = i_inCommon->print_format;


  D_orig_hsz = CALLOC(1,Data_orig);
  D_orig_hsz->nFish = i_D_orig->nHaul;
  D_orig_hsz->nHaul = i_D_orig->nHaul;
  D_orig_hsz->nBoat = i_D_orig->nBoat;
  D_orig_hsz->boat = CALLOC(i_D_orig->nHaul,int);
  D_orig_hsz->nFishBoat = CALLOC(D_orig_hsz->nHaul,int);
  D_orig_hsz->totweight = CALLOC(D_orig_hsz->nHaul,double);
  D_orig_hsz->totlength = CALLOC(D_orig_hsz->nHaul,double);
  D_orig_hsz->replength = CALLOC(D_orig_hsz->nHaul,int);
  for(i=0;i<D_orig_hsz->nHaul;i++)
    {
      D_orig_hsz->boat[i] = i_D_orig->boat[i];
      D_orig_hsz->nFishBoat[i] = 1;
      D_orig_hsz->totlength[i] = 0.0;
      D_orig_hsz->replength[i] = 1;
      D_orig_hsz->totweight[i] = i_D_orig->haulweight[i];
    }
  D_orig_hsz->coastal_cod = 0;

  #ifdef DEBUG_INPUT
  FILE    *fp;
  fprintf(stderr,"Write input haulsize model to file\n");
  sprintf(buffer,"%s/caa_input_hsz.txt",i_inCommon->inputfolder);
  fprintf(stderr,"Write input model1 to file: %s\n",buffer);
  fp = fopen(buffer,"w");   
  fprintf(fp,"nFish %d\n",D_orig_hsz->nFish);
  fprintf(fp,"nBoat %d\n",D_orig_hsz->nBoat);
  fprintf(fp,"nHaul %d\n",D_orig_hsz->nHaul);
  fprintf(fp,"coastal_cod %d\n",D_orig_hsz->coastal_cod);
  fprintf(fp,"haul boat nFishBoat totlength haulweight replength\n");
  for(i=0;i<D_orig_hsz->nHaul;i++)
    {
      fprintf(fp,"%d %d %d %f %f %d \n",i,D_orig_hsz->boat[i],D_orig_hsz->nFishBoat[i],
	      D_orig_hsz->totlength[i],D_orig_hsz->totweight[i],D_orig_hsz->replength[i]);
    }
  fclose(fp);
  #endif

  err = main_model2(inCommon_hsz, i_inHsz, i_inPrior, D_orig_hsz, 1);
  if(err)
    {
      write_warning("main_model_hsz:Error calling main_model2\n");
      return(err);
    }
  fprintf(stderr,"main_model_hsz: main_model2 ended\n\n");

  FREE(D_orig_hsz->nFishBoat);
  FREE(D_orig_hsz->totweight);
  FREE(D_orig_hsz->totlength);
  FREE(D_orig_hsz->replength);
  FREE(D_orig_hsz->boat);
  FREE(D_orig_hsz);
  FREE(inCommon_hsz->num_par2);
  FREE(inCommon_hsz);
  //NB Don't delete inHsz now because it's linked to inLga!!


  return(0);
}		/* end of main_model_hsz */
