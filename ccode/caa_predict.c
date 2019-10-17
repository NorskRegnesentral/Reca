/*!
  \file caa_predict.c
  \brief Containing the main routine for prection catch-at-age
  \author Geir Storvik and Hanne Rognebakke

  The prediction is performed by the ::predict routine.
*/
#include "caa.h"
#include "caa_routines.h"
#include "caa_read_write.h"
#include "caa_init.h"
#include "caa_cell_constr.h"
#include "caa_predict.h"
#include "caa_util.h"
#include "caa_input.h"
#include "caa_sample_g_a.h"

static int      s_N_gauher;  /*!< Number of nodes in Gauss Hermite quadrature */
static double  *s_gauher_x;  /*!< Abscissas in Gauss Hermite quadrature */
static double  *s_gauher_w;  /*!< Weights in Gauss Hermite quadrature */


#ifdef LOG_FILE
extern FILE     *g_caa_log; 
#endif
extern FILE     *g_caa_mcmc1;
extern FILE     *g_caa_mcmc2;
extern FILE     *g_caa_mcmc_hsz;
 

/*!
  \author Geir Storvik and Hanne Rognebakke
  \brief Estimates catch-at-age for simulated parameters and random effects.

  The routine therefore starts to convert input data into 
  approperiate c-structures (as defined in caa.h).

  Thereafter a loop through all simulations is performed with the main
  calculations performed in the ::find_catch_at_age routine.

  The predictions are stored in a int vector, o_mcmc_totcatch with all
  predictions from one iteration in one sequential block. See the
  ::write_it_totcatch routine for the format of this block.

  In addition mean of length-given-age and mean of weight-given-age for each
  simulation is stored in o_mcmc_mean_l and o_mcmc_mean_w. Again the
  ::write_it_totcatch routine described the format.
*/
int predict(Input_predict *i_inPredict, Input_totcatch *i_inCatch, Input_cell *i_inCell, Data_COST *i_D_COST)
{
  /* Data structures */
  Age_struct       *age;             /* Current simulations of age-parameters */ 
  LW_struct        *length;          /* Current simulations of lga-parameters */ 
  LW_struct        *weight;          /* Current simulations of wgl-parameters */ 
  LW_struct        *hsz;          /* Current simulations of wgl-parameters */ 
  TC_struct        *totcatch;        /* Contains all MCMC simulations */
  TC_struct        *totcatch_disc=NULL;
  Data_age         *D_age;           /* Age data */
  Data_lin         *D_lga;           /* Lga data */
  Data_g_a         *D_g_a;           /* g_a parameters and simulations */
  Data_lin         *D_wgl;           /* Wgl data */
  Data_lin         *D_hsz=NULL;      /* Haulsize data */
  Data_totcatch    *D_totcatch;      /**/
  LW_struct        *length_CC=NULL;  /* Current simulations of lga-parameters for coastal cod */ 
  LW_struct        *weight_CC=NULL;  /* Current simulations of wgl-parameters for coastal cod */ 
  Data_lin         *D_lga_CC=NULL;   /* Lga data for coastal cod */
  Data_lin         *D_wgl_CC=NULL;   /* Wgl data for coastal cod */
  Data_g_a         *D_g_a_CC=NULL;   /* g_a parameters and simulations for coastal cod */

  int       err=0;
  long      time_now, time_start;
  char      filename[MAX_STR];
  char      buffer[MAX_STR];

  int       lga_cens_model=0;  // Not implemented in model now. Must write parameters in write_mcmc1 and read parameters in read_mcmc1
  double   *lga_cens_par=NULL;

  time_start = clock();

  #ifdef LOG_FILE
  fprintf(g_caa_log,"Initializing predict\n"); 
  #endif
   /* Calculate Gauss hermite abscissas and weights */
  s_N_gauher = 30;
  s_gauher_x = CALLOC(s_N_gauher+1,double);  //Free ok
  s_gauher_w = CALLOC(s_N_gauher+1,double);  //Free ok
  gauher(s_gauher_x,s_gauher_w,s_N_gauher);

  /* Read input data from binary file */
  /* Open file for reading in binary format */
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Read input predict from binary files %s and %s\n",
	  i_inPredict->filename_mcmc1,i_inPredict->filename_mcmc2);
  #endif
  
  fprintf(stderr,"predict: Open %s for reading age and lga parameters\n",i_inPredict->filename_mcmc1);
  if(!(g_caa_mcmc1 = fopen(i_inPredict->filename_mcmc1, "rb")))
    {
      sprintf(buffer,"predict: Couldn't open file for reading: %s\n",i_inPredict->filename_mcmc1);
      write_warning(buffer);
      return(1);
    }
  err = read_mcmc1(i_inPredict, &D_age, &D_lga, &D_g_a, &D_lga_CC, &D_g_a_CC);
  if(err)
    {
      write_warning("predict:Error calling predict\n");
      return(err);
    } 

  fprintf(stderr,"predict: Open %s for reading wgl parameters\n",i_inPredict->filename_mcmc2);
  if(!(g_caa_mcmc2 = fopen(i_inPredict->filename_mcmc2, "rb")))
    {
      sprintf(buffer,"predict: Couldn't open file for reading: %s\n",i_inPredict->filename_mcmc2);
      write_warning(buffer);
      return(1);
    }
  err = read_mcmc2(i_inPredict, &D_wgl, &D_wgl_CC);
  if(err)
    {
      write_warning("predict:Error calling predict\n");
      return(err);
    } 

  if(i_inPredict->inc_hsz)
    {
      #ifdef DEBUG_PREDICT
      printf("Read haulsize parameters from binary file\n");
      #endif
      fprintf(stderr,"predict: Open %s for reading hsz parameters\n",i_inPredict->filename_hsz_mcmc2);
      if(!(g_caa_mcmc_hsz = fopen(i_inPredict->filename_hsz_mcmc2, "rb")))
	{
	  sprintf(buffer,"predict: Couldn't open file for reading: %s\n",i_inPredict->filename_hsz_mcmc2);
	  write_warning(buffer);
	  return(1);
	}
      err = read_hsz(i_inPredict, &D_hsz);
      if(err)
	{
	  write_warning("predict:Error calling read_hsz\n");
	  return(err);
	}
    }

  #ifdef DEBUG_PREDICT
  printf("Initialize model\n");
  #endif
  err = initialize_predict(D_age,&age,D_lga,&length,D_g_a,
			   D_lga_CC,&length_CC,D_g_a_CC,
			   D_wgl,&weight,D_wgl_CC,&weight_CC,D_hsz,&hsz,
			   lga_cens_model,lga_cens_par,i_inPredict->coastal_cod);
  if(err)
    {
      write_warning("predict:Error calling initialize\n");
      return(err);
    }

  /* Assign mcmc samples to correct structures */
  age->par->num_var = i_inPredict->num_par1[0];
  length->par->num_var = i_inPredict->num_par1[1];
  if(i_inPredict->coastal_cod)
    {
      length_CC->par->num_var = i_inPredict->num_par1[3];
    }
  weight->par->num_var = i_inPredict->num_par2[0];
  if(i_inPredict->coastal_cod)
    {
      weight_CC->par->num_var = i_inPredict->num_par2[1];
    }
  if(i_inPredict->inc_hsz)
    {
      hsz->par->num_var = D_hsz->glm->numpar;
    }

  #ifdef DEBUG_PREDICT
  printf("Make data totcatch\n");
  #endif
  err = makedata_totcatch(D_age,D_lga,D_wgl,D_hsz,
			  i_inCatch,i_inPredict->inc_hsz,&D_totcatch);
  if(err)
    {
      write_warning("predict:Error calling makedata_totcatch\n");
      return(err);
    }
  D_totcatch->nlint = i_inPredict->N_l_int;
  D_totcatch->l_int = i_inPredict->l_int;

  #ifdef DEBUG_PREDICT
  printf("Make data mcmc\n");
  #endif
  err = makedata_mcmc(i_inPredict,D_age->glm->ncat,D_totcatch->nlint,&totcatch);
  if(err)
    {
      write_warning("predict:Error calling makedata_mcmc\n");
      return(err);
    }
 
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log," CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"\nStart run_predict\n");
  #endif
  /* Run predict */
  #ifdef DEBUG_PREDICT
  printf("Run predict\n");
  #endif
  err = run_predict(D_age, age, D_lga, length, D_lga_CC, length_CC, D_g_a, D_g_a_CC,
		    D_wgl, weight, D_wgl_CC, weight_CC, D_hsz, hsz,
		    D_totcatch, totcatch, totcatch_disc, i_inPredict);
  if(err)
    {
      write_warning("predict:Error calling run_predict\n");
      return(err);
    }




  /* Clean up */
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Cleaning up\n");
  #endif
  #ifdef DEBUG_PREDICT
  printf("Cleaning up\n");
  #endif

  fclose(g_caa_mcmc1);
  fclose(g_caa_mcmc2);
  if(i_inPredict->inc_hsz)
    fclose(g_caa_mcmc_hsz);

  FREE(s_gauher_x);
  FREE(s_gauher_w);


  #ifdef DEBUG_PROG
  printf("NB!! Exit predict before cleaning up!!");
  #endif
  if(0){

  err = re_makedata_mcmc(&totcatch);
  if(err)
    {
      write_warning("predict:Error calling re_makedata_mcmc\n");
      return(err);
    }

  err = re_initialize_predict(D_age, age, D_lga, length, D_lga_CC, length_CC, D_g_a,
			      D_wgl, weight, D_wgl_CC, weight_CC, i_inPredict->coastal_cod,
			      D_hsz, hsz);
  if(err)
    {
      write_warning("predict:Error calling re_initialize_predict\n");
      return(err);
    }

  err = re_makedata_totcatch(D_age,D_lga,D_wgl,D_hsz,
			     i_inCatch,i_inPredict->inc_hsz,&D_totcatch);
  if(err)
    {
      write_warning("predict:Error calling re_makedata_totcatch\n");
      return(err);
    }
  err = re_makedata_lin1(NULL,NULL,&D_wgl,0);
  if(err)
    {
      write_warning("predict:Error calling re_makedata_lin1\n");
      return(err);
    }
  if(i_inPredict->coastal_cod)
    {
      err = re_makedata_lin1(NULL,NULL,&D_wgl_CC,0);
      if(err)
	{
	  write_warning("predict:Error calling re_makedata_lin1\n");
	  return(err);
	}
    }
  
  err = re_makedata_lin1(NULL,NULL,&D_lga,0);
  if(err)
    {
      write_warning("predict:Error calling re_makedata_lin1\n");
      return(err);
    }
  if(i_inPredict->coastal_cod)
    {
      err = re_makedata_lin1(NULL,NULL,&D_lga_CC,0);
      if(err)
	{
	  write_warning("predict:Error calling re_makedata_lin1\n");
	  return(err);
	}
    }
  if(i_inPredict->inc_hsz)
    {
      err = re_makedata_lin1(NULL,NULL,&D_hsz,0);
      if(err)
	{
	  write_warning("predict:Error calling re_makedata_lin1\n");
	  return(err);
	}
    }

  err = re_makedata_age1(NULL,NULL,&D_age,0);
  if(err)
    {
      write_warning("predict:Error calling re_makedata_age\n");
      return(err);
    }
  }
  printf("Cleaning up\n");
  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log," CPU time used: %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  #endif


  return(0);
}		/* end of predict */



int run_predict(Data_age *i_D_age, Age_struct *i_age, Data_lin *i_D_lga, LW_struct *i_length, 
		Data_lin *i_D_lga_CC, LW_struct *i_length_CC, Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC,
		Data_lin *i_D_wgl, LW_struct *i_weight, Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, 
		Data_lin *i_D_hsz, LW_struct *i_hsz, 
		Data_totcatch *i_D_totcatch, TC_struct *i_totcatch, TC_struct *i_totcatch_disc,
		Input_predict *i_inPredict)
{
  int      nMCMC,nMCMC_sim,it,burnin_hsz,err=0;
  int      burnin,coastal_cod,class_error,ret;
  double  *mean_l,*mean_w,*dtemp;
  int     *num_par1;
  char     filename[MAX_STR];
  char     buffer[MAX_STR];
  FILE    *caa_pred;

  int      time_now, time_start;

  time_start = clock();

  mean_l = CALLOC(i_D_age->glm->ncat,double); //FREE OK
  mean_w = CALLOC(i_D_age->glm->ncat,double); //FREE OK
  dtemp = CALLOC(2,double); //FREE OK

  nMCMC = i_inPredict->nMCMC;
  burnin = i_inPredict->burnin;
  num_par1 = i_inPredict->num_par1;
  coastal_cod = i_inPredict->coastal_cod;
  if(num_par1[3]>num_par1[1]) //classification error
    class_error = 1;
  else
    class_error = 0;

  /* Open file for printing in binary format */
  sprintf(filename,"%s",i_inPredict->filename_predict);
  fprintf(stderr,"print parameters to file %s\n",filename);
  #ifdef LOG_FILE
  fprintf(g_caa_log,"print parameters to file %s\n",filename);
  #endif
  if(!(caa_pred = fopen(filename, "wb")))
    {
      sprintf(buffer,"run_predict: Couldn't open file for writing: %s\n",filename);
      write_warning(buffer);
      return(1);
    }
  nMCMC_sim = nMCMC-burnin;
  fwrite(&nMCMC_sim,sizeof(int),1,caa_pred);
  fwrite(&i_D_age->glm->ncat,sizeof(int),1,caa_pred);
  fwrite(i_D_age->a_vec,sizeof(int),i_D_age->glm->ncat,caa_pred);
  fwrite(&i_D_totcatch->nlint,sizeof(int),1,caa_pred);
  fwrite(i_D_totcatch->l_int,sizeof(double),i_D_totcatch->nlint,caa_pred);

  #ifdef LOG_FILE
  time_now = clock();
  fprintf(g_caa_log,"CPU time used %fs\n", (double)(time_now-time_start)/CLOCKS_PER_SEC);
  fprintf(g_caa_log,"Perform prediction\n");
  #endif

  #ifdef LOG_FILE
  fprintf(g_caa_log,"Start MCMC iterations\n");
  #endif

  /* Haulsize file includes mcmc samples for use in fitting model1 - including burnin! */
  burnin_hsz = i_inPredict->nMCMC_hsz-i_inPredict->nMCMC;
  //fprintf(stderr,"run_predict:burnin_hsz=%d, nMCMC_hsz=%d\n",burnin_hsz,i_inPredict->nMCMC_hsz);
  for(it=0;it<burnin_hsz;it++)
    {
      if(i_D_age->glm->inc_hsz)
        {
          err = read_it(g_caa_mcmc_hsz,i_D_hsz->glm,i_hsz->par,0,i_inPredict->read_boat);
          if(err)
            {
              write_warning("run_predict:Error calling read_it for haulsize regression\n");
              return(err);
            }
        }
    }

  //printf("run_predict:burnin=%d, nMCMC=%d\n",burnin,nMCMC);

  for(it=0;it<burnin;it++)
    {
      //printf("run_predict: burnin: it=%d\n",it);
      i_D_age->glm->xcov[0]->n_cov--;
      err= read_it(g_caa_mcmc1,i_D_age->glm,i_age->par,0,i_inPredict->read_boat);
      if(err)
        {
          write_warning("run_predict:Error calling read_it for age model\n");
	  return(err);
	}
      i_D_age->glm->xcov[0]->n_cov++;
      err = read_it(g_caa_mcmc1,i_D_lga->glm,i_length->par,0,i_inPredict->read_boat);
      if(err)
        {
          write_warning("run_predict:Error calling read_it for lga model\n");
	  return(err);
	}
      if(i_D_g_a->g_a_npar>0) 
	ret = fread(i_D_g_a->g_a_par,sizeof(double),i_D_g_a->g_a_npar,g_caa_mcmc1);
      if(coastal_cod)
	{
	  err = read_it(g_caa_mcmc1,i_D_lga_CC->glm,i_length_CC->par,class_error,i_inPredict->read_boat);
	  if(err)
	    {
	      write_warning("run_predict:Error calling read_it for lga model\n");
	      return(err);
	    }
	  if(i_D_g_a_CC->g_a_npar>0) 
	    ret = fread(i_D_g_a_CC->g_a_par,sizeof(double),i_D_g_a_CC->g_a_npar,g_caa_mcmc1);
	  if(class_error)
	    {
	      ret = fread(dtemp,sizeof(double),2,g_caa_mcmc1);
	    }
	}
      err = read_it(g_caa_mcmc2,i_D_wgl->glm,i_weight->par,0,i_inPredict->read_boat);
      if(err)
        {
          write_warning("run_predict:Error calling read_it for wgl model\n");
	  return(err);
	}
      if(coastal_cod)
	{
	  err = read_it(g_caa_mcmc2,i_D_wgl_CC->glm,i_weight_CC->par,0,i_inPredict->read_boat);
	  if(err)
	    {
	      write_warning("run_predict:Error calling read_it for wgl model\n");
	      return(err);
	    }
	}
      if(i_inPredict->inc_hsz)
	{
	  err = read_it(g_caa_mcmc_hsz,i_D_hsz->glm,i_hsz->par,0,i_inPredict->read_boat);
	  if(err)
	    {
	      write_warning("run_predict:Error calling read_it for haulsize regression\n");
	      return(err);
	    }
	}
    }

  for(it=burnin;it<nMCMC;it++)
    {
      //printf("\nIt %d:\n",it);
      #ifdef DEBUG_PREDICT
      printf("Read mcmc samples\n");
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"\nIt %d:\nRead mcmc samples\n",it);
      #endif
      i_D_age->glm->xcov[0]->n_cov--;
      //fprintf(stderr,"run_predict:read_it age model\n");
      err= read_it(g_caa_mcmc1,i_D_age->glm,i_age->par,0,i_inPredict->read_boat);
      if(err)
        {
          write_warning("run_predict:Error calling read_it for age model\n");
	  return(err);
	}
      i_D_age->glm->xcov[0]->n_cov++;
      //fprintf(stderr,"run_predict:read_it lga model\n");
      err = read_it(g_caa_mcmc1,i_D_lga->glm,i_length->par,0,i_inPredict->read_boat);
      if(err)
        {
          write_warning("run_predict:Error calling read_it for lga model\n");
	  return(err);
	}
      if(i_D_g_a->g_a_npar>0) 
	{
	  ret = fread(i_D_g_a->g_a_par,sizeof(double),i_D_g_a->g_a_npar,g_caa_mcmc1);
	  err = calc_g_a_S_R2(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,i_D_g_a->g_a);
	}
      if(coastal_cod)
	{
	  err = read_it(g_caa_mcmc1,i_D_lga_CC->glm,i_length_CC->par,class_error,i_inPredict->read_boat);
	  if(err)
	    {
	      write_warning("run_predict:Error calling read_it for lga model\n");
	      return(err);
	    }
	  if(i_D_g_a_CC->g_a_npar>0) 
	    {
	      ret = fread(i_D_g_a_CC->g_a_par,sizeof(double),i_D_g_a_CC->g_a_npar,g_caa_mcmc1);
	      err = calc_g_a_S_R2(i_D_g_a_CC->ncat,i_D_g_a_CC->a_vec,i_D_g_a_CC->g_a_par,i_D_g_a_CC->g_a);
	    }
	  if(class_error)
	    {
	      ret = fread(dtemp,sizeof(double),2,g_caa_mcmc1);
	    }
	}
      //fprintf(stderr,"run_predict:read_it wgl model\n");
      err = read_it(g_caa_mcmc2,i_D_wgl->glm,i_weight->par,0,i_inPredict->read_boat);
      if(err)
        {
          write_warning("run_predict:Error calling read_it for wgl model\n");
	  return(err);
	}
      if(coastal_cod)
	{
	  err = read_it(g_caa_mcmc2,i_D_wgl_CC->glm,i_weight_CC->par,0,i_inPredict->read_boat);
	  if(err)
	    {
	      write_warning("run_predict:Error calling read_it for wgl model\n");
	      return(err);
	    }
	}
      //fprintf(stderr,"run_predict:read_it hsz model\n");
      if(i_inPredict->inc_hsz)
	{
	  err = read_it(g_caa_mcmc_hsz,i_D_hsz->glm,i_hsz->par,0,i_inPredict->read_boat);
	  if(err)
	    {
	      write_warning("run_predict:Error calling read_it for haulsize regression\n");
	      return(err);
	    }
	}
      #ifdef DEBUG_PREDICT
      printf("find_catch_at_age\n");
      #endif
      time_start = clock();
      err = find_catch_at_age(i_D_age,i_age,i_D_lga,i_length,i_D_lga_CC,i_length_CC,
			      i_D_g_a,i_D_g_a_CC,
			      i_D_wgl,i_weight,i_D_wgl_CC,i_weight_CC,
			      i_D_hsz,i_hsz,
			      i_D_totcatch,
			      i_inPredict->n_MC,it,coastal_cod,
			      i_totcatch,mean_l,mean_w);
      if(err)
	{
	  write_warning("run_predict:Error calling find_catch_at_age\n");
	  return(err);
	}
      
      #ifdef LOG_FILE
      time_now = clock();
      fprintf(g_caa_log,"CPU time used in find_catch_at_age in iteration %d: %fs\n", it,
	      (double)(time_now-time_start)/CLOCKS_PER_SEC);
      #endif

      #ifdef DEBUG_PREDICT
      printf("write samples\n");
      #endif
      err = write_it_totcatch(caa_pred, i_D_age->glm->ncat, i_D_totcatch->nlint,
			      i_totcatch, mean_l, mean_w);
      if(err)
        {
          write_warning("run_predict:Error calling write_totcatch\n");
	  return(err);
	}
    }//end for(it=s_burnin;it<nMCMC;it++)

  printf("End run_predict\n");
  /* Close files */    
  fclose(caa_pred);


  FREE(mean_l);
  FREE(mean_w);

  FREE(dtemp);

  return(0);
}           /* end of run_predict */



/*!
  \author Geir Storvik
  \brief Estimates catch-at-age for current parameters and random effects.

  The catches are assumed to be sampled from other hauls than those used in
  the fitting. This means that the haul-effects are unknown. For the lga and
  wgl models this effect can be integrated out. For the age model, however,
  the haul-effects are included through Monte Carlo estimation of the
  age-probablities. Using this both expected age-probablities and expected
  weights can be calculated from which catch-at-age can be predicted.

  For the age model the haul effect is assumed always included and is used in the
  Monte Carlo estimate of the age probabilities.
  For the lga and wgl models, haul effects in the intercept are merged together with
  observation effects in the estimation of the conditional means.
*/
int find_catch_at_age(Data_age *i_D_age, Age_struct *i_age, 
		      Data_lin *i_D_lga, LW_struct *i_length, Data_lin *i_D_lga_CC, LW_struct *i_length_CC,
		      Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC,
		      Data_lin *i_D_wgl, LW_struct *i_weight, Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC,
		      Data_lin *i_D_hsz, LW_struct *i_hsz,
		      Data_totcatch *i_D_totcatch,
		      int *i_nMC, int iter, int i_coastal_cod,
		      TC_struct *i_totcatch, double *i_mean_l, double *i_mean_w)
{
  int     a,c,h,i,inc_hsz,nMC,l,err,ncat_age,season;
  double  A,B,T,sum_T;
  double  mu_lga[2],mu_wgl[2],var_lga,var_wgl;
  double  mu_hsz[0],sd_hsz,hsz,eff_hsz,sum_hsz;
  double  mu_lga_CC[2],mu_wgl_CC[2],var_lga_CC,var_wgl_CC;
  double  eff,sd_age,var,sd_l_fish,sd_l_CC_fish;
  double  pa_sum,sum,mean_w,mean_l,E_a,E_a_tot;
  double  *E_lga,*E_wga, *mu, *pa, *E_pa, *mu_age, **E_pl_a;


  E_lga = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  E_wga = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  mu_age = CALLOC(i_D_age->glm->ncat,double);     // Free ok
  mu = CALLOC(i_D_age->glm->ncat,double);          // Free ok
  pa = CALLOC(i_D_age->glm->ncat,double);          // Free ok
  E_pa = CALLOC(i_D_age->glm->ncat,double);        // Free ok
  E_pl_a = Mmatrix_2d(0,i_D_age->glm->ncat-1,0,i_D_totcatch->nlint-1,sizeof(double),1);        // Free ok


  
  /* Initialize variables */
  sum_T = G_ZERO;
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      i_mean_l[a] = G_ZERO;
      i_mean_w[a] = G_ZERO;
      for(l=0;l<i_D_totcatch->nlint;l++)
	i_totcatch->catch_at_age[a][l] = G_ZERO;
    }
  for(i=0;i<2;i++)
    {
      mu_lga[i] = G_ZERO;
      mu_wgl[i] = G_ZERO;
    }

  inc_hsz = i_D_age->glm->inc_hsz;
  if(inc_hsz)
    {
      sd_hsz = sqrt(G_ONE/i_hsz->par->tau_obs);
    }
    
  /* Age standard deviation */
  if(i_D_age->glm->xcov[0]->iboat>-1) //Boat effect is included here
    {
      sd_age = sqrt(G_ONE/i_age->par->tau_obs + G_ONE/i_age->par->tau[0][i_D_age->glm->xcov[0]->iboat]);
    }
  else					     
    sd_age = sqrt(G_ONE/i_age->par->tau_obs);

  /* Lga variance */
  var_lga = G_ONE/i_length->par->tau_obs;
  if(i_D_lga->glm->xcov[0]->ihaul>0)// Include haul effect in variance and not in effect    
    var_lga += G_ONE/i_length->par->tau[0][i_D_lga->glm->xcov[0]->ihaul];
  if(i_D_lga->glm->xcov[0]->iboat>0)// Include boat effect in variance and not in effect
    var_lga += G_ONE/i_length->par->tau[0][i_D_lga->glm->xcov[0]->iboat];
  sd_l_fish = G_ONE/i_length->par->tau_obs;

  if(i_coastal_cod)
    {
      var_lga_CC = G_ONE/i_length_CC->par->tau_obs;
      if(i_D_lga_CC->glm->xcov[0]->ihaul>0)
	var_lga_CC += G_ONE/i_length_CC->par->tau[0][i_D_lga_CC->glm->xcov[0]->ihaul];
      if(i_D_lga_CC->glm->xcov[0]->iboat>0)
	var_lga_CC += G_ONE/i_length_CC->par->tau[0][i_D_lga_CC->glm->xcov[0]->iboat];
      sd_l_CC_fish = G_ONE/i_length_CC->par->tau_obs;
    }

  /* Wgl variance */
  var_wgl = G_ONE/i_weight->par->tau_obs;
  if(i_D_wgl->glm->xcov[0]->ihaul>0)// Include haul effect in variance and not in effect
    var_wgl += G_ONE/i_weight->par->tau[0][i_D_wgl->glm->xcov[0]->ihaul];
  if(i_D_wgl->glm->xcov[0]->iboat>0)// Include boat effect in variance and not in effect
    var_wgl += G_ONE/i_weight->par->tau[0][i_D_wgl->glm->xcov[0]->iboat];
      
  if(i_coastal_cod)
    {
      var_wgl_CC = G_ONE/i_weight_CC->par->tau_obs;
      if(i_D_wgl_CC->glm->xcov[0]->ihaul>0)
	var_wgl_CC += G_ONE/i_weight_CC->par->tau[0][i_D_wgl_CC->glm->xcov[0]->ihaul];
      if(i_D_wgl_CC->glm->xcov[0]->iboat>0)
	var_wgl_CC += G_ONE/i_weight_CC->par->tau[0][i_D_wgl_CC->glm->xcov[0]->iboat];
    }

  if(i_coastal_cod)
    ncat_age = (int) i_D_age->glm->ncat/2;
  else
    ncat_age = i_D_age->glm->ncat;

 
  /* Start prediction over all cells */
  #ifdef DEBUG_PREDICT
  printf("Start prediction over all cells\n");
  #endif
  for(c=0;c<i_D_totcatch->nCell;c++)
    {
      season = i_D_totcatch->season[c];
      #ifdef DEBUG_PREDICT
      printf("Calculate means for age\n");
      #endif
      err = calculate_means_age(mu_age,c,i_age,i_D_age,i_D_totcatch);
      if(err)
	{
	  write_warning("find_catch_at_age:Error calling calculate_means_age\n");
	  return(err);
	}

      if(inc_hsz)
	{
          #ifdef DEBUG_PREDICT
	  printf("Calculate means for haulsize\n");
          #endif
	  err = calculate_means_hsz(mu_hsz,c,i_hsz,i_D_hsz,i_D_totcatch);
	  if(err)
	    {
	      write_warning("find_catch_at_age:Error calling calculate_means_hsz\n");
	      return(err);
	    }
	}

      #ifdef DEBUG_PREDICT
      printf("Calculate means for lga\n");
      #endif
      /* Calculate means for lga */
      err = calculate_means_lga(mu_lga,c,i_length,i_D_lga,i_D_totcatch);
      if(err)
	{
	  write_warning("find_catch_at_age:Error calling calculate_means_lga\n");
	  return(err);
	}
      if(i_coastal_cod)
	{
          #ifdef DEBUG_PREDICT
	  printf("Calculate means for lga - coastal cod\n");
          #endif
	  err = calculate_means_lga(mu_lga_CC,c,i_length_CC,i_D_lga_CC,i_D_totcatch);
	  if(err)
	    {
	      write_warning("find_catch_at_age:Error calling calculate_means_lga\n");
	      return(err);
	    }
	}

      #ifdef DEBUG_PREDICT
      printf("Calculate means for wgl\n");
      #endif
      /* Calculate means for wgl */
      err = calculate_means_wgl(mu_wgl,c,i_weight,i_D_wgl,i_D_totcatch);
      if(err)
	{
	  write_warning("find_catch_at_age:Error calling calculate_means_wgl\n");
	  return(err);
	}
      if(i_coastal_cod)
	{
          #ifdef DEBUG_PREDICT
	  printf("Calculate means for wgl- coastal cod\n");
          #endif
	  err = calculate_means_wgl(mu_wgl_CC,c,i_weight_CC,i_D_wgl_CC,i_D_totcatch);
	  if(err)
	    {
	      write_warning("find_catch_at_age:Error calling calculate_means_wgl\n");
	      return(err);
	    }
	}

      #ifdef DEBUG_PREDICT
      printf("Monte Carlo estimation\n");
      #endif
      /* Monte Carlo estimation  */
      for(a=0;a<i_D_age->glm->ncat;a++)
	E_pa[a] = G_ZERO;      
      nMC = i_nMC[c];
      sum_hsz = G_ZERO;
      for(h=0; h<nMC; h++)
	{
	  if(inc_hsz)
	    {
	      eff_hsz = gennor(G_ZERO,sd_hsz);
	      //hsz = calc_eff(i_D_hsz->glm->xcov[0],i_hsz->par->eff[0][0],0) + eff_hsz;
	      hsz = mu_hsz[0] + eff_hsz;
	      //if(c==1 && h<3)
	      //	printf("mu_hsz[0]=%f,hsz_eff=%f,hsz=%f\n",mu_hsz[0],calc_eff(i_D_hsz->glm->xcov[0],i_hsz->par->eff[0][0],0),hsz);
	      if(hsz<0.0001)
		hsz = 0.0001;
	      hsz = exp(hsz);
	    }
	  else
	    hsz = 1;
	  sum_hsz += hsz;
          pa_sum = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      if(inc_hsz)
		{
		  eff = mu_age[a] + eff_hsz*i_age->par->eff[a][0][i_D_age->glm->xcov[0]->ihaulsize][0] + gennor(G_ZERO,sd_age);
		}
	      else
		{
		  eff = mu_age[a] + gennor(G_ZERO,sd_age);
		}
	      pa[a] = exp(eff);
	      pa_sum += pa[a];
	    }
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    pa[a] /= pa_sum;
	  if(i_age->delta_age>0) //subtract delta_age from age probabilities
	    {
	      pa_sum = G_ZERO;
	      for(a=0;a<i_D_age->glm->ncat;a++)
		{
		  pa[a]=pa[a]-i_age->delta_age;
		  if(pa[a]<0)
		    pa[a]=0;
		  pa_sum += pa[a];
		} 
	      for(a=0;a<i_D_age->glm->ncat;a++)
		pa[a] /= pa_sum;
	    }
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    E_pa[a] += pa[a]*hsz;
	}
      for(a=0;a<i_D_age->glm->ncat;a++)
	E_pa[a] /= sum_hsz;


      #ifdef DEBUG_PREDICT
      Data_cov *xcov;
      printf("mu_lga=(%lf,%lf),mu_wgl=(%lf,%lf)\n",
	     mu_lga[0],mu_lga[1],mu_wgl[0],mu_wgl[1]);
      xcov = i_D_lga->glm->xcov[0];
      printf("prec_lga_obs=%lf,prec_lga_haul=%lf,var_lga=%lf\n",
	     i_length->par->tau_obs,i_length->par->tau[0][xcov->ihaul],var_lga);
      xcov = i_D_wgl->glm->xcov[0];
      printf("prec_wgl_obs=%lf,prec_wgl_haul=%lf,var_wgl=%lf\n",
	     i_weight->par->tau_obs,i_weight->par->tau[0][xcov->ihaul],var_wgl);
      if(i_coastal_cod)
	{
	  printf("mu_lga_CC=(%lf,%lf),mu_wgl_CC=(%lf,%lf)\n",
		 mu_lga_CC[0],mu_lga_CC[1],mu_wgl_CC[0],mu_wgl_CC[1]);
	  xcov = i_D_lga_CC->glm->xcov[0];
	  printf("prec_lga_obs=%lf,prec_lga_haul=%lf,var_lga=%lf\n",
		 i_length_CC->par->tau_obs,i_length_CC->par->tau[0][xcov->ihaul],var_lga_CC);
	  xcov = i_D_wgl_CC->glm->xcov[0];
	  printf("prec_wgl_obs=%lf,prec_wgl_haul=%lf,var_wgl=%lf\n",
		 i_weight_CC->par->tau_obs,i_weight_CC->par->tau[0][xcov->ihaul],var_wgl_CC);
	}
      printf("mu_age E[p(a)] E_lga E_wga\n");
      #endif

      mean_w = G_ZERO;
      mean_l = G_ZERO;   
      A = mu_wgl[0] + mu_wgl[1]*mu_lga[0];
      B = mu_wgl[1]*mu_lga[1];
      var = mu_wgl[1]*mu_wgl[1]*var_lga + var_wgl;
      for(a=0;a<ncat_age;a++)
	{
	  E_lga[a] = exp(mu_lga[0]+mu_lga[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]+season-1]+G_HALF*var_lga);
	  E_wga[a] = exp(A+B*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]+season-1]+G_HALF*var);
	  mean_l += E_pa[a] * E_lga[a];
	  mean_w += E_pa[a] * E_wga[a];
          #ifdef DEBUG_PREDICT
	  if(a<5)
	  printf("%d %lf %lf %lf %lf\n",a,mu_age[a],E_pa[a],E_lga[a],E_wga[a]);
          #endif
	  mu[a] = mu_lga[0]+mu_lga[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]+season-1];
	  if(i_D_totcatch->nlint>1)
	    {
	      E_pl_a[a][0] = pnorm(i_D_totcatch->l_int[0],mu[a],sd_l_fish);
	      sum = E_pl_a[a][0];
	      for(l=1;l<(i_D_totcatch->nlint-1);l++)
		{
		  E_pl_a[a][l] = pnorm(i_D_totcatch->l_int[l],mu[a],sd_l_fish);
		  E_pl_a[a][l] -= sum;
		  sum += E_pl_a[a][l];
		}
	      E_pl_a[a][i_D_totcatch->nlint-1] = G_ONE-sum;
	    }
	  else
	    E_pl_a[a][0] = G_ONE;
	}
      if(i_coastal_cod)
	{
	  A = mu_wgl_CC[0] + mu_wgl_CC[1]*mu_lga_CC[0];
	  B = mu_wgl_CC[1]*mu_lga_CC[1];
	  var = mu_wgl_CC[1]*mu_wgl_CC[1]*var_lga_CC + var_wgl_CC;
	  for(a=ncat_age;a<i_D_age->glm->ncat;a++)
	    {
	      E_lga[a] = exp(mu_lga_CC[0]+mu_lga_CC[1]*i_D_g_a_CC->g_a[i_D_g_a_CC->a2Age_vec[a]+season-1]+G_HALF*var_lga_CC);
	      E_wga[a] = exp(A+B*i_D_g_a_CC->g_a[i_D_g_a_CC->a2Age_vec[a]+season-1]+G_HALF*var);
	      mean_l += E_pa[a] * E_lga[a];
	      mean_w += E_pa[a] * E_wga[a];
              #ifdef DEBUG_PREDICT
	      printf("%lf %lf %lf %lf\n",mu_age[a],E_pa[a],E_lga[a],E_wga[a]);
              #endif
	      mu[a] = mu_lga_CC[0]+mu_lga_CC[1]*i_D_g_a_CC->g_a[i_D_g_a_CC->a2Age_vec[a]+season-1];
	      if(i_D_totcatch->nlint>1)
		{
		  E_pl_a[a][0] = pnorm(i_D_totcatch->l_int[0],mu[a],sd_l_CC_fish);
		  sum = E_pl_a[a][0];
		  for(l=1;l<(i_D_totcatch->nlint-1);l++)
		    {
		      E_pl_a[a][l] = pnorm(i_D_totcatch->l_int[l],mu[a],sd_l_CC_fish);
		      E_pl_a[a][l] -= sum;
		      sum += E_pl_a[a][l];
		    }
		  E_pl_a[a][i_D_totcatch->nlint-1] = G_ONE-sum;
		}
	      else
		E_pl_a[a][0] = G_ONE;
	    }
	}//end if(i_coastal_cod)

      T = i_D_totcatch->catch[c]/mean_w;
      if(!(T> -9999.0  && T < 999999999999999999.99))
	{
	  write_warning("find_catch_at_age:Something is wrong 2\n");
	  fprintf(stderr,"age mu_age E_p(a) E_lga E_wga\n");
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      fprintf(stderr,"%d %lf %lf %lf %lf\n",a,mu_age[a],E_pa[a],E_lga[a],E_wga[a]);
	    }
	  fprintf(stderr,"mean_w=%lf,sd_age=%lf\n",mean_w,sd_age);
	  fprintf(stderr,"mu_lga = ");
	  for(i=0;i<i_D_lga->glm->nxcov;i++)
	    fprintf(stderr,"%lf ",mu_lga[i]);
	  fprintf(stderr,"mu_wgl = ");
	  for(i=0;i<i_D_wgl->glm->nxcov;i++)
	    fprintf(stderr,"%lf ",mu_wgl[i]);
	  fprintf(stderr,"A=%lf, B=%lf,var=%lf\n",A,B,var);
	  return(1);
	}
      #ifdef DEBUG_PREDICT
      printf("c=%d,catch=%lf,A=%lf,B=%lf,var=%lf,mean_l=%lf,mean_w=%lf,T=%lf\n",
	     c,i_D_totcatch->catch[c],A,B,var,mean_l,mean_w,T);
      #endif
 
      sum_T += T;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  E_a = 0.0;
	  sum = G_ZERO;
	  for(l=0;l<i_D_totcatch->nlint;l++)
	    {
	      i_totcatch->catch_at_age[a][l] += T * E_pa[a]*E_pl_a[a][l];
	      E_a += T * E_pa[a]*E_pl_a[a][l];
	      sum += i_totcatch->catch_at_age[a][l];
	      //      #ifdef DEBUG_PREDICT
	      //if(iter<3)
	      //printf("T=%lf,E_pa=%lf,E_pl_a=%lf,catch_at_age=%lf\n",
	      //	     T,E_pa[a],E_pl_a[a][l],i_totcatch->catch_at_age[a][l]);
              //#endif
	    }
	  i_mean_l[a] += E_lga[a]*E_a;
	  i_mean_w[a] += E_wga[a]*E_a;
	}
      
    }// end for(c=0;c<i_D_totcatch->nCell;c++)
  
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      E_a_tot = 0;
      for(l=0;l<i_D_totcatch->nlint;l++)
	E_a_tot += i_totcatch->catch_at_age[a][l];

      if(E_a_tot<1E-12)
	{
	  E_a_tot = G_ZERO;
	}
      else
	{
	  i_mean_l[a] /= E_a_tot;
	  i_mean_w[a] /= E_a_tot;
	}
    }


  #ifdef DEBUG_PREDICT
  sum = G_ZERO;
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      printf("a=%d,p_a=%lf\n",a,E_pa[a]);
      for(l=0;l<i_D_totcatch->nlint;l++)
	{
	  sum += i_totcatch->catch_at_age[a][l];
	  printf("%lf ",i_totcatch->catch_at_age[a][l]);
	}
      printf("\n");
    }
  #endif

  FREE(E_lga);
  FREE(E_wga);
  FREE(mu_age);
  FREE(mu);
  FREE(pa);
  FREE(E_pa);
  Fmatrix_2d(&E_pl_a[0][0],&E_pl_a[0]);

  return(0);
}           /* end of find_catch_at_age */


/*!
  \author Geir Storvik
  \brief Calculates intercept and slope for age model
*/
int calculate_means_age(double *mu_age,int c,Age_struct *i_age,Data_age *i_D_age,
			Data_totcatch *i_D_totcatch)
{
  int a,j,j2,k,ncov;
  double eff;
  Data_cov *xcov;
  
  xcov = i_D_age->glm->xcov[0];
  ncov = xcov->n_cov;
  ncov -= 1;  // Excluding haul effects in intercept since this is treated later
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      mu_age[a] = G_ZERO;
      for(j=0;j<ncov;j++)
	{
	  k = i_D_totcatch->age_xcov[0]->c_cov[c][j];
	  if(j!=xcov->ihaulsize && k>-1)
	    {
	      eff = i_age->par->eff[a][0][j][k];
	    }
	  else if(j==xcov->ihaulsize || j==xcov->iboat)
	    {
	      eff = 0; // Excluding boat and haulsize effects in intercept since this is treated later
	    }
	  else if(j==xcov->ispat)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,
		      "calculate_means_age:Missing area c=%d not implemented yet\n",c);
	      fprintf(g_caa_log,"Age model, a=%d,j=%d,k=%d\n",a,j,k);
	      fprintf(g_caa_log,"c=%d",c);
	      for(k=0;k<xcov->n_cov;k++)
		fprintf(g_caa_log,"%d ",i_D_totcatch->age_xcov[0]->c_cov[c][k]);
	      fprintf(g_caa_log,"\n");
              #endif
	      write_warning("calculate_means_age:Something is wrong\n");
	      return(1);
	    }
	  else if(i_D_age->glm->xcov[0]->fix[j]==0 && j!=xcov->ispat && j!=xcov->iboat)
	    {
	      printf("boat? j=%d: eff[%d][%d][%d][%d]=%f\n",j,a,0,j,k,i_age->par->eff[a][0][j][k]);
	      eff = gennor(G_ZERO,G_ONE/sqrt(i_age->par->tau[0][j]));   
	    }
	  else
	    {
              #ifdef LOG_FILE
	      for(j2=0;j2<xcov->n_cov;j2++)
		fprintf(g_caa_log,"%d ",i_D_totcatch->age_xcov[0]->c_cov[c][j2]);
	      fprintf(g_caa_log,"\n");
	      fprintf(g_caa_log,"Missing fixed effect in age model (c=%d,j=%d) not allowed\n",c,j);
              #endif
	      write_warning("calculate_means_age:Something is wrong\n");
	      return(1);
	    }
	  mu_age[a] += eff;
	}
    }
     
  return(0);
}            /* End of calculate_means_age */



/*!
  \author Geir Storvik
  \brief Calculates intercept and slope for lga model
*/
int calculate_means_lga(double *mu_lga,int c,LW_struct *i_length,Data_lin *i_D_lga,
			Data_totcatch *i_D_totcatch)
{
  int       i,j,k,ncov;
  double    eff;
  Data_cov *xcov;
  
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      xcov = i_D_lga->glm->xcov[i];
      mu_lga[i] = G_ZERO;
      ncov = xcov->n_cov;
      if(xcov->ihaul>0)// Include haul effect in variance and not in effect
	{
	  ncov -= 1;   
	}
      if(xcov->iboat>0)// Include haul effect in variance and not in effect
	{
	  ncov -= 1;   
	}
      for(j=0;j<ncov;j++)
	{
	  k=i_D_totcatch->lga_xcov[i]->c_cov[c][j];
	  if(k > -1)
	    eff = i_length->par->eff[0][i][j][k];
	  else if(i_D_lga->glm->xcov[i]->fix[j]==0&&j != xcov->ispat) 
	    {
	      eff = gennor(G_ZERO,sqrt(G_ONE/i_length->par->tau[i][j]));
	    }
	  else if(j==xcov->ispat)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"calculate_means_lga:Missing area c=%d not implemented yet\n",c);
	      for(k=0;k<xcov->n_cov;k++)
		fprintf(g_caa_log,"%d ",i_D_totcatch->lga_xcov[i]->c_cov[c][k]);
	      fprintf(g_caa_log,"\n");
              #endif
	      write_warning("calculate_means_lga:Missing area effect in lga model is not allowed\n");
	      return(1);
	    }
	  else
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"Missing fixed effect in lga model (c=%d,j=%d) not allowed\n",c,j);
              #endif
	      write_warning("calculate_means_lga:Something is wrong\n");
	      return(1);
	    }
	  mu_lga[i] += eff;
	}
    }

  return(0);
}            /* End of calculate_means_lga */

/*!
  \author Geir Storvik
  \brief Calculates intercept and slope for wgl model
*/
int calculate_means_wgl(double *mu_wgl,int c,LW_struct *i_weight,Data_lin *i_D_wgl,
			Data_totcatch *i_D_totcatch)
{
  int       i,j,k,ncov;
  double    eff;
  Data_cov *xcov;

  for(i=0;i<i_D_wgl->glm->nxcov;i++)
    {
      xcov = i_D_wgl->glm->xcov[i];
      ncov = xcov->n_cov;
      mu_wgl[i] = G_ZERO;
      if(xcov->ihaul>0)// Include haul effect in variance and not in effect
	{
	  ncov -= 1;   
	}
      if(xcov->iboat>0)// Include boat effect in variance and not in effect
	{
	  ncov -= 1;   
	}
      for(j=0;j<ncov;j++)
	{
	  k=i_D_totcatch->wgl_xcov[i]->c_cov[c][j];
	  if(k > -1)
	    eff = i_weight->par->eff[0][i][j][k];
	  else if(i_D_wgl->glm->xcov[i]->fix[j]==0&&j != xcov->ispat) 
	    eff = gennor(G_ZERO,sqrt(G_ONE/i_weight->par->tau[i][j]));
	  else if(j==xcov->ispat)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"calculate_means_wgl:Missing area c=%d not implemented yet\n",c);
	      fprintf(g_caa_log,"c=%d",c);
	      for(k=0;k<xcov->n_cov;k++)
		fprintf(g_caa_log,"%d ",i_D_totcatch->wgl_xcov[i]->c_cov[c][k]);
	      fprintf(g_caa_log,"\n");
              #endif
	      write_warning("calculate_means_wgl:Something is wrong\n");
	      return(1);
	    }
	  else
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"Missing fixed effect in wgl model (c=%d,j=%d) not allowed\n",c,j);
              #endif
	      write_warning("calculate_means_wgl:Something is wrong\n");
	      return(1);
	    }
	  mu_wgl[i] += eff;
	}
    }

  return(0);
}            /* End of calculate_means_wgl */

/*!
  \author Hanne Rognebakke
  \brief Calculates intercept and slope for haulsize model
*/
int calculate_means_hsz(double *mu_hsz,int c,LW_struct *i_hsz,Data_lin *i_D_hsz,
			Data_totcatch *i_D_totcatch)
{
  int       j,k,ncov;
  double    eff;
  Data_cov *xcov;
  
  xcov = i_D_hsz->glm->xcov[0];
  mu_hsz[0] = G_ZERO;
  ncov = xcov->n_cov;
  if(xcov->ihaul>0)// Include haul effect in variance and not in effect
    {
      ncov -= 1;   
    }
  for(j=0;j<ncov;j++)
    {
      k=i_D_totcatch->hsz_xcov[0]->c_cov[c][j];
      //if(j!= xcov->icell && k > -1)
      if(k > -1)
	eff = i_hsz->par->eff[0][0][j][k];
      //else if (j==xcov->icell)
      //	eff = i_hsz->par->cell[0][k];
      else if(i_D_hsz->glm->xcov[0]->fix[j]==0&&j != xcov->ispat) 
	{
	  eff = gennor(G_ZERO,sqrt(G_ONE/i_hsz->par->tau[0][j]));
	}
      else if(j==xcov->ispat)
	{
	  write_warning("calculate_means_hsz:Missing area effect in hsz model is not allowed\n");
	  return(1);
	}
      else
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"Missing fixed effect in hsz model (c=%d,j=%d) not allowed\n",c,j);
          #endif
	  write_warning("calculate_means_hsz:Something is wrong\n");
	  return(1);
	}
      mu_hsz[0] += eff;
    }

  return(0);
}            /* End of calculate_means_hsz */

