/*!
  \file caa_main.c
  \brief Containing the main routines for MCMC simulation
  \author Geir Storvik, Hanne Rognebakke

*/

/*Include Files:*/
#include "caa.h"
#include "caa_mcmc.h"
#include "caa_read_write.h"
#include "caa_sample_multi.h"
#include "caa_sample_gauss.h"
#include "caa_sample_g_a.h"
#include "caa_routines.h"
#include "caa_evaluate.h"
//#include "caa_COST.h"

#ifdef LOG_FILE
extern FILE     *g_caa_log;
#endif

extern FILE     *g_caa_mcmc_hsz; 
  
  

/*!
  \author Geir Storvik
  \brief  Initialization for MCMC simulations of age and lga model

  Start with initialization of parameters. This is mainly based on using data
  which is complete (that is amigo-type data and age-stratified by length data)
  and running a few simulations on this.

  Note 
*/
int MCMC_model1_init(Data_orig *i_D_orig, Data_age *i_D_age, Age_struct *i_age,
		     Data_lin *i_D_lga, LW_struct *i_length, 
		     Data_lin *i_D_lga_CC, LW_struct *i_length_CC,
		     Data_g_a *i_D_g_a,Data_g_a *i_D_g_a_CC,Data_CC *i_D_CC,
		     Data_COST *i_D_COST,int i_use_debug)
{
  int        err=0,save_sim=0;

  if(i_use_debug)
    save_sim = 1;
  #ifdef DEBUG_PROG
  printf("MCMC_model1_init:Sample ages init\n");
  save_sim = 1;
  #endif


  err = sample_ages_init(i_D_orig,i_D_CC,i_D_age,i_D_lga,i_D_g_a,i_D_lga_CC,i_D_g_a_CC,save_sim);
  if(err)
    {
      write_warning("MCMC_model1_init:Error calling sample_ages_init\n");
      return(err);
    }


  // Initialization alpha-values for full data
  //Using a small precision when optimization gives alpha-values close
  //to the observations and seems to result in realisations with high
  //posterior densities. It seems however that convergence of the haul-effects
  //goes very slowly and that it is better to not do this. This choice
  //is therefore commented out.
  //i_age->par->tau_obs = 0.0001;
  #ifdef DEBUG_PROG
  printf("MCMC_model1: call age_haul_modes\n");
  #endif
  err = age_haul_modes(0,i_D_age->glm->nHaul,i_age,i_D_age);
  //i_age->par->tau_obs = 1;

  /* Initial fit */ 
  #ifdef DEBUG_PROG
  int a,h,sum,suma;
  double *p,psum;
  FILE *fp;
  fp = fopen("data.txt","w");
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      for(a=0;a<i_D_age->glm->ncat;a++)
	fprintf(fp,"%d ",(int) i_D_age->Ages[h][a]);
      fprintf(fp,"\n");
    }
  fclose(fp);
  p = CALLOC(i_D_age->glm->ncat,double);
  sum = 0;
  printf("MCMC_model1_init: Num_age= ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      suma = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	suma += (int) i_D_age->Ages[h][a];
      printf("%d ",suma);
      sum += suma;
    }
  printf("\n");
  printf("Total aged fish=%d\n",sum);
  psum = G_ZERO;
  printf("Mean age alpha over all hauls, a: 1/h*sum_h exp(i_age->alpha[h][a]):\n");

  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      p[a] = G_ZERO;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	{
	  p[a] += exp(i_age->alpha[h][a]);
	}
      p[a] /= i_D_age->glm->nHaul;
      psum += p[a];
    }
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      printf("const=%lf ; ",i_age->par->eff[a][0][0][0]);
      printf("p(%d)=%lf\n",a,p[a]/psum);
    }
  printf("sigma_haul=%lf\n",sqrt(G_ONE/i_age->par->tau_obs));
  printf("Initial fit of lga model:\n");
  printf("fixed model=%d\n",i_length->fixed_model);
  FREE(p);
  #endif


  if(i_length->fixed_model == 0)
    {
      err = MCMC_it_lga(i_D_lga,i_D_g_a,i_length,0,0,i_D_lga->glm->nHaul,i_age->age_errors,i_D_orig->coastal_cod);
      if(err)
	{
	  write_warning("MCMC_model1_init:Error calling MCMC_it_lga\n");
	  return(err);
	}
      if(i_D_orig->coastal_cod)
	{
	  err = MCMC_it_lga(i_D_lga_CC,i_D_g_a_CC,i_length_CC,0,0,i_D_lga->glm->nHaul,i_age->age_errors,i_D_orig->coastal_cod);
	  if(err)
	    {
	      write_warning("MCMC_model1_init:Error calling MCMC_it_lga\n");
	      return(err);
	    }
	}
    }
  else /* use fixed lga model */
    {
      err = MCMC_it_lga_fixed(i_D_lga,i_D_g_a,i_length,0,0);
      if(err)
	{
	  write_warning("MCMC_model1_init:Error calling MCMC_it_lga_fixed\n");
	  return(err);
	}
      if(i_D_orig->coastal_cod)
	{
	  err = MCMC_it_lga_fixed(i_D_lga_CC,i_D_g_a_CC,i_length_CC,0,0);
	  if(err)
	    {
	      write_warning("MCMC_model1_init:Error calling MCMC_it_lga_fixed\n");
	      return(err);
	    }
	}
    }

  return(0);
}



/*!
  \author Geir Storvik
  \brief  Perform MCMC simulations of age and lga model

  The simulations are divided into burn-in and simulations after burnin in which
  num_it_inner simulations are performed for each num_it_outer simulation.
  All inner simulations are performed through the ::MCMC_model1_it routine.

  Simulations are saved for each num_it_outer iteration using the ::write_it
  routine.
*/
int MCMC_model1(Input_common *i_inCommon,Data_orig *i_D_orig, 
		Data_age *i_D_age, Age_struct *i_age, Age_struct *i_age_mean,
		Data_lin *i_D_lga, LW_struct *i_length, LW_struct *i_length_mean,
		Data_lin *i_D_lga_CC, LW_struct *i_length_CC, LW_struct *i_length_CC_mean,
		Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, Data_CC *i_D_CC,
		double *i_mod1_mean_inv_lik, double *i_age_mean_inv_lik, double *i_lga_mean_inv_lik,
		Data_COST *i_D_COST)
{
  int        err=0,force_acc,ret;
  int        it,iti,ito,it_tot,a,h;
  int       *acc_h=NULL;
  int        save_sim=0,write_alpha=0;
  FILE      *fp;
  char       filename[MAX_STR];
  char       buffer[MAX_STR];

  #ifdef LOG_FILE
  fprintf(g_caa_log,"Start MCMC\n");
  #endif

  #ifdef DEBUG_PROG
  int printAgeDist=1;
  if(printAgeDist)
    {
      fp = fopen("ageDist.dat","w");
      printf("MCMC_model1: print simulated age probabilities to file ageDist.dat\n");
      fprintf(fp,"iteration haul");
      for(a=0;a<i_D_age->glm->ncat;a++)
	fprintf(fp," age%d",i_D_age->a_vec[a]);
      fprintf(fp,"\n");
    }
  #endif
  #ifdef WRITE_ALPHAS
  FILE      *fp_alpha;
  fp_alpha = fopen("alphas.bin","wb");
  fclose(fp_alpha);
  #endif

  force_acc = 0;
  it_tot = 0;
  acc_h = CALLOC(i_D_age->glm->nHaul,int);     // Free ok

  if(i_inCommon->inc_hsz)
    {
      sprintf(filename,"%s",i_inCommon->filename_hsz_it);
      if(!(g_caa_mcmc_hsz = fopen(filename, "rb")))
	{
	  sprintf(buffer,"MCMC_model1: Couldn't open file for reading: %s\n",filename);
	  write_warning(buffer);
	  return(1);
	}
    }
  
  #ifdef LOG_FILE
  if(i_inCommon->burn_in>0)
    fprintf(g_caa_log,"Start MCMC burnin iterations\n");
  #endif
  if(i_inCommon->burn_in>0)
    fprintf(stderr,"Start MCMC burnin iterations\n");
  for(iti=0;iti<i_inCommon->burn_in;iti++)
    {
      if(iti%100 == 0) // print each 100 iteration 
	fprintf(stderr,"\nIteration %d\n",iti);
      if(i_inCommon->inc_hsz)
	{
	  ret = fread(i_age->par->eff_hsz,sizeof(double),i_D_age->glm->nHaul,g_caa_mcmc_hsz);
	  ret = fread(&i_age->par->tau_hsz,sizeof(double),1,g_caa_mcmc_hsz);
	}
      err = MCMC_model1_it(0,force_acc,1,iti,acc_h,it_tot,
			   i_D_orig,i_D_age,i_age,i_D_lga,i_length,i_D_lga_CC,i_length_CC,
			   i_D_g_a,i_D_g_a_CC,i_D_CC,i_D_COST,save_sim,write_alpha);
      it_tot++;
      if(err)
	{
	  write_warning("MCMC_model1:Error calling MCMC_model1_it\n");
	  return(err);
	}
      force_acc = 0;
    }

  /* Start full simulation */
  #ifdef WRITE_ALPHAS
  write_alpha = 1;
  #endif
  it = 0;
  force_acc = 0;
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Start MCMC iterations\n");
  #endif
  fprintf(stderr,"Start MCMC iterations\n");
  for(ito=0;ito<i_inCommon->num_it_outer;ito++)
    {
      //printf("\n\nouter it=%d\n",ito);
      for(iti=0;iti<i_inCommon->num_it_inner;iti++)
	{
          #ifdef DEBUG_PROG
	  printf("\n\nIteration %d\n",it);
          #endif
	  if(it%100 == 0) // print each 100 iteration 
	    fprintf(stderr,"\nIteration %d\n",it);
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\n\nIteration %d\n",it);
          #endif
	  if(i_inCommon->inc_hsz)
	    {
	      ret = fread(i_age->par->eff_hsz,sizeof(double),i_D_age->glm->nHaul,g_caa_mcmc_hsz);
	      ret = fread(&i_age->par->tau_hsz,sizeof(double),1,g_caa_mcmc_hsz);
	      //printf("caa_mcmc.c: \n\nIteration %d\n",it);
	      //printf("MCMC_model1: eff_hsz=%f %f %f, tau=%f\n",i_age->par->eff_hsz[0],i_age->par->eff_hsz[1],i_age->par->eff_hsz[2],i_age->par->tau_hsz);
	    }
	  if((it_tot==(i_inCommon->burn_in+i_inCommon->num_it_inner*i_inCommon->num_it_outer-1))&&i_inCommon->use_debug)
	    save_sim=1;
	  err = MCMC_model1_it(0,force_acc,1,it,acc_h,it_tot,
			       i_D_orig,i_D_age,i_age,i_D_lga,i_length,i_D_lga_CC,i_length_CC,
			       i_D_g_a,i_D_g_a_CC,i_D_CC,i_D_COST,save_sim,write_alpha);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling MCMC_model1_it\n");
	      return(err);
	    }
	  err = update_average_age(it,i_age,i_D_age,i_age_mean);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling update_average_age\n");
	      return(err);
	    }

	  err = update_average_g_a(it,i_D_g_a);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling update_average_g_a\n");
	      return(err);
	    }
	  if(i_D_orig->coastal_cod)
	    {
	      err = update_average_g_a(it,i_D_g_a_CC);
	      if(err)
		{
		  write_warning("MCMC_model1:Error calling update_average_g_a\n");
		  return(err);
		}
	    }

	  err = update_average_lin(it,i_length,i_D_lga,i_length_mean);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling update_average_lin\n");
	      return(err);
	    }
	  if(i_D_orig->coastal_cod)
	    {
	      err = update_average_lin(it,i_length_CC,i_D_lga_CC,i_length_CC_mean);
	      if(err)
		{
		  write_warning("MCMC_model1:Error calling update_average_lin\n");
		  return(err);
		}
	    }

          err = Bayes_CV_model1(it,i_age,i_D_age,i_length,i_D_lga,
			      i_mod1_mean_inv_lik);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling Bayes_CV_age\n");
	      return(err);
	    }

          err = Bayes_CV_age(it,i_age,i_D_age,i_age_mean_inv_lik);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling Bayes_CV_age\n");
	      return(err);
	    }

          err = Bayes_CV_lin(it,i_length,i_D_lga,i_lga_mean_inv_lik);
	  if(err)
	    {
	      write_warning("MCMC_model1:Error calling Bayes_CV_age\n");
	      return(err);
	    }

	  if(i_D_orig->coastal_cod)
	    {
              #ifdef DEBUG_PROG
	      printf("Calculate s_lga_mean_inv_lik for coastal_cod\n");
	      #endif
	    }

          it++;
	  it_tot++;
          force_acc = 0;
	} // end for(iti=0;iti<i_inCommon->num_it_inner;iti++)
      
      /* Write samples to file */

      err = write_samples_model1(i_D_age,i_age,i_D_lga,i_length,i_D_g_a,
				 i_D_lga_CC,i_length_CC,i_D_g_a_CC,
				 i_D_CC,i_D_orig->coastal_cod,i_inCommon->print_boat,
				 i_inCommon->print_format,i_D_COST);
      if(err)
	{
	  write_warning("MCMC_model1:Error calling save_samples_model1\n");
	  return(err);
	}

      #ifdef DEBUG_PROG
      if(printAgeDist)
	{
	  double psum;
	  int pHaul=0,i,f;
	  int sumN=0;
	  double *p;
	  p = CALLOC(i_D_age->glm->ncat,double);
	  for(h=0;h<i_D_age->glm->nHaul;h++)
	    {
	      psum = G_ZERO;
	      for(a=0;a<i_D_age->glm->ncat;a++)
		{
		  p[a] = exp(i_age->alpha[h][a]);
		  psum += p[a];
		}
	      fprintf(fp,"%d %d",it_tot,h);
	      for(a=0;a<i_D_age->glm->ncat;a++)
		fprintf(fp," %f",p[a]/psum);
	      fprintf(fp,"\n");
	    }
	  FREE(p);
	}
      #endif

    } // end for(ito=0;ito<i_inCommon->num_it_outer;ito++)
  FREE(acc_h);
  if(i_inCommon->inc_hsz)
    {
      fclose(g_caa_mcmc_hsz);
    }
 
  #ifdef DEBUG_PROG
  if(printAgeDist)
    {
      fclose(fp);
    }
  #endif

  return(0);
}		/* end of MCMC_model1 */


/*!
  \author Geir Storvik
  \brief Main simulation steps inside each MCMC iteration for model1

  Simulations are performed by switching between the following steps:
  - Simulating age parameters
  - Simulating parameters of non-linear function in lga model
  - Simulating lga parameters

  If fixed_model, then simulate only age parameters.
*/
int MCMC_model1_it(int start_h,int i_force_acc,int i_len_only,
		   int i_it,int *o_acc_h,int i_it_tot,
		   Data_orig *i_D_orig, Data_age *i_D_age, Age_struct *i_age, 
		   Data_lin *i_D_lga, LW_struct *i_length, 
		   Data_lin *i_D_lga_CC, LW_struct *i_length_CC,
		   Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, Data_CC *i_D_CC,
		   Data_COST *i_D_COST, int i_save_sim, int i_write_alpha)
{
  int err=0;
  int nHaul;

  #ifdef DEBUG_PROG 
    printf("\nSampling Age-parameters\n");
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"\nSampling Age-parameters\n");
  #endif
  

  nHaul = i_D_age->glm->nHaul;
      
  err = MCMC_it_age(start_h,i_force_acc,i_len_only,i_it,i_it_tot,o_acc_h,
		    i_D_orig,i_D_age,i_age,i_D_lga,i_length,
		    i_D_lga_CC,i_length_CC,i_D_g_a,i_D_g_a_CC,i_D_CC,
		    i_D_COST,i_save_sim,nHaul,i_write_alpha);
  if(err)
    {
      write_warning("MCMC_model1_it:Error calling MCMC_it_age\n");
      return(err);
    }

  if(i_length->fixed_model == 0) /* Sample lga model */
    {
      #ifdef DEBUG_PROG
      printf("\nSampling Length-parameters\n");
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"\nSampling Length-parameters\n");
      #endif

      if(i_D_g_a->g_a_model>0)
	{
	  err = MCMC_it_g_a(i_D_orig,i_D_COST,i_D_age,i_D_lga,i_D_g_a,i_length,start_h,nHaul,i_it);
	  if(err)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_g_a\n");
              #endif
	      write_warning("MCMC_model1_it:Error calling MCMC_it_g_a\n");
	      return(err);
	    }
	  if(i_D_orig->coastal_cod)
	    {
	      err = MCMC_it_g_a(i_D_orig,i_D_COST,i_D_age,i_D_lga_CC,i_D_g_a_CC,i_length_CC,start_h,nHaul,i_it);
	      if(err)
		{
                  #ifdef LOG_FILE
		  fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_g_a\n");
                  #endif
		  write_warning("MCMC_model1_it:Error calling MCMC_it_g_a\n");
		  return(err);
		}
	    }
	}
      //if(i_D_COST->model)
      //err = MCMC_it_lga_COST(i_D_orig,i_D_COST,i_D_lga,i_D_g_a,i_length,start_h,i_it,nHaul,i_age->age_errors);
      //else
      err = MCMC_it_lga(i_D_lga,i_D_g_a,i_length,start_h,i_it,nHaul,i_age->age_errors,i_D_orig->coastal_cod);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_lga\n");
          #endif
	  write_warning("MCMC_model1_it:Error calling MCMC_it_lga\n");
	  return(err);
	}
      if(i_D_orig->coastal_cod)
	{
          #ifdef DEBUG_PROG
	  printf("\n Coastal cod\n");
	  printf("tau[%d][%d]=%lf\n",0,1,i_length_CC->par->tau[0][1]);
	  printf("tau_obs=%lf\n",i_length_CC->par->tau_obs);
          #endif
	  //i_length_CC->par->tau[0][1] = i_length->par->tau[0][1];
	  //i_length_CC->par->tau_obs = i_length->par->tau_obs;
	  err = MCMC_it_lga(i_D_lga_CC,i_D_g_a_CC,i_length_CC,start_h,i_it,nHaul,i_age->age_errors,i_D_orig->coastal_cod);
	  if(err)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_lga\n");
              #endif
	      write_warning("MCMC_model1_it:Error calling MCMC_it_lga\n");
	      return(err);
	    }
	}
    }
  else /* Use fixed lga model */
    {
      #ifdef DEBUG_PROG
      printf("\nSampling fixed Length-parameters\n");
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"\nSampling fixed Length-parameters\n");
      #endif
      // fixed ga parameters in MCMC_it_lga_fixed routine
      err = MCMC_it_lga_fixed(i_D_lga,i_D_g_a,i_length,start_h,i_it_tot);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_lga_fixed\n");
          #endif
	  write_warning("MCMC_model1_it:Error calling MCMC_it_lga_fixed\n");
	  return(err);
	}
      if(i_D_orig->coastal_cod)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\nSampling fixed Length-parameters\n");
          #endif
	  err = MCMC_it_lga_fixed(i_D_lga_CC,i_D_g_a_CC,i_length_CC,start_h,i_it_tot);
	  if(err)
	    {
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"\nMCMC_model1_it:Error calling MCMC_it_lga_fixed\n");
              #endif
	      write_warning("MCMC_model1_it:Error calling MCMC_it_lga_fixed\n");
	      return(err);
	    }
	}
    }

  return(0);
}		/* end of MCMC_model1_it */



/*!
  \author Geir Storvik
  \brief Main part for age model inside each MCMC iteration

  Simulations are performed by the following steps:
  - Simulating missing ages using the ::sample_ages routine
  - Simulating linear structure using the ::sample_gauss routine
  - Simulating ages if errors in age-readings ::sample_ages_age_error
  - Simulating haul effects using the ::sample_age_haul routine
  - Simulating precision parameters for haul effects
  - Calculating likelihood

  All data are now of the Amigo-type structure.
*/
int MCMC_it_age(int start_h,int i_force_acc,int i_len_only,
		int i_it,int i_it_tot,int *o_acc_h,
		Data_orig *i_D_orig, Data_age *i_D_age, Age_struct *i_age,
		Data_lin *i_D_lga, LW_struct *i_length, 
		Data_lin *i_D_lga_CC, LW_struct *i_length_CC,
		Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, Data_CC *i_D_CC,
		Data_COST *i_D_COST,
		int i_save_sim, int i_nHaul,int i_write_alpha)
{
  int err=0;

  /* Sample missing ages */
  #ifdef DEBUG_PROG
  printf("Sample missing ages\n");
  #endif
  err = sample_ages(i_D_orig,i_D_CC,i_age,i_D_age,i_length,i_D_lga,i_D_g_a,
		    i_length_CC,i_D_lga_CC,i_D_g_a_CC,i_save_sim,i_it);  
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_ages\n");
      return(err);
    }

  /* Sample discarded fish - ages and lengths */
  if(i_length->cens_model) 
    {
      /*if(i_D_COST->model)
	{
          #ifdef DEBUG_PROG
	  printf("\nSampling censoring parameters\n");
          #endif
	  err = sample_cens_par(i_D_lga,i_D_orig,i_D_COST);
	  if(err)
	    {
	      write_warning("MCMC_it_age:Error calling sample_cens_par\n");
	      return(err);
	    }
          #ifdef DEBUG_PROG
	  printf("\nSampling discarded fish\n");
          #endif
	  err = sample_discard_COST(i_age,i_D_age,i_length,i_D_lga,i_D_g_a,i_D_orig,i_D_COST,i_save_sim);
	  if(err)
	    {
	      write_warning("MCMC_it_age:Error calling sample_discarded_COST\n");
	      return(err);
	    }
	    }*/
      //      else /* Sample discarded fish, old version not changed yet */
	{
	  printf("MCMC_it_age: sample_discarded: must change to continuous age\n");
	}	    
    }  
  
  /* Calculate sufficient statistics */
  #ifdef DEBUG_PROG
  printf("Calculate sufficient statistics\n");
  #endif
  err = make_suff_age(i_D_age->glm->ncat,i_age,i_D_age,i_D_lga->haulweight,start_h);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_age:Error calling make_suff_age\n");
      #endif
      write_warning("MCMC_it_age:Error calling make_suff_age\n");
      return(err);
    }

  // Remove first haul effect from sampling
  /* Sample fixed and remaining random effects */
  #ifdef DEBUG_PROG
  printf("Sample fixed and random effects (not haul effect)\n");
  #endif
  i_D_age->glm->xcov[0]->n_cov--;
  err = sample_gauss_eff(i_age->gr_str_f,i_age->par,i_D_age->glm,start_h,i_D_age->glm->nHaul,0);
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_gauss_eff\n");
      return(err);
    }

  // Add haul effect again
  i_D_age->glm->xcov[0]->n_cov++;

  /* Sample complete alpha's */
  #ifdef DEBUG_PROG
  printf("Sample haul effect\n");
  #endif
  err = sample_age_alpha(i_age,i_D_age,start_h,i_force_acc,i_it,o_acc_h,i_write_alpha);
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_age_alpha\n");
      return(err);
    }

  /* Sample precision parameters */
  #ifdef DEBUG_PROG
  printf("Sample precision parameters\n");
  #endif
  err = sample_precision_age_haul(start_h,i_age->par,i_age->alpha,i_D_age->glm,i_nHaul,i_age->gr_str_f->constr->nc);
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_precision_haul\n");
      return(err);
    }

  i_D_age->glm->xcov[0]->n_cov--;
  err = sample_precision(start_h,i_age->par,i_D_age->glm,i_D_age->glm->nHaul,0);
  i_D_age->glm->xcov[0]->n_cov++;
  if(err)
    {
      write_warning("MCMC_it_age:Error calling sample_precision\n");
      return(err);
    }

  /* Calculate likelihoods */
  err = calc_lik_age(i_age,i_D_age,&(i_age->par->loglik));
  if(err)
    {
      write_warning("MCMC_model1:Error calling calc_lik_age\n");
      return(err);
    }

  #ifdef DEBUG_PROG
  int a,h,i,ind_a,sum,suma;
  double *p,psum;

  p = CALLOC(i_D_age->glm->ncat,double);
  printf("MCMC_it_age:\n");
  psum = G_ZERO;
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      p[a] = exp(i_age->par->eff[a][0][0][0]);
      psum += p[a];
    }
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      printf("alpha_const=");
      for(i=0;i<i_D_age->glm->nxcov;i++)
	printf("%lf ",i_age->par->eff[a][i][0][0]);
      printf("p=%lf\n",p[a]/psum);
    }
  printf("sigma_haul=%lf\n",sqrt(G_ONE/i_age->par->tau_obs));
  sum = 0;
  printf("Num_age= ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      suma = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	suma += (int) i_D_age->Ages[h][a];
      printf("%d ",suma);
      sum += suma;
    }
  printf("\n");
  printf("Total aged fish=%d\n",sum);
  printf("Loglik_age=%f\n",i_age->par->loglik);
  printf("tau_obs=%f\n",i_age->par->tau_obs);
  FREE(p);
  #endif

  #ifdef LOG_FILE
  fprintf(g_caa_log,"\n");
  fprintf(g_caa_log,"Loglik_age=%f\n",i_age->par->loglik);
  fprintf(g_caa_log,"tau_obs=%f\n",i_age->par->tau_obs);
  #endif

  return(0);
}		/* end of MCMC_it_age */



/*!
 \fn int MCMC_it_lga(Data_lin *i_D_lga,Data_g_a *i_D_g_a,LW_struct *i_length,int i_start_h)
 \author Geir Storvik
  \brief Main part for lga model inside each MCMC iteration

  Simulations are performed by the following steps:
  - Calculating sufficient statistics for the lga model
  - Simulating linear structure using the ::sample_gauss routine
  - Simulating precision parameters
  - Calculating likelihood
*/
int MCMC_it_lga(Data_lin *i_D_lga,Data_g_a *i_D_g_a,LW_struct *i_length,int i_start_h,
		int i_it,int i_nHaul,int i_age_error,int i_coastal_cod)
{
  int err;

  /* Sample effects except haul effect */
  /* based on observed age-length data */

  if(i_age_error || i_D_g_a->g_a_model>0 || i_coastal_cod) 
    {
      /* If ages observed with errors or non-linear model, new sufficient statistics */
      /* Calculate sufficient statistics - using observed ages */
      #ifdef DEBUG_PROG
      printf("MCMC_it_lga:make new sufficient statistics\n");
      #endif
      err = make_suff_lga(i_D_lga,i_D_g_a,i_start_h,0);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"MCMC_it_lga:Error calling make_suff_lga\n");
          #endif
	  write_warning("MCMC_it_lga:Error calling make_suff_lga\n");
	  return(err);
	}
    }
  else //ages observed without errors
    {
      /* Copy sufficient statistics - using only observed ages */
      #ifdef DEBUG_PROG
      printf("MCMC_it_lga:copy sufficient statistics\n");
      #endif
      err = copy_suff_lga_fix(i_D_lga,i_start_h);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"MCMC_it_lga:Error calling copy_suff_lga\n");
          #endif
	  write_warning("MCMC_it_lga:Error calling copy_suff_lga\n");
	  return(err); 
	}
    }
  err = sample_gauss_eff(i_length->gr_str,i_length->par,i_D_lga->glm,i_start_h,i_nHaul,1);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_lga:Error calling sample_gauss_eff\n");
      #endif
      write_warning("MCMC_it_lga:Error calling sample_gauss_eff\n");
      return(err);
    }
  #ifdef DEBUG_PROG
     printf("MCMC_it_lga:\n");
     printf("Int=%lf\n",i_length->par->eff[0][0][0][0]);
     printf("Slope=%lf\n",i_length->par->eff[0][1][0][0]);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Int=%lf\n",i_length->par->eff[0][0][0][0]);
  fprintf(g_caa_log,"Slope=%lf\n",i_length->par->eff[0][1][0][0]);
  #endif

  /* Sample haul effects using observed and simulated data */
  if(i_D_lga->glm->xcov[0]->ihaul>0)
    {
      /* Calculate sufficient statistics - using simulated ages */
      err = make_suff_lga(i_D_lga,i_D_g_a,i_start_h,1);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"MCMC_it_lga:Error calling make_suff_lga\n");
          #endif
	  write_warning("MCMC_it_lga:Error calling make_suff_lga\n");
	  return(err);
	}
      err = sample_lga_haul(i_start_h,i_length->par,i_D_lga->glm);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"MCMC_it_lga:Error calling sample_lga_haul\n");
          #endif
	  write_warning("MCMC_it_lga:Error calling sample_lga_haul\n");
	  return(err);
	}

      /* Calculate sufficient statistics again - using observed ages */
      if(i_age_error || i_D_g_a->g_a_model>0) 
	{
	  /* If ages observed with errors or non-linear model, new sufficient statistics */
	  err = make_suff_lga(i_D_lga,i_D_g_a,i_start_h,0);
	  if(err)
	    {
          #ifdef LOG_FILE
	      fprintf(g_caa_log,"MCMC_it_lga:Error calling make_suff_lga\n");
          #endif
	      write_warning("MCMC_it_lga:Error calling make_suff_lga\n");
	      return(err);
	    }
	}
      else //ages observed without errors
	{
	  /* Copy sufficient statistics - using only observed ages */
	  err = copy_suff_lga_fix(i_D_lga,i_start_h);
	  if(err)
	    {
          #ifdef LOG_FILE
	      fprintf(g_caa_log,"MCMC_it_lga:Error calling copy_suff_lga\n");
          #endif
	      write_warning("MCMC_it_lga:Error calling copy_suff_lga\n");
	      return(err); 
	    }
	}
    }

  /* Sample precision parameter tau_obs */
  err = sample_precision_lin(i_start_h,i_length->par,i_D_lga->glm,i_nHaul);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_lga:Error calling sample_precision_lin\n");
      #endif
      write_warning("MCMC_it_lga:Error calling sample_precision_lin\n");
      return(err);
    } 

  /* Sample rest of precision parameters */
  if(i_it>0)
    {
      err = sample_precision(i_start_h,i_length->par,i_D_lga->glm,i_nHaul,1);
      if(err)
	{
	  write_warning("MCMC_it_lga:Error calling sample_precision\n");
	  return(err);
	}
    }

  #ifdef DEBUG_PROG
  int i,j,isp;
  Data_cov *xcov;
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      xcov = i_D_lga->glm->xcov[i];
      isp = xcov->ispat;
      for(j=0;j<xcov->n_cov;j++)
	if(!xcov->fix[j] && j != isp && j != xcov->icell)
	  printf("tau[%d][%d]=%lf\n",i,j,i_length->par->tau[i][j]);
    }
  printf("tau_obs=%lf\n",i_length->par->tau_obs);
  #endif

  /* Calculate likelihood */
  err = calc_lik_lin(i_length,i_D_lga,&(i_length->par->loglik));
  if(err)
    {
      write_warning("MCMC_model1:Error calling calc_lik_lin\n");
      return(err);
    }
  #ifdef DEBUG_PROG
  printf("Loglik_lga=%lf\n",i_length->par->loglik);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Loglik_lga=%lf\n",i_length->par->loglik);
  #endif



  return(0);
}		/* end of MCMC_it_lga */



/*!
  \author Hanne Rognebakke
  \brief Nonlinear connection in the lga model inside each MCMC iteration for COST program
*/
int MCMC_it_lga_COST(Data_orig *i_D_orig, Data_COST *i_D_COST, Data_lin *i_D_lga, Data_g_a *i_D_g_a,
		     LW_struct *i_length, int i_start_h, int i_it, int i_nHaul, int i_age_error)
{
  int err;

  /* Calculate sufficient statistics */
  /*err = make_suff_lga_COST(i_D_orig,i_D_COST,i_D_lga,i_D_g_a,i_start_h,0);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_lga:Error calling make_suff_lga\n");
      #endif
      write_warning("MCMC_it_lga:Error calling make_suff_lga\n");
      return(err);
      }*/

  err = sample_gauss_eff(i_length->gr_str,i_length->par,i_D_lga->glm,i_start_h,i_nHaul,1);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_lga:Error calling sample_gauss_eff\n");
      #endif
      write_warning("MCMC_it_lga:Error calling sample_gauss_eff\n");
      return(err);
    }
  #ifdef DEBUG_PROG
     printf("MCMC_it_lga:\n");
     printf("Int=%lf\n",i_length->par->eff[0][0][0][0]);
     printf("Slope=%lf\n",i_length->par->eff[0][1][0][0]);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Int=%lf\n",i_length->par->eff[0][0][0][0]);
  fprintf(g_caa_log,"Slope=%lf\n",i_length->par->eff[0][1][0][0]);
  #endif

  /* Sample haul effects */
  if(i_D_lga->glm->xcov[0]->ihaul>0)
    {
      err = sample_lga_haul(i_start_h,i_length->par,i_D_lga->glm);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"MCMC_it_lga:Error calling sample_lga_haul\n");
          #endif
	  write_warning("MCMC_it_lga:Error calling sample_lga_haul\n");
	  return(err);
	}
    }
  
  /* Sample precision parameter tau_obs */
  err = sample_precision_lin(i_start_h,i_length->par,i_D_lga->glm,i_nHaul);
  if(err)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"MCMC_it_lga:Error calling sample_precision_lin\n");
      #endif
      write_warning("MCMC_it_lga:Error calling sample_precision_lin\n");
      return(err);
    } 

  /* Sample rest of precision parameters */
  if(i_it>0)
    {
      err = sample_precision(i_start_h,i_length->par,i_D_lga->glm,i_nHaul,1);
      if(err)
	{
	  write_warning("MCMC_it_lga:Error calling sample_precision\n");
	  return(err);
	}
    }

  #ifdef DEBUG_PROG
  int i,j,isp;
  Data_cov *xcov;
  for(i=0;i<i_D_lga->glm->nxcov;i++)
    {
      xcov = i_D_lga->glm->xcov[i];
      isp = xcov->ispat;
      for(j=0;j<xcov->n_cov;j++)
	if(!xcov->fix[j] && j != isp && j != xcov->icell)
	  printf("tau[%d][%d]=%lf\n",i,j,i_length->par->tau[i][j]);
    }
  printf("tau_obs=%lf\n",i_length->par->tau_obs);
  #endif

  /* Calculate likelihood */
  err = calc_lik_lin(i_length,i_D_lga,&(i_length->par->loglik));
  if(err)
    {
      write_warning("MCMC_model1:Error calling calc_lik_lin\n");
      return(err);
    }
  #ifdef DEBUG_PROG
  printf("Loglik_lga=%lf\n",i_length->par->loglik);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Loglik_lga=%lf\n",i_length->par->loglik);
  #endif


  return(0);
}		/* end of MCMC_it_lga_COST */



/*!
  \author Geir Storvik
  \brief Nonlinear connection in the lga model inside each MCMC iteration

  Simulations are performed by the following steps:
  - Calculate sufficiednt statistics for the non-linear model
  - Simulating parameters in g-function
*/
int MCMC_it_g_a(Data_orig *i_D_orig,Data_COST *i_D_COST,Data_age *i_D_age,Data_lin *i_D_lga,
		Data_g_a *i_D_g_a,LW_struct *i_length,int i_start_h,int i_nHaul,int i_it)
{
  int err=0; 

  if(i_D_g_a->g_a_model==1)/* Schnute Richards model */
    {
      /* Calculate sufficient statistics */
      #ifdef DEBUG_PROG
      printf("\nCalculate sufficient parameters for sampling g_a-parameters\n");
      #endif
      //if(i_D_COST->model)
      //err = suff_g_a_COST(i_D_orig,i_D_COST,i_length,i_D_age,i_D_lga,i_D_g_a,i_start_h,i_nHaul,i_D_g_a->suff);
      //else
      err = suff_g_a(i_length,i_D_age,i_D_lga,i_D_g_a,i_start_h,i_nHaul,i_D_g_a->suff);
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"MCMC_it_g_a:Error calling suff_g_a\n");
          #endif
	  write_warning("MCMC_it_g_a:Error calling suff_g_a\n");
	  return(err);
	}
      
      #ifdef DEBUG_PROG
      printf("\nSampling g_a-parameters\n");
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"\nSampling g_a-parameters\n");
      #endif
      /* Sample nonlinear function g(a) */
      err = sample_g_a(i_length,i_D_age,i_D_lga,i_D_g_a,i_D_g_a->suff,i_it,i_nHaul);
      #ifdef LOG_FILE
      fprintf(g_caa_log,"\nSampling g_a-parameters: %d \n",err);
      #endif
      if(err)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"MCMC_it_g_a:Error calling sample_g_a\n");
          #endif
	  write_warning("MCMC_it_g_a:Error calling sample_g_a\n");
	  return(err);
	}
    }
  else if(i_D_g_a->g_a_model==2)/* polyn3 model */
    {
      printf("MCMC_it_g_a:Update age-length relation: polyn3 model\n");
      return(err);
    }
  else 
    {
      printf("MCMC_it_g_a:Unknown age-length relation\n");
      return(err);
    }

  return(0);
}		/* end of MCMC_it_g_a */



/*!
  \author Geir Storvik
  \brief  Perform MCMC simulations of wgl model

  Start with initialization of parameters.

  The simulations are divided into burn-in and simulations after burnin in which
  num_it_inner simulations are performed for each num_it_outer simulation.
  All inner simulations are performed through the ::MCMC_model1_it routine.

  Simulations are saved for each num_it_outer iteration using the ::write_it
  routine.
*/
int MCMC_model2(Input_common *i_inCommon,
		Data_lin *i_D_wgl, LW_struct *i_weight, LW_struct *i_weight_mean,
		Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, LW_struct *i_weight_CC_mean,
		int i_coastal_cod,
		double *i_wgl_mean_inv_lik)
{
  int        err;
  int        it,iti,ito,it_tot;

 
  it_tot=0;
  /* Burn in */
  for(iti=0;iti<i_inCommon->burn_in;iti++)
    {
      err = MCMC_model2_it(0,it_tot,i_D_wgl,i_weight,
			   i_D_wgl_CC,i_weight_CC,i_coastal_cod);
      it_tot++;
      if(err)
	{
	  write_warning("MCMC_model2:Error calling MCMC_model2_it\n");
	  return(err);
	}
    }

  /* Start full simulation */
  it = 0;
  for(ito=0;ito<i_inCommon->num_it_outer;ito++)
    {
      for(iti=0;iti<i_inCommon->num_it_inner;iti++)
	{
          #ifdef DEBUG_PROG
	  printf("\n\nIteration %d\n",it);
          #endif
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"\n\nIteration %d\n",it);
          #endif
	  err = MCMC_model2_it(0,it_tot,i_D_wgl,i_weight,
			       i_D_wgl_CC,i_weight_CC,i_coastal_cod);
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling MCMC_model2_it\n");
	      return(err);
	    }

	  err = update_average_lin(it,i_weight,i_D_wgl,i_weight_mean);
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling update_average_lin\n");
	      return(err);
	    }

	  if(i_coastal_cod)
	    {
	      err = update_average_lin(it,i_weight_CC,i_D_wgl_CC,i_weight_CC_mean);
	      if(err)
		{
		  write_warning("MCMC_model2:Error calling update_average_lin\n");
		  return(err);
		}
	    }

          err = Bayes_CV_lin(it,i_weight,i_D_wgl,i_wgl_mean_inv_lik);
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling Bayes_CV_age\n");
	      return(err);
	    }

          it++;
	  it_tot++;
	}
      if(i_inCommon->inc_hsz)
	{
	  i_weight->par->eff[0][1][0][0] = 0.0; // This is undetermined in the haulsize model
	}

      /* Calculate likelihoods */
      err = calc_lik_lin(i_weight,i_D_wgl,&(i_weight->par->loglik));
      if(err)
	{
	  write_warning("MCMC_model2:Error calling calc_lik_lin\n");
	  return(err);
	}
      if(i_coastal_cod)
	{
	  err = calc_lik_lin(i_weight_CC,i_D_wgl_CC,&(i_weight_CC->par->loglik));
	  if(err)
	    {
	      write_warning("MCMC_model2:Error calling calc_lik_lin\n");
	      return(err);
	    }
	}
      #ifdef DEBUG_PROG
      printf("Loglik_wgl=%lf\n",i_weight->par->loglik);
      printf("Int=%lf\n",i_weight->par->eff[0][0][0][0]);
      printf("Slope=%lf\n",i_weight->par->eff[0][1][0][0]);
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"Loglik_wgl=%lf\n",i_weight->par->loglik);
      fprintf(g_caa_log,"Int=%lf\n",i_weight->par->eff[0][0][0][0]);
      fprintf(g_caa_log,"Slope=%lf\n",i_weight->par->eff[0][1][0][0]);
      #endif

      /* Save samples */
      err = write_samples_model2(i_D_wgl,i_weight,i_D_wgl_CC,i_weight_CC,
				 i_coastal_cod,i_inCommon->print_boat,i_inCommon->print_format);
      if(err)
	{
	  write_warning("MCMC_model2:Error calling write_samples_model2\n");
	  return(err);
	}
      if(i_inCommon->inc_hsz) //write only parameters that are needed for the estimation in a separate file
      	{
	  err = write_samples_haulsize(g_caa_mcmc_hsz,i_D_wgl->glm,i_weight->par);
	}

    }  // end for(ito=0;ito<i_inCommon->num_it_outer;ito++)

  return(0);
}		/* end of MCMC_model2 */




/*!
  \author Geir Storvik
  \brief Main simulation steps inside each MCMC iteration for model2

  Simulations are performed by the following step:
  - Simulating wgl parameters
*/
int MCMC_model2_it(int start_h, int it_tot, 		
		   Data_lin *i_D_wgl, LW_struct *i_weight, 
		   Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC,
		   int i_coastal_cod)
{
  int err;

  if(i_weight->fixed_model == 0) /* Sample wgl model */
    {
      #ifdef DEBUG_PROG
      printf("Sampling Weight-parameters\n");
      #endif
      #ifdef LOG_FILE
      fprintf(g_caa_log,"Sampling Weight-parameters\n");
      #endif
      err = MCMC_it_wgl(i_D_wgl,i_weight,start_h);
      if(err)
	{
	  write_warning("MCMC_model2_it:Error calling MCMC_it_wgl\n");
	  return(err);
	}
      if(i_coastal_cod)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"Sampling Weight-parameters: Coastal cod\n");
          #endif
	  err = MCMC_it_wgl(i_D_wgl_CC,i_weight_CC,start_h);
	  if(err)
	    {
	      write_warning("MCMC_model2_it:Error calling MCMC_it_wgl\n");
	      return(err);
	    }
	}
    }  
  else
    {
      i_weight->par->tau_obs = i_weight->fixed_tau[it_tot];
      i_weight->par->eff[0][0][0][0] = i_weight->fixed_int[it_tot]; /* Intercept */
      i_weight->par->eff[0][1][0][0] = i_weight->fixed_slp[it_tot]; /* Slope */
      if(i_coastal_cod)
	{
	  i_weight_CC->par->tau_obs = i_weight_CC->fixed_tau[it_tot];
	  i_weight_CC->par->eff[0][0][0][0] = i_weight_CC->fixed_int[it_tot]; /* Intercept */
	  i_weight_CC->par->eff[0][1][0][0] = i_weight_CC->fixed_slp[it_tot]; /* Slope */
	}
    }


  return(0);
}		/* end of MCMC_model2_it */



/*!
  \author Geir Storvik
  \brief Main part for lga model inside each MCMC iteration

  Simulations are performed by the following steps:
  - Simulating linear structure using the ::sample_gauss routine
  - Simulating precision parameters
*/
int MCMC_it_wgl(Data_lin *i_D_wgl,LW_struct *i_weight,int start_h)
{
  int err;

  /* Sample effects */
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Sample effects\n");
  #endif
  err = sample_gauss_eff(i_weight->gr_str,i_weight->par,i_D_wgl->glm,start_h,i_D_wgl->glm->nHaul,0);
  if(err)
    {
      write_warning("MCMC_it_wgl:Error calling sample_gauss_eff\n");
      return(err);
    }
 
  /* Sample precision parameter tau_obs */
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Sample precision tau_obs\n");
  #endif
  err = sample_precision_lin(start_h,i_weight->par,i_D_wgl->glm,i_D_wgl->glm->nHaul);
  if(err)
    {
      write_warning("MCMC_it_wgl:Error calling sample_precision_lin\n");
      return(err);
    }

  /* Sample rest of precision parameters */	  
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Sample rest of precisions\n");
  #endif
  err = sample_precision(start_h,i_weight->par,i_D_wgl->glm,i_D_wgl->glm->nHaul,0);
  if(err)
    {
      write_warning("MCMC_it_wgl:Error calling sample_precision\n");
      return(err);
    }

  return(0);
}		/* end of MCMC_it_wgl */



/*!
  \author 
  \brief Main part for lga model inside each MCMC iteration if fixed lga model

  Simulations are performed by the following steps:
  - Calculating sufficient statistics for the lga model
  - Simulating lga parameters (Int, Slp, precision(fish)) using the ::sample_ routine 
  - Calculating likelihood
*/
int MCMC_it_lga_fixed(Data_lin *i_D_lga,Data_g_a *i_D_g_a,LW_struct *i_length,int start_h,int i_it)
{
  int err;

  /* "sample" int, slp and tau_obs from a set of realisations */ 
  i_length->par->eff[0][0][0][0] = i_length->fixed_int[i_it];     
  i_length->par->eff[0][1][0][0] = i_length->fixed_slp[i_it];     
  i_length->par->tau_obs = i_length->fixed_tau[i_it]; 

  if(i_D_g_a->g_a_model==1)
    {
      i_D_g_a->g_a_par[0] = i_D_g_a->fixed_c[i_it];
      i_D_g_a->g_a_par[1] = i_D_g_a->fixed_theta[i_it];
      i_D_g_a->g_a_par[2] = i_D_g_a->fixed_gamma[i_it];
    }

  #ifdef DEBUG_PROG
  printf("MCMC_it_lga:\n");
  printf("Int=%lf\n",i_length->par->eff[0][0][0][0]);
  printf("Slope=%lf\n",i_length->par->eff[0][1][0][0]);
  printf("Precision=%lf\n",i_length->par->tau_obs);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Int=%lf\n",i_length->par->eff[0][0][0][0]);
  fprintf(g_caa_log,"Slope=%lf\n",i_length->par->eff[0][1][0][0]);
  fprintf(g_caa_log,"Precision=%lf\n",i_length->par->tau_obs);
  #endif

  err = calc_lik_lin(i_length,i_D_lga,&(i_length->par->loglik));
  if(err)
    {
      write_warning("MCMC_model1:Error calling calc_lik_lin\n");
      return(err);
    }
  #ifdef DEBUG_PROG
     printf("Loglik_lga=%lf\n",i_length->par->loglik);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Loglik_lga=%lf\n",i_length->par->loglik);
  #endif

  return(0);
}		/* end of MCMC_it_lga_fixed */



