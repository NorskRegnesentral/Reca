/*!
  \file caa_init.c
  \brief Routines for initialization of variables and structs
  \author Geir Storvik

*/
#include "caa.h"
#include "caa_init.h"
#include "caa_read_write.h"
#include "caa_sample_g_a.h"
#include "caa_sample_multi.h"
#include "caa_sample_gauss.h"
#include "caa_evaluate.h"
#include "caa_routines.h"


double Qfunc_age_fix(int node, int nnode, char *arg);
double Qfunc_age_ran(int node, int nnode, char *arg);
double Qfunc_length(int node, int nnode, char *arg);
double Qfunc_weight(int node, int nnode,char *arg);

#ifdef LOG_FILE
extern FILE     *g_caa_log;
#endif


/*!
  \brief Memory allocation and initial values for age and lga models
  \author Geir Storvik
*/
int initialize_model1(int i_seed, int *i_num_par,
		      Data_age *i_D_age, Age_struct **o_age, Age_struct **o_age_mean,
		      Data_lin *i_D_lga, LW_struct **o_length, LW_struct **o_length_mean,
		      Data_lin *i_D_lga_CC, LW_struct **o_length_CC, LW_struct **o_length_CC_mean,
		      Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, 
		      int i_age_errors, double *i_A2A, double i_delta_age,
		      int i_coastal_cod,int i_inc_hsz,int i_sim_ar,
		      double *i_pri_age_eff_mean, double *i_pri_age_eff_prec,
		      double *i_pri_age_prec_par, double *i_pri_age_ar,
		      double *i_pri_lga_eff_mean, double *i_pri_lga_eff_prec,
		      double *i_pri_lga_prec_par, double *i_pri_lga_ar,
		      int i_lga_fixed_model, double *i_lga_fixed_int, 
		      double *i_lga_fixed_slp, double *i_lga_fixed_tau,
		      double *i_lga_fixed_g_a_c,double *i_lga_fixed_g_a_theta, double *i_lga_fixed_g_a_gamma,
		      double *i_lga_fixed_int_CC, double *i_lga_fixed_slp_CC, double *i_lga_fixed_tau_CC,
		      double *i_lga_fixed_g_a_c_CC,double *i_lga_fixed_g_a_theta_CC, double *i_lga_fixed_g_a_gamma_CC,
		      int i_lga_cens_model, double *i_lga_cens_par)
{
  int err=0;
  Age_struct *age, *age_mean;
  LW_struct *length, *length_mean, *length_CC, *length_CC_mean;

  (*GMRFLib_uniform_init)(i_seed);
  #ifdef DEBUG_PROG
  printf("Allocating memory for age structs\n");
  #endif
  /* Allocating memory for age structs */
  err = init_age(i_D_age,i_age_errors,i_A2A,i_coastal_cod,i_inc_hsz,i_sim_ar,
		 i_pri_age_eff_mean,i_pri_age_eff_prec,i_pri_age_prec_par,i_pri_age_ar,
		 &age,&age_mean);
  if(err)
    {
      write_warning("initialize_model1:Error calling init_age\n");
      return(err);
    }
  age->par->num_var = i_num_par[0];
  age->gr_str_f->Qfunc = Qfunc_age_fix;
  age->gr_str_r->Qfunc = Qfunc_age_ran;
  age->delta_age = i_delta_age;


  #ifdef DEBUG_PROG
  printf("Allocating memory for length structs\n");
  #endif
  /* Allocating memory for length structs */
  err = init_lin(i_D_lga,i_sim_ar,i_pri_lga_eff_mean,i_pri_lga_eff_prec,
		 i_pri_lga_prec_par,i_pri_lga_ar,i_lga_fixed_model,
		 &length,&length_mean);
  if(err)
    {
      write_warning("initialize_model1:Error calling init_lin\n");
      return(err);
    }
  length->par->num_var = i_num_par[1];
  length->gr_str->Qfunc = Qfunc_length;

  length->fixed_model = i_lga_fixed_model;
  if(length->fixed_model)
    {
      length->fixed_int = i_lga_fixed_int;
      length->fixed_slp = i_lga_fixed_slp;
      length->fixed_tau = i_lga_fixed_tau;
      if(i_D_g_a->g_a_model==1)
	{
	  i_D_g_a->fixed_c = i_lga_fixed_g_a_c;
	  i_D_g_a->fixed_theta = i_lga_fixed_g_a_theta;
	  i_D_g_a->fixed_gamma = i_lga_fixed_g_a_gamma;
	  i_D_g_a->g_a_par[0] = i_lga_fixed_g_a_c[0];
	  i_D_g_a->g_a_par[1] = i_lga_fixed_g_a_theta[0];
	  i_D_g_a->g_a_par[2] = i_lga_fixed_g_a_gamma[0];
	}
      else if(i_D_g_a->g_a_model>1)
	{
	  write_warning("initialize_model1:Error using fixed model. Unknown g-function.\n");
	  return(1);
	}
    }
 
  /* Allocating memory for extra length structs if coastal cod */
  if(i_coastal_cod)
    {
      #ifdef DEBUG_PROG
      printf("Allocating memory for length structs - coastal cod\n");
      #endif
      err = init_lin(i_D_lga_CC,i_sim_ar,i_pri_lga_eff_mean,i_pri_lga_eff_prec,
		     i_pri_lga_prec_par,i_pri_lga_ar,i_lga_fixed_model,
		     &length_CC,&length_CC_mean);
      if(err)
	{
	  write_warning("initialize_model1:Error calling init_lin\n");
	  return(err);
	}
      length_CC->par->num_var = i_num_par[3];
      length_CC->gr_str->Qfunc = Qfunc_length;
      length_CC->fixed_model = i_lga_fixed_model;
      if(length_CC->fixed_model)
	{
	  length_CC->fixed_int = i_lga_fixed_int_CC;
	  length_CC->fixed_slp = i_lga_fixed_slp_CC;
	  length_CC->fixed_tau = i_lga_fixed_tau_CC;
	  if(i_D_g_a_CC->g_a_model==1)
	    {
	      i_D_g_a_CC->fixed_c = i_lga_fixed_g_a_c_CC;
	      i_D_g_a_CC->fixed_theta = i_lga_fixed_g_a_theta_CC;
	      i_D_g_a_CC->fixed_gamma = i_lga_fixed_g_a_gamma_CC;
	      i_D_g_a_CC->g_a_par[0] = i_lga_fixed_g_a_c[0];
	      i_D_g_a_CC->g_a_par[1] = i_lga_fixed_g_a_theta[0];
	      i_D_g_a_CC->g_a_par[2] = i_lga_fixed_g_a_gamma[0];
	    }
	}
    }

  /* Initial values if sampling discards */
  length->cens_model = i_lga_cens_model;
  if(length->cens_model)
    {
      #ifdef DEBUG_PROG
      printf("Initial values if sampling discards\n");
      #endif
      length->cens_k = i_lga_cens_par[0];
      length->cens_m = i_lga_cens_par[1];
      length->cens_r = i_lga_cens_par[2];
      length->cens_Nlim = i_lga_cens_par[3];
    }
  if(i_coastal_cod)
    {
      length_CC->cens_model = i_lga_cens_model;
      if(length_CC->cens_model)
	{
	  length_CC->cens_k = i_lga_cens_par[0];
	  length_CC->cens_m = i_lga_cens_par[1];
	  length_CC->cens_r = i_lga_cens_par[2];
	  length_CC->cens_Nlim = i_lga_cens_par[3];
	}
    }

  /* Initial values for length given age model, needed sampling ages */
  #ifdef DEBUG_PROG
  printf("Initial values for length given age model\n");
  #endif
  err = init_lin_par(length,i_D_lga);
  if(err)
    {
      write_warning("initialize_model1:Error calling init_lin_par\n");
      return(err);
    }
  if(i_coastal_cod)
    {
      err = init_lin_par(length_CC,i_D_lga_CC);
      if(err)
	{
	  write_warning("initialize_model1:Error calling init_lin_par\n");
	  return(err);
	}
    }

  err = init_evaluate(i_D_age->glm->nHaul,i_D_age->glm->ncat);
  if(err)
    {
      write_warning("initialize_model1:Error calling init_evaluate\n");
      return(err);
    }

  // Initialization for routines in caa_sample_multi.c
  #ifdef DEBUG_PROG
  printf("Initialization for routines in caa_sample_multi.c\n");
  #endif
  err = sample_multi_initialize(i_D_age->glm->ncat);
  if(err)
    {
      write_warning("initialize_model1:Error calling sampling_multi_initialize\n");
      return(err);
    }

  // Initialization for routines in caa_sample_g_a.c
  if(i_D_g_a->g_a_model>0)
    {
      err = sample_g_a_initialize(i_D_g_a->ncat,i_D_g_a->g_a_model);
      if(err)
	{
	  write_warning("initialize_model1:Error calling sampling_g_a_initialize\n");
	  return(err);
	}
    }

  *o_age = age;
  *o_age_mean = age_mean;
  *o_length = length;
  *o_length_mean = length_mean;
  *o_length_CC = length_CC;
  *o_length_CC_mean = length_CC_mean;

  return(0);
}		/* end of initialize_model1 */



/*!
  \brief Reallocate memory allocated in initialize_model1
  \author Geir Storvik
*/
int re_initialize_model1(Data_age *i_D_age, Age_struct *i_age, Age_struct *i_age_mean, 
			 Data_lin *i_D_lga, LW_struct *i_length, LW_struct *i_length_mean,
			 Data_lin *i_D_lga_CC, LW_struct *i_length_CC, LW_struct *i_length_CC_mean,
			 Data_g_a *i_D_g_a, int i_coastal_cod, int i_age_errors, double *i_A2A, int i_inc_hsz)
{
  int  err;

  err = sample_multi_re_initialize();
  if(err)
    {
      write_warning("re_initialize_model1:Error calling sampling_multi_re_initialize\n");
      return(err);
    }
  if(i_D_g_a->g_a_model>0)
    {
      err = sample_g_a_re_initialize();
      if(err)
	{
	  write_warning("re_initialize_model1:Error calling sampling_g_a_re_initialize\n");
	  return(err);
	}
    }

  err = re_init_age(i_D_age,i_age_errors,i_A2A,i_inc_hsz,&i_age,&i_age_mean);
  if(err)
    {
      write_warning("re_initialize_model1:Error calling re_init_age\n");
      return(err);
    }

  err = re_init_lin(i_D_lga,&i_length,&i_length_mean);
  if(err)
    {
      write_warning("re_initialize_model1:Error calling re_init_lin\n");
      return(err);
    }

  if(i_coastal_cod)
    {
      err = re_init_lin(i_D_lga_CC,&i_length_CC,&i_length_CC_mean);
      if(err)
	{
	  write_warning("re_initialize_model1:Error calling re_init_lin\n");
	  return(err);
	}
    }

  err = re_init_evaluate(i_D_age->glm->ncat);
  if(err)
    {
      write_warning("re_initialize_model1:Error calling init_evaluate\n");
      return(err);
    }


  return(0);
}		/* end of re_initialize_model1 */



/*!
  \brief Memory allocation and initial values for wgl model
  \author Geir Storvik
*/
int initialize_model2(int i_seed, int *i_num_par,
		      Data_lin *i_D_wgl, LW_struct **o_weight, LW_struct **o_weight_mean,
		      Data_lin *i_D_wgl_CC, LW_struct **o_weight_CC, LW_struct **o_weight_CC_mean,
		      int i_sim_ar,
		      double *i_pri_wgl_eff_mean,double *i_pri_wgl_eff_prec,
		      double *i_pri_wgl_prec_par,double *i_pri_wgl_ar,
		      int i_wgl_fixed_model, double *i_wgl_fixed_int, 
		      double *i_wgl_fixed_slp, double *i_wgl_fixed_tau, 
		      double *i_wgl_fixed_int_CC, double *i_wgl_fixed_slp_CC, 
		      double *i_wgl_fixed_tau_CC, int i_coastal_cod)
{
  int err=0;
  LW_struct *weight, *weight_mean, *weight_CC=NULL, *weight_CC_mean=NULL;
  
  (*GMRFLib_uniform_init)(i_seed);
  
  /* Allocating memory for weight structs */
  err = init_lin(i_D_wgl,i_sim_ar,i_pri_wgl_eff_mean,i_pri_wgl_eff_prec,
		 i_pri_wgl_prec_par,i_pri_wgl_ar,i_wgl_fixed_model,
		 &weight,&weight_mean);
  if(err)
    {
      write_warning("initialize_model2:Error calling init_lin\n");
      return(err);
    }

  weight->par->num_var = i_num_par[0];
  weight->gr_str->Qfunc = Qfunc_weight;

  weight->fixed_model = i_wgl_fixed_model;
  if(weight->fixed_model)
    {
      weight->fixed_int = i_wgl_fixed_int;
      weight->fixed_slp = i_wgl_fixed_slp;
      weight->fixed_tau = i_wgl_fixed_tau;
    }

  /* Allocating memory for extra weight structs if coastal cod */
  if(i_coastal_cod)
    {
      err = init_lin(i_D_wgl_CC,i_sim_ar,i_pri_wgl_eff_mean,i_pri_wgl_eff_prec,
		     i_pri_wgl_prec_par,i_pri_wgl_ar,i_wgl_fixed_model,
		     &weight_CC,&weight_CC_mean);
      if(err)
	{
	  write_warning("initialize_model2:Error calling init_lin\n");
	  return(err);
	}
      weight_CC->par->num_var = i_num_par[1];
      weight_CC->gr_str->Qfunc = Qfunc_weight;
      weight_CC->fixed_model = i_wgl_fixed_model;
      if(weight_CC->fixed_model)
	{
	  weight_CC->fixed_int = i_wgl_fixed_int_CC;
	  weight_CC->fixed_slp = i_wgl_fixed_slp_CC;
	  weight_CC->fixed_tau = i_wgl_fixed_tau_CC;
	}
    }
  
  /* Initial values for wgl model */
  err = init_lin_par(weight,i_D_wgl);
  if(err)
    {
      write_warning("initialize_model2:Error calling init_lin_par\n");
      return(err);
    }
  if(i_coastal_cod)
    {
      err = init_lin_par(weight_CC,i_D_wgl_CC);
      if(err)
	{
	  write_warning("initialize_model2:Error calling init_lin_par\n");
	  return(err);
	}
    }

  err = init_evaluate(i_D_wgl->glm->nHaul,1);
  if(err)
    {
      write_warning("initialize_model2:Error calling init_evaluate\n");
      return(err);
    }

  *o_weight = weight;
  *o_weight_mean = weight_mean;
  *o_weight_CC = weight_CC;
  *o_weight_CC_mean = weight_CC_mean;
  
  return(0);
}		/* end of initialize_model2 */




/*!
  \brief Reallocate memory allocated in initialize_model2
  \author Geir Storvik
*/
int re_initialize_model2(Data_lin *i_D_wgl, LW_struct *i_weight, LW_struct *i_weight_mean, 
			 Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, LW_struct *i_weight_CC_mean, 
			 int i_coastal_cod)
{
  int   err;

  err = re_init_lin(i_D_wgl,&i_weight,&i_weight_mean);
  if(err)
    {
      write_warning("re_initialize_model2:Error calling re_init_lin\n");
      return(err);
    }
  if(i_coastal_cod)
    {
      err = re_init_lin(i_D_wgl_CC,&i_weight_CC,&i_weight_CC_mean);
      if(err)
	{
	  write_warning("re_initialize_model2:Error calling re_init_lin\n");
	  return(err);
	}
    }


  err = re_init_evaluate(1);
  if(err)
    {
      write_warning("re_initialize_model2:Error calling re_init_evaluate\n");
      return(err);
    }

  return(0);
}		/* end of re_initialize_model2 */



/*!
  \brief Memory allocation and initial values for the sub-models
  \author Geir Storvik
*/
int initialize_predict(Data_age *i_D_age, Age_struct **o_age, 
		       Data_lin *i_D_lga, LW_struct **o_length, Data_g_a *i_D_g_a,
		       Data_lin *i_D_lga_CC, LW_struct **o_length_CC, Data_g_a *i_D_g_a_CC,
		       Data_lin *i_D_wgl, LW_struct **o_weight, Data_lin *i_D_wgl_CC, LW_struct **o_weight_CC,
		       Data_lin *i_D_hsz, LW_struct **o_hsz, 
		       int i_lga_cens_model, double *i_lga_cens_par, int i_coastal_cod)
{
  int err;
  Age_struct *age;
  LW_struct *length, *length_CC=NULL, *weight, *weight_CC=NULL, *hsz=NULL;

  /* Allocating memory for age structs */
  err = alloc_age(i_D_age,&age);
  if(err)
    {
      write_warning("initialize_predict:Error calling alloc_age\n");
      return(err);
    }
  age->delta_age = i_D_age->delta_age;

  /* Allocating memory for length structs */
  err = alloc_lin(i_D_lga,&length);
  if(err)
    {
      write_warning("initialize_predict:Error calling alloc_lin\n");
      return(err);
    }
 if(i_coastal_cod)
    {
      err = alloc_lin(i_D_lga_CC,&length_CC);
      if(err)
	{
	  write_warning("initialize_predict:Error calling alloc_lin\n");
	  return(err);
	}
    }

  /* Initial values if sampling discards */
  length->cens_model = i_lga_cens_model;
  if(length->cens_model)
    {
      length->cens_k = i_lga_cens_par[0];
      length->cens_m = i_lga_cens_par[1];
      length->cens_r = i_lga_cens_par[2];
      length->cens_Nlim = i_lga_cens_par[3];
    }

  /* Initial values if coastal cod */
  if(i_coastal_cod)
    {
      length_CC->cens_model = i_lga_cens_model;
      if(length_CC->cens_model)
	{
	  length_CC->cens_k = i_lga_cens_par[0];
	  length_CC->cens_m = i_lga_cens_par[1];
	  length_CC->cens_r = i_lga_cens_par[2];
	  length_CC->cens_Nlim = i_lga_cens_par[3];
	}
    }
 
  /* Allocating memory for weight structs */
  err = alloc_lin(i_D_wgl,&weight);
  if(err)
    {
      write_warning("initialize_predict:Error calling init_lin\n");
      return(err);
    }
  if(i_coastal_cod)
    {
      err = alloc_lin(i_D_wgl_CC,&weight_CC);
      if(err)
	{
	  write_warning("initialize_predict:Error calling init_lin\n");
	  return(err);
	}
    }

  /* Allocating memory for haulsize structs */
  if(i_D_age->glm->inc_hsz)
    {
      err = alloc_lin(i_D_hsz,&hsz);
      if(err)
	{
	  write_warning("initialize_predict:Error calling init_lin\n");
	  return(err);
	}
    }

  if(i_D_g_a->g_a_model>0)
    {
      err = sample_g_a_initialize(i_D_g_a->ncat,i_D_g_a->g_a_model);
      if(err)
	{
	  write_warning("initialize:Error calling sampling_g_a_initialize\n");
	  return(err);
	}
    }

  *o_age = age;
  *o_length = length;
  *o_length_CC = length_CC;
  *o_weight = weight;
  *o_weight_CC = weight_CC;
  *o_hsz = hsz;

  return(0);
}            /* End of initialize_predict */

    
/*!
  \brief Reallocate memory allocated in initialize_predict
  \author Geir Storvik
*/
int re_initialize_predict(Data_age *i_D_age, Age_struct *i_age, 
			  Data_lin *i_D_lga, LW_struct *i_length,
			  Data_lin *i_D_lga_CC, LW_struct *i_length_CC, Data_g_a *i_D_g_a, 
			  Data_lin *i_D_wgl, LW_struct *i_weight, 
			  Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, int i_coastal_cod,
			  Data_lin *i_D_hsz, LW_struct *i_hsz)
{
  int err;

  if(i_D_age->glm->inc_hsz)
    {
      err = re_alloc_lin(i_D_hsz,&i_hsz);
      if(err)
	{
	  write_warning("re_initialize_predict:Error calling re_init_lin\n");
	  return(err);
	}
    }

  if(i_D_g_a->g_a_model>0)
    {
      err = sample_g_a_re_initialize();
      if(err)
	{
	  write_warning("re_initialize_predict:Error calling sampling_g_a_re_initialize\n");
	  return(err);
	}
    }
  /* Allocating memory for age structs */
  err = re_alloc_age(i_D_age,&i_age);
  if(err)
    {
      write_warning("re_initialize_predict:Error calling re_alloc_age\n");
      return(err);
    }
  
  /* Allocating memory for length structs */
  err = re_alloc_lin(i_D_lga,&i_length);
  if(err)
    {
      write_warning("re_initialize_predict:Error calling re_alloc_lin\n");
      return(err);
    }
  if(i_coastal_cod)
    {
      err = re_alloc_lin(i_D_lga_CC,&i_length_CC);
      if(err)
	{
	  write_warning("re_initialize_predict:Error calling re_alloc_lin\n");
	  return(err);
	}
    }

  /* Allocating memory for weight structs */
  err = re_alloc_lin(i_D_wgl,&i_weight);
  if(err)
    {
      write_warning("re_initialize_predict:Error calling re_init_lin\n");
      return(err);
    }
  if(i_coastal_cod)
    {
      err = re_alloc_lin(i_D_wgl_CC,&i_weight_CC);
      if(err)
	{
	  write_warning("re_initialize_predict:Error calling re_init_lin\n");
	  return(err);
	}
    }

  return(0);
}            /* End of re_initialize_predict */


/*!
  \brief Initialize graph structures for age and lga model (for GMRFLib)
  \author Geir Storvik

  For both structures, Precision matrix calculated, Constraints are defined
  and Graph structure initiated by GMRLFib_init_problem is made
*/
int init_graph_model1(Data_age *i_D_age, Age_struct *i_age, Data_lin *i_D_lga, LW_struct *i_length,
		      Data_lin *i_D_lga_CC, LW_struct *i_length_CC, int i_coastal_cod, int i_constr)
{
  int err;

  #ifdef DEBUG_PROG
  printf("Initializing for age\n");
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Initializing for age\n");
  #endif

  // Remove haul effect before making graph
  i_D_age->glm->xcov[0]->n_cov--;
  /* Initialization for graph of fixed effects */
  err = make_graph_gauss(i_D_age->glm,i_age->gr_str_f,i_age->par,0,i_constr);
  if(err)
    {
      write_warning("init_graph_model1:Error calling make_graph_gauss\n");
      return(err);
    }
  // Add haul effect again
  i_D_age->glm->xcov[0]->n_cov++;
  i_age->gr_str_f->x_new = CALLOC(i_age->gr_str_f->graph->n,double);  // Free ok

  #ifdef DEBUG_PROG
  printf("Initializing for length\n");
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Initializing for length\n");
  #endif
  err = make_graph_gauss(i_D_lga->glm,i_length->gr_str,i_length->par,0,i_constr);
  if(err)
    {
      write_warning("init_graph_model1:Error calling make_graph_gauss\n");
      return(err);
    }
  if(i_coastal_cod)
    {
      err = make_graph_gauss(i_D_lga_CC->glm,i_length_CC->gr_str,i_length_CC->par,0,i_constr);
      if(err)
	{
	  write_warning("init_graph_model1:Error calling make_graph_gauss\n");
	  return(err);
	}
    }

  return(0);
}		/* end of init_graph_model1 */



/*!
  \brief  Reallocate space allocated by ::init_graph_model1
  \author Geir Storvik
*/
int re_init_graph_model1(Data_age *i_D_age, Age_struct *i_age, LW_struct *i_length, 
			 LW_struct *i_length_CC, int i_coastal_cod)
{
  int err;

  // Clean up
  // Remove haul effect before making graph
  i_D_age->glm->xcov[0]->n_cov--;
  /* Initialization for graph of fixed effects */
  err = re_make_graph_gauss(i_age->gr_str_f);
  if(err)
    {
      write_warning("re_init_graph_model1:Error calling re_make_graph_gauss\n");
      return(err);
    }
  // Add haul effect again
  i_D_age->glm->xcov[0]->n_cov++;
  FREE(i_age->gr_str_f->x_new);
  
  err = re_make_graph_gauss(i_length->gr_str);
  if(err)
    {
      write_warning("re_init_graph_model1:Error calling re_make_graph_gauss\n");
      return(err);
    }
  if(i_coastal_cod)
    {
      err = re_make_graph_gauss(i_length_CC->gr_str);
      if(err)
	{
	  write_warning("re_init_graph_model1:Error calling re_make_graph_gauss\n");
	  return(err);
	}
    }
  
  return(0);
}		/* end of re_init_graph_model1 */




/*!
  \author Geir Storvik
  \brief Allocate space and initialize age struct and struct for age-data

  Memory allocated in this routine is reallocated in ::re_init_age.
*/
int init_age(Data_age *i_D_age,int i_age_errors,double *i_A2A,int i_coastal_cod,int i_inc_hsz,
	     int i_sim_ar,double *i_pri_age_eff_mean,double *i_pri_age_eff_prec,
	     double *i_pri_age_prec_par,double *i_pri_age_ar,
             Age_struct **o_age,Age_struct **o_age_mean)
{
  int         a,a2,err,h,i,j,ind,N,ind1,ind2,ind3,ind4,ncov,k;
  double      sum;
  Age_struct *age;

  err = alloc_age(i_D_age,o_age);
  if(err)
    {
      write_warning("init_age:Error calling alloc_age\n");
      return(err);
    }
  age = *o_age;

  /* Include all effects except haul effect in linear graph */
  for(i=0;i<i_D_age->glm->nxcov;i++)  
    {
      for(j=0;j<i_D_age->glm->xcov[i]->n_cov;j++)
	age->gr_str_f->in_gr[i][j] = 1;
    }

  for(i=0;i<i_D_age->glm->nxcov;i++)  /* Only fixed effects are included in linear graph */
    {
      for(j=0;j<i_D_age->glm->xcov[i]->n_cov;j++)
	age->gr_str_r->in_gr[i][j] = 1-i_D_age->glm->xcov[i]->fix[j];
    }

  /* Initialize/allocate memory for haulsize */
  if(i_inc_hsz)
    {
      age->par->tau_hsz = 1;
      age->par->eff_hsz = CALLOC(i_D_age->glm->nHaul,double);
      for(i=0;i<i_D_age->glm->nHaul;i++)
	age->par->eff_hsz[i] = 1.0;
    }

  /* Specify that one observation is available on alpha */
  /* Initial values age */
  ind1 = 0;
  ind2 = 0;
  ind3 = 0;
  ind4 = 0;
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      if(i==0)
	ncov = i_D_age->glm->xcov[i]->n_cov-1;
      else
	ncov = i_D_age->glm->xcov[i]->n_cov;
      #ifdef LOG_FILE
      fprintf(g_caa_log,"init_age: ncov=%d\n",ncov);
      #endif
      for(j=0;j<ncov;j++)  //Not haul-effect here
	{
	  if(i_D_age->glm->xcov[i]->fix[j])
	    {
	      age->par->tau[i][j] = i_pri_age_eff_prec[ind1];
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"init_age: prior: fixed effect: tau[%d][%d]=%f,eff_prec[%d]=%f\n",i,j,age->par->tau[i][j],ind1,i_pri_age_eff_prec[ind1]);
	      #endif
	      ind1++;
	      for(k=0;k<i_D_age->glm->xcov[i]->n_fac[j];k++)
		{
		  for(a=0;a<i_D_age->glm->ncat;a++)
		    {
		      age->par->prior_mean[a][i][j][k] = i_pri_age_eff_mean[ind2];
		      ind2++;
		    }
		}
	    }
	  else
	    {
	      age->par->prior_prec[i][j][0] = i_pri_age_prec_par[ind3];
	      age->par->prior_prec[i][j][1] = i_pri_age_prec_par[ind3+1];
	      age->par->tau[i][j] = age->par->prior_prec[i][j][0]/age->par->prior_prec[i][j][1];
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"init_age: prior: random effect: cov %d: tau[%d][%d]=%f, prior_prec1=%f, prior_prec2=%f, prec_par[%d]=%f,prec_par[%d]=%f\n",
		      j,i,j,age->par->tau[i][j],age->par->prior_prec[i][j][0],age->par->prior_prec[i][j][1],ind3,i_pri_age_prec_par[ind3],ind3+1,i_pri_age_prec_par[ind3+1]);
	      #endif
	      ind3 += 2;
	      for(k=0;k<i_D_age->glm->xcov[i]->n_fac[j];k++)
		for(a=0;a<i_D_age->glm->ncat;a++)
		  age->par->prior_mean[a][i][j][k] = G_ZERO;
	    }
	}
      if(i_D_age->glm->xcov[i]->ispat >= 0)
	{
	  if(i_sim_ar)
	    {
	      age->par->sim_ar = i_sim_ar;
	      age->par->prior_ar[i][0] = i_pri_age_ar[ind4];
	      age->par->prior_ar[i][1] = i_pri_age_ar[ind4+1];
	      age->par->ar[i] = age->par->prior_ar[i][0]/(age->par->prior_ar[i][0]+age->par->prior_ar[i][1]);
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"init_age: prior: ar[%d]=%f, pri_age_ar[%d]=%f,pri_age_ar[%d]=%f\n",i,age->par->ar[i],ind4,i_pri_age_ar[ind4],ind4+1,i_pri_age_ar[ind4+1]);
              #endif
	      ind4 +=2;
	    }
	  else
	    age->par->ar[i] = G_ZERO;
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"init_age: prior: ar[%d]=%f\n",i,age->par->ar[i]);
	  #endif
	}      
    }
  age->par->prior_prec_obs[0] = i_pri_age_prec_par[ind3];
  age->par->prior_prec_obs[1] = i_pri_age_prec_par[ind3+1];
  age->par->tau_obs = age->par->prior_prec_obs[0]/age->par->prior_prec_obs[1];
  #ifdef LOG_FILE
  fprintf(g_caa_log,"init_age: prior: tau_obs=%f,age_prec_par[%d]=%f,age_prec_par[%d]=%f\n\n",age->par->tau_obs,ind3,i_pri_age_prec_par[ind3],ind3+1,i_pri_age_prec_par[ind3+1]);
  #endif
  age->par->loglik = G_ZERO;

  i_D_age->glm->suff = Mmatrix_3d(0,i_D_age->glm->nHaul-1,   // Free ok
                                  0,i_D_age->glm->nxcov-1,
                                  0,i_D_age->glm->nxcov-1,sizeof(double),1);

  i_D_age->glm->beta_hat = Mmatrix_3d(0,i_D_age->glm->nHaul-1,0,i_D_age->glm->ncat-1,0,i_D_age->glm->nxcov-1,
				      sizeof(double),1); // Free ok

  i_D_age->glm->ssq = CALLOC(i_D_age->glm->nHaul,double);  // Free ok

  /* Specify that one observation is available on alpha */
  for(h=0;h<i_D_age->glm->nHaul;h++)  
      i_D_age->glm->suff[h][0][0] = G_ONE;


  for(h=0;h<i_D_age->glm->nHaul;h++)  
      i_D_age->glm->ssq[h] = G_ZERO;

  err = alloc_age(i_D_age,o_age_mean);
  if(err)
    {
      write_warning("init_age:Error calling alloc_age\n");
      return(err);
    }

  int ncat;
  i_D_age->type_age = CALLOC(i_D_age->glm->nHaul,int);   // Free ok
  age->age_errors = i_age_errors;
  #ifdef LOG_FILE
  fprintf(g_caa_log,"init_age:age_errors=%d\n\n",age->age_errors);
  #endif
  if(age->age_errors)
    {
      age->A2A = Mmatrix_2d(0,i_D_age->glm->ncat-1,      
                            0,i_D_age->glm->ncat-1,sizeof(double),1);// Free ok
      ind = 0;
      if(i_coastal_cod)//Coastal cod
	ncat = i_D_age->glm->ncat/2;
      else
	ncat = i_D_age->glm->ncat;
      for(a=0;a<ncat;a++)
	{
	  sum = G_ZERO;
	  for(a2=0;a2<ncat;a2++)
	    {
	      age->A2A[a2][a] = i_A2A[ind];
	      sum += i_A2A[ind];
	      if(i_coastal_cod)//Coastal cod
		age->A2A[ncat+a2][ncat+a] = i_A2A[ind];
	      ind++;
	    }
	  if(sum>1.000001 || sum < 0.99999)
	    {
	      write_warning("init_age:Error reading age error matrix\n");
	      printf("Age error matrix: columns must sum to 1, sum=%f\n",sum);
	      return(1);
	    }
	}
      #ifdef LOG_FILE
      fprintf(g_caa_log,"init_age: A2A=\n");
      for(a=0;a<ncat;a++)
	{
	  for(a2=0;a2<ncat;a2++)
	    fprintf(g_caa_log,"%f ",age->A2A[a][a2]);
	  fprintf(g_caa_log,"\n");
	}
      #endif
      /* Find "neighbor" ages */
      age->A_Nneigh = CALLOC(i_D_age->glm->ncat,int);   // Free ok
      age->A_neigh = CALLOC(i_D_age->glm->ncat,int *);  // Free ok
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
          /* First count number of neighbors */
	  N = 0;
	  for(a2=0;a2<i_D_age->glm->ncat;a2++)
	    {
	      if(age->A2A[a][a2]> 0.000000001)
		N++;
	    }
	  age->A_Nneigh[a] = N;
	  age->A_neigh[a] = CALLOC(N,int);              // Free ok
          N = 0;
	  for(a2=0;a2<i_D_age->glm->ncat;a2++)
	    {
	      if(age->A2A[a][a2]> 0.000000001)
		{
		  age->A_neigh[a][N] = a2;
		  N++;
		}
	    }
	}

      for(h=0;h<i_D_age->glm->nHaul;h++)
	i_D_age->type_age[h] = 1;
    }
  else
    {
      for(h=0;h<i_D_age->glm->nHaul;h++)
	i_D_age->type_age[h] = 0;
    }


  age = *o_age_mean;
  age->par->loglik = G_ZERO;



  return(0);
}		/* end of init_age */


/*!
  \author Geir Storvik
  \brief Reallocate space allocated by ::init_age
*/
int re_init_age(Data_age *i_D_age,int i_age_errors,double *i_A2A,int i_inc_hsz,
             Age_struct **o_age,Age_struct **o_age_mean)
{
  int         a,err;
  Age_struct *age;

  if(i_inc_hsz)
    {
      age = *o_age;
      FREE(age->par->eff_hsz);
    }
  err = re_alloc_age(i_D_age,o_age);
  if(err)
    {
      write_warning("re_init_age:Error calling re_alloc_age\n");
      return(err);
    }
  age = *o_age;

  Fmatrix_3d(&i_D_age->glm->suff[0][0][0],&i_D_age->glm->suff[0][0],
             &i_D_age->glm->suff[0]);

  Fmatrix_3d(&i_D_age->glm->beta_hat[0][0][0],&i_D_age->glm->beta_hat[0][0],
             &i_D_age->glm->beta_hat[0]);
  FREE(i_D_age->glm->ssq);

  err = re_alloc_age(i_D_age,o_age_mean);
  if(err)
    {
      write_warning("re_init_age:Error calling re_alloc_age\n");
      return(err);
    }

  FREE(i_D_age->type_age);
  if(i_age_errors)
    {
      Fmatrix_2d(&age->A2A[0][0],&age->A2A[0]);
      /* Find "neighbor" ages */
      FREE(age->A_Nneigh);
      for(a=0;a<i_D_age->glm->ncat;a++)
	  FREE(age->A_neigh[a]);
      FREE(age->A_neigh);

    }

  return(0);
}		/* end of re_init_age */


/*!
  \author Geir Storvik
  \brief Allocating space for age-struct
*/
int alloc_age(Data_age *i_D_age,Age_struct **o_age)
{
  int         err,i,j;
  Age_struct *age;

  age = CALLOC(1,Age_struct);                  // Free ok
  age->gr_str_f = CALLOC(1,Graph_str);         // Free ok
  age->gr_str_f->in_gr = CALLOC(i_D_age->glm->nxcov,int *);  //Free ok

  for(i=0;i<i_D_age->glm->nxcov;i++)  
      age->gr_str_f->in_gr[i] = 
	CALLOC(i_D_age->glm->xcov[i]->n_cov,int);  // Free ok

  age->gr_str_r = CALLOC(1,Graph2_str);            // Free ok
  age->gr_str_r->in_gr = CALLOC(i_D_age->glm->nxcov,int *); // Free ok
  for(i=0;i<i_D_age->glm->nxcov;i++)  /* Only fixed effects are included in linear graph */
    {
      age->gr_str_r->in_gr[i] = 
	CALLOC(i_D_age->glm->xcov[i]->n_cov,int);   // Free ok
      for(j=0;j<i_D_age->glm->xcov[i]->n_cov;j++)
	age->gr_str_r->in_gr[i][j] = 1-i_D_age->glm->xcov[i]->fix[j];
    }
  age->alpha = Mmatrix_2d(0,i_D_age->glm->nHaul-1,  // Free ok
		   	    0,i_D_age->glm->ncat-1,sizeof(double),1);

  age->mu = CALLOC(i_D_age->glm->nHaul,double);       // Free ok
  err = alloc_Eff_str(i_D_age->glm->ncat,i_D_age->glm->nxcov,
                        i_D_age->glm->xcov,&age->par);
  if(err)
    {
      write_warning("alloc_age:Error calling alloc_Eff_str\n");
      return(err);
    }

  *o_age = age;
  return(0);
}		/* end of alloc_age */


/*!
  \author Geir Storvik
  \brief Re-allocating space allocated by alloc_age
*/
int re_alloc_age(Data_age *i_D_age,Age_struct **o_age)
{
  int         err,i;
  Age_struct *age;

  age = *o_age;

  for(i=0;i<i_D_age->glm->nxcov;i++)  
      FREE(age->gr_str_f->in_gr[i]);
  FREE(age->gr_str_f->in_gr);
  FREE(age->gr_str_f);


  for(i=0;i<i_D_age->glm->nxcov;i++)  
      FREE(age->gr_str_r->in_gr[i]);
  FREE(age->gr_str_r->in_gr);
  FREE(age->gr_str_r);

  Fmatrix_2d(&age->alpha[0][0],&age->alpha[0]);

  FREE(age->mu);
  err = re_alloc_Eff_str(i_D_age->glm->ncat,i_D_age->glm->nxcov,
                        i_D_age->glm->xcov,&age->par);
  if(err)
    {
      write_warning("re_alloc_age:Error calling re_alloc_Eff_str\n");
      return(err);
    }

  FREE(age);

  return(0);
}		/* end of re_alloc_age */


/*!
  \author Geir Storvik
  \brief Allocate space and initialize linear struct and data for a linear model
*/
int init_lin(Data_lin *i_D_lin,int i_sim_ar,double *i_pri_lin_eff_mean,double *i_pri_lin_eff_prec,
	     double *i_pri_lin_prec_par,double *i_pri_lin_ar,int i_fixed_model,
	     LW_struct **o_lin,LW_struct **o_lin_mean)
{
  int        err,a,i,j,ind1,ind2,ind3,ind4,k;
  LW_struct *lin;

  err = alloc_lin(i_D_lin,o_lin);
  if(err)
    {
      write_warning("init_lin:Error calling alloc_lin\n");
      return(err);
    }
  /* Initial values */
  ind1 = 0;
  ind2 = 0;
  ind3 = 0;
  ind4 = 0;
  lin = *o_lin;
  for(i=0;i<i_D_lin->glm->nxcov;i++)
    {
      #ifdef LOG_FILE
      fprintf(g_caa_log,"init_lin: ncov=%d\n",i_D_lin->glm->xcov[i]->n_cov);
      #endif
      for(j=0;j<i_D_lin->glm->xcov[i]->n_cov;j++)
	{
	  if(i_D_lin->glm->xcov[i]->fix[j])
	    {
	      lin->par->tau[i][j] = i_pri_lin_eff_prec[ind1];
	      #ifdef LOG_FILE
	      fprintf(g_caa_log,"init_lin: prior: fixed effect: tau[%d][%d]=%f\n",i,j,lin->par->tau[i][j]);
	      #endif
	      ind1++;
	      for(k=0;k<i_D_lin->glm->xcov[i]->n_fac[j];k++)
		{
		  for(a=0;a<i_D_lin->glm->ncat;a++)
		    {
		      lin->par->prior_mean[a][i][j][k] = i_pri_lin_eff_mean[ind2];
                      #ifdef LOG_FILE
		      fprintf(g_caa_log,"prior:fixed effect: factor %d: prior_mean[%d][%d][%d][%d]=%f, eff_mean[%d]=%f\n",
			      k,a,i,j,k,lin->par->prior_mean[a][i][j][k],ind2,i_pri_lin_eff_mean[ind2]);
     	              #endif
		      //fprintf(stderr,"prior:fixed effect: factor%d ind=%d prior_mean=%f\n",k,ind2,i_pri_lin_eff_mean[ind2]);
		      ind2++;
		    }
		}
	    }
	  else
	    {
	      lin->par->prior_prec[i][j][0] = i_pri_lin_prec_par[ind3];
	      lin->par->prior_prec[i][j][1] = i_pri_lin_prec_par[ind3+1];
	      if(i_fixed_model==0)
		{
		  lin->par->tau[i][j] = lin->par->prior_prec[i][j][0]/lin->par->prior_prec[i][j][1];
		}
	      else
		{
		  lin->par->tau[i][j] = 1000000.0;
		}
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"init_lin: prior: random effect: cov%d: prior_prec1=%f, prior_prec2=%f, prior_tau=%f\n",
		      j,lin->par->prior_prec[i][j][0],lin->par->prior_prec[i][j][1],lin->par->tau[i][j]);
	      #endif
	      ind3 += 2;
	      for(k=0;k<i_D_lin->glm->xcov[i]->n_fac[j];k++)
		for(a=0;a<i_D_lin->glm->ncat;a++)
		  lin->par->prior_mean[a][i][j][k] = G_ZERO;
	    }
	}
      if(i_D_lin->glm->xcov[i]->ispat >= 0)
	{
	  if(i_sim_ar)
	    {
	      lin->par->sim_ar = i_sim_ar;
	      lin->par->prior_ar[i][0] = i_pri_lin_ar[ind4];
	      lin->par->prior_ar[i][1] = i_pri_lin_ar[ind4+1];
	      //printf("ispat=%f\n",i_D_lin->glm->xcov[i]->ispat);
	      //printf("i=%d,ar[0]=%f,ar[1]=%f\n",i,i_pri_lin_ar[ind4],i_pri_lin_ar[ind4+1]);
	      ind4 +=2;
	      lin->par->ar[i] = lin->par->prior_ar[i][0]/
		(lin->par->prior_ar[i][0]+lin->par->prior_ar[i][1]);
	    }
	  else
	    lin->par->ar[i] = G_ZERO;
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"init_lin: prior: ar[%d]=%f\n",i,lin->par->ar[i]);
	  #endif
	}
    }
  lin->par->prior_prec_obs[0] = i_pri_lin_prec_par[ind3];
  lin->par->prior_prec_obs[1] = i_pri_lin_prec_par[ind3+1];
  lin->par->tau_obs = lin->par->prior_prec_obs[0]/lin->par->prior_prec_obs[1];
  lin->par->loglik   = G_ZERO;
  #ifdef LOG_FILE
  fprintf(g_caa_log,"init_lin: prior: tau_obs=%f\n",lin->par->tau_obs);
  #endif

  err = alloc_lin(i_D_lin,o_lin_mean);
  if(err)
    {
      write_warning("init_lin:Error calling alloc_lin\n");
      return(err);
    }
  lin = *o_lin_mean;

  lin->par->loglik   = G_ZERO;


  return(0);
}		/* end of init_lin */


/*!
  \author Geir Storvik
  \brief Reallocate space allocated by init_lin
*/
int re_init_lin(Data_lin *i_D_lin,
             LW_struct **o_lin,LW_struct **o_lin_mean)
{
  int        err;

  err = re_alloc_lin(i_D_lin,o_lin);
  if(err)
    {
      write_warning("re_init_lin:Error calling re_alloc_lin\n");
      return(err);
    }

  err = re_alloc_lin(i_D_lin,o_lin_mean);
  if(err)
    {
      write_warning("re_init_lin:Error calling alloc_lin\n");
      return(err);
    }

  return(0);
}		/* end of re_init_lin */


/*!
  \author Geir Storvik
  \brief Allocating space for linear struct 
*/
int alloc_lin(Data_lin *i_D_lin,LW_struct **o_lin)
{
  int        err,i,j;
  LW_struct *lin;

  lin = CALLOC(1,LW_struct);             // Free ok
  lin->gr_str = CALLOC(1,Graph_str);     // Free ok
  err = alloc_Eff_str(1,i_D_lin->glm->nxcov,i_D_lin->glm->xcov,&lin->par);
  if(err)
    {
      write_warning("alloc_lin:Error calling alloc_Eff_str\n");
      return(err);
    }
  
  lin->gr_str->in_gr = CALLOC(i_D_lin->glm->nxcov,int *);    // Free ok
  for(i=0;i<i_D_lin->glm->nxcov;i++)  
    {
      lin->gr_str->in_gr[i] = CALLOC(i_D_lin->glm->xcov[i]->n_cov,int);   // Free ok
      for(j=0;j<i_D_lin->glm->xcov[i]->n_cov;j++)
	lin->gr_str->in_gr[i][j] = 1;/* All effects are included in linear graph */
    }

  *o_lin = lin;

  return(0);     /* end of alloc_lin */
}



/*!
  \author Geir Storvik
  \brief Allocating space for linear struct 
*/
int re_alloc_lin(Data_lin *i_D_lin,LW_struct **o_lin)
{
  int        err,i;
  LW_struct *lin;

  lin = *o_lin;

  err = re_alloc_Eff_str(1,i_D_lin->glm->nxcov,i_D_lin->glm->xcov,&lin->par);
  if(err)
    {
      write_warning("re_alloc_lin:Error calling alloc_Eff_str\n");
      return(err);
    }
  
  for(i=0;i<i_D_lin->glm->nxcov;i++)  
     FREE(lin->gr_str->in_gr[i]);
  FREE(lin->gr_str->in_gr);

  FREE(lin->gr_str);
  FREE(lin);

  return(0);     /* end of re_alloc_lin */
}


/*!
  \author Geir Storvik
  \brief Initialize parameters in lga and wgl model
*/
int init_lin_par(LW_struct *i_lin,Data_lin *i_D_lin)
{
  int     h;
  double  Int,Slp,W;

  Int = G_ZERO;
  Slp = G_ZERO;
  W = G_ZERO;

  if(i_lin->fixed_model == 0)
    {
      for(h=0;h<i_D_lin->glm->nHaul;h++)
	{
	  Int += i_D_lin->glm->beta_hat[h][0][0]*i_D_lin->glm->suff[h][0][0];
	  Slp += i_D_lin->glm->beta_hat[h][0][1]*i_D_lin->glm->suff[h][0][0];
	  W += i_D_lin->glm->suff[h][0][0];
	}
      Int /= W;
      Slp /= W;
    }
  else
    {
      Int = i_lin->fixed_int[0];
      Slp = i_lin->fixed_slp[0];
      i_lin->par->tau_obs = i_lin->fixed_tau[0];
    }

  /* Assume first covariate is constant term */
  if(i_D_lin->glm->xcov[0]->n_fac[0]==1)
      i_lin->par->eff[0][0][0][0] = Int;
  else
    {
      write_warning("init_lin_par:First covariate should be constant term\n");
      return(1);
    }

  if(i_D_lin->glm->xcov[1]->n_fac[0]==1)
      i_lin->par->eff[0][1][0][0] = Slp;
  else
    {
      write_warning("init_lin_par:First covariate should be constant term\n");
      return(1);
    }
  #ifdef DEBUG_PROG
  printf("Init_lin_par:\n");
  printf("Int = %f\n",i_lin->par->eff[0][0][0][0]);
  printf("Slp = %f\n",i_lin->par->eff[0][1][0][0]);
  #endif


  return(0);    /* End of init_lin_par */
}




/*!
  \author Geir Storvik
  \brief Pick out an element of the precision matrix for the age model

  For use when sampling using the GMRFLib library.
  Requires that make_Q_age_fix is called first.
*/
double Qfunc_age_fix(int node, int nnode,char *arg)
{
  double **Q;
  Q = (double **) arg;

  return(Q[node][nnode]);
}		/* end of Qfunc_age_fix */




/*!
  \author Geir Storvik
  \brief Pick out an element of the precision matrix for the random effects of the age model

  For use when sampling using the GMRFLib library.
  Requires that make_Q_age_ran is called first. From an earlier version where the
  GMRFLib library was used for simulating haul effects in the age model. Now this is done
  in new routines.
*/
double Qfunc_age_ran(int node, int nnode,char *arg)
{
  fprintf(stderr,"Must change Qfunc_age_ran\n");
  return(1);
    //  return(i_age->par->tau_obs*(node==nnode));
}		/* end of Qfunc_age_ran */



/*!
  \author Geir Storvik
  \brief Pick out an element of the precision matrix for the lga model

  For use when sampling using the GMRFLib library.
  Requires that make_Q_length is called first.
*/
double Qfunc_length(int node, int nnode,char *arg)
{
  double **Q;
  Q = (double **) arg;

  return(Q[node][nnode]);
}		/* end of Qfunc_length */



/*!
  \author Geir Storvik
  \brief Pick out an element of the precision matrix for the wgl model

  For use when sampling using the GMRFLib library.
  Requires that ::make_Q_weight is called first.
*/
double Qfunc_weight(int node, int nnode,char *arg)
{
  double **Q;
  Q = (double **) arg;

  return(Q[node][nnode]);
}		/* end of Qfunc_weight */


