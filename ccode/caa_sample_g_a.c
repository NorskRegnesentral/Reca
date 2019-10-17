/*!
  \file caa_sample_g_a.c
  \author Geir Storvik
  \brief Routines for sampling nonlinear part of lga model

*/
#include "caa.h"
#include "caa_sample_g_a.h"
#include "caa_read_write.h"
#include "caa_chol.h"
#include "caa_routines.h"
#include "caa_util.h"


static int calc_g_a_S_R(int i_ncat,double *i_a_vec,double *i_par,double *o_g);
static int calc_g_a_polyn3(int i_ncat,double *i_a_vec,double *i_par,double *o_g);
static int (*s_calc_g_a)(int i_ncat,double *i_a_vec,double *i_par,double *o_g);

static double   *s_g_a;

#ifdef LOG_FILE
extern FILE     *g_caa_log;
#endif


/*!
  \author Geir Storvik
  \brief Allocates space for routines used to sample non-linear part of lga-model

  This routine must be called before using routines in the caa_sample_g_a.c file.
*/
int sample_g_a_initialize(int i_ncat,int i_g_a_model)
{
  if(i_g_a_model==1)
    {
      s_calc_g_a = calc_g_a_S_R;
    }
  else if(i_g_a_model==2)
    {
      s_calc_g_a = calc_g_a_polyn3;
    }
  s_g_a = CALLOC(i_ncat,double);    // Free ok
  if(!s_g_a)
    {
      write_warning("sample_g_a_initialize:Error allocating s_g_a\n");
      return(1);
    }


  return(0);
}		/* end of sample_g_a_initialize */




/*!
  \author Geir Storvik
  \brief Reallocate space allocated by sample_g_a_initialize
*/
int sample_g_a_re_initialize()
{
  FREE(s_g_a);

  return(0);
}		/* end of sample_g_a_re_initialize */




/*L:suff_g_a*
________________________________________________________________

		suff_g_a
________________________________________________________________

Name:		suff_g_a
Syntax:		
Description:    Calculates sufficient statistics for sampling the function g(a) 
                describing the non-linear
                link between age and length
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________ 
*/
int suff_g_a(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
	     Data_g_a *i_D_g_a,int i_start_h,int i_nHaul,double **o_suff)
{
  int         a,h,i;
  double      sum_l,sum_l2,N_h_a,weight;
  double     *beta;
  Data_glm   *glm;

  glm = i_D_lga->glm;
  beta = CALLOC(glm->nxcov,double); //Free ok


  if(glm->nxcov==2)
    {
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  o_suff[0][a] = G_ZERO;
	  o_suff[1][a] = G_ZERO;
	  o_suff[2][a] = G_ZERO;
	  for(h=i_start_h;h<i_nHaul;h++)
	    {
	      sum_l = i_D_lga->sum_by_cat_fix[h][a];
	      sum_l2 = i_D_lga->sqsum_by_cat_fix[h][a];
	      N_h_a = (double) i_D_lga->Ages_fix[h][a];
	      for(i=0;i<glm->nxcov;i++)
		beta[i] = calc_eff(glm->xcov[i],i_length->par->eff[0][i],h);
	      o_suff[0][a] += sum_l2 - 2*beta[0]*sum_l + N_h_a*beta[0]*beta[0];
	      o_suff[1][a] += -2*(sum_l*beta[1] - N_h_a*beta[0]*beta[1]);
	      o_suff[2][a] += beta[1]*beta[1]*N_h_a;
	    }
	}
    }
  else if(glm->nxcov==3)
    {
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  o_suff[0][a] = G_ZERO;
	  o_suff[1][a] = G_ZERO;
	  o_suff[2][a] = G_ZERO;
	  for(h=i_start_h;h<i_nHaul;h++)
	    {
	      sum_l = i_D_lga->sum_by_cat_fix[h][a];
	      sum_l2 = i_D_lga->sqsum_by_cat_fix[h][a];
	      weight = i_D_lga->haulweight[h];
	      N_h_a = (double) i_D_lga->Ages_fix[h][a];
	      for(i=0;i<glm->nxcov;i++)
		beta[i] = calc_eff(glm->xcov[i],i_length->par->eff[0][i],h);
	      o_suff[0][1] += sum_l2 - 2*beta[0]*sum_l + N_h_a*beta[0]*beta[0] 
		              - 2*beta[2]*weight*sum_l + 2*beta[0]*beta[2]*weight
                              + beta[2]*weight*beta[2]*weight;
	      o_suff[1][a] += -2*(sum_l*beta[1] - N_h_a*beta[0]*beta[1]
                              -N_h_a*beta[2]*weight*beta[1]);
	      o_suff[2][a] += beta[1]*beta[1]*N_h_a;
	    }
	}
    }

  FREE(beta);

  return(0);
}               /* end of suff_g_a */



/*!
  \brief Calculates the log-likelihood (apart from precision, and prior)
*/
int calc_lik_g_a(double *o_log_lik,double **i_suff,int i_ncat)
{
  int a;
  double log_lik;
  
  log_lik = G_ZERO;
  for(a=0;a<i_ncat;a++)
    {
      log_lik += i_suff[0][a] + i_suff[1][a]*s_g_a[a] + i_suff[2][a]*s_g_a[a]*s_g_a[a];
    }

  *o_log_lik = log_lik;
  
  return(0);
}               /* end of calc_lik_g_a */



/*!
  \brief Samples the parameters in the function g(a) describing the non-linear
  link between age and length, assuming g(a)=log(1-exp(-K(a-a_0)))
  
  The algorithm is a systematic scan Metropolis-Hastings algorithm
  where for c and theta a scale-proposal
  
  Because ths sampling takes very little time, num=1000 iteration are
  performe in order to get reasonable estimates from the 
  conditional posterior. This could probably be made more
  efficient. 

  \author Geir Storvik
*/
int sample_g_a(LW_struct *i_length,Data_age *i_D_age,Data_lin *i_D_lga,
               Data_g_a *i_D_g_a,double **i_suff,int i_it,int i_nHaul)
{
  int       a,i;
  int       err,num=10;
  int      *acc_par;
  double    fac,u,d,tau,l_new,l_cur;
  double   *par_new;

  par_new = CALLOC(i_D_g_a->g_a_npar,double);// Free ok
  acc_par = CALLOC(3,int);

  tau = i_length->par->tau_obs;

  /* Calculate current log-likelihood */
  err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,s_g_a);
  if(err)
    {
      write_warning("sample_g_a:Error calling s_calc_g_a\n");
      return(err);
    }
  err = calc_lik_g_a(&l_cur,i_suff,i_D_g_a->ncat);
  if(err)
    {
      write_warning("sample_g_a:Error calling calc_lik_g_a\n");
      return(err);
    }

  #ifdef LOG_FILE
  fprintf(g_caa_log,"sample_g_a: Starting sample_g_a\n");
  #endif

  /* Run Metropolis-Hastings for parameters */
  fac = 1.01;
  for(i=0;i<num;i++)
    {
      // Sample parameters i g-function
      #ifdef DEBUG_G_A
      printf("old:c=%lg,theta=%lg,gamma=%lg,l_cur=%lf\n",
	     i_D_g_a->g_a_par[0],i_D_g_a->g_a_par[1],i_D_g_a->g_a_par[2],l_cur);
      #endif
      par_new[0] = i_D_g_a->g_a_par[0];
      par_new[1] = i_D_g_a->g_a_par[1];
      par_new[2] = i_D_g_a->g_a_par[2];
      if(i_D_g_a->sample_c)
	{
	  par_new[0] = scale_proposal(i_D_g_a->g_a_par[0],fac,NULL);
	  if(par_new[0]>0 && par_new[0]<5)
	    {
	      err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,par_new,s_g_a);
	      err = calc_lik_g_a(&l_new,i_suff,i_D_g_a->ncat);
	      u = genunf(G_ZERO,G_ONE);
	      d = -G_HALF*tau*(l_new-l_cur);
	      //add prior
	      d += (P_GAMMA_A-G_ONE)*(log(par_new[0])-log(i_D_g_a->g_a_par[0]))-P_GAMMA_B*(par_new[0]-i_D_g_a->g_a_par[0]);
	      if(d > -1.0e32 && d < 1.0e32 && log(u) < d)
		{
		  i_D_g_a->g_a_par[0] = par_new[0];
		  acc_par[0]++;
		  l_cur = l_new;
		}
	    }
	}
      if(i_D_g_a->sample_theta)
	{
	  par_new[1] = scale_proposal(i_D_g_a->g_a_par[1],fac,NULL);
	  if(par_new[1]>0 && par_new[1]<5)
	    {
	      err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,par_new,s_g_a);
	      err = calc_lik_g_a(&l_new,i_suff,i_D_g_a->ncat);
	      u = genunf(G_ZERO,G_ONE);
	      d = -G_HALF*tau*(l_new-l_cur);
	      d += (P_GAMMA_A-G_ONE)*(log(par_new[1])-log(i_D_g_a->g_a_par[1]))-P_GAMMA_B*(par_new[1]-i_D_g_a->g_a_par[1]);
	      if(d > -1.0e32 && d < 1.0e32 && log(u) < d)
		{
		  i_D_g_a->g_a_par[1] = par_new[1];
		  acc_par[1]++;
		  l_cur = l_new;
		}
	    }
	}
      if(i_D_g_a->sample_gamma)
	{
	  par_new[2] = scale_proposal(i_D_g_a->g_a_par[2],fac,NULL);
	  if(par_new[2]>0 && par_new[2]<5)
	    {
	      err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,par_new,s_g_a);
	      err = calc_lik_g_a(&l_new,i_suff,i_D_g_a->ncat);
	      u = genunf(G_ZERO,G_ONE);
	      d = -G_HALF*tau*(l_new-l_cur);
	      d += (P_GAMMA_A-G_ONE)*(log(par_new[2])-log(i_D_g_a->g_a_par[2]))-P_GAMMA_B*(par_new[2]-i_D_g_a->g_a_par[2]);
	      if(d > -1.0e32 && d < 1.0e32 && log(u) < d)
		{
		  i_D_g_a->g_a_par[2] = par_new[2];
		  acc_par[2]++;
		  l_cur = l_new;
		}
	    }
	}
      #ifdef DEBUG_G_A
      printf("new:c=%lg,theta=%lg,gamma=%lg,l_new=%lf\n",
	     par_new[0],par_new[1],par_new[2],l_new);
      printf("g_a=");
      for(a=0;a<i_D_g_a->ncat;a++)
	printf("%lf, ",s_g_a[a]);
      printf("\n");
      #endif

    }

  /* Calculate the g-function again for current parameters */
  err = s_calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,s_g_a);
  for(a=0;a<i_D_g_a->ncat;a++)
    i_D_g_a->g_a[a] = s_g_a[a];


  #ifdef LOG_FILE
  fprintf(g_caa_log,"sample_g_a: c=%lg,theta=%lg,gamma=%lg\n",i_D_g_a->g_a_par[0],i_D_g_a->g_a_par[1],i_D_g_a->g_a_par[2]);
  #endif

  #ifdef DEBUG_G_A
  printf("c=%lg,theta=%lg,gamma=%lg, acc=(%lf %lf %lf)\n",
         i_D_g_a->g_a_par[0],i_D_g_a->g_a_par[1],
         i_D_g_a->g_a_par[2],(double) acc_par[0]/(double) num,
         (double) acc_par[1]/(double) num,(double) acc_par[2]/(double) num);
  printf("g_a=");
  for(a=0;a<i_D_g_a->ncat;a++)
    printf("%lf, ",s_g_a[a]);
  printf("\n");
  #endif

  FREE(par_new);
  FREE(acc_par);


  return(0);
}               /* end of sample_g_a */




/*L:calc_g_a*
________________________________________________________________

		calc_g_a
________________________________________________________________

Name:		calc_g_a
Syntax:		
Description:    Calculates g_a function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
int calc_g_a(int i_ncat,double *i_a_vec,double *i_par,double *o_g)
{
  return(s_calc_g_a(i_ncat,i_a_vec,i_par,o_g));

}		/* end of calc_g_a */



/*L:calc_g_a_S_R*
________________________________________________________________

		calc_g_a_S_R
________________________________________________________________

Name:		calc_g_a_S_R
Syntax:		
Description:    Calculates S_R function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int calc_g_a_S_R(int i_ncat,double *i_a_vec,double *i_par,double *o_g_a)
{
  int    i;
  double a,amin,amax,g,g_min,g_max,r;
  double theta,gamma,c;

  amin = i_a_vec[0];
  amax = i_a_vec[i_ncat-1];


  c = i_par[0];
  theta = i_par[1];
  gamma = i_par[2];

  if(c!=1)
    fprintf(stderr,"NOTE! Change calc_g_a_S_R when c!=1\n");
  g_min = log(1-theta*exp(-gamma * amin));
  g_max = log(1-theta*exp(-gamma * amax));
  r = g_max-g_min;
  for(i=0;i<i_ncat;i++)
    {
      a = i_a_vec[i];
      g = log(1-theta*exp(-gamma * a));
      //o_g_a[i] = log(1-theta*exp(-gamma * a));
      o_g_a[i] = (g-g_min)/r;
    }

  return(0);
}		/* end of calc_g_a_S_R */

int calc_g_a_S_R2(int i_ncat,double *i_a_vec,double *i_par,double *o_g_a)
{
  int    i;
  double a,amin,amax,g,g_min,g_max,r;
  double theta,gamma,c;

  amin = i_a_vec[0];
  amax = i_a_vec[i_ncat-1];


  c = i_par[0];
  theta = i_par[1];
  gamma = i_par[2];

  if(c!=1)
    fprintf(stderr,"NOTE! Change calc_g_a_S_R when c!=1\n");
  g_min = log(1-theta*exp(-gamma * amin));
  g_max = log(1-theta*exp(-gamma * amax));
  r = g_max-g_min;
  for(i=0;i<i_ncat;i++)
    {
      a = i_a_vec[i];
      g = log(1-theta*exp(-gamma * a));
      //o_g_a[i] = log(1-theta*exp(-gamma * a));
      o_g_a[i] = (g-g_min)/r;
    }

  return(0);
}		/* end of calc_g_a_S_R */

/*!
  \brief Calculates the non-linear g-function 
  \author Hanne Rognebakke
*/
int calc_g_a_log_lin(int i_ncat,double *i_a_vec,double *i_par,double *o_g_a)
{
  int    i;
  double a,g_min,r;


  //g_min = log(i_a_vec[0]);
  //r = log(i_a_vec[i_ncat-1])-log(i_a_vec[0]);
  for(i=0;i<i_ncat;i++)
    {
      a = i_a_vec[i];
      //o_g_a[i] = (log(a)-g_min)/r;
      o_g_a[i] = log(a);
    }
  return(0);
}		/* end of calc_g_a_log_lin */



/*L:calc_g_a_polyn3*
________________________________________________________________

		calc_g_a_polyn3
________________________________________________________________

Name:		calc_g_a_polyn3
Syntax:		
Description:    Calculates polyn3 function
Side effects:   
Return value:   None
Global or static variables used: 
Example:
Linking:
Bugs:
Author:		Geir Storvik, UiO
Date:
Source file: Cvs: $
________________________________________________________________
*/
static int calc_g_a_polyn3(int i_ncat,double *i_a_vec,double *i_par,double *o_g)
{
  int    i;
  double a,a_norm,g;
  double beta1,beta2,beta3;

  beta1 = exp(i_par[0]);
  beta3 = i_par[1];
  beta2 = G_ONE+beta1-beta3;
  for(i=0;i<i_ncat;i++)
    {
      a = i_a_vec[i];
      a_norm = (log(a)-log(i_a_vec[0]))/
               (log(i_a_vec[i_ncat-1])-log(i_a_vec[0]));
      g = ((beta3*a_norm+beta2)*a_norm+beta3)*a_norm;
      o_g[i] = g;
    }

  return(0);
}		/* end of calc_g_a_polyn3 */



