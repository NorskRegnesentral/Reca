/*!
  \file caa_sample_multi.c
  \author Geir Storvik
  \brief Routines for sampling nonlinear part of age model

  This file contains a lot of routines that currently is not used but have
  been tried out in order to obtain better convergence and/or easier implementation
*/
#include "caa.h"
#include "caa_sample_multi.h"
#include "caa_read_write.h"
#include "caa_chol.h"
#include "caa_routines.h"
#include "caa_util.h"
//#include "caa_COST.h"

static GMRFLib_hidden_param_tp *s_hidden_par=NULL;
static GMRFLib_optimize_param_tp *s_optimize_param;
static int Loglik_poisson_age(double *logll,double *x,
			      int m,int idx,double *x_vec,char *arg);
static int sample_age_alpha_ages_given(int i_h,Age_struct *i_age,Data_age *i_D_age,
				       int i_force_acc,int *o_acc,int i_write_alpha);
static double age_haul_calc_posterior(int i_ncat,double *i_alpha,double *i_prob,
				      double *i_Ages,double *i_mu,double i_tau);
static int age_haul_find_mode(int i_ncat,double *i_Ages,double i_N,
			      double *i_mu,double i_tau,
			      double *x_alpha_opt,double *o_prob,double *o_log_opt,
			      double *o_grad,double **o_Hess);
static int age_haul_calc_prob(int i_ncat,double *i_alpha,double *o_prob);



static double  *s_mu;        /*!< Prior means for alpha's */
static double  *s_prob;      /*!< Probabilities for multinomial distribution */
static double  *s_alpha_opt; /*!< Optimal values of alpha's */
static double  *s_alpha_prop;/*!< Proposal values of alpha's */
static double  *s_grad;      /*!< Gradient vector for optimization of haul likelihood */
static double **s_Hess;      /*!< Hesse matrix for optimization of haul likelihood */ 
static double  *s_eps;       /*!< Working vector for sampling alpha's */
static double  *s_eps2;      /*!< Working vector for sampling alpha's */
static double  *s_delta;     /*!< Change in Newton-Raphson for finding mode */
static int      s_N_gauher;  /*!< Number of nodes in Gauss Hermite quadrature */
static double  *s_gauher_x;  /*!< Abscissas in Gauss Hermite quadrature */
static double  *s_gauher_w;  /*!< Weights in Gauss Hermite quadrature */

#define NGIBBS  10        /*! Number of Gibbs iterations for sampling alpha's */

//#define AGE_TEST 1

#ifdef LOG_FILE
extern FILE     *g_caa_log;
extern FILE     *g_caa_alpha;
#endif

/*!
  \author Geir Storvik
  \brief Allocates space for routines used to sample non-linear part of age-model

  This routine must be called before using routines in the caa_sample_multi.c file.
*/
int sample_multi_initialize(int i_ncat)
{
  int    i;
  double sum;

  if(0){
  GMRFLib_default_optimize_param(&s_optimize_param);

  /* sett default verdier */
  if (!s_hidden_par) GMRFLib_default_hidden_par(&s_hidden_par);
  s_hidden_par->cmeanmode   = GMRFLib_COND_MEAN;
  s_hidden_par->gaussapprox = 0;
  /* er denne lik 1, faas approximasjon A1 */
  s_hidden_par->modeoption  = GMRFLib_MODEOPTION_MODE;
  s_hidden_par->nantithetic = 4;
  s_hidden_par->neightype   = GMRFLib_NEIGHTYPE_GRAPH;
  s_hidden_par->norder      = 2;
  s_hidden_par->range       = 6.;
  s_hidden_par->nresolution = 10;
  /* evnt sett lik 10, = anntall regioner i spline */
  s_hidden_par->neighpar    = 1;
  /* 0 betyr ingen korreksjon for integral-leddet.
     1 betyr korreksjon for naboene i grafen, 
     2 naboene og naboenes nabo, etc...   */

  s_hidden_par->nsample     = 10;
  /* hvor mange sample som skal brukes.  */
  if(0)
    {/* approximasjon A2, finnes ved */
      s_hidden_par->neighpar = 0;
      s_hidden_par->nsample  = 0;
    }
  if(1)
    {/* approximasjon A3, billig	*/
      s_hidden_par->neighpar    = 2;
      s_hidden_par->nsample     = 1; /* kun forventningen brukes */
      s_hidden_par->nantithetic = 0; /* og da maa denne vaere null */
    }
  if(0)
    {/* approximasjon A3, bedre*/
      s_hidden_par->nantithetic = 4;  /* default verdi */
      s_hidden_par->neighpar    = 2;
      s_hidden_par->nsample     = 10;  
    }
  }

  s_mu = CALLOC(i_ncat,double);        // Free ok
  if(!s_mu)
    {
      write_warning("sample_multi_initialize:Error allocating s_mu\n");
      return(1);
    }
  s_prob = CALLOC(i_ncat,double);      // Free ok
  if(!s_prob)
    {
      write_warning("sample_multi_initialize:Error allocating s_prob\n");
      return(1);
    }
  s_alpha_opt = CALLOC(i_ncat,double); // Free ok
  if(!s_alpha_opt)
    {
      write_warning("sample_multi_initialize:Error allocating s_alpha_opt\n");
      return(1);
    }
  s_alpha_prop = CALLOC(i_ncat,double);// Free ok
  if(!s_alpha_prop)
    {
      write_warning("sample_multi_initialize:Error allocating s_alpha_prop\n");
      return(1);
    }
  s_grad = CALLOC(i_ncat,double);     // Free ok
  if(!s_grad)
    {
      write_warning("sample_multi_initialize:Error allocating s_grad\n");
      return(1);
    }
  s_Hess = Mmatrix_2d(0,i_ncat-1,0,i_ncat-1,sizeof(double),1); // Free ok
  if(!s_Hess)
    {
      write_warning("sample_multi_initialize:Error allocating s_Hess\n");
      return(1);
    }
  s_eps = CALLOC(i_ncat,double);      // Free ok
  if(!s_eps)
    {
      write_warning("sample_multi_initialize:Error allocating s_eps\n");
      return(1);
    }
  s_eps2 = CALLOC(i_ncat,double);     // Free ok
  if(!s_eps2)
    {
      write_warning("sample_multi_initialize:Error allocating s_eps2\n");
      return(1);
    }
  s_delta = CALLOC(i_ncat,double);    // Free ok
  if(!s_delta)
    {
      write_warning("sample_multi_initialize:Error allocating s_delta\n");
      return(1);
    }

  //Gauss hermite weights
  s_N_gauher = 30;
  s_gauher_x = CALLOC(s_N_gauher+1,double);  //Free ok
  s_gauher_w = CALLOC(s_N_gauher+1,double);  //Free ok
  gauher(s_gauher_x,s_gauher_w,s_N_gauher);

  if(0) 
    {
      sum = G_ZERO;
      for(i=1;i<=s_N_gauher;i++)
	sum += s_gauher_w[i];
      for(i=1;i<=s_N_gauher;i++)
	s_gauher_w[i] /= sum;
    }


  #ifdef DEBUG_MULTI_ALPHA
  g_caa_alpha = fopen("multi_alpha.dat","w");
  #endif

  return(0);
}		/* end of sample_multi_initialize */

/*!
  \author Geir Storvik
  \brief Reallocate space allocated by sample_multi_initialize
*/
int sample_multi_re_initialize()
{
  FREE(s_mu);
  FREE(s_prob);
  FREE(s_alpha_opt);
  FREE(s_alpha_prop);
  FREE(s_grad);
  Fmatrix_2d(&s_Hess[0][0],&s_Hess[0]);
  FREE(s_eps);
  FREE(s_eps2);
  FREE(s_delta);
  FREE(s_gauher_x);
  FREE(s_gauher_w);

  #ifdef DEBUG_MULTI_ALPHA
  fclose(g_caa_alpha);
  #endif
  return(0);
}		/* end of sample_multi_re_initialize */



/*!
  \author Geir Storvik
  \brief Find optimal haul values for h's given
*/
int age_haul_modes(int i_start_h,int i_stop_h,Age_struct *i_age,Data_age *i_D_age)
{
  int     a,err,h;
  double  N_h,tau;
  double  log_opt;

  double *Ages;
  Ages=CALLOC(i_D_age->glm->ncat,double);

  tau = i_age->par->tau_obs;

  /* Find mode */
  for(h=i_start_h;h<i_stop_h;h++)
    {
      N_h = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  N_h += i_D_age->Ages[h][a];
	  s_mu[a] = G_ZERO;
	  i_age->alpha[h][a] = log(i_D_age->Ages[h][a] + 0.01);
	}
      if(N_h > 0)
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    Ages[a] = i_D_age->Ages[h][a];

	  err = age_haul_find_mode(i_D_age->glm->ncat,Ages,N_h,s_mu,tau,
				   i_age->alpha[h],s_prob,&log_opt,s_grad,s_Hess);	
	  if(err)
	    {
	      write_warning("age_haul_modes:Error calling age_haul_find_mode\n");
	      return(err);
	    }
	}
    }
  FREE(Ages);
          
  return(0);
}		/* end of age_haul_modes */


/*!
  \author Geir Storvik
  \brief Sample ages based on empirical age-given-length distribution

  Assumes lengths are given by categories. For each length category posterior 
  probabilites are calculated using prior and length, and ages are sampled 
  from multinomial distribution.

  Stored in ages per haul.

  This routine is only used for initialization of the parameters to be simulated.
*/
int sample_ages_init(Data_orig *i_D_orig,Data_CC *i_D_CC,Data_age *i_D_age,
		     Data_lin *i_D_lga,Data_g_a *i_D_g_a,
		     Data_lin *i_D_lga_CC,Data_g_a *i_D_g_a_CC,int i_saveSim)
{
  int            a,f,h,ind_f,ncat,na,n_CCc,n_CCu,n_Sc,n_Su;
  int            a2,i,s,l_int,ind_a,N_int,season,nSeason;
  long          *ages;
  double       **P_al,***P_al_s;
  double         sum,lobs;
  FILE          *fp;

  if(i_saveSim)
    {
      fp = fopen("ages_miss_start.dat","w");
      printf("sample_ages_init: print simulated missing ages to file: ages_miss_start.dat\n");
      fprintf(fp,"haul age_ind season age_real lobs rep\n");
    }

  ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  nSeason = i_D_g_a->nSeason;

  ncat = i_D_age->glm->ncat;
  if(i_D_orig->coastal_cod)
    {
      na = ncat/2;
      n_CCc = 0;
      n_CCu = 0;
      n_Sc = 0;
      n_Su = 0;
    }
  else
    na = ncat;

  /* Initialize */
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      for(a=0;a<ncat;a++)
	{
	  i_D_age->Ages[h][a] = i_D_age->Ages_fix[h][a];
	}
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  i_D_lga->Ages[h][a] = i_D_lga->Ages_fix[h][a];
          i_D_lga->sum_by_cat[h][a] = i_D_lga->sum_by_cat_fix[h][a];
          i_D_lga->sqsum_by_cat[h][a] = i_D_lga->sqsum_by_cat_fix[h][a];
	}
    }
  if(i_D_orig->coastal_cod)
    {
      for(h=0;h<i_D_age->glm->nHaul;h++) 
	{
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      n_Sc += i_D_lga->Ages_fix[h][a];
	    }
	  for(a=0;a<i_D_g_a_CC->ncat;a++)
	    {
	      n_CCc += i_D_lga_CC->Ages_fix[h][a];
	    }
	}
      if(i_D_CC->class_error)
	{
	  i_D_CC->k1 = genbet(P_BETA_A+n_CCc, P_BETA_B);
	  i_D_CC->k2 = genbet(P_BETA_A+n_Sc, P_BETA_B);
	} // end if(class_error )
    }

  /* Simulate missing ages*/
  /* First making a transition matrix P(a|l) for a finite number of lengths */
  N_int = i_D_orig->n_int_len;
  /* Data combined for all seasons, for use if no observations */
  P_al = Mmatrix_2d(0,N_int,0,i_D_age->glm->ncat-1,sizeof(double),1);
  /* Data for each season */
  P_al_s = Mmatrix_3d(0,nSeason-1,0,N_int,0,i_D_age->glm->ncat-1,sizeof(double),1);// Free ok

  for(i=0;i<N_int+1;i++)
    for(a=0;a<i_D_age->glm->ncat;a++)
      P_al[i][a] = G_ZERO;

  for(s=0;s<nSeason;s++)
    for(i=0;i<N_int+1;i++)
      for(a=0;a<i_D_age->glm->ncat;a++)
	P_al_s[s][i][a] = G_ZERO;

  ind_f = 0;  
  for(h=0;h<i_D_age->glm->nHaul;h++)
    for(f=0;f<i_D_orig->nFishBoat[h];f++)
      {
	a = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
	lobs = i_D_orig->totlength[ind_f];
	season = i_D_orig->season[ind_f];
	if(a > -1000 && lobs > -1000.0)
	  {
	    l_int = 0;
	    while(lobs > i_D_orig->int_len_lim[l_int])
	      l_int++;
	    P_al[l_int][a] += (double) i_D_orig->replength[ind_f];
	    P_al_s[season-1][l_int][a] += (double) i_D_orig->replength[ind_f];
	  }
	ind_f++;
      }

  // Convert to probabilities
  for(i=0;i<N_int+1;i++)
    {
      sum = G_ZERO;
      for(a=0;a<i_D_age->glm->ncat;a++)
	sum += P_al[i][a];
      if(sum < 0.0001)
	{
	  if(i==0)
	    P_al[i][0] = G_ONE;
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al[i][a] = P_al[i-1][a];
	    }
	}
      else
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    P_al[i][a] /= sum;
	}
    }
  for(s=0;s<nSeason;s++)
    {
      for(i=0;i<N_int+1;i++)
	{
	  sum = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    sum += P_al_s[s][i][a];
	  if(sum<0.0001)//use P_al (over all seasons) if no observations in this season
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_s[s][i][a] = P_al[i][a];
	    }
	  else
	    {
	      for(a=0;a<i_D_age->glm->ncat;a++)
		P_al_s[s][i][a] /= sum;
	    }
	}
    }

  /* Start simulation of missing ages */
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      ind_f = i_D_orig->start_noAge[h];
      for(f=0;f<i_D_orig->num_noAge[h];f++)
	{
	  a = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
	  season = i_D_orig->season[ind_f];
	  lobs = i_D_orig->totlength[ind_f];
	  if(a < -1000 && lobs > -1000.0)
	    {
	      //Find right interval
	      l_int = 0;
	      while(lobs > i_D_orig->int_len_lim[l_int])
		l_int++;
	      //Simulate age
	      my_genmul(i_D_orig->replength[ind_f],P_al_s[season-1][l_int],i_D_age->glm->ncat,ages);
	      for(a2=0;a2<i_D_age->glm->ncat;a2++)
		{
		  if(ages[a2]>0)
		    {
		      i_D_age->Ages[h][a2] += (double) ages[a2];
		      ind_a = i_D_g_a->a2Age_vec[a2]+(season-1);
		      i_D_lga->Ages[h][ind_a] += (int) ages[a2];
		      i_D_lga->sum_by_cat[h][ind_a] += lobs*ages[a2];
		      i_D_lga->sqsum_by_cat[h][ind_a] += lobs*lobs*ages[a2];
		      if(i_saveSim)
			fprintf(fp,"%d %d %d %f %f %d\n",h,a2,season,
				i_D_g_a->a_vec[ind_a],lobs,(int) ages[a2]);
		    }
		}
	    }
	  ind_f++;
	}
    }
    
  #ifdef DEBUG_PROG
  int N=0;
  fprintf(stderr,"sample_ages_init:\nage:");
  for(a=0;a<i_D_age->glm->ncat;a++){
    fprintf(stderr,"%d ",(int) i_D_age->Ages[0][a]);
    N+=(int) i_D_age->Ages[0][a];
  }
  fprintf(stderr,"N=%d\nlga:",N);
  N=0;
  for(a=0;a<i_D_g_a->ncat;a++){
    fprintf(stderr,"%d ",i_D_lga->Ages[0][a]);
    N+=(int) i_D_lga->Ages[0][a];
  }
  fprintf(stderr,"N=%d\n",N);
  if(i_D_orig->coastal_cod){
    fprintf(stderr,"lga:");
    for(a=0;a<i_D_g_a->ncat;a++){
      fprintf(stderr,"%d ",i_D_lga_CC->Ages[0][a]);
      N+=(int) i_D_lga_CC->Ages[0][a];
    }
    fprintf(stderr,"N=%d\n",N);
  }
  #endif

  
  if(i_saveSim)
    fclose(fp);

  // Free allocated memory
  FREE(ages);
  Fmatrix_2d(&P_al[0][0],&P_al[0]);
  Fmatrix_3d(&P_al_s[0][0][0],&P_al_s[0][0],&P_al_s[0]);

  return(0);
}		/* end of sample_ages_init */



/*!
  \author Geir Storvik and Hanne Rognebakke
  \brief Samples missing ages.

  Here all data are assumed to be long strings of  age and length. 

  In order to speed up computation, the fish are assumed ordered in hauls. 
  Further, inside each haul, the lenghts are assumed ordered in increasing values. 
  Then the length-range is divided into a finite number of intervals in which the 
  age-probabilities are assumed constant for all length-values inside an interval.
*/
int sample_ages(Data_orig *i_D_orig,Data_CC *i_D_CC,
		    Age_struct *i_age,Data_age *i_D_age,
		    LW_struct *i_length,Data_lin *i_D_lga,Data_g_a *i_D_g_a, 
		    LW_struct *i_length_CC,Data_lin *i_D_lga_CC,Data_g_a *i_D_g_a_CC, 
		    int i_saveSim, int i_it)
{
  int            a,f,h,ind_a,ind_f,n_sim,cum_fish,aobs,season;
  int            n_notsim,sim,cc,ncat_age,type=0;
  int            n_CCc,n_CCu,n_Sc,n_Su,err=0;
  double         lobs,lobs_prev,r,sum_p,lstart,lend;
  double        *p,*p2,*mu,*beta,*sigma;
  long          *ages;
  FILE          *fp;
  char           buffer[MAX_STR];


  p = CALLOC(i_D_age->glm->ncat,double);       // Free ok
  p2 = CALLOC(i_D_age->glm->ncat,double);       // Free ok
  ages = CALLOC(i_D_age->glm->ncat,long);      // Free ok
  mu = CALLOC(i_D_age->glm->ncat,double);      // Free ok
  beta = CALLOC(i_D_lga->glm->nxcov,double);      // Free ok      
  sigma = CALLOC(i_D_age->glm->ncat,double);      // Free ok

  ncat_age = i_D_age->glm->ncat;
  cc = i_D_orig->coastal_cod;
  sim = 0;
  if(i_age->age_errors)
    {
      sim = 1;
    }

  if(i_saveSim)
    {
      sprintf(buffer,"ages_miss_%d.dat",i_it);
      fp = fopen(buffer,"w");
      fprintf(stderr,"sample_ages: print simulated ages to file %s\n",buffer);
      fprintf(fp,"haul age_ind season age_real lobs rep\n");
    }

  for(a=0;a<ncat_age;a++)
    sigma[a] = G_ONE/sqrt(i_length->par->tau_obs);
  if(cc)
    {
      ncat_age = (int) i_D_age->glm->ncat/2;
      for(a=ncat_age;a<i_D_age->glm->ncat;a++)
	sigma[a] = G_ONE/sqrt(i_length_CC->par->tau_obs);
      if(i_D_CC->class_error)
	{
	  sim = 1;
	}
      n_CCc = 0;
      n_CCu = 0;
      n_Sc = 0;
      n_Su = 0;
    }

  cum_fish = 0;
  n_sim = 0;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    {
      /* Only simulate p(a) if ages are missing or ages observed with error */
      if(i_D_orig->num_noAge[h]>0 || sim==1)
	{
	  /* Start by initializing sufficient statistics */ 
	  err = init_suff_stat_sim_age(h,i_D_age,i_D_lga,i_D_lga_CC,i_D_g_a->ncat,sim,cc);
	  if(err)
	    {
	      write_warning("sample_ages:Error calling init_suff_stat_sim_age\n");
	      return(err);
	    }
	  if(cc)
	    {
	      for(a=0;a<ncat_age;a++)
		n_Sc += i_D_age->Ages[h][a];
	      for(a=ncat_age;a<i_D_age->glm->ncat;a++)
		n_CCc += i_D_age->Ages[h][a];		
	    }
	  /* find prior age-probabilities */
	  sum_p = G_ZERO;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      p[a] = exp(i_age->alpha[h][a]);
	      sum_p += p[a];
	    }
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    {
	      p[a] /= sum_p;
	    }

	  if(sim) /* Simulate new ages */
	    {
	      /* Both age and length observed */
	      for(f=0;f<(i_D_orig->nFishBoat[h]-i_D_orig->num_noAge[h]);f++)
		{
		  ind_f = cum_fish+f;
		  aobs = i_D_orig->totage[ind_f]-i_D_age->a_vec[0];
		  lobs = i_D_orig->totlength[ind_f];
		  lstart = i_D_orig->lstart[ind_f];
		  lend = i_D_orig->lend[ind_f];
		  if(cc)
		    {
		      type = i_D_orig->tottype[ind_f] ;
		    }
		  season = i_D_orig->season[ind_f];
		  if(cc)
		    err = calc_int_slp(h,season,ncat_age,i_D_lga->glm,i_length->par,i_D_g_a,
				       cc,i_D_lga_CC->glm,i_length_CC->par,i_D_g_a_CC,mu);
		  else
		    err = calc_int_slp(h,season,ncat_age,i_D_lga->glm,i_length->par,i_D_g_a,
				       cc,NULL,NULL,NULL,mu);
		  /* Calculate probabilities */
		  err = calc_prob_age(i_D_CC,cc,i_age,aobs,lobs,lstart,lend,type,ncat_age,mu,sigma,p,p2,&sum_p);
		  if(err)
		    {
		      write_warning("sample_ages:Error calling calc_prob_age\n");
		      return(err);
		    }

		  if(sum_p<=0)
		    {
		      fprintf(stderr,"WARNING: p=0 when simulating new ages \n");
		      n_notsim += i_D_orig->replength[ind_f];
		    }
		  else
		    {
		      my_genmul(i_D_orig->replength[ind_f],p2,i_D_age->glm->ncat,ages);
		      n_sim+=i_D_orig->replength[ind_f];
		      for(a=0;a<i_D_age->glm->ncat;a++)
			{
			  r = (double) ages[a];
			  if(r > 0)
			    {
			      i_D_age->Ages[h][a] += r;
			      if(cc & (a>=ncat_age))
				{
				  ind_a = i_D_g_a_CC->a2Age_vec[a]+(season-1);
				  i_D_lga_CC->Ages[h][ind_a] += ages[a];
				  i_D_lga_CC->sum_by_cat[h][ind_a] += r*lobs;
				  i_D_lga_CC->sqsum_by_cat[h][ind_a] += r*lobs*lobs;
				  if(type==1)
				    n_CCc += ages[a];
				  else if(type==2)
				    n_CCu += ages[a];
				}
			      else
				{
				  ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
				  i_D_lga->Ages[h][ind_a] += ages[a];
				  i_D_lga->sum_by_cat[h][ind_a] += r*lobs;
				  i_D_lga->sqsum_by_cat[h][ind_a] += r*lobs*lobs;
				  if(cc){
				    if(type==5)
				      n_Sc += ages[a];
				    else if(type==4)
				      n_Su += ages[a];
				  }
				}
			      if(i_saveSim)
				fprintf(fp,"%d %d %d %f %f %d\n",h,a,season,i_D_g_a->a_vec[ind_a],lobs,(int) r);
			      //if(h==0)
			      //fprintf(stderr,"%d %d age=%d type=%d,ind_a=%d,lga=%d,lgaCC=%d,lobs=%f r=%d,nCCc=%d,nSc=%d\n",h,a,(int) i_D_age->Ages[h][a],
			      //					type,ind_a,i_D_lga->Ages[h][ind_a],i_D_lga_CC->Ages[h][ind_a],lobs,(int) r,n_CCc,n_Sc);
			    }
			  else if(r < 0)
			    {
			      write_warning("sample_ages: Negative ages generated");
			      fprintf(stderr,"h=%d,f=%d:ages=%f,replength=%d,p2=%f,A_Nneigh=%d,aobs=%d\n",
				      h,f,r,i_D_orig->replength[ind_f],p2[(int)r],i_age->A_Nneigh[aobs],aobs);
			      return(1);
			    }
			}
		    }
		} // end  for(f=0;f<(i_D_orig->nFishBoat[h]-i_D_orig->num_noAge[h]);f++)
	      
	      /* Copy to fixed version - to be used in lga model */
	      for(a=0;a<i_D_age->glm->ncat;a++)
		i_D_age->Ages_fix[h][a] = i_D_age->Ages[h][a];
	      for(a=0;a<i_D_g_a->ncat;a++)
		{
		  i_D_lga->Ages_fix[h][a] = i_D_lga->Ages[h][a];
		  i_D_lga->sum_by_cat_fix[h][a] = i_D_lga->sum_by_cat[h][a];
		  i_D_lga->sqsum_by_cat_fix[h][a] = i_D_lga->sqsum_by_cat[h][a];
		}
	      if(cc)
		{
		  for(a=0;a<i_D_g_a_CC->ncat;a++)
		    {
		      i_D_lga_CC->Ages_fix[h][a] = i_D_lga_CC->Ages[h][a];
		      i_D_lga_CC->sum_by_cat_fix[h][a] = i_D_lga_CC->sum_by_cat[h][a];
		      i_D_lga_CC->sqsum_by_cat_fix[h][a] = i_D_lga_CC->sqsum_by_cat[h][a];
		    }
		}
	    } //end if(sim)

	  ind_f = i_D_orig->start_noAge[h];
	  lobs_prev = -1;
	  n_notsim = 0;

	  
	  /* Loop through all non-aged fish */
	  for(f=0;f<i_D_orig->num_noAge[h];f++)
	    {
	      //printf("h=%d,ind_f=%d,f=%d,rep[%d]=%d\n",h,ind_f,f,ind_f+f,i_D_orig->replength[ind_f+f]);
	      if(i_D_orig->replength[ind_f+f]>0)
		{
		  lobs = i_D_orig->totlength[ind_f+f];
		  season = i_D_orig->season[ind_f+f];
		  if(fabs(lobs-lobs_prev)>1e-8)
		    {
		      //New calculation of age-probabilities needs to be performed
		      /* Calculate probabilities */
		      if(cc)
			err = calc_int_slp(h,season,ncat_age,i_D_lga->glm,i_length->par,i_D_g_a,
					   cc,i_D_lga_CC->glm,i_length_CC->par,i_D_g_a_CC,mu);
		      else
			err = calc_int_slp(h,season,ncat_age,i_D_lga->glm,i_length->par,i_D_g_a,
					   cc,NULL,NULL,NULL,mu);
		      lstart = i_D_orig->lstart[ind_f+f];
		      lend = i_D_orig->lend[ind_f+f];
		      sum_p = G_ZERO;
		      for(a=0;a<i_D_age->glm->ncat;a++)
			{
			  p2[a] = p[a]*(pnorm(lend,mu[a],sigma[a])-pnorm(lstart,mu[a],sigma[a]));
			  sum_p += p2[a];
			}
		      for(a=0;a<i_D_age->glm->ncat;a++)
			p2[a] /= sum_p;
		      lobs_prev = lobs;
		    }
		  if(sum_p<=0)
		    {
		      fprintf(stderr,"WARNING: prob=0: h=%d,f=%d: lobs=%f, mu(0)=%f, mu(age_max)=%f, sum_p=%f, rep=%d\n",
			     h,(ind_f+f),lobs,mu[0],mu[i_D_age->glm->ncat-1],sum_p,i_D_orig->replength[ind_f+f]);
		      n_notsim += i_D_orig->replength[ind_f+f];
		    }
		  else
		    {
		      my_genmul(i_D_orig->replength[ind_f+f],p2,i_D_age->glm->ncat,ages);
		      n_sim+=i_D_orig->replength[ind_f+f];
		      for(a=0;a<i_D_age->glm->ncat;a++)
			{
			  r = (double) ages[a];
			  if(r>0)
			    {
			      if(i_it>0)
				i_D_age->Ages[h][a] += r;
			      if(cc & (a>=ncat_age))
				{
				  ind_a = i_D_g_a_CC->a2Age_vec[a]+(season-1);
				  i_D_lga_CC->Ages[h][ind_a] += (int) r;
				  i_D_lga_CC->sum_by_cat[h][ind_a] += r*lobs;
				  i_D_lga_CC->sqsum_by_cat[h][ind_a] += r*lobs*lobs;
				}
			      else
				{
				  ind_a = i_D_g_a->a2Age_vec[a]+(season-1);
				  i_D_lga->Ages[h][ind_a] += (int) r;
				  i_D_lga->sum_by_cat[h][ind_a] += r*lobs;
				  i_D_lga->sqsum_by_cat[h][ind_a] += r*lobs*lobs;
				}
			      if(i_saveSim)
				fprintf(fp,"%d %d %d %f %f %d\n",h,a,season,i_D_g_a->a_vec[ind_a],lobs,(int)r);
			    }
			}
		    }		  
		}
	    } //end for(f=0...)
	  if(n_notsim>0)
	    fprintf(stderr,"h=%d: number of fish where age not simulated = %d\n",h,n_notsim);	  
	}//end if(i_D_orig->num_noAge[h]>0 || i_age->age_errors==1)
      else
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_D_age->Ages[h][a] = i_D_age->Ages_fix[h][a]+(double)i_D_age->n_h[h]*i_age->delta_age;
	}

      cum_fish += i_D_orig->nFishBoat[h];
    } //end for(h=0...)


  #ifdef DEBUG_PROG
  printf("sample_ages: n_sim=%d\n",n_sim);
  #endif

  if(cc)
    {
      if(i_D_CC->class_error)
	{
	  i_D_CC->k1 = genbet(P_BETA_A+n_CCc, P_BETA_B+n_CCu);
	  i_D_CC->k2 = genbet(P_BETA_A+n_Sc, P_BETA_B+n_Su);
          #ifdef DEBUG_PROG
	  printf("k1=%f, Beta_A=%f,n_CCc=%d, Beta_B=%f,n_CCu=%d\n",i_D_CC->k1,P_BETA_A,n_CCc,P_BETA_B,n_CCu);
	  printf("k2=%f, Beta_A=%f,n_Sc=%d, Beta_B=%f,n_Su=%d\n",i_D_CC->k2,P_BETA_A,n_Sc,P_BETA_B,n_Su);
	  #endif
	}
    }

  if(i_saveSim)
    fclose(fp);


  // Free memory allocated in this routine
  FREE(p);
  FREE(p2);
  FREE(mu);
  FREE(ages);
  FREE(beta);
  FREE(sigma);

  return(0);
}		/* end of sample_ages */


/*!
  \author Hanne Rognebakke
  \brief  Initializing sufficient statistics to be used in ::sample_ages
*/
int init_suff_stat_sim_age(int i_h,Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_lga_CC,
			   int i_ga_ncat,int i_sim,int i_cc)
{
  int a;
  
  if(i_sim) /* Age errors: All simulated */
    {
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  i_D_age->Ages[i_h][a] = G_ZERO;
	}
      for(a=0;a<i_ga_ncat;a++)
	{
	  i_D_lga->Ages[i_h][a] = 0;
	  i_D_lga->sum_by_cat[i_h][a] = G_ZERO;
	  i_D_lga->sqsum_by_cat[i_h][a] = G_ZERO;
	}
      if(i_cc)
	{
	  for(a=0;a<i_ga_ncat;a++)
	    {
	      i_D_lga_CC->Ages[i_h][a] = 0;
	      i_D_lga_CC->sum_by_cat[i_h][a] = G_ZERO;
	      i_D_lga_CC->sqsum_by_cat[i_h][a] = G_ZERO;
	    }
	}
    }
  else  /* No age errors: Keep aged fish */
    {
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  i_D_age->Ages[i_h][a] = i_D_age->Ages_fix[i_h][a];
	}
      for(a=0;a<i_ga_ncat;a++)
	{
	  i_D_lga->Ages[i_h][a] = i_D_lga->Ages_fix[i_h][a];
	  i_D_lga->sum_by_cat[i_h][a] = i_D_lga->sum_by_cat_fix[i_h][a];
	  i_D_lga->sqsum_by_cat[i_h][a] = i_D_lga->sqsum_by_cat_fix[i_h][a];
	}
      if(i_cc)
	{
	  for(a=0;a<i_ga_ncat;a++)
	    {
	      i_D_lga_CC->Ages[i_h][a] = i_D_lga_CC->Ages_fix[i_h][a];
	      i_D_lga_CC->sum_by_cat[i_h][a] = i_D_lga_CC->sum_by_cat_fix[i_h][a];
	      i_D_lga_CC->sqsum_by_cat[i_h][a] = i_D_lga_CC->sqsum_by_cat_fix[i_h][a];
	    }
	}
    }

  return(0);
}		/* end of init_suff_stat_sim_age */


/*!
  \author Hanne Rognebakke
  \brief Calculate new age probability for simulating ages observed with error or classification error
*/
int calc_prob_age(Data_CC *i_D_CC,int i_cc,Age_struct *i_age,int i_aobs,
		  double i_lobs,double i_lstart,double i_lend,int i_type,int i_ncat_age,
		  double *i_mu,double *i_sigma,double *i_p,double *o_p2,double *o_sum_p)
{
  int      i,a,A_Nneigh,ncat_tot;
  double   sum_p,A2A;
  double   pobs,ptrue,pobs2,ptrue2,prob;
  char     buffer[MAX_STR];

  if(i_cc)
    ncat_tot = 2*i_ncat_age;
  else
    ncat_tot = i_ncat_age;

  for(i=0;i<ncat_tot;i++)
    o_p2[i] = G_ZERO;
  sum_p = G_ZERO;

  if(i_age->age_errors)
    A_Nneigh = i_age->A_Nneigh[i_aobs];
  else
    A_Nneigh = 1;

  for(i=0;i<A_Nneigh;i++)
    {
      if(i_age->age_errors)
	{
	  a = i_age->A_neigh[i_aobs][i];
	  A2A = i_age->A2A[i_aobs][a];
	}
      else
	{
	  a = i_aobs;
	  A2A = 1;
	}
      if(i_cc & i_D_CC->class_error)
	{
	  if(i_type == 1) //coastal cod
	    {
	      pobs = i_D_CC->ptype1_CC1[a];
	      ptrue = i_p[a+i_ncat_age]*i_D_CC->k1;
	      pobs2 = i_D_CC->ptype1_S5[a];
	      ptrue2 = i_p[a]*i_D_CC->k2;
	      prob = pobs*ptrue/(pobs*ptrue+pobs2*ptrue2);
	    }
	  else if(i_type == 2) //uncertain coastal cod
	    {
	      pobs = i_D_CC->ptype2_CC2[a];
	      ptrue = i_p[a+i_ncat_age]*(1-i_D_CC->k1);
	      pobs2 = i_D_CC->ptype2_S4[a];
	      ptrue2 = i_p[a]*(1-i_D_CC->k2);
	      prob = pobs*ptrue/(pobs*ptrue+pobs2*ptrue2);
	    }
	  else if(i_type == 4) //uncertain skrei
	    {
	      pobs = i_D_CC->ptype4_S4[a];
	      ptrue = i_p[a]*(1-i_D_CC->k2);
	      pobs2 = i_D_CC->ptype4_CC2[a];
	      ptrue2 = i_p[a+i_ncat_age]*(1-i_D_CC->k1);
	      prob = 1-pobs*ptrue/(pobs*ptrue+pobs2*ptrue2);
	    }
	  else if(i_type == 5) //skrei
	    {
	      pobs = i_D_CC->ptype5_S5[a];
	      ptrue = i_p[a]*i_D_CC->k2;
	      pobs2 = i_D_CC->ptype5_CC1[a];
	      ptrue2 = i_p[a+i_ncat_age]*i_D_CC->k1;
	      prob = 1-pobs*ptrue/(pobs*ptrue+pobs2*ptrue2);
	    }
	  else
	    {
	      fprintf(stderr,"calc_prob_age: Something is wrong, type=%d\n",i_type);
	      prob = 0.5;
	    }
	  if(prob<0 || prob>1)
	    {
	      sprintf(buffer,"calc_prob_age: Something is wrong, type=%d,prob=%f\n",i_type,prob);
	      write_warning(buffer);
	      return(1);
	    }
	}
      else if(i_cc)
	{
	  if(i_type==1 || i_type == 2)
	    {
	      prob = 1.0;
	    }
	  else if(i_type==4 || i_type==5)
	    {
	      prob = 0.0;
	    }
	  else
	    fprintf(stderr,"calc_prob_age: Something is wrong, type=%d\n",i_type);
	}
      else
	{
	  prob = 0;
	}

      o_p2[a] = (1-prob)*i_p[a]*A2A*((pnorm(i_lend,i_mu[a],i_sigma[a])-pnorm(i_lstart,i_mu[a],i_sigma[a]))+0.00001);
      sum_p += o_p2[a];
      if(i_cc)
	{
	  o_p2[a+i_ncat_age] = prob*i_p[a+i_ncat_age]*A2A*
	    (pnorm(i_lend,i_mu[a+i_ncat_age],i_sigma[a+i_ncat_age])-pnorm(i_lstart,i_mu[a+i_ncat_age],i_sigma[a+i_ncat_age]));
	  sum_p += o_p2[a+i_ncat_age];
	}
    }
  if(sum_p<0)
    {
      sprintf(buffer,"calc_prob_age: Something is wrong, sum_p=%f\n",sum_p);
      write_warning(buffer);
      return(1);
    }
  for(i=0;i<ncat_tot;i++)
    {
      o_p2[i] /= sum_p;
    }
  *o_sum_p = sum_p;

  return(0);
}		/* end of calc_prob_age */

/*!
  \author Hanne Rognebakke
  \brief Calculate new age probability for simulating ages observed with error or classification error
*/
int calc_int_slp(int i_h, int i_season, int i_ncat_age, Data_glm *i_glm, Eff_str *i_par, Data_g_a *i_D_g_a,
		 int i_cc, Data_glm *i_glmCC, Eff_str *i_parCC, Data_g_a *i_D_g_a_CC, double *o_mu)
{
  int a,i;
  double *beta;

  beta = CALLOC(i_glm->nxcov,double);      // Free ok      

  //  calc_int_slp(h,season,i_D_lga->glm,i_length->par,i_D_lga_CC,i_length_CC->par);
  
  /* Find intercept and slope for lga model */
  for(i=0;i<i_glm->nxcov;i++)
    beta[i] = calc_eff(i_glm->xcov[i],i_par->eff[0][i],i_h);
  for(a=0;a<i_ncat_age;a++)
    {
      o_mu[a] = beta[0] + beta[1]*i_D_g_a->g_a[i_D_g_a->a2Age_vec[a]+i_season-1];
    }
  if(i_cc)/* Find intercept and slope for lga model CC*/
    {
      for(i=0;i<i_glmCC->nxcov;i++)
	beta[i] = calc_eff(i_glmCC->xcov[i],i_parCC->eff[0][i],i_h);
      for(a=0;a<i_ncat_age;a++)
	o_mu[a+i_ncat_age] = beta[0] + beta[1]*i_D_g_a_CC->g_a[i_D_g_a_CC->a2Age_vec[a]+i_season-1];
    }
  
   return(0);
}		/* end of calc_prob_age */


/*!
  \author Geir Storvik
  \brief Calculate sufficient statistics for age model

*/
int make_suff_age(int i_ncat,Age_struct *i_age,Data_age *i_D_age,double *i_haulweight,int i_start_h)
{
  int      a,h;

  if(i_D_age->glm->nxcov==1)
    {
      for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
	{
	  for(a=0;a<i_ncat;a++)
	    {
	      i_D_age->glm->beta_hat[h][a][0] = i_age->alpha[h][a];
	      i_D_age->glm->suff[h][0][0] = G_ONE;	      
	    }
	}
    }
  else
    {
      write_warning("age->nxcov different from 1 is not implemented\n");
      return(1);
    }


  return(0);
}               /* end of make_suff_age */



/*!
  \author Geir Storvik
  \brief Calculate sufficient statistics for length given age model

  The sufficient statistics are conditional on the \f$g(a)\f$ function.
  If i_use_sim_age=1, the data were ages are simulated is included, else
  if 0, only observed age data is used.

  Sufficient statistics are for each haul
  - the number of fish per age
  - the sum of $g(a)$
  - the square sum of \f$g(a)\f$
  - Least squares intercept, slope and sum of squares
*/
int make_suff_lga(Data_lin *i_D_lga,Data_g_a *i_D_g_a,int i_start_h,int i_use_sim_age)
{
  int      a,h;
  double   g_a,maxAges,N_h,N_h_a,sum_g,sum_g2,sum_l,sum_l2,sum_gl,beta0,beta1,ssq;
  int     *Ages;
  double  *sum_by_cat,*sqsum_by_cat;

  Ages = CALLOC(i_D_g_a->ncat,int);
  sum_by_cat = CALLOC(i_D_g_a->ncat,double);
  sqsum_by_cat = CALLOC(i_D_g_a->ncat,double);

  for(h=i_start_h;h<i_D_lga->glm->nHaul;h++)
    {
      if(i_use_sim_age)
	{
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      Ages[a] = i_D_lga->Ages[h][a];
	      //if(Ages[a]<0)
	      //printf("simage=%d,h=%d,Ages[%d]=%d \n",i_use_sim_age,h,a,Ages[a]); 
	      sum_by_cat[a] = i_D_lga->sum_by_cat[h][a];
	      sqsum_by_cat[a] = i_D_lga->sqsum_by_cat[h][a];
	    }
	}
      else // use only observed ages
	{
	  for(a=0;a<i_D_g_a->ncat;a++)
	    {
	      Ages[a] = i_D_lga->Ages_fix[h][a]; 
	      //  if(Ages[a]<0)
	      //printf("simage=%d,h=%d,Ages[%d]=%d \n",i_use_sim_age,h,a,Ages[a]); 
	      sum_by_cat[a] = i_D_lga->sum_by_cat_fix[h][a];
	      sqsum_by_cat[a] = i_D_lga->sqsum_by_cat_fix[h][a];
	    }
	}
      sum_g = G_ZERO;
      sum_g2 = G_ZERO;
      sum_gl = G_ZERO;
      sum_l = G_ZERO;
      sum_l2 = G_ZERO;
      maxAges = G_ZERO;
      N_h = G_ZERO;
      for(a=0;a<i_D_g_a->ncat;a++)
	{
	  g_a = i_D_g_a->g_a[a];
          N_h_a = (double) Ages[a];
          sum_g +=  N_h_a * g_a;
          maxAges = max(maxAges,N_h_a);
          sum_g2 +=  N_h_a * g_a * g_a;
          sum_gl += g_a * sum_by_cat[a];
	  sum_l += sum_by_cat[a];
	  sum_l2 += sqsum_by_cat[a];
          N_h += N_h_a;
	}
      if(N_h < 0.0001)
	{
	  /* No data */
          beta0 = G_ZERO;
          beta1 = G_ZERO;
          ssq = G_ZERO;
	}
      else if(fabs(maxAges-N_h)<0.001)
	{ /* Only one age group sampled */
          beta0 = sum_l/N_h;
          beta1 = G_ZERO;
	  ssq = sum_l2+beta0*beta0*N_h-G_TWO*beta0*sum_l+0.00000001;
	}
      else if(fabs(sum_g2-sum_g)<0.0001)
	{  // Only samples where g-function is equal
          beta0 = sum_l/N_h;
          beta1 = G_ZERO;
	  ssq = sum_l2+beta0*beta0*N_h-G_TWO*beta0*sum_l+0.00000001;
	}
      else if(fabs(sum_g*sum_g-N_h*sum_g2)<0.0000000000001)
	{  // Only samples where g-function is equal
          beta0 = sum_l/N_h;
          beta1 = G_ZERO;
	  ssq = sum_l2+beta0*beta0*N_h-G_TWO*beta0*sum_l+0.00000001;
	}
      else
	{
	  beta1 = (sum_l*sum_g-N_h*sum_gl)/(sum_g*sum_g-N_h*sum_g2);
	  beta0 = (sum_l-beta1*sum_g)/N_h;
	  ssq = sum_l2+beta0*beta0*N_h+beta1*beta1*sum_g2-
	    G_TWO*(beta0*sum_l+beta1*sum_gl-beta0*beta1*sum_g)+0.00000001;
	}
      if(ssq < G_ZERO || !(beta0 > -9999999.99 && beta0 < 999999999.99) ||
                         !(beta0 > -9999999.99 && beta0 < 999999999.99))
	{
	  printf("h=%d,beta0=%lf,beta1=%lf,ssq=%lf\n",h,beta0,beta1,ssq);
          for(a=0;a<i_D_g_a->ncat;a++)
	    printf("%d %f %f %f\n",Ages[a],g_a,sum_by_cat[a],sqsum_by_cat[a]);
	  write_warning("make_suff_lga:Something is wrong\n");
	  return(1);
	}
      i_D_lga->glm->beta_hat[h][0][1] = beta1;
      i_D_lga->glm->beta_hat[h][0][0] = beta0;
      i_D_lga->glm->ssq[h] = ssq;
      i_D_lga->glm->suff[h][0][0] = N_h;
      i_D_lga->glm->suff[h][0][1] = sum_g;
      i_D_lga->glm->suff[h][1][0] = sum_g;
      i_D_lga->glm->suff[h][1][1] = sum_g2;
if(0){      
      if(i_D_lga->glm->suff[h][1][0]<0)
	{
	  printf("suff[%d][1][0]=%f, %f\n",h,i_D_lga->glm->suff[h][1][0],sum_g);
	  for(a=0;a<i_D_g_a->ncat;a++) 
	    {
	      sum_g +=  N_h_a * g_a;
	      printf("h=%d,a=%d,N_h_a=%f,g=%f\n",h,a,N_h_a,i_D_g_a->g_a[a]);
	    }
	}
 }
      //if(h<10)
      //printf("h=%d, beta0=%12.10f, beta1=%12.10f, suff00=%10.8f, suff01=%10.8f, suff11=%10.8f\n",
      //       h,beta0,beta1,N_h,sum_g,sum_g2);

      if(i_D_lga->glm->nxcov==3)
	{
	  i_D_lga->glm->suff[h][0][2] = N_h * i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][1][2] = i_D_lga->glm->suff[h][0][1]*i_D_lga->haulweight[h];
	  i_D_lga->glm->suff[h][2][2] = N_h * i_D_lga->haulweight[h]*i_D_lga->haulweight[h];
	  
	  i_D_lga->glm->suff[h][2][0] = i_D_lga->glm->suff[h][0][2];
	  i_D_lga->glm->suff[h][2][1] = i_D_lga->glm->suff[h][1][2];
	  i_D_lga->glm->beta_hat[h][0][2] = G_ZERO;
	}
    }
  FREE(Ages);
  FREE(sum_by_cat);
  FREE(sqsum_by_cat);

  return(0);
}               /* end of make_suff_lga */


/*!
  \author Hanne Rognebakke
  \brief Copy sufficient statistics for length given age model for fixed data
*/
int copy_suff_lga_fix(Data_lin *i_D_lga,int i_start_h)
{
  int h;

  for(h=i_start_h;h<i_D_lga->glm->nHaul;h++)
    {
      i_D_lga->glm->beta_hat[h][0][0] = i_D_lga->glm->beta_hat_fix[h][0][0];
      i_D_lga->glm->beta_hat[h][0][1] = i_D_lga->glm->beta_hat_fix[h][0][1];
      i_D_lga->glm->ssq[h] = i_D_lga->glm->ssq_fix[h];
      i_D_lga->glm->suff[h][0][0] = i_D_lga->glm->suff_fix[h][0][0];
      i_D_lga->glm->suff[h][0][1] = i_D_lga->glm->suff_fix[h][0][1];
      i_D_lga->glm->suff[h][1][0] = i_D_lga->glm->suff_fix[h][1][0];
      i_D_lga->glm->suff[h][1][1] = i_D_lga->glm->suff_fix[h][1][1];
      if(i_D_lga->glm->nxcov==3)
	{
	  i_D_lga->glm->suff[h][0][2] = i_D_lga->glm->suff_fix[h][0][2];
	  i_D_lga->glm->suff[h][1][2] = i_D_lga->glm->suff_fix[h][1][2];
	  i_D_lga->glm->suff[h][2][2] = i_D_lga->glm->suff_fix[h][2][2];
	  i_D_lga->glm->suff[h][2][0] = i_D_lga->glm->suff_fix[h][2][0];
	  i_D_lga->glm->suff[h][2][1] = i_D_lga->glm->suff_fix[h][2][1];
	  i_D_lga->glm->beta_hat[h][0][2] = G_ZERO;
	}
    }
  return(0);
}               /* end of copy_suff_lga_fix */



/*!
  \author Geir Storvik
  \brief Samples random effects in age model.

  Include a haul effect as random effect. Starts sampling for h=i_start_h

  The main simulation is performed per haul through the 
  sample_age_alpha_ages_given routine.
*/
int sample_age_alpha(Age_struct *i_age,Data_age *i_D_age,int i_start_h,int i_acc,
		     int i_it,int *acc_h,int i_write_alpha)
{
  int     err,h,acc,acc_tot,acc_min,acc_max;


  acc_min = 1000000;
  acc_max = 0;
  acc_tot = G_ZERO; 
  if(i_it==0)
    {
      for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
	acc_h[h] = 0;
    }
  for(h=i_start_h;h<i_D_age->glm->nHaul;h++)
    {
      // Alternative using Gibbs sampling on each alpha's separately
      // Used for testing the sample_age_alpha_ages_given routine
      //err = sample_age_alpha_ages_given_gibbs(h,i_age,i_D_age,i_acc,&acc);

      err = sample_age_alpha_ages_given(h,i_age,i_D_age,i_acc,&acc,i_write_alpha);
      if(err)
	{
	  write_warning("sample_age_alpha:Error calling sample_age_alpha_ages_given\n");
	  return(err);
	}

      acc_tot += acc;
      acc_h[h] += acc;
      acc_min = MIN(acc_min,acc_h[h]);
      acc_max = MAX(acc_max,acc_h[h]);
    }
  
  return(0);
}		/* end of sample_age_alpha */


/*!
  \author Geir Storvik
  \brief Routine for sampling alpha's inside a haul.

  This routine samples \f$\alpha_{h,a},a=1,...,A\f$ conditioned on 
  all fixed and random effects (except the haul effect), other hyper-parameters
  and observed ages. Note that simulating the \f$\alpha_{h,a}\f$'s is
  equivalent to sampling the haul effects.

  The simulation is performed through Metropolis-Hastings independence sampling
  steps by construting a proposal distribution through optimization of the
  conditional distribution. A Gaussian proposal with mean in the mode and
  covariance given from the hessian is used.

  Since the optimization requires some computer time,
  NGIBBS iterations of the Metropolis-Hastings independence sampling is performed.

  A problem with the independence sampling can be that the density for the
  backwards proposal is very small, making acceptance very unlikely. A possibility
  in that case is to scale the covariance matrix such that both jump forwards
  and backwards becomes more likely. This is implemented through the sc_old
  and sc_new variables. Currently these are set to 1, though.

  080401: Modified the routine to sample only the i_ncat-1 first alpha's under the
  constraint that all sum to zero. The age_haul_find_mode routine is changed accordingly.
*/
static int sample_age_alpha_ages_given(int i_h,Age_struct *i_age,Data_age *i_D_age,
				       int i_force_acc,int *o_acc,int i_write_alpha)
{
  int     a,a2,it,noacc,err=0;
  double  N_h,sum,u,tau,ssq_old,ssq_new;
  double  log_q_old,log_q_new,log_d_old,log_d_new,log_d_opt,log_old,log_new;
  double  sc_old,sc_new;

  double *Ages;
  Ages=CALLOC(i_D_age->glm->ncat,double);
  for(a=0;a<i_D_age->glm->ncat;a++)
    Ages[a] = i_D_age->Ages[i_h][a];

  i_D_age->glm->xcov[0]->n_cov--; //Not haul effect
  tau = i_age->par->tau_obs;
  /* Initializing */
  N_h = G_ZERO;
  sum = G_ZERO;
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      N_h += Ages[a];
      s_mu[a] = calc_eff_suff_age(i_D_age->glm->xcov[0],i_age->par,i_h,a);
      //    for(i=0;i<i_D_age->glm->nxcov;i++)
      //	s_mu[a] += calc_eff(i_D_age->glm->xcov[i],i_age->par->eff[a][i],i_h)*i_D_age->glm->suff[i_h][0][i];
    }
  i_D_age->glm->xcov[0]->n_cov++;

  /* Find mode */
  err = age_haul_find_mode(i_D_age->glm->ncat,Ages,N_h,s_mu,tau,
			   s_alpha_opt,s_prob,&log_d_opt,s_grad,s_Hess);
  if(err)
    {
      write_warning("sample_age_alpha_ages_given:Error calling age_haul_find_mode\n");
      return(err);
    }
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      if(!(s_alpha_opt[a] > -99999999999.99 && s_alpha_opt[a] < 999999999.99))
	{
	  printf("h=%d,tau=%lf\n",i_h,tau);
	  for(a2=0;a2<i_D_age->glm->ncat;a2++)
	    printf("a=%d,Ages=%f, mu=%lf,s_alpha_opt=%lf\n",
		   a2,Ages[a2],s_mu[a2],s_alpha_opt[a2]);
	  write_warning("sample_age_alpha_ages_given:something is wrong\n");
	  return(1);
	}
    }

  /* Calculate ssq for old sample */ 
  for(a=0;a<(i_D_age->glm->ncat-1);a++)
    s_eps2[a] = i_age->alpha[i_h][a]-s_alpha_opt[a];
  ssq_old = cholssq0(s_Hess,i_D_age->glm->ncat-1,s_eps2);

  /* Calculate proposal log-likelihood for old sample */
  err = age_haul_calc_prob(i_D_age->glm->ncat,i_age->alpha[i_h],s_prob);
  if(err)
    {
      write_warning("sample_age_alpha_ages_given:Error calling age_haul_calc_prob\n");
      return(err);
    }
  log_d_old = age_haul_calc_posterior(i_D_age->glm->ncat,i_age->alpha[i_h],s_prob,
				    Ages,s_mu,tau);

  // Make a proper scaling so that jump backwards become more likely
  if(log_d_opt > log_d_old)
    sc_old = ssq_old/(G_TWO*(log_d_opt-log_d_old));
  else
    sc_old = G_ONE;
  // Turning of the scaling...
  sc_old = G_ONE;
  log_q_old = -G_HALF*ssq_old/sc_old;
  log_old = log_d_old - log_q_old;

  /* Sample */
  
  noacc = 1;
  it = 0;
  while(it < NGIBBS)
    {
      /* Draw from proposal */
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	s_eps[a] = gennor(G_ZERO,G_ONE/sqrt(sc_old));
      chollTl0(s_Hess,i_D_age->glm->ncat-1,s_eps,s_eps2);  //chollTl?
      #ifdef DEBUG_HAUL_OPT
      int a2;
      FILE  *unit;
      unit = fopen("Hchol.dat","w");
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	{
	  for(a2=0;a2<(i_D_age->glm->ncat-1);a2++)
	    fprintf(unit,"%lf ",s_Hess[a][a2]);
	  fprintf(unit,"\n");
	}
      fclose(unit);
      unit = fopen("eps_eps2.dat","w");
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	fprintf(unit,"%lf %lf\n",s_eps[a],s_eps2[a]);
      fclose(unit);
      #endif       
      
      sum = G_ZERO;
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	{
	  s_alpha_prop[a] = s_alpha_opt[a] + s_eps2[a];
          if(!(s_alpha_prop[a] > -99999999999.99 && s_alpha_prop[a] < 999999999.99))
	    {
              printf("h=%d\n",i_h);
              for(a2=0;a2<i_D_age->glm->ncat;a2++)
		printf("Ages[%d]=%f, s_alpha_opt[%d]=%lf\n",
		       a2,Ages[a2],a2,s_alpha_opt[a]);
	      write_warning("sample_age_alpha_ages_given:something is wrong\n");
	      return(1);
	    }
	  sum += s_alpha_prop[a];
	}
      s_alpha_prop[i_D_age->glm->ncat-1] = -sum; 
      // Calculate age-probabilities for new sample
      err = age_haul_calc_prob(i_D_age->glm->ncat,s_alpha_prop,s_prob);
      log_d_new = age_haul_calc_posterior(i_D_age->glm->ncat,s_alpha_prop,s_prob,
					Ages,s_mu,tau);

      /* Calculate ssq for new sampl */
      ssq_new = G_ZERO;
      for(a=0;a<(i_D_age->glm->ncat-1);a++)
	ssq_new += s_eps[a]*s_eps[a];
      // Scaling for forward jump
      //sc_new = ssq_new/(G_TWO*(log_d_opt-log_d_new));
      // Scaling turned off....
      sc_new = G_ONE;
      /* Calculate proposal likelihoods */
      log_q_new = -G_HALF*ssq_new/sc_new;
      log_new = log_d_new - log_q_new;

      // M_H acceptance
      u = genunf(G_ZERO,G_ONE);
      if(i_force_acc || u < exp(log_new-log_old))
	{
	  /* Accept sample */
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_age->alpha[i_h][a] = s_alpha_prop[a];
	  log_old = log_new;
	  noacc = 0;
	}
      it++;
    }  

  for(a=0;a<i_D_age->glm->ncat;a++)
    i_age->par->eff[a][0][i_D_age->glm->xcov[0]->n_cov-1][i_h] = i_age->alpha[i_h][a];

  #ifdef WRITE_ALPHAS 
  if(i_write_alpha)
    {
      FILE   *unit;
      unit = fopen("alphas.bin","ab");
      fwrite(i_age->alpha[i_h],sizeof(double),i_D_age->glm->ncat,unit);
      fwrite(s_alpha_opt,sizeof(double),i_D_age->glm->ncat,unit);
      fclose(unit);
    }
  #endif

  err = age_haul_calc_prob(i_D_age->glm->ncat,i_age->alpha[i_h],s_prob);
  *o_acc = 1-noacc;  

  FREE(Ages);        
  return(0);
}		/* end of sample_age_alpha_ages_given */


/*!
  \author Geir Storvik
  \brief Finds optimal values of alpha for a given haul.

  Optimization is performed through Newton-Raphson steps with initial
  values based on observed ages and hyper-parameters (i.e not on old value
  which is important in order to make it possible to use in an independence sampler)

  080401: Optimizing on first A-1 alpha's, using that the last one is minus the sum 
  of the   other. Note that the covariance in the prior for the first A-1 alpha's 
  conditional on that the sum of all are zero is I+11^T.
  
  Returns both mode, gradient, hessian, the age-probabilities and 
  conditional log-posterior density
*/
static int age_haul_find_mode(int i_ncat,double *i_Ages,double i_N,
			      double *i_mu,double i_tau,
			      double *x_alpha_opt,double *o_prob,double *o_log_opt,
			      double *o_grad,double **o_Hess)
{
  int      it,more,err,a,a2;
  double   log_opt_old,log_opt,scale,mean_alpha,sum;

  // Initial values
  mean_alpha = G_ZERO;
  sum = G_ZERO;
  for(a=0;a<i_ncat;a++)
    {
      o_prob[a]=i_Ages[a]+0.00001;
      sum +=o_prob[a];
    }
  for(a=0;a<i_ncat;a++)
    {
      o_prob[a]/=sum;
      //x_alpha_opt[a] = i_mu[a]+i_Ages[a]*(G_ONE-o_prob[a])/i_tau;  // Old initial value
      x_alpha_opt[a] = i_mu[a]+(i_Ages[a]-i_N*o_prob[a])/i_tau;
      mean_alpha += x_alpha_opt[a];
    }
  mean_alpha /= (double) i_ncat;
  for(a=0;a<i_ncat;a++)
    x_alpha_opt[a] -= mean_alpha;
  err = age_haul_calc_prob(i_ncat,x_alpha_opt,o_prob);
  log_opt_old = age_haul_calc_posterior(i_ncat,x_alpha_opt,o_prob,i_Ages,i_mu,i_tau);

  it=0;
  more = 1;
  // Run maximum 40 iterations with NR on first i_ncat-1 alpha's
  // NB! The number of iterations must be large enough to find the proper alpha_opt value
  while(more && it < 40)
     {
       for(a=0;a<(i_ncat-1);a++)
	 {
	   o_grad[a] = -i_tau*(x_alpha_opt[a]-i_mu[a]-x_alpha_opt[i_ncat-1]+i_mu[i_ncat-1]) + 
	     (double) (i_Ages[a]-i_Ages[i_ncat-1])-i_N*(o_prob[a]-o_prob[i_ncat-1]);
	   o_Hess[a][a] = G_TWO*i_tau + i_N*o_prob[a]*(G_ONE-o_prob[a]) +
	                  i_N*o_prob[i_ncat-1]*(G_ONE-o_prob[i_ncat-1])+
                          G_TWO*i_N*o_prob[a]*o_prob[i_ncat-1];
	   for(a2=a+1;a2<(i_ncat-1);a2++)
	     {
	       o_Hess[a][a2] = i_tau + i_N*o_prob[i_ncat-1]-
                               i_N*(o_prob[a]-o_prob[i_ncat-1])*
                                   (o_prob[a2]-o_prob[i_ncat-1]);
	       o_Hess[a2][a] = o_Hess[a][a2];
	     }
	 }
       #ifdef DEBUG_HAUL_OPT
       FILE *unit;
       unit = fopen("Hess.dat","w");
       for(a=0;a<(i_ncat-1);a++)
	 {
	   for(a2=0;a2<(i_ncat-1);a2++)
	     fprintf(unit,"%lf ",o_Hess[a][a2]);
	   fprintf(unit,"\n");
	 }
       fclose(unit);
       #endif
       err = choldc0(o_Hess,i_ncat-1);
       if(err)
	 {
	   write_warning("sample_age_find_mode:Error calling choldc0\n");
	   return(err);
	 }
       cholsl0_old(o_Hess,i_ncat-1,o_grad,s_delta);
       #ifdef DEBUG_HAUL_OPT
       unit = fopen("grad_detaold_delta.dat","w");
       for(a=0;a<(i_ncat-1);a++)
	 fprintf(unit,"%lf %lf\n",o_grad[a],s_delta[a]);
       fclose(unit);
       #endif       

       sum = G_ZERO;
       for(a=0;a<(i_ncat-1);a++)
	 {
	   x_alpha_opt[a] = x_alpha_opt[a]+s_delta[a];
	   sum += x_alpha_opt[a];
	 }
       x_alpha_opt[i_ncat-1] = -sum;

       err = age_haul_calc_prob(i_ncat,x_alpha_opt,o_prob);

       log_opt = age_haul_calc_posterior(i_ncat,x_alpha_opt,o_prob,i_Ages,i_mu,i_tau);
       scale = G_ONE;
       // If new value do not give a better value, use a smaller Newton-Raphson step
       while(log_opt < (log_opt_old-0.00000001))
	 {
	   scale /= G_TWO;
	   sum = G_ZERO;
	   for(a=0;a<(i_ncat-1);a++)
	     {
	       x_alpha_opt[a] = x_alpha_opt[a]-scale*s_delta[a];
	       sum += x_alpha_opt[a];
	     }
	   x_alpha_opt[i_ncat-1] = -sum;

	   err = age_haul_calc_prob(i_ncat,x_alpha_opt,o_prob);

	   log_opt = age_haul_calc_posterior(i_ncat,x_alpha_opt,o_prob,i_Ages,i_mu,i_tau);
	 }
       more = fabs(log_opt_old-log_opt) > 0.001;
       log_opt_old = log_opt;
       it++;
    }

  *o_log_opt = log_opt;

  return(0);
}		/* end of age_haul_find_mode */


/*!
  \author Geir Storvik
  \brief Calculates posterior for alpha given observed ages and prior expectations
*/
static double age_haul_calc_posterior(int i_ncat,double *i_alpha,double *i_prob,
				      double *i_Ages,double *i_mu,double i_tau)
{
  int       a;
  double    log_pr,res,log_d;

  /* Prior */
  log_pr = G_ZERO;
  for(a=0;a<i_ncat;a++)
    {
      res = i_alpha[a] - i_mu[a];
      log_pr -= res*res;
    }
  log_pr *= G_HALF*i_tau;
  /* Likelihood */
  log_d = G_ZERO;
  for(a=0;a<i_ncat;a++)
    log_d += log(i_prob[a]) * i_Ages[a];

  return(log_pr+log_d);
}		/* end of age_haul_calc_posterior */


/*!
  \author Geir Storvik
  \brief Calculates probabilities for age-groups given alpha's
*/
static int age_haul_calc_prob(int i_ncat,double *i_alpha,double *o_prob)
{
  int    a;
  double sum;
  
  sum = G_ZERO;
  for(a=0;a<i_ncat;a++)
    {
      o_prob[a] = exp(i_alpha[a]);
      sum += o_prob[a];
    }
  for(a=0;a<i_ncat;a++)
    {
      o_prob[a] /= sum;
    }
  
  return(0);
}


/*!
  \author Geir Storvik
  \brief Calculates loglikelihood of poisson distribution for age-data

  This routine is used in sample_age_ran where the multinomial distribution 
  is extended by a intensity variable in order to make the age-data
  poisson distributed.
*/
static int Loglik_poisson_age(double *logll,double *x,
		       int m,int idx,double *x_vec,char *arg)
{
  int i;
  double *N, *mu;
  char **args;

  args = (char **)arg;
  N = (double *)args[0];
  mu = (double *)args[1];

  for(i=0;i<m;i++) logll[i] = N[idx]*x[i]-mu[idx]*exp(x[i]);

  return GMRFLib_SUCCESS;
}		/* end of Loglik_poisson_age */



/*!
  \author Geir Storvik
  \brief Samples precision parameters for haul effect in age model
  080401:  Modified in order to take into account that sum_a alpha_{h,a}=0
*/
int sample_precision_age_haul(int i_start_h,Eff_str *i_par,
                              double **i_alpha,Data_glm *i_glm,int i_nHaul)
{
  int    n,a,h;
  double mu,ssq;

  i_glm->xcov[0]->n_cov--; //Not haul effect
  n = (i_nHaul-i_start_h)*(i_glm->ncat-1);
  #ifdef DEBUG_PROG
  fprintf(stderr,"caa_sample_multi.c: sample_precision_age_haul: Use n = (i_nHaul-i_start_h-1)*(i_glm->ncat-1)?\n");
  //n = (i_nHaul-i_start_h-1)*(i_glm->ncat-1);
  #endif
  ssq = 0;
  for(h=i_start_h;h<i_nHaul;h++)
    {
      for(a=0;a<i_glm->ncat;a++)
	{
	  mu = calc_eff_suff_age(i_glm->xcov[0],i_par,h,a);
	  //	    mu += calc_eff(i_glm->xcov[i],i_par->eff[a][i],h)*i_glm->suff[h][0][i];
	  ssq += pow(i_alpha[h][a] - mu,2);
	}
      if(!(ssq >= 0 && ssq < 9999999999999.99))
	{
	  fprintf(stderr,"ssq=%lf\n",ssq);
	  write_warning("sample_precision_age_haul:Something is wrong\n");
	  return(1);
	}
    }
  if(!(ssq > 0 && ssq < 9999999999999.99))
    {
      fprintf(stderr,"ssq=%lg\n",ssq);
      write_warning("sample_precision_age_haul:Something is wrong\n");
      return(1);
    }
  i_par->tau_obs = gengam(i_par->prior_prec_obs[0]+G_HALF * ssq,
			  i_par->prior_prec_obs[1]+G_HALF * (double) n);

  i_glm->xcov[0]->n_cov++;

  return(0);
}		/* end of sample_precision_age_haul */


