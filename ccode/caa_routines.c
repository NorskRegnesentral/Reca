/*!
  \file caa_routines.c
  \brief Routines for converting data into approperiate structures and some other general
  routines (including writing input data to files)
*/
#include "caa.h"
#include "caa_routines.h"
#include "caa_read_write.h"
#include "caa_utl.h"
#include "caa_util.h"
#include "caa_sample_g_a.h"
//#include "caa_COST.h"

static int make_c_cov(int i_n_cov,int *i_n_lev,int i_ispat,int i_iboat,int i_ihaulsize,int *i_fix,int *i_interaction,
		      int *i_in_landings,int **i_x_cov,int i_nHaul,Data_cov **o_xcov,int i_fit,int i_incl_haul);
static int re_make_c_cov(int i_nHaul,Data_cov **o_xcov);
static int make_spat_struct(int *i_num,int *i_adj_area,Data_cov *x_xcov);
static int re_make_spat_struct(int *i_num,int *i_adj_area,Data_cov *x_xcov);
static int convert_cov(int i_nHaul,Data_cov *i_xcov);
static int re_convert_cov(Data_cov *i_xcov);
static int convert_tot_cov(Data_totcatch *i_D_totcatch,int *i_fac,Data_cov *i_xcov,
			   Data_cov *o_xcov,int i_neff);
static int re_convert_tot_cov(Data_totcatch *i_D_totcatch,int *i_fac,Data_cov *i_xcov,
			      Data_cov *o_xcov,int i_neff);
int compare(int *i_x,int *i_y);
static int update_average_par(int i_n,int i_ncat,Eff_str *i_par,Eff_str *x_par_mean,
                              int i_nxcov,Data_cov **i_xcov);

#ifdef LOG_FILE
extern FILE     *g_caa_log;
#endif

extern FILE     *g_caa_mcmc1;
extern FILE     *g_caa_mcmc2;

extern FILE     *g_caa_mcmc_hsz_eff;


/*!
  \author Geir Storvik
  \brief Makes a struct of type containing model and covariates for age model

  Makes a struct of type Data_age 
  describing the model structure, covariates and number of hauls
  with and without age-data.

  Space allocated in this routine is reallocated in ::re_makedata_age1  
*/
int makedata_age1(int i_nHaul,int i_nAges,int *i_a_vec,
		  Input_cov *i_cov, int *i_num_adj_area,int *i_adj_area,
		  Data_age **o_D_age,int i_fit)
{
  int         a,i,ind,ncov,err;
  int         ispat,iboat,ihaul,ihaulsize,icell,n_fac_cell,n_ispat,n_iboat;
  int        *fixed,*interaction,*in_landings,*n_fac;
  Data_cov   *xcov;
  Data_age   *D_age;


  D_age = CALLOC(1,Data_age);       // Free ok
  D_age->glm = CALLOC(1,Data_glm);  // Free ok


  /* Sizes */
  D_age->glm->nHaul = i_nHaul;
  D_age->glm->ncat=i_nAges;
  ncov = i_cov->n_cov_i+i_cov->n_cov_d+1; // Include haul as covariate here

  D_age->a_vec = CALLOC(D_age->glm->ncat,int);  // Free ok
  for(a=0;a<D_age->glm->ncat;a++)
    {
      D_age->a_vec[a] = i_a_vec[a];
    }

  /* Covariates */
  D_age->glm->nxcov = 1;  /* Only Intercept */
  D_age->glm->xcov = CALLOC(D_age->glm->nxcov,Data_cov *);   // Free ok

  // create new vectors - have different length than original because ncov is larger
  fixed = CALLOC(ncov,int);  // Free ok
  interaction = CALLOC(ncov,int);  // Free ok
  in_landings = CALLOC(ncov,int);  // Free ok
  n_fac = CALLOC(ncov,int);  // Free ok
  n_fac_cell = 1;
  ispat = -1;
  iboat = -1;
  icell = i_cov->icell; // use input from R
  n_ispat = 0;
  n_iboat = 0;
  ihaulsize = i_cov->ihaulsize;
  ind = 0;
  for(i=0;i<i_cov->n_cov;i++)
    {
      #ifdef DEBUG_PROG
      fprintf(stderr,"i=%d,nlev=%d,random=%d,interaction=%d,in.landings=%d,CAR=%d,continuous=%d\n",i,i_cov->n_lev[i],i_cov->random[i],i_cov->interaction[i],i_cov->in_landings[i],i_cov->spatial[i],i_cov->continuous[i]);
      #endif
      interaction[i] = i_cov->interaction[i];
      in_landings[i] = i_cov->in_landings[i];
      if(i_cov->random[i]==1)
	fixed[ind] = 0;
      else
	fixed[ind] = 1;
      n_fac[ind] = i_cov->n_lev[i];
      ind++;
      if(i_cov->interaction[i]==1) // Use input from R instead?
	n_fac_cell = n_fac_cell*(i_cov->n_lev[i]-1);
      if(i_cov->spatial[i]==1)
	{
	  ispat = i;
	  n_ispat++;
	}
      if((i_cov->in_landings[i]==0)&&(i_cov->random[i]==1)&&(i!=icell))
	{
	  iboat = i;
	  n_iboat++;
	}
    }
  fixed[ncov-1] = 0;
  interaction[ncov-1] = 0;
  in_landings[ncov-1] = 0;
  n_fac[ncov-1] = D_age->glm->nHaul;
  ihaul = ncov-1;
  if(n_ispat>1)
    {
      write_warning("makedata_age1:Number of CAR variables can only be 0 or 1\n");
      return(1);
    }
  if(n_iboat>1)
    {
      write_warning("makedata_age1:Number of random variables not in landings can only be 0 or 1\n");
      return(1);
    }

  #ifdef DEBUG_PROG
  fprintf(stderr,"ncov=%d,ispat=%d,icell=%d,iboat=%d,ihaulsize=%d,ihaul=%d,(nfac_cell-1)=%d\n",i_cov->n_cov,ispat,icell,iboat,ihaulsize,ihaul,n_fac_cell);
  #endif
  
  // Intercept
  err = make_c_cov(ncov,n_fac,ispat,iboat,ihaulsize,fixed,interaction,in_landings,i_cov->c_cov_i,D_age->glm->nHaul,&xcov,i_fit,1);
  if(err)
    {
      write_warning("makedata_age1:Error calling make_c_cov\n");
      return(err);
    }
  D_age->glm->xcov[0] = xcov;
  D_age->glm->xcov[0]->icell = icell;
  D_age->glm->xcov[0]->iboat = iboat;
  if(i_cov->icell)
    D_age->glm->xcov[0]->n_fac_cell = n_fac_cell;
  D_age->glm->xcov[0]->ihaul = ihaul;
  //D_age->glm->xcov[0]->ihaulsize = i_int_ihaulsize-1;    /* index for haulsize boat term */
  D_age->glm->xcov[0]->ihaulsize = ihaulsize;    /* index for haulsize boat term */

  if(D_age->glm->xcov[0]->ihaulsize>0)
    D_age->glm->inc_hsz = 1;
  else 
    D_age->glm->inc_hsz = 0;

  #ifdef DEBUG_PROG
  fprintf(stderr,"ncov=%d\n",D_age->glm->xcov[0]->n_cov);
  for(i=0;i<D_age->glm->xcov[0]->n_cov;i++){
    fprintf(stderr,"nfac=%d\n",D_age->glm->xcov[0]->n_fac[i]);
  }
  fprintf(stderr,"ihaul=%d,ispat=%d,iboat=%d,icell=%d\n",D_age->glm->xcov[0]->ihaul,
	  D_age->glm->xcov[0]->ispat,D_age->glm->xcov[0]->iboat,D_age->glm->xcov[0]->icell);
  #endif
  
  
  /* Convert covariates */
  if(i_fit) //not needed in predict
    {
      err = convert_cov(D_age->glm->nHaul,D_age->glm->xcov[0]);
      if(err)
	{
	  write_warning("makedata_age1:Error calling convert_cov\n");
	  return(err);
	}
    }
  
  /* spatial structure */
  if(i_fit) //not needed in predict
    {
      if(D_age->glm->xcov[0]->ispat > -1)
	make_spat_struct(i_num_adj_area,i_adj_area,D_age->glm->xcov[0]);
    }
  
  FREE(interaction);
  FREE(in_landings);
  FREE(fixed);
  FREE(n_fac);

  *o_D_age = D_age;

  return(0);
}		/* end of makedata_age1 */

    

/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in ::makedata_age1
*/
int re_makedata_age1(int *i_num_adj_area,int *i_adj_area,Data_age **o_D_age,int i_fit)
{
  int         err;
  Data_cov   *xcov;
  Data_age   *D_age;

  D_age = *o_D_age;


  xcov = D_age->glm->xcov[0];

  if(i_fit)
    {
      if(D_age->glm->xcov[0]->ispat > -1)
	re_make_spat_struct(i_num_adj_area,i_adj_area,xcov);
      err = re_convert_cov(D_age->glm->xcov[0]);
      if(err)
	{
	  write_warning("re_makedata_age1:Error calling re_convert_cov\n");
	  return(err);
	}
    }

  re_make_c_cov(D_age->glm->nHaul,&xcov);

  FREE(D_age->glm->xcov);
  FREE(D_age->glm);
  FREE(D_age->a_vec);
  FREE(D_age);

  return(0);
}		/* end of re_makedata_age1 */

    

/*!
  \author Geir Storvik
  \brief Put observed ages for different hauls into D_age->Ages for Amigo data.

  Note that the first hauls are for length-only data
  Memory allocated in this routine is reallocated in ::re_makedata_age2
*/
int makedata_age2(Data_orig *i_D_orig,Data_age *i_D_age,int class_error)
{
  int         a,a2,f,h,n,ind;
  int         n_type_miss;

  /* Ages */
  i_D_age->n_h = CALLOC(i_D_age->glm->nHaul,int);    // Free ok
 
  i_D_age->Ages = CALLOC2_d(i_D_age->glm->nHaul,i_D_age->glm->ncat);

  i_D_age->Ages_fix = CALLOC2_d(i_D_age->glm->nHaul,i_D_age->glm->ncat);

  i_D_age->Ages_disc = Mmatrix_2d(0,i_D_age->glm->nHaul-1,0,i_D_age->glm->ncat,
				   sizeof(double),1); //Free ok
  i_D_age->Ages_land = Mmatrix_2d(0,i_D_age->glm->nHaul-1,0,i_D_age->glm->ncat,
				   sizeof(double),1); //Free ok


  if(i_D_orig->coastal_cod) //difference between coastal cod and skrei
    {
      a2 = (int) i_D_age->glm->ncat/2;
      ind = 0;
      n_type_miss = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_D_age->Ages[h][a] = G_ZERO;
	  for(f=0;f<i_D_orig->nFishBoat[h];f++)
	    {
	      if(i_D_orig->totage[ind] > -1000)
		{
		  a = i_D_orig->totage[ind]-i_D_age->a_vec[0];
		  if(a < 0 || a >= i_D_age->a_vec[i_D_age->glm->ncat-1])
		    {
		      fprintf(stderr,"h=%d,f=%d,ind=%d,totage=%d,a_vec0=%d,a_vecA=%d\n",
			      h,f,ind,i_D_orig->totage[ind],i_D_age->a_vec[0],i_D_age->a_vec[i_D_age->glm->ncat-1]);
		      write_warning("makedata_age2:Something is wrong\n");
		      return(1);
		    }
		  /* Use certain observations for starting values, even when classification error */
		  if(i_D_orig->tottype[ind] == 1) //certain coastal cod
		    i_D_age->Ages[h][a2+a] += (double) i_D_orig->replength[ind];
		  else if(i_D_orig->tottype[ind] == 2 && class_error == 0)
		    i_D_age->Ages[h][a2+a] += (double) i_D_orig->replength[ind];
		  else if(i_D_orig->tottype[ind] == 4 && class_error == 0)
		    i_D_age->Ages[h][a] += (double) i_D_orig->replength[ind];
		  else if(i_D_orig->tottype[ind] == 5) //certain skrei
		    i_D_age->Ages[h][a] += (double) i_D_orig->replength[ind];
		  else if(i_D_orig->tottype[ind]<0)
		    n_type_miss++;
		  else if((i_D_orig->tottype[ind] !=2) && (i_D_orig->tottype[ind] !=4))
		    fprintf(stderr,"makedata_age2: WARNING: Wrong type included: h=%d,f=%d,type=%d\n",
			    h,f,i_D_orig->tottype[ind]); 
		}
	      ind++;
	    }
	  n = 0;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    n += (int) i_D_age->Ages[h][a];
	  i_D_age->n_h[h] = n;
	}
      #ifdef DEBUG_PROG
      printf("%d fish with age observed and type missing\n",n_type_miss);
      #endif
    }
  else  //original version
    {
      ind = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	{
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    i_D_age->Ages[h][a] = G_ZERO;
	  for(f=0;f<i_D_orig->nFishBoat[h];f++)
	    {
	      if(i_D_orig->totage[ind] > -1000)
		{
		  a = i_D_orig->totage[ind]-i_D_age->a_vec[0];
		  if(a < 0 || a >= i_D_age->glm->ncat)
		    {
		      fprintf(stderr,"h=%d,f=%d,ind=%d,totage=%d,a_vec0=%d\n",
			      h,f,ind,i_D_orig->totage[ind],i_D_age->a_vec[0]);
		      write_warning("makedata_age2:Something is wrong\n");
		      return(1);
		    }
		  i_D_age->Ages[h][i_D_orig->totage[ind]-i_D_age->a_vec[0]] += (double) i_D_orig->replength[ind];
		}
	      ind++;
	    }
	  n = 0;
	  for(a=0;a<i_D_age->glm->ncat;a++)
	    n += (int) i_D_age->Ages[h][a];
	  i_D_age->n_h[h] = n;
	}
    } 
  for(h=0;h<i_D_age->glm->nHaul;h++)
    for(a=0;a<i_D_age->glm->ncat;a++)
      {
	i_D_age->Ages_fix[h][a] = i_D_age->Ages[h][a];
      }

  int sum,suma;
  #ifdef DEBUG_PROG
  sum = 0;
  fprintf(stderr,"makedata_age2: Num_age= ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      suma = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	suma += (int) i_D_age->Ages[h][a];
      fprintf(stderr,"%d ",suma);
      sum += suma;
    }
  fprintf(stderr,"\n");
  fprintf(stderr,"Total aged fish=%d\n",sum);
  #endif
  #ifdef LOG_FILE
  sum = 0;
  fprintf(g_caa_log,"makedata_age2: Num_age= ");
  for(a=0;a<i_D_age->glm->ncat;a++)
    {
      suma = 0;
      for(h=0;h<i_D_age->glm->nHaul;h++)
	suma += (int) i_D_age->Ages[h][a];
      fprintf(g_caa_log,"%d ",suma);
      sum += suma;
    }
  fprintf(g_caa_log,"\n");
  fprintf(g_caa_log,"Total aged fish=%d\n",sum);
  #endif

  return(0);
}		/* end of makedata_age2 */

    

/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in ::makedata_age2
*/
int re_makedata_age2(Data_age *i_D_age)
{
  /* Ages */
  FREE(i_D_age->n_h);
  FREE2_d(i_D_age->Ages,i_D_age->glm->nHaul);
  FREE2_d(i_D_age->Ages_fix,i_D_age->glm->nHaul);
  // Fmatrix_2d(&i_D_age->Ages[0][0],&i_D_age->Ages[0]);
  Fmatrix_2d(&i_D_age->Ages_disc[0][0],&i_D_age->Ages_disc[0]);
  Fmatrix_2d(&i_D_age->Ages_land[0][0],&i_D_age->Ages_land[0]);


  return(0);
}		/* end of re_makedata_age2 */



/*!
  \author Hanne Rognebakke
  \brief Allocate space and initialize struct for data in g(a)-function

  Memory allocated in this routine is reallocated in ::re_makedata_g_a
*/
int makedata_g_a(int i_age_ncat, int i_g_a_ncat, int i_g_a_nSeason, int *i_a_vec, int i_coastal_cod,
		 int i_g_a_model, double *i_par_init, int i_sample_c, int i_sample_theta,int i_sample_gamma,
		 Data_g_a **o_D_g_a)
{
  int        a,i,A,err;
  Data_g_a  *D_g_a;
 

  /* Allocating space */
  D_g_a = CALLOC(1,Data_g_a);                  // Free ok
  D_g_a->ncat = i_g_a_ncat;
  D_g_a->nSeason = i_g_a_nSeason;
  D_g_a->a_vec = CALLOC(D_g_a->ncat,double);  // Free ok
  D_g_a->a_vec[0] = i_a_vec[0]+(double) 1/D_g_a->nSeason;
  //printf("ga avec[0]=%f\n",D_g_a->a_vec[0]);
  for(a=1;a<D_g_a->ncat;a++)
    {
      D_g_a->a_vec[a] = D_g_a->a_vec[a-1]+ (double) 1/D_g_a->nSeason;
      //printf("ga avec[%d]=%f\n",a,D_g_a->a_vec[a]);
    }

  D_g_a->g_a = CALLOC(D_g_a->ncat,double);   // Free ok
  D_g_a->g_a_mean = CALLOC(D_g_a->ncat,double);   // Free ok

  if(i_coastal_cod)
    {
      A = (int)i_age_ncat/2;
      D_g_a->a2Age_vec = CALLOC(i_age_ncat,int);  // Free ok
      D_g_a->a2Age_vec[0] = 0;
      for(a=1;a<A;a++)
	{ 
	  D_g_a->a2Age_vec[a] = D_g_a->a2Age_vec[a-1]+D_g_a->nSeason;
	}
      D_g_a->a2Age_vec[A] = 0;
      for(a=A+1;a<i_age_ncat;a++)
	{ 
	  D_g_a->a2Age_vec[a] = D_g_a->a2Age_vec[a-1]+D_g_a->nSeason;
	}
      //for(a=0;a<i_age_ncat;a++)
      //printf("a2Age_vec[%d]=%d\n",a,D_g_a->a2Age_vec[a]);
    }
  else
    {
      D_g_a->a2Age_vec = CALLOC(i_age_ncat,int);  // Free ok
      D_g_a->a2Age_vec[0] = 0;
      //printf("a2Age_vec[0]=%d\n",D_g_a->a2Age_vec[0]);
      for(a=1;a<i_age_ncat;a++)
	{ 
	  D_g_a->a2Age_vec[a] = D_g_a->a2Age_vec[a-1]+D_g_a->nSeason;
	  //printf("a2Age_vec[%d]=%d\n",a,D_g_a->a2Age_vec[a]);
	}
    }

  D_g_a->g_a_model = i_g_a_model;
  if(D_g_a->g_a_model==0) /* log-linear model */
    {
      D_g_a->g_a_npar = 0;
      err = calc_g_a_log_lin(D_g_a->ncat,D_g_a->a_vec,D_g_a->g_a_par,D_g_a->g_a);
      if(err)
	{
	  write_warning("makedata_g_a:Error calling calc_g_a_log_lin\n");
	  return(err);
	}
    }
  else if(D_g_a->g_a_model==1) /* Schnute-Richards model */
    {
      D_g_a->g_a_npar = 3;
      D_g_a->g_a_par = CALLOC(3,double);            // Free ok
      D_g_a->g_a_par_mean = CALLOC(3,double);            // Free ok
      for(i=0;i<D_g_a->g_a_npar;i++)
	{
	  D_g_a->g_a_par[i] = i_par_init[i];
	  //printf("g_a_par[%d]=%f\n",i,D_g_a->g_a_par[i]);
	}
      D_g_a->sample_c = i_sample_c;
      D_g_a->sample_theta = i_sample_theta;
      D_g_a->sample_gamma = i_sample_gamma;
      err = calc_g_a_S_R2(D_g_a->ncat,D_g_a->a_vec,D_g_a->g_a_par,D_g_a->g_a);
      if(err)
	{
	  write_warning("makedata_g_a:Error calling calc_g_a_S_R\n");
	  return(err);
	}
    }
  else 
    {
      write_warning("makedata_g_a:Unknown g_a_model\n");
      return(1);
    }
  //fprintf(stderr,"makedata_g_a:\n");
  //for(a=0;a<D_g_a->ncat;a++)
  //fprintf(stderr,"g_a[%d]=%f\n",a,D_g_a->g_a[a]); 

  /* Allocating space for sufficient statistics for g(a) */
  D_g_a->suff = Mmatrix_2d(0,2,0,D_g_a->ncat-1,sizeof(double),1);  // Free ok

  *o_D_g_a = D_g_a;

  return(0);
}		/* end of makedata_g_a */



/*!
  \author Hanne Rognebakke
  \brief Re-allocating space allocated by ::makedata_g_a
*/
int re_makedata_g_a(Data_g_a **o_D_g_a)
{
  Data_g_a  *D_g_a;

  D_g_a = *o_D_g_a;
  FREE(D_g_a->a_vec);
  FREE(D_g_a->a2Age_vec);
  FREE(D_g_a->g_a);
  FREE(D_g_a->g_a_mean);
  if(D_g_a->g_a_model>0)
    {
      FREE(D_g_a->g_a_par);
      FREE(D_g_a->g_a_par_mean);
    }
  Fmatrix_2d(&D_g_a->suff[0][0],&D_g_a->suff[0]);

  FREE(D_g_a);

  return(0);
}		/* end of re_makedata_g_a */



/*!
  \author Geir Storvik
  \brief Construct a struct describing the model structure and covariates for linear models

  With and without covariates.
  Note that effects are transformed to start at zero.

  Memory allocated in this routines is reallocated in ::re_makedata_lin1
*/
int makedata_lin1(int i_nHaul,double *i_haulweight,
		  Input_cov *i_int_cov,Input_cov *i_slp_cov,
		  int *i_num_adj_area,int *i_adj_area,
		  Data_lin **o_D_lin,int i_fit,int i_inc_hsz)
{
  int        i,ind,ncov,err;
  int        ispat,iboat,ihaul,icell,ihaulsize,n_fac_cell;
  int       *fixed,*interaction,*in_landings,*n_fac;
  Data_cov  *xcov;
  Data_lin  *D_lin;


  D_lin = CALLOC(1,Data_lin);         // Free ok
  D_lin->glm = CALLOC(1,Data_glm);    // Free ok
  D_lin->glm->ncat = 1;

  /* Sizes */
  D_lin->glm->nHaul = i_nHaul;
  D_lin->haulweight = i_haulweight;

 
  /* Covariates */
  D_lin->glm->nxcov = 2;  /* Intercept and slope */
  D_lin->glm->xcov = CALLOC(D_lin->glm->nxcov,Data_cov *);   // Free ok
  D_lin->glm->inc_hsz = i_inc_hsz;

  /* Intercept */
  if(i_inc_hsz==0){
    ncov = i_int_cov->n_cov_i+1; // Include haul as covariate here for lga and wgl mod
  } else {
    ncov = i_int_cov->n_cov_i;
  }
  fixed = CALLOC(ncov,int);  // Free ok
  n_fac = CALLOC(ncov,int);  // Free ok
  interaction = CALLOC(ncov,int); //Free ok  - Must extend the current pointer since haul can be included as covariate
  in_landings = CALLOC(ncov,int); //Free ok  - Must extend the current pointer since haul can be included as covariate
  n_fac_cell = 1;
  ispat = -1;
  iboat = -1;
  icell = i_int_cov->icell;
  ind = 0;
  for(i=0;i<i_int_cov->n_cov;i++)
    {
      #ifdef DEBUG_PROG
      fprintf(stderr,"i=%d,nlev=%d,random=%d,interaction=%d,in.landings=%d,CAR=%d,continuous=%d\n",i,i_int_cov->n_lev[i],i_int_cov->random[i],i_int_cov->interaction[i],i_int_cov->in_landings[i],i_int_cov->spatial[i],i_int_cov->continuous[i]);
      #endif
      interaction[i] = i_int_cov->interaction[i];
      in_landings[i] = i_int_cov->in_landings[i];
      if(i_int_cov->random[i] == 1)
	fixed[ind] = 0;
      else
	fixed[ind] = 1;
      n_fac[ind] = i_int_cov->n_lev[i];
      ind++;
      if(i_int_cov->interaction[i]==1)
	n_fac_cell = n_fac_cell*(i_int_cov->n_lev[i]-1);
      if(i_int_cov->spatial[i]==1)
	ispat = i;
      if((i_int_cov->in_landings[i]==0)&&(i_int_cov->random[i]==1)&&(i!=icell))
	iboat = i;
    }
  if(i_inc_hsz==0){
    fixed[ncov-1] = 0;
    n_fac[ncov-1] = D_lin->glm->nHaul;
    interaction[ncov-1] = 0;
    in_landings[ncov-1] = 0;
    ihaul = ncov-1;
  } else {
    ihaul = -1; // Haul not included in haulsize model
  }  
  ihaulsize = -1;  // Not included in lga or wgl model

  #ifdef DEBUG_PROG
  fprintf(stderr,"makedata_lin1: Include haul covariate as an option??\n");
  fprintf(stderr,"ncov=%d,ispat=%d,icell=%d,iboat=%d,ihaulsize=%d,ihaul=%d,(nfac_cell-1)=%d\n",i_int_cov->n_cov,ispat,icell,iboat,ihaulsize,ihaul,n_fac_cell);
  #endif

  if(i_inc_hsz==0)
    err = make_c_cov(ncov,n_fac,ispat,iboat,ihaulsize,fixed,interaction,in_landings,i_int_cov->c_cov_i,D_lin->glm->nHaul,&xcov,i_fit,1);
  else
    err = make_c_cov(ncov,n_fac,ispat,iboat,ihaulsize,fixed,interaction,in_landings,i_int_cov->c_cov_i,D_lin->glm->nHaul,&xcov,i_fit,0);
  if(err)
    {
      write_warning("makedata_lin1:Error calling make_c_cov for intercept\n");
      return(err);
    }
  D_lin->glm->xcov[0] = xcov;
  D_lin->glm->xcov[0]->ihaul = ihaul;
  D_lin->glm->xcov[0]->ihaulsize = ihaulsize;
  D_lin->glm->xcov[0]->icell = i_int_cov->icell;
  D_lin->glm->xcov[0]->iboat = iboat;
  if(i_int_cov->icell)
    D_lin->glm->xcov[0]->n_fac_cell = n_fac_cell;
  FREE(fixed);
  FREE(interaction);
  FREE(in_landings);
  FREE(n_fac);
  
  /* Convert covariates */
  if(i_fit) //not needed in predict
    {
      err = convert_cov(D_lin->glm->nHaul,D_lin->glm->xcov[0]);
      if(err)
	{
	  write_warning("makedata_lin1:Error calling convert_cov for intercept\n");
	  return(err);
	}
    }

  /* Make spatial structure */      
  if(i_fit) //not needed in predict
    {
      if(D_lin->glm->xcov[0]->ispat > -1)
	{
	  err = make_spat_struct(i_num_adj_area,i_adj_area,D_lin->glm->xcov[0]);
	  if(err)
	    {
	      write_warning("makedata_lin1:Error calling make_spat_struct\n");
	      return(err);
	    }
	}
    }

  /* Slope */
  fixed = CALLOC(i_slp_cov->n_cov,int);  // Free ok
  n_fac = CALLOC(i_slp_cov->n_cov,int);  // Free ok
  ispat = -1;
  iboat = -1;
  ihaul = -1;
  icell = -1;
  ihaulsize = -1;
  for(i=0;i<i_slp_cov->n_cov;i++)
    {
      if(i_slp_cov->random[i] == 1)
	fixed[i] = 0;
      else
	fixed[i] = 1;
      if(i_slp_cov->spatial[i] == 1)
	ispat = i;
      n_fac[i] = 1;
    }
  err = make_c_cov(i_slp_cov->n_cov,n_fac,ispat,iboat,ihaulsize,fixed,i_slp_cov->interaction,i_slp_cov->in_landings,i_slp_cov->c_cov_i,D_lin->glm->nHaul,&xcov,i_fit,0);
  if(err)
    {
      write_warning("makedata_lin1:Error calling make_c_cov\n");
      return(err);
    }
  D_lin->glm->xcov[1] = xcov;
  D_lin->glm->xcov[1]->ihaul = ihaul;
  D_lin->glm->xcov[1]->icell = icell;
  D_lin->glm->xcov[1]->ihaulsize = ihaulsize;
  FREE(fixed);
  FREE(n_fac);

  /* Convert covariates */
  if(i_fit) //not needed in predict
    {
      err = convert_cov(D_lin->glm->nHaul,D_lin->glm->xcov[1]);
      if(err)
	{
	  write_warning("makedata_lin1:Error calling convert_cov for slope\n");
	  return(err);
	}
    }

  /* Make spatial structure */
  if(i_fit)
    {
      if(D_lin->glm->xcov[1]->ispat > -1)
	{
	  err = make_spat_struct(i_num_adj_area,i_adj_area,D_lin->glm->xcov[1]);
	  if(err)
	    {
	      write_warning("makedata_lin1:Error calling make_spat_struct\n");
	      return(err);
	    }
	}
    }


  *o_D_lin = D_lin;

  return(0);
}		/* end of makedata_lin1 */



/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in ::makedata_lin1
*/
int re_makedata_lin1(int *i_num_adj_area,int *i_adj_area,
		     Data_lin **o_D_lin,int i_fit)
{
  int        err;
  Data_cov  *xcov;
  Data_lin  *D_lin;

  D_lin = *o_D_lin;

  /* Int */
  xcov = D_lin->glm->xcov[0];
  if(i_fit)
    {
      if(xcov->ispat > -1)
	{
	  err = re_make_spat_struct(i_num_adj_area,i_adj_area,xcov);
	  if(err)
	    {
	      write_warning("re_makedata_lin1:Error calling re_make_spat_struct\n");
	      return(err);
	    }
	}
      err = re_convert_cov(xcov);
      if(err)
	{
	  write_warning("re_makedata_lin1:Error calling re_convert_cov\n");
	  return(err);
	}
    }

  err = re_make_c_cov(D_lin->glm->nHaul,&xcov);
  if(err)
    {
      write_warning("re_makedata_lin1:Error calling re_make_c_cov\n");
      return(err);
    }


  /* Slope */
  xcov = D_lin->glm->xcov[1];
  if(i_fit)
    {
      if(xcov->ispat > -1)
	{
	  err = re_make_spat_struct(i_num_adj_area,i_adj_area,xcov);
	  if(err)
	    {
	      write_warning("re_makedata_lin1:Error calling re_ make_spat_struct\n");
	      return(err);
	    }
	}    
      err = re_convert_cov(xcov);
      if(err)
	{
	  write_warning("re_makedata_lin1:Error calling re_convert_cov\n");
	  return(err);
	}
    }
  err = re_make_c_cov(D_lin->glm->nHaul,&xcov);
  if(err)
    {
      write_warning("re_makedata_lin1:Error calling re_make_c_cov\n");
      return(err);
    }


  FREE(D_lin->glm->xcov);
  FREE(D_lin->glm);
  FREE(D_lin);

  return(0);
}		/* end of re_makedata_lin1 */



/*!
  \author Geir Storvik
  \brief Calculates sufficient statistics for the lga model 

  Use amigo data and age-stratified by length data.

  This routine is only used in initialization. 
  It calculates summary statistics based on aged fish.
  Updating sufficient statistics based on length-only data is performed
  by the ::make_suff_lga routine.

  \todo The routine copy summaries to the fix-versions. This should be made
  more clean later on.

  Memory allocated in this routine is reallocated in ::re_makedata_lga_suff.
*/
int makedata_lga_suff(Data_lin *i_D_lga,Data_orig *i_D_orig,Data_g_a *i_D_g_a)
{
  int        a,b,f,h,season,n,ncat,N,ind_f,ind_a,err=0;
  double    *x,*y,beta0,beta1,ssq,length,r;


  ncat = i_D_g_a->ncat;

  /* Sufficient statistics, observed age-length data + simulated missing ages (added in sample_ages) */
  i_D_lga->Ages = Mmatrix_2d(0,i_D_lga->glm->nHaul,0,ncat-1,sizeof(int),1); // Free ok
  i_D_lga->glm->haul_obs = CALLOC(i_D_lga->glm->nHaul,int);   // Free ok
  i_D_lga->glm->boat_obs = CALLOC(i_D_orig->nBoat,int);   // Free ok

  i_D_lga->sum_by_cat = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,
				   0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->sqsum_by_cat = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,
				     0,ncat-1,sizeof(double),1);    // Free ok
  /* Need size of beta_hat equal to nxcov when finding b-vector in sample_gauss_eff */
  i_D_lga->glm->beta_hat = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,0,
				      0,i_D_lga->glm->nxcov-1,sizeof(double),1);    // Free ok
  i_D_lga->glm->ssq = CALLOC(i_D_lga->glm->nHaul,double);                 // Free ok
  i_D_lga->glm->suff = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,i_D_lga->glm->nxcov-1,
				  0,i_D_lga->glm->nxcov-1,sizeof(double),1); // Free ok

  /* Sufficient statistics based on observed age-length, constant during the simulations */
  i_D_lga->Ages_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul,0,ncat-1,sizeof(int),1); // Free ok
  i_D_lga->sum_by_cat_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,
				       0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->sqsum_by_cat_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,
					 0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->glm->beta_hat_fix = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,0,
				      0,i_D_lga->glm->nxcov-1,sizeof(double),1);    // Free ok
  i_D_lga->glm->ssq_fix = CALLOC(i_D_lga->glm->nHaul,double);                 // Free ok
  i_D_lga->glm->suff_fix = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,i_D_lga->glm->nxcov-1,
				  0,i_D_lga->glm->nxcov-1,sizeof(double),1); // Free ok



  N = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    N = max(N,i_D_orig->nFishBoat[h]);
  x = CALLOC(N,double);                  // Free ok
  y = CALLOC(N,double);                  // Free ok

  for(b=0;b<i_D_orig->nBoat;b++)
    i_D_lga->glm->boat_obs[b] = 0;

  ind_f = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      i_D_lga->glm->haul_obs[h] = 0;
      for(a=0;a<ncat;a++)
	{ 
	  i_D_lga->Ages[h][a] = 0;
	  i_D_lga->sum_by_cat[h][a] = G_ZERO;
	  i_D_lga->sqsum_by_cat[h][a] = G_ZERO;
	}
      i_D_lga->glm->suff[h][0][1] = G_ZERO;
      i_D_lga->glm->suff[h][1][1] = G_ZERO;
      n = 0;
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
          a = i_D_orig->totage[ind_f];
          length = i_D_orig->totlength[ind_f];
	  if(a > -9000 && length > -9000)
	    {
	      i_D_lga->glm->haul_obs[h] = 1;
	      season = i_D_orig->season[ind_f];
	      ind_a = i_D_g_a->a2Age_vec[a-(int) i_D_g_a->a_vec[0]]+(season-1);
              x[n] = i_D_g_a->g_a[ind_a];
              y[n] = length;
	      r = (double) i_D_orig->replength[ind_f];
              i_D_lga->Ages[h][ind_a] += i_D_orig->replength[ind_f];
	      i_D_lga->sum_by_cat[h][ind_a] += r*y[n];
	      i_D_lga->sqsum_by_cat[h][ind_a] += r*y[n]*y[n];
 	      i_D_lga->glm->suff[h][0][0] += r;
	      i_D_lga->glm->suff[h][0][1] += r*x[n];
	      i_D_lga->glm->suff[h][1][1] += r*x[n]*x[n];
	      n++;
	    }
          ind_f ++;
	}
      if(i_D_lga->glm->haul_obs[h]==1)
        {
          b = i_D_orig->boat[h]-1;
          i_D_lga->glm->boat_obs[b] = 1;
        }
      i_D_lga->glm->suff[h][1][0] = i_D_lga->glm->suff[h][0][1];

      if(i_D_lga->glm->suff[h][0][0]>G_ZERO)
	{
	  err = lm_fit_suff(ncat,i_D_g_a->g_a,i_D_lga->sum_by_cat[h],i_D_lga->sqsum_by_cat[h],
                            i_D_lga->glm->suff[h],&beta0,&beta1,&ssq);
	  if(err)
	    {
	      write_warning("makedata_lga_suff:Error calling lm_fit_suff\n");
	      return(err);
	    }
	  i_D_lga->glm->beta_hat[h][0][0] = beta0;
	  i_D_lga->glm->beta_hat[h][0][1] = beta1;
	  i_D_lga->glm->ssq[h] = ssq;
	  //if(h<10)
	  // fprintf(stderr,"makedata_lga_suff: h=%d: beta_0=%f, beta_1=%f, ssq=%f\n",h,beta0,beta1,ssq);
	}

      // Copying summaries to the fix-versions. 
      for(a=0;a<ncat;a++)
	{
	  i_D_lga->Ages_fix[h][a] = i_D_lga->Ages[h][a];
	  i_D_lga->sum_by_cat_fix[h][a] = i_D_lga->sum_by_cat[h][a];
	  i_D_lga->sqsum_by_cat_fix[h][a] = i_D_lga->sqsum_by_cat[h][a];
	}
      i_D_lga->glm->suff_fix[h][0][0] = i_D_lga->glm->suff[h][0][0];
      i_D_lga->glm->suff_fix[h][0][1] = i_D_lga->glm->suff[h][0][1];
      i_D_lga->glm->suff_fix[h][1][0] = i_D_lga->glm->suff[h][1][0];
      i_D_lga->glm->suff_fix[h][1][1] = i_D_lga->glm->suff[h][1][1];
      if(i_D_lga->glm->suff[h][0][0]>G_ZERO)
	{
	  i_D_lga->glm->beta_hat_fix[h][0][0] = beta0;
	  i_D_lga->glm->beta_hat_fix[h][0][1] = beta1;
	  i_D_lga->glm->ssq_fix[h] = ssq;
	}
      else
	{
	  i_D_lga->glm->beta_hat_fix[h][0][0] = G_ZERO;
	  i_D_lga->glm->beta_hat_fix[h][0][1] = G_ZERO;
	  i_D_lga->glm->ssq_fix[h] = G_ZERO;
	}

    } 


  FREE(x);
  FREE(y);

  return(0);
}		/* end of makedata_lga_suff */



/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in ::makedata_lga_suff
*/
int re_makedata_lga_suff(Data_lin *i_D_lga)
{
  Fmatrix_2d(&i_D_lga->Ages[0][0],&i_D_lga->Ages[0]);
  FREE(i_D_lga->glm->haul_obs);
  FREE(i_D_lga->glm->boat_obs);
  Fmatrix_2d(&i_D_lga->sum_by_cat[0][0],&i_D_lga->sum_by_cat[0]);
  Fmatrix_2d(&i_D_lga->sqsum_by_cat[0][0],&i_D_lga->sqsum_by_cat[0]);
  Fmatrix_3d(&i_D_lga->glm->beta_hat[0][0][0],&i_D_lga->glm->beta_hat[0][0],&i_D_lga->glm->beta_hat[0]);
  FREE(i_D_lga->glm->ssq);
  Fmatrix_3d(&i_D_lga->glm->suff[0][0][0],&i_D_lga->glm->suff[0][0],&i_D_lga->glm->suff[0]);

  Fmatrix_2d(&i_D_lga->Ages_fix[0][0],&i_D_lga->Ages_fix[0]);
  Fmatrix_2d(&i_D_lga->sum_by_cat_fix[0][0],&i_D_lga->sum_by_cat_fix[0]);
  Fmatrix_2d(&i_D_lga->sqsum_by_cat_fix[0][0],&i_D_lga->sqsum_by_cat_fix[0]);
  Fmatrix_3d(&i_D_lga->glm->beta_hat_fix[0][0][0],&i_D_lga->glm->beta_hat_fix[0][0],&i_D_lga->glm->beta_hat_fix[0]);
  FREE(i_D_lga->glm->ssq_fix);
  Fmatrix_3d(&i_D_lga->glm->suff_fix[0][0][0],&i_D_lga->glm->suff_fix[0][0],&i_D_lga->glm->suff_fix[0]);
  
  return(0);
}		/* end of re_makedata_lga_suff */



/*!
  \author Hanne Rognebakke
  \brief Calculates sufficient statistics for the lga model when coastal cod 

  Use amigo data and age-stratified by length data.
  Use certain observations for starting values even if classification error

  This routine is only used in initialization. 
  It calculates summary statistics based on aged fish.
  Updating sufficient statistics based on length-only data is performed
  by the ::make_suff_lga routine.

  \todo The routine copy summaries to the fix-versions. This should be made
  more clean later on.

  Memory allocated in this routine is reallocated in ::re_makedata_lga_suff.
*/
int makedata_lga_suff_CC(Data_lin *i_D_lga,Data_lin *i_D_lga_CC,Data_orig *i_D_orig,
			 Data_g_a *i_D_g_a,Data_g_a *i_D_g_a_CC,int i_class_error)
{
  int        ncat,ncat_CC,a,b,f,h,season,N,ind_f,ind_a,cum_fish,cCC,cS,err=0;
  double     beta0,beta1,ssq,length,r;


  ncat = i_D_g_a->ncat;
  ncat_CC = i_D_g_a_CC->ncat;

  /* Sufficient Statistics */
  i_D_lga->Ages = Mmatrix_2d(0,i_D_lga->glm->nHaul,0,ncat-1,sizeof(int),1); // Free ok
  i_D_lga->glm->haul_obs = CALLOC(i_D_lga->glm->nHaul,int);   // Free ok
  i_D_lga->glm->boat_obs = CALLOC(i_D_orig->nBoat,int);   // Free ok
  i_D_lga->sum_by_cat = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->sqsum_by_cat = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,0,ncat-1,sizeof(double),1);    // Free ok
  /* Need size of beta_hat equal to nxcov when finding b-vector in sample_gauss_eff */
  i_D_lga->glm->beta_hat = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,0,
				      0,i_D_lga->glm->nxcov-1,sizeof(double),1);    // Free ok
  i_D_lga->glm->ssq = CALLOC(i_D_lga->glm->nHaul,double);                 // Free ok
  i_D_lga->glm->suff = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,i_D_lga->glm->nxcov-1,0,
				  i_D_lga->glm->nxcov-1,sizeof(double),1); // Free ok
 
  /* Sufficient statistics based on observed age-length, constant during the simulations */
  i_D_lga->Ages_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul,0,ncat-1,sizeof(int),1); // Free ok
  i_D_lga->sum_by_cat_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->sqsum_by_cat_fix = Mmatrix_2d(0,i_D_lga->glm->nHaul-1,0,ncat-1,sizeof(double),1);    // Free ok
  i_D_lga->glm->beta_hat_fix = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,0,
					  0,i_D_lga->glm->nxcov-1,sizeof(double),1);    // Free ok
  i_D_lga->glm->ssq_fix = CALLOC(i_D_lga->glm->nHaul,double);                 // Free ok
  i_D_lga->glm->suff_fix = Mmatrix_3d(0,i_D_lga->glm->nHaul-1,0,i_D_lga->glm->nxcov-1,
				      0,i_D_lga->glm->nxcov-1,sizeof(double),1); // Free ok
 
  /* Coastal cod */
  i_D_lga_CC->Ages = Mmatrix_2d(0,i_D_lga_CC->glm->nHaul,0,ncat_CC-1,sizeof(int),1); // Free ok
  i_D_lga_CC->glm->haul_obs = CALLOC(i_D_lga->glm->nHaul,int);   // Free ok
  i_D_lga_CC->glm->boat_obs = CALLOC(i_D_orig->nBoat,int);   // Free ok
  i_D_lga_CC->sum_by_cat = Mmatrix_2d(0,i_D_lga_CC->glm->nHaul-1,0,ncat_CC-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->sqsum_by_cat =  Mmatrix_2d(0,i_D_lga_CC->glm->nHaul-1,0,ncat_CC-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->glm->beta_hat = Mmatrix_3d(0,i_D_lga_CC->glm->nHaul-1,0,0,
					 0,i_D_lga_CC->glm->nxcov-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->glm->ssq = CALLOC(i_D_lga_CC->glm->nHaul,double);             // Free ok
  i_D_lga_CC->glm->suff = Mmatrix_3d(0,i_D_lga_CC->glm->nHaul-1,0,i_D_lga_CC->glm->nxcov-1,
				     0,i_D_lga_CC->glm->nxcov-1,sizeof(double),1); // Free ok
  i_D_lga_CC->Ages_fix = Mmatrix_2d(0,i_D_lga_CC->glm->nHaul,0,ncat_CC-1,sizeof(int),1); // Free ok
  i_D_lga_CC->sum_by_cat_fix = Mmatrix_2d(0,i_D_lga_CC->glm->nHaul-1,0,ncat_CC-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->sqsum_by_cat_fix = Mmatrix_2d(0,i_D_lga_CC->glm->nHaul-1,0,ncat_CC-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->glm->beta_hat_fix = Mmatrix_3d(0,i_D_lga_CC->glm->nHaul-1,0,0,
					     0,i_D_lga_CC->glm->nxcov-1,sizeof(double),1);    // Free ok
  i_D_lga_CC->glm->ssq_fix = CALLOC(i_D_lga_CC->glm->nHaul,double);                 // Free ok
  i_D_lga_CC->glm->suff_fix = Mmatrix_3d(0,i_D_lga_CC->glm->nHaul-1,0,i_D_lga_CC->glm->nxcov-1,
					 0,i_D_lga_CC->glm->nxcov-1,sizeof(double),1); // Free ok

  N = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    N = max(N,i_D_orig->nFishBoat[h]);
  ind_f = 0;
  if(i_D_lga->glm->nHaul != i_D_lga_CC->glm->nHaul)
    write_warning("makedata_lga_suff_CC: something is wrong\n");

  for(b=0;b<i_D_orig->nBoat;b++)
    {
      i_D_lga->glm->boat_obs[b] = 0;
      i_D_lga_CC->glm->boat_obs[b] = 0;
    }

  cum_fish = 0;
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      i_D_lga->glm->haul_obs[h] = 0;
      i_D_lga_CC->glm->haul_obs[h] = 0;
      for(a=0;a<ncat;a++)
	{
	  i_D_lga->Ages[h][a] = 0;
	  i_D_lga->sum_by_cat[h][a] = G_ZERO;
	  i_D_lga->sqsum_by_cat[h][a] = G_ZERO;
	  i_D_lga->Ages_fix[h][a] = 0;
	  i_D_lga->sum_by_cat_fix[h][a] = G_ZERO;
	  i_D_lga->sqsum_by_cat_fix[h][a] = G_ZERO;
	}
      for(a=0;a<ncat_CC;a++)
	{
	  i_D_lga_CC->Ages[h][a] = 0;
	  i_D_lga_CC->sum_by_cat[h][a] = G_ZERO;
	  i_D_lga_CC->sqsum_by_cat[h][a] = G_ZERO;
	  i_D_lga_CC->Ages_fix[h][a] = 0;
	  i_D_lga_CC->sum_by_cat_fix[h][a] = G_ZERO;
	  i_D_lga_CC->sqsum_by_cat_fix[h][a] = G_ZERO;
	}
      i_D_lga->glm->suff[h][0][1] = G_ZERO;
      i_D_lga->glm->suff[h][1][1] = G_ZERO;
      i_D_lga_CC->glm->suff[h][0][1] = G_ZERO;
      i_D_lga_CC->glm->suff[h][1][1] = G_ZERO;
      for(f=0;f<i_D_orig->nFishBoat[h]-i_D_orig->num_noAge[h];f++) //both age and length observed
	{
	  ind_f = cum_fish+f;
	  a = i_D_orig->totage[ind_f];
	  season = i_D_orig->season[ind_f];
	  length = i_D_orig->totlength[ind_f];
	  r = (double) i_D_orig->replength[ind_f];
	  assert(a > -1000);
	  assert(length > -1000);
	  cCC = 0;
	  cS = 0;
	  /* Use certain observations for starting values even if classification error */
	  if(i_D_orig->tottype[ind_f]==1) 
	    cCC = 1;
	  else if((i_D_orig->tottype[ind_f]==2) && (i_class_error==0))
	    cCC = 1;
	  else if((i_D_orig->tottype[ind_f]==4) && (i_class_error==0))
	    cS = 1;
	  else if(i_D_orig->tottype[ind_f]==5) 
	    cS = 1;

	  if((i_D_orig->tottype[ind_f]==1) | (i_D_orig->tottype[ind_f]==2))
	    i_D_lga_CC->glm->haul_obs[h] = 1;
	  if((i_D_orig->tottype[ind_f]==5) | (i_D_orig->tottype[ind_f]==4) )
	    i_D_lga->glm->haul_obs[h] = 1;

	  if(cCC) //certain coastal cod
	    {
	      ind_a = i_D_g_a_CC->a2Age_vec[a-(int) i_D_g_a_CC->a_vec[0]]+(season-1);
	      i_D_lga_CC->Ages[h][ind_a] += i_D_orig->replength[ind_f];
	      i_D_lga_CC->sum_by_cat[h][ind_a] += r*length;
	      i_D_lga_CC->sqsum_by_cat[h][ind_a] += r*length*length;
	      i_D_lga_CC->Ages_fix[h][ind_a] += i_D_orig->replength[ind_f];
	      i_D_lga_CC->sum_by_cat_fix[h][ind_a] += r*length;
	      i_D_lga_CC->sqsum_by_cat_fix[h][ind_a] += r*length*length;
	      i_D_lga_CC->glm->suff[h][0][0] += r;
	      i_D_lga_CC->glm->suff[h][0][1] += r*i_D_g_a_CC->g_a[ind_a];
	      i_D_lga_CC->glm->suff[h][1][1] += r*i_D_g_a_CC->g_a[ind_a]*i_D_g_a_CC->g_a[ind_a];
	    }
	  else if(cS) //certain skrei
	    {
	      ind_a = i_D_g_a->a2Age_vec[a-(int) i_D_g_a->a_vec[0]]+(season-1);
	      i_D_lga->Ages[h][ind_a] += i_D_orig->replength[ind_f];
	      i_D_lga->sum_by_cat[h][ind_a] += r*length;
	      i_D_lga->sqsum_by_cat[h][ind_a] += r*length*length;
	      i_D_lga->Ages_fix[h][ind_a] += i_D_orig->replength[ind_f];
	      i_D_lga->sum_by_cat_fix[h][ind_a] += r*length;
	      i_D_lga->sqsum_by_cat_fix[h][ind_a] += r*length*length;
	      i_D_lga->glm->suff[h][0][0] += r;
	      i_D_lga->glm->suff[h][0][1] += r*i_D_g_a->g_a[ind_a];
	      i_D_lga->glm->suff[h][1][1] += r*i_D_g_a->g_a[ind_a]*i_D_g_a->g_a[ind_a];
	    }
	}
      cum_fish += i_D_orig->nFishBoat[h];
      if(i_D_lga->glm->haul_obs[h]==1)
        {
          b = i_D_orig->boat[h]-1;
          i_D_lga->glm->boat_obs[b] = 1;
        }
      if(i_D_lga_CC->glm->haul_obs[h]==1)
        {
          b = i_D_orig->boat[h]-1;
          i_D_lga_CC->glm->boat_obs[b] = 1;
        }

      i_D_lga->glm->suff[h][1][0] = i_D_lga->glm->suff[h][0][1];
      i_D_lga_CC->glm->suff[h][1][0] = i_D_lga_CC->glm->suff[h][0][1];
      
      if(i_D_lga->glm->suff[h][0][0]>G_ZERO)
	{
	  err = lm_fit_suff(ncat,i_D_g_a->g_a,i_D_lga->sum_by_cat[h],i_D_lga->sqsum_by_cat[h],
			    i_D_lga->glm->suff[h],&beta0,&beta1,&ssq);
	  if(err)
	    {
	      write_warning("makedata_lga_suff_CC:Error calling lm_fit_suff\n");
	      return(err);
	    }
	  i_D_lga->glm->beta_hat[h][0][0] = beta0;
	  i_D_lga->glm->beta_hat[h][0][1] = beta1;
	  i_D_lga->glm->ssq[h] = ssq;
	}
      //if(h<10)
      //fprintf(stderr,"h=%d:ncat=%d,sum_by_cat[0]=%f,sum_by_cat[1]=%f,sum_by_cat[2]=%f,sqsum_by_cat[0]=%f,sqsum_by_cat[1]=%f,sqsum_by_cat[2]=%f,suff[0][0]=%f,suff[0][1]=%f,suff[1][1]=%f,beta0=%f,beta1=%f\n",
      //	h,ncat,i_D_lga->sum_by_cat[h][0],i_D_lga->sum_by_cat[h][1],i_D_lga->sum_by_cat[h][2],
      //	i_D_lga->sqsum_by_cat[h][0],i_D_lga->sqsum_by_cat[h][1],i_D_lga->sqsum_by_cat[h][2],
      //	i_D_lga->glm->suff[h][0][0],i_D_lga->glm->suff[h][0][1],i_D_lga->glm->suff[h][1][1],
      //	beta0,beta1);
      if(i_D_lga_CC->glm->suff[h][0][0]>G_ZERO)
	{
	  err = lm_fit_suff(ncat_CC,i_D_g_a_CC->g_a,i_D_lga_CC->sum_by_cat[h],
			    i_D_lga_CC->sqsum_by_cat[h],
			    i_D_lga_CC->glm->suff[h],&beta0,&beta1,&ssq);
	  if(err)
	    {
	      write_warning("makedata_lga_suff_CC:Error calling lm_fit_suff\n");
	      return(err);
	    }
	  i_D_lga_CC->glm->beta_hat[h][0][0] = beta0;
	  i_D_lga_CC->glm->beta_hat[h][0][1] = beta1;
	  i_D_lga_CC->glm->ssq[h] = ssq;
	}


      // Copying summaries to the fix-versions. 
      // This should probably be changed to put all directly into the fix-versions
      i_D_lga->glm->suff_fix[h][0][0] = i_D_lga->glm->suff[h][0][0];
      i_D_lga->glm->suff_fix[h][0][1] = i_D_lga->glm->suff[h][0][1];
      i_D_lga->glm->suff_fix[h][1][0] = i_D_lga->glm->suff[h][1][0];
      i_D_lga->glm->suff_fix[h][1][1] = i_D_lga->glm->suff[h][1][1];
      if(i_D_lga->glm->suff[h][0][0]>G_ZERO)
	{
	  i_D_lga->glm->beta_hat_fix[h][0][0] = i_D_lga->glm->beta_hat[h][0][0];
	  i_D_lga->glm->beta_hat_fix[h][0][1] = i_D_lga->glm->beta_hat[h][0][1];
	  i_D_lga->glm->ssq_fix[h] = i_D_lga->glm->ssq[h];
	}
      else
	{
	  i_D_lga->glm->beta_hat_fix[h][0][0] = G_ZERO;
	  i_D_lga->glm->beta_hat_fix[h][0][1] = G_ZERO;
	  i_D_lga->glm->ssq_fix[h] = G_ZERO;
	}
      i_D_lga_CC->glm->suff_fix[h][0][0] = i_D_lga_CC->glm->suff[h][0][0];
      i_D_lga_CC->glm->suff_fix[h][0][1] = i_D_lga_CC->glm->suff[h][0][1];
      i_D_lga_CC->glm->suff_fix[h][1][0] = i_D_lga_CC->glm->suff[h][1][0];
      i_D_lga_CC->glm->suff_fix[h][1][1] = i_D_lga_CC->glm->suff[h][1][1];
      if(i_D_lga_CC->glm->suff[h][0][0]>G_ZERO)
	{
	  i_D_lga_CC->glm->beta_hat_fix[h][0][0] = i_D_lga_CC->glm->beta_hat[h][0][0];
	  i_D_lga_CC->glm->beta_hat_fix[h][0][1] = i_D_lga_CC->glm->beta_hat[h][0][1];
	  i_D_lga_CC->glm->ssq_fix[h] = i_D_lga_CC->glm->ssq[h];
	}
      else
	{
	  i_D_lga_CC->glm->beta_hat_fix[h][0][0] = G_ZERO;
	  i_D_lga_CC->glm->beta_hat_fix[h][0][1] = G_ZERO;
	  i_D_lga_CC->glm->ssq_fix[h] = G_ZERO;
	}
    }

  return(0);
}		/* end of makedata_lga_suff_CC */



/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in ::makedata_lga_suff_CC
*/
int re_makedata_lga_suff_CC(Data_lin *i_D_lga,Data_lin *i_D_lga_CC)
{
  Fmatrix_2d(&i_D_lga->Ages[0][0],&i_D_lga->Ages[0]);
  FREE(i_D_lga->glm->haul_obs);
  FREE(i_D_lga->glm->boat_obs);
  Fmatrix_2d(&i_D_lga->sum_by_cat[0][0],&i_D_lga->sum_by_cat[0]);
  Fmatrix_2d(&i_D_lga->sqsum_by_cat[0][0],&i_D_lga->sqsum_by_cat[0]);
  Fmatrix_3d(&i_D_lga->glm->beta_hat[0][0][0],&i_D_lga->glm->beta_hat[0][0],&i_D_lga->glm->beta_hat[0]);
  FREE(i_D_lga->glm->ssq);
  Fmatrix_3d(&i_D_lga->glm->suff[0][0][0],&i_D_lga->glm->suff[0][0],&i_D_lga->glm->suff[0]);
  Fmatrix_2d(&i_D_lga->Ages_fix[0][0],&i_D_lga->Ages_fix[0]);
  Fmatrix_2d(&i_D_lga->sum_by_cat_fix[0][0],&i_D_lga->sum_by_cat_fix[0]);
  Fmatrix_2d(&i_D_lga->sqsum_by_cat_fix[0][0],&i_D_lga->sqsum_by_cat_fix[0]);
  Fmatrix_3d(&i_D_lga->glm->beta_hat_fix[0][0][0],&i_D_lga->glm->beta_hat_fix[0][0],&i_D_lga->glm->beta_hat_fix[0]);
  FREE(i_D_lga->glm->ssq_fix);
  Fmatrix_3d(&i_D_lga->glm->suff_fix[0][0][0],&i_D_lga->glm->suff_fix[0][0],&i_D_lga->glm->suff_fix[0]);
  
  Fmatrix_2d(&i_D_lga_CC->Ages[0][0],&i_D_lga_CC->Ages[0]);
  FREE(i_D_lga_CC->glm->haul_obs);
  FREE(i_D_lga_CC->glm->boat_obs);
   Fmatrix_2d(&i_D_lga_CC->sum_by_cat[0][0],&i_D_lga_CC->sum_by_cat[0]);
  Fmatrix_2d(&i_D_lga_CC->sqsum_by_cat[0][0],&i_D_lga_CC->sqsum_by_cat[0]);
  Fmatrix_3d(&i_D_lga_CC->glm->beta_hat[0][0][0],&i_D_lga_CC->glm->beta_hat[0][0],&i_D_lga_CC->glm->beta_hat[0]);
  FREE(i_D_lga_CC->glm->ssq);
  Fmatrix_3d(&i_D_lga_CC->glm->suff[0][0][0],&i_D_lga_CC->glm->suff[0][0],&i_D_lga_CC->glm->suff[0]);
  Fmatrix_2d(&i_D_lga_CC->Ages_fix[0][0],&i_D_lga_CC->Ages_fix[0]);
  Fmatrix_2d(&i_D_lga_CC->sum_by_cat_fix[0][0],&i_D_lga_CC->sum_by_cat_fix[0]);
  Fmatrix_2d(&i_D_lga_CC->sqsum_by_cat_fix[0][0],&i_D_lga_CC->sqsum_by_cat_fix[0]);
  Fmatrix_3d(&i_D_lga_CC->glm->beta_hat_fix[0][0][0],&i_D_lga_CC->glm->beta_hat_fix[0][0],&i_D_lga_CC->glm->beta_hat_fix[0]);
  FREE(i_D_lga_CC->glm->ssq_fix);
  Fmatrix_3d(&i_D_lga_CC->glm->suff_fix[0][0][0],&i_D_lga_CC->glm->suff_fix[0][0],&i_D_lga_CC->glm->suff_fix[0]);
  

  return(0);
}		/* end of re_makedata_lga_suff_CC */


/*!
  \author Geir Storvik
  \brief Calculates sufficient statistics for the wgl model 

  Memory allocated in this routine is reallocated in ::re_makedata_wgl_suff.
*/
int makedata_wgl_suff(Data_lin *i_D_lin,Data_orig *i_D_orig)
{
  int        f,h,n,N,ind,err=0;
  double    *x,*y,beta0,beta1,ssq,l,w,r;


  /* Sufficient Statistics */
  i_D_lin->glm->beta_hat = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,0,
				      0,i_D_lin->glm->nxcov-1,sizeof(double),1);  // Free ok
  i_D_lin->glm->ssq = CALLOC(i_D_lin->glm->nHaul,double);      // Free ok
  i_D_lin->glm->suff = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,i_D_lin->glm->nxcov-1,
				  0,i_D_lin->glm->nxcov-1,sizeof(double),1); // Free ok

  N = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    N = max(N,i_D_orig->nFishBoat[h]);

  x = CALLOC(N,double);               // Free ok
  y = CALLOC(N,double);               // Free ok

  ind = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    {
      i_D_lin->glm->suff[h][0][0] = G_ZERO;
      i_D_lin->glm->suff[h][0][1] = G_ZERO;
      i_D_lin->glm->suff[h][1][1] = G_ZERO;
      n = 0;
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
	  l = i_D_orig->totlength[ind];
	  w = i_D_orig->totweight[ind];
	  if(l > -9000 && w > -9000)
	    {
	      r = (double) i_D_orig->replength[ind];
	      x[n] = l;
	      y[n] = w;
	      i_D_lin->glm->suff[h][0][0] += r;
	      i_D_lin->glm->suff[h][0][1] += r*x[n];
	      i_D_lin->glm->suff[h][1][1] += r*x[n]*x[n];
	      n += 1;
	    }
	  ind++;
	}
      if(n > 0)
	{
	  err = lm_fit(n,x,y,&beta0,&beta1,&ssq); 
	  if(err)
	    {
	      write_warning("makedata_wgl_suff:Error calling lm_fit_suff\n");
	      return(err);
	    }
          #ifdef DEBUG_PROG
	  if(h<5)
	    printf("makedata_wgl_suff: h=%d,suff00=%f,suff01=%f,suff11=%f,beta0=%f,beta1=%f\n",
		   h,i_D_lin->glm->suff[h][0][0],
	  	   i_D_lin->glm->suff[h][0][1],i_D_lin->glm->suff[h][1][1],beta0,beta1);
	  #endif
          i_D_lin->glm->beta_hat[h][0][0] = beta0;
          i_D_lin->glm->beta_hat[h][0][1] = beta1;
          i_D_lin->glm->ssq[h] = ssq;
	}
      else
	{
	  fprintf(stderr,"Warning: makedata_wgl_suff: h=%d, n=%d, l=%f, w=%f\n",h,n,l,w);
	}
      i_D_lin->glm->suff[h][1][0] = i_D_lin->glm->suff[h][0][1];
    }
  FREE(x);
  FREE(y);

  return(0);
}		/* end of makedata_wgl_suff */


/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in ::makedata_wgl_suff
*/
int re_makedata_wgl_suff(Data_lin *i_D_lin)
{
  Fmatrix_3d(&i_D_lin->glm->beta_hat[0][0][0],&i_D_lin->glm->beta_hat[0][0],&i_D_lin->glm->beta_hat[0]);
  FREE(i_D_lin->glm->ssq);
  Fmatrix_3d(&i_D_lin->glm->suff[0][0][0],&i_D_lin->glm->suff[0][0],&i_D_lin->glm->suff[0]);

  return(0);
}		/* end of re_makedata_wgl_suff */


/*!
  \author Hanne Rognebakke
  \brief Calculates sufficient statistics for the wgl model when coastal cod

  Memory allocated in this routine is reallocated in ::re_makedata_wgl_suff_CC.
*/
int makedata_wgl_suff_CC(Data_lin *i_D_lin,Data_lin *i_D_lin_CC,Data_orig *i_D_orig)
{
  int        f,h,N,ind,n_CC,n,err=0;
  double     beta0,beta1,ssq,l,w,r;
  double    *x_CC,*y_CC,*x,*y;



  /* Sufficient Statistics */
  i_D_lin->glm->beta_hat = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,0,
				      0,i_D_lin->glm->nxcov-1,sizeof(double),1);  // Free ok
  i_D_lin->glm->ssq = CALLOC(i_D_lin->glm->nHaul,double);      // Free ok
  i_D_lin->glm->suff = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,i_D_lin->glm->nxcov-1,
				  0,i_D_lin->glm->nxcov-1,sizeof(double),1); // Free ok

  i_D_lin_CC->glm->beta_hat = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,0,
				      0,i_D_lin->glm->nxcov-1,sizeof(double),1);  // Free ok
  i_D_lin_CC->glm->ssq = CALLOC(i_D_lin->glm->nHaul,double);      // Free ok
  i_D_lin_CC->glm->suff = Mmatrix_3d(0,i_D_lin->glm->nHaul-1,0,i_D_lin->glm->nxcov-1,
				  0,i_D_lin->glm->nxcov-1,sizeof(double),1); // Free ok

  N = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    N = max(N,i_D_orig->nFishBoat[h]);

  x_CC = CALLOC(N,double);               // Free ok
  y_CC = CALLOC(N,double);               // Free ok
  x = CALLOC(N,double);               // Free ok
  y = CALLOC(N,double);               // Free ok

  ind = 0;
  for(h=0;h<i_D_lin->glm->nHaul;h++)
    {
      i_D_lin->glm->suff[h][0][0] = G_ZERO;
      i_D_lin->glm->suff[h][0][1] = G_ZERO;
      i_D_lin->glm->suff[h][1][1] = G_ZERO;
      n_CC = 0;
      n = 0;
      for(f=0;f<i_D_orig->nFishBoat[h];f++)
	{
	  l = i_D_orig->totlength[ind];
	  w = i_D_orig->totweight[ind];
	  if(l > -9000 && w > -9000)
	    {
	      if(i_D_orig->tottype[ind]==1) //certain coastal cod
		{
		  r = (double) i_D_orig->replength[ind];
		  x_CC[n_CC] = l;
		  y_CC[n_CC] = w;
		  i_D_lin_CC->glm->suff[h][0][0] += r;
		  i_D_lin_CC->glm->suff[h][0][1] += r*l;
		  i_D_lin_CC->glm->suff[h][1][1] += r*l*l;
		  n_CC += 1;
		}
	      else if(i_D_orig->tottype[ind]==5) //certain skrei
		{
		  r = (double) i_D_orig->replength[ind];
		  x[n] = l;
		  y[n] = w;
		  i_D_lin->glm->suff[h][0][0] += r;
		  i_D_lin->glm->suff[h][0][1] += r*l;
		  i_D_lin->glm->suff[h][1][1] += r*l*l;
		  n += 1;
		}
	    }
	  ind++;
	}

      if(n_CC > 0)
	{
	  err = lm_fit(n_CC,x_CC,y_CC,&beta0,&beta1,&ssq); 
	  if(err)
	    {
	      write_warning("makedata_wgl_suff_CC:Error calling lm_fit_suff\n");
	      return(err);
	    }
          i_D_lin_CC->glm->beta_hat[h][0][0] = beta0;
          i_D_lin_CC->glm->beta_hat[h][0][1] = beta1;
          i_D_lin_CC->glm->ssq[h] = ssq;
	}
      i_D_lin_CC->glm->suff[h][1][0] = i_D_lin_CC->glm->suff[h][0][1];

      if(n > 0)
	{
	  err = lm_fit(n,x,y,&beta0,&beta1,&ssq); 
	  if(err)
	    {
	      write_warning("makedata_wgl_suff_CC:Error calling lm_fit_suff\n");
	      return(err);
	    }
          i_D_lin->glm->beta_hat[h][0][0] = beta0;
          i_D_lin->glm->beta_hat[h][0][1] = beta1;
          i_D_lin->glm->ssq[h] = ssq;
	}
      i_D_lin->glm->suff[h][1][0] = i_D_lin->glm->suff[h][0][1];


      if(i_D_lin_CC->glm->nxcov==3)
	{
	  i_D_lin_CC->glm->suff[h][0][2] = i_D_lin_CC->glm->suff[h][0][0] * i_D_orig->haulweight[h];
	  i_D_lin_CC->glm->suff[h][1][2] = i_D_lin_CC->glm->suff[h][0][1]*i_D_orig->haulweight[h];
	  i_D_lin_CC->glm->suff[h][2][2] = i_D_lin_CC->glm->suff[h][0][0] * i_D_orig->haulweight[h]*i_D_orig->haulweight[h];
	  i_D_lin_CC->glm->suff[h][2][0] = i_D_lin_CC->glm->suff[h][0][2];
	  i_D_lin_CC->glm->suff[h][2][1] = i_D_lin_CC->glm->suff[h][1][2];
          i_D_lin_CC->glm->beta_hat[h][0][2] = G_ZERO;    /* needed in find b-vector in sample_gauss_eff() */
	}
      if(i_D_lin->glm->nxcov==3)
	{
	  i_D_lin->glm->suff[h][0][2] = i_D_lin->glm->suff[h][0][0] * i_D_orig->haulweight[h];
	  i_D_lin->glm->suff[h][1][2] = i_D_lin->glm->suff[h][0][1]*i_D_orig->haulweight[h];
	  i_D_lin->glm->suff[h][2][2] = i_D_lin->glm->suff[h][0][0] * i_D_orig->haulweight[h]*i_D_orig->haulweight[h];
	  i_D_lin->glm->suff[h][2][0] = i_D_lin->glm->suff[h][0][2];
	  i_D_lin->glm->suff[h][2][1] = i_D_lin->glm->suff[h][1][2];
          i_D_lin->glm->beta_hat[h][0][2] = G_ZERO;    /* needed in find b-vector in sample_gauss_eff() */
	}
    }
  FREE(x_CC);
  FREE(y_CC);
  FREE(x);
  FREE(y);

  return(0);
}		/* end of makedata_wgl_suff_CC */


/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in ::makedata_wgl_suff_CC
*/
int re_makedata_wgl_suff_CC(Data_lin *i_D_lin,Data_lin *i_D_lin_CC)
{
  Fmatrix_3d(&i_D_lin->glm->beta_hat[0][0][0],&i_D_lin->glm->beta_hat[0][0],&i_D_lin->glm->beta_hat[0]);
  FREE(i_D_lin->glm->ssq);
  Fmatrix_3d(&i_D_lin->glm->suff[0][0][0],&i_D_lin->glm->suff[0][0],&i_D_lin->glm->suff[0]);

  Fmatrix_3d(&i_D_lin_CC->glm->beta_hat[0][0][0],&i_D_lin_CC->glm->beta_hat[0][0],&i_D_lin_CC->glm->beta_hat[0]);
  FREE(i_D_lin_CC->glm->ssq);
  Fmatrix_3d(&i_D_lin_CC->glm->suff[0][0][0],&i_D_lin_CC->glm->suff[0][0],&i_D_lin_CC->glm->suff[0]);

  return(0);
}		/* end of re_makedata_wgl_suff_CC */


/*!
  \author Geir Storvik
  \brief Make a struct of type Data_cov describing the covariate structure. 

  This structure will typically be a part of a Data_glm struct.

  Memory allocated by this routine is reallocated by ::re_make_c_cov
*/
static int make_c_cov(int i_n_cov,int *i_n_lev,int i_ispat,int i_iboat,int i_ihaulsize,int *i_fix,int *i_interaction,
		      int *i_in_landings,int **i_x_cov,int i_nHaul,Data_cov **o_xcov,int i_fit,int i_incl_haul)
{
  int       f,h,ncol,ind_col;
  Data_cov *xcov;

  xcov = CALLOC(1,Data_cov);    // Free ok
  /*number of factors */
  xcov->n_cov = i_n_cov;
  xcov->n_fac = CALLOC(xcov->n_cov,int);  // Free ok
  xcov->fix = CALLOC(xcov->n_cov,int);    // Free ok
  xcov->interaction = CALLOC(xcov->n_cov,int);    // Free ok
  xcov->in_landings = CALLOC(xcov->n_cov,int);    // Free ok
  /* index for spatial term */
  xcov->ispat = i_ispat;   
  /* index for boat term */
  xcov->iboat = i_iboat;   
  /* fixed effect */
  for(f=0;f<i_n_cov;f++)
    {
      xcov->n_fac[f] = i_n_lev[f];
      xcov->fix[f] = i_fix[f];
      xcov->interaction[f] = i_interaction[f];
      xcov->in_landings[f] = i_in_landings[f];
      //fprintf(stderr,"n_fac[%d]=%d, fix[%d]=%d, interaction[%d]=%d\n",
      //	      f,xcov->n_fac[f],f,xcov->fix[f],f,xcov->interaction[f]);
    }
  xcov->c_cov = CALLOC(i_nHaul,int *);            // Free ok
  for(h=0;h<i_nHaul;h++)
    xcov->c_cov[h] = CALLOC(xcov->n_cov,int);   // Free ok
  if(i_incl_haul)
    ncol = xcov->n_cov-1;
  else
    ncol = xcov->n_cov;
  if(i_fit)
    {
      ind_col = 0;
      for(f=0;f<ncol;f++)
	{
	  if(f!=i_ihaulsize)
	    {
	      for(h=0;h<i_nHaul;h++)
		{
		  xcov->c_cov[h][f] = i_x_cov[ind_col][h];
		}
	      ind_col++;
	    }
	}
      if(i_incl_haul)  // insert haul covariate
	{
	  for(h=0;h<i_nHaul;h++)
	    xcov->c_cov[h][ncol] = (h+1);
	}
      if(i_ihaulsize>0) // insert haulsize covariate
	{
	  for(h=0;h<i_nHaul;h++)
	    xcov->c_cov[h][i_ihaulsize] = 1;
	}
    }

  #ifdef DEBUG_PROG
  fprintf(stderr,"make_c_cov:");
  fprintf(stderr,"nHaul=%d,ncov=%d\n",i_nHaul,i_n_cov);
  fprintf(stderr,"i_incl_haul=%d\n",i_incl_haul);
  for(h=0;h<MIN(10,i_nHaul);h++)
    {
      for(f=0;f<i_n_cov;f++)
	fprintf(stderr,"%d ",xcov->c_cov[h][f]);
      fprintf(stderr,"\n");
    }
  fprintf(stderr,"\n\n"); 
  #endif
  
  *o_xcov = xcov;
  return(0);
}		/* end of make_c_cov */


 
/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in ::make_c_cov
*/
static int re_make_c_cov(int i_nHaul,Data_cov **o_xcov)
{
  int       h;
  Data_cov *xcov;

  xcov = *o_xcov;

  FREE(xcov->n_fac);
  FREE(xcov->fix);
  FREE(xcov->interaction);
  FREE(xcov->in_landings);

  for(h=0;h<i_nHaul;h++)
    FREE(xcov->c_cov[h]);
  FREE(xcov->c_cov);
      
  FREE(xcov);

  return(0);
}		/* end of re_make_c_cov */

    

/*!
  \author Geir Storvik
  \brief Make spatial structure 

  Memory allocated in this routine is reallocated by ::re_make_spat_struct
*/
static int make_spat_struct(int *i_num,int *i_adj_area,Data_cov *x_xcov)
{
  char  string[150];
  int   ind,k,r,rr,nArea,isp;

  isp = x_xcov->ispat;
  nArea = x_xcov->n_fac[isp];
  x_xcov->num_adj_area = CALLOC(nArea,int);          // Free ok
  x_xcov->adj_area = CALLOC(nArea,int *);            // Free ok

  for(r=0;r<nArea;r++)
    {
      x_xcov->num_adj_area[r] = i_num[r];
      x_xcov->adj_area[r] = CALLOC(x_xcov->num_adj_area[r],int);  // Free ok
    }

  ind = 0;
  for(r=0;r<nArea;r++)
    {
      for(rr=0;rr<x_xcov->num_adj_area[r];rr++)
	{
       
	  x_xcov->adj_area[r][rr] = i_adj_area[ind]-1;
          ind++;
	  /* convert neighbor areas  */
          k = 0;
          while(k < (nArea-1)  && x_xcov->adj_area[r][rr]!=x_xcov->conv[isp][k])
	    k++;
          if(x_xcov->adj_area[r][rr]==x_xcov->conv[isp][k])
	    x_xcov->adj_area[r][rr] = k;
          else
	    {
              printf("make_spat_struct:Missing neighbor area %d %d %d\n",
		     r,rr,(int) i_adj_area[ind]);
              sprintf(string,"make_spat_struct:Missing neighbor area %d %d %d\n",
		      r,rr,(int) i_adj_area[ind]);
              write_warning(string);
	      return(1);
	    }
	}
    }

  return(0);
}		/* end of make_spat_struct */

    
/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in ::make_spat_struct
*/
static int re_make_spat_struct(int *i_num,int *i_adj_area,Data_cov *x_xcov)
{
  int   r,nArea,isp;

  isp = x_xcov->ispat;
  nArea = x_xcov->n_fac[isp];

  FREE(x_xcov->num_adj_area);

  for(r=0;r<nArea;r++)
    FREE(x_xcov->adj_area[r]);
  FREE(x_xcov->adj_area);


  return(0);
}		/* end of re_make_spat_struct */



/*!
  \author Geir Storvik
  \brief make a struct of type Data_l containing length only and age-stratified-by-length data.

  Memory allocated by this routine is reallocated by ::re_makedata_only_length
*/
int makedata_only_length(int i_nLengths,int *i_lengthCount,double *i_length,int *i_journey,
                         int i_lga_nAgeLengths,double *i_lga_ageLength,
                         int *i_lga_ageLengthCount,int *i_lga_ageJourney,int i_nAges,
                         Data_lin *i_D_lga,Data_l **o_D_l)
{
  int         a,i,ind,l;
  Data_l     *D_l;

  D_l = CALLOC(1,Data_l);     // Free ok

  D_l->nLengths = (int) i_nLengths;

  if(D_l->nLengths>0)
    {
      D_l->count = CALLOC(D_l->nLengths,int);   // Free ok
      for(l=0;l<D_l->nLengths;l++)
	D_l->count[l] = (int) i_lengthCount[l];

      D_l->length = CALLOC(D_l->nLengths,double); // Free ok
      for(l=0;l<D_l->nLengths;l++)
	D_l->length[l] = i_length[l];

      D_l->journey = CALLOC(D_l->nLengths,int);   // Free ok
      for(l=0;l<D_l->nLengths;l++)
	{
	  D_l->journey[l] = (int) i_journey[l];
	  D_l->journey[l]--;    /* Change for hauls starting at zero */
	}
    }
  /* Aged fish */
  D_l->nAgeLengths = i_lga_nAgeLengths;
  if(D_l->nAgeLengths>0)
    {
      D_l->ageLength = i_lga_ageLength;
      D_l->ageLengthCount = 
	Mmatrix_2d(0,i_lga_nAgeLengths-1,0,i_nAges-1,sizeof(int),1); // Free ok
      ind = 0;
      for(i=0;i<i_lga_nAgeLengths;i++)
	for(a=0;a<i_nAges;a++)
	  {
	    D_l->ageLengthCount[i][a] = i_lga_ageLengthCount[ind];
	    ind++;
	  }
      D_l->ageJourney = CALLOC(i_lga_nAgeLengths,int);   // Free ok
      for(i=0;i<i_lga_nAgeLengths;i++)
	{
	  D_l->ageJourney[i] = i_lga_ageJourney[i];
          D_l->ageJourney[i]--;
	}
    }
  *o_D_l = D_l;

  return(0);
}		/* end of makedata_only_length */


/*!
  \author Geir Storvik
  \brief Reallocate memory allocated in ::makedata_only_length
*/
int re_makedata_only_length(int i_nLengths,int *i_lengthCount,double *i_length,
			    int *i_journey,
			    int i_lga_nAgeLengths,double *i_lga_ageLength,
			    int *i_lga_ageLengthCount,int *i_lga_ageJourney,int i_nAges,
			    Data_lin *i_D_lga,Data_l **o_D_l)
{
  Data_l     *D_l;

  D_l = *o_D_l;

  if(D_l->nLengths>0)
    {
      FREE(D_l->count);
      FREE(D_l->length);
      FREE(D_l->journey);
    }
  /* Aged fish */
  if(D_l->nAgeLengths>0)
    {
      Fmatrix_2d(&D_l->ageLengthCount[0][0],&D_l->ageLengthCount[0]);
      FREE(D_l->ageJourney);
    }

  FREE(D_l);

  return(0);
}		/* end of re_makedata_only_length */


/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in ::makedata_only_length
*/
int makedata_mcmc(Input_predict *i_inPredict, int i_nAges, int i_nlint, TC_struct **o_totcatch)
{
  int n;
  TC_struct *totcatch;

  totcatch = CALLOC(1,TC_struct);     

  totcatch->catch_at_age = Mmatrix_2d(0,i_nAges-1,0,i_nlint-1,sizeof(double),1); 

  n = i_nAges*i_inPredict->N_l_int*(i_inPredict->nMCMC - i_inPredict->burnin);
  totcatch->mcmc = CALLOC(n,double);

  n = i_nAges*(i_inPredict->nMCMC - i_inPredict->burnin);
  totcatch->mean_l = CALLOC(n,double);
  totcatch->mean_w = CALLOC(n,double);


  *o_totcatch = totcatch;
 
  return(0);
}		/* end of makedata_mcmc */


/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in ::makedata_mcmc
*/
int re_makedata_mcmc(TC_struct **o_totcatch)
{
  TC_struct *totcatch;

  totcatch = *o_totcatch;
  Fmatrix_2d(&totcatch->catch_at_age[0][0],&totcatch->catch_at_age[0]);

  FREE(totcatch->mcmc);
  FREE(totcatch->mean_l);
  FREE(totcatch->mean_w);

  FREE(totcatch);
 
  return(0);
}		/* end of re_makedata_mcmc */


/*!
  \author Geir Storvik
  \brief Make a struct of type Data_totcatch containing  total catch (to be used for prediction).

  Memory allocated in this routine is reallocated by ::re_makedata_totcatch
*/
int makedata_totcatch(Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_wgl,Data_lin *i_D_hsz,
		      Input_totcatch *i_inCatch,int i_inc_hsz,Data_totcatch **o_D_totcatch)
{
  int        err,c,i,ii,ind,ncov;
  Data_totcatch  *D_totcatch;
 
  D_totcatch = CALLOC(1,Data_totcatch);        // Free ok
  
  /* Number of cells and number of factors */
  D_totcatch->nCell = i_inCatch->nCell;
  D_totcatch->nFactors = i_inCatch->nFactors;

  ncov = 1;
  // Factors for hauleffect is not part of fac_age_int anymore
  D_totcatch->fac_age = Mmatrix_2d(0,ncov-1,0,i_D_age->glm->xcov[0]->n_cov-1,sizeof(int),1);    // Free ok
  #ifdef DEBUG_PREDICT
  printf("Factors corresponding to age int model\n");
  #endif
  /* Factors corresponding to age int model */
  ind = 0;
  for(i=0;i<(i_D_age->glm->xcov[0]->n_cov-1);i++)
    {
      if(i_D_age->glm->xcov[0]->in_landings[i]==1)
	{
	  D_totcatch->fac_age[0][i] = i_inCatch->fac_age_int[ind]-1; /* Assume start on zero */
	  ind++;
	}
      else
	D_totcatch->fac_age[0][i] = -1;
    }
  if(i_D_age->glm->xcov[0]->icell>0) /* The factor for cell effects is added as the last column */
    {
      D_totcatch->fac_age[0][i_D_age->glm->xcov[0]->icell] = i_inCatch->fac_age_int[ind]-1;
    }
  #ifdef DEBUG_PREDICT
  for(i=0;i<(i_D_age->glm->xcov[0]->n_cov-1);i++)
    fprintf(stderr,"fac_age[0][%d]=%d\n",i,D_totcatch->fac_age[0][i]);
  #endif
    
  if(i_inc_hsz)
    {
      D_totcatch->fac_hsz = Mmatrix_2d(0,ncov-1,0,i_D_hsz->glm->xcov[0]->n_cov-1,sizeof(int),1);    // Free ok
      /* Factors corresponding to hsz int model */
      ind = 0;
      for(i=0;i<(i_D_hsz->glm->xcov[0]->n_cov-1);i++)
	{
	  if(i_D_age->glm->xcov[0]->in_landings[i]==1)
	    {
	      D_totcatch->fac_hsz[0][i] = i_inCatch->fac_hsz_int[ind]-1; /* Assume start on zero */
	      ind++;
	    }
	  else
	    D_totcatch->fac_hsz[0][i] = -1;
	}
      if(i_D_hsz->glm->xcov[0]->icell>0)
	{
	  D_totcatch->fac_hsz[0][i_D_hsz->glm->xcov[0]->icell] = i_inCatch->fac_hsz_int[ind]-1;
    	}
      #ifdef DEBUG_PREDICT
      for(i=0;i<(i_D_hsz->glm->xcov[0]->n_cov-1);i++)
	fprintf(stderr,"fac_hsz[0][%d]=%d\n",i,D_totcatch->fac_hsz[0][i]);
      #endif
    }

  ncov = 2;
  D_totcatch->fac_lga = Mmatrix_2d(0,ncov-1,0,i_D_lga->glm->xcov[0]->n_cov-1,sizeof(int),1);  // Free ok

  #ifdef DEBUG_PREDICT
  printf("Factors corresponding to lga intercept model\n");
  #endif
  /* Factors corresponding to lga intercept model */
  ind = 0;
  for(i=0;i<i_D_lga->glm->xcov[0]->n_cov;i++)
    {
      if(i_D_lga->glm->xcov[0]->in_landings[i]==1)
	{
	  D_totcatch->fac_lga[0][i] = i_inCatch->fac_lga_int[ind]-1;  /* Assume start on zero */
	  ind++;
	}
      else
	D_totcatch->fac_lga[0][i] = -1;
    }
  if(i_D_lga->glm->xcov[0]->icell>0)
    {
      D_totcatch->fac_lga[0][i_D_lga->glm->xcov[0]->icell] = i_inCatch->fac_lga_int[ind]-1;
    }

  #ifdef DEBUG_PREDICT
  printf("Factors corresponding to lga slope model\n");
  #endif
  /* Factors corresponding to lga slope model */
  ind = 0;
  for(i=0;i<i_D_lga->glm->xcov[1]->n_cov;i++)
    {
      if(i_D_lga->glm->xcov[1]->in_landings[i]==1)
	{
	  D_totcatch->fac_lga[1][i] = i_inCatch->fac_lga_slp[ind]-1; /* Assume start on zero */
	  ind++;
	}
      else 
	D_totcatch->fac_lga[1][i] = -1;
    }
  if(i_D_lga->glm->xcov[1]->icell>0)
    {
      D_totcatch->fac_lga[1][i_D_lga->glm->xcov[0]->icell] = i_inCatch->fac_lga_int[ind]-1;
    }
  #ifdef DEBUG_PREDICT
  for(ii=0;ii<i_D_lga->glm->nxcov;ii++)
    {
      for(i=0;i<i_D_lga->glm->xcov[ii]->n_cov;i++)
	fprintf(stderr,"fac_lga[%d][%d]=%d\n",ii,i,D_totcatch->fac_lga[ii][i]);
    }
  #endif

  ncov = 2;
  D_totcatch->fac_wgl = Mmatrix_2d(0,ncov-1,0,i_D_wgl->glm->xcov[0]->n_cov-1,sizeof(int),1);  // Free ok

  #ifdef DEBUG_PREDICT
  printf("Factors corresponding to wgl intercept model\n");
  #endif
  /* Factors corresponding to wgl intercept model */
  ind = 0;
  for(i=0;i<i_D_wgl->glm->xcov[0]->n_cov;i++) 
    {
      if(i_D_wgl->glm->xcov[0]->in_landings[i]==1)
	{
	  D_totcatch->fac_wgl[0][i] = i_inCatch->fac_wgl_int[ind]-1; /* Assume start on zero */
	  ind++;
	}
      else
	D_totcatch->fac_wgl[0][i] = -1;
    }
  if(i_D_wgl->glm->xcov[0]->icell>0)
    {
      D_totcatch->fac_wgl[0][i_D_wgl->glm->xcov[0]->icell] = i_inCatch->fac_wgl_int[ind]-1;
    }

  #ifdef DEBUG_PREDICT
  printf("Factors corresponding to wgl slope model\n");
  #endif
  /* Factors corresponding to wgl slope model */
  ind = 0;
  for(i=0;i<i_D_wgl->glm->xcov[1]->n_cov;i++)
    {
      if(i_D_wgl->glm->xcov[1]->in_landings[i]==1)
	{
	  D_totcatch->fac_wgl[1][i] = i_inCatch->fac_wgl_slp[ind]-1; /* Assume start on zero */
	  ind++;
	}
      else
	D_totcatch->fac_wgl[1][i] = -1;
    }
  if(i_D_wgl->glm->xcov[1]->icell>0)
    {
      D_totcatch->fac_wgl[1][i_D_wgl->glm->xcov[0]->icell] = i_inCatch->fac_wgl_int[ind]-1;
    }
  #ifdef DEBUG_PREDICT
  for(ii=0;ii<i_D_wgl->glm->nxcov;ii++)
    {
      for(i=0;i<i_D_wgl->glm->xcov[ii]->n_cov;i++)
	fprintf(stderr,"fac_wgl[%d][%d]=%d\n",ii,i,D_totcatch->fac_wgl[ii][i]);
    }
  #endif

  #ifdef DEBUG_PREDICT
  printf("Factors\n");
  #endif
  /* factors */
  D_totcatch->factors = Mmatrix_2d(0,D_totcatch->nCell-1,         // Free ok
				    0,D_totcatch->nFactors,sizeof(int),1);
  for(c=0;c<D_totcatch->nCell;c++)
    {
      for(i=0;i<D_totcatch->nFactors;i++)
	{
	  D_totcatch->factors[c][i] = (int) i_inCatch->factors[i][c]-1;
	}
      /* New hauls */
      D_totcatch->factors[c][D_totcatch->nFactors] = -1;
     }

  #ifdef DEBUG_PREDICT
  printf("total weight\n");
  #endif
  /* total weight */
  D_totcatch->catch = CALLOC(D_totcatch->nCell,double);           // Free ok
  for(c=0;c<D_totcatch->nCell;c++)
      D_totcatch->catch[c] = i_inCatch->catch[c];

  /* season */
  D_totcatch->season = CALLOC(D_totcatch->nCell,int);           // Free ok
  for(c=0;c<D_totcatch->nCell;c++)
    {
      D_totcatch->season[c] = i_inCatch->season[c];
    }

  #ifdef DEBUG_PREDICT
  printf("Convert factors for age model\n");
  #endif
  /* Convert factors for age model */
  D_totcatch->age_xcov = CALLOC(i_D_age->glm->nxcov,Data_cov *);   // Free ok
  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      D_totcatch->age_xcov[i] = CALLOC(i_D_age->glm->nxcov,Data_cov);   // Free ok
      err = convert_tot_cov(D_totcatch,D_totcatch->fac_age[0],i_D_age->glm->xcov[i],D_totcatch->age_xcov[i],1);
      if(err)
	{
	  write_warning("makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
    }
  // Must correct for haul effect - not in age model
  for(c=0;c<D_totcatch->nCell;c++)
    D_totcatch->age_xcov[0]->c_cov[c][i_D_age->glm->xcov[0]->ihaul]=-1;
  
  #ifdef DEBUG_PREDICT
  int j;
  fprintf(stderr,"D_totcatch->age_xcov[0]->c_cov:\n");
  for(c=0;c<MIN(D_totcatch->nCell,20);c++)
    {
      for(j=0;j<i_D_age->glm->xcov[0]->n_cov;j++)
	{
	  fprintf(stderr,"%d ",D_totcatch->age_xcov[0]->c_cov[c][j]);
	}
      fprintf(stderr,"\n");
    }
  #endif


  if(i_inc_hsz)
    {
      /* Convert factors for age model */
      D_totcatch->hsz_xcov = CALLOC(i_D_hsz->glm->nxcov,Data_cov *);   // Free ok
      for(i=0;i<i_D_hsz->glm->nxcov;i++)
	{
	  D_totcatch->hsz_xcov[i] = CALLOC(i_D_hsz->glm->nxcov,Data_cov);   // Free ok
	  err = convert_tot_cov(D_totcatch,D_totcatch->fac_hsz[0],i_D_hsz->glm->xcov[i],D_totcatch->hsz_xcov[i],0);
	  if(err)
	    {
	      write_warning("makedata_totcatch:Error calling convert_tot_cov\n");
	      return(err);
	    }
	}
    }

  #ifdef DEBUG_PREDICT
  printf("Convert factors for lga model\n");
  #endif
  /* Convert factors for lga model */
  D_totcatch->lga_xcov = CALLOC(i_D_lga->glm->nxcov,Data_cov *);    // Free ok
  for(i=0;i<i_D_lga->glm->nxcov;i++)  
    {
      D_totcatch->lga_xcov[i] = CALLOC(i_D_lga->glm->nxcov,Data_cov);    // Free ok
      if(i==0)
        err = convert_tot_cov(D_totcatch,D_totcatch->fac_lga[i],i_D_lga->glm->xcov[i],D_totcatch->lga_xcov[i],1);
      else
        err = convert_tot_cov(D_totcatch,D_totcatch->fac_lga[i],i_D_lga->glm->xcov[i],D_totcatch->lga_xcov[i],0);
      if(err)
	{
	  write_warning("makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
    }

  #ifdef DEBUG_PREDICT
  printf("Convert factors for wgl model\n");
  #endif
  /* Convert factors for wgl model */
  D_totcatch->wgl_xcov = CALLOC(i_D_wgl->glm->nxcov,Data_cov *);     // Free ok
  for(i=0;i<i_D_wgl->glm->nxcov;i++) 
    {
      D_totcatch->wgl_xcov[i] = CALLOC(i_D_wgl->glm->nxcov,Data_cov);     // Free ok
      if(i==0)
        err = convert_tot_cov(D_totcatch,D_totcatch->fac_wgl[i],i_D_wgl->glm->xcov[i],D_totcatch->wgl_xcov[i],1);
      else
        err = convert_tot_cov(D_totcatch,D_totcatch->fac_wgl[i],i_D_wgl->glm->xcov[i],D_totcatch->wgl_xcov[i],0);
      if(err)
	{
	  write_warning("makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
    }

  *o_D_totcatch = D_totcatch;

  return(0);
}		/* end of makedata_totcatch */


/*!
  \author Geir Storvik
  \brief Reallocate memory allocated by ::makedata_totcatch
*/
int re_makedata_totcatch(Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_wgl,Data_lin *i_D_hsz,
			 Input_totcatch *i_inCatch,int i_inc_hsz,
			 Data_totcatch **o_D_totcatch)
{
  int        err,i;
  Data_totcatch  *D_totcatch;
 
  D_totcatch = *o_D_totcatch;

  for(i=0;i<i_D_age->glm->nxcov;i++)
    {
      err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_age[0],i_D_age->glm->xcov[i],D_totcatch->age_xcov[i],0);
      if(err)
	{
	  write_warning("re_makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
      FREE(D_totcatch->age_xcov[i]);
    }
  FREE(D_totcatch->age_xcov);
  Fmatrix_2d(&D_totcatch->fac_age[0][0],&D_totcatch->fac_age[0]);

  if(i_inc_hsz)
    {
      err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_hsz[0],i_D_hsz->glm->xcov[i],D_totcatch->hsz_xcov[i],0);
      if(err)
	{
	  write_warning("re_makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
      FREE(D_totcatch->hsz_xcov[0]);
      FREE(D_totcatch->hsz_xcov);
      Fmatrix_2d(&D_totcatch->fac_hsz[0][0],&D_totcatch->fac_hsz[0]);
    }

  for(i=0;i<i_D_lga->glm->nxcov;i++)  
    {
      if(i==0)
        err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_lga[i],i_D_lga->glm->xcov[i],D_totcatch->lga_xcov[i],1);
      else
        err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_lga[i],i_D_lga->glm->xcov[i],D_totcatch->lga_xcov[i],0);
      if(err)
	{
	  write_warning("re_makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
      FREE(D_totcatch->lga_xcov[i]);
    }
  FREE(D_totcatch->lga_xcov);         
  Fmatrix_2d(&D_totcatch->fac_lga[0][0],&D_totcatch->fac_lga[0]);


  for(i=0;i<i_D_wgl->glm->nxcov;i++) 
    {
      if(i==0)
        err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_wgl[i],i_D_wgl->glm->xcov[i],D_totcatch->wgl_xcov[i],1);
      else
        err = re_convert_tot_cov(D_totcatch,D_totcatch->fac_wgl[i],i_D_wgl->glm->xcov[i],D_totcatch->wgl_xcov[i],0);
      if(err)
	{
	  write_warning("re_makedata_totcatch:Error calling convert_tot_cov\n");
	  return(err);
	}
      FREE(D_totcatch->wgl_xcov[i]);
    }
  FREE(D_totcatch->wgl_xcov);
  Fmatrix_2d(&D_totcatch->fac_wgl[0][0],&D_totcatch->fac_wgl[0]);


  Fmatrix_2d(&D_totcatch->factors[0][0],&D_totcatch->factors[0]);
  FREE(D_totcatch->catch);
  FREE(D_totcatch->season);

  FREE(D_totcatch);
  

  return(0);
}		/* end of re_makedata_totcatch */

    
/*!
  \author Geir Storvik
  \brief Converts covariates for total catch in the approperiate order.

  i_neff is the number of effects not being converted 
  (in order to not convert haul effects for linear models)

  Memory allocated by this routine is reallocated by ::re_convert_tot_cov
*/
static int convert_tot_cov(Data_totcatch *i_D_totcatch,int *i_fac,Data_cov *i_xcov,
			   Data_cov *o_xcov,int i_neff)
{
  int   c,j,k;

  /* Copy from i_xcov */
  o_xcov->n_cov = i_xcov->n_cov;
  o_xcov->n_fac = i_xcov->n_fac;
  o_xcov->fix = i_xcov->fix;
  o_xcov->ispat = i_xcov->ispat;
  o_xcov->icell = i_xcov->icell;
  o_xcov->iboat = i_xcov->iboat;
  o_xcov->num_adj_area = i_xcov->num_adj_area;
  o_xcov->adj_area = i_xcov->adj_area;

  /* Convert factors */
  o_xcov->c_cov = Mmatrix_2d(0,i_D_totcatch->nCell-1,     // Free ok
			     0,i_xcov->n_cov-1,sizeof(int),1);
  for(j=0;j<(i_xcov->n_cov-i_neff);j++)
    {
      #ifdef DEBUG_PREDICT
      fprintf(stderr,"j=%d, ncov=%d, i_neff=%d, icell=%d,ncell=%d,ifac=%d,ihaul=%d\n",j,i_xcov->n_cov,i_neff,i_xcov->icell,i_D_totcatch->nCell,i_fac[j],i_xcov->ihaul);
      #endif
      for(c=0;c<i_D_totcatch->nCell;c++)
	{
	  if(i_fac[j]==-1)
	    {
	      o_xcov->c_cov[c][j] = -1;
	    }
	  else
	    {
	      k = 0;
	      while(k<(i_xcov->n_fac[j]-1) && i_D_totcatch->factors[c][i_fac[j]]!=k)
		k++;
	      if(i_D_totcatch->factors[c][i_fac[j]]==k)
		{
		  o_xcov->c_cov[c][j] = k;
		}
	      else
		{
		  o_xcov->c_cov[c][j] = -1;
		}
	    }
	}
    }
  
  #ifdef DEBUG_PREDICT
  fprintf(stderr,"icell=%d,ispat=%d,ihaul=%d,iboat=%d\n",i_xcov->icell,i_xcov->ispat,i_xcov->ihaul,i_xcov->iboat);
  fprintf(stderr,"factors:\n");
  for(c=0;c<MIN(i_D_totcatch->nCell,20);c++)
    {
      for(j=0;j<i_xcov->n_cov;j++)
	{
	  fprintf(stderr,"%d ",o_xcov->c_cov[c][j]);
	}
      fprintf(stderr,"\n");
    }
  #endif
  
  return(0);
}		/* end of convert_tot_cov */

    
/*!
  \author Geir Storvik
  \brief Reallocate memory allocated bo ::convert_tot_cov
*/
static int re_convert_tot_cov(Data_totcatch *i_D_totcatch,int *i_fac,Data_cov *i_xcov,
			   Data_cov *o_xcov,int i_neff)
{
  Fmatrix_2d(&o_xcov->c_cov[0][0],&o_xcov->c_cov[0]);

  return(0);
}		/* end of re_convert_tot_cov */


/*!
  \author Geir Storvik
  \brief Convert covariates to be 0,1,...,

  Memory allocated by this routine is reallocated by ::re_convert_cov
*/
static int convert_cov(int i_nHaul,Data_cov *i_xcov)
{
  char    string[150];
  int     h,j,n;

  i_xcov->conv = CALLOC(i_xcov->n_cov,int *);   // Free ok

  for(j=0;j<i_xcov->n_cov;j++)
    {
      n = 0;
      for(h=0;h<i_nHaul;h++)
	{
	  n = MAX(n,i_xcov->c_cov[h][j]);
          i_xcov->c_cov[h][j]--;
	}
      if((i_xcov->fix[j]==0) || (j==i_xcov->ihaulsize)) // changed to include all random effects, not just area effect
	n = i_xcov->n_fac[j];
      if(n!= i_xcov->n_fac[j])
	{
	  write_warning("convert_cov:Something is wrong\n");
	  sprintf(string,"j=%d,n=%d,nFac=%d\n",j,n,(int) i_xcov->n_fac[j]);
          write_warning(string);
          for(h=0;h<i_nHaul;h++)
	    {
              sprintf(string,"c_cov[%d][%d]=%d\n",h,j,i_xcov->c_cov[h][j]);
              write_warning(string);
	    }
	  return(1);
	}
      i_xcov->conv[j] = CALLOC(n,int);   // Free ok
      for(n=0;n<i_xcov->n_fac[j];n++)
	{
	  i_xcov->conv[j][n] = n;
	}
    }

  return(0);
}		/* end of convert_cov */


/*!
  \author Geir Storvik
  \brief Reallocate memory allocated by ::convert_cov
*/
static int re_convert_cov(Data_cov *i_xcov)
{
  int j;

  for(j=0;j<i_xcov->n_cov;j++)
    FREE(i_xcov->conv[j]);
  FREE(i_xcov->conv);

  return(0);
}		/* end of re_convert_cov */


/*!
  \author Geir Storvik
  \brief  Compare two integers givin -1 if x<y and 1 otherwise

   Routine used for qsort.
*/
int compare(int *i_x,int *i_y)
{
 int r=0;
 if(*i_x < *i_y)
   r= -1;
 else if(*i_x > *i_y)
   r=1;
 return(r);
}		/* end of compare */


/*!
  \author Geir Storvik
  \brief  Compare two integers givin -1 if x<y and 1 otherwise

   Routine used for qsort.
*/
int compd(double *i_x,double *i_y)
{
 int r=0;
 if(*i_x < *i_y)
   r= -1;
 else if(*i_x > *i_y)
   r=1;
 return(r);
}		/* end of compare */



/*!
  \author Geir Storvik
  \brief Define node number of effects in graph structure

  \param i_xcov Struct containing covariates for the model
  \param i_nxcov The number of main terms in model 
         (currently 1 for age model, 2 for lga and wgl model)
  \param i_c_in_gr A matrix indicating if factor j of main term i
         is to be included in the graph
  \param o_node An array giving the node index in the graph
      - First index is for the main term (intercept or slope)
      - Second index is for the covariates (constant term, year, seas,...)
      - Third index is for the factor-level inside the covariate

  The node structure is as follows:
  - list For each main terms
     - sub Node for constant term
     - sub Node for year (if present in model)
     - sub Node for seas (if present in model)
     - sub Node for gear (if present in model)
     - sub Node for area (if present in model)
     - sub Node for cell (if present in model)
  For the age model there is only one main term. On the other hand the above
  structure is inside each age-category with all effects for the second age group
  comming after all effects for the first age group and so on.

  For the lga and wgl models there are two main terms, intercept and slope.

  Memory allocated by this routine is reallocated by ::re_find_node_effect
*/
int find_node_effect(Data_cov **i_xcov,int i_nxcov,int **i_c_in_gr,int ****o_node)
{
  int       i,j,k,ind,n_cov;
  int    ***node;
  Data_cov *xcov;

  node = CALLOC(i_nxcov,int **);     // Free ok
  ind = 0;
  for(i=0;i<i_nxcov;i++)
    {
      xcov = i_xcov[i];
      n_cov = xcov->n_cov;
      node[i] = CALLOC(n_cov,int *);     // Free ok

      for(j=0;j<xcov->n_cov;j++)
	{
          if(i_c_in_gr[i][j])
	    {
  	      node[i][j] = CALLOC(xcov->n_fac[j],int);   // Free ok
	      for(k=0;k<xcov->n_fac[j];k++)
		{
		  node[i][j][k] = ind;
		  ind++;
		}
	    }
	}
    }

  *o_node = node;

  return(0);
}		/* end of find_node_effect */


/*!
  \author Geir Storvik
  \brief Reallocate space allocated in ::find_node_effect
*/
int re_find_node_effect(Data_cov **i_xcov,int i_nxcov,int **i_c_in_gr,int ****o_node)
{
  int       i,j;
  int    ***node;
  Data_cov *xcov;

  node = *o_node;

  for(i=0;i<i_nxcov;i++)
    {
      xcov = i_xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          if(i_c_in_gr[i][j])
	    {
  	      FREE(node[i][j]);
	    }
	}
      FREE(node[i]);
    }
  FREE(node);


  return(0);
}		/* end of re_find_node_effect */


/*!
  \author Geir Storvik
  \brief Allocates space for Eff_str

  Memory allocated by this routine is reallocated by ::re_alloc_Eff_str
*/
int alloc_Eff_str(int i_ncat,int i_nxcov,Data_cov **i_xcov,Eff_str **o_par)
{
  int      a,i,j,k;
  Eff_str *par;

  par = CALLOC(1,Eff_str);         // Free ok
  par->eff = Mmatrix_2d(0,i_ncat-1,0,i_nxcov-1,sizeof(double **),1); //Free ok
  par->ssq = Mmatrix_2d(0,i_ncat-1,0,i_nxcov-1,sizeof(double *),1); //Free ok
  par->n_ssq = Mmatrix_2d(0,i_ncat-1,0,i_nxcov-1,sizeof(int *),1); //Free ok
  for(a=0;a<i_ncat;a++)
    for(i=0;i<i_nxcov;i++)
      {
	par->eff[a][i] = CALLOC(i_xcov[i]->n_cov,double *);  // Free ok
	for(j=0;j<i_xcov[i]->n_cov;j++)
	  {
	    par->eff[a][i][j] = CALLOC(i_xcov[i]->n_fac[j],double); // Free ok
	    for(k=0;k<i_xcov[i]->n_fac[j];k++)
	      par->eff[a][i][j][k] = G_ZERO;
	  }
	par->ssq[a][i] = CALLOC(i_xcov[i]->n_cov,double);  // Free ok
	par->n_ssq[a][i] = CALLOC(i_xcov[i]->n_cov,int);  // Free ok
      }
  
  par->ar = CALLOC(i_nxcov,double);    // Free ok
  par->prior_mean = Mmatrix_2d(0,i_ncat-1,0,i_nxcov-1,sizeof(double **),1); //Free ok
  par->prior_prec = CALLOC(i_nxcov,double **); 
  par->prior_ar = CALLOC(i_nxcov,double *); 
  par->tau = CALLOC(i_nxcov,double *); // Free ok
  for(a=0;a<i_ncat;a++)
    for(i=0;i<i_nxcov;i++)
      {
	par->prior_mean[a][i] = CALLOC(i_xcov[i]->n_cov,double *);  // Free ok
	for(j=0;j<i_xcov[i]->n_cov;j++)
	  {
	    par->prior_mean[a][i][j] = CALLOC(i_xcov[i]->n_fac[j],double); // Free ok
	  }
      }
  for(i=0;i<i_nxcov;i++)
    {
      par->prior_prec[i] = CALLOC(i_xcov[i]->n_cov,double *); 
      par->prior_ar[i] = CALLOC(2,double);  
      for(j=0;j<i_xcov[i]->n_cov;j++)
	{
	  par->prior_prec[i][j] = CALLOC(2,double); 
	}
      par->tau[i] = CALLOC(i_xcov[i]->n_cov,double);  // Free ok
    }
  par->prior_prec_obs = CALLOC(2,double);

  *o_par = par;
  return(0);
}		/* end of alloc_Eff_str */



/*!
  \author Geir Storvik
  \brief Reallocates memory allocated by ::alloc_Eff_str
*/
int re_alloc_Eff_str(int i_ncat,int i_nxcov,Data_cov **i_xcov,Eff_str **o_par)
{
  int      a,i,j;
  Eff_str *par;

  par = *o_par;

  for(a=0;a<i_ncat;a++)
    for(i=0;i<i_nxcov;i++)
      {
	for(j=0;j<i_xcov[i]->n_cov;j++)
	  FREE(par->eff[a][i][j]);
	FREE(par->eff[a][i]);
	FREE(par->ssq[a][i]);
	FREE(par->n_ssq[a][i]);
      }
  Fmatrix_2d(&par->eff[0][0],&par->eff[0]);
  Fmatrix_2d(&par->ssq[0][0],&par->ssq[0]);
  Fmatrix_2d(&par->n_ssq[0][0],&par->n_ssq[0]);
  
  FREE(par->ar);
  for(a=0;a<i_ncat;a++)
    for(i=0;i<i_nxcov;i++)
      for(j=0;j<i_xcov[i]->n_cov;j++)
	FREE(par->prior_mean[a][i][j]);
  for(i=0;i<i_nxcov;i++)
    {
      for(j=0;j<i_xcov[i]->n_cov;j++)
	FREE(par->prior_prec[i][j]); 
      FREE(par->prior_prec[i]);
      FREE(par->prior_ar[i]);
      FREE(par->tau[i]);
    }
  Fmatrix_2d(&par->prior_mean[0][0],&par->prior_mean[0]);
  FREE(par->prior_prec);
  FREE(par->prior_ar);
  FREE(par->prior_prec_obs);
  FREE(par->tau);
  FREE(par);

  return(0);
}		/* end of re_alloc_Eff_str */



/*!
  \author Geir Storvik
  \brief Calculated sum of effects for a given haul
*/
double calc_eff(Data_cov *i_xcov,double **i_eff,int i_h)
{
  int    j,k;
  double mu;

  mu = G_ZERO;
  for(j=0;j<i_xcov->n_cov;j++)
    {
      k = i_xcov->c_cov[i_h][j];
      mu += i_eff[j][k];
    }
  return(mu);
}


/*!
  \author Geir Storvik
  \brief Calculated sum of effects for a given haul
*/
double calc_eff_no_haul(Data_cov *i_xcov,double **i_eff,int i_h)
{
  int    j,k;
  double mu;

  mu = G_ZERO;
  for(j=0;j<(i_xcov->n_cov-1);j++)
    {
      k = i_xcov->c_cov[i_h][j];
      mu += i_eff[j][k];
    }
  return(mu);
}


/*!
  \author Hanne Rognebakke
  \brief Calculated sum of effects plus sufficient statistics for a given haul in the age model
*/
double calc_eff_suff_age(Data_cov *i_xcov,Eff_str *i_par,int i_h,int i_a)
{
  int    j,k;
  double mu;

  mu = G_ZERO;
  for(j=0;j<i_xcov->n_cov;j++)
    {
      k = i_xcov->c_cov[i_h][j];
      if(j==i_xcov->ihaulsize)
	{
	  mu += i_par->eff[i_a][0][j][k]*i_par->eff_hsz[i_h];
	}
      else
	mu += i_par->eff[i_a][0][j][k];
    }

  return(mu);
}


/*!
  \author Geir Storvik
  \brief Calculated log-density of univariate normal
*/
double ldnorm(double x,double mu,double sigma,double logsigma)
{
  double d,res;

  res = (x-mu)/sigma;

  d = -G_HALF*res*res-logsigma;

  return(d);
}		/* end of ldnorm */


/*!
  \author Geir Storvik
  \brief Updates average of simulations for age-parameters
*/
int update_average_age(int i_n,Age_struct *i_age,Data_age *i_D_age,
                       Age_struct *x_age_mean)
{
  int    a,h,err;
  
  err =  update_average_par(i_n,i_D_age->glm->ncat,i_age->par,x_age_mean->par,
                            i_D_age->glm->nxcov,i_D_age->glm->xcov);
  if(err)
    {
      write_warning("update_average_age:Error calling update_average_par\n");
      return(err);
    }

  /* Update haul effects */
  for(h=0;h<i_D_age->glm->nHaul;h++)
  for(a=0;a<i_D_age->glm->ncat;a++)
    update_mean(&(x_age_mean->alpha[h][a]),i_age->alpha[h][a],i_n);

  /* Update haul precision */
  update_mean(&(x_age_mean->par->tau_obs),i_age->par->tau_obs,i_n);


  return(err);
}		/* end of update_average_age */


/*!
  \author Hanne Rognebakke
  \brief Updates average of simulations for parameters in g_a function
*/
int update_average_g_a(int i_n,Data_g_a *i_D_g_a)
{
  int    a;
  /* Update g_a */
  for(a=0;a<i_D_g_a->ncat;a++)
    update_mean(&(i_D_g_a->g_a_mean[a]),i_D_g_a->g_a[a],i_n);
  for(a=0;a<i_D_g_a->g_a_npar;a++)
    update_mean(&(i_D_g_a->g_a_par_mean[a]),i_D_g_a->g_a_par[a],i_n);

  return(0);
}		/* end of update_average_g_a */


/*!
  \author Geir Storvik
  \brief Updates average of simulations for parameters in lga or wgl model
*/
int update_average_lin(int i_n,LW_struct *i_lin,Data_lin *i_D_lin,
                       LW_struct *x_lin_mean)
{
  int    err;
  
  err =  update_average_par(i_n,i_D_lin->glm->ncat,i_lin->par,x_lin_mean->par,
                            i_D_lin->glm->nxcov,i_D_lin->glm->xcov);
  if(err)
    {
      write_warning("update_average_lin:Error calling update_average_par\n");
      return(err);
    }

  /* Update fish precision */
  update_mean(&(x_lin_mean->par->tau_obs),i_lin->par->tau_obs,i_n);

  return(0);
}		/* end of update_average_lin */


/*!
  \author Geir Storvik
  \brief  Updates average for effects, precisions and ar-coefficients
*/
static int update_average_par(int i_n,int i_ncat,Eff_str *i_par,Eff_str *x_par_mean,
                              int i_nxcov,Data_cov **i_xcov)
{
  int    a,i,j,k;
  Data_cov *xcov;
  
  /* Update linear effects */
  for(a=0;a<i_ncat;a++)
  for(i=0;i<i_nxcov;i++)
    {
      xcov = i_xcov[i];
      for(j=0;j<xcov->n_cov;j++)
      for(k=0;k<xcov->n_fac[j];k++)
	update_mean(&(x_par_mean->eff[a][i][j][k]),i_par->eff[a][i][j][k],i_n);
    }

  /* Update precisions */
  for(i=0;i<i_nxcov;i++)
    {
      xcov = i_xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	update_mean(&(x_par_mean->tau[i][j]),i_par->tau[i][j],i_n);
    }
  /* Update ar-coef */
  for(i=0;i<i_nxcov;i++)
    update_mean(&(x_par_mean->ar[i]),i_par->ar[i],i_n);

  return(0);
}		/* end of update_average_par */


/*!
  \author Geir Storvik
  \brief Update the mean of a univariate quantity
*/
int update_mean(double *i_mean,double i_x,int i_n)
{
  double n,n1;
  n = (double) i_n;
  n1 = (double) (i_n+1);

  *i_mean = (n*(*i_mean)+i_x)/n1;

  return(0);
}		/* end of update_mean */


/*!
  \author Hanne Rognebakke
  \brief Writes the current MCMC samples of model1 to the binary file
*/
int write_samples_model1(Data_age *i_D_age, Age_struct *i_age,
			 Data_lin *i_D_lga, LW_struct *i_length, Data_g_a *i_D_g_a,
			 Data_lin *i_D_lga_CC, LW_struct *i_length_CC, Data_g_a *i_D_g_a_CC, 
			 Data_CC *i_D_CC, int i_coastal_cod, 
			 int i_print_boat, int i_print_format, Data_COST *i_D_COST)
{
  int i,err;

  /* Write age model */
  /* Don't write haul effect to mcmc-vector */
  i_D_age->glm->xcov[0]->n_cov--;
  if(i_print_format==0)
    err= write_it(g_caa_mcmc1,i_D_age->glm,i_age->par,0,i_print_boat);
  else
    err= write_it_ascii(g_caa_mcmc1,i_D_age->glm,i_age->par,0,i_print_boat);
  if(err)
    {
      write_warning("write_samples_model1:Error calling write_it, age parameters\n");
      return(err);
    }
  i_D_age->glm->xcov[0]->n_cov++;

  i_D_lga->glm->inc_hsz=0;
    
  /* Write lga model */
  if(i_print_format==0)
    err = write_it(g_caa_mcmc1,i_D_lga->glm,i_length->par,0,i_print_boat);
  else
    err = write_it_ascii(g_caa_mcmc1,i_D_lga->glm,i_length->par,0,i_print_boat);
  if(err)
    {
      write_warning("write_samples_model1:Error calling write_it, lga parameters\n");
      return(err);
    }
  /* Write g-function */
  if(i_D_g_a->g_a_npar>0)
    {
      if(i_print_format==0)
	fwrite(i_D_g_a->g_a_par,sizeof(double),i_D_g_a->g_a_npar,g_caa_mcmc1);
      else
	{
	  for(i=0;i<i_D_g_a->g_a_npar;i++)
	    fprintf(g_caa_mcmc1,"%f ",i_D_g_a->g_a_par[i]);
	  fprintf(g_caa_mcmc1,"\n");
	}
    }
  
  /* Write coastal cod */
  if(i_coastal_cod)
    {
      if(i_print_format==0)
	err = write_it(g_caa_mcmc1,i_D_lga_CC->glm,i_length_CC->par,i_D_CC->class_error,i_print_boat);
      else
	err = write_it_ascii(g_caa_mcmc1,i_D_lga_CC->glm,i_length_CC->par,i_D_CC->class_error,i_print_boat);
      if(err)
	{
	  write_warning("write_samples_model1:Error calling write_it, lga parameters\n");
	  return(err);
	}
      if(i_D_g_a_CC->g_a_npar>0)
	{
	  if(i_print_format==0)
	    fwrite(i_D_g_a_CC->g_a_par,sizeof(double),i_D_g_a_CC->g_a_npar,g_caa_mcmc1);
	  else
	    {
	      for(i=0;i<i_D_g_a->g_a_npar;i++)
		fprintf(g_caa_mcmc1,"%f ",i_D_g_a->g_a_par[i]);
	      fprintf(g_caa_mcmc1,"\n");
	    }
	}
      if(i_D_CC->class_error)
	{
	  if(i_print_format==0)
	    {
	      fwrite(&i_D_CC->k1,sizeof(double),1,g_caa_mcmc1);
	      fwrite(&i_D_CC->k2,sizeof(double),1,g_caa_mcmc1);
	    }
	  else
	    {
	      fprintf(g_caa_mcmc1,"%f %f\n",i_D_CC->k1,i_D_CC->k2);
	    }
	}
    }
  
  return(0);
}		/* end of write_samples_model1 */


/*!
  \author Hanne Rognebakke
  \brief Writes the current MCMC samples of model2 to the binary file 
*/
int write_samples_model2(Data_lin *i_D_wgl, LW_struct *i_weight,
			 Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC,
			 int i_coastal_cod, int i_print_boat, int i_print_format)
{
  int err=0;


  /* Write wgl model */
  if(i_print_format==0) // print binary format
    {
      err = write_it(g_caa_mcmc2,i_D_wgl->glm,i_weight->par,0,i_print_boat);
      if(err)
	{
	  write_warning("write_samples_model2:Error calling write_it, wgl parameters\n");
	  return(err);
	}
      if(i_coastal_cod)
	{
	  err = write_it(g_caa_mcmc2,i_D_wgl_CC->glm,i_weight_CC->par,0,i_print_boat);
	  if(err)
	    {
	      write_warning("write_samples_model2:Error calling write_it, wgl parameters\n");
	      return(err);
	    }
	}
    }
  else
    {
      err = write_it_ascii(g_caa_mcmc2,i_D_wgl->glm,i_weight->par,0,i_print_boat);
      if(err)
	{
	  write_warning("write_samples_model2:Error calling write_it, wgl parameters\n");
	  return(err);
	}
      if(i_coastal_cod)
	{
	  err = write_it_ascii(g_caa_mcmc2,i_D_wgl_CC->glm,i_weight_CC->par,0,i_print_boat);
	  if(err)
	    {
	      write_warning("write_samples_model2:Error calling write_it, wgl parameters\n");
	      return(err);
	    }
	}
    }

  
  return(0);
}		/* end of write_samples_model2 */



/*!
  \author Geir Storvik
  \brief Writes the current MCMC samples of a model to binary file

  For a specific model (age, lga or wgl) all parameters involved are written as
  a int string just after the previous iteration. 

  The format is
  - Effects: categories: main terms:covariates:factors
  - AR-coefficients: main terms
  - Precisions: Random effects and then observation precision
  - Likelihood

  For lga and wgl, there is just one category. For age there is only one main term, for
  lga and wgl there are two, intercept and slope. For age, observation precision correspond to
  haul precision.
*/
int write_it(FILE *fp, Data_glm *i_glm, Eff_str *i_par, int i_class_error, int i_print_boat)
{
  int       a,h,i,j,n,ncov;
  double   *eff_hsz=NULL;
  Data_cov *xcov;

  #ifdef WRITE_HAUL_EFF
  FILE   *unit;
  //unit = fopen("haul_lin.dat","w");
  unit = stderr;
  #endif

  #ifdef DEBUG_PROG
  fprintf(stderr,"nxcov=%d, ncat=%d, inc_hsz=%d\n",i_glm->nxcov,i_glm->ncat,i_glm->inc_hsz);
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
	  fprintf(stderr,"fac=%d\n",xcov->n_fac[j]);
	}
    }
  #endif

  
  
  //fprintf(stderr,"write_it:\n"); 
  if(i_glm->inc_hsz && (i_glm->ncat==1))
    {
      eff_hsz = CALLOC(i_glm->nHaul,double);    //FREE OK
      for(h=0;h<i_glm->nHaul;h++)
	{
	  eff_hsz[h] = i_glm->beta_hat[h][0][0]-calc_eff(i_glm->xcov[0],i_par->eff[0][0],h);
	}
      fwrite(eff_hsz,sizeof(double),i_glm->nHaul,g_caa_mcmc_hsz_eff);
      //fprintf(stderr,"eff_hsz=%f %f %f\n",eff_hsz[0],eff_hsz[1],eff_hsz[2]);
    }

  /* Write linear effects */
  #ifdef WRITE_HAUL_EFF
  fprintf(unit,"Effects\n");
  #endif
  n = 0;
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(a=0;a<i_glm->ncat;a++)
	{
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      if(j==xcov->iboat)
		{
		  if(i_print_boat)
		    {
		      fwrite(i_par->eff[a][i][j],sizeof(double),xcov->n_fac[j],fp);
		    }
		}
	      else
		fwrite(i_par->eff[a][i][j],sizeof(double),xcov->n_fac[j],fp);
	      //if(i_glm->ncat==1)
	      //fprintf(stderr,"write_it: eff[%d][%d][%d][%d]=%f\n",a,i,j,0,i_par->eff[a][i][j][0]);
	      n += xcov->n_fac[j];
              #ifdef WRITE_HAUL_EFF
	      for(k=0;k<xcov->n_fac[j];k++)
		fprintf(unit,"%d %d %d %d %f\n",i,a,j,k,i_par->eff[a][i][j][k]);
              #endif
	    }
	}
      if(i==0 && i_glm->inc_hsz && i_glm->ncat==1) // haulsize effects only for intercept wgl model (i=0 & ncat=1) when running main_model_hsz
	{
	  fwrite(eff_hsz,sizeof(double),i_glm->nHaul,fp);
	}
    }
  #ifdef DEBUG_PROG
  if(i_glm->ncat==1)
    fprintf(stderr,"print lin model:\n");
  else
    fprintf(stderr,"print age model:\n");
  fprintf(stderr,"printed %d effects\n",n);
  #endif

  /* Write ar-coef */
  #ifdef WRITE_HAUL_EFF
  fprintf(unit,"Ar\n");
  #endif
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      if(xcov->ispat > -1)
	{
	  fwrite(&i_par->ar[i],sizeof(double),1,fp);
	  n++;
          #ifdef WRITE_HAUL_EFF
          fprintf(unit,"%d %f\n",i_par->ar[i]);
          #endif
	}
    }
  #ifdef DEBUG_PROG
  fprintf(stderr,"printed %d param (incl ar-parameters)\n",n);
  #endif

  /* Write precisions for random effects*/
  #ifdef WRITE_HAUL_EFF
  fprintf(unit,"Precisions\n");
  #endif
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      ncov = xcov->n_cov;
      for(j=0;j<ncov;j++)
	{
          if(!xcov->fix[j])
	    {
	      fwrite(&i_par->tau[i][j],sizeof(double),1,fp);
	      n++;
              #ifdef WRITE_HAUL_EFF
	      fprintf(unit,"%d %d %f\n",i,j,i_par->tau[i][j]);
              #endif
	    }
	}
    }
  #ifdef DEBUG_PROG
  fprintf(stderr,"printed %d param (incl precisions)\n",n); 
  #endif

  /* Write observation precision */
  fwrite(&i_par->tau_obs,sizeof(double),1,fp);
  n++;
  #ifdef DEBUG_PROG
  fprintf(stderr,"printed %d param (incl observation precision)\n",n);
  fprintf(stderr,"Observation precisions\n");
  fprintf(stderr,"%f \n",i_par->tau_obs);
  #endif

  /* Write loglikelihood */
  fwrite(&i_par->loglik,sizeof(double),1,fp);
  n++;
  #ifdef DEBUG_PROG
  fprintf(stderr,"printed %d param (incl likelihood)\n",n);
  fprintf(stderr,"Loglikelihood\n");
  fprintf(stderr,"%f\n",i_par->loglik);
  #endif

  if(i_class_error)
    {
      if(n!=(i_par->num_var-2)) //parameters k1 and k2 written outside this routine
	{      
	  printf("write_it:n=%d != num_var=%d\n",n,i_par->num_var-2);
	  write_warning("write_it:Something is wrong\n");
	  return(1);
	}
    }
  else
    {
      if(n!=i_par->num_var)
	{      
	  printf("write_it:n=%d != num_var=%d\n",n,i_par->num_var);
	  write_warning("write_it:Something is wrong\n");
	  return(1);
	}
    }
  if(i_glm->inc_hsz && (i_glm->ncat==1))
    {
      FREE(eff_hsz);
    }

  #ifdef WRITE_HAUL_EFF
  //fclose(unit);
  #endif

  return(0);
}		/* end of write_it */



/*!
  \author Geir Storvik
  \brief Picks out parameters for iteration it from mcmc-vectors
*/
int read_it(FILE *fp, Data_glm *i_glm, Eff_str *i_par, int i_class_error, int i_read_boat)
{
  int       a,i,j,k,n,ncov,ret;
  double   *eff_hsz=NULL;
  Data_cov *xcov;

  if(i_glm->inc_hsz)
    eff_hsz = CALLOC(i_glm->nHaul,double);    //FREE OK
  

  /* Read linear effects */
  n = 0;
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(a=0;a<i_glm->ncat;a++)
	{
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      if(j==xcov->iboat)
		{
		  if(i_read_boat)
		    {
		      ret = fread(i_par->eff[a][i][j],sizeof(double),xcov->n_fac[j],fp);
		    }
		}
	      else
		ret = fread(i_par->eff[a][i][j],sizeof(double),xcov->n_fac[j],fp);
	      n += xcov->n_fac[j];
	      //if(a==0||i==2)
	      //fprintf(stderr,"read_it: eff[%d][%d][%d][%d]=%f\n",a,i,j,0,i_par->eff[a][i][j][0]);
	    }
	}
      if(i==0 && i_glm->inc_hsz && i_glm->ncat==1) // haulsize effects only for intercept wgl model (i=0 & ncat=1) when running main_model_hsz
	{
	  ret = fread(eff_hsz,sizeof(double),i_glm->nHaul,fp);
	  //fprintf(stderr,"read_it hsz: eff=%f %f %f\n",eff_hsz[0],eff_hsz[1],eff_hsz[2]);
	}
    }

  /* Read ar-coef */
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      if(xcov->ispat > -1)
	{
	  ret = fread(&i_par->ar[i],sizeof(double),1,fp);
	  n++;
	}
    }

  /* Read precisions for random effects*/
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      ncov = xcov->n_cov;
      for(j=0;j<ncov;j++)
	{
          if(!xcov->fix[j])
	    {
	      ret = fread(&i_par->tau[i][j],sizeof(double),1,fp);
	      n++;
	    }
	}
    }

  /* Read observation precision */
  ret = fread(&i_par->tau_obs,sizeof(double),1,fp);
  n++;

  /* Read loglikelihood */
  ret = fread(&i_par->loglik,sizeof(double),1,fp);
  n++;
  #ifdef DEBUG_PREDICT
  fprintf(stderr,"tau_obs=%f\n",i_par->tau_obs);
  fprintf(stderr,"loglik=%f\n",i_par->loglik);
  #endif
  
  if(i_class_error)
    {
      if(n!=(i_par->num_var-2)) //parameters k1 and k2 written outside this routine
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"read_it:n=%d != num_par=%d\n",n,i_par->num_var);
          #endif
	  printf("read_it:n=%d i_par->num_var=%d\n",n,i_par->num_var);
	  write_warning("read_it:Something is wrong\n");
	  return(1);
	}
    }
  else if(i_glm->inc_hsz)
    {
      //      printf("read_it:n=%d i_par->num_var=%d\n",n,i_par->num_var);
    }
  else
    {
      if(n!=i_par->num_var)
	{
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"read_it:n=%d != num_par=%d\n",n,i_par->num_var);
          #endif
	  printf("read_it:n=%d i_par->num_var=%d\n",n,i_par->num_var);
	  write_warning("read_it:Something is wrong\n");
	  return(1);
	}
    }
  if(i_glm->inc_hsz)
    {
      FREE(eff_hsz);
    }

  return(0);
}            /* End of read_it */



/*!
  \author Geir Storvik
  \brief Writes the current MCMC samples of a model to binary file

  For a specific model (age, lga or wgl) all parameters involved are written as
  a int string just after the previous iteration. 

  The format is
  - Effects: categories: main terms:covariates:factors
  - AR-coefficients: main terms
  - Precisions: Random effects and then observation precision
  - Likelihood

  For lga and wgl, there is just one category. For age there is only one main term, for
  lga and wgl there are two, intercept and slope. For age, observation precision correspond to
  haul precision.
*/
int write_it_ascii(FILE *fp, Data_glm *i_glm, Eff_str *i_par, int i_class_error, int i_print_boat)
{
  int       a,h,i,j,k,n,ncov;
  double   *eff_hsz=NULL;
  Data_cov *xcov;

  // Only for lin model (intercept wgl model (i=0 & ncat=1) when running main_model_hsz)
  if(i_glm->inc_hsz && i_glm->ncat==1)
    {
      eff_hsz = CALLOC(i_glm->nHaul,double);    //FREE OK
      for(h=0;h<i_glm->nHaul;h++)
	eff_hsz[h] = i_glm->beta_hat[h][0][0]-calc_eff(i_glm->xcov[0],i_par->eff[0][0],h);
    }

  if(i_glm->ncat==1)
    fprintf(fp,"\nPrint MCMC sample (lin model):\n");
  else
    fprintf(fp,"\nPrint MCMC sample (age model):\n");

  /* Write linear effects */
  n = 0;
  for(i=0;i<i_glm->nxcov;i++)
    {
      if(i_glm->ncat==1)
	{
	  if(i==0)
	    fprintf(fp,"Intercept:\n");
	  else 
	    fprintf(fp,"Slope:\n");
	}
      xcov = i_glm->xcov[i];
      for(a=0;a<i_glm->ncat;a++)
	{
	  if(i_glm->ncat>1)
	    fprintf(fp,"age group %d\n",a);
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      fprintf(fp,"Covariate %d:",j);
	      if(j==xcov->iboat)
		fprintf(fp," (boat)");
	      else if(j==xcov->ihaul)
		fprintf(fp," (haul)");
	      else if(j==xcov->ihaulsize)
		fprintf(fp," (haulsize)");
	      else if(j==xcov->icell)
		fprintf(fp," (cell)");
	      else if(j==xcov->ispat)
		fprintf(fp," (spatial)");
	      fprintf(fp,"\n");
	      if(j==xcov->iboat)
		{
		  if(i_print_boat)
		    {
		      for(k=0;k<xcov->n_fac[j];k++)
			{
			  fprintf(fp,"%f ",i_par->eff[a][i][j][k]);
			}
		      fprintf(fp,"\n");
		    }
		}
	      else
		{
		  for(k=0;k<xcov->n_fac[j];k++)
		    {
		      fprintf(fp,"%f ",i_par->eff[a][i][j][k]);
		    }		      
		  fprintf(fp,"\n");
		}
	      n += xcov->n_fac[j];
	    }
	}
      if(i==0 && i_glm->inc_hsz && i_glm->ncat==1) // haulsize effects only for intercept wgl model (i=0 & ncat=1) when running main_model_hsz
	{
	  fprintf(fp,"haulsize:\n");
	  for(h=0;h<i_glm->nHaul;h++)
	    {
	      fprintf(fp,"%f ",eff_hsz[h]);
	    }
	  fprintf(fp,"\n");
	}
    }
  #ifdef DEBUG_PROG
  if(i_glm->ncat==1)
    fprintf(stderr,"print lin model:\n");
  else
    fprintf(stderr,"print age model:\n");
  #endif
  fprintf(fp,"printed %d effects\n",n);

  /* Write ar-coef */
  fprintf(fp,"ar-coefficient:\n");
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      if(xcov->ispat > -1)
	{
	  fprintf(fp,"%f\n",i_par->ar[i]);
	  n++;
	}
    }

  /* Write precisions for random effects*/
  fprintf(fp,"precisions:\n");
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      ncov = xcov->n_cov;
      for(j=0;j<ncov;j++)
	{
          if(!xcov->fix[j])
	    {
	      fprintf(fp,"Covariate %d: ",j);
	      fprintf(fp,"%f ",i_par->tau[i][j]);
	      fprintf(fp,"\n");
	      n++;
	    }
	}
      //      fprintf(fp,"\n",i_par->tau[i][j]);
    }

  /* Write observation precision */
  fprintf(fp,"observation precision:\n");
  fprintf(fp,"%f\n",i_par->tau_obs);
  n++;

  /* Write loglikelihood */
  fprintf(fp,"log-likelihood:\n");
  fprintf(fp,"%f\n",i_par->loglik);
  n++;
  #ifdef DEBUG_PROG
  fprintf(stderr,"printed %d param (incl likelihood)\n",n);
  #endif

  if(i_class_error)
    {
      if(n!=(i_par->num_var-2)) //parameters k1 and k2 written outside this routine
	{      
	  printf("write_it:n=%d != num_var=%d\n",n,i_par->num_var-2);
	  write_warning("write_it:Something is wrong\n");
	  return(1);
	}
    }
  else
    {
      if(n!=i_par->num_var)
	{      
	  printf("write_it:n=%d != num_var=%d\n",n,i_par->num_var);
	  write_warning("write_it:Something is wrong\n");
	  return(1);
	}
    }

  if(i_glm->inc_hsz && i_glm->ncat==1)
    {
      FREE(eff_hsz);
    }

  return(0);
}		/* end of write_it_ascii */



/*!
  \author Hanne Rognebakke
  \brief Writes the current MCMC samples of model1 to the binary file
*/
int write_samples_haulsize(FILE *fp, Data_glm *i_glm, Eff_str *i_par)
{
  int       a,h,i,j,n;
  double   *eff_hsz;
  Data_cov *xcov;

  n = 0;
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(a=0;a<i_glm->ncat;a++)
	{
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      if(j==xcov->iboat)
		{
		  fwrite(i_par->eff[a][i][j],sizeof(double),xcov->n_fac[j],fp);
		  n += xcov->n_fac[j];
		  //printf("write_samples_haulsize: eff_hsz[0]=%f\n",i_par->eff[a][i][j][0]);
		}
	    }
	}
    }

  eff_hsz = CALLOC(i_glm->nHaul,double);    //FREE OK
  for(h=0;h<i_glm->nHaul;h++)
    {
      eff_hsz[h] = i_glm->beta_hat[h][0][0]-calc_eff(i_glm->xcov[0],i_par->eff[0][0],h);
    }
  fwrite(eff_hsz,sizeof(double),i_glm->nHaul,fp);
  //  printf("write_samples_haulsize: eff=%f %f %f\n",eff_hsz[0],eff_hsz[1],eff_hsz[2]);
  n += i_glm->nHaul;
  FREE(eff_hsz);

  fwrite(&i_par->tau_obs,sizeof(double),1,fp);
  n += 1;
  //  printf("write_samples_haulsize: tau=%f \n n=%d\n",i_par->tau_obs,n);


  return(0);
}		/* end of write_samples_haulsize*/




/*!
  \author Geir Storvik
  \brief Pick out parameters on the non-linear lga relation for iteration it
*/
int read_it_g_a(int i_it,Data_g_a *i_D_g_a)
{
  int  err=0,a,A;
  
  if(i_D_g_a->g_a_model==0)
    {
      A = i_D_g_a->ncat-1;
      for(a=0;a<i_D_g_a->ncat;a++)
	i_D_g_a->g_a[a] = 
	  (log(i_D_g_a->a_vec[a])-log(i_D_g_a->a_vec[0]))/
	  (log(i_D_g_a->a_vec[A])-log(i_D_g_a->a_vec[0]));
    }
  else if(i_D_g_a->g_a_model==1)
    {
      err = calc_g_a(i_D_g_a->ncat,i_D_g_a->a_vec,i_D_g_a->g_a_par,i_D_g_a->g_a);
      if(err){
	write_warning("read_it_g_a:Error calling calc_g_a\n");
	return(1);
      }
    }
  return(0);
}		/* end of read_it_g_a */



/*!
  \author Geir Storvik
  \brief Writes the current MCMC samples of catch at age to the mcmc-vector.

  The format is
  - Total catch: age-categories:length-categories
  - Mean lga and mean wgl for each age-category
*/
int write_it_totcatch(FILE *fp,int i_ncat,int i_nlen,
		      TC_struct *i_totcatch,double *i_mean_l,double *i_mean_w)
{
  int  a;

  for(a=0;a<i_ncat;a++)
    {
      fwrite(i_totcatch->catch_at_age[a],sizeof(double),i_nlen,fp);
    }

  fwrite(i_mean_l,sizeof(double),i_ncat,fp);
  fwrite(i_mean_w,sizeof(double),i_ncat,fp);
  
  return(0);
}		/* end of write_it_totcatch */


/*!
  \author Geir Storvik/Haavard Rue
  \brief propose a new value, scale, on the interval [1/f, f].

  Density of proposal is \f$\propto 1+1/x\f$.  This choice makes \f$q(x,x')/q(x',x) = 1\f$,
  in the acceptance ratio, when \f$x' = scale*x\f$.
*/
double scale_proposal(double x, double f, double *la)
{
    double len = f - 1/f;
    if (la) *la = 0.0;
    if (f == 1.0) return x;
    if ((*GMRFLib_uniform)() < len/(len+2*log(f)))
        return (1/f + len*(*GMRFLib_uniform)())*x;
    else
        return pow(f, 2.0*(*GMRFLib_uniform)()-1.0)*x;
}



/*!
  \author Geir Storvik
  \brief Generate an observation from the multinomial distribution 
  \param n Number of events that will be classified into one of
           the categories 1...ncat
  \param p Vector of probabilities.  p(i) is the probability that
           an event will be classified into category i.  Thus, p(i)
           must be [0,1]. Only the first ncat-1 p(i) must be defined
           since P(NCAT) is 1.0 minus the sum of the first
  \param ncat Number of categories.  Length of p and ix.
  \param ix Observation from multinomial distribution.  All ix(i)
            will be nonnegative and their sum will be n.

   Method: Algorithm from page 559 of Devroye, Luc, Non-Uniform Random Variate Generation.  
   Springer-Verlag, New York, 1986.

   This routine is a slight modification of the genmul routine of the ranlib library
*/
void my_genmul(long n,double *p,long ncat,long *ix)
{
static double prob,ptot,sum;
static long i,icat,ntot;
    ptot = 0.0F;
    for(i=0; i<ncat; i++) {
        ptot += *(p+i);
    }
/*
     Initialize variables
*/
    ntot = n;
    sum = 1.0F;
    for(i=0; i<ncat; i++) ix[i] = 0;
/*
     Generate the observation
*/
    for(icat=0; icat<(ncat-1); icat++) {
        prob = *(p+icat)/sum;
        *(ix+icat) = ignbin(ntot,prob);
        ntot -= *(ix+icat);
	if(ntot <= 0) return;
        sum -= *(p+icat);
    }
    *(ix+ncat-1) = ntot;
/*
     Finished
*/
    return;
}

