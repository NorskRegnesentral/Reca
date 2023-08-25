/*!
  \file caa.h
  \brief Main header file for the caa system
  \author Geir Storvik
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "GMRFLib.h"
#include "utl.h"
#include "ranlib.h"
#include "numrec.h"
#include "caa_substruct.h"

/* Universal constants */

/*
#define  REAL double
*/
#define  G_PI                3.14159265  /* ? */
#define  G_E                 2.718281    /* ? */
#define  G_ZERO              0.00000000
#define  G_ONE_THIRD         0.33333333 
#define  G_ONE_FOURTH        0.25000000 
#define  G_HALF              0.50000000
#define  G_ONE               1.00000000
#define  G_TWO               2.00000000
#define  G_THREE             3.00000000
#define  G_FOUR              4.00000000
#define  TRUE              1
#define  FALSE             0
#define  CALLOC(n,type) (type *)calloc((size_t)(n),sizeof(type))
#define  FREE(ptr) {if(ptr)free((void *)(ptr));ptr=NULL;}

//Priors or initial values for hyper parameters
#define P_TAU_CON            0.0001
#define P_TAU_FIX            0.001
#define P_TAU_RAN            100.0
#define P_TAU_OBS            1.0
#define P_AR                 0.5
#define P_GAMMA_A            0.01 //parameters in Gamma priors for precisions
#define P_GAMMA_B            0.01 //parameters in Gamma priors for precisions
#define P_BETA_A             1.0 //parameters in Beta priors
#define P_BETA_B             1.0 //parameters in Beta priors


//#define  LOG_FILE          1 //Writes output to log-file
//#define  DEBUG             1 //Different output for each iteration
//#define  DEBUG_TEST        1 //Testing under programming
//#define  DEBUG_INPUT       1 //Writes input from all main routines to files
//#define  DEBUG_GAUSS       0 //For debugging simulation of linear parts of model
//#define  DEBUG_GAUSS_FILE  0 //For debugging simulation of linear parts of model
//#define  DEBUG_GAUSS_CONSTR  0 //For debugging simulation of linear parts of model


FILE       *g_caa_log;    /*!< file-unit for log-file */
FILE       *g_caa_alpha;  /*!< file-unit for writing alpha's to file, for debugging */

FILE       *g_caa_mcmc1;  /*!< file-unit for fit parameters from age and lga model */
FILE       *g_caa_mcmc2;  /*!< file-unit for fit parameters from wgl model */
FILE       *g_caa_mcmc_hsz; /*!< file-unit for part of fit parameters from haulsize regression */
FILE       *g_caa_mcmc_hsz_eff; /*!< file-unit for haul effects from haulsize regression */

/*
#define MEMCHK(ptr) GMRFLib_ASSERT(ptr,GMRFLib_EMEMORY)
*/
/* Application specific constants: */

/* ran_unif01() */
/*
#define MULTIPLIER   69069
#define SHIFT            1
#define MODULUS      256*256*256*128
#define INVMOD       ( (REAL) 1 / ((REAL) MODULUS )) / ( (REAL) 2 )
*/

/* Macro definitions */
#define  MAX(A,B)  ((A) > (B) ? (A) : (B))
#define  MIN(A,B)  ((A) < (B) ? (A) : (B))
typedef  long int BOOL;

#define SWAP(a,b) tempr=a;a=b;b=tempr
#define  MAX_STR      1000


/*!
  \struct Age_struct caa.h
  \brief Define structure for age parameters

  A struct containing the model-description and the parameters in
  the age model. Typically names x_age where x_ is s_ if a global 
  variable, i_ if an input variable and o_ if an output variable.
*/
typedef struct
{
  Eff_str   *par;           /*!< \brief Parameters to be simulated  */
  double   **alpha;         /*!< \brief alpha[a][h] in age model */
  double    *mu;            /*!< \brief mu in Poison approximation of multinomial model */
  Graph_str *gr_str_f;      /*!< \brief Graph struct for effect */ 
  Graph2_str *gr_str_r;     /*!< \brief Graph struct for random effects, not in use now */
  int        age_errors;    /*!< \brief =1 if errors in age-readings, 0 otherwise */
  double   **A2A;           /*!< \brief Transition matrix for errors in age-readings */
  int       *A_Nneigh;      /*!< \brief Number of neighbors in A2A */
  int      **A_neigh;       /*!< \brief Neighbors in A2A */
  double     delta_age;     /*!< \brief Value to be added to all age groups when estimating alpha parameters */
  int       *model;         /*!< \brief indicators for model */
} Age_struct;

/*! 
  \struct LW_struct caa.h
  \brief Define structure for linear models lga and wgl

  A struct containing the model-description and the parameters in
  the  linear models lga or wgl. Typically named x_length and x_weight.
*/
typedef struct
{
  Eff_str   *par;           /*!< \brief Parameters to be simulated */
  Graph_str *gr_str;        /*!< \brief Graph struct for effects */ 
  int       *model;         /*!< \brief indicators for model */
  int        cens_model;    /*!< \brief = 1 if sampling discards, = 0 otherwise */
  double     cens_k;        /*!< \brief k-parameter in censoring */
  double     cens_m;        /*!< \brief m-parameter in censoring */
  double     cens_r;        /*!< \brief r-parameter in censoring */
  double     cens_Nlim;     /*!< \brief limit for number of sampled discards */
  int        fixed_model;   /*!< \brief = 0 for simulated parameters, = 1 for fixed parameters */
  double     *fixed_int;    /*!< \brief Values of fixed intercept parameter */
  double     *fixed_slp;    /*!< \brief Values of fixed slope parameter */
  double     *fixed_tau;    /*!< \brief Values of tau_obs parameters */
} LW_struct;

/*!
  \struct TC_struct caa.h
  \brief Define structure for catch-at-age

  A struct containing catch-by-age for the current simulation as well as
  vectors containing simulations from all iterations
*/
typedef struct
{
  double     **catch_at_age;/*! \brief catch_at_age for current simulation */
  double      *mcmc;        /*! \brief catch at age for all iterations */   
  double      *mean_l;      /*! \brief Mean in length given age for current simulation */
  double      *mean_w;      /*! \brief Mean in weight given age for current simulation */
  double      *planded;     /*!< \brief probability of landed (if COST) */
} TC_struct;

/*! \struct Data_orig
    \brief Contain the original age and length data of amigo type. 

    Needed when a non-linear lga model is used and sufficient statistics are to be 
    recalculated when non-linear model changes.
    This struct could perhaps be merged together with Data_l.
*/
typedef struct
{
  int          nFish;       /*!< \brief total number of fish (=length(totage)) */
  int          nHaul;       /*!< \brief Number of units (=length(nFishBoat)) */
  int         *samplingID;  /*!< \brief unit id of fish (=length(totage)) */
  int          coastal_cod; /*!< \brief =1 if including coastal cod, =0 otherwise */
  int          nBoat;       /*!< \brief Number of boats */
  int         *boat;        /*!< \brief Which boat the hauls correspond to */
  int         *nFishBoat;   /*!< \brief Number of fish pr unit (haul/trip) */
  int         *totage;      /*!< \brief ages of fish (age group) */     
  double      *part_year;   /*!< \brief season (part of year) of fish */     
  double      *totlength;   /*!< \brief lengths of fish */ 
  double      *lstart;      /*!< \brief lower limits of length-intervals  corresponding to totlength*/
  double      *lend;        /*!< \brief upper limits of length-intervals  corresponding to totlength*/
  double      *totweight;   /*!< \brief weights of fish */ 
  int         *replength;   /*!< \brief repetitions fish with equal lengths */        
  int         *start_noAge; /*!< \brief First position of missing age for haul */
  int         *start_Age;   /*!< \brief First position of observation for haul - only used for discards*/
  int         *num_noAge;   /*!< \brief Number of noAge fish in haul */
  int         *season;      /*!< \brief season for each fish */
  int         *mean_season; /*!< \brief mean season in a haul */
  double      *haulweight;  /*!< \brief haul sizes */
  int         *tottype;     /*!< \brief if cod, type of cod (coastal or skrei) */
  int          n_int_len;   /*!< \brief Number of intervals for lengths */
  double      *int_len_lim; /*!< \brief lower limits on intervals for lengths */
  double      *int_len;     /*!< \brief length value for intervals */
  int         *discard;     /*!< \brief repetitions with discarded */ 
  int         *n_discard;   /*!< \brief number of discards in each haul */
  int         *landed;      /*!< \brief repetitions with landed */
  int         *n_landed;    /*!< \brief number of landed in each haul */
} Data_orig;

/*!
  \struct Data_CC caa.h
  \brief Contain data for coastal cod
*/
typedef struct
{
  int          class_error;  /*!< \brief =1 if classification error, =0 otherwise */
  double      *ptype1_CC1;   /*!< \brief p(type==1|Coastal cod type1) */
  double      *ptype1_S5;    /*!< \brief p(type==1|Skrei type5) */
  double      *ptype2_CC2;   /*!< \brief p(type==2|Coastal cod type2) */
  double      *ptype2_S4;    /*!< \brief p(type==2|Skrei type4) */
  double      *ptype4_CC2;   /*!< \brief p(type==4|Coastal cod type2) */
  double      *ptype4_S4;    /*!< \brief p(type==4|Skrei type4) */
  double      *ptype5_CC1;   /*!< \brief p(type==5|Coastal cod type1) */
  double      *ptype5_S5;    /*!< \brief p(type==5|Skrei type5) */
  double       k1;           /*!< \brief constant defining the proportion of type1 CC */
  double       k2;           /*!< \brief constant defining the proportion of type5 Skrei */
} Data_CC;

/*! 
  \struct Data_age caa.h
  \brief Contain data and covariates for age model
*/
typedef struct
{
  Data_glm  *glm;          /*!< \brief covariates */
  int       *n_h;          /*!< \brief Number of fish in haul */
  double   **Ages;         /*!< \brief Number in each age group for each haul */
  double   **Ages_fix;     /*!< \brief Numbers fixed by observed ages */
  double   **Ages_disc;    /*!< \brief Number of discards in each age group for each haul */
  double   **Ages_land;    /*!< \brief Number of landed in each age group for each haul */
  int       *a_vec;        /*!< \brief Vector defining sequence of ages */
  double     delta_age;
  int       *type_age;        
  /*!< \brief Type of age observations
  = 0 if ages are given without errors,
  = 1 if ages are given with errors,
  = 2 if only lengths are given,
  = 3 if only lengths with some ages are given,
  = 4 if only lengths with some ages with errors are given 
  */
} Data_age;

/*! 
  \struct Data_lin caa.h
  \brief Contain data and covariates for linear model (lga and wgl)
*/
typedef struct
{
  
  Data_glm  *glm;              /*!< \brief covariates */
  int      **Ages;             /*!< \brief Number in each age group for each haul
			            Ages used in lga model, could be smaller than Ages 
			            in Data_age.  */
  int      **Ages_fix;         /*!< \brief Numbers fixed by observed ages */
  double   **Ages_miss_sim;    /*!< \brief Matrix to store simulated missing ages 
                                    and corresponding observed lengths*/
  double    *haulweight;       /*!< \brief weight in hauls */
  double   **sum_by_cat;       /*!< \brief For lga model, sum of length by age */
  double   **sum_by_cat_fix;   /*!< \brief Fixed by the aged observations */
  double   **sqsum_by_cat;     /*!< \brief For lga model, square sum of length by age */
  double   **sqsum_by_cat_fix; /*!< \brief Fixed by the aged observations */
} Data_lin;

/*!
  \struct Data_g_a caa.h
  \brief Contain parameters for g-function in lga model
*/
typedef struct
{
  double    *g_a;            /*!< \brief Function of a giving linear model between log-length and age */
  double    *g_a_mean;       /*!< \brief Mean values of g-function */
  int        g_a_model;      /*!< \brief  = 0 for log-linear model, = 1 for Schute-Richards model */
  int        ncat;           /*!< \brief Number of age intervals */
  int        nSeason;        /*!< \brief Number of seasons */
  double    *a_vec;          /*!< \brief Vector defining sequence of ages */   
  int       *a2Age_vec;      /*!< \brief Vector defining a_vec corresponding to a_vec in Data_age */     
  int        g_a_npar;       /*!< \brief Number of parameters describing g_a function */
  double    *g_a_par;        /*!< \brief Parameters describing g_a function */
  double    *g_a_par_mean;   /*!< \brief Mean values of g_a parameters */
  double   **suff;           /*!< \brief Sufficient statistics */  
  int        sample_c;       /*!< \brief =1 if c is sampled, 0 otherwise */
  int        sample_theta;   /*!< \brief =1 if theta is sampled, 0 otherwise */
  int        sample_gamma;   /*!< \brief =1 if gamma is sampled, 0 otherwise */
  double    *fixed_c;        /*!< \brief Values of fixed c parameter */
  double    *fixed_theta;    /*!< \brief Values of fixed theta parameter */
  double    *fixed_gamma;    /*!< \brief Values of fixed gamma parameter */
} Data_g_a;

/*! 
  \struct Data_l caa.h
  \brief Contain data and covariates for length-only and stratified by length data
*/
typedef struct
{
  /*! \brief Number of lengths in data, for non-Amigo boats */
  int     nLengths;      
  /*! \brief Integers corresponding to different hauls */
  int    *journey;       
  /*! \brief Number of fish in each length category */
  int    *count;         
  /*! \brief Length for each category  */
  double *cat;           
  /*! \brief Length for length category */
  double *length;        
  /*! \brief Number of length-categories for which fish are aged */
  int     nAgeLengths;   
  /*! \brief Lengths for which fish are aged */
  double *ageLength;     
  /*! \brief Number of aged fish in each length/age combination */
  int   **ageLengthCount;
  /*! \brief Corresponding journey */
  int    *ageJourney;    
  /*! \brief Number in each age category corresponding to length */
  int   **Ages;          
} Data_l;                     /* Only length data */


/*! 
  \struct Data_totcatch caa.h
  \brief Contain data and covariates for total catch

  To be used for prediction of catch-at-age
*/
typedef struct
{
  /*! \brief Number of cells with total catch */
  int          nCell;    
  /*! \brief Number of length intervals to divide predictions on */
  int          nlint;    
  /*! \brief Limits of length intervals */
  double      *l_int;    
  /*! \brief Number of factors */
  int          nFactors; 
  /*! \brief Factors corresponding to age */
  int        **fac_age;  
  /*! \brief Factors corresponding to lga */
  int        **fac_lga;  
  /*! \brief Factors corresponding to wgl */
  int        **fac_wgl;  
  /*! \brief Factors corresponding to hsz */
  int        **fac_hsz;  
  /*! \brief Factors corresponding to catch */
  int        **factors;  
  /*! \brief Catch  */
  double      *catch;
  /*! \brief  Season */
  int         *season;
  /*! \brief Factors for age model */
  Data_cov   **age_xcov; 
  /*! \brief Factors for lga model */
  Data_cov   **lga_xcov; 
  /*! \brief Factors for wgl model */
  Data_cov   **wgl_xcov; 
  /*! \brief Factors for hsz model */
  Data_cov   **hsz_xcov; 
} Data_totcatch;    /* Data for total catch */


/*!
  \struct Input_common caa.h
  \brief Contain common input parameters from R objects.
*/
typedef struct
{
  int          seed;        /*!< \brief Seed value */
  int          burn_in;     /*!< \brief Burn in */
  int          num_it_inner;/*!< \brief Number of iterations inside main loop */
  int          num_it_outer;/*!< \brief Main loop, number of iterations stored */
  int         *num_par1;    /*!< \brief Number of parameters in model1 (age,lga,ga,lga_CC,ga_CC)*/
  int         *num_par2;    /*!< \brief Number of parameters in model2 (wgl,wgl_CC)*/
  int          constr;      /*!< \brief Constr type, 1 for sum, 2 for treatment */
  int          sim_ar;
  int          use_debug;   /*!< \brief Print debug information if 1, not if 0 */
  char        *inputfolder; /*!< \brief Name of folder where input files are located */
  char        *filename_mcmc1;/*!< \brief Filename for printing binary files with age+lga mcmc results */
  char        *filename_mcmc2;/*!< \brief Filename for printing binary files with wgl mcmc results */
  char        *filename_hsz_hauleff;/*!< \brief Filename for printing binary files with haulsize effects */
  char        *filename_hsz_mcmc2;/*!< \brief Filename for printing binary files with hsz parameters mcmc results */
  char        *filename_hsz_it;/*!< \brief Filename for printing binary files with specific haulsize parameters */
  int          print_boat;  /*!< \brief 1 if print boat effects, 0 otherwise */
  int          inc_hsz;     /*!< \brief 1 if haulsize is included in the model, 0 otherwise */
  int          print_format;/*!< \brief format of output files: 0=binary, 1=ascii*/
} Input_common;

/*!
  \struct Input_cov caa.h
  \brief Contain covariate input parameters from R objects.
*/
typedef struct
{
  int          n_cov; /*!< \brief number of covariates */
  int         *n_lev; /*!< \brief  number of levels for each covariate */
  int         *random;
  int         *spatial;
  int         *continuous;
  int         *interaction;
  int         *in_landings;
  int         *in_slopeModel;
  int          n_cov_i; /*!< \brief number of integer covariates */
  int          n_cov_d; /*!< \brief number of double covariates */
  int        **c_cov_i; /*!< \brief  matrix of integer covariates */
  double     **c_cov_d; /*!< \brief  matrix of double covariates */
  double      *Sigma_cell;
  double      *constr_cell;
  int          nconstr_cell;
  int          ispat;
  int          iboat;
  int          icell;
  int          ihaul;  
  int          ihaulsize;
} Input_cov;


/*!
  \struct Input_age caa.h
  \brief Contain age input parameters from R objects.
*/
typedef struct
{
  int          nAges;
  int         *a_vec;
  double      *A2A;
  int          errors;
  int         *season;
  Input_cov   *cov;
  double       delta_age;
  int         *num_adj_area;
  int         *adj_area;
} Input_age;

/*!
  \struct Input_lga caa.h
  \brief Contain length-given-age input parameters from R objects.
*/
typedef struct
{
  int          g_a_model;
  int          g_a_ncat;
  int          g_a_nSeason;
  double      *g_a_par_init;
  int          g_a_sample_c;
  int          g_a_sample_theta;
  int          g_a_sample_gamma;
  int          fixed_model;
  double      *fixed_int;
  double      *fixed_slp;
  double      *fixed_tau;
  double      *fixed_g_a_c;
  double      *fixed_g_a_theta;
  double      *fixed_g_a_gamma;
  double      *fixed_int_CC;
  double      *fixed_slp_CC;
  double      *fixed_tau_CC;
  double      *fixed_g_a_c_CC;
  double      *fixed_g_a_theta_CC;
  double      *fixed_g_a_gamma_CC;
  int          cens_model;
  double      *cens_par;
  Input_cov   *int_cov;
  Input_cov   *slp_cov;
  int         *num_adj_area;
  int         *adj_area;  
} Input_lga;

/*!
  \struct Input_wgl caa.h
  \brief Contain weight-given-length input parameters from R objects.
*/
typedef struct
{
  int          fixed_model;
  double      *fixed_int;
  double      *fixed_slp;
  double      *fixed_tau;
  double      *fixed_int_CC;
  double      *fixed_slp_CC;
  double      *fixed_tau_CC;
  Input_cov   *int_cov;
  Input_cov   *slp_cov;
  int         *num_adj_area;
  int         *adj_area;  
} Input_wgl;

/*!
  \struct Input_prior caa.h
  \brief Contain prior input parameters from R objects.
*/
typedef struct
{
  double      *age_eff_mean;
  double      *age_eff_prec;
  double      *age_prec_par;
  double      *age_ar;
  double      *lga_eff_mean;
  double      *lga_eff_prec;
  double      *lga_prec_par;
  double      *lga_ar;
  double      *wgl_eff_mean;
  double      *wgl_eff_prec;
  double      *wgl_prec_par;
  double      *wgl_ar;
} Input_prior;

/*!
  \struct Input_cell caa.h
  \brief Contain cell input parameters from R objects.
*/
typedef struct
{
  int         *num_cell_o;
  int         *num_cell_u;
  int          age_int_nC;
  double      *age_int_E;
  double      *age_int_C;
  int          hsz_int_nC;
  double      *hsz_int_E;
  double      *hsz_int_C;
  int          lga_int_nC;
  double      *lga_int_E;
  double      *lga_int_C;
  int          lga_slp_nC;
  double      *lga_slp_E;
  double      *lga_slp_C;
  int          wgl_int_nC;
  double      *wgl_int_E;
  double      *wgl_int_C;
  int          wgl_slp_nC;
  double      *wgl_slp_E;
  double      *wgl_slp_C;
} Input_cell;

/*!
  \struct Input_totcatch caa.h
  \brief Contain total catch input parameters from R objects.
*/
typedef struct
{
  int         *season; 
  double      *midseason; 
  int          nCell;
  int          nFactors;
  int         *fac_age_int;
  int         *fac_hsz_int;
  int         *fac_lga_int;
  int         *fac_lga_slp;
  int         *fac_wgl_int;
  int         *fac_wgl_slp;
  int        **factors;
  double      *catch;
} Input_totcatch;

/*!
  \struct Input_predict caa.h
  \brief Contain common input parameters to predict from R objects.
*/
typedef struct
{
  int          nMCMC;
  int          nMCMC_hsz;
  int          burnin;
  int         *num_par1;
  int         *num_par2;
  int          nHaul;
  int          N_l_int;
  int         *n_MC;
  double      *l_int;
  int          coastal_cod;
  char        *inputfolder; /*!< \brief Name of folder where input files are located */
  char        *filename_mcmc1;/*!< \brief Filename for reading binary files with age+lga mcmc results */
  char        *filename_mcmc2;/*!< \brief Filename for reading binary files with wgl mcmc results */
  char        *filename_hsz_mcmc2;/*!< \brief Filename for reading binary files with hsz mcmc results */
  char        *filename_hsz_it;/*!< \brief Filename for printing binary files with specific haulsize parameters */
  char        *filename_predict;/*!< \brief Filename for printing binary files with predict mcmc results */
  int         *Npar1;
  int         *Npar2;
  int         *Npar3;
  int          split;
  int          read_boat;        /*!< \brief 1 if read boat effects, 0 otherwise */
  int          inc_hsz;     /*!< \brief 1 if haulsize is included in the model, 0 otherwise */
} Input_predict;

/*! 
  \struct Data_obs caa.h
  \brief Contain the original data for observer data

  Currently used in COST project
*/
typedef struct
{
  /*! \brief number of trips */
  int          n_trip;
  /*! \brief number of hauls pr trip */
  int         *num_trip;
  /*! \brief number of measured discarded fish pr haul */
  int         *num_haul_disc;
  /*! \brief observed month */
  int         *season;
  /*! \brief length categories for discard samples */
  double      *l_disc;
  /*! \brief number at length for discards */
  int         *lfreq_disc;
  /*! \brief number of discards in haul */
  double      *haulsize_disc;
  /*! \brief number of discards sampled */
  double      *sampsize_disc;
  /*! \brief numbers of discard age-length data in trip */
  int         *num_alk;
  /*! \brief ages for discard age-length data */
  int         *alk_a;
  /*! \brief lengths for discard age-length data */
  double      *alk_l;
  /*! \brief number at length for discard age-length data */
  int         *alk_lfreq;
  /*! \brief number of size classes pr trip with landings */
  int         *num_trip_land;
  /*! \brief number of measured landed fish pr size class */
  int         *num_size_land;
  /*! \brief length categories for landing samples */
  double      *l_land;
  /*! \brief number at length for landings */
  int         *lfreq_land;
  /*! \brief total weight landed in size class */
  double      *totsize_land;
  /*! \brief weight of landings sampled for lengths in size class */
  double     *sampsize_land;
} Data_obs; /* Observer data */

/*! 
  \struct cens_struct caa.h
  \brief Contain the parameters in the censoring model

  Currently used in COST project
*/
typedef struct
{
  int         ncat;          /*!< \brief number of categories for censoring parameters, e.g. n_trip */
  double      k;             /*!< \brief k-parameter in censoring (same for all hauls) */
  double      m;             /*!< \brief m-parameter in censoring (same for all hauls) */
  double     *r;             /*!< \brief r-parameter in censoring */
  double     *mu;            /*!< \brief expected value for censoring parameters */
  double     *tau;           /*!< \brief precision for censoring parameters */
  double      a_prior;       /*!< \brief parameter in gamma prior */ 
  double      b_prior;       /*!< \brief parameter in gamma prior */ 
  double      mu_prior_mean; /*!< \brief prior mean for expected value */
  double      mu_prior_prec; /*!< \brief prior precision for expected value */
} cens_struct;  /* Parameters in censoring model */

/*!
  \struct Data_mland caa.h
  \brief Contain the original data for market landings

  Currently used in COST project
*/
typedef struct
{
  /*! \brief number of trips */
  int          n_trip;
  /*! \brief number of size classes pr trip*/
  int         *num_trip;
  /*! \brief observed month */
  int         *season;
  /*! \brief number of measured fish pr size class */
  int         *num_size;
  /*! \brief length categories for landing samples */
  double      *l;
  /*! \brief number at length */
  int         *lfreq;
  /*! \brief total weight landed in size class */
  double      *totsize;
  /*! \brief weight of landings sampled for lengths in size class */
  double     *sampsize;
  /*! \brief numbers of age-length data in trip */
  int         *num_alk;
  /*! \brief ages for age-length data */
  int         *alk_a;
  /*! \brief lengths for age-length data */
  double      *alk_l;
  /*! \brief number at length for age-length data */
  int         *alk_lfreq;
  /*! \brief number of length categories for simulated discards */
  int          N_int_disc;
  /*! \brief length categories for simulated discards */
  double      *l_disc;
  /*! \brief number at length for simulated discards */
  int         *lfreq_disc;
  /*! \brief upper limits on intervals for lengths */
  double      *int_len_lim;
  /*! \brief parameter in Poisson distribution for number of fish */
  double      *lambda;
  /*! \brief hyperparameter for lambda */
  double       c;
  /*! \brief hyperparameter for lambda */
  double       d;
} Data_mland; /* Market landing data */

/*! 
  \struct Data_COST caa.h
  \brief Contain the original data with both observer data and market landing data
  and the simulated parameters in censoring model

  Currently used in COST project
*/
typedef struct
{
  /*! \brief =1 if COST model, =0 otherwise */
  int          model;
  /*! \brief number of iterations when resampling the length data */
  int          resample;
  /*! \brief Observer data */
  Data_obs    *obs;
  /*! \brief Market landing data */
  Data_mland  *mland;
  /*! \brief parameters in censoring model */
  cens_struct *cens;
  /*! \brief Number of COST specific variables to be simulated */
  int        num_var;
  /*! \brief MCMC samples of COST specific variables from all iterations */
  double      *mcmc; 
} Data_COST;  /* Original observer and market landing data */




