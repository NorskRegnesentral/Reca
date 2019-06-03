/*!
  \file caa_substruct.h
  \brief Structs included in other structs defined in caa.h
  \author Geir Storvik
*/

/*!
  \struct Eff_str caa.h
  \brief Struct for intercept or slope parameters

  A struct containing the current simulations of the effects (either in intercept 
  or slope coefficient) as well as precisions and spatial parameters. In addition,
  values from all simulations are stored in mcmc, see the routine write_it for 
  format of this.
  Typically a part of an Age_struct or a LW_struct object.
*/
typedef struct
{
  /*! \brief Array over all effects
    categories (first index) 
    x-covar (second index)
    c-covar (third index) 
    individual effects (forth index) 
  */
  double ****eff;
  double *eff_hsz;
  /*! \brief matrix of cell effects
    categories (first index) 
    individual effects (second index) 
  */
  double   **cell;
  /*! \brief Sum of squares for effects  */
  double  ***ssq;
  /*! \brief Number in sum of squares for effects  */
  int     ***n_ssq;
  /*! \brief Ar coefficient for each x-covar */
  double    *ar;
  int        sim_ar;
  /*! \brief Relative precisions for factors
    x-covar (first index)
    c-covar (second index)
    (assumed to be equal for different categories) 
  */
  double   **tau;
  /*! \brief Precision for observation - 
    corresponds to precision of haul effect in age model
   */
  double     tau_obs;
  double     tau_hsz;
  /*! \brief Prior mean of fixed effects */
  double ****prior_mean; 
  /*! \brief Parameters for gamma-prior on precisions */
  double  ***prior_prec; 
  /*! \brief Parameters for gamma-prior on obs(haul for age)-precisions */
  double    *prior_prec_obs;  
  /*! \brief Parameters for beta-prior on ar-coeff */
  double   **prior_ar;  
  /*! \brief Loglikelihood */
  double     loglik;
  /*! \brief Number of variables to be simulated */
  int        num_var;
} Eff_str;

/*!
  \struct Graph_str
  \brief A struct containing the graph structure for the linear part of a model

  To be used for simulation using the GMRFLib library 
  Typically a part of an Age_struct or a LW_struct object
*/
typedef struct
{
  int      **in_gr;          /*!< \brief =1 if factor is to be included in graph */
  int     ***node;           /*!< \brief Node number for effect */
  double   **Q;              /*!< \brief Q_matrix for parameters */
  double    *b;              /*!< \brief b vector in distribution for graph */
  GMRFLib_graph_tp *graph;   /*!< \brief Graph for model */
  GMRFLib_problem_tp *model; /*!< \brief Model for graph */
  GMRFLib_constr_tp *constr; /*!< \brief Constraints for model */
  double    *init;           /*!< \brief Initial values to GMRFLib_sample */
  double    *x_new;          
  GMRFLib_Qfunc_tp *Qfunc;   /*!< \brief Pick out Q-element */
  GMRFLib_Qfunc_tp *Qfunc_new;/*!< \brief Pick out Q-element for new hyperpar */
  double     logdens_x;     /*!< \brief log-dens of x given hyperpar */
  double     logdens_x_y;    /*!< \brief log-dens of x given y and hyperpar */
  double     logdens_y_x;    /*!< \brief log-dens of y given x and hyperpar */
} Graph_str;  /* For eff. in lin model and fix. eff. in age model */

/*!
  \struct Graph2_str caa.h
  \brief Not in use now.

  A struct for constructing graph for simulation of non-linear part of age
  model. 
*/
typedef struct
{
  
  int      **in_gr;   /*!< \brief =1 if factor is to be included in graph */
  
  int       *acc;     /*!< \brief accept rate */
  double    *N;       
  double    *mean;
  double    *mu_old;
  double    *mu_new;
  double    *x_new;
  double    *d;
  GMRFLib_graph_tp *graph;    /*!< \brief Graph for model for GMRFLib */
  GMRFLib_Qfunc_tp *Qfunc;    /*!< \brief Qfunction for GMRFLib Graph */
} Graph2_str;  /* For ran. eff. in age model */

/*!
  \struct Cell_dist
  \brief specifies conditional distribution for unobserved cells given
  observed cells

  Typically part of a Data_age or a Data_lin object. */
typedef struct
{
  /*! \brief Number of observed cells */
  int        n_o;            
  /*! \brief Number of unobserved cells */
  int        n_u;           
  /*! \brief Rank in conditional distribution */
  int        n_C;
  /*! \brief Expectation matrix */
  double   **E;
  /*! \brief "Cholesky" decomposition of covariance matrix  */
  double   **C; 
} Cell_dist;

/*!
  \struct Data_cov caa.h
  \brief Covariates for either the intercept or the slope structure of the linear models. 

  Typically a part of a Data_glm object. */
typedef struct
{
  /*! \brief Number of covariates */
  int        n_cov;            
  /*! \brief For each factor, number of categories */
  int       *n_fac; 
  /*! \brief If cell, number of categories */
  int        n_fac_cell; 
  /*! \brief =1 if factor is a fixed effect (i.e. presision is not estimated), 0 otherwise */
  int       *fix; 
  /*! \brief =1 if factor is an interaction effect (i.e. part of cell effect), 0 otherwise */
  int       *interaction; 
  /*! \brief =1 if factor exists in landings data (i.e. not boat effect), 0 otherwise */
  int       *in_landings; 
  /*! \brief Index indicating which factor has spatial structure */
  int        ispat;
  /*! \brief Index indicating which factor has haul structure */
  int        ihaul;
  /*! \brief Index indicating which factor has cell structure */
  int        icell;
  /*! \brief Index indicating which factor has boat structure */
  int        iboat;
  /*! \brief Index indicating which factor has haulsize haul structure */
  int        ihaulsize;
  /*! \brief Prior precision matrix for observed cell effects */
  double   **Q_cell;
  /*! \brief Chol of Q_cell */
  double   **Qchol_cell;
  /*! \brief Vector of cell effects */
  double    *cell_vec;
  /*! \brief Number of constraints on observed cell effects */
  int        n_constr_cell;
  /*! \brief Constraints on observed cell effects */
  double   **constr_cell;
  /*! \brief categorial covariates */
  int      **c_cov;            
  /*! \brief Table for converting factors to 0,1,... */
  int      **conv;             
  /*! \brief Number of neighbors for each area */
  int       *num_adj_area;     
  /*! \brief List of neighbors */
  int      **adj_area;  
  /*! \brief cell distribution for submodels */
  Cell_dist *cell_dist;
} Data_cov;

/*!
  \struct Data_glm
  \brief Organize covariates and containing sufficient statistics for the linear models. 

  Typically part of a Data_age or a Data_lin object. */
typedef struct
{
  /*! \brief Number of parameters in model */
  int        numpar;
  /*! \brief Number of categories, nAges for age, 1 for length and weight */
  int        ncat;            
  /*! \brief Number of nHaul */
  int        nHaul;           
  /*! \brief Number of nHaul for presplit runs */
  int        nHaulOld;           
  /*!< \brief Indicator if haul includes observed ages, 1 if observed, 0 otherwise */
  int       *haul_obs;      
  /*!< \brief Indicator if boat includes observed ages, 1 if observed, 0 otherwise */  
  int       *boat_obs;        
  /*!< \brief Indicator if haulsize included in age model, 1 if included, 0 otherwise */  
  int        inc_hsz;
  /*! \brief Number indexes */
  int        nindex;
  /*! \brief index for each haul */
  int       *index;
  /*! \brief Number of continuous covariates (including intercept) */
  int        nxcov; 
  /*! \brief xcov[0] contains intercept covariates, xcov[1] contains slope covariates */
  Data_cov **xcov;
  /*! \brief Sufficient statistics

  suff[h][0][0] = number of observations 
  suff[h][0][1] = sum, sum of log(predictor)
  suff[h][1][1] = sum2, sum of log(predictor)^2 */
  double  ***suff;        
  /*! \brief beta_hat[h][a][0] = intercept of regression, beta_hat[h][a][1] = slope */
  double  ***beta_hat;
  /*! \brief Square error of regression */
  double    *ssq;    
  /*! \brief Sufficient statistics based on observed ages - fixed during simulations */
  double  ***suff_fix;        
  /*! \brief Regression parameters based on observed ages - fixed during simulations */
  double  ***beta_hat_fix;
  /*! \brief Square error of regression based on observed ages - fixed during simulations */
  double    *ssq_fix;    
} Data_glm;

