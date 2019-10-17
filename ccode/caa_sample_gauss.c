/*!
  \file caa_sample_gauss.c
  \author Geir Storvik
  \brief Routines for sampling linear models using the GMRFLib library

  For the age model this corresponds to simulating all linear coefficients
  and the precisions for the random terms given the \f$\alpha_{h,a}\f$'s.

  For the lga model this corresponds to simulating all linear coefficients
  and the precisions for the random terms given the \f$\alpha_{h,a}\f$'s.

  The variables to be simulated are stored in a graph with a structure as defined
  in the find_node_effect routine
*/
#include "caa.h"
#include "caa_sample_gauss.h"
#include "caa_read_write.h"
#include "caa_chol.h"
#include "caa_routines.h"


static int make_Q_inc_gauss(int i_nNodes,int i_ncat,int i_nxcov,int i_start_h,int i_nHaul,
		            Data_cov **i_xcov,int ***i_node,int **i_in_gr,int **o_Q_inc);
static double draw_spatial_ar(Data_cov *i_xcov,Eff_str *i_par,int i_ncat,int i_i);

#ifdef LOG_FILE
extern FILE     *g_caa_log;
#endif

#define NPROB 100 /*!< Number of possible values for ar-coef, see draw_spatial_ar */



/*!
  \brief Making graph for simulation from multivariate gaussian distribution

  The dependence structure is described by first making a Q_inc matrix (with
  the make_Q_inc_gauss routine) which has ones if nodes are neighbors and
  zero otherwise.

  Also include constraints. Only the constrast sum has been properly tested.
  \author Geir Storvik, UiO
*/
int make_graph_gauss(Data_glm *i_glm,Graph_str *i_gr,Eff_str *i_par,
                     int i_start_h,int i_constr)
{
  int  *no_data;
  int   a,err,i,j,n;
  int  *ind;
  int **Q_inc;
  Data_cov         *xcov;
  GMRFLib_graph_tp *graph;

  #ifdef DEBUG_PROG
  printf("Create graph\n");
  #endif
  err = GMRFLib_create_graph(&graph);
  if(err)
    {
      write_warning("make_graph_gauss:Problems with call to  GMRFLib_create_graph\n");
      return(err);
    }
  graph->n = 0;
  for(a=0;a<i_glm->ncat;a++)
    for(i=0;i<i_glm->nxcov;i++)
      {
	xcov = i_glm->xcov[i];
	for(j=0;j<xcov->n_cov;j++)
	  {
 	    if(i_gr->in_gr[i][j])
	      graph->n += xcov->n_fac[j];
	  }
      }
  ind = CALLOC(graph->n,int);    // Free ok

  //Q_inc = Mmatrix_2d(0,graph->n-1,0,graph->n-1,sizeof(int),1);   // Free ok
  //if(Q_inc==NULL)
  //  {
  //    write_warning("make_graph_gass:Error allocating memory for Q_inc\n");
  //    return(1);
  //  }
  Q_inc = CALLOC(graph->n,int *);
  for(i=0;i<graph->n;i++)
    Q_inc[i] = CALLOC(graph->n,int);
  err = make_Q_inc_gauss(graph->n,i_glm->ncat,i_glm->nxcov,i_start_h,i_glm->nHaul,
			 i_glm->xcov,i_gr->node,i_gr->in_gr,Q_inc);
  if(err)
    {
      write_warning("make_graph_gauss:Error calling make_Q_inc_gauss\n");
      return(err);
    }
  

  graph->nnbs = CALLOC(graph->n,int);      // Free ok
  graph->nbs = CALLOC(graph->n, int *);    // Free ok
  no_data = CALLOC(graph->n,int);      // Free ok

  for(i=0;i<graph->n;i++)
    {
      n = 0;
      for(j=0;j<graph->n;j++)
	{
          if(Q_inc[i][j])
	    {
	      ind[n] = j;
              n++;
	    }
	}
      graph->nnbs[i] = n;
      if(n>0)
	{
	  graph->nbs[i] = CALLOC(n,int);   // Free ok
	  for(j=0;j<n;j++)
	    graph->nbs[i][j] = ind[j];
	}
    }
  #ifdef DEBUG_PROG
  printf("make_graph_gauss:prepare graph\n");
  #endif
  err = GMRFLib_prepare_graph(graph);
  if(err)
    {
      write_warning("make_graph_gauss:Problems with call to  GMRFLib_prepare_graph\n");
      return(err);
    }
  i_gr->graph = graph;

  i_gr->b = CALLOC(i_gr->graph->n,double);     // Free ok
  i_gr->init = CALLOC(i_gr->graph->n,double);  // Free ok
  for(i=0;i<i_gr->graph->n;i++)
    i_gr->init[i] = GMRFLib_uniform();

  #ifdef DEBUG_PROG
  printf("make_graph_gauss: make constraints\n");
  #endif
  if(i_constr==1) // Currently used today
    {
      err = make_constr(i_gr,i_glm);
      if(err)
	{
	  write_warning("make_graph_gauss:Error calling make_constr\n");
	  return(err);
	}
    }
  else if(i_constr==2)
    {
      fprintf(stderr,"Warning: make_graph_gauss: constr=2 not tested in recent years\n");
      err = make_constr2(i_gr,i_glm);
      err = 1;
      if(err)
	{
	  write_warning("make_graph_gauss:Error calling make_constr2\n");
	  return(err);
	}
    }
  else
    {
      write_warning("make_graph_gauss:Unknown constraint\n");
      return(1);
    }

  i_gr->Q = Mmatrix_2d(0,i_gr->graph->n-1,     // Free ok
                       0,i_gr->graph->n-1,sizeof(double),1);

  #ifdef DEBUG_PROG
  printf("make_graph_gauss: make_Q_gauss\n");
  #endif
  err = make_Q_gauss(i_gr,i_glm,i_par,0,i_glm->nHaul,i_gr->Q,no_data);
  if(err)
    {
      write_warning("make_graph_gauss:Error calling make_Q_gauss\n");
      return(err);
    }
  #ifdef DEBUG_GAUSS_FILE
  FILE    *unit;
  unit = fopen("Q_gauss.dat","w");
  for(i=0;i<i_gr->graph->n;i++)
    {
      for(j=0;j<i_gr->graph->n;j++)
	fprintf(unit,"%lf ",i_gr->Q[i][j]);
      fprintf(unit,"\n");
    }
  fclose(unit);
  #endif

  #ifdef DEBUG_PROG
  fprintf(stderr,"make_graph_gauss: GMRFLib_init_problem\n");
  #endif
  err = GMRFLib_init_problem(
          &i_gr->model,i_gr->init,NULL,NULL,NULL,
          i_gr->graph,i_gr->Qfunc,(char *) i_gr->Q,NULL,i_gr->constr,0);
  if(err)
    {
      write_warning("make_graph_gauss:Error calling GMRFLib_init_problem\n");
      return(err);
    }
  #ifdef DEBUG_PROG
  fprintf(stderr,"make_graph_gauss: GMRFLib_init_problem OK!\n");
  #endif

  FREE(ind);
  //Fmatrix_2d(&Q_inc[0][0],&Q_inc[0]);
  for(i=0;i<graph->n;i++)
    FREE(Q_inc[i]);
  FREE(Q_inc);
  
  FREE(no_data);

  return(0);
}		/* end of make_graph_gauss */


/*!
  \author Geir Storvik
  \brief Reallocate memory allocated by make_graph_gauss
*/
int re_make_graph_gauss(Graph_str *i_gr)
{
  int   err,i,n;
  GMRFLib_graph_tp *graph;

  re_make_constr(i_gr);

  graph = i_gr->graph;

  for(i=0;i<graph->n;i++)
    {
      n = graph->nnbs[i];
      if(n>0)
	  FREE(graph->nbs[i]);
    }
  FREE(i_gr->b);
  FREE(i_gr->init);

  err = GMRFLib_free_graph(graph);
  if(err)
    {
      write_warning("re_make_graph_gauss:Problems with call to  GMRFLib_free_graph\n");
      return(err);
    }

  Fmatrix_2d(&i_gr->Q[0][0],&i_gr->Q[0]);

  err = GMRFLib_free_problem(i_gr->model);


  return(0);
}		/* end of re_make_graph_gauss */


/*!
  \author Geir Storvik
  \brief Makes a matrix Q_inc which is one for nodes which are neigbors.

  The neighbor structure is taking data into account.
*/
static int make_Q_inc_gauss(int i_nNodes,int i_ncat,int i_nxcov,int i_start_h,int i_nHaul,
			    Data_cov **i_xcov,int ***i_node,int **i_in_gr,int **o_Q_inc)
{
  int         a,h,i,i2,j,j2,k,ind,ind2,nArea,r1,r2,Nsub;
  Data_cov   *xcov, *xcov2;

  Nsub = i_nNodes/i_ncat;
  for(i=0;i<i_nNodes;i++)
    for(j=0;j<i_nNodes;j++)
      o_Q_inc[i][j] = 0;

  for(a=0;a<i_ncat;a++)
    for(i=0;i<i_nxcov;i++)
      {  
	xcov = i_xcov[i];
	/* Area effects */
	ind = xcov->ispat; 
	if(ind >= 0 && i_in_gr[i][ind])
	  {
	    nArea = xcov->n_fac[ind];
	    for(j=0;j<nArea;j++)
	      {
		r1 = a*Nsub + i_node[i][ind][j];
		for(k=0;k<xcov->num_adj_area[j];k++)
		  {
		    r2 = a*Nsub + i_node[i][ind][xcov->adj_area[j][k]];
		    o_Q_inc[r1][r2] = 1;
		  }
	      }
	  }
      }

  /* Interaction between effects because of observations */
     
  for(a=0;a<i_ncat;a++)
  for(i=0;i<i_nxcov;i++)
  for(i2=0;i2<i_nxcov;i2++)
    {
      xcov = i_xcov[i];
      xcov2 = i_xcov[i2];
      for(j=0;j<xcov->n_cov;j++)
      for(j2=0;j2<xcov2->n_cov;j2++)
	{
	  if(i_in_gr[i][j] && i_in_gr[i2][j2])
	    {
	      for(h=0;h<i_nHaul;h++)
		{
		  ind = a*Nsub + i_node[i][j][xcov->c_cov[h][j]];
		  ind2 = a*Nsub + i_node[i2][j2][xcov2->c_cov[h][j2]];
		  if(ind>i_nNodes || ind2>i_nNodes)
		    write_warning("make_Q_inc_gauss:Something is wrong\n");
		  o_Q_inc[ind][ind2] = 1;
		}
	    }
	}
    }
  for(i=0;i<i_nNodes;i++)
    o_Q_inc[i][i] = 0;
  /*
  unit = fopen("Q_inc_gauss.dat","w");
  for(i=0;i<i_nNodes;i++)
  for(j=0;j<i_nNodes;j++)
    fprintf(unit,"%d ",o_Q_inc[i][j]);
  fclose(unit);
  */
  return(0);
}		/* end of make_Q_inc_gauss */


/*!
  \author Geir Storvik
  \brief Calculates precision matrix for linear model 
*/
int make_Q_gauss(Graph_str *i_gr,Data_glm *i_glm,Eff_str *i_par,
		 int i_start_h,int i_nHaul,double **o_Q,int *o_no_data)
{
  int        a,i,i1,i2,j,j1,j2,k,isp,ind,ind1,ind2,h,Nsub;
  int        ihaulsize;
  double     ar,n_j,tau;
  Data_cov  *xcov, *xcov1, *xcov2;

  Nsub = i_gr->graph->n/i_glm->ncat;
  /* Initialize */  
  for(i=0;i<i_gr->graph->n;i++)
    {
      o_Q[i][i] = G_ZERO;
      for(j=0;j<i_gr->graph->nnbs[i];j++)
        o_Q[i][i_gr->graph->nbs[i][j]] = G_ZERO;
    }
  /* FIRST PRIOR */
  for(a=0;a<i_glm->ncat;a++)
    for(i=0;i<i_glm->nxcov;i++)
      {
	xcov = i_glm->xcov[i];
	isp = xcov->ispat;	
	for(j=0;j<xcov->n_cov;j++)
	  {
	    //if(i_gr->in_gr[i][j] && j!= isp && j!=xcov->icell)  (CHANGED BY HANNE 120118 )
	    if(i_gr->in_gr[i][j] && j!= isp)
	      {
		for(k=0;k<xcov->n_fac[j];k++)
		  {
		    ind = a*Nsub + i_gr->node[i][j][k];
		    o_Q[ind][ind] = i_par->tau[i][j];
		  }
	      }
	  }
	/* Spatial interaction */
	if(isp>=0 && i_gr->in_gr[i][isp])
	  {
	    ar = i_par->ar[i];
	    tau = i_par->tau[i][isp];
	    for(j=0;j<xcov->n_fac[isp];j++)
	      {
		n_j = (double) xcov->num_adj_area[j];
		ind = a*Nsub + i_gr->node[i][isp][j];
		//o_Q[ind][ind] = tau*n_j;
		o_Q[ind][ind] = tau*(ar*n_j+(G_ONE-ar));
		for(k=0;k<xcov->num_adj_area[j];k++)
		  {
		    ind2 = a*Nsub + i_gr->node[i][isp][xcov->adj_area[j][k]];
		    o_Q[ind][ind2] = -tau*ar;
		  }
	      }
	  }
	if(0)
	  {
	    printf("iboat=%d,i_gr=%d \n",xcov->iboat,i_gr->in_gr[i][xcov->iboat]);
	    tau = i_par->tau[i][xcov->iboat];	    
	    for(j=0;j<xcov->n_fac[xcov->iboat];j++)
	      {
		ind = a*Nsub + i_gr->node[i][xcov->iboat][j];
		if(j<3)
		  printf("Q[%d]=%f ",ind,o_Q[ind][ind]);
		o_Q[ind][ind] = tau;
		if(j<3)
		  printf("Q[%d]=%f ",ind,o_Q[ind][ind]);
	      }
	    printf("\n");
	  }
      }

  /* INCLUDE DATA */
  if(i_glm->xcov[0]->ihaulsize<0) //  Haulsize not included in the model
    {
      for(i1=0;i1<i_glm->nxcov;i1++)
	{
	  xcov1 = i_glm->xcov[i1];
	  for(j1=0;j1<xcov1->n_cov;j1++)
	    if(i_gr->in_gr[i1][j1])
	      {
		for(a=0;a<i_glm->ncat;a++)
		  for(h=i_start_h;h<i_nHaul;h++)
		    {
		      ind1 = a*Nsub + i_gr->node[i1][j1][xcov1->c_cov[h][j1]];
		      o_no_data[ind1] = 0;
		      for(i2=0;i2<i_glm->nxcov;i2++)
			{
			  xcov2 = i_glm->xcov[i2];
			  for(j2=0;j2<xcov2->n_cov;j2++)
			    {
			      if(i_gr->in_gr[i2][j2])
				{
				  ind2 = a*Nsub + i_gr->node[i2][j2][xcov2->c_cov[h][j2]];
				  o_Q[ind1][ind2] += i_par->tau_obs * i_glm->suff[h][i1][i2];
				}
			    }
			}
		    }
	      }
	}
    }
  else
    {
      ihaulsize = i_glm->xcov[0]->ihaulsize;
      for(i1=0;i1<i_glm->nxcov;i1++)
	{
	  xcov1 = i_glm->xcov[i1];
	  for(j1=0;j1<xcov1->n_cov;j1++)
	    if(i_gr->in_gr[i1][j1])
	      {
		for(a=0;a<i_glm->ncat;a++)
		  for(h=i_start_h;h<i_nHaul;h++)
		    {
		      ind1 = a*Nsub + i_gr->node[i1][j1][xcov1->c_cov[h][j1]];
		      o_no_data[ind1] = 0;
		      for(i2=0;i2<i_glm->nxcov;i2++)
			{
			  xcov2 = i_glm->xcov[i2];
			  for(j2=0;j2<xcov2->n_cov;j2++)
			    {
			      if(i_gr->in_gr[i2][j2])
				{
				  ind2 = a*Nsub + i_gr->node[i2][j2][xcov2->c_cov[h][j2]];
				  if(j1!=ihaulsize && j2!=ihaulsize)
				    {
				      o_Q[ind1][ind2] += i_par->tau_obs * i_glm->suff[h][i1][i2];
				    }
				  else if(j1==ihaulsize && j2!=ihaulsize)
				    {
				      o_Q[ind1][ind2] += i_par->tau_obs * i_glm->suff[h][i1][i2]*i_par->eff_hsz[h];
				    }
				  else if(j1!=ihaulsize && j2==ihaulsize)
				    {
				      o_Q[ind1][ind2] += i_par->tau_obs * i_glm->suff[h][i1][i2]*i_par->eff_hsz[h];
				    }
				  else if(j1==ihaulsize && j2==ihaulsize)   
				    {
				      o_Q[ind1][ind2] += i_par->tau_obs * i_glm->suff[h][i1][i2]*i_par->eff_hsz[h]*i_par->eff_hsz[h];
				    }
				  else
				    {
				      write_warning("make_Q_gauss:Something is wrong\n");				      
				    }
				  //	  if(a==1 && h<2)
				  // printf("a=%d,h=%d: j1=%d,j2=%d,tau*suff=%f,Q=%f,eff_hsz=%f\n",a,h,j1,j2,i_par->tau_obs * i_glm->suff[h][i1][i2],o_Q[ind1][ind2],i_par->eff_hsz[h]);
				}
			    }
			}
		    }
	      }
	}
    }

  return(0);
}		/* end of make_Q_gauss */


/*!
  \author Geir Storvik
  \brief  Construct sum constraints on model for GMRFLib

  Making constraints on all fixed effects to sum to zero. 

  In the current version also random effects are constrained to sum to zero.
  This is taken into account when sampling the precisions in the sample_precision
  routine.
  
*/
int make_constr(Graph_str *i_gr,Data_glm *i_glm)
{
  int                a,i,j,k,ind,n,n_constr,err,Nsub,ncat_max;
  int               *n_interact,**fac_interact;
  int             ***cov;
  Data_cov          *xcov;
  GMRFLib_constr_tp *constr;

  Nsub = i_gr->graph->n/i_glm->ncat;
  err = GMRFLib_create_constr(&constr);
  if(err)
    {
      write_warning("make_constr:Error calling GMRFLib_create_constr\n");
      return(err);
    }
  constr->nc = 0;

  
  /* Count constraints inside categories */
  ncat_max = max(1,i_glm->ncat-1);
  n_interact = CALLOC(i_glm->nxcov,int);
  fac_interact = CALLOC(i_glm->nxcov,int *);
  for(a=0;a<ncat_max;a++)
    {
      for(i=0;i<i_glm->nxcov;i++)
	{
	  xcov = i_glm->xcov[i];	  
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      //Constraints on fixed and random effects
	      //if(j==xcov->ispat && a==0)// if not include, also change below,see L644
	      //fprintf(stderr,"caa_sample_gauss.c: make_constr: Include constraint on spatial effect??\n");
	      //if(j!=xcov->ihaul && j!=xcov->iboat && j!=xcov->ispat)   // Not including haul effects, spatial or boat effects
	      if(j!=xcov->ihaul && j!=xcov->iboat)   // Not including haul effects and boat effects
		{
		  if(i_gr->in_gr[i][j] && xcov->n_fac[j]>1)
		    {
		      if(j==xcov->icell)
			{
			  n_interact[i]=0;
			  fac_interact[i] = CALLOC(xcov->n_cov,int);
			  for(k=0;k<xcov->n_cov;k++)
			    {
			      if(xcov->interaction[k]==1)
				{
				  fac_interact[i][n_interact[i]] = xcov->n_fac[k];
				  n_interact[i]++;
				  constr->nc = constr->nc+xcov->n_fac[k];
				}
			    }
			}
		      else
			{
			  constr->nc++;
			}
		    }
		}
	    }
	}
    }
  #ifdef DEBUG_PROG
  printf("Number of groups: %d\n",ncat_max);
  printf("Number of categories: %d\n",i_glm->nxcov);
  for(i=0;i<i_glm->nxcov;i++)
    {
      printf("\tNumber of covariates: %d\n",i_glm->xcov[i]->n_cov);
      for(j=0;j<i_glm->xcov[i]->n_cov;j++)
	printf("\tNumber of factors for covariate %d: %d\n",j+1,i_glm->xcov[i]->n_fac[j]);
    }
  printf("Number of constraints inside categories:%d\n",constr->nc);
  #endif

  if(i_glm->ncat>1) // age model
    {
      n=0;
      /* Count constraints between categories */
      for(i=0;i<i_glm->nxcov;i++)
	{
	  xcov = i_glm->xcov[i];
	  for(j=0;j<xcov->n_cov;j++) 
	    {
	      // Constraints on fixed and random effects - also including cell effects
	      if(i_gr->in_gr[i][j]) 
		{ 
		  constr->nc += xcov->n_fac[j];
		  n += xcov->n_fac[j];
		}
	    }
	}
      #ifdef DEBUG_PROG
      printf("Number of constraints between categories: %d\n",n);
      #endif
    }

  #ifdef DEBUG_PROG
  printf("Total number of constraints: %d\n",constr->nc);
  #endif

  constr->a_matrix = CALLOC(constr->nc * i_gr->graph->n,double);  // Free ok
  constr->e_vector = CALLOC(constr->nc,double);  // Free ok

  /* Find covariate matrix for cell effects */
  int ii,jj,ind_r,ind_f,nfac;
  cov = CALLOC(i_glm->nxcov,int **);
  for(i=0;i<i_glm->nxcov;i++) {
    xcov = i_glm->xcov[i];
    cov[i] = CALLOC(n_interact[i],int *);
    for(j=0;j<xcov->n_cov;j++) {
      if(j==xcov->icell)
	{
	  for(ii=0;ii<n_interact[i];ii++)
	    cov[i][ii] = CALLOC(xcov->n_fac[j],int);
	  ind_r = 0;
	  while(ind_r<xcov->n_fac[j])
	    {
	      for(k=0;k<fac_interact[i][0];k++)
		{
		  cov[i][0][ind_r] = k;
		  ind_r++;
		}
	    }
	  nfac = fac_interact[i][0];
	}
    }
    for(ii=1;ii<n_interact[i];ii++)
      {
	ind_r = 0;
	ind_f = 0;
	for(jj=0;jj<fac_interact[i][ii];jj++)
	  {
	    for(k=0;k<nfac;k++)
	      {
		cov[i][ii][ind_r] = ind_f;
		ind_r++;                    
	      }
	    ind_f++;
	  }
	nfac = nfac*fac_interact[i][ii];
      }
  }
  
  n_constr = 0;
  /* Constraints inside categories */
  for(a=0;a<ncat_max;a++) {
    for(i=0;i<i_glm->nxcov;i++)	{
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++) {
	//Constraints on fixed and random effects
	//if(j!=xcov->ihaul && j!=xcov->iboat && j!=xcov->ispat) {  // Not including haul effects - but include cell effects
	if(j!=xcov->ihaul && j!=xcov->iboat) {  // Not including haul effects - but include cell effects
	  if(i_gr->in_gr[i][j] && xcov->n_fac[j]>1) {
	    if(j==xcov->icell) {
	      for(ii=0;ii<n_interact[i];ii++) {
		for(jj=0;jj<fac_interact[i][ii];jj++)
		  {
		    for(k=0;k<xcov->n_fac[j];k++)
		      {
			if(cov[i][ii][k]==jj)
			  {
			    ind = a*Nsub + i_gr->node[i][j][k];
			    constr->a_matrix[constr->nc*ind+n_constr] = G_ONE;
			    //if(a==0)
			    // fprintf(stderr,"i=%d,j=%d,k=%d,ind=%d,ind_constr=%d,n_constr=%d\n",
			    //	      i,j,k,ind,constr->nc*ind+n_constr,n_constr);
			  }
		      }
		    constr->e_vector[n_constr] = G_ZERO;
		    n_constr++;
		  }
	      }
	    } else { // all other covariates
	      for(k=0;k<xcov->n_fac[j];k++) {
		ind = a*Nsub + i_gr->node[i][j][k];
		constr->a_matrix[constr->nc*ind+n_constr] = G_ONE;
		//if(a==0)
		// fprintf(stderr,"i=%d,j=%d,k=%d,ind=%d,ind_constr=%d,n_constr=%d\n",i,j,k,ind,constr->nc*ind+n_constr,n_constr);
	      }
	      constr->e_vector[n_constr] = G_ZERO;
	      n_constr++;
	    }
	  }
	}
      }
    }// for(i=0;i<i_glm->nxcov;i++)
  }
  for(i=0;i<i_glm->nxcov;i++)
    {
      for(ii=0;ii<n_interact[i];ii++)
	FREE(cov[i][ii]);
      FREE(cov[i]);
    }
  FREE(cov);
  
  if(i_glm->ncat>1) { // age model
    /* Count constraints between categories */
    for(i=0;i<i_glm->nxcov;i++)	{
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++) {
	// Constraints on fixed and random effects
	if(i_gr->in_gr[i][j]) {
	  for(k=0;k<xcov->n_fac[j];k++) {
	    for(a=0;a<i_glm->ncat;a++)
	      {
		ind = a*Nsub + i_gr->node[i][j][k];
		constr->a_matrix[constr->nc*ind+n_constr] = G_ONE;			      
	      }
	    constr->e_vector[n_constr] = G_ZERO;
	    n_constr++;
	  }
	}
      }
    }
  }


  if(n_constr != constr->nc)
    {
      fprintf(stderr,"n_const=%d,constr->nc=%d\n",n_constr,constr->nc);
      write_warning("make_constr:Not consistent number of constraints\n");
      return(1);
    }

  #ifdef DEBUG_GAUSS_CONSTR
  FILE   *unit; 
  if(i_glm->ncat>1)
    {
      fprintf(stderr,"Writing constraints\n");
      unit = fopen("constr_e.dat","w");
      for(i=0;i<constr->nc;i++)
	fprintf(unit,"%lf\n",constr->e_vector[i]);
      fclose(unit);
      unit = fopen("constr_a.dat","w");
      for(i=0;i<(constr->nc*i_gr->graph->n);i++)
	fprintf(unit,"%lf\n",constr->a_matrix[i]);
      fclose(unit);
      fprintf(stderr,"Writing constraints OK\n");
    }
  #endif

  err = GMRFLib_prepare_constr(constr,i_gr->graph,0);
  if(err)
    {
      write_warning("make_constr:Error calling GMRFLib_prepare_constr\n");
      return(err);
    }

  for(i=0;i<i_glm->nxcov;i++)
    FREE(fac_interact[i]);
  FREE(fac_interact);
  FREE(n_interact);
  
  i_gr->constr = constr;
  return(0);
}		/* end of make_constr */


/*!
  \author Geir Storvik
  \brief Reallocate space allocated by make_constr
*/
int re_make_constr(Graph_str *i_gr)
{
  int                err;
  GMRFLib_constr_tp *constr;

  constr = i_gr->constr;

  err = GMRFLib_free_constr(constr);

  return(0);
}		/* end of re_make_constr */


/*!
  \author Geir Storvik
  \brief Construct treatment constraints on model for GMRFLib

  Making constraints on first effect and last category to be zero.
  This constraint has not been properly tested!
*/
int make_constr2(Graph_str *i_gr,Data_glm *i_glm)
{
  int                a,i,j,k,ind,n_constr,err,Nsub;
  Data_cov          *xcov;
  GMRFLib_constr_tp *constr;

  Nsub = i_gr->graph->n/i_glm->ncat;
  err = GMRFLib_create_constr(&constr);
  if(err)
    {
      write_warning("make_constr2:Error calling GMRFLib_create_constr\n");
      return(err);
    }
  constr->nc = 0;

  /* Count constraints inside categories */
  for(a=0;a<i_glm->ncat;a++)
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          /* Only constraints on fixed effects with more than one category*/
          if(xcov->fix[j] && i_gr->in_gr[i][j] && xcov->n_fac[j]>1)  
	    constr->nc++;
	}
    }
    if(i_glm->ncat>1)
    {
      /* Count constraints between categories */
      for(i=0;i<i_glm->nxcov;i++)
	{
	  xcov = i_glm->xcov[i];
	  for(j=0;j<xcov->n_cov;j++)
	    if(xcov->fix[j] && i_gr->in_gr[i][j])
	      constr->nc += xcov->n_fac[j];
	}
    }
    constr->a_matrix = CALLOC(constr->nc * i_gr->graph->n,double);  // Free ok
    constr->e_vector = CALLOC(constr->nc,double);    // Free ok
  
  n_constr = 0;
  /* Constraints inside categories */
  for(a=0;a<i_glm->ncat;a++)
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          if(xcov->fix[j] && i_gr->in_gr[i][j] && xcov->n_fac[j]>1)  
	    { /* Only constraints on fixed effects with more than one category */
              k = 0;
	      ind = a*Nsub + i_gr->node[i][j][k];
	      constr->a_matrix[constr->nc*ind+n_constr] = G_ONE;
	      constr->e_vector[n_constr] = G_ZERO;
	      n_constr++;
	    }
	}
    }
  if(i_glm->ncat>1)
    {
      /* Count constraints between categories */
      for(i=0;i<i_glm->nxcov;i++)
	{
	  xcov = i_glm->xcov[i];
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      if(xcov->fix[j] && i_gr->in_gr[i][j])
		{
		  for(k=0;k<xcov->n_fac[j];k++)
		    {
                      a = i_glm->ncat-1;
		      ind = a*Nsub + i_gr->node[i][j][k];
		      constr->a_matrix[constr->nc*ind+n_constr] = G_ONE;
		      constr->e_vector[n_constr] = G_ZERO;
		      n_constr++;
		    }
		}
	    }
	}
    }

  if(n_constr != constr->nc)
    {
      write_warning("make_constr2:Not consistent number of constraints\n");
      return(1);
    }

  err = GMRFLib_prepare_constr(constr,i_gr->graph,0);
  if(err)
    {
      write_warning("make_constr2:Error calling GMRFLib_prepare_constr\n");
      return(err);
    }

  i_gr->constr = constr;
  return(0);
}		/* end of make_constr2 */



/*!
  \author Geir Storvik
  \brief Samples all effects in graph using GMRFLib_sample.

  First calculates precision matrix \f$Q\f$ through make_Q_gauss. Then calculates
  the \f$b\f$ vector defining the expectation through \f$\mu=Q^{-1}b\f$.
  The simulation is performed through cals to the GMRFLib_init_problem and
  GMRFLib_sample routines. Finally the samples are extracted from the graphs
  and stored in the approperiate structures.

  If LOG_FILE is defined, some of the simulated effects are written to the log-file
*/
int sample_gauss_eff(Graph_str *i_gr,Eff_str *i_par,Data_glm *i_glm,int i_start_h,int i_nHaul,int i_excl_haul)
{
  int        a,h,i,i2,j,k,keep,ind,err,Nsub;
  double     mu;
  Data_cov  *xcov;
  int       *no_data;

  no_data = CALLOC(i_gr->graph->n,int);
  for(i=0;i<i_gr->graph->n;i++)
    no_data[i] = 1;

  Nsub = i_gr->graph->n/i_glm->ncat;
  /* Calculate Q */
  err = make_Q_gauss(i_gr,i_glm,i_par,i_start_h,i_nHaul,i_gr->Q,no_data);
  if(err)
    {
      write_warning("sample_gauss_eff:Error calling make_Q_gauss\n");
      return(err);
    }


  /* Find b-vector */
  for(i=0;i<i_gr->graph->n;i++)
    i_gr->b[i] = 0;
  //Include prior mean on fixed effects
  //It is here assumed that fixed effects are independent a priori
  for(a=0;a<i_glm->ncat;a++)
    for(i=0;i<i_glm->nxcov;i++)
      {
	xcov = i_glm->xcov[i];
	for(j=0;j<xcov->n_cov;j++)
	  {
	    if(i_gr->in_gr[i][j] && xcov->fix[j])
	      {
		for(k=0;k<xcov->n_fac[j];k++)
		  {
		    ind = a*Nsub + i_gr->node[i][j][k];
		    i_gr->b[ind] = i_par->prior_mean[a][i][j][k] * i_par->tau[i][j];
		  }
	      }
	  }
      }
  //Data contribution
  if(i_glm->xcov[0]->ihaulsize<0) // Haulsize not included in model
    {
      for(a=0;a<i_glm->ncat;a++)
	for(h=i_start_h;h<i_nHaul;h++)
	  for(i=0;i<i_glm->nxcov;i++)
	    {
	      xcov = i_glm->xcov[i];
	      mu = G_ZERO;
	      for(i2=0;i2<i_glm->nxcov;i2++)
		{
		  mu += i_glm->beta_hat[h][a][i2]*i_glm->suff[h][i][i2];		
		}
	      for(j=0;j<xcov->n_cov;j++)
		{
		  if(i_gr->in_gr[i][j])
		    {
		      ind = a*Nsub + i_gr->node[i][j][xcov->c_cov[h][j]];
		      i_gr->b[ind] += i_par->tau_obs*mu;
		    }
		}
	    }
    }
  else
    {
      for(a=0;a<i_glm->ncat;a++)
	for(h=i_start_h;h<i_nHaul;h++)
	  for(i=0;i<i_glm->nxcov;i++)
	    {
	      xcov = i_glm->xcov[i];
	      mu = G_ZERO;
	      for(i2=0;i2<i_glm->nxcov;i2++)
		{
		  mu += i_glm->beta_hat[h][a][i2]*i_glm->suff[h][i][i2];		
		}
	      for(j=0;j<xcov->n_cov;j++)
		{
		  if(i_gr->in_gr[i][j])
		    {
		      ind = a*Nsub + i_gr->node[i][j][xcov->c_cov[h][j]];
		      if(j==xcov->ihaulsize)
			{
			  i_gr->b[ind] += i_par->tau_obs*mu*i_par->eff_hsz[h];
			}
		      else
			{
			  i_gr->b[ind] += i_par->tau_obs*mu;
			}
		    }
		}
	    }
    } 
  
  #ifdef DEBUG_GAUSS
  FILE  *unit;
  //  if(i_glm->ncat==1)
  if(i_glm->ncat>1)
    {
      printf("print debug Gauss\n");
      unit = fopen("gauss_b.dat","w");
      for(i=0;i<i_gr->model->n;i++)
	fprintf(unit,"%lf\n",i_gr->b[i]);
      fclose(unit);
      unit = fopen("gauss_Q.dat","w");
      printf("size Q graph=%d\n",i_gr->model->n);
      for(i=0;i<i_gr->model->n;i++)
	{
	  for(j=0;j<i_gr->model->n;j++)
	    fprintf(unit,"%lf ",i_gr->Q[i][j]);
	  fprintf(unit,"\n");
	}
      fclose(unit);
      unit = fopen("gauss_alpha.dat","w");
      for(h=0;h<i_nHaul;h++)
	{
	  for(a=0;a<i_glm->ncat;a++)
	    fprintf(unit,"%lf ",i_glm->beta_hat[h][a][0]);
	  fprintf(unit,"\n");
	}
      fclose(unit);
    }
  #endif

  /* Initialize problem */
  keep = GMRFLib_KEEP_graph;
  keep = GMRFLib_KEEP_graph+GMRFLib_KEEP_constr;
  err = GMRFLib_init_problem(&i_gr->model,i_gr->init,i_gr->b,NULL,NULL,i_gr->graph,
  			     i_gr->Qfunc,(char *) i_gr->Q,NULL,i_gr->constr,keep);
  if(err)
    {
      write_warning("sample_gauss_eff:Error calling GMRFLib_init_problem\n");
      return(err);
    }
      
  /* Perform sampling */
  err = GMRFLib_sample(i_gr->model);
  if(err)
    {
      write_warning("sample_gauss_eff:Error calling GMRFLib_sample\n");
      return(err);
    }
  if(!(i_gr->model->sample[0]> -999999999999.99 && i_gr->model->sample[0] < 99999999.99))
    {
      write_warning("sample_gauss_eff:Error calling GMRFLib_sample\n");
      err = 1;
      return(err);
    }

  /* Extract samples */
  if(i_excl_haul==1)
    {
      a=0; //only for lga model
      i_glm->xcov[0]->n_cov--;//haul effect stored last - will be sampled separately in sample_lga_haul()
      for(i=0;i<i_glm->nxcov;i++)
      {
	xcov = i_glm->xcov[i];
	for(j=0;j<xcov->n_cov;j++)
	  {
	    i_par->ssq[a][i][j] = G_ZERO;
	    i_par->n_ssq[a][i][j] = 0;
	    if(i_gr->in_gr[i][j])
	      {
		for(k=0;k<xcov->n_fac[j];k++)
		  {
		    ind = a*Nsub + i_gr->node[i][j][k];
		    if(i_start_h == 0 || no_data[ind]==0)
		      {
			i_par->eff[a][i][j][k] = i_gr->model->sample[ind];
			i_par->ssq[a][i][j] +=  i_par->eff[a][i][j][k]*i_par->eff[a][i][j][k];
			i_par->n_ssq[a][i][j] += 1;
		      }
		  }
	      }
	  }
      }
      i_glm->xcov[0]->n_cov++;
    }
  else
    {
      for(a=0;a<i_glm->ncat;a++)
      for(i=0;i<i_glm->nxcov;i++)
      {
	xcov = i_glm->xcov[i];
	for(j=0;j<xcov->n_cov;j++)
	  {
	    i_par->ssq[a][i][j] = G_ZERO;
	    i_par->n_ssq[a][i][j] = 0;
	    if(i_gr->in_gr[i][j])
	      {
		for(k=0;k<xcov->n_fac[j];k++)
		  {
		    ind = a*Nsub + i_gr->node[i][j][k];
		    if(i_start_h == 0 || no_data[ind]==0)
		      {
			i_par->eff[a][i][j][k] = i_gr->model->sample[ind];
			i_par->ssq[a][i][j] +=  i_par->eff[a][i][j][k]*i_par->eff[a][i][j][k];
			i_par->n_ssq[a][i][j] += 1;
		      }
		  }
	      }
	  }
      }
    }

  #ifdef DEBUG_CELL
  FILE   *unit_cell, *unit_eff;
  if(i_glm->ncat>1)
    unit_cell = fopen("cell_age.dat","w");
  else
    unit_cell = fopen("cell_lin.dat","w");
  i = 0;
  xcov = i_glm->xcov[i];
  j = xcov->icell;
  if(j>=0)
    {
      for(a=0;a<i_glm->ncat;a++)
      for(k=0;k<xcov->n_fac[j];k++)
	fprintf(unit_cell,"%lf\n",i_par->eff[a][i][j][k]);
    }
  fclose(unit_cell);
  if(i_glm->ncat>1)
    unit_eff = fopen("eff_age.dat","w");
  else
    unit_eff = fopen("eff_lin.dat","w");
  for(i=0;i<i_gr->graph->n;i++)
    fprintf(unit_eff,"%lf\n",i_gr->model->sample[i]);
  fclose(unit_eff);
  #endif

  FREE(no_data);

  return(0);
}		/* end of sample_gauss_eff */


/*! 
  \author Geir Storvik and Hanne Rognebakke
  \brief Samples precision parameters for random effects.

  Also samples ar-coefficients in spatial effects by calling draw_spatial_ar

  If LOG_FILE is defined, writes simulates values to log-file
*/
int sample_precision(int i_start_h,Eff_str *i_par,Data_glm *i_glm,int i_nHaul,int i_lgamod)
{
  int    i,isp,n,a,j,k,nfac,incl_cov;
  double ssq,ar,x,n_j;
  Data_cov *xcov;


  #ifdef LOG_FILE
  fprintf(g_caa_log,"sample_precision:nHaul=%d\n",i_nHaul);
  #endif

  for(i=0;i<i_glm->nxcov;i++)
    {  
      xcov = i_glm->xcov[i];
      isp = xcov->ispat;
      for(j=0;j<xcov->n_cov;j++)
	{
	  //if(!xcov->fix[j] && j != isp && j != xcov->icell) // cell and area sampled separately (CHANGED BY HANNE 120118 )
	  if(!xcov->fix[j] && j != isp) // area sampled separately
	    {
	      if(j != xcov->ihaul)
		nfac = xcov->n_fac[j];
	      else
		nfac = i_nHaul;
	      
	      if(i_glm->ncat>1){
		if(j== xcov->icell)// NB! n_fac_cell already subtracted 1!
		  n = xcov->n_fac_cell*(i_glm->ncat-1);
		else
		  n = (nfac-1)*(i_glm->ncat-1);
		//n = nfac*i_glm->ncat-(nfac+i_glm->ncat-1);
	      } 
	      else {
		if(j== xcov->icell)
		  n = xcov->n_fac_cell;  // NB! n_fac_cell already subtracted 1
		else
		  n = nfac-1;
	      }

	      ssq = G_ZERO;
	      for(a=0;a<i_glm->ncat;a++)
		{
		  for(k=0;k<nfac;k++)
		    {
		      incl_cov = 1;
		      if(i_lgamod==1 && j==xcov->ihaul) // if lga model and haul effect
                        {
                          if(i_glm->haul_obs[k]==0)
                            {
                              incl_cov = 0;  // not include haul effects for lga model if no observed values in that haul
                              n--;
                            }
                        }
		      if(i_lgamod==1 && j==xcov->iboat) // if lga model and boat effect
                        {
                          if(i_glm->boat_obs[k]==0)
                            {
                              incl_cov = 0;  // not include boat effects for lga model if no observed values in that boat
                              n--;
                            }
                        }
		      ssq += incl_cov*pow(i_par->eff[a][i][j][k],2);		    
		    }
		}
	      i_par->tau[i][j] = gengam(i_par->prior_prec[i][j][0]+G_HALF * ssq,
					i_par->prior_prec[i][j][1]+G_HALF * (double) n);//with prior
	      //if(j == xcov->icell)
	      //fprintf(stderr,"sample_precision cell:tau[%d][%d]=%f, ssq=%f,n=%d,nfac=%d,ncat=%d,prior[0]=%f,prior[1]=%f\n",
	      //	      i,j,i_par->tau[i][j],ssq,n,nfac,i_glm->ncat,i_par->prior_prec[i][j][0],i_par->prior_prec[i][j][1]);
	      //if(j == xcov->ihaul)
	      //fprintf(stderr,"sample_precision haul:tau[%d][%d]=%f, ssq=%f,n=%d\n",
	      //	i,j,i_par->tau[i][j],ssq,n);
              #ifdef LOG_FILE
	      fprintf(g_caa_log,"sample_precision:tau[%d][%d]=%f, ssq=%f,n=%d\n",
		      i,j,i_par->tau[i][j],ssq,n);
              #endif
	    } 
	} // end for(j=0;j<xcov->n_cov;j++)
      if(isp>=0 && !xcov->fix[isp]) 
	{
	  i_par->ar[i] = G_ZERO;
	  if(i_par->sim_ar)
	    i_par->ar[i] = draw_spatial_ar(xcov,i_par,i_glm->ncat,i);
	  if(i_glm->ncat>1)
	    n = xcov->n_fac[isp]*i_glm->ncat-(xcov->n_fac[isp]+i_glm->ncat-1);
	  else
	    n = xcov->n_fac[isp]-1;
  //	  printf("n=%d\n",n);
          /* Calculate ssq */
	  ssq = 0;
	  ar = i_par->ar[i];
	  for(a=0;a<i_glm->ncat;a++)
	    for(j=0;j<xcov->n_fac[isp];j++)
	      {
		x = i_par->eff[a][i][isp][j];
                n_j = (double) xcov->num_adj_area[j];
		//ssq += x*x*n_j;
		ssq += x*x*(n_j*ar+G_ONE-ar);
		for(k=0;k<xcov->num_adj_area[j];k++)
		    ssq -= x*i_par->eff[a][i][isp][xcov->adj_area[j][k]]*ar;
	      }
	  i_par->tau[i][isp] = gengam(i_par->prior_prec[i][isp][0]+G_HALF * ssq,
				      i_par->prior_prec[i][isp][1]+G_HALF * (double) n);//with prior
	  // fprintf(stderr,"sample_precision ar:tau_spat[%d]=%f,ar[%d]=%f: ssq=%f,n=%d,n_fac=%d,ncat=%d\n",
	  //  i,i_par->tau[i][isp],i,i_par->ar[i],ssq,n,xcov->n_fac[isp],i_glm->ncat);
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"sample_precision ar:tau_spat[%d]=%f,ar[%d]=%f: ssq=%f,n=%d,n_fac=%d,ncat=%d\n",
		  i,i_par->tau[i][isp],i,i_par->ar[i],ssq,n,xcov->n_fac[isp],i_glm->ncat);
          #endif
	}
      if(0){
      if(xcov->icell>=0)
	{
	  j = xcov->icell;
	  ssq = G_ZERO;
	  n = 0;
	  for(a=0;a<i_glm->ncat;a++)
	  for(k=0;k<xcov->n_fac[xcov->icell];k++)
	    {
	      xcov->cell_vec[n] = i_par->eff[a][i][xcov->icell][k];
	      n++;
	    }

	  //ssq = cholssq0(xcov->Qchol_cell,n,xcov->cell_vec);
	  ssq = symssq0(xcov->Q_cell,n,xcov->cell_vec);

	  n -= xcov->n_constr_cell;
 
	  i_par->tau[i][j] = gengam(i_par->prior_prec[i][j][0]+G_HALF * ssq,
				    i_par->prior_prec[i][j][1]+G_HALF * (double) n);//with prior
          #ifdef LOG_FILE
	  fprintf(g_caa_log,"sample_precision cell: tau_cell[%d][%d]=%f: ssq=%f,n=%d,n_constr=%d\n",
		 i,j,i_par->tau[i][j],ssq,n,xcov->n_constr_cell);
	  #endif
          #ifdef DEBUG_CELL
	  int c;
	  double sum; 
	  for(c=0;c<xcov->n_constr_cell;c++)
	    {
	      sum = G_ZERO;
	      for(k=0;k<n;k++)
		sum += xcov->constr_cell[c][k] * xcov->cell_vec[k];
	    }
	  if(fabs(sum)>0.00001) 
	    write_warning("sample_precision: Something is wrong\n");
	  fprintf(stderr,"cell, n_constr_cell=%d,sum_constr=%f,ssq=%f,n=%d\n",xcov->n_constr_cell,sum,ssq,n);
          #endif
	}
      }
    }

  return(0);
}		/* end of sample_precision */


/*!
  \author Geir Storvik
  \brief Samples ar-coef in spatial model

  Sampling is performed by calculating probabilities on a discrete number of outcomes
  between -1 and 1 specified by NPROB
*/
static double draw_spatial_ar(Data_cov *i_xcov,Eff_str *i_par,int i_ncat,int i_i)
{
  int    a,j,k,l,isp,err;
  double ar,prob[NPROB], **Q0,ssq,log_det_half,u,n,n_j,tau,x,scale_down=0.99;

  isp = i_xcov->ispat;
  tau = i_par->tau[i_i][isp];
  n = (double) (i_xcov->n_fac[isp]*i_ncat);
  Q0 = Mmatrix_2d(0,i_xcov->n_fac[isp]-1,                   // Free ok
		  0,i_xcov->n_fac[isp]-1,sizeof(double),1);
  
  /* Calculate probabilities */
  for(l=0;l<NPROB;l++)
    {
      ar = (double) ((l+1))/(double) (NPROB+1);
      ar *= scale_down;
      /* Calculate Q0 */
      for(j=0;j<i_xcov->n_fac[isp];j++)
	{
	  //Q0[j][j] = (double) i_xcov->num_adj_area[j];
	  Q0[j][j] = (double) i_xcov->num_adj_area[j]*ar + G_ONE-ar;
	  for(k=0;k<i_xcov->num_adj_area[j];k++)
	    Q0[j][i_xcov->adj_area[j][k]] = -ar;
	}
      /* Calculate determinant */
      err = choldc0(Q0,i_xcov->n_fac[isp]);
      if(err)
	{
	  write_warning("draw_spatial_ar:Error calling choldc0\n");
	  return(err);
	}
      log_det_half = G_ZERO;
      for(j=0;j<i_xcov->n_fac[isp];j++)
	{
	  log_det_half += log(Q0[j][j]);
	}
      /* Calculate ssq */
      ssq = 0;
      for(a=0;a<i_ncat;a++)
	for(j=0;j<i_xcov->n_fac[isp];j++)
	  {
	    x = i_par->eff[a][i_i][isp][j];
	    n_j = (double) i_xcov->num_adj_area[j];
	    ssq += x*x*(n_j*ar+G_ONE-ar);
	    for(k=0;k<i_xcov->num_adj_area[j];k++)
	      ssq -= x*i_par->eff[a][i_i][isp][i_xcov->adj_area[j][k]]*ar;
	  }
      prob[l] = exp(log_det_half * (double) i_ncat - G_HALF*ssq*tau); 
      if(l>0)
        prob[l] += prob[l-1];
    }
  /* Sample a */
  u = genunf(G_ZERO,prob[NPROB-1]);
  l = 0;
  while(l < (NPROB-1) && u>prob[l])
    l++;
  ar = (double) (l)/(double) (NPROB-1);
  ar = (double) (l+1)/(double) (NPROB+1);
  ar *= scale_down;
  

  Fmatrix_2d(&Q0[0][0],&Q0[0]);

  return(ar);
}		/* end of draw_spatial_ar */


/*!
  \author Geir Storvik
  \brief Samples precision parameter for fish (tau_obs) in lga or wgl model

  If LOG_FILE is defined, writes simulates values to log-file
*/
int sample_precision_lin(int i_start_h,Eff_str *i_par,Data_glm *i_glm,int i_nHaul)
{
  int    i,i2,n,h;
  double ssq,tot_ssq;
  Data_cov *xcov;
  double *mu;

  mu = CALLOC(i_glm->nxcov,double);    // Free ok

  #ifdef DEBUG_PROG
  printf("sample_precision_lin: use nHaul=%d\n",i_nHaul);
  #endif

  /* sample par->tau_obs */
  n = 0;
  tot_ssq = 0;
  for(h=i_start_h;h<i_nHaul;h++)
    {
      if(i_glm->suff[h][0][0]>0)
	{
	  n += i_glm->suff[h][0][0];
	  ssq = i_glm->ssq[h];
	  for(i=0;i<i_glm->nxcov;i++)
	    {
	      xcov = i_glm->xcov[i];
	      mu[i] = calc_eff(xcov,i_par->eff[0][i],h);
	    }
	  for(i=0;i<i_glm->nxcov;i++)
	    for(i2=0;i2<i_glm->nxcov;i2++)
	      {
		ssq += (mu[i]-i_glm->beta_hat[h][0][i])*
		  (mu[i2]-i_glm->beta_hat[h][0][i2])*i_glm->suff[h][i][i2];
		//		if(h<3)
		//printf("h=%d:i=%d,i2=%d,glm_ssq=%f,mu[%d]=%f,beta_hat=%f,mu[%d]=%f,beta_hat=%f,suff=%f,add=%f\n",
		//	h,i,i2,i_glm->ssq[h],i,mu[i],i_glm->beta_hat[h][0][i],i2,mu[i2],i_glm->beta_hat[h][0][i2],i_glm->suff[h][i][i2],
		//	 (mu[i]-i_glm->beta_hat[h][0][i])*(mu[i2]-i_glm->beta_hat[h][0][i2])*i_glm->suff[h][i][i2]);
	      }
	  //if(!(ssq > 0 && ssq < 9999999999999999999999.0))
	  if(!(ssq < 9999999999999999999999.0))
	    {	      
	      write_warning("sample_precision_lin:Something is wrong\n");
              #ifdef LOG_FILE
	      double dtmp;
	      fprintf(g_caa_log,"glm->ssq[%d]=%f\n",h,i_glm->ssq[h]);
	      for(i=0;i<i_glm->nxcov;i++)
		fprintf(g_caa_log,"mu[%d]=%f\n",i,mu[i]);
              fprintf(g_caa_log,"ssq = %lf\n",ssq);
	      for(i=0;i<i_glm->nxcov;i++)
		for(i2=0;i2<i_glm->nxcov;i2++)
		  {
		    fprintf(g_caa_log,"beta_hat[%d][0][%d]=%f, ",h,i,i_glm->beta_hat[h][0][i]);
		    fprintf(g_caa_log,"beta_hat[%d][0][%d]=%f, ",h,i2,i_glm->beta_hat[h][0][i2]);
		    fprintf(g_caa_log,"suff[%d][%d][%d]=%f\n",h,i,i2,i_glm->suff[h][i][i2]);
		    fprintf(g_caa_log,"%f  ",(mu[i]-i_glm->beta_hat[h][0][i]));
		    fprintf(g_caa_log,"%f\n",(mu[i2]-i_glm->beta_hat[h][0][i2]));
		    dtmp = (mu[i]-i_glm->beta_hat[h][0][i])*
		      (mu[i2]-i_glm->beta_hat[h][0][i2])*i_glm->suff[h][i][i2];
		    ssq += dtmp;
		    fprintf(g_caa_log,"ssq=%f,dtmp=%f\n",ssq,dtmp);
		  }
              #endif 
	      return(1);
	    }
	  if(!(ssq > 0))
	    {
              #ifdef LOG_FILE
	      double dtmp;
	      fprintf(g_caa_log,"glm->ssq[%d]=%f\n",h,i_glm->ssq[h]);
	      for(i=0;i<i_glm->nxcov;i++)
		fprintf(g_caa_log,"mu[%d]=%f\n",i,mu[i]);
              fprintf(g_caa_log,"ssq = %lf\n",ssq);
	      for(i=0;i<i_glm->nxcov;i++)
		for(i2=0;i2<i_glm->nxcov;i2++)
		  {
		    fprintf(g_caa_log,"beta_hat[%d][0][%d]=%f, ",h,i,i_glm->beta_hat[h][0][i]);
		    fprintf(g_caa_log,"beta_hat[%d][0][%d]=%f, ",h,i2,i_glm->beta_hat[h][0][i2]);
		    fprintf(g_caa_log,"suff[%d][%d][%d]=%f\n",h,i,i2,i_glm->suff[h][i][i2]);
		    fprintf(g_caa_log,"%f  ",(mu[i]-i_glm->beta_hat[h][0][i]));
		    fprintf(g_caa_log,"%f\n",(mu[i2]-i_glm->beta_hat[h][0][i2]));
		    dtmp = (mu[i]-i_glm->beta_hat[h][0][i])*
		      (mu[i2]-i_glm->beta_hat[h][0][i2])*i_glm->suff[h][i][i2];
		    ssq += dtmp;
		    fprintf(g_caa_log,"ssq=%f,dtmp=%f\n",ssq,dtmp);
		  }
              #endif
	    }
	  tot_ssq += ssq;
	}
    }
  i_par->tau_obs = gengam(i_par->prior_prec_obs[0]+G_HALF * tot_ssq,
			  i_par->prior_prec_obs[1]+G_HALF * (double) n);
  if(!(i_par->tau_obs > G_ZERO && i_par->tau_obs < 99999999999999.99))
    {
      write_warning("sample_precision_lin:Wrong tau_obs value\n");
      fprintf(stderr,"sample_precision_lin:tau_obs=%lf\n",i_par->tau_obs);
      #ifdef LOG_FILE
      fprintf(g_caa_log,"sample_precision_lin:tot_ssq=%lf,n=%d\n",tot_ssq,n);
      fprintf(g_caa_log,"sample_precision_lin:tau_obs=%lf\n",i_par->tau_obs);
      #endif
      return(1);
    }
  #ifdef LOG_FILE
  fprintf(g_caa_log,"sample_precision_lin:tot_ssq=%lf,n=%d\n",tot_ssq,n);
  fprintf(g_caa_log,"sample_precision_lin:tau_obs=%lf\n",i_par->tau_obs);
  #endif
  #ifdef DEBUG_PROG
  fprintf(stderr,"sample_precision_lin:tot_ssq=%lf,n=%d\n",tot_ssq,n);
  fprintf(stderr,"sample_precision_lin:tau_obs=%lf\n",i_par->tau_obs);
  #endif

  FREE(mu);
  return(0);
}		/* end of sample_precision_lin */



/*!
  \author Hanne Rognebakke
  \brief Samples lga haul effects 

*/
int sample_lga_haul(int i_start_h,Eff_str *i_par,Data_glm *i_glm)
{
  int        h,j,k,ixcov;
  double     mu,mu_h,b_h,var_h;
  Data_cov  *xcov;

  ixcov = i_glm->xcov[0]->ihaul; //index for which number the haul effect is stored
  i_par->ssq[0][0][ixcov] = G_ZERO;
  i_par->n_ssq[0][0][ixcov] = 0;

  for(h=i_start_h;h<i_glm->nHaul;h++)
    {
      // intercept
      xcov = i_glm->xcov[0];
      mu = G_ZERO;
      for(j=0;j<xcov->n_cov-1;j++) //not haul effect
	{
	  k = xcov->c_cov[h][j];
	  mu += i_par->eff[0][0][j][k];
	}
      b_h = (i_glm->beta_hat[h][0][0] - mu) * i_glm->suff[h][0][0];

      // slope
      mu = i_par->eff[0][1][0][0];
      b_h += (i_glm->beta_hat[h][0][1] - mu) * i_glm->suff[h][0][1];

      b_h = b_h * i_par->tau_obs;

      var_h = G_ONE/(i_par->tau[0][ixcov]+i_par->tau_obs*i_glm->suff[h][0][0]);
      mu_h = b_h * var_h;
      
      i_par->eff[0][0][ixcov][h] = gennor(mu_h,sqrt(var_h));
      i_par->ssq[0][0][ixcov] +=  i_par->eff[0][0][ixcov][h]*i_par->eff[0][0][ixcov][h];
      i_par->n_ssq[0][0][ixcov] += 1;
	
    }
  return(0);
}		/* end of sample_lga_haul */
