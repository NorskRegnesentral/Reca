/*!
  \file caa_read_binary.c
  \brief Routines for reading binary files and put data into structures.
*/
#include "caa.h"
#include "caa_read_write.h"


static int read_bin_var(FILE *fp, char **o_varname, int *o_nvec, int **o_ivec,
			double **o_dvec, char **o_cvec, int *o_nchar, int *o_type);

/*!
\author Hanne Rognebakke
\brief Reads a variable in the data list and the value.
*/
static int read_bin_var(FILE *fp, char **o_varname, int *o_nvec, int **o_ivec,
			double **o_dvec, char **o_cvec, int *o_nchar, int *o_type)
{
  int         i,nchar,nvec,type;
  int        *ivec=NULL;
  double     *dvec=NULL;
  char        buffer[MAX_STR];
  char        varname[MAX_STR];
  char        cvec[MAX_STR];
  size_t      res;

  
  res = fread(&nchar,sizeof(int),1,fp);
  res = fread(&buffer,sizeof(char),nchar,fp);
  sprintf(varname,"%.*s",nchar,buffer); // Note! Copy into varname to get the right number of characters
  res = fread(&nvec,sizeof(int),1,fp);
  res = fread(&type,sizeof(int),1,fp);
  #ifdef DEBUG_PROG
  fprintf(stderr,"nchar=%d, nvec=%d, type=%d, varname=%s\n",nchar,nvec,type,varname);
  #endif
  if(type==0)
    {
      ivec = CALLOC(nvec,int);
      res = fread(ivec,sizeof(int),nvec,fp);
      #ifdef DEBUG_PROG
      for(i=0;i<MIN(nvec,10);i++)
      	fprintf(stderr,"%d ",ivec[i]);
      fprintf(stderr,"\n");
      #endif
    }
  else if(type==1)
    {
      dvec = CALLOC(nvec,double);
      res = fread(dvec,sizeof(double),nvec,fp);
      #ifdef DEBUG_PROG
      for(i=0;i<MIN(nvec,10);i++)
      	fprintf(stderr,"%f ",dvec[i]);
      fprintf(stderr,"\n");
      #endif
    }
  else if(type==2)
    {
      for(i=0;i<nvec;i++)
	{
	  res = fread(&nchar,sizeof(int),1,fp);
	  res = fread(&cvec,sizeof(char),nchar,fp);
	}
      if(nvec>1)
	fprintf(stderr,"read_bin_var: Error!: Not implemented for reading character strings in vectors !!\n");
    }
  else
    {
      write_warning("read_bin_var: Something is wrong: type not recognized\n");
    }
  *o_varname = varname;
  *o_nvec = nvec;
  *o_ivec = ivec;
  *o_dvec = dvec;
  *o_cvec = cvec;
  *o_nchar = nchar;
  *o_type = type;

  return(0);
}		/* end of read_bin_var */
    

/*!
  \author Hanne Rognebakke
  \brief Makes structs of type Input_common, Data_orig, Input_age, Data_CC and Input_lga

  Space allocated in this routine is reallocated in ::re_alloc_objects
*/
int alloc_objects_age_lga(Input_common **o_inCommon, Data_orig **o_D_orig, Input_age **o_inAge,
			  Input_lga **o_inLga, Data_CC **o_D_CC, Input_wgl **o_inHsz, Input_prior **o_inPrior)
{
  Input_common  *inCommon;
  Data_orig     *D_orig;
  Input_age     *inAge;
  Data_CC       *D_CC;
  Input_lga     *inLga;
  Input_prior   *inPrior;
  Input_wgl     *inHsz;

  D_orig = CALLOC(1,Data_orig);     // Free ok

  inCommon = CALLOC(1,Input_common);     // Free ok
  
  inAge = CALLOC(1,Input_age);     // Free ok
  inAge->cov = CALLOC(1,Input_cov);      // Free ok

  inLga = CALLOC(1,Input_lga);     // Free ok
  inLga->int_cov = CALLOC(1,Input_cov);      // Free ok
  inLga->slp_cov = CALLOC(1,Input_cov);      // Free ok

  inHsz = CALLOC(1,Input_wgl);     // Free ok
  inHsz->int_cov = CALLOC(1,Input_cov);      // Free ok
  inHsz->slp_cov = CALLOC(1,Input_cov);      // Free ok

  D_CC = CALLOC(1,Data_CC);     // Free ok

  inPrior = CALLOC(1,Input_prior);       // Free ok

  *o_inCommon = inCommon;
  *o_D_orig = D_orig;
  *o_inAge = inAge;
  *o_D_CC = D_CC;
  *o_inLga = inLga;
  *o_inHsz = inHsz;
  *o_inPrior = inPrior;

  return(0);
}		/* end of alloc_objects_age_lga */
     
    

/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in ::alloc_objects
*/
int re_alloc_objects_age_lga(Input_common **o_inCommon,Data_orig **o_D_orig, Input_age **o_inAge,
			     Input_lga **o_inLga, Data_CC **o_D_CC, Input_wgl **o_inHsz, Input_prior **o_inPrior)
{
  Input_common   *inCommon;
  Data_orig      *D_orig;
  Input_age      *inAge;
  Data_CC        *D_CC;
  Input_lga      *inLga;
  Input_wgl      *inHsz;
  Input_prior    *inPrior;

  inCommon = *o_inCommon;
  D_orig = *o_D_orig;
  inAge = *o_inAge;
  D_CC = *o_D_CC;
  inLga = *o_inLga;
  inHsz = *o_inHsz;
  inPrior = *o_inPrior;
 
  if(D_orig->coastal_cod)
    {
      FREE(D_CC->ptype1_CC1);
      FREE(D_CC->ptype1_S5);
      FREE(D_CC->ptype2_CC2);
      FREE(D_CC->ptype2_S4);
      FREE(D_CC->ptype4_CC2);
      FREE(D_CC->ptype4_S4);
      FREE(D_CC->ptype5_CC1);
      FREE(D_CC->ptype5_S5);
    }
  FREE(D_CC);
  
  FREE(inAge->a_vec);
  FREE(inAge->cov->n_lev);
  FREE(inAge->cov->random);
  FREE(inAge->cov->spatial);
  FREE(inAge->cov->interaction);
  FREE(inAge->cov->continuous);
  FREE(inAge->cov->c_cov_i);
  FREE(inAge->cov->c_cov_d);
  FREE(inAge->cov->Sigma_cell);
  FREE(inAge->cov->constr_cell);
  FREE(inAge->cov);
  FREE(inAge->num_adj_area);
  FREE(inAge->adj_area);
  if(inAge->errors)
    FREE(inAge->A2A);
  FREE(inAge);
  
  FREE(inLga->int_cov->n_lev);
  FREE(inLga->int_cov->random);
  FREE(inLga->int_cov->spatial);
  FREE(inLga->int_cov->interaction);
  FREE(inLga->int_cov->continuous);
  FREE(inLga->int_cov->c_cov_i);
  FREE(inLga->int_cov->c_cov_d);
  FREE(inLga->int_cov->Sigma_cell);
  FREE(inLga->int_cov->constr_cell);
  FREE(inLga->int_cov);
  
  FREE(inLga->slp_cov->n_lev);
  FREE(inLga->slp_cov->random);
  FREE(inLga->slp_cov->spatial);
  FREE(inLga->slp_cov->interaction);
  FREE(inLga->slp_cov->continuous);
  FREE(inLga->slp_cov->c_cov_i);
  FREE(inLga->slp_cov->Sigma_cell);
  FREE(inLga->slp_cov->constr_cell);
  FREE(inLga->slp_cov);
  
  FREE(inLga->num_adj_area);
  FREE(inLga->adj_area);

  if(inLga->g_a_model)
    FREE(inLga->g_a_par_init);
  if(inLga->fixed_model)
    {
      FREE(inLga->fixed_int);
      FREE(inLga->fixed_slp);
      FREE(inLga->fixed_tau);
      if(inLga->g_a_model)
	{
	  FREE(inLga->fixed_g_a_c);
	  FREE(inLga->fixed_g_a_theta);
	  FREE(inLga->fixed_g_a_gamma);
	}
      if(D_orig->coastal_cod)
	{
	  FREE(inLga->fixed_int_CC);
	  FREE(inLga->fixed_slp_CC);
	  FREE(inLga->fixed_tau_CC);
	  if(inLga->g_a_model)
	    {
	      FREE(inLga->fixed_g_a_c_CC);
	      FREE(inLga->fixed_g_a_theta_CC);
	      FREE(inLga->fixed_g_a_gamma_CC);
	    }
	}
    }
  
  FREE(inLga);

  if(inCommon->inc_hsz)
    {
      FREE(inHsz->int_cov->n_lev);
      FREE(inHsz->int_cov->random);
      FREE(inHsz->int_cov->spatial);
      FREE(inHsz->int_cov->interaction);
      FREE(inHsz->int_cov->continuous);
      FREE(inHsz->int_cov->c_cov_i);
      FREE(inHsz->int_cov->c_cov_d);
      FREE(inHsz->int_cov->Sigma_cell);
      FREE(inHsz->int_cov->constr_cell);
      FREE(inHsz->int_cov);
      
      FREE(inHsz->slp_cov->n_lev);
      FREE(inHsz->slp_cov->random);
      FREE(inHsz->slp_cov->spatial);
      FREE(inHsz->slp_cov->interaction);
      FREE(inHsz->slp_cov->continuous);
      FREE(inHsz->slp_cov->c_cov_i);
      FREE(inHsz->slp_cov->c_cov_d);
      FREE(inHsz->slp_cov->Sigma_cell);
      FREE(inHsz->slp_cov->constr_cell);
      FREE(inHsz->slp_cov);
      
      FREE(inHsz->num_adj_area);
      FREE(inHsz->adj_area);
      
      if(inHsz->fixed_model)
	{
	  FREE(inHsz->fixed_int);
	  FREE(inHsz->fixed_slp);
	  FREE(inHsz->fixed_tau);
	}
      FREE(inHsz);
      FREE(D_orig->haulweight);
    }
  
  FREE(D_orig->boat);
  FREE(D_orig->nFishBoat);
  FREE(D_orig->totage);
  FREE(D_orig->totlength);
  FREE(D_orig->lstart);
  FREE(D_orig->lend);
  FREE(D_orig->replength);
  FREE(D_orig->season);
  FREE(D_orig->start_noAge);
  FREE(D_orig->num_noAge);
  FREE(D_orig->int_len_lim);
  FREE(D_orig->int_len);
  FREE(D_orig->tottype);
  FREE(D_orig);

  FREE(inCommon->num_par1);
  FREE(inCommon->num_par2);
  FREE(inCommon->filename_mcmc1);
  FREE(inCommon->filename_mcmc2);
  FREE(inCommon);
  
  
  FREE(inPrior->age_eff_mean);
  FREE(inPrior->age_eff_prec);
  FREE(inPrior->age_prec_par);
  FREE(inPrior->age_ar);
  FREE(inPrior->lga_eff_mean);
  FREE(inPrior->lga_eff_prec);
  FREE(inPrior->lga_prec_par);
  FREE(inPrior->lga_ar);
  FREE(inPrior->wgl_eff_mean);
  FREE(inPrior->wgl_eff_prec);
  FREE(inPrior->wgl_prec_par);
  FREE(inPrior->wgl_ar);
  FREE(inPrior);

  return(0);
}		/* end of re_alloc_objects_age_lga */

    

/*!
  \author Hanne Rognebakke
  \brief Makes structs of type Input_common, Data_orig, Input_wgl and Data_CC

  Space allocated in this routine is reallocated in ::re_alloc_objects_wgl
*/
int alloc_objects_wgl(Input_common **o_inCommon, Data_orig **o_D_orig, Input_wgl **o_inWgl,
		      Data_CC **o_D_CC, Input_prior **o_inPrior)
{
  Input_common  *inCommon;
  Data_orig     *D_orig;
  Data_CC       *D_CC;
  Input_prior   *inPrior;
  Input_wgl     *inWgl;

  D_orig = CALLOC(1,Data_orig);     // Free ok

  inCommon = CALLOC(1,Input_common);     // Free ok
  
  inWgl = CALLOC(1,Input_wgl);     // Free ok
  inWgl->int_cov = CALLOC(1,Input_cov);      // Free ok
  inWgl->slp_cov = CALLOC(1,Input_cov);      // Free ok

  D_CC = CALLOC(1,Data_CC);     // Free ok

  inPrior = CALLOC(1,Input_prior);       // Free ok

  *o_inCommon = inCommon;
  *o_D_orig = D_orig;
  *o_inWgl = inWgl;
  *o_D_CC = D_CC;
  *o_inPrior = inPrior;

  return(0);
}		/* end of alloc_objects_wgl */
     
    

/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in ::alloc_objects_wgl
*/
int re_alloc_objects_wgl(Input_common **o_inCommon,Data_orig **o_D_orig, 
			 Input_wgl **o_inWgl, Data_CC **o_D_CC, Input_prior **o_inPrior)
{
  Input_common   *inCommon;
  Data_orig      *D_orig;
  Input_wgl      *inWgl;
  Data_CC        *D_CC;
  Input_prior    *inPrior;

  inCommon = *o_inCommon;
  D_orig = *o_D_orig;
  inWgl = *o_inWgl;
  D_CC = *o_D_CC;
  inPrior = *o_inPrior;
 
  FREE(inCommon->num_par1);
  FREE(inCommon->num_par2);
  FREE(inCommon->filename_mcmc1);
  FREE(inCommon->filename_mcmc2);
  FREE(inCommon);
  
  if(D_orig->coastal_cod)
    {
      FREE(D_CC->ptype1_CC1);
      FREE(D_CC->ptype1_S5);
      FREE(D_CC->ptype2_CC2);
      FREE(D_CC->ptype2_S4);
      FREE(D_CC->ptype4_CC2);
      FREE(D_CC->ptype4_S4);
      FREE(D_CC->ptype5_CC1);
      FREE(D_CC->ptype5_S5);
    }
  FREE(D_CC);
  
  FREE(inWgl->int_cov->n_lev);
  FREE(inWgl->int_cov->random);
  FREE(inWgl->int_cov->spatial);
  FREE(inWgl->int_cov->interaction);
  FREE(inWgl->int_cov->continuous);
  FREE(inWgl->int_cov->c_cov_i);
  FREE(inWgl->int_cov->c_cov_d);
  FREE(inWgl->int_cov->Sigma_cell);
  FREE(inWgl->int_cov->constr_cell);
  FREE(inWgl->int_cov);
  
  FREE(inWgl->slp_cov->n_lev);
  FREE(inWgl->slp_cov->random);
  FREE(inWgl->slp_cov->spatial);
  FREE(inWgl->slp_cov->interaction);
  FREE(inWgl->slp_cov->continuous);
  FREE(inWgl->slp_cov->c_cov_i);
  FREE(inWgl->slp_cov->c_cov_d);
  FREE(inWgl->slp_cov->Sigma_cell);
  FREE(inWgl->slp_cov->constr_cell);
  FREE(inWgl->slp_cov);
  
  FREE(inWgl->num_adj_area);
  FREE(inWgl->adj_area);

  if(inWgl->fixed_model)
    {
      FREE(inWgl->fixed_int);
      FREE(inWgl->fixed_slp);
      FREE(inWgl->fixed_tau);
      if(D_orig->coastal_cod)
	{
	  FREE(inWgl->fixed_int_CC);
	  FREE(inWgl->fixed_slp_CC);
	  FREE(inWgl->fixed_tau_CC);
	}
    }
  FREE(inWgl);
  
  FREE(D_orig->boat);
  FREE(D_orig->nFishBoat);
  FREE(D_orig->totlength);
  FREE(D_orig->lstart);
  FREE(D_orig->lend);
  FREE(D_orig->totweight);
  FREE(D_orig->replength);
  FREE(D_orig->haulweight);
  FREE(D_orig->start_noAge);
  FREE(D_orig->num_noAge);
  FREE(D_orig->int_len_lim);
  FREE(D_orig->int_len);
  FREE(D_orig->tottype);
  FREE(D_orig);
  
  FREE(inPrior->age_eff_mean);
  FREE(inPrior->age_eff_prec);
  FREE(inPrior->age_prec_par);
  FREE(inPrior->age_ar);
  FREE(inPrior->lga_eff_mean);
  FREE(inPrior->lga_eff_prec);
  FREE(inPrior->lga_prec_par);
  FREE(inPrior->lga_ar);
  FREE(inPrior->wgl_eff_mean);
  FREE(inPrior->wgl_eff_prec);
  FREE(inPrior->wgl_prec_par);
  FREE(inPrior->wgl_ar);
  FREE(inPrior);

  return(0);
}		/* end of re_alloc_objects_wgl */

    

/*!
  \author Hanne Rognebakke
  \brief Makes structs of type Input_common

  Space allocated in this routine is reallocated in ::re_alloc_objects
*/
int readdata_common_ascii(char *i_filename, Input_common *i_inCommon, Data_orig *i_D_orig, Data_CC *i_D_CC)
{
  int            i,j,nvar,nvec;
  char           buffer[MAX_STR],varname[MAX_STR];
  char          *inputfolder,*fname_mcmc1,*fname_mcmc2,*fname_hsz_mcmc2,*fname_hsz_it,*fname_hsz_hauleff;
  FILE          *fp;
  size_t         res;

  i_inCommon->constr = 1;
  i_inCommon->print_format = 0;  //print_format =0 (binary), 1 (ascii)
  i_inCommon->old_version = 0;
  
  /* Open file for reading in ascii format */
  if(!(fp = fopen(i_filename, "r")))
    {
      sprintf(buffer,"readdata_common_ascii: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }
  res = fscanf(fp,"%d",&nvar);

  for(i=0;i<nvar;i++)
    {
      res = fscanf(fp,"%s %d",varname,&nvec);
      #ifdef DEBUG_PROG
      fprintf(stderr,"nvec=%d, varname=%s\n",nvec,varname);
      #endif
      
      if (strcmp(varname, "seed") == 0)
	{
	  res = fscanf(fp,"%d",&i_inCommon->seed);
	}
      else if(strcmp(varname, "mcmc.par") == 0)
	{
	  if(nvec != 3)
	    {
	      sprintf(buffer,"readdata_common: wrong number of parameters in mcmc.par (%d)\n",nvec);
	      write_warning(buffer);
	      return(1);
	    }
	  res = fscanf(fp,"%d %d %d",&i_inCommon->burn_in,&i_inCommon->num_it_inner,
		       &i_inCommon->num_it_outer);
	}
      else if(strcmp(varname, "num.par.model1") == 0)
	{
	  if(nvec != 5)
	    {
	      sprintf(buffer,"readdata_common: wrong number of parameters in num.par.model1 (%d must be 5)\n",nvec);
	      write_warning(buffer);
	      return(1);
	    }
	  i_inCommon->num_par1 = CALLOC(5,int);
	  for(j=0;j<5;j++)
	    res = fscanf(fp,"%d",&i_inCommon->num_par1[j]);
	}
      else if(strcmp(varname, "num.par.model2") == 0)
	{
	  if(nvec != 2)
	    {
	      sprintf(buffer,"readdata_common: wrong number of parameters in num.par.model1 (%d must be 2)\n",nvec);
	      write_warning(buffer);
	      return(1);
	    }
	  i_inCommon->num_par2 = CALLOC(2,int);
	  for(j=0;j<2;j++)
	    res = fscanf(fp,"%d",&i_inCommon->num_par2[j]);
	}
      else if(strcmp(varname, "inputfolder") == 0)
	{
	  inputfolder = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",inputfolder);
	  i_inCommon->inputfolder = inputfolder;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"inputfolder = %s\n",i_inCommon->inputfolder);
	  #endif
	}
      else if(strcmp(varname, "filename.mcmc1") == 0)
	{
	  fname_mcmc1 = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_mcmc1);
	  i_inCommon->filename_mcmc1 = fname_mcmc1;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename mcmc1 = %s\n",i_inCommon->filename_mcmc1);
	  #endif
	}
      else if(strcmp(varname, "filename.mcmc2") == 0)
	{
	  fname_mcmc2 = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_mcmc2);
	  i_inCommon->filename_mcmc2 = fname_mcmc2;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename mcmc2 = %s\n",i_inCommon->filename_mcmc2);
	  #endif
	}
      else if(strcmp(varname, "filename.hsz.mcmc2") == 0)
	{
	  fname_hsz_mcmc2 = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_hsz_mcmc2);
	  i_inCommon->filename_hsz_mcmc2 = fname_hsz_mcmc2;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename hsz_mcmc2 = %s\n",i_inCommon->filename_hsz_mcmc2);
	  #endif
	}
      else if(strcmp(varname, "filename.hsz.it") == 0)
	{
	  fname_hsz_it = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_hsz_it);
	  i_inCommon->filename_hsz_it = fname_hsz_it;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename hsz_it = %s\n",i_inCommon->filename_hsz_it);
	  #endif
	}
      else if(strcmp(varname, "filename.hsz.hauleff") == 0)
	{
	  fname_hsz_hauleff = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_hsz_hauleff);
	  i_inCommon->filename_hsz_hauleff = fname_hsz_hauleff;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename hsz_hauleff = %s\n",i_inCommon->filename_hsz_hauleff);
	  #endif
	}
      else if(strcmp(varname, "sim.ar") == 0)
	{
	  res = fscanf(fp,"%d",&i_inCommon->sim_ar);
	}
      else if(strcmp(varname, "usedebug") == 0)
	{
	  res = fscanf(fp,"%d",&i_inCommon->use_debug);
	}
      else if(strcmp(varname, "print.boat") == 0)
	{
	  res = fscanf(fp,"%d",&i_inCommon->print_boat);
	}
      else if(strcmp(varname, "print.format") == 0)
	{
	  res = fscanf(fp,"%d",&i_inCommon->print_format);
	}
      else if(strcmp(varname, "old.version") == 0)
	{
	  res = fscanf(fp,"%d",&i_inCommon->old_version);
	}
      else if(strcmp(varname, "inc.haulsize") == 0)
	{
	  res = fscanf(fp,"%d",&i_inCommon->inc_hsz);
	}
      else if (strcmp(varname, "n.int.len") == 0)
	{
	  res = fscanf(fp,"%d",&i_D_orig->n_int_len);
	}
      else if(strcmp(varname, "int.len.lim") == 0)
	{
	  i_D_orig->int_len_lim = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_orig->int_len_lim[j]);	    
	}
      else if(strcmp(varname, "int.len.vec") == 0)
	{
	  i_D_orig->int_len = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_orig->int_len[j]);	    
	}
      else if (strcmp(varname, "coastal.cod") == 0)
	{
	  res = fscanf(fp,"%d",&i_D_orig->coastal_cod);
	}
      else if (strcmp(varname, "CCerror") == 0)
	{
	  res = fscanf(fp,"%d",&i_D_CC->class_error);
	}
      else if(strcmp(varname, "ptype1.CC") == 0)
	{
	  i_D_CC->ptype1_CC1 = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_CC->ptype1_CC1[j]);	    
	}
      else if(strcmp(varname, "ptype1.S") == 0)
	{
	  i_D_CC->ptype1_S5 = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_CC->ptype1_S5[j]);	    
	}
      else if(strcmp(varname, "ptype2.CC") == 0)
	{
	  i_D_CC->ptype2_CC2 = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_CC->ptype2_CC2[j]);	    
	}
      else if(strcmp(varname, "ptype2.S") == 0)
	{
	  i_D_CC->ptype2_S4 = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_CC->ptype2_S4[j]);	    
	}
      else if(strcmp(varname, "ptype4.CC") == 0)
	{
	  i_D_CC->ptype4_CC2 = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_CC->ptype4_CC2[j]);	    
	}
      else if(strcmp(varname, "ptype4.S") == 0)
	{
	  i_D_CC->ptype4_S4 = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_CC->ptype4_S4[j]);	    
	}
      else if(strcmp(varname, "ptype5.CC") == 0)
	{
	  i_D_CC->ptype5_CC1 = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_CC->ptype5_CC1[j]);	    
	}
      else if(strcmp(varname, "ptype5.S") == 0)
	{
	  i_D_CC->ptype5_S5 = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&i_D_CC->ptype5_S5[j]);	    
	}
      else
	{
	  sprintf(buffer,"readdata_common_ascii: Unknown variable read: %s\n",varname);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"%s\n",buffer);
	  #endif
	  write_warning(buffer);
	  return(1);
	}
    }
  fclose(fp);
  
  return(0);
}		/* end of readdata_common_ascii*/

    

/*!
  \author Hanne Rognebakke
  \brief Makes structs of type Input_common

  Space allocated in this routine is reallocated in ::re_alloc_objects
*/
int readdata_common(char *i_filename, Input_common *i_inCommon, Data_orig *i_D_orig, Data_CC *i_D_CC)
{
  int            nvar,i,nchar,nvec,type,err;
  int           *ivec;
  double        *dvec;
  char           buffer[MAX_STR];
  char          *varname,*cvec;
  FILE          *fp;
  size_t         res;

  i_inCommon->constr = 1;
  i_inCommon->print_format = 0;  //print_format =0 (binary), 1 (ascii)
  i_inCommon->old_version = 0;
  
  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"readdata_common: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }
  
  res = fread(&nvar,sizeof(int),1,fp);

  for(i=0;i<nvar;i++)
    {
      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type); 
      if(err)
	{
	  write_warning("readdata_common:Error calling read_bin_var\n");
	  return(1);
	}

      if (strcmp(varname, "seed") == 0)
	{
	  if(type != 0)
	    {
	      sprintf(buffer,"readdata_common: wrong type of parameters in seed (%d must be 0)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  i_inCommon->seed = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "mcmc.par") == 0)
	{
	  if(nvec != 3)
	    {
	      sprintf(buffer,"readdata_common: wrong number of parameters in mcmc.par (%d)\n",nvec);
	      write_warning(buffer);
	      return(1);
	    }
	  i_inCommon->burn_in = ivec[0];
	  i_inCommon->num_it_inner = ivec[1];
	  i_inCommon->num_it_outer = ivec[2];
	  FREE(ivec);
	}
      else if(strcmp(varname, "num.par.model1") == 0)
	{
	  if(nvec != 5)
	    {
	      sprintf(buffer,"readdata_common: wrong number of parameters in num.par.model1 (%d must be 5)\n",nvec);
	      write_warning(buffer);
	      return(1);
	    }
	  i_inCommon->num_par1 = ivec;
	}
      else if(strcmp(varname, "num.par.model2") == 0)
	{
	  if(nvec != 2)
	    {
	      sprintf(buffer,"readdata_common: wrong number of parameters in num.par.model1 (%d must be 2)\n",nvec);
	      write_warning(buffer);
	      return(1);
	    }
	  i_inCommon->num_par2 = ivec;
	}
      else if(strcmp(varname, "filename.mcmc1") == 0)
	{
	  i_inCommon->filename_mcmc1 = CALLOC(nchar,char);
	  strncpy(i_inCommon->filename_mcmc1,cvec,nchar);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"common filename mcmc1 = %s nchar=%d\n",i_inCommon->filename_mcmc1,nchar);
	  #endif
	}
      else if(strcmp(varname, "filename.mcmc2") == 0)
	{
	  i_inCommon->filename_mcmc2 = CALLOC(nchar,char);
	  strncpy(i_inCommon->filename_mcmc2,cvec,nchar);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"common filename mcmc2 = %s nchar=%d\n",i_inCommon->filename_mcmc2,nchar);
	  #endif
	}
      else if(strcmp(varname, "filename.hsz.mcmc2") == 0)
	{
	  i_inCommon->filename_hsz_mcmc2 = CALLOC(nchar,char);
	  strncpy(i_inCommon->filename_hsz_mcmc2,cvec,nchar);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"common filename hsz_mcmc2 = %s nchar=%d\n",i_inCommon->filename_hsz_mcmc2,nchar);
	  #endif
	}
      else if(strcmp(varname, "filename.hsz.it") == 0)
	{
	  i_inCommon->filename_hsz_it = CALLOC(nchar,char);
	  strncpy(i_inCommon->filename_hsz_it,cvec,nchar);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"common filename hsz_it = %s nchar=%d\n",i_inCommon->filename_hsz_it,nchar);
	  #endif
	}
      else if(strcmp(varname, "filename.hsz.hauleff") == 0)
	{
	  i_inCommon->filename_hsz_hauleff = CALLOC(nchar,char);
	  strncpy(i_inCommon->filename_hsz_hauleff,cvec,nchar);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"common filename hsz_hauleff = %s nchar=%d\n",i_inCommon->filename_hsz_hauleff,nchar);
	  #endif
	}
      else if(strcmp(varname, "sim.ar") == 0)
	{
	  i_inCommon->sim_ar = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "usedebug") == 0)
	{
	  i_inCommon->use_debug = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "print.boat") == 0)
	{
	  i_inCommon->print_boat = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "print.format") == 0)
	{
	  i_inCommon->print_format = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "old.version") == 0)
	{
	  i_inCommon->old_version = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "inc.haulsize") == 0)
	{
	  i_inCommon->inc_hsz = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "n.int.len") == 0)
	{
	  i_D_orig->n_int_len = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "int.len.lim") == 0)
	{
	  i_D_orig->int_len_lim = dvec;
	}
      else if(strcmp(varname, "int.len.vec") == 0)
	{
	  i_D_orig->int_len = dvec;
	}
      else if (strcmp(varname, "coastal.cod") == 0)
	{
	  i_D_orig->coastal_cod = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "ptype1.CC") == 0)
	{
	  i_D_CC->ptype1_CC1 = dvec;
	}
      else if(strcmp(varname, "ptype1.S") == 0)
	{
	  i_D_CC->ptype1_S5 = dvec;
	}
      else if(strcmp(varname, "ptype2.CC") == 0)
	{
	  i_D_CC->ptype2_CC2 = dvec;
	}
      else if(strcmp(varname, "ptype2.S") == 0)
	{
	  i_D_CC->ptype2_S4 = dvec;
	}
      else if(strcmp(varname, "ptype4.CC") == 0)
	{
	  i_D_CC->ptype4_CC2 = dvec;
	}
      else if(strcmp(varname, "ptype4.S") == 0)
	{
	  i_D_CC->ptype4_S4 = dvec;
	}
      else if(strcmp(varname, "ptype5.CC") == 0)
	{
	  i_D_CC->ptype5_CC1 = dvec;
	}
      else if(strcmp(varname, "ptype5.S") == 0)
	{
	  i_D_CC->ptype5_S5 = dvec;
	}
      else if (strcmp(varname, "CCerror") == 0)
	{
	  i_D_CC->class_error = ivec[0];
	  FREE(ivec);
	}
      else
	{
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"readdata_common: Unknown variable read: %s\n",varname);
	  #endif
	}
    }

  fclose(fp);
  
  return(0);
}		/* end of readdata_common */
     
    

    
/*!
  \author Hanne Rognebakke
  \brief Makes a struct of type containing input age data and length data

  Space allocated in this routine is reallocated in ::re_alloc_objects
*/
int readdata_age_lga(char *i_filename, Data_orig *i_D_orig, Input_age *i_inAge, Input_lga *i_inLga)
{
  int           nvar,var_count,i,j,nchar,nvec,type,ncov,ncov_i,ncov_d,err;
  int          *ivec;
  double       *dvec;
  char          buffer[MAX_STR];
  char         *varname,*cvec;
  FILE         *fp;
  size_t        res;

  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"readdata_age_lga: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

  res = fread(&nvar,sizeof(int),1,fp);

  var_count = 0;
  while(var_count<nvar)
    {
      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type);
      var_count++;
      if(err)
	{
	  write_warning("readdata_age:Error calling read_bin_var\n");
	}

      if (strcmp(varname, "a.vec") == 0)
	{
	  i_inAge->nAges = nvec;
	  i_inAge->a_vec = ivec;
	}
      else if(strcmp(varname, "age") == 0)
	{
	  i_D_orig->nFish = nvec;
	  i_D_orig->totage = ivec;
	}
      else if(strcmp(varname, "part.year") == 0)
	{
	  i_D_orig->part_year = dvec;
	}
      else if(strcmp(varname, "samplingID") == 0)
	{
	  i_D_orig->samplingID = ivec;
	}
      else if(strcmp(varname, "nFishBoat") == 0)
	{
	  i_D_orig->nHaul = nvec;
	  i_D_orig->nFishBoat = ivec;
	}
      else if(strcmp(varname, "lengthCM") == 0)
	{
	  i_D_orig->totlength = dvec;
	}
      else if(strcmp(varname, "start.noAge") == 0)
	{
	  i_D_orig->start_noAge = ivec;
	}
      else if(strcmp(varname, "num.noAge") == 0)
	{
	  i_D_orig->num_noAge = ivec;
	}
      else if(strcmp(varname, "otolithtype") == 0)
	{
	  i_D_orig->tottype = ivec;
	}
      else if (strcmp(varname, "nlev") == 0)
	{
	  i_inAge->cov->n_cov = nvec;
	  i_inAge->cov->n_lev = ivec;
	  //ivec = CALLOC(nvec,int); // Create and copy new vector for Lga object
	  //for(i=0;i<nvec;i++)
	  //ivec[i] = i_inAge->cov->n_lev[i];
	  //i_inLga->int_cov->n_cov = nvec;
	  //i_inLga->int_cov->n_lev = ivec;
	}
      else if (strcmp(varname, "n.col.cov") == 0)  // covariates
	{
	  ncov = ivec[0];
	  ncov_i = 0;
	  ncov_d = 0;
	  i_inAge->cov->c_cov_i = CALLOC(ncov,int *);
	  i_inAge->cov->c_cov_d = CALLOC(ncov,double *);
	  i_inAge->cov->ispat = -1;
	  i_inAge->cov->icell = -1;
	  i_inAge->cov->iboat = -1;
	  i_inAge->cov->ihaulsize = -1;
	  //i_inLga->int_cov->c_cov_i = CALLOC(ncov,int *);
	  //i_inLga->int_cov->c_cov_d = CALLOC(ncov,double *);
	  //i_inLga->int_cov->ispat = -1;
	  for(j=0;j<ncov;j++)
	    {
	      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type);
	      var_count++;
	      if(type==0)
		{
		  i_inAge->cov->c_cov_i[ncov_i] = ivec;
		  //ivec = CALLOC(nvec,int); // Create and copy new vector for Lga object
		  //for(i=0;i<nvec;i++)
		  //ivec[i] = i_inAge->cov->c_cov_i[ncov_i][i];
		  //i_inLga->int_cov->c_cov_i[ncov_i] = ivec;
		  ncov_i++;
		  if(strcmp(varname, "cell") ==0)
		    i_inAge->cov->icell = j;
		}
	      else if(type==1)
		{
		  i_inAge->cov->c_cov_d[ncov_d] = dvec;
		  //dvec = CALLOC(nvec,double); // Create and copy new vector for Lga object
		  //for(i=0;i<nvec;i++)
		  //dvec[i] = i_inAge->cov->c_cov_d[ncov_d][i];
		  //i_inLga->int_cov->c_cov_d[ncov_d] = dvec;
		  ncov_d++;
		  if(strcmp(varname, "haulcount") ==0) // FIND ihaulsize FROM "interaction" INSTEAD??
		    i_inAge->cov->ihaulsize = j;
		  else{
		    write_warning("readdata_age_lga:Unknown continuous covariate\n");
		    return(1);
		  }
		}
	      else
		{
		  write_warning("readdata_age_lga:Unknown type for covariate\n");
		  return(1);
		}
	    }
	  //i_inLga->int_cov->icell = i_inAge->cov->icell;
	  //i_inLga->int_cov->ihaulsize = i_inAge->cov->ihaulsize;
	  i_inAge->cov->n_cov_i = ncov_i;
	  i_inAge->cov->n_cov_d = ncov_d;
	  //i_inLga->int_cov->n_cov_i = ncov_i;
	  //i_inLga->int_cov->n_cov_d = ncov_d;
	}
      else if (strcmp(varname, "random") == 0)   //int.fix
	{
	  i_inAge->cov->random = ivec;
	  //ivec = CALLOC(nvec,int); // Create and copy new vector for Lga object
	  //for(i=0;i<nvec;i++)
	  //ivec[i] = i_inAge->cov->random[i];
	  //i_inLga->int_cov->random = ivec;
	}
      else if (strcmp(varname, "CAR") == 0)
	{
	  i_inAge->cov->spatial = ivec;
	  //ivec = CALLOC(nvec,int); // Create and copy new vector for Lga object
	  //for(i=0;i<nvec;i++)
	  //ivec[i] = i_inAge->cov->spatial[i];
	  //i_inLga->int_cov->spatial = ivec;
	}
      else if (strcmp(varname, "continuous") == 0)
	{
	  i_inAge->cov->continuous = ivec;
	  //ivec = CALLOC(nvec,int); // Create and copy new vector for Lga object
	  //for(i=0;i<nvec;i++)
	  //ivec[i] = i_inAge->cov->continuous[i];
	  //i_inLga->int_cov->continuous = ivec;
	}
      else if (strcmp(varname, "interaction") == 0)
	{
	  i_inAge->cov->interaction = ivec;
	  //ivec = CALLOC(nvec,int); // Create and copy new vector for Lga object
	  //for(i=0;i<nvec;i++)
	  //ivec[i] = i_inAge->cov->interaction[i];
	  //i_inLga->int_cov->interaction = ivec;
	}
      else if(strcmp(varname, "in.slopeModel") == 0)
	{
	  i_inAge->cov->in_slopeModel = ivec;
	}
      else if(strcmp(varname, "in.landings") == 0)
	{
	  i_inAge->cov->in_landings = ivec;
	  //ivec = CALLOC(nvec,int); // Create and copy new vector for Lga object
	  //for(i=0;i<nvec;i++)
	  //ivec[i] = i_inAge->cov->in_landings[i];
	  //i_inLga->int_cov->in_landings = ivec;
	}
      else if (strcmp(varname, "Sigma") == 0)  //int.Sigma.cell
	{
	  i_inAge->cov->Sigma_cell = dvec;
	}
      else if (strcmp(varname, "constr") == 0)   //int.constr.cell
	{
	  i_inAge->cov->constr_cell = dvec;
	}
      else if (strcmp(varname, "n.constr.cell") == 0)
	{
	  i_inAge->cov->nconstr_cell = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "num.adj.area") == 0)
	{
	  i_inAge->num_adj_area = ivec;
	  ivec = CALLOC(nvec,int); // Create and copy new vector for Lga object
	  for(i=0;i<nvec;i++)
	    ivec[i] = i_inAge->num_adj_area[i];
	  i_inLga->num_adj_area = ivec;
	}
      else if (strcmp(varname, "adj.area") == 0)
	{
	  i_inAge->adj_area = ivec;
	  ivec = CALLOC(nvec,int); // Create and copy new vector for Lga object
	  for(i=0;i<nvec;i++)
	    ivec[i] = i_inAge->adj_area[i];
	  i_inLga->adj_area = ivec;
	}
      else if (strcmp(varname, "age.errors") == 0)
	{
	  i_inAge->errors = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "A2A") == 0)
	{
	  i_inAge->A2A = dvec;
	}
      else if (strcmp(varname, "delta.age") == 0)
	{
	  i_inAge->delta_age = dvec[0];
	  FREE(dvec);
	}
      else
	{
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"readdata_age_lga: Unknown variable read: %s\n",varname);
	  #endif
	}
    }


  fclose(fp);
  
  return(0);
}		/* end of readdata_age_lga */
     
    


/*!
  \author Hanne Rognebakke
  \brief Makes a struct of type containing input lga data

  Space allocated in this routine is reallocated in ::re_readdata_lga
*/
int readdata_lga(char *i_filename, Data_orig *i_D_orig, Input_lga *i_inLga)
{
  int           nvar,var_count,i,j,nchar,nvec,type,ncov,ncov_i,ncov_d,err;
  int          *ivec;
  double       *dvec;
  char          buffer[MAX_STR];
  char         *varname,*cvec;
  FILE         *fp;
  size_t        res;

  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"readdata_lga: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

  res = fread(&nvar,sizeof(int),1,fp);
 

  // default values - if not input from R
  i_inLga->g_a_nSeason = 0;
  i_inLga->g_a_sample_c = 0;
  i_inLga->g_a_sample_theta = 0;
  i_inLga->g_a_sample_gamma = 1;
  i_inLga->fixed_model = 0;

  var_count = 0;
  while(var_count<nvar)
    {
      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type); 
      var_count++;
      if(err)
	{
	  write_warning("readdata_lga:Error calling read_bin_var\n");
	}

      if (strcmp(varname, "Sigma") == 0)
	{
	  i_inLga->int_cov->Sigma_cell = dvec;
	}
      else if (strcmp(varname, "constr") == 0) 
	{
	  i_inLga->int_cov->constr_cell = dvec;
	}
      else if (strcmp(varname, "n.constr.cell") == 0)
	{
	  i_inLga->int_cov->nconstr_cell = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "nlev") == 0)
	{
	  i_inLga->int_cov->n_cov = nvec;
	  i_inLga->int_cov->n_lev = ivec;
	}
      else if (strcmp(varname, "n.col.cov") == 0)  // covariates
	{
	  ncov = ivec[0];
	  ncov_i = 0;
	  ncov_d = 0;
	  i_inLga->int_cov->c_cov_i = CALLOC(ncov,int *);
	  i_inLga->int_cov->c_cov_d = CALLOC(ncov,double *);
	  i_inLga->int_cov->ispat = -1;
	  i_inLga->int_cov->icell = -1;
	  i_inLga->int_cov->iboat = -1;
	  i_inLga->int_cov->ihaulsize = -1;
	  for(j=0;j<ncov;j++)
	    {
	      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type);
	      var_count++;
	      if(type==0)
		{
		  i_inLga->int_cov->c_cov_i[ncov_i]= ivec;
		  ncov_i++;
		  if(strcmp(varname, "cell") ==0)
		    i_inLga->int_cov->icell = j;
		}
	      else if(type==1)
		{
		  i_inLga->int_cov->c_cov_d[ncov_d] = dvec;
		  ncov_d++;
		  if(strcmp(varname, "haulcount") ==0) // FIND ihaulsize FROM "interaction" INSTEAD??
		    i_inLga->int_cov->ihaulsize = j;
		  else{
		    write_warning("readdata_lga:Unknown continuous covariate\n");
		    return(1);
		  }
		}
	      else
		{
		  write_warning("readdata_lga:Unknown type for covariate\n");
		  return(1);
		}
	    }
	  i_inLga->int_cov->n_cov_i = ncov_i;
	  i_inLga->int_cov->n_cov_d = ncov_d;
	}
      else if (strcmp(varname, "random") == 0)   //int.fix
	{
	  i_inLga->int_cov->random = ivec;
	}
      else if (strcmp(varname, "CAR") == 0)
	{
	  i_inLga->int_cov->spatial = ivec;
	}
      else if (strcmp(varname, "continuous") == 0)
	{
	  i_inLga->int_cov->continuous = ivec;
	}
      else if (strcmp(varname, "interaction") == 0)
	{
	  i_inLga->int_cov->interaction = ivec;
	}
      else if(strcmp(varname, "in.slopeModel") == 0)
	{
	  i_inLga->int_cov->in_slopeModel = ivec;
	}
      else if(strcmp(varname, "in.landings") == 0)
	{
	  i_inLga->int_cov->in_landings = ivec;
	}
      else if (strcmp(varname, "g.a.model") == 0)
	{
	  i_inLga->g_a_model = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "g.a.nSeason") == 0)
	{
	  i_inLga->g_a_nSeason = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "g.a.sample.c") == 0)
	{
	  i_inLga->g_a_sample_c = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "g.a.sample.theta") == 0)
	{
	  i_inLga->g_a_sample_theta = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "g.a.sample.gamma") == 0)
	{
	  i_inLga->g_a_sample_gamma = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "g.a.par.init") == 0)
	{
	  i_inLga->g_a_par_init = dvec;
	}
      else if (strcmp(varname, "fixed.model") == 0)
	{
	  i_inLga->fixed_model = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "fixed.int") == 0)
	{
	  i_inLga->fixed_int = dvec;
	}
      else if (strcmp(varname, "fixed.slp") == 0)
	{
	  i_inLga->fixed_slp = dvec;
	}
      else if (strcmp(varname, "fixed.tau") == 0)
	{
	  i_inLga->fixed_tau = dvec;
	}
      else if (strcmp(varname, "fixed.g.a.c") == 0)
	{
	  i_inLga->fixed_g_a_c = dvec;
	}
      else if (strcmp(varname, "fixed.g.a.theta") == 0)
	{
	  i_inLga->fixed_g_a_theta = dvec;
	}
      else if (strcmp(varname, "fixed.g.a.gamma") == 0)
	{
	  i_inLga->fixed_g_a_gamma = dvec;
	}
      else if (strcmp(varname, "fixed.int.cc") == 0)
	{
	  i_inLga->fixed_int_CC = dvec;
	}
      else if (strcmp(varname, "fixed.slp.cc") == 0)
	{
	  i_inLga->fixed_slp_CC = dvec;
	}
      else if (strcmp(varname, "fixed.tau.cc") == 0)
	{
	  i_inLga->fixed_tau_CC = dvec;
	}
      else if (strcmp(varname, "fixed.g.a.c.cc") == 0)
	{
	  i_inLga->fixed_g_a_c_CC = dvec;
	}
      else if (strcmp(varname, "fixed.g.a.theta.cc") == 0)
	{
	  i_inLga->fixed_g_a_theta_CC = dvec;
	}
      else if (strcmp(varname, "fixed.g.a.gamma.cc") == 0)
	{
	  i_inLga->fixed_g_a_gamma_CC = dvec;
	}
      else
	{
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"readdata_lga: Unknown variable read: %s\n",varname);
	  #endif
	}
    }

  fclose(fp);

  return(0);
}		/* end of readdata_lga */
     
    
 
/*!
  \author Hanne Rognebakke
  \brief Makes a struct of type containing input wgl data

  Space allocated in this routine is reallocated in ::re_readdata_wgl
*/
int readdata_wgl(char *i_filename, Data_orig *i_D_orig, Input_wgl *i_inWgl)
{ 
  int           nvar,var_count,j,nchar,nvec,type,ncov,ncov_i,ncov_d,err;
  int          *ivec;
  double       *dvec;
  char          buffer[MAX_STR];
  char         *varname,*cvec;
  FILE         *fp;
  size_t        res;

  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"readdata_wgl: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

  res = fread(&nvar,sizeof(int),1,fp);

  var_count = 0;
  while(var_count<nvar)
    {
      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type);
      var_count++;
      if(err)
	{
	  write_warning("readdata_wgl:Error calling read_bin_var\n");
	}

      if(strcmp(varname, "weightKG") == 0)
	{
	  i_D_orig->nFish = nvec;
	  i_D_orig->totweight = dvec;
	}
      else if(strcmp(varname, "lengthCM") == 0)
	{
	  i_D_orig->totlength = dvec;
	}
      else if(strcmp(varname, "samplingID") == 0)
	{
	  i_D_orig->samplingID = ivec;
	}
      else if(strcmp(varname, "otolithtype") == 0)
	{
	  i_D_orig->tottype = ivec;
	}
      else if(strcmp(varname, "nFishBoat") == 0)
	{
	  i_D_orig->nHaul = nvec;
	  i_D_orig->nFishBoat = ivec;
	}
      else if (strcmp(varname, "nlev") == 0)
	{
	  i_inWgl->int_cov->n_cov = nvec;
	  i_inWgl->int_cov->n_lev = ivec;
	}
      else if (strcmp(varname, "n.col.cov") == 0)  // covariates
	{
	  ncov = ivec[0];
	  ncov_i = 0;
	  ncov_d = 0;
	  i_inWgl->int_cov->c_cov_i = CALLOC(ncov,int *);
	  i_inWgl->int_cov->c_cov_d = CALLOC(ncov,double *);
	  i_inWgl->int_cov->ispat = -1;
	  i_inWgl->int_cov->icell = -1;
	  i_inWgl->int_cov->iboat = -1;
	  i_inWgl->int_cov->ihaulsize = -1;
	  for(j=0;j<ncov;j++)
	    {
	      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type);
	      var_count++;
	      if(type==0)
		{
		  i_inWgl->int_cov->c_cov_i[ncov_i] = ivec;
		  if(strcmp(varname, "cell") ==0)
		    i_inWgl->int_cov->icell = j;
		  ncov_i++;
		}
	      else if(type==1)
		{
		  i_inWgl->int_cov->c_cov_d[ncov_d] = dvec;
		  if(strcmp(varname, "haulcount") ==0)
		    i_inWgl->int_cov->ihaulsize = j;
		  ncov_d++;
		}
	      else
		{
		  write_warning("readdata_wgl:Unknown type for covariate\n");
		  return(1);
		}
	    }
	  i_inWgl->int_cov->n_cov_i = ncov_i;
	  i_inWgl->int_cov->n_cov_d = ncov_d;
	}
      else if (strcmp(varname, "random") == 0)
	{
	  i_inWgl->int_cov->random = ivec;
	}
      else if (strcmp(varname, "CAR") == 0)
	{
	  i_inWgl->int_cov->spatial = ivec;
	  for(j=0;j<nvec;j++){
	    if(ivec[j]==1){
	      i_inWgl->int_cov->ispat = j;
	    }
	  }
	}
      else if (strcmp(varname, "interaction") == 0)
	{
	  i_inWgl->int_cov->interaction = ivec;
	}
      else if (strcmp(varname, "continuous") == 0)
	{
	  i_inWgl->int_cov->continuous = ivec;
	}
      else if(strcmp(varname, "in.landings") == 0)
	{
	  i_inWgl->int_cov->in_landings = ivec;
	}
      else if (strcmp(varname, "Sigma") == 0)
	{
	  i_inWgl->int_cov->Sigma_cell = dvec;
	}
      else if (strcmp(varname, "constr") == 0) 
	{
	  i_inWgl->int_cov->constr_cell = dvec;
	}
      else if (strcmp(varname, "n.constr.cell") == 0)
	{
	  i_inWgl->int_cov->nconstr_cell = ivec[0];
	  FREE(ivec);
	}
     else if (strcmp(varname, "num.adj.area") == 0)
	{
	  i_inWgl->num_adj_area = ivec;
	}
      else if (strcmp(varname, "adj.area") == 0)
	{
	  i_inWgl->adj_area = ivec;
	}
      else if (strcmp(varname, "fixed.model") == 0)
	{
	  i_inWgl->fixed_model = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "fixed.int") == 0)
	{
	  i_inWgl->fixed_int = dvec;
	}
      else if (strcmp(varname, "fixed.slp") == 0)
	{
	  i_inWgl->fixed_slp = dvec;
	}
      else if (strcmp(varname, "fixed.tau") == 0)
	{
	  i_inWgl->fixed_tau = dvec;
	}
      else if (strcmp(varname, "fixed.int.cc") == 0)
	{
	  i_inWgl->fixed_int_CC = dvec;
	}
      else if (strcmp(varname, "fixed.slp.cc") == 0)
	{
	  i_inWgl->fixed_slp_CC = dvec;
	}
      else if (strcmp(varname, "fixed.tau.cc") == 0)
	{
	  i_inWgl->fixed_tau_CC = dvec;
	}
      else
	{
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"readdata_wgl: Unknown variable read: %s\n",varname);
	  #endif
	}
    }



  fclose(fp);
  
  return(0);
}		/* end of readdata_wgl */
     
    
 
/*!
  \author Hanne Rognebakke
  \brief Makes a struct of type containing input wgl data, but contains only haulsize data

  Space allocated in this routine is reallocated in ::re_readdata_hsz
*/
int readdata_hsz(char *i_filename, Data_orig *i_D_orig, Input_wgl *i_inHsz)
{ 
  int           nvar,var_count,j,nchar,nvec,type,ncov,ncov_i,ncov_d,err;
  int          *ivec;
  double       *dvec;
  char          buffer[MAX_STR];
  char         *varname,*cvec;
  FILE         *fp;
  size_t        res;

  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"readdata_hsz: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

  res = fread(&nvar,sizeof(int),1,fp);

  var_count = 0;
  while(var_count<nvar)
    {
      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type);
      var_count++;
      if(err)
	{
	  write_warning("readdata_hsz:Error calling read_bin_var\n");
	}

      if (strcmp(varname, "hsz") == 0)
	{
	  i_D_orig->haulweight = dvec;
	}
      else if (strcmp(varname, "nlev") == 0)
	{
	  i_inHsz->int_cov->n_cov = nvec;
	  i_inHsz->int_cov->n_lev = ivec;
	}
      else if (strcmp(varname, "n.col.cov") == 0)  // covariates
	{
	  ncov = ivec[0];
	  ncov_i = 0;
	  ncov_d = 0;
	  i_inHsz->int_cov->c_cov_i = CALLOC(ncov,int *);
	  i_inHsz->int_cov->c_cov_d = CALLOC(ncov,double *);
	  i_inHsz->int_cov->ispat = -1;
	  i_inHsz->int_cov->icell = -1;
	  i_inHsz->int_cov->iboat = -1;
	  i_inHsz->int_cov->ihaulsize = -1;
	  for(j=0;j<ncov;j++)
	    {
	      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type);
	      var_count++;
	      if(type==0)
		{
		  i_inHsz->int_cov->c_cov_i[ncov_i] = ivec;
		  if(strcmp(varname, "cell") ==0)
		    i_inHsz->int_cov->icell = j;
		  ncov_i++;
		}
	      else if(type==1)
		{
		  i_inHsz->int_cov->c_cov_d[ncov_d] = dvec;
		  if(strcmp(varname, "haulcount") ==0)
		    i_inHsz->int_cov->ihaulsize = j;
		  ncov_d++;
		}
	      else
		{
		  write_warning("readdata_hsz:Unknown type for covariate\n");
		  return(err);
		}
	    }
	  i_inHsz->int_cov->n_cov_i = ncov_i;
	  i_inHsz->int_cov->n_cov_d = ncov_d;
	}
      else if (strcmp(varname, "random") == 0)
	{
	  i_inHsz->int_cov->random = ivec;
	}
      else if (strcmp(varname, "CAR") == 0)
	{
	  i_inHsz->int_cov->spatial = ivec;
	  for(j=0;j<nvec;j++){
	    if(ivec[j]==1){
	      i_inHsz->int_cov->ispat = j;
	    }
	  }
	}
      else if (strcmp(varname, "interaction") == 0)
	{
	  i_inHsz->int_cov->interaction = ivec;
	}
      else if (strcmp(varname, "continuous") == 0)
	{
	  i_inHsz->int_cov->continuous = ivec;
	}
      else if(strcmp(varname, "in.landings") == 0)
	{
	  i_inHsz->int_cov->in_landings = ivec;
	}
      else if (strcmp(varname, "Sigma") == 0)
	{
	  i_inHsz->int_cov->Sigma_cell = dvec;
	}
      else if (strcmp(varname, "constr") == 0) 
	{
	  i_inHsz->int_cov->constr_cell = dvec;
	}
      else if (strcmp(varname, "n.constr.cell") == 0)
	{
	  i_inHsz->int_cov->nconstr_cell = ivec[0];
	  FREE(ivec);
	}
     else if (strcmp(varname, "num.adj.area") == 0)
	{
	  i_inHsz->num_adj_area = ivec;
	}
      else if (strcmp(varname, "adj.area") == 0)
	{
	  i_inHsz->adj_area = ivec;
	}
      else if (strcmp(varname, "fixed.model") == 0)
	{
	  i_inHsz->fixed_model = ivec[0];
	  FREE(ivec);
	}
      else if (strcmp(varname, "fixed.int") == 0)
	{
	  i_inHsz->fixed_int = dvec;
	}
      else if (strcmp(varname, "fixed.slp") == 0)
	{
	  i_inHsz->fixed_slp = dvec;
	}
      else if (strcmp(varname, "fixed.tau") == 0)
	{
	  i_inHsz->fixed_tau = dvec;
	}
      else
	{
	  fprintf(stderr,"readdata_hsz: Unknown variable read: %s\n",varname);
	}
    }



  fclose(fp);
  
  return(0);
}		/* end of readdata_hsz*/
     
    

/*!
  \author Hanne Rognebakke
  \brief Makes a struct of type containing input prior data for age and lga model

  Space allocated in this routine is reallocated in ::re_readdata_prior
*/
int readdata_prior(char *i_filename, Input_prior *i_inPrior)
{
  int           nvar,i,nchar,nvec,type,err;
  int          *ivec;
  double       *dvec;
  char          buffer[MAX_STR];
  char         *varname,*cvec;
  FILE         *fp;
  size_t        res;


  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"readdata_prior: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

  res = fread(&nvar,sizeof(int),1,fp);

  for(i=0;i<nvar;i++)
    {
      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type); 
      if(err)
	{
	  write_warning("readdata_prior:Error calling read_bin_var\n");
	}

      if (strcmp(varname, "age.eff.mean") == 0)
	{
	  i_inPrior->age_eff_mean = dvec;
	}
      else if (strcmp(varname, "age.eff.prec") == 0)
	{
	  i_inPrior->age_eff_prec = dvec;
	}
      else if (strcmp(varname, "age.prec.par") == 0)
	{
	  i_inPrior->age_prec_par = dvec;
	}
      else if (strcmp(varname, "age.ar") == 0)
	{
	  i_inPrior->age_ar = dvec;
	}
      else if (strcmp(varname, "lga.eff.mean") == 0)
	{
	  i_inPrior->lga_eff_mean = dvec;
	}
      else if (strcmp(varname, "lga.eff.prec") == 0)
	{
	  i_inPrior->lga_eff_prec = dvec;
	}
      else if (strcmp(varname, "lga.prec.par") == 0)
	{
	  i_inPrior->lga_prec_par = dvec;
	}
      else if (strcmp(varname, "lga.ar") == 0)
	{
	  i_inPrior->lga_ar = dvec;
	}
      else if (strcmp(varname, "wgl.eff.mean") == 0)
	{
	  i_inPrior->wgl_eff_mean = dvec;
	}
      else if (strcmp(varname, "wgl.eff.prec") == 0)
	{
	  i_inPrior->wgl_eff_prec = dvec;
	}
      else if (strcmp(varname, "wgl.prec.par") == 0)
	{
	  i_inPrior->wgl_prec_par = dvec;
	}
      else if (strcmp(varname, "wgl.ar") == 0)
	{
	  i_inPrior->wgl_ar = dvec;
	}
      else
	{
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"readdata_prior: Unknown variable read: %s\n",varname);
	  #endif
	}
    }


  fclose(fp);
  
  return(0);
}		/* end of readdata_prior */
     
    


/*!
  \author Hanne Rognebakke
  \brief Makes structs of type Input_predict

  Space allocated in this routine is reallocated in ::re_readdata_predict
*/
int readdata_predict_ascii(char *i_filename, Input_predict **o_inPredict)
{
  Input_predict *inPredict;
  int            i,j,nvar,nvec;
  char           buffer[MAX_STR],varname[MAX_STR];
  char          *inputfolder,*fname_mcmc1,*fname_mcmc2,*fname_predict,*fname_hsz_mcmc2,*fname_hsz_it;
  FILE          *fp;
  size_t         res;

  inPredict = CALLOC(1,Input_predict);     // Free ok

  /* Open file for reading in ascii format */
  if(!(fp = fopen(i_filename, "r")))
    {
      sprintf(buffer,"readdata_predict_ascii: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }
  res = fscanf(fp,"%d",&nvar);

  for(i=0;i<nvar;i++)
    {
      res = fscanf(fp,"%s %d",varname,&nvec);
      //fprintf(stderr,"readdata_predict_ascii: var=%s, nvec=%d\n",varname,nvec);
      
      if (strcmp(varname, "burnin") == 0)
	{
	  res = fscanf(fp,"%d",&inPredict->burnin);
	}
      else if(strcmp(varname, "int") == 0)
	{
	  res = fscanf(fp,"%d",&j);
	}
      else if(strcmp(varname, "N.lint") == 0)
	{
	  res = fscanf(fp,"%d",&inPredict->N_l_int);
	}
      else if(strcmp(varname, "l.int") == 0)
	{
 	  inPredict->l_int = CALLOC(nvec,double);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%lf",&inPredict->l_int[j]);
	}
      else if(strcmp(varname, "nMCvec") == 0)
	{
 	  inPredict->n_MC = CALLOC(nvec,int);
	  for(j=0;j<nvec;j++)
	    res = fscanf(fp,"%d",&inPredict->n_MC[j]);
	}
      else if(strcmp(varname, "inputfolder") == 0)
	{
	  inputfolder = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",inputfolder);
	  inPredict->inputfolder = inputfolder;
	}
      else if(strcmp(varname, "filename.mcmc1") == 0)
	{
	  fname_mcmc1 = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_mcmc1);
	  inPredict->filename_mcmc1 = fname_mcmc1;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename mcmc1 = %s\n",inPredict->filename_mcmc1);
	  #endif
	}
      else if(strcmp(varname, "filename.mcmc2") == 0)
	{
	  fname_mcmc2 = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_mcmc2);
	  inPredict->filename_mcmc2 = fname_mcmc2;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename mcmc2 = %s\n",inPredict->filename_mcmc2);
	  #endif
	}
      else if(strcmp(varname, "filename.predict") == 0)
	{
	  fname_predict = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_predict);
	  inPredict->filename_predict = fname_predict;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename predict = %s\n",inPredict->filename_predict);
	  #endif
	}
      else if(strcmp(varname, "filename.hsz.mcmc2") == 0)
	{
	  fname_hsz_mcmc2 = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_hsz_mcmc2);
	  inPredict->filename_hsz_mcmc2 = fname_hsz_mcmc2;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename hsz_mcmc2 = %s\n",inPredict->filename_hsz_mcmc2);
	  #endif
	}
      else if(strcmp(varname, "filename.hsz.it") == 0)
	{
	  fname_hsz_it = malloc(MAX_STR*sizeof(char));
	  res = fscanf(fp,"%s",fname_hsz_it);
	  inPredict->filename_hsz_it = fname_hsz_it;
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"filename hsz_it = %s\n",inPredict->filename_hsz_it);
	  #endif
	}
      else if(strcmp(varname, "inc.haulsize") == 0)
	{
	  res = fscanf(fp,"%d",&inPredict->inc_hsz);
	}
     else
	{
	  sprintf(buffer,"readdata_predict_ascii: Unknown variable read: %s\n",varname);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"%s\n",buffer);
	  #endif
	  write_warning(buffer);
	  return(1);
	}
    }
  fclose(fp);

  *o_inPredict = inPredict;

  return(0);
}		/* end of readdata_predict_ascii */
    


/*!
  \author Hanne Rognebakke
  \brief Makes structs of type Input_predict

  Space allocated in this routine is reallocated in ::re_readdata_predict
*/
int readdata_predict(char *i_filename, Input_predict **o_inPredict)
{
  Input_predict *inPredict;
  int            nvar,i,nchar,nvec,type,err;
  int           *ivec;
  double        *dvec;
  char           buffer[MAX_STR];
  char          *varname,*cvec;
  FILE          *fp;
  size_t         res;

  inPredict = CALLOC(1,Input_predict);     // Free ok

  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"readdata_predict: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

  
  res = fread(&nvar,sizeof(int),1,fp);

  for(i=0;i<nvar;i++)
    {
      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type); 
      if(err)
	{
	  write_warning("readdata_predict:Error calling read_bin_var\n");
	}

      if (strcmp(varname, "burnin") == 0)
	{
	  if(type != 0)
	    {
	      sprintf(buffer,"readdata_predict: wrong type of parameters in burnin (%d must be 0)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  inPredict->burnin = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "N.lint") == 0)
	{
	  if(type != 0)
	    {
	      sprintf(buffer,"readdata_predict: wrong type of parameters in N.lint (%d must be 0)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  inPredict->N_l_int= ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "l.int") == 0)
	{
	  if(type != 1)
	    {
	      sprintf(buffer,"readdata_predict: wrong type of parameters in l.int (%d must be 1)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  inPredict->l_int = dvec;
	}
      else if(strcmp(varname, "nMCvec") == 0)
	{
	  if(type != 0)
	    {
	      sprintf(buffer,"readdata_predict: wrong type of parameters in nMCvec (%d must be 0)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  inPredict->n_MC = ivec;
	}
      else if(strcmp(varname, "filename.mcmc1") == 0)
	{
	  inPredict->filename_mcmc1 = CALLOC(nchar,char);
	  strncpy(inPredict->filename_mcmc1,cvec,nchar);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"common filename mcmc1 = %s nchar=%d\n",inPredict->filename_mcmc1,nchar);
	  #endif
	}
      else if(strcmp(varname, "filename.mcmc2") == 0)
	{
	  inPredict->filename_mcmc2 = CALLOC(nchar,char);
	  strncpy(inPredict->filename_mcmc2,cvec,nchar);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"common filename mcmc2 = %s nchar=%d\n",inPredict->filename_mcmc2,nchar);
	  #endif
	}
      else if(strcmp(varname, "filename.predict") == 0)
	{
	  inPredict->filename_predict = CALLOC(nchar,char);
	  strncpy(inPredict->filename_predict,cvec,nchar);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"common filename predict = %s nchar=%d\n",inPredict->filename_predict,nchar);
	  #endif
	}
      else if(strcmp(varname, "filename.hsz.mcmc2") == 0)
	{
	  inPredict->filename_hsz_mcmc2 = CALLOC(nchar,char);
	  strncpy(inPredict->filename_hsz_mcmc2,cvec,nchar);
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"common filename hsz_mcmc2 = %s nchar=%d\n",inPredict->filename_hsz_mcmc2,nchar);
	  #endif
	}
      else if(strcmp(varname, "inc.haulsize") == 0)
	{
	  inPredict->inc_hsz = ivec[0];
	  FREE(ivec);
	}
      else
	{
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"readdata_predict: Unknown variable read: %s\n",varname);
	  #endif
	}
    }
 
  *o_inPredict = inPredict;

  fclose(fp);
  
  return(0);
}		/* end of readdata_predict */

    

/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in ::readdata_predict
*/
int re_readdata_predict(Input_predict **o_inPredict)
{
  Input_predict   *inPredict;

  inPredict = *o_inPredict;

  FREE(inPredict->num_par1);
  FREE(inPredict->num_par2);
  FREE(inPredict->l_int);
  FREE(inPredict->n_MC);
  FREE(inPredict->filename_mcmc1);
  FREE(inPredict->filename_mcmc2);
  FREE(inPredict->filename_predict);
  FREE(inPredict);
  
  return(0);
}		/* end of re_readdata_predict */


/*!
  \author Hanne Rognebakke
  \brief Makes structs of type Input_catch

  Space allocated in this routine is reallocated in ::re_readdata_catch
*/
int readdata_catch(char *i_filename, Input_totcatch **o_inCatch)
{
  Input_totcatch  *inCatch;
  int              nvar,var_count,j,nfac,nchar,nvec,type,tmp,err;
  int             *ivec;
  double          *dvec;
  char             buffer[MAX_STR];
  char            *varname,*cvec;
  FILE            *fp;
  size_t           res;

  inCatch = CALLOC(1,Input_totcatch);     // Free ok

  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"readdata_catch: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

 
  res = fread(&nvar,sizeof(int),1,fp);
  
  var_count = 0;
  while(var_count<nvar)
    {
      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type); 
      var_count++;
      if(err)
	{
	  write_warning("readdata_catch:Error calling read_bin_var\n");
	}

      if(strcmp(varname, "nCell") == 0)
	{
	  if(type != 0)
	    {
	      sprintf(buffer,"readdata_catch: wrong type of parameters in nCell (%d must be 0)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  inCatch->nCell= ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "nfactors") == 0)
	{
	  if(type != 0)
	    {
	      sprintf(buffer,"readdata_catch: wrong type of parameters in nfactors (%d must be 0)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  inCatch->nFactors = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "fac.age.int") == 0)
	{
	  inCatch->fac_age_int = ivec;
	}
      else if(strcmp(varname, "fac.hsz.int") == 0)
	{
	  inCatch->fac_hsz_int = ivec;
	}
      else if(strcmp(varname, "fac.lga.int") == 0)
	{
	  inCatch->fac_lga_int = ivec;
	}
      else if(strcmp(varname, "fac.lga.slp") == 0)
	{
	  inCatch->fac_lga_slp = ivec;
	}
      else if(strcmp(varname, "fac.wgl.int") == 0)
	{
	  inCatch->fac_wgl_int = ivec;
	}
      else if(strcmp(varname, "fac.wgl.slp") == 0)
	{
	  inCatch->fac_wgl_slp = ivec;
	}
      else if (strcmp(varname, "n.col.cov") == 0)  // covariates
	{
	  nfac = ivec[0];
	  inCatch->factors = CALLOC(nfac,int *);
	  for(j=0;j<nfac;j++)
	    {
	      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type);
	      var_count++;
	      if(type==0)
		{
		  inCatch->factors[j] = ivec;
		}
	      else
		{
		  write_warning("readdata_catch:Unknown type for factors\n");
		  return(err);
		}
	    }
	}
      else if(strcmp(varname, "catch") == 0)
	{
	  inCatch->catch = dvec;
	}
      else if(strcmp(varname, "season") == 0)
	{
	  inCatch->season = ivec;
	  write_warning("readdata_catch:Season should be replaced with input 'midseason'\n");
	}
      else if(strcmp(varname, "midseason") == 0)
	{
	  inCatch->midseason = dvec;
          inCatch->season = CALLOC(nvec,int);
	  for(j=0;j<nvec;j++)
	    {
	      tmp = (int) ceil(inCatch->midseason[j]*12);
              inCatch->season[j] = tmp;
	    }
	}
      else
	{
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"readdata_catch: Unknown variable read: %s\n",varname);
	  #endif
	}
    }
 
  *o_inCatch = inCatch;

  fclose(fp);
  
  return(0);
}		/* end of readdata_catch */

    
/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in ::readdata_catch
*/
int re_readdata_catch(Input_totcatch **o_inCatch)
{
  Input_totcatch   *inCatch;

  inCatch = *o_inCatch;

  FREE(inCatch->season);
  FREE(inCatch->midseason);
  FREE(inCatch->fac_age_int);
  FREE(inCatch->fac_lga_int);
  FREE(inCatch->fac_lga_slp);
  FREE(inCatch->fac_wgl_int);
  FREE(inCatch->fac_wgl_slp);
  FREE(inCatch->factors);
  FREE(inCatch->catch);
  FREE(inCatch);
  
  return(0);
}		/* end of re_readdata_catch */


/*!
  \author Hanne Rognebakke
  \brief Makes structs of type Input_cell

  Space allocated in this routine is reallocated in ::re_readdata_cell
*/
int readdata_cell(char *i_filename, Input_cell **o_inCell)
{
  Input_cell  *inCell;
  int          nvar,i,nchar,nvec,type,err;
  int         *ivec;
  double      *dvec;
  char         buffer[MAX_STR];
  char        *varname,*cvec;
  FILE        *fp;
  size_t       res;

  inCell = CALLOC(1,Input_cell);     // Free ok

  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"readdata_cell: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

 
  res = fread(&nvar,sizeof(int),1,fp);

   for(i=0;i<nvar;i++)
    {
      err = read_bin_var(fp, &varname, &nvec, &ivec, &dvec, &cvec, &nchar, &type); 
      if(err)
	{
	  write_warning("readdata_cell:Error calling read_bin_var\n");
	}
 
      if (strcmp(varname, "num.cell.o") == 0)
	{
	  if(type != 0)
	    {
	      sprintf(buffer,"readdata_cell: wrong type of parameters in num.cell.o (%d must be 0)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  inCell->num_cell_o = ivec;
	}
      else if (strcmp(varname, "num.cell.u") == 0)
	{
	  if(type != 0)
	    {
	      sprintf(buffer,"readdata_cell: wrong type of parameters in num.cell.u (%d must be 0)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  inCell->num_cell_u = ivec;
	}
      else if(strcmp(varname, "age.int.nC") == 0)
	{
	  if(type != 0)
	    {
	      sprintf(buffer,"readdata_cell: wrong type of parameters in age.int.nC (%d must be 0)\n",type);
	      write_warning(buffer);
	      return(1);
	    }
	  inCell->age_int_nC = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "age.int.E") == 0)
	{
	  inCell->age_int_E = dvec;
	}
      else if(strcmp(varname, "age.int.C") == 0)
	{
	  inCell->age_int_C = dvec;
	}
      else if(strcmp(varname, "hsz.int.nC") == 0)
	{
	  inCell->hsz_int_nC = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "hsz.int.E") == 0)
	{
	  inCell->hsz_int_E = dvec;
	}
      else if(strcmp(varname, "hsz.int.C") == 0)
	{
	  inCell->hsz_int_C = dvec;
	}
      else if(strcmp(varname, "lga.int.nC") == 0)
	{
	  inCell->lga_int_nC = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "lga.int.C") == 0)
	{
	  inCell->lga_int_C = dvec;
	}
      else if(strcmp(varname, "lga.int.E") == 0)
	{
	  inCell->lga_int_E = dvec;
	}
      else if(strcmp(varname, "lga.slp.nC") == 0)
	{
	  inCell->lga_slp_nC = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "lga.slp.C") == 0)
	{
	  inCell->lga_slp_C = dvec;
	}
      else if(strcmp(varname, "lga.slp.E") == 0)
	{
	  inCell->lga_slp_E = dvec;
	}
      else if(strcmp(varname, "wgl.int.nC") == 0)
	{
	  inCell->wgl_int_nC = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "wgl.int.C") == 0)
	{
	  inCell->wgl_int_C = dvec;
	}
      else if(strcmp(varname, "wgl.int.E") == 0)
	{
	  inCell->wgl_int_E = dvec;
	}
      else if(strcmp(varname, "wgl.slp.nC") == 0)
	{
	  inCell->wgl_slp_nC = ivec[0];
	  FREE(ivec);
	}
      else if(strcmp(varname, "wgl.slp.C") == 0)
	{
	  inCell->wgl_slp_C = dvec;
	}
      else if(strcmp(varname, "wgl.slp.E") == 0)
	{
	  inCell->wgl_slp_E = dvec;
	}
      else
	{
	  #ifdef DEBUG_PROG
	  fprintf(stderr,"readdata_cell: Unknown variable read: %s\n",varname);
	  #endif
	}
    }
 
  *o_inCell = inCell;

  fclose(fp);
  
  return(0);
}		/* end of readdata_cell */

    
/*!
  \author Hanne Rognebakke
  \brief Reallocate memory allocated in ::readdata_cell
*/
int re_readdata_cell(Input_cell **o_inCell)
{
  Input_cell   *inCell;

  inCell = *o_inCell;

  FREE(inCell->num_cell_o);
  FREE(inCell->num_cell_u);
  FREE(inCell->age_int_E);
  FREE(inCell->age_int_C);
  FREE(inCell->lga_int_E);
  FREE(inCell->lga_int_C);
  FREE(inCell->lga_slp_E);
  FREE(inCell->lga_slp_C);
  FREE(inCell->wgl_int_E);
  FREE(inCell->wgl_int_C);
  FREE(inCell->wgl_slp_E);
  FREE(inCell->wgl_slp_C);
  FREE(inCell);
  
  return(0);
}		/* end of re_readdata_cell */




int add_object_info_age_lga(Input_common *i_inCommon, Data_orig *i_D_orig, Input_age *i_inAge,
			    Input_lga *i_inLga, Input_wgl *i_inHsz)
{
  int i,nHaul,nFish,n_cov,max_boat;

  int old_version = 0;
  if(i_inCommon->old_version)
    old_version = 1;
  
  i_inCommon->constr = 1;
  //  i_inCommon->print_format = 1;  //print_format =0 (binary), 1 (ascii)

  nHaul = i_D_orig->nHaul;
  nFish = i_D_orig->nFish;

  if(i_inLga->g_a_nSeason==0)
    {
      i_inLga->g_a_nSeason = 12;  //Use 12 seasons for continuous age
      if(old_version){
	i_inLga->g_a_nSeason = 4;  //Use 12 seasons for continuous age
	fprintf(stderr,"add_object_info: NOTE! USING nSeason=4 !!!\n");
      }
    }

  if(i_D_orig->coastal_cod)
    i_inLga->g_a_ncat = (int) i_inAge->nAges/2 *i_inLga->g_a_nSeason;
  else 
    i_inLga->g_a_ncat = i_inAge->nAges*i_inLga->g_a_nSeason;
	  
  i_D_orig->boat = CALLOC(nHaul,int);
  i_D_orig->season = CALLOC(nFish,int);
  
  i_D_orig->replength = CALLOC(nFish,int);
  #ifdef DEBUG_PROG
  fprintf(stderr,"add_object_info: Takes log of totlength\n");
  #endif
  int tmp;
  for(i=0;i<nFish;i++)
    {
      i_D_orig->replength[i] = 1;
      i_D_orig->totlength[i] = log(i_D_orig->totlength[i]);
      //NOTE: CHANGE TO USING ONLY PART_YEAR WHEN CHANGING CONTINUOUS AGE!!
      tmp = (int) ceil(i_D_orig->part_year[i]*i_inLga->g_a_nSeason);
      if((tmp<1)|(tmp>i_inLga->g_a_nSeason))
      	fprintf(stderr,"WARNING: Season=%d\n",tmp);
      i_D_orig->season[i] = tmp;
    }

  if(i_inCommon->inc_hsz)
    {
      #ifdef DEBUG_PROG
      fprintf(stderr,"add_object_info: Takes log of haulweight\n");
      #endif
      for(i=0;i<nHaul;i++)
	{
	  i_D_orig->haulweight[i] = log(i_D_orig->haulweight[i]);
	}
      i_inHsz->num_adj_area = i_inAge->num_adj_area;
      i_inHsz->adj_area = i_inAge->adj_area;      
      i_inHsz->slp_cov->n_cov = 1;
      n_cov = 1;
      i_inHsz->slp_cov->n_lev = CALLOC(n_cov,int);
      i_inHsz->slp_cov->n_lev[0] = 1;
      i_inHsz->slp_cov->random = CALLOC(n_cov,int);
      i_inHsz->slp_cov->random[0] = 0;
      i_inHsz->slp_cov->spatial = CALLOC(n_cov,int);
      i_inHsz->slp_cov->spatial[0] = 0;
      i_inHsz->slp_cov->continuous = CALLOC(n_cov,int);
      i_inHsz->slp_cov->continuous[0] = 0;
      i_inHsz->slp_cov->interaction = CALLOC(n_cov,int);
      i_inHsz->slp_cov->interaction[0] = 0;
      i_inHsz->slp_cov->in_landings = CALLOC(n_cov,int);
      i_inHsz->slp_cov->in_landings[0] = 1;
      i_inHsz->slp_cov->c_cov_i = CALLOC(1,int *);
      i_inHsz->slp_cov->c_cov_i[0] = CALLOC(n_cov*i_D_orig->nHaul,int);
      for(i=0;i<(n_cov*i_D_orig->nHaul);i++)
	{
	  i_inHsz->slp_cov->c_cov_i[0][i] = 1;
	}
    }

  i_inLga->slp_cov->n_cov = 1;
  n_cov = 1;
  i_inLga->slp_cov->n_lev = CALLOC(n_cov,int);
  i_inLga->slp_cov->n_lev[0] = 1;
  i_inLga->slp_cov->random = CALLOC(n_cov,int);
  i_inLga->slp_cov->random[0] = 0;
  i_inLga->slp_cov->spatial = CALLOC(n_cov,int);
  i_inLga->slp_cov->spatial[0] = 0;
  i_inLga->slp_cov->continuous = CALLOC(n_cov,int);
  i_inLga->slp_cov->continuous[0] = 0;
  i_inLga->slp_cov->interaction = CALLOC(n_cov,int);
  i_inLga->slp_cov->interaction[0] = 0;
  i_inLga->slp_cov->in_landings = CALLOC(n_cov,int);
  i_inLga->slp_cov->in_landings[0] = 1;
  i_inLga->slp_cov->c_cov_i = CALLOC(1,int *);
  i_inLga->slp_cov->c_cov_i[0] = CALLOC(n_cov*i_D_orig->nHaul,int);
  for(i=0;i<(n_cov*i_D_orig->nHaul);i++)
    {
      i_inLga->slp_cov->c_cov_i[0][i] = 1;
    }

  return(0);
}		/* end of add_object_info */



int add_object_info_wgl(Input_common *i_inCommon, Data_orig *i_D_orig, Input_wgl *i_inWgl)
{
  int i,nHaul,nFish,n_cov;
  
  i_inCommon->constr = 1;
  //i_inCommon->print_format = 1; //print_format =0 (binary), 1 (ascii)

  nHaul = i_D_orig->nHaul;
  nFish = i_D_orig->nFish;

  i_D_orig->nBoat = 1;
  i_D_orig->boat = CALLOC(nHaul,int);
  i_D_orig->haulweight = CALLOC(nHaul,double);
  i_D_orig->replength = CALLOC(nFish,int);
  //  fprintf(stderr,"add_object_info_wgl: CHANGE boat, haulweight\n");
  for(i=0;i<nHaul;i++)
    {
      i_D_orig->boat[i] = 0;
      i_D_orig->haulweight[i] = 1.0;
    }
  for(i=0;i<nFish;i++)
    {
      i_D_orig->replength[i] = 1;
      i_D_orig->totweight[i] = log(i_D_orig->totweight[i]);
      i_D_orig->totlength[i] = log(i_D_orig->totlength[i]);
    }

  i_inWgl->fixed_model = 0;
  
  i_inWgl->slp_cov->n_cov = 1;
  n_cov = 1;
  i_inWgl->slp_cov->n_lev = CALLOC(n_cov,int);
  i_inWgl->slp_cov->n_lev[0] = 1;
  i_inWgl->slp_cov->random = CALLOC(n_cov,int);
  i_inWgl->slp_cov->random[0] = 0;
  i_inWgl->slp_cov->spatial = CALLOC(n_cov,int);
  i_inWgl->slp_cov->spatial[0] = 0;
  i_inWgl->slp_cov->continuous = CALLOC(n_cov,int);
  i_inWgl->slp_cov->continuous[0] = 0;
  i_inWgl->slp_cov->interaction = CALLOC(n_cov,int);
  i_inWgl->slp_cov->interaction[0] = 0;
  i_inWgl->slp_cov->in_landings = CALLOC(n_cov,int);
  i_inWgl->slp_cov->in_landings[0] = 1;
  i_inWgl->slp_cov->c_cov_i = CALLOC(1,int *);
  i_inWgl->slp_cov->c_cov_i[0] = CALLOC(n_cov*i_D_orig->nHaul,int);
  for(i=0;i<(n_cov*i_D_orig->nHaul);i++)
    {
      i_inWgl->slp_cov->c_cov_i[0][i] = 1;
    }

	  
  return(0);
}		/* end of add_object_info_wgl */


