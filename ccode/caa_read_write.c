#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "caa.h"
#include "caa_read_write.h"
#include "caa_routines.h"

#ifdef LOG_FILE
extern FILE     *g_caa_log; 
#endif
extern FILE     *g_caa_mcmc1;
extern FILE     *g_caa_mcmc2;


/*!
  \author Geir Storvik
  \brief Writes a warning message to screen
*/
void write_warning(char *i_text)
{
  #ifdef LOG_FILE
  fprintf(g_caa_log,"%s\n",i_text);
  #endif
  printf("%s\n",i_text);

  return;
}		/* end of write_warning */

void write_output(char *filename, char *i_text)
{
  FILE *fp;
  fp = fopen(filename,"a");
  fprintf(fp,"%s\n",i_text);
  fclose(fp);

  return;
}		/* end of write_output */

int printWarning(char *i_text)
{
  #ifdef LOG_FILE
  fprintf(g_caa_log,"%s\n",i_text);
  fprintf(g_caa_log,"\n");
  #endif
  printf("Warning:  ");
  printf("%s\n",i_text);
  printf("\n");

  return(1);
}		/* end of printWarning */

/*!
  \brief Print error message and exit from program.
  Copyright:    Norwegian Computing Center, SAMBA, 2008
*/
void printError(char *i_text, char *i_filename)
{
  printf("%s\n",i_text);
  printf(": %s ", i_filename);
  printf(" !!!\n");
  exit(0);
}



/*!
  \brief Write parameters from an iteration (typically last) in model1 to file

  Parameters can be used as starting values in a later run.
  \author Hanne Rognebakke
*/
int write_par_model1(char *i_filename, Data_age *i_D_age, Age_struct *i_age,
		     Data_lin *i_D_lga, LW_struct *i_length, 
		     Data_lin *i_D_lga_CC, LW_struct *i_length_CC,  
		     Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, int i_coastal_cod, Data_CC *i_D_CC)
{
  char   buffer[MAX_STR];
  int    err,i;
  FILE  *fp;

  /* Open file for printing in binary format */
  if(!(fp = fopen(i_filename, "wb")))
    {
      sprintf(buffer,"write_par_model1: Couldn't open file for writing: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }
  
  /* Age model */
  err = write_par_bin(fp,i_D_age->glm,i_age->par);
  if(err)
    {
      write_warning("write_par_model1:Error calling write_par_bin\n");
      return(err);
    }
  /* Write complete alphas for age model */  
  int h;
  for(h=0;h<i_D_age->glm->nHaul;h++)
    fwrite(i_age->alpha[h],sizeof(double),i_D_age->glm->ncat,fp);
  
  /* Lga model */
  err = write_par_bin(fp,i_D_lga->glm,i_length->par);
  if(err)
    {
      write_warning("write_par_model1:Error calling write_par_bin\n");
      return(err);
    }
  
  /* g-function */
  if(i_D_g_a->g_a_model==1)
    fwrite(i_D_g_a->g_a_par,sizeof(double),i_D_g_a->g_a_npar,fp);

  /* Coastal cod */  
  if(i_coastal_cod)
    {
      err = write_par_bin(fp,i_D_lga_CC->glm,i_length_CC->par);
      if(err)
	{
	  write_warning("write_par_model1:Error calling write_par_bin\n");
	  return(err);
	}
      /* g-function */
      if(i_D_g_a_CC->g_a_model==1)
	fwrite(i_D_g_a_CC->g_a_par,sizeof(double),i_D_g_a_CC->g_a_npar,fp);
      /* Classification error */
      fwrite(&i_D_CC->k1,sizeof(int),1,fp);
      fwrite(&i_D_CC->k2,sizeof(int),1,fp);
    }
  
  fclose(fp);
 

  if(0)//use ascii file when testing
    {
      fp = fopen(i_filename, "w");
      err = write_par_ascii(fp,i_D_age->glm,i_age->par);
      /* Write complete alphas for age model */  
      int a,h;
      for(a=0;a<i_D_age->glm->ncat;a++)
	{
	  for(h=0;h<i_D_age->glm->nHaul;h++)
	    fprintf(fp,"%f ",i_age->alpha[h][a]);
	  fprintf(fp,"\n");
	}
      err = write_par_ascii(fp,i_D_lga->glm,i_length->par);
      
      /* g-function */
      if(i_D_g_a->g_a_model==1)
	{
	  for(i=0;i<i_D_g_a->g_a_npar;i++)
	    fprintf(fp,"%f ",i_D_g_a->g_a_par[i]);
	  fprintf(fp,"\n");
	}
      
      /* Coastal cod */
      if(i_coastal_cod)
	{
	  err = write_par_ascii(fp,i_D_lga_CC->glm,i_length_CC->par);
	  /* g-function */
	  if(i_D_g_a_CC->g_a_model==1)
	    {
	      for(i=0;i<i_D_g_a_CC->g_a_npar;i++)
		fprintf(fp,"%f ",i_D_g_a_CC->g_a_par[i]);
	      fprintf(fp,"\n");
	    }
	  /* Classification error */
	  fprintf(fp,"%f %f\n",i_D_CC->k1,i_D_CC->k2);
	}
      
      fclose(fp);
    }

  
  return(0);
}               /* end of write_par_model1 */



/*!
  \brief Write parameters in binary format
  \author Hanne Rognebakke
*/
int write_par_bin(FILE *fp,Data_glm *i_glm,Eff_str *i_par)
{
  int       a,i,j;
  Data_cov *xcov;

  /* Write number of categories and factors */
  fwrite(&i_glm->nxcov,sizeof(int),1,fp);
  fwrite(&i_glm->ncat,sizeof(int),1,fp);
  for(i=0;i<i_glm->nxcov;i++)
    {
      fwrite(&i_glm->xcov[i]->n_cov,sizeof(int),1,fp);
      fwrite(i_glm->xcov[i]->n_fac,sizeof(int),i_glm->xcov[i]->n_cov,fp);
    }	

  /* Write linear effects */
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(a=0;a<i_glm->ncat;a++)
	{
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      fwrite(i_par->eff[a][i][j],sizeof(double),xcov->n_fac[j],fp);
	    }
	}
    }

  /* Write precisions for random effects*/  
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov--; //haul precision in age model stored in tau_obs
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          if(!xcov->fix[j])
	    {
	      fwrite(&i_par->tau[i][j],sizeof(double),1,fp);
	    }
	}
    }
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov++;

  /* Write observation precision */
  fwrite(&i_par->tau_obs,sizeof(double),1,fp);

  return(0);
}               /* end of write_par_bin */



/*!
  \brief Write parameters in ascii format
  \author Hanne Rognebakke
*/
int write_par_ascii(FILE *fp,Data_glm *i_glm,Eff_str *i_par)
{
  int       a,i,j,k;
  Data_cov *xcov;

  /* Write number of categories and factors */
  fprintf(fp,"%d %d\n",i_glm->nxcov,i_glm->ncat);
  for(i=0;i<i_glm->nxcov;i++)
    {
      fprintf(fp,"%d ",i_glm->xcov[i]->n_cov);
      for(j=0;j<i_glm->xcov[i]->n_cov;j++)
	fprintf(fp,"%d ",i_glm->xcov[i]->n_fac[j]);
      fprintf(fp,"\n");
    }

  /* Write linear effects */
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(a=0;a<i_glm->ncat;a++)
	{
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      for(k=0;k<xcov->n_fac[j];k++)
		{
		  fprintf(fp,"%f ",i_par->eff[a][i][j][k]);
		}
	      fprintf(fp,"\n");
	    }
	}
    }

  /* Write precisions for random effects*/  
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov--; //haul precision in age model stored in tau_obs
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          if(!xcov->fix[j])
	    {
	      fprintf(fp,"%f\n",i_par->tau[i][j]);
	    }
	}
    }
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov++;

  /* Write observation precision */
  fprintf(fp,"%f\n",i_par->tau_obs);

  return(0);
}               /* end of write_par_ascii */



/*!
  \brief Write parameters from an iteration (typically last) in model2 to file

  Parameters can be used as starting values in a later run.
  \author Hanne Rognebakke
*/
int write_par_model2(char *i_filename, Data_lin *i_D_wgl, LW_struct *i_weight, 
		     Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, int i_coastal_cod)
{
  char   buffer[MAX_STR];
  int    err;
  FILE  *fp;

  /* Open file for printing in binary format */
  if(!(fp = fopen(i_filename, "wb")))
    {
      sprintf(buffer,"write_par_model2: Couldn't open file for writing: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

  /* Wgl model */
  err = write_par_bin(fp,i_D_wgl->glm,i_weight->par);
  if(err)
    {
      write_warning("write_par_model2:Error calling write_par_bin\n");
      return(err);
    }

  if(i_coastal_cod)
    {
      err = write_par_bin(fp,i_D_wgl_CC->glm,i_weight_CC->par);
      if(err)
	{
	  write_warning("write_par_model2:Error calling write_par_bin\n");
	  return(err);
	}
    }

  fclose(fp);

  if(0) // when testing
    {
      fp = fopen(i_filename, "w");
      err = write_par_ascii(fp,i_D_wgl->glm,i_weight->par);
      if(i_coastal_cod)
	err = write_par_ascii(fp,i_D_wgl_CC->glm,i_weight_CC->par);
      
      fclose(fp);
    }

  
  return(0);
}               /* end of write_par_model2 */


/*!
  \brief Read parameters from binary file and use as starting values
  \author Hanne Rognebakke
*/
int read_par_model1(char *i_filename, Data_age *i_D_age, Age_struct *i_age,
		    Data_lin *i_D_lga, LW_struct *i_length, 
		    Data_lin *i_D_lga_CC, LW_struct *i_length_CC, 
		    Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, int i_coastal_cod, Data_CC *i_D_CC)
{
  int    i,ret,nFac;
  char   buffer[MAX_STR];
  int    err;
  double dtemp;
  FILE  *fp;

  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"read_par_model1: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }

  /* Age model */
  err = read_par_bin(fp,i_D_age->glm,i_age->par,&nFac);
  if(err)
    {
      write_warning("read_par_model1:Error calling read_par_bin\n");
      return(err);
    }
  i_D_age->glm->nHaulOld = nFac;

  err = read_complete_alpha_bin(fp,i_D_age->glm,i_age);
  if(err)
    {
      write_warning("read_par_model1:Error calling read_complete_alpha_bin\n");
      return(err);
    }

  /* Lga model */
  err = read_par_bin(fp,i_D_lga->glm,i_length->par,&nFac);
  if(err)
    {
      write_warning("read_par_model1:Error calling read_par_bin\n");
      return(err);
    }
    
  /* g-function */
  if(i_D_g_a->g_a_model==1)
    {
      ret = fread(i_D_g_a->g_a_par,sizeof(double),i_D_g_a->g_a_npar,fp);
    }

  /* Coastal cod */  
  if(i_coastal_cod)
    {
      err = read_par_bin(fp,i_D_lga_CC->glm,i_length_CC->par,&nFac);
      if(err)
	{
	  write_warning("read_par_model1:Error calling read_par_bin\n");
	  return(err);
	}
      /* g-function */
      if(i_D_g_a_CC->g_a_model==1)
	{
	  ret = fread(i_D_g_a_CC->g_a_par,sizeof(double),i_D_g_a_CC->g_a_npar,fp);
	}
      /* Classification error */
      ret = fread(&i_D_CC->k1,sizeof(int),1,fp);
      ret = fread(&i_D_CC->k2,sizeof(int),1,fp);
    }

  fclose(fp);

  if(0) // If testing
    {
      fp = fopen(i_filename, "r");
      err = read_par_ascii(fp,i_D_age->glm,i_age->par,i_age->alpha);
    
      err = read_par_ascii_lga(fp,i_D_lga->glm,i_length->par);

      /* g-function */
      if(i_D_g_a->g_a_model==1)
	{
	  for(i=0;i<i_D_g_a->g_a_npar;i++)
	    {
	      ret = fscanf(fp,"%lf",&dtemp);
	      i_D_g_a->g_a_par[i] = dtemp;
	    }
	}
      
      /* Coastal cod */
      if(i_coastal_cod)
	{
	  err = read_par_ascii_lga(fp,i_D_lga_CC->glm,i_length_CC->par);
	  if(i_D_g_a_CC->g_a_model==1)
	    {
	      for(i=0;i<i_D_g_a_CC->g_a_npar;i++)
		{
		  ret = fscanf(fp,"%lf",&dtemp);
		  i_D_g_a_CC->g_a_par[i] = dtemp;
		}
	    }
	  ret = fscanf(fp,"%lf",&dtemp);
	  i_D_CC->k1 = dtemp;
	  ret = fscanf(fp,"%lf",&dtemp);
	  i_D_CC->k2 = dtemp;
	}
      
      fclose(fp);
    }

  return(0);
}               /* end of read_par_model1 */


/*!
  \brief Read parameters in binary format and put into struct Eff_str
  \author Hanne Rognebakke
*/
int read_par_bin(FILE *fp,Data_glm *i_glm,Eff_str *i_par,int *i_nFac)
{
  int       a,i,j,ret;
  int       nxcov,ncat;
  int      *n_cov,**n_fac;
  Data_cov *xcov;

  /* Read number of categories and factors */
  ret = fread(&nxcov,sizeof(int),1,fp);
  ret = fread(&ncat,sizeof(int),1,fp);
  n_cov = CALLOC(nxcov,int);
  n_fac = CALLOC(nxcov,int *);
  for(i=0;i<nxcov;i++)
    {
      ret = fread(&n_cov[i],sizeof(int),1,fp);
      n_fac[i] = CALLOC(n_cov[i],int);
      ret = fread(n_fac[i],sizeof(int),n_cov[i],fp);
    }
  *i_nFac = n_fac[0][n_cov[0]-1];

  /* Read linear effects */
  for(i=0;i<nxcov;i++)
    {
      for(a=0;a<ncat;a++)
	{
	  for(j=0;j<n_cov[i];j++)
	    {
	      ret = fread(i_par->eff[a][i][j],sizeof(double),n_fac[i][j],fp);
	    }
	}
    }

  /* Read precisions for random effects*/  
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov--; //haul precision in age model stored in tau_obs
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          if(!xcov->fix[j])
	    {
	      ret = fread(&i_par->tau[i][j],sizeof(double),1,fp);
	    }
	}
    }
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov++;

  /* Read observation precision */
  ret = fread(&i_par->tau_obs,sizeof(double),1,fp);

  FREE(n_cov);
  for(i=0;i<nxcov;i++)
    FREE(n_fac[i]);
  FREE(n_fac);

  return(0);
}               /* end of read_par_bin */

/*!
  \brief Read parameters for complete alphas in binary format 
  \author Hanne Rognebakke
*/
int read_complete_alpha_bin(FILE *fp,Data_glm *i_glm,Age_struct *i_age)
{
  int a,h,nHaul,ret;
  
  nHaul = i_glm->nHaulOld;
  
  for(h=0;h<nHaul;h++)
    ret = fread(i_age->alpha[h],sizeof(double),i_glm->ncat,fp);
  for(h=nHaul;h<i_glm->nHaul;h++)
    for(a=0;a<i_glm->ncat;a++)
      i_age->alpha[h][a] = calc_eff_no_haul(i_glm->xcov[0],i_age->par->eff[a][0],0); //haul effect = 0
  
  return(0);
}               /* end of read_par_bin */

/*!
  \brief Read parameters in ascii format and put into struct Eff_str
  \author Hanne Rognebakke
*/
int read_par_ascii_lga(FILE *fp,Data_glm *i_glm,Eff_str *i_par)
{
  int       a,i,j,k,ret;
  int       nxcov,ncat;
  double    eff;
  int      *n_cov,**n_fac;
  Data_cov *xcov;

  /* Read number of categories and factors */
  ret = fscanf(fp,"%d %d",&nxcov,&ncat);
  n_cov = CALLOC(nxcov,int);
  n_fac = CALLOC(nxcov,int *);
  for(i=0;i<nxcov;i++)
    {
      ret = fscanf(fp,"%d",&n_cov[i]);
      n_fac[i] = CALLOC(n_cov[i],int);
      for(j=0;j<n_cov[i];j++)
	ret = fscanf(fp,"%d",&n_fac[i][j]);
    }

  /* Read linear effects */
  for(i=0;i<nxcov;i++)
    {
      for(a=0;a<ncat;a++)
	{
	  for(j=0;j<n_cov[i];j++)
	    {
	      for(k=0;k<n_fac[i][j];k++)
		{
		  ret = fscanf(fp,"%lf",&eff);
		  i_par->eff[a][i][j][k]=eff;
		}
 	    }
	}
    }  

  /* Read precisions for random effects*/  
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov--; //haul precision in age model stored in tau_obs
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          if(!xcov->fix[j])
	    {
	      ret = fscanf(fp,"%lf",&eff);
	      i_par->tau[i][j] = eff;
	    }
	}
    }
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov++;

  /* Read observation precision */
  ret = fscanf(fp,"%lf",&i_par->tau_obs);

  FREE(n_cov);
  for(i=0;i<nxcov;i++)
    FREE(n_fac[i]);
  FREE(n_fac);

  return(0);
}               /* end of read_par_ascii_lga */

/*!
  \brief Read parameters in ascii format and put into struct Eff_str
  \author Hanne Rognebakke
*/
 int read_par_ascii(FILE *fp,Data_glm *i_glm,Eff_str *i_par,double **i_alpha)
{
  int       a,i,j,k,ret;
  int       nxcov,ncat;
  int      *n_cov,**n_fac;
  Data_cov *xcov;

  /* Read number of categories and factors */
  ret = fscanf(fp,"%d %d",&nxcov,&ncat);
  n_cov = CALLOC(nxcov,int);
  n_fac = CALLOC(nxcov,int *);
  for(i=0;i<nxcov;i++)
    {
      ret = fscanf(fp,"%d",&n_cov[i]);
      n_fac[i] = CALLOC(n_cov[i],int);
      for(j=0;j<n_cov[i];j++)
	{
	  ret = fscanf(fp,"%d",&n_fac[i][j]);
	}
    }

  /* Read linear effects */
  for(i=0;i<nxcov;i++)
    {
      for(a=0;a<ncat;a++)
	{
	  for(j=0;j<n_cov[i];j++)
	    {
	      for(k=0;k<n_fac[i][j];k++)
		{
		  ret = fscanf(fp,"%lf",&i_par->eff[a][i][j][k]);
		}
 	    }
	}
    }

  /* Read precisions for random effects*/  
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov--; //haul precision in age model stored in tau_obs
  for(i=0;i<i_glm->nxcov;i++)
    {
      xcov = i_glm->xcov[i];
      for(j=0;j<xcov->n_cov;j++)
	{
          if(!xcov->fix[j])
	    {
	      ret = fscanf(fp,"%lf",&i_par->tau[i][j]);
	    }
	}
    }
  if(i_glm->ncat>1)
    i_glm->xcov[0]->n_cov++;

  /* Read observation precision */
  ret = fscanf(fp,"%lf",&i_par->tau_obs);
  #ifdef LOG_FILE
  fprintf(g_caa_log,"read_par_ascii: tau_obs=%f\n",i_par->tau_obs);
  #endif
 
  int h,nHaul;
  double dtemp;
  nHaul = n_fac[0][n_cov[0]-1];
  i_glm->nHaulOld = nHaul;
  /* Read complete alphas for age model */  
  if(i_glm->ncat>1)
    {
      for(a=0;a<i_glm->ncat;a++)
	{
	  for(h=0;h<nHaul;h++)
	    {
	      ret = fscanf(fp,"%lf",&dtemp);
	      i_alpha[h][a] = dtemp;
	    }
	  for(h=nHaul;h<i_glm->nHaul;h++)
	    {
	      i_alpha[h][a] = calc_eff_no_haul(i_glm->xcov[0],i_par->eff[a][0],0); //haul effect = 0
	    }
	}
    }
  
  FREE(n_cov);
  for(i=0;i<nxcov;i++)
    FREE(n_fac[i]);
  FREE(n_fac);

  return(0);
}               /* end of read_par_ascii */

/*!
  \brief Read parameters for model2 from binary file and use as starting values
  \author Hanne Rognebakke
*/
int read_par_model2(char *i_filename, Data_lin *i_D_wgl, LW_struct *i_weight, 
		    Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, int i_coastal_cod)
{
  char   buffer[MAX_STR];
  int    err,nFac;
  FILE  *fp;

  /* Open file for reading in binary format */
  if(!(fp = fopen(i_filename, "rb")))
    {
      sprintf(buffer,"read_par_model2: Couldn't open file for reading: %s\n",i_filename);
      write_warning(buffer);
      return(1);
    }
  
  /* Wgl model */
  err = read_par_bin(fp,i_D_wgl->glm,i_weight->par,&nFac);
  if(err)
    {
      write_warning("read_par_model2:Error calling read_par_bin\n");
      return(err);
    }
  if(i_coastal_cod)
    {
      err = read_par_bin(fp,i_D_wgl_CC->glm,i_weight_CC->par,&nFac);
      if(err)
	{
	  write_warning("read_par_model2:Error calling read_par_bin\n");
	  return(err);
	}
    }
  
  fclose(fp);

  if(0)
    {
      fp = fopen(i_filename, "r");
      err = read_par_ascii(fp,i_D_wgl->glm,i_weight->par,NULL);
      if(i_coastal_cod)
	err = read_par_ascii(fp,i_D_wgl_CC->glm,i_weight_CC->par,NULL);
      
      fclose(fp);
    }

  return(0);
}               /* end of read_par_model2 */


/*!
  \brief Write model parameters from estimation (age+lga model) to file, binary or ascii. 
  The output file is to be used in prediction.
  \author Hanne Rognebakke
*/
int write_mcmc1(Data_age *i_D_age, Age_struct *i_age, 
		Data_lin *i_D_lga, LW_struct *i_length, Data_lin *i_D_lga_CC,  
		LW_struct *i_length_CC, Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, 
		int i_nMCMC, int *i_num_par1, int i_coastal_cod, int i_print_boat,
		double i_delta_age, int i_print_format)
{
  int    i;

  if(i_print_format==0) // Binary format
    {
      /* Write common parameters */    
      fwrite(&i_nMCMC,sizeof(int),1,g_caa_mcmc1);
      fwrite(i_num_par1,sizeof(int),5,g_caa_mcmc1);
      fwrite(&i_print_boat,sizeof(int),1,g_caa_mcmc1);

      /* Write age model */
      write_glm_object(g_caa_mcmc1, i_D_age->glm, i_print_format);
      fwrite(i_D_age->a_vec,sizeof(int),i_D_age->glm->ncat,g_caa_mcmc1);
      fwrite(&i_delta_age,sizeof(double),1,g_caa_mcmc1);

      /* Write lga model */
      write_glm_object(g_caa_mcmc1, i_D_lga->glm, i_print_format);

      /* Write ga model */
      fwrite(&i_D_g_a->g_a_model,sizeof(int),1,g_caa_mcmc1);

      /* Write coastal cod parameters */
      fwrite(&i_coastal_cod,sizeof(int),1,g_caa_mcmc1);
      if(i_coastal_cod)
	{
	  write_glm_object(g_caa_mcmc1, i_D_lga_CC->glm, i_print_format);
	  /* Write ga model */
	  fwrite(&i_D_g_a_CC->g_a_model,sizeof(int),1,g_caa_mcmc1);
	}
    }
  else //Ascii format
    {
      /* Write common parameters */    
      fprintf(g_caa_mcmc1,"%d\n",i_nMCMC);
      fprintf(g_caa_mcmc1,"%d %d %d %d %d\n",i_num_par1[0],i_num_par1[1],i_num_par1[2],i_num_par1[3],i_num_par1[4]);
      fprintf(g_caa_mcmc1,"%d\n",i_print_boat);

      /* Write age model */
      write_glm_object(g_caa_mcmc1, i_D_age->glm, i_print_format);
      for(i=0;i<i_D_age->glm->ncat;i++)
	fprintf(g_caa_mcmc1,"%d ",i_D_age->a_vec[i]);
      fprintf(g_caa_mcmc1,"\n");
      fprintf(g_caa_mcmc1,"%f\n",i_delta_age);

      /* Write lga model */
      write_glm_object(g_caa_mcmc1, i_D_lga->glm, i_print_format);

      /* Write ga model */
      fprintf(g_caa_mcmc1,"%d\n",i_D_g_a->g_a_model);
      
      /* Write coastal cod parameters */
      fprintf(g_caa_mcmc1,"%d\n",i_coastal_cod);
      if(i_coastal_cod)
	{
	  write_glm_object(g_caa_mcmc1, i_D_lga_CC->glm, i_print_format);
	  /* Write ga model */
	  fprintf(g_caa_mcmc1,"%d\n",i_D_g_a_CC->g_a_model);
	}
    }

  return(0);
}               /* end of write_mcmc1 */


int write_glm_object(FILE *fp, Data_glm *i_glm, int i_print_format)
{
  int    i,j;
  
  if(i_print_format==0) // Binary format
    {
      fwrite(&i_glm->ncat,sizeof(int),1,fp);
      fwrite(&i_glm->nHaul,sizeof(int),1,fp);
      fwrite(&i_glm->nxcov,sizeof(int),1,fp);
      for(i=0;i<i_glm->nxcov;i++)
	{
	  fwrite(&i_glm->xcov[i]->n_cov,sizeof(int),1,fp);
	  fwrite(i_glm->xcov[i]->n_fac,sizeof(int),i_glm->xcov[i]->n_cov,fp);
	  fwrite(i_glm->xcov[i]->fix,sizeof(int),i_glm->xcov[i]->n_cov,fp);
	  fwrite(i_glm->xcov[i]->interaction,sizeof(int),i_glm->xcov[i]->n_cov,fp);
	  fwrite(i_glm->xcov[i]->in_landings,sizeof(int),i_glm->xcov[i]->n_cov,fp);
	}
      fwrite(&i_glm->xcov[0]->ispat,sizeof(int),1,fp);
      fwrite(&i_glm->xcov[0]->ihaul,sizeof(int),1,fp);
      fwrite(&i_glm->xcov[0]->icell,sizeof(int),1,fp);
      fwrite(&i_glm->xcov[0]->iboat,sizeof(int),1,fp);
      fwrite(&i_glm->xcov[0]->ihaulsize,sizeof(int),1,fp);
    }
  else //Ascii format
    {
      fprintf(fp,"%d\n",i_glm->ncat);
      fprintf(fp,"%d\n",i_glm->nHaul);
      fprintf(fp,"%d\n",i_glm->nxcov);
      for(i=0;i<i_glm->nxcov;i++)
	{
	  fprintf(fp,"%d\n",i_glm->xcov[i]->n_cov);
	  for(j=0;j<i_glm->xcov[i]->n_cov;j++)
	    {
	      fprintf(fp,"%d ",i_glm->xcov[i]->n_fac[j]);
	    }
	  fprintf(fp,"\n");
	  for(j=0;j<i_glm->xcov[i]->n_cov;j++)
	    {
	      fprintf(fp,"%d ",i_glm->xcov[i]->fix[j]);
	    }
	  fprintf(fp,"\n");
	  for(j=0;j<i_glm->xcov[i]->n_cov;j++)
	    {
	      fprintf(fp,"%d ",i_glm->xcov[i]->interaction[j]);
	    }
	  fprintf(fp,"\n");
	  for(j=0;j<i_glm->xcov[i]->n_cov;j++)
	    {
	      fprintf(fp,"%d ",i_glm->xcov[i]->in_landings[j]);
	    }
	  fprintf(fp,"\n");
	}
      fprintf(fp,"%d\n",i_glm->xcov[0]->ispat);
      fprintf(fp,"%d\n",i_glm->xcov[0]->ihaul);
      fprintf(fp,"%d\n",i_glm->xcov[0]->icell);
      fprintf(fp,"%d\n",i_glm->xcov[0]->iboat);
      fprintf(fp,"%d\n",i_glm->xcov[0]->ihaulsize);
    }
  #ifdef LOG_FILE
  fprintf(g_caa_log,"write_glm_object: ncat=%d, nHaul=%d, nxcov=%d\n",i_glm->ncat,i_glm->nHaul,i_glm->nxcov);
  for(i=0;i<i_glm->nxcov;i++)
    {
      fprintf(g_caa_log,"n_cov=%d\n",i_glm->xcov[i]->n_cov);
      for(j=0;j<i_glm->xcov[i]->n_cov;j++)
	fprintf(g_caa_log,"fac[%d]=%d,fix[%d]=%d,interaction[%d]=%d,in.landings[%d]=%d\n",
		j,i_glm->xcov[i]->n_fac[j],j,i_glm->xcov[i]->fix[j],j,i_glm->xcov[i]->interaction[j],j,i_glm->xcov[i]->in_landings[j]);
    }
  #endif
  
  return(0);
}         /* end of write_glm_object */


int read_glm_object(FILE *fp, Data_glm **o_glm)
{
  int        i,j,ret;
  Data_glm  *glm;
  Data_cov  *xcov;
  
  glm = CALLOC(1,Data_glm);

  ret = fread(&glm->ncat,sizeof(int),1,fp);
  ret = fread(&glm->nHaul,sizeof(int),1,fp);
  ret = fread(&glm->nxcov,sizeof(int),1,fp);
  #ifdef DEBUG_PREDICT
  fprintf(stderr,"read_glm_object: ncat=%d, nHaul=%d, nxcov=%d\n",glm->ncat,glm->nHaul,glm->nxcov);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"read_glm_object: ncat=%d, nHaul=%d, nxcov=%d\n",glm->ncat,glm->nHaul,glm->nxcov);
  #endif
  /* Covariates */
  glm->xcov = CALLOC(glm->nxcov,Data_cov *); 
  for(i=0;i<glm->nxcov;i++)
    {
      xcov = CALLOC(1,Data_cov); 
      ret = fread(&xcov->n_cov,sizeof(int),1,fp);
      xcov->n_fac = CALLOC(xcov->n_cov,int); 
      xcov->fix = CALLOC(xcov->n_cov,int);
      xcov->interaction = CALLOC(xcov->n_cov,int);
      xcov->in_landings = CALLOC(xcov->n_cov,int);
      ret = fread(xcov->n_fac,sizeof(int),xcov->n_cov,fp);
      ret = fread(xcov->fix,sizeof(int),xcov->n_cov,fp);
      ret = fread(xcov->interaction,sizeof(int),xcov->n_cov,fp);
      ret = fread(xcov->in_landings,sizeof(int),xcov->n_cov,fp);
      glm->xcov[i] = xcov;
      #ifdef DEBUG_PREDICT
      fprintf(stderr,"n_cov=%d\n",xcov->n_cov);
      for(j=0;j<xcov->n_cov;j++)
	fprintf(stderr,"fac[%d]=%d,fix[%d]=%d\n",j,xcov->n_fac[j],j,xcov->fix[j]);
      #endif
    }
  ret = fread(&glm->xcov[0]->ispat,sizeof(int),1,fp);
  ret = fread(&glm->xcov[0]->ihaul,sizeof(int),1,fp);
  ret = fread(&glm->xcov[0]->icell,sizeof(int),1,fp);
  ret = fread(&glm->xcov[0]->iboat,sizeof(int),1,fp);
  ret = fread(&glm->xcov[0]->ihaulsize,sizeof(int),1,fp);
  #ifdef DEBUG_PREDICT
  fprintf(stderr,"ispat=%d,ihaul=%d,icell=%d,iboat=%d,ihaulsize=%d\n",glm->xcov[0]->ispat,glm->xcov[0]->ihaul,glm->xcov[0]->icell,glm->xcov[0]->iboat,glm->xcov[0]->ihaulsize);
  #endif
  #ifdef LOG_FILE
  fprintf(g_caa_log,"ispat=%d,ihaul=%d,icell=%d,iboat=%d,ihaulsize=%d\n",glm->xcov[0]->ispat,glm->xcov[0]->ihaul,glm->xcov[0]->icell,glm->xcov[0]->iboat,glm->xcov[0]->ihaulsize);
  #endif
  if(glm->xcov[0]->ihaulsize>0)
    glm->inc_hsz = 1;
  else 
    glm->inc_hsz = 0;

  *o_glm = glm;
  
  return(0);
}         /* end of read_glm_object */




/*!
  \brief Read binary file from estimation (age+lga) and put data in structures. 
  Memory allocated in this routine is reallocated in ::re_makedata_predict.
  \author Hanne Rognebakke
*/
int read_mcmc1(Input_predict *i_inPredict, Data_age **o_D_age, Data_lin **o_D_lga, Data_g_a **o_D_g_a,
	       Data_lin **o_D_lga_CC, Data_g_a **o_D_g_a_CC)
{
  int        i,ret,g_a_model,g_a_nSeason,g_a_ncat;
  double    *par_init;
  Data_age  *D_age;
  Data_lin  *D_lga;
  Data_g_a  *D_g_a;
  Data_lin  *D_lga_CC=NULL;
  Data_g_a  *D_g_a_CC=NULL;
 
  /* Allocate space and initialize */
  i_inPredict->num_par1 = CALLOC(5,int);
  
  D_age = CALLOC(1,Data_age);      
  D_lga = CALLOC(1,Data_lin);      
  
  /* Read input data from binary file */
  
  /* Read common parameters */
  ret = fread(&i_inPredict->nMCMC,sizeof(int),1,g_caa_mcmc1);
  ret = fread(i_inPredict->num_par1,sizeof(int),5,g_caa_mcmc1);
  ret = fread(&i_inPredict->read_boat,sizeof(int),1,g_caa_mcmc1);
  if(i_inPredict->num_par1[3]>0)
    i_inPredict->coastal_cod = 1;
  else
    i_inPredict->coastal_cod = 0;
  #ifdef DEBUG_PREDICT
  fprintf(stderr,"\nnMCMC=%d\n",i_inPredict->nMCMC);
  fprintf(stderr,"num_par1=");
  for(i=0;i<5;i++)
    fprintf(stderr," %d",i_inPredict->num_par1[i]);
  fprintf(stderr,"\ncoastal cod=%d\n",i_inPredict->coastal_cod);
  fprintf(stderr,"read_boat=%d\n",i_inPredict->read_boat);
  #endif
	  
  /* Read age parameters */
  read_glm_object(g_caa_mcmc1, &D_age->glm);
  D_age->glm->numpar = i_inPredict->num_par1[0];
  D_age->a_vec = CALLOC(D_age->glm->ncat,int);
  ret = fread(D_age->a_vec,sizeof(int),D_age->glm->ncat,g_caa_mcmc1); 
  ret = fread(&D_age->delta_age,sizeof(double),1,g_caa_mcmc1); // delta_age
  
  /* Read lga parameters */
  read_glm_object(g_caa_mcmc1, &D_lga->glm);
  D_lga->glm->numpar = i_inPredict->num_par1[1];
  D_lga->glm->xcov[1]->ispat = -1;
  D_lga->glm->xcov[1]->ihaul = -1;
  D_lga->glm->xcov[1]->icell = -1;
  D_lga->glm->xcov[1]->iboat = -1;
  D_lga->glm->xcov[1]->ihaulsize = -1;

  /* Read ga parameters */
  ret = fread(&g_a_model,sizeof(int),1,g_caa_mcmc1);
  #ifdef DEBUG_PREDICT
  fprintf(stderr,"ga_model=%d\n",g_a_model);
  #endif
  if(g_a_model==1)
    {
      par_init = CALLOC(3,double); // Only for not getting seg.fault. Not used in prediction.
      for(i=0;i<3;i++)
	par_init[i] = 0.0;
    }
  g_a_nSeason = 12;  //Use 12 seasons for continuous age - Not needed in prediction
  g_a_ncat = D_age->glm->ncat*g_a_nSeason;

   /* Read coastal cod parameters */
  ret = fread(&i_inPredict->coastal_cod,sizeof(int),1,g_caa_mcmc1);
  #ifdef DEBUG_PREDICT
  fprintf(stderr,"coastal cod=%d\n",i_inPredict->coastal_cod);
  #endif
  ret = makedata_g_a(D_age->glm->ncat,g_a_ncat,g_a_nSeason,D_age->a_vec,i_inPredict->coastal_cod,
		     g_a_model, par_init, 0, 0, 0, &D_g_a);
  if(i_inPredict->coastal_cod)
    {
      D_lga_CC = CALLOC(1,Data_lin);      
      D_g_a_CC = CALLOC(1,Data_g_a);          
      read_glm_object(g_caa_mcmc1, &D_lga_CC->glm);
      D_lga_CC->glm->numpar = i_inPredict->num_par1[3];
      D_lga_CC->glm->xcov[1]->ispat = -1;
      D_lga_CC->glm->xcov[1]->ihaul = -1;
      D_lga_CC->glm->xcov[1]->icell = -1;
      D_lga_CC->glm->xcov[1]->iboat = -1;
      D_lga_CC->glm->xcov[1]->ihaulsize = -1;
      ret = fread(&g_a_model,sizeof(int),1,g_caa_mcmc1);
      ret = makedata_g_a(D_age->glm->ncat,g_a_ncat,g_a_nSeason,D_age->a_vec,i_inPredict->coastal_cod,
			 g_a_model, par_init, 0, 0, 0, &D_g_a_CC);
    }

  if(g_a_model==1)
    FREE(par_init);
  *o_D_age = D_age;
  *o_D_lga = D_lga;
  *o_D_g_a = D_g_a;
  *o_D_lga_CC = D_lga_CC;
  *o_D_g_a_CC = D_g_a_CC;

  return(0);
}               /* end of read_mcmc1 */



/*!
  \brief Read binary file from estimation (age+lga) and put data in structures. 
  Memory allocated in this routine is reallocated in ::re_makedata_predict.
  \author Hanne Rognebakke
*/
int read_mcmc2(Input_predict *i_inPredict, Data_lin **o_D_wgl, Data_lin **o_D_wgl_CC)
{
  int        ret,itmp;
  char       buffer[MAX_STR];
  Data_lin  *D_wgl;
  Data_lin  *D_wgl_CC=NULL;

  /* Allocate space and initialize */
  i_inPredict->num_par2 = CALLOC(2,int);
  D_wgl = CALLOC(1,Data_lin);      
  
  /* Read common parameters */
  ret = fread(&itmp,sizeof(int),1,g_caa_mcmc2);
  if(itmp != i_inPredict->nMCMC)
    {
      sprintf(buffer,"read_mcmc2: Number of MCMC samples doesn't match in the two files\n");
      write_warning(buffer);
      return(1);
    }
  ret = fread(i_inPredict->num_par2,sizeof(int),2,g_caa_mcmc2);
  ret = fread(&itmp,sizeof(int),1,g_caa_mcmc2);
  if(itmp!=i_inPredict->read_boat)
    {
      sprintf(buffer,"read_mcmc2: read_boat variable doesn't match in the two files\n");
      write_warning(buffer);
      return(1);
    }
  #ifdef DEBUG_PREDICT
  fprintf(stderr,"npar_wgl=%d %d\n",i_inPredict->num_par2[0],i_inPredict->num_par2[1]);
  #endif
  
  /* Read wgl parameters */
  read_glm_object(g_caa_mcmc2, &D_wgl->glm);
  D_wgl->glm->numpar = i_inPredict->num_par2[0];
  D_wgl->glm->xcov[1]->ispat = -1;
  D_wgl->glm->xcov[1]->ihaul = -1;
  D_wgl->glm->xcov[1]->icell = -1;
  D_wgl->glm->xcov[1]->iboat = -1;
  D_wgl->glm->xcov[1]->ihaulsize = -1;

   /* Read coastal cod parameters */
  ret = fread(&itmp,sizeof(int),1,g_caa_mcmc2);
  if(itmp != i_inPredict->coastal_cod)
    {
      sprintf(buffer,"read_mcmc2: Coastal cod parameter not equal in the two files\n");
      write_warning(buffer);
      return(1);
    }
  if(i_inPredict->coastal_cod)
    {
      D_wgl_CC = CALLOC(1,Data_lin);      
      read_glm_object(g_caa_mcmc2, &D_wgl_CC->glm);
      D_wgl_CC->glm->numpar = i_inPredict->num_par2[1];
      D_wgl_CC->glm->xcov[1]->ispat = -1;
      D_wgl_CC->glm->xcov[1]->ihaul = -1;
      D_wgl_CC->glm->xcov[1]->icell = -1;
      D_wgl_CC->glm->xcov[1]->iboat = -1;
      D_wgl_CC->glm->xcov[1]->ihaulsize = -1;
    }

  *o_D_wgl = D_wgl;
  *o_D_wgl_CC = D_wgl_CC;

  return(0);
}               /* end of read_mcmc2 */


int read_hsz(Input_predict *i_inPredict, Data_lin **o_D_wgl)
{
  int       ret,itmp;
  int       ivec[2];
  char      buffer[MAX_STR];
  Data_lin *D_wgl;

  D_wgl = CALLOC(1,Data_lin);      

  /* Read common parameters */
  ret = fread(&itmp,sizeof(int),1,g_caa_mcmc_hsz);
  i_inPredict->nMCMC_hsz = itmp;
  //fprintf(stderr,"nMCMC_hsz=%d\n",i_inPredict->nMCMC_hsz);
  ret = fread(ivec,sizeof(int),2,g_caa_mcmc_hsz);
  //fprintf(stderr,"numpar=%d,%d\n",ivec[0],ivec[1]);
  ret = fread(&itmp,sizeof(int),1,g_caa_mcmc_hsz); // print_boat parameter, included since printed from wgl model
  
  /* Read wgl parameters */
  read_glm_object(g_caa_mcmc_hsz, &D_wgl->glm);
  D_wgl->glm->numpar = ivec[0];
  D_wgl->glm->xcov[1]->ispat = -1;
  D_wgl->glm->xcov[1]->ihaul = -1;
  D_wgl->glm->xcov[1]->icell = -1;
  D_wgl->glm->xcov[1]->iboat = -1;
  D_wgl->glm->xcov[1]->ihaulsize = -1;

  D_wgl->glm->inc_hsz = 1;

  /* Read coastal cod parameters (=0, but included since printed from wgl model)*/
  ret = fread(&itmp,sizeof(int),1,g_caa_mcmc_hsz);

  *o_D_wgl = D_wgl;
  
  
  return(0);
}               /* end of read_hsz */



/*!
  \brief Write mcmc parameters from estimation (wgl model) to binary file. 
  The output file is to be used in prediction.
  \author Hanne Rognebakke
*/
int write_mcmc2(Data_lin *i_D_wgl, LW_struct *i_weight, 
		Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC,
		int i_nMCMC, int *i_num_par2, int i_coastal_cod, int i_print_boat, int i_print_format)
{
  if(i_print_format==0) // Binary format
    {
      /* Write common parameters */    
      fwrite(&i_nMCMC,sizeof(int),1,g_caa_mcmc2);
      fwrite(i_num_par2,sizeof(int),2,g_caa_mcmc2);
      fwrite(&i_print_boat,sizeof(int),1,g_caa_mcmc2);

      /* Write wgl model */
      write_glm_object(g_caa_mcmc2, i_D_wgl->glm, i_print_format);

      /* Write coastal cod parameters */
      fwrite(&i_coastal_cod,sizeof(int),1,g_caa_mcmc2);
      if(i_coastal_cod)
 	  write_glm_object(g_caa_mcmc2, i_D_wgl_CC->glm, i_print_format);
    }
  else //Ascii format
    {
      /* Write common parameters */    
      fprintf(g_caa_mcmc2,"%d\n",i_nMCMC);
      fprintf(g_caa_mcmc2,"%d %d\n",i_num_par2[0],i_num_par2[1]);
      fprintf(g_caa_mcmc2,"%d\n",i_print_boat);
      /* Write wgl model */
      write_glm_object(g_caa_mcmc2, i_D_wgl->glm, i_print_format);
      /* Write coastal cod parameters */
      fprintf(g_caa_mcmc2,"%d\n",i_coastal_cod);
      if(i_coastal_cod)
	write_glm_object(g_caa_mcmc2, i_D_wgl_CC->glm, i_print_format);
    }
  
  return(0);
}               /* end of write_mcmc2 */



/*!
  \brief Writes the containts of an LW-struct to file "filename" in the current directory.

  Only to be used for testing.
  \author Hanne Rognebakke
*/
int write_LW_struct(LW_struct *i_lin,Data_glm *i_glm,char *i_filename)
{
  int i,j,k,a;
  int err;
  Data_cov   *xcov;
  FILE *fp;
  
  fp = fopen(i_filename,"w");

  fprintf(fp,"LW_struct\n\n");
  fprintf(fp,"Eff_str:\n");  
  fprintf(fp,"ncat=%d\n",i_glm->ncat);
  for(a=0;a<i_glm->ncat;a++)
    {
      fprintf(fp,"nxcov=%d\n",i_glm->nxcov);
      for(i=0;i<i_glm->nxcov;i++)
	{
	  xcov = i_glm->xcov[i];
	  fprintf(fp,"  n_cov=%d\n",xcov->n_cov);
	  for(j=0;j<xcov->n_cov;j++)
	    {
	      for(k=0;k<xcov->n_fac[j];k++)
		fprintf(fp,"  eff[%d][%d][%d][%d]=%f\n",a,i,j,k,i_lin->par->eff[a][i][j][k]);
	      fprintf(fp,"  ssq[%d][%d][%d]=%f\n",a,i,j,i_lin->par->ssq[a][i][j]);
	      fprintf(fp,"  n_ssq[%d][%d][%d]=%d\n",a,i,j,i_lin->par->n_ssq[a][i][j]);
	    }
	}
    }
  for(i=0;i<i_glm->nxcov;i++)
    {
      fprintf(fp,"ar[%d]=%f\n",i,i_lin->par->ar[i]);
      xcov = i_glm->xcov[i];
      for(j=0;j<2;j++)
	fprintf(fp,"prior_ar[%d][%d]=%f\n",i,j,i_lin->par->prior_ar[i][j]);
      for(j=0;j<xcov->n_cov;j++)
	{
	  for(k=0;k<xcov->n_fac[j];k++)
	    for(a=0;a<i_glm->ncat;a++)
	      fprintf(fp,"prior_mean[%d][%d][%d][%d]=%f\n",a,i,j,k,i_lin->par->prior_mean[a][i][j][k]);
	  for(k=0;k<2;k++)
	    fprintf(fp,"prior_prec[%d][%d][%d]=%f\n",i,j,k,i_lin->par->prior_prec[i][j][k]);
	  fprintf(fp,"tau[%d][%d]=%f\n",i,j,i_lin->par->tau[i][j]);
	}
    }
  for(i=0;i<2;i++)
    fprintf(fp,"prior_prec_obs[%d]=%f\n",i,i_lin->par->prior_prec_obs[i]);

  fprintf(fp,"tau_obs=%f\n",i_lin->par->tau_obs);
  fprintf(fp,"loglik=%f\n",i_lin->par->loglik);
  fprintf(fp,"num_var=%d\n",i_lin->par->num_var);

  fprintf(fp,"\nGraph_str:\n");  
  for(i=0;i<i_glm->nxcov;i++)  
    {
      for(j=0;j<i_glm->xcov[i]->n_cov;j++)
	fprintf(fp,"in_gr[%d][%d]=%d\n",i,j,i_lin->gr_str->in_gr[i][j]);
    }
  err = GMRFLib_print_graph(fp,i_lin->gr_str->graph);


  fprintf(fp,"\ncens_model=%d\n",i_lin->cens_model);
  if(i_lin->cens_model)
    {
      fprintf(fp,"cens_k=%f,cens_m=%f,cens_r=%f,cens_Nlim=%f\n",
	      i_lin->cens_k,i_lin->cens_m,i_lin->cens_r,i_lin->cens_Nlim);
    }
  fprintf(fp,"\nfixed_model=%d\n",i_lin->fixed_model);

  fclose(fp);

  return(0);
}           /* end of write_LW_struct */

/*!
  \brief Writes the containts of Data_lin to file "filename" in the current directory.

  Only to be used for testing.
  \author Hanne Rognebakke
*/
int write_Data_lin(Data_lin *i_D_lga,int i_ncat,char *i_filename)
{
  int h,a;
  FILE *fp;
  fp=fopen(i_filename,"w");

  fprintf(fp,"Data_lin\n\n");

  fprintf(fp,"Data_glm:\n");
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      if(i_D_lga->glm->nxcov==3)
	fprintf(fp,"beta_hat[%d][0][0]=%f,beta_hat[%d][0][1]=%f,beta_hat[%d][0][2]=%f\n",
		h,i_D_lga->glm->beta_hat[h][0][0],h,i_D_lga->glm->beta_hat[h][0][1],
		h,i_D_lga->glm->beta_hat[h][0][2]);
      else
	fprintf(fp,"beta_hat[%d][0][0]=%f,beta_hat[%d][0][1]=%f\n",
		h,i_D_lga->glm->beta_hat[h][0][0],h,i_D_lga->glm->beta_hat[h][0][1]);
    }
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      if(i_D_lga->glm->nxcov==3)
	{
	  fprintf(fp,"suff[%d][0][0]=%f,suff[%d][0][1]=%f,suff[%d][0][2]=%f\n",
		  h,i_D_lga->glm->suff[h][0][0],h,i_D_lga->glm->suff[h][0][1],
		  h,i_D_lga->glm->suff[h][0][2]);
	  fprintf(fp,"suff[%d][1][0]=%f,suff[%d][1][1]=%f,suff[%d][1][2]=%f\n",
		  h,i_D_lga->glm->suff[h][1][0],h,i_D_lga->glm->suff[h][1][1],
		  h,i_D_lga->glm->suff[h][1][2]);
	  fprintf(fp,"suff[%d][2][0]=%f,suff[%d][2][1]=%f,suff[%d][2][2]=%f\n",
		  h,i_D_lga->glm->suff[h][2][0],h,i_D_lga->glm->suff[h][2][1],
		  h,i_D_lga->glm->suff[h][2][2]);
	}
      else
	{
	  fprintf(fp,"suff[%d][0][0]=%f,suff[%d][0][1]=%f\n",
		  h,i_D_lga->glm->suff[h][0][0],h,i_D_lga->glm->suff[h][0][1]);
	  fprintf(fp,"suff[%d][1][0]=%f,suff[%d][1][1]=%f\n",
		  h,i_D_lga->glm->suff[h][1][0],h,i_D_lga->glm->suff[h][1][1]);
	}
      fprintf(fp,"ssq[%d]=%f\n",h,i_D_lga->glm->ssq[h]);
    }

  fprintf(fp,"\n");
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      for(a=0;a<i_ncat;a++)
	{
	  fprintf(fp,"Ages[%d][%d]=%d,sum_by_cat[%d][%d]=%f,sqsum_by_cat[%d][%d]=%f\n",
		  h,a,i_D_lga->Ages[h][a],h,a,i_D_lga->sum_by_cat[h][a],
		  h,a,i_D_lga->sqsum_by_cat[h][a]);
	}
    }
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    fprintf(fp,"haulweight[%d]=%f\n",h,i_D_lga->haulweight[h]);
  for(h=0;h<i_D_lga->glm->nHaul;h++)
    {
      for(a=0;a<i_ncat;a++)
	{
	  fprintf(fp,"Ages_fix[%d][%d]=%d,sum_by_cat_fix[%d][%d]=%f,sqsum_by_cat_fix[%d][%d]=%f\n",
		  h,a,i_D_lga->Ages_fix[h][a],h,a,i_D_lga->sum_by_cat_fix[h][a],
		  h,a,i_D_lga->sqsum_by_cat_fix[h][a]);
	}
    }

  fclose(fp);

  return(0);
}            /* end of write_Data_lin */
