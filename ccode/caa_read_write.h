void write_warning(char *i_text);
void write_output(char *filename, char *i_text);
int printWarning(char *i_text);
void printError(char *i_text, char *i_filename);
int write_par_model1(char *i_filename, Data_age *i_D_age, Age_struct *i_age,
		     Data_lin *i_D_lga, LW_struct *i_length, 
		     Data_lin *i_D_lga_CC, LW_struct *i_length_CC,  
		     Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, int i_coastal_cod, Data_CC *i_D_CC);
int write_par_bin(FILE *fp,Data_glm *i_glm,Eff_str *i_par);
int write_par_ascii(FILE *fp,Data_glm *i_glm,Eff_str *i_par);
int write_par_model2(char *i_filename, Data_lin *i_D_wgl, LW_struct *i_weight, 
		     Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, int i_coastal_cod);
int read_par_model1(char *i_filename, Data_age *i_D_age, Age_struct *i_age,
		    Data_lin *i_D_lga, LW_struct *i_length, 
		    Data_lin *i_D_lga_CC, LW_struct *i_length_CC, 
		    Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, int i_coastal_cod, Data_CC *i_D_CC);
int read_par_bin(FILE *fp,Data_glm *i_glm,Eff_str *i_par,int *i_nFac);
int read_complete_alpha_bin(FILE *fp,Data_glm *i_glm,Age_struct *i_age);
int read_par_ascii(FILE *fp,Data_glm *i_glm,Eff_str *i_par,double **i_alpha);
int read_par_model2(char *i_filename, Data_lin *i_D_wgl, LW_struct *i_weight, 
		    Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC, int i_coastal_cod);
int write_mcmc1(Data_age *i_D_age, Age_struct *i_age, 
		Data_lin *i_D_lga, LW_struct *i_length, Data_lin *i_D_lga_CC,  
		LW_struct *i_length_CC, Data_g_a *i_D_g_a, Data_g_a *i_D_g_a_CC, 
		int i_nMCMC, int *i_num_par1, int i_coastal_cod, int i_print_boat,
		double i_delta_age, int i_print_format);
int write_glm_object(FILE *fp, Data_glm *i_glm, int i_print_format);
int read_glm_object(FILE *fp, Data_glm **o_glm);
int read_mcmc1(Input_predict *i_inPredict, Data_age **o_D_age, Data_lin **o_D_lga, Data_g_a **o_D_g_a,
	       Data_lin **o_D_lga_CC, Data_g_a **o_D_g_a_CC);
int read_mcmc2(Input_predict *i_inPredict, Data_lin **o_D_wgl, Data_lin **o_D_wgl_CC);
int read_hsz(Input_predict *i_inPredict,  Data_lin **o_D_wgl);
int write_mcmc2(Data_lin *i_D_wgl, LW_struct *i_weight, 
		Data_lin *i_D_wgl_CC, LW_struct *i_weight_CC,
		int i_nMCMC, int *i_num_par2, int i_coastal_cod, int i_print_boat, int i_print_format);
int write_LW_struct(LW_struct *i_lin, Data_glm *i_glm, char *i_filename);
int write_Data_lin(Data_lin *i_D_lga, int i_ncat, char *i_filename);
int read_par_ascii_lga(FILE *fp,Data_glm *i_glm,Eff_str *i_par);
