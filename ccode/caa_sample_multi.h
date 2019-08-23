int sample_multi_initialize(int i_ncat);
int sample_multi_re_initialize();
int age_haul_modes(int i_start_h,int i_stop_h,Age_struct *i_age,Data_age *i_D_age);
int sample_ages_init(Data_orig *i_D_orig,Data_CC *i_D_CC,Data_age *i_D_age,
		     Data_lin *i_D_lga,Data_g_a *i_D_g_a,
		     Data_lin *i_D_lga_CC,Data_g_a *i_D_g_a_CC,int i_saveSim);
int sample_ages(Data_orig *i_D_orig,Data_CC *i_D_CC,
		Age_struct *i_age,Data_age *i_D_age,
		LW_struct *i_length,Data_lin *i_D_lga,Data_g_a *i_D_g_a, 
		LW_struct *i_length_CC,Data_lin *i_D_lga_CC,Data_g_a *i_D_g_a_CC, 
		int i_saveSim, int i_it);
int init_suff_stat_sim_age(int i_h,Data_age *i_D_age,Data_lin *i_D_lga,Data_lin *i_D_lga_CC,
			   int i_ga_ncat,int i_sim,int i_cc);
int calc_prob_age(Data_CC *i_D_CC,int i_cc,Age_struct *i_age,int i_aobs,
		  double i_lobs,double i_lstart,double i_lend,int i_type,int i_ncat_age,
		  double *i_mu,double *i_sigma,double *i_p,double *o_p2,double *o_sum_p);
int calc_int_slp(int i_h, int i_season, int i_ncat_age, Data_glm *i_glm, Eff_str *i_par, Data_g_a *i_D_g_a,
		 int i_cc, Data_glm *i_glmCC, Eff_str *parCC, Data_g_a *i_D_g_a_CC, double *o_mu);
int make_suff_age(int i_ncat,Age_struct *i_age,Data_age *i_D_age,double *i_haulweight,
		  int i_start_h);
int make_suff_lga(Data_lin *i_D_lga,Data_g_a *i_g_a,int i_start_h,int i_use_sim_age);
int copy_suff_lga_fix(Data_lin *i_D_lga,int i_start_h);
int sample_age_alpha(Age_struct *i_age,Data_age *i_D_age,int i_start_h,int i_acc,
		     int i_it,int *acc_h,int i_write_alpha);
int sample_precision_age_haul(int i_start_h,Eff_str *i_par,
			      double **i_alpha,Data_glm *i_glm,int i_nHaul,int i_nconstr);
