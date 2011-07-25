as2t_su3 *allocate_as2t_su3(int length,char *tag);
char *allocate_vector(int length,char *tag);
colorspinspin *allocate_colorspinspin(int length,char *tag);
double calculate_weighted_residue_RL(spincolor *source,spincolor *sol,quad_su3 *conf,double kappa,double m,spincolor *s,spincolor *t,int dinf,int RL);
double global_plaquette(quad_su3 *conf);
double max_double(double a,double b);
double min_double(double a,double b);
double ran2(int loc_ind);
double su3_normq(su3 U);
double take_time();
FILE* open_text_file_for_output(char *outfile);
int bordlx_of_coord(int *x,int idir);
int bordlx_of_coord_list(int x0,int x1,int x2,int x3,int idir);
int check_cgmms_residue_RL(int *run_flag,double *residue_mass,int nrun,double rr,double *zfs,int st_crit,double st_res,double st_res2,int iter,spincolor **sol,int nmass,double *m,spincolor *source,quad_su3 *conf,double kappa,spincolor *s,spincolor *t,int RL);
int edgelx_of_coord(int *x,int idir,int jdir);
int glblx_of_coord(int *x);
int loclx_of_coord(int *x);
int loclx_of_coord_list(int x0,int x1,int x2,int x3);
int lx_of_coord(int *x,int *s);
int main();
int max_int(int a,int b);
int min_int(int a,int b);
int pm_one(int loc_site);
int read_binary_blob(char *data_out,char *path,const char *expected_record,int max_nbytes_per_site);
int trivialize_dirac_matr(short int *gt,dirac_matr *g);
quad_su3 *allocate_quad_su3(int length,char *tag);
redspincolor *allocate_redspincolor(int length,char *tag);
spincolor *allocate_spincolor(int length,char *tag);
su3 *allocate_su3(int length,char *tag);
su3spinspin *allocate_su3spinspin(int length,char *tag);
void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis);
void adapt_theta(quad_su3 *conf,double *old_theta,double *put_theta,int putonbords,int putonedges);
void ape_smearing(quad_su3 *smear_conf,quad_su3 *origi_conf,double alpha,int nstep);
void apply_lxspincolor_border_projector(redspincolor *p,spincolor *s);
void apply_Q2_left(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,spincolor *temp,redspincolor *tin,redspincolor *tout);
void apply_Q2_RL(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,spincolor *temp,redspincolor *tin,redspincolor *tout,int RL);
void apply_Q2(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,spincolor *temp,redspincolor *tin,redspincolor *tout);
void apply_Q_left(spincolor *out,spincolor *in,quad_su3 *conf,double kappac,double mu);
void apply_Q_RL(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,double mu,int RL);
void apply_Q(spincolor *out,spincolor *in,quad_su3 *conf,double kappac,double mu);
void apply_Q_v0(spincolor *out,spincolor *in,quad_su3 *conf,double kappac,double mu);
void apply_Q_v1(spincolor *out,spincolor *in,quad_su3 *conf,double kappac,double mu);
void as2t_saturate(complex out,as2t a,as2t b);
void as2t_su3_put_to_zero(as2t_su3 m);
void assign_complex_prod_i(complex a);
void assign_complex_prod_minus_i(complex a);
void assign_spincolor_prod_real(spincolor out,double factor);
void check_endianess();
void check_free(void *a);
void close_appretto();
void close_input();
void close_random();
void color_copy(color b,color a);
void color_isubt(color a,color b,color c);
void color_isumm(color a,color b,color c);
void color_prod_real(color a,color b,double c);
void color_put_to_zero(color m);
void color_subt(color a,color b,color c);
void color_summ(color a,color b,color c);
void communicate_gauge_borders(quad_su3 *conf);
void communicate_gauge_edges(quad_su3 *conf);
void communicate_lx_borders(char *data,MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *MPI_BORD_RECE,int nbytes_per_site);
void communicate_lx_edges(char *data,MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,int nbytes_per_site);
void communicate_lx_spincolor_borders(spincolor *s);
void communicate_su3_borders(su3 *u);
void complex_pow(complex res,complex base,double exp);
void complex_prod_real(complex a,complex b,double c);
void complex_prod_with_real(complex a,complex b,double c);
void complex_reciprocal(complex rec,complex c);
void complex_subtassign(complex a,complex b);
void complex_subt(complex a,complex b,complex c);
void complex_subt_the_conj1_prod(complex a,complex b,complex c);
void complex_subt_the_conj1_prod_i(complex a,complex b,complex c);
void complex_subt_the_conj2_prod(complex a,complex b,complex c);
void complex_subt_the_conj2_prod_i(complex a,complex b,complex c);
void complex_subt_the_conj_conj_prod(complex a,complex b,complex c);
void complex_subt_the_conj_conj_prod_i(complex a,complex b,complex c);
void complex_subt_the_prod(complex a,complex b,complex c);
void complex_subt_the_prod_i(complex a,complex b,complex c);
void complex_summassign(complex a,complex b);
void complex_summ(complex a,complex b,complex c);
void complex_summ_the_conj1_prod(complex a,complex b,complex c);
void complex_summ_the_conj1_prod_i(complex a,complex b,complex c);
void complex_summ_the_conj2_prod(complex a,complex b,complex c);
void complex_summ_the_conj2_prod_i(complex a,complex b,complex c);
void complex_summ_the_conj_conj_prod(complex a,complex b,complex c);
void complex_summ_the_conj_conj_prod_i(complex a,complex b,complex c);
void complex_summ_the_prod(complex a,complex b,complex c);
void complex_summ_the_prod_i(complex a,complex b,complex c);
void density_profile(double *glb_rho,spincolor *sp,int *or_pos);
void dirac_prod(dirac_matr *out,dirac_matr *in1,dirac_matr *in2);
void dirac_summ(dirac_matr *out,dirac_matr *in1,dirac_matr *in2);
void doubles_to_doubles_changing_endianess(double *dest,double *sour,int ndoubles);
void doubles_to_floats_changing_endianess(float *dest,double *sour,int n);
void doubles_to_floats_same_endianess(float *dest,double *sour,int n);
void expect_str(const char *exp_str);
void find_temporal_gauge_fixing_matr(su3 *fixm,quad_su3 *u);
void floats_to_doubles_changing_endianess(double *dest,float *sour,int n);
void floats_to_doubles_same_endianess(double *dest,float *sour,int n);
void free_null(void *a);
void gauge_transform_conf(quad_su3 *uout,su3 *g,quad_su3 *uin);
void get_spincolor_from_colorspinspin(spincolor out,colorspinspin in,int id_source);
void get_spincolor_from_su3spinspin(spincolor out,su3spinspin in,int id_source,int ic_source);
void init_appretto();
void init_base_gamma();
void init_dirac(dirac_matr *out,int pos0,double rea0,double ima0,int pos1,double rea1,double ima1,int pos2,double rea2,double ima2,int pos3,double rea3,double ima3);
void init_grid();
void initialize_e_bord_senders_of_kind(MPI_Datatype *MPI_E_BORD_SEND,MPI_Datatype *base);
void initialize_lx_bord_receivers_of_kind(MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base);
void initialize_lx_bord_senders_of_kind(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *base);
void initialize_lx_edge_receivers_of_kind(MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
void initialize_lx_edge_senders_of_kind(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *base);
void init_random(int seed);
void inv_Q2_cg_left(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue);
void inv_Q2_cgmms_left(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit);
void inv_Q2_cgmms_RL(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit,int RL);
void inv_Q2_cgmms(spincolor **sol,spincolor *source,spincolor **guess,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit);
void inv_Q2_cg_RL(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue,int RL);
void inv_Q2_cg_sorc_RL(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue,int RL);
void inv_Q2_cg(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue);
void jacobi_smearing(spincolor *smear_sc,spincolor *origi_sc,quad_su3 *conf,double kappa,int niter);
void Momentum(int **iP,double *bc,double *P2,double *SinP2,double **P,double **SinP,double *SinP4,int nmom);
void open_input(char *input_path);
void Pline_backward(su3 *Pline, quad_su3 *conf);
void Pline_forward(su3 *Pline, quad_su3 *conf);
void Pline(su3 *Pline,quad_su3 *conf);
void Pmunu_term(as2t_su3 *Pmunu,quad_su3 *conf);
void print_contractions_to_file(FILE *fout,int ncontr,int *op1,int *op2,complex *contr,int twall,const char *tag);
void print_contraction_to_file(FILE *fout,int op1,int op2,complex *contr,int twall,const char *tag,double norm);
void print_dirac(dirac_matr *in);
void print_spinspin(spinspin s);
void put_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
void put_spincolor_into_colorspinspin(colorspinspin out,spincolor in,int id_source);
void put_spincolor_into_su3spinspin(su3spinspin out,spincolor in,int id_source,int ic_source);
void quad_su3_copy(quad_su3 b,quad_su3 a);
void read_colorspinspin(colorspinspin *css,char *base_path,char *end_path);
void read_colorspinspin_reconstructing(colorspinspin **css,char *base_path,char *end_path,quad_su3 *conf,double kappa,double mu);
void read_double(double *in);
void read_int(int *in);
void read_list_of_doubles(char *tag,int *nentries,double **list);
void read_list_of_ints(char *tag,int *nentries,int **list);
void read_list_of_var(char *tag,int *nentries,char **list,int size_of_el,const char *par);
void read_local_gauge_conf(quad_su3 *out,char *path);
void read_real_vector(double *out,char *path,const char *header,int nreals_per_site);
void read_spincolor_reconstructing(spincolor **out,spincolor *temp,char *path,quad_su3 *conf,double kappa,double mu);
void read_spincolor(spincolor *out,char *path);
void read_str(char *str,int length);
void read_str_double(const char *exp_str,double *in);
void read_str_int(const char *exp_str,int *in);
void read_str_str(const char *exp_str,char *in,int length);
void read_var(char *in,const char *par,int size_of);
void reconstruct_doublet(spincolor *outminus,spincolor *outplus,spincolor *in,quad_su3 *conf,double kappac,double mu);
void rem_boundaries_conditions(quad_su3 *conf,double *theta_in_pi,int putonbords,int putonedges);
void rotate_spinspin_to_physical_basis(spinspin s,int rsi,int rso);
void rotate_vol_colorspinspin_to_physical_basis(colorspinspin *s,int rsi,int rso);
void safe_complex_conj1_prod(complex a,complex b,complex c);
void safe_complex_conj2_prod(complex a,complex b,complex c);
void safe_complex_prod(complex a,complex b,complex c);
void safe_complex_prod_i(complex a,complex b);
void safe_complex_prod_minus_i(complex a,complex b);
void safe_complex_prod_spincolor(spincolor out,complex a,spincolor in);
void safe_dirac_compl_prod(dirac_matr *out,dirac_matr *in,complex c);
void safe_dirac_prod_spincolor(spincolor out,dirac_matr *m,spincolor in);
void safe_spincolor_prod_dirac(spincolor out,spincolor in,dirac_matr *m);
void safe_spincolor_summ_with_cfactor(spincolor a,spincolor b,spincolor c,complex factor);
void safe_spincolor_summ_with_rfactor(spincolor a,spincolor b,spincolor c,double factor);
void safe_spinspin_prod_dirac(spinspin out,spinspin in,dirac_matr *m);
void safe_su3_dag_prod_color(color a,su3 b,color c);
void safe_su3_hermitian(su3 out,su3 in);
void safe_su3_prod_color(color a,su3 b,color c);
void safe_su3_prod_complex(su3 a,su3 b,complex c);
void set_eo_geometry();
void set_lx_bord_senders_and_receivers(MPI_Datatype *MPI_BORD_SEND,MPI_Datatype *MPI_BORD_RECE,MPI_Datatype *base);
void set_lx_edge_senders_and_receivers(MPI_Datatype *MPI_EDGE_SEND,MPI_Datatype *MPI_EDGE_RECE,MPI_Datatype *base);
void set_lx_geometry();
void shift_gauge_conf_down(quad_su3 *conf,int *amount);
void site_trace_g_sdag_g_s(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2);
void site_trace_g_sdag_g_s_g_sdag_g_s(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2, dirac_matr *g3,spinspin s3,dirac_matr *g4,spinspin s4);
void smearing_apply_kappa_H(spincolor *H,double kappa,quad_su3 *conf,spincolor *smear_sc);
void spincolor_copy(spincolor b,spincolor a);
void spincolor_FT(spincolor *S,spincolor *FT,double *theta,int **iP,int nmom);
void spincolor_put_to_zero(spincolor m);
void spincolor_subt(spincolor a,spincolor b,spincolor c);
void spincolor_summ(spincolor a,spincolor b,spincolor c);
void spincolor_summ_the_prod_complex(spincolor out,spincolor in,complex factor);
void spinspin_dirac_spinspindag_prod(spinspin out,dirac_matr *m,spinspin in);
void spinspin_dirac_spinspin_prod(spinspin out,dirac_matr *m,spinspin in);
void spinspin_spinspindag_prod(spinspin out,spinspin a,spinspin b);
void spinspin_spinspin_prod(spinspin out,spinspin a,spinspin b);
void squared_path(su3 square,quad_su3 *conf,int A,int mu,int nu);
void start_communicate_lx_redspincolor_borders(redspincolor *in,redspincolor *out,MPI_Request *request);
void su3_copy(su3 b,su3 a);
void su3_dag_prod_su3_dag(su3 a,su3 b,su3 c);
void su3_dag_prod_su3(su3 a,su3 b,su3 c);
void su3_dag_summ_the_color_prod(color a,su3 b,color c);
void su3_det(complex d,su3 U);
void su3_explicit_inverse(su3 invU,su3 U);
void su3_print(su3 U);
void su3_prod_real(su3 a,su3 b,double r);
void su3_prod_su3_dag(su3 a,su3 b,su3 c);
void su3_prod_su3(su3 a,su3 b,su3 c);
void su3_put_to_id(su3 m);
void su3_put_to_zero(su3 m);
void su3spinspin_put_to_zero(su3spinspin m);
void su3_subt_complex(su3 a,su3 b,complex c);
void su3_subt(su3 a,su3 b,su3 c);
void su3_summ_real(su3 a,su3 b,double c);
void su3_summ(su3 a,su3 b,su3 c);
void su3_summ_the_color_prod(color a,su3 b,color c);
void su3_trace(complex tr,su3 m);
void su3_unitarize(su3 new_link,su3 prop_link);
void subtassign_color(color a,color b);
void subtassign_icolor(color a,color b);
void subtassign_spincolor(spincolor out,spincolor in);
void summassign_color(color a,color b);
void summassign_icolor(color a,color b);
void summassign_spincolor(spincolor out,spincolor in);
void sum_trace_g_sdag_g_s_times_trace_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr);
void sum_trace_id_sdag_g_s_times_trace_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr);
void swap_doubles(double *d1,double *d2);
void swap_eo_to_lx(char *out_lx,char *in_e,char *in_o,int nbytes_per_site,int bord);
void swap_lx_to_eo(char *out_e,char *out_o,char *in_lx,int nbytes_per_site,int bord);
void swap_lx_to_eo_or_eo_to_lx(char *vect_e,char *vect_o,char *vect_lx,int nbytes_per_site,int bord,int eotolx_lxtoeo);
void swap_spincolor_eo_to_lx(spincolor *out_lx,spincolor *in_e,spincolor *in_o,int bord);
void swap_spincolor_lx_to_eo(spincolor *out_e,spincolor *out_o,spincolor *in_lx,int bord);
void trace_g_sdag_g_s(complex *glb_c,dirac_matr *g1,colorspinspin *s1,dirac_matr *g2,colorspinspin *s2,const int ncontr);
void trace_g_sdag_g_s_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr);
void trace_id_sdag_g_s_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr);
void trace_prod_dirac_spinspin(complex c,dirac_matr *a,spinspin b);
void trace_prod_spinspins(complex c,spinspin a,spinspin b);
void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,as2t_su3 *Pmunu,colorspinspin *in);
void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,as2t_su3 *Pmunu,spincolor *in);
void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,as2t_su3 *Pmunu,su3spinspin *in);
void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,as2t_su3 Pmunu,spincolor in);
void unsafe_color_prod_su3(color a,color b,su3 c);
void unsafe_color_prod_su3_dag(color a,color b,su3 c);
void unsafe_complex_conj1_prod(complex a,complex b,complex c);
void unsafe_complex_conj2_prod(complex a,complex b,complex c);
void unsafe_complex_conj_conj_prod(complex a,complex b,complex c);
void unsafe_complex_prod(complex a,complex b,complex c);
void unsafe_dirac_compl_prod(dirac_matr *out,dirac_matr *in,complex c);
void unsafe_dirac_prod_spincolor(spincolor out,dirac_matr *m,spincolor in);
void unsafe_spincolor_prod_complex(spincolor out,spincolor in,complex factor);
void unsafe_spincolor_prod_dirac(spincolor out,spincolor in,dirac_matr *m);
void unsafe_spincolor_summ_with_ifactor(spincolor out,spincolor a,spincolor b,double factor);
void unsafe_su3_dag_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in);
void unsafe_su3_dag_prod_color(color a,su3 b,color c);
void unsafe_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in);
void unsafe_su3_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in);
void unsafe_su3_hermitian(su3 out,su3 in);
void unsafe_su3_prod_color(color a,su3 b,color c);
void unsafe_su3_prod_spincolor(spincolor out,su3 U,spincolor in);
void unsafe_su3spinspin_prod_complex(su3spinspin out,su3spinspin in,complex factor);
void unsafe_subt_su3_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in);
void unsafe_subt_su3_prod_color(color a,su3 b,color c);
void unsafe_summassign_spincolor_prod_ireal(spincolor out,spincolor in,double factor);
void unsafe_summassign_spincolor_prod_real(spincolor out,spincolor in,double factor);
void unsafe_summ_su3_dag_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in);
void unsafe_summ_su3_dag_prod_color(color a,su3 b,color c);
void unsafe_summ_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in);
void unsafe_summ_su3_prod_color(color a,su3 b,color c);
void unsafe_summ_su3_prod_spincolor(spincolor out,su3 U,spincolor in);
void vol_assign_spincolor_prod_real(spincolor *sc,double c);
void vol_spincolor_summassign(spincolor *smear_sc,spincolor *H);
void write_double_vector(LemonWriter *writer,char *data,char *header_message,int nreals_per_site,int nbits);
void write_header(LemonWriter *writer,char *header,uint64_t record_bytes);
void write_local_gauge_conf(char *path,quad_su3 *in);
void write_spincolor(char *path,spincolor *spinor,int prec);
void write_text_record(LemonWriter *writer,char *header,char *message);
