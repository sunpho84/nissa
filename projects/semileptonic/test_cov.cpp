#include <nissa.h>

int T=8,L=4;
int muS=1,tsource=0;
double kappa=0.177,mass=0.50,residue=1e-24;
char conf_path[]="test/data/L4T8conf";

void apply_nabla(su3spinspin *out,su3spinspin *in,quad_su3 *conf,int mu)
{
  spincolor *temp_out=nissa_malloc("temp_out",loc_vol,spincolor);
  spincolor *temp_in=nissa_malloc("temp_in",loc_vol,spincolor);
  
  for(int id=0;id<4;id++)
    for(int ic=0;ic<3;ic++)
      {
	get_spincolor_from_su3spinspin(temp_in,in,id,ic);
	apply_nabla_i(temp_out,temp_in,conf,mu);
	put_spincolor_into_su3spinspin(out,temp_out,id,ic);
      }
  
  nissa_free(temp_out);
  nissa_free(temp_in);
}

//invert the prop
void get_prop(su3spinspin *out,quad_su3 *conf,int cov)
{
  //source
  su3spinspin *source=nissa_malloc("source",loc_vol+bord_vol,su3spinspin);
  coords O={tsource,0,0,0};
  generate_delta_source(source,O);
  
  //allocate solution
  su3spinspin **solution[2]={nissa_malloc("a",1,su3spinspin*),nissa_malloc("a",1,su3spinspin*)};
  for(int r=0;r<2;r++) solution[r][0]=nissa_malloc("solution",loc_vol+bord_vol,su3spinspin);
  
  //prepare the source applying der if needed
  if(cov) apply_nabla_i(source,source,conf,muS);

  //invert and put in place
  int r=1;
  compute_su3spinspin_tm_propagators_multi_mass((su3spinspin***)solution,conf,kappa,&mass,1,10000,&residue,source);
  vector_copy(out,solution[r][0]);
  rotate_vol_su3spinspin_to_physical_basis(out,!r,!r);
  
  nissa_free(source);
  nissa_free(solution);
}

//compute correlators and print to stdout
void compute_correlators(int *op_source,int *op_sink,su3spinspin *prop1,su3spinspin *prop2,const char *tag)
{
  int ncontr=1;
  complex corr[T];
  meson_two_points_Wilson_prop(corr,op_source,prop1,op_sink,prop2,1);
  print_contractions_to_file(stdout,ncontr,op_source,op_sink,corr,tsource,tag,1.0);
}

void in_main(int narg,char **arg)
{
  init_grid(T,L); 
  
  const int LOCAL=0,WITH_DER=1;
  
  //conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  read_ildg_gauge_conf(conf,conf_path);
  /*
  start_loc_rnd_gen(100);
  perform_random_gauge_transform(conf,conf);  
  */
  
  master_printf("plaq: %.18g\n",global_plaquette_lx_conf(conf));
  double put_theta[4]={1,0,0,0},old_theta[4]={0,0,0,0};
  adapt_theta(conf,old_theta,put_theta,1,1);
  
  //prop from local source
  su3spinspin *loc_prop=nissa_malloc("loc_prop",loc_vol+bord_vol,su3spinspin);
  get_prop(loc_prop,conf,LOCAL);
  
  //propagator from der
  su3spinspin *der_prop=nissa_malloc("der_prop",loc_vol+bord_vol,su3spinspin);
  get_prop(der_prop,conf,WITH_DER);
  
  //derivative in the sink for the local source
  su3spinspin *loc_prop_der=nissa_malloc("loc_prop_der",loc_vol,su3spinspin);
  apply_nabla_i(loc_prop_der,loc_prop,conf,muS);
	
  //derivative in the sink for the covariant source
  su3spinspin *der_prop_der=nissa_malloc("der_prop_der",loc_vol,su3spinspin);
  apply_nabla_i(der_prop_der,der_prop,conf,muS);
  
  //operators
  int op_l[1]={5};
  int op_d[1]={6};
  
  compute_correlators(op_l,op_l,loc_prop,loc_prop,"LSOU-LSIN-");
  compute_correlators(op_d,op_l,loc_prop,der_prop,"DSOU-LSIN-");
  compute_correlators(op_l,op_d,loc_prop,loc_prop_der,"LSOU-DSIN-");
  compute_correlators(op_d,op_d,loc_prop,der_prop_der,"DSOU-DSIN-");
  
  nissa_free(conf);
  nissa_free(loc_prop);
  nissa_free(der_prop);
  nissa_free(loc_prop_der);
  nissa_free(der_prop_der);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
    
  return 0;
}
