#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/macros.hpp"
#include "base/thread_macros.hpp"
#include "geometry/geometry_mix.hpp"
#include "io/input.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "operations/gauge_fixing.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "arbitrary.hpp"

#include "flux_tube.hpp"

namespace nissa
{
  //compute 1/4 of the watusso-path
  void watusso_quarter(path_drawing_t *c,su3 *out,quad_su3 *conf,int L,int D,int mu,int nu,int mu_orie,int nu_orie,int rh,int si)
  {
    int steps[2*10]={
      nu,-nu_orie*L/2,
      si,-D,
      mu,-mu_orie*L/2,
      nu,+nu_orie*L/2,
      si,+D,
      mu,+mu_orie*L/2,
      rh,-nu_orie,
      mu,-mu_orie,
      rh,+nu_orie,
      mu,+mu_orie
    };
    
    elong_su3_path(c,out,conf,steps,10);
  }
  
  //compute the flux-tube in the watusso way
  void watusso_path(path_drawing_t *c,su3 *out,quad_su3 *conf,int L,int D,int mu,int nu,int rh,int si)
  {
    init_su3_path(c,out);
    
    //perform the four sub-paths
    watusso_quarter(c,out,conf,L,D,mu,nu,+1,+1,rh,si);
    watusso_quarter(c,out,conf,L,D,nu,mu,+1,-1,mu,si);
    watusso_quarter(c,out,conf,L,D,mu,nu,-1,-1,rh,si);
    watusso_quarter(c,out,conf,L,D,nu,mu,-1,+1,mu,si);
  }
  
  //compute the flux tube for a fixed size and various distances along a fixed direction
  //[mu,nu] define the plane on which the large Wilson loop is lying
  //[mu,rho] define the probe-plaquette
  //sigma defines the axis (that should be perpendicular to mu,nu and rho)
  THREADABLE_FUNCTION_8ARG(compute_flux_tube, double*,out, quad_su3*,conf, int,L, int,D, int,mu, int,nu, int,rh, int,si)
  {
    GET_THREAD_ID();

    //check that the size is even
    if(L%2) crash("implemented only for even sizes, but L=%d",L);

    //for debug purpose
    path_drawing_t *c=new path_drawing_t;
    
    //compute the path
    su3 *Wloop=nissa_malloc("Wloop",loc_vol+bord_vol,su3);
    watusso_path(c,Wloop,conf,L,D,mu,nu,rh,si);

    //trace it
    complex *point_diag=nissa_malloc("point_diag",loc_vol,complex);
    complex temp_out={0,0};
    for(int ic=0;ic<3;ic++)
      {
	//copy
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  complex_copy(point_diag[ivol],Wloop[ivol][ic][ic]);
	THREAD_BARRIER();

	//reduce
	complex diag_ic;
	complex_vector_glb_collapse(diag_ic,point_diag,loc_vol);
	complex_summassign(temp_out,diag_ic);
      }
    
    //free
    nissa_free(point_diag);
    nissa_free(Wloop);
    
    //copy out
    complex_copy(out,temp_out);

    //print the path
    FILE *fout=open_file("path_draw","w");
    master_fprintf(fout,"set terminal x11\n"
		   "set pm3d explicit at s depthorder hidden3d 1\n"
		   "set hidden3d front\n"
		   "set style fill transparent solid 0.65\n"
		   "set palette rgb 9,9,3\n"
		   "set xrange [-5:5]\n"
		   "set yrange [-5:5]\n"
		   "set zrange [-2:5]\n"
		   "set isosamples 2,2\n"
		   "unset key\n"
		   "unset tics\n"
		   "unset border\n"
		   "unset colorbox\n"
		   "set view 60,30,1.,1\n");
    for(int i=1;i<(int)c->size();i++)
      {
	master_fprintf(fout,"pause 1\n");
	master_fprintf(fout,"set arrow %d from %d,%d,%d to %d,%d,%d lw 1 head front\n",i,
		       (*c)[i-1][mu],(*c)[i-1][nu],(*c)[i-1][si],
		       (*c)[i][mu],(*c)[i][nu],(*c)[i][si]);
	master_fprintf(fout,"splot 0 w pm3d\n");
      }
    
    delete c;
  }
  THREADABLE_FUNCTION_END

  //compute using the pars specified in the class
  void compute_flux_tube(quad_su3 *conf,flux_tube_meas_pars_t &ft)
  {
    //compute 
    complex out;
    compute_flux_tube((double*)out,conf,ft.L,ft.D,ft.mu,ft.nu,ft.rh,ft.si);
    master_printf("%lg %lg\n",out[0],out[1]);

    //make a random gauge transofrmation to check gauge invariance
    complex out_gt;
    if(ft.gauge_check)
      {
	//perform a random gauge transformation
	quad_su3 *conf_gt=nissa_malloc("conf_gt",loc_vol,quad_su3);
	perform_random_gauge_transform(conf_gt,conf);

	//recompute
	compute_flux_tube((double*)out_gt,conf_gt,ft.L,ft.D,ft.mu,ft.nu,ft.rh,ft.si);
	master_printf("%lg %lg\n",out_gt[0],out_gt[1]);
	nissa_free(conf_gt);
      }
    
    if(fabs(out_gt[RE]-out[RE])>1.e-12||fabs(out_gt[IM]-out[IM])>1.e-12) crash("flux: (%lg,%lg) differs after gauge transformation: (%lg,%lg)",out[RE],out[IM],out_gt[RE],out_gt[IM]);
  }
  void compute_flux_tube(quad_su3 **eo_conf,flux_tube_meas_pars_t &ft)
  {
    quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol,quad_su3);
    paste_eo_parts_into_lx_conf(lx_conf,eo_conf);
    compute_flux_tube(lx_conf,ft);
    nissa_free(lx_conf);
  }
  
  //read pars
  void read_flux_tube_meas_pars(flux_tube_meas_pars_t &ft,bool flag)
  {
    if(flag==true) ft.flag=true;
    else read_str_int("MeasureFluxTube",&ft.flag);
    if(ft.flag)
      {
	read_str_str("Path",ft.path,100);
	read_str_int("L",&ft.L);
	read_str_int("D",&ft.D);
        ft.mu=0;
	ft.nu=1;
        ft.rh=1;
	ft.si=2;
	ft.gauge_check=true;
      }
  }
}
