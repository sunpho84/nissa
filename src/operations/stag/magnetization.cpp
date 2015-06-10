#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "inverters/staggered/cg_invert_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //compute the magnetization starting from chi and rnd
  //please note that the conf must hold backfield and stagphases
  THREADABLE_FUNCTION_10ARG(magnetization, complex*,magn, complex*,magn_proj_x, quad_su3**,conf, quark_content_t*,quark, color**,rnd, color**,chi, complex*,point_magn, coords*,arg, int,mu, int,nu)
  {
    GET_THREAD_ID();
    
    communicate_ev_and_od_color_borders(chi);
    vector_reset(point_magn);
    
    //allocate a thread-local reduction
    complex thr_magn_proj_x[glb_size[1]];
    for(int i=0;i<glb_size[1];i++) thr_magn_proj_x[i][RE]=thr_magn_proj_x[i][IM]=0;
    
    //summ the results of the derivative
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
        {
          int ivol=loclx_of_loceo[par][ieo];
          
          //summ the contribution of the derivative in mu and nu directions
          int rho_list[2]={mu,nu};
          for(int irho=0;irho<2;irho++)
            {
              int rho=rho_list[irho];
              
              int iup_eo=loceo_neighup[par][ieo][rho];
              int idw_eo=loceo_neighdw[par][ieo][rho];
              int idw_lx=loclx_neighdw[ivol][rho];
              
              color v;
              complex t;
              
              //mu component of x
              int ix=glb_coord_of_loclx[ivol][1];
              
              //forward derivative
              unsafe_su3_prod_color(v,conf[par][ieo][rho],chi[!par][iup_eo]);
              color_scalar_prod(t,v,rnd[par][ieo]);
              complex_summ_the_prod_double(point_magn[ivol],t,arg[ivol][rho]);
              //compute also the projected current
              complex_summ_the_prod_double(thr_magn_proj_x[ix],t,arg[ivol][rho]);
              
              //backward derivative: note that we should multiply for -arg*(-U^+)
              unsafe_su3_dag_prod_color(v,conf[!par][idw_eo][rho],chi[!par][idw_eo]);
              color_scalar_prod(t,v,rnd[par][ieo]);
	      complex_summ_the_prod_double(point_magn[ivol],t,arg[idw_lx][rho]);
              //compute also the projected current
              complex_summ_the_prod_double(thr_magn_proj_x[ix],t,arg[idw_lx][rho]);
            }
        }
    THREAD_BARRIER();
    
    //reduce the projected magnetization
    complex temp_proj_x[glb_size[1]];
    for(int x=0;x<glb_size[1];x++) glb_reduce_complex(temp_proj_x[x],thr_magn_proj_x[x]);
    
    //reduce across all nodes and threads
    complex temp;
    complex_vector_glb_collapse(temp,point_magn,loc_vol);
    
    //add normalization, corresponding to all factors relative to derivative with respects to "b": 
    //-quark_deg/4 coming from the determinant
    //-1/vol coming from stochastic trace
    //-1/2 coming from dirac operator
    //-i*2*quark_charge*M_PI/glb_size[mu]/glb_size[nu] coming EM potential prefactor in front of "b"
    //and a minus because F=-logZ
    if(IS_MASTER_THREAD)
      {
        double coeff=-quark->deg*2*M_PI*quark->charge/(4.0*glb_vol*2*glb_size[mu]*glb_size[nu]);
        unsafe_complex_prod_idouble(*magn,temp,coeff);
        for(int x=0;x<glb_size[1];x++) unsafe_complex_prod_idouble(magn_proj_x[x],temp_proj_x[x],coeff);
      }
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //compute the magnetization
  THREADABLE_FUNCTION_8ARG(magnetization, complex*,magn, complex*,magn_proj_x, quad_su3**,conf, int,quantization, quad_u1**,u1b, quark_content_t*,quark, double,residue, color**,rnd)
  {
    GET_THREAD_ID();
    
    //fixed to Z magnetization
    int mu=1,nu=2;
    
    //allocate source and propagator
    color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
    
    //we need to store phases
    coords *arg=nissa_malloc("arg",loc_vol+bord_vol,coords);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol+bord_vol)
      get_args_of_quantization[quantization](arg[ivol],ivol,mu,nu);
    
    //array to store magnetization on single site (actually storing backward contrib at displaced site)
    complex *point_magn=nissa_malloc("app",loc_vol,complex);
    
    //we add stagphases and backfield externally because we need them for derivative
    addrem_stagphases_to_eo_conf(conf);
    add_backfield_to_conf(conf,u1b);
    
    //invert
    inv_stD_cg(chi,conf,quark->mass,100000,residue,rnd);
    
    //compute mag
    magnetization(magn,magn_proj_x,conf,quark,rnd,chi,point_magn,arg,mu,nu);
    
    //remove stag phases and u1 field
    rem_backfield_from_conf(conf,u1b);
    addrem_stagphases_to_eo_conf(conf);
    
    //free
    for(int par=0;par<2;par++) nissa_free(chi[par]);
    nissa_free(point_magn);
    nissa_free(arg);
  }
  THREADABLE_FUNCTION_END
  
  //compute the magnetization
  void magnetization(complex *magn,complex *magn_proj_x,quad_su3 **conf,int quantization,quad_u1 **u1b,quark_content_t *quark,double residue)
  {
    //allocate source and generate it
    color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
    generate_fully_undiluted_eo_source(rnd,RND_GAUSS,-1);
    
    //call inner function
    magnetization(magn,magn_proj_x,conf,quantization,u1b,quark,residue,rnd);

    for(int par=0;par<2;par++) nissa_free(rnd[par]);
  }
  
  //measure magnetization
  void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(theory_pars.magnetization_meas_pars.path,conf_created?"w":"a");
    FILE *file_proj=open_file(combine("%s_proj_x",theory_pars.magnetization_meas_pars.path).c_str(),conf_created?"w":"a");
    
    int ncopies=theory_pars.magnetization_meas_pars.ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      {
        master_fprintf(file,"%d",iconf);
        
        //measure magnetization for each quark
        for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
          {
            complex magn={0,0};
            complex magn_proj_x[glb_size[1]]; //this makes pair and pact with "1" and "2" upstairs
            for(int i=0;i<glb_size[1];i++) magn_proj_x[i][RE]=magn_proj_x[i][IM]=0;
            
            //loop over hits
            int nhits=theory_pars.magnetization_meas_pars.nhits;
            for(int hit=0;hit<nhits;hit++)
              {
                verbosity_lv2_master_printf("Evaluating magnetization for flavor %d/%d, ncopies %d/%d nhits %d/%d\n",
                                            iflav+1,theory_pars.nflavs,icopy+1,ncopies,hit+1,nhits);
            
                //compute and summ
                complex temp,temp_magn_proj_x[glb_size[1]];
                magnetization(&temp,temp_magn_proj_x,conf,theory_pars.em_field_pars.flag,theory_pars.backfield[iflav],theory_pars.quark_content+iflav,
                              theory_pars.magnetization_meas_pars.residue); //flag holds quantization
                
                //normalize
                complex_summ_the_prod_double(magn,temp,1.0/nhits);
                for(int x=0;x<glb_size[1];x++) complex_summ_the_prod_double(magn_proj_x[x],temp_magn_proj_x[x],1.0/nhits);
              }
            
            //output
            master_fprintf(file,"\t%+016.16lg \t%+016.16lg",magn[RE],magn[IM]);
            for(int x=0;x<glb_size[1];x++)
              master_fprintf(file_proj,"%d\t%d\t%d\t%d\t%+016.16lg \t%+016.16lg\n",iconf,icopy,iflav,x,magn_proj_x[x][RE],magn_proj_x[x][IM]);
          }
        
        master_fprintf(file,"\n");
      }
    
    close_file(file);
    close_file(file_proj);
  }
}
