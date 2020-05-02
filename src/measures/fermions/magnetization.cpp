#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "inverters/staggered/cg_invert_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "operations/gaugeconf.hpp"
#include "routines/mpi_routines.hpp"

#include "magnetization.hpp"

namespace nissa
{
  //compute the magnetization starting from chi and rnd
  //please note that the conf must hold backfield
  THREADABLE_FUNCTION_9ARG(magnetization, complex*,magn, quad_su3**,conf, quark_content_t*,quark, color**,rnd, color**,chi, complex*,point_magn, coords*,arg, int,mu, int,nu)
  {
    GET_THREAD_ID();
    
    communicate_ev_and_od_color_borders(chi);
    vector_reset(point_magn);
    
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
              
	      //forward derivative
              unsafe_su3_prod_color(v,conf[par][ieo][rho],chi[!par][iup_eo]);
              color_scalar_prod(t,rnd[par][ieo],v);
              complex_summ_the_prod_double(point_magn[ivol],t,arg[ivol][rho]);
              
              //backward derivative: note that we should multiply for -arg*(-U^+)
              unsafe_su3_dag_prod_color(v,conf[!par][idw_eo][rho],chi[!par][idw_eo]);
              color_scalar_prod(t,rnd[par][ieo],v);
	      complex_summ_the_prod_double(point_magn[ivol],t,arg[idw_lx][rho]);
            }
        }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
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
        complex_prod_idouble(*magn,temp,coeff);
      }
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END
  
  //compute the magnetization
  THREADABLE_FUNCTION_7ARG(magnetization, complex*,magn, quad_su3**,conf, int,quantization, quad_u1**,u1b, quark_content_t*,quark, double,residue, color**,rnd)
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
    NISSA_PARALLEL_LOOP_END;
    
    //array to store magnetization on single site (actually storing backward contrib at displaced site)
    complex *point_magn=nissa_malloc("app",loc_vol,complex);
    
    //we add backfield externally because we need them for derivative
    add_backfield_with_stagphases_to_conf(conf,u1b);
    
    //invert
    inv_stD_cg(chi,conf,quark->mass,1000000,residue,rnd);
    
    //compute mag
    magnetization(magn,conf,quark,rnd,chi,point_magn,arg,mu,nu);
    
    //remove backfield
    rem_backfield_with_stagphases_from_conf(conf,u1b);
    
    //free
    for(int par=0;par<2;par++) nissa_free(chi[par]);
    nissa_free(point_magn);
    nissa_free(arg);
  }
  THREADABLE_FUNCTION_END
  
  //compute the magnetization
  void magnetization(complex *magn,complex *magn_free,rnd_t rnd_type,quad_su3 **conf,quad_su3 **free_conf,int quantization,quad_u1 **u1b,quark_content_t *quark,double residue)
  {
    //allocate source and generate it
    color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
    generate_fully_undiluted_eo_source(rnd,rnd_type,-1);
    
    //call inner function
    magnetization(magn,conf,quantization,u1b,quark,residue,rnd);
    magnetization(magn_free,free_conf,quantization,u1b,quark,residue,rnd);
    
    for(int par=0;par<2;par++) nissa_free(rnd[par]);
  }
  
  //measure magnetization
  void measure_magnetization(quad_su3 **conf,theory_pars_t &theory_pars,magnetization_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    FILE *file_free=open_file(meas_pars.path+"_free",conf_created?"w":"a");
    
    quad_su3 *free_conf[2];
    for(int eo=0;eo<2;eo++)
      free_conf[eo]=nissa_malloc("free_conf",loc_volh+bord_volh,quad_su3);
    generate_cold_eo_conf(free_conf);
    
    int ncopies=meas_pars.ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      {
        master_fprintf(file,"%d",iconf);
        master_fprintf(file_free,"%d",iconf);
        
        //measure magnetization for each quark
        for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
          {
	    if(theory_pars.quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	    
            complex magn={0,0};
            complex magn_free={0,0};
            
            //loop over hits
            int nhits=meas_pars.nhits;
            for(int hit=0;hit<nhits;hit++)
              {
                verbosity_lv2_master_printf("Evaluating magnetization for flavor %d/%d, ncopies %d/%d nhits %d/%d\n",
                                            iflav+1,theory_pars.nflavs(),icopy+1,ncopies,hit+1,nhits);
            
                //compute and summ
                complex temp,temp_free;
                magnetization(&temp,&temp_free,meas_pars.rnd_type,conf,free_conf,theory_pars.em_field_pars.flag,theory_pars.backfield[iflav],&theory_pars.quarks[iflav],meas_pars.residue); //flag holds quantization
                
                //normalize
                complex_summ_the_prod_double(magn,temp,1.0/nhits);
                complex_summ_the_prod_double(magn_free,temp_free,1.0/nhits);
              }
            
            //output
            master_fprintf(file,"\t%+016.16lg \t%+016.16lg",magn[RE],magn[IM]);
            master_fprintf(file_free,"\t%+016.16lg \t%+016.16lg",magn_free[RE],magn_free[IM]);
          }
        
        master_fprintf(file,"\n");
        master_fprintf(file_free,"\n");
      }
    
    for(int eo=0;eo<2;eo++)
      nissa_free(free_conf[eo]);
    
    close_file(file);
    close_file(file_free);
  }
  
  //print
  std::string magnetization_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasMagnetiz\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
