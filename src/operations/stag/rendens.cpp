#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "io/input.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"

#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "rendens.hpp"
#include "stag.hpp"

namespace nissa
{
  namespace rende
  {
    //multpiply by the derivative of M w.r.t mu
    THREADABLE_FUNCTION_6ARG(mult_dMdmu, color**,out, theory_pars_t*,pars, quad_su3**,conf, int,iflav, int,ord, color**,in)
    {
      GET_THREAD_ID();
      
      if(ord==0) crash("makes no sense to call with order zero");
      
      add_backfield_to_conf(conf,pars->backfield[iflav]);
      
      for(int par=0;par<2;par++)
	{
	  communicate_ev_or_od_color_borders(in[!par],!par);
	  
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      color temp;
	      unsafe_su3_prod_color(temp,conf[par][ieo][0],in[!par][loceo_neighup[par][ieo][0]]);
	      int idw=loceo_neighdw[par][ieo][0];
	      if(ord%2==0) su3_dag_subt_the_prod_color(temp,conf[!par][idw][0],in[!par][idw]);
	      else         su3_dag_summ_the_prod_color(temp,conf[!par][idw][0],in[!par][idw]);
	      color_prod_double(out[par][ieo],temp,0.5);
	      set_borders_invalid(out[par]);
	    }
	}
      
      rem_backfield_from_conf(conf,pars->backfield[iflav]);
    }
    THREADABLE_FUNCTION_END
	
    //take the trace between A^dag and B
    THREADABLE_FUNCTION_4ARG(summ_the_trace, double*,out, complex*,point_result, color**, A, color**, B)
    {
      GET_THREAD_ID();
      
      //compute results for single points
      vector_reset(point_result);
      for(int par=0;par<2;par++)
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  for(int ic=0;ic<3;ic++)
	    complex_summ_the_conj1_prod(point_result[loclx_of_loceo[par][ieo]],A[par][ieo][ic],B[par][ieo][ic]);
      THREAD_BARRIER();
      
      //final reduction
      complex temp;
      complex_vector_glb_collapse(temp,point_result,loc_vol);
      complex_summassign(out,temp);
    }
    THREADABLE_FUNCTION_END
    
    //fill a source
    void fill_rende_source(color **src)
    {generate_fully_undiluted_eo_source(src,RND_GAUSS,-1);}
    
    typedef color* rende_t[2];
#define NEW_RENDE_T(A)					\
    rende_t A;						\
    A[0]=nissa_malloc(#A,loc_volh+bord_volh,color);	\
    A[1]=nissa_malloc(#A,loc_volh+bord_volh,color)
#define DELETE_RENDE_T(A)				\
    nissa_free(A[0]);					\
    nissa_free(A[1]);
#define MINV(out,iflav,in)						\
    mult_Minv(out,conf,&pars,iflav,pars.quark_rendens_meas_pars.residue,in,true)
#define NEW_MINV(out,iflav,in)			\
    NEW_RENDE_T(out);				\
    MINV(out,iflav,in)
#define DM(out,iflav,ord,in)			\
    mult_dMdmu(out,&pars,conf,iflav,ord,in)
#define NEW_DM(out,iflav,ord,in)		\
    NEW_RENDE_T(out);				\
    DM(out,iflav,ord,in)
#define NEW_TRACE_RES(o)			\
    complex o={0,0}
#define SUMM_THE_TRACE(A,B,C)				\
      summ_the_trace((double*)A,point_result,B,C)
#define PRINT(A)							\
      master_fprintf(file,"%+016.016lg %+016.016lg\t",A[0]/pars.quark_rendens_meas_pars.nhits,A[1]/pars.quark_rendens_meas_pars.nhits)
      }
  
  using namespace rende;
  
  //measure the quark number and its derivative w.r.t mu
  void measure_quark_rendens(quad_su3 **conf,theory_pars_t &pars,int iconf,int conf_created)
  {
    //open the file, allocate point result and source
    FILE *file=open_file(pars.quark_rendens_meas_pars.path,conf_created?"w":"a");
    complex *point_result=nissa_malloc("point_result",loc_vol,complex);
    NEW_RENDE_T(source);
    addrem_stagphases_to_eo_conf(conf);
    
    //vectors for calculation
    NEW_RENDE_T(dM);
    NEW_RENDE_T(d2M);
    NEW_RENDE_T(M_dM);
    NEW_RENDE_T(M_d2M);
    
    for(int icopy=0;icopy<pars.quark_rendens_meas_pars.ncopies;icopy++)
      {
	//print conf id
	master_fprintf(file,"%d\t",iconf);
	
	//loop over flavor
	for(int iflav=0;iflav<pars.nflavs;iflav++)
	  {
	    //vectors for output
	    NEW_TRACE_RES(Tr_M_dM);
	    NEW_TRACE_RES(Tr_M_d2M);
	    
	    //loop over hits
	    for(int ihit=0;ihit<pars.quark_rendens_meas_pars.nhits;ihit++)
	      {
		//fill the source
		fill_rende_source(source);
		
		//compute M^-1*dM
		DM(dM,iflav,1,source);
		MINV(M_dM,iflav,dM);
		
		//compute M^-1*d2M
		DM(d2M,iflav,2,source);
		MINV(M_d2M,iflav,d2M);
		
		//trace
		SUMM_THE_TRACE(Tr_M_dM,source,M_dM);
		SUMM_THE_TRACE(Tr_M_d2M,source,M_d2M);
	      }
	    
	    //print out (automatical normalisation for nhits)
	    PRINT(Tr_M_dM);
	    PRINT(Tr_M_d2M);
	  }
	
	master_fprintf(file,"\n");
      }
    
    DELETE_RENDE_T(dM);
    DELETE_RENDE_T(d2M);
    DELETE_RENDE_T(M_dM);
    DELETE_RENDE_T(M_d2M);
	
    //close and deallocate
    addrem_stagphases_to_eo_conf(conf);
    close_file(file);
    nissa_free(point_result);
    DELETE_RENDE_T(source);
  }
}
