#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "communicate/communicate.hpp"
#include "new_types/su3.hpp"

#include "stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif


namespace nissa
{
  //check it is an hypercube origin
  inline int is_hypercube_origin(int ivol)
  {
    int is=1;
    for(int mu=0;mu<NDIM;mu++) is&=(loc_coord_of_loclx[ivol][mu]%2==0);
    return is;
  }
  
  void prepare_source(int idest,quad_su3 **conf,int nto_av,int *dir_to_shift,int *map_sources,color *ori[8][2])
  {
    GET_THREAD_ID();
    
    color **dest=ori[idest];
    //communicate all sources
    for(int isource=0;isource<nto_av;isource++)
      communicate_ev_and_od_color_borders(ori[isource]);
    
    for(int eo=0;eo<2;eo++)
      {
	vector_reset(dest[eo]);
	for(int isource=0;isource<nto_av;isource++)
	  {
	    int mu=dir_to_shift[isource];
	    NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	      {
		int iave=loceo_neighdw[eo][ieo][mu];
		su3_dag_summ_the_prod_color(dest[eo][ieo],conf[!eo][iave][mu],ori[map_sources[isource]][!eo][iave]);
	      }
	  }
	set_borders_invalid(dest[eo]);
      }
  }
  
  THREADABLE_FUNCTION_4ARG(staggered_meson_corr, quad_su3**,conf, theory_pars_t*,tp, int,iconf, int,conf_created)
  {
    GET_THREAD_ID();
    
    int nflavs=tp->nflavs;
    int ncombo=nflavs*(nflavs+1)/2;
    int nhits=1;
    
    //allocate
    color *source[8][2],*sol[nflavs][8][2];
    for(int isite=0;isite<8;isite++)
      for(int eo=0;eo<2;eo++)
	{
	  source[isite][eo]=nissa_malloc("source",loc_volh+bord_volh,color);
	  for(int iflav=0;iflav<nflavs;isite++) sol[iflav][isite][eo]=nissa_malloc("sol",loc_volh+bord_volh,color);
	}
    complex *glb_corr=nissa_malloc("glb_corr",glb_size[0]*ncombo,complex);
    complex *thread_corr=new complex[glb_size[0]*ncombo];
    memset(thread_corr,0,sizeof(complex)*glb_size[0]*ncombo);
    
    for(int ihit=0;ihit<nhits;ihit++)
      {
	//prepare the source
	for(int eo=0;eo<2;eo++) vector_reset(source[0][eo]);
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  {
	    int ivol=loclx_of_loceo[EVN][ieo];
	    
	    //fill only if it is an hypercube origin
	    if(is_hypercube_origin(ivol)&&glb_coord_of_loclx[ivol][0]==0)
	      for(int icol=0;icol<NCOL;icol++)
		rnd_get_Z4(source[0][EVN][ieo][icol],loc_rnd_gen+ivol);
	  }
	THREAD_BARRIER();
	
	//perform all shifts
	///////////////////////////// 0  X   Y   XY    Z   XZ    YZ    XYZ
	int nsource_to_average[8]=   {0, 1,  1,  2,    1,  2,    2,    3};
	int sources_to_average[8][3]={{},{0},{0},{1,2},{0},{1,4},{2,4},{3,5,6}};
	int dir_to_shift      [8][3]={{},{1},{2},{2,1},{3},{3,1},{3,2},{3,2,1}};
	for(int isource=1;isource<8;isource++)
	  prepare_source(isource,conf,nsource_to_average[isource],dir_to_shift[isource],sources_to_average[isource],source);
	
	//compute all props
	for(int iflav=0;iflav<nflavs;iflav++)
	  for(int isource=0;isource<8;isource++)
	    mult_Minv(sol[iflav][isource],conf,tp,iflav,1e-12,source[isource]);
      }
    
    //free
    for(int isite=0;isite<8;isite++)
      for(int eo=0;eo<2;eo++)
	{
	  nissa_free(source[isite][eo]);
	  for(int iflav=0;iflav<nflavs;iflav++) nissa_free(sol[iflav][isite][eo]);
	}
    delete[] thread_corr;
    nissa_free(glb_corr);
  }
  THREADABLE_FUNCTION_END
}
