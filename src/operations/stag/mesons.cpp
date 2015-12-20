#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_mix.hpp"
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
  
  //perform a shift inside the hypercube
  void shift_in_hypercube(color *dest,quad_su3 *conf,int nto_av,int *dir_to_shift,int *map_sources,color **ori)
  {
    GET_THREAD_ID();
    
    //communicate all sources and the conf
    communicate_lx_quad_su3_borders(conf);
    
    vector_reset(dest);
    for(int isource=0;isource<nto_av;isource++)
      {
	int mu=dir_to_shift[isource];
	int iori=map_sources[isource];
	communicate_lx_color_borders(ori[iori]);
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  if(glb_coord_of_loclx[ivol][mu]%2==1)
	    {
	      int idw=loclx_neighdw[ivol][mu];
	      su3_dag_summ_the_prod_color(dest[ivol],conf[idw][mu],ori[iori][idw]);
	    }
	  else
	    {
	      int iup=loclx_neighup[ivol][mu];
	      su3_summ_the_prod_color(dest[ivol],conf[ivol][mu],ori[iori][iup]);
	    }
      }
    THREAD_BARRIER();
    
    //final normalization
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      color_prod_double(dest[ivol],dest[ivol],1.0/nto_av);
    set_borders_invalid(dest);
  }
  
  //compute correlation functions for staggered mesons, arbitary taste and spin
  THREADABLE_FUNCTION_4ARG(staggered_meson_corr, quad_su3**,conf, theory_pars_t*,tp, int,iconf, int,conf_created)
  {
    GET_THREAD_ID();
    
    std::vector<std::pair<int,int> > mesons;
    int nop=mesons.size();
    
    int nflavs=tp->nflavs;
    int ncombo=nop*nflavs*(nflavs+1)/2;
    int nhits=1;
    
    //allocate
    const int nshift=8;
    color *source[nshift],*sol[nflavs][nshift*nshift];
    for(int ishift=0;ishift<nshift;ishift++)
      source[ishift]=nissa_malloc("source",loc_vol+bord_vol,color);
    for(int iflav=0;iflav<nflavs;iflav++)
      for(int ijshift=0;ijshift<nshift*nshift;ijshift++)
	    sol[iflav][ijshift]=nissa_malloc("sol",loc_vol+bord_vol,color);
    complex *glb_corr=nissa_malloc("glb_corr",glb_size[0]*ncombo,complex);
    complex *thread_corr=new complex[glb_size[0]*ncombo];
    quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol,quad_su3);
    memset(thread_corr,0,sizeof(complex)*glb_size[0]*ncombo);
    color *temp_source[2],*temp_sol[2];
    for(int eo=0;eo<2;eo++)
      {
	temp_source[eo]=nissa_malloc("temp_source",loc_volh+bord_volh,color);
	temp_sol[eo]=nissa_malloc("temp_sol",loc_volh+bord_volh,color);
      }
    
    //convert the conf and smear
    paste_eo_parts_into_lx_vector(lx_conf,conf);
    
    for(int ihit=0;ihit<nhits;ihit++)
      {
	//prepare the source
	//filling only if it is an hypercube origin
	vector_reset(source[0]);
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  if(is_hypercube_origin(ivol)&&glb_coord_of_loclx[ivol][0]==0)
	    for(int icol=0;icol<NCOL;icol++)
	      rnd_get_Z4(source[0][ivol][icol],loc_rnd_gen+ivol);
      	THREAD_BARRIER();
	
	//perform all shifts
	///////////////////////////// 0  X   Y   XY    Z   XZ    YZ    XYZ
	int nsource_to_ave[nshift]=   {0, 1,  1,  2,    1,  2,    2,    3};
	int sources_to_ave[nshift][3]={{},{0},{0},{1,2},{0},{1,4},{2,4},{3,5,6}};
	int dir_to_shift      [nshift][3]={{},{1},{2},{2,1},{3},{3,1},{3,2},{3,2,1}};
	for(int ishift=1;ishift<nshift;ishift++)
	  shift_in_hypercube(source[ishift],lx_conf,nsource_to_ave[ishift],dir_to_shift[ishift],sources_to_ave[ishift],source);
	
	//compute all props
	for(int iflav=0;iflav<nflavs;iflav++)
	  for(int so_shift=0;so_shift<nshift;so_shift++)
	    {
	      int ba_shift=so_shift*nshift;
	      
	      //convert
	      split_lx_vector_into_eo_parts(temp_source,source[so_shift]);
	      mult_Minv(temp_sol,conf,tp,iflav,1e-12,temp_source);
	      paste_eo_parts_into_lx_vector(sol[iflav][0+ba_shift],temp_sol);
	      
	      //shift the sink
	      for(int si_shift=1;si_shift<nshift;si_shift++)
		shift_in_hypercube(sol[iflav][si_shift+ba_shift],lx_conf,nsource_to_ave[si_shift],dir_to_shift[si_shift],sources_to_ave[si_shift],sol[iflav]+ba_shift);
	    }
	

	//int itaste=5,ispin=5;
	//for(int ia=0;ia<8;ia++)
	//{
	    //int ib=ia^(itaste^ispin);
	//}
      }
    
    for(int ishift=0;ishift<nshift;ishift++) nissa_free(source[ishift]);
    for(int iflav=0;iflav<nflavs;iflav++) for(int ijshift=0;ijshift<nshift*nshift;ijshift++) nissa_free(sol[iflav][ijshift]);
    for(int eo=0;eo<2;eo++)
      {
	nissa_free(temp_source[eo]);
	nissa_free(temp_sol[eo]);
      }
    delete[] thread_corr;
    nissa_free(glb_corr);
    nissa_free(lx_conf);
  }
  THREADABLE_FUNCTION_END
}
