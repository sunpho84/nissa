#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"
#include "operations/gauge_fixing.hpp"

#include "stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

//in this formalism, shift=t+2*(x+2*(y+2*z))

namespace nissa
{
  namespace
  {
    const int nshift=8;
    
    //return final combo
    inline int ind_combo(int iop,int iflav,int jflav,int nflavs)
    {return jflav+nflavs*(iflav+nflavs*iop);}
    
    //return the index of the sink-summed
    inline int ind_summ(int iflav,int sto_so_shift,int sto_si_shift,int t,int so_col,int si_col)
    {return si_col+NCOL*(so_col+NCOL*(t+glb_size[0]*(sto_si_shift+nshift*(sto_so_shift+nshift*iflav))));}
  }
  
  //check it is an hypercube origin
  inline int is_hypercube_shift(int ivol,int ishift)
  {
    int is=1;
    for(int mu=0;mu<NDIM;mu++) is&=((glb_coord_of_loclx[ivol][mu]%2)==((ishift>>mu)==1));
    return is;
  }
  
  //form the mask for x (-1)^[x*(s^<+n^>)](-1)^[n*(s+n)^<]
  inline int form_stag_meson_pattern(int ispin,int itaste)
  {
    //add g5*g5
    ispin^=15;
    itaste^=15;
    
    int res=0;
    for(int mu=NDIM-1;mu>=0;mu--)
      {
	int p=0;
	for(int nu=0;nu<mu;nu++) p+=(itaste>>nu)&1;
	for(int nu=mu+1;nu<NDIM;nu++) p+=(ispin>>nu)&1;
	p&=1;
	
	res=res*2+p;
      }
    
    return res;
  }
  
  //prepare the source
  //filling only if it is an hypercube origin
  void prepare_source_in_hypercube_origin(color *source[2],int ishift,int icol,int t)
  {
    GET_THREAD_ID();
    
    for(int eo=0;eo<2;eo++)
      {
	vector_reset(source[eo]);
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  {
	    int ivol=loclx_of_loceo[eo][ieo];
	    if(is_hypercube_shift(ivol,ishift)&&glb_coord_of_loclx[ivol][0]==t)
	      complex_put_to_real(source[eo][ieo][icol],1);
	  }
	set_borders_invalid(source[eo]);
      }
  }
  
  // //perform a shift inside the hypercube
  // void shift_in_hypercube(color *dest,quad_su3 *conf,int nto_av,int *dir_to_shift,int *map_sources,color **ori)
  // {
  //   GET_THREAD_ID();
  //
  //   //communicate all sources and the conf
  //   communicate_lx_quad_su3_borders(conf);
  //
  //   vector_reset(dest);
  //   for(int isource=0;isource<nto_av;isource++)
  //     {
  // 	int mu=dir_to_shift[isource];
  // 	int iori=map_sources[isource];
  // 	communicate_lx_color_borders(ori[iori]);
  //
  // 	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
  // 	  if(glb_coord_of_loclx[ivol][mu]%2==1)
  // 	    {
  // 	      int idw=loclx_neighdw[ivol][mu];
  // 	      su3_dag_summ_the_prod_color(dest[ivol],conf[idw][mu],ori[iori][idw]);
  // 	    }
  // 	  else
  // 	    {
  // 	      int iup=loclx_neighup[ivol][mu];
  // 	      su3_summ_the_prod_color(dest[ivol],conf[ivol][mu],ori[iori][iup]);
  // 	    }
  //     }
  //   THREAD_BARRIER();
  //
  //   //final normalization
  //   NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
  //     color_prod_double(dest[ivol],dest[ivol],1.0/nto_av);
  //   set_borders_invalid(dest);
  // }
  // //perform all shifts
  // ///////////////////////////// 0  X   Y   XY    Z   XZ    YZ    XYZ
  // int nsource_to_ave[nshift]=   {0, 1,  1,  2,    1,  2,    2,    3};
  // int sources_to_ave[nshift][3]={{},{0},{0},{1,2},{0},{1,4},{2,4},{3,5,6}};
  // int dir_to_shift      [nshift][3]={{},{1},{2},{2,1},{3},{3,1},{3,2},{3,2,1}};
  // for(int ishift=1;ishift<nshift;ishift++)
  //   shift_in_hypercube(source[ishift],lx_conf,nsource_to_ave[ishift],dir_to_shift[ishift],sources_to_ave[ishift],source);
  
  
  //summ on the sink with a given shift
  inline void store_sol(complex *sink_summed,int iflav,int so_col,int sto_so_shift,color **sol)
  {
    GET_THREAD_ID();
    
    for(int eo=0;eo<2;eo++)
      NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	{
	  //find lx and time
	  int ivol=loclx_of_loceo[eo][ieo];
	  int t=glb_coord_of_loclx[ivol][0];
	  
	  //compute the shift index of the sink
	  int si_shift=0;
	  for(int mu=NDIM-1;mu>=0;mu--) si_shift=2*si_shift+glb_coord_of_loclx[ivol][mu]%2;
	  int sto_si_shift=si_shift/2;
	  
	  for(int si_col=0;si_col<NCOL;si_col++)
	    complex_summassign(sink_summed[ind_summ(iflav,sto_so_shift,sto_si_shift,t,so_col,si_col)],sol[eo][ieo][si_col]);
	}
    THREAD_BARRIER();
  }
  
  //return the other end of the pair
  inline int get_other_sto_shift(int sto_shift,int mask)
  {return ((sto_shift*2)^mask)/2;}
  
  //compute correlation functions for staggered mesons, arbitary taste and spin
  THREADABLE_FUNCTION_4ARG(staggered_meson_corr, quad_su3**,ext_conf, theory_pars_t*,tp, int,iconf, int,conf_created)
  {
    GET_THREAD_ID();
    
    //list of operators
    std::vector<std::pair<int,int> > mesons;
    mesons.push_back(std::make_pair(15,15));
    
    int nop=mesons.size();
    int nflavs=tp->nflavs;
    int ncombo=ind_combo(nop-1,nflavs-1,nflavs-1, nflavs)+1;
    int nhits=1;
    
    //allocate
    color *source[2],*sol[2];
    for(int eo=0;eo<2;eo++)
      {
	sol[eo]=nissa_malloc("sol",loc_volh+bord_volh,color);
	source[eo]=nissa_malloc("source",loc_volh+bord_volh,color);
      }
    complex *glb_corr=nissa_malloc("glb_corr",glb_size[0]*ncombo,complex);
    quad_su3 *gf_conf[2];
    for(int eo=0;eo<2;eo++) gf_conf[eo]=nissa_malloc("gf_conf",loc_volh+bord_volh,quad_su3);
    int sink_summed_size=ind_summ(nflavs-1,nshift-1,nshift-1,glb_size[0]-1,NCOL-1,NCOL-1)+1;
    complex *sink_summed=new complex[sink_summed_size];
    for(int i=0;i<sink_summed_size;i++) complex_put_to_zero(sink_summed[i]);
    
    //perform gauge fixing
    quad_su3 *temp_conf=nissa_malloc("temp_conf",loc_vol+bord_vol,quad_su3);
    paste_eo_parts_into_lx_vector(temp_conf,ext_conf);
    coulomb_gauge_fix(temp_conf,temp_conf,1e-14);
    split_lx_vector_into_eo_parts(gf_conf,temp_conf);
    nissa_free(temp_conf);
    
    //form the masks
    int mask[nop];
    for(int iop=0;iop<nop;iop++) mask[iop]=form_stag_meson_pattern(mesons[iop].first,mesons[iop].second);
    
    for(int ihit=0;ihit<nhits;ihit++)
      {
	//compute all props and summ over sink
	for(int so_col=0;so_col<NCOL;so_col++)
	  for(int sto_so_shift=0;sto_so_shift<nshift;sto_so_shift++)
	    {
	      int so_shift=sto_so_shift*2;
	      prepare_source_in_hypercube_origin(source,so_shift,so_col,ihit*glb_size[0]/nhits);
	      for(int iflav=0;iflav<nflavs;iflav++)
		{
		  mult_Minv(sol,gf_conf,tp,iflav,1e-12,source);
		  store_sol(sink_summed,iflav,so_col,sto_so_shift,sol);
		}
	    }
	
	//reduce
	glb_threads_reduce_double_vect((double*)sink_summed,2*sink_summed_size);
	if(IS_MASTER_THREAD) glb_nodes_reduce_complex_vect(sink_summed,sink_summed_size);
	THREAD_BARRIER();
	
	//bind them
	if(IS_MASTER_THREAD)
	  for(int iop=0;iop<nop;iop++)
	    {
	      int ispin=mesons[iop].first;
	      int itaste=mesons[iop].second;
	      
	      for(int iflav=0;iflav<nflavs;iflav++)
		for(int jflav=0;jflav<nflavs;jflav++)
		  {
		    int ic=ind_combo(iop,iflav,jflav, nflavs);
		    for(int sto_so_ishift=0;sto_so_ishift<nshift;sto_so_ishift++)
		      for(int sto_si_ishift=0;sto_si_ishift<nshift;sto_si_ishift++)
			{
			  int sto_so_jshift=get_other_sto_shift(sto_so_ishift,ispin^itaste);
			  int sto_si_jshift=get_other_sto_shift(sto_si_ishift,ispin^itaste);

			  int expo=((sto_so_ishift*2)&mask[iop])^((sto_si_ishift*2)&mask[iop]);
			  int sign=1;
			  for(int temp=expo;temp;temp/=2) sign*=1-2*(temp&1);
			  
			  for(int t=0;t<glb_size[0];t++)
			    for(int so_col=0;so_col<NCOL;so_col++)
			      for(int si_col=0;si_col<NCOL;si_col++)
				{
				  complex temp;
				  unsafe_complex_conj2_prod(temp,
							    sink_summed[ind_summ(iflav,sto_so_ishift,sto_si_ishift,t,so_col,si_col)],
							    sink_summed[ind_summ(iflav,sto_so_jshift,sto_si_jshift,t,so_col,si_col)]);

				  complex_summ_the_prod_double(glb_corr[ic],temp,sign);
				}
			}
		  }
	      }
      }
    
    for(int eo=0;eo<2;eo++) nissa_free(sol[eo]);
    for(int eo=0;eo<2;eo++) nissa_free(source[eo]);
    for(int eo=0;eo<2;eo++) nissa_free(gf_conf[eo]);
    nissa_free(glb_corr);
    delete[] sink_summed;
  }
  THREADABLE_FUNCTION_END
}
