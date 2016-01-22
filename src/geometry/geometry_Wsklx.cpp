#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "new_types/su3.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //initialize the vector explaining how to store data when applying Wilson
  //hopping matrix on a sink-based way, without E/O preconditioning
  void define_Wsklx_hopping_matrix_output_pos()
  {
    two_stage_computation_pos_t *out=&Wsklx_hopping_matrix_output_pos; //for brevity
    out->inter_fr_in_pos=nissa_malloc("inter_fr_in_pos",2*NDIM*loc_vol,int); //offset for intermediate result
    out->final_fr_inter_pos=NULL; //we are not overlapping communication with intermediate->final,
    //so final_fr_inter_pos=inter, basically
    out->inter_fr_recv_pos=nissa_malloc("inter_fr_recv_pos",bord_vol,int); //offset for intermediate from nissa_recv_buf
    
    //separate local from ext border
    for(int iWsklx=0;iWsklx<loc_vol;iWsklx++)
      {
	int iloclx=loclx_of_Wsklx[iWsklx];
	
	for(int idir=0;idir<NDIM;idir++)
	  {
	    int bw=loclx_neighdw[iloclx][idir]; //we should subtract loc_vol from bw and put 2*NDIM-1 in place of 2*NDIM
	    out->inter_fr_in_pos[iWsklx*2*NDIM+idir]=(bw>=loc_vol)?bw+(2*NDIM-1)*loc_vol:(2*NDIM)*Wsklx_of_loclx[bw]+idir;
	    
	    int fw=loclx_neighup[iloclx][idir]; //idem for fw
	    out->inter_fr_in_pos[iWsklx*2*NDIM+NDIM+idir]=(fw>=loc_vol)?fw+(2*NDIM-1)*loc_vol:(2*NDIM)*Wsklx_of_loclx[fw]+NDIM+idir;
	  }
      }
    
    //final dest of incoming border
    for(int mu=0;mu<NDIM;mu++)
      for(int base_src=0;base_src<bord_dir_vol[mu];base_src++)
	{
	  int bw_src=bord_offset[mu]+base_src;
	  out->inter_fr_recv_pos[bw_src]=2*NDIM*Wsklx_of_loclx[surflx_of_bordlx[bw_src]]+NDIM+mu;
	  
	  int fw_src=bord_volh+bord_offset[mu]+base_src;
	  out->inter_fr_recv_pos[fw_src]=2*NDIM*Wsklx_of_loclx[surflx_of_bordlx[fw_src]]+mu;
	}
  }
  
  /*
  //Define output index of sink-applied hopping matrix for non-eo preco operator.
  void define_Wsklx_hopping_matrix_lx_output_pointers()
  {
  //allocate hopping matrix output pointer
  Wsklx_hopping_matrix_output_pointer=nissa_malloc("Wsklx_hopping_matrix_output_pointer",8*loc_vol,int);
  Wsklx_hopping_matrix_final_output=nissa_malloc("Wsklx_hopping_matrix_final_output",bord_vol,int);
  
  //separate local from ext border
  for(int iWsklx=0;iWsklx<loc_vol;iWsklx++)
  {
  int iloclx=loclx_of_Wsklx[iWsklx];
  
  for(int idir=0;idir<4;idir++)
  {
  int bw=loclx_neighdw[iloclx][idir];
  Wsklx_hopping_matrix_output_pointer[iWsklx*8+idir]=(bw>=loc_vol)?bw+7*loc_vol:8*Wsklx_of_loclx[bw]+idir;
  
  int fw=loclx_neighup[iloclx][idir];
  Wsklx_hopping_matrix_output_pointer[iWsklx*8+4+idir]=(fw>=loc_vol)?fw+7*loc_vol:8*Wsklx_of_loclx[fw]+4+idir;
  }
  }
  
  //final dest of incoming border
  for(int mu=0;mu<4;mu++)
  for(int base_src=0;base_src<bord_dir_vol[mu];base_src++)
  {
  int bw_src=bord_offset[mu]+base_src;
  Wsklx_hopping_matrix_final_output[bw_src]=8*Wsklx_of_loclx[surflx_of_bordlx[bw_src]]+4+mu;
  
  int fw_src=bord_volh+bord_offset[mu]+base_src;
  Wsklx_hopping_matrix_final_output[fw_src]=8*Wsklx_of_loclx[surflx_of_bordlx[fw_src]]+mu;
  }
  }
  */
  
  //Define Wsk ordering: surface sites comes first, then bulk sites
  void set_Wsklx_order()
  {
    Wsklx_order_inited=1;
    
    Wsklx_of_loclx=nissa_malloc("Wsklx_of_loclx",loc_vol,int);
    loclx_of_Wsklx=nissa_malloc("loclx_of_Wsklx",loc_vol,int);
    
    //reset Wsklx index
    int Wsklx=0;
    
    //surface
    for(int isurflx=0;isurflx<surf_vol;isurflx++)
      {
	int loclx=loclx_of_surflx[isurflx];
	
	loclx_of_Wsklx[Wsklx]=loclx;
	Wsklx_of_loclx[loclx]=Wsklx;
	
	Wsklx++;
      }
    
    //bulk
    for(int ibulklx=0;ibulklx<bulk_vol;ibulklx++)
      {
	int loclx=loclx_of_bulklx[ibulklx];
	
	loclx_of_Wsklx[Wsklx]=loclx;
	Wsklx_of_loclx[loclx]=Wsklx;
	
	Wsklx++;
      }
    
    if(Wsklx!=loc_vol) crash("defining Wsk_lx ordering: %d!=%d",Wsklx,loc_vol);
    
    define_Wsklx_hopping_matrix_output_pos();
  }
  
  //free
  void unset_Wsklx_order()
  {
    nissa_free(Wsklx_of_loclx);
    nissa_free(loclx_of_Wsklx);
    
    nissa_free(Wsklx_hopping_matrix_output_pos.inter_fr_in_pos);
    nissa_free(Wsklx_hopping_matrix_output_pos.inter_fr_recv_pos);
    
    //nissa_free(Wsklx_hopping_matrix_output_pointer);
    //nissa_free(Wsklx_hopping_matrix_final_output);
  }
  
  /*
    put together the 2*NDIM links to be applied to a single point
    first comes the links needed to scatter backward the signal (not to be daggered)
    then those needed to scatter it forward (to be daggered)
  */
  THREADABLE_FUNCTION_2ARG(lx_conf_remap_to_Wsklx, oct_su3*,out, quad_su3*,in)
  {
    GET_THREAD_ID();
    
    //communicate conf border so we can accede it
    communicate_lx_quad_su3_borders(in);
    
    NISSA_PARALLEL_LOOP(isrc_lx,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  //catch links needed to scatter signal forward
	  su3_copy(out[Wsklx_of_loclx[isrc_lx]][NDIM+mu],in[isrc_lx][mu]);
	  
	  //copy links also where they are needed to scatter the signal backward, if 
	  //sites that need them are not in the border (that would mean that computation must be 
	  //done in another node)
	  int idst_lx=loclx_neighup[isrc_lx][mu];
	  if(idst_lx<loc_vol) su3_copy(out[Wsklx_of_loclx[idst_lx]][mu],in[isrc_lx][mu]);
	}
    
    //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward 
    for(int mu=0;mu<NDIM;mu++) //border and link direction
      if(paral_dir[mu])
	NISSA_PARALLEL_LOOP(ibord,loc_vol+bord_offset[mu],loc_vol+bord_offset[mu]+bord_dir_vol[mu])
	  su3_copy(out[Wsklx_of_loclx[loclx_neighup[ibord][mu]]][mu],in[ibord][mu]);
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END

#define DEFINE_lx_remap_to_Wsklx_or_reverse(T)				\
  /*remap a type T from lx to Wsklx layout or reverse*/			\
  THREADABLE_FUNCTION_3ARG(lx_remap_to_Wsklx_or_reverse, T*,ext_out, T*,in, int,lx_to_Wsklx) \
  {									\
    GET_THREAD_ID();							\
									\
    /*bufferize if needed*/						\
    int bufferize=(void*)ext_out==(void*)in;				\
    T *out=bufferize?nissa_malloc("out",loc_vol,T):ext_out;		\
									\
    if(lx_to_Wsklx) NISSA_PARALLEL_LOOP(ivol_lx,0,loc_vol) memcpy(out[Wsklx_of_loclx[ivol_lx]],in[ivol_lx],sizeof(T)); \
    else            NISSA_PARALLEL_LOOP(iWsk_lx,0,loc_vol) memcpy(out[loclx_of_Wsklx[iWsk_lx]],in[iWsk_lx],sizeof(T)); \
									\
    /*wait filling*/							\
    set_borders_invalid(out);						\
									\
    /*unbuffer if needed*/						\
    if(bufferize)							\
      {									\
	vector_copy(ext_out,out);					\
	nissa_free(out);						\
      }									\
  }									\
  THREADABLE_FUNCTION_END						\
									\
    void lx_remap_to_Wsklx(T *out,T *in)				\
    {lx_remap_to_Wsklx_or_reverse(out,in,1);}				\
    void Wsklx_remap_to_lx(T *out,T *in)				\
    {lx_remap_to_Wsklx_or_reverse(out,in,0);}

DEFINE_lx_remap_to_Wsklx_or_reverse(spincolor)
DEFINE_lx_remap_to_Wsklx_or_reverse(colorspinspin)
DEFINE_lx_remap_to_Wsklx_or_reverse(su3spinspin)
}
