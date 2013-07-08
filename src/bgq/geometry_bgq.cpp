#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "bgq_macros.h"

#include "../communicate/communicate.h"
#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/thread_macros.h"
#include "../base/vectors.h"
#include "../new_types/complex.h"
#ifdef USE_THREADS
 #include "../routines/thread.h"
#endif

/*
  Define bgq ordering: surface sites comes first, then bulk sites
  Each site at coord c is matched with site at coord c+T/2
  Since t is the slowest running coord, the two sites pair is given
  always by ivol and ivol+loc_volh, so mapping array are defined 
  only for first half of the lattice
  
  We first scan T bw and fw surface, other surfaces (if present), then bulk.
*/
void define_bgq_lx_ordering()
{
  bgqlx_of_loclx=nissa_malloc("bgqlx_of_loclx",loc_volh,int);
  loclx_of_bgqlx=nissa_malloc("loclx_of_bgqlx",loc_volh,int);
  
  //reset bgqlx index
  int bgqlx=0;
  
  //scan T surface
  for(int loclx=0;loclx<loc_volh;loclx++)
    //check that we are on T==0 or T==loc_size[0]/2-1
    if(loc_coord_of_loclx[loclx][0]==0||loc_coord_of_loclx[loclx][0]==loc_size[0]/2-1)
      {
	loclx_of_bgqlx[bgqlx]=loclx;
	bgqlx_of_loclx[loclx]=bgqlx;
	
	bgqlx++;
      }
  
  //scan non-T surface
  for(int isurflx=0;isurflx<surf_vol;isurflx++)
    {
      int loclx=loclx_of_surflx[isurflx];
      
      //check that we are on first half of time direction, and that we are not on T==0 or T==loc_size[0]/2-1
      if(loc_coord_of_loclx[loclx][0]!=0 && loc_coord_of_loclx[loclx][0]<loc_size[0]/2-1)
	{
	  loclx_of_bgqlx[bgqlx]=loclx;
	  bgqlx_of_loclx[loclx]=bgqlx;
	  
	  bgqlx++;
	}
    }

  //take not of the total virtual and non-virtual surface
  bgq_vsurf_vol=bgqlx;
  
  //scan bulk
  for(int ibulklx=0;ibulklx<bulk_vol;ibulklx++)
    {
      int loclx=loclx_of_bulklx[ibulklx];
      
      //check that we are on first half of time direction, and that we are not on T==0 or T==loc_size[0]/2-1
      if(loc_coord_of_loclx[loclx][0]!=0 && loc_coord_of_loclx[loclx][0]<loc_size[0]/2-1)
	{
	  loclx_of_bgqlx[bgqlx]=loclx;
	  bgqlx_of_loclx[loclx]=bgqlx;
	  
	  bgqlx++;
	}
    }
  
  if(bgqlx!=loc_volh) crash("defining bgq_lx ordering: %d!=%d",bgqlx,loc_volh);
}

/*
  Define output index of sink-applied hopping matrix for non-eo preco operator.
  For T border we must store two different indices, one pointing to the communication buffer and the other to
  other virtual node buffer. In T buffer data is ordered as for buffered communication send buffer, that is,
  in the same order as explained in "communicate.h"
  Other three borders points directly to communication borders.
*/
void define_bgq_hopping_matrix_lx_output_pointers_and_T_buffers()
{
 //t dir is at least virtually parallelized
  bgqlx_t_vbord_vol=2*loc_vol/loc_size[0];
  
  //allocate hopping matrix output pointer
  bgq_hopping_matrix_output_pointer=nissa_malloc("bgq_hopping_matrix_output_pointer",8*loc_volh,bi_halfspincolor*);
  bgq_hopping_matrix_output_T_buffer=nissa_malloc("bgq_hopping_matrix_output_T_buff",bgqlx_t_vbord_vol,bi_halfspincolor);
  //T dir contains one pointer per true lattice site, others only half
  bgq_hopping_matrix_final_output=nissa_malloc("bgq_hopping_matrix_final_output",2*(bord_vol/2+bord_dir_vol[0]),bi_halfspincolor*);

  for(int ibgqlx=0;ibgqlx<loc_volh;ibgqlx++)
    {
      int iloclx=loclx_of_bgqlx[ibgqlx];
      
      //t direction bw scattering
      bgq_hopping_matrix_output_pointer[ibgqlx*8+0]=(loc_coord_of_loclx[iloclx][0]==0)?
	//we move in other vnode
	bgq_hopping_matrix_output_T_buffer+iloclx:
	//we are still in the same vnode
	bgq_hopping_matrix_output_data+8*bgqlx_of_loclx[loclx_neighdw[iloclx][0]]+0;
      
      //t direction fw scattering
      bgq_hopping_matrix_output_pointer[ibgqlx*8+4]=(loc_coord_of_loclx[iloclx][0]==loc_size[0]/2-1)?
	//we moved in another vnode
	bgq_hopping_matrix_output_T_buffer+bgqlx_t_vbord_vol/2+loclx_neighup[iloclx][0]-loc_volh:
	//we are still in the same vnode
	bgq_hopping_matrix_output_data+8*bgqlx_of_loclx[loclx_neighup[iloclx][0]]+4;
      
      //other direction derivatives
      for(int idir=1;idir<4;idir++)
	{
	  int bw=loclx_neighdw[iloclx][idir];
	  bgq_hopping_matrix_output_pointer[ibgqlx*8+0+idir]=(bw>=loc_vol)?
	    //we moved to another vnode
	    (bi_halfspincolor*)nissa_send_buf+bw-loc_vol-bord_offset[idir]/2:
	    //we are still in local vnode
	    bgq_hopping_matrix_output_data+8*bgqlx_of_loclx[bw]+0+idir;
	  
	  int fw=loclx_neighup[iloclx][idir];
	  bgq_hopping_matrix_output_pointer[ibgqlx*8+4+idir]=(fw>=loc_vol)?
	    //we moved in another vnode
	    (bi_halfspincolor*)nissa_send_buf+fw-loc_vol-bord_volh/2-bord_offset[idir]/2:              
	    //we are still in local vnode
	    bgq_hopping_matrix_output_data+8*bgqlx_of_loclx[fw]+4+idir;
	}
    }
  
  for(int base_src=0;base_src<bord_dir_vol[0];base_src++)
    {
      //T bw border (bw derivative): contributions start from 4
      int bw_src=base_src;
      int bw_dst=8*bgqlx_of_loclx[base_src]+4;
      bgq_hopping_matrix_final_output[bw_src]=bgq_hopping_matrix_output_data+bw_dst;
      
      //T fw border (fw derivative): contributions start from 0
      int fw_src=bord_vol/4+bord_dir_vol[0]/2+base_src; 
      int fw_dst=8*bgqlx_of_loclx[loclx_neighdw[base_src+loc_volh][0]]; 
      bgq_hopping_matrix_final_output[fw_src]=bgq_hopping_matrix_output_data+fw_dst;
    }
      
  for(int mu=1;mu<4;mu++)
    for(int base_src=0;base_src<bord_dir_vol[mu]/2;base_src++)
      {
	//other 3 bw borders
	int bw_src=bord_dir_vol[0]/2+bord_offset[mu]/2+base_src;
	int bw_dst=8*bgqlx_of_loclx[surflx_of_bordlx[bord_offset[mu]+base_src]]+4+mu;
	bgq_hopping_matrix_final_output[bw_src]=bgq_hopping_matrix_output_data+bw_dst;
	
	//other 3 fw borders
	int fw_src=bord_vol/4+bord_dir_vol[0]+bord_offset[mu]/2+base_src;
	int fw_dst=8*bgqlx_of_loclx[surflx_of_bordlx[bord_volh+bord_offset[mu]+base_src]]+mu;
	bgq_hopping_matrix_final_output[fw_src]=bgq_hopping_matrix_output_data+fw_dst;
      }
}

/*
  put together the 8 links to be applied to a single point
  first comes the links needed to scatter backward the signal (not to be daggered)
  then those needed to scatter it forward (to be daggered)
*/
THREADABLE_FUNCTION_2ARG(lx_conf_remap_to_bgqlx, bi_oct_su3*,out, quad_su3*,in)
{
  GET_THREAD_ID();
  
  //communicate conf border so we can accede it
  communicate_lx_quad_su3_borders(in);
  
  //scan the two virtual nodes
  for(int vn_fw_dst_bgqlx=0;vn_fw_dst_bgqlx<2;vn_fw_dst_bgqlx++)
    NISSA_PARALLEL_LOOP(ibase,0,loc_volh)
      {
	int isrc_lx=ibase+vn_fw_dst_bgqlx*loc_volh;
	int ifw_dst_bgqlx=bgqlx_of_loclx[ibase];
	
	for(int mu=0;mu<4;mu++)
	  {
	    //catch links needed to scatter signal forward
	    SU3_TO_BI_SU3(out[ifw_dst_bgqlx][4+mu],in[isrc_lx][mu],vn_fw_dst_bgqlx);

	    //copy links also where they are needed to scatter the signal backward, if 
	    //sites that need them are not in the border (that would mean that computation must be 
	    //done in another node
	    int idst_lx=loclx_neighup[isrc_lx][mu];
	    //1 only if idst_lx is in the first node, 2 if in boder
	    int vn_bw_dst_bgqlx=(idst_lx/loc_volh);
	    if(vn_bw_dst_bgqlx<2)
	      {
		int idst_bgqlx=bgqlx_of_loclx[idst_lx%loc_volh];
		SU3_TO_BI_SU3(out[idst_bgqlx][mu],in[isrc_lx][mu],vn_bw_dst_bgqlx);
	      }
	  }
      }
  
  //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward 
  for(int mu=0;mu<4;mu++) //border and link direction
    if(paral_dir[mu])
      NISSA_PARALLEL_LOOP(ibord,loc_vol+bord_offset[mu],loc_vol+bord_offset[mu]+bord_dir_vol[mu])
	{
	  int idst_lx=loclx_neighup[ibord][mu];
	  int vn_dst_bgqlx=(idst_lx>=loc_volh);
	  int idst_bgqlx=bgqlx_of_loclx[idst_lx%loc_volh];
	  
	  SU3_TO_BI_SU3(out[idst_bgqlx][mu],in[ibord][mu],vn_dst_bgqlx);
	}
  
  set_borders_invalid(out);
}}

//remap a spincolor from lx to bgqlx layout
THREADABLE_FUNCTION_2ARG(lx_spincolor_remap_to_bgqlx, bi_spincolor*,ext_out, spincolor*,in)
{
  GET_THREAD_ID();
    
  //bufferize if needed
  int bufferize=(void*)ext_out==(void*)in;
  bi_spincolor *out=bufferize?nissa_malloc("out",loc_volh,bi_spincolor):ext_out;
  
  //copy the two VN
  for(int vn=0;vn<2;vn++)
    NISSA_PARALLEL_LOOP(ivol_lx,0,loc_volh)
      SPINCOLOR_TO_BI_SPINCOLOR(out[bgqlx_of_loclx[ivol_lx]],in[ivol_lx+vn*loc_volh],vn);
  
  //wait filling
  set_borders_invalid(out);
  
  //unbuffer if needed
  if(bufferize)
    {
      vector_copy(ext_out,out);
      nissa_free(out);
    }
}}
//reverse
THREADABLE_FUNCTION_2ARG(bgqlx_spincolor_remap_to_lx, spincolor*,ext_out, bi_spincolor*,in)
{
  GET_THREAD_ID();
  
  //buffer if needed
  int bufferize=(void*)ext_out==(void*)in;
  spincolor *out=bufferize?nissa_malloc("out",loc_vol,spincolor):ext_out;

  //split to the two VN
  NISSA_PARALLEL_LOOP(ivol_bgqlx,0,loc_volh)
    BI_SPINCOLOR_TO_SPINCOLOR(out[loclx_of_bgqlx[ivol_bgqlx]],out[loc_volh+loclx_of_bgqlx[ivol_bgqlx]],in[ivol_bgqlx]);
  
  //wait filling
  set_borders_invalid(out);

  //unbuffer if needed
  if(bufferize)
    {
      vector_copy(ext_out,out);
      nissa_free(out);
    }
}}

//set bgq geometry
void set_bgq_geometry()
{
  if(!nissa_eo_geom_inited)
    {
      if(!nissa_use_eo_geom) crash("eo geometry must be enabled in order to use bgq one");
      crash("initialize eo_geometry before bgq one");
    }
  
  define_bgq_lx_ordering();

  //allocate a temporary vector to apply hopping matrix
  bgq_hopping_matrix_output_data=nissa_malloc("bgq_hopping_matrix_output_data",8*loc_volh,bi_halfspincolor);

  //define output index pointers
  define_bgq_hopping_matrix_lx_output_pointers_and_T_buffers();
}

//unset it
void unset_bgq_geometry()
{
  nissa_free(bgqlx_of_loclx);
  nissa_free(loclx_of_bgqlx);
  
  nissa_free(bgq_hopping_matrix_output_pointer);
  nissa_free(bgq_hopping_matrix_output_T_buffer);
  nissa_free(bgq_hopping_matrix_output_data);
  nissa_free(bgq_hopping_matrix_final_output);
}
