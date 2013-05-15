#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "nissa.h"

#include "bgq_macros.h"
#include "new_vars_and_types.h"

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
  in reversed order with respect to border geometry.
  The other three borders points directly to communication borders.
*/
void define_bgq_hopping_matrix_lx_output_pointers_and_T_buffers(bi_halfspincolor *binded)
{
  bgqlx_t_vbord_vol=2*loc_vol/loc_size[0]; //t dir is at least virtually parallelized
  
  //bind
  bgq_hopping_matrix_output_binded=binded;
  
  //allocate hopping matrix output pointer
  bgq_hopping_matrix_output_pointer=nissa_malloc("bgq_hopping_matrix_output_pointer",8*loc_volh,bi_halfspincolor*);
  bgq_hopping_matrix_output_T_buffer=nissa_malloc("bgq_hopping_matrix_output_T_buffr",bgqlx_t_vbord_vol,bi_halfspincolor);
  
  for(int ibgqlx=0;ibgqlx<loc_volh;ibgqlx++)
    {
      int iloclx=loclx_of_bgqlx[ibgqlx];
      
      //t direction bw scattering
      bgq_hopping_matrix_output_pointer[ibgqlx*8+0]=(loc_coord_of_loclx[iloclx][0]==0)?
	bgq_hopping_matrix_output_T_buffer+iloclx:  //we move in other vnode
	binded+8*bgqlx_of_loclx[loclx_neighdw[iloclx][0]]+0;            //we are still in the same vnode
      
      //t direction fw scattering
      bgq_hopping_matrix_output_pointer[ibgqlx*8+4]=(loc_coord_of_loclx[iloclx][0]==loc_size[0]/2-1)?
	//we moved in another vnode
	bgq_hopping_matrix_output_T_buffer+bgqlx_t_vbord_vol/2+loclx_neighup[iloclx][0]-loc_volh:
	binded+8*bgqlx_of_loclx[loclx_neighup[iloclx][0]]+4;  //we are still in the same vnode
      
      //other direction derivatives
      for(int idir=1;idir<4;idir++)
	{
	  int bw=loclx_neighdw[iloclx][idir];
	  bgq_hopping_matrix_output_pointer[ibgqlx*8+0+idir]=(bw>=loc_vol)?
	    (bi_halfspincolor*)nissa_send_buf+bw-loc_vol-bord_offset[idir]/2:   //we moved to another vnode
	    binded+8*bgqlx_of_loclx[bw]+0+idir;                        //we are still in local vnode
	  
	  int fw=loclx_neighup[iloclx][idir];
	  bgq_hopping_matrix_output_pointer[ibgqlx*8+4+idir]=(fw>=loc_vol)?
	    (bi_halfspincolor*)nissa_send_buf+fw-loc_vol-bord_volh/2-bord_offset[idir]/2:              //we moved in another vnode
	    binded+8*bgqlx_of_loclx[fw]+4+idir;                        //we are still in local vnode
	}
    }
}

/*
  put together the 8 links to be applied to a single point
  first comes the links needed to scatter backward the signal (not to be daggered)
  then the ones needed to scatter it forward (to be daggered)
*/
THREADABLE_FUNCTION_2ARG(lx_conf_remap_to_bgqlx, bi_oct_su3*,out, quad_su3*,in)
{
  GET_THREAD_ID();
  
  //communicate conf border so we can accede it
  communicate_lx_quad_su3_borders(in);
  
  /*debug
    int tot_fixed=0;*/
  
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

	    /*debug
	      master_printf("%04d (%d,%d) remapA to (%d,%d,%d)\n",ifw_dst_bgqlx*16+(4+mu)*2+vn_fw_dst_bgqlx,isrc_lx,mu,
	      ifw_dst_bgqlx,4+mu,vn_fw_dst_bgqlx);
	      tot_fixed++;*/
	    
	    //catch links needed to scatter signal backward 
	    int idst_lx=loclx_neighup[isrc_lx][mu];
	    int vn_bw_dst_bgqlx=(idst_lx/loc_volh); //1 only if idst_lx is in the first node, 2 if in boder
	    
	    //discard border out
	    if(vn_bw_dst_bgqlx<2)
	      {
		int idst_bgqlx=bgqlx_of_loclx[idst_lx%loc_volh];
		SU3_TO_BI_SU3(out[idst_bgqlx][mu],in[isrc_lx][mu],vn_bw_dst_bgqlx);
		
		/*debug
		  master_printf("%04d (%d,%d) remapB to (%d,%d,%d)\n",idst_bgqlx*16+mu*2+vn_bw_dst_bgqlx,isrc_lx,mu,
		  idst_bgqlx,mu,vn_bw_dst_bgqlx);
		  tot_fixed++;*/
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
	  
	  /*debug
	    master_printf("%04d (%d,%d) remapC to (%d,%d,%d)\n",idst_bgqlx*16+mu*2+vn_dst_bgqlx,ibord,mu,
	    idst_bgqlx,mu,vn_dst_bgqlx);
	    tot_fixed++;*/
	}
  
  THREAD_BARRIER();
  
  /*debug
    if(tot_fixed!=8*loc_vol) crash("fixed only %d when expectd %d",tot_fixed,8*loc_vol);
    
    if(IS_MASTER_THREAD)
    for(int ivol=0;ivol<loc_volh;ivol++)
    for(int mu=0;mu<8;mu++)
    CHECK_BI_SU3(out[ivol][mu]);*/
}}

//remap a spincolor from lx to bgqlx layout
THREADABLE_FUNCTION_2ARG(lx_spincolor_remap_to_bgqlx, bi_spincolor*,out, spincolor*,in)
{
  GET_THREAD_ID();
  
  for(int vn=0;vn<2;vn++)
  NISSA_PARALLEL_LOOP(ivol_lx,0,loc_volh)
    SPINCOLOR_TO_BI_SPINCOLOR(out[bgqlx_of_loclx[ivol_lx]],in[ivol_lx+vn*loc_volh],vn);

  THREAD_BARRIER();
}}
//reverse
THREADABLE_FUNCTION_2ARG(bgqlx_spincolor_remap_to_lx, spincolor*,out, bi_spincolor*,in)
{
  GET_THREAD_ID();
  
  NISSA_PARALLEL_LOOP(ivol_bgqlx,0,loc_volh)
    BI_SPINCOLOR_TO_SPINCOLOR(out[loclx_of_bgqlx[ivol_bgqlx]],out[loc_volh+loclx_of_bgqlx[ivol_bgqlx]],in[ivol_bgqlx]);
  
  THREAD_BARRIER();
}}

//set bgq geometry
void set_bgq_geometry()
{
  define_bgq_lx_ordering();
}

//unset it
void unset_bgq_geometry()
{
  nissa_free(bgqlx_of_loclx);
  nissa_free(loclx_of_bgqlx);
}
