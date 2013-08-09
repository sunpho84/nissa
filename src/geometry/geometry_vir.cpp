#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

//for the moment we still call "bgq"
#include "../bgq/bgq_macros.h"

#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/thread_macros.h"
#include "../base/vectors.h"
#include "../communicate/communicate.h"
#include "../new_types/complex.h"

#ifdef USE_THREADS
 #include "../routines/thread.h"
#endif

inline int vnode_of_loclx(int lx)
{return nvnodes*loc_coord_of_loclx[lx][nissa_vnode_paral_dir]/loc_size[nissa_vnode_paral_dir];}

/*
  Define virtual nodes ordering: surface sites comes first, then bulk sites.
  Let's call v the direction where we have virtual nodes.
  Each site at coord c is matched with site at coord c+L_v/2
  Mapping array are defined only for first virtual node, other sites are
  separated by the first (in the lx layout) by vnode_lx_offset.
  We first scan v bw and fw surface, then others, then bulk.
  Putting v border in front makes easier putting back in place vhalo data
*/
void define_virlx_ordering()
{
  int v=nissa_vnode_paral_dir;
  
  virlx_of_loclx=nissa_malloc("virlx_of_loclx",loc_vol+bord_vol,int);
  loclx_of_virlx=nissa_malloc("loclx_of_virlx",(loc_vol+bord_vol)/nvnodes,int);
  
  //reset virtual index
  int virlx=0;
  
  //scan v- and + surface, in order
  int coord_to_compare[2]={0,loc_size[v]/nvnodes-1};
  for(int iter=0;iter<2;iter++)
    for(int loclx=0;loclx<loc_vol;loclx++)
      if(loc_coord_of_loclx[loclx][v]==coord_to_compare[iter])
	{
	  loclx_of_virlx[virlx]=loclx;
	  for(int inode=0;inode<nvnodes;inode++) virlx_of_loclx[loclx+inode*vnode_lx_offset]=virlx;
	  virlx++;
	}
  
  //scan non-v surface
  for(int isurflx=0;isurflx<surf_vol;isurflx++)
    {
      int loclx=loclx_of_surflx[isurflx];
      
      //check that we are on first half of v direction, and that we are not on x_v==0 or x_v==loc_size[v]/nvnodes-1
      if(loc_coord_of_loclx[loclx][v]!=0 && loc_coord_of_loclx[loclx][v]<loc_size[v]/nvnodes-1)
        {
          loclx_of_virlx[virlx]=loclx;
	  for(int inode=0;inode<nvnodes;inode++) virlx_of_loclx[loclx+inode*vnode_lx_offset]=virlx;

          virlx++;
        }
    }

  //scan bulk
  for(int ibulklx=0;ibulklx<bulk_vol;ibulklx++)
    {
      int loclx=loclx_of_bulklx[ibulklx];
      
      //check that we are on first half and not on the surf, even virtual
      if(loc_coord_of_loclx[loclx][v]!=0 && loc_coord_of_loclx[loclx][v]<loc_size[v]/nvnodes-1)
	{
	  loclx_of_virlx[virlx]=loclx;
	  for(int inode=0;inode<nvnodes;inode++) virlx_of_loclx[loclx+inode*vnode_lx_offset]=virlx;

	  virlx++;
	}
    }
  
  if(virlx!=loc_vol/nvnodes) crash("defining virlx ordering: %d!=%d",virlx,loc_vol/nvnodes);
  
  //define virtual size
  for(int mu=0;mu<4;mu++) vir_loc_size[mu]=loc_size[mu];
  vir_loc_size[v]/=nvnodes;
  
  //scan bw and fw borders
  for(int bf=0;bf<2;bf++)
    for(int mu0=0;mu0<4;mu0++)
      {
	//take dirs
	int mu1=perp_dir[mu0][0],mu2=perp_dir[mu0][1],mu3=perp_dir[mu0][2];
	int l1=vir_loc_size[mu1],l2=vir_loc_size[mu2],l3=vir_loc_size[mu3];
	for(int base_bordlx=0;base_bordlx<bord_dir_vol[mu0];base_bordlx++)
	  {
	    int loclx=loc_vol+bord_offset[mu0]+bord_volh*bf+base_bordlx;
	    
	    //take coord of bordlx including nvnodes factor in v dir
	    int surflx=surflx_of_bordlx[bord_offset[mu0]+bord_volh*bf+base_bordlx];
	    coords c;
	    for(int nu=0;nu<4;nu++) c[nu]=loc_coord_of_loclx[surflx][nu];
	    c[v]%=vir_loc_size[v];
	    
	    //define destination in vir geometry, pairing sites in v dir
	    int temp=c[mu3]+l3*(c[mu2]+l2*c[mu1]);
	    if(mu0==v) temp/=nvnodes;
	    int virlx=(loc_vol+bord_volh*bf+bord_offset[mu0])/nvnodes+temp;
	    virlx_of_loclx[loclx]=virlx;
	    if(virlx<(loc_vol+bord_vol)/nvnodes) loclx_of_virlx[virlx]=loclx;
	  }
      }
}

/*
  Define output index of sink-applied hopping matrix for non-eo preco operator, using virtual nodes.
  See "two_stage_computations" doc for more explenations.
  In v virtual halo data is ordered as for buffered communication send buffer, that is,
  in the same order as explained in "communicate.h"
*/
void define_virlx_hopping_matrix_output_pos()
{
  //order will be:     bord_vol/nvnodes | 8*loc_vol/nvnodes | vbord_vol
  int loc_data_start=bord_vol/nvnodes;
  int vbord_start=loc_data_start+8*loc_vol/nvnodes;
  int req_size=vbord_start+vbord_vol; //note that this is in unity of the vparallelized structure
  if(nissa_send_buf_size<req_size*sizeof(bi_halfspincolor)) crash("we need larger nissa_send_buf"); //to be moved
  
  two_stage_computation_pos_t *out=&virlx_hopping_matrix_output_pos; //for brevity
  out->inter_fr_in_pos=nissa_malloc("inter_fr_in_pos",8*loc_vol/nvnodes,int);
  out->final_fr_inter_pos=NULL; //we are not overlapping communication with intermediate->final transfer, for the moment
  out->inter_fr_recv_pos=nissa_malloc("inter_fr_recv_pos",bord_vol/nvnodes,int);
  
  int v=nissa_vnode_paral_dir;
  int mu1=perp_dir[v][0],mu2=perp_dir[v][1],mu3=perp_dir[v][2];
  
  for(int ivirlx=0;ivirlx<loc_vol/nvnodes;ivirlx++)
    {
      int iloclx=loclx_of_virlx[ivirlx];
      int *c=loc_coord_of_loclx[iloclx];
      //v direction bw scattering (fw derivative)
      out->inter_fr_in_pos[ivirlx*8+0+v]=(loc_coord_of_loclx[iloclx][v]==0)?
	//we move in other vnode, so we must point to local position in the v plane
	vbord_start+c[mu3]+loc_size[mu3]*(c[mu2]+loc_size[mu2]*c[mu1]):
	//we are still in the same vnode
	loc_data_start+8*virlx_of_loclx[loclx_neighdw[iloclx][v]]+0+v;

      //v direction fw scattering (bw derivative)
      out->inter_fr_in_pos[ivirlx*8+4+v]=(loc_coord_of_loclx[iloclx][v]==loc_size[v]/nvnodes-1)?
	//idem, but we shift of half vbord because this is + bord
	vbord_start+vbord_vol/2+c[mu3]+loc_size[mu3]*(c[mu2]+loc_size[mu2]*c[mu1]):
	//we are still in the same vnode, but we shift of 4 (this is + contr)
	loc_data_start+8*virlx_of_loclx[loclx_neighup[iloclx][v]]+4+v;
      
      //other direction derivatives
      for(int imu=0;imu<3;imu++)
	{
	  int mu=perp_dir[v][imu];

	  int bw=loclx_neighdw[iloclx][mu];
	  out->inter_fr_in_pos[ivirlx*8+0+mu]=(bw>=loc_vol)?
	    //we moved to another node
	    virlx_of_loclx[bw]-loc_vol/nvnodes:
	    //we are still in local vnode
	    loc_data_start+8*virlx_of_loclx[bw]+0+mu;
	    
	    int fw=loclx_neighup[iloclx][mu];
	    out->inter_fr_in_pos[ivirlx*8+4+mu]=(fw>=loc_vol)?
	      //we moved in another node
	      virlx_of_loclx[fw]-loc_vol/nvnodes:
	      //we are still in local vnode
	      loc_data_start+8*virlx_of_loclx[fw]+4+mu;
	}
    }
  
  //"v" border is consecutive, in the vhalo as in the border, so we just need to define other three borders
  for(int imu=0;imu<3;imu++)
    {
      int mu=perp_dir[v][imu];
      for(int base_src=0;base_src<bord_dir_vol[mu]/nvnodes;base_src++)
	{
	  //other 3 bw borders
	  int bw_vir_src=bord_offset[mu]/nvnodes+base_src;
	  int bw_bordlx_src=loclx_of_virlx[bw_vir_src+loc_vol/nvnodes]-loc_vol;
	  int bw_dst=8*virlx_of_loclx[surflx_of_bordlx[bw_bordlx_src]]+4+mu;
	  out->inter_fr_recv_pos[bw_vir_src]=bw_dst;
	  
	  //other 3 fw borders
	  int fw_vir_src=bord_volh/nvnodes+bord_offset[mu]/nvnodes+base_src;
	  int fw_bordlx_src=loclx_of_virlx[fw_vir_src+loc_vol/nvnodes]-loc_vol;
	  int fw_dst=8*virlx_of_loclx[surflx_of_bordlx[fw_bordlx_src]]+mu;
	  out->inter_fr_recv_pos[fw_vir_src]=fw_dst;
	}
    }
}

/*
  put together the 8 links to be applied to a single point
  first comes the links needed to scatter backward the signal (not to be daggered)
  then those needed to scatter it forward (to be daggered)
*/
THREADABLE_FUNCTION_2ARG(lx_conf_remap_to_virlx, bi_oct_su3*,out, quad_su3*,in)
{
  GET_THREAD_ID();
  
  int v=nissa_vnode_paral_dir;
  
  //communicate conf border so we can accede it
  communicate_lx_quad_su3_borders(in);
  
  //scan the two virtual nodes
  NISSA_PARALLEL_LOOP(isrc_lx,0,loc_vol)
    {
      int ifw_dst_virlx=virlx_of_loclx[isrc_lx];
      
      for(int mu=0;mu<4;mu++)
	{
	  //catch links needed to scatter signal forward
	  SU3_TO_BI_SU3(out[ifw_dst_virlx][4+mu],in[isrc_lx][mu],vnode_of_loclx(isrc_lx));
	  
	  //copy links also where they are needed to scatter the signal backward, if 
	  //sites that need them are not in the border (that would mean that computation must be 
	  //done in another node
	  int idst_lx=loclx_neighup[isrc_lx][mu];
	  if(idst_lx<loc_vol)
	    {
	      int idst_virlx=virlx_of_loclx[idst_lx];
	      SU3_TO_BI_SU3(out[idst_virlx][mu],in[isrc_lx][mu],vnode_of_loclx(idst_lx));
	    }
	  }
      }
  
  //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward 
  for(int mu=0;mu<4;mu++) //border and link direction
    if(paral_dir[mu])
      NISSA_PARALLEL_LOOP(ibord,loc_vol+bord_offset[mu],loc_vol+bord_offset[mu]+bord_dir_vol[mu])
	{
	  int idst_lx=loclx_neighup[ibord][mu];
	  int vn_dst_virlx=vnode_of_loclx(idst_lx);
	  int idst_virlx=virlx_of_loclx[idst_lx];
	  
	  SU3_TO_BI_SU3(out[idst_virlx][mu],in[ibord][mu],vn_dst_virlx);
	}
  
  set_borders_invalid(out);
}}

//remap a spincolor from lx to virlx layout
THREADABLE_FUNCTION_2ARG(lx_spincolor_remap_to_virlx, bi_spincolor*,ext_out, spincolor*,in)
{
  GET_THREAD_ID();
    
  //bufferize if needed
  int bufferize=(void*)ext_out==(void*)in;
  bi_spincolor *out=bufferize?nissa_malloc("out",loc_vol/nvnodes,bi_spincolor):ext_out;
  
  //copy the various VN
  NISSA_PARALLEL_LOOP(ivol_lx,0,loc_vol)
    SPINCOLOR_TO_BI_SPINCOLOR(out[virlx_of_loclx[ivol_lx]],in[ivol_lx],vnode_of_loclx(ivol_lx));
  
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
THREADABLE_FUNCTION_2ARG(virlx_spincolor_remap_to_lx, spincolor*,ext_out, bi_spincolor*,in)
{
  GET_THREAD_ID();
  
  //buffer if needed
  int bufferize=(void*)ext_out==(void*)in;
  spincolor *out=bufferize?nissa_malloc("out",loc_vol,spincolor):ext_out;

  //split to the two VN
  NISSA_PARALLEL_LOOP(ivol_virlx,0,loc_vol/nvnodes)
    BI_SPINCOLOR_TO_SPINCOLOR(out[loclx_of_virlx[ivol_virlx]],out[loclx_of_virlx[ivol_virlx]+vnode_lx_offset],
			      in[ivol_virlx]);
  
  //wait filling
  set_borders_invalid(out);

  //unbuffer if needed
  if(bufferize)
    {
      vector_copy(ext_out,out);
      nissa_free(out);
    }
}}

//set virtual geometry
void set_vir_geometry()
{
  if(!nissa_eo_geom_inited)
    {
      if(!nissa_use_eo_geom) crash("eo geometry must be enabled in order to use vir one");
      crash("initialize eo_geometry before vir one");
      
      if(loc_size[nissa_vnode_paral_dir]%4) crash("virtual parallelized dir must have local size multiple of 4");
    }
  
  define_virlx_ordering();

  //define output index pointers
  define_virlx_hopping_matrix_output_pos();
}

//unset it
void unset_vir_geometry()
{
  virlx_hopping_matrix_output_pos.free();
  
  nissa_free(virlx_of_loclx);
  nissa_free(loclx_of_virlx);
}
