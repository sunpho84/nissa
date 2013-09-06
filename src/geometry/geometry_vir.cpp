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

//routine to add site in the list
void mark_vir_of_loclx(int &virlx,int *vireo,int &loclx)
{
  //take parity
  int par=loclx_parity[loclx];
  
  //mark loclx of virlx and vireo
  loclx_of_virlx[virlx]=loclx;
  loclx_of_vireo[par][vireo[par]]=loclx;
  
  //and virlx and vireo for all vnodes loclx
  for(int inode=0;inode<nvnodes;inode++)
    {
      virlx_of_loclx[loclx+inode*vnode_lx_offset]=virlx;
      vireo_of_loclx[loclx+inode*vnode_lx_offset]=vireo[par];
    }
  
  //increment virlx and vireo
  virlx++;
  vireo[par]++;
}

/*
  Define virtual nodes ordering: surface sites comes first, then bulk sites.
  Let's call v the direction where we have virtual nodes.
  Each site at coord c is matched with site at coord c+L_v/2
  Mapping array are defined only for first virtual node, other sites are
  separated by the first (in the lx layout) by vnode_lx_offset.
  We first scan v bw and fw surface, then others, then bulk.
  Putting v border in front makes easier putting back in place vhalo data
*/
void define_vir_ordering()
{
  //define virtual size
  int v=nissa_vnode_paral_dir;
  for(int mu=0;mu<4;mu++) vir_loc_size[mu]=loc_size[mu];
  vir_loc_size[v]/=nvnodes;
  
  //allocate
  virlx_of_loclx=nissa_malloc("virlx_of_loclx",loc_vol+bord_vol,int);
  vireo_of_loclx=nissa_malloc("vireo_of_loclx",loc_vol+bord_vol,int);
  loclx_of_virlx=nissa_malloc("loclx_of_virlx",(loc_vol+bord_vol)/nvnodes,int);
  for(int par=0;par<2;par++) loclx_of_vireo[par]=nissa_malloc("loceo_of_virlx",(loc_vol+bord_vol)/nvnodes/2,int);
  
  //reset virtual index
  int virlx=0;
  int vireo[2]={0,0};
  
  //scan v- and + surface, in order
  int coord_to_compare[2]={0,loc_size[v]/nvnodes-1};
  for(int iter=0;iter<2;iter++)
    for(int loclx=0;loclx<loc_vol;loclx++)
      if(loc_coord_of_loclx[loclx][v]==coord_to_compare[iter]) mark_vir_of_loclx(virlx,vireo,loclx);
  
  //scan non-v surface
  for(int isurflx=0;isurflx<surf_vol;isurflx++)
    {
      int loclx=loclx_of_surflx[isurflx];
      
      //check that we are on first half of v direction, and that we are not on x_v==0 or x_v==loc_size[v]/nvnodes-1
      int c=loc_coord_of_loclx[loclx][v];
      if(c!=0 && c<loc_size[v]/nvnodes-1) mark_vir_of_loclx(virlx,vireo,loclx);
    }
  
  //take note of vsurf_vol
  vsurf_vol=virlx;
  vsurf_volh=vireo[0];
  
  //scan bulk
  for(int ibulklx=0;ibulklx<bulk_vol;ibulklx++)
    {
      int loclx=loclx_of_bulklx[ibulklx];
      
      //check that we are on first half and not on the surf, even virtual
      int c=loc_coord_of_loclx[loclx][v];
      if(c!=0 && c<loc_size[v]/nvnodes-1) mark_vir_of_loclx(virlx,vireo,loclx);
    }
  
  //trivial checks
  if(virlx!=loc_vol/nvnodes) crash("defining virlx ordering: %d!=%d",virlx,loc_vol/nvnodes);
  if(vireo[EVN]!=loc_vol/nvnodes/2) crash("defining vireo[EVN] ordering: %d!=%d",vireo[EVN],loc_vol/nvnodes/2);
  if(vireo[ODD]!=loc_vol/nvnodes/2) crash("defining vireo[EVN] ordering: %d!=%d",vireo[ODD],loc_vol/nvnodes/2);
  
  //scan bw and fw borders
  for(int bf=0;bf<2;bf++)
    for(int mu0=0;mu0<4;mu0++)
      {
	//take dirs
	int mu1=perp_dir[mu0][0],mu2=perp_dir[mu0][1],mu3=perp_dir[mu0][2];
	int l2=vir_loc_size[mu2],l3=vir_loc_size[mu3];//l1=vir_loc_size[mu1];
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
	    int vireo=virlx/2;
	    virlx_of_loclx[loclx]=virlx;
	    vireo_of_loclx[loclx]=vireo;
	    if(virlx<(loc_vol+bord_vol)/nvnodes)
	      {
		loclx_of_virlx[virlx]=loclx;
		loclx_of_vireo[loclx_parity[loclx]][vireo]=loclx;
	      }
	  }
      }
}

/*
  Define output index of sink-applied hopping matrix for hoppint matrix, using virtual nodes.
  See "two_stage_computations" doc for more explenations.
  In v virtual halo data is ordered as for buffered communication send buffer, that is,
  in the same order as explained in "communicate.h"
  if par==2 it makes the computation for lx case
  if par==1 for the eo case
  if par==0 for the oe case
*/
void define_vir_hopping_matrix_output_pos()
{
  for(int par=0;par<=2;par++)
    {
      int fact=(par==2)?1:2;
      
      //order will be:     bord_vol/nvnodes | 8*loc_vol/nvnodes | vbord_vol
      int loc_data_start=bord_vol/nvnodes/fact;
      int vbord_start=loc_data_start+8*loc_vol/nvnodes/fact;
      int req_size=vbord_start+vbord_vol/fact; //note that this is in unity of the vparallelized structure
      if(nissa_send_buf_size<req_size*sizeof(bi_halfspincolor)) //to be moved
	crash("we need larger nissa_send_buf: %d > %d",req_size*sizeof(bi_halfspincolor),nissa_send_buf_size);
      
      int *loclx_of_vir=(par==2)?loclx_of_virlx:loclx_of_vireo[par];
      int *vir_of_loclx=(par==2)?virlx_of_loclx:vireo_of_loclx;
      
      two_stage_computation_pos_t *out=(par==2)?&virlx_hopping_matrix_output_pos:
	viroe_or_vireo_hopping_matrix_output_pos+par;
      out->inter_fr_in_pos=nissa_malloc("inter_fr_in_pos",8*loc_vol/nvnodes/fact,int);
      out->final_fr_inter_pos=NULL; //we are not overlapping communication with intermediate->final transfer
      out->inter_fr_recv_pos=nissa_malloc("inter_fr_recv_pos",bord_vol/nvnodes/fact,int);
  
      int v=nissa_vnode_paral_dir;
      int mu1=perp_dir[v][0],mu2=perp_dir[v][1],mu3=perp_dir[v][2];
      
      for(int ivir=0;ivir<loc_vol/nvnodes/fact;ivir++)
	{
	  int iloclx=loclx_of_vir[ivir];
	  int *c=loc_coord_of_loclx[iloclx];
	  //v direction bw scattering (fw derivative)
	  out->inter_fr_in_pos[ivir*8+0+v]=(loc_coord_of_loclx[iloclx][v]==0)?
	    //we move in other vnode, so we must point to local position in the v plane
	    vbord_start+(vbord_volh*0+c[mu3]+loc_size[mu3]*(c[mu2]+loc_size[mu2]*c[mu1]))/fact:
	    //we are still in the same vnode
	    loc_data_start+8*vir_of_loclx[loclx_neighdw[iloclx][v]]+0+v;
	  /*
	  if(0 && par==0 && rank==0)
	    printf("ANNA_FWDER_BWSCAT_%d(v) (loc %d, %d %d %d %d glb %d, bw %d, %d) %d %d\n",
		   v,iloclx,loc_coord_of_loclx[iloclx][0],loc_coord_of_loclx[iloclx][1],
		   loc_coord_of_loclx[iloclx][2],loc_coord_of_loclx[iloclx][3],
		   glblx_of_loclx[iloclx],loclx_neighdw[iloclx][v],
		   (loc_coord_of_loclx[iloclx][v]==0),
		   ivir*8+0+v,out->inter_fr_in_pos[ivir*8+0+v]);
	  */
	  
	  //v direction fw scattering (bw derivative)
	  out->inter_fr_in_pos[ivir*8+4+v]=(loc_coord_of_loclx[iloclx][v]==loc_size[v]/nvnodes-1)?
	    //idem, but we shift of half vbord because this is + bord
	    vbord_start+(vbord_volh*1+c[mu3]+loc_size[mu3]*(c[mu2]+loc_size[mu2]*c[mu1]))/fact:
	    //we are still in the same vnode, but we shift of 4 (this is + contr)
	    loc_data_start+8*vir_of_loclx[loclx_neighup[iloclx][v]]+4+v;
	  /*
	  if(0 && par==0 && rank==0)
	    printf("ANNA_BWDER_FWSCAT_%d(v) (loc %d, %d %d %d %d, glb %d, fw %d, %d) %d %d\n",
		   v,iloclx,loc_coord_of_loclx[iloclx][0],loc_coord_of_loclx[iloclx][1],
		   loc_coord_of_loclx[iloclx][2],loc_coord_of_loclx[iloclx][3],
		   glblx_of_loclx[iloclx],loclx_neighup[iloclx][v],
		   (loc_coord_of_loclx[iloclx][v]==loc_size[v]/nvnodes-1),
		   ivir*8+4+v,out->inter_fr_in_pos[ivir*8+4+v]);
	  */
	  
	  //other direction derivatives
	  for(int imu=0;imu<3;imu++)
	    {
	      int mu=perp_dir[v][imu];
	      
	      int bw=loclx_neighdw[iloclx][mu];
	      out->inter_fr_in_pos[ivir*8+0+mu]=(bw>=loc_vol)?
		//we moved to another node
		vir_of_loclx[bw]-loc_vol/nvnodes/fact: //SEEMS THAT IT IS MAKING A MESS
		//we are still in local vnode
		loc_data_start+8*vir_of_loclx[bw]+0+mu;
	      /*
	      if(0 && par==0 && rank==0)
		printf("ANNA_FWDER_BWSCAT_%d (loc %d, %d %d %d %d, glb %d, bw %d, %d) %d %d\n",
		       mu,iloclx,loc_coord_of_loclx[iloclx][0],loc_coord_of_loclx[iloclx][1],
		       loc_coord_of_loclx[iloclx][2],loc_coord_of_loclx[iloclx][3],
		       glblx_of_loclx[iloclx],bw,
		       (bw>=loc_vol),
		       ivir*8+0+mu,out->inter_fr_in_pos[ivir*8+0+mu]);
	      */
	      
	      int fw=loclx_neighup[iloclx][mu];
	      out->inter_fr_in_pos[ivir*8+4+mu]=(fw>=loc_vol)?
		//we moved in another node
		vir_of_loclx[fw]-loc_vol/nvnodes/fact:
		//we are still in local vnode
		loc_data_start+8*vir_of_loclx[fw]+4+mu;
	      /*
	      if(0 && par==0 && rank==0)
		printf("ANNA_BWDER_FWSCAT_%d (loc %d, %d %d %d %d, glb %d, fw %d, %d) %d %d\n",
		       mu,iloclx,loc_coord_of_loclx[iloclx][0],loc_coord_of_loclx[iloclx][1],
		       loc_coord_of_loclx[iloclx][2],loc_coord_of_loclx[iloclx][3],
		       glblx_of_loclx[iloclx],fw,
		       (fw>=loc_vol),
		       ivir*8+4+mu,out->inter_fr_in_pos[ivir*8+4+mu]);
	      */
	    }
	}
      
      //"v" border is consecutive, in the vhalo as in the border, so we just need to define other three borders
      for(int imu=0;imu<3;imu++)
	{
	  int mu=perp_dir[v][imu];
	  for(int base_src=0;base_src<bord_dir_vol[mu]/nvnodes/fact;base_src++)
	    {
	      //other 3 bw borders
	      int bw_vir_src=(0*bord_volh+bord_offset[mu])/nvnodes/fact+base_src;
	      int bw_bordlx_src=loclx_of_vir[bw_vir_src+loc_vol/nvnodes/fact]-loc_vol;
	      int bw_dst=8*vir_of_loclx[surflx_of_bordlx[bw_bordlx_src]]+4+mu;
	      out->inter_fr_recv_pos[bw_vir_src]=bw_dst;
	      
	      //other 3 fw borders
	      int fw_vir_src=(1*bord_volh+bord_offset[mu])/nvnodes/fact+base_src;
	      int fw_bordlx_src=loclx_of_vir[fw_vir_src+loc_vol/nvnodes/fact]-loc_vol;
	      int fw_dst=8*vir_of_loclx[surflx_of_bordlx[fw_bordlx_src]]+0+mu;
	      out->inter_fr_recv_pos[fw_vir_src]=fw_dst;
	    }
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
  
  //communicate conf border so we can accede it
  communicate_lx_quad_su3_borders(in);
  
  //scan the two virtual nodes
  NISSA_PARALLEL_LOOP(isrc_lx,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      {
	//catch links needed to scatter signal forward
	SU3_TO_BI_SU3(out[virlx_of_loclx[isrc_lx]][4+mu],in[isrc_lx][mu],vnode_of_loclx(isrc_lx));
	
	//copy links also where they are needed to scatter the signal backward, if 
	//sites that need them are not in the border (that would mean that computation must be 
	//done in another node
	int idst_lx=loclx_neighup[isrc_lx][mu];
	if(idst_lx<loc_vol) SU3_TO_BI_SU3(out[virlx_of_loclx[idst_lx]][mu],in[isrc_lx][mu],vnode_of_loclx(idst_lx));
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

//similar for eo
THREADABLE_FUNCTION_2ARG(lx_conf_remap_to_vireo, bi_oct_su3**,out, quad_su3*,in)
{
  GET_THREAD_ID();
  
  const int debug=0;
  
  //communicate conf border so we can accede it
  communicate_lx_quad_su3_borders(in);
  
  //scan the two virtual nodes
  NISSA_PARALLEL_LOOP(isrc_lx,0,loc_vol)
    for(int mu=0;mu<4;mu++)
      {
	//catch links needed to scatter signal forward
	if(debug) in[isrc_lx][mu][0][0][1]=4+mu;
	SU3_TO_BI_SU3(out[!loclx_parity[isrc_lx]][vireo_of_loclx[isrc_lx]][4+mu],in[isrc_lx][mu],vnode_of_loclx(isrc_lx));
	if(debug) in[isrc_lx][mu][0][0][1]=0;

	//copy links also where they are needed to scatter the signal backward, if 
	//sites that need them are not in the border (that would mean that computation must be 
	//done in another node)
	int idst_lx=loclx_neighup[isrc_lx][mu],vn=vnode_of_loclx(idst_lx);
	if(idst_lx<loc_vol)
	  {
	    if(debug) in[isrc_lx][mu][0][0][1]=mu;
	    SU3_TO_BI_SU3(out[!loclx_parity[idst_lx]][vireo_of_loclx[idst_lx]][mu],in[isrc_lx][mu],vn);
	    if(debug) in[isrc_lx][mu][0][0][1]=0;
	  }
      }
  
  //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward 
  for(int mu=0;mu<4;mu++) //border and link direction
    if(paral_dir[mu])
      NISSA_PARALLEL_LOOP(ibord,loc_vol+bord_offset[mu],loc_vol+bord_offset[mu]+bord_dir_vol[mu])
	{
	  int idst_lx=loclx_neighup[ibord][mu],vn=vnode_of_loclx(idst_lx);
	  
	  if(debug) in[ibord][mu][0][0][1]=mu;
	  SU3_TO_BI_SU3(out[!loclx_parity[idst_lx]][vireo_of_loclx[idst_lx]][mu],in[ibord][mu],vn);
	  if(debug) in[ibord][mu][0][0][1]=0;
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

//remap a spincolor from lx to vireo layout
THREADABLE_FUNCTION_2ARG(lx_spincolor_remap_to_vireo, bi_spincolor**,out, spincolor*,in)
{
  GET_THREAD_ID();
    
  //copy the various VN
  NISSA_PARALLEL_LOOP(ivol_lx,0,loc_vol)
    SPINCOLOR_TO_BI_SPINCOLOR(out[loclx_parity[ivol_lx]][vireo_of_loclx[ivol_lx]],in[ivol_lx],vnode_of_loclx(ivol_lx));
  
  //wait filling
  set_borders_invalid(out[EVN]);
  set_borders_invalid(out[ODD]);
}}
//reverse
THREADABLE_FUNCTION_2ARG(vireo_spincolor_remap_to_lx, spincolor*,out, bi_spincolor**,in)
{
  GET_THREAD_ID();

  //split to the two VN
  for(int par=0;par<2;par++)
    NISSA_PARALLEL_LOOP(ivol_vireo,0,loc_volh/nvnodes)
      BI_SPINCOLOR_TO_SPINCOLOR(out[loclx_of_vireo[par][ivol_vireo]],
				out[loclx_of_vireo[par][ivol_vireo]+vnode_lx_offset],
				in[par][ivol_vireo]);
  
  //wait filling
  set_borders_invalid(out);
}}

//remap a color from lx to vireo layout
THREADABLE_FUNCTION_2ARG(lx_color_remap_to_vireo, bi_color**,out, color*,in)
{
  GET_THREAD_ID();
    
  //copy the various VN
  NISSA_PARALLEL_LOOP(ivol_lx,0,loc_vol)
    COLOR_TO_BI_COLOR(out[loclx_parity[ivol_lx]][vireo_of_loclx[ivol_lx]],in[ivol_lx],vnode_of_loclx(ivol_lx));
  
  //wait filling
  set_borders_invalid(out[EVN]);
  set_borders_invalid(out[ODD]);
}}
//reverse
THREADABLE_FUNCTION_2ARG(vireo_color_remap_to_lx, color*,out, bi_color**,in)
{
  GET_THREAD_ID();

  //split to the two VN
  for(int par=0;par<2;par++)
    NISSA_PARALLEL_LOOP(ivol_vireo,0,loc_volh/nvnodes)
      BI_COLOR_TO_COLOR(out[loclx_of_vireo[par][ivol_vireo]],out[loclx_of_vireo[par][ivol_vireo]+vnode_lx_offset],
				in[par][ivol_vireo]);
  
  //wait filling
  set_borders_invalid(out);
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
  
  define_vir_ordering();

  //define output index pointers
  define_vir_hopping_matrix_output_pos();
}

//unset it
void unset_vir_geometry()
{
  virlx_hopping_matrix_output_pos.free();
  
  nissa_free(loclx_of_virlx);
  nissa_free(virlx_of_loclx);

  for(int par=0;par<2;par++)
    {
      viroe_or_vireo_hopping_matrix_output_pos[par].free();
      nissa_free(loclx_of_vireo[par]);
    }
  nissa_free(vireo_of_loclx);
}
