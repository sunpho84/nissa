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
  
  We first scan surface, then bulk.
  Backward links are (0,1,2,3), forward links are (4,5,6,7)
*/
void define_bgq_lx_ordering()
{
  bgqlx_of_loclx=nissa_malloc("bgqlx_of_loclx",loc_volh,int);
  loclx_of_bgqlx=nissa_malloc("loclx_of_bgqlx",loc_volh,int);
  
  //scan surface
  int bgqlx=0;
  for(int isurflx=0;isurflx<surf_vol;isurflx++)
    {
      int loclx=loclx_of_surflx[isurflx];
      
      //check that we are on first half of time direction
      if(loc_coord_of_loclx[loclx][0]<loc_size[0]/2)
	{
	  loclx_of_bgqlx[bgqlx]=loclx;
	  bgqlx_of_loclx[loclx]=bgqlx;
	  
	  bgqlx++;
	}
    }

  //scan bulk
  for(int ibulklx=0;ibulklx<bulk_vol;ibulklx++)
    {
      int loclx=loclx_of_bulklx[ibulklx];
      
      //check that we are on first half of time direction
      if(loc_coord_of_loclx[loclx][0]<loc_size[0]/2)
	{
	  loclx_of_bgqlx[bgqlx]=loclx;
	  bgqlx_of_loclx[loclx]=bgqlx;
	  
	  bgqlx++;
	}
    }
  
  if(bgqlx!=loc_volh) crash("defining bgq_lx ordering");
}

/*
  Define output index of sink-applied hopping matrix for non-eo preco operator.
  First we put T border, which is always present and of the same size of the 
  true node T border, then other three borders which are only half the size of the
  true node border. As always, forward border comes after backward one.
  Also backward derivative output comes before forward.
*/
void define_bgq_hopping_matrix_output_index()
{
  bgqlx_t_vbord_vol=loc_vol/loc_size[0]; //t dir is at least virtually parallelized
  bgqlx_vbord_vol=2*(bgqlx_t_vbord_vol+(bord_dir_vol[1]+bord_dir_vol[2]+bord_dir_vol[3])/2);
  
  /*debug
    master_printf("bgqlx_t_vbord_vol %d, bgqlx_vbord_vol %d\n",bgqlx_t_vbord_vol,bgqlx_vbord_vol);*/
  
  //allocate out hopping matrix index
  bgq_hopping_matrix_output_index=nissa_malloc("bgq_hopping_matrix_output_index",8*loc_volh,int);
  
  //for debug
  for(int i=0;i<8*loc_volh;i++)
    bgq_hopping_matrix_output_index[i]=-1;
  
  for(int ibgqlx=0;ibgqlx<loc_volh;ibgqlx++)
    {
      int iloclx=loclx_of_bgqlx[ibgqlx];
      
      //t direction bw derivative
      bgq_hopping_matrix_output_index[ibgqlx*8+0]=(loc_coord_of_loclx[iloclx][0]==0)?
	iloclx:                                                       //we move in other vnode
	bgqlx_vbord_vol+8*bgqlx_of_loclx[loclx_neighdw[iloclx][0]];   //we are still in the same vnode
      /*debug
	master_printf("A (%d) bgq_hopping_matrix_output_index[%d]: %d\n",loc_coord_of_loclx[iloclx][0]==0,ibgqlx*8+0,
	bgq_hopping_matrix_output_index[ibgqlx*8+0]);*/
      
      //t direction fw derivative
      bgq_hopping_matrix_output_index[ibgqlx*8+4]=(loc_coord_of_loclx[iloclx][0]==loc_size[0]/2-1)?
	bgqlx_vbord_vol/2+loclx_neighup[iloclx][0]-loc_volh:           //we moved in another vnode
	bgqlx_vbord_vol+8*bgqlx_of_loclx[loclx_neighup[iloclx][0]]+4;  //we are still in the same vnode
      /*debug
	master_printf("B (%d) bgq_hopping_matrix_output_index[%d]: %d\n",loc_coord_of_loclx[iloclx][0]==loc_size[0]/2-1,
	ibgqlx*8+4,
	bgq_hopping_matrix_output_index[ibgqlx*8+4]);*/
      
      //other direction derivatives
      for(int idir=1;idir<4;idir++)
	{
	  int bw=loclx_neighdw[iloclx][idir];
	  bgq_hopping_matrix_output_index[ibgqlx*8+0+idir]=(bw>=loc_vol)?
	    bgqlx_t_vbord_vol+(bord_offset[idir]-bord_offset[1])/2+
	    (loc_coord_of_loclx[iloclx][perp_dir[idir][2]]+loc_size[perp_dir[idir][2]]*
	     (loc_coord_of_loclx[iloclx][perp_dir[idir][1]]+loc_size[perp_dir[idir][1]]*
	      (loc_coord_of_loclx[iloclx][0]))):           //we moved in another vnode
	    bgqlx_vbord_vol+bgqlx_of_loclx[bw]*8+idir;     //we are still in local vnode
	  
	  /*debug
	    master_printf("C (%d  %d dir, %d %d %d) bgq_hopping_matrix_output_index[%d]: %d\n",bw>=loc_vol,
	    idir,
	    loc_coord_of_loclx[iloclx][0],
	    loc_coord_of_loclx[iloclx][perp_dir[idir][1]],
	    loc_coord_of_loclx[iloclx][perp_dir[idir][2]],
	    ibgqlx*8+0+idir,
	    bgq_hopping_matrix_output_index[ibgqlx*8+0+idir]);*/
	  
	  int fw=loclx_neighup[iloclx][idir];
	  bgq_hopping_matrix_output_index[ibgqlx*8+4+idir]=(fw>=loc_vol)?
	    bgqlx_vbord_vol/2+bgqlx_t_vbord_vol+(bord_offset[idir]-bord_offset[1])/2+
	    (loc_coord_of_loclx[iloclx][perp_dir[idir][2]]+loc_size[perp_dir[idir][2]]*
	     (loc_coord_of_loclx[iloclx][perp_dir[idir][1]]+loc_size[perp_dir[idir][1]]*
	      (loc_coord_of_loclx[iloclx][0]))):  //we moved in another vnode
	    bgqlx_vbord_vol+bgqlx_of_loclx[fw]*8+4+idir;                               //we are still in local vnode
	  
	  /*debug
	    master_printf("D (%d) bgq_hopping_matrix_output_index[%d]: %d\n",fw>=loc_vol,ibgqlx*8+4+idir,
	    bgq_hopping_matrix_output_index[ibgqlx*8+4+idir]);*/
	}
    }
  
  /*debug
    for(int i=0;i<8*loc_volh;i++)
    if(bgq_hopping_matrix_output_index[i]==-1) crash("index %d was not properly defined");
    
    for(int i=0;i<8*loc_volh;i++)
    master_printf("site*8+dirversed %d will output to index %d\n",i,bgq_hopping_matrix_output_index[i]);*/    
}

//we must put together the 8 links to be applied to a single point
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
	
	//we put together the four forward links
	for(int mu=0;mu<4;mu++)
	  {
	    //catch fw links
	    SU3_TO_BI_SU3(out[ifw_dst_bgqlx][4+mu],in[isrc_lx][mu],vn_fw_dst_bgqlx);

	    /*debug
	      master_printf("%04d (%d,%d) remapA to (%d,%d,%d)\n",ifw_dst_bgqlx*16+(4+mu)*2+vn_fw_dst_bgqlx,isrc_lx,mu,
	      ifw_dst_bgqlx,4+mu,vn_fw_dst_bgqlx);
	      tot_fixed++;*/
	    
	    //catch bw links
	    int idst_lx=loclx_neighup[isrc_lx][mu]; //bw links are used by fw sites bw der
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
  
  //scan the backward borders (first half of lx border)
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
  
  thread_barrier(REMAP_BARRIER);
  
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

  thread_barrier(REMAP_BARRIER);
}}
//reverse
THREADABLE_FUNCTION_2ARG(bgqlx_spincolor_remap_to_lx, spincolor*,out, bi_spincolor*,in)
{
  GET_THREAD_ID();
  
  NISSA_PARALLEL_LOOP(ivol_bgqlx,0,loc_volh)
    BI_SPINCOLOR_TO_SPINCOLOR(out[loclx_of_bgqlx[ivol_bgqlx]],out[loc_volh+loclx_of_bgqlx[ivol_bgqlx]],in[ivol_bgqlx]);
  
  thread_barrier(REMAP_BARRIER);
}}
