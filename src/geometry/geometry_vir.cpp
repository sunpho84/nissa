#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_GEOMETRY_VIR
#include "geometry_vir.hpp"

#include <functional>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
//#include "bgq/Wilson_hopping_matrix_lx_bgq.hpp"
#include "communicate/borders.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_vir.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "new_types/two_stage_computation.hpp"
#include "operations/su3_paths/topological_charge.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  /*
    Define output index of sink-applied hopping matrix for hoppint matrix, using virtual nodes.
    See "two_stage_computations" doc for more explenations.
    In v virtual halo data is ordered as for buffered communication send buffer, that is,
    in the same order as explained in "communicate.h"
    if par==2 it makes the computation for lx case
    if par==1 for the eo case
    if par==0 for the oe case
  */
  // void define_vir_hopping_matrix_output_pos()
  // {
  //   for(int par=0;par<=2;par++)
  //     {
  // 	int fact=(par==2)?1:2;
	
  // 	//order will be:     bord_vol/nvranks | 8*vloc_vol/nvranks | vbord_vol
  // 	int loc_data_start=bord_vol/nvranks/fact;
  // 	int vbord_start=loc_data_start+8*loc_vol/nvranks/fact;
	
  // 	//note that this is in unity of the vparallelized structure, here assumed to be halfspincolor
  // 	size_t req_size=(vbord_start+vbord_vol/fact)*sizeof(vir_halfspincolor);
  // 	recv_buf_size=std::max(recv_buf_size,req_size);
  // 	send_buf_size=std::max(send_buf_size,req_size);
	
  // 	int *loclx_of_vir=(par==2)?loclx_of_virlx:loclx_of_vireo[par];
  // 	int *vir_of_loclx=(par==2)?virlx_of_loclx:vireo_of_loclx;
	
  // 	two_stage_computation_pos_t *out=(par==2)?&virlx_hopping_matrix_output_pos:
  // 	  viroe_or_vireo_hopping_matrix_output_pos+par;
  // 	out->inter_fr_in_pos=nissa_malloc("inter_fr_in_pos",8*loc_vol/nvranks/fact,int);
  // 	out->final_fr_inter_pos=NULL; //we are not overlapping communication with intermediate->final transfer
  // 	out->inter_fr_recv_pos=nissa_malloc("inter_fr_recv_pos",bord_vol/nvranks/fact,int);
	
  // 	int v=vrank_paral_dir;
  // 	int mu1=perp_dir[v][0],mu2=perp_dir[v][1],mu3=perp_dir[v][2];
	
  // 	for(int ivir=0;ivir<loc_vol/nvranks/fact;ivir++)
  // 	  {
  // 	    int iloclx=loclx_of_vir[ivir];
  // 	    int *c=loc_coord_of_loclx[iloclx];
  // 	    //v direction bw scattering (fw derivative)
  // 	    out->inter_fr_in_pos[ivir*8+0+v]=(loc_coord_of_loclx[iloclx][v]==0)?
  // 	      //we move in other vrank, so we must point to local position in the v plane
  // 	      vbord_start+(vbord_volh*0+c[mu3]+loc_size[mu3]*(c[mu2]+loc_size[mu2]*c[mu1]))/fact:
  // 	      //we are still in the same vrank
  // 	      loc_data_start+8*vir_of_loclx[loclx_neighdw[iloclx][v]]+0+v;
  // 	    /*
  // 	      if(0 && par==0 && rank==0)
  // 	      printf("ANNA_FWDER_BWSCAT_%d(v) (loc %d, %d %d %d %d glb %d, bw %d, %d) %d %d\n",
  // 	      v,iloclx,loc_coord_of_loclx[iloclx][0],loc_coord_of_loclx[iloclx][1],
  // 	      loc_coord_of_loclx[iloclx][2],loc_coord_of_loclx[iloclx][3],
  // 	      glblx_of_loclx[iloclx],loclx_neighdw[iloclx][v],
  // 	      (loc_coord_of_loclx[iloclx][v]==0),
  // 	      ivir*8+0+v,out->inter_fr_in_pos[ivir*8+0+v]);
  // 	    */
	    
  // 	    //v direction fw scattering (bw derivative)
  // 	    out->inter_fr_in_pos[ivir*8+4+v]=(loc_coord_of_loclx[iloclx][v]==loc_size[v]/nvranks-1)?
  // 	      //idem, but we shift of half vbord because this is + bord
  // 	      vbord_start+(vbord_volh*1+c[mu3]+loc_size[mu3]*(c[mu2]+loc_size[mu2]*c[mu1]))/fact:
  // 	      //we are still in the same vrank, but we shift of 4 (this is + contr)
  // 	      loc_data_start+8*vir_of_loclx[loclx_neighup[iloclx][v]]+4+v;
  // 	    /*
  // 	      if(0 && par==0 && rank==0)
  // 	      printf("ANNA_BWDER_FWSCAT_%d(v) (loc %d, %d %d %d %d, glb %d, fw %d, %d) %d %d\n",
  // 	      v,iloclx,loc_coord_of_loclx[iloclx][0],loc_coord_of_loclx[iloclx][1],
  // 	      loc_coord_of_loclx[iloclx][2],loc_coord_of_loclx[iloclx][3],
  // 	      glblx_of_loclx[iloclx],loclx_neighup[iloclx][v],
  // 	      (loc_coord_of_loclx[iloclx][v]==loc_size[v]/nvranks-1),
  // 	      ivir*8+4+v,out->inter_fr_in_pos[ivir*8+4+v]);
  // 	    */
	    
  // 	    //other direction derivatives
  // 	    for(int imu=0;imu<3;imu++)
  // 	      {
  // 		int mu=perp_dir[v][imu];
		
  // 		int bw=loclx_neighdw[iloclx][mu];
  // 		out->inter_fr_in_pos[ivir*8+0+mu]=(bw>=loc_vol)?
  // 		  //we moved to another node
  // 		  vir_of_loclx[bw]-loc_vol/nvranks/fact: //SEEMS THAT IT IS MAKING A MESS
  // 		  //we are still in local vrank
  // 		  loc_data_start+8*vir_of_loclx[bw]+0+mu;
  // 		/*
  // 		  if(0 && par==0 && rank==0)
  // 		  printf("ANNA_FWDER_BWSCAT_%d (loc %d, %d %d %d %d, glb %d, bw %d, %d) %d %d\n",
  // 		  mu,iloclx,loc_coord_of_loclx[iloclx][0],loc_coord_of_loclx[iloclx][1],
  // 		  loc_coord_of_loclx[iloclx][2],loc_coord_of_loclx[iloclx][3],
  // 		  glblx_of_loclx[iloclx],bw,
  // 		  (bw>=loc_vol),
  // 		  ivir*8+0+mu,out->inter_fr_in_pos[ivir*8+0+mu]);
  // 		*/
		
  // 		int fw=loclx_neighup[iloclx][mu];
  // 		out->inter_fr_in_pos[ivir*8+4+mu]=(fw>=loc_vol)?
  // 		  //we moved in another node
  // 		  vir_of_loclx[fw]-loc_vol/nvranks/fact:
  // 		  //we are still in local vrank
  // 		  loc_data_start+8*vir_of_loclx[fw]+4+mu;
  // 		/*
  // 		  if(0 && par==0 && rank==0)
  // 		  printf("ANNA_BWDER_FWSCAT_%d (loc %d, %d %d %d %d, glb %d, fw %d, %d) %d %d\n",
  // 		  mu,iloclx,loc_coord_of_loclx[iloclx][0],loc_coord_of_loclx[iloclx][1],
  // 		  loc_coord_of_loclx[iloclx][2],loc_coord_of_loclx[iloclx][3],
  // 		  glblx_of_loclx[iloclx],fw,
  // 		  (fw>=loc_vol),
  // 		  ivir*8+4+mu,out->inter_fr_in_pos[ivir*8+4+mu]);
  // 		*/
  // 	      }
  // 	  }
	
  // 	//"v" border is consecutive, in the vhalo as in the border, so we just need to define other three borders
  // 	for(int imu=0;imu<3;imu++)
  // 	  {
  // 	    int mu=perp_dir[v][imu];
  // 	    for(int base_src=0;base_src<bord_dir_vol[mu]/nvranks/fact;base_src++)
  // 	      {
  // 		//other 3 bw borders
  // 		int bw_vir_src=(0*bord_volh+bord_offset[mu])/nvranks/fact+base_src;
  // 		int bw_bordlx_src=loclx_of_vir[bw_vir_src+loc_vol/nvranks/fact]-loc_vol;
  // 		int bw_dst=8*vir_of_loclx[surflx_of_bordlx[bw_bordlx_src]]+4+mu;
  // 		out->inter_fr_recv_pos[bw_vir_src]=bw_dst;
		
  // 		//other 3 fw borders
  // 		int fw_vir_src=(1*bord_volh+bord_offset[mu])/nvranks/fact+base_src;
  // 		int fw_bordlx_src=loclx_of_vir[fw_vir_src+loc_vol/nvranks/fact]-loc_vol;
  // 		int fw_dst=8*vir_of_loclx[surflx_of_bordlx[fw_bordlx_src]]+0+mu;
  // 		out->inter_fr_recv_pos[fw_vir_src]=fw_dst;
  // 	      }
  // 	  }
  //     }
  // }
  
  // /*
  //   put together the 8 links to be applied to a single point
  //   first comes the links needed to scatter backward the signal (not to be daggered)
  //   then those needed to scatter it forward (to be daggered)
  // */
  // THREADABLE_FUNCTION_2ARG(lx_conf_remap_to_virlx, vir_oct_su3*,out, quad_su3*,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //communicate conf border so we can accede it
  //   communicate_lx_quad_su3_borders(in);
    
  //   //scan the two virtual nodes
  //   NISSA_PARALLEL_LOOP(isrc_lx,0,loc_vol)
  //     for(int mu=0;mu<NDIM;mu++)
  // 	{
  // 	  //catch links needed to scatter signal forward
  // 	  SU3_TO_VIR_SU3(out[virlx_of_loclx[isrc_lx]][NDIM+mu],in[isrc_lx][mu],vrank_of_loclx(isrc_lx));
	  
  // 	  //copy links also where they are needed to scatter the signal backward, if
  // 	  //sites that need them are not in the border (that would mean that computation must be
  // 	  //done in another node
  // 	  int idst_lx=loclx_neighup[isrc_lx][mu];
  // 	  if(idst_lx<loc_vol) SU3_TO_VIR_SU3(out[virlx_of_loclx[idst_lx]][mu],in[isrc_lx][mu],vrank_of_loclx(idst_lx));
  // 	}
    
  //   //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward
  //   for(int mu=0;mu<NDIM;mu++) //border and link direction
  //     if(paral_dir[mu])
  // 	NISSA_PARALLEL_LOOP(ibord,loc_vol+bord_offset[mu],loc_vol+bord_offset[mu]+bord_dir_vol[mu])
  // 	  {
  // 	    int idst_lx=loclx_neighup[ibord][mu];
  // 	    int vn_dst_virlx=vrank_of_loclx(idst_lx);
  // 	    int idst_virlx=virlx_of_loclx[idst_lx];
	    
  // 	    SU3_TO_VIR_SU3(out[idst_virlx][mu],in[ibord][mu],vn_dst_virlx);
  // 	  }
    
  //   set_borders_invalid(out);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // /*
  //   similar, but the various links are in different blocks
  // */
  // THREADABLE_FUNCTION_2ARG(lx_conf_remap_to_virlx_blocked, vir_su3*,out, quad_su3*,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //communicate conf border so we can accede it
  //   communicate_lx_quad_su3_borders(in);
    
  //   //scan the two virtual nodes
  //   NISSA_PARALLEL_LOOP(isrc_lx,0,loc_vol)
  //     for(int mu=0;mu<NDIM;mu++)
  // 	{
  // 	  //catch links needed to scatter signal forward
  // 	  SU3_TO_VIR_SU3(out[(NDIM+mu)*loc_vol/nvranks+virlx_of_loclx[isrc_lx]],in[isrc_lx][mu],vrank_of_loclx(isrc_lx));
	  
  // 	  //copy links also where they are needed to scatter the signal backward, if
  // 	  //sites that need them are not in the border (that would mean that computation must be
  // 	  //done in another node
  // 	  int idst_lx=loclx_neighup[isrc_lx][mu];
  // 	  if(idst_lx<loc_vol) SU3_TO_VIR_SU3_DAG(out[mu*loc_vol/nvranks+virlx_of_loclx[idst_lx]],in[isrc_lx][mu],vrank_of_loclx(idst_lx));
  // 	}
    
  //   //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward
  //   for(int mu=0;mu<NDIM;mu++) //border and link direction
  //     if(paral_dir[mu])
  // 	NISSA_PARALLEL_LOOP(ibord,loc_vol+bord_offset[mu],loc_vol+bord_offset[mu]+bord_dir_vol[mu])
  // 	  {
  // 	    int idst_lx=loclx_neighup[ibord][mu];
  // 	    int vn_dst_virlx=vrank_of_loclx(idst_lx);
  // 	    int idst_virlx=mu*loc_vol/nvranks+virlx_of_loclx[idst_lx];
	    
  // 	    SU3_TO_VIR_SU3_DAG(out[idst_virlx],in[ibord][mu],vn_dst_virlx);
  // 	  }
    
  //   set_borders_invalid(out);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // THREADABLE_FUNCTION_2ARG(virlx_conf_remap_to_lx, quad_su3*,ext_out, vir_oct_su3*,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //buffer if needed
  //   int bufferize=((void*)ext_out==(void*)in);
  //   quad_su3 *out=bufferize?nissa_malloc("out",loc_vol,quad_su3):ext_out;
    
  //   //split to the two VN
  //   NISSA_PARALLEL_LOOP(ivol_virlx,0,loc_vol/nvranks)
  //     for(int mu=0;mu<NDIM;mu++)
  // 	VIR_SU3_TO_SU3(out[loclx_of_virlx[ivol_virlx]][mu],out[loclx_of_virlx[ivol_virlx]+vrank_lx_offset][mu],
  // 		      in[ivol_virlx][4+mu]);
    
  //   //wait filling
  //   set_borders_invalid(out);
    
  //   //unbuffer if needed
  //   if(bufferize)
  //     {
  // 	vector_copy(ext_out,out);
  // 	nissa_free(out);
  //     }
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // //similar for eo
  // THREADABLE_FUNCTION_2ARG(lx_conf_remap_to_vireo, vir_oct_su3**,out, quad_su3*,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //communicate conf border so we can accede it
  //   communicate_lx_quad_su3_borders(in);
    
  //   //scan the two virtual nodes
  //   NISSA_PARALLEL_LOOP(isrc_lx,0,loc_vol)
  //     for(int mu=0;mu<4;mu++)
  // 	{
  // 	  //catch links needed to scatter signal forward
  // 	  SU3_TO_VIR_SU3(out[loclx_parity[isrc_lx]][vireo_of_loclx[isrc_lx]][4+mu],in[isrc_lx][mu],vrank_of_loclx(isrc_lx));
	  
  // 	  //copy links also where they are needed to scatter the signal backward, if 
  // 	  //sites that need them are not in the border (that would mean that computation must be 
  // 	  //done in another node)
  // 	  int idst_lx=loclx_neighup[isrc_lx][mu],vn=vrank_of_loclx(idst_lx);
  // 	  if(idst_lx<loc_vol) SU3_TO_VIR_SU3(out[loclx_parity[idst_lx]][vireo_of_loclx[idst_lx]][mu],in[isrc_lx][mu],vn);
  // 	}
    
  //   //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward 
  //   for(int mu=0;mu<4;mu++) //border and link direction
  //     if(paral_dir[mu])
  // 	NISSA_PARALLEL_LOOP(ibord,loc_vol+bord_offset[mu],loc_vol+bord_offset[mu]+bord_dir_vol[mu])
  // 	  {
  // 	    int idst_lx=loclx_neighup[ibord][mu],vn=vrank_of_loclx(idst_lx);
  // 	    SU3_TO_VIR_SU3(out[loclx_parity[idst_lx]][vireo_of_loclx[idst_lx]][mu],in[ibord][mu],vn);
  // 	  }
    
  //   for(int eo=0;eo<2;eo++) set_borders_invalid(out[eo]);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // //similar for eo
  // THREADABLE_FUNCTION_2ARG(lx_conf_remap_to_single_vireo, vir_single_oct_su3**,out, quad_su3*,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //communicate conf border so we can accede it
  //   communicate_lx_quad_su3_borders(in);
    
  //   //scan the two virtual nodes
  //   NISSA_PARALLEL_LOOP(isrc_lx,0,loc_vol)
  //     for(int mu=0;mu<4;mu++)
  // 	{
  // 	  //catch links needed to scatter signal forward
  // 	  SU3_TO_VIR_SINGLE_SU3(out[loclx_parity[isrc_lx]][vireo_of_loclx[isrc_lx]][4+mu],in[isrc_lx][mu],vrank_of_loclx(isrc_lx));
	  
  // 	  //copy links also where they are needed to scatter the signal backward, if
  // 	  //sites that need them are not in the border (that would mean that computation must be
  // 	  //done in another node)
  // 	  int idst_lx=loclx_neighup[isrc_lx][mu],vn=vrank_of_loclx(idst_lx);
  // 	  if(idst_lx<loc_vol) SU3_TO_VIR_SINGLE_SU3(out[loclx_parity[idst_lx]][vireo_of_loclx[idst_lx]][mu],in[isrc_lx][mu],vn);
  // 	}
    
  //   //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward
  //   for(int mu=0;mu<4;mu++) //border and link direction
  //     if(paral_dir[mu])
  // 	NISSA_PARALLEL_LOOP(ibord,loc_vol+bord_offset[mu],loc_vol+bord_offset[mu]+bord_dir_vol[mu])
  // 	  {
  // 	    int idst_lx=loclx_neighup[ibord][mu],vn=vrank_of_loclx(idst_lx);
  // 	    SU3_TO_VIR_SINGLE_SU3(out[loclx_parity[idst_lx]][vireo_of_loclx[idst_lx]][mu],in[ibord][mu],vn);
  // 	  }
    
  //   for(int eo=0;eo<2;eo++) set_borders_invalid(out[eo]);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // //similar for evn or odd
  // THREADABLE_FUNCTION_2ARG(eo_conf_remap_to_vireo, vir_oct_su3**,out, quad_su3**,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //communicate conf border so we can accede it
  //   communicate_ev_and_od_quad_su3_borders(in);
    
  //   //scan the two virtual nodes
  //   for(int par=0;par<2;par++)
  //     NISSA_PARALLEL_LOOP(isrc_eo,0,loc_volh)
  // 	for(int mu=0;mu<4;mu++)
  // 	  {
  // 	    //catch links needed to scatter signal forward
  // 	    int vn=vrank_of_loceo(par,isrc_eo);
  // 	    SU3_TO_VIR_SU3(out[par][vireo_of_loceo[par][isrc_eo]][4+mu],in[par][isrc_eo][mu],vn);
	    
  // 	    //copy links also where they are needed to scatter the signal backward, if
  // 	    //sites that need them are not in the border (that would mean that computation must be
  // 	    //done in another node)
  // 	    int idst_eo=loceo_neighup[par][isrc_eo][mu];
  // 	    vn=vrank_of_loceo(!par,idst_eo);
  // 	    if(idst_eo<loc_volh) SU3_TO_VIR_SU3(out[!par][vireo_of_loceo[!par][idst_eo]][mu],in[par][isrc_eo][mu],vn);
  // 	  }
        
  //   //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward
  //   for(int par=0;par<2;par++)
  //     for(int mu=0;mu<4;mu++) //border and link direction
  // 	if(paral_dir[mu])
  // 	  NISSA_PARALLEL_LOOP(ibord,(loc_vol+bord_offset[mu])/2,(loc_vol+bord_offset[mu]+bord_dir_vol[mu])/2)
  // 	    {
  // 	      int idst_eo=loceo_neighup[par][ibord][mu],vn=vrank_of_loceo(!par,idst_eo);
  // 	      SU3_TO_VIR_SU3(out[!par][vireo_of_loceo[!par][idst_eo]][mu],in[par][ibord][mu],vn);
  // 	    }
    
  //   set_borders_invalid(out[EVN]);
  //   set_borders_invalid(out[ODD]);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // //similar for evn or odd
  // THREADABLE_FUNCTION_2ARG(eo_conf_remap_to_single_vireo, vir_single_oct_su3**,out, quad_su3**,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //communicate conf border so we can accede it
  //   communicate_ev_and_od_quad_su3_borders(in);
    
  //   //scan the two virtual nodes
  //   for(int par=0;par<2;par++)
  //     NISSA_PARALLEL_LOOP(isrc_eo,0,loc_volh)
  // 	for(int mu=0;mu<4;mu++)
  // 	  {
  // 	    //catch links needed to scatter signal forward
  // 	    int vn=vrank_of_loceo(par,isrc_eo);
  // 	    SU3_TO_VIR_SINGLE_SU3(out[par][vireo_of_loceo[par][isrc_eo]][4+mu],in[par][isrc_eo][mu],vn);
	    
  // 	    //copy links also where they are needed to scatter the signal backward, if
  // 	    //sites that need them are not in the border (that would mean that computation must be
  // 	    //done in another node)
  // 	    int idst_eo=loceo_neighup[par][isrc_eo][mu];
  // 	    vn=vrank_of_loceo(!par,idst_eo);
  // 	    if(idst_eo<loc_volh) SU3_TO_VIR_SINGLE_SU3(out[!par][vireo_of_loceo[!par][idst_eo]][mu],
  // 						      in[par][isrc_eo][mu],vn);
  // 	  }
    
  //   //scan the backward borders (first half of lx border) to finish catching links needed to scatter signal backward
  //   for(int par=0;par<2;par++)
  //     for(int mu=0;mu<4;mu++) //border and link direction
  // 	if(paral_dir[mu])
  // 	  NISSA_PARALLEL_LOOP(ibord,(loc_vol+bord_offset[mu])/2,(loc_vol+bord_offset[mu]+bord_dir_vol[mu])/2)
  // 	    {
  // 	      int idst_eo=loceo_neighup[par][ibord][mu],vn=vrank_of_loceo(!par,idst_eo);
  // 	      SU3_TO_VIR_SINGLE_SU3(out[!par][vireo_of_loceo[!par][idst_eo]][mu],in[par][ibord][mu],vn);
  // 	    }
    
  //   set_borders_invalid(out[EVN]);
  //   set_borders_invalid(out[ODD]);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  //remap an lx vector to vir[some] layout
  THREADABLE_FUNCTION_8ARG(something_remap_to_virsome_internal, void*,out, void*,in, int,vol, int,size_per_site, int,nel_per_site, int,nvranks, int*,vrank_of_locsite, int*,idx_out)
  {
    GET_THREAD_ID();
    START_TIMING(remap_time,nremap);
    
    if(out==in) CRASH("cannot use with out==in");
    
    //copy the various virtual ranks
    NISSA_PARALLEL_LOOP(isite,0,vol)
      for(int iel=0;iel<nel_per_site;iel++)
	memcpy((char*)out+size_per_site*(vrank_of_locsite[isite]+nvranks*(iel+nel_per_site*idx_out[isite]))
	       ,
	       (char*)in+size_per_site*(iel+nel_per_site*isite)
	       ,
	       size_per_site);
    
    //wait filling
    set_borders_invalid(out);
    
    STOP_TIMING(remap_time);
  }
  THREADABLE_FUNCTION_END
  //reverse
  THREADABLE_FUNCTION_9ARG(virsome_remap_to_something_internal, void*,out, void*,in, int,vol, int,size_per_site, int,nel_per_site, int,nvranks, int*,vrank_of_locsite, int*,vrank_locsite_offset, int*,idx_out)
  {
    GET_THREAD_ID();
    START_TIMING(remap_time,nremap);
    
    if(out==in) CRASH("cannot use with out==in");
    
    //split the virtual ranks
    NISSA_PARALLEL_LOOP(virsome,0,vol/nvranks)
      for(int iel=0;iel<nel_per_site;iel++)
  	for(int vrank=0;vrank<nvranks;vrank++)
	  memcpy((char*)out+size_per_site*(iel+nel_per_site*(idx_out[virsome]+vrank_locsite_offset[vrank]))
		 ,
		 (char*)in+size_per_site*(vrank+nvranks*(iel+nel_per_site*virsome))
		 ,
		 size_per_site);
    
    //wait filling
    set_borders_invalid(out);
    
    STOP_TIMING(remap_time);
  }
  THREADABLE_FUNCTION_END
  
  // //only even or odd
  // THREADABLE_FUNCTION_3ARG(evn_or_odd_spincolor_remap_to_virevn_or_odd, vir_spincolor*,out, spincolor*,in, int,par)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //split to the two VN
  //   NISSA_PARALLEL_LOOP(ivol_eo,0,loc_volh)
  //     SPINCOLOR_TO_VIR_SPINCOLOR(out[vireo_of_loceo[par][ivol_eo]],in[ivol_eo],vrank_of_loceo(EVN,ivol_eo));
    
  //   //wait filling
  //   set_borders_invalid(out);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  // //reverse
  // THREADABLE_FUNCTION_3ARG(virevn_or_odd_spincolor_remap_to_evn_or_odd, spincolor*,out, vir_spincolor*,in, int,par)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //split to the two VN
  //   NISSA_PARALLEL_LOOP(ivol_vireo,0,loc_volh/nvranks)
  //     VIR_SPINCOLOR_TO_SPINCOLOR(out[loceo_of_vireo[par][ivol_vireo]],out[loceo_of_vireo[par][ivol_vireo]+vrank_eo_offset],
  // 			in[ivol_vireo]);
    
  //   //wait filling
  //   set_borders_invalid(out);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // //remap a spincolor_128 from lx to virlx layout
  // THREADABLE_FUNCTION_2ARG(lx_spincolor_128_remap_to_virlx, vir_spincolor_128*,ext_out, spincolor_128*,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //bufferize if needed
  //   int bufferize=(void*)ext_out==(void*)in;
  //   vir_spincolor_128 *out=bufferize?nissa_malloc("out",loc_vol/nvranks,vir_spincolor_128):ext_out;
    
  //   //copy the various VN
  //   NISSA_PARALLEL_LOOP(ivol_lx,0,loc_vol)
  //     SPINCOLOR_128_TO_VIR_SPINCOLOR_128(out[virlx_of_loclx[ivol_lx]],in[ivol_lx],vrank_of_loclx(ivol_lx));
    
  //   //wait filling
  //   set_borders_invalid(out);
    
  //   //unbuffer if needed
  //   if(bufferize)
  //     {
  // 	vector_copy(ext_out,out);
  // 	nissa_free(out);
  //     }
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  // //reverse
  // THREADABLE_FUNCTION_2ARG(virlx_spincolor_128_remap_to_lx, spincolor_128*,ext_out, vir_spincolor_128*,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //buffer if needed
  //   int bufferize=(void*)ext_out==(void*)in;
  //   spincolor_128 *out=bufferize?nissa_malloc("out",loc_vol,spincolor_128):ext_out;
    
  //   //split to the two VN
  //   NISSA_PARALLEL_LOOP(ivol_virlx,0,loc_vol/nvranks)
  //     VIR_SPINCOLOR_128_TO_SPINCOLOR_128(out[loclx_of_virlx[ivol_virlx]],out[loclx_of_virlx[ivol_virlx]+vrank_lx_offset],
  // 					in[ivol_virlx]);
    
  //   //wait filling
  //   set_borders_invalid(out);
    
  //   //unbuffer if needed
  //   if(bufferize)
  //     {
  // 	vector_copy(ext_out,out);
  // 	nissa_free(out);
  //     }
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // //remap a spincolor from lx to vireo layout
  // THREADABLE_FUNCTION_2ARG(lx_spincolor_remap_to_vireo, vir_spincolor**,out, spincolor*,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //copy the various VN
  //   NISSA_PARALLEL_LOOP(ivol_lx,0,loc_vol)
  //     SPINCOLOR_TO_VIR_SPINCOLOR(out[loclx_parity[ivol_lx]][vireo_of_loclx[ivol_lx]],in[ivol_lx],vrank_of_loclx(ivol_lx));
    
  //   //wait filling
  //   set_borders_invalid(out[EVN]);
  //   set_borders_invalid(out[ODD]);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  // //reverse
  // THREADABLE_FUNCTION_2ARG(vireo_spincolor_remap_to_lx, spincolor*,out, vir_spincolor**,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //split to the two VN
  //   for(int par=0;par<2;par++)
  //     NISSA_PARALLEL_LOOP(ivol_vireo,0,loc_volh/nvranks)
  // 	VIR_SPINCOLOR_TO_SPINCOLOR(out[loclx_of_vireo[par][ivol_vireo]],
  // 				  out[loclx_of_vireo[par][ivol_vireo]+vrank_lx_offset],
  // 				  in[par][ivol_vireo]);
    
  //   //wait filling
  //   set_borders_invalid(out);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // ////////////////////////////////////////////////////////////////////
  
  // //quad_su3 to vir_quad_su3
  // THREADABLE_FUNCTION_2ARG(lx_quad_su3_remap_to_virlx, vir_quad_su3*,out, quad_su3*,in)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   NISSA_PARALLEL_LOOP(ivol_lx,0,loc_vol)
  //     for(int mu=0;mu<NDIM;mu++)
  // 	SU3_TO_VIR_SU3(out[virlx_of_loclx[ivol_lx]][mu],in[ivol_lx][mu],vrank_of_loclx(ivol_lx));
    
  //   set_borders_invalid(out);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // //only even or odd
  // THREADABLE_FUNCTION_3ARG(evn_or_odd_quad_su3_remap_to_virevn_or_odd, vir_quad_su3*,out, quad_su3*,in, int,par)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //split to the two VN
  //   NISSA_PARALLEL_LOOP(ivol_eo,0,loc_volh)
  //     for(int mu=0;mu<NDIM;mu++)
  // 	SU3_TO_VIR_SU3(out[vireo_of_loceo[par][ivol_eo]][mu],in[ivol_eo][mu],vrank_of_loceo(EVN,ivol_eo));
    
  //   //wait filling
  //   set_borders_invalid(out);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  // //reverse
  // THREADABLE_FUNCTION_3ARG(virevn_or_odd_quad_su3_remap_to_evn_or_odd, quad_su3*,out, vir_quad_su3*,in, int,par)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //split to the two VN
  //   NISSA_PARALLEL_LOOP(ivol_vireo,0,loc_volh/nvranks)
  //     for(int mu=0;mu<NDIM;mu++)
  // 	VIR_SU3_TO_SU3(out[loceo_of_vireo[par][ivol_vireo]][mu],out[loceo_of_vireo[par][ivol_vireo]+vrank_eo_offset][mu],
  // 			   in[ivol_vireo][mu]);
    
  //   //wait filling
  //   set_borders_invalid(out);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  // ///////////////// GENERAL COMPLEX TO BE CHECKED //////////////////
  
  // //only even or odd
  // THREADABLE_FUNCTION_4ARG(evn_or_odd_complex_vect_remap_to_virevn_or_odd, vir_complex*,out, complex*,in, int,par, int,vl)
  // {
  //   GET_THREAD_ID();
  //   START_TIMING(remap_time,nremap);
    
  //   //split to the two VN
  //   NISSA_PARALLEL_LOOP(ivol_eo,0,loc_volh)
  //     for(int i=0;i<vl;i++)
  // 	COMPLEX_TO_VIR_COMPLEX(((vir_complex*)out)[i+vl*vireo_of_loceo[par][ivol_eo]],
  // 			       ((complex*)in)[i+vl*ivol_eo],vrank_of_loceo(EVN,ivol_eo));
    
  //   //wait filling
  //   set_borders_invalid(out);
  //   STOP_TIMING(remap_time);
  // }
  // THREADABLE_FUNCTION_END
  
  //! assign vlx to all loclx
  template <class T> void fill_vlx_index(vranks_ord_t<T> *vgeo)
  {
    vranks_grid_t<T> *vg=vgeo->vg;
    for(int loclx=0;loclx<loc_vol;loclx++)
      {
	coords vloc_coords; //< coordinates inside the virtual rank
	for(int mu=0;mu<NDIM;mu++) vloc_coords[mu]=loc_coord_of_loclx[loclx][mu]-vg->vrank_coord[vg->vrank_of_loclx[loclx]][mu]*vg->vloc_size[mu];
	vgeo->vloc_of_loclx[loclx]=lx_of_coord(vloc_coords,vg->vloc_size);
      }
  }
  
  //! assign vsf to all loclx
  template <class T> void fill_vsf_index(vranks_ord_t<T> *vgeo_sf,vranks_ord_t<T> *vgeo_lx)
  {
    vranks_grid_t<T> *vg=vgeo_lx->vg;
    NISSA_LOC_VOL_LOOP(loclx) vgeo_sf->vloc_of_loclx[loclx]=-1;
    int vsf=0;
    for(int is_on_surf_loop=1;is_on_surf_loop>=0;is_on_surf_loop--)
      for(int vloclx=0;vloclx<vg->vloc_vol;vloclx++)
	{
	  int base_loclx=vgeo_lx->loclx_of_vloc[vloclx];
	  bool is_on_surf=false;
	  for(int mu=0;mu<NDIM;mu++) //either 0 or vloc_size[mu]-1
	    is_on_surf|=((loc_coord_of_loclx[base_loclx][mu]%vg->vloc_size[mu])%(vg->vloc_size[mu]-1))==0;
	  //if it is part of current loop, mark it and sign itbasetype::
	  if(is_on_surf==is_on_surf_loop && vgeo_sf->vloc_of_loclx[base_loclx]==-1)
	    {
	      for(int vrank=0;vrank<vg->nvranks;vrank++) vgeo_sf->vloc_of_loclx[base_loclx+vg->vrank_loclx_offset[vrank]]=vsf;
	      vsf++;
	    }
	}
	
	//checks
	if(vsf!=vg->vloc_vol) CRASH("vsuflx arrived at %d, expecting %d",vsf,vg->vloc_vol);
  }
  
  //set virtual geometry
  void set_vranks_geometry()
  {
    vlx_double_geom.init(vdouble_grid,fill_vlx_index<double>);
    vlx_float_geom.init(vfloat_grid,fill_vlx_index<float>);
    vsf_double_geom.init(vdouble_grid,std::bind(fill_vsf_index<double>,std::placeholders::_1,&vlx_double_geom));
    vsf_float_geom.init(vfloat_grid,std::bind(fill_vsf_index<float>,std::placeholders::_1,&vlx_float_geom));
  }
  
  //unset it
  void unset_vranks_geometry()
  {
    //destroy the geometries
    vlx_double_geom.destroy();
    vlx_float_geom.destroy();
    vsf_double_geom.destroy();
    vsf_float_geom.destroy();
    //destroy also the grids
    vdouble_grid.destroy();
    vfloat_grid.destroy();
  }
}
