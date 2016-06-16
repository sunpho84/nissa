#include "nissa.hpp"

#include <math.h>

#include "bgq/spi.hpp"

#include "bgq/staggered_hopping_matrix_eo_or_oe_bgq.hpp"
#include "bgq/Wilson_hopping_matrix_lx_bgq.hpp"

#if (defined BGQ && !defined BGQ_EMU)
 #include <bgpm/include/bgpm.h>
 #include <firmware/include/personality.h>
#endif

using namespace nissa;

const int nbench=10,nbench_port=10;

int seed=100;
double mu=0.03,kappa=0.137;
coords is_closed;

#ifndef SPI
  coords_5D spi_dir_is_torus;
#endif

int comp_nhop(coords_5D n,coords_5D tcoords,coords_5D tdims,bool p=0)
{
  int nhop=0;
  for(int idir=0;idir<5;idir++)
    {
      if(p) master_printf(" net_dir %d, loc_coord %d, neigh_coord %d, is_torus %d\n",
			  idir,tcoords[idir],n[idir],spi_dir_is_torus[idir]);
      
      int off=abs(n[idir]-tcoords[idir]);
      if(spi_dir_is_torus[idir]) off=std::min(off,std::min(tdims[idir]+n[idir]-tcoords[idir],
							   tdims[idir]+tcoords[idir]-n[idir]));
      nhop+=off;
    }
  
  return nhop;
}

void check_torus()
{
#ifdef SPI
  //get the personality
  Personality_t pers;
  Kernel_GetPersonality(&pers,sizeof(pers));
      
  for(int idir=0;idir<5;idir++) spi_dir_is_torus[idir]=ND_GET_TORUS(idir,pers.Network_Config.NetFlags);
  
  coords_5D tcoords={pers.Network_Config.Acoord,pers.Network_Config.Bcoord,pers.Network_Config.Ccoord,
		     pers.Network_Config.Dcoord,pers.Network_Config.Ecoord};
  
  coords_5D tdims={pers.Network_Config.Anodes,pers.Network_Config.Bnodes,pers.Network_Config.Cnodes,
		   pers.Network_Config.Dnodes,pers.Network_Config.Enodes};
  
  //write ordered coordinates
  if(0)
    for(int irank=0;irank<nranks;irank++)
      {
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==irank)
	  {
	    printf("%d   ",irank);
	    for(int mu=0;mu<4;mu++) printf("%d ",rank_coord[mu]);
	    printf("   ");
	    for(int mu=0;mu<5;mu++) printf("%d ",tcoords[mu]);
	    printf("\n");
	  }
	fflush(stdout);
      }
  
  for(int idir=0;idir<5;idir++)
    master_printf("Dir %d dim: %d, is torus: %d\n",idir,tdims[idir],spi_dir_is_torus[idir]);
  
  //copy the coords
  coords_5D spi_rank_coord;
  spi_rank_coord[0]=pers.Network_Config.Acoord;
  spi_rank_coord[1]=pers.Network_Config.Bcoord;
  spi_rank_coord[2]=pers.Network_Config.Ccoord;
  spi_rank_coord[3]=pers.Network_Config.Dcoord;
  spi_rank_coord[4]=pers.Network_Config.Ecoord;

  coords_5D n[4][2];
  for(int mu=0;mu<4;mu++)
    for(int bf=0;bf<2;bf++)
      {
        //send to one dir, receive from the opposite
        MPI_Sendrecv(spi_rank_coord,sizeof(coords_5D),MPI_BYTE,rank_neigh[bf][mu],0,
                     n[mu][bf],sizeof(coords_5D),MPI_BYTE,rank_neigh[!bf][mu],0,
                     cart_comm,MPI_STATUS_IGNORE);
	
	int nhop=comp_nhop(n[mu][bf],tcoords,tdims,1);
	master_printf("Dir %d bf %d nhop: %d\n",mu,bf,nhop);
      }
  
  for(int mu=0;mu<4;mu++)
    {
      is_closed[mu]=0;
      for(int idir=0;idir<5;idir++) is_closed[mu]|=(n[mu][0][idir]!=n[mu][1][idir]);
    }
#endif
}

void time_mpi_timing()
{
  double ave_timing=0,var_timing=0;
  for(int ibench=0;ibench<nbench;ibench++)
    {
      double iter_timing=-take_time();
      iter_timing+=take_time();
      
      ave_timing+=iter_timing;
      var_timing+=iter_timing*iter_timing;
    }
  ave_timing/=nbench;
  var_timing/=nbench;
  var_timing-=ave_timing*ave_timing;
  var_timing=sqrt(var_timing);
  master_printf("timing overherad: %lg +- %lg sec per timing\n",ave_timing,var_timing);
}

//benchmark thread_barrier
THREADABLE_FUNCTION_0ARG(bench_thread_barrier)
{
  int nbench_barr=1;
  double barr_time=-take_time();
  for(int ibench=0;ibench<nbench_barr;ibench++)
    THREAD_BARRIER();
  barr_time+=take_time();
  barr_time/=nbench_barr;
  
  //and pure omp barrier
  double omp_barr_time=-take_time();
  for(int ibench=0;ibench<nbench_barr;ibench++)
    {
#pragma omp barrier
    }
  omp_barr_time+=take_time();
  omp_barr_time/=nbench_barr;
  master_printf("thread barrier: %lg sec/barr, %lg sec/omp_barr\n",barr_time,omp_barr_time);
}}

THREADABLE_FUNCTION_0ARG(bench_scalar_prod)
{
  GET_THREAD_ID();
  spincolor *in=nissa_malloc("in",loc_vol,spincolor);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  in[ivol][id][ic][ri]=glblx_of_loclx[ivol];
  set_borders_invalid(in);
  
  double sc_time=-take_time();
  double sc;
  for(int ibench=0;ibench<nbench;ibench++)
    double_vector_glb_scalar_prod(&sc,(double*)in,(double*)in,loc_vol*sizeof(spincolor)/sizeof(double));
  sc_time+=take_time();
  sc_time/=nbench;
  
  int nflops=2*loc_vol*sizeof(spincolor)/sizeof(double);
  master_printf("time to make a global scalar product (%lg): %lg (%lg Mflops)\n",sc,sc_time,nflops/sc_time*1.e-6);
  
  nissa_free(in);
}}

THREADABLE_FUNCTION_1ARG(bench_vector_copy, int,mode)
{
  GET_THREAD_ID();
  spincolor *in=nissa_malloc("in",loc_vol,spincolor);
  spincolor *in1=nissa_malloc("in1",loc_vol,spincolor);
  spincolor *out=nissa_malloc("out",loc_vol,spincolor);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
        for(int ri=0;ri<2;ri++)
          in[ivol][id][ic][ri]=glblx_of_loclx[ivol];
  set_borders_invalid(in);
  int nbytes=sizeof(spincolor)*loc_vol;
  vector_reset(out);
  vector_copy(in1,in);
  
  double c_time=-take_time();
  switch(mode)
    {
    case 0:
      for(int ibench=0;ibench<nbench;ibench++)
	{
	  int morso=nbytes/nthreads;
	  memcpy((char*)out+thread_id*morso,(char*)in+thread_id*morso,morso);
	  THREAD_BARRIER();
	}
      break;
    case 1:
#ifndef BGQ_EMU
      for(int ibench=0;ibench<nbench;ibench++)
	{
	  DECLARE_REG_VIR_SPINCOLOR(reg);
	  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	    {
	      REG_LOAD_VIR_SPINCOLOR(reg,in[ivol]);
	      STORE_REG_VIR_SPINCOLOR(out[ivol],reg);
	    }
	}
#endif
      break;
    case 2:
#ifndef BGQ_EMU
      for(int ibench=0;ibench<nbench;ibench++)
        {
          int morso=loc_vol*sizeof(spincolor)/sizeof(double)/nthreads;
          double *x=(double*)in,*y=(double*)out;
          int is=thread_id*morso,ie=is+morso;
	  
          for(int i=is;i<ie;i+=4) vec_st(vec_ld(0,x+i),0,y+i);
        }
#endif
      break;
    case 3:
#ifndef BGQ_EMU
      for(int ibench=0;ibench<nbench;ibench++)
        {
          int morso=loc_vol*sizeof(spincolor)/sizeof(double)/nthreads;
          double *x=(double*)in,*y=(double*)out,*z=(double*)in1;
          int is=thread_id*morso,ie=is+morso;
	  
          for(int i=is;i<ie;i+=4)
	    {
	      vector4double reg;
	      vec_st(vec_madd(vec_ld(0,x+i),vec_splats(1.0),vec_mul(vec_ld(0,z+i),vec_splats(1.0))),0,y+i);
	    }
        }
#endif
      break;
    default:
      crash("unknown method");
    }
  c_time+=take_time();
  c_time/=nbench;
  
  master_printf("time to copy %d bytes: %lg, %lg Mbyte/s\n",nbytes,c_time,nbytes/c_time*1.e-6);
		  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
        for(int ri=0;ri<2;ri++)
          if(in[ivol][id][ic][ri]*(1+(mode==3))!=out[ivol][id][ic][ri])
	    crash("%d %d %d %d expcted %lg obtained %lg",ivol,id,ic,ri,in[ivol][id][ic][ri],out[ivol][id][ic][ri]);
  
  nissa_free(in1);
  nissa_free(in);
  nissa_free(out);
}}

THREADABLE_FUNCTION_0ARG(bench_linear_comb)
{
  GET_THREAD_ID();
  spincolor *in1=nissa_malloc("in",loc_vol,spincolor);
  spincolor *in2=nissa_malloc("in",loc_vol,spincolor);
  spincolor *in3=nissa_malloc("in",loc_vol,spincolor);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  in1[ivol][id][ic][ri]=glblx_of_loclx[ivol];
  set_borders_invalid(in1);
  
  vector_copy(in2,in1);
  vector_copy(in3,in1);
  
  double lc_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    double_vector_linear_comb((double*)in3,(double*)in1,1.0,(double*)in2,2.0,loc_vol*sizeof(spincolor)/sizeof(double));
  lc_time+=take_time();
  lc_time/=nbench;
  
  double sc;
  double_vector_glb_scalar_prod(&sc,(double*)in3,(double*)in3,loc_vol*sizeof(spincolor)/sizeof(double));
  
  int nflops=3*loc_vol*sizeof(spincolor)/sizeof(double);
  master_printf("time to make a linear combination (%lg): %lg (%lg Mflops)\n",sc,lc_time,nflops/lc_time*1.e-6);
  
  nissa_free(in1);
  nissa_free(in2);
  nissa_free(in3);
}}

void bench_comm()
{
  double saved_time[8];
  //benchmark communications
  for(int itest=0;itest<32;itest++)
    {
      //set all dir to 0
      int comm_dir[8];
      for(int kdir=0;kdir<8;kdir++) comm_dir[kdir]=0;
      
      //first 8 tests: 8 dirs separately
      if(itest<8) comm_dir[itest]=1;
      
      //test 8: all dirs
      if(itest==8) for(int kdir=0;kdir<8;kdir++) comm_dir[kdir]=1;
      
      //test 9: up and dw 0
      if(itest==9) comm_dir[0]=comm_dir[4]=1;
      
      //test 10: up 1 and 2 dirs
      if(itest==10) comm_dir[1]=comm_dir[2]=1;
      
      //test 11: up 1, 2 and 3 dirs
      if(itest==11) comm_dir[1]=comm_dir[2]=comm_dir[3]=1;
      
      //test 12: up and dw 1, 2 and 3 dirs
      if(itest==12) comm_dir[1]=comm_dir[2]=comm_dir[3]=comm_dir[5]=comm_dir[6]=comm_dir[7]=1;
      
      //test 13: up and dw 1 dir
      if(itest==13) comm_dir[1]=comm_dir[5]=1;
      
      //test 14: up and dw 2 dir
      if(itest==14) comm_dir[2]=comm_dir[6]=1;
      
      //test 15: up and dw 3 dir
      if(itest==15) comm_dir[3]=comm_dir[7]=1;
      
      //test 16: up 0 and 1 dir
      if(itest==16) comm_dir[0]=comm_dir[1]=1;
      
      //test 17: all the up dir
      if(itest==17) comm_dir[0]=comm_dir[1]=comm_dir[2]=comm_dir[3]=1;
      
      //test 18: x+ and y-
      if(itest==18) comm_dir[1]=comm_dir[6]=1;
      
      //test 19: t+ and y+
      if(itest==19) comm_dir[0]=comm_dir[2]=1;
      
      //test 20-23: all up and one dw
      if(itest>=20 && itest<24) comm_dir[0]=comm_dir[1]=comm_dir[2]=comm_dir[3]=comm_dir[itest-16]=1;
      
      //test 24-27: all but one of the up dirs
      if(itest>=24 && itest<28)
	{
	  comm_dir[0]=comm_dir[1]=comm_dir[2]=comm_dir[3]=1;
	  comm_dir[itest-24]=0;
	}
      
      //test 28-31: all but one of the up and corresponding dirs
      if(itest>=28 && itest<32)
	{
	  for(int kdir=0;kdir<8;kdir++) comm_dir[kdir]=1;
	  comm_dir[itest-24]=comm_dir[itest-28]=0;
	}
      
      //compute buff size
      int buff_size=0;
      for(int kdir=0;kdir<8;kdir++) buff_size+=bord_dir_vol[kdir%4]*comm_dir[kdir]*sizeof(spincolor);
      
      //take time of n comms
      double time=-take_time();
      for(int ibench=0;ibench<nbench;ibench++)
	{
	  comm_start(lx_spincolor_comm,comm_dir,buff_size);
	  comm_wait(lx_spincolor_comm);
	}
      time+=take_time();
      
      //print out
      double size=buff_size/1024.0/1024*nbench;
      char tag[8][3]={"t-","x-","y-","z-","t+","x+","y+","z+"};
      master_printf("Rank %d, speed in %d test ( ",rank,itest);
      for(int kdir=0;kdir<8;kdir++) if(comm_dir[kdir]) master_printf("%s ",tag[kdir]);
      master_printf("): %lg Mb / %lg sec = %lg Mb/sec",size,time,size/time);
      if(itest<8) master_printf("\n");
      else
	{
	  double max_time=0;
	  for(int kdir=0;kdir<8;kdir++)
	    if(comm_dir[kdir])
	      {
		max_time=std::max(max_time,saved_time[kdir]);
		if(kdir>=4 && comm_dir[kdir-4] && !is_closed[kdir-4])
		  max_time=std::max(max_time,saved_time[kdir]*2);
	      }
	  master_printf(" [%lg Mb/sec exp]\n",max_time);
	}
      
      
      if(itest<8) saved_time[itest]=time;
    }
}

void debug_apply_stDoe_or_eo(int oe_or_eo)
{
  //create random conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      su3_put_to_id(conf[ivol][mu]);
  set_borders_invalid(conf);
  
  //create random in vector
  color *in=nissa_malloc("in",loc_vol+bord_vol,color);
  vector_reset(in);
  NISSA_LOC_VOL_LOOP(ivol)
    in[ivol][0][0]=glblx_of_loclx[ivol];
  set_borders_invalid(in);
  communicate_lx_color_borders(in);
  
  //remap conf to bgq
  vir_oct_su3 *vir_conf[2]={nissa_malloc("vir_conf_evn",loc_volh/2,vir_oct_su3),
			  nissa_malloc("vir_conf_odd",loc_volh/2,vir_oct_su3)};
  lx_conf_remap_to_vireo(vir_conf,conf);
  
  //remap in to bgq
  vir_color *vir_in_eo[2]={nissa_malloc("vir_in",loc_volh/2,vir_color),
			 nissa_malloc("vir_in",loc_volh/2,vir_color)};
  lx_color_remap_to_vireo(vir_in_eo,in);
  
  {
    //precheck: unmapping
    color *un_out=nissa_malloc("un_out",loc_vol+bord_vol,color);
    vireo_color_remap_to_lx(un_out,vir_in_eo);
    
    //compute average diff
    double diff;
    double_vector_subt((double*)un_out,(double*)un_out,(double*)in,loc_vol*sizeof(color)/sizeof(double));
    double_vector_glb_scalar_prod(&diff,(double*)un_out,(double*)un_out,loc_vol*sizeof(color)/sizeof(double));
    char tag_oe_or_eo[2][3]={"oe","eo"};
    master_printf("Map unmap %s diff: %lg\n",tag_oe_or_eo[oe_or_eo],diff);
    
    nissa_free(un_out);
  }
  
  memset(send_buf,0,send_buf_size);
  memset(recv_buf,0,recv_buf_size);
  
  //apply stag bgq
  apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm(vir_conf,0,vsurf_volh,vir_in_eo[oe_or_eo],oe_or_eo);
  THREAD_BARRIER();
  start_double_staggered_hopping_matrix_oe_or_eo_bgq_communications();
  THREAD_BARRIER();
  
  //compute on the bulk and finish communications
  apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm(vir_conf,vsurf_volh,loc_volh/2,vir_in_eo[oe_or_eo],oe_or_eo);
  THREAD_BARRIER();
  finish_double_staggered_hopping_matrix_oe_or_eo_bgq_communications(oe_or_eo);
  THREAD_BARRIER();
  
  const char EVN_ODD_TAG[2][4]={"EVN","ODD"};
  vir_color *bgq_hopping_matrix_output_data=(vir_color*)send_buf+bord_volh/2;
  for(int ivol=0;ivol<loc_vol;ivol++)
    if(loclx_parity[ivol]==!oe_or_eo)
      {
	int vn=(loc_coord_of_loclx[ivol][vnode_paral_dir]>=loc_size[vnode_paral_dir]/2);
	
	for(int mu=0;mu<4;mu++)
	  {
	    int idw_att=(int)(in[loclx_neighdw[ivol][mu]][0][0]);
	    int mu_att=4+mu;
	    int idw_ott=bgq_hopping_matrix_output_data[vireo_of_loclx[ivol]*8+4+mu][0][vn][0];
	    int mu_ott=-bgq_hopping_matrix_output_data[vireo_of_loclx[ivol]*8+4+mu][0][vn][1]/idw_ott;
	    coords catt,cott;
	    glb_coord_of_glblx(catt,idw_att);
	    glb_coord_of_glblx(cott,idw_ott);
	    if(idw_att!=idw_ott)
	      {
		char tag[2][10]={"corre","WRONG"};
		master_printf("%s_st ivol %d (%d %d %d %d) mu %d dw att %d %d (%d %d %d %d) ott %d (%s) %d (%d %d %d %d)\n",
			      tag[(idw_att!=idw_ott)],ivol,
			      loc_coord_of_loclx[ivol][0],loc_coord_of_loclx[ivol][1],
			      loc_coord_of_loclx[ivol][2],loc_coord_of_loclx[ivol][3],
			      mu,
			      idw_att,mu_att,catt[0],catt[1],catt[2],catt[3],
			      idw_ott,EVN_ODD_TAG[loclx_parity[idw_ott]],mu_ott,cott[0],cott[1],cott[2],cott[3]);
	    }
	}
    }
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    if(loclx_parity[ivol]==!oe_or_eo)
      {
	int vn=(loc_coord_of_loclx[ivol][vnode_paral_dir]>=loc_size[vnode_paral_dir]/2);
	
	for(int mu=0;mu<4;mu++)
	  {
	    int idw_att=(int)(in[loclx_neighup[ivol][mu]][0][0]);
	    int mu_att=mu;
	    int idw_ott=bgq_hopping_matrix_output_data[vireo_of_loclx[ivol]*8+mu][0][vn][0];
	    int mu_ott=+bgq_hopping_matrix_output_data[vireo_of_loclx[ivol]*8+mu][0][vn][1]/idw_ott; //+because dag
	    coords catt,cott;
	    glb_coord_of_glblx(catt,idw_att);
	    glb_coord_of_glblx(cott,idw_ott);
	    if(idw_att!=idw_ott)
	      {
		char tag[2][10]={"corre","WRONG"};
		master_printf("%s_st ivol %d (%d %d %d %d) mu %d fw att %d %d (%d %d %d %d) ott %d (%s) %d (%d %d %d %d)\n",
			      tag[(idw_att!=idw_ott)],ivol,
			      loc_coord_of_loclx[ivol][0],loc_coord_of_loclx[ivol][1],
			      loc_coord_of_loclx[ivol][2],loc_coord_of_loclx[ivol][3],
			      mu,
			      idw_att,mu_att,catt[0],catt[1],catt[2],catt[3],
			      idw_ott,EVN_ODD_TAG[loclx_parity[idw_ott]],mu_ott,cott[0],cott[1],cott[2],cott[3]);
	      }
	  }
      }
  
  nissa_free(in);
  nissa_free(conf);
  nissa_free(vir_in_eo[EVN]);
  nissa_free(vir_in_eo[ODD]);
  nissa_free(vir_conf[EVN]);
  nissa_free(vir_conf[ODD]);
}

void debug_apply_stDeo()
{debug_apply_stDoe_or_eo(1);}

void debug_apply_stDoe()
{debug_apply_stDoe_or_eo(0);}

void debug_apply_tmQ()
{
  //create random conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      su3_put_to_id(conf[ivol][mu]);
  set_borders_invalid(conf);
  
  //create random in vector
  spincolor *in=nissa_malloc("in",loc_vol+bord_vol,spincolor);
  vector_reset(in);
  NISSA_LOC_VOL_LOOP(ivol)
    in[ivol][0][0][0]=glblx_of_loclx[ivol];
  set_borders_invalid(in);
  communicate_lx_spincolor_borders(in);
  
  //remap conf to bgq
  vir_oct_su3 *vir_conf=nissa_malloc("vir_conf",loc_vol+bord_vol,vir_oct_su3);
  lx_conf_remap_to_virlx(vir_conf,conf);
  
  //remap in to bgq
  vir_spincolor *vir_in=nissa_malloc("vir_in",loc_vol+bord_vol,vir_spincolor);
  lx_spincolor_remap_to_virlx(vir_in,in);
  
  {
    //precheck: unmapping
    spincolor *un_out=nissa_malloc("un_out",loc_vol+bord_vol,spincolor);
    virlx_spincolor_remap_to_lx(un_out,vir_in);
    //compute average diff
    double diff;
    double_vector_subt((double*)un_out,(double*)un_out,(double*)in,loc_vol*sizeof(spincolor)/sizeof(double));
    double_vector_glb_scalar_prod(&diff,(double*)un_out,(double*)un_out,loc_vol*sizeof(spincolor)/sizeof(double));
    master_printf("Map unmap lx diff: %lg\n",diff);
    
    nissa_free(un_out);
  }
  
  memset(send_buf,0,send_buf_size);
  memset(recv_buf,0,recv_buf_size);
  
  //apply tm bgq
  master_printf("Applying on the vsurface: %d-%d\n",0,vsurf_vol);
  apply_Wilson_hopping_matrix_lx_bgq_nocomm(vir_conf,0,vsurf_vol,vir_in);
  THREAD_BARRIER();
  start_Wilson_hopping_matrix_lx_bgq_communications();
  THREAD_BARRIER();
  
  //compute on the bulk and finish communications
  master_printf("Applying on the bulk: %d-%d\n",vsurf_vol,loc_volh);
  apply_Wilson_hopping_matrix_lx_bgq_nocomm(vir_conf,vsurf_vol,loc_volh,vir_in);
  THREAD_BARRIER();
  finish_Wilson_hopping_matrix_lx_bgq_communications();
  THREAD_BARRIER();
  
  vir_halfspincolor *bgq_hopping_matrix_output_data=(vir_halfspincolor*)send_buf+bord_volh;
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      int vn=(loc_coord_of_loclx[ivol][vnode_paral_dir]>=loc_size[vnode_paral_dir]/2);
      
      for(int mu=0;mu<4;mu++)
	{
	  int idw_att=(int)(in[loclx_neighdw[ivol][mu]][0][0][0]);
	  int mu_att=4+mu;
	  int idw_ott=bgq_hopping_matrix_output_data[virlx_of_loclx[ivol]*8+4+mu][0][0][vn][0];
	  int mu_ott=-bgq_hopping_matrix_output_data[virlx_of_loclx[ivol]*8+4+mu][0][0][vn][1]/idw_ott;
	  coords catt,cott;
	  glb_coord_of_glblx(catt,idw_att);
	  glb_coord_of_glblx(cott,idw_ott);
	  if(idw_att!=idw_ott)
	    {
	      char tag[2][10]={"corre","WRONG"};
	      master_printf("%s_tm ivol %d (%d %d %d %d) mu %d dw att %d %d (%d %d %d %d) ott %d %d (%d %d %d %d)\n",
			    tag[(idw_att!=idw_ott)],ivol,
			    loc_coord_of_loclx[ivol][0],loc_coord_of_loclx[ivol][1],
			    loc_coord_of_loclx[ivol][2],loc_coord_of_loclx[ivol][3],
			    mu,
			    idw_att,mu_att,catt[0],catt[1],catt[2],catt[3],
			    idw_ott,mu_ott,cott[0],cott[1],cott[2],cott[3]);
	    }
	}
    }
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      int vn=(loc_coord_of_loclx[ivol][vnode_paral_dir]>=loc_size[vnode_paral_dir]/2);
      
      for(int mu=0;mu<4;mu++)
	{
	  int idw_att=(int)(in[loclx_neighup[ivol][mu]][0][0][0]);
	  int mu_att=mu;
	  int idw_ott=bgq_hopping_matrix_output_data[virlx_of_loclx[ivol]*8+mu][0][0][vn][0];
	  int mu_ott=+bgq_hopping_matrix_output_data[virlx_of_loclx[ivol]*8+mu][0][0][vn][1]/idw_ott; //+because dag
	  coords catt,cott;
	  glb_coord_of_glblx(catt,idw_att);
	  glb_coord_of_glblx(cott,idw_ott);
	  if(idw_att!=idw_ott)
	    {
	      char tag[2][10]={"corre","WRONG"};
	      master_printf("%s_tm ivol %d (%d %d %d %d) mu %d fw att %d %d (%d %d %d %d) ott %d %d (%d %d %d %d)\n",
			    tag[(idw_att!=idw_ott)],ivol,
			    loc_coord_of_loclx[ivol][0],loc_coord_of_loclx[ivol][1],
			    loc_coord_of_loclx[ivol][2],loc_coord_of_loclx[ivol][3],
			    mu,
			    idw_att,mu_att,catt[0],catt[1],catt[2],catt[3],
			    idw_ott,mu_ott,cott[0],cott[1],cott[2],cott[3]);
	    }
	}
    }
  
  nissa_free(in);
  nissa_free(conf);
  nissa_free(vir_in);
  nissa_free(vir_conf);
}

void debug2_tm()
{
  master_printf("normal way\n");
  
  //create random conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      su3_put_to_rnd(conf[ivol][mu],loc_rnd_gen[ivol]);
  set_borders_invalid(conf);
  
  //create random in vector
  spincolor *in=nissa_malloc("in",loc_vol+bord_vol,spincolor);
  generate_undiluted_source(in,RND_Z4,-1);
  
  //apply a fixed number of time
  spincolor *out=nissa_malloc("out",loc_vol+bord_vol,spincolor);
  spincolor *tmp=nissa_malloc("tmp",loc_vol+bord_vol,spincolor);
  double port_time=-take_time();
  for(int ibench=0;ibench<nbench_port;ibench++)
    apply_tmQ2(out,conf,kappa,tmp,mu,in);
  port_time+=take_time();
  port_time/=2*nbench_port;
  
  master_printf("bgq way\n");
  
  //remap conf to bgq
  vir_oct_su3 *vir_conf=nissa_malloc("vir_conf",loc_vol+bord_vol,vir_oct_su3);
  lx_conf_remap_to_virlx(vir_conf,conf);
  
  //remap in to bgq
  vir_spincolor *vir_in=nissa_malloc("vir_in",loc_vol+bord_vol,vir_spincolor);
  lx_spincolor_remap_to_virlx(vir_in,in);
  
  //apply bgq
  vir_spincolor *vir_out=nissa_malloc("vir_out",loc_vol+bord_vol,vir_spincolor);
  double bgq_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    apply_tmQ2_bgq(vir_out,vir_conf,kappa,mu,vir_in);
  bgq_time+=take_time();
  bgq_time/=2*nbench;
  
  //unmap to compare
  spincolor *un_out=nissa_malloc("un_out",loc_vol+bord_vol,spincolor);
  virlx_spincolor_remap_to_lx(un_out,vir_out);
  
  //compute average diff
  double diff;
  double_vector_subt((double*)un_out,(double*)un_out,(double*)out,loc_vol*sizeof(spincolor)/sizeof(double));
  double_vector_glb_scalar_prod(&diff,(double*)un_out,(double*)un_out,loc_vol*sizeof(spincolor)/sizeof(double));
  master_printf("TM application diff: %lg\n",diff);
  
  //print benchmark
  master_printf("In+conf+addr size: %lg Mbytes (L2 cache: 32 Mb)\n",
		(loc_vol*(sizeof(quad_su3)+sizeof(spincolor))+8*loc_volh*sizeof(vir_halfspincolor*))/1024.0/1024.0);
  master_printf("Time to apply %d time:\n",nbench);
  master_printf(" %lg sec in port mode\n",port_time);
  master_printf(" %lg sec in bgq mode\n",bgq_time);
  
  //benchmark pure hopping matrix application
  double hop_bgq_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    apply_Wilson_hopping_matrix_lx_bgq_nocomm(vir_conf,0,loc_volh,vir_in);
  hop_bgq_time+=take_time();
  hop_bgq_time/=nbench;
  int nflops_hop=1152*loc_vol;
  master_printf("hop_bgq_time: %lg sec, %d flops, %lg Mflops\n",hop_bgq_time,nflops_hop,nflops_hop*1e-6/hop_bgq_time);
  
  //benchmark expansion
  double exp_bgq_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    hopping_matrix_lx_expand_to_Q_bgq(vir_out);
  exp_bgq_time+=take_time();
  exp_bgq_time/=nbench;
  int nflops_exp=312*loc_vol;
  master_printf("exp_bgq_time: %lg sec, %d flops, %lg Mflops\n",exp_bgq_time,nflops_exp,nflops_exp*1e-6/exp_bgq_time);
  master_printf("(hop+exp)_bgq_time: %lg sec, %d flops, %lg Mflops\n",
		hop_bgq_time+exp_bgq_time,nflops_exp+nflops_hop,(nflops_exp+nflops_hop)*
		1e-6/(hop_bgq_time+exp_bgq_time));
  
  //benchmark buff filling
  double buff_filling_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    bgq_Wilson_hopping_matrix_lx_vdir_VN_comm_and_buff_fill();
  buff_filling_time+=take_time();
  buff_filling_time/=nbench;
  master_printf("buff_filling_time: %lg sec\n",buff_filling_time);
  
  //benchmark bulk computation
  double bulk_computation_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    apply_Wilson_hopping_matrix_lx_bgq_nocomm(vir_conf,vsurf_vol,loc_volh,vir_in);
  bulk_computation_time+=take_time();
  bulk_computation_time/=nbench;
  master_printf("bulk_computation_time: %lg sec\n",bulk_computation_time);
  
  //benchmark surf computation
  double surf_compuation_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    apply_Wilson_hopping_matrix_lx_bgq_nocomm(vir_conf,0,vsurf_vol,vir_in);
  surf_compuation_time+=take_time();
  surf_compuation_time/=nbench;
  master_printf("surf_compuation_time: %lg sec\n",surf_compuation_time);
  
  //benchmark buff unfilling
  double buff_unfilling_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    finish_Wilson_hopping_matrix_lx_bgq_communications();
  buff_unfilling_time+=take_time();
  buff_unfilling_time/=nbench;
  master_printf("buff_unfilling_time: %lg sec\n",buff_unfilling_time);
  
  //benchmark pure spi communication without data moving
  double pure_spi_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    {
      comm_start(lx_halfspincolor_comm);
      comm_wait(lx_halfspincolor_comm);
    }
  pure_spi_time+=take_time();
  pure_spi_time/=nbench;
  double data_exch=lx_halfspincolor_comm.tot_mess_size/1024.0/1024.0;
  master_printf("pure_comm_time: %lg sec, %lg Mbytes, %lg Mb/sec\n",pure_spi_time,data_exch,data_exch/pure_spi_time);
  
  nissa_free(conf);
  nissa_free(in);
  nissa_free(tmp);
  nissa_free(out);
  
  nissa_free(vir_conf);
  nissa_free(vir_in);
  nissa_free(vir_out);
  nissa_free(un_out);
}

THREADABLE_FUNCTION_0ARG(bull)
{
}}

void debug2_st()
{
  master_printf("------staggered------\n");
  
  double mass2=0.8;
  
  //create random conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  NISSA_LOC_VOL_LOOP(ivol)
    for(int mu=0;mu<4;mu++)
      {
	su3_put_to_rnd(conf[ivol][mu],loc_rnd_gen[ivol]);
	//int fw=loclx_neighup[ivol][mu];
	//if(fw<loc_vol) su3_put_to_diag(conf[ivol][mu],glblx_of_loclx[fw]);
	//else           su3_put_to_diag(conf[ivol][mu],glblx_of_bordlx[fw-loc_vol]);
      }
  set_borders_invalid(conf);
  
  //create random in vector
  color *in=nissa_malloc("in",loc_vol+bord_vol,color);
  vector_reset(in);
  NISSA_LOC_VOL_LOOP(ivol)
    //in[ivol][0][0]=1;
    in[ivol][0][0]=glblx_of_loclx[ivol];
  set_borders_invalid(in);
  communicate_lx_color_borders(in);
  
  ///////////////////////// normal way //////////////////////
  
  master_printf("normal way\n");
  
  //split into eo the conf
  quad_su3 *conf_eo[2]={nissa_malloc("conf",loc_volh+bord_volh,quad_su3),
			nissa_malloc("conf",loc_volh+bord_volh,quad_su3)};
  split_lx_vector_into_eo_parts(conf_eo,conf);
  
  color *in_eo[2]={nissa_malloc("in",loc_volh+bord_volh,color),nissa_malloc("in",loc_volh+bord_volh,color)};
  split_lx_vector_into_eo_parts(in_eo,in);
  
  double port_time=-take_time();
  color *out[2]={nissa_malloc("out_e",loc_volh+bord_volh,color),nissa_malloc("out_o",loc_volh+bord_volh,color)};
  for(int ibench=0;ibench<nbench_port;ibench++)
    apply_stD2ee_m2(out[EVN],conf_eo,in_eo[ODD]/*used as tmp*/,mass2,in_eo[EVN]);
    //apply_st2Doe(out[ODD],conf_eo,in_eo[EVN]);
  port_time+=take_time();
  port_time/=nbench_port;
  
  nissa_free(conf_eo[EVN]);
  nissa_free(conf_eo[ODD]);
  nissa_free(in_eo[EVN]);
  nissa_free(in_eo[ODD]);
  
  /////////////////////////// bgq way //////////////////////
  
  master_printf("bgq way\n");
  
  //remap conf to bgq
  vir_oct_su3 *vir_conf_eo[2]={nissa_malloc("vir_conf_evn",loc_volh/2,vir_oct_su3),
			     nissa_malloc("vir_conf_odd",loc_volh/2,vir_oct_su3)};
  lx_conf_remap_to_vireo(vir_conf_eo,conf);
  
  //remap in to bgq
  vir_color *vir_in_eo[2]={nissa_malloc("vir_in",loc_volh/2,vir_color),
			 nissa_malloc("vir_in",loc_volh/2,vir_color)};
  lx_color_remap_to_vireo(vir_in_eo,in);
  
  //apply bgq
  vir_color *vir_out_eo[2]={nissa_malloc("vir_out",loc_volh/2,vir_color),
			  nissa_malloc("vir_out",loc_volh/2,vir_color)};
  
  vir_single_oct_su3 *vir_single_conf_eo[2]={nissa_malloc("vir_single_conf_evn",loc_volh/2,vir_single_oct_su3),
					   nissa_malloc("vir_single_conf_odd",loc_volh/2,vir_single_oct_su3)};
  vir_single_color *vir_single_in_eo[2]={nissa_malloc("vir_single_in",loc_volh/2,vir_single_color),
				       nissa_malloc("vir_single_in",loc_volh/2,vir_single_color)};
  vir_single_color *vir_single_out_eo[2]={nissa_malloc("vir_single_out",loc_volh/2,vir_single_color),
					nissa_malloc("vir_single_out",loc_volh/2,vir_single_color)};
  
  //remap
  for(int eo=0;eo<2;eo++)
    {
      for(int i=0;i<(int)(sizeof(vir_oct_su3)/sizeof(double)*loc_volh/2);i++)
	((float*)(vir_single_conf_eo[eo]))[i]=((double*)(vir_conf_eo[eo]))[i];
      set_borders_invalid(vir_single_conf_eo[eo]);
      
      for(int i=0;i<(int)(sizeof(vir_color)/sizeof(double)*loc_volh/2);i++)
	((float*)(vir_single_in_eo[eo]))[i]=((double*)(vir_in_eo[eo]))[i];
      set_borders_invalid(vir_single_in_eo[eo]);
    }
  
  double bgq_double_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    apply_stD2ee_m2_bgq(vir_out_eo[EVN],vir_conf_eo,mass2,vir_in_eo[EVN]);
  bgq_double_time+=take_time();
  bgq_double_time/=nbench;
  
  double bgq_single_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    apply_single_stD2ee_m2_bgq(vir_single_out_eo[EVN],vir_single_conf_eo,mass2,vir_single_in_eo[EVN]);
  bgq_single_time+=take_time();
  bgq_single_time/=nbench;
  
  //unamp to compare
  color *un_out=nissa_malloc("un_out",loc_volh,color);
  color *un_single_out=nissa_malloc("un_single_out",loc_volh,color);
  virevn_or_odd_color_remap_to_evn_or_odd(un_out,vir_out_eo[EVN],EVN);
  virevn_or_odd_single_color_remap_to_evn_or_odd(un_single_out,vir_single_out_eo[EVN],EVN);
  
  //compute average diff
  if(rank==0)
  for(int i=0;i<loc_volh;i++)
    {
      int ex=0;
      int lx=loclx_of_loceo[EVN][i];
      for(int mu=0;mu<4;mu++)
	{
	  int fw=loclx_neighup[lx][mu];
	  int bw=loclx_neighdw[lx][mu];
	  int cf=(fw<loc_vol)?glblx_of_loclx[fw]:glblx_of_bordlx[fw-loc_vol];
	  int cb=(bw<loc_vol)?glblx_of_loclx[bw]:glblx_of_bordlx[bw-loc_vol];
	  ex+=cf;
	  ex-=cb;
	  //printf("  %d %d, %d %d\n",lx,mu,cf,cb);
	}
      //for(int ic=0;ic<3;ic++)
      //for(int ri=0;ri<2;ri++)
    }
  double diff2,diff_single2,norm2;
  double_vector_subt((double*)un_out,(double*)un_out,(double*)out[EVN],loc_volh*sizeof(color)/sizeof(double));
  double_vector_subt((double*)un_single_out,(double*)un_single_out,(double*)out[EVN],loc_volh*sizeof(color)/sizeof(double));
  double_vector_glb_scalar_prod(&diff2,(double*)un_out,(double*)un_out,loc_volh*sizeof(color)/sizeof(double));
  double_vector_glb_scalar_prod(&diff_single2,(double*)un_single_out,(double*)un_single_out,loc_volh*sizeof(color)/sizeof(double));
  double_vector_glb_scalar_prod(&norm2,(double*)out[EVN],(double*)out[EVN],loc_volh*sizeof(color)/sizeof(double));
  master_printf("ST application relative diff: %lg\n",sqrt(diff2/norm2));
  master_printf("ST application relative diff single: %lg\n",sqrt(diff_single2/norm2));
  master_printf("Time to apply %d time (on average):\n",nbench);
  master_printf(" %lg sec in port mode, %lg MFlops\n",port_time,1156e-6*loc_volh/port_time);
  master_printf(" %lg sec in bgq double mode, %lg MFlops\n",bgq_double_time,1156e-6*loc_volh/bgq_double_time);
  master_printf(" %lg sec in bgq single mode, %lg MFlops\n",bgq_single_time,1156e-6*loc_volh/bgq_single_time);
  
  //benchmark pure hopping matrix application
  double hop_bgq_time[2];
  int nflops_hop=480*loc_volh;
  for(int eo_or_oe=0;eo_or_oe<2;eo_or_oe++)
    {
      hop_bgq_time[eo_or_oe]=-take_time();
      for(int ibench=0;ibench<nbench;ibench++)
	apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm(vir_conf_eo,0,loc_volh/2,vir_in_eo[EVN],eo_or_oe);
      hop_bgq_time[eo_or_oe]+=take_time();
      hop_bgq_time[eo_or_oe]/=nbench;
      master_printf("hop_bgq_time_%s: %lg sec, %d flops, %lg Mflops\n",
		    (eo_or_oe==0)?"eo":"oe",
		    hop_bgq_time[eo_or_oe],nflops_hop,nflops_hop*1e-6/hop_bgq_time[eo_or_oe]);
    }
  
  //benchmark expansion
  double exp_bgq_time[2][2];
  int nflops_exp=NCOL*(9*flops_per_complex_summ)*loc_volh;
  for(int single_double=0;single_double<2;single_double++)
    for(int eo_or_oe=0;eo_or_oe<2;eo_or_oe++)
      {
	exp_bgq_time[single_double][eo_or_oe]=-take_time();
	if(single_double==0)
	  {
	    if(eo_or_oe==0)
	      for(int ibench=0;ibench<nbench;ibench++)
		hopping_matrix_oe_or_eo_expand_to_single_staggered_D_bgq((vir_single_color*)(vir_out_eo[EVN]));
	    else
	      for(int ibench=0;ibench<nbench;ibench++)
		hopping_matrix_oe_or_eo_expand_to_single_staggered_D_subtract_from_mass2_times_in_bgq((vir_single_color*)(vir_out_eo[ODD]),mass2,(vir_single_color*)(vir_in_eo[EVN]));
	  }
	else
	  {
	    if(eo_or_oe==0)
	      for(int ibench=0;ibench<nbench;ibench++)
		hopping_matrix_oe_or_eo_expand_to_double_staggered_D_bgq(vir_out_eo[EVN]);
	    else
	      for(int ibench=0;ibench<nbench;ibench++)
		hopping_matrix_oe_or_eo_expand_to_double_staggered_D_subtract_from_mass2_times_in_bgq(vir_out_eo[ODD],mass2,vir_in_eo[EVN]);
	  }
	exp_bgq_time[single_double][eo_or_oe]+=take_time();
	exp_bgq_time[single_double][eo_or_oe]/=nbench;
	master_printf("exp_bgq_time_%s (%s): %lg sec, %d flops, %lg Mflops\n",
		      (eo_or_oe==0)?"eo":"oe",(single_double==0?"single":"double"),
		      exp_bgq_time[single_double][eo_or_oe],nflops_exp,nflops_exp*1e-6/exp_bgq_time[single_double][eo_or_oe]);
      }
  
  //total
  double hop_plus_exp_bgq_time=hop_bgq_time[0]+hop_bgq_time[1]+exp_bgq_time[0][0]+exp_bgq_time[0][1];
  int nflops_hop_plus_exp=2*(nflops_exp+nflops_hop);
  master_printf("(hop+exp)_bgq_time: %lg sec, %d flops, %lg Mflops\n",
		hop_plus_exp_bgq_time,nflops_hop_plus_exp,nflops_hop_plus_exp*
		1e-6/hop_plus_exp_bgq_time);

  int ncas=11;
  double tie[ncas];
  for(int cas=0;cas<ncas;cas++) tie[cas]=0;

  for(int ibench=0;ibench<nbench;ibench++)
    for(int cas=0;cas<ncas;cas++)
      {
	const int OE=0;
	const int EO=1;
	tie[cas]-=take_time();
	switch(cas)
	  {
	  case 0:apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm(vir_conf_eo,0,vsurf_volh,vir_in_eo[EVN],OE);break;
	  case 1:start_double_staggered_hopping_matrix_oe_or_eo_bgq_communications();break;
	  case 2:apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm(vir_conf_eo,vsurf_volh,loc_volh/2,vir_in_eo[EVN],OE);break;
	  case 3:finish_double_staggered_hopping_matrix_oe_or_eo_bgq_communications(OE);break;
	  case 4:hopping_matrix_oe_or_eo_expand_to_double_staggered_D_bgq(vir_out_eo[EVN]);break;
	  case 5:apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm(vir_conf_eo,0,vsurf_volh,vir_out_eo[EVN],EO);break;
	  case 6:start_double_staggered_hopping_matrix_oe_or_eo_bgq_communications();break;
	  case 7:apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm(vir_conf_eo,vsurf_volh,loc_volh/2,vir_out_eo[EVN],EO);break;
	  case 8:finish_double_staggered_hopping_matrix_oe_or_eo_bgq_communications(EO);break;
	  case 9:if(mass2!=0) hopping_matrix_oe_or_eo_expand_to_double_staggered_D_subtract_from_mass2_times_in_bgq
				(vir_out_eo[EVN],mass2,vir_in_eo[EVN]);
	    hopping_matrix_oe_or_eo_expand_to_double_staggered_D_bgq(vir_out_eo[EVN],-1);
	    break;
	  default:bull();break;
	  }
	tie[cas]+=take_time();
      }
  
  double tot_tie=0;
  const char name[11][100]={"apply_oe_0_vsurfvolh",
		       "start",
		       "apply_oe_vsurfvolh_loc_volh/2",
		       "finish",
		       "expand",
		       "apply_eo_0_vsurfvolh",
		       "start",
		       "apply_eo_vsurfvolh_loc_volh/2",
		       "finish",
		       "expand",
		       "bench"};

  for(int cas=0;cas<ncas;cas++)
    {
      tie[cas]/=nbench;
      tot_tie+=tie[cas];
      master_printf("Piece %d (%s) time: %lg\n",cas,name[cas],tie[cas]);
    }
  master_printf("Total time: %lg\n",tot_tie);
  
  //benchmark buff filling
  double buff_filling_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    bgq_double_staggered_hopping_matrix_oe_or_eo_vdir_VN_comm_and_buff_fill();
  buff_filling_time+=take_time();
  buff_filling_time/=nbench;
  master_printf("buff_filling_time: %lg sec\n",buff_filling_time);
  
  //benchmark bulk computation
  double bulk_computation_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm(vir_conf_eo,0,vsurf_volh,vir_in_eo[EVN],0);
  bulk_computation_time+=take_time();
  bulk_computation_time/=nbench;
  master_printf("bulk_computation_time: %lg sec\n",bulk_computation_time);
  
  //benchmark surf computation
  double surf_compuation_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    apply_double_staggered_hopping_matrix_oe_or_eo_bgq_nocomm(vir_conf_eo,vsurf_volh,loc_volh/2,vir_in_eo[EVN],0);
  surf_compuation_time+=take_time();
  surf_compuation_time/=nbench;
  master_printf("surf_compuation_time: %lg sec\n",surf_compuation_time);
  
  //benchmark buff unfilling
  double buff_unfilling_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    finish_double_staggered_hopping_matrix_oe_or_eo_bgq_communications(0);
  buff_unfilling_time+=take_time();
  buff_unfilling_time/=nbench;
  master_printf("buff_unfilling_time: %lg sec\n",buff_unfilling_time);
  
  //benchmark pure spi communication without data moving
  for(int single_double=0;single_double<2;single_double++)
    {
      double pure_spi_time=-take_time();
      for(int ibench=0;ibench<nbench;ibench++)
	if(single_double==0)
	  {
	    comm_start(eo_single_color_comm);
	    comm_wait(eo_single_color_comm);
	  }
      else
	{
	  comm_start(eo_color_comm);
	  comm_wait(eo_color_comm);
	}
      
      pure_spi_time+=take_time();
      pure_spi_time/=nbench;
      
      double data_exch=(single_double==0?eo_single_color_comm:eo_color_comm).tot_mess_size/1024.0/1024.0;
      master_printf("pure_comm_time (%s): %lg sec, %lg Mbytes, %lg Mb/sec\n",(single_double==0)?"single":"double",pure_spi_time,data_exch,data_exch/pure_spi_time);
    }
  
  nissa_free(conf);
  nissa_free(vir_conf_eo[0]);
  nissa_free(vir_conf_eo[1]);
  nissa_free(vir_single_conf_eo[0]);
  nissa_free(vir_single_conf_eo[1]);
  nissa_free(in);
  nissa_free(vir_in_eo[0]);
  nissa_free(vir_in_eo[1]);
  nissa_free(vir_single_in_eo[0]);
  nissa_free(vir_single_in_eo[1]);
  nissa_free(out[0]);
  nissa_free(out[1]);
  nissa_free(un_out);
  nissa_free(un_single_out);
  nissa_free(vir_out_eo[0]);
  nissa_free(vir_out_eo[1]);
  nissa_free(vir_single_out_eo[0]);
  nissa_free(vir_single_out_eo[1]);
}

THREADABLE_FUNCTION_0ARG(bench_su3_path_prod)
{
  GET_THREAD_ID();
  
  su3 *conf=nissa_malloc("conf",loc_vol,su3);
  su3 *path_in=nissa_malloc("path_in",loc_vol,su3);
  su3 *path_out=nissa_malloc("path_out",loc_vol,su3);
  
  vector_reset(conf);
  vector_reset(path_in);
  vector_reset(path_out);
  
  double ti=-take_time();
  
  for(int ibench=0;ibench<nbench;ibench++)
    {
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	unsafe_su3_prod_su3(path_out[ivol],conf[ivol],path_in[ivol]);
      THREAD_BARRIER();
    }
  
  ti+=take_time();
  ti/=nbench;
  
  master_printf("Time to take portably a (local) path product: %lg s, %lg MFlops\n",ti,loc_vol*1e-6*flops_per_su3_prod/ti);
  
  nissa_free(path_in);
  nissa_free(path_out);
  nissa_free(conf);
  
  //////////////////////////////////////////////////////////
  
  vir_su3 *vir_conf=nissa_malloc("vir_conf",loc_vol/2,vir_su3);
  vir_su3 *vir_path_in=nissa_malloc("vir_path_in",loc_vol/2,vir_su3);
  vir_su3 *vir_path_out=nissa_malloc("vir_path_out",loc_vol/2,vir_su3);
  
  vector_reset(vir_conf);
  vector_reset(vir_path_in);
  vector_reset(vir_path_out);
  
  ti=-take_time();
  
  for(int ibench=0;ibench<nbench;ibench++)
    {
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol/2)
	{
	  DECLARE_REG_VIR_SU3(REG_VIR_CONF);
	  REG_LOAD_VIR_SU3(REG_VIR_CONF,vir_conf[ivol]);
	  VIR_SU3_PREFETCH_NEXT(vir_conf[ivol]);
	  
	  DECLARE_REG_VIR_SU3(REG_VIR_PATH_IN);
	  REG_LOAD_VIR_SU3(REG_VIR_PATH_IN,vir_path_in[ivol]);
	  VIR_SU3_PREFETCH_NEXT(vir_path_in[ivol]);
	  
	  DECLARE_REG_VIR_SU3(REG_VIR_PATH_OUT);
  	  REG_VIR_SU3_DAG_PROD_VIR_SU3(REG_VIR_PATH_OUT,REG_VIR_CONF,REG_VIR_PATH_IN);
	  STORE_REG_VIR_SU3(vir_path_out[ivol],REG_VIR_PATH_OUT);
	}
      THREAD_BARRIER();
    }
  
  ti+=take_time();
  ti/=nbench;
  
  master_printf("Time to take vectorially a (local) path product with dag: %lg s, %lg MFlops\n",ti,loc_vol*1e-6*flops_per_su3_prod/ti);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ti=-take_time();
  
  for(int ibench=0;ibench<nbench;ibench++)
    {
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol/2)
	{
	  DECLARE_REG_VIR_SU3(REG_VIR_CONF);
	  REG_LOAD_VIR_SU3(REG_VIR_CONF,vir_conf[ivol]);
	  VIR_SU3_PREFETCH_NEXT(vir_conf[ivol]);
	  
	  DECLARE_REG_VIR_SU3(REG_VIR_PATH_IN);
	  REG_LOAD_VIR_SU3(REG_VIR_PATH_IN,vir_path_in[ivol]);
	  VIR_SU3_PREFETCH_NEXT(vir_path_in[ivol]);
	  
	  DECLARE_REG_VIR_SU3(REG_VIR_PATH_OUT);
  	  REG_VIR_SU3_PROD_VIR_SU3(REG_VIR_PATH_OUT,REG_VIR_CONF,REG_VIR_PATH_IN);
	  STORE_REG_VIR_SU3(vir_path_out[ivol],REG_VIR_PATH_OUT);
	}
      THREAD_BARRIER();
    }
  
  ti+=take_time();
  ti/=nbench;
  
  master_printf("Time to take vectorially a (local) path non-daggered product: %lg s, %lg MFlops\n",ti,loc_vol*1e-6*flops_per_su3_prod/ti);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ti=-take_time();
  
  for(int ibench=0;ibench<nbench;ibench++)
    {
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol/2)
	{
	  DECLARE_REG_VIR_PARTIAL_SU3(REG_VIR_CONF);
	  REG_LOAD_VIR_PARTIAL_SU3(REG_VIR_CONF,vir_conf[ivol]);
	  VIR_PARTIAL_SU3_PREFETCH_NEXT(vir_conf[ivol]);
	  
	  DECLARE_REG_VIR_SU3(REG_VIR_PATH_IN);
	  REG_LOAD_VIR_SU3(REG_VIR_PATH_IN,vir_path_in[ivol]);
	  VIR_SU3_PREFETCH_NEXT(vir_path_in[ivol]);
	  
	  DECLARE_REG_VIR_PARTIAL_SU3(REG_VIR_PATH_OUT);
  	  REG_VIR_PARTIAL_SU3_PROD_VIR_SU3(REG_VIR_PATH_OUT,REG_VIR_CONF,REG_VIR_PATH_IN);
	  STORE_REG_VIR_PARTIAL_SU3(vir_path_out[ivol],REG_VIR_PATH_OUT);
	}
      THREAD_BARRIER();
    }
  
  ti+=take_time();
  ti/=nbench;
  
  master_printf("Time to take vectorially a (local) path non-daggered partial product: %lg s, %lg MFlops\n",ti,loc_vol*1e-6*flops_per_su3_prod/ti);
  
  nissa_free(vir_path_in);
  nissa_free(vir_path_out);
  nissa_free(vir_conf);
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  
  for(int mu=0;mu<8;mu++)
    {
      //set all dir to 0 but mu
      int comm_dir[8];
      memset(comm_dir,0,sizeof(int)*8);
      comm_dir[mu]=1;
      
      int buff_size=bord_dir_vol[mu%4]*sizeof(su3);
      
      //take time of n comms
      double time=-take_time();
      for(int ibench=0;ibench<nbench;ibench++)
	{
	  comm_start(lx_su3_comm,comm_dir,buff_size);
	  comm_wait(lx_su3_comm);
	}
      time+=take_time();
      time/=nbench;
      
      //print out
      double size=buff_size/1024.0/1024;
      master_printf("Communications dir %d: %lg Mb / %lg sec = %lg Mb/sec\n",mu,size,time,size/time);
    }
}
THREADABLE_FUNCTION_END

void in_main(int narg,char **arg)
{
  if(narg<3) crash("use %s L T",arg[0]);
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  
  //init
  init_grid(T,L);
  
  start_loc_rnd_gen(seed);
  
  check_torus();
  
  time_mpi_timing();
  bench_thread_barrier();
  
  bench_su3_path_prod();
  
  //master_printf("debugging tmQ\n");
  //debug_apply_tmQ();
  
  //master_printf("debugging stDeo\n");
  //debug_apply_stDeo();
  
  //master_printf("debugging stDoe\n");
  //debug_apply_stDoe();
  
  //prebench
  comm_start(lx_spincolor_comm);
  comm_wait(lx_spincolor_comm);
  
  bench_comm();
  
  bench_scalar_prod();
  bench_linear_comb();
#ifdef BGQ_EMU
  int nmodes=1;
#else
  int nmodes=4;
#endif
  for(int mode=0;mode<nmodes;mode++) bench_vector_copy(mode);
  
  //////////////////////////////////////
  
  debug2_tm();
  debug2_st();
  
  /////////////////////////////////////////////
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
