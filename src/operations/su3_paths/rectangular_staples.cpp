#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "rectangular_staples.hpp"
#include "threads/threads.hpp"

namespace nissa
{
#define COMPUTE_RECT_FW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_prod_su3(TEMP,A,B);		\
  unsafe_su3_prod_su3_dag(OUT,TEMP,C);
#define SUMM_RECT_FW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_prod_su3(TEMP,A,B);		\
  su3_summ_the_prod_su3_dag(OUT,TEMP,C);
  
#define COMPUTE_POINT_RECT_FW_STAPLES(out,conf,sq_staples,A,B,F,imu,mu,inu,nu,temp) \
  COMPUTE_RECT_FW_STAPLE(out[A][mu][3+inu],sq_staples[A][nu][imu],conf[B][mu],conf[F][nu],temp); /*bw sq staple*/ \
  SUMM_RECT_FW_STAPLE(out[A][mu][3+inu],conf[A][nu],sq_staples[B][mu][3+inu],conf[F][nu],temp);  /*fw sq staple*/ \
  SUMM_RECT_FW_STAPLE(out[A][mu][3+inu],conf[A][nu],conf[B][mu],sq_staples[F][nu][3+imu],temp);  /*fw sq staple*/
  
#define COMPUTE_RECT_BW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_dag_prod_su3(TEMP,A,B);		\
  unsafe_su3_prod_su3(OUT,TEMP,C);
#define SUMM_RECT_BW_STAPLE(OUT,A,B,C,TEMP)	\
  unsafe_su3_dag_prod_su3(TEMP,A,B);		\
  su3_summ_the_prod_su3(OUT,TEMP,C);
  
#define COMPUTE_POINT_RECT_BW_STAPLES(out,conf,sq_staples,A,D,E,imu,mu,inu,nu,temp) \
  COMPUTE_RECT_BW_STAPLE(out[A][mu][inu],sq_staples[D][nu][imu],conf[D][mu],conf[E][nu],temp); /*bw sq staple*/ \
  SUMM_RECT_BW_STAPLE(out[A][mu][inu],conf[D][nu],sq_staples[D][mu][inu],conf[E][nu],temp);    /*bw sq staple*/ \
  SUMM_RECT_BW_STAPLE(out[A][mu][inu],conf[D][nu],conf[D][mu],sq_staples[E][nu][3+imu],temp);  /*fw sq staple*/
  
  // 1) start communicating lower surface forward staples
  std::vector<MPI_Request> rectangular_staples_lx_conf_start_communicating_lower_surface_fw_squared_staples(const LxField<squared_staples_t>& sq_staples)
  {
    //copy lower surface into sending buf to be sent to dw nodes
    //obtained scanning on first half of the border, and storing them
    //in the first half of sending buf
    for(int nu=0;nu<4;nu++) //border and staple direction
      if(isDirParallel[nu])
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    const int mu=perpDirs[nu][imu];
	    const int inu=(nu<mu)?nu:nu-1;
	    
	    PAR(bordOffset[nu],
		bordOffset[nu]+bordDirVol[nu],
		CAPTURE(mu,
			inu,
			TO_READ(sq_staples)),
		ibord,
		{
		  su3_copy(((quad_su3*)send_buf)[ibord][mu],sq_staples[surflxOfBordlx[ibord]][mu][3+inu]); //one contribution per link in the border
		});
	  }
    
    //start communication of lower surf to backward nodes
    START_TIMING(tot_comm_time,ntot_comm);
    const std::vector<std::pair<int,int>> dir_comm={{0,0},{0,1},{0,2},{0,3}};
    
    return startBufHaloNeighExchange<quad_su3>(1,dir_comm);
  }
  
  // 2) compute non_fwsurf fw staples that are always local
  void rectangular_staples_lx_conf_compute_non_fw_surf_fw_staples(LxField<rectangular_staples_t>& out,
								  const LxField<quad_su3>& conf,
								  const LxField<squared_staples_t>& sq_staples)
  {
    for(int mu=0;mu<4;mu++) //link direction
      for(int inu=0;inu<3;inu++) //staple direction
	{
	  const int nu=perpDirs[mu][inu];
	  const int imu=(mu<nu)?mu:mu-1;
	  
	  PAR(0,
	      nonFwSurfVol,
	      CAPTURE(mu,
		      imu,
		      inu,
		      nu,
		      TO_WRITE(out),
		      TO_READ(conf),
		      TO_READ(sq_staples)),
	      ibulk,
	      {
		su3 temp; //three staples in clocwise order
		int A=loclxOfNonFwSurflx[ibulk],B=loclxNeighup[A][nu],F=loclxNeighup[A][mu];
		COMPUTE_POINT_RECT_FW_STAPLES(out,conf,sq_staples,A,B,F,imu,mu,inu,nu,temp);
	      });
	}
  }
  
  // 3) finish communication of lower surface fw squared staples
  void rectangular_staples_lx_conf_finish_communicating_lower_surface_fw_squared_staples(LxField<squared_staples_t>& sq_staples,
											 std::vector<MPI_Request>& requests)
  {
    waitAsyncCommsFinish(requests);
    STOP_TIMING(tot_comm_time);
    
    //copy the received forward border (stored in the second half of receiving buf) to its destination
    for(int nu=0;nu<4;nu++) //border and staple direction
      if(isDirParallel[nu])
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    int mu=perpDirs[nu][imu];
	    int inu=(nu<mu)?nu:nu-1;
	    
	    PAR(bordVolh+bordOffset[nu],
		bordVolh+bordOffset[nu]+bordDirVol[nu],
		CAPTURE(mu,
			inu,
			TO_WRITE(sq_staples)),
		ibord,
		{
		  su3_copy(sq_staples[locVol+ibord][mu][3+inu],((quad_su3*)recv_buf)[ibord][mu]); //one contribution per link in the border
		});
	  }
  }
  
  // 4) compute backward staples to be sent to up nodes and send them
  std::vector<MPI_Request> rectangular_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(LxField<rectangular_staples_t>& out,
													  const LxField<quad_su3>& conf,
													  const LxField<squared_staples_t>& sq_staples)
  {
    //compute backward staples to be sent to up nodes
    //obtained scanning D on fw_surf and storing data as they come
    for(int inu=0;inu<3;inu++) //staple direction
      for(int mu=0;mu<4;mu++) //link direction
	{
	  int nu=perpDirs[mu][inu];
	  int imu=(mu<nu)?mu:mu-1;
	  
	  PAR(0,
	      fwSurfVol,
	      CAPTURE(mu,
		      imu,
		      nu,
		      inu,
		      TO_WRITE(out),
		      TO_READ(conf),
		      TO_READ(sq_staples)),
	      ifw_surf,
	      {
		su3 temp;
		int D=loclxOfFwSurflx[ifw_surf],A=loclxNeighup[D][nu],E=loclxNeighup[D][mu];
		COMPUTE_POINT_RECT_BW_STAPLES(out,conf,sq_staples,A,D,E,imu,mu,inu,nu,temp);
	      });
	}
    
    //copy in send buf, obtained scanning second half of each parallelized direction external border and
    //copying the three perpendicular links staple
    for(int nu=0;nu<4;nu++) //border and staple direction
      if(isDirParallel[nu])
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    const int mu=perpDirs[nu][imu];
	    const int inu=(nu<mu)?nu:nu-1;
	    
	    PAR(bordVolh+bordOffset[nu],
		bordVolh+bordOffset[nu]+bordDirVol[nu],
		CAPTURE(mu,
			inu,
			TO_WRITE(out)),
		ibord,
		{
		  su3_copy(((quad_su3*)send_buf)[ibord][mu],out[locVol+ibord][mu][inu]); //one contribution per link in the border
		});
	  }
    
    //start communication of fw surf backward staples to forward nodes
    START_TIMING(tot_comm_time,ntot_comm);
    const std::vector<std::pair<int,int>> dir_comm={{1,0},{1,1},{1,2},{1,3}};
    
    return startBufHaloNeighExchange<quad_su3>(1,dir_comm);
  }
  
  // 5) compute non_fw_surf bw staples
  void rectangular_staples_lx_conf_compute_non_fw_surf_bw_staples(LxField<rectangular_staples_t>& out,
								  const LxField<quad_su3>& conf,
								  const LxField<squared_staples_t>& sq_staples)
  {
    for(int mu=0;mu<4;mu++) //link direction
      for(int inu=0;inu<3;inu++) //staple direction
	{
	  int nu=perpDirs[mu][inu];
	  int imu=(mu<nu)?mu:mu-1;
	  
	  //obtained scanning D on fw_surf
	  PAR(0,
	      nonFwSurfVol,
	      CAPTURE(nu,
		      inu,
		      mu,
		      imu,
		      TO_READ(conf),
		      TO_READ(sq_staples),
		      TO_WRITE(out)),
	      inon_fw_surf,
	      {
		su3 temp;
		const int D=loclxOfNonFwSurflx[inon_fw_surf],A=loclxNeighup[D][nu],E=loclxNeighup[D][mu];
		COMPUTE_POINT_RECT_BW_STAPLES(out,conf,sq_staples,A,D,E,imu,mu,inu,nu,temp);
	      });
	}
  }
  
  // 6) compute fw_surf fw staples
  void rectangular_staples_lx_conf_compute_fw_surf_fw_staples(LxField<rectangular_staples_t>& out,
							      const LxField<quad_su3>& conf,
							      const LxField<squared_staples_t>& sq_staples)
  {
    for(int mu=0;mu<4;mu++) //link direction
      for(int inu=0;inu<3;inu++) //staple direction
	{
	  int nu=perpDirs[mu][inu];
	  int imu=(mu<nu)?mu:mu-1;
	  
	  //obtained looping A on forward surface
	  PAR(0,
	      fwSurfVol,
	      CAPTURE(inu,
		      nu,
		      imu,
		      mu,
		      TO_WRITE(out),
		      TO_READ(conf),
		      TO_READ(sq_staples)),
	      ifw_surf,
	      {
		su3 temp;
		int A=loclxOfFwSurflx[ifw_surf],B=loclxNeighup[A][nu],F=loclxNeighup[A][mu];
		COMPUTE_POINT_RECT_FW_STAPLES(out,conf,sq_staples,A,B,F,imu,mu,inu,nu,temp);
	      });
	}
  }
  
  // 7) finish communication of fw_surf bw staples
  void rectangular_staples_lx_conf_finish_communicating_fw_surf_bw_staples(LxField<rectangular_staples_t>& out,
									   std::vector<MPI_Request> requests)
  {
    waitAsyncCommsFinish(requests);
    STOP_TIMING(tot_comm_time);
    
    //copy the received backward staples (stored on first half of receiving buf) on bw_surf sites
    for(int nu=0;nu<4;nu++) //staple and fw bord direction
      if(isDirParallel[nu])
	for(int imu=0;imu<3;imu++) //link direction
	  {
	    const int mu=perpDirs[nu][imu];
	    const int inu=(nu<mu)?nu:nu-1;
	    
	    PAR(bordOffset[nu],
		bordOffset[nu]+bordDirVol[nu],
		CAPTURE(mu,
			inu,
			TO_WRITE(out)),
		ibord,
		{
		  su3_copy(out[surflxOfBordlx[ibord]][mu][inu],((quad_su3*)recv_buf)[ibord][mu]);//one contribution per link in the border
		});
	  }
  }
  
  //compute rectangular staples overlapping computation and communications, and avoiding using edges
  void compute_rectangular_staples_lx_conf(LxField<rectangular_staples_t>& out,
					   const LxField<quad_su3>& conf,
					   const LxField<squared_staples_t>& sq_staples)
  {
    //compute non_fw_surf fw staples
    std::vector<MPI_Request> bwRequests=
      rectangular_staples_lx_conf_start_communicating_lower_surface_fw_squared_staples(sq_staples);
    rectangular_staples_lx_conf_compute_non_fw_surf_fw_staples(out,conf,sq_staples); //nasty
    rectangular_staples_lx_conf_finish_communicating_lower_surface_fw_squared_staples(const_cast<LxField<squared_staples_t>&>(sq_staples),bwRequests);
    
    //compute fw_surf bw staples, non_fw_surf bw staples and fw_surf fw staples
    std::vector<MPI_Request> fwRequests=
      rectangular_staples_lx_conf_compute_and_start_communicating_fw_surf_bw_staples(out,conf,sq_staples);
    rectangular_staples_lx_conf_compute_non_fw_surf_bw_staples(out,conf,sq_staples);
    rectangular_staples_lx_conf_compute_fw_surf_fw_staples(out,conf,sq_staples);
    rectangular_staples_lx_conf_finish_communicating_fw_surf_bw_staples(out,fwRequests);
  }
  
  //summ everything together
  void compute_summed_rectangular_staples_lx_conf(LxField<quad_su3>& out,
						  const LxField<quad_su3>& conf,
						  const LxField<squared_staples_t>& sq_staples)
  {
    LxField<rectangular_staples_t> rectangular_staples("rectangular_staples",WITH_HALO);
    
    compute_rectangular_staples_lx_conf(rectangular_staples,conf,sq_staples);
    
    //summ
    PAR(0,locVol,
	CAPTURE(TO_WRITE(out),
		TO_READ(rectangular_staples)),
	ivol,
	{
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      su3_copy(out[ivol][mu],rectangular_staples[ivol][mu][0]);
	      for(int iterm=1;iterm<2*(NDIM-1);iterm++)
		su3_summassign(out[ivol][mu],rectangular_staples[ivol][mu][iterm]);
	    }
	});
  }
}
