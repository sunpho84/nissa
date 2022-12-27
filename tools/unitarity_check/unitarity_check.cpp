#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <math.h>

#include "nissa.hpp"

using namespace nissa;

template <typename Conf>
void test_unitarity(FILE *fout,
		    LxField<quad_su3>& conf,
		    char *filename)
{
  double loc_max=0,loc_avg=0;
  
  read_ildg_gauge_conf(conf,filename);
  master_printf("conf has valid halo: %d\n",conf.haloIsValid);
  conf.invalidateHalo();
  conf.updateHalo();
  master_printf("forced update, conf has valid halo: %d\n",conf.haloIsValid);
  PAR(0,1,
      CAPTURE(TO_WRITE(conf)),
      ivol,
      {
      });
  master_printf("writing, conf has valid halo: %d\n",conf.haloIsValid);
  master_printf("Plaquette: %16.16lg\n",global_plaquette_lx_conf(conf));
  
  // OddField<oct_su3> test("test");
  // conf.updateHalo();
  
  // NISSA_PARALLEL_LOOP(iev,0,locVolh)
  //   {
  //     for(int mu=0;mu<NDIM;mu++)
  // 	{
  // 	  const int ilx=loclx_of_loceo[EVN][iev];
  // 	  su3_copy(test[iev][0][mu],conf[loclxNeighdw[ilx][mu]][mu]);
  // 	  su3_copy(test[iev][1][mu],conf[ilx][mu]);
  // 	}
  //   }
  // NISSA_PARALLEL_LOOP_END;
  
  // using HalfStaple=su3[2][NDIM][2][NDIM-1];
  // using ComprHalfStaple=su3[2][NDIM-1];
  
  // OddField<HalfStaple,WITH_HALO> stapPart("stapPart");
  // NISSA_PARALLEL_LOOP(iev,0,locVolh)
  //   {
  //     for(int mu=0;mu<NDIM;mu++)
  // 	for(int inu=0;inu<NDIM-1;inu++)
  // 	  {
  // 	    const int nu=perp_dir[mu][inu];
  // 	    unsafe_su3_prod_su3(stapPart[loceo_neighup[EVN][iev][nu]][1][mu][1][inu],
  // 				test[iev][1][mu],
  // 				test[iev][1][nu]);
  // 		     }
  // 	    }
  // NISSA_PARALLEL_LOOP_END;
  
  // stapPart.fillSendingBufWithHalo<ComprHalfStaple>([](ComprHalfStaple& out,
  // 						      const auto& in,
  // 						      const int& bf,
  // 						      const int& mu)
  // {
  //   for(int bf2=0;bf2<2;bf2++)
  //     for(int inu=0;inu<NDIM-1;inu++)
  // 	su3_copy(out[bf2][inu],in[bf][mu][bf2][inu]);
  // });
  
  // exchangeNeighBuf<ComprHalfStaple>(/*half vol*/ 2);
  
  // stapPart.fillSurfaceWithReceivingBuf<ComprHalfStaple>([](auto&& out,
  // 							   const ComprHalfStaple& in,
  // 							   const int& bf,
  // 							   const int& mu)
  // {
  //   for(int bf2=0;bf2<2;bf2++)
  //     for(int inu=0;inu<NDIM-1;inu++)
  // 	su3_copy(out[bf][mu][bf2][inu],in[bf2][inu]);
  // });
  
  NISSA_LOC_VOL_LOOP(ivol)
    for(int idir=0;idir<4;idir++)
      {
	su3 zero;
	su3_put_to_id(zero);
	su3_subt_the_prod_su3_dag(zero,conf[ivol][idir],conf[ivol][idir]);
	
	double r=real_part_of_trace_su3_prod_su3_dag(zero,zero)/18;
	
	loc_avg+=r;
	if(loc_max<r) loc_max=r;
	
	if(0)
	  if(r>1.e-30)
	    {
	      master_printf("diff %d %d %d %d   %d   %lg\n",glbCoordOfLoclx[ivol][0],glbCoordOfLoclx[ivol][1],
			    glbCoordOfLoclx[ivol][2],glbCoordOfLoclx[ivol][3],idir,r);
	      su3_print(conf[ivol][idir]);
	      for(int i=0;i<3;i++)
		for(int j=i;j<3;j++)
		  {
		    complex t;
		    color_scalar_prod(t,conf[ivol][idir][i],conf[ivol][idir][j]);
		    
		    // if(fabs(t[0])>1.e-15 && fabs(t[0]-1)>1.e-15)
		    //   {
		    //     printf(" %d%d prod: %lg %lg\n",i,j,t[0],t[1]);
		    
		    //     //search for orthogonals
		    //     for(int jvol=0;jvol<locVol*4*9-18;jvol++)
		    // 	{
		    // 	  color_scalar_prod(t,conf[ivol][idir][i],(*((color*)((complex*)conf+jvol))));
		    // 	  if(fabs(t[0])<=1.e-15)
		    // 	    printf(" %d orth to %d (%d %d %d), prod: %lg %lg\n",
		    // 		   ivol,jvol,jvol/36,(jvol%36)/9,jvol%9,t[0],t[1]);
		    // 	}
		    //   }
		  }
	    }
      }
  
  double glb_max=0,glb_avg=0;
  MPI_Reduce(&loc_avg,&glb_avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&loc_max,&glb_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  glb_avg/=4*glbVol;
  
  glb_avg=sqrt(glb_avg);
  glb_max=sqrt(glb_max);
  
  master_fprintf(fout,"%s, Max: %16.16lg, Avg: %16.16lg\n",filename,glb_max,glb_avg);
}

int main(int narg,char **arg)
{
  char filename[1024];
  
  //basic mpi initialization
  init_nissa(narg,arg);
  
  {
    if(narg<2) crash("Use: %s input_file",arg[0]);
    
    open_input(arg[1]);
    
    //grid sizes
    int L,T;
    read_str_int("L",&L);
    read_str_int("T",&T);
    
    //init the MPI grid
    init_grid(T,L);
    
    testLxHaloExchange();
    testEoHaloExchange();
    testLxEdgesExchange();
    testEoEdgesExchange();
    
    //summary
    char output[1024];
    read_str_str("SummaryFile",output,1024);
    
    FILE *fout=open_text_file_for_output(output);
    
    /// Number of configurations
    int nconf;
    read_str_int("NGaugeConf",&nconf);
    
    LxField<quad_su3> conf("conf",WITH_HALO);
    for(int iconf=0;iconf<nconf;iconf++)
      {
	read_str(filename,1024);
	test_unitarity(fout,conf,filename);
      }
    
    close_input();
    
    ///////////////////////////////////////////
    
    if(rank==0) fclose(fout);
  }
  
  close_nissa();
  
  return 0;
}
