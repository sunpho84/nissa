#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <math.h>

#include "nissa.hpp"

using namespace nissa;

template <typename Conf>
void test_unitarity(FILE *fout,
		    Conf& conf,
		    char *filename)
{
  double loc_max=0,loc_avg=0;
  
  read_ildg_gauge_conf(conf,filename);
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

std::string siteAsString(const int& n)
{
  const auto c=glb_coord_of_glblx(n);
  
  std::string res=
    std::to_string(n)+
    "("+std::to_string(c[0]);
  
  for(int nu=1;nu<NDIM;nu++)
    res+=","+std::to_string(c[nu]);
  res+=")";
  
  return res;
}

void taintTheCommBuffers()
{
  int* r=(int*)recv_buf;;
  NISSA_PARALLEL_LOOP(i,0,recv_buf_size/sizeof(int))
    r[i]=-3;
  NISSA_PARALLEL_LOOP_END;
  
  int* s=(int*)send_buf;
  NISSA_PARALLEL_LOOP(i,0,send_buf_size/sizeof(int))
    s[i]=-4;
  NISSA_PARALLEL_LOOP_END;
}

void testLxHaloExchange()
{
  LxField<int> test("testHalo",WITH_HALO);
  NISSA_PARALLEL_LOOP(i,0,locVol)
    test[i]=glblxOfLoclx[i];
  NISSA_PARALLEL_LOOP_END;
  
  test.invalidateHalo();
  
  taintTheCommBuffers();
  
  test.updateHalo();
  
  NISSA_PARALLEL_LOOP(site,0,locVol)
    {
      for(int ori=0;ori<2;ori++)
	for(int mu=0;mu<NDIM;mu++)
	  {
	    const int ln=loclx_neigh[ori][site][mu];
	    const int gn=(ln<locVol)?glblxOfLoclx[ln]:glblxOfBordlx[ln-locVol];
	    const int neighVal=test[ln];
	    
	    if(neighVal!=gn)
	      master_printf("site %s ori %d dir %d has neigh %s with val %s\n",
			    siteAsString(glblxOfLoclx[site]).c_str(),
			    ori,mu,
			    siteAsString(gn).c_str(),
			    siteAsString(neighVal).c_str());
	  }
    }
  NISSA_PARALLEL_LOOP_END;
}

void testEoHaloExchange()
{
  EoField<int> test("testEoHalo",WITH_HALO);
  forBothParities([&test](const auto& par)
    {
      NISSA_PARALLEL_LOOP(i,0,locVolh)
	test[par][i]=glblxOfLoclx[loclx_of_loceo[par][i]];
      NISSA_PARALLEL_LOOP_END;
    });
  test.invalidateHalo();
  
  taintTheCommBuffers();
  
  test.updateHalo();
  
  forBothParities([&test](const auto& par)
    {
      NISSA_PARALLEL_LOOP(eoSite,0,locVolh)
	{
	  for(int ori=0;ori<2;ori++)
	    for(int mu=0;mu<NDIM;mu++)
	      {
		const int ln=((ori==0)?loceo_neighdw:loceo_neighup)[par][eoSite][mu];
		const int lx=loclx_of_loceo[!par][ln];
		const int gn=(lx<locVol)?glblxOfLoclx[lx]:glblxOfBordlx[lx-locVol];
		const int neighVal=test[!par][ln];
		
		if(neighVal!=gn)
		  master_printf("site %s ori %d dir %d has neigh %s with val %s\n",
				siteAsString(glblxOfLoclx[loclx_of_loceo[!par][eoSite]]).c_str(),
				ori,mu,
				siteAsString(gn).c_str(),
				siteAsString(neighVal).c_str());
	      }
	}
      NISSA_PARALLEL_LOOP_END;
    });
}

void testLxEdgesExchange()
{
  LxField<int> test("testEdge",WITH_HALO_EDGES);
  NISSA_PARALLEL_LOOP(i,0,locVol)
    test[i]=glblxOfLoclx[i];
  NISSA_PARALLEL_LOOP_END;
  
  NISSA_PARALLEL_LOOP(i,0,bord_vol)
    test[i+locVol]=-1;
  NISSA_PARALLEL_LOOP_END;
  
  NISSA_PARALLEL_LOOP(i,0,edge_vol)
    test[i+locVol+bord_vol]=-2;
  NISSA_PARALLEL_LOOP_END;
  
  test.invalidateHalo();
  test.invalidateEdges();
  test.updateHalo();
  
  taintTheCommBuffers();
  
  test.updateEdges();
  
  NISSA_PARALLEL_LOOP(site,0,locVol)
    {
      for(int ori1=0;ori1<2;ori1++)
	for(int ori2=0;ori2<2;ori2++)
	  for(int iEdge=0;iEdge<nEdges;iEdge++)
	    {
	      const auto [mu,nu]=edge_dirs[iEdge];
	      
	      const int l1n=loclx_neigh[ori1][site][mu];
	      const int ln=loclx_neigh[ori2][l1n][nu];
	      const int gn=(ln<locVol)?glblxOfLoclx[ln]:((ln<locVol+bord_vol)?glblxOfBordlx[ln-locVol]:glblxOfEdgelx[ln-locVol-bord_vol]);
	      const int neighVal=test[ln];
	      
	      if(neighVal!=gn)
		master_printf("site %s ori (%d,%d) dir (%d,%d) has edgelx neigh %s with val %s\n",
			      siteAsString(glblxOfLoclx[site]).c_str(),
			      ori1,ori2,mu,nu,
			      siteAsString(gn).c_str(),
			      siteAsString(neighVal).c_str());
	    }
    }
  NISSA_PARALLEL_LOOP_END;
}

void testEoEdgesExchange()
{
  EoField<int> test("testEoEdge",WITH_HALO_EDGES);
  forBothParities([&test](const auto& par)
    {
      NISSA_PARALLEL_LOOP(i,0,locVolh)
	test[par][i]=glblxOfLoclx[loclx_of_loceo[par][i]];
      NISSA_PARALLEL_LOOP_END;
      
      NISSA_PARALLEL_LOOP(i,0,bord_volh)
	test[par][i+locVolh]=-1;
      NISSA_PARALLEL_LOOP_END;
      
      NISSA_PARALLEL_LOOP(i,0,edge_volh)
	test[par][i+locVolh+bord_volh]=-2;
      NISSA_PARALLEL_LOOP_END;
    });
  
  int* r=(int*)recv_buf;
  int* s=(int*)send_buf;
  const int n=recv_buf_size/sizeof(int);
  NISSA_PARALLEL_LOOP(i,0,n)
    {
      r[i]=-3;
      s[i]=-4;
    }
  NISSA_PARALLEL_LOOP_END;
  
  test.invalidateHalo();
  test.invalidateEdges();
  test.updateHalo();
  test.updateEdges();
  
  forBothParities([&test](const auto& par)
  {
    NISSA_PARALLEL_LOOP(site,0,locVolh)
      {
	for(int ori1=0;ori1<2;ori1++)
	  for(int ori2=0;ori2<2;ori2++)
	    for(int iEdge=0;iEdge<nEdges;iEdge++)
	      {
		const auto [mu,nu]=edge_dirs[iEdge];
		
		const int l1n=((ori1==0)?loceo_neighdw:loceo_neighup)[par][site][mu];
		const int ln=((ori2==0)?loceo_neighdw:loceo_neighup)[!par][l1n][nu];
		const int lx=loclx_of_loceo[par][ln];
		const int gn=(lx<locVol)?glblxOfLoclx[lx]:((lx<locVol+bord_vol)?glblxOfBordlx[lx-locVol]:glblxOfEdgelx[lx-locVol-bord_vol]);
		const int neighVal=test[par][ln];
	      
		if(neighVal!=gn)
		  master_printf("par %d site %s ori (%d,%d) dir (%d,%d) neigh %s with val %s\n",
				par(),
				siteAsString(glblxOfLoclx[loclx_of_loceo[par][site]]).c_str(),
				ori1,ori2,mu,nu,
				siteAsString(gn).c_str(),
				siteAsString(neighVal).c_str());
	      }
      }
    NISSA_PARALLEL_LOOP_END;
  });
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
