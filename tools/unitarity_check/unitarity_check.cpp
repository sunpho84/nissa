#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <math.h>

#include "nissa.hpp"

using namespace nissa;

void test_unitarity(FILE *fout,quad_su3 *conf,char *filename)
{
  double loc_max=0,loc_avg=0;
  
  read_ildg_gauge_conf(conf,filename);
  master_printf("Plaquette: %16.16lg\n",global_plaquette_lx_conf(conf));
  
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
		  
		  if(fabs(t[0])>1.e-15 && fabs(t[0]-1)>1.e-15)
		    {
		      printf(" %d%d prod: %lg %lg\n",i,j,t[0],t[1]);
		      
		      //search for orthogonals
		      for(int jvol=0;jvol<locVol*4*9-18;jvol++)
			{
			  color_scalar_prod(t,conf[ivol][idir][i],(*((color*)((complex*)conf+jvol))));
			  if(fabs(t[0])<=1.e-15)
			    printf(" %d orth to %d (%d %d %d), prod: %lg %lg\n",
				   ivol,jvol,jvol/36,(jvol%36)/9,jvol%9,t[0],t[1]);
			}
		    }
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

using C=Field<quad_su3,nissa::CPU_LAYOUT>;
using G=Field<quad_su3,nissa::GPU_LAYOUT>;

//The product of two complex number
template <typename A,
	  typename B,
	  typename C>
CUDA_HOST_AND_DEVICE INLINE_FUNCTION
void unsafe_complex_prod(A&& a,const B& b,const C& c)
{
  a[0]=b[0]*c[0]-b[1]*c[1];
  a[1]=b[0]*c[1]+b[1]*c[0];
}

//Summ to the output the product of two complex number
template <typename A,
	  typename B,
	  typename C>
CUDA_HOST_AND_DEVICE INLINE_FUNCTION
void complex_summ_the_prod(A&& a,const B& b,const C& c)
{
  const auto t=b[0]*c[0]-b[1]*c[1];
  a[1]+=b[0]*c[1]+b[1]*c[0];
  a[0]+=t;
}

//Product of two su3 matrixes
template <typename A,
	  typename B,
	  typename C>
CUDA_HOST_AND_DEVICE INLINE_FUNCTION
void unsafe_su3_prod_su3(A&& a,const B& b,const C& c)
{
  for(int ir_out=0;ir_out<NCOL;ir_out++)
    for(int ic_out=0;ic_out<NCOL;ic_out++)
      {
	unsafe_complex_prod(a[ir_out][ic_out],b[ir_out][0],c[0][ic_out]);
	for(int itemp=1;itemp<NCOL;itemp++)
	  complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
      }
}

void f(C& c,const int ivol)
{
  ASM_BOOKMARK_BEGIN("CIAO");
  unsafe_su3_prod_su3(c[ivol+1][1],c[ivol][2],c[ivol][0]);
  ASM_BOOKMARK_END("CIAO");
}

void f(quad_su3* c,const int ivol)
{
  ASM_BOOKMARK_BEGIN("GIAO");
  unsafe_su3_prod_su3(c[ivol+1][1],c[ivol][2],c[ivol][0]);
  ASM_BOOKMARK_END("GIAO");
}

int main(int narg,char **arg)
{
  char filename[1024];
  
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  open_input(arg[1]);
  
  //grid sizes
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //init the MPI grid 
  init_grid(T,L);
  //summary
  char output[1024];
  read_str_str("SummaryFile",output,1024);
  FILE *fout=open_text_file_for_output(output);
  //nconf
  int nconf;
  read_str_int("NGaugeConf",&nconf);
  
  quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol,quad_su3);
  
  for(int iconf=0;iconf<nconf;iconf++)
    {
      read_str(filename,1024);
      test_unitarity(fout,conf,filename);
    }
  
  close_input();

    G g("g",locVol);
    C c("c",locVol);
  
    //f(g,10);
  {
    
    
    master_printf("%d %d\n",&g[5][3][2][1][1],&c[5][3][2][1][1]);
    
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	auto& r=f[ivol][1][0][0][1];
	
    // 	r=0;
    //   }
    // NISSA_PARALLEL_LOOP_END;
  }
  
  ///////////////////////////////////////////
  
  if(rank==0) fclose(fout);
  
  nissa_free(conf);

  close_nissa();

  return 0;
}
