#include "nissa.h"
const double rad2=1.414213562373095048801688724209;

void apply_D(spincolor *out,quad_su3 *conf,spincolor *in,double div)
{
  int Xup,Xdw;

  for(int X=0;X<loc_vol;X++)
  {
    color temp_c0,temp_c1,temp_c2,temp_c3;

    //Forward 0
    Xup=loclx_neighup[X][0];
    color_summ(temp_c0,in[Xup][0],in[Xup][2]);
    color_summ(temp_c1,in[Xup][1],in[Xup][3]);
    unsafe_su3_prod_color(out[X][0],conf[X][0],temp_c0);
    unsafe_su3_prod_color(out[X][1],conf[X][0],temp_c1);
    color_copy(out[X][2],out[X][0]);
    color_copy(out[X][3],out[X][1]);
        
    //Backward 0
    Xdw=loclx_neighdw[X][0];
    color_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
    color_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
    unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][0],temp_c0);
    unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][0],temp_c1);
    summassign_color(out[X][0],temp_c2);
    summassign_color(out[X][1],temp_c3);
    subtassign_color(out[X][2],temp_c2);
    subtassign_color(out[X][3],temp_c3);

    //Forward 1
    Xup=loclx_neighup[X][1];
    color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
    color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
    unsafe_su3_prod_color(temp_c2,conf[X][1],temp_c0);
    unsafe_su3_prod_color(temp_c3,conf[X][1],temp_c1);
    summassign_color(out[X][0],temp_c2);
    summassign_color(out[X][1],temp_c3);
    subtassign_icolor(out[X][2],temp_c3);
    subtassign_icolor(out[X][3],temp_c2);
        
    //Backward 1
    Xdw=loclx_neighdw[X][1];
    color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
    color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
    unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][1],temp_c0);
    unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][1],temp_c1);
    summassign_color(out[X][0],temp_c2);
    summassign_color(out[X][1],temp_c3);
    summassign_icolor(out[X][2],temp_c3);
    summassign_icolor(out[X][3],temp_c2);
    
    //Forward 2
    Xup=loclx_neighup[X][2];
    color_summ(temp_c0,in[Xup][0],in[Xup][3]);
    color_subt(temp_c1,in[Xup][1],in[Xup][2]);
    unsafe_su3_prod_color(temp_c2,conf[X][2],temp_c0);
    unsafe_su3_prod_color(temp_c3,conf[X][2],temp_c1);
    summassign_color(out[X][0],temp_c2);
    summassign_color(out[X][1],temp_c3);
    subtassign_color(out[X][2],temp_c3);
    summassign_color(out[X][3],temp_c2);
        
    //Backward 2
    Xdw=loclx_neighdw[X][2];
    color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
    color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
    unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][2],temp_c0);
    unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][2],temp_c1);
    summassign_color(out[X][0],temp_c2);
    summassign_color(out[X][1],temp_c3);
    summassign_color(out[X][2],temp_c3);
    subtassign_color(out[X][3],temp_c2);
    
    //Forward 3
    Xup=loclx_neighup[X][3];
    color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
    color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
    unsafe_su3_prod_color(temp_c2,conf[X][3],temp_c0);
    unsafe_su3_prod_color(temp_c3,conf[X][3],temp_c1);
    summassign_color(out[X][0],temp_c2);
    summassign_color(out[X][1],temp_c3);
    subtassign_icolor(out[X][2],temp_c2);
    summassign_icolor(out[X][3],temp_c3);
        
    //Backward 3
    Xdw=loclx_neighdw[X][3];
    color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
    color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
    unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][3],temp_c0);
    unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][3],temp_c1);
    summassign_color(out[X][0],temp_c2);
    summassign_color(out[X][1],temp_c3);
    summassign_icolor(out[X][2],temp_c2);
    subtassign_icolor(out[X][3],temp_c3);
    
    //Put the -1/2 factor on derivative, the gamma5, and the imu
    //ok this is horrible, but fast
    for(int d=0;d<4;d++)
      for(int c=0;c<3;c++)
	for(int r=0;r<2;r++)
	  out[X][d][c][r]*=-0.5/div;
  }
}

void stoc_source(spincolor *spinore)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
        {
	  spinore[ivol][id][ic][0]=pm_one(ivol)/rad2;
	  spinore[ivol][id][ic][1]=pm_one(ivol)/rad2;
	}
}

void scal_prod(complex c,spincolor *a,spincolor *b)
{
  complex loc_c={0,0};
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	complex_summ_the_conj1_prod(loc_c,a[ivol][id][ic],b[ivol][id][ic]);
  MPI_Allreduce(loc_c,c,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

int main(int narg,char **arg)
{
  int seed=0;
  char filename[1024];

  //basic mpi initialization
  init_nissa();
  
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  open_input(arg[1]);

  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //init the MPI grid 
  init_grid(T,L);

  read_str_int("Seed",&seed);
  read_str_str("Filename",filename,1024);
  
  close_input();

  //Initialize the random generator
  init_random(seed);

  ///////////////////////////////////////////

  quad_su3 *conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol));
  read_ildg_gauge_conf(conf,filename);
  
  spincolor *eta=allocate_spincolor(loc_vol+bord_vol,"eta");
  spincolor *source=allocate_spincolor(loc_vol+bord_vol,"source");
  spincolor *phi=allocate_spincolor(loc_vol+bord_vol,"phi");
  
  int ns=100;
  int nt=4;
  complex tr[nt][ns+1];memset(tr,0,sizeof(complex)*nt*(ns+1));

  for(int is=0;is<ns;is++)
    {
      stoc_source(eta);
      memcpy(source,eta,loc_vol*sizeof(spincolor));
      
      for(int it=0;it<nt;it++)
	{
	  double div=(it>0)?it:1;

	  apply_D(phi,conf,source,div);
	  scal_prod(tr[it][is],eta,phi);
	  memcpy(source,phi,loc_vol*sizeof(spincolor));
	}
    }
  
  for(int it=1;it<nt;it+=2)
    {
      for(int is=0;is<ns;is++) complex_summassign(tr[it][ns],tr[it][is]);
      for(int is=0;is<ns;is++) complex_subt(tr[it][is],tr[it][ns],tr[it][is]);
    printf("%d %lg %lg\n",it+1,-tr[it][is][0]/ns,-tr[it][is][1]/ns);
  
  return 0;
}
