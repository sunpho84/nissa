#include "nissa.h"

int main(int narg,char **arg)
{
  int or_pos[4]={0,0,0,0};
  char filename[1024];
  
  //basic mpi initialization
  initNissa();

  if(narg<2) CRASH("Use: %s input_file",arg[0]);

  open_input(arg[1]);

  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  //Init the MPI grid 
  initGrid(T,L);
  
  //Smearing parameters
  double ape_alpha;
  read_str_double("ApeAlpha",&ape_alpha);
  int ape_niters;
  read_str_int("ApeNiter",&ape_niters);
  double jacobi_kappa;
  read_str_double("JacobiKappa",&jacobi_kappa);
  int jacobi_niters;
  read_str_int("JacobiNiters",&jacobi_niters);
  
  ///////////////////////////////////////////
  
  //allocate and read conf
  read_str_str("GaugeConf",filename,1024);
  quad_su3 *orig_conf=nissa_malloc("Ori_Conf",loc_vol+bord_vol+edge_vol,quad_su3);
  read_ildg_gauge_conf(orig_conf,filename);
  
  //smear the conf
  quad_su3 *smea_conf=nissa_malloc("Smea_Conf",loc_vol+bord_vol+edge_vol,quad_su3);
  ape_spatial_smear_conf(smea_conf,orig_conf,ape_alpha,ape_niters);
  MASTER_PRINTF("gauge conf smeared\n");
  
  //allocate and generate the source
  spincolor *sp=nissa_malloc("orig_spincolor",loc_vol+bord_vol,spincolor);
  memset(sp,0,sizeof(spincolor)*loc_vol);

  //put the source
  int local=1;
  int lx[4];
  for(int mu=0;mu<4;mu++)
    {
      lx[mu]=or_pos[mu]-rank_coord[mu]*loc_size[mu];
      local=local && lx[mu]>=0 && lx[mu]<loc_size[mu];
    }
  if(local==1) sp[loclxOfCoord(lx)][0][0][0]=1;
  
  //read the output base
  char base_out[1024];
  read_str_str("BaseOut",base_out,1024);
  
  close_input();
  
  /////////////////////////////////////////////
  
  for(int ismear=0;ismear<jacobi_niters;ismear++)
    {
      //comput the density profile
      int r=sqrt(glb_size[0]*glb_size[0]+3*glb_size[1]*glb_size[1]);
      double rho[r];
      density_profile(rho,sp,or_pos);
      
      //write the profile
      char path[1024];
      sprintf(path,"%s_%02d",base_out,ismear);
      FILE *fout=open_text_file_for_output(path);
      for(int d=0;d<r;d++) master_fprintf(fout,"%d %g\n",d,rho[d]);
      if(rank==0) fclose(fout);
      
      jacobi_smearing(sp,sp,smea_conf,jacobi_kappa,1);
    }
  
  ///////////////////////////////////////////
  
  nissa_free(smea_conf);
  nissa_free(orig_conf);
  nissa_free(sp);
  
  ///////////////////////////////////////////
  
  closeNissa();

  return 0;
}
