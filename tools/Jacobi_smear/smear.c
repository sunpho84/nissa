#include "nissa.h"

int main(int narg,char **arg)
{
  int or_pos[4]={31,10,19,13};
  char filename[1024];

  //basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  read_str_str("Filename",filename,1024);

  close_input();
  
  //Init the MPI grid 
  init_grid(T,L);
  
  ///////////////////////////////////////////
  
  //allocate and read conf
  quad_su3 *orig_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"orig_conf");
  quad_su3 *smea_conf=allocate_quad_su3(loc_vol+loc_bord+loc_edge,"smea_conf");
  read_gauge_conf(orig_conf,filename);

  ape_smearing(smea_conf,orig_conf,0.5,20);
  if(rank==0) printf("gauge conf smeared\n");
  //memcpy(smea_conf,orig_conf,sizeof(quad_su3)*loc_vol);

  //allocate and generate the source
  spincolor *origi_sp=allocate_spincolor(loc_vol+loc_bord,"orig_spincolor");
  memset(origi_sp,0,sizeof(spincolor)*loc_vol);

  //put the source
  int local=1;
  int lx[4];
  for(int mu=0;mu<4;mu++)
    {
      lx[mu]=or_pos[mu]-proc_coord[mu]*loc_size[mu];
      local=local && lx[mu]>=0 && lx[mu]<loc_size[mu];
    }
  if(local==1) origi_sp[loclx_of_coord(lx)][0][0][0]=1;
  
  //allocate and smeard
  spincolor *smear_sp=allocate_spincolor(loc_vol+loc_bord,"smear_spincolor");
  dina_smearing(smear_sp,origi_sp,smea_conf,4,50,or_pos[0]);
  
  int L=glb_size[1];
  double n_or[L],n_sm[L];
  density_profile(n_or,origi_sp,or_pos);
  density_profile(n_sm,smear_sp,or_pos);

  FILE *fout=open_text_file_for_output("profile");
  for(int d=0;d<L;d++)
    {
      n_or[d]=sqrt(n_or[d]);
      n_sm[d]=sqrt(n_sm[d]);
      if(rank==0) fprintf(fout,"%d %g %g\n",d,n_or[d],n_sm[d]);
    }
  if(rank==0) fclose(fout);
    
  ///////////////////////////////////////////
  
  close_nissa();

  return 0;
}
