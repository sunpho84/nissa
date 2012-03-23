#include "nissa.h"

//compute the static quark pair potential for time going from 1 to maxT and distance from 1 to maxD
//writing the output to the passed path 
void compute_static_quark_pair_potential(char *outpath,quad_su3 *conf,int maxT,int maxD)
{
  //open outpath
  FILE *fout=open_text_file_for_output(outpath);
  
  //allocate product of links
  su3 *u=nissa_malloc("U",loc_vol+loc_bord,su3);

  //fix temporal direction
  int mu=0;
  
  //loop over the distance
  for(int nstep_D=1;nstep_D<=maxD;nstep_D++)
    {
      master_fprintf(fout,"# potential at distance %d between 1 and %d\n",nstep_D,maxT);
      
      //loop over time
      for(int nstep_T=1;nstep_T<=maxT;nstep_T++)
	{
	  //loops averaging on orthogonal dirs
	  double rec=0;
	  for(int nu=1;nu<4;nu++)
	    rec+=average_real_part_of_trace_of_rectangle_path(conf,mu,nu,nstep_T,nstep_D,u);
	  master_fprintf(fout,"%.18lg\n",rec/3/3);
	  if(rank==0) fflush(fout);
	}
      master_fprintf(fout,"\n");
    }
  
  //close output
  if(rank==0) fclose(fout);
  
  //free product
  nissa_free(u);
}

int main(int narg,char **arg)
{
  init_nissa();

  //open input
  if(narg<2) crash("Use: %s input_file",arg[0]);
  open_input(arg[1]);
  
  //Init the MPI grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  char conf_path[1024];
  read_str_str("GaugePath",conf_path,1024);
  char out_path[1024];
  read_str_str("OutPath",out_path,1024);
  double alpha0,alpha1,alpha2;
  read_str_double("Alpha0",&alpha0);
  read_str_double("Alpha1",&alpha1);
  read_str_double("Alpha2",&alpha2);
  int Tmax,Dmax;
  read_str_int("Dmax",&Dmax);
  read_str_int("Tmax",&Tmax);
  
  close_input();
  
  /////////////////////////////////
  
  //read conf
  quad_su3 *read_conf=nissa_malloc("read_conf",loc_vol+loc_bord+loc_edge,quad_su3);
  read_ildg_gauge_conf(read_conf,conf_path);
  communicate_lx_quad_su3_borders(read_conf);
  master_printf(" read conf plaq: %.18g\n",global_plaquette_lx_conf(read_conf));
  
  //hyp conf
  quad_su3 *hyp_smeared_conf=nissa_malloc("hyp_smeared_conf",loc_vol+loc_bord+loc_edge,quad_su3);
  hyp_smear_conf_dir(hyp_smeared_conf,read_conf,alpha0,alpha1,alpha2,0);
  communicate_lx_quad_su3_borders(hyp_smeared_conf);
  master_printf(" hyp smeared conf plaq: %.18g\n",global_plaquette_lx_conf(hyp_smeared_conf));
  
  //smear conf
  ape_spatial_smear_conf(read_conf,read_conf,0.25,8);
  communicate_lx_quad_su3_borders(read_conf);
  master_printf(" ape smeared conf plaq: %.18g\n",global_plaquette_lx_conf(read_conf));
  
  //overwrite spatial links
  nissa_loc_vol_loop(ivol)
    for(int mu=1;mu<3;mu++)
      su3_copy(hyp_smeared_conf[ivol][mu],read_conf[ivol][mu]);
  
  compute_static_quark_pair_potential(out_path,hyp_smeared_conf,Tmax,Dmax);
  
  //free
  nissa_free(hyp_smeared_conf);
  nissa_free(read_conf);
  
  close_nissa();
  
  return 0;
}
