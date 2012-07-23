#include <string.h>

#include "nissa.h"

//compute the static quark pair potential
void compute_static_quark_pair_potential(char *out_path,quad_su3 *conf,double hyp_alpha0,double hyp_alpha1,double hyp_alpha2,double ape_alpha,int nlev_spat_smear,
					 int *niter_lev_spat_smear,int Tmax,int Dmax)
{
  //open output
  FILE *fout=open_text_file_for_output(out_path);
      
  //hyp conf
  hyp_smear_conf_dir(conf,conf,hyp_alpha0,hyp_alpha1,hyp_alpha2,0);
  communicate_lx_quad_su3_borders(conf);
  master_printf(" hyp smeared conf plaq: %.18g\n",global_plaquette_lx_conf(conf));
  
  //allocate su3 line in time
  su3 *tline=nissa_malloc("tline",loc_vol*Tmax,su3);
  
  //allocate su3 line in space
  su3 *xline=nissa_malloc("dline",loc_vol*Dmax,su3);
  
  /////////////////////////////////
  
  //compute the tline: this is done once forever since we are using only one level of hyp
  master_printf("Computing tline\n");
  nissa_loc_vol_loop(ivol)
    {
      const int mu=0;
      
      su3_copy(tline[ivol*Tmax],conf[ivol][mu]);
      int jvol=ivol;
      for(int t=1;t<Tmax;t++)
	{
	  jvol=loclx_neighup[jvol][mu];
	  unsafe_su3_prod_su3(tline[ivol*Tmax+t],tline[ivol*Tmax+t-1],conf[jvol][mu]);
	}
    }
  
  //loop over spatial smearing
  int niter_prec_lev_spat_smear=0;
  for(int ilev_spat_smear=0;ilev_spat_smear<nlev_spat_smear;ilev_spat_smear++)
    {
      int niter_to_apply=niter_lev_spat_smear[ilev_spat_smear]-niter_prec_lev_spat_smear;
      master_printf("\nComputing potential for level smearing: %d, niter: %d.\n",ilev_spat_smear,niter_lev_spat_smear[ilev_spat_smear]);
      
      //smear conf
      ape_spatial_smear_conf(conf,conf,ape_alpha,niter_to_apply);
      niter_prec_lev_spat_smear=niter_lev_spat_smear[ilev_spat_smear];
      master_printf(" hyp + ape smeared conf plaq: %.18g\n",global_plaquette_lx_conf(conf));

      //reset the potential
      double U[Dmax][Tmax];
      memset(U,0,sizeof(double)*Tmax*Dmax);
      
      //loop over the three spatial directions
      for(int mu=1;mu<4;mu++)
	{
	  //compute the xline
	  master_printf("Computing xline for direction %d\n",mu);
	  nissa_loc_vol_loop(ivol)
	    {
	      su3_copy(xline[ivol*Dmax],conf[ivol][mu]);
	      int jvol=ivol;
	      for(int d=1;d<Dmax;d++)
		{
		  jvol=loclx_neighup[jvol][mu];
		  unsafe_su3_prod_su3(xline[ivol*Dmax+d],xline[ivol*Dmax+d-1],conf[jvol][mu]);
		}
	    }
	  
	  //loop over the whole lattice
	  nissa_loc_vol_loop(ivol)
	    {
	      int A=ivol;
	      
	      //loop over time distance
	      for(int t=0;t<Tmax;t++)
		{
		  //find D
		  int D=A;
		  for(int i=0;i<=t;i++) D=loclx_neighup[D][0];
		  
		  //loop over space distance
		  for(int d=0;d<Dmax;d++)
		    {
		      //find B
		      int B=A;
		      for(int i=0;i<=d;i++) B=loclx_neighup[B][mu];

		      // |
		      // |__   part 
		      su3 p1;
		      unsafe_su3_dag_prod_su3(p1,tline[A*Tmax+t],xline[A*Dmax+d]);
		      
		      // |  |
		      // |__|  part 
		      su3 p2;
		      unsafe_su3_prod_su3(p2,p1,tline[B*Tmax+t]);
		      
		      //close with horizontal line DC
		      U[d][t]+=real_part_of_trace_su3_prod_su3_dag(p2,xline[D*Dmax+d]);
		    }
		}
	    }
	}
    
      //normalize and print
      for(int d=0;d<Dmax;d++)
	{
	  for(int t=0;t<Tmax;t++)  master_fprintf(fout,"%d  %.18lg\n",t+1,U[d][t]/glb_vol/3/3);
	  master_fprintf(fout,"&\n");
	}
    } 
  
  //close out
  if(rank==0) fclose(fout);
      
  nissa_free(xline);
  nissa_free(tline);
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
  
  //check that we are not working in parallel
  if(nissa_nranks>1) crash("scalar only code");
  
  char conf_path[1024];
  read_str_str("GaugePath",conf_path,1024);
  char out_path[1024];
  read_str_str("OutPath",out_path,1024);
  double hyp_alpha0,hyp_alpha1,hyp_alpha2;
  read_str_double("HypAlpha0",&hyp_alpha0);
  read_str_double("HypAlpha1",&hyp_alpha1);
  read_str_double("HypAlpha2",&hyp_alpha2);
  double ape_alpha;
  read_str_double("ApeAlpha",&ape_alpha);
  int nlev_spat_smear,*niter_lev_spat_smear;
  read_list_of_ints("ApeNiter",&nlev_spat_smear,&niter_lev_spat_smear);
  int Tmax,Dmax;
  read_str_int("Dmax",&Dmax);
  read_str_int("Tmax",&Tmax);
  
  /////////////////////////////////
  
  //read conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  read_ildg_gauge_conf(conf,conf_path);
  communicate_lx_quad_su3_borders(conf);
  master_printf(" read conf plaq: %.18g\n",global_plaquette_lx_conf(conf));
  
  compute_static_quark_pair_potential(out_path,conf,hyp_alpha0,hyp_alpha1,hyp_alpha2,ape_alpha,nlev_spat_smear,niter_lev_spat_smear,Tmax,Dmax);

  ///////////////////////////////////
  
  close_input();
  
  nissa_free(conf);
  
  close_nissa();
  
  return 0;
}
