#include <nissa.hpp>

using namespace nissa;

THREADABLE_FUNCTION_4ARG(compute_gaussianity, double*,x, color*,source, int,maxpow, coords*,source_pos)
{
  GET_THREAD_ID();
  
  //reset local pows
  double locx[glb_size[0]][maxpow];
  for(int t=0;t<glb_size[0];t++)
    for(int ipow=0;ipow<maxpow;ipow++)
      locx[t][ipow]=0.0;
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      int t=glb_coord_of_loclx[ivol][0];
      
      //get site norm
      double n=0.0;
      for(int ic=0;ic<NCOL;ic++)
	for(int ri=0;ri<2;ri++)
	  n+=sqr(source[ivol][ic][ri]);
      
      //loop over all powers to be computed
      for(int ipow=0;ipow<maxpow;ipow++)
	{
	  //compute distance
	  double xpow=0.0;
	  for(int mu=1;mu<NDIM;mu++)
	    {
	      int xmu=(glb_coord_of_loclx[ivol][mu]-source_pos[t][mu]+glb_size[mu])%glb_size[mu];
	      if(xmu>=glb_size[mu]/2) xmu-=glb_size[mu];
	      xpow+=pow(xmu,ipow*2);
	    }
	  
	  locx[t][ipow]+=n*xpow;
	}
    }
  THREAD_BARRIER();
  
  //reduce
  for(int t=0;t<glb_size[0];t++)
    for(int ipow=0;ipow<maxpow;ipow++)
      x[t*maxpow+ipow]=glb_reduce_double(locx[t][ipow]);
}
THREADABLE_FUNCTION_END

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //open input file
  open_input(arg[1]);
  
  //init the grid
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //read conf
  char conf_path[1024];
  read_str_str("Conf",conf_path,1024);
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  read_ildg_gauge_conf(conf,conf_path);
  
  //read APE smearing pars
  int ape_smearing_niters;
  double ape_smearing_alpha;
  read_str_double("ApeSmearingAlpha",&ape_smearing_alpha);
  read_str_int("ApeSmearingNiters",&ape_smearing_niters);
  ape_spatial_smear_conf(conf,conf,ape_smearing_alpha,ape_smearing_niters);
  
  //read Gaussian smearing pars
  int nlevels;
  double kappa;
  read_str_double("Kappa",&kappa);
  read_str_int("NLevels",&nlevels);
  
  //set the source
  color *source=nissa_malloc("source",loc_vol+bord_vol,color);
  vector_reset(source);
  coords source_pos[glb_size[0]];
  for(int t=0;t<glb_size[0];t++)
    {
      //generate coords and fix t
      generate_random_coord(source_pos[t]);
      source_pos[t][0]=t;
      
      //get loclx and rank
      int l,r;
      get_loclx_and_rank_of_coord(&l,&r,source_pos[t]);
      
      if(rank==r) source[l][0][0]=1;
    }
    
  for(int ilev=0;ilev<=nlevels;ilev++)
    {
      compute_gaussianity(x,source,maxpow,source_pos);
      
      int maxpow=4;
      double x[maxpow*glb_size[0]];
      
      master_printf("smear %d \n",ilev);
      for(int t=0;t<glb_size[0];t++)
	{
	  master_printf("%d\t%lg\t",t,x[t*maxpow+0]);
	  for(int ipow=1;ipow<maxpow;ipow++) master_printf("%lg\t",x[t*maxpow+ipow]/x[t*maxpow+0]);
	  master_printf("\n");
	}
      master_printf("\n");
      
      if(ilev<nlevels) gaussian_smearing(source,source,conf,kappa,1);
    }
  
  nissa_free(source);
  nissa_free(conf);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}

