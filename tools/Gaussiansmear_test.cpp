#include <nissa.hpp>

using namespace nissa;

THREADABLE_FUNCTION_4ARG(compute_gaussianity_pars, double*,x, color*,source, int,maxpow, coords*,source_pos)
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

//get average and error of gaussianity pars
void process_gaussianity(double *a,double *e,double *x,int maxpow)
{
  //reset summ and errors
  for(int ipow=0;ipow<maxpow;ipow++)
    a[ipow]=e[ipow]=0.0;
  
  for(int t=0;t<glb_size[0];t++)
    for(int ipow=0;ipow<maxpow;ipow++)
      {
	double s=0;
	if(ipow==0) s=x[t*maxpow+0];
	if(ipow==1) s=sqrt(x[t*maxpow+1]/x[t*maxpow+0]);
	
	//increment
	a[ipow]+=s;
	e[ipow]+=s*s;
      }
  
  //build averages and errors
  for(int ipow=0;ipow<maxpow;ipow++)
    {
      a[ipow]/=glb_size[0];
      e[ipow]/=glb_size[0];
      e[ipow]-=sqr(a[ipow]);
      e[ipow]=sqrt(fabs(e[ipow]));
    }
}

//according to Bali
double expected_radius(double kappa,int nlevels,double plaq)
{
  kappa*=pow(plaq,0.25);
  return sqrt(nlevels*kappa/(1+2*(NDIM-1)*kappa));
}

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
  
  //print spatial plaquette
  double plaqs[2];
  global_plaquette_lx_conf(plaqs,conf);
  master_printf("Plaquettes: %16.16lg, %16.16lg\n",plaqs[0],plaqs[1]);
  
  //read Gaussian smearing pars
  int nlevels,meas_each;
  double kappa;
  read_str_double("Kappa",&kappa);
  read_str_int("NLevels",&nlevels);
  read_str_int("MeasEach",&meas_each);
  
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
      
      //put the source only if on correct rank
      if(rank==r) source[l][0][0]=1;
    }
    
  for(int ilev=0;ilev<=nlevels;ilev+=meas_each)
    {
      //compute gaussianity
      int maxpow=2;
      double x[maxpow*glb_size[0]];
      compute_gaussianity_pars(x,source,maxpow,source_pos);
      
      //take averages
      double a[maxpow],e[maxpow];
      process_gaussianity(a,e,x,maxpow);
      
      //write
      master_printf("Smearing level %d\n",ilev);
      master_printf(" - average norm:   %lg +- %lg\n",a[0],e[0]);
      master_printf(" - average radius: %lg +- %lg\n",a[1],e[1]);
      master_printf("   expected:       %lg\n",expected_radius(kappa,ilev,plaqs[1]));
      master_printf("\n");
      
      //smear
      if(ilev<nlevels) gaussian_smearing(source,source,conf,kappa,meas_each);
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

