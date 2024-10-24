#include <nissa.hpp>

using namespace nissa;

void compute_gaussianity_pars(double* x,color* source,int maxpow,coords_t* source_pos)
{
  //reset local pows
  double locx[glbSize[0]][maxpow];
  for(int t=0;t<glbSize[0];t++)
    for(int ipow=0;ipow<maxpow;ipow++)
      locx[t][ipow]=0.0;
  
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    {
      int t=glbCoordOfLoclx[ivol][0];
      
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
	      int xmu=(glbCoordOfLoclx[ivol][mu]-source_pos[t][mu]+glbSize[mu])%glbSize[mu];
	      if(xmu>=glbSize[mu]/2) xmu-=glbSize[mu];
	      xpow+=pow(xmu,ipow*2);
	    }
	  
	  locx[t][ipow]+=n*xpow;
	}
    }
  NISSA_PARALLEL_LOOP_END;
  THREAD_BARRIER();
  
  //reduce
  for(int t=0;t<glbSize[0];t++)
    for(int ipow=0;ipow<maxpow;ipow++)
      non_loc_reduce(&x[t*maxpow+ipow],&locx[t][ipow]);
}

//get average and error of gaussianity pars
void process_gaussianity(double *a,double *e,double *x,int maxpow)
{
  //reset summ and errors
  for(int ipow=0;ipow<maxpow;ipow++)
    a[ipow]=e[ipow]=0.0;
  
  for(int t=0;t<glbSize[0];t++)
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
      a[ipow]/=glbSize[0];
      e[ipow]/=glbSize[0];
      e[ipow]-=sqr(a[ipow]);
      e[ipow]=sqrt(fabs(e[ipow])/glbSize[0]);
    }
}

//according to Bali
double expected_radius(double kappa,int nlevels,double plaq)
{
  kappa*=pow(plaq,0.25);
  return sqrt(nlevels*kappa/(1+2*(NDIM-1)*kappa));
}

//hold a tern to keep density
struct dens_t
{
  int n;
  double s;
  dens_t() : n(0),s(0.0) {}
};

typedef std::map<int,dens_t> mapdens_t;

//compute the density distribution
void compute_density(FILE *fout,color *source,coords_t *source_pos)
{
  mapdens_t density[glbSize[0]];
  
  NISSA_LOC_VOL_LOOP(ivol)
    {
      int t=glbCoordOfLoclx[ivol][0];
      
      //get site norm
      double n=0.0;
      for(int ic=0;ic<NCOL;ic++)
	for(int ri=0;ri<2;ri++)
	  n+=sqr(source[ivol][ic][ri]);
      
      //compute distance
      int rho=0.0;
      for(int mu=1;mu<NDIM;mu++)
	{
	  int xmu=(glbCoordOfLoclx[ivol][mu]-source_pos[t][mu]+glbSize[mu])%glbSize[mu];
	  if(xmu>=glbSize[mu]/2) xmu-=glbSize[mu];
	  rho+=sqr(xmu);
	}
      
      //increment
      dens_t &it=density[t][rho];
      it.n++;
      it.s+=n;
    }
  
  //reduce and print
  master_fprintf(fout," NDists %d\n",(int)density[0].size());

  for(int t=0;t<glbSize[0];t++)
    {
      master_fprintf(fout," t %d\n",t);
      
      for(mapdens_t::iterator it=density[t].begin();it!=density[t].end();it++)
	{
	  int r2=it->first;
	  dens_t d=it->second;
	  non_loc_reduce(&d.n);
	  non_loc_reduce(&d.s);
	  double n=d.n;
	  double s=d.s/n;
	  
	  master_fprintf(fout,"%d %lg\n",r2,s);
	}
      master_fprintf(fout,"\n");
    }
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
  quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol,quad_su3);
  read_ildg_gauge_conf(conf,conf_path);
  
  //read APE smearing pars
  int ape_smearing_niters;
  double ape_smearing_alpha;
  read_str_double("ApeSmearingAlpha",&ape_smearing_alpha);
  read_str_int("ApeSmearingNiters",&ape_smearing_niters);
  ape_spatial_smear_conf(conf,conf,ape_smearing_alpha,ape_smearing_niters);
  
  //read Gaussian smearing pars
  int ncoeffs;
  double kappa;
  read_str_double("Kappa",&kappa);
  read_str_int("NCoeffs",&ncoeffs);
  std::map<size_t,double> coeffs;
  for(size_t iCoeff=0;iCoeff<ncoeffs;iCoeff++)
    {
      double w;
      read_double(&w);
      int n;
      read_int(&n);
      
      coeffs[n]+=w;
    }
  
  //output file
  char out_path[1024];
  read_str_str("Output",out_path,1024);
  FILE *fout=open_file(out_path,"w");
  master_fprintf(fout,"T %d\n",glbSize[0]);
  master_fprintf(fout,"Kappa %lg\n",kappa);
  for(const auto [n,w] : coeffs)
    master_fprintf(fout,"%lg*H^%zu",w,n);
  master_fprintf(fout,"\n");
  
  //print spatial plaquette
  double plaqs[2];
  global_plaquette_lx_conf(plaqs,conf);
  master_fprintf(fout,"TimePlaquette %16.16lg\n",plaqs[0]);
  master_fprintf(fout,"SpatPlaquette %16.16lg\n",plaqs[1]);
  master_fprintf(fout,"\n");
  
  //set the source
  color *source=nissa_malloc("source",locVol+bord_vol,color);
  vector_reset(source);
  coords_t source_pos[glbSize[0]];
  for(int t=0;t<glbSize[0];t++)
    {
      //generate coords and fix t
      source_pos[t]=generate_random_coord();
      source_pos[t][0]=t;
      
      //get loclx and rank
      int l,r;
      get_loclx_and_rank_of_coord(l,r,source_pos[t]);
      
      //put the source only if on correct rank
      if(rank==r) source[l][0][0]=1;
    }
    
  color *tot=nissa_malloc("tot",locVol+bord_vol,color);
  vector_reset(tot);
  size_t p=0;
  for(const auto [n,w] : coeffs)
    {
      gaussian_smearing(source,source,conf,kappa,n-p);
      p=n;
      double_vector_summassign_double_vector_prod_double(tot[0][0],source[0][0],w,locVol*sizeof(color)/sizeof(double));
    }
  
  // //compute gaussianity
  // int maxpow=2;
  // double x[maxpow*glbSize[0]];
  // compute_gaussianity_pars(x,source,maxpow,source_pos);
  
  // //take averages
  // double a[maxpow],e[maxpow];
  // process_gaussianity(a,e,x,maxpow);
  
  // //write
  // master_printf(" - average norm:   %lg +- %lg\n",a[0],e[0]);
  // master_printf(" - average radius: %lg +- %lg\n",a[1],e[1]);
  // master_printf("   expected:       %lg\n",expected_radius(kappa,ilev,plaqs[1]));
  // master_printf("\n");
      
  //     master_fprintf(fout," Smearlevel %d\n",ilev);
  compute_density(fout,tot,source_pos);
  
  //smear
  
  close_file(fout);
  
  nissa_free(source);
  nissa_free(conf);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}

