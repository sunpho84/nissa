#include <set>

#include <nissa.hpp>

using namespace nissa;

//hold a tern to keep density
struct dens_t
{
  int64_t n{};
  double s{};
  double s2{};
};

typedef std::map<int,dens_t> mapdens_t;


//Compute the density distribution
void compute_density(FILE *fout,
		     color* vec)
{
  mapdens_t density;
  
  static std::set<int> dists;
  if(dists.size()==0)
    for(int64_t gvol=0;gvol<glbVol/glbSize[0];gvol++)
      {
	coords_t g=glb_coord_of_glblx(gvol);
	int rho=0;
	for(int mu=1;mu<NDIM;mu++)
	  {
	    if(g[mu]>=glbSize[mu]/2)
	      g[mu]=glbSize[mu]-g[mu];
	    rho+=sqr(g[mu]);
	  }
	dists.insert(rho);
      }
  
  for(int64_t ivol=0;ivol<locVol;ivol++)
    {
      //get site norm
      double s2=0.0;
      for(int ic=0;ic<NCOL;ic++)
	for(int ri=0;ri<2;ri++)
	  s2+=sqr(vec[ivol][ic][ri]);
      
      //compute distance
      int rho=0.0;
      for(int mu=1;mu<NDIM;mu++)
	{
	  int xmu=glbCoordOfLoclx[ivol][mu];
	  if(xmu>=glbSize[mu]/2) xmu-=glbSize[mu];
	  rho+=sqr(xmu);
	}
      
      //increment
      dens_t &it=density[rho];
      it.n++;
      it.s+=s2;
      it.s2+=s2*s2;
    }
  
  //reduce and print
  // master_fprintf(fout," NDists %d\n",(int)density.size());
  
  for(auto& dist : dists)
    {
      dens_t& d=density[dist];
      
      non_loc_reduce(&d.n);
      non_loc_reduce(&d.s);
      non_loc_reduce(&d.s2);
      
      const int64_t n=d.n;
      const double s=d.s/n;
      const double s2=d.s2/n-s*s;
      
      master_fprintf(fout,"%d %lg %lg\n",dist,s,sqrt(s2/n));
    }
  
  master_fprintf(fout,"\n");
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
  double kappa;
  read_str_double("Kappa",&kappa);
  size_t nPoly;
  {
    int n;
    read_str_int("NPoly",&n);
    nPoly=n;
  }
  std::vector<std::map<size_t,double>> coeffs(nPoly);
  
  std::vector<color*> tot(nPoly);
  for(auto& t : tot)
    {
      t=nissa_malloc("tot",locVol,color);
      vector_reset(t);
    }
  
  std::set<size_t> nList;
  std::vector<std::string> outPaths(nPoly);
  for(size_t iPoly=0;iPoly<nPoly;iPoly++)
    {
      char outPath[1024];
      read_str_str("Output",outPath,1024);
      outPaths[iPoly]=outPath;
      
      int nCoeffs;
      read_str_int("NCoeffs",&nCoeffs);
      for(int iCoeff=0;iCoeff<nCoeffs;iCoeff++)
	{
	  double w;
	  read_double(&w);
	  int n;
	  read_int(&n);
	  
	  coeffs[iPoly][n]+=w;
	  
	  nList.insert(n);
	}
    }
  
  //print spatial plaquette
  double plaqs[2];
  global_plaquette_lx_conf(plaqs,conf);
  master_printf("TimePlaquette %16.16lg\n",plaqs[0]);
  master_printf("SpatPlaquette %16.16lg\n",plaqs[1]);
  
  //set the source
  color *source=nissa_malloc("source",locVol+bord_vol,color);
  vector_reset(source);
  
  for(coords_t c{};c[0]<glbSize[0];c[0]++)
    {
      int ivol;
      int r;
      get_loclx_and_rank_of_coord(ivol,r,c);
      
      if(r==rank)
	source[ivol][0][0]=1;
    }
  set_borders_invalid(source);
  
  size_t p=0;
  for(const auto n : nList)
    {
      gaussian_smearing(source,source,conf,kappa,n-p);
      
      for(size_t iPoly=0;iPoly<nPoly;iPoly++)
	if(const auto nw=coeffs[iPoly].find(n);nw!=coeffs[iPoly].end())
	  {
	    const auto& [n,w]=*nw;
	    master_printf("Adding %lg*H^%zu (computed with %zu new steps) to poly %zu\n",w,n,n-p,iPoly);
	    
	    double_vector_summassign_double_vector_prod_double(&tot[iPoly][0][0][0],&source[0][0][0],w,locVol*sizeof(color)/sizeof(double));
	  }
      
      p=n;
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
  //output file
  for(size_t iPoly=0;iPoly<nPoly;iPoly++)
    {
      FILE* f=open_file(outPaths[iPoly],"w");
      
      // master_fprintf(f,"T %d\n",glbSize[0]);
      // master_fprintf(f,"Kappa %lg\n",kappa);
      // for(const auto [n,w] : coeffs[iPoly])
      // 	master_fprintf(f,"%+lg*H^%zu",w,n);
      // master_fprintf(f,"\n");
      
      master_fprintf(f,"\n");
      
      compute_density(f,tot[iPoly]);
      
      close_file(f);
    }
  
  for(size_t iPoly=0;iPoly<nPoly;iPoly++)
    nissa_free(tot[iPoly]);
  nissa_free(source);
  nissa_free(conf);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}

