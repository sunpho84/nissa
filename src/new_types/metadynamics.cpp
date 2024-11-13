#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "threads/threads.hpp"
#include "io/input.hpp"

#include "metadynamics.hpp"

namespace nissa
{
  //print all pars
  std::string meta_pars_t::get_str(const bool& full) const
  {
    std::ostringstream os;
    os<<"After\t\t=\t"<<after<<"\n";
    os<<"Each\t\t=\t"<<each<<"\n";
    os<<"Coeff\t\t=\t"<<coeff<<"\n";
    os<<"Width\t\t=\t"<<width<<"\n";
    os<<"Barr\t\t=\t"<<barr<<"\n";
    os<<"ForceOut\t=\t"<<force_out<<"\n";
    os<<"WellTempering\t=\t"<<well_tempering<<"\n";
    os<<"Bend\t\t=\t"<<bend<<"\n";
    
    return os.str();
  }
  
  //update the history-dependent potential
  void meta_pars_t::update(const int& isweep,
			   const double& Q)
  {
    if(isweep>=after && (isweep-after)%each==0)
      {
	const int igrid=floor(Q/width)+ngrid/2.0;
	
	double alpha=Q/width;
	alpha=alpha-floor(alpha);
	
	if(igrid>=0 and igrid<=ngrid)
	  grid[igrid]+=(1-alpha)*coeff;
	
	if(igrid+1>=0 and igrid+1<=ngrid)
	  grid[igrid+1]+=alpha*coeff;
      }
  }
  
  //compute the derivative of the potential
  double meta_pars_t::compute_pot_der(const double& x) const
  {
    //take igrid
    const int igrid=floor((x+barr)/width);
    
    //inside the barriers
    if(igrid>=0 and igrid<ngrid)
      return (grid[igrid+1]-grid[igrid])/width;
    else
      if(igrid<0)
	return -force_out*(-x-barr);
      else
	return +force_out*(+x-barr);
  }
  
  //compute the potential using past history
  double meta_pars_t::compute_pot(const double& x) const
  {
    //take igrid
    const int igrid=floor((x+barr)/width);
    
    //inside the barriers
    if(igrid>=0 and igrid<ngrid)
      {
	const double x0=igrid*width-barr;
	const double m=(grid[igrid+1]-grid[igrid])/width;
	const double q=grid[igrid]-m*x0;
	
	return q+m*x;
      }
    else
      if(igrid<0)
	return force_out*sqr(-x-barr)/2+grid[0];
      else
	return force_out*sqr(+x-barr)/2+grid[ngrid];
  }
  
  //write
  void meta_pars_t::save(const char *path) const
  {
    FILE *fout=open_file(path,"w");
    
    for(int i=0;i<=ngrid;i++)
      nissa::master_fprintf(fout,"%lg %16.16lg\n",-barr+i*width,grid[i]);
    
    close_file(fout);
  }
  
  //read
  void meta_pars_t::load(const char *path)
  {
    //to be sure, resize
    grid.resize(ngrid+1);
    
    FILE *fin=open_file(path,"r");
    if(rank==0)
      for(int igrid=0;igrid<=ngrid;igrid++)
	{
	  double xread;
	  int rc=fscanf(fin,"%lg %lg",&xread,&grid[igrid]);
	  if(rc!=2) crash("reading line %d of \"%s\"",igrid,path);
	  int jgrid=floor((xread+barr+width/2)/width);
	  if(igrid!=jgrid) crash("found %d (%lg) when expecting %d",jgrid,xread,igrid);
	}
    close_file(fin);
    
    //broadcast
    for(int igrid=0;igrid<=ngrid;igrid++) MPI_Bcast(&grid[igrid],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  
  //draw the chronological force
  void meta_pars_t::draw_force(const char *force_path) const
  {
    const double x_min=-barr*1.1;
    const double x_max=+barr*1.1;
    const double x_diff=x_max-x_min;
    int n=ceil(x_diff/width*10);
    if(n==0) n=1;
    double dx=x_diff/n;
    
    //compute
    double *xy=new double[n+1];
    double *xz=new double[n+1];
    for(int i=0;i<=n;i++)
      {
	xy[i]=compute_pot_der(x_min+i*dx);
	xz[i]=(compute_pot(x_min+i*dx+dx/10)-compute_pot(x_min+i*dx-dx/10))/(dx/5);
      }
    
    //write
    FILE *fout=open_file(force_path,"w");
    for(int i=0;i<=n;i++) nissa::master_fprintf(fout,"%16.16lg %16.16lg\n",x_min+i*dx,xy[i]);
    nissa::master_fprintf(fout,"&\n");
    for(int i=0;i<=n;i++) nissa::master_fprintf(fout,"%16.16lg %16.16lg\n",x_min+i*dx,xz[i]);
    close_file(fout);
    
    delete[] xy;
    delete[] xz;
  }
  
  //initialize
  void meta_pars_t::init()
  {
    ngrid=(2*barr+width/2)/width;
    grid.resize(ngrid+1);
    for(int igrid=0;igrid<=ngrid;igrid++)
      grid[igrid]=0;
  }
  
  //read from a file all the parameters
  void meta_pars_t::read_pars()
  {
    read_str_int("MetaAfter",&after);
    read_str_int("MetaEach",&each);
    read_str_double("MetaCoeff",&coeff);
    read_str_double("MetaWidth",&width);
    read_str_double("MetaBarr",&barr);
    read_str_double("MetaForceOut",&force_out);
    read_str_double("MetaBend",&bend);
    read_str_double("MetaWellTempering",&well_tempering);
    
    init();
  }
}
