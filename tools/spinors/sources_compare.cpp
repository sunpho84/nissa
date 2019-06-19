#include <nissa.hpp>

using namespace nissa;

//! integrate
double integrate_corr(const std::vector<double> &c,const std::vector<double> &d,double powd=0.0)
{
  double out;
  out=0.0;
  for(size_t iel=0;iel<c.size()-1;iel++)
    {
      double x1=d[iel+1];
      double x0=d[iel];
      double dx=x1-x0;
      double y1=c[iel+1]*pow(x1,powd);
      double y0=c[iel]*pow(x0,powd);
      
      out+=dx*(y1+y0)/2.0;
    }
  
  return 2.0*out;
}

void in_main(int narg,char **arg)
{
  GET_THREAD_ID();
  
  if(narg<6) crash("use: %s L T file_in file_out output_path tag",arg[0]);
  
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  char *source_path=arg[3];
  char *smeared_source_path=arg[4];
  char *output_path=arg[5];
  char *tag=arg[6];
  
  //Init the MPI grid
  init_grid(T,L);
  
  ///////////////////////////////////////////
  
  spincolor *source=nissa_malloc("source",loc_vol,spincolor);
  spincolor *smeared_source=nissa_malloc("smeared_source",loc_vol,spincolor);
  read_real_vector(source,source_path,tag);
  read_real_vector(smeared_source,smeared_source_path,tag);
  
  int iglb_max=0;
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    if(spincolor_norm2(source[ivol])>1e-10)
      iglb_max=glblx_of_loclx[ivol];
  MPI_Allreduce(MPI_IN_PLACE,&iglb_max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  coords g;
  glb_coord_of_glblx(g,iglb_max);
  master_printf("Source location: %d %d %d %d\n",g[0],g[1],g[2],g[3]);
  
  //check the norm
  double source_norm=double_vector_glb_norm2(source,loc_vol);
  if(fabs(source_norm-1.0)>1e-10)
    crash("Norm %lg, needs to be 1, excess of %lg",source_norm,source_norm-1.0);
  
  std::map<int,std::pair<double,int>> rho;
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    if(glb_coord_of_loclx[ivol][0]==g[0])
      {
	int r2=0;
	for(int mu=1;mu<NDIM;mu++)
	  {
	    int c=(glb_size[mu]+glb_coord_of_loclx[ivol][mu]-g[mu])%glb_size[mu];
	    r2+=sqr(std::min(c,glb_size[mu]-c));
	  }
	rho[r2].first+=spincolor_norm2(smeared_source[ivol]);
	rho[r2].second++;
      }
  THREAD_BARRIER();
  
  for(int i=1;i<nranks;i++)
    {
      if(rank==i)
	{
	  int nel=rho.size();
	  MPI_Send(&nel,1,MPI_INT,0,909,MPI_COMM_WORLD);
	  for(auto &r : rho)
	    {
	      int r2=r.first;
	      double p=r.second.first;
	      int n=r.second.second;
	      MPI_Send(&r2,1,MPI_INT,0,910,MPI_COMM_WORLD);
	      MPI_Send(&p,1,MPI_DOUBLE,0,911,MPI_COMM_WORLD);
	      MPI_Send(&n,1,MPI_INT,0,912,MPI_COMM_WORLD);
	    }
	}
      
      if(rank==0)
	{
	  verbosity_lv2_master_printf("Communicating with %d\n",i);
	  
	  int nel;
	  MPI_Recv(&nel,1,MPI_INT,i,909,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  verbosity_lv2_master_printf("Needs to receive %d elements from %d\n",nel,i);
	  for(int iel=0;iel<nel;iel++)
	    {
	      verbosity_lv2_master_printf("Receiving element %d from %d\n",iel,i);
	      
	      int r2;
	      double p;
	      int n;
	      MPI_Recv(&r2,1,MPI_INT,i,910,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      MPI_Recv(&p,1,MPI_DOUBLE,i,911,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      MPI_Recv(&n,1,MPI_INT,i,912,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      
	      rho[r2].first+=p;
	      rho[r2].second+=n;
	    }
	}
      
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  if(rank==0)
    {
      FILE *fout=open_file(output_path,"w");
      std::vector<double> c,d;
      for(auto &p : rho)
	{
	  double r2=p.first;
	  double w=p.second.first;
	  int n=p.second.second;
	  master_fprintf(fout,"%lg" "\t" "%lg" "\t" "%d" "\n",sqrt(r2),w,n);
	  
	  d.push_back(sqrt(r2));
	  c.push_back(w/n);
	}
      
      double x2=integrate_corr(c,d,2);
      double n=integrate_corr(c,d,0.0);
      master_printf("Radius: %lg\n",sqrt(x2/n));
      
      close_file(fout);
    }
  
  nissa_free(source);
  nissa_free(smeared_source);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  close_nissa();
  
  return 0;
}
