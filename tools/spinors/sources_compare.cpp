#include <nissa.hpp>

using namespace nissa;

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

  fft4d(source,source,+1,false);
  fft4d(smeared_source,smeared_source,+1,false);
  
  complex *prod=nissa_malloc("prod",loc_vol,complex);
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    spincolor_scalar_prod(prod[ivol],source[ivol],smeared_source[ivol]);
  
  fft4d(prod,prod,-1, true);
  
  std::map<int,std::pair<double,int>> rho;
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      int r2=0;
      for(int mu=0;mu<NDIM;mu++)
	{
	  int c=glb_coord_of_loclx[ivol][mu];
	  r2+=sqr(std::min(c,glb_size[mu]/2-c));
	}
      rho[r2].first+=prod[ivol][RE];
      rho[r2].second++;
    }
  THREAD_BARRIER();
  
  for(int i=1;i<nranks;i++)
    {
      if(i==rank)
	{
	  int nel=rho.size();
	  MPI_Send(&nel,1,MPI_INT,0,909,MPI_COMM_WORLD);
	  for(auto &r : rho)
	    {
	      MPI_Send(&r.first,1,MPI_INT,0,910,MPI_COMM_WORLD);
	      MPI_Send(&r.second.first,1,MPI_DOUBLE,0,911,MPI_COMM_WORLD);
	    }
	}
      
      if(i==0)
	{
	  int nel;
	  MPI_Recv(&nel,1,MPI_INT,1,909, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  for(int i=0;i<nel;i++)
	    {
	      int r2;
	      double p;
	      MPI_Recv(&r2,1,MPI_INT,0,910,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      MPI_Recv(&p,1,MPI_DOUBLE,0,911,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      
	      rho[r2].first+=p;
	      rho[r2].second++;
	    }
	}
      
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  FILE *fout=open_file(output_path,"w");
  for(auto &r : rho)
    master_fprintf(fout,"%lg" "\t" "%lg" "\n",sqrt(r.first),r.second.first/r.second.second);
  close_file(fout);
  
  nissa_free(prod);
  nissa_free(source);
  nissa_free(smeared_source);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  close_nissa();
  
  return 0;
}
