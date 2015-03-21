#include <nissa.hpp>

using namespace nissa;

double tot_prog_time=0;
char path_in[200];
char path_out[200];

int nr;
int nconfs,nterm,njacks,each,clust_size;

double *in_buffer,*out_buffer;

//read the input
void init_clusterize(const char *path)
{
  open_input(path);
  
  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid 
  init_grid(T,L); 
  //Pars
  read_str_int("Nr",&nr);
  //Input and output paths
  read_str_str("PathIn",path_in,200);
  read_str_str("PathOut",path_out,200);
  //Statistics
  read_str_int("NConfs",&nconfs);
  read_str_int("Each",&each);
  read_str_int("NTerm",&nterm);
  read_str_int("NJacks",&njacks);
  
  close_input();
  
  //////////////////////////////
  
  //find cluster sizes
  clust_size=(nconfs-nterm)/each/njacks;
  nconfs=nterm+clust_size*each*njacks;
  master_printf("Adapted nconfs to %d\n",nconfs);
  if(njacks<=1) crash("cannot use njacks %d (at least 2)");
  if(each==0) crash("cannot use each 0");
  if(clust_size==0) crash("cannot use clust_size 0");
  
  master_printf("Nconfs to consider: %d\n",nconfs);
  master_printf("Cluster size: %d\n",clust_size);
  
  //allocate in buffer
  in_buffer=nissa_malloc("in_buffer",loc_vol,double);
  out_buffer=nissa_malloc("out_bufffer",loc_vol*(njacks+1)*nr,double);
}

//add it to the cluster
THREADABLE_FUNCTION_4ARG(add_cluster, double*,out_buffer, double*,in_buffer, int,iconf, int,r)
{
  GET_THREAD_ID();
  
  int iclust;
  if(iconf==0) iclust=njacks;
  else
    {
      iclust=(iconf-nterm)/each/clust_size;
      if(iclust>=njacks) crash("iclust: %d >= njacks: %d, iconf: %d",iclust,njacks,iconf);
    }
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    out_buffer[iclust+(njacks+1)*(ivol+loc_vol*r)]+=in_buffer[ivol];
  THREAD_BARRIER();
}
THREADABLE_FUNCTION_END

//load the correlation funcions
void load_data(const char *path)
{
  ILDG_File file=ILDG_File_open_for_read(path);
  
  //reset output
  vector_reset(out_buffer);
  
  while(!ILDG_File_reached_EOF(file))
    {
      //read header
      ILDG_header header;
      header=ILDG_File_get_next_record_header(file);
      
      //parse the header
      master_printf("Found header: %s\n",header.type);
      int iconf,r;
      int rc=sscanf(header.type,"%d_%d",&iconf,&r);
      if(rc!=2) crash("returned %d",rc); 
      if(r>=nr) crash("loaded r=%d while nr=%d",r,nr);
      
      //prompt found conf and r
      verbosity_lv2_master_printf("Conf: %d, r: %d\n",iconf,r);
      
      //if we passed termalization we read otherwise we skip it
      if((iconf<nterm || iconf>=nconfs) && (iconf!=0))
	{
	  ILDG_File_skip_record(file,header);
	  //discard checksum
	  checksum check_read;
	  ILDG_File_read_checksum(check_read,file);
	}
      else
	{
	  //read data
	  read_real_vector(in_buffer,file,header,1);
	  add_cluster(out_buffer,in_buffer,iconf,r);
	}
    }
  
  ILDG_File_close(file);
}

//clusterize the single jackninfe
void clusterize(double *data)
{
  double ze=data[njacks];
  
  data[njacks]=0;
  for(int ijack=0;ijack<njacks;ijack++) data[njacks]+=data[ijack];
  for(int ijack=0;ijack<njacks;ijack++) data[ijack]=(data[njacks]-data[ijack])/((njacks-1)*clust_size);
  data[njacks]/=njacks*clust_size;
  
  for(int ijack=0;ijack<=njacks;ijack++) data[ijack]=ze-data[ijack];
}

//clusterize
THREADABLE_FUNCTION_0ARG(clusterize)
{
  GET_THREAD_ID();

  for(int r=0;r<nr;r++)
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      clusterize(out_buffer+0+(njacks+1)*(ivol+loc_vol*r));
  THREAD_BARRIER();
}
THREADABLE_FUNCTION_END

//save the clusterized data
void save_data(const char *path)
{
  FILE *fout=open_file(path,"w");
  
  if(!little_endian) change_endianness(out_buffer,out_buffer,loc_vol*nr*(njacks+1));
  int rc=fwrite(out_buffer,sizeof(double),loc_vol*nr*(njacks+1),fout);
  if(rc!=loc_vol*nr*(njacks+1)) crash("returned %d instead of %d",rc,loc_vol*nr*(njacks+1));
     
  close_file(fout);
}

//close it
void close_clusterize()
{
  nissa_free(in_buffer);
  nissa_free(out_buffer);
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //check arguments
  if(narg<2) crash("Use %s input",arg[0]);
  init_clusterize(arg[1]);
  
  load_data(path_in);
  clusterize();
  save_data(path_out);
  
  tot_prog_time+=take_time();
  close_clusterize();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
    
  return 0;
}
