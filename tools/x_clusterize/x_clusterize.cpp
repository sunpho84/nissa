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
  initGrid(T,L); 
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
  MASTER_PRINTF("Adapted nconfs to %d\n",nconfs);
  if(njacks<=1) CRASH("cannot use njacks %d (at least 2)");
  if(each==0) CRASH("cannot use each 0");
  if(clust_size==0) CRASH("cannot use clust_size 0");
  
  MASTER_PRINTF("Nconfs to consider: %d\n",nconfs);
  MASTER_PRINTF("Cluster size: %d\n",clust_size);
  
  //allocate in buffer
  in_buffer=nissa_malloc("in_buffer",loc_vol,double);
  out_buffer=nissa_malloc("out_bufffer",loc_vol*(njacks+1)*nr,double);
}

//add it to the cluster
void add_cluster(double* out_buffer,double* in_buffer,int iconf,int r)
{
  
  int iclust;
  if(iconf==0) iclust=njacks;
  else
    {
      iclust=(iconf-nterm)/each/clust_size;
      if(iclust>=njacks) CRASH("iclust: %d >= njacks: %d, iconf: %d",iclust,njacks,iconf);
    }
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    out_buffer[iclust+(njacks+1)*(ivol+loc_vol*r)]+=in_buffer[ivol];
  NISSA_PARALLEL_LOOP_END;
  THREAD_BARRIER();
}

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
      MASTER_PRINTF("Found header: %s\n",header.type);
      int iconf,r;
      int rc=sscanf(header.type,"%d_%d",&iconf,&r);
      if(rc!=2) CRASH("returned %d",rc); 
      if(r>=nr) CRASH("loaded r=%d while nr=%d",r,nr);
      
      //prompt found conf and r
      VERBOSITY_LV2_MASTER_PRINTF("Conf: %d, r: %d\n",iconf,r);
      
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
void clusterize()
{
  
  for(int r=0;r<nr;r++)
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      clusterize(out_buffer+0+(njacks+1)*(ivol+loc_vol*r));
  NISSA_PARALLEL_LOOP_END;
  THREAD_BARRIER();
}

//save the clusterized data
void save_data(const char *path)
{
  FILE *fout=open_file(path,"w");
  
  if(!little_endian) change_endianness(out_buffer,out_buffer,loc_vol*nr*(njacks+1));
  int rc=fwrite(out_buffer,sizeof(double),loc_vol*nr*(njacks+1),fout);
  if(rc!=loc_vol*nr*(njacks+1)) CRASH("returned %d instead of %d",rc,loc_vol*nr*(njacks+1));
     
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
  if(narg<2) CRASH("Use %s input",arg[0]);
  init_clusterize(arg[1]);
  
  load_data(path_in);
  clusterize();
  save_data(path_out);
  
  tot_prog_time+=take_time();
  close_clusterize();
}

int main(int narg,char **arg)
{
  initNissa_threaded(narg,arg,in_main);
  closeNissa();
    
  return 0;
}
