#include "nissa.h"

//use at most this memory
int max_mem_usable=128000000;

typedef char sa_string[300];

int T,ncorr,ncombo,ncorr_type;
int ncombo_per_block;
int njack,clust_size;
int start_file_id,nfile,nfile_teo;
char base_path_in[1024],base_path_out[1024];
sa_string *corr_name,*outpath;
double *data;
int mem_asked,nblock;
int *pos;

void parse_input(char *path)
{
  //Open input
  open_input(path);

  //in pars
  read_str_int("T",&T);
  read_str_str("BasePathIn",base_path_in,1024);
  read_str_int("StartFileId",&start_file_id);

  //nfiles
  read_str_int("NFile",&nfile_teo);
  read_str_int("NJack",&njack);
  
  //n corrs
  read_str_str("BasePathOut",base_path_out,1024);
  
  close_input();
  
  //compute real numb of files
  clust_size=nfile_teo/njack;
  nfile=clust_size*njack;
}

//count the number of combo
void count_corr(char *path)
{
  ncombo=ncorr=0;
  
  master_printf("Considering file: %s\n",path);
  
  FILE *fin=open_file(path,"r");
  
  char line[1024];
  while(fgets(line,1024,fin)==line)
    if(line[1]=='#')
      {
	char test1[100],test2[100];
	int ntest1=sscanf(line," # %s ",test1);
	int ntest2=sscanf(line+strlen(test1)+3,"%s",test2);
	
	if(ntest1>0)
	  if(ntest2<=0) ncorr++;
	  else ncombo++;
      }
  
  ncorr_type=ncorr/ncombo;
  
  printf("ncombo: %d\n",ncombo);
  printf("ncorr: %d\n",ncorr);
  printf("ncorr_type: %d\n",ncorr_type);
  
  //alloc space for the name of corrs and outpath
  corr_name=nissa_malloc("corr_name",ncorr_type,sa_string);
  outpath=nissa_malloc("outpath",ncorr,sa_string);
  pos=nissa_malloc("pos",nfile,int);
  memset(pos,0,sizeof(int)*nfile);
  
  //go to the file start and load all corr type names
  fseek(fin,0,SEEK_SET);  
  uint64_t icorr_type=0;
  while(fgets(line,1024,fin)==line && icorr_type<ncorr_type)
    if(line[1]=='#')
      {
	char test[100];
	int ntest1=sscanf(line," # %s ",corr_name[icorr_type]);
	int ntest2=sscanf(line+strlen(corr_name[icorr_type])+3,"%s",test);
	
	if(ntest1>0 && ntest2<=0)
	  {
	    printf("Corr %llu: %s\n",icorr_type,corr_name[icorr_type]);
	    icorr_type++;
	  }
      }
    
  fclose(fin);
  
  //compute size of each combo
  uint64_t combo_size=ncorr_type*2*T*(njack+1)*sizeof(double);
  
  //compute the total amount of memory needed
  uint64_t max_mem_needed=ncombo*combo_size;
  
  printf("Max memory needed: %llu\n",max_mem_needed);
  printf("Max memory usable: %llu\n",max_mem_usable);
  
  //count number of blocks needed
  nblock=max_mem_needed/max_mem_usable;
  
  //ceil
  if(nblock*max_mem_usable<max_mem_needed) nblock++;
  
  //if ncombo not multiple of ncombo_per_block increase nblock
  while(nblock*ncombo_per_block!=ncombo)
    {
      nblock++;
      ncombo_per_block=ncombo/nblock;
    }
  if(nblock>ncombo) crash("boh");
  
  //compute the number of combo in each block and the memory asked
  ncombo_per_block=ncombo/nblock;
  mem_asked=ncombo_per_block*combo_size;
  
  printf("nblock: %d\n",nblock);
  printf("ncombo per block: %d\n",ncombo_per_block);
  
  //allocate room for data
  data=nissa_malloc("data",mem_asked/sizeof(double),double);
  
  //prepare the output names
  for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
    sprintf(outpath[icorr_type],base_path_out,corr_name[icorr_type]);
  
  //check that paths differs
  for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
    for(int jcorr_type=icorr_type+1;jcorr_type<ncorr_type;jcorr_type++)
      if(strcmp(outpath[icorr_type],outpath[jcorr_type])==0)
	crash("%d corr outpath %s is equal to %d outpath, %s",icorr_type,outpath[icorr_type],jcorr_type,outpath[jcorr_type]);
}

void close_clust()
{
  nissa_free(corr_name);
  nissa_free(outpath);
  nissa_free(data);
  nissa_free(pos);
  
  close_nissa();
}

void parse_file(int ifile,char *path,int &start)
{
  int iclust=ifile/clust_size;
  printf("clust: %d/%d\n",iclust,njack);
  
  //open and seek
  FILE *file=open_file(path,"r");
  if(fseek(file,start,SEEK_SET)) crash("while seeking on %s",path);
  
  char line[1024];
  
  int nread_corr=0,nread_line=0;
  
  for(int icombo=0;icombo<ncombo_per_block;icombo++)
    for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
      {
	//search the corr
	do
	  {
	    char *succ=fgets(line,1024,file);
	    nread_line++;
	    if(succ!=line) crash("file ended before finishing reading");
	  }
	while(line[1]=='#' || strlen(line)<=1);
	
	//advance the corr read
	nread_corr++;
	
	//read the corr
	for(int t=0;t<T;t++)
	  {
	    //if not t==0, read line
	    if(t>0)
	      do if(line!=fgets(line,1024,file)) crash("error reading line, obtained: %s",line);
	      while(strlen(line)<=1);
	    
	    //scanning line
	    double t1,t2;
	    int n=sscanf(line,"%lg %lg",&t1,&t2);
	    if(n!=2) crash("scanning line %s obtained only %d numbers",line,n);
	    nread_line++;
	    
	    //find output place
	    int i1=iclust+(njack+1)*(t+T*(0+2*(icombo+ncombo_per_block*icorr_type)));
	    int i2=iclust+(njack+1)*(t+T*(1+2*(icombo+ncombo_per_block*icorr_type)));
	    
	    //summ into the cluster
	    data[i1]+=t1;
	    data[i2]+=t2;
	  }
      }  
  
  //check to have finished the file
  //while(line==fgets(line,1024,file))
  //if(strlen(line)>1) crash("should have reached end and instead got: %s",line);
  
  //get pos
  start=ftell(file);
  
  fclose(file);
}

int main(int narg,char **arg)
{
  init_nissa();
  
  if(narg<2) crash("use %s input",arg[0]);
  
  //read input and prepare structures
  parse_input(arg[1]);
    
  //create first file path
  char path[1024];
  sprintf(path,base_path_in,start_file_id);
  count_corr(path);
  
  //looping over blocks
  for(int iblock=0;iblock<nblock;iblock++)
    {
      //reset data
      memset(data,0,mem_asked);
      
      //loop over all the confs
      int targ_file=0;
      for(int ifile=0;ifile<nfile;ifile++)
	{
	  printf("Considering file %d\n",ifile);
	  
	  //search existing file
	  int found;
	  do
	    {
	      sprintf(path,base_path_in,targ_file+start_file_id);
	      found=file_exists(path);
	      if(!found)
		{
		  master_printf("file %s not available, skipping",path);
		  targ_file++;
		  if(targ_file>=nfile_teo) crash("finished all available files");
		}
	    }
	  while(!found);
	  
	  //parse the file into data clusters
	  parse_file(ifile,path,pos[targ_file]);
	  
	  //advance
	  targ_file++;
	}
      
      //make the jacknives
      for(int icorr=0;icorr<ncombo_per_block*2*ncorr_type*T;icorr++)
	{
	  double *d=data+icorr*(njack+1);
	  
	  //compute average
	  for(int iclust=0;iclust<njack;iclust++) d[njack]+=d[iclust];
	  for(int ijack=0;ijack<njack;ijack++) d[ijack]=(d[njack]-d[ijack])/(nfile-nfile/njack);
	  d[njack]/=nfile;
	}
      
      //if necessary change the endianess
      if(!little_endian)
	doubles_to_doubles_changing_endianess(data,data,mem_asked/sizeof(double));
      
      //write different files for each corr
      for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
	{
	  FILE *file;
	  if(iblock==0) file=open_file(outpath[icorr_type],"w");
	  else file=open_file(outpath[icorr_type],"a");
	  
	  int n=fwrite((void*)(data+ncombo_per_block*2*T*(njack+1)*icorr_type),sizeof(double),ncombo_per_block*2*T*(njack+1),file);
	  fclose(file);
	}
    }
  
  //exiting
  close_clust();
  
  return 0;
}
