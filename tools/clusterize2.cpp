#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>

#include "common_tools.cpp"

//use at most this memory
uint64_t max_mem_usable=128000000;

typedef char sa_string[300];

int T,ncorr,ncombo,ncorr_type;
int nfile_names,ncombo_per_block;
int njack,clust_size,nblock;;
int start_conf_id,nconfs,nconfs_teo;
char base_path_in[1024],base_path_out[1024];
sa_string *corr_name,*outpath;
double *data;
uint64_t mem_asked;
int *pos;

void parse_input(char *path)
{
  //Open input
  open_input(path);
  
  //in pars
  read_str_int("T",&T);
  read_str_int("StartConfId",&start_conf_id);
  
  //nconfs
  read_str_int("Nconfs",&nconfs_teo);
  read_str_int("NJack",&njack);
  
  //n files
  read_str_int("NFileNames",&nfile_names);
  
  //compute real numb of confs
  clust_size=nconfs_teo/njack;
  nconfs=clust_size*njack;
}

//count the number of combo
void count_corr(char *path)
{
  ncombo=ncorr=0;
  
  printf("Considering conf: %s\n",path);
  
  FILE *fin=open_file(path,"r");
  
  char line[1024];
  while(fgets(line,1024,fin)==line)
    if(line[1]=='#')
      {
	char test1[100],test2[100];
	int ntest1=sscanf(line," # %s ",test1);
	int ntest2=sscanf(line+strlen(test1)+3,"%s",test2);
	
	if(ntest1>0)
	  {
	    if(ntest2<=0) ncorr++;
	    else ncombo++;
	  }
      }
  
  ncorr_type=ncorr/ncombo;
  
  printf("ncombo: %d\n",ncombo);
  printf("ncorr: %d\n",ncorr);
  printf("ncorr_type: %d\n",ncorr_type);
  
  //alloc space for the name of corrs and outpath
  corr_name=(sa_string*)malloc(ncorr_type*sizeof(sa_string));
  outpath=(sa_string*)malloc(ncorr*sizeof(sa_string));
  pos=(int*)malloc(nconfs_teo*sizeof(int));
  memset(pos,0,sizeof(int)*nconfs_teo);
  
  //go to the file start and load all corr type names
  fseek(fin,0,SEEK_SET);  
  int icorr_type=0;
  while(fgets(line,1024,fin)==line && icorr_type<ncorr_type)
    if(line[1]=='#')
      {
	char test[100];
	int ntest1=sscanf(line," # %s ",corr_name[icorr_type]);
	int ntest2=sscanf(line+strlen(corr_name[icorr_type])+3,"%s",test);
	
	if(ntest1>0 && ntest2<=0)
	  {
	    printf("Corr %d: %s\n",icorr_type,corr_name[icorr_type]);
	    icorr_type++;
	  }
      }
    
  fclose(fin);
  
  //compute size of each combo
  int combo_size=ncorr_type*2*T*(njack+1)*sizeof(double);
  
  //compute the total amount of memory needed
  uint64_t max_mem_needed=(uint64_t)ncombo*combo_size;
  
  printf("max memory needed: %ju\n",max_mem_needed);
  printf("max memory usable: %ju\n",max_mem_usable);
  
  //count number of blocks needed
  if(max_mem_needed<max_mem_usable) nblock=1;
  else
    {
      nblock=max_mem_needed/max_mem_usable;
      if(nblock*max_mem_usable<max_mem_needed) nblock++;
    }
  
  //if ncombo not multiple of ncombo_per_block increase nblock
  ncombo_per_block=ncombo/nblock;
  while(nblock*ncombo_per_block!=ncombo)
    {
      nblock++;
      ncombo_per_block=ncombo/nblock;
      if(ncombo_per_block<1) crash("not enough memory for a single combo!");
    }
  if(nblock>ncombo) crash("something went wrong when computing nblock");
  
  //compute the number of combo in each block and the memory asked
  mem_asked=ncombo_per_block*combo_size;
  
  printf("\n");
  printf(" memory asked: %ju\n",mem_asked);
  printf(" nblock: %d\n",nblock);
  printf(" ncombo per block: %d\n",ncombo_per_block);
  printf("\n");
  
  //allocate room for data
  data=(double*)malloc(mem_asked);
  
  //prepare the output names
  for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
    sprintf(outpath[icorr_type],base_path_out,corr_name[icorr_type]);
  
  //check that paths differs
  for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
    for(int jcorr_type=icorr_type+1;jcorr_type<ncorr_type;jcorr_type++)
      if(strcmp(outpath[icorr_type],outpath[jcorr_type])==0)
	crash("%d corr outpath %s is equal to %d outpath, %s",icorr_type,outpath[icorr_type],jcorr_type,outpath[jcorr_type]);
}

void parse_conf(int iconf,char *path,int &start)
{
  int iclust=iconf/clust_size;
  
  //open and seek
  FILE *file=open_file(path,"r");
  if(fseek(file,start,SEEK_SET)) crash("while seeking on %s",path);
  
  char line[1024];
  
  int nread_corr=0,nread_line=0;
  
  for(int icombo=0;icombo<ncombo_per_block;icombo++)
    for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
      {
	//search the corr
	int nonblank;
	do
	  {
	    char *succ=fgets(line,1024,file);
	    nread_line++;
	    if(succ!=line) crash("file ended before finishing reading");

	    nonblank=0;
	    for(size_t i=0;i<strlen(line)-1;i++) nonblank+=(line[i]!=' ');
	  }
	while(line[1]=='#' || strlen(line)<= 1||!nonblank);
	
	//advance the corr read
	nread_corr++;
	
	//read the corr
	for(int t=0;t<T;t++)
	  {
	    //if not t==0, read line
	    if(t>0)
	      {
		int nonblank;
		do
		  {
		    if(line!=fgets(line,1024,file)) crash("error reading line, obtained: %s",line);
		    nonblank=0;
		    for(size_t i=0;i<strlen(line)-1;i++) nonblank+=(line[i]!=' ');
		  }
		while(strlen(line)<=1||!nonblank);
	      }
	    
	    //scanning line
	    double t1,t2;
	    int n=sscanf(line,"%lg %lg",&t1,&t2);
	    if(n!=2) crash("scanning line '%s' obtained only %d numbers",line,n);
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
  if(narg<2) crash("use %s input",arg[0]);
  
  check_endianess();
  
  //read input and prepare structures
  parse_input(arg[1]);
  
  for(int ifile_name=0;ifile_name<nfile_names;ifile_name++)
    {
      read_str_str("BasePathIn",base_path_in,1024);
      read_str_str("BasePathOut",base_path_out,1024);
      
      //create first file path
      char path[1024];
      sprintf(path,base_path_in,start_conf_id);
      count_corr(path);
      
      //looping over blocks
      for(int iblock=0;iblock<nblock;iblock++)
	{
	  if(iblock>1) printf("Running on block %d/%d\n",iblock+1,nblock);
	  
	  //reset data
	  memset(data,0,mem_asked);
	  
	  //loop over all the confs
	  int targ_conf=0;
	  for(int iconf=0;iconf<nconfs;iconf++)
	    {
	      printf("Considering conf %d/%d\n",iconf+1,nconfs);
	      
	      //search existing conf
	      int found;
	      do
		{
		  sprintf(path,base_path_in,targ_conf+start_conf_id);
		  found=file_exists(path);
		  if(!found)
		    {
		      printf("conf %s not available, skipping",path);
		      targ_conf++;
		      if(targ_conf>=nconfs_teo) crash("finished all available confs");
		    }
		}
	      while(!found);
	      
	      //parse the conf into data clusters
	      parse_conf(iconf,path,pos[targ_conf]);
	      
	      //advance
	      targ_conf++;
	    }
	  
	  //make the jacknives
	  for(int icorr=0;icorr<ncombo_per_block*2*ncorr_type*T;icorr++)
	    {
	      double *d=data+icorr*(njack+1);
	      
	      //compute average
	      for(int iclust=0;iclust<njack;iclust++) d[njack]+=d[iclust];
	      for(int ijack=0;ijack<njack;ijack++) d[ijack]=(d[njack]-d[ijack])/(nconfs-nconfs/njack);
	      d[njack]/=nconfs;
	    }
	  
	  //if necessary change the endianess
	  if(!little_endian)
	    doubles_to_doubles_changing_endianess(data,data,mem_asked/sizeof(double));
	  
	  printf("\n");
	  
	  //write different files for each corr
	  for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
	    {
	      FILE *file;
	      if(iblock==0)
		{
		  printf("Writing corr %s (%d/%d)\n",corr_name[icorr_type],icorr_type+1,ncorr_type);
		  file=open_file(outpath[icorr_type],"w");
		}
	      else
		{
		  printf("Appending corr %s (%d/%d)\n",corr_name[icorr_type],icorr_type+1,ncorr_type);
		  file=open_file(outpath[icorr_type],"a");
		}
	      
	      int n=fwrite((void*)(data+ncombo_per_block*2*T*(njack+1)*icorr_type),sizeof(double),ncombo_per_block*2*T*(njack+1),file);
	      if(n!=ncombo_per_block*2*T*(njack+1)) crash("obtanied %d instead of %d",n,ncombo_per_block*2*T*(njack+1));
	      fclose(file);
	    }
	  
	  printf("\n");
	}
      
      //free
      free((void*)data);
      free((void*)corr_name);
      free((void*)outpath);
      free((void*)pos);
    }
  
  return 0;
}
