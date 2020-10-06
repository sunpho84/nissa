#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>

#include "common_tools.cpp"

uint64_t mem_needed;

typedef char sa_string[300];

int REIM,use_new_contraction_layout,binary_out;
int T,ncorr,ncombo,ncorr_type,ntriple;
int nfile_names;
int njack,clust_size;
int start_conf_id,nconfs,nconfs_teo,confs_each;
int start_copy_id,ncopies;
char base_path_in[1024],base_path_out[1024];
sa_string *corr_name,*outpath;
double *data;

void parse_input(char *path)
{
  //Open input
  open_input(path);
  
  //in pars
  read_str_int("T",&T);
  read_str_int("StartConfId",&start_conf_id);
  
  //nconfs
  read_str_int("Nconfs",&nconfs_teo);
  read_str_int("ConfsEach",&confs_each);
  read_str_int("StartCopyId",&start_copy_id);
  read_str_int("Ncopies",&ncopies);
  read_str_int("NJack",&njack);
  
  //n files
  read_str_int("NFileNames",&nfile_names);
  
  //REIM
  read_str_int("UseNewContractionLayout",&use_new_contraction_layout);
  REIM=use_new_contraction_layout?1:2;
  
  //binary or ascii
  read_str_int("BinaryOut",&binary_out);
  
  //compute real numb of confs
  clust_size=nconfs_teo/njack;
  nconfs=clust_size*njack;
}

//count the number of combo
void count_corr(char *path)
{
  ncombo=ncorr=ntriple=0;
  
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
    else
      {
	int t;
	double dre,dim;
	int ntest3=sscanf(line,"%d %lg %lg",&t,&dre,&dim);
	if(ntest3>0) ntriple++;
      }
  
  //if no header was present
  if(ncombo==0) ncombo=1;
  
  if(ncorr==0) ncorr=ncombo=ntriple/T;
  
  ncorr_type=ncorr/ncombo;
  
  printf("ncombo: %d\n",ncombo);
  printf("ncorr: %d\n",ncorr);
  printf("ncorr_type: %d\n",ncorr_type);
  
  //alloc space for the name of corrs and outpath
  corr_name=(sa_string*)malloc(ncorr_type*sizeof(sa_string));
  outpath=(sa_string*)malloc(ncorr*sizeof(sa_string));
  
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
  int combo_size=ncorr_type*REIM*T*(njack+1)*sizeof(double);
  
  //compute the total amount of memory needed
  mem_needed=(uint64_t)ncombo*combo_size;
  
  printf("memory needed: %lu\n",mem_needed);
  printf("\n");
  
  //allocate room for data
  data=(double*)malloc(mem_needed);
  
  //prepare the output names
  for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
    sprintf(outpath[icorr_type],base_path_out,corr_name[icorr_type]);
  
  //check that paths differs
  for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
    for(int jcorr_type=icorr_type+1;jcorr_type<ncorr_type;jcorr_type++)
      if(strcmp(outpath[icorr_type],outpath[jcorr_type])==0)
	crash("%d corr outpath %s is equal to %d outpath, %s",icorr_type,outpath[icorr_type],jcorr_type,outpath[jcorr_type]);
}

void parse_conf(int iconf,char *path)
{
  int iclust=iconf/clust_size;
  printf("Considering conf %d/%d (%s) on thread %d\n",iconf+1,nconfs,path,omp_get_thread_num());

  //open
  FILE *file=open_file(path,"r");
  
  char line[1024];
  
  int nread_corr=0,nread_line=0;
  
  for(int icombo=0;icombo<ncombo;icombo++)
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
	    int n;
	    if(use_new_contraction_layout)
	      {
		n=sscanf(line,"%lg",&t1);
		if(n!=1) crash("scanning line '%s' obtained only %d numbers",line,n);
	      }
	    else
	      {
		n=sscanf(line,"%lg %lg",&t1,&t2);
		if(n!=2) crash("scanning line '%s' obtained only %d numbers",line,n);
	      }
	    nread_line++;
	    
	    //find output place
	    int i1=iclust+(njack+1)*(t+T*(0+REIM*(icombo+ncombo*icorr_type)));
	    int i2=iclust+(njack+1)*(t+T*(1+REIM*(icombo+ncombo*icorr_type)));
	    
	    //summ into the cluster
	    data[i1]+=t1;
	    if(REIM) data[i2]+=t2;
	  }

	//advance the corr read
	nread_corr++;
      }  
  
  //check to have finished the file
  while(line==fgets(line,1024,file))
  if(strlen(line)>1) crash("should have reached end and instead got: %s",line);
  
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
      char path[nconfs][ncopies][100];
      for(int icopy=0;icopy<ncopies;icopy++)
	{
	  int curr_conf=0;
	  for(int iconf=0;iconf<nconfs;iconf++)
	    {
	      //search existing conf
	      int found;
	      do
		{
		  sprintf(path[iconf][icopy],base_path_in,curr_conf+start_conf_id,icopy+start_copy_id);
		  found=file_exists(path[iconf][icopy]);
		  if(!found)
		    {
		      printf("conf %s not available, skipping",path[iconf][icopy]);
		      curr_conf+=confs_each;
		      if(curr_conf>=nconfs_teo*confs_each) crash("finished all available confs");
		    }
		}
	      while(!found);
	      
	      //increment
	      curr_conf+=confs_each;
	    }
	}
      
      //count corrs and reset data
      count_corr(path[0][0]);
      memset(data,0,mem_needed);
      
      //parse the conf into data clusters
#pragma omp parallel for
      for(int ijack=0;ijack<njack;ijack++)
	for(int iconf=ijack*clust_size;iconf<(ijack+1)*clust_size;iconf++)
	  for(int icopy=0;icopy<ncopies;icopy++)
	    parse_conf(iconf,path[iconf][icopy]);
      
      //make the jacknives
      for(int icorr=0;icorr<ncombo*REIM*ncorr_type*T;icorr++)
	{
	  double *d=data+icorr*(njack+1);
	  
	  //compute average
	  for(int iclust=0;iclust<njack;iclust++) d[njack]+=d[iclust];
	  for(int ijack=0;ijack<njack;ijack++) d[ijack]=(d[njack]-d[ijack])/(nconfs-nconfs/njack)/ncopies;
	  d[njack]/=nconfs*ncopies;
	}
      
      //if necessary change the endianess
      if(binary_out && !little_endian)
	doubles_to_doubles_changing_endianess(data,data,mem_needed/sizeof(double));
      
      printf("\n");
      
      //write different files for each corr
      for(int icorr_type=0;icorr_type<ncorr_type;icorr_type++)
	{
	  FILE *file;
	  printf("Writing corr %s (%d/%d)\n",corr_name[icorr_type],icorr_type+1,ncorr_type);
	  file=open_file(outpath[icorr_type],"w");
	  
	  if(binary_out)
	    {
	      int n=fwrite((void*)(data+ncombo*REIM*T*(njack+1)*icorr_type),sizeof(double),ncombo*REIM*T*(njack+1),file);
	      if(n!=ncombo*REIM*T*(njack+1)) crash("obtained %d instead of %d",n,ncombo*REIM*T*(njack+1));
	    }
	  else
	    for(int i=0;i<ncombo*REIM*T*(njack+1);i++) fprintf(file,"%+016.16lg\n",
							       data[i+ncombo*REIM*T*(njack+1)*icorr_type]);
	  fclose(file);
	}
      
      printf("\n");
    }
  
  //free
  free((void*)data);
  free((void*)corr_name);
  free((void*)outpath);
  
  return 0;
}
