#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>

//use at most this memory
uint64_t max_mem_usable=128000000;

typedef char sa_string[300];

FILE *input_global;

int T,ncorr,ncombo,ncorr_type;
int nfile_names,ncombo_per_block;
int njack,clust_size,nblock;;
int start_conf_id,nconfs,nconfs_teo;
char base_path_in[1024],base_path_out[1024];
sa_string *corr_name,*outpath;
double *data;
uint64_t mem_asked;
int *pos;
int little_endian;

//check the endianess of the machine
void check_endianess()
{
  little_endian=1;
  little_endian=(int)(*(char*)(&little_endian));
}

//fucking tool to revert the endianess of doubles
void doubles_to_doubles_changing_endianess(double *dest,double *sour,int ndoubles)
{
  char *cdest,*csour;
  char temp;
  
  if(dest==sour)
    for(int idouble=0;idouble<ndoubles;idouble++)
      {
        cdest=(char*)(dest+idouble);
        csour=(char*)(sour+idouble);
        
        temp=csour[7];
        csour[7]=cdest[0];
        cdest[0]=temp;

        temp=csour[6];
        csour[6]=cdest[1];
        cdest[1]=temp;

        temp=csour[5];
        csour[5]=cdest[2];
        cdest[2]=temp;

        temp=csour[4];
        csour[4]=cdest[3];
        cdest[3]=temp;
      }
  else
    for(int idouble=0;idouble<ndoubles;idouble++)
      {
        cdest=(char*)(dest+idouble);
        csour=(char*)(sour+idouble);
        
        cdest[0]=csour[7];
        cdest[1]=csour[6];
        cdest[2]=csour[5];
        cdest[3]=csour[4];
        cdest[4]=csour[3];
        cdest[5]=csour[2];
        cdest[6]=csour[1];
        cdest[7]=csour[0];

      }
}

void crash(const char *templ,...)
{
  va_list ap;
  va_start(ap,templ);
  
  char mess[1024];
  vsprintf(mess,templ,ap);
  va_end(ap);
  
  fprintf(stderr,"%s\n",mess);  
  
  exit(1);
}

//read a token from file
int read_next_token(char *tok)
{
  int ok=fscanf(input_global,"%s",tok);

  return ok;
}

//check whether a token starts or end a comment
int check_tok_starts_comment(char *tok)
{return strncasecmp(tok,"/*",2)==0;}
int check_tok_ends_comment(char *tok)
{return strncasecmp(tok+strlen(tok)-2,"*/",2)==0;}

//read up to "*/"
void read_up_to_end_of_comment()
{
  char tok[1024];
  
  do
    {
      int ok=read_next_token(tok);
      if(ok!=1) crash("reached end of file without finding end of comment");
    }
  while(check_tok_ends_comment(tok)==0);
}

int read_var_catcherr(char *out,const char *par,int size_of)
{
  char tok[1024];
  int ok,comment_found;
  
  do
    {
      //read the token
      ok=read_next_token(tok);
      
      //check if token comment
      if(check_tok_starts_comment(tok))
        {
          comment_found=1;
          //if the comments do not ends itself, read up to finding its end
          if(!check_tok_ends_comment(tok)) read_up_to_end_of_comment();
        }
      else comment_found=0;
    }
  while(comment_found);
  
  //parse the token
  if(ok==1)
    ok=(sscanf(tok,par,out)>0);
  
  return ok;
}

void read_var(char *out,const char *par,int size_of)
{if(!read_var_catcherr(out,par,size_of)) crash("Couldn't read from input file!!!");}

//Read an integer from the file
void read_int(int *out)
{read_var((char*)out,"%d",sizeof(int));}

//Read a double from the file
void read_double(double *out)
{read_var((char*)out,"%lg",sizeof(double));}

//Read a string from the file
void read_str(char *str,int length)
{read_var(str,"%s",length);}

//Read a string from the file and check against the argument
void expect_str(const char *exp_str)
{
  char obt_str[1024];
  
  read_str(obt_str,1024);
  
  if(strcasecmp(exp_str,obt_str)!=0) crash("Error, expexcted '%s' in input file, obtained: '%s'",exp_str,obt_str);
}

//Read an int checking the tag
void read_str_int(const char *exp_str,int *in)
{
  expect_str(exp_str);
  read_int(in);
}

//Read a double checking the tag
void read_str_double(const char *exp_str,double *in)
{
  expect_str(exp_str);
  read_double(in);
}

//Read a string checking the tag
void read_str_str(const char *exp_str,char *in,int length)
{
  expect_str(exp_str);
  read_str(in,length);
}

FILE* open_file(const char *outfile,const char *mode)
{
  FILE *fout=NULL;
  
  fout=fopen(outfile,mode);
  if(fout==NULL) crash("Couldn't open the file: %s for mode: %s",outfile,mode);
  
  return fout;
}

//check if a file exists
int file_exists(const char *path)
{
  int status=1;
  
  FILE *f=fopen(path,"r");
  if(f!=NULL)
    {
      status=1;
      fclose(f);
    }
  else status=0;

  return status;
}

void open_input(char *input_path)
{
  input_global=fopen(input_path,"r");
  if(input_global==NULL) crash("File '%s' not found",input_path);
}

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
	  if(ntest2<=0) ncorr++;
	  else ncombo++;
      }
  
  ncorr_type=ncorr/ncombo;
  
  printf("ncombo: %d\n",ncombo);
  printf("ncorr: %d\n",ncorr);
  printf("ncorr_type: %d\n",ncorr_type);
  
  //alloc space for the name of corrs and outpath
  corr_name=(sa_string*)malloc(ncorr_type*sizeof(sa_string));
  outpath=(sa_string*)malloc(ncorr*sizeof(sa_string));
  pos=(int*)malloc(nconfs*sizeof(int));
  memset(pos,0,sizeof(int)*nconfs);
  
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
  
  printf("max memory needed: %llu\n",max_mem_needed);
  printf("max memory usable: %llu\n",max_mem_usable);
  
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
  printf(" memory asked: %llu\n",mem_asked);
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
