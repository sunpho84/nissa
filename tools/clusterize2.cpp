#include "nissa.h"

typedef char sa_string[100];

int T,ncorr,ncombo,ncorr_type;
int njack,clust_size;
int start_file_id,nfile;
char base_path_in[1024],base_path_out[1024];
sa_string *corr_name,*out_path;
double *data;

void parse_input(char *path)
{
  //Open input
  open_input(path);

  //in pars
  read_str_int("T",&T);
  read_str_str("BasePathIn",base_path_in,1024);
  read_str_int("StartFileId",&start_file_id);

  //nfiles
  int nfile_teo;
  read_str_int("NFile",&nfile_teo);
  read_str_int("NJack",&njack);
  
  //n corrs
  read_str_str("BasePathOut",base_path_out,1024);
  
  close_input();
  
  //compute real numb of files
  clust_size=nfile_teo/njack;
  nfile=clust_size*njack;

  //allocate
  corr_name=nissa_malloc("corr_name",ncorr,sa_string);
  out_path=nissa_malloc("out_path",ncorr,sa_string);
  
  //allocate room for data
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
  /*
      if(strcmp(line,line_to_find[icorr])==0)
	ncombo[icorr]++;
  */
  fclose(fin);
  
  ncorr_type=ncorr/ncombo;
  
  /*
  for(int icorr=0;icorr<ncorr;icorr++)
    master_printf("%s %d\n",line_to_find[icorr],ncombo[icorr]);
  */
  
  printf("ncombo: %d\n",ncombo);
  printf("ncorr: %d\n",ncorr);
  printf("ncorr_type: %d\n",ncorr_type);
  
}

void close_clust()
{
  nissa_free(corr_name);
  nissa_free(out_path);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_nissa();
  
  if(narg<2) crash("uese %s input",arg[0]);
  
  //read input and prepare structures
  parse_input(arg[1]);
    
  //create first file path
  char path[1024];
  sprintf(path,base_path_in,start_file_id);
  count_corr(path);
  
  /*
  //loop over all the confs
  int targ_conf=0;
  for(int iconf=0;iconf<nconf;iconf++)
    {
      FILE *fin;
      
      printf("%d\n",iconf);
        
      //open file
      do
        {
          sprintf(path,arg[1],targ_conf+base_gauge);
          fin=fopen(path,"r");
          if(fin==NULL)
            {
              fprintf(stderr,"Error, couldn't open file: %s, skipping\n",path);
              targ_conf++;
              if(targ_conf>=nconf_teo)
                {
                  fprintf(stderr,"Error, finished the available confs\n");
                  exit(1);
                }
            }
          
        }
      while(fin==NULL);
      targ_conf++;
      
      dwqd
//search the 
      for(int icorr=0;icorr<ncorr;icorr++)
        {
          char line[1024],*out;
          do out=fgets(line,1024,fin);
          while(out==line && (strcmp(line,line_to_find)!=0));
          
          //read
          for(int t=0;t<T;t++)
            {
              int r=fscanf(fin,"%lg %lg",&(data[icorr][0][t][iconf]),&(data[icorr][1][t][iconf]));
              if(r!=2)
                {
                  fprintf(stderr,"Error reading from file %d time: %d\n",iconf,t);
                  exit(1);
                }
            }
        }
      
      //close file
      fclose(fin);
    }
  
  FILE *fout=fopen(arg[5],"w");
  if(fout==NULL)
    {
      fprintf(stderr,"Error, couldn't open for output file: %s\n",arg[5]);
      exit(1);
    }
  */
  //exiting
  close_clust();
  
  return 0;
}
