#include "nissa.h"

typedef char sa_string[100];

int T,ncorr;
int njack,clust_size;
int start_file_id,nfile;
char base_path[1024];
sa_string *corr_name,*out_path;

void parse_input(char *path)
{
  //Open input
  open_input(path);

  //in pars
  read_str_int("T",&T);
  read_str_str("BasePath",base_path,1024);
  read_str_int("StartFileId",&start_file_id);

  //nfiles
  int nfile_teo;
  read_str_int("NFile",&nfile_teo);
  read_str_int("NJack",&njack);
  clust_size=nfile_teo/njack;
  nfile=clust_size*njack;
  
  //n corrs
  read_str_int("NCorr",&ncorr);
  corr_name=nissa_malloc("corr_name",ncorr,sa_string);
  out_path=nissa_malloc("out_path",ncorr,sa_string);
  
  //load corr names and out
  for(int icorr=0;icorr<ncorr;icorr++)
    {
      read_str(corr_name[icorr],100);
      read_str(out_path[icorr],100);
    }
  
  close_input();
}

//count the number of combo
void count_corr(int *ncombo,char *path,const char **line_to_find)
{
  memset(ncombo,0,sizeof(int)*ncorr);
  
  FILE *fin=open_file(path,"r");
  
  char line[1024];
  while(fgets(line,1024,fin)==line)
    for(int icorr=0;icorr<ncorr;icorr++)
      if(strcmp(line,line_to_find[icorr])==0)
	ncombo[icorr]++;
  
  fclose(fin);
  
  for(int icorr=0;icorr<ncorr;icorr++)
    printf("%s %d\n",line_to_find[icorr],ncombo[icorr]);
}


int main(int narg,char **arg)
{
  init_nissa();
  
  if(narg<2) crash("uese %s input",arg[0]);
  
  parse_input(arg[1]);
  
  close_nissa();
  
  return 0;
}
