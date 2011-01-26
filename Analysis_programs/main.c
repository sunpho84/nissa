#include "input.c"
#include <string.h>

void fuf()
{
  char line[1024];
  char *out;

  do out=fgets(line,1024,input_global);
  while(out!=NULL && (line[0]=='#' || strlen(line)<=1));
  
  double dou,dod;
  if(out!=NULL)
    {
      dou=strtod(line,&out);
      dod=strtod(out,NULL);
      printf("%g %g\n",dou,dod);
    }
}

int main(int narg,char **arg)
{
  int ok;

  if(narg<2)
    {
      fprintf(stderr,"Error, use %s input_file\n",arg[0]);
      exit(1);
    }
  
  open_input(arg[1]);
  
  int nconf;
  read_str_int("NConf",&nconf);
  fuf();
  fuf();
  
  return 0;
}
