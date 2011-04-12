#pragma once
#include <stdio.h>
#include <stdlib.h>

FILE *open_file(const char* path,const char* mod)
{
  FILE *out=fopen(path,mod);

  if(out==NULL)
    {
      fprintf(stderr,"Error opening file '%s' in mode '%s'\n",path,mod);
      exit(1);
    }
  
  return out;
}
