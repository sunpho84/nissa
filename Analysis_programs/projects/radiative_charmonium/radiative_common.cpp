#include "include.h"
#include "../nf2/common_pars.cpp"

void read_ensemble_pars(char *base_path,int &T,int &ibeta,int &nmass,double *&mass,double &sea_mass,const char *data_list_file)
{
  FILE *input=open_file(data_list_file,"r");
  
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","njack"); //dum
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");

  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");
  read_formatted_from_file_expecting((char*)&sea_mass,input,"%lg","sea_mass");
  
  fclose(input);
}

void boot_from_jack(boot &out,jack in,int *iboot_jack)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
  out.data[nboot]=in.data[njack];
}
