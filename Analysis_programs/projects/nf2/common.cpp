#pragma once

#include "include.h"
#include "common_pars.cpp"

void read_ensemble_pars(char *base_path,int &T,int &ibeta,int &nmass,double *&mass,int &iml_un,int &nlights,const char *data_list_file)
{
  FILE *input=open_file(data_list_file,"r");
  
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");

  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");
  read_formatted_from_file_expecting((char*)&iml_un,input,"%d","iml_un");
  read_formatted_from_file_expecting((char*)&nlights,input,"%d","nlights");
  
  fclose(input);
}

int icombo(int im1,int im2,int nmass,int nlights,int mode)
{
  int imc,ims,ic;
  
  switch(mode)
    {
    case 0:
      imc=max(im1,im2);
      ims=min(im1,im2);
  
      if((im1<0||im1>=nmass)||(im2<0||im2>=nmass))
	{
	  cerr<<"Error, im1="<<im1<<", im2="<<im2<<" has to be in the interval: [0,"<<nmass-1<<"]"<<endl;
	  exit(1);
	}
  
      ic=ims*nmass-(ims*(ims-1))/2+(imc-ims);
      break;
    case 1:
      if(im1>=nlights)
	{
	  cerr<<"Error, requiring uncompute combo: im1="<<im1<<" > nlights="<<nlights<<"!"<<endl;
	  exit(1);
	}
      
      ic=im1*nmass+im2;
      break;
    case 2:
      ic=0;
      if(im2<im1)
	{
	  int imt=im1;
	  im1=im2;
	  im2=imt;
	}
      if(im1>=nlights||im2<im1)
	{
	  fprintf(stderr,"Error, im1=%d im2=%d nli=%d combo mode 2 not possible\n",im1,im2,nlights);
	  exit(1);
	}
      for(int ims=0;ims<im1;ims++)
	for(int imc=ims;imc<nmass;imc++)
	  ic++;
	for(int imc=im1;imc<im2;imc++)
	  ic++;
      printf("%d %d %d\n",im1,im2,ic);
      break;
    default:
      cerr<<"Error, unkwnown mode"<<endl;
      ic=0;
      exit(1);
      break;
    }
  return ic;
}

//load all data
void load_ensembles_list(char **&base_corrs_path,char **&ens_name,int &nens,int *&T,int *&ibeta,int *&nmass,double **&mass,int *&iml_un,int *&nlights,const char *ens_list_path)
{
  FILE *ens_list_file=open_file(ens_list_path,"r");
  
  //load the number of ensembles
  read_formatted_from_file_expecting((char*)(&nens),ens_list_file,"%d","nens");
  
  //base address of all correlations
  char base_ens_path[1024];
  read_formatted_from_file_expecting((char*)(&base_ens_path),ens_list_file,"%s","base_ens_path");
  
  //ensemble parameters
  base_corrs_path=(char**)malloc(nens*sizeof(char*));
  ens_name=(char**)malloc(nens*sizeof(char*));
  T=(int*)malloc(nens*sizeof(int));
  ibeta=(int*)malloc(nens*sizeof(int));
  nmass=(int*)malloc(nens*sizeof(int));
  mass=(double**)malloc(nens*sizeof(double*));
  iml_un=(int*)malloc(nens*sizeof(int));
  nlights=(int*)malloc(nens*sizeof(int));
  
  //Loop over all ensemble. Ensemble specification is supposed to be stored in a file named [base_ens_path]/[ens_path]/data_list
  for(int iens=0;iens<nens;iens++)
    {
      //allocate ensemble name and correlations path
      base_corrs_path[iens]=(char*)malloc(1024);
      ens_name[iens]=(char*)malloc(1024);
      
      //load ensemble name and parameters
      read_formatted_from_file(ens_name[iens],ens_list_file,"%s","ens_name");
      read_ensemble_pars(base_corrs_path[iens],T[iens],ibeta[iens],nmass[iens],mass[iens],iml_un[iens],nlights[iens],combine("%s/%s/data_list",base_ens_path,ens_name[iens]).c_str());
    }
  
  fclose(ens_list_file);
}

void boot_from_jack(boot &out,jack in,int *iboot_jack)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
  out.data[nboot]=in.data[njack];
}
