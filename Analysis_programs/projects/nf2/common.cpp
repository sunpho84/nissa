#pragma once

#include "include.h"

const int nboot=100;
int njack=16;

const double hc=0.19733;

//results taken by arxiv:1004.1115
double Za_med[4]={0.746,0.746,0.772,0.780};
double Za_err[4]={0.011,0.006,0.006,0.006};
double Zv_med[4]={0.5816,0.6103,0.6451};
double Zv_err[4]={0.0002,0.0003,0.0003};
double Zp_fr_Zs_med[4]={0.580,0.626,0.686, 0.746};
double Zp_fr_Zs_err[4]={0.017,0.013,0.012, 0.011};

bool latpars_initted=false;

boot lat[4],Zp[4],Za[4],Zv[4],Zp_fr_Zs[4];
boot ml_phys,ms_phys,mc_phys;
boot aml_phys[4],ams_phys[4],amc_phys[4];
boot f0,db0;

double mPi_phys=0.1350; //neutral
double fPi_phys=0.1307;

double mK_phys=0.493667;
double delta_mK_phys_med=-0.006;
double delta_mK_phys_err=0.0006;

double mD_phys=1.8696;
double mD_s_phys=1.9685;

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

boot get_latpar(FILE *fin)
{
  boot out(nboot,njack);
  double buf[101];
  int nr=fread(buf,sizeof(double),101,fin);
  if(nr!=101)
    {
      perror("Error reading from latpars file");
      exit(1);
    }
  out.put(buf);
  return out;
}

void init_latpars()
{
  if(latpars_initted==true) return;
  
  FILE *input_latpars=fopen("/Users/francesco/QCD/LAVORI/NF2/latpars_E","r");
  
  for(int ib=0;ib<4;ib++) lat[ib]=get_latpar(input_latpars);
  for(int ib=0;ib<4;ib++) Zp[ib]=get_latpar(input_latpars);
  
  f0=get_latpar(input_latpars);
  db0=get_latpar(input_latpars);
  
  ml_phys=get_latpar(input_latpars);
  ms_phys=get_latpar(input_latpars);
  mc_phys=get_latpar(input_latpars);
  
  for(int ib=0;ib<4;ib++)
    {
      aml_phys[ib]=ml_phys*lat[ib]*Zp[ib];
      ams_phys[ib]=ms_phys*lat[ib]*Zp[ib];
      amc_phys[ib]=mc_phys*lat[ib]*Zp[ib];
    }

  for(int ib=0;ib<4;ib++)
    {
      Za[ib]=boot(nboot,njack);
      Zv[ib]=boot(nboot,njack);
      Zp_fr_Zs[ib]=boot(nboot,njack);
      Za[ib].fill_gauss(Za_med[ib],Za_err[ib],2873246+ib);
      Zv[ib].fill_gauss(Zv_med[ib],Zv_err[ib],4334943+ib);
      Zp_fr_Zs[ib].fill_gauss(Zp_fr_Zs_med[ib],Zp_fr_Zs_err[ib],5486357249+ib);
    }      
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
