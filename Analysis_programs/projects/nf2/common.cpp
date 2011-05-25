#pragma once

#include "include.h"

const int nboot=100;
int njack=16;

const double hc=0.19733;
double lat_med[4]={1/2.0198,1/2.3286,1/2.9419,1/3.6800};
double lat_err[4]={lat_med[0]/31.5,lat_med[1]/36.8,lat_med[2]/41.9,lat_med[3]/44.7};
double lat_med_fm[4]={lat_med[0]*hc,lat_med[1]*hc,lat_med[2]*hc,lat_med[3]*hc};

double lat_rat_med[4]={1,0.868996,0.79319,0.800388};
double lat_rat_err[4]={0,0.00620553,0.00533782,0.00512307};

double Zp_med[4]={0.411,0.437,0.477,0.501};
double Zp_err[4]={0.012,0.007,0.006,0.020};

//results taken by arxiv:1004.1115
double Za_med[4]={0.746,0.746,0.772,0.780};
double Za_err[4]={0.011,0.006,0.006,0.006};

double ml_phys_med=3.6/1000,ms_phys_med=95.0/1000,mc_phys_med=1140.0/1000;
double ml_phys_err=0.2/1000,ms_phys_err= 6.0/1000,mc_phys_err=  40.0/1000;

double f0_med=121.9019/1000,db0_med=5.2756;
double f0_err=0.1668/1000,db0_err=0.2078;

bool latpars_initted=false;

boot lat[4],Zp[4],Za[4];
boot ml_phys,ms_phys,mc_phys;
boot aml_phys[4],ams_phys[4],amc_phys[4];
boot f0,db0;
int iboot_jack[4][nboot+1];

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

int icombo(int im1,int im2,int nmass)
{
  int imc=max(im1,im2);
  int ims=min(im1,im2);
  
  if((im1<0||im1>=nmass)||(im2<0||im2>=nmass))
    {
      cerr<<"Error, im1="<<im1<<", im2="<<im2<<" has to be in the interval: [0,"<<nmass-1<<"]"<<endl;
      exit(1);
    }
  
  return ims*nmass-(ims*(ims-1))/2+(imc-ims);
}

void init_latpars()
{
  if(latpars_initted==true) return;
  
  for(int ib=0;ib<4;ib++)
    {
      lat[ib]=boot(nboot,njack);
      if(ib==0) lat[ib].fill_gauss(lat_med[ib],lat_err[ib],156010+ib);
      else
	{
	  lat[ib].fill_gauss(lat_rat_med[ib],lat_rat_err[ib],156010+ib);
	  lat[ib]*=lat[ib-1];
	}
      cout<<"Lattice spacing "<<ib<<" "<<lat[ib]*hc<<endl;
      
      Zp[ib]=boot(nboot,njack);
      Zp[ib].fill_gauss(Zp_med[ib],Zp_err[ib],1534572+ib);
    }
  
  ml_phys=boot(nboot,njack);
  ms_phys=boot(nboot,njack);
  mc_phys=boot(nboot,njack);
  
  ml_phys.fill_gauss(ml_phys_med,ml_phys_err,2478463);
  ms_phys.fill_gauss(ms_phys_med,ms_phys_err,6764732);
  mc_phys.fill_gauss(mc_phys_med,mc_phys_err,7623647);
  
  for(int ib=0;ib<4;ib++)
    {
      ran_gen estr(22141425+ib);
      for(int iboot=0;iboot<nboot;iboot++)
	iboot_jack[ib][iboot]=estr.get_int(njack);
      iboot_jack[ib][nboot]=njack;

      aml_phys[ib]=ml_phys*lat[ib]*Zp[ib];
      ams_phys[ib]=ms_phys*lat[ib]*Zp[ib];
      amc_phys[ib]=mc_phys*lat[ib]*Zp[ib];
    }

  for(int ib=0;ib<4;ib++)
    {
      Za[ib]=boot(nboot,njack);
      Za[ib].fill_gauss(Za_med[ib],Za_err[ib],2873246+ib);
    }      
  
  f0=boot(nboot,njack);
  db0=boot(nboot,njack);

  f0.fill_gauss(f0_med,f0_err,855438428);
  db0.fill_gauss(db0_med,db0_err,76547898);
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

void boot_from_jack(boot &out,jack in,int ibeta)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[ibeta][iboot]];
  out.data[nboot]=in.data[njack];
}
