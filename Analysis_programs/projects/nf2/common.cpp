#pragma once

#include "include.h"

const int nboot=100;
int njack=16;

const double hc=0.19733;
double lat_med[4]={1/2.0198,1/2.3286,1/2.9419,1/3.6800};
double lat_err[4]={lat_med[0]/31.5,lat_med[1]/36.8,lat_med[2]/41.9,lat_med[3]/44.7};
double lat_med_fm[4]={lat_med[0]*hc,lat_med[1]*hc,lat_med[2]*hc,lat_med[3]*hc};

double Zp_med[4]={0.411,0.437,0.477,0.501};
double Zp_err[4]={0.012,0.007,0.006,0.020};

//results taken by arxiv:1004.1115
double Za_med[4]={0.746,0.746,0.772,0.780};
double Za_err[4]={0.011,0.006,0.006,0.006};

bool latpars_initted=false;

double ml_phys_med=3.6/1000,ms_phys_med=95.0/1000,mc_phys_med=1140.0/1000;
double ml_phys_err=0.2/1000,ms_phys_err= 6.0/1000,mc_phys_err=  40.0/1000;

boot lat[4],Zp[4],Za[4];
boot ml_phys,ms_phys,mc_phys;
boot aml_phys[4],ams_phys[4],amc_phys[4];
int iboot_jack[4][nboot+1];

void read_ensamble_pars(char *base_path,int &T,int &ibeta,int &nmass,double *&mass,int &iml_un,int &nlights,const char *data_list_file)
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
      lat[ib].fill_gauss(lat_med[ib],lat_err[ib],1010+ib);

      Zp[ib]=boot(nboot,njack);
      Zp[ib].fill_gauss(Zp_med[ib],Zp_err[ib],1532+ib);
    }
  
  ml_phys=boot(nboot,njack);
  ms_phys=boot(nboot,njack);
  mc_phys=boot(nboot,njack);

  ml_phys.fill_gauss(ml_phys_med,ml_phys_err,2478463);
  ms_phys.fill_gauss(ms_phys_med,ms_phys_err,6764732);
  mc_phys.fill_gauss(mc_phys_med,mc_phys_err,7623647);
  
  for(int ibeta=0;ibeta<4;ibeta++)
    {
      ran_gen estr(22141425);
      for(int iboot=0;iboot<nboot;iboot++)
	iboot_jack[ibeta][iboot]=estr.get_int(njack);
      iboot_jack[ibeta][nboot]=njack;

      aml_phys[ibeta]=ml_phys*lat[ibeta]*Zp[ibeta];
      ams_phys[ibeta]=ms_phys*lat[ibeta]*Zp[ibeta];
      amc_phys[ibeta]=mc_phys*lat[ibeta]*Zp[ibeta];
    }

  for(int ib=0;ib<4;ib++)
    {
      Za[ib]=boot(nboot,njack);
      Za[ib].fill_gauss(Za_med[ib],Za_err[ib],2323456);
    }      

}

void boot_from_jack(boot &out,jack in,int ibeta)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[ibeta][iboot]];
  out.data[nboot]=in.data[njack];
}
