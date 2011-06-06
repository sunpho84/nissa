#pragma once

#include "include.h"

const double hc=0.19733;

//results taken by arxiv:1004.1115
double Za_med[4]={0.746,0.746,0.772,0.780};
double Za_err[4]={0.011,0.006,0.006,0.006};

bool latpars_initted=false;

double lat_med[4],lat_med_fm[4];
boot lat[4],Zp[4],Za[4];
boot ml_phys,ms_phys,mc_phys;
boot aml_phys[4],ams_phys[4],amc_phys[4];
boot f0,db0;

double fPi_phys=0.1307;

double mK_phys=0.493667;
double delta_mK_phys_med=-0.006;
double delta_mK_phys_err=0.0006;

double mD_phys=1.8696;
double mD_s_phys=1.9685;

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
  
  FILE *input_latpars=fopen("/home/francesco/QCD/LAVORI/NF2/DATA/latpars_E","r");
  
  for(int ib=0;ib<4;ib++) lat[ib]=get_latpar(input_latpars);
  for(int ib=0;ib<4;ib++) Zp[ib]=get_latpar(input_latpars);
  
  f0=get_latpar(input_latpars);
  db0=get_latpar(input_latpars);
  
  ml_phys=get_latpar(input_latpars);
  ms_phys=get_latpar(input_latpars);
  mc_phys=get_latpar(input_latpars);
  
  for(int ib=0;ib<4;ib++)
    {
      lat_med[ib]=lat[ib].med();
      lat_med_fm[ib]=lat_med[ib]/hc;
      
      aml_phys[ib]=ml_phys*lat[ib]*Zp[ib];
      ams_phys[ib]=ms_phys*lat[ib]*Zp[ib];
      amc_phys[ib]=mc_phys*lat[ib]*Zp[ib];
    }

  for(int ib=0;ib<4;ib++)
    {
      Za[ib]=boot(nboot,njack);
      Za[ib].fill_gauss(Za_med[ib],Za_err[ib],2873246+ib);
    }      
}
