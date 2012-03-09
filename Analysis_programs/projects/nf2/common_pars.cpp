#pragma once

#include "include.h"

const int nboot=100;
int njack=16;

const double hc=0.19733;

//results taken by arxiv:1004.1115
double Zp_med[4]={0.411,0.437,0.477,0.501};
double Zp_err[4]={0.012,0.007,0.006,0.020};
double Za_med[4]={0.746,0.746,0.772,0.780};
double Za_err[4]={0.011,0.006,0.006,0.006};
double Zv_med[4]={0.5816,0.6103,0.6451, 0.6751};
double Zv_err[4]={0.0002,0.0003,0.0003, 0.0003};
double Zp_fr_Zs_med[4]={0.580,0.626,0.686, 0.746};
double Zp_fr_Zs_err[4]={0.017,0.013,0.012, 0.011};

double ml_phys_med=3.6e-03,ml_phys_err=0.2e-03;
double ms_phys_med=95e-03,ms_phys_err=6e-03;
double mc_phys_med=1.14,mc_phys_err=0.04;

double lat_med[4]={0.486508,0.422773,0.335339,0.268402};
double lat_err[4]={0.0151488,0.0112503,0.00766779,0.00570585};

double f0_med=0.121806,f0_err=0.000162696;
double db0_med=5.35229,db0_err=0.205398;

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

int nref_hmass=11;
bvec ref_hmass(11,nboot,njack);
double ref_ahmass[11]={0.1828,0.2150,0.2529,0.2974,0.3498,0.4114,0.4839,0.5691,0.6694,0.7873,0.9260};

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
  
  /*
  ml_phys=boot(nboot,njack);
  ms_phys=boot(nboot,njack);
  mc_phys=boot(nboot,njack);
  ml_phys.fill_gauss(ml_phys_med,ml_phys_err,83748245);
  ms_phys.fill_gauss(ms_phys_med,ms_phys_err,847329);
  mc_phys.fill_gauss(mc_phys_med,mc_phys_err,9834654);

  f0=boot(nboot,njack);
  db0=boot(nboot,njack);
  f0.fill_gauss(f0_med,f0_err,86546);
  db0.fill_gauss(db0_med,db0_err,89349721);
  
  for(int ib=0;ib<4;ib++)
    {
      lat[ib]=boot(nboot,njack);
      lat[ib].fill_gauss(lat_med[ib],lat_err[ib],4938724+ib);
      
      Zp[ib]=boot(nboot,njack);
      Za[ib]=boot(nboot,njack);
      Zv[ib]=boot(nboot,njack);
      Zp_fr_Zs[ib]=boot(nboot,njack);
      Zp[ib].fill_gauss(Zp_med[ib],Zp_err[ib],4329852+ib);
      Za[ib].fill_gauss(Za_med[ib],Za_err[ib],2873246+ib);
      Zv[ib].fill_gauss(Zv_med[ib],Zv_err[ib],4334943+ib);
      Zp_fr_Zs[ib].fill_gauss(Zp_fr_Zs_med[ib],Zp_fr_Zs_err[ib],5486357249+ib);
      
      aml_phys[ib]=ml_phys*lat[ib]*Zp[ib];
      ams_phys[ib]=ms_phys*lat[ib]*Zp[ib];
      amc_phys[ib]=mc_phys*lat[ib]*Zp[ib];
    }      
  
  for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++)
    ref_hmass[iref_hmass]=ref_ahmass[iref_hmass]/lat[1].med()/Zp[1].med();
  */
    
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

  for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++)
    ref_hmass[iref_hmass]=ref_ahmass[iref_hmass]/lat[1].med()/Zp[1].med();
  
  latpars_initted=true;
}
