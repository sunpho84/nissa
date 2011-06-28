#pragma once

#include "include.h"

const double hc=0.19733;

bool latpars_initted=false;

//renormalization constants and lattice spacings
int **ijack_boot;
double lat_med[4],lat_med_fm[4];
boot lat[4],Zp[4],Za[4];
//results taken by arxiv:1004.1115
double Za_med[4]={0.746,0.746,0.772,0.780};
double Za_err[4]={0.011,0.006,0.006,0.006};
//l physics
boot ml_phys,aml_phys[4];
boot f0,db0;
double fPi_phys=0.1307;
//s physics
boot ms_phys,ams_phys[4];
double mK_phys=0.493667;
double dmK_phys_med=-0.006;
double dmK_phys_err=0.0006/1000;
//c physics
boot mc_phys,amc_phys[4];
double mD_phys=1.8696;
double mDs_phys=1.9685;
//lattice tags
const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

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
  
  //now it's time to shrink the variance!
  for(int iboot=0;iboot<nboot;iboot++)
    out.data[iboot]=out.data[nboot]+(out.data[iboot]-out.data[nboot])*sqrt(15.0/9);
  
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
      lat_med_fm[ib]=lat_med[ib]*hc;
      
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

void init_ijack_boot(int nel,int nboot,int njack,int seed)
{
  ijack_boot=(int**)malloc(sizeof(int)*nel);
  for(int iel=0;iel<nel;iel++)
    {
      ijack_boot[iel]=(int*)malloc(sizeof(int)*(nboot+1));
      ran_gen estr(seed+iel);
      for(int iboot=0;iboot<nboot;iboot++) ijack_boot[iel][iboot]=estr.get_int(njack);
      ijack_boot[iel][nboot]=njack;
    }   
}

void boot_from_jack(boot &out,jack &in,int iel)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++)
    out.data[iboot]=in.data[ijack_boot[iel][iboot]];
  out.data[nboot]=in.data[njack];
}

