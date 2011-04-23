#pragma once

#include "include.h"

int T,TH,L;
int nmoms,nmass;
int njack;
int ibeta;

char base_path[1024];
double *mass,*theta;
int *iopp_th;

double lat[4]={1/2.0198,1/2.3286,1/2.9419,1/3.6800};
double Zp[4]={0.411,0.437,0.477,0.501};

void read_input()
{
  FILE *input=open_file("input","r");
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  L=TH=T/2;
  
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");
  
  read_formatted_from_file_expecting((char*)&nmoms,input,"%d","nmoms");
  expect_string_from_file(input,"theta_list");
  
  theta=(double*)malloc(sizeof(double)*nmoms);
  for(int imom=0;imom<nmoms;imom++) read_formatted_from_file((char*)&(theta[imom]),input,"%lg","theta");
  
  //find opposite theta
  iopp_th=(int*)malloc(sizeof(int)*nmoms);
  for(int ik=0;ik<nmoms;ik++)
    {
      int ik_opp=0;
      while(theta[ik_opp]!=-theta[ik]) ik_opp++;
      iopp_th[ik]=ik_opp;
    }
  
  fclose(input);
}

void check_interval(int var,int min,int max)
{
  if(var>=max || var<min)
    {
      cerr<<"Asked for correlation of impossible combination, "<<var<<" not in the interval: ["<<min<<","<<max<<")"<<endl;
      exit(1);
    }
}

//read a particular two point passed as argument
jvec read_two_points(const char *in,int im1,int im2,int ik1,int ik2,int r1,int r2,int ri)
{
  check_interval(ik1,0,nmoms);check_interval(im1,0,nmass);check_interval(r1,0,2);
  check_interval(ik2,0,nmoms);check_interval(im2,0,nmass);check_interval(r2,0,2);
  check_interval(ri,0,2);
  
  int icorr=2*(ik1+nmoms*(ik2+nmoms*(im1*2+r1+2*nmass*(im2*2+r2))))+ri;
  return -jvec_load_naz(in,T,njack,icorr)/(L*L*L);
}
//read a particular three point passed as argument
jvec read_three_points(const char *in,int im1,int im2,int im3,int ik1,int ik2,int r1,int r2,int r3,int mu,int ri)
{
  jvec c(T,njack);

  check_interval(ik1,0,nmoms);check_interval(im1,0,nmass);check_interval(r1,0,2);
  check_interval(ik2,0,nmoms);check_interval(im2,0,nmass);check_interval(r2,0,2);
  check_interval(im3,0,nmass);check_interval(r3,0,2);
  check_interval(ri,0,2);
  check_interval(mu,0,4);

  int icorr=ik1+nmoms*(ik2+nmoms*(im1*2+r1+2*nmass*(im2*2+r2+2*nmass*(im3*2+r3))));
  
  FILE *file=open_file(in,"r");

  if(fseek(file,icorr*2*4*sizeof(double)*T*(njack+1),SEEK_SET))
    {
      cerr<<"Error while searching for correlation"<<icorr<<endl;
      exit(1);
    }

  double data_in[T][4][2][njack+1];
  int stat=fread(data_in,sizeof(double)*4*2*T*(njack+1),1,file);
  if(stat!=1)
    {
      if(stat==EOF)
	{
	  cerr<<"Error, reached EOF while reading data!"<<endl;
	  exit(1);
	}
      else
        {
          perror("Error while reading data");
          exit(1);
        }
    }
  
  for(int iel=0;iel<T;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c.data[iel].data[ijack]=data_in[iel][mu][ri][ijack];
  
  fclose(file);
  
  return c/(L*L*L);
}

jvec read_P5_Vmu_P5(int im1,int im2,int im3,int ik1,int ik2,int r1,int r2,int r3,int mu)
{
  char path[1024];sprintf(path,"%s/oPVmuPo-sss_conf.1.dat",base_path);
  
  //real or imag part, according to mu
  int ri[4]={0,1,1,1};
  
  return read_three_points(path,im1,im2,im3,ik1,ik2,r1,r2,r3,mu,ri[mu]);
}

jvec read_P5_P5(int im1,int im2,int ik1,int ik2,int r1,int r2)
{
  char path[1024];sprintf(path,"%s/oPPo-ss_conf.1.dat",base_path);
  //read real part
  int ri=0;return read_two_points(path,im1,im2,ik1,ik2,r1,r2,ri);
}

jvec read_A0_P5(int im1,int im2,int ik1,int ik2,int r1,int r2)
{
  char path[1024];sprintf(path,"%s/oA0Po-ss_conf.1.dat",base_path);
  //read real part
  int ri=0;return read_two_points(path,im1,im2,ik1,ik2,r1,r2,ri);
}

jvec read_thimpr_P5_P5(int im1,int im2,int ik1,int ik2,int r1,int r2)
{return (read_P5_P5(im1,im2,ik1,ik2,r1,r2)+read_P5_P5(im1,im2,iopp_th[ik1],iopp_th[ik2],r1,r2))/2;}

jvec read_ch_thimpr_P5_P5(int im1,int im2,int ik1,int r1)
{return read_thimpr_P5_P5(im1,im2,ik1,0,r1,r1);}

jvec read_thimpr_P5_Vmu_P5(int im1,int im2,int im3,int ik1,int ik2,int r1,int r2,int r3,int mu)
{
  //summ or subtract, according to mu
  int si[4]={+1,-1,-1,-1};
  return (read_P5_Vmu_P5(im1,im2,im3,ik1,ik2,r1,r2,r3,mu)+si[mu]*
	  read_P5_Vmu_P5(im1,im2,im3,iopp_th[ik1],iopp_th[ik2],r1,r2,r3,mu))/2;
}

jvec read_ch_thimpr_P5_Vmu_P5(int im1,int im2,int im3,int ik1,int ik2,int r,int mu)
{return read_thimpr_P5_Vmu_P5(im1,im2,im3,ik1,ik2,r,!r,r,mu);}
  
jvec read_deg_ch_thimpr_P5_Vmu_P5(int im_spec,int im_val,int ik1,int ik2,int r,int mu)
{return read_ch_thimpr_P5_Vmu_P5(im_spec,im_val,im_val,ik1,ik2,r,mu);}

jvec read_deg_ch_thimpr_P5_V0_P5(int im_spec,int im_val,int ik1,int ik2,int r)
{return read_deg_ch_thimpr_P5_Vmu_P5(im_spec,im_val,ik1,ik2,r,0);}

  
jack calculate_Q2(jack E1,double theta1,jack E2,double theta2)
{return pow(E1-E2,2)-3*pow(2*M_PI*(theta1-theta2)/L,2);}

void calculate_Q_P(jack *Q,jack *P,jack E1,double theta1,jack E2,double theta2)
{
  Q[0]=E1-E2;
  P[0]=E1+E2;

  for(int i=1;i<4;i++)
    {
      Q[i]=jack(njack);Q[i]=2*M_PI/L*(theta1-theta2);
      P[i]=jack(njack);P[i]=2*M_PI/L*(theta1+theta2);
    }
}

jack quad_jack_prod_quad_jack(jack *v,jack *p)
{
  jack res=v[0]*p[0];
  for(int mu=1;mu<4;mu++) res-=v[mu]*p[mu];
  return res;
}
