#pragma once

char base_path[1024];

int njack=15;
int T,L,TH;

double *mass,*theta;
int *iopp_th;
int ith_spec,im_spec,r_spec;

int nmass,nmoms;

int n2pts=30;
const char name_2pts[30][20]={"P5P5","A0P5","P5A0","V0V0","A0A0","V1V1","V2V2","V3V3","A1A1","A2A2","A3A3","S0S0","T1T1","T2T2","T3T3","B1B1","B2B2","B3B3","V0S0","S0V0","S0P5","P5S0","V0P5","P5V0","A0S0","S0A0","CHROMO-S0P5","CHROMO-P5P5","CHROMO-S0S0","CHROMO-P5S0"};
FILE *file_2pts;

int n3pts=43;
const char name_3pts[43][20]={"V0P5","V1P5","V2P5","V3P5","S0P5","T1P5","T2P5","T3P5","P5S0","V1V2","V2V3","V3V1","A0V0","A0V1","A0V2","A0V3","A1V0","A2V0","A3V0","A1V1","A2V2","A3V3","T1V2","T2V3","T3V1","B1V1","B2V2","B3V3","B1V2","B2V3","B3V1","A0P5","A1P5","A2P5","A3P5","P5P5","B1P5","B2P5","B3P5","CHROMO-S0P5","CHROMO-P5P5","CHROMO-S0S0","CHROMO-P5S0"};
FILE *file_3pts;

int check_interval(int a,int b,int c)
{
  if(a>=c||a<b)
    {
      fprintf(stderr,"Error, %d outside interval: [%d,%d]\n",a,b,c);
      exit(1);
    }
  return a;
}

void read_data_list()
{
  FILE *input=open_file("data_list","r");
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  L=TH=T/2;
  
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");
  
  read_formatted_from_file_expecting((char*)&nmoms,input,"%d","nmoms");
  expect_string_from_file(input,"theta_list");
  theta=(double*)malloc(sizeof(double)*nmoms);
  for(int imom=0;imom<nmoms;imom++) read_formatted_from_file((char*)&(theta[imom]),input,"%lg","theta");

  read_formatted_from_file_expecting((char*)&im_spec,input,"%d","im_spec");
  read_formatted_from_file_expecting((char*)&ith_spec,input,"%d","ith_spec");
  read_formatted_from_file_expecting((char*)&r_spec,input,"%d","r_spec");
  
  //find opposite theta
  iopp_th=(int*)malloc(sizeof(int)*nmoms);
  for(int ik=0;ik<nmoms;ik++)
    {
      int ik_opp=0;
      while(theta[ik_opp]!=-theta[ik]) ik_opp++;
      iopp_th[ik]=ik_opp;
    }
  
  fclose(input);
  
  file_2pts=open_file(combine("%s/two_points",base_path).c_str(),"r");
  file_3pts=open_file(combine("%s/three_points",base_path).c_str(),"r");
}

jvec read_two_points(int ith2,int im2,int r2,int im1,int r1,int iobs,int ri=0)
{
  int icorr=check_interval(ith2,0,nmoms);
  icorr=icorr*nmass+check_interval(im2,0,nmass);
  icorr=icorr*2+check_interval(r2,0,2);
  icorr=icorr*nmass+check_interval(im1,0,nmass);
  icorr=icorr*2+check_interval(r1,0,2);
  icorr=icorr*n2pts+check_interval(iobs,0,n2pts);
  icorr=icorr*2+check_interval(ri,0,2);

  jvec out(T,njack);
  out.load(file_2pts,icorr);
  
  return out;
}

jvec read_ch_thimpr_P5_P5(int ith2,int im2,int r2,int im1)
{return (read_two_points(ith2,im2,r2,im1,r2,0,0)+read_two_points(iopp_th[ith2],im2,r2,im1,r2,0,0))/2;}

jvec read_three_points(int ith2,int im2,int ith1,int im1,int iobs,int ri=0)
{
  int icorr=check_interval(ith2,0,nmoms);
  icorr=icorr*nmass+check_interval(im2,0,nmass);
  icorr=icorr*nmoms+check_interval(ith1,0,nmoms);
  icorr=icorr*nmass+check_interval(im1,0,nmass);
  icorr=icorr*n3pts+check_interval(iobs,0,n2pts);
  icorr=icorr*2+check_interval(ri,0,2);
  
  jvec out(T,njack);
  out.load(file_3pts,icorr);
  
  return out;
}

jvec read_ch_thimpr_P5_Vmu_P5(int ith2,int im2,int ith1,int im1,int mu)
{
  //real or imag part, according to mu
  int ri[4]={0,1,1,1};
  int iobs=mu;
  
  return (read_three_points(ith2,im2,ith1,im1,iobs,ri[mu])+read_three_points(iopp_th[ith2],im2,iopp_th[ith1],im1,iobs,ri[mu]))/2;
}

void calculate_Q_P(jack *Q,jack *P,jack E1,double theta1,jack E2,double theta2)
{
  Q[0]=E1-E2;
  P[0]=E1+E2;
  
  for(int i=1;i<4;i++)
    {
      Q[i]=jack(njack);Q[i]=M_PI/L*(theta1-theta2);
      P[i]=jack(njack);P[i]=M_PI/L*(theta1+theta2);
    }
}

jack quad_jack_prod_quad_jack(jack *v,jack *p)
{
  jack res=v[0]*p[0];
  for(int mu=1;mu<4;mu++) res-=v[mu]*p[mu];
  return res;
}
