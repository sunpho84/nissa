#include <include.h>
#include <iostream>
#include <sstream>

using namespace std;

int nmass;
int nsm_lev,njack;
double *mass;
int *sm_steps;
int T,TH;
char base_path[1024];

jack f(jack Z,jack M,double m1,double m2)
{
  return Z/(M*sinh(M))*(m1+m2);
}

double fun_fit(double Z1,double Z2,double M,int t)
{
  return Z1*Z2*exp(-M*TH)*cosh(M*(TH-t))/M;
}

int icombo(int im1,int im2,int r1,int r2,int ri)
{
  return ri+2*(r1+2*(r2+2*(im2+nmass*im1)));
}

jvec load_corr(int sm_lev_so,int sm_lev_si,int im1,int im2,int ri)
{
  char path[1024];
  sprintf(path,"%s_%02d_%02d",base_path,sm_steps[sm_lev_so],sm_steps[sm_lev_si]);
  int ic1=icombo(im1,im2,0,0,ri);
  int ic2=icombo(im1,im2,1,1,ri);
  
  //printf("Path: %s, im1=%d,im2=%d,sm_so=%d,sm_si=%d, combos: %d %d\n",path,im1,im2,sm_lev_so,sm_lev_si,ic1,ic2);
  
  return ((jvec_load(path,T,njack,ic1)+jvec_load(path,T,njack,ic2))/2).simmetrized(1);
}

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  
  read_formatted_from_file_expecting(base_path,fin,"%s","base_path");
  printf("Base_Path: %s\n",base_path);
  read_formatted_from_file_expecting((char*)&nmass,fin,"%d","nmass");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)(&mass[imass]),fin,"%lg","mass");
  read_formatted_from_file_expecting((char*)&njack,fin,"%d","njack");
  read_formatted_from_file_expecting((char*)&nsm_lev,fin,"%d","nsm_lev");
  sm_steps=(int*)malloc(sizeof(int)*nsm_lev);
  for(int ism_lev=0;ism_lev<nsm_lev;ism_lev++) read_formatted_from_file((char*)(&sm_steps[ism_lev]),fin,"%d","sm_steps");
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  TH=T/2;
}

int main()
{
  read_pars("input");
  
  int im1=0;
  for(int im2=0;im2<nmass;im2++)
    for(int sm_lev_so=0;sm_lev_so<nsm_lev;sm_lev_so++)
      for(int sm_lev_si=0;sm_lev_si<nsm_lev;sm_lev_si++)
	{
	  //load the current combo
	  jvec c=load_corr(sm_lev_so,sm_lev_si,im1,im2,0);
	  c.print_to_file(combine("out_%02d_%02d_%02d_%02d.xmg",im1,im2,sm_steps[sm_lev_so],sm_steps[sm_lev_si]).c_str());
	  effective_mass(c).print_to_file(combine("eff_%02d_%02d_%02d_%02d.xmg",im1,im2,sm_steps[sm_lev_so],sm_steps[sm_lev_si]).c_str());
	}
  
  return 0;
}
