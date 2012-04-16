#include "include.h"

char data_list_file[1024];
char corr_name[1024];
char out_file[1024];

char base_path[1024];
int T,L,TH;
int ntheta,nmass;
int njack=16;
int tmin,tmax;

void read_ensemble_pars(char *base_path,int &T,int &nmass,int &ntheta,const char *data_list_file)
{
  FILE *input=open_file(data_list_file,"r");
  
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  L=TH=T/2;
  int ibeta;
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  expect_string_from_file(input,"mass_list");
  double mass;
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass),input,"%lg","mass");
  
  read_formatted_from_file_expecting((char*)&ntheta,input,"%d","ntheta");
  expect_string_from_file(input,"theta_list");
  double theta;
  for(int itheta=0;itheta<ntheta;itheta++) read_formatted_from_file((char*)&(theta),input,"%lg","theta");
  
  fclose(input);
}

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  read_formatted_from_file_expecting(data_list_file,fin,"%s","data_list_file");
  
  read_formatted_from_file_expecting(corr_name,fin,"%s","corr_name");
  expect_string_from_file(fin,"tint");
  read_formatted_from_file((char*)&tmin,fin,"%d","tmin");
  read_formatted_from_file((char*)&tmax,fin,"%d","tmax");
  
  read_formatted_from_file_expecting(out_file,fin,"%s","out_file");
    
  fclose(fin);
}

int icombo(int ith2,int im1,int im2,int r1,int r2,int ri)
{return ri+2*(r1+2*(im1+nmass*(r2+2*(im2+nmass*ith2))));}

jvec load_standing_tm_pion()
{
  int ic1=icombo(3,0,0,0,0,0),ic2=icombo(3,0,0,1,1,0);
  
  cout<<"Combo: "<<ic1<<" "<<ic2<<endl;
  return (jvec_load(combine("%s/%s",base_path,corr_name).c_str(),T,njack,ic1)+
	  jvec_load(combine("%s/%s",base_path,corr_name).c_str(),T,njack,ic2))/2;
}

//function to fit
double fun_fit(double Z2,double M,int t)
{return Z2*exp(-M*TH)*cosh(M*(TH-t))/sinh(M);}

//chi2 calculation
double *corr_fit,*corr_err;
void chi2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  for(int t=tmin;t<=min(tmax,TH);t++)
    ch+=sqr((corr_fit[t]-fun_fit(p[0],p[1],t))/corr_err[t]);
}

void two_pts_fit_minuit(jack &M,jack &Z,jvec corr)
{
  //define minuit staff
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(chi2);
  corr_fit=new double[TH+1];
  corr_err=new double[TH+1];
  
  minu.DefineParameter(0,"Z2",Z[0],Z.err(),0,2*Z[0]);
  minu.DefineParameter(1,"M",M[0],M.err(),0,2*M[0]);
  for(int t=tmin;t<=tmax;t++) corr_err[t]=corr.data[t].err();
  
  //jacknife analysis
  for(int ijack=0;ijack<njack+1;ijack++)
    {
      //copy data so that glob function may access it
      for(int t=tmin;t<=tmax;t++) corr_fit[t]=corr.data[t].data[ijack];
      
      //fit
      double dum;
      minu.Migrad();          
      minu.GetParameter(0,Z.data[ijack],dum);
      minu.GetParameter(1,M.data[ijack],dum);
    }
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s input",arg[0]);
  
  read_pars(arg[1]);
  read_ensemble_pars(base_path,T,nmass,ntheta,data_list_file);
  
  jvec corr=load_standing_tm_pion().simmetrized(1);
  
  jack M,Z2;
  two_pts_fit(M,Z2,corr,tmin,tmax,"M_2steps.xmg","Z2_2steps.xmg");
  cout<<"M="<<M<<", Z2="<<Z2<<endl;
  
  //two_pts_fit_minuit(M,Z2,corr);
  //cout<<"M="<<M<<", Z2="<<Z2<<endl;
  
  M.write_to_binfile(out_file);
  Z2.append_to_binfile(out_file);
  
  return 0;
}
