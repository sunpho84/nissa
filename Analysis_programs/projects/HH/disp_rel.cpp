#include <include.h>
#include <iostream>
#include <sstream>

using namespace std;

int nmass,ntheta;
int njack;
double *mass;
double *theta;
int T,TH,L;
char base_path[1024];
char corr_name[1024],out_file[1024];
int tmin[2],tmax[2];

double *X_fit;
double *Y_fit;
double *err_Y_fit;
int nc_fit;

jack fixed_A;

int *use;

double parab(double a,double b,double c,double x)
{
  return a*x*x+b*x+c;
}

//calculate the chi square
double chi2(double A,double B,double C,double D,double *a)
{
  double ch2=0;
  
  for(int iel=0;iel<nc_fit;iel++)
    if(use[iel])
      ch2+=pow((Y_fit[iel]-parab(A,B,C,X_fit[iel]))/err_Y_fit[iel],2);
  
  return ch2;
}

//wrapper for the calculation of the chi2
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0],B=p[1],C=p[2],D=p[3];
  double *a=p+4;
  ch=chi2(A,B,C,D,a);
}


void parab_fit(double ch2,jack &a,jack &b,jack &c,double *X,jvec Y,const char *path,bool fl_parab=1)
{
  ofstream out(path);
  out<<"@type xydy\n";
  
  double MX=0;
  //copy X
  nc_fit=Y.nel;
  X_fit=new double[nc_fit];
  for(int iel=0;iel<nc_fit;iel++)
    if(use[iel])
    {
      X_fit[iel]=X[iel];
      out<<X[iel]<<" "<<Y[iel]<<endl;
      if(X[iel]>MX) MX=X[iel];
    }
  Y_fit=new double[nc_fit];
  err_Y_fit=new double[nc_fit];
  
  TMinuit minu;
  minu.SetPrintLevel(0);
  
  int npars=3;
  minu.DefineParameter(0,"A",0.0,0.0001,0,0);
  minu.DefineParameter(1,"B",0.0,0.0001,0,0);
  minu.DefineParameter(2,"C",0.0,0.0001,0,0);
  if(!fl_parab)
    {
      minu.FixParameter(0);
      npars--;
    }
  minu.SetFCN(chi2_wr);
  
  for(int ijack=0;ijack<=njack;ijack++)
    {
      minu.DefineParameter(2,"C",fixed_A[ijack],0.0001,0,0);
      minu.FixParameter(2);
      if(ijack>0)
        minu.SetPrintLevel(-1);
      
      for(int iel=0;iel<nc_fit;iel++)
        {
          Y_fit[iel]=Y.data[iel].data[ijack];
          err_Y_fit[iel]=Y[iel].err();
        }
      
      //minimize
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,a.data[ijack],dum);
      minu.GetParameter(1,b.data[ijack],dum);
      minu.GetParameter(2,c.data[ijack],dum);
    }
  
  out<<"&\n@type xy\n";
  for(double t=0;t<MX;t+=MX/100)
    out<<t<<" "<<parab(a[njack],b[njack],c[njack],t)<<endl;
  
  out.close();
}

int icombo(int ith1,int ith2,int im1,int im2,int r1,int r2,int ri)
{
  int ic=ri+2*(r1+2*(im1+nmass*(r2+2*(im2+nmass*(ith2+ntheta*ith1)))));

  return ic;
}

jvec load_corr(char *base_path,int sm_lev,int ith1,int ith2,int im1,int im2,int r1,int r2,const char *obs_name)
{
  char path[1024];
  int tag[2]={0,30};
  sprintf(path,"%s/2pts_%s_30_%02d",base_path,obs_name,tag[sm_lev]);
  cout<<endl<<"ith1:"<<ith1<<" ith2:"<<ith2<<" im1:"<<im1<<" im2:"<<im2<<endl;
  jvec l=(jvec_load(path,T,njack,icombo(ith1,ith2,im1,im2,!r1,!r2,0))+
	  jvec_load(path,T,njack,icombo(ith1,ith2,im1,im2, r1, r2,0)))/2;
  
  jvec r=(jvec_load(path,T,njack,icombo(ith2,ith1,im2,im1,!r1,!r2,0))+
	  jvec_load(path,T,njack,icombo(ith2,ith1,im2,im1, r1, r2,0)))/2;
  
  return (l+r)/2;
}

void read_ensemble_pars(const char *data_list_file)
{
  FILE *input=open_file(data_list_file,"r");
  
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&njack,input,"%d","njack");
  int ibeta;
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");
  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");

  read_formatted_from_file_expecting((char*)&ntheta,input,"%d","ntheta");
  expect_string_from_file(input,"theta_list");
  theta=(double*)malloc(sizeof(double)*ntheta);
  for(int itheta=0;itheta<ntheta;itheta++) read_formatted_from_file((char*)&(theta[itheta]),input,"%lg","theta");
  
  fclose(input);
}

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  char data_list_file[256];
  read_formatted_from_file_expecting(data_list_file,fin,"%s","data_list_file");

  read_formatted_from_file_expecting(corr_name,fin,"%s","corr_name");
  
  read_formatted_from_file_expecting((char*)(&tmin[0]),fin,"%d","tint_sl");
  read_formatted_from_file((char*)(&tmax[0]),fin,"%d","tint_sl");
  
  read_formatted_from_file_expecting((char*)(&tmin[1]),fin,"%d","tint_ss");
  read_formatted_from_file((char*)(&tmax[1]),fin,"%d","tint_ss");
  
  read_formatted_from_file_expecting(out_file,fin,"%s","out_file");
  
  fclose(fin);
  
  read_ensemble_pars(data_list_file);
  
  TH=L=T/2;
}

int icombo_int(int ith1,int ith2,int im1,int im2)
{
  return im1+nmass*(im2+nmass*(ith2+ntheta*ith1));
}

int main(int narg,char **arg)
{
  read_pars("input");
  use=new int[nc_fit];
  jvec c;
  //loop over heavier mass
  int ncombo=nmass*nmass*ntheta*ntheta;
  jvec aM(ncombo,njack);
  for(int ith1=0;ith1<ntheta;ith1++)
    for(int ith2=0;ith2<ntheta;ith2++)
      for(int im2=0;im2<nmass;im2++)
	for(int im1=0;im1<nmass;im1++)
	  {
	    int ic=icombo_int(ith1,ith2,im1,im2);
	    if(ith1>=ith2) use[ic]=0;
	    else use[ic]=1;
	    if(ith1==ntheta-1 && ic<ntheta/2 ) use[ic]=0;
	    //load the corrs
	    c=load_corr(base_path,0,ith1,ith2,im1,im2,0,0,corr_name); //TM
	    aM.data[ic]=constant_fit(effective_mass(c.simmetrized(1)),tmin[0],tmax[0],combine("eff_plot_%02d_%02d_%02d_%02d.xmg",ith1,ith2,im1,im2).c_str());
	    
	    ofstream out(combine("plot_%02d_%02d_%02d_%02d.xmg",ith1,ith2,im1,im2).c_str());
	    out<<c<<endl;
	    cout<<"ith1:"<<ith1<<" ith2:"<<ith2<<" im1:"<<im1<<" im2:"<<im2<<" ic:"<<ic<<"/"<<ncombo<<" aM:"<<aM[ic]<<
	      " c[0]/c["<<TH<<"]="<<c[0]/c[TH]<<
	      endl;
	  }
  
  aM.write_to_binfile("M");
  
  for(int im2=0;im2<nmass;im2++)
    for(int im1=0;im1<nmass;im1++)
      {
	
	/////////////////////// check lattice free boson dispertion relation ///////////////////////
	{
	  ofstream out(combine("free_disp_rel_%02d_%02d.xmg",im1,im2).c_str());
	  out<<"@type xydy"<<endl;
	  out<<"@s0 line type 0"<<endl;      
	  out<<"@s0 symbol color 2"<<endl;
	  out<<"@s0 errorbar color 2"<<endl;
	  out<<"@s0 symbol 1"<<endl;
	  
	  double x[ntheta*ntheta];
	  jvec y(ntheta*ntheta,njack);
	  int ic=0;
	  for(int ith1=0;ith1<ntheta;ith1++)
	    for(int ith2=0;ith2<ntheta;ith2++)
	      {
		double qi=(theta[ith1]-theta[ith2])*M_PI/L;
		double q2=3*qi*qi;
		
		cout<<"th1:"<<theta[ith1]<<" th2:"<<theta[ith2]<<endl;
		
		x[ic]=4*3*sqr(sin(qi/2));
		y[ic]=4*sqr(sinh(aM[icombo_int(ith1,ith2,im1,im2)]/2));
		out<<x[ic]<<" "<<y[ic]<<endl;
		
		ic++;
	      }
	  
	  double ch2;
	  jack a(njack),b(njack),c(njack);
	  fixed_A=y[0];
	  parab_fit(ch2,a,b,c,x,y,"test_free_lattice_disp_rel.xmg",0);
	  cout<<"Parameters: "<<endl;
	  cout<<" "<<a<<" should be compatible with 0"<<endl;
	  cout<<" "<<b<<" should be compatible with 1"<<endl;
	  cout<<" "<<c<<" should be equal to "<<y[0]<<endl;
	  out.close();
	}
	
	////////////////////// check heavy meson on the lattice //////////////////////////////////
	{
	  ofstream out(combine("heavy_disp_rel_%02d_%02d.xmg",im1,im2).c_str());
	  out<<"@type xydy"<<endl;
	  out<<"@s0 line type 0"<<endl;      
	  out<<"@s0 symbol color 2"<<endl;
	  out<<"@s0 errorbar color 2"<<endl;
	  out<<"@s0 symbol 1"<<endl;
	  
	  double x[ntheta*ntheta];
	  jvec y(ntheta*ntheta,njack);
	  int ic=0;
	  for(int ith1=0;ith1<ntheta;ith1++)
	    for(int ith2=0;ith2<ntheta;ith2++)
	      {
		double qi=(theta[ith1]-theta[ith2])*M_PI/L;
		double q2=3*qi*qi;
		
		cout<<"th1:"<<theta[ith1]<<" th2:"<<theta[ith2]<<endl;
		
		x[ic]=q2;
		y[ic]=aM[icombo_int(ith1,ith2,im1,im2)];
		out<<x[ic]<<" "<<y[ic]<<endl;
		
		ic++;
	      }
	  
	  double ch2;
	  jack a(njack),b(njack),c(njack);
	  fixed_A=y[0];
	  parab_fit(ch2,a,b,c,x,y,"test_heavy_disp_rel.xmg",0);
	  cout<<"Parameters: "<<endl;
	  cout<<" "<<a<<" should be compatible with 0"<<endl;
	  cout<<" "<<0.5/b<<" should be compatible (in the continuum) with m (se below)"<<endl;
	  cout<<" "<<c<<" should be fixed to m="<<y[0]<<endl;
	  out.close();
	}
	
      }
  
  //check AKAK
  jvec csl=load_corr(base_path,0,0,0,0,0,0,0,"AKAK");
  jvec css=load_corr(base_path,1,0,0,0,0,0,0,"AKAK");
  
  {
    ofstream out("AKAK_comparison.xmg");
    out<<"@type xydy"<<endl;
    out<<"@s0 line type 0"<<endl;      
    out<<"@s0 symbol color 2"<<endl;
    out<<"@s0 errorbar color 2"<<endl;
    out<<"@s0 symbol 1"<<endl;
    out<<csl<<endl;
    out<<"&"<<endl;
    out<<"@type xydy"<<endl;
    out<<"@s1 line type 0"<<endl;      
    out<<"@s1 symbol color 3"<<endl;
    out<<"@s1 errorbar color 3"<<endl;
    out<<"@s1 symbol 2"<<endl;
    out<<css<<endl;
  }
  
  {
    ofstream out("AKAK_effmass_comparison.xmg");
    out<<"@type xydy"<<endl;
    out<<"@s0 line type 0"<<endl;      
    out<<"@s0 symbol color 2"<<endl;
    out<<"@s0 errorbar color 2"<<endl;
    out<<"@s0 symbol 1"<<endl;
    out<<effective_mass(csl.simmetrized(1))<<endl;
    out<<"&"<<endl;
    out<<"@type xydy"<<endl;
    out<<"@s1 line type 0"<<endl;      
    out<<"@s1 symbol color 3"<<endl;
    out<<"@s1 errorbar color 3"<<endl;
    out<<"@s1 symbol 2"<<endl;
    out<<effective_mass(css.simmetrized(1))<<endl;
  }
  
  return 0;
}
