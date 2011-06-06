#include <include.h>
#include <iostream>

using namespace std;

int nmass,nm_low;
int nsm_lev,njack;
int tmin_list[3][2];
double *mass;
jvec **c;
int T,TH;
char base_path[1024];

double fun_fit(double Z1,double Z2,double M,int t)
{
  return Z1*Z2*exp(-M*TH)*cosh(M*(TH-t))/M;
}

int icombo(int im1,int im2,int r1,int r2)
{
  int imin=min(im1,im2);
  int imax=max(im1,im2);
  
  return r1+2*imax+nmass*2*(imin*2+r2);
}

jvec load_corr(int sm_lev_so,int sm_lev_si,int im1,int im2)
{
  char path[1024];
  sprintf(path,"%save_%d%d",base_path,sm_lev_so,sm_lev_si);
  
  return ((jvec_load(path,T,njack,icombo(im1,im2,0,0))+jvec_load(path,T,njack,icombo(im1,im2,1,1)))/2).simmetrized(1);
}

int ijack_fit=0;
void ch2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double ZL0=p[0];
  double M0=p[1];
  
  double *Z=p+2;
  double M=p[2+nsm_lev];
  
  for(int t=0;t<=TH;t++)
    {
      if(t>=16) ch+=sqr((c[0][0][t].data[ijack_fit]-fun_fit(ZL0,ZL0,M0,t))/c[0][0][t].err());
      
      for(int sm_lev_so=0;sm_lev_so<nsm_lev;sm_lev_so++)
	for(int sm_lev_si=0;sm_lev_si<nsm_lev;sm_lev_si++)
	  if(t>=6) ch+=sqr((c[sm_lev_so][sm_lev_si][t].data[ijack_fit]-fun_fit(Z[sm_lev_so],Z[sm_lev_si],M,t))/c[sm_lev_so][sm_lev_si][t].err());
    }
}

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  
  read_formatted_from_file_expecting(base_path,fin,"%s","base_path");
  read_formatted_from_file_expecting((char*)&nmass,fin,"%d","nmass");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)(&mass[imass]),fin,"%lg","mass");
  read_formatted_from_file_expecting((char*)&nm_low,fin,"%d","nmass_low");
  read_formatted_from_file_expecting((char*)&njack,fin,"%d","njack");
  read_formatted_from_file_expecting((char*)&nsm_lev,fin,"%d","nsm_lev");
  read_formatted_from_file_expecting((char*)&T,fin,"%d","T");
  TH=T/2;
  
  read_formatted_from_file_expecting((char*)(&tmin_list[0][0]),fin,"%d","tmin_ll");
  read_formatted_from_file((char*)(&tmin_list[0][1]),fin,"%d","tmin_ll_smeared");

  read_formatted_from_file_expecting((char*)(&tmin_list[1][0]),fin,"%d","tmin_lh");
  read_formatted_from_file((char*)(&tmin_list[1][1]),fin,"%d","tmin_lh_smeared");
  
  read_formatted_from_file_expecting((char*)(&tmin_list[2][0]),fin,"%d","tmin_hh");
  read_formatted_from_file((char*)(&tmin_list[2][1]),fin,"%d","tmin_hh_smeared");
}

int main()
{
  read_pars("input");
  
  jack mass_estim;

  //allocate correlation function
  c=new jvec*[nsm_lev];
  for(int sm_lev_so=0;sm_lev_so<nsm_lev;sm_lev_so++) c[sm_lev_so]=new jvec[nsm_lev];

  //loop over heavier mass
  int icombo=0,ncombo=nmass*nmass;
  jvec Z0(ncombo,njack),M0(ncombo,njack);
  jvec Z(ncombo,njack),M(ncombo,njack);
  for(int im1=0;im1<nmass;im1++)
    for(int im2=0;im2<nmass;im2++)
      {
	int combo_type=(int)(im1>=nm_low)+(int)(im2>=nm_low);
	
	for(int sm_lev_so=0;sm_lev_so<nsm_lev;sm_lev_so++)
	  for(int sm_lev_si=0;sm_lev_si<nsm_lev;sm_lev_si++)
	    {
	      int smear_type=(sm_lev_so or sm_lev_si);
	      int tmin=tmin_list[combo_type][smear_type];
	      if(sm_lev_so>1||sm_lev_si>1) tmin=13;
	      
	      //load the current combo
	      c[sm_lev_so][sm_lev_si]=load_corr(sm_lev_so,sm_lev_si,im1,im2);
	      c[sm_lev_so][sm_lev_si].print_to_file(combine("out_%02d_%02d_%d%d.xmg",im1,im2,sm_lev_so,sm_lev_si).c_str());
	      
	      //print effective mass
	      jvec eff=effective_mass(c[sm_lev_so][sm_lev_si]);
	      
	      //fit mass
	      mass_estim=constant_fit(eff,tmin,TH);
	      
	      //prepare the plot of the mass fit
	      ofstream out(combine("plot_%02d_%02d_%d%d.xmg",im1,im2,sm_lev_so,sm_lev_si).c_str());
	      out<<"@type xydy"<<endl<<"@s0 line type 0"<<endl;
	      out<<eff<<endl;
	      out<<"&"<<endl<<"@type xy"<<endl;
	      double av_mass=mass_estim.med(),er_mass=mass_estim.err();
	      out<<tmin<<" "<<av_mass-er_mass<<endl<<TH<<" "<<av_mass-er_mass<<endl;
	      out<<TH<<" "<<av_mass+er_mass<<endl<<tmin<<" "<<av_mass+er_mass<<endl;
	      out<<tmin<<" "<<av_mass-er_mass<<endl;	    
	    }
      
	cout<<"Fitting combo "<<im1<<" "<<im2<<" mass estimate: "<<mass_estim<<endl;
	cout<<"---------------------------------"<<endl;
      
	TMinuit minu(1+nsm_lev);
	minu.SetFCN(ch2);
	minu.SetPrintLevel(-1);
      
	int t_LOC=tmin_list[combo_type][0];
	int t_SME=tmin_list[combo_type][1];
      
	for(ijack_fit=0;ijack_fit<njack+1;ijack_fit++)
	  {
	    jack Z_estim=sqrt(c[0][0][t_LOC]/fun_fit(1,1,mass_estim[ijack_fit],t_LOC));
	    minu.DefineParameter(0,"ZL0",Z_estim[ijack_fit],0.001,0,100);
	    minu.DefineParameter(1,"M0",mass_estim[ijack_fit],0.001,0,0);
	    for(int sm_lev=0;sm_lev<nsm_lev;sm_lev++)
	      {
		jack Z_estim=sqrt(c[sm_lev][sm_lev][t_SME]/fun_fit(1,1,mass_estim[ijack_fit],t_SME));
		minu.DefineParameter(2+sm_lev,"Z",Z_estim[ijack_fit],0.001,0,100);
	      }
	    minu.DefineParameter(2+nsm_lev,"M",mass_estim[ijack_fit],0.001,0,0);
	  
	    double dum;
	  
	    //minu.mnseek();
	    //minu.mnsimp();
	    minu.Migrad();
	    
	    minu.GetParameter(0,Z0[icombo].data[ijack_fit],dum);
	    minu.GetParameter(1,M0[icombo].data[ijack_fit],dum);
	    minu.GetParameter(2,Z[icombo].data[ijack_fit],dum);
	    minu.GetParameter(2+nsm_lev,M[icombo].data[ijack_fit],dum);
	  }
      
	cout<<Z0[icombo]<<endl;
	cout<<M0[icombo]<<endl;
	cout<<Z[icombo]<<endl;
	cout<<M[icombo]<<endl;
	
	icombo++;
      }
  
  ofstream Mout("M.xmg");
  icombo=0;
  Mout<<"@type xydy"<<endl;
  for(int im1=0;im1<nmass;im1++)
    for(int im2=0;im2<nmass;im2++)
      Mout<<mass[im1]<<" "<<mass[im2]<<" "<<M0[icombo++].err()<<endl;
  Mout<<"&"<<endl;
  icombo=0;
  for(int im1=0;im1<nmass;im1++)
    for(int im2=0;im2<nmass;im2++)
      Mout<<mass[im1]<<" "<<mass[im2]<<" "<<M[icombo++].err()<<endl;
  
  return 0;
}
