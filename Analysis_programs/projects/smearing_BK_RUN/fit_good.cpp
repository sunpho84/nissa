#include <include.h>
#include <iostream>
#include <sstream>

using namespace std;

int nmass,nlights,iml_un;
int njack=16;
double *mass;
jvec c[2];
int T,TH;
char base_path[1024];
int tSEP;
int ibeta;
int nwalls;
double afm[4]={0.0977,0.0847,0.0671,0.0536};
int ijack_fit;

char corr_name[1024],out_file[1024];
int tmin[3],tmax[3];
int ifit_int,parity;

double fun_fit(double Z1,double Z2,double M,int t)
{return Z1*Z2*exp(-M*TH)*cosh(M*(TH-t))/M;}

void ch2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M=p[0];
  double ZL=p[1];
  double ZS=p[2];
  
  for(int t=tmin[ifit_int];t<=min(tmax[ifit_int],TH);t++)
    {
      ch+=sqr((c[0][t].data[ijack_fit]-fun_fit(ZL,ZS,M,t))/c[0][t].err());
      ch+=sqr((c[1][t].data[ijack_fit]-fun_fit(ZS,ZS,M,t))/c[1][t].err());
    }
}

int icombo(int im1,int im2,int r1,int r2,int iwall,int ri)
{
  if(im1>=nlights)
    {
      cerr<<"Error, requiring uncomputed combo!"<<endl;
      exit(1);
    }
  
  int ic=ri+2*(iwall+nwalls*(r1+2*(im1+nlights*(r2+2*im2))));

  return ic;
}

jvec load_corr(char *base_path,int sm_lev,int im1,int im2,int r1,int r2,const char *obs_name)
{
  char path[1024];
  int tag[2]={0,30};
  sprintf(path,"%s/%s_30_%02d",base_path,obs_name,tag[sm_lev]);

  jvec l=(jvec_load(path,T,njack,icombo(im1,im2,!r1,!r2,0,0))+
	  jvec_load(path,T,njack,icombo(im1,im2, r1, r2,0,0)))/2;
  
  jvec r=(jvec_load(path,T,njack,icombo(im1,im2,!r1,!r2,1,0))+
	  jvec_load(path,T,njack,icombo(im1,im2, r1, r2,1,0)))/2;
  

  for(int tL=0;tL<T;tL++)
    {
      int tR=tL+tSEP;
      if(tR>=T) tR-=T;
      
      l.data[tL]+=r.data[tR];
    }
  
  return l/2;
}

void read_ensemble_pars(char *base_path,int &T,int &ibeta,int &nmass,double *&mass,int &nlights,const char *data_list_file)
{
  FILE *input=open_file(data_list_file,"r");
  
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  nwalls=2;
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");

  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");
  read_formatted_from_file_expecting((char*)&iml_un,input,"%d","iml_un");
  read_formatted_from_file_expecting((char*)&nlights,input,"%d","nlights");
  read_formatted_from_file_expecting((char*)&tSEP,input,"%d","tSEP");
  
  fclose(input);
}

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  char data_list_file[256];
  read_formatted_from_file_expecting(data_list_file,fin,"%s","data_list_file");

  read_formatted_from_file_expecting(corr_name,fin,"%s","corr_name");
  read_formatted_from_file_expecting((char*)(&parity),fin,"%d","parity");
  
  read_formatted_from_file_expecting((char*)(&tmin[0]),fin,"%d","tint_ll");
  read_formatted_from_file((char*)(&tmax[0]),fin,"%d","tint_ll");
  
  read_formatted_from_file_expecting((char*)(&tmin[1]),fin,"%d","tint_lh");
  read_formatted_from_file((char*)(&tmax[1]),fin,"%d","tint_lh");
  
  read_formatted_from_file_expecting((char*)(&tmin[2]),fin,"%d","tint_hh");
  read_formatted_from_file((char*)(&tmax[2]),fin,"%d","tint_hh");
  
  read_formatted_from_file_expecting(out_file,fin,"%s","out_file");
  
  fclose(fin);
  
  read_ensemble_pars(base_path,T,ibeta,nmass,mass,nlights,data_list_file);
  
  TH=T/2;
}

string rectangle(int a,int b,jack c)
{
  ostringstream out;
  
  out<<a<<" "<<c.med()-c.err()<<endl<<b<<" "<<c.med()-c.err()<<endl;
  out<<b<<" "<<c.med()+c.err()<<endl<<a<<" "<<c.med()+c.err()<<endl;
  out<<a<<" "<<c.med()-c.err()<<endl;	    
  
  return out.str();
}

int main(int narg,char **arg)
{
  read_pars("input");
  
  double a=afm[ibeta]/0.197;
  jack mass_estim[2];
  jack Z_estim[2];
  
  //loop over heavier mass
  int ncombo=nlights*nmass;
  jvec Zl(ncombo,njack),Zs(ncombo,njack),aM(ncombo,njack);
  TMinuit minu(3);
  minu.SetFCN(ch2);
  minu.SetPrintLevel(-1);
  //WARNING!!!! here order reverted!!!!
  for(int im1=0;im1<nlights;im1++)
    for(int im2=0;im2<nmass;im2++)
      {
	int ic=im1*nmass+im2;
	
	if(im2>=nlights) ifit_int=1;
	else ifit_int=0;
	
	jvec eff[2];
	for(int sm_lev=0;sm_lev<2;sm_lev++)
	  {
	    //load the current combo
	    c[sm_lev]=load_corr(base_path,sm_lev,im1,im2,0,0,corr_name).simmetrized(1); //TM
	    if(corr_name==string("S0S0")) c[sm_lev]=-c[sm_lev];

	    //compute effective mass
	    eff[sm_lev]=effective_mass(c[sm_lev]);
	  }

	for(int sm_lev=0;sm_lev<2;sm_lev++)
	  {	    
	    //check intervals
	    int tin=tmin[ifit_int];
	    int tfi=tin;
	    bool adv=1;
	    while(tfi<tmax[ifit_int] && adv)
	      {
		for(int ijack=0;ijack<=njack;ijack++)
		  adv=adv && (eff[0][tfi][ijack]>0 && eff[1][tfi][ijack]>0);
		if(adv) tfi++;
	      }	
	    if(tfi!=tmax[ifit_int]) cout<<"Warning, combo "<<ic<<" restricting to "<<tfi<<" instead than "<<tmax[ifit_int]<<endl;
	    
	    //fit mass
	    mass_estim[sm_lev]=constant_fit(eff[sm_lev],tmin[ifit_int],tmax[ifit_int]);
	    
	    //compute the corr func removed of the temporal dependency,
	    //fit Z and print it
	    {
	      ofstream out(combine("rem_tdep_%02d_%02d_%d.xmg",im1,im2,sm_lev).c_str());
	      out<<"@type xydy"<<endl<<"@s0 line type 0"<<endl;
	      jvec temp=c[sm_lev];
	      for(int iel=0;iel<TH;iel++)
		for(int ijack=0;ijack<njack+1;ijack++)
		  temp[iel].data[ijack]=temp[iel].data[ijack]/fun_fit(1,1,mass_estim[sm_lev][ijack],iel);
	      Z_estim[sm_lev]=constant_fit(temp,tmin[ifit_int],tmax[ifit_int]);
	      out<<temp<<"&\n@type xy\n";
	      out<<rectangle(tmin[ifit_int],tmax[ifit_int],Z_estim[sm_lev])<<endl;
	    }
	  }
	
	Z_estim[1]=sqrt(Z_estim[1]);
	Z_estim[0]/=Z_estim[1];
	
	//compute Zl and M
	aM[ic]=(mass_estim[0]+mass_estim[1])/2;
	Zl[ic]=Z_estim[0];
	
	//print M and f
	jack M=aM[ic]/a;
	jack af=Zl[ic]*(mass[im1]+mass[im2])/aM[ic]/sinh(aM[ic]);
	jack f=af/a;
	
	for(ijack_fit=0;ijack_fit<=njack;ijack_fit++)
	  {
	    minu.DefineParameter(0,"aM",aM[ic][ijack_fit],aM[ic].err(),0,0);
	    minu.DefineParameter(1,"ZL",Z_estim[0][ijack_fit],Z_estim[0].err(),0,0);
	    minu.DefineParameter(2,"ZS",Z_estim[1][ijack_fit],Z_estim[1].err(),0,0);
	    minu.Migrad();
	    double dum;
	    minu.GetParameter(0,aM[ic].data[ijack_fit],dum);
	    minu.GetParameter(1,Zl[ic].data[ijack_fit],dum);
	    minu.GetParameter(2,Zs[ic].data[ijack_fit],dum);
	  }
	af=Zl[ic]*(mass[im1]+mass[im2])/aM[ic]/sinh(aM[ic]);
	f=af/a;
	cout<<"im1:"<<im1<<" im2:"<<im2<<" ic:"<<ic<<"/"<<ncombo<<" M: "<<aM[ic]<<", f:"<<f<<" Zl:"<<Zl[ic]<<" Zs:"<<Zs[ic]<<" Zest:"<<Z_estim[0]<<", "<<Z_estim[1]<<endl;
	for(int sm_lev=0;sm_lev<2;sm_lev++)
	  {
	    ofstream fout(combine("out_%02d_%02d_%d.xmg",im1,im2,sm_lev).c_str());
	    fout<<"@type xydy"<<endl<<c[sm_lev]<<endl<<"&\n@type xy\n";
	    for(int t=tmin[ifit_int];t<=tmax[ifit_int];t++)
	      if(sm_lev==0) fout<<t<<" "<<fun_fit(Zl[ic][njack],Zs[ic][njack],aM[ic][njack],t)<<endl;
	      else          fout<<t<<" "<<fun_fit(Zs[ic][njack],Zs[ic][njack],aM[ic][njack],t)<<endl;
	    fout.close();
	  }
	
	//fit mass and print it
	for(int sm_lev=0;sm_lev<2;sm_lev++)
	  {
	    ofstream out(combine("eff_%02d_%02d_%d.xmg",im1,im2,sm_lev).c_str());
	    out<<"@type xydy"<<endl<<"@s0 line type 0"<<endl;
	    out<<eff[sm_lev];
	    out<<"&\n@type xy\n";
	    out<<rectangle(tmin[ifit_int],tmax[ifit_int],aM[ic])<<endl;
	  }
	    

      }
  
  aM.write_to_binfile(out_file);
  (Zl*Zl).append_to_binfile(out_file);

  return 0;
}
