#include <include.h>
#include <iostream>
#include <sstream>

using namespace std;

int nmass,nlights,iml_un;
int njack=16;
double *mass;
int T,TH;
char base_path[1024];
int tSEP;
int ibeta;
int nwalls;
int ijack_fit;

char out_file[1024];
int tmin[3],tmax[3];
int ifit_int;

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
  jvec l=(jvec_load(path,T,njack,icombo(im1,im2,r1,r2,0,0))+
	  jvec_load(path,T,njack,icombo(im1,im2,r1,r2,0,0)))/2;
  
  jvec r=(jvec_load(path,T,njack,icombo(im1,im2,r1,r2,1,0))+
	  jvec_load(path,T,njack,icombo(im1,im2,r1,r2,1,0)))/2;
  
  
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
  
  //loop over heavier mass
  int ncombo=nlights*nmass;
  jvec Z_TfrV(ncombo,njack);

  jvec c[ncombo];
  
  //WARNING!!!! here order reverted!!!!
  for(int im1=0;im1<nlights;im1++)
    for(int im2=0;im2<nmass;im2++)
      {
	int ic=im1*nmass+im2;
	
	if(im2>=nlights) ifit_int=1;
	else ifit_int=0;

	jvec ct[2];
	//load the current combo
	for(int sm_lev=0;sm_lev<2;sm_lev++)
	  ct[sm_lev]=(-load_corr(base_path,sm_lev,im1,im2,0,0,"TKTK")/load_corr(base_path,sm_lev,im1,im2,0,0,"VKVK")).simmetrized(1); //TM
	
	c[ic]=ct[0]/sqrt(ct[1]);
	//compute Z_TfrV_local
	Z_TfrV[ic]=constant_fit(c[ic],tmin[ifit_int],tmax[ifit_int]);
	
	cout<<"im1:"<<im1<<" im2:"<<im2<<" ic:"<<ic<<"/"<<ncombo<<" Z:"<<Z_TfrV[ic]<<endl;
	
	ofstream fout(combine("out_%02d_%02d.xmg",im1,im2).c_str());
	fout<<"@type xydy"<<endl<<c[ic]<<endl<<"&\n@type xy\n";
	fout<<rectangle(tmin[ifit_int],tmax[ifit_int],Z_TfrV[ic]);
	fout.close();		    
      }  
  
  Z_TfrV.write_to_binfile(out_file);

  return 0;
}
