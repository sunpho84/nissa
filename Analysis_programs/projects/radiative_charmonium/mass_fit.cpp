#include <include.h>

#define REAL 0
#define IMAG 1

int T,TH,L,tsep;
int njack,ntheta=2;

char base_path[1024];
int ibeta;
double lmass,cmass,theta;
char corr_name[1024];

int tminL,tmaxL;
int tminS,tmaxS;

int icombo_2pts(int r1,int r2,int ith1,int reim)
{
  return reim+2*(r1+2*(r2+2*ith1));
}

jvec load_2pts(const char *name,int ith1,int r1,int r2,int reim,const char *sl)
{
  return (
	  jvec_load(combine("%s/2pts_%s_%s",base_path,name,sl).c_str(),T,njack,icombo_2pts(r1,r2,ith1,reim))
	  +jvec_load(combine("%s/2pts_%s_%s",base_path,name,sl).c_str(),T,njack,icombo_2pts(!r1,!r2,ith1,reim))
	  )/2;
}

void read_data_list(const char *path)
{
  FILE *input_file=open_file(path,"r");
  
  read_formatted_from_file_expecting(base_path,input_file,"%s","base_path");
  read_formatted_from_file_expecting((char*)(&T),input_file,"%d","T");
  L=TH=T/2;
  read_formatted_from_file_expecting((char*)(&njack),input_file,"%d","njack");
  read_formatted_from_file_expecting((char*)(&ibeta),input_file,"%d","ibeta");
  read_formatted_from_file_expecting((char*)(&theta),input_file,"%lg","theta");
  read_formatted_from_file_expecting((char*)(&tsep),input_file,"%d","tsep");
  read_formatted_from_file_expecting((char*)(&lmass),input_file,"%lg","lmass");
  read_formatted_from_file_expecting((char*)(&cmass),input_file,"%lg","cmass");
}

void read_input()
{
  char data_list_path[1024];
  FILE *input_file=open_file("input","r");
  
  read_formatted_from_file_expecting(data_list_path,input_file,"%s","data_list_file");
  read_formatted_from_file_expecting(corr_name,input_file,"%s","corr_name");
  
  read_formatted_from_file_expecting((char*)(&tminL),input_file,"%d","tminL");
  read_formatted_from_file_expecting((char*)(&tmaxL),input_file,"%d","tmaxL");
  read_formatted_from_file_expecting((char*)(&tminS),input_file,"%d","tminS");
  read_formatted_from_file_expecting((char*)(&tmaxS),input_file,"%d","tmaxS");
  
  read_data_list(data_list_path);
    
  fclose(input_file);
}

int main()
{
  read_input();
  
  ///////////////////////////// Load two points for standing D and D* //////////////////////////
  
  //load ss
  jvec C_ss=load_2pts(corr_name,0, 0,0, 0, "30_30");
  //load sl
  jvec C_sl=load_2pts(corr_name,0, 0,0, 0, "30_00");
  
  
  //////////////////////////////////// Fit masses and Z for standing D and D* ////////////////////////////////////////
  
  //compute mass
  jack M,ZL,ZS;
  two_pts_SL_fit(M,ZL,ZS,C_sl.simmetrized(1),C_ss.simmetrized(1),tminL,tmaxL,tminS,tmaxS,"MSL.xmg","MSS.xmg");
  cout<<"mass: "<<M<<", Z: "<<ZL<<endl;
  
  M.write_to_binfile("M");
  
  return 0;
}
