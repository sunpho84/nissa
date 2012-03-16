#include "../HH_common.cpp"

const int nbeta=4;

//calculate the chi square
double chi2(double A,double B,double *a)
{
  double ch2=0;
  
  for(int iens=0;iens<nens;iens++)
    if(include_380 || ibeta[iens]!=0)
      ch2+=pow((Y_fit[iens]-fun_fit_Mh(A,B,C,D,X_fit[iens],a[ibeta[iens]]))/err_Y_fit[iens],2);
  
  return ch2;
}

//wrapper for the calculation of the chi2
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0],B=p[1],C=p[2],D=p[3];
  double *a=p+4;
  ch=chi2(A,B,C,D,a);
}

void fit(boot &A,boot &B,boot &C,boot &D,bvec &X,bvec &Y)
{
  //copy X
  X_fit=new double[nens];
  for(int iens=0;iens<nens;iens++) X_fit[iens]=X[iens].med();
  Y_fit=new double[nens];
  err_Y_fit=new double[nens];
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  
  int npars=4;
  minu.DefineParameter(0,"A",0.0,0.0001,0,0);
  minu.DefineParameter(1,"B",0.0,0.0001,0,0);
  minu.DefineParameter(2,"C",0.0,0.0001,0,0);
  minu.DefineParameter(3,"D",0.0,0.0001,0,0);
  if(!include_a4)
    {
      minu.FixParameter(3);
      npars--;
    }
  minu.SetFCN(chi2_wr);
  
  double C2;
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      if(iboot>0)
        minu.SetPrintLevel(-1);
      
      minu.DefineParameter(4,"a380",lat[0][iboot],0.0001,0,0);
      minu.DefineParameter(5,"a390",lat[1][iboot],0.0001,0,0);
      minu.DefineParameter(6,"a405",lat[2][iboot],0.0001,0,0);
      minu.DefineParameter(7,"a420",lat[3][iboot],0.0001,0,0);
      minu.FixParameter(4);
      minu.FixParameter(5);
      minu.FixParameter(6);
      minu.FixParameter(7);
      
      for(int iens=0;iens<nens;iens++)
        {
          Y_fit[iens]=Y.data[iens].data[iboot];
          err_Y_fit[iens]=Y.data[iens].err();
        }
      
      //minimize
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,A.data[iboot],dum);
      minu.GetParameter(1,B.data[iboot],dum);
      minu.GetParameter(2,C.data[iboot],dum);
      minu.GetParameter(3,D.data[iboot],dum);
      
      double lat_med[4]={lat[0].med(),lat[1].med(),lat[2].med(),lat[3].med()};
      if(iboot==0) C2=chi2(A.data[iboot],B[iboot],C[iboot],D[iboot],lat_med);
    }
  
  int ninc_ens=0;
  for(int iens=0;iens<nens;iens++)
    if(ibeta[iens]!=0 || include_380) ninc_ens++;
  
  //calculate the chi2
  cout<<"A=("<<A<<"), B=("<<B<<"), C=("<<C<<"), D=("<<D<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<ninc_ens-npars<<" = "<<C2/(ninc_ens-npars)<<endl;
  
  delete[] X_fit;
  delete[] Y_fit;
  delete[] err_Y_fit;
}

double aM0[4]={1.264,1.131,0.942,0.780};
double aM1[4]={1.71,1.42,1.08,0.84};

double eM0=0.001;
double eM1[4]={0.01,0.02,0.01,0.01};

int main(int narg,char **arg)
{
  init_latpars();
  
  bvec M0(4,nboot,njack);
  bvec M1(4,nboot,njack);
  
  for(int ibeta=0;ibeta<nbeta;ibeta++)
    {
      M0[ibeta].fill_gauss(aM0[ibeta],eM0,67832854+ibeta);
      M1[ibeta].fill_gauss(aM1[ibeta],eM1[ibeta],832854+ibeta);
      
      M0[ibeta]/=lat[ibeta];
      M1[ibeta]/=lat[ibeta];
    }
  
  
  ofstream out("/tmp/check_fermilab.xmg");
  out<<"@type xydy"<<endl;
  
  for(int ibeta=0;ibeta<nbeta;ibeta++)
    out<<(lat[ibeta]*lat[ibeta]).med()<<" "<<M0[ibeta]<<endl;
  
  
  out<<"&"<<endl;
  
  for(int ibeta=0;ibeta<nbeta;ibeta++)
    out<<(lat[ibeta]*lat[ibeta]).med()<<" "<<M1[ibeta]<<endl;
  
  out.close();

  //Mh_chir_cont.write_to_binfile("results");
  
  return 0;
}
