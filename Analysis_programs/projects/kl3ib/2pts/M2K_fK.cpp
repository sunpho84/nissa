#include <include.h>
#include <iostream>
#include <math.h>

using namespace std;

//parameters
const int nboot=100;
const int njack=10;
const int nbeta=4;

int nel;
int *ibeta;
int ref_ml_beta[nbeta];
double *aml_bare;
boot *ml;

#include "../latpars.cpp"
#include "plot_funs.cpp"
const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (GeV)";
const char tag_a2[1024]="a\\S2\\N (fm)";

//physical input to fix dem
boot dmK_phys;
boot dm2K_phys;
//read data
bvec dam2K_fr_dl_bare,dfK_fr_afK_dl_bare;
bvec amK,afK;
//derived data
bvec dm2K_fr_dl,dmK_fr_dl,dl;
bvec dfK,dfK_fr_fK_dl,dfK_fr_fK;
bvec dfK_fr_dl;
bvec m2K,fK;
//output and plot quantities
bvec dl_chir;
boot dl_chir_cont;
//resulting parameters of the fits
bvec par_res_fit_dfK_fr_dl;
bvec par_res_fit_dfK_fr_fK;
bvec par_res_fit_dm2K_fr_dl;

double fun_fit_m2K(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{return A + B*ml + C*a*a + D*ml*ml + E*a*a*ml;}

double fun_fit_fK(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{
  double xill=db0*ml/(16*M_PI*M_PI*f0*f0);
  return A*(1-3.0/4*xill*log(xill)+B*xill+(C+E*xill)*a*a+D*xill*xill);
}

double fun_fit_dm2K_fr_dl(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{return B + E*a*a;}

double fun_fit_dfK_fr_dl(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{
  double dxill_dml=db0/(16*M_PI*M_PI*f0*f0);
  double xill=dxill_dml*ml;
  return A+B*xill*log(xill)+C*xill+D*a*a;
}

double fun_fit_dfK_fr_dfK(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{return fun_fit_dfK_fr_dl(A,B,C,D,E,f0,db0,ml,a);}

double fun_fit_dfK_fr_fK(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{
  double dxill_dml=db0/(16*M_PI*M_PI*f0*f0);
  double xill=dxill_dml*ml;
  return A+B*xill*log(xill)+C*xill+D*a*a;
}

double fun_fit_dl(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{return dmK_phys[plot_iboot]/(fun_fit_dm2K_fr_dl(A,B,C,D,E,f0,db0,ml,a)/2/pow(fun_fit_m2K(A,B,C,D,E,f0,db0,ml,a),0.5));}


double fun_fit_dfK(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{
  return fun_fit_dfK_fr_dl(A,B,C,D,E,f0,db0,ml,a)*fun_fit_dl(par_res_fit_dm2K_fr_dl[0][plot_iboot],par_res_fit_dm2K_fr_dl[1][plot_iboot],par_res_fit_dm2K_fr_dl[2][plot_iboot],par_res_fit_dm2K_fr_dl[3][plot_iboot],par_res_fit_dm2K_fr_dl[4][plot_iboot],f0,db0,ml,a);
}

void read_input(const char *path)
{
  ifstream input(path);
  if(input.good()!=1)
    {
      cerr<<"Erorr, unable to open file: "<<path<<endl;
      exit(1);
    }
  double njack; //dummy njack
  input>>nel>>njack;
  ibeta=(int*)malloc(sizeof(int)*nel);
  init_ijack_boot(nel,nboot,njack,11352637);
  
  aml_bare=(double*)malloc(sizeof(double)*nel);
  ml=(boot*)malloc(sizeof(boot)*nel);
  dam2K_fr_dl_bare=bvec(nel,nboot,njack);
  dm2K_fr_dl=bvec(nel,nboot,njack);
  dmK_fr_dl=bvec(nel,nboot,njack);
  dl=bvec(nel,nboot,njack);
  dfK=bvec(nel,nboot,njack);
  dfK_fr_fK=bvec(nel,nboot,njack);
  dfK_fr_fK_dl=bvec(nel,nboot,njack);
  dfK_fr_afK_dl_bare=bvec(nel,nboot,njack);
  dfK_fr_dl=bvec(nel,nboot,njack);
  amK=bvec(nel,nboot,njack);
  afK=bvec(nel,nboot,njack);
  m2K=bvec(nel,nboot,njack);
  fK=bvec(nel,nboot,njack);
  
  dmK_phys=boot(nboot,njack);
  dmK_phys.fill_gauss(dmK_phys_med,dmK_phys_err,2525245);
  dm2K_phys=dmK_phys*2*mK_phys;

  //reset index list of lighter ml
  for(int ib=0;ib<nbeta;ib++) ref_ml_beta[ib]=-1;
  
  for(int iel=0;iel<nel;iel++)
    {
      //read input
      char str[1024];
      input>>str>>ibeta[iel]>>aml_bare[iel];
      
      //read data
      jvec A(7,njack);
      A.load(str,0);
      
      //prepare boot
      ml[iel].create(nboot,njack);
      ml[iel]=aml_bare[iel]/lat[ibeta[iel]]/Zp[ibeta[iel]];
      boot_from_jack(dam2K_fr_dl_bare[iel],A[1],iel);
      boot_from_jack(dfK_fr_afK_dl_bare[iel],A[4],iel);
      boot_from_jack(amK[iel],A[5],iel);
      boot_from_jack(afK[iel],A[6],iel);
      
      //prepare renormalized quantities
      dm2K_fr_dl[iel]=dam2K_fr_dl_bare[iel]/lat[ibeta[iel]]*Zp[ibeta[iel]];
      dfK_fr_fK_dl[iel]=dfK_fr_afK_dl_bare[iel]*lat[ibeta[iel]]*Zp[ibeta[iel]]; //propagate Zp error
      dfK_fr_fK[iel]=dfK_fr_afK_dl_bare[iel]/dam2K_fr_dl_bare[iel]*sqr(lat[ibeta[iel]])*dm2K_phys;
      //calculate the physical quantity
      dmK_fr_dl[iel]=dm2K_fr_dl[iel]/(2*mK_phys);
      dl[iel]=dmK_phys/dmK_fr_dl[iel];
      m2K[iel]=sqr(amK[iel]/lat[ibeta[iel]]);
      fK[iel]=afK[iel]/lat[ibeta[iel]];
      dfK_fr_dl[iel]=dfK_fr_afK_dl_bare[iel]*afK[iel]*Zp[ibeta[iel]];
      dfK[iel]=dfK_fr_dl[iel]*dl[iel];
      
      cout<<ml[iel]<<" "<<dm2K_fr_dl[iel]<<" "<<fK[iel]<<endl;
      
      //set the lighter mass
      int b=ibeta[iel],r=ref_ml_beta[b];
      if(r==-1||ml[r].med()>ml[iel].med()) ref_ml_beta[b]=iel;
    }
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++) cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<ml[ref_ml_beta[ib]]<<" meV"<<endl;
  cout<<"---"<<endl;
}

//calculate the chi square
double chi2(double *p,int iboot,double (*fun)(double,double,double,double,double,double,double,double,double),double *X,double *Y,double *err_Y,double *apow)
{
  double ch2=0;

  for(int iel=0;iel<nel;iel++)
    {
      double f=fun(p[0],p[1],p[2],p[3],p[4],f0[iboot],db0[iboot],X[iel],lat[ibeta[iel]][iboot])*apow[ibeta[iel]];
      ch2+=sqr((Y[iel]-f)/err_Y[iel]);
    }
  
  return ch2;
}

//wrapper for the calculation of the chi2
double *fit_X;
double *fit_Y,*fit_err_Y;
double (*fit_fun)(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a);
int fit_iboot;
double fit_apow[4];
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{ch=chi2(p,fit_iboot,fit_fun,fit_X,fit_Y,fit_err_Y,fit_apow);}

void fit(boot &A,boot &B,boot &C,boot &D,boot &E,boot *X,bvec &Ya,double (*fun)(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a),bool *fix_par,double power)
{
  fit_fun=fun;
  bvec apow(4,nboot,njack);
  for(int ib=0;ib<4;ib++) apow[ib]=pow(lat[ib],power);
  bvec Y(nel,nboot,njack);
  for(int iel=0;iel<nel;iel++) Y[iel]=Ya[iel]*apow[ibeta[iel]];
  //allocate X, Y and err_Y and copy X
  fit_X=new double[nel];
  fit_Y=new double[nel];
  fit_err_Y=new double[nel];
  for(int iel=0;iel<nel;iel++)
    {
      fit_X[iel]=X[iel].med();
      fit_err_Y[iel]=Y[iel].err();
    }
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  
  minu.DefineParameter(0,"A",0.0,0.0001,0,0);
  minu.DefineParameter(1,"B",0.01,0.0001,0,0);
  minu.DefineParameter(2,"C",0,0.0001,0,0);
  minu.DefineParameter(3,"D",0,0.0001,0,0);
  minu.DefineParameter(4,"E",0,0.0001,0,0);

  int nfree=0;
  for(int ipar=0;ipar<5;ipar++)
    {
      if(fix_par[ipar]) minu.FixParameter(ipar);
      else nfree++;
    }

  minu.SetFCN(chi2_wr);
  
  double C2;
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      for(int ib=0;ib<4;ib++) fit_apow[ib]=apow[ib][iboot];
      fit_iboot=iboot;
      //copy current Y
      for(int iel=0;iel<nel;iel++)
	fit_Y[iel]=Y.data[iel].data[iboot];
	      
      //minimize
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,A.data[iboot],dum);
      minu.GetParameter(1,B.data[iboot],dum);
      minu.GetParameter(2,C.data[iboot],dum);
      minu.GetParameter(3,D.data[iboot],dum);
      minu.GetParameter(4,E.data[iboot],dum);
      
      if(iboot==0)
	{
	  double p[5]={A[iboot],B[iboot],C[iboot],D[iboot],E[iboot]};
	  C2=chi2(p,iboot,fit_fun,fit_X,fit_Y,fit_err_Y,fit_apow);
	}
    }
  
  //calculate the chi2
  cout<<"A = ("<<A<<"), B=("<<B<<"), C=("<<C<<"), D=("<<D<<"), E=("<<E<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<nel-nfree<<" = "<<C2/(nel-nfree)<<endl;
  
  delete[] fit_X;
  delete[] fit_Y;
  delete[] fit_err_Y;
}

void analysis_dm2K_fr_dl()
{
  const char tag_dl[1024]="\\xd\\0m\\sl\\N\\SMS,2GeV\\N (GeV)";
  const char tag_dm2K_fr_dl[1024]="\\xd\\0M\\S2\\N\\sK";
  
  //perform the fit
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack),D(nboot,njack),E(nboot,njack);
  bool fixing[5]={1,0,1,1,0};
  fit(A,B,C,D,E,ml,dm2K_fr_dl,fun_fit_dm2K_fr_dl,fixing,1);
  cout<<endl;  
  
  //chiral extrapolation
  boot dm2K_fr_dl_chir_cont(nboot,njack);
  bvec dm2K_fr_dl_estr_ml(nbeta,nboot,njack);
  for(int ib=0;ib<nbeta;ib++)
    for(int iboot=0;iboot <nboot+1;iboot++)
      {
	int r=ref_ml_beta[ib];
	dm2K_fr_dl_chir_cont.data[iboot]=fun_fit_dm2K_fr_dl(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0);
	dm2K_fr_dl_estr_ml.data[ib].data[iboot]=dm2K_fr_dl[r][iboot]*fun_fit_dm2K_fr_dl(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0)/fun_fit_dm2K_fr_dl(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml[r][iboot],0);
      }
  
  //calculate md-mu
  boot dmK_chir_cont=dm2K_fr_dl_chir_cont/(2*mK_phys);
  dl_chir_cont=dmK_phys/dmK_chir_cont;
  boot mu_chir_cont=ml_phys-dl_chir_cont/2;
  boot md_chir_cont=ml_phys+dl_chir_cont/2;
  boot mu_fr_md_chir_cont=mu_chir_cont/md_chir_cont;
  
  //chiral and continuum
  cout<<"dm2K = "<<dm2K_fr_dl_chir_cont<<endl;
  cout<<"dmK = "<<dmK_chir_cont<<endl;
  cout<<"md - mu = ("<<dl_chir_cont*1000<<") Mev"<<endl;
  cout<<"mu / md = ("<<mu_chir_cont*1000<<") / ("<<md_chir_cont*1000<<") = "<<mu_fr_md_chir_cont<<endl;
  
  par_res_fit_dm2K_fr_dl=bvec(5,nboot,njack);
  
  par_res_fit_dm2K_fr_dl.data[0]=A;
  par_res_fit_dm2K_fr_dl.data[1]=B;
  par_res_fit_dm2K_fr_dl.data[2]=C;
  par_res_fit_dm2K_fr_dl.data[3]=D;
  par_res_fit_dm2K_fr_dl.data[4]=E;
  
  plot_funz_ml("dm2K_fr_dl_funz_ml.xmg",tag_dm2K_fr_dl,tag_ml,tag_dm2K_fr_dl,ml,dm2K_fr_dl,par_res_fit_dm2K_fr_dl,ml_phys.med(),fun_fit_dm2K_fr_dl,dm2K_fr_dl_chir_cont);
  plot_funz_ml("dl_funz_ml.xmg",tag_dl,tag_ml,tag_dl,ml,dl,par_res_fit_dm2K_fr_dl,ml_phys.med(),fun_fit_dl,dl_chir_cont);
  plot_funz_a2("dm2K_fr_dl_funz_a2.xmg",tag_dm2K_fr_dl,tag_a2,tag_dm2K_fr_dl,lat_med_fm,dm2K_fr_dl_estr_ml,par_res_fit_dm2K_fr_dl,fun_fit_dm2K_fr_dl,dm2K_fr_dl_chir_cont);
}

void analysis_dfK_fr_dl()
{
  const char tag_dfK[1024]="\\xd\\0f\\sK (GeV)";
  const char tag_dfK_fr_dl[1024]="\\xd\\0f\\sK\\N/\\xd\\0m";
  
  //perform the fit
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack),D(nboot,njack),E(nboot,njack);
  bool fixing[5]={0,0,0,0,1};
  fit(A,B,C,D,E,ml,dfK_fr_dl,fun_fit_dfK_fr_dl,fixing,0);
  cout<<endl;  
  
  //chiral extrapolation
  boot dfK_fr_dl_chir_cont(nboot,njack);
  bvec dfK_fr_dl_estr_ml(nbeta,nboot,njack);
  for(int ib=0;ib<nbeta;ib++)
    for(int iboot=0;iboot <nboot+1;iboot++)
      {
	int r=ref_ml_beta[ib];
	dfK_fr_dl_chir_cont.data[iboot]=fun_fit_dfK_fr_dl(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0);
	dfK_fr_dl_estr_ml.data[ib].data[iboot]=dfK_fr_dl[r][iboot]*fun_fit_dfK_fr_dl(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0)/fun_fit_dfK_fr_dl(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml[r][iboot],0);
      }
  
  //chiral and continuum
  boot dfK_chir_cont=dfK_fr_dl_chir_cont*dl_chir_cont;
  cout<<"dfK = "<<dfK_fr_dl_chir_cont<<endl;
  cout<<"fK+ - fK0 = ("<<dfK_chir_cont*1000<<") MeV"<<endl;
  cout<<"fK+ / fK - 1 = ("<<dfK_chir_cont/0.158/2<<")"<<endl;
  
  par_res_fit_dfK_fr_dl=bvec(5,nboot,njack);
  par_res_fit_dfK_fr_dl.data[0]=A;
  par_res_fit_dfK_fr_dl.data[1]=B;
  par_res_fit_dfK_fr_dl.data[2]=C;
  par_res_fit_dfK_fr_dl.data[3]=D;
  par_res_fit_dfK_fr_dl.data[4]=E;
  
  plot_funz_ml("dfK_fr_dl_funz_ml.xmg",tag_dfK_fr_dl,tag_ml,tag_dfK_fr_dl,ml,dfK_fr_dl,par_res_fit_dfK_fr_dl,ml_phys.med(),fun_fit_dfK_fr_dl,dfK_fr_dl_chir_cont);
  plot_funz_ml("dfK_funz_ml.xmg",tag_dfK,tag_ml,tag_dfK,ml,dfK,par_res_fit_dfK_fr_dl,ml_phys.med(),fun_fit_dfK,dfK_chir_cont);
  plot_funz_a2("dfK_fr_dl_funz_a2.xmg",tag_dfK_fr_dl,tag_a2,tag_dfK_fr_dl,lat_med_fm,dfK_fr_dl_estr_ml,par_res_fit_dfK_fr_dl,fun_fit_dfK_fr_dl,dfK_fr_dl_chir_cont);
}

void analysis_dfK_fr_fK()
{
  const char tag_dfK_fr_fK[1024]="\\xd\\0f\\sK\\N/\\xd\\0m";

  //perform the fit
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack),D(nboot,njack),E(nboot,njack);
  bool fixing[5]={0,0,0,0,1};
  fit(A,B,C,D,E,ml,dfK_fr_fK,fun_fit_dfK_fr_fK,fixing,0);
  cout<<endl;  
  
  //chiral extrapolation
  boot dfK_fr_fK_chir_cont(nboot,njack);
  bvec dfK_fr_fK_estr_ml(nbeta,nboot,njack);
  for(int ib=0;ib<nbeta;ib++)
    for(int iboot=0;iboot <nboot+1;iboot++)
      {
	int r=ref_ml_beta[ib];
	dfK_fr_fK_chir_cont.data[iboot]=fun_fit_dfK_fr_fK(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0);
	dfK_fr_fK_estr_ml.data[ib].data[iboot]=dfK_fr_fK[r][iboot]*fun_fit_dfK_fr_fK(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0)/fun_fit_dfK_fr_dfK(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml[r][iboot],0);
      }
  
  //chiral and continuum
  cout<<"(fK+ - fK0) / fK = ("<<dfK_fr_fK_chir_cont<<")"<<endl;
  cout<<"fK+/fK - 1 = ("<<0.5*dfK_fr_fK_chir_cont<<")"<<endl;

  par_res_fit_dfK_fr_fK=bvec(5,nboot,njack);
  par_res_fit_dfK_fr_fK.data[0]=A;
  par_res_fit_dfK_fr_fK.data[1]=B;
  par_res_fit_dfK_fr_fK.data[2]=C;
  par_res_fit_dfK_fr_fK.data[3]=D;
  par_res_fit_dfK_fr_fK.data[4]=E;
  
  plot_funz_ml("dfK_fr_fK_funz_ml.xmg",tag_dfK_fr_fK,tag_ml,tag_dfK_fr_fK,ml,dfK_fr_fK,par_res_fit_dfK_fr_fK,ml_phys.med(),fun_fit_dfK_fr_fK,dfK_fr_fK_chir_cont);
  plot_funz_a2("dfK_fr_fK_funz_a2.xmg",tag_dfK_fr_fK,tag_a2,tag_dfK_fr_fK,lat_med_fm,dfK_fr_fK_estr_ml,par_res_fit_dfK_fr_fK,fun_fit_dfK_fr_fK,dfK_fr_fK_chir_cont);
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      cerr<<"Error, use: "<<arg[0]<<" input"<<endl;
      exit(1);
    }
  
  init_latpars();
  read_input(arg[1]);
  
  analysis_dm2K_fr_dl();
  cout<<"---"<<endl;
  analysis_dfK_fr_dl();
  cout<<"---"<<endl;
  analysis_dfK_fr_fK();
  
  return 0;
}
