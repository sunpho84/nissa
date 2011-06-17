#include <include.h>
#include <iostream>
#include <math.h>

using namespace std;

int nel;
int *ibeta;
const int nboot=100;
int njack=10;
const int nbeta=4;

#include "latpars.cpp"

int **ijack_boot;
double *aml_bare;
boot *ml;
int ref_ml_beta[nbeta];

//derived data
bvec delta_fK_fr_fK;
boot delta_m2K_phys;
bvec par_res_fit;

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

int plot_iboot;

double fun_fit(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{
  double dxill_dml=db0/(16*M_PI*M_PI*f0*f0);
  double xill=dxill_dml*ml;
  
  return A+B*xill*log(xill)+C*xill+D*a*a;
  //return A+B*xill*(xill)+C*xill+D*a*a;
}

int rni(int n)
{return rand()%n;}

void boot_from_jack(boot &out,jack &in,int iel)
{
  int nboot=out.nboot;
  for(int iboot=0;iboot<nboot+1;iboot++) out.data[iboot]=in.data[ijack_boot[iel][iboot]];
}

void read_input(const char *path)
{
  ifstream input(path);
  if(input.good()!=1)
    {
      cerr<<"Erorr, unable to open file: "<<path<<endl;
      exit(1);
    }
  input>>nel>>njack;
  ibeta=(int*)malloc(sizeof(int)*nel);
  
  aml_bare=(double*)malloc(sizeof(double)*nel);
  ml=(boot*)malloc(sizeof(boot)*nel);
  ijack_boot=(int**)malloc(sizeof(int)*nel);
  for(int iel=0;iel<nel;iel++)
    {
      ijack_boot[iel]=(int*)malloc(sizeof(int)*(nboot+1));
      for(int iboot=0;iboot<nboot;iboot++) ijack_boot[iel][iboot]=rni(njack);
      ijack_boot[iel][nboot]=njack;
    }	
  
  boot delta_mK_phys=boot(nboot,njack);
  delta_mK_phys.fill_gauss(delta_mK_phys_med,delta_mK_phys_err,2525245);
  delta_m2K_phys=delta_mK_phys*2*mK_phys;
  cout<<"delta_MK_phys="<<delta_mK_phys<<endl;
  cout<<"delta_M2K_phys="<<delta_m2K_phys<<endl;
  delta_fK_fr_fK=bvec(nel,nboot,njack);
  
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
      
      boot delta_am2K_fr_delta_ml_bare(nboot,njack);
      boot delta_fK_fr_afK_delta_ml_bare(nboot,njack);
      boot afK(nboot,njack);
      boot_from_jack(delta_fK_fr_afK_delta_ml_bare,A[4],iel);
      boot_from_jack(delta_am2K_fr_delta_ml_bare,A[1],iel);
      boot_from_jack(afK,A[6],iel);
      
      //prepare renormalized quantities
      delta_fK_fr_fK[iel]=delta_fK_fr_afK_delta_ml_bare/delta_am2K_fr_delta_ml_bare*sqr(lat[ibeta[iel]]);
      //jack temp=A[3]/A[0]*sqr(A[5]);
      //set the lighter mass
      int b=ibeta[iel],r=ref_ml_beta[b];
      if(r==-1||ml[r].med()>ml[iel].med()) ref_ml_beta[b]=iel;
    }
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++) cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<ml[ref_ml_beta[ib]]<<" meV"<<endl;
  cout<<"---"<<endl;
}

double *X_fit;
double *Y_fit,*err_Y_fit;

//calculate the chi square
double chi2(double A,double B,double C,double D,double E,double f0,double db0,double *a)
{
  double ch2=0;

  for(int iel=0;iel<nel;iel++)
    {
      double ch2_term=pow((Y_fit[iel]-fun_fit(A,B,C,D,E,f0,db0,X_fit[iel],a[ibeta[iel]]))/err_Y_fit[iel],2);
      ch2+=ch2_term;
    }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0],B=p[1],C=p[2],D=p[3],E=p[4],f0=p[5],db0=p[6];
  double *a=p+7;
  ch=chi2(A,B,C,D,E,f0,db0,a);
}

void fit(boot &A,boot &B,boot &C,boot &D,boot &E)
{
  //copy X
  X_fit=new double[nel];
  for(int iel=0;iel<nel;iel++) X_fit[iel]=ml[iel].med();
  Y_fit=new double[nel];
  err_Y_fit=new double[nel];
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  
  minu.DefineParameter(0,"A",0.0,0.0001,0,0);
  minu.DefineParameter(1,"B",0.0,0.0001,0,0);
  minu.DefineParameter(2,"C",0,0.0001,0,0);
  minu.DefineParameter(3,"D",0,0.0001,0,0);
  minu.DefineParameter(4,"E",0,0.0001,0,0);
  //minu.FixParameter(1);
  minu.FixParameter(4);

  minu.SetFCN(chi2_wr);
  
  double C2;
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      if(iboot>0)
	minu.SetPrintLevel(-1);
      
      minu.DefineParameter(5,"f0",f0[iboot],0.0001,0,0);
      minu.DefineParameter(6,"db0",db0[iboot],0.0001,0,0);
      minu.DefineParameter(7,"a380",lat[0][iboot],0.0001,0,0);
      minu.DefineParameter(8,"a390",lat[1][iboot],0.0001,0,0);
      minu.DefineParameter(9,"a405",lat[2][iboot],0.0001,0,0);
      minu.DefineParameter(10,"a420",lat[3][iboot],0.0001,0,0);
      minu.FixParameter(5);
      minu.FixParameter(6);
      minu.FixParameter(7);
      minu.FixParameter(8);
      minu.FixParameter(9);
      minu.FixParameter(10);
      
      for(int iel=0;iel<nel;iel++)
	{
	  Y_fit[iel]=delta_fK_fr_fK.data[iel].data[iboot];
	  err_Y_fit[iel]=delta_fK_fr_fK.data[iel].err();
	}
      
      //minimize
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,A.data[iboot],dum);
      minu.GetParameter(1,B.data[iboot],dum);
      minu.GetParameter(2,C.data[iboot],dum);
      minu.GetParameter(3,D.data[iboot],dum);
      minu.GetParameter(4,E.data[iboot],dum);
      
      if(iboot==0) C2=chi2(A.data[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],lat_med);
    }
  
  //calculate the chi2
  cout<<"A = ("<<A<<"), B=("<<B<<"), C=("<<C<<"), D=("<<D<<"), E=("<<E<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<nel-4<<" = "<<C2/(nel-4)<<endl;
  
  delete[] X_fit;
  delete[] Y_fit;
  delete[] err_Y_fit;
}

void plot_funz_ml(const char *out_path,const char *title,const char *xlab,const char *ylab,boot *X,bvec &Y,bvec &par,double X_phys,double (*fun)(double,double,double,double,double,double,double,double,double),boot &chiral_extrap_cont)
{
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Chiral extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  if(fun!=NULL)
    {
      //plot the function with error
      int npoint=100;
      double X_pol[npoint];
      bvec Y_pol(npoint,nboot,njack);
      for(int ipoint=0;ipoint<npoint;ipoint++) X_pol[ipoint]=0.0599/(npoint-1)*ipoint+0.0001;
      for(int ib=0;ib<nbeta;ib++)
	{
	  bvec Y_pol(npoint,nboot,njack);
	  for(int ipoint=0;ipoint<npoint;ipoint++)
	    for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
	      Y_pol.data[ipoint].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],par[4][iboot],f0[iboot],db0[iboot],X_pol[ipoint],lat[ib][iboot]);
	  
	  out.set(1,set_fill_color[ib]);
	  out.polygon(X_pol,Y_pol);
	  out.new_set();
	  out.set(1,set_color[ib]);
	  out.set_line_size(2);
	  out.ave_line(X_pol,Y_pol);
	  out.new_set();
	}
      //plot continuum curve
      for(int ipoint=0;ipoint<npoint;ipoint++)
	for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
	  Y_pol.data[ipoint]=fun(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],par[4][iboot],f0[iboot],db0[iboot],X_pol[ipoint],0);
      //out.set(1,"magenta");
      //out.polygon(X_pol,Y_pol);
      //out.new_set();
      out.set(1,"magenta");
      out.set_line_size(3);
      out.ave_line(X_pol,Y_pol);
      out.new_set();
    }
  
  //plot the original data with error  
  for(int ib=0;ib<nbeta;ib++)
    {
      out.set(4,"none",set_symbol[ib],set_color[ib],"filled");
      out.set_legend(set_legend_fm[ib]);
      out.set_line_size(2);
      out.fout<<"@type xydy"<<endl;
      for(int iel=0;iel<nel;iel++) if(ibeta[iel]==ib) out.fout<<X[iel].med()<<" "<<Y.data[iel]<<endl;
      out.new_set();
    }
  
  //plot the extrapolated point with error
  out.set(4,"none","circle","indigo","filled");
  out.set_line_size(3);
  out.set_legend("Physical point");
  out.print_graph(X_phys,chiral_extrap_cont);
  out.new_set();
}

void plot_funz_a2(const char *out_path,const char *title,const char *xlab,const char *ylab,double *X,bvec &Y,bvec &par,double (*fun)(double,double,double,double,double,double,double,double,double),boot &chiral_extrap_cont)
{
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Continuum extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  
  //plot the function with error
  double X_pol[100],X2_pol[100];
  bvec Y_pol(100,nboot,njack);
  for(int iel=0;iel<100;iel++)
    {
      X_pol[iel]=0.1/99*iel;
      X2_pol[iel]=X_pol[iel]*X_pol[iel];
	for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
	Y_pol.data[iel].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],par[4][iboot],f0[iboot],db0[iboot],ml_phys[iboot],X_pol[iel]/hc);
    }
  out.set(1,"yellow");
  out.polygon(X2_pol,Y_pol);
  out.new_set();
  out.set(1,"red");
  out.set_line_size(2);
  out.ave_line(X2_pol,Y_pol);
  out.new_set();

  //plot the data with error
  for(int ib=0;ib<nbeta;ib++)
    {
      out.set(4,"none",set_symbol[ib],set_color[ib],"filled");
      out.set_legend(set_legend_fm[ib]);
      out.set_line_size(2);
      out.fout<<"@type xydy"<<endl;
      out.fout<<X[ib]*X[ib]<<" "<<Y.data[ib]<<endl;
      out.new_set();
    }
  
  //plot the extrapolated point with error
  out.set(4,"none","circle","indigo","filled");
  out.set_line_size(3);
  out.set_legend("Physical point");
  out.print_graph(0,chiral_extrap_cont);
  out.new_set();
}

void analysis()
{
  const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (GeV)";
  const char tag_delta_fK_fr_fK[1024]="\\xd\\0f\\sK\\N/f\\sK";
  const char tag_a2[1024]="a\\S2\\N (fm)";
  
  //perform the fit
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack),D(nboot,njack),E(nboot,njack);
  fit(A,B,C,D,E);
  cout<<endl;  
  
  //chiral extrapolation
  boot delta_fK_fr_fK_chir[nbeta];
  for(int ib=0;ib<nbeta;ib++)
    {
      //int r=ref_ml_beta[ib];
      
      delta_fK_fr_fK_chir[ib]=boot(nboot,njack);
      for(int iboot=0;iboot<nboot+1;iboot++)
	delta_fK_fr_fK_chir[ib].data[iboot]=fun_fit(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],lat[ib][iboot])*delta_m2K_phys[iboot];
    }
  
  //chiral and continuum extrapolation
  boot delta_fK_fr_fK_chir_cont(nboot,njack);
  for(int iboot=0;iboot<nboot+1;iboot++)
    delta_fK_fr_fK_chir_cont.data[iboot]=fun_fit(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0)*delta_m2K_phys[iboot];
  
  cout<<"(fK+ - fK0) / fK = ("<<delta_fK_fr_fK_chir_cont<<") "<<endl;
  cout<<"fK+/fK - 1 = ("<<0.5*delta_fK_fr_fK_chir_cont<<") "<<endl;
  
  par_res_fit=bvec(11,nboot,njack);
  par_res_fit.data[0]=A;
  par_res_fit.data[1]=B;
  par_res_fit.data[2]=C;
  par_res_fit.data[3]=D;
  par_res_fit.data[4]=E;
  par_res_fit.data[5]=f0;
  par_res_fit.data[6]=db0;
  par_res_fit.data[7]=lat[0];
  par_res_fit.data[8]=lat[1];
  par_res_fit.data[9]=lat[2];
  par_res_fit.data[10]=lat[3];
  
  plot_funz_ml("delta_fK_fr_fK_funz_ml.xmg",tag_delta_fK_fr_fK,tag_ml,tag_delta_fK_fr_fK,ml,delta_fK_fr_fK,par_res_fit,ml_phys.med(),fun_fit,delta_fK_fr_fK_chir_cont);
  plot_funz_a2("delta_fK_fr_fK_funz_a2.xmg",tag_delta_fK_fr_fK,tag_a2,tag_delta_fK_fr_fK,lat_med,delta_fK_fr_fK,par_res_fit,fun_fit,delta_fK_fr_fK_chir_cont);
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
  
  analysis();
  
  return 0;
}
