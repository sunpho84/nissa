#include <include.h>
#include <iostream>
#include <math.h>

using namespace std;

int nel;
int *ibeta;
const int nboot=100;
int njack=10;
const int nbeta=4;

double *aml_bare;
boot *ml;
int ref_ml_beta[nbeta];

const double hc=0.19733;
double lat_med[nbeta]={1/2.0198,1/2.3286,1/2.9419,1/3.6800};
double lat_err[nbeta]={lat_med[0]/31.5,lat_med[1]/36.8,lat_med[2]/41.9,lat_med[3]/44.7};
double lat_med_fm[nbeta]={lat_med[0]*hc,lat_med[1]*hc,lat_med[2]*hc,lat_med[3]*hc};
double Zp_med[nbeta]={0.411,0.437,0.477,0.501};
double Zp_err[nbeta]={0.012,0.007,0.006,0.020};

double ml_phys_med=3.6/1000;
double ml_phys_err=0.2/1000;
double MK_phys=0.493667;
double delta_MK_phys_med=-0.006;
double delta_MK_phys_err=0.0006;
//double MPi_phys=0.135,M2Pi_phys=MPi_phys*MPi_phys;

//read data
bvec da_M2K,delta_fK_fr_afK_2dm_bare;
bvec aMK,afK;
//derived data
bvec dM2K,dMK,delta_ml;
bvec delta_fK,delta_fK_fr_fK_2dm;
bvec dfK;
bvec M2K,fK;

bvec delta_ml_chir;
boot delta_ml_chir_cont;

boot ml_phys,f0,db0,lat[nbeta],Zp[nbeta];
boot delta_MK_phys;

bvec par_res_fit_fK;
bvec par_res_fit_M2K;

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

int plot_iboot;

double fun_fit_fK(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{
  double xill=db0*ml/(16*M_PI*M_PI*f0*f0);
  return A*(1-3.0/4*xill*log(xill)+B*xill+(C+E*xill)*a*a+D*xill*xill);
}

double fun_fit_dfK(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{
  double dxill_dml=db0/(16*M_PI*M_PI*f0*f0);
  double xill=dxill_dml*ml;
  return A+B*xill*log(xill)+C*xill+D*a*a;
}


double fun_fit_M2K(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{return A + B*ml + C*a*a + D*ml*ml + E*a*a*ml;}

double fun_fit_dM2K(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{return B + E*a*a;}

double fun_fit_delta_ml(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{return delta_MK_phys[plot_iboot]/(fun_fit_dM2K(A,B,C,D,E,f0,db0,ml,a)/2/pow(fun_fit_M2K(A,B,C,D,E,f0,db0,ml,a),0.5));}

double fun_fit_delta_fK(double A,double B,double C,double D,double E,double f0,double db0,double ml,double a)
{
  return fun_fit_dfK(A,B,C,D,E,f0,db0,ml,a)*fun_fit_delta_ml(par_res_fit_M2K[0][plot_iboot],par_res_fit_M2K[1][plot_iboot],par_res_fit_M2K[2][plot_iboot],par_res_fit_M2K[3][plot_iboot],par_res_fit_M2K[4][plot_iboot],f0,db0,ml,a);
}


int rni(int n)
{return rand()%n;}

void boot_from_jack(boot &out,jack &in)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[rni(njack)];
  out.data[nboot]=in.data[njack];
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
  da_M2K=bvec(nel,nboot,njack);
  dM2K=bvec(nel,nboot,njack);
  dMK=bvec(nel,nboot,njack);
  delta_ml=bvec(nel,nboot,njack);
  delta_fK=bvec(nel,nboot,njack);
  delta_fK_fr_fK_2dm=bvec(nel,nboot,njack);
  delta_fK_fr_afK_2dm_bare=bvec(nel,nboot,njack);
  dfK=bvec(nel,nboot,njack);
  aMK=bvec(nel,nboot,njack);
  afK=bvec(nel,nboot,njack);
  M2K=bvec(nel,nboot,njack);
  fK=bvec(nel,nboot,njack);
  for(int ib=0;ib<nbeta;ib++)
    {
      lat[ib]=boot(nboot,njack);
      lat[ib].fill_gauss(lat_med[ib],lat_err[ib],235235+ib);
      Zp[ib]=boot(nboot,njack);
      Zp[ib].fill_gauss(Zp_med[ib],Zp_err[ib],87466+ib);
    }
  delta_MK_phys=boot(nboot,njack);
  ml_phys=boot(nboot,njack);
  f0=boot(nboot,njack);
  db0=boot(nboot,njack);
  
  delta_MK_phys.fill_gauss(delta_MK_phys_med,delta_MK_phys_err,678456);
  ml_phys.fill_gauss(ml_phys_med,ml_phys_err,97635);
  f0.fill_gauss(121.9019/1000,0.1668/1000,765473457);
  db0.fill_gauss(5.2756,0.2078,6732576);
  

  //reset index list of lighter ml
  for(int ib=0;ib<nbeta;ib++) ref_ml_beta[ib]=-1;
  
  for(int i=0;i<nel;i++)
    {
      //read input
      char str[1024];
      input>>str>>ibeta[i]>>aml_bare[i];
      
      //read data
      jvec A(7,njack);
      A.load(str,0);
      
      //prepare boot
      ml[i]=aml_bare[i]/lat[ibeta[i]]/Zp[ibeta[i]];
      boot_from_jack(da_M2K[i],A[1]);
      boot_from_jack(delta_fK_fr_afK_2dm_bare[i],A[4]);
      boot_from_jack(aMK[i],A[5]);
      boot_from_jack(afK[i],A[6]);
      
      //prepare renormalized quantities
      dM2K[i]=da_M2K[i]/lat[ibeta[i]]*Zp[ibeta[i]];
      delta_fK_fr_fK_2dm[i]=delta_fK_fr_afK_2dm_bare[i]*lat[ibeta[i]]*Zp[ibeta[i]]; //propagate Zp error
      //calculate the physical quantity
      dMK[i]=dM2K[i]/(2*MK_phys);
      delta_ml[i]=delta_MK_phys/dMK[i];
      M2K[i]=sqr(aMK[i]/lat[ibeta[i]]);
      fK[i]=afK[i]/lat[ibeta[i]];
      dfK[i]=delta_fK_fr_afK_2dm_bare[i]*afK[i]*Zp[ibeta[i]];
      delta_fK[i]=dfK[i]*delta_ml[i];
      
      cout<<ml[i]<<" "<<dM2K[i]<<" "<<fK[i]<<endl;
      
      //set the lighter mass
      int b=ibeta[i],r=ref_ml_beta[b];
      if(r==-1||ml[r].med()>ml[i].med()) ref_ml_beta[b]=i;
    }
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++) cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<ml[ref_ml_beta[ib]]<<" MeV"<<endl;
  cout<<"---"<<endl;
}

bool fitting_fK;
double *X_fit;
double *Y_fit,*err_Y_fit;
double *Z_fit,*err_Z_fit;

//calculate the chi square
double chi2(double A,double B,double C,double D,double E,double f0,double db0,double *a)
{
  double ch2=0;

  for(int iel=0;iel<nel;iel++)
    {
      if(fitting_fK==1)
	{
	  //double ch2_fK_term=pow((Y_fit[iel]-fun_fit_fK(A,B,C,D,E,f0,db0,X_fit[iel],a[ibeta[iel]]))/err_Y_fit[iel],2);
	  double ch2_dfK_term=pow((Z_fit[iel]-fun_fit_dfK(A,B,C,D,E,f0,db0,X_fit[iel],a[ibeta[iel]]))/err_Z_fit[iel],2);
	  ch2+=ch2_dfK_term;//+ch2_fK_term;
	}
      else
	{
	  //double ch2_M2K_term=pow((Y_fit[iel]-fun_fit_M2K(A,B,C,D,E,f0,db0,X_fit[iel],a[ibeta[iel]]))/err_Y_fit[iel],2);
	  double ch2_dM2K_term=pow((Z_fit[iel]-fun_fit_dM2K(A,B,C,D,E,f0,db0,X_fit[iel],a[ibeta[iel]]))/err_Z_fit[iel],2);
	  ch2+=ch2_dM2K_term;//+ch2_M2K_term;
	}
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

void fit(boot &A,boot &B,boot &C,boot &D,boot &E,boot *X,bvec &Y,bvec &Z)
{
  //copy X
  X_fit=new double[nel];
  for(int iel=0;iel<nel;iel++) X_fit[iel]=X[iel].med();
  Y_fit=new double[nel];
  err_Y_fit=new double[nel];
  Z_fit=new double[nel];
  err_Z_fit=new double[nel];
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  
  minu.DefineParameter(0,"A",0.0,0.0001,0,0);
  minu.DefineParameter(1,"B",0.01,0.0001,0,0);
  minu.DefineParameter(2,"C",0,0.0001,0,0);
  minu.DefineParameter(3,"D",0,0.0001,0,0);
  minu.DefineParameter(4,"E",0,0.0001,0,0);
  if(fitting_fK==0)
    {
      minu.FixParameter(0);
      minu.FixParameter(2);
      minu.FixParameter(3);
    }
  else
    {
      minu.FixParameter(1);
      minu.FixParameter(4);
    }
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
	  Y_fit[iel]=Y.data[iel].data[iboot];
	  err_Y_fit[iel]=Y.data[iel].err();
	  Z_fit[iel]=Z.data[iel].data[iboot];
	  err_Z_fit[iel]=Z.data[iel].err();
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
  cout<<"Chi2 = "<<C2<<" / "<<2*nel-5<<" = "<<C2/(2*nel-5)<<endl;
  
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

void analysis_M2K()
{
  const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (GeV)";
  const char tag_delta_ml[1024]="\\xd\\0m\\sl\\N\\SMS,2GeV\\N (GeV)";
  const char tag_M2K[1024]="M\\S2\\N\\sK";
  const char tag_dM2K[1024]="\\xd\\0M\\S2\\N\\sK";
  const char tag_a2[1024]="a\\S2\\N (fm)";
  
  //perform the fit
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack),D(nboot,njack),E(nboot,njack);
  fitting_fK=0;
  fit(A,B,C,D,E,ml,M2K,dM2K);
  cout<<endl;  
  
  //chiral extrapolation
  boot M2K_chir[nbeta],M2K_chir_cont(nboot,njack),dM2K_chir_cont(nboot,njack);
  bvec dM2K_estr_ml(nbeta,nboot,njack);
  for(int ib=0;ib<nbeta;ib++)
    {
      M2K_chir[ib]=boot(nboot,njack);
      for(int iboot=0;iboot <nboot+1;iboot++)
	{
	  int r=ref_ml_beta[ib];
	  M2K_chir_cont.data[iboot]=fun_fit_M2K(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0);
	  dM2K_chir_cont.data[iboot]=fun_fit_dM2K(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0);
	  M2K_chir[ib].data[iboot]=fun_fit_M2K(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],lat[ib][iboot]);
	  dM2K_estr_ml.data[ib].data[iboot]=dM2K[r][iboot]*fun_fit_dM2K(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0)/fun_fit_dM2K(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml[r][iboot],0);
	}
    }
  
  //calculate md-mu
  //boot dMK_chir_cont=dM2K_chir_cont/(2*sqrt(M2K_chir_cont));
  boot dMK_chir_cont=dM2K_chir_cont/(2*MK_phys);
  delta_ml_chir_cont=delta_MK_phys/dMK_chir_cont;
  boot mu_chir_cont=ml_phys-delta_ml_chir_cont/2;
  boot md_chir_cont=ml_phys+delta_ml_chir_cont/2;
  boot mu_fr_md_chir_cont=mu_chir_cont/md_chir_cont;
  
  //chiral and continuum
  cout<<"MK = ("<<sqrt(M2K_chir_cont)*1000<<") MeV"<<endl;
  cout<<"md - mu = ("<<delta_ml_chir_cont*1000<<") Mev"<<endl;
  cout<<"mu / md = ("<<mu_chir_cont*1000<<") / ("<<md_chir_cont*1000<<") = "<<mu_fr_md_chir_cont<<endl;
    
  par_res_fit_M2K=bvec(11,nboot,njack);
  
  par_res_fit_M2K.data[0]=A;
  par_res_fit_M2K.data[1]=B;
  par_res_fit_M2K.data[2]=C;
  par_res_fit_M2K.data[3]=D;
  par_res_fit_M2K.data[4]=E;
  par_res_fit_M2K.data[5]=f0;
  par_res_fit_M2K.data[6]=db0;
  par_res_fit_M2K.data[7]=lat[0];
  par_res_fit_M2K.data[8]=lat[1];
  par_res_fit_M2K.data[9]=lat[2];
  par_res_fit_M2K.data[10]=lat[3];
  
  plot_funz_ml("M2K_funz_ml.xmg",tag_M2K,tag_ml,tag_M2K,ml,M2K,par_res_fit_M2K,ml_phys.med(),fun_fit_M2K,M2K_chir_cont);
  plot_funz_ml("dM2K_funz_ml.xmg",tag_dM2K,tag_ml,tag_dM2K,ml,dM2K,par_res_fit_M2K,ml_phys.med(),fun_fit_dM2K,dM2K_chir_cont);
  plot_funz_ml("delta_ml_funz_ml.xmg",tag_delta_ml,tag_ml,tag_delta_ml,ml,delta_ml,par_res_fit_M2K,ml_phys.med(),fun_fit_delta_ml,delta_ml_chir_cont);
  plot_funz_a2("dM2K_funz_a2.xmg",tag_dM2K,tag_a2,tag_dM2K,lat_med_fm,dM2K_estr_ml,par_res_fit_M2K,fun_fit_dM2K,dM2K_chir_cont);
}

void analysis_fK()
{
  const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (GeV)";
  const char tag_fK[1024]="f\\sK (GeV)";
  const char tag_delta_fK[1024]="\\xd\\0f\\sK (GeV)";
  const char tag_dfK[1024]="\\xd\\0f\\sK\\N/\\xd\\0m";
  const char tag_a2[1024]="a\\S2\\N (fm)";
  
  //perform the fit
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack),D(nboot,njack),E(nboot,njack);
  fitting_fK=1;
  fit(A,B,C,D,E,ml,fK,dfK);
  cout<<endl;  
  
  //chiral extrapolation
  boot fK_chir[nbeta],fK_chir_cont(nboot,njack),dfK_chir_cont(nboot,njack);
  bvec dfK_estr_ml(nbeta,nboot,njack);
  for(int ib=0;ib<nbeta;ib++)
    {
      int r=ref_ml_beta[ib];
      
      fK_chir[ib]=boot(nboot,njack);
      for(int iboot=0;iboot <nboot+1;iboot++)
	{
	  fK_chir_cont.data[iboot]=fun_fit_fK(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0);
	  dfK_chir_cont.data[iboot]=fun_fit_dfK(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0);
	  fK_chir[ib].data[iboot]=fun_fit_fK(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],lat[ib][iboot]);
	  dfK_estr_ml.data[ib].data[iboot]=dfK[r][iboot]*fun_fit_dfK(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml_phys[iboot],0)/fun_fit_dfK(A[iboot],B[iboot],C[iboot],D[iboot],E[iboot],f0[iboot],db0[iboot],ml[r][iboot],0);
	}
    }
  
  //chiral and continuum
  boot delta_fK_chir_cont=dfK_chir_cont*delta_ml_chir_cont;
  cout<<"fK = ("<<fK_chir_cont*1000<<") MeV"<<endl;
  cout<<"fK+ - fK0 = ("<<delta_fK_chir_cont*1000<<") MeV"<<endl;
  cout<<"(fK+/fK)-1 = ("<<delta_fK_chir_cont/fK_chir_cont/2<<")"<<endl;
  
  par_res_fit_fK=bvec(11,nboot,njack);
  par_res_fit_fK.data[0]=A;
  par_res_fit_fK.data[1]=B;
  par_res_fit_fK.data[2]=C;
  par_res_fit_fK.data[3]=D;
  par_res_fit_fK.data[4]=E;
  par_res_fit_fK.data[5]=f0;
  par_res_fit_fK.data[6]=db0;
  par_res_fit_fK.data[7]=lat[0];
  par_res_fit_fK.data[8]=lat[1];
  par_res_fit_fK.data[9]=lat[2];
  par_res_fit_fK.data[10]=lat[3];
  
  plot_funz_ml("fK_funz_ml.xmg",tag_fK,tag_ml,tag_fK,ml,fK,par_res_fit_fK,ml_phys.med(),fun_fit_fK,fK_chir_cont);
  plot_funz_ml("dfK_funz_ml.xmg",tag_dfK,tag_ml,tag_dfK,ml,dfK,par_res_fit_fK,ml_phys.med(),fun_fit_dfK,dfK_chir_cont);
  plot_funz_ml("delta_fK_funz_ml.xmg",tag_delta_fK,tag_ml,tag_delta_fK,ml,delta_fK,par_res_fit_fK,ml_phys.med(),fun_fit_delta_fK,delta_fK_chir_cont);
  plot_funz_a2("dfK_funz_a2.xmg",tag_dfK,tag_a2,tag_dfK,lat_med_fm,dfK_estr_ml,par_res_fit_fK,fun_fit_dfK,dfK_chir_cont);
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      cerr<<"Error, use: "<<arg[0]<<" input"<<endl;
      exit(1);
    }
  
  read_input(arg[1]);
  
  analysis_M2K();
  cout<<"---"<<endl;
  analysis_fK();
  
  return 0;
}
