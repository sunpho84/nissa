#include <include.h>
#include <iostream>

using namespace std;

int nel;
int *ibeta;
int nboot=100;
int njack=10;
const int nbeta=4;

double *aml_bare;
boot *ml;
int ref_ml_beta[nbeta];

const double hc=0.19733;
double lat[nbeta]={1/2.0198,1/2.3286,1/2.9419,1/3.6800};
double err_lat[nbeta]={lat[0]/31.5,lat[1]/36.8,lat[2]/41.9,lat[3]/44.7};
double lat_fm[nbeta]={lat[0]*hc,lat[1]*hc,lat[2]*hc,lat[3]*hc};
double Zp[nbeta]={0.411,0.437,0.477,0.501};
double err_Zp[nbeta]={0.012,0.007,0.006,0.020};

double ml_phys=3.6;
double MK_phys=0.493667;
double DMK_phys=-0.006;
double err_DMK_phys=0.0006;
//double MPi_phys=0.135,M2Pi_phys=MPi_phys*MPi_phys;

//read data
bvec da_M2K,delta_fK_fr_afK_2dm_bare;
bvec aMK,afK;
//derived data
bvec dM2K,dMK,delta_ml;
bvec delta_fK_fr_fK_2dm;
bvec M2K,fK;

bvec delta_ml_chir;
boot delta_ml_chir_cont;

bvec par_res_fit_ml,par_res_fit_fK;

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

double fun_fit(double m,double q,double k,double x,double a)
{return m*x+q+k*a*a;}

boot fun_plot_fit(boot m,boot q,boot k,double x,double a)
{return m*x+q+k*a*a;}

boot fun_plot_delta_ml(boot m,boot q,boot k,double x,double a)
{return DMK_phys/(fun_plot_fit(m,q,k,x,a)/(2*MK_phys))*1000;}

boot fun_plot_delta_fK_fr_fK(boot m,boot q,boot k,double x,double a)
{return fun_plot_fit(m,q,k,x,a)/fun_plot_delta_ml(par_res_fit_ml[0],par_res_fit_ml[1],par_res_fit_ml[2],x,a);}

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
  input>>nel>>nboot;
  ibeta=(int*)malloc(sizeof(int)*nel);
  
  aml_bare=(double*)malloc(sizeof(double)*nel);
  ml=(boot*)malloc(sizeof(boot)*nel);
  da_M2K=bvec(nel,nboot,njack);
  dM2K=bvec(nel,nboot,njack);
  dMK=bvec(nel,nboot,njack);
  delta_ml=bvec(nel,nboot,njack);
  delta_fK_fr_fK_2dm=bvec(nel,nboot,njack);
  delta_fK_fr_afK_2dm_bare=bvec(nel,nboot,njack);
  aMK=bvec(nel,nboot,njack);
  afK=bvec(nel,nboot,njack);
  M2K=bvec(nel,nboot,njack);
  fK=bvec(nel,nboot,njack);
  
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
      boot temp_Zp(nboot,njack),temp_lat(nboot,njack);
      temp_Zp.fill_gauss(Zp[ibeta[i]],err_Zp[ibeta[i]]);
      temp_lat.fill_gauss(lat[ibeta[i]],err_lat[ibeta[i]]);
      ml[i]=aml_bare[i]/temp_lat/temp_Zp*1000;
      boot_from_jack(da_M2K[i],A[1]);
      boot_from_jack(delta_fK_fr_afK_2dm_bare[i],A[4]);
      boot_from_jack(aMK[i],A[5]);
      boot_from_jack(afK[i],A[6]);
      
      //prepare renormalized quantities
      dM2K[i]=da_M2K[i]/temp_lat*temp_Zp;
      delta_fK_fr_fK_2dm[i]=delta_fK_fr_afK_2dm_bare[i]/temp_lat*temp_Zp; //propagate Zp error
      
      //calculate the physical quantity
      boot temp_DMK_phys(nboot,njack);
      temp_DMK_phys.fill_gauss(DMK_phys,err_DMK_phys);
      dMK[i]=dM2K[i]/(2*MK_phys);
      delta_ml[i]=temp_DMK_phys/dMK[i]*1000;
      M2K[i]=sqr(aMK[i]/temp_lat);
      fK[i]=afK[i]/temp_lat;
      
      cout<<ml[i]<<" "<<dM2K[i]<<endl;
      
      //set the lighter mass
      int b=ibeta[i],r=ref_ml_beta[b];
      if(r==-1||ml[r].med()>ml[i].med()) ref_ml_beta[b]=i;
    }
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++) cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<ml[ref_ml_beta[ib]]<<" MeV"<<endl;
  cout<<"---"<<endl;
}

double *X_fit,*Y_fit,*err_Y_fit;

//calculate the chi square
double chi2(double m,double q,double k)
{
  double ch2=0;

  for(int iel=0;iel<nel;iel++)
    {
      double ch2_term=pow((Y_fit[iel]-fun_fit(m,q,k,X_fit[iel],lat[ibeta[iel]]))/err_Y_fit[iel],2);
      ch2+=ch2_term;
    }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double m=p[0];double q=p[1];double k=p[2];
  ch=chi2(m,q,k);
}

void fit(boot &m,boot &q,boot &k,boot *X,bvec &Y,bool fix=0)
{
  //copy X
  X_fit=new double[nel];
  for(int iel=0;iel<nel;iel++) X_fit[iel]=X[iel].med();
  Y_fit=new double[nel];
  err_Y_fit=new double[nel];
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.DefineParameter(0,"m",0,0.0001,0,0);
  minu.DefineParameter(1,"q",0,0.0001,0,0);
  minu.DefineParameter(2,"k",0,0.0001,0,0);
  if(fix==1) minu.FixParameter(0);
  minu.SetFCN(chi2_wr);

  double C2;
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      for(int iel=0;iel<nel;iel++)
	{
	  Y_fit[iel]=Y.data[iel].data[iboot];
	  err_Y_fit[iel]=Y.data[iel].err();
	}
      
      //if(iboot==0) minu.SetPrintLevel(+1);
      //else
      minu.SetPrintLevel(-1);
      
      //minimize
      minu.Migrad();
      
      //get back parameters
      double dum;
      minu.GetParameter(0,m.data[iboot],dum);
      minu.GetParameter(1,q.data[iboot],dum);
      minu.GetParameter(2,k.data[iboot],dum);
      
      if(iboot==0) C2=chi2(m.data[iboot],q[iboot],k[iboot]);
    }
  
  //calculate the chi2
  cout<<"M = ("<<m<<"), Q=("<<q<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<nel-3<<" = "<<C2/(nel-3)<<endl;
  
  delete[] X_fit;
  delete[] Y_fit;
  delete[] err_Y_fit;
}

void plot_funz_ml(const char *out_path,const char *title,const char *xlab,const char *ylab,boot *X,bvec &Y,bvec &par,boot (*plot_fun)(boot,boot,boot,double,double),double X_phys,boot &chiral_extrap_cont)
{
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Chiral extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  //plot the function with error
  int npoint=100;
  double X_pol[npoint];
  bvec Y_pol(npoint,nboot,njack);
  for(int ipoint=0;ipoint<npoint;ipoint++) X_pol[ipoint]=60.0/(npoint-1)*ipoint;
  for(int ib=0;ib<nbeta;ib++)
    {
      bvec Y_pol(npoint,nboot,njack);
      for(int ipoint=0;ipoint<npoint;ipoint++) Y_pol.data[ipoint]=plot_fun(par[0],par[1],par[2],X_pol[ipoint],lat[ib]);
      out.set(1,set_fill_color[ib]);
      out.polygon(X_pol,Y_pol);
      out.new_set();
      out.set(1,set_color[ib]);
      out.set_line_size(2);
      out.ave_line(X_pol,Y_pol);
      out.new_set();
    }
  //plot continuum curve
  for(int ipoint=0;ipoint<npoint;ipoint++) Y_pol.data[ipoint]=plot_fun(par[0],par[1],par[2],X_pol[ipoint],lat[0]*0);
  //out.set(1,"magenta");
  //out.polygon(X_pol,Y_pol);
  //out.new_set();
  out.set(1,"magenta");
  out.set_line_size(3);
  out.ave_line(X_pol,Y_pol);
  out.new_set();
  
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

void plot_funz_a2(const char *out_path,const char *title,const char *xlab,const char *ylab,double *X,bvec &Y,bvec &par,boot (*plot_fun)(boot,boot,boot,double,double),boot &chiral_extrap_cont)
{
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Continuum extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  
  //plot the function with error
  if(plot_fun!=NULL)
    {
      double X_pol[100],X2_pol[100];
      bvec Y_pol(100,nboot,njack);
      for(int iel=0;iel<100;iel++)
	{
	  X_pol[iel]=0.1/99*iel;
	  X2_pol[iel]=X_pol[iel]*X_pol[iel];
	  Y_pol.data[iel]=plot_fun(par[0],par[1],par[2],ml_phys/1000,X_pol[iel]/0.197);
	}
      out.set(1,"yellow");
      out.polygon(X2_pol,Y_pol);
      out.new_set();
      out.set(1,"red");
      out.set_line_size(2);
      out.ave_line(X2_pol,Y_pol);
      out.new_set();
    }
  
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

void analysis_ml()
{
  const char tag_dM2K[1024]="\\xd\\0M\\S2\\N\\sK";
  const char tag_delta_ml[1024]="\\xD\\0m\\sl\\N\\SMS,2GeV\\N (MeV)";
  
  const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (MeV)";
  const char tag_a2[1024]="a\\S2\\N (fm)";
  
  //perform the fit
  boot m(nboot,njack),q(nboot,njack),k(nboot,njack);
  fit(m,q,k,ml,dM2K,1);
  cout<<endl;  
  
  //compute m_d - m_u
  boot dMK_chir[nbeta],dMK_estr_ml[nbeta];
  bvec dM2K_estr_ml(nbeta,nboot,njack),dM2K_chir(nbeta,nboot,njack);
  bvec delta_ml_estr_ml(nbeta,nboot,njack);
  delta_ml_chir=bvec(nbeta,nboot,njack);
  
  for(int ib=0;ib<nbeta;ib++)
    {
      int r=ref_ml_beta[ib];
      
      dM2K_chir[ib]=m*ml_phys+q+k*lat[ib]*lat[ib];
      dM2K_estr_ml[ib]=dM2K_chir[ib]*(m*ml_phys+q+k*lat[ib]*lat[ib])/(m*ml[r]+q+k*lat[ib]*lat[ib]);
      
      dMK_chir[ib]=dM2K_chir[ib]/(2*MK_phys);
      dMK_estr_ml[ib]=dM2K_estr_ml[ib]/(2*MK_phys);
      
      delta_ml_chir[ib]=DMK_phys/dMK_chir[ib]*1000;
      delta_ml_estr_ml[ib]=DMK_phys/dMK_estr_ml[ib]*1000;
    }
  
  boot dM2K_chir_cont=m*ml_phys+q;
  boot dMK_chir_cont=dM2K_chir_cont/(2*MK_phys);
  delta_ml_chir_cont=DMK_phys/dMK_chir_cont*1000;
  
  boot mu_chir_cont=ml_phys-delta_ml_chir_cont;
  boot md_chir_cont=ml_phys+delta_ml_chir_cont;
  boot ratio_m_chir_cont=mu_chir_cont/md_chir_cont;
  
  //plot of the fitted function
  par_res_fit_ml=bvec(3,nboot,njack);
  par_res_fit_ml.data[0]=m;par_res_fit_ml.data[1]=q;par_res_fit_ml.data[2]=k;
  plot_funz_ml("dM2K_funz_ml.xmg",tag_dM2K,tag_ml,tag_dM2K,ml,dM2K,par_res_fit_ml,fun_plot_fit,ml_phys,dM2K_chir_cont);
  plot_funz_a2("dM2K_funz_a2.xmg",tag_dM2K,tag_a2,tag_dM2K,lat_fm,dM2K_estr_ml,par_res_fit_ml,fun_plot_fit,dM2K_chir_cont);
  
  //plot the md-mu plot in physical units vs ml and vs a2
  plot_funz_ml("delta_ml_funz_ml.xmg",tag_delta_ml,tag_ml,tag_delta_ml,ml,delta_ml,par_res_fit_ml,fun_plot_delta_ml,ml_phys,delta_ml_chir_cont);
  plot_funz_a2("delta_ml_funz_a2.xmg",tag_delta_ml,tag_a2,tag_delta_ml,lat_fm,delta_ml_estr_ml,par_res_fit_ml,fun_plot_delta_ml,delta_ml_chir_cont);
  
  cout<<"md - mu = "<<delta_ml_chir_cont<<" MeV"<<endl;
  cout<<"mu / md = ( "<<mu_chir_cont<<" MeV ) / ( "<<md_chir_cont<<" MeV ) = "<<ratio_m_chir_cont<<endl;
  
  cout<<"---"<<endl;
  
  plot_funz_ml("M2K_funz_ml.xmg",tag_dM2K,tag_ml,tag_dM2K,ml,M2K,par_res_fit_ml,NULL,ml_phys,dM2K_chir_cont);
}

void analysis_fK()
{
  const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (MeV)";
  const char tag_delta_fK_fr_fK[1024]="\\xd\\0f\\sK\\N/\\f\\sK";
  const char tag_delta_fK_fr_fK_2dm[1024]="\\xd\\0f\\sK\\N/\\f\\sK aggiustare";
  
  //perform the fit
  boot m(nboot,njack),q(nboot,njack),k(nboot,njack);
  fit(m,q,k,ml,delta_fK_fr_fK_2dm,0);
  cout<<endl;  
  
  //compute m_d - m_u
  boot delta_fK_fr_fK_2dm_chir[nbeta];
  boot delta_fK_fr_fK_chir[nbeta];
  for(int ib=0;ib<nbeta;ib++)
    {
      delta_fK_fr_fK_2dm_chir[ib]=m*ml_phys+q+k*lat[ib]*lat[ib];
      delta_fK_fr_fK_chir[ib]=delta_fK_fr_fK_2dm_chir[ib]*delta_ml_chir[ib]/1000;
    }
  
  //extrapolate
  boot delta_fK_fr_fK_2dm_chir_cont=m*ml_phys+q;
  boot delta_fK_fr_fK_chir_cont=delta_fK_fr_fK_2dm_chir_cont*delta_ml_chir_cont/1000;
  
  //plot of the fitted function
  bvec delta_fK_fr_fK=delta_fK_fr_fK_2dm*delta_ml/1000;
  par_res_fit_fK=bvec(3,nboot,njack);
  par_res_fit_fK.data[0]=m;par_res_fit_fK.data[1]=q;par_res_fit_fK.data[2]=k;
  plot_funz_ml("delta_fK_fr_fK_2dm_chiral_extrap.xmg",tag_delta_fK_fr_fK_2dm,tag_ml,tag_delta_fK_fr_fK_2dm,ml,delta_fK_fr_fK_2dm,par_res_fit_fK,fun_plot_fit,ml_phys,delta_fK_fr_fK_2dm_chir_cont);
  plot_funz_ml("delta_fK_fr_fK_chiral_extrap.xmg",tag_delta_fK_fr_fK,tag_ml,tag_delta_fK_fr_fK,ml,delta_fK_fr_fK,par_res_fit_fK,fun_plot_delta_fK_fr_fK,ml_phys,delta_fK_fr_fK_chir_cont);
  
  cout<<"( fK+ - fK0 ) / fK = "<<delta_fK_fr_fK_chir_cont<<endl;
  cout<<"---"<<endl;
  
  plot_funz_ml("fK_funz_ml",tag_delta_fK_fr_fK_2dm,tag_ml,tag_delta_fK_fr_fK_2dm,ml,fK,par_res_fit_fK,NULL,ml_phys,delta_fK_fr_fK_2dm_chir_cont);
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      cerr<<"Error, use: "<<arg[0]<<" input"<<endl;
      exit(1);
    }
  
  read_input(arg[1]);
  
  analysis_ml();
  analysis_fK();
  
  return 0;
}
