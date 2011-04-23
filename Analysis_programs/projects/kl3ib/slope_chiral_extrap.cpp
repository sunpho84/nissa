#include <include.h>
#include <iostream>

using namespace std;

int nel;
int ibeta;
int nboot=100;
int njack=10;
double *a2_M2Pi,*M2Pi;

double lat[4]={1/2.0198,1/2.3286,1/2.9419,1/3.6800};
double Zp[4]={0.411,0.437,0.477,0.501};

double clat;
double ml_chir=3.6;
double MK_phys=0.493667;
double DMK_phys=-0.006;
double MPi_phys=0.135,M2Pi_phys=MPi_phys*MPi_phys;
double a2_M2Pi_phys;

bvec da_M2K_P5P5,da_M2K_A0P5,delta_fK_fr_fK_P5P5,delta_fK_fr_fK_A0P5;

double delta_m_plot_fun(double const x,double*p)
{return DMK_phys/((p[0]*x*clat*clat+p[1])/(clat*2*MK_phys))/Zp[ibeta]*1000;}

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
  input>>nel>>nboot>>ibeta;
  clat=lat[ibeta];
  a2_M2Pi_phys=M2Pi_phys*clat*clat;
  
  a2_M2Pi=(double*)malloc(sizeof(double)*nel);
  M2Pi=(double*)malloc(sizeof(double)*nel);
  da_M2K_P5P5=bvec(nel,nboot,njack);
  da_M2K_A0P5=bvec(nel,nboot,njack);
  delta_fK_fr_fK_P5P5=bvec(nel,nboot,njack);
  delta_fK_fr_fK_A0P5=bvec(nel,nboot,njack);

  for(int i=0;i<nel;i++)
    {
      char str[1024];
      input>>str;
      jvec A(5,njack);
      A.load(str,0);
      a2_M2Pi[i]=A[0].med();
      M2Pi[i]=A[0].med()/clat/clat;
      boot_from_jack(da_M2K_P5P5[i],A[1]);
      boot_from_jack(da_M2K_A0P5[i],A[2]);
      boot_from_jack(delta_fK_fr_fK_P5P5[i],A[3]);
      boot_from_jack(delta_fK_fr_fK_A0P5[i],A[4]);
      
      cout<<a2_M2Pi[i]<<" "<<da_M2K_P5P5[i]<<endl;
    }
  cout<<"---"<<endl;
}

void fit(boot &m,boot &q,double *X,bvec &Y)
{
  //perform the fit
  linear_fit(m,q,X,Y);
  //calculate the chi2
  double C2=0;
  for(int i=0;i<nel;i++)
    {
      boot C=Y.data[i]-m*X[i]-q;
      C2+=pow(C.med()/C.err(),2);
    }
  cout<<"M = ("<<m<<"), Q=("<<q<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<nel-2<<" = "<<C2/(nel-2)<<endl;
}

void plot(const char *out_path,const char *title,const char *xlab,const char *ylab,double *X,bvec &Y,double (*fun)(double,double*),bvec &par,double X_phys,boot &chiral_extrap)
{
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Chiral extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  //plot the function with error
  out.set(1,"orange");
  out.polygon(fun,0,X[nel-1]*1.1,100,par);
  out.new_set();
  out.set(1,"red");
  out.set_line_size(2);
  out.contour(fun,0,X[nel-1]*1.1,100,par);
  //plot the original data with error  
  out.new_set();
  out.set(4,"none","square","black","filled");
  out.set_legend("Lattice data");
  out.set_line_size(2);
  out.print_graph(X,Y);
  //plot the extrapolated point with error
  out.new_set();
  out.set(4,"none","circle","green4","filled");
  out.set_line_size(3);
  out.set_legend("Physical point");
  out.print_graph(X_phys,chiral_extrap);
}

void analysis_ml(bvec &da_M2K)
{
  //perform the fit
  boot m,q;
  fit(m,q,a2_M2Pi,da_M2K);
  cout<<endl;
  
  //compute m_d - m_u
  boot da_M2K_chir=m*a2_M2Pi_phys+q;
  boot dMK_chir=da_M2K_chir/(clat*2*MK_phys);
  boot delta_m_chir=DMK_phys/dMK_chir/Zp[ibeta]*1000;
  boot mu_chir=ml_chir-delta_m_chir;
  boot md_chir=ml_chir+delta_m_chir;
  boot ratio_m_chir=mu_chir/md_chir;
  
  //plot of the fitted function
  bvec par(2,nboot,njack);par.data[0]=m;par.data[1]=q;
  plot("d_a2_M2K_chiral_extrap.xmg","\\xd\\0m=m\\sd\\N-m\\su","M\\S2\\N\\s\\xp\\N\\0 (GeV\\S2\\N)","\\xd\\0m (MeV)",a2_M2Pi,da_M2K,lin_fun,par,a2_M2Pi_phys,da_M2K_chir);
  
  //plot the md-mu plot in physical units
  bvec dMK=da_M2K/(clat*2*MK_phys);
  bvec delta_m=DMK_phys/dMK/Zp[ibeta]*1000;
  plot("delta_m_chiral_extrap.xmg","\\xd\\0m=m\\sd\\N-m\\su","M\\S2\\N\\s\\xp\\N\\0 (GeV\\S2\\N)","\\xd\\0m (MeV)",M2Pi,delta_m,delta_m_plot_fun,par,M2Pi_phys,delta_m_chir);
  
  cout<<"md - mu = "<<delta_m_chir<<" MeV"<<endl;
  cout<<"mu / md = ( "<<mu_chir<<" MeV ) / ( "<<md_chir<<" MeV ) = "<<ratio_m_chir<<endl;
  
  cout<<"---"<<endl;
}

void analysis_fK(bvec &delta_fK_fr_fK)
{
  //perform the fit
  boot m,q;
  fit(m,q,M2Pi,delta_fK_fr_fK);
  cout<<endl;
  
  //extrapolate
  boot delta_fK_fr_fK_chir=m*M2Pi_phys+q;
  
  //plot of the fitted function
  bvec par(2,nboot,njack);par.data[0]=m;par.data[1]=q;
  plot("delta_fK_fr_fK_chiral_extrap.xmg","\\xd\\0f\\sK\\N/f\\sK","M\\S2\\N\\s\\xp\\N\\0 (GeV\\S2\\N)","a\\xD\\0M\\S2\\N\\sK\\N",M2Pi,delta_fK_fr_fK,lin_fun,par,M2Pi_phys,delta_fK_fr_fK_chir);
  
  cout<<"( fK+ - fK0 ) / fK = "<<delta_fK_fr_fK_chir<<endl;
  cout<<"---"<<endl;
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      cerr<<"Error, use: "<<arg[0]<<" input"<<endl;
      exit(1);
    }
  
  read_input(arg[1]);
  
  analysis_ml(da_M2K_P5P5);
  analysis_fK(delta_fK_fr_fK_P5P5);
  
  return 0;
}
