#include "radiative_common.cpp"

const int nbeta=4;
char **base_corrs_path,**ens_name;
int nens,*T,*ibeta,*nmass;
double *cmass,*lmass;

//read data
bvec F;

bvec ml;
int ref_ml_beta[4]={-1,-1,-1,-1};

bvec par_res_fit_F;

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

int plot_iboot;
int include_a4;
int include_380;

double fun_fit_F(double A,double B,double C,double D,double ml,double a)
{
  double a1=a/lat[1][plot_iboot];
  return A*(1 + B*ml + C*a1*a1 + D*a1*a1*a1*a1);
}

void plot_funz_ml(const char *out_path,const char *title,const char *xlab,const char *ylab,bvec &X,bvec &Y,bvec &par,double X_phys,double (*fun)(double,double,double,double,double,double),boot &chiral_extrap_cont)
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
	      Y_pol.data[ipoint].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],X_pol[ipoint],lat[ib][iboot]);
	  
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
	  Y_pol.data[ipoint]=fun(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],X_pol[ipoint],0);
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
      for(int iens=0;iens<nens;iens++) if(ibeta[iens]==ib) out.fout<<X[iens].med()<<" "<<Y.data[iens]<<endl;
      out.new_set();
    }
  
  //plot the extrapolated point with error
  out.set(4,"none","circle","indigo","filled");
  out.set_line_size(3);
  out.set_legend("Physical point");
  out.print_graph(X_phys,chiral_extrap_cont);
  out.new_set();
}

bool fitting_fK;
double *X_fit;
double *Y_fit,*err_Y_fit;

//calculate the chi square
double chi2(double A,double B,double C,double D,double *a)
{
  double ch2=0;
  
  for(int iens=0;iens<nens;iens++)
    if(include_380 || ibeta[iens]!=0)
      ch2+=pow((Y_fit[iens]-fun_fit_F(A,B,C,D,X_fit[iens],a[ibeta[iens]]))/err_Y_fit[iens],2);
  
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


void plot_funz_a2(const char *out_path,const char *title,const char *xlab,const char *ylab,double *X,bvec &Y,bvec &par,double (*fun)(double,double,double,double,double,double),boot &chiral_extrap_cont)
{
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Continuum extrapolation of %s",title).c_str());
  out.axis_label(xlab,ylab);
  
  //plot the function with error
  double X_pol[100],X2_pol[100];
  bvec Y_pol(100,nboot,njack);
  for(int iens=0;iens<100;iens++)
    {
      X_pol[iens]=0.1/99*iens;
      X2_pol[iens]=X_pol[iens]*X_pol[iens];
	for(int iboot=plot_iboot=0;iboot<nboot+1;plot_iboot=iboot++)
	  Y_pol.data[iens].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],par[3][iboot],ml_phys[iboot],X_pol[iens]/hc);
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

void load_iboot(int *iboot_jack,char *ens_name)
{
  char path[1024];
  sprintf(path,"../../DATA1/%s/iboot",ens_name);
  
  FILE *fiboot=fopen(path,"r");
  if(fiboot==NULL)
    {
      perror(combine("Error opening file iboot for ensamble %s",path).c_str());
      exit(1);
    }
  int nr=fread(iboot_jack,sizeof(int),100,fiboot);
  if(nr!=100)
    {
      perror(combine("Error loading iboot data for ensamble %s",path).c_str());
      exit(1);
    }
  fclose(fiboot);
}

int main(int narg,char **arg)
{
  init_latpars();
  
  //read ensemble list, meson masses and meson name
  FILE *an_input_file=open_file("analysis_pars","r");
  char meson_name[1024];
  read_formatted_from_file_expecting((char*)&include_a4,an_input_file,"%d","include_a4");
  read_formatted_from_file_expecting((char*)&include_380,an_input_file,"%d","include_380");
  read_formatted_from_file_expecting((char*)&nens,an_input_file,"%d","nens");
  read_formatted_from_file_expecting(meson_name,an_input_file,"%s","meson_name");
  lmass=new double[nens];
  ibeta=new int[nens];
  F=bvec(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      char path[1024];
      int L;
      read_formatted_from_file((char*)&(ibeta[iens]),an_input_file,"%d","ibeta");
      read_formatted_from_file((char*)&(lmass[iens]),an_input_file,"%lg","lmass");
      read_formatted_from_file((char*)&L,an_input_file,"%d","L");
      read_formatted_from_file(path,an_input_file,"%s","path");
      
      jack BKBK(njack);
      BKBK.load(combine("../%s/M",path).c_str());

      jack VKVK(njack);
      VKVK.load(combine("../../VEC_JPSI_ETAC/%s/M_VK",path).c_str());
      
      jack P5P5(njack);
      P5P5.load(combine("../../VEC_JPSI_ETAC/%s/M_P5",path).c_str());
      
      jack k=(BKBK*BKBK-P5P5*P5P5)/(2*BKBK);
      //jack k(VKVK*VKVK-P5P5*P5P5)/(2*VKVK);
      cout<<"theta="<<k*L/sqrt(3)/M_PI<<endl;
      //load iboot
      int iboot_jack[100];
      load_iboot(iboot_jack,path);
      
      jack temp=P5P5;
      boot_from_jack(F.data[iens],temp,iboot_jack);
      F[iens]/=lat[ibeta[iens]];
    }
  fclose(an_input_file);
  
  //define ml and ref ml
  ml=bvec(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens],r=ref_ml_beta[b];
      //define ml
      cout<<iens<<" "<<b<<" "<<lmass[iens]<<endl;
      ml[iens]=lmass[iens]/lat[b]/Zp[b];
      //set the lighter mass
      if(r==-1||fabs(ml[r].med()-0.050)>fabs(ml[iens].med()-0.050)) ref_ml_beta[b]=iens;
    }
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++) if(ref_ml_beta[ib]!=-1) cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<ml[ref_ml_beta[ib]]<<" MeV"<<endl;
  cout<<"---"<<endl;

  //perform the fit
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack),D(nboot,njack);
  fit(A,B,C,D,ml,F);
  
  //chiral extrapolation
  bvec F_chir(nbeta,nboot,njack);
  boot F_chir_cont(nboot,njack);
  bvec F_estr_ml(nbeta,nboot,njack);
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      F_chir_cont.data[iboot]=fun_fit_F(A[iboot],B[iboot],C[iboot],D[iboot],ml_phys[iboot],0);
      for(int ib=0;ib<nbeta;ib++)
	{
	  int r=ref_ml_beta[ib];
	  F_chir[ib].data[iboot]=fun_fit_F(A[iboot],B[iboot],C[iboot],D[iboot],ml_phys[iboot],lat[ib][iboot]);
	  if(r!=-1)
	    F_estr_ml.data[ib].data[iboot]=F[r][iboot]*fun_fit_F(A[nboot],B[nboot],C[nboot],D[nboot],ml_phys[nboot],0)/fun_fit_F(A[nboot],B[nboot],C[nboot],D[nboot],ml[r][nboot],0);
	}
    }
  
  //chiral and continuum
  cout<<"F = "<<F_chir_cont<<" GeV"<<endl;
  cout<<endl;
  
  par_res_fit_F=bvec(4,nboot,njack);
  
  par_res_fit_F.data[0]=A;
  par_res_fit_F.data[1]=B;
  par_res_fit_F.data[2]=C;
  par_res_fit_F.data[3]=D;
  
  const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (GeV)";
  const char tag_a2[1024]="a\\S2\\N (fm)";
  double lat_med_fm[4]={lat[0].med()*hc,lat[1].med()*hc,lat[2].med()*hc,lat[3].med()*hc};
  
  plot_funz_ml("F_funz_ml.xmg",meson_name,tag_ml,meson_name,ml,F,par_res_fit_F,ml_phys.med(),fun_fit_F,F_chir_cont);
  plot_funz_a2("F_funz_a2.xmg",meson_name,tag_a2,meson_name,lat_med_fm,F_estr_ml,par_res_fit_F,fun_fit_F,F_chir_cont);

  F_chir_cont.write_to_binfile("results");
  
  return 0;
}
