#include "../common.cpp"

const int nbeta=4;
char **base_corrs_path,**ens_name;
int nens,*T,*ibeta,*iml_un,*nlights,*nmass;
double **mass;

//read data
bvec phiD_s;

bvec ml;
int ref_ml_beta[4]={-1,-1,-1,-1};

bvec par_res_fit_phiD_s;

const char set_color[nbeta][1024]={"black","blue","red","green4"};
const char set_fill_color[nbeta][1024]={"grey","turquoise","yellow","green"};
const char set_symbol[nbeta][1024]={"square","circle","triangle up","triangle left"};
const char set_legend[nbeta][1024]={"\\xb\\0=3.80","\\xb\\0=3.90","\\xb\\0=4.05","\\xb\\0=4.20"};
const char set_legend_fm[nbeta][1024]={"a = 0.098 fm","a = 0.085 fm","a = 0.067 fm","a = 0.054 fm"};

int glb_iboot;

double fun_fit_phiD_s(double A,double B,double C,double ml,double a)
{return A * ( 1 + B*ml + C*a*a );}

void plot_funz_ml(const char *out_path,const char *title,const char *xlab,const char *ylab,bvec &X,bvec &Y,bvec &par,double X_phys,double (*fun)(double,double,double,double,double),boot &chiral_extrap_cont)
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
	    for(int iboot=glb_iboot=0;iboot<nboot+1;glb_iboot=iboot++)
	      Y_pol.data[ipoint].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],X_pol[ipoint],lat[ib][iboot]);
	  
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
	for(int iboot=glb_iboot=0;iboot<nboot+1;glb_iboot=iboot++)
	  Y_pol.data[ipoint]=fun(par[0][iboot],par[1][iboot],par[2][iboot],X_pol[ipoint],0);
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
double chi2(double A,double B,double C,double *a)
{
  double ch2=0;

  for(int iens=0;iens<nens;iens++)
    {
      double ch2_phiD_s_term=pow((Y_fit[iens]-fun_fit_phiD_s(A,B,C,X_fit[iens],a[ibeta[iens]]))/err_Y_fit[iens],2);
      ch2+=ch2_phiD_s_term;
    }
  
  return ch2;
}

//wrapper for the calculation of the chi2
void chi2_wr(int &npar,double *fuf,double &ch,double *p,int flag)
{
  double A=p[0],B=p[1],C=p[2];
  double *a=p+3;
  ch=chi2(A,B,C,a);
}

void fit(boot &A,boot &B,boot &C,bvec &X,bvec &Y)
{
  //copy X
  X_fit=new double[nens];
  for(int iens=0;iens<nens;iens++) X_fit[iens]=X[iens].med();
  Y_fit=new double[nens];
  err_Y_fit=new double[nens];
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  
  minu.DefineParameter(0,"A",0.0,0.0001,0,0);
  minu.DefineParameter(1,"B",0.0,0.0001,0,0);
  minu.DefineParameter(2,"C",0.0,0.0001,0,0);
  
  minu.SetFCN(chi2_wr);
  
  double C2;
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      if(iboot>0)
        minu.SetPrintLevel(-1);
      
      minu.DefineParameter(3,"a380",lat[0][iboot],0.0001,0,0);
      minu.DefineParameter(4,"a390",lat[1][iboot],0.0001,0,0);
      minu.DefineParameter(5,"a405",lat[2][iboot],0.0001,0,0);
      minu.DefineParameter(6,"a420",lat[3][iboot],0.0001,0,0);
      minu.FixParameter(3);
      minu.FixParameter(4);
      minu.FixParameter(5);
      minu.FixParameter(6);
      
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
      
      if(iboot==0) C2=chi2(A.data[iboot],B[iboot],C[iboot],lat_med);
    }
  
  //calculate the chi2
  cout<<"A = ("<<A<<"), B=("<<B<<"), C=("<<C<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<nens-3<<" = "<<C2/(nens-3)<<endl;
  
  delete[] X_fit;
  delete[] Y_fit;
  delete[] err_Y_fit;
}


void plot_funz_a2(const char *out_path,const char *title,const char *xlab,const char *ylab,double *X,bvec &Y,bvec &par,double (*fun)(double,double,double,double,double),boot &chiral_extrap_cont)
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
	for(int iboot=glb_iboot=0;iboot<nboot+1;glb_iboot=iboot++)
	  Y_pol.data[iens].data[iboot]=fun(par[0][iboot],par[1][iboot],par[2][iboot],ml_phys[iboot],X_pol[iens]/hc);
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

int main(int narg,char **arg)
{
  init_latpars();
  
  //read ensemble list, meson phiD_s and meson name
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],meson_phi_file[1024],meson_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(meson_phi_file,an_input_file,"%s","meson_phi_file");
  read_formatted_from_file_expecting(meson_name,an_input_file,"%s","meson_name");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,iml_un,nlights,ens_list_path);
  
  //define ml and ref ml, and chiral logs
  ml=bvec(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens],r=ref_ml_beta[b];
      //define ml
      cout<<iens<<" "<<b<<" "<<iml_un[iens]<<" "<<mass[iens][iml_un[iens]]<<endl;
      ml[iens]=mass[iens][iml_un[iens]]/lat[b]/Zp[b];
      
      //set the lighter mass
      if(r==-1||ml[r].med()>ml[iens].med()) ref_ml_beta[b]=iens;
    }
  cout<<"---"<<endl;
  for(int ib=0;ib<nbeta;ib++) if(ref_ml_beta[ib]!=-1) cout<<"Ref "<<ib<<" = "<<ref_ml_beta[ib]<<", "<<ml[ref_ml_beta[ib]]<<" MeV"<<endl;
  cout<<"---"<<endl;

  //load data
  phiD_s=bvec(nens,nboot,njack);
  phiD_s.load(meson_phi_file,0);  

  //perform the fit
  boot A(nboot,njack),B(nboot,njack),C(nboot,njack);
  fit(A,B,C,ml,phiD_s);
  cout<<endl;
  
  //chiral extrapolation
  boot phiD_s_chir[nbeta],phiD_s_chir_cont(nboot,njack);
  bvec phiD_s_estr_ml(nbeta,nboot,njack),fD_s_estr_ml(nbeta,nboot,njack);
  for(int ib=0;ib<nbeta;ib++)
    {
      phiD_s_chir[ib]=boot(nboot,njack);
      for(int iboot=0;iboot<nboot+1;iboot++)
	{
	  int r=ref_ml_beta[ib];
	  phiD_s_chir_cont.data[iboot]=fun_fit_phiD_s(A[iboot],B[iboot],C[iboot],ml_phys[iboot],0);
	  phiD_s_chir[ib].data[iboot]=fun_fit_phiD_s(A[iboot],B[iboot],C[iboot],ml_phys[iboot],lat[ib][iboot]);
	  
	  if(r!=-1) phiD_s_estr_ml.data[ib].data[iboot]=phiD_s[r][iboot]*fun_fit_phiD_s(A[iboot],B[iboot],C[iboot],ml_phys[iboot],0)/fun_fit_phiD_s(A[iboot],B[iboot],C[iboot],ml[r][iboot],0);
	}
    }
  boot fD_s_chir_cont=phiD_s_chir_cont/sqrt(mD_s_phys);
  
  //chiral and continuum
  cout<<"phiD_s = ("<<phiD_s_chir_cont*1000<<") MeV"<<endl;
  cout<<"fD_s = ("<<fD_s_chir_cont*1000<<") MeV"<<endl;
  phiD_s_chir_cont.write_to_binfile("results_phiD_s");
  fD_s_chir_cont.write_to_binfile("results_fD_s");
  
  par_res_fit_phiD_s=bvec(3,nboot,njack);
  
  par_res_fit_phiD_s.data[0]=A;
  par_res_fit_phiD_s.data[1]=B;
  par_res_fit_phiD_s.data[2]=C;
  
  const char tag_ml[1024]="m\\sl\\N\\SMS,2GeV\\N (GeV)";
  const char tag_a2[1024]="a\\S2\\N (fm)";
  plot_funz_ml("phiD_s_funz_ml.xmg",meson_name,tag_ml,meson_name,ml,phiD_s,par_res_fit_phiD_s,ml_phys.med(),fun_fit_phiD_s,phiD_s_chir_cont);
  plot_funz_a2("phiD_s_funz_a2.xmg",meson_name,tag_a2,meson_name,lat_med_fm,phiD_s_estr_ml,par_res_fit_phiD_s,fun_fit_phiD_s,phiD_s_chir_cont);
  
  return 0;
}
