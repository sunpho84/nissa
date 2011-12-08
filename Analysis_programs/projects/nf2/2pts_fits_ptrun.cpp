#include "common.cpp"

char data_list_file[1024];
char corr_name[1024],out_file[1024];

int tmin[3],tmax[3];
int nlights,ifit_int,parity;
double *corr_fit,*corr_err;

int T,TH,L;
int iml_un,nmass,ibeta;
char base_path[1024];
double *mass;

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  read_formatted_from_file_expecting(data_list_file,fin,"%s","data_list_file");
  
  read_formatted_from_file_expecting(corr_name,fin,"%s","corr_name");
  read_formatted_from_file_expecting((char*)(&parity),fin,"%d","parity");
  
  read_formatted_from_file_expecting((char*)(&tmin[0]),fin,"%d","tint_ll");
  read_formatted_from_file((char*)(&tmax[0]),fin,"%d","tint_ll");
  
  read_formatted_from_file_expecting((char*)(&tmin[1]),fin,"%d","tint_lh");
  read_formatted_from_file((char*)(&tmax[1]),fin,"%d","tint_lh");
  
  read_formatted_from_file_expecting((char*)(&tmin[2]),fin,"%d","tint_hh");
  read_formatted_from_file((char*)(&tmax[2]),fin,"%d","tint_hh");
  
  read_formatted_from_file_expecting(out_file,fin,"%s","out_file");
  
  fclose(fin);
}

//function to fit
double fun_fit(double Z2,double M,int t)
{
  if(parity==1) return Z2*exp(-M*TH)*cosh(M*(TH-t))/sinh(M);
  else return Z2*exp(-M*TH)*sin(M*(TH-t))/sinh(M);
}

//chi2 calculation
void chi2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  for(int t=tmin[ifit_int];t<=min(tmax[ifit_int],TH);t++)
    ch+=sqr((corr_fit[t]-fun_fit(p[0],p[1],t))/corr_err[t]);
}

int main()
{
  read_pars("input");
  read_ensemble_pars(base_path,T,ibeta,nmass,mass,iml_un,nlights,data_list_file);
  TH=L=T/2;
  
  init_latpars();
  
  int ncombo=nmass*nmass;
  
  //load all the corrs
  double *buf=new double[4*ncombo*T*(njack+1)];
  FILE *fin=open_file(combine("%s/%s",base_path,corr_name).c_str(),"r");
  int stat=fread(buf,sizeof(double),4*ncombo*(njack+1)*T,fin);
  if(stat!=4*ncombo*(njack+1)*T)
    {
      cerr<<"Error loading data!"<<endl;
      exit(1);
    }
  
  jvec M(ncombo,njack);
  jvec Z2(ncombo,njack);
  
  //define minuit staff
  TMinuit minu(2);
  minu.SetPrintLevel(-1);
  minu.SetFCN(chi2);
  corr_fit=new double[TH+1];
  corr_err=new double[TH+1];
  
  //fit each combo
  int ic=0;
  for(int ims=0;ims<nmass;ims++)
    for(int imc=ims;imc<nmass;imc++)
      {
	//take into account corr
	jvec corr1(T,njack),corr2(T,njack),corr;
	int ic1=0+2*(ims+nmass*(0+2*imc));
	int ic2=1+2*(ims+nmass*(1+2*imc));
	cout<<ims<<" "<<imc<<" "<<ic1<<" "<<ic2<<endl;
	corr1.put(buf+ic1*T*(njack+1));
	corr2.put(buf+ic2*T*(njack+1));
	corr=(corr1+corr2)/2;
	cout<<corr1[0]<<" "<<corr2[0]<<endl;
	
	//choose the index of the fitting interval
	if(ims>=nlights) ifit_int=2;
	else
	  if(imc>=nlights) ifit_int=1;
	  else ifit_int=0;
	
	//simmetrize
	corr=corr.simmetrized(parity);
	int ttmin=tmin[ifit_int];
	int ttmax=tmax[ifit_int];
	jvec Mcor=effective_mass(corr),Z2cor(TH+1,njack);
	jack Meff=constant_fit(Mcor,ttmin,ttmax);
	for(int t=0;t<=TH;t++)
	  for(int ijack=0;ijack<=njack;ijack++)
	    Z2cor[t].data[ijack]=corr[t].data[ijack]/fun_fit(1,Meff[ijack],t);
	jack Z2eff=constant_fit(Z2cor,ttmin,ttmax);
	
	if(!isnan(Z2eff[0])) minu.DefineParameter(0,"Z2",Z2eff[0],Z2eff.err(),0,2*Z2eff[0]);
	if(!isnan(Meff[0])) minu.DefineParameter(1,"M",Meff[0],Meff.err(),0,2*Meff[0]);
	for(int t=tmin[ifit_int];t<=tmax[ifit_int];t++) corr_err[t]=corr.data[t].err();
	
	//jacknife analysis
	for(int ijack=0;ijack<njack+1;ijack++)
	  {
	    //copy data so that glob function may access it
	    for(int t=tmin[ifit_int];t<=tmax[ifit_int];t++) corr_fit[t]=corr.data[t].data[ijack];
	  
	    //fit
	    double dum;
	    minu.Migrad();	    
	    minu.GetParameter(0,Z2.data[ic].data[ijack],dum);
	    minu.GetParameter(1,M.data[ic].data[ijack],dum);
	  }
	
	//if((ims==iml_un||ims==nlights-1||ims==nlights||ims==nmass-1)&&
	//(imc==iml_un||imc==nlights-1||imc==nlights||imc==nmass-1))
	  {
	    //plot eff mass
	    {
	      ofstream out(combine("eff_mass_plot_%02d_%02d.xmg",ims,imc).c_str());
	      out<<"@type xydy"<<endl;
	      out<<"@s0 line type 0"<<endl;
	      out<<Mcor<<endl;
	      out<<"&"<<endl;
	      out<<"@type xy"<<endl;
	      double av_mass=M[ic].med();
	      double er_mass=M[ic].err();
	      out<<tmin[ifit_int]<<" "<<av_mass-er_mass<<endl;
	      out<<tmax[ifit_int]<<" "<<av_mass-er_mass<<endl;
	      out<<tmax[ifit_int]<<" "<<av_mass+er_mass<<endl;
	      out<<tmin[ifit_int]<<" "<<av_mass+er_mass<<endl;
	      out<<tmin[ifit_int]<<" "<<av_mass-er_mass<<endl;
	    }
	    //plot fun
	    {
	      ofstream out(combine("fun_plot_%02d_%02d.xmg",ims,imc).c_str());
	      out<<"@type xydy"<<endl;
	      out<<"@s0 line type 0"<<endl;
	      out<<corr<<endl;
	      out<<"&"<<endl;
	      out<<"@type xy"<<endl;
	      for(int t=tmin[ifit_int];t<tmax[ifit_int];t++)
		out<<t<<" "<<fun_fit(Z2[ic][njack],M[ic][njack],t)<<endl;
	    }
	  }
	
	  cout<<mass[ims]<<" "<<mass[imc]<<"  "<<M[ic]<<" "<<Z2[ic]<<" f:"<<sqrt(Z2[ic])/(sinh(M[ic])*M[ic])*(mass[ims]+mass[imc])<<" fV:"<<Za_med[ibeta]*sqrt(Z2[ic])/M[ic]/lat[ibeta].med()<<endl;
	  ic++;
      }
  
  ofstream out("fitted_mass.xmg");
  out<<"@type xydy"<<endl;
  for(int ims=0;ims<nmass;ims++)
    {
      //out<<"s0 line type 0"<<endl;
      for(int imc=0;imc<nmass;imc++) out<<mass[imc]<<" "<<M[icombo(ims,imc,nmass,nlights,0)]<<endl;
      out<<"&"<<endl;
    }

  M.write_to_binfile(out_file);
  Z2.append_to_binfile(out_file);
  
  return 0;
}
