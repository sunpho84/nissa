#include "common.cpp"

char data_list_file[1024];
char corr_name[1024],out_file[1024];

int tmin,tmax;
int nlights,parity;
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
  
  read_formatted_from_file_expecting((char*)(&tmin),fin,"%d","tint");
  read_formatted_from_file((char*)(&tmax),fin,"%d","tint");
  
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
  for(int t=tmin;t<=min(tmax,TH);t++)
    ch+=sqr((corr_fit[t]-fun_fit(p[0],p[1],t))/corr_err[t]);
}

int main()
{
  read_pars("input");
  read_ensemble_pars(base_path,T,ibeta,nmass,mass,iml_un,nlights,data_list_file);
  TH=L=T/2;
  
  int ncombo=nlights*(nmass-(nlights-1)/2);
  
  //load all the corrs
  double *buf=new double[ncombo*T*(njack+1)];
  FILE *fin=open_file(combine("%s/%s",base_path,corr_name).c_str(),"r");
  int stat=fread(buf,sizeof(double),ncombo*(njack+1)*T,fin);
  if(stat!=ncombo*(njack+1)*T)
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
  
  int ic=0;
  //fit each combo
  for(int ims=0;ims<nlights;ims++)
    for(int imc=ims;imc<nmass;imc++)
      {
	cout<<ims<<" "<<imc<<endl;
	//take into account corr
	jvec corr(T,njack);
	corr.put(buf+ic*T*(njack+1));
	if((string)corr_name=="VKVK") corr*=-1;
	
	//simmetrize
	corr=corr.simmetrized(parity);
	jvec Mcor=effective_mass(corr),Z2cor(TH+1,njack);
	jack Meff=constant_fit(Mcor,tmin,tmax);
	for(int t=0;t<=TH;t++)
	  for(int ijack=0;ijack<=njack;ijack++)
	    Z2cor[t].data[ijack]=corr[t].data[ijack]/fun_fit(1,Meff[ijack],t);
	jack Z2eff=constant_fit(Z2cor,tmin,tmax);
	
	if(!isnan(Z2eff[0])) minu.DefineParameter(0,"Z2",Z2eff[0],Z2eff.err(),0,2*Z2eff[0]);
	if(!isnan(Meff[0])) minu.DefineParameter(1,"M",Meff[0],Meff.err(),0,2*Meff[0]);
	for(int t=tmin;t<=tmax;t++) corr_err[t]=corr.data[t].err();
	
	//jacknife analysis
	for(int ijack=0;ijack<njack+1;ijack++)
	  {
	    //copy data so that glob function may access it
	    for(int t=tmin;t<=tmax;t++) corr_fit[t]=corr.data[t].data[ijack];
	    
	    //fit
	    double dum;
	    minu.Migrad();	    
	    minu.GetParameter(0,Z2.data[ic].data[ijack],dum);
	    minu.GetParameter(1,M.data[ic].data[ijack],dum);
	  }
	
	//plot eff mass
	{
	  ofstream out(combine("eff_mass_plot_%02d_%02d.xmg",ims,imc).c_str());
	  out<<"@type xydy"<<endl;
	  out<<"@s0 line type 0"<<endl;
	  out<<Mcor<<endl;
	  out<<"&"<<endl;
	  out<<"@type xy"<<endl;
	  double av_mass=Meff.med();
	  double er_mass=Meff.err();
	  out<<tmin<<" "<<av_mass-er_mass<<endl;
	  out<<tmax<<" "<<av_mass-er_mass<<endl;
	  out<<tmax<<" "<<av_mass+er_mass<<endl;
	  out<<tmin<<" "<<av_mass+er_mass<<endl;
	  out<<tmin<<" "<<av_mass-er_mass<<endl;
	}
	//plot fun
	{
	  ofstream out(combine("fun_plot_%02d_%02d.xmg",ims,imc).c_str());
	  out<<"@yaxes scale Logarithmic"<<endl;
	  out<<"@type xydy"<<endl;
	  out<<"@s0 line type 0"<<endl;
	  out<<corr<<endl;
	  out<<"&"<<endl;
	  out<<"@type xy"<<endl;
	  jvec temp(tmax-tmin+1,njack);
	  for(int t=tmin;t<=tmax;t++)
	    for(int ijack=0;ijack<=njack;ijack++) temp[t-tmin].data[ijack]=fun_fit(Z2[ic][ijack],M[ic][ijack],t);
	  for(int t=tmin;t<=tmax;t++)
	    out<<t<<" "<<temp[t-tmin].med()-temp[t-tmin].err()<<endl;
	  for(int t=tmax;t>=tmin;t--)
	    out<<t<<" "<<temp[t-tmin].med()+temp[t-tmin].err()<<endl;
	}
	
	cout<<mass[ims]<<" "<<mass[imc]<<"  "<<M[ic]<<" ("<<Meff<<") "<<Z2[ic]<<" "<<endl;
	ic++;
      }
  
  ic=0;
  ofstream out("fitted_mass.xmg");
  out<<"@type xydy"<<endl;
  for(int ims=0;ims<nlights;ims++)
    {
      //out<<"s0 line type 0"<<endl;
      for(int imc=ims;imc<nmass;imc++)
	{
	  out<<mass[imc]<<" "<<M[ic]<<endl;
	  ic++;
	}
      out<<"&"<<endl;
    }

  M.write_to_binfile(out_file);
  Z2.append_to_binfile(out_file);
  
  return 0;
}
