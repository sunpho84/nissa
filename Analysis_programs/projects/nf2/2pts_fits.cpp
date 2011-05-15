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
  
  read_formatted_from_file_expecting((char*)(&nlights),fin,"%d","nlights");
  
  read_formatted_from_file_expecting(out_file,fin,"%s","out_file");
  
  fclose(fin);
}

//function to fit
double fun_fit(double Z2,double M,int t)
{
  if(parity==1) return Z2*exp(-M*TH)*cosh(M*(TH-t))/M;
  else return Z2*exp(-M*TH)*sin(M*(TH-t))/M;
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
  read_setup_pars(base_path,T,ibeta,nmass,iml_un,mass,data_list_file);
  TH=L=T/2;
  
  init_latpars(100);
  
  int ncombo=nmass*(nmass+1)/2;
  
  //load all the corrs
  double *buf=new double[ncombo*T*(njack+1)];
  FILE *fin=open_file(combine("%s/%s",base_path,"P5P5").c_str(),"r");
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
  
  //fit each combo
  int icombo=0;
  for(int ims=0;ims<nmass;ims++)
    for(int imc=ims;imc<nmass;imc++)
      {
	//take into account corr
	jvec corr(T,njack);
	corr.put(buf+icombo*T*(njack+1));
	
	//choose the index of the fitting interval
	if(ims>=nlights) ifit_int=2;
	else
	  if(imc>=nlights) ifit_int=1;
	  else ifit_int=0;
	
	//simmetrize
	corr=corr.simmetrized(parity);
	jvec mcor=effective_mass(corr);
	jack meff=constant_fit(mcor,tmin[ifit_int]+1,tmax[ifit_int]);
	jvec temp(corr.nel,corr.njack);
	for(int t=tmin[ifit_int];t<=tmax[ifit_int];t++)
	  for(int ijack=0;ijack<njack+1;ijack++)
	    temp[t].data[ijack]=corr[t].data[ijack]/fun_fit(1,meff.data[ijack],t);
	jack cestim=constant_fit(temp,tmin[ifit_int],tmax[ifit_int]);
	
	//jacknife analysis
	for(int ijack=0;ijack<njack+1;ijack++)
	  {
	    //copy data so that glob function may access it
	    for(int t=tmin[ifit_int];t<=tmax[ifit_int];t++)
	      {
		corr_fit[t]=corr.data[t].data[ijack];	
		corr_err[t]=corr.data[t].err();
		//cout<<t<<" "<<corr_fit[t]<<" "<<corr_err[t]<<endl;
	      }
	  
	    //fit
	    double dum;
	    minu.DefineParameter(0,"Z2",meff.data[ijack],0.001,0,0);
	    minu.DefineParameter(1,"M",cestim.data[ijack],0.001,0,0);
	    minu.Migrad();
	    minu.GetParameter(0,Z2.data[icombo].data[ijack],dum);
	    minu.GetParameter(1,M.data[icombo].data[ijack],dum);
	    //cout<<ims<<" "<<imc<<"  "<<mass[ims]<<" "<<mass[imc]<<"  "<<icombo<<"  "<<Z2.data[icombo].data[ijack]<<" "<<M.data[icombo].data[ijack]<<endl;
	  }
	
	if((ims==iml_un||ims==nlights-1||ims==nlights||ims==nmass-1)&&
	   (imc==iml_un||imc==nlights-1||imc==nlights||imc==nmass-1))
	  {
	    
	    //plot fit
	    ofstream out(combine("fit_plot_%02d_%02d.xmg",ims,imc).c_str());
	    out<<"@type xydy"<<endl;
	    out<<"@s0 line type 0"<<endl;
	    out<<mcor<<endl;
	    out<<"&"<<endl;
	    out<<"@type xy"<<endl;
	    double av_mass=M[icombo].med();
	    double er_mass=M[icombo].err();
	    out<<tmin[ifit_int]+1<<" "<<av_mass-er_mass<<endl;
	    out<<tmax[ifit_int]+1<<" "<<av_mass-er_mass<<endl;
	    out<<tmax[ifit_int]+1<<" "<<av_mass+er_mass<<endl;
	    out<<tmin[ifit_int]+1<<" "<<av_mass+er_mass<<endl;
	    out<<tmin[ifit_int]+1<<" "<<av_mass-er_mass<<endl;
	  }
	
	cout<<ims<<" "<<imc<<"  "<<mass[ims]<<" "<<mass[imc]<<"  "<<icombo<<"  "<<M[icombo]<<"  "<<Z2[icombo]<<endl;
	
	icombo++;
      }
  
  M.write_to_binfile(out_file);
  Z2.append_to_binfile(out_file);
  
  return 0;
}
