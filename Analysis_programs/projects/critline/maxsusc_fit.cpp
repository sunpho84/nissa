#include "include.h"

int nx,nx_rew,nazione,noss,npow,nboot,njack=2;
double *x,*y;
double *x_rew,*y_rew;
const int vol=16384;

int iof(int is,int io,int ipow){return ipow+npow*(io+noss*is);}

bvec lin_solve(double *A,bvec b)
{
  int d=b.nel;
  int nboot=b.nboot;
  int njack=b.njack;

  bvec x(d,nboot,njack);
  
  for(int i=0;i<d;i++)
    {
      double C=A[i*d+i];
      for(int j=i;j<d;j++) A[i*d+j]/=C;
      b[i]/=C;
      
      for(int k=i+1;k<d;k++)
	{
	  double C=A[k*d+i];
	  for(int j=i;j<d;j++) A[k*d+j]-=A[i*d+j]*C;
	  b[k]-=C*b[i];
	}
    }
  
  for(int k=d-1;k>=0;k--)
    {
      boot S(nboot,njack);
      S=0;
      for(int i=k+1;i<d;i++) S+=A[k*d+i]*x[i];
      x[k]=b[k]-S;
    }
  
  return x;
}

bvec poly_fit(double *x,bvec y,int d)
{
  int np=y.nel;
  int nboot=y.nboot;
  int njack=y.njack;
  
  double Al[2*d+1];memset(Al,0,sizeof(double)*(2*d+1));
  bvec c(d+1,nboot,njack);c=0;

  for(int p=0;p<np;p++)
    {
      //calculate the weight
      double w=pow(y[p].err(),-2);
      //compute Al and c
      for(int f=0;f<=2*d;f++)
	{
	  Al[f]+=w;
	  if(f<=d) c[f]+=y[p]*w;
	  w*=x[p];
	}
    }
  
  double A[(d+1)*(d+1)];
  for(int i=0;i<=d;i++)
    for(int j=0;j<=d;j++)
      A[i*(d+1)+j]=Al[i+j];
  
  return lin_solve(A,c);
}      

boot pol(bvec p,double x)
{
  boot out=p[0];
  double xx=x;
  for(int i=1;i<p.nel;i++)
    {
      out+=xx*p[i];
      xx*=x;
    }
  
  return out;
}

bvec lorentzian_fit(double *x,bvec y)
{
  bvec d=poly_fit(x,1/y,2);
  bvec out(3,y.nboot,y.njack);
  //for(int i=0;i<3;i++) cout<<i<<" = "<<d[i]<<endl;
  out[1]=-0.5*d[1]/d[2];            //x0
  out[0]=1/(d[0]-sqr(out[1])*d[2]); //N
  out[2]=1/sqrt(d[2]*out[0]);       //sigma
  //cout<<"arg:"<<d[2]*out[0]<<endl;
  //for(int i=0;i<3;i++) cout<<i<<" out= "<<out[i]<<endl;

  return out;
}

boot find_xmax(double *x,bvec y)
{
  boot xmax(y.nboot,y.njack);
  
  for(int iboot=0;iboot<=y.nboot;iboot++)
    {
      int imax=0;
      for(int iel=0;iel<y.nel;iel++)
	if(y[iel][iboot]>=y[imax][iboot]) imax=iel;
      xmax.data[iboot]=x[imax];
    }
  
  return xmax;
}

boot get_xmax(const char *path,const int det)
{
  FILE *f=open_file(combine("%s/data.raw",path).c_str(),"r");
  fread(&nx,sizeof(int),1,f);
  fread(&nx_rew,sizeof(int),1,f);
  fread(&nazione,sizeof(int),1,f);
  fread(&noss,sizeof(int),1,f);
  fread(&npow,sizeof(int),1,f);
  fread(&nboot,sizeof(int),1,f);
  nboot-=1;
  
  x=new double[nazione*nx];
  y=new double[nazione*nx*(nboot+1)*noss*npow];
  
  x_rew=new double[nazione*nx_rew];
  y_rew=new double[nazione*nx_rew*(nboot+1)*noss*npow];
  
  fread(x,sizeof(double),nazione*nx,f);
  fread(y,sizeof(double),nx*(nboot+1)*noss*npow,f);
  fread(x_rew,sizeof(double),nazione*nx_rew,f);
  fread(y_rew,sizeof(double),nx_rew*(nboot+1)*noss*npow,f);

  fclose(f);
  
  bvec da(nx,nboot,njack);
  bvec da_rew(nx_rew,nboot,njack);
  int ipowo=1,io=3;
  for(int is=0;is<nx;is++) da[is].put(y+(nboot+1)*iof(is,io,ipowo));
  for(int ip=0;ip<nx_rew;ip++) da_rew[ip].put(y_rew+(nboot+1)*iof(ip,io,ipowo));

  da*=vol;
  da_rew*=vol;
  
  bvec par_rew=lorentzian_fit(x_rew,da_rew);
  //for(int i=0;i<3;i++) cout<<par[i]<<" "<<par[i].err()/par[i].med()<<endl;
  //for(int i=0;i<3;i++) cout<<par_rew[i]<<" "<<par_rew[i].err()/par_rew[i].med()<<endl;
  
  //filter
  int nrew_filt=0,nfilt=0;
  double tresh=par_rew[0].med()/4;
  double tresh_rew=par_rew[0].med()*3/4;
  for(int i=0;i<da_rew.nel;i++) nrew_filt+=da_rew[i].med() > tresh_rew ? 1 : 0;
  for(int i=0;i<da.nel;i++) nfilt+=da[i].med() > tresh ? 1 : 0;
  double x_f[nfilt];
  double x_rew_f[nrew_filt];
  bvec da_f(nfilt,nboot,njack);
  bvec da_rew_f(nrew_filt,nboot,njack);
  int ifilt=0;
  int irew_filt=0;
  for(int i=0;i<da_rew.nel;i++)
    if(da_rew[i].med()>tresh_rew)
      {
	x_rew_f[irew_filt]=x_rew[i];
	da_rew_f[irew_filt]=da_rew[i];
	//cout<<ifilt<<" "<<x_rew_f[irew_filt]<<" "<<1/da_rew_f[irew_filt]<<endl;
	irew_filt++;
      }
  for(int i=0;i<da.nel;i++)
    if(da[i].med()>tresh)
      {
	x_f[ifilt]=x[i];
	da_f[ifilt]=da[i];
	//cout<<ifilt<<" "<<x_f[ifilt]<<" "<<1/da_f[ifilt]<<endl;
	ifilt++;
      }
  
  //cout<<"Det: "<<det<<endl;
  //cout<<"NFilt: "<<nfilt<<endl;
  
  bvec par_f=lorentzian_fit(x_f,da_f);
  bvec par_rew_f=lorentzian_fit(x_rew_f,da_rew_f);
  for(int i=0;i<3;i++)
    for(int iboot=0;iboot<nboot;iboot++)
      {
	if(isnan(par_f[i][iboot]))
	  {
	    par_f[i].data[iboot]=par_f[i].data[(iboot-1+nboot)%nboot];
	    //cout<<"par: "<<par_f[i][iboot]<<" "<<par_f[i].err()/par_f[i].med()<<endl;
	  }
	if(isnan(par_rew_f[i][iboot]))
	  {
	    par_rew_f[i].data[iboot]=par_rew_f[i].data[(iboot-1+nboot)%nboot];
	    //cout<<"par_: "<<par_rew_f[i][iboot]<<" "<<par_rew_f[i].err()/par_rew_f[i].med()<<endl;
	  }
      }
  bvec re_f(da_rew_f.nel,nboot,njack);
  bvec re_rew_f(da_rew_f.nel,nboot,njack);
  for(int i=0;i<da_rew_f.nel;i++) re_f[i]=par_f[0]/(1+sqr((x_rew_f[i]-par_f[1])/par_f[2]));
  for(int i=0;i<da_rew_f.nel;i++) re_rew_f[i]=par_rew_f[0]/(1+sqr((x_rew_f[i]-par_rew_f[1])/par_rew_f[2]));

  grace out("%s/fit_susc_%d.xmg",path,det);
  out.set(2,"none","square");
  if(det==0) out.print_graph(x_rew_f,da_rew_f);
  if(det==2) out.print_graph(x_f,da_f);
  //out.print_graph(x_rew,da_rew);
  out.new_set();
  out.set(1,"green");
  if(det==0) out.polygon(x_rew_f,da_rew_f);
  if(det==2) out.polygon(x_f,da_f);
  out.new_set();
  out.set(1,"blue");
  if(det==0) out.polygon(x_rew_f,re_rew_f);
  if(det==2) out.polygon(x_rew_f,re_f);
  
  boot xmax_rew=find_xmax(x_rew_f,da_rew_f);
  //cout<<xmax<<endl;
  
  switch(det)
    {
    case 0:return par_rew_f[1];break;
    case 1:return xmax_rew;break;
    case 2:return par_f[1];break;
    }
  
  cerr<<"Unknown value of 'det'"<<endl;
  exit(0);
  return par_f[1];
}

double chi2_pol(bvec par,double *x,bvec y)
{
  double c2=0;
  for(int i=0;i<y.nel;i++) c2+=sqr((y[i]-pol(par,x[i])).med()/y[i].err());
  
  return c2;
}

int nxfit;
double *xfit,*erry,*yfit;
double ratfun21(double x,double *p){return (p[0]+x*p[1]+x*x*p[2])/(1+p[3]*x);}
//double ratfun21(double x,double *p){return (p[0]+x*p[1]+x*x*p[2])/(1-(0.41-p[1])/p[0]*x);}
double chi2_rat21(double *p)
{
  double ch2=0;
  for(int i=0;i<nxfit;i++)
    ch2+=sqr((yfit[i]-ratfun21(xfit[i],p))/erry[i]);
  return ch2;
}

void ch2_wrap(int &npar,double *fuf,double &ch,double *p,int flag)
{ch=chi2_rat21(p);}

bvec rat21_fit(bvec par_fit4,double *muI2,bvec xmax)
{
  int nmu=xmax.nel;
  TMinuit minu(4);
  minu.SetFCN(ch2_wrap);
  minu.DefineParameter(0,"A0",par_fit4[0].med(),par_fit4[0].err(),0,0);
  minu.DefineParameter(1,"A1",par_fit4[1].med(),par_fit4[1].err(),0,0);
  minu.DefineParameter(2,"A2",par_fit4[2].med(),par_fit4[2].err(),0,0);
  minu.DefineParameter(3,"A3",0,0.001,0,0);
  //minu.FixParameter(3);
  nxfit=nmu;
  xfit=muI2;
  yfit=new double[nmu];
  erry=new double[nmu];
  for(int ifit=0;ifit<nmu;ifit++) erry[ifit]=xmax[ifit].err();
  bvec par_fitR(4,nboot,njack);
  for(int iboot=0;iboot<=nboot;iboot++)
    {
      for(int ifit=0;ifit<nmu;ifit++) yfit[ifit]=xmax[ifit][iboot];
      minu.Migrad();
      double dum;
      for(int ipar=0;ipar<4;ipar++) minu.GetParameter(ipar,par_fitR[ipar].data[iboot],dum);
    }
  
  return par_fitR;
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      fprintf(stderr,"Error, use: %s file\n",arg[0]);
      exit(1);
    }
  
  grace out("fit_xmax.xmg");
  FILE *fin=open_file(arg[1],"r");
  const char color[3][1024]={"blue","green","red"};
  int nmu;
  read_formatted_from_file_expecting((char*)&nmu,fin,"%d","nmu");

  double muI[nmu],muI2[nmu];
  bvec xmax[3];xmax[0]=xmax[1]=xmax[2]=bvec(nmu,99,2);
  for(int imu=0;imu<nmu;imu++)
    {
      char path[1024];
      read_formatted_from_file((char*)&(muI[imu]),fin,"%lg","path");
      read_formatted_from_file(path,fin,"%s","path");
      for(int det=0;det<3;det++) xmax[det][imu]=get_xmax(path,det);
      muI2[imu]=sqr(muI[imu]);
      if(muI[imu]<0) muI2[imu]*=-1;
    }
  
  for(int det=0;det<3;det++)
    {
      bvec par_fit4=poly_fit(muI2,xmax[det],2);
      bvec par_fit6=poly_fit(muI2,xmax[det],3);
      bvec par_fitR=rat21_fit(par_fit4,muI2,xmax[det]);
      
      cout<<chi2_pol(par_fit4,muI2,xmax[det])/(nmu-3)<<endl;
      cout<<chi2_pol(par_fit6,muI2,xmax[det])/(nmu-4)<<endl;
      {
	double p[4]={par_fitR[0][nboot],par_fitR[1][nboot],par_fitR[2][nboot],par_fitR[3][nboot]};
	cout<<chi2_rat21(p)/(nmu-4)<<endl;
      }
      
      int nfit=300;
      double mu2F[nfit+1],mu2F_min=-0.5,mu2F_max=muI2[nmu-1];
      
      bvec xmaxF4(nfit+1,99,2);
      bvec xmaxF6(nfit+1,99,2);
      bvec xmaxFR(nfit+1,99,2);
      for(int ifit=0;ifit<=nfit;ifit++)
	{
	  double x=mu2F[ifit]=mu2F_min+(mu2F_max-mu2F_min)/nfit*ifit;
	  xmaxF4[ifit]=pol(par_fit4,x);
	  xmaxF6[ifit]=pol(par_fit6,x);
	  for(int iboot=0;iboot<=nboot;iboot++)
	    {
	      double p[4]={par_fitR[0][iboot],par_fitR[1][iboot],par_fitR[2][iboot],par_fitR[3][iboot]};
	      xmaxFR[ifit].data[iboot]=ratfun21(x,p);
	    }
	}
      
      out.set(2,"none","square");
      out.print_graph(muI2,xmax[det]);
      out.new_set();
      out.set(1,color[det]);
      out.polygon(mu2F,xmaxF4);
      out.new_set();
      out.set(1,color[det]);
      out.polygon(mu2F,xmaxF6);
      out.new_set();
      out.set(1,color[det]);
      out.polygon(mu2F,xmaxFR);
      out.new_set();
      
      cout<<par_fit4[1]<<endl;
      cout<<par_fit6[1]<<endl;
      cout<<par_fitR[1]-par_fitR[3]*par_fitR[0]<<endl;
    }
  
  return 0;
}
