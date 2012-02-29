#pragma once

//jack-vec version
jvec effective_mass(jvec a)
{
  int TH=a.nel-1;
  int njack=a.njack;
  
  jvec b(TH,a.njack);
  
  for(int t=0;t<TH;t++)
    {
      jack temp=-log(a[t+1]/a[t]);
      double miniz=temp.med();
      double einiz=temp.err();
      for(int ijack=0;ijack<=njack;ijack++)
	{
	  double m=miniz;
	  double e=einiz;
	  
	  double targ=a[t+1][ijack]/a[t][ijack];
	  
	  double yl;
	  double yr;
	  
	  //increment the range up to reaching opposite sign
	  int q;
	  do
	    {
	      yl=cosh((m-e)*(TH-(t+1)))/cosh((m-e)*(TH-t))-targ;
	      yr=cosh((m+e)*(TH-(t+1)))/cosh((m+e)*(TH-t))-targ;
	      q=((yl<0 && yr<0) || (yl>=0 && yr>=0));
	      if(q) e*=2;
	    }
	  while(q);
	  
	  //bisect
	  double xl=m-e,xr=m+e,ym;
	  do
	    {
	      m=(xl+xr)/2;
	      ym=cosh(m*(TH-(t+1)))/cosh(m*(TH-t))-targ;
	      if(yl<0 && ym<0 || yl>0 && ym>0)
		xl=m;
	      else
		xr=m;
	    }
	  while(fabs(ym)>1.e-13);
	  
	  b.data[t].data[ijack]=m;
	}
    }
  return b;
}

//jack-vec
jvec aperiodic_effective_mass(const jvec a)
{
  int TH=a.nel-1;
  int njack=a.njack;
  
  jvec b(TH,njack);
  
  for(int t=0;t<TH;t++) b.data[t]=log(a.data[t]/a.data[t+1]);
  
  return b;
}

//fit the mass
jack mass_fit(jvec corr,int tmin,int tmax,const char *path=NULL)
{
  jvec effe=effective_mass(corr);
  jack mass=constant_fit(effe,tmin,tmax,path);
  
  return mass;
}

//fit the mass and the matrix element
void two_pts_fit(jack &E,jack &Z,jvec corr,int tmin,int tmax,const char *path1=NULL,const char *path2=NULL)
{
  E=mass_fit(corr,tmin,tmax,path1);
  jvec temp(corr.nel,corr.njack);
  int TH=temp.nel-1;
  for(int t=0;t<=TH;t++) temp[t]=corr[t]/exp(-E*TH)/cosh(E*(TH-t))*E;
  Z=constant_fit(temp,tmin,tmax,path2);
}

//fit the mass and the matrix element in SS and SL combo
double *c_two_pts_SL_fit[2],*e_two_pts_SL_fit[2];
int TH_two_pts_SL_fit;
int tmin_two_pts_SL_fit[2];
int tmax_two_pts_SL_fit[2];

double fun_two_pts_SL_fit(double Z1,double Z2,double M,double t)
{return Z1*Z2*exp(-M*TH_two_pts_SL_fit)*cosh(M*(TH_two_pts_SL_fit-t))/M;}

void ch2_two_pts_SL_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  double M=p[0];
  double ZL=p[1];
  double ZS=p[2];
  
  for(int t=tmin_two_pts_SL_fit[0];t<=min(tmax_two_pts_SL_fit[0],TH_two_pts_SL_fit);t++)
    {
      double diff=c_two_pts_SL_fit[0][t]-fun_two_pts_SL_fit(ZL,ZS,M,t);
      double err=e_two_pts_SL_fit[0][t];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3) cout<<"SL, t="<<t<<", diff="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
  for(int t=tmin_two_pts_SL_fit[1];t<=min(tmax_two_pts_SL_fit[1],TH_two_pts_SL_fit);t++)
    {
      double diff=c_two_pts_SL_fit[1][t]-fun_two_pts_SL_fit(ZS,ZS,M,t);
      double err=e_two_pts_SL_fit[1][t];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3) cout<<"SS, t="<<t<<", diff="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void two_pts_SL_fit(jack &M,jack &ZL,jack &ZS,jvec corrSL,jvec corrSS,int tminL,int tmaxL,int tminS,int tmaxS,const char *path1=NULL,const char *path2=NULL)
{
  jack ML=constant_fit(effective_mass(corrSL),tminL,tmaxL,path1);
  jack MS=constant_fit(effective_mass(corrSS),tminS,tmaxS,path2);
  M=jack_weighted_average(ML,MS);
  jvec tempSL(corrSS.nel,corrSS.njack),tempSS(corrSS.nel,corrSS.njack);
  int TH=tempSS.nel-1;
  for(int t=0;t<=TH;t++)
    {
      tempSL[t]=corrSL[t]/exp(-M*TH)/cosh(M*(TH-t))*M;
      tempSS[t]=corrSS[t]/exp(-M*TH)/cosh(M*(TH-t))*M;
    }
  ZS=sqrt(constant_fit(tempSS,tminS,tmaxS,NULL));
  ZL=constant_fit(tempSL,tminL,tmaxL,NULL)/ZS;
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_two_pts_SL_fit);
  minu.DefineParameter(0,"M",M.med(),0.001,0,0);
  minu.DefineParameter(1,"ZL",ZL.med(),0.001,0,0);
  minu.DefineParameter(2,"ZS",ZS.med(),0.001,0,0);
  
  int njack=ML.njack;
  c_two_pts_SL_fit[0]=new double[TH+1];
  c_two_pts_SL_fit[1]=new double[TH+1];
  e_two_pts_SL_fit[0]=new double[TH+1];
  e_two_pts_SL_fit[1]=new double[TH+1];
  
  TH_two_pts_SL_fit=TH;
  tmin_two_pts_SL_fit[0]=tminL;
  tmin_two_pts_SL_fit[1]=tminS;
  tmax_two_pts_SL_fit[0]=tmaxL;
  tmax_two_pts_SL_fit[1]=tmaxS;
  
  for(int iel=0;iel<=TH;iel++)
    {
      e_two_pts_SL_fit[0][iel]=corrSL[iel].err();
      e_two_pts_SL_fit[1][iel]=corrSS[iel].err();
    }
  
  for(int ijack_fit=0;ijack_fit<=njack;ijack_fit++)
    {
      for(int iel=0;iel<=TH;iel++)
	{
	  c_two_pts_SL_fit[0][iel]=corrSL[iel][ijack_fit];
	  c_two_pts_SL_fit[1][iel]=corrSS[iel][ijack_fit];
	}
      minu.Migrad();
      double dum;
      minu.GetParameter(0,M.data[ijack_fit],dum);
      minu.GetParameter(1,ZL.data[ijack_fit],dum);
      minu.GetParameter(2,ZS.data[ijack_fit],dum);
    }
  
  double ch2,grad[3],par[3]={M[njack],ZL[njack],ZS[njack]};
  minu.Eval(3,grad,ch2,par,3);
  cout<<"ch2: "<<ch2<<endl;
  
  /*
  if(path1!=NULL)
    {
      double x[100];
      double dx=(tmaxL-tminL)/99.0;
      jvec y(100,njack);
      jvec z=corrSL;
      for(int i=0;i<100;i++)
	{
	  x[i]=tminL+dx*i;
	  for(int ijack=0;ijack<=njack;ijack++) y[i].data[ijack]=fun_two_pts_SL_fit(ZS[ijack],ZL[ijack],M[ijack],x[i]);
	  for(int ijack=0;ijack<=njack;ijack++) y[i].data[ijack]/=fun_two_pts_SL_fit(ZS[njack],ZL[njack],M[njack],x[i]);
	}
      for(int t=0;t<z.nel;t++) z[t]/=fun_two_pts_SL_fit(ZS[njack],ZL[njack],M[njack],t);
      grace out(path1);
      //out.fout<<"@    yaxes scale Logarithmic"<<endl;
      out.fout<<"@type xydy"<<endl;
      out.fout<<"@s0 line type 0"<<endl;
      out<<z<<endl;
      out.new_set();
      out.contour(x,y);
    }
  if(path2!=NULL)
    {
      double x[100];
      double dx=(tmaxS-tminS)/99.0;
      jvec y(100,njack);
      jvec z=corrSS;
      for(int i=0;i<100;i++)
	{
	  x[i]=tminS+dx*i;
	  for(int ijack=0;ijack<=njack;ijack++) y[i].data[ijack]=fun_two_pts_SL_fit(ZS[ijack],ZS[ijack],M[ijack],x[i])/fun_two_pts_SL_fit(ZS[njack],ZS[njack],M[njack],x[i]);
	}
      for(int t=0;t<z.nel;t++) z[t]/=fun_two_pts_SL_fit(ZS[njack],ZS[njack],M[njack],t);
      grace out(path2);
      //out.fout<<"@    yaxes scale Logarithmic"<<endl;
      out.fout<<"@type xydy"<<endl;
      out.fout<<"@s0 line type 0"<<endl;
      out<<z<<endl;
      out.new_set();
      out.contour(x,y);
    }
  */
}

void linear_fit(jack &m,jack &q,jvec corr,int tmin,int tmax)
{
  tmin=max(0,tmin);
  tmax=min(tmax+1,corr.nel);
  int njack=corr.njack;
  jack Y(njack),XY(njack),Y2(njack);
  double X=0,W=0,X2=0;

  Y=XY=0;
  for(int t=tmin;t<tmax;t++)
    {
      double err=corr[t].err();
      double w=1/sqr(err);
      
      W+=w;
      
      X+=t*w;
      X2+=t*t*w;
      
      Y+=corr[t]*w;
      Y2+=sqr(corr[t])*w;
      
      XY+=t*corr[t]*w;
    }

  XY-=X*Y/W;
  Y2-=Y*Y/W;
  X2-=X*X/W;
  
  m=XY/X2;
  q=(Y-m*X)/W;
}

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
