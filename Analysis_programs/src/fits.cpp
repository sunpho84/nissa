#pragma once

//return effective mass
double effective_mass(double y,double yp,int t,int TH)
{
  double rapp=yp/y;
  double meff=-log(rapp);

  double res;
  int iter=0;
  do
    {
      double rapr=cosh(meff*(TH-(t+1)))/cosh(meff*(TH-t));
      double scale=rapr/rapp;
      res=fabs(1-scale)/meff;
      meff*=scale;
      iter++;
    }
  while(res>1.e-8 && iter<1000);
  
  return meff;
}

//jack version of effective mass
jack effective_mass(const jack y,const jack yp,int t,int TH)
{
  int njack=y.njack;
  jack m(njack);
  
  for(int ij=0;ij<njack+1;ij++) m.data[ij]=effective_mass(y.data[ij],yp.data[ij],t,TH);
  
  return m;
}

//jack-vec version
jvec effective_mass(const jvec a)
{
  int TH=a.nel-1;
  int njack=a.njack;
  
  jvec b(TH,njack);
  
  for(int t=0;t<TH;t++) b.data[t]=effective_mass(a.data[t],a.data[t+1],t,TH);
  
  return b;
}

//fit the mass
jack mass_fit(jvec corr,int tmin,int tmax,const char *path=NULL)
{
  jvec effe=effective_mass(corr);
  jack mass=constant_fit(effe,tmin+1,tmax,path);
  
  return mass;
}

//fit the mass and the matrix element
void P5P5_fit(jack &E,jack &Z,jvec corr,int tmin,int tmax,const char *path1=NULL,const char *path2=NULL)
{
  E=mass_fit(corr,tmin,tmax,path1);
  jvec temp(corr.nel,corr.njack);
  int TH=temp.nel-1;
  for(int t=0;t<=TH;t++) temp[t]=corr[t]/exp(-E*TH)/cosh(E*(TH-t))*E;
  Z=constant_fit(temp,tmin,tmax,path2);
}
