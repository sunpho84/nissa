#include "common.cpp"

char data_list_file[1024];
int tmin_l,tmax_l;
int tmin_h,tmax_h;

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  read_formatted_from_file_expecting(data_list_file,fin,"%s","data_list_file");
  read_formatted_from_file_expecting((char*)(&tmin_l),fin,"%d","tmin_l");
  read_formatted_from_file_expecting((char*)(&tmax_l),fin,"%d","tmax_l");
  read_formatted_from_file_expecting((char*)(&tmin_h),fin,"%d","tmin_h");
  read_formatted_from_file_expecting((char*)(&tmax_h),fin,"%d","tmax_h");
  fclose(fin);
}

int icombo(int im1,int im2)
{
  int ims=min(im1,im2);
  int imc=max(im1,im2);
  
  return ims*nmass-ims*(ims-1)/2+(imc-ims);
}

double det3(double *a,double *b,double *c)
{
  double d;
  d= a[0]*(b[1]*c[2]-b[2]*c[1]);
  d+=b[0]*(c[1]*a[2]-c[2]*a[1]);
  d+=c[0]*(a[1]*b[2]-a[2]*b[1]);
  
  return d;
}
void parabol_int(double &a,double &b,double &c,double *xin,double *yd)
{
  double ya[3],yb[3],yc[3];
  for(int i=0;i<3;i++)
    {
      ya[i]=sqr(xin[i]);
      yb[i]=xin[i];
      yc[i]=1;
    }
  
  double D=det3(ya,yb,yc);
  double Da=det3(yd,yb,yc);
  double Db=det3(ya,yd,yc);
  double Dc=det3(ya,yb,yd);
  
  a=Da/D;
  b=Db/D;
  c=Dc/D;
}

boot interpolate_single(bvec vec,double *x,boot xf)
{
  int nx=vec.nel;
  boot out(nboot,njack);

  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      //check that the needed x is not below the lower one or above the heavier one
      double dmin=(x[1]-x[0])/2;
      double dmax=(x[nx-1]-x[nx-2])/2;
      bool linbelow=(xf[iboot]-x[0]<dmin);
      bool linabove=(x[nx-1]-xf[iboot]<dmax);
      if(nx==2) linabove=1;
      
      //if linear
      if(linabove||linbelow)
	{
	  int sup;
	  if(linabove) sup=nx-1;
	  else sup=1;
	  
	  double m=(vec[sup][iboot]-vec[sup-1][iboot])/(x[sup]-x[sup-1]);
	  double q=vec[sup-1][iboot]-m*x[sup-1];
	  
	  out.data[iboot]=m*xf[iboot]+q;
	  
	  if(iboot==nboot)
	    {
	      cout<<" "<<sup-1<<" "<<x[sup-1]<<" "<<vec[sup-1][iboot]<<" "<<m*x[sup-1]+q<<endl;
	      cout<<" "<<sup<<" "<<x[sup]<<" "<<vec[sup][iboot]<<" "<<m*x[sup]+q<<endl;
	      cout<<"  int: "<<xf[iboot]<<" "<<out.data[iboot]<<endl;
	    }	  
	}
      else
	{
	  //if parabolic interpolation, find the nearest point
	  int nearx=0;
	  for(int ix=0;ix<nx;ix++) if(fabs(x[ix]-xf[iboot])<fabs(x[nearx]-xf[iboot])) nearx=ix;
	  
	  //copy the x and y for current bootstrap
	  double xin[3]={x[nearx-1],x[nearx],x[nearx+1]};
	  double yin[3]={vec[nearx-1][iboot],vec[nearx][iboot],vec[nearx+1][iboot]};
	  
	  //find the spline
	  double a,b,c;
	  parabol_int(a,b,c,xin,yin);
	  
	  //interpolate
	  out.data[iboot]=a*sqr(xf[iboot])+b*xf[iboot]+c;

	  if(iboot==nboot)
	    {
	      cout<<" "<<nearx-1<<" "<<x[nearx-1]<<" "<<vec[nearx-1][iboot]<<endl;
	      cout<<" "<<nearx<<" "<<x[nearx]<<" "<<vec[nearx][iboot]<<endl;
	      cout<<" "<<nearx+1<<" "<<x[nearx+1]<<" "<<vec[nearx+1][iboot]<<endl;
	      cout<<"  int: "<<xf[iboot]<<" "<<out.data[iboot]<<endl;
	    }	  
	}
    }
  
  return out;
}

boot interpolate_double(bvec vec)
{
  bvec strange_inted(nstrange,nboot,njack);
  boot out(nboot,njack);

  //first of all interpolate in the charm
  for(int _ims=0;_ims<nstrange;_ims++)
    {
      int ims=_ims+1;
      
      //interpolate around the charm
      int ncharm=nmass-nstrange-1;
      bvec temp_vec(ncharm,nboot,njack);
      double x[ncharm];
      cout<<" istrange: "<<ims<<endl;
      for(int _imc=0;_imc<ncharm;_imc++)
	{
	  int imc=_imc+nstrange+1;
	  x[_imc]=mass[imc];
	  temp_vec[_imc]=vec[icombo(ims,imc)];
	  cout<<" "<<ims<<" "<<imc<<" "<<icombo(ims,imc)<<" "<<temp_vec[_imc]<<endl;
	}
      
      strange_inted[_ims]=interpolate_single(temp_vec,x,amc_phys[ibeta]);
      cout<<"mc_int_c "<<mass[ims]<<" "<<strange_inted[_ims]<<endl;
    }
  
  //now interpolate around strange
  double x[nstrange];
  
  for(int _ims=0;_ims<nstrange;_ims++) x[_ims]=mass[_ims+1]; //take away l
  out=interpolate_single(strange_inted,x,ams_phys[ibeta]);
  
  cout<<"mc_int_s "<<out<<endl;
  
  return out;
}
  
void fit(jack &mfit,jack &ffit,jvec corr,const char *out_path,int tmin,int tmax)
{
  jvec mcor=effective_mass(corr);
  
  mfit=constant_fit(mcor,tmin+1,tmax);
  jvec temp(corr.nel,corr.njack);
  for(int t=0;t<=TH;t++) temp[t]=sqrt(corr[t]/exp(-mfit*TH)/cosh(mfit*(TH-t))*mfit);
  ffit=constant_fit(temp,tmin,tmax);
  
  ofstream out(out_path);
  out<<"@type xydy"<<endl;
  //out<<"@line type none"<<endl;
  out<<mcor<<endl;
  out<<"&"<<endl;
  out<<"@type xy"<<endl;
  double av_mass=mfit.med();
  double er_mass=mfit.err();
  out<<tmin<<" "<<av_mass-er_mass<<endl;
  out<<tmax<<" "<<av_mass-er_mass<<endl;
  out<<tmax<<" "<<av_mass+er_mass<<endl;
  out<<tmin<<" "<<av_mass+er_mass<<endl;
}

int main()
{
  read_pars("input");
  read_input(data_list_file);
  init_latpars(100);
  
  int ncombo=nmass*nmass;
  bvec amD12(ncombo,nboot,njack);
  bvec afD12(ncombo,nboot,njack);
  
  //load all the corrs
  int ncombo_red=nmass*(nmass+1)/2;
  double *buf=new double[ncombo_red*T*(njack+1)];
  FILE *fin=open_file(combine("%s/%s",base_path,"P5P5").c_str(),"r");
  int stat=fread(buf,sizeof(double),ncombo_red*(njack+1)*T,fin);
  if(stat!=ncombo_red*(njack+1)*T)
    {
      cerr<<"Error loading data!"<<endl;
      exit(1);
    }
  
  int icombo=0;
  for(int ims=0;ims<nmass;ims++)
    for(int imc=ims;imc<nmass;imc++)
      {
	jvec corr(T,njack);
	corr.put(buf+icombo*T*(njack+1));
	
	int tmin=tmin_l;
	int tmax=tmax_l;
	if((ims>=1+nstrange)||(imc>=1+nstrange))
	  {
	    tmin=tmin_h;
	    tmax=tmax_h;
	  }
	
	jack tempm(njack),tempf(njack);
	fit(tempm,tempf,corr.simmetrized(1),combine("Ds_%d_%d",ims,imc).c_str(),tmin,tmax);

	boot_from_jack(amD12.data[icombo],tempm,ibeta);
	boot_from_jack(afD12.data[icombo],(mass[imc]+mass[ims])*tempf/(tempm*sinh(tempm)),ibeta);

	cout<<ims<<" "<<imc<<" "<<mass[ims]<<" "<<mass[imc]<<" "<<icombo<<" "<<amD12[icombo]<<" "<<tempf<<endl;
	
	icombo++;
      }
  
  boot mDs=interpolate_double(amD12)/lat[ibeta];
  boot fDs=interpolate_double(afD12)/lat[ibeta];
  
  mDs.write_to_binfile("results");
  fDs.append_to_binfile("results_fDs");

  cout<<mDs<<endl;
  cout<<fDs<<endl;

  return 0;
}
