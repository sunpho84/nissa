#include <include.h>
#include <sstream>
#include <iostream>

using namespace std;

int tmin,tmax;
int T=48,TH=24;
int njack;
int nmass=2;
int nslso,nslsi;

string print_error_band(jvec a,double *x)
{
  int nt=a.nel;
  ostringstream out;
  out<<"@type xy"<<endl;
  for(int t=0;t<nt;t++) out<<x[t]<<" "<<a[t].med()-a[t].err()<<endl;
  for(int t=nt-1;t>=0;t--) out<<x[t]<<" "<<a[t].med()+a[t].err()<<endl;
  out<<x[0]<<" "<<a[0].med()+a[0].err()<<endl;
  out<<"&\n@type xy"<<endl;
  for(int t=0;t<nt;t++) out<<x[t]<<" "<<a[t].med()<<endl;
  
  return out.str();
}

int iprop(int im3,int im2,int im1,int parity,int islso,int islsi)
{return 2*(islso+nslso*(islsi+nslsi*(parity+2*(im1+nmass*(im2+nmass*im3)))));}

jvec load_naz(const char *path,int ip)
{
  jvec temp(T,njack);
  temp.load_naz(path,ip);
  
  return temp;
}

jvec load_improved(const char *path,int r3,int r2,int r1,int islso,int islsi,int parity=0)
{
  return (
	  +load_naz(path,iprop( r3, r2, r1,parity,islso,islsi))
	  +load_naz(path,iprop(!r3,!r2,!r1,parity,islso,islsi))
	  
	  -
	  
	  (load_naz(path,iprop( r3, r2, r1,!parity,islso,islsi))+
	   load_naz(path,iprop(!r3,!r2,!r1,!parity,islso,islsi)))
	  .simmetric()
	  ).first_half()/4;
}

void read_input()
{
  FILE *fin=open_file("input","r");
  read_formatted_from_file_expecting((char*)&njack,fin,"%d","njack");
  read_formatted_from_file_expecting((char*)&nslso,fin,"%d","nslso");
  read_formatted_from_file_expecting((char*)&nslsi,fin,"%d","nslsi");
  read_formatted_from_file_expecting((char*)&tmin,fin,"%d","tmin");
  read_formatted_from_file_expecting((char*)&tmax,fin,"%d","tmax");
  fclose(fin);
}

int main()
{
  read_input();
  
  int r3=0;
  for(int r2=0;r2<2;r2++)
    for(int r1=0;r1<2;r1++)
      for(int islsi=0;islsi<nslsi;islsi++)
	for(int islso=0;islso<nslso;islso++)
	  {
	    cout<<"============="<<endl;
	    cout<<"Combo: "<<r3<<r2<<r1<<" "<<islso<<islsi<<endl;
	    cout<<"============="<<endl;
	    jvec a=load_improved("oC5C5o-sss_conf.1.dat",r3,r2,r1,islso,islsi);
	    jvec b=load_improved("oC5C5o-dss_conf.1.dat",r3,r2,r1,islso,islsi);
	    jvec c=load_improved("oC5C5o-sds_conf.1.dat",r3,r2,r1,islso,islsi);
	    jvec d=load_improved("oC5C5o-ssd_conf.1.dat",r3,r2,r1,islso,islsi);
	    
	    a.print_to_file(combine("a%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    b.print_to_file(combine("b%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    c.print_to_file(combine("c%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    d.print_to_file(combine("d%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    
	    jvec e=aperiodic_effective_mass(a);
	    e.print_to_file(combine("e%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    cout<<"MASS: "<<constant_fit(e,tmin,tmax)<<endl;
	    
	    (b/a).print_to_file(combine("r%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    (c/a).print_to_file(combine("s%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    (d/a).print_to_file(combine("t%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());

	    jvec u=(d+c-b)/a;
	    jvec v=(c-b)/a;

	    u.print_to_file(combine("u%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    v.print_to_file(combine("v%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    
	    jack SLOPE(njack),C(njack);
	    cout<<tmin<<" "<<tmax<<endl;
	    linear_fit(SLOPE,C,u,tmin,tmax);
	    
	    ofstream out_ratio(combine("ud_pres_%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    int nt=101;
	    int Dt=tmax-tmin;
	    double dt=Dt*1.1/(nt-1);
	    double t[nt];
	    jvec ratio_fitted(nt,njack);
	    for(int it=0;it<nt;it++)
	      {
		t[it]=tmin-0.05*Dt+dt*it;
		for(int ijack=0;ijack<=njack;ijack++)
		  ratio_fitted[it].data[ijack]=C[ijack]+SLOPE[ijack]*t[it];
	      }
	    out_ratio<<print_error_band(ratio_fitted,t);
	    out_ratio<<"@type xydy"<<endl
		     <<u<<"&"<<endl;

	    jvec ud(TH-1,njack);
	    for(int t=0;t<TH-1;t++) ud.data[t]=u[t+1]-u[t];
	    ud.print_to_file(combine("ud%d%d%d_%d%d",r3,r2,r1,islso,islsi).c_str());
	    cout<<"SLOPE: "<<SLOPE<<endl;
	    cout<<"SLOPE: "<<constant_fit(ud,tmin,tmax)<<endl;
	    linear_fit(SLOPE,C,d/a,tmin,tmax);
	    cout<<"SLOPE2: "<<SLOPE<<endl;
	    //cout<<"INTER: "<<C<<endl;
	  }

  
  
  return 0;
}
