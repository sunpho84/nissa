#include <include.h>
#include <iostream>

using namespace std;

int tmin,tmax;
int T=48,TH=24;
int njack;
int nmass=2;
int nsl;

int iprop(int im3,int im2,int im1,int parity,int islsink,int islsrc)
{return 2*(islsrc+nsl*(islsink+nsl*(parity+2*(im1+nmass*(im2+nmass*im3)))));}

jvec load_naz(const char *path,int ip)
{
  jvec temp(T,njack);
  temp.load_naz(path,ip);
  
  return temp;
}

jvec load_improved(const char *path,int r3,int r2,int r1,int islsink,int islsrc,int parity=0)
{
  return (
	  +load_naz(path,iprop( r3, r2, r1,parity,islsrc,islsink))
	  +load_naz(path,iprop(!r3,!r2,!r1,parity,islsrc,islsink))
	  
	  -
	  
	  (load_naz(path,iprop( r3, r2, r1,!parity,islsrc,islsink))+
	  load_naz(path,iprop(!r3,!r2,!r1,!parity,islsrc,islsink))).simmetric()
	  ).first_half()/4;
}

void read_input()
{
  FILE *fin=open_file("input","r");
  read_formatted_from_file_expecting((char*)&njack,fin,"%d","njack");
  read_formatted_from_file_expecting((char*)&nsl,fin,"%d","nsl");
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
      for(int islsrc=0;islsrc<nsl;islsrc++)
	for(int islsink=0;islsink<nsl;islsink++)
	  {
	    cout<<"============="<<endl;
	    cout<<"Combo: "<<r3<<r2<<r1<<" "<<islsrc<<islsink<<endl;
	    cout<<"============="<<endl;
	    jvec a=load_improved("oC5C5o-sss_conf.1.dat",r3,r2,r1,islsrc,islsink);
	    jvec b=load_improved("oC5C5o-dss_conf.1.dat",r3,r2,r1,islsrc,islsink);
	    jvec c=load_improved("oC5C5o-sds_conf.1.dat",r3,r2,r1,islsrc,islsink);
	    jvec d=load_improved("oC5C5o-ssd_conf.1.dat",r3,r2,r1,islsrc,islsink);
	    
	    a.print_to_file(combine("a%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    b.print_to_file(combine("b%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    c.print_to_file(combine("c%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    d.print_to_file(combine("d%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    
	    jvec e=aperiodic_effective_mass(a);
	    e.print_to_file(combine("e%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    cout<<"MASS: "<<constant_fit(e,tmin,tmax)<<endl;
	    
	    (b/a).print_to_file(combine("r%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    (c/a).print_to_file(combine("s%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    (d/a).print_to_file(combine("t%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());

	    jvec u=(d+c-b)/a;
	    jvec v=(c-b)/a;

	    u.print_to_file(combine("u%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    v.print_to_file(combine("v%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    
	    jack SLOPE(njack),C(njack);
	    cout<<tmin<<" "<<tmax<<endl;
	    linear_fit(SLOPE,C,u,tmin,tmax);
	    jvec ud(TH-1,njack);
	    for(int t=0;t<TH-1;t++) ud.data[t]=u[t+1]-u[t];
	    ud.print_to_file(combine("ud%d%d%d_%d%d",r3,r2,r1,islsrc,islsink).c_str());
	    cout<<"SLOPE: "<<SLOPE<<endl;
	    cout<<"SLOPE: "<<constant_fit(ud,tmin,tmax)<<endl;
	    linear_fit(SLOPE,C,d/a,tmin,tmax);
	    cout<<"SLOPE2: "<<SLOPE<<endl;
	    //cout<<"INTER: "<<C<<endl;
	  }

  
  
  return 0;
}
