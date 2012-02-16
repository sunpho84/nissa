#include "include.h"

int nboot=100;
int njack=16;
int nm=11;
double mass[11]={0.984873,1.15836,1.36255,1.6023,1.88462,2.2165,2.60711,3.06615,3.60653,4.24174,4.98902};
bvec Mh(nm,nboot,njack);

int main()
{
  FILE *an_input_file=open_file("analysis_pars","r");
  char chiral_data[1024];
  read_formatted_from_file_expecting(chiral_data,an_input_file,"%s","chiral_data");
  fclose(an_input_file);
  Mh.load(chiral_data,0);
  
  bvec par(3,nboot,njack);
  int last=0;
  while(!isnan(Mh[last].err()) && last<nm)
    last++;

  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      double y[3];
      double *x=mass+last-3;
      y[0]=Mh[last-3][iboot];
      y[1]=Mh[last-2][iboot];
      y[2]=Mh[last-1][iboot];
      parabolic_spline(par.data[0].data[iboot],par.data[1].data[iboot],par.data[2].data[iboot],x,y);
    }
  
  double MB_phys=5.279;
  boot a=par[0],b=par[1],c=par[2];
  boot d=b*b-4*a*(c-MB_phys);
  boot Mb_phys=(-b+sqrt(d))/(2*a);
  cout<<Mb_phys<<endl;
  
  grace out("M_vs_m.xmg");
  out.fout<<"@type xy"<<endl;
  int nx=100;
  double x[nx],x0=mass[last-3],x1=Mb_phys.med()+2*Mb_phys.err(),dx=(x1-x0)/(nx-1);
  bvec y(nx,nboot,njack);
  for(int i=0;i<nx;i++)
    {
      x[i]=x0+dx*i;
      y[i]=a*x[i]*x[i]+b*x[i]+c;
    }
  out.polygon(x,y);
  out.fout<<"@type xydy"<<endl;
  for(int im=0;im<nm;im++) out.fout<<mass[im]<<" "<<Mh[im]<<endl;

  out.fout<<"@type xydx"<<endl;
  out.fout<<Mb_phys.med()<<" "<<MB_phys<<" "<<Mb_phys.err()<<endl;
  
  out.fout<<"&\n@type xy"<<endl;
  out.ave_line(x,y);
  ofstream err_perc("err_perc.xmg");
  for(int im=0;im<nm;im++)
    err_perc<<mass[im]<<" "<<Mh[im].err()/Mh[im].med()*100<<endl;
  
  out.fout<<"&\n0 "<<MB_phys<<"\n 10 "<<MB_phys<<endl;
  
  cout<<"pusc"<<endl;
  
  return 0;
}
