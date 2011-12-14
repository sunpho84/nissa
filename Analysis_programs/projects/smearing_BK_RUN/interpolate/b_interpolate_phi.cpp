#include "include.h"

int nboot=100;
int njack=16;
int nm=11;
double mass[11]={0.984873,1.15836,1.36255,1.6023,1.88462,2.2165,2.60711,3.06615,3.60653,4.24174,4.98902};
double inv_mass[11];
bvec phih(nm,nboot,njack);

void linear_fit(boot &m,boot &q,bvec corr,double *x)
{
  int nel=corr.nel;
  int nboot=corr.nboot;
  int njack=corr.njack;
  boot Y(nboot,njack),XY(nboot,njack),Y2(nboot,njack);
  double X=0,W=0,X2=0;

  Y=XY=0;
  for(int i=0;i<nel;i++)
    {
      double err=corr[i].err();
      double w=1/sqr(err);
      
      W+=w;
      
      X+=x[i]*w;
      X2+=x[i]*x[i]*w;
      
      Y+=corr[i]*w;
      Y2+=sqr(corr[i])*w;
      
      XY+=x[i]*corr[i]*w;
    }

  XY-=X*Y/W;
  Y2-=Y*Y/W;
  X2-=X*X/W;
  
  m=XY/X2;
  q=(Y-m*X)/W;
}


int main()
{
  FILE *an_input_file=open_file("analysis_pars","r");
  char chiral_data[1024];
  read_formatted_from_file_expecting(chiral_data,an_input_file,"%s","chiral_data");
  fclose(an_input_file);
  
  phih.load(chiral_data,0);
  
  for(int imass=0;imass<nm;imass++) inv_mass[imass]=1/mass[imass];
  
  bvec par(3,nboot,njack);
  for(int iboot=0;iboot<nboot+1;iboot++)
    {
      double y[3];
      double *x=inv_mass+8;
      y[0]=phih[8][iboot];
      y[1]=phih[9][iboot];
      y[2]=phih[10][iboot];
      parabolic_spline(par.data[0].data[iboot],par.data[1].data[iboot],par.data[2].data[iboot],x,y);
    }
  
  boot Mb_phys(nboot,njack);
  Mb_phys.fill_gauss(4.91,0.14,124234);
  cout<<Mb_phys<<endl;
  
  grace out("phi_vs_m.xmg");
  
  out.fout<<"@type xy"<<endl;
  int nx=100;
  double x[nx],x0=1/(Mb_phys.med()+2*Mb_phys.err()),x1=1/mass[8],dx=(x1-x0)/(nx-1);
  bvec y(nx,nboot,njack);
  for(int i=0;i<nx;i++)
    {
      x[i]=x0+dx*i;
      y[i]=par[0]*x[i]*x[i]+par[1]*x[i]+par[2];
    }
  out.polygon(x,y);
  
  out.fout<<"@type xydy"<<endl;
  for(int im=0;im<nm;im++) out.fout<<1/mass[im]<<" "<<phih[im]<<endl;
  
  out.fout<<"&\n@type xy"<<endl;
  out.ave_line(x,y);
  ofstream err_perc("err_perc.xmg");
  for(int im=0;im<nm;im++)
    err_perc<<mass[im]<<" "<<phih[im].err()/phih[im].med()*100<<endl;
  
  cout<<"pusc"<<endl;
  
  return 0;
}
