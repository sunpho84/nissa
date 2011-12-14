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
  
  cout<<"pusc"<<endl;
  
  double im[5];
  bvec t(5,nboot,njack);
  for(int iel=0;iel<5;iel++)
    {
      t[iel]=phih[iel];
      im[iel]=1/mass[iel];
    }
  boot m,q;
  linear_fit(m,q,t,im);
  cout<<"m: "<<m<<" q: "<<q<<endl;

  {
    bvec yg(100,nboot,njack);
    double xg[100];
    for(int i=0;i<100;i++)
      {
	xg[i]=1.3/100*i;
	yg[i]=m*xg[i]+q;
      }
    
    grace out("phi_vs_m.xmg");
    out.fout<<"@type xy"<<endl;
    out.polygon(xg,yg);
    out.fout<<"&"<<endl;
    out.ave_line(xg,yg);
    out.fout<<"@type xydy"<<endl;
    for(int im=0;im<nm;im++)
      {
	out.fout<<1/mass[im]<<" "<<phih[im]<<endl;
	if(im==4) out.fout<<"&"<<endl;
      }
  }
  
  return 0;
}
