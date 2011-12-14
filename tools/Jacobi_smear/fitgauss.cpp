#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <cstdlib>

#include <fstream>
#include <iostream>

using namespace std;

int main(int narg,char **arg)
{
  int i;
  double x,y;
  ifstream file;
  TApplication myapp("App",NULL,NULL);
  TCanvas tela;
  TGraphErrors graf;
  TF1 fun("fun","[0]*exp(-x*x/(2*[1]*[1]))");

  if(narg<2)
    {
      cerr<<"usa: "<<arg[0]<<" file"<<endl;
      exit(0);
    }

  file.open(arg[1]);

  if(file.good()==0)
    {
      cout<<"File: "<<arg[1]<<" inesistente!"<<endl;
      exit(0);
    }

  i=0;
  double dum,N=1;
  while(file>>x>>dum>>y)
    {
      if(i==0) N=y;
      graf.SetPoint(i,x,y);
      graf.SetPointError(i,0,0.0001);
      i++;
    }
  
  fun.SetParameter(0,N);
  fun.SetParameter(1,1);

  fun.SetParName(0,"N");
  fun.SetParName(1,"S");

  graf.Draw("AP");

  graf.Fit("fun");

  //graf.GetXaxis()->SetLimits(0,x);

  tela.Modified();
  tela.Update();

  myapp.Run(true);

  return 0;
}
