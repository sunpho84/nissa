#include <include.h>
#include <iostream>

using namespace std;

int nel;
int njack;
double *X;

void fit(const char *out_path,const char *title,jvec Y)
{
  //perform the fit
  jack m,q;
  linear_fit(m,q,X,Y);
  //calculate the chi2
  double C2=0;
  for(int i=0;i<nel;i++)
    {
      jack C=Y.data[i]-m*X[i]-q;
      C2+=pow(C.med()/C.err(),2);
    }
  cout<<"M = ("<<m<<"), Q=("<<q<<")"<<endl;
  cout<<"Chi2 = "<<C2<<" / "<<nel-2<<" = "<<C2/(nel-2)<<endl;
  
  //setup the plot
  grace out(out_path);
  out.plot_size(800,600);
  out.plot_title(combine("Chiral extrapolation of %s",title).c_str());
  out.axis_label("(aM\\s\\xp\\N)\\S2",title);
  //plot the function with error
  out.set(1,"orange");
  jvec par(2,njack);par.data[0]=-m;par.data[1]=-q;
  out.polygon(lin_fun,0,X[nel-1]*1.1,100,par);
  //plot the original data with error  
  out.new_set();
  out.set(3,"none","square","black");
  out.set_line_size(2);
  out.print_graph(X,-Y);
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      cerr<<"Error, use: "<<arg[0]<<" input"<<endl;
      exit(1);
    }
  
  ifstream input(arg[1]);
  
  input>>nel>>njack;
  
  X=new double[nel];
  jvec dM2K_P5P5(nel,njack);
  jvec dM2K_A0P5(nel,njack);
  jvec dfK_P5P5(nel,njack);
  jvec dfK_A0P5(nel,njack);
  
  for(int i=0;i<nel;i++)
    {
      char str[1024];
      input>>str;
      jvec A(5,njack);
      A.load(str,0);
      X[i]=A[0].med();
      dM2K_P5P5[i]=A[1];
      dM2K_A0P5[i]=A[2];
      dfK_P5P5[i]=A[3];
      dfK_A0P5[i]=A[4];
      
      cout<<X[i]<<" "<<dM2K_P5P5[i]<<endl;
    }
  
  fit("M2K_chiral_extrap.xmg","a\\xD\\0M\\S2\\N\\sK\\N",dM2K_P5P5);
  fit("fK_chiral_extrap.xmg","\\xD\\0f\\sK\\N/f\\sK",dfK_P5P5);
  
  return 0;
}
