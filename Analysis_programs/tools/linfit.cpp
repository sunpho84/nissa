#include <include.h>
#include <iostream>

using namespace std;

int main(int narg,char **arg)
{
  if(narg<2)
    {
      cerr<<"Error, use: "<<arg[0]<<" input"<<endl;
      exit(1);
    }
  
  ifstream input(arg[1]);
  int nel;
  int njack;
  
  input>>nel>>njack;
  
  double X[nel];
  jvec Y(nel,njack);
  
  for(int i=0;i<nel;i++)
    {
      char str[1024];
      input>>str;
      jvec A(2,njack);
      A.load(str,0);
      X[i]=A[0].med();
      Y[i]=A[1];
      
      cout<<X[i]<<" "<<Y[i]<<endl;
    }
  
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
  grace outh("/tmp/res");
  outh.plot_size(800,600);
  outh.plot_title("Linear fit");
  //plot the function with error
  outh.set(1,"orange");
  jvec par(2,njack);par.data[0]=m;par.data[1]=q;  
  outh.polygon(lin_fun,X[0],X[nel-1],100,par);
  //plot the original data with error  
  outh.new_set();
  outh.set(3,"none","square","black");
  outh.set_line_size(2);
  outh.print_graph(X,Y);
  
  return 0;
}
