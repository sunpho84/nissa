#include <include.h>
#include <iostream>

#include "prace_common.cpp"

using namespace std;

int main()
{
  read_data_list();
  
  jvec E(nmoms,njack);
  jvec d[nmoms];
  
  for(int ith2=0;ith2<nmoms;ith2++)
    {
      int im2=0;
      int ru=0;
      
      jvec c(T,njack);
      c=(read_ch_thimpr_P5_P5(ith2,im2,ru,im_spec)+read_ch_thimpr_P5_P5(ith2,im2,!ru,im_spec))/2;
      d[ith2]=effective_mass(c.simmetrized(1));
      E.data[ith2]=constant_fit(d[ith2],9,19);
      
      //output plot
      grace out("/tmp/out_%d",ith2);
      out.plot_size(800,600);
      //function with error
      out.set(1,"orange");
      jvec par(1,njack);par.data[0]=E[ith2];
      out.polygon(const_fun,9,19,100,par);
      //plot the original data with error  
      out.new_set();
      out.set(3,"none","square","black");
      out.set_line_size(2);
      out.print_graph(d[ith2]);
    }
  
  return 0;
}
