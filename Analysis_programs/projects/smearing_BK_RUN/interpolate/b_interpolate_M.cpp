#include "include.h"
#include "../../nf2/common_pars.cpp"

bvec Mh(nref_hmass,nboot,njack);

int main()
{
  double mb_phys_med=4.98;
  double mb_phys_err=0.2;
  boot mb_phys(nboot,njack);
  mb_phys.fill_gauss(mb_phys_med,mb_phys_err,2432112);

  init_latpars();
  
  FILE *an_input_file=open_file("analysis_pars","r");
  char chiral_data[1024];
  read_formatted_from_file_expecting(chiral_data,an_input_file,"%s","chiral_data");
  fclose(an_input_file);
  Mh.load(chiral_data,0);
  
  double MB_phys=5.279;
  
  boot yout=interpolate_single(ref_hmass,Mh,mb_phys,"M_vs_m.xmg");
  
  cout<<"C_interpolated: "<<interpolate_single(ref_hmass,Mh,mc_phys,"M_vs_m.xmg")<<endl;
  cout<<"B_interpolated: "<<yout<<endl;
  
  return 0;
}
