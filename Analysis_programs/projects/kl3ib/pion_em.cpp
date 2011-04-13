#include "include.h"
#include "kl3_common.cpp"

int main()
{
  njack=10;
  read_input();
  
  //allocate vector where to load data
  jvec P5_V0_P5[nmoms][nmoms],P5_P5[nmoms];
  
  //load two and three point for Zv
  int r=0,im_spec=0,im_valence=0;
  for(int ik1=0;ik1<nmoms;ik1++)
    {
      P5_P5[ik1]=read_ch_thimpr_P5_P5(im_spec,im_valence,ik1,r);
      for(int ik2=0;ik2<nmoms;ik2++) P5_V0_P5[ik1][ik2]=read_deg_ch_thimpr_P5_V0_P5(im_spec,im_valence,ik1,ik2,r);
    }
  
  //calculate Zv
  jvec two_Zv=P5_P5[0][TH]/P5_V0_P5[0][0];
  
  //fit mass
  jack M[nmoms];
  for(int ik=0;ik<nmoms;ik++) M[ik]=constant_fit(effective_mass(P5_P5[ik].simmetrized(1)),12,23);
  
  for(int ik1=0;ik1<nmoms;ik1++)
    for(int ik2=0;ik2<nmoms;ik2++)
      if(theta[ik2]>=0)
        {
          //calculate the double ratio simmetrized and print
          jvec drs=(P5_V0_P5[ik1][ik2]/sqrt(P5_P5[ik1][TH]*P5_P5[ik2][TH])*two_Zv).simmetrized(1);
	  drs.print_to_file("EM_ff_mom_%d_%d.out",ik1,ik2);
	  
          //calculate the transferred impulse and the form factor
          jack Q2=calculate_Q2(M[ik1],theta[ik1],M[ik2],theta[ik2]);
          jack C=constant_fit(drs,10,14);
          
	  //print Q2 and form factor
	  cout<<Q2<<" "<<C<<endl;
	}

  //print results
  for(int ik1=0;ik1<nmoms;ik1++)
    {
      int ik2=iopp_th[ik1];
      P5_P5[ik1].print_to_file("P5_P5_mom%d.out",ik1);
      P5_V0_P5[ik1][ik2].print_to_file("P5_V0_P5_mom_%d_%d.out",ik1,ik2);
    }

  return 0;
}
