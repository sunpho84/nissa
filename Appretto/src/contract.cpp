#pragma once

//Trace of the product of gamma1 * spinspin1^dag * gamma2 * spinspin2
void trace_g_sdag_g_s(complex c,dirac_matr &g1,spinspin &s1,dirac_matr &g2,spinspin &s2)
{
  spinspin t1,t2;
  
  spinspin_dirac_spinspindag_prod(t1,g1,s1);
  if(rank==0)
    {
      print_spinspin(t1);
      cout<<endl;
    }

  spinspin_dirac_spinspin_prod(t2,g2,s2);
  if(rank==0)
    {
      print_spinspin(t2);
      cout<<endl;
    }

  trace_prod_spinspins(c,t1,t2);
}

