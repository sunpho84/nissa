#pragma once

//load the gauge configuration, compute plaquette and verify it to be
//correct within expected precision
int test_plaquette_computation()
{
  const double exp_plaq=0.6350055722288719;
  const double tolerance=1.e-14;
  
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //read conf, compute plaquette, print it
  read_ildg_gauge_conf(conf,"../data/L4T8conf");
  
  double plaq=global_plaquette_lx_conf(conf);
  master_printf("Loaded plaquette: %16.16g, expected: %16.16g\n",plaq,exp_plaq);
  
  //check precision
  double rel_diff=(plaq-exp_plaq)/exp_plaq;
  int test_passed=fabs(rel_diff)<=tolerance;  
  master_printf("Relative difference: %16.16g\n",rel_diff);
  master_printf("Tolerance: %16.16g\n",tolerance);
  
  nissa_free(conf);
  
  return test_passed;
}
