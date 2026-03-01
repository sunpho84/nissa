#include "nissa.hpp"

using namespace nissa;

void inMain(int narg,char **arg)
{
  if(narg<2)
    CRASH("Use: %s input_file",arg[0]);
  
  open_input(arg[1]);
  
  //Init the MPI grid
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  initGrid(T,L);
  
  char inPath[1024];
  read_str_str("InGaugePath",inPath,1024);
  
  char outPath[1024];
  read_str_str("OutGaugePath",outPath,1024);
  
  LC_gauge_fixing_pars_t pars;
  pars.gauge=LC_gauge_fixing_pars_t::LANDAU;
  read_str_double("Precision",&pars.targetPrecision);
  read_str_int("UseFftAcc",&pars.useFftAcc);
  read_str_int("UseAdaptativeSearch",&pars.useAdaptativeSearch);
  read_str_int("UseGeneralizedCG",&pars.useGeneralizedCg);
  
  close_input();
  
  //set pars
  field_rng_stream.init(3472291050);
  
  ///////////////////////////////////////////
  
  LxField<quad_su3> conf("Conf",WITH_HALO);
  LxField<quad_su3> fixedConf("FixedConf",WITH_HALO);
  
  read_ildg_gauge_conf(conf,inPath);
  
  Landau_or_Coulomb_gauge_fix(fixedConf,pars,conf);
  
  write_ildg_gauge_conf(outPath,fixedConf);
  
  MASTER_PRINTF("plaq before: %16.16lg\n",global_plaquette_lx_conf(conf));
  MASTER_PRINTF("plaq after: %16.16lg\n",global_plaquette_lx_conf(fixedConf));
}

int main(int narg,char** arg)
{
  initNissa(narg,arg);
  inMain(narg,arg);
  closeNissa();
  
  return 0;
}
