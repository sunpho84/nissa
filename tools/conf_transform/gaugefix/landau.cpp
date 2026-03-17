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
  read_str_double("Alpha",&pars.alphaExp);
  read_str_int("NMaxIters",&pars.nmaxIterations);
  read_str_int("UseFftAcc",&pars.useFftAcc);
  read_str_int("UseAdaptativeSearch",&pars.useAdaptativeSearch);
  read_str_int("UseGeneralizedCG",&pars.useGeneralizedCg);
  int preliminaryRandomGaugeTransformation;
  read_str_int("PreliminaryRandomGaugeTransformation",&preliminaryRandomGaugeTransformation);
  int seed{};
  if(preliminaryRandomGaugeTransformation)
    read_str_int("Seed",&seed);

  int computeProfile;
  read_str_int("ComputeProfile",&computeProfile);
  char profilePath[1024];
  if(computeProfile)
    read_str_str("ProfilePath",profilePath,1024);
  
  close_input();
  
  //set pars
  field_rng_stream.init(seed);
  
  ///////////////////////////////////////////
  
  LxField<quad_su3> conf("Conf",WITH_HALO);
  LxField<quad_su3> fixedConf("FixedConf",WITH_HALO);
  
  read_ildg_gauge_conf(conf,inPath);
  if(preliminaryRandomGaugeTransformation)
    perform_random_gauge_transform(conf,conf);
  
  Landau_or_Coulomb_gauge_fix(fixedConf,pars,conf);
  
  write_ildg_gauge_conf(outPath,fixedConf);
  
  MASTER_PRINTF("plaq before: %16.16lg\n",global_plaquette_lx_conf(conf));
  MASTER_PRINTF("plaq after: %16.16lg\n",global_plaquette_lx_conf(fixedConf));
  
  if(computeProfile)
    {
      FILE* profileFile=open_file(profilePath,"w");
      
      for(int mu=0;mu<NDIM;mu++)
	{
	  LxField<su3> path("path",WITH_HALO);
	  PAR(0,
	      locVol,
	      CAPTURE(TO_WRITE(path)),
		  iVol,
	      {
		su3_put_to_id(path[iVol]);
	      });
	  
	  for(int x=0;x<=glbSize[mu]/2;x++)
	    {
	      LxField<double> t("t");
	      PAR(0,
		  locVol,
		  CAPTURE(TO_READ(path),
			  TO_WRITE(t)),
		  iVol,
		  {
		    t[iVol]=su3_real_trace(path[iVol]);
		  });
	      
	      double p;
	      t.reduce(p);
	      p/=glbVol;
	      
	      master_fprintf(profileFile,"%d %d %.16lg\n",mu,x,p);
	      
	      path.updateHalo();
	      LxField<su3> tmp("tmp");
	      PAR(0,
		  locVol,
		  CAPTURE(TO_READ(path),
			  TO_READ(fixedConf),
			  TO_WRITE(tmp),
			  mu),
		  iVol,
		  {
		    unsafe_su3_prod_su3(tmp[iVol],fixedConf[iVol][mu],path[loclxNeighup[iVol][mu]]);
		  });
	      
	      path=tmp;
	    }
	  
	  master_fprintf(profileFile,"\n");
	}
    }
}

int main(int narg,char** arg)
{
  initNissa(narg,arg);
  inMain(narg,arg);
  closeNissa();
  
  return 0;
}
