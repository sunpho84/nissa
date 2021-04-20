#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "random/randomGenerate.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/linalgs.hpp"
#include "measures/fermions/tm_corr_op.hpp"
#include "new_types/su3.hpp"
#include "nucleon.hpp"
#include "operations/gauge_fixing.hpp"
#include "operations/smearing/APE.hpp"
#include "operations/smearing/gaussian.hpp"
#include "routines/mpi_routines.hpp"
#include "stag.hpp"

namespace nissa
{
  using namespace stag;
  
  void measure_nucleon_corr(eo_ptr<quad_su3> conf,theory_pars_t theory_pars,nucleon_corr_meas_pars_t meas_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    
    /// Number of flavors
    const int nflavs=theory_pars.nflavs();
    
    /// Source
    spincolor** source=nissa_malloc("source",NDIRAC*NCOL,spincolor*);
    for(int idc=0;idc<NDIRAC*NCOL;idc++)
      source[idc]=nissa_malloc("source[]",locVolWithBord.nastyConvert(),spincolor);
    
    /// Smearing conf
    quad_su3* smearingConf=nissa_malloc("smearingConf",locVolWithBord.nastyConvert(),quad_su3);
    paste_eo_parts_into_lx_vector(smearingConf,conf);
    ape_smear_conf(smearingConf,smearingConf,meas_pars.apeSmeAlpha,meas_pars.apeSmeNSteps,all_other_dirs[0],1);
    
    /// Propagators
    spincolor*** prop;
    prop=nissa_malloc("prop",nflavs,spincolor**);
    for(int iflav=0;iflav<nflavs;iflav++)
      {
	prop[iflav]=nissa_malloc("prop[iflav]",NDIRAC*NCOL,spincolor*);
	for(int idc=0;idc<NDIRAC*NCOL;idc++)
	  prop[iflav][idc]=nissa_malloc("prop[iflav][idc]",locVolWithBord.nastyConvert(),spincolor);
      }
    
    /// Operations for the propagator
    tm_corr_op tmCorrOp(conf,meas_pars.residue,theory_pars);
    
    /// Correlation function
    complex* corr=nissa_malloc("corr",glbTimeSize.nastyConvert()*nflavs*nflavs,complex);
    
    for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
      {
	vector_reset(corr);
	
	for(int ihit=0;ihit<meas_pars.nhits;ihit++)
	  {
	    /// Position of the source
	    GlbCoords glbSourceCoords;
	    generate_random_coord(glbSourceCoords);
	    
	    /// Which rank hosts the source
	    Rank whichRank;
	    
	    /// Local site
	    LocLxSite locSourcePos;
	    
	    get_loclx_and_rank_of_coord(locSourcePos,whichRank,glbSourceCoords);
	    
	    for(int idirac=0;idirac<NDIRAC;idirac++)
	      for(int icol=0;icol<NCOL;icol++)
		{
		  spincolor* &s=source[icol+NCOL*idirac];
		  vector_reset(s);
		  
		  if(rank==whichRank)
		    s[locSourcePos.nastyConvert()][idirac][icol][RE]=1;
		  set_borders_invalid(s);
		  
		  gaussian_smearing(s,s,smearingConf,meas_pars.gaussSmeKappa,meas_pars.gaussSmeNSteps);
		}
	    
	    for(int iflav=0;iflav<nflavs;iflav++)
	      for(int idc=0;idc<NDIRAC*NCOL;idc++)
		{
		  spincolor* p=prop[iflav][idc];
		  tmCorrOp.inv(p,source[idc],iflav);
		  
		  gaussian_smearing(p,p,smearingConf,meas_pars.gaussSmeKappa,meas_pars.gaussSmeNSteps);
		}
	    
	    for(int ilikeFlav=0;ilikeFlav<nflavs;ilikeFlav++)
	      for(int idislikeFlav=0;idislikeFlav<nflavs;idislikeFlav++)
		{
		  master_printf("Computing %d %d\n",ilikeFlav,idislikeFlav);
		  complex tempCorr[glbTimeSize.nastyConvert()];
		  tm_corr_op::compute_nucleon_2pts_contr(tempCorr,
							 prop[ilikeFlav],
							 prop[idislikeFlav],
							 glbSourceCoords(tDir),-1);
		  
		  master_printf("Summing %d %d\n",ilikeFlav,idislikeFlav);
		  FOR_ALL_GLB_TIMES(t)
		    complex_summassign(corr[t.nastyConvert()+glbTimeSize.nastyConvert()*(ilikeFlav+nflavs*idislikeFlav)],tempCorr[t.nastyConvert()]);
		}
	  }
	
	for(int ilikeFlav=0;ilikeFlav<nflavs;ilikeFlav++)
	  for(int idislikeFlav=0;idislikeFlav<nflavs;idislikeFlav++)
	    {
	      master_printf("Printing %d %d\n",ilikeFlav,idislikeFlav);
	      master_fprintf(file," # conf %d ; like1 = %d ; dislike = %d ; like2 = %d\n",
			     iconf,ilikeFlav,idislikeFlav,ilikeFlav);
	      
	      FOR_ALL_GLB_TIMES(t)
		{
		  complex c;
		  complex_prod_double(c,corr[t.nastyConvert()+glbTimeSize.nastyConvert()*(ilikeFlav+nflavs*idislikeFlav)],1.0/meas_pars.nhits);
		  master_fprintf(file,"%d %+.16lg %+.16lg\n",t,c[RE],c[IM]);
		}
	    }
      }
    
    for(int idc=0;idc<NDIRAC*NCOL;idc++)
      nissa_free(source[idc]);
    nissa_free(source);
    
    for(int iflav=0;iflav<nflavs;iflav++)
      {
	for(int idc=0;idc<NDIRAC*NCOL;idc++)
	  nissa_free(prop[iflav][idc]);
	nissa_free(prop[iflav]);
      }
    nissa_free(prop);
    
    nissa_free(corr);
    
    nissa_free(smearingConf);
    
    close_file(file);
  }
  
  //nucleon correlators
  std::string nucleon_corr_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasNucleonCorrs\n";
    
    if(is_nonstandard() or full) os<<base_fermionic_meas_t::get_str(full);
    if(gaussSmeKappa!=def_gaussSmeKappa() or full) os<<"GaussSmeKappa\t=\t"<<gaussSmeKappa<<"\n";
    if(gaussSmeNSteps!=def_gaussSmeNSteps() or full) os<<"GaussSmeNSteps\t=\t"<<gaussSmeNSteps<<"\n";
    if(apeSmeAlpha!=def_apeSmeAlpha() or full) os<<"ApeSmeAlpha\t=\t"<<apeSmeAlpha<<"\n";
    if(apeSmeNSteps!=def_apeSmeNSteps() or full) os<<"ApeSmeNSteps\t=\t"<<apeSmeNSteps<<"\n";
    
    return os.str();
  }
}
