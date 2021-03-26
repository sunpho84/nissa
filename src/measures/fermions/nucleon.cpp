#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/linalgs.hpp"
#include "measures/fermions/tm_corr_op.hpp"
#include "new_types/su3.hpp"
#include "nucleon.hpp"
#include "operations/gauge_fixing.hpp"
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
    spincolor* source=nissa_malloc("source",locVol+bord_vol,spincolor);
    
    /// Propagators
    spincolor*** prop;
    prop=nissa_malloc("prop",nflavs,spincolor**);
    for(int iflav=0;iflav<nflavs;iflav++)
      {
	prop[iflav]=nissa_malloc("prop[iflav]",NDIRAC*NCOL,spincolor*);
	for(int idc=0;idc<NDIRAC*NCOL;idc++)
	  prop[iflav][idc]=nissa_malloc("prop[iflav][idc]",locVol+bord_vol,spincolor);
      }
    
    /// Operations for the propagator
    tm_corr_op tmCorrOp(conf,meas_pars.residue,theory_pars);
    
    /// Correlation function
    complex* corr=nissa_malloc("corr",glbSize[0]*nflavs*nflavs,complex);
    
    for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
      {
	vector_reset(corr);
	
	for(int ihit=0;ihit<meas_pars.nhits;ihit++)
	  {
	    /// Position of the source
	    coords glbSourceCoords;
	    generate_random_coord(glbSourceCoords);
	    
	    for(int iflav=0;iflav<nflavs;iflav++)
	      for(int idirac=0;idirac<NDIRAC;idirac++)
		for(int icol=0;icol<NCOL;icol++)
		  {
		    vector_reset(source);
		    
		    /// Which rank hosts the source
		    int whichRank;
		    
		    /// Local site
		    int locSourcePos;
		    
		    get_loclx_and_rank_of_coord(&locSourcePos,&whichRank,glbSourceCoords);
		    
		    if(rank==whichRank)
		      source[locSourcePos][idirac][icol][RE]=1;
		    
		    tmCorrOp.inv(prop[iflav][icol+NCOL*idirac],source,iflav);
		  }
	    
	    for(int ilikeFlav=0;ilikeFlav<nflavs;ilikeFlav++)
	      for(int idislikeFlav=0;idislikeFlav<nflavs;idislikeFlav++)
		{
		  master_printf("Computing %d %d\n",ilikeFlav,idislikeFlav);
		  complex tempCorr[glbSize[0]];
		  tm_corr_op::compute_nucleon_2pts_contr(tempCorr,
							 prop[ilikeFlav],
							 prop[idislikeFlav],
							 glbSourceCoords[0],-1);
		  
		  master_printf("Summing %d %d\n",ilikeFlav,idislikeFlav);
		  for(int t=0;t<glbSize[0];t++)
		    complex_summassign(corr[t+glbSize[0]*(ilikeFlav+nflavs*idislikeFlav)],tempCorr[t]);
		}
	  }
	
	for(int ilikeFlav=0;ilikeFlav<nflavs;ilikeFlav++)
	  for(int idislikeFlav=0;idislikeFlav<nflavs;idislikeFlav++)
	    {
	      master_printf("Printing %d %d\n",ilikeFlav,idislikeFlav);
	      master_fprintf(file," # conf %d ; like1 = %d ; dislike = %d ; like2 = %d\n",
			     iconf,ilikeFlav,idislikeFlav,ilikeFlav);
	      
	      for(int t=0;t<glbSize[0];t++)
		{
		  complex c;
		  complex_prod_double(c,corr[t+glbSize[0]*(ilikeFlav+nflavs*idislikeFlav)],1.0/meas_pars.nhits);
		  master_fprintf(file,"%d %+.16lg %+.16lg\n",t,c[RE],c[IM]);
		}
	    }
      }
    
    nissa_free(source);
    
    for(int iflav=0;iflav<nflavs;iflav++)
      {
	for(int idc=0;idc<NDIRAC*NCOL;idc++)
	  nissa_free(prop[iflav][idc]);
	nissa_free(prop[iflav]);
      }
    nissa_free(prop);
    
    nissa_free(corr);
    
    close_file(file);
  }
  
  //nucleon correlators
  std::string nucleon_corr_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasNucleonCorrs\n";
    if(is_nonstandard() or full) os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
