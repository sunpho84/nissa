#ifndef _APE_HPP
#define _APE_HPP

#include <sstream>

#include <base/old_field.hpp>
#include <geometry/geometry_lx.hpp>
#include <new_types/su3_op.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  //parameters to ape smear
  struct ape_pars_t
  {
    int nlevels;
    
    double alpha;
    
    int def_nlevels() const
    {
      return 20;
    }
    
    double def_alpha() const
    {
      return 0.5;
    }
    
    int master_fprintf(FILE *fout,
		       const bool full) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      os<<"Ape\n";
      if(full or is_nonstandard())
	{
	  if(full or alpha!=def_alpha()) os<<" Alpha\t\t=\t"<<alpha<<"\n";
	  if(full or nlevels!=def_nlevels()) os<<" NLevels\t=\t"<<nlevels<<"\n";
	}
      
      return os.str();
    }
    
    int is_nonstandard() const
    {
      return
	nlevels!=def_nlevels() or
	alpha!=def_alpha();
    }
    
    ape_pars_t() :
      nlevels(def_nlevels()),
      alpha(def_alpha())
    {
    }
  };
  
  inline void ape_smear_conf(LxField<quad_su3>& smearConf,
			     LxField<quad_su3> conf,
			     const double& alpha,
			     const int& nSteps,
			     const which_dir_t& dirs=all_dirs,
			     const int& minStapleDir=0)
  {
    char listed_dirs[21]="";
    for(int mu=0;mu<NDIM;mu++)
      if(dirs[mu])
	{
	  char temp[3];
	  snprintf(temp,3,"%d ",mu);
	  strncat(listed_dirs,temp,20);
	}
    verbosity_lv1_master_printf("APE { %s} smearing with alpha=%g, %d iterations\n",listed_dirs,alpha,nSteps);
    
    for(int istep=1;istep<=nSteps;istep++)
      {
	verbosity_lv3_master_printf("APE spatial smearing with alpha=%g iteration %d of %d\n",alpha,istep,nSteps);
	
	conf.updateEdges();
	
	PAR(0,locVol,
	    CAPTURE(alpha,dirs,minStapleDir,
		    TO_READ(conf),
		    TO_WRITE(smearConf)),
	    ivol,
	    {
	      for(int mu=0;mu<NDIM;mu++)
		if(dirs[mu])
		  {
		    //calculate staples
		    su3 stap,temp1,temp2;
		    su3_put_to_zero(stap);
		    
		    for(int nu=minStapleDir;nu<NDIM;nu++)     //  E---F---C
		      if(nu!=mu)                              //  |   |   | mu
			{                                     //  D---A---B
			  const int A=ivol;                   //   nu
			  const int B=loclxNeighup[A][nu];
			  const int F=loclxNeighup[A][mu];
			  unsafe_su3_prod_su3(temp1,conf[A][nu],conf[B][mu]);
			  unsafe_su3_prod_su3_dag(temp2,temp1,conf[F][nu]);
			  su3_summ(stap,stap,temp2);
			  
			  const int D=loclxNeighdw[A][nu];
			  const int E=loclxNeighup[D][mu];
			  unsafe_su3_dag_prod_su3(temp1,conf[D][nu],conf[D][mu]);
			  unsafe_su3_prod_su3(temp2,temp1,conf[E][nu]);
			  su3_summ(stap,stap,temp2);
			}
		    
		    //create new link to be reunitarized
		    su3 prop_link;
		    for(int icol1=0;icol1<NCOL;icol1++)
		      for(int icol2=0;icol2<NCOL;icol2++)
			for(int ri=0;ri<2;ri++)
			  //prop_link[icol1][icol2][ri]=(1-alpha)*temp_conf[ivol][mu][icol1][icol2][ri]+alpha/6*stap[icol1][icol2][ri];
			  prop_link[icol1][icol2][ri]=conf[ivol][mu][icol1][icol2][ri]+alpha*stap[icol1][icol2][ri];
		    
		    su3_unitarize_maximal_trace_projecting(prop_link);
		    
		    su3_copy(smearConf[ivol][mu],prop_link);
		  }
	    });
	
	if(istep!=nSteps)
	  conf=smearConf;
	else
	  smearConf.invalidateHalo();
      }
  }
  
  inline void ape_single_dir_smear_conf(LxField<quad_su3>& smear_conf,
					const LxField<quad_su3>& origi_conf,
					const double& alpha,
					const int& nstep,
					const int& mu,
					const int& min_staple_dir=0)
  {
    ape_smear_conf(smear_conf,origi_conf,alpha,nstep,only_dir[mu],min_staple_dir);
  }
  
  inline void ape_perp_dir_smear_conf(LxField<quad_su3>& smear_conf,
				      const LxField<quad_su3>& origi_conf,
				      const double& alpha,
				      const int& nstep,
				      const int& mu,
				      const int& min_staple_dir=0)
  {
    ape_smear_conf(smear_conf,origi_conf,alpha,nstep,all_other_dirs[mu]);
  }
  
  inline void ape_temporal_smear_conf(LxField<quad_su3>& smear_conf,
				      const LxField<quad_su3>& origi_conf,
				      const double& alpha,
				      const int& nstep)
  {
    ape_single_dir_smear_conf(smear_conf,origi_conf,alpha,nstep,0);
  }
  
  inline void ape_spatial_smear_conf(LxField<quad_su3>& smear_conf,
				     const LxField<quad_su3>& origi_conf,
				     const double& alpha,
				     const int& nstep)
  {
    ape_perp_dir_smear_conf(smear_conf,origi_conf,alpha,nstep,0,1);
  }
}

#endif
