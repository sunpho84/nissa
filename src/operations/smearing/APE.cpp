#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/field.hpp>
#include <base/vectors.hpp>
#include <communicate/borders.hpp>
#include <communicate/edges.hpp>
#include <geometry/geometry_lx.hpp>
#include <linalgs/linalgs.hpp>
#include <new_types/su3_op.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  //perform ape smearing
  //be sure not to have border condition added
  void ape_smear_conf(LxField<quad_su3>& smear_conf,
		      LxField<quad_su3> conf,
		      const double& alpha,
		      const int& nSteps,
		      const which_dir_t& dirs,
		      const int& min_staple_dir)
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
	
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  for(int mu=0;mu<NDIM;mu++)
	    if(dirs[mu])
	      {
		//calculate staples
		su3 stap,temp1,temp2;
		su3_put_to_zero(stap);
		
		for(int nu=min_staple_dir;nu<NDIM;nu++)   //  E---F---C
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
		
		su3_copy(smear_conf[ivol][mu],prop_link);
	      }
	NISSA_PARALLEL_LOOP_END;
	
	if(istep!=nSteps)
	  conf=smear_conf;
	else
	  smear_conf.invalidateHalo();
      }
  }
}
