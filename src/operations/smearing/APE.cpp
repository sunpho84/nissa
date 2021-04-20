#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //perform ape smearing
  //be sure not to have border condition added
  void ape_smear_conf(quad_su3* smear_conf,quad_su3* origi_conf,double alpha,int nstep,const Coords<bool>& dirs,const Dir& min_staple_dir)
  {
    
    quad_su3 *temp_conf=nissa_malloc("temp_conf",locVolWithBordAndEdge.nastyConvert(),quad_su3);
    if(origi_conf!=smear_conf) double_vector_copy((double*)smear_conf,(double*)origi_conf,locVol.nastyConvert()*sizeof(quad_su3)/sizeof(double));
    
    char listed_dirs[21]="";
    FOR_ALL_DIRS(mu)
      if(dirs(mu))
	{
	  char temp[3];
	  snprintf(temp,3,"%d ",mu());
	  strncat(listed_dirs,temp,20);
	}
    verbosity_lv1_master_printf("APE { %s} smearing with alpha=%g, %d iterations\n",listed_dirs,alpha,nstep);
    
    for(int istep=0;istep<nstep;istep++)
      {
	verbosity_lv3_master_printf("APE spatial smearing with alpha=%g iteration %d of %d\n",alpha,istep,nstep);
	double_vector_copy((double*)temp_conf,(double*)smear_conf,locVol.nastyConvert()*sizeof(quad_su3)/sizeof(double));
	
	//communicate the borders
	communicate_lx_quad_su3_edges(temp_conf);
	
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  FOR_ALL_DIRS(mu)
	    if(dirs(mu))
	      {
		//calculate staples
		su3 stap,temp1,temp2;
		su3_put_to_zero(stap);
		for(Dir nu=min_staple_dir;nu<NDIM;nu++)    //  E---F---C
		  if(nu!=mu)                                     //  |   |   | mu
		    {                                            //  D---A---B
		      const LocLxSite& A=ivol;                   //   nu
		      const LocLxSite& B=loclxNeighup(A,nu);
		      const LocLxSite& F=loclxNeighup(A,mu);
		      unsafe_su3_prod_su3(temp1,temp_conf[A.nastyConvert()][nu.nastyConvert()],temp_conf[B.nastyConvert()][mu.nastyConvert()]);
		      unsafe_su3_prod_su3_dag(temp2,temp1,temp_conf[F.nastyConvert()][nu.nastyConvert()]);
		      su3_summ(stap,stap,temp2);
		      
		      const LocLxSite& D=loclxNeighdw(A,nu);
		      const LocLxSite& E=loclxNeighup(D,mu);
		      unsafe_su3_dag_prod_su3(temp1,temp_conf[D.nastyConvert()][nu.nastyConvert()],temp_conf[D.nastyConvert()][mu.nastyConvert()]);
		      unsafe_su3_prod_su3(temp2,temp1,temp_conf[E.nastyConvert()][nu.nastyConvert()]);
		      su3_summ(stap,stap,temp2);
		    }
		
		//create new link to be reunitarized
		su3 prop_link;
		for(int icol1=0;icol1<NCOL;icol1++)
		  for(int icol2=0;icol2<NCOL;icol2++)
		    for(int ri=0;ri<2;ri++)
		      //prop_link[icol1][icol2][ri]=(1-alpha)*temp_conf[ivol.nastyConvert()][mu][icol1][icol2][ri]+alpha/6*stap[icol1][icol2][ri];
		      prop_link[icol1][icol2][ri]=temp_conf[ivol.nastyConvert()][mu.nastyConvert()][icol1][icol2][ri]+alpha*stap[icol1][icol2][ri];
		
		su3_unitarize_maximal_trace_projecting(smear_conf[ivol.nastyConvert()][mu.nastyConvert()],prop_link);
	      }
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(smear_conf);
      }
    
    nissa_free(temp_conf);
  }
}
