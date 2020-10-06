#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //perform ape smearing
  //be sure not to have border condition added
  THREADABLE_FUNCTION_6ARG(ape_smear_conf, quad_su3*,smear_conf, quad_su3*,origi_conf, double,alpha, int,nstep, int*,dirs, int,min_staple_dir)
  {
    GET_THREAD_ID();
    
    quad_su3 *temp_conf=nissa_malloc("temp_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    if(origi_conf!=smear_conf) double_vector_copy((double*)smear_conf,(double*)origi_conf,loc_vol*sizeof(quad_su3)/sizeof(double));
    
    char listed_dirs[21]="";
    for(int mu=0;mu<NDIM;mu++)
      if(dirs[mu])
	{
	  char temp[3];
	  snprintf(temp,3,"%d ",mu);
	  strncat(listed_dirs,temp,20);
	}
    verbosity_lv1_master_printf("APE { %s} smearing with alpha=%g, %d iterations\n",listed_dirs,alpha,nstep);
    
    for(int istep=0;istep<nstep;istep++)
      {
	verbosity_lv3_master_printf("APE spatial smearing with alpha=%g iteration %d of %d\n",alpha,istep,nstep);
	double_vector_copy((double*)temp_conf,(double*)smear_conf,loc_vol*sizeof(quad_su3)/sizeof(double));
	
	//communicate the borders
	communicate_lx_quad_su3_edges(temp_conf);
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  for(int mu=0;mu<NDIM;mu++)
	    if(dirs[mu])
	      {
		//calculate staples
		su3 stap,temp1,temp2;
		su3_put_to_zero(stap);
		for(int nu=min_staple_dir;nu<NDIM;nu++)   //  E---F---C
		  if(nu!=mu)                              //  |   |   | mu
		    {                                     //  D---A---B
		      int A=ivol;                         //   nu
		      int B=loclx_neighup[A][nu];
		      int F=loclx_neighup[A][mu];
		      unsafe_su3_prod_su3(temp1,temp_conf[A][nu],temp_conf[B][mu]);
		      unsafe_su3_prod_su3_dag(temp2,temp1,temp_conf[F][nu]);
		      su3_summ(stap,stap,temp2);
		      
		      int D=loclx_neighdw[A][nu];
		      int E=loclx_neighup[D][mu];
		      unsafe_su3_dag_prod_su3(temp1,temp_conf[D][nu],temp_conf[D][mu]);
		      unsafe_su3_prod_su3(temp2,temp1,temp_conf[E][nu]);
		      su3_summ(stap,stap,temp2);
		    }
		
		//create new link to be reunitarized
		su3 prop_link;
		for(int icol1=0;icol1<NCOL;icol1++)
		  for(int icol2=0;icol2<NCOL;icol2++)
		    for(int ri=0;ri<2;ri++)
		      //prop_link[icol1][icol2][ri]=(1-alpha)*temp_conf[ivol][mu][icol1][icol2][ri]+alpha/6*stap[icol1][icol2][ri];
		      prop_link[icol1][icol2][ri]=temp_conf[ivol][mu][icol1][icol2][ri]+alpha*stap[icol1][icol2][ri];
		
		su3_unitarize_maximal_trace_projecting(smear_conf[ivol][mu],prop_link);
	      }
	
	set_borders_invalid(smear_conf);
      }
    
    nissa_free(temp_conf);
  }
  THREADABLE_FUNCTION_END
}
