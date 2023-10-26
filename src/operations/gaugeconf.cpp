#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/random.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_eo.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/su3_op.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"
#include "measures/gauge/topological_charge.hpp"
#include "routines/mpi_routines.hpp"

/*
  rotate a field anti-clockwise by 90 degrees
  
   0---1---2---0        0---6---3---0
   |   |   |   |        |   |   |   |
   6---7---8---6        2---8---5---2
   |   |   |   |        |   |   |   |
   3---4---5---3        1---7---4---1
   |   |   |   |        |   |   |   |
   O---1---2---0        O---6---3---0
   
   d2
   O d1
   
   where d1=axis+1
   and   d2=d1+1
   
*/

namespace nissa
{
  void ac_rotate_vector(void *out,void *in,int axis,size_t bps)
  {
    //find the two swapping direction
    int d1=1+(axis-1+1)%3;
    int d2=1+(axis-1+2)%3;
    
    //check that the two directions have the same size and that we are not asking 0 as axis
    if(glbSizes[d1]!=glbSizes[d2]) crash("Rotation works only if dir %d and %d have the same size!",glbSizes[d1],glbSizes[d2]);
    if(axis==0) crash("Error, only spatial rotations implemented");
    int L=glbSizes[d1];
    
    //allocate destinations and sources
    coords_t *xto=nissa_malloc("xto",locVol,coords_t);
    coords_t *xfr=nissa_malloc("xfr",locVol,coords_t);
    
    //scan all local sites to see where to send and from where to expect data
    NISSA_LOC_VOL_LOOP(ivol)
    {
      //copy 0 and axis coord to "to" and "from" sites
      xto[ivol][0]=xfr[ivol][0]=glbCoordOfLoclx[ivol][0];
      xto[ivol][axis]=xfr[ivol][axis]=glbCoordOfLoclx[ivol][axis];
      
      //find reamining coord of "to" site
      xto[ivol][d1]=(L-glbCoordOfLoclx[ivol][d2])%L;
      xto[ivol][d2]=glbCoordOfLoclx[ivol][d1];
      
      //find remaining coord of "from" site
      xfr[ivol][d1]=glbCoordOfLoclx[ivol][d2];
      xfr[ivol][d2]=(L-glbCoordOfLoclx[ivol][d1])%L;
    }
    
    //call the remapping
    //remap_vector((char*)out,(char*)in,xto,xfr,bps);
    crash("to be reimplemented");
    
    //free vectors
    nissa_free(xfr);
    nissa_free(xto);
  }
  
  /*
    rotate the gauge configuration anti-clockwise by 90 degrees
    this is more complicated than a single vector because of link swaps
    therefore the rotation is accomplished through 2 separates steps
    
    .---.---.---.     .---.---.---.       .---.---.---.
    |           |     |           |       |           |
    .   B 3 C   .     .   B 2'C   .       .   C 4'D   .
    |   2   4   |     |   1   3   |       |   3   1   |
    .   A 1 D   .     .   A 4'D   .       .   B 2'A   .
    |           |     |           |       |           |
    O---.---.---.     O---.---.---.       O---.---.---.
    
    d2
    O d1
    
  */
  
  void ac_rotate_gauge_conf(quad_su3 *out,quad_su3 *in,int axis)
  {
    crash("reimplement");
    
    // int d0=0;
    // int d1=1+(axis-1+1)%3;
    // int d2=1+(axis-1+2)%3;
    // int d3=axis;
    
    // //allocate a temporary conf with borders
    // quad_su3 *temp_conf=nissa_malloc("temp_conf",locVol+bord_vol,quad_su3);
    // memcpy(temp_conf,in,locVol*sizeof(quad_su3));
    // communicate_lx_quad_su3_borders(temp_conf);
    
    // //now reorder links
    // NISSA_LOC_VOL_LOOP(ivol)
    // {
    //   //copy temporal direction and axis
    //   memcpy(out[ivol][d0],temp_conf[ivol][d0],sizeof(su3));
    //   memcpy(out[ivol][d3],temp_conf[ivol][d3],sizeof(su3));
    //   //swap the other two
    //   unsafe_su3_hermitian(out[ivol][d1],temp_conf[loclxNeighdw[ivol][d2]][d2]);
    //   memcpy(out[ivol][d2],temp_conf[ivol][d1],sizeof(su3));
    // }
    
    // //rotate rigidly
    // ac_rotate_vector(out,out,axis,sizeof(quad_su3));
  }
  
  void put_boundaries_conditions(LxField<quad_su3>& conf,
				 const momentum_t& theta_in_pi,
				 const int& putonbords,
				 const int& putonedges)
  {
    complex theta[NDIM];
    for(int idir=0;idir<NDIM;idir++)
      {
	theta[idir][0]=cos(theta_in_pi[idir]*M_PI/glbSizes[idir]);
	theta[idir][1]=sin(theta_in_pi[idir]*M_PI/glbSizes[idir]);
      }
    
    int nsite=locVol;
    if(putonbords) nsite+=bord_vol;
    if(putonedges) nsite+=edge_vol;
    
    PAR(0,nsite,
	CAPTURE(theta,
		TO_WRITE(conf)),
	ivol,
      {
	for(int idir=0;idir<NDIM;idir++)
	  safe_su3_prod_complex(conf[ivol][idir],conf[ivol][idir],theta[idir]);
      });
  }
  
  void rem_boundaries_conditions(LxField<quad_su3>& conf,
				 const momentum_t& theta_in_pi,
				 const int& putonbords,
				 const int& putonedges)
  {
    const momentum_t minus_theta_in_pi={-theta_in_pi[0],-theta_in_pi[1],-theta_in_pi[2],-theta_in_pi[3]};
    put_boundaries_conditions(conf,minus_theta_in_pi,putonbords,putonedges);
  }
  
  void adapt_theta(LxField<quad_su3>& conf,
		   momentum_t& old_theta,
		   const momentum_t& put_theta,
		   const int& putonbords,
		   const int& putonedges)
  {
    momentum_t diff_theta;
    int adapt=0;
    
    for(int idir=0;idir<NDIM;idir++)
      {
	adapt=adapt or (old_theta[idir]!=put_theta[idir]);
	diff_theta[idir]=put_theta[idir]-old_theta[idir];
	old_theta[idir]=put_theta[idir];
      }
    
    if(adapt)
      {
	master_printf("Necessary to add boundary condition: %f %f %f %f\n",diff_theta[0],diff_theta[1],diff_theta[2],diff_theta[3]);
	put_boundaries_conditions(conf,diff_theta,putonbords,putonedges);
      }
  }
  
  //generate an identical conf
  void generate_cold_eo_conf(EoField<quad_su3>& conf)
  {
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(par,
		    TO_WRITE(conf)),
	    ieo,
	    {
	      for(int mu=0;mu<NDIM;mu++)
		su3_put_to_id(conf[par][ieo][mu]);
	    });
      }
  }
  
  //generate a random conf
  void generate_hot_eo_conf(EoField<quad_su3>& conf)
  {
    if(loc_rnd_gen_inited==0) crash("random number generator not inited");
    
    for(int par=0;par<2;par++)
      {
	PAR(0,locVolh,
	    CAPTURE(par,
		    TO_WRITE(conf)),
	    ieo,
	    {
	      int ilx=loclx_of_loceo[par][ieo];
	      for(int mu=0;mu<NDIM;mu++)
		su3_put_to_rnd(conf[par][ieo][mu],loc_rnd_gen[ilx]);
	    });
      }
  }
  
  //unitarize an a lx conf
  void unitarize_lx_conf_orthonormalizing(quad_su3* conf)
  {
    crash("reimplement");
    
    // START_TIMING(unitarize_time,nunitarize);
    
    // PAR(0,locVol,
    // 	CAPTURE(TO_WRITE(conf)),
    // 	ivol,
    // 	{
    // 	  for(int idir=0;idir<NDIM;idir++)
    // 	    {
    // 	      su3 t;
    // 	      su3_unitarize_orthonormalizing(t,conf[ivol][idir]);
    // 	      su3_copy(conf[ivol][idir],t);
    // 	    }
    // 	});
    // STOP_TIMING(unitarize_time);
  }
  
  //unitarize the conf by explicitly by projecting it maximally to su3
  void unitarize_lx_conf_maximal_trace_projecting(quad_su3* conf)
  {
    crash("reimplement");
    
    // START_TIMING(unitarize_time,nunitarize);
    
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   for(int mu=0;mu<NDIM;mu++)
    //     su3_unitarize_maximal_trace_projecting(conf[ivol][mu],conf[ivol][mu]);
    // NISSA_PARALLEL_LOOP_END;
    
    // set_borders_invalid(conf);
    // STOP_TIMING(unitarize_time);
  }
  
  //eo version
  void unitarize_eo_conf_maximal_trace_projecting(eo_ptr<quad_su3> conf)
  {
    crash("reimplement");
    
    // START_TIMING(unitarize_time,nunitarize);
    
    // for(int par=0;par<2;par++)
    //   {
    //     NISSA_PARALLEL_LOOP(ivol,0,locVolh)
    //       for(int mu=0;mu<NDIM;mu++)
    //         su3_unitarize_maximal_trace_projecting(conf[par][ivol][mu],conf[par][ivol][mu]);
    // 	NISSA_PARALLEL_LOOP_END;
        
    //     set_borders_invalid(conf[par]);
    //   }
    
    // STOP_TIMING(unitarize_time);
  }
  
  //overrelax an lx configuration
  void overrelax_lx_conf_handle(su3 out,su3 staple,int ivol,int mu,void *pars)
  {
    su3_find_overrelaxed(out,out,staple,((int*)pars)[0]);
  }
  
  void overrelax_lx_conf(LxField<quad_su3>& conf,gauge_sweeper_t* sweeper,int nhits)
  {
    crash("reimplement");
    
    //(sweeper->sweep_conf(conf,overrelax_lx_conf_handle,(void*)&nhits);
  }
  
  //same for heatbath
  namespace heatbath_lx_conf_ns
  {
    struct pars_t
    {
      double beta;
      int nhits;
      pars_t(double beta,int nhits) : beta(beta),nhits(nhits){}
    };
    void handle(su3 out,su3 staple,int ivol,int mu,void *pars)
    {su3_find_heatbath(out,out,staple,((pars_t*)pars)->beta,((pars_t*)pars)->nhits,loc_rnd_gen+ivol);}
  }
  
  void heatbath_lx_conf(LxField<quad_su3>& conf,gauge_sweeper_t* sweeper,const double& beta,const int& nhits)
  {
    crash("reimplement");
    // heatbath_lx_conf_ns::pars_t pars(beta,nhits);sweeper->sweep_conf(conf,heatbath_lx_conf_ns::handle,&pars);
  }
  
  //same for cooling
  void cool_lx_conf_handle(su3 out,su3 staple,int ivol,int mu,void *pars)
  {su3_unitarize_maximal_trace_projecting_iteration(out,staple);}
  
  void cool_lx_conf(LxField<quad_su3>& conf,gauge_sweeper_t* sweeper)
  {
    sweeper->sweep_conf(conf,cool_lx_conf_handle,NULL);}
  
  //measure the average gauge energy
  double average_gauge_energy(const LxField<quad_su3>& conf)
  {
    double energy;
    
    conf.updateEdges();
    LxField<double> loc_energy("energy");
    
    PAR(0,locVol,
	CAPTURE(TO_WRITE(loc_energy),
		TO_READ(conf)),
	ivol,
      {
	//compute the clover-shape paths
	as2t_su3 leaves;
	four_leaves_point(leaves,conf,ivol);
	
	loc_energy[ivol]=0.0;
	for(int i=0;i<NDIM*(NDIM-1)/2;i++)
	  {
	    su3 A;
	    unsafe_su3_subt_su3_dag(A,leaves[i],leaves[i]);
	    su3_prodassign_double(A,1.0/8.0); //factor 1/2 for the antihermitian, 1/4 for average leave
	    complex temp;
	    trace_su3_prod_su3(temp,A,A);
	    loc_energy[ivol]-=temp[RE];
	  }
	loc_energy[ivol]/=glbVol;
      });
    
    glb_reduce(&energy,loc_energy,locVol);
    
    return energy;
  }
}
