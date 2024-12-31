#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <algorithm>

#include "geometry/geometry_mix.hpp"
#include "io/endianness.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "operations/remap_vector.hpp"
#include "routines/ios.hpp"
#include "topological_charge.hpp"

namespace nissa
{
  /// Computes the four leaves on all sites
  void four_leaves(LxField<as2t_su3>& leavesSumm,
		   const LxField<quad_su3>& conf)
  {
    conf.updateEdges();
    
    PAR(0,locVol,
	CAPTURE(TO_READ(conf),
		TO_WRITE(leavesSumm)),
	ivol,
	{
	  four_leaves_point(leavesSumm[ivol],conf,ivol);
	});
  }
  
  /// Measure the topological charge site by site
  void local_topological_charge(LxField<double>& charge,
				const LxField<quad_su3>& conf)
  {
    const double norm_fact=1/(128*M_PI*M_PI);
    
    LxField<as2t_su3> leaves("leaves");
    
    charge.reset();
    
    //compute the clover-shape paths
    four_leaves(leaves,conf);
    
    //list the three combinations of plans
    constexpr int plan_id[3][2]={{0,5},{1,4},{2,3}};
    
    //loop on the three different combinations of plans
    for(int iperm=0;iperm<3;iperm++)
      {
	//take the index of the two plans
	const int ip0=plan_id[iperm][0];
	const int ip1=plan_id[iperm][1];
	
	PAR(0,locVol,
	    CAPTURE(iperm,
		    norm_fact,
		    ip0,ip1,
		    TO_READ(conf),
		    TO_READ(leaves),
		    TO_WRITE(charge)),
	    ivol,
	    {
	      constexpr int sign[3]={1,-1,1};
	     
	     //products
	     su3 clock,aclock;
	     unsafe_su3_prod_su3_dag(clock,leaves[ivol][ip0],leaves[ivol][ip1]);
	     unsafe_su3_prod_su3(aclock,leaves[ivol][ip0],leaves[ivol][ip1]);
	     
	     //take the trace
	     complex tclock,taclock;
	     su3_trace(tclock,clock);
	     su3_trace(taclock,aclock);
	     
	     //takes the combination with appropriate sign
	     charge[ivol]+=sign[iperm]*(tclock[RE]-taclock[RE])*norm_fact;
	    });
      }
    
    charge.invalidateHalo();
  }
  
  /// total topological charge
  double total_topological_charge_lx_conf(const LxField<quad_su3>& conf)
  {
    LxField<double> charge("charge");
    local_topological_charge(charge,conf);
    
    double totCharge;
    glb_reduce(&totCharge,charge,locVol);
    
    return totCharge;
  }
  
  /// Wrapper for eos case
  double total_topological_charge_eo_conf(const EoField<quad_su3>& eoConf)
  {
    //convert to lx
    LxField<quad_su3> lxConf("lx_conf",WITH_HALO);
    paste_eo_parts_into_lx_vector(lxConf,eoConf);
    
    return total_topological_charge_lx_conf(lxConf);
  }
  
  //compute the correlator between topological charge
  void compute_topo_corr(double* charge)
  {
    CRASH("reimplement");
    // //pass to complex
    // complex *ccharge=nissa_malloc("ccharge",locVol,complex);
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   complex_put_to_real(ccharge[ivol],charge[ivol]);
    // NISSA_PARALLEL_LOOP_END;
    // THREAD_BARRIER();
    
    // //transform
    // fft4d(ccharge,ccharge,all_dirs,1/*complex per site*/,+1,true/*normalize*/);
    
    // //multiply to build correlators
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   safe_complex_prod(ccharge[ivol],ccharge[ivol],ccharge[ivol]);
    // NISSA_PARALLEL_LOOP_END;
    // THREAD_BARRIER();
    
    // //transform back
    // fft4d(ccharge,ccharge,all_dirs,1/*complex per site*/,-1,false/*do not normalize*/);
    
    // //return to double
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   charge[ivol]=ccharge[ivol][RE];
    // NISSA_PARALLEL_LOOP_END;
    // nissa_free(ccharge);
  }
  
  //finding the index to put only 1/16 of the data
  int index_to_topo_corr_remapping(int iloc_lx)
  {
    int subcube=0,subcube_el=0;
    int subcube_size[NDIM][2],subcube_coord[NDIM],subcube_el_coord[NDIM];
    for(int mu=0;mu<NDIM;mu++)
      {
	subcube_size[mu][0]=glbSize[mu]/2+1;
	subcube_size[mu][1]=glbSize[mu]/2-1;
	
	//take global coord and identify subcube
	int glx_mu=glbCoordOfLoclx[iloc_lx][mu];
	subcube_coord[mu]=(glx_mu>=subcube_size[mu][0]);
	subcube=subcube*2+subcube_coord[mu];
	
	//identify also the local coord
	subcube_el_coord[mu]=glx_mu-subcube_coord[mu]*subcube_size[mu][0];
	subcube_el=subcube_el*subcube_size[mu][subcube_coord[mu]]+subcube_el_coord[mu];
      }
    
    //summ the smaller-index cubes
    Coords nsubcubes_per_dir;
    for(int mu=0;mu<NDIM;mu++) nsubcubes_per_dir[mu]=2;
    int minind_cube_vol=0;
    for(int isubcube=0;isubcube<subcube;isubcube++)
      {
	//get coords
	Coords c=coordOfLx(isubcube,nsubcubes_per_dir);
	//compute vol
	int subcube_vol=1;
	for(int mu=0;mu<NDIM;mu++) subcube_vol*=subcube_size[mu][c[mu]];
	minind_cube_vol+=subcube_vol;
      }
    
    return subcube_el+minind_cube_vol;
  }
  
  //wrapper
  void index_to_topo_corr_remapping(int &irank,int &iloc,int iloc_lx,void *pars)
  {
    int iglb=index_to_topo_corr_remapping(iloc_lx);
    
    //find rank and loclx
    irank=iglb/locVol;
    iloc=iglb%locVol;
  }
  
  //store only 1/16 of the file
  void store_topo_corr(FILE *file,double *corr,int itraj,double top,vector_remap_t *topo_corr_rem)
  {
    //remap
    topo_corr_rem->remap(corr,corr,sizeof(double));
    
    //change endianness to little
    if(not LittleEndian)
      {
	CRASH("reimplement");
	// change_endianness((int*)&itraj,(int*)&itraj,1);
	// change_endianness(corr,corr,locVol);
	// change_endianness(&top,&top,1);
      }
    
    //offset to mantain 16 byte alignement
    if(fseek(file,3*sizeof(int),SEEK_CUR)) CRASH("seeking to align");
    MPI_Barrier(MPI_COMM_WORLD);
    
    //write conf id and polyakov
    if(rank==0)
      {
	off_t nwr=fwrite(&itraj,sizeof(int),1,file);
	if(nwr!=1) CRASH("wrote %ld int instead of 1",nwr);
	nwr=fwrite(&top,sizeof(double),1,file);
	if(nwr!=1) CRASH("wrote %ld doubles instead of 1",nwr);
      }
    else
      if(fseek(file,sizeof(int)+sizeof(double),SEEK_CUR)) CRASH("seeking");
    MPI_Barrier(MPI_COMM_WORLD);
    
    //find which piece has to write data
    int64_t tot_data=1;
    for(int mu=0;mu<NDIM;mu++) tot_data*=glbSize[mu]/2+1;
    
    //fix possible exceding boundary
    int64_t istart=std::min(tot_data,(int64_t)locVol*rank);
    int64_t iend=std::min(tot_data,istart+locVol);
    int64_t loc_data=iend-istart;
    
    //take original position of the file
    off_t ori=ftell(file);
    
    //jump to the correct point in the file
    if(fseek(file,ori+istart*sizeof(double),SEEK_SET)) CRASH("seeking");
    MPI_Barrier(MPI_COMM_WORLD);
    
    //write if something has to be written
    if(loc_data!=0)
      {
	int nbytes_to_write=loc_data*sizeof(double);
	off_t nbytes_wrote=fwrite(corr,1,nbytes_to_write,file);
	if(nbytes_wrote!=nbytes_to_write) CRASH("wrote %ld bytes instead of %d",nbytes_wrote,nbytes_to_write);
      }
    
    //point to after the data
    fseek(file,ori+tot_data*sizeof(double),SEEK_SET);
  }
  
  //measure the topological charge
  void measure_topology_lx_conf(const top_meas_pars_t &pars,
				const LxField<quad_su3>& unsmoothed_conf,
				const int& iconf,
				const bool& conf_created,
				const bool& preserve_unsmoothed)
  {
    CRASH("reimplement");
    
    // //open the file and allocate remapper
    // FILE *file=open_file(pars.path,conf_created?"w":"a"),*corr_file=NULL;
    // vector_remap_t *topo_corr_rem=NULL;
    // if(pars.meas_corr)
    //   {
    // 	corr_file=fopen(pars.corr_path.c_str(),(conf_created or !file_exists(pars.corr_path))?"w":"r+");
    // 	if(corr_file==NULL) CRASH("opening %s",pars.corr_path.c_str());
    // 	if(fseek(corr_file,0,SEEK_END)) CRASH("seeking to the end");
    // 	topo_corr_rem=new vector_remap_t(locVol,index_to_topo_corr_remapping,NULL);
    //   }
    
    // //allocate a temorary conf to be smoothed
    // LxField<double> charge("charge");
    // LxField<quad_su3> smoothed_conf("smoothed_conf",WITH_HALO_EDGES);
    // smoothed_conf=unsmoothed_conf;
    
    // int nsmooth=0;
    // bool finished;
    // do
    //   {
    // 	//plaquette and local charge
    // 	const double plaq=global_plaquette_lx_conf(smoothed_conf);
    // 	local_topological_charge(charge,smoothed_conf);
	
    // 	//total charge
    // 	double tot_charge;
    // 	glb_reduce(&tot_charge,charge,locVol);
    // 	total_topological_charge_lx_conf(&tot_charge,smoothed_conf);
    // 	master_fprintf(file,"%d %d %+16.16lg %16.16lg\n",iconf,nsmooth,tot_charge,plaq);
    // 	finished=smooth_lx_conf_until_next_meas(smoothed_conf,pars.smooth_pars,nsmooth);
	
    // 	//correlators if asked
    // 	if(pars.meas_corr)
    // 	  {
    // 	    CRASH("reimplement");
    // 	    // compute_topo_corr(charge);
    // 	    // store_topo_corr(corr_file,charge,iconf,tot_charge,topo_corr_rem);
    // 	  }
    //   }
    // while(not finished);
    
    // //discard smoothed conf
    // if(preserve_unsmoothed) nissa_free(smoothed_conf);
    // nissa_free(charge);
    
    // close_file(file);
    // if(pars.meas_corr)
    //   {
    // 	fclose(corr_file);
    // 	delete topo_corr_rem;
    //   }
  }
  
  void measure_topology_eo_conf(const top_meas_pars_t &pars,
				const EoField<quad_su3>& unsmoothed_conf_eo,
				const int& iconf,
				const bool& conf_created)
  {
	    CRASH("reimplement");
    // quad_su3 *unsmoothed_conf_lx=nissa_malloc("unsmoothed_conf_lx",locVol+bord_vol+edge_vol,quad_su3);
    // paste_eo_parts_into_lx_vector(unsmoothed_conf_lx,unsmoothed_conf_eo);
    // measure_topology_lx_conf(pars,unsmoothed_conf_lx,iconf,conf_created,false);
    // nissa_free(unsmoothed_conf_lx);
  }
  
  //compute the topological staples site by site
  void topological_staples(LxField<quad_su3>& staples,
			   const LxField<quad_su3>& conf)
  {
    LxField<as2t_su3> leaves("leaves",WITH_HALO_EDGES);
    
    //compute the clover-shape paths
    four_leaves(leaves,conf);
    
    //takes the anti-symmetric part (apart from a factor 2), in an horrendous way
    PAR(0,locVol,
	CAPTURE(TO_WRITE(leaves)),
	ivol,
	{
	  for(int imunu=0;imunu<6;imunu++)
	    {
	      auto u=leaves[ivol][imunu];
	      for(int ic1=0;ic1<NCOL;ic1++)
		for(int ic2=ic1;ic2<NCOL;ic2++)
		  { //do not look here please, it is better to put a carpet on this uglyness
		    u[ic2][ic1][0]=-(u[ic1][ic2][0]=u[ic1][ic2][0]-u[ic2][ic1][0]);
		    u[ic2][ic1][1]=+(u[ic1][ic2][1]=u[ic1][ic2][1]+u[ic2][ic1][1]);
		  }
	    }
	});
    
    leaves.invalidateHalo();
    leaves.updateEdges();
    
    //loop on the three different combinations of plans
    staples.reset();
    PAR(0,locVol,
	CAPTURE(TO_READ(conf),
		TO_READ(leaves),
		TO_WRITE(staples)),
	A,
	{
	  //list the plan and coefficients for each staples
	  const int plan_perp[4][3]={{ 5, 4, 3},{ 5, 2, 1},{ 4, 2, 0},{ 3, 1, 0}};
	  const int plan_sign[4][3]={{+1,-1,+1},{-1,+1,-1},{+1,-1,+1},{-1,+1,-1}};
	  
	  for(int mu=0;mu<NDIM;mu++) //link direction
	    for(int inu=0;inu<NDIM-1;inu++)              //  E---F---C
	      {                                          //  |   |   | mu
		int nu=perpDirs[mu][inu];                //  D---A---B
		//this gives the other pair element      //        nu
		int iplan=plan_perp[mu][inu];
		
		//takes neighbours
		int B=loclxNeighup[A][nu];
		int C=loclxNeighup[B][mu];
		int D=loclxNeighdw[A][nu];
		int E=loclxNeighup[D][mu];
		int F=loclxNeighup[A][mu];
		
		//compute ABC, BCF and the full SU(3) staple, ABCF
		su3 ABC,BCF,ABCF;
		unsafe_su3_prod_su3(ABC,conf[A][nu],conf[B][mu]);
		unsafe_su3_prod_su3_dag(BCF,conf[B][mu],conf[F][nu]);
		unsafe_su3_prod_su3_dag(ABCF,ABC,conf[F][nu]);
		
		//compute ADE, DEF and the full SU(3) staple, ADEF
		su3 ADE,DEF,ADEF;
		unsafe_su3_dag_prod_su3(ADE,conf[D][nu],conf[D][mu]);
		unsafe_su3_prod_su3(DEF,conf[D][mu],conf[E][nu]);
		unsafe_su3_prod_su3(ADEF,ADE,conf[E][nu]);
		
		//local summ and temp
		su3 loc_staples,temp;
		su3_put_to_zero(loc_staples);
		//insert the leave in the four possible forward positions
		
		unsafe_su3_prod_su3(loc_staples,leaves[A][iplan],ABCF);     //insertion on A
		unsafe_su3_prod_su3(temp,conf[A][nu],leaves[B][iplan]);
		su3_summ_the_prod_su3(loc_staples,temp,BCF);                //insertion on B
		unsafe_su3_prod_su3(temp,ABC,leaves[C][iplan]);
		su3_summ_the_prod_su3_dag(loc_staples,temp,conf[F][nu]);    //insertion on C
		su3_summ_the_prod_su3(loc_staples,ABCF,leaves[F][iplan]);   //insertion on F
	      
		//insert the leave in the four possible backward positions
		su3_summ_the_dag_prod_su3(loc_staples,leaves[A][iplan],ADEF);    //insertion on A
		unsafe_su3_dag_prod_su3_dag(temp,conf[D][nu],leaves[D][iplan]);
		su3_summ_the_prod_su3(loc_staples,temp,DEF);                     //insertion on D
		unsafe_su3_prod_su3_dag(temp,ADE,leaves[E][iplan]);
		su3_summ_the_prod_su3(loc_staples,temp,conf[E][nu]);             //insertion on E
		su3_summ_the_prod_su3_dag(loc_staples,ADEF,leaves[F][iplan]);    //insertion on F
		
		//summ or subtract, according to the coefficient
		if(plan_sign[mu][inu]==+1) su3_summassign(staples[A][mu],loc_staples);
		else                       su3_subtassign(staples[A][mu],loc_staples);
	      }
	});
    
    staples.invalidateHalo();
  }
}
