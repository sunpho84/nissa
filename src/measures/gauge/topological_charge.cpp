#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_mix.hpp"
#include "io/input.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/float_128.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "operations/fft.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/remap_vector.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/smearing/Wflow.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "topological_charge.hpp"

namespace nissa
{
  //This will calculate the six independent components of
  //              2*a^2*ig*P_{mu,nu}
  //for a single point. Note that P_{mu,nu} is still not anti-symmetric
  //please ensure to have communicated the edges outside!
  /*
    ^                   C--<-- B --<--Y
    |                   |  2  | |  1  |
    n                   |     | |     |
    u                   D-->--\X/-->--A
    |                   D--<--/X\--<--A
    -----mu---->        |  3  | |  4  |
    .                   |     | |     |
    .                   E-->-- F -->--G
    in order to have the anti-symmetric part, use
    the routine inside "clover_term"
  */
  CUDA_HOST_AND_DEVICE void four_leaves_point(as2t_su3 leaves_summ,quad_su3 *conf,int X)
  {
    if(!check_edges_valid(conf[0])) crash("communicate edges externally");
    
    for(int mu=0;mu<NDIM;mu++)
      {
	int A=loclx_neighup[X][mu];
	int D=loclx_neighdw[X][mu];
        
	for(int nu=mu+1;nu<NDIM;nu++)
	  {
	    int munu=edge_numb[mu][nu];
	    
	    int B=loclx_neighup[X][nu];
	    int F=loclx_neighdw[X][nu];
            
	    int C=loclx_neighup[D][nu];
	    int E=loclx_neighdw[D][nu];
            
	    int G=loclx_neighdw[A][nu];
            
	    su3 temp1,temp2;
	    
	    //Leaf 1
	    unsafe_su3_prod_su3(temp1,conf[X][mu],conf[A][nu]);           //    B--<--Y
	    unsafe_su3_prod_su3_dag(temp2,temp1,conf[B][mu]);             //    |  1  |
	    unsafe_su3_prod_su3_dag(leaves_summ[munu],temp2,conf[X][nu]); //    |     |
	    /*                                                 */         //    X-->--A
            
	    //Leaf 2
	    unsafe_su3_prod_su3_dag(temp1,conf[X][nu],conf[C][mu]);       //    C--<--B
	    unsafe_su3_prod_su3_dag(temp2,temp1,conf[D][nu]);             //    |  2  |
	    unsafe_su3_prod_su3(temp1,temp2,conf[D][mu]);                 //    |     |
	    su3_summ(leaves_summ[munu],leaves_summ[munu],temp1);          //    D-->--X
            
	    //Leaf 3
	    unsafe_su3_dag_prod_su3_dag(temp1,conf[D][mu],conf[E][nu]);    //   D--<--X
	    unsafe_su3_prod_su3(temp2,temp1,conf[E][mu]);                  //   |  3  |
	    unsafe_su3_prod_su3(temp1,temp2,conf[F][nu]);                  //   |     |
	    su3_summ(leaves_summ[munu],leaves_summ[munu],temp1);           //   E-->--F
            
	    //Leaf 4
	    unsafe_su3_dag_prod_su3(temp1,conf[F][nu],conf[F][mu]);         //  X--<--A
	    unsafe_su3_prod_su3(temp2,temp1,conf[G][nu]);                   //  |  4  |
	    unsafe_su3_prod_su3_dag(temp1,temp2,conf[X][mu]);               //  |     |
	    su3_summ(leaves_summ[munu],leaves_summ[munu],temp1);            //  F-->--G
            
	    munu++;
	  }
      }
  }
  THREADABLE_FUNCTION_2ARG(four_leaves, as2t_su3*,leaves_summ, quad_su3*,conf)
  {
    GET_THREAD_ID();
    communicate_lx_quad_su3_edges(conf);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      four_leaves_point(leaves_summ[ivol],conf,ivol);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(leaves_summ);
  }
  THREADABLE_FUNCTION_END
  
  //measure the topological charge site by site
  THREADABLE_FUNCTION_2ARG(local_topological_charge, double*,charge, quad_su3*,conf)
  {
    double norm_fact=1/(128*M_PI*M_PI);
    
    as2t_su3 *leaves=nissa_malloc("leaves",loc_vol,as2t_su3);
    
    vector_reset(charge);
    
    //compute the clover-shape paths
    four_leaves(leaves,conf);
    
    //list the three combinations of plans
    int plan_id[3][2]={{0,5},{1,4},{2,3}};
    
    //loop on the three different combinations of plans
    GET_THREAD_ID();
    for(int iperm=0;iperm<3;iperm++)
      {
	//take the index of the two plans
	int ip0=plan_id[iperm][0];
	int ip1=plan_id[iperm][1];
	
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  {
	    const int sign[3]={1,-1,1};
	    
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
	  }
	NISSA_PARALLEL_LOOP_END;
      }
    
    set_borders_invalid(charge);
    
    nissa_free(leaves);
  }
  THREADABLE_FUNCTION_END
  
  //total topological charge
  THREADABLE_FUNCTION_2ARG(total_topological_charge_lx_conf, double*,tot_charge, quad_su3*,conf)
  {
    double *charge=nissa_malloc("charge",loc_vol,double);
    local_topological_charge(charge,conf);
    double_vector_glb_collapse(tot_charge,charge,loc_vol);
    nissa_free(charge);
  }
  THREADABLE_FUNCTION_END
  
  //wrapper for eos case
  THREADABLE_FUNCTION_2ARG(total_topological_charge_eo_conf, double*,tot_charge, eo_ptr<quad_su3>,eo_conf)
  {
    //convert to lx
    quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    
    total_topological_charge_lx_conf(tot_charge,lx_conf);
    
    nissa_free(lx_conf);
  }
  THREADABLE_FUNCTION_END
  
  //compute the correlator between topological charge
  THREADABLE_FUNCTION_1ARG(compute_topo_corr, double*,charge)
  {
    GET_THREAD_ID();
    
    //pass to complex
    complex *ccharge=nissa_malloc("ccharge",loc_vol,complex);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      complex_put_to_real(ccharge[ivol],charge[ivol]);
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    //transform
    fft4d(ccharge,ccharge,all_dirs,1/*complex per site*/,+1,true/*normalize*/);
    
    //multiply to build correlators
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      safe_complex_prod(ccharge[ivol],ccharge[ivol],ccharge[ivol]);
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    //transform back
    fft4d(ccharge,ccharge,all_dirs,1/*complex per site*/,-1,false/*do not normalize*/);
    
    //return to double
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      charge[ivol]=ccharge[ivol][RE];
    NISSA_PARALLEL_LOOP_END;
    nissa_free(ccharge);
  }
  THREADABLE_FUNCTION_END
  
  //finding the index to put only 1/16 of the data
  int index_to_topo_corr_remapping(int iloc_lx)
  {
    int subcube=0,subcube_el=0;
    int subcube_size[NDIM][2],subcube_coord[NDIM],subcube_el_coord[NDIM];
    for(int mu=0;mu<NDIM;mu++)
      {
	subcube_size[mu][0]=glb_size[mu]/2+1;
	subcube_size[mu][1]=glb_size[mu]/2-1;
	
	//take global coord and identify subcube
	int glx_mu=glb_coord_of_loclx[iloc_lx][mu];
	subcube_coord[mu]=(glx_mu>=subcube_size[mu][0]);
	subcube=subcube*2+subcube_coord[mu];
	
	//identify also the local coord
	subcube_el_coord[mu]=glx_mu-subcube_coord[mu]*subcube_size[mu][0];
	subcube_el=subcube_el*subcube_size[mu][subcube_coord[mu]]+subcube_el_coord[mu];
      }
    
    //summ the smaller-index cubes
    coords nsubcubes_per_dir;
    for(int mu=0;mu<NDIM;mu++) nsubcubes_per_dir[mu]=2;
    int minind_cube_vol=0;
    for(int isubcube=0;isubcube<subcube;isubcube++)
      {
	//get coords
	coords c;
	coord_of_lx(c,isubcube,nsubcubes_per_dir);
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
    irank=iglb/loc_vol;
    iloc=iglb%loc_vol;
  }
  
  //store only 1/16 of the file
  void store_topo_corr(FILE *file,double *corr,int itraj,double top,vector_remap_t *topo_corr_rem)
  {
    if(IS_PARALLEL) crash("cannot work threaded!");
    
    //remap
    topo_corr_rem->remap(corr,corr,sizeof(double));
    
    //change endianness to little
    if(!little_endian)
      {
	change_endianness((int*)&itraj,(int*)&itraj,1);
	change_endianness(corr,corr,loc_vol);
	change_endianness(&top,&top,1);
      }
    
    //offset to mantain 16 byte alignement
    if(fseek(file,3*sizeof(int),SEEK_CUR)) crash("seeking to align");
    MPI_Barrier(MPI_COMM_WORLD);
    
    //write conf id and polyakov
    if(rank==0)
      {
	off_t nwr=fwrite(&itraj,sizeof(int),1,file);
	if(nwr!=1) crash("wrote %d int instead of 1",nwr);
	nwr=fwrite(&top,sizeof(double),1,file);
	if(nwr!=1) crash("wrote %d doubles instead of 1",nwr);
      }
    else
      if(fseek(file,sizeof(int)+sizeof(double),SEEK_CUR)) crash("seeking");
    MPI_Barrier(MPI_COMM_WORLD);
    
    //find which piece has to write data
    int64_t tot_data=1;
    for(int mu=0;mu<NDIM;mu++) tot_data*=glb_size[mu]/2+1;
    
    //fix possible exceding boundary
    int64_t istart=std::min(tot_data,loc_vol*rank);
    int64_t iend=std::min(tot_data,istart+loc_vol);
    int64_t loc_data=iend-istart;
    
    //take original position of the file
    off_t ori=ftell(file);
    
    //jump to the correct point in the file
    if(fseek(file,ori+istart*sizeof(double),SEEK_SET)) crash("seeking");
    MPI_Barrier(MPI_COMM_WORLD);
    
    //write if something has to be written
    if(loc_data!=0)
      {
	int nbytes_to_write=loc_data*sizeof(double);
	off_t nbytes_wrote=fwrite(corr,1,nbytes_to_write,file);
	if(nbytes_wrote!=nbytes_to_write) crash("wrote %d bytes instead of %d",nbytes_wrote,nbytes_to_write);
      }
    
    //point to after the data
    fseek(file,ori+tot_data*sizeof(double),SEEK_SET);
  }
  
  //measure the topological charge
  void measure_topology_lx_conf(top_meas_pars_t &pars,quad_su3 *unsmoothed_conf,int iconf,bool conf_created,bool preserve_unsmoothed)
  {
    //open the file and allocate remapper
    FILE *file=open_file(pars.path,conf_created?"w":"a"),*corr_file=NULL;
    vector_remap_t *topo_corr_rem=NULL;
    if(pars.meas_corr)
      {
	corr_file=fopen(pars.corr_path.c_str(),(conf_created or !file_exists(pars.corr_path))?"w":"r+");
	if(corr_file==NULL) crash("opening %s",pars.corr_path.c_str());
	if(fseek(corr_file,0,SEEK_END)) crash("seeking to the end");
	topo_corr_rem=new vector_remap_t(loc_vol,index_to_topo_corr_remapping,NULL);
      }
    
    //allocate a temorary conf to be smoothed
    double *charge=nissa_malloc("charge",loc_vol,double);
    quad_su3 *smoothed_conf;
    if(preserve_unsmoothed)
      {
	smoothed_conf=nissa_malloc("smoothed_conf",loc_vol+bord_vol+edge_vol,quad_su3);
	vector_copy(smoothed_conf,unsmoothed_conf);
      }
    else smoothed_conf=unsmoothed_conf;
    
    int nsmooth=0;
    bool finished;
    do
      {
	//plaquette and local charge
	double plaq=global_plaquette_lx_conf(smoothed_conf);
	local_topological_charge(charge,smoothed_conf);
	//total charge
	double tot_charge;
	double_vector_glb_collapse(&tot_charge,charge,loc_vol);
	total_topological_charge_lx_conf(&tot_charge,smoothed_conf);
	master_fprintf(file,"%d %d %+16.16lg %16.16lg\n",iconf,nsmooth,tot_charge,plaq);
	finished=smooth_lx_conf_until_next_meas(smoothed_conf,pars.smooth_pars,nsmooth);
	
	//correlators if asked
	if(pars.meas_corr)
	  {
	    compute_topo_corr(charge);
	    store_topo_corr(corr_file,charge,iconf,tot_charge,topo_corr_rem);
	  }
      }
    while(!finished);
    
    //discard smoothed conf
    if(preserve_unsmoothed) nissa_free(smoothed_conf);
    nissa_free(charge);
    
    close_file(file);
    if(pars.meas_corr)
      {
	fclose(corr_file);
	delete topo_corr_rem;
      }
  }
  
  void measure_topology_eo_conf(top_meas_pars_t &pars,eo_ptr<quad_su3> unsmoothed_conf_eo,int iconf,bool conf_created)
  {
    quad_su3 *unsmoothed_conf_lx=nissa_malloc("unsmoothed_conf_lx",loc_vol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_vector(unsmoothed_conf_lx,unsmoothed_conf_eo);
    measure_topology_lx_conf(pars,unsmoothed_conf_lx,iconf,conf_created,false);
    nissa_free(unsmoothed_conf_lx);
  }
  
  //compute the topological staples site by site
  THREADABLE_FUNCTION_2ARG(topological_staples, quad_su3*,staples, quad_su3*,conf)
  {
    GET_THREAD_ID();
    as2t_su3 *leaves=nissa_malloc("leaves",loc_vol+bord_vol+edge_vol,as2t_su3);
    
    //compute the clover-shape paths
    four_leaves(leaves,conf);
    //takes the anti-symmetric part (apart from a factor 2), in an horrendous way
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int imunu=0;imunu<6;imunu++)
	{
	  color *u=leaves[ivol][imunu];
	  for(int ic1=0;ic1<NCOL;ic1++)
	    for(int ic2=ic1;ic2<NCOL;ic2++)
	      { //do not look here please, it is better to put a carpet on this uglyness
		u[ic2][ic1][0]=-(u[ic1][ic2][0]=u[ic1][ic2][0]-u[ic2][ic1][0]);
		u[ic2][ic1][1]=+(u[ic1][ic2][1]=u[ic1][ic2][1]+u[ic2][ic1][1]);
	      }
	}
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    set_borders_invalid(leaves);
    communicate_lx_as2t_su3_edges(leaves);
    
    //loop on the three different combinations of plans
    vector_reset(staples);
    NISSA_PARALLEL_LOOP(A,0,loc_vol)
      {
	//list the plan and coefficients for each staples
	const int plan_perp[4][3]={{ 5, 4, 3},{ 5, 2, 1},{ 4, 2, 0},{ 3, 1, 0}};
	const int plan_sign[4][3]={{+1,-1,+1},{-1,+1,-1},{+1,-1,+1},{-1,+1,-1}};
	
	for(int mu=0;mu<NDIM;mu++) //link direction
	  for(int inu=0;inu<NDIM-1;inu++)              //  E---F---C
	    {                                          //  |   |   | mu
	      int nu=perp_dir[mu][inu];                //  D---A---B
	      //this gives the other pair element      //        nu
	      int iplan=plan_perp[mu][inu];
	      
	      //takes neighbours
	      int B=loclx_neighup[A][nu];
	      int C=loclx_neighup[B][mu];
	      int D=loclx_neighdw[A][nu];
	      int E=loclx_neighup[D][mu];
	      int F=loclx_neighup[A][mu];
	      
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
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(staples);
    
    nissa_free(leaves);
  }
  THREADABLE_FUNCTION_END
  
  //store the topological charge if needed
  void topotential_pars_t::store_if_needed(eo_ptr<quad_su3> ext_conf,int iconf)
  {
    if(flag==2 and iconf%each==0 and iconf>=after)
      {
	double charge;
	eo_ptr<quad_su3> conf;
	if(stout_pars.nlevels==0)
	  {
	    conf[0]=ext_conf[0];
	    conf[1]=ext_conf[1];
	  }
	else
	  {
	    conf[0]=nissa_malloc("stout_conf_e",loc_volh+bord_volh+edge_volh,quad_su3);
	    conf[1]=nissa_malloc("stout_conf_o",loc_volh+bord_volh+edge_volh,quad_su3);
	    stout_smear(conf,ext_conf,&stout_pars);
	  }
	
	//compute topocharge
	total_topological_charge_eo_conf(&charge,conf);
	master_printf("Topological charge to be stored: %lg\n",charge);
	update(iconf,charge);
	
	//free if needed
	if(stout_pars.nlevels!=0)
	  {
	    nissa_free(conf[0]);
	    nissa_free(conf[1]);
	  }
      }
  }
  
  //print pars
  std::string top_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasTop\n";
    if(each!=def_each() or full) os<<" Each\t\t=\t"<<each<<"\n";
    if(after!=def_after() or full) os<<" After\t\t=\t"<<after<<"\n";
    if(path!=def_path() or full) os<<" Path\t\t=\t\""<<path.c_str()<<"\"\n";
    if(meas_corr!=def_meas_corr() or full) os<<" MeasCorr\t=\t"<<meas_corr<<"\n";
    if(corr_path!=def_corr_path() or full) os<<" CorrPath\t=\t\""<<corr_path<<"\"\n";
    os<<smooth_pars.get_str(full);
    
    return os.str();
  }
}
