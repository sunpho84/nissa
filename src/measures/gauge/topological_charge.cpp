#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include <base/vectors.hpp>
#include <communicate/borders.hpp>
#include <communicate/edges.hpp>
#include <geometry/geometry_mix.hpp>
#include <io/input.hpp>
#include <linalgs/linalgs.hpp>
#include <linalgs/reduce.hpp>
#include <measures/gauge/topological_charge.hpp>
#include <new_types/complex.hpp>
#include <new_types/float_128.hpp>
#include <new_types/spin.hpp>
#include <new_types/su3.hpp>
#include <operations/fft.hpp>
#include <operations/gaugeconf.hpp>
#include <operations/remap_vector.hpp>
#include <operations/su3_paths/gauge_sweeper.hpp>
#include <operations/smearing/stout.hpp>
#include <operations/smearing/Wflow.hpp>
#include <operations/su3_paths/plaquette.hpp>
#include <routines/ios.hpp>
#include <routines/mpi_routines.hpp>

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
  CUDA_HOST_DEVICE void four_leaves_point(as2t_su3 leaves_summ,quad_su3 *conf,const LocLxSite& X)
  {
    //if(!check_edges_valid(conf[0])) crash("communicate edges externally");
    
    FOR_ALL_DIRECTIONS(mu)
      {
	const LocLxSite& A=loclxNeighup(X,mu);
	const LocLxSite& D=loclxNeighdw(X,mu);
        
	for(Direction nu=mu+1;nu<NDIM;nu++)
	  {
	    const int munu=edge_numb[mu.nastyConvert()][nu.nastyConvert()];
	    
	    const LocLxSite& B=loclxNeighup(X,nu);
	    const LocLxSite& F=loclxNeighdw(X,nu);
            
	    const LocLxSite& C=loclxNeighup(D,nu);
	    const LocLxSite& E=loclxNeighdw(D,nu);
            
	    const LocLxSite& G=loclxNeighdw(A,nu);
            
	    su3 temp1,temp2;
	    
	    //Leaf 1
	    unsafe_su3_prod_su3(temp1,conf[X.nastyConvert()][mu.nastyConvert()],conf[A.nastyConvert()][nu.nastyConvert()]);           //    B--<--Y
	    unsafe_su3_prod_su3_dag(temp2,temp1,conf[B.nastyConvert()][mu.nastyConvert()]);             //    |  1  |
	    unsafe_su3_prod_su3_dag(leaves_summ[munu],temp2,conf[X.nastyConvert()][nu.nastyConvert()]); //    |     |
	    /*                                                 */         //    X-->--A
            
	    //Leaf 2
	    unsafe_su3_prod_su3_dag(temp1,conf[X.nastyConvert()][nu.nastyConvert()],conf[C.nastyConvert()][mu.nastyConvert()]);       //    C--<--B
	    unsafe_su3_prod_su3_dag(temp2,temp1,conf[D.nastyConvert()][nu.nastyConvert()]);             //    |  2  |
	    unsafe_su3_prod_su3(temp1,temp2,conf[D.nastyConvert()][mu.nastyConvert()]);                 //    |     |
	    su3_summ(leaves_summ[munu],leaves_summ[munu],temp1);          //    D-->--X
            
	    //Leaf 3
	    unsafe_su3_dag_prod_su3_dag(temp1,conf[D.nastyConvert()][mu.nastyConvert()],conf[E.nastyConvert()][nu.nastyConvert()]);    //   D--<--X
	    unsafe_su3_prod_su3(temp2,temp1,conf[E.nastyConvert()][mu.nastyConvert()]);                  //   |  3  |
	    unsafe_su3_prod_su3(temp1,temp2,conf[F.nastyConvert()][nu.nastyConvert()]);                  //   |     |
	    su3_summ(leaves_summ[munu],leaves_summ[munu],temp1);           //   E-->--F
            
	    //Leaf 4
	    unsafe_su3_dag_prod_su3(temp1,conf[F.nastyConvert()][nu.nastyConvert()],conf[F.nastyConvert()][mu.nastyConvert()]);         //  X--<--A
	    unsafe_su3_prod_su3(temp2,temp1,conf[G.nastyConvert()][nu.nastyConvert()]);                   //  |  4  |
	    unsafe_su3_prod_su3_dag(temp1,temp2,conf[X.nastyConvert()][mu.nastyConvert()]);               //  |     |
	    su3_summ(leaves_summ[munu],leaves_summ[munu],temp1);            //  F-->--G
	  }
      }
  }
  void four_leaves(as2t_su3* leaves_summ,quad_su3* conf)
  {
    communicate_lx_quad_su3_edges(conf);
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      four_leaves_point(leaves_summ[ivol.nastyConvert()],conf,ivol);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(leaves_summ);
  }
  
  //measure the topological charge site by site
  void local_topological_charge(double* charge,quad_su3* conf)
  {
    double norm_fact=1/(128*M_PI*M_PI);
    
    as2t_su3 *leaves=nissa_malloc("leaves",locVol.nastyConvert(),as2t_su3);
    
    vector_reset(charge);
    
    //compute the clover-shape paths
    four_leaves(leaves,conf);
    
    //list the three combinations of plans
    int plan_id[3][2]={{0,5},{1,4},{2,3}};
    
    //loop on the three different combinations of plans
    for(int iperm=0;iperm<3;iperm++)
      {
	//take the index of the two plans
	int ip0=plan_id[iperm][0];
	int ip1=plan_id[iperm][1];
	
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  {
	    const int sign[3]={1,-1,1};
	    
	    //products
	    su3 clock,aclock;
	    unsafe_su3_prod_su3_dag(clock,leaves[ivol.nastyConvert()][ip0],leaves[ivol.nastyConvert()][ip1]);
	    unsafe_su3_prod_su3(aclock,leaves[ivol.nastyConvert()][ip0],leaves[ivol.nastyConvert()][ip1]);
	    
	    //take the trace
	    complex tclock,taclock;
	    su3_trace(tclock,clock);
	    su3_trace(taclock,aclock);
	    
	    //takes the combination with appropriate sign
	    charge[ivol.nastyConvert()]+=sign[iperm]*(tclock[RE]-taclock[RE])*norm_fact;
	  }
	NISSA_PARALLEL_LOOP_END;
      }
    
    set_borders_invalid(charge);
    
    nissa_free(leaves);
  }
  
  //total topological charge
  void total_topological_charge_lx_conf(double* tot_charge,quad_su3* conf)
  {
    double *charge=nissa_malloc("charge",locVol.nastyConvert(),double);
    local_topological_charge(charge,conf);
    glb_reduce(tot_charge,charge,locVol.nastyConvert());
    nissa_free(charge);
  }
  
  //wrapper for eos case
  void total_topological_charge_eo_conf(double* tot_charge,eo_ptr<quad_su3> eo_conf)
  {
    //convert to lx
    quad_su3 *lx_conf=nissa_malloc("lx_conf",locVolWithBordAndEdge.nastyConvert(),quad_su3);
    paste_eo_parts_into_lx_vector(lx_conf,eo_conf);
    
    total_topological_charge_lx_conf(tot_charge,lx_conf);
    
    nissa_free(lx_conf);
  }
  
  //compute the correlator between topological charge
  void compute_topo_corr(double* charge)
  {
    
    //pass to complex
    complex *ccharge=nissa_malloc("ccharge",locVol.nastyConvert(),complex);
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      complex_put_to_real(ccharge[ivol.nastyConvert()],charge[ivol.nastyConvert()]);
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    //transform
    fft4d(ccharge,ccharge,all_dirs,1/*complex per site*/,+1,true/*normalize*/);
    
    //multiply to build correlators
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      safe_complex_prod(ccharge[ivol.nastyConvert()],ccharge[ivol.nastyConvert()],ccharge[ivol.nastyConvert()]);
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    //transform back
    fft4d(ccharge,ccharge,all_dirs,1/*complex per site*/,-1,false/*do not normalize*/);
    
    //return to double
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      charge[ivol.nastyConvert()]=ccharge[ivol.nastyConvert()][RE];
    NISSA_PARALLEL_LOOP_END;
    nissa_free(ccharge);
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
    irank=(iglb/locVol).nastyConvert();
    iloc=(iglb%locVol).nastyConvert();
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
	change_endianness(corr,corr,locVol.nastyConvert());
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
    for(int mu=0;mu<NDIM;mu++) tot_data*=glbSize[mu]/2+1;
    
    //fix possible exceding boundary
    int64_t istart=std::min(tot_data,(locVol*rank).nastyConvert());
    int64_t iend=std::min(tot_data,(istart+locVol).nastyConvert());
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
	topo_corr_rem=new vector_remap_t(locVol(),index_to_topo_corr_remapping,NULL);
      }
    
    //allocate a temorary conf to be smoothed
    double *charge=nissa_malloc("charge",locVol.nastyConvert(),double);
    quad_su3 *smoothed_conf;
    if(preserve_unsmoothed)
      {
	smoothed_conf=nissa_malloc("smoothed_conf",locVolWithBordAndEdge.nastyConvert(),quad_su3);
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
	glb_reduce(&tot_charge,charge,locVol.nastyConvert());
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
    quad_su3 *unsmoothed_conf_lx=nissa_malloc("unsmoothed_conf_lx",locVolWithBordAndEdge.nastyConvert(),quad_su3);
    paste_eo_parts_into_lx_vector(unsmoothed_conf_lx,unsmoothed_conf_eo);
    measure_topology_lx_conf(pars,unsmoothed_conf_lx,iconf,conf_created,false);
    nissa_free(unsmoothed_conf_lx);
  }
  
  //compute the topological staples site by site
  void topological_staples(quad_su3* staples,quad_su3* conf)
  {
    as2t_su3 *leaves=nissa_malloc("leaves",locVolWithBordAndEdge.nastyConvert(),as2t_su3);
    
    //compute the clover-shape paths
    four_leaves(leaves,conf);
    //takes the anti-symmetric part (apart from a factor 2), in an horrendous way
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int imunu=0;imunu<6;imunu++)
	{
	  color *u=leaves[ivol.nastyConvert()][imunu];
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
    NISSA_PARALLEL_LOOP(A,0,locVol)
      {
	//list the plan and coefficients for each staples
	const int plan_perp[4][3]={{ 5, 4, 3},{ 5, 2, 1},{ 4, 2, 0},{ 3, 1, 0}};
	const int plan_sign[4][3]={{+1,-1,+1},{-1,+1,-1},{+1,-1,+1},{-1,+1,-1}};
	
	FOR_ALL_DIRECTIONS(mu)                         //link direction
	  for(int inu=0;inu<NDIM-1;inu++)              //  E---F---C
	    {                                          //  |   |   | mu
	      const Direction nu=perp_dir[mu.nastyConvert()][inu];                //  D---A---B //nasty
	      //this gives the other pair element      //        nu
	      int iplan=plan_perp[mu.nastyConvert()][inu];
	      
	      //takes neighbours
	      const LocLxSite& B=loclxNeighup(A,nu);
	      const LocLxSite& C=loclxNeighup(B,mu);
	      const LocLxSite& D=loclxNeighdw(A,nu);
	      const LocLxSite& E=loclxNeighup(D,mu);
	      const LocLxSite& F=loclxNeighup(A,mu);
	      
	      //compute ABC, BCF and the full SU(3) staple, ABCF
	      su3 ABC,BCF,ABCF;
	      unsafe_su3_prod_su3(ABC,conf[A.nastyConvert()][nu.nastyConvert()],conf[B.nastyConvert()][mu.nastyConvert()]);
	      unsafe_su3_prod_su3_dag(BCF,conf[B.nastyConvert()][mu.nastyConvert()],conf[F.nastyConvert()][nu.nastyConvert()]);
	      unsafe_su3_prod_su3_dag(ABCF,ABC,conf[F.nastyConvert()][nu.nastyConvert()]);
	      
	      //compute ADE, DEF and the full SU(3) staple, ADEF
	      su3 ADE,DEF,ADEF;
	      unsafe_su3_dag_prod_su3(ADE,conf[D.nastyConvert()][nu.nastyConvert()],conf[D.nastyConvert()][mu.nastyConvert()]);
	      unsafe_su3_prod_su3(DEF,conf[D.nastyConvert()][mu.nastyConvert()],conf[E.nastyConvert()][nu.nastyConvert()]);
	      unsafe_su3_prod_su3(ADEF,ADE,conf[E.nastyConvert()][nu.nastyConvert()]);
	      
	      //local summ and temp
	      su3 loc_staples,temp;
	      su3_put_to_zero(loc_staples);
	      //insert the leave in the four possible forward positions
	      
	      unsafe_su3_prod_su3(loc_staples,leaves[A.nastyConvert()][iplan],ABCF);     //insertion on A
	      unsafe_su3_prod_su3(temp,conf[A.nastyConvert()][nu.nastyConvert()],leaves[B.nastyConvert()][iplan]);
	      su3_summ_the_prod_su3(loc_staples,temp,BCF);                //insertion on B
	      unsafe_su3_prod_su3(temp,ABC,leaves[C.nastyConvert()][iplan]);
	      su3_summ_the_prod_su3_dag(loc_staples,temp,conf[F.nastyConvert()][nu.nastyConvert()]);    //insertion on C
	      su3_summ_the_prod_su3(loc_staples,ABCF,leaves[F.nastyConvert()][iplan]);   //insertion on F
	      
	      //insert the leave in the four possible backward positions
	      su3_summ_the_dag_prod_su3(loc_staples,leaves[A.nastyConvert()][iplan],ADEF);    //insertion on A
	      unsafe_su3_dag_prod_su3_dag(temp,conf[D.nastyConvert()][nu.nastyConvert()],leaves[D.nastyConvert()][iplan]);
	      su3_summ_the_prod_su3(loc_staples,temp,DEF);                     //insertion on D
	      unsafe_su3_prod_su3_dag(temp,ADE,leaves[E.nastyConvert()][iplan]);
	      su3_summ_the_prod_su3(loc_staples,temp,conf[E.nastyConvert()][nu.nastyConvert()]);             //insertion on E
	      su3_summ_the_prod_su3_dag(loc_staples,ADEF,leaves[F.nastyConvert()][iplan]);    //insertion on F
	      
	      //summ or subtract, according to the coefficient
	      if(plan_sign[mu.nastyConvert()][inu]==+1)
		su3_summassign(staples[A.nastyConvert()][mu.nastyConvert()],loc_staples);
	      else
		su3_subtassign(staples[A.nastyConvert()][mu.nastyConvert()],loc_staples);
	    }
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(staples);
    
    nissa_free(leaves);
  }
  
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
	    conf[0]=nissa_malloc("stout_conf_e",(locVolh+bord_volh+edge_volh).nastyConvert(),quad_su3);
	    conf[1]=nissa_malloc("stout_conf_o",(locVolh+bord_volh+edge_volh).nastyConvert(),quad_su3);
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
