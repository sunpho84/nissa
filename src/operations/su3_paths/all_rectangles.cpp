#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <map>
#include <vector>

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "operations/shift.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "smearing/APE.hpp"
#include "smearing/HYP.hpp"

namespace nissa
{
  //compute all possible rectangular paths among a defined interval
  THREADABLE_FUNCTION_4ARG(measure_all_rectangular_paths_new, all_rect_meas_pars_t*,pars, quad_su3*,ori_conf, int,iconf, int,create_output_file)
  {
    GET_THREAD_ID();
    
    FILE *fout=NULL;
    if(rank==0 && IS_MASTER_THREAD) fout=open_file(pars->path,create_output_file?"w":"a");
    
    //running conf
    quad_su3 *sme_conf=nissa_malloc("sme_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    
    //compute competing volume
    int prp_vol=(int)ceil((double)glb_size[1]*glb_size[1]/nranks);
    int min_vol=prp_vol*rank,max_vol=min_vol+prp_vol;
    if(max_vol>glb_size[1]*glb_size[1]) max_vol=glb_size[1]*glb_size[1];
    int cmp_vol=max_vol-min_vol;
    
    //saved spatial links for competing volume, for the three dirs and different smearing level, 
    //and for all length up to Dmax
    int nape_spat_levls=pars->nape_spat_levls;
    int Dmax=pars->Dmax;
    //su3 *spat_links=nissa_malloc("spat_links",3*cmp_vol*nape_spat_levls*Dmax,su3);
    
    //hyp or temporal APE smear the conf
    if(pars->use_hyp_or_ape_temp==0) hyp_smear_conf_dir(sme_conf,ori_conf,pars->hyp_temp_alpha0,pars->hyp_temp_alpha1,pars->hyp_temp_alpha2,0);
    else ape_temporal_smear_conf(sme_conf,ori_conf,pars->ape_temp_alpha,pars->nape_temp_iters);
    
    /*
    //list of destinations for each link and dir
    std::map<int,std::vector<int> > out_per_rank[3];
    
    //first find where to scatter data
    for(int ii=0;ii<3;ii++)
      {
	int i=ii+1;
	int j=perp2_dir[0][i][0],k=perp2_dir[0][i][1];
	for(int ivol=0;ivol<loc_vol;ivol++)
	  {
	    //find dest in the global indexing
	    int glb_prp_dest_site=glb_coord_of_loclx[ivol][k]+glb_size[k]*glb_coord_of_loclx[ivol][j];
	    int dest_rank=glb_prp_dest_site/prp_vol;
	    int loc_prp_dest_site=glb_prp_dest_site-dest_rank*prp_vol;
	    
	    //find size on the remote rank
	    int dest_prp_vol=prp_vol;
	    if(dest_rank==nranks-1) dest_prp_vol=glb_size[1]*glb_size[1]-dest_rank*prp_vol;
	    
	    //find local site in the new layout
	    int dest_site=glb_coord_of_loclx[ivol][i]+glb_size[i]*
	      (glb_coord_of_loclx[ivol][0]+glb_size[0]*(loc_prp_dest_site+dest_prp_vol*ii));
	    
	    //mark it in the list
	    out_per_rank[ii][dest_rank].push_back(dest_site);
	  }
      }
    
    //send info to corresponding rank
    int *nin_per_rank=nissa_malloc("nin_per_rank",nranks*3,int);
    int *nout_per_rank=nissa_malloc("nout_per_rank",nranks*3,int);
    for(int delta_rank=0;delta_rank<nranks;delta_rank++)
      {
	int dest_rank=(rank+delta_rank)%nranks;
	int recv_rank=(rank+nranks-delta_rank)%nranks;
	for(int ii=0;ii<3;ii++)
	  {
	    //check if is different from 0
	    int temp_nout_per_rank=0;
	    std::map<int,std::vector<int> >::iterator it=out_per_rank[ii].find(dest_rank);
	    if(it!=out_per_rank[ii].end()) temp_nout_per_rank=it->second.size();
	    nout_per_rank[ii+3*dest_rank]=temp_nout_per_rank;
	  }
	MPI_Sendrecv(nout_per_rank+dest_rank*3,3,MPI_INT,dest_rank,0,
		     nin_per_rank+recv_rank*3 ,3,MPI_INT,recv_rank,0,
		     MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    
    for(int dest_rank=0;dest_rank<nranks;dest_rank++)
      master_printf("Send to %d: %d %d %d\n",dest_rank,nout_per_rank[dest_rank*3+0],nout_per_rank[dest_rank*3+1],
		    nout_per_rank[dest_rank*3+2]);
    
    printf("%d\n",928341);
    THREAD_BARRIER();
    */
    /*
    //buffer out and in
    su3 *buffer_out=nissa_mallo("buffer_out",comp_vol,su3);
    su3 *buffer_in=nissa_mallo("buffer_in",comp_vol,su3);
    
    //loop over spatial APE smeared levels
    for(int iape=0;iape<pars->nape_spat_levls;iape++)
      {
        //spatial APE smearing
        ape_spatial_smear_conf(sme_conf,sme_conf,pars->ape_spat_alpha,
	       (iape==0)?pars->nape_spat_iters[0]:(pars->nape_spat_iters[iape]-pars->nape_spat_iters[iape-1]));
	
	//scatter all links
	for(int i=1;i<=3;i++)
	  {
	    
	  }
      }
    */
    
    //printf("%d\n",928341);
    nissa_free(sme_conf);
    //printf("%d\n",928342);
    if(rank==0 && IS_MASTER_THREAD) fclose(fout);
    //printf("%d\n",928343);
  }}

  //compute all possible rectangular paths among a defined interval
  THREADABLE_FUNCTION_4ARG(measure_all_rectangular_paths, all_rect_meas_pars_t*,pars, quad_su3*,ori_conf, int,iconf, int,create_output_file)
  {
    GET_THREAD_ID();
    
    FILE *fout=NULL;
    if(rank==0 && IS_MASTER_THREAD) fout=open_file(pars->path,create_output_file?"w":"a");
    
    //hypped and APE spatial smeared conf
    quad_su3 *sme_conf=nissa_malloc("sme_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    
    //allocate time path, time-spatial paths, closing paths and point contribution
    su3 *T_path=nissa_malloc("T_path",loc_vol+bord_vol,su3);
    su3 *TS_path=nissa_malloc("TS_path",loc_vol+bord_vol,su3);
    su3 *closed_path=nissa_malloc("closed_path",loc_vol+bord_vol,su3);
    double *point_path=nissa_malloc("point_path",loc_vol,double);
    
    //hyp or temporal APE smear the conf
    if(pars->use_hyp_or_ape_temp==0) hyp_smear_conf_dir(sme_conf,ori_conf,pars->hyp_temp_alpha0,pars->hyp_temp_alpha1,pars->hyp_temp_alpha2,0);
    else ape_temporal_smear_conf(sme_conf,ori_conf,pars->ape_temp_alpha,pars->nape_temp_iters);
    
    //loop over APE smeared levels
    for(int iape=0;iape<pars->nape_spat_levls;iape++)
      {
	//APE smearing
	ape_spatial_smear_conf(sme_conf,sme_conf,pars->ape_spat_alpha,
		  (iape==0)?pars->nape_spat_iters[0]:(pars->nape_spat_iters[iape]-pars->nape_spat_iters[iape-1]));
	
	//reset the Tpath link product
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  su3_put_to_id(T_path[ivol]);
	
	//move along T up to Tmax
	for(int t=1;t<=pars->Tmax;t++)
	  {
	    //take the product
	    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	      safe_su3_prod_su3(T_path[ivol],T_path[ivol],sme_conf[ivol][0]);
	    set_borders_invalid(T_path);
	    
	    //push up the vector along T
	    su3_vec_single_shift(T_path,0,+1);
	    
	    //results to be printed, averaged along the three dirs
	    double paths[pars->Dmax+1][3];
	    for(int ii=0;ii<3;ii++) for(int d=0;d<=pars->Dmax;d++) paths[d][ii]=0;
	    
	    //if T_path is long enough we move along spatial dirs
	    if(t>=pars->Tmin)
	      for(int ii=0;ii<3;ii++)
		{
		  int i=ii+1;
		  
		  //copy T_path
		  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
		    su3_copy(TS_path[ivol],T_path[ivol]);
		  
		  //move along i up to Dmax
		  for(int d=1;d<=pars->Dmax;d++)
		    {
		      //take the product
		      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
			safe_su3_prod_su3(TS_path[ivol],TS_path[ivol],sme_conf[ivol][i]);
		      set_borders_invalid(TS_path);
		      
		      //push up the vector along i
		      su3_vec_single_shift(TS_path,i,+1);
		      
		      //if TS_path is long enough we close the path
		      if(d>=pars->Dmin)
			{
			  //copy TS_path
			  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
			    su3_copy(closed_path[ivol],TS_path[ivol]);
			  set_borders_invalid(closed_path);
			  
			  //move back along time
			  for(int tp=1;tp<=t;tp++)
			    {
			      //push dw the vector along 0
			      su3_vec_single_shift(closed_path,0,-1);
			      
			      //take the product with dag
			      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
				safe_su3_prod_su3_dag(closed_path[ivol],closed_path[ivol],sme_conf[ivol][0]);
			      set_borders_invalid(closed_path);
			    }
			  
			  //move back along space
			  for(int dp=1;dp<=d;dp++)
			    {
			      //push dw the vector along i
			      su3_vec_single_shift(closed_path,i,-1);
			      
			      //take the product with dag
			      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
				safe_su3_prod_su3_dag(closed_path[ivol],closed_path[ivol],sme_conf[ivol][i]);
			      set_borders_invalid(closed_path);
			    }
			  
			  //take the trace and store it in the point contribution
			  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
			    point_path[ivol]=
			    closed_path[ivol][0][0][RE]+closed_path[ivol][1][1][RE]+closed_path[ivol][2][2][RE];
			  
			  //reduce among all threads and ranks and summ it
			  double temp;
			  double_vector_glb_collapse(&temp,point_path,loc_vol);
			  paths[d][ii]+=temp;
			}
		    }
		}
	    
	    //print all the Dmax contributions, with ncol*nspat_dir*glb_vol normalization
	    if(t>=pars->Tmin && rank==0 && IS_MASTER_THREAD)
	      for(int d=pars->Dmin;d<=pars->Dmax;d++)
		{
		  fprintf(fout,"%d %d  %d %d",iconf,iape,t,d);
		  for(int ii=0;ii<3;ii++) fprintf(fout,"\t%16.16lg",paths[d][ii]/(3*glb_vol));
		  fprintf(fout,"\n");
		}
	    THREAD_BARRIER();
	  }
      }
    
    //close file
    if(rank==0 && IS_MASTER_THREAD) fclose(fout);
    
    //free stuff
    nissa_free(sme_conf);
    nissa_free(T_path);
    nissa_free(TS_path);
    nissa_free(closed_path);
    nissa_free(point_path);
  }}
  
  void measure_all_rectangular_paths(all_rect_meas_pars_t *pars,quad_su3 **conf_eo,int iconf,int create_output_file)
  {
    quad_su3 *conf_lx=nissa_malloc("conf_lx",loc_vol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_conf(conf_lx,conf_eo);
    
    measure_all_rectangular_paths(pars,conf_lx,iconf,create_output_file);
    
    nissa_free(conf_lx);
  }
}
