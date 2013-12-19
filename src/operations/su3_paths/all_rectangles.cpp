#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <map>
#include <vector>
#include <stdlib.h>

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
#include "su3_paths/plaquette.hpp"

namespace nissa
{
  typedef int tricoords_t[3];
  
  //return the three coords of site in the transposed space
  void get_tricoords_of_site(tricoords_t c,int icmp,tricoords_t L)
  {
    for(int mu=2;mu>=0;mu--)
      {
        c[mu]=icmp%L[mu];
        icmp/=L[mu];
      }
  }

  //return the three coords of site in the transposed space
  int get_site_of_tricoords(tricoords_t c,tricoords_t L)
  {
    int icmp=0;
    
    for(int mu=0;mu<3;mu++)
      icmp=icmp*L[mu]+c[mu];
    
    return icmp;
  }
  
  //copy
  void tricoords_copy(tricoords_t out,tricoords_t in)
  {
    out[0]=in[0];
    out[1]=in[1];
    out[2]=in[2];
  }
  
  //shift in a single dir
  int site_shift(int icmp,tricoords_t L,int mu,int shift)
  {
    tricoords_t in;
    get_tricoords_of_site(in,icmp,L);
    
    in[mu]=in[mu]+shift;
    while(in[mu]<0) in[mu]+=L[mu];
    while(in[mu]>=L[mu]) in[mu]-=L[mu];
    
    return get_site_of_tricoords(in,L);
  }
  
  //compute the transposed from lx index
  void index_transp(int &irank_transp,int &iloc_transp,int iloc_lx,void *pars)
  {
    int ii=((int*)pars)[0],prp_vol=((int*)pars)[1];
    int i=ii+1;
    
    //directions perpendicular to 0 and i
    int j=perp2_dir[0][ii][0],k=perp2_dir[0][ii][1];

    //find dest in the global indexing
    int *g=glb_coord_of_loclx[iloc_lx];
    int glb_dest_site=g[i]+glb_size[i]*(g[0]+glb_size[0]*(g[k]+glb_size[k]*g[j]));
    irank_transp=glb_dest_site/prp_vol;
    iloc_transp=glb_dest_site-irank_transp*prp_vol;
  }

  //compute all possible rectangular paths among a defined interval
  THREADABLE_FUNCTION_4ARG(measure_all_rectangular_paths_new, all_rect_meas_pars_t*,pars, quad_su3*,ori_conf, int,iconf, int,create_output_file)
  {
    GET_THREAD_ID();
    
    //running conf
    quad_su3 *sme_conf=nissa_malloc("sme_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    
    //compute competing volume
    int prp_vol=glb_size[0]*glb_size[1]*((int)ceil((double)glb_size[1]*glb_size[1]/nranks));
    int min_vol=prp_vol*rank,max_vol=min_vol+prp_vol;
    if(max_vol>=glb_vol) max_vol=glb_vol;
    int cmp_vol=max_vol-min_vol;
    
    //define the three remapper
    vector_remap_t *remap[3];
    for(int ii=0;ii<3;ii++)
      {
	int pars[2]={ii,prp_vol};
	remap[ii]=new vector_remap_t(loc_vol,index_transp,pars);
	if(remap[ii]->nel_in!=cmp_vol) crash("expected %d obtained %d",cmp_vol,remap[ii]->nel_in);
      }
    
    //transposed configurations
    //we need 3 copies, each holding 1 smeared temporal links and 3*nape_spat_levls spatial links per site
    int nape_spat_levls=pars->nape_spat_levls;
    int ntot_sme=1+nape_spat_levls;
    su3 *transp_conf=nissa_malloc("transp_confs",3*cmp_vol*ntot_sme,su3);
    //local conf holders pre-transposing
    su3 *pre_transp_conf_holder=nissa_malloc("pre_transp_conf_holder",loc_vol,su3);
    su3 *post_transp_conf_holder=nissa_malloc("post_transp_conf_holder",cmp_vol,su3);

    //hyp or temporal APE smear the conf
    if(pars->use_hyp_or_ape_temp==0) hyp_smear_conf_dir(sme_conf,ori_conf,pars->hyp_temp_alpha0,
							pars->hyp_temp_alpha1,pars->hyp_temp_alpha2,0);
    else ape_temporal_smear_conf(sme_conf,ori_conf,pars->ape_temp_alpha,pars->nape_temp_iters);
    master_printf("Plaquette after temp smear: %16.16lg\n",global_plaquette_lx_conf(sme_conf));
    
    //store temporal links and send them
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      su3_copy(pre_transp_conf_holder[ivol],sme_conf[ivol][0]);
    THREAD_BARRIER();
    for(int ii=0;ii<3;ii++)
      {
	remap[ii]->remap(post_transp_conf_holder,pre_transp_conf_holder,sizeof(su3));
	NISSA_PARALLEL_LOOP(icmp,0,cmp_vol)
	  su3_copy(transp_conf[icmp+cmp_vol*(0+ntot_sme*ii)],post_transp_conf_holder[icmp]);
	THREAD_BARRIER();
      }
    
    //spatial APE smearing
    for(int iape=0;iape<nape_spat_levls;iape++)
      {
        ape_spatial_smear_conf(sme_conf,sme_conf,pars->ape_spat_alpha,
	       (iape==0)?pars->nape_spat_iters[0]:(pars->nape_spat_iters[iape]-pars->nape_spat_iters[iape-1]));
	master_printf("Plaquette after %d ape smears: %16.16lg\n",iape+1,global_plaquette_lx_conf(sme_conf));
	
	//store spatial links and send them
	for(int ii=0;ii<3;ii++)
	  {
	    int i=ii+1;
	    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	      su3_copy(pre_transp_conf_holder[ivol],sme_conf[ivol][i]);
	    THREAD_BARRIER();

	    remap[ii]->remap(post_transp_conf_holder,pre_transp_conf_holder,sizeof(su3));
	    NISSA_PARALLEL_LOOP(icmp,0,cmp_vol)
	      su3_copy(transp_conf[icmp+cmp_vol*((1+iape)+ntot_sme*ii)],post_transp_conf_holder[icmp]);
	    THREAD_BARRIER();
	  }
      }

    //free smeared conf, pre-post buffers and remappers
    nissa_free(post_transp_conf_holder);
    nissa_free(pre_transp_conf_holder);
    nissa_free(sme_conf);
    for(int ii=0;ii<3;ii++) delete remap[ii];
    
    ////////////////////////////////////////////////////////////////////////////////////
    
    //all the rectangles, for each thread
    int dD=pars->Dmax+1-pars->Dmin;
    int dT=pars->Tmax+1-pars->Tmin;
    double *all_rectangles=nissa_malloc("all_rectangles",dD*dT*3*nape_spat_levls*nthreads,double);
    vector_reset(all_rectangles);
    double *all_rectangles_loc_thread=all_rectangles+dD*dT*3*nape_spat_levls*THREAD_ID;
    
    //all time-lines for all distances dT, and running space-lines
    su3 *Tline=nissa_malloc("Tline",cmp_vol*dT,su3);
    su3 *Dline=nissa_malloc("Dline",cmp_vol,su3);
    
    for(int ii=0;ii<3;ii++)
      {
	int i=ii+1;
	tricoords_t L={cmp_vol/glb_size[0]/glb_size[i],glb_size[0],glb_size[i]};
	
	//create all Tline
	NISSA_PARALLEL_LOOP(icmp,0,cmp_vol)
	  {
	    su3 U;
	    su3_copy(U,transp_conf[icmp+cmp_vol*(0+ntot_sme*ii)]);
	    
	    for(int t=1;t<pars->Tmin;t++)
	      safe_su3_prod_su3(U,U,transp_conf[site_shift(icmp,L,1,t)+cmp_vol*(0+ntot_sme*ii)]);
	    for(int dt=0;dt<dT;dt++)
	      {
		su3_copy(Tline[icmp*dT+dt],U);
		safe_su3_prod_su3(U,U,transp_conf[site_shift(icmp,L,1,dt+pars->Tmin)+cmp_vol*(0+ntot_sme*ii)]);
	      }
	  }
	THREAD_BARRIER();
	
	for(int iape=0;iape<nape_spat_levls;iape++)
	  {
	    //create Dlines up to Dmin
	    NISSA_PARALLEL_LOOP(icmp,0,cmp_vol)
	      {
		su3_copy(Dline[icmp],transp_conf[icmp+cmp_vol*((1+iape)+ntot_sme*ii)]);
		
		for(int d=1;d<pars->Dmin;d++)
		  safe_su3_prod_su3(Dline[icmp],Dline[icmp],
				    transp_conf[site_shift(icmp,L,2,d)+cmp_vol*((1+iape)+ntot_sme*ii)]);
	      }
	    THREAD_BARRIER();
	    
	    for(int dd=0;dd<dD;dd++)
	      {
		int d=dd+pars->Dmin;

		//closes
		NISSA_PARALLEL_LOOP(icmp,0,cmp_vol)
		  {
		    su3 part1,part2;
		    for(int dt=0;dt<dT;dt++)
		      {
			int t=dt+pars->Tmin;
			unsafe_su3_prod_su3(part1,Dline[icmp],Tline[site_shift(icmp,L,2,d)*dT+dt]);
			unsafe_su3_prod_su3(part2,Tline[icmp*dT+dt],Dline[site_shift(icmp,L,1,t)]);
			all_rectangles_loc_thread[dd+dD*(dt+dT*(ii+3*iape))]+=
			  real_part_of_trace_su3_prod_su3_dag(part1,part2);
		      }
		  }
		
		//prolong
		NISSA_PARALLEL_LOOP(icmp,0,cmp_vol)
		  safe_su3_prod_su3(Dline[icmp],Dline[icmp],
				    transp_conf[site_shift(icmp,L,2,d)+cmp_vol*((1+iape)+ntot_sme*ii)]);
		THREAD_BARRIER();
	      }
	  }
      }
    
#ifdef USE_THREADS
    if(nthreads>1)
      if(IS_MASTER_THREAD)
	for(unsigned int other_thread=1;other_thread<nthreads;other_thread++)
	  {
	    double *all_rectangles_other_thread=all_rectangles+dD*dT*3*nape_spat_levls*other_thread;
	    for(int i=0;i<dD*dT*3*nape_spat_levls;i++)
	      all_rectangles_loc_thread[i]+=all_rectangles_other_thread[i];		
	  }
    THREAD_BARRIER();
#endif
    
    //perform all rank reduction and print
    double *all_rectangles_glb=nissa_malloc("all_rectangles",dD*dT*3*nape_spat_levls,double);    
    if(IS_MASTER_THREAD)
      {
	decript_MPI_error(MPI_Reduce(all_rectangles,all_rectangles_glb,dD*dT*3*nape_spat_levls,
				     MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD),"reducing");
	
	//open file
	if(rank==0)
	  {
	    FILE *fout=open_file(pars->path,create_output_file?"w":"a");
	    for(int iape=0;iape<nape_spat_levls;iape++)
	      for(int dt=0;dt<dT;dt++)
		for(int dd=0;dd<dD;dd++)
		  {
		    fprintf(fout,"%d %d  %d %d",iconf,iape,dt+pars->Tmin,dd+pars->Dmin);
		    for(int ii=0;ii<3;ii++) fprintf(fout,"\t%16.16lg",
						    all_rectangles_glb[dd+dD*(dt+dT*(ii+3*iape))]/(3*glb_vol));
		    fprintf(fout,"\n");
		  }
	    fclose(fout);
	  }
      }
    
    nissa_free(all_rectangles_glb);
    nissa_free(all_rectangles);
    nissa_free(Tline);
    nissa_free(Dline);
    nissa_free(transp_conf);
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
	    
	    //print all the Dmax contributions, with ncol*glb_vol normalization
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
