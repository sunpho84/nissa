#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "operations/shift.hpp"
#include "operations/remap_vector.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "operations/smearing/APE.hpp"
#include "operations/smearing/HYP.hpp"
#include "operations/su3_paths/plaquette.hpp"

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
    int mu0=((int*)pars)[0],imu1=((int*)pars)[1],prp_vol=((int*)pars)[2];
    int mu1=perp_dir[mu0][imu1],mu2=perp2_dir[mu0][imu1][0],mu3=perp2_dir[mu0][imu1][1];

    //find dest in the global indexing
    int *g=glb_coord_of_loclx[iloc_lx];
    int glb_dest_site=g[mu1]+glb_size[mu1]*(g[mu0]+glb_size[mu0]*(g[mu2]+glb_size[mu2]*g[mu3]));
    irank_transp=glb_dest_site/prp_vol;
    iloc_transp=glb_dest_site-irank_transp*prp_vol;
  }

  //compute all possible rectangular paths among a defined interval
  THREADABLE_FUNCTION_4ARG(measure_all_rectangular_paths, all_rect_meas_pars_t*,pars, quad_su3*,ori_conf, int,iconf, int,create_output_file)
  {
    GET_THREAD_ID();
    
    //remapping
    int nape_spat_levls=pars->nape_spat_levls,ntot_sme=1+nape_spat_levls;
    int prp_vol[12],cmp_vol[12],imu01=0,mu0_l[12],mu1_l[12],cmp_vol_max=0;
    vector_remap_t *remap[12];
    su3 *transp_conf[12];
    su3 *pre_transp_conf_holder=nissa_malloc("pre_transp_conf_holder",loc_vol,su3);
    for(int mu0=0;mu0<4;mu0++)
      for(int imu1=0;imu1<3;imu1++)
	{
	  //find dirs
	  int mu1=perp_dir[mu0][imu1],mu2=perp2_dir[mu0][imu1][0],mu3=perp2_dir[mu0][imu1][1];
	  mu0_l[imu01]=mu0;mu1_l[imu01]=mu1;
	  
	  //compute competing volume
	  prp_vol[imu01]=glb_size[mu0]*glb_size[mu1]*((int)ceil((double)glb_size[mu2]*glb_size[mu3]/nranks));
	  int min_vol=prp_vol[imu01]*rank,max_vol=min_vol+prp_vol[imu01];
	  if(max_vol>=glb_vol) max_vol=glb_vol;
	  cmp_vol[imu01]=max_vol-min_vol;
	  cmp_vol_max=std::max(cmp_vol_max,cmp_vol[imu01]);
	  
	  //define the six remapper
	  int pars[3]={mu0,imu1,prp_vol[imu01]};
	  remap[imu01]=new vector_remap_t(loc_vol,index_transp,pars);
	  if(remap[imu01]->nel_in!=cmp_vol[imu01]) crash("expected %d obtained %d",cmp_vol[imu01],remap[imu01]->nel_in);
	  
	  //allocate transp conf
	  transp_conf[imu01]=nissa_malloc("transp_conf",cmp_vol[imu01]*ntot_sme,su3);
	  
	  imu01++;
	}
    
    //local conf holders post-transposing
    su3 *post_transp_conf_holder=nissa_malloc("post_transp_conf_holder",cmp_vol_max,su3);
	  
    //hyp or APE smear the conf
    quad_su3 *sme_conf=nissa_malloc("sme_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    for(int mu0=0;mu0<4;mu0++)
      {
	if(pars->use_hyp_or_ape_temp==0) hyp_smear_conf_dir(sme_conf,ori_conf,pars->hyp_temp_alpha0,
							    pars->hyp_temp_alpha1,pars->hyp_temp_alpha2,mu0);
	else ape_single_dir_smear_conf(sme_conf,ori_conf,pars->ape_temp_alpha,pars->nape_temp_iters,mu0);
	verbosity_lv1_master_printf("Plaquette after \"temp\" (%d) smear: %16.16lg\n",mu0,global_plaquette_lx_conf(sme_conf));
	
	//store temporal links and send them
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  su3_copy(pre_transp_conf_holder[ivol],sme_conf[ivol][mu0]);
	THREAD_BARRIER();	
	for(int imu1=0;imu1<3;imu1++)
	  {
	    int imu01=mu0*3+imu1;
	    remap[imu01]->remap(post_transp_conf_holder,pre_transp_conf_holder,sizeof(su3));
	    NISSA_PARALLEL_LOOP(icmp,0,cmp_vol[imu01])
	      su3_copy(transp_conf[imu01][icmp+cmp_vol[imu01]*0],post_transp_conf_holder[icmp]);
	    THREAD_BARRIER();
	  }

	//spatial APE smearing
	for(int iape=0;iape<nape_spat_levls;iape++)
	  {
	    ape_perp_dir_smear_conf(sme_conf,sme_conf,pars->ape_spat_alpha,
		(iape==0)?pars->nape_spat_iters[0]:(pars->nape_spat_iters[iape]-pars->nape_spat_iters[iape-1]),mu0);
	    verbosity_lv1_master_printf("Plaquette after %d perp %d ape smears: %16.16lg\n",
					iape+1,mu0,global_plaquette_lx_conf(sme_conf));
	    
	    //store "spatial" links and send them
	    for(int imu1=0;imu1<3;imu1++)
	      {
		int mu1=perp_dir[mu0][imu1];
		int imu01=mu0*3+imu1;
		NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
		  su3_copy(pre_transp_conf_holder[ivol],sme_conf[ivol][mu1]);
		THREAD_BARRIER();
		remap[imu01]->remap(post_transp_conf_holder,pre_transp_conf_holder,sizeof(su3));
		NISSA_PARALLEL_LOOP(icmp,0,cmp_vol[imu01])
		  su3_copy(transp_conf[imu01][icmp+cmp_vol[imu01]*(1+iape)],post_transp_conf_holder[icmp]);
		THREAD_BARRIER();
	      }
	  }
      }
    
    //free smeared conf, pre-post buffers and remappers
    nissa_free(post_transp_conf_holder);
    nissa_free(pre_transp_conf_holder);
    nissa_free(sme_conf);
    for(int imu01=0;imu01<12;imu01++) delete remap[imu01];

    ////////////////////////////////////////////////////////////////////////////////////
    
    //all the rectangles
    int dD=pars->Dmax+1-pars->Dmin;
    int dT=pars->Tmax+1-pars->Tmin;
    int nrect=dD*dT*12*nape_spat_levls;
    double *all_rectangles=nissa_malloc("all_rectangles",nrect*NACTIVE_THREADS,double);
    vector_reset(all_rectangles);
    double *all_rectangles_loc_thread=all_rectangles+nrect*THREAD_ID;
    
    //all time-lines for all distances dT, and running space-lines
    su3 *Tline=nissa_malloc("Tline",cmp_vol_max*dT,su3);
    su3 *Dline=nissa_malloc("Dline",cmp_vol_max,su3);
    
    int irect=0;
    for(int imu01=0;imu01<12;imu01++)
      {
	tricoords_t L={cmp_vol[imu01]/glb_size[mu0_l[imu01]]/glb_size[mu1_l[imu01]],glb_size[mu0_l[imu01]],
		       glb_size[mu1_l[imu01]]};
	
	//create all Tline
	NISSA_PARALLEL_LOOP(icmp,0,cmp_vol[imu01])
	  {
	    su3 U;
	    su3_copy(U,transp_conf[imu01][icmp+cmp_vol[imu01]*0]);
	    
	    for(int t=1;t<pars->Tmin;t++)
	      safe_su3_prod_su3(U,U,transp_conf[imu01][site_shift(icmp,L,1,t)+cmp_vol[imu01]*0]);
	    for(int dt=0;dt<dT;dt++)
	      {
		su3_copy(Tline[icmp*dT+dt],U);
		safe_su3_prod_su3(U,U,transp_conf[imu01][site_shift(icmp,L,1,dt+pars->Tmin)+cmp_vol[imu01]*0]);
	      }
	  }
	
	for(int iape=0;iape<nape_spat_levls;iape++)
	  {
	    //create all D'line
	    NISSA_PARALLEL_LOOP(icmp,0,cmp_vol[imu01])
	      {
		su3 U;
		su3_copy(U,transp_conf[imu01][icmp+cmp_vol[imu01]*iape]);
		
		for(int d=1;d<pars->Dmin;d++)
		  safe_su3_prod_su3(U,U,transp_conf[imu01][site_shift(icmp,L,1,d)+cmp_vol[imu01]*iape]);
		for(int dd=0;dd<dD;dd++)
		  {
		    su3_copy(Tline[icmp*dD+dd],U);
		    safe_su3_prod_su3(U,U,transp_conf[imu01][site_shift(icmp,L,1,dd+pars->Dmin)+cmp_vol[imu01]*iape]);
		  }
	      }
	    
	    //create Dlines up to Dmin
	    NISSA_PARALLEL_LOOP(icmp,0,cmp_vol[imu01])
	      {
		su3_copy(Dline[icmp],transp_conf[imu01][icmp+cmp_vol[imu01]*(1+iape)]);
		
		for(int d=1;d<pars->Dmin;d++)
		  safe_su3_prod_su3(Dline[icmp],Dline[icmp],
				    transp_conf[imu01][site_shift(icmp,L,2,d)+cmp_vol[imu01]*(1+iape)]);
	      }
	    THREAD_BARRIER();
	    
	    for(int dd=0;dd<dD;dd++)
	      {
		int d=dd+pars->Dmin;
		
		//closes
		NISSA_PARALLEL_LOOP(icmp,0,cmp_vol[imu01])
		  {
		    su3 part1,part2;
		    for(int dt=0;dt<dT;dt++)
		      {
			int t=dt+pars->Tmin;
			unsafe_su3_prod_su3(part1,Dline[icmp],Tline[site_shift(icmp,L,2,d)*dT+dt]);
			unsafe_su3_prod_su3(part2,Tline[icmp*dT+dt],Dline[site_shift(icmp,L,1,t)]);
			all_rectangles_loc_thread[irect+dt]+=real_part_of_trace_su3_prod_su3_dag(part1,part2);
		      }
		  }
		irect+=dT;
		THREAD_BARRIER();
		
		//prolong
		NISSA_PARALLEL_LOOP(icmp,0,cmp_vol[imu01])
		  safe_su3_prod_su3(Dline[icmp],Dline[icmp],
				    transp_conf[imu01][site_shift(icmp,L,2,d)+cmp_vol[imu01]*(1+iape)]);
		THREAD_BARRIER();
	      }
	  }
      }
    if(nrect!=irect) crash("expected %d rects, obtained %d",nrect,irect);
    for(int imu01=0;imu01<12;imu01++) nissa_free(transp_conf[imu01]);
    
#ifdef USE_THREADS
    if(nthreads>1)
      if(IS_MASTER_THREAD)
	for(unsigned int other_thread=1;other_thread<nthreads;other_thread++)
	  {
	    double *all_rectangles_other_thread=all_rectangles+nrect*other_thread;
	    for(int irect=0;irect<nrect;irect++) all_rectangles_loc_thread[irect]+=all_rectangles_other_thread[irect];
	  }
    THREAD_BARRIER();
#endif
    
    //perform all rank reduction and print
    double *all_rectangles_glb=nissa_malloc("all_rectangles",nrect,double);
    if(IS_MASTER_THREAD)
      {
	decript_MPI_error(MPI_Reduce(all_rectangles,all_rectangles_glb,nrect,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD),"red.");
	
	//open file
	if(rank==0)
	  {
	    FILE *fout=open_file(pars->path,create_output_file?"w":"a");
	    
	    irect=0;
	    char dir_name[5]="txyz";
	    for(int imu01=0;imu01<12;imu01++)
	      for(int iape=0;iape<nape_spat_levls;iape++)
		for(int dd=0;dd<dD;dd++)
		  for(int dt=0;dt<dT;dt++)
		    {
		      fprintf(fout,"cnf %d ism %d %c %d %c %d %16.16lg\n",
			      iconf,
			      iape,
			      dir_name[mu0_l[imu01]],dd+pars->Dmin,
			      dir_name[mu1_l[imu01]],dt+pars->Tmin,
			      all_rectangles_glb[irect++]/(3*glb_vol));
		    }
	    fclose(fout);
	  }
      }
  
    nissa_free(all_rectangles_glb);
    nissa_free(all_rectangles);
    nissa_free(Tline);
    nissa_free(Dline);
  }
  THREADABLE_FUNCTION_END
  
  //compute all possible rectangular paths among a defined interval
  THREADABLE_FUNCTION_4ARG(measure_all_rectangular_paths_old, all_rect_meas_pars_t*,pars, quad_su3*,ori_conf, int,iconf, int,create_output_file)
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
  }
  THREADABLE_FUNCTION_END
  
  void measure_all_rectangular_paths(all_rect_meas_pars_t *pars,quad_su3 **conf_eo,int iconf,int create_output_file)
  {
    quad_su3 *conf_lx=nissa_malloc("conf_lx",loc_vol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_conf(conf_lx,conf_eo);
    
    measure_all_rectangular_paths(pars,conf_lx,iconf,create_output_file);
    
    nissa_free(conf_lx);
  }
}
