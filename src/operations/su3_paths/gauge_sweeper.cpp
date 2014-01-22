#include "base/global_variables.hpp"
#include "base/debug.hpp"
#include "base/thread_macros.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#ifdef BGQ
 #include "bgq/bgq_macros.hpp"
#endif

#include <stdlib.h>

namespace nissa
{  
  //constructor
  gauge_sweeper_t::gauge_sweeper_t()
  {
    //mark not to have inited geom and staples
    comm_init_time=comp_time=comm_time=0;
    staples_inited=par_geom_inited=packing_inited=false;
    max_cached_link=max_sending_link=0;
  }
  
  //destructor
  gauge_sweeper_t::~gauge_sweeper_t()
  {
    //if we inited really, undef staples and also site counters
    if(par_geom_inited)
      {
	if(!packing_inited) nissa_free(ilink_per_staples);
	nissa_free(nsite_per_box_dir_par);
      }
    
    //deallocate staples computer
    if(staples_inited)
      {
	nissa_free(buf_out);
	nissa_free(buf_in);
	nissa_free(ivol_of_box_dir_par);
	for(int ibox=0;ibox<16;ibox++) delete box_comm[ibox];
      }
    
    //and packer
    if(packing_inited)
      {
	nissa_free(packing_link_source_dest);
	nissa_free(packing_link_buf);
      }
  }
  
  //initialize the geometry of the boxes subdirpar sets
  void gauge_sweeper_t::init_box_dir_par_geometry(int ext_gpar,int(*par_comp)(coords ivol_coord,int dir))
  {
    comm_init_time-=take_time();
    
    //check and mark as inited, store parity and allocate geometry and subbox volume
    if(par_geom_inited) crash("parity checkboard already initialized");
    par_geom_inited=true;
    gpar=ext_gpar;
    ivol_of_box_dir_par=nissa_malloc("ivol_of_box_dir_par",4*loc_vol,int);
    nsite_per_box_dir_par=nissa_malloc("nsite_per_box_dir_par",16*4*gpar,int);
    
    //find the elements of all boxes
    int ibox_dir_par=0; //order in the parity-split order
    for(int ibox=0;ibox<16;ibox++)
      for(int dir=0;dir<4;dir++)
	for(int par=0;par<gpar;par++)
	  {
	    nsite_per_box_dir_par[par+gpar*(dir+4*ibox)]=0;
	    for(int isub=0;isub<nsite_per_box[ibox];isub++)
	      {
		//get coordinates of site
		coords isub_coord;
		coord_of_lx(isub_coord,isub,box_size[ibox]);
		
		//get coords in the local size, and parity
		coords ivol_coord;
		for(int mu=0;mu<4;mu++) ivol_coord[mu]=box_size[0][mu]*box_coord[ibox][mu]+isub_coord[mu];
		int site_par=par_comp(ivol_coord,dir);
		if(site_par>=gpar||site_par<0) crash("obtained par %d while expecting in the range [0,%d]",par,gpar-1);
		
		//map sites in current parity
		if(site_par==par)
		  {
		    ivol_of_box_dir_par[ibox_dir_par++]=lx_of_coord(ivol_coord,loc_size);
		    nsite_per_box_dir_par[par+gpar*(dir+4*ibox)]++;
		  }
	      }
	  }  
    comm_init_time+=take_time();
  }
  
  //check that everything is hit once and only once
  void gauge_sweeper_t::check_hit_exactly_once()
  {
    coords *hit=nissa_malloc("hit",loc_vol,coords);
    vector_reset(hit);
    
    //mark what we hit
    int ibase=0;
    for(int ibox=0;ibox<16;ibox++)
      for(int dir=0;dir<4;dir++)
	for(int par=0;par<gpar;par++)
	  {
	    for(int ibox_dir_par=0;ibox_dir_par<nsite_per_box_dir_par[par+gpar*(dir+4*ibox)];ibox_dir_par++)
	      {
		int ivol=ivol_of_box_dir_par[ibox_dir_par+ibase];
		if(ivol>=loc_vol) crash("ibox %d ibox_par %d ibase %d ivol %d",ibox,ibox_dir_par,ibase,ivol);
		hit[ivol][dir]++;
	      }
	    ibase+=nsite_per_box_dir_par[par+gpar*(dir+4*ibox)];
	  }
    
    //check to have hit everything
    for(int ivol=0;ivol<loc_vol;ivol++)
      for(int mu=0;mu<4;mu++)
	if(hit[ivol][mu]!=1) crash("missing hit ivol %d mu %d",ivol,mu);
    
    nissa_free(hit);
  }
  
  //check that everything is hit without overlapping
  void gauge_sweeper_t::check_hit_in_the_exact_order()
  {
    int *hit=nissa_malloc("hit",4*loc_vol+max_cached_link,int);
    int ibase=0;
    for(int ibox=0;ibox<16;ibox++)
      for(int dir=0;dir<4;dir++)
        for(int par=0;par<gpar;par++)
          {
            vector_reset(hit);
            for(int iter=0;iter<2;iter++)
              //at first iter mark all site updating
              //at second iter check that none of them is internally used
              for(int ibox_dir_par=0;ibox_dir_par<nsite_per_box_dir_par[par+gpar*(dir+4*ibox)];ibox_dir_par++)
                {
                  int ivol=ivol_of_box_dir_par[ibox_dir_par+ibase];
                  
                  if(iter==0)
                    {
                      hit[dir+4*ivol]=par+1;

                      if(0)
                        printf("ibox %d dir %d par %d hitting %d, ivol %d{%d,%d,%d,%d};%d\n",
                               ibox,dir,par,
                               dir+4*ivol,ivol,
                               loc_coord_of_loclx[ivol][0],
                               loc_coord_of_loclx[ivol][1],
                               loc_coord_of_loclx[ivol][2],
                               loc_coord_of_loclx[ivol][3],
                               dir);
                    }
                  else
                    for(int ihit=0;ihit<nlinks_per_staples_of_link;ihit++)
                      {
                        int link=ilink_per_staples[ihit+nlinks_per_staples_of_link*(ibox_dir_par+ibase)];
                        if(link>=4*loc_vol+max_cached_link)
                          crash("ibox %d ibox_dir_par %d ibase %d ihit %d, link %d, max_cached_link %d",
                                ibox,ibox_dir_par,ibase,ihit,link,max_cached_link);
                        
                        if(0)
			  printf("ibox %d dir %d par %d link[%d], ivol %d{%d,%d,%d,%d}: %d,%d,%d,%d;%d\n",
				 ibox,dir,par,ihit,
				 ivol,
				 loc_coord_of_loclx[ivol][0],
				 loc_coord_of_loclx[ivol][1],
				 loc_coord_of_loclx[ivol][2],
				 loc_coord_of_loclx[ivol][3],
				 link>=4*loc_vol?-1:loc_coord_of_loclx[link/4][0],
				 link>=4*loc_vol?-1:loc_coord_of_loclx[link/4][1],
				 link>=4*loc_vol?-1:loc_coord_of_loclx[link/4][2],
				 link>=4*loc_vol?-1:loc_coord_of_loclx[link/4][3],
				 link>=4*loc_vol?-1:link%4);
                        
                        if(hit[link]!=0) crash("ivol %d:{%d,%d,%d,%d} ibox %d dir %d par %d ihit %d link %d" 
                                               "[site %d:{%d,%d,%d,%d},dir %d]: par %d",
                                               ivol,
                                               loc_coord_of_loclx[ivol][0],
                                               loc_coord_of_loclx[ivol][1],
                                               loc_coord_of_loclx[ivol][2],
                                               loc_coord_of_loclx[ivol][3],
                                               ibox,dir,par,ihit,link,link/4,
                                               link>=4*loc_vol?-1:loc_coord_of_loclx[link/4][0],
                                               link>=4*loc_vol?-1:loc_coord_of_loclx[link/4][1],
                                               link>=4*loc_vol?-1:loc_coord_of_loclx[link/4][2],
                                               link>=4*loc_vol?-1:loc_coord_of_loclx[link/4][3],
                                               link%4,hit[link]-1);
                      }
                }
            ibase+=nsite_per_box_dir_par[par+gpar*(dir+4*ibox)];
          }
    
    nissa_free(hit);
  }
  
  //init the list of staples
  void gauge_sweeper_t::init_staples(int ext_nlinks_per_staples_of_link,void(*ext_add_staples_per_link)(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu),void (*ext_compute_staples)(su3 staples,su3 *links,int *ilinks))
  {
    //take external nlinks and mark
    if(!par_geom_inited) crash("call geom initializer before");
    if(staples_inited) crash("staples already initialized");
    staples_inited=true;
    nlinks_per_staples_of_link=ext_nlinks_per_staples_of_link;
    add_staples_per_link=ext_add_staples_per_link;
    compute_staples=ext_compute_staples;
    
    //allocate
    ilink_per_staples=nissa_malloc("ilink_per_staples",nlinks_per_staples_of_link*4*loc_vol,int);
    all_to_all_gathering_list_t *gl[16];
    for(int ibox=0;ibox<16;ibox++) gl[ibox]=new all_to_all_gathering_list_t;
    add_staples_required_links(gl);
    
    //initialize the communicator
    for(int ibox=0;ibox<16;ibox++)
      {
	box_comm[ibox]=new all_to_all_comm_t(*(gl[ibox]));
	delete gl[ibox];
      }
    
    //compute the maximum number of link to send and receive and allocate buffers
    for(int ibox=0;ibox<16;ibox++)
      {
	max_cached_link=std::max(max_cached_link,box_comm[ibox]->nel_in);
	max_sending_link=std::max(max_sending_link,box_comm[ibox]->nel_out);
      }
    buf_out=nissa_malloc("buf_out",max_sending_link,su3);
    buf_in=nissa_malloc("buf_in",max_cached_link,su3);
    
    //check cached
    verbosity_lv3_master_printf("Max cached links: %d\n",max_cached_link);
    if(max_cached_link>bord_vol+edge_vol) crash("larger buffer needed [really? recheck this]");

    //perform two checks
    check_hit_exactly_once();
    check_hit_in_the_exact_order();
  }
  
  //ordering for link_source_dest
  int compare_link_source_dest(const void *a,const void *b)
  {return ((int*)a)[0]-((int*)b)[0];}
  
  //find the place where each link must be copied to access it sequentially
  void gauge_sweeper_t::find_packing_index(void (*ext_compute_staples_packed)(su3 staples,su3 *links))
  {
    if(!packing_inited)
      {	
	//mark packing to be have been inited and allocate
	packing_inited=true;
	compute_staples_packed=ext_compute_staples_packed;
	packing_link_source_dest=nissa_malloc("packing_link_source_dest",2*nlinks_per_staples_of_link*4*loc_vol,int);
	
	int ibase=0;
	int max_packing_link_nel=0;
	for(int ibox=0;ibox<16;ibox++)
	  for(int dir=0;dir<4;dir++)
	    for(int par=0;par<gpar;par++)
	      {
		//find the packing size
		int ns=nsite_per_box_dir_par[par+gpar*(dir+4*ibox)];
		max_packing_link_nel=std::max(max_packing_link_nel,nlinks_per_staples_of_link*ns);
		
		//scan the destination
		for(int ibox_dir_par=0;ibox_dir_par<ns;ibox_dir_par++)
		  for(int ilink=0;ilink<nlinks_per_staples_of_link;ilink++)
		    {		  
		      packing_link_source_dest[2*(ilink+nlinks_per_staples_of_link*(ibox_dir_par+ibase))+0]=
			ilink_per_staples[(ibox_dir_par+ibase)*nlinks_per_staples_of_link+ilink];
		      packing_link_source_dest[2*(ilink+nlinks_per_staples_of_link*(ibox_dir_par+ibase))+1]=
#ifdef BGQ
			//if using BGQ, pack info on vnode too, as last bit
			((ilink+nlinks_per_staples_of_link*(ibox_dir_par%(ns/2)))<<1)+
			 (2*ibox_dir_par>=ns) //vnode bit
#else
			ilink+nlinks_per_staples_of_link*ibox_dir_par
#endif
			;
		    }
		
	       //sort and increase the base
	       qsort(packing_link_source_dest+2*nlinks_per_staples_of_link*ibase,
		     nsite_per_box_dir_par[par+gpar*(dir+4*ibox)],2*sizeof(int),compare_link_source_dest);
	       ibase+=nsite_per_box_dir_par[par+gpar*(dir+4*ibox)];	       
	      }
	
	//allocate packing link and deallocate ilink_per_staples
	packing_link_buf=nissa_malloc("packing_link_buf",max_packing_link_nel,su3);
	nissa_free(ilink_per_staples);
      }
  }
  
  //add all the links needed to compute staples separately for each box
  THREADABLE_FUNCTION_2ARG(add_staples_required_links_to_gauge_sweep, gauge_sweeper_t*,gs, all_to_all_gathering_list_t**,gl)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ibox,0,16)
      {
	if(rank==0) printf("box %d being initialized by rank %d\n",ibox,thread_id);
	//find base for curr box
	int ibase=0;
	for(int jbox=0;jbox<ibox;jbox++) ibase+=nsite_per_box[jbox];
	ibase*=4;
	
	//scan all the elements of sub-box, selecting only those with the good parity
	for(int dir=0;dir<4;dir++)
	  for(int par=0;par<gs->gpar;par++)
	    {	  
	      for(int ibox_dir_par=ibase;ibox_dir_par<ibase+gs->nsite_per_box_dir_par[par+gs->gpar*(dir+4*ibox)];
		  ibox_dir_par++)
		{
		  int ivol=gs->ivol_of_box_dir_par[ibox_dir_par];
		  gs->add_staples_per_link(gs->ilink_per_staples+ibox_dir_par*gs->nlinks_per_staples_of_link,*(gl[ibox]),
					   ivol,dir);
		}
	      ibase+=gs->nsite_per_box_dir_par[par+gs->gpar*(dir+4*ibox)];
	    }
    }
    THREAD_BARRIER();
  }
  THREADABLE_FUNCTION_END

  //wrapper to use threads
  void gauge_sweeper_t::add_staples_required_links(all_to_all_gathering_list_t **gl)
  {
    add_staples_required_links_to_gauge_sweep(this,gl);
  }

  //update a single link using either heat-bath or overrelaxation
  void gauge_sweeper_t::update_link_using_staples(quad_su3 *conf,int ivol,int dir,su3 staples,quenched_update_alg_t update_alg,double beta,int nhits)
  {
    switch(update_alg)
      {
      case HEATBATH: su3_find_heatbath(conf[ivol][dir],conf[ivol][dir],staples,beta,nhits,loc_rnd_gen+ivol);break;
      case OVERRELAX: su3_find_overrelaxed(conf[ivol][dir],conf[ivol][dir],staples,nhits);break;
      case COOL_FULLY: su3_unitarize_maximal_trace_projecting(conf[ivol][dir],staples);break;
      case COOL_PARTLY: su3_unitarize_maximal_trace_projecting_iteration(conf[ivol][dir],staples);break;
      default: crash("unknown update algorithm");
      }
  }
  
#ifdef BGQ
  void compute_tlSym_staples_packed_bgq(su3 staples1,su3 staples2,bi_su3 *links);
#endif
  
  //sweep the conf using the appropriate routine to compute staples
  THREADABLE_FUNCTION_5ARG(sweep_conf_fun, quad_su3*,conf, gauge_sweeper_t*,gs, quenched_update_alg_t,update_alg, double,beta, int,nhits)
  {
    GET_THREAD_ID();
  
    int ibase=0;
    for(int ibox=0;ibox<16;ibox++)
      {
	//communicate needed links
	if(IS_MASTER_THREAD) gs->comm_time-=take_time();
	gs->box_comm[ibox]->communicate(conf,conf,sizeof(su3),NULL,NULL,ibox+100);
	master_printf("Barriering\n");
	THREAD_BARRIER();
	MPI_Barrier(MPI_COMM_WORLD);
	if(IS_MASTER_THREAD) gs->comm_time+=take_time();
      
	if(IS_MASTER_THREAD) gs->comp_time-=take_time();
	for(int dir=0;dir<4;dir++)
	  for(int par=0;par<gs->gpar;par++)
	    {
	      int nbox_dir_par=gs->nsite_per_box_dir_par[par+gs->gpar*(dir+4*ibox)];
	      
	      //pack
	      if(gs->packing_inited)
		{
		  int *packing_link_source_dest=gs->packing_link_source_dest+2*gs->nlinks_per_staples_of_link*ibase;
		  NISSA_PARALLEL_LOOP(ilink_to_ship,0,gs->nlinks_per_staples_of_link*nbox_dir_par)
		    {
		      int isource=packing_link_source_dest[2*ilink_to_ship+0];
		      int idest=packing_link_source_dest[2*ilink_to_ship+1];
#ifdef BGQ
		      int true_dest=(idest>>1);
		      int vnode=idest&1;
		      //master_printf("source %d (%lg) dest: %d/%d %d\n",isource,((su3*)conf)[isource][0][0][0],
		      //true_dest,get_vect(gs->packing_link_buf)->nel/2,vnode);
		      SU3_TO_BI_SU3(((bi_su3*)gs->packing_link_buf)[true_dest],((su3*)conf)[isource],vnode);
#else
		      su3_copy(gs->packing_link_buf[idest],((su3*)conf)[isource]);
		      //master_printf("source %d (%lg) dest: %d/%d\n",
		      //isource,((su3*)conf)[isource][0][0][0],idest,get_vect(gs->packing_link_buf)->nel);
#endif	      
		    }
		  THREAD_BARRIER();
		}
	      
#ifdef BGQ	      
	      su3 *staples_list=nissa_malloc("staples_list",nbox_dir_par,su3);
	      if(gs->packing_inited)
		NISSA_PARALLEL_LOOP(ibox_dir_par,0,nbox_dir_par/2)
		  {
		   //master_printf("ibox_dir_par: %d\n");
		   compute_tlSym_staples_packed_bgq(staples_list[ibox_dir_par],staples_list[ibox_dir_par+nbox_dir_par/2],
						    ((bi_su3*)gs->packing_link_buf)+ibox_dir_par*gs->nlinks_per_staples_of_link);
		  }
	      THREAD_BARRIER();
#endif	      
	      
	      //scan all the box
	      NISSA_PARALLEL_LOOP(ibox_dir_par,ibase,ibase+nbox_dir_par)
		{
		  //compute the staples
		  su3 staples;
		  
		  if(gs->packing_inited) 
		    {
#ifdef BGQ
		      //su3 staples_temp;
		      su3_copy(staples,staples_list[ibox_dir_par-ibase]);
#else
		      gs->compute_staples_packed(staples,
					      gs->packing_link_buf+(ibox_dir_par-ibase)*gs->nlinks_per_staples_of_link);
#endif
		    }
		  else
		    gs->compute_staples(staples,
					   (su3*)conf,gs->ilink_per_staples+gs->nlinks_per_staples_of_link*ibox_dir_par);

		  //su3 staples_temp;
		  //gs->compute_staples(staples_temp,
		  //(su3*)conf,gs->ilink_per_staples+gs->nlinks_per_staples_of_link*ibox_dir_par);
		  
		  //master_printf("ibox: %d\n",ibox_dir_par-ibase);
#ifdef BGQ
		  if(gs->packing_inited)
		    {
		      //su3 temp;
		      //su3_subt(temp,staples,staples_temp);
		      //master_printf("%lg\n",real_part_of_trace_su3_prod_su3_dag(temp,temp));
		      //master_printf("bgq:\n");
		      //su3_print(staples_temp);
		    }
#endif	  
		  //master_printf("pc:\n");
		  //su3_print(staples);
		  
		  //find new link
		  int ivol=gs->ivol_of_box_dir_par[ibox_dir_par];
		  gs->update_link_using_staples(conf,ivol,dir,staples,update_alg,beta,nhits);
		}
	      
#ifdef BGQ
	      nissa_free(staples_list);
#endif
	      //increment the box-dir-par subset
	      ibase+=nbox_dir_par;
	      THREAD_BARRIER();
	    }
	if(IS_MASTER_THREAD) gs->comp_time+=take_time();
      }
  
    set_borders_invalid(conf);
  }
  THREADABLE_FUNCTION_END

  //wrapper to use threads
  void gauge_sweeper_t::sweep_conf(quad_su3 *conf,quenched_update_alg_t up_alg,double beta,int nhits)
  {sweep_conf_fun(conf,this,up_alg,beta,nhits);}
  
  //return the appropriate sweeper
  gauge_sweeper_t *get_sweeper(gauge_action_name_t gauge_action_name)
  {
    gauge_sweeper_t *sweeper=NULL;
    switch(gauge_action_name)
      {
      case WILSON_GAUGE_ACTION:sweeper=Wilson_sweeper;break;
      case TLSYM_GAUGE_ACTION:sweeper=tlSym_sweeper;break;
      case UNSPEC_GAUGE_ACTION:crash("unspecified action");break;
      default: crash("not implemented action");break;
      }
    
    return sweeper;
  }

  ///////////////////////////////////// tlSym action //////////////////////////
  
  //compute the parity according to the tlSym requirements
  int tlSym_par(coords ivol_coord,int dir)
  {
    int site_par=0;
    for(int mu=0;mu<4;mu++) site_par+=((mu==dir)?2:1)*ivol_coord[mu];

    site_par=site_par%4;
  
    return site_par;
  }

  //add all links needed for a certain site
  void add_tlSym_staples(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu)
  {
    int *c=glb_coord_of_loclx[ivol];      
    coords A={c[0],c[1],c[2],c[3]};                         //       P---O---N    
    coords B,C,D,E,F,G,H,I,J,K,L,M,N,O,P;                   //       |   |   |    
    //find coord mu                                         //   H---G---F---E---D
    K[mu]=L[mu]=M[mu]=(A[mu]-1+glb_size[mu])%glb_size[mu];  //   |   |   |   |   |
    I[mu]=J[mu]=B[mu]=C[mu]=A[mu];                          //   I---J---A---B---C
    D[mu]=E[mu]=F[mu]=G[mu]=H[mu]=(A[mu]+1)%glb_size[mu];   //       |   |   |    
    N[mu]=O[mu]=P[mu]=(A[mu]+2)%glb_size[mu];               //       K---L---M    
    for(int inu=0;inu<3;inu++)
      {            
	int &nu=perp_dir[mu][inu];
	int &rh=perp2_dir[mu][inu][0];
	int &si=perp2_dir[mu][inu][1];
      
	//find coord rho
	B[rh]=C[rh]=D[rh]=E[rh]=F[rh]=G[rh]=H[rh]=I[rh]=J[rh]=K[rh]=L[rh]=M[rh]=N[rh]=O[rh]=P[rh]=A[rh];
	//find coord sigma
	B[si]=C[si]=D[si]=E[si]=F[si]=G[si]=H[si]=I[si]=J[si]=K[si]=L[si]=M[si]=N[si]=O[si]=P[si]=A[si];
	//find coord nu
	H[nu]=I[nu]=(A[nu]-2+glb_size[nu])%glb_size[nu];
	K[nu]=J[nu]=G[nu]=P[nu]=(I[nu]+1)%glb_size[nu];
	L[nu]=F[nu]=O[nu]=A[nu];
	M[nu]=B[nu]=E[nu]=N[nu]=(A[nu]+1)%glb_size[nu];
	C[nu]=D[nu]=(A[nu]+2)%glb_size[nu];
      
	//backward square staple
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(J,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(J,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(G,nu);
	//forward square staple
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(A,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(B,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(F,nu);
	//backward dw rectangle
	//*(ilink_to_be_used++)=gat.add_conf_link_for_paths(L,mu); //vertical common link dw (see below)
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(K,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(K,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(J,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(G,nu);
	//backward backward rectangle
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(J,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(I,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(I,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(H,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(G,nu);
	//backward up rectangle
	//*(ilink_to_be_used++)=gat.add_conf_link_for_paths(J,nu); //already computed at...
	//*(ilink_to_be_used++)=gat.add_conf_link_for_paths(J,mu); //...backward square staple
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(G,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(P,nu);
	//*(ilink_to_be_used++)=gat.add_conf_link_for_paths(F,mu); //vertical common link up (see below)
	//forward dw rectangle
	//*(ilink_to_be_used++)=gat.add_conf_link_for_paths(L,mu); //vertical common link dw (see below)
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(L,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(M,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(B,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(F,nu);
	//forward forward rectangle
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(A,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(B,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(C,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(E,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(F,nu);
	//forward up rectangle
	//*(ilink_to_be_used++)=gat.add_conf_link_for_paths(A,nu); //already computed at...
	//*(ilink_to_be_used++)=gat.add_conf_link_for_paths(B,mu); //...forward square staple
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(E,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(O,nu);
	//*(ilink_to_be_used++)=gat.add_conf_link_for_paths(F,mu); //vertical common link up (see below)
      }
    *(ilink_to_be_used++)=gat.add_conf_link_for_paths(F,mu); //vertical common link up
    *(ilink_to_be_used++)=gat.add_conf_link_for_paths(L,mu); //verical common link dw
  
    //8*3 link missing, 2 readded = -22 links
  }

  //compute the summ of the staples pointed by "ilinks"
  void compute_tlSym_staples(su3 staples,su3 *links,int *ilinks)
  {
    su3 squares,rectangles,up_rectangles,dw_rectangles;
    su3_put_to_zero(squares);
    su3_put_to_zero(rectangles);
    su3_put_to_zero(up_rectangles);
    su3_put_to_zero(dw_rectangles);
    
    const int PARTIAL=2;
    for(int inu=0;inu<3;inu++)
      {  
	su3 hb;
	//backward square staple
	unsafe_su3_dag_prod_su3(hb,links[ilinks[ 0]],links[ilinks[ 1]]);
	su3_summ_the_prod_su3(squares,hb,links[ilinks[ 2]]);
	//forward square staple
	su3 hf;
	unsafe_su3_prod_su3(hf,links[ilinks[ 3]],links[ilinks[ 4]]);
	su3_summ_the_prod_su3_dag(squares,hf,links[ilinks[ 5]]);
      
	su3 temp1,temp2;
	//backward dw rectangle
	unsafe_su3_dag_prod_su3(temp2,links[ilinks[ 6]],links[ilinks[ 7]],PARTIAL);
	unsafe_su3_prod_su3(temp1,temp2,links[ilinks[ 8]],PARTIAL);
	su3_build_third_row(temp1);
	su3_summ_the_prod_su3(dw_rectangles,temp1,links[ilinks[9]]);
	//backward backward rectangle
	unsafe_su3_dag_prod_su3_dag(temp1,links[ilinks[10]],links[ilinks[11]],PARTIAL);
	unsafe_su3_prod_su3(temp2,temp1,links[ilinks[12]],PARTIAL);
	unsafe_su3_prod_su3(temp1,temp2,links[ilinks[13]],PARTIAL);
	su3_build_third_row(temp1);
	su3_summ_the_prod_su3(rectangles,temp1,links[ilinks[14]]);
	//backward up rectangle
	unsafe_su3_prod_su3(temp2,hb,links[ilinks[15]]);
	su3_summ_the_prod_su3(up_rectangles,temp2,links[ilinks[16]]);
	//forward dw rectangle
	unsafe_su3_prod_su3(temp2,links[ilinks[17]],links[ilinks[18]],PARTIAL);
	unsafe_su3_prod_su3(temp1,temp2,links[ilinks[19]],PARTIAL);
	su3_build_third_row(temp1);
	su3_summ_the_prod_su3_dag(dw_rectangles,temp1,links[ilinks[20]]);
	//forward forward rectangle
	unsafe_su3_prod_su3(temp1,links[ilinks[21]],links[ilinks[22]],PARTIAL);
	unsafe_su3_prod_su3(temp2,temp1,links[ilinks[23]],PARTIAL);
	unsafe_su3_prod_su3_dag(temp1,temp2,links[ilinks[24]],PARTIAL);
	su3_build_third_row(temp1);
	su3_summ_the_prod_su3_dag(rectangles,temp1,links[ilinks[25]]);
	//forward up rectangle
	unsafe_su3_prod_su3(temp2,hf,links[ilinks[26]]);
	su3_summ_the_prod_su3_dag(up_rectangles,temp2,links[ilinks[27]]);

	ilinks+=28;
      }
    
    //close the two partial rectangles
    su3_summ_the_prod_su3_dag(rectangles,up_rectangles,links[ilinks[ 0]]);
    su3_summ_the_dag_prod_su3(rectangles,links[ilinks[ 1]],dw_rectangles);
  
    //compute the summed staples
    double b1=-1.0/12,b0=1-8*b1;
    su3_linear_comb(staples,squares,b0,rectangles,b1);
  }
  
  //compute the summ of the staples pointed by "ilinks"
  void compute_tlSym_staples_packed(su3 staples,su3 *links)
  {
    su3 squares,rectangles,up_rectangles,dw_rectangles;
    su3_put_to_zero(squares);
    su3_put_to_zero(rectangles);
    su3_put_to_zero(up_rectangles);
    su3_put_to_zero(dw_rectangles);
    
    const int PARTIAL=2;
    for(int inu=0;inu<3;inu++)
      {  
	su3 hb;
	//backward square staple
	unsafe_su3_dag_prod_su3(hb,links[ 0],links[ 1]);
	su3_summ_the_prod_su3(squares,hb,links[ 2]);
	//forward square staple
	su3 hf;
	unsafe_su3_prod_su3(hf,links[ 3],links[ 4]);
	su3_summ_the_prod_su3_dag(squares,hf,links[ 5]);
      
	su3 temp1,temp2;
	//backward dw rectangle
	unsafe_su3_dag_prod_su3(temp2,links[ 6],links[ 7],PARTIAL);
	unsafe_su3_prod_su3(temp1,temp2,links[ 8],PARTIAL);
	su3_build_third_row(temp1);
	su3_summ_the_prod_su3(dw_rectangles,temp1,links[9]);
	//backward backward rectangle
	unsafe_su3_dag_prod_su3_dag(temp1,links[10],links[11],PARTIAL);
	unsafe_su3_prod_su3(temp2,temp1,links[12],PARTIAL);
	unsafe_su3_prod_su3(temp1,temp2,links[13],PARTIAL);
	su3_build_third_row(temp1);
	su3_summ_the_prod_su3(rectangles,temp1,links[14]);
	//backward up rectangle
	unsafe_su3_prod_su3(temp2,hb,links[15]);
	su3_summ_the_prod_su3(up_rectangles,temp2,links[16]);
	//forward dw rectangle
	unsafe_su3_prod_su3(temp2,links[17],links[18],PARTIAL);
	unsafe_su3_prod_su3(temp1,temp2,links[19],PARTIAL);
	su3_build_third_row(temp1);
	su3_summ_the_prod_su3_dag(dw_rectangles,temp1,links[20]);
	//forward forward rectangle
	unsafe_su3_prod_su3(temp1,links[21],links[22],PARTIAL);
	unsafe_su3_prod_su3(temp2,temp1,links[23],PARTIAL);
	unsafe_su3_prod_su3_dag(temp1,temp2,links[24],PARTIAL);
	su3_build_third_row(temp1);
	su3_summ_the_prod_su3_dag(rectangles,temp1,links[25]);
	//forward up rectangle
	unsafe_su3_prod_su3(temp2,hf,links[26]);
	su3_summ_the_prod_su3_dag(up_rectangles,temp2,links[27]);

	links+=28;
      }
    
    //close the two partial rectangles
    su3_summ_the_prod_su3_dag(rectangles,up_rectangles,links[ 0]);
    su3_summ_the_dag_prod_su3(rectangles,links[ 1],dw_rectangles);
  
    //compute the summed staples
    double b1=-1.0/12,b0=1-8*b1;
    su3_linear_comb(staples,squares,b0,rectangles,b1);
  }
  
#ifdef BGQ
  //compute the summ of the staples pointed by "ilinks"
  void compute_tlSym_staples_packed_bgq(su3 staples1,su3 staples2,bi_su3 *links)
  {
    bi_su3 squares,rectangles,up_rectangles,dw_rectangles;
    
    BI_SU3_PUT_TO_ZERO(squares);
    BI_SU3_PUT_TO_ZERO(rectangles);
    BI_SU3_PUT_TO_ZERO(up_rectangles);
    BI_SU3_PUT_TO_ZERO(dw_rectangles);
  
    DECLARE_REG_BI_SU3(REG_U1);
    DECLARE_REG_BI_SU3(REG_U2);
    DECLARE_REG_BI_SU3(REG_U3);
    
    for(int inu=0;inu<3;inu++)
      {  
	//backward square staple
	bi_su3 hb;
	REG_LOAD_BI_SU3(REG_U2,links[0]);
	REG_LOAD_BI_SU3(REG_U3,links[1]);
	REG_BI_SU3_DAG_PROD_BI_SU3(REG_U1,REG_U2,REG_U3);
	STORE_REG_BI_SU3(hb,REG_U1);
	REG_LOAD_BI_SU3(REG_U2,squares);
	REG_LOAD_BI_SU3(REG_U3,links[2]);
	REG_BI_SU3_SUMM_THE_PROD_BI_SU3(REG_U2,REG_U1,REG_U3);
	STORE_REG_BI_SU3(squares,REG_U2);
	//forward square staple
	bi_su3 hf;
	REG_LOAD_BI_SU3(REG_U2,links[3]);
	REG_LOAD_BI_SU3(REG_U3,links[4]);
	REG_BI_SU3_PROD_BI_SU3(REG_U1,REG_U2,REG_U3);
	STORE_REG_BI_SU3(hf,REG_U1);
	REG_LOAD_BI_SU3(REG_U2,squares);
	REG_LOAD_BI_SU3(REG_U3,links[5]);
	REG_BI_SU3_SUMM_THE_PROD_BI_SU3_DAG(REG_U2,REG_U1,REG_U3);
	STORE_REG_BI_SU3(squares,REG_U2);
      
	//backward dw rectangle
	REG_LOAD_BI_SU3(REG_U2,links[6]);
	REG_LOAD_BI_SU3(REG_U3,links[7]);
	REG_BI_SU3_DAG_PROD_BI_SU3(REG_U1,REG_U2,REG_U3);
	REG_LOAD_BI_SU3(REG_U3,links[8]);
	REG_BI_SU3_PROD_BI_SU3(REG_U2,REG_U1,REG_U3);
	REG_LOAD_BI_SU3(REG_U1,dw_rectangles);
	REG_LOAD_BI_SU3(REG_U3,links[9]);
	REG_BI_SU3_SUMM_THE_PROD_BI_SU3(REG_U1,REG_U2,REG_U3);
	STORE_REG_BI_SU3(dw_rectangles,REG_U1);
	//backward backward rectangle
	REG_LOAD_BI_SU3(REG_U2,links[10]);
	REG_LOAD_BI_SU3(REG_U3,links[11]);
	REG_BI_SU3_PROD_BI_SU3(REG_U1,REG_U3,REG_U2); //not dag dag, so reverted
	REG_LOAD_BI_SU3(REG_U3,links[12]);
	REG_BI_SU3_DAG_PROD_BI_SU3(REG_U2,REG_U1,REG_U3); //daggering the first instead
	REG_LOAD_BI_SU3(REG_U3,links[13]);
	REG_BI_SU3_PROD_BI_SU3(REG_U1,REG_U2,REG_U3);
	REG_LOAD_BI_SU3(REG_U3,links[14]);
	REG_LOAD_BI_SU3(REG_U2,rectangles);
	REG_BI_SU3_SUMM_THE_PROD_BI_SU3(REG_U2,REG_U1,REG_U3);
	STORE_REG_BI_SU3(rectangles,REG_U2);
	//backward up rectangle
	REG_LOAD_BI_SU3(REG_U2,hb);
	REG_LOAD_BI_SU3(REG_U3,links[15]);
        REG_BI_SU3_PROD_BI_SU3(REG_U1,REG_U2,REG_U3);
	REG_LOAD_BI_SU3(REG_U3,links[16]);
	REG_LOAD_BI_SU3(REG_U2,up_rectangles);
	REG_BI_SU3_SUMM_THE_PROD_BI_SU3(REG_U2,REG_U1,REG_U3);
	STORE_REG_BI_SU3(up_rectangles,REG_U2);
	//forward dw rectangle
	REG_LOAD_BI_SU3(REG_U2,links[17]);
	REG_LOAD_BI_SU3(REG_U3,links[18]);
	REG_BI_SU3_PROD_BI_SU3(REG_U1,REG_U2,REG_U3);
	REG_LOAD_BI_SU3(REG_U3,links[19]);
	REG_BI_SU3_PROD_BI_SU3(REG_U2,REG_U1,REG_U3);
	REG_LOAD_BI_SU3(REG_U1,dw_rectangles);
	REG_LOAD_BI_SU3(REG_U3,links[20]);
	REG_BI_SU3_SUMM_THE_PROD_BI_SU3_DAG(REG_U1,REG_U2,REG_U3);
	STORE_REG_BI_SU3(dw_rectangles,REG_U1);
	//forward forward rectangle
	REG_LOAD_BI_SU3(REG_U2,links[21]);
	REG_LOAD_BI_SU3(REG_U3,links[22]);
	REG_BI_SU3_PROD_BI_SU3(REG_U1,REG_U2,REG_U3);
	REG_LOAD_BI_SU3(REG_U3,links[23]);
	REG_BI_SU3_PROD_BI_SU3(REG_U2,REG_U1,REG_U3);
	REG_LOAD_BI_SU3(REG_U3,links[24]);
	REG_BI_SU3_PROD_BI_SU3_DAG(REG_U1,REG_U2,REG_U3);
	REG_LOAD_BI_SU3(REG_U3,links[25]);
	REG_LOAD_BI_SU3(REG_U2,rectangles);
	REG_BI_SU3_SUMM_THE_PROD_BI_SU3_DAG(REG_U2,REG_U1,REG_U3);
	STORE_REG_BI_SU3(rectangles,REG_U2);
	//forward up rectangle
	REG_LOAD_BI_SU3(REG_U2,hf);
	REG_LOAD_BI_SU3(REG_U3,links[26]);
        REG_BI_SU3_PROD_BI_SU3(REG_U1,REG_U2,REG_U3);
	REG_LOAD_BI_SU3(REG_U3,links[27]);
	REG_LOAD_BI_SU3(REG_U2,up_rectangles);
	REG_BI_SU3_SUMM_THE_PROD_BI_SU3_DAG(REG_U2,REG_U1,REG_U3);
	STORE_REG_BI_SU3(up_rectangles,REG_U2);

	links+=28;
      }
    
    //close the two partial rectangles
    REG_LOAD_BI_SU3(REG_U1,rectangles);
    REG_LOAD_BI_SU3(REG_U2,up_rectangles);
    REG_LOAD_BI_SU3(REG_U3,links[0]);
    REG_BI_SU3_SUMM_THE_PROD_BI_SU3_DAG(REG_U1,REG_U2,REG_U3);
    REG_LOAD_BI_SU3(REG_U2,dw_rectangles);
    REG_LOAD_BI_SU3(REG_U3,links[1]);
    REG_BI_SU3_DAG_SUMM_THE_PROD_BI_SU3(REG_U1,REG_U3,REG_U2);
    STORE_REG_BI_SU3(rectangles,REG_U1);
    
    //compute the summed staples
    double b1=-1.0/12,b0=1-8*b1;
    DECLARE_REG_BI_COMPLEX(reg_b0);
    DECLARE_REG_BI_COMPLEX(reg_b1);
    REG_SPLAT_BI_COMPLEX(reg_b0,b0);
    REG_SPLAT_BI_COMPLEX(reg_b1,b1);
    REG_LOAD_BI_SU3(REG_U2,squares);
    REG_BI_SU3_PROD_4DOUBLE(REG_U3,REG_U2,reg_b0);
    REG_BI_SU3_SUMM_THE_PROD_4DOUBLE(REG_U3,REG_U3,REG_U1,reg_b1);
    
    //split staples
    bi_su3 bi_staples;
    STORE_REG_BI_SU3(bi_staples,REG_U3);
    BI_SU3_TO_SU3(staples1,staples2,bi_staples);
  }
#endif

  //initialize the tlSym sweeper using the above defined routines
  void init_tlSym_sweeper()
  {
    if(!tlSym_sweeper->staples_inited)
      {
	verbosity_lv3_master_printf("Initializing tlSym sweeper\n");
	//checking consistency for gauge_sweeper initialization
	for(int mu=0;mu<4;mu++) if(loc_size[mu]<4) crash("loc_size[%d]=%d must be at least 4",mu,loc_size[mu]);
	//initialize the tlSym sweeper
	const int nlinks_per_tlSym_staples_of_link=3*2*(3+5*3)-3*8+2;
	tlSym_sweeper->init_box_dir_par_geometry(4,tlSym_par);
	tlSym_sweeper->init_staples(nlinks_per_tlSym_staples_of_link,add_tlSym_staples,compute_tlSym_staples);
#ifdef BGQ
	tlSym_sweeper->find_packing_index(compute_tlSym_staples_packed);
#endif
      }
  }

  ///////////////////////////////////////// Wilson ////////////////////////////////////////
  
  //compute the parity according to the Wilson requirements
  int Wilson_par(coords ivol_coord,int dir)
  {
    int site_par=0;
    for(int mu=0;mu<4;mu++) site_par+=ivol_coord[mu];

    site_par=site_par%2;
  
    return site_par;
  }

  //add all links needed for a certain site
  void add_Wilson_staples(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu)
  {
    int *c=glb_coord_of_loclx[ivol];      
    coords A={c[0],c[1],c[2],c[3]};
    coords B,F,G,J;                                         //       G---F---E
    //find coord mu                                         //       |   |   |
    J[mu]=B[mu]=A[mu];                                      //       J---A---B
    F[mu]=G[mu]=(A[mu]+1)%glb_size[mu];
    for(int inu=0;inu<3;inu++)
      {            
	int &nu=perp_dir[mu][inu];
	int &rh=perp2_dir[mu][inu][0];
	int &si=perp2_dir[mu][inu][1];
	
	//find coord rho
	B[rh]=F[rh]=G[rh]=J[rh]=A[rh];
	//find coord sigma
	B[si]=F[si]=G[si]=J[si]=A[si];
	//find coord nu
	J[nu]=G[nu]=(A[nu]-1+glb_size[nu])%glb_size[nu];
	F[nu]=A[nu];
	B[nu]=(A[nu]+1)%glb_size[nu];
      
	//backward square staple
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(J,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(J,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(G,nu);
	//forward square staple
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(A,nu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(B,mu);
	*(ilink_to_be_used++)=gat.add_conf_link_for_paths(F,nu);
      }
  }

  //compute the summ of the staples pointed by "ilinks"
  void compute_Wilson_staples(su3 staples,su3 *links,int *ilinks)
  {
    su3_put_to_zero(staples);
  
    for(int inu=0;inu<3;inu++)
      {  
	su3 hb;
	//backward square staple
	unsafe_su3_dag_prod_su3(hb,links[ilinks[ 0]],links[ilinks[ 1]]);
	su3_summ_the_prod_su3(staples,hb,links[ilinks[ 2]]);
	//forward square staple
	su3 hf;
	unsafe_su3_prod_su3(hf,links[ilinks[ 3]],links[ilinks[ 4]]);
	su3_summ_the_prod_su3_dag(staples,hf,links[ilinks[ 5]]);

	ilinks+=6;
      }
  }
  
  //initialize the Wilson sweeper using the above defined routines
  void init_Wilson_sweeper()
  {
    if(!Wilson_sweeper->staples_inited)
      {
	verbosity_lv3_master_printf("Initializing Wilson sweeper\n");
	//checking consistency for gauge_sweeper initialization
	for(int mu=0;mu<4;mu++) if(loc_size[mu]<2) crash("loc_size[%d]=%d must be at least 2",mu,loc_size[mu]);
	//initialize the Wilson sweeper
	Wilson_sweeper->init_box_dir_par_geometry(2,Wilson_par);
	const int nlinks_per_Wilson_staples_of_link=18;
	Wilson_sweeper->init_staples(nlinks_per_Wilson_staples_of_link,add_Wilson_staples,compute_Wilson_staples);
      }
  }
  
  //call the appropriate sweeper intializator
  void init_sweeper(gauge_action_name_t gauge_action_name)
  {
#ifdef USE_THREADS
    if(!thread_pool_locked) crash("call from non-parallel environment");
#endif
    switch(gauge_action_name)
      {
      case WILSON_GAUGE_ACTION:if(!Wilson_sweeper->staples_inited) init_Wilson_sweeper();break;
      case TLSYM_GAUGE_ACTION:if(!tlSym_sweeper->staples_inited) init_tlSym_sweeper();break;
      case UNSPEC_GAUGE_ACTION:crash("unspecified action");break;
      default: crash("not implemented action");break;
      }
  }
}
