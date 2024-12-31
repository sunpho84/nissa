#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/debug.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/gauge/Symanzik_action.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"

#include <stdlib.h>

#define EXTERN
#include "gauge_sweeper.hpp"

namespace nissa
{
  void add_Symanzik_staples(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu);
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
	for(int ibox=0;ibox<(1<<NDIM);ibox++) delete box_comm[ibox];
      }
    
    //and packer
    if(packing_inited)
      {
	nissa_free(packing_link_source_dest);
	nissa_free(packing_link_buf);
      }
  }
  
  //initialize the geometry of the boxes subdirpar sets
  void gauge_sweeper_t::init_box_dir_par_geometry(int ext_gpar,int(*par_comp)(const Coords& ivol_coord,const int& dir))
  {
    comm_init_time-=take_time();
    
    //check and mark as inited, store parity and allocate geometry and subbox volume
    if(par_geom_inited) CRASH("parity checkboard already initialized");
    par_geom_inited=true;
    gpar=ext_gpar;
    ivol_of_box_dir_par=nissa_malloc("ivol_of_box_dir_par",NDIM*locVol,int);
    nsite_per_box_dir_par=nissa_malloc("nsite_per_box_dir_par",(1<<NDIM)*NDIM*gpar,int);
    
    //find the elements of all boxes
    int ibox_dir_par=0; //order in the parity-split order
    for(int ibox=0;ibox<(1<<NDIM);ibox++)
      for(int dir=0;dir<NDIM;dir++)
	for(int par=0;par<gpar;par++)
	  {
	    nsite_per_box_dir_par[par+gpar*(dir+NDIM*ibox)]=0;
	    for(int isub=0;isub<nsite_per_box[ibox];isub++)
	      {
		//get coordinates of site
		Coords isub_coord=coordOfLx(isub,boxSize[ibox]);
		
		//get coords in the local size, and parity
		Coords ivol_coord;
		for(int mu=0;mu<NDIM;mu++)
		  ivol_coord[mu]=boxSize[0][mu]*boxCoord[ibox][mu]+isub_coord[mu];
		
		int site_par=par_comp(ivol_coord,dir);
		if(site_par>=gpar||site_par<0) CRASH("obtained par %d while expecting in the range [0,%d]",par,gpar-1);
		
		//map sites in current parity
		if(site_par==par)
		  {
		    ivol_of_box_dir_par[ibox_dir_par++]=lxOfCoord(ivol_coord,locSize);
		    nsite_per_box_dir_par[par+gpar*(dir+NDIM*ibox)]++;
		  }
	      }
	  }
    comm_init_time+=take_time();
  }
  
  //check that everything is hit once and only once
  void gauge_sweeper_t::check_hit_exactly_once()
  {
    Coords *hit=nissa_malloc("hit",locVol,Coords);
    vector_reset(hit);
    
    //mark what we hit
    int ibase=0;
    for(int ibox=0;ibox<(1<<NDIM);ibox++)
      for(int dir=0;dir<NDIM;dir++)
	for(int par=0;par<gpar;par++)
	  {
	    for(int ibox_dir_par=0;ibox_dir_par<nsite_per_box_dir_par[par+gpar*(dir+NDIM*ibox)];ibox_dir_par++)
	      {
		int ivol=ivol_of_box_dir_par[ibox_dir_par+ibase];
		if(ivol>=locVol) CRASH("ibox %d ibox_par %d ibase %d ivol %d",ibox,ibox_dir_par,ibase,ivol);
		hit[ivol][dir]++;
	      }
	    ibase+=nsite_per_box_dir_par[par+gpar*(dir+NDIM*ibox)];
	  }
    
    //check to have hit everything
    for(int ivol=0;ivol<locVol;ivol++)
      for(int mu=0;mu<NDIM;mu++)
	if(hit[ivol][mu]!=1) CRASH("missing hit ivol %d mu %d",ivol,mu);
    
    nissa_free(hit);
  }
  
  //check that everything is hit without overlapping
  void gauge_sweeper_t::check_hit_in_the_exact_order()
  {
    int *hit=nissa_malloc("hit",NDIM*locVol+max_cached_link,int);
    int ibase=0;
    for(int ibox=0;ibox<(1<<NDIM);ibox++)
      for(int dir=0;dir<NDIM;dir++)
        for(int par=0;par<gpar;par++)
          {
            vector_reset(hit);
            for(int iter=0;iter<2;iter++)
              //at first iter mark all sites updating
              //at second iter check that none of them is internally used
              for(int ibox_dir_par=0;ibox_dir_par<nsite_per_box_dir_par[par+gpar*(dir+NDIM*ibox)];ibox_dir_par++)
                {
                  int ivol=ivol_of_box_dir_par[ibox_dir_par+ibase];
                  
                  if(iter==0)
                    {
                      hit[dir+NDIM*ivol]=par+1;
		      
                      if(0)
			{
			  printf("ibox %d dir %d par %d hitting %d, ivol %d{%d",
				 ibox,dir,par,
				 dir+NDIM*ivol,ivol,
				 locCoordOfLoclx[ivol][0]);
			  for(int mu=1;mu<NDIM;mu++) printf(",%d",locCoordOfLoclx[ivol][mu]);
			  printf("};%d\n",dir);
			}
                    }
                  else
                    for(int ihit=0;ihit<nlinks_per_staples_of_link;ihit++)
                      {
                        int link=ilink_per_staples[ihit+nlinks_per_staples_of_link*(ibox_dir_par+ibase)];
                        if(link>=NDIM*locVol+max_cached_link)
                          CRASH("ibox %d ibox_dir_par %d ibase %d ihit %d, link %d, max_cached_link %d",
                                ibox,ibox_dir_par,ibase,ihit,link,max_cached_link);
                        
			if(0)
			  {
			    printf("ibox %d dir %d par %d link[%d], ivol %d{%d",ibox,dir,par,ihit,ivol,locCoordOfLoclx[ivol][0]);
			    for(int mu=1;mu<NDIM;mu++) printf(",%d",locCoordOfLoclx[ivol][mu]);
			    printf("}: %d",link>=NDIM*locVol?-1:locCoordOfLoclx[link/NDIM][0]);
			    for(int mu=0;mu<NDIM;mu++) printf(",%d",link>=NDIM*locVol?-1:locCoordOfLoclx[link/NDIM][mu]);
			    printf(";%d\n",link>=NDIM*locVol?-1:link%NDIM);
			  }
			
                        if(hit[link]!=0)
			  {
			    char message[1024],*ap=message;
			    ap+=sprintf(message,"ivol %d:{%d",ivol,locCoordOfLoclx[ivol][0]);
			    for(int mu=1;mu<NDIM;mu++) ap+=sprintf(ap,",%d",locCoordOfLoclx[ivol][mu]);
			    ap+=sprintf(ap,"} ibox %d dir %d par %d got hit by %d link %d [site %d: {%d",
					ibox,dir,par,ihit,link,link/NDIM,
					link>=NDIM*locVol?-1:locCoordOfLoclx[link/NDIM][0]);
			    for(int mu=1;mu<NDIM;mu++) ap+=sprintf(ap,",%d",link>=NDIM*locVol?-1:locCoordOfLoclx[link/NDIM][mu]);
			    ap+=sprintf(ap,"},dir %d]: par %d",link%NDIM,hit[link]-1);
			    CRASH("%s",message);
			  }
		      }
                }
            ibase+=nsite_per_box_dir_par[par+gpar*(dir+NDIM*ibox)];
          }
    
    nissa_free(hit);
  }
  
  //init the list of staples
  void gauge_sweeper_t::init_staples(int ext_nlinks_per_staples_of_link,void(*ext_add_staples_per_link)(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu),void (*ext_compute_staples)(su3 staples,su3 *links,int *ilinks,double C1))
  {
    //take external nlinks and mark
    if(!par_geom_inited) CRASH("call geom initializer before");
    if(staples_inited) CRASH("staples already initialized");
    staples_inited=true;
    nlinks_per_staples_of_link=ext_nlinks_per_staples_of_link;
    add_staples_per_link=ext_add_staples_per_link;
    compute_staples=ext_compute_staples;
    
    //allocate
    ilink_per_staples=nissa_malloc("ilink_per_staples",nlinks_per_staples_of_link*NDIM*locVol,int);
    all_to_all_gathering_list_t *gl[1<<NDIM];
    for(int ibox=0;ibox<(1<<NDIM);ibox++) gl[ibox]=new all_to_all_gathering_list_t;
    add_staples_required_links(gl);
    
    //initialize the communicator
    for(int ibox=0;ibox<(1<<NDIM);ibox++)
      {
	box_comm[ibox]=new all_to_all_comm_t(*(gl[ibox]));
	delete gl[ibox];
      }
    
    //compute the maximum number of links to send and receive and allocate buffers
    for(int ibox=0;ibox<(1<<NDIM);ibox++)
      {
	max_cached_link=std::max(max_cached_link,box_comm[ibox]->nel_in);
	max_sending_link=std::max(max_sending_link,box_comm[ibox]->nel_out);
      }
    buf_out=nissa_malloc("buf_out",max_sending_link,su3);
    buf_in=nissa_malloc("buf_in",max_cached_link,su3);
    
    //check cached
    VERBOSITY_LV3_MASTER_PRINTF("Max cached links: %d\n",max_cached_link);
    if(max_cached_link>bordVol+edgeVol) CRASH("larger buffer needed [really? recheck this]");
    
    //perform two checks
    check_hit_exactly_once();
    check_hit_in_the_exact_order();
  }
  
  //ordering for link_source_dest
  int compare_link_source_dest(const void *a,const void *b)
  {return ((int*)a)[0]-((int*)b)[0];}
  
  //reorder the packer
  void reorder_packing_link_source_dest(gauge_sweeper_t* gs)
  {
    CRASH("reimplement");
    // //split workload and find starting point
    // NISSA_CHUNK_WORKLOAD(bdp_start,chunk_load,bdp_end,0,(1<<NDIM)*NDIM*gs->gpar,THREAD_ID,nthreads);
    // int ibase=0;
    // for(int bdp=0;bdp<bdp_start;bdp++) ibase+=gs->nsite_per_box_dir_par[bdp];
    
    // for(int bdp=bdp_start;bdp<bdp_end;bdp++)
    //   {
    // 	//sort and increase the base
    // 	qsort(gs->packing_link_source_dest+2*gs->nlinks_per_staples_of_link*ibase,
    // 	      gs->nlinks_per_staples_of_link*gs->nsite_per_box_dir_par[bdp],
    // 	      2*sizeof(int),
    // 	      compare_link_source_dest);
    // 	ibase+=gs->nsite_per_box_dir_par[bdp];
    //   }
  }
  
  //find the place where each link must be copied to access it sequentially
  void gauge_sweeper_t::find_packing_index(void (*ext_compute_staples_packed)(su3 staples,su3 *links,double C1))
  {
    if(!packing_inited)
      {
	//mark packing to be have been inited and allocate
	packing_inited=true;
	compute_staples_packed=ext_compute_staples_packed;
	packing_link_source_dest=nissa_malloc("packing_link_source_dest",2*(nlinks_per_staples_of_link*NDIM*locVol+1),int);
	
	int ibase=0;
	int max_packing_link_nel=0;
	for(int ibox=0;ibox<(1<<NDIM);ibox++)
	  for(int dir=0;dir<NDIM;dir++)
	    for(int par=0;par<gpar;par++)
	      {
		//find the packing size
		int ns=nsite_per_box_dir_par[par+gpar*(dir+NDIM*ibox)],nsh=ns/2; //only for bg/q
		if(nsh*2!=ns) nsh++;
		max_packing_link_nel=std::max(max_packing_link_nel,nlinks_per_staples_of_link*2*nsh);
		
		//scan the destination
		for(int ibox_dir_par=0;ibox_dir_par<ns;ibox_dir_par++)
		  for(int ilink=0;ilink<nlinks_per_staples_of_link;ilink++)
		    {
		      packing_link_source_dest[2*(ilink+nlinks_per_staples_of_link*(ibox_dir_par+ibase))+0]=
			ilink_per_staples[(ibox_dir_par+ibase)*nlinks_per_staples_of_link+ilink];
		      packing_link_source_dest[2*(ilink+nlinks_per_staples_of_link*(ibox_dir_par+ibase))+1]=
			ilink+nlinks_per_staples_of_link*ibox_dir_par
			;
		    }
		
		//increase the base
		ibase+=nsite_per_box_dir_par[par+gpar*(dir+NDIM*ibox)];
	      }
	
	reorder_packing_link_source_dest(this);
	
	//allocate packing link and deallocate ilink_per_staples
	packing_link_buf=nissa_malloc("packing_link_buf",max_packing_link_nel,su3);
	nissa_free(ilink_per_staples);
      }
  }
  
  //add all the links needed to compute staples separately for each box
  void add_staples_required_links_to_gauge_sweep(gauge_sweeper_t *gs,all_to_all_gathering_list_t **gl)
  {
    for(int ibox=0;ibox<(1<<NDIM);ibox++)
      {
	//find base for curr box
	int ibase=0;
	for(int jbox=0;jbox<ibox;jbox++) ibase+=nsite_per_box[jbox];
	ibase*=NDIM;
	
	//scan all the elements of sub-box, selecting only those with the good parity
	for(int dir=0;dir<NDIM;dir++)
	  for(int par=0;par<gs->gpar;par++)
	    {
	      for(int ibox_dir_par=ibase;ibox_dir_par<ibase+gs->nsite_per_box_dir_par[par+gs->gpar*(dir+NDIM*ibox)];ibox_dir_par++)
		{
		  int ivol=gs->ivol_of_box_dir_par[ibox_dir_par];
		  gs->add_staples_per_link(gs->ilink_per_staples+ibox_dir_par*gs->nlinks_per_staples_of_link,*(gl[ibox]),ivol,dir);
		}
	      ibase+=gs->nsite_per_box_dir_par[par+gs->gpar*(dir+NDIM*ibox)];
	    }
    }
  }
  
  //wrapper to use threads
  void gauge_sweeper_t::add_staples_required_links(all_to_all_gathering_list_t **gl)
  {add_staples_required_links_to_gauge_sweep(this,gl);}
  
  //pack all the links required to compute staples
  void gauge_sweeper_t::pack_links(quad_su3 *conf,int ibase,int nbox_dir_par)
  {
    CRASH("redo");
    // //prepare the chunk load
    // NISSA_CHUNK_WORKLOAD(start,chunk_load,end,0,nlinks_per_staples_of_link*nbox_dir_par,THREAD_ID,NACTIVE_THREADS);
    // int *source_dest=packing_link_source_dest+2*(nlinks_per_staples_of_link*ibase+start);
    
    // for(int ilink_to_ship=start;ilink_to_ship<end;ilink_to_ship++)
    //   {
    // 	int isource=*(source_dest++);
    // 	int idest=*(source_dest++);
    // 	su3_copy(packing_link_buf[idest],((su3*)conf)[isource]);
    //   }
    // THREAD_BARRIER();
  }
  
  //compute the parity according to the Symanzik requirements
  int Symanzik_par(const Coords& ivol_coord,const int& dir)
  {
    int site_par=0;
    for(int mu=0;mu<NDIM;mu++) site_par+=((mu==dir)?2:1)*ivol_coord[mu];
    
    site_par=site_par%NDIM;
    
    return site_par;
  }
  
  //add all links needed for a certain site
  void add_Symanzik_staples(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu)
  {
    Coords& A=glbCoordOfLoclx[ivol];                            //       P---O---N
    Coords B,C,/*D,*/E,F,G,H,I,J,K,L,M,/*N,*/O,P;               //       |   |   |
    //find coord mu                                             //   H---G---F---E---D
    K[mu]=L[mu]=M[mu]=(A[mu]-1+glbSize[mu])%glbSize[mu];      //   |   |   |   |   |
    I[mu]=J[mu]=B[mu]=C[mu]=A[mu];                              //   I---J---A---B---C
    /*D[mu]=*/E[mu]=F[mu]=G[mu]=H[mu]=(A[mu]+1)%glbSize[mu];   //       |   |   |
    /*N[mu]=*/O[mu]=P[mu]=(A[mu]+2)%glbSize[mu];               //       K---L---M
    for(int inu=0;inu<NDIM-1;inu++)
      {
	int nu=perpDirs[mu][inu];
	
	//copy orthogonal coords
#if NDIM>=3
	for(int irh=0;irh<NDIM-2;irh++)
	  {
	    int rh=perp2Dirs[mu][inu][irh];
	    B[rh]=C[rh]=/*D[rh]=*/E[rh]=F[rh]=G[rh]=H[rh]=I[rh]=J[rh]=K[rh]=L[rh]=M[rh]=/*N[rh]=*/O[rh]=P[rh]=A[rh];
	  }
#endif
	
	//find coord nu
	H[nu]=I[nu]=(A[nu]-2+glbSize[nu])%glbSize[nu];
	K[nu]=J[nu]=G[nu]=P[nu]=(I[nu]+1)%glbSize[nu];
	L[nu]=F[nu]=O[nu]=A[nu];
	M[nu]=B[nu]=E[nu]=/*N[nu]=*/(A[nu]+1)%glbSize[nu];
	C[nu]=/*D[nu]=*/(A[nu]+2)%glbSize[nu];
	
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
    
    //8*(NDIM-1) link missing, 2 readded = -22 links (in 4d)
  }
  
  //compute the summ of the staples pointed by "ilinks"
  void compute_Symanzik_staples(su3 staples,su3 *links,int *ilinks,double C1)
  {
    su3 squares,rectangles,up_rectangles,dw_rectangles;
    su3_put_to_zero(squares);
    su3_put_to_zero(rectangles);
    su3_put_to_zero(up_rectangles);
    su3_put_to_zero(dw_rectangles);
    
    const int PARTIAL=2;
    for(int inu=0;inu<NDIM-1;inu++)
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
    su3_linear_comb(staples,squares,get_C0(C1),rectangles,C1);
  }
  
  //compute the summ of the staples pointed by "ilinks"
  void compute_Symanzik_staples_packed(su3 staples,su3 *links,double C1)
  {
    double C0=get_C0(C1);
    
    su3 squares,rectangles,up_rectangles,dw_rectangles;
    su3_put_to_zero(squares);
    su3_put_to_zero(rectangles);
    su3_put_to_zero(up_rectangles);
    su3_put_to_zero(dw_rectangles);
    
    const int PARTIAL=2;
    for(int inu=0;inu<NDIM-1;inu++)
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
    su3_linear_comb(staples,squares,C0,rectangles,C1);
  }
  
  //initialize the tlSym sweeper using the above defined routines
  void init_Symanzik_sweeper()
  {
    if(!Symanzik_sweeper->staples_inited)
      {
	VERBOSITY_LV3_MASTER_PRINTF("Initializing Symanzik sweeper\n");
	//checking consistency for gauge_sweeper initialization
	for(int mu=0;mu<NDIM;mu++) if(locSize[mu]<4) CRASH("loc_size[%d]=%d must be at least 4",mu,locSize[mu]);
	//initialize the Symanzik sweeper
	const int nlinks_per_Symanzik_staples_of_link=(NDIM-1)*2*(3+5*3)-(NDIM-1)*8+2;
	Symanzik_sweeper->init_box_dir_par_geometry(4,Symanzik_par);
	Symanzik_sweeper->init_staples(nlinks_per_Symanzik_staples_of_link,add_Symanzik_staples,compute_Symanzik_staples);
      }
  }
  
  ///////////////////////////////////////// Wilson ////////////////////////////////////////
  
  //compute the parity according to the Wilson requirements
  int Wilson_par(const Coords& ivol_coord,const int& dir)
  {
    int site_par=0;
    for(int mu=0;mu<NDIM;mu++) site_par+=ivol_coord[mu];
    
    site_par=site_par%2;
    
    return site_par;
  }
  
  //add all links needed for a certain site
  void add_Wilson_staples(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu)
  {
    Coords& A=glbCoordOfLoclx[ivol];
    Coords B,F,G,J;                                       //       G---F---E
    //find coord mu                                         //       |   |   |
    J[mu]=B[mu]=A[mu];                                      //       J---A---B
    F[mu]=G[mu]=(A[mu]+1)%glbSize[mu];
    for(int inu=0;inu<NDIM-1;inu++)
      {
	int nu=perpDirs[mu][inu];
	
	//copy orthogonal coords
#if NDIM>=3
	for(int irh=0;irh<NDIM-2;irh++)
	  {
	    int rh=perp2Dirs[mu][inu][irh];
	    B[rh]=F[rh]=G[rh]=J[rh]=A[rh];
	  }
#endif
	
	//find coord nu
	J[nu]=G[nu]=(A[nu]-1+glbSize[nu])%glbSize[nu];
	F[nu]=A[nu];
	B[nu]=(A[nu]+1)%glbSize[nu];
	
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
  void compute_Wilson_staples(su3 staples,su3 *links,int *ilinks,double C1)
  {
    su3_put_to_zero(staples);
    
    for(int inu=0;inu<NDIM-1;inu++)
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
  
  void compute_Wilson_staples_packed(su3 staples,su3 *links,double C1)
  {
    su3_put_to_zero(staples);
    
    for(int inu=0;inu<NDIM-1;inu++)
      {
	su3 hb;
	//backward square staple
	unsafe_su3_dag_prod_su3(hb,links[ 0],links[ 1]);
	su3_summ_the_prod_su3(staples,hb,links[ 2]);
	//forward square staple
	su3 hf;
	unsafe_su3_prod_su3(hf,links[ 3],links[ 4]);
	su3_summ_the_prod_su3_dag(staples,hf,links[ 5]);
	
	links+=6;
      }
  }
  
  //initialize the Wilson sweeper using the above defined routines
  void init_Wilson_sweeper()
  {
    if(!Wilson_sweeper->staples_inited)
      {
	VERBOSITY_LV3_MASTER_PRINTF("Initializing Wilson sweeper\n");
	//checking consistency for gauge_sweeper initialization
	for(int mu=0;mu<NDIM;mu++) if(locSize[mu]<2) CRASH("loc_size[%d]=%d must be at least 2",mu,locSize[mu]);
	//initialize the Wilson sweeper
	Wilson_sweeper->init_box_dir_par_geometry(2,Wilson_par);
	const int nlinks_per_Wilson_staples_of_link=6*(NDIM-1);
	Wilson_sweeper->init_staples(nlinks_per_Wilson_staples_of_link,add_Wilson_staples,compute_Wilson_staples);
      }
  }
  
  //call the appropriate sweeper intializator
  void init_sweeper(gauge_action_name_t gauge_action_name)
  {
    switch(gauge_action_name)
      {
      case WILSON_GAUGE_ACTION:
	if(!Wilson_sweeper->staples_inited)
	  init_Wilson_sweeper();
	break;
      case TLSYM_GAUGE_ACTION:
      case IWASAKI_GAUGE_ACTION:
	if(!Symanzik_sweeper->staples_inited)
	  init_Symanzik_sweeper();
	break;
      case UNSPEC_GAUGE_ACTION:
	CRASH("unspecified action");
	break;
      default:
	CRASH("not implemented action");
	break;
      }
  }
}
