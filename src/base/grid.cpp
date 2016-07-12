#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_MPI
 #include <mpi.h>
#endif

#include "base/grid.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_vir.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //compute the variance of the border
  int compute_border_variance(coords L)
  {
    int S2B=0,SB=0;
    for(int ib=0;ib<NDIM;ib++)
      {
	int B=1;
	for(int mu=0;mu<NDIM;mu++) if(mu!=ib) B*=L[mu];
	SB+=B;
	S2B+=B*B;
      }
    SB/=NDIM;
    S2B/=NDIM;
    S2B-=SB*SB;
    
    return S2B;
  }
  
  //if nfact_V>=nfact_R factorize the number of ranks, otherwise the volume
  //in the first case we find the best way to assign the ranks to different directions
  //in the second case we find how many sites per direction to assign to each ranks
  partitioning_t::partitioning_t(long long int V,int R)
  {
    std::vector<int> list_fact_V=factorize(V);
    std::vector<int> list_fact_R=factorize(R);
    factorize_R=(list_fact_V.size()>=list_fact_R.size());
    if(factorize_R) list_fact=list_fact_R;
    else            list_fact=list_fact_V;
    
    //compute the number of combinations: this is given by NDIM^nfact
    if(list_fact.size())
      {
	ncombo=1;
	for(size_t ifact=0;ifact<list_fact.size();ifact++) ncombo*=NDIM;
      }
    else ncombo=0;
    
    //restart combo
    restart();
  }
  
  int partitioning_t::partitioning_t::decrypt_and_validate_partition(coords R_per_dir,coords grid_size,coords min_grid_size,coords fix_R_per_dir)
  {
    //reset loc
    for(int mu=0;mu<NDIM;mu++) R_per_dir[mu]=1;
    
    //special case
    if(ncombo==0) return 0;
    
    //compute mask factor
    int mask=1;
    for(int jfact=0;jfact<(int)list_fact.size()-1;jfact++) mask*=NDIM;
    
    //find the partioning corresponding to icombo
    int ifact=list_fact.size()-1;
    int valid_partitioning=1;
    while(valid_partitioning && ifact>=0)
      {
	//find the direction: this is given by the ifact digit of icombo wrote in base NDIM
	int mu=(icombo/mask)%NDIM;
	
	//if we are factorizing the number of ranks, factor is given by list_fact, otherwise from glb[mu]/list_fact
	if(factorize_R) R_per_dir[mu]*=list_fact[ifact];
	else            R_per_dir[mu]*=grid_size[mu]/list_fact[ifact];
	
	//check that the total volume is a multiple of what trying to divide it for
	if(valid_partitioning && grid_size[mu]%R_per_dir[mu]!=0)
	  {
	    valid_partitioning=false;
	    master_printf("grid_size(%d)%%R_per_dir(%d) mu=%d !=0\n",grid_size[mu],R_per_dir[mu],mu);
	  }
	//check that it is larger than the minimal size per this dir
	if(valid_partitioning && grid_size[mu]<R_per_dir[mu]*min_grid_size[mu])
	  {
	    valid_partitioning=false;
	    master_printf("grid_size(%d)<R_per_dir*min_grid_size(%d*%d) mu=%d !=0\n",grid_size[mu],R_per_dir[mu],min_grid_size[mu],mu);
	  }
	
	//if fixed R per this dir
	if(valid_partitioning && fix_R_per_dir[mu]>0)
	  {
	    //check that it is not larger than the agreed size
	    valid_partitioning&=(R_per_dir[mu]<=fix_R_per_dir[mu]);
	    //check that it divides fix_R_per_dir
	    valid_partitioning&=(fix_R_per_dir[mu]%R_per_dir[mu]==0);
	  }
	
	//if partition is valid move to next factor
	if(valid_partitioning)
	  {
	    ifact--;
	    mask/=NDIM;
	  }
      }
    
    return ifact;
  }
  
  //
  void partitioning_t::skip_combo_of_factors(int ifact)
  {
    int mask=1;
    for(int jfact=0;jfact<ifact-1;jfact++) mask*=NDIM;
    icombo+=mask;
  }
  
  //
  bool partitioning_t::find_next_valid_partition(coords R_per_dir,coords grid_size,coords min_grid_size,coords fix_R_per_dir)
  {
    icombo++;
    
    int nunassign_factors;
    if(!finished())
      do
	{
	  nunassign_factors=decrypt_and_validate_partition(R_per_dir,grid_size,min_grid_size,fix_R_per_dir);
	  master_printf("icombo %d unassigned %d ncombo %d\n",icombo,nunassign_factors,ncombo);
	  if(nunassign_factors>0) skip_combo_of_factors(nunassign_factors);
	}
      while(nunassign_factors>0 && !finished());
    
    return !finished();
  }
  
  //find the grid minimizing the surface
  void find_grid()
  {
    //check that we did not ask to fix ranks and vranks in an impossible way
    int res_nranks=nranks;
    int res_nvranks=nvranks_max;
    coords min_loc_size;
    coords min_vir_loc_size;
    for(int mu=0;mu<NDIM;mu++)
      {
	min_loc_size[mu]=1;
	if(use_eo_geom) min_loc_size[mu]*=2;
	min_vir_loc_size[mu]=min_loc_size[mu];
	if(fix_nvranks_max[mu]) min_vir_loc_size[mu]*=fix_nvranks_max[mu];
	
	if(fix_nranks[mu])
	  {
	    if(glb_size[mu]%fix_nranks[mu]) CRASH("asked to fix nranks to %d in dir %d in an impossible way, global size %d min_loc_size %d",fix_nranks[mu],mu,glb_size[mu],min_loc_size[mu]);
	    if(fix_nranks[mu]>nranks) CRASH("asked to fix nranks in dir %d to a number %d larger than the total ranks, %d",mu,fix_nranks[mu],nranks);
	    res_nranks/=fix_nranks[mu];
	  }
	
	if(fix_nvranks_max[mu])
	  {
	    if(glb_size[mu]%fix_nvranks_max[mu]||fix_nvranks_max[mu]>glb_size[mu]) CRASH("asked to fix nvranks to %d in dir %d in an impossible way, global size %d",fix_nvranks_max[mu],mu,glb_size[mu]);
	    if(fix_nvranks_max[mu]>nvranks_max) CRASH("asked to fix nvranks in dir %d to a number %d larger than the maximal vnranks, %d",mu,fix_nvranks_max[mu],nvranks_max);
	    res_nvranks/=fix_nvranks_max[mu];
	  }
      }
    if(res_nranks<1) CRASH("overfixed the number of ranks per direction");
    if(res_nvranks<1) CRASH("overfixed the number of vranks per direction");
    
    //////////////////// find the partitioning which minimize the surface /////////////////////
    
    int something_found=0;
    
    //find the partition which minimize the surface and the surface variance
    double min_rel_surf=-1;
    
    partitioning_t ranks_partitioning(glb_vol,nranks);
    master_printf("partitioning glb_vol %d in %d ranks obtained %u possible combos\n",glb_vol,nranks,ranks_partitioning.ncombo);
    partitioning_t vranks_partitioning(loc_vol,nvranks_max);
    master_printf("partitioning lc_vol %d in %d vranks obtained %u possible combos\n",loc_vol,nvranks_max,vranks_partitioning.ncombo);
    
    coords RPD;
    //use min_vloc_size because we want to further partition that dir
    while(ranks_partitioning.find_next_valid_partition(RPD,glb_size,min_vir_loc_size,fix_nranks))
      {
	//set the local size
	master_printf("%d %d %d %d\n",RPD[0],RPD[1],RPD[2],RPD[3]);
	coords LS;
	for(int mu=0;mu<NDIM;mu++) LS[mu]=glb_size[mu]/RPD[mu];
	vranks_partitioning.restart();
	
	coords VPD;
	coords fix_nvranks;
	memset(fix_nvranks,0,sizeof(coords));
	while(vranks_partitioning.find_next_valid_partition(VPD,LS,min_loc_size,fix_nvranks))
	  {
	    //set the vir local size
	    coords VLS;
	    for(int mu=0;mu<NDIM;mu++) VLS[mu]=LS[mu]/VPD[mu];
	    
	    //compute the rank surface=loc_vol-bulk_volume
	    coords loc_bulk_size,vloc_bulk_size;
	    int loc_bulk_vol=1,vloc_bulk_vol=1;
	    for(int mu=0;mu<NDIM;mu++)
	      {
		loc_bulk_size[mu]=LS[mu];
		vloc_bulk_size[mu]=VLS[mu];
		if(RPD[mu]>1)
		  {
		    loc_bulk_size[mu]-=2;
		    vloc_bulk_size[mu]-=2;
		  }
		else if(VPD[mu]>1) vloc_bulk_size[mu]-=2;
		loc_bulk_vol*=loc_bulk_size[mu];
		vloc_bulk_vol*=vloc_bulk_size[mu];
	      }
	    loc_bulk_vol=std::max(0,loc_bulk_vol);
	    vloc_bulk_vol=std::max(0,vloc_bulk_vol);
	    
	    //compute relative surface
	    double bulk_rel_surf=(double)loc_bulk_vol/loc_vol;
	    double vbulk_rel_surf=(double)vloc_bulk_vol/vloc_min_vol;
	    double rel_surf=sqrt(sqr(bulk_rel_surf)+sqr(vbulk_rel_surf));
	    
	    //if it is the minimal surface (or first valid combo) copy it and compute the border size
	    if(rel_surf<min_rel_surf||min_rel_surf<0)
	      {
		min_rel_surf=rel_surf;
		
		for(int mu=0;mu<NDIM;mu++)
		  {
		    nranks_per_dir[mu]=RPD[mu];
		    vloc_min_size[mu]=loc_size[mu]/VPD[mu];
		  }
		
		something_found=1;
	      }
	  }
      }
    
    if(!something_found) CRASH("no valid partitioning found");
  }
  
  //define boxes
  void init_boxes()
  {
    //get the size of box 0
    for(int mu=0;mu<NDIM;mu++)
      {
	if(loc_size[mu]<2) CRASH("loc_size[%d]=%d must be at least 2",mu,loc_size[mu]);
	box_size[0][mu]=loc_size[mu]/2;
      }
    
    //get coords of cube ans box size
    coords nboxes;
    for(int mu=0;mu<NDIM;mu++) nboxes[mu]=2;
    for(int ibox=0;ibox<(1<<NDIM);ibox++)
      {
	//coords
	verbosity_lv3_master_printf("Box %d coord [ ",ibox);
	coord_of_lx(box_coord[ibox],ibox,nboxes);
	for(int mu=0;mu<NDIM;mu++) verbosity_lv3_master_printf("%d ",box_coord[ibox][mu]);
      
	//size
	verbosity_lv3_master_printf("] size [ ",ibox);
	nsite_per_box[ibox]=1;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    if(ibox!=0) box_size[ibox][mu]=((box_coord[ibox][mu]==0)?
					    (box_size[0][mu]):(loc_size[mu]-box_size[0][mu]));
	    nsite_per_box[ibox]*=box_size[ibox][mu];
	    verbosity_lv3_master_printf("%d ",box_size[ibox][mu]);
	  }
	verbosity_lv3_master_printf("], nsites: %d\n",nsite_per_box[ibox]);
      }
  }
  
  //initialize MPI grid
  //if you need non-homogeneus glb_size[i] pass L=T=0 and
  //set glb_size before calling the routine
  void init_grid(int T,int L)
  {
    //take initial time
    double time_init=-take_time();
    master_printf("\nInitializing grid, geometry and communications\n");
    
    if(grid_inited==1) CRASH("grid already intialized!");
    grid_inited=1;
    
    //set the volume
    if(T!=0 && L!=0)
      {
	glb_size[0]=T;
	for(int mu=1;mu<NDIM;mu++) glb_size[mu]=L;
      }
    
    //broadcast the global sizes
    coords_broadcast(glb_size);
    
    //calculate global volume
    glb_vol=1;
    for(int mu=0;mu<NDIM;mu++) glb_vol*=glb_size[mu];
    glb_volh=glb_vol/2;
    glb_spat_vol=glb_vol/glb_size[0];
    glb_vol2=(double)glb_vol*glb_vol;
    if(use_eo_geom) if(glb_vol%(1<<NDIM)!=0) CRASH("in order to use eo geometry, global volume %d must be a multiple of (1<<NDIM)=%d",glb_vol,(1<<NDIM));
    
    //print info on global lattice
    master_printf("Global lattice:\t%d",glb_size[0]);
    for(int mu=1;mu<NDIM;mu++) master_printf("x%d",glb_size[mu]);
    master_printf(" = %d\n",glb_vol);
    
    //initialize local volume and check it
    master_printf("Number of MPI ranks: %d\n",nranks);
    loc_vol=glb_vol/nranks;
    loc_volh=loc_vol/2;
    if(loc_vol*nranks!=glb_vol) CRASH("global volume %d cannot be divided by nranks %d",glb_vol,nranks);
    if(use_eo_geom) if(loc_vol%(1<<NDIM)!=0) CRASH("in order to use eo geometry, local volume %d must be a multiple of (1<<NDIM)=%d",loc_vol,(1<<NDIM));
    
    //initialize virtual local volume and check it
    master_printf("Number of maximal virtual ranks: %d\n",nvranks_max);
    vloc_min_vol=loc_vol/nvranks_max;
    if(vloc_min_vol*nvranks_max!=loc_vol) CRASH("local volume %d cannot be divided by nvranks %d",loc_vol,nvranks_max);
    if(use_eo_geom) if(vloc_min_vol%(1<<NDIM)) CRASH("in order to use eo geometry and virtual ranks, local virtual volume %d must be a multiple of (1<<NDIM)=%d",vloc_min_vol,(1<<NDIM));
    
    //find the grid minimizing the surface
    find_grid();
    
    //check whether each direction dir is parallelized
    for(int mu=0;mu<NDIM;mu++)
      {
	paral_dir[mu]=(nranks_per_dir[mu]>1);
	nparal_dir+=paral_dir[mu];
      }
    
    //creates the grid
    master_printf("Creating MPI cartesian grid:\t%d",nranks_per_dir[0]);
    for(int mu=1;mu<NDIM;mu++) master_printf("x%d",nranks_per_dir[mu]);
    master_printf("\n");
    create_MPI_cartesian_grid();
    
    //compute the local volume
    for(int mu=0;mu<NDIM;mu++) loc_size[mu]=glb_size[mu]/nranks_per_dir[mu];
    loc_spat_vol=loc_vol/loc_size[0];
    loc_vol2=(double)loc_vol*loc_vol;
    
    //calculate bulk size
    bulk_vol=non_bw_surf_vol=1;
    for(int mu=0;mu<NDIM;mu++)
      if(paral_dir[mu])
	{
	  bulk_vol*=loc_size[mu]-2;
	  non_bw_surf_vol*=loc_size[mu]-1;
	}
      else
	{
	  bulk_vol*=loc_size[mu];
	  non_bw_surf_vol*=loc_size[mu];
	}
    non_fw_surf_vol=non_bw_surf_vol;
    fw_surf_vol=bw_surf_vol=loc_vol-non_bw_surf_vol;
    surf_vol=loc_vol-bulk_vol;
    
    //calculate the border size
    bord_volh=0;
    bord_offset[0]=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	//bord size along the mu dir
	if(paral_dir[mu]) bord_dir_vol[mu]=loc_vol/loc_size[mu];
	else bord_dir_vol[mu]=0;
	
	//total bord
	bord_volh+=bord_dir_vol[mu];
	
	//summ of the border extent up to dir mu
	if(mu>0) bord_offset[mu]=bord_offset[mu-1]+bord_dir_vol[mu-1];
      }
    bord_vol=2*bord_volh;
    
    init_boxes();
    
#ifdef USE_VNODES
    //two times the size of vnode_paral_dir face
    vbord_vol=2*loc_vol/loc_size[vnode_paral_dir]; //so is not counting all vsites
    vbord_volh=vbord_vol/2;
    vdir_bord_vol=vbord_vol/2;
    vdir_bord_volh=vbord_volh/2;
    //compute the offset between sites of different vnodes
    //this amount to the product of the local size of the direction running faster than
    //vnode_paral_dir, times half the local size along vnode_paral_dir
    vnode_lx_offset=loc_size[vnode_paral_dir]/NVNODES;
    for(int mu=vnode_paral_dir+1;mu<NDIM;mu++) vnode_lx_offset*=loc_size[mu];
    vnode_eo_offset=vnode_lx_offset/2;
#endif
    
    //calculate the egdes size
    edge_vol=0;
    edge_offset[0]=0;
    int iedge=0;
    for(int mu=0;mu<NDIM;mu++)
      for(int nu=mu+1;nu<NDIM;nu++)
	{
	  //edge among the i and j dir
	  if(paral_dir[mu] && paral_dir[nu]) edge_dir_vol[iedge]=bord_dir_vol[mu]/loc_size[nu];
	  else edge_dir_vol[iedge]=0;
	  
	  //total edge
	  edge_vol+=edge_dir_vol[iedge];
	  
	  //summ of the border extent up to dir i
	  if(iedge>0)
	  edge_offset[iedge]=edge_offset[iedge-1]+edge_dir_vol[iedge-1];
	  iedge++;
	}
    edge_vol*=4;
    edge_volh=edge_vol/2;
    master_printf("Edge vol: %d\n",edge_vol);
    
    //set edge numb
    {
      int iedge=0;
      for(int mu=0;mu<NDIM;mu++)
	{
	  edge_numb[mu][mu]=-1;
	  for(int nu=mu+1;nu<NDIM;nu++)
	    {
	      edge_numb[mu][nu]=edge_numb[nu][mu]=iedge;
	      iedge++;
	    }
	}
    }
    
    //print information on local volume
    master_printf("Local volume\t%d",loc_size[0]);
    for(int mu=1;mu<NDIM;mu++) master_printf("x%d",loc_size[mu]);
    master_printf(" = %d\n",loc_vol);
    master_printf("List of parallelized dirs:\t");
    for(int mu=0;mu<NDIM;mu++) if(paral_dir[mu]) master_printf("%d ",mu);
    if(nparal_dir==0) master_printf("(none)");
    master_printf("\n");
    master_printf("Border size: %d\n",bord_vol);
    for(int mu=0;mu<NDIM;mu++)
      verbosity_lv3_master_printf("Border offset for dir %d: %d\n",mu,bord_offset[mu]);
    
    //print orderd list of the rank names
    if(VERBOSITY_LV3)
      {
	char proc_name[1024];
	int proc_name_length;
	MPI_Get_processor_name(proc_name,&proc_name_length);
	
	for(int irank=0;irank<nranks;irank++)
	  {
	    if(irank==rank)
	      {
		printf("Rank %d of %d running on processor %s: %d (%d",rank,nranks,proc_name,cart_rank,rank_coord[0]);
		for(int mu=1;mu<NDIM;mu++) printf(" %d",rank_coord[mu]);
		printf(")\n");
	      }
	    fflush(stdout);
	    ranks_barrier();
	    MPI_Barrier(MPI_COMM_WORLD);
	  }
      }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    
    //set the cartesian and e/o and vir geometry
    set_lx_geometry();
    if(use_eo_geom) set_eo_geometry();
    if(use_vranks) set_vranks_geometry();
    
    ///////////////////////////////////// start communicators /////////////////////////////////
    
    ncomm_allocated=0;
    
    //allocate only now buffers, so we should have finalized its size
    recv_buf=nissa_malloc("recv_buf",recv_buf_size,char);
    send_buf=nissa_malloc("send_buf",send_buf_size,char);
    
    //setup all lx borders communicators
    set_lx_comm(lx_su3_comm,sizeof(su3));
    set_lx_comm(lx_quad_su3_comm,sizeof(quad_su3));
    set_lx_comm(lx_single_quad_su3_comm,sizeof(single_quad_su3));
    set_lx_comm(lx_as2t_su3_comm,sizeof(as2t_su3));
    set_lx_comm(lx_spin_comm,sizeof(spin));
    set_lx_comm(lx_color_comm,sizeof(color));
    set_lx_comm(lx_single_color_comm,sizeof(single_color));
    set_lx_comm(lx_spinspin_comm,sizeof(spinspin));
    set_lx_comm(lx_spin1field_comm,sizeof(spin1field));
    set_lx_comm(lx_spincolor_comm,sizeof(spincolor));
    set_lx_comm(lx_single_halfspincolor_comm,sizeof(single_halfspincolor));
    set_lx_comm(lx_spincolor_128_comm,sizeof(spincolor_128));
    set_lx_comm(lx_halfspincolor_comm,sizeof(halfspincolor));
    set_lx_comm(lx_colorspinspin_comm,sizeof(colorspinspin));
    set_lx_comm(lx_su3spinspin_comm,sizeof(su3spinspin));
    
    //setup all lx edges communicators
#ifdef USE_MPI
    set_lx_edge_senders_and_receivers(MPI_LX_SU3_EDGES_SEND,MPI_LX_SU3_EDGES_RECE,&MPI_SU3);
    set_lx_edge_senders_and_receivers(MPI_LX_QUAD_SU3_EDGES_SEND,MPI_LX_QUAD_SU3_EDGES_RECE,&MPI_QUAD_SU3);
    set_lx_edge_senders_and_receivers(MPI_LX_AS2T_SU3_EDGES_SEND,MPI_LX_AS2T_SU3_EDGES_RECE,&MPI_AS2T_SU3);
#endif
    
    if(use_eo_geom)
      {
	set_eo_comm(eo_spin_comm,sizeof(spin));
	set_eo_comm(eo_spincolor_comm,sizeof(spincolor));
	set_eo_comm(eo_spincolor_128_comm,sizeof(spincolor_128));
	set_eo_comm(eo_color_comm,sizeof(color));
	set_eo_comm(eo_single_color_comm,sizeof(single_color));
	set_eo_comm(eo_halfspincolor_comm,sizeof(halfspincolor));
	set_eo_comm(eo_single_halfspincolor_comm,sizeof(single_halfspincolor));
	set_eo_comm(eo_quad_su3_comm,sizeof(quad_su3));
	set_eo_comm(eo_su3_comm,sizeof(su3));
	
#ifdef USE_MPI
	set_eo_edge_senders_and_receivers(MPI_EO_QUAD_SU3_EDGES_SEND,MPI_EO_QUAD_SU3_EDGES_RECE,&MPI_QUAD_SU3);
#endif
      }
    
    //take final time
    master_printf("Time elapsed for grid inizialization: %f s\n",time_init+take_time());
  }
}
