/*
  This program can be used to generate gauge configurations
  according to the rooted-staggered action in presence of 
  electromagnetic fields and/or imaginary chemical potentials.
  
  The molecular dynamic routines are in the file:
   ../../src/hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.cpp
*/

#include <math.h>

#include "nissa.hpp"

using namespace nissa;

//observables
char gauge_obs_path[1024];

//input and output path for confs
char conf_path[1024];
char store_conf_path[1024];

//conf and staples
quad_su3 *conf;
su3 *buf_out,*buf_in;

//evol pars
double beta;
gauge_action_name_t gauge_action_name;
pure_gauge_evol_pars_t evol_pars;

//confs
int nprod_confs;
int iconf,max_nconfs;
int store_conf_each;
int store_running_temp_conf;
int seed;

//communications and building
int max_cached_link=0,max_sending_link=0;
int nlinks_per_paths_site=3*2*(3+5*3);
coords box_size[16];
int nsite_per_box[16];
int nsite_per_box_dir_par[16*4*4];
int *ilink_per_paths;
int *ivol_of_box_dir_par;
all_to_all_comm_t *box_comm[16];

//bench
double comm_time=0,comp_time=0,meas_time=0,read_time=0,write_time=0,base_init_time=0;

enum start_conf_cond_t{UNSPEC_COND,HOT,COLD};
enum update_alg_t{UNSPEC_UP,HEAT,OVER};

void measure_gauge_obs(bool );

//write a conf adding info
void write_conf()
{
  write_time-=take_time();
  
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  char text[1024];

  //conf id
  sprintf(text,"%d",iconf);
  ILDG_string_message_append_to_last(&mess,"ConfID",text);
  
  //skip 10 random numbers
  for(int iskip=0;iskip<10;iskip++) rnd_get_unif(&glb_rnd_gen,0,1);

  //glb_rnd_gen status
  convert_rnd_gen_to_text(text,&glb_rnd_gen);
  ILDG_string_message_append_to_last(&mess,"RND_gen_status",text);
  
  //write the conf
  write_ildg_gauge_conf(conf_path,conf,64,&mess);
  
  //free messages
  ILDG_message_free_all(&mess);

  write_time+=take_time();
}

//read conf
void read_conf()
{
  read_time-=take_time();
  
  //init messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);

  //read the conf
  read_ildg_gauge_conf(conf,conf_path,&mess);
  
  //scan messages
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    {  
      if(strcasecmp(cur_mess->name,"ConfID")==0) sscanf(cur_mess->data,"%d",&iconf);
      if(strcasecmp(cur_mess->name,"RND_gen_status")==0) start_loc_rnd_gen(cur_mess->data);
    }
  
  //if message with string not found start from input seed
  if(loc_rnd_gen_inited==0)
    {
      master_printf("RND_gen status not found inside conf, starting from input passed seed\n");
      start_loc_rnd_gen(seed);
    }
  
  //free all messages
  ILDG_message_free_all(&mess);

  read_time+=take_time();
}

//add link to the map if needed
int add_link(all_to_all_gathering_list_t &gat,coords g,int mu)
{
  //find rank and local position
  int ivol,irank;
  get_loclx_and_rank_of_coord(&ivol,&irank,g);
  int ilink_asked=4*ivol+mu;
  
  //if it is local, return local position
  if(irank==rank) return ilink_asked;
  else
    {
      int irank_link_asked=ilink_asked*nranks+irank;
      
      //if it is non local search it in the list of to-be-gathered
      all_to_all_gathering_list_t::iterator it=gat.find(irank_link_asked);
      
      //if it is already in the list, return its position
      if(it!=gat.end()) return it->second;
      else
        {
          //otherwise add it to the list of to-be-gathered
          int nel_gathered=gat.size();
          int igathered=4*loc_vol+nel_gathered;
          gat[irank_link_asked]=igathered;
          
          return igathered;
        }
    }
}

//add all links needed for a certain site
void add_tlSym_paths(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu)
{
  int *c=glb_coord_of_loclx[ivol];      
  coords A={c[0],c[1],c[2],c[3]};                         //       P---O---N    
  coords B,C,D,E,F,G,H,I,J,K,L,M,N,O,P;      		  //       |   |   |    
  //find coord mu					  //   H---G---F---E---D
  K[mu]=L[mu]=M[mu]=(A[mu]-1+glb_size[mu])%glb_size[mu];  //   |   |   |   |   |
  I[mu]=J[mu]=B[mu]=C[mu]=A[mu];			  //   I---J---A---B---C
  D[mu]=E[mu]=F[mu]=G[mu]=H[mu]=(A[mu]+1)%glb_size[mu];	  //       |   |   |    
  N[mu]=O[mu]=P[mu]=(A[mu]+2)%glb_size[mu];		  //       K---L---M    
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
      *(ilink_to_be_used++)=add_link(gat,J,nu);
      *(ilink_to_be_used++)=add_link(gat,J,mu);
      *(ilink_to_be_used++)=add_link(gat,G,nu);
      //forward square staple
      *(ilink_to_be_used++)=add_link(gat,A,nu);
      *(ilink_to_be_used++)=add_link(gat,B,mu);
      *(ilink_to_be_used++)=add_link(gat,F,nu);
      //backward dw rectangle
      *(ilink_to_be_used++)=add_link(gat,L,mu);
      *(ilink_to_be_used++)=add_link(gat,K,nu);
      *(ilink_to_be_used++)=add_link(gat,K,mu);
      *(ilink_to_be_used++)=add_link(gat,J,mu);
      *(ilink_to_be_used++)=add_link(gat,G,nu);
      //backward backward rectangle
      *(ilink_to_be_used++)=add_link(gat,J,nu);
      *(ilink_to_be_used++)=add_link(gat,I,nu);
      *(ilink_to_be_used++)=add_link(gat,I,mu);
      *(ilink_to_be_used++)=add_link(gat,H,nu);
      *(ilink_to_be_used++)=add_link(gat,G,nu);
      //backward up rectangle
      *(ilink_to_be_used++)=add_link(gat,J,nu);
      *(ilink_to_be_used++)=add_link(gat,J,mu);
      *(ilink_to_be_used++)=add_link(gat,G,mu);
      *(ilink_to_be_used++)=add_link(gat,P,nu);
      *(ilink_to_be_used++)=add_link(gat,F,mu);
      //forward dw rectangle
      *(ilink_to_be_used++)=add_link(gat,L,mu);
      *(ilink_to_be_used++)=add_link(gat,L,nu);
      *(ilink_to_be_used++)=add_link(gat,M,mu);
      *(ilink_to_be_used++)=add_link(gat,B,mu);
      *(ilink_to_be_used++)=add_link(gat,F,nu);
      //forward forward rectangle
      *(ilink_to_be_used++)=add_link(gat,A,nu);
      *(ilink_to_be_used++)=add_link(gat,B,nu);
      *(ilink_to_be_used++)=add_link(gat,C,mu);
      *(ilink_to_be_used++)=add_link(gat,E,nu);
      *(ilink_to_be_used++)=add_link(gat,F,nu);
      //forward up rectangle
      *(ilink_to_be_used++)=add_link(gat,A,nu);
      *(ilink_to_be_used++)=add_link(gat,B,mu);
      *(ilink_to_be_used++)=add_link(gat,E,mu);
      *(ilink_to_be_used++)=add_link(gat,O,nu);
      *(ilink_to_be_used++)=add_link(gat,F,mu);
    }
}

//add all the links needed to compute staple separataly for each box
THREADABLE_FUNCTION_1ARG(add_links_to_paths, all_to_all_gathering_list_t**,gl)
{
  GET_THREAD_ID();
  
  //add the links to paths
  NISSA_PARALLEL_LOOP(ibox,0,16)
    {
      //find base for curr box
      int ibase=0;
      for(int jbox=0;jbox<ibox;jbox++) ibase+=nsite_per_box[jbox];
      ibase*=4;

      //scan all the elements of sub-box, selecting only those with the good parity
      for(int dir=0;dir<4;dir++)
	for(int par=0;par<4;par++)
	  {	  
	    for(int ibox_dir_par=ibase;ibox_dir_par<ibase+nsite_per_box_dir_par[par+4*(dir+4*ibox)];ibox_dir_par++)
	      {
		int ivol=ivol_of_box_dir_par[ibox_dir_par];
		add_tlSym_paths(ilink_per_paths+ibox_dir_par*nlinks_per_paths_site,*(gl[ibox]),ivol,dir);
	      }
	    ibase+=nsite_per_box_dir_par[par+4*(dir+4*ibox)];
	  }
    }
  THREAD_BARRIER();
}}

//initialize the geometry of the sub-boxes
void init_box_geometry()
{
  ivol_of_box_dir_par=nissa_malloc("ivol_of_box_dir_par",4*loc_vol,int);
  ilink_per_paths=nissa_malloc("ilink_per_paths",nlinks_per_paths_site*4*loc_vol,int);
  
  //get the size of box 0
  for(int mu=0;mu<4;mu++)
    {
      if(loc_size[mu]<4) crash("loc_size[%d]=%d must be at least 4",mu,loc_size[mu]);
      box_size[0][mu]=loc_size[mu]/2;
    }

  //get coords of cube ans box size
  coords box_coord[16];
  coords nboxes={2,2,2,2};
  for(int ibox=0;ibox<16;ibox++)
    {
      //coords
      verbosity_lv2_master_printf("Box %d coord [ ",ibox);
      coord_of_lx(box_coord[ibox],ibox,nboxes);
      for(int mu=0;mu<4;mu++) verbosity_lv2_master_printf("%d ",box_coord[ibox][mu]);
      
      //size
      verbosity_lv2_master_printf("] size [ ",ibox);
      nsite_per_box[ibox]=1;
      for(int mu=0;mu<4;mu++)
	{
	  if(ibox!=0) box_size[ibox][mu]=((box_coord[ibox][mu]==0)?(box_size[0][mu]):(loc_size[mu]-box_size[0][mu]));
	  nsite_per_box[ibox]*=box_size[ibox][mu];
	  verbosity_lv2_master_printf("%d ",box_size[ibox][mu]);
	}
      verbosity_lv2_master_printf("], nsites: %d\n",nsite_per_box[ibox]);
    }

  //find the elements of all boxes
  int ibox_dir_par=0; //order in the parity-split order
  for(int ibox=0;ibox<16;ibox++)
    for(int dir=0;dir<4;dir++)
      for(int par=0;par<4;par++)
	{
	  nsite_per_box_dir_par[par+4*(dir+4*ibox)]=0;
	  for(int isub=0;isub<nsite_per_box[ibox];isub++)
	    {
	      //get coordinates of site
	      coords isub_coord;
	      coord_of_lx(isub_coord,isub,box_size[ibox]);
	      
	      //get coords in the local size, and parity
	      coords ivol_coord;
	      int site_par=0;
	      for(int mu=0;mu<4;mu++)
		{
		  ivol_coord[mu]=box_size[0][mu]*box_coord[ibox][mu]+isub_coord[mu];
		  site_par+=((mu==dir)?2:1)*ivol_coord[mu];
		}
	      site_par=site_par%4;
	      
	      //map sites in current parity
	      if(site_par==par)
		{
		  ivol_of_box_dir_par[ibox_dir_par++]=lx_of_coord(ivol_coord,loc_size);
		  nsite_per_box_dir_par[par+4*(dir+4*ibox)]++;
		}
	    }
	}  
  
  //init the links
  all_to_all_gathering_list_t *gl[16];
  for(int ibox=0;ibox<16;ibox++) gl[ibox]=new all_to_all_gathering_list_t;
  add_links_to_paths(gl);
  
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
  verbosity_lv2_master_printf("Max cached links: %d\n",max_cached_link);
  if(max_cached_link>bord_vol+edge_vol) crash("larger buffer needed");
  
  base_init_time+=take_time();
}

//initialize the simulation
void init_simulation(char *path)
{
  base_init_time-=take_time();
  
  //////////////////////////// read the input /////////////////////////
  
  //open input file
  open_input(path);
  
  //init the grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);  
  
  read_str_str("GaugeObsPath",gauge_obs_path,1024); //gauge observables path
  read_str_int("MaxNConfs",&max_nconfs); //number of confs to produce
  read_str_int("Seed",&seed); //seed

  //kind of action and evolution pars
  char gauge_action_name_str[1024];
  read_str_str("GaugeAction",gauge_action_name_str,1024);
  if(strcmp(gauge_action_name_str,"Wilson")==0) gauge_action_name=Wilson_action;
  else
    if(strcmp(gauge_action_name_str,"tlSym")==0) gauge_action_name=tlSym_action;
    else crash("unknown gauge action: %s",gauge_action_name_str);
  read_str_double("Beta",&beta);
  read_pure_gauge_evol_pars(evol_pars);
  
  //read in and out conf path
  read_str_str("ConfPath",conf_path,1024);
  read_str_str("StoreConfPath",store_conf_path,1024);
  read_str_int("StoreConfEach",&store_conf_each);
  read_str_int("StoreRunningTempConf",&store_running_temp_conf);
  
  //read if configuration must be generated cold or hot
  char start_conf_cond_str[1024];
  read_str_str("StartConfCond",start_conf_cond_str,1024);
  start_conf_cond_t start_conf_cond=UNSPEC_COND;
  if(strcasecmp(start_conf_cond_str,"HOT")==0) start_conf_cond=HOT;
  if(strcasecmp(start_conf_cond_str,"COLD")==0) start_conf_cond=COLD;
  if(start_conf_cond==UNSPEC_COND)
    crash("unknown starting condition cond %s, expected 'HOT' or 'COLD'",start_conf_cond_str);
  
  close_input();
  
  ////////////////////////// allocate stuff ////////////////////////
   
  //allocate conf and staples
  conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);

  //load conf or generate it
  if(file_exists(conf_path))
    {
      master_printf("File %s found, loading\n",conf_path);
      read_conf();
    }
  else
    {
      //start the random generator using passed seed
      start_loc_rnd_gen(seed);
      
      //generate hot or cold conf
      if(start_conf_cond==HOT)
	{
	  master_printf("File %s not found, generating hot conf\n",conf_path);
	  generate_hot_lx_conf(conf);
	}
      else
	{
	  master_printf("File %s not found, generating cold conf\n",conf_path);
	  generate_cold_lx_conf(conf);
	}
      
      //reset conf id
      iconf=0;
      
      //write initial measures
      measure_gauge_obs(true);
    }  
  
  init_box_geometry();
}

//finalize everything
void close_simulation()
{
  master_printf("========== Performance report ===========\n");
  master_printf("Basic initialization time: %lg sec\n",base_init_time);
  master_printf("Communication time: %lg sec\n",comm_time);
  master_printf("Link update time: %lg sec\n",comp_time);
  master_printf("Measurement time: %lg sec\n",meas_time);
  master_printf("Read conf time: %lg sec\n",read_time);
  master_printf("Write conf time: %lg sec\n",write_time);
  master_printf("=========================================\n");
  master_printf("\n");
  
  if(!store_running_temp_conf) write_conf();
  nissa_free(conf);
  nissa_free(buf_out);
  nissa_free(buf_in);
  
  nissa_free(ivol_of_box_dir_par);
  nissa_free(ilink_per_paths);
  for(int ibox=0;ibox<16;ibox++) delete box_comm[ibox];
}

void check()
{  
  //////////////////////////// make sure everything is hit ///////////////////////
  
  {
    coords *hit=nissa_malloc("hit",loc_vol,coords);
    vector_reset(hit);
    
    int ibase=0;
    for(int ibox=0;ibox<16;ibox++)
      for(int dir=0;dir<4;dir++)
	for(int par=0;par<4;par++)
	  {
	    for(int ibox_dir_par=0;ibox_dir_par<nsite_per_box_dir_par[par+4*(dir+4*ibox)];ibox_dir_par++)
	      {
		int ivol=ivol_of_box_dir_par[ibox_dir_par+ibase];
		if(ivol>=loc_vol) crash("ibox %d ibox_par %d ibase %d ivol %d",ibox,ibox_dir_par,ibase,ivol);
		hit[ivol][dir]++;
	      }
	    ibase+=nsite_per_box_dir_par[par+4*(dir+4*ibox)];
	  }
    
    for(int ivol=0;ivol<loc_vol;ivol++)
      for(int mu=0;mu<4;mu++)
	if(hit[ivol][mu]!=1) crash("missing hit ivol %d mu %d",ivol,mu);
    
    nissa_free(hit);
  }
  
  //////////////////////////// make sure everything is hit in the appropriate order ///////////////////////

  {
    int *hit=nissa_malloc("hit",4*loc_vol+max_cached_link,int);
    int ibase=0;
    for(int ibox=0;ibox<16;ibox++)
      for(int dir=0;dir<4;dir++)
	for(int par=0;par<4;par++)
	  {
	    vector_reset(hit);
	    for(int iter=0;iter<2;iter++)
	      //at first iter mark all site updating
	      //at second iter check that none of them is internally used
	      for(int ibox_dir_par=0;ibox_dir_par<nsite_per_box_dir_par[par+4*(dir+4*ibox)];ibox_dir_par++)
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
		    for(int ihit=0;ihit<nlinks_per_paths_site;ihit++)
		      {
			int link=ilink_per_paths[ihit+nlinks_per_paths_site*(ibox_dir_par+ibase)];
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
            ibase+=nsite_per_box_dir_par[par+4*(dir+4*ibox)];
	  }
    
    nissa_free(hit);
  }
}

//compute the summ of the staples pointed by "ilinks"
void compute_tlSym_staples(su3 staples,su3 *links,int *ilinks)
{
  su3 squares,rectangles;
  su3_put_to_zero(squares);
  su3_put_to_zero(rectangles);
  
  for(int inu=0;inu<3;inu++)
    {  
      su3 temp1,temp2;
      //backward square staple
      unsafe_su3_dag_prod_su3(temp1,links[ilinks[ 0]],links[ilinks[ 1]]);
      su3_summ_the_prod_su3(squares,temp1,links[ilinks[ 2]]);
      //forward square staple
      unsafe_su3_prod_su3(temp1,links[ilinks[ 3]],links[ilinks[ 4]]);
      su3_summ_the_prod_su3_dag(squares,temp1,links[ilinks[ 5]]);
      //backward dw rectangle
      unsafe_su3_dag_prod_su3_dag(temp1,links[ilinks[ 6]],links[ilinks[ 7]]);
      unsafe_su3_prod_su3(temp2,temp1,links[ilinks[ 8]]);
      unsafe_su3_prod_su3(temp1,temp2,links[ilinks[ 9]]);
      su3_summ_the_prod_su3(rectangles,temp1,links[ilinks[10]]);
      //backward backward rectangle
      unsafe_su3_dag_prod_su3_dag(temp1,links[ilinks[11]],links[ilinks[12]]);
      unsafe_su3_prod_su3(temp2,temp1,links[ilinks[13]]);
      unsafe_su3_prod_su3(temp1,temp2,links[ilinks[14]]);
      su3_summ_the_prod_su3(rectangles,temp1,links[ilinks[15]]);
      //backward up rectangle
      unsafe_su3_dag_prod_su3(temp1,links[ilinks[16]],links[ilinks[17]]);
      unsafe_su3_prod_su3(temp2,temp1,links[ilinks[18]]);
      unsafe_su3_prod_su3(temp1,temp2,links[ilinks[19]]);
      su3_summ_the_prod_su3_dag(rectangles,temp1,links[ilinks[20]]);
      //forward dw rectangle
      unsafe_su3_dag_prod_su3(temp1,links[ilinks[21]],links[ilinks[22]]);
      unsafe_su3_prod_su3(temp2,temp1,links[ilinks[23]]);
      unsafe_su3_prod_su3(temp1,temp2,links[ilinks[24]]);
      su3_summ_the_prod_su3_dag(rectangles,temp1,links[ilinks[25]]);
      //forward forward rectangle
      unsafe_su3_prod_su3(temp1,links[ilinks[26]],links[ilinks[27]]);
      unsafe_su3_prod_su3(temp2,temp1,links[ilinks[28]]);
      unsafe_su3_prod_su3_dag(temp1,temp2,links[ilinks[29]]);
      su3_summ_the_prod_su3_dag(rectangles,temp1,links[ilinks[30]]);
      //forward up rectangle
      unsafe_su3_prod_su3(temp1,links[ilinks[31]],links[ilinks[32]]);
      unsafe_su3_prod_su3(temp2,temp1,links[ilinks[33]]);
      unsafe_su3_prod_su3_dag(temp1,temp2,links[ilinks[34]]);
      su3_summ_the_prod_su3_dag(rectangles,temp1,links[ilinks[35]]);

      ilinks+=36;
    }
  
  //compute the summed staples
  double b1=-1.0/12,b0=1-8*b1;
  su3_linear_comb(staples,squares,b0,rectangles,b1);
}

//compute action
double compute_tlSym_action(complex paths)
{
  //coefficient of rectangles and squares
  double b1=-1.0/12,b0=1-8*b1;
  
  //compute the total action
  global_plaquette_and_rectangles_lx_conf(paths,conf);
  return (b0*6*glb_vol*(1-paths[RE])+b1*12*glb_vol*(1-paths[IM]))*beta;
}

//updated all sites with heat-bath or overrelaxation
THREADABLE_FUNCTION_4ARG(sweep_conf, update_alg_t,update_alg, quad_su3*,conf, void*,buf_out, void*,buf_in)
{
  GET_THREAD_ID();
  
  int ibase=0;
  for(int ibox=0;ibox<16;ibox++)
    {
      //communicate needed links
      if(IS_MASTER_THREAD) comm_time-=take_time();
      box_comm[ibox]->communicate(conf,conf,sizeof(su3),buf_out,buf_in);
      if(IS_MASTER_THREAD) comm_time+=take_time();
      
      if(IS_MASTER_THREAD) comp_time-=take_time();
      for(int dir=0;dir<4;dir++)
	for(int par=0;par<4;par++)
	  {
	    //scan all the box
	    int nbox_dir_par=nsite_per_box_dir_par[par+4*(dir+4*ibox)];
	    NISSA_PARALLEL_LOOP(ibox_dir_par,ibase,ibase+nbox_dir_par)
	      {
		//compute the staples
		su3 staples;
		compute_tlSym_staples(staples,(su3*)conf,ilink_per_paths+nlinks_per_paths_site*ibox_dir_par);
		
		//find new link
		int ivol=ivol_of_box_dir_par[ibox_dir_par];
		su3 new_link;
		if(update_alg==HEAT)
		  su3_find_heatbath(new_link,conf[ivol][dir],staples,beta,evol_pars.nhb_hits,loc_rnd_gen+ivol);
		else
		  su3_find_overrelaxed(new_link,conf[ivol][dir],staples,evol_pars.nov_hits);
		
		//copy new link
		su3_copy(conf[ivol][dir],new_link);
	      }
	    
	    //increment the box-dir-par subset
	    ibase+=nbox_dir_par;
	    THREAD_BARRIER();
	  }
      if(IS_MASTER_THREAD) comp_time+=take_time();
    }
  
  set_borders_invalid(conf);
}}

//heatbath or overrelax algorithm for the quenched simulation case, Wilson action
void generate_new_conf(quad_su3 *conf,void *buf_out,void *buf_in)
{
  //number of hb sweeps
  for(int ihb_sweep=0;ihb_sweep<evol_pars.nhb_sweeps;ihb_sweep++) sweep_conf(HEAT,conf,buf_out,buf_in);
  //numer of overrelax sweeps
  for(int iov_sweep=0;iov_sweep<evol_pars.nov_sweeps;iov_sweep++) sweep_conf(OVER,conf,buf_out,buf_in);
}

//measure plaquette and polyakov loop
void measure_gauge_obs(bool conf_created=false)
{
  meas_time-=take_time();
  
  //open creating or appending
  FILE *file=open_file(gauge_obs_path,conf_created?"w":"a");

  //compute action
  double time_action=-take_time();
  double paths[2];
  double action=compute_tlSym_action(paths);
  master_printf("Action: %015.15lg measured in %lg sec\n",action,time_action+take_time());

  master_fprintf(file,"%6d\t%015.15lg\t%015.15lg\t%015.15lg\n",iconf,action,paths[0],paths[1]);
  
  if(rank==0) fclose(file);
  meas_time+=take_time();
}

//store conf when appropriate
void store_conf_if_necessary()
{
  if(store_conf_each!=0 && iconf%store_conf_each==0)
    {
      char path[1024];
      sprintf(path,"%s.%05d",store_conf_path,iconf);
      write_conf();
    }
}

void in_main(int narg,char **arg)
{
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  ///////////////////////////////////////
  
  //generate the required amount of confs
  nprod_confs=0;
  master_printf("\n");
  do
    {
      master_printf("--------Configuration %d--------\n",iconf);
      
      // 1) produce new conf
      if(max_nconfs!=0)
	{
	  double gen_time=-take_time();
	  generate_new_conf(conf,buf_out,buf_in);
	  gen_time+=take_time();
	  master_printf("Generate new conf in %lg sec\n",gen_time);
	  nprod_confs++;
	  iconf++;
	}
      
      // 2) measure
      measure_gauge_obs();
      
      // 3) increment id and write conf
      if(store_running_temp_conf) write_conf();
      
      // 4) if conf is multiple of store_conf_each copy it
      store_conf_if_necessary();
      
      // 5) spacing between output
      master_printf("\n");
    }
  while(nprod_confs<max_nconfs && !file_exists("stop") && !file_exists("restart"));
  
  /////////////////////////////////////// timings /////////////////////////////////
  
  close_simulation();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
