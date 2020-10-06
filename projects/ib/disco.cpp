#include "nissa.hpp"
#include <unistd.h>

#define PROP_TYPE colorspinspin

#include "conf.hpp"
#include "pars.hpp"
#include "prop.hpp"

using namespace nissa;

spincolor *eta,*phi,*temp;
spin1field **bubble,*temp_bubble;
spin1field **tot_bubble;
storable_vector_t<double> tadpole_val,scalop_val,pseudop_val;

const int version=1;
int nread_sources;
int silv_par=2;

int nB;
int iB(int par,int iquark)
{return iquark+nquarks*par;}

//initialize the simulation
void init_simulation(const char *path)
{
  open_input(path);
  
  read_input_preamble();
  read_photon_pars();
  read_use_photon_field();
  read_seed_start_random();
  read_noise_type();
  read_free_theory_flag();
  read_random_gauge_transform();
  read_loc_hadr_curr();
  read_nsources();
  read_ngauge_conf();
  
  //allocate
  nB=iB(silv_par-1,nquarks-1)+1;
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  eta=nissa_malloc("eta",loc_vol+bord_vol,spincolor);
  phi=nissa_malloc("phi",loc_vol+bord_vol,spincolor);
  temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  temp_bubble=nissa_malloc("temp_bubble",loc_vol+bord_vol,spin1field);
  bubble=nissa_malloc("bubble*",nB,spin1field*);
  tot_bubble=nissa_malloc("tot_bubble*",nB,spin1field*);
  for(int iB=0;iB<nB;iB++)
    {
      bubble[iB]=nissa_malloc("bubble",loc_vol+bord_vol,spin1field);
      tot_bubble[iB]=nissa_malloc("tot_bubble",loc_vol+bord_vol,spin1field);
    }
}

//read the data
void read_data()
{
  //open the file for reading
  std::string path=combine("%s/bubbles",outfolder);
  if(!file_exists(path.c_str()))
    {
      nread_sources=0;
      master_printf("No source ever computed\n");
    }
  else
    {
      ILDG_File fin=ILDG_File_open_for_read(path.c_str());
      
      //read info and nsources
      ILDG_header head;
      ILDG_File_search_record(head,fin,"Info");
      char info_data[head.data_length+1];
      ILDG_File_read_all(info_data,fin,head.data_length);
      
      //parse
      int read_version,read_par,read_nquarks;
      int nread_info=sscanf(info_data," Version %d Par %d NQuarks %d NSources %d",&read_version,&read_par,&read_nquarks,&nread_sources);
      if(nread_info!=5) crash("could not parse info");
      
      //check the info record
      if(version!=read_version) crash("read version %d while program is version %d",read_version,version);
      if(silv_par!=read_par) crash("read parity %d while input is %d",read_par,silv_par);
      if(nquarks!=read_nquarks) crash("read nquarks %d while input is %d",read_nquarks,nquarks);
      
      //load what have been done
      if(nread_sources<nsources)
	{
	  tadpole_val.read_from_ILDG_file(fin,"Tadpole");
	  pseudop_val.read_from_ILDG_file(fin,"Pseudoscalar");
	  scalop_val.read_from_ILDG_file(fin,"Scalar");
	  
	  //read
	  for(int par=0;par<silv_par;par++)
	    for(int iquark=0;iquark<nquarks;iquark++)
	      {
		ILDG_header head=ILDG_File_get_next_record_header(fin);
		std::string exp_type=combine("par_%d_mass_%d_r_%d",par,iquark);
		if(exp_type!=head.type) crash("expecting \"%s\", obtained \"%s\"",exp_type.c_str(),head.type);
		ILDG_File_read_ildg_data_all(bubble[par][iquark],fin,head);
	      }
	}
      
      ILDG_File_close(fin);
    }
}

//write the data
void write_data()
{
  //open the file for output
  ILDG_File fout=ILDG_File_open_for_write(combine("%s/bubbles",outfolder).c_str());
  
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  
  //write the info
  std::ostringstream info;
  info<<" Version "<<version
      <<" Par "<<silv_par
      <<" Nquarks "<<nquarks
      <<" NSources "<<nsources;
  ILDG_string_message_append_to_last(&mess,"Info",info.str().c_str());
  
  //write tadpole and the rest
  tadpole_val.append_to_message_with_name(mess,"Tadpole");
  pseudop_val.append_to_message_with_name(mess,"Pseudoscalar");
  scalop_val.append_to_message_with_name(mess,"Scalar");
  ILDG_File_write_all_messages(fout,&mess);
  
  //store the different parity, masses
  for(int par=0;par<silv_par;par++)
    for(int iquark=0;iquark<nquarks;iquark++)
      ILDG_File_write_ildg_data_all(fout,bubble[par][iquark],sizeof(spin1field),combine("par_%d_quark_%d_r",par,iquark).c_str());
	
  ILDG_File_close(fout);
}

//close deallocating everything
void close()
{
  //free
  nissa_free(eta);
  nissa_free(phi);
  nissa_free(temp);
  nissa_free(conf);
  for(int iB=0;iB<nB;iB++)
    {
      nissa_free(tot_bubble[iB]);
      nissa_free(bubble[iB]);
    }
  nissa_free(temp_bubble);
  nissa_free(tot_bubble);
  nissa_free(bubble);
  
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g, of which:\n",tot_prog_time);
}

//tadpole
THREADABLE_FUNCTION_3ARG(compute_tadpole, spincolor*,eta, spincolor*,phi, int,r)
{
  GET_THREAD_ID();
  
  if(twisted_run) insert_tm_tadpole(temp,conf,phi,r,tadpole,-1);
  else            insert_Wilson_tadpole(temp,conf,phi,tadpole,-1);
  
  //trace
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<NCOL;ic++)
	safe_complex_conj1_prod(temp[ivol][id][ic],eta[ivol][id][ic],temp[ivol][id][ic]);
  THREAD_BARRIER();
  
  //reduce and store
  complex tadpole_res;
  complex_vector_glb_collapse(tadpole_res,(complex*)temp,sizeof(spincolor)/sizeof(complex)*loc_vol);
  if(IS_MASTER_THREAD)
    {
      tadpole_val.push_back(tadpole_res[RE]);
      tadpole_val.push_back(tadpole_res[IM]);
    }
  
  master_printf("Tadpole: %+16.16lg %+16.16lg\n",tadpole_res[0],tadpole_res[1]);
}
THREADABLE_FUNCTION_END

//scalar and pseudoscalar
THREADABLE_FUNCTION_2ARG(compute_scalar_pseudoscalar, spincolor*,phi, spincolor*,eta)
{
  GET_THREAD_ID();
  
  //trace
  complex loc_scalar={0,0},loc_pseudo={0,0};
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      spincolor temp;
      unsafe_dirac_prod_spincolor(temp,base_gamma+5,eta[ivol]);
      spincolor_scalar_prod(loc_scalar,eta[ivol],phi[ivol]);
      spincolor_scalar_prod(loc_pseudo,eta[ivol],temp);
    }
  THREAD_BARRIER();
  
  //reduce and store
  complex glb_scalar,glb_pseudo;
  glb_reduce_complex(glb_scalar,loc_scalar);
  glb_reduce_complex(glb_pseudo,loc_pseudo);
  if(IS_MASTER_THREAD)
    for(int ri=0;ri<2;ri++)
      {
	scalop_val.push_back(glb_scalar[ri]);
	pseudop_val.push_back(glb_pseudo[ri]);
      }
  
  master_printf("Pseudo: %+16.16lg %+16.16lg\n",glb_pseudo[0],glb_pseudo[1]);
  master_printf("Scalar: %+16.16lg %+16.16lg\n",glb_scalar[0],glb_scalar[1]);
}
THREADABLE_FUNCTION_END

//matrix element of the vector current as a function of the photon position
THREADABLE_FUNCTION_5ARG(vector_matrix_element, spin1field*,out, spincolor*,sink, quad_su3*,conf, int,r, spincolor*,source)
{
  GET_THREAD_ID();
  
  //set parameters
  dirac_matr GAMMA;
  if(twisted_run) dirac_prod_idouble(&GAMMA,base_gamma+5,-tau3[r]);
  else            GAMMA=base_gamma[0];
  
  //reset the output and communicate borders
  vector_reset(out);
  communicate_lx_spincolor_borders(sink);
  communicate_lx_spincolor_borders(source);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<NDIM;mu++)
      {
	//find neighbors
	int ifw=loclx_neighup[ivol][mu];
	
	//transport down and up
	spincolor fw,bw;
	unsafe_su3_prod_spincolor(fw,conf[ivol][mu],source[ifw]);
	unsafe_su3_dag_prod_spincolor(bw,conf[ivol][mu],source[ivol]);
	
	//put GAMMA-gmu on the forward piece
	spincolor GA_fw,gmu_fw,GA_M_gmu_fw;
	unsafe_dirac_prod_spincolor(GA_fw,&GAMMA,fw);
	unsafe_dirac_prod_spincolor(gmu_fw,base_gamma+map_mu[mu],fw);
	spincolor_subt(GA_M_gmu_fw,GA_fw,gmu_fw);
	//put GAMMA+gmu on the backward piece
	spincolor GA_bw,gmu_bw,GA_P_gmu_bw;
	unsafe_dirac_prod_spincolor(GA_bw,&GAMMA,bw);
	unsafe_dirac_prod_spincolor(gmu_bw,base_gamma+map_mu[mu],bw);
	spincolor_summ(GA_P_gmu_bw,GA_bw,gmu_bw);
	
	//take the scalar product
	complex c_fw,c_bw;
	spincolor_scalar_prod(c_fw,sink[ivol],GA_M_gmu_fw);
	spincolor_scalar_prod(c_bw,sink[ifw],GA_P_gmu_bw);
	
	//summ the result
	complex fact_fw={0,-0.5},fact_bw={0,+0.5};
	complex_summ_the_prod(out[ivol][mu],c_fw,fact_fw);
	complex_summ_the_prod(out[ivol][mu],c_bw,fact_bw);
      }
  
  set_borders_invalid(out);
}
THREADABLE_FUNCTION_END

//compute everything with a single mass
void single_mass(int par,int iquark)
{
  master_printf(" iquark %d/%d\n",iquark+1,nquarks);
  
  //solve for phi for each quark
  get_qprop(phi,eta,iquark);
  
  //compute the various operator
  compute_tadpole(eta,phi,qr[iquark]);
  compute_scalar_pseudoscalar(eta,phi);
  
  //compute the single bubble
  int ib=iB(par,iquark);
  vector_matrix_element(bubble[ib],eta,conf,qr[iquark],phi);
  
  //summ to the total
  double_vector_summassign((double*)(tot_bubble[ib]),(double*)(bubble[ib]),sizeof(spin1field)/sizeof(double)*loc_vol);
}

//compute everything with a single source
void single_source(int isource)
{
  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
  
  //generate the source
  generate_undiluted_source(eta,noise_type,-1);
  
  for(int iquark=0;iquark<nquarks;iquark++)
    single_mass(isource%silv_par,iquark);
}

//read the number of sources in the file and check that it works
bool sources_missing()
{
  read_data();
  if(nread_sources<nsources)
    {
      master_printf("Computed %d sources, we need to do %d\n",nread_sources,nsources);
      return true;
    }
  else
    {
      master_printf("Computed all %d sources, skipping\n",nsources);
      return false;
    }
}

void finish_conf()
{
  write_data();
  nanalyzed_conf++;
  unlink(combine("%s/running",outfolder).c_str());
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  //loop over the configs
  int iconf=0,enough_time=1;
  while(iconf<ngauge_conf && enough_time && !file_exists(stop_path) && read_conf_parameters(iconf,sources_missing))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      //loop over sources
      for(int isource=nread_sources;isource<nsources;isource++)
	single_source(isource);
      
      finish_conf();
    }
  
  //close the simulation
  tot_prog_time+=take_time();
  close();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
