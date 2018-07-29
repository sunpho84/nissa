#include "nissa.hpp"

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

using namespace nissa;

int hits_done_so_far;

/////////////////////////////////////////////////////////////////

namespace curr
{
  //fileds J,j and xi
  spin1field **fields;
  
  //fields to store the average, to get stoch estimate of the current, or to convolve it
  enum field_t{J,j,xi};
  int nfields_type=3;
  
  //get the index of the field
  int ifield_idx(int iquark,int f)
  {
    return f+nfields_type*iquark;
  }
  
  //allocate all currents
  int ntot_fields;
  void allocate_all_fields()
  {
    ntot_fields=ifield_idx(qprop_name_list.size()-1,nfields_type-1)+1;
    
    fields=nissa_malloc("fields",ntot_fields,spin1field*);
    for(int ifield=0;ifield<ntot_fields;ifield++) fields[ifield]=nissa_malloc("field",loc_vol+bord_vol, spin1field);
  }
  
  //free all currents
  void free_all_fields()
  {
    for(int ifield=0;ifield<ntot_fields;ifield++) nissa_free(fields[ifield]);
    nissa_free(fields);
  }
  
  //form the path of the current
  std::string path(int iquark)
  {
    return combine("%s/%s_current_%d_hits",outfolder,qprop_name_list[iquark].c_str(),iquark);
  }
  
  //read all currents
  void load_all_or_reset()
  {
    for(int iquark=0;iquark<(int)qprop_name_list.size();iquark++)
      {
	spin1field *c=fields[ifield_idx(iquark,J)];
	if(hits_done_so_far) read_real_vector(c,path(hits_done_so_far),"ildg-binary-data");
	else vector_reset(c);
      }
  }
  
  //write all currents
  void store_all()
  {
    for(int iquark=0;iquark<(int)qprop_name_list.size();iquark++)
      write_real_vector(path(hits_done_so_far),fields[ifield_idx(iquark,J)],64,"ildg-binary-data");
  }
}

/////////////////////////////////////////////////////////////////

//path of the number of hits
std::string hits_done_so_far_path()
{
  return combine("%s/hits_done_so_far",outfolder);
}

//gets the number of hits done
void get_hits_done_so_far()
{
  hits_done_so_far=
    file_exists(hits_done_so_far_path())
    ?
    master_fscan_int(hits_done_so_far_path())
    :
    0;
}

//write the number of hits done so far
void write_hits_done_so_far()
{
  FILE *fout=open_file(hits_done_so_far_path(),"w");
  master_fprintf(fout,"%d\n",hits_done_so_far);
  close_file(fout);
}

//skip the first hits
void skip_hits_done_so_far()
{
  get_hits_done_so_far();
  skip_nhits(0,hits_done_so_far);
}

//give a name to the propagator
std::string get_prop_name(int iquark,insertion_t ins)
{
  return combine("Q%d_%s",iquark,ins_name[ins]);
}

//compute the matrix element of the conserved current between two propagators
THREADABLE_FUNCTION_4ARG(calc_cur, spin1field*,cur, spincolor*,source, quad_su3*,conf, spincolor*,prop)
{
  GET_THREAD_ID();
  
  vector_reset(cur);
  
  communicate_lx_spincolor_borders(source);
  communicate_lx_spincolor_borders(prop);
  communicate_lx_quad_su3_borders(conf);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<NDIM;mu++)
      {
	// int ifw=loclx_neighup[ivol][mu];
	// spincolor f;
	
	// //piece psi_ivol U_ivol psi_fw
	// unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ifw]);
	// spincolor_scalar_prod(cur[ivol][mu],source[ivol],f);
	
	int ifw=loclx_neighup[ivol][mu];
	spincolor f,Gf;
	complex c;
	
	//piece psi_ivol U_ivol psi_fw
	unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ifw]);
	spincolor_copy(Gf,f);
	dirac_subt_the_prod_spincolor(Gf,base_gamma+igamma_of_mu[mu],f);
	spincolor_scalar_prod(c,source[ivol],Gf);
	complex_summ_the_prod_idouble(cur[ivol][mu],c,-0.5);
	
	//piece psi_fw U_ivol^dag psi_ivol
	unsafe_su3_dag_prod_spincolor(f,conf[ivol][mu],prop[ivol]);
	spincolor_copy(Gf,f);
	dirac_summ_the_prod_spincolor(Gf,base_gamma+igamma_of_mu[mu],f);
	spincolor_scalar_prod(c,source[ifw],Gf);
	complex_summ_the_prod_idouble(cur[ivol][mu],c,+0.5);
      }
  set_borders_invalid(cur);
}
THREADABLE_FUNCTION_END

void in_main_Marcus(int narg,char **arg)
{
  init_grid(16,8);
  
  spincolor *source=nissa_malloc("source",loc_vol,spincolor);
  spincolor *prop=nissa_malloc("prop",loc_vol,spincolor);
  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  
  read_ildg_gauge_conf(conf,"conf");
  momentum_t th={1,0,0,0};
  put_boundaries_conditions(conf,th,0,0);
  
  read_real_vector(source,"source.0000.00000","scidac-binary-data");
  read_real_vector(prop,"source.0000.00000.inverted","scidac-binary-data");
  
  spin1field *cur=nissa_malloc("cur",loc_vol,spin1field);
  calc_cur(cur,source,conf,prop);
  
  FILE *fout=open_file("cur.txt","w");
  master_fprintf(fout,"# [loops_em]  norm square of source     = %.16lg\n",double_vector_glb_norm2(source,loc_vol));
  master_fprintf(fout,"# [loops_em]  norm square of propagator = %.16lg\n",double_vector_glb_norm2(prop,loc_vol));
  NISSA_LOC_VOL_LOOP(ivol)
  {
    master_fprintf(fout,"# [loops_em] x\t");
    for(int mu=0;mu<NDIM;mu++) master_fprintf(fout,"%d ",loc_coord_of_loclx[ivol][mu]);
    master_fprintf(fout,"\n");
    
    for(int mu=0;mu<NDIM;mu++) master_fprintf(fout,"%d" "\t" "%.16lg" "\t" "%.16lg\n",mu,cur[ivol][mu][RE],cur[ivol][mu][IM]);
  }
  close_file(fout);
  
  nissa_free(cur);
  nissa_free(conf);
  nissa_free(prop);
  nissa_free(source);
}

void init_simulation(int narg,char **arg)
{
  //set the prefix for contr file
  mes2pts_prefix="disco";
  
  std::string input_path="";
  
  //parse opts
  int c;
  while((c=getopt(narg,arg,"i:"))!= -1)
    switch (c)
      {
      case 'i': input_path=optarg; break;
      default: crash("Unknown option -%c",optopt);
      }
  
  if(input_path=="") crash("Please specify -i");
  open_input(input_path);
  
  //init the grid
  read_init_grid();
  
  //Wall time
  read_str_double("WallTime",&wall_time);
  
  //Seed
  read_seed_start_random();
  
  //Set stochastic source
  stoch_source=true;
  
  //Read the noise type
  char str_noise_type[20];
  read_str_str("NoiseType",str_noise_type,20);
  rnd_t noise_type=convert_str_to_rnd_t(str_noise_type);
  
  //Everything diluted so far
  set_diluted_spin(1);
  set_diluted_color(1);
  set_diluted_space(1);
  
  //Number of hits
  read_nhits();
  
  //put the source in the list
  std::string source_name="source";
  int store_source=false;
  ori_source_name_list.push_back(source_name);
  Q[source_name].init_as_source(noise_type,ALL_TIMES,0,store_source);
  
  //Read whether this is a twisted and/or clover run
  read_twisted_run();
  read_clover_run();
  
  /////////////////////////////////////////////////////////////////
  
  //Boundary conditions
  momentum_t theta;
  theta[0]=temporal_bc;
  for(int mu=1;mu<NDIM;mu++) theta[mu]=0;
  int store_prop=false;
  
  //Read kappa
  double kappa=0.125;
  if(twisted_run) read_str_double("Kappa",&kappa);
  
  //Read the number of quarks
  int nquarks;
  read_str_int("NQuarks",&nquarks);
  
  //Placeholder for mass, kappa and resiude
  if(twisted_run) expect_str("MassRResidue");
  else            expect_str("KappaResidue");
  
  for(int iquark=0;iquark<nquarks;iquark++)
    {
      //Read mass and r
      double mass=0.0;
      int r=0;
      if(twisted_run)
	{
	  read_double(&mass);
	  master_printf("Read variable 'Mass' with value: %lg\n",mass);
	  read_int(&r);
	  master_printf("Read variable 'R' with value: %d\n",r);
	  
	  //include tau in the mass
	  mass*=tau3[r];
	}
      else
	read_double(&kappa);
      
      //Residue for solver
      double residue;
      read_double(&residue);
      master_printf("Read variable 'Residue' with value: %lg\n",residue);
      
      //Set charge to zero
      double charge=0.0;
      
      //Put the propagator in the list
      std::string prop_name=get_prop_name(iquark,PROP);
      qprop_name_list.push_back(prop_name);
      Q[prop_name].init_as_propagator(PROP,{{source_name,{1.0,0.0}}},ALL_TIMES,residue,kappa,mass,r,charge,theta,store_prop);
      
      //Add insertion of S, P and T
      for(auto ins : {SCALAR,PSEUDO,TADPOLE})
	{
	  std::string prop_name_ins=get_prop_name(iquark,ins);
	  qprop_name_list.push_back(prop_name_ins);
	  Q[prop_name_ins].init_as_propagator(ins,{{prop_name,{1.0,0.0}}},ALL_TIMES,residue,kappa,mass,r,charge,theta,store_prop);
	  
	  mes2pts_contr_map.push_back(mes_contr_map_t(prop_name_ins,source_name,prop_name_ins));
	}
    }
  
  //set g5 on the source if not diluted
  int id_so=(nso_spi==NDIRAC ? 0 : 5);
  mes_gamma_list.push_back(idirac_pair_t(id_so,0));
  
  read_photon_pars();
  
  read_free_theory_flag();
  read_random_gauge_transform();
  read_Landau_gauge_fix();
  read_store_conf();
  
  read_loc_hadr_curr();
  
  read_ngauge_conf();
  
  /////////////////////////////////////////////////////////////////
  
  allocate_loop_source();
  curr::allocate_all_fields();
  allocate_mes2pts_contr();
  glb_conf=nissa_malloc("glb_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  inner_conf=nissa_malloc("inner_conf",loc_vol+bord_vol+edge_vol,quad_su3);
}

void close()
{
  close_input();
  
  print_statistics();
  
  Q.clear();

  curr::free_all_fields();
  free_loop_source();
  free_mes2pts_contr();
  nissa_free(glb_conf);
  nissa_free(inner_conf);
  
  if(clover_run)
    {
      nissa_free(Cl);
      nissa_free(invCl);
    }
}

void compute_disco_PST(int ihit)
{
  //compute "2pts"
  vector_reset(mes2pts_contr);
  compute_mes2pts_contr();
  
  //print correlators
  int force_append(ihit>0);
  int skip_inner_header=true;
  std::string hit_header=combine("\n # hit %d\n\n",ihit);
  print_mes2pts_contr(1.0,force_append,skip_inner_header,hit_header);
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //init simulation according to input file
  init_simulation(narg,arg);
  
  //loop over the configs
  int iconf=0;
  while(read_conf_parameters(iconf,finish_file_present))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      //start at hits_done_so_far
      skip_hits_done_so_far();
      curr::load_all_or_reset();
      
      while(hits_done_so_far<nhits)
	{
	  start_hit(hits_done_so_far);
	  generate_propagators(hits_done_so_far);
	  compute_disco_PST(hits_done_so_far);
	  
	  hits_done_so_far++;
	}
      
      curr::store_all();
      write_hits_done_so_far();
      mark_finished();
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
