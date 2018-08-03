#include "nissa.hpp"

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

using namespace nissa;

int nquarks;
int nhits_done_so_far;
std::string source_name="source";
int compute_EU6_alt;

/////////////////////////////////////////////////////////////////

namespace curr
{
  //fields J,j and xi
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
    ntot_fields=ifield_idx(nquarks-1,nfields_type-1)+1;
    
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
  std::string path(int iquark,int nhits)
  {
    return combine("%s/current_%d_hits_quark_%d",outfolder,nhits,iquark);
  }
  
  //read all currents
  void load_all_or_reset()
  {
    for(int iquark=0;iquark<nquarks;iquark++)
      {
	spin1field *c=fields[ifield_idx(iquark,J)];
	if(nhits_done_so_far) read_real_vector(c,path(iquark,nhits_done_so_far),"ildg-binary-data");
	else vector_reset(c);
      }
  }
  
  //write all currents
  void store_all()
  {
    for(int iquark=0;iquark<nquarks;iquark++)
      write_real_vector(path(iquark,nhits_done_so_far),fields[ifield_idx(iquark,J)],64,"ildg-binary-data");
  }
}

/////////////////////////////////////////////////////////////////

namespace hits
{
  complex *field;
  complex *f1_f2;
  int nsingle_or_all=2;
  enum{SINGLE,ALL};
  
  //index to access the correct combo
  int idx(int f1,int f2,int ihit,int single_or_all)
  {
    return f1+nquarks*(f2+nquarks*(ihit+nhits*single_or_all));
  }
  
  //allocate all fields
  void allocate()
  {
    field=nissa_malloc("Eff",loc_vol,complex);
    f1_f2=nissa_malloc("f1_f2",idx(nquarks-1,nquarks-1,nhits-1,nsingle_or_all-1)+1,complex);
  }
  
  //free all fields
  void free()
  {
    nissa_free(field);
    nissa_free(f1_f2);
  }
  
  //forms the path to all hits
  std::string path()
  {
    return combine("%s/f1_f2",outfolder);
  }
  
  //store or load the hits
  enum{STORE,LOAD};
  void store_or_load(int nh,int store_or_load)
  {
    FILE *file=open_file(path(),store_or_load==STORE?"w":"r");
    for(int f1=0;f1<nquarks;f1++)
      for(int f2=0;f2<nquarks;f2++)
	for(int ihit=0;ihit<nh;ihit++)
	  for(int single_or_all=0;single_or_all<nsingle_or_all;single_or_all++)
	    {
	      complex &c=f1_f2[idx(f1,f2,ihit,single_or_all)];
	      for(int RI=0;RI<2;RI++)
		if(store_or_load==STORE) master_fprintf(file,"%.16lg%s",c[RI],RI?"\n":" ");
		else  		         c[RI]=master_fscan_double(file);
	    }
    
    close_file(file);
  }
  
  //store the hits
  void store(int nh)
  {
    store_or_load(nh,STORE);
  }
  
  //load the hits if non zero
  void load(int nh)
  {
    if(nh)
      store_or_load(nh,LOAD);
  }
}

/////////////////////////////////////////////////////////////////

//Compute diagrams 1,2,4,and 6alt
THREADABLE_FUNCTION_0ARG(compute_EU1_EU2_EU4_EU6alt)
{
  //suffix and complex weight for EU1,2,4,6alt
  std::vector<std::tuple<std::string,double,double,std::string>> suff_weight_id={{"PSEUDO",0.0,1.0,"1"},{"SCALAR",1.0,0.0,"2"},{"TADPOLE",1.0,0.0,"4"}};
  if(compute_EU6_alt) suff_weight_id.push_back({"PROP_F_PROP_F",1.0,0.0,"6alt"});
  
  for(auto &i : suff_weight_id)
    {
      //decompose pars
      std::string suff=std::get<0>(i);
      complex w={std::get<1>(i),std::get<2>(i)};
      std::string id=std::get<3>(i);
      
      //open the file
      FILE *fout=open_file(combine("%s/EU%s",outfolder,id.c_str()),nhits_done_so_far?"a":"w");
      
      for(int iquark=0;iquark<nquarks;iquark++)
	{
	  //takes the scalar product
	  complex p;
	  compute_prop_scalprod(p,source_name,combine("Q%d_%s",iquark,suff.c_str()));
	  
	  //put the weight and write
	  complex_prodassign(p,w);
	  master_fprintf(fout,"%.16lg %.16lg\t",p[RE],p[IM]);
	}
      
      master_fprintf(fout,"\n");
      close_file(fout);
    }
}
THREADABLE_FUNCTION_END

//Compute EU5 and EU6
void compute_EU5_EU6()
{
  //tag to access files and data
  enum{iEU5,iEU6};
  int eu_tag[2]={5,6};
  
  //open the file
  FILE *fout_EU[2];
  for(int ieu=0;ieu<2;ieu++)
    fout_EU[ieu]=open_file(combine("%s/EU%d",outfolder,eu_tag[ieu]),"w");
  
  //store the summ
  complex E_f1_f2[nquarks*nquarks];
  memset(E_f1_f2,0,sizeof(nquarks*nquarks));
  
  for(int ih=0;ih<nhits_done_so_far;ih++)
    {
      for(int iquark=0;iquark<nquarks;iquark++)
	for(int jquark=0;jquark<nquarks;jquark++)
	  {
	    //update the sum of single hit, E, up to ih
	    int i=jquark+nquarks*iquark;
	    complex_summassign(E_f1_f2[i],hits::f1_f2[hits::idx(iquark,jquark,ih,hits::SINGLE)]);
	    
	    //compute the difference of J*Xi with E
	    complex H_f1_f2;
	    complex &JXi_f1_f2=hits::f1_f2[hits::idx(iquark,jquark,ih,hits::ALL)];
	    complex_subt(H_f1_f2,JXi_f1_f2,E_f1_f2[i]);
	    
	    //normalize: aovoid dividing by 0, replacing 0 with 1
	    int norm=(ih+1)*std::max(1,ih);
	    complex_prodassign_double(H_f1_f2,1.0/norm);
	    
	    //print
	    master_fprintf(fout_EU[iEU5],"%.16lg %.16lg\t",H_f1_f2[RE],H_f1_f2[IM]);
	    
	    //compute EU6
	    if(iquark==jquark)
	      {
		complex F_f;
		complex_prod_double(F_f,E_f1_f2[jquark+nquarks*iquark],1.0/(ih+1));
		complex_subtassign(F_f,H_f1_f2);
		
		//include obscure nomralization
		complex_prodassign_double(F_f,glb_vol/6.0);
		
		//print
		master_fprintf(fout_EU[iEU6],"%.16lg %.16lg\t",F_f[RE],F_f[IM]);
	      }
	  }
      
      //send to new line
      for(int ieu=0;ieu<2;ieu++)
	master_fprintf(fout_EU[ieu],"\n");
    }
  
  //close both
  for(int ieu=0;ieu<2;ieu++)
    close_file(fout_EU[ieu]);
}

/////////////////////////////////////////////////////////////////

//path of the number of hits
std::string nhits_done_so_far_path()
{
  return combine("%s/nhits_done_so_far",outfolder);
}

//gets the number of hits done
void get_nhits_done_so_far()
{
  nhits_done_so_far=
    file_exists(nhits_done_so_far_path())
    ?
    master_fscan_int(nhits_done_so_far_path())
    :
    0;
}

//write the number of hits done so far
void write_nhits_done_so_far()
{
  FILE *fout=open_file(nhits_done_so_far_path(),"w");
  master_fprintf(fout,"%d\n",nhits_done_so_far);
  close_file(fout);
}

//skip the first hits
void skip_hits_done_so_far()
{
  get_nhits_done_so_far();
  skip_nhits(0,nhits_done_so_far);
}

//give a name to the propagator
std::string get_prop_name(int iquark,insertion_t ins)
{
  return combine("Q%d_%s",iquark,ins_name[ins]);
}

// //compute the matrix element of the conserved current between two propagators
// THREADABLE_FUNCTION_4ARG(calc_cur, spin1field*,cur, spincolor*,source, quad_su3*,conf, spincolor*,prop)
// {
//   GET_THREAD_ID();
  
//   vector_reset(cur);
  
//   communicate_lx_spincolor_borders(source);
//   communicate_lx_spincolor_borders(prop);
//   communicate_lx_quad_su3_borders(conf);
  
//   NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
//     for(int mu=0;mu<NDIM;mu++)
//       {
// 	// int ifw=loclx_neighup[ivol][mu];
// 	// spincolor f;
	
// 	// //piece psi_ivol U_ivol psi_fw
// 	// unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ifw]);
// 	// spincolor_scalar_prod(cur[ivol][mu],source[ivol],f);
	
// 	int ifw=loclx_neighup[ivol][mu];
// 	spincolor f,Gf;
// 	complex c;
	
// 	//piece psi_ivol U_ivol psi_fw
// 	unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ifw]);
// 	spincolor_copy(Gf,f);
// 	dirac_subt_the_prod_spincolor(Gf,base_gamma+igamma_of_mu[mu],f);
// 	spincolor_scalar_prod(c,source[ivol],Gf);
// 	complex_summ_the_prod_idouble(cur[ivol][mu],c,-0.5);
	
// 	//piece psi_fw U_ivol^dag psi_ivol
// 	unsafe_su3_dag_prod_spincolor(f,conf[ivol][mu],prop[ivol]);
// 	spincolor_copy(Gf,f);
// 	dirac_summ_the_prod_spincolor(Gf,base_gamma+igamma_of_mu[mu],f);
// 	spincolor_scalar_prod(c,source[ifw],Gf);
// 	complex_summ_the_prod_idouble(cur[ivol][mu],c,+0.5);
//       }
//   set_borders_invalid(cur);
// }
// THREADABLE_FUNCTION_END

// void in_main_Marcus(int narg,char **arg)
// {
//   init_grid(16,8);
  
//   spincolor *source=nissa_malloc("source",loc_vol,spincolor);
//   spincolor *prop=nissa_malloc("prop",loc_vol,spincolor);
//   quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  
//   read_ildg_gauge_conf(conf,"conf");
//   momentum_t th={1,0,0,0};
//   put_boundaries_conditions(conf,th,0,0);
  
//   read_real_vector(source,"source.0000.00000","scidac-binary-data");
//   read_real_vector(prop,"source.0000.00000.inverted","scidac-binary-data");
  
//   spin1field *cur=nissa_malloc("cur",loc_vol,spin1field);
//   calc_cur(cur,source,conf,prop);
  
//   FILE *fout=open_file("cur.txt","w");
//   master_fprintf(fout,"# [loops_em]  norm square of source     = %.16lg\n",double_vector_glb_norm2(source,loc_vol));
//   master_fprintf(fout,"# [loops_em]  norm square of propagator = %.16lg\n",double_vector_glb_norm2(prop,loc_vol));
//   NISSA_LOC_VOL_LOOP(ivol)
//   {
//     master_fprintf(fout,"# [loops_em] x\t");
//     for(int mu=0;mu<NDIM;mu++) master_fprintf(fout,"%d ",loc_coord_of_loclx[ivol][mu]);
//     master_fprintf(fout,"\n");
    
//     for(int mu=0;mu<NDIM;mu++) master_fprintf(fout,"%d" "\t" "%.16lg" "\t" "%.16lg\n",mu,cur[ivol][mu][RE],cur[ivol][mu][IM]);
//   }
//   close_file(fout);
  
//   nissa_free(cur);
//   nissa_free(conf);
//   nissa_free(prop);
//   nissa_free(source);
// }

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
  set_diluted_spin(0);
  set_diluted_color(0);
  set_diluted_space(1);
  
  //Number of hits
  read_nhits();
  
  //put the source in the list
  int store_source=false;
  ori_source_name_list.push_back(source_name);
  Q[source_name].init_as_source(noise_type,ALL_TIMES,0,store_source);
  
  //Read whether this is a twisted and/or clover run
  read_twisted_run();
  read_clover_run();
  
  //Read if to compute EU6 with two inversions
  read_str_int("ComputeEU6Alt",&compute_EU6_alt);
  
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
      
      //includes the propagators needed to compute EU6 with two inversions
      if(compute_EU6_alt)
	{
	  std::string prop_name_F=combine("Q%d_PROP_F",iquark);
	  std::string prop_name_F_PROP=combine("Q%d_PROP_F_PROP",iquark);
	  std::string prop_name_F_PROP_F=combine("Q%d_PROP_F_PROP_F",iquark);
	  for(auto& i : std::vector<std::tuple<std::string,std::string,insertion_t>>{
	      {prop_name_F,prop_name,PHOTON},
	      {prop_name_F_PROP,prop_name_F,PROP},
	      {prop_name_F_PROP_F,prop_name_F_PROP,PHOTON}})
	    {
	      std::string out=std::get<0>(i);
	      std::string in=std::get<1>(i);
	      insertion_t ins=std::get<2>(i);
	      qprop_name_list.push_back(out);
	      Q[out].init_as_propagator(ins,{{in,{1.0,0.0}}},ALL_TIMES,residue,kappa,mass,r,charge,theta,store_prop);
	    }
	  mes2pts_contr_map.push_back(mes_contr_map_t("EU6",source_name,prop_name_F_PROP_F));
	}
    }
  
  mes_gamma_list.push_back(idirac_pair_t(5,5));
  
  read_photon_pars();
  
  read_free_theory_flag();
  read_random_gauge_transform();
  read_Landau_gauge_fix();
  read_store_conf();
  
  read_ngauge_conf();
  
  /////////////////////////////////////////////////////////////////
  
  allocate_loop_source();
  curr::allocate_all_fields();
  hits::allocate();
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
  
  hits::free();
  
  if(clover_run)
    {
      nissa_free(Cl);
      nissa_free(invCl);
    }
}

//compute pseudoscalar, scalar and tadpoles (which can be used for EU1, EU2, EU4)
void compute_disco_PST(int ihit)
{
  //compute "2pts"
  vector_reset(mes2pts_contr);
  compute_mes2pts_contr(false);
  
  //print correlators
  int force_append(ihit>0);
  int skip_inner_header=true;
  std::string hit_header=combine("\n # hit %d\n\n",ihit);
  print_mes2pts_contr(1.0,force_append,skip_inner_header,hit_header);
}

//compute j_{f,mu}^i
THREADABLE_FUNCTION_0ARG(compute_all_quark_currents)
{
  for(int iquark=0;iquark<nquarks;iquark++)
    {
      //select name and field
      std::string quark_name=get_prop_name(iquark,PROP);
      spin1field *j=curr::fields[ifield_idx(iquark,curr::j)];
      spin1field *J=curr::fields[ifield_idx(iquark,curr::J)];
      
      //compute and summ
      int revert=false;
      local_or_conserved_vector_current_mel(j,base_gamma[0],source_name,quark_name,revert);
      double_vector_summassign((double*)J,(double*)j,loc_vol*sizeof(spin1field)/sizeof(double));
    }
}
THREADABLE_FUNCTION_END

//take the scalar propduct of j or J with xi to get hits::f,f
THREADABLE_FUNCTION_0ARG(compute_all_E_f1_f2)
{
  GET_THREAD_ID();
  
  for(int isingle_or_all=0;isingle_or_all<hits::nsingle_or_all;isingle_or_all++)
    {
      curr::field_t ic=(isingle_or_all?curr::J:curr::j);
      
      //convolve
      for(int iquark=0;iquark<nquarks;iquark++)
	{
	  spin1field *c=curr::fields[ifield_idx(iquark,ic)];
	  spin1field *xi=curr::fields[ifield_idx(iquark,curr::xi)];
	  
	  multiply_by_tlSym_gauge_propagator(xi,c,photon);
	}
      
      //scalar product
      for(int iquark=0;iquark<nquarks;iquark++)
	for(int jquark=iquark;jquark<nquarks;jquark++)
	  {
	    spin1field *c=curr::fields[ifield_idx(iquark,ic)];
	    spin1field *xi=curr::fields[ifield_idx(jquark,curr::xi)];
	    
	    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	      {
		complex_put_to_zero(hits::field[ivol]);
		for(int mu=0;mu<NDIM;mu++) complex_summ_the_prod(hits::field[ivol],c[ivol][mu],xi[ivol][mu]);
	      }
	    THREAD_BARRIER();
	    
	    //collapse
	    complex glb_E_f1_f2;
	    complex_vector_glb_collapse(glb_E_f1_f2,hits::field,loc_vol);
	    for(int isw=0;isw<2;isw++)
	      complex_copy(hits::f1_f2[hits::idx(isw?iquark:jquark,
						isw?jquark:iquark,
						nhits_done_so_far,isingle_or_all)],glb_E_f1_f2);
	    //master_printf("E_%d_%d single_or_all %d= ( %lg , %lg )\n",iquark,jquark,isingle_or_all,glb_E_f1_f2[RE],glb_E_f1_f2[IM]);
	  }
    }
}
THREADABLE_FUNCTION_END

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
      
      //start at nhits_done_so_far
      skip_hits_done_so_far();
      curr::load_all_or_reset();
      hits::load(nhits_done_so_far);
      
      while(nhits_done_so_far<nhits)
	{
	  start_hit(nhits_done_so_far);
	  generate_propagators(nhits_done_so_far);
	  compute_disco_PST(nhits_done_so_far);
	  
	  compute_all_quark_currents();
	  compute_all_E_f1_f2();
	  
	  compute_EU1_EU2_EU4_EU6alt();
	  
	  nhits_done_so_far++;
	}
      
      curr::store_all();
      hits::store(nhits_done_so_far);
      write_nhits_done_so_far();
      compute_EU5_EU6();
      
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
