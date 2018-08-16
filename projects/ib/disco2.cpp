#include "nissa.hpp"

using namespace nissa;

const int ALL_TIMES=-1;
momentum_t theta={-1,0,0,0};

namespace mel
{
  //buffer for reduction
  complex *buffer;
  
  //compute the local matrix element between source and prop of gamma[igamma]
  THREADABLE_FUNCTION_4ARG(local_mel, double*,out, spincolor*,source, int,igamma, spincolor*,prop)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	spincolor t;
	unsafe_dirac_prod_spincolor(t,base_gamma+igamma,prop[ivol]);
	spincolor_scalar_prod(buffer[ivol],source[ivol],t);
      }
    THREAD_BARRIER();
    
    complex_vector_glb_collapse(out,buffer,loc_vol);
  }
  THREADABLE_FUNCTION_END
  
  //compute the matrix element of the conserved current between two propagators. If asking to revert, g5 is inserted between the two propagators
  THREADABLE_FUNCTION_5ARG(conserved_vector_current_mel, spin1field*,out, spincolor*,source, quad_su3*,conf, int,r, spincolor*,prop)
  {
    GET_THREAD_ID();
    
    vector_reset(out);
    
    //compute the gammas
    dirac_matr GAMMA[5];
    dirac_prod_idouble(GAMMA+4,base_gamma+5,-tau3[r]);
    for(int mu=0;mu<NDIM;mu++) GAMMA[mu]=base_gamma[igamma_of_mu[mu]];
    
    communicate_lx_spincolor_borders(source);
    communicate_lx_spincolor_borders(prop);
    communicate_lx_quad_su3_borders(conf);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  int ivol_fw=loclx_neighup[ivol][mu];
	  spincolor f,Gf;
	  complex c;
	  
	  //piece psi_ivol U_ivol psi_fw
	  unsafe_su3_prod_spincolor(f,conf[ivol][mu],prop[ivol_fw]);
	  unsafe_dirac_prod_spincolor(Gf,GAMMA+4,f);
	  dirac_subt_the_prod_spincolor(Gf,GAMMA+mu,f);
	  spincolor_scalar_prod(c,source[ivol],Gf);
	  complex_summ_the_prod_idouble(out[ivol][mu],c,-0.5);
		
	  //piece psi_fw U_ivol^dag psi_ivol
	  unsafe_su3_dag_prod_spincolor(f,conf[ivol][mu],prop[ivol]);
	  unsafe_dirac_prod_spincolor(Gf,GAMMA+4,f);
	  dirac_summ_the_prod_spincolor(Gf,GAMMA+mu,f);
	  spincolor_scalar_prod(c,source[ivol_fw],Gf);
	  complex_summ_the_prod_idouble(out[ivol][mu],c,+0.5);
	}
  }
  THREADABLE_FUNCTION_END
  
  //compute the summ of the product of the two vectors
  THREADABLE_FUNCTION_3ARG(global_product, double*,out, spin1field*,a, spin1field*,b)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	complex_put_to_zero(buffer[ivol]);
	for(int mu=0;mu<NDIM;mu++)
	  complex_summ_the_prod(buffer[ivol],a[ivol][mu],b[ivol][mu]);
      }
    THREAD_BARRIER();
    
    complex_vector_glb_collapse(out,buffer,loc_vol);
  }
  THREADABLE_FUNCTION_END
}

void in_main(int narg,char **arg)
{
  std::string input_path;
  
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
  
  //geometry
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //Wall time
  double wall_time;
  read_str_double("WallTime",&wall_time);
  
  //local random generator
  int seed;
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);
  
  //fermion
  int r=0;
  //Read kappa
  double kappa;
  read_str_double("Kappa",&kappa);
  double am;
  read_str_double("Mass",&am);
  double residue;
  read_str_double("Residue",&residue);
  
  //allocate the source and prop
  int nhits;
  read_str_int("NHits",&nhits);
  spincolor *eta[nhits];
  spincolor *phi[nhits];
  for(int ihit=0;ihit<nhits;ihit++)
    {
      eta[ihit]=nissa_malloc("eta",loc_vol+bord_vol,spincolor);
      phi[ihit]=nissa_malloc("phi",loc_vol+bord_vol,spincolor);
    }
  
  //photon
  gauge_info photon_pars;
  photon_pars.alpha=FEYNMAN_ALPHA;
  photon_pars.c1=WILSON_C1;
  photon_pars.zms=UNNO_ALEMANNA;
  
  //compute the tadpole coefficient
  momentum_t tadpole_coeff;
  compute_tadpole(tadpole_coeff,photon_pars);
  
  //free theory
  int free_theory;
  read_str_int("FreeTheory",&free_theory);
  
  //conf
  int nconfs;
  read_str_int("NGaugeConfs",&nconfs);
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  
  //currents
  spin1field *J[nhits];
  for(int ihit=0;ihit<nhits;ihit++) J[ihit]=nissa_malloc("J",loc_vol+bord_vol,spin1field);
  spin1field *xi=nissa_malloc("xi",loc_vol+bord_vol,spin1field);
  
  //source and solution for the solver
  spincolor *source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  spincolor *solution=nissa_malloc("solution",loc_vol+bord_vol,spincolor);
  
  //buffer for local matrix element
  mel::buffer=nissa_malloc("loc_mel::buffer",loc_vol,complex);
  
  //propagator used for tadpole
  spincolor *tadpole_prop=nissa_malloc("tadpole_prop",loc_vol+bord_vol,spincolor);
  
  /////////////////////////////////////////////////////////////////
  
  int iconf=0;
  do
    {
      //input conf and output dir
      char conf_path[1024];
      read_str(conf_path,1024);
      char outfolder[1024];
      read_str(outfolder,1024);
      
      //generate the source
      for(int ihit=0;ihit<nhits;ihit++)
	generate_undiluted_source(eta[ihit],RND_Z4,ALL_TIMES);
      
      if(file_exists(combine("%s/running",outfolder))) master_printf(" Skipping %s\n",conf_path);
      else
	{
	  int ris=create_dir(outfolder);
	  if(ris==0) master_printf(" Output path \"%s\" not present, created.\n",outfolder);
	  else       crash(" Failed to create the output \"%s\" for conf \"%s\".",outfolder,conf_path);
	  
	  //read the configuration and put phases
	  if(free_theory) generate_cold_lx_conf(conf);
	  else            read_ildg_gauge_conf(conf,conf_path);
	  
	  momentum_t old_theta;
	  old_theta[0]=0;old_theta[1]=old_theta[2]=old_theta[3]=0;
	  adapt_theta(conf,old_theta,theta,0,0);
	  
	  //compute all propagators
	  for(int ihit=0;ihit<nhits;ihit++)
	    {
	      master_printf("Prop Hit %d\n",ihit);
	      
	      safe_dirac_prod_spincolor(source,(tau3[r]==-1)?&Pminus:&Pplus,eta[ihit]);
	      
	      if(free_theory)
		{
		  tm_quark_info qu(kappa,am,r,theta);
		  tm_basis_t basis=WILSON_BASE;
		  multiply_from_left_by_x_space_twisted_propagator_by_fft(solution,source,qu,basis,false);
		}
	      else
		inv_tmD_cg_eoprec(solution,NULL,conf,kappa,am*tau3[r],1000000,residue,source);
	      
	      safe_dirac_prod_spincolor(phi[ihit],(tau3[r]==-1)?&Pminus:&Pplus,solution);
	    }
	  
	  //compute all currents
	  for(int ihit=0;ihit<nhits;ihit++)
	    {
	      //master_printf("Cur Hit %d\n",ihit);
	      mel::conserved_vector_current_mel(J[ihit],eta[ihit],conf,r,phi[ihit]);
	    }
	  
	  //compute diagrams EU1, EU2 and EU4
	  complex EU1={0.0,0.0},EU2={0.0,0.0},EU4={0.0,0.0};
	  
	  //open the output files
	  FILE *fout_EU1=open_file(combine("%s/EU1",outfolder),"w");
	  FILE *fout_EU2=open_file(combine("%s/EU2",outfolder),"w");
	  FILE *fout_EU4=open_file(combine("%s/EU4",outfolder),"w");
	  FILE *fout_EU5=open_file(combine("%s/EU5",outfolder),"w");
	  FILE *fout_EU6=open_file(combine("%s/EU6",outfolder),"w");
	  
	  for(int ihit=0;ihit<nhits;ihit++)
	    {
	      complex temp;
	      
	      //Pseudo
	      mel::local_mel(temp,eta[ihit],5,phi[ihit]);
	      complex_summ_the_prod_idouble(EU1,temp,-1.0);
	      master_fprintf(fout_EU1,"%.16lg %.16lg\n",EU1[RE]/(ihit+1),EU1[IM]/(ihit+1));
	      
	      //Scalar
	      mel::local_mel(temp,eta[ihit],0,phi[ihit]);
	      complex_summassign(EU2,temp);
	      master_fprintf(fout_EU2,"%.16lg %.16lg\n",EU2[RE]/(ihit+1),EU2[IM]/(ihit+1));
	      
	      //Tadpole
	      insert_tm_tadpole(tadpole_prop,conf,phi[ihit],r,tadpole_coeff,ALL_TIMES);
	      mel::local_mel(temp,eta[ihit],0,tadpole_prop);
	      complex_summassign(EU4,temp);
	      master_fprintf(fout_EU4,"%.16lg %.16lg\n",EU4[RE]/(ihit+1),EU4[IM]/(ihit+1));
	    }
	  
	  //Compute diagram EU5
	  complex EU5={0.0,0.0};
	  int nEU5=0;
	  for(int ihit=0;ihit<nhits;ihit++)
	    {
	      multiply_by_tlSym_gauge_propagator(xi,J[ihit],photon_pars);
	      
	      for(int jhit=0;jhit<ihit;jhit++)
		{
		  complex temp;
		  mel::global_product(temp,xi,J[jhit]);
		  complex_summassign(EU5,temp);
		  nEU5++;
		}
	      master_fprintf(fout_EU5,"%.16lg %.16lg %d %d\n",EU5[RE]/nEU5,EU5[IM]/nEU5,ihit,nEU5);
	    }
	  
	  //Compute diagram EU6
	  complex EU6={0.0,0.0};
	  int nEU6=0;
	  for(int ihit=0;ihit<nhits;ihit++)
	    {
	      for(int jhit=0;jhit<ihit;jhit++)
		{
		  mel::conserved_vector_current_mel(J[ihit],eta[ihit],conf,r,phi[jhit]);
		  mel::conserved_vector_current_mel(J[jhit],eta[jhit],conf,r,phi[ihit]);
		  multiply_by_tlSym_gauge_propagator(xi,J[ihit],photon_pars);
		  complex temp;
		  mel::global_product(temp,J[jhit],xi);
		  complex_summassign(EU6,temp);
		  nEU6++;
		}
	      master_fprintf(fout_EU6,"%.16lg %.16lg %d %d\n",EU6[RE]/nEU6,EU6[IM]/nEU6,ihit,nEU6);
	    }
	  
	  close_file(fout_EU1);
	  close_file(fout_EU2);
	  close_file(fout_EU4);
	  close_file(fout_EU5);
	  close_file(fout_EU6);
	}
      
      iconf++;
    }
  while(iconf<nconfs);
  
  /////////////////////////////////////////////////////////////////
  
  nissa_free(tadpole_prop);
  nissa_free(xi);
  nissa_free(solution);
  nissa_free(source);
  nissa_free(mel::buffer);
  for(int ihit=0;ihit<nhits;ihit++) nissa_free(J[ihit]);
  nissa_free(conf);
  
  //free the source and prop
  for(int ihit=0;ihit<nhits;ihit++)
    {
      nissa_free(eta[ihit]);
      nissa_free(phi[ihit]);
    }
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
