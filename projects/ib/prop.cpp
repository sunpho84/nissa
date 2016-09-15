#include <nissa.hpp>

#define EXTERN_PROP
 #include "prop.hpp"

namespace nissa
{
  //get a propagator inverting on "in"
  void get_qprop(spincolor *out,spincolor *in,double kappa,double mass,int r,double residue,double theta)
  {
    GET_THREAD_ID();
    
    //put theta on the boundaries
    adapt_spatial_theta(conf,theta);
    
    //rotate the source index - the propagator rotate AS the sign of mass term
    if(twisted_run) safe_dirac_prod_spincolor(in,(tau3[r]==-1)?&Pminus:&Pplus,in);
    
    //invert
    START_TIMING(inv_time,ninv_tot);
    
    if(clover_run) inv_tmclovD_cg_eoprec(out,NULL,conf,kappa,Cl,invCl,mass,1000000,residue,in);
    else inv_tmD_cg_eoprec(out,NULL,conf,kappa,mass,1000000,residue,in);
    
    STOP_TIMING(inv_time);
    
    //rotate the sink index
    if(twisted_run) safe_dirac_prod_spincolor(out,(tau3[r]==-1)?&Pminus:&Pplus,out);
    
    //put back no theta on the boundaries
    put_spatial_theta_periodic(conf);
  }
  
  //generate a source, wither a wall or a point in the origin
  THREADABLE_FUNCTION_1ARG(generate_original_source, qprop_t*,sou)
  {
    GET_THREAD_ID();
    
    //consistency check
    if(!stoch_source&&(!diluted_spi_source||!diluted_col_source)) crash("for a non-stochastic source, spin and color must be diluted");
    
    //reset all to begin
    for(int i=0;i<nso_spi*nso_col;i++) vector_reset(sou->sp[i]);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	spincolor c;
	spincolor_put_to_zero(c);
	
	//fill colour and spin index 0
	for(int id_si=0;id_si<(diluted_spi_source?1:NDIRAC);id_si++)
	  for(int ic_si=0;ic_si<(diluted_col_source?1:NCOL);ic_si++)
	    {
	      if(stoch_source && (sou->tins==-1||glb_coord_of_loclx[ivol][0]==sou->tins)) comp_get_rnd(c[id_si][ic_si],&(loc_rnd_gen[ivol]),sou->noise_type);
	      else if(glblx_of_loclx[ivol]==0) complex_put_to_real(c[id_si][ic_si],1);
	  }
	
	//fill other spin indices
	for(int id_so=0;id_so<nso_spi;id_so++)
	  for(int ic_so=0;ic_so<nso_col;ic_so++)
	    for(int id_si=0;id_si<NDIRAC;id_si++)
	      for(int ic_si=0;ic_si<NCOL;ic_si++)
		  if((!diluted_spi_source||(id_so==id_si))&&(!diluted_col_source||(ic_so==ic_si)))
		    complex_copy(sou->sp[so_sp_col_ind(id_so,ic_so)][ivol][id_si][ic_si],c[diluted_spi_source?0:id_si][diluted_col_source?0:ic_si]);
      }
    
    //compute the norm2, set borders invalid and if required, smear the conf
    double loc_ori_source_norm2=0;
    for(int id_so=0;id_so<nso_spi;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	{
	  spincolor *s=sou->sp[so_sp_col_ind(id_so,ic_so)];
	  set_borders_invalid(s);
	  loc_ori_source_norm2+=double_vector_glb_norm2(s,loc_vol);
	  if(sou->sme) gaussian_smearing(s,s,ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
	}
    if(IS_MASTER_THREAD) sou->ori_source_norm2=loc_ori_source_norm2;
  }
  THREADABLE_FUNCTION_END
  
  //insert the photon on the source side
  void insert_external_loc_source(spincolor *out,spin1field *curr,coords dirs,spincolor *in,int t)
  {
    GET_THREAD_ID();
    
    if(in==out) crash("in==out");
    
    vector_reset(out);
    
    for(int mu=0;mu<NDIM;mu++)
      if(dirs[mu])
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  if(t==-1||glb_coord_of_loclx[ivol][0]==t)
	    {
	      spincolor temp1,temp2;
	      unsafe_dirac_prod_spincolor(temp1,base_gamma+map_mu[mu],in[ivol]);
	      unsafe_spincolor_prod_complex(temp2,temp1,curr[ivol][mu]);
	      spincolor_summ_the_prod_idouble(out[ivol],temp2,1);
	    }
    
    set_borders_invalid(out);
  }
  
  //insert the photon on the source
  void insert_external_loc_source(spincolor *out,spin1field *curr,spincolor *in,int t)
  {insert_external_loc_source(out,curr,all_dirs,in,t);}
  
  void insert_external_source(spincolor *out,spin1field *curr,spincolor *ori,int t,int r,int loc)
  {
    if(loc) insert_external_loc_source(loop_source,curr,ori,t);
    else
      if(twisted_run) insert_tm_external_source(loop_source,conf,curr,ori,r,t);
      else            insert_Wilson_external_source(loop_source,conf,curr,ori,t);
  }
  
  //generate a sequential source
  void generate_source(insertion_t inser,int r,spincolor *ori,int t)
  {
    source_time-=take_time();
    
    switch(inser)
      {
      case ORI_SOU:prop_multiply_with_gamma(loop_source,0,ori,t);break;
      case SCALAR:prop_multiply_with_gamma(loop_source,0,ori,t);break;
      case PSEUDO:prop_multiply_with_gamma(loop_source,5,ori,t);break;
      case PHOTON:insert_external_source(loop_source,photon_field,ori,t,r,loc_hadr_curr);break;
      case PHOTON_PHI:insert_external_source(loop_source,photon_phi,ori,t,r,loc_hadr_curr);break;
      case PHOTON_ETA:insert_external_source(loop_source,photon_eta,ori,t,r,loc_hadr_curr);break;
      case TADPOLE:
	if(twisted_run) insert_tm_tadpole(loop_source,conf,ori,r,tadpole,t);
	else            insert_Wilson_tadpole(loop_source,conf,ori,tadpole,t);
	break;
	//case VECTOR:insert_external_source(source,NULL,ori,t,r,loc_pion_curr);break;
      }
    
    source_time+=take_time();
    nsource_tot++;
  }
  
  //generate all the quark propagators
  void generate_quark_propagators(int ihit)
  {
    GET_THREAD_ID();
    for(size_t i=0;i<qprop_name_list.size();i++)
      {
	//get names
	std::string name=qprop_name_list[i];
	qprop_t &q=Q[name];
	std::string source_name=q.source_name;
	qprop_t &qsource=Q[source_name];
	
	//copy norm
	q.ori_source_norm2=qsource.ori_source_norm2;
	
	//write info on mass and r
	if(twisted_run) master_printf(" mass[%d]=%lg, r=%d, theta=%lg\n",i,q.mass,q.r,q.theta);
	else            master_printf(" kappa[%d]=%lg, theta=%lg\n",i,q.kappa,q.theta);
	
	//compute the inverse clover term, if needed
	if(clover_run) invert_twisted_clover_term(invCl,q.mass,q.kappa,Cl);
	
	insertion_t insertion=q.insertion;
	master_printf("Generating propagator %s inserting %s on source %s\n",name.c_str(),ins_name[insertion],ins_name[(int)insertion],source_name.c_str());
	for(int id_so=0;id_so<nso_spi;id_so++)
	  for(int ic_so=0;ic_so<nso_col;ic_so++)
	    {
	      int isou=so_sp_col_ind(id_so,ic_so);
	      generate_source(insertion,q.r,qsource[isou],q.tins);
	      spincolor *sol=q[isou];
	      
	      //combine the filename
	      std::string path=combine("%s/hit%d_prop%s_idso%d_icso%d",outfolder,ihit,name.c_str(),id_so,ic_so);
	      
	      //if the prop exists read it
	      if(file_exists(path))
		{
		  master_printf("  loading the inversion dirac index %d, color %d\n",id_so,ic_so);
		  START_TIMING(read_prop_time,nread_prop);
		  read_real_vector(sol,path,"prop");
		  STOP_TIMING(read_prop_time);
		}
	      else
		{
		  //otherwise compute it and possibly store it
		  get_qprop(sol,loop_source,q.kappa,q.mass,q.r,q.residue,q.theta);
		  if(q.sme) gaussian_smearing(sol,sol,ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
		  
		  if(q.store)
		    {
		      START_TIMING(store_prop_time,nstore_prop);
		      write_double_vector(path,sol,64,"prop");
		      STOP_TIMING(store_prop_time);
		    }
		  master_printf("  finished the inversion dirac index %d, color %d\n",id_so,ic_so);
		}
	    }
      }
  }
  
  /////////////////////////////////////////////// lepton propagators ///////////////////////////////////////////
  
  //allocate all leptonic propagators
  void allocate_L_prop()
  {
    nlprop=ilprop(nquark_lep_combos-1,nlins-1,norie-1,nr_lep-1)+1;
    L=nissa_malloc("L*",nlprop,spinspin*);
    for(int iprop=0;iprop<nlprop;iprop++) L[iprop]=nissa_malloc("L",loc_vol+bord_vol,spinspin);
  }
  
  //free all leptonic propagators
  void free_L_prop()
  {
    for(int iprop=0;iprop<nlprop;iprop++) nissa_free(L[iprop]);
    nissa_free(L);
    nissa_free(temp_lep);
  }
  
  //return appropriately modified info
  tm_quark_info get_lepton_info(int ilepton,int orie,int r)
  {
    tm_quark_info le=leps[ilepton];
    le.r=r;
    le.bc[0]=QUARK_BOUND_COND;
    for(int i=1;i<NDIM;i++) le.bc[i]*=sign_orie[orie];
    
    return le;
  }
  
  //compute phase exponent for space part: vec{p}*\vec{x}
  double get_space_arg(int ivol,momentum_t bc)
  {
    double arg=0;
    for(int mu=1;mu<NDIM;mu++)
      {
	double step=bc[mu]*M_PI/glb_size[mu];
	arg+=step*glb_coord_of_loclx[ivol][mu];
      }
    return arg;
  }
  
  //compute the phase for lepton on its sink
  void get_lepton_sink_phase_factor(complex out,int ivol,int ilepton,tm_quark_info le)
  {
    //compute space and time factor
    double arg=get_space_arg(ivol,le.bc);
    int t=glb_coord_of_loclx[ivol][0];
    if(follow_chris_or_nazario==follow_nazario && t>=glb_size[0]/2) t=glb_size[0]-t;
    double ext=exp(t*lep_energy[ilepton]);
    
    //compute full exponential (notice the factor -1)
    out[RE]=cos(-arg)*ext;
    out[IM]=sin(-arg)*ext;
  }
  
  //compute the phase for antineutrino - the orientation is that of the muon (as above)
  void get_antineutrino_source_phase_factor(complex out,int ivol,int ilepton,momentum_t bc)
  {
    //compute space and time factor
    double arg=get_space_arg(ivol,bc);
    int t=glb_coord_of_loclx[ivol][0];
    if(follow_chris_or_nazario==follow_nazario && t>=glb_size[0]/2) t=glb_size[0]-t;
    double ext=exp(t*neu_energy[ilepton]);
    
    //compute full exponential (notice the factor +1)
    out[RE]=cos(+arg)*ext;
    out[IM]=sin(+arg)*ext;
  }
  
  //set everything to a phase factor
  void set_to_lepton_sink_phase_factor(spinspin *prop,int ilepton,tm_quark_info &le)
  {
    GET_THREAD_ID();
    
    vector_reset(prop);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	complex ph;
	get_lepton_sink_phase_factor(ph,ivol,ilepton,le);
	spinspin_put_to_diag(prop[ivol],ph);
      }
    set_borders_invalid(prop);
  }
  
  //insert the photon on the source side
  THREADABLE_FUNCTION_5ARG(insert_photon_on_the_source, spinspin*,prop, spin1field*,A, int*,dirs, tm_quark_info,le, int,twall)
  {
    GET_THREAD_ID();
    
    //select A
    communicate_lx_spin1field_borders(A);
    
    //copy on the temporary and communicate borders
    vector_copy(temp_lep,prop);
    communicate_lx_spinspin_borders(temp_lep);
    vector_reset(prop);
    
    if(!loc_muon_curr)
      {
	dirac_matr GAMMA;
	if(twisted_run) dirac_prod_double(&GAMMA,base_gamma+0,1);
	else dirac_prod_idouble(&GAMMA,base_gamma+5,-tau3[le.r]);
	
	//phases
	quad_u1 phases;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    phases[mu][0]=cos(le.bc[mu]*M_PI);
	    phases[mu][1]=sin(le.bc[mu]*M_PI);
	  }
	
	//prepare each propagator for a single lepton
	//by computing i(phi(x-mu)A_mu(x-mu)(-i t3 g5-gmu)/2-phi(x+mu)A_mu(x)(-i t3 g5+gmu)/2)=
	//(ph0 A_mu(x-mu)g[r][0][mu]-ph0 A_mu(x)g[r][1][mu])=
	for(int mu=0;mu<NDIM;mu++)
	  if(dirs[mu])
	    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	      if(twall==-1||glb_coord_of_loclx[ivol][0]==twall)
		{
		  //find neighbors
		  int ifw=loclx_neighup[ivol][mu];
		  int ibw=loclx_neighdw[ivol][mu];
		  
		  //compute phase factor
		  spinspin ph_bw,ph_fw;
		  
		  //transport down and up
		  if(glb_coord_of_loclx[ivol][mu]==glb_size[mu]-1) unsafe_spinspin_prod_complex_conj2(ph_fw,temp_lep[ifw],phases[mu]);
		  else spinspin_copy(ph_fw,temp_lep[ifw]);
		  if(glb_coord_of_loclx[ivol][mu]==0) unsafe_spinspin_prod_complex(ph_bw,temp_lep[ibw],phases[mu]);
		  else spinspin_copy(ph_bw,temp_lep[ibw]);
		  
		  //fix coefficients - i is inserted here!
		  //also dir selection is made here
		  spinspin_prodassign_idouble(ph_fw,-0.5*dirs[mu]);
		  spinspin_prodassign_idouble(ph_bw,+0.5*dirs[mu]);
		  
		  //fix insertion of the current
		  safe_spinspin_prod_complex(ph_fw,ph_fw,A[ivol][mu]);
		  safe_spinspin_prod_complex(ph_bw,ph_bw,A[ibw][mu]);
		  
		  //summ and subtract the two
		  spinspin fw_M_bw,fw_P_bw;
		  spinspin_subt(fw_M_bw,ph_fw,ph_bw);
		  spinspin_summ(fw_P_bw,ph_fw,ph_bw);
		  
		  //put GAMMA on the summ
		  spinspin temp_P;
		  unsafe_spinspin_prod_dirac(temp_P,fw_P_bw,&GAMMA);
		  spinspin_summassign(prop[ivol],temp_P);
		  
		  //put gmu on the diff
		  spinspin temp_M;
		  unsafe_spinspin_prod_dirac(temp_M,fw_M_bw,base_gamma+map_mu[mu]);
		  spinspin_summassign(prop[ivol],temp_M);
		}
      }
    else
      {
	for(int mu=0;mu<NDIM;mu++)
	  if(dirs[mu])
	    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	      if(twall==-1||glb_coord_of_loclx[ivol][0]==twall)
		{
		  spinspin temp1,temp2;
		  unsafe_spinspin_prod_dirac(temp1,temp_lep[ivol],base_gamma+map_mu[mu]);
		  unsafe_spinspin_prod_complex(temp2,temp1,A[ivol][mu]);
		  spinspin_summ_the_prod_idouble(prop[ivol],temp2,1);
		}
	
      }
    
    set_borders_invalid(prop);
  }
  THREADABLE_FUNCTION_END
  
  //generate all the lepton propagators, pointing outward
  //the computations is done by:
  // 1)putting the correct phase in x space, given by exp(E_mu*t-i*vec(p)*vec(x))
  // 2)multiplying it by the conserved current inserting eta or phi
  // 3)going to momentum space
  // 4)multiplying by the lepton propagator in momentum space
  // 5)coming back to x space
  THREADABLE_FUNCTION_0ARG(generate_lepton_propagators)
  {
    GET_THREAD_ID();
    
    if(IS_MASTER_THREAD) lepton_prop_time-=take_time();
    
    for(int ilepton=0;ilepton<nquark_lep_combos;ilepton++)
      for(int ilins=0;ilins<nlins;ilins++)
	for(int ori=0;ori<norie;ori++)
	  for(int r=0;r<nr_lep;r++)
	    {
	      //set the properties of the meson
	      //time boundaries are anti-periodic, space are as for external line
	      tm_quark_info le=get_lepton_info(ilepton,ori,r);
	      
	      //select the propagator
	      int iprop=ilprop(ilepton,ilins,ori,r);
	      spinspin *prop=L[iprop];
	      
	      //put it to a phase
	      set_to_lepton_sink_phase_factor(prop,ilepton,le);
	      
	      //if we are doing Nazario's way (with the external line) add it
	      if(follow_chris_or_nazario==follow_nazario)
		{
		  //select only the wall
		  int tmiddle=glb_size[0]/2;
		  select_propagator_timeslice(prop,prop,tmiddle);
		  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le,base);
		}
	      
	      //insert or not photon
	      if(ilins)
		{
		  //insert photon and prolong
		  insert_photon_on_the_source(prop,photon_field,all_dirs,le,-1); //all times
		  multiply_from_right_by_x_space_twisted_propagator_by_fft(prop,prop,le,base);
		}
	    }
    
    if(IS_MASTER_THREAD) lepton_prop_time+=take_time();
  }
  THREADABLE_FUNCTION_END
  
  /////////////////////////////////////////////// photon propagators ///////////////////////////////////////////
  
  //allocate the photon fields
  void allocate_photon_fields()
  {
    photon_eta=nissa_malloc("photon_eta",loc_vol+bord_vol,spin1field);
    photon_field=nissa_malloc("photon_field",loc_vol+bord_vol,spin1field);
    photon_phi=nissa_malloc("photon_phi",loc_vol+bord_vol,spin1field);
  }
  
  //free the photon fields
  void free_photon_fields()
  {
    nissa_free(photon_eta);
    nissa_free(photon_phi);
    nissa_free(photon_field);
  }
  
  //wrapper to generate a stochastic propagator
  void generate_photon_stochastic_propagator()
  {
    photon_prop_time-=take_time();
    
    //generate source and stochastich propagator
    generate_stochastic_tlSym_gauge_propagator(photon_phi,photon_eta,photon);
    
    //generate the photon field
    multiply_by_sqrt_tlSym_gauge_propagator(photon_field,photon_eta,photon);
    
    photon_prop_time+=take_time();
    nphoton_prop_tot++;
  }
}
