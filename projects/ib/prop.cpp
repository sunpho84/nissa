#include <nissa.hpp>

#define EXTERN_PROP
 #include "prop.hpp"

#include <set>

namespace nissa
{
  //multiply the configuration for an additional u(1) field, defined as exp(-i e q A /3)
  THREADABLE_FUNCTION_2ARG(add_photon_field_to_conf, quad_su3*,conf, double,charge)
  {
    const double alpha_em=1/137.04;
    const double e2=4*M_PI*alpha_em;
    const double e=sqrt(e2);
    verbosity_lv2_master_printf("Adding backfield, for a quark of charge q=e*%lg/3\n",charge);
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  complex ph;
	  complex_iexp(ph,-e*photon_field[ivol][mu][RE]*charge/3.0);
	  safe_su3_prod_complex(conf[ivol][mu],conf[ivol][mu],ph);
	}
    set_borders_invalid(conf);
  }
  THREADABLE_FUNCTION_END
  
  //remove the field
  void rem_photon_field_to_conf(quad_su3 *conf,double q)
  {add_photon_field_to_conf(conf,-q);}
  
  //get a propagator inverting on "in"
  void get_qprop(spincolor *out,spincolor *in,double kappa,double mass,int r,double charge,double residue,double theta)
  {
    GET_THREAD_ID();
    
    //rotate the source index - the propagator rotate AS the sign of mass term
    if(twisted_run) safe_dirac_prod_spincolor(in,(tau3[r]==-1)?&Pminus:&Pplus,in);
    
    //invert
    START_TIMING(inv_time,ninv_tot);
    
    quad_su3 *conf=get_updated_conf(charge,QUARK_BOUND_COND,theta,glb_conf);
    
    if(clover_run) inv_tmclovD_cg_eoprec(out,NULL,conf,kappa,Cl,invCl,glb_cSW,mass,1000000,residue,in);
    else inv_tmD_cg_eoprec(out,NULL,conf,kappa,mass,1000000,residue,in);
    
    STOP_TIMING(inv_time);
    
    //rotate the sink index
    if(twisted_run) safe_dirac_prod_spincolor(out,(tau3[r]==-1)?&Pminus:&Pplus,out);
  }
  
  //generate a source, wither a wall or a point in the origin
  THREADABLE_FUNCTION_1ARG(generate_original_source, qprop_t*,sou)
  {
    GET_THREAD_ID();
    
    //consistency check
    if(!stoch_source and (!diluted_spi_source or !diluted_col_source)) crash("for a non-stochastic source, spin and color must be diluted");
    
    //reset all to begin
    for(int i=0;i<nso_spi*nso_col;i++) vector_reset(sou->sp[i]);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	spincolor c;
	spincolor_put_to_zero(c);
	
	//compute relative coords
	bool is_spat_orig=true;
	coords rel_c;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    rel_c[mu]=rel_coord_of_loclx(ivol,mu);
	    if(mu) is_spat_orig&=(rel_c[mu]==0);
	  }
	
	//dilute in space
	int mask=1;
	for(int mu=0;mu<NDIM;mu++) mask&=(rel_c[mu]%diluted_spat_source==0);
	
	//fill colour and spin index 0
	for(int id_si=0;id_si<(diluted_spi_source?1:NDIRAC);id_si++)
	  for(int ic_si=0;ic_si<(diluted_col_source?1:NCOL);ic_si++)
	    {
	      if(stoch_source and mask and (sou->tins==-1 or rel_c[0]==sou->tins)) comp_get_rnd(c[id_si][ic_si],&(loc_rnd_gen[ivol]),sou->noise_type);
	      if(!stoch_source and is_spat_orig and (sou->tins==-1 or rel_c[0]==sou->tins)) complex_put_to_real(c[id_si][ic_si],1);
	    }
	
	//fill other spin indices
	for(int id_so=0;id_so<nso_spi;id_so++)
	  for(int ic_so=0;ic_so<nso_col;ic_so++)
	    for(int id_si=0;id_si<NDIRAC;id_si++)
	      for(int ic_si=0;ic_si<NCOL;ic_si++)
		  if((!diluted_spi_source or (id_so==id_si)) and (!diluted_col_source or (ic_so==ic_si)))
		    complex_copy(sou->sp[so_sp_col_ind(id_so,ic_so)][ivol][id_si][ic_si],c[diluted_spi_source?0:id_si][diluted_col_source?0:ic_si]);
      }
    
    //compute the norm2, set borders invalid
    double ori_source_norm2=0;
    for(int id_so=0;id_so<nso_spi;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	{
	  spincolor *s=sou->sp[so_sp_col_ind(id_so,ic_so)];
	  set_borders_invalid(s);
	  ori_source_norm2+=double_vector_glb_norm2(s,loc_vol);
	}
    if(IS_MASTER_THREAD) sou->ori_source_norm2=ori_source_norm2;
  }
  THREADABLE_FUNCTION_END
  
  //insert the photon on the source side
  void insert_external_loc_source(spincolor *out,spin1field *curr,spincolor *in,int t,coords dirs)
  {
    GET_THREAD_ID();
    
    if(in==out) crash("in==out");
    
    vector_reset(out);
    
    for(int mu=0;mu<NDIM;mu++)
      if(dirs[mu])
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  if(t==-1 or glb_coord_of_loclx[ivol][0]==t)
	    {
	      spincolor temp1,temp2;
	      unsafe_dirac_prod_spincolor(temp1,base_gamma+igamma_of_mu[mu],in[ivol]);
	      unsafe_spincolor_prod_complex(temp2,temp1,curr[ivol][mu]);
	      spincolor_summ_the_prod_idouble(out[ivol],temp2,1);
	    }
    
    set_borders_invalid(out);
  }
  
  //insert the external source
  void insert_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *ori,int t,int r,coords dirs,int loc)
  {
    if(loc) insert_external_loc_source(out,curr,ori,t,dirs);
    else
      if(twisted_run) insert_tm_external_source(out,conf,curr,ori,r,dirs,t);
      else            insert_Wilson_external_source(out,conf,curr,ori,dirs,t);
  }
  
  //insert the tadpole
  void insert_tadpole(spincolor *out,quad_su3 *conf,spincolor *ori,int t,int r)
  {
    if(twisted_run) insert_tm_tadpole(loop_source,conf,ori,r,tadpole,t);
    else            insert_Wilson_tadpole(loop_source,conf,ori,tadpole,t);
  }
  
  //insert the conserved current
  void insert_conserved_current(spincolor *out,quad_su3 *conf,spincolor *ori,int t,int r,coords dirs)
  {
    if(twisted_run) insert_tm_conserved_current(loop_source,conf,ori,r,dirs,t);
    else            insert_Wilson_conserved_current(loop_source,conf,ori,dirs,t);
  }
  
  //smear the propagator
  void smear_prop(spincolor *out,quad_su3 *conf,spincolor *ori,int t,double kappa,int nlevels)
  {
    //nb: the smearing radius is given by
    //a*sqrt(2*n*kappa/(1+6*kappa))
    
    gaussian_smearing(out,ori,conf,kappa,nlevels);
    select_propagator_timeslice(out,out,t);
  }
  
  //phase the propagator
  THREADABLE_FUNCTION_4ARG(phase_prop, spincolor*,out, spincolor*,ori, int,t, double,th)
  {
    GET_THREAD_ID();
    
    if(fabs((int)(th/2)-th/2)>1e-10) crash("Error: phase %lg must be an even integer",th);
    
    vector_reset(out);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute x*p
	double arg=0.0;
	for(int mu=1;mu<NDIM;mu++) arg+=2*M_PI*th*glb_coord_of_loclx[ivol][mu]/glb_size[mu]; //N.B: valid only if source is on origin...
	
	//compute exp(ip)
	complex factor;
	complex_iexp(factor,arg);
	
	//put the phase
	if(t==-1 or glb_coord_of_loclx[ivol][0]==t)
	  unsafe_spincolor_prod_complex(out[ivol],ori[ivol],factor);
      }
    
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //generate a sequential source
  void generate_source(insertion_t inser,int r,double charge,double kappa,double theta,spincolor *ori,int t)
  {
    source_time-=take_time();
    
    int rel_t=t;
    if(rel_t!=-1) rel_t=(t+source_coord[0])%glb_size[0];
    coords dirs[NDIM]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    
    quad_su3 *conf;
    if(inser!=SMEARING) conf=get_updated_conf(charge,QUARK_BOUND_COND,theta,glb_conf);
    else
      {
	quad_su3 *ext_conf;
	if(ape_smeared_conf) ext_conf=ape_smeared_conf;
	else                 ext_conf=glb_conf;
	conf=get_updated_conf(0.0,PERIODIC,theta,ext_conf);
      }
    
    master_printf("Inserting r: %d\n",r);
    switch(inser)
      {
      case PROP:prop_multiply_with_gamma(loop_source,0,ori,rel_t);break;
      case SCALAR:prop_multiply_with_gamma(loop_source,0,ori,rel_t);break;
      case PSEUDO:prop_multiply_with_gamma(loop_source,5,ori,rel_t);break;
      case PHOTON:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,all_dirs,loc_hadr_curr);break;
      case PHOTON0:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,dirs[0],loc_hadr_curr);break;
      case PHOTON1:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,dirs[1],loc_hadr_curr);break;
      case PHOTON2:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,dirs[2],loc_hadr_curr);break;
      case PHOTON3:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,dirs[3],loc_hadr_curr);break;
      case PHOTON_PHI:insert_external_source(loop_source,conf,photon_phi,ori,rel_t,r,all_dirs,loc_hadr_curr);break;
      case PHOTON_ETA:insert_external_source(loop_source,conf,photon_eta,ori,rel_t,r,all_dirs,loc_hadr_curr);break;
      case TADPOLE:insert_tadpole(loop_source,conf,ori,rel_t,r);break;
      case CVEC0:insert_conserved_current(loop_source,conf,ori,rel_t,r,dirs[0]);break;
      case CVEC1:insert_conserved_current(loop_source,conf,ori,rel_t,r,dirs[1]);break;
      case CVEC2:insert_conserved_current(loop_source,conf,ori,rel_t,r,dirs[2]);break;
      case CVEC3:insert_conserved_current(loop_source,conf,ori,rel_t,r,dirs[3]);break;
      case SMEARING:smear_prop(loop_source,conf,ori,rel_t,kappa,r);break;
      case PHASING:phase_prop(loop_source,ori,rel_t,theta);break;
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
	master_printf("Generating propagator %s inserting %s on source %s\n",name.c_str(),ins_name[insertion],source_name.c_str());
	for(int id_so=0;id_so<nso_spi;id_so++)
	  for(int ic_so=0;ic_so<nso_col;ic_so++)
	    {
	      int isou=so_sp_col_ind(id_so,ic_so);
	      generate_source(insertion,q.r,q.charge,q.kappa,q.theta,qsource[isou],q.tins);
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
		  //otherwise compute it
		  if(q.insertion==PROP) get_qprop(sol,loop_source,q.kappa,q.mass,q.r,q.charge,q.residue,q.theta);
		  else                  vector_copy(sol,loop_source);
		  
		  //and store if needed
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
	arg+=step*rel_coord_of_loclx(ivol,mu);
      }
    return arg;
  }
  
  //compute the phase for lepton on its sink
  void get_lepton_sink_phase_factor(complex out,int ivol,int ilepton,tm_quark_info le)
  {
    //compute space and time factor
    double arg=get_space_arg(ivol,le.bc);
    int t=rel_time_of_loclx(ivol);
    if(follow_chris_or_nazario==follow_nazario and t>=glb_size[0]/2) t=glb_size[0]-t;
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
    int t=rel_time_of_loclx(ivol);
    if(follow_chris_or_nazario==follow_nazario and t>=glb_size[0]/2) t=glb_size[0]-t;
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
	      if(twall==-1 or rel_time_of_loclx(ivol)==twall)
		{
		  //find neighbors
		  int ifw=loclx_neighup[ivol][mu];
		  int ibw=loclx_neighdw[ivol][mu];
		  
		  //compute phase factor
		  spinspin ph_bw,ph_fw;
		  
		  //transport down and up
		  if(rel_coord_of_loclx(ivol,mu)==glb_size[mu]-1) unsafe_spinspin_prod_complex_conj2(ph_fw,temp_lep[ifw],phases[mu]);
		  else spinspin_copy(ph_fw,temp_lep[ifw]);
		  if(rel_coord_of_loclx(ivol,mu)==0) unsafe_spinspin_prod_complex(ph_bw,temp_lep[ibw],phases[mu]);
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
		  unsafe_spinspin_prod_dirac(temp_M,fw_M_bw,base_gamma+igamma_of_mu[mu]);
		  spinspin_summassign(prop[ivol],temp_M);
		}
      }
    else
      {
	for(int mu=0;mu<NDIM;mu++)
	  if(dirs[mu])
	    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	      if(twall==-1 or rel_time_of_loclx(ivol)==twall)
		{
		  spinspin temp1,temp2;
		  unsafe_spinspin_prod_dirac(temp1,temp_lep[ivol],base_gamma+igamma_of_mu[mu]);
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
    
    //invalidate internal conf
    inner_conf_valid=false;
    
    photon_prop_time+=take_time();
    nphoton_prop_tot++;
  }
  
  //put the phase of the source due to missing e(iky)
  THREADABLE_FUNCTION_2ARG(put_fft_source_phase, spincolor*,qtilde, double,fft_sign)
  {
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
      {
	//sum_mu -sign*2*pi*p_mu*y_mu/L_mu
	double arg=0;
	for(int mu=0;mu<NDIM;mu++) arg+=-fft_sign*2*M_PI*glb_coord_of_loclx[imom][mu]*source_coord[mu]/glb_size[mu];
	
	complex f={cos(arg),sin(arg)};
	spincolor_prodassign_complex(qtilde[imom],f);
	
	// spincolor_put_to_zero(qtilde[imom]);
	// for(int mu=0;mu<4;mu++) qtilde[imom][mu][0][0]=glb_coord_of_loclx[imom][mu];
      }
    
    set_borders_invalid(qtilde);
  }
  THREADABLE_FUNCTION_END
  
  //initialize the fft filter, once forever
  void init_fft_filter()
  {
    master_printf("Initializing fft filter\n");
    
    //file where to store output
    FILE *fout=NULL;
    const char path_list[]="mom_list.txt";
    if(not file_exists("mom_list.txt")) fout=open_file(path_list,"w");
    
    //store the list of filtered
    std::set<int> list_of_filtered;;
    
    //scattering list
    all_to_all_scattering_list_t sl;
    for(std::vector<fft_mom_range_t>::iterator f=fft_mom_range_list.begin();f!=fft_mom_range_list.end();f++)
      for(int vol=vol_of_lx(f->width),ifilt=0;ifilt<vol;ifilt++)
	{
	  //gets the coordinate in the filtering volume
	  coords c;
	  coord_of_lx(c,ifilt,f->width);
	  coord_summassign(c,f->offs,glb_size);
	  
	  //compute p~4/p~2^2
	  double pt2=0,pt4=0;
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      double pmu=M_PI*(2*c[mu]+(mu==0)*QUARK_BOUND_COND)/glb_size[mu];
	      double ptmu=sin(pmu);
	      pt2+=sqr(ptmu);
	      pt4+=pow(ptmu,4.0);
	    }
	  
	  if(pt4/sqr(pt2)<p4_fr_p22_max)
	    for(int imir=0;imir<pow(2,NDIM);imir++)
	      {
		//get mirrorized
		coords cmir;
		for(int mu=0;mu<NDIM;mu++)
		cmir[mu]=get_mirrorized_site_coord(c[mu]+(mu==0 and get_bit(imir,0) and QUARK_BOUND_COND==1),mu,get_bit(imir,mu));
		
		//check if not already collected
		int iglb=glblx_of_coord(cmir);
		if(list_of_filtered.find(iglb)==list_of_filtered.end())
		  {
		    //print momentum coordinates
		    if(fout)
		      {
			for(int mu=0;mu<NDIM;mu++)
			  if(cmir[mu]<glb_size[mu]/2) master_fprintf(fout,"%d ",cmir[mu]);
			  else                        master_fprintf(fout,"%d ",cmir[mu]-glb_size[mu]);
			master_fprintf(fout,"\n");
		      }
		    
		    //search where data is stored
		    int wrank,iloc;
		    get_loclx_and_rank_of_coord(&iloc,&wrank,cmir); //the remapper will leave holes
		    if(rank==wrank) sl.push_back(std::make_pair(iloc,list_of_filtered.size()*nranks*nso_spi*nso_col+0));
		    
		    list_of_filtered.insert(iglb);
		  }
	      }
	}
    
    //close file if opened
    if(fout) close_file(fout);
    
    //setup
    nfft_filtered=list_of_filtered.size();
    fft_filter_remap.setup_knowing_where_to_send(sl);
  }
  
  //perform fft and store the propagators
  void propagators_fft(int ihit)
  {
    GET_THREAD_ID();
    
    spincolor *qtilde=nissa_malloc("qtilde",loc_vol+bord_vol,spincolor);
    spincolor *qfilt=nissa_malloc("qfilt",nfft_filtered*nso_spi*nso_col,spincolor);
    
    double fft_sign=-1;
    for(size_t iprop=0;iprop<fft_prop_list.size();iprop++)
      {
	const std::string tag=fft_prop_list[iprop];
	master_printf("Fourier transforming propagator %s\n",tag.c_str());
	
	//loop on dirac and color source index
	for(int id_so=0;id_so<nso_spi;id_so++)
	  for(int ic_so=0;ic_so<nso_col;ic_so++)
	    {
	      START_TIMING(fft_time,nfft_tot);
	      
	      //perform fft
	      spincolor *q=Q[tag][so_sp_col_ind(id_so,ic_so)];
	      fft4d((complex*)qtilde,(complex*)q,sizeof(spincolor)/sizeof(complex),fft_sign,1);
	      put_fft_source_phase(qtilde,fft_sign);
	      
	      //gather - check the rewriting pattern above!
	      fft_filter_remap.communicate(qfilt+so_sp_col_ind(id_so,ic_so),qtilde,sizeof(spincolor));
	      
	      STOP_TIMING(fft_time);
	    }
	
	//create filename
	std::string filename=combine("%s/fft_%s",outfolder,tag.c_str());
	if(nhits>1) filename+=combine("_hit_%d",ihit);
	
	//open and write
	FILE *fout=open_file(filename,"w");
	if(rank==0) fwrite(qfilt,sizeof(spincolor)*nso_spi*nso_col,nfft_filtered,fout);
	close_file(fout);
      }
    
    nissa_free(qtilde);
    nissa_free(qfilt);
  }
}
