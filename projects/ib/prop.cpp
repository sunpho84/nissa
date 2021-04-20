#include <nissa.hpp>

#define EXTERN_PROP
 #include "prop.hpp"

#include <set>
#include <tuple>

namespace nissa
{
  //multiply the configuration for an additional u(1) field, defined as exp(-i e q A /3)
  void add_photon_field_to_conf(quad_su3* conf,double charge)
  {
    const double alpha_em=1/137.04;
    const double e2=4*M_PI*alpha_em;
    const double e=sqrt(e2);
    verbosity_lv2_master_printf("Adding backfield, for a quark of charge q=e*%lg/3\n",charge);
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  complex ph;
	  complex_iexp(ph,-e*photon_field[ivol.nastyConvert()][mu][RE]*charge/3.0);
	  safe_su3_prod_complex(conf[ivol.nastyConvert()][mu],conf[ivol.nastyConvert()][mu],ph);
	}
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(conf);
  }
  
  //remove the field
  void rem_photon_field_to_conf(quad_su3 *conf,double q)
  {
    add_photon_field_to_conf(conf,-q);
  }
  
  //get a propagator inverting on "in"
  void get_qprop(spincolor *out,spincolor *in,double kappa,double mass,int r,double charge,double residue,const Momentum& theta)
  {
    
    //rotate the source index - the propagator rotate AS the sign of mass term
    if(twisted_run>0) safe_dirac_prod_spincolor(in,(tau3[r]==-1)?&Pminus:&Pplus,in);
    
    //invert
    START_TIMING(inv_time,ninv_tot);
    
    //get an intenral time
    double tin=take_time();
    
    if(free_theory and charge==0)
      {
	master_printf("   working in FT\n");
	
	tm_quark_info qu(kappa,fabs(mass),r,theta);
	tm_basis_t basis=WILSON_BASE;
	multiply_from_left_by_x_space_twisted_propagator_by_fft(out,in,qu,basis,false);
      }
    else
      {
	quad_su3 *conf=get_updated_conf(charge,theta,glb_conf);
	
	master_printf("   inverting explicitly\n");
	if(clover_run) inv_tmclovD_cg_eoprec(out,NULL,conf,kappa,Cl,invCl,glb_cSW,mass,1000000,residue,in);
	else inv_tmD_cg_eoprec(out,NULL,conf,kappa,mass,1000000,residue,in);
      }
    
    verbosity_lv1_master_printf("Solving time: %lg s\n",take_time()-tin);
    
    STOP_TIMING(inv_time);
    
    //rotate the sink index
    if(twisted_run>0) safe_dirac_prod_spincolor(out,(tau3[r]==-1)?&Pminus:&Pplus,out);
  }
  
  //generate a source, wither a wall or a point in the origin
  void generate_original_source(qprop_t* sou)
  {
    
    //consistency check
    if(!stoch_source and (!diluted_spi_source or !diluted_col_source)) crash("for a non-stochastic source, spin and color must be diluted");
    
    //reset all to begin
    for(int i=0;i<nso_spi*nso_col;i++) vector_reset(sou->sp[i]);
    
    spincolor **sou_proxy=nissa_malloc("sou_proxy",nso_spi*nso_col,spincolor*);
    for(int id_so=0;id_so<nso_spi;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	sou_proxy[so_sp_col_ind(id_so,ic_so)]=sou->sp[so_sp_col_ind(id_so,ic_so)];
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	spincolor c;
	spincolor_put_to_zero(c);
	
	//compute relative coords
	bool is_spat_orig=true;
	GlbCoords rel_c;
	FOR_ALL_SPATIAL_DIRS(mu)
	  {
	    rel_c(mu)=rel_coord_of_loclx(ivol,mu);
	    if(mu()) is_spat_orig&=(rel_c(mu)==0);
	  }
	
	//dilute in space
	int mask=1;
	FOR_ALL_SPATIAL_DIRS(mu) mask&=(rel_c(mu)%diluted_spat_source==0);
	
	//fill colour and spin index 0
	for(int id_si=0;id_si<(diluted_spi_source?1:NDIRAC);id_si++)
	  for(int ic_si=0;ic_si<(diluted_col_source?1:NCOL);ic_si++)
	    {
	      if(stoch_source and mask and (sou->tins==-1 or rel_c(tDir)==sou->tins)) comp_get_rnd(c[id_si][ic_si],&(loc_rnd_gen[ivol.nastyConvert()]),sou->noise_type);
	      if(!stoch_source and is_spat_orig and (sou->tins==-1 or rel_c(tDir)==sou->tins)) complex_put_to_real(c[id_si][ic_si],1);
	    }
	
	//fill other spin indices
	for(int id_so=0;id_so<nso_spi;id_so++)
	  for(int ic_so=0;ic_so<nso_col;ic_so++)
	    for(int id_si=0;id_si<NDIRAC;id_si++)
	      for(int ic_si=0;ic_si<NCOL;ic_si++)
		  if((!diluted_spi_source or (id_so==id_si)) and (!diluted_col_source or (ic_so==ic_si)))
		    complex_copy(sou_proxy[so_sp_col_ind(id_so,ic_so)][ivol.nastyConvert()][id_si][ic_si],c[diluted_spi_source?0:id_si][diluted_col_source?0:ic_si]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    //compute the norm2, set borders invalid
    double ori_source_norm2=0;
    for(int id_so=0;id_so<nso_spi;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	{
	  spincolor *s=sou->sp[so_sp_col_ind(id_so,ic_so)];
	  set_borders_invalid(s);
	  ori_source_norm2+=double_vector_glb_norm2(s,locVol.nastyConvert());
	}
    if(IS_MASTER_THREAD) sou->ori_source_norm2=ori_source_norm2;
    
    nissa_free(sou_proxy);
  }
  
  //insert the photon on the source side
  void insert_external_loc_source(spincolor *out,spin1field *curr,spincolor *in,const GlbCoord& t,const Coords<bool>& dirs)
  {
    if(in==out) crash("in==out");
    
    vector_reset(out);
    
    FOR_ALL_SPATIAL_DIRS(mu)
      if(dirs(mu))
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  if(t==-1 or glbCoordOfLoclx(ivol,tDir)==t)
	    {
	      spincolor temp1,temp2;
	      unsafe_dirac_prod_spincolor(temp1,base_gamma+igamma_of_mu(mu),in[ivol.nastyConvert()]);
	      unsafe_spincolor_prod_complex(temp2,temp1,curr[ivol.nastyConvert()][mu.nastyConvert()]);
	      spincolor_summ_the_prod_idouble(out[ivol.nastyConvert()],temp2,1);
	    }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //insert the external source
  void insert_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *ori,const GlbCoord& t,int r,const Coords<bool>& dirs,int loc)
  {
    if(loc) insert_external_loc_source(out,curr,ori,t,dirs);
    else
      if(twisted_run>0) insert_tm_external_source(out,conf,curr,ori,r,dirs,t);
      else              insert_Wilson_external_source(out,conf,curr,ori,dirs,t);
  }
  
  //insert the tadpole
  void insert_tadpole(spincolor *out,quad_su3 *conf,spincolor *ori,const GlbCoord& t,int r)
  {
    if(twisted_run>0) insert_tm_tadpole(loop_source,conf,ori,r,tadpole,t);
    else              insert_Wilson_tadpole(loop_source,conf,ori,tadpole,t);
  }
  
  //insert the conserved current
  void insert_conserved_current(spincolor *out,quad_su3 *conf,spincolor *ori,const GlbCoord& t,int r,const Coords<bool>& dirs)
  {
    if(twisted_run>0) insert_tm_conserved_current(loop_source,conf,ori,r,dirs,t);
    else              insert_Wilson_conserved_current(loop_source,conf,ori,dirs,t);
  }
  
  //smear the propagator
  template <typename T>
  void smear_prop(spincolor *out,quad_su3 *conf,spincolor *ori,const GlbCoord& t,const T& kappa,int nlevels)
  {
    
    //nb: the smearing radius is given by
    //a*sqrt(2*n*kappa/(1+6*kappa))
    
    START_TIMING(sme_time,nsme_tot);
    
    gaussian_smearing(out,ori,conf,kappa,nlevels);
    select_propagator_timeslice(out,out,t);
    
    STOP_TIMING(sme_time);
  }
  
  //phase the propagator
  void phase_prop(spincolor* out,spincolor* ori,const GlbCoord& t,const Momentum& th)
  {
    FOR_ALL_SPATIAL_DIRS(mu)
      if(fabs((int)(th(mu)/2)-th(mu)/2)>1e-10)
	crash("Error: phase %lg must be an even integer",th(mu));
    
    vector_reset(out);
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	//compute x*p
	double arg=0.0;
	FOR_ALL_SPATIAL_DIRS(mu)
	  arg+=M_PI*th(mu)*rel_coord_of_loclx(ivol,mu)()/glbSize(mu)(); //N.B: valid only if source is on origin...
	
	//compute exp(ip)
	complex factor;
	complex_iexp(factor,arg);
	
	//put the phase
	if(t==-1 or glbCoordOfLoclx(ivol,tDir)==t)
	  unsafe_spincolor_prod_complex(out[ivol.nastyConvert()],ori[ivol.nastyConvert()],factor);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //backward flow the propagator
  void back_flow_prop(spincolor* out,quad_su3* conf,spincolor* ori,const GlbCoord& t,double dt,int nflows)
  {
    
    START_TIMING(bflw_time,nbflw_tot);
    
    //the flown conf
    quad_su3 *flown_conf=nissa_malloc("flown_conf",locVolWithBordAndEdge.nastyConvert(),quad_su3);
    vector_copy(flown_conf,conf);
    
    //the recursive flower, need to cache backward integration
    Wflow_pars_t Wf;
    Wf.nflows=nflows;
    Wf.dt=dt;
    recursive_Wflower_t recu(Wf,flown_conf);
    
    //the adjoint flower needed for fermionic source
    fermion_adjoint_flower_t<spincolor> adj_ferm_flower(dt,all_other_dirs[0]);
    
    //at each step it goes from iflow+1 to iflow
    select_propagator_timeslice(out,ori,t);
    for(int iflow=nflows-1;iflow>=0;iflow--)
      {
	//update conf to iflow
	double t=dt*iflow;
	master_printf(" flow back to %d/%d, t %lg\n",iflow,nflows,t);
	recu.update(iflow);
	
	//make the flower generate the intermediate step between iflow and iflow+1
	adj_ferm_flower.generate_intermediate_steps(flown_conf);
	
	adj_ferm_flower.flow_fermion(out);
      }
    
    STOP_TIMING(bflw_time);
    
    nissa_free(flown_conf);
  }
  
  //flow the propagator
  void flow_prop(spincolor* out,quad_su3* conf,spincolor* ori,const GlbCoord& t,double dt,int nflows)
  {
    
    START_TIMING(flw_time,nflw_tot);
    
    //the flown conf
    quad_su3 *flown_conf=nissa_malloc("flown_conf",locVolWithBordAndEdge.nastyConvert(),quad_su3);
    vector_copy(flown_conf,conf);
    
    //the flower, need to cache integration
    fermion_flower_t<spincolor,4> ferm_flower(dt,all_other_dirs[0]);
    
    select_propagator_timeslice(out,ori,t);
    for(int iflow=0;iflow<=nflows;iflow++)
      {
	//update conf to iflow
	double t=dt*iflow;
	master_printf(" flow forward to %d/%d, t %lg, initial plaquette: %.16lg\n",iflow,nflows,t,global_plaquette_lx_conf(flown_conf));
	
	//make the flower generate the intermediate step between iflow-1 and iflow
	ferm_flower.generate_intermediate_steps(flown_conf);
	ferm_flower.flow_fermion(out);
	ferm_flower.prepare_for_next_flow(flown_conf);
      }
    
    STOP_TIMING(flw_time);
    
    nissa_free(flown_conf);
  }
  
  void build_source(spincolor* out,std::vector<source_term_t>* source_terms,int isou)
  {
    
    vector_reset(out);
    
    for(auto& c : *source_terms)
      {
	complex coef={c.second.first,c.second.second};
	spincolor *p=Q[c.first][isou];
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  spincolor_summ_the_prod_complex(out[ivol.nastyConvert()],p[ivol.nastyConvert()],coef);
	NISSA_PARALLEL_LOOP_END;
      }
    set_borders_invalid(out);
  }
  
  //generate a sequential source
  void generate_source(insertion_t inser,char *ext_field_path,int r,double charge,double kappa,const Momentum& kappa_asymm,const Momentum& theta,std::vector<source_term_t>& source_terms,int isou,const GlbCoord& t)
  {
    source_time-=take_time();
    
    GlbCoord rel_t=t;
    if(rel_t!=-1) rel_t=(t+source_coord(tDir))%glbSize(tDir);
    
    quad_su3 *conf;
    if(not is_smearing_ins(inser)) conf=get_updated_conf(charge,theta,glb_conf);
    else
      {
	quad_su3 *ext_conf;
	if(ape_smeared_conf) ext_conf=ape_smeared_conf;
	else                 ext_conf=glb_conf;
	conf=get_updated_conf(0.0,theta,ext_conf);
      }
    
    spincolor* ori=nissa_malloc("ori",locVolWithBord.nastyConvert(),spincolor);
    build_source(ori,&source_terms,isou);
    
    spin1field *ext_field=nullptr;
    if(inser==EXT_FIELD)
      {
	ext_field=nissa_malloc("ext_field",locVolWithBord.nastyConvert(),spin1field);
	read_real_vector(ext_field,combine("%s/%s",outfolder,ext_field_path),"Current");
      }
    
    master_printf("Inserting r: %d\n",r);
    switch(inser)
      {
      case PROP:prop_multiply_with_gamma(loop_source,0,ori,rel_t);break;
      case SCALAR:prop_multiply_with_gamma(loop_source,0,ori,rel_t);break;
      case PSEUDO:prop_multiply_with_gamma(loop_source,5,ori,rel_t);break;
      case GAMMA:prop_multiply_with_gamma(loop_source,r,ori,rel_t);break;
      case PHOTON:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,all_dirs,loc_hadr_curr);break;
      case PHOTON0:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,only_dir[0],loc_hadr_curr);break;
      case PHOTON1:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,only_dir[1],loc_hadr_curr);break;
      case PHOTON2:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,only_dir[2],loc_hadr_curr);break;
      case PHOTON3:insert_external_source(loop_source,conf,photon_field,ori,rel_t,r,only_dir[3],loc_hadr_curr);break;
      case PHOTON_PHI:insert_external_source(loop_source,conf,photon_phi,ori,rel_t,r,all_dirs,loc_hadr_curr);break;
      case PHOTON_ETA:insert_external_source(loop_source,conf,photon_eta,ori,rel_t,r,all_dirs,loc_hadr_curr);break;
      case TADPOLE:insert_tadpole(loop_source,conf,ori,rel_t,r);break;
      case CVEC:insert_conserved_current(loop_source,conf,ori,rel_t,r,all_dirs);break;
      case CVEC0:insert_conserved_current(loop_source,conf,ori,rel_t,r,only_dir[0]);break;
      case CVEC1:insert_conserved_current(loop_source,conf,ori,rel_t,r,only_dir[1]);break;
      case CVEC2:insert_conserved_current(loop_source,conf,ori,rel_t,r,only_dir[2]);break;
      case CVEC3:insert_conserved_current(loop_source,conf,ori,rel_t,r,only_dir[3]);break;
      case EXT_FIELD:insert_external_source(loop_source,conf,ext_field,ori,rel_t,r,all_dirs,loc_hadr_curr);break;
      case SMEARING:smear_prop(loop_source,conf,ori,rel_t,kappa,r);break;
      case ANYSM:smear_prop(loop_source,conf,ori,rel_t,kappa_asymm,r);break;
      case WFLOW:flow_prop(loop_source,conf,ori,rel_t,kappa,r);break;
      case BACK_WFLOW:back_flow_prop(loop_source,conf,ori,rel_t,kappa,r);break;
      case PHASING:phase_prop(loop_source,ori,rel_t,theta);break;
      }
    
    nissa_free(ori);
    if(ext_field)
      nissa_free(ext_field);
    
    source_time+=take_time();
    nsource_tot++;
  }
  
  //Generate all the original sources
  void generate_original_sources(int ihit,bool skip_io)
  {
    
    for(size_t i=0;i<ori_source_name_list.size();i++)
      {
	std::string &name=ori_source_name_list[i];
	master_printf("Generating source \"%s\"\n",name.c_str());
	qprop_t *q=&Q[name];
	generate_original_source(q);
	
	if(not skip_io)
	  for(int id_so=0;id_so<nso_spi;id_so++)
	    for(int ic_so=0;ic_so<nso_col;ic_so++)
	      {
		//combine the filename
		std::string path=combine("%s/hit%d_source%s_idso%d_icso%d",outfolder,ihit,name.c_str(),id_so,ic_so);
		
		int isou=so_sp_col_ind(id_so,ic_so);
		spincolor *sou=(*q)[isou];
		
		//if the prop exists read it
		if(file_exists(path))
		  {
		    master_printf("  loading the source dirac index %d, color %d\n",id_so,ic_so);
		    START_TIMING(read_prop_time,nread_prop);
		    read_real_vector(sou,path,"scidac-binary-data");
		    STOP_TIMING(read_prop_time);
		  }
		else
		  {
		    master_printf("  file %s not available, skipping loading\n",path.c_str());
		    
		    //and store if needed
		    if(q->store)
		      {
			master_printf("  writing the source dirac index %d, color %d\n",id_so,ic_so);
			START_TIMING(store_prop_time,nstore_prop);
			write_real_vector(path,sou,64,"scidac-binary-data");
			STOP_TIMING(store_prop_time);
		      }
		  }
	      }
      }
  }
  
  //generate all the quark propagators
  void generate_quark_propagators(int ihit)
  {
    for(size_t i=0;i<qprop_name_list.size();i++)
      {
	//get names
	std::string name=qprop_name_list[i];
	qprop_t &q=Q[name];
	
	//get ori_source norm2
	const std::string& first_source=q.source_terms.front().first;
	const double ori_source_norm2=q.ori_source_norm2=Q[first_source].ori_source_norm2;
	for(auto& n : q.source_terms)
	  {
	    double this_source_norm2=Q[n.first].ori_source_norm2;
	    if(ori_source_norm2!=this_source_norm2)
	      crash("first source %s has different norm2 %lg than %s, %lg",first_source.c_str(),ori_source_norm2,n.first.c_str(),this_source_norm2);
	  }
	
	//write info on mass and r
	if(twisted_run) master_printf(" mass[%d]=%lg, r=%d, theta={%lg,%lg,%lg}\n",i,q.mass,q.r,q.theta(xDir),q.theta(yDir),q.theta(zDir));
	else            master_printf(" kappa[%d]=%lg, theta={%lg,%lg,%lg}\n",i,q.kappa,q.theta(xDir),q.theta(yDir),q.theta(zDir));
	
	//compute the inverse clover term, if needed
	if(clover_run) invert_twisted_clover_term(invCl,q.mass,q.kappa,Cl);
	
	//create the description of the source
	std::string source_descr;
	if(q.source_terms.size()==1)
	  source_descr=first_source;
	else
	  {
	    source_descr="(";
	    for(int i=0;i<(int)q.source_terms.size();i++)
	      {
		source_term_t& this_source=q.source_terms[i];
		complex c={this_source.second.first,this_source.second.second};
		if(i>0) source_descr+="+";
		source_descr+=this_source.first+"*("+std::to_string(c[RE])+","+std::to_string(c[IM])+")";
	      }
	    source_descr+=")";
	  }
	
	insertion_t insertion=q.insertion;
	master_printf("Generating propagator %s inserting %s on source %s\n",name.c_str(),ins_name[insertion],source_descr.c_str());
	for(int id_so=0;id_so<nso_spi;id_so++)
	  for(int ic_so=0;ic_so<nso_col;ic_so++)
	    {
	      int isou=so_sp_col_ind(id_so,ic_so);
	      generate_source(insertion,q.ext_field_path,q.r,q.charge,q.kappa,q.kappa_asymm,q.theta,q.source_terms,isou,q.tins);
	      spincolor *sol=q[isou];
	      
	      //combine the filename
	      std::string path=combine("%s/hit%d_prop%s_idso%d_icso%d",outfolder,ihit,name.c_str(),id_so,ic_so);
	      
	      //if the prop exists read it
	      if(file_exists(path))
		{
		  master_printf("  loading the solution, dirac index %d, color %d\n",id_so,ic_so);
		  START_TIMING(read_prop_time,nread_prop);
		  read_real_vector(sol,path,"scidac-binary-data");
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
		      write_real_vector(path,sol,64,"scidac-binary-data");
		      STOP_TIMING(store_prop_time);
		    }
		  master_printf("  finished the calculation of dirac index %d, color %d\n",id_so,ic_so);
		}
	    }
      }
  }
  
  /////////////////////////////////////////////// photon propagators ///////////////////////////////////////////
  
  //allocate the photon fields
  void allocate_photon_fields()
  {
    photon_eta=nissa_malloc("photon_eta",locVolWithBord.nastyConvert(),spin1field);
    photon_field=nissa_malloc("photon_field",locVolWithBord.nastyConvert(),spin1field);
    photon_phi=nissa_malloc("photon_phi",locVolWithBord.nastyConvert(),spin1field);
  }
  
  //free the photon fields
  void free_photon_fields()
  {
    nissa_free(photon_eta);
    nissa_free(photon_phi);
    nissa_free(photon_field);
  }
  
  //wrapper to generate a stochastic propagator
  void generate_photon_stochastic_propagator(int ihit)
  {
    
    photon_prop_time-=take_time();
    
    //generate source and stochastich propagator
    master_printf("Generating photon stochastic propagator\n");
    generate_stochastic_tlSym_gauge_propagator(photon_phi,photon_eta,photon);
    
    //generate the photon field
    multiply_by_sqrt_tlSym_gauge_propagator(photon_field,photon_eta,photon);
    
    std::vector<std::pair<std::string,spin1field*> > name_field;
    name_field.push_back(std::make_pair("eta",photon_eta));
    name_field.push_back(std::make_pair("phi",photon_phi));
    name_field.push_back(std::make_pair("A",photon_field));
    for(std::vector<std::pair<std::string,spin1field*> >::iterator nf=name_field.begin();nf!=name_field.end();nf++)
      {
	std::string &name=(*nf).first;
	spin1field *ph=(*nf).second;
	
	//combine the filename
	std::string path=combine("%s/hit%d_field%s",outfolder,ihit,name.c_str());
	
	//if asked and the file exists read it
	if(load_photons)
	  {
	    if(file_exists(path))
	      {
		master_printf("  loading the photon field %s\n",name.c_str());
		START_TIMING(read_prop_time,nread_prop);
		read_real_vector(ph,path,"scidac-binary-data");
		STOP_TIMING(read_prop_time);
	      }
	    else master_printf("  file %s not available, skipping loading\n",path.c_str());
	  }
	
	//if asked, write it
	if(store_photons)
	  {
	    master_printf("  storing the photon field %s\n",name.c_str());
	    START_TIMING(store_prop_time,nstore_prop);
	    write_real_vector(path,ph,64,"scidac-binary-data");
	    STOP_TIMING(store_prop_time);
	  }
      }
    
    //invalidate internal conf
    inner_conf_valid=false;
    
    photon_prop_time+=take_time();
    nphoton_prop_tot++;
  }
  
  //put the phase of the source due to missing e(iky)
  void put_fft_source_phase(spincolor* qtilde,double fft_sign)
  {
    
    NISSA_PARALLEL_LOOP(imom,0,locVol)
      {
	//sum_mu -sign*2*pi*p_mu*y_mu/L_mu
	double arg=0;
	FOR_ALL_SPATIAL_DIRS(mu)
	  arg+=-fft_sign*2*M_PI*glbCoordOfLoclx(imom,mu)()*source_coord(mu)()/glbSize(mu)();
	
	complex f={cos(arg),sin(arg)};
	spincolor_prodassign_complex(qtilde[imom.nastyConvert()],f);
	
	// spincolor_put_to_zero(qtilde[imom]);
	// for(int mu=0;mu<4;mu++) qtilde[imom](mu)(timeDirection)(timeDirection)=glb_coord_of_loclx[imom](mu);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(qtilde);
  }
  
  //initialize the fft filter, once forever
  void init_fft_filter_from_range(std::vector<std::pair<fft_mom_range_t,double>>& fft_mom_range_list)
  {
    master_printf("Initializing fft filter\n");
    
    //file where to store output
    FILE *fout=NULL;
    const char path_list[]="mom_list.txt";
    if(not file_exists("mom_list.txt")) fout=open_file(path_list,"w");
    
    //store the list of filtered
    std::set<GlbLxSite> list_of_filtered;;
    
    //scattering list
    all_to_all_scattering_list_t sl;
    for(auto &f : fft_mom_range_list)
      for(GlbLxSite vol=vol_of_lx(f.first.width),ifilt=0;ifilt<vol;ifilt++)
	{
	  //gets the coordinate in the filtering volume
	  GlbCoords c;
	  coord_of_lx(c,ifilt,f.first.width);
	  coord_summassign(c,f.first.offs,glbSize);
	  
	  //compute p~4/p~2^2
	  double pt2=0,pt4=0;
	  FOR_ALL_SPATIAL_DIRS(mu)
	    {
	      double pmu=M_PI*(2*c(mu)()+(mu==0)*temporal_bc)/glbSize(mu)();
	      double ptmu=sin(pmu);
	      pt2+=sqr(ptmu);
	      pt4+=pow(ptmu,4.0);
	    }
	  
	  double p4_fr_p22_max=f.second;
	  if(pt4/sqr(pt2)<p4_fr_p22_max)
	    for(int imir=0;imir<pow(2,NDIM);imir++)
	      {
		//get mirrorized
		GlbCoords cmir;
		FOR_ALL_SPATIAL_DIRS(mu)
		  cmir(mu)=get_mirrorized_site_coord(c(mu)+(mu==0 and get_bit(imir,0) and temporal_bc==ANTIPERIODIC_BC),mu,get_bit(imir,mu()));
		
		//check if not already collected
		const GlbLxSite iglb=glblx_of_coord(cmir);
		if(list_of_filtered.find(iglb)==list_of_filtered.end())
		  {
		    //print momentum coordinates
		    if(fout)
		      {
			FOR_ALL_SPATIAL_DIRS(mu)
			  if(cmir(mu)<glbSize(mu)/2) master_fprintf(fout,"%d ",cmir(mu));
			  else                        master_fprintf(fout,"%d ",cmir(mu)-glbSize(mu));
			master_fprintf(fout,"\n");
		      }
		    
		    //search where data is stored
		    Rank wrank;
		    LocLxSite iloc;
		    get_loclx_and_rank_of_coord(iloc,wrank,cmir); //the remapper will leave holes
		    if(rank==wrank)
		      sl.push_back(std::make_pair(iloc(),list_of_filtered.size()*nranks));
		    
		    list_of_filtered.insert(iglb);
		  }
	      }
	}
    
    //close file if opened
    if(fout) close_file(fout);
    
    //setup
    fft_filterer.emplace_back((int)list_of_filtered.size(),sl,std::string(""));
  }
  
  ///helper to do on master
  template <typename F,
	    typename...Args>
  auto do_on_master(F f,Args&&...args)
  {
    using R=decltype(f(std::forward<Args>(args)...));
    
    R res;
    memset((void*)&res,0,sizeof(R));
    
    if(rank==0)
      res=f(std::forward<Args>(args)...);
    MPI_Bcast(&res,sizeof(R),MPI_CHAR,0,MPI_COMM_WORLD);
    
    return res;
  }
  
  //read the list of momenta to write
  void init_fft_filterer_from_file(const char *fileout_suff,const char *filein_name)
  {
    //file where to read the list
    FILE *fin=open_file(filein_name,"r");
    auto read_int=[fin]()
      {
	int res=0;
	bool eof;
	
	eof=(fscanf(fin,"%d",&res)!=1);
	
	if(eof) res=0;
	
	return std::make_tuple(res,eof);
      };
    
    std::set<GlbLxSite> list_of_filtered;;
    
    all_to_all_scattering_list_t sl;
    bool ended=false;
    while(not ended)
      {
	GlbCoords c;
	FOR_ALL_SPATIAL_DIRS(mu)
	  {
	    int cmu;
	    bool emu;
	    
	    std::tie(cmu,emu)=do_on_master(read_int);
	    
	    c(mu)=(cmu+glbSize(mu))%glbSize(mu);
	    
	    ended|=emu;
	  }
	
	if(not ended)
	  {
	    //search where data is stored
	    Rank wrank;
	    LocLxSite iloc;
	    get_loclx_and_rank_of_coord(iloc,wrank,c);
	    if(rank==wrank)
	      sl.push_back(std::make_pair(iloc(),list_of_filtered.size()*nranks));
	    
	    const GlbLxSite iglb=glblx_of_coord(c);
	    list_of_filtered.insert(iglb);
	  }
      }
    
    //close file if opened
    close_file(fin);
    
    fft_filterer.emplace_back(list_of_filtered.size(),sl,fileout_suff);
  }
  
  //perform fft and store the propagators
  void propagators_fft(int ihit)
  {
    
    spincolor *qtilde=nissa_malloc("qtilde",locVolWithBord.nastyConvert(),spincolor);
    
    int nf=fft_filterer.size();
    spincolor *qfilt[nf];
    spincolor *qfilt_temp[nf];
    for(int i=0;i<nf;i++)
      {
	qfilt[i]=nissa_malloc("qfilt",fft_filterer[i].nfft_filtered*nso_spi*nso_col,spincolor);
	qfilt_temp[i]=nissa_malloc("qfilt_temp",fft_filterer[i].nfft_filtered,spincolor);
      }
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
	      for(int i=0;i<nf;i++)
		{
		  master_printf("Filtering %d/%d\n",i,nf);
		  fft_filterer[i].fft_filter_remap.communicate(qfilt_temp[i],qtilde,sizeof(spincolor));
		  crash("#warning reimplement");
		  // NISSA_PARALLEL_LOOP(imom,0,fft_filterer[i].nfft_filtered)
		  //   spincolor_copy(qfilt[i][imom*nso_spi*nso_col+so_sp_col_ind(id_so,ic_so)],qfilt_temp[i][imom]);
		  // NISSA_PARALLEL_LOOP_END;
		  set_borders_invalid(qfilt[i]);
		}
	      
	      STOP_TIMING(fft_time);
	    }
	
	//create filename
	for(int i=0;i<nf;i++)
	  {
	    std::string filename=combine("%s/fft_%s",outfolder,tag.c_str());
	    if(fft_filterer[i].file_suffix!="") filename+=combine("_%s",fft_filterer[i].file_suffix.c_str());
	    if(nhits>1) filename+=combine("_hit_%d",ihit);
	    
	    //open and write
	    FILE *fout=open_file(filename,"w");
	    if(rank==0) fwrite(qfilt[i],sizeof(spincolor)*nso_spi*nso_col,fft_filterer[i].nfft_filtered,fout);
	    close_file(fout);
	  }
      }
    
    for(int i=0;i<nf;i++)
      {
	nissa_free(qfilt[i]);
	nissa_free(qfilt_temp[i]);
      }
    
    nissa_free(qtilde);
  }
}
