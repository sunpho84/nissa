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
  void get_qprop(spincolor *out,spincolor *in,double kappa,double mass,int r,double charge,double residue,double *theta)
  {
    GET_THREAD_ID();
    
    //rotate the source index - the propagator rotate AS the sign of mass term
    if(twisted_run) safe_dirac_prod_spincolor(in,(tau3[r]==-1)?&Pminus:&Pplus,in);
    
    //invert
    START_TIMING(inv_time,ninv_tot);
    
    quad_su3 *conf=get_updated_conf(charge,theta,glb_conf);
    
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
  THREADABLE_FUNCTION_4ARG(phase_prop, spincolor*,out, spincolor*,ori, int,t, double*,th)
  {
    GET_THREAD_ID();
    
    for(int mu=1;mu<NDIM;mu++) if(fabs((int)(th[mu]/2)-th[mu]/2)>1e-10) crash("Error: phase %lg must be an even integer",th);
    
    vector_reset(out);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	//compute x*p
	double arg=0.0;
	for(int mu=1;mu<NDIM;mu++) arg+=M_PI*th[mu]*rel_coord_of_loclx(ivol,mu)/glb_size[mu]; //N.B: valid only if source is on origin...
	
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
  void generate_source(insertion_t inser,int r,double charge,double kappa,double *theta,spincolor *ori,int t)
  {
    source_time-=take_time();
    
    int rel_t=t;
    if(rel_t!=-1) rel_t=(t+source_coord[0])%glb_size[0];
    coords dirs[NDIM]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    
    quad_su3 *conf;
    if(inser!=SMEARING) conf=get_updated_conf(charge,theta,glb_conf);
    else
      {
	quad_su3 *ext_conf;
	if(ape_smeared_conf) ext_conf=ape_smeared_conf;
	else                 ext_conf=glb_conf;
	conf=get_updated_conf(0.0,theta,ext_conf);
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
  
  //Generate all the original sources
  void generate_original_sources(int ihit)
  {
    GET_THREAD_ID();
    
    for(size_t i=0;i<ori_source_name_list.size();i++)
      {
	std::string &name=ori_source_name_list[i];
	master_printf("Generating source \"%s\"\n",name.c_str());
	qprop_t *q=&Q[name];
	generate_original_source(q);
	
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
	      else master_printf("  file %s not available, skipping loading\n",path.c_str());
	      
	      //and store if needed
	      if(q->store)
		{
		  master_printf("  writing the source dirac index %d, color %d\n",id_so,ic_so);
		  START_TIMING(store_prop_time,nstore_prop);
		  write_double_vector(path,sou,64,"scidac-binary-data");
		  STOP_TIMING(store_prop_time);
		}
	    }
      }
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
	if(twisted_run) master_printf(" mass[%d]=%lg, r=%d, theta={%lg,%lg,%lg}\n",i,q.mass,q.r,q.theta[1],q.theta[2],q.theta[3]);
	else            master_printf(" kappa[%d]=%lg, theta={%lg,%lg,%lg}\n",i,q.kappa,q.theta[1],q.theta[2],q.theta[3]);
	
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
		      write_double_vector(path,sol,64,"scidac-binary-data");
		      STOP_TIMING(store_prop_time);
		    }
		  master_printf("  finished the inversion dirac index %d, color %d\n",id_so,ic_so);
		}
	    }
      }
  }
  
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
