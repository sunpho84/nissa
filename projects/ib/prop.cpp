#include <nissa.hpp>

#define EXTERN_PROP
# include "prop.hpp"

#include <tuple>

namespace nissa
{
  void add_photon_field_to_conf(LxField<quad_su3>& conf,
				const double& charge)
  {
    const double alpha_em=1/137.04;
    const double e2=4*M_PI*alpha_em;
    const double e=sqrt(e2);
    
    VERBOSITY_LV2_MASTER_PRINTF("Adding backfield, for a quark of charge q=e*%lg/3\n",charge);
    
    PAR(0,locVol,
	CAPTURE(charge,e,
		TO_WRITE(conf)),
	ivol,
	{
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      complex ph;
	      complex_iexp(ph,-e*(*photon_field)[ivol][mu][RE]*charge/3.0);
	      safe_su3_prod_complex(conf[ivol][mu],conf[ivol][mu],ph);
	    }
	});
  }
  
  /// Remove the field
  void rem_photon_field_to_conf(LxField<quad_su3>& conf,
				const double& charge)
  {
    add_photon_field_to_conf(conf,-charge);
  }
  
  //get a propagator inverting on "in"
  void get_qprop(LxField<spincolor>& out,
		 LxField<spincolor>& in,
		 const double& kappa,
		 const double& mass,
		 const int& r,
		 const double charge,
		 const double& residue,
		 const Momentum& theta)
  {
    //rotate the source index - the propagator rotate AS the sign of mass term
    if(twisted_run>0)
      PAR(0,locVol,
	  CAPTURE(r,
		  TO_WRITE(in)),
	  ivol,
	  {
	    auto s=in[ivol];
	    safe_dirac_prod_spincolor(s,(tau3[r]==-1)?Pminus:Pplus,s);
	  });
    
    //invert
    START_TIMING(inv_time,ninv_tot);
    
    //get an intenral time
    double tin=take_time();
    
    if(free_theory and charge==0)
      {
	MASTER_PRINTF("   working in FT\n");
	
	const TmQuarkInfo qu(kappa,fabs(mass),r,theta);
	const tm_basis_t basis=WILSON_BASE;
	multiply_from_left_by_x_space_twisted_propagator_by_fft(out,in,qu,basis,false);
      }
    else
      {
	MASTER_PRINTF("   inverting explicitly\n");
	
#ifdef USE_EXTERNAL_SOLVER
	
	std::ostringstream os;
	os.precision(16);
	os<<conf_path<<"ch"<<charge;
	for(int mu=0;mu<NDIM;mu++)
	  os<<"th["<<mu<<"]"<<theta[mu];
	export_conf::confTag=os.str();
	export_conf::relyOnTag=true;
#endif
	
	const LxField<quad_su3>* conf=get_updated_conf(charge,theta,*glb_conf);
	
	if(clover_run)
	  inv_tmclovD_cg_eoprec(out,std::nullopt,*conf,kappa,*Cl,invCl,glb_cSW,mass,1000000,residue,in);
	else
	  inv_tmD_cg_eoprec(out,std::nullopt,*conf,kappa,mass,1000000,residue,in);
      }
    
    VERBOSITY_LV1_MASTER_PRINTF("Solving time: %lg s\n",take_time()-tin);
    
    STOP_TIMING(inv_time);
    
    //rotate the sink index
    PAR(0,locVol,
	CAPTURE(r,
		TO_WRITE(out)),
	ivol,
	{
	  auto s=out[ivol];
	  safe_dirac_prod_spincolor(s,(tau3[r]==-1)?Pminus:Pplus,s);
	});
  }
  
  //insert the photon on the source side
  void insert_external_loc_source(LxField<spincolor>& out,
				  const LxField<spin1field>& curr,
				  const LxField<spincolor>& in,
				  const int& t,
				  const WhichDirs& dirs)
  {
    
    if(in==out) CRASH("in==out");
    
    out.reset();
    
    for(int mu=0;mu<NDIM;mu++)
      if(dirs[mu])
	PAR(0,locVol,
	    CAPTURE(t,mu,
		    TO_WRITE(out),
		    TO_READ(in),
		    TO_READ(curr)),
	    ivol,
	    {
	      if(t==-1 or glbCoordOfLoclx[ivol][0]==t)
		{
		  spincolor temp1,temp2;
		  unsafe_dirac_prod_spincolor(temp1,base_gamma[iGammaOfMu(mu)],in[ivol]);
		  unsafe_spincolor_prod_complex(temp2,temp1,curr[ivol][mu]);
		  spincolor_summ_the_prod_idouble(out[ivol],temp2,1);
		}
	    });
  }
  
  //insert a field in the current
  template <typename F>
  void insert_external_loc_source(LxField<spincolor>& out,
				  F&& currCalc,
				  const LxField<spincolor>& in,
				  const int& t)
  {
    if(in==out) CRASH("in==out");
    
    out.reset();
    
    for(int mu=0;mu<NDIM;mu++)
      PAR(0,locVol,
	  CAPTURE(t,mu,currCalc,
		  TO_WRITE(out),
		  TO_READ(in)),
	  ivol,
	  {
	    if(t==-1 or glbCoordOfLoclx[ivol][0]==t)
	      {
		spincolor temp1;
		unsafe_dirac_prod_spincolor(temp1,base_gamma[iGammaOfMu(mu)],in[ivol]);
		
		complex curr;
		currCalc(curr,ivol,mu,0.0);
		
		spincolor temp2;
		unsafe_spincolor_prod_complex(temp2,temp1,curr);
		spincolor_summ_the_prod_idouble(out[ivol],temp2,1);
	      }
	  });
  }
  
  //insert the external source
  void insert_external_source(LxField<spincolor>& out,
			      const LxField<quad_su3>& conf,
			      const LxField<spin1field>& curr,
			      const LxField<spincolor>& ori,
			      const int& t,
			      const int& r,
			      const WhichDirs& dirs,
			      const int loc)
  {
    if(loc) insert_external_loc_source(out,curr,ori,t,dirs);
    else
      if(twisted_run>0) insert_tm_external_source(out,conf,curr,ori,r,dirs,t);
      else              insert_Wilson_external_source(out,conf,curr,ori,dirs,t);
  }
  
  //insert the external source
  template <typename F>
  void insert_external_source(LxField<spincolor>& out,
			      const LxField<quad_su3>& conf,
			      F&& currCalc,
			      const LxField<spincolor>& ori,
			      const int& t,
			      const int& r,
			      const int& loc)
  {
    if(loc) insert_external_loc_source(out,currCalc,ori,t);
    else
      if(twisted_run>0) insert_tm_external_source(out,conf,currCalc,ori,r,t);
      else              insert_Wilson_external_source(out,conf,currCalc,ori,t);
  }
  
  //insert the tadpole
  void insert_tadpole(LxField<spincolor>& out,
		      const LxField<quad_su3>& conf,
		      const LxField<spincolor>& ori,
		      const int& t,
		      const int& r)
  {
    if(twisted_run>0) insert_tm_tadpole(*loop_source,conf,ori,r,tadpole,t);
    else              insert_Wilson_tadpole(*loop_source,conf,ori,tadpole,t);
  }
  
  /// Insert the lepton loop
  LxField<spin1field> get_lepton_loop(const double& mass,
				      const Momentum& theta)
  {
    LxField<spin1field> lepton_loop("leptonLoop",WITH_HALO);
    
    return lepton_loop;
  }
  
  //insert the conserved current
  void insert_conserved_current(LxField<spincolor>& out,
				const LxField<quad_su3>& conf,
				const LxField<spincolor>& ori,
				const int& t,
				const int& r,
				const WhichDirs& dirs)
  {
    if(twisted_run>0) insert_tm_conserved_current(*loop_source,conf,ori,r,dirs,t);
    else              insert_Wilson_conserved_current(*loop_source,conf,ori,dirs,t);
  }
  
  //smear the propagator
  template <typename T>
  void smear_prop(LxField<spincolor>& out,
		  const LxField<quad_su3>& conf,
		  const LxField<spincolor>& ori,
		  const int& t,
		  const T& kappa,
		  const int& nlevels)
  {
    
    //nb: the smearing radius is given by
    //a*sqrt(2*n*kappa/(1+6*kappa))
    
    START_TIMING(sme_time,nsme_tot);
    
    gaussian_smearing(out,ori,conf,kappa,nlevels);
    select_propagator_timeslice(out,out,t);
    
    STOP_TIMING(sme_time);
  }
  
  //phase the propagator
  void phase_prop(LxField<spincolor>& out,
		  const LxField<spincolor>& ori,
		  const int& t,
		  const Momentum& th)
  {
    
    for(int mu=1;mu<NDIM;mu++)
      if(fabs((int)(th[mu]/2)-th[mu]/2)>1e-10)
	CRASH("Error: phase %lg must be an even integer",th[mu]);
    
    out.reset();
    PAR(0,locVol,
	CAPTURE(t,th,
		TO_WRITE(out),
		TO_READ(ori)),
	ivol,
	{
	  //compute x*p
	  double arg=0.0;
	  for(int mu=1;mu<NDIM;mu++)
	    arg+=M_PI*th[mu]*rel_coord_of_loclx(ivol,mu)/glbSize[mu]; //N.B: valid only if source is on origin...
	  
	  //compute exp(ip)
	  complex factor;
	  complex_iexp(factor,arg);
	  
	  //put the phase
	  if(t==-1 or glbCoordOfLoclx[ivol][0]==t)
	    unsafe_spincolor_prod_complex(out[ivol],ori[ivol],factor);
	});
  }
  
  //backward flow the propagator
  void back_flow_prop(LxField<spincolor>& out,
		      const LxField<quad_su3>& conf,
		      const LxField<spincolor>& ori,
		      const int& t,
		      const double& dt,
		      const int& nflows)
  {
    
    START_TIMING(bflw_time,nbflw_tot);
    
    //the flown conf
    LxField<quad_su3> flown_conf("flown_conf",WITH_HALO_EDGES);
    flown_conf=conf;
    
    //the recursive flower, need to cache backward integration
    Wflow_pars_t Wf;
    Wf.nflows=nflows;
    Wf.dt=dt;
    recursive_Wflower_t recu(Wf,flown_conf);
    
    //the adjoint flower needed for fermionic source
    fermion_adjoint_flower_t<spincolor> adj_ferm_flower(dt,allOtherDirs[0]);
    
    //at each step it goes from iflow+1 to iflow
    select_propagator_timeslice(out,ori,t);
    for(int iflow=nflows-1;iflow>=0;iflow--)
      {
	//update conf to iflow
	double t=dt*iflow;
	MASTER_PRINTF(" flow back to %d/%d, t %lg\n",iflow,nflows,t);
	recu.update(iflow);
	
	//make the flower generate the intermediate step between iflow and iflow+1
	adj_ferm_flower.generate_intermediate_steps(flown_conf);
	
	adj_ferm_flower.flow_fermion(out);
      }
    
    STOP_TIMING(bflw_time);
  }
  
  //flow the propagator
  void flow_prop(LxField<spincolor>& out,
		 const LxField<quad_su3>& conf,
		 const LxField<spincolor>& ori,
		 const int& t,
		 const double& dt,
		 const int& nflows)
  {
    
    START_TIMING(flw_time,nflw_tot);
    
    //the flown conf
    LxField<quad_su3> flown_conf("flown_conf",WITH_HALO_EDGES);
    flown_conf=conf;
    
    //the flower, need to cache integration
    fermion_flower_t<spincolor,4> ferm_flower(dt,allOtherDirs[0]);
    
    select_propagator_timeslice(out,ori,t);
    for(int iflow=0;iflow<=nflows;iflow++)
      {
	//update conf to iflow
	double t=dt*iflow;
	MASTER_PRINTF(" flow forward to %d/%d, t %lg, initial plaquette: %.16lg\n",iflow,nflows,t,global_plaquette_lx_conf(flown_conf));
	
	//make the flower generate the intermediate step between iflow-1 and iflow
	ferm_flower.generate_intermediate_steps(flown_conf);
	ferm_flower.flow_fermion(out);
	ferm_flower.prepare_for_next_flow(flown_conf);
      }
    
    STOP_TIMING(flw_time);
  }
  
  void build_source(LxField<spincolor>& out,
		    const std::vector<source_term_t>& source_terms,
		    const int& isou)
  {
    MASTER_PRINTF("Creating the source\n");
    
    out.reset();
    
    for(auto& c : source_terms)
      {
	complex coef={c.second.first,c.second.second};
	decltype(auto) p=Q[c.first][isou].getSurelyReadableOn<defaultMemorySpace>();
	
	
	PAR(0,locVol,
	    CAPTURE(coef,
		    TO_WRITE(out),
		    TO_READ(p)),
	    ivol,
	    {
	      spincolor_summ_the_prod_complex(out[ivol],p[ivol],coef);
	    });
      }
    
    MASTER_PRINTF("Source created\n");
  }
  
  void generate_photon_source(LxField<spin1field>& photon_eta)
  {
    auto source_filler=field_rng_stream.getDrawer<spin1field>();
    source_filler.fillField(photon_eta);
    
    PAR(0,locVol,
	CAPTURE(TO_WRITE(photon_eta)),
	ivol,
	{
	  for(int mu=0;mu<NDIM;mu++)
	    z4Transform(photon_eta[ivol][mu]);
	});
  }
  
  //multiply for the Dirac operator
  void mult_by_Dop(LxField<spincolor>& out,
		   const LxField<spincolor>& in,
		   const double kappa,
		   const double mass,
		   const int r,
		   const double& charge,
		   const Momentum& theta)
  {
    LxField<spincolor> tmp("tmp",WITH_HALO);
    
    //rotate the source index - the dirac operator rotate opposite of the mass term
    if(twisted_run>0)
      PAR(0,locVol,
	  CAPTURE(r,
		  TO_WRITE(tmp),
		  TO_READ(in)),
	  ivol,
	  {
	    safe_dirac_prod_spincolor(tmp[ivol],(tau3[r]==+1)?Pminus:Pplus,in[ivol]);
	  });
    else tmp=in;
    
    LxField<quad_su3>& conf=
      *get_updated_conf(charge,theta,*glb_conf);
    
    if(clover_run)
      apply_tmclovQ(out,conf,kappa,*Cl,mass,tmp);
    else
      apply_tmQ(out,conf,kappa,mass,tmp);
    
    PAR(0,locVol,
	CAPTURE(TO_WRITE(out)),
	ivol,
	{
	  safe_dirac_prod_spincolor(out[ivol],base_gamma[5],out[ivol]);
	});
    
    //rotate the sink index
    if(twisted_run>0)
      PAR(0,locVol,
	  CAPTURE(r,
		  TO_WRITE(out)),
	  ivol,
	  {
	    safe_dirac_prod_spincolor(out[ivol],(tau3[r]==+1)?Pminus:Pplus,out[ivol]);
	  });
  }
  
  /// Choose a single position
  void select_position(LxField<spincolor>& out,
		       const LxField<spincolor>& ori,
		       const int& pos)
  {
    const Coords g=
      glbCoordOfGlblx(pos);
    
    PAR(0,locVol,
	CAPTURE(g,
		TO_WRITE(out),
		TO_READ(ori)),
	ivol,
	{
	  bool mask=true;
	  for(int mu=0;mu<NDIM;mu++)
	    mask&=rel_coord_of_loclx(ivol,mu)==g[mu];
	  
	  spincolor_prod_double(out[ivol],ori[ivol],mask);
	});
  }
  
  /// Choose a single spin index s
  void select_spin(LxField<spincolor>& out,
		   const LxField<spincolor>& ori,
		   const int& s)
  {
    PAR(0,locVol,
	CAPTURE(s,
		TO_WRITE(out),
		TO_READ(ori)),
	ivol,
	{
	  for(int id=0;id<NDIRAC;id++)
	    color_prod_double(out[ivol][id],ori[ivol][id],id==s);
	});
  }
  
  /// Choose a single color index c
  void select_color(LxField<spincolor>& out,
		    const LxField<spincolor>& ori,
		    const int& c)
  {
    PAR(0,locVol,
	CAPTURE(c,
		TO_WRITE(out),
		TO_READ(ori)),
	ivol,
	{
	  for(int id=0;id<NDIRAC;id++)
	    for(int ic=0;ic<NCOL;ic++)
	      complex_prod_double(out[ivol][id][ic],ori[ivol][id][ic],ic==c);
	});
  }
  
  enum class BwFw {BW,FW};
  
  CUDA_HOST_AND_DEVICE
  double HeavyTheta(const int x)
  {
    return ((x>=0)+(x>0))/2.0;
  }
  
  //generate a sequential source
  void generate_source(insertion_t inser,char *ext_field_path,double mass,int r,double charge,double kappa,const Momentum& kappa_asymm,const Momentum& theta,std::vector<source_term_t>& source_terms,int isou,int t)
  {
    source_time-=take_time();
    
    int rel_t=t;
    if(rel_t!=-1) rel_t=(t+source_coord[0])%glbSize[0];
    
    LxField<quad_su3> *conf;
    if(not is_smearing_ins(inser)) conf=get_updated_conf(charge,theta,*glb_conf);
    else
      {
	LxField<quad_su3> *ext_conf;
	if(ape_smeared_conf) ext_conf=ape_smeared_conf;
	else                 ext_conf=glb_conf;
	conf=get_updated_conf(0.0,theta,*ext_conf);
      }
    
    LxField<spincolor> ori("ori",WITH_HALO);
    build_source(ori,source_terms,isou);
    
    LxField<spin1field> *ext_field=nullptr;
    if(inser==EXT_FIELD)
      {
	const std::string path=combine("%s/%s",outfolder,ext_field_path);
	ext_field=new LxField<spin1field>("ext_field",WITH_HALO);
	ReadWriteRealVector<spin1field,defaultSpaceTimeLayout,defaultMemorySpace> r(*ext_field,ext_field_path);
	r.read();
      }
    
    /// Function to insert the virtual photon emission projection
    auto vphotonInsertCurr=
      [mass,theta](const BwFw bwFw,const int nu)
    {
      gauge_info insPhoton;
      insPhoton.which_gauge=photon.which_gauge;
      insPhoton.bc[0]=0;
      for(int mu=1;mu<NDIM;mu++)
	insPhoton.bc[mu]=theta[mu];
      insPhoton.c1=WILSON_C1;
      insPhoton.zms=photon.zms;
      
      const double Eg=gluon_energy(insPhoton,mass,0);
      
      return
	[bwFw,nu,Eg,theta] CUDA_HOST_AND_DEVICE(complex ph,
						const int ivol,
						const int mu,
						const double fwbw_phase)
      {
	const double a=-3*0.5*fwbw_phase*M_PI*theta[mu]/glbSize[mu];
	
	complex_iexp(ph,a);
	complex_prodassign_idouble(ph,-1.0);
	
	if(mu==nu)
	  {
	    const int TH=glbSize[0]/2;
	    const int t=(glbSize[0]+glbCoordOfLoclx[ivol][0]-source_coord[0])%glbSize[0];
	    
	    const double f1= //eq.3.8 of reph.pdf
	      (bwFw==BwFw::BW)?
	      exp(-(TH+t)*Eg):
	      exp(-(TH-t)*Eg);
	    
	    const double f2=
	      (bwFw==BwFw::BW)?
	      exp((TH-t)*Eg):
	      exp(-(3*TH-t)*Eg);
	    
	    const double h1=HeavyTheta(TH-t);
	    const double h2=HeavyTheta(t-TH);
	    
	    complex_prodassign_double(ph,h1*f1+h2*f2);
	  }
	else
	  complex_put_to_zero(ph);
      };
    };
    
    MASTER_PRINTF("Inserting r: %d\n",r);
    LxField<spincolor>& loop_source=
      *nissa::loop_source;
    
    switch(inser)
      {
      case PROP:
      case SCALAR:prop_multiply_with_gamma(loop_source,0,ori,rel_t);break;
      case PSEUDO:prop_multiply_with_gamma(loop_source,5,ori,rel_t);break;
      case GAMMA:prop_multiply_with_gamma(loop_source,r,ori,rel_t);break;
      case COLOR:prop_multiply_with_color_delta(loop_source,r,ori,rel_t);break;
      case PHOTON:insert_external_source(loop_source,*conf,*photon_field,ori,rel_t,r,allDirs,loc_hadr_curr);break;
      case PHOTON0:insert_external_source(loop_source,*conf,*photon_field,ori,rel_t,r,onlyDir[0],loc_hadr_curr);break;
      case PHOTON1:insert_external_source(loop_source,*conf,*photon_field,ori,rel_t,r,onlyDir[1],loc_hadr_curr);break;
      case PHOTON2:insert_external_source(loop_source,*conf,*photon_field,ori,rel_t,r,onlyDir[2],loc_hadr_curr);break;
      case PHOTON3:insert_external_source(loop_source,*conf,*photon_field,ori,rel_t,r,onlyDir[3],loc_hadr_curr);break;
      case VBHOTON0:insert_external_source(loop_source,*conf,vphotonInsertCurr(BwFw::BW,0),ori,rel_t,r,loc_hadr_curr);break;
      case VBHOTON1:insert_external_source(loop_source,*conf,vphotonInsertCurr(BwFw::BW,1),ori,rel_t,r,loc_hadr_curr);break;
      case VBHOTON2:insert_external_source(loop_source,*conf,vphotonInsertCurr(BwFw::BW,2),ori,rel_t,r,loc_hadr_curr);break;
      case VBHOTON3:insert_external_source(loop_source,*conf,vphotonInsertCurr(BwFw::BW,3),ori,rel_t,r,loc_hadr_curr);break;
      case VPHOTON0:insert_external_source(loop_source,*conf,vphotonInsertCurr(BwFw::FW,0),ori,rel_t,r,loc_hadr_curr);break;
      case VPHOTON1:insert_external_source(loop_source,*conf,vphotonInsertCurr(BwFw::FW,1),ori,rel_t,r,loc_hadr_curr);break;
      case VPHOTON2:insert_external_source(loop_source,*conf,vphotonInsertCurr(BwFw::FW,2),ori,rel_t,r,loc_hadr_curr);break;
      case VPHOTON3:insert_external_source(loop_source,*conf,vphotonInsertCurr(BwFw::FW,3),ori,rel_t,r,loc_hadr_curr);break;
      case PHOTON_PHI:insert_external_source(loop_source,*conf,*photon_phi,ori,rel_t,r,allDirs,loc_hadr_curr);break;
      case PHOTON_ETA:insert_external_source(loop_source,*conf,*photon_eta,ori,rel_t,r,allDirs,loc_hadr_curr);break;
      case TADPOLE:insert_tadpole(loop_source,*conf,ori,rel_t,r);break;
      case CVEC:insert_conserved_current(loop_source,*conf,ori,rel_t,r,allDirs);break;
      case CVEC0:insert_conserved_current(loop_source,*conf,ori,rel_t,r,onlyDir[0]);break;
      case CVEC1:insert_conserved_current(loop_source,*conf,ori,rel_t,r,onlyDir[1]);break;
      case CVEC2:insert_conserved_current(loop_source,*conf,ori,rel_t,r,onlyDir[2]);break;
      case CVEC3:insert_conserved_current(loop_source,*conf,ori,rel_t,r,onlyDir[3]);break;
      case EXT_FIELD:insert_external_source(loop_source,*conf,*ext_field,ori,rel_t,r,allDirs,loc_hadr_curr);break;
      case SMEARING:smear_prop(loop_source,*conf,ori,rel_t,kappa,r);break;
      case ANYSM:smear_prop(loop_source,*conf,ori,rel_t,kappa_asymm,r);break;
      case WFLOW:flow_prop(loop_source,*conf,ori,rel_t,kappa,r);break;
      case BACK_WFLOW:back_flow_prop(loop_source,*conf,ori,rel_t,kappa,r);break;
      case PHASING:phase_prop(loop_source,ori,rel_t,theta);break;
      case DIROP:mult_by_Dop(loop_source,ori,kappa,mass,r,charge,theta);break;
      case DEL_POS:select_position(loop_source,ori,r);break;
      case DEL_SPIN:select_spin(loop_source,ori,r);break;
      case DEL_COL:select_color(loop_source,ori,r);break;
      case LEP_LOOP:insert_external_source(loop_source,*conf,get_lepton_loop(mass,theta),ori,rel_t,r,allDirs,loc_hadr_curr);break;
      }
    
    if(ext_field)
      delete ext_field;
    
    source_time+=take_time();
    nsource_tot++;
  }
  
  //generate all the quark propagators
  void generate_quark_propagator(std::string& name,
				 qprop_t& q,
				 int ihit)
  {
    q.alloc_storage();
    
    //get ori_source norm2
    const std::string& first_source=q.source_terms.front().first;
    const double ori_source_norm2=q.ori_source_norm2=Q[first_source].ori_source_norm2;
    for(auto& n : q.source_terms)
      {
   	double this_source_norm2=Q[n.first].ori_source_norm2;
	if(ori_source_norm2!=this_source_norm2)
	  CRASH("first source %s has different norm2 %lg than %s, %lg",first_source.c_str(),ori_source_norm2,n.first.c_str(),this_source_norm2);
      }
    
    //write info on mass and r
    if(twisted_run) MASTER_PRINTF(" mass=%lg, r=%d, theta={%lg,%lg,%lg}\n",q.mass,q.r,q.theta[1],q.theta[2],q.theta[3]);
    else            MASTER_PRINTF(" kappa=%lg, theta={%lg,%lg,%lg}\n",q.kappa,q.theta[1],q.theta[2],q.theta[3]);
    
    //compute the inverse clover term, if needed
    if(clover_run and q.insertion==PROP)
      {
	static double m=0,k=0;
	if(m!=q.mass or k!=q.kappa)
	  {
	    m=q.mass;
	    k=q.kappa;
	    const double init_time=take_time();
	    MASTER_PRINTF("Inverting clover\n");
	    invert_twisted_clover_term(*invCl,q.mass,q.kappa,*Cl);
	    MASTER_PRINTF("Clover inverted in %lg s\n",take_time()-init_time);
	  }
      }
    
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
    MASTER_PRINTF("Generating propagator %s inserting %s on source %s\n",name.c_str(),ins_name[insertion],source_descr.c_str());
    for(int id_so=0;id_so<nso_spi;id_so++)
      for(int ic_so=0;ic_so<nso_col;ic_so++)
	{
	  int isou=so_sp_col_ind(id_so,ic_so);
	  generate_source(insertion,q.ext_field_path,q.mass,q.r,q.charge,q.kappa,q.kappa_asymm,q.theta,q.source_terms,isou,q.tins);
	  q[isou].initOn<defaultMemorySpace>([&name,
					      ihit,
					      id_so,
					      ic_so,
					      &q](LxField<spincolor>& sol)
	  {
	    //combine the filename
	    const std::string path=combine("%s/hit%d_prop%s_idso%d_icso%d",outfolder,ihit,name.c_str(),id_so,ic_so);
	    
	    ReadWriteRealVector<spincolor,defaultSpaceTimeLayout,defaultMemorySpace> rw(sol,path);
	    
	    //if the prop exists read it
	    if(rw.canLoad())
	      {
		MASTER_PRINTF("  loading the solution, dirac index %d, color %d\n",id_so,ic_so);
		START_TIMING(read_prop_time,nread_prop);
		rw.read();
		STOP_TIMING(read_prop_time);
	      }
	    else
	      {
		//otherwise compute it
		if(q.insertion==PROP) get_qprop(sol,*loop_source,q.kappa,q.mass,q.r,q.charge,q.residue,q.theta);
		else    *sol=*loop_source;
		
		//and store if needed
		if(q.store)
		  {
		    START_TIMING(store_prop_time,nstore_prop);
		    rw.write();
		    STOP_TIMING(store_prop_time);
		  }
		MASTER_PRINTF("  finished the calculation of dirac index %d, color %d\n",id_so,ic_so);
	      }
	  });
	}
  }
  
  /////////////////////////////////////////////// photon propagators ///////////////////////////////////////////
  
  //allocate the photon fields
  void allocate_photon_fields()
  {
    photon_eta=new LxField<spin1field>("photon_eta",WITH_HALO);
    photon_field=new LxField<spin1field>("photon_field",WITH_HALO);
    photon_phi=new LxField<spin1field>("photon_phi",WITH_HALO);
  }
  
  //free the photon fields
  void free_photon_fields()
  {
    delete photon_eta;
    delete photon_phi;
    delete photon_field;
  }
  
  //wrapper to generate a stochastic propagator
  void generate_photon_stochastic_propagator(const int& ihit)
  {
    photon_prop_time-=take_time();
    
    generate_photon_source(*photon_eta);
    
    //generate source and stochastich propagator
    MASTER_PRINTF("Generating photon stochastic propagator\n");
    multiply_by_tlSym_gauge_propagator(*photon_phi,*photon_eta,photon);
    multiply_by_sqrt_tlSym_gauge_propagator(*photon_field,*photon_eta,photon);
    
    std::vector<std::pair<std::string,LxField<spin1field>*>> name_field;
    name_field.push_back(std::make_pair("eta",photon_eta));
    name_field.push_back(std::make_pair("phi",photon_phi));
    name_field.push_back(std::make_pair("A",photon_field));
    for(auto nf=name_field.begin();nf!=name_field.end();nf++)
      {
	std::string &name=(*nf).first;
	LxField<spin1field> *ph=(*nf).second;
	
	//combine the filename
	const std::string path=
	  combine("%s/hit%d_field%s",outfolder,ihit,name.c_str());
	
	ReadWriteRealVector<spin1field,defaultSpaceTimeLayout,defaultMemorySpace> rw(*ph,path);
	
	//if asked and the file exists read it
	if(load_photons)
	  {
	    if(rw.canLoad())
	      {
		MASTER_PRINTF("  loading the photon field %s\n",name.c_str());
		START_TIMING(read_prop_time,nread_prop);
		rw.read();
		STOP_TIMING(read_prop_time);
	      }
	    else MASTER_PRINTF("  file %s not available, skipping loading\n",path.c_str());
	  }
	
	//if asked, write it
	if(store_photons)
	  {
	    MASTER_PRINTF("  storing the photon field %s\n",name.c_str());
	    START_TIMING(store_prop_time,nstore_prop);
	    rw.write();
	    STOP_TIMING(store_prop_time);
	  }
      }
    
    //invalidate internal conf
    inner_conf_valid=false;
    
    photon_prop_time+=take_time();
    nphoton_prop_tot++;
  }
  
  //put the phase of the source due to missing e(iky)
  void put_fft_source_phase(LxField<spincolor>& qtilde,
			    const double& fft_sign)
  {
    PAR(0,locVol,
	CAPTURE(fft_sign,
		TO_WRITE(qtilde)),
	imom,
	{
	  //sum_mu -sign*2*pi*p_mu*y_mu/L_mu
	  double arg=0;
	  for(int mu=0;mu<NDIM;mu++)
	    arg+=-fft_sign*2*M_PI*glbCoordOfLoclx[imom][mu]*source_coord[mu]/glbSize[mu];
	  
	  complex f={cos(arg),sin(arg)};
	  spincolor_prodassign_complex(qtilde[imom],f);
	  
	  // spincolor_put_to_zero(qtilde[imom]);
	  // for(int mu=0;mu<4;mu++) qtilde[imom][mu][0][0]=glb_coord_of_loclx[imom][mu];
	});
  }
  
  //initialize the fft filter, once forever
  void init_fft_filter_from_range(std::vector<std::pair<fft_mom_range_t,double>>& fft_mom_range_list)
  {
    MASTER_PRINTF("Initializing fft filter\n");
    
    //file where to store output
    FILE *fout=NULL;
    const char path_list[]="mom_list.txt";
    if(not file_exists("mom_list.txt")) fout=open_file(path_list,"w");
    
    //store the list of filtered
    std::set<int> list_of_filtered;
    
    //scattering list
    all_to_all_scattering_list_t sl;
    for(auto &f : fft_mom_range_list)
      for(int vol=volOfLx(f.first.width),ifilt=0;ifilt<vol;ifilt++)
	{
	  //gets the coordinate in the filtering volume
	  Coords c=coordOfLx(ifilt,f.first.width);
	  coordsSummassign(c,f.first.offs,glbSize);
	  
	  //compute p~4/p~2^2
	  double pt2=0,pt4=0;
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      double pmu=M_PI*(2*c[mu]+(mu==0)*temporal_bc)/glbSize[mu];
	      double ptmu=sin(pmu);
	      pt2+=sqr(ptmu);
	      pt4+=pow(ptmu,4.0);
	    }
	  
	  double p4_fr_p22_max=f.second;
	  if(pt4/sqr(pt2)<p4_fr_p22_max)
	    for(int imir=0;imir<pow(2,NDIM);imir++)
	      {
		//get mirrorized
		Coords cmir;
		for(int mu=0;mu<NDIM;mu++)
		  cmir[mu]=get_mirrorized_site_coord(c[mu]+(mu==0 and get_bit(imir,0) and temporal_bc==ANTIPERIODIC_BC),mu,get_bit(imir,mu));
		
		//check if not already collected
		int iglb=glblxOfCoord(cmir);
		if(list_of_filtered.find(iglb)==list_of_filtered.end())
		  {
		    //print momentum coordinates
		    if(fout)
		      {
			for(int mu=0;mu<NDIM;mu++)
			  if(cmir[mu]<glbSize[mu]/2) master_fprintf(fout,"%d ",cmir[mu]);
			  else                        master_fprintf(fout,"%d ",cmir[mu]-glbSize[mu]);
			master_fprintf(fout,"\n");
		      }
		    
		    //search where data is stored
		    const auto [wrank,iloc]=
		      getLoclxAndRankOfCoords(cmir); //the remapper will leave holes
		    if(rank==wrank) sl.push_back(std::make_pair(iloc,list_of_filtered.size()*nranks));
		    
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
    
    std::set<int> list_of_filtered;;
    
    all_to_all_scattering_list_t sl;
    bool ended=false;
    while(not ended)
      {
	Coords c;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    int cmu;
	    bool emu;
	    
	    std::tie(cmu,emu)=do_on_master(read_int);
	    
	    c[mu]=(cmu+glbSize[mu])%glbSize[mu];
	    
	    ended|=emu;
	  }
	
	if(not ended)
	  {
	    //search where data is stored
	    const auto [wrank,iloc]=
	      getLoclxAndRankOfCoords(c);
	    if(rank==wrank) sl.push_back(std::make_pair(iloc,list_of_filtered.size()*nranks));
	    
	    int iglb=glblxOfCoord(c);
	    list_of_filtered.insert(iglb);
	  }
      }
    
    //close file if opened
    close_file(fin);
    
    fft_filterer.emplace_back(list_of_filtered.size(),sl,fileout_suff);
  }
  
  //perform fft and store the propagators
  void propagators_fft(const int& ihit)
  {
    LxField<spincolor> qtilde("qtilde",WITH_HALO);
    
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
	CRASH("dependencies are broken");
	const std::string tag=fft_prop_list[iprop];
	MASTER_PRINTF("Fourier transforming propagator %s\n",tag.c_str());
	
	//loop on dirac and color source index
	for(int id_so=0;id_so<nso_spi;id_so++)
	  for(int ic_so=0;ic_so<nso_col;ic_so++)
	    {
	      START_TIMING(fft_time,nfft_tot);
	      
	      //perform fft
	      // LxField<spincolor>& q=Q[tag][so_sp_col_ind(id_so,ic_so)];
	      
	      CRASH("reimplement");
	      //fft4d((complex*)qtilde,(complex*)q,sizeof(spincolor)/sizeof(complex),fft_sign,1);
	      put_fft_source_phase(qtilde,fft_sign);
	      
	      //gather - check the rewriting pattern above!
	      for(int i=0;i<nf;i++)
		{
		  MASTER_PRINTF("Filtering %d/%d\n",i,nf);
		  //fft_filterer[i].fft_filter_remap.communicate(qfilt_temp[i],qtilde,sizeof(spincolor));
		  CRASH("#warning reimplement");
		  // NISSA_PARALLEL_LOOP(imom,0,fft_filterer[i].nfft_filtered)
		  //   spincolor_copy(qfilt[i][imom*nso_spi*nso_col+so_sp_col_ind(id_so,ic_so)],qfilt_temp[i][imom]);
		  // NISSA_PARALLEL_LOOP_END;
		  // set_borders_invalid(qfilt[i]);
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
  }
}
