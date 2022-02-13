#include <nissa.hpp>

#define EXTERN_CONF
 #include "conf.hpp"

#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

namespace nissa
{
  //init the MPI grid
  void read_init_grid()
  {
    int L,T;
    read_str_int("L",&L);
    read_str_int("T",&T);
    
    init_grid(T,L);
  }
  
  //needed to avoid any check
  bool finish_file_present()
  {
    return not file_exists(combine("%s/finished",outfolder).c_str());
  }
  
  //allocate confs needed by the program
  void allocate_confs()
  {
    if(not conf_allocated)
      {
	master_printf("Allocating confs\n");
	
	glb_conf=nissa_malloc("glb_conf",locVolWithBordAndEdge.nastyConvert(),quad_su3);
	inner_conf=nissa_malloc("inner_conf",locVolWithBordAndEdge.nastyConvert(),quad_su3);
	ape_smeared_conf=nissa_malloc("ape_smeared_conf",locVolWithBordAndEdge.nastyConvert(),quad_su3);
      }
    else
      master_printf("Skipping allocating confs\n");
    
    conf_allocated=true;
  }
  
  //freeing the confs
  void free_confs()
  {
    if(conf_allocated)
      {
	master_printf("Freeing confs\n");
	
	nissa_free(glb_conf);
	nissa_free(inner_conf);
	nissa_free(ape_smeared_conf);
      }
    else
      master_printf("Skipping freeing confs\n");
    
    conf_allocated=false;
  }
  
  using Su3=Tensor<OfComps<ColorRow,ColorCln,ComplId>,double>;
  
  using ComplDouble=Tensor<OfComps<ComplId>,double>;
  
  //The product of two complex number
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void unsafeComplexProd(A&& a,
			 const B& b,
			 const C& c)
  {
    double t=b(Re)*c(Re)-b(Im)*c(Im);
    a(Im)=b(Re)*c(Im)+b(Im)*c(Re);
    a(Re)=t;
  }
  
  //The product of a complex number by the conjugate of the second
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void unsafeComplexConj2Prod(A&& a,
			      const B& b,
			      const C& c)
  {
    //ASM_BOOKMARK_BEGIN("unsafeComplexConj2Prod");
    unsafeComplexProd(a,b,conj(c));
    // a(Re)=+b(Re)*c(Re)+b(Im)*c(Im);
    // a(Im)=-b(Re)*c(Im)+b(Im)*c(Re);
    //ASM_BOOKMARK_END("unsafeComplexConj2Prod");
  }
  
  void u(complex a,
	 const complex b,
	 const complex c)
  {
    ASM_BOOKMARK_BEGIN("unsafe_complex_conj2_prod");
    unsafe_complex_conj2_prod(a,b,c);
    ASM_BOOKMARK_END("unsafe_complex_conj2_prod");
  }
  
  void u(su3 a,
	 const su3 b,
	 const su3 c)
  {
    ASM_BOOKMARK_BEGIN("safe_su3_prod_su3_dag");
    safe_su3_prod_su3_dag(a,b,c);
    ASM_BOOKMARK_END("safe_su3_prod_su3_dag");
  }
  
  void u(Tensor<OfComps<ComplId>,double>& a,
	 const Tensor<OfComps<ComplId>,double>& b,
	 const Tensor<OfComps<ComplId>,double>& c)
  {
    ASM_BOOKMARK_BEGIN("unsafeComplexConj2Prod");
    unsafeComplexConj2Prod(a,b,c);
    ASM_BOOKMARK_END("unsafeComplexConj2Prod");
  }
  
  //Summ to the output the product of two complex number
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void complexSummTheProd(A&& a,
						  const B& b,
						  const C& c)
  {
    const double t=b(Re)*c(Re)-b(Im)*c(Im);
    a(Im)+=b(Re)*c(Im)+b(Im)*c(Re);
    a(Re)+=t;
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void complexSummTheConj2Prod(A&& a,
			       const B& b,
			       const C& c)
  {
    const double t=+b(Re)*c(Re)+b(Im)*c(Im);
    a(Im)+=-b(Re)*c(Im)+b(Im)*c(Re);
    a(Re)+=t;
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void complexCopy(A&& a,
		   const B& b)
  {
    UNROLL_FOR_ALL_COMPONENT_VALUES(ComplId,reim)
      a(reim)=b(reim);
    UNROLL_FOR_END;
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_DEVICE inline void complexSummassign(A&& a,
						 const B& b)
  {
    UNROLL_FOR_ALL_COMPONENT_VALUES(ComplId,reim)
      a(reim)+=b(reim);
    UNROLL_FOR_END;
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void su3Copy(A&& a,
	       const B& b)
  {
    UNROLL_FOR_ALL_COMPONENT_VALUES(ColorRow,ir)
      UNROLL_FOR_ALL_COMPONENT_VALUES(ColorCln,ic)
      complexCopy(a(ir,ic),b(ir,ic));
    UNROLL_FOR_END;
    UNROLL_FOR_END;
  }
  
  //Product of two su3 matrixes
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void unsafeSu3ProdSu3(A&& a,
			const B& b,
			const C& c)
  {
    // unrollFor<NCOL>([&](const ColorRow& colRowOut) INLINE_ATTRIBUTE
    // {
    //   unrollFor<NCOL>([&](const ColorCln& colClnOut) INLINE_ATTRIBUTE
    //   {
    // 	  unsafeComplexProd(a(colRowOut,colClnOut),b(colRowOut,ColorCln(0)),c(ColorRow(0),colClnOut));
    // 	  for(int itemp=1;itemp<NCOL;itemp++)
    // 	    complexSummTheProd(a(colRowOut,colClnOut),b(colRowOut,ColorCln(itemp)),c(ColorRow(itemp),colClnOut));
    //   });
    // });
    UNROLL_FOR_ALL_COMPONENT_VALUES(ColorRow,colRowOut)
      UNROLL_FOR_ALL_COMPONENT_VALUES(ColorCln,colClnOut)
      {
	unsafeComplexProd(a(colRowOut,colClnOut),b(colRowOut,ColorCln(0)),c(ColorRow(0),colClnOut));
	for(int itemp=1;itemp<NCOL;itemp++)
	  complexSummTheProd(a(colRowOut,colClnOut),b(colRowOut,ColorCln(itemp)),c(ColorRow(itemp),colClnOut));
      }
    UNROLL_FOR_END;
    UNROLL_FOR_END;
  }
  
  //Product of two su3 matrixes
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void unsafeSu3ProdSu3Dag(A&& a,
			   const B& b,
			   const C& c)
  {
    unsafeSu3ProdSu3(a,b,dag(c));
    // FOR_ALL_COMPONENT_VALUES(ColorRow,colRowOut)
    //   FOR_ALL_COMPONENT_VALUES(ColorCln,colClnOut)
    // 	{
    // 	  unsafeComplexConj2Prod(a(colRowOut,colClnOut),b(colRowOut,ColorCln(0)),c(colClnOut.transp(),ColorCln(0)));
    // 	  for(int jc=1;jc<NCOL;jc++)
    // 	    complexSummTheConj2Prod(a(colRowOut,colClnOut),b(colRowOut,ColorCln(jc)),c(colClnOut.transp(),ColorCln(jc)));
    // 	}
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void safeSu3ProdSu3Dag(A&& a,
			 const B& b,
			 const C& c)
  {
    Su3 d;
    unsafeSu3ProdSu3Dag(d,b,c);
    su3Copy(a,d);
  }
  
  //return the trace of an su3 matrix
  template <typename A,
	    typename B>
  CUDA_HOST_DEVICE INLINE_FUNCTION
  void su3Trace(A&& tr,
		const B& m)
  {
    complexCopy(tr,m(ColorRow(0),ColorCln(0)));
    for(int ic=1;ic<NCOL;ic++)
      complexSummassign(tr,m(ColorRow(ic),ColorCln(ic)));
  }
  
  //read the conf and setup it
  void setup_conf(quad_su3 *conf,const char *conf_path,int rnd_gauge_transform,int free_theory)
  {
    //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
    if(not free_theory)
      {
	START_TIMING(conf_load_time,nconf_load);
	read_ildg_gauge_conf(conf,conf_path);
	STOP_TIMING(conf_load_time);
	master_printf("plaq: %+16.16g\n",global_plaquette_lx_conf(conf));
      }
    else generate_cold_lx_conf(conf);
    
    //if asked, randomly transform the configurations
    if(rnd_gauge_transform) perform_random_gauge_transform(conf,conf);
    if(Landau_gauge_fix_flag) Landau_or_Coulomb_gauge_fix(conf,&gauge_fixing_pars,conf);
    if(store_conf) write_ildg_gauge_conf(combine("%s/conf",outfolder),conf,64);
    
    //if clover term is included, compute it
    if(clover_run) clover_term(Cl,glb_cSW,conf);
    
    //// DEBUG
    master_printf("DEBUG plaquette: %+16.16lg\n",global_plaquette_lx_conf(conf));
    
    auto newConf=lxField<OfComps<Dir,ColorRow,ColorCln,ComplId>,AllocateBord::YES>();
    
    for(LocLxSite site=0;site<locVol;site++)
      FOR_ALL_DIRS(dir)
	FOR_ALL_ROW_COLORS(rowCol)
          FOR_ALL_CLN_COLORS(clnCol)
	    FOR_REIM_PARTS(reIm)
	      newConf(site,dir,rowCol,clnCol,reIm)=conf[site()][dir()][rowCol()][clnCol()][reIm()];
    
  // using ToBeFiltered = std::tuple<TensorComp<LocLxSiteSignature, ANY, 0>, TensorComp<ColorSignature, ROW, 0>, TensorComp<ColorSignature, CLN, 0>, TensorComp<ComplIdSignature, ANY, 0> >;
  // using Filter = std::tuple<std::tuple<TensorComp<ComplIdSignature, ANY, 0>, TensorComp<LocLxSiteSignature, ANY, 0> >, TensorComp<ColorSignature, ROW, 0>, TensorComp<ColorSignature, CLN, 0> >;
  // using aa=typename TupleFilterAllTypes<ToBeFiltered,Filter>::type;
    
    auto plaquettes=lxField<OfComps<>>();
    for(LocLxSite site=0;site<locVol;site++)
      plaquettes(site)=0.0;
    
    FOR_ALL_DIRS(dir)
      for(Dir otherDir=dir+1;otherDir<NDIM;otherDir++)
	{
	  // auto r=std::decay<decltype(std::get<1>(t.nestedExprs))>::Comps{};
	  // auto shifted=shiftDw(newConf(otherDir),dir).close();
	  auto lowerPart=(newConf(dir)*shiftDw(newConf(otherDir),dir)).close();
	  auto upperPart=(newConf(otherDir)*shiftDw(newConf(dir),otherDir)).close();
	  
	  // const double i=lowerPart(LocLxSite(0),ColorRow(0),ColorCln(0),ComplId(0));
	  // su3 temp;
	  // unsafe_su3_prod_su3(temp,conf[0][dir()],conf[loclxNeighup(LocLxSite(0),dir)()][otherDir()]);
	  // const double j=temp[0][0][0];
	  // master_printf("ori: %lg %lg\n",
	  // 		conf[0][dir()][0][0][0],
	  // 		newConf(LocLxSite(0),dir,ColorRow(0),ColorCln(0),ComplId(0)));
	  // master_printf("shifted: %lg %lg\n",
	  // 		conf[loclxNeighup(LocLxSite(0),dir)()][otherDir()][0][0][0],
	  // 		shifted(LocLxSite(0),ColorRow(0),ColorCln(0),ComplId(0)));
	  
	  // master_printf("%lg %lg\n",i,j);
	  
	  auto c=(lowerPart*dag(upperPart));
	  
	  
	  ASM_BOOKMARK_BEGIN("ciccione");
	  plaquettes=plaquettes+real(trace(c));
	  //printf("%lg\n",glbPlaq);
	  ASM_BOOKMARK_END("ciccione");
	}
    const double glbPlaq=plaquettes.globalReduce()();
    printf("%.16lg\n",glbPlaq/glbVol()/6/3);
    
    Tensor<OfComps<Dir,ColorRow>> rt;
    for (Dir mu = 0; mu < Dir ::sizeAtCompileTimeAssertingNotDynamic(); mu++)
      for(ColorRow cr=0;cr<3;cr++)
	rt(mu,cr)=cr()+3*mu();
    Tensor<OfComps<Dir,ColorCln>> ct;
    FOR_ALL_DIRS(mu)
      for(ColorCln cc=0;cc<3;cc++)
	ct(mu,cc)=cc()+3*mu();
    
    Tensor<OfComps<Dir>> res;
    ASM_BOOKMARK_BEGIN("prod");
    res=prod(ct,rt);
    ASM_BOOKMARK_END("prod");
    master_printf("%lg\n",res(Dir(0)));
    master_printf("%lg\n",res(Dir(1)));
    
    crash("");
    // newConf=dag(oldConf);
    
    // typename decltype(newConf)::Fund a;

    
    // loopOnAllComponents<typename decltype(newConf)::Comps>(newConf.data.indexComputer.dynamicSizes,[](auto...){});
      //decltype(prod(newConf,newField)) b="ciao";


    
    //TransposeTensorComps<decltype(newConf)::Comps> t="ciao";
    
    NISSA_LOC_VOL_LOOP(ivol)
      FOR_ALL_DIRS(mu)
      FOR_ALL_COMPONENT_VALUES(ColorRow,colRow)
      FOR_ALL_COMPONENT_VALUES(ColorCln,colCol)
      FOR_ALL_COMPONENT_VALUES(ComplId,reIm)
      {
	newConf(ivol,mu,colRow,colCol,reIm)=
	  conf[ivol.nastyConvert()][mu.nastyConvert()][colRow.nastyConvert()][colCol.nastyConvert()][reIm.nastyConvert()];
      }
    
    ComplDouble r;
    FOR_ALL_COMPONENT_VALUES(ComplId,reim)
      r(ComplId(reim))=0;
    
    NISSA_LOC_VOL_LOOP(ivol)
      FOR_ALL_DIRS(mu)
      FOR_ALL_DIRS(nu)
      if(mu>nu)
	{
        Su3 p;
	
	unsafeSu3ProdSu3(p,newConf(ivol,mu),newConf(loclxNeighup(ivol,mu),nu));
	safeSu3ProdSu3Dag(p,p,newConf(loclxNeighup(ivol,nu),mu));
	
	ASM_BOOKMARK_BEGIN("safeSu3ProdSu3Dag");
	safeSu3ProdSu3Dag(p,p,newConf(ivol,nu));
	ASM_BOOKMARK_END("safeSu3ProdSu3Dag");
	
	ComplDouble c;
	su3Trace(c,p);
	complexSummassign(r,c);
      }
        //// DEBUG
    master_printf("DEBUG2 plaquette: %+16.16lg\n",r(Re)/18/glbVol());
    
    //if the copied conf exists, ape smear
    if(ape_smeared_conf)
      {
	ape_spatial_smear_conf(ape_smeared_conf,conf,ape_smearing_alpha,ape_smearing_niters);
	master_printf("Smeared plaquette: %+16.16lg\n",global_plaquette_lx_conf(ape_smeared_conf));
      }
    
    //invalidate internal conf
    inner_conf_valid=false;
  }
  
  //take a set of theta, charge and photon field, and update the conf
  quad_su3* get_updated_conf(double charge,const Momentum& theta,quad_su3 *in_conf)
  {
    //check if the inner conf is valid or not
    static quad_su3 *stored_conf=NULL;
    static double stored_charge=0;
    static Momentum stored_theta;
    
    if(not inner_conf_valid) master_printf("Inner conf is invalid (loaded new conf, or new photon generated)\n");
    
    //check ref conf
    if(stored_conf!=in_conf)
      {
	master_printf("Inner conf is invalid (ref conf from %p to %p)\n",stored_conf,in_conf);
	inner_conf_valid=false;
      }
    
    //check charge
    if(charge!=stored_charge)
      {
	master_printf("Inner conf is invalid (charge changed from %lg to %lg)\n",stored_charge,charge);
	inner_conf_valid=false;
      }
    //check theta
    bool same_theta=true;
    FOR_ALL_DIRS(mu)
      same_theta&=(theta(mu)==stored_theta(mu));
    
    if(not same_theta)
      {
	master_printf("Inner conf is invalid (theta changed from {%lg,%lg,%lg,%lg} to {%lg,%lg,%lg,%lg}\n",
		      stored_theta(Dir(0)),stored_theta(xDir),stored_theta(yDir),stored_theta(zDir),theta(Dir(0)),theta(xDir),theta(yDir),theta(zDir));
	inner_conf_valid=false;
      }
    
    if(not inner_conf_valid)
      {
	master_printf("Inner conf not valid: updating it\n");
	
	//copy
	vector_copy(inner_conf,in_conf);
	
	//put momentum
	Momentum old_theta;
	FOR_ALL_DIRS(mu)
	  old_theta(mu)=0.0;
	adapt_theta(inner_conf,old_theta,theta,0,0);
	
	//include the photon field, with correct charge
	if(charge) add_photon_field_to_conf(inner_conf,charge);
      }
    
    //update value and set valid
    stored_conf=in_conf;
    stored_charge=charge;
    FOR_ALL_DIRS(mu)
      stored_theta(mu)=theta(mu);
    inner_conf_valid=true;
    
    return inner_conf;
  }
  
  //check if the time is enough
  int check_remaining_time()
  {
    if(nanalyzed_conf)
      {
	//check remaining time
	double temp_time=take_time()+tot_prog_time;
	double ave_time=temp_time/nanalyzed_conf;
	double left_time=wall_time-temp_time;
	int enough_time=left_time>(ave_time*1.1);
	
	master_printf("\nRemaining time: %lg sec\n",left_time);
	master_printf("Average time per conf: %lg sec, pessimistically: %lg\n",ave_time,ave_time*1.1);
	if(enough_time) master_printf("Time is enough to go on!\n");
	else master_printf("Not enough time, exiting!\n");
	
	return enough_time;
      }
    else return true;
  }
  
  //init a new conf
  void start_new_conf()
  {
    setup_conf(glb_conf,conf_path,rnd_gauge_transform,free_theory);
    
    //reset contractions
    if(mes2pts_contr_size) vector_reset(mes2pts_contr);
    if(handcuffs_contr_size) vector_reset(handcuffs_contr);
    if(bar2pts_contr_size) vector_reset(bar2pts_contr);
    if(bar2pts_alt_contr_size) vector_reset(bar2pts_alt_contr);
    if(nmeslep_corr) vector_reset(meslep_contr);
  }
  
  //skip a number of hits [a,b)
  void skip_nhits(int a,int b)
  {
    for(int ihit=a;ihit<b;ihit++)
      {
	GlbCoords coord;
	generate_random_coord(coord);
	if(need_photon) generate_stochastic_tlSym_gauge_propagator_source(photon_eta);
	generate_original_sources(ihit,true);
      }
  }
  
  //handle to discard the source
  void skip_conf()
  {
    skip_nhits(0,nhits);
  }
  
  //find a new conf
  int read_conf_parameters(int &iconf,bool(*external_condition)())
  {
    //Check if asked to stop or restart
    int asked_stop=file_exists(stop_path);
    verbosity_lv2_master_printf("Asked to stop: %d\n",asked_stop);
    int asked_restart=file_exists("restart");
    verbosity_lv2_master_printf("Asked to restart: %d\n",asked_restart);
    //check if enough time
    int enough_time=check_remaining_time();
    verbosity_lv2_master_printf("Enough time: %d\n",enough_time);
    //check that there are still conf to go
    int still_conf=iconf<ngauge_conf;
    verbosity_lv2_master_printf("Still conf: %d\n",still_conf);
    
    allocate_confs();
    
    int ok_conf=false;
    if(!asked_stop and !asked_restart and enough_time and still_conf)
      do
	{
	  //Gauge path
	  read_str(conf_path,1024);
	  
	  //Out folder
	  read_str(outfolder,1024);
	  
	  //Check if the conf has been finished or is already running
	  master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
	  char run_file[1024];
	  if(snprintf(run_file,1024,"%s/running",outfolder)<0) crash("witing %s",run_file);
	  ok_conf=!(file_exists(run_file)) and external_condition();
	  
	  //if not finished
	  if(ok_conf)
	    {
	      master_printf(" Configuration \"%s\" not yet analyzed, starting\n",conf_path);
	      if(!dir_exists(outfolder))
		{
		  int ris=create_dir(outfolder);
		  if(ris==0) master_printf(" Output path \"%s\" not present, created.\n",outfolder);
		  else
		    {
		      master_printf(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
		      ok_conf=0;
		      skip_conf();
		    }
		}
	      if(ok_conf)
		{
		  //try to lock the running file
		  lock_file.try_lock(run_file);
		  
		  //setup the conf and generate the source
		  start_new_conf();
		  
		  //verify that nobody took the lock
		  if(not lock_file.check_lock())
		    {
		      ok_conf=false;
		      master_printf("Somebody acquired the lock on %s\n",run_file);
		    }
		}
	    }
	  else
	    {
	      //skipping conf
	      master_printf("\"%s\" finished or running, skipping configuration \"%s\"\n",outfolder,conf_path);
	      skip_conf();
	    }
	  iconf++;
	  
	  still_conf=(iconf<ngauge_conf);
	}
      while(!ok_conf and still_conf);
    
    master_printf("\n");
    
    //write if it was asked to stop or restart
    if(asked_stop) master_printf("Asked to stop\n");
    if(asked_restart) master_printf("Asked to restart\n");
    
    //writing that all confs have been measured and write it
    if(!ok_conf and iconf>=ngauge_conf)
      {
	master_printf("Analyzed all confs, exiting\n\n");
	file_touch(stop_path);
      }
    
    return ok_conf;
  }
  
  //mark a conf as finished
  void mark_finished()
  {
    char fin_file[1024];
    if(snprintf(fin_file,1024,"%s/finished",outfolder)<0) crash("writing %s",fin_file);
    file_touch(fin_file);
    nanalyzed_conf++;
  }
  
  inline void print_single_statistic(double frac_time,double tot_time,int niter,const char *tag)
  {if(niter) master_printf(" - %02.2f%% for %d %s (%2.2gs avg)\n",frac_time/tot_time*100,niter,tag,frac_time/niter);}
  
  //print all statisticd
  void print_statistics()
  {
    if(nanalyzed_conf)
      {
	master_printf("\n");
	master_printf("Inverted %d configurations.\n",nanalyzed_conf);
	master_printf("Total time: %g, of which:\n",tot_prog_time);
	print_single_statistic(conf_load_time,tot_prog_time,nconf_load,"loading conf");
	print_single_statistic(smear_oper_time,tot_prog_time,nsmear_oper,"smearing");
	print_single_statistic(lepton_prop_time,tot_prog_time,nlprop,"preparation of lepton propagators");
	print_single_statistic(source_time,tot_prog_time,nsource_tot,"preparation of generalized sources");
	print_single_statistic(inv_time,tot_prog_time,ninv_tot,"calculation of quark propagator");
	if(ninv_tot) master_printf("    of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",
				   cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
	print_single_statistic(store_prop_time,tot_prog_time,nstore_prop,"storing propagators");
	print_single_statistic(read_prop_time,tot_prog_time,nread_prop,"reading propagators");
	print_single_statistic(mes2pts_contr_time,tot_prog_time,nmes2pts_contr_made,"calculation of mesonic 2pts_contractions");
	print_single_statistic(handcuffs_contr_time,tot_prog_time,nhandcuffs_contr_made,"calculation of handcuff 2pts_contractions");
	print_single_statistic(bar2pts_alt_contr_time,tot_prog_time,nbar2pts_alt_contr_made,"calculation of barionic 2pts alt contractions");
	print_single_statistic(bar2pts_contr_time,tot_prog_time,nbar2pts_contr_made,"calculation of barionic 2pts contractions");
	print_single_statistic(meslep_contr_time,tot_prog_time,nmeslep_contr_made,"calculation of hadro-leptonic contractions");
	print_single_statistic(contr_print_time,tot_prog_time,nmeslep_contr_made,"printing contractions");
      	print_single_statistic(fft_time,tot_prog_time,nfft_tot,"Fourier transforming and writing fft-propagators");
      	print_single_statistic(sme_time,tot_prog_time,nsme_tot,"Gaussian smearing");
      	print_single_statistic(flw_time,tot_prog_time,nflw_tot,"Flowing forward");
      	print_single_statistic(bflw_time,tot_prog_time,nbflw_tot,"Flowing backward");
      }
  }
}
