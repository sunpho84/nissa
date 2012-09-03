// Template to invert using c.g.
// macro to be defined:
//  -apply_operator
//  -cg_invert
//  -cg_parameters_proto
//  -cg_inner_solver
//  -size_of_bulk, size_of_bord
//  -basetype
//  -ndoubles_per_site

void cg_invert(basetype *sol,basetype *guess,cg_additional_parameters_proto,int niter,int rniter,double residue,basetype *source)
{
  int riter=0;
  basetype *s=nissa_malloc("s",loc_vol,basetype);
  basetype *p=nissa_malloc("p",loc_vol+bord_vol,basetype);
  basetype *r=nissa_malloc("r",loc_vol,basetype);

  //macro to be defined externally, allocating all the required additional vectors
  cg_additional_vectors_allocation();
  
  if(guess==NULL) vector_reset(sol);
  else vector_copy(sol,guess);
  
  //external loop, used if the internal exceed the maximal number of iterations
  double source_norm,lambda;
  do
    {
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      apply_operator(s,cg_operator_parameters,sol);
      
      double_vector_subt_double_vector_prod_double((double*)r,(double*)source,(double*)s,1,bulk_vol*ndoubles_per_site);
      double_vector_copy((double*)p,(double*)r,bulk_vol*ndoubles_per_site);
      source_norm=double_vector_glb_scalar_prod((double*)source,(double*)source,bulk_vol*ndoubles_per_site);
      double delta=double_vector_glb_scalar_prod((double*)r,(double*)r,bulk_vol*ndoubles_per_site);
      
      if(riter==0) verbosity_lv2_master_printf("Source norm: %lg\n",source_norm);
      if(source_norm==0 || isnan(source_norm)) crash("invalid norm: %lg",source_norm);
      
      verbosity_lv2_master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
      
      //main loop
      int iter=0;
      do
	{	  
	  //(r_k,r_k)/(p_k*DD*p_k)
	  apply_operator(s,cg_operator_parameters,p);
	  
	  double alpha=double_vector_glb_scalar_prod((double*)s,(double*)p,bulk_vol*ndoubles_per_site);
	  double omega=delta/alpha;
	  
	  //sol_(k+1)=x_k+omega*p_k
	  double_vector_summ_double_vector_prod_double((double*)sol,(double*)sol,(double*)p,omega,bulk_vol*ndoubles_per_site);
	  //r_(k+1)=x_k-omega*p_k
	  double_vector_summ_double_vector_prod_double((double*)r,(double*)r,(double*)s,-omega,bulk_vol*ndoubles_per_site);
	  //(r_(k+1),r_(k+1))
	  lambda=double_vector_glb_scalar_prod((double*)r,(double*)r,bulk_vol*ndoubles_per_site);

	  //(r_(k+1),r_(k+1))/(r_k,r_k)
	  double gammag=lambda/delta;
	  delta=lambda;
	  
	  //p_(k+1)=r_(k+1)+gammag*p_k
	  double_vector_summ_double_vector_prod_double((double*)p,(double*)r,(double*)p,gammag,bulk_vol*ndoubles_per_site);
	  
	  iter++;

	  if(iter%10==0) verbosity_lv2_master_printf("iter %d relative residue: %lg\n",iter,lambda/source_norm);
	}
      while(lambda>(residue*source_norm) && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      apply_operator(s,cg_operator_parameters,sol);
      double_vector_subt_double_vector_prod_double((double*)r,(double*)source,(double*)s,1,bulk_vol*ndoubles_per_site);
      lambda=double_vector_glb_scalar_prod((double*)r,(double*)r,bulk_vol*ndoubles_per_site);
      
      verbosity_lv1_master_printf("\nfinal relative residue (after %d iters): %lg where %lg was required\n",iter,lambda/source_norm,residue);
      
      riter++;
    }
  while(lambda>(residue*source_norm) && riter<rniter);
  
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  
  //macro to be defined externally
  cg_additional_vectors_free();
}

#undef basetype
#undef ndoubles_per_site
#undef size_of_bulk
#undef size_of_bord

#undef apply_operator
#undef cg_operator_parameters
#undef cg_invert
#undef cg_additional_parameters_proto
#undef cg_additional_vectors_free
#undef cg_additional_vectors_allocation

