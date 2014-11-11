  THREADABLE_FUNCTION_5ARG(APPLY_STD2EE_M2_BGQ, BI_32_64_COLOR*,bi_out, BI_32_64_OCT_SU3**,bi_conf, BI_32_64_COLOR*,bi_temp, PREC_TYPE,mass2, BI_32_64_COLOR*,bi_in)
  {
    GET_THREAD_ID();
    
#ifdef BENCH
    if(IS_MASTER_THREAD) bgq_stdD_app_time-=take_time();
#endif
    
    //----------------------looping on E--------------------
    const int OE=0;

    //compute on the surface and start communications
    NAME3(apply,PREC,staggered_hopping_matrix_oe_or_eo_bgq_nocomm)(bi_conf,0,vsurf_volh,bi_in,OE);
    NAME3(start,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)();
    
    //compute on the bulk and finish communications
    NAME3(apply,PREC,staggered_hopping_matrix_oe_or_eo_bgq_nocomm)(bi_conf,vsurf_volh,loc_volh/2,bi_in,OE);
    NAME3(finish,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)(OE);
    
    //put the eight pieces together
    NAME3(hopping_matrix_eo_or_eo_expand_to,PREC,staggered_D)(bi_out);
    
    //----------------------looping on O--------------------
    const int EO=1;  
    
    //compute on the surface and start communications
    NAME3(apply,PREC,staggered_hopping_matrix_oe_or_eo_bgq_nocomm)(bi_conf,0,vsurf_volh,bi_out,EO);
    NAME3(start,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)();
    
    //compute on the bulk and finish communications
    NAME3(apply,PREC,staggered_hopping_matrix_oe_or_eo_bgq_nocomm)(bi_conf,vsurf_volh,loc_volh/2,bi_out,EO);
    NAME3(finish,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)(EO);
    
    //put the eight pieces subtracting them from diag (in fact one of the two D is daggered)
    if(mass2!=0) NAME3(hopping_matrix_eo_or_eo_expand_to,PREC,staggered_D_subtract_from_mass2_times_in)
		   (bi_out,mass2,bi_in);
    else         NAME3(hopping_matrix_eo_or_eo_expand_to,PREC,staggered_D)(bi_out,-1);

#ifdef BENCH
    if(IS_MASTER_THREAD)
      {
	bgq_stdD_app_time+=take_time();
	bgq_stdD_napp++;
      }
#endif
  }
  THREADABLE_FUNCTION_END

#undef PREC
#undef PREC_TYPE
#undef BI_32_64_OCT_SU3
#undef BI_32_64_COLOR
#undef APPLY_STD2EE_M2_BGQ
