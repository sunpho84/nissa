  THREADABLE_FUNCTION_4ARG(APPLY_STD2EE_M2_BGQ, VIR_32_64_COLOR*,vir_out, VIR_32_64_OCT_SU3**,vir_conf, PREC_TYPE,mass2, VIR_32_64_COLOR*,vir_in)
  {
    GET_THREAD_ID();
    
    START_TIMING(bgq_stdD_app_time,nbgq_stdD_app);
    
    //----------------------looping on E--------------------
    const int OE=0;
    
    //compute on the surface and start communications
    NAME3(apply,PREC,staggered_hopping_matrix_oe_or_eo_bgq_nocomm)(vir_conf,0,vsurf_volh,vir_in,OE);
    NAME3(start,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)();
    
    //compute on the bulk and finish communications
    NAME3(apply,PREC,staggered_hopping_matrix_oe_or_eo_bgq_nocomm)(vir_conf,vsurf_volh,loc_volh/2,vir_in,OE);
    NAME3(finish,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)(OE);
    
    //put the eight pieces together
    NAME3(hopping_matrix_oe_or_eo_expand_to,PREC,staggered_D_bgq)(vir_out);
    
    //----------------------looping on O--------------------
    const int EO=1;
    
    //compute on the surface and start communications
    NAME3(apply,PREC,staggered_hopping_matrix_oe_or_eo_bgq_nocomm)(vir_conf,0,vsurf_volh,vir_out,EO);
    NAME3(start,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)();
    
    //compute on the bulk and finish communications
    NAME3(apply,PREC,staggered_hopping_matrix_oe_or_eo_bgq_nocomm)(vir_conf,vsurf_volh,loc_volh/2,vir_out,EO);
    NAME3(finish,PREC,staggered_hopping_matrix_oe_or_eo_bgq_communications)(EO);
    
    //put the eight pieces subtracting them from diag (in fact one of the two D is daggered)
    if(mass2!=0) NAME3(hopping_matrix_oe_or_eo_expand_to,PREC,staggered_D_subtract_from_mass2_times_in_bgq)
		   (vir_out,mass2,vir_in);
    else         NAME3(hopping_matrix_oe_or_eo_expand_to,PREC,staggered_D_bgq)(vir_out,-1);
    
    STOP_TIMING(bgq_stdD_app_time);
  }
  THREADABLE_FUNCTION_END

#undef PREC
#undef PREC_TYPE
#undef VIR_32_64_OCT_SU3
#undef VIR_32_64_COLOR
#undef APPLY_STD2EE_M2_BGQ
