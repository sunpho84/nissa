
//return all the shifts summed together
//placed here because common to 64 and 32 bits template
#if CGM_NARG == 0
void SUMM_SRC_AND_ALL_INV_CGM(BASETYPE* sol,rat_approx_t* appr,int niter_max,double req_res,BASETYPE* source)
#elif CGM_NARG == 1
void SUMM_SRC_AND_ALL_INV_CGM(BASETYPE* sol,AT1 A1,rat_approx_t* appr,int niter_max,double req_res,BASETYPE* source)
#elif CGM_NARG == 2
void SUMM_SRC_AND_ALL_INV_CGM(BASETYPE* sol,AT1 A1,AT2 A2,rat_approx_t* appr,int niter_max,double req_res,BASETYPE* source)
#elif CGM_NARG == 3
void SUMM_SRC_AND_ALL_INV_CGM(BASETYPE* sol,AT1 A1,AT2 A2,AT3 A3,rat_approx_t* appr,int niter_max,double req_res,BASETYPE* source)
#elif CGM_NARG == 4
void SUMM_SRC_AND_ALL_INV_CGM(BASETYPE* sol,AT1 A1,AT2 A2,AT3 A3,AT4 A4,rat_approx_t* appr,int niter_max,double req_res,BASETYPE* source)
#elif CGM_NARG == 5
void SUMM_SRC_AND_ALL_INV_CGM(BASETYPE* sol,AT1 A1,AT2 A2,AT3 A3,AT4 A4,AT5 A5,rat_approx_t* appr,int niter_max,double req_res,BASETYPE* source)
#endif
{
  CRASH("reimplement");
  // const int nterms=appr->degree();
  
  // //allocate temporary single solutions
  // BASETYPE **temp=nissa_malloc("temp",nterms,BASETYPE*);
  // for(int iterm=0;iterm<nterms;iterm++)
  //   temp[iterm]=nissa_malloc(combine("temp%d",iterm).c_str(),BULK_VOL+BORD_VOL,BASETYPE);
  
  // //call multi-shift solver
  // CGM_INVERT_RUN_HM_UP_TO_COMM_PREC(temp,CGM_ADDITIONAL_PARAMETERS_CALL appr->poles.data(),appr->degree(),niter_max,req_res,source);

  // double *weights=nissa_malloc("weights",nterms,double);
  // for(int iterm=0;iterm<nterms;iterm++)
  //   weights[iterm]=appr->weights[iterm];
  // const double cons=appr->cons;
  
  // //summ all the shifts
  // NISSA_PARALLEL_LOOP(i,0,BULK_VOL*NDOUBLES_PER_SITE)
  //   {
  //     ((double*)sol)[i]=cons*((double*)source)[i];
  //     for(int iterm=0;iterm<nterms;iterm++)
  // 	((double*)sol)[i]+=weights[iterm]*((double*)(temp[iterm]))[i];
  //   }
  // NISSA_PARALLEL_LOOP_END;
  
  // set_borders_invalid(sol);
  
  // //free temp vectors
  // for(int iterm=0;iterm<appr->degree();iterm++)
  //   nissa_free(temp[iterm]);
  
  // nissa_free(temp);
  // nissa_free(weights);
}
