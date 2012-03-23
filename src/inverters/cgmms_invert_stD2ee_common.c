#pragma once

double calculate_weighted_residue_stD2ee(color *source,quad_su3 **conf,double m2,color *s,color *t,int dinf,color *sol)
{
  apply_stD2ee(s,conf,t,sqrt(m2),sol);
  
  double loc_weighted_residue=0,tot_weighted_residue;
  double loc_weight=0,tot_weight;
  double *ds=(double*)s,*dsource=(double*)source,*dsol=(double*)sol;
  
  for(int i=0;i<loc_volh*3;i++)
    {
      double nr=(*ds)-(*dsource);
      double ni=(*(ds+1))-(*(dsource+1));
      
      double weight=1/((*dsol)*(*dsol)+(*dsol+1)*(*dsol+1));
      double contrib=(nr*nr+ni*ni)*weight;
      
      if(dinf==2)
        {
          loc_weighted_residue+=contrib;
          loc_weight+=weight;
        }
      else if(contrib>loc_weighted_residue) loc_weighted_residue=contrib;
        
      ds+=2;dsource+=2;dsol+=2;
    }
  
    if(dinf==2)
      {
        MPI_Allreduce(&loc_weighted_residue,&tot_weighted_residue,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&loc_weight,&tot_weight,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
    else MPI_Allreduce(&loc_weighted_residue,&tot_weighted_residue,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  
  if(dinf==2) tot_weighted_residue/=loc_vol*tot_weight/2;
  
  return tot_weighted_residue;
}

int check_cgmm2s_residue_stD2ee(int *run_flag,double *residue_mass,int nrun,double rr,double *zfs,int st_crit,double st_res,double st_res2,int iter,int nmass,double *m2,color *source,quad_su3 **conf,color *s,color *t,double source_norm,color **sol)
{
  const int each=10;

  for(int imass=0;imass<nmass;imass++)
    if(run_flag[imass])
      {
	int fini=0;

	if(st_crit==sc_standard||st_crit==sc_unilevel)
	  {
	    residue_mass[imass]=rr*zfs[imass]*zfs[imass]/source_norm;
	    
	    if(st_crit==sc_standard) fini=(residue_mass[imass]<st_res2||residue_mass[0]<st_res);
	    else fini=residue_mass[imass]<st_res;
	  }
	else if(st_crit==sc_weighted_norm2||st_crit==sc_weighted_norm_inf)
	  if(iter%each==0)
	    {  //locally weighted norm
	      if(st_crit==sc_weighted_norm2)
		residue_mass[imass]=calculate_weighted_residue_stD2ee(source,conf,sqrt(m2[imass]),s,t,2,sol[imass]);
	      else residue_mass[imass]=calculate_weighted_residue_stD2ee(source,conf,sqrt(m2[imass]),s,t,-1,sol[imass]);
	    }
	
	if(fini)
	  {
	    run_flag[imass]=0;
	    nrun--;
	  }
      } 
  
  if(iter%each==0)
    {    
      master_printf(" cgmms iter %d rel. residues: ",iter);
      for(int imass=0;imass<nmass;imass++)
	if(run_flag[imass])
	  master_printf("%1.4e  ",residue_mass[imass]);
      master_printf("\n");
    }
  
  return nrun;
}

//return all the masses summed together
void summ_src_and_all_inv_stD2ee_cgmm2s(color *sol,quad_su3 **conf,rat_approx *appr,int niter,double st_res,double st_minres,int st_crit,color *source)
{
  //allocate temporary single solutions
  color *temp[appr->nterms];
  for(int iterm=0;iterm<appr->nterms;iterm++)
    temp[iterm]=nissa_malloc("temp",loc_volh+loc_bordh,color);
  
  //call multi-mass solver
  inv_stD2ee_cgmm2s(temp,conf,appr->poles,appr->nterms,niter,st_res,st_minres,st_crit,source);
  
  //summ all the masses
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	{
	  sol[ivol][ic][ri]=appr->cons*source[ivol][ic][ri];
	  for(int iterm=0;iterm<appr->nterms;iterm++)
	    sol[ivol][ic][ri]+=appr->weights[iterm]*temp[iterm][ivol][ic][ri];
	}

  set_borders_invalid(sol);
  
  //free temp vectors
  for(int iterm=0;iterm<appr->nterms;iterm++)
    nissa_free(temp[iterm]);
}
