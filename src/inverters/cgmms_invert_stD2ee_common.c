#pragma once

double calculate_weighted_residue_stD2ee(color *source,color *sol,quad_su3 **conf,double m2,color *s,color *t,int dinf)
{
  apply_stD2ee(s,conf,t,sqrt(m2),sol);
  
  double loc_weighted_residue=0,tot_weighted_residue;
  double loc_weight=0,tot_weight;
  double *ds=(double*)s,*dsource=(double*)source,*dsol=(double*)sol;
  
  for(int i=0;i<loc_vol*3/2;i++)
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

int check_cgmm2s_residue_stD2ee(int *run_flag,double *residue_mass,int nrun,double rr,double *zfs,int st_crit,double st_res,double st_res2,int iter,color **sol,int nmass,double *m2,color *source,quad_su3 **conf,color *s,color *t,double source_norm)
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
	      communicate_ev_color_borders(sol[imass]);
	      if(st_crit==sc_weighted_norm2)
		residue_mass[imass]=calculate_weighted_residue_stD2ee(source,sol[imass],conf,sqrt(m2[imass]),s,t,2);
	      else residue_mass[imass]=calculate_weighted_residue_stD2ee(source,sol[imass],conf,sqrt(m2[imass]),s,t,-1);
	    }
	
	if(fini)
	  {
	    run_flag[imass]=0;
	    nrun--;
	  }
      } 
  
  if(iter%each==0)
    {    
      master_printf("cgmms iter %d rel. residues: ",iter);
      for(int imass=0;imass<nmass;imass++)
	if(run_flag[imass])
	  master_printf("%1.4e  ",residue_mass[imass]);
      master_printf("\n");
    }
  
  return nrun;
}

//return all the masses summed together
void summ_src_and_all_inv_stD2ee_cgmm2s(color *sol,color *source,quad_su3 **conf,double con,double *m2,double *coef,int nmass,int niter,double st_res,double st_minres,int st_crit)
{
  //allocate temporary single solutions
  color *temp[nmass];
  temp[0]=nissa_malloc("temp",(loc_vol+loc_bord)/2*nmass,color);
  for(int imass=1;imass<nmass;imass++)
    temp[imass]=temp[imass-1]+(loc_vol+loc_bord)/2;
  
  //coll multi-mass solver
  inv_stD2ee_cgmm2s(temp,source,conf,m2,nmass,niter,st_res,st_minres,st_crit);
  
  //summ all the masses
  for(int ivol=0;ivol<loc_vol/2;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	{
	  sol[ivol][ic][ri]=con*source[ivol][ic][ri];
	  for(int imass=0;imass<nmass;imass++)
	    sol[ivol][ic][ri]+=coef[imass]*temp[imass][ivol][ic][ri];
	}
  
  nissa_free(temp[0]);
}
