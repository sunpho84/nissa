#pragma once

double calculate_weighted_residue_RL(spincolor *source,spincolor *sol,quad_su3 *conf,double kappa,double m,spincolor *s,spincolor *t,int dinf,int RL)
{
  apply_Q2_RL(s,sol,conf,kappa,m,t,RL);
  
  double loc_weighted_residue=0,tot_weighted_residue;
  double loc_weight=0,tot_weight;
  double *ds=(double*)s,*dsource=(double*)source,*dsol=(double*)sol;
  
  for(int i=0;i<loc_vol*3*4;i++)
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
  
  if(dinf==2) tot_weighted_residue/=loc_vol*tot_weight;
  
  return tot_weighted_residue;
}

int check_cgmms_residue_RL(int *run_flag,double *residue_mass,int nrun,double rr,double *zfs,int st_crit,double st_res,double st_res2,int iter,spincolor **sol,int nmass,double *m,spincolor *source,quad_su3 *conf,double kappa,spincolor *s,spincolor *t,int RL)
{
  const int each=10;

  for(int imass=0;imass<nmass;imass++)
    if(run_flag[imass])
      {
	int fini=0;

	if(st_crit==sc_standard||st_crit==sc_unilevel)
	  {
	    residue_mass[imass]=rr*zfs[imass]*zfs[imass];
	    
	    if(st_crit==sc_standard) fini=(residue_mass[imass]<st_res2||residue_mass[0]<st_res);
	    else fini=residue_mass[imass]<st_res;
	  }
	else if(st_crit==sc_weighted_norm2||st_crit==sc_weighted_norm_inf)
	  if(iter%each==0)
	    {  //locally weighted norm
	      communicate_lx_spincolor_borders(sol[imass]);
	      if(st_crit==sc_weighted_norm2)
		residue_mass[imass]=calculate_weighted_residue_RL(source,sol[imass],conf,kappa,m[imass],s,t,2,RL);
	      else residue_mass[imass]=calculate_weighted_residue_RL(source,sol[imass],conf,kappa,m[imass],s,t,-1,RL);
	
	      if(residue_mass[imass]<st_res) fini=1;
	    }
	
	if(fini)
	  {
	    run_flag[imass]=0;
	    nrun--;
	  }
      } 
  
  if(iter%each==0 && rank==0)
    {    
      printf("cgmms iter %d residues: ",iter);
      for(int imass=0;imass<nmass;imass++) if(run_flag[imass]) printf("%1.4e  ",residue_mass[imass]);
      printf("\n");
    }
  
  return nrun;
}
