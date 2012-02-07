#pragma once

//Refers to the doc: "doc/eo_inverter.lyx" for explenations

//invert Koo defined in equation (7)
void inv_tmDkern_eoprec_square_eos(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 **conf,double kappa,double mu,int nitermax,double residue)
{
  int niter=nitermax;
  int riter=0;
  int rniter=5;
  spincolor *p=nissa_malloc("p",loc_volh+loc_bordh,spincolor);
  spincolor *r=nissa_malloc("r",loc_volh,spincolor);
  spincolor *s=nissa_malloc("s",loc_volh,spincolor);
  spincolor *temp1=nissa_malloc("temp1",loc_volh+loc_bordh,spincolor);
  spincolor *temp2=nissa_malloc("temp2",loc_volh+loc_bordh,spincolor);

  ///////////////// prepare the internal source /////////////////
  
  if(guess==NULL) memset(sol,0,sizeof(spincolor)*(loc_volh+loc_bordh));
  else
    {
      communicate_od_spincolor_borders(guess);
      memcpy(sol,guess,sizeof(spincolor)*(loc_volh+loc_bordh));
    }

  //external loop, used if the internal exceed the maximal number of iterations
  double lambda; //(r_(k+1),r_(k+1))
  double source_norm;
  do
    {
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      double delta;
      {
	tmDkern_eoprec_square_eos(s,temp1,temp2,sol,conf,kappa,mu);
	
        double loc_delta=0,loc_source_norm=0;
	for(int X=0;X<loc_volh;X++)
	  for(int id=0;id<4;id++)
	    for(int ic=0;ic<3;ic++)
	      for(int ri=0;ri<2;ri++)
		{
		  double c1=source[X][id][ic][ri]-s[X][id][ic][ri];
		  p[X][id][ic][ri]=r[X][id][ic][ri]=c1;
		  if(riter==0) loc_source_norm+=source[X][id][ic][ri]*source[X][id][ic][ri];
		  loc_delta+=c1*c1;
		}
        MPI_Allreduce(&loc_delta,&delta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if(riter==0)
	  {
	    MPI_Allreduce(&loc_source_norm,&source_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	    master_printf("iter 0 relative residue: %lg\n",delta/source_norm);
	  }
      }
      
      //main loop
      int iter=0;
      do
	{
	  double omega; //(r_k,r_k)/(p_k*DD*p_k)
	  double alpha;
	  
	  tmDkern_eoprec_square_eos(s,temp1,temp2,p,conf,kappa,mu);
	  
	  double loc_alpha=0; //real part of the scalar product
	  for(int X=0;X<loc_volh;X++)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		for(int ri=0;ri<2;ri++)
		  loc_alpha+=s[X][id][ic][ri]*p[X][id][ic][ri];
	  MPI_Allreduce(&loc_alpha,&alpha,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  omega=delta/alpha;
	  
	  double loc_lambda=0;
	  for(int X=0;X<loc_volh;X++)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		for(int ri=0;ri<2;ri++)
		  {
		    sol[X][id][ic][ri]+=omega*p[X][id][ic][ri];    //sol_(k+1)=x_k+omega*p_k
		    double c1=r[X][id][ic][ri]-omega*s[X][id][ic][ri];//r_(k+1)=x_k-omega*pk
		    r[X][id][ic][ri]=c1;
		    loc_lambda+=c1*c1;
		  }
	  MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  
	  double gammag=lambda/delta; //(r_(k+1),r_(k+1))/(r_k,r_k)
	  delta=lambda;
	  
	  //p_(k+1)=r_(k+1)+gammag*p_k
	  for(int X=0;X<loc_volh;X++)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		for(int ri=0;ri<2;ri++)
		  p[X][id][ic][ri]=r[X][id][ic][ri]+gammag*p[X][id][ic][ri];
	  
	  iter++;

          if(iter%10==0) master_printf("iter %d relative residue: %lg\n",iter,lambda/source_norm);
	}
      while(lambda>(residue*source_norm) && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      tmDkern_eoprec_square_eos(s,temp1,temp2,sol,conf,kappa,mu);
      {
        double loc_lambda=0;
	for(int X=0;X<loc_volh;X++)
	  for(int id=0;id<4;id++)
	    for(int ic=0;ic<3;ic++)
	      for(int ri=0;ri<2;ri++)
		{
		  double c1=source[X][id][ic][ri]-s[X][id][ic][ri];
		  loc_lambda+=c1*c1;
		}
	if(rank_tot>0) MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else lambda=loc_lambda;
      }
      master_printf("\nfinal relative residue (after %d iters): %lg where %lg was required\n",iter,lambda/source_norm,residue);

      riter++;
    }
  while(lambda>(residue*source_norm) && riter<rniter);

  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
  nissa_free(temp1);
  nissa_free(temp2);
}

//Invert twisted mass operator using e/o preconditioning.
void inv_tmD_cg_eoprec_eos(spincolor *solution_lx,spincolor *source_lx,spincolor *guess_Koo,quad_su3 *conf_lx,double kappa,double mu,int nitermax,double residue)
{
  //prepare the e/o split version of the source
  spincolor *source_eos[2];
  source_eos[0]=nissa_malloc("source_eos",loc_vol+loc_bord,spincolor);
  source_eos[1]=source_eos[0]+loc_volh+loc_bordh;
  split_lx_spincolor_into_eo_parts(source_eos,source_lx,loc_vol);
  
  //prepare the e/o split version of the solution
  spincolor *solution_eos[2];
  solution_eos[0]=nissa_malloc("solution_eos",loc_vol+loc_bord,spincolor);
  solution_eos[1]=solution_eos[0]+loc_volh+loc_bordh;
  
  //prepare the e/o split version of the conf
  quad_su3 *conf_eos[2];
  conf_eos[0]=nissa_malloc("conf_eos",loc_vol+loc_bord,quad_su3);
  conf_eos[1]=conf_eos[0]+loc_volh+loc_bordh;
  communicate_lx_gauge_borders(conf_lx);
  split_lx_conf_into_eo_parts(conf_eos,conf_lx,loc_vol+loc_bord);
  
  ///////////////////////////////////// invert with e/o improvement ///////////////////////////////////
  
  spincolor *varphi=nissa_malloc("varphi",loc_volh+loc_bordh,spincolor);
  
  //Equation (8.a)
  spincolor *temp=nissa_malloc("temp",loc_volh+loc_bordh,spincolor);
  inv_tmDee_or_oo_eos(temp,source_eos[EVN],kappa,mu);
  
  //Equation (8.b)
  tmn2Doe_eos(varphi,temp,conf_eos);
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int id=0;id<2;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  { //gamma5 is explicitely wrote
	    varphi[ivol][id  ][ic][ri]=+source_eos[ODD][ivol][id  ][ic][ri]+varphi[ivol][id  ][ic][ri]*0.5;
	    varphi[ivol][id+2][ic][ri]=-source_eos[ODD][ivol][id+2][ic][ri]-varphi[ivol][id+2][ic][ri]*0.5;
	  }
  
  //Equation (9) using solution_eos[EVN] as temporary vector
  inv_tmDkern_eoprec_square_eos(temp,varphi,guess_Koo,conf_eos,kappa,mu,nitermax,residue);
  tmDkern_eoprec_eos(solution_eos[ODD],solution_eos[EVN],temp,conf_eos,kappa,-mu);
  if(guess_Koo!=NULL) memcpy(guess_Koo,temp,sizeof(spincolor)*loc_volh); //if a guess was passed, return new one
  nissa_free(temp);

  //Equation (10)
  tmn2Deo_eos(varphi,solution_eos[ODD],conf_eos);
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  varphi[ivol][id][ic][ri]=source_eos[EVN][ivol][id][ic][ri]+varphi[ivol][id][ic][ri]*0.5;
  inv_tmDee_or_oo_eos(solution_eos[EVN],varphi,kappa,mu);

  nissa_free(varphi);
  
  /////////////////////////// paste the e/o parts of the solution together and free ///////////////////
  
  paste_eo_parts_into_lx_spincolor(solution_lx,solution_eos[EVN],solution_eos[ODD],loc_vol);
  
  nissa_free(conf_eos[0]);
  nissa_free(source_eos[0]);
  nissa_free(solution_eos[0]);
}
