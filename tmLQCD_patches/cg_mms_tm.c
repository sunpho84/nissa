/***********************************************************************
 *
 * $Id: cg_mms_tm.c 1460 2010-02-03 13:20:45Z urbach $
 *
 * Copyright (C) 2004 Andrea Shindler
 *               2009 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Author: Andrea Shindler <shindler@ifh.de> Jan 2004
 * 
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include <io/spinor.h>
#include <io/params.h>
#include "gamma.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver/matrix_mult_typedef.h"
#include "cg_mms_tm.h"
#include <io/params.h>

static spinor * xs_qmms;
static spinor * ps_qmms;
static spinor ** xs_mms_solver;
static spinor ** ps_mms_solver;
static double * sigma;
static double * zitam1, * zita;
static double * alphas, * betas;


void init_mms_tm(const int nr);


/* P output = solution , Q input = source */
int cg_mms_tm(spinor * const P, spinor * const Q, const int max_iter, 
	      double eps_sq, const int rel_prec, const int N, matrix_mult f) {

  static double normsq, pro, err, alpha_cg = 1., beta_cg = 0., normsp, squarenorm;
  int iteration, im;
  char filename[100];
  static double gamma,alpham1;
  
  double tmp_mu = g_mu;
  WRITER * writer;
  paramsPropagatorFormat *propagatorFormat = NULL;
  
  init_mms_tm(g_no_extra_masses);

  /* currently only implemented for P=0 */
  zero_spinor_field(P, N);
  /*  Value of the bare MMS-masses (\mu^2 - \mu_0^2) */
  for(im = 0; im < g_no_extra_masses; im++) {
    sigma[im] = g_extra_masses[im]*g_extra_masses[im] - g_mu*g_mu;
    assign(xs_mms_solver[im], P, N);
    assign(ps_mms_solver[im], Q, N);
    zitam1[im] = 1.0;
    zita[im] = 1.0;
    alphas[im] = 1.0;
    betas[im] = 0.0;
  }
 
  squarenorm = square_norm(Q, N, 1);
  assign(g_spinor_field[DUM_SOLVER], P, N);
  normsp = square_norm(P, N, 1);
  assign(g_spinor_field[DUM_SOLVER+5], Q, N);

  /* initialize residue r and search vector p */
/*   if(normsp == 0){ */
  /* currently only implemented for P=0 */
  if(1) {
    /* if a starting solution vector equal to zero is chosen */
    assign(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+5], N);
    assign(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+5], N);
    normsq = square_norm(Q, N, 1);
  }
  else{
    /* if a starting solution vector different from zero is chosen */
    f(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER]);
   
    diff(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+5], g_spinor_field[DUM_SOLVER+3], N);
    assign(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], N);
    normsq = square_norm(g_spinor_field[DUM_SOLVER+2], N, 1);
  }

  /* main loop */
  for(iteration = 0; iteration < max_iter; iteration++) {
    
    /*   Q^2*p and then (p,Q^2*p)  */
    f(g_spinor_field[DUM_SOLVER+4], g_spinor_field[DUM_SOLVER+2]);
    pro = scalar_prod_r(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+4], N, 1);
    
    /* For the update of the coeff. of the shifted pol. we need alpha_cg(i-1) and alpha_cg(i).
       This is the reason why we need this double definition of alpha */
    alpham1 = alpha_cg;

    /* Compute alpha_cg(i+1) */
    alpha_cg = normsq/pro;
    for(im = 0; im < g_no_extra_masses; im++) {
      
      /* Now gamma is a temp variable that corresponds to zita(i+1) */ 
      gamma = zita[im]*alpham1/(alpha_cg*beta_cg*(1.-zita[im]/zitam1[im]) 
				+ alpham1*(1.+sigma[im]*alpha_cg));
      
      /* Now zita(i-1) is put equal to the old zita(i) */
      zitam1[im] = zita[im];
      /* Now zita(i+1) is updated */
      zita[im] = gamma;
      /* Update of alphas(i) = alpha_cg(i)*zita(i+1)/zita(i) */ 
      alphas[im] = alpha_cg*zita[im]/zitam1[im];
      /* Compute xs(i+1) = xs(i) + alphas(i)*ps(i) */
      assign_add_mul_r(xs_mms_solver[im], ps_mms_solver[im], alphas[im], N); 
    } 
    
    /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
    assign_add_mul_r(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+2],  alpha_cg, N);
    /*  Compute r_(i+1) = r_i - alpha_cg(i+1) Qp_i   */
    assign_add_mul_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+4], -alpha_cg, N);

    /* Check whether the precision eps_sq is reached */
    
    err = square_norm(g_spinor_field[DUM_SOLVER+1], N, 1);
    if(g_debug_level > 0 && g_proc_id == g_stdio_proc) {
      printf("CG MMS %d\t%g\n", iteration, err); fflush( stdout );
    }
    
    if( ((err <= eps_sq) && (rel_prec == 0)) ||
	((err <= eps_sq*squarenorm) && (rel_prec == 1)) ) {

      assign(P, g_spinor_field[DUM_SOLVER], N);
      f(g_spinor_field[DUM_SOLVER+2], P);
      diff(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], Q, N);
      err = square_norm(g_spinor_field[DUM_SOLVER+3], N, 1);
      if(g_debug_level > 0 && g_proc_id == g_stdio_proc) {
	printf("true residue %d\t%g\t\n",iteration, err); 
	fflush( stdout);
      }
      g_sloppy_precision = 0;
      g_mu = tmp_mu;

      //Try to open the file which instruct what to save
      FILE* instruct;
      int saveDD,save0,save1,nbits;
      save0=save1=1;
      nbits=32;
      if(g_proc_id==0)
	{
	  instruct=fopen("cg_mms_instructions","r");
	  if(instruct!=NULL)
	    {
	      int nintread=0;
	      nintread+=fscanf(instruct,"%d",&save0);
	      nintread+=fscanf(instruct,"%d",&save1);
	      nintread+=fscanf(instruct,"%d",&nbits);
	      nintread+=fscanf(instruct,"%d",&saveDD);
	      if(nintread!=4)
		{
		  printf("Error in reading instructions for the cgmms save\n");
		  saveDD=1;
		  save0=1;
		  save1=1;
		  nbits=32;
		}
	      if(nbits!=32) nbits=64;
	      fclose(instruct);
	    }
	  else printf("No instruction for the cgmms save, saving both props in prec 32 bits\n");
	}
      
#ifdef MPI
      MPI_Bcast(&save0,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&save1,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nbits,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
      
      /* save all the results of (Q^dagger Q)^(-1) \gamma_5 \phi */
      /* here ... */
      if(saveDD==1)
	{
	  sprintf(filename,"%s.%.4d.%.2d.%.2d.cgmms.%.2d.inverted", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, 0);
	  if(g_kappa != 0) {
	    mul_r(g_spinor_field[DUM_SOLVER+1], (2*g_kappa)*(2*g_kappa), g_spinor_field[DUM_SOLVER], N);
	  }
	  
	  // the 0 is for appending 
	  construct_writer(&writer, filename, 0);
	  
	  propagatorFormat = construct_paramsPropagatorFormat(nbits, 1);
	  write_propagator_format(writer, propagatorFormat);
	  free(propagatorFormat);
	  
	  convert_lexic_to_eo(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], 
			      g_spinor_field[DUM_SOLVER+1]);
	  write_spinor(writer, &g_spinor_field[DUM_SOLVER+3], &g_spinor_field[DUM_SOLVER+2], 1, nbits);
	  destruct_writer(writer);
	  
	}
      
      //Now we save the 'up' and 'down' propagators

      mul_r(g_spinor_field[DUM_SOLVER], 2*g_kappa, g_spinor_field[DUM_SOLVER], N);
      
      //up
      if(save0==1)
	{
	  Q_plus_psi(g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER]);
	  convert_lexic_to_eo(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], 
			      g_spinor_field[DUM_SOLVER+1]);
	  sprintf(filename,"%s.%.4d.%.2d.%.2d.cgmms.%.2d.inverted.0", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, 0);
	  construct_writer(&writer, filename, 1);
	  propagatorFormat = construct_paramsPropagatorFormat(nbits, 1);
	  write_propagator_format(writer, propagatorFormat);
	  free(propagatorFormat);
	  write_spinor(writer, &g_spinor_field[DUM_SOLVER+3], &g_spinor_field[DUM_SOLVER+2], 1, nbits);
	  destruct_writer(writer);
	}
      //down
      if(save1==1)
	{
	  Q_minus_psi(g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER]);
	  convert_lexic_to_eo(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2],
			      g_spinor_field[DUM_SOLVER+1]);
	  sprintf(filename,"%s.%.4d.%.2d.%.2d.cgmms.%.2d.inverted.1", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, 0);
	  construct_writer(&writer, filename, 1);
	  propagatorFormat = construct_paramsPropagatorFormat(nbits, 1);
	  write_propagator_format(writer, propagatorFormat);
	  free(propagatorFormat);
	  write_spinor(writer, &g_spinor_field[DUM_SOLVER+3], &g_spinor_field[DUM_SOLVER+2], 1, nbits);
	  destruct_writer(writer);
	}

      for(im = 0; im < g_no_extra_masses; im++) {
	
	if(saveDD)
	  {
	    g_mu = g_extra_masses[im];
	    
	    
	    sprintf(filename,"%s.%.4d.%.2d.%.2d.cgmms.%.2d.inverted", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, im+1);
	    
	    construct_writer(&writer, filename, 0);
	    
	    propagatorFormat = construct_paramsPropagatorFormat(nbits, 1);
	    write_propagator_format(writer, propagatorFormat);
	    free(propagatorFormat);
	    
	    if(g_kappa != 0) {
	      mul_r(g_spinor_field[DUM_SOLVER+1], (2*g_kappa)*(2*g_kappa), xs_mms_solver[im], N);
	    }
	    convert_lexic_to_eo(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1]);
	    
	    write_spinor(writer, &g_spinor_field[DUM_SOLVER+3], &g_spinor_field[DUM_SOLVER+2], 1, nbits);
	    destruct_writer(writer);
	  }
	
	//Now we save the 'up' and 'down' propagators
	mul_r(xs_mms_solver[im], 2*g_kappa,xs_mms_solver[im] , N);
	
	//up
	if(save0==1)
	  {
	    Q_plus_psi(g_spinor_field[DUM_SOLVER+1],xs_mms_solver[im]);
	    convert_lexic_to_eo(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], 
				g_spinor_field[DUM_SOLVER+1]);
	    sprintf(filename,"%s.%.4d.%.2d.%.2d.cgmms.%.2d.inverted.0", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, im+1);
	    construct_writer(&writer, filename, 1);
	    propagatorFormat = construct_paramsPropagatorFormat(nbits, 1);
	    write_propagator_format(writer, propagatorFormat);
	    free(propagatorFormat);
	    write_spinor(writer, &g_spinor_field[DUM_SOLVER+3], &g_spinor_field[DUM_SOLVER+2], 1, nbits);
	    destruct_writer(writer);
	  }
	//down
	if(save1==1)
	  {
	    Q_minus_psi(g_spinor_field[DUM_SOLVER+1],xs_mms_solver[im]);
	    convert_lexic_to_eo(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2], 
				g_spinor_field[DUM_SOLVER+1]);
	    sprintf(filename,"%s.%.4d.%.2d.%.2d.cgmms.%.2d.inverted.1", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, im+1);
	    construct_writer(&writer, filename, 1);
	    propagatorFormat = construct_paramsPropagatorFormat(nbits, 1);
	    write_propagator_format(writer, propagatorFormat);
	    free(propagatorFormat);
	    write_spinor(writer, &g_spinor_field[DUM_SOLVER+3], &g_spinor_field[DUM_SOLVER+2], 1, nbits);
	    destruct_writer(writer);	
	  }
      }

      g_mu = tmp_mu;

      return(iteration+1);
    }
    
    /* Compute beta_cg(i+1) = (r(i+1),r(i+1))/(r(i),r(i))
       Compute p(i+1) = r(i+1) + beta(i+1)*p(i)  */
    beta_cg = err/normsq;
    assign_mul_add_r(g_spinor_field[DUM_SOLVER+2], beta_cg, g_spinor_field[DUM_SOLVER+1], N);
    normsq = err;
    
    /* Compute betas(i+1) = beta_cg(i)*(zita(i+1)*alphas(i))/(zita(i)*alpha_cg(i))
       Compute ps(i+1) = zita(i+1)*r(i+1) + betas(i+1)*ps(i)  */
    for(im = 0; im < g_no_extra_masses; im++) {
      betas[im] = beta_cg*zita[im]*alphas[im]/(zitam1[im]*alpha_cg);
      assign_mul_add_mul_r(ps_mms_solver[im], g_spinor_field[DUM_SOLVER+1], betas[im], zita[im], N);
    }
  }
  assign(P, g_spinor_field[DUM_SOLVER], N);
  g_sloppy_precision = 0;
  return(-1);
}


void init_mms_tm(const int nr) {
  static int ini = 0;
  int i;
  if(ini == 0) {

    sigma = (double*)calloc((nr), sizeof(double));
    zitam1 = (double*)calloc((nr), sizeof(double));
    zita = (double*)calloc((nr), sizeof(double));
    alphas = (double*)calloc((nr), sizeof(double));
    betas = (double*)calloc((nr), sizeof(double));

#if (defined SSE2 || defined SSE)
    xs_qmms = (spinor*)calloc(VOLUMEPLUSRAND*(nr)+1,sizeof(spinor));
    xs_mms_solver = (spinor**)calloc((nr)+1,sizeof(spinor*));
    
    ps_qmms = (spinor*)calloc(VOLUMEPLUSRAND*(nr)+1,sizeof(spinor));
    ps_mms_solver = (spinor**)calloc((nr)+1,sizeof(spinor*));

    for(i = 0; i < nr; i++) {
      xs_mms_solver[i]=(spinor*)(((unsigned long int)(xs_qmms)+ALIGN_BASE)&~ALIGN_BASE) + i*VOLUMEPLUSRAND;
      ps_mms_solver[i]=(spinor*)(((unsigned long int)(ps_qmms)+ALIGN_BASE)&~ALIGN_BASE) + i*VOLUMEPLUSRAND;
    }
#else
    xs_qmms = (spinor*)calloc(VOLUMEPLUSRAND*(nr),sizeof(spinor));
    xs_mms_solver = (spinor**)calloc((nr),sizeof(spinor*));

    ps_qmms = (spinor*)calloc(VOLUMEPLUSRAND*(nr),sizeof(spinor));
    ps_mms_solver = (spinor**)calloc((nr),sizeof(spinor*));

    for(i = 0; i < nr; i++) {
      xs_mms_solver[i] = xs_qmms + i*VOLUMEPLUSRAND;
      ps_mms_solver[i] = ps_qmms + i*VOLUMEPLUSRAND;
    }
#endif
    ini=1;
  }
}
