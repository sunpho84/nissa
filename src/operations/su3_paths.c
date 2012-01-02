#pragma once

/////////////////////////////////////// Complicated things /////////////////////

//square (the proto-plaquette)
/*
     
  C------D
n |      | 
u |      | 
  A--mu--B
  
The square path P_{mu,nu} is defined as U(A,mu)U(B,nu)U^(C,mu)U^(A,nu)=
=U(A,mu)U(B,nu)(U(A,nu)U^(C,mu))^=U(AB,munu)*U^(AC,numu)
*/

void squared_path(su3 square,quad_su3 *conf,int A,int mu,int nu)
{
  int B=loclx_neighup[A][mu];
  int C=loclx_neighup[A][nu];

  su3 AB,AC;

  su3_prod_su3(AB,conf[A][mu],conf[B][nu]);
  su3_prod_su3(AC,conf[A][nu],conf[C][mu]);
  su3_prod_su3_dag(square,AB,AC);
}

//This calculate the global plaquette. It's not done in a very
//efficient way, but it's ok for our scope.
double global_plaquette(quad_su3 *conf)
{
  communicate_gauge_borders(conf);
  
  su3 square;
  complex pl;
  double totlocplaq=0;
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int idir=0;idir<4;idir++)
      for(int jdir=idir+1;jdir<4;jdir++)
	{
	  squared_path(square,conf,ivol,idir,jdir);
	  su3_trace(pl,square);
	  totlocplaq+=pl[0]/3;
	}
  
  double totplaq;
  MPI_Reduce(&totlocplaq,&totplaq,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  return totplaq/glb_vol/6;
}

void Pline(su3 *Pline,quad_su3 *conf)
{
  int X,iX[4],T=loc_size[0];
  
  su3 *U0=nissa_malloc("U0",loc_size[0]+1,su3);
  su3 *plf=nissa_malloc("plf",loc_size[0]+1,su3);
  su3 *plb=nissa_malloc("plb",loc_size[0]+1,su3);
  
  MPI_Status status_f,status_b;
  int tag_f=0,tag_b=1;
  su3 aux,plf_dw_buffer,plb_up_buffer;
  
  //local loop
  for(iX[1]=0;iX[1]<loc_size[1];iX[1]++)
    for(iX[2]=0;iX[2]<loc_size[2];iX[2]++)
      for(iX[3]=0;iX[3]<loc_size[3];iX[3]++)
	{
	  //forward line
	  for(iX[0]=0;iX[0]<T;iX[0]++)
	    {
	      X=loclx_of_coord(iX);
	      unsafe_su3_hermitian(U0[iX[0]+1],conf[X][0]);
             }
            su3_copy(plf[1],U0[1]);
            for(int t=2;t<=T;t++) su3_prod_su3(plf[t],U0[t],plf[t-1]);

	    if(rank_coord[0]<nrank_dir[0]-1)
	      {
		MPI_Send(plf[T],1,MPI_SU3,rank_neighup[0],tag_f,MPI_COMM_WORLD);
		MPI_Recv(plb_up_buffer,1,MPI_SU3,rank_neighup[0],tag_b,MPI_COMM_WORLD,&status_b);
	      }

	    //backward line
            for(iX[0]=0;iX[0]<T;iX[0]++)
	      {
		X=loclx_of_coord(iX);
		su3_copy(U0[iX[0]+1],conf[X][0]);
	      }
            su3_copy(plb[T],U0[T]);
            for(int t=1;t<T;t++) su3_prod_su3(plb[T-t],U0[T-t],plb[T-t+1]);

	   if(rank_coord[0]>0)
	     {
	       MPI_Send(plb[1],1,MPI_SU3,rank_neighdw[0],tag_b,MPI_COMM_WORLD);
	       MPI_Recv(plf_dw_buffer,1,MPI_SU3,rank_neighdw[0],tag_f,MPI_COMM_WORLD,&status_f);
	     }

	    //forward P-line
            iX[0]=0;
            X=loclx_of_coord(iX);
            su3_put_to_id(Pline[X]); //if glb_t=0 -> P-line=Identity
            if(glb_coord_of_loclx[X][0]>0 && glb_coord_of_loclx[X][0]<=(glb_size[0]/2-1))
	      { //if 0<t<=T/2-1 
		su3_prod_su3(aux,Pline[X],plf_dw_buffer);
		su3_copy(Pline[X],aux);
	      }
            for(iX[0]=1;iX[0]<T;iX[0]++)
	      {
                int t=iX[0];
                X=loclx_of_coord(iX);
                //This is a "local" P-line
                if(glb_coord_of_loclx[X][0]<=(glb_size[0]/2-1))
		  {
		    su3_copy(Pline[X],plf[t]);
		    if(rank_coord[0]>0)
		      {
			su3_prod_su3(aux,Pline[X],plf_dw_buffer);
			su3_copy(Pline[X],aux);
		      }
		  }
	      }

	    //backward P-line
	   for(iX[0]=T-1;iX[0]>=0;iX[0]--)
	     {
	       int t=iX[0];
	       X=loclx_of_coord(iX);
	       if(glb_coord_of_loclx[X][0]>(glb_size[0]/2-1))
		 {
		   su3_copy(Pline[X],plb[t+1]);
		   if(rank_coord[0]<(nrank_dir[0]-1))
		     {
		       su3_prod_su3(aux,Pline[X],plb_up_buffer);
		       su3_copy(Pline[X],aux);
		     }
		 }
	     }	
	}
  
  nissa_free(U0);nissa_free(plf);nissa_free(plb);
}

void Pline_forward(su3 *Pline, quad_su3 *conf)
{
  int X,iX[4],T=loc_size[0];
  
  su3 *U0=nissa_malloc("U0",loc_size[0]+1,su3);
  su3 *plf=nissa_malloc("plf",loc_size[0]+1,su3);
  
  MPI_Status status_f;
  int tag_f=0;
  su3 aux,plf_dw_buffer;
  
  //local loop
  for(iX[1]=0;iX[1]<loc_size[1];iX[1]++)
    for(iX[2]=0;iX[2]<loc_size[2];iX[2]++)
      for(iX[3]=0;iX[3]<loc_size[3];iX[3]++)
	{
	  for(iX[0]=0;iX[0]<T;iX[0]++)
	    {
	      X=loclx_of_coord(iX);
	      unsafe_su3_hermitian(U0[iX[0]+1],conf[X][0]);
	    }
           su3_copy(plf[1],U0[1]);
           for(int t=2;t<=T;t++) su3_prod_su3(plf[t],U0[t],plf[t-1]);

           if(rank_coord[0]<nrank_dir[0]-1)
	     {
	       su3_copy(plf_dw_buffer,plf[T]);
	       MPI_Send(plf_dw_buffer,1,MPI_SU3,rank_neighup[0],tag_f,MPI_COMM_WORLD);
	     }
           if(rank_coord[0]>0) MPI_Recv(plf_dw_buffer,1,MPI_SU3,rank_neighdw[0],tag_f,MPI_COMM_WORLD,&status_f);

	   iX[0]=0;
	   X=loclx_of_coord(iX);
	   su3_put_to_id(Pline[X]);
	   if(rank_coord[0]>0)
	     {
	       su3_prod_su3(aux,Pline[X],plf_dw_buffer);
	       su3_copy(Pline[X],aux);
	     }

           for(iX[0]=1;iX[0]<T;iX[0]++)
	     {
	       X=loclx_of_coord(iX);
	       //This is a "local" P-line
	       su3_copy(Pline[X],plf[iX[0]]);
	       //Here I acummulate the global P-line 
	       if(rank_coord[0]>0)
		 {
		   su3_prod_su3(aux,Pline[X],plf_dw_buffer);
		   su3_copy(Pline[X],aux);
		 }
	     }
        }

  nissa_free(U0);nissa_free(plf);
}

void Pline_backward(su3 *Pline, quad_su3 *conf)
{
  int X,iX[4],T=loc_size[0];
  
  su3 *U0=nissa_malloc("U0",loc_size[0]+1,su3);
  su3 *plb=nissa_malloc("plb",loc_size[0]+1,su3);
  
  MPI_Status status_b;
  int tag_b=1;
  su3 aux,plb_up_buffer;
  
  //local loop
  for(iX[1]=0;iX[1]<loc_size[1];iX[1]++)
    for(iX[2]=0;iX[2]<loc_size[2];iX[2]++)
      for(iX[3]=0;iX[3]<loc_size[3];iX[3]++)
	{
	  for(iX[0]=0;iX[0]<T;iX[0]++)
	    {
	      X=loclx_of_coord(iX);
	      su3_copy(U0[iX[0]+1],conf[X][0]);
            }
            su3_copy(plb[T],U0[T]);
            for(int t=1;t<T;t++) su3_prod_su3(plb[T-t],U0[T-t],plb[T-t+1]);

           if(rank_coord[0]<nrank_dir[0]-1)
	     MPI_Recv(plb_up_buffer,1,MPI_SU3,rank_neighup[0],tag_b,MPI_COMM_WORLD,&status_b);
	   
           if(rank_coord[0]>0)
	     {
	       su3_copy(plb_up_buffer,plb[1]);
	       MPI_Send(plb_up_buffer,1,MPI_SU3,rank_neighdw[0],tag_b,MPI_COMM_WORLD);
	     }

            //backward P-line
	   for(iX[0]=T-1;iX[0]>=0;iX[0]--)
	     {
	       int t=iX[0];
	       X=loclx_of_coord(iX);
	       su3_copy(Pline[X],plb[t+1]);
	       if(rank_coord[0]<nrank_dir[0]-1)
		 {
		   su3_prod_su3(aux,Pline[X],plb_up_buffer);
		   su3_copy(Pline[X],aux);
		 }
	     }
	}

  nissa_free(U0);nissa_free(plb);
}

//This will calculate 2*a^2*ig*P_{mu,nu}
/*
  ^                   C--<-- B --<--Y 
  |                   |  2  | |  1  | 
  n                   |     | |     | 
  u                   D-->--\X/-->--A 
  |                   D--<--/X\--<--A 
  -----mu---->        |  3  | |  4  | 
  		      |     | |     | 
		      E-->-- F -->--G 
*/
void Pmunu_term(as2t_su3 *Pmunu,quad_su3 *conf)
{
  int A,B,C,D,E,F,G;
  int munu;

  su3 temp1,temp2,leaves_summ;

  for(int X=0;X<loc_vol;X++)
    {
      as2t_su3_put_to_zero(Pmunu[X]);

      munu=0;
      for(int mu=0;mu<4;mu++)
	{
	  A=loclx_neighup[X][mu];
	  D=loclx_neighdw[X][mu];
	  
	  for(int nu=mu+1;nu<4;nu++)
	    {
	      B=loclx_neighup[X][nu];
	      F=loclx_neighdw[X][nu];
	      
	      C=loclx_neighup[D][nu];
	      E=loclx_neighdw[D][nu];

	      G=loclx_neighdw[A][nu];

	      //Put to 0 the summ of the leaves
	      su3_put_to_zero(leaves_summ);
	      
	      //Leaf 1
	      su3_prod_su3(temp1,conf[X][mu],conf[A][nu]);         //    B--<--Y 
	      su3_prod_su3_dag(temp2,temp1,conf[B][mu]);           //    |  1  | 
	      su3_prod_su3_dag(temp1,temp2,conf[X][nu]);           //    |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);	           //    X-->--A 

	      //Leaf 2
	      su3_prod_su3_dag(temp1,conf[X][nu],conf[C][mu]);      //   C--<--B
	      su3_prod_su3_dag(temp2,temp1,conf[D][nu]);            //   |  2  | 
	      su3_prod_su3(temp1,temp2,conf[D][mu]);		    //   |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		    //   D-->--X
	      
	      //Leaf 3
	      su3_dag_prod_su3_dag(temp1,conf[D][mu],conf[E][nu]);  //   D--<--X
	      su3_prod_su3(temp2,temp1,conf[E][mu]);		    //   |  3  | 
	      su3_prod_su3(temp1,temp2,conf[F][nu]);		    //   |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		    //   E-->--F
	      
	      //Leaf 4
	      su3_dag_prod_su3(temp1,conf[F][nu],conf[F][mu]);       //  X--<--A 
	      su3_prod_su3(temp2,temp1,conf[G][nu]);                 //  |  4  | 
	      su3_prod_su3_dag(temp1,temp2,conf[X][mu]);             //  |     |  
	      su3_summ(leaves_summ,leaves_summ,temp1);               //  F-->--G 

	      //calculate U-U^dagger
	      for(int ic1=0;ic1<3;ic1++)
		for(int ic2=0;ic2<3;ic2++)
		  {
		    Pmunu[X][munu][ic1][ic2][0]=(leaves_summ[ic1][ic2][0]-leaves_summ[ic2][ic1][0])/4;
		    Pmunu[X][munu][ic1][ic2][1]=(leaves_summ[ic1][ic2][1]+leaves_summ[ic2][ic1][1])/4;
		  }

	      munu++;
	    }
	}
    }
}

//apply the chromo operator to the passed spinor site by site (not yet fully optimized)
void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,as2t_su3 Pmunu,spincolor in)
{
  color temp_d1;
  
  for(int d1=0;d1<4;d1++)
    {
      color_put_to_zero(out[d1]);
      for(int imunu=0;imunu<6;imunu++)
	{
	  unsafe_su3_prod_color(temp_d1,Pmunu[imunu],in[smunu_pos[d1][imunu]]);
	  for(int c=0;c<3;c++) complex_summ_the_prod(out[d1][c],smunu_entr[d1][imunu],temp_d1[c]);
	}
    }
}

//apply the chromo operator to the passed spinor to the whole volume
void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,as2t_su3 *Pmunu,spincolor *in)
{
  for(int ivol=0;ivol<loc_vol;ivol++) unsafe_apply_point_chromo_operator_to_spincolor(out[ivol],Pmunu[ivol],in[ivol]);
}

//apply the chromo operator to the passed colorspinspin
//normalization as in ape next
void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,as2t_su3 *Pmunu,colorspinspin *in)
{
  spincolor temp1,temp2;
  
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<4;id_source++) //dirac index of source
	{
	  //Switch the color_spinspin into the spincolor.
	  get_spincolor_from_colorspinspin(temp1,in[loc_site],id_source);
	  
	  unsafe_apply_point_chromo_operator_to_spincolor(temp2,Pmunu[loc_site],temp1);
	  
	  //Switch back the spincolor into the colorspinspin
	  put_spincolor_into_colorspinspin(out[loc_site],temp2,id_source);
	}
    }
}

//apply the chromo operator to the passed su3spinspin
//normalization as in ape next
void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin *out,as2t_su3 *Pmunu,su3spinspin *in)
{
  spincolor temp1,temp2;
  
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<4;id_source++) //dirac index of source
	for(int ic_source=0;ic_source<3;ic_source++) //color index of source
	  {
	    //Switch the su3spinspin into the spincolor.
	    get_spincolor_from_su3spinspin(temp1,in[loc_site],id_source,ic_source);
	    
	    unsafe_apply_point_chromo_operator_to_spincolor(temp2,Pmunu[loc_site],temp1);
	    
	    //Switch back the spincolor into the colorspinspin
	    put_spincolor_into_su3spinspin(out[loc_site],temp2,id_source,ic_source);
	  }
    }
}
