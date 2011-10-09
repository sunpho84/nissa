#pragma once

//Trace the product of gamma1 * su3spinspin1^dag * gamma2 * su3spinspin2,
void trace_g_ccss_dag_g_ccss(complex *glb_c,dirac_matr *g1,su3spinspin *s1,dirac_matr *g2,su3spinspin *s2,const int ncontr)
{
  //Allocate a contiguous memory area where to store local node results
  complex *loc_c=appretto_malloc("loc_c",ncontr*glb_size[0],complex);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;
  
  for(int icontr=0;icontr<ncontr;icontr++) 
    {
      if(debug_lvl>1 && rank==0) printf("Contraction %d/%d\n",icontr+1,ncontr);

      //Local loop
      for(int ivol=0;ivol<loc_vol;ivol++)
        {
          int glb_t=glb_coord_of_loclx[ivol][0];
          //Color loop
          for(int ic1=0;ic1<3;ic1++)
	    for(int ic2=0;ic2<3;ic2++)
	      {
		complex ctemp;
		site_trace_g_sdag_g_s(ctemp,&(g1[icontr]),s1[ivol][ic2][ic1],&(g2[icontr]),s2[ivol][ic2][ic1]);
		complex_summ(loc_c[icontr*glb_size[0]+glb_t],loc_c[icontr*glb_size[0]+glb_t],ctemp);
	      }
        }
    }
  
  //Final reduction
  if(debug_lvl>1 && rank==0) printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(debug_lvl>1 && rank==0) printf("Reduction done\n");
  
  appretto_free(loc_c);
}

//print all the passed contractions to the file
void print_contractions_to_file_ptv(FILE *fout,int ncontr,int *op1,int *op2,complex *contr,int twall,const char *tag)
{
  double norm=1;
  
  if(rank==0)
    for(int icontr=0;icontr<ncontr;icontr++)
      {
        fprintf(fout,"\n");
        print_contraction_to_file(fout,op1[icontr],op2[icontr],contr+icontr*glb_size[0],twall,tag,norm);
      }
}

