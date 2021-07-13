#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_GEOMETRY_LEB
 #include "geometry_Leb.hpp"

#include "geometry_eo.hpp"

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"

namespace nissa
{
  //return the loclx coordinate of a Leblx
  coords_t loclx_coord_of_Leblx(const Leb_factors_t &factors,int Leblx)
  {
    coords_t c;
    
    int nfactors=factors.size();
    coords_t t1,t2,t3;
    
    //init
    for(int i=0;i<NDIM;i++)
      {
 	t1[i]=t2[i]=c[i]=0;
	t3[i]=(NDIM*nfactors-1)%nfactors;
      }
    
    //convert to mixed base
    int x_mixed_base[NDIM][NDIM*nfactors];
    memset(x_mixed_base,0,sizeof(x_mixed_base));
    for(int i=0;i<NDIM*nfactors;i++)
      {
	int mu=NDIM-1-(i%NDIM);
	int f=factors[t1[mu]][mu];
	
	x_mixed_base[mu][t2[mu]]=Leblx%f;
	Leblx/=f;
	
	t1[mu]=(t1[mu]+1)%nfactors;
	t2[mu]++;
      }
    
    //build coordinate in lx format
    for(int mu=0;mu<NDIM;mu++)
      for(int j=NDIM*nfactors-1;j>=0;j--)
	{
	  c[mu]=x_mixed_base[mu][j]+factors[t3[mu]][mu]*c[mu];
	  t3[mu]=(t3[mu]+nfactors-1)%nfactors;
	}
    
    return c;
  }
  
  //init the Lebesgue geometry
  void set_Leb_geometry()
  {
    if(not use_Leb_geom) crash("Lebesgue Geometry was not to be used!");
    if(Leb_geom_inited) crash("Lebesgue Geometry already initialized!");
    Leb_geom_inited=true;
    
    loclx_of_Leblx=nissa_malloc("loclx_of_Leblx",locVol+bord_vol+edge_vol,int);
    Leblx_of_loclx=nissa_malloc("Leblx_of_loclx",locVol+bord_vol+edge_vol,int);
    Leblx_neighup=nissa_malloc("Leblx_neighup",locVol,coords_t);
    Leblx_neighdw=nissa_malloc("Leblx_neighdw",locVol,coords_t);
    Leblx_parity=nissa_malloc("Leblx_parity",locVol+bord_vol+edge_vol,int);
    surfLeblx_of_bordLeblx=nissa_malloc("surfLeblx_of_Bordlx",bord_vol,int);
    
    //get nmax_fact
    int nmax_facts=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	int list_fact_mu[log2N(locSize[mu])];
	nmax_facts=std::max(nmax_facts,factorize(list_fact_mu,locSize[mu]));
      }
    
    //set all factors to 1
    Leb_factors_t factors(nmax_facts,std::vector<int>(NDIM));
    for(int i=0;i<nmax_facts;i++) for(int mu=0;mu<NDIM;mu++) factors[i][mu]=1;
    
    //put all the non-1 factors
    for(int mu=0;mu<NDIM;mu++)
      {
	int list_fact_mu[log2N(locSize[mu])];
	int nfacts=factorize(list_fact_mu,locSize[mu]);
	int nfacts1=nmax_facts-nfacts;
	for(int ifact=0;ifact<nfacts;ifact++) factors[nfacts1+ifact][mu]=list_fact_mu[ifact];
      }
    
    //print the factors
    verbosity_lv3_master_printf("Leb factors\n");
    for(int mu=0;mu<NDIM;mu++)
      {
	for(int ifact=0;ifact<nmax_facts;ifact++) verbosity_lv3_master_printf("%d ",factors[ifact][mu]);
	verbosity_lv3_master_printf("\n");
      }
    
    //on a first basis, everything's the same
    for(int i=0;i<locVol+bord_vol+edge_vol;i++)
      Leblx_of_loclx[i]=loclx_of_Leblx[i]=i;
    
    //fill Leblx of loclx and the opposite
    for(int Leblx=0;Leblx<locVol;Leblx++)
      {
	const coords_t c=loclx_coord_of_Leblx(factors,Leblx);
	int loclx=loclx_of_coord(c);
	loclx_of_Leblx[Leblx]=loclx;
	Leblx_of_loclx[loclx]=Leblx;
	
	verbosity_lv3_master_printf("%d %d\n",Leblx_of_loclx[loclx],loclx_of_Leblx[Leblx]);
      }
    
    //set movements
    for(int Leblx=0;Leblx<locVol;Leblx++)
      {
	int loclx=loclx_of_Leblx[Leblx];
	for(int mu=0;mu<NDIM;mu++)
	  {
	    Leblx_neighup[Leblx][mu]=Leblx_of_loclx[loclxNeighup[loclx][mu]];
	    Leblx_neighdw[Leblx][mu]=Leblx_of_loclx[loclxNeighdw[loclx][mu]];
	  }
      }
    
    //set surf filler
    for(int ibord=0;ibord<bord_vol;ibord++) surfLeblx_of_bordLeblx[ibord]=Leblx_of_loclx[surflxOfBordlx[ibord]];
    
    //set e/o
    if(use_eo_geom)
      {
	for(int eo=0;eo<2;eo++)
	  {
	    loceo_of_Lebeo[eo]=nissa_malloc("loceo_of_Lebeo",locVolh+bord_volh+edge_volh,int);
	    Lebeo_of_loceo[eo]=nissa_malloc("Lebeo_of_loceo",locVolh+bord_volh+edge_volh,int);
	    Lebeo_neighup[eo]=nissa_malloc("Lebeo_neighup",locVolh,coords_t);
	    Lebeo_neighdw[eo]=nissa_malloc("Lebeo_neighdw",locVolh,coords_t);
	  }
	
	int iLebeo[2]={0,0};
	for(int Leblx=0;Leblx<locVol+bord_vol+edge_vol;Leblx++)
	  {
	    int loclx=loclx_of_Leblx[Leblx];
	    int par=loclx_parity[loclx];
	    Leblx_parity[Leblx]=par;
	    int loceo=loceo_of_loclx[loclx];
	    
	    //set Lebeo of loceo and opposite
	    int Lebeo=iLebeo[par]++;
	    Lebeo_of_loceo[par][loceo]=Lebeo;
	    loceo_of_Lebeo[par][Lebeo]=loceo;
	  }
	
	//set e/o movements
	for(int Leblx=0;Leblx<locVol;Leblx++)
	  {
	    int loclx=loclx_of_Leblx[Leblx];
	    int par=loclx_parity[loclx];
	    int loceo=loceo_of_loclx[loclx];
	    int Lebeo=Lebeo_of_loceo[par][loceo];
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		Lebeo_neighup[par][Lebeo][mu]=Lebeo_of_loceo[!par][loceo_neighup[par][loceo][mu]];
		Lebeo_neighdw[par][Lebeo][mu]=Lebeo_of_loceo[!par][loceo_neighdw[par][loceo][mu]];
	      }
	  }
      }
  }
  
  //unset the Lebesgue geometry
  void unset_Leb_geometry()
  {
    if(not Leb_geom_inited) crash("asking to unset never initialized Lebesgue Geometry!");
    
    nissa_free(loclx_of_Leblx);
    nissa_free(Leblx_of_loclx);
    nissa_free(Leblx_neighup);
    nissa_free(Leblx_neighdw);
    nissa_free(Leblx_parity);
    nissa_free(surfLeblx_of_bordLeblx);
    
    if(use_eo_geom)
      for(int eo=0;eo<2;eo++)
	{
	  nissa_free(loceo_of_Lebeo[eo]);
	  nissa_free(Lebeo_of_loceo[eo]);
	  nissa_free(Lebeo_neighup[eo]);
	  nissa_free(Lebeo_neighdw[eo]);
	}
  }
}
