#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/complex.hpp"
#include "new_types/new_types_definitions.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{  
  double rnd_get_unif(rnd_gen *gen,double min,double max);
  
  //initialize a random number generator
  void start_rnd_gen(rnd_gen *out,int seed)
  {
    const int im1=2147483563,ia1=40014;
    const int iq1=53668,ir1=12211;
    int j,k;
    
    //initialization
    out->idum=seed;
    out->idum=std::max(out->idum+1,1);
    out->idum2=out->idum;
    for(j=RAN2_NTAB+7;j>=0;j--)
      {
	k=out->idum/iq1;
	out->idum=ia1*(out->idum-k*iq1)-k*ir1;
	if(out->idum<0) out->idum+=im1;
	if(j<RAN2_NTAB) out->iv[j]=out->idum;
      }
    out->iy=out->iv[0];
  }
  
  //print all the entries of the random generator into a string
  void convert_rnd_gen_to_text(char *text,rnd_gen *gen)
  {
    int *ptr=(int*)gen;
    
    //init output text
    sprintf(text,"%d",ptr[0]);
    
    //append all the elements
    for(int i=1;i<RAN2_NTAB+3;i++)
      {
	char temp[20];
	sprintf(temp," %d",ptr[i]);
	strcat(text,temp);
      }
  }
  
  //read all the entries of the random generator from a string
  void convert_text_to_rnd_gen(rnd_gen *gen,char *text)
  {
    int *ptr=(int*)gen;
    
    for(int i=0;i<RAN2_NTAB+3;i++)
      {
	char temp[20];
	if(sscanf(text,"%s",temp)!=1) crash("while reading element %d from %s",i,text);
	text+=1+strlen(temp);
	if(sscanf(temp,"%d",ptr+i)!=1) crash("while converting to int %s",temp);
      }
  }
  
  //initialize the global random generator
  void start_glb_rnd_gen(int seed)
  {
    if(glb_rnd_gen_inited==1) crash("global random generator already initialized");
    start_rnd_gen(&(glb_rnd_gen),seed);
    
    glb_rnd_gen_inited=1;
    master_printf("Global random generators initialized with seed: %d\n",seed);
  }
  
  //init from text
  void start_glb_rnd_gen(char *text)
  {
    if(glb_rnd_gen_inited==1) crash("global random generator already initialized");
    convert_text_to_rnd_gen(&(glb_rnd_gen),text);
    
    glb_rnd_gen_inited=1;
    master_printf("Global random generators initialized from text\n");
  }
  
  //initialize the grid of local random number generator
  void start_loc_rnd_gen(int seed)
  {
    if(loc_rnd_gen_inited==1) crash("local random generator already initialized!");
    
    //check the grid to be initiaized
    if(loc_vol==0) crash("grid not initalized!");
    
    //Generate the true seed
    if(glb_rnd_gen_inited==0) start_glb_rnd_gen(seed);
    int internal_seed=(int)rnd_get_unif(&glb_rnd_gen,0,RAND_MAX);
    
    //allocate the grid of random generator, one for site
    loc_rnd_gen=nissa_malloc("Loc_rnd_gen",loc_vol,rnd_gen);
    for(int ivol=0;ivol<loc_vol;ivol++) start_rnd_gen(&(loc_rnd_gen[ivol]),internal_seed+glblx_of_loclx[ivol]);
    
    loc_rnd_gen_inited=1;
    master_printf("Grid of local random generators initialized with internal seed: %d\n",internal_seed);
  }
  
  //init from text
  void start_loc_rnd_gen(char *text)
  {
    start_glb_rnd_gen(text);
    start_loc_rnd_gen(0); //glb grid already started so does not matter
  }
  
  //stop grid of local random generators
  void stop_loc_rnd_gen()
  {
    if(loc_rnd_gen_inited==0) crash("local random generator not initialized");
    master_printf("Stopping local random generators\n");
    
    nissa_free(loc_rnd_gen);
    loc_rnd_gen_inited=0;
  }
  
  //standard ran2 from numerical recipes
  double rnd_get_unif(rnd_gen *gen,double min,double max)
  {
    const int im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,ia2=40692;
    const int iq1=53668,iq2=52774,ir1=12211,ir2=3791,ndiv=1+imm1/RAN2_NTAB;
    const double am=1.0/im1,eps=1.2e-7,rnmx=1-eps;
    int j,k;
    double out;
    
    k=gen->idum/iq1;
    gen->idum=ia1*(gen->idum-k*iq1)-k*ir1;
    if(gen->idum<0) gen->idum+=im1;
    
    k=gen->idum2/iq2;
    gen->idum2=ia2*(gen->idum2-k*iq2)-k*ir2;
    if(gen->idum2<0) gen->idum2+=im2;
    
    j=gen->iy/ndiv;
    gen->iy=gen->iv[j]-gen->idum2;
    gen->iv[j]=gen->idum;
    if(gen->iy<0) gen->iy+=imm1;
    
    out=std::min(am*gen->iy,rnmx);
    
    return out*(max-min)+min;
  }
  
  //return a numer between 0 and 1
  int rnd_get_pm_one(rnd_gen *gen)
  {
    double r=rnd_get_unif(gen,0,1);
    if(r>0.5) return 1;
    else return -1;
  }
  
  //return a Z2 complex
  void rnd_get_Z2(complex out,rnd_gen *gen)
  {
    out[0]=rnd_get_pm_one(gen);
    out[1]=0;
  }
  
  //return a Z4 complex
  void rnd_get_Z4(complex out,rnd_gen *gen)
  {
    out[0]=rnd_get_pm_one(gen)/RAD2;
    out[1]=rnd_get_pm_one(gen)/RAD2;
  }
  
  //return a gaussian complex
  void rnd_get_gauss_complex(complex out,rnd_gen *gen,complex ave,double sig)
  {
    const double one_by_sqrt2=0.707106781186547;
    double norm=sig*one_by_sqrt2;
    double q,r;
    
    r=sqrt(-2*log(1-rnd_get_unif(gen,0,1)));
    q=2*M_PI*rnd_get_unif(gen,0,1);
    
    out[0]=r*cos(q)*norm+ave[0];
    out[1]=r*sin(q)*norm+ave[1];
  }
  
  //return a complex number of appropriate type
  void comp_get_rnd(complex out,rnd_gen *gen,enum rnd_t rtype)
  {
    complex z={0,0};
    switch(rtype)
      {
      case RND_ALL_PLUS_ONE:  out[0]=1;                     out[1]=0;break;
      case RND_ALL_MINUS_ONE: out[0]=-1;                    out[1]=0;break;
      case RND_UNIF:          out[0]=rnd_get_unif(gen,0,1); out[1]=0;break;
      case RND_Z2:            rnd_get_Z2(out,gen);                   break;
      case RND_Z4:            rnd_get_Z4(out,gen);                   break;
      case RND_GAUSS:         rnd_get_gauss_complex(out,gen,z,1);    break;
      }
  }
  
  //fill a grid of vectors with numbers between 0 and 1
  THREADABLE_FUNCTION_4ARG(rnd_fill_unif_loc_vector, double*,v, int,dps, double,min, double,max)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int i=0;i<dps;i++)
	v[ivol*dps+i]=rnd_get_unif(&(loc_rnd_gen[ivol]),min,max);
    
    set_borders_invalid(v);
  }
  THREADABLE_FUNCTION_END
  
  //return a grid of +-x numbers
  THREADABLE_FUNCTION_2ARG(rnd_fill_pm_one_loc_vector, double*,v, int,nps)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int i=0;i<nps;i++)
	v[ivol*nps+i]=rnd_get_pm_one(&(loc_rnd_gen[ivol]));
    
    set_borders_invalid(v);
  }
  THREADABLE_FUNCTION_END
  
  //generate a spindiluted vector according to the passed type
  THREADABLE_FUNCTION_3ARG(generate_spindiluted_source, colorspinspin*,source, enum rnd_t,rtype, int,twall)
  {
    //reset
    vector_reset(source);
    
    //compute normalization norm: spat vol if twall>=0, glb_vol else
    int norm2=glb_vol;
    if(twall>=0) norm2/=glb_size[0];
    double inv_sqrt_norm2=1.0/sqrt(norm2);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      if(glb_coord_of_loclx[ivol][0]==twall||twall<0)
	for(int ic=0;ic<3;ic++)
	  {
	    //generate the id_sink==id_source==0 entry
	    comp_get_rnd(source[ivol][ic][0][0],&(loc_rnd_gen[ivol]),rtype);
	    complex_prod_double(source[ivol][ic][0][0],source[ivol][ic][0][0],inv_sqrt_norm2);
	    //copy the other three dirac indexes
	    for(int d=1;d<4;d++)
	      memcpy(source[ivol][ic][d][d],source[ivol][ic][0][0],sizeof(complex));
	  }
    
    set_borders_invalid(source);
  }
  THREADABLE_FUNCTION_END
  
  //generate an undiluted vector according to the passed type
  THREADABLE_FUNCTION_3ARG(generate_undiluted_source, spincolor*,source, enum rnd_t,rtype, int,twall)
  {
    //reset
    vector_reset(source);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      if(glb_coord_of_loclx[ivol][0]==twall||twall<0)
	for(int id=0;id<4;id++)
	  for(int ic=0;ic<3;ic++)
	    comp_get_rnd(source[ivol][id][ic],&(loc_rnd_gen[ivol]),rtype);
    
    set_borders_invalid(source);
  }
  THREADABLE_FUNCTION_END
  
  //generate a fully undiluted source
  THREADABLE_FUNCTION_5ARG(generate_fully_undiluted_eo_source, color*,source, enum rnd_t,rtype, int,twall, int,par, int,dir)
  {
    vector_reset(source);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
      {
	int ilx=loclx_of_loceo[par][ieo];
	if(glb_coord_of_loclx[ilx][dir]==twall||twall<0)
	  for(int ic=0;ic<3;ic++)
	    comp_get_rnd(source[ieo][ic],&(loc_rnd_gen[ilx]),rtype);
      }
    
    set_borders_invalid(source);
  }
  THREADABLE_FUNCTION_END
  void generate_fully_undiluted_eo_source(color **source,enum rnd_t rtype,int twall,int dir=0)
  {for(int par=0;par<2;par++) generate_fully_undiluted_eo_source(source[par],rtype,twall,par,dir);}
  
  //generate a delta source
  THREADABLE_FUNCTION_2ARG(generate_delta_source, su3spinspin*,source, int*,x)
  {
    //reset
    vector_reset(source);
    
    int islocal=1;
    coords lx;
    for(int mu=0;mu<NDIM;mu++)
      {
	lx[mu]=x[mu]-rank_coord[mu]*loc_size[mu];
	islocal&=(lx[mu]>=0);
	islocal&=(lx[mu]<loc_size[mu]);
      }
    
    if(islocal)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  source[loclx_of_coord(lx)][ic][ic][id][id][0]=1;
    
    set_borders_invalid(source);
  }
  THREADABLE_FUNCTION_END

  //generate a delta source
  THREADABLE_FUNCTION_2ARG(generate_delta_eo_source, color**,source, int*,x)
  {
    //reset
    for(int par=0;par<2;par++) vector_reset(source[par]);
    
    int islocal=1;
    coords lx;
    for(int mu=0;mu<NDIM;mu++)
      {
        lx[mu]=x[mu]-rank_coord[mu]*loc_size[mu];
        islocal&=(lx[mu]>=0);
        islocal&=(lx[mu]<loc_size[mu]);
      }
    
    if(islocal)
      {
	int ivol=loclx_of_coord(lx);
	for(int ic=0;ic<3;ic++) source[loclx_parity[ivol]][loceo_of_loclx[ivol]][ic][0]=1;
      }
    
    for(int par=0;par<2;par++) set_borders_invalid(source[par]);
  }
  THREADABLE_FUNCTION_END
}
