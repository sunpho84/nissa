#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#define EXTERN_RANDOM
# include <random/randomGenerate.hpp>

#include <base/vectors.hpp>
#include <new_types/su3_op.hpp>

namespace nissa
{
  //print all the entries of the random generator into a string
  void convert_rnd_gen_to_text(char *text,rnd_gen *gen,int size)
  {
    int *ptr=(int*)gen;
    memset(text,0,size);
    
    //init output text
    snprintf(text,size,"%d",ptr[0]);
    
    //append all the elements
    for(int i=1;i<RAN2_NTAB+3;i++)
      {
	char temp[20];
	snprintf(temp,20," %d",ptr[i]);
	strncat(text,temp,size);
      }
    if(strlen(text)==1024) crash("use larger text");
  }
  
  //read all the entries of the random generator from a string
  void convert_text_to_rnd_gen(rnd_gen *gen,const char *text)
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
  void start_glb_rnd_gen(const char *text)
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
    if(locVol==0l) crash("grid not initalized!"); //nasty
    
    //Generate the true seed
    if(glb_rnd_gen_inited==0) start_glb_rnd_gen(seed);
    int internal_seed=(int)rnd_get_unif(&glb_rnd_gen,0,RAND_MAX-glbVol());
    
    //allocate the grid of random generator, one for site
    loc_rnd_gen=nissa_malloc("Loc_rnd_gen",locVol(),rnd_gen);
    for(LocLxSite ivol=0;ivol<locVol;ivol++)
      {
	/// Reference to global site
	const GlbLxSite& glbLx=glblxOfLoclx(ivol);
	
	int loc_seed=internal_seed+glbLx();
	if(loc_seed<0) crash("something went wrong with local seed: %d + %d = %d",internal_seed,glbLx,loc_seed);
	start_rnd_gen(&(loc_rnd_gen[ivol.nastyConvert()]),loc_seed);
      }
    loc_rnd_gen_inited=1;
    master_printf("Grid of local random generators initialized with internal seed: %d\n",internal_seed);
  }
  
  //init from text
  void start_loc_rnd_gen(const char *text)
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
  
  //generate a random postion
  void generate_random_coord(GlbCoords& c)
  {
    FOR_ALL_DIRS(mu)
      c(mu)=(int)(rnd_get_unif(&glb_rnd_gen,0,glbSize(mu)()));
  }
  
  //return a numer between 0 and 1
  CUDA_HOST_AND_DEVICE int rnd_get_pm_one(rnd_gen *gen)
  {
    double r=rnd_get_unif(gen,0,1);
    if(r>0.5) return 1;
    else return -1;
  }
  
  //return a Z2 complex
  CUDA_HOST_AND_DEVICE void rnd_get_Z2(complex out,rnd_gen *gen)
  {
    out[0]=rnd_get_pm_one(gen);
    out[1]=0;
  }
  
  //return a Z4 complex
  CUDA_HOST_AND_DEVICE void rnd_get_Z4(complex out,rnd_gen *gen)
  {
    out[0]=rnd_get_pm_one(gen)/(double)RAD2;
    out[1]=rnd_get_pm_one(gen)/(double)RAD2;
  }
  
  //return a ZN complex
  CUDA_HOST_AND_DEVICE void rnd_get_ZN(complex out,rnd_gen *gen,int N)
  {complex_iexp(out,2*M_PI*(int)rnd_get_unif(gen,0,N)/N);}
  
  //return a gaussian double
  double rnd_get_gauss_double(rnd_gen *gen,double ave,double sig)
  {
    double q,r;
    
    r=sqrt(-2*log(1-rnd_get_unif(gen,0,1)));
    q=2*M_PI*rnd_get_unif(gen,0,1);
    
    return r*cos(q)*sig+ave;
  }
  
  //return a gaussian complex with sigma=sig/sqrt(2)
  CUDA_HOST_AND_DEVICE void rnd_get_gauss_complex(complex out,rnd_gen *gen,complex ave,double sig)
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
  CUDA_HOST_AND_DEVICE void comp_get_rnd(complex out,rnd_gen *gen,enum rnd_t rtype)
  {
    complex z={0,0};
    switch(rtype)
      {
      case RND_ALL_PLUS_ONE: complex_put_to_real(out,+1);                   break;
      case RND_ALL_MINUS_ONE:complex_put_to_real(out,-1);                   break;
      case RND_UNIF:         complex_put_to_real(out,rnd_get_unif(gen,0,1));break;
      case RND_Z2:           rnd_get_Z2(out,gen);                           break;
      case RND_Z3:           rnd_get_Z3(out,gen);                           break;
      case RND_Z4:           rnd_get_Z4(out,gen);                           break;
      case RND_GAUSS:        rnd_get_gauss_complex(out,gen,z,1);            break;
      }
  }
  
  //get the type of random from a string
  rnd_t convert_str_to_rnd_t(const char *str)
  {
    for(int i=0;i<nrnd_type;i++)
      if(strcasecmp(str,rnd_t_str[i])==0)
	return (rnd_t)i;
    
    master_fprintf(stderr,"Error, unknown random string %s,known ones:\n",str);
    
    for(int i=0;i<nrnd_type;i++)
      master_fprintf(stderr," %s\n",rnd_t_str[i]);
    crash("Choose one of them");
    
    return (rnd_t)RND_ALL_MINUS_ONE;
  }
  
  //fill a grid of vectors with numbers between 0 and 1
  void rnd_fill_unif_loc_vector(double* v,int dps,double min,double max)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int i=0;i<dps;i++)
	v[(ivol*dps+i).nastyConvert()]=rnd_get_unif(&(loc_rnd_gen[ivol.nastyConvert()]),min,max);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(v);
  }
  
  //return a grid of +-x numbers
  void rnd_fill_pm_one_loc_vector(double* v,int nps)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int i=0;i<nps;i++)
	v[(ivol*nps+i).nastyConvert()]=rnd_get_pm_one(&(loc_rnd_gen[ivol.nastyConvert()]));
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(v);
  }
  
  //generate a spindiluted vector according to the passed type
  void generate_colorspindiluted_source(su3spinspin* source,enum rnd_t rtype,const GlbCoord& tWall)
  {
    //reset
    vector_reset(source);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      if(glbCoordOfLoclx(ivol,tDir)==tWall or tWall<0)
	{
	  comp_get_rnd(source[ivol.nastyConvert()][0][0][0][0],&(loc_rnd_gen[ivol.nastyConvert()]),rtype);
	  for(int c=0;c<NCOL;c++)
	    for(int d=0;d<NDIRAC;d++)
	      if(c or d)
		memcpy(source[ivol.nastyConvert()][c][c][d][d],source[ivol.nastyConvert()][0][0][0][0],sizeof(complex));
	  }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(source);
  }
  
  //generate a spindiluted vector according to the passed type
  void generate_spindiluted_source(colorspinspin* source,enum rnd_t rtype,const GlbCoord& tWall)
  {
    //reset
    vector_reset(source);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      if(glbCoordOfLoclx(ivol,tDir)==tWall or tWall<0)
	for(int ic=0;ic<NCOL;ic++)
	  {
	    comp_get_rnd(source[ivol.nastyConvert()][ic][0][0],&(loc_rnd_gen[ivol.nastyConvert()]),rtype);
	    for(int d=1;d<NDIRAC;d++)
	      memcpy(source[ivol.nastyConvert()][ic][d][d],source[ivol.nastyConvert()][ic][0][0],sizeof(complex));
	  }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(source);
  }
  
  //generate an undiluted vector according to the passed type
  void generate_undiluted_source(spincolor* source,enum rnd_t rtype,const GlbCoord& tWall)
  {
    vector_reset(source);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      if(glbCoordOfLoclx(ivol,tDir)==tWall or tWall<0)
	for(int id=0;id<NDIRAC;id++)
	  for(int ic=0;ic<NCOL;ic++)
	    comp_get_rnd(source[ivol.nastyConvert()][id][ic],&(loc_rnd_gen[ivol.nastyConvert()]),rtype);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(source);
  }
  
  //generate a fully undiluted source
  void generate_fully_undiluted_lx_source(color* source,enum rnd_t rtype,const GlbCoord& tWall,const Dir& dir)
  {
    vector_reset(source);
    
    NISSA_PARALLEL_LOOP(ilx,0,locVol)
      if(tWall<0 or glbCoordOfLoclx(ilx,dir)==tWall)
	for(int ic=0;ic<NCOL;ic++)
	  comp_get_rnd(source[ilx.nastyConvert()][ic],&(loc_rnd_gen[ilx.nastyConvert()]),rtype);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(source);
  }
  //eo version
  void generate_fully_undiluted_eo_source(color* source,enum rnd_t rtype,const GlbCoord& tWall,const Parity& par,const Dir& dir)
  {
    vector_reset(source);
    
    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
      {
	const LocLxSite& ilx=loclx_of_loceo(par,ieo);
	if(tWall<0 or glbCoordOfLoclx(ilx,dir)==tWall)
	  for(int ic=0;ic<NCOL;ic++)
	    comp_get_rnd(source[ieo.nastyConvert()][ic],&(loc_rnd_gen[ilx.nastyConvert()]),rtype);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(source);
  }
  
  void generate_fully_undiluted_eo_source(eo_ptr<color> source,enum rnd_t rtype,const GlbCoord& tWall,const Dir& dir)
  {
    for(Parity par=0;par<2;par++)
      generate_fully_undiluted_eo_source(source[par],rtype,tWall,par,dir);
  }
  
  //same for spincolor
  void generate_fully_undiluted_eo_source(spincolor* source,enum rnd_t rtype,const GlbCoord& tWall,const Parity& par,const Dir& dir)
  {
    vector_reset(source);
    
    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
      {
	const LocLxSite& ilx=loclx_of_loceo(par,ieo);
	if(tWall<0 or glbCoordOfLoclx(ilx,dir)==tWall)
	  for(int id=0;id<NDIRAC;id++)
	    for(int ic=0;ic<NCOL;ic++)
	    comp_get_rnd(source[ieo.nastyConvert()][id][ic],&(loc_rnd_gen[ilx.nastyConvert()]),rtype);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(source);
  }
  
  void generate_fully_undiluted_eo_source(eo_ptr<spincolor> source,enum rnd_t rtype,const GlbCoord& tWall,const Dir& dir)
  {
    FOR_BOTH_PARITIES(par)
      generate_fully_undiluted_eo_source(source[par],rtype,tWall,par,dir);
  }
  
  //generate a delta source
  void generate_delta_source(su3spinspin* source,const GlbCoords& x)
  {
    //reset
    vector_reset(source);
    
    LocLxSite ivol;
    Rank where;
    
    get_loclx_and_rank_of_coord(ivol,where,x);
    
    if(rank==where)
      for(int id=0;id<NDIRAC;id++)
	for(int ic=0;ic<NCOL;ic++)
	  source[ivol.nastyConvert()][ic][ic][id][id][0]=1;
    
    set_borders_invalid(source);
  }
  
  //generate a delta source
  void generate_delta_eo_source(eo_ptr<su3> source,const GlbCoords& x)
  {
    //reset
    for(int par=0;par<2;par++)
      vector_reset(source[par]);
    
    LocLxSite ivol;
    Rank where;
    
    get_loclx_and_rank_of_coord(ivol,where,x);
    
    if(rank==where)
      su3_put_to_id(source[loclx_parity(ivol).nastyConvert()][loceo_of_loclx(ivol).nastyConvert()]);
    
    for(int par=0;par<2;par++)
      set_borders_invalid(source[par]);
  }
  
    //Taken from M.D'Elia
#if NCOL == 3
  CUDA_HOST_AND_DEVICE void herm_put_to_gauss(su3 H,rnd_gen *gen,double sigma)
  {
    const double one_by_sqrt3=0.577350269189626;
    const double two_by_sqrt3=1.15470053837925;
    
    double r[8];
    for(size_t ir=0;ir<4;ir++)
      {
	complex rc,ave={0,0};
	rnd_get_gauss_complex(rc,gen,ave,sigma);
	r[ir*2+0]=rc[0];
	r[ir*2+1]=rc[1];
      }
    
    //real part of diagonal elements
    H[0][0][0]= r[2]+one_by_sqrt3*r[7];
    H[1][1][0]=-r[2]+one_by_sqrt3*r[7];
    H[2][2][0]=     -two_by_sqrt3*r[7];
    
    //put immaginary part of diagonal elements to 0
    H[0][0][1]=H[1][1][1]=H[2][2][1]=0;
    
    //remaining
    H[0][1][0]=H[1][0][0]=r[0];
    H[0][1][1]=-(H[1][0][1]=r[1]);
    H[0][2][0]=H[2][0][0]=r[3];
    H[0][2][1]=-(H[2][0][1]=r[4]);
    H[1][2][0]=H[2][1][0]=r[5];
    H[1][2][1]=-(H[2][1][1]=r[6]);
  }
#endif
  
  /////////////////////////////// Generate an hermitian matrix ///////////////////////
  
  // A gauss vector has complex components z which are gaussian distributed
  // with <z~ z> = sigma
  void color_put_to_gauss(color H,rnd_gen *gen,double sigma)
  {
    complex ave={0,0};
    for(size_t ic=0;ic<NCOL;ic++) rnd_get_gauss_complex(H[ic],gen,ave,sigma);
  }
  
  //put a matrix to random used passed random generator
  CUDA_HOST_AND_DEVICE void su3_put_to_rnd(su3 u_ran,rnd_gen &rnd)
  {
    su3_put_to_id(u_ran);
    
    for(size_t i1=0;i1<NCOL;i1++)
      for(size_t i2=i1+1;i2<NCOL;i2++)
	{
	  //generate u0,u1,u2,u3 random on the four dim. sphere
	  double u0=rnd_get_unif(&rnd,-1,1);
	  double alpha=sqrt(1-u0*u0);
	  double phi=rnd_get_unif(&rnd,0,2*M_PI);
	  double costheta=rnd_get_unif(&rnd,-1,1);
	  double sintheta=sqrt(1-costheta*costheta);
	  double u3=alpha*costheta;
	  double u1=alpha*sintheta*cos(phi);
	  double u2=alpha*sintheta*sin(phi);
	  
	  //define u_l as unit matrix ...
	  su3 u_l;
	  su3_put_to_id(u_l);
	  
	  //... and then modify the elements in the chosen su(2) subgroup
	  u_l[i1][i1][RE]=u0;
	  u_l[i1][i1][IM]=u3;
	  u_l[i1][i2][RE]=u2;
	  u_l[i1][i2][IM]=u1;
	  u_l[i2][i1][RE]=-u2;
	  u_l[i2][i1][IM]=u1;
	  u_l[i2][i2][RE]=u0;
	  u_l[i2][i2][IM]=-u3;
	  
	  safe_su3_prod_su3(u_ran,u_l,u_ran);
	}
  }
  
  //return a single link after the heatbath procedure
  //routines to be shrunk!
  void su3_find_heatbath(su3 out,su3 in,su3 staple,double beta,int nhb_hits,rnd_gen *gen)
  {
    //compute the original contribution to the action due to the given link 
    su3 prod;
    unsafe_su3_prod_su3_dag(prod,in,staple);
    
    //copy in link to out
    if(out!=in) su3_copy(out,in);
    
    //iterate over heatbath hits
    for(int ihit=0;ihit<nhb_hits;ihit++)
      //scan all the three possible subgroups
      for(int isub_gr=0;isub_gr<3;isub_gr++)
	{
	  //take the part of the su3 matrix
	  double r0,r1,r2,r3;
	  double smod=su2_part_of_su3(r0,r1,r2,r3,prod,isub_gr);
	  
	  //omega is the coefficient of the plaquette, divided by the module of the su2 submatrix for normalization
	  double omega_f=beta/(3*smod);
	  double z_norm=exp(-2*omega_f);
	  omega_f=1/omega_f;
	  
	  double temp_f,z_f,a0;
	  do
	    {
	      double z_temp=(z_norm-1)*rnd_get_unif(gen,0,1)+1;
	      a0     = 1+omega_f*log(z_temp);
	      z_f    = 1-a0*a0;
	      temp_f = sqr(rnd_get_unif(gen,0,1))-z_f;
	    }
	  while(temp_f>0);
	  
	  double x_rat=sqrt(z_f);
	  
	  //generate an su2 matrix
	  double fi=rnd_get_unif(gen,0,2*M_PI);
	  double cteta=rnd_get_unif(gen,-1,1);
	  double steta=sqrt(1-cteta*cteta);
	  
	  double a1=steta*cos(fi)*x_rat;
	  double a2=steta*sin(fi)*x_rat;
	  double a3=cteta*x_rat;
	  
	  double x0 = a0*r0 + a1*r1 + a2*r2 + a3*r3;
	  double x1 = r0*a1 - a0*r1 + a2*r3 - r2*a3;
	  double x2 = r0*a2 - a0*r2 + a3*r1 - r3*a1;
	  double x3 = r0*a3 - a0*r3 + a1*r2 - r1*a2;
	  
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,prod);
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,out);
	}
  }
}
