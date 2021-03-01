#include "nissa.hpp"

#include <random>

using namespace nissa;

int nhits;
int seed;
double wall_time;

double mass;
double kappa;
double residue;

int nanalyzed_conf;
int ngauge_conf;

quad_su3 *glb_conf;
spincolor *source,*temp;
spincolor *prop[2];

int iconf=0;
char conf_path[1024];
char outfolder[1024];
char run_file[1024];
lock_file_t<uint64_t> lock_file;
double init_time;

/*

The SitmoRng produces a stream of uint32, which must be rectified into
a stream of double. The stream is associated to lattice sites, such
that the i-th random number of the j-th site is the k-th double
number, with k=j+glb_vol*i. In this way, to draw n random numbers r_i
for each lattice site, one needs to move with glb_vol stride along the
stream. It this convenient to define a structure FieldRngStream which
wraps the SitmoRng, providing a view on the stream, in which the
number of random numbers drawn is kept homoegeneous across the
lattice. The user must draw N numbers in turns, requiring the kind of
local structure T to be be filled. The request produces a proxy,
FieldRngOf<T>, which keeps track of the request. When the request is
made, the status of the FieldRngStream is advanced. The user can pick
up the requested numbers from the FieldRngOf<T>, requiring a specific
site. Only at this moment the actual Sitmo is accessed, and the numbers generated.

SitmoRng
--------
-initialized on a given seed
-given an uint64, it produces an uint32 in the range 0-max

FieldRngStream
--------------
-embeds the SitmoRng
-keeps track of the number of real random real generated across the
lattice, which in turns correspond to a certain number of calls to the
SitmoRng
-accessed requiring to draw a certain number of double
-this call returns a proxy and advances the number of random reals

FieldRngOf<T>
--------------------
-keeps reference to the Sitmo of the FieldRngStream
-when passed a T, and the local lattice site, it fill T with the required elements
 */

namespace Sitmo
{
  /// Type of the key
  using Key=
    std::array<uint64_t,5>;
  
  /// Encrypted word
  using Word=
    std::array<uint64_t,4>;
  
  /// Encrypts the input
  inline Word encrypt(const Key& key,  ///< Key to encrypt
		      Word x)          ///< Input to encrypt
  {
    for(int j=0;j<5;j++)
      {
	for(int i=0;i<2;i++)
	  {
	    constexpr uint64_t mk[2][2]=
	      {{14,16},{25,33}};
	    
	    uint64_t& x0=x[2*i];
	    uint64_t& x1=x[(2*i+1)%4];
	    const uint64_t& rx=mk[j%2][i];
	    
	    x1+=key[(2*i+1+j)%5]+((i==1)?j:0);
	    x0+=x1+key[(2*i+j)%5];
	    x1=(x1<<rx)|(x1>>(64-rx));
	    x1^=x0;
	  }
	
	for(int l=0;l<3;l++)
	  for(int i=0;i<2;i++)
	    {
	      constexpr uint64_t m[2][3][2]=
		{{{52,57},{23,40},{5,37}},
		 {{46,12},{58,22},{32,32}}};
	      
	      uint64_t& x0=x[2*i];
	      uint64_t& x1=x[(2*i+((l%2==0)?3:1))%4];
	      const uint64_t& rx=m[j%2][l][i];
	      
	      x0+=x1;
	      x1=(x1<<rx)|(x1>>(64-rx));
	      x1^=x0;
	    }
      }
    
    // Increment entry 3
    x[3]+=5;
    
    return x;
  }
  
  /// Build a key from a word
  inline Key buildKey(const Word& word)
  {
    /// Output
    Key key;
    
    // Copy the first 4
    for(int i=0;i<4;i++)
      key[i]=word[i];
    
    // Set the fifth
    key[4]=0x1BD11BDAA9FC1A22^key[0]^key[1]^key[2]^key[3];
    
    return key;
  }
  
  /// Class encapsulating the Sitmo random generator
  struct Rng
  {
    /// Default seed
    static constexpr uint32_t DEFAULT_SEED()
    {
      return 3472291050;
    }
    
    /// Holds the state of the generator
    struct State : public Sitmo::Word
    {
      /// Increment of a certain amount
      State operator+(const uint64_t& z) const
      {
	/// Result
	State out=(*this);
	
	// Increment
	out[0]=(*this)[0]+z;
	
	// Overflow check
	if(out[0]<=(*this)[0])
	  {
	    /// Digit to check for overflow
	    int iDigit=1;
	    
	    // Carry over
	    do out[iDigit]++;
	    while(out[iDigit++]==0 and iDigit<4);
	  }
	
	return out;
      }
      
      /// Self increment
      State operator+=(const uint64_t& z)
      {
	return(*this)=(*this)+z;
      }
      
      /// Unitary self-increment
      State operator++()
      {
	return (*this)+=1;
      }
    };
    
    /// Seed, used to encrypt
    Key key;
    
    /// State of the random number generator
    State state{};
    
    /// Type returned
    using ResultType=uint32_t;
    
    static constexpr ResultType max_value=(~(ResultType)0);
    // /// Draw a seed from a random generator
    // template <typename F>
    // void seed(F& f)
    // {
    //   /// Result type of the callabale generator
    //   using FRes=typename F::ResultType;
      
    //   /// Number of calls to be issued to fill the state
    //   constexpr int nCall=sizeof(Word)/sizeof(FRes);
      
    //   union
    //   {
    // 	/// Partial key
    // 	Word partKey;
	
    // 	/// Access to the key in a way which allows to fill with the generator
    // 	FRes rawKey[nCall];
    //   };
      
    //   /// Fill the key
    //   for(int i=0;i<nCall;i++)
    // 	rawKey[i]=f();
      
    //   key=buildKey(partKey);
    // }
    
    /// Seed from a seed
    void seed(const uint32_t& s)
    {
      /// Partial word to be used
      const Word& partKey={s};
      
      key=buildKey(partKey);
    }
    
    /// Generate a number with a given offset w.r.t current state
    ResultType generateNth(const uint64_t& offset)
    {
      union
      {
	/// Temporary result of the encpter
	Word temp;
	
	/// Output to be returned
	const ResultType out{};
      };
      
      /// Shifted state
      State shiftedState=state+offset;
      
      temp=encrypt(key,static_cast<Word>(shiftedState));
      
      return out;
    }
    
    /// Default constructor
    explicit Rng(const uint32_t& s=DEFAULT_SEED())
    {
      seed(s);
    }
  };
}

//Embeds the Sitmo rng at a certain point in the stream
struct RngView
{
  /// Field random number generator
  Sitmo::Rng& ref;
  
  //Number of uint32_t generated so far
  uint64_t irng;
  
  //Minimal value
  constexpr uint32_t min()
  {
    return 0;
  }
  
  //Maximal value
  constexpr uint32_t max()
  {
    return ~0;
  }
  
  //Constructor
  RngView(Sitmo::Rng& ref,const uint64_t& irng) :
    ref(ref),irng(irng)
  {
  }
  
  /// Returns an uint32
  uint32_t operator()()
  {
    return ref.generateNth(irng++);
  }
};


template <typename T>
struct FieldRngOf
{
  //Mark if the generator can be used or not
  bool used;
  
  //Reference to the Sitmo
  Sitmo::Rng& rng;
  
  //Distribution [0,1)
  std::uniform_real_distribution<> distr;
  
  //Number of reals to be shifted, in units of global volume
  const uint64_t offsetReal;
  
  //Number of reals for each site
  static constexpr int nRealsPerSite=
    sizeof(T)/sizeof(double);
    
  //Constructor
  FieldRngOf(Sitmo::Rng& rng,const uint64_t& offsetReal) :
    used(false),rng(rng),offsetReal(offsetReal)
  {
    master_printf("All entries of FieldRng initialized\n");
  }
  
  /// Returns a view on a specific site and real number
  RngView getRngViewOnGlbSiteIRndReal(const int& glblx,const int& irnd_real_per_site)
  {
    //Computes the number in the stream of reals
    const uint64_t irnd_double=offsetReal+glblx+glb_vol*irnd_real_per_site;
    
    //Computes the offset in the rng stream
    const uint64_t irnd_uint32_t=2*irnd_double;
    
    return RngView(rng,irnd_uint32_t);
  }
  
  //Fill a specific site
  void _fillSite(T out,const int glblx)
  {
    //Flatten access to out
    double* reals=(double*)out;
    
    for(int irnd_real=0;irnd_real<nRealsPerSite;irnd_real++)
      {
	auto view=getRngViewOnGlbSiteIRndReal(glblx,irnd_real);
	
	reals[irnd_real]=distr(view);
      }
  }
  
  void enforce_single_usage()
  {
    if(used)
      crash("cannot use a rng field twice");
    used=true;
  }
  
  //Fill with origin
  void fillWithOrigin(T out)
  {
    enforce_single_usage();
    
    _fillSite(out,0);
  }
  
  //Fill all sites
  void fillField(T* out)
  {
    enforce_single_usage();
    
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(loclx,0,loc_vol)
      {
	//Finds the global site of local one
	const int& glblx=glblx_of_loclx[loclx];
	_fillSite(out[loclx],glblx);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
};

struct FieldRngStream
{
  //Embeds the random number generator
  Sitmo::Rng rng;
  
  //Number of double precision numbers generated per site
  uint64_t ndouble_gen;
  
  //Skip n drawers
  template <typename T>
  void skipDrawers(const uint64_t& nDrawers)
  {
    master_printf("Skipping drawer\n");
    ndouble_gen+=FieldRngOf<T>::nRealsPerSite*nDrawers;
  }
  
  //Return a proxy from which to fetch the random numbers
  template <typename T>
  auto getDrawer()
  {
    static_assert(sizeof(T)%sizeof(double)==0,"Type not multiple in size of double");
    
    master_printf("Creating the drawer\n");
    
    //Result
    FieldRngOf<T> res(rng,ndouble_gen);
    
    skipDrawers<T>(1);
    
    return res;
  }
  
  //Draw a single instance of T
  template <typename T>
  void drawScalar(T& out)
  {
    static_assert(sizeof(T)%sizeof(double)==0,"Type not multiple in size of double");
    
    auto drawer=getDrawer<T>();
    
    drawer.fillWithOrigin(out);
  }
  
  //Initilizes with a seed
  void init(const uint32_t& seed)
  {
    rng.seed(seed);
    ndouble_gen=0;
  }
};

FieldRngStream field_rng_stream;

static constexpr complex zero_complex={0.0,0.0};
static constexpr complex sqrt_2_half_complex{M_SQRT1_2,M_SQRT1_2};
void BoxMullerTransform(complex out,const complex ave=zero_complex,const complex sig=sqrt_2_half_complex)
{
  const double r=sqrt(-2*log(1-out[RE]));
  const double q=2*M_PI*out[IM];
  
  out[RE]=r*cos(q)*sig[RE]+ave[RE];
  out[IM]=r*sin(q)*sig[IM]+ave[IM];
}

void init_simulation(int narg,char **arg)
{
  lock_file.init();
  
  open_input(arg[1]);
  
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  init_grid(T,L);
  
  read_str_double("WallTime",&wall_time);
  
  read_str_int("Seed",&seed);
  
  read_str_int("NHits",&nhits);
  
  read_str_double("Kappa",&kappa);
  
  read_str_double("Mass",&mass);
  
  read_str_double("Residue",&residue);
  
  read_str_int("NGaugeConf",&ngauge_conf);
  
  glb_conf=nissa_malloc("glb_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  
  temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
  
  for(int r=0;r<2;r++)
    prop[r]=nissa_malloc("prop",loc_vol+bord_vol,spincolor);
}

//close the program
void close()
{
  close_input();
  
  master_printf("\n");
  master_printf("Nanalyzed confs: %d\n",nanalyzed_conf);
  master_printf("Total time: %lg s\n",take_time()-init_time);
  if(nanalyzed_conf)
    master_printf("Time per conf: %lg s\n",(take_time()-init_time)/nanalyzed_conf);
  master_printf("\n");
  
  nissa_free(glb_conf);
  for(int r=0;r<2;r++)
    nissa_free(prop[r]);
  nissa_free(source);
  nissa_free(temp);
}

//check if asked to stop or restart
bool asked_to_stop_or_restart()
{
  const bool asked_stop=file_exists("stop");
  master_printf("Asked to stop: %d\n",asked_stop);
  const int asked_restart=file_exists("restart");
  master_printf("Asked to restart: %d\n",asked_restart);
  
  return asked_stop or asked_restart;
}

//check that enough time is available
bool enough_time()
{
  const double passed_time=take_time()-init_time;
  const double remaining_time=wall_time-passed_time;
  const double time_per_conf=passed_time/(nanalyzed_conf+1e-10);
  const double time_per_conf_with_tol=time_per_conf*1.1;
  const bool no_analyzed_conf=(nanalyzed_conf==0);
  if(not no_analyzed_conf)
    {
      master_printf("Time per conf: %lg s (%d confs)\n",time_per_conf,nanalyzed_conf);
      master_printf("Time per conf with tol: %lg s\n",time_per_conf_with_tol);
      master_printf("Remaining time: %lg s\n",remaining_time);
    }
  
  const bool enough_remaining_time=no_analyzed_conf or (remaining_time>time_per_conf_with_tol);
  if(no_analyzed_conf)
    master_printf("No configuration analyzed yet, proceeding in any case\n");
  else
    if(enough_remaining_time)
      master_printf("Time is enough, going on\n");
    else
      master_printf("Not enough time\n");
  
  return enough_remaining_time;
}

//check if all confs has been analyzed
bool finished_confs()
{
  return iconf>=ngauge_conf;
}

bool read_conf_path_and_check_outpath_not_exists()
{
  read_str(conf_path,1024);
  
  read_str(outfolder,1024);
  
  master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
  
  const bool has_to_be_created=not dir_exists(outfolder);
  if(has_to_be_created)
      master_printf(" Output path \"%s\" not present.\n",outfolder);
  else
    master_printf(" Output path \"%s\" already present.\n",outfolder);
  
  return has_to_be_created;
}

bool create_outpath()
{
  const bool created=not create_dir(outfolder);
  
  if(created)
    master_printf(" Output path \"%s\" for conf \"%s\" created\n",outfolder,conf_path);
  else
    master_printf(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
  
  return created;
}

bool create_run_file()
{
  if(snprintf(run_file,1024,"%s/running",outfolder)<0) crash("writing %s",run_file);
  
  return lock_file.try_lock(run_file);
}

bool read_conf()
{
  read_ildg_gauge_conf(glb_conf,conf_path);
  master_printf("plaq: %+16.16g\n",global_plaquette_lx_conf(glb_conf));
  
  momentum_t old_theta={0,0,0,0};
  momentum_t theta={1,0,0,0};
  adapt_theta(glb_conf,old_theta,theta,0,0);
  
  return true;
}

bool check_lock_file()
{
  const bool lock_file_valid=lock_file.check_lock();
  
  if(not lock_file_valid)
    master_printf("Somebody acquired the lock on %s\n",run_file);
  
  return lock_file_valid;
}

bool check_if_next_conf_has_to_be_analyzed()
{
  return
    ((not asked_to_stop_or_restart()) and
     enough_time() and
     read_conf_path_and_check_outpath_not_exists() and
     create_outpath() and
     create_run_file() and
     read_conf() and
     check_lock_file());
}

void skip_conf()
{
  field_rng_stream.skipDrawers<spincolor>(nhits*glb_size[0]);
}

bool find_next_conf_not_analyzed()
{
  bool valid_conf=false;
  
  while((not finished_confs()) and not valid_conf)
    {
      valid_conf=check_if_next_conf_has_to_be_analyzed();
      
      if(not valid_conf)
	{
	  master_printf(" Configuration \"%s\" not to be processed\n",conf_path);
	  iconf++;
	  
	  skip_conf();
	}
    }
  
  if(valid_conf)
    master_printf(" Configuration \"%s\" valid, starting\n",conf_path);
  
  return valid_conf and not finished_confs();
}

//mark a conf as finished
void mark_finished()
{
  char fin_file[1024];
  if(snprintf(fin_file,1024,"%s/finished",outfolder)<0) crash("writing %s",fin_file);
  file_touch(fin_file);
  
  iconf++;
  nanalyzed_conf++;
}

void fill_source(const int glbT)
{
  GET_THREAD_ID();
  
  // double tFrT[1];
  // field_rng_stream.drawScalar(tFrT);
  // const int glbT=tFrT[0]*glb_size[0];
  master_printf("Source position: %d\n",glbT);
  
  auto source_filler=field_rng_stream.getDrawer<spincolor>();
  master_printf("Drawer initialized\n");
  source_filler.fillField(source);
  master_printf("Source filled\n");
  
  NISSA_PARALLEL_LOOP(loclx,0,loc_vol)
    {
      if(glbT==glb_coord_of_loclx[loclx][0])
	for(int id=0;id<NDIRAC;id++)
	  for(int ic=0;ic<NCOL;ic++)
	    BoxMullerTransform(source[loclx][id][ic]);
      else
	spincolor_put_to_zero(source[loclx]);
    }
  NISSA_PARALLEL_LOOP_END;
  
  master_printf("Box-Muller transformation performed\n");
  
  set_borders_invalid(source);
}

void get_prop(const int& r)
{
  safe_dirac_prod_spincolor(temp,(tau3[r]==-1)?&Pminus:&Pplus,source);
  inv_tmD_cg_eoprec(prop[r],NULL,glb_conf,kappa,mass*tau3[r],1000000,residue,temp);
  safe_dirac_prod_spincolor(prop[r],(tau3[r]==-1)?&Pminus:&Pplus,prop[r]);
}

void in_main(int narg,char **arg)
{
  const char stop_path[]="stop";
  
  init_time=take_time();
  
  //init simulation according to input file
  init_simulation(narg,arg);
  
  field_rng_stream.init(seed);
  
  //loop over the configs
  while(find_next_conf_not_analyzed())
    {
      FILE* disco_contr_file=open_file(combine("%s/disco_contr",outfolder),"w");
      FILE* conn_contr_file=open_file(combine("%s/conn_contr",outfolder),"w");
      
      for(int ihit=0;ihit<nhits;ihit++)
	{
	  complex disco_contr[2][glb_size[0]];
	  
	  for(int glbT=0;glbT<glb_size[0];glbT++)
	    {
	      fill_source(glbT);
	      
	      for(int r=0;r<2;r++)
		{
		  get_prop(r);
		  // const double n=double_vector_glb_norm2(source,loc_vol);
		  // master_printf("n: %lg\n",n);
		  
		  prop_multiply_with_gamma(temp,5,prop[r],-1);
		  
		  complex_vector_glb_scalar_prod(disco_contr[r][glbT],(complex*)source,(complex*)temp,loc_vol*sizeof(spincolor)/sizeof(complex));
		  
		  //master_printf("Prop: %lg %lg\n",prop[0][0][0][0],prop[0][0][0][1]);
		}
	      
	      for(int r1=0;r1<2;r1++)
		for(int r2=0;r2<2;r2++)
		  {
		    complex conn_contr[glb_size[0]];
		    for(int locT=0;locT<loc_size[0];locT++)
		      {
			const int cubeOrigin=loc_spat_vol*locT;
			const int glbTshifted=(glb_coord_of_loclx[cubeOrigin][0]-glbT+glb_size[0])%glb_size[0];
			complex* slice1=(complex*)(prop[r1][cubeOrigin]);
			complex* slice2=(complex*)(prop[r2][cubeOrigin]);
			
			complex_vector_glb_scalar_prod(conn_contr[glbTshifted],slice1,slice2,loc_spat_vol*sizeof(spincolor)/sizeof(complex));
		      }
		    
		    MPI_Allreduce(MPI_IN_PLACE,conn_contr,2*glb_size[0],MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
		    master_fprintf(conn_contr_file,"\n# hit %d , r1 %d , r2 %d\n\n",ihit,r1,r2);
		    for(int t=0;t<glb_size[0];t++)
		      master_fprintf(conn_contr_file,"%.16lg %.16lg\n",conn_contr[t][RE],conn_contr[t][IM]);
		  }
	    }
	  
	  for(int r=0;r<2;r++)
	    {
	      master_fprintf(disco_contr_file,"\n# hit %d , r %d\n\n",ihit,r);
	      for(int glbT=0;glbT<glb_size[0];glbT++)
		master_fprintf(disco_contr_file,"%.16lg %.16lg\n",disco_contr[r][glbT][RE],disco_contr[r][glbT][IM]);
	    }
	}
      
      mark_finished();
    }
  
  if(iconf>=ngauge_conf)
    {
      master_printf("Analyzed all confs, exiting\n\n");
      file_touch(stop_path);
    }
  
  //close the simulation
  close();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
