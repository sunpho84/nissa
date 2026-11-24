#include <nissa.hpp>
#include <operations/smearing/gaussian.hpp>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <limits>
#include <string_view>
#include <utility>

namespace
{
  using namespace nissa;
  
  constexpr double smearingKappa=0.2;
  
  /// Parses a positive integer without accepting trailing characters
  int positiveArgument(const char* text)
  {
    const std::string_view input{text};
    int value{};
    const auto result=std::from_chars(input.data(),input.data()+input.size(),value);
    
    if(result.ec!=std::errc{} or result.ptr!=input.data()+input.size() or value<=0)
      CRASH("Expected a positive integer, received '%s'",text);
    
    return value;
  }
  
  /// Generates bounded values from global coordinates, independent of MPI partition
  CUDA_HOST_AND_DEVICE double sample(const Coords& coordinate,
                                     const int component)
  {
    int value=component*13;
    
    for(int mu=0;mu<NDIM;mu++)
      value+=(mu+1)*coordinate[mu];
    
    return double(value%101-50)/101;
  }
  
  /// Initializes nonzero complex spinors and nearly diagonal synthetic gauge links
  void initialize(LxField<spincolor>& input,
                  LxField<quad_su3>& gauge)
  {
    PAR(0,locVol,CAPTURE(TO_WRITE(input),TO_WRITE(gauge)),site,
	{
	  const Coords coordinate=glbCoordOfLoclx[site];
	  
	  for(int spin=0;spin<4;spin++)
	    for(int color=0;color<3;color++)
	      for(int ri=0;ri<2;ri++)
		input[site][spin][color][ri]=sample(coordinate,2*(3*spin+color)+ri);
	  
	  for(int mu=0;mu<NDIM;mu++)
	    for(int row=0;row<3;row++)
	      for(int col=0;col<3;col++)
		for(int ri=0;ri<2;ri++)
		  gauge[site][mu][row][col][ri]=
		    (row==col and ri==0 ? 1.0 : 0.0)+0.03*sample(coordinate,2*(9*mu+3*row+col)+ri);
	});
  }
  
  /// Evaluates table neighbors or arithmetic interior neighbors with halo fallback
  template <bool computedNeighbors,bool fused>
  void applyTransport(LxField<spincolor>& output,
                      const LxField<spincolor>& input,
                      const LxField<quad_su3>& gauge)
  {
    input.updateHalo();
    gauge.updateHalo();
    const Coords sizes=locSize;
    const double kappa=smearingKappa;
    
    // PAR provides the same tuner, field backup, and synchronization as production
    PAR(0,locVol,CAPTURE(TO_WRITE(output),TO_READ(input),TO_READ(gauge),sizes,kappa),site,
	{
	  spincolor sum;
	  spincolor_put_to_zero(sum);
	  int64_t stride=int64_t(sizes[2])*sizes[3];
	  
	  for(int mu=1;mu<NDIM;mu++)
	    {
	      int64_t up,down;
	      
	      if constexpr(computedNeighbors)
		{
		  const int coordinate=(site/stride)%sizes[mu];
		  up=coordinate+1<sizes[mu] ? site+stride : loclxNeighup[site][mu];
		  down=coordinate>0 ? site-stride : loclxNeighdw[site][mu];
		}
	      else
		{
		  up=loclxNeighup[site][mu];
		  down=loclxNeighdw[site][mu];
		}
	      
	      spincolor value;
	      unsafe_su3_prod_spincolor(value,gauge[site][mu],input[up]);
	      su3_dag_summ_the_prod_spincolor(value,gauge[down][mu],input[down]);
	      spincolor_summ_the_prod_double(sum,value,kappa);
	      
	      if(mu+1<NDIM)
		stride/=sizes[mu+1];
	    }
	  
	  if constexpr(fused)
	    {
	      const double normalization=1/(1+6*kappa);
	      
	      for(int spin=0;spin<4;spin++)
		for(int color=0;color<3;color++)
		  for(int ri=0;ri<2;ri++)
		    output[site][spin][color][ri]=
		      normalization*(input[site][spin][color][ri]+sum[spin][color][ri]);
	    }
	  else
	    spincolor_copy(output[site],sum);
	});
  }
  
  /// Runs the variants with real evolving fields and halo exchanges at every step
  template <bool computedNeighbors,bool fused>
  void smear(LxField<spincolor>& output,
             const LxField<spincolor>& input,
             const LxField<quad_su3>& gauge,
             LxField<spincolor>& temporary,
             LxField<spincolor>& hopping,
             const int iterations)
  {
    temporary=input;
    
    if constexpr(fused)
      {
        auto* current=&temporary;
        auto* next=&hopping;
	
        for(int iteration=0;iteration<iterations;iteration++)
          {
            applyTransport<computedNeighbors,true>(*next,*current,gauge);
            std::swap(current,next);
          }
	
        output=*current;
      }
    else
      {
        const double normalization=1/(1+6*smearingKappa);
	
        for(int iteration=0;iteration<iterations;iteration++)
          {
            hopping.reset();
            applyTransport<computedNeighbors,false>(hopping,temporary,gauge);
            temporary+=hopping;
            temporary*=normalization;
          }
	
        output=temporary;
      }
  }
  
  /// Returns the maximum elapsed wall time across MPI ranks, including GPU work
  template <typename F>
  double measure(const F& operation)
  {
#ifdef USE_CUDA
    DECRYPT_CUDA_ERROR(cudaDeviceSynchronize(),"Before smearing benchmark");
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    const double start=MPI_Wtime();
    operation();
#ifdef USE_CUDA
    DECRYPT_CUDA_ERROR(cudaDeviceSynchronize(),"After smearing benchmark");
#endif
    const double local=MPI_Wtime()-start;
    double maximum;
    MPI_Allreduce(&local,&maximum,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    
    return maximum;
  }
  
  /// Compares every physical field component with the current production routine
  template <typename Reference>
  void validate(const LxField<spincolor>& output,
                const Reference& reference)
  {
    const auto host=output.template getSurelyReadableOn<MemorySpace::CPU>();
    double error=0;
    
    for(int64_t site=0;site<locVol;site++)
      for(int component=0;component<24;component++)
        {
          const double expected=reference(site,component);
          const double actual=host(site,component);
	  
          if(not std::isfinite(actual))
            error=std::numeric_limits<double>::infinity();
          else
            error=std::max(error,std::abs(actual-expected)/(1+std::abs(expected)));
        }
    
    MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    
    if(error>1e-12)
      CRASH("Smearing comparison failed: maximum scaled error %.17g",error);
    
    MASTER_PRINTF("SMEAR validation maximum_scaled_error=%.17g\n",error);
  }
  
  /// Warms up tuning before timing repeated complete smearing calls
  template <typename F,typename Reference>
  void benchmark(const char* name,
                 const F& operation,
                 const LxField<spincolor>& output,
                 const Reference& reference,
                 const int iterations,
                 const int batches)
  {
    MASTER_PRINTF("SMEAR warmup case=%s (kernel tuning excluded from timings)\n",name);
    operation();
    validate(output,reference);
    double total=0;
    
    for(int batch=0;batch<batches;batch++)
      {
        const double seconds=measure(operation);
        total+=seconds;
        MASTER_PRINTF("SMEAR sample case=%s batch=%d ms_per_iteration=%.6f\n",name,batch,seconds*1000/iterations);
      }
    
    const double milliseconds=1000*total/(double(iterations)*batches);
    MASTER_PRINTF("SMEAR result case=%s mean_ms_per_iteration=%.6f estimated_global_GFLOP_per_s=%.3f\n",
                  name,milliseconds,1848.0*glbVol/(milliseconds*1e6));
  }
}

/// Benchmarks production smearing and localized variants using the configured backend
int main(int narg,
         char** arg)
{
  using namespace nissa;
  initNissa(narg,arg);
  
  if(narg!=1 and narg!=5 and narg!=9)
    CRASH("Usage: gaussian_smearing [L T iterations batches [rT rX rY rZ]]");
  
  const int length=narg>1 ? positiveArgument(arg[1]) : 16;
  const int time=narg>1 ? positiveArgument(arg[2]) : 32;
  const int iterations=narg>1 ? positiveArgument(arg[3]) : 20;
  const int batches=narg>1 ? positiveArgument(arg[4]) : 3;
  
  if(narg==9)
    for(int mu=0;mu<NDIM;mu++)
      fix_nranks[mu]=positiveArgument(arg[5+mu]);
  
  initGrid(time,length);
  MASTER_PRINTF("SMEAR L=%d T=%d ranks=%d rank_grid=%d,%d,%d,%d local=%d,%d,%d,%d iterations=%d batches=%d layout=%s\n",
                length,time,nranks,nRanksDir[0],nRanksDir[1],nRanksDir[2],nRanksDir[3],
                locSize[0],locSize[1],locSize[2],locSize[3],iterations,batches,spaceTimeLayoutName(defaultSpaceTimeLayout));
  {
    LxField<quad_su3> gauge("gauge",WITH_HALO);
    LxField<spincolor> input("input",WITH_HALO);
    LxField<spincolor> output("output",WITH_HALO);
    LxField<spincolor> temporary("temporary",WITH_HALO);
    LxField<spincolor> hopping("hopping",WITH_HALO);
    initialize(input,gauge);
    const Momentum kappa{0,smearingKappa,smearingKappa,smearingKappa};
    const auto baseline=[&]()
    {
      gaussian_smearing(output,input,gauge,kappa,iterations,&temporary,&hopping);
    };
    
    // Preserve the production result on the CPU, outside all measurement regions
    baseline();
    LxField<spincolor,defaultSpaceTimeLayout,MemorySpace::CPU> reference("reference",WITH_HALO);
    reference=output;
    benchmark("production",baseline,output,reference,iterations,batches);
    benchmark("four_pass_table",[&]()
    {
      smear<false,false>(output,input,gauge,temporary,hopping,iterations);
    },output,reference,iterations,batches);
    benchmark("four_pass_computed",[&]()
    {
      smear<true,false>(output,input,gauge,temporary,hopping,iterations);
    },output,reference,iterations,batches);
    benchmark("fused_table",[&]()
    {
      smear<false,true>(output,input,gauge,temporary,hopping,iterations);
    },output,reference,iterations,batches);
    benchmark("fused_computed",[&]()
    {
      smear<true,true>(output,input,gauge,temporary,hopping,iterations);
    },output,reference,iterations,batches);
    
    // Force fresh spinor exchanges to expose the unchanged input-halo cost
    temporary=input;
    temporary.updateHalo();
    const double haloSeconds=measure([&]()
    {
      for(int iteration=0;iteration<iterations;iteration++)
        {
          temporary.invalidateHalo();
          temporary.updateHalo();
        }
    });
    MASTER_PRINTF("SMEAR halo_only ms_per_exchange=%.6f (all library halo directions)\n",1000*haloSeconds/iterations);
  }
  closeNissa();
  
  return 0;
}
