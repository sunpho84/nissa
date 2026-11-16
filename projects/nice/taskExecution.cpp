#include "taskExecution.hpp"

#include "base/debug.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
#include "inverters/twisted_mass/cg_invert_tmD_eoprec.hpp"
#include "inverters/twisted_clover/cg_invert_tmclovD_eoprec.hpp"
#include "io/reader.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/dirac.hpp"
#include "new_types/su3_op.hpp"
#include "operations/smearing/gaussian.hpp"
#include "routines/math_routines.hpp"
#include "threads/threads.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <span>
#include <vector>

namespace nissa::nice
{
  namespace
  {
    /// Return outer coordinates in the SField's physical component order
    std::vector<size_t> outerCoordinates(const LineValue& line,
                                         const size_t iCopy,
                                         const size_t iSpin,
                                         const size_t iColor)
    {
      std::vector<size_t> coordinates(line.field.outer.size());
      
      for(size_t i=0;i<line.field.outer.size();i++)
        {
          const Comp& component=line.field.outer[i];
	  
          if(component.kind==StandardComp::Copy)
            coordinates[i]=iCopy;
          else if(component.kind==StandardComp::Spin and
                  component.side==CompSide::Column)
            coordinates[i]=iSpin;
          else if(component.kind==StandardComp::Color and
                  component.side==CompSide::Column)
            coordinates[i]=iColor;
          else
            CRASH("NICE line has an unsupported outer component");
        }
      
      return coordinates;
    }
    
    /// Load every outer block using IB's hit, copy, and dilution naming convention
    void readILDGSource(LineValue& line,
                        const std::string& path,
                        const size_t iHit)
    {
      const std::filesystem::path stem(path);
      
      if(stem.filename().empty())
        CRASH("NICE readILDGLine requires a file stem after the directory");
      
      for(size_t iCopy=0;iCopy<line.nCopies;iCopy++)
        for(size_t iSpin=0;iSpin<(line.diluteSpin?NDIRAC:1);iSpin++)
          for(size_t iColor=0;iColor<(line.diluteColor?NCOL:1);iColor++)
            {
              std::string name=std::format("hit{}_{}",iHit,stem.filename().string());
              
              if(line.nCopies>1)
                name+=std::format("_copy{}",iCopy);
              
              name+=std::format("_idso{}_icso{}",iSpin,iColor);
              const std::string slotPath=(stem.parent_path()/name).string();
              const auto coordinates=outerCoordinates(line,iCopy,iSpin,iColor);
              auto field=line.field.template shadow<spincolor,
                defaultSpaceTimeLayout,defaultMemorySpace>(
                  "NICE ILDG source",std::span<const size_t>{coordinates});
              MASTER_PRINTF("NICE: reading line slot '%s'\n",slotPath.c_str());
              read_real_vector(field,slotPath,"scidac-binary-data");
              set_borders_invalid(field);
            }
    }
    
    /// Fill a source column with a point source at this hit's shared origin
    void fillPointSource(LineValue& line,
                         const Source& source,
                         const Coords& origin)
    {
      if(source.tSelect>=glbSize[0])
        CRASH("NICE source time %d is outside temporal extent %d",
              source.tSelect,glbSize[0]);
      
      const int sourceTime=source.tSelect<0?origin[0]:source.tSelect;
      const size_t nSourceSpins=line.diluteSpin?NDIRAC:1;
      const size_t nSourceColors=line.diluteColor?NCOL:1;
      
      for(size_t iCopy=0;iCopy<line.nCopies;iCopy++)
        for(size_t iSpin=0;iSpin<nSourceSpins;iSpin++)
          for(size_t iColor=0;iColor<nSourceColors;iColor++)
            {
              const std::vector<size_t> coordinates=
                outerCoordinates(line,iCopy,iSpin,iColor);
              auto field=line.field.template shadow<spincolor,
						    defaultSpaceTimeLayout,defaultMemorySpace>("NICE point source",std::span<const size_t>{coordinates});
	      
              field.template initOn<defaultMemorySpace>([origin,sourceTime,iSpin,iColor]
                                                        (LxField<spincolor>& output)
              {
                PAR(0,locVol,
                    CAPTURE(origin,sourceTime,iSpin,iColor,
                            TO_WRITE(output)),
                    site,
		    {
		      spincolor value{};
		      bool isOrigin=glbCoordOfLoclx[site][0]==sourceTime;
		      
		      for(int mu=1;mu<NDIM;mu++)
			isOrigin=isOrigin and glbCoordOfLoclx[site][mu]==origin[mu];
		      
		      if(isOrigin)
			value[iSpin][iColor][RE]=1.0;
		      
		      spincolor_copy(output[site],value);
		    });
		
                set_borders_invalid(output);
              });
            }
    }
    
    /// Add weighted line inputs into an allocated line field
    void sumLineInputs(LineValue& output,
                       const std::span<const WeightedLineInput> inputs)
    {
      if(inputs.empty())
        CRASH("NICE line extension has no inputs");
      
      for(size_t iCopy=0;iCopy<output.nCopies;iCopy++)
        for(size_t iSpin=0;iSpin<(output.diluteSpin?NDIRAC:1);iSpin++)
          for(size_t iColor=0;iColor<(output.diluteColor?NCOL:1);iColor++)
            {
              const std::vector<size_t> coordinates=
                outerCoordinates(output,iCopy,iSpin,iColor);
              auto destination=output.field.template shadow<spincolor,
							    defaultSpaceTimeLayout,defaultMemorySpace>(
												       "NICE weighted line",std::span<const size_t>{coordinates});
              destination.reset();
	      
              for(const WeightedLineInput& input : inputs)
                {
                  if(input.value==nullptr)
                    CRASH("NICE extension received a null line input");
		  
                  const std::vector<size_t> inputCoordinates=
                    outerCoordinates(*input.value,iCopy,iSpin,iColor);
                  auto source=input.value->field.template shadow<spincolor,
								 defaultSpaceTimeLayout,defaultMemorySpace>(
													    "NICE weighted input",std::span<const size_t>{inputCoordinates});
		  
                  if(input.weight!=0.0)
                    {
                      if(input.weight==1.0)
                        destination+=source;
                      else
                        {
                          LxField<spincolor> weighted("NICE weighted temporary");
                          weighted=source;
                          weighted*=input.weight;
                          destination+=weighted;
                        }
                    }
                }
	      
              set_borders_invalid(destination);
            }
    }
    
    /// Apply a gamma matrix to every sink spin index
    void applyGamma(LineValue& output,
                    const LineValue& input,
                    const Gamma& gamma)
    {
      gamma.validate();
      
      for(size_t iCopy=0;iCopy<input.nCopies;iCopy++)
        for(size_t iSpin=0;iSpin<(input.diluteSpin?NDIRAC:1);iSpin++)
          for(size_t iColor=0;iColor<(input.diluteColor?NCOL:1);iColor++)
            {
              const std::vector<size_t> coordinates=
                outerCoordinates(input,iCopy,iSpin,iColor);
              auto source=input.field.template shadow<spincolor,
						      defaultSpaceTimeLayout,defaultMemorySpace>(
												 "NICE gamma input",std::span<const size_t>{coordinates});
              auto destination=output.field.template shadow<spincolor,
							    defaultSpaceTimeLayout,defaultMemorySpace>(
												       "NICE gamma output",std::span<const size_t>{coordinates});
	      
              destination=source;
              PAR(0,locVol,
                  CAPTURE(iGamma=gamma.iGamma,TO_WRITE(destination),TO_READ(source)),
                  site,
		  {
		    safe_dirac_prod_spincolor(destination[site],base_gamma[iGamma],source[site]);
		  });
              set_borders_invalid(destination);
            }
    }
    
    /// Keep only sites at the selected relative time from the source origin
    void selectRelativeTime(LineValue& line,
                            const int deltaTime,
                            const Coords& origin)
    {
      if(deltaTime<0)
        return;
      
      if(deltaTime>=glbSize[0])
        CRASH("NICE time selection %d is outside temporal extent %d",
              deltaTime,glbSize[0]);
      
      const int selectedTime=(origin[0]+deltaTime)%glbSize[0];
      
      for(size_t iCopy=0;iCopy<line.nCopies;iCopy++)
        for(size_t iSpin=0;iSpin<(line.diluteSpin?NDIRAC:1);iSpin++)
          for(size_t iColor=0;iColor<(line.diluteColor?NCOL:1);iColor++)
            {
              const std::vector<size_t> coordinates=
                outerCoordinates(line,iCopy,iSpin,iColor);
              auto field=line.field.template shadow<spincolor,
						    defaultSpaceTimeLayout,defaultMemorySpace>(
											       "NICE time selection",std::span<const size_t>{coordinates});
	      
              PAR(0,locVol,
                  CAPTURE(selectedTime,TO_WRITE(field)),
                  site,
		  {
		    if(glbCoordOfLoclx[site][0]!=selectedTime)
		      spincolor_put_to_zero(field[site]);
		  });
              set_borders_invalid(field);
            }
    }
    
    /// Solve in the library's twisted basis while exposing physical-basis fields to NICE
    void solvePhysicalPropagator(LxField<spincolor>& solution,
                                 const LxField<spincolor>& source,
                                 const GaugeField& gauge,
                                 const Prop& prop)
    {
      const dirac_matr rotation=prop.r==1?Pplus:Pminus;
      LxField<spincolor> rotatedSource("NICE twisted-basis source",WITH_HALO);
      
      PAR(0,locVol,
          CAPTURE(rotation,TO_READ(source),TO_WRITE(rotatedSource)),
          site,
          {
            safe_dirac_prod_spincolor(rotatedSource[site],rotation,source[site]);
          });
      set_borders_invalid(rotatedSource);
      
      const AnisDopPars anis{false,1.0,false};
      // Use IB's solver tolerance when the NICE residue default is unset
      const double residue=prop.residue==0.0?1e-16:prop.residue;
      const double twistedMass=prop.r*prop.mass;
      
      if(prop.csw==0.0)
        inv_tmD_cg_eoprec(solution,std::nullopt,gauge,prop.kappa,
                          anis,twistedMass,1000000,residue,rotatedSource);
      else
        inv_tmclovD_cg_eoprec(solution,std::nullopt,gauge,prop.kappa,
                              anis,prop.csw,twistedMass,1000000,residue,rotatedSource);
      
      PAR(0,locVol,
          CAPTURE(rotation,TO_WRITE(solution)),
          site,
          {
            safe_dirac_prod_spincolor(solution[site],rotation,solution[site]);
          });
      set_borders_invalid(solution);
    }
    
    /// Apply a physical-basis propagator to every outer component of an input line
    void applyPropagator(LineValue& output,
                         const LineValue& input,
                         const GaugeRecipes& gaugeRecipes,
                         const Prop& prop)
    {
      prop.validate();
      
      if(prop.charge!=0.0 or prop.mom!=SpatMomentum{})
        CRASH("NICE propagator execution does not yet support charge or momentum");
      
      if(not std::isfinite(prop.residue) or prop.residue<0)
        CRASH("NICE propagator residue must be finite and nonnegative");
      
      const GaugeField& gauge=gaugeField(gaugeRecipes,prop.gaugeRecipe);
      
      for(size_t iCopy=0;iCopy<input.nCopies;iCopy++)
        for(size_t iSpin=0;iSpin<(input.diluteSpin?NDIRAC:1);iSpin++)
          for(size_t iColor=0;iColor<(input.diluteColor?NCOL:1);iColor++)
            {
              const std::vector<size_t> coordinates=
                outerCoordinates(input,iCopy,iSpin,iColor);
              auto source=input.field.template shadow<spincolor,
                defaultSpaceTimeLayout,defaultMemorySpace>(
                  "NICE inversion source",std::span<const size_t>{coordinates});
              auto solution=output.field.template shadow<spincolor,
                defaultSpaceTimeLayout,defaultMemorySpace>(
                  "NICE inversion result",std::span<const size_t>{coordinates});
	      
              solvePhysicalPropagator(solution,source,gauge,prop);
            }
    }
    
    /// Apply spatial Gaussian smearing to every outer component of an input line
    void applySmearing(LineValue& output,
                       const LineValue& input,
                       const GaugeRecipes& gaugeRecipes,
                       const Smear& smear)
    {
      smear.validate();
      
      if(smear.mom!=SpatMomentum{})
        CRASH("NICE smearing execution does not yet support momentum phases");
      
      const GaugeField& gauge=gaugeField(gaugeRecipes,smear.gaugeRecipe);
      const Momentum kappa{0.0,smear.kappa,smear.kappa,smear.kappa};
      
      for(size_t iCopy=0;iCopy<input.nCopies;iCopy++)
        for(size_t iSpin=0;iSpin<(input.diluteSpin?NDIRAC:1);iSpin++)
          for(size_t iColor=0;iColor<(input.diluteColor?NCOL:1);iColor++)
            {
              const std::vector<size_t> coordinates=
                outerCoordinates(input,iCopy,iSpin,iColor);
              auto source=input.field.template shadow<spincolor,
                defaultSpaceTimeLayout,defaultMemorySpace>(
                  "NICE smearing input",std::span<const size_t>{coordinates});
              auto destination=output.field.template shadow<spincolor,
                defaultSpaceTimeLayout,defaultMemorySpace>(
                  "NICE smearing result",std::span<const size_t>{coordinates});
	      
              gaussian_smearing(destination,source,gauge,kappa,
                                static_cast<int>(smear.n));
            }
    }
    
    /// Apply the relative momentum phase to every outer component of an input line
    void applyPhase(LineValue& output,
                    const LineValue& input,
                    const Phase& phase,
                    const Coords& origin)
    {
      for(size_t iCopy=0;iCopy<input.nCopies;iCopy++)
        for(size_t iSpin=0;iSpin<(input.diluteSpin?NDIRAC:1);iSpin++)
          for(size_t iColor=0;iColor<(input.diluteColor?NCOL:1);iColor++)
            {
              const std::vector<size_t> coordinates=
                outerCoordinates(input,iCopy,iSpin,iColor);
              auto source=input.field.template shadow<spincolor,
                defaultSpaceTimeLayout,defaultMemorySpace>(
                  "NICE phase input",std::span<const size_t>{coordinates});
              auto destination=output.field.template shadow<spincolor,
                defaultSpaceTimeLayout,defaultMemorySpace>(
                  "NICE phase result",std::span<const size_t>{coordinates});
	      
              PAR(0,locVol,
                  CAPTURE(momentumX=phase.mom.x,
                          momentumY=phase.mom.y,
                          momentumZ=phase.mom.z,
                          origin,
                          TO_WRITE(destination),TO_READ(source)),
                  site,
              {
                double argument=0.0;
		
                for(int mu=1;mu<NDIM;mu++)
                  {
                    const int relativeCoordinate=
                      (glbCoordOfLoclx[site][mu]-origin[mu]+glbSize[mu])%
                      glbSize[mu];
                    const double momentumComponent=mu==1?momentumX:
                      mu==2?momentumY:momentumZ;
                    argument+=M_PI*momentumComponent*relativeCoordinate/glbSize[mu];
                  }
		
                complex factor;
                complex_iexp(factor,argument);
                safe_spincolor_prod_complex(destination[site],source[site],factor);
              });
              set_borders_invalid(destination);
            }
    }
    
    /// Dispatch one operator over a line and its prepared gauge fields
    void applyExtensionOperator(LineValue& output,
                                const LineValue& input,
                                const Pars& pars,
                                const GaugeRecipes& gaugeRecipes,
                                const Coords& origin)
    {
      if(const auto* gamma=std::get_if<Gamma>(&pars))
        applyGamma(output,input,*gamma);
      else if(const auto* prop=std::get_if<Prop>(&pars))
        applyPropagator(output,input,gaugeRecipes,*prop);
      else if(const auto* smear=std::get_if<Smear>(&pars))
        applySmearing(output,input,gaugeRecipes,*smear);
      else if(const auto* phase=std::get_if<Phase>(&pars))
        applyPhase(output,input,*phase,origin);
      else
        CRASH("NICE line extension does not support this parameter type");
    }
    
    /// Dispatch one line extension over its weighted dependencies
    LineValue applyExtension(const RunOptions& options,
                             const GaugeRecipes& gaugeRecipes,
                             const ExtendSelectingTPars& extension,
                             const std::span<const WeightedLineInput> inputs,
                             const Coords& origin)
    {
      LineValue input(options);
      sumLineInputs(input,inputs);
      LineValue output(options);
      
      applyExtensionOperator(output,input,extension.pars,gaugeRecipes,origin);
      
      selectRelativeTime(output,extension.tSelect,origin);
      
      return output;
    }
    
    /// Calculate one IB-compatible meson gamma-pair contraction
    Tensor contractLines(const TraceContraction& contraction,
                         const LineValue& backward,
                         const LineValue& forward,
                         const Coords& origin)
    {
      if(not backward.sameShapeAs(forward))
        CRASH("NICE trace inputs have incompatible line shapes");
      
      if(not backward.diluteSpin)
        for(const auto& [sink,source] : contraction.gammaList)
          if(source!=5)
            CRASH("NICE non-spin-diluted contractions require source gamma 5");
      
      const size_t nGammaPairs=contraction.gammaList.size();
      Tensor result({ScalarKind::Complex,
                     {{StandardComp::GammaPair,CompSide::None,nGammaPairs},
                      {StandardComp::Time,CompSide::None,
                       static_cast<size_t>(glbSize[0])},
                      {StandardComp::Copy,CompSide::None,backward.nCopies}}},
                    MemorySpace::CPU);
      const size_t nSourceSpins=backward.diluteSpin?NDIRAC:1;
      const size_t nSourceColors=backward.diluteColor?NCOL:1;
      const double normalization=12.0/static_cast<double>(nSourceSpins*nSourceColors);
      size_t iGammaPair=0;
      
      for(const auto& [sinkGamma,sourceGamma] : contraction.gammaList)
        {
          const dirac_matr gammaSink=base_gamma[sinkGamma];
          const dirac_matr gammaSource=base_gamma[sourceGamma];
          const dirac_matr gammaFive=base_gamma[5];
	  
          for(size_t iCopy=0;iCopy<backward.nCopies;iCopy++)
            {
              LxField<complex> local("NICE local contraction",WITHOUT_HALO);
              local.reset();
	      
              for(size_t iSourceSpin=0;iSourceSpin<nSourceSpins;iSourceSpin++)
                for(size_t iSourceColor=0;iSourceColor<nSourceColors;iSourceColor++)
                  {
                    const size_t backwardSpin=backward.diluteSpin?
                      gammaSource.pos[iSourceSpin]:0;
                    const size_t forwardSpin=backward.diluteSpin?iSourceSpin:0;
                    const std::vector<size_t> backwardCoordinates=
                      outerCoordinates(backward,iCopy,backwardSpin,iSourceColor);
                    const std::vector<size_t> forwardCoordinates=
                      outerCoordinates(forward,iCopy,forwardSpin,iSourceColor);
                    auto backwardField=backward.field.template shadow<spincolor,
                      defaultSpaceTimeLayout,defaultMemorySpace>(
                        "NICE backward contraction input",
                        std::span<const size_t>{backwardCoordinates});
                    auto forwardField=forward.field.template shadow<spincolor,
                      defaultSpaceTimeLayout,defaultMemorySpace>(
                        "NICE forward contraction input",
                        std::span<const size_t>{forwardCoordinates});
		    
                    PAR(0,locVol,
                        CAPTURE(gammaSink,gammaSource,gammaFive,
                                iSourceSpin,normalization,
                                TO_WRITE(local),
                                TO_READ(backwardField),TO_READ(forwardField)),
                        site,
                    {
                      const int pairedSpin=gammaSource.pos[iSourceSpin];
                      complex sourceFactor;
                      unsafe_complex_prod(sourceFactor,
                                          gammaSource.entr[iSourceSpin],
                                          gammaFive.entr[pairedSpin]);
		      
                      for(int iSinkSpin=0;iSinkSpin<NDIRAC;iSinkSpin++)
                        {
                          const int forwardSinkSpin=gammaSink.pos[iSinkSpin];
                          complex sinkFactor;
                          complex factor;
                          unsafe_complex_prod(sinkFactor,
                                              gammaFive.entr[iSinkSpin],
                                              gammaSink.entr[iSinkSpin]);
                          unsafe_complex_prod(factor,sourceFactor,sinkFactor);
                          complex_prodassign_double(factor,normalization);
			  
                          complex colorSum{};
			  
                          UNROLL_FOR_ALL_COLS(iColor)
                            complex_summ_the_conj1_prod(
                              colorSum,
                              backwardField[site][iSinkSpin][iColor],
                              forwardField[site][forwardSinkSpin][iColor]);
			  
                          complex_summ_the_prod(local[site],colorSum,factor);
                        }
                    });
                  }
	      
              std::vector<complex> global(static_cast<size_t>(glbSize[0]));
              glb_reduce(global.data(),local,locVol,glbSize[0],locSize[0],
                         glbCoordOfLoclx[0][0]);
	      
              for(size_t iTime=0;iTime<global.size();iTime++)
                {
                  const size_t relativeTime=
                    (iTime+global.size()-static_cast<size_t>(origin[0]))%global.size();
                  const size_t index=
                    (iGammaPair*global.size()+relativeTime)*backward.nCopies+iCopy;
                  result.data[index*2]=global[iTime][RE];
                  result.data[index*2+1]=global[iTime][IM];
                }
            }
	  
          iGammaPair++;
        }
      
      return result;
    }
  }
  
  TaskHandlers makeNissaTaskHandlers(const RunOptions& options,
				     const GaugeRecipes& gaugeRecipes,
				     const Sitmo::Rng& rng,
				     const RngLayout& rngLayout,
				     const size_t iGaugeConf,
				     const size_t iHit)
  {
    if(iHit>=static_cast<size_t>(options.policy.nHits))
      CRASH("NICE RNG hit index %zu exceeds configured hit count %d",
            iHit,options.policy.nHits);
    
    const OperationRngRange originRange=
      rngLayout.operationRange(rng,rngLayout.origin,iGaugeConf,iHit);
    const Coords origin=
      options.physics.drawGlobalOriginTask.execute(originRange);
    
    MASTER_PRINTF("NICE: hit %zu source origin (%d,%d,%d,%d)\n",
                  iHit,origin[0],origin[1],origin[2],origin[3]);
    
    return TaskHandlers{
      .source=[&options,origin,iHit](const Source& source)
      {
        LineValue line(options);
        
        if(source.ildgPath)
          readILDGSource(line,*source.ildgPath,iHit);
        else
          fillPointSource(line,source,origin);
	
        return line;
      },
      .extend=[&options,&gaugeRecipes,origin](const ExtendSelectingTPars& extension,
                                             const std::span<const WeightedLineInput> inputs)
      {
        return applyExtension(options,gaugeRecipes,extension,inputs,origin);
      },
      .trace=[origin](const TraceContraction& contraction,
                      const LineValue& backward,
                      const LineValue& forward)
      {
        return contractLines(contraction,backward,forward,origin);
      }};
  }
}
