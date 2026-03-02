#ifndef _GAUGE_FIXING_HPP
#define _GAUGE_FIXING_HPP

#include <sstream>

#include <base/field.hpp>
#include <io/input.hpp>
#include <new_types/su3.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  struct LC_gauge_fixing_pars_t
  {
    /// At the same time this defines the method and the direction from which to start... maybe disentangle?
    enum gauge_t{LANDAU=0,COULOMB=1};
    
    /// Minimization method
    enum method_t{OVERRELAX,EXPONENTIATE};
    
    /// Converts the gauge to a string
    static inline std::string gaugeTag(const gauge_t& gauge)
    {
      std::string res;
      
      switch(gauge)
	{
	case LANDAU:
	  res="Landau";
	  break;
	case COULOMB:
	  res="Coulomb";
	  break;
	default:
	  CRASH("unknown gauge %d",gauge);
	}
      
      return res;
    }
    
    /// Converts a string to a enum
    static inline gauge_t gaugeFromTag(const std::string_view& tag)
    {
      if(tag=="Landau")
	return LANDAU;
      
      if(tag=="Coulomb")
	return COULOMB;
      
      CRASH("Unknown gauge %s, use \"Landau\" or \"Coulomb\"",tag.begin());
      
      return {};
    }
    
    /// Precision for minimization
    static constexpr double defTargetPrecision=1e-14;
    
    /// Default gauge
    static constexpr gauge_t defGauge=LANDAU;
    
    /// Max number of iterations
    static constexpr int defNmaxIterations=100000;
    
    /// Iterations between reunitarization
    static constexpr int defUnitarizeEach=100;
    
    /// Default method
    static constexpr method_t defMethod=EXPONENTIATE;
    
    // Default robability to overrelax
    static constexpr double defOverrelaxProb=0.9;
    
    /// Alpha in exp(-i alpha der /2)
    static constexpr double defAlphaExp=0.16;
    
    /// Use or not the adaptative search of 1405.5812
    static constexpr int defUseAdaptativeSearch=true;
    
    /// Use or not the generalized cg
    static constexpr int defUseGeneralizedCg=true;
    
    /// Use or not fft to accelerate
    static constexpr int defUseFftAcc=true;
    
    gauge_t gauge;
    
    double targetPrecision;
    
    int nmaxIterations;
    
    int unitarizeEach;
    
    method_t method;
    
    double overrelaxProb;
    
    double alphaExp;
    
    int useAdaptativeSearch;
    
    int useGeneralizedCg;
    
    int useFftAcc;
    
    static inline std::string methodTag(const method_t& method)
    {
      std::string res;
      
      switch(method)
	{
	case OVERRELAX:
	  res="Overrelax";
	  break;
	case EXPONENTIATE:
	  res="Exponentiate";
	  break;
	default:
	  CRASH("unknwon method %d",method);
	}
      
      return res;
    }
    
    static inline method_t methodFromTag(const std::string_view& tag)
    {
      if(tag=="Overrelax")
	return OVERRELAX;
      
      if(tag=="Exponentiate")
	return EXPONENTIATE;
      
      CRASH("Unknown method %s, use \"Overrelax\" or \"Exponentiate\"",tag.begin());
      
      return {};
    }
    
    int master_fprintf(FILE* fout,
		       const int full=false) const 
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      if(full or is_nonstandard())
	{
	  os<<"GaugeFix\n";
	  if(full or gauge!=defGauge) os<<" Gauge\t=\t"<<gaugeTag(gauge)<<"\n";
	  if(full or targetPrecision!=defTargetPrecision) os<<" TargetPrecision\t=\t"<<targetPrecision<<"\n";
	  if(full or nmaxIterations!=defNmaxIterations) os<<" TargetPrecision\t=\t"<<nmaxIterations<<"\n";
	  if(full or unitarizeEach!=defUnitarizeEach) os<<" TargetPrecision\t=\t"<<unitarizeEach<<"\n";
	  if(full or method!=defMethod) os<<" Method\t=\t"<<methodTag(method)<<"\n";
	  if(full or overrelaxProb!=defOverrelaxProb) os<<" OverrelaxProb\t=\t"<<overrelaxProb<<"\n";
	  if(full or alphaExp!=defAlphaExp) os<<" AlphaExp\t=\t"<<alphaExp<<"\n";
	  if(full or useAdaptativeSearch!=defUseAdaptativeSearch) os<<" UseAdaptativeSearch\t=\t"<<useAdaptativeSearch<<"\n";
	  if(full or useGeneralizedCg!=defUseGeneralizedCg) os<<" UseGeneralizedCg\t=\t"<<useGeneralizedCg<<"\n";
	  if(full or useFftAcc!=defUseFftAcc) os<<" UseFFTacc\t=\t"<<useFftAcc<<"\n";
	}
      return os.str();
    }
    
    bool is_nonstandard() const
    {
      return
	gauge!=defGauge or
	targetPrecision!=defTargetPrecision or
	nmaxIterations!=defNmaxIterations or
	unitarizeEach!=defUnitarizeEach or
	method!=defMethod or
	overrelaxProb!=defOverrelaxProb or
	alphaExp!=defAlphaExp or
	useAdaptativeSearch!=defUseAdaptativeSearch or
	useGeneralizedCg!=defUseGeneralizedCg or
	useFftAcc!=defUseFftAcc;
    }
    
    LC_gauge_fixing_pars_t() :
      gauge(defGauge),
      targetPrecision(defTargetPrecision),
      nmaxIterations(defNmaxIterations),
      unitarizeEach(defUnitarizeEach),
      method(defMethod),
      overrelaxProb(defOverrelaxProb),
      alphaExp(defAlphaExp),
      useAdaptativeSearch(defUseAdaptativeSearch),
      useGeneralizedCg(defUseGeneralizedCg),
      useFftAcc(defUseFftAcc)
    {
    }
  };
  
  //read all Landau or Coulomb gauge fixing pars
  inline void read_LC_gauge_fixing_pars(LC_gauge_fixing_pars_t &pars)
  {
    char gauge_tag[200];
    read_str_str("Gauge",gauge_tag,200);
    pars.gauge=LC_gauge_fixing_pars_t::gaugeFromTag(gauge_tag);
    read_str_double("TargetPrecision",&pars.targetPrecision);
    read_str_int("NMaxIterations",&pars.nmaxIterations);
    read_str_int("UnitarizeEach",&pars.unitarizeEach);
    char method_tag[200];
    read_str_str("Method",method_tag,200);
    pars.method=LC_gauge_fixing_pars_t::methodFromTag(method_tag);
    switch(pars.method)
      {
      case LC_gauge_fixing_pars_t::OVERRELAX:
	read_str_double("OverrelaxProb",&pars.overrelaxProb);
	break;
      case LC_gauge_fixing_pars_t::EXPONENTIATE:
	read_str_double("AlphaExp",&pars.alphaExp);
	read_str_int("UseAdaptativeSearch",&pars.useAdaptativeSearch);
	read_str_int("UseGeneralizedCg",&pars.useGeneralizedCg);
	read_str_int("UseFFTacc",&pars.useFftAcc);
	break;
      default:
	CRASH("unkown method %d",pars.method);
      }
  }
  
  /// Apply a gauge transformation to the conf
  void gauge_transform_conf(LxField<quad_su3>& uout,
			    const LxField<su3>& g,
			    const LxField<quad_su3>& uin);
  
  void gauge_transform_conf(eo_ptr<quad_su3> uout,eo_ptr<su3> g,eo_ptr<quad_su3> uin);
  
  void gauge_transform_color(eo_ptr<color> out,eo_ptr<su3> g,eo_ptr<color> in);
  
  void Landau_or_Coulomb_gauge_fix(LxField<quad_su3>& fixed_conf,
				   const LC_gauge_fixing_pars_t& pars,
				   const LxField<quad_su3>& ext_conf);
  
  /// Perform a random gauge transformation
  void perform_random_gauge_transform(LxField<quad_su3>& conf_out,
				      const LxField<quad_su3>& conf_in);
  
  void perform_random_gauge_transform(eo_ptr<quad_su3> conf_out,eo_ptr<quad_su3> conf_in);
}

#endif
