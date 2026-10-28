#ifndef _CORR_PARTS_HPP
#define _CORR_PARTS_HPP

#include "base/debug.hpp"

#include <array>
#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace nissa::nice
{
  /// Three spatial momentum components; temporal boundary conditions are separate
  struct SpatMomentum
  {
    /// Spatial momentum component along x
    double x{};
    
    /// Spatial momentum component along y
    double y{};
    
    /// Spatial momentum component along z
    double z{};
    
    /// Access x, y or z by axis number; reject any other axis
    double& operator[](const std::size_t axis)
    {
      switch(axis)
        {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: CRASH("SpatMomentum axis must be 0, 1 or 2");
        }
    }
    
    /// Access a component of a const momentum
    const double& operator[](const std::size_t axis) const
    {
      switch(axis)
        {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: CRASH("SpatMomentum axis must be 0, 1 or 2");
        }
    }
    
    /// Return a value copy suitable for array-based consumers
    std::array<double,3> asArray() const
    {
      return {x,y,z};
    }
    
    /// Divide each momentum component by the same scalar
    SpatMomentum operator/(const double divisor) const
    {
      return {x/divisor,y/divisor,z/divisor};
    }
    
    /// Compare momentum components in order
    auto operator<=>(const SpatMomentum&) const=default;
  };
  
  /// Print a momentum in the form {x,y,z}
  std::ostream& operator<<(std::ostream& out,
                           const SpatMomentum& momentum);
  
  /// Negate a momentum without creating negative zero components
  SpatMomentum operator-(const SpatMomentum& momentum);
  
  /// Scale each component of a momentum
  SpatMomentum operator*(const double weight,
                         const SpatMomentum& momentum);
  
  /// Scale each component of a momentum with the scalar on the right
  SpatMomentum operator*(const SpatMomentum& momentum,
                         const double weight);
  
  /// Smearing parameters, including their momentum phase
  struct Smear
  {
    /// Smearing strength
    double kappa{};
    
    /// Smearing iteration count
    double n{};
    
    /// Spatial momentum inserted by the smearing
    SpatMomentum mom{};
    
    /// Smearing parameters are unchanged by conjugation
    Smear dag() const
    {
      return *this;
    }
    
    /// Describe the smearing parameters
    std::string describe() const
    {
      std::ostringstream out;
      out<<"Smear(kappa="<<kappa<<",n="<<n<<",mom="<<mom<<")";
      
      return out.str();
    }
    
    /// Compare smearing parameters in member order
    auto operator<=>(const Smear&) const=default;
  };
  
  /// Spatial momentum phase; conjugation reverses its momentum
  struct Phase
  {
    /// Spatial momentum inserted by the phase
    SpatMomentum mom{};
    
    /// Reverse the momentum under conjugation
    Phase dag() const
    {
      Phase result=
        *this;
      result.mom=
        -mom;
      
      return result;
    }
    
    /// Describe the momentum phase
    std::string describe() const
    {
      std::ostringstream out;
      out<<"Phase(mom="<<mom<<")";
      
      return out.str();
    }
    
    /// Compare phase momenta
    auto operator<=>(const Phase&) const=default;
  };
  
  /// Gamma matrix identified by its index; zero denotes the scalar
  struct Gamma
  {
    /// Gamma matrix index; zero is the identity/scalar operator
    std::size_t iGamma{};
    
    /// Accept only the 16 standard gamma matrices
    void validate() const
    {
      if(iGamma>=16)
        CRASH("invalid gamma matrix index %zu: expected 0 through 15",iGamma);
    }
    
    /// Gamma parameters are unchanged by conjugation
    Gamma dag() const
    {
      return *this;
    }
    
    /// Describe the gamma index
    std::string describe() const
    {
      return "Gamma(iGamma="+std::to_string(iGamma)+")";
    }
    
    /// Compare gamma indices
    auto operator<=>(const Gamma&) const=default;
  };
  
  /// Select a time separation in a line expression
  struct DeltaT
  {
    /// Selected nonnegative time separation
    std::size_t t{};
    
    /// Time selection is unchanged by conjugation
    DeltaT dag() const
    {
      return *this;
    }
    
    /// Describe the selected separation
    std::string describe() const
    {
      return "DeltaT(t="+std::to_string(t)+")";
    }
    
    /// Compare selected separations
    auto operator<=>(const DeltaT&) const=default;
  };
  
  /// Propagator parameters; residue belongs to the propagator alone
  struct Prop
  {
    /// Hopping parameter
    double kappa{};
    
    /// Bare mass parameter
    double mass{};
    
    /// Wilson sign parameter
    int r{};
    
    /// Electric charge
    double charge{};
    
    /// Propagator momentum
    SpatMomentum mom{};
    
    /// Residue used by the propagator
    double residue{};
    
    /// Conjugation reverses r and momentum; r must be +1 or -1
    Prop dag() const
    {
      if(r!=-1 and r!=1)
        CRASH("Prop::dag requires r to be +1 or -1");
      
      Prop result=
        *this;
      result.r=
        -r;
      result.mom=
        -mom;
      
      return result;
    }
    
    /// Describe the propagator parameters
    std::string describe() const
    {
      std::ostringstream out;
      out<<"Prop(kappa="<<kappa<<",mass="<<mass<<",r="<<r
         <<",charge="<<charge<<",mom="<<mom<<",residue="<<residue<<")";
      
      return out.str();
    }
    
    /// Compare propagator parameters in member order
    auto operator<=>(const Prop&) const=default;
  };
  
  /// Parameter alternatives accepted by an operator or line extender
  using Pars=
    std::variant<Smear,Phase,Gamma,DeltaT,Prop>;
  
  /// Conjugate whichever parameter alternative is held
  Pars dag(const Pars& pars);
  
  /// Describe whichever parameter alternative is held
  std::string describe(const Pars& pars);
  
  /// Ordered product of parameters, with a scalar weight on the product
  struct Oper
  {
    /// Leftmost parameter in the product
    Pars pars;
    
    /// Remaining product, if any
    std::unique_ptr<Oper> rhs;
    
    /// Scalar factor applying to the whole product
    double weight{1.0};
    
    /// Build a one-parameter operator product
    template <typename T>
    requires std::constructible_from<Pars,T>
    Oper(const T& value) :
      pars(value)
    {
    }
    
    /// Copy the complete product chain rather than sharing its tail
    Oper(const Oper& other) :
      pars(other.pars),
      rhs(other.rhs?std::make_unique<Oper>(*other.rhs):nullptr),
      weight(other.weight)
    {
    }
    
    /// Replace this product with a deep copy of another
    Oper& operator=(const Oper& other)
    {
      if(this!=&other)
        {
          pars=
            other.pars;
          rhs=
            other.rhs?std::make_unique<Oper>(*other.rhs):nullptr;
          weight=
            other.weight;
        }
      
      return *this;
    }
    
    /// Move an operator product
    Oper(Oper&&) noexcept=default;
    
    /// Move-assign an operator product
    Oper& operator=(Oper&&) noexcept=default;
    
    /// Conjugate the parameters and reverse their product order
    Oper dag() const;
    
    /// Describe the weighted product in parameter order
    std::string describe() const
    {
      std::ostringstream out;
      
      if(weight!=1.0)
        out<<weight<<"*";
      
      out<<nissa::nice::describe(pars);
      
      if(rhs)
        out<<"*"<<rhs->describe();
      
      return out.str();
    }
    
    /// Compare complete product structure and weight
    bool operator==(const Oper& other) const
    {
      return weight==other.weight and pars==other.pars and
        ((not rhs and not other.rhs) or
         (rhs and other.rhs and *rhs==*other.rhs));
    }
    
    /// Order products by parameters, weight and remaining chain
    std::partial_ordering operator<=>(const Oper& other) const
    {
      if(auto cmp=
         pars<=>other.pars;cmp!=0)
        return cmp;
      
      if(auto cmp=
         weight<=>other.weight;cmp!=0)
        return cmp;
      
      if(not rhs and not other.rhs)
        return std::partial_ordering::equivalent;
      
      if(not rhs)
        return std::partial_ordering::less;
      
      if(not other.rhs)
        return std::partial_ordering::greater;
      
      return *rhs<=>*other.rhs;
    }
  };
  
  /// Concatenate two operator products, retaining their order
  Oper operator*(const Oper& lhs,
                 const Oper& rhs);
  
  /// Scale an operator product on the left
  Oper operator*(const double weight,
                 const Oper& oper);
  
  /// Scale an operator product on the right
  Oper operator*(const Oper& oper,
                 const double weight);
  
  /// Unique source node, optionally fixed at a selected time
  struct Source
  {
    /// Process-wide identity assigned at construction
    std::uint64_t count;
    
    /// Selected time, or -1 when unrestricted
    int tSelect;
    
    /// Allocate a new source identity with an optional time selection
    explicit Source(const int tSelect=-1);
    
    /// Describe this source's identity and time selection
    std::string describe() const
    {
      return "Source(count="+std::to_string(count)+",tSelect="+
        std::to_string(tSelect)+")";
    }
    
    /// Compare source identities and time selections
    auto operator<=>(const Source&) const=default;
  };
  
  struct Line;
  
  /// Parameters applied to a line, with an optional time selection
  struct ExtendSelectingTPars
  {
    /// Parameter applied to each input line
    Pars pars;
    
    /// Selected time, or -1 when unrestricted
    int tSelect{-1};
    
    /// Describe the parameter and optional time selection
    std::string describe() const
    {
      std::string result=
        nissa::nice::describe(pars);
      
      if(tSelect!=-1)
        result+=
          "deltaT("+std::to_string(tSelect)+")";
      
      return result;
    }
    
    /// Compare parameters and time selections
    auto operator<=>(const ExtendSelectingTPars&) const=default;
  };
  
  /// Apply an operator to one line or combine weighted input lines
  struct Extender
  {
    /// Operation applied after its weighted inputs have been combined
    ExtendSelectingTPars op;
    
    /// Weighted input lines; a single input need not have unit weight
    std::vector<std::pair<double,Line>> rhs;
    
    /// Whether this extender applies the scalar gamma parameter
    bool isScalar() const
    {
      const auto* gamma=
        std::get_if<Gamma>(&op.pars);
      
      return gamma and gamma->iGamma==0;
    }
    
    /// Whether the inputs form a sum or have a nonunit weight
    bool isCombiner() const;
    
    /// Whether the first input is a source
    bool isExtendingSource() const;
    
    /// Whether the extender selects a time
    bool isSelectingTime() const
    {
      return op.tSelect!=-1;
    }
    
    /// Describe the operation and all weighted input lines
    std::string describe() const;
    
    /// Compare operations and weighted inputs
    bool operator==(const Extender& other) const;
    
    /// Order operations and weighted inputs
    std::partial_ordering operator<=>(const Extender& other) const;
  };
  
  /// Recursive correlation line rooted at a source or extender
  struct Line
  {
    /// A source leaf or an operation over earlier lines
    std::variant<Extender,Source> data;
    
    /// Optional label; not part of structural comparison
    std::string name;
    
    /// Build a line rooted at the given source
    Line(const Source& source) :
      data(source)
    {
    }
    
    /// Copy an extender and simplify the resulting line
    Line(const Extender& extender) :
      data(extender)
    {
      simplify();
    }
    
    /// Move an extender and simplify the resulting line
    Line(Extender&& extender) :
      data(std::move(extender))
    {
      simplify();
    }
    
    /// Build an extender and normalize time selection and redundant scalars
    Line(const Pars& pars,
         const std::vector<std::pair<double,Line>>& rhs,
         const int tSelect=-1) :
      data(Extender{ExtendSelectingTPars{pars,tSelect},rhs})
    {
      simplify();
    }
    
    /// Apply a complete operator product to a line
    Line(const Oper& oper,
         const Line& line) :
      data(Extender{ExtendSelectingTPars{oper.pars},{}})
    {
      auto& extender=
        std::get<Extender>(data);
      
      if(oper.rhs)
        extender.rhs.emplace_back(oper.weight,Line(*oper.rhs,line));
      else
        extender.rhs.emplace_back(oper.weight,line);
      
      simplify();
    }
    
    /// Copy a line and assign a new display name
    Line(const Line& other,
         const std::string& name) :
      data(other.data),
      name(name)
    {
    }
    
    /// Copy a line while retaining its source identity
    Line(const Line&)=default;
    
    /// Move a line
    Line(Line&&) noexcept=default;
    
    /// Copy-assign a line
    Line& operator=(const Line&)=default;
    
    /// Move-assign a line
    Line& operator=(Line&&) noexcept=default;
    
    /// Whether the root is an extender rather than a source
    bool isExtender() const
    {
      return std::holds_alternative<Extender>(data);
    }
    
    /// Describe this line recursively
    std::string describe() const
    {
      const auto describeValue=
        [](const auto& value)
      {
        return value.describe();
      };
      
      return std::visit(describeValue,data);
    }
    
    /// The name labels a line but does not affect its structural identity
    bool operator==(const Line& other) const
    {
      return data==other.data;
    }
    
    /// Order lines by structure, ignoring their display names
    std::partial_ordering operator<=>(const Line& other) const
    {
      return data<=>other.data;
    }
  
  private:
    
    /// Normalize time selection and remove redundant scalar operations
    void simplify()
    {
      removeDeltaT();
      removeUselessScalar();
    }
    
    /// Fold DeltaT parameters into the extender's time selector
    void removeDeltaT()
    {
      auto* extender=
        std::get_if<Extender>(&data);
      if(not extender)
        return;
      
      if(const auto* dt=
         std::get_if<DeltaT>(&extender->op.pars))
        {
          if(dt->t>static_cast<std::size_t>(std::numeric_limits<int>::max()))
            CRASH("DeltaT exceeds the supported time range");
          
          extender->op.tSelect=
            static_cast<int>(dt->t);
          extender->op.pars=
            Gamma{.iGamma=0};
        }
      
      for(auto& [weight,line] : extender->rhs)
        {
          (void)weight;
          line.removeDeltaT();
        }
    }
    
    /// Elide scalar extenders when they carry no additional operation
    void removeUselessScalar()
    {
      auto* extender=
        std::get_if<Extender>(&data);
      if(not extender)
        return;
      
      for(auto& [weight,line] : extender->rhs)
        {
          (void)weight;
          line.removeUselessScalar();
        }
      
      if(extender->isCombiner() or extender->rhs.empty() or not extender->isScalar())
        return;
      
      const Line& inner=
        extender->rhs.front().second;
      
      if(const auto* source=
         std::get_if<Source>(&inner.data))
        {
          if(extender->op.tSelect==-1 or source->tSelect==extender->op.tSelect)
            {
              Line replacement(inner);
              *this=
                std::move(replacement);
            }
          
          return;
        }
      
      if(const auto* nested=
         std::get_if<Extender>(&inner.data))
        if(nested->isScalar() and extender->op.tSelect==-1)
          {
            auto replacement=
              nested->rhs;
            extender->op.tSelect=
              nested->op.tSelect;
            extender->rhs=
              std::move(replacement);
          }
    }
  };
  
  inline bool Extender::isCombiner() const
  {
    return rhs.size()>1 or (rhs.size()==1 and rhs.front().first!=1.0);
  }
  
  inline bool Extender::isExtendingSource() const
  {
    return not rhs.empty() and not rhs.front().second.isExtender();
  }
  
  inline std::string Extender::describe() const
  {
    if(rhs.empty())
      CRASH("cannot describe an extender without an input line");
    
    std::ostringstream out;
    out<<op.describe();
    
    if(op.tSelect!=-1)
      out<<"*"<<DeltaT{.t=static_cast<std::size_t>(op.tSelect)}.describe();
    
    out<<"*";
    
    if(rhs.size()>1)
      out<<"(";
    
    for(std::size_t i=0;i<rhs.size();++i)
      {
        const auto& [weight,line]=
          rhs[i];
        
        if(weight<0)
          out<<weight<<"*";
        else
          {
            if(i)
              out<<"+";
            
            if(weight!=1.0)
              out<<weight<<"*";
          }
        
        out<<line.describe();
      }
    
    if(rhs.size()>1)
      out<<")";
    
    return out.str();
  }
  
  inline bool Extender::operator==(const Extender& other) const
  {
    return op==other.op and rhs==other.rhs;
  }
  
  inline std::partial_ordering Extender::operator<=>(const Extender& other) const
  {
    if(auto cmp=
       op<=>other.op;cmp!=0)
      return cmp;
    
    return rhs<=>other.rhs;
  }
  
  /// Apply an operator product to a line
  Line operator*(const Oper& oper,
                 const Line& line);
  
  /// Scale a line on the left
  Line operator*(const double weight,
                 const Line& line);
  
  /// Scale a line on the right
  Line operator*(const Line& line,
                 const double weight);
  
  /// Sum lines, merging equal outer operations when requested
  Line addLines(const Line& lhs,
                const Line& rhs,
                const bool merge=true);
  
  /// Sum lines with compatible outer operations merged
  Line operator+(const Line& lhs,
                 const Line& rhs);
}

#endif
