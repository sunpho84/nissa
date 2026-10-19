#include "corrParts.hpp"

#include <atomic>

namespace nissa::nice
{
  namespace
  {
    std::atomic<std::uint64_t> nextSourceId{};
  }
  
  std::ostream& operator<<(std::ostream& out,
                           const Momentum& momentum)
  {
    return out<<"{"<<momentum.x<<","<<momentum.y<<","<<momentum.z<<"}";
  }
  
  Momentum operator-(const Momentum& momentum)
  {
    // Avoid turning a zero component into a negative zero
    const auto negate=
      [](const double value)
    {
      return value==0.0?value:-value;
    };
    
    return {negate(momentum.x),negate(momentum.y),negate(momentum.z)};
  }
  
  Momentum operator*(const double weight,
                     const Momentum& momentum)
  {
    return {weight*momentum.x,weight*momentum.y,weight*momentum.z};
  }
  
  Momentum operator*(const Momentum& momentum,
                     const double weight)
  {
    return weight*momentum;
  }
  
  Pars dag(const Pars& pars)
  {
    const auto conjugate=
      [](const auto& value)->Pars
    {
      return value.dag();
    };
    
    return std::visit(conjugate,pars);
  }
  
  std::string describe(const Pars& pars)
  {
    const auto describeValue=
      [](const auto& value)
    {
      return value.describe();
    };
    
    return std::visit(describeValue,pars);
  }
  
  Oper Oper::dag() const
  {
    Oper result=
      rhs?rhs->dag()*Oper(nissa::nice::dag(pars)):
      Oper(nissa::nice::dag(pars));
    result.weight*=
      weight;
    
    return result;
  }
  
  Oper operator*(const Oper& lhs,
                 const Oper& rhs)
  {
    Oper result(lhs);
    Oper* tail=
      &result;
    while(tail->rhs)
      tail=
        tail->rhs.get();
    tail->rhs=
      std::make_unique<Oper>(rhs);
    result.weight*=
      rhs.weight;
    tail->rhs->weight=1.0;
    
    return result;
  }
  
  Oper operator*(const double weight,
                 const Oper& oper)
  {
    Oper result(oper);
    result.weight*=
      weight;
    
    return result;
  }
  
  Oper operator*(const Oper& oper,
                 const double weight)
  {
    return weight*oper;
  }
  
  Source::Source(const int tSelect) :
    count(nextSourceId.fetch_add(1,std::memory_order_relaxed)),
    tSelect(tSelect)
  {
  }
  
  Line operator*(const Oper& oper,
                 const Line& line)
  {
    return Line(oper,line);
  }
  
  Line operator*(const double weight,
                 const Line& line)
  {
    if(const auto* source=
       std::get_if<Source>(&line.data))
      return Line(Pars(Gamma{.iGamma=0}),{{weight,Line(*source)}});
    
    Extender result=
      std::get<Extender>(line.data);
    for(auto& [factor,input] : result.rhs)
      {
        (void)input;
        factor*=
          weight;
      }
    
    return Line(std::move(result));
  }
  
  Line operator*(const Line& line,
                 const double weight)
  {
    return weight*line;
  }
  
  Line addLines(const Line& lhs,
                const Line& rhs,
                const bool merge)
  {
    const auto* left=
      std::get_if<Extender>(&lhs.data);
    const auto* right=
      std::get_if<Extender>(&rhs.data);
    
    if(merge and left and right and left->op==right->op)
      {
        Extender result=
          *left;
        result.rhs.insert(result.rhs.end(),right->rhs.begin(),right->rhs.end());
        
        return Line(std::move(result));
      }
    
    return Line(Extender{ExtendSelectingTPars{Pars(Gamma{.iGamma=0})},
                         {{1.0,lhs},{1.0,rhs}}});
  }
  
  Line operator+(const Line& lhs,
                 const Line& rhs)
  {
    return addLines(lhs,rhs);
  }
}
