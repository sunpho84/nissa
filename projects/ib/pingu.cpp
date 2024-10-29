#include <nissa.hpp>

using namespace nissa;

namespace nissa
{
  /// One of the fourth roots of unity
  enum class IdFourthRoot{One,I,MinusOne,MinusI};
  
  /// Product of two fourth roots of unity
  constexpr IdFourthRoot operator*(const IdFourthRoot& first,
				   const IdFourthRoot& second)
  {
    return IdFourthRoot(((int)first+(int)second)%4);
  }
  
  /// Complex conjugate of one of the fourth roots of unity
  constexpr IdFourthRoot operator~(const IdFourthRoot& i)
  {
    return IdFourthRoot(((int)i+((int)i%2)*2)%4);
  }
  
  /// Negation of one of the fourth roots of unity
  constexpr IdFourthRoot operator-(const IdFourthRoot& i)
  {
    return IdFourthRoot(((int)i+2)%4);
  }
  
  /// Dirac Gamma matrix, can contain only +/-1 or +/-i, and only one
  /// nonzero column per row
  struct DiracGamma
  {
    /// Type defining a single entry
    using Entry=
      std::pair<SpinCln,IdFourthRoot>;
    
    /// Type defining all the entries
    using Entries=
      StackTens<CompsList<SpinRow>,Entry>;
    
    /// Entries of the gamma matrix
    const Entries data;
    
    /// Product with another gamma
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr DiracGamma operator*(const DiracGamma& oth) const
    {
      ///Result gamma data
      Entries outData;
      
      for(SpinRow ig1=0;ig1<nDirac;ig1++)
	{
	  //This is the line to be taken on the second matrix
	  const SpinRow ig2=transp(data(ig1).first);
	  
	  //The entries of the output is, on each line, the complex
	  //product of the entries of the first matrix on that line, for
	  //the entries of the second matrix on the line with the index
	  //equal to the column of the first matrix which is different
	  //from 0 (which again is ig2)
	  const IdFourthRoot r=
	    data(ig1).second*oth.data(ig2).second;
	  
	  //For each line, the column of the output matrix which is
	  //different from 0 is the column of the second matrix different
	  //from 0 on the line with index equal to the column of the first
	  //matrix which is different from 0 (that is, ig2)
	  outData(ig1)={oth.data(ig2).first,r};
	};
      
      return {outData};
    }
    
    /// Negated Gamma
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr DiracGamma operator-() const
    {
      /// Result gamma data
      Entries outData;
      
      for(SpinRow i=0;i<nDirac;i++)
	outData(i)={outData(i).first,-outData(i).second};
      
      return {outData};
    }
    
    /// Hermitian matrix, transposed and conjugated
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr DiracGamma dag() const
    {
      /// Result gamma data
      Entries outData;
      
      for(SpinRow sr=0;sr<nDirac;sr++)
	{
	  const auto& [sc,c]=this->data(sr);
	  outData(transp(sc))={sc,~c};
	}
      
      return {outData};
    }
  };
  
  //static constexpr Gamma Id{{{0,IdFourthRoot::One}}};
}

namespace GammaBasis
{
#define E(P,V) Gamma::Entry{P,IdFourthRoot::V}
#define GAMMA(NAME,P0,V0,P1,V1,P2,V2,P3,V3)				\
  constexpr Gamma NAME{{E(P0,V0),E(P1,V1),E(P2,V2),E(P3,V3)}}
  
  GAMMA(X,3,MinusI,2,MinusI,1,I,0,I);
  GAMMA(Y,3,MinusOne,2,One,1,One,0,MinusOne);
  GAMMA(Z,2,MinusI,3,I,0,I,1,MinusI);
  GAMMA(T,2,MinusOne,3,MinusOne,0,MinusOne,1,MinusOne);
  
#undef GAMMA
#undef E
  constexpr DiracGamma G5=T*X*Y*Z;
  [[maybe_unused]]
  constexpr DiracGamma e=-G5.dag();
  [[maybe_unused]]
  constexpr std::array<DiracGamma,16> basis{X*X,X,Y,Z,T,G5,G5*X,G5*Y,G5*Z,G5*T,T*X,T*Y,T*Z,Y*Z,Z*X,X*Y};
}

DECLARE_TRANSPOSABLE_COMP(SpinTcc,int,4,spinTcc);
DECLARE_TRANSPOSABLE_COMP(SpinTbs,int,4,spinTbs);
DECLARE_TRANSPOSABLE_COMP(ColTcc,int,3,colTcc);
DECLARE_TRANSPOSABLE_COMP(ColTbs,int,3,colTbs);
DECLARE_DYNAMIC_COMP(TimeTcc);
DECLARE_DYNAMIC_COMP(TimeTbs);

void in_main(int narg,char **arg)
{
  const int T=64;
  TimeTbs tBS(T);
  TimeTcc tCC(T);
  
  DynamicTens<OfComps<TimeTcc,SpinTccRow,ColTccRow,SpinTccCln,ColTccCln,ComplId>,double> CC(tCC);
  DynamicTens<OfComps<TimeTbs,SpinTbsRow,ColTbsRow,SpinTbsCln,ColTbsCln,ComplId>,double> BS(tBS);
  
  CC(spinTcc(0));
  
  auto CCc=traceOver<SpinTcc>(traceOver<ColTcc>(CC));
  auto BSc=traceOver<SpinTbs>(traceOver<ColTbs>(BS));
  
  DynamicTens<OfComps<TimeTcc,TimeTbs,ComplId>,double> c;
  c=CCc*BSc;
  
}

int main(int narg,char **arg)
{
  initNissa(narg,arg);
  in_main(narg,arg);
  closeNissa();
  
  return 0;
}
