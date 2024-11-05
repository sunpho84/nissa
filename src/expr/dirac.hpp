#ifndef _DIRAC_HPP
#define _DIRAC_HPP

#include <array>

#include <expr/comp.hpp>
#include <expr/conj.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/subExprs.hpp>

namespace nissa
{
#define NSPIN 4
  
  DECLARE_TRANSPOSABLE_COMP(Spin,int,NSPIN,spin);
  
  constexpr Spin nDirac=NSPIN;
  
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
  
  /// Unary plus
  constexpr IdFourthRoot operator+(const IdFourthRoot& i)
  {
    return i;
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Dirac Gamma matrix, can contain only +/-1 or +/-i, and only one
  /// nonzero column per row
  struct DiracGammaData
  {
    /// Non null column
    std::array<int,4> nnc;
    
    /// Value of the non-null entry
    std::array<IdFourthRoot,4> val;
    
    /// Negative of the present Gamma
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr DiracGammaData operator-() const
    {
      /// Result non null column
      std::array<int,4> outNnc;
      
      /// Value of non null entry
      std::array<IdFourthRoot,4> outVal;
      
      for(int i=0;i<4;i++)
	{
	  outNnc[i]=nnc[i];
	  outVal[i]=-val[i];
	}
      
      return {outNnc,outVal};
    }
    
    /// Hermitian matrix, transposed and conjugated
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr DiracGammaData dag() const
    {
      /// Result non null column
      std::array<int,4> outNnc;
      
      /// Value of non null entry
      std::array<IdFourthRoot,4> outVal;
      
      for(int sr=0;sr<4;sr++)
	{
	  const int sc=nnc[sr];
	  outNnc[sc]=sr;
	  outVal[sc]=~val[sr];
	}
      
      return {outNnc,outVal};
    }
    
    /// Constructor
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr DiracGammaData(const std::array<int,4>& nnc,
			 const std::array<IdFourthRoot,4>& val) :
      nnc{nnc},val{val}
    {
    }
    
    /// Product of two gammas
    INLINE_FUNCTION HOST_DEVICE_ATTRIB
    constexpr DiracGammaData operator*(const DiracGammaData& b) const
    {
      /// Result non null column
      std::array<int,4> outNnc;
      
      /// Value of non null entry
      std::array<IdFourthRoot,4> outVal;
      
      for(int ig1=0;ig1<4;ig1++)
	{
	  //This is the line to be taken on the second matrix
	  const int ig2=nnc[ig1];
	  
	  //The entries of the output is, on each line, the complex
	  //product of the entries of the first matrix on that line, for
	  //the entries of the second matrix on the line with the index
	  //equal to the column of the first matrix which is different
	  //from 0 (which again is ig2)
	  outVal[ig1]=val[ig1]*b.val[ig2];
	  
	  //For each line, the column of the output matrix which is
	  //different from 0 is the column of the second matrix different
	  //from 0 on the line with index equal to the column of the first
	  //matrix which is different from 0 (that is, ig2)
	  outNnc[ig1]=b.nnc[ig2];
	}
      
      return {outNnc,outVal};
    }
  };
  
  namespace SpinorialBasis
  {
    using enum IdFourthRoot;
    
    constexpr inline DiracGammaData GX{{3,2,1,0},{-I,-I,+I,+I}};
    
    constexpr inline DiracGammaData GY{{3,2,1,0},{-One,+One,+One,-One}};
    
    constexpr inline DiracGammaData GZ{{2,3,0,1},{-I,+I,+I,-I}};
    
    constexpr inline DiracGammaData G0{{2,3,0,1},{-One,-One,-One,-One}};
    
    constexpr inline DiracGammaData G5=G0*GX*GY*GZ;
    
    constexpr inline DiracGammaData ID=G0*G0;
    
    static constexpr DiracGammaData AX=GX*G5;
    
    static constexpr DiracGammaData AY=GY*G5;
    
    static constexpr DiracGammaData AZ=GZ*G5;
    
    static constexpr DiracGammaData A0=G0*G5;
    
    static constexpr DiracGammaData TX=G0*GX;
    
    static constexpr DiracGammaData TY=G0*GY;
    
    static constexpr DiracGammaData TZ=G0*GZ;
    
    static constexpr DiracGammaData BX=GY*GZ;
    
    static constexpr DiracGammaData BY=GZ*GX;
    
    static constexpr DiracGammaData BZ=GX*GY;
    
    [[maybe_unused]]
    static constexpr std::array<DiracGammaData,16> basis{ID,GX,GY,GZ,G0,G5,AX,AY,AZ,A0,TX,TY,TZ,BX,BY,BZ};
  }
  
  /////////////////////////////////////////////////////////////////
  
  /// Dirac Gamma matrix, can contain only +/-1 or +/-i, and only one
  /// nonzero column per row
  ///
  /// Forward declaration
  template <DerivedFromTransposableComp _CSo=Spin,
	    DerivedFromTransposableComp _CSi=_CSo,
	    typename _Fund=int>
  struct DiracGammaOver;
  
  PROVIDE_FEATURE(DiracGamma);
  
#define THIS					\
  DiracGammaOver<_CSo,_CSi,_Fund>
  
#define COMPS					\
  CompsList<RowOf<_CSo>,ClnOf<_CSi>,ReIm>
  
#define BASE					\
  Node<THIS,COMPS>
  
  /// Dirac Gamma matrix, can contain only +/-1 or +/-i, and only one
  /// nonzero column per row
  template <DerivedFromTransposableComp _CSo,
	    DerivedFromTransposableComp _CSi,
	    typename _Fund>
  struct DiracGammaOver :
    NoSubExprs,
    ProvideNoDynamicsComp,
    DiracGammaFeat,
    BASE
  {
#undef BASE
    
    using This=
      THIS;
    
#undef THIS
    
    const DiracGammaData data;
    
    using CSo=RowOf<_CSo>;
    
    using CSi=ClnOf<_CSi>;
    
    /// Components
    using Comps=
      COMPS;
    
    /// Fundamentl type
    using Fund=
      _Fund;
    
#undef COMPS
    
    /// Do not need to store by ref
    static constexpr bool storeByRef=
      false;
    
    /// Can run on both GPU and CPU as it is trivially copyable
    static constexpr ExecSpace execSpace=
      execOnCPUAndGPU;
    
    /// Evaluate
    HOST_DEVICE_ATTRIB INLINE_FUNCTION constexpr
    Fund eval(const DerivedFromComp auto&...cs) const
    {
      const auto c=
	std::make_tuple(cs...);
      
      const int cSo=
	std::get<CSo>(c)();
      
      const int cSi=
	std::get<CSi>(c)();
      
      const ReIm ri=
	std::get<ReIm>(c);
      
      const int nonNullCSi=data.nnc[cSo];
      
      const IdFourthRoot idR=data.val[cSo];
      
      using enum IdFourthRoot;
      
      if(nonNullCSi==cSi)
	switch(idR)
	  {
	  case One:
	    if(ri==0)
	      return 1;
	    break;
	  case MinusOne:
	    if(ri==0)
	      return -1;
	    break;
	  case I:
	    if(ri==1)
	      return 1;
	    break;
	  case MinusI:
	    if(ri==1)
	      return -1;
	    break;
	  }
      
      return 0;
    }
    
    /// Take a copy as a reference
    INLINE_FUNCTION DiracGammaOver getRef() const
    {
      return *this;
    }
    
    /// Type obtained reinterpreting the fund
    template <typename NFund>
    using ReinterpretFund=
      DiracGammaOver<CSo,CSi,NFund>;
    
    /// Default constructor
    INLINE_FUNCTION constexpr HOST_DEVICE_ATTRIB
    DiracGammaOver(const DiracGammaData& data) :
      data(data)
    {
    }
  };
  
  /// Defaulting the Dirac gamma
  using DiracGamma=
    DiracGammaOver<>;
}

#endif
