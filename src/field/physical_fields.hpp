#ifndef _PHYSICAL_FIELDS_HPP
#define _PHYSICAL_FIELDS_HPP

#include "field.hpp"

namespace nissa
{
  /// Define components
  ///
  /// Name is obtained combining NAME and Comps. Components are
  /// obtained from variadic args.
#define DEFINE_COMPS(NAME,...)			\
  /*! Components carried by a NAME */			\
  using NAME ## Comps=TensComps<__VA_ARGS__>
  
  /// Define a tensor
  ///
  /// Name is obtained combining NAME and SHORT_FUND.
  /// Components are obtained combining NAME and Comps.
  /// Fundamental type is FUND.
#define DEFINE_TENS(NAME,SHORT_FUND,FUND)		\
  /*! NAME tensor of FUND */				\
  using NAME ## SHORT_FUND=Tens<NAME ## Comps,FUND>
  
  /// Define a tensor with many fundamental type
#define DEFINE_TENS_WITH_DOUBLE_FLOAT(NAME)		\
  DEFINE_TENS(NAME,D,double);				\
  DEFINE_TENS(NAME,F,float);				\
  using NAME = NAME ## D
  
  /// Define a field
  ///
  /// Name is obtained combining NAME, SHORT_LOC and SHORT_FUND.
  /// Field run over LOC sites, has components combination of NAME
  /// and COMPS, and fundamental type FUND.
#define DEFINE_FIELD(NAME,SHORT_LOC,LOC,SHORT_FUND,FUND)		\
  /*! Field of type NAME over LOC volume with FUND type*/		\
  using NAME ## SHORT_LOC ## SHORT_FUND = Field<LOC,NAME ## Comps,FUND>
  
  /// Define a field with many fundamental type
#define DEFINE_FIELD_WITH_DOUBLE_FLOAT(NAME,SHORT_LOC,LOC)		\
  DEFINE_FIELD(NAME,SHORT_LOC,LOC,D,double);				\
  DEFINE_FIELD(NAME,SHORT_LOC,LOC,F,float);				\
  using NAME ## SHORT_LOC = NAME ## SHORT_LOC ## D
  
  /// Define lx, even and odd, double and float version of the field
#define DEFINE_FIELD_WITH_LX_EVN_ODD_DOUBLE_FLOAT(NAME)		\
  DEFINE_FIELD_WITH_DOUBLE_FLOAT(NAME,Lx,LocVolIdx);		\
  DEFINE_FIELD_WITH_DOUBLE_FLOAT(NAME,Evn,LocVolEvnIdx);	\
  DEFINE_FIELD_WITH_DOUBLE_FLOAT(NAME,Odd,LocVolOddIdx);	\
  /*! Alias for the lx double precision version of the field */\
  using NAME ## Field = NAME ## Lx
  
  /// Define a "physical field" with name NAME based on the variadic components
#define DEFINE_PHYSICAL_FIELD(NAME,...)				     \
  DEFINE_COMPS(NAME,__VA_ARGS__);				     \
  DEFINE_TENS_WITH_DOUBLE_FLOAT(NAME);				     \
  DEFINE_FIELD_WITH_LX_EVN_ODD_DOUBLE_FLOAT(NAME)
  
  DEFINE_PHYSICAL_FIELD(Scalar,);
  DEFINE_PHYSICAL_FIELD(Compl,ComplIdx);
  DEFINE_PHYSICAL_FIELD(ColorCompl,ColorIdx<>,ComplIdx);
  DEFINE_PHYSICAL_FIELD(SpinCompl,SpinIdx<>,ComplIdx);
  DEFINE_PHYSICAL_FIELD(SpinColorCompl,SpinIdx<>,ColorIdx<>,ComplIdx);
  DEFINE_PHYSICAL_FIELD(ColorColorCompl,ColorIdx<ROW>,ColorIdx<COL>,ComplIdx);
  DEFINE_PHYSICAL_FIELD(SpinSpinCompl,SpinIdx<ROW>,SpinIdx<COL>,ComplIdx);
  DEFINE_PHYSICAL_FIELD(LorentzColorColorCompl,LorentzIdx<ROW>,ColorIdx<ROW>,ColorIdx<COL>,ComplIdx);
  
  /// Alias for LorentzColorColorCompl
  using QuadSU3=LorentzColorColorCompl;
  
  /// Alias for LorentzColorColorComplField
  using GaugeConf=LorentzColorColorComplField;
  
  /// Alias for even part of the gauge conf (LorentzColorColorComplField
  using EvnGaugeConf=LorentzColorColorComplEvnD;
  
  /// Alias for odd part of the gauge conf (LorentzColorColorComplField
  using OddGaugeConf=LorentzColorColorComplOddD;
  
  /// Real part access to complex
  constexpr ComplIdx re{0};
  
  /// Imaginary part access to complex
  constexpr ComplIdx im{1};
  
}

#endif
