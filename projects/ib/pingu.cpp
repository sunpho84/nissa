#include <nissa.hpp>

using namespace nissa;

namespace nissa
{
  DECLARE_TRANSPOSABLE_COMP(SpinCc,int,4,spinCc);
  DECLARE_TRANSPOSABLE_COMP(SpinHsl,int,4,spinHsl);
  DECLARE_TRANSPOSABLE_COMP(ColCc,int,3,colCc);
  DECLARE_TRANSPOSABLE_COMP(ColHsl,int,3,colHsl);
}

double res(const StackTens<CompsList<SpinRow,SpinCln>>& in,
	   const DiracGamma& G)
{
  ASM_BOOKMARK_BEGIN("res");
  const double e=real(trace(in*G));
  ASM_BOOKMARK_END("res");
  return e;
}

double res2(const StackTens<CompsList<SpinRow,SpinCln>>& in)
{
  using namespace SpinorialBasis;
  
  ASM_BOOKMARK_BEGIN("res2");
  const double e=real(trace(in*DiracGamma(ID)));
  ASM_BOOKMARK_END("res2");
  return e;
}

double res3(const StackTens<CompsList<SpinRow,SpinCln>>& in)
{
  using namespace SpinorialBasis;
  
  ASM_BOOKMARK_BEGIN("res3");
  double e=0;
  for(SpinRow srow=0;srow<4;srow++)
    {
      const SpinCln scln=DiracGamma(ID).data.nnc[srow()];
      const auto v=DiracGamma(ID).data.val[srow()];
      
      switch(v)
	{
	case IdFourthRoot::One:
	  e+=in[srow][scln];
	break;
	case IdFourthRoot::I:
	  break;
	case IdFourthRoot::MinusOne:
	  e-=in[srow][scln];
	  break;
	case IdFourthRoot::MinusI:
	  break;
      }
    }
  ASM_BOOKMARK_END("res3");
  
  return e;
}

double ys{78};

double res4(double x)
{
 ASM_BOOKMARK_BEGIN("res4");

 constexpr int i=0;
 double y;
 if(not std::isfinite(x))
   y=0;
 else y=i*x;
 ASM_BOOKMARK_END("res4");
 
 return y;
}

constexpr int nHits=12,T=128,nContr=9;

namespace nissa
{
  DECLARE_UNTRANSPOSABLE_COMP(Contr,int,nContr,contr);
  DECLARE_UNTRANSPOSABLE_COMP(Hit,int,nHits,hit);
  DECLARE_DYNAMIC_COMP(TimeCc);
  DECLARE_DYNAMIC_COMP(TimeHsl);
}

template <DerivedFromComp Time,
	  DerivedFromComp SpiRow,
	  DerivedFromComp ColRow,
	  DerivedFromComp SpiCln,
	  DerivedFromComp ColCln>
auto load(const std::string& tag)
{
  DynamicTens<OfComps<Contr,Hit,Time,SpiCln,ColCln,SpiRow,ColRow,ReIm>> res(std::make_tuple(Time{T}));
  
  FILE *fin=fopen(("/home/francesco/QCD/LAVORI/PINGU/data/mes_contr_"+tag).c_str(),"r");
  if(fin==nullptr)
    CRASH("Unable to open");
  
  compsLoop<CompsList<SpiRow,
		      ColRow,
		      SpiCln,
		      ColCln>>([fin,
				&res](const DerivedFromComp auto...cs)
		      {
			for(Hit iHit=0;iHit<nHits;iHit++)
			  {
			    auto exp=
			      [fin](const auto...str)
			      {
				char a[100];
				for(auto e : {str...})
				  if(int rc=fscanf(fin,"%s",a)!=1;rc or strcmp(e,a))
				    CRASH("rc=%d, or %s!=%s",rc,e,a);
			      };
			    
			    auto get=
			      [fin](const char* fmt,auto out)
			      {
				if(int rc=fscanf(fin,fmt,out)!=1)
				  CRASH("rc=%d for %s",rc,fmt);
			      };
			    
			    exp("#","Contraction","of");
			    char a[100];
			    get("%s",a);
			    exp("^","\\dag","and");
			    char b[100];
			    get("%s",b);
			    exp("origin","located","at","txyz","=");
			    char c[100];
			    get("%s",c);
			    
			    char g[100];
			    for(Contr iContr=0;iContr<nContr;iContr++)
			      {
				exp("#");
				get("%s",g);
				// masterPrintf("%s\n",g);
				
				for(Time t=0;t<T;t++)
				  {
				    for(ReIm reIm=0;reIm<2;reIm++)
				      if(fscanf(fin,"%lg",&res(iHit,iContr,cs...,t,reIm))!=1)
					CRASH("Unable to get reim for t %ld",t());
				  }
			      }
			  }
		      },{});
  
  char dum[100];
  if(fscanf(fin,"%s",dum)==1)
    CRASH("Unexpected to be able to read something, obtained %s",dum);
  
  fclose(fin);
  
  return res;
}

void gammaTest()
{
  using namespace SpinorialBasis;
  
  for(const auto& [name,G] : std::vector<std::pair<const char*,DiracGammaData>>{std::make_pair("Id",ID),{"GX",GX},{"GY",GY},{"GZ",GZ},{"G0",G0},{"G5",G5}})
    {
      masterPrintf("%s\n",name);
      
      const DiracGamma DG(G);
      
      for(SpinRow sr=0;sr<4;sr++)
	{
	  for(SpinCln sc=0;sc<4;sc++)
	    masterPrintf("%+02.2f,%+02.2f\t",(double)DG(sr,sc,ReIm(0)),(double)DG(sr,sc,ReIm(1)));
	  masterPrintf("\n");
	}
      
      masterPrintf("\n\n");
    }
}

void in_main(int narg,char **arg)
{
  //gammaTest();
  
  using namespace SpinorialBasis;
  
  const auto hsl=load<TimeHsl,SpinHslCln,ColHslCln,SpinHslRow,ColHslRow>("HSL").contr(0).close();
  const auto cc=load<TimeCc,SpinCcCln,ColCcCln,SpinCcRow,ColCcRow>("CC").contr(0).close();
  const auto v0=(compAverage<Hit>(trace(DiracGammaOver<SpinHsl>(G5)*hsl))).close();
  for(TimeHsl tHsl=0;tHsl<T;tHsl++)
    {
      masterPrintf("%ld ",tHsl());
      for(ReIm reIm=0;reIm<2;reIm++)
	masterPrintf("%lg ",v0(tHsl,reIm));
          masterPrintf("\n");
    }
  
  const auto dccd=(DiracGammaOver<SpinHsl,SpinCc>(G5*GZ)*cc*DiracGammaOver<SpinCc,SpinHsl>(G5*GZ)).close();
  
 // ASM_BOOKMARK_BEGIN("otto");
 //  auto otto=
 //    (compAverage<Hit>(trace(hsl*kronDelta<ColHslRow,ColCcCln>*dccd*kronDelta<ColCcRow,ColHslCln>))).close();
 // ASM_BOOKMARK_END("otto");
  
  for(TimeCc tCc=0;tCc<T;tCc++)
    for(TimeHsl tHsl=0;tHsl<T;tHsl++)
    {
      masterPrintf("%ld %ld ",tCc(),tHsl());
      for(ReIm reIm=0;reIm<2;reIm++)
	//masterPrintf("%lg ",otto(tCc,tHsl,reIm));
	masterPrintf("%lg ",dccd(tCc,reIm).hit(0).colCcRow(0).colCcCln(0).spinHslCln(0).spinHslRow(0));
      masterPrintf("\n");
    }
  
  // constexpr auto a=
  //   []()
  // {
  //   StackTens<CompsList<SpinRow,SpinCln>> res;
    
  //   for(SpinRow i=0;i<4;i++)
  //     for(SpinCln j=0;j<4;j++)
  // 	res(i,j)=i()+4*j();
    
  //   return res;
  // }();

  // const double e2=res2(a);
  // const double e3=res3(a);
  
#if 0
  const int T=64;
  TimeTbs tBS(T);
  TimeCc tCC(T);
  
  DynamicTens<OfComps<TimeCc,SpinCcRow,ColCcRow,SpinCcCln,ColCcCln,ReIm>,double> CC(tCC);
  DynamicTens<OfComps<TimeTbs,SpinTbsRow,ColTbsRow,SpinTbsCln,ColTbsCln,ReIm>,double> BS(tBS);
  
  CC(spinCc(0));
  
  auto CCc=traceOver<SpinCc>(traceOver<ColCc>(CC));
  auto BSc=traceOver<SpinTbs>(traceOver<ColTbs>(BS));
  
  DynamicTens<OfComps<TimeCc,TimeTbs,ReIm>,double> c;

  constexpr int a=decltype(BS*DiracGammas<SpinTbs,SpinCc>::GX*CC*DiracGammas<SpinCc,SpinTbs>::GX*getKronDelta<ColCcRow,ColTbsCln>()
			   )::Comps{};
  
  
#endif
  
  // constexpr auto i=
  // 	      DiracGamma(SpinorialBasis::GX).close();

  
}

int main(int narg,char **arg)
{
  initNissa(narg,arg);
  in_main(narg,arg);
  closeNissa();
  
  return 0;
}
