#include <nissa.hpp>

using namespace nissa;

namespace nissa
{
  
  
    {
      
	{
    }
  ASM_BOOKMARK_END("res3");
  
}

double ys{78};
{


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
