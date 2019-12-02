#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <base/field.hpp>
#include <routines/ios.hpp>

#include<complex>
namespace nissa
{
  
  void new_index_test()
  {

    typedef double comple[2];
    typedef comple color[NCOL];
    typedef color spincolor[NDIRAC];
    
    field<spincolor> a;
    
    asm("#here1");
    for(int ispin=0;ispin<4;ispin++)
      for(int icol=0;icol<3;icol++)
	{
	  a(ispin,icol)=a.index(ispin,icol);
	}
    asm("#here2");
    
    master_printf("%lg\n",a(0,0,0));
    
    crash("");
  }
}
