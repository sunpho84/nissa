#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_REDUCE
 #include <linalgs/reduce.hpp>

namespace nissa
{
  /// Releases the reduction buffer if allocated
  void deallocate_reduction_buffer()
  {
    if(reducing_buffer!=nullptr)
      {
	master_printf("Freeing reduction buffer, used size: %ld\n",reducing_buffer_size);
	reducing_buffer_size=0;
	
	nissa_free(reducing_buffer);
      }
  }
}
