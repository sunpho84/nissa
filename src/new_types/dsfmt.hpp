#ifndef _DSFMT_HPP
#define _DSFMT_HPP

#define DSFMT_N           191
#define DSFMT_N64         382
#define DSFMT_POS1        117

#include <stdint.h>
#include <string.h>

#include "io/endianness.hpp"

namespace nissa
{
  //new generator
  class dsfmt_t
  {
    //to access internal binary part
    union w128_t
    {
      uint64_t u[2];
      uint32_t u32[4];
      double d[2];
    };
    
  private:
    dsfmt_t(){memset(this,0,sizeof(dsfmt_t));}
    
    //hold internal status
    w128_t status[DSFMT_N+1];
    uint64_t idx;
    uint32_t seed;
    uint64_t nextr;
    
    //fix endianness
    int idxof(int i)
    {
      if(!little_endian) return i^1;
      else return i;
    }
    
    //this is the cripitic part that update internal status
    void do_recursion(w128_t &r,w128_t &a,w128_t &b,w128_t &lung)
    {
      w128_t L=lung;
      lung.u[0]=(a.u[0]<<19)^(L.u[1]>>32)^(L.u[1]<<32)^b.u[0];
      lung.u[1]=(a.u[1]<<19)^(L.u[0]>>32)^(L.u[0]<<32)^b.u[1];
      r.u[0]=(lung.u[0]>>12)^(lung.u[0]&0x000ffafffffffb3fULL)^a.u[0];
      r.u[1]=(lung.u[1]>>12)^(lung.u[1]&0x000ffdfffc90fffdULL)^a.u[1];
    }
    
    //generate new set of number
    void generate_all_new()
    {
      w128_t lung=status[DSFMT_N];
      
      for(int i=0;i<DSFMT_N-DSFMT_POS1;i++) do_recursion(status[i],status[i],status[i+DSFMT_POS1],lung);
      for(int i=DSFMT_N-DSFMT_POS1;i<DSFMT_N;i++) do_recursion(status[i],status[i],status[i+DSFMT_POS1-DSFMT_N],lung);
      status[DSFMT_N]=lung;
      
      //reset index
      idx=0;
    }

    //internally used to skip
    void next_state()
    {
      int in_idx=(idx/2)%DSFMT_N;
      w128_t *pstate=&status[0];
      
      w128_t *lung=&pstate[DSFMT_N];
      do_recursion(pstate[in_idx],pstate[in_idx],pstate[(in_idx + DSFMT_POS1)%DSFMT_N],lung[0]);
      idx=(idx+2)%DSFMT_N64;
    }
    
    //internally used to skip
    void add(dsfmt_t &src)
    {
      int dp=idx/2;
      int sp=src.idx/2;
      int diff=(sp-dp+DSFMT_N)%DSFMT_N;
      for (int i=0;i<DSFMT_N-diff;i++)
	{
	  int p=i+diff;
	  status[i].u[0]^=src.status[p].u[0];
	  status[i].u[1]^=src.status[p].u[1];
	}
      for(int i=DSFMT_N-diff;i<DSFMT_N;i++)
	{
	  int p=i+diff-DSFMT_N;
	  status[i].u[0]^=src.status[p].u[0];
	  status[i].u[1]^=src.status[p].u[1];
	}
      status[DSFMT_N].u[0]^=src.status[DSFMT_N].u[0];
      status[DSFMT_N].u[1]^=src.status[DSFMT_N].u[1];
    }
    
  public:  
    //initialise from seed
    void init(uint32_t ext_seed,uint64_t nextr=0)
    {
      seed=ext_seed;
      nextr=0;
      
      //cryptic initialisation
      (&status[0].u32[0])[idxof(0)]=seed;
      for(int i=1;i<(DSFMT_N+1)*4;i++) (&status[0].u32[0])[idxof(i)]=1812433253UL*((&status[0].u32[0])[idxof(i-1)]^((&status[0].u32[0])[idxof(i-1)]>>30))+i;
      for(int i=0;i<DSFMT_N*2;i++) (&status[0].u[0])[i]=((&status[0].u[0])[i]&0x000FFFFFFFFFFFFFULL)|0x3FF0000000000000ULL;
      
      //certificate the period
      uint64_t inner=((status[DSFMT_N].u[0]^0x90014964b32f4329ULL)&0x3d84e1ac0dc82880ULL)^((status[DSFMT_N].u[1]^0x3b8d12ac548a7c7aULL)&0x0000000000000001ULL);
      for(int i=32;i>0;i>>=1) inner^=(inner>>i);
      if((inner&1)!=1) status[DSFMT_N].u[1]^=1;
      
      //set generated number
      idx=DSFMT_N64;
      
      if(nextr) (*this)=this->get_skipped(nextr);
    }
    
    //get a double distribute between 0 and 1
    double get()
    {
      if(idx==DSFMT_N64) generate_all_new();
      nextr++;
      return 2-(&status[0].d[0])[idx++];
    }
    
    //jump forward
    dsfmt_t get_skipped(const char *skip_string)
    {
      dsfmt_t work=(*this),out;
      int index=idx;
      work.idx=DSFMT_N64;
      
      for(int i=0;skip_string[i]!='\0';i++)
	{
	  int bits=skip_string[i];
	  if(bits>='a'&&bits<='f') bits+=10-'a';
	  else bits-='0';
	  bits&=0x0f;
	  
	  for(int j=0;j<4;j++)
	    {
	      if((bits&1)!=0) out.add(work);
	      
	      work.next_state();
	      bits>>=1;
	    }
	}
      out.idx=index;
      
      return out;
    }
    
    //skip by an integer amount
    dsfmt_t get_skipped(uint64_t skip)
    {
      uint64_t ori_skip=skip;

      //prepare output
      dsfmt_t out=(*this);
      while(out.idx+skip>=DSFMT_N64)
	{
	  int incr=DSFMT_N64-out.idx;
	  out.idx+=incr;
	  skip-=incr;
	  out.generate_all_new();
	}
      out.idx+=skip;

      //fix the number of extractions
      out.nextr+=ori_skip;
      
      return out;
    }

    //constructor
    dsfmt_t(uint32_t seed){init(seed);}
  };
  
}

#endif
