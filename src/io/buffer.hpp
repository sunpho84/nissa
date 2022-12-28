#ifndef _BUFFER_HPP
#define _BUFFER_HPP

#include <sstream>

#include "io/endianness.hpp"

namespace nissa
{
  struct buffer_t
  {
    //return size
    int size()
    {
      int cur=st.tellg();
      st.seekg(0,std::ios::end);
      int s=st.tellg();
      st.seekg(cur,std::ios::beg);
      return s;
    }
    std::ostream& write(const char *s,std::streamsize n){return st.write(s,n);}
    std::istream& read(char *s,std::streamsize n){return st.read(s,n);}
    bool operator!(){return !st;}
  private:
    std::stringstream st;
  };
  
  template <class T> buffer_t &operator<<(buffer_t &out,
					  T in)
  {
    fixFromNativeEndianness<LittleEndian>(in);
    out.write((char*)&in,sizeof(T));
    
    return out;
  }
  
  template <class T>
  buffer_t &operator>>(buffer_t &in,T &out)
  {
    in.read((char*)&out,sizeof(T));
    fixToNativeEndianness<LittleEndian>(out);
    
    return in;
  }
}

#endif
