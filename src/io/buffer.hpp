#ifndef _BUFFER_HPP
#define _BUFFER_HPP

#include <sstream>

#include "io/endianness.hpp"

namespace nissa
{
  struct buffer_t : public std::stringstream
  {
    //return size
    int size()
    {
      int cur=tellg();
      seekg(0,std::ios::end);
      int s=tellg();
      seekg(cur,std::ios::beg);
      
      return s;
    }
  };
  template <class T> buffer_t &operator<<(buffer_t &out,T &in)
  {
    if(!little_endian)
      {
        T temp=in;
        change_endianness(temp);
        out.write((char*)&temp,sizeof(T));
      }
    else out.write((char*)&in,sizeof(T));
    
    return out;
  }
  
  template <class T> buffer_t &operator>>(buffer_t &in,T &out)
  {
    in.read((char*)&out,sizeof(T));
    if(!little_endian) change_endianness(out);
    
    return in;
  }
}

#endif
