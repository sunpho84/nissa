#ifndef _NICE_HOST_MEMBERS_HPP
#define _NICE_HOST_MEMBERS_HPP

#include <memory>
#include <string>
#include <tuple>

struct HObjectInfo;

/// Name and member pointer shared by script registration and parameter inspection
template <typename T,
          typename M>
struct HostMemberInfo
{
  /// Script-visible member name
  std::string name;
  
  /// Pointer to the concrete member
  M T::* pointer;
  
  /// Optional metadata for a nested host object
  std::shared_ptr<const HObjectInfo> memberType;
};

/// Describes a member with the same name in C++ and the script
/// Use inside PROVIDE_NICE_MEMBERS, which supplies the local Self alias
#define MEMB(NAME)				\
  HostMemberInfo{#NAME,&Self::NAME}

/// Describes a nested host member using its registered type metadata
/// Requires the same Self alias as MEMB
#define HMEMB(NAME,INFO)			\
  HostMemberInfo{#NAME,&Self::NAME,(INFO)}

/// Defines static members() from descriptors such as MEMB(field)
/// HMEMB and explicit HostMemberInfo descriptors may also be included in the list
#define PROVIDE_NICE_MEMBERS(TYPE,...)	 \
  static auto members()			 \
  {					 \
    using Self=TYPE;			 \
    return std::make_tuple(__VA_ARGS__); \
  }

#endif
