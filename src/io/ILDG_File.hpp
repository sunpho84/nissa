#ifndef _ILDG_FILE_HPP
#define _ILDG_FILE_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdio.h>
#include <stdint.h>

#include <string>
#include <sstream>
#include <vector>

#include "checksum.hpp"
#include "base/debug.hpp"
#include "geometry/geometry_lx.hpp"
#include "operations/remap_vector.hpp"

#ifndef EXTERN_ILDG
# define EXTERN_ILDG extern
# define INIT_TO(A)
#else
# define INIT_TO(A) =A
#endif

#define ILDG_MAGIC_NO                   0x456789ab
#define ILDG_MB_MASK                    ((uint16_t)0x80)
#define ILDG_ME_MASK                    ((uint16_t)0x40)

namespace nissa
{
#ifdef USE_MPI
# ifdef USE_MPI_IO
  typedef MPI_Offset ILDG_Offset;
  typedef MPI_File ILDG_File;
# else
  typedef off64_t ILDG_Offset;
  typedef FILE* ILDG_File;
# endif
#endif
  EXTERN_ILDG int ignore_ILDG_magic_number INIT_TO(false);
  EXTERN_ILDG int fast_read_write_vectors INIT_TO(false);
  
  //ILDG header
  struct ILDG_header
  {
    uint32_t magic_no;
    uint16_t version;
    uint16_t mbme_flag;
    uint64_t data_length;
    char type[128]={};
  };
  
  //store messages
  struct ILDG_message
  {
    bool is_last;
    char *data;
    char *name;
    uint64_t data_length;
    ILDG_message *next;
  };
  
  //ILDG file view
  struct ILDG_File_view
  {
    char format[100];
  };
  
  ILDG_File ILDG_File_open(const std::string &path,const char *mode);
  ILDG_File ILDG_File_open_for_read(const std::string &path);
  ILDG_File ILDG_File_open_for_write(const std::string &path);
  ILDG_File_view ILDG_File_create_scidac_mapped_view(ILDG_File &file,ILDG_Offset nbytes_per_site);
  ILDG_File_view ILDG_File_get_current_view(ILDG_File &file);
  ILDG_Offset ILDG_File_get_position(ILDG_File &file);
  ILDG_Offset ILDG_File_get_size(ILDG_File &file);
  ILDG_header ILDG_File_build_record_header(int MB_flag,int ME_flag,const char *type,uint64_t data_length);
  ILDG_header ILDG_File_get_next_record_header(ILDG_File &file);
  void ILDG_message_init_to_last(ILDG_message *mess);
  ILDG_message *ILDG_message_find_last(ILDG_message *mess);
  ILDG_message* ILDG_bin_message_append_to_last(ILDG_message *mess,const char *name,const char *data,uint64_t length);
  ILDG_message* ILDG_string_message_append_to_last(ILDG_message *mess,const char *name,const char *data);
  
  void ILDG_File_write_all_messages(ILDG_File& file,
				    const ILDG_message* mess);
  
  void ILDG_message_free_all(ILDG_message *mess);
  bool ILDG_File_reached_EOF(ILDG_File &file);
  bool get_MB_flag(ILDG_header &header);
  bool get_ME_flag(ILDG_header &header);
  int ILDG_File_search_record(ILDG_header &header,ILDG_File &file,const char *record_name,ILDG_message *mess=NULL);
  void ILDG_File_close(ILDG_File &file);
  void ILDG_File_master_write(ILDG_File &file,void *data,size_t nbytes_req);
  void ILDG_File_read_all(void *data,ILDG_File &file,size_t nbytes_req);
  Checksum ILDG_File_read_checksum(ILDG_File &file);
  
  void ILDG_File_seek_to_next_eight_multiple(ILDG_File &file);
  void ILDG_File_set_position(ILDG_File &file,ILDG_Offset pos,int amode);
  void ILDG_File_set_view(ILDG_File &file,ILDG_File_view &view);
  void ILDG_File_skip_nbytes(ILDG_File &file,ILDG_Offset nbytes);
  void ILDG_File_skip_record(ILDG_File &file,ILDG_header header);
  void ILDG_File_write_checksum(ILDG_File &file,const Checksum& check);
  
  void ILDG_File_write_record_header(ILDG_File &file,
				     const ILDG_header &header_to_write);
  
  /// Define the remapping from the layout having in each rank a
  /// consecutive block of data holding a consecutive piece of ildg
  /// data to canonical lx
  inline std::pair<int,int> index_from_ILDG_remapping(const int& iloc_ILDG)
  {
    int iglb_ILDG=rank*locVol+iloc_ILDG;
    
    //find global coords in ildg ordering
    Coords xto;
    for(int mu=NDIM-1;mu>=0;mu--)
      {
	int nu=scidacMapping[mu];
	xto[nu]=iglb_ILDG%glbSize[nu];
	iglb_ILDG/=glbSize[nu];
      }
    
    return getLoclxAndRankOfCoords(xto);
  }
  
  /// Defines the reampping from lx in order to have in each rank a
  /// consecutive block of data holding a consecutive piece of ildg
  /// data
  inline std::pair<int,int> index_to_ILDG_remapping(const int& iloc_lx)
  {
    //find global index in ildg transposed ordering
    int iglb_ILDG=0;
    for(int mu=0;mu<NDIM;mu++)
      {
	const int nu=scidacMapping[mu];
	iglb_ILDG=iglb_ILDG*glbSize[nu]+glbCoordOfLoclx[iloc_lx][nu];
      }
    
    const int irank_ILDG=iglb_ILDG/locVol;
    const int iloc_ILDG=iglb_ILDG%locVol;
    
    return {irank_ILDG,iloc_ILDG};
  }
  
  /// Write the data according to the ILDG mapping
  template <typename T>
  void ILDG_File_write_ildg_data_all(ILDG_File &file,
				     LxField<T> in,
				     const char* header_message)
  {
    /// Take note the number of reals per site
    constexpr int nrealsPerSite=
      LxField<T>::nInternalDegs;
    
    /// Take notes of the number of bytes per site
    constexpr uint64_t nBytesPerSite=
      nrealsPerSite*sizeof(typename LxField<T>::Fund);
    
    /// Computes the length of the data to be written
    const uint64_t data_length=
      nBytesPerSite*glbVol;
    
    /// Temporary buffer with CPU spacetime layout
    LxField<T,SpaceTimeLayout::CPU> temp("temp");
    
    //reorder data to the appropriate place
    vector_remap_t(locVol,index_to_ILDG_remapping).
      remap(temp.template getPtr<defaultMemorySpace>(),
	    in.template getSurelyWithSpaceTimeLayout<SpaceTimeLayout::CPU>().template getPtr<defaultMemorySpace>(),
	    data_length/glbVol);
    
    FOR_EACH_SITE_DEG_OF_FIELD(temp,
			       CAPTURE(TO_WRITE(temp)),
			       site,
			       iDeg,
			       {
				 fixFromNativeEndianness<BigEndian>(temp(site,iDeg));
			       });
    
    /// Prepare the header
    const ILDG_header header=
      ILDG_File_build_record_header(0,0,header_message,data_length);
    
    ILDG_File_write_record_header(file,header);
    
    //allocate the buffer
    ILDG_Offset nbytes_per_rank=
      header.data_length/nranks;
    
    // Take note of the original position
    const ILDG_Offset ori_pos=
      ILDG_File_get_position(file);
    
    /// Find the starting point
    const ILDG_Offset new_pos=
      ori_pos+rank*nbytes_per_rank;
    ILDG_File_set_position(file,new_pos,SEEK_SET);
    
    /// Writes
    const ILDG_Offset nbytes_wrote=
      fwrite(temp.template getSurelyReadableOn<MemorySpace::CPU>().template getPtr<MemorySpace::CPU>(),1,nbytes_per_rank,file);
    if(nbytes_wrote!=nbytes_per_rank)
      CRASH("wrote %ld bytes instead of %ld",nbytes_wrote,nbytes_per_rank);
    
    // Place at the end of the record, including padding
    ILDG_File_set_position(file,ori_pos+ceil_to_next_eight_multiple(header.data_length),SEEK_SET);
    
    // Pad if necessary
    if(const size_t pad_diff=
       header.data_length%8;pad_diff!=0)
      {
	char buf[8];
	memset(buf,0,8);
	ILDG_File_master_write(file,(void*)buf,8-pad_diff);
      }
  }
  
  void ILDG_File_write_record(ILDG_File &file,const char *type,const char *buf,uint64_t len);
  void ILDG_File_write_text_record(ILDG_File &file,const char *type,const char *text);
  
  //Writes a field to a file (data is a vector of loc_vol) with no frill
  template <typename T>
  void write_lattice_field(ILDG_File &file,T *data)
  {
    ILDG_File_write_ildg_data_all_raw(file,data,locVol*sizeof(T));
  }
  
  //Writes a field opening the file with given path (data is a vector of loc_vol) with no frill
  template <typename T>
  void write_lattice_field(const char *path,T *data)
  {
    ILDG_File file=ILDG_File_open_for_write(path);
    
    write_lattice_field(file,data);
    
    ILDG_File_close(file);
  }
  
  //storable vector
  ILDG_message* ILDG_string_message_append_to_last(ILDG_message *mess,const char *name,const char *data);
  template<class T> struct storable_vector_t : std::vector<T>
  {
    //append to last message
    ILDG_message *append_to_message_with_name(ILDG_message &mess,const char *name)
    {
      std::ostringstream os;
      os.precision(16);
      for(typename std::vector<T>::iterator it=this->begin();it!=this->end();it++) os<<*it<<" ";
      return ILDG_string_message_append_to_last(&mess,name,os.str().c_str());
    }
    //convert from a text message
    void convert_from_text(const char *data)
    {
      std::istringstream is(data);
      T temp;
      while(is>>temp) this->push_back(temp);
    }
    void convert_from_message(ILDG_message &mess)
    {convert_from_text(mess.data);}
    
    //read it from file
    void read_from_ILDG_file(ILDG_File fin, const char *tag)
    {
      ILDG_header head;
      head=ILDG_File_get_next_record_header(fin);
      if(strcasecmp(tag,head.type)==0)
	{
	  char *data=new char[head.data_length+1];
	  ILDG_File_read_all(data,fin,head.data_length);
	  this->convert_from_text(data);
	  delete[] data;
	}
      else CRASH("Unable to convert, tag %d while expecting %d",head.type,tag);
    }
  };
  
  /// Read the data according to ILDG mapping
  template <typename T>
  void ILDG_File_read_ildg_data_all(LxField<T>& out,
				    ILDG_File &file,
				    const ILDG_header &header)
  {
    // Take note of the original position
    const ILDG_Offset ori_pos=
      ILDG_File_get_position(file);
    
    /// Number of bytes expected per rank
    const ILDG_Offset nbytes_per_rank_exp=
      header.data_length/nranks;
    
    // Find the starting point and moves there
    const ILDG_Offset new_pos=
      ori_pos+rank*nbytes_per_rank_exp;
    ILDG_File_set_position(file,new_pos,SEEK_SET);
    
    // Allocated buffer on CPU with the CPU layout
    LxField<T,SpaceTimeLayout::CPU,MemorySpace::CPU> buf("buf");
    
    // Take not of the initial read time
    const double beg=
      take_time();
    
    /// Reads taking note of the number of read bytes
    const ILDG_Offset nbytes_read=
      fread(buf.template getPtr<MemorySpace::CPU>(),1,nbytes_per_rank_exp,file);
    
    // Report and verify full reading
    MASTER_PRINTF("Bare reading %zu bytes took %lg s\n",nbytes_per_rank_exp,take_time()-beg);
    if(nbytes_read!=nbytes_per_rank_exp)
      CRASH("read %zu bytes instead of %ld",nbytes_read,nbytes_per_rank_exp);
    
    // Place at the end of the record, including padding
    ILDG_File_set_position(file,ori_pos+ceil_to_next_eight_multiple(header.data_length),SEEK_SET);
    
    /// Ensures that the buffer is on the default memory space
    decltype(auto) bufOnDefaultMemorySpace=
      buf.template getSurelyReadableOn<defaultMemorySpace>();
    
    /// Incapsulates the action of remapping, which is carried out on the default memory space
    auto act=
      [&bufOnDefaultMemorySpace,
       rem=vector_remap_t(locVol,index_from_ILDG_remapping),
       &header](auto& f)
      {
	rem.remap(f.template getPtr<defaultMemorySpace>(),
		  bufOnDefaultMemorySpace.template getPtr<defaultMemorySpace>(),
		  header.data_length/glbVol);
      };
    
    if constexpr(defaultSpaceTimeLayout==SpaceTimeLayout::CPU)
      act(out);
    else
      {
	LxField<T,SpaceTimeLayout::CPU> tmp("tmp");
	act(tmp);
	out=tmp;
      }
      
      FOR_EACH_SITE_DEG_OF_FIELD(out,
				 CAPTURE(TO_WRITE(out)),
				 site,
				 iDeg,
				 {
				   fixToNativeEndianness<BigEndian>(out(site,iDeg));
				 });
      
      VERBOSITY_LV3_MASTER_PRINTF("ildg data record read: %lu bytes\n",header.data_length);
  }
}

#undef EXTERN_ILDG
#undef INIT_TO

#endif
