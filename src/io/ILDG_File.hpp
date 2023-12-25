#ifndef _ILDG_FILE_HPP
#define _ILDG_FILE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <stdio.h>
#include <stdint.h>

#include <optional>
#include <string>
#include <sstream>
#include <vector>

#include <io/checksum.hpp>
#include <base/debug.hpp>
#include <base/lattice.hpp>
#include <geometry/geometry_lx.hpp>

namespace nissa
{
  typedef off_t ILDGOffset;
  
  inline bool ignoreIldgMagicNumber{false};
  
  PROVIDE_RESOURCE(ildgToNissaRemapper,AllToAllComm<LocLxSite,LocLxSite>);
  
  PROVIDE_RESOURCE(nissaToIldgRemapper,AllToAllComm<LocLxSite,LocLxSite>);
  
  /// ILDG header
  struct ILDGHeader
  {
    static constexpr uint16_t defaultVersion=1;
    
    static constexpr uint32_t defaultMagicNo=0x456789ab;
    
    static constexpr uint16_t MB_MASK=0x80;
    
    static constexpr uint16_t ME_MASK=0x40;
    
    /// Magic number
    uint32_t magicNo;
    
    /// Version number
    uint16_t version;
    
    /// Obscure flag
    uint16_t mbmeFlag;
    
    /// Length of the data
    uint64_t dataLength;
    
    char type[128];
    
    /////////////////////////////////////////////////////////////////
    
    /// Take the MB flag
    bool getMBFlag() const
    {
      return mbmeFlag&MB_MASK;
    }
    
    /// Take the ME flag
    bool getMEFlag() const
    {
      return mbmeFlag&ME_MASK;
    }
    
    /// Returns a tuple of reference to the data which needs to be adjusted of endianness
    auto tieEndiannessVariableData()
    {
      return std::tie(dataLength,magicNo,version);
    }
    
    /// Fix the endianness to the native
    void fixToNativeEndianness()
    {
      std::apply([](auto&...args)
      {
	(nissa::fixToNativeEndianness<BigEndian>(args),...);
      },tieEndiannessVariableData());
    }
    
    /// Fix the endianness from the native
    void fixFromNativeEndianness()
    {
      std::apply([](auto&...args)
      {
	(nissa::fixFromNativeEndianness<BigEndian>(args),...);
      },tieEndiannessVariableData());
    }
    
    /// Build record header
    ILDGHeader(const char *extType,
	       const uint64_t& dataLength) :
      magicNo{defaultMagicNo},
      version{defaultVersion},
      mbmeFlag{},
      dataLength{dataLength}
    {
      strncpy(type,extType,sizeof(type)-1);
    }
    
    ILDGHeader()
    {
    }
  };
  
  /// ILDG messages
  using ILDGMessages=
    std::map<std::string,std::vector<char>>;
  
  /// Structure to manipulate Ildg files
  struct ILDGFile
  {
    static constexpr char scidacChecksumRecordName[]=
      "scidac-checksum";
    
    /// Internal handle
    FILE* file;
    
    /// Open the given file
    void open(const std::string &path,
	      const char *mode)
    {
      file=fopen(path.c_str(),mode);
      if(file==nullptr)
	crash("while opening file %s",path.c_str());
    }
    
    /// Constructor
    ILDGFile(const std::string &path,
	      const char *mode)
    {
      open(path,mode);
    }
    
    /// Default constructor
    ILDGFile() :
      file{nullptr}
    {
    }
    
    /// Gets current position
    ILDGOffset getPosition() const
    {
      return ftell(file);
    }
    
    /// Set position
    void setPosition(const ILDGOffset& pos,
		     const int& amode)
    {
      crash_printing_error(fseek(file,pos,amode),"while seeking");
      
      mpiRanksBarrier();
    }
    
    /// Gets the total size
    ILDGOffset getSize()
    {
      const ILDGOffset oriPos=
	getPosition();
      
      setPosition(0,SEEK_END);
      
      const ILDGOffset size=
	getPosition();
      
      setPosition(oriPos,SEEK_SET);
      
      return size;
    }
    
    /// Close an open file
    void close()
    {
      crash_printing_error(fclose(file),"while closing file");
      
      mpiRanksBarrier();
      
      file=nullptr;
    }
    
    /// Destructor
    ~ILDGFile()
    {
      if(file)
	close();
    }
    
    /// Skip the passed amount of bytes starting from curr position
    void skipNbytes(const ILDGOffset& nBytes)
    {
      if(nBytes)
	crash_printing_error(fseek(file,nBytes,SEEK_CUR),"while seeking ahead %ld bytes from current position",nBytes);
    
      mpiRanksBarrier();
    }
    
    /// Seek to position corresponding to next multiple of eight
    void seekToNextEightMultiple()
    {
      skipNbytes(diffWithNextMultipleOf<8>(getPosition()));
    }
    
    /// Check if end of file reached
    bool reachedEOF()
    {
      return getSize()==getPosition();
    }
    
    /// Simultaneous read from all node
    template <typename T>
    void readAll(T& data,
		 const size_t& nbytesReq=sizeof(T))
    {
      seekToNextEightMultiple();
      
      const size_t nbytesRead=
	fread(&data,1,nbytesReq,file);
      
      if(nbytesRead!=nbytesReq)
	crash("read %u bytes instead of %u required",nbytesRead,nbytesReq);
      
      //padding
      seekToNextEightMultiple();
      
      verbosity_lv3_master_printf("record read: %u bytes\n",nbytesReq);
    }
    
    /// Search next record
    ILDGHeader getNextRecordHeader()
    {
      //padding
      seekToNextEightMultiple();
      
      ILDGHeader header;
      readAll(header);
      header.fixToNativeEndianness();
      
      verbosity_lv3_master_printf("record %s contains: %ld bytes\n",header.type,header.dataLength);
      
      if(header.magicNo!=header.defaultMagicNo)
	{
	  char buf[1024];
	  snprintf(buf,1024,"wrong magic number, expected %x and obtained %x",header.defaultMagicNo,header.magicNo);
	  
	  if(ignoreIldgMagicNumber)
	    master_printf("Warning, %s\n",buf);
	  else
	    crash(buf);
	}
      
      return header;
    }
    
    /// Skip the current record
    void skipRecord(const ILDGHeader& header)
    {
      seekToNextEightMultiple();
      skipNbytes(header.dataLength);
      seekToNextEightMultiple();
    }
    
    /// Write from first node
    template <typename T>
    void masterWrite(const T& data,
		     const size_t nBytesReq=sizeof(T))
    {
      if(isMasterRank())
	{
	  if(const size_t nBytesWritten=fwrite(&data,1,nBytesReq,file);
	     nBytesWritten!=nBytesReq)
	    crash("wrote %lu bytes instead of %lu required",nBytesWritten,nBytesReq);
	  
	  //this is a blocking routine!
	  skipNbytes(0);
	}
      else
	skipNbytes(nBytesReq);
      
      //sync
      fflush(file);
      mpiRanksBarrier();
    }
    
    /// Write the header taking into account endianness
    void writeRecordHeader(ILDGHeader header)
    {
      header.fixFromNativeEndianness();
      
      masterWrite(header);
    }
    
    /// Write a record
    template <typename T>
    void writeRecord(const char *type,
		     const T& data,
		     const uint64_t& len=sizeof(T))
    {
      writeRecordHeader({type,len});
      
      masterWrite(data,len);
      
      seekToNextEightMultiple();
    }
    
    /// Special write for strings
    void writeTextRecord(const char *type,
			 const char *text)
    {
      writeRecord(type,text,strlen(text)+1);
    }
    
    /// Write the checksum
    void writeChecksum(const Checksum& check)
    {
      /// String to be written
      char mess[160];
      
      snprintf(mess,
	       159,
	       "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
	       "<scidacChecksum>"
	       "<version>1.0</version>"
	       "<suma>%#010x</suma>"
	       "<sumb>%#010x</sumb>"
	       "</scidacChecksum>",
	       check[0],
	       check[1]);
      
      writeTextRecord(scidacChecksumRecordName,mess);
    }
    
    static void initIldgNissaRemapper()
    {
      resources::_ildgToNissaRemapper.init(lat->getLocVol(),
					   [lat=lat->getRef(),
					    nRanksPerDir=nRanksPerDir.getReadable()](const LocLxSite ildgChunkEl)
					   {
					     const LocLxSite ildgEl=
					       thisRank()*lat.getLocVol()+ildgChunkEl;
					     
					     const GlbCoords ildgSizes=
					       Lattice::scidacRemap(lat.getGlbSizes());
					     
					     const GlbCoords ildgGlbCoords=
					       decomposeLxToCoords(ildgEl,ildgSizes);
					     
					     const GlbCoords nissaGlbCoords=
					       Lattice::scidacRemap(ildgGlbCoords);
					     
					     const MpiRankCoords rankCoords=
					       (nissaGlbCoords.template reinterpretFund<LocCoord>()/lat.getLocSizes()).reinterpretFund<MpiRankCoord>();
					     
					     const MpiRank rank=
					       lxOfCoords<MpiRank>(rankCoords,nRanksPerDir);
					     
					     const LocCoords nissaLocCoords=
					       nissaGlbCoords.template reinterpretFund<LocCoord>()-
					       rankCoords.template reinterpretFund<LocCoord>()*lat.getLocSizes();
					     
					     const LocLxSite locNissaSite=
					       lxOfCoords<LocLxSite>(nissaLocCoords,lat.getLocSizes());
					     
					     return std::make_tuple(rank,locNissaSite);
					   });
    }
    
    /// Read the data in the ILDG data format, according to the header
    template <DerivedFromComp...C,
	      typename Fund>
    void readFieldAndReorder(Field<OfComps<C...>,Fund,FieldLayout::CPU,MemoryType::CPU>& data,
			     const ILDGHeader& header)
    {
      if(const size_t& nF=data.nElements,
	 &nH=header.dataLength;
	 nF<header.dataLength)
	crash("Field has %zu elements, smaller than needed %zu",nF,nH);
      
      const ILDGOffset nbytesPerRankExp=
	header.dataLength/nRanks();
      
      /// Original position
      const ILDGOffset oriPos=
	getPosition();
      
      /// Starting point
      const ILDGOffset newPos=
	oriPos+thisRank()*nbytesPerRankExp;
      setPosition(newPos,SEEK_SET);
      
      //read
      if(const ILDGOffset nbytesRead=
	 fread(data,1,nbytesPerRankExp,file);
	 nbytesRead!=nbytesPerRankExp)
	crash("read %ld bytes instead of %ld",nbytesRead,nbytesPerRankExp);
      
      //place at the end of the record, including padding
      setPosition(oriPos+ceilToNextMultipleOf<8>(header.dataLength),SEEK_SET);
      
      crash("");
      //reorder data to the appropriate place
      // vector_remap_t *rem=new vector_remap_t(locVol,index_from_ILDG_remapping);
      // rem->remap(data,buf,header.data_length/glbVol);
      // delete rem;
      
      verbosity_lv3_master_printf("ildg data record read: %ld bytes\n",header.dataLength);
    }
    
    /// Search a particular record in a file
    std::optional<ILDGHeader> searchRecord(const std::string& recordName,
					   ILDGMessages* mess=nullptr)
    {
      while(not reachedEOF())
	{
	  ILDGHeader header=
	    getNextRecordHeader();
	  
	  verbosity_lv3_master_printf("found record: %s\n",header.type);
	  
	  if(strcmp(recordName.c_str(),header.type)==0)
	    return header;
	  else
	    if(mess==nullptr)
	      skipRecord(header); //ignore message
	    else
	      {
		//load the message and pass to next
		auto [data,inserted]=
		  mess->emplace(header.type,header.dataLength);
		readAll(*data,header.dataLength);
	      }
	}
      
      return {};
    }
    
    /// Read the checksum
    Checksum readScidacChecksum()
    {
      /// Setup as non found as search it
      Checksum checkRead{};
      
      if(const auto header=
	 searchRecord(scidacChecksumRecordName);header)
	{
	  const uint64_t nBytes=
	    header->dataLength;
	  
	  /// Message read
	  char mess[nBytes+1];
	  readAll(*mess,nBytes);
	  mess[nBytes]='\0';
	  
	  for(const auto& [name,id] :
		{std::make_pair("<suma>",0),
		 std::make_pair("<sumb>",1)})
	    if(char *handle=
	       strstr(mess,name);
	       handle==nullptr or
	       sscanf(handle+6,"%x",&checkRead[0])==0)
	      master_printf("WARNING: Broken checksum\n");
	}
      
      return checkRead;
    }
    
    // /// Remap to ildg
    // void remapToWriteIldgData(char* buf,char* data,int nbytes_per_site)
    // {
    //   crash("reimplement");
      // PAR(0,locVol,
      // 	CAPTURE(),
      // 	ivol,
      // 	{
      // 	  int64_t idest=0;
      // 	  for(int mu=0;mu<NDIM;mu++)
      // 	    {
      // 	      int nu=scidac_mapping[mu];
      // 	      idest=idest*locSize[nu]+locCoordOfLoclx[isour][nu];
      // 	    }
      // 	  memcpy(buf+nbytes_per_site*idest,data+nbytes_per_site*isour,nbytes_per_site);
      // 	});
    // }
    
    // /// Bare write data in the ILDG order
    // void writeIldgDataAllRaw(void *data,
    // 			     const uint64_t& dataLength)
    // {
    //   crash("reimplement");
      
      // //allocate the buffer
      // ILDG_Offset nbytes_per_rank=data_length/nranks;
      // char *buf=nissa_malloc("buf",nbytes_per_rank,char);
      
      // //take original position
      // ILDG_Offset ori_pos=ILDG_File_get_position(file);
      
      // //find starting point
      // ILDG_Offset new_pos=ori_pos+rank*nbytes_per_rank;
      // ILDG_File_set_position(file,new_pos,SEEK_SET);
      
      // //reorder data to the appropriate place
      // vector_remap_t *rem=new vector_remap_t(locVol,index_to_ILDG_remapping,NULL);
      // rem->remap(buf,data,data_length/glbVol);
      // delete rem;
      
      // //write
      // ILDG_Offset nbytes_wrote=fwrite(buf,1,nbytes_per_rank,file);
      // if(nbytes_wrote!=nbytes_per_rank) crash("wrote %d bytes instead of %d",nbytes_wrote,nbytes_per_rank);
      
      // //place at the end of the record, including padding
      // ILDG_File_set_position(file,ori_pos+ceil_to_next_eight_multiple(data_length),SEEK_SET);
      
      // //free buf and ord
      // nissa_free(buf);
    // }
    
    // /// Read the data according to ILDG mapping
    // void writeIldgDataAll(void *data,
    // 			  const ILDGOffset& nbytesPerSite,
    // 			  const char *type)
    // {
    //   //prepare the header and write it
    //   const uint64_t dataLength=
    // 	nbytesPerSite*lat->getGlbVol()();
      
    //   writeRecordHeader({type,dataLength});
      
    //   writeIldgDataAllRaw(data,dataLength);
    // }
    
    // void writeAllMessages(const ILDGMessages& messages)
    // {
    //   for(const auto& [name,data] : messages)
    // 	writeRecord(name.c_str(),data[0],data.size());
    // }
  };
  
  /// Define the remapping from the layout having in each rank a
  /// consecutive block of data holding a consecutive piece of ildg
  /// data to canonical lx
  inline std::pair<int,int> index_from_ILDG_remapping(const int& iloc_ILDG)
  {
    crash("");
    // int iglb_ILDG=thisRank()*locVol+iloc_ILDG;
    
    // //find global coords in ildg ordering
    // coords_t xto;
    // for(int mu=NDIM-1;mu>=0;mu--)
    //   {
    // 	int nu=scidac_mapping[mu];
    // 	xto[nu]=iglb_ILDG%glbSizes[nu];
    // 	iglb_ILDG/=glbSizes[nu];
    //   }
    
    // return get_loclx_and_rank_of_coord(xto);
    return {};
  }
  
  /// Defines the reampping from lx in order to have in each rank a
  /// consecutive block of data holding a consecutive piece of ildg
  /// data
  inline std::pair<int,int> index_to_ILDG_remapping(const int& iloc_lx)
  {
    crash("");
    // //find global index in ildg transposed ordering
    // int iglb_ILDG=0;
    // for(int mu=0;mu<NDIM;mu++)
    //   {
    // 	const int nu=scidac_mapping[mu];
    // 	iglb_ILDG=iglb_ILDG*glbSizes[nu]+glbCoordOfLoclx[iloc_lx][nu];
    //   }
    
    // const int irank_ILDG=iglb_ILDG/locVol;
    // const int iloc_ILDG=iglb_ILDG%locVol;
    
    // return {irank_ILDG,iloc_ILDG};
    
    return {};
  }
  
  // //Writes a field to a file (data is a vector of loc_vol) with no frill
  // template <typename T>
  // void write_lattice_field(ILDGFile &file,T *data)
  // {
  //   ILDGFile_write_ildg_data_all_raw(file,data,lat->getLocVol()()*sizeof(T));
  // }
  
  // //Writes a field opening the file with given path (data is a vector of loc_vol) with no frill
  // template <typename T>
  // void write_lattice_field(const char *path,T *data)
  // {
  //   ILDGFile file(path,"w");
    
  //   write_lattice_field(file,data);
    
  //   ILDGFile_close(file);
  // }
  
  // //storable vector
  // ILDG_message* ILDG_string_message_append_to_last(ILDG_message *mess,const char *name,const char *data);
  // template<class T> struct storable_vector_t : std::vector<T>
  // {
  //   //append to last message
  //   ILDG_message *append_to_message_with_name(ILDG_message &mess,const char *name)
  //   {
  //     std::ostringstream os;
  //     os.precision(16);
  //     for(typename std::vector<T>::iterator it=this->begin();it!=this->end();it++) os<<*it<<" ";
  //     return ILDG_string_message_append_to_last(&mess,name,os.str().c_str());
  //   }
  //   //convert from a text message
  //   void convert_from_text(const char *data)
  //   {
  //     std::istringstream is(data);
  //     T temp;
  //     while(is>>temp) this->push_back(temp);
  //   }
  //   void convert_from_message(ILDG_message &mess)
  //   {convert_from_text(mess.data);}
    
  //   //read it from file
  //   void read_from_ILDG_file(ILDGFile fin, const char *tag)
  //   {
  //     ILDG_header head;
  //     head=ILDGFile_get_next_record_header(fin);
  //     if(strcasecmp(tag,head.type)==0)
  // 	{
  // 	  char *data=new char[head.data_length+1];
  // 	  ILDGFile_read_all(data,fin,head.data_length);
  // 	  this->convert_from_text(data);
  // 	  delete[] data;
  // 	}
  //     else crash("Unable to convert, tag %d while expecting %d",head.type,tag);
  //   }
  // };
}

#endif
