#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#define EXTERN_LOCALIZER
# include <operations/localizer.hpp>

namespace nissa::localizer
{
  std::tuple<FullLocDirCoord,OrthoSpaceTime> computeDimensions(const Dir& dir)
  {
    /// Fully local size in a given direction is the global size
    const FullLocDirCoord fcSize=
      lat->getGlbSizes(dir)();
    
    /// Perpendicular size across the whole lattice
    const OrthoSpaceTime glbOsdSize=
      lat->getGlbVol()()/lat->getGlbSizes(dir)();
    
    /// Portion of the perpendicular size relative to each lattice
    const OrthoSpaceTime locOsdSize=
      (glbOsdSize+nRanks()-1)/nRanks();
    
    return {fcSize,locOsdSize};
  }
  
  LocDirMaker getLocDirMaker(const Dir& dir)
  {
    return {LocLxSite(lat->getLocVol()),
	    [&dir](const LocLxSite& locLxSite)
	    {
	      /// Dimensions of the current direction
	      const auto& [fcSize,locOsdSize]=
		dimensions[dir];
	      
	      /// Coordinate in the current direction of the requires site
	      const FullLocDirCoord fc=
		lat->getGlbCoordsOfLocLx(locLxSite)[dir]();
	      
	      /// Global index in the space perpendicular to the current direction
	      int64_t glbOsd=0;
	      for(PerpDir pDir=0;pDir<NDIM-1;pDir++)
		{
		  const Dir jDir=perpDirOf[dir][pDir];
		  glbOsd=lat->getGlbSizes(jDir)()*glbOsd+lat->getGlbCoordsOfLocLx(locLxSite)(jDir)();
		}
	      
	      /// Rank hosting global site
	      const MpiRank orthoRank=
		(int)(glbOsd/locOsdSize());
	      
	      /// Local index in the rank
	      const OrthoSpaceTime locOsd=
		glbOsd-orthoRank()*locOsdSize();
	      
	      /// Merges the components
	      const auto mc=
		MC::merge(std::make_tuple(fcSize,locOsdSize),fc,OrthoSpaceTime(locOsd));
	      
	      return std::make_tuple(orthoRank,mc);
	    }};
  }
  
  void init()
  {
    if(initialized)
      crash("already initialized");
    
    for(Dir dir=0;dir<NDIM;dir++)
      dimensions[dir]=computeDimensions(dir);
    
    /// Communicators which make a given direction local
    std::vector<LocDirMaker> locDirMaker;
    for(Dir dir=0;dir<NDIM;dir++)
      locDirMaker.push_back(getLocDirMaker(dir));
    
    /// Inverse localizer
    std::vector<LocDirUnmaker> locDirUnmaker;
    for(Dir dir=0;dir<NDIM;dir++)
      locDirUnmaker.push_back(locDirMaker[dir()].inverse());
    
    for(Dir dir=0;dir<nDim-1;dir++)
      locDirChanger.push_back(locDirMaker[dir()]*locDirUnmaker[dir()+1]);
    
    firstLocDirMaker=new LocDirMaker(std::move(locDirMaker.back()));
    
    lastLocDirUnmaker=new LocDirUnmaker(std::move(locDirUnmaker.front()));
    
    initialized=true;
  }
  
  void dealloc()
  {
    if(initialized)
      {
	delete firstLocDirMaker;
	locDirChanger.clear();
	delete lastLocDirUnmaker;
      }
    else
      crash("localizer not initialized");
    
    initialized=false;
  }
}
