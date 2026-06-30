#ifndef _CONF_HPP
#define _CONF_HPP

#include <base/field.hpp>
#include <geometry/geometry_lx.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  inline int nAnalyzedConfs;
  inline int nTotHitsDone;
  inline double tot_prog_time,wall_time;
  
  inline double conf_load_time;
  inline int nconf_load;
  
  inline double ape_time;
  inline int nape_tot;
  
  inline char conf_path[1024],outfolder[1024];
  inline int ngauge_conf;
  inline int inner_conf_valid;
  inline bool conf_allocated{false};
  inline LxField<quad_su3>* glb_conf;
  inline LxField<quad_su3>* inner_conf;
  inline LxField<quad_su3> *ape_smeared_conf;
  
  inline int lock_fd;
  
  void allocate_confs();
  void free_confs();
  void read_init_grid();
  Coords generate_random_coord();
  
  LxField<quad_su3>* get_updated_conf(const double& charge,
				      const Momentum& theta,
				      const LxField<quad_su3>& in_conf);
  
  void start_new_conf();
  
  void setup_conf(LxField<quad_su3>& conf,
		  const char *conf_path,
		  const int& rnd_gauge_transform,
		  const int& free_theory);
  
  /// Check if the time is enough
  inline bool isRemainingTimeEnough()
  {
    if(nTotHitsDone)
      {
	MASTER_PRINTF("Total number of hits done: %d\n",nTotHitsDone);
	
	const double passedTime=
	  take_time()+tot_prog_time;
	MASTER_PRINTF("Total elapsed time: %lg s\n",passedTime);
	
	const double aveTimePerHit=
	  passedTime/nTotHitsDone;
	MASTER_PRINTF("Average time per hit: %lg sec, pessimistically: %lg\n",aveTimePerHit,aveTimePerHit*1.1);
	
	const double remainingTime=
	  wall_time-passedTime;
	MASTER_PRINTF("Remaining time: %lg sec\n",remainingTime);
	
	const bool isRemainingTimeEnough=
	  broadcast(remainingTime>(aveTimePerHit*1.1));
	
	if(isRemainingTimeEnough)
	  MASTER_PRINTF("Time is enough to go on!\n");
	else
	  MASTER_PRINTF("Not enough time, exiting!\n");
	
	return isRemainingTimeEnough;
      }
    else
      return true;
  }
  
  int read_conf_parameters(int &iconf);
  
  /// Forward declaration of the hit looper
  struct HitLooper;
  
  void releaseConf();
  
  /// Cleanup the conf
  void finalizeConf(const HitLooper& hitLooper);
  
  void print_statistics();
  void skip_nhits(int a,int b);
}

#endif
