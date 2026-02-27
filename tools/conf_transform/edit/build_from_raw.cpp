#include <nissa.hpp>
#include <inttypes.h>

using namespace nissa;

//open the file and scan it until an su3 link is found
FILE* find_conf_beginning(std::string path)
{
  FILE *fin=open_file(path,"r");
  su3 link;
  
  long long int off=0;
  int rc=0,link_found=0;
  double non_su3=0;
  do
    {
      VERBOSITY_LV3_MASTER_PRINTF("==============================\n");
      //seek
      rc=fseeko64(fin,off,SEEK_SET);
      
      if(rc) MASTER_PRINTF("Error seeking to %lld\n",off);
      else
	{
	  //check feof
	  if(feof(fin))
	    {
	      CRASH("Reached EOF");
	      rc=1;
	    }
	  //read
	  int nr=fread(link,sizeof(su3),1,fin);
      CRASH("reimplement");
	  // if(little_endian) change_endianness((double*)link,(double*)link,sizeof(su3)/sizeof(double));
	  
	  if(nr!=1)
	    {
	      CRASH("Error reading: %d",nr);
	      rc=1;
	    }
	  else
	    {
	      non_su3=su3_get_non_unitariness(link);
	      if(fabs(non_su3)<1e-13)
		{
		  link_found=1;
		  VERBOSITY_LV2_MASTER_PRINTF("Link found at %lld, non-su3ness: %lg\n",off,non_su3);
		}
	      else
		{
		  VERBOSITY_LV3_MASTER_PRINTF("Link not yet found, deviation: %lg at %lld\n",non_su3,off);
		  off++;
		}
	    }
	  // double d=link[NCOL-1][NCOL-1][IM];
	  // if(fabs(d)<1 and fabs(d)>1e-6) VERBOSITY_LV2_MASTER_PRINTF("Suspect!\n");
	}
    }
  while(rc==0 and (not link_found));
  
  //set to the position
  if(link_found)
    {
      VERBOSITY_LV2_MASTER_PRINTF("Putting offset to beginning of found link: %lld\n",off);
      rc=fseeko64(fin,off,SEEK_SET);
    }
  
  return fin;
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  initNissa(narg,arg);
  
  if(narg<5) CRASH("use: %s L T file_in file_out",arg[0]);
  
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  
  //Init the MPI grid
  initGrid(T,L);
  //////////////////////////// read the conf /////////////////////////////
  
      CRASH("reimplement");
  // quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol,quad_su3);
  
  // FILE *fin=find_conf_beginning(arg[3]);
  // int rc=fread(conf,sizeof(quad_su3),glbVol,fin);
  // if(rc!=glbVol) CRASH("Unable to read, returned: %d",rc);
  // close_file(fin);
  
  // //convert and reorder
  // if(little_endian) change_endianness((double*)conf,(double*)conf,glbVol*sizeof(quad_su3)/sizeof(double));
  // vector_remap_t(locVol,index_from_ILDG_remapping).remap(conf,conf,sizeof(quad_su3));
  // NISSA_PARALLEL_LOOP(ivol,0,locVol)
  //   quad_su3_ildg_to_nissa_reord(conf[ivol],conf[ivol]);
  // NISSA_PARALLEL_LOOP_END;
  
  // //perform unitarity test
  // unitarity_check_result_t unitarity_check_result;
  // unitarity_check_lx_conf(unitarity_check_result,conf);
  
      CRASH("reimplement");
  // VERBOSITY_LV1_MASTER_PRINTF("Plaquette of read conf: %16.16lg\n",global_plaquette_lx_conf(conf));
  // VERBOSITY_LV1_MASTER_PRINTF("Deviation from unitarity: %lg average, %lg max\n",unitarity_check_result.average_diff,unitarity_check_result.max_diff);
  
  //////////////////////////// write the conf ////////////////////////////
  
      CRASH("reimplement");
  // write_ildg_gauge_conf(arg[4],conf,64);
  // nissa_free(conf);
  
  ///////////////////////////////////////////
  
  closeNissa();
  
  return 0;
}
