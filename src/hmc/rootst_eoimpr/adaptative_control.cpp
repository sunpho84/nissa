#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <sstream>

#include "../../base/global_variables.h"
#include "../../base/random.h"
#include "../../base/thread_macros.h"
#include "../../new_types/new_types_definitions.h"
#include "../../routines/ios.h"
#ifdef USE_THREADS
 #include "../../routines/thread.h"
#endif

//copy the key into hmc_evol_pars
void copy_key(hmc_evol_pars_t &out,const std::pair<int,int> &in)
{
  out.nmd_steps=in.first;
  out.ngauge_substeps=in.second;
}

//increase or decrase
void random_increase(int &x)
{
  if(x==1) x=2;
  else x+=2*(int)rnd_get_unif(&glb_rnd_gen,0,2)-1;
}

//check that it is last
template <typename Iter> Iter next(Iter iter) {return ++iter;}
template <typename Iter, typename Cont> bool is_last(Iter iter,const Cont& cont)
{return (iter!=cont.end())&&(next(iter)==cont.end());}

//compute derivative around x
void compute_second_derivative_around(std::pair<int,std::pair<double,double> > *der,std::map<int,std::pair<double,double> > &slice,int x)
{
  std::map<int,std::pair<double,double> >::iterator it=slice.find(x);
  if(it==slice.begin())
    {
      der[0]=*(it++);
      der[1]=*(it++);
      der[2]=*(it++);
    }
  else
    if(is_last(it,slice))
      {
	der[2]=*(it--);
	der[1]=*(it--);
	der[0]=*(it--);
      }
    else
      {
	der[1]=*(it--);
	der[0]=*(it++);
	it++;
	der[2]=*it;
      }
}

//move around
void move_around(std::map<int,std::pair<double,double> > slice,int &x,const char *tag)
{
  //header
  master_printf("Checking %s timeslice around %d\n",tag,x);
  for(std::map<int,std::pair<double,double> >::iterator it=slice.begin();it!=slice.end();it++)
    master_printf(" %d %lg %lg is on the same %s timelsice\n",it->first,it->second.first,it->second.second,tag);
  
  // take a random step in one direction
  if(slice.size()<=1)
    {
      master_printf(" No info on how to move on this slice, just making a random choise\n");
      random_increase(x);
    }
  else
    {
      if(slice.size()>=3) //compute symmetric derivative
	{
	  std::pair<int,std::pair<double,double> > der[3];
	  compute_second_derivative_around(der,slice,x);
	  
	  int x0=der[0].first,x1=der[1].first,x2=der[2].first;
	  double y0=der[0].second.first,y1=der[1].second.first,y2=der[2].second.first;
	  double a=((y0-y1)/(x0-x1)-(y1-y2)/(x1-x2))/(x0-x2);
	  double b=((y0-a*x0*x0)-(y1-a*x1*x1))/(x0-x1);
	  //double c=y0-a*x0*x0-b*x0;
	  double first_der=2*a*x0+b;
	  
	  master_printf(" Second order derivative: %d %d %d %lg %lg %lg %lg\n",x0,x1,x2,y0,y1,y2,a);
	  master_printf(" First order derivative: %lg\n",first_der);
	  
	  if(a>0) 
	    {
	      //we are far away from the maximum: just use first order derivative
	      master_printf(" Positive parabola, moving along first order derivative\n");
	      if(first_der>0) x++;
	      else x=std::max(1,x-1);
	    }
	  else
	    {
	      //locate vertex
	      double xv=-b/(2*a);
	      master_printf(" Vertex: %lg\n",xv);
	      if(xv>x) x=std::min((int)(xv+0.5),x+1);
	      else     x=std::max(1,std::max((int)(xv+0.5),x-1));
	    }
	}
      else //compute first order derivative 
	{
	  int x0=slice.begin()->first,x1=slice.rbegin()->first;
	  double y0=slice.begin()->second.first,y1=slice.rbegin()->second.first;
	  double der=(y0-y1)/(x0-x1);
	  master_printf(" First order derivative: (%lg-%lg)/(%d-%d)=%lg\n",y0,y1,x0,x1,der);
	  if(der>0) x++;
	  else x=std::max(1,x-1);
	}
    }
}

//choose the best parameters for next trajectory
void choose_hmc_traj_pars(hmc_evol_pars_t &in,adaptative_algorithm_pars_t &pars,int itraj)
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD)
    {
      master_printf("----------- Self adapting HMC algorithm -----------\n");
      hmc_traj_stat_t *stat=&(pars.stat);
      
      //check consistency
      for(hmc_traj_stat_t::iterator it=stat->begin();it!=stat->end();it++)
	  {
	    int n=it->second.n;
	    double sx=it->second.sx;
	    double s2x=it->second.s2x;
	    if(s2x*n<sx*sx)
	      {
		master_printf("Error found in entry %d %d: %d %lg %lg\n",it->first.first,it->first.second,n,sx,s2x);
		stat->erase(it->first);
	      }
	    else master_printf("Checking %d %d: %d %lg %lg\n",it->first.first,it->first.second,n,sx,s2x);
	  }
      
      int nentries=stat->size();
      copy_key(in,pars.current);
      
      if(itraj<0) master_printf("Not yet thermalized, skipping adaptation\n");
      if(itraj>=pars.use_for) master_printf("Thermalized traj %d passed adaptation time (%d),"
					    " fixing to best choice\n",itraj,pars.use_for);
      if(nentries==0 && itraj>=0 && itraj<pars.use_for)
	master_printf("No previous history recorderd, using default settings");
      /*
	1) find if there is at least one entry
	   otherwise, we cannot but use default settings
	2) compute averages and error for all the data, 
	   finding the entry with the smallest number of iters
	3) if we skipped beyond use_for interaval, fix to max
	4) if there is an entry with less than three iterations use it
	5) find out the largest entry, taking into account errors
	6) make a step in a random direction
      */
      
      // 1) if there is at least one entry we can start discussing it
      if(nentries>=1 && itraj>=0)
	{
	  master_printf("Adaptative hmc parameters, nentries: %d\n",nentries);
	  
	  // 2) compute averages and errors finding the entry with the smaller number of iterations
	  double ave[nentries],err[nentries];
	  double rate_max=0;
	  int nmd[nentries],ngss[nentries];
	  int niters_min=0,ientry=0;
	  hmc_traj_stat_t::iterator it_min,it_max;
	  for(hmc_traj_stat_t::iterator it=stat->begin();it!=stat->end();it++)
	    {
	      nmd[ientry]=it->first.first;
	      ngss[ientry]=it->first.second;
	      it->second.ave_err(ave[ientry],err[ientry]);
	      
	      //print the ave end err
	      master_printf(" %d %d: (%d) %lg %lg\n",nmd[ientry],ngss[ientry],it->second.n,ave[ientry],err[ientry]);
	      
	      //update lower entry
	      int niters=it->second.n;
	      if(niters<niters_min||ientry==0)
		{
		  niters_min=niters;
		  it_min=it;
		}
	      
	      //update max rate
	      if(it->second.n>=3 && (ave[ientry]>rate_max||ientry==0))
		{
		  rate_max=ave[ientry];
		  it_max=it;
		}
	      
	      ientry++;
	    }
	  
	  // 3) if we skipped beyond max larger adapted traj, fix to the one with the largest rate
	  if(itraj>=pars.use_for)
	    {
	      copy_key(in,it_max->first);
	      master_printf("Finished adaptation time, sticking to to largest rate\n");
	    }
	  else
	    {
	      // 4) if there is an entry with an almost undefined error, use it
	      if(niters_min<3)
		{
		  copy_key(in,it_min->first);
		  master_printf("Need to finish estimating one entry\n");
		}
	      else
		{
		  // 5) find the entry with the larger rate of production taking into account error
		  complex extracted_rate,zero={0,0};
		  double max_rate=0;
		  int ientry=0;
		  hmc_traj_stat_t::iterator it_max_rate;
		  for(hmc_traj_stat_t::iterator it=stat->begin();it!=stat->end();it++)
		    {
		      //extract a gaussian number
		      if(ientry%2==0) rnd_get_gauss_complex(extracted_rate,&glb_rnd_gen,zero,1);
		      double rate=extracted_rate[ientry%2]*err[ientry]+ave[ientry];
		      master_printf(" %d ientry: %lg\n",ientry,rate);
		      if(ientry==0||max_rate<rate)
			{
			  max_rate=rate;
			  it_max_rate=it;
			}
		      
		      ientry++;
		    }
		  
		  //copy max rate
		  copy_key(in,it_max_rate->first);
		  master_printf("max rate: %d %d\n",in.nmd_steps,in.ngauge_substeps);
		  
		  // 6) explore neighbourhood and update
		  std::map<int,std::pair<double,double> > nmd_slice,ngss_slice;
		  for(int ientry=0;ientry<nentries;ientry++)
		    {
		      if(nmd[ientry]==in.nmd_steps) ngss_slice[ngss[ientry]]=std::make_pair(ave[ientry],err[ientry]);
		      if(ngss[ientry]==in.ngauge_substeps)
			nmd_slice[nmd[ientry]]=std::make_pair(ave[ientry],err[ientry]);
		    }        
		  move_around(nmd_slice,in.nmd_steps,"nmd");
		  move_around(ngss_slice,in.ngauge_substeps,"ngss");
		}
	    }
	}
      
      //copy back
      pars.current.first=in.nmd_steps;
      pars.current.second=in.ngauge_substeps;
      
      //write final output
      master_printf("Using: %d %d\n\n",in.nmd_steps,in.ngauge_substeps);
    }
  
  THREAD_BARRIER();
}

//read the cureent pars
void adaptative_algorithm_pars_t::init_from_text(char *text)
{
  //init with current pars
  std::istringstream is(text);
  is>>current.first>>current.second;
  
  //read numer of entries
  int nentries;
  is>>nentries;
  
  //read all entries
  for(int ientry=0;ientry<nentries;ientry++)
    {
      int nmd,ngss,n;
      double sx,s2x;
      is>>nmd>>ngss>>n>>sx>>s2x;
      set(nmd,ngss,n,sx,s2x);
    }
}

//print the current pars
std::string adaptative_algorithm_pars_t::save_to_text()
{
  //init with current pars
  std::ostringstream os;
  os<<current.first<<" "<<current.second<<" ";
  
  //print number of entries
  int nentries=stat.size();
  os<<nentries<<" ";
  
  //save all the stat
  for(hmc_traj_stat_t::iterator it=stat.begin();it!=stat.end();it++)
    os<<it->first.first<<" "<<it->first.second<<" "<<it->second.n<<" "<<it->second.sx<<" "<<it->second.s2x<<" ";
      
  return os.str();
}
