#include <nissa.hpp>
#define N 21

using namespace nissa;

//ON type
typedef complex ON_t[N];
comm_t lx_ON_t_comm,eo_ON_t_comm;
comm_t lx_quad_u1_comm,eo_quad_u1_comm;
DEFINE_BORDERS_ROUTINES(ON_t)
DEFINE_BORDERS_ROUTINES(quad_u1)

//configuration
ON_t *zeta,*zeta_old,*pi,*fpi;
quad_u1 *lambda,*lambda_old;
momentum_t *omega,*fomega;

bool conf_created=true;
double beta,g;
int nsweep,nterm;
int nhmc_steps,nacc;

//write a conf adding info
void write_conf(const char *path,int iconf)
{
  //open for writing
  ILDG_File file=ILDG_File_open_for_write(path);
  
  //write data
  char data_iconf[30],data_rnd[1000];
  sprintf(data_iconf,"%d",iconf);
  convert_rnd_gen_to_text(data_rnd,&glb_rnd_gen);
  ILDG_File_write_record(file,"ConfID",data_iconf,strlen(data_iconf));
  ILDG_File_write_record(file,"RND_gen_status",data_rnd,strlen(data_rnd));
  write_double_vector(file,(double*)lambda,sizeof(quad_u1)/sizeof(double),64,"lambda");
  write_double_vector(file,(double*)zeta,sizeof(ON_t)/sizeof(double),64,"zeta");
  
  ILDG_File_close(file);
}

//load the conf
void read_conf(const char *path,int &iconf)
{
  ILDG_File file=ILDG_File_open_for_read(path);
  
  //read
  ILDG_header header;
  
  //iconf
  header=ILDG_File_get_next_record_header(file);
  if(strcasecmp(header.type,"ConfID")) CRASH("expecting ConfID record, found %s",header.type);
  char data_iconf[30];
  ILDG_File_read_all(data_iconf,file,header.data_length);
  sscanf(data_iconf,"%d",&iconf);
  
  //rnd
  header=ILDG_File_get_next_record_header(file);
  if(strcasecmp(header.type,"RND_gen_status")) CRASH("expecting RND_gen_status record, found %s",header.type);
  char data_rnd[1000];
  ILDG_File_read_all(data_rnd,file,header.data_length);
  start_loc_rnd_gen(data_rnd);
  
  //lambda
  header=ILDG_File_get_next_record_header(file);
  if(strcasecmp(header.type,"lambda")) CRASH("expecting lambda record, found %s",header.type);
  read_real_vector((double*)lambda,file,header,sizeof(quad_u1)/sizeof(double));

  //zeta
  header=ILDG_File_get_next_record_header(file);
  if(strcasecmp(header.type,"zeta")) CRASH("expecting zeta record, found %s",header.type);
  read_real_vector((double*)zeta,file,header,sizeof(ON_t)/sizeof(double));
  
  ILDG_File_close(file);
}


//initialize cpn simulation
void read_input(int &seed,const char *path)
{
  open_input(path);
  
  //read geometry
  int L;
  read_str_int("L",&L);
  init_grid(L,L);
  
  //read beta and compute g
  read_str_double("Beta",&beta);
  g=1/(N*beta);

  //read seed and initialize generator
  read_str_int("Seed",&seed);
  
  //read the number of sweeps and substeps
  read_str_int("NSweep",&nsweep);
  read_str_int("NTerm",&nterm);
  read_str_int("NHmcSteps",&nhmc_steps);

  close_input();
}

//initialize the system to id
void init_system_to_cold()
{
  
  //set the configuration to cold
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      memset(zeta[ivol],0,sizeof(ON_t));
      complex_put_to_real(zeta[ivol][0],1);
      for(int mu=0;mu<NDIM;mu++) complex_put_to_real(lambda[ivol][mu],1);
    }
  
  set_borders_invalid(zeta);
  set_borders_invalid(lambda);
}

//initialize the code
void init(int seed,int &base_isweep)
{
  set_lx_comm(lx_ON_t_comm,sizeof(ON_t));
  set_lx_comm(lx_quad_u1_comm,sizeof(quad_u1)); 
  
  //allocate configuration
  zeta=nissa_malloc("zeta",loc_vol+bord_vol,ON_t);
  lambda=nissa_malloc("lambda",loc_vol+bord_vol,quad_u1);
  
  //allocate configuration
  zeta_old=nissa_malloc("zeta_old",loc_vol+bord_vol,ON_t);
  lambda_old=nissa_malloc("lambda_old",loc_vol+bord_vol,quad_u1);
  
  //allocate momenta
  pi=nissa_malloc("pi",loc_vol,ON_t);
  omega=nissa_malloc("omega",loc_vol,momentum_t);
  
  //allocate force
  fpi=nissa_malloc("fpi",loc_vol,ON_t);
  fomega=nissa_malloc("fomega",loc_vol,momentum_t);

  init_system_to_cold();
  nacc=0;
  base_isweep=0;
  
  if(file_exists("conf"))
    {
      read_conf("conf",base_isweep);
      conf_created=false;
    }
  else
    {
      start_loc_rnd_gen(seed);
      conf_created=true;
    }
}

//close everything
void close_cpn(int isweep)
{
  write_conf("conf",isweep);
  
  MASTER_PRINTF("Closing CPN\n");
  
  nissa_free(zeta);
  nissa_free(lambda);
  nissa_free(zeta_old);
  nissa_free(lambda_old);
  nissa_free(pi);
  nissa_free(omega);
  nissa_free(fpi);
  nissa_free(fomega);
}

//compute the total energy or action
void energy(double* en,ON_t* z,quad_u1* l)
{

  communicate_lx_ON_t_borders(z);

  double loc_thread_res=0;
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<NDIM;mu++)
      {
        int iup=loclx_neighup[ivol][mu];
        for(int n=0;n<N;n++)
	  {
	    complex pr;
	    unsafe_complex_conj1_prod(pr,z[iup][n],z[ivol][n]);
	    loc_thread_res+=real_part_of_complex_prod(pr,l[ivol][mu]);
	  }
      }

  (*en)=-(2*glb_reduce_double(loc_thread_res)-2*glb_vol*NDIM);
}
double energy(ON_t *z,quad_u1 *l)
{
  double en;
  energy(&en,z,l);
  return en;
}
double action(ON_t *z,quad_u1 *l)
{return energy(z,l)/g;}

//return a scalar product between two zetas
void get_zeta_complex_scalar_prod(complex res,ON_t a,ON_t b)
{
  res[RE]=res[IM]=0;
  for(int n=0;n<N;n++) complex_summ_the_conj1_prod(res,a[n],b[n]);
}

//only real part
double get_zeta_real_scalar_prod(const ON_t a,const ON_t b)
{
  double res=0;
  for(int n=0;n<N;n++) res+=real_part_of_complex_scalar_prod(a[n],b[n]);
  return res;
}

//return the norm squared of a zeta
inline double get_zeta_norm2(const ON_t z)
{return get_zeta_real_scalar_prod(z,z);}
inline double get_zeta_norm(const ON_t z)
{return sqrt(get_zeta_norm2(z));}

//orthogonalize
void zeta_orthogonalize_with(ON_t z,const ON_t w)
{
  double norm_with=get_zeta_real_scalar_prod(w,z)/get_zeta_norm2(w);
  for(int n=0;n<N;n++) complex_subt_the_prod_double(z[n],w[n],norm_with);
}

//reunitarize a zeta
void zeta_unitarize(ON_t z)
{
  double zeta_norm_inv=1/get_zeta_norm(z);
  for(int n=0;n<N;n++) complex_prodassign_double(z[n],zeta_norm_inv);
}

//reunitarize a lambda
double lambda_unitarize(complex l)
{
  double n=sqrt(squared_complex_norm(l));
  complex_prodassign_double(l,1/n);
  return n;
}

//draw momenta
void generate_momenta()
{
  
  complex ave={0,0};
  const double sigma=sqrt(2);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      //generate the lambda momenta, gaussianly
      for(int mu=0;mu<NDIM;mu++) omega[ivol][mu]=rnd_get_gauss_double(loc_rnd_gen+ivol);
      
      //generate zeta momenta and orthogonalize them
      for(int n=0;n<N;n++) rnd_get_gauss_complex(pi[ivol][n],loc_rnd_gen+ivol,ave,sigma);
      //orthogonalize
      zeta_orthogonalize_with(pi[ivol],zeta[ivol]);
    }
  THREAD_BARRIER();
}

//compute the action for zeta momenta
void zeta_momenta_action(double* act)
{
  
  double loc_thread_act=0;
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    loc_thread_act+=get_zeta_real_scalar_prod(pi[ivol],pi[ivol]);
  loc_thread_act/=2;
  (*act)=glb_reduce_double(loc_thread_act);
}

//compute the action for lambda momenta
void lambda_momenta_action(double* act)
{
  
  double loc_thread_act=0;
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<NDIM;mu++)
      loc_thread_act+=sqr(omega[ivol][mu]);
  loc_thread_act/=2;
  (*act)=glb_reduce_double(loc_thread_act);
}

//compute the action of momenta
double momenta_action()
{
  double zm,lm;
  zeta_momenta_action(&zm);
  lambda_momenta_action(&lm);
  
  return zm+lm;
}

//compute lambda force
void compute_lambda_forces()
{
  
  communicate_lx_ON_t_borders(zeta);
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<NDIM;mu++)
      {
        //reset
        fomega[ivol][mu]=0;
	int iup=loclx_neighup[ivol][mu];
        for(int n=0;n<N;n++)
	  fomega[ivol][mu]+=(-lambda[ivol][mu][RE]*zeta[iup][n][IM]+lambda[ivol][mu][IM]*zeta[iup][n][RE])*zeta[ivol][n][RE]+
	    (+lambda[ivol][mu][RE]*zeta[iup][n][RE]+lambda[ivol][mu][IM]*zeta[iup][n][IM])*zeta[ivol][n][IM];
        fomega[ivol][mu]*=-2*beta*N;
      }
  
#if 0
  int site=0;
  double eps=1.e-6;
  double pre_act=action(zeta,lambda);
  for(int mu=0;mu<NDIM;mu++)
    {
      double pre_val=atan2(lambda[site][mu][IM],lambda[site][mu][RE]);
      
      lambda[site][mu][RE]=cos(pre_val+eps);
      lambda[site][mu][IM]=sin(pre_val+eps);
      
      double post_act=action(zeta,lambda);
      lambda[site][mu][RE]=cos(pre_val);
      lambda[site][mu][IM]=sin(pre_val);
      
      double f=-(post_act-pre_act)/eps;
      MASTER_PRINTF("mu: %d f: %lg %lg\n",mu,f,f/fomega[site][mu]);
    }
#endif
}

//compute the staple of zeta
inline void site_staple(ON_t staple,ON_t *z,quad_u1 *l,int ivol)
{
  for(int n=0;n<N;n++) staple[n][RE]=staple[n][IM]=0;

  for(int mu=0;mu<NDIM;mu++)
    {
      int iup=loclx_neighup[ivol][mu];
      for(int n=0;n<N;n++) complex_summ_the_conj2_prod(staple[n],z[iup][n],l[ivol][mu]);
      int idw=loclx_neighdw[ivol][mu];
      for(int n=0;n<N;n++) complex_summ_the_prod(staple[n],z[idw][n],l[idw][mu]);
    }
  for(int n=0;n<N;n++) complex_prodassign_double(staple[n],2);
}

//compute zeta forces
void compute_zeta_forces()
{
  
  communicate_lx_ON_t_borders(zeta);
  communicate_lx_quad_u1_borders(lambda);
  
  //zeta orthogonalized spin projection
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      site_staple(fpi[ivol],zeta,lambda,ivol);
      zeta_orthogonalize_with(fpi[ivol],zeta[ivol]);
      for(int n=0;n<N;n++) complex_prodassign_double(fpi[ivol][n],beta*N);
    }

#if 0
  int site=0;
  double eps=1.e-8;
  double pre_act=action(zeta,lambda);
  ON_t pre_val;
  for(int m=0;m<N;m++) complex_copy(pre_val[m],zeta[site][m]);
  for(int n=0;n<N;n++)
    {
      zeta[site][n][IM]=pre_val[n][IM]+eps;
      zeta_unitarize(zeta[site]);
      double post_act=action(zeta,lambda);
      
      double f=-(post_act-pre_act)/eps;
      
      MASTER_PRINTF("n: %d f: %lg %lg\n",n,f,f/fpi[site][n][IM]);
      for(int m=0;m<N;m++) complex_copy(zeta[site][m],pre_val[m]);
    }
#endif
}

//update pi momenta
inline void update_zeta_momenta(double eps)
{
  
  //compute the zeta forces
  compute_zeta_forces();

  //update zeta momenta
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol) for(int n=0;n<N;n++) complex_summ_the_prod_double(pi[ivol][n],fpi[ivol][n],eps);
  THREAD_BARRIER();
}

//update omega momenta
inline void update_lambda_momenta(double eps)
{
  
  //compute the lambda forces
  compute_lambda_forces();

  //update momenta
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol) for(int mu=0;mu<NDIM;mu++) omega[ivol][mu]+=fomega[ivol][mu]*eps;
  THREAD_BARRIER();
}

//update both momenta
inline void update_momenta(double eps)
{
  update_zeta_momenta(eps);
  update_lambda_momenta(eps);
}

//update the zetas
inline void update_zeta_positions(double eps)
{
  
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      //compute the norm of pi, to form the matrix
      double pi_norm=get_zeta_norm(pi[ivol]);
      
      //compute parameters of rotating matrix
      double al=eps*pi_norm;
      double cal=cos(al),sal=sin(al);
  
      //update zeta according to ortho-bound
      for(int n=0;n<N;n++)
        {
          //get old values of coord and momenta for z
          complex x,p;
	  complex_copy(x,zeta[ivol][n]);
          complex_copy(p,pi[ivol][n]);
          
          //rotate
	  complex_prodassign_double(zeta[ivol][n],cal);
	  complex_summ_the_prod_double(zeta[ivol][n],p,sal/pi_norm);
          complex_prodassign_double(pi[ivol][n],cal);
	  complex_summ_the_prod_double(pi[ivol][n],x,-pi_norm*sal);
        }
    }
  set_borders_invalid(zeta);
}

//update the lambdas
inline void update_lambda_positions(double eps)
{
  
  //update lambda
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    for(int mu=0;mu<NDIM;mu++)
      {
	complex t={cos(eps*omega[ivol][mu]),sin(eps*omega[ivol][mu])};
	safe_complex_prod(lambda[ivol][mu],lambda[ivol][mu],t);
      }
  set_borders_invalid(lambda);
}

//update zetas and lambdas
void update_positions(double eps)
{
  update_zeta_positions(eps);
  update_lambda_positions(eps);
}

//integrate equation of motion
void hmc_integrate(double tl)
{
  
  double dt=tl/nhmc_steps/2,dth=dt/2,ldt=dt*OMELYAN_LAMBDA,l2dt=2*OMELYAN_LAMBDA*dt,uml2dt=(1-2*OMELYAN_LAMBDA)*dt;

  //     Compute H(t+lambda*dt) i.e. v1=v(t)+a[r(t)]*lambda*dt (first half step)
  update_momenta(ldt);
  
  //         Main loop
  for(int istep=0;istep<nhmc_steps;istep++)
    {
      //decide if last step is final or not
      double last_dt=(istep==(nhmc_steps-1)) ? ldt : l2dt;
        
      //     Compute U(t+dt/2) i.e. r1=r(t)+v1*dt/2
      update_positions(dth);
      //     Compute H(t+(1-2*lambda)*dt) i.e. v2=v1+a[r1]*(1-2*lambda)*dt
      update_momenta(uml2dt);
      //     Compute U(t+dt/2) i.e. r(t+dt)=r1+v2*dt/2
      update_positions(dth);
      //     Compute H(t+dt) i.e. v(t+dt)=v2+a[r(t+dt)]*lambda*dt (at last step) or *2*lambda*dt
      update_momenta(last_dt);
    }
  
  //check_lambda_conf_unitarity(lambda);
  //check_zeta_conf_unitarity(zeta);
  //normalize the configuration
  NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
    {
      zeta_unitarize(zeta[ivol]);
      for(int mu=0;mu<NDIM;mu++)
        lambda_unitarize(lambda[ivol][mu]);
    }
  set_borders_invalid(zeta);
  set_borders_invalid(lambda);
}

//perform a hybid monte carlo update
void hmc_update(bool skip_test=false)
{
  //copy configuration
  vector_copy(zeta_old,zeta);
  vector_copy(lambda_old,lambda);
  
  //generate momenta and compute action
  generate_momenta();
  double start_mom_action=momenta_action();
  double start_theo_action=action(zeta,lambda);
  double start_action=start_mom_action+start_theo_action;
  VERBOSITY_LV2_MASTER_PRINTF(" action: mom=%lg, coord=%lg\n",start_mom_action,start_theo_action);
  
  //integrate for unitary length
  hmc_integrate(1.0);
  
  //compute final action
  double final_mom_action=momenta_action();
  double final_theo_action=action(zeta,lambda);
  double final_action=final_mom_action+final_theo_action;
  MASTER_PRINTF(" action: mom=%lg, coord=%lg\n",final_mom_action,final_theo_action);

  //compute difference of action and print it
  double diff_action=final_action-start_action;
  
  //make metropolis test
  bool acc=metro_test(diff_action);
  
  //copy back old configuration
  if(!skip_test && !acc)
    {
      vector_copy(zeta,zeta_old);
      vector_copy(lambda,lambda_old);
    }
  
  const char acc_flag[2][4]={"rej","acc"};
  MASTER_PRINTF("diff_action: %lg-%lg=%lg, %s %s\n",final_action,start_action,diff_action,acc_flag[acc],(skip_test?" (skip)":""));
  
  if(acc) nacc++;
}

//return the geometric definition of topology
void geometric_topology_simplified(double* topo,ON_t* z)
{
  
  communicate_lx_ON_t_borders(z);
  
  double loc_topo_thread=0;
  int mu=0,nu=1;
  NISSA_PARALLEL_LOOP(n,0,loc_vol)
    {
      int nPmu=loclx_neighup[n][mu];
      int nPnu=loclx_neighup[n][nu];
      int nMmu=loclx_neighdw[n][mu];
      int nMnu=loclx_neighdw[n][nu];
      
      complex a,b,c,d,e,f;
      get_zeta_complex_scalar_prod(a,z[nPnu],z[nMmu]);
      get_zeta_complex_scalar_prod(b,z[n],z[nPnu]);
      get_zeta_complex_scalar_prod(c,z[nMmu],z[n]);
      get_zeta_complex_scalar_prod(d,z[n],z[nMnu]);
      get_zeta_complex_scalar_prod(e,z[nPmu],z[n]);
      get_zeta_complex_scalar_prod(f,z[nMnu],z[nPmu]);
      complex ab;
      unsafe_complex_prod(ab,a,b);
      complex abc;
      unsafe_complex_prod(abc,ab,c);
      complex de;
      unsafe_complex_prod(de,d,e);
      complex def;
      unsafe_complex_prod(def,de,f);
      loc_topo_thread+=
        atan2(abc[IM],abc[RE])+
	atan2(def[IM],def[RE]);
    }
  (*topo)=glb_reduce_double(loc_topo_thread)/(2*M_PI);
}
double geometric_topology_simplified(ON_t *z)
{double topo;geometric_topology_simplified(&topo,z);return topo;}

void in_main(int narg,char **arg)
{
  //init
  if(narg<2) CRASH("Use %s input",arg[0]);
  int seed;
  read_input(seed,arg[1]);
  int base_isweep;
  init(seed,base_isweep);
  
  //open observable files
  FILE *energy_file=open_file("energy_file",conf_created?"w":"a");
  FILE *topology_file=open_file("topology",conf_created?"w":"a");
  
  int isweep;
  for(isweep=base_isweep;isweep<base_isweep+nsweep;isweep++)
    {
      //compute energy
      master_fprintf(energy_file,"%d %12.12lg\n",isweep,energy(zeta,lambda)/glb_vol/NDIM);
      master_fprintf(topology_file,"%d %12.12lg\n",isweep,geometric_topology_simplified(zeta));

      hmc_update(isweep<nterm);
    }
  
  //close observables file
  close_file(energy_file);
  close_file(topology_file);
  
  //close
  close_cpn(isweep);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  close_nissa();
    
  return 0;
}
