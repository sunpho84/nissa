#include <nissa.hpp>

using namespace nissa;

int T,L;
vector_remap_t *remapper;

void index_from_lx_to_Neo(Rank &rank_out,LocLxSite &iel_out,const LocLxSite& iel_in,void *pars)
{
  Coords<Dir> mu_ord;
  mu_ord(tDir)=0;
  mu_ord(xDir)=zDir;
  mu_ord(yDir)=yDir;
  mu_ord(zDir)=xDir;
  const LocLxSite ilx=iel_in/NDIM;
  const Dir mu=(int)(iel_in()%NDIM);
  
  //odd sites goes with themseleves
  const bool shift_comp=(loclx_parity(ilx)==0);
  
  GlbCoords g;
  FOR_ALL_DIRS(nu)
    g(mu_ord(nu))=glbCoordOfLoclx(ilx,nu);
  if(shift_comp)
    g(mu_ord(mu))=(g(mu_ord(mu))+1)%glbSize(mu_ord(mu));
  
  int glb_site_dest=(glblx_of_coord(g)()/2)*8+mu_ord(mu)()*2+shift_comp;
  rank_out=glb_site_dest/(8*locVolh());
  iel_out=glb_site_dest-rank_out()*8*locVolh();
}

void conf_convert(char *outpath,char *inpath)
{
  quad_su3 *conf=nissa_malloc("conf",locVolWithBord.nastyConvert(),quad_su3);
  read_ildg_gauge_conf(conf,inpath);
  
  //compute and convert the plaquette
  double plaq=global_plaquette_lx_conf(conf)*3;
  if(not little_endian) change_endianness(&plaq,&plaq,1);
  
  //convert the lattice size
  GlbCoords temp;
  FOR_ALL_DIRS(mu)
    if(not little_endian)
      change_endianness(temp(mu)(),glbSize(mu)());
    else
      temp(mu)=glbSize(mu);
  
  //write the header
  FILE *fout=fopen(outpath,"w");
  if(rank==0)
    {
      int nw;
      FOR_ALL_DIRS(mu)
	{
	  const int c=temp(mu)();
	  nw=fwrite(&c,sizeof(int),1,fout);
	  if(nw!=1) crash("did not success in writing");
	}
      nw=fwrite(&plaq,sizeof(double),1,fout);
      if(nw!=1) crash("did not success in writing");
    }
  MPI_Barrier(MPI_COMM_WORLD);
  
  //if needed convert the endianess of the conf
  if(!little_endian) change_endianness((double*)conf,(double*)conf,locVol()*4*18);
  
  //reorder
  char *buf=nissa_malloc("buf",locVol.nastyConvert()*sizeof(quad_su3),char);
  remapper->remap(buf,conf,sizeof(su3));
  
  //seek to the correct point
  fseek(fout,24+rank*sizeof(quad_su3)*locVol(),SEEK_SET);
  MPI_Barrier(MPI_COMM_WORLD);
  
  //write
  int nw=fwrite(buf,sizeof(quad_su3),locVol(),fout);
  if(nw!=locVol) crash("did not success in writing");
  MPI_Barrier(MPI_COMM_WORLD);
  
  //flush
  fflush(fout);
  MPI_Barrier(MPI_COMM_WORLD);
  
  //close
  fclose(fout);
  MPI_Barrier(MPI_COMM_WORLD);
  
  nissa_free(buf);
  nissa_free(conf);
}

void in_main(int narg,char **arg)
{
  if(narg<2) crash("Use: %s input",arg[0]);
  
  open_input(arg[1]);
  
  // 1) Read information about the gauge conf
  
  //Read the volume
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  //Init the MPI grid
  init_grid(T,L);
  
  //init the remapper
  remapper=new vector_remap_t(4*locVol(),index_from_lx_to_Neo,NULL);
  
  //read the number of gauge configurations
  int N;
  read_str_int("NGaugeConf",&N);
  
  ///////////////////////////////////////////
  
  for(int i=0;i<N;i++)
  {
      char in[1024],out[1024];
      read_str(in,1024);
      read_str(out,1024);
      
      if(not file_exists(out))
	 conf_convert(out,in);
      else
	master_printf("Skipping conf %s\n",in);
  }
  
  close_input();
  
  ///////////////////////////////////////////
  
  delete remapper;
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
