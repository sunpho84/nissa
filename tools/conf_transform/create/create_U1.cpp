#include <nissa.hpp>

using namespace nissa;

typedef double realspin1field[NDIM];

typedef realspin1field realspin1prop[NDIM];

void write_u1_field(ILDG_File file,
		    const LxField<realspin1field>& field)
{
  LxField<realspin1field> permuted("permuted");
  
  PAR(0,
      locVol,
      CAPTURE(TO_READ(field),
	      TO_WRITE(permuted)),
      loclx,
      {
	for(int mu=0;mu<NDIM;mu++)
	  permuted[loclx][(mu+NDIM-1)%NDIM]=field[loclx][mu];
      });
  
  write_real_vector(file,permuted,"ildg-binary-data");
}

int main(int narg,char **arg)
{
  initNissa(narg,arg);
  
  //check argument
  if(narg<7) CRASH("Use: %s L T gauge[L=LANDAU,F=FEYNMAN,C=COULOMB] seed nconfs pattern [e for plaquette] [su3conf input] [u3conf output]",arg[0]);
  
  const int L=atoi(arg[1]);
  const int T=atoi(arg[2]);
  const char gauge=arg[3][0];
  const int seed=atoi(arg[4]);
  const int nConfs=atoi(arg[5]);
  const char* pattern=arg[6];
  const double e=(narg>7)?strtod(arg[7],nullptr):0;
  const char* su3Path=(narg>8)?arg[8]:nullptr;
  const char* u3Path=(narg>9)?arg[9]:nullptr;
  
  initGrid(T,L);
  
  field_rng_stream.init(seed);
  
  gauge_info gl;
  gl.zms=UNNO_ALEMANNA;
  switch(gauge)
    {
    case 'L':
      gl.which_gauge=gauge_info::LANDAU;
      break;
    case 'F':
      gl.which_gauge=gauge_info::FEYNMAN;
      break;
    case 'C':
      gl.which_gauge=gauge_info::COULOMB;
      break;
    default:
      CRASH("unknown gauge '%c'",gauge);
    }
  gl.c1=WILSON_C1;
  
  LxField<spin1field> photonEta("photonEta");
  LxField<spin1field> photonField("photonField");
  LxField<realspin1field> realPhotonField("realPhotonField");
  // spin1prop *propRecoMom=nissa_malloc("propRecoMom",locVol,spin1prop);
  // spin1prop *propRecoMom2=nissa_malloc("propRecoMom2",locVol,spin1prop);
  // spin1prop *propReco=nissa_malloc("propReco",locVol,spin1prop);
  // vector_reset(propRecoMom);
  // vector_reset(propRecoMom2);
  
  for(int iConf=0;iConf<nConfs;iConf++)
    {
      auto sourceFiller=
	field_rng_stream.getDrawer<spin1field>();
      
      sourceFiller.fillField(photonEta);
      
      PAR(0,
	  locVol,
	  CAPTURE(TO_WRITE(photonEta)),
	  loclx,
	  {
	    for(int mu=0;mu<NDIM;mu++)
	      z2Transform(photonEta[loclx][mu]);
	  });
      
      multiply_by_sqrt_tlSym_gauge_propagator(photonField,photonEta,gl);
      
      PAR(0,
	  locVol,
	  CAPTURE(TO_WRITE(realPhotonField),
		  TO_READ(photonField)),
	  loclx,
	  {
	    for(int mu=0;mu<NDIM;mu++)
	      realPhotonField[loclx][mu]=photonField[loclx][mu][RE];
	  });
      
      char outPath[128];
      snprintf(outPath,128,pattern,iConf);
      if(strncasecmp(outPath,pattern,128)==0)
	CRASH("pattern %s cannot be used to creat the conf name, try something like: conf%%04d.dat",pattern);
      
      ILDG_File file=ILDG_File_open_for_write(outPath);
      char ildg_format_message[1024];
      snprintf(ildg_format_message,1024,
	      "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	      "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n"
	      "            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
	      "            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n"
	      "  <version>1.0</version>\n"
	      "  <field>photonField</field>\n"
	      "  <precision>%zu</precision>\n"
	      "  <dof>4</dof>\n"
	      "  <lx>%d</lx>\n"
	      "  <ly>%d</ly>\n"
	      "  <lz>%d</lz>\n"
	      "  <lt>%d</lt>\n"
	      "</ildgFormat>",
	      64lu,glbSize[3],glbSize[2],glbSize[1],glbSize[0]);
      ILDG_File_write_text_record(file,"ildg-format",ildg_format_message);
      
      write_u1_field(file,realPhotonField);
      ILDG_File_close(file);
      
      // pass_spin1field_from_x_to_mom_space(photonEta,photonField,gl.bc,true,true);
      
      // NISSA_PARALLEL_LOOP(loclx,0,locVol)
      // 	{
      // 	  for(int mu=0;mu<NDIM;mu++)
      // 	    for(int nu=0;nu<NDIM;nu++)
      // 	      {
      // 		complex t;
      // 		unsafe_complex_conj1_prod(t,photonEta[loclx][mu],photonEta[loclx][nu]);
      // 		if(loclx==1 and mu==0 and nu==0)
      // 		  MASTER_PRINTF("%lg\n",t[RE]);
      // 		complex_summassign(propRecoMom[loclx][mu][nu],t);
      // 		for(int ri=0;ri<2;ri++)
      // 		  propRecoMom2[loclx][mu][nu][ri]+=sqr(t[ri]);
      // 	      }
      // 	}
      // NISSA_PARALLEL_LOOP_END;
      // set_borders_invalid(propRecoMom);
    }
  
  if(e)
    {
      LxField<quad_su3> u1("u1",WITH_HALO);
      
      PAR(0,
	  locVol,
	  CAPTURE(e,
		  TO_READ(photonField),
		  TO_WRITE(u1)),
	  loclx,
	  {
	    for(int mu=0;mu<NDIM;mu++)
	      {
		complex c;
		complex_iexp(c,e*photonField[loclx][mu][RE]);
		
		su3_put_to_diag_complex(u1[loclx][mu],c);
	      }
	  });
      
      const double pU1=
	global_plaquette_lx_conf(u1);
      
      if(su3Path)
	{
	  LxField<quad_su3> conf("conf",WITH_HALO);
	  read_ildg_gauge_conf(conf,su3Path);
	  
	  const double pSU3=
	    global_plaquette_lx_conf(conf);
	  
	  MASTER_PRINTF("\nPlaquette of the conf without the u1 phase: %.16lg\n",pSU3);
	  
	  PAR(0,
	      locVol,
	      CAPTURE(TO_READ(u1),
		      TO_WRITE(conf)),
	      loclx,
	      {
		for(int mu=0;mu<NDIM;mu++)
		  su3_prodassign_su3(conf[loclx][mu],u1[loclx][mu]);
	      });
	  
	  const double pU3=
	    global_plaquette_lx_conf(conf);
	  
	  MASTER_PRINTF("Plaquette of the conf with the u1 phase: %.16lg\n",pU3);
	  
	  if(u3Path)
	    write_ildg_gauge_conf(u3Path,conf);
	}
      
      MASTER_PRINTF("Plaquette of the u1 phase: %.16lg\n\n",pU1);
    }

  // spin1prop *prop=nissa_malloc("prop",locVol,spin1prop);
  // compute_x_space_tlSym_gauge_propagator_by_fft(prop,gl);
  // for(int d=0;d<2;d++)
  //   for(int mu=0;mu<NDIM;mu++)
  //     for(int nu=0;nu<NDIM;nu++)
  // 	{
  // 	  int ivol,r;
  // 	  get_loclx_and_rank_of_coord(ivol,r,{0,0,d,0});
  // 	  double f=prop[d][mu][nu][RE];
  // 	  if(fabs(f)<1e-15) f=0;
  // 	  if(r==rank) printf("A_mu=%d_nu=%d(t=0,x=0,y=%d,z=0); %.16lg\n",mu,nu,d,f);
  // 	}
  // nissa_free(prop);
  // NISSA_PARALLEL_LOOP(loclx,0,locVol)
  //   {
  //     for(int mu=0;mu<NDIM;mu++)
  // 	for(int nu=0;nu<NDIM;nu++)
  // 	  {
  // 	    complex_prodassign_double(propRecoMom[loclx][mu][nu],1.0/nConfs);
  // 	    complex_prodassign_double(propRecoMom2[loclx][mu][nu],1.0/nConfs);
  // 	    for(int ri=0;ri<2;ri++)
  // 	      {
  // 		propRecoMom2[loclx][mu][nu][ri]-=sqr(propRecoMom[loclx][mu][nu][ri]);
  // 		propRecoMom2[loclx][mu][nu][ri]=sqrt(propRecoMom2[loclx][mu][nu][ri]/(nConfs-1));
  // 	      }
  // 	  }
  //   }
  // NISSA_PARALLEL_LOOP_END;
  // set_borders_invalid(propRecoMom);
  
  // // pass_spin1prop_from_mom_to_x_space(propReco,propRecoMom,gl.bc,true,true);
  
  // // spin1prop *prop=nissa_malloc("prop",locVol,spin1prop);
  // spin1prop *propMom=nissa_malloc("propMom",locVol,spin1prop);
  // compute_mom_space_tlSym_gauge_propagator(propMom,gl);
  // // compute_x_space_tlSym_gauge_propagator_by_fft(prop,gl);
  
  // for(int site=0;site<2;site++)
  //   for(int mu=0;mu<NDIM;mu++)
  //     for(int nu=0;nu<NDIM;nu++)
  // 	MASTER_PRINTF("%d %d %d %lg %lg %lg\n",site,mu,nu,propRecoMom[site][mu][nu][RE],propRecoMom2[site][mu][nu][RE],propMom[site][mu][nu][RE]);
      
  // nissa_free(prop);
  // nissa_free(propMom);
  // nissa_free(propReco);
  // nissa_free(propRecoMom);
  closeNissa();
  
  return 0;
}
