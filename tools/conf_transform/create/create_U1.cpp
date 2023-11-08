#include <nissa.hpp>

using namespace nissa;

typedef double realspin1field[NDIM];
typedef realspin1field realspin1prop[NDIM];

void inMain(int narg,char **arg)
{
  //check argument
  if(narg<7) crash("Use: %s L T gauge[L=LANDAU,F=FEYNMAN] seed nconfs pattern",arg[0]);
  
  const int L=atoi(arg[1]);
  const int T=atoi(arg[2]);
  const char gauge=arg[3][0];
  const int seed=atoi(arg[4]);
  const int nConfs=atoi(arg[5]);
  const char* pattern=arg[6];
  
  init_grid(T,L);
  
  FieldRngStream fieldRngStream;
  fieldRngStream.init(seed);
  
  gauge_info gl;
  gl.zms=UNNO_ALEMANNA;
  switch(gauge)
    {
    case 'L':
      gl.alpha=LANDAU_ALPHA;
      break;
    case 'F':
      gl.alpha=FEYNMAN_ALPHA;
      break;
    default:
      crash("unknown gauge '%c'",gauge);
    }
  gl.c1=WILSON_C1;
  
  spin1field *photonEta=nissa_malloc("photonEta",locVol,spin1field);
  spin1field *photonField=nissa_malloc("photonField",locVol,spin1field);
  realspin1field *realPhotonField=nissa_malloc("realPhotonField",locVol,realspin1field);
  // spin1prop *propRecoMom=nissa_malloc("propRecoMom",locVol,spin1prop);
  // spin1prop *propRecoMom2=nissa_malloc("propRecoMom2",locVol,spin1prop);
  // spin1prop *propReco=nissa_malloc("propReco",locVol,spin1prop);
  // vector_reset(propRecoMom);
  // vector_reset(propRecoMom2);
  
  for(int iConf=0;iConf<nConfs;iConf++)
    {
      auto sourceFiller=
	fieldRngStream.getDrawer<spin1field>();
      
      sourceFiller.fillField(photonEta);
      
      NISSA_PARALLEL_LOOP(loclx,0,locVol)
	{
	  for(int mu=0;mu<NDIM;mu++)
	    z2Transform(photonEta[loclx][mu]);
	}
      NISSA_PARALLEL_LOOP_END;
      set_borders_invalid(photonEta);
      
      multiply_by_sqrt_tlSym_gauge_propagator(photonField,photonEta,gl);
      
      NISSA_PARALLEL_LOOP(loclx,0,locVol)
	{
	  for(int mu=0;mu<NDIM;mu++)
	    realPhotonField[loclx][mu]=photonField[loclx][mu][RE];
	}
      NISSA_PARALLEL_LOOP_END;
      set_borders_invalid(realPhotonField);
      
      char outPath[128];
      snprintf(outPath,128,pattern,iConf);
      if(strncasecmp(outPath,pattern,128)==0)
	crash("pattern %s cannot be used to creat the conf name, try something like: conf%%04d.dat",pattern);
      
      ILDG_File file=ILDG_File_open_for_write(outPath);
      char ildg_format_message[1024];
      sprintf(ildg_format_message,
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
      
      write_real_vector(file,realPhotonField,64,"ildg-binary-data");
      
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
      // 		  master_printf("%lg\n",t[RE]);
      // 		complex_summassign(propRecoMom[loclx][mu][nu],t);
      // 		for(int ri=0;ri<2;ri++)
      // 		  propRecoMom2[loclx][mu][nu][ri]+=sqr(t[ri]);
      // 	      }
      // 	}
      // NISSA_PARALLEL_LOOP_END;
      // set_borders_invalid(propRecoMom);
    }
  
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
  // 	master_printf("%d %d %d %lg %lg %lg\n",site,mu,nu,propRecoMom[site][mu][nu][RE],propRecoMom2[site][mu][nu][RE],propMom[site][mu][nu][RE]);
      
  // nissa_free(prop);
  // nissa_free(propMom);
  // nissa_free(propReco);
  // nissa_free(propRecoMom);
  nissa_free(realPhotonField);
  nissa_free(photonField);
  nissa_free(photonEta);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,inMain);
  close_nissa();
  
  return 0;
}
