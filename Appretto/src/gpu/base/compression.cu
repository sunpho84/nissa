#pragma once

// ----------------------- double to float and vice-versa conversion ----------------------

void double_to_float_vec(float *out,double *in,int nel)
{for(int iel=0;iel<nel;iel++) out[iel]=(float)(in[iel]);}

void float_to_double_vec(double *out,float *in,int nel)
{for(int iel=0;iel<nel;iel++) out[iel]=(double)(in[iel]);}

inline void double_to_sploat_vec(float *high,float *low,double *in,int nel)
{
  for(int iel=0;iel<nel;iel++)
    {
      high[iel]=(float)(in[iel]);
      low[iel]=(float)(in[iel]-high[iel]);
    }
}

#define sploat_to_double_vec(out,high,low,nel) for(int iel=0;iel<nel;iel++) (out)[iel]=(double)((high)[iel])+(low)[iel];


// -------------------------- su3 matrix compression amd decompression ---------------------

//compress an su3 matrix to 8 doubles
void su3d_to_su3g(su3g out,su3d in)
{
  out[0]=in[0][1][0]; //a2
  out[1]=in[0][1][1];
  out[2]=in[0][2][0]; //a3
  out[3]=in[0][2][1];
  out[4]=atan2(in[0][0][1],in[0][0][0]); //theta_a1
  out[5]=atan2(in[2][0][1],in[2][0][0]); //theta_a2
  out[6]=in[1][0][0]; //b1
  out[7]=in[1][0][1];
}

//compress an su3 matrix to 8 doubles, separating each into 2 floats 
void su3d_to_su3cc(su3c high,su3c low,su3d in)
{
  su3g temp;
  su3d_to_su3g(temp,in);
  double_to_sploat_vec(high,low,temp,8);
}

//compress an su3 matrix to 8 floats
void su3d_to_su3c(su3c out,su3d in)
{
  su3g temp;
  su3d_to_su3g(temp,in);
  double_to_float_vec(out,temp,8);
}

//decompress 8 float parameters to an su3 matrix
void su3c_to_su3f(su3f out,su3c in)
{
#include "gauge_decompression_float.c"
} 

//gpu version
__device__ void gpu_su3c_to_su3f(su3f out,su3c in)
{
#include "gauge_decompression_float.c"
} 

//decompress 8 split double parameters to an su3 double matrix
void su3cc_to_su3d(su3d out,su3c in_high,su3c in_low)
{
#include "gauge_decompression_double.c"
} 

//gpu version
__device__ void gpu_su3cc_to_su3d(su3d out,su3c in_high,su3c in_low)
{
#include "gauge_decompression_double.c"
} 


// ----------------------- full gaugeconf compression and decompression ----------------------

//compress an array of su3
void gaugeconf_compress(su3c *high,su3c *low,quad_su3d *in,int vol)
{for(int i=0;i<vol;i++) for(int mu=0;mu<4;mu++) su3d_to_su3cc(high[vol*mu+i],low[vol*mu+i],in[i][mu]);}

//reconstruction of a full gauge conf in float
__global__ void gpu_gaugeconf_reconstruct_float(quad_su3f *out,su3c *in,int vol)
{
  int pos=threadIdx.x+blockDim.x*blockIdx.x;
  int mu=pos/vol,ivol=pos-vol*mu;
  
  if(pos<vol*4) gpu_su3c_to_su3f(out[ivol][mu],in[pos]);
}

//reconstruction of a full gauge conf in double
__global__ void gpu_gaugeconf_reconstruct_double(quad_su3d *out,su3c *high,su3c *low,int vol)
{
  int pos=threadIdx.x+blockDim.x*blockIdx.x;
  int mu=pos/vol,ivol=pos-vol*mu;
  
  if(pos<vol*4) gpu_su3cc_to_su3d(out[ivol][mu],high[pos],low[pos]);
}

//Compress the configuration on cpu and send it to the gpu
void gaugeconf_compress_and_move(su3c *gpu_comp_conf_high,su3c *gpu_comp_conf_low,quad_su3d *cpu_uncomp_conf,int vol)
{
  //allocate room for the compressed conf on cpu
  su3c *cpu_comp_conf_high=cpu_allocate_su3c(2*4*vol,"cpu_comp_conf");
  su3c *cpu_comp_conf_low=cpu_comp_conf_high+4*vol;
  
  //compress the configuration and send it
  gaugeconf_compress(cpu_comp_conf_high,cpu_comp_conf_low,cpu_uncomp_conf,vol);
  if(gpu_comp_conf_low==gpu_comp_conf_high+4*vol) //single bunch
    memcpy_cpu_to_gpu(gpu_comp_conf_high,cpu_comp_conf_high,2*4*vol*sizeof(su3c));
  else
    {
      memcpy_cpu_to_gpu(gpu_comp_conf_high,cpu_comp_conf_high,4*vol*sizeof(su3c));
      if(gpu_comp_conf_low!=NULL)
	memcpy_cpu_to_gpu(gpu_comp_conf_low,cpu_comp_conf_low,4*vol*sizeof(su3c));
    }
  
  //free cpu compressed conf
  free(cpu_comp_conf_high);
}

//Decompress the configuration on gpu and send it to the cpu in float
void gaugeconf_decompress_and_move_float(quad_su3f *cpu_uncomp_conf,su3c *gpu_comp_conf,int vol,int blocksize)
{
  //allocate room for the decompressed conf on gpu
  quad_su3f *gpu_uncomp_conf=gpu_allocate_quad_su3f(vol,"gpu_uncomp_conf");
  
  //expand the configuration and move it
  int gridsize=(int)ceil(vol*4.0/blocksize);
  gpu_gaugeconf_reconstruct_float<<<gridsize,blocksize>>>(gpu_uncomp_conf,gpu_comp_conf,vol);
  memcpy_gpu_to_cpu(cpu_uncomp_conf,gpu_uncomp_conf,vol*sizeof(quad_su3f));

  //free the gpu decompressed conf
  gpu_free(gpu_uncomp_conf,"gpu_uncomp_conf");
}

//Decompress the configuration on gpu and send it to the cpu in double
void gaugeconf_decompress_and_move_double(quad_su3d *cpu_uncomp_conf,su3c *gpu_comp_conf_high,su3c *gpu_comp_conf_low,int vol,int blocksize)
{
  //allocate room for the decompressed conf on gpu
  quad_su3d *gpu_uncomp_conf=gpu_allocate_quad_su3d(vol,"gpu_uncomp_conf");
  
  //expand the configuration and move it
  int gridsize=(int)ceil(vol*4.0/blocksize);
  gpu_gaugeconf_reconstruct_double<<<gridsize,blocksize>>>(gpu_uncomp_conf,gpu_comp_conf_high,gpu_comp_conf_low,vol);
  memcpy_gpu_to_cpu(cpu_uncomp_conf,gpu_uncomp_conf,vol*sizeof(quad_su3d));

  //free the gpu decompressed conf
  gpu_free(gpu_uncomp_conf,"gpu_uncomp_conf");
}


// ------------------------------------ testing routines -----------------------------------

extern "C" void test_gaugeconf_compression(quad_su3d *conf,int vol)
{
  int blocksize=192;
  
  //Compress the configuration on cpu and send it to the gpu
  su3c *gpu_comp_conf_high=gpu_allocate_su3c(2*4*vol,"gpu_comp_conf");
  su3c *gpu_comp_conf_low=gpu_comp_conf_high+4*vol;
  gaugeconf_compress_and_move(gpu_comp_conf_high,gpu_comp_conf_low,conf,vol);
  
  //Decompress the configuration on the gpu and get it back as a float
  quad_su3f *cpu_uncomp_conf_float=cpu_allocate_quad_su3f(vol,"cpu_uncomp_conf_float");
  gaugeconf_decompress_and_move_float(cpu_uncomp_conf_float,gpu_comp_conf_high,vol,blocksize);
  
  //Decompress the configuration on the gpu and get it back as a double
  quad_su3d *cpu_uncomp_conf_double=cpu_allocate_quad_su3d(vol,"cpu_uncomp_conf_double");
  gaugeconf_decompress_and_move_double(cpu_uncomp_conf_double,gpu_comp_conf_high,gpu_comp_conf_low,vol,blocksize);
  
  //Now compare the decompressed conf with the original one
  for(int ivol=0;ivol<vol;ivol++)
    for(int idir=0;idir<4;idir++)
      {
	double err_float=0;
	double err_double=0;
	
	for(int i=0;i<3;i++)
	  for(int j=0;j<3;j++)
	    for(int ri=0;ri<2;ri++)
	      {
		double diff_float=conf[ivol][idir][i][j][ri]-cpu_uncomp_conf_float[ivol][idir][i][j][ri];
		double diff_double=conf[ivol][idir][i][j][ri]-cpu_uncomp_conf_double[ivol][idir][i][j][ri];
		err_float+=diff_float*diff_float;
		err_double+=diff_double*diff_double;
	      }
	
	printf("%d %d %lg %lg\n",ivol,idir,err_float,err_double);
      }
  
  //free the gpu compressed conf and cpu decompressed conf
  gpu_free(gpu_comp_conf_high,"gpu_comp_conf");
  free(cpu_uncomp_conf_float);
  free(cpu_uncomp_conf_double);
}

