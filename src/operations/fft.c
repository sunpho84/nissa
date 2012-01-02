#pragma once

//return the bit inverse of an int
int bitrev(int in,int l2n)
{
  int out=0;
  
  for(int i=0;i<l2n;i++) if(in & (1<<i)) out+=(1<<(l2n-i-1));
  
  return out;
}

//return the powers of n contained in the input
int find_max_pow2(int a)
{
  int nl=0;
  while((a&0x1)==0)
    {
      nl++;
      a>>=1;
    }
  
  return nl;
}

//transpose the data shifting coordinate order by 1
//if normal ordering TXYZ is present and you pass mu0=0,
//it will put data in XYZT order. Than if called with mu0=1
//it will put data in YZTX order, etc
void data_coordinate_order_shift(complex *data,int ncpp,int mu0)
{
  int *pos=nissa_malloc("Pos",loc_vol,int);
  
  //order of directions
  int in_mu[4] ={mu0,(mu0+1)%4,(mu0+2)%4,(mu0+3)%4};
  int out_mu[4]={(mu0+1)%4,(mu0+2)%4,(mu0+3)%4,mu0};
  
  //loop over all sites
  int x[4];
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    //find in and out position
	    int in=0,out=0;
	    for(int mu=0;mu<4;mu++)
	      {
		in=in*loc_size[in_mu[mu]]+x[in_mu[mu]];
		out=out*loc_size[out_mu[mu]]+x[out_mu[mu]];
	      }
	   
	    //mark its final order
	    pos[in]=out;
	  }
  
  //apply reordering
  reorder_vector((char*)data,pos,loc_vol,sizeof(complex)*ncpp);
  
  nissa_free(pos);
}

// Perform the 1d transform (mu is used to choose the communicator and n)
// The fft consist of four step:
//  1) data of is rearranged across the ranks
//  2) ft is done on the odd length blocks (not needed if data length is a power of 2)
//  3) Lanczos lemma is implemented on the local data
//  4) Lanczos lemma is applied on non-local data (if needed)
void fft1d(complex *out,complex *in,int ncpp,int mu,double sign,int normalize)
{
  int log2glb_nblk=find_max_pow2(glb_size[mu]); //number of powers of 2 contained in n
  int glb_nblk=1<<log2glb_nblk;                 //remaining odd block
  int blk_size=glb_size[mu]/glb_nblk;           //block size
  int loc_nblk=glb_nblk/nrank_dir[mu];          //number of local blocks
  
  //check if the loc_size is a multiple of block_size
  if(loc_size[mu]%blk_size!=0) crash("Error, FFT implemented only for power of 2 grids!");
  
  //allocate the buffer to send data
  complex *buf=nissa_malloc("buf",loc_size[mu]*ncpp,complex);
  
  ////////////////////////////// first part /////////////////////////////////
  
  //reorder data across the various rank: glb_iel_in=iel_blk_out*glb_nblk+glb_iblk_out_rev
  int nrequest0=0,nexp_request0=loc_size[mu]*2;
  MPI_Request request0[nexp_request0];
  for(int loc_iblk_in=0;loc_iblk_in<loc_nblk;loc_iblk_in++)
    for(int iel_blk_in=0;iel_blk_in<blk_size;iel_blk_in++)
      {
	//find the full index
	int loc_iel_in=loc_iblk_in*blk_size+iel_blk_in;
	int glb_iel_in=glb_coord_of_loclx[0][mu]+loc_iel_in;
	int glb_iblk_in=rank_coord[mu]*loc_nblk+loc_iblk_in;
	
	//look at the output
	int iel_blk_out=glb_iel_in/glb_nblk;
	int glb_iblk_out=bitrev(glb_iel_in-iel_blk_out*glb_nblk,log2glb_nblk);
	int rank_coord_out=glb_iblk_out/loc_nblk;
	int loc_iblk_out=glb_iblk_out-rank_coord_out*loc_nblk;
	int loc_iel_out=loc_iblk_out*blk_size+iel_blk_out;
	
	//now, if the destination is local, put it on place
	if(rank_coord_out==rank_coord[mu]) memcpy(buf+ncpp*loc_iel_out,in+ncpp*loc_iel_in,sizeof(complex)*ncpp);
	else //send the data to the destination
	  MPI_Isend((void*)(in+loc_iel_in*ncpp),ncpp*2,MPI_DOUBLE,rank_coord_out,1241+loc_iel_out,
		    line_comm[mu],&request0[nrequest0++]);

	//if necessary receive data: glb_iel_send=iel_blk_in*glb_nblk+glb_iblk_in_rev
	int glb_iel_send=iel_blk_in*glb_nblk+bitrev(glb_iblk_in,log2glb_nblk);
	int rank_coord_send=glb_iel_send/loc_size[mu];
	if(rank_coord_send!=line_rank[mu])
	  MPI_Irecv((void*)(buf+loc_iel_in*ncpp),ncpp*2,MPI_DOUBLE,rank_coord_send,1241+loc_iel_in, 
		    line_comm[mu],&request0[nrequest0++]);
      }
  
  //wait for scramble to finish
  MPI_Status status0[nexp_request0];
  if(nrequest0>0) MPI_Waitall(nrequest0,request0,status0);
  
  /////////////////////////// second part ////////////////////////////////////
  
  //perform the block fourier transform if needed
  if(blk_size==1) memcpy(out,buf,loc_size[mu]*ncpp*sizeof(complex));
  else
    {
      //initialize the output
      memset(out,0,loc_size[mu]*ncpp*sizeof(complex));

      //loop over local blocks
      for(int loc_iblk=0;loc_iblk<loc_nblk;loc_iblk++)
	{
	  //take initial position of the output and input block
	  complex *in_blk=buf+loc_iblk*blk_size*ncpp;
	  complex *out_blk=out+loc_iblk*blk_size*ncpp;
	  
	  //loop over out elements of out block
	  for(int iel_out=0;iel_out<blk_size;iel_out++)
	    {
	      //take initial position of out local elements
	      complex *out_pad=out_blk+iel_out*ncpp;
	      
	      //incrementing factor
	      double theta=iel_out*sign*2*M_PI/blk_size;
	      double wtemp=sin(0.5*theta);
	      double wpr=-2*wtemp*wtemp;
	      double wpi=sin(theta);
	      
	      //fourier factor
	      double wr=1;
	      double wi=0;
	      
	      //loop over elements of in blocks
	      for(int iel_in=0;iel_in<blk_size;iel_in++)	  
		{
		  //take initial position of in local elements
		  complex *in_pad=in_blk+iel_in*ncpp;
		  
		  //loop over local elements
		  for(int ic=0;ic<ncpp;ic++)
		    {
		      out_pad[ic][0]+=wr*in_pad[ic][0]-wi*in_pad[ic][1];
		      out_pad[ic][1]+=wr*in_pad[ic][1]+wi*in_pad[ic][0];
		    }
		  
		  //increment the twiddle
		  wtemp=wr;
		  wr+=wr*wpr-   wi*wpi;
		  wi+=wi*wpr+wtemp*wpi;
		}
	    }
	}
    }
  
  /////////////////////////////  third part ////////////////////////////
  
  //now perform the lanczos procedure up to when it does not need communications
  for(int delta=blk_size;delta<loc_size[mu];delta*=2)
    {
      //incrementing factor
      double theta=sign*2*M_PI/(2*delta);
      double wtemp=sin(0.5*theta);
      double wpr=-2*wtemp*wtemp;
      double wpi=sin(theta);
      
      //fourier coefficient
      double wr=1;
      double wi=0;

      //loop over the delta length (each m will correspond to increasing twiddle)
      for(int m=0;m<delta;m++)
        {
	  //loop on the first addend
          for(int i=m;i<loc_size[mu];i+=2*delta)
            {
	      //second addend
              int j=i+delta;
              
	      //site data multiplication
	      for(int ic=0;ic<ncpp;ic++)
		{
		  double tempr=wr*out[j*ncpp+ic][0]-wi*out[j*ncpp+ic][1];
		  double tempi=wr*out[j*ncpp+ic][1]+wi*out[j*ncpp+ic][0];
		  
		  out[j*ncpp+ic][0]=out[i*ncpp+ic][0]-tempr;
		  out[j*ncpp+ic][1]=out[i*ncpp+ic][1]-tempi;

		  out[i*ncpp+ic][0]+=tempr;
		  out[i*ncpp+ic][1]+=tempi;
		}
            }
          wtemp=wr;
          wr+=wr*wpr-   wi*wpi;
          wi+=wi*wpr+wtemp*wpi;
	}
    }
  
  ///////////////////////////////////// fourth part //////////////////////////////
  
  //now perform the lanczos procedure up to the end
  for(int delta_rank=1;delta_rank<nrank_dir[mu];delta_rank*=2) //this is rank width of delta
    {
      int delta=delta_rank*loc_size[mu];    //block extent
      int idelta=rank_coord[mu]/delta_rank; //identify the delta of the current rank
      
      //find if the current rank holding the first or second block
      complex *first,*second;
      MPI_Request request2[2];
      MPI_Status status2[2];
      
      if(idelta%2==0) //first: so it has to receive second block and send first
	{
	  first=out;
	  second=buf;
	  
	  MPI_Irecv((void*)buf,2*ncpp*loc_size[mu],MPI_DOUBLE,line_rank[mu]+delta_rank,113+line_rank[mu],
		    line_comm[mu],&(request2[0]));
	  MPI_Isend((void*)out,2*ncpp*loc_size[mu],MPI_DOUBLE,line_rank[mu]+delta_rank,113+line_rank[mu]+delta_rank,
		    line_comm[mu],&(request2[1]));
	}
      else           //second: so it has to receive first block and send second
	{
	  first=buf;
	  second=out;
	  
	  MPI_Irecv((void*)buf,2*ncpp*loc_size[mu],MPI_DOUBLE,line_rank[mu]-delta_rank,113+line_rank[mu],
		    line_comm[mu],&(request2[0]));
	  MPI_Isend((void*)out,2*ncpp*loc_size[mu],MPI_DOUBLE,line_rank[mu]-delta_rank,113+line_rank[mu]-delta_rank,
		    line_comm[mu],&(request2[1]));
	}
      
      //incrementing factor
      double theta=sign*2*M_PI/(2*delta);
      double wtemp=sin(0.5*theta);
      double wpr=-2*wtemp*wtemp;
      double wpi=sin(theta);
      
      int rank_pos_delta=rank_coord[mu]%delta_rank;  //position of the rank inside the delta
      int pos_delta=rank_pos_delta*loc_size[mu];     //starting coord of local delta inside the delta

      //fourier coefficient
      double wr=cos(pos_delta*theta);
      double wi=sin(pos_delta*theta);      
      
      //wait for communications to finish
      MPI_Waitall(2,request2,status2);

      //loop over the delta length (each m will correspond to increasing twiddle)
      for(int loc_m=0;loc_m<loc_size[mu];loc_m++)
        {
	  //site data multiplication
	  for(int ic=0;ic<ncpp;ic++)
	    {
	      double tempr=wr*second[loc_m*ncpp+ic][0]-wi*second[loc_m*ncpp+ic][1];
	      double tempi=wr*second[loc_m*ncpp+ic][1]+wi*second[loc_m*ncpp+ic][0];
	      
	      if(idelta%2==0)
		{
		  first[loc_m*ncpp+ic][0]+=tempr;
		  first[loc_m*ncpp+ic][1]+=tempi;
		}
	      else
		{
		  second[loc_m*ncpp+ic][0]=first[loc_m*ncpp+ic][0]-tempr;
		  second[loc_m*ncpp+ic][1]=first[loc_m*ncpp+ic][1]-tempi;
		}
	    }
          wtemp=wr;
          wr+=wr*wpr-   wi*wpi;
          wi+=wi*wpr+wtemp*wpi;
	}
    }
  
  if(normalize==1)
    for(int i=0;i<loc_size[mu]*ncpp;i++)
      for(int ri=0;ri<2;ri++) out[i][ri]/=glb_size[mu];
  
  nissa_free(buf);
}

//perform the fft in all directions
void fft4d(complex *out,complex *in,int ncpp,double sign,int normalize)
{
  //copy input in the output (if they differ!)
  if(out!=in) memcpy(out,in,ncpp*sizeof(complex)*loc_vol);

  //perform the fft in each direction
  for(int mu=0;mu<4;mu++)
    {
      //perform the 1d fft (slower dir)
      fft1d(out,out,ncpp*loc_vol/loc_size[mu],mu,sign,normalize);
      
      //for the time being we stick to transpose the data
      data_coordinate_order_shift(out,ncpp,mu);
    }
}

