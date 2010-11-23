#pragma once

#include <iostream>

#include "global.cpp"

using namespace std;

//The structure for gamma matrix
struct gamma
{
public:
  int pos[4];
  complex entr[4];
};

//The base of the 16 gamma matrixes and the two rotators
gamma base_gamma[16];
gamma Pplus,Pminus;

//Initialize a gamma with outside entries
void init_gamma(gamma &out,
		int pos0,double rea0,double ima0,
		int pos1,double rea1,double ima1,
		int pos2,double rea2,double ima2,
		int pos3,double rea3,double ima3)
{
  out.pos[0]=pos0;
  out.pos[1]=pos1;
  out.pos[2]=pos2;
  out.pos[3]=pos3;

  out.entr[0][0]=rea0;
  out.entr[1][0]=rea1;
  out.entr[2][0]=rea2;
  out.entr[3][0]=rea3;

  out.entr[0][1]=ima0;
  out.entr[1][1]=ima1;
  out.entr[2][1]=ima2;
  out.entr[3][1]=ima3;
}

//If the two gamma in1 and in2 have the same position structure sum
//them, otherwise it crashes
void gamma_summ(gamma &out,gamma &in1,gamma &in2)
{
  for(int ig=0;ig<4;ig++)
    if(in1.pos[ig]==in2.pos[ig])
      {
	out.pos[ig]=in1.pos[ig];
	complex_summ(out.entr[ig],in1.entr[ig],in2.entr[ig]);
      }
    else 
      {
	if(rank==0) cerr<<"The two matrix passed to sum have different positions"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
}

//Assign to the first gamma the product of the second and the third
void gamma_prod(gamma &out,gamma &in1,gamma &in2)
{
  //This is the line on the first matrix
  for(int ig1=0;ig1<4;ig1++)
    {
      //This is the line to be taken on the second matrix
      int ig2=in1.pos[ig1];

      //For each line, then column of the output matrix which is
      //different from 0 is the column of the second matrix different
      //from 0 on the line with index equal to the column of the first
      //matrix which is different from 0 (that is, ig2)
      out.pos[ig1]=in2.pos[ig2];

      //The entries of the output is, on each line, the complex
      //product of the entries of the first matrix on that line, for
      //the entries of the second matrix on the line with the index
      //equal to the column of the first matrix which is different
      //from 0 (which again is ig2)
      complex_prod(out.entr[ig1],in1.entr[ig1],in2.entr[ig2]);
    }
}

//Print the gamma on node 0
void print_gamma(gamma &in)
{
  for(int ir=0;ir<4;ir++)
    {
      int pos=in.pos[ir];
      for(int ic=0;ic<pos;ic++) cout<<"0,0\t";
      cout<<in.entr[ir][0]<<","<<in.entr[ir][1]<<"\t";
      for(int ic=pos+1;ic<4;ic++) cout<<"0,0\t";
      cout<<endl;
    }
}

//Initialize the gamma matrix base and the rotators
// base_gamma[ 0] = identity
// base_gamma[ 1] = gamma_1
// base_gamma[ 2] = gamma_2
// base_gamma[ 3] = gamma_3
// base_gamma[ 4] = gamma_0
// base_gamma[ 5] = gamma_5
// base_gamma[ 6] = gamma_15
// base_gamma[ 7] = gamma_25
// base_gamma[ 8] = gamma_35
// base_gamma[ 9] = gamma_05
// base_gamma[10] = gamma_01
// base_gamma[11] = gamma_02
// base_gamma[12] = gamma_03
// base_gamma[13] = gamma_23
// base_gamma[14] = gamma_31
// base_gamma[15] = gamma_12

void init_base_gamma()
{
  init_gamma(base_gamma[ 0],  0,1,0  , 1,1,0  , 2,1,0  , 3,1,0 );
  init_gamma(base_gamma[ 1],  3,0,-1 , 2,0,-1 , 1,0,1  , 0,0,1 );
  init_gamma(base_gamma[ 2],  3,-1,0 , 2,1,0  , 1,1,0  , 0,-1,0);
  init_gamma(base_gamma[ 3],  2,0,-1 , 3,0,1  , 0,0,1  , 1,0,-1);
  init_gamma(base_gamma[ 4],  2,-1,0 , 3,-1,0 , 0,-1,0 , 1,-1,0);
  init_gamma(base_gamma[ 5],  0,1,0  , 1,1,0  , 2,-1,0 , 3,-1,0);
  init_gamma(base_gamma[ 6],  3,0,1  , 2,0,1  , 1,0,1  , 0,0,1 );
  init_gamma(base_gamma[ 7],  3,1,0  , 2,-1,0 , 1,1,0  , 0,-1,0);
  init_gamma(base_gamma[ 8],  2,0,1  , 3,0,-1 , 0,0,1  , 1,0,-1);
  init_gamma(base_gamma[ 9],  2,1,0  , 3,1,0  , 0,-1,0 , 1,-1,0);
  init_gamma(base_gamma[10],  1,0,-1 , 0,0,-1 , 3,0,1  , 2,0,1 );
  init_gamma(base_gamma[11],  1,-1,0 , 0,1,0  , 3,1,0  , 2,-1,0);
  init_gamma(base_gamma[12],  0,0,-1 , 1,0,1  , 2,0,1  , 3,0,-1);
  init_gamma(base_gamma[13],  1,0,1  , 0,0,1  , 3,0,1  , 2,0,1 );
  init_gamma(base_gamma[14],  1,1,0  , 0,-1,0 , 3,1,0  , 2,-1,0);
  init_gamma(base_gamma[15],  0,0,1  , 1,0,-1 , 2,0,1  , 3,0,-1);
}
