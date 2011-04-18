#include "include.h"
#include "prace_common.cpp"

int main()
{
  //read the input file
  read_data_list();
  //allocate vector where to load data
  jvec P5_Vmu_P5[nmass][nmass][nmoms][nmoms][4];
  jvec P5_P5[nmass][nmoms][2];
  
  /////////////////////////////////////////////////// data loading ////////////////////////////////////////////////////////
  
  //load two and three points
  int im_s=2,im_l=im_spec;
  for(int im2=0;im2<nmass;im2++)
    for(int ik2=0;ik2<nmoms;ik2++)
      for(int r=0;r<2;r++)
	P5_P5[im2][ik2][r]=read_two_points(ik2,im2,r_spec,im_spec,r,0,0);
  
  for(int im2=0;im2<nmass;im2++)
    for(int ik2=0;ik2<nmoms;ik2++)
      for(int mu=0;mu<4;mu++) //average k-pi with pi-k
	for(int im1=0;im1<nmass;im1++)
	  for(int ik1=0;ik1<nmoms;ik1++)
	    {
	      int ri[4]={0,1,1,1};
	      P5_Vmu_P5[im2][im1][ik2][ik1][mu]=read_three_points(ik2,im2,ik1,im1,mu,ri[mu]);
	      //P5_Vmu_P5[im2][im1][ik2][ik1][mu].print_to_file("/tmp/out_P5_Vmu_P5_%d_%d_%d_%d_%d",im2,ik2,mu,im1,ik1);
	    }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for(int im2=0;im2<nmass;im2++)
    for(int ik2=0;ik2<nmoms;ik2++)
      {
	jvec b=(P5_P5[im2][ik2][0]+P5_P5[im2][ik2][1])/2;
	b.print_to_file("/tmp/P5_P5_%d_%d",im2,ik2);
	jvec c=effective_mass(b.simmetrized(1));
	jack m=constant_fit(c,10,21);
	cout<<im2<<" "<<ik2<<" "<<m<<endl;
      }
  
  //calculate Zv
  jvec Zv[nmass];
  for(int imass=0;imass<nmass;imass++)
    Zv[imass]=P5_P5[imass][ith_spec][r_spec].data[TH]/(2*P5_Vmu_P5[imass][imass][ith_spec][ith_spec][0]);
  
  cout<<Zv[0]<<endl;
  exit(1);
  
  //fit mass from 2pts
  jack E2pt[nmass][nmoms],Z;
  for(int imass=0;imass<nmass;imass++)
    for(int ik=0;ik<nmoms;ik++)
      E2pt[imass][ik]=constant_fit(effective_mass(P5_P5[imass][ik][r_spec].simmetrized(1)),12,23);
  
  //fit mass from 3 pts
  jack E3pt[nmass][nmoms];
  for(int imass=0;imass<nmass;imass++)
    for(int ik=0;ik<nmoms;ik++)
      {
	jack X=constant_fit((P5_Vmu_P5[imass][imass][ik][ik][1]+P5_Vmu_P5[imass][imass][ik][ik][2]+
			     P5_Vmu_P5[imass][imass][ik][ik][3])/P5_Vmu_P5[imass][imass][ik][ik][0],10,14)/sqrt(3);
	double p=sqrt(3)*2*M_PI*theta[ik]/L;
	E3pt[imass][ik]=p/X;//sqrt(1.0/(X*X)-1);
	cout<<imass<<" "<<ik<<" "<<E2pt[imass][ik]<<" "<<E3pt[imass][ik]<<endl;
      }
  
  //calculate M2_K-M2_Pi
  jack DM2=pow(E2pt[im_s][0],2)-pow(E2pt[im_l][0],2);
  
  //Form factor
  ofstream file_fp("/tmp/fp"),file_fm("/tmp/fm"),file_f0("/tmp/f0");
  for(int ik_l=0;ik_l<nmoms;ik_l++)
    for(int ik_s=0;ik_s<nmoms;ik_s++)
      {
	cout<<"-------"<<endl;
	
	jvec cV[2];
	//temporal direction: simmetrize t with T-t
	cV[0]=P5_Vmu_P5[im_s][im_l][ik_s][ik_l][0].simmetrized(-1)*P5_Vmu_P5[im_l][im_s][ik_l][ik_s][0].simmetrized(-1);
	//spatial direction
	cV[1]=cV[0]*0;
	for(int idir=1;idir<4;idir++) cV[1]+=P5_Vmu_P5[im_s][im_l][ik_s][ik_l][idir].simmetrized(1)*
					P5_Vmu_P5[im_l][im_s][ik_l][ik_s][idir].simmetrized(1)/3;
	//denominator
	jvec Den=P5_Vmu_P5[im_l][im_l][ik_l][ik_l][0].simmetrized(-1)*P5_Vmu_P5[im_s][im_s][ik_s][ik_s][0].simmetrized(-1);
	
	//calculate 2*sqrt(E_K*E_Pi)
	jack A=2*sqrt(E2pt[im_s][ik_s]*E2pt[im_l][ik_l]);
	//calculate the two ratios
	for(int mu=0;mu<2;mu++)
	  {
	    cV[mu]=sqrt(cV[mu]/Den);
            //reassign the correct sign
	    if(P5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu][T/4][njack]<0) cV[mu]=-cV[mu];
	    //multiply by 2*sqrt(E_K*E_Pi)
	    cV[mu]*=A;
	    cV[mu].print_to_file("/tmp/out_double_ratio_simm_P5_V%d_P5_%d_%d",mu,ik_s,ik_l);
	  }
	
	jack V0=constant_fit(cV[0],10,14);
	jack Vi=constant_fit(cV[1],10,14);
	
	cout<<"V0: "<<V0/A<<" "<<A<<endl;
	cout<<"Vi: "<<Vi/A<<endl;
	
	//calculate the impulses
	jack Q[4],P[4],Q2;
	calculate_Q_P(Q,P,E2pt[im_s][ik_s],theta[ik_s],E2pt[im_l][ik_l],theta[ik_l]);
	Q2=quad_jack_prod_quad_jack(Q,Q);
	//calculate_Q_P(Q,P,E3pt[im_s][ik_s],theta[ik_s],E3pt[im_l][ik_l],theta[ik_l]);

	jack fp,fm,f0;
	if(theta[ik_s]!=-theta[ik_l] && ik_s!=ik_l)
	  {
	    jack DET=P[0]*Q[1]-Q[0]*P[1];
	    //calculate (f+) and (f-)*Q2/(M_K^2-M_Pi^2)
	    fp=(V0*Q[1]-Vi*Q[0])/DET;
	    fm=(P[0]*Vi-P[1]*V0)/DET;
	    f0=fp+Q2/DM2*fm;
	  }
	if(theta[ik_s]==-theta[ik_l] && ik_s!=ik_l)
	  {
	    //calculate (f+) and (f-)*Q2/(M_K^2-M_Pi^2)
	    fm=Vi/Q[1];
	    fp=(V0-Q[0]*fm)/P[0];
	    f0=fp+Q2/DM2*fm;
	  }
	if(ik_s!=ik_l)
	  {
	    cout<<ik_s<<" "<<ik_l<<" Q0="<<Q[0].med()<<" "<<"Qi="<<Q[1].med()<<" P0="<<P[0].med()<<" P1="<<P[1].med()<<endl;
	    file_fp<<Q2.med()<<" "<<fp<<" "<<ik_s<<" "<<ik_l<<endl;
	    file_fm<<Q2.med()<<" "<<fm<<" "<<ik_s<<" "<<ik_l<<endl;
	    file_f0<<Q2.med()<<" "<<f0<<" "<<ik_s<<" "<<ik_l<<endl;
	  }
      }
  
  return 0;
}
