#include "include.h"
#include "prace_common.cpp"

int main()
{
  //read the input file
  read_data_list();
  //allocate vector where to load data
  jvec P5_Vmu_P5[nmass][nmass][nmoms][nmoms][4];
  jvec P5_P5[nmass][nmoms];
  
  int im_s=3,im_l=im_spec;
  
  /////////////////////////////////////////////////// data loading ////////////////////////////////////////////////////////
  
  //load two points
  int ntwo_to_load=nmass*nmoms;
  cout<<"Now loading "<<ntwo_to_load<<" 2pts correlations functions"<<endl;
  for(int im2=0;im2<nmass;im2++)
    for(int ik2=0;ik2<nmoms;ik2++)
      P5_P5[im2][ik2]=read_two_points(ik2,im2,r_spec,im_spec,r_spec,0,0).simmetrized(1);

  //average theta with -theta
  for(int im2=0;im2<nmass;im2++)
    for(int ik2=0;ik2<nmoms;ik2++)
      P5_P5[im2][ik2]=P5_P5[im2][iopp_th[ik2]]=(P5_P5[im2][ik2]+P5_P5[im2][iopp_th[ik2]])/2;
  
  //load three points
  FILE *cache_3pts=fopen("3pts_cache","r");
  bool use_cache=(cache_3pts!=NULL);
  if(use_cache) cout<<"Cache file existing!"<<endl;
  else 
    {
      cache_3pts=fopen("3pts_cache","w");
      cout<<"Creating cache!"<<endl;
    }
  int ithree=0,nthree_to_load=nmass*nmoms*nmass*nmoms*4;
  double cache_buf[(TH+1)*(njack+1)];
  cout<<"Now loading "<<nthree_to_load<<" 3pts correlations functions"<<endl;
  for(int im2=0;im2<nmass;im2++)
    for(int ik2=0;ik2<nmoms;ik2++)
      for(int mu=0;mu<4;mu++)
	for(int im1=0;im1<nmass;im1++)
	  for(int ik1=0;ik1<nmoms;ik1++)
	    {
	      int pa[4]={-1,1,1,1};
	      int ri[4]={0,1,1,1};
	      if(!use_cache)
		{
		  P5_Vmu_P5[im2][im1][ik2][ik1][mu]=-read_three_points(ik2,im2,ik1,im1,mu,ri[mu]).simmetrized(pa[mu]);
		  P5_Vmu_P5[im2][im1][ik2][ik1][mu].get(cache_buf);
		  fwrite(cache_buf,sizeof(double),(TH+1)*(njack+1),cache_3pts);
		}
	      else
		{
		  fread(cache_buf,sizeof(double),(TH+1)*(njack+1),cache_3pts);
		  P5_Vmu_P5[im2][im1][ik2][ik1][mu]=jvec(TH+1,njack);
		  P5_Vmu_P5[im2][im1][ik2][ik1][mu].put(cache_buf);
		}
	      //P5_Vmu_P5[im2][im1][ik2][ik1][mu].print_to_file("/tmp/out_P5_Vmu_P5_%d_%d_%d_%d_%d",im2,ik2,mu,im1,ik1);
	      
	      ithree++;
	      if(ithree%(nthree_to_load/10)==0) cout<<"Loading completed at "<<ithree*100/nthree_to_load<<"%"<<endl;
	    }
  
  /////////////////////////////////// perform improvements on three points //////////////////////////////////////
  
  //average K with Pi
  cout<<"Performing K-Pi and Pi-K improvemnt"<<endl;
  for(int im2=0;im2<nmass;im2++)
    for(int ik2=0;ik2<nmoms;ik2++)
      for(int mu=0;mu<4;mu++)
	for(int im1=0;im1<nmass;im1++)
	  for(int ik1=0;ik1<nmoms;ik1++)
	    {
	      jvec a=(P5_Vmu_P5[im2][im1][ik2][ik1][mu]+
		      P5_Vmu_P5[im1][im2][ik1][ik2][mu].simmetric())/2;
	      
	      P5_Vmu_P5[im2][im1][ik2][ik1][mu]=a;
	      P5_Vmu_P5[im1][im2][ik1][ik2][mu]=a.simmetric();
	    }
  
  //perform theta improvement
  cout<<"Performing theta improvemnt"<<endl;
  for(int im2=0;im2<nmass;im2++)
    for(int ik2=0;ik2<nmoms;ik2++)
      for(int mu=0;mu<4;mu++)
	for(int im1=0;im1<nmass;im1++)
	  for(int ik1=0;ik1<nmoms;ik1++)
	    {
	      int sig[4]={+1,-1,-1,-1};
	      jvec a=(P5_Vmu_P5[im2][im1][ik2][ik1][mu]+sig[mu]*
		      P5_Vmu_P5[im2][im1][iopp_th[ik2]][iopp_th[ik1]][mu])/2;
	      
	      P5_Vmu_P5[im2][im1][ik2][ik1][mu]=a;
	      P5_Vmu_P5[im2][im1][iopp_th[ik2]][iopp_th[ik1]][mu]=sig[mu]*a;
	    }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //calculate Zv
  jvec Zv[nmass];
  for(int imass=0;imass<nmass;imass++)
    Zv[imass]=P5_P5[imass][ith_spec].data[TH]/(2*P5_Vmu_P5[imass][imass][ith_spec][ith_spec][0]);
  
  P5_P5[0][ith_spec].print_to_file("/tmp/P5P5Zv0");
  P5_Vmu_P5[0][0][ith_spec][ith_spec][0].print_to_file("/tmp/P5V0P5Zv0");
  Zv[0].print_to_file("Zv0");
  
  //fit mass from 2pts
  jack E2pt[nmass][nmoms],Z;
  for(int imass=0;imass<nmass;imass++)
    for(int ik=0;ik<nmoms;ik++)
      {
	jvec temp=effective_mass(P5_P5[imass][ik]);
	
	int tmin;
	if(imass==0) tmin=11;
	else if(imass<4) tmin=13;
	else tmin=15;
	
	E2pt[imass][ik]=constant_fit(temp,tmin,23);
	
	//create P5P5 graph
	grace out("P5_P5_%d_%d",imass,ik);
	//create polygon
	const int npoint=2;
	out.set(1,"grey");
	double X_pol[npoint]={tmin,23};
	jvec Y_pol(npoint,njack);
	Y_pol=E2pt[imass][ik];
	out.polygon(X_pol,Y_pol);
	//average
	out.new_set();
	out.set(1,"black");
	out.ave_line(X_pol,Y_pol);
	//create data
	out.new_set();
	out.set(2,"none","black");
	out.fout<<"@type xydy"<<endl;
	out<<temp;
      }
  
  //fit mass from 3 pts
  jack E3pt[nmass][nmoms];
  ofstream out_mass("mass");
  for(int imass=0;imass<nmass;imass++)
    for(int ik=0;ik<nmoms;ik++)
      {
	jvec a=(P5_Vmu_P5[imass][imass][ik][ik][1]+P5_Vmu_P5[imass][imass][ik][ik][2]+
		P5_Vmu_P5[imass][imass][ik][ik][3])/P5_Vmu_P5[imass][imass][ik][ik][0]/sqrt(3);
	for(int mu=0;mu<4;mu++) P5_Vmu_P5[imass][imass][ik][ik][mu].print_to_file("/tmp/3pts_%d_%d_%d",imass,ik,mu);
	jack X=constant_fit(a,10,14);
	double p=sqrt(3)*M_PI*theta[ik]/L;
	
	if(p!=0) E3pt[imass][ik]=p/X;
	else     E3pt[imass][ik]=E2pt[imass][ik];
	
	out_mass<<imass<<" "<<ik<<" "<<E2pt[imass][ik]<<" "<<E3pt[imass][ik]<<endl;
      }
  out_mass.close();
  
  //calculate M2_K-M2_Pi
  jack DM2=pow(E2pt[im_s][ith_spec],2)-pow(E2pt[im_l][ith_spec],2);
  
  //Form factor
  ofstream file_fp("/tmp/fp"),file_fm("/tmp/fm"),file_f0("/tmp/f0");
  for(int ik_l=0;ik_l<nmoms;ik_l++)
    for(int ik_s=0;ik_s<nmoms;ik_s++)
      {
	cout<<"-------"<<endl;
	
	jvec cV[2];
	//temporal direction
	cV[0]=P5_Vmu_P5[im_s][im_l][ik_s][ik_l][0]*P5_Vmu_P5[im_l][im_s][ik_l][ik_s][0];
	//spatial direction
	cV[1]=cV[0]*0;
	for(int idir=1;idir<4;idir++) cV[1]+=P5_Vmu_P5[im_s][im_l][ik_s][ik_l][idir]*
					P5_Vmu_P5[im_l][im_s][ik_l][ik_s][idir]/3;
	//denominator
	jvec Den=P5_Vmu_P5[im_l][im_l][ik_l][ik_l][0]*P5_Vmu_P5[im_s][im_s][ik_s][ik_s][0];
	
	//calculate 2*sqrt(E_K*E_Pi)
	jack A=2*sqrt(E2pt[im_s][ik_s]*E2pt[im_l][ik_l]);
	A=2*sqrt(E3pt[im_s][ik_s]*E3pt[im_l][ik_l]);
	//calculate the two ratios
	for(int mu=0;mu<2;mu++)
	  {
	    cV[mu]=sqrt(cV[mu]/Den);
            //reassign the correct sign
	    if(P5_Vmu_P5[im_l][im_s][ik_l][ik_s][mu][T/4][njack]<0) cV[mu]=-cV[mu];
	    //multiply by 2*sqrt(E_K*E_Pi)
	    cV[mu]*=A;
	    //cV[mu].print_to_file("/tmp/out_double_ratio_simm_P5_V%d_P5_%d_%d",mu,ik_s,ik_l);
	  }
	
	jack V0=constant_fit(cV[0],10,14);
	jack Vi=constant_fit(cV[1],10,14);
	
	cout<<"V0: "<<V0/A<<" "<<A<<endl;
	cout<<"Vi: "<<Vi/A<<endl;
	
	//calculate the impulses
	jack Q[4],P[4],Q2,P2;
	calculate_Q_P(Q,P,E2pt[im_s][ik_s],theta[ik_s],E2pt[im_l][ik_l],theta[ik_l]);
	Q2=quad_jack_prod_quad_jack(Q,Q);
	P2=quad_jack_prod_quad_jack(P,P);
	calculate_Q_P(Q,P,E3pt[im_s][ik_s],theta[ik_s],E3pt[im_l][ik_l],theta[ik_l]);
	
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
	    
	    cout<<"Solved in the b case"<<endl;
	  }
	if(ik_s!=ik_l)
	  {
	    cout<<ik_s<<" "<<ik_l<<" Q0="<<Q[0].med()<<" "<<"Qi="<<Q[1].med()<<" P0="<<P[0].med()<<" P1="<<P[1].med()<<endl;
	    file_fp<<Q2.med()<<" "<<fp<<" "<<ik_s<<" "<<ik_l<<" "<<P[0].med()<<endl;
	    file_fm<<Q2.med()<<" "<<fm<<" "<<ik_s<<" "<<ik_l<<" "<<P2.med()<<endl;
	    file_f0<<Q2.med()<<" "<<f0<<" "<<ik_s<<" "<<ik_l<<" "<<P2.med()<<endl;
	  }
      }
  
  return 0;
}
